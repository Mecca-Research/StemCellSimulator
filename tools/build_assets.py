#!/usr/bin/env python3
"""Derive the simulator's committed assets from the raw reference micrographs.

Run ``python tools/fetch_micrographs.py`` first, then::

    python tools/build_assets.py            # all targets
    python tools/build_assets.py hipsc npc  # selected targets

Outputs (assets/derived/):

hipsc_volume.jpg     8x8 mosaic of the 60 z-slices of the 3D stem-cell stack
                     (R = membrane, G = DNA; JPEG 4:4:4 so channels stay separate).
hipsc_labels.png     same mosaic, segmentation code (cell id in bits 0-6,
                     +128 inside the nucleus); lossless.
hipsc_meshes.bin     little-endian mesh buffer: uint16-quantised vertices and
                     uint16 indices (see hipsc_cells.json "meshes").
hipsc_cells.json     per-cell measurements + mesh offsets + calibration summary.
kidney_volume.jpg    4x4 mosaic of the 16-slice kidney stack (RGB channels).
kidney_volume.json   volume metadata.
mitosis.png/.json    DNA image + segmented nuclei, mitotic classification.
skin.jpg/.json       H&E reference + hematoxylin nuclei positions + layer depth.
npc_timelapse.png    15-frame x 2-channel mosaic of the nuclear-targeting movie.
npc_kinetics.json    rim-enrichment time course and first-order fit.
cell_qpi.png         single suspended cell (quantitative phase).
rbc_*.jpg/rbc.json   blood smears; RBC internal ruler and absorbance profile
                     compared with the Evans-Fung 3D thickness.
organoids.jpg/.json  intestinal organoids with Hough-detected circles.
calibration.json     every measured number the JS models consume.
"""
from __future__ import annotations

import json
import struct
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
from skimage import color, feature, filters, measure, morphology, segmentation, transform

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "assets" / "raw"
OUT = ROOT / "assets" / "derived"

HIPSC_SPACING = (0.29, 0.26, 0.26)  # z, y, x in micrometres (scikit-image cells3d docs)
KIDNEY_SPACING = (1.25, 1.24, 1.24)


# --------------------------------------------------------------------------- utils
def pnorm(a: np.ndarray, lo: float = 0.5, hi: float = 99.8) -> np.ndarray:
    """Percentile-normalise to [0, 1]."""
    a = a.astype(np.float32)
    p0, p1 = np.percentile(a, [lo, hi])
    return np.clip((a - p0) / max(p1 - p0, 1e-9), 0.0, 1.0)


def to_u8(a: np.ndarray) -> np.ndarray:
    return np.round(np.clip(a, 0, 1) * 255).astype(np.uint8)


def mosaic(slices: np.ndarray, cols: int) -> np.ndarray:
    """Tile (Z, H, W, C) slices row-major into a (rows*H, cols*W, C) image."""
    z, h, w = slices.shape[:3]
    rows = -(-z // cols)
    out = np.zeros((rows * h, cols * w) + slices.shape[3:], slices.dtype)
    for i in range(z):
        r, c = divmod(i, cols)
        out[r * h:(r + 1) * h, c * w:(c + 1) * w] = slices[i]
    return out


def save_jpeg(path: Path, arr: np.ndarray, quality: int = 88) -> None:
    # 4:4:4 (no chroma subsampling) keeps each fluorescence channel independent.
    Image.fromarray(arr).save(path, quality=quality, subsampling=0, optimize=True)


def save_png(path: Path, arr: np.ndarray) -> None:
    # optimize=True keeps the files small; no alpha channel so browsers never
    # premultiply the data channels on decode.
    Image.fromarray(arr).save(path, optimize=True)


def dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n")


def r3(x) -> float:
    return float(round(float(x), 3))


# ------------------------------------------------------------------ segmentation
def segment_nuclei_3d(dna: np.ndarray, spacing=HIPSC_SPACING) -> np.ndarray:
    """Smoothed Otsu mask, distance transform, h-maxima seeded watershed."""
    smooth = filters.gaussian(dna, sigma=(1, 2, 2))
    mask = smooth > filters.threshold_otsu(smooth)
    mask = morphology.remove_small_objects(mask, max_size=2000)
    mask = ndi.binary_fill_holes(mask)
    dist = ndi.distance_transform_edt(mask, sampling=spacing)
    dist = filters.gaussian(dist, sigma=(1, 3, 3))
    markers = measure.label(morphology.h_maxima(dist, 1.0))
    return segmentation.watershed(-dist, markers, mask=mask).astype(np.int32)


def cell_layer(nuclei: np.ndarray) -> tuple[int, int]:
    area = (nuclei > 0).sum(axis=(1, 2))
    peak = area.max()
    z_lo = int(np.argmax(area > 0.01 * peak))
    z_hi = int(len(area) - 1 - np.argmax(area[::-1] > 0.05 * peak))
    return max(z_lo - 4, 0), min(z_hi + 2, len(area) - 1)


def segment_cells_3d(membrane: np.ndarray, nuclei: np.ndarray) -> np.ndarray:
    """Membrane-intensity watershed seeded by the nuclei, restricted to the monolayer."""
    z0, z1 = cell_layer(nuclei)
    layer = np.zeros(membrane.shape, bool)
    layer[z0:z1 + 1] = True
    smooth = filters.gaussian(membrane, sigma=(1.0, 2, 2))
    return segmentation.watershed(smooth, nuclei, mask=layer, compactness=0.001).astype(np.int32)


# ----------------------------------------------------------------------- meshing
def taubin(verts: np.ndarray, faces: np.ndarray, iters: int = 12,
           lam: float = 0.5, mu: float = -0.53) -> np.ndarray:
    """Taubin (lambda|mu) smoothing: removes voxel stair-steps without shrinking."""
    n = len(verts)
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges = np.concatenate([edges, edges[:, ::-1]])
    deg = np.bincount(edges[:, 0], minlength=n).astype(np.float64)
    deg[deg == 0] = 1
    v = verts.astype(np.float64).copy()
    for _ in range(iters):
        for k in (lam, mu):
            nb = np.zeros_like(v)
            np.add.at(nb, edges[:, 0], v[edges[:, 1]])
            v += k * (nb / deg[:, None] - v)
    return v


def mesh_label(mask: np.ndarray, spacing, step: int = 2, sigma: float = 1.0):
    """Marching cubes on a smoothed binary mask -> (verts[x,y,z] um, faces)."""
    # pad so that (a) the surface is closed at the image border and (b) every
    # axis length is 1 mod step, otherwise a coarse marching-cubes grid skips
    # the last (zero) samples and leaves the mesh open
    lo = step + 1
    hi = [lo + (-(n + 2 * lo - 1)) % step for n in mask.shape]
    pad = np.pad(mask.astype(np.float32), [(lo, h) for h in hi])
    soft = filters.gaussian(pad, sigma=sigma)
    verts, faces, _, _ = measure.marching_cubes(soft, level=0.5, spacing=spacing,
                                                step_size=step, allow_degenerate=False)
    verts -= np.array(spacing) * lo  # undo padding
    verts = taubin(verts, faces)[:, ::-1]  # zyx -> xyz
    # orient every closed surface outward (positive signed volume)
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    if np.einsum("ij,ij->i", v0, np.cross(v1, v2)).sum() < 0:
        faces = faces[:, [0, 2, 1]]
    return verts.astype(np.float32), faces.astype(np.uint32)


def ellipsoid_axes(region) -> dict:
    """Principal semi-axes (um) and orientation of a regionprops object."""
    evals, evecs = np.linalg.eigh(region.inertia_tensor)
    # For a solid ellipsoid, inertia eigenvalue I_a = m/5 (b^2 + c^2); regionprops
    # normalises by mass, so semi-axis^2 = 5/2 (I_b + I_c - I_a).
    i1, i2, i3 = evals
    semi = np.sqrt(np.maximum(2.5 * np.array([i2 + i3 - i1, i1 + i3 - i2, i1 + i2 - i3]), 1e-9))
    return {"semi_axes_um": [r3(s) for s in semi],
            "axes_zyx": [[r3(c) for c in evecs[:, k]] for k in range(3)]}


# --------------------------------------------------------------------- targets
def build_hipsc() -> dict:
    stack = tifffile.imread(RAW / "cells3d.tif")
    mem, dna = pnorm(stack[:, 0]), pnorm(stack[:, 1])
    raw_dna = stack[:, 1].astype(np.float32)
    nuclei = segment_nuclei_3d(dna)
    cells = segment_cells_3d(mem, nuclei)
    assert nuclei.max() < 128, "segmentation code packs ids into 7 bits"

    # --- volume mosaic -----------------------------------------------------
    code = cells.astype(np.uint8) | np.where(nuclei > 0, 128, 0).astype(np.uint8)
    vol = np.stack([to_u8(mem), to_u8(dna), np.zeros_like(code)], axis=-1)
    save_jpeg(OUT / "hipsc_volume.jpg", mosaic(vol, cols=8))
    save_png(OUT / "hipsc_labels.png", mosaic(code, cols=8))

    # --- per-cell measurements and meshes ---------------------------------
    z_dim, y_dim, x_dim = cells.shape
    extent = [x_dim * HIPSC_SPACING[2], y_dim * HIPSC_SPACING[1], z_dim * HIPSC_SPACING[0]]
    nuc_props = {p.label: p for p in measure.regionprops(nuclei, intensity_image=raw_dna, spacing=HIPSC_SPACING)}
    cell_props = {p.label: p for p in measure.regionprops(cells, spacing=HIPSC_SPACING)}
    median_dna = float(np.median([p.intensity_mean for p in nuc_props.values()]))

    # contact graph: labels sharing a face in the 3D cell segmentation
    contacts: dict[int, set[int]] = {k: set() for k in cell_props}
    for axis in range(3):
        a = np.moveaxis(cells, axis, 0)
        left, right = a[:-1].ravel(), a[1:].ravel()
        sel = (left != right) & (left > 0) & (right > 0)
        for i, j in set(zip(left[sel].tolist(), right[sel].tolist())):
            contacts[i].add(j)
            contacts[j].add(i)

    blobs: list[bytes] = []
    offset = 0
    records = []
    # vertices are quantised to uint16 over a fixed box that contains the stack
    qbox = np.array(extent, np.float64) + 4.0  # small margin for smoothing overshoot
    qorigin = np.array([-2.0, -2.0, -2.0])
    for lab in sorted(cell_props):
        cp, npp = cell_props[lab], nuc_props[lab]
        bbox_touches = (cp.bbox[1] == 0 or cp.bbox[2] == 0 or cp.bbox[4] == y_dim or cp.bbox[5] == x_dim)
        nuc_touches = (npp.bbox[1] == 0 or npp.bbox[2] == 0 or npp.bbox[4] == y_dim or npp.bbox[5] == x_dim)
        mesh_meta = {}
        for kind, labels in (("cell", cells), ("nucleus", nuclei)):
            # ~1 um sampling for cells (18 um across), ~0.8 um for nuclei (11 um)
            verts, faces = mesh_label(labels == lab, HIPSC_SPACING,
                                      step=4 if kind == "cell" else 3,
                                      sigma=1.4 if kind == "cell" else 1.0)
            assert len(verts) < 65536
            q = np.round((verts - qorigin) / qbox * 65535).clip(0, 65535).astype("<u2")
            vb, fb = q.tobytes(), faces.astype("<u2").tobytes()
            mesh_meta[kind] = {"vertex_offset": offset, "vertex_count": int(len(verts)),
                               "index_offset": offset + len(vb), "index_count": int(faces.size)}
            pad = (-(len(vb) + len(fb))) % 4  # keep every block 4-byte aligned
            blobs += [vb, fb, b"\0" * pad]
            offset += len(vb) + len(fb) + pad
        brightness = npp.intensity_mean / median_dna
        records.append({
            "id": int(lab),
            "centroid_um": [r3(c) for c in npp.centroid[::-1]],  # x, y, z
            "cell_centroid_um": [r3(c) for c in cp.centroid[::-1]],
            "nucleus_volume_um3": r3(npp.area),
            "cell_volume_um3": r3(cp.area),
            "nc_ratio": r3(npp.area / max(cp.area - npp.area, 1e-9)),
            "nucleus": ellipsoid_axes(npp),
            "cell_height_um": r3((cp.bbox[3] - cp.bbox[0]) * HIPSC_SPACING[0]),
            "dna_brightness": r3(brightness),
            # condensed chromatin is >1.5x brighter than interphase nuclei
            "mitotic": bool(brightness > 1.5 and npp.solidity < 0.85),
            "solidity": r3(npp.solidity),
            "touches_border": bool(bbox_touches),
            "nucleus_touches_border": bool(nuc_touches),
            "neighbours": sorted(int(n) for n in contacts[lab]),
            "mesh": mesh_meta,
        })
    (OUT / "hipsc_meshes.bin").write_bytes(b"".join(blobs))

    interior = [r for r in records if not r["touches_border"]] or records
    whole_nuclei = [r for r in records if not r["nucleus_touches_border"] and not r["mitotic"]]
    summary = {
        "n_cells": len(records),
        "n_interior": len(interior),
        "n_whole_nuclei": len(whole_nuclei),
        "n_mitotic": sum(r["mitotic"] for r in records),
        "nucleus_volume_um3_median": r3(np.median([r["nucleus_volume_um3"] for r in whole_nuclei])),
        "cell_volume_um3_median": r3(np.median([r["cell_volume_um3"] for r in interior])),
        "nc_ratio_median": r3(np.median([r["nc_ratio"] for r in interior])),
        "cell_height_um_median": r3(np.median([r["cell_height_um"] for r in interior])),
        "neighbours_mean_interior": r3(np.mean([len(r["neighbours"]) for r in interior])),
        "nucleus_radius_um_equiv": r3(np.cbrt(3 * np.median([r["nucleus_volume_um3"] for r in whole_nuclei]) / (4 * np.pi))),
        "nucleus_semi_axes_um_median": [r3(np.median([r["nucleus"]["semi_axes_um"][k] for r in whole_nuclei])) for k in range(3)],
        "cell_radius_um_equiv": r3(np.cbrt(3 * np.median([r["cell_volume_um3"] for r in interior]) / (4 * np.pi))),
        "areal_density_per_100um2": r3(100 * len(records) / (extent[0] * extent[1])),
        # in-plane nearest-neighbour spacing of nucleus centroids (whole nuclei only)
        "nn_distance_um_median": r3(np.median([
            min(np.hypot(r["centroid_um"][0] - q["centroid_um"][0], r["centroid_um"][1] - q["centroid_um"][1])
                for q in records if q is not r)
            for r in records if not r["nucleus_touches_border"]])),
    }
    meta = {
        "source": "cells3d.tif (Allen Institute for Cell Science, CC0)",
        "shape_zyx": list(cells.shape),
        "spacing_zyx_um": list(HIPSC_SPACING),
        "extent_xyz_um": [r3(e) for e in extent],
        "mosaic": {"file": "hipsc_volume.jpg", "labels": "hipsc_labels.png", "cols": 8, "slices": z_dim,
                   "tile": [x_dim, y_dim], "channels": {"r": "membrane", "g": "dna"},
                   "label_code": "cell id | 128 * inside-nucleus"},
        "meshes": {"file": "hipsc_meshes.bin", "vertex": "uint16 xyz, um = origin + q / 65535 * box",
                   "quant_origin_um": [r3(v) for v in qorigin], "quant_box_um": [r3(v) for v in qbox],
                   "index": "uint16 triangle list"},
        "summary": summary,
        "cells": records,
    }
    dump(OUT / "hipsc_cells.json", meta)
    return {"hipsc": summary}


def build_kidney() -> dict:
    stack = tifffile.imread(RAW / "kidney-tissue-fluorescence.tif").astype(np.float32)  # Z C Y X
    z, c, h, w = stack.shape
    small = stack.reshape(z, c, h // 2, 2, w // 2, 2).mean(axis=(3, 5))
    chans = [pnorm(small[:, k], 1, 99.7) for k in range(c)]
    vol = np.stack([to_u8(ch) for ch in chans], axis=-1)  # Z Y X 3
    save_jpeg(OUT / "kidney_volume.jpg", mosaic(vol, cols=4))
    spacing = (KIDNEY_SPACING[0], KIDNEY_SPACING[1] * 2, KIDNEY_SPACING[2] * 2)
    dump(OUT / "kidney_volume.json", {
        "source": "kidney-tissue-fluorescence.tif (G. Buckley, CC0)",
        "shape_zyx": [z, h // 2, w // 2], "spacing_zyx_um": list(spacing),
        "extent_xyz_um": [r3(w // 2 * spacing[2]), r3(h // 2 * spacing[1]), r3(z * spacing[0])],
        "mosaic": {"file": "kidney_volume.jpg", "cols": 4, "slices": z, "tile": [w // 2, h // 2],
                   "channels": {"r": "450 nm", "g": "515 nm", "b": "605 nm"}},
    })
    return {"kidney": {"slices": z}}


def build_mitosis() -> dict:
    img = tifffile.imread(RAW / "AS_09125_050116030001_D03f00d0.tif").astype(np.float32)
    # Three-class Otsu separates background from nuclei. A nucleus is scored as
    # mitotic (condensed chromatin) when its 90th-percentile DNA intensity is an
    # outlier: > median + 4 * MAD over all nuclei in the field.
    t_lo, t_hi = filters.threshold_multiotsu(img, classes=3)
    nuclei_mask = ndi.binary_fill_holes(img > t_lo)
    nuclei_mask = morphology.remove_small_objects(nuclei_mask, max_size=20)
    dist = ndi.distance_transform_edt(nuclei_mask)
    markers = measure.label(morphology.h_maxima(filters.gaussian(dist, 1), 1.0))
    labels = segmentation.watershed(-dist, markers, mask=nuclei_mask)
    regions = measure.regionprops(labels, intensity_image=img)
    q90 = np.array([np.percentile(img[labels == p.label], 90) for p in regions])
    med = float(np.median(q90))
    mad = float(np.median(np.abs(q90 - med)) * 1.4826)
    cut = med + 4 * mad
    div_ids = {p.label for p, q in zip(regions, q90) if q > cut}
    nuclei = []
    for p in regions:
        y, x = p.centroid
        nuclei.append({"x": r3(x), "y": r3(y), "r": r3(np.sqrt(p.area / np.pi)),
                       "orientation": r3(p.orientation), "ecc": r3(p.eccentricity),
                       "mitotic": p.label in div_ids})
    n_div = sum(n["mitotic"] for n in nuclei)
    mi = n_div / max(len(nuclei), 1)
    save_png(OUT / "mitosis.png", to_u8(pnorm(img, 0.5, 99.9)))
    dump(OUT / "mitosis.json", {"source": "CellProfiler ExampleHuman (Moffat et al. 2006), CC0",
                                "size": list(img.shape[::-1]), "otsu": [r3(t_lo), r3(t_hi)],
                                "mitotic_rule": "q90 > median + 4*MAD", "q90_cut": r3(cut),
                                "n_nuclei": len(nuclei), "n_mitotic": n_div,
                                "mitotic_index": r3(mi), "nuclei": nuclei})
    return {"mitosis": {"n_nuclei": len(nuclei), "n_mitotic": n_div, "mitotic_index": r3(mi)}}


def build_skin() -> dict:
    im = Image.open(RAW / "Normal_Epidermis_and_Dermis_with_Intradermal_Nevus_10x.JPG").convert("RGB")
    im = im.resize((1024, 768), Image.LANCZOS)
    im.save(OUT / "skin.jpg", quality=88, optimize=True)
    rgb = np.asarray(im).astype(np.float32) / 255
    hed = color.rgb2hed(rgb)
    hema = filters.gaussian(hed[..., 0], 1.5)
    t = filters.threshold_otsu(hema)
    blobs = measure.label(morphology.remove_small_objects(hema > t, max_size=12))
    nuclei = [{"x": r3(p.centroid[1]), "y": r3(p.centroid[0]), "r": r3(np.sqrt(p.area / np.pi))}
              for p in measure.regionprops(blobs) if 12 < p.area < 600]
    # epidermis = the hematoxylin-dense band; report its mean depth profile per column
    dens = filters.gaussian(hed[..., 0] > t, 6)
    col_profile = [int(np.argmax(dens[:, x] > 0.3)) if (dens[:, x] > 0.3).any() else -1
                   for x in range(0, 1024, 16)]
    dump(OUT / "skin.json", {"source": "Kilbad, Wikimedia Commons, public domain",
                             "size": [1024, 768], "n_nuclei": len(nuclei), "nuclei": nuclei,
                             "epidermis_top_px_every16": col_profile})
    return {"skin": {"n_nuclei": len(nuclei)}}


def first_order(t, e_inf, e0, tau):
    return e_inf - (e_inf - e0) * np.exp(-t / tau)


def build_npc() -> dict:
    movie = tifffile.imread(RAW / "NPCsingleNucleus.tif").astype(np.float32)  # T C Y X
    frames = movie.shape[0]
    enrich = []
    for t in range(frames):
        dna, marker = movie[t, 0], movie[t, 1]
        sm = filters.gaussian(dna, 2)
        lab = measure.label(sm > filters.threshold_otsu(sm))
        nucleus = lab == max(measure.regionprops(lab), key=lambda p: p.area).label
        nucleus = ndi.binary_fill_holes(nucleus)
        grow = lambda r: ndi.binary_dilation(nucleus, morphology.disk(r))
        rim = grow(3) & ~ndi.binary_erosion(nucleus, morphology.disk(3))
        cyto = grow(14) & ~grow(5)
        enrich.append(float(marker[rim].mean() / marker[cyto].mean()))
    t = np.arange(frames, dtype=float)
    e = np.array(enrich)
    popt, _ = curve_fit(first_order, t, e, p0=[e.max(), e[0], frames / 3], maxfev=20000)
    pred = first_order(t, *popt)
    r2 = 1 - np.sum((e - pred) ** 2) / np.sum((e - e.mean()) ** 2)
    tiles = []
    for k in range(frames):
        tiles.append(np.concatenate([pnorm(movie[k, 0], 0.5, 99.9), pnorm(movie[k, 1], 0.5, 99.9)], axis=1))
    save_png(OUT / "npc_timelapse.png", to_u8(mosaic(np.stack(tiles), cols=5)))
    result = {"frames": frames, "rim_over_cytoplasm": [r3(x) for x in e],
              "fit": {"model": "E(t) = E_inf - (E_inf - E0) exp(-t / tau)",
                      "E_inf": r3(popt[0]), "E0": r3(popt[1]), "tau_frames": r3(popt[2]), "r2": r3(r2)},
              "note": "The frame interval is not recorded in the TIFF; tau is expressed in frames and "
                      "used as a relative (shape) calibration of the nuclear-targeting step.",
              "tile": [movie.shape[3] * 2, movie.shape[2]], "cols": 5,
              "source": "Boni et al. 2015 J Cell Biol (Ellenberg lab), CC0"}
    dump(OUT / "npc_kinetics.json", result)
    return {"npc": result["fit"]}


def build_cell() -> dict:
    src = RAW / "cell.png"
    im = Image.open(src).convert("L")
    im.save(OUT / "cell_qpi.png", optimize=True)
    a = np.asarray(im).astype(np.float32)
    mask = a > filters.threshold_otsu(a)
    mask = morphology.remove_small_objects(ndi.binary_fill_holes(mask), max_size=500)
    p = max(measure.regionprops(measure.label(mask)), key=lambda q: q.area)
    diameter_um = 2 * np.sqrt(p.area / np.pi) * 0.107
    return {"single_cell": {"equivalent_diameter_um": r3(diameter_um)}}


EVANS_FUNG = (7.82, 0.207, 2.003, -1.123)  # D (um), c0, c1, c2


def evans_fung_thickness(x: np.ndarray) -> np.ndarray:
    """Full thickness of the Evans & Fung (1972) biconcave erythrocyte at x = r / R."""
    d, c0, c1, c2 = EVANS_FUNG
    inside = np.clip(1 - x ** 2, 0, 1)
    return np.where(x < 1, d * np.sqrt(inside) * (c0 + c1 * x ** 2 + c2 * x ** 4), 0.0)


def build_rbc() -> dict:
    """Erythrocytes as an internal ruler, and absorbance profile vs 3D model thickness."""
    out = {}
    bins = np.linspace(0, 1.2, 25)
    mids = 0.5 * (bins[1:] + bins[:-1])
    for tag, stem in (("normal", "chula_rbc_1"), ("lymphocyte", "chula_rbc_2")):
        rgb = np.asarray(Image.open(RAW / f"{stem}.jpg").convert("RGB")).astype(np.float64)
        green = rgb[..., 1]                        # eosin/hemoglobin absorb strongly in green
        absorb = -np.log(np.clip(green / np.percentile(green, 95), 1e-3, 1))
        smooth = filters.gaussian(absorb, 1)
        mask = ndi.binary_fill_holes(smooth > 0.8 * filters.threshold_otsu(smooth))
        lab = measure.label(morphology.remove_small_objects(mask, max_size=60))
        pts = np.loadtxt(RAW / f"{stem}.txt", ndmin=2).astype(int)
        yy, xx = np.mgrid[0:lab.shape[0], 0:lab.shape[1]]
        diam, profiles = [], []
        for x, y, cls in pts:
            if cls != 0:
                continue
            region = lab[min(y, lab.shape[0] - 1), min(x, lab.shape[1] - 1)]
            if region == 0:
                continue
            p = measure.regionprops((lab == region).astype(int))[0]
            if p.eccentricity > 0.6 or p.solidity < 0.93:  # isolated, round cells only
                continue
            diam.append(p.equivalent_diameter_area)
            cy, cx = p.centroid
            rr = np.hypot(yy - cy, xx - cx) / (p.equivalent_diameter_area / 2)
            profiles.append([absorb[(rr >= bins[i]) & (rr < bins[i + 1])].mean() for i in range(24)])
        d_px = float(np.median(diam))
        entry = {"n_cells": len(diam), "diameter_px_median": r3(d_px),
                 "diameter_px_iqr": [r3(v) for v in np.percentile(diam, [25, 75])],
                 "um_per_px": r3(EVANS_FUNG[0] / d_px)}
        if tag == "normal":
            prof = np.nanmedian(np.array(profiles), axis=0)
            prof /= prof.max()
            model = evans_fung_thickness(mids)
            model /= model.max()
            inside = mids < 0.95
            rmse = float(np.sqrt(np.mean((prof[inside] - model[inside]) ** 2)))
            corr = float(np.corrcoef(prof[inside], model[inside])[0, 1])
            entry.update({"r_over_R": [r3(v) for v in mids], "absorbance_profile": [r3(v) for v in prof],
                          "evans_fung_thickness": [r3(v) for v in model],
                          "profile_rmse": r3(rmse), "profile_correlation": r3(corr),
                          "pallor_min_over_max": r3(prof[0])})
        else:
            # the one leukocyte in this field (a lymphocyte) is not annotated; its
            # hematoxylin-dense nucleus is the only large, round, dark-violet blob
            # (the thin dark diagonal is a slide scratch and is rejected by shape)
            dark = filters.gaussian(rgb[..., 1], 1.5) < np.percentile(rgb[..., 1], 1.5)
            blobs = [b for b in measure.regionprops(measure.label(ndi.binary_fill_holes(dark)))
                     if b.area > 150 and b.eccentricity < 0.8 and b.solidity > 0.9]
            if blobs:
                big = max(blobs, key=lambda b: b.area)
                entry["lymphocyte_nucleus_diameter_um"] = r3(big.equivalent_diameter_area * EVANS_FUNG[0] / d_px)
                entry["lymphocyte_centre_px"] = [r3(big.centroid[1]), r3(big.centroid[0])]
        out[tag] = entry
        Image.fromarray(rgb.astype(np.uint8)).save(OUT / f"rbc_{tag}.jpg", quality=90)
    dump(OUT / "rbc.json", {"source": "Chula-RBC-12-Dataset (MIT)", **out,
                            "method": "RBC diameter = 7.82 um internal ruler; absorbance -ln(I/I0) in green; "
                                      "median radial profile of isolated normal cells vs Evans-Fung thickness"})
    return {"rbc": {k: {kk: vv for kk, vv in v.items() if not isinstance(vv, list)} for k, v in out.items()}}


def build_organoids() -> dict:
    im = Image.open(RAW / "orgaquant_3.jpg").convert("L")
    width = 750
    im = im.resize((width, round(width * im.height / im.width)), Image.LANCZOS)
    a = np.asarray(im).astype(np.float64) / 255
    flat = a / filters.gaussian(a, 40)           # remove vignetting
    edges = feature.canny(flat, sigma=2.0, low_threshold=0.03, high_threshold=0.08)
    radii = np.arange(9, 46)
    acc, cx, cy, rad = transform.hough_circle_peaks(transform.hough_circle(edges, radii), radii,
                                                   min_xdistance=8, min_ydistance=8,
                                                   threshold=0.45, total_num_peaks=400)
    centre = np.array([a.shape[1] / 2, a.shape[0] / 2])
    keep = []
    for i in np.argsort(-acc):
        p = np.array([cx[i], cy[i]])
        if np.linalg.norm(p - centre) > 0.47 * a.shape[1]:
            continue
        if any(np.linalg.norm(p - np.array([cx[j], cy[j]])) < 0.8 * max(rad[i], rad[j]) for j in keep):
            continue
        keep.append(i)
    circles = [{"x": int(cx[i]), "y": int(cy[i]), "r": int(rad[i])} for i in keep]
    r = np.array([c["r"] for c in circles], float)
    Image.fromarray(to_u8(pnorm(flat, 0.5, 99.5))).save(OUT / "organoids.jpg", quality=88)
    stats = {"n_organoids": len(circles), "radius_px_median": r3(np.median(r)),
             "radius_cv": r3(r.std() / r.mean()), "radius_p90_over_median": r3(np.percentile(r, 90) / np.median(r))}
    dump(OUT / "organoids.json", {"source": "OrgaQuant test data (MIT)", "size": [a.shape[1], a.shape[0]],
                                  "method": "vignetting-flattened Canny edges + Hough circles (r 9-45 px)",
                                  **stats, "circles": circles})
    return {"organoids": stats}


TARGETS = {"hipsc": build_hipsc, "kidney": build_kidney, "mitosis": build_mitosis,
           "skin": build_skin, "npc": build_npc, "cell": build_cell, "rbc": build_rbc,
           "organoids": build_organoids}


def main(argv: list[str]) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    selected = argv or list(TARGETS)
    cal_path = OUT / "calibration.json"
    calibration = json.loads(cal_path.read_text()) if cal_path.exists() else {}
    for name in selected:
        print(f"building {name} ...", flush=True)
        calibration.update(TARGETS[name]())
    dump(cal_path, calibration)
    print(json.dumps(calibration, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
