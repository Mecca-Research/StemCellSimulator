"""Tests for the micrograph pipeline (tools/build_assets.py, fetch_micrographs.py).

Synthetic tests always run; checks of the committed derived assets run against
assets/derived/; raw-data checks are skipped when assets/raw/ is absent.
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
DERIVED = ROOT / "assets" / "derived"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ba = load("build_assets")
fm = load("fetch_micrographs")


# ----------------------------------------------------------------- catalogue
def test_catalog_entries_are_pinned_and_open():
    entries = fm.load_catalog()
    assert len(entries) >= 6
    ids = set()
    for e in entries:
        assert e["id"] not in ids
        ids.add(e["id"])
        assert re.fullmatch(r"[0-9a-f]{64}", e["sha256"]), e["id"]
        assert e["url"].startswith("https://")
        assert e["license"] in {"CC0-1.0", "Public domain", "MIT"}, e["license"]
        for key in ("title", "credit", "used_for"):
            assert e[key]


def test_raw_files_match_their_pins():
    raw = ROOT / "assets" / "raw"
    present = [e for e in fm.load_catalog() if (raw / e["file"]).exists()]
    if not present:
        pytest.skip("raw micrographs not fetched")
    for e in present:
        assert fm.sha256(raw / e["file"]) == e["sha256"], e["id"]


# ----------------------------------------------------------------- geometry
def ball(shape, centre, radius):
    zz, yy, xx = np.indices(shape)
    return ((zz - centre[0]) ** 2 + (yy - centre[1]) ** 2 + (xx - centre[2]) ** 2) <= radius ** 2


def signed_volume(verts, faces):
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(np.einsum("ij,ij->i", v0, np.cross(v1, v2)).sum() / 6)


@pytest.mark.parametrize("step", [1, 2, 4])
def test_mesh_of_a_ball_is_closed_outward_and_has_the_right_volume(step):
    mask = ball((40, 44, 48), (20, 22, 24), 12)
    verts, faces = ba.mesh_label(mask, (1.0, 1.0, 1.0), step=step, sigma=1.0)
    vol = signed_volume(verts, faces)
    assert vol > 0, "normals must point outward"
    assert abs(vol - mask.sum()) / mask.sum() < 0.12
    # a closed surface: every edge is shared by exactly two triangles
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert np.all(counts == 2)


def test_mesh_touching_the_image_border_is_still_closed():
    mask = ball((30, 30, 30), (15, 15, 2), 9)  # cut by the x = 0 face
    verts, faces = ba.mesh_label(mask, (1.0, 1.0, 1.0), step=4, sigma=1.2)
    assert signed_volume(verts, faces) > 0
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert np.all(counts == 2)


def test_ellipsoid_axes_recovers_semi_axes():
    from skimage import measure
    zz, yy, xx = np.indices((40, 60, 80))
    mask = ((zz - 20) / 6) ** 2 + ((yy - 30) / 10) ** 2 + ((xx - 40) / 18) ** 2 <= 1
    props = measure.regionprops(mask.astype(int))[0]
    semi = sorted(ba.ellipsoid_axes(props)["semi_axes_um"], reverse=True)
    assert np.allclose(semi, [18, 10, 6], rtol=0.06)


def test_nuclei_segmentation_separates_touching_blobs():
    img = np.zeros((24, 64, 96), np.float32)
    for cx in (28, 50, 72):
        img[ball(img.shape, (12, 32, cx), 9)] = 1.0
    labels = ba.segment_nuclei_3d(img, spacing=(1.0, 1.0, 1.0))
    assert labels.max() == 3


def test_mosaic_layout():
    s = np.arange(5 * 2 * 3).reshape(5, 2, 3)
    m = ba.mosaic(s, cols=2)
    assert m.shape == (6, 6)
    assert np.array_equal(m[2:4, 3:6], s[3])


def test_evans_fung_profile():
    x = np.linspace(0, 0.999, 400)
    t = ba.evans_fung_thickness(x)
    assert abs(t[0] - 7.82 * 0.207) < 1e-9          # central thickness D * c0
    assert 0.6 < x[np.argmax(t)] < 0.9               # thickest at the rim torus
    assert ba.evans_fung_thickness(np.array([1.0]))[0] == 0


def test_first_order_fit_recovers_parameters():
    from scipy.optimize import curve_fit
    t = np.arange(15, dtype=float)
    y = ba.first_order(t, 3.0, 1.2, 6.0)
    popt, _ = curve_fit(ba.first_order, t, y, p0=[2.5, 1.0, 4.0])
    assert np.allclose(popt, [3.0, 1.2, 6.0], rtol=1e-4)


# ----------------------------------------------------------------- committed assets
def test_hipsc_mesh_buffer_matches_metadata():
    meta = json.loads((DERIVED / "hipsc_cells.json").read_text())
    size = (DERIVED / "hipsc_meshes.bin").stat().st_size
    assert meta["summary"]["n_cells"] == len(meta["cells"]) > 10
    for c in meta["cells"]:
        for kind in ("cell", "nucleus"):
            m = c["mesh"][kind]
            assert m["vertex_offset"] % 2 == 0 and m["index_offset"] % 2 == 0
            assert m["index_offset"] + 2 * m["index_count"] <= size
            assert m["index_count"] % 3 == 0
        assert c["nucleus_volume_um3"] < c["cell_volume_um3"]


def test_calibration_values_are_plausible():
    cal = json.loads((DERIVED / "calibration.json").read_text())
    h = cal["hipsc"]
    assert 300 < h["nucleus_volume_um3_median"] < 2000
    assert 0.1 < h["nc_ratio_median"] < 1.5
    assert 5 < h["cell_height_um_median"] < 20
    assert 0.01 < cal["mitosis"]["mitotic_index"] < 0.2
    assert cal["npc"]["r2"] > 0.9
    rbc = cal["rbc"]["normal"]
    assert rbc["profile_correlation"] > 0.9            # 3D red-cell model vs real smear
    assert 5 < cal["rbc"]["lymphocyte"]["lymphocyte_nucleus_diameter_um"] < 12
