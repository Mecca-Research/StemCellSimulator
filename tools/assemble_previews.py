#!/usr/bin/env python3
"""Assemble captured frames (tests/e2e/out/frames/<name>/*.png) into animated
WebP previews in docs/images/<name>.webp plus a still <name>.jpg."""
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
FRAMES = ROOT / "tests" / "e2e" / "out" / "frames"
OUT = ROOT / "docs" / "images"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for d in sorted(p for p in FRAMES.iterdir() if p.is_dir()):
        files = sorted(d.glob("*.png"))
        if not files:
            continue
        frames = []
        for f in files:
            im = Image.open(f).convert("RGB")
            w = 640
            frames.append(im.resize((w, round(im.height * w / im.width)), Image.LANCZOS))
        frames[0].save(OUT / f"{d.name}.webp", save_all=True, append_images=frames[1:], duration=90, loop=0, quality=72, method=6)
        frames[len(frames) // 2].save(OUT / f"{d.name}.jpg", quality=85, optimize=True)
        print(f"{d.name}: {len(frames)} frames -> {(OUT / f'{d.name}.webp').stat().st_size // 1024} KiB")


if __name__ == "__main__":
    main()
