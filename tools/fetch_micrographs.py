#!/usr/bin/env python3
"""Download the reference micrographs listed in assets/micrographs/catalog.json.

Files land in assets/raw/ (git-ignored). Every download is verified against the
SHA-256 pin in the catalog, so the derived assets are reproducible bit-for-bit.

    python tools/fetch_micrographs.py            # fetch everything missing
    python tools/fetch_micrographs.py hipsc-3d   # fetch selected ids
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "assets" / "micrographs" / "catalog.json"
RAW = ROOT / "assets" / "raw"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_catalog() -> list[dict]:
    return json.loads(CATALOG.read_text())["entries"]


def fetch(entry: dict) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / entry["file"]
    if dest.exists() and sha256(dest) == entry["sha256"]:
        print(f"ok      {entry['id']:<16} {dest.name}")
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"fetch   {entry['id']:<16} {entry['url']}")
    with urllib.request.urlopen(entry["url"], timeout=600) as resp, tmp.open("wb") as out:
        while chunk := resp.read(1 << 20):
            out.write(chunk)
    digest = sha256(tmp)
    if digest != entry["sha256"]:
        tmp.unlink()
        raise SystemExit(f"checksum mismatch for {entry['id']}: {digest} != {entry['sha256']}")
    tmp.replace(dest)
    return dest


def main(argv: list[str]) -> int:
    entries = load_catalog()
    wanted = set(argv)
    unknown = wanted - {e["id"] for e in entries}
    if unknown:
        raise SystemExit(f"unknown catalog ids: {sorted(unknown)}")
    for entry in entries:
        if not wanted or entry["id"] in wanted:
            fetch(entry)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
