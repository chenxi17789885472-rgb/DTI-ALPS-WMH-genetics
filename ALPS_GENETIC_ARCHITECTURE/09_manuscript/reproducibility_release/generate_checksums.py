#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "SHA256SUMS.tsv"

TARGETS = [
    PROJECT / "06_scripts",
    PROJECT / "09_manuscript" / "scripts",
    PROJECT / "09_manuscript" / "submission_assets" / "figure_data",
    PROJECT / "09_manuscript" / "submission_assets" / "tables",
]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


rows: list[tuple[str, str]] = []
for target in TARGETS:
    for path in sorted(target.rglob("*")):
        if path.is_file() and not path.name.startswith("."):
            rows.append((str(path.relative_to(PROJECT)), digest(path)))

OUTPUT.write_text(
    "relative_path\tsha256\n"
    + "\n".join(f"{path}\t{sha}" for path, sha in rows)
    + "\n",
    encoding="utf-8",
)
print(f"Wrote {len(rows)} checksums to {OUTPUT}")
