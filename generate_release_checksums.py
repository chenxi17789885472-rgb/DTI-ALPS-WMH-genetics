#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "SHA256SUMS.txt"
EXCLUDED_NAMES = {"SHA256SUMS.txt", ".DS_Store"}
RENV_LOCAL_DIRS = {"library", "local", "cellar", "lock", "python", "sandbox", "staging"}


def is_local_environment_file(path: Path) -> bool:
    parts = path.relative_to(ROOT).parts
    return any(
        parts[index] == "renv" and parts[index + 1] in RENV_LOCAL_DIRS
        for index in range(len(parts) - 1)
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


files = [
    path
    for path in ROOT.rglob("*")
    if path.is_file()
    and ".git" not in path.parts
    and path.name not in EXCLUDED_NAMES
    and not is_local_environment_file(path)
]

OUTPUT.write_text(
    "".join(
        f"{sha256(path)}  {path.relative_to(ROOT)}\n"
        for path in sorted(files)
    ),
    encoding="utf-8",
)
print(f"Wrote checksums for {len(files)} release files to {OUTPUT}")
