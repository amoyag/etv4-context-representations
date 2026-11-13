#!/usr/bin/env python3
"""Compile the existing manuscript, keeping TeX caches inside the workspace."""
import argparse
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", default="logs/paper-build.txt")
    args = parser.parse_args()
    log = (ROOT / args.log).resolve()
    if not log.is_relative_to(ROOT):
        parser.error("Log must stay inside the pilot workspace")
    cache = ROOT / ".cache/texmf"
    cache.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TEXMFVAR=str(cache), TEXMFCACHE=str(cache))
    with log.open("w") as stream:
        result = subprocess.run(["latexmk", "-g", "-lualatex", "-interaction=nonstopmode",
                                 "-halt-on-error", "-outdir=build", "main.tex"],
                                cwd=ROOT / "paper", env=env, stdout=stream, stderr=subprocess.STDOUT)
    print(f"LaTeX exit {result.returncode}; log: {log}")
    return result.returncode

if __name__ == "__main__":
    raise SystemExit(main())
