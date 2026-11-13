#!/bin/sh
set -eu
cd "$(dirname "$0")"
tex_cache="/tmp/etv4-preprint-tex-cache"
mkdir -p build "$tex_cache/bin"
export TEXMFVAR="$tex_cache"
export TEXMFCACHE="$tex_cache"
# TeX Live's universal Biber launcher can fail with Apple's newer lipo.
# Extract the installed arm64 binary locally; Overleaf uses standard Biber.
if [ "$(uname -s)" = Darwin ] && [ "$(uname -m)" = arm64 ]; then
    if ! biber --version >/dev/null 2>&1; then
        /usr/bin/lipo "$(command -v biber)" -thin arm64 -output "$tex_cache/bin/biber"
        chmod +x "$tex_cache/bin/biber"
        export PATH="$tex_cache/bin:$PATH"
    fi
fi
latexmk -lualatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
