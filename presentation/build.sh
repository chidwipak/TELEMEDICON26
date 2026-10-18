#!/usr/bin/env bash
# Build the seven-slide Beamer deck -> TELEMEDICON26_CliffForge_Pitch.pdf (works from any directory).
# Needs: pdflatex + biber + latexmk with the beamer, biblatex, pgf/tikz packages (TeX Live / TinyTeX).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
# TinyTeX keeps biber in its own bin dir, which is not always on PATH
if ! command -v biber >/dev/null 2>&1; then
  for d in "$HOME"/.TinyTeX/bin/*/ /usr/local/texlive/*/bin/*/; do [ -x "${d}biber" ] && export PATH="$d:$PATH" && break; done
fi
command -v latexmk >/dev/null || { echo "latexmk not found" >&2; exit 1; }
command -v biber   >/dev/null || { echo "biber not found (tlmgr install biber)" >&2; exit 1; }
mkdir -p build
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
cp build/main.pdf TELEMEDICON26_CliffForge_Pitch.pdf
if command -v pdfinfo >/dev/null; then echo "pages: $(pdfinfo TELEMEDICON26_CliffForge_Pitch.pdf | awk '/^Pages:/{print $2}')"; fi
echo "wrote $HERE/TELEMEDICON26_CliffForge_Pitch.pdf"
