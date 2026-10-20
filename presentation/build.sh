#!/usr/bin/env bash
# Build the eleven-slide Beamer deck -> TELEMEDICON26_CliffForge_Pitch.pdf (works from any directory).
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

# Some macOS TeX installations ship a universal Biber that fails while unpacking
# its own architecture. Use a temporary, locally thinned copy only in that case.
if ! biber --version >/dev/null 2>&1; then
  command -v lipo >/dev/null || { echo "biber cannot run and lipo is unavailable" >&2; exit 1; }
  biber_source="$(command -v biber)"
  biber_arch="$(uname -m)"
  biber_fallback_dir="$(mktemp -d "${TMPDIR:-/private/tmp}/cliffforge-biber.XXXXXX")"
  trap 'rm -f "$biber_fallback_dir/biber"; rmdir "$biber_fallback_dir" 2>/dev/null || true' EXIT
  lipo -thin "$biber_arch" "$biber_source" -output "$biber_fallback_dir/biber" || {
    echo "biber cannot run and its $biber_arch binary could not be extracted" >&2
    exit 1
  }
  chmod +x "$biber_fallback_dir/biber"
  "$biber_fallback_dir/biber" --version >/dev/null || {
    echo "temporary $biber_arch biber cannot run" >&2
    exit 1
  }
  export PATH="$biber_fallback_dir:$PATH"
  echo "Using temporary $biber_arch biber for this build"
fi

mkdir -p build
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex
if command -v pdfinfo >/dev/null; then
  pages="$(pdfinfo build/main.pdf | awk '/^Pages:/{print $2}')"
  echo "pages: $pages"
  [ "$pages" = 11 ] || { echo "expected 11 slides, found $pages" >&2; exit 1; }
fi
cp build/main.pdf TELEMEDICON26_CliffForge_Pitch.pdf
cp build/main.pdf main.pdf
echo "wrote $HERE/TELEMEDICON26_CliffForge_Pitch.pdf"
