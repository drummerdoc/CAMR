#!/usr/bin/env bash
# Render every docs/*.md to PDF (pandoc + xelatex):  ./build_pdf.sh  ->  docs/pdf/*.pdf
# Needs pandoc and a TeX distribution with xelatex (MacTeX / BasicTeX + fontspec;
# TeX Live on Linux).  Override the fonts with MAINFONT / MONOFONT if the defaults
# are not installed.
set -euo pipefail
cd "$(dirname "$0")" && mkdir -p pdf
command -v pandoc  >/dev/null || { echo "pandoc not found";  exit 1; }
command -v xelatex >/dev/null || { echo "xelatex not found (install MacTeX or BasicTeX)"; exit 1; }
if [ -z "${MAINFONT:-}" ]; then
    case "$(uname -s)" in
        Darwin) MAINFONT="Times New Roman"; MONOFONT="${MONOFONT:-Menlo}" ;;
        *)      MAINFONT="DejaVu Serif";    MONOFONT="${MONOFONT:-DejaVu Sans Mono}" ;;
    esac
fi
for f in *.md; do
    echo "  $f"
    pandoc "$f" -f gfm+tex_math_dollars --pdf-engine=xelatex \
        -V mainfont="$MAINFONT" -V monofont="$MONOFONT" \
        -V geometry:margin=2.2cm -V fontsize=10pt -V colorlinks=true \
        --toc --toc-depth=2 --metadata title="${f%.md}" \
        -o "pdf/${f%.md}.pdf" 2>&1 | grep -v '^\[WARNING\]' || true
done
echo "wrote docs/pdf/"
