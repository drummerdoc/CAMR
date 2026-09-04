#!/usr/bin/env bash
# Render the docs to standalone HTML (math via MathJax) for viewers without
# Markdown+LaTeX support:  ./build_html.sh  ->  docs/html/*.html
cd "$(dirname "$0")" && mkdir -p html
for f in *.md; do
    pandoc "$f" -f gfm+tex_math_dollars -t html5 -s --mathjax --toc --toc-depth=3 \
           --metadata title="${f%.md}" -o "html/${f%.md}.html"
done
echo "wrote docs/html/"
