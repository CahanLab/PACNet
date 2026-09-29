#!/usr/bin/env bash
# Builds the static PACNet web app: a shinylive export of webApp/ plus the
# reference data in data/. Usage (from repo root): tools/build_site.sh [outdir]
set -euo pipefail
out="${1:-site}"
rm -rf "$out"
Rscript -e "options(timeout = 900)
shinylive::export('webApp', '$out', template_params = list(
  title = 'PACNet',
  include_after_body = paste(readLines('tools/site_loading.html'), collapse = '\\n')))"
cp -R data "$out/data"
touch "$out/.nojekyll"
echo "Built $out ($(du -sh "$out" | cut -f1))"
