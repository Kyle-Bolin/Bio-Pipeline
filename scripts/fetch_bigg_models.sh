#!/usr/bin/env bash
# Download BiGG models as plain SBML into a folder, for BENCHMARK_MODEL_DIR.
# Usage: scripts/fetch_bigg_models.sh DEST_DIR MODEL_ID [MODEL_ID ...]
set -euo pipefail

dest=$1
shift
mkdir -p "$dest"
for id in "$@"; do
  if [[ -s "$dest/$id.xml" ]]; then
    echo "$id: already present"
    continue
  fi
  # -L follows BiGG's http -> https redirect, which COBRApy's downloader does not.
  curl -fsSL --retry 3 --retry-delay 5 --max-time 120 \
    "https://bigg.ucsd.edu/static/models/$id.xml.gz" | gunzip > "$dest/$id.xml.tmp"
  mv "$dest/$id.xml.tmp" "$dest/$id.xml"
  echo "$id: downloaded ($(du -h "$dest/$id.xml" | cut -f1))"
done
