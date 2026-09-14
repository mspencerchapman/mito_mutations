#!/usr/bin/env bash
# Stages the files listed in data/zenodo_manifest.csv into a flat directory for
# upload to Zenodo (Zenodo stores files by name only, with no directories).
#
#   bash prepare_zenodo_upload.sh [output_dir]
#
# Then upload the contents of that directory to the Zenodo deposition, publish,
# and set ZENODO_RECORD_ID at the top of download_data.R to the published
# record id.
set -euo pipefail

OUT="${1:-zenodo_upload}"
MANIFEST="data/zenodo_manifest.csv"

[ -f "$MANIFEST" ] || { echo "Run from the repository root ($MANIFEST not found)" >&2; exit 1; }
mkdir -p "$OUT"

n=0
# Skip the header; fields are quoted CSV: filename,path,bytes,md5,description
while IFS= read -r line; do
  fname=$(printf '%s' "$line" | cut -d',' -f1 | tr -d '"')
  fpath=$(printf '%s' "$line" | cut -d',' -f2 | tr -d '"')
  [ -f "$fpath" ] || { echo "MISSING: $fpath" >&2; exit 1; }
  cp "$fpath" "$OUT/$fname"
  n=$((n+1))
done < <(tail -n +2 "$MANIFEST")

echo "Staged $n files into $OUT/ ($(du -sh "$OUT" | cut -f1))"
echo "Upload the contents of $OUT/ to Zenodo, publish, then set ZENODO_RECORD_ID in download_data.R"
