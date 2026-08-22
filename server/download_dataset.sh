#!/usr/bin/env bash
set -euo pipefail

ROOT=${A617_SERVER_ROOT:-$HOME/a617}
RAW="$ROOT/raw"
API=https://nrds.inl.gov/api/3/action/package_show
mkdir -p "$RAW"

download_package() {
  local package=$1
  curl -fsSL "$API?id=$package" |
    jq -r '.result.resources[]?.url | select(test("/download/"))' |
    while IFS= read -r url; do printf '%s\n' "$url"; done
}

urls=$(mktemp)
trap 'rm -f "$urls"' EXIT
download_package a617_test6-7_images_ebsd_spd_format >>"$urls"
download_package a617_test6-7_images_ebsd_dat_format >>"$urls"
download_package a617_test6-7_images_ebsd_ang_format_a617_test6-7_ebsd >>"$urls"
for element in al cl co cr f fe mg mn mo na ni o si ti w; do download_package "a617_test6-7_images_ebsd___eds_${element}" >>"$urls"; done

export RAW
download_one() {
  local url=$1 file="$RAW/${1##*/}" tmp="$RAW/.${1##*/}.part"
  [ -s "$file" ] || { curl -fL --retry 3 -o "$tmp" "$url" && mv -f "$tmp" "$file"; }
}
export -f download_one
xargs -P 4 -n 1 -I {} bash -c 'download_one "$1"' _ {} <"$urls"

printf 'raw files: '; find "$RAW" -type f | wc -l
