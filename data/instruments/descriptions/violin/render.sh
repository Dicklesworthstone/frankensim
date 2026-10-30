#!/usr/bin/env bash
# Four material strings coupled to the sourced il Cannone shell assembly and
# enclosed air, using the existing material_strings engine. The old .performance
# cards supply ONLY bow speed/force/articulation and duration; their one-way body,
# EI=0 string and constant modal damping are no longer the rendering route.
# usage: render.sh MATERIAL_STRINGS_BIN NEW_OUT_DIR [launcher options]
# See README.md for sources, limitations, material/gas variants and lab export.
set -euo pipefail
bin="${1:?usage: render.sh MATERIAL_STRINGS_BIN NEW_OUT_DIR [launcher options]}"
out="${2:?usage: render.sh MATERIAL_STRINGS_BIN NEW_OUT_DIR [launcher options]}"
shift 2
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo="$(cd "$here/../../../.." && pwd -P)"
materials="$repo/crates/fs-couple/examples/material_strings"
assembly_options="${ASSEMBLY_OPTIONS:-$materials/cannone-h8.assembly}"
performance="${MATERIAL_STRINGS_PERFORMANCE:-$materials/cannone-dominant-g3.material-strings}"
geometry="$repo/data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026"
for note in ${NOTES:-G3 D4 A4 E5}; do
    case "$note" in
        G3) string=0 ;; D4) string=1 ;; A4) string=2 ;; E5) string=3 ;;
        *) echo "unsupported open string: $note" >&2; exit 2 ;;
    esac
    python3 "$repo/apps/instrument-lab/tools/render_material_strings.py" \
        --renderer "$bin" --out "$out/$note" \
        --performance "$performance" \
        --pack-root "$materials/packs" \
        --shell-table "$geometry/shell-assembly-h8.tsv" \
        --assembly-options "$assembly_options" \
        --bow-card "$here/arco-open-$note.performance" --string-index "$string" \
        --bow-distance-from-bridge 0.1 --family violin --note "$note" \
        --observer "material shell-volume and cavity pressure at 1 m; compact receiver approximation" \
        --recipe data/instruments/descriptions/violin/render.sh "$@"
done
