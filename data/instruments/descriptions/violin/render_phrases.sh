#!/usr/bin/env bash
# The violin phrases of the instrument-tests site (the G major scale up and down, the first phrase of Ode to Joy at
# 100 bpm, four double stops in G major) on the route of render.sh's open notes.
# gen_violin_phrase_notes.py writes one note per rendered note from the same h6 card (its tensions, its
# residual-flexibility record, and the bow cards' force and speed), each with its own bow schedule and a 0.45 s ring
# after the bow lifts. A stopped note adds the finger as a preloaded spring and dashpot pair (--stop 1e6 1.9:
# K = 1e6 N/m, a stated unknown; C = 1.9 N s/m). Every note renders through render.sh's launcher, packs, shell table
# and assembly options, with its stop profile and no --bow-card. One core per note.
# apps/instrument-explorer/tools/phrases/phrasekit.py joins the renders (a note stops its own string).
# usage: render_phrases.sh MATERIAL_STRINGS_BIN NEW_OUT_DIR [NOTES_AT_ONCE]
#   NEW_OUT_DIR/notes/violin/NAME.material-strings and NAME.stops   the 21 written notes
#   NEW_OUT_DIR/violin/NAME/render.wav                              each note's render, 144 kHz
# NOTES_AT_ONCE (default 1) notes render in parallel; each takes one core and about 1.5 GB.
set -euo pipefail
usage="usage: render_phrases.sh MATERIAL_STRINGS_BIN NEW_OUT_DIR [NOTES_AT_ONCE]"
bin="${1:?$usage}"
out="${2:?$usage}"
parallel="${3:-1}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo="$(cd "$here/../../../.." && pwd -P)"
materials="$repo/crates/fs-couple/examples/material_strings"
performance="${MATERIAL_STRINGS_PERFORMANCE:-$materials/cannone-dominant-g3-h6.material-strings}"
assembly_options="${ASSEMBLY_OPTIONS:-$materials/cannone-h6.assembly}"
geometry="$repo/data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026"
shell_table="${SHELL_TABLE:-$geometry/shell-assembly-h6.tsv}"
mkdir "$out"
out="$(cd "$out" && pwd -P)"
python3 -I "$materials/gen_violin_phrase_notes.py" "$repo" "$performance" "$out/notes" --stop 1e6 1.9
mkdir "$out/violin"

render_note() {
    local note="$1" name stops
    name="$(basename "$note" .material-strings)"
    stops="${note%.material-strings}.stops"
    set -- --renderer "$bin" --out "$out/violin/$name" --performance "$note" \
        --pack-root "$materials/packs" --shell-table "$shell_table" --assembly-options "$assembly_options"
    if [ -e "$stops" ]; then
        set -- "$@" --stop-profile "$stops"
    fi
    python3 "$repo/apps/instrument-lab/tools/render_material_strings.py" "$@" --family violin --note "$name" \
        --observer "material shell-volume and cavity pressure at 1 m; compact receiver approximation" \
        --recipe data/instruments/descriptions/violin/render_phrases.sh
}
export -f render_note
export bin out repo materials assembly_options shell_table
find "$out/notes/violin" -name '*.material-strings' | sort | xargs -P "$parallel" -I{} bash -c 'render_note "$1"' _ {}
