#!/usr/bin/env bash
# Gt-S1's plain nylon treble strings, one at a time, on the assembled body of the 1964 Esteso guitar, using the existing
# material_strings engine (pexfp.377). G3, B3 and E4 are open; F#4 (Gb4) and F5 are the B string stopped at frets 7
# and 18. The launcher assembles each string's body profile from the bridged h12 shell table and that string's assembly
# options (esteso-h12-{g,b,e}.assembly), then renders the card. F5's card runs at 96 kHz. The wound E2 takes Woodhouse
# 2004's whole-string coefficients through the measured section (wound-basses.tsv, its own assembly options
# esteso-h12-e2.assembly); the wound strings' core and winding construction stays UNKNOWN (see sources.tsv).
# The flat cedar-plate cards (gt-s1-*) stay in this directory; their recipe is this file's history.
# usage: render.sh MATERIAL_STRINGS_BIN XTASK_BIN NEW_OUT_DIR [launcher options]
# Reference notes in the Iowa library (48 kHz mono, mf): G3 Guitar.mf.sulG.G3B3, B3 Guitar.mf.sulB.B3,
# E4 Guitar.mf.sul_E.E4B4, F#4 Guitar.mf.sulB.C4B4, F5 Guitar.mf.sulB.C5Gb5, E2 Guitar.mf.sulE.E2B2 (sources.tsv).
set -euo pipefail
bin="${1:?usage: render.sh MATERIAL_STRINGS_BIN XTASK_BIN NEW_OUT_DIR [launcher options]}"
xtask="${2:?usage: render.sh MATERIAL_STRINGS_BIN XTASK_BIN NEW_OUT_DIR [launcher options]}"
out="${3:?usage: render.sh MATERIAL_STRINGS_BIN XTASK_BIN NEW_OUT_DIR [launcher options]}"
shift 3
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo="$(cd "$here/../../../.." && pwd -P)"
materials="$repo/crates/fs-couple/examples/material_strings"
geometry="$repo/data/instruments/geometry/guitar-esteso-1964-zimmer2024"
mkdir "$out" "$out/packs"
for pack in nylon-g3 nylon-b3 nylon-e4 wound-basses; do
    cp "$here/$pack.tsv" "$here/$pack-manifest.tsv" "$out/packs/"
    (cd "$repo" && "$xtask" matdb-pack --manifest "$out/packs/$pack-manifest.tsv" --out "$out/packs/$pack.fsmatpk")
done
for note in ${NOTES:-G3 B3 E4 Gb4 F5 E2}; do
    case "$note" in
        G3) card=esteso-g3; string=g; pack=nylon-g3 ;;
        B3) card=esteso-b3; string=b; pack=nylon-b3 ;;
        E4) card=esteso-e4; string=e; pack=nylon-e4 ;;
        Gb4) card=esteso-gb4; string=b; pack=nylon-b3 ;;
        F5) card=esteso-f5; string=b; pack=nylon-b3 ;;
        E2) card=esteso-e2; string=e2; pack=wound-basses ;;
        *) echo "unsupported note: $note" >&2; exit 2 ;;
    esac
    python3 "$repo/apps/instrument-lab/tools/render_material_strings.py" \
        --renderer "$bin" --out "$out/$note" \
        --performance "$here/$card.material-strings" \
        --pack-root "$materials/packs" --pack "$pack.fsmatpk=$out/packs/$pack.fsmatpk" \
        --shell-table "$geometry/shell-assembly-h12-bridge.tsv" \
        --assembly-options "$here/esteso-h12-$string.assembly" \
        --family guitar --note "$note" \
        --observer "material shell-volume and cavity pressure at 1 m; compact receiver approximation" \
        --recipe data/instruments/descriptions/guitar/render.sh "$@"
done
