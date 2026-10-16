#!/usr/bin/env bash
# Oracle determinism matrix leg (bead frankensim-uf7cw).
#
# Builds ONE leg of the G5 determinism matrix: an isolated snapshot tree of
# the committed frankensim source plus every constellation sibling frozen at
# either the LOCK PIN (`MODE=pinned`) or the checkout's live HEAD
# (`MODE=live`), then runs the gb_003 / gb_004 goldens plus the
# determinism_probe witness inside it. Nothing here mutates a shared
# checkout: snapshots are produced with read-only `git archive` extractions,
# so a drifted sibling can be measured without being moved.
#
# Outputs under OUT_DIR (created):
#   probe.jsonl        determinism probe payload (one header + case/error lines)
#   bore_goldens.log   golden-test stdout tail
#   leg_receipt.json   machine-readable provenance for the comparator
#
# Environment: FS_DET_ENV_LABEL overrides the emitted environment tag
# (default "<mode>-<arch>"). Exit codes: 0 green; 2 usage; 3 observation
# refusal (a pinned sibling object missing locally is NO-DATA for this leg,
# never silently skipped).
set -euo pipefail

usage() {
    echo "usage: MODE=pinned|live OUT_DIR=<dir> [$0]" >&2
    exit 2
}

MODE="${MODE:-}"
OUT_DIR="${OUT_DIR:-}"
[ -n "$MODE" ] && [ -n "$OUT_DIR" ] || usage
case "$MODE" in
    pinned | live) ;;
    *) usage ;;
esac

REPO="$(pwd)"
if [ ! -f "$REPO/constellation.lock" ] || [ ! -f "$REPO/Cargo.toml" ]; then
    echo "refusal: run from the frankensim repository root" >&2
    exit 3
fi

SNAPSHOT="$OUT_DIR/snapshot"
if [ -e "$SNAPSHOT" ]; then
    echo "refusal: $SNAPSHOT already exists; use a fresh OUT_DIR per leg" >&2
    exit 3
fi
mkdir -p "$OUT_DIR" "$SNAPSHOT"

# ---- frankensim: committed source only (immune to shared-tree churn) -----
git -C "$REPO" archive --prefix=frankensim/ HEAD | tar -x -C "$OUT_DIR"

# ---- siblings: archive rows out of constellation.lock --------------------
# Rows are single-line JSON objects; collapse spacing so a plain sed pair
# extraction works without jq/python assumptions on workers.
LOCK_ROWS=$(tr -d ' \t\r' < "$REPO/constellation.lock" |
    sed -n 's/^.*"lib":"\([^"]*\)","version":"[^"]*","git_head":"\([^"]*\)".*$/\1 \2/p')

CONSTELLATION_ROOT="$(cd "$REPO/.." && pwd)"
SIBLING_RECEIPTS=""
while read -r lib pin; do
    [ -n "$lib" ] || continue
    CHECKOUT="$CONSTELLATION_ROOT/$lib"
    if [ ! -d "$CHECKOUT/.git" ]; then
        echo "{\"record\":\"refusal\",\"reason\":\"sibling checkout missing\",\"lib\":\"$lib\"}" >&2
        exit 3
    fi
    if [ "$MODE" = "pinned" ]; then
        REF="$pin"
        if ! git -C "$CHECKOUT" cat-file -e "$REF^{commit}" 2>/dev/null; then
            echo "{\"record\":\"refusal\",\"reason\":\"pinned object absent locally\",\"lib\":\"$lib\",\"pin\":\"$pin\"}" >&2
            exit 3
        fi
    else
        REF="$(git -C "$CHECKOUT" rev-parse HEAD)"
    fi
    git -C "$CHECKOUT" archive --prefix="$lib/" "$REF" | tar -x -C "$OUT_DIR"
    SIBLING_RECEIPTS="$SIBLING_RECEIPTS,\"$lib\":\"$REF\""
done <<ROWS
$LOCK_ROWS
ROWS
SIBLING_RECEIPTS="${SIBLING_RECEIPTS#,}"

# ---- build + run inside the frozen tree ----------------------------------
export CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/oracle-det-target"
SNAP_FS="$SNAPSHOT/frankensim"
cd "$SNAP_FS"

GOLDENS_PASS=true
cargo test -p fs-query --test bore -- gb_003 gb_004 2>&1 |
    tee "$OUT_DIR/bore_goldens.log" || GOLDENS_PASS=false

LABEL="${FS_DET_ENV_LABEL:-$MODE-$(uname -m)}"
FS_DET_ENV="$LABEL" cargo run -q -p fs-query --example determinism_probe \
    >"$OUT_DIR/probe.jsonl"

cat >"$OUT_DIR/leg_receipt.json" <<RECEIPT
{"record":"leg_receipt","schema":"frankensim-oracle-leg-v1","mode":"$MODE",
 "arch":"$(uname -m)","os":"$(uname -s)","env_label":"$LABEL",
 "frankensim_head":"$FS_HEAD","goldens_pass":$GOLDENS_PASS,
 "siblings":{$SIBLING_RECEIPTS}}
RECEIPT
if [ "${FS_DET_STDOUT:-0}" = "1" ]; then
    echo "=====FS_DET_BEGIN leg_receipt.json====="
    cat "$OUT_DIR/leg_receipt.json"
    echo "=====FS_DET_END====="
    echo "=====FS_DET_BEGIN probe.jsonl====="
    cat "$OUT_DIR/probe.jsonl"
    echo "=====FS_DET_END====="
fi

"$GOLDENS_PASS" || {
    echo "leg failed: goldens red inside frozen tree; receipts emitted above" >&2
    exit 1
}
echo "leg complete: $OUT_DIR"
"$GOLDENS_PASS" || {
    echo "leg failed: goldens red inside frozen tree; receipts retained" >&2
    exit 1
}
echo "leg complete: $OUT_DIR"

