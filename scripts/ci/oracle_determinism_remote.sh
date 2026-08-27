#!/usr/bin/env bash
# Remote EXECUTE half of an oracle determinism matrix leg (bead
# frankensim-uf7cw). Runs on a worker whose /data/tmp/det-<label>/snapshot
# was assembled and streamed from the AUTHORITATIVE checkout (mirrors carry
# no .git, so archive-based freezing must happen on the origin host).
#
# Contract (env):
#   MODE            pinned|live          (receipt metadata)
#   LABEL           unique matrix label
#   DET_DIR         directory containing snapshot/ (receives outputs too)
#   FS_DET_STDOUT   1 => emit receipt + payload blocks on stdout
#
# The EXIT trap self-reports the terminal status, so any transport that can
# pipe stdout gets the verdict regardless of its own quoting behavior.
# Exit codes: 2 usage; 3 refusal; 1 goldens red (receipts still emitted).
set -euo pipefail

MODE="${MODE:-}"
LABEL="${LABEL:-}"
DET_DIR="${DET_DIR:-}"
[ -n "$MODE" ] && [ -n "$LABEL" ] && [ -n "$DET_DIR" ] || {
    echo "usage: MODE=pinned|live LABEL=<label> DET_DIR=<dir> [$0]" >&2
    exit 2
}

trap 'echo "LEG_RC=$?"' EXIT

SNAPSHOT="$DET_DIR/snapshot"
FS_TREE="$SNAPSHOT/frankensim"
[ -f "$FS_TREE/Cargo.toml" ] || {
    echo "refusal: snapshot missing at $SNAPSHOT" >&2
    exit 3
}

export CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/oracle-det-target"

cd "$FS_TREE"

GOLDENS_PASS=true
cargo test -p fs-query --test bore -- gb_003 gb_004 2>&1 |
    tail -200 >"$DET_DIR/bore_goldens.log" || GOLDENS_PASS=false

FS_DET_ENV="$LABEL" cargo run -q -p fs-query --example determinism_probe \
    >"$DET_DIR/probe.jsonl"

cat >"$DET_DIR/leg_receipt.json" <<RECEIPT
{"record":"leg_receipt","schema":"frankensim-oracle-leg-v1","mode":"$MODE",
 "arch":"$(uname -m)","os":"$(uname -s)","env_label":"$LABEL",
 "goldens_pass":$GOLDENS_PASS,
 "siblings_note":"snapshotted on origin host before transfer"}
RECEIPT

if [ "${FS_DET_STDOUT:-0}" = "1" ]; then
    echo "=====FS_DET_BEGIN leg_receipt.json====="
    cat "$DET_DIR/leg_receipt.json"
    echo "=====FS_DET_END====="
    echo "=====FS_DET_BEGIN probe.jsonl====="
    cat "$DET_DIR/probe.jsonl"
    echo "=====FS_DET_END====="
fi

if [ "$GOLDENS_PASS" != "true" ]; then
    echo "leg failed: goldens red inside frozen tree; receipts emitted above" >&2
    exit 1
fi
echo "leg complete: $DET_DIR"
