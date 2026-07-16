#!/usr/bin/env bash
# depgraph_receipt_e2e.sh (bead ongz) — the production dependency-receipt
# chain, end to end, with REAL invocations at every layer (no synthetic
# receipt authority anywhere):
#
#   xtask depgraph-receipt mints a locked single-root receipt
#     -> the nested Cargo invocation exports it to fs-la/build.rs
#     -> fs-session folds the resulting build identity into the GEMM
#        tune key and dispatch binding
#     -> fs-roofline retains the exact receipt artifact in the run ledger
#     -> a fresh process revalidates the retained evidence to Fresh
#     -> a warm process reuses the exact admitted tune row without resweep
#
# and then mutates every trust input, requiring a DIFFERENT build/tune key
# or a NAMED fail-closed refusal:
#
#   receipt-byte tamper, receipt/salt export omission, development salt
#   (never citable), unified feature drift, target drift, local path
#   package content drift, and cargo executable identity drift.
#
# The lane binds one before/after tree snapshot: the content-drift arm
# transiently adds ONE new file under crates/fs-la (removed by trap; the
# terminal seal proves the tree state is restored). Run this lane on BOTH
# reference ISA families (aarch64 + x86_64); every verdict row records the
# host ISA. The attested stages need a machine that passes the roofline
# plausibility floors — a saturated shared host is a named refusal, not a
# skip; rerun in a quiet window (see bead o582).
set -euo pipefail

cd "$(dirname "$0")/../.."

absolute_path() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s\n' "$PWD/$1" ;;
  esac
}

TOOLCHAIN="${FSIM_DEPGRAPH_E2E_TOOLCHAIN:-nightly-2026-07-06}"
CARGO_BIN="${FSIM_DEPGRAPH_E2E_CARGO:-$HOME/.cargo/bin/cargo}"
GEMM_N="${FSIM_DEPGRAPH_E2E_N:-4096}"
PROMOTE_RETRIES="${FSIM_DEPGRAPH_E2E_PROMOTE_RETRIES:-3}"
ROOT_PACKAGE="fs-roofline"
HOST_ISA=$(uname -m)
HEAD_SHA=$(git rev-parse HEAD)
LOG_ROOT=$(absolute_path "${FSIM_DEPGRAPH_E2E_LOG_DIR:-target/depgraph-receipt-e2e}")
mkdir -p "$LOG_ROOT"
LOG_DIR=$(mktemp -d "$LOG_ROOT/${HEAD_SHA:0:12}-XXXXXXXX")
TARGET_ROOT=$(absolute_path "${FSIM_DEPGRAPH_E2E_TARGET_DIR:-target/depgraph-receipt-e2e-builds}")
TARGET_DIR="$TARGET_ROOT/$(basename "$LOG_DIR")"
mkdir -p "$TARGET_DIR"
VERDICTS="$LOG_DIR/verdicts.jsonl"
: >"$VERDICTS"
CONTENT_PROBE="crates/fs-la/.depgraph-e2e-content-probe"
FAILURES=0

cleanup() { rm -f "$CONTENT_PROBE"; }
trap cleanup EXIT

tree_snapshot() {
  {
    git rev-parse HEAD
    git -c core.excludesFile=/dev/null status --porcelain --untracked-files=all
  } | shasum -a 256 | cut -d' ' -f1
}

row() { # lane status detail log
  python3 - "$1" "$2" "$3" "$4" "$HEAD_SHA" "$HOST_ISA" <<'PY' >>"$VERDICTS"
import json, sys
lane, status, detail, log, head, isa = sys.argv[1:]
print(json.dumps({
    "check": "depgraph-receipt-e2e",
    "lane": lane,
    "status": status,
    "head": head,
    "host_isa": isa,
    "detail": detail,
    "log": log,
}, separators=(",", ":")))
PY
  if [[ "$2" != "ok" ]]; then
    FAILURES=$((FAILURES + 1))
  fi
  printf '%-34s %-4s %s\n' "$1" "$2" "$3"
}

run_cargo() { # log-name args...
  local log="$LOG_DIR/$1.log"
  shift
  CARGO_TARGET_DIR="$TARGET_DIR" "$CARGO_BIN" "+$TOOLCHAIN" "$@" \
    >"$log" 2>&1
}

json_field() { # file key -> value of last occurrence at any nesting depth
  python3 - "$1" "$2" <<'PY'
import json, sys
path, key = sys.argv[1:]
value = None

def find(node):
    global value
    if isinstance(node, dict):
        if key in node:
            value = node[key]
        for child in node.values():
            find(child)
    elif isinstance(node, list):
        for child in node:
            find(child)

with open(path, encoding="utf-8") as handle:
    for line in handle:
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            find(json.loads(line))
        except json.JSONDecodeError:
            continue
print(json.dumps(value) if isinstance(value, (dict, list)) else value)
PY
}

SNAPSHOT_BEFORE=$(tree_snapshot)

# ---- S1: mint the locked single-root receipt (twice-derived internally),
# then prove verify-mode round-trips the exact bytes.
RECEIPT="$LOG_DIR/receipt.json"
if run_cargo s1-mint run -q -p xtask -- depgraph-receipt -- --package "$ROOT_PACKAGE" \
  && mv "$LOG_DIR/s1-mint.log" "$RECEIPT" && [[ -s "$RECEIPT" ]]; then
  row "s1-receipt-mint" ok "locked single-root receipt: $(wc -c <"$RECEIPT" | tr -d ' ') bytes" "$RECEIPT"
else
  row "s1-receipt-mint" fail "xtask depgraph-receipt refused" "$LOG_DIR/s1-mint.log"
  echo "cannot continue without a receipt" >&2
  exit 1
fi
if FRANKENSIM_DEPGRAPH_RECEIPT="$(cat "$RECEIPT")" \
   run_cargo s1-verify run -q -p xtask -- depgraph-receipt --verify -- --package "$ROOT_PACKAGE"; then
  row "s1-receipt-verify" ok "verify mode round-trips the exact receipt bytes" "$LOG_DIR/s1-verify.log"
else
  row "s1-receipt-verify" fail "verify refused the just-minted receipt" "$LOG_DIR/s1-verify.log"
fi

# ---- S2: the selected Cargo invocation supplies the receipt to fs-la
# build.rs (receipt-class build of the production binary + ceremony helper).
if FRANKENSIM_DEPGRAPH_RECEIPT="$(cat "$RECEIPT")" \
   run_cargo s2-build build -p fs-roofline --bin roofline --example attest_baseline; then
  row "s2-receipt-build" ok "roofline + attest_baseline built with the receipt exported" "$LOG_DIR/s2-build.log"
else
  row "s2-receipt-build" fail "receipt-class build refused" "$LOG_DIR/s2-build.log"
  echo "cannot continue without the receipt-class binary" >&2
  exit 1
fi
ROOFLINE="$TARGET_DIR/debug/roofline"
ATTEST="$TARGET_DIR/debug/examples/attest_baseline"

# ---- S3: operator ceremony (real machine probes; the plausibility floors
# are a hard gate) and the attested clean run reaching Fresh.
STORE="$LOG_DIR/plain_store.jsonl"
PROMOTED=false
for attempt in $(seq 1 "$PROMOTE_RETRIES"); do
  if "$ROOFLINE" promote --store "$STORE" --firmware depgraph-e2e-fw-v1 \
       --operator depgraph-e2e --justification "ongz production receipt e2e" \
       >"$LOG_DIR/s3-promote.log" 2>&1; then
    PROMOTED=true
    break
  fi
  sleep 15
done
if [[ "$PROMOTED" == true ]]; then
  row "s3-promote" ok "baseline promoted on attempt $attempt/$PROMOTE_RETRIES" "$LOG_DIR/s3-promote.log"
else
  row "s3-promote" fail "promotion refused $PROMOTE_RETRIES times (host too noisy for plausibility floors?)" "$LOG_DIR/s3-promote.log"
fi

LEDGER="$LOG_DIR/production.db"
RUN1="$LOG_DIR/s3-run1.jsonl"
if [[ "$PROMOTED" == true ]]; then
  FP=$(json_field "$LOG_DIR/s3-promote.log" fingerprint)
  if "$ATTEST" --store "$STORE" --fingerprint "$FP" \
       --key-id ops/depgraph-e2e --signature depgraph-e2e-signature-v1 \
       --out-attested "$LOG_DIR/attested.jsonl" \
       --out-authority "$LOG_DIR/authority.tsv" \
       --out-receipts "$LOG_DIR/receipts.txt" \
       >"$LOG_DIR/s3-attest.log" 2>&1; then
    row "s3-attest" ok "attested store + authority + retained receipts emitted" "$LOG_DIR/s3-attest.log"
  else
    row "s3-attest" fail "ceremony helper refused" "$LOG_DIR/s3-attest.log"
  fi
  : >"$LOG_DIR/dep_authority.txt" # empty canonical policy: no revocations
  attested_run() { # out-file
    "$ROOFLINE" --n "$GEMM_N" --warmup 1 --reps 2 --ledger "$LEDGER" \
      --baseline "$LOG_DIR/attested.jsonl" --firmware depgraph-e2e-fw-v1 \
      --authority-policy "$LOG_DIR/authority.tsv" \
      --retained-receipts "$LOG_DIR/receipts.txt" \
      --dependency-authority-policy "$LOG_DIR/dep_authority.txt" \
      >"$1" 2>&1
  }
  if attested_run "$RUN1" \
     && [[ "$(json_field "$RUN1" revalidated_fresh)" == "True" ]] \
     && [[ "$(json_field "$RUN1" citable)" == "True" ]]; then
    row "s3-fresh-run" ok "attested run recorded, revalidated Fresh, citable" "$RUN1"
  else
    row "s3-fresh-run" fail "attested run did not reach citable+Fresh" "$RUN1"
  fi

  # The exact receipt artifact is retained in the ledger with an In edge.
  EXPECTED_LEN=$(wc -c <"$RECEIPT" | tr -d ' ')
  RETAINED=$(sqlite3 "$LEDGER" \
    "SELECT count(*) FROM artifacts a JOIN edges e ON e.artifact = a.hash \
     WHERE a.kind = 'fs-la-depgraph-receipt' AND a.len = $EXPECTED_LEN AND e.role = 'in';" \
    2>>"$LOG_DIR/s3-ledger.log" || echo query-failed)
  if [[ "$RETAINED" =~ ^[1-9] ]]; then
    row "s3-retained-artifact" ok "ledger retains the exact $EXPECTED_LEN-byte receipt with an In edge" "$LEDGER"
  else
    row "s3-retained-artifact" fail "receipt artifact/edge missing (count=$RETAINED)" "$LEDGER"
  fi

  # ---- S4: a warm process reuses the exact admitted row without resweep.
  ROW1=$(json_field "$RUN1" tune_row_identity)
  TUNE_COUNT1=$(sqlite3 "$LEDGER" "SELECT count(*) FROM tune;" 2>/dev/null || echo query-failed)
  RUN2="$LOG_DIR/s4-run2.jsonl"
  if attested_run "$RUN2"; then
    ROW2=$(json_field "$RUN2" tune_row_identity)
    TUNE_COUNT2=$(sqlite3 "$LEDGER" "SELECT count(*) FROM tune;" 2>/dev/null || echo query-failed)
    if [[ -n "$ROW1" && "$ROW1" == "$ROW2" && "$TUNE_COUNT1" == "$TUNE_COUNT2" \
          && "$(json_field "$RUN2" citable)" == "True" ]]; then
      row "s4-warm-reuse" ok "warm dispatch reused admitted tune row $ROW1 (rows: $TUNE_COUNT2, no resweep)" "$RUN2"
    else
      row "s4-warm-reuse" fail "row1=$ROW1 row2=$ROW2 tune_rows=$TUNE_COUNT1->$TUNE_COUNT2" "$RUN2"
    fi
  else
    row "s4-warm-reuse" fail "warm attested run refused" "$RUN2"
  fi
else
  row "s3-fresh-run" fail "skipped: no promoted baseline" "$LOG_DIR/s3-promote.log"
  row "s3-retained-artifact" fail "skipped: no promoted baseline" "$LOG_DIR/s3-promote.log"
  row "s4-warm-reuse" fail "skipped: no promoted baseline" "$LOG_DIR/s3-promote.log"
fi

# ---- S5: mutations. Every trust input either moves the build/tune key or
# refuses with a named error.

# (a) receipt-byte tamper -> fs-la build script refuses the non-canonical
# receipt before any identity is minted.
TAMPERED=$(python3 - "$RECEIPT" <<'PY'
import sys
text = open(sys.argv[1], encoding="utf-8").read()
needle = '"schema":"fs-la-depgraph-receipt-v1"'
assert needle in text
print(text.replace(needle, '"schema":"fs-la-depgraph-receipt-v2"', 1), end="")
PY
)
if FRANKENSIM_DEPGRAPH_RECEIPT="$TAMPERED" \
   run_cargo s5-tamper check -p fs-la; then
  row "s5-receipt-tamper" fail "tampered receipt bytes were accepted by the build" "$LOG_DIR/s5-tamper.log"
elif grep -q "FRANKENSIM_DEPGRAPH_RECEIPT" "$LOG_DIR/s5-tamper.log"; then
  row "s5-receipt-tamper" ok "tampered receipt refused by fs-la build.rs (named)" "$LOG_DIR/s5-tamper.log"
else
  row "s5-receipt-tamper" fail "build failed without naming the receipt refusal" "$LOG_DIR/s5-tamper.log"
fi

# (b) export omission -> outside the workspace (no .cargo/config.toml salt)
# a build with neither receipt nor salt is a named refusal.
OMIT_DIR=$(mktemp -d "${TMPDIR:-/tmp}/depgraph-e2e-omit-XXXXXXXX")
cat >"$OMIT_DIR/Cargo.toml" <<TOML
[package]
name = "depgraph-e2e-omission-probe"
version = "0.0.0"
edition = "2021"

[dependencies]
fs-la = { path = "$PWD/crates/fs-la" }

[workspace]
TOML
mkdir -p "$OMIT_DIR/src" && echo "pub fn probe() {}" >"$OMIT_DIR/src/lib.rs"
if (cd "$OMIT_DIR" && env -u FRANKENSIM_DEPGRAPH_RECEIPT -u FRANKENSIM_DEPGRAPH_SALT \
     CARGO_TARGET_DIR="$TARGET_DIR/omission" RCH_DISABLE=1 \
     "$CARGO_BIN" "+$TOOLCHAIN" check --offline >"$LOG_DIR/s5-omission.log" 2>&1); then
  row "s5-export-omission" fail "build without receipt or salt was accepted" "$LOG_DIR/s5-omission.log"
elif grep -q "requires dependency-graph evidence" "$LOG_DIR/s5-omission.log"; then
  row "s5-export-omission" ok "omitted evidence refused by fs-la build.rs (named, bead fz2.6)" "$LOG_DIR/s5-omission.log"
else
  row "s5-export-omission" fail "build failed without naming the evidence refusal" "$LOG_DIR/s5-omission.log"
fi
rm -rf "$OMIT_DIR"

# (c) development salt builds, but can never become citable: the same
# attested ceremony against a salt-class binary must refuse freshness.
if [[ "$PROMOTED" == true ]] \
   && run_cargo s5-salt-build build -p fs-roofline --bin roofline; then
  SALT_LEDGER="$LOG_DIR/salt.db"
  SALT_RUN="$LOG_DIR/s5-salt-run.jsonl"
  if "$TARGET_DIR/debug/roofline" --n "$GEMM_N" --warmup 1 --reps 2 \
       --ledger "$SALT_LEDGER" --baseline "$LOG_DIR/attested.jsonl" \
       --firmware depgraph-e2e-fw-v1 --authority-policy "$LOG_DIR/authority.tsv" \
       --retained-receipts "$LOG_DIR/receipts.txt" \
       --dependency-authority-policy "$LOG_DIR/dep_authority.txt" \
       >"$SALT_RUN" 2>&1 \
     && [[ "$(json_field "$SALT_RUN" citable)" == "False" ]] \
     && grep -q "development equivalence salt" "$SALT_RUN"; then
    row "s5-salt-never-citable" ok "salt-class build refused citation (named)" "$SALT_RUN"
  else
    row "s5-salt-never-citable" fail "salt-class run did not refuse citation by name" "$SALT_RUN"
  fi
else
  row "s5-salt-never-citable" fail "skipped: no promoted baseline or salt build refused" "$LOG_DIR/s5-salt-build.log"
fi

verify_expect_mismatch() { # lane log-name detail extra-args...
  local lane="$1" log_name="$2" detail="$3"
  shift 3
  if FRANKENSIM_DEPGRAPH_RECEIPT="$(cat "$RECEIPT")" \
     run_cargo "$log_name" run -q -p xtask -- depgraph-receipt --verify -- \
       --package "$ROOT_PACKAGE" "$@"; then
    row "$lane" fail "verify accepted a drifted derivation" "$LOG_DIR/$log_name.log"
  elif grep -q "depgraph receipt mismatch" "$LOG_DIR/$log_name.log"; then
    row "$lane" ok "$detail" "$LOG_DIR/$log_name.log"
  else
    row "$lane" fail "verify failed without the named mismatch refusal" "$LOG_DIR/$log_name.log"
  fi
}

# (d) unified feature drift and (e) target drift move the derivation.
verify_expect_mismatch s5-feature-drift s5-features \
  "feature-selection drift is a named mismatch refusal" --all-features
if [[ "$HOST_ISA" == "arm64" || "$HOST_ISA" == "aarch64" ]]; then
  OTHER_TARGET="x86_64-unknown-linux-gnu"
else
  OTHER_TARGET="aarch64-apple-darwin"
fi
verify_expect_mismatch s5-target-drift s5-target \
  "target drift ($OTHER_TARGET) is a named mismatch refusal" --target "$OTHER_TARGET"

# (f) local path package content drift: one ADDED file inside the fs-la
# package moves the path digest (removed by trap; sealed below).
echo "depgraph-e2e content probe" >"$CONTENT_PROBE"
verify_expect_mismatch s5-content-drift s5-content \
  "local package content drift is a named mismatch refusal"
rm -f "$CONTENT_PROBE"

# (g) cargo executable identity drift: deriving under a different cargo
# changes the receipt's captured executable identity.
if STABLE_CARGO=$(rustup which cargo --toolchain stable 2>/dev/null) && [[ -x "$STABLE_CARGO" ]]; then
  if FRANKENSIM_DEPGRAPH_RECEIPT="$(cat "$RECEIPT")" CARGO="$STABLE_CARGO" \
     run_cargo s5-executable run -q -p xtask -- depgraph-receipt --verify -- \
       --package "$ROOT_PACKAGE"; then
    row "s5-executable-drift" fail "verify accepted a different cargo executable" "$LOG_DIR/s5-executable.log"
  elif grep -Eq "depgraph receipt mismatch|cargo executable" "$LOG_DIR/s5-executable.log"; then
    row "s5-executable-drift" ok "cargo executable drift is a named refusal" "$LOG_DIR/s5-executable.log"
  else
    row "s5-executable-drift" fail "verify failed without a named refusal" "$LOG_DIR/s5-executable.log"
  fi
else
  row "s5-executable-drift" fail "stable toolchain unavailable for executable-identity drift" "$LOG_DIR/s5-executable.log"
fi

# ---- S6: terminal seal — the tree snapshot is unchanged by the lane.
SNAPSHOT_AFTER=$(tree_snapshot)
if [[ "$SNAPSHOT_BEFORE" == "$SNAPSHOT_AFTER" ]]; then
  row "s6-tree-seal" ok "before/after tree snapshots match ($SNAPSHOT_BEFORE)" "$VERDICTS"
else
  row "s6-tree-seal" fail "tree moved during the lane: $SNAPSHOT_BEFORE -> $SNAPSHOT_AFTER" "$VERDICTS"
fi

echo
echo "verdicts: $VERDICTS"
if [[ "$FAILURES" -ne 0 ]]; then
  echo "depgraph-receipt-e2e: $FAILURES lane(s) failed" >&2
  exit 1
fi
echo "depgraph-receipt-e2e: all lanes green ($HOST_ISA)"
