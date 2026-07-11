#!/usr/bin/env bash
# x86_runtime_sweep.sh (bead yuyy) — the RUNTIME-class x86 lane.
#
# The ebro cross-check firewall catches what fails to COMPILE for
# x86_64; this lane catches what compiles everywhere but RUNS
# differently on x86 — per-ISA schedule engagement (the obq0 fs-cheb
# golden divergence), storage/page-size behavior (u8og), capsule
# dispatch, fsqlite MVCC differences. It picks the first REACHABLE and
# QUIET Threadripper, fast-forwards its clone to origin/main, and runs
# the full workspace test suite there (--no-fail-fast), reporting one
# quality-lane JSONL verdict bound to the tested HEAD and machine.
#
# Operational lessons encoded (from the hand runs this lane replaces):
# the Threadrippers double as rch workers whose background lanes can
# kill long test runs and saturate memory bandwidth, so each host is
# load-gated with bounded retries, and the remote run is wrapped in
# setsid so an ssh drop cannot orphan-kill it. A dirty or diverged
# clone is a named diagnostic, never silently tested or modified. Exactly one
# terminal quality-lane verdict is emitted for the overall run.
#
# Policy (owner decision recorded in the bead): when no quiet host is
# reachable this emits status=skipped with the reason — the nightly
# flow treats staleness as its own alarm; this lane is not a required
# pre-commit gate.
set -euo pipefail

cd "$(dirname "$0")/../.."

HOSTS="${FSIM_RUNTIME_SWEEP_HOSTS:-ts1:/home/ubuntu/frankensim ts2:/home/ubuntu/frankensim}"
LOAD_MAX="${FSIM_RUNTIME_SWEEP_LOAD_MAX:-12}"
RETRIES="${FSIM_RUNTIME_SWEEP_RETRIES:-3}"
RETRY_SLEEP="${FSIM_RUNTIME_SWEEP_RETRY_SLEEP:-120}"
LOG_DIR="${FSIM_RUNTIME_SWEEP_LOG_DIR:-$PWD/target/x86-runtime-logs}"
mkdir -p "$LOG_DIR"
RUN_ID="$(date +%s)-$$"
LOCAL_HEAD=$(git rev-parse HEAD)

repo_dirt() { # repository
  local repo="$1" tracked untracked index_flags hidden_index
  tracked=$(git -C "$repo" -c core.fileMode=true -c core.excludesFile=/dev/null \
    status --porcelain=v1 --untracked-files=all) || return 1
  untracked=$(git -C "$repo" -c core.excludesFile=/dev/null \
    ls-files --others --exclude-per-directory=.gitignore) || return 1
  index_flags=$(git -C "$repo" ls-files -v) || return 1
  hidden_index=$(printf '%s\n' "$index_flags" \
    | LC_ALL=C awk 'substr($0, 1, 1) == "S" || substr($0, 1, 1) ~ /[a-z]/ { print; exit }')
  printf '%s%s%s' "$tracked" "$untracked" \
    "${hidden_index:+index flag hides worktree state: $hidden_index}"
}

row() { # status host head detail log
  python3 - "$1" "$2" "$3" "$4" "$5" <<'PY'
import json
import sys

status, host, head, detail, log = sys.argv[1:]
print(json.dumps({
    "check": "quality-lane",
    "lane": "x86-runtime-workspace-tests",
    "status": status,
    "host": host,
    "head": head,
    "detail": detail,
    "log": log,
}, separators=(",", ":")))
PY
}

if ! local_dirt_before=$(repo_dirt .); then
  row "fail" "none" "$LOCAL_HEAD" "could not inspect local working-tree state" "-"
  exit 1
fi
if [ -n "$local_dirt_before" ]; then
  row "fail" "none" "$LOCAL_HEAD" \
    "local working tree is dirty; a remote origin/main run cannot prove this source state" "-"
  exit 1
fi

for attempt in $(seq 1 "$RETRIES"); do
  for pair in $HOSTS; do
    host="${pair%%:*}"
    clone="${pair#*:}"
    # Reachable?
    if ! ssh -o ConnectTimeout=10 -o BatchMode=yes "$host" true 2>/dev/null; then
      continue
    fi
    # Quiet? (1-minute load; the rch worker traffic this gates against
    # barely moves CPU load, so the in-suite roofline guards remain the
    # deeper filter — this check just avoids obviously doomed windows.)
    load=$(ssh -o ConnectTimeout=10 "$host" "cut -d' ' -f1 /proc/loadavg" 2>/dev/null || echo 999)
    if ! awk -v l="$load" -v m="$LOAD_MAX" 'BEGIN{exit !(l<m)}'; then
      continue
    fi
    # Clone present, clean, and fast-forwardable to origin/main.
    if ! state=$(ssh -o ConnectTimeout=10 "$host" bash -s -- "$clone" 2>/dev/null <<'REMOTE_STATE'
clone=$1
cd "$clone" 2>/dev/null || { echo no-clone; exit 0; }
repo_dirt() {
  tracked=$(git -c core.fileMode=true -c core.excludesFile=/dev/null \
    status --porcelain=v1 --untracked-files=all) || return 1
  untracked=$(git -c core.excludesFile=/dev/null \
    ls-files --others --exclude-per-directory=.gitignore) || return 1
  index_flags=$(git ls-files -v) || return 1
  hidden_index=$(printf '%s\n' "$index_flags" \
    | LC_ALL=C awk 'substr($0, 1, 1) == "S" || substr($0, 1, 1) ~ /[a-z]/ { print; exit }')
  printf '%s%s%s' "$tracked" "$untracked" \
    "${hidden_index:+index flag hides worktree state: $hidden_index}"
}
if ! dirt=$(repo_dirt); then
  echo inspect-failed
  exit 0
fi
if test -n "$dirt"; then
  echo dirty
  exit 0
fi
git fetch -q origin main 2>/dev/null || { echo fetch-failed; exit 0; }
git merge-base --is-ancestor HEAD origin/main || { echo diverged; exit 0; }
git merge --ff-only -q origin/main 2>/dev/null || { echo ff-failed; exit 0; }
if ! dirt=$(repo_dirt); then
  echo inspect-after-ff-failed
  exit 0
fi
if test -n "$dirt"; then
  echo dirty-after-ff
  exit 0
fi
echo ok
REMOTE_STATE
); then
      printf 'x86 runtime sweep admission probe failed for %s:%s\n' \
        "$host" "$clone" >&2
      continue
    fi
    state=$(printf '%s\n' "$state" | tail -1)
    if [ "$state" != "ok" ]; then
      printf 'x86 runtime sweep rejected %s:%s: %s\n' "$host" "$clone" "$state" >&2
      continue
    fi
    if ! head=$(ssh -o ConnectTimeout=10 "$host" bash -s -- "$clone" <<'REMOTE_HEAD'
cd "$1" && git rev-parse HEAD
REMOTE_HEAD
); then
      printf 'x86 runtime sweep could not capture admitted HEAD from %s:%s\n' \
        "$host" "$clone" >&2
      continue
    fi
    if [ "$head" != "$LOCAL_HEAD" ]; then
      row "fail" "$host" "$head" \
        "remote origin/main head does not equal admitted local head $LOCAL_HEAD" "-"
      exit 1
    fi
    log="$LOG_DIR/${head:0:12}-${host}.log"
    remote_log="/tmp/frankensim-x86-runtime-${head:0:12}-${RUN_ID}.log"
    remote_status=0
    ssh -o ConnectTimeout=10 "$host" bash -s -- "$clone" "$head" "$remote_log" <<'REMOTE_RUN' \
      || remote_status=$?
clone=$1
expected_head=$2
log=$3
cd "$clone" || exit 120
repo_dirt() {
  tracked=$(git -c core.fileMode=true -c core.excludesFile=/dev/null \
    status --porcelain=v1 --untracked-files=all) || return 1
  untracked=$(git -c core.excludesFile=/dev/null \
    ls-files --others --exclude-per-directory=.gitignore) || return 1
  index_flags=$(git ls-files -v) || return 1
  hidden_index=$(printf '%s\n' "$index_flags" \
    | LC_ALL=C awk 'substr($0, 1, 1) == "S" || substr($0, 1, 1) ~ /[a-z]/ { print; exit }')
  printf '%s%s%s' "$tracked" "$untracked" \
    "${hidden_index:+index flag hides worktree state: $hidden_index}"
}
: >"$log" || exit 120
test "$(git rev-parse HEAD)" = "$expected_head" || exit 121
initial_worktree=$(repo_dirt) || exit 122
if test -n "$initial_worktree"; then
  printf 'FSIM_WORKTREE_BEFORE_DIAGNOSTIC_BEGIN\n%s\nFSIM_WORKTREE_BEFORE_DIAGNOSTIC_END\n' \
    "$initial_worktree" >>"$log"
  exit 122
fi
set +e
snapshot_before_output=$(scripts/ci/checkout_constellation.sh --snapshot 2>&1)
snapshot_before_status=$?
set -e
printf 'FSIM_SNAPSHOT_BEFORE_OUTPUT_BEGIN\n%s\nFSIM_SNAPSHOT_BEFORE_OUTPUT_END\n' \
  "$snapshot_before_output" >>"$log"
snapshot_before=$(printf '%s\n' "$snapshot_before_output" | tail -1)
printf 'FSIM_SNAPSHOT_BEFORE_STATUS=%s\nFSIM_SNAPSHOT_BEFORE=%s\n' \
  "$snapshot_before_status" "$snapshot_before" >>"$log"
if test "$snapshot_before_status" -ne 0 \
    || ! printf '%s\n' "$snapshot_before" | grep -Eq '^[0-9a-f]{64}$'; then
  exit 126
fi
set +e
# HERMETIC target dir: host-global CARGO_TARGET_DIR values point at
# shared caches where other toolchains' artifacts cause E0514
# cross-rustc contamination (seen live on ts2) and concurrent-agent
# artifact races. The sweep owns its own warm dir inside the clone.
setsid env PATH="$HOME/.cargo/bin:$PATH" \
  CARGO_TARGET_DIR="$clone/target-x86-runtime-sweep" \
  cargo test --locked --workspace --no-fail-fast >>"$log" 2>&1
status=$?
set -e
final_head=$(git rev-parse HEAD) || exit 123
final_worktree=$(repo_dirt) || exit 124
set +e
snapshot_after_output=$(scripts/ci/checkout_constellation.sh --snapshot 2>&1)
snapshot_after_status=$?
set -e
printf 'FSIM_SNAPSHOT_AFTER_OUTPUT_BEGIN\n%s\nFSIM_SNAPSHOT_AFTER_OUTPUT_END\n' \
  "$snapshot_after_output" >>"$log"
snapshot_after=$(printf '%s\n' "$snapshot_after_output" | tail -1)
if test -n "$final_worktree"; then
  worktree_state=dirty
else
  worktree_state=clean
fi
printf '\nFSIM_CARGO_EXIT=%s\nFSIM_HEAD_BEFORE=%s\nFSIM_HEAD_AFTER=%s\nFSIM_WORKTREE_AFTER=%s\nFSIM_SNAPSHOT_AFTER_STATUS=%s\nFSIM_SNAPSHOT_AFTER=%s\n' \
  "$status" "$expected_head" "$final_head" "$worktree_state" \
  "$snapshot_after_status" "$snapshot_after" >>"$log"
if test -n "$final_worktree"; then
  printf 'FSIM_WORKTREE_PORCELAIN_BEGIN\n%s\nFSIM_WORKTREE_PORCELAIN_END\n' \
    "$final_worktree" >>"$log"
fi
test "$final_head" = "$expected_head" || exit 124
test "$worktree_state" = clean || exit 125
if test "$snapshot_after_status" -ne 0 \
    || ! printf '%s\n' "$snapshot_after" | grep -Eq '^[0-9a-f]{64}$'; then
  exit 126
fi
test "$snapshot_after" = "$snapshot_before" || exit 127
exit "$status"
REMOTE_RUN
    if ! scp -q "$host:$remote_log" "$log"; then
      row "fail" "$host" "$head" \
        "could not retain the remote runtime log (ssh status ${remote_status})" "$log"
      exit 1
    fi
    ok=$(grep -cE 'test result: ok' "$log" || true)
    failed=$(grep -cE 'test result: FAILED' "$log" || true)
    recorded_status=$(sed -n 's/^FSIM_CARGO_EXIT=//p' "$log" | tail -1)
    before=$(sed -n 's/^FSIM_HEAD_BEFORE=//p' "$log" | tail -1)
    after=$(sed -n 's/^FSIM_HEAD_AFTER=//p' "$log" | tail -1)
    worktree=$(sed -n 's/^FSIM_WORKTREE_AFTER=//p' "$log" | tail -1)
    snapshot_before_status=$(sed -n 's/^FSIM_SNAPSHOT_BEFORE_STATUS=//p' "$log" | tail -1)
    snapshot_before=$(sed -n 's/^FSIM_SNAPSHOT_BEFORE=//p' "$log" | tail -1)
    snapshot_after_status=$(sed -n 's/^FSIM_SNAPSHOT_AFTER_STATUS=//p' "$log" | tail -1)
    snapshot_after=$(sed -n 's/^FSIM_SNAPSHOT_AFTER=//p' "$log" | tail -1)
    local_after=$(git rev-parse HEAD)
    if ! local_dirt_after=$(repo_dirt .); then
      local_dirt_after="<inspection-failed>"
    fi
    if [ "$remote_status" -eq 0 ] \
        && [ "$recorded_status" = "0" ] \
        && [ "$before" = "$head" ] \
        && [ "$after" = "$head" ] \
        && [ "$worktree" = "clean" ] \
        && [ "$snapshot_before_status" = "0" ] \
        && [ "$snapshot_after_status" = "0" ] \
        && [ -n "$snapshot_before" ] \
        && [ "$snapshot_before" = "$snapshot_after" ] \
        && [ "$local_after" = "$LOCAL_HEAD" ] \
        && [ -z "$local_dirt_after" ] \
        && [ "$failed" -eq 0 ] \
        && [ "$ok" -gt 0 ]; then
      row "pass" "$host" "$head" "workspace runtime tests green on x86: ${ok} suites" "$log"
      exit 0
    fi
    row "fail" "$host" "$head" \
      "x86 runtime command status=${remote_status}/${recorded_status:-missing}, heads=${before:-missing}/${after:-missing}, worktree=${worktree:-missing}, snapshots=${snapshot_before_status:-missing}/${snapshot_after_status:-missing}:${snapshot_before:-missing}/${snapshot_after:-missing}, local=${LOCAL_HEAD}/${local_after}:${local_dirt_after:+dirty}, suites=${failed} red/${ok} green" "$log"
    exit 1
  done
  if [ "$attempt" -lt "$RETRIES" ]; then
    printf 'x86 runtime sweep attempt %s/%s found no admissible host; retrying after %ss\n' \
      "$attempt" "$RETRIES" "$RETRY_SLEEP" >&2
    sleep "$RETRY_SLEEP"
  fi
done

row "skipped" "none" "unknown" "no quiet reachable Threadripper within ${RETRIES} attempts (load gate ${LOAD_MAX}); staleness is the nightly flow's alarm" "-"
exit 0
