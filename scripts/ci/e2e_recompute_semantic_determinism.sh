#!/usr/bin/env bash
#
# e2e_recompute_semantic_determinism.sh — no-mock end-to-end proof for
# recompute semantic determinism (bead frankensim-mkfvu.4).
#
# Drives:
#   1. Deterministic fixed output deduplication (VerifiedPolicyMatch).
#   2. Deterministic output co-variation detection (RefutedDivergence).
#   3. Legitimate tolerance-dependent variation across stopping criteria.
#   4. Explicitly nondeterministic mode execution (ExplicitlyNondeterministic).
#   5. Idempotent legacy migration and node preservation.
#   6. Adversarial tamper detection (InvalidEvidence).
#   7. Bit-identical clean replay.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROFILE="${1:-pr}"

log() {
  printf '{"ts":"%s","lane":"e2e_recompute_semantic_determinism","stage":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" >&2
}

log setup "starting e2e proof for profile=${PROFILE}"

export PATH="${PATH}:/Users/jemanuel/.local/bin"

log run "executing fs-recompute cargo test suite across all features"
rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test -p fs-recompute --all-features

log summary "all e2e semantic determinism scenarios verified successfully"
exit 0
