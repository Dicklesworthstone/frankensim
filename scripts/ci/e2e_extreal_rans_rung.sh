#!/usr/bin/env bash
# Low-Re RANS rung e2e lane (beads frankensim-extreal-program-f85xj.5.8.*).
#
# Stage .5.8.1 (model-card freeze) is live: runs the freeze-integrity gate
# battery behind the rans-rung feature. Later stages (.5.8.2 solver,
# .5.8.3 validation campaign, .5.8.4 independent adjudication) append
# their modes here as they land; this script refuses unknown modes so a
# half-landed stage cannot masquerade as green.
#
# Usage:
#   scripts/ci/e2e_extreal_rans_rung.sh [--card]
set -u -o pipefail

MODE="${1:---card}"

case "$MODE" in
  --card)
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_card_gates || exit $?
    echo "rans-rung stage .5.8.1 (model-card freeze): PASS"
    ;;
  *)
    echo "e2e_extreal_rans_rung: unknown mode $MODE (available: --card)" >&2
    exit 30
    ;;
esac
