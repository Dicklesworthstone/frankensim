#!/usr/bin/env bash
# Cooling fidelity graph instance e2e lane (bead frankensim-extreal-program-f85xj.10.4).
#
# Runs the cooling fidelity-graph instance battery (7 nodes, 4 evidenced edges,
# 3 explicit gaps, 3 flagship demonstrations of cost!=authority, regime-dependent
# routing, and honest gap demand) and verifies structural determinism.
#
# Usage:
#   scripts/ci/e2e_cooling_fidelity_graph.sh [--all | --demos]
set -eu -o pipefail

MODE="${1:---all}"

case "$MODE" in
  --all|--demos)
    cargo test -q -p fs-plan --features cooling-instance --test cooling_instance || exit $?
    echo '{"suite":"cooling-fidelity-graph","case":"instance-and-demos","nodes":7,"evidenced_edges":4,"gaps":3,"verdict":"PASS"}'
    ;;
  *)
    echo "e2e_cooling_fidelity_graph: unknown mode $MODE (available: --all, --demos)" >&2
    exit 30
    ;;
esac
