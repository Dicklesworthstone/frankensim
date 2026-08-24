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

MODE="${1:---all}"

case "$MODE" in
  --card)
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_card_gates || exit $?
    echo "rans-rung stage .5.8.1 (model-card freeze): PASS"
    ;;
  --solver)
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_solver || exit $?
    echo "rans-rung stage .5.8.2 (solver rung): PASS"
    ;;
  --validation)
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_validation || exit $?
    echo "rans-rung stage .5.8.3 (validation matrix): PASS"
    ;;
  --adjudicate)
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_adjudication || exit $?
    echo "rans-rung stage .5.8.4 (independent adjudication): PASS"
    ;;
  --all)
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_card_gates || exit $?
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_solver || exit $?
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_validation || exit $?
    cargo test -q -p fs-scenario --features rans-rung \
      --test rans_adjudication || exit $?

    REPORT_DIR="target/reports"
    mkdir -p "$REPORT_DIR"
    RECEIPT_FILE="${REPORT_DIR}/extreal_rans_rung_receipt.json"

    cat > "${RECEIPT_FILE}" << 'EOF'
{
  "receipt_version": 1,
  "domain": "org.frankensim.fs-scenario.rans-adjudication.v1",
  "node": {
    "node_id": "e10-low-re-rans",
    "governing_regime": "steady forced convection, attached / channel / fin array flow",
    "cost_tier": "Moderate (O(N_cells))",
    "authority_class": "Estimate"
  },
  "edge": {
    "source_node": "e10-low-re-rans",
    "target_qoi": "temperature_and_thermal_resistance",
    "evidence_tier": "ContextualValidatedEstimate",
    "admitted_contexts": [
      "attached_channel_flow",
      "heatsink_fin_array",
      "mild_buoyancy_forced_convection"
    ],
    "refused_contexts": [
      "massive_unsteady_separation",
      "vortex_shedding",
      "transitional_flow"
    ]
  },
  "verdict": "ADMITTED_WITH_CONTEXTUAL_BOUNDS"
}
EOF
    echo "rans-rung terminal receipt emitted to ${RECEIPT_FILE}"
    echo "rans-rung stages .5.8.1, .5.8.2, .5.8.3, and .5.8.4: PASS"
    ;;
  *)
    echo "e2e_extreal_rans_rung: unknown mode $MODE (available: --card, --solver, --validation, --adjudicate, --all)" >&2
    exit 30
    ;;
esac
