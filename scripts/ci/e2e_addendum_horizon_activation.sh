#!/usr/bin/env bash
# Horizon trigger adjudication lane (bead frankensim-epic-addendum-xpck.5.4,
# feeding xpck.5.8's quarterly adjudication).
#
# Executes the splitting-error demand-gate battery and emits today's honest
# population receipt. The gate engine is live and tested; the controller is
# NOT activated here — activation requires bound paying coupled-transient
# workloads with complete error budgets, which do not exist yet (NoData).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly SCRIPT_DIR
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
readonly REPO_ROOT
readonly CARGO_BIN="${CARGO_BIN:-cargo}"

echo "[horizon-activation] running fs-govern horizon_splitting battery"
"$CARGO_BIN" test -p fs-govern --test horizon_splitting

echo "[horizon-activation] running fs-govern horizon_metrology battery (Trigger 11)"
"$CARGO_BIN" test -p fs-govern --test horizon_metrology

echo "[horizon-activation] running fs-govern horizon_symmetry battery (Trigger 13b)"
"$CARGO_BIN" test -p fs-govern --test horizon_symmetry

echo "[horizon-activation] running fs-govern horizon_explanation battery (Trigger B)"
"$CARGO_BIN" test -p fs-govern --test horizon_explanation

echo "[horizon-activation] running fs-govern horizon_goodhart battery (Trigger D)"
"$CARGO_BIN" test -p fs-govern --test horizon_goodhart

echo "[horizon-activation] emitting current population receipts"
cat <<'JSON'
{"suite":"addendum-horizon","trigger":4,"proposal":"4",
 "disposition":"NoData",
 "reason":"no paying coupled-transient workload exists in the program yet",
 "controller_state":"off-instrumented-only"}
{"suite":"addendum-horizon","trigger":11,"proposal":"11",
 "disposition":"NoData",
 "reason":"no authorized metrology partnership or retained scan data exists in the program yet",
 "fallback_state":"point-sensor-assimilation-active"}
{"suite":"addendum-horizon","trigger":13,"proposal":"13b",
 "disposition":"NoData",
 "reason":"no representative real-workload symmetry census exists in the program yet",
 "solver_state":"detection-only"}
{"suite":"addendum-horizon","trigger":"B","proposal":"B",
 "disposition":"Rule4Defer",
 "reason":"Rule 4: human-driven operator mode defers explanation-object activation",
 "explanation_state":"deferred"}
{"suite":"addendum-horizon","trigger":"D","proposal":"D",
 "disposition":"Rule4Defer",
 "reason":"Rule 4: human-driven operator mode defers Goodhart guard activation",
 "goodhart_state":"deferred"}
JSON

echo "[horizon-activation] PASS (gate batteries green; disposition NoData/Rule4Defer — activation requires bound empirical evidence)"

