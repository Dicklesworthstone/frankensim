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
echo "[horizon-activation] emitting current population receipt"
cat <<'JSON'
{"suite":"addendum-horizon","trigger":4,"proposal":"4",
 "disposition":"NoData",
 "reason":"no paying coupled-transient workload exists in the program yet",
 "controller_state":"off-instrumented-only"}
JSON

echo "[horizon-activation] PASS (gate battery green; disposition NoData — activation requires bound paying workloads)"
