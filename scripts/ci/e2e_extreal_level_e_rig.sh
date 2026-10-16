#!/usr/bin/env bash
# Level-E rig and metrology specification CI verification lane
# (bead frankensim-extreal-program-f85xj.4.5.2).
#
# Validates Level-E rig specification, instrument calibration chains,
# energy-balance interval gate, and blind-holdout sealing.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly SCRIPT_DIR
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
readonly REPO_ROOT
readonly CARGO_BIN="${CARGO_BIN:-cargo}"

echo "[level-e-rig] running fs-vvreg rig test battery"
"$CARGO_BIN" test -p fs-vvreg --test rig

echo "[level-e-rig] running fs-vvreg campaign test battery"
"$CARGO_BIN" test -p fs-vvreg --test campaign

echo "[level-e-rig] PASS (rig specification, campaign operating matrix, and ingest tests green; energy-balance gate verified)"
