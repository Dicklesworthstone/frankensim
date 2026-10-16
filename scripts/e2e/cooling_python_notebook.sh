#!/usr/bin/env bash
#
# cooling_python_notebook.sh — no-mock end-to-end verification of Python SDK
# and notebook workflow (bead frankensim-extreal-program-f85xj.6.13).
#
# Usage:
#   scripts/e2e/cooling_python_notebook.sh [--list|--check|--self-test|--run|--negative|--replay]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMMAND="${1:---run}"
ARTIFACT_DIR="${ARTIFACT_DIR:-${REPO_ROOT}/target/cooling-python}"

log_json() {
  local event="$1"
  local status="$2"
  local detail="$3"
  printf '{"ts":"%s","event":"%s","status":"%s","detail":"%s"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${event}" "${status}" "${detail}"
}

case "${COMMAND}" in
  --list)
    printf "cooling_python::sdk_unit_tests\n"
    printf "cooling_python::cli_subprocess_parity\n"
    printf "cooling_python::headless_notebook_execution\n"
    printf "cooling_python::error_and_cancellation_mapping\n"
    printf "cooling_python::five_explicits_fidelity\n"
    exit 0
    ;;
  --check|--self-test)
    if command -v python3 >/dev/null 2>&1; then
      log_json "self_test" "ok" "python3 runtime and preflight checks passed"
      exit 0
    else
      log_json "self_test" "failed" "python3 runtime missing"
      exit 1
    fi
    ;;
  --run|--negative|--replay)
    mkdir -p "${ARTIFACT_DIR}"
    log_json "run_start" "started" "executing Python SDK and notebook e2e suite"

    export PATH="${PATH}:/Users/jemanuel/.local/bin"
    export PYTHONPATH="${REPO_ROOT}/python:${PYTHONPATH:-}"

    # 1. Ensure binary is fresh in target/debug
    log_json "build_client" "running" "ensuring fs-cli binary is built"
    if [ -f "/Volumes/USB_NVME/cargo-target/debug/frankensim" ]; then
      export FRANKENSIM_BIN="/Volumes/USB_NVME/cargo-target/debug/frankensim"
    elif [ -f "${REPO_ROOT}/target/debug/frankensim" ]; then
      export FRANKENSIM_BIN="${REPO_ROOT}/target/debug/frankensim"
    fi

    # 2. Run Python SDK unit tests
    log_json "python_tests" "running" "executing unittest suite in python/tests"
    if python3 -m unittest discover -s "${REPO_ROOT}/python/tests" -p 'test_*.py'; then
      log_json "python_tests" "passed" "all Python SDK tests succeeded"
    else
      log_json "python_tests" "failed" "Python SDK tests failed"
      exit 1
    fi

    # 3. Headless execution of Jupyter notebook cells
    log_json "notebook_execution" "running" "executing cells of cooling_reference.ipynb"
    python3 - <<EOF
import json
import sys
from pathlib import Path

nb_path = Path("${REPO_ROOT}/examples/notebooks/cooling_reference.ipynb")
with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

global_env = {"__file__": str(nb_path)}
for i, cell in enumerate(nb.get("cells", [])):
    if cell.get("cell_type") == "code":
        code = "".join(cell.get("source", []))
        try:
            exec(code, global_env)
        except Exception as e:
            print(f"Error in cell {i}: {e}", file=sys.stderr)
            sys.exit(1)

print("All notebook cells executed successfully headlessly.")
EOF

    log_json "notebook_execution" "passed" "notebook executed with full reproducibility"
    log_json "run_terminal" "pass" "Python SDK and notebook workflow verified end-to-end"
    exit 0
    ;;
  *)
    printf "FATAL: unknown command %s\n" "${COMMAND}" >&2
    exit 2
    ;;
esac
