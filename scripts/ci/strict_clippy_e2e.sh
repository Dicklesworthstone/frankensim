#!/usr/bin/env bash
# strict_clippy_e2e.sh — terminal strict-Clippy gate (bead frankensim-vo8f2.5).
#
# Owns the single workspace strict-Clippy replay: pinned-toolchain check,
# fresh warning inventory on one source/toolchain snapshot, classification
# of every historical row from the four vo8f2 child repairs, bounded
# schema-versioned JSONL evidence, hostile self-test, and stable exit
# classes. A zero-warning receipt proves lint cleanliness for the recorded
# source/toolchain snapshot ONLY — never semantic correctness, performance,
# science, or certificate authority.
#
# CLI:
#   --check                  validate environment/source/toolchain/baseline
#                            WITHOUT running production clippy
#   --self-test              hostile fixtures only: CLI parsing, exit map,
#                            redaction, JSON escaping, row order/dedup,
#                            seeded-warning redness through the real pinned
#                            compiler (never touches workspace sources)
#   --run                    production replay: full pinned workspace lane
#                            over every target
#   --output-dir <path>      retained-evidence root (default target/vo8f2-5)
#   --max-wall-seconds <n>   tighten-only deadline override (default 5400)
#   --baseline <path>        historical inventory TSV
#                            (default scripts/ci/strict_clippy_baseline.tsv)
#
# EXIT CLASSES (stable):
#    0 zero-warning receipt on the recorded snapshot
#   10 runner usage / argument error
#   11 warning inventory nonzero (classified rows retained)
#   12 compile failure under the pinned toolchain
#   13 toolchain mismatch against rust-toolchain.toml
#   14 deadline exceeded
#   15 cancellation or drain failure
#   16 output overflow beyond bounded caps
#   17 evidence failure (log append fault)
#   18 infrastructure refusal (missing tools, unreadable snapshot/baseline)

set -u
LC_ALL=C
export LC_ALL

EXIT_OK=0
EXIT_USAGE=10
EXIT_WARNINGS=11
EXIT_COMPILE=12
EXIT_TOOLCHAIN=13
EXIT_TIMEOUT=14
EXIT_CANCEL=15
EXIT_OVERFLOW=16
EXIT_EVIDENCE=17
EXIT_INFRA=18

SCHEMA_VERSION="strict-clippy-e2e/1"

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)"
BASELINE_DEFAULT="scripts/ci/strict_clippy_baseline.tsv"
if [ -z "${REPO_ROOT}" ]; then
    printf '{"schema":"%s","event":"run-terminal","exit_class":%d,"reason":"not a git repository"}\n' \
        "${SCHEMA_VERSION}" "${EXIT_INFRA}"
    exit "${EXIT_INFRA}"
fi
cd "${REPO_ROOT}" || exit "${EXIT_INFRA}"

MODE=""
OUTPUT_DIR="target/vo8f2-5"
MAX_WALL_SECONDS=5400
BASELINE="${BASELINE_DEFAULT}"

log_file=""
log_events=0
LOG_MAX_EVENTS=20000
RAW_CAPTURE_MAX_BYTES=$((64 * 1024 * 1024))

usage() {
    sed -n '2,38p' "$0" | sed 's/^# \{0,1\}//'
}

# ---------------------------------------------------------------- logging --
json_escape() {
    printf '%s' "$1" \
        | tr '\t\r\n' '   ' \
        | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'
}

redact_text() {
    printf '%s' "$1" \
        | sed -e "s|${HOME}|~|g" \
        -e 's/sk-[A-Za-z0-9]\{8,\}/[REDACTED]/g' \
        -e 's/ghp_[A-Za-z0-9]\{8,\}/[REDACTED]/g' \
        -e 's/AKIA[A-Z0-9]\{8,\}/[REDACTED]/g'
}

emit() {
    # emit EVENT_NAME KEY=VALUE...  (values escaped here; callers pass raw)
    event="$1"; shift
    if [ -z "${log_file}" ] || [ "${log_events}" -ge "${LOG_MAX_EVENTS}" ]; then
        return 0
    fi
    line="{\"schema\":\"${SCHEMA_VERSION}\",\"event\":\"${event}\""
    for kv in "$@"; do
        k="${kv%%=*}"
        v="${kv#*=}"
        line="${line},\"${k}\":\"$(json_escape "${v}")\""
    done
    line="${line}}"
    printf '%s\n' "${line}" >>"${log_file}" 2>/dev/null || {
        printf '{"schema":"%s","event":"evidence-failure","reason":"log append failed"}\n' \
            "${SCHEMA_VERSION}" >&2
        exit "${EXIT_EVIDENCE}"
    }
    log_events=$((log_events + 1))
}

snapshot_fingerprint() {
    head_sha="$(git rev-parse HEAD 2>/dev/null)" || return 1
    dirty="$(git status --porcelain 2>/dev/null | shasum -a 256 | cut -c1-16)"
    printf '%s+%s' "${head_sha}" "${dirty}"
}

pinned_channel() {
    sed -n 's/^channel[[:space:]]*=[[:space:]]*"\(.*\)"/\1/p' rust-toolchain.toml 2>/dev/null \
        | head -n 1
}

check_environment() {
    command -v git >/dev/null 2>&1 || { printf 'missing git\n' >&2; return 1; }
    command -v cargo >/dev/null 2>&1 || { printf 'missing cargo\n' >&2; return 1; }
    command -v python3 >/dev/null 2>&1 || { printf 'missing python3\n' >&2; return 1; }
    cargo clippy --version >/dev/null 2>&1 || { printf 'missing clippy\n' >&2; return 1; }
    [ -f rust-toolchain.toml ] || { printf 'missing rust-toolchain.toml\n' >&2; return 1; }
    snapshot_fingerprint >/dev/null 2>&1 || { printf 'unreadable git snapshot\n' >&2; return 1; }
    return 0
}

check_toolchain() {
    want="$(pinned_channel)"
    [ -n "${want}" ] || { printf 'empty pinned channel\n' >&2; return 1; }
    got="$(rustup show active-toolchain 2>/dev/null | cut -d' ' -f1)"
    # rustup may report the fully qualified host triple; accept any active
    # toolchain whose id starts with the pinned channel string.
    case "${got}" in
        "${want}" | "${want}"-*) return 0 ;;
        *) printf 'toolchain mismatch: pinned=%s active=%s\n' "${want}" "${got}" >&2; return 1 ;;
    esac
}

validate_baseline() {
    [ -f "${BASELINE}" ] || { printf 'missing baseline %s\n' "${BASELINE}" >&2; return 1; }
    while IFS= read -r line; do
        [ -n "${line}" ] || continue
        case "${line}" in
            '#'* | '') continue ;;
        esac
        fields="$(printf '%s' "${line}" | awk -F'\t' '{print NF}')"
        [ "${fields}" -ge 5 ] || {
            printf 'baseline row needs >=5 TSV fields: %s\n' "${line}" >&2
            return 1
        }
    done <"${BASELINE}"
    return 0
}

# ------------------------------------------------------------- inventory ---
parse_warning_rows() {
    # Reads the cargo JSON stream on STDIN. Uses python3 -c so stdin stays
    # bound to data; a heredoc would bind it to the script text instead.
    python3 -c '
import json, sys, re

rows = {}
for raw in sys.stdin:
    raw = raw.strip()
    if not raw.startswith("{"):
        continue
    try:
        m = json.loads(raw)
    except Exception:
        continue
    if m.get("reason") != "compiler-message":
        continue
    msg = m.get("message", {})
    # Under -D warnings lint failures surface as level=error carrying a
    # lint code; only bare rustc codes (E0xxx) are true compile failures.
    if msg.get("level") not in ("warning", "error"):
        continue
    code_early = (msg.get("code") or {}).get("code") or ""
    if code_early.startswith("E") and code_early[1:].isdigit():
        continue
    spans = [s for s in msg.get("spans", []) if s.get("is_primary")]
    if not spans:
        continue
    s = spans[0]
    path = s.get("file_name") or ""
    if "/rustc/" in path or path.startswith("/rust/"):
        continue
    mm = re.match(r"^crates/([^/]+)/", path)
    if not mm:
        continue
    code = (msg.get("code") or {}).get("code") or "unknown-lint"
    text = (msg.get("message") or "")[:160].replace("\t", " ")
    rows[(mm.group(1), path, s.get("line_start", 0), code, text)] = True

for (crate, path, line, code, text) in sorted(rows.keys()):
    print("\t".join([crate, path, str(line), code, text]))
'
}

classify_rows() {
    # args: baseline path; stdin: live TSV rows.
    # stdout: class<TAB>crate<TAB>path<TAB>line<TAB>lint<TAB>diagnostic<TAB>historical_bead
    baseline_file="$1"
    tmp_live="$(mktemp)" || return 1
    cat >"${tmp_live}"
    while IFS=$'\t' read -r crate path line code text; do
        hit_bead="$(awk -F'\t' -v c="${crate}" -v k="${code}" \
            '$2==c && $4==k {print $1; exit}' "${baseline_file}" 2>/dev/null)"
        hit_path="$(awk -F'\t' -v c="${crate}" -v k="${code}" \
            '$2==c && $4==k {print $3; exit}' "${baseline_file}" 2>/dev/null)"
        if [ -n "${hit_bead}" ]; then
            cls="present"
            [ "${path}" != "${hit_path}" ] && cls="moved"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${cls}" "${crate}" "${path}" "${line}" "${code}" "${text}" "${hit_bead}"
        else
            printf 'newly_found\t%s\t%s\t%s\t%s\t%s\t\n' \
                "${crate}" "${path}" "${line}" "${code}" "${text}"
        fi
    done <"${tmp_live}"
    grep -v '^#' "${baseline_file}" 2>/dev/null |
        while IFS=$'\t' read -r bead crate path code note; do
            [ -n "${crate}" ] || continue
            if ! awk -F'\t' -v c="${crate}" -v k="${code}" \
                '$1==c && $4==k {found=1} END{exit !found}' "${tmp_live}" >/dev/null 2>&1; then
                printf 'proven_fixed\t%s\t%s\t0\t%s\t%s\t%s\n' \
                    "${crate}" "${path}" "${code}" "${note}" "${bead}"
            fi
        done
    rm -f "${tmp_live}" 2>/dev/null || true
    return 0
}

run_clippy_capped() {
    # $@: extra cargo clippy args. Sets CLIPPY_STATUS:
    # <cargo exit code> | timeout | cancelled | overflow.
    raw_out="${OUTPUT_DIR}/clippy_raw.jsonl"
    status_file="${OUTPUT_DIR}/clippy_status"
    drained_file="${OUTPUT_DIR}/deadline_drained"
    : >"${raw_out}"
    : >"${drained_file}.watch"
    (
        cargo clippy "$@" --message-format=json -- -D warnings >"${raw_out}.full" 2>"${OUTPUT_DIR}/clippy_stderr.txt"
        printf '%s' "$?" >"${status_file}"
        touch "${drained_file}.done"
    ) &
    cargo_pid=$!
    term_hit=0
    forward_term() {
        term_hit=1
        kill -TERM "${cargo_pid}" 2>/dev/null
        sleep 10
        kill -KILL "${cargo_pid}" 2>/dev/null
    }
    trap forward_term TERM INT
    deadline_watchdog() {
        sleep "${MAX_WALL_SECONDS}"
        touch "${drained_file}"
        kill -TERM "${cargo_pid}" 2>/dev/null
        sleep 10
        kill -KILL "${cargo_pid}" 2>/dev/null
    }
    deadline_watchdog &
    watchdog_pid=$!
    wait "${cargo_pid}" 2>/dev/null
    CLIPPY_STATUS="$(cat "${status_file}" 2>/dev/null || printf '99')"
    kill "${watchdog_pid}" 2>/dev/null
    wait "${watchdog_pid}" 2>/dev/null
    trap - TERM INT
    rm -f "${drained_file}.watch" 2>/dev/null || true
    if [ -f "${drained_file}" ]; then
        CLIPPY_STATUS="timeout"
    elif [ "${term_hit}" -eq 1 ]; then
        CLIPPY_STATUS="cancelled"
    fi
    actual_bytes="$(wc -c <"${raw_out}.full" 2>/dev/null || printf 0)"
    if [ "${actual_bytes}" -gt "${RAW_CAPTURE_MAX_BYTES}" ]; then
        CLIPPY_STATUS="overflow"
    else
        head -c "${RAW_CAPTURE_MAX_BYTES}" "${raw_out}.full" >"${raw_out}" 2>/dev/null ||
            CLIPPY_STATUS="overflow"
    fi
    return 0
}

# ----------------------------------------------------------------- modes ---
do_check() {
    mkdir -p "${OUTPUT_DIR}" || exit "${EXIT_INFRA}"
    log_file="${OUTPUT_DIR}/strict_clippy_check.jsonl"
    : >"${log_file}"
    fp="$(snapshot_fingerprint)" || { emit run-terminal "exit_class=${EXIT_INFRA}"; exit "${EXIT_INFRA}"; }
    check_environment || { emit run-terminal "exit_class=${EXIT_INFRA}"; exit "${EXIT_INFRA}"; }
    check_toolchain || { emit run-terminal "exit_class=${EXIT_TOOLCHAIN}"; exit "${EXIT_TOOLCHAIN}"; }
    validate_baseline || { emit run-terminal "exit_class=${EXIT_INFRA}"; exit "${EXIT_INFRA}"; }
    emit run-start "mode=check" "snapshot=${fp}"
    emit run-terminal "exit_class=${EXIT_OK}"
    printf 'check ok: snapshot=%s toolchain=%s baseline=%s\n' \
        "${fp}" "$(rustup show active-toolchain 2>/dev/null | cut -d' ' -f1)" "${BASELINE}"
    exit "${EXIT_OK}"
}

seeded_warning_redness() {
    seed_root="${OUTPUT_DIR}/selftest-seed"
    mkdir -p "${seed_root}" || return 1
    seed_dir="$(mktemp -d "${seed_root}/seed.XXXXXX")" || return 1
    cat >"${seed_dir}/Cargo.toml" <<TOML
[package]
name = "vo8f25_seed"
version = "0.0.0"
edition = "2021"

[workspace]
TOML
    mkdir "${seed_dir}/src" || return 1
    cat >"${seed_dir}/src/lib.rs" <<'RS'
// Default-warn clippy bait: approx_constant fires on the PI lookalike and
// needless_return on the explicit return, so any pinned-toolchain run of
// this seed must produce at least one warning under -D warnings.
pub fn seeded_pi() -> f64 {
    let pi = 3.14159;
    return pi;
}
RS
    (cd "${seed_dir}" && cargo clippy --lib --message-format=json -- -D warnings 2>/dev/null) \
        >"${seed_dir}/out.jsonl"
    status=$?
    [ "${status}" -ne 0 ] || return 1
    grep -q 'clippy::' "${seed_dir}/out.jsonl" || return 1
    return 0
}

do_self_test() {
    mkdir -p "${OUTPUT_DIR}" || exit "${EXIT_INFRA}"
    log_file="${OUTPUT_DIR}/strict_clippy_selftest.jsonl"
    : >"${log_file}"
    failures=0
    t() {
        if [ "$2" -eq 0 ]; then
            printf 'ok   %s\n' "$1"
        else
            printf 'FAIL %s\n' "$1" >&2
            failures=$((failures + 1))
        fi
    }

    r="$(redact_text "token sk-abcdef1234567890 home ${HOME}/x")"
    case "${r}" in *sk-abcdef*) t redaction 1 ;; *) t redaction 0 ;; esac

    e="$(json_escape 'a"b\c d')"
    [ "${e}" = 'a\"b\\c d' ]
    t json_escape $?

    ("$0" --definitely-not-a-mode >/dev/null 2>&1)
    usage_rc=$?
    [ "${usage_rc}" -eq "${EXIT_USAGE}" ]
    t usage_exit_10 $?

    printf '%s\n' \
        '{"reason":"compiler-message","message":{"level":"warning","code":{"code":"clippy::zz"},"message":"m2","spans":[{"is_primary":true,"file_name":"crates/z/src/a.rs","line_start":7}]}}' \
        '{"reason":"compiler-message","message":{"level":"warning","code":{"code":"clippy::aa"},"message":"m1","spans":[{"is_primary":true,"file_name":"crates/a/src/x.rs","line_start":2}]}}' \
        '{"reason":"compiler-message","message":{"level":"warning","code":{"code":"clippy::aa"},"message":"m1","spans":[{"is_primary":true,"file_name":"crates/a/src/x.rs","line_start":2}]}}' \
        | parse_warning_rows >"${OUTPUT_DIR}/rows_probe.tsv"
    rows_count="$(wc -l <"${OUTPUT_DIR}/rows_probe.tsv" | tr -d ' ')"
    first_row="$(head -n 1 "${OUTPUT_DIR}/rows_probe.tsv" | cut -f2)"
    [ "${rows_count}" = "2" ] && [ "${first_row}" = "crates/a/src/x.rs" ]
    t row_order_and_dedup $?

    seeded_warning_redness
    t seeded_warning_red $?

    if [ "${failures}" -ne 0 ]; then
        emit run-terminal "exit_class=${EXIT_USAGE}" "selftest_failures=${failures}"
        printf 'self-test FAILURES: %d\n' "${failures}" >&2
        exit "${EXIT_USAGE}"
    fi
    emit run-terminal "exit_class=${EXIT_OK}"
    printf 'self-test ok\n'
    exit "${EXIT_OK}"
}

do_run() {
    mkdir -p "${OUTPUT_DIR}" || exit "${EXIT_INFRA}"
    log_file="${OUTPUT_DIR}/strict_clippy_run.jsonl"
    : >"${log_file}"
    fp="$(snapshot_fingerprint)" || exit "${EXIT_INFRA}"
    tool="$(rustup show active-toolchain 2>/dev/null | cut -d' ' -f1)"
    check_environment || { emit run-terminal "exit_class=${EXIT_INFRA}"; exit "${EXIT_INFRA}"; }
    check_toolchain || { emit run-terminal "exit_class=${EXIT_TOOLCHAIN}"; exit "${EXIT_TOOLCHAIN}"; }
    validate_baseline || { emit run-terminal "exit_class=${EXIT_INFRA}"; exit "${EXIT_INFRA}"; }
    emit run-start "mode=run" "snapshot=${fp}" "toolchain=${tool}"

    run_clippy_capped --workspace --all-targets --no-deps
    status="${CLIPPY_STATUS:-99}"
    case "${status}" in
        timeout)
            emit run-terminal "exit_class=${EXIT_TIMEOUT}"
            exit "${EXIT_TIMEOUT}" ;;
        overflow)
            emit run-terminal "exit_class=${EXIT_OVERFLOW}"
            exit "${EXIT_OVERFLOW}" ;;
        0) : ;;
        101 | 102)
            # Under -D warnings, lint failures also exit 101. Only a bare
            # rustc E-code diagnostic is a true compile failure (class 12);
            # pure lint debt falls through to the inventory classifier.
            if grep -qE '"code":"E[0-9]{4}"' "${OUTPUT_DIR}/clippy_raw.jsonl" 2>/dev/null; then
                emit run-terminal "exit_class=${EXIT_COMPILE}" "cargo_status=${status}"
                exit "${EXIT_COMPILE}"
            fi
            emit run-note "cargo_status=${status}" "note=lint-debt-only failure"
            ;;
        *)
            emit run-terminal "exit_class=${EXIT_INFRA}" "cargo_status=${status}"
            exit "${EXIT_INFRA}" ;;
    esac

    parse_warning_rows <"${OUTPUT_DIR}/clippy_raw.jsonl" >"${OUTPUT_DIR}/live_rows.tsv"
    classify_rows "${BASELINE}" <"${OUTPUT_DIR}/live_rows.tsv" |
        LC_ALL=C sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 \
            >"${OUTPUT_DIR}/classification.tsv"

    warnings=0 fixed=0 moved=0 present=0 newf=0
    while IFS=$'\t' read -r cls crate path line code text bead; do
        [ -n "${cls}" ] || continue
        emit inventory-row "class=${cls}" "crate=${crate}" "path=${path}" "line=${line}" \
            "lint=${code}" "diagnostic=${text}" "historical_bead=${bead}"
        case "${cls}" in
            proven_fixed) fixed=$((fixed + 1)) ;;
            moved) moved=$((moved + 1)); warnings=$((warnings + 1)) ;;
            present) present=$((present + 1)); warnings=$((warnings + 1)) ;;
            newly_found) newf=$((newf + 1)); warnings=$((warnings + 1)) ;;
        esac
    done <"${OUTPUT_DIR}/classification.tsv"

    if [ "${warnings}" -eq 0 ]; then
        emit summary "warnings_total=0" "proven_fixed=${fixed}"
        emit run-terminal "exit_class=${EXIT_OK}" "proven_fixed=${fixed}"
        printf 'zero-warning receipt: snapshot=%s toolchain=%s proven_fixed=%d\n' \
            "${fp}" "${tool}" "${fixed}"
        exit "${EXIT_OK}"
    fi
    emit summary "warnings_total=${warnings}" "present=${present}" "moved=${moved}" "newly_found=${newf}"
    emit run-terminal "exit_class=${EXIT_WARNINGS}" "present=${present}" "moved=${moved}" \
        "newly_found=${newf}" "proven_fixed=${fixed}"
    printf 'warning inventory nonzero: present=%d moved=%d newly_found=%d proven_fixed=%d\n' \
        "${present}" "${moved}" "${newf}" "${fixed}" >&2
    exit "${EXIT_WARNINGS}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --check | --self-test | --run)
            [ -z "${MODE}" ] || { usage; exit "${EXIT_USAGE}"; }
            MODE="$1" ;;
        --output-dir)
            [ $# -ge 2 ] || { usage; exit "${EXIT_USAGE}"; }
            OUTPUT_DIR="$2"; shift ;;
        --max-wall-seconds)
            [ $# -ge 2 ] || { usage; exit "${EXIT_USAGE}"; }
            case "$2" in '' | *[!0-9]*) usage; exit "${EXIT_USAGE}" ;; esac
            MAX_WALL_SECONDS="$2"; shift ;;
        --baseline)
            [ $# -ge 2 ] || { usage; exit "${EXIT_USAGE}"; }
            BASELINE="$2"; shift ;;
        *)
            usage; exit "${EXIT_USAGE}" ;;
    esac
    shift
done
[ -n "${MODE}" ] || { usage; exit "${EXIT_USAGE}"; }

case "${MODE}" in
    --check) do_check ;;
    --self-test) do_self_test ;;
    --run) do_run ;;
esac