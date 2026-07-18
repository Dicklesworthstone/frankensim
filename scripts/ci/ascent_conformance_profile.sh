#!/usr/bin/env bash
# ASCENT conformance profile runner.
#
# The PR profile is a deliberately small, family-representative selector. The
# nightly profile is the complete fs-ascent package with Cargo fail-fast
# disabled. Both profiles are lockfile-closed and enforce an aggregate monotonic
# budget that intentionally includes compilation. This is an execution guard,
# not a machine-independent performance claim.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly SCRIPT_DIR
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
readonly REPO_ROOT
readonly CARGO_BIN="${CARGO_BIN:-cargo}"
readonly PYTHON_BIN="${PYTHON_BIN:-python3}"
readonly PR_BUDGET_SECONDS="${FS_ASCENT_PR_BUDGET_SECONDS:-900}"
readonly NIGHTLY_BUDGET_SECONDS="${FS_ASCENT_NIGHTLY_BUDGET_SECONDS:-7200}"

usage() {
  cat >&2 <<'EOF'
usage:
  scripts/ci/ascent_conformance_profile.sh pr
  scripts/ci/ascent_conformance_profile.sh nightly
  scripts/ci/ascent_conformance_profile.sh --list pr
  scripts/ci/ascent_conformance_profile.sh --list nightly
  scripts/ci/ascent_conformance_profile.sh --self-test

environment:
  CARGO_BIN                         Cargo executable path (default: cargo)
  PYTHON_BIN                        Python 3 executable path (default: python3)
  FS_ASCENT_PR_BUDGET_SECONDS       Aggregate PR wall budget (default: 900)
  FS_ASCENT_NIGHTLY_BUDGET_SECONDS  Aggregate nightly wall budget (default: 7200)
  FS_ASCENT_PROFILE_LOG_DIR         Retained-run root (default: target/ascent-conformance-profile)
  FS_ASCENT_TERMINATION_GRACE_SECONDS  Grace before KILL (default: 5)
  FS_ASCENT_KILL_DRAIN_SECONDS      Bounded KILL drain wait (default: 5)

The aggregate monotonic budget includes compilation and test execution. Callers may
override it for a declared host/target policy; the emitted receipt retains the
effective value. The self-test uses internal fake processes and never invokes
real Cargo.
EOF
}

run_python_supervisor() {
  local mode="$1"
  local profile="$2"
  exec "$PYTHON_BIN" -I - "$mode" "$profile" "$REPO_ROOT" \
    "$SCRIPT_DIR/ascent_conformance_profile.sh" "$CARGO_BIN" \
    "$PR_BUDGET_SECONDS" "$NIGHTLY_BUDGET_SECONDS" <<'PY'
import atexit
import hashlib
import json
import os
import pathlib
import signal
import stat
import subprocess
import sys
import tempfile
import time
import traceback

(
    MODE,
    PROFILE,
    REPO_ROOT_TEXT,
    SCRIPT_PATH_TEXT,
    CARGO_BIN,
    PR_BUDGET_TEXT,
    NIGHTLY_BUDGET_TEXT,
) = sys.argv[1:]

REPO_ROOT = pathlib.Path(REPO_ROOT_TEXT).resolve()
SCRIPT_PATH = pathlib.Path(SCRIPT_PATH_TEXT).resolve()
PR_TARGETS = (
    "ascent_battery",
    "constrained_battery",
    "pareto_battery",
    "runner_battery",
    "budget_trend_manifest",
    "bbob_budget_ledger",
)
MAX_BUDGET_SECONDS = 7 * 24 * 60 * 60
MAX_DRAIN_SECONDS = 60.0
REQUESTED_SIGNAL = None


class RunInterrupted(Exception):
    def __init__(self, signum):
        super().__init__(signum)
        self.signum = signum


class IdentityError(RuntimeError):
    pass


def canonical_bytes(value):
    return (
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        + b"\n"
    )


def emit_stdout(value):
    payload = canonical_bytes(value)
    try:
        sys.stdout.buffer.write(payload)
        sys.stdout.buffer.flush()
        return True
    except (BrokenPipeError, OSError):
        return False


def parse_positive_integer(label, text):
    try:
        value = int(text, 10)
    except ValueError as error:
        raise ValueError(f"{label} must be a positive integer, got {text!r}") from error
    if value <= 0 or value > MAX_BUDGET_SECONDS:
        raise ValueError(
            f"{label} must be in 1..={MAX_BUDGET_SECONDS}, got {value}"
        )
    return value


def parse_positive_seconds(label, text):
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError(f"{label} must be a positive number, got {text!r}") from error
    if not (0.0 < value <= MAX_DRAIN_SECONDS):
        raise ValueError(
            f"{label} must be in (0,{MAX_DRAIN_SECONDS}], got {value}"
        )
    return value


def budget_for_profile(profile):
    if profile == "pr":
        return parse_positive_integer("PR budget", PR_BUDGET_TEXT)
    if profile == "nightly":
        return parse_positive_integer("nightly budget", NIGHTLY_BUDGET_TEXT)
    raise ValueError(f"unknown profile {profile!r}")


def command_for_target(profile, target):
    command = [CARGO_BIN, "test", "--locked", "-p", "fs-ascent"]
    if profile == "pr":
        command.extend(["--test", target])
    else:
        command.append("--no-fail-fast")
    command.extend(["--", "--nocapture"])
    return command


def list_profile(profile):
    budget = budget_for_profile(profile)
    emit_stdout(
        {
            "schema": "frankensim-ascent-conformance-profile-v2",
            "profile": profile,
            "budget_seconds": budget,
            "build_time_included": True,
            "deadline_enforced": True,
            "deadline_clock": "monotonic",
            "receipt_artifact": "verdicts.jsonl",
            "child_output": "per-target retained logs",
            "provenance_status": "unsealed",
        }
    )
    targets = PR_TARGETS if profile == "pr" else ("all",)
    for target in targets:
        emit_stdout(
            {
                "schema": "frankensim-ascent-conformance-selector-v2",
                "profile": profile,
                "package": "fs-ascent",
                "target": target,
                "command": command_for_target(profile, target),
                "provenance_status": "unsealed",
            }
        )
    return 0


def sanitized_git_environment():
    environment = os.environ.copy()
    for name in tuple(environment):
        if name.startswith("GIT_"):
            environment.pop(name, None)
    environment["LC_ALL"] = "C"
    return environment


def git_bytes(*arguments):
    command = ["git", "-c", "core.excludesFile=/dev/null", *arguments]
    try:
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=sanitized_git_environment(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=20,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise IdentityError(
            f"Git observation failed for {arguments!r}: "
            f"{type(error).__name__}: {error}"
        ) from error
    if result.returncode != 0:
        detail = result.stderr.decode("utf-8", "backslashreplace").strip()
        raise IdentityError(
            f"Git observation {arguments!r} exited {result.returncode}: {detail}"
        )
    return result.stdout


def hash_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def feed_record(digest, label, payload):
    label_bytes = label.encode("ascii")
    digest.update(len(label_bytes).to_bytes(4, "big"))
    digest.update(label_bytes)
    digest.update(len(payload).to_bytes(8, "big"))
    digest.update(payload)


def root_materialized_tree_sha256():
    inventory = git_bytes(
        "-c",
        "core.fileMode=true",
        "ls-files",
        "-z",
        "--cached",
        "--others",
        "--exclude-per-directory=.gitignore",
    )
    relative_paths = sorted(set(path for path in inventory.split(b"\0") if path))
    root_bytes = os.fsencode(str(REPO_ROOT))
    digest = hashlib.sha256()
    feed_record(digest, "domain", b"frankensim-root-materialized-tree-v1")
    for relative_path in relative_paths:
        feed_record(digest, "path", relative_path)
        full_path = os.path.join(root_bytes, relative_path)
        try:
            before = os.lstat(full_path)
        except FileNotFoundError:
            feed_record(digest, "kind", b"missing")
            continue
        mode = stat.S_IFMT(before.st_mode) | stat.S_IMODE(before.st_mode)
        feed_record(digest, "mode", str(mode).encode("ascii"))
        if stat.S_ISREG(before.st_mode):
            content = hashlib.sha256()
            with open(full_path, "rb") as handle:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    content.update(chunk)
            after = os.lstat(full_path)
            stable_fields = (
                before.st_dev,
                before.st_ino,
                before.st_mode,
                before.st_size,
                before.st_mtime_ns,
            )
            observed_fields = (
                after.st_dev,
                after.st_ino,
                after.st_mode,
                after.st_size,
                after.st_mtime_ns,
            )
            if stable_fields != observed_fields:
                raise IdentityError(
                    f"source file moved while hashing: "
                    f"{os.fsdecode(relative_path)!r}"
                )
            feed_record(digest, "size", str(before.st_size).encode("ascii"))
            feed_record(digest, "content-sha256", content.hexdigest().encode("ascii"))
        elif stat.S_ISLNK(before.st_mode):
            target = os.readlink(full_path)
            target_bytes = target if isinstance(target, bytes) else os.fsencode(target)
            feed_record(digest, "symlink-target", target_bytes)
        elif stat.S_ISDIR(before.st_mode):
            feed_record(digest, "kind", b"directory")
        else:
            feed_record(digest, "kind", b"special")
    return digest.hexdigest()


def repo_identity():
    head = git_bytes("rev-parse", "--verify", "HEAD").strip().decode("ascii")
    head_tree = (
        git_bytes("rev-parse", "--verify", "HEAD^{tree}").strip().decode("ascii")
    )
    index_rows = git_bytes(
        "-c", "core.fileMode=true", "ls-files", "-s", "-z"
    )
    status_rows = git_bytes(
        "-c",
        "core.fileMode=true",
        "status",
        "--porcelain=v1",
        "-z",
        "--untracked-files=all",
    )
    cargo_lock = REPO_ROOT / "Cargo.lock"
    if not cargo_lock.is_file():
        raise IdentityError("required Cargo.lock is missing")
    constellation_lock = REPO_ROOT / "constellation.lock"
    return {
        "head_sha": head,
        "head_tree_sha": head_tree,
        "index_sha256": hashlib.sha256(index_rows).hexdigest(),
        "root_tree_sha256": root_materialized_tree_sha256(),
        "git_status_sha256": hashlib.sha256(status_rows).hexdigest(),
        "dirty": bool(status_rows),
        "cargo_lock_sha256": hash_file(cargo_lock),
        "constellation_lock_sha256": (
            hash_file(constellation_lock) if constellation_lock.is_file() else None
        ),
    }


class VerdictWriter:
    def __init__(self, run_dir, profile, runner_argv):
        self.run_dir = pathlib.Path(run_dir)
        self.path = self.run_dir / "verdicts.jsonl"
        self.profile = profile
        self.runner_argv = tuple(runner_argv)
        self.handle = open(self.path, "xb", buffering=0)
        self.rows = 0
        self.sealed = False
        self.stdout_mirror_ok = True
        self.source_before = None

    def _write(self, row):
        payload = canonical_bytes(row)
        self.handle.write(payload)
        os.fsync(self.handle.fileno())
        self.rows += 1
        self.stdout_mirror_ok = emit_stdout(row) and self.stdout_mirror_ok

    def append(self, row):
        if self.sealed:
            raise RuntimeError("cannot append after proof seal")
        self._write(row)

    def seal_once(
        self,
        *,
        status,
        terminal_exit_code,
        provenance_state,
        source_after,
        detail,
    ):
        if self.sealed:
            return
        os.fsync(self.handle.fileno())
        prefix_hash = hash_file(self.path)
        seal = {
            "identity_domain": "org.frankensim.ci.ascent-conformance-proof.v1",
            "identity_version": 1,
            "schema": "frankensim-ascent-conformance-proof-seal-v1",
            "event": "proof-seal",
            "profile": self.profile,
            "status": status,
            "terminal_exit_code": terminal_exit_code,
            "provenance_state": provenance_state,
            "source_before": self.source_before,
            "source_after": source_after,
            "runner_argv": list(self.runner_argv),
            "verdicts_prefix_sha256": prefix_hash,
            "prefix_rows": self.rows,
            "run_dir": str(self.run_dir),
            "verdicts_path": str(self.path),
            "stdout_mirror_ok": self.stdout_mirror_ok,
            "detail": detail,
        }
        self.sealed = True
        self._write(seal)

    def close(self):
        try:
            self.handle.close()
        except OSError:
            pass


def remember_signal(signum, _frame):
    global REQUESTED_SIGNAL
    if REQUESTED_SIGNAL is None:
        REQUESTED_SIGNAL = signum


def install_signal_handlers():
    for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, remember_signal)


def ignore_control_signals():
    for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, signal.SIG_IGN)


def close_control_signal_admission():
    global REQUESTED_SIGNAL
    controlled = (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)
    try:
        signal.pthread_sigmask(signal.SIG_BLOCK, controlled)
        if REQUESTED_SIGNAL is None:
            pending = signal.sigpending()
            for signum in controlled:
                if signum in pending:
                    REQUESTED_SIGNAL = signum
                    break
    except AttributeError:
        ignore_control_signals()


def check_interrupted():
    if REQUESTED_SIGNAL is not None:
        raise RunInterrupted(REQUESTED_SIGNAL)


def leader_status_without_reaping(process, inspection_errors):
    try:
        info = os.waitid(
            os.P_PID,
            process.pid,
            os.WEXITED | os.WNOHANG | os.WNOWAIT,
        )
    except (ChildProcessError, OSError) as error:
        inspection_errors.append(
            f"leader waitid failed: {type(error).__name__}: {error}"
        )
        return None
    if info is None:
        return None
    if info.si_code == os.CLD_EXITED:
        return int(info.si_status)
    return -int(info.si_status)


def live_group_descendants(process, inspection_errors, observed_descendants):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    try:
        result = subprocess.run(
            ["/bin/ps", "-axo", "pid=,pgid=,stat="],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
            env=environment,
        )
    except (OSError, subprocess.SubprocessError) as error:
        inspection_errors.append(
            f"process inspection failed: {type(error).__name__}: {error}"
        )
        return None
    if result.returncode != 0:
        inspection_errors.append(
            f"process inspection exited {result.returncode}: "
            f"{result.stderr.strip()}"
        )
        return None
    live = set()
    for line in result.stdout.splitlines():
        fields = line.split(None, 2)
        if len(fields) != 3:
            continue
        pid_text, pgid_text, state = fields
        try:
            pid = int(pid_text)
            pgid = int(pgid_text)
        except ValueError:
            continue
        try:
            sid = os.getsid(pid)
        except ProcessLookupError:
            continue
        except PermissionError as error:
            inspection_errors.append(
                f"cannot inspect session for pid {pid}: "
                f"{type(error).__name__}: {error}"
            )
            return None
        if (
            pgid == process.pid
            and sid == process.pid
            and pid != process.pid
            and not state.startswith("Z")
        ):
            live.add(pid)
    observed_descendants.update(live)
    return live


def signal_owned_group(process, signum, inspection_errors):
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        return False
    except PermissionError as error:
        inspection_errors.append(
            f"group signal denied: {type(error).__name__}: {error}"
        )
        return False
    return True


def wait_for_drain(
    process,
    deadline_ns,
    inspection_errors,
    observed_descendants,
):
    leader_returncode = None
    while True:
        leader_returncode = leader_status_without_reaping(
            process, inspection_errors
        )
        live = live_group_descendants(
            process, inspection_errors, observed_descendants
        )
        if leader_returncode is not None and live == set():
            return True, leader_returncode
        now_ns = time.monotonic_ns()
        if now_ns >= deadline_ns:
            return False, leader_returncode
        time.sleep(min(0.05, (deadline_ns - now_ns) / 1_000_000_000))


def reap_if_exited(process, leader_returncode, inspection_errors):
    if leader_returncode is None:
        return False
    try:
        observed = process.wait(timeout=0.2)
    except subprocess.TimeoutExpired:
        inspection_errors.append(
            "leader was observable as exited but bounded reap timed out"
        )
        return False
    if observed != leader_returncode:
        inspection_errors.append(
            f"leader status changed between waitid and reap: "
            f"{leader_returncode} -> {observed}"
        )
        return False
    return True


def drain_owned_group(
    process,
    *,
    grace_signal,
    term_grace_seconds,
    kill_grace_seconds,
    inspection_errors,
    observed_descendants,
):
    leader_returncode = leader_status_without_reaping(process, inspection_errors)
    live = live_group_descendants(
        process, inspection_errors, observed_descendants
    )
    if leader_returncode is not None and live == set():
        reaped = reap_if_exited(process, leader_returncode, inspection_errors)
        return {
            "drain_status": "not_needed" if reaped and not inspection_errors else "incomplete",
            "leader_returncode": leader_returncode,
            "grace_signal_sent": False,
            "kill_signal_sent": False,
        }

    grace_sent = signal_owned_group(process, grace_signal, inspection_errors)
    complete, observed_returncode = wait_for_drain(
        process,
        time.monotonic_ns() + int(term_grace_seconds * 1_000_000_000),
        inspection_errors,
        observed_descendants,
    )
    kill_sent = False
    if not complete:
        kill_sent = signal_owned_group(process, signal.SIGKILL, inspection_errors)
        complete, observed_returncode = wait_for_drain(
            process,
            time.monotonic_ns() + int(kill_grace_seconds * 1_000_000_000),
            inspection_errors,
            observed_descendants,
        )
    reaped = reap_if_exited(process, observed_returncode, inspection_errors)
    complete = complete and reaped and not inspection_errors
    return {
        "drain_status": "complete" if complete else "incomplete",
        "leader_returncode": observed_returncode,
        "grace_signal_sent": grace_sent,
        "kill_signal_sent": kill_sent,
    }


def wait_for_leader_until(process, deadline_ns, inspection_errors):
    while True:
        check_interrupted()
        leader_returncode = leader_status_without_reaping(
            process, inspection_errors
        )
        if leader_returncode is not None:
            return leader_returncode
        now_ns = time.monotonic_ns()
        if now_ns >= deadline_ns:
            return None
        time.sleep(min(0.05, (deadline_ns - now_ns) / 1_000_000_000))


def leader_classification(returncode):
    if returncode is None:
        return {
            "leader_returncode": None,
            "leader_exit_code": None,
            "leader_signal": None,
            "leader_exit_kind": "unknown",
        }
    if returncode >= 0:
        return {
            "leader_returncode": returncode,
            "leader_exit_code": returncode,
            "leader_signal": None,
            "leader_exit_kind": "code",
        }
    signum = -returncode
    try:
        signal_name = signal.Signals(signum).name
    except ValueError:
        signal_name = f"SIG{signum}"
    return {
        "leader_returncode": returncode,
        "leader_exit_code": None,
        "leader_signal": signal_name,
        "leader_exit_kind": "signal",
    }


def log_identity(log_path):
    return {
        "log_path": str(log_path),
        "log_sha256": hash_file(log_path),
        "log_bytes": log_path.stat().st_size,
    }


def supervise_target(
    *,
    profile,
    target,
    command,
    deadline_ns,
    term_grace_seconds,
    kill_grace_seconds,
    log_path,
    source_identity,
):
    started_ns = time.monotonic_ns()
    inspection_errors = []
    observed_descendants = set()
    process = None
    log_handle = open(log_path, "xb", buffering=0)
    log_handle.write(
        canonical_bytes(
            {
                "schema": "frankensim-ascent-target-log-header-v1",
                "profile": profile,
                "target": target,
                "command": command,
            }
        )
    )
    try:
        check_interrupted()
        if time.monotonic_ns() >= deadline_ns:
            log_handle.close()
            receipt = {
                "schema": "frankensim-ascent-conformance-target-v2",
                "event": "target-result",
                "profile": profile,
                "package": "fs-ascent",
                "target": target,
                "command": command,
                "source": source_identity,
                "status": "budget_exceeded",
                "launched": False,
                "terminal_exit_code": 124,
                "elapsed_seconds": 0,
                "budget_status": "exceeded",
                "deadline_clock": "monotonic",
                "drain_status": "not_applicable",
                "drain_trigger": "deadline_before_spawn",
                "containment_scope": "new-session-process-group",
                "escaped_session_descendants_claimed": False,
                "process_group_identity_pinned_until_drain": False,
                "observed_descendant_count": 0,
                "inspection_errors": [],
                "forwarded_signal": None,
                "grace_signal_sent": False,
                "kill_signal_sent": False,
                **leader_classification(None),
                **log_identity(log_path),
            }
            return receipt, 124
        try:
            process = subprocess.Popen(
                command,
                cwd=REPO_ROOT,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        except OSError as error:
            log_handle.write(
                (
                    f"spawn failure: {type(error).__name__}: {error}\n"
                ).encode("utf-8", "backslashreplace")
            )
            log_handle.close()
            receipt = {
                "schema": "frankensim-ascent-conformance-target-v2",
                "event": "target-result",
                "profile": profile,
                "package": "fs-ascent",
                "target": target,
                "command": command,
                "source": source_identity,
                "status": "launch_failed",
                "launched": False,
                "terminal_exit_code": 126,
                "elapsed_seconds": max(
                    0, (time.monotonic_ns() - started_ns) // 1_000_000_000
                ),
                "budget_status": "within",
                "deadline_clock": "monotonic",
                "drain_status": "not_applicable",
                "drain_trigger": "spawn_failure",
                "containment_scope": "new-session-process-group",
                "escaped_session_descendants_claimed": False,
                "process_group_identity_pinned_until_drain": False,
                "observed_descendant_count": 0,
                "inspection_errors": [],
                "spawn_error": f"{type(error).__name__}: {error}",
                "forwarded_signal": None,
                "grace_signal_sent": False,
                "kill_signal_sent": False,
                **leader_classification(None),
                **log_identity(log_path),
            }
            return receipt, 126

        try:
            leader_returncode = wait_for_leader_until(
                process, deadline_ns, inspection_errors
            )
            if leader_returncode is None:
                drain = drain_owned_group(
                    process,
                    grace_signal=signal.SIGTERM,
                    term_grace_seconds=term_grace_seconds,
                    kill_grace_seconds=kill_grace_seconds,
                    inspection_errors=inspection_errors,
                    observed_descendants=observed_descendants,
                )
                wrapper_code = (
                    124
                    if drain["drain_status"] in ("complete", "not_needed")
                    else 125
                )
                status = "budget_exceeded"
                budget_status = "exceeded"
                drain_trigger = "deadline"
                forwarded_signal = None
            else:
                drain = drain_owned_group(
                    process,
                    grace_signal=signal.SIGTERM,
                    term_grace_seconds=term_grace_seconds,
                    kill_grace_seconds=kill_grace_seconds,
                    inspection_errors=inspection_errors,
                    observed_descendants=observed_descendants,
                )
                leader_returncode = drain["leader_returncode"]
                if drain["drain_status"] == "not_needed":
                    status = "pass" if leader_returncode == 0 else "fail"
                    wrapper_code = 0 if leader_returncode == 0 else 1
                    drain_trigger = "none"
                else:
                    status = "fail"
                    wrapper_code = (
                        1 if drain["drain_status"] == "complete" else 127
                    )
                    drain_trigger = "leader_exit_with_live_group"
                budget_status = "within"
                forwarded_signal = None
        except RunInterrupted as interruption:
            ignore_control_signals()
            drain = drain_owned_group(
                process,
                grace_signal=interruption.signum,
                term_grace_seconds=term_grace_seconds,
                kill_grace_seconds=kill_grace_seconds,
                inspection_errors=inspection_errors,
                observed_descendants=observed_descendants,
            )
            status = "interrupted"
            budget_status = "within"
            drain_trigger = (
                f"signal_{signal.Signals(interruption.signum).name}"
            )
            forwarded_signal = signal.Signals(interruption.signum).name
            wrapper_code = (
                128 + interruption.signum
                if drain["drain_status"] in ("complete", "not_needed")
                else 127
            )

        log_handle.close()
        leader_returncode = drain["leader_returncode"]
        receipt = {
            "schema": "frankensim-ascent-conformance-target-v2",
            "event": "target-result",
            "profile": profile,
            "package": "fs-ascent",
            "target": target,
            "command": command,
            "source": source_identity,
            "status": status,
            "launched": True,
            "terminal_exit_code": wrapper_code,
            "elapsed_seconds": max(
                0, (time.monotonic_ns() - started_ns) // 1_000_000_000
            ),
            "budget_status": budget_status,
            "deadline_clock": "monotonic",
            "drain_status": drain["drain_status"],
            "drain_trigger": drain_trigger,
            "containment_scope": "new-session-process-group",
            "escaped_session_descendants_claimed": False,
            "process_group_identity_pinned_until_drain": True,
            "observed_descendant_count": len(observed_descendants),
            "inspection_errors": sorted(set(inspection_errors)),
            "forwarded_signal": forwarded_signal,
            "grace_signal_sent": drain["grace_signal_sent"],
            "kill_signal_sent": drain["kill_signal_sent"],
            **leader_classification(leader_returncode),
            **log_identity(log_path),
        }
        return receipt, wrapper_code
    except BaseException:
        if process is not None:
            ignore_control_signals()
            drain_owned_group(
                process,
                grace_signal=signal.SIGTERM,
                term_grace_seconds=term_grace_seconds,
                kill_grace_seconds=kill_grace_seconds,
                inspection_errors=inspection_errors,
                observed_descendants=observed_descendants,
            )
        raise
    finally:
        if not log_handle.closed:
            log_handle.close()


def create_run_dir(profile):
    try:
        head_hint = (
            git_bytes("rev-parse", "--verify", "HEAD").strip().decode("ascii")[:12]
        )
    except Exception:
        head_hint = "unknown"
    log_root = pathlib.Path(
        os.environ.get(
            "FS_ASCENT_PROFILE_LOG_DIR",
            str(REPO_ROOT / "target" / "ascent-conformance-profile"),
        )
    )
    if not log_root.is_absolute():
        log_root = REPO_ROOT / log_root
    log_root.mkdir(parents=True, exist_ok=True)
    return pathlib.Path(
        tempfile.mkdtemp(prefix=f"{head_hint}-{profile}-", dir=log_root)
    ).resolve()


def safe_repo_identity():
    if os.environ.get("FS_ASCENT_INTERNAL_SELF_TEST") == "1":
        return {
            "head_sha": "self-test-head",
            "head_tree_sha": "self-test-head-tree",
            "index_sha256": "self-test-index",
            "root_tree_sha256": "self-test-root-tree",
            "git_status_sha256": "self-test-status",
            "dirty": False,
            "cargo_lock_sha256": "self-test-cargo-lock",
            "constellation_lock_sha256": "self-test-constellation-lock",
            "identity_fixture": "frankensim-ascent-self-test-source-v1",
        }, None
    try:
        return repo_identity(), None
    except Exception as error:
        return None, f"{type(error).__name__}: {error}"


def run_profile(profile):
    global REQUESTED_SIGNAL
    REQUESTED_SIGNAL = None
    install_signal_handlers()
    budget_seconds = budget_for_profile(profile)
    term_grace_seconds = parse_positive_seconds(
        "termination grace",
        os.environ.get("FS_ASCENT_TERMINATION_GRACE_SECONDS", "5"),
    )
    kill_grace_seconds = parse_positive_seconds(
        "kill drain",
        os.environ.get("FS_ASCENT_KILL_DRAIN_SECONDS", "5"),
    )
    run_dir = create_run_dir(profile)
    runner_argv = [str(SCRIPT_PATH), profile]
    writer = VerdictWriter(run_dir, profile, runner_argv)

    def incomplete_exit_guard():
        if not writer.sealed:
            writer.seal_once(
                status="internal_error",
                terminal_exit_code=127,
                provenance_state="incomplete",
                source_after=None,
                detail="process exited before the normal terminal seal",
            )
            writer.close()

    atexit.register(incomplete_exit_guard)
    source_before, identity_error = safe_repo_identity()
    writer.source_before = source_before
    if identity_error is not None:
        writer.append(
            {
                "schema": "frankensim-ascent-conformance-run-error-v1",
                "event": "run-error",
                "profile": profile,
                "status": "admission_failed",
                "detail": identity_error,
                "run_dir": str(run_dir),
            }
        )
        writer.seal_once(
            status="admission_failed",
            terminal_exit_code=127,
            provenance_state="incomplete",
            source_after=None,
            detail=identity_error,
        )
        writer.close()
        return 127

    started_ns = time.monotonic_ns()
    deadline_ns = started_ns + budget_seconds * 1_000_000_000
    targets = PR_TARGETS if profile == "pr" else ("all",)
    writer.append(
        {
            "schema": "frankensim-ascent-conformance-run-v2",
            "event": "run-start",
            "profile": profile,
            "package": "fs-ascent",
            "budget_seconds": budget_seconds,
            "deadline_clock": "monotonic",
            "build_time_included": True,
            "deadline_enforced": True,
            "total_targets": len(targets),
            "runner_argv": runner_argv,
            "run_dir": str(run_dir),
            "verdicts_path": str(writer.path),
            "child_output_channel": "per-target-log",
            "receipt_channel": "stdout-and-verdicts-jsonl",
            "containment_scope": "new-session-process-group",
            "escaped_session_descendants_claimed": False,
            "source": source_before,
        }
    )
    if (
        os.environ.get("FS_ASCENT_INTERNAL_SELF_TEST") == "1"
        and os.environ.get("FS_ASCENT_SELF_TEST_INJECT") == "after-run-start"
    ):
        raise RuntimeError("injected failure after run-start")

    attempted = 0
    passed = 0
    failed = 0
    budget_exceeded = 0
    run_status = "pass"
    terminal_exit_code = 0
    detail = "all selected targets passed"
    try:
        for index, target in enumerate(targets, start=1):
            check_interrupted()
            command = command_for_target(profile, target)
            log_path = run_dir / f"target-{index:02d}-{target}.log"
            receipt, wrapper_code = supervise_target(
                profile=profile,
                target=target,
                command=command,
                deadline_ns=deadline_ns,
                term_grace_seconds=term_grace_seconds,
                kill_grace_seconds=kill_grace_seconds,
                log_path=log_path,
                source_identity=source_before,
            )
            attempted += 1
            writer.append(receipt)
            if wrapper_code == 0:
                passed += 1
                continue
            if wrapper_code == 1:
                failed += 1
                run_status = "fail"
                terminal_exit_code = 1
                detail = "one or more selected targets failed"
                continue
            if wrapper_code in (123, 124, 125):
                budget_exceeded += 1
                run_status = "budget_exceeded"
                terminal_exit_code = 124 if wrapper_code in (123, 124) else 125
                detail = "aggregate monotonic deadline was exceeded"
                break
            if wrapper_code in (129, 130, 143):
                failed += 1
                run_status = "interrupted"
                terminal_exit_code = wrapper_code
                detail = f"runner received {receipt['forwarded_signal']}"
                break
            failed += 1
            run_status = (
                "launch_failed" if wrapper_code == 126 else "containment_failed"
            )
            terminal_exit_code = wrapper_code
            detail = receipt.get("spawn_error", "target containment failed")
            break
    except RunInterrupted as interruption:
        run_status = "interrupted"
        terminal_exit_code = 128 + interruption.signum
        detail = f"runner received {signal.Signals(interruption.signum).name}"

    source_after, final_identity_error = safe_repo_identity()
    close_control_signal_admission()
    if final_identity_error is not None or source_after != source_before:
        provenance_state = "incomplete"
    elif source_before.get("identity_fixture") is not None:
        provenance_state = "self_test"
    else:
        provenance_state = "sealed"
    if provenance_state == "incomplete":
        run_status = "provenance_failed"
        terminal_exit_code = 127
        detail = (
            final_identity_error
            or "HEAD, tree, index, status, or lock identity moved during the run"
        )
    elif REQUESTED_SIGNAL is not None and run_status != "interrupted":
        run_status = "interrupted"
        terminal_exit_code = 128 + REQUESTED_SIGNAL
        detail = f"runner received {signal.Signals(REQUESTED_SIGNAL).name}"
    unattempted = len(targets) - attempted
    writer.append(
        {
            "schema": "frankensim-ascent-conformance-run-v2",
            "event": "run-summary",
            "profile": profile,
            "status": run_status,
            "terminal_exit_code": terminal_exit_code,
            "budget_status": (
                "exceeded" if budget_exceeded else "within"
            ),
            "budget_seconds": budget_seconds,
            "elapsed_seconds": max(
                0, (time.monotonic_ns() - started_ns) // 1_000_000_000
            ),
            "total_targets": len(targets),
            "attempted_targets": attempted,
            "passed_targets": passed,
            "failed_targets": failed,
            "budget_exceeded_targets": budget_exceeded,
            "unattempted_targets": unattempted,
            "deadline_clock": "monotonic",
            "source_before": source_before,
            "source_after": source_after,
            "provenance_state": provenance_state,
            "detail": detail,
        }
    )
    writer.seal_once(
        status=run_status,
        terminal_exit_code=terminal_exit_code,
        provenance_state=provenance_state,
        source_after=source_after,
        detail=detail,
    )
    writer.close()
    return terminal_exit_code


FAKE_CARGO_SOURCE = r'''#!/usr/bin/env python3
import json
import os
import pathlib
import signal
import subprocess
import sys
import time

scenario = os.environ.get("FS_ASCENT_FAKE_SCENARIO", "pass")
counter = os.environ.get("FS_ASCENT_FAKE_COUNTER")
ready = os.environ.get("FS_ASCENT_FAKE_READY")
if counter:
    with open(counter, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(sys.argv, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
if scenario == "pass":
    print("FAKE_CHILD_STDOUT")
    print("FAKE_CHILD_STDERR", file=sys.stderr)
    raise SystemExit(0)
if scenario == "exit7":
    raise SystemExit(7)
if scenario == "grandchild":
    for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, signal.SIG_IGN)
    time.sleep(30)
    raise SystemExit(0)
if scenario == "ignore-tree":
    for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, signal.SIG_IGN)
    environment = os.environ.copy()
    environment["FS_ASCENT_FAKE_SCENARIO"] = "grandchild"
    environment.pop("FS_ASCENT_FAKE_COUNTER", None)
    environment.pop("FS_ASCENT_FAKE_READY", None)
    child = subprocess.Popen([sys.executable, __file__], env=environment)
    if ready:
        pathlib.Path(ready).write_text(
            json.dumps({"leader": os.getpid(), "grandchild": child.pid}),
            encoding="utf-8",
        )
    time.sleep(30)
    raise SystemExit(0)
raise SystemExit(f"unknown fake scenario: {scenario}")
'''


def parse_json_lines(text):
    rows = []
    for line in text.splitlines():
        if line:
            rows.append(json.loads(line))
    return rows


def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)


def verify_case_receipts(case):
    stdout_rows = parse_json_lines(case["stdout"])
    starts = [
        row
        for row in stdout_rows
        if row.get("event") == "run-start"
    ]
    assert_true(len(starts) == 1, f"{case['name']}: expected one run-start")
    run_dir = pathlib.Path(starts[0]["run_dir"])
    verdict_path = run_dir / "verdicts.jsonl"
    payload = verdict_path.read_bytes()
    lines = payload.splitlines(keepends=True)
    file_rows = [json.loads(line) for line in lines]
    seals = [row for row in file_rows if row.get("event") == "proof-seal"]
    assert_true(len(seals) == 1, f"{case['name']}: expected exactly one seal")
    assert_true(
        file_rows[-1].get("event") == "proof-seal",
        f"{case['name']}: seal was not terminal",
    )
    assert_true(
        hashlib.sha256(b"".join(lines[:-1])).hexdigest()
        == seals[0]["verdicts_prefix_sha256"],
        f"{case['name']}: prefix hash mismatch",
    )
    assert_true(
        stdout_rows == file_rows,
        f"{case['name']}: stdout mirror differs from authoritative JSONL",
    )
    case["rows"] = file_rows
    case["seal"] = seals[0]
    case["run_dir"] = run_dir
    return case


def bounded_self_test_process(command, *, environment, signal_to_send=None, ready=None):
    process = subprocess.Popen(
        command,
        cwd=REPO_ROOT,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    if signal_to_send is not None:
        deadline = time.monotonic() + 5
        while ready is not None and not ready.exists():
            if process.poll() is not None:
                break
            if time.monotonic() >= deadline:
                break
            time.sleep(0.02)
        if process.poll() is None:
            os.kill(process.pid, signal_to_send)
    try:
        stdout, stderr = process.communicate(timeout=15)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        if ready is not None and ready.exists():
            try:
                leader = json.loads(ready.read_text(encoding="utf-8"))["leader"]
                os.killpg(int(leader), signal.SIGKILL)
            except (OSError, KeyError, ValueError, json.JSONDecodeError):
                pass
        stdout, stderr = process.communicate(timeout=2)
        raise AssertionError(
            f"self-test runner exceeded bounded wait: {command!r}"
        )
    return process.returncode, stdout, stderr


def run_self_tests():
    root = pathlib.Path(
        tempfile.mkdtemp(
            prefix="frankensim ascent self test ",
            dir=os.environ.get("TMPDIR", "/tmp"),
        )
    ).resolve()
    fake_dir = root / "path with spaces"
    fake_dir.mkdir(parents=True)
    fake_cargo = fake_dir / "fake cargo.py"
    fake_cargo.write_text(FAKE_CARGO_SOURCE, encoding="utf-8")
    fake_cargo.chmod(0o755)
    cases = []

    def run_case(
        name,
        *,
        profile,
        scenario,
        expected_returncode,
        cargo_path=fake_cargo,
        budget=5,
        signal_to_send=None,
        inject=None,
    ):
        case_root = root / name
        case_root.mkdir()
        ready = case_root / "ready.json"
        counter = case_root / "launches.jsonl"
        environment = os.environ.copy()
        environment.update(
            {
                "PYTHON_BIN": sys.executable,
                "CARGO_BIN": str(cargo_path),
                "FS_ASCENT_PR_BUDGET_SECONDS": str(budget),
                "FS_ASCENT_NIGHTLY_BUDGET_SECONDS": str(budget),
                "FS_ASCENT_TERMINATION_GRACE_SECONDS": "0.10",
                "FS_ASCENT_KILL_DRAIN_SECONDS": "0.75",
                "FS_ASCENT_PROFILE_LOG_DIR": str(
                    case_root / "retained logs with spaces"
                ),
                "FS_ASCENT_FAKE_SCENARIO": scenario,
                "FS_ASCENT_FAKE_READY": str(ready),
                "FS_ASCENT_FAKE_COUNTER": str(counter),
                "FS_ASCENT_INTERNAL_SELF_TEST": "1",
            }
        )
        if inject is not None:
            environment["FS_ASCENT_SELF_TEST_INJECT"] = inject
        returncode, stdout, stderr = bounded_self_test_process(
            [str(SCRIPT_PATH), profile],
            environment=environment,
            signal_to_send=signal_to_send,
            ready=ready if signal_to_send is not None else None,
        )
        assert_true(
            returncode == expected_returncode,
            f"{name}: expected exit {expected_returncode}, got "
            f"{returncode}; stderr={stderr!r}",
        )
        case = verify_case_receipts(
            {
                "name": name,
                "returncode": returncode,
                "stdout": stdout,
                "stderr": stderr,
                "ready": ready,
                "counter": counter,
            }
        )
        cases.append(case)
        return case

    passed = run_case(
        "pass and path spaces",
        profile="nightly",
        scenario="pass",
        expected_returncode=0,
    )
    target_rows = [
        row for row in passed["rows"] if row.get("event") == "target-result"
    ]
    assert_true(len(target_rows) == 1, "pass: target receipt cardinality")
    assert_true(target_rows[0]["leader_exit_code"] == 0, "pass: exit mapping")
    assert_true(
        "FAKE_CHILD_STDOUT" not in passed["stdout"],
        "pass: child stdout contaminated receipt channel",
    )
    assert_true(
        "FAKE_CHILD_STDOUT"
        in pathlib.Path(target_rows[0]["log_path"]).read_text(
            encoding="utf-8"
        ),
        "pass: child stdout missing from retained target log",
    )
    assert_true(
        " " in str(passed["run_dir"]),
        "pass: run directory did not exercise path spaces",
    )
    assert_true(
        passed["seal"]["provenance_state"] == "self_test",
        "pass: self-test source identity was not labeled as a fixture",
    )

    exited = run_case(
        "exit mapping",
        profile="nightly",
        scenario="exit7",
        expected_returncode=1,
    )
    exit_target = next(
        row for row in exited["rows"] if row.get("event") == "target-result"
    )
    assert_true(
        exit_target["leader_exit_code"] == 7
        and exit_target["leader_exit_kind"] == "code",
        "exit mapping: exact child exit was not retained",
    )

    missing = root / "missing path with spaces" / "cargo"
    launch_failed = run_case(
        "launch failure",
        profile="nightly",
        scenario="pass",
        expected_returncode=126,
        cargo_path=missing,
    )
    launch_target = next(
        row
        for row in launch_failed["rows"]
        if row.get("event") == "target-result"
    )
    assert_true(
        launch_target["status"] == "launch_failed"
        and launch_target["launched"] is False,
        "launch failure: classification mismatch",
    )

    timed_out = run_case(
        "timeout child grandchild",
        profile="pr",
        scenario="ignore-tree",
        expected_returncode=124,
        budget=1,
    )
    timeout_targets = [
        row
        for row in timed_out["rows"]
        if row.get("event") == "target-result"
    ]
    assert_true(len(timeout_targets) == 1, "timeout: later target launched")
    timeout_target = timeout_targets[0]
    assert_true(
        timeout_target["status"] == "budget_exceeded"
        and timeout_target["drain_status"] == "complete"
        and timeout_target["kill_signal_sent"] is True
        and timeout_target["observed_descendant_count"] >= 1,
        "timeout: stubborn process tree was not bounded and drained",
    )
    launch_count = (
        len(timed_out["counter"].read_text(encoding="utf-8").splitlines())
        if timed_out["counter"].exists()
        else 0
    )
    assert_true(launch_count == 1, "timeout: later Cargo target was launched")

    for signum, expected in (
        (signal.SIGHUP, 129),
        (signal.SIGINT, 130),
        (signal.SIGTERM, 143),
    ):
        interrupted = run_case(
            f"forward {signal.Signals(signum).name}",
            profile="nightly",
            scenario="ignore-tree",
            expected_returncode=expected,
            budget=30,
            signal_to_send=signum,
        )
        interrupt_target = next(
            row
            for row in interrupted["rows"]
            if row.get("event") == "target-result"
        )
        assert_true(
            interrupt_target["status"] == "interrupted"
            and interrupt_target["forwarded_signal"]
            == signal.Signals(signum).name
            and interrupt_target["drain_status"] == "complete",
            f"signal forwarding failed for {signal.Signals(signum).name}",
        )

    injected = run_case(
        "exit guard seal",
        profile="nightly",
        scenario="pass",
        expected_returncode=127,
        inject="after-run-start",
    )
    assert_true(
        injected["seal"]["provenance_state"] == "incomplete"
        and injected["seal"]["status"] == "internal_error",
        "exit guard: incomplete terminal seal missing",
    )

    emit_stdout(
        {
            "schema": "frankensim-ascent-conformance-self-test-v2",
            "status": "pass",
            "cases": len(cases),
            "cargo_invocations": 0,
            "retained_artifact_root": str(root),
            "temporary_files_deleted": 0,
        }
    )
    return 0


def main():
    if MODE == "list":
        return list_profile(PROFILE)
    if MODE == "self-test":
        try:
            return run_self_tests()
        except BaseException as error:
            emit_stdout(
                {
                    "schema": "frankensim-ascent-conformance-self-test-v2",
                    "status": "fail",
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc().splitlines(),
                    "cargo_invocations": 0,
                }
            )
            return 1
    if MODE != "run":
        raise ValueError(f"unknown mode {MODE!r}")
    try:
        return run_profile(PROFILE)
    except BaseException as error:
        print(
            f"ascent profile internal error: {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        return 127


try:
    raise SystemExit(main())
except ValueError as error:
    emit_stdout(
        {
            "schema": "frankensim-ascent-conformance-invocation-error-v1",
            "status": "fail",
            "error": str(error),
        }
    )
    raise SystemExit(2)
PY
}

if (( $# == 2 )) && [[ "$1" == "--list" ]]; then
  run_python_supervisor list "$2"
fi

if (( $# == 1 )) && [[ "$1" == "--self-test" ]]; then
  run_python_supervisor self-test ""
fi

if (( $# == 1 )) && [[ "$1" == "pr" || "$1" == "nightly" ]]; then
  run_python_supervisor run "$1"
fi

usage
exit 2
