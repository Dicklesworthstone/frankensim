//! Workspace-suite green receipt (bead frankensim-fluz9).
//!
//! Why this exists. Nothing in the repository establishes that the workspace
//! test suite passes; the README's test-file count is inventory, not a green
//! claim, and L3 cannot be admitted over a suite whose green/red state is
//! unknown. This module produces the retained, dated, reproducible receipt:
//! exact command, toolchain, host fingerprint, HEAD and constellation-lock
//! identity, per-crate pass/fail/ignored counts, the totals, and an EXPLICIT
//! known-red set whose every entry carries an owning bead and a disposition.
//!
//! How the run is made (no format guessing, no coverage laundering):
//!
//! 1. `cargo metadata` enumerates every test target of every native-workspace
//!    package and of the attested fs-wasm nested workspace.
//! 2. `cargo test --workspace --no-run --message-format json --keep-going`
//!    builds every test target it can; targets with no emitted executable
//!    are recorded as BUILD FAILURES by name (a broken integration-test
//!    target must not silently reduce the receipt's coverage).
//! 3. Each built test executable is run directly with `-Z unstable-options
//!    --format json`, cwd at its package root — exactly how cargo itself
//!    runs it — giving per-test events with certain package attribution.
//!    Doc tests run per package through `cargo test --doc` for the same
//!    reason. Executables run sequentially, in sorted order, so the receipt
//!    is a deterministic function of the same run inputs.
//! 4. Failures are partitioned by the tracked `suite-known-red.json`
//!    registry: every entry carries an owning bead and a disposition;
//!    anything else is UNEXPECTED RED. The receipt refuses to claim
//!    `green` or `green-with-known-red` while unexpected red exists — the
//!    negative control the bead names as its test.
//!
//! The check gate validates the receipt's schema and internal coherence,
//! refuses green claims over red facts, refuses a registered known-red whose
//! owner bead is closed while the test still fails, and renders staleness
//! (HEAD moved since the run) as a visible note, never a wedge.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation, fnv1a64};

pub(crate) const CHECK: &str = "suite-receipt";
const RECEIPT_PATH: &str = "suite-receipt.json";
const REGISTRY_PATH: &str = "suite-known-red.json";
const ISSUES_PATH: &str = ".beads/issues.jsonl";
const SCHEMA: &str = "frankensim-suite-receipt-v2";
const REGISTRY_SCHEMA: &str = "frankensim-suite-known-red-v1";
const MAX_RUN_BYTES: u64 = 512 * 1024 * 1024;
const MAX_JSON_STRING: usize = 8 * 1024 * 1024;
/// One test executable may not run longer than this before it is recorded
/// as a timeout rather than allowed to hang the receipt forever.
const TARGET_TIMEOUT_SECS: u64 = 900;
/// Browser flagships are intentionally a separate Cargo workspace, but their
/// native test result is still part of the repository suite receipt.
const ATTESTED_NESTED_WORKSPACES: &[&str] = &["crates/fs-wasm"];

/// A failing or passing-again registered known-red test.
#[derive(Debug, Clone, PartialEq, Eq)]
struct KnownRedEntry {
    test: String,
    krate: String,
    workspace: String,
    owner_bead: String,
    disposition: String,
}

/// Outcome of one test target (one executable).
#[derive(Debug, Clone, Default)]
struct TargetOutcome {
    passed: usize,
    failed: usize,
    ignored: usize,
    failures: Vec<String>,
    /// Set when the target could not run at all (crash, timeout, no JSON).
    target_error: Option<String>,
}

/// The full run model, pre-render.
#[derive(Debug)]
struct RunModel {
    command: String,
    executed_at: String,
    host: String,
    toolchain: String,
    head_sha: String,
    head_dirty: bool,
    lock_hash: String,
    /// Stable workspace identity -> that workspace's Cargo.lock fingerprint.
    /// A nested workspace can resolve a different graph than the native root.
    workspace_locks: BTreeMap<String, String>,
    target_triple: String,
    crates: BTreeMap<String, TargetOutcome>,
    build_failures: Vec<String>,
    known_red: Vec<(KnownRedEntry, bool)>, // entry, observed_failed
    unexpected_red: Vec<(String, String)>, // (crate, test)
    excluded: Vec<String>,
}

impl RunModel {
    fn totals(&self) -> (usize, usize, usize) {
        let mut totals = (0usize, 0usize, 0usize);
        for outcome in self.crates.values() {
            totals.0 += outcome.passed;
            totals.1 += outcome.failed;
            totals.2 += outcome.ignored;
        }
        totals
    }

    fn status(&self) -> &'static str {
        let (_, failed, _) = self.totals();
        if !self.build_failures.is_empty()
            || !self.unexpected_red.is_empty()
            || self
                .crates
                .values()
                .any(|outcome| outcome.target_error.is_some())
        {
            "not-green"
        } else if failed == 0 {
            "green"
        } else {
            "green-with-known-red"
        }
    }
}

fn violation(detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: RECEIPT_PATH.to_string(),
        detail: detail.into(),
    }
}

fn json_obj_field<'a>(map: &'a BTreeMap<String, JsonValue>, key: &str) -> Option<&'a JsonValue> {
    map.get(key)
}

fn json_string(map: &BTreeMap<String, JsonValue>, key: &str) -> Option<String> {
    match map.get(key) {
        Some(JsonValue::String(value)) => Some(value.clone()),
        _ => None,
    }
}

fn json_count(map: &BTreeMap<String, JsonValue>, key: &str) -> Result<usize, String> {
    match map.get(key) {
        Some(JsonValue::Number(raw)) => raw
            .parse::<usize>()
            .map_err(|error| format!("`{key}` is not a count: {error}")),
        _ => Err(format!("missing count `{key}`")),
    }
}

fn parse_json_lines(text: &str, context: &str) -> Vec<BTreeMap<String, JsonValue>> {
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if !line.starts_with('{') {
            continue;
        }
        if let Ok(JsonValue::Object(map)) =
            JsonParser::with_string_limit(line, MAX_JSON_STRING).finish()
        {
            out.push(map);
        } else {
            // A non-JSON line inside a JSON channel is noise, not evidence;
            // skip it rather than fabricate a parse failure for output that
            // may carry panic text.
            let _ = context;
        }
    }
    out
}

/// Load and validate the known-red registry.
fn load_registry(root: &Path) -> Result<Vec<KnownRedEntry>, String> {
    let text = std::fs::read_to_string(root.join(REGISTRY_PATH))
        .map_err(|error| format!("{REGISTRY_PATH} unreadable: {error}"))?;
    let parsed = JsonParser::with_string_limit(&text, MAX_JSON_STRING)
        .finish()
        .map_err(|error| format!("{REGISTRY_PATH} is not valid JSON: {error}"))?;
    let JsonValue::Object(map) = &parsed else {
        return Err(format!("{REGISTRY_PATH} is not a JSON object"));
    };
    match map.get("schema") {
        Some(JsonValue::String(schema)) if schema == REGISTRY_SCHEMA => {}
        _ => {
            return Err(format!(
                "{REGISTRY_PATH} schema is not `{REGISTRY_SCHEMA}`; refusing to read a \
                 foreign registry"
            ));
        }
    }
    let mut entries = Vec::new();
    match map.get("known_red") {
        Some(JsonValue::Array(rows)) => {
            for row in rows {
                let JsonValue::Object(row) = row else {
                    return Err(format!("{REGISTRY_PATH} has a non-object entry"));
                };
                let get = |key: &str| {
                    json_string(row, key)
                        .ok_or_else(|| format!("{REGISTRY_PATH} entry is missing `{key}`"))
                };
                let disposition = get("disposition")?;
                if disposition != "blocked-upstream" && disposition != "repaired" {
                    return Err(format!(
                        "{REGISTRY_PATH} entry has unknown disposition `{disposition}`; \
                         dispositions are blocked-upstream or repaired"
                    ));
                }
                entries.push(KnownRedEntry {
                    test: get("test")?,
                    krate: get("crate")?,
                    workspace: json_string(row, "workspace")
                        .unwrap_or_else(|| "native".to_string()),
                    owner_bead: get("owner_bead")?,
                    disposition,
                });
            }
        }
        _ => return Err(format!("{REGISTRY_PATH} has no `known_red` array")),
    }
    Ok(entries)
}

/// Which bead ids are closed in the tracker (for owner-liveness checks).
fn closed_bead_ids(root: &Path) -> Result<std::collections::BTreeSet<String>, String> {
    let text = std::fs::read_to_string(root.join(ISSUES_PATH))
        .map_err(|error| format!("{ISSUES_PATH} unreadable: {error}"))?;
    let mut closed = std::collections::BTreeSet::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let Ok(JsonValue::Object(map)) =
            JsonParser::with_string_limit(line, MAX_JSON_STRING).finish()
        else {
            continue;
        };
        if let (Some(id), Some(status)) = (json_string(&map, "id"), json_string(&map, "status"))
            && status == "closed"
        {
            closed.insert(id);
        }
    }
    Ok(closed)
}

fn run_command(command: &str, args: &[&str], cwd: &Path) -> Result<(bool, String, String), String> {
    let output = Command::new(command)
        .args(args)
        .current_dir(cwd)
        .env("CARGO_TERM_COLOR", "never")
        // The receipt's cargo children must use the HOST's own default target
        // dir. An inherited CARGO_TARGET_DIR travels badly across machines
        // (a Mac NVMe path on a Linux worker refuses every target build),
        // and the target dir is ephemeral build state, never receipt content.
        .env_remove("CARGO_TARGET_DIR")
        .env_remove("CARGO_BUILD_TARGET_DIR")
        .output()
        .map_err(|error| format!("cannot spawn {command} {args:?}: {error}"))?;
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    Ok((output.status.success(), stdout, stderr))
}

/// The cargo package name owning an artifact, from its package_id.
///
/// Cargo emits three shapes: `name version (source)`, `name@version`, and
/// URL specs — `path+file:///dir/name#version` for path packages and
/// `registry+url#name@version` for registry packages. All three resolve to
/// the bare package name; anything unrecognizable is returned whole rather
/// than silently renamed.
fn package_name(package_id: &str) -> String {
    let token = package_id.split_whitespace().next().unwrap_or(package_id);
    if token.contains("+file://") || token.contains("+http") {
        let (source, fragment) = match token.split_once('#') {
            Some((source, fragment)) => (source, Some(fragment)),
            None => (token, None),
        };
        // `registry+url#name@version`: the name is the fragment, but only
        // when the fragment actually carries a name@ (a bare version, as in
        // `path+file:///dir/name#version`, is not a name).
        if let Some(fragment) = fragment
            && fragment.contains('@')
            && let Some(name) = fragment.split('@').next()
            && !name.is_empty()
        {
            return name.to_string();
        }
        // `path+file:///dir/name#version`: the name is the last path segment.
        let name = source.rsplit('/').next().unwrap_or(source);
        return name.split('@').next().unwrap_or(name).to_string();
    }
    token.split('@').next().unwrap_or(token).to_string()
}

/// Enumerate workspace test targets via cargo metadata, separating
/// feature-gated targets (which `cargo test --workspace` itself skips under
/// the default feature set) from the expected set. A gated target is not a
/// build failure; it is recorded as excluded, on the record.
fn metadata_targets(root: &Path) -> Result<(Vec<(String, String, PathBuf)>, Vec<String>), String> {
    let (ok, stdout, stderr) = run_command(
        "cargo",
        &["metadata", "--no-deps", "--format-version", "1"],
        root,
    )?;
    if !ok {
        return Err(format!("cargo metadata failed: {}", stderr.trim()));
    }
    let parsed = JsonParser::with_string_limit(&stdout, MAX_JSON_STRING * 16)
        .finish()
        .map_err(|error| format!("cargo metadata is not valid JSON: {error}"))?;
    let JsonValue::Object(map) = &parsed else {
        return Err("cargo metadata is not an object".to_string());
    };
    let mut targets = Vec::new();
    let mut feature_excluded = Vec::new();
    match json_obj_field(map, "packages") {
        Some(JsonValue::Array(packages)) => {
            for package in packages {
                let JsonValue::Object(package) = package else {
                    continue;
                };
                let name = json_string(package, "name")
                    .ok_or_else(|| "metadata package without name".to_string())?;
                let manifest = json_string(package, "manifest_path")
                    .ok_or_else(|| format!("metadata package {name} without manifest_path"))?;
                let package_root = Path::new(&manifest)
                    .parent()
                    .map(Path::to_path_buf)
                    .ok_or_else(|| format!("metadata package {name} manifest has no parent"))?;
                if let Some(JsonValue::Array(targets_json)) = json_obj_field(package, "targets") {
                    for target in targets_json {
                        let JsonValue::Object(target) = target else {
                            continue;
                        };
                        // `test: true` marks targets compiled in test mode,
                        // but examples and build scripts also carry it:
                        // cargo test does not RUN those as tests. Only lib
                        // (unit tests), test (integration), and bin (binary
                        // unit tests) targets produce runnable test
                        // executables; anything else is enumeration noise.
                        let runnable_kind = match json_obj_field(target, "kind") {
                            Some(JsonValue::Array(kinds)) => kinds.iter().any(|kind| {
                                matches!(kind, JsonValue::String(k) if k == "lib" || k == "test" || k == "bin")
                            }),
                            _ => false,
                        };
                        if !runnable_kind {
                            continue;
                        }
                        if !matches!(json_obj_field(target, "test"), Some(JsonValue::Bool(true))) {
                            continue;
                        }
                        let target_name = json_string(target, "name")
                            .ok_or_else(|| format!("metadata {name} target without name"))?;
                        if let Some(JsonValue::Array(required)) =
                            json_obj_field(target, "required-features")
                        {
                            let features: Vec<String> = required
                                .iter()
                                .filter_map(|feature| match feature {
                                    JsonValue::String(feature) => Some(feature.clone()),
                                    _ => None,
                                })
                                .collect();
                            if !features.is_empty() {
                                feature_excluded.push(format!(
                                    "{name}::{target_name} (feature-gated: {})",
                                    features.join(",")
                                ));
                                continue;
                            }
                        }
                        targets.push((name.clone(), target_name, package_root.clone()));
                    }
                }
            }
        }
        _ => return Err("cargo metadata has no packages array".to_string()),
    }
    if targets.is_empty() {
        return Err("cargo metadata enumerated zero test targets".to_string());
    }
    Ok((targets, feature_excluded))
}

/// Run the workspace suite and build the model. This is the long pole: a
/// full test-target build plus every test executable, sequentially.
fn crate_key(workspace: &str, package: &str) -> String {
    format!("{workspace}::{package}")
}

fn run_workspace_targets(
    cargo: &str,
    workspace_id: &str,
    workspace_root: &Path,
    targets: Vec<(String, String, PathBuf)>,
    crates: &mut BTreeMap<String, TargetOutcome>,
    build_failures: &mut Vec<String>,
) -> Result<(), String> {
    let build_failures_before = build_failures.len();
    let (build_ok, build_stdout, build_stderr) = run_command(
        cargo,
        &[
            "build",
            "--tests",
            "--workspace",
            "--message-format",
            "json",
            "--keep-going",
        ],
        workspace_root,
    )?;
    // Map (package, target-name) -> executable for successfully built targets.
    let mut executables: BTreeMap<(String, String), (PathBuf, PathBuf)> = BTreeMap::new();
    let mut package_roots: BTreeMap<String, PathBuf> = BTreeMap::new();
    for (name, _, package_root) in &targets {
        package_roots.insert(name.clone(), package_root.clone());
    }
    for map in parse_json_lines(&build_stdout, "compiler-artifact") {
        if json_string(&map, "reason").as_deref() != Some("compiler-artifact") {
            continue;
        }
        let Some(JsonValue::Object(target)) = json_obj_field(&map, "target") else {
            continue;
        };
        let Some(executable) = json_string(&map, "executable") else {
            continue;
        };
        let Some(package_id) = json_string(&map, "package_id") else {
            continue;
        };
        let package = package_name(&package_id);
        let Some(target_name) = json_string(target, "name") else {
            continue;
        };
        let package_root = package_roots
            .get(&package)
            .cloned()
            .unwrap_or_else(|| workspace_root.to_path_buf());
        executables.insert(
            (package, target_name),
            (PathBuf::from(executable), package_root),
        );
    }
    if !build_ok {
        for (package, target_name, _) in &targets {
            // Doc-test targets have no executable artifact by construction;
            // they run through rustdoc, not the linker.
            if !executables.contains_key(&(package.clone(), target_name.clone())) {
                build_failures.push(format!(
                    "{} (test target {target_name})",
                    crate_key(workspace_id, package)
                ));
            }
        }
        if build_failures.len() == build_failures_before {
            build_failures.push(format!(
                "unidentified target; cargo said: {}",
                stderr_tail(&build_stderr)
            ));
        }
    }
    let mut sorted_executables: Vec<_> = executables.into_iter().collect();
    sorted_executables.sort();
    for ((package, target_name), (executable, package_root)) in sorted_executables {
        let key = crate_key(workspace_id, &package);
        eprintln!("suite-receipt: running {key}::{target_name}");
        let outcome = crates.entry(key).or_default();
        run_test_executable(&executable, &package_root, outcome);
    }
    // Doc tests per package (rustdoc has no executable artifact).
    for (package, package_root) in &package_roots {
        let key = crate_key(workspace_id, package);
        eprintln!("suite-receipt: running {key} doc tests");
        let outcome = crates.entry(key).or_default();
        if cargo == "cargo" {
            run_doc_tests(workspace_root, package, package_root, outcome);
        } else {
            run_doc_tests_with_cargo(cargo, workspace_root, package, package_root, outcome);
        }
    }
    Ok(())
}

/// Run the native workspace plus each explicitly attested nested workspace.
/// A nested target-build failure is recorded as a build failure, never an
/// exclusion; metadata admission failures instead prevent a receipt.
pub(crate) fn run_suite(root: &Path) -> Result<RunModel, String> {
    let (targets, mut feature_excluded) = metadata_targets(root)?;
    let mut crates = BTreeMap::new();
    let mut build_failures = Vec::new();
    let mut workspace_locks = BTreeMap::new();
    workspace_locks.insert("native".to_string(), cargo_lock_hash(root));
    run_workspace_targets(
        "cargo",
        "native",
        root,
        targets,
        &mut crates,
        &mut build_failures,
    )?;

    for nested in ATTESTED_NESTED_WORKSPACES {
        let nested_root = root.join(nested);
        workspace_locks.insert((*nested).to_string(), cargo_lock_hash(&nested_root));
        let (targets, nested_excluded) = metadata_targets(&nested_root)
            .map_err(|error| format!("nested workspace {nested}: {error}"))?;
        feature_excluded.extend(
            nested_excluded
                .into_iter()
                .map(|target| format!("{nested}/{target}")),
        );
        run_workspace_targets(
            "cargo",
            nested,
            &nested_root,
            targets,
            &mut crates,
            &mut build_failures,
        )
        .map_err(|error| format!("nested workspace {nested}: {error}"))?;
    }
    if crates.is_empty() {
        return Err("suite receipt recorded no runnable crate outcomes".to_string());
    }

    Ok(RunModel {
        command: "cargo build --tests --workspace --message-format json --keep-going in the \
                  native workspace and each attested nested workspace; then each test executable \
                  with -Z unstable-options --format json"
            .to_string(),
        executed_at: utc_now(),
        host: host_fingerprint(),
        toolchain: toolchain_string(root),
        head_sha: head_sha(root),
        head_dirty: head_dirty(root),
        lock_hash: constellation_lock_hash(root),
        workspace_locks,
        target_triple: target_triple(),
        crates,
        build_failures,
        known_red: Vec::new(),
        unexpected_red: Vec::new(),
        excluded: feature_excluded,
    })
}

fn stderr_tail(stderr: &str) -> String {
    stderr
        .lines()
        .filter(|line| line.contains("error"))
        .take(3)
        .collect::<Vec<_>>()
        .join(" | ")
}

/// Run one test executable with JSON output into the crate outcome.
///
/// The suite must survive a hung target unattended (measured live: a storm
/// lane deadlocked all-threads-futex for 30+ minutes, bead frankensim-kh5tf).
/// std offers no child timeout, so a watchdog thread reaps the child after
/// [`TARGET_TIMEOUT_SECS`]; the outcome records the kill as a timeout, never
/// as a pass and never as a hang.
fn run_test_executable(executable: &Path, cwd: &Path, outcome: &mut TargetOutcome) {
    let child = Command::new(executable)
        .args(["-Z", "unstable-options", "--format", "json"])
        .current_dir(cwd)
        .env("CARGO_TERM_COLOR", "never")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn();
    let child = match child {
        Ok(child) => child,
        Err(error) => {
            outcome.target_error = Some(format!("cannot spawn {}: {error}", executable.display()));
            return;
        }
    };
    let pid = child.id();
    let watchdog = std::thread::spawn(move || {
        std::thread::sleep(std::time::Duration::from_secs(TARGET_TIMEOUT_SECS));
        // std has no portable kill-by-pid; the project is unix-only, and a
        // kill on an already-reaped pid is a harmless ignored error.
        let _ = Command::new("kill").arg(pid.to_string()).status();
    });
    let start = std::time::Instant::now();
    match child.wait_with_output() {
        Ok(output) => {
            let elapsed = start.elapsed().as_secs();
            let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
            let mut saw_suite = false;
            for map in parse_json_lines(&stdout, "libtest") {
                match json_string(&map, "type").as_deref() {
                    Some("test") => match json_string(&map, "event").as_deref() {
                        Some("ok") => outcome.passed += 1,
                        Some("failed") => {
                            outcome.failed += 1;
                            if let Some(name) = json_string(&map, "name") {
                                outcome.failures.push(name);
                            }
                        }
                        Some("ignored") => outcome.ignored += 1,
                        _ => {}
                    },
                    Some("suite") => saw_suite = true,
                    _ => {}
                }
            }
            if !saw_suite {
                outcome.target_error = Some(format!(
                    "no suite summary in output (exit {:?})",
                    output.status.code()
                ));
            }
            if elapsed >= TARGET_TIMEOUT_SECS {
                outcome.target_error = Some(format!(
                    "target exceeded the {TARGET_TIMEOUT_SECS}s bound (killed at {elapsed}s)"
                ));
            }
        }
        Err(error) => {
            outcome.target_error = Some(format!("cannot await {}: {error}", executable.display()));
        }
    }
    drop(watchdog);
}

fn run_doc_tests(root: &Path, package: &str, package_root: &Path, outcome: &mut TargetOutcome) {
    run_doc_tests_with_cargo("cargo", root, package, package_root, outcome);
}

fn run_doc_tests_with_cargo(
    cargo: &str,
    root: &Path,
    package: &str,
    package_root: &Path,
    outcome: &mut TargetOutcome,
) {
    let args = [
        "test",
        "--doc",
        "-p",
        package,
        "--",
        "-Z",
        "unstable-options",
        "--format",
        "json",
    ];
    let result = run_command(cargo, &args, root);
    match result {
        Ok((ok, stdout, stderr)) => {
            for map in parse_json_lines(&stdout, "rustdoc") {
                if json_string(&map, "type").as_deref() == Some("test") {
                    match json_string(&map, "event").as_deref() {
                        Some("ok") => outcome.passed += 1,
                        Some("failed") => {
                            outcome.failed += 1;
                            if let Some(name) = json_string(&map, "name") {
                                outcome.failures.push(format!("(doc) {name}"));
                            }
                        }
                        Some("ignored") => outcome.ignored += 1,
                        _ => {}
                    }
                }
            }
            if !ok {
                outcome.target_error = Some(format!(
                    "doc tests for {package} at {} failed: {}",
                    package_root.display(),
                    stderr_tail(&stderr)
                ));
            }
        }
        Err(error) => {
            outcome.target_error = Some(format!(
                "doc tests for {package} at {} failed to spawn: {error}",
                package_root.display()
            ));
        }
    }
}

fn utc_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_secs())
        .unwrap_or(0);
    // Civil-from-days, stable and dependency-free.
    let days = secs / 86_400;
    let secs_of_day = secs % 86_400;
    let (year, month, day) = civil_from_days(days as i64);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        secs_of_day / 3600,
        (secs_of_day % 3600) / 60,
        secs_of_day % 60
    )
}

fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    (
        if month <= 2 { year + 1 } else { year },
        month as u32,
        day as u32,
    )
}

fn host_fingerprint() -> String {
    let uname = run_command("uname", &["-srm"], Path::new("/")).map_or_else(
        |_| "unknown-uname".to_string(),
        |(_, out, _)| out.trim().to_string(),
    );
    let cpus = run_command("nproc", &[], Path::new("/"))
        .ok()
        .and_then(|(ok, out, _)| ok.then(|| out.trim().to_string()))
        .or_else(|| {
            run_command("sysctl", &["-n", "hw.ncpu"], Path::new("/"))
                .ok()
                .and_then(|(ok, out, _)| ok.then(|| out.trim().to_string()))
        })
        .unwrap_or_else(|| "?".to_string());
    format!("{uname}, {cpus} cpus")
}

fn toolchain_string(root: &Path) -> String {
    let channel = std::fs::read_to_string(root.join("rust-toolchain.toml"))
        .ok()
        .and_then(|text| {
            text.lines()
                .find(|line| line.trim_start().starts_with("channel"))
                .map(str::to_string)
        })
        .unwrap_or_else(|| "channel unpinned?".to_string());
    let rustc = run_command("rustc", &["--version"], root).map_or_else(
        |_| "rustc unknown".to_string(),
        |(_, out, _)| out.trim().to_string(),
    );
    format!("{}; {rustc}", channel.trim())
}

fn head_sha(root: &Path) -> String {
    crate::git_out(root, &["rev-parse", "HEAD"]).unwrap_or_else(|_| "no-git".to_string())
}

fn head_dirty(root: &Path) -> bool {
    crate::git_out(root, &["status", "--porcelain"]).map_or(true, |status| !status.is_empty())
}

fn constellation_lock_hash(root: &Path) -> String {
    std::fs::read(root.join("constellation.lock"))
        .map(|bytes| format!("{:016x}", fnv1a64(&bytes)))
        .unwrap_or_else(|_| "no-lock".to_string())
}

fn cargo_lock_hash(workspace_root: &Path) -> String {
    std::fs::read(workspace_root.join("Cargo.lock"))
        .map(|bytes| format!("{:016x}", fnv1a64(&bytes)))
        .unwrap_or_else(|_| "no-cargo-lock".to_string())
}

fn target_triple() -> String {
    run_command("rustc", &["-vV"], Path::new("/"))
        .ok()
        .and_then(|(_, out, _)| {
            out.lines()
                .find(|line| line.starts_with("host:"))
                .map(|line| line.trim_start_matches("host:").trim().to_string())
        })
        .unwrap_or_else(|| "unknown-target".to_string())
}

/// Attach the known-red partition to a completed run.
fn partition_failures(root: &Path, model: &mut RunModel) -> Result<(), String> {
    let registry = load_registry(root)?;
    let mut known_red = Vec::new();
    let mut unexpected = Vec::new();
    for entry in &registry {
        let observed = model
            .crates
            .get(&crate_key(&entry.workspace, &entry.krate))
            .is_some_and(|outcome| outcome.failures.iter().any(|name| name == &entry.test));
        known_red.push((entry.clone(), observed));
    }
    for (krate, outcome) in &model.crates {
        for failure in &outcome.failures {
            let registered = registry.iter().any(|entry| {
                crate_key(&entry.workspace, &entry.krate) == *krate && entry.test == *failure
            });
            if !registered {
                unexpected.push((krate.clone(), failure.clone()));
            }
        }
    }
    model.known_red = known_red;
    model.unexpected_red = unexpected;
    Ok(())
}

fn escape(text: &str) -> String {
    crate::json_escape(text)
}

/// Render the canonical receipt bytes.
fn render(model: &RunModel) -> String {
    let (passed, failed, ignored) = model.totals();
    let mut out = String::new();
    out.push_str("{\n");
    writeln!(out, "  \"schema\": \"{SCHEMA}\"").ok();
    out.push_str("  \"bead\": \"frankensim-fluz9\",\n");
    out.push_str("  \"run\": {\n");
    writeln!(
        out,
        "    \"command\": \"{}\",\n    \"executed_at\": \"{}\",\n    \"host_fingerprint\": \"{}\",\n    \
         \"toolchain\": \"{}\",\n    \"head_sha\": \"{}\",\n    \"head_dirty\": {},\n    \
         \"constellation_lock_fnv1a64\": \"{}\",\n    \"target\": \"{}\"\n  }},",
        escape(&model.command),
        model.executed_at,
        escape(&model.host),
        escape(&model.toolchain),
        model.head_sha,
        model.head_dirty,
        model.lock_hash,
        model.target_triple
    )
    .ok();
    out.push_str("  \"workspaces\": {\n");
    for (index, (workspace, lock_hash)) in model.workspace_locks.iter().enumerate() {
        let comma = if index + 1 == model.workspace_locks.len() {
            ""
        } else {
            ","
        };
        writeln!(
            out,
            "    \"{}\": {{\"cargo_lock_fnv1a64\": \"{}\"}}{comma}",
            escape(workspace),
            escape(lock_hash)
        )
        .ok();
    }
    out.push_str("  },\n");
    writeln!(out, "  \"status\": \"{}\",", model.status()).ok();
    writeln!(
        out,
        "  \"totals\": {{\"passed\": {passed}, \"failed\": {failed}, \"ignored\": {ignored}}},"
    )
    .ok();
    out.push_str("  \"crates\": {\n");
    let crate_count = model.crates.len();
    for (index, (krate, outcome)) in model.crates.iter().enumerate() {
        let comma = if index + 1 == crate_count { "" } else { "," };
        let mut line = format!(
            "    \"{krate}\": {{\"passed\": {}, \"failed\": {}, \"ignored\": {}",
            outcome.passed, outcome.failed, outcome.ignored
        );
        if outcome.target_error.is_some() {
            line.push_str(", \"target_error\": true");
        }
        writeln!(line, "}}{comma}").ok();
        out.push_str(&line);
    }
    out.push_str("  },\n");
    out.push_str("  \"build_failures\": [");
    for (index, failure) in model.build_failures.iter().enumerate() {
        if index > 0 {
            out.push_str(", ");
        }
        write!(out, "\"{}\"", escape(failure)).ok();
    }
    out.push_str("],\n");
    out.push_str("  \"known_red\": [\n");
    for (index, (entry, observed)) in model.known_red.iter().enumerate() {
        let comma = if index + 1 == model.known_red.len() {
            ""
        } else {
            ","
        };
        writeln!(
            out,
            "    {{\"test\": \"{}\", \"crate\": \"{}\", \"workspace\": \"{}\", \"owner_bead\": \"{}\", \"disposition\": \
             \"{}\", \"observed_failed\": {observed}}}{comma}",
            escape(&entry.test),
            escape(&entry.krate),
            escape(&entry.workspace),
            escape(&entry.owner_bead),
            entry.disposition
        )
        .ok();
    }
    out.push_str("  ],\n");
    out.push_str("  \"unexpected_red\": [\n");
    for (index, (krate, test)) in model.unexpected_red.iter().enumerate() {
        let comma = if index + 1 == model.unexpected_red.len() {
            ""
        } else {
            ","
        };
        writeln!(
            out,
            "    {{\"crate\": \"{}\", \"test\": \"{}\"}}{comma}",
            escape(krate),
            escape(test)
        )
        .ok();
    }
    out.push_str("  ],\n");
    out.push_str("  \"excluded\": [");
    for (index, excluded) in model.excluded.iter().enumerate() {
        if index > 0 {
            out.push_str(", ");
        }
        write!(out, "\"{}\"", escape(excluded)).ok();
    }
    out.push_str("],\n");
    out.push_str(
        "  \"no_claim\": \"one run on one host against the recorded HEAD and lock state; green \
         means these tests passed here, not that the suite is exhaustive, flake-free on other \
         hosts, or that known-red tests are anything but tracked upstream debt\"\n",
    );
    out.push_str("}\n");
    out
}

/// The validated receipt, for the gate.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SuiteReceipt {
    status: String,
    head_sha: String,
    passed: usize,
    failed: usize,
    ignored: usize,
    crate_count: usize,
    build_failures: usize,
    known_red: Vec<(KnownRedEntry, bool)>,
    unexpected_red: usize,
    /// Crates whose run ended with a target error (watchdog timeout or a
    /// missing suite summary): a known-red test on one of these crates is
    /// UNKNOWN, never "repaired".
    target_error_crates: std::collections::BTreeSet<String>,
}

fn parse_receipt(text: &str) -> Result<SuiteReceipt, String> {
    let parsed = JsonParser::with_string_limit(text, MAX_JSON_STRING)
        .finish()
        .map_err(|error| format!("{RECEIPT_PATH} is not valid JSON: {error}"))?;
    let JsonValue::Object(map) = &parsed else {
        return Err(format!("{RECEIPT_PATH} is not a JSON object"));
    };
    match map.get("schema") {
        Some(JsonValue::String(schema)) if schema == SCHEMA => {}
        Some(JsonValue::String(schema)) => {
            return Err(format!(
                "{RECEIPT_PATH} schema is `{schema}`, expected `{SCHEMA}`"
            ));
        }
        _ => return Err(format!("{RECEIPT_PATH} has no schema string")),
    }
    let status =
        json_string(map, "status").ok_or_else(|| format!("{RECEIPT_PATH} has no status"))?;
    if !matches!(
        status.as_str(),
        "green" | "green-with-known-red" | "not-green"
    ) {
        return Err(format!("{RECEIPT_PATH} has unknown status `{status}`"));
    }
    let Some(JsonValue::Object(run)) = json_obj_field(map, "run") else {
        return Err(format!("{RECEIPT_PATH} has no run object"));
    };
    let head_sha = json_string(run, "head_sha")
        .ok_or_else(|| format!("{RECEIPT_PATH} run has no head_sha"))?;
    let Some(JsonValue::Object(totals)) = json_obj_field(map, "totals") else {
        return Err(format!("{RECEIPT_PATH} has no totals object"));
    };
    let passed = json_count(totals, "passed")?;
    let failed = json_count(totals, "failed")?;
    let ignored = json_count(totals, "ignored")?;
    let mut crate_count = 0usize;
    let mut crate_sum = (0usize, 0usize, 0usize);
    let mut target_error_crates = std::collections::BTreeSet::new();
    match json_obj_field(map, "crates") {
        Some(JsonValue::Object(crates)) => {
            for (krate, value) in crates {
                let JsonValue::Object(value) = value else {
                    return Err(format!("{RECEIPT_PATH} crate row is not an object"));
                };
                crate_count += 1;
                crate_sum.0 += json_count(value, "passed")?;
                crate_sum.1 += json_count(value, "failed")?;
                crate_sum.2 += json_count(value, "ignored")?;
                if matches!(
                    json_obj_field(value, "target_error"),
                    Some(JsonValue::Bool(true))
                ) {
                    target_error_crates.insert(krate.clone());
                }
            }
        }
        _ => return Err(format!("{RECEIPT_PATH} has no crates object")),
    }
    if crate_sum != (passed, failed, ignored) {
        return Err(format!(
            "{RECEIPT_PATH} is internally inconsistent: totals ({passed}/{failed}/{ignored}) \
             do not equal the crate sums ({}/{}/{})",
            crate_sum.0, crate_sum.1, crate_sum.2
        ));
    }
    let build_failures = match json_obj_field(map, "build_failures") {
        Some(JsonValue::Array(failures)) => failures.len(),
        _ => return Err(format!("{RECEIPT_PATH} has no build_failures array")),
    };
    let mut known_red = Vec::new();
    match json_obj_field(map, "known_red") {
        Some(JsonValue::Array(rows)) => {
            for row in rows {
                let JsonValue::Object(row) = row else {
                    return Err(format!("{RECEIPT_PATH} known_red row is not an object"));
                };
                let observed = match json_obj_field(row, "observed_failed") {
                    Some(JsonValue::Bool(observed)) => *observed,
                    _ => return Err("known_red row without observed_failed boolean".to_string()),
                };
                known_red.push((
                    KnownRedEntry {
                        test: json_string(row, "test")
                            .ok_or_else(|| "known_red row without test".to_string())?,
                        krate: json_string(row, "crate")
                            .ok_or_else(|| "known_red row without crate".to_string())?,
                        workspace: json_string(row, "workspace")
                            .unwrap_or_else(|| "native".to_string()),
                        owner_bead: json_string(row, "owner_bead")
                            .ok_or_else(|| "known_red row without owner_bead".to_string())?,
                        disposition: json_string(row, "disposition")
                            .ok_or_else(|| "known_red row without disposition".to_string())?,
                    },
                    observed,
                ));
            }
        }
        _ => return Err(format!("{RECEIPT_PATH} has no known_red array")),
    }
    let unexpected_red = match json_obj_field(map, "unexpected_red") {
        Some(JsonValue::Array(rows)) => rows.len(),
        _ => return Err(format!("{RECEIPT_PATH} has no unexpected_red array")),
    };
    // The negative control, as a parse-time law: a receipt that claims any
    // green flavor while carrying unexpected red or build failures is
    // self-contradicting and must not validate.
    if unexpected_red > 0 && status != "not-green" {
        return Err(format!(
            "{RECEIPT_PATH} claims `{status}` while carrying {unexpected_red} unexpected red \
             test(s); a receipt may not launder unowned failures"
        ));
    }
    if build_failures > 0 && status != "not-green" {
        return Err(format!(
            "{RECEIPT_PATH} claims `{status}` while carrying {build_failures} build failure(s)"
        ));
    }
    if !target_error_crates.is_empty() && status != "not-green" {
        return Err(format!(
            "{RECEIPT_PATH} claims `{status}` while carrying target errors for {} crate(s)",
            target_error_crates.len()
        ));
    }
    let Some(JsonValue::Object(workspaces)) = json_obj_field(map, "workspaces") else {
        return Err(format!("{RECEIPT_PATH} has no workspaces object"));
    };
    let Some(JsonValue::Object(fs_wasm)) = workspaces.get("crates/fs-wasm") else {
        return Err(format!(
            "{RECEIPT_PATH} omits the attested crates/fs-wasm workspace"
        ));
    };
    match json_string(fs_wasm, "cargo_lock_fnv1a64") {
        Some(lock) if lock != "no-cargo-lock" => {}
        _ => {
            return Err(format!(
                "{RECEIPT_PATH} crates/fs-wasm workspace has no Cargo.lock identity"
            ));
        }
    }
    Ok(SuiteReceipt {
        status,
        head_sha,
        passed,
        failed,
        ignored,
        crate_count,
        build_failures,
        known_red,
        unexpected_red,
        target_error_crates,
    })
}

/// Generate the receipt. Long: builds every test target, then runs every
/// test executable sequentially.
pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let mut model = run_suite(root)?;
    partition_failures(root, &mut model)?;
    std::fs::write(root.join(RECEIPT_PATH), render(&model))
        .map_err(|error| format!("cannot write {RECEIPT_PATH}: {error}"))
}

/// The standing gate: validate the receipt, cross-check the known-red
/// registry against the tracker, and render staleness as a note.
pub(crate) fn check(root: &Path) -> (Vec<Violation>, Vec<PolicyNote>) {
    let text = match std::fs::read_to_string(root.join(RECEIPT_PATH)) {
        Ok(text) => text,
        Err(error) => {
            return (
                vec![violation(format!(
                    "{RECEIPT_PATH} unreadable: {error}; a missing receipt is not a green \
                     suite — run the documented regeneration command"
                ))],
                Vec::new(),
            );
        }
    };
    let receipt = match parse_receipt(&text) {
        Ok(receipt) => receipt,
        Err(error) => return (vec![violation(error)], Vec::new()),
    };
    let mut violations = Vec::new();
    let mut notes = Vec::new();
    // The receipt and the registry must cover the same tests: an entry in
    // one but not the other is a disposition drifting out of sync.
    match load_registry(root) {
        Ok(registry) => {
            for entry in &registry {
                let recorded = receipt.known_red.iter().any(|(receipt_entry, _)| {
                    receipt_entry.test == entry.test
                        && receipt_entry.krate == entry.krate
                        && receipt_entry.workspace == entry.workspace
                });
                if !recorded {
                    violations.push(violation(format!(
                        "registry entry {}::{} is absent from the receipt; either the suite \
                         regenerated without it or the registry drifted — reconcile deliberately",
                        entry.krate, entry.test
                    )));
                }
            }
            // An owner bead that closed while its test still fails is a
            // contradiction: closure says resolved, the suite says red.
            if let Ok(closed) = closed_bead_ids(root) {
                for (entry, observed) in &receipt.known_red {
                    if *observed && closed.contains(&entry.owner_bead) {
                        violations.push(violation(format!(
                            "known-red test {}::{} still fails but its owner bead {} is \
                             CLOSED; a closed owner with a red test repeats the es6pt error",
                            entry.krate, entry.test, entry.owner_bead
                        )));
                    }
                }
            }
        }
        Err(error) => violations.push(violation(error)),
    }
    // A registered test that passes again is repair evidence: visible, and
    // the registry row must leave deliberately. But a crate that ended in a
    // target error (watchdog timeout, no suite summary) proves NOTHING about
    // its known-red tests — "unknown", never "repaired".
    for (entry, observed) in &receipt.known_red {
        if receipt
            .target_error_crates
            .contains(&crate_key(&entry.workspace, &entry.krate))
        {
            notes.push(PolicyNote {
                check: CHECK,
                crate_name: RECEIPT_PATH.to_string(),
                verdict: "known-red-unknown",
                detail: format!(
                    "{}::{} state is UNKNOWN: the crate's run ended in a target error \
                     (watchdog timeout or no suite summary); not a pass, not a failure (owner {})",
                    entry.krate, entry.test, entry.owner_bead
                ),
            });
        } else if !observed {
            notes.push(PolicyNote {
                check: CHECK,
                crate_name: RECEIPT_PATH.to_string(),
                verdict: "known-red-repaired",
                detail: format!(
                    "{}::{} no longer fails; remove its registry row deliberately (owner {})",
                    entry.krate, entry.test, entry.owner_bead
                ),
            });
        }
    }
    let current_head = head_sha(root);
    if current_head != "no-git" && current_head != receipt.head_sha {
        notes.push(PolicyNote {
            check: CHECK,
            crate_name: RECEIPT_PATH.to_string(),
            verdict: "stale-receipt",
            detail: format!(
                "HEAD moved since the suite ran (receipt {}, live {current_head}); the \
                 receipt is a point-in-time fact — regenerate when the suite state should be \
                 re-attested",
                receipt.head_sha
            ),
        });
    }
    notes.push(PolicyNote {
        check: CHECK,
        crate_name: RECEIPT_PATH.to_string(),
        verdict: "receipt-status",
        detail: format!(
            "status {}: {}/{} passed/failed, {} ignored over {} crates, {} build failure(s), \
             {} unexpected red",
            receipt.status,
            receipt.passed,
            receipt.failed,
            receipt.ignored,
            receipt.crate_count,
            receipt.build_failures,
            receipt.unexpected_red
        ),
    });
    (violations, notes)
}

#[cfg(test)]
mod tests {
    //! G0: partition laws, the negative control, registry validation, and
    //! coherence refusals. The receipt's whole point is that it cannot be
    //! talked into claiming green over red facts.

    use super::*;

    fn registry_text(entries: &[(&str, &str, &str)]) -> String {
        let rows: Vec<String> = entries
            .iter()
            .map(|(krate, test, owner)| {
                format!(
                    "{{\"test\":\"{test}\",\"crate\":\"{krate}\",\"owner_bead\":\"{owner}\",\"disposition\":\"blocked-upstream\"}}"
                )
            })
            .collect();
        format!(
            "{{\"schema\":\"{REGISTRY_SCHEMA}\",\"known_red\":[{}]}}",
            rows.join(",")
        )
    }

    #[test]
    fn g0_civil_from_days_matches_known_dates() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(19_744), (2024, 1, 22));
        assert_eq!(civil_from_days(20_368), (2025, 10, 7));
    }

    #[test]
    fn g0_registry_refuses_foreign_schema_and_bad_disposition() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root")
            .to_path_buf();
        let scratch =
            std::env::temp_dir().join(format!("fsim-suite-receipt-test-{}", std::process::id()));
        std::fs::create_dir_all(&scratch).expect("scratch");
        std::fs::write(scratch.join(REGISTRY_PATH), "{\"schema\":\"other\"}").expect("write");
        assert!(load_registry(&scratch).is_err());
        std::fs::write(
            scratch.join(REGISTRY_PATH),
            "{\"schema\":\"frankensim-suite-known-red-v1\",\"known_red\":[{\"test\":\"t\",\"crate\":\"c\",\"owner_bead\":\"o\",\"disposition\":\"invented\"}]}",
        )
        .expect("write");
        let error = load_registry(&scratch).expect_err("bad disposition must refuse");
        assert!(error.contains("disposition"), "{error}");
        drop(root);
        let _ = std::fs::remove_dir_all(&scratch);
    }

    #[test]
    fn g0_receipt_refuses_totals_that_do_not_match_crate_sums() {
        let receipt = "{
  \"schema\": \"frankensim-suite-receipt-v2\",
  \"run\": {\"head_sha\": \"abc\"},
  \"status\": \"green\",
  \"totals\": {\"passed\": 2, \"failed\": 0, \"ignored\": 0},
  \"crates\": {\"fs-x\": {\"passed\": 1, \"failed\": 0, \"ignored\": 0}},
  \"build_failures\": [],
  \"known_red\": [],
  \"unexpected_red\": []
}\n";
        let error = parse_receipt(receipt).expect_err("sum mismatch must refuse");
        assert!(error.contains("inconsistent"), "{error}");
    }

    #[test]
    fn g0_the_negative_control_no_green_over_unexpected_red() {
        let receipt = "{
  \"schema\": \"frankensim-suite-receipt-v2\",
  \"run\": {\"head_sha\": \"abc\"},
  \"status\": \"green-with-known-red\",
  \"totals\": {\"passed\": 1, \"failed\": 1, \"ignored\": 0},
  \"crates\": {\"fs-x\": {\"passed\": 1, \"failed\": 1, \"ignored\": 0}},
  \"build_failures\": [],
  \"known_red\": [],
  \"unexpected_red\": [{\"crate\": \"fs-x\", \"test\": \"t\"}],
  \"workspaces\": {\"crates/fs-wasm\": {\"cargo_lock_fnv1a64\": \"nested-lock\"}}
}\n";
        let error = parse_receipt(receipt).expect_err("green over unexpected red must refuse");
        assert!(error.contains("launder"), "{error}");
    }

    #[test]
    fn g0_no_green_over_build_failures_either() {
        let receipt = "{
  \"schema\": \"frankensim-suite-receipt-v2\",
  \"run\": {\"head_sha\": \"abc\"},
  \"status\": \"green\",
  \"totals\": {\"passed\": 1, \"failed\": 0, \"ignored\": 0},
  \"crates\": {\"fs-x\": {\"passed\": 1, \"failed\": 0, \"ignored\": 0}},
  \"build_failures\": [\"fs-y (test target z)\"],
  \"known_red\": [],
  \"unexpected_red\": []
}\n";
        let error = parse_receipt(receipt).expect_err("green over build failure must refuse");
        assert!(error.contains("build failure"), "{error}");
    }

    #[test]
    fn g0_target_error_forces_not_green_status() {
        let model = RunModel {
            command: "test".to_string(),
            executed_at: "now".to_string(),
            host: "host".to_string(),
            toolchain: "toolchain".to_string(),
            head_sha: "head".to_string(),
            head_dirty: false,
            lock_hash: "lock".to_string(),
            workspace_locks: BTreeMap::new(),
            target_triple: "target".to_string(),
            crates: BTreeMap::from([(
                "native::fs-x".to_string(),
                TargetOutcome {
                    target_error: Some("no suite summary".to_string()),
                    ..TargetOutcome::default()
                },
            )]),
            build_failures: Vec::new(),
            known_red: Vec::new(),
            unexpected_red: Vec::new(),
            excluded: Vec::new(),
        };

        assert_eq!(model.status(), "not-green");
        let rendered = render(&model);
        assert!(rendered.contains("\"status\": \"not-green\""));
        assert!(rendered.contains("\"target_error\": true"));
    }

    #[test]
    fn g0_receipt_refuses_green_over_target_error() {
        let receipt = "{
  \"schema\": \"frankensim-suite-receipt-v2\",
  \"run\": {\"head_sha\": \"abc\"},
  \"status\": \"green\",
  \"totals\": {\"passed\": 0, \"failed\": 0, \"ignored\": 0},
  \"crates\": {\"native::fs-x\": {\"passed\": 0, \"failed\": 0, \"ignored\": 0, \"target_error\": true}},
  \"build_failures\": [],
  \"known_red\": [],
  \"unexpected_red\": [],
  \"workspaces\": {\"crates/fs-wasm\": {\"cargo_lock_fnv1a64\": \"nested-lock\"}}
}\n";

        let error = parse_receipt(receipt).expect_err("green over target error must refuse");
        assert!(error.contains("target errors"), "{error}");
    }

    #[test]
    fn g0_honest_not_green_receipt_validates() {
        let receipt = "{
  \"schema\": \"frankensim-suite-receipt-v2\",
  \"run\": {\"head_sha\": \"abc\"},
  \"status\": \"not-green\",
  \"totals\": {\"passed\": 1, \"failed\": 1, \"ignored\": 0},
  \"crates\": {\"fs-x\": {\"passed\": 1, \"failed\": 1, \"ignored\": 0}},
  \"build_failures\": [],
  \"known_red\": [],
  \"unexpected_red\": [{\"crate\": \"fs-x\", \"test\": \"t\"}],
  \"workspaces\": {\"crates/fs-wasm\": {\"cargo_lock_fnv1a64\": \"nested-lock\"}}
}\n";
        let receipt = parse_receipt(receipt).expect("an honest not-green receipt validates");
        assert_eq!(receipt.status, "not-green");
        assert_eq!(receipt.unexpected_red, 1);
    }

    #[test]
    fn g0_receipt_preserves_observed_failed_false() {
        let receipt = "{
  \"schema\": \"frankensim-suite-receipt-v2\",
  \"run\": {\"head_sha\": \"abc\"},
  \"status\": \"green\",
  \"totals\": {\"passed\": 1, \"failed\": 0, \"ignored\": 0},
  \"crates\": {\"native::fs-x\": {\"passed\": 1, \"failed\": 0, \"ignored\": 0}},
  \"build_failures\": [],
  \"known_red\": [{\"test\": \"t\", \"crate\": \"fs-x\", \"workspace\": \"native\", \"owner_bead\": \"b\", \"disposition\": \"blocked-upstream\", \"observed_failed\": false}],
  \"unexpected_red\": [],
  \"workspaces\": {\"crates/fs-wasm\": {\"cargo_lock_fnv1a64\": \"nested-lock\"}}
}\n";
        let receipt = parse_receipt(receipt).expect("receipt parses");
        assert!(!receipt.known_red[0].1, "false must not mean observed");
    }

    #[test]
    fn g0_package_name_parses_both_id_shapes() {
        assert_eq!(
            package_name("fs-ledger 0.1.0 (path+file:///x)"),
            "fs-ledger"
        );
        assert_eq!(package_name("fs-ledger@0.1.0"), "fs-ledger");
        assert_eq!(
            package_name("path+file:///data/projects/frankensim/crates/fs-cli#0.0.1"),
            "fs-cli"
        );
        assert_eq!(
            package_name("registry+https://github.com/rust-lang/crates.io-index#serde@1.0.0"),
            "serde"
        );
    }

    #[test]
    fn g0_fs_wasm_nested_workspace_is_attested() {
        assert_eq!(ATTESTED_NESTED_WORKSPACES, ["crates/fs-wasm"]);
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root")
            .to_path_buf();
        let nested = root.join("crates/fs-wasm");
        let (targets, _) = metadata_targets(&nested).expect("fs-wasm metadata");
        assert!(
            targets
                .iter()
                .any(|(krate, target, _)| krate == "fs-wasm" && target == "fs_wasm"),
            "the receipt must enumerate the fs-wasm test target"
        );
    }

    #[cfg(unix)]
    #[test]
    fn g0_nested_workspace_target_survives_render_and_known_red_partition() {
        use std::os::unix::fs::PermissionsExt;

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let scratch = std::env::temp_dir().join(format!("fsim-receipt-{nonce}"));
        std::fs::create_dir_all(&scratch).expect("scratch");
        let executable = scratch.join("fs-wasm-test");
        std::fs::write(
            &executable,
            "#!/bin/sh\nprintf '%s\\n' '{\"type\":\"test\",\"event\":\"failed\",\"name\":\"campaigns::tests::sensorforge_defaults\"}' '{\"type\":\"suite\",\"event\":\"ok\"}'\n",
        )
        .expect("test executable");
        let mut permissions = std::fs::metadata(&executable)
            .expect("metadata")
            .permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(&executable, permissions).expect("executable mode");

        let cargo = scratch.join("cargo");
        std::fs::write(
            &cargo,
            format!(
                "#!/bin/sh\nif [ \"$1\" = build ]; then printf '%s\\n' '{{\"reason\":\"compiler-artifact\",\"package_id\":\"fs-wasm@0.0.1\",\"target\":{{\"name\":\"fs_wasm\"}},\"executable\":\"{}\"}}'; elif [ \"$1\" = test ]; then echo doc-test-failure >&2; exit 1; fi\n",
                executable.display()
            ),
        )
        .expect("cargo shim");
        let mut permissions = std::fs::metadata(&cargo).expect("metadata").permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(&cargo, permissions).expect("cargo mode");
        std::fs::write(
            scratch.join(REGISTRY_PATH),
            registry_text(&[("fs-wasm", "campaigns::tests::sensorforge_defaults", "vxnd3")])
                .replace(
                    "\"crate\":\"fs-wasm\"",
                    "\"crate\":\"fs-wasm\",\"workspace\":\"crates/fs-wasm\"",
                ),
        )
        .expect("registry");

        let mut crates = BTreeMap::new();
        let mut build_failures = Vec::new();
        run_workspace_targets(
            cargo.to_str().expect("utf8 path"),
            "crates/fs-wasm",
            &scratch,
            vec![(
                "fs-wasm".to_string(),
                "fs_wasm".to_string(),
                scratch.clone(),
            )],
            &mut crates,
            &mut build_failures,
        )
        .expect("target run");
        assert!(
            crates
                .get("crates/fs-wasm::fs-wasm")
                .and_then(|outcome| outcome.target_error.as_ref())
                .is_some(),
            "a failing cargo doc-test command is a target error"
        );
        // The same package name in the native workspace must not overwrite
        // the exercised nested-workspace result.
        crates.insert(
            crate_key("native", "fs-wasm"),
            TargetOutcome {
                passed: 1,
                ..TargetOutcome::default()
            },
        );
        assert_eq!(crates.len(), 2, "workspace identity prevents key collision");
        let mut model = RunModel {
            command: "test".to_string(),
            executed_at: "now".to_string(),
            host: "host".to_string(),
            toolchain: "toolchain".to_string(),
            head_sha: "head".to_string(),
            head_dirty: false,
            lock_hash: "lock".to_string(),
            workspace_locks: BTreeMap::from([
                ("native".to_string(), "native-lock".to_string()),
                ("crates/fs-wasm".to_string(), "nested-lock".to_string()),
            ]),
            target_triple: "target".to_string(),
            crates,
            build_failures,
            known_red: Vec::new(),
            unexpected_red: Vec::new(),
            excluded: Vec::new(),
        };
        partition_failures(&scratch, &mut model).expect("partition");
        assert!(model.known_red[0].1, "the fs-wasm known red was observed");
        let rendered = render(&model);
        assert!(
            parse_receipt(&rendered).is_ok(),
            "nested lock identity validates"
        );
        assert!(rendered.contains("\"native::fs-wasm\""));
        assert!(rendered.contains("crates/fs-wasm::fs-wasm"));
        assert!(rendered.contains("\"crates/fs-wasm\": {\"cargo_lock_fnv1a64\": \"nested-lock\"}"));
    }

    #[test]
    fn g0_registry_text_helper_parses() {
        let scratch =
            std::env::temp_dir().join(format!("fsim-suite-registry-{}", std::process::id()));
        std::fs::create_dir_all(&scratch).expect("scratch");
        std::fs::write(
            scratch.join(REGISTRY_PATH),
            registry_text(&[("fs-ledger", "t1", "frankensim-x")]),
        )
        .expect("write");
        let entries = load_registry(&scratch).expect("registry parses");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].test, "t1");
        let _ = std::fs::remove_dir_all(&scratch);
    }
}
