#![cfg_attr(windows, feature(windows_by_handle))]

//! FrankenSim repository policy checks (`cargo run -p xtask -- <command>`).
//!
//! Commands:
//! - `check-layers`   — enforce the L0..L6 layer dependency direction (plan §4, AGENTS.md).
//! - `check-deps`     — enforce the Franken-only runtime dependency policy (Decalogue P1).
//! - `check-contracts`— every workspace `fs-*` crate ships a CONTRACT.md with required sections.
//! - `check-unsafe`   — unsafe code only in registered capsules (<300 lines, SAFETY.md).
//! - `check-powi`     — no build-mode-dependent `f64::powi` in deterministic paths (bead 4xnt).
//! - `check-libm`     — cross-ISA-claiming crates route transcendentals via fs_math::det (bead lyms).
//! - `check-color-admission` — no new positive Color literals outside admission authorities (bead 6pf9).
//! - `check-obs-events` — committed *.events.jsonl fixtures are schema-valid fs-obs lines (bead huq.16).
//! - `check-casual-print` — core libraries cannot create untyped stdout/stderr truth paths (bead i94v.7.3.3).
//! - `check-terminology` — enforce the repository's sole-branch vocabulary policy.
//! - `check-goldens`  — golden hashes declare upstream couplings; drift re-freezes deliberately (bead y4pt).
//! - `check-identities` — identity schemas classify fields and link mutation coverage (bead iu5l).
//! - `check-manifest-fixture` — admit only declared new-domain Cargo edges and an acyclic same-layer order.
//! - `check-claims`   — README hashes/crates/sentinels must exist in code (bead 06yc).
//! - `check-closures` — closed bug beads must cite regression evidence or a disposition (bead hx4p).
//! - `check-citable-producers` — exhaustively inventory authority-gated `citation_eligible` sinks.
//! - `check-all`      — all of the above; non-zero exit on any violation.
//! - `lock-constellation` / `check-constellation` — pin/verify the Franken library states.
//! - `matdb-pack`     — compile licensed material TSV or NASA-9 sources into normalized packs.
//!
//! Output is JSON-lines (one verdict object per check per crate) so agents parse
//! outcomes without scraping; a human-readable summary goes to stderr.
//!
//! Parsing note: this tool intentionally hand-parses the *subset* of TOML that this
//! repository's own generated manifests use (see docs/CONVENTIONS.md "Manifest
//! conventions": one dependency per line, `[section]` headers on their own line).
//! It fails loudly on shapes it does not understand rather than guessing — these are
//! our files, and the conventions are enforced, not inferred.

mod bootstrap_provenance;
mod claims;
mod closures;
pub mod constellation_admission;
mod constellation_cleanliness;
mod depgraph;
mod identities;
mod manifest_fixture;
mod matdb_pack;

use bootstrap_provenance::{
    BootstrapProvenanceRow, bootstrap_provenance_support_preflight, provenance_path_text,
    write_bootstrap_provenance,
};
use constellation_cleanliness::{
    is_redirecting_entry, pinned_repository_worktree_status, repository_worktree_status,
    sanitized_git_command, verify_two_complete_passes,
};

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::io::{Read as _, Write as _};
use std::path::{Component, Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_CONSTELLATION_LOCK_TEMP_SUFFIX: AtomicU64 = AtomicU64::new(0);

/// Layers, in the plan's order. `Util` crates (fs-qty, fs-obs) are usable by every
/// layer; `Tool` crates (xtask) are outside the shipped dependency graph entirely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layer {
    Util,
    L0,
    L1,
    L2,
    L3,
    L4,
    L5,
    L6,
    Tool,
}

impl Layer {
    fn parse(s: &str) -> Option<Layer> {
        Some(match s {
            "UTIL" => Layer::Util,
            "L0" => Layer::L0,
            "L1" => Layer::L1,
            "L2" => Layer::L2,
            "L3" => Layer::L3,
            "L4" => Layer::L4,
            "L5" => Layer::L5,
            "L6" => Layer::L6,
            "TOOL" => Layer::Tool,
            _ => return None,
        })
    }

    fn name(self) -> &'static str {
        match self {
            Layer::Util => "UTIL",
            Layer::L0 => "L0",
            Layer::L1 => "L1",
            Layer::L2 => "L2",
            Layer::L3 => "L3",
            Layer::L4 => "L4",
            Layer::L5 => "L5",
            Layer::L6 => "L6",
            Layer::Tool => "TOOL",
        }
    }

    /// May a crate in layer `self` depend on a crate in layer `dep`?
    ///
    /// This is NOT a linear order: ASCENT (L4) and LUMEN (L5) are siblings —
    /// LUMEN may not depend on ASCENT and vice versa (AGENTS.md dependency
    /// direction). HELM (L6) may depend on everything; everything may depend
    /// on UTIL; TOOL may depend on anything but nothing may depend on TOOL.
    fn may_depend_on(self, dep: Layer) -> bool {
        use Layer::{L0, L1, L2, L3, L4, L5, L6, Tool, Util};
        if dep == Tool {
            return false;
        }
        match self {
            Util => matches!(dep, Util),
            L0 => matches!(dep, Util | L0),
            L1 => matches!(dep, Util | L0 | L1),
            L2 => matches!(dep, Util | L0 | L1 | L2),
            L3 => matches!(dep, Util | L0 | L1 | L2 | L3),
            L4 => matches!(dep, Util | L0 | L1 | L2 | L3 | L4),
            // LUMEN: MORPH, FLUX field abstractions, BEDROCK, SUBSTRATE — not ASCENT.
            L5 => matches!(dep, Util | L0 | L1 | L2 | L3 | L5),
            L6 | Tool => true,
        }
    }
}

/// The Franken constellation: the ONLY permitted external runtime dependencies
/// besides `std` (Decalogue P1). The libraries are workspaces whose member
/// crates use these prefixes (probed from the actual sibling repos):
/// asupersync (single crate), fsqlite* (FrankenSQLite), fnx-* (FrankenNetworkx),
/// fnp-* (FrankenNumpy), ft-* (FrankenTorch), fsci-* (FrankenScipy),
/// fp-* (FrankenPandas).
const CONSTELLATION_PREFIXES: &[&str] = &[
    "asupersync",
    "fsqlite",
    "fnx-",
    "fnp-",
    "ft-",
    "fsci-",
    "fp-",
];

/// Is this dependency name part of the Franken constellation?
fn is_constellation_dep(name: &str) -> bool {
    CONSTELLATION_PREFIXES
        .iter()
        .any(|p| name == p.trim_end_matches('-') || name.starts_with(p))
}

/// The constellation repositories (sibling directories of this workspace):
/// (canonical library name, directory name under `~/projects`).
const CONSTELLATION_REPOS: &[(&str, &str)] = &[
    ("asupersync", "asupersync"),
    ("frankensqlite", "frankensqlite"),
    ("franken_numpy", "franken_numpy"),
    ("frankentorch", "frankentorch"),
    ("frankenscipy", "frankenscipy"),
    ("frankenpandas", "frankenpandas"),
    ("franken_networkx", "franken_networkx"),
];

/// Required CONTRACT.md sections (contract-conformance discipline, plan §13.3).
const CONTRACT_SECTIONS: &[&str] = &[
    "## Purpose and layer",
    "## Public types and semantics",
    "## Invariants",
    "## Error model",
    "## Determinism class",
    "## Cancellation behavior",
    "## Unsafe boundary",
    "## Feature flags",
    "## Conformance tests",
    "## No-claim boundaries",
];

#[derive(Debug)]
struct Manifest {
    name: String,
    layer: Layer,
    /// `[dependencies]` keys only — dev-dependencies are exempt from the layer
    /// rule and the constellation policy (dev-only oracles are permitted when
    /// isolated; they are still REPORTED so the exemption stays visible).
    runtime_deps: Vec<String>,
    dev_deps: Vec<String>,
    dir: PathBuf,
}

#[derive(Debug)]
struct Violation {
    check: &'static str,
    crate_name: String,
    detail: String,
}

#[derive(Debug)]
struct PolicyNote {
    check: &'static str,
    crate_name: String,
    verdict: &'static str,
    detail: String,
}

fn parse_manifest(path: &Path, text: &str) -> Result<Manifest, String> {
    fn value_delimiter_delta(value: &str) -> i32 {
        let mut quote = None;
        let mut escaped = false;
        let mut delta = 0i32;
        for byte in value.bytes() {
            if let Some(delimiter) = quote {
                if escaped {
                    escaped = false;
                } else if byte == b'\\' && delimiter == b'"' {
                    escaped = true;
                } else if byte == delimiter {
                    quote = None;
                }
                continue;
            }
            match byte {
                b'\'' | b'"' => quote = Some(byte),
                b'[' | b'{' => delta += 1,
                b']' | b'}' => delta -= 1,
                b'#' => break,
                _ => {}
            }
        }
        delta
    }

    let mut section = String::new();
    let mut name = None;
    let mut layer = None;
    let mut runtime_deps = Vec::new();
    let mut dev_deps = Vec::new();
    let mut multiline_value_depth = 0i32;

    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if multiline_value_depth > 0 {
            multiline_value_depth += value_delimiter_delta(line);
            if multiline_value_depth < 0 {
                return Err(format!(
                    "{}:{}: multiline value closes more delimiters than it opens",
                    path.display(),
                    lineno + 1
                ));
            }
            continue;
        }
        if line.starts_with('[') {
            section = line.to_string();
            continue;
        }
        let Some(eq) = line.find('=') else {
            return Err(format!(
                "{}:{}: unparseable line {line:?}",
                path.display(),
                lineno + 1
            ));
        };
        let key = line[..eq].trim().trim_matches('"').to_string();
        let value = line[eq + 1..].trim();
        multiline_value_depth = value_delimiter_delta(value);
        if multiline_value_depth < 0 {
            return Err(format!(
                "{}:{}: value closes more delimiters than it opens",
                path.display(),
                lineno + 1
            ));
        }
        match section.as_str() {
            "[package]" if key == "name" => name = Some(value.trim_matches('"').to_string()),
            "[package.metadata.frankensim]" if key == "layer" => {
                let v = value.trim_matches('"');
                layer = Some(
                    Layer::parse(v)
                        .ok_or_else(|| format!("{}: unknown layer {v:?}", path.display()))?,
                );
            }
            "[dependencies]" => runtime_deps.push(key),
            "[dev-dependencies]" => dev_deps.push(key),
            _ => {}
        }
    }
    if multiline_value_depth != 0 {
        return Err(format!(
            "{}: unterminated multiline manifest value",
            path.display()
        ));
    }

    Ok(Manifest {
        name: name.ok_or_else(|| format!("{}: missing package.name", path.display()))?,
        layer: layer.ok_or_else(|| {
            format!(
                "{}: missing [package.metadata.frankensim] layer — every workspace crate must \
                 declare its layer",
                path.display()
            )
        })?,
        runtime_deps,
        dev_deps,
        dir: path.parent().unwrap_or(Path::new(".")).to_path_buf(),
    })
}

fn load_workspace(root: &Path) -> Result<Vec<Manifest>, String> {
    let crates_dir = root.join("crates");
    let mut out = Vec::new();
    let entries = std::fs::read_dir(&crates_dir)
        .map_err(|e| format!("cannot read {}: {e}", crates_dir.display()))?;
    for entry in entries {
        let entry = entry.map_err(|e| e.to_string())?;
        let manifest_path = entry.path().join("Cargo.toml");
        if manifest_path.is_file() {
            let text = std::fs::read_to_string(&manifest_path).map_err(|e| e.to_string())?;
            out.push(parse_manifest(&manifest_path, &text)?);
        }
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    if out.is_empty() {
        return Err(format!(
            "no crate manifests found under {}",
            crates_dir.display()
        ));
    }
    Ok(out)
}

fn check_layers(manifests: &[Manifest]) -> Vec<Violation> {
    let layer_of: BTreeMap<&str, Layer> = manifests
        .iter()
        .map(|m| (m.name.as_str(), m.layer))
        .collect();
    let mut violations = Vec::new();
    for m in manifests {
        for dep in &m.runtime_deps {
            if let Some(&dep_layer) = layer_of.get(dep.as_str())
                && !m.layer.may_depend_on(dep_layer)
            {
                violations.push(Violation {
                    check: "layers",
                    crate_name: m.name.clone(),
                    detail: format!(
                        "{} ({}) must not depend on {} ({}); lower layers never know \
                         about higher ones, and L4/L5 are siblings (plan §4)",
                        m.name,
                        m.layer.name(),
                        dep,
                        dep_layer.name()
                    ),
                });
            }
        }
    }
    violations
}

fn check_deps(manifests: &[Manifest]) -> Vec<Violation> {
    let workspace: BTreeMap<&str, Layer> = manifests
        .iter()
        .map(|m| (m.name.as_str(), m.layer))
        .collect();
    let mut violations = Vec::new();
    for m in manifests {
        if m.layer == Layer::Tool {
            continue; // xtask itself is outside the shipped graph.
        }
        for dep in &m.runtime_deps {
            let is_workspace = workspace.contains_key(dep.as_str());
            let is_constellation = is_constellation_dep(dep);
            if !is_workspace && !is_constellation {
                violations.push(Violation {
                    check: "dependency-policy",
                    crate_name: m.name.clone(),
                    detail: format!(
                        "{} depends on {dep:?}, which is neither a workspace fs-* crate nor a \
                         Franken-constellation library (Decalogue P1: std + constellation only)",
                        m.name
                    ),
                });
            }
        }
    }
    violations
}

fn check_contracts(manifests: &[Manifest]) -> Vec<Violation> {
    let mut violations = Vec::new();
    for m in manifests {
        if m.layer == Layer::Tool {
            continue;
        }
        let contract = m.dir.join("CONTRACT.md");
        match std::fs::read_to_string(&contract) {
            Err(_) => violations.push(Violation {
                check: "contracts",
                crate_name: m.name.clone(),
                detail: format!("{} is missing CONTRACT.md", m.name),
            }),
            Ok(text) => {
                for section in CONTRACT_SECTIONS {
                    if !text.contains(section) {
                        violations.push(Violation {
                            check: "contracts",
                            crate_name: m.name.clone(),
                            detail: format!("{} CONTRACT.md missing section {section:?}", m.name),
                        });
                    }
                }
            }
        }
    }
    violations
}

/// Dev-dependencies are exempt from the constellation policy (dev-only oracles
/// are permitted when isolated) but the exemption must stay VISIBLE: emit a
/// JSON note per external dev-dependency so review sees every use of the escape
/// hatch (docs/CONVENTIONS.md "Dependency policy").
fn dev_dep_notes(manifests: &[Manifest]) -> Vec<(String, String)> {
    let workspace: Vec<&str> = manifests.iter().map(|m| m.name.as_str()).collect();
    let mut notes = Vec::new();
    for m in manifests {
        for dep in &m.dev_deps {
            if !workspace.contains(&dep.as_str()) && !is_constellation_dep(dep) {
                notes.push((m.name.clone(), dep.clone()));
            }
        }
    }
    notes
}

// ---------------------------------------------------------------------------
// Constellation lockfile: pins the exact library states (version + git head)
// this workspace builds against. The lock HASH covers only the portable
// identity (name, version, head) — never local paths — and becomes part of
// every ledger op's Five Explicits once fs-ledger lands (plan §12).
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct ConstellationEntry {
    lib: String,
    dir: PathBuf,
    version: String,
    git_head: String,
    remote: String,
}

/// FNV-1a 64 (mirrors fs-obs::fnv1a64; duplicated here because TOOL crates
/// keep zero workspace deps so policy tooling never blocks on library builds).
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn first_version_in(manifest_text: &str) -> Option<String> {
    for line in manifest_text.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("version") {
            let rest = rest.trim_start();
            if let Some(v) = rest.strip_prefix('=') {
                return Some(v.trim().trim_matches('"').to_string());
            }
        }
    }
    None
}

fn constellation_entries(workspace_root: &Path) -> Result<Vec<ConstellationEntry>, String> {
    let projects = workspace_root
        .parent()
        .ok_or_else(|| "workspace root has no parent".to_string())?;
    let mut out = Vec::new();
    for (lib, dirname) in CONSTELLATION_REPOS {
        let dir = projects.join(dirname);
        let manifest = dir.join("Cargo.toml");
        let text = std::fs::read_to_string(&manifest)
            .map_err(|e| format!("constellation repo {lib} missing at {}: {e}", dir.display()))?;
        let version = first_version_in(&text).unwrap_or_else(|| "unversioned-workspace".into());
        let head = git_out(&dir, &["rev-parse", "HEAD"]).unwrap_or_else(|_| "no-git".to_string());
        let remote = git_out(&dir, &["remote", "get-url", "origin"])
            .unwrap_or_else(|_| "no-remote".to_string());
        out.push(ConstellationEntry {
            lib: (*lib).to_string(),
            dir,
            version,
            git_head: head,
            remote,
        });
    }
    out.sort_by(|a, b| a.lib.cmp(&b.lib));
    Ok(out)
}

/// Canonical portable content the lock hash covers (no local paths).
fn lock_identity(entries: &[ConstellationEntry]) -> String {
    let mut s = String::new();
    for e in entries {
        let _ = writeln!(s, "{}={}@{}", e.lib, e.version, e.git_head);
    }
    s
}

fn render_lock_rows(rows: &[LockRow], lock_hash: &str) -> String {
    let mut s = String::new();
    s.push_str("{\n  \"schema\": \"");
    s.push_str(CONSTELLATION_LOCK_SCHEMA);
    s.push_str("\",\n  \"identity_domain\": \"");
    s.push_str(CONSTELLATION_LOCK_IDENTITY_DOMAIN);
    let _ = write!(
        s,
        "\",\n  \"identity_version\": {CONSTELLATION_LOCK_IDENTITY_VERSION}"
    );
    s.push_str(",\n  \"lock_hash\": \"");
    s.push_str(lock_hash);
    s.push_str("\",\n  \"note\": \"");
    s.push_str(&json_escape(CONSTELLATION_LOCK_NOTE));
    s.push_str("\",\n");
    s.push_str("  \"libraries\": [\n");
    for (index, row) in rows.iter().enumerate() {
        let comma = if index + 1 == rows.len() { "" } else { "," };
        let _ = writeln!(
            s,
            "    {{\"lib\": \"{}\", \"version\": \"{}\", \"git_head\": \"{}\", \"remote\": \"{}\", \"path\": \"{}\"}}{comma}",
            json_escape(&row.lib),
            json_escape(&row.version),
            json_escape(&row.git_head),
            json_escape(&row.remote),
            json_escape(&row.path),
        );
    }
    s.push_str("  ]\n}\n");
    s
}

fn render_lockfile(entries: &[ConstellationEntry], hash: u64) -> String {
    let rows = entries
        .iter()
        .map(|entry| LockRow {
            lib: entry.lib.clone(),
            version: entry.version.clone(),
            git_head: entry.git_head.clone(),
            remote: entry.remote.clone(),
            path: entry.dir.display().to_string(),
        })
        .collect::<Vec<_>>();
    render_lock_rows(&rows, &format!("{hash:016x}"))
}

fn write_constellation_lock(
    path: &Path,
    entries: &[ConstellationEntry],
    hash: u64,
) -> std::io::Result<()> {
    let identity_domain = CONSTELLATION_LOCK_WRITER_IDENTITY_DOMAIN;
    let identity_version = CONSTELLATION_LOCK_WRITER_IDENTITY_VERSION;
    let document = render_lockfile(entries, hash);
    debug_assert!(!document.contains(&format!("\"identity_domain\": \"{identity_domain}\"")));
    let file_name = path.file_name().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "constellation lock path has no file name",
        )
    })?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let (temporary, mut staging) = (0..128)
        .find_map(|_| {
            let suffix = NEXT_CONSTELLATION_LOCK_TEMP_SUFFIX.fetch_add(1, Ordering::Relaxed);
            let mut temporary_name = file_name.to_os_string();
            temporary_name.push(format!(".tmp.{}.{suffix}", std::process::id()));
            let temporary = parent.join(temporary_name);
            match std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
            {
                Ok(file) => Some(Ok((temporary, file))),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => None,
                Err(error) => Some(Err(error)),
            }
        })
        .transpose()?
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::AlreadyExists,
                format!(
                    "constellation lock writer {identity_domain}@{identity_version} could not reserve a staging file beside {}",
                    path.display()
                ),
            )
        })?;
    staging.write_all(document.as_bytes()).map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!(
                "constellation lock writer {identity_domain}@{identity_version} could not stage {} in retained file {}: {error}",
                path.display(),
                temporary.display()
            ),
        )
    })?;
    staging.sync_all().map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!(
                "constellation lock writer {identity_domain}@{identity_version} could not make retained staging file {} durable for {}: {error}",
                temporary.display(),
                path.display()
            ),
        )
    })?;
    drop(staging);
    std::fs::rename(&temporary, path).map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!(
                "constellation lock writer {identity_domain}@{identity_version} could not atomically replace {} from retained staging file {}: {error}",
                path.display(),
                temporary.display()
            ),
        )
    })?;
    if let Ok(directory) = std::fs::File::open(parent) {
        let _ = directory.sync_all();
    }
    Ok(())
}

/// `lock-constellation` writes the lockfile; `check-constellation` verifies
/// the recorded hash still matches the live repos (CI drift gate).
fn cmd_constellation(root: &Path, check: bool) -> ExitCode {
    let entries = match constellation_entries(root) {
        Ok(e) => e,
        Err(msg) => {
            eprintln!("error: {msg}");
            return ExitCode::FAILURE;
        }
    };
    let hash = fnv1a64(lock_identity(&entries).as_bytes());
    let lock_path = root.join("constellation.lock");
    if check {
        let existing = match read_constellation_lock(&lock_path) {
            Ok(existing) => existing,
            Err(error) => {
                eprintln!(
                    "error: {error}; run `cargo run -p xtask -- lock-constellation` if the lock is missing"
                );
                return ExitCode::FAILURE;
            }
        };
        let (recorded, rows) = match parse_lock_rows(&existing) {
            Ok(parsed) => parsed,
            Err(detail) => {
                println!(
                    "{{\"check\":\"constellation-lock\",\"verdict\":\"invalid\",\"detail\":\"{}\"}}",
                    json_escape(&detail)
                );
                eprintln!("constellation lock INVALID: {detail}");
                return ExitCode::FAILURE;
            }
        };
        let declared_identity = match lock_rows_identity(&rows) {
            Ok(identity) => identity,
            Err(detail) => {
                println!(
                    "{{\"check\":\"constellation-lock\",\"verdict\":\"invalid\",\"detail\":\"{}\"}}",
                    json_escape(&detail)
                );
                eprintln!("constellation lock INVALID: {detail}");
                return ExitCode::FAILURE;
            }
        };
        let declared_hash = fnv1a64(declared_identity.as_bytes());
        if recorded != format!("{declared_hash:016x}") {
            println!(
                "{{\"check\":\"constellation-lock\",\"verdict\":\"invalid\",\"recorded\":\"{}\",\"declared\":\"{declared_hash:016x}\"}}",
                json_escape(&recorded)
            );
            eprintln!(
                "constellation lock SELF-DRIFT: recorded {recorded}, declared rows hash to \
                 {declared_hash:016x}"
            );
            return ExitCode::FAILURE;
        }
        if recorded != format!("{hash:016x}") {
            println!(
                "{{\"check\":\"constellation-lock\",\"verdict\":\"drift\",\"recorded\":\"{recorded}\",\"live\":\"{hash:016x}\"}}"
            );
            eprintln!(
                "constellation DRIFT: recorded {recorded}, live {hash:016x}; a constellation \
                 repo moved — re-lock deliberately with lock-constellation and note why"
            );
            return ExitCode::FAILURE;
        }
        if let Err(detail) = verify_constellation_rows(root, &rows) {
            println!(
                "{{\"check\":\"constellation-lock\",\"verdict\":\"refused\",\"hash\":\"{hash:016x}\",\"detail\":\"{}\"}}",
                json_escape(&detail)
            );
            eprintln!("constellation lock REFUSED: {detail}");
            return ExitCode::FAILURE;
        }
        println!(
            "{{\"check\":\"constellation-lock\",\"verdict\":\"ok\",\"hash\":\"{hash:016x}\",\"clean\":true}}"
        );
        eprintln!("constellation lock OK and all sibling trees clean ({hash:016x})");
        ExitCode::SUCCESS
    } else {
        if let Err(detail) = verify_live_entries_clean(&entries) {
            eprintln!("refusing to lock a dirty or unreadable constellation: {detail}");
            return ExitCode::FAILURE;
        }
        match write_constellation_lock(&lock_path, &entries, hash) {
            Ok(()) => {
                println!(
                    "{{\"check\":\"constellation-lock\",\"verdict\":\"written\",\"hash\":\"{hash:016x}\"}}"
                );
                eprintln!("wrote {} (hash {hash:016x})", lock_path.display());
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("error writing lockfile: {e}");
                ExitCode::FAILURE
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Unsafe-capsule enforcement (patch Rev P): every module containing `unsafe`
// must be registered in unsafe-capsules.json, stay under 300 lines, and ship
// a SAFETY.md (docs/SAFETY_TEMPLATE.md). Decalogue P1 made mechanical.
// ---------------------------------------------------------------------------

/// Strip `//` line comments (incl. doc comments) so prose ABOUT unsafe does
/// not trip the scanner. String literals are NOT stripped: a string containing
/// the token forces a visible rename or registration — loud beats silent.
fn strip_line_comments(line: &str) -> &str {
    line.find("//").map_or(line, |i| &line[..i])
}

fn contains_unsafe_token(text: &str) -> bool {
    for raw in text.lines() {
        let line = strip_line_comments(raw);
        for (i, _) in line.match_indices("unsafe") {
            let before_ok = i == 0
                || !line.as_bytes()[i - 1].is_ascii_alphanumeric()
                    && line.as_bytes()[i - 1] != b'_';
            let after = line.as_bytes().get(i + 6);
            let after_ok = after.is_none_or(|&b| !b.is_ascii_alphanumeric() && b != b'_');
            if before_ok && after_ok {
                return true;
            }
        }
    }
    false
}

/// Registered capsule file paths (repo-relative), parsed from the registry's
/// one-capsule-per-line convention: {"crate": "...", "module": "...", ...}.
fn registry_modules(registry_text: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in registry_text.lines() {
        let t = line.trim();
        if t.starts_with("{\"crate\"")
            && let Some(pos) = t.find("\"module\":")
        {
            let rest = &t[pos + 9..];
            if let Some(start) = rest.find('"') {
                let rest = &rest[start + 1..];
                if let Some(end) = rest.find('"') {
                    out.push(rest[..end].to_string());
                }
            }
        }
    }
    out
}

fn check_unsafe(root: &Path) -> Vec<Violation> {
    let registry_text =
        std::fs::read_to_string(root.join("unsafe-capsules.json")).unwrap_or_default();
    let mut violations = Vec::new();
    if registry_text.is_empty() {
        violations.push(Violation {
            check: "unsafe-capsules",
            crate_name: "<repo>".to_string(),
            detail: "unsafe-capsules.json missing at workspace root".to_string(),
        });
        return violations;
    }
    let registered = registry_modules(&registry_text);
    let crates_dir = root.join("crates");
    let Ok(entries) = std::fs::read_dir(&crates_dir) else {
        return violations;
    };
    let mut stack: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path().join("src"))
        .filter(|p| p.is_dir())
        .collect();
    while let Some(dir) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in rd.flatten() {
            let p = entry.path();
            if p.is_dir() {
                stack.push(p);
                continue;
            }
            if p.extension().is_none_or(|e| e != "rs") {
                continue;
            }
            let Ok(text) = std::fs::read_to_string(&p) else {
                continue;
            };
            let has_unsafe = contains_unsafe_token(&text);
            let has_allow =
                text.contains("#[allow(unsafe_code)]") || text.contains("#![allow(unsafe_code)]");
            if !(has_unsafe || has_allow) {
                continue;
            }
            let rel = p.strip_prefix(root).unwrap_or(&p).display().to_string();
            let is_registered = registered.iter().any(|m| rel.ends_with(m) || m == &rel);
            if !is_registered {
                violations.push(Violation {
                    check: "unsafe-capsules",
                    crate_name: rel.clone(),
                    detail: format!(
                        "{rel} contains `unsafe`/allow(unsafe_code) but is not a registered \
                         capsule; register it in unsafe-capsules.json with a SAFETY.md \
                         (docs/SAFETY_TEMPLATE.md) or remove the unsafe"
                    ),
                });
                continue;
            }
            let loc = text.lines().count();
            if loc >= 300 {
                violations.push(Violation {
                    check: "unsafe-capsules",
                    crate_name: rel.clone(),
                    detail: format!(
                        "capsule {rel} is {loc} lines; capsules must stay under 300 (split the \
                         safe parts out from behind the facade)"
                    ),
                });
            }
            let safety = p.parent().map(|d| d.join("SAFETY.md"));
            if safety.as_deref().is_none_or(|s| !s.is_file()) {
                violations.push(Violation {
                    check: "unsafe-capsules",
                    crate_name: rel.clone(),
                    detail: format!("capsule {rel} has no SAFETY.md beside it"),
                });
            }
        }
    }
    violations
}

/// `f64::powi`/`f32::powi` take exactly ONE argument, so a `.powi(` call
/// whose argument list contains a top-level comma is some other method
/// (e.g. the fs-opt expression builders) and is skipped. Returns the
/// single argument text when the call could be the float intrinsic; on
/// an unbalanced (multi-line) call it conservatively returns the rest of
/// the line so hidden float powi cannot slip through.
fn powi_single_arg(after_open: &str) -> Option<&str> {
    let mut depth = 0usize;
    for (i, c) in after_open.char_indices() {
        match c {
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => {
                if depth == 0 {
                    return if c == ')' {
                        Some(&after_open[..i])
                    } else {
                        None
                    };
                }
                depth -= 1;
            }
            ',' if depth == 0 => return None,
            _ => {}
        }
    }
    Some(after_open)
}

/// No optimization-level-dependent integer powers in deterministic paths
/// (bead 4xnt): `f64::powi`'s rounding differs between debug and release
/// from exponent 4 upward (llvm.powi has no pinned operation order), so
/// any `.powi(arg)` in crate sources, tests, examples, or benches must
/// either use a literal exponent
/// in -3..=3 (where all lowerings agree), be migrated to
/// `fs_math::det::powi` (or fs-qty's `powi_pinned`), or carry a
/// `// det-ok: <reason>` annotation on the same or preceding line.
/// Explicit primitive UFCS calls are always flagged because their receiver
/// and exponent need separate parsing; migrate or annotate them locally.
fn check_powi(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();
    let crates_dir = root.join("crates");
    let Ok(entries) = std::fs::read_dir(&crates_dir) else {
        return violations;
    };
    let mut stack: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_dir())
        .collect();
    while let Some(dir) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in rd.flatten() {
            let p = entry.path();
            if p.is_dir() {
                if p.file_name().is_none_or(|name| name != "target") {
                    stack.push(p);
                }
                continue;
            }
            if p.extension().is_none_or(|e| e != "rs") {
                continue;
            }
            let rel = p.strip_prefix(root).unwrap_or(&p).display().to_string();
            let Ok(text) = std::fs::read_to_string(&p) else {
                continue;
            };
            let mut prev_raw = "";
            for (idx, raw) in text.lines().enumerate() {
                let code = strip_line_comments(raw);
                let annotated = raw.contains("det-ok:") || prev_raw.contains("det-ok:");
                let mut flagged = false;

                let mut rest = code;
                while let Some(pos) = rest.find(".powi") {
                    let after_name = &rest[pos + ".powi".len()..];
                    if let Some(after) = after_name.trim_start().strip_prefix('(')
                        && let Some(arg) = powi_single_arg(after)
                    {
                        let literal_ok = arg
                            .trim()
                            .parse::<i64>()
                            .is_ok_and(|v| (-3..=3).contains(&v));
                        if !literal_ok && !annotated && !flagged {
                            violations.push(Violation {
                                check: "powi-determinism",
                                crate_name: rel.clone(),
                                detail: format!(
                                    "{rel}:{}: `.powi({})` — f64::powi rounding is \
                                     optimization-level-dependent (1-ULP debug/release \
                                     divergence from exponent 4; bead 4xnt); use \
                                     fs_math::det::powi, or a literal exponent in -3..=3, \
                                     or annotate `// det-ok: <reason>` on this or the \
                                     preceding line",
                                    idx + 1,
                                    arg.trim()
                                ),
                            });
                            flagged = true;
                        }
                    }
                    rest = after_name;
                }

                if !annotated
                    && !flagged
                    && ["f32::powi", "f64::powi", "<f32>::powi", "<f64>::powi"]
                        .iter()
                        .any(|needle| code.contains(needle))
                {
                    violations.push(Violation {
                        check: "powi-determinism",
                        crate_name: rel.clone(),
                        detail: format!(
                            "{rel}:{}: explicit primitive `powi` call bypasses the pinned deterministic primitive; use fs_math::det::powi or annotate `// det-ok: <reason>` on this or the preceding line",
                            idx + 1
                        ),
                    });
                }
                prev_raw = raw;
            }
        }
    }
    violations
}

/// Crates whose CONTRACT claims CROSS-ISA bitwise determinism and have been
/// audited onto the det:: routing doctrine (bead frankensim-lyms). Tiered
/// policy from that bead: a crate claiming cross-ISA/bitwise determinism
/// must route transcendentals through `fs_math::det` (platform libm differs
/// by ≥1 ULP across ISAs/libm versions); crates claiming only same-ISA
/// determinism are exempt and are NOT listed here. Grow this list as crates
/// are audited — adding a name turns the doctrine on for it.
const LIBM_DOCTRINE_CRATES: &[&str] = &["fs-geocon", "fs-toleralloc", "fs-uq"];

/// Method-call transcendentals that platform libm implements without a
/// correct-rounding guarantee. `sqrt` is deliberately absent (IEEE-754
/// requires correct rounding, so it is bit-stable everywhere), as are
/// exact ops (`abs`, `rem_euclid`, `floor`, ...).
const LIBM_METHODS: &[&str] = &[
    "cbrt", "ln", "ln_1p", "log2", "log10", "exp", "exp2", "exp_m1", "expm1", "sin", "cos",
    "sin_cos", "tan", "asin", "acos", "atan", "atan2", "sinh", "cosh", "tanh", "asinh", "acosh",
    "atanh", "hypot", "powf",
];

/// check-libm (bead frankensim-lyms): inside [`LIBM_DOCTRINE_CRATES`],
/// every raw libm transcendental method call must be migrated to its
/// `fs_math::det` equivalent or carry a `// det-ok: <reason>` annotation
/// on the same or preceding line (the escape hatch for dev-only oracle
/// comparisons). Same annotation discipline as `check-powi`.
fn check_libm(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();
    for crate_name in LIBM_DOCTRINE_CRATES {
        let crate_dir = root.join("crates").join(crate_name);
        let mut stack = vec![crate_dir];
        while let Some(dir) = stack.pop() {
            let Ok(rd) = std::fs::read_dir(&dir) else {
                continue;
            };
            for entry in rd.flatten() {
                let p = entry.path();
                if p.is_dir() {
                    if p.file_name().is_none_or(|name| name != "target") {
                        stack.push(p);
                    }
                    continue;
                }
                if p.extension().is_none_or(|e| e != "rs") {
                    continue;
                }
                let rel = p.strip_prefix(root).unwrap_or(&p).display().to_string();
                let Ok(text) = std::fs::read_to_string(&p) else {
                    continue;
                };
                let mut prev_raw = "";
                for (idx, raw) in text.lines().enumerate() {
                    let code = strip_line_comments(raw);
                    let annotated = raw.contains("det-ok:") || prev_raw.contains("det-ok:");
                    prev_raw = raw;
                    if annotated {
                        continue;
                    }
                    // `det::sin(x)` is a path call and never matches the
                    // `.sin(` method needle, so routed code cannot false-flag.
                    for method in LIBM_METHODS {
                        let needle = format!(".{method}(");
                        if code.contains(&needle) {
                            violations.push(Violation {
                                check: "libm-determinism",
                                crate_name: rel.clone(),
                                detail: format!(
                                    "{rel}:{}: raw `.{method}(` — platform libm is not \
                                     correctly rounded and differs across ISAs; this crate's \
                                     CONTRACT claims cross-ISA bitwise determinism (bead \
                                     lyms), so use fs_math::det::{method} (or det::pow for \
                                     cbrt/powf), or annotate `// det-ok: <reason>`",
                                    idx + 1
                                ),
                            });
                        }
                    }
                }
            }
        }
    }
    violations
}

/// The two color-admission AUTHORITY crates (bead 6pf9): fs-evidence owns
/// the Color/AdmittedColor types themselves and fs-ledger's ColorGraph is
/// the replay-audited admission oracle. Positive literals inside them are
/// the mechanism, not a bypass.
const COLOR_AUTHORITY_CRATES: &[&str] = &["fs-evidence", "fs-ledger"];

/// Grandfathered direct positive-Color constructors (bead 6pf9 stage S3
/// worklist, surveyed 2026-07-16 at d5a2da0). This list is a RATCHET: it
/// only shrinks. Remove a crate here once its owner group migrates its
/// positive-evidence surfaces to `fs_evidence::AdmittedColor` (or renames
/// remaining construction sites under explicit declared/unverified
/// discipline); never add a crate — new code routes through an admission
/// authority (`ColorGraph::admission_receipt`,
/// `VerifiedPackage::claim_admission_receipt`,
/// `LoopReport::admitted_headline`) or carries a
/// `// declared-color-ok: <reason>` annotation.
const COLOR_DECLARED_CRATES: &[&str] = &[];

/// check-color-admission (bead 6pf9, stage S4 lint, check-powi/check-libm
/// precedent): a direct positive `Color::Verified`/`Color::Validated`
/// literal outside the admission authorities is capability fabrication —
/// structural validation cannot tell an admitted certificate from a
/// fabricated literal. Inside `crates/*/src`, every such literal must live
/// in an authority crate, a grandfathered [`COLOR_DECLARED_CRATES`] entry,
/// test code (a `tests/` path or below a `#[cfg(test)]` marker), or carry a
/// `// declared-color-ok: <reason>` annotation on the same or preceding
/// line.
fn check_color_admission(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();
    let crates_dir = root.join("crates");
    let Ok(crates) = std::fs::read_dir(&crates_dir) else {
        return violations;
    };
    for krate in crates.flatten() {
        let crate_path = krate.path();
        let Some(crate_name) = crate_path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if COLOR_AUTHORITY_CRATES.contains(&crate_name)
            || COLOR_DECLARED_CRATES.contains(&crate_name)
        {
            continue;
        }
        let mut stack = vec![crate_path.join("src")];
        while let Some(dir) = stack.pop() {
            let Ok(rd) = std::fs::read_dir(&dir) else {
                continue;
            };
            for entry in rd.flatten() {
                let p = entry.path();
                if p.is_dir() {
                    stack.push(p);
                    continue;
                }
                if p.extension().is_none_or(|e| e != "rs") {
                    continue;
                }
                let rel = p.strip_prefix(root).unwrap_or(&p).display().to_string();
                let Ok(text) = std::fs::read_to_string(&p) else {
                    continue;
                };
                let mut prev_raw = "";
                let mut in_test_code = false;
                for (idx, raw) in text.lines().enumerate() {
                    // Test fixtures legitimately declare colors; skip the
                    // conventional trailing test module once its cfg marker
                    // appears.
                    if raw.contains("#[cfg(test)]") {
                        in_test_code = true;
                    }
                    if in_test_code {
                        continue;
                    }
                    let code = strip_line_comments(raw);
                    let annotated = raw.contains("declared-color-ok:")
                        || prev_raw.contains("declared-color-ok:");
                    prev_raw = raw;
                    if annotated {
                        continue;
                    }
                    // Pattern READS are not fabrication: destructuring a
                    // positive color in a match arm, `matches!`, or `if let`
                    // only branches on rank. Skip the three lexical read
                    // shapes; a construction hidden on the same line as an
                    // arm body is a tolerated false negative — this check is
                    // a shrink-only ratchet, not a proof.
                    if code.contains("=>") || code.contains("matches!") || code.contains("if let ")
                    {
                        continue;
                    }
                    for variant in ["Color::Verified", "Color::Validated"] {
                        if code.contains(variant) {
                            violations.push(Violation {
                                check: "color-admission",
                                crate_name: rel.clone(),
                                detail: format!(
                                    "{rel}:{}: direct positive `{variant}` literal — a \
                                     fabricated literal is indistinguishable from admitted \
                                     evidence (bead 6pf9); mint through an admission authority \
                                     (ColorGraph::admission_receipt, \
                                     VerifiedPackage::claim_admission_receipt) and consume \
                                     fs_evidence::AdmittedColor, or annotate \
                                     `// declared-color-ok: <reason>`",
                                    idx + 1
                                ),
                            });
                        }
                    }
                }
            }
        }
    }
    violations
}

/// check-obs-events (bead huq.16): every committed `*.events.jsonl` file must
/// contain only schema-valid fs-obs event lines. The schema authority is
/// `fs_obs::validate_line` itself — this check deliberately depends on the
/// crate instead of re-implementing the wire format, so xtask can never
/// drift into a second dialect. No fixtures exist yet: the check arms the
/// test-log convention and enforces from the first committed fixture
/// (citable-scanner precedent).
fn check_obs_events(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();
    let mut stack = vec![root.join("crates"), root.join("data")];
    while let Some(dir) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in rd.flatten() {
            let p = entry.path();
            if p.is_dir() {
                if p.file_name().is_none_or(|name| name != "target") {
                    stack.push(p);
                }
                continue;
            }
            let Some(name) = p.file_name().and_then(|n| n.to_str()) else {
                continue;
            };
            if !name.ends_with(".events.jsonl") {
                continue;
            }
            let rel = p.strip_prefix(root).unwrap_or(&p).display().to_string();
            let Ok(text) = std::fs::read_to_string(&p) else {
                violations.push(Violation {
                    check: "obs-events",
                    crate_name: rel.clone(),
                    detail: format!("{rel}: unreadable event fixture"),
                });
                continue;
            };
            for (idx, line) in text.lines().enumerate() {
                if line.is_empty() {
                    continue;
                }
                if let Err(error) = fs_obs::validate_line(line) {
                    violations.push(Violation {
                        check: "obs-events",
                        crate_name: rel.clone(),
                        detail: format!(
                            "{rel}:{}: invalid fs-obs event line at byte {}: {}",
                            idx + 1,
                            error.at,
                            error.message
                        ),
                    });
                }
            }
        }
    }
    violations
}

// ---------------------------------------------------------------------------
// Sole-branch terminology policy (bead sj31i.54). The prohibited token is
// assembled from fragments so this guard and its mutation fixtures do not
// exempt themselves from the same tracked-source scan they enforce.
// ---------------------------------------------------------------------------

const TERMINOLOGY_CHECK: &str = "branch-terminology";
const TERMINOLOGY_ALLOWLIST_VERSION: u32 = 2;

fn legacy_branch_word() -> String {
    ["mas", "ter"].concat()
}

/// Compound technical terms in which the legacy branch word is NOT a
/// branch reference and cannot be renamed: SQLite's system-catalog
/// table (external API), the Maurer-Cartan/BV "<word> equation" of
/// mathematics (spaced or hyphenated), and DMA/bus device-arbitration
/// roles (hardware architecture). Allowlist v2; each compound is
/// stripped before the branch-term scan so any residual bare
/// occurrence on the same line still trips.
fn terminology_allowed_compounds(word: &str) -> [String; 5] {
    [
        format!("sqlite_{word}"),
        format!("{word} equation"),
        format!("{word}-equation"),
        format!("dma {word}"),
        format!("bus {word}"),
    ]
}

fn terminology_policy_line_is_allowlisted(path: &str, line: &str, word: &str) -> bool {
    if path != "AGENTS.md" {
        return false;
    }
    [
        format!("## Git Branch: ONLY Use `main`, NEVER `{word}`"),
        format!("- Never reference `{word}` in code or docs. If you see it, treat it as a bug."),
        format!("- If the remote also needs a legacy `{word}` ref, synchronize it from `main`"),
    ]
    .iter()
    .any(|allowed| line == allowed)
}

fn scan_terminology_sources(sources: &BTreeMap<String, String>) -> Vec<Violation> {
    let word = legacy_branch_word();
    let needle = word.to_ascii_lowercase();
    let compounds = terminology_allowed_compounds(&needle);
    let mut violations = Vec::new();
    for (path, source) in sources {
        for (line_index, line) in source.lines().enumerate() {
            let mut scrubbed = line.to_ascii_lowercase();
            for compound in &compounds {
                scrubbed = scrubbed.replace(compound.as_str(), "");
            }
            if scrubbed.contains(&needle)
                && !terminology_policy_line_is_allowlisted(path, line, &word)
            {
                violations.push(Violation {
                    check: TERMINOLOGY_CHECK,
                    crate_name: path.clone(),
                    detail: format!(
                        "{path}:{}: prohibited legacy branch term {word:?}; use main for the branch, root/study for seeds, or representative/constrained for domain reductions (allowlist v{TERMINOLOGY_ALLOWLIST_VERSION})",
                        line_index + 1
                    ),
                });
            }
        }
    }
    violations
}

fn terminology_scans_path(path: &str) -> bool {
    if path.starts_with(".beads/") {
        return false;
    }
    Path::new(path)
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| {
            matches!(
                extension,
                "c" | "cc"
                    | "cpp"
                    | "h"
                    | "hpp"
                    | "js"
                    | "json"
                    | "jsonl"
                    | "lean"
                    | "md"
                    | "py"
                    | "rs"
                    | "sh"
                    | "toml"
                    | "ts"
                    | "txt"
                    | "yaml"
                    | "yml"
            )
        })
}

fn tracked_terminology_sources(root: &Path) -> Result<BTreeMap<String, String>, String> {
    let output = std::process::Command::new("git")
        .args(["ls-files", "-z"])
        .current_dir(root)
        .output()
        .map_err(|error| format!("cannot inventory tracked files: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "git ls-files failed with status {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    if !output.stdout.is_empty() && !output.stdout.ends_with(&[0]) {
        return Err("git ls-files inventory is not NUL-terminated".to_string());
    }

    let mut sources = BTreeMap::new();
    for raw_path in output.stdout.split(|byte| *byte == 0) {
        if raw_path.is_empty() {
            continue;
        }
        let path = std::str::from_utf8(raw_path)
            .map_err(|_| "git ls-files emitted a non-UTF-8 path".to_string())?;
        if !terminology_scans_path(path) {
            continue;
        }
        let bytes = std::fs::read(root.join(path))
            .map_err(|error| format!("cannot read tracked source {path:?}: {error}"))?;
        let source = String::from_utf8(bytes)
            .map_err(|_| format!("tracked source {path:?} is not valid UTF-8"))?;
        sources.insert(path.to_string(), source);
    }
    Ok(sources)
}

fn check_terminology(root: &Path) -> Vec<Violation> {
    match tracked_terminology_sources(root) {
        Ok(sources) => scan_terminology_sources(&sources),
        Err(detail) => vec![Violation {
            check: TERMINOLOGY_CHECK,
            crate_name: "<repo>".to_string(),
            detail,
        }],
    }
}

// ---------------------------------------------------------------------------
// Citable-producer inventory (bead pd16): a positive citation bit is a
// security boundary. Every Rust string/debug field that can emit that bit is
// discovered mechanically and must remain at one audited, authority-gated
// sink. A fixed-false report-only serializer is classified separately: it
// cannot make a positive claim, but it still needs an explicit allowlist row
// so a later edit cannot silently turn it into a positive-capable sink.
// ---------------------------------------------------------------------------

const CITABLE_PRODUCER_CHECK: &str = "citable-producer-inventory";
const CITATION_FIELD_NAME: &str = concat!("citation_", "eligible");
const CITATION_FIELD_SEAM: usize = 9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CitationFieldMode {
    PositiveCapable,
    ReportOnlyFalse,
}

impl CitationFieldMode {
    fn name(self) -> &'static str {
        match self {
            Self::PositiveCapable => "positive-capable",
            Self::ReportOnlyFalse => "report-only-fixed-false",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CitationFieldOccurrence {
    path: String,
    owner: String,
    anchor_text: String,
    mode: CitationFieldMode,
    line: usize,
}

#[derive(Debug, Clone, Copy)]
struct CitableProducerSpec {
    path: &'static str,
    owner: &'static str,
    /// A stable schema/metric marker in the same string literal. The Debug
    /// field has no schema marker, so path + owner + exact multiplicity is its
    /// stable identity.
    anchor: &'static str,
    mode: CitationFieldMode,
    /// File-local source-shape guards that tie a positive-capable serializer
    /// to the typed authority path. These are deliberately semantic markers,
    /// not line numbers.
    authority_markers: &'static [&'static str],
    /// Markers proving the same entry point has an explicit fixed-false
    /// report-only branch (or a separate candidate namespace).
    report_only_markers: &'static [&'static str],
}

const CITABLE_PRODUCER_ALLOWLIST: &[CitableProducerSpec] = &[
    CitableProducerSpec {
        path: "crates/fs-roofline/src/production.rs",
        owner: "fmt",
        anchor: "",
        mode: CitationFieldMode::PositiveCapable,
        authority_markers: &[
            concat!("&self.citation_", "eligible()"),
            "self.admission_error().is_none()",
            "CitationAuthority::Receipt(binding)",
        ],
        report_only_markers: &[
            "pub struct ReportOnlyProductionRun",
            "crate::EvidenceNamespace::Custom",
        ],
    },
    CitableProducerSpec {
        path: "crates/fs-roofline/src/bin/roofline.rs",
        owner: "evidence_admission_json",
        anchor: "fs-roofline-evidence-admission-v2",
        mode: CitationFieldMode::PositiveCapable,
        authority_markers: &[
            concat!(
                "let (citation_",
                "eligible, admission_error) = run.evidence_admission();"
            ),
            "Self::Attested(run) => {",
            "(refusal.is_none(), refusal)",
        ],
        report_only_markers: &["Self::ReportOnly(run) => (false, run.admission_error())"],
    },
    CitableProducerSpec {
        path: "crates/fs-roofline/src/bin/roofline.rs",
        owner: "main",
        anchor: "fs-roofline-recorded-evidence-v2",
        mode: CitationFieldMode::PositiveCapable,
        authority_markers: &[
            concat!(
                "let (citation_",
                "eligible, admission_error) = run.evidence_admission();"
            ),
            "let recorded = match run.record(&ledger)",
            "let dependency_authority = load_dependency_authority(&args);",
            ".revalidate(&ledger, &current)",
            "dependency_authority_policy_receipt",
            concat!("let citable = citation_", "eligible && revalidated_fresh;"),
        ],
        report_only_markers: &["Self::ReportOnly(run) => (false, run.admission_error())"],
    },
    CitableProducerSpec {
        path: "crates/fs-feec/tests/perf_lane.rs",
        owner: "sum_factorized_attainment",
        anchor: "feec-gate",
        mode: CitationFieldMode::PositiveCapable,
        authority_markers: &[
            "configured_citable_ledger(",
            "classify_gate_admission(&snapshot, configuration_refusal)",
            "GateAdmission::Citable => (true, None)",
            "Ok(path) if !path.is_empty() && path != \":memory:\"",
            "record_external_perf_gate_at_path(",
            concat!("\\\"recorded\\\":{citation_", "eligible}"),
        ],
        report_only_markers: &["GateAdmission::ReportOnly(reason) => (false, Some(reason))"],
    },
    CitableProducerSpec {
        path: "crates/fs-fft/tests/perf_lane.rs",
        owner: "fft_attainment",
        anchor: "fft-gate",
        mode: CitationFieldMode::PositiveCapable,
        authority_markers: &[
            "configured_citable_ledger(",
            "classify_gate_admission(&snapshot, configuration_refusal)",
            "GateAdmission::Citable => (true, None)",
            "Ok(path) if !path.is_empty() && path != \":memory:\"",
            "record_external_perf_gate_at_path(",
            concat!("\\\"recorded\\\":{citation_", "eligible}"),
        ],
        report_only_markers: &["GateAdmission::ReportOnly(reason) => (false, Some(reason))"],
    },
];

#[derive(Debug, Clone, Copy)]
struct RustStringLiteral<'a> {
    start: usize,
    end: usize,
    content: &'a str,
}

fn char_literal_end(bytes: &[u8], start: usize) -> Option<usize> {
    let tail = core::str::from_utf8(bytes.get(start + 1..)?).ok()?;
    let first = tail.chars().next()?;
    if first != '\\' {
        let close = start + 1 + first.len_utf8();
        // A non-escaped Rust char contains exactly one Unicode scalar. If the
        // next byte is not the closing quote, this apostrophe begins a lifetime
        // or label (including Unicode XID spellings), not an opaque literal.
        return (bytes.get(close) == Some(&b'\'')).then_some(close + 1);
    }

    // Escaped chars have variable-width spellings (`\\xNN`, `\\u{...}`, and
    // escaped quotes). Retain the existing fail-closed newline boundary while
    // finding the first unescaped closing quote.
    let mut escaped = false;
    for (offset, byte) in bytes[start + 1..].iter().copied().enumerate() {
        let index = start + 1 + offset;
        if escaped {
            escaped = false;
        } else if byte == b'\\' {
            escaped = true;
        } else if byte == b'\'' {
            return Some(index + 1);
        } else if byte == b'\n' {
            return None;
        }
    }
    None
}

fn raw_string_open(bytes: &[u8], start: usize) -> Option<(usize, usize)> {
    if start > 0
        && bytes
            .get(start - 1)
            .is_some_and(|byte| byte.is_ascii_alphanumeric() || *byte == b'_')
    {
        return None;
    }
    let mut cursor = start;
    if matches!(bytes.get(cursor), Some(b'b' | b'c')) {
        cursor += 1;
    }
    if bytes.get(cursor) != Some(&b'r') {
        return None;
    }
    cursor += 1;
    let hash_start = cursor;
    while bytes.get(cursor) == Some(&b'#') {
        cursor += 1;
    }
    (bytes.get(cursor) == Some(&b'"')).then_some((cursor, cursor - hash_start))
}

/// Extract Rust string literals while ignoring line/block comments and char
/// literals. Both cooked and raw strings are supported because a raw JSON
/// literal is just as capable of creating a new evidence producer.
fn rust_string_literals(source: &str) -> Result<Vec<RustStringLiteral<'_>>, String> {
    let bytes = source.as_bytes();
    let mut literals = Vec::new();
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        if bytes.get(cursor..cursor + 2) == Some(b"//") {
            cursor = source[cursor..]
                .find('\n')
                .map_or(bytes.len(), |offset| cursor + offset + 1);
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"/*") {
            let mut depth = 1usize;
            cursor += 2;
            while cursor < bytes.len() && depth > 0 {
                if bytes.get(cursor..cursor + 2) == Some(b"/*") {
                    depth += 1;
                    cursor += 2;
                } else if bytes.get(cursor..cursor + 2) == Some(b"*/") {
                    depth -= 1;
                    cursor += 2;
                } else {
                    cursor += 1;
                }
            }
            if depth != 0 {
                return Err("unterminated block comment while scanning Rust strings".to_string());
            }
            continue;
        }
        if bytes[cursor] == b'\''
            && let Some(end) = char_literal_end(bytes, cursor)
        {
            cursor = end;
            continue;
        }
        if let Some((quote, hashes)) = raw_string_open(bytes, cursor) {
            let content_start = quote + 1;
            let mut close = content_start;
            let end = loop {
                let Some(relative) = source[close..].find('"') else {
                    return Err("unterminated raw string while scanning producers".to_string());
                };
                let quote_end = close + relative;
                let suffix_end = quote_end + 1 + hashes;
                if suffix_end <= bytes.len()
                    && bytes[quote_end + 1..suffix_end]
                        .iter()
                        .all(|byte| *byte == b'#')
                {
                    break suffix_end;
                }
                close = quote_end + 1;
            };
            literals.push(RustStringLiteral {
                start: cursor,
                end,
                content: &source[content_start..end - 1 - hashes],
            });
            cursor = end;
            continue;
        }
        if bytes[cursor] == b'"' {
            let content_start = cursor + 1;
            let mut escaped = false;
            let mut end = None;
            for (offset, byte) in bytes[content_start..].iter().copied().enumerate() {
                let index = content_start + offset;
                if escaped {
                    escaped = false;
                } else if byte == b'\\' {
                    escaped = true;
                } else if byte == b'"' {
                    end = Some(index + 1);
                    break;
                }
            }
            let Some(end) = end else {
                return Err("unterminated cooked string while scanning producers".to_string());
            };
            literals.push(RustStringLiteral {
                start: cursor,
                end,
                content: &source[content_start..end - 1],
            });
            cursor = end;
            continue;
        }
        cursor += 1;
    }
    Ok(literals)
}

fn rust_offset_is_code(source: &str, offset: usize, literals: &[RustStringLiteral<'_>]) -> bool {
    let bytes = source.as_bytes();
    let mut cursor = 0usize;
    let mut literal_index = 0usize;
    while cursor <= offset && cursor < bytes.len() {
        while literals
            .get(literal_index)
            .is_some_and(|literal| literal.end <= cursor)
        {
            literal_index += 1;
        }
        if let Some(literal) = literals.get(literal_index)
            && literal.start == cursor
        {
            if offset < literal.end {
                return false;
            }
            cursor = literal.end;
            literal_index += 1;
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"//") {
            let end = source[cursor..]
                .find('\n')
                .map_or(bytes.len(), |relative| cursor + relative + 1);
            if offset < end {
                return false;
            }
            cursor = end;
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"/*") {
            let mut depth = 1usize;
            cursor += 2;
            while cursor < bytes.len() && depth > 0 {
                if bytes.get(cursor..cursor + 2) == Some(b"/*") {
                    depth += 1;
                    cursor += 2;
                } else if bytes.get(cursor..cursor + 2) == Some(b"*/") {
                    depth -= 1;
                    cursor += 2;
                } else {
                    cursor += 1;
                }
            }
            if offset < cursor || depth != 0 {
                return false;
            }
            continue;
        }
        if bytes[cursor] == b'\''
            && let Some(end) = char_literal_end(bytes, cursor)
        {
            if offset < end {
                return false;
            }
            cursor = end;
            continue;
        }
        cursor += 1;
    }
    cursor > offset
}

fn rust_declaration_start(
    source: &str,
    literals: &[RustStringLiteral<'_>],
    declaration: &str,
) -> Option<usize> {
    let mut line_start = 0usize;
    for line in source.split_inclusive('\n') {
        let leading = line
            .bytes()
            .take_while(|byte| matches!(*byte, b' ' | b'\t'))
            .count();
        let candidate = line_start + leading;
        let trimmed = &line[leading..];
        if let Some(after) = trimmed.strip_prefix(declaration)
            && after.as_bytes().first().is_some_and(|byte| {
                matches!(
                    *byte,
                    b':' | b'<' | b';' | b'{' | b' ' | b'\t' | b'\r' | b'\n'
                )
            })
            && rust_offset_is_code(source, candidate, literals)
        {
            return Some(candidate);
        }
        line_start += line.len();
    }
    None
}

fn declared_function_name(raw: &str) -> Option<String> {
    let mut line = raw.trim_start();
    if line.starts_with("//") || line.starts_with("/*") || line.starts_with('*') {
        return None;
    }
    if line.starts_with("pub(")
        && let Some(close) = line.find(')')
    {
        line = line[close + 1..].trim_start();
    } else if let Some(rest) = line.strip_prefix("pub ") {
        line = rest.trim_start();
    }
    for qualifier in ["async ", "const ", "unsafe "] {
        if let Some(rest) = line.strip_prefix(qualifier) {
            line = rest.trim_start();
        }
    }
    let rest = line.strip_prefix("fn ")?;
    let name: String = rest
        .chars()
        .take_while(|ch| ch.is_ascii_alphanumeric() || *ch == '_')
        .collect();
    (!name.is_empty()).then_some(name)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum RustOwnerScope {
    Module,
    Function(usize),
}

fn enclosing_function(source: &str, offset: usize) -> (String, RustOwnerScope) {
    let mut line_end = offset;
    loop {
        let line_start = source[..line_end].rfind('\n').map_or(0, |index| index + 1);
        if let Some(name) = declared_function_name(&source[line_start..line_end])
            && let Some(relative_open) = source[line_start..offset].find('{')
        {
            let open = line_start + relative_open;
            if matches!(rust_scope_contains(source, open, offset), Ok(true)) {
                return (name, RustOwnerScope::Function(line_start));
            }
        }
        if line_start == 0 {
            break;
        }
        line_end = line_start - 1;
    }
    ("<module>".to_string(), RustOwnerScope::Module)
}

fn matcher_only_literal(source: &str, offset: usize) -> bool {
    let line_start = source[..offset].rfind('\n').map_or(0, |index| index + 1);
    let prefix = source[line_start..offset].trim_end();
    [".contains(", ".starts_with(", ".ends_with("]
        .iter()
        .any(|matcher| prefix.ends_with(matcher))
}

fn parser_key_literal(source: &str, literal: RustStringLiteral<'_>) -> bool {
    literal.content == CITATION_FIELD_NAME && source[literal.end..].trim_start().starts_with("=>")
}

fn literal_is_citation_field(literal: RustStringLiteral<'_>) -> bool {
    if literal.content == CITATION_FIELD_NAME {
        return true;
    }
    let compact: String = literal
        .content
        .chars()
        .filter(|ch| !ch.is_whitespace())
        .collect();
    let escaped_key = format!("\\\"{CITATION_FIELD_NAME}\\\":");
    let raw_key = format!("\"{CITATION_FIELD_NAME}\":");
    compact.contains(&escaped_key) || compact.contains(&raw_key)
}

/// The scanner's own protected-name lexicon and allowlist deliberately split
/// the field name so the guard does not inventory its policy data as an output
/// sink. Only this declaration region is exempt; executable xtask functions
/// remain subject to the same producer scan as every other Rust source.
#[derive(Debug, Clone, Copy)]
struct CitableGuardMetadataRanges {
    field: (usize, usize),
    allowlist: (usize, usize),
}

fn citable_guard_metadata_ranges(
    path: &str,
    source: &str,
    literals: &[RustStringLiteral<'_>],
) -> Option<CitableGuardMetadataRanges> {
    if path != "xtask/src/main.rs" {
        return None;
    }
    let field_start = rust_declaration_start(source, literals, "const CITATION_FIELD_NAME")?;
    let allowlist_start =
        rust_declaration_start(source, literals, "const CITABLE_PRODUCER_ALLOWLIST")?;
    let struct_start = rust_declaration_start(source, literals, "struct RustStringLiteral")?;
    if !(field_start < allowlist_start && allowlist_start < struct_start) {
        return None;
    }
    let field_end = source[field_start..allowlist_start]
        .find(';')
        .map(|relative| field_start + relative + 1)?;
    let allowlist_end = source[allowlist_start..struct_start]
        .find("];")
        .map(|relative| allowlist_start + relative + 2)?;
    Some(CitableGuardMetadataRanges {
        field: (field_start, field_end),
        allowlist: (allowlist_start, allowlist_end),
    })
}

fn citable_guard_metadata_literal(
    source: &str,
    literal: RustStringLiteral<'_>,
    ranges: Option<CitableGuardMetadataRanges>,
) -> bool {
    let Some(ranges) = ranges else {
        return false;
    };
    ((ranges.field.0..ranges.field.1).contains(&literal.start)
        || (ranges.allowlist.0..ranges.allowlist.1).contains(&literal.start))
        && matches!(
            enclosing_function(source, literal.start).1,
            RustOwnerScope::Module
        )
}

/// Conservative per-function state for the field's fixed semantic seam:
/// `citation_` + `eligible`. Seeing both halves in distinct non-consumer
/// literals is enough to inventory a positive-capable dynamic construction;
/// declaration order does not matter.
#[derive(Debug, Clone, Copy, Default)]
struct CitationFragmentState<'a> {
    prefix: Option<RustStringLiteral<'a>>,
    suffix: Option<RustStringLiteral<'a>>,
}

fn literal_has_dynamic_citation_fragments(
    content: &str,
    citation_prefix: &str,
    citation_suffix: &str,
) -> bool {
    if !content.contains(citation_prefix) || !content.contains(citation_suffix) {
        return false;
    }
    if content.contains(CITATION_FIELD_NAME) {
        // Remove every complete spelling first. This keeps ordinary diagnostic
        // repetitions and format captures such as `{citation_eligible}` from
        // masquerading as split producers, while an unmatched extra half
        // remains visible and fail-closed.
        let residual = content.replace(CITATION_FIELD_NAME, "");
        residual.contains(citation_prefix) || residual.contains(citation_suffix)
    } else {
        true
    }
}

/// Whether `offset` is still inside the braced scope beginning at `open`.
/// The small lexer mirrors `rust_string_literals`: braces in comments,
/// strings, raw strings, and char literals do not affect the scope depth.
fn rust_scope_contains(source: &str, open: usize, offset: usize) -> Result<bool, String> {
    let bytes = source.as_bytes();
    let mut cursor = open;
    let mut depth = 0usize;
    while cursor < offset {
        if bytes.get(cursor..cursor + 2) == Some(b"//") {
            cursor = source[cursor..]
                .find('\n')
                .map_or(offset, |relative| (cursor + relative + 1).min(offset));
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"/*") {
            let mut comment_depth = 1usize;
            cursor += 2;
            while cursor < offset && comment_depth > 0 {
                if bytes.get(cursor..cursor + 2) == Some(b"/*") {
                    comment_depth += 1;
                    cursor += 2;
                } else if bytes.get(cursor..cursor + 2) == Some(b"*/") {
                    comment_depth -= 1;
                    cursor += 2;
                } else {
                    cursor += 1;
                }
            }
            if comment_depth != 0 {
                return Err("unterminated block comment while finding test scope".to_string());
            }
            continue;
        }
        if bytes[cursor] == b'\''
            && let Some(end) = char_literal_end(bytes, cursor)
        {
            cursor = end.min(offset);
            continue;
        }
        if let Some((quote, hashes)) = raw_string_open(bytes, cursor) {
            let mut close = quote + 1;
            let end = loop {
                let Some(relative) = source[close..].find('"') else {
                    return Err("unterminated raw string while finding test scope".to_string());
                };
                let quote_end = close + relative;
                let suffix_end = quote_end + 1 + hashes;
                if suffix_end <= bytes.len()
                    && bytes[quote_end + 1..suffix_end]
                        .iter()
                        .all(|byte| *byte == b'#')
                {
                    break suffix_end;
                }
                close = quote_end + 1;
            };
            cursor = end.min(offset);
            continue;
        }
        if bytes[cursor] == b'"' {
            let mut escaped = false;
            let mut end = None;
            for (relative, byte) in bytes[cursor + 1..].iter().copied().enumerate() {
                let index = cursor + 1 + relative;
                if escaped {
                    escaped = false;
                } else if byte == b'\\' {
                    escaped = true;
                } else if byte == b'"' {
                    end = Some(index + 1);
                    break;
                }
            }
            let Some(end) = end else {
                return Err("unterminated cooked string while finding test scope".to_string());
            };
            cursor = end.min(offset);
            continue;
        }
        match bytes[cursor] {
            b'{' => depth += 1,
            b'}' => {
                let Some(next) = depth.checked_sub(1) else {
                    return Ok(false);
                };
                depth = next;
                if depth == 0 {
                    return Ok(false);
                }
            }
            _ => {}
        }
        cursor += 1;
    }
    Ok(depth > 0)
}

fn cfg_test_only_literal(source: &str, offset: usize) -> Result<bool, String> {
    let mut search_end = offset;
    while let Some(attribute) = source[..search_end].rfind("#[cfg(test)]") {
        let after_attribute = &source[attribute + "#[cfg(test)]".len()..offset];
        let Some(module) = after_attribute.find("mod tests") else {
            search_end = attribute;
            continue;
        };
        let after_module = attribute + "#[cfg(test)]".len() + module + "mod tests".len();
        let Some(relative_open) = source[after_module..offset].find('{') else {
            search_end = attribute;
            continue;
        };
        let open = after_module + relative_open;
        if rust_scope_contains(source, open, offset)? {
            return Ok(true);
        }
        search_end = attribute;
    }
    Ok(false)
}

fn citation_field_mode(source: &str, literal: RustStringLiteral<'_>) -> CitationFieldMode {
    let compact: String = literal
        .content
        .chars()
        .filter(|ch| !ch.is_whitespace())
        .collect();
    let escaped_false = format!("\\\"{CITATION_FIELD_NAME}\\\":false");
    let raw_false = format!("\"{CITATION_FIELD_NAME}\":false");
    let field_call_is_false = literal.content == CITATION_FIELD_NAME
        && source[literal.end..].trim_start().starts_with(", &false");
    if compact.contains(&escaped_false) || compact.contains(&raw_false) || field_call_is_false {
        CitationFieldMode::ReportOnlyFalse
    } else {
        CitationFieldMode::PositiveCapable
    }
}

fn scan_citable_producer_source(
    path: &str,
    source: &str,
) -> Result<Vec<CitationFieldOccurrence>, String> {
    let mut occurrences = Vec::new();
    let mut fragment_states = BTreeMap::<RustOwnerScope, CitationFragmentState<'_>>::new();
    let (citation_prefix, citation_suffix) = CITATION_FIELD_NAME.split_at(CITATION_FIELD_SEAM);
    let literals = rust_string_literals(source)?;
    let guard_metadata = citable_guard_metadata_ranges(path, source, &literals);
    for literal in literals {
        if matcher_only_literal(source, literal.start)
            || parser_key_literal(source, literal)
            || cfg_test_only_literal(source, literal.start)?
            || citable_guard_metadata_literal(source, literal, guard_metadata)
        {
            continue;
        }
        let (owner, owner_scope) = enclosing_function(source, literal.start);
        let dynamic_fragments = literal_has_dynamic_citation_fragments(
            literal.content,
            citation_prefix,
            citation_suffix,
        );
        if literal_is_citation_field(literal) {
            occurrences.push(CitationFieldOccurrence {
                path: path.to_string(),
                owner: owner.clone(),
                anchor_text: literal.content.to_string(),
                mode: citation_field_mode(source, literal),
                line: source[..literal.start]
                    .bytes()
                    .filter(|byte| *byte == b'\n')
                    .count()
                    + 1,
            });
            if dynamic_fragments {
                occurrences.push(CitationFieldOccurrence {
                    path: path.to_string(),
                    owner,
                    anchor_text: CITATION_FIELD_NAME.to_string(),
                    mode: CitationFieldMode::PositiveCapable,
                    line: source[..literal.start]
                        .bytes()
                        .filter(|byte| *byte == b'\n')
                        .count()
                        + 1,
                });
            }
            continue;
        }

        let state = fragment_states.entry(owner_scope).or_default();
        // A diagnostic or parser message may spell the whole vocabulary in
        // one CONTIGUOUS literal without emitting it as a field. Noncontiguous
        // halves inside one literal can be split and rejoined dynamically, so
        // they are themselves a positive-capable occurrence.
        if dynamic_fragments {
            occurrences.push(CitationFieldOccurrence {
                path: path.to_string(),
                owner,
                anchor_text: CITATION_FIELD_NAME.to_string(),
                mode: CitationFieldMode::PositiveCapable,
                line: source[..literal.start]
                    .bytes()
                    .filter(|byte| *byte == b'\n')
                    .count()
                    + 1,
            });
            continue;
        }
        let has_prefix = literal.content.contains(citation_prefix);
        let has_suffix = literal.content.contains(citation_suffix);
        if has_prefix && has_suffix {
            continue;
        }
        if state.prefix.is_none() && has_prefix {
            state.prefix = Some(literal);
        }
        if state.suffix.is_none() && has_suffix {
            state.suffix = Some(literal);
        }
        if let (Some(prefix), Some(suffix)) = (state.prefix, state.suffix) {
            state.prefix = None;
            state.suffix = None;
            let start = if prefix.start <= suffix.start {
                prefix
            } else {
                suffix
            };
            occurrences.push(CitationFieldOccurrence {
                path: path.to_string(),
                owner,
                anchor_text: CITATION_FIELD_NAME.to_string(),
                // A fragmented construction cannot earn the direct
                // literal's fixed-false classification: dynamic use is
                // conservatively positive-capable.
                mode: CitationFieldMode::PositiveCapable,
                line: source[..start.start]
                    .bytes()
                    .filter(|byte| *byte == b'\n')
                    .count()
                    + 1,
            });
        }
    }
    Ok(occurrences)
}

#[allow(clippy::too_many_lines)] // exact five-sink inventory + authority/report-only proofs
fn audit_citable_producer_sources(sources: &BTreeMap<String, String>) -> Vec<Violation> {
    let mut violations = Vec::new();
    let mut occurrences = Vec::new();
    for (path, source) in sources {
        match scan_citable_producer_source(path, source) {
            Ok(found) => occurrences.extend(found),
            Err(error) => violations.push(Violation {
                check: CITABLE_PRODUCER_CHECK,
                crate_name: path.clone(),
                detail: format!("cannot exhaustively scan {path}: {error}"),
            }),
        }
    }

    let mut matched = vec![0usize; occurrences.len()];
    for spec in CITABLE_PRODUCER_ALLOWLIST {
        let matching: Vec<usize> = occurrences
            .iter()
            .enumerate()
            .filter(|(_, occurrence)| {
                occurrence.path == spec.path
                    && occurrence.owner == spec.owner
                    && occurrence.mode == spec.mode
                    && occurrence.anchor_text.contains(spec.anchor)
            })
            .map(|(index, _)| index)
            .collect();
        match matching.as_slice() {
            [_] => {}
            [] => violations.push(Violation {
                check: CITABLE_PRODUCER_CHECK,
                crate_name: spec.path.to_string(),
                detail: format!(
                    "audited {} sink {}#{} (anchor {:?}) is missing, moved, or changed mode; update the authority proof and allowlist together",
                    spec.mode.name(),
                    spec.path,
                    spec.owner,
                    spec.anchor,
                ),
            }),
            _ => violations.push(Violation {
                check: CITABLE_PRODUCER_CHECK,
                crate_name: spec.path.to_string(),
                detail: format!(
                    "audited sink {}#{} (anchor {:?}) occurs {} times; every producer needs a distinct reviewed allowlist row",
                    spec.path,
                    spec.owner,
                    spec.anchor,
                    matching.len(),
                ),
            }),
        }
        for index in matching {
            matched[index] += 1;
        }

        let Some(source) = sources.get(spec.path) else {
            continue;
        };
        for marker in spec.authority_markers {
            if !source.contains(marker) {
                violations.push(Violation {
                    check: CITABLE_PRODUCER_CHECK,
                    crate_name: spec.path.to_string(),
                    detail: format!(
                        "audited sink {}#{} lost authority marker {:?}; positive citation output must remain downstream of admitted authority and durable recording",
                        spec.path, spec.owner, marker,
                    ),
                });
            }
        }
        for marker in spec.report_only_markers {
            if !source.contains(marker) {
                violations.push(Violation {
                    check: CITABLE_PRODUCER_CHECK,
                    crate_name: spec.path.to_string(),
                    detail: format!(
                        "audited sink {}#{} lost report-only marker {:?}; candidate/refused paths must serialize a fixed false citation state",
                        spec.path, spec.owner, marker,
                    ),
                });
            }
        }
    }

    for (occurrence, match_count) in occurrences.iter().zip(matched) {
        if match_count == 1 {
            continue;
        }
        let classification = if occurrence.mode == CitationFieldMode::ReportOnlyFalse {
            "this fixed-false serializer cannot make a positive claim, but must still be explicitly audited"
        } else {
            "a positive-capable producer must be downstream of admitted authority and durable recording"
        };
        violations.push(Violation {
            check: CITABLE_PRODUCER_CHECK,
            crate_name: occurrence.path.clone(),
            detail: format!(
                "unexpected {} citation field at {}:{} in fn {}: {classification}; add a precise path/function/schema allowlist row only after reviewing both positive and report-only branches",
                occurrence.mode.name(),
                occurrence.path,
                occurrence.line,
                occurrence.owner,
            ),
        });
    }
    violations
}

/// Scan every repository-owned Rust source root directly. This deliberately
/// does not ask Git for a file list: a newly-created source file enters the
/// inventory immediately, before it can be staged or committed.
fn workspace_rust_sources(root: &Path) -> Result<BTreeMap<String, String>, String> {
    let roots = ["crates", "tools", "xtask"].map(|directory| root.join(directory));
    let mut stack = Vec::from(roots);
    let mut paths = Vec::new();
    while let Some(directory) = stack.pop() {
        let metadata = std::fs::symlink_metadata(&directory)
            .map_err(|error| format!("cannot inspect {}: {error}", directory.display()))?;
        if metadata.file_type().is_symlink() {
            return Err(format!(
                "refusing symlinked Rust-source inventory root {}",
                directory.display()
            ));
        }
        let entries = std::fs::read_dir(&directory)
            .map_err(|error| format!("cannot read {}: {error}", directory.display()))?;
        for entry in entries {
            let entry = entry
                .map_err(|error| format!("cannot enumerate {}: {error}", directory.display()))?;
            let path = entry.path();
            let file_type = entry
                .file_type()
                .map_err(|error| format!("cannot inspect {}: {error}", path.display()))?;
            if file_type.is_symlink() {
                return Err(format!(
                    "refusing symlink under Rust-source inventory roots: {}",
                    path.display()
                ));
            }
            if file_type.is_dir() {
                // The only excluded build root is the repository-level
                // `target/`, which is outside these three explicitly-owned
                // source roots. A directory merely named `target` below a
                // source root is ordinary source and must be inventoried.
                stack.push(path);
            } else if file_type.is_file() && path.extension().is_some_and(|ext| ext == "rs") {
                paths.push(path);
                if paths.len() > CASUAL_MAX_SOURCE_FILES {
                    return Err(format!(
                        "Rust-source inventory exceeds the {CASUAL_MAX_SOURCE_FILES}-file audit cap"
                    ));
                }
            }
        }
    }
    paths.sort();

    let mut sources = BTreeMap::new();
    let mut total_bytes = 0usize;
    for path in paths {
        let relative = casual_inventory_relative(root, &path)?;
        let remaining = CASUAL_MAX_SOURCE_BYTES
            .checked_sub(total_bytes)
            .ok_or_else(|| "Rust-source inventory byte count overflowed".to_string())?;
        let source = casual_read_bounded_utf8(&path, remaining, "Rust source")
            .map_err(|error| format!("cannot read Rust source {relative}: {error}"))?;
        total_bytes = total_bytes
            .checked_add(source.len())
            .ok_or_else(|| "Rust-source inventory byte count overflowed".to_string())?;
        if total_bytes > CASUAL_MAX_SOURCE_BYTES {
            return Err(format!(
                "Rust-source inventory exceeds the {CASUAL_MAX_SOURCE_BYTES}-byte audit cap"
            ));
        }
        if sources.insert(relative.clone(), source).is_some() {
            return Err(format!(
                "multiple filesystem paths collapse to Rust-source inventory key {relative}"
            ));
        }
    }
    Ok(sources)
}

fn casual_inventory_relative(root: &Path, path: &Path) -> Result<String, String> {
    let relative = path
        .strip_prefix(root)
        .map_err(|error| format!("{} escaped workspace root: {error}", path.display()))?;
    let mut components = Vec::new();
    for component in relative.components() {
        match component {
            Component::CurDir => {}
            Component::Normal(component) => {
                let component = component.to_str().ok_or_else(|| {
                    format!(
                        "refusing non-UTF-8 inventory path component in {}",
                        path.display()
                    )
                })?;
                if component.contains('\\') {
                    return Err(format!(
                        "refusing backslash-bearing inventory path component {component:?} in {}; portable graph identity would be ambiguous",
                        path.display()
                    ));
                }
                components.push(component);
            }
            Component::ParentDir | Component::RootDir | Component::Prefix(_) => {
                return Err(format!(
                    "{} is not a normalized workspace-relative inventory path",
                    path.display()
                ));
            }
        }
    }
    if components.is_empty() {
        return Err(format!(
            "{} does not name an entry below the workspace root",
            path.display()
        ));
    }
    Ok(components.join("/"))
}

fn casual_read_bounded_utf8(path: &Path, limit: usize, kind: &str) -> Result<String, String> {
    let file = std::fs::File::open(path)
        .map_err(|error| format!("cannot open {}: {error}", path.display()))?;
    let metadata = file
        .metadata()
        .map_err(|error| format!("cannot inspect {}: {error}", path.display()))?;
    if !metadata.is_file() {
        return Err(format!("{} is not a regular file", path.display()));
    }
    let advertised = usize::try_from(metadata.len()).map_err(|_| {
        format!(
            "{kind} {} length does not fit the host address space",
            path.display()
        )
    })?;
    if advertised > limit {
        return Err(format!(
            "{kind} {} exceeds the remaining {limit}-byte audit budget",
            path.display()
        ));
    }

    // Metadata is only a preflight: the file may grow between stat and read.
    // Limit the reader itself to one byte beyond the remaining budget so that
    // a concurrent replacement/growth race still refuses before unbounded
    // allocation or I/O.
    let read_limit = u64::try_from(limit).unwrap_or(u64::MAX).saturating_add(1);
    let mut source = String::with_capacity(advertised);
    file.take(read_limit)
        .read_to_string(&mut source)
        .map_err(|error| format!("cannot read {} as UTF-8: {error}", path.display()))?;
    if source.len() > limit {
        return Err(format!(
            "{kind} {} grew beyond the remaining {limit}-byte audit budget while being read",
            path.display()
        ));
    }
    Ok(source)
}

fn check_citable_producers(root: &Path) -> Vec<Violation> {
    match workspace_rust_sources(root) {
        Ok(sources) => audit_citable_producer_sources(&sources),
        Err(detail) => vec![Violation {
            check: CITABLE_PRODUCER_CHECK,
            crate_name: "<repo>".to_string(),
            detail,
        }],
    }
}

// ---------------------------------------------------------------------------
// No-casual-print policy (bead i94v.7.3.3). Process output is data: core
// libraries return typed values/events while CLI binaries own presentation.
// ---------------------------------------------------------------------------

const CASUAL_PRINT_CHECK: &str = "casual-print";
const CASUAL_MAX_SOURCE_FILES: usize = 8_192;
const CASUAL_MAX_SOURCE_BYTES: usize = 256 * 1024 * 1024;
const CASUAL_MAX_MANIFEST_BYTES: usize = 1024 * 1024;
const CASUAL_MAX_TOKENS_PER_SOURCE: usize = 2_000_000;
const CASUAL_MAX_TOTAL_TOKENS: usize = 32_000_000;
const CASUAL_MAX_REACHABLE_SOURCES: usize = 16_384;
const CASUAL_MAX_MODULE_EDGES: usize = 65_536;
const CASUAL_MAX_MODULE_DEPTH: usize = 1_024;
const CASUAL_MAX_PORTABLE_PATH_BYTES: usize = 16 * 1_024;
const CASUAL_MAX_DIAGNOSTICS: usize = 256;

#[derive(Debug, Clone, Copy)]
struct CasualPrintAllowance {
    path: &'static str,
    owner: &'static str,
    invocation_anchors: &'static [&'static str],
    reason: &'static str,
}

/// Pre-policy structured emitters. Each allowance names one unique function
/// body and its exact normalized macro invocations. A duplicate same-named
/// owner, added invocation, changed macro, or changed argument token stream
/// invalidates the whole allowance instead of inheriting it.
const CASUAL_PRINT_ALLOWLIST: &[CasualPrintAllowance] = &[
    CasualPrintAllowance {
        path: "crates/fs-casebook/src/lib.rs",
        owner: "run",
        invocation_anchors: &[
            r#"::std::println!("{}",record.json_line())"#,
            r#"::std::println!("{}",replay_record.json_line())"#,
            r#"::std::println!("{}",disagreement.json_line())"#,
        ],
        reason: "legacy casebook JSONL harness emitter",
    },
    CasualPrintAllowance {
        path: "crates/fs-propcheck/src/lib.rs",
        owner: "check_structured",
        invocation_anchors: &[r#"::std::println!("{failure_row}")"#],
        reason: "legacy failing-counterexample JSONL harness emitter",
    },
    CasualPrintAllowance {
        path: "crates/fs-vskeleton/src/lib.rs",
        owner: "emit",
        invocation_anchors: &[r#"::std::eprintln!("{}",e.to_jsonl())"#],
        reason: "legacy fs-obs JSONL application adapter",
    },
];

#[derive(Debug, Clone, PartialEq, Eq)]
struct CasualPrintOccurrence {
    path: String,
    owner: String,
    owner_start: Option<usize>,
    macro_name: String,
    invocation_anchor: String,
    alias_of: Option<String>,
    in_macro_tokens: bool,
    line: usize,
}

#[derive(Debug, Clone, Default)]
struct CasualPrintScan {
    occurrences: Vec<CasualPrintOccurrence>,
    function_starts: BTreeMap<String, Vec<usize>>,
    authority_hazards: Vec<(usize, String)>,
}

#[derive(Debug, Clone, Copy)]
struct CasualRustToken<'a> {
    start: usize,
    text: &'a str,
}

#[derive(Debug, Clone)]
struct CasualFunctionScope {
    name: String,
    declaration_start: usize,
    open: usize,
    close: usize,
}

fn casual_path_literal(token: &str) -> Result<String, String> {
    let bytes = token.as_bytes();
    if token.starts_with('r')
        && let Some((quote, hashes)) = raw_string_open(bytes, 0)
    {
        let close = bytes
            .len()
            .checked_sub(hashes + 1)
            .ok_or_else(|| "truncated raw #[path] literal".to_string())?;
        if bytes.get(close) != Some(&b'"') || !bytes[close + 1..].iter().all(|byte| *byte == b'#') {
            return Err("malformed raw #[path] literal".to_string());
        }
        if close.saturating_sub(quote + 1) > CASUAL_MAX_PORTABLE_PATH_BYTES {
            return Err(format!(
                "#[path] literal exceeds the {CASUAL_MAX_PORTABLE_PATH_BYTES}-byte portable-path cap"
            ));
        }
        return Ok(token[quote + 1..close].to_string());
    }
    if bytes.first() != Some(&b'"') || bytes.last() != Some(&b'"') {
        return Err("#[path] requires a cooked or raw UTF-8 string literal".to_string());
    }

    let mut output = String::new();
    let mut characters = token[1..token.len() - 1].chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            output.push(character);
            if output.len() > CASUAL_MAX_PORTABLE_PATH_BYTES {
                return Err(format!(
                    "#[path] literal exceeds the {CASUAL_MAX_PORTABLE_PATH_BYTES}-byte portable-path cap"
                ));
            }
            continue;
        }
        let escaped = characters
            .next()
            .ok_or_else(|| "truncated escape in #[path] literal".to_string())?;
        match escaped {
            '\\' => output.push('\\'),
            '"' => output.push('"'),
            '\'' => output.push('\''),
            'n' => output.push('\n'),
            'r' => output.push('\r'),
            't' => output.push('\t'),
            '0' => output.push('\0'),
            'x' => {
                let high = characters
                    .next()
                    .and_then(|digit| digit.to_digit(16))
                    .ok_or_else(|| "invalid \\x escape in #[path] literal".to_string())?;
                let low = characters
                    .next()
                    .and_then(|digit| digit.to_digit(16))
                    .ok_or_else(|| "invalid \\x escape in #[path] literal".to_string())?;
                let value = high * 16 + low;
                if value > 0x7f {
                    return Err("non-ASCII \\x escape in #[path] literal".to_string());
                }
                output.push(char::from_u32(value).expect("ASCII escape is a scalar"));
            }
            'u' => {
                if characters.next() != Some('{') {
                    return Err("invalid \\u escape in #[path] literal".to_string());
                }
                let mut value = 0_u32;
                let mut digits = 0_u8;
                loop {
                    let digit = characters
                        .next()
                        .ok_or_else(|| "unterminated \\u escape in #[path] literal".to_string())?;
                    if digit == '}' {
                        break;
                    }
                    let digit = digit
                        .to_digit(16)
                        .ok_or_else(|| "invalid \\u escape in #[path] literal".to_string())?;
                    value = value
                        .checked_mul(16)
                        .and_then(|value| value.checked_add(digit))
                        .ok_or_else(|| "overflowing \\u escape in #[path] literal".to_string())?;
                    digits = digits.saturating_add(1);
                    if digits > 6 {
                        return Err("overlong \\u escape in #[path] literal".to_string());
                    }
                }
                if digits == 0 {
                    return Err("empty \\u escape in #[path] literal".to_string());
                }
                output.push(
                    char::from_u32(value)
                        .ok_or_else(|| "non-scalar \\u escape in #[path] literal".to_string())?,
                );
            }
            _ => {
                return Err(format!("unsupported escape \\{escaped} in #[path] literal"));
            }
        }
        if output.len() > CASUAL_MAX_PORTABLE_PATH_BYTES {
            return Err(format!(
                "#[path] literal exceeds the {CASUAL_MAX_PORTABLE_PATH_BYTES}-byte portable-path cap"
            ));
        }
    }
    Ok(output)
}

fn casual_normalized_inclusion(source_path: &str, literal: &str) -> Result<String, String> {
    if literal.contains('\\') {
        return Err(format!(
            "#[path] `{literal}` contains a backslash; portable source identity requires `/` separators"
        ));
    }
    if literal.starts_with('/') {
        return Err(format!("#[path] `{literal}` must be workspace-relative"));
    }

    let mut normalized = source_path
        .rsplit_once('/')
        .map_or_else(Vec::new, |(parent, _)| {
            parent.split('/').map(str::to_string).collect()
        });
    for component in literal.split('/') {
        match component {
            "" | "." => {}
            ".." => {
                if normalized.pop().is_none() {
                    return Err(format!("#[path] `{literal}` escapes the workspace root"));
                }
            }
            component => normalized.push(component.to_string()),
        }
    }
    let normalized = normalized.join("/");
    if normalized.len() > CASUAL_MAX_PORTABLE_PATH_BYTES {
        return Err(format!(
            "#[path] `{literal}` resolves beyond the {CASUAL_MAX_PORTABLE_PATH_BYTES}-byte portable-path cap"
        ));
    }
    Ok(normalized)
}

fn casual_identifier_start(character: char) -> bool {
    character == '_'
        || character.is_alphabetic()
        // Unicode's Other_ID_Start set is part of XID_Start but is not
        // uniformly covered by `char::is_alphabetic`. Keep the complete
        // stable set explicit so valid Rust identifiers are not split.
        || matches!(character, '\u{1885}' | '\u{1886}' | '\u{2118}' | '\u{212e}' | '\u{309b}' | '\u{309c}')
}

fn casual_identifier_continue(character: char) -> bool {
    character == '_'
        || character.is_alphanumeric()
        // Rust uses Unicode XID_Continue. `std` does not expose that table, so
        // conservatively retain non-ASCII continuation scalars (combining
        // marks and join controls included) inside one token. ASCII Rust
        // punctuation remains a hard token boundary.
        || (!character.is_ascii() && !character.is_whitespace() && !character.is_control())
}

fn casual_token_is_identifier(token: &str) -> bool {
    let mut characters = token.chars();
    characters.next().is_some_and(casual_identifier_start)
        && characters.all(casual_identifier_continue)
}

fn casual_identifier_at<'source>(
    tokens: &[CasualRustToken<'source>],
    index: usize,
) -> Option<(&'source str, usize)> {
    if tokens.get(index)?.text == "r"
        && tokens.get(index + 1).is_some_and(|token| token.text == "#")
        && tokens
            .get(index + 2)
            .is_some_and(|token| casual_token_is_identifier(token.text))
    {
        Some((tokens[index + 2].text, index + 3))
    } else {
        casual_token_is_identifier(tokens[index].text).then_some((tokens[index].text, index + 1))
    }
}

#[allow(clippy::too_many_lines)] // one fail-closed lexer keeps every trivia/string boundary aligned
fn casual_rust_tokens(source: &str) -> Result<Vec<CasualRustToken<'_>>, String> {
    let bytes = source.as_bytes();
    let mut tokens = Vec::new();
    macro_rules! push_token {
        ($token:expr) => {{
            if tokens.len() == CASUAL_MAX_TOKENS_PER_SOURCE {
                return Err(format!(
                    "source exceeds the {CASUAL_MAX_TOKENS_PER_SOURCE}-token audit cap"
                ));
            }
            tokens.push($token);
        }};
    }
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        if bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"//") {
            cursor = source[cursor..]
                .find('\n')
                .map_or(bytes.len(), |relative| cursor + relative + 1);
            continue;
        }
        if bytes.get(cursor..cursor + 2) == Some(b"/*") {
            let mut depth = 1usize;
            cursor += 2;
            while cursor < bytes.len() && depth > 0 {
                if bytes.get(cursor..cursor + 2) == Some(b"/*") {
                    depth += 1;
                    cursor += 2;
                } else if bytes.get(cursor..cursor + 2) == Some(b"*/") {
                    depth -= 1;
                    cursor += 2;
                } else {
                    cursor += 1;
                }
            }
            if depth != 0 {
                return Err("unterminated block comment while scanning casual output".to_string());
            }
            continue;
        }
        if let Some((quote, hashes)) = raw_string_open(bytes, cursor) {
            let mut close = quote + 1;
            let end = loop {
                let Some(relative) = source[close..].find('"') else {
                    return Err("unterminated raw string while scanning casual output".to_string());
                };
                let quote_end = close + relative;
                let suffix_end = quote_end + 1 + hashes;
                if suffix_end <= bytes.len()
                    && bytes[quote_end + 1..suffix_end]
                        .iter()
                        .all(|byte| *byte == b'#')
                {
                    break suffix_end;
                }
                close = quote_end + 1;
            };
            push_token!(CasualRustToken {
                start: cursor,
                text: &source[cursor..end],
            });
            cursor = end;
            continue;
        }
        if bytes[cursor] == b'"' {
            let mut escaped = false;
            let mut end = None;
            for (relative, byte) in bytes[cursor + 1..].iter().copied().enumerate() {
                let index = cursor + 1 + relative;
                if escaped {
                    escaped = false;
                } else if byte == b'\\' {
                    escaped = true;
                } else if byte == b'"' {
                    end = Some(index + 1);
                    break;
                }
            }
            let Some(end) = end else {
                return Err("unterminated cooked string while scanning casual output".to_string());
            };
            push_token!(CasualRustToken {
                start: cursor,
                text: &source[cursor..end],
            });
            cursor = end;
            continue;
        }
        if bytes[cursor] == b'\''
            && let Some(end) = char_literal_end(bytes, cursor)
        {
            push_token!(CasualRustToken {
                start: cursor,
                text: &source[cursor..end],
            });
            cursor = end;
            continue;
        }
        let first = source[cursor..]
            .chars()
            .next()
            .expect("cursor remains inside source");
        if casual_identifier_start(first) {
            let start = cursor;
            cursor += first.len_utf8();
            while cursor < bytes.len() {
                let character = source[cursor..]
                    .chars()
                    .next()
                    .expect("cursor remains inside source");
                if !casual_identifier_continue(character) {
                    break;
                }
                cursor += character.len_utf8();
            }
            push_token!(CasualRustToken {
                start,
                text: &source[start..cursor],
            });
            continue;
        }
        let width = source[cursor..]
            .chars()
            .next()
            .expect("cursor remains inside source")
            .len_utf8();
        push_token!(CasualRustToken {
            start: cursor,
            text: &source[cursor..cursor + width],
        });
        cursor += width;
    }
    Ok(tokens)
}

fn casual_delimiter_pairs(tokens: &[CasualRustToken<'_>]) -> Result<Vec<Option<usize>>, String> {
    let mut pairs = vec![None; tokens.len()];
    let mut stack = Vec::<(&str, usize)>::new();
    for (index, token) in tokens.iter().enumerate() {
        match token.text {
            "(" | "[" | "{" => stack.push((token.text, index)),
            ")" | "]" | "}" => {
                let Some((open, open_index)) = stack.pop() else {
                    return Err(format!(
                        "unmatched closing delimiter {} at byte {}",
                        token.text, token.start
                    ));
                };
                let expected = match open {
                    "(" => ")",
                    "[" => "]",
                    "{" => "}",
                    _ => unreachable!("only opening delimiters enter the stack"),
                };
                if token.text != expected {
                    return Err(format!(
                        "mismatched delimiter {open} ... {} at byte {}",
                        token.text, token.start
                    ));
                }
                pairs[open_index] = Some(index);
                pairs[index] = Some(open_index);
            }
            _ => {}
        }
    }
    if let Some((open, index)) = stack.pop() {
        return Err(format!(
            "unterminated delimiter {open} at byte {}",
            tokens[index].start
        ));
    }
    Ok(pairs)
}

fn casual_macro_token_spans(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
) -> Vec<(usize, usize)> {
    let mut spans = Vec::new();
    for (open, token) in tokens.iter().enumerate() {
        if !matches!(token.text, "(" | "[" | "{") {
            continue;
        }
        let invocation = open >= 2
            && tokens[open - 1].text == "!"
            && casual_token_is_identifier(tokens[open - 2].text);
        let macro_rules_body = (open >= 3
            && tokens[open - 3].text == "macro_rules"
            && tokens[open - 2].text == "!"
            && casual_token_is_identifier(tokens[open - 1].text))
            || (open >= 5
                && tokens[open - 5].text == "macro_rules"
                && tokens[open - 4].text == "!"
                && tokens[open - 3].text == "r"
                && tokens[open - 2].text == "#"
                && casual_token_is_identifier(tokens[open - 1].text));
        if (invocation || macro_rules_body)
            && let Some(close) = pairs[open]
        {
            spans.push((open, close));
        }
    }
    spans
}

#[derive(Debug, Clone, Copy)]
struct CasualAttribute {
    hash: usize,
    bracket: usize,
    close: usize,
    inner: bool,
}

fn casual_attributes(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
) -> Vec<CasualAttribute> {
    let mut attributes = Vec::new();
    for hash in 0..tokens.len() {
        if tokens[hash].text != "#" {
            continue;
        }
        let inner = tokens.get(hash + 1).is_some_and(|token| token.text == "!");
        let bracket = hash + if inner { 2 } else { 1 };
        if tokens.get(bracket).is_some_and(|token| token.text == "[")
            && let Some(close) = pairs[bracket]
        {
            attributes.push(CasualAttribute {
                hash,
                bracket,
                close,
                inner,
            });
        }
    }
    attributes
}

fn casual_span_mask(token_count: usize, spans: &[(usize, usize)]) -> Vec<bool> {
    let mut deltas = vec![0_i32; token_count.saturating_add(1)];
    for &(start, end) in spans {
        if start >= token_count || end < start {
            continue;
        }
        deltas[start] += 1;
        if end + 1 < deltas.len() {
            deltas[end + 1] -= 1;
        }
    }
    let mut depth = 0_i32;
    deltas
        .into_iter()
        .take(token_count)
        .map(|delta| {
            depth += delta;
            depth > 0
        })
        .collect()
}

fn casual_span_end_map(token_count: usize, spans: &[(usize, usize)]) -> Vec<Option<usize>> {
    let mut ends = vec![None; token_count];
    for &(start, end) in spans {
        if start < token_count {
            ends[start] = Some(ends[start].map_or(end, |existing| existing.max(end)));
        }
    }
    ends
}

fn casual_structural_token_mask(
    token_count: usize,
    macro_spans: &[(usize, usize)],
    attributes: &[CasualAttribute],
) -> (Vec<bool>, Vec<bool>) {
    let macro_mask = casual_span_mask(token_count, macro_spans);
    let attribute_interiors = attributes
        .iter()
        .filter_map(|attribute| {
            (attribute.bracket + 1 < attribute.close)
                .then_some((attribute.bracket + 1, attribute.close - 1))
        })
        .collect::<Vec<_>>();
    let attribute_mask = casual_span_mask(token_count, &attribute_interiors);
    let structural_mask = macro_mask
        .iter()
        .zip(&attribute_mask)
        .map(|(macro_data, attribute_data)| *macro_data || *attribute_data)
        .collect();
    (macro_mask, structural_mask)
}

fn casual_delimiter_containers(tokens: &[CasualRustToken<'_>]) -> Vec<Option<usize>> {
    let mut containers = vec![None; tokens.len()];
    let mut stack = Vec::new();
    for (index, token) in tokens.iter().enumerate() {
        containers[index] = stack.last().copied();
        match token.text {
            "(" | "[" | "{" => stack.push(index),
            ")" | "]" | "}" => {
                stack.pop();
            }
            _ => {}
        }
    }
    containers
}

fn casual_cfg_test_attribute(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
    index: usize,
) -> Option<(usize, bool)> {
    if tokens.get(index)?.text != "#" {
        return None;
    }
    let inner = tokens.get(index + 1).is_some_and(|token| token.text == "!");
    let bracket = index + if inner { 2 } else { 1 };
    if tokens.get(bracket)?.text != "[" {
        return None;
    }
    let close = pairs.get(bracket).copied().flatten()?;
    let (name, cursor) = casual_identifier_at(tokens, bracket + 1)?;
    let (predicate, predicate_end) = casual_identifier_at(tokens, cursor + 1)?;
    (name == "cfg"
        && tokens.get(cursor)?.text == "("
        && predicate == "test"
        && tokens.get(predicate_end)?.text == ")"
        && close == predicate_end + 1)
        .then_some((close, inner))
}

fn casual_cfg_test_spans(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
    structural_mask: &[bool],
    macro_span_ends: &[Option<usize>],
    containers: &[Option<usize>],
) -> Result<(bool, Vec<(usize, usize)>), String> {
    let mut file_test_only = false;
    let mut spans = Vec::new();
    let mut covered_until = 0usize;
    for index in 0..tokens.len() {
        if index < covered_until {
            continue;
        }
        // Macro inputs/bodies and every arbitrary attribute interior are token
        // data, not structural Rust syntax. Only the outer shell of this exact
        // recognized cfg attribute is permitted to confer authority.
        if structural_mask[index] {
            continue;
        }
        let Some((attribute_close, inner)) = casual_cfg_test_attribute(tokens, pairs, index) else {
            continue;
        };
        if inner {
            match containers[index] {
                Some(open) if tokens[open].text == "{" => {
                    let close = pairs[open].expect("validated containing brace has a close");
                    spans.push((open, close));
                    covered_until = close.saturating_add(1);
                }
                None => {
                    file_test_only = true;
                    break;
                }
                Some(_) => {
                    // A `#![cfg(test)]` token inside a parenthesized or bracketed
                    // macro input is not an inner attribute on the surrounding
                    // function/module. Treating the nearest outer brace as its
                    // owner would let inert macro input hide later production
                    // output in that brace.
                }
            }
            continue;
        }

        let mut cursor = attribute_close + 1;
        let span_end = loop {
            let Some(token) = tokens.get(cursor) else {
                return Err(format!(
                    "#[cfg(test)] at byte {} has no following item",
                    tokens[index].start
                ));
            };
            if let Some(close) = macro_span_ends[cursor] {
                cursor = close + 1;
                continue;
            }
            match token.text {
                "(" | "[" => {
                    cursor = pairs[cursor]
                        .ok_or_else(|| "unpaired item-header delimiter".to_string())?
                        + 1;
                }
                "{" => {
                    break pairs[cursor]
                        .ok_or_else(|| "unpaired cfg(test) item body".to_string())?;
                }
                ";" => break cursor,
                "}" => {
                    return Err(format!(
                        "#[cfg(test)] at byte {} is not attached to a complete item",
                        tokens[index].start
                    ));
                }
                _ => cursor += 1,
            }
        };
        // `cfg(test)` applies to the item, not merely to attributes which
        // happen to follow it. Include every contiguous preceding outer
        // attribute so `#[path = "main.rs"] #[cfg(test)] mod diagnostics;`
        // cannot make the path attribute look production-reachable merely by
        // reversing otherwise-equivalent attribute order.
        let mut span_start = index;
        while span_start >= 2 && tokens[span_start - 1].text == "]" {
            let previous_close = span_start - 1;
            let Some(previous_open) = pairs[previous_close] else {
                break;
            };
            let Some(previous_hash) = previous_open.checked_sub(1) else {
                break;
            };
            if tokens[previous_open].text != "["
                || tokens[previous_hash].text != "#"
                || structural_mask[previous_hash]
            {
                break;
            }
            span_start = previous_hash;
        }
        spans.push((span_start, span_end));
        // Multiple contiguous cfg(test) attributes all govern the same item.
        // The first recognized shell covers the complete item, so revisiting
        // later shells would only rescan the same suffix quadratically.
        covered_until = span_end.saturating_add(1);
    }
    Ok((file_test_only, spans))
}

fn casual_function_scopes(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
    structural_mask: &[bool],
    test_mask: &[bool],
    macro_span_ends: &[Option<usize>],
) -> Vec<CasualFunctionScope> {
    let mut functions = Vec::new();
    for index in 0..tokens.len() {
        if structural_mask[index] || test_mask[index] || tokens[index].text != "fn" {
            continue;
        }
        let Some((name, mut cursor)) = tokens.get(index + 1).and_then(|token| {
            if token.text == "r"
                && tokens.get(index + 2).is_some_and(|token| token.text == "#")
                && tokens
                    .get(index + 3)
                    .is_some_and(|token| casual_token_is_identifier(token.text))
            {
                Some((tokens[index + 3].text, index + 4))
            } else {
                casual_token_is_identifier(token.text).then_some((token.text, index + 2))
            }
        }) else {
            continue;
        };
        while let Some(token) = tokens.get(cursor) {
            if let Some(close) = macro_span_ends[cursor] {
                // A braced macro invocation in a return type or where-clause
                // is token data, not the function body. Parenthesized and
                // bracketed invocations were already skipped by delimiter
                // pairing; handle the braced form explicitly as well.
                cursor = close + 1;
                continue;
            }
            match token.text {
                "(" | "[" => {
                    let Some(close) = pairs[cursor] else {
                        break;
                    };
                    cursor = close + 1;
                }
                "{" => {
                    let Some(close) = pairs[cursor] else {
                        break;
                    };
                    functions.push(CasualFunctionScope {
                        name: name.to_string(),
                        declaration_start: tokens[index].start,
                        open: cursor,
                        close,
                    });
                    break;
                }
                ";" | "}" => break,
                _ => cursor += 1,
            }
        }
    }
    functions
}

fn casual_test_mask(
    token_count: usize,
    file_test_only: bool,
    test_spans: &[(usize, usize)],
) -> Vec<bool> {
    if file_test_only {
        vec![true; token_count]
    } else {
        casual_span_mask(token_count, test_spans)
    }
}

fn casual_function_owner_map(
    token_count: usize,
    functions: &[CasualFunctionScope],
) -> Vec<Option<usize>> {
    let mut order = (0..functions.len()).collect::<Vec<_>>();
    order.sort_unstable_by_key(|index| functions[*index].open);
    let mut owners = vec![None; token_count];
    let mut stack = Vec::new();
    let mut next = 0usize;
    for token in 0..token_count {
        while stack
            .last()
            .is_some_and(|index| functions[*index].close <= token)
        {
            stack.pop();
        }
        while next < order.len() && functions[order[next]].open < token {
            stack.push(order[next]);
            next += 1;
        }
        owners[token] = stack.last().copied();
    }
    owners
}

fn casual_line_starts(source: &str) -> Vec<usize> {
    let mut starts = vec![0];
    starts.extend(
        source
            .bytes()
            .enumerate()
            .filter_map(|(index, byte)| (byte == b'\n').then_some(index + 1)),
    );
    starts
}

fn casual_line_for_byte(line_starts: &[usize], byte: usize) -> usize {
    line_starts.partition_point(|start| *start <= byte)
}

#[derive(Debug, Clone)]
struct CasualModuleEdge {
    target: String,
    module_dir: String,
    declaration: String,
}

#[derive(Debug, Clone)]
struct CasualModuleContext {
    path: String,
    package_root: String,
    module_dir: String,
    depth: usize,
}

struct CasualModuleFrame {
    context: CasualModuleContext,
    edges: Vec<CasualModuleEdge>,
}

#[derive(Debug)]
struct CasualCoreGraph {
    paths: BTreeSet<String>,
}

fn casual_join_portable(base: &str, child: &str) -> String {
    if base.is_empty() {
        child.to_string()
    } else {
        format!("{base}/{child}")
    }
}

fn casual_path_is_inside(path: &str, directory: &str) -> bool {
    path == directory
        || path
            .strip_prefix(directory)
            .is_some_and(|tail| tail.starts_with('/'))
}

fn casual_module_name(tokens: &[CasualRustToken<'_>], index: usize) -> Option<(String, usize)> {
    let token = tokens.get(index)?;
    if token.text == "r"
        && tokens.get(index + 1).is_some_and(|token| token.text == "#")
        && tokens
            .get(index + 2)
            .is_some_and(|token| casual_token_is_identifier(token.text))
    {
        Some((tokens[index + 2].text.to_string(), index + 3))
    } else {
        casual_token_is_identifier(token.text).then(|| (token.text.to_string(), index + 1))
    }
}

fn casual_module_header_start(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
    module: usize,
) -> usize {
    if module > 0 && tokens[module - 1].text == "pub" {
        return module - 1;
    }
    if module > 0
        && tokens[module - 1].text == ")"
        && let Some(open) = pairs[module - 1]
        && open > 0
        && tokens[open - 1].text == "pub"
    {
        return open - 1;
    }
    module
}

fn casual_attached_attributes(
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
    item_start: usize,
) -> Vec<usize> {
    let mut hashes = Vec::new();
    let mut cursor = item_start;
    while cursor > 0 && tokens[cursor - 1].text == "]" {
        let close = cursor - 1;
        let Some(bracket) = pairs[close] else {
            break;
        };
        let Some(hash) = bracket.checked_sub(1) else {
            break;
        };
        if tokens[bracket].text != "[" || tokens[hash].text != "#" {
            break;
        }
        hashes.push(hash);
        cursor = hash;
    }
    hashes.reverse();
    hashes
}

fn casual_module_path_attribute(
    path: &str,
    tokens: &[CasualRustToken<'_>],
    pairs: &[Option<usize>],
    hash: usize,
) -> Result<Option<String>, String> {
    let bracket = hash + 1;
    let close = pairs[bracket].ok_or_else(|| format!("unpaired attribute in {path}"))?;
    let Some((name, cursor)) = casual_identifier_at(tokens, bracket + 1) else {
        return Ok(None);
    };
    if name == "cfg_attr"
        && tokens[bracket + 2..close]
            .iter()
            .any(|token| token.text == "path")
    {
        return Err(format!(
            "{path} uses cfg_attr to select a module path; conditional path resolution is unsupported and therefore refused"
        ));
    }
    if name != "path" {
        return Ok(None);
    }
    if close != cursor + 2 || tokens.get(cursor).is_none_or(|token| token.text != "=") {
        return Err(format!("malformed #[path] module attribute in {path}"));
    }
    casual_path_literal(tokens[cursor + 1].text)
        .map(Some)
        .map_err(|error| format!("cannot decode #[path] in {path}: {error}"))
}

fn casual_module_edges(
    path: &str,
    package_root: &str,
    module_dir: &str,
    source: &str,
    sources: &BTreeMap<String, String>,
) -> Result<(Vec<CasualModuleEdge>, usize), String> {
    let tokens = casual_rust_tokens(source)?;
    let pairs = casual_delimiter_pairs(&tokens)?;
    let macro_spans = casual_macro_token_spans(&tokens, &pairs);
    let attributes = casual_attributes(&tokens, &pairs);
    let (_, structural_mask) =
        casual_structural_token_mask(tokens.len(), &macro_spans, &attributes);
    let macro_span_ends = casual_span_end_map(tokens.len(), &macro_spans);
    let containers = casual_delimiter_containers(&tokens);
    let (file_test_only, test_spans) = casual_cfg_test_spans(
        &tokens,
        &pairs,
        &structural_mask,
        &macro_span_ends,
        &containers,
    )?;
    let test_mask = casual_test_mask(tokens.len(), file_test_only, &test_spans);
    if file_test_only {
        return Ok((Vec::new(), tokens.len()));
    }

    let mut edges = Vec::new();
    let mut inline_modules = Vec::<(usize, String)>::new();
    for index in 0..tokens.len() {
        while inline_modules
            .last()
            .is_some_and(|(close, _)| *close < index)
        {
            inline_modules.pop();
        }
        if structural_mask[index] || test_mask[index] || tokens[index].text != "mod" {
            continue;
        }
        let Some((name, cursor)) = casual_module_name(&tokens, index + 1) else {
            continue;
        };
        let Some(body) = tokens.get(cursor).map(|token| token.text) else {
            return Err(format!("truncated module declaration in {path}"));
        };
        if !matches!(body, ";" | "{") {
            continue;
        }
        let header_start = casual_module_header_start(&tokens, &pairs, index);
        let attached = casual_attached_attributes(&tokens, &pairs, header_start);
        let mut direct_path = None;
        for hash in attached {
            if let Some(literal) = casual_module_path_attribute(path, &tokens, &pairs, hash)? {
                if direct_path.replace(literal).is_some() {
                    return Err(format!(
                        "multiple #[path] attributes on module {name} in {path}"
                    ));
                }
            }
        }
        let current_dir = inline_modules
            .last()
            .map_or(module_dir, |(_, directory)| directory.as_str());
        let child_identity_bytes = current_dir
            .len()
            .checked_add(if current_dir.is_empty() { 0 } else { 1 })
            .and_then(|length| length.checked_add(name.len()))
            .ok_or_else(|| format!("module identity length overflowed in {path}"))?;
        if child_identity_bytes > CASUAL_MAX_PORTABLE_PATH_BYTES {
            return Err(format!(
                "module {name} in {path} exceeds the {CASUAL_MAX_PORTABLE_PATH_BYTES}-byte portable-path cap"
            ));
        }
        let child_module_dir = casual_join_portable(current_dir, &name);
        if direct_path.is_some() && !inline_modules.is_empty() {
            return Err(format!(
                "module {name} in inline-module context in {path} carries #[path]; this path base is unsupported and refused"
            ));
        }
        if body == "{" {
            if direct_path.is_some() {
                return Err(format!(
                    "inline module {name} in {path} carries #[path]; this path-context combination is unsupported and refused"
                ));
            }
            let close =
                pairs[cursor].ok_or_else(|| format!("unpaired inline module {name} in {path}"))?;
            inline_modules.push((close, child_module_dir));
            continue;
        }

        let target = if let Some(literal) = direct_path {
            casual_normalized_inclusion(path, &literal)
                .map_err(|error| format!("cannot resolve #[path] in {path}: {error}"))?
        } else {
            let flat = format!("{child_module_dir}.rs");
            let directory = casual_join_portable(&child_module_dir, "mod.rs");
            match (
                sources.contains_key(&flat),
                sources.contains_key(&directory),
            ) {
                (true, false) => flat,
                (false, true) => directory,
                (true, true) => {
                    return Err(format!(
                        "ambiguous module {name} in {path}: both {flat} and {directory} exist"
                    ));
                }
                (false, false) => {
                    return Err(format!(
                        "missing module {name} in {path}: neither {flat} nor {directory} exists"
                    ));
                }
            }
        };
        if !target.ends_with(".rs") || !casual_path_is_inside(&target, package_root) {
            return Err(format!(
                "module {name} in {path} resolves outside package root {package_root}: {target}"
            ));
        }
        if !sources.contains_key(&target) {
            return Err(format!(
                "module {name} in {path} resolves to absent Rust source {target}"
            ));
        }
        edges.push(CasualModuleEdge {
            target,
            module_dir: child_module_dir,
            declaration: format!("{path}:mod {name}"),
        });
        if edges.len() > CASUAL_MAX_MODULE_EDGES {
            return Err(format!(
                "{path} exceeds the {CASUAL_MAX_MODULE_EDGES}-edge module audit cap"
            ));
        }
    }
    Ok((edges, tokens.len()))
}

fn casual_core_graph(
    sources: &BTreeMap<String, String>,
    roots: &[(String, String)],
) -> Result<CasualCoreGraph, String> {
    let mut paths = BTreeSet::new();
    let mut owners = BTreeMap::<String, String>::new();
    let mut root_contexts = Vec::new();
    for (path, package_root) in roots {
        if !sources.contains_key(path) {
            return Err(format!(
                "library root {path} is absent from the Rust-source inventory"
            ));
        }
        if owners
            .insert(
                path.clone(),
                format!("manifest library root {package_root}"),
            )
            .is_some()
        {
            return Err(format!("library root {path} is declared more than once"));
        }
        root_contexts.push(CasualModuleContext {
            path: path.clone(),
            package_root: package_root.clone(),
            module_dir: path
                .rsplit_once('/')
                .map_or_else(String::new, |(parent, _)| parent.to_string()),
            depth: 0,
        });
    }

    let mut total_tokens = 0usize;
    let mut edge_count = 0usize;
    let mut active = BTreeMap::<String, bool>::new();
    let mut frames = Vec::<CasualModuleFrame>::new();
    for root in root_contexts {
        active.insert(root.path.clone(), true);
        frames.push(casual_module_frame(
            root,
            sources,
            &mut paths,
            &mut total_tokens,
            &mut edge_count,
        )?);

        while !frames.is_empty() {
            let edge = frames
                .last_mut()
                .expect("non-empty graph stack has a current frame")
                .edges
                .pop();
            let Some(edge) = edge else {
                let completed = frames
                    .pop()
                    .expect("the observed graph frame remains present");
                active.insert(completed.context.path, false);
                continue;
            };
            match active.get(&edge.target).copied() {
                Some(true) => {
                    return Err(format!(
                        "cyclic module inclusion through {} reaches active ancestor {}",
                        edge.declaration, edge.target
                    ));
                }
                Some(false) => {
                    let previous = owners
                        .get(&edge.target)
                        .expect("completed module targets retain their owner");
                    return Err(format!(
                        "aliased module target {} is owned by both {previous} and {}",
                        edge.target, edge.declaration
                    ));
                }
                None => {}
            }
            if let Some(previous) = owners.insert(edge.target.clone(), edge.declaration.clone()) {
                return Err(format!(
                    "aliased module target {} is owned by both {previous} and {}",
                    edge.target, edge.declaration
                ));
            }
            let parent = frames
                .last()
                .expect("an outgoing edge retains its parent graph frame");
            let child = CasualModuleContext {
                path: edge.target,
                package_root: parent.context.package_root.clone(),
                module_dir: edge.module_dir,
                depth: parent.context.depth.saturating_add(1),
            };
            active.insert(child.path.clone(), true);
            frames.push(casual_module_frame(
                child,
                sources,
                &mut paths,
                &mut total_tokens,
                &mut edge_count,
            )?);
        }
    }
    Ok(CasualCoreGraph { paths })
}

fn casual_module_frame(
    context: CasualModuleContext,
    sources: &BTreeMap<String, String>,
    paths: &mut BTreeSet<String>,
    total_tokens: &mut usize,
    edge_count: &mut usize,
) -> Result<CasualModuleFrame, String> {
    if context.depth > CASUAL_MAX_MODULE_DEPTH {
        return Err(format!(
            "library graph exceeds the {CASUAL_MAX_MODULE_DEPTH}-level module-depth cap at {}",
            context.path
        ));
    }
    if context.path.len() > CASUAL_MAX_PORTABLE_PATH_BYTES
        || context.module_dir.len() > CASUAL_MAX_PORTABLE_PATH_BYTES
    {
        return Err(format!(
            "module identity for {} exceeds the {CASUAL_MAX_PORTABLE_PATH_BYTES}-byte portable-path cap",
            context.path
        ));
    }
    if !paths.insert(context.path.clone()) {
        return Err(format!(
            "module target {} entered the reachability walk more than once",
            context.path
        ));
    }
    if paths.len() > CASUAL_MAX_REACHABLE_SOURCES {
        return Err(format!(
            "library graph exceeds the {CASUAL_MAX_REACHABLE_SOURCES}-source reachability cap"
        ));
    }
    let source = sources
        .get(&context.path)
        .expect("module targets were inventory-checked before graph entry");
    let (edges, token_count) = casual_module_edges(
        &context.path,
        &context.package_root,
        &context.module_dir,
        source,
        sources,
    )?;
    *total_tokens = total_tokens
        .checked_add(token_count)
        .ok_or_else(|| "library token count overflowed".to_string())?;
    if *total_tokens > CASUAL_MAX_TOTAL_TOKENS {
        return Err(format!(
            "library graph exceeds the {CASUAL_MAX_TOTAL_TOKENS}-token audit cap"
        ));
    }
    *edge_count = edge_count
        .checked_add(edges.len())
        .ok_or_else(|| "module edge count overflowed".to_string())?;
    if *edge_count > CASUAL_MAX_MODULE_EDGES {
        return Err(format!(
            "library graph exceeds the {CASUAL_MAX_MODULE_EDGES}-edge audit cap"
        ));
    }
    Ok(CasualModuleFrame { context, edges })
}

fn casual_manifest_string(value: &str) -> Result<String, String> {
    let value = value.trim();
    let Some(delimiter) = value.as_bytes().first().copied() else {
        return Err("missing TOML string value".to_string());
    };
    if delimiter == b'\'' {
        let close = value[1..]
            .find('\'')
            .map(|relative| relative + 1)
            .ok_or_else(|| "unterminated literal TOML string".to_string())?;
        let tail = value[close + 1..].trim();
        if !tail.is_empty() && !tail.starts_with('#') {
            return Err("unexpected tokens after TOML string".to_string());
        }
        return Ok(value[1..close].to_string());
    }
    if delimiter != b'"' {
        return Err("expected a quoted TOML string".to_string());
    }
    let mut escaped = false;
    let mut close = None;
    for (relative, byte) in value.as_bytes()[1..].iter().copied().enumerate() {
        if escaped {
            escaped = false;
        } else if byte == b'\\' {
            escaped = true;
        } else if byte == b'"' {
            close = Some(relative + 1);
            break;
        }
    }
    let close = close.ok_or_else(|| "unterminated basic TOML string".to_string())?;
    let tail = value[close + 1..].trim();
    if !tail.is_empty() && !tail.starts_with('#') {
        return Err("unexpected tokens after TOML string".to_string());
    }
    casual_path_literal(&value[..=close])
}

fn casual_manifest_key_path(input: &str) -> Result<Vec<String>, String> {
    let mut cursor = 0usize;
    let bytes = input.as_bytes();
    let mut keys = Vec::new();
    loop {
        while bytes
            .get(cursor)
            .is_some_and(|byte| byte.is_ascii_whitespace())
        {
            cursor += 1;
        }
        let Some(first) = bytes.get(cursor).copied() else {
            return Err("empty TOML key path".to_string());
        };
        let key = if matches!(first, b'\'' | b'"') {
            let delimiter = first;
            let start = cursor;
            cursor += 1;
            let mut escaped = false;
            let close = loop {
                let Some(byte) = bytes.get(cursor).copied() else {
                    return Err("unterminated quoted TOML key".to_string());
                };
                if delimiter == b'"' && escaped {
                    escaped = false;
                } else if delimiter == b'"' && byte == b'\\' {
                    escaped = true;
                } else if byte == delimiter {
                    break cursor;
                }
                cursor += 1;
            };
            cursor += 1;
            if delimiter == b'\'' {
                input[start + 1..close].to_string()
            } else {
                casual_path_literal(&input[start..cursor])
                    .map_err(|error| format!("invalid quoted TOML key: {error}"))?
            }
        } else {
            let start = cursor;
            while bytes
                .get(cursor)
                .is_some_and(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-'))
            {
                cursor += 1;
            }
            if cursor == start {
                return Err(format!(
                    "unsupported TOML key token starting at {:?}",
                    &input[start..]
                ));
            }
            input[start..cursor].to_string()
        };
        keys.push(key);

        while bytes
            .get(cursor)
            .is_some_and(|byte| byte.is_ascii_whitespace())
        {
            cursor += 1;
        }
        match bytes.get(cursor) {
            None => return Ok(keys),
            Some(b'.') => cursor += 1,
            Some(_) => {
                return Err(format!(
                    "unexpected TOML key-path suffix {:?}",
                    &input[cursor..]
                ));
            }
        }
    }
}

fn casual_manifest_key_is(actual: &[String], expected: &[&str]) -> bool {
    actual.len() == expected.len()
        && actual
            .iter()
            .zip(expected)
            .all(|(actual, expected)| actual.as_str() == *expected)
}

fn casual_manifest_comment_free(input: &str) -> Result<&str, String> {
    let bytes = input.as_bytes();
    let mut quote = None;
    let mut escaped = false;
    for (index, byte) in bytes.iter().copied().enumerate() {
        match quote {
            Some(b'"') if escaped => escaped = false,
            Some(b'"') if byte == b'\\' => escaped = true,
            Some(delimiter) if byte == delimiter => quote = None,
            Some(_) => {}
            None if matches!(byte, b'\'' | b'"') => quote = Some(byte),
            None if byte == b'#' => return Ok(&input[..index]),
            None => {}
        }
    }
    if quote.is_some() {
        return Err("unterminated quote in TOML key/table syntax".to_string());
    }
    Ok(input)
}

fn casual_manifest_table_path(line: &str) -> Result<Option<Vec<String>>, String> {
    let header = casual_manifest_comment_free(line)?.trim();
    if !header.starts_with('[') {
        return Ok(None);
    }
    let (inner, array) = if header.starts_with("[[") {
        if !header.ends_with("]]") {
            return Err("unterminated TOML array-table header".to_string());
        }
        (&header[2..header.len() - 2], true)
    } else {
        if !header.ends_with(']') {
            return Err("unterminated TOML table header".to_string());
        }
        (&header[1..header.len() - 1], false)
    };
    let path = casual_manifest_key_path(inner)?;
    if array
        && path
            .first()
            .is_some_and(|key| matches!(key.as_str(), "package" | "lib"))
    {
        return Err("package/lib must be TOML tables, not arrays of tables".to_string());
    }
    Ok(Some(path))
}

fn casual_manifest_assignment(line: &str) -> Result<Option<(Vec<String>, &str)>, String> {
    let bytes = line.as_bytes();
    let mut quote = None;
    let mut escaped = false;
    for (index, byte) in bytes.iter().copied().enumerate() {
        match quote {
            Some(b'"') if escaped => escaped = false,
            Some(b'"') if byte == b'\\' => escaped = true,
            Some(delimiter) if byte == delimiter => quote = None,
            Some(_) => {}
            None if matches!(byte, b'\'' | b'"') => quote = Some(byte),
            None if byte == b'#' => return Ok(None),
            None if byte == b'=' => {
                let keys = casual_manifest_key_path(line[..index].trim())?;
                return Ok(Some((keys, &line[index + 1..])));
            }
            None => {}
        }
    }
    if quote.is_some() {
        return Err("unterminated quote in TOML assignment key".to_string());
    }
    Ok(None)
}

fn casual_manifest_library_root(
    package_root: &str,
    manifest: &str,
    sources: &BTreeMap<String, String>,
) -> Result<Option<String>, String> {
    let mut section = Vec::<String>::new();
    let mut package_seen = false;
    let mut autolib = true;
    let mut lib_seen = false;
    let mut lib_path = None;
    for (line_number, raw) in manifest.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(table) = casual_manifest_table_path(line).map_err(|error| {
            format!(
                "cannot parse {package_root}/Cargo.toml:{} table: {error}",
                line_number + 1
            )
        })? {
            section = table;
            if casual_manifest_key_is(&section, &["package"]) {
                package_seen = true;
            } else if casual_manifest_key_is(&section, &["lib"]) {
                if lib_seen {
                    return Err(format!(
                        "{package_root}/Cargo.toml declares [lib] more than once"
                    ));
                }
                lib_seen = true;
            }
            continue;
        }
        let Some((key, value)) = casual_manifest_assignment(line).map_err(|error| {
            format!(
                "cannot parse {package_root}/Cargo.toml:{} assignment: {error}",
                line_number + 1
            )
        })?
        else {
            continue;
        };
        let mut full_key = section.clone();
        full_key.extend(key);
        if full_key.first().is_some_and(|key| key == "package") {
            package_seen = true;
        }
        if full_key.first().is_some_and(|key| key == "lib") {
            lib_seen = true;
        }
        if casual_manifest_key_is(&full_key, &["package"])
            || casual_manifest_key_is(&full_key, &["lib"])
        {
            return Err(format!(
                "{package_root}/Cargo.toml:{} uses an inline package/lib table; this audit requires an explicit or dotted table so target authority stays reviewable",
                line_number + 1
            ));
        }
        if casual_manifest_key_is(&full_key, &["package", "autolib"]) {
            autolib = match value
                .split_once('#')
                .map_or(value, |(value, _)| value)
                .trim()
            {
                "true" => true,
                "false" => false,
                _ => {
                    return Err(format!(
                        "{package_root}/Cargo.toml:{} has non-boolean package.autolib",
                        line_number + 1
                    ));
                }
            };
        } else if casual_manifest_key_is(&full_key, &["lib", "path"]) {
            if lib_path.is_some() {
                return Err(format!(
                    "{package_root}/Cargo.toml declares lib.path more than once"
                ));
            }
            lib_path = Some(casual_manifest_string(value).map_err(|error| {
                format!(
                    "cannot parse {package_root}/Cargo.toml:{} lib.path: {error}",
                    line_number + 1
                )
            })?);
        }
    }
    if !package_seen {
        return Ok(None);
    }
    let default = casual_join_portable(package_root, "src/lib.rs");
    let target = if lib_seen {
        if let Some(literal) = lib_path {
            casual_normalized_inclusion(&format!("{package_root}/Cargo.toml"), &literal)
                .map_err(|error| format!("cannot resolve {package_root} lib.path: {error}"))?
        } else {
            default
        }
    } else if autolib && sources.contains_key(&default) {
        default
    } else {
        return Ok(None);
    };
    if !target.ends_with(".rs") || !casual_path_is_inside(&target, package_root) {
        return Err(format!(
            "{package_root} library target escapes its package or is not Rust source: {target}"
        ));
    }
    if !sources.contains_key(&target) {
        return Err(format!(
            "{package_root} library target {target} is absent from the Rust-source inventory"
        ));
    }
    Ok(Some(target))
}

fn casual_workspace_library_roots(
    root: &Path,
    sources: &BTreeMap<String, String>,
) -> Result<Vec<(String, String)>, String> {
    let crates = root.join("crates");
    let mut packages = std::fs::read_dir(&crates)
        .map_err(|error| format!("cannot read {}: {error}", crates.display()))?
        .map(|entry| {
            entry
                .map(|entry| entry.path())
                .map_err(|error| format!("cannot enumerate {}: {error}", crates.display()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    packages.sort();
    let mut roots = Vec::new();
    for package in packages {
        if !package.is_dir() || !package.join("Cargo.toml").is_file() {
            continue;
        }
        let package_root = casual_inventory_relative(root, &package)?;
        let manifest_path = package.join("Cargo.toml");
        let manifest =
            casual_read_bounded_utf8(&manifest_path, CASUAL_MAX_MANIFEST_BYTES, "Cargo manifest")
                .map_err(|error| format!("cannot read {package_root}/Cargo.toml: {error}"))?;
        if let Some(target) = casual_manifest_library_root(&package_root, &manifest, sources)? {
            roots.push((target, package_root));
        }
    }
    Ok(roots)
}

fn casual_fixture_library_roots(sources: &BTreeMap<String, String>) -> Vec<(String, String)> {
    sources
        .keys()
        .filter(|path| path.starts_with("crates/") && path.ends_with("/src/lib.rs"))
        .filter_map(|path| {
            path.strip_suffix("/src/lib.rs")
                .map(|package| (path.clone(), package.to_string()))
        })
        .collect()
}

fn protected_print_macro(name: &str) -> bool {
    matches!(name, "print" | "println" | "eprint" | "eprintln" | "dbg")
}

fn casual_raw_identifier_start(tokens: &[CasualRustToken<'_>], identifier: usize) -> usize {
    if identifier >= 2 && tokens[identifier - 2].text == "r" && tokens[identifier - 1].text == "#" {
        identifier - 2
    } else {
        identifier
    }
}

fn casual_macro_path_start(tokens: &[CasualRustToken<'_>], identifier: usize) -> usize {
    let mut start = casual_raw_identifier_start(tokens, identifier);
    loop {
        if start < 2 || tokens[start - 2].text != ":" || tokens[start - 1].text != ":" {
            break;
        }
        let separator = start - 2;
        if separator == 0 {
            start = 0;
            break;
        }
        let previous = separator - 1;
        if !casual_token_is_identifier(tokens[previous].text) {
            start = separator;
            break;
        }
        start = casual_raw_identifier_start(tokens, previous);
    }
    start
}

fn casual_protected_name_binding(
    tokens: &[CasualRustToken<'_>],
    use_statement_mask: &[bool],
    index: usize,
) -> Option<&'static str> {
    let name_start = casual_raw_identifier_start(tokens, index);
    if name_start > 0 && tokens[name_start - 1].text == "as" {
        return Some("import bound to protected name");
    }
    if tokens
        .get(index + 1)
        .is_some_and(|token| token.text == "as")
    {
        return Some("protected name renamed by import");
    }
    if name_start >= 2
        && tokens[name_start - 2].text == "macro_rules"
        && tokens[name_start - 1].text == "!"
    {
        return Some("local macro declaration uses protected name");
    }
    if name_start >= 1 && tokens[name_start - 1].text == "macro" {
        return Some("local macro declaration uses protected name");
    }
    if name_start >= 2
        && tokens[name_start - 2].text == "extern"
        && tokens[name_start - 1].text == "crate"
    {
        return Some("extern crate binds a protected macro-capable name");
    }

    if use_statement_mask[index] {
        return Some("import binds or resolves through protected name");
    }
    None
}

fn casual_use_statement_mask(
    tokens: &[CasualRustToken<'_>],
    structural_mask: &[bool],
) -> Vec<bool> {
    let mut mask = vec![false; tokens.len()];
    let mut in_use = false;
    for index in 0..tokens.len() {
        if !structural_mask[index] && tokens[index].text == "use" {
            in_use = true;
        }
        mask[index] = in_use;
        if in_use && tokens[index].text == ";" {
            in_use = false;
        }
    }
    mask
}

fn casual_macro_authority_hazards(
    tokens: &[CasualRustToken<'_>],
    attributes: &[CasualAttribute],
    structural_mask: &[bool],
    test_mask: &[bool],
    line_starts: &[usize],
) -> Result<Vec<(usize, String)>, String> {
    let mut hazards = Vec::new();
    for attribute in attributes {
        if structural_mask[attribute.hash] || test_mask[attribute.hash] {
            continue;
        }
        let interior = &tokens[attribute.bracket + 1..attribute.close];
        let attribute_name = casual_identifier_at(interior, 0).map(|(name, _)| name);
        if attribute.inner {
            let can_remove_std = attribute_name
                .is_some_and(|name| matches!(name, "no_std" | "no_core"))
                || (attribute_name == Some("cfg_attr")
                    && interior
                        .iter()
                        .any(|token| matches!(token.text, "no_std" | "no_core")));
            if can_remove_std {
                casual_push_authority_hazard(
                    &mut hazards,
                    casual_line_for_byte(line_starts, tokens[attribute.hash].start),
                    "an allowlisted absolute ::std macro requires the compiler-provided std extern-prelude binding; no_std/no_core can rebind that authority".to_string(),
                )?;
            }
            continue;
        }
        let imports_macros = attribute_name == Some("macro_use")
            || (attribute_name == Some("cfg_attr")
                && interior.iter().any(|token| token.text == "macro_use"));
        if imports_macros {
            casual_push_authority_hazard(
                &mut hazards,
                casual_line_for_byte(line_starts, tokens[attribute.hash].start),
                "legacy or conditional #[macro_use] can inject a protected macro name".to_string(),
            )?;
        }
    }

    let mut index = 0usize;
    while index < tokens.len() {
        if structural_mask[index] || test_mask[index] || tokens[index].text != "use" {
            index += 1;
            continue;
        }
        let mut cursor = index + 1;
        let mut has_glob = false;
        while let Some(token) = tokens.get(cursor) {
            if token.text == ";" {
                break;
            }
            has_glob |= token.text == "*";
            cursor += 1;
        }
        if tokens.get(cursor).is_none_or(|token| token.text != ";") {
            return Err(format!(
                "unterminated use declaration at byte {}",
                tokens[index].start
            ));
        }
        if !has_glob {
            index = cursor + 1;
            continue;
        }
        let local_root = tokens
            .get(index + 1)
            .is_some_and(|token| matches!(token.text, "crate" | "self" | "super"));
        if !local_root {
            casual_push_authority_hazard(
                &mut hazards,
                casual_line_for_byte(line_starts, tokens[index].start),
                "potentially external glob import can inject a protected macro name".to_string(),
            )?;
        }
        index = cursor + 1;
    }
    Ok(hazards)
}

fn casual_push_authority_hazard(
    hazards: &mut Vec<(usize, String)>,
    line: usize,
    detail: String,
) -> Result<(), String> {
    if hazards.len() == CASUAL_MAX_DIAGNOSTICS {
        return Err(format!(
            "macro authority analysis exceeds the {CASUAL_MAX_DIAGNOSTICS}-hazard cap"
        ));
    }
    hazards.push((line, detail));
    Ok(())
}

fn scan_casual_print_source(
    path: &str,
    source: &str,
    explicitly_core: bool,
) -> Result<CasualPrintScan, String> {
    if !explicitly_core {
        return Ok(CasualPrintScan::default());
    }

    let tokens = casual_rust_tokens(source)?;
    let pairs = casual_delimiter_pairs(&tokens)?;
    let macro_token_spans = casual_macro_token_spans(&tokens, &pairs);
    let attributes = casual_attributes(&tokens, &pairs);
    let (macro_mask, structural_mask) =
        casual_structural_token_mask(tokens.len(), &macro_token_spans, &attributes);
    let macro_span_ends = casual_span_end_map(tokens.len(), &macro_token_spans);
    let containers = casual_delimiter_containers(&tokens);
    let (file_test_only, test_spans) = casual_cfg_test_spans(
        &tokens,
        &pairs,
        &structural_mask,
        &macro_span_ends,
        &containers,
    )?;
    let test_mask = casual_test_mask(tokens.len(), file_test_only, &test_spans);
    let functions = casual_function_scopes(
        &tokens,
        &pairs,
        &structural_mask,
        &test_mask,
        &macro_span_ends,
    );
    let function_owners = casual_function_owner_map(tokens.len(), &functions);
    let use_statement_mask = casual_use_statement_mask(&tokens, &structural_mask);
    let line_starts = casual_line_starts(source);
    let authority_hazards = casual_macro_authority_hazards(
        &tokens,
        &attributes,
        &structural_mask,
        &test_mask,
        &line_starts,
    )?;
    let mut occurrences = Vec::new();
    for (index, token) in tokens.iter().enumerate() {
        if test_mask[index] {
            continue;
        }

        let owner = function_owners[index].map(|owner| &functions[owner]);
        let owner_name = owner.map_or("<module>", |function| function.name.as_str());
        let owner_start = owner.map(|function| function.declaration_start);
        let line = casual_line_for_byte(&line_starts, token.start);

        if protected_print_macro(token.text)
            && let Some(binding) =
                casual_protected_name_binding(&tokens, &use_statement_mask, index)
        {
            if occurrences.len() == CASUAL_MAX_DIAGNOSTICS {
                return Err(format!(
                    "protected-output occurrence count exceeds the {CASUAL_MAX_DIAGNOSTICS}-occurrence cap"
                ));
            }
            occurrences.push(CasualPrintOccurrence {
                path: path.to_string(),
                owner: owner_name.to_string(),
                owner_start,
                macro_name: token.text.to_string(),
                invocation_anchor: format!("protected-binding:{binding}:{}", token.text),
                alias_of: Some(binding.to_string()),
                in_macro_tokens: structural_mask[index],
                line,
            });
        }

        if !protected_print_macro(token.text)
            || !tokens.get(index + 1).is_some_and(|next| next.text == "!")
        {
            continue;
        }
        let open = index + 2;
        if !tokens
            .get(open)
            .is_some_and(|candidate| matches!(candidate.text, "(" | "[" | "{"))
        {
            return Err(format!(
                "protected macro {}! at byte {} lacks a balanced token tree",
                token.text, token.start
            ));
        }
        let close = pairs[open]
            .ok_or_else(|| format!("protected macro {}! has an unpaired body", token.text))?;
        let invocation_start = casual_macro_path_start(&tokens, index);
        let invocation_anchor = tokens[invocation_start..=close]
            .iter()
            .map(|candidate| candidate.text)
            .collect::<String>();
        if occurrences.len() == CASUAL_MAX_DIAGNOSTICS {
            return Err(format!(
                "protected-output occurrence count exceeds the {CASUAL_MAX_DIAGNOSTICS}-occurrence cap"
            ));
        }
        occurrences.push(CasualPrintOccurrence {
            path: path.to_string(),
            owner: owner_name.to_string(),
            owner_start,
            macro_name: token.text.to_string(),
            invocation_anchor,
            alias_of: None,
            in_macro_tokens: macro_mask[index] || structural_mask[index],
            line,
        });
    }
    let mut function_starts = BTreeMap::<String, Vec<usize>>::new();
    for function in functions {
        function_starts
            .entry(function.name)
            .or_default()
            .push(function.declaration_start);
    }
    Ok(CasualPrintScan {
        occurrences,
        function_starts,
        authority_hazards,
    })
}

fn casual_push_violation(violations: &mut Vec<Violation>, violation: Violation) -> bool {
    if violations.len() >= CASUAL_MAX_DIAGNOSTICS - 1 {
        if violations.len() == CASUAL_MAX_DIAGNOSTICS - 1 {
            violations.push(Violation {
                check: CASUAL_PRINT_CHECK,
                crate_name: "<repo>".to_string(),
                detail: format!(
                    "casual-print diagnostics reached the {CASUAL_MAX_DIAGNOSTICS}-record cap; fix the reported prefix and rerun"
                ),
            });
        }
        false
    } else {
        violations.push(violation);
        true
    }
}

fn audit_casual_print_sources_with_roots(
    sources: &BTreeMap<String, String>,
    roots: &[(String, String)],
) -> Vec<Violation> {
    let mut violations = Vec::new();
    let core_graph = match casual_core_graph(sources, roots) {
        Ok(graph) => graph,
        Err(error) => {
            violations.push(Violation {
                check: CASUAL_PRINT_CHECK,
                crate_name: "<repo>".to_string(),
                detail: format!("cannot exhaustively resolve library module graph: {error}"),
            });
            return violations;
        }
    };
    for path in &core_graph.paths {
        let source = sources
            .get(path)
            .expect("core graph contains inventory-owned sources");
        let scan = match scan_casual_print_source(path, source, true) {
            Ok(scan) => scan,
            Err(error) => {
                if !casual_push_violation(
                    &mut violations,
                    Violation {
                        check: CASUAL_PRINT_CHECK,
                        crate_name: path.clone(),
                        detail: format!("cannot exhaustively scan {path}: {error}"),
                    },
                ) {
                    return violations;
                }
                continue;
            }
        };
        let CasualPrintScan {
            occurrences,
            function_starts,
            authority_hazards,
        } = scan;
        let mut allowed = vec![false; occurrences.len()];
        for allowance in CASUAL_PRINT_ALLOWLIST
            .iter()
            .filter(|allowance| allowance.path == path)
        {
            let candidates = occurrences
                .iter()
                .enumerate()
                .filter(|(_, occurrence)| {
                    occurrence.owner == allowance.owner
                        && occurrence.alias_of.is_none()
                        && !occurrence.in_macro_tokens
                })
                .collect::<Vec<_>>();
            let mut owner_starts = function_starts
                .get(allowance.owner)
                .cloned()
                .unwrap_or_default();
            owner_starts.sort_unstable();
            owner_starts.dedup();
            let actual = candidates
                .iter()
                .map(|(_, occurrence)| occurrence.invocation_anchor.as_str())
                .collect::<Vec<_>>();
            let expected = allowance.invocation_anchors.to_vec();
            let candidates_belong_to_unique_owner = owner_starts.len() == 1
                && candidates
                    .iter()
                    .all(|(_, occurrence)| occurrence.owner_start == owner_starts.first().copied());
            let authority_is_unambiguous = authority_hazards.is_empty();
            if candidates_belong_to_unique_owner && authority_is_unambiguous && actual == expected {
                for (index, _) in candidates {
                    allowed[index] = true;
                }
            } else {
                if !casual_push_violation(
                    &mut violations,
                    Violation {
                        check: CASUAL_PRINT_CHECK,
                        crate_name: path.clone(),
                        detail: format!(
                            "{}: allowance for unique fn {} ({}) expected exact invocations {expected:?}, observed {} owner body/bodies, {actual:?}, and {} macro-authority hazard(s); update code and ratchet together",
                            allowance.path,
                            allowance.owner,
                            allowance.reason,
                            owner_starts.len(),
                            authority_hazards.len(),
                        ),
                    },
                ) {
                    return violations;
                }
            }
        }
        if CASUAL_PRINT_ALLOWLIST
            .iter()
            .any(|allowance| allowance.path == path)
        {
            for (line, hazard) in authority_hazards {
                if !casual_push_violation(
                    &mut violations,
                    Violation {
                        check: CASUAL_PRINT_CHECK,
                        crate_name: path.clone(),
                        detail: format!(
                            "{path}:{line}: {hazard}; an allowlisted emitter may not inherit ambiguous macro authority"
                        ),
                    },
                ) {
                    return violations;
                }
            }
        }
        for (index, occurrence) in occurrences.into_iter().enumerate() {
            if allowed[index] {
                continue;
            }
            let detail = if let Some(binding) = occurrence.alias_of {
                format!(
                    "{}:{}: {binding} `{}` can change protected macro identity and bypass the typed-output ratchet; core libraries may not import, rename, or declare protected output names",
                    occurrence.path, occurrence.line, occurrence.macro_name,
                )
            } else if occurrence.in_macro_tokens {
                format!(
                    "{}:{}: {}! is spelled inside macro token input/body or an attribute interior, where nested attributes and apparent function owners are non-authoritative; return a typed fs-obs event/record outside token generation",
                    occurrence.path, occurrence.line, occurrence.macro_name,
                )
            } else {
                format!(
                    "{}:{}: {}! in fn {} creates an untyped process-output path; return a typed fs-obs event/record and let a CLI or process runner render it",
                    occurrence.path, occurrence.line, occurrence.macro_name, occurrence.owner,
                )
            };
            if !casual_push_violation(
                &mut violations,
                Violation {
                    check: CASUAL_PRINT_CHECK,
                    crate_name: occurrence.path.clone(),
                    detail,
                },
            ) {
                return violations;
            }
        }
    }
    violations
}

fn audit_casual_print_sources(sources: &BTreeMap<String, String>) -> Vec<Violation> {
    let roots = casual_fixture_library_roots(sources);
    audit_casual_print_sources_with_roots(sources, &roots)
}

fn check_casual_print(root: &Path) -> Vec<Violation> {
    let sources = match workspace_rust_sources(root) {
        Ok(sources) => sources,
        Err(detail) => {
            return vec![Violation {
                check: CASUAL_PRINT_CHECK,
                crate_name: "<repo>".to_string(),
                detail,
            }];
        }
    };
    let roots = match casual_workspace_library_roots(root, &sources) {
        Ok(roots) => roots,
        Err(detail) => {
            return vec![Violation {
                check: CASUAL_PRINT_CHECK,
                crate_name: "<repo>".to_string(),
                detail: format!("cannot determine package library targets: {detail}"),
            }];
        }
    };
    if roots.len() > CASUAL_MAX_REACHABLE_SOURCES {
        return vec![Violation {
            check: CASUAL_PRINT_CHECK,
            crate_name: "<repo>".to_string(),
            detail: format!(
                "library-root count exceeds the {CASUAL_MAX_REACHABLE_SOURCES}-root audit cap"
            ),
        }];
    }
    audit_casual_print_sources_with_roots(&sources, &roots)
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

fn emit(
    violations: &[Violation],
    dev_notes: &[(String, String)],
    policy_notes: &[PolicyNote],
    checks_run: &[&str],
    crates_checked: usize,
) -> ExitCode {
    // JSON-lines verdicts on stdout (agent-facing).
    for (krate, dep) in dev_notes {
        println!(
            "{{\"check\":\"dev-dependency-note\",\"crate\":\"{}\",\"verdict\":\"note\",\"detail\":\"external dev-dependency {} (exempt; must be an isolated, documented oracle)\"}}",
            json_escape(krate),
            json_escape(dep)
        );
    }
    for note in policy_notes {
        println!(
            "{{\"check\":\"{}\",\"crate\":\"{}\",\"verdict\":\"{}\",\"detail\":\"{}\"}}",
            json_escape(note.check),
            json_escape(&note.crate_name),
            json_escape(note.verdict),
            json_escape(&note.detail)
        );
    }
    for v in violations {
        println!(
            "{{\"check\":\"{}\",\"crate\":\"{}\",\"verdict\":\"violation\",\"detail\":\"{}\"}}",
            json_escape(v.check),
            json_escape(&v.crate_name),
            json_escape(&v.detail)
        );
    }
    println!(
        "{{\"check\":\"summary\",\"checks\":\"{}\",\"crates\":{},\"violations\":{}}}",
        json_escape(&checks_run.join("+")),
        crates_checked,
        violations.len()
    );
    // Human-facing summary on stderr.
    if violations.is_empty() {
        eprintln!(
            "policy OK: {} crates, checks: {}",
            crates_checked,
            checks_run.join(", ")
        );
        ExitCode::SUCCESS
    } else {
        for v in violations {
            eprintln!("VIOLATION [{}] {}: {}", v.check, v.crate_name, v.detail);
        }
        eprintln!("{} violation(s)", violations.len());
        ExitCode::FAILURE
    }
}

/// Golden-coupling discipline (bead y4pt): every golden hash declares
/// the upstream semantic surfaces it was frozen against in
/// golden-couplings.json; a surface whose source version const drifts
/// from the registry — or a golden pinned against a stale surface
/// version — fails with a pointer to every row that must be
/// deliberately re-frozen (per docs/GOLDEN_POLICY.md).
#[allow(clippy::too_many_lines)] // registry parse + two rule passes: the check IS the semantics
fn check_goldens(root: &Path) -> Vec<Violation> {
    // One-line-per-entry registry; extract string/number fields by key.
    fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
        let tag = format!("\"{key}\"");
        let start = line.find(&tag)? + tag.len();
        let rest = line[start..].trim_start().strip_prefix(':')?.trim_start();
        if let Some(stripped) = rest.strip_prefix('"') {
            stripped.split('"').next()
        } else {
            rest.split([',', '}']).next().map(str::trim)
        }
    }
    let mut violations = Vec::new();
    let bail = |detail: String| Violation {
        check: "golden-couplings",
        crate_name: "<repo>".to_string(),
        detail,
    };
    let Ok(registry) = std::fs::read_to_string(root.join("golden-couplings.json")) else {
        violations.push(bail(
            "golden-couplings.json missing at workspace root".to_string(),
        ));
        return violations;
    };
    let external_versions = if root.join("identity-authorities.json").is_file() {
        match identities::external_coupling_versions(root) {
            Ok(versions) => Some(versions),
            Err(errors) => {
                violations.extend(errors.into_iter().map(|error| {
                    bail(format!(
                        "external identity coupling {}: {}",
                        error.crate_name, error.detail
                    ))
                }));
                None
            }
        }
    } else {
        None
    };
    // Surfaces: id -> (registry version, dependents filled later).
    let mut surface_versions: Vec<(String, u32)> = Vec::new();
    let mut in_surfaces = false;
    let mut in_goldens = false;
    let mut goldens: Vec<(String, String, String, String)> = Vec::new();
    for line in registry.lines() {
        if line.starts_with("\"surfaces\"") {
            in_surfaces = true;
            in_goldens = false;
            continue;
        }
        if line.starts_with("\"goldens\"") {
            in_goldens = true;
            in_surfaces = false;
            continue;
        }
        if in_surfaces && line.trim_start().starts_with('{') {
            let (Some(id), Some(file), Some(ver)) = (
                field(line, "id"),
                field(line, "file"),
                field(line, "version"),
            ) else {
                violations.push(bail(format!("malformed surface row: {line}")));
                continue;
            };
            let Ok(reg_ver) = ver.parse::<u32>() else {
                violations.push(bail(format!("surface {id}: bad version {ver:?}")));
                continue;
            };
            if let Some(symbol) = field(line, "symbol") {
                if let Some(versions) = &external_versions
                    && versions.get(id) != Some(&reg_ver)
                {
                    violations.push(bail(format!(
                        "external surface {id}: {file}#{symbol} v{reg_ver} is not an exact validated identity authority coupling"
                    )));
                }
                surface_versions.push((id.to_string(), reg_ver));
                continue;
            }
            let Some(name) = field(line, "const") else {
                violations.push(bail(format!("malformed surface row: {line}")));
                continue;
            };
            let needle = format!("pub const {name}: u32 = ");
            let src = std::fs::read_to_string(root.join(file)).unwrap_or_default();
            let Some(actual) = src
                .find(&needle)
                .and_then(|at| src[at + needle.len()..].split(';').next())
                .and_then(|v| v.trim().parse::<u32>().ok())
            else {
                violations.push(bail(format!(
                    "surface {id}: {file} does not declare `{needle}<version>;`"
                )));
                continue;
            };
            if actual != reg_ver {
                let dependents: Vec<&str> = registry
                    .lines()
                    .filter(|l| l.contains("\"golden\"") && l.contains(id))
                    .filter_map(|l| field(l, "golden"))
                    .collect();
                violations.push(bail(format!(
                    "surface {id} version drifted: source declares {actual}, registry pins \
                    {reg_ver} — an upstream semantic change must deliberately re-freeze its \
                     dependents {dependents:?} (docs/GOLDEN_POLICY.md), then update both pins"
                )));
            }
            if let Some(domain_const) = field(line, "domain_const") {
                let Some(expected_domain) = field(line, "domain") else {
                    violations.push(bail(format!(
                        "surface {id}: domain_const {domain_const} requires a domain literal"
                    )));
                    continue;
                };
                let (domain_file, domain_symbol) =
                    if let Some((path, symbol)) = domain_const.split_once('#') {
                        let path = Path::new(path);
                        if symbol.is_empty()
                            || path.is_absolute()
                            || path.components().any(|component| {
                                !matches!(component, std::path::Component::Normal(_))
                            })
                        {
                            violations.push(bail(format!(
                                "surface {id}: domain_const {domain_const:?} is not a safe \
                                 repo-relative path#symbol reference"
                            )));
                            continue;
                        }
                        (path.to_path_buf(), symbol)
                    } else {
                        (PathBuf::from(file), domain_const)
                    };
                let domain_src =
                    std::fs::read_to_string(root.join(&domain_file)).unwrap_or_default();
                let domain_needle = format!("const {domain_symbol}: &str =");
                let actual_domain = domain_src
                    .find(&domain_needle)
                    .and_then(|at| {
                        domain_src[at + domain_needle.len()..]
                            .trim_start()
                            .strip_prefix('"')
                    })
                    .and_then(|rest| rest.split('"').next());
                if actual_domain != Some(expected_domain) {
                    violations.push(bail(format!(
                        "surface {id} domain drifted: {} must declare \
                         `{domain_needle} \"{expected_domain}\";` — domain rotations require a \
                         semantic version bump and deliberate golden-coupling annotation",
                        domain_file.display(),
                    )));
                }
                let explicit_version = expected_domain == "fsid"
                    || expected_domain
                        .split(['.', ':', '-'])
                        .any(|segment| segment == format!("v{reg_ver}"));
                if !explicit_version {
                    violations.push(bail(format!(
                        "surface {id}: domain {expected_domain:?} must carry an exact \
                         v{reg_ver} dot/colon segment; schema/domain \
                         versions may not drift independently"
                    )));
                }
            }
            surface_versions.push((id.to_string(), reg_ver));
        }
        if in_goldens && line.trim_start().starts_with('{') {
            let (Some(g), Some(file), Some(name), Some(deps)) = (
                field(line, "golden"),
                field(line, "file"),
                field(line, "const"),
                field(line, "depends_on"),
            ) else {
                violations.push(bail(format!("malformed golden row: {line}")));
                continue;
            };
            if field(line, "justification").is_none_or(|j| j.len() < 20) {
                violations.push(bail(format!(
                    "golden {g}: missing/thin justification (the protocol requires the \
                     committed-tree, two-mode evidence trail)"
                )));
            }
            goldens.push((
                g.to_string(),
                file.to_string(),
                name.to_string(),
                deps.to_string(),
            ));
        }
    }
    for (g, file, name, deps) in &goldens {
        let src = std::fs::read_to_string(root.join(file)).unwrap_or_default();
        if !src.contains(&format!("const {name}")) {
            violations.push(bail(format!(
                "golden {g}: {file} no longer declares `const {name}` — update or retire \
                 the registry row"
            )));
        }
        for pair in deps.split(',') {
            let Some((sid, pinned)) = pair.split_once('=') else {
                violations.push(bail(format!("golden {g}: malformed dependency {pair:?}")));
                continue;
            };
            let Ok(pinned) = pinned.trim().parse::<u32>() else {
                violations.push(bail(format!("golden {g}: bad pinned version {pair:?}")));
                continue;
            };
            match surface_versions.iter().find(|(id, _)| id == sid.trim()) {
                None => violations.push(bail(format!(
                    "golden {g} depends on unknown surface {sid:?} — declare the surface row"
                ))),
                Some((_, current)) if *current != pinned => violations.push(bail(format!(
                    "golden {g} was frozen against {sid}={pinned} but the surface is now \
                     {current} — re-freeze deliberately per docs/GOLDEN_POLICY.md, then \
                     update the pin"
                ))),
                Some(_) => {}
            }
        }
    }
    violations
}

// ---------------------------------------------------------------------------
// Constellation verification/fetch (beads huq.17, 1t8i): once this binary can
// run, obtain pinned sources FROM THE LOCK, verify content identity (the commit
// hash), and never silently substitute a nearby working tree that does not
// match. Cargo cannot build this in-workspace binary while a required sibling
// path is absent; the pre-Cargo clean-host entry point remains tracked by 1t8i.
// Drift, dirt, missing revisions, and identity mismatches fail with structured
// diagnostics.
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
struct LockRow {
    lib: String,
    version: String,
    git_head: String,
    remote: String,
    path: String,
}

const CONSTELLATION_LOCK_SCHEMA: &str = "frankensim-constellation-lock-v2";
const CONSTELLATION_LOCK_IDENTITY_VERSION: u32 = 1;
const CONSTELLATION_LOCK_IDENTITY_DOMAIN: &str = "org.frankensim.xtask.constellation-lock.v1";
const CONSTELLATION_LOCK_WRITER_IDENTITY_VERSION: u32 = 2;
const CONSTELLATION_LOCK_WRITER_IDENTITY_DOMAIN: &str =
    "org.frankensim.xtask.constellation-lock-writer.v2";
const CONSTELLATION_LOCK_NOTE: &str = "lock_hash covers (lib, version, git_head) only — paths are per-machine; remote is transport for bootstrap-constellation (content identity is the git head)";
const MAX_CONSTELLATION_LOCK_BYTES: usize = 1_048_576;

fn decode_constellation_lock(bytes: Vec<u8>) -> Result<String, String> {
    if bytes.len() > MAX_CONSTELLATION_LOCK_BYTES {
        return Err(format!(
            "constellation.lock exceeds the {MAX_CONSTELLATION_LOCK_BYTES}-byte parser bound"
        ));
    }
    String::from_utf8(bytes).map_err(|error| format!("constellation.lock is not UTF-8: {error}"))
}

fn read_constellation_lock(path: &Path) -> Result<String, String> {
    let file = std::fs::File::open(path)
        .map_err(|error| format!("{} unreadable: {error}", path.display()))?;
    let limit = u64::try_from(MAX_CONSTELLATION_LOCK_BYTES + 1)
        .map_err(|_| "constellation lock read bound does not fit u64".to_string())?;
    let mut bytes = Vec::new();
    file.take(limit)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("{} unreadable: {error}", path.display()))?;
    decode_constellation_lock(bytes)
}

struct CanonicalJsonParser<'a> {
    input: &'a str,
    cursor: usize,
}

impl<'a> CanonicalJsonParser<'a> {
    const fn new(input: &'a str) -> Self {
        Self { input, cursor: 0 }
    }

    fn rest(&self) -> &'a str {
        &self.input[self.cursor..]
    }

    fn expect(&mut self, literal: &str) -> Result<(), String> {
        if self.rest().starts_with(literal) {
            self.cursor += literal.len();
            Ok(())
        } else {
            Err(format!(
                "expected canonical token {literal:?} at byte {}",
                self.cursor
            ))
        }
    }

    fn consume(&mut self, literal: &str) -> bool {
        if self.rest().starts_with(literal) {
            self.cursor += literal.len();
            true
        } else {
            false
        }
    }

    fn hex_quad(&mut self) -> Result<u16, String> {
        let bytes = self
            .input
            .as_bytes()
            .get(self.cursor..self.cursor + 4)
            .ok_or_else(|| "truncated JSON Unicode escape".to_string())?;
        let mut value = 0_u16;
        for byte in bytes {
            let digit = match byte {
                b'0'..=b'9' => u16::from(byte - b'0'),
                b'a'..=b'f' => u16::from(byte - b'a' + 10),
                b'A'..=b'F' => u16::from(byte - b'A' + 10),
                _ => return Err("invalid JSON Unicode escape".to_string()),
            };
            value = (value << 4) | digit;
        }
        self.cursor += 4;
        Ok(value)
    }

    fn unicode_escape(&mut self) -> Result<char, String> {
        let high = self.hex_quad()?;
        let scalar = if (0xd800..=0xdbff).contains(&high) {
            self.expect("\\u")?;
            let low = self.hex_quad()?;
            if !(0xdc00..=0xdfff).contains(&low) {
                return Err("high surrogate is not followed by a low surrogate".to_string());
            }
            0x1_0000 + ((u32::from(high) - 0xd800) << 10) + (u32::from(low) - 0xdc00)
        } else if (0xdc00..=0xdfff).contains(&high) {
            return Err("unpaired low surrogate in JSON string".to_string());
        } else {
            u32::from(high)
        };
        char::from_u32(scalar).ok_or_else(|| "invalid JSON Unicode scalar".to_string())
    }

    fn string(&mut self) -> Result<String, String> {
        self.expect("\"")?;
        let mut value = String::new();
        loop {
            let Some(character) = self.rest().chars().next() else {
                return Err("unterminated JSON string".to_string());
            };
            self.cursor += character.len_utf8();
            match character {
                '"' => return Ok(value),
                '\\' => {
                    let Some(escaped) = self.rest().chars().next() else {
                        return Err("truncated JSON escape".to_string());
                    };
                    self.cursor += escaped.len_utf8();
                    match escaped {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        '/' => value.push('/'),
                        'b' => value.push('\u{0008}'),
                        'f' => value.push('\u{000c}'),
                        'n' => value.push('\n'),
                        'r' => value.push('\r'),
                        't' => value.push('\t'),
                        'u' => value.push(self.unicode_escape()?),
                        _ => return Err(format!("invalid JSON escape \\{escaped}")),
                    }
                }
                control if control.is_control() => {
                    return Err("unescaped control character in JSON string".to_string());
                }
                other => value.push(other),
            }
        }
    }

    fn finish(self) -> Result<(), String> {
        if self.cursor == self.input.len() {
            Ok(())
        } else {
            Err(format!("trailing JSON data at byte {}", self.cursor))
        }
    }
}

fn validate_lock_field(label: &str, value: &str) -> Result<(), String> {
    if value.is_empty() || value.chars().any(char::is_control) {
        return Err(format!(
            "constellation lock {label} must be non-empty and control-free"
        ));
    }
    Ok(())
}

fn parse_lock_rows(text: &str) -> Result<(String, Vec<LockRow>), String> {
    if text.is_empty() || text.len() > MAX_CONSTELLATION_LOCK_BYTES {
        return Err(format!(
            "constellation.lock must contain 1..={MAX_CONSTELLATION_LOCK_BYTES} UTF-8 bytes"
        ));
    }
    let mut parser = CanonicalJsonParser::new(text);
    parser.expect("{\n  \"schema\": ")?;
    let schema = parser.string()?;
    parser.expect(",\n  \"identity_domain\": ")?;
    let identity_domain = parser.string()?;
    parser.expect(&format!(
        ",\n  \"identity_version\": {CONSTELLATION_LOCK_IDENTITY_VERSION}"
    ))?;
    parser.expect(",\n  \"lock_hash\": ")?;
    let lock_hash = parser.string()?;
    parser.expect(",\n  \"note\": ")?;
    let note = parser.string()?;
    parser.expect(",\n  \"libraries\": [\n")?;

    let mut rows = Vec::new();
    loop {
        parser.expect("    {\"lib\": ")?;
        let lib = parser.string()?;
        parser.expect(", \"version\": ")?;
        let version = parser.string()?;
        parser.expect(", \"git_head\": ")?;
        let git_head = parser.string()?;
        parser.expect(", \"remote\": ")?;
        let remote = parser.string()?;
        parser.expect(", \"path\": ")?;
        let path = parser.string()?;
        parser.expect("}")?;

        for (label, value) in [
            ("library", lib.as_str()),
            ("version", version.as_str()),
            ("remote", remote.as_str()),
            ("path", path.as_str()),
        ] {
            validate_lock_field(label, value)?;
        }
        if !matches!(git_head.len(), 40 | 64)
            || !git_head
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return Err(format!("non-canonical git head for {lib}"));
        }
        rows.push(LockRow {
            lib,
            version,
            git_head,
            remote,
            path,
        });
        if parser.consume(",\n") {
            continue;
        }
        parser.expect("\n  ]\n}\n")?;
        break;
    }
    parser.finish()?;

    if schema != CONSTELLATION_LOCK_SCHEMA {
        return Err(format!("unsupported constellation lock schema {schema:?}"));
    }
    if identity_domain != CONSTELLATION_LOCK_IDENTITY_DOMAIN {
        return Err(format!(
            "unsupported constellation lock identity domain {identity_domain:?}"
        ));
    }
    if note != CONSTELLATION_LOCK_NOTE {
        return Err("constellation.lock carries a non-canonical identity note".to_string());
    }
    if lock_hash.len() != 16
        || !lock_hash
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        || rows.is_empty()
    {
        return Err("constellation.lock has no canonical hash or no libraries".to_string());
    }
    let mut expected: Vec<&str> = CONSTELLATION_REPOS.iter().map(|(lib, _)| *lib).collect();
    expected.sort_unstable();
    let mut declared: Vec<&str> = rows.iter().map(|row| row.lib.as_str()).collect();
    declared.sort_unstable();
    if declared != expected {
        return Err(format!(
            "constellation library set mismatch: declared={declared:?}, expected={expected:?}"
        ));
    }
    if render_lock_rows(&rows, &lock_hash) != text {
        return Err("constellation.lock is valid JSON but not canonical".to_string());
    }
    let identity = lock_rows_identity(&rows)?;
    let expected_hash = format!("{:016x}", fnv1a64(identity.as_bytes()));
    if lock_hash != expected_hash {
        return Err(format!(
            "constellation.lock hash {lock_hash} disagrees with declared rows {expected_hash}"
        ));
    }
    Ok((lock_hash, rows))
}

fn lock_rows_identity(rows: &[LockRow]) -> Result<String, String> {
    let mut ordered = BTreeMap::new();
    for row in rows {
        if ordered.insert(row.lib.as_str(), row).is_some() {
            return Err(format!("duplicate constellation library {:?}", row.lib));
        }
    }
    let mut identity = String::new();
    for row in ordered.values() {
        let _ = writeln!(identity, "{}={}@{}", row.lib, row.version, row.git_head);
    }
    Ok(identity)
}

fn git_out(dir: &Path, args: &[&str]) -> Result<String, String> {
    let output = sanitized_git_command(dir, args)
        .output()
        .map_err(|e| format!("git {args:?} failed to spawn: {e}"))?;
    if output.status.success() {
        let stdout = std::str::from_utf8(&output.stdout).map_err(|error| {
            format!(
                "git {args:?} in {} returned non-UTF-8 stdout: {error}",
                dir.display()
            )
        })?;
        Ok(stdout.trim().to_string())
    } else {
        let stderr = std::str::from_utf8(&output.stderr).map_err(|error| {
            format!(
                "git {args:?} in {} returned non-UTF-8 stderr: {error}",
                dir.display()
            )
        })?;
        Err(format!(
            "git {args:?} in {} failed: {}",
            dir.display(),
            stderr.trim()
        ))
    }
}

fn git_run(dir: &Path, args: &[&str]) -> Result<(), String> {
    let output = sanitized_git_command(dir, args)
        .output()
        .map_err(|error| format!("git {args:?} failed to spawn: {error}"))?;
    if output.status.success() {
        std::str::from_utf8(&output.stdout).map_err(|error| {
            format!(
                "git {args:?} in {} returned non-UTF-8 stdout: {error}",
                dir.display()
            )
        })?;
        Ok(())
    } else {
        let stderr = std::str::from_utf8(&output.stderr).map_err(|error| {
            format!(
                "git {args:?} in {} returned non-UTF-8 stderr: {error}",
                dir.display()
            )
        })?;
        Err(format!(
            "git {args:?} in {} failed: {}",
            dir.display(),
            stderr.trim()
        ))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RepositoryObservation {
    head: String,
    status: String,
}

fn repository_observation(
    target: &Path,
    expected_head: &str,
) -> Result<RepositoryObservation, String> {
    let head_before = git_out(target, &["rev-parse", "HEAD"])
        .map_err(|error| format!("{}: {error}", target.display()))?;
    let status = pinned_repository_worktree_status(target, expected_head)?;
    let head_after = git_out(target, &["rev-parse", "HEAD"])
        .map_err(|error| format!("{}: {error}", target.display()))?;
    coherent_repository_observation(target, &head_before, status, head_after)
}

fn coherent_repository_observation(
    target: &Path,
    head_before: &str,
    status: String,
    head_after: String,
) -> Result<RepositoryObservation, String> {
    if head_before != head_after {
        return Err(format!(
            "{} moved while its worktree state was being observed: before={head_before}, after={head_after}",
            target.display()
        ));
    }
    Ok(RepositoryObservation {
        head: head_after,
        status,
    })
}

const BOOTSTRAP_INCOMPLETE_KEY: &str = "frankensim.bootstrapIncomplete";

fn clear_bootstrap_marker(target: &Path) -> Result<(), String> {
    git_run(
        target,
        &["config", "--local", "--unset-all", BOOTSTRAP_INCOMPLETE_KEY],
    )
}

fn directory_is_empty(path: &Path) -> Result<bool, String> {
    let mut entries = std::fs::read_dir(path)
        .map_err(|error| format!("cannot inspect {}: {error}", path.display()))?;
    Ok(entries.next().is_none())
}

fn is_repository_root(target: &Path) -> Result<bool, String> {
    match std::fs::symlink_metadata(target.join(".git")) {
        Ok(_) => repository_worktree_status(target).map(|_| true),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(error) => Err(format!(
            "cannot inspect repository marker in {}: {error}",
            target.display()
        )),
    }
}

/// Admit a destination without deleting or repurposing existing content.
/// Returns true only when this invocation initialized the repository.
fn ensure_bootstrap_repository(target: &Path, offline: bool) -> Result<bool, String> {
    let existed = target.exists();
    if existed {
        let metadata = target
            .symlink_metadata()
            .map_err(|error| format!("cannot inspect {}: {error}", target.display()))?;
        if is_redirecting_entry(&metadata) || !metadata.is_dir() {
            return Err(format!(
                "{} exists but is not an ordinary directory",
                target.display()
            ));
        }
    }
    if !existed {
        if offline {
            return Err(format!(
                "{} missing from the source cache in --offline mode",
                target.display()
            ));
        }
        std::fs::create_dir_all(target)
            .map_err(|error| format!("cannot create {}: {error}", target.display()))?;
    }

    if is_repository_root(target)? {
        return Ok(false);
    }
    if existed && !directory_is_empty(target)? {
        return Err(format!(
            "{} is a non-empty non-git directory; refusing to repurpose it",
            target.display()
        ));
    }
    if offline {
        return Err(format!(
            "{} is not a usable cached git checkout in --offline mode",
            target.display()
        ));
    }
    git_run(target, &["init", "--quiet", "--template="])?;
    git_run(
        target,
        &["config", "--local", BOOTSTRAP_INCOMPLETE_KEY, "true"],
    )?;
    git_run(target, &["config", "--local", "core.autocrlf", "false"])?;
    Ok(true)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BootstrapOutcome {
    state: &'static str,
    transport_used: bool,
}

impl BootstrapOutcome {
    const fn pinned(marker_present: bool) -> Self {
        Self {
            state: if marker_present {
                "resumed"
            } else {
                "verified"
            },
            transport_used: false,
        }
    }

    const fn materialized(initialized: bool, offline: bool) -> Self {
        Self {
            state: if initialized { "cloned" } else { "resumed" },
            transport_used: !offline,
        }
    }
}

/// Materialize or resume one exact-pin checkout without deleting or silently
/// repurposing an ordinary existing repository. New checkouts carry a local
/// incomplete marker before fetch; an existing wrong-head checkout is resumable
/// only when that marker is present. An older unborn checkout is also safe to
/// resume when it has no worktree content and its origin matches exactly.
fn bootstrap_checkout(
    row: &LockRow,
    target: &Path,
    url: &str,
    offline: bool,
) -> Result<BootstrapOutcome, String> {
    let initialized = ensure_bootstrap_repository(target, offline)?;

    let marker_present = git_out(
        target,
        &["config", "--local", "--get", BOOTSTRAP_INCOMPLETE_KEY],
    )
    .is_ok_and(|value| value == "true");
    let current_head = git_out(target, &["rev-parse", "HEAD"]);
    if let Ok(head) = &current_head {
        if head.as_str() == row.git_head {
            verify_pinned_clean(row, target)?;
            if marker_present {
                clear_bootstrap_marker(target)?;
                return Ok(BootstrapOutcome::pinned(true));
            }
            return Ok(BootstrapOutcome::pinned(false));
        }
        if !marker_present {
            return check_pinned_clean_observation(row, target, head, "")
                .map(|()| BootstrapOutcome::pinned(false));
        }
    }

    if current_head.is_err() {
        let existing_origin = git_out(target, &["remote", "get-url", "origin"]).ok();
        if !may_resume_unborn_checkout(marker_present, existing_origin.as_deref(), url) {
            return Err(format!(
                "{} is an unmarked unborn checkout without the exact locked origin {url:?}; \
                 refusing to adopt it",
                target.display()
            ));
        }
    }

    let status = repository_worktree_status(target)?;
    if !status.is_empty() {
        return Err(format!(
            "{} is an incomplete bootstrap with worktree changes; refusing to overwrite it:\n{status}",
            target.display(),
        ));
    }
    if !offline && url == "no-remote" {
        return Err(format!(
            "lock declares no remote for {} — re-lock on a host that has one",
            row.lib
        ));
    }
    git_run(
        target,
        &["config", "--local", BOOTSTRAP_INCOMPLETE_KEY, "true"],
    )?;
    match git_out(target, &["remote", "get-url", "origin"]) {
        Ok(existing) if existing != url => {
            return Err(format!(
                "{} incomplete bootstrap origin is {existing:?}, expected {url:?}",
                target.display()
            ));
        }
        Ok(_) => {}
        Err(_) if offline => {}
        Err(_) => git_run(target, &["remote", "add", "origin", url])?,
    }
    if !offline {
        git_run(
            target,
            &[
                "fetch",
                "--no-auto-maintenance",
                "--no-recurse-submodules",
                "--quiet",
                "--depth",
                "1",
                "origin",
                &row.git_head,
            ],
        )?;
    }
    git_run(
        target,
        &[
            "checkout",
            "--no-recurse-submodules",
            "--quiet",
            "--no-overwrite-ignore",
            "--detach",
            &row.git_head,
        ],
    )?;
    verify_pinned_clean(row, target)?;
    clear_bootstrap_marker(target)?;
    Ok(BootstrapOutcome::materialized(initialized, offline))
}

fn may_resume_unborn_checkout(
    marker_present: bool,
    existing_origin: Option<&str>,
    expected_origin: &str,
) -> bool {
    marker_present || existing_origin == Some(expected_origin)
}

fn check_pinned_clean_observation(
    row: &LockRow,
    target: &Path,
    head: &str,
    status: &str,
) -> Result<(), String> {
    if head != row.git_head {
        return Err(format!(
            "{} is at {head}, lock pins {} — refusing to silently substitute a nearby \
             working tree; align or replace that sibling deliberately",
            target.display(),
            row.git_head
        ));
    }
    if !status.is_empty() {
        return Err(format!(
            "{} is DIRTY at the locked head — a modified working tree is not the pinned \
             source (a case-folding checkout collision also surfaces here); restore or \
             replace that sibling deliberately:\n{status}",
            target.display(),
        ));
    }
    Ok(())
}

fn verify_pinned_clean_with<F>(row: &LockRow, target: &Path, mut observe: F) -> Result<(), String>
where
    F: FnMut() -> Result<RepositoryObservation, String>,
{
    let first = observe()?;
    check_pinned_clean_observation(row, target, &first.head, &first.status)?;
    let confirmed = observe()?;
    check_pinned_clean_observation(row, target, &confirmed.head, &confirmed.status)?;
    if confirmed != first {
        return Err(format!(
            "{} changed between two complete pinned-state observations; refusing incoherent provenance",
            target.display()
        ));
    }
    Ok(())
}

fn verify_pinned_clean(row: &LockRow, target: &Path) -> Result<(), String> {
    verify_pinned_clean_with(row, target, || {
        repository_observation(target, &row.git_head)
    })
}

fn verify_constellation_rows(root: &Path, rows: &[LockRow]) -> Result<(), String> {
    if rows.len() != CONSTELLATION_REPOS.len() {
        return Err(format!(
            "lock declares {} libraries, expected {}",
            rows.len(),
            CONSTELLATION_REPOS.len()
        ));
    }
    let projects = root
        .parent()
        .ok_or_else(|| "workspace root has no parent".to_string())?;
    for &(expected_lib, dirname) in CONSTELLATION_REPOS {
        let row = rows
            .iter()
            .find(|row| row.lib == expected_lib)
            .ok_or_else(|| format!("lock is missing constellation library {expected_lib}"))?;
        verify_pinned_clean(row, &projects.join(dirname))?;
    }
    Ok(())
}

fn verify_live_entries_clean(entries: &[ConstellationEntry]) -> Result<(), String> {
    for entry in entries {
        let row = LockRow {
            lib: entry.lib.clone(),
            version: entry.version.clone(),
            git_head: entry.git_head.clone(),
            remote: entry.remote.clone(),
            path: entry.dir.display().to_string(),
        };
        verify_pinned_clean(&row, &entry.dir)?;
    }
    Ok(())
}

fn dirname_of(lib: &str) -> &str {
    CONSTELLATION_REPOS
        .iter()
        .find(|(l, _)| *l == lib)
        .map_or(lib, |(_, d)| d)
}

#[derive(Debug, PartialEq, Eq)]
struct BootstrapOptions {
    dest: Option<PathBuf>,
    offline: bool,
    from: Option<String>,
}

fn required_bootstrap_value<'a>(flag: &str, value: Option<&'a str>) -> Result<&'a str, String> {
    match value {
        Some(value) if !value.is_empty() && !value.starts_with('-') => Ok(value),
        _ => Err(format!("{flag} requires a non-empty value")),
    }
}

fn parse_bootstrap_options(
    default_dest: Option<PathBuf>,
    args: &[String],
) -> Result<BootstrapOptions, String> {
    let mut options = BootstrapOptions {
        dest: default_dest,
        offline: false,
        from: None,
    };
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--dest" => {
                options.dest = Some(PathBuf::from(required_bootstrap_value(
                    "--dest",
                    it.next().map(String::as_str),
                )?));
            }
            "--offline" => options.offline = true,
            "--from" => {
                options.from = Some(
                    required_bootstrap_value("--from", it.next().map(String::as_str))?.to_string(),
                );
            }
            other => return Err(format!("unknown flag {other:?}")),
        }
    }
    Ok(options)
}

fn bootstrap_provenance_row(
    row: &LockRow,
    selected_transport: &str,
    transport_used: bool,
    state: &str,
) -> BootstrapProvenanceRow {
    BootstrapProvenanceRow::new(
        &row.lib,
        &row.git_head,
        &row.remote,
        selected_transport,
        transport_used,
        state,
    )
}

fn require_unchanged_bootstrap_lock(
    lock_path: &Path,
    original_lock_text: &str,
) -> Result<(), String> {
    let observed = read_constellation_lock(lock_path)?;
    if observed == original_lock_text {
        Ok(())
    } else {
        Err(format!(
            "{} changed after bootstrap admission; refusing to publish provenance for a mixed lock epoch",
            lock_path.display()
        ))
    }
}

fn verify_provenance_publication_barrier(
    lock_path: &Path,
    original_lock_text: &str,
    dest: &Path,
    rows: &[LockRow],
) -> Result<(), String> {
    require_unchanged_bootstrap_lock(lock_path, original_lock_text)?;
    verify_two_complete_passes(rows, |row| {
        require_unchanged_bootstrap_lock(lock_path, original_lock_text)?;
        verify_pinned_clean(row, &dest.join(dirname_of(&row.lib)))?;
        require_unchanged_bootstrap_lock(lock_path, original_lock_text)
    })?;
    require_unchanged_bootstrap_lock(lock_path, original_lock_text)
}

/// `bootstrap-constellation [--offline] [--from <base>]`.
/// dest defaults to the workspace parent (where the manifests' relative
/// paths point). `--from <base>` overrides every remote with
/// `<base>/<dirname>` — the air-gapped-mirror / local-test transport.
// One protocol: fetch/verify/refuse per library + provenance; splitting
// would scatter the fail-closed invariants the acceptance audits as one.
#[allow(clippy::too_many_lines)]
fn cmd_bootstrap(root: &Path) -> ExitCode {
    let args: Vec<String> = std::env::args().skip(2).collect();
    let options = match parse_bootstrap_options(root.parent().map(Path::to_path_buf), &args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("bootstrap-constellation: {error}");
            return ExitCode::FAILURE;
        }
    };
    let Some(dest) = options.dest else {
        eprintln!("bootstrap-constellation: no destination directory");
        return ExitCode::FAILURE;
    };
    let Some(workspace_parent) = root.parent() else {
        eprintln!("bootstrap-constellation: workspace has no parent directory");
        return ExitCode::FAILURE;
    };
    if dest != workspace_parent {
        eprintln!(
            "bootstrap-constellation: --dest {} cannot satisfy this workspace's fixed sibling \
             path dependencies; use {}",
            dest.display(),
            workspace_parent.display()
        );
        return ExitCode::FAILURE;
    }
    if let Err(error) = bootstrap_provenance_support_preflight() {
        eprintln!("bootstrap-constellation: {error}");
        return ExitCode::FAILURE;
    }
    let dest_text = match provenance_path_text(&dest) {
        Ok(dest_text) => dest_text,
        Err(error) => {
            eprintln!("bootstrap-constellation: {error}");
            return ExitCode::FAILURE;
        }
    };
    let lock_path = root.join("constellation.lock");
    let lock_text = match read_constellation_lock(&lock_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("error: constellation.lock unreadable: {e} — the lock IS the input");
            return ExitCode::FAILURE;
        }
    };
    let (lock_hash, rows) = match parse_lock_rows(&lock_text) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mut provenance = Vec::new();
    let mut failures = 0usize;
    for row in &rows {
        let dirname = dirname_of(&row.lib);
        let target = dest.join(dirname);
        let url = options
            .from
            .as_ref()
            .map_or_else(|| row.remote.clone(), |base| format!("{base}/{dirname}"));
        let outcome = bootstrap_checkout(row, &target, &url, options.offline);
        match outcome {
            Ok(outcome) => {
                let state = outcome.state;
                println!(
                    "{{\"check\":\"constellation-bootstrap\",\"lib\":\"{}\",\"state\":\"{state}\",\"head\":\"{}\"}}",
                    json_escape(&row.lib),
                    json_escape(&row.git_head),
                );
                provenance.push(bootstrap_provenance_row(
                    row,
                    &url,
                    outcome.transport_used,
                    state,
                ));
            }
            Err(why) => {
                println!(
                    "{{\"check\":\"constellation-bootstrap\",\"lib\":\"{}\",\"state\":\"failed\",\"why\":\"{}\"}}",
                    json_escape(&row.lib),
                    json_escape(&why),
                );
                eprintln!("bootstrap FAILED for {}: {why}", row.lib);
                failures += 1;
            }
        }
    }
    if failures > 0 {
        eprintln!(
            "constellation bootstrap failed for {failures}/{} libraries (fail closed)",
            rows.len()
        );
        return ExitCode::FAILURE;
    }
    let prov_path = dest.join("constellation-bootstrap.json");
    if let Err(e) =
        write_bootstrap_provenance(&prov_path, &lock_hash, dest_text, &provenance, || {
            verify_provenance_publication_barrier(&lock_path, &lock_text, &dest, &rows)
        })
    {
        eprintln!("error writing bootstrap provenance: {e}");
        return ExitCode::FAILURE;
    }
    eprintln!(
        "constellation bootstrap OK: {} libraries at their locked heads under {} \
         (provenance: {})",
        rows.len(),
        dest.display(),
        prov_path.display()
    );
    ExitCode::SUCCESS
}

fn main() -> ExitCode {
    let cmd = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "check-all".to_string());
    let root = std::env::var("CARGO_WORKSPACE_DIR").map_or_else(
        |_| {
            // xtask runs from the workspace root or from xtask/; find the root by Cargo.toml.
            let cwd = std::env::current_dir().expect("cwd");
            if cwd.join("crates").is_dir() {
                cwd
            } else {
                cwd.parent().map(Path::to_path_buf).unwrap_or(cwd)
            }
        },
        PathBuf::from,
    );
    // Lockfile commands don't need the workspace manifests.
    match cmd.as_str() {
        "lock-constellation" => return cmd_constellation(&root, false),
        "check-constellation" => return cmd_constellation(&root, true),
        "bootstrap-constellation" => return cmd_bootstrap(&root),
        "generate-identities" => return identities::generate_identities(&root),
        "depgraph-receipt" => {
            let rest: Vec<String> = std::env::args().skip(2).collect();
            return match depgraph::cmd_depgraph_receipt(&root, &rest) {
                Ok(()) => ExitCode::SUCCESS,
                Err(error) => {
                    eprintln!("error: {error}");
                    ExitCode::FAILURE
                }
            };
        }
        "matdb-pack" => {
            let rest: Vec<String> = std::env::args().skip(2).collect();
            return match matdb_pack::cmd_matdb_pack(&root, &rest) {
                Ok(()) => ExitCode::SUCCESS,
                Err(error) => {
                    eprintln!("error: {error}");
                    ExitCode::FAILURE
                }
            };
        }
        _ => {}
    }
    let manifests = match load_workspace(&root) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mut policy_notes = Vec::new();
    let (violations, checks): (Vec<Violation>, Vec<&str>) = match cmd.as_str() {
        "check-layers" => (check_layers(&manifests), vec!["layers"]),
        "check-deps" => (check_deps(&manifests), vec!["dependency-policy"]),
        "check-contracts" => (check_contracts(&manifests), vec!["contracts"]),
        "check-unsafe" => (check_unsafe(&root), vec!["unsafe-capsules"]),
        "check-powi" => (check_powi(&root), vec!["powi-determinism"]),
        "check-libm" => (check_libm(&root), vec!["libm-determinism"]),
        "check-color-admission" => (check_color_admission(&root), vec!["color-admission"]),
        "check-obs-events" => (check_obs_events(&root), vec!["obs-events"]),
        "check-casual-print" => (check_casual_print(&root), vec![CASUAL_PRINT_CHECK]),
        "check-terminology" => (check_terminology(&root), vec![TERMINOLOGY_CHECK]),
        "check-goldens" => (check_goldens(&root), vec!["golden-couplings"]),
        "check-identities" => (
            identities::check_identities(&root),
            vec!["semantic-identities"],
        ),
        "check-manifest-fixture" => {
            let report = manifest_fixture::check_manifest_fixture(&root);
            policy_notes = report.decisions;
            (report.violations, vec!["manifest-fixture"])
        }
        "check-claims" => (claims::check_claims(&root), vec!["claim-state"]),
        "check-closures" => (closures::check_closures(&root), vec!["closure-evidence"]),
        "check-citable-producers" => (check_citable_producers(&root), vec![CITABLE_PRODUCER_CHECK]),
        "check-all" => {
            let mut v = check_layers(&manifests);
            v.extend(check_deps(&manifests));
            v.extend(check_contracts(&manifests));
            v.extend(check_unsafe(&root));
            v.extend(check_powi(&root));
            v.extend(check_libm(&root));
            v.extend(check_color_admission(&root));
            v.extend(check_obs_events(&root));
            v.extend(check_casual_print(&root));
            v.extend(check_terminology(&root));
            v.extend(check_goldens(&root));
            v.extend(identities::check_identities(&root));
            let manifest_report = manifest_fixture::check_manifest_fixture(&root);
            v.extend(manifest_report.violations);
            policy_notes = manifest_report.decisions;
            v.extend(claims::check_claims(&root));
            v.extend(closures::check_closures(&root));
            v.extend(check_citable_producers(&root));
            (
                v,
                vec![
                    "layers",
                    "dependency-policy",
                    "contracts",
                    "unsafe-capsules",
                    "powi-determinism",
                    "libm-determinism",
                    "color-admission",
                    "obs-events",
                    CASUAL_PRINT_CHECK,
                    TERMINOLOGY_CHECK,
                    "golden-couplings",
                    "semantic-identities",
                    "manifest-fixture",
                    "claim-state",
                    "closure-evidence",
                    CITABLE_PRODUCER_CHECK,
                ],
            )
        }
        other => {
            eprintln!(
                "unknown command {other:?}; use check-layers|check-deps|check-contracts|\
                 check-unsafe|check-powi|check-obs-events|check-casual-print|check-terminology|\
                 check-goldens|check-claims|check-closures|\
                 check-identities|check-manifest-fixture|check-citable-producers|check-all|generate-identities|\
                 lock-constellation|check-constellation|depgraph-receipt|matdb-pack"
            );
            return ExitCode::FAILURE;
        }
    };
    emit(
        &violations,
        &dev_dep_notes(&manifests),
        &policy_notes,
        &checks,
        manifests.len(),
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn casual_print_scanner_rejects_each_untyped_macro_in_core_libraries() {
        let mut sources = BTreeMap::new();
        sources.insert(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "pub fn leak() {\n",
                "    print!(\"a\"); println!(\"b\");\n",
                "    eprint!(\"c\"); eprintln!(\"d\"); dbg!(5);\n",
                "}\n",
            )
            .to_string(),
        );
        let violations = audit_casual_print_sources(&sources);
        assert_eq!(violations.len(), 5, "every casual output macro must fail");
        for name in ["print!", "println!", "eprint!", "eprintln!", "dbg!"] {
            assert!(
                violations
                    .iter()
                    .any(|violation| violation.detail.contains(name)),
                "missing mutation witness for {name}: {violations:?}"
            );
        }
    }

    #[test]
    fn casual_print_scanner_ignores_tests_strings_and_cli_boundaries() {
        let fixture = concat!(
            "pub fn typed() { let _ = \"println!(not code)\"; }\n",
            "#[cfg(test)] fn direct_test_helper() { eprintln!(\"test-only fn\"); }\n",
            "#[cfg(test)] mod diagnostics {\n",
            "    #[test] fn diagnostic() { println!(\"test-only\"); }\n",
            "}\n",
            "#[cfg(test)] impl Probe {\n",
            "    fn diagnostic(&self) { dbg!(self); }\n",
            "}\n",
        );
        let sources = [
            ("crates/fs-new/src/lib.rs", fixture),
            (
                "crates/fs-new/src/main.rs",
                "fn main() { println!(\"cli\"); }",
            ),
            (
                "crates/fs-new/src/bin/tool.rs",
                "fn main() { eprintln!(\"cli diagnostic\"); }",
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        assert!(
            audit_casual_print_sources(&sources).is_empty(),
            "test diagnostics, literal text, and CLI rendering are outside the core-library ban"
        );
    }

    #[test]
    fn casual_print_allowlist_is_exact_and_owner_scoped() {
        assert_eq!(CASUAL_PRINT_ALLOWLIST.len(), 3);
        let allowed = [
            (
                "crates/fs-casebook/src/lib.rs",
                concat!(
                    "pub fn run() {\n",
                    "::std::println!(\"{}\", record.json_line());\n",
                    "::std::println!(\"{}\", replay_record.json_line());\n",
                    "::std::println!(\"{}\", disagreement.json_line());\n",
                    "}",
                ),
            ),
            (
                "crates/fs-propcheck/src/lib.rs",
                "pub fn check_structured() { ::std::println!(\"{failure_row}\"); }",
            ),
            (
                "crates/fs-vskeleton/src/lib.rs",
                "fn emit() { ::std::eprintln!(\"{}\", e.to_jsonl()); }",
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        assert!(audit_casual_print_sources(&allowed).is_empty());

        let escaped = [(
            "crates/fs-casebook/src/lib.rs".to_string(),
            "pub fn unrelated() { println!(\"new path\"); }".to_string(),
        )]
        .into_iter()
        .collect();
        let violations = audit_casual_print_sources(&escaped);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("unrelated")),
            "an unrelated function cannot inherit the allowance: {violations:?}"
        );
    }

    #[test]
    fn casual_print_scanner_treats_comments_as_trivia_but_not_as_cfg_code() {
        let sources = [(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "const FAKE: &str = \"#[cfg(test)] mod diagnostics {\";\n",
                "// #[cfg(test)] mod diagnostics {\n",
                "pub fn leak() { println /* nested /* trivia */ remains */ ! (\"x\"); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let violations = audit_casual_print_sources(&sources);
        assert_eq!(violations.len(), 1, "comment trivia cannot hide the bang");
        assert!(violations[0].detail.contains("println!"));
    }

    #[test]
    fn casual_print_scanner_does_not_promote_macro_tokens_to_enclosing_cfg() {
        let sources = [(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "macro_rules! swallow { ($($tokens:tt)*) => {}; }\n",
                "pub fn leak() {\n",
                "    swallow!(#![cfg(test)]);\n",
                "    println!(\"production output remains visible\");\n",
                "}\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let violations = audit_casual_print_sources(&sources);
        assert_eq!(
            violations.len(),
            1,
            "macro token input cannot make the enclosing production function test-only"
        );
        assert!(violations[0].detail.contains("println!"));
    }

    #[test]
    fn casual_print_scanner_denies_cfg_and_owner_authority_to_macro_tokens() {
        let cfg_spoofs = [(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "macro_rules! strip_outer { (#[$meta:meta] $($body:tt)*) => { $($body)* }; }\n",
                "macro_rules! strip_inner { ({ #![$meta:meta] $($body:tt)* }) => {{ $($body)* }}; }\n",
                "pub fn first() { strip_outer!(#[cfg(test)] println!(\"outer production\")); }\n",
                "pub fn second() { strip_inner!({ #![cfg(test)] eprintln!(\"inner production\"); }); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let cfg_violations = audit_casual_print_sources(&cfg_spoofs);
        assert_eq!(cfg_violations.len(), 2, "{cfg_violations:?}");
        assert!(
            cfg_violations
                .iter()
                .all(|violation| violation.detail.contains("inside macro token input/body"))
        );

        let fake_owner = [(
            "crates/fs-casebook/src/lib.rs".to_string(),
            concat!(
                "struct Row; impl Row { fn json_line(&self) -> &'static str { \"{}\" } }\n",
                "macro_rules! relocate {\n",
                "(fn $claimed:ident() { $($body:tt)* }) => {\n",
                "pub fn leak() {\n",
                "let record = Row; let replay_record = Row; let disagreement = Row;\n",
                "$($body)*\n",
                "}\n",
                "};\n",
                "}\n",
                "relocate!(fn run() {\n",
                "println!(\"{}\", record.json_line());\n",
                "println!(\"{}\", replay_record.json_line());\n",
                "println!(\"{}\", disagreement.json_line());\n",
                "});\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let owner_violations = audit_casual_print_sources(&fake_owner);
        assert!(
            owner_violations
                .iter()
                .any(|violation| violation.detail.contains("expected exact invocations"))
        );
        assert_eq!(
            owner_violations
                .iter()
                .filter(|violation| violation.detail.contains("inside macro token input/body"))
                .count(),
            3,
            "macro-input `fn run` tokens cannot mint a real allowance owner: {owner_violations:?}"
        );

        let braced_header_macro = [(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "macro_rules! output_ty { () => { Result<(), ()> }; }\n",
                "pub fn genuine() -> output_ty!{} { println!(\"real body\"); Ok(()) }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let braced_header_violations = audit_casual_print_sources(&braced_header_macro);
        assert_eq!(
            braced_header_violations.len(),
            1,
            "a braced type macro is not the apparent function body: {braced_header_violations:?}"
        );
        assert!(
            braced_header_violations[0]
                .detail
                .contains("println! in fn genuine")
        );
    }

    #[test]
    fn casual_print_scanner_denies_structural_authority_to_attribute_interiors() {
        let sources = [(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "#[derive(fake::Wrap { nested: #[cfg(test)] fn forged() { println!(\"attribute output\"); } })]\n",
                "#[tool(parenthesized(#![cfg(test)] fn also_forged() {}))]\n",
                "pub fn genuine() { eprintln!(\"production output\"); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let violations = audit_casual_print_sources(&sources);
        assert_eq!(violations.len(), 2, "{violations:?}");
        assert!(violations.iter().any(|violation| {
            violation.detail.contains("attribute interior") && violation.detail.contains("println!")
        }));
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("eprintln! in fn genuine"))
        );
        assert!(
            violations
                .iter()
                .all(|violation| !violation.detail.contains("fn forged"))
        );
    }

    #[test]
    fn casual_print_scanner_keeps_unicode_lifetimes_and_identifiers_as_code() {
        let sources = [(
            "crates/fs-new/src/lib.rs".to_string(),
            "pub fn ℘<'α>() { println!(\"hidden before repair\"); } const C: char = 'x';"
                .to_string(),
        )]
        .into_iter()
        .collect();
        let violations = audit_casual_print_sources(&sources);
        assert_eq!(violations.len(), 1, "{violations:?}");
        assert!(violations[0].detail.contains("println! in fn ℘"));

        let literals = [(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "const A: &str = r#\"println!(raw)\"#;\n",
                "const B: &[u8] = b\"eprintln!(bytes)\";\n",
                "const C: &[u8] = br#\"dbg!(raw bytes)\"#;\n",
                "const D: &core::ffi::CStr = c\"print!(c string)\";\n",
                "const E: &core::ffi::CStr = cr#\"eprint!(raw c string)\"#;\n",
                "const F: char = 'λ';\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        assert!(
            audit_casual_print_sources(&literals).is_empty(),
            "protected spellings inside every Rust string/char family remain opaque"
        );
    }

    #[test]
    fn casual_print_allowlist_rejects_duplicate_owners_and_count_expansion() {
        let duplicate_owner = [(
            "crates/fs-casebook/src/lib.rs".to_string(),
            concat!(
                "fn run() {\n",
                "::std::println!(\"{}\", record.json_line());\n",
                "::std::println!(\"{}\", replay_record.json_line());\n",
                "::std::println!(\"{}\", disagreement.json_line());\n",
                "}\n",
                "mod second { fn run() {} }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let duplicate_violations = audit_casual_print_sources(&duplicate_owner);
        assert!(
            duplicate_violations
                .iter()
                .any(|violation| violation.detail.contains("2 owner body/bodies")),
            "same-name owners cannot pool their allowed count: {duplicate_violations:?}"
        );

        let expanded = [(
            "crates/fs-vskeleton/src/lib.rs".to_string(),
            concat!(
                "fn emit() {\n",
                "::std::eprintln!(\"{}\", e.to_jsonl());\n",
                "::std::eprintln!(\"{}\", extra.to_jsonl());\n",
                "}\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let expanded_violations = audit_casual_print_sources(&expanded);
        assert!(
            expanded_violations
                .iter()
                .any(|violation| violation.detail.contains("expected exact invocations")),
            "new output inside an allowed body must invalidate the ratchet: {expanded_violations:?}"
        );
    }

    #[test]
    fn casual_print_allowlist_pins_order_qualification_and_protected_bindings() {
        let reordered = [(
            "crates/fs-casebook/src/lib.rs".to_string(),
            concat!(
                "fn run() {\n",
                "::std::println!(\"{}\", replay_record.json_line());\n",
                "::std::println!(\"{}\", record.json_line());\n",
                "::std::println!(\"{}\", disagreement.json_line());\n",
                "}\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let reordered_violations = audit_casual_print_sources(&reordered);
        assert!(
            reordered_violations
                .iter()
                .any(|violation| violation.detail.contains("expected exact invocations"))
        );

        let qualified = [(
            "crates/fs-propcheck/src/lib.rs".to_string(),
            "fn check_structured() { evil::println!(\"{failure_row}\"); }".to_string(),
        )]
        .into_iter()
        .collect();
        let qualified_violations = audit_casual_print_sources(&qualified);
        assert!(
            qualified_violations
                .iter()
                .any(|violation| violation.detail.contains("evil::println"))
        );

        let rebound = [(
            "crates/fs-propcheck/src/lib.rs".to_string(),
            concat!(
                "use evil::emit as println;\n",
                "fn check_structured() { println!(\"{failure_row}\"); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let rebound_violations = audit_casual_print_sources(&rebound);
        assert!(
            rebound_violations
                .iter()
                .any(|violation| violation.detail.contains("import bound to protected name"))
        );

        let raw_declaration = [(
            "crates/fs-new/src/lib.rs".to_string(),
            "macro_rules! r#println { () => {}; }".to_string(),
        )]
        .into_iter()
        .collect();
        let raw_declaration_violations = audit_casual_print_sources(&raw_declaration);
        assert!(raw_declaration_violations.iter().any(|violation| {
            violation
                .detail
                .contains("local macro declaration uses protected name")
        }));
    }

    #[test]
    fn casual_print_allowlist_rejects_ambient_macro_import_authority() {
        let selective_macro_use = [(
            "crates/fs-propcheck/src/lib.rs".to_string(),
            concat!(
                "#[macro_use(println)] extern crate evil;\n",
                "fn check_structured() { ::std::println!(\"{failure_row}\"); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let selective_violations = audit_casual_print_sources(&selective_macro_use);
        assert!(selective_violations.iter().any(|violation| {
            violation.detail.contains("#[macro_use]")
                || violation.detail.contains("macro-authority hazard")
        }));

        let external_glob = [(
            "crates/fs-propcheck/src/lib.rs".to_string(),
            concat!(
                "use evil::*;\n",
                "fn check_structured() { ::std::println!(\"{failure_row}\"); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        let glob_violations = audit_casual_print_sources(&external_glob);
        assert!(
            glob_violations
                .iter()
                .any(|violation| violation.detail.contains("external glob import"))
        );

        let exhaustively_local_glob = [(
            "crates/fs-propcheck/src/lib.rs".to_string(),
            concat!(
                "mod local {}\n",
                "use crate::local::*;\n",
                "fn check_structured() { ::std::println!(\"{failure_row}\"); }\n",
            )
            .to_string(),
        )]
        .into_iter()
        .collect();
        assert!(
            audit_casual_print_sources(&exhaustively_local_glob).is_empty(),
            "a crate-rooted glob cannot change the explicitly pinned ::std macro path"
        );

        for crate_attribute in ["#![no_std]", "#![cfg_attr(all(), no_std)]"] {
            let rebound_std = [(
                "crates/fs-propcheck/src/lib.rs".to_string(),
                format!(
                    "{crate_attribute}\nextern crate evil as std;\nfn check_structured() {{ ::std::println!(\"{{failure_row}}\"); }}\n"
                ),
            )]
            .into_iter()
            .collect();
            let rebound_std_violations = audit_casual_print_sources(&rebound_std);
            assert!(rebound_std_violations.iter().any(|violation| {
                violation.detail.contains("compiler-provided std")
                    || violation.detail.contains("macro-authority hazard")
            }));
        }
    }

    #[test]
    fn casual_print_macro_authority_hazards_are_strictly_capped() {
        let mut source = String::new();
        for index in 0..=CASUAL_MAX_DIAGNOSTICS {
            let _ = writeln!(source, "#[macro_use] mod imported_{index} {{}}");
        }
        source.push_str("fn check_structured() { ::std::println!(\"{failure_row}\"); }");
        let sources = [("crates/fs-propcheck/src/lib.rs".to_string(), source)]
            .into_iter()
            .collect();
        let violations = audit_casual_print_sources(&sources);
        assert_eq!(violations.len(), 1, "{violations:?}");
        assert!(violations[0].detail.contains("hazard cap"));
    }

    #[test]
    fn casual_print_scanner_rejects_import_aliases_and_malformed_input() {
        let aliased = [(
            "crates/fs-new/src/lib.rs".to_string(),
            "pub fn leak() { use std::println as alias; alias!(\"hidden\"); }".to_string(),
        )]
        .into_iter()
        .collect();
        let alias_violations = audit_casual_print_sources(&aliased);
        assert!(
            alias_violations.iter().any(|violation| violation
                .detail
                .contains("protected name renamed by import")),
            "protected macro aliases must fail closed: {alias_violations:?}"
        );

        let malformed = [(
            "crates/fs-new/src/lib.rs".to_string(),
            "pub fn leak() { println /* never closed".to_string(),
        )]
        .into_iter()
        .collect();
        let malformed_violations = audit_casual_print_sources(&malformed);
        assert_eq!(malformed_violations.len(), 1);
        assert!(
            malformed_violations[0]
                .detail
                .contains("cannot exhaustively scan")
        );
    }

    #[test]
    fn casual_print_inventory_scans_nested_target_dirs_and_graph_ignores_decoys() {
        let root = std::env::temp_dir().join(format!(
            "xtask-casual-print-inventory-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(root.join("crates/mini/src/target")).expect("source target dir");
        std::fs::create_dir_all(root.join("crates/mini/target")).expect("nested target dir");
        std::fs::create_dir_all(root.join("target")).expect("workspace build root");
        std::fs::create_dir_all(root.join("tools")).expect("tools root");
        std::fs::create_dir_all(root.join("xtask")).expect("xtask root");
        std::fs::write(
            root.join("crates/mini/Cargo.toml"),
            "[package]\nname='mini'\n",
        )
        .expect("mini manifest");
        std::fs::write(
            root.join("crates/mini/src/target/mod.rs"),
            "pub fn leak() { println!(\"must be inventoried\"); }\n",
        )
        .expect("source target module");
        std::fs::write(root.join("crates/mini/src/lib.rs"), "mod target;\n").expect("library root");
        std::fs::write(
            root.join("crates/mini/target/generated.rs"),
            "pub fn decoy() { println!(\"unreachable nested target decoy\"); }\n",
        )
        .expect("nested target fixture");
        std::fs::write(
            root.join("target/generated.rs"),
            "compile_error!(\"workspace build output is outside owned roots\");\n",
        )
        .expect("workspace build-root fixture");

        let sources = workspace_rust_sources(&root).expect("workspace inventory");
        assert!(sources.contains_key("crates/mini/src/target/mod.rs"));
        assert!(sources.contains_key("crates/mini/target/generated.rs"));
        assert!(!sources.contains_key("target/generated.rs"));
        assert_eq!(audit_casual_print_sources(&sources).len(), 1);
    }

    #[test]
    fn casual_print_inventory_enforces_bytes_before_unbounded_reading() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "xtask-casual-print-byte-cap-{}-{nonce}.rs",
            std::process::id()
        ));
        std::fs::write(&path, "123456789").expect("bounded-read fixture");
        let error = casual_read_bounded_utf8(&path, 8, "fixture source")
            .expect_err("metadata must refuse the oversized file before reading it");
        assert!(error.contains("exceeds the remaining 8-byte"), "{error}");
    }

    #[cfg(unix)]
    #[test]
    fn casual_print_inventory_refuses_lossy_backslash_path_identity() {
        let root = Path::new("/repo");
        let error =
            casual_inventory_relative(root, Path::new("/repo/crates/fs-new/src/decoy\\name.rs"))
                .expect_err("a Unix filename backslash must not collapse into a separator");
        assert!(error.contains("backslash-bearing"), "{error}");
        assert!(
            casual_normalized_inclusion("crates/fs-new/src/lib.rs", "decoy\\name.rs").is_err(),
            "module literals use slash-only portable identity"
        );
    }

    #[cfg(unix)]
    #[test]
    fn casual_print_inventory_refuses_symlinks_under_owned_roots() {
        use std::os::unix::fs::symlink;

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "xtask-casual-print-symlink-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir_all(root.join("crates/mini/src")).expect("crate source root");
        std::fs::create_dir_all(root.join("tools")).expect("tools root");
        std::fs::create_dir_all(root.join("xtask")).expect("xtask root");
        std::fs::write(root.join("crates/mini/src/lib.rs"), "").expect("real source");
        symlink("lib.rs", root.join("crates/mini/src/indirect.rs")).expect("source symlink");
        let error = workspace_rust_sources(&root).expect_err("symlink must refuse inventory");
        assert!(error.contains("refusing symlink"), "{error}");
    }

    #[test]
    fn casual_print_scanner_reclassifies_path_included_cli_sources() {
        let included = [
            (
                "crates/fs-new/src/lib.rs",
                concat!(
                    "#[path = \"main.rs\"] mod embedded_main;\n",
                    "#[path = r#\"bin/tool.rs\"#] mod embedded_bin;\n",
                ),
            ),
            (
                "crates/fs-new/src/main.rs",
                "pub fn leak() { println!(\"library-reachable main\"); }",
            ),
            (
                "crates/fs-new/src/bin/tool.rs",
                "pub fn leak() { eprintln!(\"library-reachable bin\"); }",
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        let included_violations = audit_casual_print_sources(&included);
        assert_eq!(included_violations.len(), 2, "{included_violations:?}");

        let test_only = [
            (
                "crates/fs-new/src/lib.rs",
                concat!(
                    "#[cfg(test)] #[path = \"main.rs\"] mod diagnostics_a;\n",
                    "#[path = \"main.rs\"] #[cfg(test)] mod diagnostics_b;\n",
                ),
            ),
            (
                "crates/fs-new/src/main.rs",
                "fn main() { println!(\"conventional CLI remains exempt\"); }",
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        assert!(audit_casual_print_sources(&test_only).is_empty());
    }

    #[test]
    fn casual_print_scanner_follows_the_complete_library_module_graph() {
        let ordinary = [
            ("crates/fs-new/src/lib.rs", "mod first;\n"),
            ("crates/fs-new/src/first.rs", "mod second;\n"),
            (
                "crates/fs-new/src/first/second.rs",
                "pub fn leak() { println!(\"transitively reachable\"); }\n",
            ),
            (
                "crates/fs-new/src/decoy.rs",
                "pub fn decoy() { eprintln!(\"unreachable decoy\"); }\n",
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        let ordinary_violations = audit_casual_print_sources(&ordinary);
        assert_eq!(ordinary_violations.len(), 1, "{ordinary_violations:?}");
        assert!(ordinary_violations[0].detail.contains("first/second.rs"));

        let direct_and_transitive = [
            (
                "crates/fs-new/src/lib.rs",
                "#[path = \"alternate/first.rs\"] mod first;\n",
            ),
            (
                "crates/fs-new/src/alternate/first.rs",
                "#[path = \"../real/second.rs\"] mod second;\n",
            ),
            (
                "crates/fs-new/src/real/second.rs",
                "pub fn leak() { dbg!(\"direct transitive path\"); }\n",
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        let path_violations = audit_casual_print_sources(&direct_and_transitive);
        assert_eq!(path_violations.len(), 1, "{path_violations:?}");
        assert!(path_violations[0].detail.contains("real/second.rs"));
    }

    #[test]
    fn casual_print_scanner_walks_deep_graphs_without_ancestry_cloning() {
        const DEPTH: usize = 512;
        let mut sources = BTreeMap::new();
        sources.insert(
            "crates/fs-new/src/lib.rs".to_string(),
            "#[path = \"node-0.rs\"] mod node_0;".to_string(),
        );
        for index in 0..DEPTH {
            let source = if index + 1 == DEPTH {
                String::new()
            } else {
                format!(
                    "#[path = \"node-{}.rs\"] mod node_{};",
                    index + 1,
                    index + 1
                )
            };
            sources.insert(format!("crates/fs-new/src/node-{index}.rs"), source);
        }
        assert!(
            audit_casual_print_sources(&sources).is_empty(),
            "the iterative active/done walk must admit a deep acyclic graph"
        );
    }

    #[test]
    fn casual_print_cfg_span_sweep_coalesces_repeated_test_attributes() {
        let mut source = "#[cfg(test)]\n".repeat(512);
        source.push_str("fn diagnostics() { println!(\"test only\"); }");
        let sources = [("crates/fs-new/src/lib.rs".to_string(), source)]
            .into_iter()
            .collect();
        assert!(audit_casual_print_sources(&sources).is_empty());
    }

    #[test]
    fn casual_print_scanner_refuses_ambiguous_missing_cyclic_escaped_or_aliased_modules() {
        let cases: [(&str, Vec<(&str, &str)>); 4] = [
            (
                "ambiguous module",
                vec![
                    ("crates/fs-new/src/lib.rs", "mod a;"),
                    ("crates/fs-new/src/a.rs", ""),
                    ("crates/fs-new/src/a/mod.rs", ""),
                ],
            ),
            (
                "missing module",
                vec![("crates/fs-new/src/lib.rs", "mod absent;")],
            ),
            (
                "cyclic module inclusion",
                vec![
                    ("crates/fs-new/src/lib.rs", "#[path = \"a.rs\"] mod a;"),
                    ("crates/fs-new/src/a.rs", "#[path = \"lib.rs\"] mod root;"),
                ],
            ),
            (
                "outside package root",
                vec![
                    (
                        "crates/fs-new/src/lib.rs",
                        "#[path = \"../../../outside.rs\"] mod escaped;",
                    ),
                    ("outside.rs", ""),
                ],
            ),
        ];
        for (expected, rows) in cases {
            let sources = rows
                .into_iter()
                .map(|(path, source)| (path.to_string(), source.to_string()))
                .collect();
            let violations = audit_casual_print_sources(&sources);
            assert_eq!(violations.len(), 1, "{expected}: {violations:?}");
            assert!(
                violations[0].detail.contains(expected),
                "expected {expected:?}: {violations:?}"
            );
        }

        let aliased = [
            (
                "crates/fs-new/src/lib.rs",
                "#[path = \"shared.rs\"] mod a; #[path = \"shared.rs\"] mod b;",
            ),
            ("crates/fs-new/src/shared.rs", ""),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        let alias_violations = audit_casual_print_sources(&aliased);
        assert!(alias_violations[0].detail.contains("aliased module target"));
    }

    #[test]
    fn casual_print_scanner_refuses_unsupported_conditional_or_inline_path_contexts() {
        let conditional = [
            (
                "crates/fs-new/src/lib.rs",
                "#[cfg_attr(unix, path = \"unix.rs\")] mod platform;",
            ),
            ("crates/fs-new/src/unix.rs", ""),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        let conditional_violations = audit_casual_print_sources(&conditional);
        assert!(
            conditional_violations[0]
                .detail
                .contains("cfg_attr to select a module path")
        );

        let inline = [
            (
                "crates/fs-new/src/lib.rs",
                "mod inline { #[path = \"nested.rs\"] mod nested; }",
            ),
            ("crates/fs-new/src/nested.rs", ""),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect();
        let inline_violations = audit_casual_print_sources(&inline);
        assert!(
            inline_violations[0]
                .detail
                .contains("inline-module context")
        );
    }

    #[test]
    fn casual_print_scanner_seeds_explicit_manifest_library_paths() {
        let sources = [(
            "crates/fs-new/src/main.rs".to_string(),
            "fn library_entry() { println!(\"not a CLI target\"); }".to_string(),
        )]
        .into_iter()
        .collect();
        let root = casual_manifest_library_root(
            "crates/fs-new",
            "[package]\nname = \"fs-new\"\nautolib = false\n[lib]\npath = \"src/main.rs\"\n",
            &sources,
        )
        .expect("manifest must parse")
        .expect("explicit [lib] target");
        let violations =
            audit_casual_print_sources_with_roots(&sources, &[(root, "crates/fs-new".to_string())]);
        assert_eq!(violations.len(), 1, "{violations:?}");
        assert!(violations[0].detail.contains("src/main.rs"));
    }

    #[test]
    fn casual_print_manifest_normalizes_quoted_and_dotted_target_keys() {
        let sources = [(
            "crates/fs-new/src/hidden.rs".to_string(),
            "fn library_entry() { println!(\"not a CLI target\"); }".to_string(),
        )]
        .into_iter()
        .collect();
        for manifest in [
            "['package']\nname = 'fs-new'\n['lib']\n'path' = 'src/hidden.rs'\n",
            "package.name = 'fs-new'\nlib.path = 'src/hidden.rs'\n",
        ] {
            let root = casual_manifest_library_root("crates/fs-new", manifest, &sources)
                .expect("Cargo-equivalent key spelling must parse")
                .expect("explicit library root");
            assert_eq!(root, "crates/fs-new/src/hidden.rs");
            let violations = audit_casual_print_sources_with_roots(
                &sources,
                &[(root, "crates/fs-new".to_string())],
            );
            assert_eq!(violations.len(), 1, "{manifest:?}: {violations:?}");
        }

        let inline = casual_manifest_library_root(
            "crates/fs-new",
            "package = { name = 'fs-new' }\nlib = { path = 'src/hidden.rs' }\n",
            &sources,
        )
        .expect_err("unsupported inline authority tables must refuse, not disappear");
        assert!(inline.contains("inline package/lib table"), "{inline}");
    }

    #[test]
    fn goldens_seeded_upstream_drift_points_at_the_coupling_row() {
        let root = std::env::temp_dir().join(format!("xtask-goldens-{}", std::process::id()));
        let src_dir = root.join("crates/mini/src");
        std::fs::create_dir_all(&src_dir).expect("fixture dirs");
        let write = |rel: &str, text: &str| {
            std::fs::write(root.join(rel), text).expect("fixture write");
        };
        write(
            "crates/mini/src/lib.rs",
            "pub const MINI_SEMANTICS_VERSION: u32 = 1;\nconst GOLDEN_HASH: u64 = 7;\n",
        );
        write(
            "golden-couplings.json",
            "{\n\"surfaces\": [\n{\"id\":\"mini:semantics\",\"file\":\"crates/mini/src/lib.rs\",\"const\":\"MINI_SEMANTICS_VERSION\",\"version\":1}\n],\n\"goldens\": [\n{\"golden\":\"mini:golden\",\"file\":\"crates/mini/src/lib.rs\",\"const\":\"GOLDEN_HASH\",\"depends_on\":\"mini:semantics=1\",\"justification\":\"recorded at fixture landing, both modes, committed tree\"}\n]\n}\n",
        );
        assert!(
            check_goldens(&root).is_empty(),
            "clean fixture must pass: {:?}",
            check_goldens(&root)
        );
        // The SEEDED UPSTREAM CHANGE: bump the source semantics const
        // without re-freezing — the checker must fail and point at the
        // dependent golden row.
        write(
            "crates/mini/src/lib.rs",
            "pub const MINI_SEMANTICS_VERSION: u32 = 2;\nconst GOLDEN_HASH: u64 = 7;\n",
        );
        let v = check_goldens(&root);
        assert_eq!(v.len(), 1, "{v:?}");
        assert!(
            v[0].detail.contains("mini:semantics")
                && v[0].detail.contains("mini:golden")
                && v[0].detail.contains("re-freeze"),
            "drift names the surface AND the dependent golden: {}",
            v[0].detail
        );
        // A stale golden pin (registry surface moved on without the
        // golden) also refuses with the protocol pointer.
        write(
            "crates/mini/src/lib.rs",
            "pub const MINI_SEMANTICS_VERSION: u32 = 1;\nconst GOLDEN_HASH: u64 = 7;\n",
        );
        write(
            "golden-couplings.json",
            "{\n\"surfaces\": [\n{\"id\": \"mini:semantics\", \"file\": \"crates/mini/src/lib.rs\", \"const\": \"MINI_SEMANTICS_VERSION\", \"version\": 1}\n],\n\"goldens\": [\n{\"golden\": \"mini:golden\", \"file\": \"crates/mini/src/lib.rs\", \"const\": \"GOLDEN_HASH\", \"depends_on\": \"mini:semantics=0\", \"justification\": \"recorded at fixture landing, both modes, committed tree\"}\n]\n}\n",
        );
        let v = check_goldens(&root);
        assert_eq!(v.len(), 1, "{v:?}");
        assert!(v[0].detail.contains("GOLDEN_POLICY"), "{}", v[0].detail);
    }

    #[test]
    fn goldens_resolve_qualified_domain_constants_fail_closed() {
        let root = std::env::temp_dir().join(format!(
            "xtask-goldens-qualified-domain-{}",
            std::process::id()
        ));
        let write = |rel: &str, text: &str| {
            let path = root.join(rel);
            std::fs::create_dir_all(path.parent().expect("fixture parent")).expect("fixture dirs");
            std::fs::write(path, text).expect("fixture write");
        };
        write(
            "crates/owner/src/lib.rs",
            "pub const OWNER_VERSION: u32 = 2;\nconst GOLDEN_HASH: u64 = 7;\n",
        );
        write(
            "crates/shared/src/lib.rs",
            "pub const SHARED_DOMAIN: &str = \"fs-recompute-node-v2\";\n",
        );
        let registry = "{\n\"surfaces\": [\n{\"id\": \"owner:identity\", \"file\": \"crates/owner/src/lib.rs\", \"const\": \"OWNER_VERSION\", \"version\": 2, \"domain_const\": \"crates/shared/src/lib.rs#SHARED_DOMAIN\", \"domain\": \"fs-recompute-node-v2\"}\n],\n\"goldens\": [\n{\"golden\": \"owner:golden\", \"file\": \"crates/owner/src/lib.rs\", \"const\": \"GOLDEN_HASH\", \"depends_on\": \"owner:identity=2\", \"justification\": \"recorded at fixture landing, both modes, committed tree\"}\n]\n}\n";
        write("golden-couplings.json", registry);
        assert!(
            check_goldens(&root).is_empty(),
            "qualified domain fixture must pass: {:?}",
            check_goldens(&root)
        );

        write(
            "crates/shared/src/lib.rs",
            "pub const SHARED_DOMAIN: &str = \"fs-recompute-node-v3\";\n",
        );
        let violations = check_goldens(&root);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("domain drifted")
                    && violation.detail.contains("crates/shared/src/lib.rs")),
            "qualified domain drift must name the external source: {violations:?}"
        );

        write(
            "golden-couplings.json",
            &registry.replace(
                "crates/shared/src/lib.rs#SHARED_DOMAIN",
                "../outside.rs#SHARED_DOMAIN",
            ),
        );
        let violations = check_goldens(&root);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("not a safe repo-relative")),
            "unsafe qualified domain reference must fail closed: {violations:?}"
        );
    }

    use super::*;

    fn audited_citable_source_fixture() -> BTreeMap<String, String> {
        [
            (
                "crates/fs-roofline/src/production.rs",
                include_str!("../../crates/fs-roofline/src/production.rs"),
            ),
            (
                "crates/fs-roofline/src/lib.rs",
                include_str!("../../crates/fs-roofline/src/lib.rs"),
            ),
            (
                "crates/fs-roofline/src/bin/roofline.rs",
                include_str!("../../crates/fs-roofline/src/bin/roofline.rs"),
            ),
            (
                "crates/fs-feec/tests/perf_lane.rs",
                include_str!("../../crates/fs-feec/tests/perf_lane.rs"),
            ),
            (
                "crates/fs-fft/tests/perf_lane.rs",
                include_str!("../../crates/fs-fft/tests/perf_lane.rs"),
            ),
        ]
        .into_iter()
        .map(|(path, source)| (path.to_string(), source.to_string()))
        .collect()
    }

    #[test]
    fn citable_producer_inventory_accepts_exact_five_audited_sinks() {
        assert_eq!(CITABLE_PRODUCER_ALLOWLIST.len(), 5);
        let violations = audit_citable_producer_sources(&audited_citable_source_fixture());
        assert!(
            violations.is_empty(),
            "the five audited authority-gated sinks must be exact: {violations:?}"
        );
    }

    #[test]
    fn citable_producer_inventory_is_clean_on_the_real_workspace_tree() {
        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let root = manifest_dir
            .parent()
            .expect("xtask lives under workspace root");
        let violations = check_citable_producers(root);
        assert!(
            violations.is_empty(),
            "the command must accept the actual workspace inventory: {violations:?}"
        );
    }

    #[test]
    fn external_gate_recorder_parser_diagnostics_and_test_fixtures_are_not_producers() {
        let source = include_str!("../../crates/fs-roofline/src/lib.rs");
        let occurrences = scan_citable_producer_source("crates/fs-roofline/src/lib.rs", source)
            .expect("recorder source is valid Rust");
        assert!(
            occurrences.is_empty(),
            "parser keys, diagnostics, and cfg(test) gate fixtures are consumers/oracles: {occurrences:?}"
        );
    }

    #[test]
    fn citable_producer_inventory_rejects_an_unknown_sixth_producer() {
        let mut sources = audited_citable_source_fixture();
        sources.insert(
            "crates/fs-new/src/lib.rs".to_string(),
            concat!(
                "\nfn serialize_sixth(citation_",
                "eligible: bool) -> String {\n    format!(\"{{\\\"schema\\\":\\\"new-producer-v1\\\",\\\"citation_",
                "eligible\\\":{citation_",
                "eligible}}}\")\n}\n",
            )
            .to_string(),
        );
        let violations = audit_citable_producer_sources(&sources);
        let sixth: Vec<_> = violations
            .iter()
            .filter(|violation| violation.crate_name == "crates/fs-new/src/lib.rs")
            .collect();
        assert_eq!(
            sixth.len(),
            1,
            "unknown sixth sink must fail once: {violations:?}"
        );
        assert!(
            sixth[0].detail.contains("positive-capable")
                && sixth[0].detail.contains("serialize_sixth"),
            "failure must name the positive producer and its owner: {}",
            sixth[0].detail
        );
    }

    #[test]
    fn citable_producer_scanner_separates_fixed_false_output_from_matchers() {
        let source = concat!(
            "\nfn report_only() {\n    println!(\"{{\\\"citation_",
            "eligible\\\":false,\\\"reason\\\":\\\"candidate\\\"}}\")\n}\n\n",
            "fn matcher_only(payload: &str) {\n    assert!(payload.contains(\"\\\"citation_",
            "eligible\\\":true\"));\n}\n",
        );
        let occurrences = scan_citable_producer_source("crates/fs-new/src/lib.rs", source)
            .expect("valid Rust-shaped fixture");
        assert_eq!(occurrences.len(), 1, "matcher strings are not producers");
        assert_eq!(occurrences[0].owner, "report_only");
        assert_eq!(occurrences[0].mode, CitationFieldMode::ReportOnlyFalse);
    }

    #[test]
    fn citable_producer_scanner_rejects_split_and_variable_backed_fields() {
        let source = concat!(
            "fn split_concat(value: bool) -> String {\n",
            "    let key = concat!(\"citation_\", \"eligible\");\n",
            "    format!(\"{{\\\"{key}\\\":{value}}}\")\n",
            "}\n\n",
            "fn variable_backed(value: bool) -> String {\n",
            "    let suffix = \"eligible\";\n",
            "    let prefix = \"citation_\";\n",
            "    let key = format!(\"{prefix}{suffix}\");\n",
            "    format!(\"{{\\\"{key}\\\":{value}}}\")\n",
            "}\n",
        );
        let occurrences = scan_citable_producer_source("crates/fs-new/src/lib.rs", source)
            .expect("valid Rust-shaped fixture");
        assert_eq!(
            occurrences.len(),
            2,
            "both dynamic protected-key constructions must enter the inventory"
        );
        assert_eq!(occurrences[0].owner, "split_concat");
        assert_eq!(occurrences[1].owner, "variable_backed");
        assert!(
            occurrences
                .iter()
                .all(|occurrence| occurrence.mode == CitationFieldMode::PositiveCapable),
            "dynamic assembly can never claim the fixed-false exemption"
        );
    }

    #[test]
    fn citable_producer_fragments_survive_an_intervening_direct_sink() {
        let source = concat!(
            "fn mixed(value: bool) -> String {\n",
            "    let prefix = \"citation_\";\n",
            "    let report = \"{\\\"citation_eligible\\\":false}\";\n",
            "    let suffix = \"eligible\";\n",
            "    let key = format!(\"{prefix}{suffix}\");\n",
            "    format!(\"{report} {{\\\"{key}\\\":{value}}}\")\n",
            "}\n",
        );
        let occurrences = scan_citable_producer_source("crates/fs-new/src/lib.rs", source)
            .expect("valid Rust-shaped fixture");
        assert_eq!(
            occurrences.len(),
            2,
            "the audited-looking direct field must not erase a hidden split producer"
        );
        assert_eq!(occurrences[0].mode, CitationFieldMode::ReportOnlyFalse);
        assert_eq!(occurrences[1].mode, CitationFieldMode::PositiveCapable);
    }

    #[test]
    fn citable_producer_scanner_rejects_one_literal_runtime_reassembly() {
        let source = concat!(
            "fn reassembled(value: bool) -> String {\n",
            "    let (prefix, suffix) = \"citation_|eligible\".split_once('|').unwrap();\n",
            "    let key = format!(\"{prefix}{suffix}\");\n",
            "    format!(\"{{\\\"{key}\\\":{value}}}\")\n",
            "}\n",
            "fn compound(value: bool) -> String {\n",
            "    let (field, suffix) = \"citation_eligible|eligible\".split_once('|').unwrap();\n",
            "    let key = format!(\"{}{suffix}\", &field[..9]);\n",
            "    format!(\"{{\\\"{key}\\\":{value}}}\")\n",
            "}\n",
            "fn direct_compound(value: bool) -> String {\n",
            "    let (field, suffix) = \"{\\\"citation_eligible\\\":false}|eligible\".split_once('|').unwrap();\n",
            "    let key = format!(\"{}{suffix}\", &field[2..11]);\n",
            "    format!(\"{field} {{\\\"{key}\\\":{value}}}\")\n",
            "}\n",
        );
        let occurrences = scan_citable_producer_source("crates/fs-new/src/lib.rs", source)
            .expect("valid Rust-shaped fixture");
        assert_eq!(occurrences.len(), 4, "one-literal split producers");
        assert_eq!(occurrences[0].owner, "reassembled");
        assert_eq!(occurrences[1].owner, "compound");
        assert_eq!(occurrences[2].owner, "direct_compound");
        assert_eq!(occurrences[3].owner, "direct_compound");
        assert_eq!(occurrences[2].mode, CitationFieldMode::ReportOnlyFalse);
        assert!(
            occurrences[0..2]
                .iter()
                .chain(&occurrences[3..4])
                .all(|occurrence| occurrence.mode == CitationFieldMode::PositiveCapable)
        );
    }

    #[test]
    fn citable_producer_complete_format_capture_remains_one_direct_sink() {
        let source = concat!(
            "fn admitted(citation_eligible: bool) -> String {\n",
            "    format!(\"{{\\\"citation_eligible\\\":{citation_eligible}}}\")\n",
            "}\n",
        );
        let occurrences = scan_citable_producer_source("crates/fs-new/src/lib.rs", source)
            .expect("valid Rust-shaped fixture");
        assert_eq!(
            occurrences.len(),
            1,
            "the value capture repeats vocabulary but emits no second key"
        );
        assert_eq!(occurrences[0].owner, "admitted");
        assert_eq!(occurrences[0].mode, CitationFieldMode::PositiveCapable);
    }

    #[test]
    fn citable_producer_fragments_do_not_join_across_functions() {
        let source = concat!(
            "fn prefix_only() { let _ = \"citation_\"; }\n",
            "const MODULE_SUFFIX: &str = \"eligible\";\n",
            "fn suffix_only() { let _ = \"eligible\"; }\n",
        );
        let occurrences = scan_citable_producer_source("crates/fs-new/src/lib.rs", source)
            .expect("valid Rust-shaped fixture");
        assert!(
            occurrences.is_empty(),
            "closed functions and module constants cannot assemble one runtime field: {occurrences:?}"
        );
    }

    #[test]
    fn citable_producer_guard_does_not_inventory_its_own_literals() {
        let source = include_str!("main.rs");
        let occurrences = scan_citable_producer_source("xtask/src/main.rs", source)
            .expect("the guard source is valid Rust");
        assert!(occurrences.is_empty(), "guard self-match: {occurrences:?}");
    }

    #[test]
    fn citable_producer_guard_metadata_exemption_never_hides_executable_code() {
        let source = concat!(
            "// const CITATION_FIELD_NAME: &str = \"spoof\";\n",
            "static COMMENT_HIDDEN: &str = \"{\\\"citation_eligible\\\":true}\";\n",
            "const CITATION_FIELD_NAME: &str = concat!(\"citation_\", \"eligible\");\n",
            "static HIDDEN: &str = \"{\\\"citation_eligible\\\":true}\";\n",
            "const CITABLE_PRODUCER_ALLOWLIST: &[()] = &[];\n",
            "static AFTER: [&str; 1] = [\"{\\\"citation_eligible\\\":true}\"];\n",
            "fn hidden(value: bool) -> String {\n",
            "    let prefix = \"citation_\";\n",
            "    let suffix = \"eligible\";\n",
            "    let key = format!(\"{prefix}{suffix}\");\n",
            "    format!(\"{{\\\"{key}\\\":{value}}}\")\n",
            "}\n",
            "struct RustStringLiteral;\n",
        );
        let occurrences = scan_citable_producer_source("xtask/src/main.rs", source)
            .expect("valid Rust-shaped guard fixture");
        assert_eq!(
            occurrences.len(),
            4,
            "module and executable guard-region producers"
        );
        assert_eq!(occurrences[0].owner, "<module>");
        assert_eq!(occurrences[1].owner, "<module>");
        assert_eq!(occurrences[2].owner, "<module>");
        assert_eq!(occurrences[3].owner, "hidden");
        assert!(
            occurrences
                .iter()
                .all(|occurrence| occurrence.mode == CitationFieldMode::PositiveCapable)
        );
    }

    fn manifest(name: &str, layer: &str, deps: &[&str]) -> Manifest {
        let mut toml = format!(
            "[package]\nname = \"{name}\"\n[package.metadata.frankensim]\nlayer = \"{layer}\"\n[dependencies]\n"
        );
        for d in deps {
            let _ = writeln!(toml, "{d} = {{ path = \"../{d}\" }}");
        }
        parse_manifest(Path::new(&format!("crates/{name}/Cargo.toml")), &toml)
            .expect("fixture parses")
    }

    #[test]
    fn parses_generated_manifest_shape() {
        let m = manifest("fs-exec", "L0", &["fs-qty", "asupersync"]);
        assert_eq!(m.name, "fs-exec");
        assert_eq!(m.layer, Layer::L0);
        assert_eq!(
            m.runtime_deps,
            vec!["fs-qty".to_string(), "asupersync".to_string()]
        );
    }

    #[test]
    fn manifest_parser_skips_multiline_feature_arrays() {
        let toml = concat!(
            "[package]\n",
            "name = \"fs-x\"\n",
            "[package.metadata.frankensim]\n",
            "layer = \"L1\"\n",
            "[dependencies]\n",
            "fs-qty = { path = \"../fs-qty\", optional = true }\n",
            "[features]\n",
            "frontier = [\n",
            "  \"dep:fs-qty\",\n",
            "]\n",
        );
        let parsed = parse_manifest(Path::new("crates/fs-x/Cargo.toml"), toml)
            .expect("multiline feature arrays are valid TOML");
        assert_eq!(parsed.runtime_deps, vec!["fs-qty"]);
    }

    #[test]
    fn missing_layer_declaration_is_an_error() {
        let toml = "[package]\nname = \"fs-x\"\n[dependencies]\n";
        let err = parse_manifest(Path::new("crates/fs-x/Cargo.toml"), toml).unwrap_err();
        assert!(
            err.contains("missing [package.metadata.frankensim] layer"),
            "got: {err}"
        );
    }

    #[test]
    fn layer_direction_allows_downward_and_util() {
        let ms = vec![
            manifest("fs-qty", "UTIL", &[]),
            manifest("fs-la", "L1", &["fs-qty", "fs-substrate"]),
            manifest("fs-substrate", "L0", &["fs-qty"]),
        ];
        assert!(check_layers(&ms).is_empty());
    }

    #[test]
    fn seeded_layer_violation_is_caught() {
        // L1 depending on L2 must fail: lower layers never know about higher ones.
        let ms = vec![
            manifest("fs-geom", "L2", &[]),
            manifest("fs-la", "L1", &["fs-geom"]),
        ];
        let v = check_layers(&ms);
        assert_eq!(v.len(), 1, "expected exactly one violation, got {v:?}");
        assert_eq!(v[0].crate_name, "fs-la");
        assert!(v[0].detail.contains("must not depend on fs-geom"));
    }

    #[test]
    fn ci_gate_self_test_injected_failure_trips() {
        // The CI property-gate meta-test (ci-self-test.yml, foundations
        // CI/CD bead): red ONLY when the workflow injects the failure, so a
        // deliberately failing test demonstrably blocks the gate. Inert in
        // every normal run.
        assert!(
            std::env::var("FS_CI_INJECT_FAILURE").is_err(),
            "FS_CI_INJECT_FAILURE is set: this red run is the proof that a \
             failing test blocks the CI gate"
        );
    }

    #[test]
    fn ascent_and_lumen_are_siblings() {
        // LUMEN (L5) must not depend on ASCENT (L4) even though 4 < 5.
        let ms = vec![
            manifest("fs-opt", "L4", &[]),
            manifest("fs-render", "L5", &["fs-opt"]),
        ];
        let v = check_layers(&ms);
        assert_eq!(v.len(), 1, "L5 -> L4 must be rejected: {v:?}");
        // And HELM (L6) may depend on both.
        let ms = vec![
            manifest("fs-opt", "L4", &[]),
            manifest("fs-render", "L5", &[]),
            manifest("fs-ir", "L6", &["fs-opt", "fs-render"]),
        ];
        assert!(check_layers(&ms).is_empty());
    }

    #[test]
    fn seeded_forbidden_external_dependency_is_caught() {
        let ms = vec![manifest("fs-la", "L1", &["ndarray"])];
        let v = check_deps(&ms);
        assert_eq!(v.len(), 1, "expected policy violation for ndarray: {v:?}");
        assert!(v[0].detail.contains("Franken-constellation"));
    }

    #[test]
    fn constellation_dependencies_are_permitted() {
        // Real member-crate names from the probed sibling workspaces.
        let ms = vec![manifest(
            "fs-exec",
            "L0",
            &[
                "asupersync",
                "fsqlite-core",
                "fnx-algorithms",
                "fnp-dtype",
                "ft-autograd",
                "fsci-cluster",
                "fp-columnar",
            ],
        )];
        assert!(check_deps(&ms).is_empty(), "{:?}", check_deps(&ms));
    }

    #[test]
    fn prefix_matching_does_not_over_accept() {
        // Names merely RESEMBLING constellation prefixes must still be caught.
        for bad in ["ftui", "fparse", "fnxlib", "fsq", "asupersync2-evil"] {
            let ms = vec![manifest("fs-la", "L1", &[bad])];
            // ft-/fp-/fnx- require the hyphen; fsqlite requires the full stem;
            // asupersync2-evil starts_with("asupersync") — document that the
            // gate is prefix-trusting WITHIN our own manifests (review catches
            // deliberate evasion; the gate catches accidents).
            if bad == "asupersync2-evil" {
                continue; // known prefix-trust boundary, documented above
            }
            assert_eq!(check_deps(&ms).len(), 1, "should reject {bad}");
        }
    }

    #[test]
    fn unsafe_scanner_finds_tokens_and_respects_comments() {
        assert!(contains_unsafe_token(
            "fn f() { unsafe { core::hint::unreachable_unchecked() } }"
        ));
        assert!(contains_unsafe_token("unsafe fn g() {}"));
        // Prose about unsafe in comments must NOT trip the scanner.
        assert!(!contains_unsafe_token("// this crate denies unsafe code"));
        assert!(!contains_unsafe_token("/// `unsafe` is forbidden here"));
        // Identifiers merely containing the substring must not trip it.
        assert!(!contains_unsafe_token(
            "let unsafety = 1; let not_unsafe_thing = 2;"
        ));
    }

    #[test]
    fn terminology_scanner_has_an_exact_policy_allowlist_and_mutation_coverage() {
        let word = legacy_branch_word();
        let title_case = format!("{}{}", word[..1].to_ascii_uppercase(), &word[1..]);
        let upper_case = word.to_ascii_uppercase();
        let mut clean = BTreeMap::from([
            (
                "AGENTS.md".to_string(),
                format!(
                    "## Git Branch: ONLY Use `main`, NEVER `{word}`\n- Never reference `{word}` in code or docs. If you see it, treat it as a bug.\n- If the remote also needs a legacy `{word}` ref, synchronize it from `main`\n"
                ),
            ),
            (
                "docs/domain.md".to_string(),
                "The main manifold uses one primary representative and a root seed.\n".to_string(),
            ),
        ]);
        assert!(
            scan_terminology_sources(&clean).is_empty(),
            "the exact controlling policy lines and unrelated domain vocabulary are allowed"
        );

        clean.insert(
            "crates/a/src/lib.rs".to_string(),
            format!("let {word} = representative;\n"),
        );
        clean.insert("docs/b.md".to_string(), format!("{title_case} seed.\n"));
        clean.insert("scripts/c.sh".to_string(), format!("role={upper_case}\n"));
        let violations = scan_terminology_sources(&clean);
        assert_eq!(violations.len(), 3, "each spelling class must trip once");
        assert_eq!(
            violations
                .iter()
                .map(|violation| violation.crate_name.as_str())
                .collect::<Vec<_>>(),
            vec!["crates/a/src/lib.rs", "docs/b.md", "scripts/c.sh"]
        );
        assert!(violations.iter().all(|violation| {
            violation.detail.contains("allowlist v2") && violation.detail.contains(":1:")
        }));

        let near_miss = BTreeMap::from([(
            "AGENTS.md".to_string(),
            format!("- Never reference `{word}` in generated docs.\n"),
        )]);
        assert_eq!(
            scan_terminology_sources(&near_miss).len(),
            1,
            "the AGENTS path does not broadly exempt non-policy occurrences"
        );

        let compounds = BTreeMap::from([(
            "crates/d/src/lib.rs".to_string(),
            format!(
                "/// sqlite_{word} rows, the {word} equation, DMA {word}s, and bus {word} roles.\n"
            ),
        )]);
        assert!(
            scan_terminology_sources(&compounds).is_empty(),
            "established compound technical terms are not branch references"
        );

        let compound_plus_bare = BTreeMap::from([(
            "crates/e/src/lib.rs".to_string(),
            format!("/// sqlite_{word} on the {word} branch.\n"),
        )]);
        assert_eq!(
            scan_terminology_sources(&compound_plus_bare).len(),
            1,
            "a bare branch reference still trips beside an allowed compound"
        );
    }

    #[test]
    fn terminology_path_scope_excludes_issue_storage_but_includes_source_and_docs() {
        assert!(!terminology_scans_path(".beads/issues.jsonl"));
        assert!(terminology_scans_path("AGENTS.md"));
        assert!(terminology_scans_path("crates/demo/src/lib.rs"));
        assert!(terminology_scans_path("scripts/ci/check.sh"));
        assert!(!terminology_scans_path("assets/plot.png"));
    }

    #[test]
    fn terminology_check_is_clean_on_the_real_workspace_tree() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("xtask lives below the workspace root");
        let violations = check_terminology(root);
        assert!(
            violations.is_empty(),
            "tracked source/docs must satisfy the sole-branch vocabulary policy: {violations:?}"
        );
    }

    #[test]
    fn unsafe_check_end_to_end_on_fixture_tree() {
        // Build a throwaway crate tree with a seeded violation and a
        // properly registered capsule; run the real filesystem check.
        let base = std::env::temp_dir().join(format!("fsim-unsafe-test-{}", std::process::id()));
        let mk = |rel: &str, content: &str| {
            let p = base.join(rel);
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            std::fs::write(p, content).unwrap();
        };
        // Unregistered unsafe -> violation.
        mk(
            "crates/fs-bad/src/lib.rs",
            "pub fn f() { unsafe { std::hint::spin_loop() } }\n",
        );
        // Registered capsule under 300 lines with SAFETY.md -> clean.
        mk(
            "crates/fs-good/src/capsule.rs",
            "#[allow(unsafe_code)]\nmod inner { }\n",
        );
        mk(
            "crates/fs-good/src/SAFETY.md",
            "# SAFETY: fs-good/capsule\n",
        );
        mk(
            "unsafe-capsules.json",
            "{\n\"capsules\": [\n{\"crate\": \"fs-good\", \"module\": \"crates/fs-good/src/capsule.rs\", \"safety_md\": \"crates/fs-good/src/SAFETY.md\"}\n]\n}\n",
        );
        let v = check_unsafe(&base);
        assert_eq!(v.len(), 1, "exactly the seeded violation expected: {v:?}");
        assert!(v[0].detail.contains("fs-bad"), "{v:?}");
        // Oversized registered capsule -> violation.
        let big = format!("#[allow(unsafe_code)]\n{}", "// pad\n".repeat(300));
        mk("crates/fs-good/src/capsule.rs", &big);
        let v = check_unsafe(&base);
        assert!(
            v.iter().any(|x| x.detail.contains("lines")),
            "LOC cap must trip: {v:?}"
        );
    }

    #[test]
    fn powi_single_arg_distinguishes_float_from_builder_calls() {
        // Single-argument calls (candidate float powi) return the arg.
        assert_eq!(powi_single_arg("n)"), Some("n"));
        assert_eq!(
            powi_single_arg("i32::from(*k) - 1)"),
            Some("i32::from(*k) - 1")
        );
        assert_eq!(powi_single_arg("-3)"), Some("-3"));
        // Two-argument calls (typed builders) are skipped.
        assert_eq!(powi_single_arg("a, n)"), None);
        assert_eq!(powi_single_arg("base, exp)"), None);
        // Nested commas inside the single argument do not split it.
        assert_eq!(
            powi_single_arg("i32::try_from(k).unwrap_or(0))"),
            Some("i32::try_from(k).unwrap_or(0)")
        );
        // Unbalanced (multi-line) stays suspect — conservative Some.
        assert_eq!(
            powi_single_arg("2 * i32::try_from("),
            Some("2 * i32::try_from(")
        );
    }

    #[test]
    fn powi_check_end_to_end_on_fixture_tree() {
        let base = std::env::temp_dir().join(format!("fsim-powi-test-{}", std::process::id()));
        let mk = |rel: &str, content: &str| {
            let p = base.join(rel);
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            std::fs::write(p, content).unwrap();
        };
        // Seeded violations: variable method exponent, whitespace before
        // the call delimiter, a test-only call, and explicit primitive UFCS.
        mk(
            "crates/fs-bad/src/lib.rs",
            "pub fn f(x: f64, n: i32) -> f64 { x.powi(n) + x.powi (n) }\n",
        );
        mk(
            "crates/fs-bad/tests/golden.rs",
            "pub fn oracle(x: f64, n: i32) -> f64 { f64::powi(x, n) }\n",
        );
        // Clean: small literals, det-ok (same and preceding line), builder
        // two-arg calls, and prose in comments.
        mk(
            "crates/fs-good/src/lib.rs",
            concat!(
                "pub fn a(x: f64) -> f64 { x.powi(2) + x.powi(-3) }\n",
                "pub fn b(x: f64, n: i32) -> f64 { x.powi(n) } // det-ok: test fixture\n",
                "// det-ok: annotated on preceding line\n",
                "pub fn c(x: f64, n: i32) -> f64 { x.powi(n) }\n",
                "// prose mentioning .powi( in a comment is not a call\n",
                "pub fn d(p: &mut B, a: u32, n: i8) { p.powi(a, n); }\n",
            ),
        );
        // The standalone WASM workspace is part of the deterministic claim.
        mk(
            "crates/fs-wasm/src/lib.rs",
            "pub fn w(x: f64, n: i32) -> f64 { x.powi(n) }\n",
        );
        let v = check_powi(&base);
        assert_eq!(v.len(), 3, "one violation per seeded file expected: {v:?}");
        assert!(v.iter().any(|x| x.detail.contains("tests/golden")), "{v:?}");
        assert!(v.iter().any(|x| x.detail.contains("fs-wasm")), "{v:?}");
        assert!(
            v.iter().all(|x| x.detail.contains("det-ok")),
            "fix hint expected: {v:?}"
        );
    }

    #[test]
    fn lock_identity_is_canonical_and_version_extraction_works() {
        let entries = vec![
            ConstellationEntry {
                remote: "https://example.invalid/x.git".to_string(),
                lib: "a".into(),
                dir: PathBuf::from("/x/a"),
                version: "1.0.0".into(),
                git_head: "abc".into(),
            },
            ConstellationEntry {
                remote: "https://example.invalid/x.git".to_string(),
                lib: "b".into(),
                dir: PathBuf::from("/elsewhere/b"),
                version: "2.0.0".into(),
                git_head: "def".into(),
            },
        ];
        // Identity excludes paths: same libs at different paths hash equal.
        let id = lock_identity(&entries);
        assert_eq!(id, "a=1.0.0@abc\nb=2.0.0@def\n");
        assert!(
            !id.contains("/x/"),
            "paths must not enter the lock identity"
        );
        assert_eq!(
            first_version_in("[workspace.package]\nversion = \"0.2.0\"\n"),
            Some("0.2.0".to_string())
        );
        assert_eq!(first_version_in("[package]\nname = \"x\"\n"), None);
    }

    #[test]
    fn bootstrap_value_flags_refuse_missing_empty_and_option_operands() {
        for flag in ["--dest", "--from"] {
            for operand in [
                None,
                Some(""),
                Some("--offline"),
                Some("--dest"),
                Some("--from"),
            ] {
                let mut args = vec![flag];
                args.extend(operand);
                let args: Vec<String> = args.into_iter().map(str::to_string).collect();
                let error = parse_bootstrap_options(Some(PathBuf::from("/default")), &args)
                    .expect_err("malformed value-taking flag must refuse");
                assert_eq!(error, format!("{flag} requires a non-empty value"));
            }
        }

        let valid: Vec<String> = ["--dest", "/constellation", "--from", "/mirror", "--offline"]
            .into_iter()
            .map(str::to_string)
            .collect();
        assert_eq!(
            parse_bootstrap_options(Some(PathBuf::from("/default")), &valid),
            Ok(BootstrapOptions {
                dest: Some(PathBuf::from("/constellation")),
                offline: true,
                from: Some("/mirror".to_string()),
            })
        );
    }

    #[test]
    fn unborn_bootstrap_requires_marker_or_exact_locked_origin() {
        let expected = "https://example.invalid/repo";
        assert!(may_resume_unborn_checkout(true, None, expected));
        assert!(may_resume_unborn_checkout(false, Some(expected), expected));
        assert!(!may_resume_unborn_checkout(false, None, expected));
        assert!(!may_resume_unborn_checkout(
            false,
            Some("https://example.invalid/other"),
            expected
        ));
    }

    #[test]
    fn constellation_lock_decode_is_bounded_and_utf8_only() {
        assert_eq!(
            decode_constellation_lock(b"canonical".to_vec()).expect("bounded UTF-8"),
            "canonical"
        );
        assert!(decode_constellation_lock(vec![b'x'; MAX_CONSTELLATION_LOCK_BYTES + 1]).is_err());
        assert!(decode_constellation_lock(vec![0xff]).is_err());
    }

    #[test]
    fn nested_directory_does_not_inherit_ancestor_repository_identity() {
        let base = std::env::temp_dir().join(format!(
            "xtask-bootstrap-root-identity-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&base);
        let ancestor = base.join("ancestor");
        let target = ancestor.join("target");
        std::fs::create_dir_all(&target).expect("fixture directories");
        let init_ancestor = std::process::Command::new("git")
            .arg("-C")
            .arg(&ancestor)
            .args(["init", "--quiet"])
            .output()
            .expect("git init ancestor");
        assert!(init_ancestor.status.success());
        assert!(
            !is_repository_root(&target).expect("inspect empty child"),
            "an empty child must not inherit its ancestor's repository identity"
        );

        let init_target = std::process::Command::new("git")
            .arg("-C")
            .arg(&target)
            .args(["init", "--quiet"])
            .output()
            .expect("git init target");
        assert!(init_target.status.success());
        assert!(is_repository_root(&target).expect("inspect initialized target"));
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // retained matrix covers every nested concealment mode
    fn pinned_cleanliness_forces_nested_submodule_visibility() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("wall clock is after the Unix epoch")
            .as_nanos();
        let base = std::env::temp_dir().join(format!(
            "xtask-submodule-cleanliness-{}-{unique}",
            std::process::id(),
        ));
        let nested_upstream = base.join("sqlite-upstream");
        let outer = base.join("frankensqlite");
        std::fs::create_dir_all(&nested_upstream).expect("create nested upstream");
        std::fs::create_dir_all(&outer).expect("create outer repository");
        for repository in [&nested_upstream, &outer] {
            git_run(repository, &["init", "--quiet", "-b", "main"])
                .expect("initialize fixture repository");
            git_run(
                repository,
                &[
                    "config",
                    "--local",
                    "user.email",
                    "submodule@frankensim.test",
                ],
            )
            .expect("configure fixture email");
            git_run(
                repository,
                &["config", "--local", "user.name", "submodule fixture"],
            )
            .expect("configure fixture name");
            git_run(
                repository,
                &["config", "--local", "commit.gpgsign", "false"],
            )
            .expect("disable fixture signing");
        }

        std::fs::write(
            nested_upstream.join("sqlite3.c"),
            "/* pinned sqlite fixture */\n",
        )
        .expect("write nested pin");
        git_run(&nested_upstream, &["add", "sqlite3.c"]).expect("stage nested pin");
        git_run(&nested_upstream, &["commit", "--quiet", "-m", "nested pin"])
            .expect("commit nested pin");

        std::fs::write(outer.join("lib.rs"), "pub fn outer() {}\n").expect("write outer fixture");
        git_run(&outer, &["add", "lib.rs"]).expect("stage outer fixture");
        git_run(&outer, &["commit", "--quiet", "-m", "outer base"]).expect("commit outer base");
        let nested_url = nested_upstream.to_str().expect("UTF-8 fixture path");
        git_run(
            &outer,
            &[
                "-c",
                "protocol.file.allow=always",
                "submodule",
                "add",
                "--quiet",
                nested_url,
                "legacy_sqlite_code/sqlite",
            ],
        )
        .expect("add nested submodule");
        git_run(
            &outer,
            &[
                "config",
                "-f",
                ".gitmodules",
                "submodule.legacy_sqlite_code/sqlite.ignore",
                "dirty",
            ],
        )
        .expect("install committed ignore=dirty policy");
        git_run(&outer, &["add", ".gitmodules", "legacy_sqlite_code/sqlite"])
            .expect("stage gitlink");
        git_run(&outer, &["commit", "--quiet", "-m", "pin nested sqlite"]).expect("commit gitlink");

        let row = LockRow {
            lib: "frankensqlite".to_string(),
            version: "fixture".to_string(),
            git_head: git_out(&outer, &["rev-parse", "HEAD"]).expect("outer head"),
            remote: nested_url.to_string(),
            path: outer.display().to_string(),
        };
        assert_eq!(repository_worktree_status(&outer), Ok(String::new()));
        verify_pinned_clean(&row, &outer).expect("clean initialized submodule verifies");
        assert_eq!(
            bootstrap_checkout(&row, &outer, "/unused", true),
            Ok(BootstrapOutcome::pinned(false)),
            "already-pinned offline bootstrap must verify the same clean nested state"
        );

        let nested_checkout = outer.join("legacy_sqlite_code/sqlite");
        std::fs::write(
            nested_checkout.join("sqlite3.c"),
            "/* concealed tracked nested dirt */\n",
        )
        .expect("write concealed nested dirt");
        assert_eq!(
            git_out(&outer, &["status", "--porcelain"]),
            Ok(String::new()),
            "ordinary status must demonstrate committed ignore=dirty concealment"
        );
        assert!(
            !repository_worktree_status(&outer)
                .expect("forced nested status")
                .is_empty()
        );
        let dirty = verify_pinned_clean(&row, &outer)
            .expect_err("tracked nested dirt must refuse despite ignore=dirty");
        assert!(dirty.contains("DIRTY"), "{dirty}");
        assert!(dirty.contains("legacy_sqlite_code/sqlite"), "{dirty}");
        std::fs::write(
            nested_checkout.join("sqlite3.c"),
            "/* pinned sqlite fixture */\n",
        )
        .expect("restore fixture bytes");
        assert_eq!(repository_worktree_status(&outer), Ok(String::new()));

        git_run(
            &nested_checkout,
            &["update-index", "--assume-unchanged", "sqlite3.c"],
        )
        .expect("install nested assume-unchanged flag");
        std::fs::write(
            nested_checkout.join("sqlite3.c"),
            "/* concealed assume-unchanged dirt */\n",
        )
        .expect("write concealed assume-unchanged dirt");
        assert_eq!(
            git_out(&nested_checkout, &["status", "--porcelain"]),
            Ok(String::new()),
            "ordinary nested status must demonstrate index-flag concealment"
        );
        let hidden_index =
            repository_worktree_status(&outer).expect("recursive hidden-index inspection");
        assert!(hidden_index.contains("legacy_sqlite_code/sqlite"));
        assert!(hidden_index.contains("index flag hides worktree state"));
        verify_pinned_clean(&row, &outer).expect_err("nested assume-unchanged dirt must refuse");
        git_run(
            &nested_checkout,
            &["update-index", "--no-assume-unchanged", "sqlite3.c"],
        )
        .expect("clear nested assume-unchanged flag");
        std::fs::write(
            nested_checkout.join("sqlite3.c"),
            "/* pinned sqlite fixture */\n",
        )
        .expect("restore fixture after hidden-index case");
        assert_eq!(repository_worktree_status(&outer), Ok(String::new()));

        std::fs::write(
            nested_upstream.join("sqlite3.c"),
            "/* drifted sqlite fixture */\n",
        )
        .expect("write nested drift");
        git_run(&nested_upstream, &["add", "sqlite3.c"]).expect("stage nested drift");
        git_run(
            &nested_upstream,
            &["commit", "--quiet", "-m", "nested drift"],
        )
        .expect("commit nested drift");
        let drift_head =
            git_out(&nested_upstream, &["rev-parse", "HEAD"]).expect("nested drift head");
        git_run(
            &outer,
            &[
                "config",
                "--local",
                "submodule.legacy_sqlite_code/sqlite.ignore",
                "all",
            ],
        )
        .expect("install repository-local ignore=all policy");
        git_run(
            &nested_checkout,
            &["fetch", "--quiet", "origin", drift_head.as_str()],
        )
        .expect("fetch nested drift");
        git_run(
            &nested_checkout,
            &["checkout", "--quiet", "--detach", drift_head.as_str()],
        )
        .expect("checkout nested drift");
        assert_eq!(
            git_out(&outer, &["status", "--porcelain"]),
            Ok(String::new()),
            "ordinary status must demonstrate repository-local ignore=all concealment"
        );
        assert!(
            !repository_worktree_status(&outer)
                .expect("forced nested HEAD status")
                .is_empty()
        );
        let drift = bootstrap_checkout(&row, &outer, "/unused", true)
            .expect_err("nested HEAD drift must refuse despite ignore=all");
        assert!(drift.contains("DIRTY"), "{drift}");

        let local_exclude = PathBuf::from(
            git_out(
                &nested_checkout,
                &["rev-parse", "--git-path", "info/exclude"],
            )
            .expect("resolve nested local exclude"),
        );
        let local_exclude = if local_exclude.is_absolute() {
            local_exclude
        } else {
            nested_checkout.join(local_exclude)
        };
        std::fs::write(&local_exclude, "hidden-local.c\n").expect("install nested local exclusion");
        std::fs::write(
            nested_checkout.join("hidden-local.c"),
            "/* locally excluded nested source */\n",
        )
        .expect("write locally excluded nested source");
        assert_eq!(
            git_out(&nested_checkout, &["status", "--porcelain"]),
            Ok(String::new()),
            "ordinary nested status must demonstrate local-exclude concealment"
        );
        let hidden_untracked =
            repository_worktree_status(&outer).expect("recursive hidden-untracked inspection");
        assert!(hidden_untracked.contains("legacy_sqlite_code/sqlite"));
        assert!(hidden_untracked.contains("hidden-local.c"));
    }

    #[test]
    fn bootstrap_provenance_distinguishes_canonical_remote_from_selected_mirror() {
        let row = LockRow {
            lib: "fixture".to_string(),
            version: "1.0.0".to_string(),
            git_head: "a".repeat(40),
            remote: "https://canonical.invalid/fixture.git".to_string(),
            path: "/constellation/fixture".to_string(),
        };
        let receipt = bootstrap_provenance::render_bootstrap_provenance_row(
            &bootstrap_provenance_row(&row, "/airgap/mirror/fixture", true, "cloned"),
        );
        assert!(receipt.contains("\"remote\": \"https://canonical.invalid/fixture.git\""));
        assert!(receipt.contains("\"selected_transport\": \"/airgap/mirror/fixture\""));
        assert!(receipt.contains("\"transport_used\": true"));
        assert_ne!(row.remote, "/airgap/mirror/fixture");
    }

    #[test]
    fn bootstrap_outcomes_match_standalone_state_and_transport_semantics() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("wall clock is after the Unix epoch")
            .as_nanos();
        let base = std::env::temp_dir().join(format!(
            "xtask-bootstrap-outcome-{}-{unique}",
            std::process::id(),
        ));
        let upstream = base.join("upstream");
        std::fs::create_dir_all(&upstream).expect("create isolated upstream");
        git_run(&upstream, &["init", "--quiet"]).expect("initialize upstream");
        git_run(
            &upstream,
            &[
                "config",
                "--local",
                "user.email",
                "bootstrap@frankensim.test",
            ],
        )
        .expect("configure upstream email");
        git_run(
            &upstream,
            &["config", "--local", "user.name", "bootstrap parity"],
        )
        .expect("configure upstream name");
        git_run(&upstream, &["config", "--local", "commit.gpgsign", "false"])
            .expect("disable fixture signing");
        std::fs::write(upstream.join("fixture.rs"), "pub fn pinned() {}\n")
            .expect("write pinned fixture");
        git_run(&upstream, &["add", "fixture.rs"]).expect("stage pinned fixture");
        git_run(&upstream, &["commit", "--quiet", "-m", "pinned"]).expect("commit pinned fixture");
        let pinned_head = git_out(&upstream, &["rev-parse", "HEAD"]).expect("read fixture head");
        let selected = format!("file://{}", upstream.display());
        let row = LockRow {
            lib: "fixture".to_string(),
            version: "1.0.0".to_string(),
            git_head: pinned_head,
            remote: selected.clone(),
            path: "/constellation/fixture".to_string(),
        };

        let marked_target = base.join("marked-at-pin");
        std::fs::create_dir_all(&marked_target).expect("create marked destination");
        git_run(&marked_target, &["init", "--quiet"]).expect("initialize marked destination");
        git_run(
            &marked_target,
            &["remote", "add", "origin", selected.as_str()],
        )
        .expect("configure marked origin");
        git_run(
            &marked_target,
            &[
                "fetch",
                "--quiet",
                "--depth",
                "1",
                "origin",
                row.git_head.as_str(),
            ],
        )
        .expect("fetch marked pin");
        git_run(
            &marked_target,
            &["checkout", "--quiet", "--detach", row.git_head.as_str()],
        )
        .expect("checkout marked pin");
        git_run(
            &marked_target,
            &["config", "--local", BOOTSTRAP_INCOMPLETE_KEY, "true"],
        )
        .expect("mark interrupted checkout");
        let marked_at_pin = bootstrap_checkout(
            &row,
            &marked_target,
            "/unreachable-transport-must-not-be-used",
            false,
        )
        .expect("already-pinned marked checkout resumes without transport");
        assert_eq!(marked_at_pin.state, "resumed");
        assert!(!marked_at_pin.transport_used);
        assert!(
            git_out(
                &marked_target,
                &["config", "--local", "--get", BOOTSTRAP_INCOMPLETE_KEY],
            )
            .is_err(),
            "successful pinned replay clears the incomplete marker",
        );
        assert_eq!(
            bootstrap_provenance::render_bootstrap_provenance_row(&bootstrap_provenance_row(
                &row,
                "/unreachable-transport-must-not-be-used",
                marked_at_pin.transport_used,
                marked_at_pin.state,
            )),
            format!(
                "{{\"lib\": \"fixture\", \"git_head\": \"{}\", \"remote\": \"{}\", \"selected_transport\": \"/unreachable-transport-must-not-be-used\", \"transport_used\": false, \"state\": \"resumed\"}}",
                row.git_head, row.remote,
            ),
            "clearing an already-pinned marker must not claim network or mirror use",
        );

        let empty_target = base.join("pre-existing-empty");
        std::fs::create_dir_all(&empty_target).expect("create empty destination");
        let initialized_empty = bootstrap_checkout(&row, &empty_target, &selected, false)
            .expect("pre-existing empty destination materializes at the pin");
        assert_eq!(initialized_empty.state, "cloned");
        assert!(initialized_empty.transport_used);
        assert_eq!(
            git_out(
                &empty_target,
                &["config", "--local", "--get", "core.autocrlf"],
            )
            .expect("initialized repository pins checkout byte policy"),
            "false",
        );
        assert_eq!(
            bootstrap_provenance::render_bootstrap_provenance_row(&bootstrap_provenance_row(
                &row,
                &selected,
                initialized_empty.transport_used,
                initialized_empty.state,
            )),
            format!(
                "{{\"lib\": \"fixture\", \"git_head\": \"{}\", \"remote\": \"{}\", \"selected_transport\": \"{}\", \"transport_used\": true, \"state\": \"cloned\"}}",
                row.git_head, row.remote, selected,
            ),
            "a pre-existing empty directory initialized by this invocation is a clone, not a resume",
        );
    }

    #[test]
    fn bootstrap_pinned_tree_observation_requires_exact_head_and_clean_status() {
        let row = LockRow {
            lib: "fixture".to_string(),
            version: "1.0.0".to_string(),
            git_head: "locked-head".to_string(),
            remote: "unused".to_string(),
            path: "/constellation/fixture".to_string(),
        };
        let target = Path::new("/constellation/fixture");

        assert!(check_pinned_clean_observation(&row, target, "locked-head", "").is_ok());

        let drift = check_pinned_clean_observation(&row, target, "other-head", "")
            .expect_err("head drift must refuse");
        assert!(drift.contains("lock pins locked-head"), "{drift}");

        let dirty = check_pinned_clean_observation(&row, target, "locked-head", " M lib.rs")
            .expect_err("checkout-time dirt must refuse");
        assert!(dirty.contains("DIRTY"), "{dirty}");
    }

    #[test]
    fn pinned_tree_acceptance_rechecks_status_after_the_first_clean_observation() {
        let row = LockRow {
            lib: "fixture".to_string(),
            version: "1.0.0".to_string(),
            git_head: "locked-head".to_string(),
            remote: "unused".to_string(),
            path: "/constellation/fixture".to_string(),
        };
        let target = Path::new("/constellation/fixture");
        let clean = RepositoryObservation {
            head: row.git_head.clone(),
            status: String::new(),
        };
        let dirty = RepositoryObservation {
            head: row.git_head.clone(),
            status: " M src/lib.rs".to_string(),
        };
        let mut stable_observations = [clean.clone(), clean.clone()].into_iter();
        verify_pinned_clean_with(&row, target, || {
            Ok(stable_observations
                .next()
                .expect("verifier must take exactly two observations"))
        })
        .expect("two identical clean pinned observations must pass");
        assert!(stable_observations.next().is_none());

        let mut observations = [clean, dirty].into_iter();
        let error = verify_pinned_clean_with(&row, target, || {
            Ok(observations
                .next()
                .expect("verifier must take exactly two observations"))
        })
        .expect_err("dirt appearing after the first observation must refuse provenance");
        assert!(error.contains("DIRTY"), "{error}");
        assert!(
            observations.next().is_none(),
            "confirmation observation must be consumed before acceptance"
        );

        let moved = coherent_repository_observation(
            target,
            "locked-head",
            String::new(),
            "other-head".to_string(),
        )
        .expect_err("HEAD movement inside one status observation must refuse");
        assert!(
            moved.contains("worktree state was being observed"),
            "{moved}"
        );
    }

    #[test]
    fn constellation_lock_identity_is_ordered_complete_and_duplicate_safe() {
        let first = LockRow {
            lib: "zeta".to_string(),
            version: "2.0.0".to_string(),
            git_head: "2222".to_string(),
            remote: "unused-zeta".to_string(),
            path: "/constellation/zeta".to_string(),
        };
        let second = LockRow {
            lib: "alpha".to_string(),
            version: "1.0.0".to_string(),
            git_head: "1111".to_string(),
            remote: "unused-alpha".to_string(),
            path: "/constellation/alpha".to_string(),
        };
        let identity = lock_rows_identity(&[first, second]).expect("unique rows are canonical");
        assert_eq!(identity, "alpha=1.0.0@1111\nzeta=2.0.0@2222\n");

        let duplicate = [
            LockRow {
                lib: "same".to_string(),
                version: "1".to_string(),
                git_head: "a".to_string(),
                remote: "one".to_string(),
                path: "/constellation/same-a".to_string(),
            },
            LockRow {
                lib: "same".to_string(),
                version: "2".to_string(),
                git_head: "b".to_string(),
                remote: "two".to_string(),
                path: "/constellation/same-b".to_string(),
            },
        ];
        assert!(lock_rows_identity(&duplicate).is_err());
    }

    #[test]
    fn tracked_constellation_lock_is_self_consistent_and_schema_bound() {
        let tracked = include_str!("../../constellation.lock");
        assert_eq!(CONSTELLATION_LOCK_IDENTITY_VERSION, 1);
        assert_eq!(CONSTELLATION_LOCK_WRITER_IDENTITY_VERSION, 2);
        assert_ne!(
            CONSTELLATION_LOCK_IDENTITY_DOMAIN,
            CONSTELLATION_LOCK_WRITER_IDENTITY_DOMAIN
        );
        assert!(!tracked.contains(CONSTELLATION_LOCK_WRITER_IDENTITY_DOMAIN));
        let (recorded, mut rows) = parse_lock_rows(tracked).expect("tracked lock is canonical");
        let identity = lock_rows_identity(&rows).expect("tracked rows are unique");
        assert_eq!(recorded, format!("{:016x}", fnv1a64(identity.as_bytes())));

        rows[0].git_head = "0".repeat(rows[0].git_head.len());
        let changed = lock_rows_identity(&rows).expect("mutated rows remain structurally valid");
        assert_ne!(recorded, format!("{:016x}", fnv1a64(changed.as_bytes())));

        let wrong_schema = tracked.replacen(
            "frankensim-constellation-lock-v2",
            "frankensim-constellation-lock-v3",
            1,
        );
        assert!(parse_lock_rows(&wrong_schema).is_err());
        let wrong_identity_domain = tracked.replacen(
            CONSTELLATION_LOCK_IDENTITY_DOMAIN,
            "org.frankensim.xtask.constellation-lock.v0",
            1,
        );
        assert!(parse_lock_rows(&wrong_identity_domain).is_err());
        let wrong_identity_version = tracked.replacen(
            &format!("\"identity_version\": {CONSTELLATION_LOCK_IDENTITY_VERSION}"),
            "\"identity_version\": 0",
            1,
        );
        assert!(parse_lock_rows(&wrong_identity_version).is_err());
        assert!(parse_lock_rows(&format!("{tracked}trailing-junk\n")).is_err());

        let escaped = render_lockfile(
            &[ConstellationEntry {
                lib: "fixture".to_string(),
                dir: PathBuf::from("/tmp/quoted\"path"),
                version: "1.0".to_string(),
                git_head: "a".repeat(40),
                remote: "https://example.invalid/quoted\"remote".to_string(),
            }],
            0,
        );
        assert!(escaped.contains("quoted\\\"path"));
        assert!(escaped.contains("quoted\\\"remote"));
    }

    #[test]
    fn constellation_lock_identity_metadata_matrix_refuses_before_pin_access() {
        let tracked = include_str!("../../constellation.lock");
        let domain_line =
            format!("  \"identity_domain\": \"{CONSTELLATION_LOCK_IDENTITY_DOMAIN}\",\n");
        let version_line =
            format!("  \"identity_version\": {CONSTELLATION_LOCK_IDENTITY_VERSION},\n");
        let malformed = [
            tracked.replacen(&domain_line, "", 1),
            tracked.replacen(&domain_line, &format!("{domain_line}{domain_line}"), 1),
            tracked.replacen(
                &format!("\"identity_domain\": \"{CONSTELLATION_LOCK_IDENTITY_DOMAIN}\""),
                "\"identity_domain\": 1",
                1,
            ),
            tracked.replacen(&version_line, "", 1),
            tracked.replacen(&version_line, &format!("{version_line}{version_line}"), 1),
            tracked.replacen(
                &format!("\"identity_version\": {CONSTELLATION_LOCK_IDENTITY_VERSION}"),
                "\"identity_version\": true",
                1,
            ),
            tracked.replacen(
                &version_line,
                &format!("{version_line}  \"identity_epoch\": 1,\n"),
                1,
            ),
        ];
        for document in malformed {
            assert!(
                parse_lock_rows(&document).is_err(),
                "identity metadata defect reached a lock row: {document}"
            );
        }
    }

    #[test]
    fn constellation_lock_json_round_trips_escaped_fields_and_refuses_malformed_documents() {
        let tracked = include_str!("../../constellation.lock");
        let (recorded, mut rows) = parse_lock_rows(tracked).expect("tracked lock is canonical");
        rows[0].remote = "https://example.invalid/quoted\"\\remote/é".to_string();
        rows[0].path = "/constellation/quoted\"\\path/é".to_string();
        let rendered = render_lock_rows(&rows, &recorded);
        assert!(rendered.contains("quoted\\\"\\\\remote/é"));
        let (reparsed_hash, reparsed_rows) =
            parse_lock_rows(&rendered).expect("escaped canonical lock must parse");
        assert_eq!(reparsed_hash, recorded);
        assert_eq!(reparsed_rows, rows);

        let malformed = [
            format!("{tracked}trailing-junk\n"),
            tracked.trim_end().to_string(),
            tracked.replacen("\"schema\":", "\"extra\": 0, \"schema\":", 1),
            tracked.replacen(
                "  \"libraries\": [\n",
                "  \"libraries\": [\n    \"junk\",\n",
                1,
            ),
            tracked.replacen(
                "\"lib\": \"asupersync\",",
                "\"lib\": \"asupersync\", \"lib\": \"asupersync\",",
                1,
            ),
            tracked.replacen("\"remote\": \"https", "\"remote\": \"bad\\qhttps", 1),
            tracked.replacen("\"remote\": \"https", "\"remote\": \"\\ud800https", 1),
            tracked.replacen("\n  ]", ",\n  ]", 1),
            tracked.replacen(
                "\"lib\": \"asupersync\", \"version\":",
                "\"version\": \"0.3.5\", \"lib\": \"asupersync\", \"version\":",
                1,
            ),
            tracked.replacen("https://", "https\\u003a//", 1),
        ];
        for document in malformed {
            assert!(
                parse_lock_rows(&document).is_err(),
                "malformed or non-canonical lock unexpectedly parsed: {document}"
            );
        }
    }

    #[test]
    fn tool_crates_are_never_dependable() {
        let ms = vec![
            manifest("xtask", "TOOL", &[]),
            manifest("fs-la", "L1", &["xtask"]),
        ];
        assert_eq!(check_layers(&ms).len(), 1);
    }

    /// A manifest whose crate dir is a real temp directory holding the
    /// given CONTRACT.md text (or none), for exercising the contract lint
    /// against complete/incomplete fixtures (bead huq.5 acceptance).
    fn contract_fixture(name: &str, layer: &str, contract: Option<&str>) -> Manifest {
        let dir = std::env::temp_dir().join(format!(
            "fs-xtask-contract-fixture-{}-{name}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).expect("fixture dir");
        match contract {
            Some(text) => std::fs::write(dir.join("CONTRACT.md"), text).expect("fixture write"),
            None => {
                let _ = std::fs::remove_file(dir.join("CONTRACT.md"));
            }
        }
        let mut m = manifest(name, layer, &[]);
        m.dir = dir;
        m
    }

    #[test]
    fn contracts_lint_accepts_complete_and_names_every_gap() {
        let complete: String = CONTRACT_SECTIONS
            .iter()
            .map(|s| format!("{s}\n\nbody\n\n"))
            .collect();
        let ok = contract_fixture("fs-ctest-ok", "L0", Some(&complete));
        assert!(
            check_contracts(&[ok]).is_empty(),
            "a complete CONTRACT.md must lint clean"
        );

        let missing_file = contract_fixture("fs-ctest-none", "L0", None);
        let v = check_contracts(&[missing_file]);
        assert_eq!(v.len(), 1);
        assert!(v[0].detail.contains("missing CONTRACT.md"));

        for dropped in CONTRACT_SECTIONS {
            let partial: String = CONTRACT_SECTIONS
                .iter()
                .filter(|s| *s != dropped)
                .map(|s| format!("{s}\n\nbody\n\n"))
                .collect();
            let m = contract_fixture("fs-ctest-partial", "L0", Some(&partial));
            let v = check_contracts(&[m]);
            assert_eq!(
                v.len(),
                1,
                "dropping exactly {dropped:?} must produce exactly one violation"
            );
            assert!(
                v[0].detail.contains(dropped),
                "the violation must name the missing section {dropped:?}"
            );
        }
    }

    #[test]
    fn contracts_lint_exempts_tool_crates_only() {
        let tool = contract_fixture("fs-ctest-tool", "TOOL", None);
        assert!(
            check_contracts(&[tool]).is_empty(),
            "TOOL crates are exempt from the contract lint"
        );
        let util = contract_fixture("fs-ctest-util", "UTIL", None);
        assert_eq!(
            check_contracts(&[util]).len(),
            1,
            "UTIL crates are NOT exempt"
        );
    }
}
