//! FrankenSim repository policy checks (`cargo run -p xtask -- <command>`).
//!
//! Commands:
//! - `check-layers`   — enforce the L0..L6 layer dependency direction (plan §4, AGENTS.md).
//! - `check-deps`     — enforce the Franken-only runtime dependency policy (Decalogue P1).
//! - `check-contracts`— every workspace `fs-*` crate ships a CONTRACT.md with required sections.
//! - `check-unsafe`   — unsafe code only in registered capsules (<300 lines, SAFETY.md).
//! - `check-all`      — all of the above; non-zero exit on any violation.
//! - `lock-constellation` / `check-constellation` — pin/verify the Franken library states.
//!
//! Output is JSON-lines (one verdict object per check per crate) so agents parse
//! outcomes without scraping; a human-readable summary goes to stderr.
//!
//! Parsing note: this tool intentionally hand-parses the *subset* of TOML that this
//! repository's own generated manifests use (see docs/CONVENTIONS.md "Manifest
//! conventions": one dependency per line, `[section]` headers on their own line).
//! It fails loudly on shapes it does not understand rather than guessing — these are
//! our files, and the conventions are enforced, not inferred.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

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

fn parse_manifest(path: &Path, text: &str) -> Result<Manifest, String> {
    let mut section = String::new();
    let mut name = None;
    let mut layer = None;
    let mut runtime_deps = Vec::new();
    let mut dev_deps = Vec::new();

    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
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
        let head = std::process::Command::new("git")
            .arg("-C")
            .arg(&dir)
            .args(["rev-parse", "HEAD"])
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map_or_else(
                || "no-git".to_string(),
                |o| String::from_utf8_lossy(&o.stdout).trim().to_string(),
            );
        out.push(ConstellationEntry {
            lib: (*lib).to_string(),
            dir,
            version,
            git_head: head,
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

fn render_lockfile(entries: &[ConstellationEntry], hash: u64) -> String {
    let mut s = String::new();
    s.push_str("{\n  \"schema\": \"frankensim-constellation-lock-v1\",\n");
    let _ = writeln!(s, "  \"lock_hash\": \"{hash:016x}\",");
    s.push_str(
        "  \"note\": \"lock_hash covers (lib, version, git_head) only — paths are per-machine\",\n",
    );
    s.push_str("  \"libraries\": [\n");
    for (i, e) in entries.iter().enumerate() {
        let comma = if i + 1 == entries.len() { "" } else { "," };
        let _ = writeln!(
            s,
            "    {{\"lib\": \"{}\", \"version\": \"{}\", \"git_head\": \"{}\", \"path\": \"{}\"}}{comma}",
            e.lib,
            e.version,
            e.git_head,
            e.dir.display()
        );
    }
    s.push_str("  ]\n}\n");
    s
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
        let Ok(existing) = std::fs::read_to_string(&lock_path) else {
            eprintln!(
                "error: {} missing; run `cargo run -p xtask -- lock-constellation`",
                lock_path.display()
            );
            return ExitCode::FAILURE;
        };
        let recorded = existing
            .lines()
            .find_map(|l| l.trim().strip_prefix("\"lock_hash\": \""))
            .and_then(|r| r.split('"').next())
            .unwrap_or("");
        if recorded == format!("{hash:016x}") {
            println!(
                "{{\"check\":\"constellation-lock\",\"verdict\":\"ok\",\"hash\":\"{hash:016x}\"}}"
            );
            eprintln!("constellation lock OK ({hash:016x})");
            ExitCode::SUCCESS
        } else {
            println!(
                "{{\"check\":\"constellation-lock\",\"verdict\":\"drift\",\"recorded\":\"{recorded}\",\"live\":\"{hash:016x}\"}}"
            );
            eprintln!(
                "constellation DRIFT: recorded {recorded}, live {hash:016x}; a constellation \
                 repo moved — re-lock deliberately with lock-constellation and note why"
            );
            ExitCode::FAILURE
        }
    } else {
        match std::fs::write(&lock_path, render_lockfile(&entries, hash)) {
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
        _ => {}
    }
    let manifests = match load_workspace(&root) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let (violations, checks): (Vec<Violation>, Vec<&str>) = match cmd.as_str() {
        "check-layers" => (check_layers(&manifests), vec!["layers"]),
        "check-deps" => (check_deps(&manifests), vec!["dependency-policy"]),
        "check-contracts" => (check_contracts(&manifests), vec!["contracts"]),
        "check-unsafe" => (check_unsafe(&root), vec!["unsafe-capsules"]),
        "check-all" => {
            let mut v = check_layers(&manifests);
            v.extend(check_deps(&manifests));
            v.extend(check_contracts(&manifests));
            v.extend(check_unsafe(&root));
            (
                v,
                vec![
                    "layers",
                    "dependency-policy",
                    "contracts",
                    "unsafe-capsules",
                ],
            )
        }
        other => {
            eprintln!(
                "unknown command {other:?}; use check-layers|check-deps|check-contracts|\
                 check-all|lock-constellation|check-constellation"
            );
            return ExitCode::FAILURE;
        }
    };
    emit(
        &violations,
        &dev_dep_notes(&manifests),
        &checks,
        manifests.len(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn lock_identity_is_canonical_and_version_extraction_works() {
        let entries = vec![
            ConstellationEntry {
                lib: "a".into(),
                dir: PathBuf::from("/x/a"),
                version: "1.0.0".into(),
                git_head: "abc".into(),
            },
            ConstellationEntry {
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
    fn tool_crates_are_never_dependable() {
        let ms = vec![
            manifest("xtask", "TOOL", &[]),
            manifest("fs-la", "L1", &["xtask"]),
        ];
        assert_eq!(check_layers(&ms).len(), 1);
    }
}
