//! FrankenSim repository policy checks (`cargo run -p xtask -- <command>`).
//!
//! Commands:
//! - `check-layers`   — enforce the L0..L6 layer dependency direction (plan §4, AGENTS.md).
//! - `check-deps`     — enforce the Franken-only runtime dependency policy (Decalogue P1).
//! - `check-contracts`— every workspace `fs-*` crate ships a CONTRACT.md with required sections.
//! - `check-unsafe`   — unsafe code only in registered capsules (<300 lines, SAFETY.md).
//! - `check-powi`     — no build-mode-dependent `f64::powi` in deterministic paths (bead 4xnt).
//! - `check-goldens`  — golden hashes declare upstream couplings; drift re-freezes deliberately (bead y4pt).
//! - `check-claims`   — README hashes/crates/sentinels must exist in code (bead 06yc).
//! - `check-closures` — closed bug beads must cite regression evidence or a disposition (bead hx4p).
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

mod claims;
mod closures;
mod depgraph;

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
        let remote = std::process::Command::new("git")
            .arg("-C")
            .arg(&dir)
            .args(["remote", "get-url", "origin"])
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map_or_else(
                || "no-remote".to_string(),
                |o| String::from_utf8_lossy(&o.stdout).trim().to_string(),
            );
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

fn render_lockfile(entries: &[ConstellationEntry], hash: u64) -> String {
    let mut s = String::new();
    s.push_str("{\n  \"schema\": \"frankensim-constellation-lock-v2\",\n");
    let _ = writeln!(s, "  \"lock_hash\": \"{hash:016x}\",");
    s.push_str(
        "  \"note\": \"lock_hash covers (lib, version, git_head) only — paths are per-machine; \
         remote is transport for bootstrap-constellation (content identity is the git head)\",\n",
    );
    s.push_str("  \"libraries\": [\n");
    for (i, e) in entries.iter().enumerate() {
        let comma = if i + 1 == entries.len() { "" } else { "," };
        let _ = writeln!(
            s,
            "    {{\"lib\": \"{}\", \"version\": \"{}\", \"git_head\": \"{}\", \"remote\": \"{}\", \"path\": \"{}\"}}{comma}",
            e.lib,
            e.version,
            e.git_head,
            e.remote,
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
        let tag = format!("\"{key}\": ");
        let start = line.find(&tag)? + tag.len();
        let rest = &line[start..];
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
            let (Some(id), Some(file), Some(name), Some(ver)) = (
                field(line, "id"),
                field(line, "file"),
                field(line, "const"),
                field(line, "version"),
            ) else {
                violations.push(bail(format!("malformed surface row: {line}")));
                continue;
            };
            let Ok(reg_ver) = ver.parse::<u32>() else {
                violations.push(bail(format!("surface {id}: bad version {ver:?}")));
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

struct LockRow {
    lib: String,
    git_head: String,
    remote: String,
}

fn parse_lock_rows(text: &str) -> Result<(String, Vec<LockRow>), String> {
    let mut rows = Vec::new();
    let mut lock_hash = String::new();
    for line in text.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("\"lock_hash\": \"") {
            lock_hash = rest.split('"').next().unwrap_or("").to_string();
        }
        if t.starts_with("{\"lib\"") {
            let field = |key: &str| -> Option<String> {
                let tag = format!("\"{key}\": \"");
                let start = t.find(&tag)? + tag.len();
                t[start..].split('"').next().map(str::to_string)
            };
            let (Some(lib), Some(git_head)) = (field("lib"), field("git_head")) else {
                return Err(format!("malformed lock row: {t}"));
            };
            rows.push(LockRow {
                lib,
                git_head,
                remote: field("remote").unwrap_or_else(|| "no-remote".to_string()),
            });
        }
    }
    if lock_hash.is_empty() || rows.is_empty() {
        return Err("constellation.lock has no hash or no libraries".to_string());
    }
    Ok((lock_hash, rows))
}

fn git_out(dir: &Path, args: &[&str]) -> Result<String, String> {
    let out = std::process::Command::new("git")
        .arg("-C")
        .arg(dir)
        .args(args)
        .output()
        .map_err(|e| format!("git {args:?} failed to spawn: {e}"))?;
    if out.status.success() {
        Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
    } else {
        Err(format!(
            "git {args:?} in {} failed: {}",
            dir.display(),
            String::from_utf8_lossy(&out.stderr).trim()
        ))
    }
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
             replace that sibling deliberately",
            target.display()
        ));
    }
    Ok(())
}

fn verify_pinned_clean(row: &LockRow, target: &Path) -> Result<(), String> {
    let head = git_out(target, &["rev-parse", "HEAD"])
        .map_err(|e| format!("{}: {e}", target.display()))?;
    if head != row.git_head {
        return check_pinned_clean_observation(row, target, &head, "");
    }
    let status = git_out(target, &["status", "--porcelain"])?;
    check_pinned_clean_observation(row, target, &head, &status)
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
    let lock_text = match std::fs::read_to_string(root.join("constellation.lock")) {
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
        let state: Result<&'static str, String> = if target.is_dir() {
            // EXISTING tree: verify identity; never silently substitute.
            verify_pinned_clean(row, &target).map(|()| "verified")
        } else if options.offline {
            Err(format!(
                "{} missing from the source cache in --offline mode",
                target.display()
            ))
        } else {
            // FETCH: clone the declared transport, check out the pinned
            // revision detached, then apply the same identity and cleanliness
            // verifier as the existing-tree path.
            let url = options
                .from
                .as_ref()
                .map_or_else(|| row.remote.clone(), |b| format!("{b}/{dirname}"));
            if url == "no-remote" {
                Err(format!(
                    "lock declares no remote for {} — re-lock on a host that has one",
                    row.lib
                ))
            } else {
                let clone = std::process::Command::new("git")
                    .args(["clone", "--no-checkout", "-c", "core.autocrlf=false", &url])
                    .arg(&target)
                    .output();
                match clone {
                    Ok(o) if o.status.success() => {
                        match git_out(&target, &["checkout", "--detach", &row.git_head]) {
                            Ok(_) => verify_pinned_clean(row, &target).map(|()| "cloned"),
                            Err(e) => Err(format!(
                                "locked revision {} unavailable from {url}: {e}",
                                row.git_head
                            )),
                        }
                    }
                    Ok(o) => Err(format!(
                        "clone of {url} failed: {}",
                        String::from_utf8_lossy(&o.stderr).trim()
                    )),
                    Err(e) => Err(format!("git clone failed to spawn: {e}")),
                }
            }
        };
        match state {
            Ok(st) => {
                println!(
                    "{{\"check\":\"constellation-bootstrap\",\"lib\":\"{}\",\"state\":\"{st}\",\"head\":\"{}\"}}",
                    row.lib, row.git_head
                );
                provenance.push(format!(
                    "{{\"lib\": \"{}\", \"git_head\": \"{}\", \"remote\": \"{}\", \"state\": \"{st}\"}}",
                    row.lib, row.git_head, row.remote
                ));
            }
            Err(why) => {
                println!(
                    "{{\"check\":\"constellation-bootstrap\",\"lib\":\"{}\",\"state\":\"failed\",\"why\":\"{}\"}}",
                    row.lib,
                    why.replace('"', "'")
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
    // Build provenance: the lock hash + every fetched/verified identity.
    let prov = format!(
        "{{\n\"schema\": \"frankensim-constellation-bootstrap-v1\",\n\"lock_hash\": \"{lock_hash}\",\n\"dest\": \"{}\",\n\"libraries\": [\n{}\n]\n}}\n",
        dest.display(),
        provenance.join(",\n")
    );
    let prov_path = dest.join("constellation-bootstrap.json");
    if let Err(e) = std::fs::write(&prov_path, prov) {
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
        "check-powi" => (check_powi(&root), vec!["powi-determinism"]),
        "check-goldens" => (check_goldens(&root), vec!["golden-couplings"]),
        "check-claims" => (claims::check_claims(&root), vec!["claim-state"]),
        "check-closures" => (closures::check_closures(&root), vec!["closure-evidence"]),
        "check-all" => {
            let mut v = check_layers(&manifests);
            v.extend(check_deps(&manifests));
            v.extend(check_contracts(&manifests));
            v.extend(check_unsafe(&root));
            v.extend(check_powi(&root));
            v.extend(check_goldens(&root));
            v.extend(claims::check_claims(&root));
            v.extend(closures::check_closures(&root));
            (
                v,
                vec![
                    "layers",
                    "dependency-policy",
                    "contracts",
                    "unsafe-capsules",
                    "powi-determinism",
                    "golden-couplings",
                    "claim-state",
                    "closure-evidence",
                ],
            )
        }
        other => {
            eprintln!(
                "unknown command {other:?}; use check-layers|check-deps|check-contracts|\
                 check-unsafe|check-powi|check-goldens|check-claims|check-closures|\
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
            "{\n\"surfaces\": [\n{\"id\": \"mini:semantics\", \"file\": \"crates/mini/src/lib.rs\", \"const\": \"MINI_SEMANTICS_VERSION\", \"version\": 1}\n],\n\"goldens\": [\n{\"golden\": \"mini:golden\", \"file\": \"crates/mini/src/lib.rs\", \"const\": \"GOLDEN_HASH\", \"depends_on\": \"mini:semantics=1\", \"justification\": \"recorded at fixture landing, both modes, committed tree\"}\n]\n}\n",
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
        let _ = std::fs::remove_dir_all(&base);
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
    fn bootstrap_pinned_tree_observation_requires_exact_head_and_clean_status() {
        let row = LockRow {
            lib: "fixture".to_string(),
            git_head: "locked-head".to_string(),
            remote: "unused".to_string(),
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
    fn tool_crates_are_never_dependable() {
        let ms = vec![
            manifest("xtask", "TOOL", &[]),
            manifest("fs-la", "L1", &["xtask"]),
        ];
        assert_eq!(check_layers(&ms).len(), 1);
    }
}
