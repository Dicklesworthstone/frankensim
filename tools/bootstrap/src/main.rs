//! frankensim-bootstrap (bead 1t8i): the clean-machine constellation
//! bootstrap.
//!
//! A fresh `git clone` of frankensim cannot build: every workspace
//! manifest declares fixed relative path dependencies on sibling
//! repositories (`../asupersync`, `../franken_numpy`, …), and Cargo
//! resolves those paths before it will build ANYTHING in the workspace —
//! including xtask, where the in-workspace verifier lives. This package
//! is therefore deliberately NOT a workspace member: it builds alone
//! (`cargo run --manifest-path tools/bootstrap/Cargo.toml`), reads
//! `constellation.lock`, and materializes every pinned sibling next to
//! the workspace so the fixed relative paths resolve — that sibling
//! layout IS the reproducible Cargo configuration; no config files are
//! generated or mutated.
//!
//! Trust rules (all fail closed):
//! - An EXISTING sibling is verified: head must equal the lock pin and
//!   the tree must be clean. Drift and dirt are refusals, never
//!   silently substituted — a case-folding checkout collision (the
//!   7n2n counterexample) surfaces here as a dirty tree and refuses.
//! - A MISSING sibling is cloned from the lock's declared remote (or
//!   `--from <base>/<dirname>` for air-gapped mirrors), checked out
//!   DETACHED at the pinned revision, then subjected to the same pinned-head
//!   and clean-tree verification as an existing sibling. No branches or
//!   worktrees are created anywhere.
//! - `--offline` never touches the network: missing siblings are
//!   structured failures (the offline-cache replay contract).
//! - Idempotent: a second run over a successful first run verifies
//!   every sibling and rewrites identical provenance.
//!
//! Output: one JSON line per library plus
//! `constellation-bootstrap.json` (schema
//! `frankensim-constellation-bootstrap-v1`) beside the siblings. The
//! logic mirrors `cargo run -p xtask -- bootstrap-constellation`, which
//! remains the in-workspace verifier once the workspace can build; this
//! binary is the pre-Cargo entry point for machines that cannot.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Library name → sibling directory name (identity mapping today; kept
/// explicit so a future rename cannot silently retarget a clone).
const CONSTELLATION_REPOS: &[(&str, &str)] = &[
    ("asupersync", "asupersync"),
    ("frankensqlite", "frankensqlite"),
    ("franken_numpy", "franken_numpy"),
    ("frankentorch", "frankentorch"),
    ("frankenscipy", "frankenscipy"),
    ("frankenpandas", "frankenpandas"),
    ("franken_networkx", "franken_networkx"),
];

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

fn dirname_of(lib: &str) -> &str {
    CONSTELLATION_REPOS
        .iter()
        .find(|(l, _)| *l == lib)
        .map_or(lib, |(_, d)| d)
}

fn usage() -> &'static str {
    "frankensim-bootstrap [--root <frankensim-checkout>] [--offline] [--from <mirror-base>]\n\
     \n\
     Reads <root>/constellation.lock and materializes every pinned sibling\n\
     repository in <root>'s PARENT directory (where the workspace's fixed\n\
     relative path dependencies point). Existing siblings are verified\n\
     (pinned head + clean tree) and never silently substituted."
}

fn required_option_value<'a>(flag: &str, value: Option<&'a String>) -> Result<&'a str, String> {
    match value {
        Some(value) if !value.is_empty() && !value.starts_with('-') => Ok(value),
        _ => Err(format!("{flag} requires a non-empty value")),
    }
}

fn verify_pinned_clean(row: &LockRow, target: &Path) -> Result<(), String> {
    let head = git_out(target, &["rev-parse", "HEAD"])
        .map_err(|e| format!("{}: {e}", target.display()))?;
    if head != row.git_head {
        return Err(format!(
            "{} is at {head}, lock pins {} — refusing to silently substitute a nearby \
             working tree; align or replace that sibling deliberately",
            target.display(),
            row.git_head
        ));
    }
    let status = git_out(target, &["status", "--porcelain"])?;
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

/// One library's bootstrap: verify an existing tree or clone a missing
/// one. Returns the terminal state name or the structured refusal.
fn bootstrap_one(
    row: &LockRow,
    dest: &Path,
    offline: bool,
    from: Option<&str>,
) -> Result<&'static str, String> {
    let dirname = dirname_of(&row.lib);
    let target = dest.join(dirname);
    if target.is_dir() {
        // EXISTING tree: verify identity; never silently substitute.
        verify_pinned_clean(row, &target)?;
        return Ok("verified");
    }
    if offline {
        return Err(format!(
            "{} missing from the source cache in --offline mode",
            target.display()
        ));
    }
    // FETCH: clone the declared transport, check out the pinned revision
    // detached, then apply the same identity and cleanliness verifier as
    // the existing-tree path.
    let url = from.map_or_else(|| row.remote.clone(), |b| format!("{b}/{dirname}"));
    if url == "no-remote" {
        return Err(format!(
            "lock declares no remote for {} — re-lock on a host that has one",
            row.lib
        ));
    }
    let clone = std::process::Command::new("git")
        .args(["clone", "--no-checkout", "-c", "core.autocrlf=false", &url])
        .arg(&target)
        .output()
        .map_err(|e| format!("git clone failed to spawn: {e}"))?;
    if !clone.status.success() {
        return Err(format!(
            "clone of {url} failed: {}",
            String::from_utf8_lossy(&clone.stderr).trim()
        ));
    }
    git_out(&target, &["checkout", "--detach", &row.git_head]).map_err(|e| {
        format!(
            "locked revision {} unavailable from {url}: {e}",
            row.git_head
        )
    })?;
    verify_pinned_clean(row, &target)?;
    Ok("cloned")
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut root: Option<PathBuf> = None;
    let mut offline = false;
    let mut from: Option<String> = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--root" => match required_option_value("--root", it.next()) {
                Ok(value) => root = Some(PathBuf::from(value)),
                Err(error) => {
                    eprintln!("frankensim-bootstrap: {error}\n\n{}", usage());
                    return ExitCode::FAILURE;
                }
            },
            "--offline" => offline = true,
            "--from" => match required_option_value("--from", it.next()) {
                Ok(value) => from = Some(value.to_string()),
                Err(error) => {
                    eprintln!("frankensim-bootstrap: {error}\n\n{}", usage());
                    return ExitCode::FAILURE;
                }
            },
            "--help" | "-h" => {
                println!("{}", usage());
                return ExitCode::SUCCESS;
            }
            other => {
                eprintln!(
                    "frankensim-bootstrap: unknown flag {other:?}\n\n{}",
                    usage()
                );
                return ExitCode::FAILURE;
            }
        }
    }
    // Default root: the frankensim checkout this binary lives in, or cwd.
    let root = root.unwrap_or_else(|| {
        let cwd = std::env::current_dir().expect("cwd");
        if cwd.join("constellation.lock").is_file() {
            cwd
        } else {
            // Manifest-path invocations run from anywhere; walk up from
            // this source file's package to the checkout root.
            let tool_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
            tool_root
                .parent()
                .and_then(Path::parent)
                .map(Path::to_path_buf)
                .unwrap_or(cwd)
        }
    });
    let lock_path = root.join("constellation.lock");
    let lock_text = match std::fs::read_to_string(&lock_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!(
                "error: {} unreadable: {e} — the lock IS the input (pass --root <checkout>)",
                lock_path.display()
            );
            return ExitCode::FAILURE;
        }
    };
    let Some(dest) = root.parent().map(Path::to_path_buf) else {
        eprintln!("error: {} has no parent directory", root.display());
        return ExitCode::FAILURE;
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
        match bootstrap_one(row, &dest, offline, from.as_deref()) {
            Ok(state) => {
                println!(
                    "{{\"check\":\"constellation-bootstrap\",\"lib\":\"{}\",\"state\":\"{state}\",\"head\":\"{}\"}}",
                    row.lib, row.git_head
                );
                provenance.push(format!(
                    "{{\"lib\": \"{}\", \"git_head\": \"{}\", \"remote\": \"{}\", \"state\": \"{state}\"}}",
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
        "constellation bootstrap OK: {} libraries at their locked heads under {} (provenance: {})",
        rows.len(),
        dest.display(),
        prov_path.display()
    );
    ExitCode::SUCCESS
}
