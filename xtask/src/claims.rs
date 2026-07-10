//! Claim-state lint (bead 06yc): README prose must not drift from code.
//!
//! Public capability prose already has machine counterparts in the tree
//! (golden-hash constants, test function names, crate directories). This
//! check verifies the three cheapest, highest-yield couplings:
//!
//! 1. Every 16-hex-digit hash cited in README.md exists verbatim
//!    (underscore-insensitive, case-insensitive) somewhere under
//!    `crates/*/src` or `crates/*/tests` — a hash quoted in prose that no
//!    longer matches any recorded golden is stale evidence language.
//! 2. Every backticked `fs-<name>` crate reference in README.md exists as
//!    `crates/fs-<name>/` (wildcards like `fs-rep-*` and paths containing
//!    `::` or `_` are skipped — they are module prose, not crate names).
//! 3. Every backticked `*_hash` symbol in README.md exists as a
//!    `fn <name>` in some crate source or test — sentinel names in prose
//!    must be real tests.
//!
//! The deeper claim-state machinery (landed flags, no-claim rows, site
//! generation from evidence packages) belongs to huq.15.1; this lint is
//! the repo-level drift stop until that exists.

use crate::Violation;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// Normalize a hash token: strip `0x`, underscores, lowercase.
fn norm_hash(tok: &str) -> String {
    tok.trim_start_matches("0x")
        .chars()
        .filter(|c| *c != '_')
        .collect::<String>()
        .to_ascii_lowercase()
}

/// Extract candidate 64-bit hash literals (16 hex digits after
/// normalization) from a text.
fn hashes_in(text: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for (idx, _) in text.match_indices("0x") {
        let tail: String = text[idx + 2..]
            .chars()
            .take_while(|c| c.is_ascii_hexdigit() || *c == '_')
            .collect();
        let norm = norm_hash(&tail);
        if norm.len() == 16 {
            out.insert(norm);
        }
    }
    out
}

/// Backticked tokens in a markdown text.
fn backticked(text: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(open) = rest.find('`') {
        let after = &rest[open + 1..];
        let Some(close) = after.find('`') else { break };
        out.push(&after[..close]);
        rest = &after[close + 1..];
    }
    out
}

/// Walk all `.rs` files under `crates/*/{src,tests}`.
fn rust_files(root: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    let Ok(entries) = std::fs::read_dir(root.join("crates")) else {
        return files;
    };
    let mut stack: Vec<PathBuf> = entries
        .flatten()
        .flat_map(|e| [e.path().join("src"), e.path().join("tests")])
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
            } else if p.extension().is_some_and(|e| e == "rs") {
                files.push(p);
            }
        }
    }
    files
}

/// README claim-state lint: see module docs for the three rules.
pub fn check_claims(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();
    let readme = match std::fs::read_to_string(root.join("README.md")) {
        Ok(t) => t,
        Err(_) => {
            violations.push(Violation {
                check: "claim-state",
                crate_name: "<repo>".to_string(),
                detail: "README.md missing at workspace root".to_string(),
            });
            return violations;
        }
    };

    // Corpus: all code text (sources + tests) for hash and fn lookups.
    let mut code_hashes: BTreeSet<String> = BTreeSet::new();
    let mut code_text = String::new();
    for f in rust_files(root) {
        if let Ok(t) = std::fs::read_to_string(&f) {
            code_hashes.extend(hashes_in(&t));
            code_text.push_str(&t);
            code_text.push('\n');
        }
    }

    // Rule 1: cited hashes exist in code.
    for h in hashes_in(&readme) {
        if !code_hashes.contains(&h) {
            violations.push(Violation {
                check: "claim-state",
                crate_name: "README.md".to_string(),
                detail: format!(
                    "README cites hash 0x{h} but no crate source/test contains it — the prose \
                     is stale relative to the recorded goldens (re-check the sentinel it \
                     describes; golden bumps must update citing prose, bead 06yc)"
                ),
            });
        }
    }

    // Rules 2 and 3 over backticked tokens.
    for tok in backticked(&readme) {
        // Rule 2: crate references.
        if let Some(name) = tok.strip_prefix("fs-") {
            let clean = name
                .chars()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-');
            if clean && !name.is_empty() && !root.join("crates").join(tok).is_dir() {
                violations.push(Violation {
                    check: "claim-state",
                    crate_name: "README.md".to_string(),
                    detail: format!(
                        "README references crate `{tok}` but crates/{tok}/ does not exist \
                         (renamed or removed crate leaves stale capability prose, bead 06yc)"
                    ),
                });
            }
        }
        // Rule 3: sentinel test symbols.
        if tok.ends_with("_hash")
            && tok
                .chars()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
            && !code_text.contains(&format!("fn {tok}"))
        {
            violations.push(Violation {
                check: "claim-state",
                crate_name: "README.md".to_string(),
                detail: format!(
                    "README names sentinel `{tok}` but no `fn {tok}` exists in any crate \
                     source/test (bead 06yc)"
                ),
            });
        }
    }
    violations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_and_backtick_extraction() {
        let hs = hashes_in("golden `0xeef1_0550_7daf_c0d5` and 0xDEAD (too short)");
        assert!(hs.contains("eef105507dafc0d5"));
        assert_eq!(hs.len(), 1);
        assert_eq!(backticked("a `b` c `d-e`"), vec!["b", "d-e"]);
    }

    #[test]
    fn claims_check_end_to_end_on_fixture_tree() {
        let base = std::env::temp_dir().join(format!("fsim-claims-test-{}", std::process::id()));
        let mk = |rel: &str, content: &str| {
            let p = base.join(rel);
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            std::fs::write(p, content).unwrap();
        };
        mk(
            "crates/fs-real/src/lib.rs",
            "pub const G: u64 = 0x1111_2222_3333_4444;\n",
        );
        mk(
            "crates/fs-real/tests/battery.rs",
            "fn real_golden_hash() {}\n",
        );
        // Seeded drift: stale hash, missing crate, missing sentinel fn.
        mk(
            "README.md",
            concat!(
                "Good: `fs-real` golden `0x1111_2222_3333_4444` via `real_golden_hash`.\n",
                "Stale hash 0xaaaa_bbbb_cccc_dddd.\n",
                "Gone crate `fs-vanished`.\n",
                "Gone sentinel `ghost_golden_hash`.\n",
            ),
        );
        let v = check_claims(&base);
        assert_eq!(v.len(), 3, "exactly the three seeded drifts: {v:?}");
        assert!(v.iter().any(|x| x.detail.contains("aaaabbbbccccdddd")));
        assert!(v.iter().any(|x| x.detail.contains("fs-vanished")));
        assert!(v.iter().any(|x| x.detail.contains("ghost_golden_hash")));
        let _ = std::fs::remove_dir_all(&base);
    }
}
