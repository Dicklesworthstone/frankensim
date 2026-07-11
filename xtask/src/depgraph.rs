//! `depgraph-receipt` — canonical resolved dependency+feature receipt for
//! the GEMM build identity (bead fz2.6).
//!
//! The fs-la build fingerprint binds current-crate `CARGO_FEATURE_*` plus
//! sources/manifests/Cargo.lock, but Cargo feature unification can compile
//! path and registry dependencies (notably asupersync and its aes-gcm/
//! rand_core closure) under a different active feature set without changing
//! any of those inputs. This command derives a canonical receipt of the
//! resolved normal-dependency closure of `fs-la` — package identity plus the
//! unified feature set per package — from `cargo tree`, OUTSIDE any build
//! script (build.rs must never invoke Cargo recursively). Build tooling
//! exports the receipt as `FRANKENSIM_DEPGRAPH_RECEIPT`; fs-la's build.rs
//! binds the bytes into the fingerprint and fails closed without either a
//! receipt or the explicit workspace salt.
//!
//! Precision note (documented contract): `cargo tree` resolves features like
//! `cargo metadata` — dev-dependencies of the selected packages participate
//! in unification, so the receipt is a conservative over-approximation of a
//! plain `cargo build` and exact for test builds of the same selection. The
//! invocation-exact `-Z unit-graph` no longer exists in Cargo; this is the
//! strongest evidence derivable from stable tooling, and two builds whose
//! normal-dependency unification differs always receive different receipts.

use std::path::Path;
use std::process::Command;

const RECEIPT_SCHEMA: &str = "fs-la-depgraph-receipt-v1";
const CLOSURE_ROOT: &str = "fs-la ";

/// One resolved package in the fs-la closure: display id and unified,
/// comma-joined feature list exactly as `cargo tree --format "{p}|{f}"`
/// reports them.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ResolvedPackage {
    id: String,
    features: String,
}

/// Run the receipt command. `args` are the Cargo selection flags of the
/// build being fingerprinted (e.g. `-p fs-roofline`); they are passed to
/// `cargo tree` verbatim so the resolve mirrors the build's selection.
///
/// Prints the canonical receipt JSON on success. With `--verify`, instead
/// compares against the `FRANKENSIM_DEPGRAPH_RECEIPT` environment variable
/// and fails on mismatch.
pub fn cmd_depgraph_receipt(root: &Path, raw_args: &[String]) -> Result<(), String> {
    let mut verify = false;
    let mut selection = Vec::new();
    let mut past_separator = false;
    for arg in raw_args {
        match (past_separator, arg.as_str()) {
            (false, "--verify") => verify = true,
            (false, "--") => past_separator = true,
            _ => selection.push(arg.clone()),
        }
    }
    let receipt = derive_receipt(root, &selection)?;
    if verify {
        let supplied = std::env::var("FRANKENSIM_DEPGRAPH_RECEIPT").map_err(|_| {
            "verify mode requires FRANKENSIM_DEPGRAPH_RECEIPT in the environment".to_string()
        })?;
        if supplied == receipt {
            println!("depgraph receipt verified: {} bytes", receipt.len());
            Ok(())
        } else {
            Err(format!(
                "depgraph receipt mismatch: environment carries {} bytes, recomputation yields {} bytes; \
                 the build graph changed under the receipt",
                supplied.len(),
                receipt.len()
            ))
        }
    } else {
        println!("{receipt}");
        Ok(())
    }
}

/// Derive the canonical receipt for the given Cargo selection flags.
fn derive_receipt(root: &Path, selection: &[String]) -> Result<String, String> {
    let mut command = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()));
    command
        .current_dir(root)
        .args([
            "tree",
            "-e",
            "normal",
            "--no-dedupe",
            "--format",
            "{p}|{f}",
            "--prefix",
            "depth",
        ])
        .args(selection);
    let output = command
        .output()
        .map_err(|error| format!("cannot execute cargo tree: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "cargo tree failed for selection {selection:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    let text = String::from_utf8(output.stdout)
        .map_err(|error| format!("cargo tree emitted non-UTF-8 output: {error}"))?;
    let packages = fs_la_closure(&text)?;
    if packages.is_empty() {
        return Err(format!(
            "selection {selection:?} does not build fs-la; no depgraph receipt applies"
        ));
    }
    Ok(canonical_receipt_json(&packages))
}

/// Extract the union of every `fs-la`-rooted subtree from `--prefix depth`
/// output. Refuses ambiguous evidence: one package resolved with two
/// different feature sets inside the closure.
fn fs_la_closure(tree: &str) -> Result<Vec<ResolvedPackage>, String> {
    let mut closure: Vec<ResolvedPackage> = Vec::new();
    let mut in_subtree_above: Option<usize> = None;
    for line in tree.lines() {
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let digits = line.bytes().take_while(u8::is_ascii_digit).count();
        if digits == 0 {
            return Err(format!("cargo tree line has no depth prefix: {line:?}"));
        }
        let depth: usize = line[..digits]
            .parse()
            .map_err(|error| format!("bad depth prefix in {line:?}: {error}"))?;
        let rest = &line[digits..];
        let (id, features) = rest
            .split_once('|')
            .ok_or_else(|| format!("cargo tree line missing feature separator: {line:?}"))?;
        // `cargo tree` suffixes repeat visits with a ` (*)` de-duplication
        // marker; it is display state, not resolution evidence.
        let features = features.strip_suffix(" (*)").unwrap_or(features);
        if let Some(root_depth) = in_subtree_above {
            if depth <= root_depth {
                in_subtree_above = None;
            }
        }
        let is_root = id.starts_with(CLOSURE_ROOT);
        if is_root && in_subtree_above.is_none() {
            in_subtree_above = Some(depth);
        }
        if in_subtree_above.is_some() {
            let entry = ResolvedPackage {
                id: id.to_string(),
                features: features.to_string(),
            };
            if let Some(existing) = closure.iter().find(|candidate| candidate.id == entry.id) {
                if existing.features != entry.features {
                    return Err(format!(
                        "ambiguous depgraph evidence: {} resolved with features {:?} and {:?} \
                         in one selection; refusing to mint a receipt",
                        entry.id, existing.features, entry.features
                    ));
                }
            } else {
                closure.push(entry);
            }
        }
    }
    closure.sort();
    Ok(closure)
}

/// Canonical, byte-stable JSON: schema tag plus id-sorted packages with
/// their unified feature CSV split into a sorted array.
fn canonical_receipt_json(packages: &[ResolvedPackage]) -> String {
    let body = packages
        .iter()
        .map(|package| {
            let mut features: Vec<&str> = package
                .features
                .split(',')
                .filter(|feature| !feature.is_empty())
                .collect();
            features.sort_unstable();
            let features = features
                .iter()
                .map(|feature| format!("\"{feature}\""))
                .collect::<Vec<_>>()
                .join(",");
            format!("{{\"id\":\"{}\",\"features\":[{features}]}}", package.id)
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{{\"schema\":\"{RECEIPT_SCHEMA}\",\"packages\":[{body}]}}")
}

#[cfg(test)]
mod tests {
    use super::*;

    const BASE_TREE: &str = "\
0fs-roofline v0.0.1 (/w/crates/fs-roofline)|
1fs-la v0.0.1 (/w/crates/fs-la)|
2fs-exec v0.0.1 (/w/crates/fs-exec)|
3asupersync v0.3.5 (/w/../asupersync)|default,proc-macros
4aes-gcm v0.10.3|aes,alloc,default
2fs-alloc v0.0.1 (/w/crates/fs-alloc)|
1fs-ledger v0.0.1 (/w/crates/fs-ledger)|
";

    #[test]
    fn closure_is_the_fs_la_subtree_only() {
        let closure = fs_la_closure(BASE_TREE).expect("closure");
        let ids: Vec<&str> = closure
            .iter()
            .map(|package| package.id.as_str())
            .collect();
        assert!(ids.iter().any(|id| id.starts_with("fs-la ")));
        assert!(ids.iter().any(|id| id.starts_with("aes-gcm ")));
        assert!(
            !ids.iter().any(|id| id.starts_with("fs-ledger ")),
            "sibling outside the fs-la subtree must not enter the receipt"
        );
    }

    #[test]
    fn different_unified_dependency_features_change_the_receipt() {
        let drifted = BASE_TREE.replace(
            "4aes-gcm v0.10.3|aes,alloc,default",
            "4aes-gcm v0.10.3|aes,alloc,default,getrandom",
        );
        let left = canonical_receipt_json(&fs_la_closure(BASE_TREE).expect("base"));
        let right = canonical_receipt_json(&fs_la_closure(&drifted).expect("drifted"));
        assert_ne!(
            left, right,
            "a dependency feature invisible to CARGO_FEATURE_* must still move the receipt"
        );
    }

    #[test]
    fn identical_graphs_replay_byte_identically() {
        let left = canonical_receipt_json(&fs_la_closure(BASE_TREE).expect("first"));
        let right = canonical_receipt_json(&fs_la_closure(BASE_TREE).expect("second"));
        assert_eq!(left, right);
    }

    #[test]
    fn feature_order_is_canonicalized() {
        let reordered = BASE_TREE.replace(
            "4aes-gcm v0.10.3|aes,alloc,default",
            "4aes-gcm v0.10.3|default,aes,alloc",
        );
        let left = canonical_receipt_json(&fs_la_closure(BASE_TREE).expect("base"));
        let right = canonical_receipt_json(&fs_la_closure(&reordered).expect("reordered"));
        assert_eq!(left, right, "feature listing order must not move the receipt");
    }

    #[test]
    fn conflicting_feature_sets_for_one_package_refuse_a_receipt() {
        let ambiguous = format!(
            "{BASE_TREE}1fs-la v0.0.1 (/w/crates/fs-la)|\n2aes-gcm v0.10.3|aes,alloc\n"
        );
        let error = fs_la_closure(&ambiguous).expect_err("ambiguity must refuse");
        assert!(error.contains("ambiguous depgraph evidence"), "{error}");
    }

    #[test]
    fn selections_without_fs_la_are_refused_upstream() {
        let closure = fs_la_closure("0fs-qty v0.0.1 (/w/crates/fs-qty)|\n").expect("closure");
        assert!(closure.is_empty());
    }
}
