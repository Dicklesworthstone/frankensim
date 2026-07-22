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
//! Rule 4 (huq.18) derives README inventory counts from the tree, and rule 5
//! (f85xj.2.1) keeps the claim-integrity defect class defined and its label
//! taxonomy documented — see the section further down for why that definition
//! is load-bearing rather than decorative.
//!
//! The deeper claim-state machinery (landed flags, no-claim rows, site
//! generation from evidence packages) belongs to huq.15.1; this lint is
//! the repo-level drift stop until that exists.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const DOC_FACTS_BEGIN: &str = "<!-- BEGIN GENERATED FRANKENSIM DOC FACTS -->";
const DOC_FACTS_END: &str = "<!-- END GENERATED FRANKENSIM DOC FACTS -->";
const DOC_FACTS_CHECK: &str = "doc-facts";
const DOC_INVENTORY_FILE: &str = "doc-facts-inventory.json";
const DOC_INVENTORY_SCHEMA: &str = "frankensim-doc-facts-inventory-v1";

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

/// The number immediately preceding `pat` on `line` (digits touching the
/// pattern), if any.
fn count_before(line: &str, pat: &str) -> Option<usize> {
    let pos = line.find(pat)?;
    let digits: String = line[..pos]
        .chars()
        .rev()
        .take_while(char::is_ascii_digit)
        .collect::<String>()
        .chars()
        .rev()
        .collect();
    digits.parse::<usize>().ok()
}

#[cfg(test)]
fn workspace_fs_member_count(manifest: &str) -> Option<usize> {
    workspace_fs_members(manifest).map(|members| members.len())
}

fn workspace_fs_members(manifest: &str) -> Option<BTreeSet<String>> {
    let mut lines = manifest.lines();
    lines.find(|line| line.trim() == "members = [")?;
    let mut members = BTreeSet::new();
    for line in lines {
        let entry = line.trim();
        if entry == "]" {
            return Some(members);
        }
        let entry = entry.strip_suffix(',').unwrap_or(entry).trim();
        let entry = entry.strip_prefix('"')?.strip_suffix('"')?;
        if entry.starts_with("crates/fs-") {
            members.insert(entry.to_string());
        }
    }
    None
}

fn git_tracked_files(root: &Path, pathspec: &str) -> Result<Option<Vec<String>>, String> {
    let repository = std::process::Command::new("git")
        .args(["-C"])
        .arg(root)
        .args(["rev-parse", "--is-inside-work-tree"])
        .output()
        .map_err(|error| format!("cannot inspect Git worktree authority: {error}"))?;
    if !repository.status.success() {
        let stderr = String::from_utf8_lossy(&repository.stderr);
        if stderr.contains("not a git repository") {
            return Ok(None);
        }
        return Err(format!(
            "git rev-parse failed while inspecting documentation inventory authority: {}",
            stderr.trim()
        ));
    }
    if String::from_utf8_lossy(&repository.stdout).trim() != "true" {
        return Ok(None);
    }

    let output = std::process::Command::new("git")
        .args(["-C"])
        .arg(root)
        .args(["ls-files", "-z", "--", pathspec])
        .output()
        .map_err(|error| format!("cannot execute git ls-files for {pathspec:?}: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "git ls-files failed for {pathspec:?}: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    let mut seen = BTreeSet::new();
    let mut files = Vec::new();
    for raw in output.stdout.split(|byte| *byte == 0) {
        if raw.is_empty() {
            continue;
        }
        let relative = std::str::from_utf8(raw).map_err(|error| {
            format!("git inventory for {pathspec:?} contains a non-UTF-8 path: {error}")
        })?;
        if relative.starts_with('/')
            || Path::new(relative)
                .components()
                .any(|component| matches!(component, std::path::Component::ParentDir))
        {
            return Err(format!(
                "git inventory for {pathspec:?} contains non-portable path {relative:?}"
            ));
        }
        if !seen.insert(relative.to_string()) {
            return Err(format!(
                "git inventory for {pathspec:?} repeats path {relative:?}"
            ));
        }
        files.push(relative.to_string());
    }
    Ok(Some(files))
}

fn json_object(value: &JsonValue) -> Option<&BTreeMap<String, JsonValue>> {
    match value {
        JsonValue::Object(map) => Some(map),
        _ => None,
    }
}

fn json_array(value: &JsonValue) -> Option<&[JsonValue]> {
    match value {
        JsonValue::Array(items) => Some(items),
        _ => None,
    }
}

fn json_text(value: &JsonValue) -> Option<&str> {
    match value {
        JsonValue::String(text) => Some(text),
        _ => None,
    }
}

fn valid_inventory_path(group: &str, relative: &str) -> bool {
    if relative.contains('\\') {
        return false;
    }
    let parts = relative.split('/').collect::<Vec<_>>();
    let crate_path = parts.first() == Some(&"crates")
        && parts
            .get(1)
            .is_some_and(|name| name.starts_with("fs-") && name.len() > 3);
    match group {
        "manifests" => crate_path && parts.len() == 3 && parts[2] == "Cargo.toml",
        "contracts" => crate_path && parts.len() == 3 && parts[2] == "CONTRACT.md",
        "integration_tests" => {
            crate_path
                && parts.len() >= 4
                && parts[2] == "tests"
                && parts.last().is_some_and(|name| name.ends_with(".rs"))
                && parts.last().is_some_and(|name| name.len() > 3)
        }
        _ => false,
    }
}

fn inventory_paths(
    root: &Path,
    object: &BTreeMap<String, JsonValue>,
    group: &str,
) -> Result<Vec<String>, String> {
    let values = object
        .get(group)
        .and_then(json_array)
        .ok_or_else(|| format!("{DOC_INVENTORY_FILE} has no array field {group:?}"))?;
    let mut paths = Vec::with_capacity(values.len());
    for (index, value) in values.iter().enumerate() {
        let relative = json_text(value)
            .ok_or_else(|| format!("{DOC_INVENTORY_FILE} {group}[{index}] is not a string"))?;
        if !valid_inventory_path(group, relative) {
            return Err(format!(
                "{DOC_INVENTORY_FILE} {group}[{index}] has invalid path {relative:?}"
            ));
        }
        let path = Path::new(relative);
        if path
            .components()
            .any(|component| !matches!(component, std::path::Component::Normal(_)))
        {
            return Err(format!(
                "{DOC_INVENTORY_FILE} {group}[{index}] has non-portable path {relative:?}"
            ));
        }
        if !root.join(path).is_file() {
            return Err(format!(
                "{DOC_INVENTORY_FILE} {group}[{index}] names missing or non-file path {relative:?}"
            ));
        }
        if paths
            .last()
            .is_some_and(|previous| previous >= &relative.to_string())
        {
            return Err(format!(
                "{DOC_INVENTORY_FILE} {group} must be strictly sorted with no duplicates; {relative:?} is out of order"
            ));
        }
        paths.push(relative.to_string());
    }
    Ok(paths)
}

struct TrackedInventory {
    manifests: Vec<String>,
    contracts: Vec<String>,
    integration_tests: Vec<String>,
    git_index_verified: bool,
}

impl TrackedInventory {
    fn read(root: &Path) -> Result<Self, String> {
        let source = std::fs::read_to_string(root.join(DOC_INVENTORY_FILE))
            .map_err(|error| format!("cannot read {DOC_INVENTORY_FILE}: {error}"))?;
        let parsed = JsonParser::new(&source)
            .finish()
            .map_err(|error| format!("cannot parse {DOC_INVENTORY_FILE}: {error}"))?;
        let object = json_object(&parsed)
            .ok_or_else(|| format!("{DOC_INVENTORY_FILE} is not a JSON object"))?;
        let expected_fields =
            BTreeSet::from(["schema", "manifests", "contracts", "integration_tests"]);
        let actual_fields = object.keys().map(String::as_str).collect::<BTreeSet<_>>();
        if actual_fields != expected_fields {
            return Err(format!(
                "{DOC_INVENTORY_FILE} fields are {actual_fields:?}, expected {expected_fields:?}"
            ));
        }
        match object.get("schema").and_then(json_text) {
            Some(DOC_INVENTORY_SCHEMA) => {}
            Some(schema) => {
                return Err(format!(
                    "{DOC_INVENTORY_FILE} schema is {schema:?}, expected {DOC_INVENTORY_SCHEMA:?}"
                ));
            }
            None => {
                return Err(format!("{DOC_INVENTORY_FILE} has no string schema field"));
            }
        }

        let manifests = inventory_paths(root, object, "manifests")?;
        let contracts = inventory_paths(root, object, "contracts")?;
        let integration_tests = inventory_paths(root, object, "integration_tests")?;
        let groups = [
            ("manifests", "crates/fs-*/Cargo.toml", &manifests),
            ("contracts", "crates/fs-*/CONTRACT.md", &contracts),
            (
                "integration_tests",
                "crates/fs-*/tests/*.rs",
                &integration_tests,
            ),
        ];
        let mut git_index_verified = true;
        for (group, pathspec, recorded) in groups {
            match git_tracked_files(root, pathspec)? {
                Some(actual) if &actual == recorded => {}
                Some(actual) => {
                    let first_difference = actual
                        .iter()
                        .zip(recorded.iter())
                        .position(|(left, right)| left != right)
                        .unwrap_or_else(|| actual.len().min(recorded.len()));
                    return Err(format!(
                        "{DOC_INVENTORY_FILE} {group} is stale against the Git index at entry {first_difference}: recorded={} git={}",
                        recorded.len(),
                        actual.len()
                    ));
                }
                None => git_index_verified = false,
            }
        }
        Ok(Self {
            manifests,
            contracts,
            integration_tests,
            git_index_verified,
        })
    }
}

#[derive(Debug, PartialEq, Eq)]
struct DocFacts {
    native_workspace_crates: usize,
    standalone_crates: Vec<String>,
    crate_directories: usize,
    contracts: usize,
    tracked_integration_tests: usize,
    layers: BTreeMap<&'static str, usize>,
}

impl DocFacts {
    fn derive(root: &Path) -> Result<(Self, bool), String> {
        let workspace_source = std::fs::read_to_string(root.join("Cargo.toml"))
            .map_err(|error| format!("cannot read Cargo.toml: {error}"))?;
        let native_members = workspace_fs_members(&workspace_source).ok_or_else(|| {
            "cannot derive native fs-* workspace members from [workspace].members".to_string()
        })?;

        let inventory = TrackedInventory::read(root)?;
        let manifest_paths = inventory
            .manifests
            .iter()
            .map(|relative| root.join(relative))
            .collect::<Vec<_>>();
        if manifest_paths.is_empty() {
            return Err("no fs-* crate manifests were found".to_string());
        }

        let mut standalone_crates = Vec::new();
        let mut layers = BTreeMap::from([
            ("UTIL", 0usize),
            ("L0", 0),
            ("L1", 0),
            ("L2", 0),
            ("L3", 0),
            ("L4", 0),
            ("L5", 0),
            ("L6", 0),
        ]);
        for path in &manifest_paths {
            let source = std::fs::read_to_string(path)
                .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
            let manifest = crate::parse_manifest(path, &source)?;
            let layer = manifest.layer.name();
            let Some(count) = layers.get_mut(layer) else {
                return Err(format!(
                    "{} declares non-product layer {layer}; fs-* documentation inventory covers UTIL and L0-L6",
                    path.display()
                ));
            };
            *count = count
                .checked_add(1)
                .ok_or_else(|| "layer inventory count overflow".to_string())?;

            let relative = path
                .strip_prefix(root)
                .map_err(|_| format!("{} is outside workspace root", path.display()))?
                .to_string_lossy()
                .replace('\\', "/");
            let member = relative
                .strip_suffix("/Cargo.toml")
                .ok_or_else(|| format!("unexpected manifest path {relative}"))?;
            if !native_members.contains(member) {
                standalone_crates.push(manifest.name);
            }
        }
        standalone_crates.sort();

        let crate_directories = manifest_paths.len();
        if native_members.len() + standalone_crates.len() != crate_directories {
            return Err(format!(
                "native ({}) plus standalone ({}) fs-* crates do not cover all {crate_directories} tracked manifests",
                native_members.len(),
                standalone_crates.len()
            ));
        }

        let contracts = inventory.contracts.len();
        let tracked_integration_tests = inventory.integration_tests.len();

        Ok((
            Self {
                native_workspace_crates: native_members.len(),
                standalone_crates,
                crate_directories,
                contracts,
                tracked_integration_tests,
                layers,
            },
            inventory.git_index_verified,
        ))
    }

    fn render(&self) -> String {
        let standalone = if self.standalone_crates.is_empty() {
            "0".to_string()
        } else {
            format!(
                "{} ({})",
                self.standalone_crates.len(),
                self.standalone_crates
                    .iter()
                    .map(|name| format!("`{name}`"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        };
        let layer_inventory = self
            .layers_ordered()
            .map(|(layer, count)| format!("`{layer}={count}`"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "{DOC_FACTS_BEGIN}\n\
| Derived repository fact | Value |\n\
|-------------------------|-------|\n\
| Native workspace `fs-*` crates | {} |\n\
| Standalone `fs-*` workspaces | {standalone} |\n\
| Tracked `fs-*` crate directories | {} |\n\
| Tracked `CONTRACT.md` files | {} of {} |\n\
| Tracked crate integration-test files | {} |\n\
| `fs-*` layer inventory | {layer_inventory} |\n\
{DOC_FACTS_END}",
            self.native_workspace_crates,
            self.crate_directories,
            self.contracts,
            self.crate_directories,
            self.tracked_integration_tests,
        )
    }

    fn layers_ordered(&self) -> impl Iterator<Item = (&'static str, usize)> + '_ {
        ["UTIL", "L0", "L1", "L2", "L3", "L4", "L5", "L6"]
            .into_iter()
            .map(|layer| (layer, self.layers.get(layer).copied().unwrap_or(0)))
    }
}

fn marked_block<'a>(source: &'a str, begin: &str, end: &str) -> Result<&'a str, String> {
    let starts: Vec<usize> = source
        .match_indices(begin)
        .map(|(index, _)| index)
        .collect();
    let ends: Vec<usize> = source.match_indices(end).map(|(index, _)| index).collect();
    if starts.len() != 1 || ends.len() != 1 {
        return Err(format!(
            "expected exactly one {begin:?} and one {end:?} marker, found {} and {}",
            starts.len(),
            ends.len()
        ));
    }
    let start = starts[0];
    let finish = ends[0]
        .checked_add(end.len())
        .ok_or_else(|| "generated-block end offset overflow".to_string())?;
    if ends[0] <= start {
        return Err(format!("{end:?} appears before {begin:?}"));
    }
    Ok(&source[start..finish])
}

fn check_doc_facts(readme: &str, facts: &DocFacts) -> Vec<Violation> {
    let expected = facts.render();
    match marked_block(readme, DOC_FACTS_BEGIN, DOC_FACTS_END) {
        Ok(actual) if actual == expected => Vec::new(),
        Ok(_) => vec![Violation {
            check: DOC_FACTS_CHECK,
            crate_name: "README.md".to_string(),
            detail: format!(
                "README generated repository-facts block is stale; replace it with the exact tree-derived block:\n{expected}"
            ),
        }],
        Err(error) => vec![Violation {
            check: DOC_FACTS_CHECK,
            crate_name: "README.md".to_string(),
            detail: format!("README generated repository-facts block is malformed: {error}"),
        }],
    }
}

pub struct DocsReport {
    pub violations: Vec<Violation>,
    pub decisions: Vec<PolicyNote>,
}

fn doc_note(source: &str, detail: String) -> PolicyNote {
    PolicyNote {
        check: DOC_FACTS_CHECK,
        crate_name: source.to_string(),
        verdict: "verified",
        detail,
    }
}

fn check_docs_with_facts(readme: &str, facts: &DocFacts) -> DocsReport {
    let mut violations = check_inventory_counts(readme, facts);
    let block_violations = check_doc_facts(readme, facts);
    let block_clean = block_violations.is_empty();
    violations.extend(block_violations);

    let mut decisions = vec![
        doc_note(
            "Cargo.toml",
            format!(
                "native_workspace_fs_members={} source=[workspace].members",
                facts.native_workspace_crates
            ),
        ),
        doc_note(
            "doc-facts-inventory.json:manifests",
            format!(
                "tracked_fs_crate_manifests={} standalone={} source=portable-checked-git-index",
                facts.crate_directories,
                facts.standalone_crates.join(",")
            ),
        ),
        doc_note(
            "doc-facts-inventory.json:contracts",
            format!(
                "tracked_contracts={} crate_manifests={}",
                facts.contracts, facts.crate_directories
            ),
        ),
        doc_note(
            "doc-facts-inventory.json:integration_tests",
            format!(
                "tracked_integration_test_files={}",
                facts.tracked_integration_tests
            ),
        ),
        doc_note(
            "crates/*/Cargo.toml:[package.metadata.frankensim].layer",
            format!(
                "declared_layer_inventory={} (metadata inventory, not dependency-validity proof)",
                facts
                    .layers_ordered()
                    .map(|(layer, count)| format!("{layer}={count}"))
                    .collect::<Vec<_>>()
                    .join(" ")
            ),
        ),
    ];
    if block_clean {
        decisions.push(doc_note(
            "README.md",
            "generated repository-facts block exactly matches its tracked sources".to_string(),
        ));
    }
    DocsReport {
        violations,
        decisions,
    }
}

pub fn check_docs(root: &Path) -> DocsReport {
    let readme = match std::fs::read_to_string(root.join("README.md")) {
        Ok(readme) => readme,
        Err(error) => {
            return DocsReport {
                violations: vec![Violation {
                    check: DOC_FACTS_CHECK,
                    crate_name: "README.md".to_string(),
                    detail: format!("cannot read README.md: {error}"),
                }],
                decisions: Vec::new(),
            };
        }
    };
    match DocFacts::derive(root) {
        Ok((facts, git_index_verified)) => {
            let mut report = check_docs_with_facts(&readme, &facts);
            report.decisions.insert(
                0,
                doc_note(
                    DOC_INVENTORY_FILE,
                    if git_index_verified {
                        "portable tracked-file registry exactly matches git ls-files for all documentation fact inputs"
                            .to_string()
                    } else {
                        "source-snapshot mode: .git is unavailable; portable tracked-file registry is schema-valid, strictly sorted, duplicate-free, and every recorded input exists (a Git worktree additionally verifies exact index equality)"
                            .to_string()
                    },
                ),
            );
            report
        }
        Err(error) => DocsReport {
            violations: vec![Violation {
                check: DOC_FACTS_CHECK,
                crate_name: "README.md".to_string(),
                detail: format!("cannot derive generated documentation facts: {error}"),
            }],
            decisions: Vec::new(),
        },
    }
}

/// huq.18: README inventory counts (crate/contract/test-file numbers in
/// the badges and the What-Exists table) must equal the tree's actual
/// counts — counts are DERIVED, never hand-promoted, so drift turns the
/// gate red instead of aging silently.
fn check_inventory_counts(readme: &str, facts: &DocFacts) -> Vec<Violation> {
    let mut violations = Vec::new();
    let crate_count = facts.crate_directories;
    let contract_count = facts.contracts;
    let test_file_count = facts.tracked_integration_tests;
    let workspace_crate_count = facts.native_workspace_crates;
    let checks = [
        (
            "%20native%20fs--%2A%20crates",
            workspace_crate_count,
            "native fs-* workspace crates (badge)",
        ),
        (
            "%20fs--%2A%20crates",
            workspace_crate_count,
            "fs-* crates (badge)",
        ),
        (
            " native `fs-*` workspace crates",
            workspace_crate_count,
            "native fs-* workspace crates",
        ),
        (
            " `fs-*` crate directories",
            crate_count,
            "fs-* crate directories",
        ),
        (" fs-* crates", crate_count, "fs-* crates (layout)"),
        (" `CONTRACT.md` files", contract_count, "CONTRACT.md files"),
        (" crate test files", test_file_count, "crate test files"),
        (
            " crate-level conformance",
            test_file_count,
            "crate test files (What Exists table)",
        ),
        (
            "%20crate%20test%20files",
            test_file_count,
            "crate test files (badge)",
        ),
        (
            "%20tracked%20integration%20test%20files",
            test_file_count,
            "tracked integration-test files (badge)",
        ),
        (
            " tracked Rust files under crate `tests/` directories",
            test_file_count,
            "tracked integration-test files",
        ),
    ];
    for line in readme.lines() {
        for (pat, actual, what) in checks {
            if let Some(claimed) = count_before(line, pat)
                && claimed != actual
            {
                violations.push(Violation {
                    check: DOC_FACTS_CHECK,
                    crate_name: "README.md".to_string(),
                    detail: format!(
                        "README claims {claimed} {what} but the tree has {actual} — counts \
                         are derived, never hand-promoted (bead huq.18)"
                    ),
                });
            }
        }
        // Contracts badge: `contracts-<n>%20of%20<m>%20crates`.
        if let Some(at) = line.find("badge/contracts-") {
            let tail = &line[at + "badge/contracts-".len()..];
            let n: String = tail.chars().take_while(char::is_ascii_digit).collect();
            let m = tail
                .find("%20of%20")
                .map(|p| &tail[p + "%20of%20".len()..])
                .map(|t| {
                    t.chars()
                        .take_while(char::is_ascii_digit)
                        .collect::<String>()
                });
            if let (Ok(n), Some(Ok(m))) = (n.parse::<usize>(), m.map(|s| s.parse::<usize>()))
                && (n != contract_count || m != crate_count)
            {
                violations.push(Violation {
                    check: DOC_FACTS_CHECK,
                    crate_name: "README.md".to_string(),
                    detail: format!(
                        "README contracts badge says {n} of {m} but the tree has \
                         {contract_count} CONTRACT.md files across {crate_count} crates \
                         (bead huq.18)"
                    ),
                });
            }
        }
        if line.contains("| Contracts |") || line.contains("`CONTRACT.md` files for") {
            let numbers: Vec<usize> = line
                .split(|ch: char| !ch.is_ascii_digit())
                .filter(|token| !token.is_empty())
                .filter_map(|token| token.parse().ok())
                .collect();
            if numbers.len() >= 2 && numbers[numbers.len() - 2..] != [contract_count, crate_count] {
                violations.push(Violation {
                    check: DOC_FACTS_CHECK,
                    crate_name: "README.md".to_string(),
                    detail: format!(
                        "README contract inventory says {} of {} but the tree has \
                         {contract_count} CONTRACT.md files across {crate_count} crates \
                         (bead huq.18)",
                        numbers[numbers.len() - 2],
                        numbers[numbers.len() - 1]
                    ),
                });
            }
        }
    }
    violations
}

// ---------------------------------------------------------------------------
// Claim-integrity defect class (bead frankensim-extreal-program-f85xj.2.1).
//
// The class definition is the input the E02 sweep and promotion gate consume
// verbatim, so it must not silently disappear, lose a decision-rule section,
// or drift out of agreement with the label taxonomy the gate queries. Code is
// the single source of truth for the canonical severity labels; the definition
// doc and the CONVENTIONS taxonomy must both name exactly these.
//
// This lint proves the definition is present and structurally intact. It does
// not, and cannot, judge whether an audit was performed honestly — that is
// what the sweep's recorded verdicts and the gate drills are for. Claiming
// otherwise here would itself be a claim-integrity defect.
// ---------------------------------------------------------------------------

/// Canonical severity labels. The gate and the inventory script accept exactly
/// these; adding a severity means changing this array and both documents.
pub const CLAIM_INTEGRITY_SEVERITY_LABELS: [&str; 3] = [
    "severity:default-path",
    "severity:gated",
    "severity:doc-only",
];

/// The mandatory class-membership label; `br list -l <label>` is the inventory.
pub const CLAIM_INTEGRITY_LABEL: &str = "claim-integrity";

const CLAIM_INTEGRITY_DOC: &str = "docs/CLAIM_INTEGRITY.md";
const CLAIM_INTEGRITY_CONVENTIONS: &str = "docs/CONVENTIONS.md";
const CLAIM_INTEGRITY_INVENTORY_SCRIPT: &str = "scripts/ci/claim_integrity_inventory.sh";

/// Sections the definition must keep. Each one is consumed by a downstream
/// bead: decision rules and audit method by the `.2.2` sweep, severity rules
/// and label taxonomy by the `.2.3` gate, known instances by both as the
/// known-answer set.
const CLAIM_INTEGRITY_REQUIRED_SECTIONS: [&str; 6] = [
    "## Definition",
    "## Decision rules",
    "## Severity rules",
    "## Label taxonomy",
    "## Audit method",
    "## Known instances",
];

/// The CONVENTIONS taxonomy section heading that must point agents at the
/// definition.
const CLAIM_INTEGRITY_CONVENTIONS_SECTION: &str = "## Claim-integrity defect class";

fn claim_integrity_violation(file: &str, detail: String) -> Violation {
    Violation {
        check: "claim-state",
        crate_name: file.to_string(),
        detail,
    }
}

/// Lint the claim-integrity definition and its taxonomy (bead f85xj.2.1).
fn check_claim_integrity_docs(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();

    let Ok(definition) = std::fs::read_to_string(root.join(CLAIM_INTEGRITY_DOC)) else {
        violations.push(claim_integrity_violation(
            CLAIM_INTEGRITY_DOC,
            format!(
                "{CLAIM_INTEGRITY_DOC} is missing — the claim-integrity defect class is the \
                 definition the E02 sweep and promotion gate consume verbatim; without it the \
                 gate counts an inventory it cannot define (bead f85xj.2.1)"
            ),
        ));
        return violations;
    };

    for section in CLAIM_INTEGRITY_REQUIRED_SECTIONS {
        if !definition.contains(section) {
            violations.push(claim_integrity_violation(
                CLAIM_INTEGRITY_DOC,
                format!(
                    "{CLAIM_INTEGRITY_DOC} lost required section {section:?} — downstream beads \
                     (.2.2 sweep, .2.3 gate) consume these sections verbatim (bead f85xj.2.1)"
                ),
            ));
        }
    }

    for label in CLAIM_INTEGRITY_SEVERITY_LABELS {
        if !definition.contains(label) {
            violations.push(claim_integrity_violation(
                CLAIM_INTEGRITY_DOC,
                format!(
                    "{CLAIM_INTEGRITY_DOC} does not name canonical severity label {label:?} — the \
                     doc and xtask must agree on the label set the gate queries (bead f85xj.2.1)"
                ),
            ));
        }
    }

    // The inventory script is named by the definition as its enforcement arm;
    // a definition citing a script that does not exist overstates enforcement.
    if definition.contains(CLAIM_INTEGRITY_INVENTORY_SCRIPT)
        && !root.join(CLAIM_INTEGRITY_INVENTORY_SCRIPT).is_file()
    {
        violations.push(claim_integrity_violation(
            CLAIM_INTEGRITY_DOC,
            format!(
                "{CLAIM_INTEGRITY_DOC} cites {CLAIM_INTEGRITY_INVENTORY_SCRIPT} as its enforcement \
                 arm but that script does not exist — documented enforcement that cannot run is \
                 itself an overstated claim (bead f85xj.2.1)"
            ),
        ));
    }

    let Ok(conventions) = std::fs::read_to_string(root.join(CLAIM_INTEGRITY_CONVENTIONS)) else {
        violations.push(claim_integrity_violation(
            CLAIM_INTEGRITY_CONVENTIONS,
            format!("{CLAIM_INTEGRITY_CONVENTIONS} is missing (bead f85xj.2.1)"),
        ));
        return violations;
    };

    if !conventions.contains(CLAIM_INTEGRITY_CONVENTIONS_SECTION) {
        violations.push(claim_integrity_violation(
            CLAIM_INTEGRITY_CONVENTIONS,
            format!(
                "{CLAIM_INTEGRITY_CONVENTIONS} lost section \
                 {CLAIM_INTEGRITY_CONVENTIONS_SECTION:?} — the label taxonomy must be discoverable \
                 where agents read conventions, not only in the definition (bead f85xj.2.1)"
            ),
        ));
    }
    if !conventions.contains(CLAIM_INTEGRITY_LABEL) {
        violations.push(claim_integrity_violation(
            CLAIM_INTEGRITY_CONVENTIONS,
            format!(
                "{CLAIM_INTEGRITY_CONVENTIONS} does not name the {CLAIM_INTEGRITY_LABEL:?} label \
                 (bead f85xj.2.1)"
            ),
        ));
    }
    for label in CLAIM_INTEGRITY_SEVERITY_LABELS {
        if !conventions.contains(label) {
            violations.push(claim_integrity_violation(
                CLAIM_INTEGRITY_CONVENTIONS,
                format!(
                    "{CLAIM_INTEGRITY_CONVENTIONS} taxonomy omits canonical severity label \
                     {label:?} — an undocumented severity is one the sweep will not apply \
                     (bead f85xj.2.1)"
                ),
            ));
        }
    }

    violations
}

/// README claim-state lint: see module docs for the three rules.
fn check_claims_from_readme(
    root: &Path,
    readme: &str,
    mut violations: Vec<Violation>,
) -> Vec<Violation> {
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

    // Rule 5 (f85xj.2.1): the claim-integrity defect class stays defined and
    // its taxonomy stays documented where agents read conventions.
    violations.extend(check_claim_integrity_docs(root));

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

pub fn check_claim_language(root: &Path) -> Vec<Violation> {
    let Ok(readme) = std::fs::read_to_string(root.join("README.md")) else {
        return vec![Violation {
            check: "claim-state",
            crate_name: "<repo>".to_string(),
            detail: "README.md missing at workspace root".to_string(),
        }];
    };
    check_claims_from_readme(root, &readme, Vec::new())
}

pub fn check_claims(root: &Path) -> Vec<Violation> {
    let mut violations = check_docs(root).violations;
    violations.extend(check_claim_language(root));
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
    fn workspace_member_count_excludes_tools_and_nested_workspaces() {
        let manifest = r#"
[workspace]
members = [
    "crates/fs-a",
    "crates/fs-b",
    "xtask",
]
"#;
        assert_eq!(workspace_fs_member_count(manifest), Some(2));
        assert_eq!(workspace_fs_member_count("[workspace]\n"), None);
    }

    #[test]
    fn generated_doc_facts_fail_on_one_seeded_stale_count() {
        let facts = DocFacts {
            native_workspace_crates: 2,
            standalone_crates: vec!["fs-wasm".to_string()],
            crate_directories: 3,
            contracts: 3,
            tracked_integration_tests: 2,
            layers: BTreeMap::from([
                ("UTIL", 0),
                ("L0", 1),
                ("L1", 1),
                ("L2", 0),
                ("L3", 0),
                ("L4", 0),
                ("L5", 0),
                ("L6", 1),
            ]),
        };
        let generated = facts.render();
        assert!(
            check_docs_with_facts(&generated, &facts)
                .violations
                .is_empty(),
            "the exact generated block must pass"
        );

        let stale = generated.replacen(
            "| Native workspace `fs-*` crates | 2 |",
            "| Native workspace `fs-*` crates | 9 |",
            1,
        );
        let violations = check_docs_with_facts(&stale, &facts).violations;
        assert_eq!(
            violations.len(),
            1,
            "one seeded stale count: {violations:?}"
        );
        assert!(violations[0].detail.contains("block is stale"));
    }

    #[test]
    fn tracked_doc_inventory_is_portable_to_source_snapshots_without_dot_git() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("xtask must live directly below the workspace root");
        let inventory = TrackedInventory::read(root).expect("live inventory must validate");
        assert!(!inventory.manifests.is_empty());
        assert_eq!(inventory.manifests.len(), inventory.contracts.len());
        assert!(!inventory.integration_tests.is_empty());

        let (facts, git_index_verified) =
            DocFacts::derive(root).expect("documentation facts must derive from the registry");
        assert_eq!(facts.crate_directories, inventory.manifests.len());
        assert_eq!(facts.contracts, inventory.contracts.len());
        assert_eq!(
            facts.tracked_integration_tests,
            inventory.integration_tests.len()
        );
        if root.join(".git").is_dir() {
            assert!(
                git_index_verified,
                "a real worktree must verify the portable registry against Git"
            );
        }
    }

    /// Build a minimal docs pair that satisfies the claim-integrity lint, so
    /// each negative case below differs from green by exactly one mutation.
    fn claim_integrity_fixture(base: &Path) {
        let mut definition = String::new();
        for section in CLAIM_INTEGRITY_REQUIRED_SECTIONS {
            definition.push_str(section);
            definition.push_str("\n\nbody\n\n");
        }
        for label in CLAIM_INTEGRITY_SEVERITY_LABELS {
            definition.push_str(&format!("- `{label}`\n"));
        }
        let mut conventions = format!("{CLAIM_INTEGRITY_CONVENTIONS_SECTION}\n\n");
        conventions.push_str(&format!("label `{CLAIM_INTEGRITY_LABEL}`\n"));
        for label in CLAIM_INTEGRITY_SEVERITY_LABELS {
            conventions.push_str(&format!("- `{label}`\n"));
        }
        let write = |rel: &str, text: &str| {
            let path = base.join(rel);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(path, text).unwrap();
        };
        write(CLAIM_INTEGRITY_DOC, &definition);
        write(CLAIM_INTEGRITY_CONVENTIONS, &conventions);
    }

    #[test]
    fn claim_integrity_lint_accepts_a_complete_definition_and_taxonomy() {
        let base = std::env::temp_dir().join(format!("fsim-ci-ok-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        claim_integrity_fixture(&base);
        let violations = check_claim_integrity_docs(&base);
        assert!(violations.is_empty(), "expected clean: {violations:?}");
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn claim_integrity_lint_fails_closed_on_each_single_mutation() {
        let base = std::env::temp_dir().join(format!("fsim-ci-mut-{}", std::process::id()));

        // A missing definition is one violation, not a silent pass: the gate
        // must never count an inventory whose class it cannot define.
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(base.join("docs")).unwrap();
        std::fs::write(base.join(CLAIM_INTEGRITY_CONVENTIONS), "irrelevant").unwrap();
        let missing = check_claim_integrity_docs(&base);
        assert_eq!(missing.len(), 1, "{missing:?}");
        assert!(missing[0].detail.contains("is missing"));

        // Dropping any one required section is caught by name.
        for dropped in CLAIM_INTEGRITY_REQUIRED_SECTIONS {
            let _ = std::fs::remove_dir_all(&base);
            claim_integrity_fixture(&base);
            let text = std::fs::read_to_string(base.join(CLAIM_INTEGRITY_DOC)).unwrap();
            std::fs::write(
                base.join(CLAIM_INTEGRITY_DOC),
                text.replace(dropped, "## Removed"),
            )
            .unwrap();
            let violations = check_claim_integrity_docs(&base);
            assert!(
                violations.iter().any(|v| v.detail.contains(dropped)),
                "dropping {dropped:?} must be caught: {violations:?}"
            );
        }

        // Dropping any one canonical severity label is caught in both files,
        // because doc and taxonomy must agree on the set the gate queries.
        for label in CLAIM_INTEGRITY_SEVERITY_LABELS {
            let _ = std::fs::remove_dir_all(&base);
            claim_integrity_fixture(&base);
            for file in [CLAIM_INTEGRITY_DOC, CLAIM_INTEGRITY_CONVENTIONS] {
                let text = std::fs::read_to_string(base.join(file)).unwrap();
                std::fs::write(base.join(file), text.replace(label, "severity:unknown")).unwrap();
            }
            let violations = check_claim_integrity_docs(&base);
            assert!(
                violations
                    .iter()
                    .filter(|v| v.detail.contains(label))
                    .count()
                    >= 2,
                "dropping {label:?} must be caught in both files: {violations:?}"
            );
        }

        // The CONVENTIONS taxonomy section must stay discoverable.
        let _ = std::fs::remove_dir_all(&base);
        claim_integrity_fixture(&base);
        let text = std::fs::read_to_string(base.join(CLAIM_INTEGRITY_CONVENTIONS)).unwrap();
        std::fs::write(
            base.join(CLAIM_INTEGRITY_CONVENTIONS),
            text.replace(CLAIM_INTEGRITY_CONVENTIONS_SECTION, "## Something else"),
        )
        .unwrap();
        let violations = check_claim_integrity_docs(&base);
        assert!(
            violations
                .iter()
                .any(|v| v.detail.contains(CLAIM_INTEGRITY_CONVENTIONS_SECTION)),
            "{violations:?}"
        );

        // Citing an enforcement script that does not exist is itself an
        // overstated claim.
        let _ = std::fs::remove_dir_all(&base);
        claim_integrity_fixture(&base);
        let text = std::fs::read_to_string(base.join(CLAIM_INTEGRITY_DOC)).unwrap();
        std::fs::write(
            base.join(CLAIM_INTEGRITY_DOC),
            format!("{text}\nrun {CLAIM_INTEGRITY_INVENTORY_SCRIPT} for the report\n"),
        )
        .unwrap();
        let violations = check_claim_integrity_docs(&base);
        assert!(
            violations
                .iter()
                .any(|v| v.detail.contains(CLAIM_INTEGRITY_INVENTORY_SCRIPT)),
            "{violations:?}"
        );
        std::fs::create_dir_all(base.join("scripts/ci")).unwrap();
        std::fs::write(base.join(CLAIM_INTEGRITY_INVENTORY_SCRIPT), "#!/bin/sh\n").unwrap();
        assert!(
            check_claim_integrity_docs(&base).is_empty(),
            "materializing the cited script must clear the violation"
        );

        let _ = std::fs::remove_dir_all(&base);
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
            "Cargo.toml",
            "[workspace]\nmembers = [\n    \"crates/fs-real\",\n]\n",
        );
        mk(
            "crates/fs-real/src/lib.rs",
            "pub const G: u64 = 0x1111_2222_3333_4444;\n",
        );
        mk(
            "crates/fs-real/Cargo.toml",
            "[package]\nname = \"fs-real\"\n[package.metadata.frankensim]\nlayer = \"L0\"\n",
        );
        mk(
            "crates/fs-real/tests/battery.rs",
            "fn real_golden_hash() {}\n",
        );
        let facts = DocFacts {
            native_workspace_crates: 1,
            standalone_crates: Vec::new(),
            crate_directories: 1,
            contracts: 0,
            tracked_integration_tests: 1,
            layers: BTreeMap::from([
                ("UTIL", 0),
                ("L0", 1),
                ("L1", 0),
                ("L2", 0),
                ("L3", 0),
                ("L4", 0),
                ("L5", 0),
                ("L6", 0),
            ]),
        };
        let doc_facts = facts.render();
        // Seeded drift: stale hash, missing crate, missing sentinel fn.
        mk(
            "README.md",
            &format!(
                "{doc_facts}\n\nGood: `fs-real` golden `0x1111_2222_3333_4444` via `real_golden_hash`.\n\
                 Stale hash 0xaaaa_bbbb_cccc_dddd.\n\
                 Gone crate `fs-vanished`.\n\
                 Gone sentinel `ghost_golden_hash`.\n"
            ),
        );
        // Rule 5's docs are present and complete so this case still isolates
        // the three seeded README drifts; the claim-integrity lint has its own
        // mutation tests above.
        claim_integrity_fixture(&base);
        let readme = std::fs::read_to_string(base.join("README.md")).unwrap();
        let docs = check_docs_with_facts(&readme, &facts);
        let v = check_claims_from_readme(&base, &readme, docs.violations);
        assert_eq!(v.len(), 3, "exactly the three seeded drifts: {v:?}");
        assert!(v.iter().any(|x| x.detail.contains("aaaabbbbccccdddd")));
        assert!(v.iter().any(|x| x.detail.contains("fs-vanished")));
        assert!(v.iter().any(|x| x.detail.contains("ghost_golden_hash")));
        let _ = std::fs::remove_dir_all(&base);
    }
}
