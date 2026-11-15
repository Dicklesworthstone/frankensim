//! Contract cut and index-plus-parts checks (bead
//! frankensim-rc-1109-root-3gc12.2.1).
//!
//! Twenty-seven crate contracts were each cut to a fraction of their size by
//! one commit, every one at 1,051,276 bytes or more (the cause is recorded on
//! the bead). Three were restored as an index plus parts: `CONTRACT.md` keeps
//! the required headings and a table that links `contracts/part-NN.md`, and
//! the parts hold the text. This module checks for the cut itself and for
//! that layout. It sets no size limit.
//!
//! - **Parts** (`contract-parts`). Every `contracts/part-NN.md` that a crate's
//!   `CONTRACT.md` links exists, and every `part-NN.md` file in the crate's
//!   `contracts/` directory is linked from its `CONTRACT.md`.
//! - **Shrink** (`contract-shrink`). A crate's `CONTRACT.md` or
//!   `contracts/part-NN.md` that is now smaller than 80% of its size at the
//!   base revision fails, unless every line removed from it is present in
//!   another contract file of the same crate, which is what moving text into a
//!   part looks like. A file that is gone counts as 0 bytes. The present bytes
//!   are the working tree's, so an uncommitted cut fails as well as a
//!   committed one. The base revision is `HEAD~1`, the parent of the checked
//!   commit, or the revision named by `--contract-base <rev>`.
//!
//! When the base revision cannot be read (no `.git`, as in a `git archive`
//! export, or a shallow clone without the parent), the shrink check does not
//! run, and a policy note with verdict `not-run` says so. When it runs, a note
//! names the base commit and the number of files compared.
//!
//! No claim: lines are compared whole after trailing whitespace is removed,
//! blank lines are skipped, and a removed line counts as kept when it appears
//! anywhere in another contract file of the crate. A deliberate rewrite that
//! shortens a file by more than a fifth and rewords what it removed fails at
//! its commit; the violation names the file for the reviewer.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::Path;
use std::process::Command;

use crate::{Manifest, PolicyNote, Violation};

pub(crate) const PARTS_CHECK: &str = "contract-parts";
pub(crate) const SHRINK_CHECK: &str = "contract-shrink";
const DEFAULT_BASE: &str = "HEAD~1";
const BASE_FLAG: &str = "--contract-base";

/// The base revision for the shrink check: `--contract-base <rev>` or
/// `--contract-base=<rev>` on the command line, else `HEAD~1`.
pub(crate) fn base_from_args() -> String {
    let args: Vec<String> = std::env::args().collect();
    for (index, arg) in args.iter().enumerate() {
        if let Some(rest) = arg.strip_prefix("--contract-base=") {
            return rest.to_string();
        }
        if arg == BASE_FLAG
            && let Some(rev) = args.get(index + 1)
        {
            return rev.clone();
        }
    }
    DEFAULT_BASE.to_string()
}

/// Both checks over the workspace at `root`.
pub(crate) fn check(
    root: &Path,
    manifests: &[Manifest],
    base: &str,
) -> (Vec<Violation>, Vec<PolicyNote>) {
    let mut violations = parts_check(manifests);
    let mut notes = Vec::new();
    match shrink_check(root, base) {
        Ok((found, compared, commit)) => {
            violations.extend(found);
            notes.push(PolicyNote {
                check: SHRINK_CHECK,
                crate_name: "workspace".to_string(),
                verdict: "note",
                detail: format!("compared {compared} contract files with base {base} = {commit}"),
            });
        }
        Err(reason) => notes.push(PolicyNote {
            check: SHRINK_CHECK,
            crate_name: "workspace".to_string(),
            verdict: "not-run",
            detail: format!("shrink check not run against base {base}: {reason}"),
        }),
    }
    (violations, notes)
}

/// Whether `name` is `part-NN.md`, with one or more ASCII digits.
fn is_part_name(name: &str) -> bool {
    name.strip_prefix("part-")
        .and_then(|rest| rest.strip_suffix(".md"))
        .is_some_and(|digits| !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()))
}

/// Whether a workspace-relative path is a crate contract file:
/// `crates/<crate>/CONTRACT.md` or `crates/<crate>/contracts/part-NN.md`.
fn contract_crate(path: &str) -> Option<&str> {
    let rest = path.strip_prefix("crates/")?;
    let (krate, tail) = rest.split_once('/')?;
    let is_contract =
        tail == "CONTRACT.md" || tail.strip_prefix("contracts/").is_some_and(is_part_name);
    (is_contract && !krate.is_empty()).then_some(krate)
}

/// The part names (`part-NN.md`) that an index links as `contracts/part-NN.md`.
fn linked_parts(index: &str) -> BTreeSet<String> {
    const LINK: &str = "contracts/part-";
    let mut linked = BTreeSet::new();
    let mut rest = index;
    while let Some(at) = rest.find(LINK) {
        let tail = &rest[at + "contracts/".len()..];
        let end = tail
            .find(|c: char| !(c.is_ascii_alphanumeric() || c == '-' || c == '.'))
            .unwrap_or(tail.len());
        let name = tail[..end].trim_end_matches('.');
        if is_part_name(name) {
            linked.insert(name.to_string());
        }
        // `tail` starts with "part-", so `end` is past at least that.
        rest = &tail[end..];
    }
    linked
}

/// The parts check for one crate: `index` is its `CONTRACT.md` text and
/// `files` the `part-NN.md` names present in its `contracts/` directory.
fn parts_violations(crate_name: &str, index: &str, files: &BTreeSet<String>) -> Vec<Violation> {
    let linked = linked_parts(index);
    let missing = linked.difference(files).map(|part| {
        format!("{crate_name} CONTRACT.md links contracts/{part}, which does not exist")
    });
    let unlisted = files
        .difference(&linked)
        .map(|part| format!("{crate_name} contracts/{part} is not linked from CONTRACT.md"));
    missing
        .chain(unlisted)
        .map(|detail| Violation {
            check: PARTS_CHECK,
            crate_name: crate_name.to_string(),
            detail,
        })
        .collect()
}

fn parts_check(manifests: &[Manifest]) -> Vec<Violation> {
    let mut violations = Vec::new();
    for m in manifests {
        // A missing CONTRACT.md is the `contracts` check's finding.
        let Ok(index) = std::fs::read_to_string(m.dir.join("CONTRACT.md")) else {
            continue;
        };
        let mut files = BTreeSet::new();
        if let Ok(entries) = std::fs::read_dir(m.dir.join("contracts")) {
            for entry in entries.flatten() {
                let name = entry.file_name().to_string_lossy().into_owned();
                if is_part_name(&name) {
                    files.insert(name);
                }
            }
        }
        violations.extend(parts_violations(&m.name, &index, &files));
    }
    violations
}

/// Non-blank lines, trailing whitespace removed.
fn content_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .map(str::trim_end)
        .filter(|line| !line.trim_start().is_empty())
}

/// The shrink verdict for one file. `base_len` and `base_text` are the file
/// at the base revision; `now` is its present text (`None` when it is gone);
/// `others` are the present texts of the crate's other contract files.
/// Returns the violation detail, or `None` when the file passes.
fn shrink_verdict(
    path: &str,
    base_len: u64,
    base_text: &str,
    now: Option<&str>,
    others: &[&str],
    base_label: &str,
) -> Option<String> {
    let now_len = now.map_or(0, |text| text.len() as u64);
    // Under 80%: now / base < 4 / 5, in integers.
    if now_len.saturating_mul(5) >= base_len.saturating_mul(4) {
        return None;
    }
    let kept: HashSet<&str> = now
        .into_iter()
        .chain(others.iter().copied())
        .flat_map(content_lines)
        .collect();
    let mut lost_lines = 0usize;
    let mut lost_bytes = 0usize;
    let mut first = None;
    let mut seen = HashSet::new();
    for line in content_lines(base_text) {
        if kept.contains(line) || !seen.insert(line) {
            continue;
        }
        lost_lines += 1;
        lost_bytes += line.len();
        first.get_or_insert(line);
    }
    let shown: String = first?.chars().take(120).collect();
    Some(format!(
        "{path} is {now_len} bytes, under 80% of its {base_len} bytes at {base_label}; \
         {lost_lines} removed lines ({lost_bytes} bytes) are in no other contract file \
         of the crate, the first: {shown:?}"
    ))
}

fn git(root: &Path, args: &[&str]) -> Result<Vec<u8>, String> {
    let output = Command::new("git")
        .args(args)
        .current_dir(root)
        .output()
        .map_err(|error| format!("cannot run git {}: {error}", args.join(" ")))?;
    if !output.status.success() {
        return Err(format!(
            "git {} failed with status {}: {}",
            args.join(" "),
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    Ok(output.stdout)
}

/// The contract files of every crate at `commit`, with their sizes in bytes,
/// keyed by crate directory name and then by path.
fn base_contracts(
    root: &Path,
    commit: &str,
) -> Result<BTreeMap<String, BTreeMap<String, u64>>, String> {
    let listing = git(root, &["ls-tree", "-r", "-l", "-z", commit, "--", "crates"])?;
    let mut crates: BTreeMap<String, BTreeMap<String, u64>> = BTreeMap::new();
    for record in listing.split(|byte| *byte == 0) {
        if record.is_empty() {
            continue;
        }
        let record = std::str::from_utf8(record)
            .map_err(|_| "git ls-tree emitted a non-UTF-8 record".to_string())?;
        let (meta, path) = record
            .split_once('\t')
            .ok_or_else(|| format!("git ls-tree record without a tab: {record:?}"))?;
        let Some(krate) = contract_crate(path) else {
            continue;
        };
        let fields: Vec<&str> = meta.split_whitespace().collect();
        let [_mode, kind, _object, size] = fields[..] else {
            return Err(format!(
                "git ls-tree record with {} fields: {record:?}",
                fields.len()
            ));
        };
        if kind != "blob" {
            continue;
        }
        let size: u64 = size
            .parse()
            .map_err(|_| format!("git ls-tree size {size:?} for {path} is not a number"))?;
        crates
            .entry(krate.to_string())
            .or_default()
            .insert(path.to_string(), size);
    }
    Ok(crates)
}

/// The present contract files of the crate in `crates/<krate>`, by path.
fn present_contracts(root: &Path, krate: &str) -> BTreeMap<String, String> {
    let mut present = BTreeMap::new();
    let dir = root.join("crates").join(krate);
    let index = format!("crates/{krate}/CONTRACT.md");
    if let Ok(text) = std::fs::read(dir.join("CONTRACT.md")) {
        present.insert(index, String::from_utf8_lossy(&text).into_owned());
    }
    if let Ok(entries) = std::fs::read_dir(dir.join("contracts")) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            if is_part_name(&name)
                && let Ok(text) = std::fs::read(entry.path())
            {
                present.insert(
                    format!("crates/{krate}/contracts/{name}"),
                    String::from_utf8_lossy(&text).into_owned(),
                );
            }
        }
    }
    present
}

/// Shrink check against `base`: the violations, the number of base contract
/// files compared and the base commit id.
fn shrink_check(root: &Path, base: &str) -> Result<(Vec<Violation>, usize, String), String> {
    let spec = format!("{base}^{{commit}}");
    let commit = git(root, &["rev-parse", "--verify", "--quiet", &spec])?;
    let commit = String::from_utf8_lossy(&commit).trim().to_string();
    if commit.is_empty() {
        return Err(format!("{base} does not name a commit"));
    }
    let short: String = commit.chars().take(10).collect();
    let base_label = format!("{base} ({short})");
    let crates = base_contracts(root, &commit)?;
    let mut violations = Vec::new();
    let mut compared = 0usize;
    for (krate, files) in &crates {
        compared += files.len();
        let shrunk: Vec<(&String, u64)> = files
            .iter()
            .filter(|(path, base_len)| {
                let now_len = std::fs::metadata(root.join(path.as_str())).map_or(0, |m| m.len());
                now_len.saturating_mul(5) < base_len.saturating_mul(4)
            })
            .map(|(path, base_len)| (path, *base_len))
            .collect();
        if shrunk.is_empty() {
            continue;
        }
        let present = present_contracts(root, krate);
        for (path, base_len) in shrunk {
            let object = format!("{commit}:{path}");
            let base_text = git(root, &["cat-file", "blob", &object])?;
            let base_text = String::from_utf8_lossy(&base_text);
            let now = present.get(path).map(String::as_str);
            let others: Vec<&str> = present
                .iter()
                .filter(|(other, _)| *other != path)
                .map(|(_, text)| text.as_str())
                .collect();
            if let Some(detail) =
                shrink_verdict(path, base_len, &base_text, now, &others, &base_label)
            {
                violations.push(Violation {
                    check: SHRINK_CHECK,
                    crate_name: krate.clone(),
                    detail,
                });
            }
        }
    }
    Ok((violations, compared, commit))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn numbered(range: std::ops::Range<usize>) -> String {
        range.map(|i| format!("contract line {i:03}\n")).collect()
    }

    fn parts(names: &[&str]) -> BTreeSet<String> {
        names.iter().map(|name| (*name).to_string()).collect()
    }

    /// `shrink_verdict` for `crates/fs-a/CONTRACT.md` against its `base` text.
    fn verdict(base: &str, now: Option<&str>, others: &[&str]) -> Option<String> {
        let (path, size) = ("crates/fs-a/CONTRACT.md", base.len() as u64);
        shrink_verdict(path, size, base, now, others, "HEAD~1")
    }

    const INDEX: &str = "## Purpose and layer\n\n\
        | [`contracts/part-02.md`](contracts/part-02.md) | 10 | 1 | a | b |\n\
        | [`contracts/part-01.md`](contracts/part-01.md) | 10 | 1 | c | d |\n\n\
        Full text: [`contracts/part-01.md`](contracts/part-01.md), under the same heading.\n";

    #[test]
    fn part_names_and_contract_paths_are_recognised_exactly() {
        assert!(is_part_name("part-01.md"));
        assert!(is_part_name("part-123.md"));
        assert!(!is_part_name("part-.md"));
        assert!(!is_part_name("part-01.md.bak"));
        assert!(!is_part_name("part-1a.md"));
        assert_eq!(contract_crate("crates/fs-a/CONTRACT.md"), Some("fs-a"));
        assert_eq!(
            contract_crate("crates/fs-a/contracts/part-07.md"),
            Some("fs-a")
        );
        assert_eq!(contract_crate("crates/fs-a/src/CONTRACT.md"), None);
        assert_eq!(contract_crate("crates/fs-a/contracts/notes.md"), None);
        assert_eq!(contract_crate("apps/x/CONTRACT.md"), None);
        assert_eq!(
            linked_parts("see contracts/part-03.md. and contracts/part-04.md)"),
            parts(&["part-03.md", "part-04.md"])
        );
    }

    #[test]
    fn an_index_whose_parts_all_exist_and_are_all_linked_passes() {
        assert!(parts_violations("fs-a", INDEX, &parts(&["part-01.md", "part-02.md"])).is_empty());
        assert!(parts_violations("fs-a", "## Purpose and layer\n", &parts(&[])).is_empty());
    }

    #[test]
    fn an_index_that_links_a_missing_part_fails() {
        let v = parts_violations("fs-a", INDEX, &parts(&["part-01.md"]));
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].check, PARTS_CHECK);
        let detail = &v[0].detail;
        assert!(detail.contains("links contracts/part-02.md, which does not exist"));
    }

    #[test]
    fn a_part_file_the_index_does_not_link_fails() {
        let files = parts(&["part-01.md", "part-02.md", "part-03.md"]);
        let v = parts_violations("fs-a", INDEX, &files);
        assert_eq!(v.len(), 1);
        let detail = &v[0].detail;
        assert!(detail.contains("contracts/part-03.md is not linked from CONTRACT.md"));
    }

    #[test]
    fn a_contract_cut_to_half_fails() {
        let base = numbered(0..100);
        let now = numbered(0..50);
        let detail = verdict(&base, Some(&now), &[]).expect("a 50% cut must fail");
        assert!(detail.contains("50 removed lines"), "{detail}");
        assert!(detail.contains("\"contract line 050\""), "{detail}");
        // Gone outright is a cut to 0 bytes.
        assert!(verdict(&base, None, &[]).is_some());
    }

    #[test]
    fn a_shrink_whose_text_moved_to_a_part_passes() {
        let base = numbered(0..100);
        let now = format!(
            "{}Full text: [`contracts/part-01.md`](contracts/part-01.md).\n",
            numbered(0..50)
        );
        let part = format!("## Moved\n\n{}", numbered(50..100));
        assert!(verdict(&base, Some(&now), &[&part]).is_none());
        // A whole file folded into another passes the same way.
        assert!(verdict(&base, None, &[&base]).is_none());
        // One line that went nowhere fails, and is the one named.
        let lossy = format!("## Moved\n\n{}", numbered(51..100));
        let detail = verdict(&base, Some(&now), &[&lossy]).expect("a lost line must fail");
        assert!(detail.contains("1 removed lines"), "{detail}");
        assert!(detail.contains("\"contract line 050\""), "{detail}");
    }

    #[test]
    fn a_shrink_to_80_percent_or_more_is_not_judged() {
        let base = numbered(0..100);
        assert!(verdict(&base, Some(&numbered(0..80)), &[]).is_none());
        assert!(verdict(&base, Some(&numbered(0..79)), &[]).is_some());
    }

    /// The git path end to end on a scratch repository: a committed cut fails
    /// against its parent, an uncommitted cut fails against HEAD, and a base
    /// that does not exist is a `not-run` note, not a pass.
    #[test]
    fn the_shrink_check_reads_its_base_from_git() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("wall clock is after the Unix epoch")
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "xtask-contract-shrink-{}-{unique}",
            std::process::id()
        ));
        let dir = root.join("crates").join("fs-a");
        std::fs::create_dir_all(dir.join("contracts")).expect("fixture directories");
        let run = |args: &[&str]| {
            let output = Command::new("git")
                .args([
                    "-c",
                    "user.name=xtask",
                    "-c",
                    "user.email=xtask@example.invalid",
                    "-c",
                    "commit.gpgsign=false",
                ])
                .args(args)
                .current_dir(&root)
                .output()
                .expect("git runs");
            assert!(output.status.success(), "git {args:?}: {output:?}");
        };
        run(&["init", "--quiet"]);
        let contract = dir.join("CONTRACT.md");
        std::fs::write(&contract, numbered(0..100)).expect("write");
        run(&["add", "-A"]);
        run(&["commit", "--quiet", "-m", "full"]);
        // Only one commit: HEAD~1 does not exist.
        let (violations, notes) = check(&root, &[], "HEAD~1");
        assert!(violations.is_empty());
        assert_eq!(notes.len(), 1);
        assert_eq!(notes[0].verdict, "not-run", "{}", notes[0].detail);

        // An uncommitted cut fails against HEAD.
        std::fs::write(&contract, numbered(0..30)).expect("write");
        let (violations, notes) = check(&root, &[], "HEAD");
        assert_eq!(violations.len(), 1, "{violations:?}");
        assert_eq!(violations[0].check, SHRINK_CHECK);
        assert_eq!(violations[0].crate_name, "fs-a");
        assert!(violations[0].detail.contains("crates/fs-a/CONTRACT.md is "));
        assert_eq!(notes[0].verdict, "note");
        let detail = &notes[0].detail;
        assert!(detail.contains("compared 1 contract files"), "{detail}");

        // Committed, it fails against the parent, the default base.
        run(&["commit", "--quiet", "-am", "cut"]);
        let (violations, _) = check(&root, &[], DEFAULT_BASE);
        assert_eq!(violations.len(), 1, "{violations:?}");

        // Moving the cut text into a part file passes.
        let part = dir.join("contracts").join("part-01.md");
        std::fs::write(part, numbered(30..100)).expect("write");
        let (violations, _) = check(&root, &[], DEFAULT_BASE);
        assert!(violations.is_empty(), "{violations:?}");
    }
}
