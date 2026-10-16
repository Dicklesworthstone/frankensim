#![cfg_attr(windows, feature(windows_by_handle))]

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
//! - Idempotent: repeated verification runs over the same terminal source
//!   state rewrite byte-identical provenance. The initial clone receipt and
//!   its first verification receipt legitimately record different states.
//!
//! Output: one JSON line per library plus
//! `constellation-bootstrap.json` (schema
//! `frankensim-constellation-bootstrap-v2`) beside the siblings. The
//! logic mirrors `cargo run -p xtask -- bootstrap-constellation`, which
//! remains the in-workspace verifier once the workspace can build; this
//! binary is the pre-Cargo entry point for machines that cannot.

#[path = "../../../xtask/src/bootstrap_provenance.rs"]
mod bootstrap_provenance;
#[path = "../../../xtask/src/constellation_admission.rs"]
pub mod constellation_admission;
#[path = "../../../xtask/src/constellation_cleanliness.rs"]
mod constellation_cleanliness;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bootstrap_provenance::BootstrapProvenanceRow;
use bootstrap_provenance::{
    bootstrap_provenance_support_preflight, provenance_path_text, write_bootstrap_provenance,
};
use constellation_cleanliness::{
    is_redirecting_entry, pinned_repository_worktree_status, repository_worktree_status,
    sanitized_git_command, verify_two_complete_passes,
};

const BOOTSTRAP_INCOMPLETE_KEY: &str = "frankensim.bootstrapIncomplete";
const CONSTELLATION_LOCK_SCHEMA: &str = "frankensim-constellation-lock-v2";
const CONSTELLATION_LOCK_IDENTITY_VERSION: u32 = 1;
const CONSTELLATION_LOCK_IDENTITY_DOMAIN: &str = "org.frankensim.xtask.constellation-lock.v1";
const CONSTELLATION_LOCK_NOTE: &str = "lock_hash covers (lib, version, git_head) only — paths are per-machine; remote is transport for bootstrap-constellation (content identity is the git head)";
const MAX_CONSTELLATION_LOCK_BYTES: usize = 1_048_576;

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
    version: String,
    git_head: String,
    remote: String,
    path: String,
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
        let mut value = 0u16;
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

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut value = 0xcbf2_9ce4_8422_2325u64;
    for &byte in bytes {
        value ^= u64::from(byte);
        value = value.wrapping_mul(0x0000_0100_0000_01b3);
    }
    value
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

fn render_lock_rows(rows: &[LockRow], lock_hash: &str) -> String {
    let mut rendered = format!(
        "{{\n  \"schema\": \"{CONSTELLATION_LOCK_SCHEMA}\",\n  \"identity_domain\": \"{CONSTELLATION_LOCK_IDENTITY_DOMAIN}\",\n  \"identity_version\": {CONSTELLATION_LOCK_IDENTITY_VERSION},\n  \"lock_hash\": \"{lock_hash}\",\n  \"note\": \"{}\",\n  \"libraries\": [\n",
        json_escape(CONSTELLATION_LOCK_NOTE)
    );
    for (index, row) in rows.iter().enumerate() {
        let comma = if index + 1 == rows.len() { "" } else { "," };
        let _ = writeln!(
            rendered,
            "    {{\"lib\": \"{}\", \"version\": \"{}\", \"git_head\": \"{}\", \"remote\": \"{}\", \"path\": \"{}\"}}{comma}",
            json_escape(&row.lib),
            json_escape(&row.version),
            json_escape(&row.git_head),
            json_escape(&row.remote),
            json_escape(&row.path),
        );
    }
    rendered.push_str("  ]\n}\n");
    rendered
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
    let mut seen = BTreeSet::new();
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
        if !seen.insert(lib.clone()) {
            return Err(format!("duplicate constellation library {lib:?}"));
        }
        rows.push(LockRow {
            lib,
            version,
            git_head,
            remote,
            path,
        });
        if rows.len() > CONSTELLATION_REPOS.len() {
            return Err("constellation lock declares too many libraries".to_string());
        }
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
    {
        return Err("constellation lock hash is not canonical lowercase hex".to_string());
    }
    let expected: BTreeSet<_> = CONSTELLATION_REPOS.iter().map(|(lib, _)| *lib).collect();
    let declared: BTreeSet<_> = rows.iter().map(|row| row.lib.as_str()).collect();
    if declared != expected || rows.len() != expected.len() {
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
            "constellation lock hash {lock_hash} disagrees with declared rows {expected_hash}"
        ));
    }
    Ok((lock_hash, rows))
}

fn read_lock(path: &Path) -> Result<String, String> {
    let file = std::fs::File::open(path)
        .map_err(|error| format!("{} unreadable: {error}", path.display()))?;
    let limit = u64::try_from(MAX_CONSTELLATION_LOCK_BYTES + 1)
        .map_err(|_| "constellation lock read bound does not fit u64".to_string())?;
    let mut text = String::new();
    file.take(limit)
        .read_to_string(&mut text)
        .map_err(|error| format!("{} is not bounded UTF-8: {error}", path.display()))?;
    if text.len() > MAX_CONSTELLATION_LOCK_BYTES {
        return Err(format!(
            "{} exceeds the {MAX_CONSTELLATION_LOCK_BYTES}-byte parser bound",
            path.display()
        ));
    }
    Ok(text)
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

fn dirname_of(lib: &str) -> &str {
    CONSTELLATION_REPOS
        .iter()
        .find(|(l, _)| *l == lib)
        .map_or(lib, |(_, d)| d)
}

fn json_escape(value: &str) -> String {
    use std::fmt::Write as _;

    let mut out = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if control.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            other => out.push(other),
        }
    }
    out
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
    let observed = read_lock(lock_path)?;
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

fn canonical_bootstrap_root(root: &Path) -> Result<PathBuf, String> {
    let metadata = std::fs::symlink_metadata(root)
        .map_err(|error| format!("cannot inspect bootstrap root {}: {error}", root.display()))?;
    if is_redirecting_entry(&metadata) || !metadata.is_dir() {
        return Err(format!(
            "bootstrap root {} must be an ordinary directory",
            root.display()
        ));
    }
    root.canonicalize().map_err(|error| {
        format!(
            "cannot canonicalize bootstrap root {}: {error}",
            root.display()
        )
    })
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
                "{} exists but is not an ordinary directory; refusing to repurpose it",
                target.display()
            ));
        }
        if is_repository_root(target)? {
            return Ok(false);
        }
        if !directory_is_empty(target)? {
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
    } else {
        if offline {
            return Err(format!(
                "{} missing from the source cache in --offline mode",
                target.display()
            ));
        }
        std::fs::create_dir(target)
            .map_err(|error| format!("cannot create {}: {error}", target.display()))?;
    }

    git_run(target, &["init", "--quiet", "--template="])?;
    git_run(
        target,
        &["config", "--local", BOOTSTRAP_INCOMPLETE_KEY, "true"],
    )?;
    git_run(target, &["config", "--local", "core.autocrlf", "false"])?;
    Ok(true)
}

fn bootstrap_marker_present(target: &Path) -> bool {
    git_out(
        target,
        &["config", "--local", "--get", BOOTSTRAP_INCOMPLETE_KEY],
    )
    .is_ok_and(|value| value == "true")
}

fn clear_bootstrap_marker(target: &Path) -> Result<(), String> {
    git_run(
        target,
        &["config", "--local", "--unset-all", BOOTSTRAP_INCOMPLETE_KEY],
    )
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
    let status = pinned_repository_worktree_status(target, &row.git_head)?;
    if !status.is_empty() {
        return Err(format!(
            "{} is DIRTY at the locked head — a modified working tree is not the pinned \
             source (a case-folding checkout collision also surfaces here); restore or \
             replace that sibling deliberately:\n{status}",
            target.display(),
        ));
    }
    let confirmed_head = git_out(target, &["rev-parse", "HEAD"])?;
    if confirmed_head != row.git_head {
        return Err(format!(
            "{} moved while its pinned state was being verified: before={head}, after={confirmed_head}",
            target.display()
        ));
    }
    Ok(())
}

/// One library's bootstrap: verify an existing tree or initialize/fetch a
/// missing tree in place. Interrupted work remains explicitly marked and may be
/// resumed only from a clean marked repository or a clean unmarked unborn
/// repository with the exact selected origin.
struct BootstrapOutcome {
    state: &'static str,
    transport_used: bool,
}

fn bootstrap_one(
    row: &LockRow,
    dest: &Path,
    offline: bool,
    selected_transport: &str,
) -> Result<BootstrapOutcome, String> {
    let dirname = dirname_of(&row.lib);
    let target = dest.join(dirname);
    let initialized = ensure_bootstrap_repository(&target, offline)?;
    let marked = bootstrap_marker_present(&target);
    let current_head = git_out(&target, &["rev-parse", "HEAD"]).ok();

    if current_head.as_deref() == Some(row.git_head.as_str()) {
        verify_pinned_clean(row, &target)?;
        if marked {
            clear_bootstrap_marker(&target)?;
            return Ok(BootstrapOutcome {
                state: "resumed",
                transport_used: false,
            });
        }
        return Ok(BootstrapOutcome {
            state: "verified",
            transport_used: false,
        });
    }
    if !marked {
        if let Some(head) = current_head.as_deref() {
            return Err(format!(
                "{} is an ordinary existing checkout at {}, but the lock pins {}; refusing to repurpose it",
                target.display(),
                head,
                row.git_head
            ));
        }
    }

    let existing_origin = git_out(&target, &["remote", "get-url", "origin"]).ok();
    if current_head.is_none() && !marked && existing_origin.as_deref() != Some(selected_transport) {
        return Err(format!(
            "{} is an unmarked unborn checkout without the exact selected origin {selected_transport:?}; refusing to adopt it",
            target.display()
        ));
    }
    if let Some(origin) = &existing_origin {
        if origin != selected_transport {
            return Err(format!(
                "{} incomplete bootstrap origin is {origin:?}, expected {selected_transport:?}",
                target.display()
            ));
        }
    }
    let status = repository_worktree_status(&target)?;
    if !status.is_empty() {
        return Err(format!(
            "{} is an incomplete bootstrap with worktree or hidden-index changes; refusing to overwrite it:\n{status}",
            target.display(),
        ));
    }
    if !offline && selected_transport == "no-remote" {
        return Err(format!(
            "lock declares no remote for {} — re-lock on a host that has one",
            row.lib
        ));
    }

    git_run(
        &target,
        &["config", "--local", BOOTSTRAP_INCOMPLETE_KEY, "true"],
    )?;
    if existing_origin.is_none() && !offline {
        git_run(&target, &["remote", "add", "origin", selected_transport])?;
    }
    if !offline {
        git_run(
            &target,
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
        &target,
        &[
            "checkout",
            "--no-recurse-submodules",
            "--quiet",
            "--no-overwrite-ignore",
            "--detach",
            &row.git_head,
        ],
    )
    .map_err(|error| {
        format!(
            "locked revision {} unavailable from {selected_transport}: {error}",
            row.git_head
        )
    })?;
    verify_pinned_clean(row, &target)?;
    clear_bootstrap_marker(&target)?;
    Ok(BootstrapOutcome {
        state: if initialized { "cloned" } else { "resumed" },
        transport_used: !offline,
    })
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
    if let Err(error) = bootstrap_provenance_support_preflight() {
        eprintln!("frankensim-bootstrap: {error}");
        return ExitCode::FAILURE;
    }
    // Default root: the frankensim checkout this binary lives in, or cwd.
    let requested_root = root.unwrap_or_else(|| {
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
    let root = match canonical_bootstrap_root(&requested_root) {
        Ok(root) => root,
        Err(error) => {
            eprintln!("error: {error}");
            return ExitCode::FAILURE;
        }
    };
    let lock_path = root.join("constellation.lock");
    let lock_text = match read_lock(&lock_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("error: {e} — the lock IS the input (pass --root <checkout>)");
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
    let Some(dest) = root.parent().map(Path::to_path_buf) else {
        eprintln!("error: {} has no parent directory", root.display());
        return ExitCode::FAILURE;
    };
    let dest_text = match provenance_path_text(&dest) {
        Ok(dest_text) => dest_text,
        Err(error) => {
            eprintln!("frankensim-bootstrap: {error}");
            return ExitCode::FAILURE;
        }
    };
    let mut provenance = Vec::new();
    let mut failures = 0usize;
    for row in &rows {
        let dirname = dirname_of(&row.lib);
        let selected_transport = from
            .as_ref()
            .map_or_else(|| row.remote.clone(), |base| format!("{base}/{dirname}"));
        match bootstrap_one(row, &dest, offline, &selected_transport) {
            Ok(outcome) => {
                println!(
                    "{{\"check\":\"constellation-bootstrap\",\"lib\":\"{}\",\"state\":\"{}\",\"head\":\"{}\"}}",
                    json_escape(&row.lib),
                    outcome.state,
                    json_escape(&row.git_head),
                );
                provenance.push(bootstrap_provenance_row(
                    row,
                    &selected_transport,
                    outcome.transport_used,
                    outcome.state,
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
        "constellation bootstrap OK: {} libraries at their locked heads under {} (provenance: {})",
        rows.len(),
        dest.display(),
        prov_path.display()
    );
    ExitCode::SUCCESS
}
