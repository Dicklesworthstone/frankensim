#![cfg_attr(windows, feature(windows_by_handle))]

//! Offline-cache replay drills (bead 1t8i): the bootstrap binary run
//! against a SYNTHETIC constellation — a temp checkout root with its own
//! lock, and a local mirror base standing in for the network — so every
//! trust rule is exercised hermetically:
//!
//! - clean machine: missing siblings clone from the mirror, detached at
//!   the pinned head, provenance written;
//! - idempotence/replay: the first replay records verification instead of
//!   cloning, then identical verification runs write byte-identical provenance;
//! - interrupted replay: clean marked and exact-origin unborn destinations
//!   resume in place, verify the exact pin, and clear the marker;
//! - adoption refusal: ordinary non-git, dirty marked, and wrong-origin unborn
//!   destinations are never repurposed;
//! - drift refusal: a sibling at the wrong head refuses;
//! - dirty refusal: a modified tree at the right head refuses;
//! - hidden-dirt refusal: local ignore rules and index flags cannot conceal
//!   untracked or modified source from the verifier;
//! - hostile-Git refusal: inherited hooks, filters, replacement objects,
//!   grafts, redirectors, hidden index flags, Python module injection, and
//!   unlisted transport helpers cannot redefine source or execute helpers;
//! - raw-index refusal: corrupted checksums, forbidden administrative paths,
//!   FSMonitor, v4, and split-index authorities all fail closed;
//! - `--offline` refusal: a missing sibling is a structured failure and
//!   the network (here: the mirror) is never touched;
//! - CLI admission: value-taking flags refuse absent operands before any
//!   lock or sibling access.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

#[allow(dead_code)] // the replay crate includes the shared module but exercises only publication seams
#[path = "../../../xtask/src/bootstrap_provenance.rs"]
mod bootstrap_provenance;

use bootstrap_provenance::{
    BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN, BOOTSTRAP_PROVENANCE_IDENTITY_VERSION,
};

const LOCK_NOTE: &str = "lock_hash covers (lib, version, git_head) only — paths are per-machine; remote is transport for bootstrap-constellation (content identity is the git head)";
const LOCK_IDENTITY_DOMAIN: &str = "org.frankensim.xtask.constellation-lock.v1";
const LOCK_IDENTITY_VERSION: u32 = 1;
const SQLITE_SUBMODULE_PATH: &str = "legacy_sqlite_code/sqlite";
const SQLITE_LEAF_SUBMODULE_PATH: &str = "leaf_sqlite";
const TRACKED_LOCK: &str = include_str!("../../../constellation.lock");
const CHECKOUT_CONSTELLATION_SH: &str =
    include_str!("../../../scripts/ci/checkout_constellation.sh");
static NEXT_FIXTURE_SUFFIX: AtomicU64 = AtomicU64::new(0);

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut value = 0xcbf2_9ce4_8422_2325u64;
    for &byte in bytes {
        value ^= u64::from(byte);
        value = value.wrapping_mul(0x0000_0100_0000_01b3);
    }
    value
}

#[cfg(windows)]
const NULL_DEVICE: &str = "NUL";
#[cfg(not(windows))]
const NULL_DEVICE: &str = "/dev/null";

const HERMETIC_GIT_CONFIG: &[(&str, &str)] = &[
    ("commit.gpgSign", "false"),
    ("tag.gpgSign", "false"),
    ("core.hooksPath", NULL_DEVICE),
    ("init.templateDir", ""),
    ("core.autocrlf", "false"),
    ("core.fsmonitor", "false"),
    ("core.untrackedCache", "false"),
    ("core.attributesFile", NULL_DEVICE),
    ("protocol.allow", "never"),
    ("protocol.file.allow", "always"),
    ("protocol.https.allow", "always"),
    ("protocol.ssh.allow", "always"),
];

fn hermetic_git_command() -> Command {
    let mut command = Command::new("git");
    hermetic_git_environment(&mut command);
    for (key, value) in HERMETIC_GIT_CONFIG {
        command.arg("-c").arg(format!("{key}={value}"));
    }
    command
}

fn hermetic_git_environment(command: &mut Command) {
    // Fixture commands must not inherit redirectors, object overlays, config
    // injection, templates, quarantine state, or tracing/prompt helpers from
    // the host. Tests that exercise a hostile variable install it explicitly
    // *after* this baseline, and replacement-object setup therefore remains
    // replacement-aware by default.
    for (key, _) in std::env::vars_os() {
        if key.to_string_lossy().starts_with("GIT_") {
            command.env_remove(key);
        }
    }
    for key in ["SSH_ASKPASS", "GCM_INTERACTIVE"] {
        command.env_remove(key);
    }
    command
        .env("GIT_CONFIG_NOSYSTEM", "1")
        .env("GIT_ATTR_NOSYSTEM", "1")
        .env("GIT_CONFIG_GLOBAL", NULL_DEVICE)
        .env("GIT_CONFIG_SYSTEM", NULL_DEVICE)
        .env("GIT_TERMINAL_PROMPT", "0")
        .env("GIT_NO_LAZY_FETCH", "1");
}

fn git(dir: &Path, args: &[&str]) -> String {
    let out = hermetic_git_command()
        .arg("-C")
        .arg(dir)
        .args(args)
        .output()
        .expect("git spawns");
    assert!(
        out.status.success(),
        "git {args:?} in {}: {}",
        dir.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).trim().to_string()
}

fn git_succeeds(dir: &Path, args: &[&str]) -> bool {
    hermetic_git_command()
        .arg("-C")
        .arg(dir)
        .args(args)
        .output()
        .expect("git spawns")
        .status
        .success()
}

fn git_path(dir: &Path, name: &str) -> PathBuf {
    let path = PathBuf::from(git(dir, &["rev-parse", "--git-path", name]));
    if path.is_absolute() {
        path
    } else {
        dir.join(path)
    }
}

fn write_git_admin_file(dir: &Path, name: &str, contents: impl AsRef<[u8]>) {
    let path = git_path(dir, name);
    std::fs::create_dir_all(path.parent().expect("Git administrative file parent"))
        .expect("create Git administrative file parent");
    std::fs::write(path, contents).expect("write Git administrative fixture");
}

fn rewrite_primary_index_entry(
    dir: &Path,
    original: &str,
    replacement: &str,
    replacement_mode: Option<&str>,
) {
    let index_path = git_path(dir, "index");
    let object_format = git(dir, &["rev-parse", "--show-object-format=storage"]);
    let script = r#"
import hashlib
import pathlib
import sys

index_path = pathlib.Path(sys.argv[1])
original = sys.argv[2].encode()
replacement = sys.argv[3].encode()
object_format = sys.argv[4]
replacement_mode = sys.argv[5]
checksum_width = {"sha1": 20, "sha256": 32}[object_format]
checksum = {"sha1": hashlib.sha1, "sha256": hashlib.sha256}[object_format]
raw = bytearray(index_path.read_bytes())
body = raw[:-checksum_width]
if body[:4] != b"DIRC" or int.from_bytes(body[4:8], "big") not in {2, 3}:
    raise SystemExit("fixture requires a v2/v3 primary index")
entry_count = int.from_bytes(body[8:12], "big")
matches = []
cursor = 12
for entry_number in range(entry_count):
    entry_start = cursor
    flags_offset = entry_start + 40 + checksum_width
    flags = int.from_bytes(body[flags_offset:flags_offset + 2], "big")
    path_start = flags_offset + 2
    if flags & 0x4000:
        path_start += 2
    encoded_length = flags & 0x0FFF
    if encoded_length < 0x0FFF:
        path_end = path_start + encoded_length
    else:
        path_end = body.find(b"\0", path_start)
        if path_end < 0:
            raise SystemExit(f"entry {entry_number} has no path terminator")
    if body[path_end:path_end + 1] != b"\0":
        raise SystemExit(f"entry {entry_number} has no bounded path terminator")
    padding_size = 8 - ((path_end - entry_start) % 8)
    entry_end = path_end + padding_size
    if body[path_start:path_end] == original:
        matches.append((entry_start, flags_offset, flags, path_start, entry_end))
    cursor = entry_end
if len(matches) != 1:
    raise SystemExit(f"expected exactly one index path {original!r}, found {matches!r}")
entry_start, flags_offset, flags, path_start, entry_end = matches[0]
if len(replacement) >= 0x0FFF:
    encoded_replacement_length = 0x0FFF
else:
    encoded_replacement_length = len(replacement)
body[flags_offset:flags_offset + 2] = (
    (flags & ~0x0FFF) | encoded_replacement_length
).to_bytes(2, "big")
if replacement_mode:
    body[entry_start + 24:entry_start + 28] = int(replacement_mode, 8).to_bytes(4, "big")
rebuilt_entry = bytearray(replacement)
rebuilt_entry.append(0)
while (path_start + len(rebuilt_entry) - entry_start) % 8:
    rebuilt_entry.append(0)
body = body[:path_start] + rebuilt_entry + body[entry_end:]
index_path.write_bytes(body + checksum(body).digest())
"#;
    let output = Command::new("python3")
        .arg("-I")
        .arg("-c")
        .arg(script)
        .arg(&index_path)
        .arg(original)
        .arg(replacement)
        .arg(&object_format)
        .arg(replacement_mode.unwrap_or(""))
        .output()
        .expect("Python index fixture editor spawns");
    assert!(
        output.status.success(),
        "rewrite primary-index path: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// A synthetic upstream repo with one committed file; returns its head.
fn make_upstream(base: &Path, name: &str, content: &str) -> String {
    let dir = base.join(name);
    std::fs::create_dir_all(&dir).expect("mkdir");
    git(&dir, &["init", "-q", "-b", "main"]);
    git(&dir, &["config", "user.email", "drill@frankensim.test"]);
    git(&dir, &["config", "user.name", "bootstrap drill"]);
    std::fs::write(dir.join("lib.rs"), content).expect("write");
    git(&dir, &["add", "lib.rs"]);
    git(&dir, &["commit", "-qm", "pinned"]);
    git(&dir, &["rev-parse", "HEAD"])
}

struct Constellation {
    /// Unique retained fixture root for this test invocation.
    base: PathBuf,
    root: PathBuf,
    mirror: PathBuf,
    heads: Vec<(String, String)>,
    lock_hash: String,
    sqlite_submodule_upstream: PathBuf,
    sqlite_submodule_drift_head: String,
}

/// Build a synthetic world: upstreams, a bare-clone mirror base, and a
/// fake frankensim checkout (parent = the bootstrap destination) whose
/// lock pins the upstream heads. Library names deliberately NOT in the
/// tool's rename map, so dirname == lib.
fn make_constellation(tag: &str) -> Constellation {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("fixture clock follows the Unix epoch")
        .as_nanos();
    let suffix = NEXT_FIXTURE_SUFFIX.fetch_add(1, Ordering::Relaxed);
    let base = std::env::temp_dir().join(format!(
        "fs-bootstrap-drill-{tag}-{}-{unique}-{suffix}",
        std::process::id()
    ));
    let upstreams = base.join("upstreams");
    let mirror = base.join("mirror");
    let parent = base.join("dest");
    let root = parent.join("frankensim");
    std::fs::create_dir_all(&root).expect("mkdir root");
    std::fs::create_dir_all(&mirror).expect("mkdir mirror");
    let sqlite_leaf_upstream = upstreams.join("sqlite-leaf-upstream");
    std::fs::create_dir_all(&sqlite_leaf_upstream).expect("mkdir sqlite leaf upstream");
    git(&sqlite_leaf_upstream, &["init", "-q", "-b", "main"]);
    git(
        &sqlite_leaf_upstream,
        &["config", "user.email", "drill@frankensim.test"],
    );
    git(
        &sqlite_leaf_upstream,
        &["config", "user.name", "bootstrap drill"],
    );
    std::fs::write(
        sqlite_leaf_upstream.join("sqlite_leaf.c"),
        "/* pinned sqlite leaf fixture */\n",
    )
    .expect("write sqlite leaf fixture");
    git(&sqlite_leaf_upstream, &["add", "sqlite_leaf.c"]);
    git(
        &sqlite_leaf_upstream,
        &["commit", "-qm", "pinned sqlite leaf"],
    );
    let sqlite_submodule_upstream = upstreams.join("sqlite-submodule-upstream");
    std::fs::create_dir_all(&sqlite_submodule_upstream).expect("mkdir sqlite submodule upstream");
    git(&sqlite_submodule_upstream, &["init", "-q", "-b", "main"]);
    git(
        &sqlite_submodule_upstream,
        &["config", "user.email", "drill@frankensim.test"],
    );
    git(
        &sqlite_submodule_upstream,
        &["config", "user.name", "bootstrap drill"],
    );
    std::fs::write(
        sqlite_submodule_upstream.join("sqlite3.c"),
        "/* pinned sqlite fixture */\n",
    )
    .expect("write pinned sqlite fixture");
    git(&sqlite_submodule_upstream, &["add", "sqlite3.c"]);
    git(
        &sqlite_submodule_upstream,
        &["commit", "-qm", "pinned sqlite gitlink"],
    );
    let out = hermetic_git_command()
        .arg("-C")
        .arg(&sqlite_submodule_upstream)
        .args([
            "-c",
            "protocol.file.allow=always",
            "submodule",
            "add",
            "--quiet",
        ])
        .arg(&sqlite_leaf_upstream)
        .arg(SQLITE_LEAF_SUBMODULE_PATH)
        .output()
        .expect("git sqlite leaf submodule add spawns");
    assert!(
        out.status.success(),
        "add sqlite leaf submodule: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    git(
        &sqlite_submodule_upstream,
        &[
            "config",
            "-f",
            ".gitmodules",
            "submodule.leaf_sqlite.ignore",
            "all",
        ],
    );
    git(
        &sqlite_submodule_upstream,
        &["add", ".gitmodules", SQLITE_LEAF_SUBMODULE_PATH],
    );
    git(
        &sqlite_submodule_upstream,
        &["commit", "-qm", "pin sqlite leaf submodule"],
    );
    let mut heads = Vec::new();
    for name in [
        "asupersync",
        "franken_networkx",
        "franken_numpy",
        "frankenpandas",
        "frankenscipy",
        "frankensqlite",
        "frankentorch",
    ] {
        let mut head = make_upstream(
            &upstreams,
            name,
            &format!("pub fn {}_fixture() {{}}\n", name.replace('-', "_")),
        );
        if name == "frankensqlite" {
            let outer = upstreams.join(name);
            let out = hermetic_git_command()
                .arg("-C")
                .arg(&outer)
                .args([
                    "-c",
                    "protocol.file.allow=always",
                    "submodule",
                    "add",
                    "--quiet",
                ])
                .arg(&sqlite_submodule_upstream)
                .arg(SQLITE_SUBMODULE_PATH)
                .output()
                .expect("git submodule add spawns");
            assert!(
                out.status.success(),
                "add sqlite-style submodule: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            git(
                &outer,
                &[
                    "config",
                    "-f",
                    ".gitmodules",
                    "submodule.legacy_sqlite_code/sqlite.ignore",
                    "dirty",
                ],
            );
            git(&outer, &["add", ".gitmodules", SQLITE_SUBMODULE_PATH]);
            git(&outer, &["commit", "-qm", "pin sqlite submodule"]);
            head = git(&outer, &["rev-parse", "HEAD"]);
        }
        // Bare mirror clone = the offline cache / air-gapped transport.
        let out = hermetic_git_command()
            .args(["clone", "-q", "--bare"])
            .arg(upstreams.join(name))
            .arg(mirror.join(name))
            .output()
            .expect("git spawns");
        assert!(out.status.success(), "bare mirror clone");
        heads.push((name.to_string(), head));
    }
    std::fs::write(
        sqlite_submodule_upstream.join("sqlite3.c"),
        "/* drifted sqlite fixture */\n",
    )
    .expect("write drifted sqlite fixture");
    git(&sqlite_submodule_upstream, &["add", "sqlite3.c"]);
    git(
        &sqlite_submodule_upstream,
        &["commit", "-qm", "drift sqlite submodule"],
    );
    let sqlite_submodule_drift_head = git(&sqlite_submodule_upstream, &["rev-parse", "HEAD"]);
    heads.sort_by(|left, right| left.0.cmp(&right.0));
    let identity = heads
        .iter()
        .map(|(name, head)| format!("{name}=0.0.0@{head}\n"))
        .collect::<String>();
    let lock_hash = format!("{:016x}", fnv1a64(identity.as_bytes()));
    let rows: Vec<String> = heads
        .iter()
        .map(|(name, head)| {
            format!(
                "    {{\"lib\": \"{name}\", \"version\": \"0.0.0\", \"git_head\": \"{head}\", \"remote\": \"no-remote\", \"path\": \"unused\"}}"
            )
        })
        .collect();
    let lock = format!(
        "{{\n  \"schema\": \"frankensim-constellation-lock-v2\",\n  \"identity_domain\": \"{LOCK_IDENTITY_DOMAIN}\",\n  \"identity_version\": {LOCK_IDENTITY_VERSION},\n  \"lock_hash\": \"{lock_hash}\",\n  \"note\": \"{LOCK_NOTE}\",\n  \"libraries\": [\n{}\n  ]\n}}\n",
        rows.join(",\n")
    );
    std::fs::write(root.join("constellation.lock"), lock).expect("write lock");
    Constellation {
        base,
        root,
        mirror,
        heads,
        lock_hash,
        sqlite_submodule_upstream,
        sqlite_submodule_drift_head,
    }
}

fn bootstrap_command(c: &Constellation, extra: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim-bootstrap"));
    command.arg("--root").arg(&c.root).args(extra);
    hermetic_git_environment(&mut command);
    command
}

fn run_bootstrap(c: &Constellation, extra: &[&str]) -> (bool, String) {
    let mut command = bootstrap_command(c, extra);
    run_command(&mut command)
}

fn run_command(command: &mut Command) -> (bool, String) {
    let out = command.output().expect("bootstrap binary spawns");
    let text = format!(
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    (out.status.success(), text)
}

fn identity_tamper_cases(canonical: &str) -> Vec<(&'static str, String)> {
    let domain_line = format!("  \"identity_domain\": \"{LOCK_IDENTITY_DOMAIN}\",\n");
    let version_line = format!("  \"identity_version\": {LOCK_IDENTITY_VERSION},\n");
    vec![
        (
            "missing identity domain",
            canonical.replacen(&domain_line, "", 1),
        ),
        (
            "wrong identity domain",
            canonical.replacen(
                LOCK_IDENTITY_DOMAIN,
                "org.frankensim.xtask.constellation-lock.v0",
                1,
            ),
        ),
        (
            "duplicate identity domain",
            canonical.replacen(&domain_line, &format!("{domain_line}{domain_line}"), 1),
        ),
        (
            "type-invalid identity domain",
            canonical.replacen(
                &format!("\"identity_domain\": \"{LOCK_IDENTITY_DOMAIN}\""),
                "\"identity_domain\": 1",
                1,
            ),
        ),
        (
            "missing identity version",
            canonical.replacen(&version_line, "", 1),
        ),
        (
            "wrong identity version",
            canonical.replacen(
                &format!("\"identity_version\": {LOCK_IDENTITY_VERSION}"),
                "\"identity_version\": 0",
                1,
            ),
        ),
        (
            "duplicate identity version",
            canonical.replacen(&version_line, &format!("{version_line}{version_line}"), 1),
        ),
        (
            "type-invalid identity version",
            canonical.replacen(
                &format!("\"identity_version\": {LOCK_IDENTITY_VERSION}"),
                "\"identity_version\": true",
                1,
            ),
        ),
        (
            "unknown identity field",
            canonical.replacen(
                &version_line,
                &format!("{version_line}  \"identity_epoch\": 1,\n"),
                1,
            ),
        ),
    ]
}

fn noncanonical_encoding_cases(canonical: &str) -> Vec<(&'static str, String)> {
    let schema_line = "  \"schema\": \"frankensim-constellation-lock-v2\",\n";
    let domain_line = format!("  \"identity_domain\": \"{LOCK_IDENTITY_DOMAIN}\",\n");
    vec![
        (
            "noncanonical whitespace",
            canonical.replacen("  \"schema\": ", "  \"schema\":", 1),
        ),
        (
            "reordered top-level keys",
            canonical.replacen(
                &format!("{schema_line}{domain_line}"),
                &format!("{domain_line}{schema_line}"),
                1,
            ),
        ),
        (
            "reordered row keys",
            canonical.replacen(
                "    {\"lib\": \"asupersync\", \"version\": \"0.0.0\"",
                "    {\"version\": \"0.0.0\", \"lib\": \"asupersync\"",
                1,
            ),
        ),
        ("noncanonical CRLF", canonical.replace('\n', "\r\n")),
    ]
}

fn install_shell_checkout(c: &Constellation) {
    let script = c.root.join("scripts/ci/checkout_constellation.sh");
    std::fs::create_dir_all(script.parent().expect("script parent")).expect("mkdir scripts/ci");
    std::fs::write(&script, CHECKOUT_CONSTELLATION_SH).expect("write checkout script");
}

fn shell_checkout_command(c: &Constellation, mode: Option<&str>) -> Command {
    let mut command = Command::new("bash");
    command.arg(c.root.join("scripts/ci/checkout_constellation.sh"));
    if let Some(mode) = mode {
        command.arg(mode);
    }
    command.arg(c.root.parent().expect("constellation destination"));
    hermetic_git_environment(&mut command);
    command
}

fn run_shell_checkout(c: &Constellation, mode: Option<&str>) -> (bool, String) {
    let mut command = shell_checkout_command(c, mode);
    run_command(&mut command)
}

fn use_local_mirror_as_locked_transport(c: &Constellation) {
    let lock_path = c.root.join("constellation.lock");
    let mut lock = std::fs::read_to_string(&lock_path).expect("fixture lock");
    for (name, _) in &c.heads {
        let remote = c.mirror.join(name).display().to_string();
        lock = lock.replacen(
            "\"remote\": \"no-remote\"",
            &format!("\"remote\": \"{remote}\""),
            1,
        );
    }
    std::fs::write(lock_path, lock).expect("write local transports");
}

fn commit_synthetic_root(c: &Constellation) {
    git(&c.root, &["init", "-q", "-b", "main"]);
    git(&c.root, &["config", "user.email", "drill@frankensim.test"]);
    git(&c.root, &["config", "user.name", "bootstrap drill"]);
    git(
        &c.root,
        &[
            "add",
            "constellation.lock",
            "scripts/ci/checkout_constellation.sh",
        ],
    );
    git(&c.root, &["commit", "-qm", "synthetic root"]);
}

fn forced_submodule_status(repository: &Path) -> String {
    git(
        repository,
        &[
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--ignore-submodules=none",
        ],
    )
}

fn initialize_sqlite_submodule(c: &Constellation) -> PathBuf {
    let outer = c
        .root
        .parent()
        .expect("constellation destination")
        .join("frankensqlite");
    let before = git(&outer, &["submodule", "status", SQLITE_SUBMODULE_PATH]);
    assert!(
        before.starts_with('-'),
        "production bootstrap must not initialize absent submodules: {before}"
    );
    git(
        &outer,
        &[
            "-c",
            "protocol.file.allow=always",
            "submodule",
            "update",
            "--init",
            "--quiet",
            "--",
            SQLITE_SUBMODULE_PATH,
        ],
    );
    let nested = outer.join(SQLITE_SUBMODULE_PATH);
    assert_eq!(forced_submodule_status(&outer), "");
    nested
}

fn initialize_sqlite_leaf(sqlite_checkout: &Path) -> PathBuf {
    let before = git(
        sqlite_checkout,
        &["submodule", "status", SQLITE_LEAF_SUBMODULE_PATH],
    );
    assert!(
        before.starts_with('-'),
        "clean-tree verification must not initialize a nested leaf: {before}"
    );
    git(
        sqlite_checkout,
        &[
            "-c",
            "protocol.file.allow=always",
            "submodule",
            "update",
            "--init",
            "--quiet",
            "--",
            SQLITE_LEAF_SUBMODULE_PATH,
        ],
    );
    sqlite_checkout.join(SQLITE_LEAF_SUBMODULE_PATH)
}

fn materialize_locked_siblings_without_provenance(c: &Constellation) {
    let dest = c.root.parent().expect("constellation destination");
    for (name, head) in &c.heads {
        let out = hermetic_git_command()
            .args(["clone", "--quiet"])
            .arg(c.mirror.join(name))
            .arg(dest.join(name))
            .output()
            .expect("git clone fixture spawns");
        assert!(
            out.status.success(),
            "clone {name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let checkout = dest.join(name);
        git(&checkout, &["checkout", "--quiet", "--detach", head]);
    }
    assert!(
        !dest.join("constellation-bootstrap.json").exists(),
        "fixture materialization must not synthesize verified provenance"
    );
}

fn install_marked_wrong_head_ignored_collision(c: &Constellation) -> (PathBuf, Vec<u8>) {
    let target = c
        .root
        .parent()
        .expect("constellation destination")
        .join("asupersync");
    let pinned = c
        .heads
        .iter()
        .find(|(name, _)| name == "asupersync")
        .map(|(_, head)| head.as_str())
        .expect("asupersync pin");
    git(&target, &["config", "user.email", "drill@frankensim.test"]);
    git(&target, &["config", "user.name", "bootstrap drill"]);

    // Construct a wrong child commit without deleting lib.rs from the retained
    // fixture: the wrong tree tracks an ignore rule and omits lib.rs, while the
    // worktree keeps lib.rs as an ignored, untracked owner-controlled file.
    std::fs::write(target.join(".gitignore"), "lib.rs\n").expect("write tracked ignore rule");
    git(&target, &["add", "--", ".gitignore"]);
    git(&target, &["update-index", "--force-remove", "--", "lib.rs"]);
    git(
        &target,
        &["commit", "-qm", "wrong head ignores pinned source"],
    );
    let wrong_head = git(&target, &["rev-parse", "HEAD"]);
    assert_ne!(wrong_head, pinned, "fixture must be at the wrong head");

    let sentinel = b"owner-controlled ignored collision bytes\n".to_vec();
    std::fs::write(target.join("lib.rs"), &sentinel).expect("write ignored collision sentinel");
    git(
        &target,
        &[
            "config",
            "--local",
            "frankensim.bootstrapIncomplete",
            "true",
        ],
    );
    let ignored = git(
        &target,
        &["status", "--ignored", "--porcelain", "--", "lib.rs"],
    );
    assert_eq!(
        ignored, "!! lib.rs",
        "fixture must expose a real ignored, untracked checkout collision"
    );
    (target, sentinel)
}

#[test]
fn value_taking_flags_refuse_missing_operands() {
    for args in [
        &["--root"][..],
        &["--from"][..],
        &["--root", "--offline"][..],
        &["--from", "--offline"][..],
    ] {
        let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim-bootstrap"));
        command.args(args);
        let (ok, text) = run_command(&mut command);
        assert!(!ok, "missing operand must refuse for {args:?}:\n{text}");
        assert!(
            text.contains(&format!("{} requires a non-empty value", args[0])),
            "{args:?} reports its admission defect:\n{text}"
        );
    }
}

#[test]
fn shell_status_probes_force_initialized_submodules_visible() {
    assert_eq!(
        CHECKOUT_CONSTELLATION_SH
            .matches("\"--ignore-submodules=all\"")
            .count(),
        0,
        "no shell status probe may suppress initialized submodule state"
    );
    assert_eq!(
        CHECKOUT_CONSTELLATION_SH
            .matches("\"--ignore-submodules=none\"")
            .count(),
        4,
        "the two cached diffs and two post-admission status probes must force visibility"
    );
}

#[test]
fn canonical_lock_tamper_and_traversal_rows_refuse_before_destination_access() {
    let c = make_constellation("lock-tamper");
    let lock_path = c.root.join("constellation.lock");
    let canonical = std::fs::read_to_string(&lock_path).expect("canonical fixture lock");
    let mut cases = identity_tamper_cases(&canonical);
    cases.extend(noncanonical_encoding_cases(&canonical));
    cases.extend([
        (
            "wrong schema",
            canonical.replacen(
                "frankensim-constellation-lock-v2",
                "frankensim-constellation-lock-v3",
                1,
            ),
        ),
        (
            "wrong row hash",
            canonical.replacen(&c.lock_hash, "0000000000000000", 1),
        ),
        (
            "duplicate library",
            canonical.replacen(
                "\"lib\": \"franken_networkx\"",
                "\"lib\": \"asupersync\"",
                1,
            ),
        ),
        (
            "path-traversing library",
            canonical.replacen("\"lib\": \"asupersync\"", "\"lib\": \"../../escaped\"", 1),
        ),
        ("trailing data", format!("{canonical}trailing-data\n")),
        ("oversized lock", "x".repeat(1_048_577)),
    ]);
    for (case, tampered) in cases {
        std::fs::write(&lock_path, tampered).expect("write tampered lock");
        let (ok, text) = run_bootstrap(&c, &["--offline"]);
        assert!(!ok, "tampered lock case {case} must refuse:\n{text}");
        assert!(
            !c.root.parent().unwrap().join("asupersync").exists(),
            "case {case} reached sibling materialization"
        );
        assert!(
            !c.base.join("escaped").exists(),
            "case {case} escaped the constellation destination"
        );
    }
}

#[test]
fn shell_lock_identity_tamper_refuses_before_repository_or_destination_access() {
    let c = make_constellation("shell-lock-identity-tamper");
    install_shell_checkout(&c);
    let lock_path = c.root.join("constellation.lock");
    let canonical = std::fs::read_to_string(&lock_path).expect("canonical fixture lock");

    let mut cases = identity_tamper_cases(&canonical);
    cases.extend(noncanonical_encoding_cases(&canonical));
    cases.push(("oversized lock", "x".repeat(1_048_577)));
    for (case, tampered) in cases {
        std::fs::write(&lock_path, tampered).expect("write tampered lock");
        let (ok, text) = run_shell_checkout(&c, None);
        assert!(!ok, "shell accepted {case}:\n{text}");
        assert!(
            text.contains("could not parse"),
            "shell must classify {case} as a lock refusal before pin checks:\n{text}"
        );
        assert!(
            !c.root.parent().unwrap().join("asupersync").exists(),
            "shell case {case} reached repository materialization"
        );
        assert!(
            !c.root
                .parent()
                .unwrap()
                .join("constellation-bootstrap.json")
                .exists(),
            "shell case {case} published destination provenance"
        );
    }
}

#[test]
fn tracked_xtask_lock_reaches_every_consumer_before_pin_checks() {
    let c = make_constellation("tracked-lock-consumers");
    std::fs::write(c.root.join("constellation.lock"), TRACKED_LOCK)
        .expect("install tracked xtask lock");
    install_shell_checkout(&c);

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(!standalone_ok, "synthetic cache is intentionally empty");
    assert!(
        standalone_text.contains("missing from the source cache in --offline mode"),
        "standalone must parse the tracked producer bytes before reporting missing pins:\n{standalone_text}"
    );
    assert!(
        !standalone_text.contains("expected canonical token")
            && !standalone_text.contains("unsupported constellation lock")
            && !standalone_text.contains("not canonical"),
        "tracked xtask lock was rejected as grammar drift:\n{standalone_text}"
    );

    let (shell_ok, shell_text) = run_shell_checkout(&c, Some("--verify-only"));
    assert!(!shell_ok, "synthetic shell cache is intentionally empty");
    assert!(
        shell_text.contains("required constellation sibling") && shell_text.contains("is missing"),
        "shell must parse the tracked producer bytes before reporting missing pins:\n{shell_text}"
    );
    assert!(
        !shell_text.contains("could not parse"),
        "tracked xtask lock was rejected by the shell grammar:\n{shell_text}"
    );
}

#[test]
fn shell_checkout_verify_and_snapshot_share_the_canonical_lock_grammar() {
    let c = make_constellation("shell-modes");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);

    let lock = std::fs::read_to_string(c.root.join("constellation.lock")).expect("fixture lock");
    assert!(
        lock.contains(&format!("\"lock_hash\": \"{}\"", c.lock_hash)),
        "transport-only fixture rewrites must preserve row-only lock_hash semantics"
    );

    let (checkout_ok, checkout_text) = run_shell_checkout(&c, None);
    assert!(
        checkout_ok,
        "synthetic shell checkout failed:\n{checkout_text}"
    );
    for (name, head) in &c.heads {
        assert!(
            checkout_text.contains(&format!(
                "\"constellation\":\"{name}\",\"verdict\":\"cloned\""
            )),
            "{name} was not cloned through the canonical shell parser:\n{checkout_text}"
        );
        let checkout = c.root.parent().unwrap().join(name);
        assert_eq!(&git(&checkout, &["rev-parse", "HEAD"]), head);
        assert_eq!(forced_submodule_status(&checkout), "");
    }

    let nested = initialize_sqlite_submodule(&c);
    assert!(nested.join("sqlite3.c").is_file());
    let (checkout_replay_ok, checkout_replay_text) = run_shell_checkout(&c, None);
    assert!(
        checkout_replay_ok,
        "clean initialized submodule must pass shell checkout replay:\n{checkout_replay_text}"
    );

    let (verify_ok, verify_text) = run_shell_checkout(&c, Some("--verify-only"));
    assert!(verify_ok, "synthetic shell verify failed:\n{verify_text}");
    for (name, _) in &c.heads {
        assert!(
            verify_text.contains(&format!(
                "\"constellation\":\"{name}\",\"verdict\":\"verified\""
            )),
            "{name} was not verified through the canonical shell parser:\n{verify_text}"
        );
    }

    let (snapshot_ok, snapshot_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        snapshot_ok,
        "synthetic shell snapshot failed:\n{snapshot_text}"
    );
    let snapshot = snapshot_text.trim();
    assert_eq!(
        snapshot.len(),
        64,
        "snapshot must be one SHA-256 hex identity"
    );
    assert!(
        snapshot
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()),
        "snapshot must be canonical lowercase hex: {snapshot}"
    );

    let (replay_ok, replay_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(replay_ok, "second shell snapshot failed:\n{replay_text}");
    assert_eq!(
        snapshot,
        replay_text.trim(),
        "unchanged shell snapshot must replay bitwise"
    );

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        standalone_ok,
        "clean initialized submodule must pass standalone replay:\n{standalone_text}"
    );
    let provenance = std::fs::read_to_string(
        c.root
            .parent()
            .unwrap()
            .join("constellation-bootstrap.json"),
    )
    .expect("standalone replay provenance");
    assert!(provenance.contains(BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN));
    assert!(provenance.contains(&format!(
        "\"identity_version\": {BOOTSTRAP_PROVENANCE_IDENTITY_VERSION}"
    )));
}

#[test]
fn shell_python_entrypoints_ignore_inherited_module_injection() {
    let c = make_constellation("pythonpath-injection");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    let hostile_modules = c.base.join("hostile-python-modules");
    std::fs::create_dir_all(&hostile_modules).expect("mkdir hostile Python module path");
    let marker = c.base.join("hostile-json-module-executed");
    std::fs::write(
        hostile_modules.join("json.py"),
        format!(
            "open({:?}, 'w').write('executed\\n')\nraise RuntimeError('hostile json module')\n",
            marker.display().to_string()
        ),
    )
    .expect("write hostile json module");

    let mut command = shell_checkout_command(&c, None);
    command.env("PYTHONPATH", &hostile_modules);
    let (ok, text) = run_command(&mut command);

    assert!(ok, "isolated shell Python bootstrap failed:\n{text}");
    assert!(
        !marker.exists(),
        "embedded Python imported a hostile PYTHONPATH module before lock admission"
    );
}

#[test]
fn fresh_init_ignores_inherited_object_ref_format_and_backend_defaults() {
    for consumer in ["standalone", "shell"] {
        let c = make_constellation(&format!("hostile-init-defaults-{consumer}"));
        let mut command = if consumer == "shell" {
            use_local_mirror_as_locked_transport(&c);
            install_shell_checkout(&c);
            shell_checkout_command(&c, None)
        } else {
            let mirror = c.mirror.to_str().expect("UTF-8 mirror");
            bootstrap_command(&c, &["--from", mirror])
        };
        let trace_marker = c.base.join("inherited-git-trace");
        let trace2_marker = c.base.join("inherited-git-trace2-event");
        command
            .env("GIT_DEFAULT_HASH", "sha256")
            .env("GIT_DEFAULT_REF_FORMAT", "reftable")
            .env("GIT_REFERENCE_BACKEND", "invalid://must-not-be-used")
            .env("GIT_TRACE", &trace_marker)
            .env("GIT_TRACE2_EVENT", &trace2_marker);
        let (ok, text) = run_command(&mut command);
        assert!(
            ok,
            "{consumer} inherited incompatible fresh-init defaults:\n{text}"
        );
        assert!(
            !trace_marker.exists() && !trace2_marker.exists(),
            "{consumer} honored inherited Git trace destinations"
        );
        for (name, head) in &c.heads {
            let checkout = c.root.parent().unwrap().join(name);
            assert_eq!(
                git(&checkout, &["rev-parse", "HEAD"]),
                *head,
                "{consumer} must initialize the lock-compatible SHA-1 repository"
            );
        }
    }
}

#[test]
fn clean_machine_clone_then_idempotent_replay_from_offline_mirror() {
    let c = make_constellation("replay");
    let mirror = c.mirror.to_str().expect("utf8").to_string();
    // CLEAN MACHINE: both siblings missing → cloned from the mirror
    // (a local path — no network exists in this test).
    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(ok, "clean-machine bootstrap succeeds:\n{text}");
    for (name, head) in &c.heads {
        assert!(
            text.contains(&format!(
                "\"lib\":\"{name}\",\"state\":\"cloned\",\"head\":\"{head}\""
            )),
            "{name} cloned at pin:\n{text}"
        );
        // Detached at the pin and clean. (The clone naturally carries
        // the upstream's default-branch ref; the contract is that HEAD
        // is DETACHED at the pin, not that refs are absent.)
        let dir = c.root.parent().unwrap().join(name);
        assert_eq!(&git(&dir, &["rev-parse", "HEAD"]), head);
        assert_eq!(forced_submodule_status(&dir), "");
        let detached = hermetic_git_command()
            .args(["-C"])
            .arg(&dir)
            .args(["symbolic-ref", "-q", "HEAD"])
            .output()
            .expect("git spawns");
        assert!(!detached.status.success(), "{name}: HEAD must be detached");
    }
    let prov_path = c
        .root
        .parent()
        .unwrap()
        .join("constellation-bootstrap.json");
    let prov1 = std::fs::read_to_string(&prov_path).expect("provenance written");
    let provenance_prefix = format!(
        "{{\n\"schema\": \"frankensim-constellation-bootstrap-v2\",\n\
         \"identity_domain\": \"{BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN}\",\n\
         \"identity_version\": {BOOTSTRAP_PROVENANCE_IDENTITY_VERSION},\n"
    );
    assert!(
        prov1.starts_with(&provenance_prefix),
        "standalone and xtask provenance headers must match exactly: {prov1}"
    );
    assert!(
        prov1.contains(&c.lock_hash),
        "canonical lock hash bound: {prov1}"
    );
    assert!(prov1.contains("\"remote\": \"no-remote\""));
    assert!(prov1.contains("\"transport_used\": true"));
    for (name, _) in &c.heads {
        assert!(
            prov1.contains(&format!("\"selected_transport\": \"{mirror}/{name}\"")),
            "mirror transport is retained for {name}: {prov1}"
        );
    }

    // REPLAY: the second run is pure verification, so transport and source-state
    // fields legitimately differ from the initial clone receipt.
    let (ok2, text2) = run_bootstrap(&c, &["--offline"]);
    assert!(
        ok2,
        "offline replay over a populated cache passes:\n{text2}"
    );
    for (name, _) in &c.heads {
        assert!(
            text2.contains(&format!("\"lib\":\"{name}\",\"state\":\"verified\"")),
            "{name} verified on replay:\n{text2}"
        );
    }
    let prov2 = std::fs::read_to_string(&prov_path).expect("provenance rewritten");
    assert!(prov2.contains("frankensim-constellation-bootstrap-v2"));
    assert!(prov2.contains("\"selected_transport\": \"no-remote\""));
    assert!(prov2.contains("\"transport_used\": false"));
    assert!(!prov2.contains(&format!("\"selected_transport\": \"{mirror}/")));

    let (ok3, text3) = run_bootstrap(&c, &["--offline"]);
    assert!(ok3, "identical offline replay succeeds:\n{text3}");
    let prov3 = std::fs::read_to_string(&prov_path).expect("provenance rewritten again");
    assert_eq!(
        prov2, prov3,
        "identical invocation and source state produce byte-identical v3 provenance"
    );
}

#[test]
fn interrupted_and_exact_origin_unborn_checkouts_resume_in_place() {
    let c = make_constellation("interrupted-resume");
    let dest = c.root.parent().expect("destination parent");
    let mirror = c.mirror.to_str().expect("utf8 mirror").to_string();

    let marked_name = "asupersync";
    let marked_head = c
        .heads
        .iter()
        .find(|(name, _)| name == marked_name)
        .map(|(_, head)| head)
        .expect("marked fixture head");
    let marked = dest.join(marked_name);
    std::fs::create_dir(&marked).expect("marked destination");
    git(&marked, &["init", "--quiet"]);
    git(
        &marked,
        &[
            "config",
            "--local",
            "frankensim.bootstrapIncomplete",
            "true",
        ],
    );
    git(
        &marked,
        &[
            "remote",
            "add",
            "origin",
            &format!("{mirror}/{marked_name}"),
        ],
    );
    git(
        &marked,
        &["fetch", "--quiet", "--depth", "1", "origin", marked_head],
    );
    assert!(
        !git_succeeds(&marked, &["rev-parse", "HEAD"]),
        "fixture must be interrupted after fetch but before checkout"
    );

    let origin_name = "franken_networkx";
    let origin_only = dest.join(origin_name);
    std::fs::create_dir(&origin_only).expect("origin-only destination");
    git(&origin_only, &["init", "--quiet"]);
    git(
        &origin_only,
        &[
            "remote",
            "add",
            "origin",
            &format!("{mirror}/{origin_name}"),
        ],
    );
    assert!(
        !git_succeeds(
            &origin_only,
            &[
                "config",
                "--local",
                "--get",
                "frankensim.bootstrapIncomplete"
            ]
        ),
        "exact-origin fixture is intentionally unmarked"
    );

    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(ok, "interrupted bootstrap resumes:\n{text}");
    for name in [marked_name, origin_name] {
        assert!(
            text.contains(&format!("\"lib\":\"{name}\",\"state\":\"resumed\"")),
            "{name} must be reported as resumed:\n{text}"
        );
        let checkout = dest.join(name);
        let expected = c
            .heads
            .iter()
            .find(|(candidate, _)| candidate == name)
            .map(|(_, head)| head)
            .expect("fixture head");
        assert_eq!(git(&checkout, &["rev-parse", "HEAD"]), *expected);
        assert_eq!(forced_submodule_status(&checkout), "");
        assert!(
            !git_succeeds(
                &checkout,
                &[
                    "config",
                    "--local",
                    "--get",
                    "frankensim.bootstrapIncomplete"
                ]
            ),
            "successful resume must clear the incomplete marker"
        );
    }

    let (replay_ok, replay_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        replay_ok,
        "completed resume replays offline:\n{replay_text}"
    );
    for name in [marked_name, origin_name] {
        assert!(
            replay_text.contains(&format!("\"lib\":\"{name}\",\"state\":\"verified\"")),
            "{name} must become an ordinary verified checkout:\n{replay_text}"
        );
    }
}

#[test]
fn marked_wrong_head_resume_never_overwrites_an_ignored_collision() {
    {
        let c = make_constellation("standalone-ignored-resume-collision");
        let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
        let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
        assert!(seeded, "seed standalone replay:\n{seed_text}");
        let provenance_path = c
            .root
            .parent()
            .expect("constellation destination")
            .join("constellation-bootstrap.json");
        let prior_provenance = std::fs::read(&provenance_path).expect("seed provenance");
        let (target, sentinel) = install_marked_wrong_head_ignored_collision(&c);

        let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
        assert_eq!(
            std::fs::read(target.join("lib.rs")).expect("retained standalone sentinel"),
            sentinel,
            "standalone resume must retain ignored owner bytes"
        );
        assert_eq!(
            std::fs::read(&provenance_path).expect("retained standalone provenance"),
            prior_provenance,
            "failed standalone resume must retain prior provenance"
        );
        assert!(
            !ok,
            "standalone resume overwrote an ignored collision:\n{text}"
        );
        assert!(
            text.contains("lib.rs") && (text.contains("checkout") || text.contains("overwrite")),
            "standalone refusal must identify the checkout collision:\n{text}"
        );
    }

    {
        let c = make_constellation("shell-ignored-resume-collision");
        use_local_mirror_as_locked_transport(&c);
        install_shell_checkout(&c);
        commit_synthetic_root(&c);
        let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
        let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
        assert!(seeded, "seed shell replay provenance:\n{seed_text}");
        let provenance_path = c
            .root
            .parent()
            .expect("constellation destination")
            .join("constellation-bootstrap.json");
        let prior_provenance = std::fs::read(&provenance_path).expect("seed provenance");
        let (target, sentinel) = install_marked_wrong_head_ignored_collision(&c);

        let (ok, text) = run_shell_checkout(&c, None);
        assert_eq!(
            std::fs::read(target.join("lib.rs")).expect("retained shell sentinel"),
            sentinel,
            "shell resume must retain ignored owner bytes"
        );
        assert_eq!(
            std::fs::read(&provenance_path).expect("retained shell provenance"),
            prior_provenance,
            "failed shell resume must retain prior provenance"
        );
        assert!(!ok, "shell resume overwrote an ignored collision:\n{text}");
        assert!(
            text.contains("lib.rs") && (text.contains("checkout") || text.contains("overwrite")),
            "shell refusal must identify the checkout collision:\n{text}"
        );
    }
}

#[test]
fn untracked_ignore_policy_cannot_hide_bytes_in_an_exact_origin_unborn_checkout() {
    let c = make_constellation("unborn-untracked-gitignore");
    let dest = c.root.parent().expect("destination parent");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let target = dest.join("asupersync");
    std::fs::create_dir(&target).expect("unborn destination");
    git(&target, &["init", "--quiet"]);
    git(
        &target,
        &["remote", "add", "origin", &format!("{mirror}/asupersync")],
    );
    std::fs::write(target.join(".gitignore"), ".gitignore\npartial.rs\n")
        .expect("write untracked ignore policy");
    std::fs::write(target.join("partial.rs"), "unverified partial source\n")
        .expect("write ignored partial source");
    assert_eq!(
        git(&target, &["status", "--porcelain"]),
        "",
        "an untracked ignore policy must demonstrate the ordinary-status gap"
    );

    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(
        !ok,
        "an exact-origin unborn checkout with untracked ignore policy must refuse:\n{text}"
    );
    assert!(
        text.contains("untracked-ignore-policy") && text.contains(".gitignore"),
        "refusal must identify the untracked ignore authority:\n{text}"
    );
    assert!(
        !dest.join("constellation-bootstrap.json").exists(),
        "refused unborn content must not publish provenance"
    );
}

#[test]
fn unsafe_partial_destinations_are_refused_without_repurposing() {
    let c = make_constellation("partial-refusal");
    let dest = c.root.parent().expect("destination parent");
    let mirror = c.mirror.to_str().expect("utf8 mirror").to_string();

    let non_git = dest.join("asupersync");
    std::fs::create_dir(&non_git).expect("non-git destination");
    std::fs::write(non_git.join("owner-data"), "do not replace\n").expect("owner data");

    let dirty_marked = dest.join("franken_networkx");
    std::fs::create_dir(&dirty_marked).expect("dirty marked destination");
    git(&dirty_marked, &["init", "--quiet"]);
    git(
        &dirty_marked,
        &[
            "config",
            "--local",
            "frankensim.bootstrapIncomplete",
            "true",
        ],
    );
    std::fs::write(dirty_marked.join("partial.rs"), "uncommitted\n").expect("dirty partial data");

    let wrong_origin = dest.join("franken_numpy");
    std::fs::create_dir(&wrong_origin).expect("wrong-origin destination");
    git(&wrong_origin, &["init", "--quiet"]);
    git(
        &wrong_origin,
        &[
            "remote",
            "add",
            "origin",
            "https://example.invalid/not-the-lock",
        ],
    );

    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(!ok, "unsafe partial destinations must refuse:\n{text}");
    assert!(text.contains("non-empty non-git directory"), "{text}");
    assert!(
        text.contains("incomplete bootstrap with worktree or hidden-index changes"),
        "{text}"
    );
    assert!(
        text.contains("unmarked unborn checkout without the exact selected origin"),
        "{text}"
    );
    assert_eq!(
        std::fs::read_to_string(non_git.join("owner-data")).expect("owner data retained"),
        "do not replace\n"
    );
    assert!(!non_git.join(".git").exists());
    assert_eq!(
        std::fs::read_to_string(dirty_marked.join("partial.rs")).expect("dirty partial retained"),
        "uncommitted\n"
    );
    assert_eq!(
        git(&wrong_origin, &["remote", "get-url", "origin"]),
        "https://example.invalid/not-the-lock"
    );
    assert!(
        !git_succeeds(
            &wrong_origin,
            &[
                "config",
                "--local",
                "--get",
                "frankensim.bootstrapIncomplete"
            ]
        ),
        "refused wrong-origin checkout must remain unmarked"
    );
}

#[test]
fn submodule_ignore_policy_cannot_conceal_constellation_dirt() {
    for case in [
        "tracked",
        "untracked",
        "head",
        "gitlink",
        "local-exclude",
        "assume-unchanged",
        "skip-worktree",
    ] {
        let c = make_constellation(&format!("submodule-{case}"));
        use_local_mirror_as_locked_transport(&c);
        install_shell_checkout(&c);
        commit_synthetic_root(&c);
        materialize_locked_siblings_without_provenance(&c);
        let nested = initialize_sqlite_submodule(&c);
        let outer = c
            .root
            .parent()
            .expect("constellation destination")
            .join("frankensqlite");

        match case {
            "tracked" => {
                std::fs::write(nested.join("sqlite3.c"), "/* tracked nested dirt */\n")
                    .expect("write tracked nested dirt");
                assert_eq!(
                    git(&outer, &["status", "--porcelain"]),
                    "",
                    "committed ignore=dirty must demonstrate the concealed tracked fixture"
                );
            }
            "untracked" => {
                std::fs::write(nested.join("untracked.c"), "/* untracked nested dirt */\n")
                    .expect("write untracked nested dirt");
                assert_eq!(
                    git(&outer, &["status", "--porcelain"]),
                    "",
                    "committed ignore=dirty must demonstrate the concealed untracked fixture"
                );
            }
            "head" => {
                git(
                    &outer,
                    &[
                        "config",
                        "--local",
                        "submodule.legacy_sqlite_code/sqlite.ignore",
                        "all",
                    ],
                );
                git(
                    &nested,
                    &[
                        "checkout",
                        "--quiet",
                        "--detach",
                        c.sqlite_submodule_drift_head.as_str(),
                    ],
                );
                assert_eq!(
                    git(&outer, &["status", "--porcelain"]),
                    "",
                    "repository-local ignore=all must demonstrate the concealed HEAD drift"
                );
            }
            "gitlink" => {
                git(
                    &nested,
                    &[
                        "checkout",
                        "--quiet",
                        "--detach",
                        c.sqlite_submodule_drift_head.as_str(),
                    ],
                );
                git(&outer, &["add", "--", SQLITE_SUBMODULE_PATH]);
            }
            "local-exclude" => {
                write_git_admin_file(&nested, "info/exclude", "hidden-local.c\n");
                std::fs::write(nested.join("hidden-local.c"), "/* locally hidden */\n")
                    .expect("write locally excluded nested source");
                assert_eq!(
                    git(&nested, &["status", "--porcelain"]),
                    "",
                    "nested repository-local excludes must demonstrate concealment"
                );
            }
            "assume-unchanged" => {
                git(
                    &nested,
                    &["update-index", "--assume-unchanged", "sqlite3.c"],
                );
                std::fs::write(nested.join("sqlite3.c"), "/* assume-unchanged dirt */\n")
                    .expect("write assume-unchanged nested source");
                assert_eq!(
                    git(&nested, &["status", "--porcelain"]),
                    "",
                    "assume-unchanged must demonstrate nested concealment"
                );
            }
            "skip-worktree" => {
                git(&nested, &["update-index", "--skip-worktree", "sqlite3.c"]);
                std::fs::write(nested.join("sqlite3.c"), "/* skip-worktree dirt */\n")
                    .expect("write skip-worktree nested source");
                assert_eq!(
                    git(&nested, &["status", "--porcelain"]),
                    "",
                    "skip-worktree must demonstrate nested concealment"
                );
            }
            _ => unreachable!("fixed mutation matrix"),
        }

        if matches!(case, "local-exclude" | "assume-unchanged" | "skip-worktree") {
            assert_eq!(
                forced_submodule_status(&outer),
                "",
                "even forced porcelain alone must demonstrate the {case} concealment gap"
            );
        } else {
            assert!(
                !forced_submodule_status(&outer).is_empty(),
                "forced status must expose {case} nested-submodule drift"
            );
        }
        for mode in [None, Some("--verify-only"), Some("--snapshot")] {
            let (ok, text) = run_shell_checkout(&c, mode);
            assert!(
                !ok,
                "shell mode {mode:?} accepted {case} nested-submodule drift:\n{text}"
            );
        }

        let (ok, text) = run_bootstrap(&c, &["--offline"]);
        assert!(
            !ok,
            "standalone replay accepted {case} nested-submodule drift:\n{text}"
        );
        assert!(
            text.contains(SQLITE_SUBMODULE_PATH),
            "standalone refusal must identify the nested repository for {case}:\n{text}"
        );
        assert!(
            !c.root
                .parent()
                .unwrap()
                .join("constellation-bootstrap.json")
                .exists(),
            "failed {case} replay must not publish verified bootstrap provenance"
        );
    }
}

#[test]
fn recursive_cleanliness_reaches_initialized_submodules_of_submodules() {
    let c = make_constellation("recursive-submodule-cleanliness");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    materialize_locked_siblings_without_provenance(&c);
    let sqlite = initialize_sqlite_submodule(&c);
    let leaf = initialize_sqlite_leaf(&sqlite);
    let outer = c.root.parent().unwrap().join("frankensqlite");

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(ok, "clean recursive shell mode {mode:?} must pass:\n{text}");
    }
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(ok, "clean recursive standalone replay must pass:\n{text}");
    let provenance_path = c
        .root
        .parent()
        .unwrap()
        .join("constellation-bootstrap.json");
    let clean_provenance =
        std::fs::read_to_string(&provenance_path).expect("clean recursive provenance");

    std::fs::write(
        leaf.join("sqlite_leaf.c"),
        "/* concealed recursive leaf dirt */\n",
    )
    .expect("write concealed recursive leaf dirt");
    assert_eq!(
        git(&sqlite, &["status", "--porcelain"]),
        "",
        "intermediate ignore=all must conceal ordinary child status"
    );
    assert_eq!(
        git(&outer, &["status", "--porcelain"]),
        "",
        "outer ignore=dirty must conceal the recursively ignored leaf"
    );

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !ok,
            "shell mode {mode:?} accepted recursive leaf dirt:\n{text}"
        );
    }
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "standalone accepted recursive leaf dirt:\n{text}");
    assert!(
        text.contains("legacy_sqlite_code/sqlite/leaf_sqlite"),
        "recursive refusal must name the complete nested scope:\n{text}"
    );
    assert_eq!(
        std::fs::read_to_string(&provenance_path).expect("prior provenance remains"),
        clean_provenance,
        "recursive validation failure must not replace verified provenance"
    );
}

#[test]
fn snapshot_refuses_untracked_root_ignore_policy_before_digesting_hidden_bytes() {
    let c = make_constellation("root-untracked-gitignore-snapshot");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let (checkout_ok, checkout_text) = run_shell_checkout(&c, None);
    assert!(checkout_ok, "seed clean siblings:\n{checkout_text}");
    let (clean_ok, clean_snapshot) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(clean_ok, "clean snapshot must succeed:\n{clean_snapshot}");

    std::fs::write(c.root.join(".gitignore"), ".gitignore\nhidden-root.rs\n")
        .expect("write untracked root ignore policy");
    std::fs::write(
        c.root.join("hidden-root.rs"),
        "untracked hidden root bytes\n",
    )
    .expect("write hidden root bytes");
    assert_eq!(
        git(&c.root, &["status", "--porcelain"]),
        "",
        "ordinary status must demonstrate the root ignore-policy gap"
    );

    let (ignored_ok, ignored_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        !ignored_ok,
        "snapshot accepted an untracked root ignore authority:\n{ignored_text}"
    );
    assert!(
        ignored_text.contains("untracked-ignore-policy")
            || ignored_text.contains("untracked .gitignore"),
        "snapshot refusal must name the root ignore authority:\n{ignored_text}"
    );
}

#[cfg(unix)]
#[test]
fn group_only_execute_bits_do_not_change_git_semantic_mode() {
    use std::os::unix::fs::PermissionsExt as _;

    let c = make_constellation("group-execute-mode-parity");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let (checkout_ok, checkout_text) = run_shell_checkout(&c, None);
    assert!(checkout_ok, "seed clean siblings:\n{checkout_text}");
    let (baseline_ok, baseline_snapshot) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(baseline_ok, "baseline snapshot:\n{baseline_snapshot}");

    let sibling = c
        .root
        .parent()
        .expect("constellation destination")
        .join("asupersync");
    let source = sibling.join("lib.rs");
    let mut permissions = std::fs::metadata(&source)
        .expect("source metadata")
        .permissions();
    permissions.set_mode(0o654);
    std::fs::set_permissions(&source, permissions).expect("set group-only execute bit");
    assert_eq!(
        forced_submodule_status(&sibling),
        "",
        "Git tracks only the owner execute bit"
    );

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        standalone_ok,
        "standalone disagreed with Git mode semantics:\n{standalone_text}"
    );
    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            ok,
            "shell mode {mode:?} rejected group-only execute:\n{text}"
        );
        if mode == Some("--snapshot") {
            assert_eq!(
                text, baseline_snapshot,
                "group-only execute must not perturb Snapshot v3"
            );
        }
    }
}

#[test]
fn fsmonitor_index_extensions_are_explicit_refusals() {
    let c = make_constellation("fsmonitor-valid-refusal");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let sibling = c
        .root
        .parent()
        .expect("constellation destination")
        .join("asupersync");

    let persist_fsmonitor_extension = || {
        git(
            &sibling,
            &["-c", "core.fsmonitor=true", "update-index", "--fsmonitor"],
        );
        // CE_FSMONITOR_VALID is an in-memory bit. Its persistent authority is
        // the FSMN extension, so exercise the public command that marks the
        // entry and then inspect the raw index without asking Git to refresh
        // (and potentially erase) the evidence.
        git(
            &sibling,
            &[
                "-c",
                "core.fsmonitor=true",
                "update-index",
                "--fsmonitor-valid",
                "lib.rs",
            ],
        );
        let raw_index = std::fs::read(git_path(&sibling, "index")).expect("read raw Git index");
        assert!(
            raw_index.windows(4).any(|window| window == b"FSMN"),
            "fixture must persist an FSMN index extension"
        );
    };

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        persist_fsmonitor_extension();
        let (shell_ok, shell_text) = run_shell_checkout(&c, mode);
        assert!(
            !shell_ok,
            "shell mode {mode:?} accepted an FSMN index extension:\n{shell_text}"
        );
        assert!(
            {
                let diagnostic = shell_text.to_ascii_lowercase();
                diagnostic.contains("fsmonitor") || diagnostic.contains("fsmn")
            },
            "{shell_text}"
        );
    }

    persist_fsmonitor_extension();
    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone accepted an FSMN index extension:\n{standalone_text}"
    );
    assert!(
        {
            let diagnostic = standalone_text.to_ascii_lowercase();
            diagnostic.contains("fsmonitor") || diagnostic.contains("fsmn")
        },
        "{standalone_text}"
    );
}

#[test]
fn corrupted_primary_index_checksums_are_explicit_refusals() {
    let c = make_constellation("primary-index-checksum-refusal");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let sibling = c
        .root
        .parent()
        .expect("constellation destination")
        .join("asupersync");
    let index_path = git_path(&sibling, "index");
    let pristine_index = std::fs::read(&index_path).expect("read pristine Git index");
    assert!(
        pristine_index.len() > 32,
        "fixture index must contain a checksum trailer"
    );
    let persist_corrupted_checksum = || {
        let mut corrupted = pristine_index.clone();
        *corrupted.last_mut().expect("nonempty index") ^= 0x01;
        std::fs::write(&index_path, corrupted).expect("write checksum-corrupted Git index");
    };

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        persist_corrupted_checksum();
        let (shell_ok, shell_text) = run_shell_checkout(&c, mode);
        assert!(
            !shell_ok,
            "shell mode {mode:?} accepted a checksum-corrupted index:\n{shell_text}"
        );
        assert!(
            shell_text.to_ascii_lowercase().contains("checksum"),
            "shell mode {mode:?} did not identify the invalid checksum:\n{shell_text}"
        );
    }

    persist_corrupted_checksum();
    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone accepted a checksum-corrupted index:\n{standalone_text}"
    );
    assert!(
        standalone_text.to_ascii_lowercase().contains("checksum"),
        "standalone did not identify the invalid checksum:\n{standalone_text}"
    );
}

#[test]
fn forbidden_primary_index_paths_are_explicit_refusals() {
    let c = make_constellation("primary-index-path-refusal");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let sibling = c
        .root
        .parent()
        .expect("constellation destination")
        .join("asupersync");
    let index_path = git_path(&sibling, "index");
    let pristine_index = std::fs::read(&index_path).expect("read pristine Git index");
    for (case, replacement, replacement_mode, expected_diagnostic) in [
        ("hfs-dotgit", ".g\u{200c}it", None, "hfs"),
        (
            "symlink-dotgitmodules",
            ".GITMODULES",
            Some("120000"),
            "gitmodules",
        ),
    ] {
        let persist_forbidden_path = || {
            std::fs::write(&index_path, &pristine_index).expect("restore pristine Git index");
            rewrite_primary_index_entry(&sibling, "lib.rs", replacement, replacement_mode);
        };

        for mode in [None, Some("--verify-only"), Some("--snapshot")] {
            persist_forbidden_path();
            let (shell_ok, shell_text) = run_shell_checkout(&c, mode);
            assert!(
                !shell_ok,
                "shell mode {mode:?} accepted {case}:\n{shell_text}"
            );
            assert!(
                shell_text
                    .to_ascii_lowercase()
                    .contains(expected_diagnostic),
                "shell mode {mode:?} did not identify {case}:\n{shell_text}"
            );
        }

        persist_forbidden_path();
        let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
        assert!(
            !standalone_ok,
            "standalone accepted {case}:\n{standalone_text}"
        );
        assert!(
            standalone_text
                .to_ascii_lowercase()
                .contains(expected_diagnostic),
            "standalone did not identify {case}:\n{standalone_text}"
        );
    }
}

#[test]
fn unsupported_primary_index_layouts_are_explicit_refusals() {
    for (case, expected_diagnostic) in [("version-4", "version 4"), ("split-index", "split-index")]
    {
        let c = make_constellation(&format!("primary-index-{case}-refusal"));
        let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
        let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
        assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
        install_shell_checkout(&c);
        commit_synthetic_root(&c);
        let sibling = c
            .root
            .parent()
            .expect("constellation destination")
            .join("asupersync");

        let install_unsupported_layout = || {
            match case {
                "version-4" => git(&sibling, &["update-index", "--index-version", "4"]),
                "split-index" => git(&sibling, &["update-index", "--split-index"]),
                _ => unreachable!("closed fixture case set"),
            };
            let raw_index = std::fs::read(git_path(&sibling, "index")).expect("read raw Git index");
            match case {
                "version-4" => assert_eq!(
                    raw_index.get(4..8),
                    Some(&[0, 0, 0, 4][..]),
                    "fixture must persist an index-v4 header"
                ),
                "split-index" => assert!(
                    raw_index.windows(4).any(|window| window == b"link"),
                    "fixture must persist a split-index link extension"
                ),
                _ => unreachable!("closed fixture case set"),
            }
        };

        for mode in [None, Some("--verify-only"), Some("--snapshot")] {
            install_unsupported_layout();
            let (shell_ok, shell_text) = run_shell_checkout(&c, mode);
            assert!(
                !shell_ok,
                "shell mode {mode:?} accepted {case}:\n{shell_text}"
            );
            assert!(
                shell_text
                    .to_ascii_lowercase()
                    .contains(expected_diagnostic),
                "shell mode {mode:?} did not identify {case}:\n{shell_text}"
            );
        }

        install_unsupported_layout();
        let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
        assert!(
            !standalone_ok,
            "standalone accepted {case}:\n{standalone_text}"
        );
        assert!(
            standalone_text
                .to_ascii_lowercase()
                .contains(expected_diagnostic),
            "standalone did not identify {case}:\n{standalone_text}"
        );
    }
}

#[cfg(unix)]
#[test]
fn symlinked_standalone_root_refuses_before_alias_destination_access() {
    use std::os::unix::fs::symlink;

    let c = make_constellation("symlinked-standalone-root");
    let alias_parent = c.base.join("alias-destination");
    std::fs::create_dir_all(&alias_parent).expect("mkdir alias destination");
    let root_alias = alias_parent.join("frankensim");
    symlink(&c.root, &root_alias).expect("install standalone root symlink");

    let mirror = c.mirror.to_str().expect("UTF-8 mirror");
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim-bootstrap"));
    command
        .arg("--root")
        .arg(&root_alias)
        .args(["--from", mirror]);
    hermetic_git_environment(&mut command);
    let (ok, text) = run_command(&mut command);

    assert!(!ok, "standalone accepted a symlinked --root:\n{text}");
    assert!(
        text.contains(&format!(
            "error: bootstrap root {} must be an ordinary directory",
            root_alias.display()
        )),
        "standalone must identify the non-ordinary root before lock access:\n{text}"
    );
    for (name, _) in &c.heads {
        assert!(
            !alias_parent.join(name).exists(),
            "symlinked --root materialized {name} beside its alias"
        );
    }
    assert!(
        !alias_parent.join("constellation-bootstrap.json").exists(),
        "symlinked --root published provenance beside its alias"
    );
}

#[cfg(unix)]
#[test]
fn symlinked_sibling_roots_refuse_in_every_consumer() {
    use std::os::unix::fs::symlink;

    let c = make_constellation("symlinked-sibling-root");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let (checkout_ok, checkout_text) = run_shell_checkout(&c, None);
    assert!(checkout_ok, "seed clean siblings:\n{checkout_text}");

    let dest = c.root.parent().expect("constellation destination");
    let sibling = dest.join("asupersync");
    let retained = dest.join("asupersync-retained-target");
    std::fs::rename(&sibling, &retained).expect("retain sibling at an ordinary directory path");
    symlink(&retained, &sibling).expect("install sibling-root symlink");

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone accepted a symlinked sibling root:\n{standalone_text}"
    );
    for mode in [Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !ok,
            "shell mode {mode:?} accepted a symlinked sibling root:\n{text}"
        );
    }
}

#[cfg(windows)]
#[test]
fn empty_junction_destinations_refuse_without_initializing_the_external_target() {
    for consumer in ["standalone", "shell"] {
        let c = make_constellation(&format!("empty-junction-destination-{consumer}"));
        use_local_mirror_as_locked_transport(&c);
        if consumer == "shell" {
            install_shell_checkout(&c);
        }
        let external = c.base.join("external-empty-destination");
        std::fs::create_dir_all(&external).expect("mkdir external junction target");
        let destination = c.root.parent().unwrap().join("asupersync");
        let junction = Command::new("cmd")
            .args(["/D", "/C", "mklink", "/J"])
            .arg(&destination)
            .arg(&external)
            .output()
            .expect("create Windows junction");
        assert!(
            junction.status.success(),
            "mklink /J failed: {}",
            String::from_utf8_lossy(&junction.stderr)
        );

        let (ok, text) = if consumer == "shell" {
            run_shell_checkout(&c, None)
        } else {
            run_bootstrap(&c, &[])
        };
        assert!(
            !ok,
            "{consumer} accepted an empty junction destination:\n{text}"
        );
        assert!(
            !external.join(".git").exists(),
            "{consumer} initialized Git through the junction into the external target"
        );
    }
}

#[test]
fn snapshot_identity_observes_root_nested_repository_dirt() {
    let c = make_constellation("root-submodule-snapshot");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let (checkout_ok, checkout_text) = run_shell_checkout(&c, None);
    assert!(checkout_ok, "seed clean siblings:\n{checkout_text}");

    let nested_path = "root_sqlite";
    let out = hermetic_git_command()
        .arg("-C")
        .arg(&c.root)
        .args([
            "-c",
            "protocol.file.allow=always",
            "submodule",
            "add",
            "--quiet",
        ])
        .arg(&c.sqlite_submodule_upstream)
        .arg(nested_path)
        .output()
        .expect("root submodule add spawns");
    assert!(
        out.status.success(),
        "add root submodule: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    git(
        &c.root,
        &[
            "config",
            "-f",
            ".gitmodules",
            "submodule.root_sqlite.ignore",
            "all",
        ],
    );
    git(&c.root, &["add", ".gitmodules", nested_path]);
    git(&c.root, &["commit", "-qm", "pin root sqlite submodule"]);

    let (clean_ok, clean_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(clean_ok, "clean root submodule snapshot:\n{clean_text}");
    let clean_snapshot = clean_text.trim().to_string();
    let root_nested = c.root.join(nested_path);
    assert!(
        git(
            &root_nested,
            &["submodule", "status", SQLITE_LEAF_SUBMODULE_PATH]
        )
        .starts_with('-'),
        "snapshot mode must not initialize an absent nested leaf"
    );
    let (clean_replay_ok, clean_replay_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        clean_replay_ok,
        "stable uninitialized root-submodule replay:\n{clean_replay_text}"
    );
    assert_eq!(
        clean_replay_text.trim(),
        clean_snapshot,
        "an uninitialized nested gitlink must replay deterministically"
    );

    std::fs::write(root_nested.join("sqlite3.c"), "/* root nested dirt A */\n")
        .expect("write root nested dirt A");
    assert_eq!(
        git(&c.root, &["status", "--porcelain"]),
        "",
        "committed ignore=all must demonstrate the concealed root fixture"
    );
    let nested_status_a = git(&root_nested, &["status", "--porcelain"]);

    let (dirty_a_ok, dirty_a_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        dirty_a_ok,
        "dirty root content is identity-bearing rather than a sibling refusal:\n{dirty_a_text}"
    );
    let dirty_snapshot_a = dirty_a_text.trim().to_string();
    assert_eq!(dirty_snapshot_a.len(), 64);
    assert_ne!(
        clean_snapshot, dirty_snapshot_a,
        "snapshot v3 must bind initialized root-submodule dirt despite ignore=all"
    );

    std::fs::write(root_nested.join("sqlite3.c"), "/* root nested dirt B */\n")
        .expect("write root nested dirt B");
    let nested_status_b = git(&root_nested, &["status", "--porcelain"]);
    assert_eq!(
        nested_status_a, nested_status_b,
        "the fixture must prove the old HEAD-plus-status collision shape"
    );
    let (dirty_b_ok, dirty_b_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(dirty_b_ok, "dirty-B root snapshot:\n{dirty_b_text}");
    let dirty_snapshot_b = dirty_b_text.trim().to_string();
    assert_ne!(
        dirty_snapshot_a, dirty_snapshot_b,
        "snapshot v3 must bind changed bytes even when porcelain status is identical"
    );
    let (dirty_b_replay_ok, dirty_b_replay_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        dirty_b_replay_ok,
        "stable dirty-B replay:\n{dirty_b_replay_text}"
    );
    assert_eq!(dirty_b_replay_text.trim(), dirty_snapshot_b);

    std::fs::write(root_nested.join("sqlite3.c"), "/* staged index A */\n")
        .expect("write staged index A");
    git(&root_nested, &["add", "sqlite3.c"]);
    std::fs::write(root_nested.join("sqlite3.c"), "/* shared worktree */\n")
        .expect("write shared worktree after index A");
    let staged_status_a = git(&root_nested, &["status", "--porcelain"]);
    let (staged_a_ok, staged_a_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(staged_a_ok, "staged-A root snapshot:\n{staged_a_text}");

    std::fs::write(root_nested.join("sqlite3.c"), "/* staged index B */\n")
        .expect("write staged index B");
    git(&root_nested, &["add", "sqlite3.c"]);
    std::fs::write(root_nested.join("sqlite3.c"), "/* shared worktree */\n")
        .expect("restore shared worktree after index B");
    let staged_status_b = git(&root_nested, &["status", "--porcelain"]);
    assert_eq!(
        staged_status_a, staged_status_b,
        "index-only fixture must retain the same porcelain shape and worktree bytes"
    );
    let (staged_b_ok, staged_b_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(staged_b_ok, "staged-B root snapshot:\n{staged_b_text}");
    assert_ne!(
        staged_a_text.trim(),
        staged_b_text.trim(),
        "snapshot v3 must bind the nested index independently of worktree bytes"
    );

    let leaf = initialize_sqlite_leaf(&root_nested);
    let (leaf_clean_ok, leaf_clean_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        leaf_clean_ok,
        "initialized leaf snapshot:\n{leaf_clean_text}"
    );
    std::fs::write(leaf.join("sqlite_leaf.c"), "/* recursive leaf dirt A */\n")
        .expect("write recursive leaf dirt A");
    let concealed_intermediate_status_a = git(&root_nested, &["status", "--porcelain"]);
    assert_eq!(
        git(&c.root, &["status", "--porcelain"]),
        "",
        "root ignore=all must conceal the depth-two dirt"
    );
    let (leaf_a_ok, leaf_a_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(leaf_a_ok, "recursive leaf-A snapshot:\n{leaf_a_text}");
    std::fs::write(leaf.join("sqlite_leaf.c"), "/* recursive leaf dirt B */\n")
        .expect("write recursive leaf dirt B");
    assert_eq!(
        git(&root_nested, &["status", "--porcelain"]),
        concealed_intermediate_status_a,
        "intermediate ignore=all must preserve the same concealed status shape"
    );
    let (leaf_b_ok, leaf_b_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(leaf_b_ok, "recursive leaf-B snapshot:\n{leaf_b_text}");
    assert_ne!(
        leaf_a_text.trim(),
        leaf_b_text.trim(),
        "snapshot v3 must recursively bind depth-two worktree bytes"
    );
}

#[test]
fn repository_local_excludes_cannot_hide_untracked_source() {
    let c = make_constellation("hidden-untracked");
    let mirror = c.mirror.to_str().expect("utf8").to_string();
    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(ok, "seed the source cache:\n{text}");
    let provenance_path = c
        .root
        .parent()
        .unwrap()
        .join("constellation-bootstrap.json");
    let admitted_provenance = std::fs::read_to_string(&provenance_path).expect("seed provenance");

    let sibling = c.root.parent().unwrap().join("asupersync");
    write_git_admin_file(&sibling, "info/exclude", "hidden-source.rs\n");
    std::fs::write(sibling.join("hidden-source.rs"), "pub fn hidden() {}\n")
        .expect("write hidden source");
    assert_eq!(
        git(&sibling, &["status", "--porcelain"]),
        "",
        "ordinary porcelain must demonstrate the hidden-file fixture"
    );

    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "repository-local excludes must not hide dirt:\n{text}");
    assert!(
        text.contains("DIRTY") && text.contains("asupersync"),
        "hidden source refusal is explicit:\n{text}"
    );
    assert_eq!(
        std::fs::read_to_string(&provenance_path).expect("prior provenance remains readable"),
        admitted_provenance,
        "failed replay must not replace an earlier verified provenance document"
    );
}

#[test]
fn inherited_pathspec_modes_cannot_hide_nested_untracked_ignore_policy() {
    let c = make_constellation("hostile-pathspec-environment");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let (seed_ok, seed_text) = run_shell_checkout(&c, None);
    assert!(seed_ok, "seed shell checkout failed:\n{seed_text}");

    let sibling = c.root.parent().unwrap().join("asupersync");
    let hidden = sibling.join("nested-ignore-policy");
    std::fs::create_dir_all(&hidden).expect("mkdir nested ignore policy");
    std::fs::write(hidden.join(".gitignore"), "*\n").expect("write self-ignoring policy");
    std::fs::write(hidden.join("hidden-source.rs"), "pub fn concealed() {}\n")
        .expect("write concealed source");

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let mut command = shell_checkout_command(&c, mode);
        command.env("GIT_LITERAL_PATHSPECS", "1");
        let (ok, text) = run_command(&mut command);
        assert!(
            !ok,
            "shell mode {mode:?} accepted a pathspec-hidden ignore authority:\n{text}"
        );
        assert!(
            text.contains("untracked-ignore-policy") || text.contains("untracked .gitignore"),
            "shell mode {mode:?} must identify the hidden ignore authority:\n{text}"
        );
    }

    let mut command = bootstrap_command(&c, &["--offline"]);
    command.env("GIT_LITERAL_PATHSPECS", "1");
    let (ok, text) = run_command(&mut command);
    assert!(
        !ok,
        "standalone accepted a pathspec-hidden ignore authority:\n{text}"
    );
    assert!(
        text.contains("untracked-ignore-policy") && text.contains(".gitignore"),
        "standalone must identify the hidden ignore authority:\n{text}"
    );
    assert!(
        !c.root
            .parent()
            .unwrap()
            .join("constellation-bootstrap.json")
            .exists(),
        "pathspec-hidden failure must not publish provenance"
    );
}

#[test]
fn index_flags_cannot_hide_modified_tracked_source() {
    let c = make_constellation("hidden-tracked");
    let mirror = c.mirror.to_str().expect("utf8").to_string();
    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(ok, "seed the source cache:\n{text}");

    let sibling = c.root.parent().unwrap().join("asupersync");
    git(&sibling, &["update-index", "--assume-unchanged", "lib.rs"]);
    std::fs::write(sibling.join("lib.rs"), "pub fn concealed_tamper() {}\n")
        .expect("modify assume-unchanged source");
    assert_eq!(
        git(&sibling, &["status", "--porcelain"]),
        "",
        "ordinary porcelain must demonstrate the index-flag fixture"
    );

    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "index flags must not hide tracked dirt:\n{text}");
    assert!(
        text.contains("DIRTY") && text.contains("asupersync"),
        "hidden tracked-source refusal is explicit:\n{text}"
    );
}

#[cfg(unix)]
#[test]
fn tracked_filesystem_aliases_refuse_in_shell_standalone_and_snapshot_paths() {
    let c = make_constellation("tracked-filesystem-alias");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    materialize_locked_siblings_without_provenance(&c);

    let sibling = c.root.parent().unwrap().join("asupersync");
    let sibling_oid = git(&sibling, &["rev-parse", ":lib.rs"]);
    std::fs::hard_link(sibling.join("lib.rs"), sibling.join("lib-alias.rs"))
        .expect("create tracked hard-link alias");
    git(
        &sibling,
        &[
            "update-index",
            "--add",
            "--cacheinfo",
            "100644",
            &sibling_oid,
            "lib-alias.rs",
        ],
    );

    for mode in [None, Some("--verify-only")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !ok,
            "shell mode {mode:?} accepted two index paths for one file:\n{text}"
        );
        assert!(
            text.contains("filesystem identity") && text.contains("hard-link aliases"),
            "shell mode {mode:?} must identify the filesystem alias:\n{text}"
        );
    }

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone accepted two index paths for one file:\n{standalone_text}"
    );
    assert!(
        standalone_text.contains("filesystem identity")
            && standalone_text.contains("hard-link aliases"),
        "standalone must identify the filesystem alias:\n{standalone_text}"
    );

    let root_oid = git(&c.root, &["rev-parse", ":constellation.lock"]);
    std::fs::hard_link(
        c.root.join("constellation.lock"),
        c.root.join("constellation-lock-alias"),
    )
    .expect("create root hard-link alias");
    git(
        &c.root,
        &[
            "update-index",
            "--add",
            "--cacheinfo",
            "100644",
            &root_oid,
            "constellation-lock-alias",
        ],
    );
    let (snapshot_ok, snapshot_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        !snapshot_ok,
        "snapshot accepted two index paths for one file:\n{snapshot_text}"
    );
    assert!(
        snapshot_text.contains("filesystem identity")
            && snapshot_text.contains("hard-link aliases"),
        "snapshot must identify the filesystem alias:\n{snapshot_text}"
    );
}

#[cfg(unix)]
#[test]
fn tracked_ancestor_symlinks_refuse_before_external_bytes_are_hashed() {
    use std::os::unix::fs::symlink;

    let c = make_constellation("tracked-ancestor-symlink");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    materialize_locked_siblings_without_provenance(&c);

    let sibling = c.root.parent().unwrap().join("asupersync");
    let sibling_oid = git(&sibling, &["rev-parse", ":lib.rs"]);
    let outside_sibling = c.base.join("outside-sibling-source");
    std::fs::create_dir_all(&outside_sibling).expect("mkdir external sibling source");
    std::fs::copy(sibling.join("lib.rs"), outside_sibling.join("lib.rs"))
        .expect("copy external sibling bytes");
    symlink(&outside_sibling, sibling.join("redirected"))
        .expect("install tracked ancestor symlink");
    git(
        &sibling,
        &[
            "update-index",
            "--add",
            "--cacheinfo",
            "100644",
            &sibling_oid,
            "redirected/lib.rs",
        ],
    );

    for mode in [None, Some("--verify-only")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !ok,
            "shell mode {mode:?} hashed through a tracked ancestor symlink:\n{text}"
        );
        assert!(
            text.contains("tracked index prefix") && text.contains("ancestor"),
            "shell mode {mode:?} must identify the ancestor redirection:\n{text}"
        );
    }

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone hashed through a tracked ancestor symlink:\n{standalone_text}"
    );
    assert!(
        standalone_text.contains("tracked index prefix") && standalone_text.contains("ancestor"),
        "standalone must identify the ancestor redirection:\n{standalone_text}"
    );

    let root_oid = git(&c.root, &["rev-parse", ":constellation.lock"]);
    let outside_root = c.base.join("outside-root-source");
    std::fs::create_dir_all(&outside_root).expect("mkdir external root source");
    std::fs::copy(
        c.root.join("constellation.lock"),
        outside_root.join("constellation.lock"),
    )
    .expect("copy external root bytes");
    symlink(&outside_root, c.root.join("redirected")).expect("install root ancestor symlink");
    git(
        &c.root,
        &[
            "update-index",
            "--add",
            "--cacheinfo",
            "100644",
            &root_oid,
            "redirected/constellation.lock",
        ],
    );
    let (snapshot_ok, snapshot_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        !snapshot_ok,
        "snapshot hashed through a tracked ancestor symlink:\n{snapshot_text}"
    );
    assert!(
        snapshot_text.contains("tracked index prefix") && snapshot_text.contains("ancestor"),
        "snapshot must identify the ancestor redirection:\n{snapshot_text}"
    );
}

#[cfg(unix)]
#[test]
fn clean_filter_cannot_normalize_modified_raw_source_into_a_false_clean_result() {
    let c = make_constellation("raw-filter-bypass");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    materialize_locked_siblings_without_provenance(&c);
    let sibling = c.root.parent().unwrap().join("asupersync");
    write_git_admin_file(&sibling, "info/attributes", "lib.rs filter=conceal\n");
    let filter_marker = c.base.join("clean-filter-executed");
    let hostile_clean = format!(
        "printf 'executed\\n' > '{}'; printf 'pub fn asupersync_fixture() {{}}\\n'",
        filter_marker.display()
    );
    git(
        &sibling,
        &["config", "--local", "filter.conceal.clean", &hostile_clean],
    );
    git(
        &sibling,
        &["config", "--local", "filter.conceal.smudge", "cat"],
    );
    std::fs::write(
        sibling.join("lib.rs"),
        "pub fn raw_bytes_are_not_the_pinned_source() {}\n",
    )
    .expect("write filter-concealed raw source");
    assert!(
        !filter_marker.exists(),
        "installing a filter authority must not execute it"
    );

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !filter_marker.exists(),
            "shell mode {mode:?} executed a refused clean-filter helper"
        );
        assert!(
            !ok,
            "shell mode {mode:?} accepted filter-concealed raw bytes:\n{text}"
        );
        assert!(
            text.contains("filter") && text.contains("authority"),
            "shell refusal must identify the executable filter authority:\n{text}"
        );
    }
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !filter_marker.exists(),
        "standalone verification executed a refused clean-filter helper"
    );
    assert!(
        !ok,
        "standalone accepted filter-concealed raw source:\n{text}"
    );
    assert!(
        text.contains("filter") && text.contains("authority"),
        "standalone refusal must identify the executable filter authority:\n{text}"
    );
    assert!(
        !c.root
            .parent()
            .unwrap()
            .join("constellation-bootstrap.json")
            .exists(),
        "filter-concealed failure must not publish provenance"
    );
}

#[cfg(unix)]
#[test]
fn worktree_config_cannot_install_an_executable_filter_authority() {
    let c = make_constellation("worktree-filter-authority");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    materialize_locked_siblings_without_provenance(&c);

    let sibling = c.root.parent().unwrap().join("asupersync");
    let filter_marker = c.base.join("worktree-filter-executed");
    let hostile_process = format!("printf 'executed\\n' > '{}'; cat", filter_marker.display());
    git(
        &sibling,
        &["config", "--local", "extensions.worktreeConfig", "true"],
    );
    git(
        &sibling,
        &[
            "config",
            "--worktree",
            "filter.worktree-hostile.process",
            &hostile_process,
        ],
    );
    write_git_admin_file(
        &sibling,
        "info/attributes",
        "lib.rs filter=worktree-hostile\n",
    );

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !filter_marker.exists(),
            "shell mode {mode:?} executed a refused worktree-config filter helper"
        );
        assert!(
            !ok,
            "shell mode {mode:?} accepted executable worktree config:\n{text}"
        );
        assert!(
            text.contains("--worktree") && text.contains("filter.worktree-hostile.process"),
            "shell refusal must identify the worktree authority:\n{text}"
        );
    }

    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !filter_marker.exists(),
        "standalone verification executed a refused worktree-config filter helper"
    );
    assert!(
        !ok,
        "standalone accepted executable worktree config:\n{text}"
    );
    assert!(
        text.contains("--worktree") && text.contains("filter.worktree-hostile.process"),
        "standalone refusal must identify the worktree authority:\n{text}"
    );
}

#[cfg(unix)]
#[test]
fn parent_status_never_executes_nested_filter_config_before_child_admission() {
    let c = make_constellation("nested-filter-admission-order");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    materialize_locked_siblings_without_provenance(&c);
    let nested = initialize_sqlite_submodule(&c);
    let filter_marker = c.base.join("nested-clean-filter-executed");
    let hostile_clean = format!(
        "printf 'executed\\n' > '{}'; printf '/* pinned sqlite fixture */\\n'",
        filter_marker.display()
    );
    write_git_admin_file(
        &nested,
        "info/attributes",
        "sqlite3.c filter=nested-conceal\n",
    );
    git(
        &nested,
        &[
            "config",
            "--local",
            "filter.nested-conceal.clean",
            &hostile_clean,
        ],
    );
    git(
        &nested,
        &["config", "--local", "filter.nested-conceal.smudge", "cat"],
    );
    std::fs::write(
        nested.join("sqlite3.c"),
        "/* raw nested source is not pinned */\n",
    )
    .expect("modify nested raw source");

    for mode in [None, Some("--verify-only"), Some("--snapshot")] {
        let (ok, text) = run_shell_checkout(&c, mode);
        assert!(
            !filter_marker.exists(),
            "shell mode {mode:?} executed nested filter config before admission"
        );
        assert!(
            !ok,
            "shell mode {mode:?} accepted nested executable config:\n{text}"
        );
        assert!(
            text.contains("--local") && text.contains("filter.nested-conceal.clean"),
            "shell refusal must identify the nested filter authority:\n{text}"
        );
    }

    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !filter_marker.exists(),
        "standalone executed nested filter config before admission"
    );
    assert!(!ok, "standalone accepted nested executable config:\n{text}");
    assert!(
        text.contains(SQLITE_SUBMODULE_PATH) && text.contains("filter.nested-conceal.clean"),
        "standalone refusal must identify the nested filter authority:\n{text}"
    );
}

#[cfg(unix)]
#[test]
fn unlisted_transport_helpers_are_never_executed() {
    use std::os::unix::fs::PermissionsExt;

    for consumer in ["standalone", "shell"] {
        let c = make_constellation(&format!("unlisted-transport-{consumer}"));
        if consumer == "shell" {
            install_shell_checkout(&c);
        }
        let lock_path = c.root.join("constellation.lock");
        let original_lock = std::fs::read_to_string(&lock_path).expect("fixture lock");
        let hostile_lock =
            original_lock.replace("\"remote\": \"no-remote\"", "\"remote\": \"evil::payload\"");
        assert_ne!(
            hostile_lock, original_lock,
            "fixture must replace transports"
        );
        std::fs::write(&lock_path, hostile_lock).expect("install hostile transports");

        let helper_dir = c.base.join("hostile-transport-bin");
        std::fs::create_dir_all(&helper_dir).expect("mkdir hostile helper directory");
        let marker = c.base.join("git-remote-evil-executed");
        let helper = helper_dir.join("git-remote-evil");
        std::fs::write(
            &helper,
            format!(
                "#!/bin/sh\nprintf 'executed\\n' > '{}'\nexit 1\n",
                marker.display()
            ),
        )
        .expect("write hostile remote helper");
        let mut permissions = std::fs::metadata(&helper)
            .expect("hostile helper metadata")
            .permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(&helper, permissions).expect("make hostile helper executable");

        let mut search_path = vec![helper_dir];
        if let Some(inherited) = std::env::var_os("PATH") {
            search_path.extend(std::env::split_paths(&inherited));
        }
        let search_path = std::env::join_paths(search_path).expect("join hostile helper PATH");
        let mut command = if consumer == "shell" {
            shell_checkout_command(&c, None)
        } else {
            bootstrap_command(&c, &[])
        };
        command.env("PATH", search_path);
        let (ok, text) = run_command(&mut command);

        assert!(
            !ok,
            "{consumer} accepted an unlisted Git transport:\n{text}"
        );
        assert!(
            !marker.exists(),
            "{consumer} executed git-remote-evil before refusing the transport"
        );
        assert!(
            text.contains("evil") && text.contains("not allowed"),
            "{consumer} refusal must identify the denied transport:\n{text}"
        );
        assert!(
            !c.root
                .parent()
                .unwrap()
                .join("constellation-bootstrap.json")
                .exists(),
            "{consumer} transport refusal must not publish provenance"
        );
    }
}

#[test]
fn local_and_worktree_reference_storage_authorities_refuse_before_bootstrap_mutation() {
    for config_scope in ["--local", "--worktree"] {
        let c = make_constellation(&format!(
            "reference-storage-authority-{}",
            config_scope.trim_start_matches('-')
        ));
        use_local_mirror_as_locked_transport(&c);
        install_shell_checkout(&c);
        let target = c.root.parent().unwrap().join("asupersync");
        std::fs::create_dir_all(&target).expect("mkdir incomplete destination");
        git(&target, &["init", "-q", "-b", "main"]);
        let exact_origin = c.mirror.join("asupersync").display().to_string();
        git(&target, &["remote", "add", "origin", &exact_origin]);
        git(
            &target,
            &["config", "--local", "core.repositoryFormatVersion", "1"],
        );
        if config_scope == "--worktree" {
            git(
                &target,
                &["config", "--local", "extensions.worktreeConfig", "true"],
            );
        }
        git(
            &target,
            &["config", config_scope, "extensions.refStorage", "files"],
        );

        for consumer in ["shell", "standalone"] {
            let mut command = if consumer == "shell" {
                shell_checkout_command(&c, None)
            } else {
                bootstrap_command(&c, &[])
            };
            let (ok, text) = run_command(&mut command);
            assert!(
                !ok,
                "{consumer} accepted {config_scope} reference-storage authority:\n{text}"
            );
            assert!(
                text.contains(config_scope)
                    && text.to_ascii_lowercase().contains("extensions.refstorage"),
                "{consumer} refusal must identify {config_scope} reference storage:\n{text}"
            );
            assert_eq!(
                git(
                    &target,
                    &["config", config_scope, "--get", "extensions.refStorage"]
                ),
                "files",
                "refusal must leave the owner-managed reference authority unchanged"
            );
        }
    }
}

#[cfg(unix)]
#[test]
fn local_and_worktree_remote_vcs_authorities_never_execute_allowlisted_helpers() {
    use std::os::unix::fs::PermissionsExt;

    for config_scope in ["--local", "--worktree"] {
        let c = make_constellation(&format!(
            "remote-vcs-authority-{}",
            config_scope.trim_start_matches('-')
        ));
        use_local_mirror_as_locked_transport(&c);
        install_shell_checkout(&c);
        let target = c.root.parent().unwrap().join("asupersync");
        std::fs::create_dir_all(&target).expect("mkdir incomplete destination");
        git(&target, &["init", "-q", "-b", "main"]);
        let exact_origin = c.mirror.join("asupersync").display().to_string();
        git(&target, &["remote", "add", "origin", &exact_origin]);
        if config_scope == "--worktree" {
            git(
                &target,
                &["config", "--local", "extensions.worktreeConfig", "true"],
            );
        }
        git(
            &target,
            &["config", config_scope, "remote.origin.vcs", "ssh"],
        );

        let helper_dir = c.base.join("hostile-ssh-helper-bin");
        std::fs::create_dir_all(&helper_dir).expect("mkdir hostile ssh helper directory");
        let marker = c.base.join("git-remote-ssh-executed");
        let helper = helper_dir.join("git-remote-ssh");
        std::fs::write(
            &helper,
            format!(
                "#!/bin/sh\nprintf 'executed\\n' > '{}'\nexit 1\n",
                marker.display()
            ),
        )
        .expect("write hostile ssh remote helper");
        let mut permissions = std::fs::metadata(&helper)
            .expect("hostile ssh helper metadata")
            .permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(&helper, permissions).expect("make hostile ssh helper executable");
        let mut search_path = vec![helper_dir];
        if let Some(inherited) = std::env::var_os("PATH") {
            search_path.extend(std::env::split_paths(&inherited));
        }
        let search_path = std::env::join_paths(search_path).expect("join hostile ssh helper PATH");

        for consumer in ["shell", "standalone"] {
            let mut command = if consumer == "shell" {
                shell_checkout_command(&c, None)
            } else {
                bootstrap_command(&c, &[])
            };
            command.env("PATH", &search_path);
            let (ok, text) = run_command(&mut command);
            assert!(
                !ok,
                "{consumer} accepted {config_scope} remote VCS authority:\n{text}"
            );
            assert!(
                !marker.exists(),
                "{consumer} executed git-remote-ssh from {config_scope} config"
            );
            assert!(
                text.contains(config_scope) && text.contains("remote.origin.vcs"),
                "{consumer} refusal must identify the {config_scope} remote VCS authority:\n{text}"
            );
        }
    }
}

#[cfg(unix)]
#[test]
fn local_and_worktree_maintenance_hooks_refuse_before_fetch() {
    use std::os::unix::fs::PermissionsExt;

    for config_scope in ["--local", "--worktree"] {
        let c = make_constellation(&format!(
            "maintenance-hook-authority-{}",
            config_scope.trim_start_matches('-')
        ));
        use_local_mirror_as_locked_transport(&c);
        install_shell_checkout(&c);
        let target = c.root.parent().unwrap().join("asupersync");
        std::fs::create_dir_all(&target).expect("mkdir incomplete destination");
        git(&target, &["init", "-q", "-b", "main"]);
        let exact_origin = c.mirror.join("asupersync").display().to_string();
        git(&target, &["remote", "add", "origin", &exact_origin]);
        if config_scope == "--worktree" {
            git(
                &target,
                &["config", "--local", "extensions.worktreeConfig", "true"],
            );
        }

        let marker = c.base.join("recent-objects-hook-executed");
        let hook = c.base.join("recent-objects-hook");
        std::fs::write(
            &hook,
            format!(
                "#!/bin/sh\nprintf 'executed\\n' > '{}'\nexit 0\n",
                marker.display()
            ),
        )
        .expect("write recent-objects hook");
        let mut permissions = std::fs::metadata(&hook)
            .expect("recent-objects hook metadata")
            .permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(&hook, permissions).expect("make maintenance hook executable");
        git(
            &target,
            &[
                "config",
                config_scope,
                "gc.recentObjectsHook",
                hook.to_str().expect("UTF-8 fixture hook path"),
            ],
        );

        for consumer in ["shell", "standalone"] {
            let mut command = if consumer == "shell" {
                shell_checkout_command(&c, None)
            } else {
                bootstrap_command(&c, &[])
            };
            let (ok, text) = run_command(&mut command);
            assert!(
                !ok,
                "{consumer} accepted {config_scope} maintenance hook authority:\n{text}"
            );
            assert!(
                !marker.exists(),
                "{consumer} executed gc.recentObjectsHook from {config_scope} config"
            );
            assert!(
                text.contains(config_scope)
                    && text.to_ascii_lowercase().contains("gc.recentobjectshook"),
                "{consumer} refusal must identify the {config_scope} maintenance hook:\n{text}"
            );
        }

        git(
            &target,
            &[
                "config",
                config_scope,
                "--unset-all",
                "gc.recentObjectsHook",
            ],
        );
        git(
            &target,
            &[
                "config",
                config_scope,
                "fetch.bundleURI",
                "file:///unrecorded-bootstrap.bundle",
            ],
        );
        for consumer in ["shell", "standalone"] {
            let mut command = if consumer == "shell" {
                shell_checkout_command(&c, None)
            } else {
                bootstrap_command(&c, &[])
            };
            let (ok, text) = run_command(&mut command);
            assert!(
                !ok,
                "{consumer} accepted {config_scope} bundle-URI authority:\n{text}"
            );
            assert!(
                text.contains(config_scope)
                    && text.to_ascii_lowercase().contains("fetch.bundleuri"),
                "{consumer} refusal must identify the {config_scope} bundle URI:\n{text}"
            );
        }

        git(
            &target,
            &["config", config_scope, "--unset-all", "fetch.bundleURI"],
        );
        git(
            &target,
            &["config", config_scope, "fetch.uriProtocols", "http"],
        );
        for consumer in ["shell", "standalone"] {
            let mut command = if consumer == "shell" {
                shell_checkout_command(&c, None)
            } else {
                bootstrap_command(&c, &[])
            };
            let (ok, text) = run_command(&mut command);
            assert!(
                !ok,
                "{consumer} accepted {config_scope} packfile-URI authority:\n{text}"
            );
            assert!(
                text.contains(config_scope)
                    && text.to_ascii_lowercase().contains("fetch.uriprotocols"),
                "{consumer} refusal must identify the {config_scope} packfile URI protocols:\n{text}"
            );
        }
    }
}

#[test]
fn bootstrap_never_recurses_into_owner_initialized_submodules() {
    for consumer in ["standalone", "shell"] {
        let c = make_constellation(&format!("no-submodule-recursion-{consumer}"));
        use_local_mirror_as_locked_transport(&c);
        if consumer == "shell" {
            install_shell_checkout(&c);
        }
        materialize_locked_siblings_without_provenance(&c);
        let nested = initialize_sqlite_submodule(&c);
        let outer = c.root.parent().unwrap().join("frankensqlite");
        git(
            &nested,
            &["fetch", "--quiet", "origin", &c.sqlite_submodule_drift_head],
        );
        git(
            &nested,
            &[
                "checkout",
                "--quiet",
                "--detach",
                &c.sqlite_submodule_drift_head,
            ],
        );
        git(&outer, &["config", "user.email", "drill@frankensim.test"]);
        git(&outer, &["config", "user.name", "bootstrap drill"]);
        git(&outer, &["add", SQLITE_SUBMODULE_PATH]);
        git(&outer, &["commit", "-qm", "marked wrong-head gitlink"]);
        git(
            &outer,
            &[
                "config",
                "--local",
                "frankensim.bootstrapIncomplete",
                "true",
            ],
        );
        git(
            &outer,
            &["config", "--local", "fetch.recurseSubmodules", "true"],
        );
        git(&outer, &["config", "--local", "submodule.recurse", "true"]);
        assert_eq!(
            forced_submodule_status(&outer),
            "",
            "wrong-head fixture must be internally clean before bootstrap"
        );

        let nested_before = git(&nested, &["rev-parse", "HEAD"]);
        let mut command = if consumer == "shell" {
            shell_checkout_command(&c, None)
        } else {
            bootstrap_command(&c, &[])
        };
        let (ok, text) = run_command(&mut command);

        assert!(
            !ok,
            "{consumer} unexpectedly repaired an owner-managed nested checkout:\n{text}"
        );
        assert!(
            text.contains(SQLITE_SUBMODULE_PATH) || text.contains("pinned but not clean"),
            "{consumer} must report the unchanged nested drift after the outer checkout:\n{text}"
        );
        assert_eq!(
            git(&nested, &["rev-parse", "HEAD"]),
            nested_before,
            "{consumer} must not update an initialized nested worktree"
        );
    }
}

#[cfg(unix)]
#[test]
fn shell_refuses_a_destination_when_enumeration_fails_before_git_init() {
    use std::os::unix::fs::PermissionsExt;

    let c = make_constellation("destination-enumeration-failure");
    use_local_mirror_as_locked_transport(&c);
    install_shell_checkout(&c);
    let target = c.root.parent().unwrap().join("asupersync");
    std::fs::create_dir_all(&target).expect("mkdir destination");
    let sentinel = target.join("sentinel");
    std::fs::write(&sentinel, "must remain byte-identical\n").expect("write sentinel");

    let helper_dir = c.base.join("failing-find-bin");
    std::fs::create_dir_all(&helper_dir).expect("mkdir failing find directory");
    let find_marker = c.base.join("failing-find-invoked");
    let find_helper = helper_dir.join("find");
    std::fs::write(
        &find_helper,
        format!(
            "#!/bin/sh\nprintf 'invoked\\n' > '{}'\nexit 73\n",
            find_marker.display()
        ),
    )
    .expect("write failing find helper");
    let mut permissions = std::fs::metadata(&find_helper)
        .expect("failing find metadata")
        .permissions();
    permissions.set_mode(0o755);
    std::fs::set_permissions(&find_helper, permissions).expect("make failing find executable");
    let mut search_path = vec![helper_dir];
    if let Some(inherited) = std::env::var_os("PATH") {
        search_path.extend(std::env::split_paths(&inherited));
    }
    let search_path = std::env::join_paths(search_path).expect("join failing find PATH");

    let mut command = shell_checkout_command(&c, None);
    command.env("PATH", search_path);
    let (ok, text) = run_command(&mut command);

    assert!(
        !ok,
        "shell accepted a destination it could not enumerate:\n{text}"
    );
    assert!(
        find_marker.exists(),
        "the enumeration failure was not injected"
    );
    assert!(
        text.contains("cannot enumerate bootstrap destination")
            && text.contains("unreadable non-git directory"),
        "shell refusal must identify failed destination enumeration:\n{text}"
    );
    assert_eq!(
        std::fs::read_to_string(&sentinel).expect("sentinel remains readable"),
        "must remain byte-identical\n"
    );
    assert!(
        !target.join(".git").exists(),
        "shell must refuse before initializing an unreadable destination"
    );
}

#[cfg(unix)]
#[test]
fn fresh_clone_ignores_hostile_global_checkout_hooks() {
    use std::os::unix::fs::PermissionsExt;

    let c = make_constellation("post-clone-dirty");
    let hooks = c.base.join("hooks");
    std::fs::create_dir_all(&hooks).expect("mkdir hooks");
    let post_checkout = hooks.join("post-checkout");
    let hook_marker = c.base.join("post-checkout-hook-executed");
    std::fs::write(
        &post_checkout,
        format!(
            "#!/bin/sh\nprintf 'executed\\n' > '{}'\nprintf 'checkout-time mutation\\n' > lib.rs\n",
            hook_marker.display()
        ),
    )
    .expect("write post-checkout hook");
    let mut permissions = std::fs::metadata(&post_checkout)
        .expect("hook metadata")
        .permissions();
    permissions.set_mode(0o755);
    std::fs::set_permissions(&post_checkout, permissions).expect("make hook executable");
    let global_config = c.base.join("post-checkout.gitconfig");
    std::fs::write(
        &global_config,
        format!("[core]\n\thooksPath = {}\n", hooks.display()),
    )
    .expect("write checkout-hook config");

    let mirror = c.mirror.to_str().expect("utf8");
    let mut command = bootstrap_command(&c, &["--from", mirror]);
    command
        .env("GIT_CONFIG_GLOBAL", &global_config)
        .env("GIT_CONFIG_NOSYSTEM", "1");
    let (ok, text) = run_command(&mut command);

    assert!(
        !hook_marker.exists(),
        "hostile inherited post-checkout hook must never execute"
    );
    assert!(
        ok,
        "bootstrap must neutralize inherited global checkout hooks:\n{text}"
    );
    assert_eq!(
        std::fs::read_to_string(c.root.parent().unwrap().join("asupersync/lib.rs"))
            .expect("pinned source remains readable"),
        "pub fn asupersync_fixture() {}\n",
        "fresh checkout must retain the exact pinned source bytes"
    );
    assert!(
        c.root
            .parent()
            .unwrap()
            .join("constellation-bootstrap.json")
            .exists(),
        "successful hermetic bootstrap must publish provenance"
    );
}

#[test]
fn inherited_git_redirection_cannot_spoof_a_pinned_sibling() {
    let c = make_constellation("hostile-git-environment");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
    install_shell_checkout(&c);
    commit_synthetic_root(&c);
    let (baseline_ok, baseline_text) = run_shell_checkout(&c, Some("--snapshot"));
    assert!(
        baseline_ok,
        "hostile-environment fixture must begin from an admissible root:\n{baseline_text}"
    );

    let dest = c.root.parent().expect("constellation destination");
    let target = dest.join("asupersync");
    git(&target, &["config", "user.email", "drill@frankensim.test"]);
    git(&target, &["config", "user.name", "bootstrap drill"]);
    std::fs::write(target.join("drift.rs"), "pub fn drift() {}\n").expect("write drift");
    git(&target, &["add", "drift.rs"]);
    git(&target, &["commit", "-qm", "wrong clean head"]);
    let wrong_head = git(&target, &["rev-parse", "HEAD"]);

    let prior_provenance = std::fs::read(dest.join("constellation-bootstrap.json"))
        .expect("seed provenance remains available");
    let decoy_git_dir = c.mirror.join("asupersync");
    let mut standalone = bootstrap_command(&c, &["--offline"]);
    standalone
        .env("GIT_DIR", &decoy_git_dir)
        .env("GIT_WORK_TREE", &target)
        .env("GIT_NAMESPACE", "hostile-namespace");
    let (standalone_ok, standalone_text) = run_command(&mut standalone);
    assert!(
        !standalone_ok,
        "standalone accepted a wrong clean head through inherited Git redirection:\n{standalone_text}"
    );
    assert!(
        standalone_text.contains("asupersync") && standalone_text.contains(&wrong_head),
        "the refusal must prove the real target repository was observed:\n{standalone_text}"
    );

    for mode in [Some("--verify-only"), Some("--snapshot")] {
        let mut shell = shell_checkout_command(&c, mode);
        shell
            .env("GIT_DIR", &decoy_git_dir)
            .env("GIT_WORK_TREE", &target)
            .env("GIT_NAMESPACE", "hostile-namespace");
        let (shell_ok, shell_text) = run_command(&mut shell);
        assert!(
            !shell_ok,
            "shell mode {mode:?} accepted a wrong clean head through inherited Git redirection:\n{shell_text}"
        );
        assert!(
            shell_text.contains("asupersync") && shell_text.contains(&wrong_head),
            "the refusal must prove the real target repository was observed:\n{shell_text}"
        );
    }
    assert_eq!(
        std::fs::read(dest.join("constellation-bootstrap.json"))
            .expect("refusal retains prior provenance"),
        prior_provenance,
    );
}

#[test]
fn replacement_objects_cannot_redefine_locked_source() {
    let c = make_constellation("replacement-object-refusal");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
    install_shell_checkout(&c);

    let dest = c.root.parent().expect("constellation destination");
    let target = dest.join("asupersync");
    let pinned = c
        .heads
        .iter()
        .find(|(name, _)| name == "asupersync")
        .map(|(_, head)| head.as_str())
        .expect("asupersync pin");
    git(&target, &["config", "user.email", "drill@frankensim.test"]);
    git(&target, &["config", "user.name", "bootstrap drill"]);
    std::fs::write(target.join("lib.rs"), "replacement source\n").expect("write replacement tree");
    git(&target, &["add", "lib.rs"]);
    let replacement_tree = git(&target, &["write-tree"]);
    let replacement = git(
        &target,
        &[
            "commit-tree",
            &replacement_tree,
            "-p",
            pinned,
            "-m",
            "replacement tree",
        ],
    );
    git(&target, &["replace", pinned, &replacement]);
    assert_eq!(git(&target, &["rev-parse", "HEAD"]), pinned);
    assert_eq!(
        git(&target, &["status", "--porcelain"]),
        "",
        "ordinary replacement-aware Git must demonstrate the false-clean gap"
    );

    let prior_provenance = std::fs::read(dest.join("constellation-bootstrap.json"))
        .expect("seed provenance remains available");
    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone accepted replacement-defined locked source:\n{standalone_text}"
    );
    assert!(
        standalone_text.contains("replace") || standalone_text.contains("tracked"),
        "replacement refusal must be actionable:\n{standalone_text}"
    );
    for mode in [Some("--verify-only"), Some("--snapshot")] {
        let (shell_ok, shell_text) = run_shell_checkout(&c, mode);
        assert!(
            !shell_ok,
            "shell mode {mode:?} accepted replacement-defined locked source:\n{shell_text}"
        );
        assert!(
            shell_text.contains("replace") || shell_text.contains("tracked"),
            "replacement refusal must be actionable:\n{shell_text}"
        );
    }
    assert_eq!(
        std::fs::read(dest.join("constellation-bootstrap.json"))
            .expect("refusal retains prior provenance"),
        prior_provenance,
    );
}

#[test]
fn graft_authority_cannot_redefine_locked_source() {
    let c = make_constellation("graft-authority-refusal");
    let mirror = c.mirror.to_str().expect("UTF-8 mirror").to_string();
    let (seeded, seed_text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(seeded, "seed exact pinned siblings:\n{seed_text}");
    install_shell_checkout(&c);

    let dest = c.root.parent().expect("constellation destination");
    let target = dest.join("asupersync");
    let pinned = c
        .heads
        .iter()
        .find(|(name, _)| name == "asupersync")
        .map(|(_, head)| head.as_str())
        .expect("asupersync pin");
    write_git_admin_file(&target, "info/grafts", format!("{pinned}\n"));
    let prior_provenance = std::fs::read(dest.join("constellation-bootstrap.json"))
        .expect("seed provenance remains available");

    let (standalone_ok, standalone_text) = run_bootstrap(&c, &["--offline"]);
    assert!(
        !standalone_ok,
        "standalone accepted a legacy graft authority:\n{standalone_text}"
    );
    assert!(
        standalone_text.contains("graft"),
        "standalone graft refusal must be actionable:\n{standalone_text}"
    );
    for mode in [Some("--verify-only"), Some("--snapshot")] {
        let (shell_ok, shell_text) = run_shell_checkout(&c, mode);
        assert!(
            !shell_ok,
            "shell mode {mode:?} accepted a legacy graft authority:\n{shell_text}"
        );
        assert!(
            shell_text.contains("graft"),
            "shell graft refusal must be actionable:\n{shell_text}"
        );
    }
    assert_eq!(
        std::fs::read(dest.join("constellation-bootstrap.json"))
            .expect("graft refusal retains prior provenance"),
        prior_provenance,
    );
}

#[test]
fn drift_dirty_and_offline_missing_all_refuse() {
    let c = make_constellation("refuse");
    let mirror = c.mirror.to_str().expect("utf8").to_string();
    let (ok, _) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(ok, "seed the cache");
    let dest = c.root.parent().unwrap().to_path_buf();

    // DIRTY: modify a file at the pinned head → refuse, name the tree.
    let alpha = dest.join("asupersync");
    let alpha_source = std::fs::read(alpha.join("lib.rs")).expect("read pinned source");
    std::fs::write(alpha.join("lib.rs"), "tampered\n").expect("tamper");
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "dirty sibling refuses");
    assert!(
        text.contains("DIRTY") && text.contains("asupersync"),
        "{text}"
    );
    std::fs::write(alpha.join("lib.rs"), alpha_source).expect("restore pinned source bytes");

    // DRIFT: advance the sibling one commit past the pin → refuse.
    let beta = dest.join("frankensqlite");
    std::fs::write(beta.join("extra.rs"), "pub fn c() {}\n").expect("write");
    git(&beta, &["add", "extra.rs"]);
    git(&beta, &["config", "user.email", "drill@frankensim.test"]);
    git(&beta, &["config", "user.name", "bootstrap drill"]);
    git(&beta, &["commit", "-qm", "drifted"]);
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "drifted sibling refuses");
    assert!(
        text.contains("refusing to repurpose") && text.contains("frankensqlite"),
        "{text}"
    );

    // OFFLINE + MISSING: retain a sibling at a different path → structured
    // refusal, and the
    // still-present drifted sibling's refusal is independent (fail
    // closed per library).
    std::fs::rename(&alpha, dest.join("asupersync-retained-for-missing"))
        .expect("retain sibling away from its admitted path");
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "missing sibling refuses offline");
    assert!(
        text.contains("missing from the source cache in --offline mode"),
        "{text}"
    );
}
