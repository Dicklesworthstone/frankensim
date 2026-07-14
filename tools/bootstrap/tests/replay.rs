//! Offline-cache replay drills (bead 1t8i): the bootstrap binary run
//! against a SYNTHETIC constellation — a temp checkout root with its own
//! lock, and a local mirror base standing in for the network — so every
//! trust rule is exercised hermetically:
//!
//! - clean machine: missing siblings clone from the mirror, detached at
//!   the pinned head, provenance written;
//! - idempotence/replay: a second run verifies (no re-clone, identical
//!   provenance);
//! - interrupted replay: clean marked and exact-origin unborn destinations
//!   resume in place, verify the exact pin, and clear the marker;
//! - adoption refusal: ordinary non-git, dirty marked, and wrong-origin unborn
//!   destinations are never repurposed;
//! - drift refusal: a sibling at the wrong head refuses;
//! - dirty refusal: a modified tree at the right head refuses;
//! - hidden-dirt refusal: local ignore rules and index flags cannot conceal
//!   untracked or modified source from the verifier;
//! - post-clone dirty refusal: checkout-time mutation is caught before a
//!   freshly cloned sibling can be accepted;
//! - `--offline` refusal: a missing sibling is a structured failure and
//!   the network (here: the mirror) is never touched;
//! - CLI admission: value-taking flags refuse absent operands before any
//!   lock or sibling access.

use std::path::{Path, PathBuf};
use std::process::Command;

const LOCK_NOTE: &str = "lock_hash covers (lib, version, git_head) only — paths are per-machine; remote is transport for bootstrap-constellation (content identity is the git head)";
const LOCK_IDENTITY_DOMAIN: &str = "org.frankensim.xtask.constellation-lock.v1";
const LOCK_IDENTITY_VERSION: u32 = 1;
const BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN: &str =
    "org.frankensim.xtask.constellation-bootstrap-provenance.v2";
const BOOTSTRAP_PROVENANCE_IDENTITY_VERSION: u32 = 2;
const TRACKED_LOCK: &str = include_str!("../../../constellation.lock");
const CHECKOUT_CONSTELLATION_SH: &str =
    include_str!("../../../scripts/ci/checkout_constellation.sh");

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut value = 0xcbf2_9ce4_8422_2325u64;
    for &byte in bytes {
        value ^= u64::from(byte);
        value = value.wrapping_mul(0x0000_0100_0000_01b3);
    }
    value
}

fn git(dir: &Path, args: &[&str]) -> String {
    let out = Command::new("git")
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
    Command::new("git")
        .arg("-C")
        .arg(dir)
        .args(args)
        .output()
        .expect("git spawns")
        .status
        .success()
}

/// A synthetic upstream repo with one committed file; returns its head.
fn make_upstream(base: &Path, name: &str, content: &str) -> String {
    let dir = base.join(name);
    std::fs::create_dir_all(&dir).expect("mkdir");
    git(&dir, &["init", "-q"]);
    git(&dir, &["config", "user.email", "drill@frankensim.test"]);
    git(&dir, &["config", "user.name", "bootstrap drill"]);
    std::fs::write(dir.join("lib.rs"), content).expect("write");
    git(&dir, &["add", "lib.rs"]);
    git(&dir, &["commit", "-qm", "pinned"]);
    git(&dir, &["rev-parse", "HEAD"])
}

struct Constellation {
    /// Keep the temp tree alive for the test's duration.
    base: PathBuf,
    root: PathBuf,
    mirror: PathBuf,
    heads: Vec<(String, String)>,
    lock_hash: String,
}

impl Drop for Constellation {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.base);
    }
}

/// Build a synthetic world: upstreams, a bare-clone mirror base, and a
/// fake frankensim checkout (parent = the bootstrap destination) whose
/// lock pins the upstream heads. Library names deliberately NOT in the
/// tool's rename map, so dirname == lib.
fn make_constellation(tag: &str) -> Constellation {
    let base =
        std::env::temp_dir().join(format!("fs-bootstrap-drill-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    let upstreams = base.join("upstreams");
    let mirror = base.join("mirror");
    let parent = base.join("dest");
    let root = parent.join("frankensim");
    std::fs::create_dir_all(&root).expect("mkdir root");
    std::fs::create_dir_all(&mirror).expect("mkdir mirror");
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
        let head = make_upstream(
            &upstreams,
            name,
            &format!("pub fn {}_fixture() {{}}\n", name.replace('-', "_")),
        );
        // Bare mirror clone = the offline cache / air-gapped transport.
        let out = Command::new("git")
            .args(["clone", "-q", "--bare"])
            .arg(upstreams.join(name))
            .arg(mirror.join(name))
            .output()
            .expect("git spawns");
        assert!(out.status.success(), "bare mirror clone");
        heads.push((name.to_string(), head));
    }
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
    }
}

fn run_bootstrap(c: &Constellation, extra: &[&str]) -> (bool, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim-bootstrap"));
    command.arg("--root").arg(&c.root).args(extra);
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

fn run_shell_checkout(c: &Constellation, mode: Option<&str>) -> (bool, String) {
    let mut command = Command::new("bash");
    command.arg(c.root.join("scripts/ci/checkout_constellation.sh"));
    if let Some(mode) = mode {
        command.arg(mode);
    }
    command.arg(c.root.parent().expect("constellation destination"));
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
        assert_eq!(git(&checkout, &["status", "--porcelain"]), "");
    }

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
        assert_eq!(git(&dir, &["status", "--porcelain"]), "");
        let detached = Command::new("git")
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

    // REPLAY: second run is pure verification — same provenance bytes
    // except the recorded state flips cloned → verified.
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
        "identical invocation and source state produce byte-identical v2 provenance"
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
        assert_eq!(git(&checkout, &["status", "--porcelain"]), "");
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
fn repository_local_excludes_cannot_hide_untracked_source() {
    let c = make_constellation("hidden-untracked");
    let mirror = c.mirror.to_str().expect("utf8").to_string();
    let (ok, text) = run_bootstrap(&c, &["--from", &mirror]);
    assert!(ok, "seed the source cache:\n{text}");

    let sibling = c.root.parent().unwrap().join("asupersync");
    std::fs::write(sibling.join(".git/info/exclude"), "hidden-source.rs\n")
        .expect("install repository-local exclusion");
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
fn fresh_clone_is_rechecked_for_checkout_time_dirt() {
    use std::os::unix::fs::PermissionsExt;

    let c = make_constellation("post-clone-dirty");
    let hooks = c.base.join("hooks");
    std::fs::create_dir_all(&hooks).expect("mkdir hooks");
    let post_checkout = hooks.join("post-checkout");
    std::fs::write(
        &post_checkout,
        "#!/bin/sh\nprintf 'checkout-time mutation\\n' > lib.rs\n",
    )
    .expect("write post-checkout hook");
    let mut permissions = std::fs::metadata(&post_checkout)
        .expect("hook metadata")
        .permissions();
    permissions.set_mode(0o755);
    std::fs::set_permissions(&post_checkout, permissions).expect("make hook executable");

    let mirror = c.mirror.to_str().expect("utf8");
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim-bootstrap"));
    command
        .arg("--root")
        .arg(&c.root)
        .args(["--from", mirror])
        .env("GIT_CONFIG_COUNT", "1")
        .env("GIT_CONFIG_KEY_0", "core.hooksPath")
        .env("GIT_CONFIG_VALUE_0", &hooks);
    let (ok, text) = run_command(&mut command);

    assert!(!ok, "a dirty post-checkout tree must refuse:\n{text}");
    assert!(
        text.contains("DIRTY") && text.contains("asupersync"),
        "fresh-clone cleanliness refusal is explicit:\n{text}"
    );
    assert!(
        !c.root
            .parent()
            .unwrap()
            .join("constellation-bootstrap.json")
            .exists(),
        "failed bootstrap must not publish success provenance"
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
    std::fs::write(alpha.join("lib.rs"), "tampered\n").expect("tamper");
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "dirty sibling refuses");
    assert!(
        text.contains("DIRTY") && text.contains("asupersync"),
        "{text}"
    );
    git(&alpha, &["checkout", "--", "lib.rs"]);

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

    // OFFLINE + MISSING: remove a sibling → structured refusal, and the
    // still-present drifted sibling's refusal is independent (fail
    // closed per library).
    std::fs::remove_dir_all(&alpha).expect("remove sibling");
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "missing sibling refuses offline");
    assert!(
        text.contains("missing from the source cache in --offline mode"),
        "{text}"
    );
}
