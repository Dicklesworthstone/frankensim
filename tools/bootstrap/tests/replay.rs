//! Offline-cache replay drills (bead 1t8i): the bootstrap binary run
//! against a SYNTHETIC constellation — a temp checkout root with its own
//! lock, and a local mirror base standing in for the network — so every
//! trust rule is exercised hermetically:
//!
//! - clean machine: missing siblings clone from the mirror, detached at
//!   the pinned head, provenance written;
//! - idempotence/replay: a second run verifies (no re-clone, identical
//!   provenance);
//! - drift refusal: a sibling at the wrong head refuses;
//! - dirty refusal: a modified tree at the right head refuses;
//! - `--offline` refusal: a missing sibling is a structured failure and
//!   the network (here: the mirror) is never touched.

use std::path::{Path, PathBuf};
use std::process::Command;

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
    for (name, content) in [
        ("drill_alpha", "pub fn a() {}\n"),
        ("drill_beta", "pub fn b() {}\n"),
    ] {
        let head = make_upstream(&upstreams, name, content);
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
    let rows: Vec<String> = heads
        .iter()
        .map(|(name, head)| {
            format!(
                "    {{\"lib\": \"{name}\", \"version\": \"0.0.0\", \"git_head\": \"{head}\", \"remote\": \"no-remote\", \"path\": \"unused\"}}"
            )
        })
        .collect();
    let lock = format!(
        "{{\n\"schema\": 2,\n\"lock_hash\": \"drill-{tag}\",\n\"libraries\": [\n{}\n]\n}}\n",
        rows.join(",\n")
    );
    std::fs::write(root.join("constellation.lock"), lock).expect("write lock");
    Constellation {
        base,
        root,
        mirror,
        heads,
    }
}

fn run_bootstrap(c: &Constellation, extra: &[&str]) -> (bool, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_frankensim-bootstrap"))
        .arg("--root")
        .arg(&c.root)
        .args(extra)
        .output()
        .expect("bootstrap binary spawns");
    let text = format!(
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    (out.status.success(), text)
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
    assert!(prov1.contains("frankensim-constellation-bootstrap-v1"));
    assert!(prov1.contains("drill-replay"), "lock hash bound: {prov1}");

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
    assert_eq!(
        prov1.replace("\"state\": \"cloned\"", "\"state\": \"verified\""),
        prov2,
        "replay provenance is byte-identical modulo the cloned→verified state"
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
    let alpha = dest.join("drill_alpha");
    std::fs::write(alpha.join("lib.rs"), "tampered\n").expect("tamper");
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "dirty sibling refuses");
    assert!(
        text.contains("DIRTY") && text.contains("drill_alpha"),
        "{text}"
    );
    git(&alpha, &["checkout", "--", "lib.rs"]);

    // DRIFT: advance the sibling one commit past the pin → refuse.
    let beta = dest.join("drill_beta");
    std::fs::write(beta.join("extra.rs"), "pub fn c() {}\n").expect("write");
    git(&beta, &["add", "extra.rs"]);
    git(&beta, &["config", "user.email", "drill@frankensim.test"]);
    git(&beta, &["config", "user.name", "bootstrap drill"]);
    git(&beta, &["commit", "-qm", "drifted"]);
    let (ok, text) = run_bootstrap(&c, &["--offline"]);
    assert!(!ok, "drifted sibling refuses");
    assert!(
        text.contains("refusing to silently substitute") && text.contains("drill_beta"),
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
