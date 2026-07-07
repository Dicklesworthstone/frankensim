//! Compile-fail battery for `kernel!` static analysis (wf9.11): every
//! forbidden construct must produce OUR structured compile error with
//! an actionable message — no silent fallbacks. In-house offline
//! harness (the fs-soa pattern): scratch cargo project, path deps,
//! `cargo check` per fixture, stderr asserted.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

const PRELUDE: &str = "use fs_tilelang::kernel;\n";

const CASES: &[(&str, &str, &str)] = &[
    (
        "alias",
        "kernel! { name: k, reads: [a], writes: [a], reduction: none, body: { a = a; } }",
        "declared twice",
    ),
    (
        "unsafe_body",
        "kernel! { name: k, reads: [a], writes: [o], reduction: none, body: { o = unsafe { a }; } }",
        "unsafe blocks are not allowed",
    ),
    (
        "user_loop",
        "kernel! { name: k, reads: [a], writes: [o], reduction: none, body: { for _ in 0..2 {} o = a; } }",
        "loops are not allowed",
    ),
    (
        "allocation",
        "kernel! { name: k, reads: [a], writes: [o], reduction: none, body: { let v = vec![a]; o = v[0]; } }",
        "allocation is not allowed",
    ),
    (
        "unknown_reduction",
        "kernel! { name: k, reads: [a], writes: [o], reduction: median, body: { o = a; } }",
        "unknown reduction",
    ),
    (
        "acc_without_reduction",
        "kernel! { name: k, reads: [a], writes: [o], reduction: none, body: { acc = a; o = a; } }",
        "reduction is `none`",
    ),
    (
        "reduction_without_acc",
        "kernel! { name: k, reads: [a], writes: [], reduction: deterministic_sum, body: { let t = a; } }",
        "never assigns `acc`",
    ),
    (
        "unassigned_write",
        "kernel! { name: k, reads: [a], writes: [o], reduction: none, body: { let t = a; } }",
        "never assigned",
    ),
    (
        "gather_undeclared",
        "kernel! { name: k, reads: [a], writes: [o], reduction: none, body: { o = gather(b, 0); } }",
        "not a declared read buffer",
    ),
    (
        "reserved_name",
        "kernel! { name: k, reads: [acc], writes: [o], reduction: none, body: { o = acc; } }",
        "reserved by the kernel grammar",
    ),
];

#[test]
fn kernel_rejections_are_ours_and_actionable() {
    let tilelang_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let root =
        std::env::temp_dir().join(format!("fs_tilelang_compile_fail_{}", std::process::id()));
    let src = root.join("src");
    fs::create_dir_all(&src).expect("scratch dirs");
    fs::write(
        root.join("Cargo.toml"),
        format!(
            "[package]\nname = \"tile-fixture\"\nversion = \"0.0.0\"\nedition = \"2024\"\n\n\
             [dependencies]\nfs-tilelang = {{ path = {tilelang_dir:?} }}\n\n[workspace]\n"
        ),
    )
    .expect("scratch manifest");
    for (case, body, expected) in CASES {
        fs::write(src.join("lib.rs"), format!("{PRELUDE}{body}\n")).expect("fixture source");
        let out = Command::new("cargo")
            .args(["check", "--offline", "--quiet"])
            .current_dir(&root)
            .env("RCH_DISABLE", "1")
            .env("CARGO_TARGET_DIR", root.join("target"))
            .output()
            .expect("cargo check runs");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            !out.status.success(),
            "fixture `{case}` compiled but must fail:\n{body}"
        );
        assert!(
            stderr.contains(expected),
            "fixture `{case}`: expected diagnostic containing {expected:?}, got:\n{stderr}"
        );
        println!(
            "{{\"suite\":\"fs-tilelang\",\"case\":\"compile-fail-{case}\",\"verdict\":\"pass\",\
             \"detail\":\"rejected with actionable diagnostic\"}}"
        );
    }
}
