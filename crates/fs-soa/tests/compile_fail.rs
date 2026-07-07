//! Compile-fail battery for `#[derive(Soa)]` diagnostics (wf9.5):
//! every unsupported shape must produce OUR structured compile error —
//! no silent fallbacks. In-house harness (no trybuild — Franken-only
//! law): a scratch cargo project with path-deps back to fs-soa is
//! `cargo check`ed offline per fixture and stderr is asserted to
//! contain the expected `compile_error!` text.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

const CASES: &[(&str, &str, &str)] = &[
    (
        "tuple_struct",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct T(f64, u32);\n",
        "tuple structs have no field names",
    ),
    (
        "unit_struct",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct U;\n",
        "unit structs have no fields",
    ),
    (
        "enum_input",
        "use fs_soa::Soa;\n#[derive(Soa)]\nenum E { A, B }\n",
        "not enums",
    ),
    (
        "lifetime_param",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct L<'a> { r: &'a f64 }\n",
        "does not support lifetime parameters",
    ),
    (
        "generic_default",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct D<T: Copy = f64> { x: T }\n",
        "does not support generic parameter defaults",
    ),
    (
        "zero_fields",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct Z {}\n",
        "requires at least one field",
    ),
    (
        "reserved_name",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct R { len: f64 }\n",
        "collides with the generated container API",
    ),
    (
        "unknown_attr",
        "use fs_soa::Soa;\n#[derive(Soa)]\nstruct A { #[soa(bogus)] x: f64 }\n",
        "unknown #[soa(",
    ),
];

#[test]
fn derive_diagnostics_are_ours_and_clear() {
    let soa_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let root = std::env::temp_dir().join(format!("fs_soa_compile_fail_{}", std::process::id()));
    let src = root.join("src");
    fs::create_dir_all(&src).expect("scratch dirs");
    fs::write(
        root.join("Cargo.toml"),
        format!(
            "[package]\nname = \"soa-fixture\"\nversion = \"0.0.0\"\nedition = \"2024\"\n\n\
             [dependencies]\nfs-soa = {{ path = {soa_dir:?} }}\n\n[workspace]\n"
        ),
    )
    .expect("scratch manifest");

    for (case, source, expected) in CASES {
        fs::write(src.join("lib.rs"), source).expect("fixture source");
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
            "fixture `{case}` compiled but must fail:\n{source}"
        );
        assert!(
            stderr.contains(expected),
            "fixture `{case}`: expected diagnostic containing {expected:?}, got:\n{stderr}"
        );
        println!(
            "{{\"suite\":\"fs-soa\",\"case\":\"compile-fail-{case}\",\"verdict\":\"pass\",\
             \"detail\":\"rejected with expected diagnostic\"}}"
        );
    }
}
