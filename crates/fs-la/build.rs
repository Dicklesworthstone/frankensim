//! Generates the fail-closed GEMM codegen identity used by durable tune keys.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

mod build_identity_support;
mod depgraph_receipt_format;

use build_identity_support::{
    ASUPERSYNC_NON_SRC_INPUTS, GEMM_BUILD_PAYLOAD_SCHEMA, GemmBuildIdentityInput,
    append_executable_identity, append_external_identity, append_source_fields,
    gemm_build_fingerprint, push_field, push_optional_field,
};

const CARGO_PROFILES: [&str; 4] = ["DEV", "RELEASE", "TEST", "BENCH"];
const PROFILE_CODEGEN_KEYS: [&str; 11] = [
    "OPT_LEVEL",
    "LTO",
    "CODEGEN_UNITS",
    "DEBUG",
    "DEBUG_ASSERTIONS",
    "OVERFLOW_CHECKS",
    "PANIC",
    "INCREMENTAL",
    "RPATH",
    "STRIP",
    "SPLIT_DEBUGINFO",
];
fn required_env(name: &str) -> String {
    println!("cargo:rerun-if-env-changed={name}");
    env::var(name)
        .unwrap_or_else(|error| panic!("{name} is required for GEMM build identity: {error}"))
}

fn optional_env(name: &str) -> Option<String> {
    println!("cargo:rerun-if-env-changed={name}");
    match env::var(name) {
        Ok(value) => Some(value),
        Err(env::VarError::NotPresent) => None,
        Err(env::VarError::NotUnicode(_)) => {
            panic!("{name} is non-Unicode and cannot be represented in GEMM build identity")
        }
    }
}

fn read_identity_file(payload: &mut Vec<u8>, name: &str, path: &Path) {
    println!("cargo:rerun-if-changed={}", path.display());
    let bytes = read_required_file(path);
    push_field(payload, name, &bytes);
}

fn read_required_file(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|error| {
        panic!(
            "cannot read required GEMM build-identity input {}: {error}",
            path.display()
        )
    })
}

fn resolve_executable(name: &str, value: &str) -> PathBuf {
    let declared = PathBuf::from(value);
    let candidate = if declared.is_absolute() || declared.components().count() > 1 {
        declared
    } else {
        env::var_os("PATH")
            .into_iter()
            .flat_map(|path| env::split_paths(&path).collect::<Vec<_>>())
            .map(|directory| directory.join(&declared))
            .find(|path| path.is_file())
            .unwrap_or_else(|| panic!("cannot resolve {name} executable {value:?} on PATH"))
    };
    println!("cargo:rerun-if-changed={}", candidate.display());
    let resolved = candidate.canonicalize().unwrap_or_else(|error| {
        panic!(
            "cannot canonicalize {name} executable {}: {error}",
            candidate.display()
        )
    });
    assert!(
        resolved.is_file(),
        "resolved {name} executable is not a regular file: {}",
        resolved.display()
    );
    println!("cargo:rerun-if-changed={}", resolved.display());
    resolved
}

fn add_executable_identity(payload: &mut Vec<u8>, name: &str, value: &str) {
    let path = resolve_executable(name, value);
    let canonical_path = path.to_str().unwrap_or_else(|| {
        panic!(
            "resolved {name} executable path is not Unicode: {}",
            path.display()
        )
    });
    let bytes = read_required_file(&path);
    append_executable_identity(payload, name, canonical_path, &bytes);
}

fn collect_rust_sources(dir: &Path, files: &mut Vec<PathBuf>) {
    let entries = std::fs::read_dir(dir).unwrap_or_else(|error| {
        panic!(
            "cannot enumerate GEMM build-identity source tree {}: {error}",
            dir.display()
        )
    });
    for entry in entries {
        let entry = entry.unwrap_or_else(|error| {
            panic!(
                "cannot read GEMM build-identity directory entry in {}: {error}",
                dir.display()
            )
        });
        let path = entry.path();
        let file_type = entry.file_type().unwrap_or_else(|error| {
            panic!(
                "cannot inspect GEMM build-identity source {}: {error}",
                path.display()
            )
        });
        if file_type.is_dir() {
            collect_rust_sources(&path, files);
        } else if file_type.is_file() && path.extension().is_some_and(|extension| extension == "rs")
        {
            files.push(path);
        }
    }
}

fn collect_regular_files(dir: &Path, files: &mut Vec<PathBuf>) {
    let entries = std::fs::read_dir(dir).unwrap_or_else(|error| {
        panic!(
            "cannot enumerate GEMM build-identity source tree {}: {error}",
            dir.display()
        )
    });
    for entry in entries {
        let entry = entry.unwrap_or_else(|error| {
            panic!(
                "cannot read GEMM build-identity directory entry in {}: {error}",
                dir.display()
            )
        });
        let path = entry.path();
        let file_type = entry.file_type().unwrap_or_else(|error| {
            panic!(
                "cannot inspect GEMM build-identity source {}: {error}",
                path.display()
            )
        });
        if file_type.is_dir() {
            collect_regular_files(&path, files);
        } else if file_type.is_file() {
            files.push(path);
        }
    }
}

fn add_source_closure(payload: &mut Vec<u8>, workspace_root: &Path) {
    let crate_roots = [
        "crates/fs-la",
        "crates/fs-simd",
        "crates/fs-exec",
        "crates/fs-alloc",
        "crates/fs-substrate",
        "crates/fs-blake3",
        "crates/fs-obs",
    ];
    let mut files = Vec::new();
    for relative_root in crate_roots {
        let crate_root = workspace_root.join(relative_root);
        println!(
            "cargo:rerun-if-changed={}",
            crate_root.join("src").display()
        );
        files.push(crate_root.join("Cargo.toml"));
        let build_script = crate_root.join("build.rs");
        println!("cargo:rerun-if-changed={}", build_script.display());
        if build_script.is_file() {
            files.push(build_script);
        }
        collect_rust_sources(&crate_root.join("src"), &mut files);
    }
    files.push(workspace_root.join("crates/fs-la/build_identity_support.rs"));
    files.push(workspace_root.join("crates/fs-la/depgraph_receipt_format.rs"));

    let mut source_fields = Vec::with_capacity(files.len());
    for path in files {
        let relative = path.strip_prefix(workspace_root).unwrap_or_else(|error| {
            panic!(
                "GEMM build-identity source {} escaped workspace root {}: {error}",
                path.display(),
                workspace_root.display()
            )
        });
        let relative = relative.to_str().unwrap_or_else(|| {
            panic!(
                "GEMM build-identity source path is not Unicode: {}",
                relative.display()
            )
        });
        let canonical_relative = relative.replace(std::path::MAIN_SEPARATOR, "/");
        println!("cargo:rerun-if-changed={}", path.display());
        source_fields.push((canonical_relative, read_required_file(&path)));
    }
    append_source_fields(payload, source_fields);
}

fn command_stdout(command: &mut Command, context: &str) -> Vec<u8> {
    let output = command
        .output()
        .unwrap_or_else(|error| panic!("cannot execute {context}: {error}"));
    assert!(
        output.status.success(),
        "{context} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn add_asupersync_identity(
    payload: &mut Vec<u8>,
    workspace_root: &Path,
    constellation_lock: &[u8],
) {
    let checkout = workspace_root.join("../asupersync");
    let git_dir_bytes = command_stdout(
        Command::new("git")
            .arg("-C")
            .arg(&checkout)
            .args(["rev-parse", "--absolute-git-dir"]),
        "git-dir discovery for the required asupersync checkout",
    );
    let git_dir_text = String::from_utf8_lossy(&git_dir_bytes);
    let git_dir = PathBuf::from(git_dir_text.trim());
    let head_path = git_dir.join("HEAD");
    let head_file = read_required_file(&head_path);
    println!("cargo:rerun-if-changed={}", head_path.display());
    println!("cargo:rerun-if-changed={}", git_dir.join("index").display());
    println!(
        "cargo:rerun-if-changed={}",
        git_dir.join("packed-refs").display()
    );
    if let Some(symbolic_ref) = head_file.strip_prefix(b"ref: ") {
        let symbolic_ref = String::from_utf8_lossy(symbolic_ref);
        println!(
            "cargo:rerun-if-changed={}",
            git_dir.join(symbolic_ref.trim()).display()
        );
    }
    let head = command_stdout(
        Command::new("git")
            .arg("-C")
            .arg(&checkout)
            .args(["rev-parse", "HEAD"]),
        "git rev-parse for the required asupersync checkout",
    );
    let head_text = String::from_utf8_lossy(&head);
    let head_text = head_text.trim();
    assert!(
        head_text.len() == 40 && head_text.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "asupersync HEAD is not a full Git object id: {head_text:?}"
    );
    let package_roots = [
        "",
        "asupersync-macros",
        "franken_decision",
        "franken_evidence",
        "franken_kernel",
    ];
    let mut files = Vec::new();
    for relative_root in package_roots {
        let package_root = checkout.join(relative_root);
        let source_root = package_root.join("src");
        println!("cargo:rerun-if-changed={}", source_root.display());
        collect_regular_files(&source_root, &mut files);
        files.push(package_root.join("Cargo.toml"));
        let build_script = package_root.join("build.rs");
        println!("cargo:rerun-if-changed={}", build_script.display());
        if build_script.is_file() {
            files.push(build_script);
        }
    }
    for relative in ASUPERSYNC_NON_SRC_INPUTS {
        files.push(checkout.join(relative));
    }

    let mut source_fields = Vec::with_capacity(files.len());
    for path in files {
        let relative = path.strip_prefix(&checkout).unwrap_or_else(|error| {
            panic!(
                "asupersync identity source {} escaped checkout {}: {error}",
                path.display(),
                checkout.display()
            )
        });
        let relative = relative.to_str().unwrap_or_else(|| {
            panic!(
                "asupersync identity source path is not Unicode: {}",
                relative.display()
            )
        });
        let canonical_relative = relative.replace(std::path::MAIN_SEPARATOR, "/");
        println!("cargo:rerun-if-changed={}", path.display());
        source_fields.push((
            format!("external/asupersync/{canonical_relative}"),
            read_required_file(&path),
        ));
    }
    append_external_identity(payload, constellation_lock, head_text, source_fields);
}

struct DepgraphEvidence {
    class_identity: String,
    kind: &'static str,
    receipt: Option<String>,
    receipt_digest: Option<String>,
}

/// Resolved dependency-graph evidence for feature unification (bead fz2.6).
///
/// `CARGO_FEATURE_*` only covers this crate's own features: Cargo can compile
/// path and registry dependencies (asupersync and its aes-gcm/rand_core
/// closure) under a different unified feature set without changing any other
/// fingerprint input, letting a tune row cross binaries whose dependency
/// codegen differs. Build tooling derives a canonical receipt of the resolved
/// normal/build closure OUTSIDE this build script
/// (`cargo run -p xtask -- depgraph-receipt -- --package <ROOT>`) and exports it as
/// `FRANKENSIM_DEPGRAPH_RECEIPT`; interactive workspace builds carry the
/// explicit `FRANKENSIM_DEPGRAPH_SALT` from `.cargo/config.toml`, which marks
/// all such builds as ONE deliberate development equivalence class, never
/// verified graph evidence. Neither present → fail closed. Receipt wins when both are set
/// (tooling overrides the workspace default).
fn add_depgraph_evidence(payload: &mut Vec<u8>) -> DepgraphEvidence {
    let receipt = optional_env("FRANKENSIM_DEPGRAPH_RECEIPT");
    let salt = optional_env("FRANKENSIM_DEPGRAPH_SALT");
    if let Some(receipt) = receipt {
        let parsed = depgraph_receipt_format::parse(&receipt).unwrap_or_else(|error| {
            panic!(
                "FRANKENSIM_DEPGRAPH_RECEIPT is not the strict canonical single-root fs-la \
                 receipt minted by xtask: {error}"
            )
        });
        push_field(payload, "depgraph-receipt", receipt.as_bytes());
        let digest = depgraph_receipt_format::content_identity(&parsed)
            .expect("a parsed canonical dependency receipt has a content identity")
            .to_hex();
        push_field(payload, "depgraph-receipt-domain-digest", digest.as_bytes());
        return DepgraphEvidence {
            class_identity: format!("receipt:{digest}"),
            kind: "operator-observed-receipt",
            receipt: Some(receipt),
            receipt_digest: Some(digest),
        };
    }
    if let Some(salt) = salt.as_deref() {
        assert!(
            !salt.is_empty()
                && salt.len() <= 128
                && salt
                    .bytes()
                    .all(|byte| byte.is_ascii_graphic() || byte == b' '),
            "FRANKENSIM_DEPGRAPH_SALT must be short printable ASCII, got {salt:?}"
        );
        push_field(payload, "depgraph-salt", salt.as_bytes());
        return DepgraphEvidence {
            class_identity: format!("salt:{salt}"),
            kind: "development-equivalence-salt",
            receipt: None,
            receipt_digest: None,
        };
    }
    panic!(
        "GEMM build identity requires dependency-graph evidence (bead fz2.6): \
         export FRANKENSIM_DEPGRAPH_RECEIPT via \
         `cargo run -p xtask -- depgraph-receipt -- --package <ROOT>` for \
         operator-observed receipt builds, \
         or build inside the workspace whose .cargo/config.toml supplies the \
         explicit FRANKENSIM_DEPGRAPH_SALT equivalence class"
    );
}

#[allow(clippy::too_many_lines)] // one ordered payload defines the complete code-generation identity
fn main() {
    let mut payload = Vec::new();
    push_field(&mut payload, "schema", GEMM_BUILD_PAYLOAD_SCHEMA.as_bytes());
    let graph_evidence = add_depgraph_evidence(&mut payload);
    println!(
        "cargo:rustc-env=FS_LA_GEMM_GRAPH_EVIDENCE={}",
        graph_evidence.class_identity
    );
    println!(
        "cargo:rustc-env=FS_LA_GEMM_GRAPH_EVIDENCE_KIND={}",
        graph_evidence.kind
    );
    let out_dir = PathBuf::from(required_env("OUT_DIR"));
    let receipt_path = out_dir.join("fs_la_depgraph_receipt.json");
    std::fs::write(
        &receipt_path,
        graph_evidence.receipt.as_deref().unwrap_or("").as_bytes(),
    )
    .unwrap_or_else(|error| {
        panic!(
            "cannot write canonical dependency receipt {}: {error}",
            receipt_path.display()
        )
    });
    if graph_evidence.receipt.is_some() {
        println!("cargo:rustc-env=FS_LA_GEMM_HAS_DEPGRAPH_RECEIPT=1");
    }
    if let Some(digest) = &graph_evidence.receipt_digest {
        println!("cargo:rustc-env=FS_LA_GEMM_DEPGRAPH_RECEIPT_DIGEST={digest}");
    }

    for name in [
        "PROFILE",
        "OPT_LEVEL",
        "DEBUG",
        "TARGET",
        "HOST",
        "CARGO_CFG_TARGET_ARCH",
        "CARGO_CFG_TARGET_ENDIAN",
        "CARGO_CFG_TARGET_ENV",
        "CARGO_CFG_TARGET_FEATURE",
        "CARGO_CFG_TARGET_OS",
        "CARGO_CFG_TARGET_POINTER_WIDTH",
        "CARGO_CFG_TARGET_VENDOR",
    ] {
        let value = required_env(name);
        push_field(&mut payload, name, value.as_bytes());
    }

    for name in [
        "CARGO_ENCODED_RUSTFLAGS",
        "CARGO_CFG_PANIC",
        "CARGO_CFG_TARGET_ABI",
        "CARGO_INCREMENTAL",
        "FRANKENSIM_GEMM_CODEGEN_ID",
    ] {
        let value = optional_env(name);
        push_optional_field(&mut payload, name, value.as_deref().map(str::as_bytes));
    }
    let rustc_wrapper = optional_env("RUSTC_WRAPPER");
    push_optional_field(
        &mut payload,
        "RUSTC_WRAPPER",
        rustc_wrapper.as_deref().map(str::as_bytes),
    );
    let rustc_workspace_wrapper = optional_env("RUSTC_WORKSPACE_WRAPPER");
    push_optional_field(
        &mut payload,
        "RUSTC_WORKSPACE_WRAPPER",
        rustc_workspace_wrapper.as_deref().map(str::as_bytes),
    );
    let path_identity = optional_env("PATH");
    push_optional_field(
        &mut payload,
        "PATH",
        path_identity.as_deref().map(str::as_bytes),
    );

    // Cargo does not export every resolved profile knob to build scripts.
    // Watch each supported override even while it is absent, then include its
    // explicit `<unset>` state in the fingerprint. This closes the same-target
    // invalidation hole when an operator later supplies one of these values.
    let mut watched_profile_vars = Vec::new();
    for profile in CARGO_PROFILES {
        for key in PROFILE_CODEGEN_KEYS {
            let name = format!("CARGO_PROFILE_{profile}_{key}");
            let value = optional_env(&name);
            push_optional_field(&mut payload, &name, value.as_deref().map(str::as_bytes));
            watched_profile_vars.push(name);
        }
    }

    let mut cargo_identity_vars: Vec<_> = env::vars()
        .filter(|(name, _)| {
            (name.starts_with("CARGO_PROFILE_") || name.starts_with("CARGO_FEATURE_"))
                && !watched_profile_vars.contains(name)
        })
        .collect();
    cargo_identity_vars.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    for (name, value) in cargo_identity_vars {
        println!("cargo:rerun-if-env-changed={name}");
        push_field(&mut payload, &name, value.as_bytes());
    }

    let rustc = required_env("RUSTC");
    add_executable_identity(&mut payload, "RUSTC", &rustc);
    for (name, value) in [
        ("RUSTC_WRAPPER", rustc_wrapper.as_deref()),
        (
            "RUSTC_WORKSPACE_WRAPPER",
            rustc_workspace_wrapper.as_deref(),
        ),
    ] {
        if let Some(value) = value.filter(|value| !value.is_empty()) {
            add_executable_identity(&mut payload, name, value);
        }
    }
    let output = Command::new(rustc)
        .arg("-vV")
        .output()
        .unwrap_or_else(|error| panic!("cannot execute RUSTC for GEMM build identity: {error}"));
    assert!(
        output.status.success(),
        "RUSTC -vV failed while constructing GEMM build identity: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    push_field(&mut payload, "rustc-vV", &output.stdout);

    let manifest_dir = PathBuf::from(required_env("CARGO_MANIFEST_DIR"));
    let workspace_root = manifest_dir.join("../..");
    read_identity_file(
        &mut payload,
        "workspace-Cargo.toml",
        &workspace_root.join("Cargo.toml"),
    );
    read_identity_file(
        &mut payload,
        "workspace-Cargo.lock",
        &workspace_root.join("Cargo.lock"),
    );
    let constellation_path = workspace_root.join("constellation.lock");
    println!("cargo:rerun-if-changed={}", constellation_path.display());
    let constellation_lock = read_required_file(&constellation_path);
    add_asupersync_identity(&mut payload, &workspace_root, &constellation_lock);
    add_source_closure(&mut payload, &workspace_root);

    let fingerprint = gemm_build_fingerprint(&GemmBuildIdentityInput {
        canonical_payload: &payload,
    })
    .to_hex();
    println!("cargo:rustc-env=FS_LA_GEMM_BUILD_FINGERPRINT={fingerprint}");
}
