//! Execute the real binary dispatch, not only the library producers.
use std::process::Command;

#[test]
fn binary_admits_a_file_driven_mixed_solid_network() {
    let request = concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/mixed-slab.json");
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", request]).output().expect("start binary");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let result = String::from_utf8(output.stdout).unwrap();
    assert!(result.contains("\"schema\":\"frankensim.cooling-network.result.v1\""));
    assert!(result.contains("\"solid_temperatures_k\":["));
    assert!(result.contains("\"dmean_dinlet_k\":["));
    assert!(output.stderr.is_empty());
}

#[test]
fn library_and_binary_execute_the_same_coupled_thermal_model() {
    let request = concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/mixed-slab.json");
    let library = fs_cli::run(["cooling-network", request, "--json"].map(str::to_owned));
    assert_eq!(library.exit_code, fs_cli::exit::SUCCESS, "{}", library.stderr);
    let binary = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["cooling-network", request, "--json"]).output().unwrap();
    assert_eq!(binary.status.code(), Some(i32::from(library.exit_code)));
    assert_eq!(binary.stdout, library.stdout.as_bytes());
    assert_eq!(binary.stderr, library.stderr.as_bytes());
    assert!(library.stdout.contains("\"solid_temperatures_k\":["));
    assert!(library.stdout.contains("\"dmean_dinlet_k\":["));
}

#[cfg(unix)]
#[test]
fn native_cooling_paths_survive_shared_dispatch_without_utf8_conversion() {
    use std::os::unix::ffi::OsStringExt;
    let directory = std::env::temp_dir().join(format!("fs-cli-native-path-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join(std::ffi::OsString::from_vec(b"model-\xff.json".to_vec()));
    std::fs::write(&path, include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/mixed-slab.json"))).unwrap();
    let arguments = ["--json".into(), "cooling-network".into(), path.into_os_string()];
    let library = fs_cli::run_os(arguments.clone());
    assert_eq!(library.exit_code, fs_cli::exit::SUCCESS, "{}", library.stderr);
    let binary = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(arguments).output().unwrap();
    assert!(binary.status.success(), "{}", String::from_utf8_lossy(&binary.stderr));
    assert_eq!(binary.stdout, library.stdout.as_bytes());
    assert_eq!(binary.stderr, library.stderr.as_bytes());
}

#[test]
fn binary_help_and_missing_request_have_stable_exit_classes() {
    let help = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["cooling-network", "--help"]).output().unwrap();
    assert!(help.status.success());
    assert!(String::from_utf8(help.stdout).unwrap().contains("cooling-network <request.json>"));
    let bad = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).output().unwrap();
    assert_eq!(bad.status.code(), Some(2));
    assert!(bad.stdout.is_empty());
}
