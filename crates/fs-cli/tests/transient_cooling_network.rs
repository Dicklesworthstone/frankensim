use std::process::Command;

#[test]
fn transient_workload_reaches_the_actual_binary_dispatch() {
    let request = concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/transient-contact-pulse.json");
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", request]).output().expect("run binary");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let text = String::from_utf8(output.stdout).expect("UTF-8 JSON");
    assert!(text.contains("\"scheme\":\"backward-euler\""));
    assert!(text.contains("\"steps\":75"));
    assert!(text.contains("\"first_sampled_violation_s\":22"));
    assert!(text.contains("\"fan\":") && text.contains("\"speed_ratio\":1.5"));
    assert!(text.contains("\"contacts\":[{"));
}
