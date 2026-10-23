//! The default drum render refuses the recording-fitted preset heads (bead .370, DONE-WHEN 1);
//! `--preset-heads` renders them on purpose and says so. Runs the binary on the CPU backend.

use std::process::Command;

#[test]
fn drums_refuse_the_fitted_preset_heads_unless_asked() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("render-gpu-preset-heads-{}-{stamp}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let hits = dir.join("hits.csv");
    std::fs::write(&hits, "time_s,speed_m_s\n0.001,2.0\n").unwrap();
    let out = dir.join("tom.wav");
    let run = |extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_render-gpu"))
            .args(["tom-12", hits.to_str().unwrap(), out.to_str().unwrap(), "0.02", "--backend", "cpu"])
            .args(extra)
            .output()
            .unwrap()
    };
    let refused = run(&[]);
    let stderr = String::from_utf8_lossy(&refused.stderr);
    assert!(!refused.status.success(), "{stderr}");
    assert!(stderr.contains("fitted to recordings") && stderr.contains("--preset-heads"), "{stderr}");
    assert!(!out.exists(), "a refused render writes nothing");

    let fitted = run(&["--preset-heads"]);
    let stderr = String::from_utf8_lossy(&fitted.stderr);
    assert!(fitted.status.success(), "{stderr}");
    assert!(stderr.contains("heads: PRESET, recording-fitted T60 (opt-in --preset-heads)"), "{stderr}");
    assert!(out.exists());
}
