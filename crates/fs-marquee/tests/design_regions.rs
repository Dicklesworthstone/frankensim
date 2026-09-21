//! Authored solid/void geometry through real projected, resume and refine commands.
use fs_topols::GridSdf;
use fs_topols::design_regions::{DesignPhase, DesignRegion, prepare_design_regions};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);
fn workspace() -> PathBuf {
    let path = std::env::temp_dir().join(format!("fs-design-regions-{}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    std::fs::write(path.join("loads.csv"), "right,0.375,0.625,0,-1,1\n").unwrap();
    std::fs::write(path.join("regions.csv"), "void,0.5,0.5,0.625,0.625,0.02\n").unwrap();
    path
}
fn binary() -> Command { Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust")) }
fn start(root: &Path, output: &str, regions: bool, pause: &str, budget: &str) -> Output {
    let mut cmd = binary();
    cmd.arg("--projected").arg(root.join(output)).arg(root.join("loads.csv"))
        .args(["3", "1", "0.7", "16", "sum", budget, "--checkpoint", "--pause-after", pause,
            "--stress-limit", "100000000"]);
    if regions { cmd.arg("--design-regions").arg(root.join("regions.csv")); }
    cmd.output().unwrap()
}
fn field(path: &Path, n: usize) -> GridSdf {
    let text = std::fs::read_to_string(path).unwrap();
    let values: Vec<f64> = text.lines().skip(1)
        .map(|line| line.split(',').nth(2).unwrap().parse().unwrap()).collect();
    let mut out = GridSdf::from_fn(n, &|_, _| 0.0);
    assert_eq!(out.nodes().len(), values.len());
    out.nodes_mut().copy_from_slice(&values);
    out
}
fn check_void(phi: &GridSdf) {
    for i in 0..=10 {
        for j in 0..=10 {
            assert!(phi.value_at([0.5 + f64::from(i) / 80.0,
                0.5 + f64::from(j) / 80.0]) > 0.0);
        }
    }
}

#[test]
fn authored_hole_survives_an_actual_update_restart_and_finer_baseline() {
    let root = workspace();
    let result = start(&root, "paused", true, "0", "17");
    assert_eq!(result.status.code(), Some(14), "{}", String::from_utf8_lossy(&result.stderr));
    let source = root.join("paused/checkpoint.fscp");
    let source_bytes = std::fs::read(&source).unwrap();
    let region_bytes = std::fs::read(root.join("regions.csv")).unwrap();
    let input = field(&root.join("paused/input-level-set.csv"), 8);
    assert!(input.node(4, 4) < 0.0);
    let baseline = field(&root.join("paused/baseline-level-set.csv"), 8);
    check_void(&baseline);
    let authored = field(&root.join("paused/design-region-level-set.csv"), 8);
    let fixed: Vec<_> = input.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8 || i / 9 == 0 || i / 9 == 8).collect();
    let expected = prepare_design_regions(&input, &fixed, &[
        DesignRegion::new(DesignPhase::Void, [0.5, 0.5], [0.625, 0.625], 0.02).unwrap(),
    ]).unwrap();
    assert_eq!(authored.nodes(), expected.geometry.nodes());
    let result = binary().args(["--projected", "--resume"]).arg(&source).arg(root.join("continued"))
        .args(["--recovery-solves", "2"]).output().unwrap();
    assert_eq!(result.status.code(), Some(0), "{}", String::from_utf8_lossy(&result.stderr));
    assert_eq!(std::fs::read_to_string(root.join("continued/trajectory.jsonl")).unwrap().lines().count(), 1,
        "must execute an accepted update, not simply re-export the paused baseline");
    check_void(&field(&root.join("continued/level-set.csv"), 8));
    let result = binary().args(["--projected", "--refine"])
        .arg(root.join("continued/checkpoint.fscp")).arg(root.join("fine"))
        .args(["--updates", "1", "--max-solves", "1", "--recovery-solves", "2", "--pause-after", "0"])
        .output().unwrap();
    assert_eq!(result.status.code(), Some(14), "{}", String::from_utf8_lossy(&result.stderr));
    check_void(&field(&root.join("fine/baseline-level-set.csv"), 16));
    assert_eq!(std::fs::read(&source).unwrap(), source_bytes);
    assert_eq!(std::fs::read(root.join("regions.csv")).unwrap(), region_bytes);
}

#[test]
fn conflicts_and_impossible_area_publish_no_baseline_or_checkpoint() {
    let root = workspace();
    for (name, csv) in [
        ("phase-conflict", "material,0.25,0.25,0.5,0.5,0.1\nvoid,0.5,0.25,0.625,0.5,0.1\n"),
        ("boundary-conflict", "void,0,0.375,0.125,0.625,0.1\n"),
        // Together with the frozen boundary, this fixes every node of the
        // original ~0.84-area strip; the requested 0.7 area is unattainable.
        ("area-infeasible", "material,0.125,0.125,0.875,0.875,0.04\n"),
    ] {
        std::fs::write(root.join("regions.csv"), csv).unwrap();
        let result = start(&root, name, true, "0", "1");
        assert_eq!(result.status.code(), Some(1));
        assert!(result.stdout.is_empty());
        assert!(!root.join(name).exists());
        assert_eq!(std::fs::read_to_string(root.join("regions.csv")).unwrap(), csv);
    }
}

#[test]
fn bounded_inputs_duplicate_options_and_existing_outputs_refuse() {
    let root = workspace();
    let result = binary().arg("--projected").arg(root.join("no-input"))
        .arg(root.join("missing-loads.csv")).args(["--design-regions", "a", "--design-regions", "b"])
        .output().unwrap();
    assert!(String::from_utf8_lossy(&result.stderr).contains("duplicate --design-regions"));
    assert!(!root.join("no-input").exists());
    std::fs::write(root.join("regions.csv"), vec![b' '; 1_048_577]).unwrap();
    let result = start(&root, "oversized", true, "0", "1");
    assert!(String::from_utf8_lossy(&result.stderr).contains("exceeds 1 MiB"));
    assert!(!root.join("oversized").exists());
    std::fs::create_dir(root.join("occupied")).unwrap();
    std::fs::write(root.join("occupied/keep"), b"original").unwrap();
    let result = start(&root, "occupied", true, "0", "1");
    assert!(String::from_utf8_lossy(&result.stderr).contains("refusing to overwrite"));
    assert_eq!(std::fs::read(root.join("occupied/keep")).unwrap(), b"original");
}

#[test]
fn absent_region_option_retains_original_outputs_and_policy_cannot_change_on_resume() {
    let root = workspace();
    let result = start(&root, "plain", false, "0", "1");
    assert_eq!(result.status.code(), Some(14));
    let summary = std::fs::read_to_string(root.join("plain/summary.json")).unwrap();
    assert!(!summary.contains("design_regions"));
    assert!(!root.join("plain/design-regions.csv").exists());
    let result = binary().args(["--projected", "--resume"]).arg(root.join("plain/checkpoint.fscp"))
        .arg(root.join("replace-policy")).args(["--recovery-solves", "2", "--design-regions"])
        .arg(root.join("regions.csv")).output().unwrap();
    assert_eq!(result.status.code(), Some(1));
    assert!(!root.join("replace-policy").exists());
}
