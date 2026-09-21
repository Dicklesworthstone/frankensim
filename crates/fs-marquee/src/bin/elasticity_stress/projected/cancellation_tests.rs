//! Exercise the command's real entry/publication and on-disk recovery paths.
use super::*;
use fs_topols::evaluated::DesignEvaluationStage;
use fs_topols::projected::ProjectedSetupStage;
use fs_topols::projected_stress::ProjectedStressSetupStage;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!("fs-stress-cancel-{}-{}",
            std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
}

fn expired() -> Instant { Instant::now().checked_sub(Duration::from_secs(2)).unwrap() }

fn optimizer() -> (ProjectedStressOptimizer, checkpoint::Policy) {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed: Vec<_> = geometry.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    let settings = OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6,
        move_cells: 0.1, nucleation_period: 0, ..OptimizeSettings::default() };
    let policy = checkpoint::Policy {
        fixed,
        projection: VolumeProjectionSettings {
            target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
        },
        controls: ProjectedSettings { max_candidates: 8, poll_iters: 1, ..ProjectedSettings::default() },
    };
    let area = ProjectedOptimizer::new(geometry, Cantilever { load: 1.0, band: 0.125 },
        settings, policy.fixed.clone(), policy.projection, policy.controls).unwrap();
    (ProjectedStressOptimizer::new(area, SampledStressLimit::new(1e12, 0.0).unwrap()).unwrap(), policy)
}

#[test]
fn expired_startup_and_resume_budgets_publish_no_study() {
    let scratch = Scratch::new();
    let fresh = scratch.0.join("fresh");
    let mut args = vec![fresh.to_string_lossy().into_owned()];
    args.extend(["1e12", "3", "2", "0.6", "8", "0", "1", "--checkpoint"].map(str::to_string));
    assert_eq!(run_at(&args, expired()).unwrap(), 6);
    assert!(!fresh.exists());

    let resumed = scratch.0.join("resumed");
    let args = vec!["--resume".to_string(), scratch.0.join("not-read.fscp").to_string_lossy().into_owned(),
        resumed.to_string_lossy().into_owned(), "--wall-seconds".to_string(), "1".to_string()];
    // The expired budget stops before even opening the missing source file.
    assert_eq!(run_at(&args, expired()).unwrap(), 6);
    assert!(!resumed.exists());
}

#[test]
fn budget_expiring_after_admission_still_prevents_first_publication() {
    let scratch = Scratch::new();
    let (optimizer, policy) = optimizer();
    let output = scratch.0.join("no-publication");
    let result = run_optimizer(&output, optimizer, policy, checkpoint::Options::default(),
        None, expired(), 1, false).unwrap();
    assert_eq!(result, 6);
    assert!(!output.exists());
}

#[test]
fn interrupted_disk_recovery_keeps_exact_source_and_retry_matches_uninterrupted_tail() {
    let scratch = Scratch::new();
    let (mut original, policy) = optimizer();
    assert!(matches!(original.advance_one().unwrap().progress, ProjectedProgress::Accepted(_)));
    // This test isolates the real disk codec and PDE replay. The CLI obtains
    // its actual executable fingerprint before calling the same functions.
    let fingerprint = "0123456789abcdef".repeat(4);
    let input = scratch.0.join("accepted.fscp");
    checkpoint::save(&input, &original, &policy, &fingerprint).unwrap();
    let bytes = std::fs::read(&input).unwrap();
    for stop in 0..4 {
        let result = checkpoint::load_controlled(&input, &fingerprint, |stage| {
            let stop_here = match stop {
                0 => matches!(stage, ProjectedStressSetupStage::Area(ProjectedSetupStage::Evaluation(DesignEvaluationStage::Solve(n))) if n > 0),
                1 => matches!(stage, ProjectedStressSetupStage::Stress(DesignEvaluationStage::Solve(n)) if n > 0),
                2 => matches!(stage, ProjectedStressSetupStage::Stress(DesignEvaluationStage::StressCell(32))),
                _ => matches!(stage, ProjectedStressSetupStage::Publish),
            };
            if stop_here { ControlFlow::Break(format!("disk stop {stop}")) }
            else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break(ref why) if *why == format!("disk stop {stop}")));
        assert_eq!(std::fs::read(&input).unwrap(), bytes);
    }
    let ControlFlow::Continue((mut resumed, restored_policy)) = checkpoint::load_controlled(
        &input, &fingerprint, |_| ControlFlow::<()>::Continue(())).unwrap()
    else { panic!("uninterrupted recovery stopped") };
    let copy = scratch.0.join("recovered.fscp");
    checkpoint::save(&copy, &resumed, &restored_policy, &fingerprint).unwrap();
    assert_eq!(std::fs::read(copy).unwrap(), bytes);
    let resumed_update = resumed.advance_one().unwrap();
    let original_update = original.advance_one().unwrap();
    assert_eq!(format!("{resumed_update:?}"), format!("{original_update:?}"));
    assert_eq!(resumed.current(), original.current());
    let a = scratch.0.join("resumed-tail.fscp");
    let b = scratch.0.join("original-tail.fscp");
    checkpoint::save(&a, &resumed, &restored_policy, &fingerprint).unwrap();
    checkpoint::save(&b, &original, &policy, &fingerprint).unwrap();
    assert_eq!(std::fs::read(a).unwrap(), std::fs::read(b).unwrap());
    assert_eq!(std::fs::read(input).unwrap(), bytes);
}
