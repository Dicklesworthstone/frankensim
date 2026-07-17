//! Fork-steering conformance (the qlvf bead, lane a): world-forking
//! round-trip — fork, diverge, and BOTH branches replay bitwise from
//! their lineage; the parent is untouched by the fork; steering
//! events are ledger-ready. Aggregate verdicts carry explicit input
//! seeds; measurements use validated fs-obs Custom events.

use fs_dfo::steer::SteeredStudy;

fn verdict(case: &str, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-dfo/steer", case);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-dfo/steer".to_string(),
            case: case.to_string(),
            pass: true,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("steering verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("steering verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn measurement(case: &str, name: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new("fs-dfo/steer", case);
    let line = emitter
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: name.to_string(),
                json,
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("steering measurement must use the fs-obs wire schema");
    println!("{line}");
}

/// A two-objective toy: f1 = ||x - a||^2, f2 = ||x + a||^2.
fn objectives(x: &[f64]) -> Vec<f64> {
    let f1: f64 = x.iter().map(|v| (v - 1.0) * (v - 1.0)).sum();
    let f2: f64 = x.iter().map(|v| (v + 1.0) * (v + 1.0)).sum();
    vec![f1, f2]
}

#[test]
fn st_001_fork_diverge_replay_bitwise() {
    const SEED: u64 = 0x5eed;
    let mut obj = objectives;
    let bounds = (-2.0, 2.0);
    // Grow the trunk 5 generations.
    let mut trunk = SteeredStudy::start(&mut obj, 3, bounds, 16, 2, SEED);
    trunk.advance(&mut obj, bounds, 5);
    let trunk_at_fork = trunk.fingerprint();
    // FORK toward objective 1; the PARENT must be untouched.
    let mut child = trunk.fork(vec![0.9, 0.1]);
    assert_eq!(
        trunk.fingerprint(),
        trunk_at_fork,
        "forking never mutates the parent (P9: world-forking)"
    );
    assert_eq!(child.lineage.len(), 1, "the steering event is recorded");
    assert!(child.lineage[0].to_json().contains("\"at\":5"));
    // Diverge: 6 more generations each under different weights.
    trunk.advance(&mut obj, bounds, 6);
    child.advance(&mut obj, bounds, 6);
    // The branches genuinely diverged toward their weights.
    let trunk_best = trunk.best_score();
    let child_best = child.best_score();
    assert_ne!(trunk.fingerprint(), child.fingerprint(), "branches differ");
    // REPLAY both branches from scratch by re-applying the lineage.
    let mut trunk_replay = SteeredStudy::start(&mut obj, 3, bounds, 16, 2, SEED);
    trunk_replay.advance(&mut obj, bounds, 5);
    let mut child_replay = trunk_replay.fork(vec![0.9, 0.1]);
    trunk_replay.advance(&mut obj, bounds, 6);
    child_replay.advance(&mut obj, bounds, 6);
    assert_eq!(
        trunk.fingerprint(),
        trunk_replay.fingerprint(),
        "the trunk replays bitwise"
    );
    assert_eq!(
        child.fingerprint(),
        child_replay.fingerprint(),
        "the fork replays bitwise from its lineage"
    );
    measurement(
        "st-001/fork",
        "dfo-fork-steering",
        format!(
            "{{\"seed\":{SEED},\"trunk_best\":{trunk_best:.4},\"child_best\":{child_best:.4}}}"
        ),
    );
    verdict(
        "st-001",
        "fork leaves the parent untouched, records the ledger-ready steering event, and \
         BOTH diverged branches replay bitwise from their lineage (G5)",
        SEED,
    );
}

#[test]
fn st_002_steering_actually_steers() {
    const SEED: u64 = 0xace;
    let mut obj = objectives;
    let bounds = (-2.0, 2.0);
    let mut study = SteeredStudy::start(&mut obj, 3, bounds, 24, 2, SEED);
    study.advance(&mut obj, bounds, 4);
    // Two forks with opposite priorities.
    let mut toward_1 = study.fork(vec![0.95, 0.05]);
    let mut toward_2 = study.fork(vec![0.05, 0.95]);
    toward_1.advance(&mut obj, bounds, 12);
    toward_2.advance(&mut obj, bounds, 12);
    // Each branch's best individual leans toward its favored optimum
    // (x -> +1 for f1, x -> -1 for f2).
    let lean = |s: &SteeredStudy| -> f64 {
        let best = s
            .state
            .population
            .iter()
            .min_by(|p, q| {
                let sp: f64 = p.f.iter().zip(&s.state.weights).map(|(f, w)| f * w).sum();
                let sq: f64 = q.f.iter().zip(&s.state.weights).map(|(f, w)| f * w).sum();
                sp.total_cmp(&sq)
            })
            .expect("pop");
        best.x.iter().sum::<f64>() / best.x.len() as f64
    };
    let (l1, l2) = (lean(&toward_1), lean(&toward_2));
    measurement(
        "st-002/steering",
        "dfo-opposite-steering",
        format!("{{\"seed\":{SEED},\"toward_f1_mean_x\":{l1:.3},\"toward_f2_mean_x\":{l2:.3}}}"),
    );
    assert!(
        l1 > 0.5 && l2 < -0.5,
        "re-weighting steers the branches to opposite optima: {l1:.2} vs {l2:.2}"
    );
    verdict(
        "st-002",
        "opposite re-weightings drive sibling forks to opposite ends of the Pareto \
         landscape — steering is real, and every step of it is a ledgered op",
        SEED,
    );
}

#[test]
#[should_panic(expected = "steered study needs a positive population")]
fn steered_study_rejects_an_empty_population() {
    // Regression: `advance` selects parents via `... % population.len()`; an
    // empty population divides by zero. Fail closed at construction.
    let mut obj = |x: &[f64]| vec![x[0], x[1]];
    let _ = SteeredStudy::start(&mut obj, 3, (0.0, 1.0), 0, 2, 0x5eed);
}
