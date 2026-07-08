//! Checkpoint tolerance-calibration probe: documents WHERE the accept
//! economy opens on the P1/24-cell elliptic family (the basis for the
//! phase-1 gate's 0.05 realistic / 0.02 hostile tolerances) and asserts
//! the acceptance rate is monotone in tolerance.
#![cfg(feature = "flywheel-e2e")]
use fs_verify::economics::{DriftGuard, EconDecision, run_speculative};
use fs_verify::fem1d::{MmsProblem, Poly};
use fs_verify::zoo::{
    CoarseRungProlongation, NeighborExtrapolation, Registry, SpeculationQuery, ZooTelemetry,
};

fn family(theta: f64, cells: usize) -> MmsProblem {
    let base = Poly(vec![0.0, theta, -theta, 0.25 * theta, -0.25 * theta]);
    #[allow(clippy::cast_precision_loss)]
    let mesh: Vec<f64> = (0..=cells).map(|k| k as f64 / cells as f64).collect();
    MmsProblem::new("elliptic-theta", base, mesh)
}

fn accepts_at(tol: f64) -> usize {
    let cache: Vec<(f64, Vec<f64>, Option<Vec<f64>>)> = [0.8, 1.0, 1.2, 1.4]
        .iter()
        .map(|&t| (t, fs_verify::fem1d::solve_p1(&family(t, 24)), None))
        .collect();
    let mut registry = Registry::new();
    registry.register(Box::new(NeighborExtrapolation { cache }));
    registry.register(Box::new(CoarseRungProlongation));
    let mut tele = ZooTelemetry::default();
    let mut guard = DriftGuard::default();
    let mut acc = 0;
    for k in 0..20 {
        let theta = 0.85 + 0.03 * f64::from(k as u8);
        let q = SpeculationQuery {
            problem: family(theta, 24),
            theta,
            tolerance: tol,
            regime: "r".to_string(),
        };
        if matches!(
            run_speculative(&q, &registry, &mut tele, &mut guard, 200),
            EconDecision::AcceptedOutright { .. }
        ) {
            acc += 1;
        }
    }
    acc
}

#[test]
fn calibration_monotone_in_tolerance() {
    let curve: Vec<usize> = [0.02, 0.05, 0.08].iter().map(|&t| accepts_at(t)).collect();
    println!(
        "{{\"metric\":\"checkpoint-calibration\",\"tols\":[0.02,0.05,0.08],\
         \"accepts\":{curve:?}}}"
    );
    assert!(curve[0] <= curve[1] && curve[1] <= curve[2], "monotone: {curve:?}");
    assert_eq!(curve[0], 0, "hostile tolerance rejects everything");
    assert!(curve[1] >= 6, "the realistic tolerance opens the economy");
}
