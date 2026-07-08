//! Planner calibration probe: uniform-mesh certified bounds vs cells on
//! the steep family (documents the kill-test tolerance choice).
#![cfg(feature = "ladder-planner")]
use fs_ir::planner::ProblemFamily;
use fs_verify::estimator::verify;
use fs_verify::fem1d::{Poly, solve_p1};

#[test]
fn uniform_bound_curve() {
    let mut c = vec![0.0; 11];
    c[1] = 0.2;
    c[2] = -0.2;
    c[9] = 1.0;
    c[10] = -1.0;
    let family = ProblemFamily {
        base: Poly(c),
        kernel: "steep".to_string(),
    };
    for cells in [12, 24, 48, 96, 192, 384] {
        let mesh: Vec<f64> = (0..=cells)
            .map(|k| f64::from(k) / f64::from(cells))
            .collect();
        let p = family.at(1.0, mesh);
        let u = solve_p1(&p);
        let rep = verify(&p, &u, 1e-9);
        println!("cells={cells} bound={:.4e}", rep.bound.hi);
    }
}

#[test]
fn trace_kill_run() {
    use fs_ir::planner::{CostTable, MemCache, PlanOutcome, plan};
    let mut c = vec![0.0; 11];
    c[1] = 0.2;
    c[2] = -0.2;
    c[9] = 1.0;
    c[10] = -1.0;
    let family = ProblemFamily {
        base: Poly(c),
        kernel: "steep".to_string(),
    };
    let out = plan(
        &family,
        1.0,
        6e-3,
        100_000.0,
        &[12, 24, 48, 96],
        &mut MemCache::default(),
        &mut CostTable::new(200.0),
    );
    if let PlanOutcome::Discharged { ops, cost, .. } = out {
        for o in &ops {
            println!(
                "op={} cost={} bound={:.3e}",
                o.op.name(),
                o.cost,
                o.bound_after
            );
        }
        println!("total={cost}");
    }
}
