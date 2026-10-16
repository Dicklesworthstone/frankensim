//! Diagnosis-calibration regression (kept from the tfz.10 probe):
//! CG's RECURSIVE residual on an inconsistent singular system
//! DIVERGES rather than plateauing — which is why the Plateau fixture
//! in the battery is restarted-GMRES stagnation (true-residual
//! measurement) and why this trajectory must keep reading as
//! BudgetExhausted/Breakdown territory, never Plateau.
use fs_feec::{element_geometry, kuhn_cube};
use fs_rand::StreamKey;
use fs_solver::{CgState, CsrOp};
use fs_sparse::precond::IdentityPrecond;

#[test]
fn probe() {
    let (complex, positions) = kuhn_cube(4);
    let geo = element_geometry(&complex, &positions);
    let k_full = fs_feec::stiffness(
        &fs_feec::incidence_to_csr(&complex.d0()),
        &fs_feec::mass_matrix(&complex, &geo, 1),
    );
    let n = positions.len();
    let op = CsrOp::symmetric(k_full);
    let mut s = StreamKey {
        seed: 21,
        kernel: 0x501E,
        tile: 81,
    }
    .stream();
    let b: Vec<f64> = (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect();
    let mut st = CgState::new(&op, &IdentityPrecond, &b);
    let rep = st.run(&op, &IdentityPrecond, 1e-12, 500);
    let h = &rep.history;
    println!("len={} last={:e}", h.len(), h.last().unwrap());
    assert!(!rep.converged);
    assert!(
        *h.last().expect("history") > 1.0,
        "inconsistent singular CG should diverge, not settle"
    );
    assert_ne!(
        rep.diagnosis,
        Some(fs_solver::StallDiagnosis::Plateau),
        "a diverging trajectory must not read as Plateau"
    );
}
