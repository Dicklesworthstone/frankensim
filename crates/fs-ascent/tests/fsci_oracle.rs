//! Bead ijil (c): FrankenScipy CROSS-VALIDATION — the §12 oracle
//! contract on shared problems. API-CHECK OUTCOME (frankenscipy 0.1.0
//! per constellation.lock): `fsci_opt::minimize`/`slsqp` accept NO
//! general constraint callbacks (bounds/penalty only), so the oracle
//! pairing is: UNCONSTRAINED parity vs `minimize(Bfgs/LBfgsB)` on
//! shared smooth fixtures (including fsci's own Rosenbrock), and
//! CONSTRAINED parity vs `differential_evolution_constrained`
//! (penalty-based, seeded) within documented tolerances.

use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{LbfgsState, StopRule, augmented_lagrangian, interior_point, sqp};
use fsci_opt::{
    DifferentialEvolutionOptions, MinimizeOptions, OptimizeMethod, differential_evolution_constrained,
    minimize, rosen, rosen_der,
};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

#[test]
fn unconstrained_parity_on_rosenbrock() {
    // fs-ascent L-BFGS vs fsci minimize on fsci's OWN rosen fixture.
    let x0 = [-1.2f64, 1.0, 0.8, -0.5];
    let mut fg = |x: &[f64]| -> (f64, Vec<f64>) { (rosen(x), rosen_der(x)) };
    let mut st = LbfgsState::new(&x0, 10, &mut fg);
    let rep = st.run(&mut fg, &StopRule::GradNorm(1e-9), 2000);
    let ours = st.x.clone();
    let _ = rep;
    for method in [OptimizeMethod::Bfgs, OptimizeMethod::LBfgsB] {
        let res = minimize(
            rosen,
            &x0,
            MinimizeOptions {
                method: Some(method),
                tol: Some(1e-12),
                maxiter: Some(5000),
                ..MinimizeOptions::default()
            },
        )
        .expect("fsci minimize runs");
        let xdev = ours
            .iter()
            .zip(&res.x)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(
            xdev < 1e-4,
            "{method:?}: fs-ascent vs fsci optimum deviate: {xdev:.2e} (ours {ours:?} vs {:?})",
            res.x
        );
    }
    let f_ours = rosen(&ours);
    verdict(
        "ijil-fsci-unconstrained",
        f_ours < 1e-12,
        &format!(
            "Rosenbrock n=4: fs-ascent L-BFGS f*={f_ours:.2e}, optimum matches fsci Bfgs AND LBfgsB within 1e-4"
        ),
    );
}

#[test]
fn constrained_parity_vs_de_oracle() {
    // Shared problem: minimize (x−2)² + (y−1)² s.t. x + y = 2,
    // x ≤ 1.2 — all three fs-ascent engines vs the fsci DE oracle.
    // DE is penalty-based and stochastic-but-seeded: documented
    // tolerance 2e-2 on x*, 1e-2 on f*.
    let mk_fg = || {
        |x: &[f64]| -> (f64, Vec<f64>) {
            (
                (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2),
                vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)],
            )
        }
    };
    let ce = |x: &[f64]| vec![x[0] + x[1] - 2.0];
    let ce_jt = |_: &[f64], w: &[f64]| vec![w[0], w[0]];
    let ci = |x: &[f64]| vec![x[0] - 1.2];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0], 0.0];

    let de = differential_evolution_constrained(
        |x: &[f64]| (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2),
        &[(-3.0, 3.0), (-3.0, 3.0)],
        |x: &[f64]| {
            // Violation: equality as |c| and inequality as max(0, c).
            (x[0] + x[1] - 2.0).abs().max(0.0) + (x[0] - 1.2).max(0.0)
        },
        DifferentialEvolutionOptions {
            seed: Some(42),
            ..DifferentialEvolutionOptions::default()
        },
    )
    .expect("fsci DE runs");

    let mut fg1 = mk_fg();
    let mut p1 = ConstrainedProblem {
        fg: &mut fg1,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };
    let al = augmented_lagrangian(&mut p1, &[0.0, 0.0], 1e-7, 40);
    let mut fg2 = mk_fg();
    let mut p2 = ConstrainedProblem {
        fg: &mut fg2,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };
    let ip = interior_point(&mut p2, &[0.0, 0.0], 1e-6, 60);
    let mut fg3 = mk_fg();
    let mut p3 = ConstrainedProblem {
        fg: &mut fg3,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };
    let sq = sqp(&mut p3, &[0.0, 0.0], 1e-7, 60);

    let mut worst_x = 0.0f64;
    let mut worst_f = 0.0f64;
    for (label, x, f) in [
        ("AL", &al.x, al.f),
        ("IP", &ip.x, ip.f),
        ("SQP", &sq.x, sq.f),
    ] {
        let xdev = x
            .iter()
            .zip(&de.x)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        let fdev = (f - de.fun).abs();
        worst_x = worst_x.max(xdev);
        worst_f = worst_f.max(fdev);
        assert!(
            xdev < 2e-2 && fdev < 1e-2,
            "{label} vs DE oracle: xdev {xdev:.3e}, fdev {fdev:.3e} (x {x:?} vs {:?})",
            de.x
        );
    }
    verdict(
        "ijil-fsci-constrained",
        true,
        &format!(
            "AL/IP/SQP vs fsci differential_evolution_constrained (seeded): worst xdev {worst_x:.2e}, worst fdev {worst_f:.2e} (documented DE tolerance 2e-2/1e-2)"
        ),
    );
}
