//! Chance-constraint conformance (the qlvf bead, lane b): the
//! Gaussian toy with a KNOWN analytic feasible boundary, solved with
//! anytime-stopped probability estimates feeding the outer loop.

use fs_uq::chance_constrained_min;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-uq/chance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// Deterministic standard-normal germ (12-uniform sum).
fn gauss(seed: u64, k: u64) -> f64 {
    let unit = |j: u64| -> f64 {
        let mut z = seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(k * 12 + j + 1);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..12).map(unit).sum::<f64>() - 6.0
}

#[test]
fn cc_001_gaussian_toy_hits_the_analytic_boundary() {
    // min x  s.t.  P(xi <= x) >= 1 - alpha,  xi ~ N(mu, sigma).
    // Analytic optimum: x* = mu + sigma * z_{1-alpha}.
    // alpha = 0.1 -> z_0.9 = 1.2816.
    let (mu, sigma, alpha) = (2.0, 0.5, 0.10);
    let x_star = mu + sigma * 1.281_551_6;
    let solution = chance_constrained_min(
        |x| x,                                      // minimize x itself
        |x, i| mu + sigma * gauss(0xc0ffee, i) - x, // g = xi - x
        alpha,
        (mu - 1.0, mu + 3.0),
        4,
        0.05,
        0.02,
        50_000,
    )
    .expect("the bounded grid contains a simultaneously admitted design");
    let (x, p_est, spent) = (solution.x, solution.probability, solution.samples);
    println!(
        "{{\"metric\":\"chance\",\"x\":{x:.4},\"x_star\":{x_star:.4},\"p_at_x\":{p_est:.3},\
         \"samples\":{spent}}}"
    );
    // The anytime CS uses its LOWER bound for feasibility, so the
    // solution sits AT or conservatively ABOVE the analytic boundary
    // (never below it — feasibility is what the validity buys).
    assert!(
        x >= x_star - 0.05,
        "never infeasible-by-optimism: {x:.4} vs boundary {x_star:.4}"
    );
    assert!(
        x <= x_star + 0.30,
        "and not absurdly conservative: {x:.4} vs {x_star:.4}"
    );
    assert!(
        solution.probability_lo >= 1.0 - alpha,
        "the simultaneous lower bound backs feasibility"
    );
    verdict(
        "cc-001",
        "the chance-constrained minimum lands at-or-conservatively-above the analytic \
         Gaussian quantile boundary, never below — anytime-stopped estimates with the \
         CS lower bound enforcing feasibility",
    );
}

#[test]
fn cc_002_multiplicity_and_work_admission_fail_closed() {
    let no_feasible =
        chance_constrained_min(|x| x, |_x, _i| 1.0, 0.1, (0.0, 1.0), 1, 0.05, 0.1, 100);
    assert!(no_feasible.is_err());
    assert!(
        std::panic::catch_unwind(|| {
            let _ = chance_constrained_min(|x| x, |_x, _i| -1.0, 0.1, (0.0, 1.0), 65, 0.05, 0.1, 1);
        })
        .is_err()
    );
    assert!(
        std::panic::catch_unwind(|| {
            let _ =
                chance_constrained_min(|x| x, |_x, _i| f64::NAN, 0.1, (0.0, 1.0), 1, 0.05, 0.1, 1);
        })
        .is_err()
    );
}
