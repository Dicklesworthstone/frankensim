//! Success-rate regression battery (grew out of the bring-up probe —
//! filename kept for history; deletion needs user permission per
//! RULE 1): large-population CMA-ES on rastrigin(5) must solve the
//! global basin in a majority of seeds, and FAILED runs must terminate
//! via the stagnation stop instead of burning their full budget (the
//! property BIPOP's restart ladder depends on).

use fs_dfo::{CmaParams, cmaes};

fn rastrigin(x: &[f64]) -> f64 {
    let a = 10.0f64;
    x.iter()
        .map(|&t| a + t.mul_add(t, -a * fs_math::det::cos(2.0 * std::f64::consts::PI * t)))
        .sum()
}

#[test]
fn large_population_success_rate_and_stagnation_stop() {
    let mut successes = 0usize;
    for seed in 1u64..=5 {
        let mut f = |x: &[f64]| rastrigin(x);
        let p = CmaParams {
            lambda: 150,
            sigma0: 3.0,
            max_evals: 120_000,
            f_target: 1e-8,
            eigen_interval: 1,
        };
        let rep = cmaes(&mut f, &[3.0; 5], &p, seed);
        if rep.converged {
            successes += 1;
        } else {
            // The stagnation stop must have fired well short of budget.
            assert!(
                rep.evals < 100_000,
                "failed run must stop on stagnation, not burn budget: {} evals",
                rep.evals
            );
        }
    }
    assert!(
        successes >= 3,
        "lambda=150 must solve rastrigin(5) in a majority of seeds: {successes}/5"
    );
    println!(
        "{{\"suite\":\"fs-dfo\",\"case\":\"success-rate\",\"verdict\":\"pass\",\"detail\":\"{successes}/5 seeds converged at lambda=150\"}}"
    );
}
