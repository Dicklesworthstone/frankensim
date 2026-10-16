//! Success-rate regression battery (grew out of the bring-up probe —
//! filename kept for history; deletion needs user permission per
//! RULE 1): large-population CMA-ES on rastrigin(5) must solve the
//! global basin in a majority of seeds, and FAILED runs must terminate
//! via the stagnation stop instead of burning their full budget (the
//! property BIPOP's restart ladder depends on).

use fs_dfo::{CmaParams, cmaes};

const SUITE: &str = "fs-dfo";
const FIXED_INPUT_SEED: u64 = 0;
const FIRST_INPUT_SEED: u64 = 1;
const LAST_INPUT_SEED: u64 = 5;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("DFO success-rate verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("DFO success-rate verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn rastrigin(x: &[f64]) -> f64 {
    let a = 10.0f64;
    x.iter()
        .map(|&t| a + t.mul_add(t, -a * fs_math::det::cos(2.0 * std::f64::consts::PI * t)))
        .sum()
}

#[test]
fn large_population_success_rate_and_stagnation_stop() {
    let mut successes = 0usize;
    for seed in FIRST_INPUT_SEED..=LAST_INPUT_SEED {
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
    verdict(
        "success-rate",
        true,
        &format!(
            "{successes}/5 seeds converged at lambda=150; optimizer input roots \
             {FIRST_INPUT_SEED}..={LAST_INPUT_SEED}; composite aggregate seed zero"
        ),
        FIXED_INPUT_SEED,
    );
}
