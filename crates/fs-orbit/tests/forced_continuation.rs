//! Original forced dynamics under independent state/equation/parameter units.

use fs_math::c64::C64;
use fs_orbit::{
    ContinuableProblem, ContinuationBudget, ContinuationScales, HbAnchor, HbBudget, HbScales,
    IslandStateJacobian, OrbitError, OrbitProblem, ScaledContinuationPoint,
    continue_branch_scaled_seeded, solve_hb_scaled_seeded,
};
use std::cell::Cell;

fn budget(harmonics: usize) -> ContinuationBudget {
    ContinuationBudget {
        max_steps: 500,
        initial_step: 0.025,
        min_step: 1e-6,
        hb: HbBudget {
            harmonics,
            max_newton: 35,
            tolerance: 2e-10,
            deflate_equilibrium: false,
        },
    }
}

struct Duffing {
    omega: f64,
    state: [f64; 2],
    equations: [f64; 2],
    parameter: f64,
    analytic: bool,
    island_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}

impl Duffing {
    fn scales(&self) -> ContinuationScales {
        ContinuationScales {
            hb: HbScales {
                states: self.state.to_vec(),
                equations: self.equations.to_vec(),
                angular_frequency: 1.0,
            },
            parameter: self.parameter,
        }
    }
    fn tangent(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        let q = x[0] / self.state[0];
        out.copy_from_slice(&[
            0.0,
            self.equations[0] / self.state[1],
            -self.equations[1] * (1.0 + 0.45 * q * q) / self.state[0],
            -0.1 * self.equations[1] / self.state[1],
        ]);
        true
    }
}

impl OrbitProblem for Duffing {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        self.island_calls.set(self.island_calls.get() + 1);
        let q = x[0] / self.state[0];
        let v = x[1] / self.state[1];
        out[0] = self.equations[0] * v;
        out[1] =
            self.equations[1] * (-q - 0.15 * q.powi(3) - 0.1 * v + 0.3 * (self.omega * t).cos());
    }
    fn port(&self, s: C64) -> Vec<C64> {
        vec![
            s.scale(self.equations[0] / self.state[0]),
            C64::ZERO,
            C64::ZERO,
            s.scale(self.equations[1] / self.state[1]),
        ]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.analytic.then_some(Self::tangent)
    }
}

impl ContinuableProblem for Duffing {
    fn set_parameter(&mut self, value: f64) {
        self.omega = value / self.parameter;
    }
}

fn family(state: [f64; 2], equations: [f64; 2], parameter: f64, analytic: bool) -> Duffing {
    Duffing {
        omega: 0.6,
        state,
        equations,
        parameter,
        analytic,
        island_calls: Cell::new(0),
        tangent_calls: Cell::new(0),
    }
}

fn branch(problem: &mut Duffing) -> Vec<ScaledContinuationPoint> {
    let budget = budget(7);
    let mut seed = vec![0.0; 2 * (2 * budget.hb.harmonics + 1)];
    seed[2] = 0.175 * problem.state[0];
    seed[5] = 0.105 * problem.state[1];
    let scales = problem.scales();
    let first = solve_hb_scaled_seeded(
        problem,
        HbAnchor::Forced { omega: 0.6 },
        seed,
        &budget.hb,
        &scales.hb,
    )
    .unwrap();
    let unit = problem.parameter;
    continue_branch_scaled_seeded(
        problem,
        |lambda| lambda / unit,
        0.6 * unit,
        1.8 * unit,
        first.orbit.packed_seed(&HbAnchor::Forced { omega: 0.6 }),
        &budget,
        &scales,
    )
    .unwrap()
}

fn folds(branch: &[ScaledContinuationPoint]) -> Vec<f64> {
    branch
        .windows(3)
        .filter_map(|p| {
            let a = p[1].point.orbit.omega - p[0].point.orbit.omega;
            let b = p[2].point.orbit.omega - p[1].point.orbit.omega;
            (a * b < 0.0).then_some(p[1].point.orbit.omega)
        })
        .collect()
}

// Independent first-harmonic balance F²=A²[(1-w²+3*beta*A²/4)²+c²*w²].
// Its folds approximate the fully coupled seven-harmonic Duffing branch.
fn approximate_folds() -> Vec<f64> {
    let mut found = Vec::new();
    let (mut previous, mut slope) = (f64::NAN, 0.0);
    for i in 200..60_000 {
        let amplitude = i as f64 * 1e-4;
        let a = 1.0 + 0.75 * 0.15 * amplitude * amplitude;
        let b = a - 0.5 * 0.1_f64.powi(2);
        let discriminant = b * b - a * a + (0.3 / amplitude).powi(2);
        if discriminant >= 0.0 {
            let omega = (b + discriminant.sqrt()).sqrt();
            if previous.is_finite() {
                let next = omega - previous;
                if next * slope < 0.0 {
                    found.push(previous);
                }
                slope = next;
            }
            previous = omega;
        }
    }
    found.reverse();
    found
}

#[test]
fn nonlinear_folds_full_waveform_and_parameter_units_are_preserved() {
    let mut base = family([1.0; 2], [1.0; 2], 1.0, true);
    let mut converted = family([1e-6, 1e3], [1e7, 1e-9], 1e4, true);
    let a = branch(&mut base);
    let b = branch(&mut converted);
    let expected = approximate_folds();
    assert_eq!(expected.len(), 2);
    for result in [&a, &b] {
        let measured = folds(result);
        assert_eq!(measured.len(), 2, "{measured:?}");
        for (actual, reference) in measured.iter().zip(&expected) {
            assert!((actual / reference - 1.0).abs() < 0.03);
        }
        assert!(result.last().unwrap().point.orbit.omega >= 1.8);
        for point in result {
            assert!(point.scaled_residual.hypot(point.arclength_residual) < budget(7).hb.tolerance);
            assert!(point.scaled_sampling_defect < budget(7).hb.tolerance);
            assert_eq!(
                point.point.orbit.residual,
                *point.point.orbit.residual_trace.last().unwrap()
            );
        }
    }
    assert_eq!(
        a.len(),
        b.len(),
        "normalized predictor should follow the same arc"
    );
    for (a, b) in a.iter().zip(&b) {
        assert!((a.point.orbit.omega - b.point.orbit.omega).abs() < 2e-6);
        assert!((a.point.lambda - b.point.lambda / converted.parameter).abs() < 2e-6);
        for state in 0..2 {
            for (x, y) in a.point.orbit.coeffs[state]
                .iter()
                .zip(&b.point.orbit.coeffs[state])
            {
                assert!((*x - y.scale(1.0 / converted.state[state])).abs() < 5e-6);
            }
        }
    }
    assert!(
        a.iter()
            .map(|p| p.point.orbit.coeffs[0][3].abs())
            .fold(0.0, f64::max)
            > 1e-3
    );

    // Restart a high-amplitude point at a frequency with multiple responses.
    // Preserve every harmonic and the forced phase, not just its amplitude.
    let selected = a
        .iter()
        .find(|p| {
            p.point.orbit.omega > 1.2
                && p.point.orbit.omega < 1.3
                && p.point.orbit.first_harmonic_amplitude(0) > 2.0
        })
        .unwrap();
    let omega = selected.point.orbit.omega;
    let seed = selected
        .point
        .orbit
        .packed_seed(&HbAnchor::Forced { omega });
    let scales = base.scales();
    let restart = continue_branch_scaled_seeded(
        &mut base,
        |w| w,
        omega,
        omega + 0.03,
        seed.clone(),
        &budget(7),
        &scales,
    )
    .unwrap();
    assert_eq!(
        restart[0]
            .point
            .orbit
            .packed_seed(&HbAnchor::Forced { omega }),
        seed
    );
    assert!(restart[0].point.orbit.coeffs[0][3].abs() > 1e-3);
    assert!(base.tangent_calls.get() > 0 && converted.tangent_calls.get() > 0);
}

#[test]
fn continuation_reuses_exact_tangents_without_coefficient_source_probes() {
    let mut analytic = family([1.0; 2], [1.0; 2], 1.0, true);
    let mut differenced = family([1.0; 2], [1.0; 2], 1.0, false);
    let a = branch(&mut analytic);
    let b = branch(&mut differenced);
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(&b) {
        assert!((a.point.orbit.omega - b.point.orbit.omega).abs() < 2e-6);
        for (a, b) in a
            .point
            .orbit
            .coeffs
            .iter()
            .flatten()
            .zip(b.point.orbit.coeffs.iter().flatten())
        {
            assert!((*a - *b).abs() < 5e-6);
        }
    }
    assert!(analytic.tangent_calls.get() > 0);
    assert_eq!(differenced.tangent_calls.get(), 0);
    assert!(
        4 * analytic.island_calls.get() < differenced.island_calls.get(),
        "analytic {}, differenced {}",
        analytic.island_calls.get(),
        differenced.island_calls.get()
    );
}

struct Septic {
    force: f64,
}
impl OrbitProblem for Septic {
    fn dim(&self) -> usize {
        1
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn port(&self, _: C64) -> Vec<C64> {
        vec![C64::from_re(1e-194)]
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = 1e-200 * (self.force * t.cos() - (x[0] / 1e-6).powi(7));
    }
}
impl ContinuableProblem for Septic {
    fn set_parameter(&mut self, value: f64) {
        self.force = value;
    }
}

#[test]
fn continuation_checks_aliases_in_tiny_physical_equation_units() {
    let scales = ContinuationScales {
        hb: HbScales {
            states: vec![1e-6],
            equations: vec![1e-200],
            angular_frequency: 1.0,
        },
        parameter: 1.0,
    };
    let mut problem = Septic { force: 0.2 };
    let points = continue_branch_scaled_seeded(
        &mut problem,
        |_| 1.0,
        0.2,
        1.0,
        vec![0.0, 1e-7, 0.0],
        &budget(1),
        &scales,
    )
    .unwrap();
    for point in points {
        let z = point.point.orbit.coeffs[0][1].scale(1e6);
        let balance = z + z.scale(35.0 * z.abs().powi(6)) - C64::from_re(0.5 * point.point.lambda);
        assert!(balance.abs() < 2e-9, "{balance:?}");
        assert!(point.point.orbit.aft_samples >= 35);
        assert!(point.scaled_sampling_defect < budget(1).hb.tolerance);
        assert!(point.point.orbit.sampling_defect.is_finite());
    }
}

#[test]
fn malformed_scaled_branch_refuses_before_source_callbacks() {
    struct Counted {
        calls: Cell<usize>,
    }
    impl OrbitProblem for Counted {
        fn dim(&self) -> usize {
            1
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, _: f64, _: &[f64], out: &mut [f64]) {
            self.calls.set(self.calls.get() + 1);
            out[0] = 0.0;
        }
    }
    impl ContinuableProblem for Counted {
        fn set_parameter(&mut self, _: f64) {
            self.calls.set(self.calls.get() + 1);
        }
    }
    let scales = ContinuationScales {
        hb: HbScales {
            states: vec![1.0],
            equations: vec![1.0],
            angular_frequency: 1.0,
        },
        parameter: 1.0,
    };
    for case in 0..9 {
        let mut problem = Counted {
            calls: Cell::new(0),
        };
        let mut scales = scales.clone();
        let mut settings = budget(1);
        let mut seed = vec![0.0; 3];
        let mut end = 2.0;
        match case {
            0 => scales.parameter = 0.0,
            1 => scales.parameter = f64::INFINITY,
            2 => scales.hb.states[0] = f64::NAN,
            3 => scales.hb.equations.clear(),
            4 => {
                seed.pop();
            }
            5 => seed[0] = f64::NAN,
            6 => settings.max_steps = 0,
            7 => settings.min_step = settings.initial_step * 2.0,
            _ => end = 1.0,
        }
        assert!(matches!(
            continue_branch_scaled_seeded(&mut problem, |w| w, 1.0, end, seed, &settings, &scales),
            Err(OrbitError::BadParameter { .. })
        ));
        assert_eq!(problem.calls.get(), 0, "case {case}");
    }
}
