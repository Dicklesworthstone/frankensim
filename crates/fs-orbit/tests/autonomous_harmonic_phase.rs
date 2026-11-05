//! Full-waveform phase, physical conditioning and original autonomous balance.

use fs_math::{c64::C64, det};
use fs_orbit::{
    HbBudget, HbRefinementBudget, HbScales, IslandStateJacobian, OrbitError, OrbitProblem,
    solve_hb_autonomous_adaptive_scaled_seeded, solve_hb_autonomous_scaled_seeded,
};
use std::cell::Cell;

const OMEGA: f64 = 2.3;
const GROWTH: f64 = 0.7;
const BIAS: f64 = 0.4;

// A stable static coordinate precedes an isolated Stuart-Landau cycle. Both
// state and equation units are independent, so the linear port is not just sI.
struct Cycle {
    states: [f64; 3],
    equations: [f64; 3],
    exact: bool,
    autonomous: bool,
    calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}

impl Cycle {
    fn new(exact: bool) -> Self {
        Self {
            states: [1.0e24, 1.0e-18, 1.0e12],
            equations: [1.0e-10, 1.0e12, 1.0e-18],
            exact,
            autonomous: true,
            calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }

    fn scales(&self) -> HbScales {
        HbScales {
            states: self.states.to_vec(),
            equations: self.equations.to_vec(),
            angular_frequency: OMEGA,
        }
    }

    fn seed(&self, n: usize, phase: f64) -> Vec<f64> {
        let mut seed = vec![0.0; 3 * (2 * n + 1) + 1];
        seed[0] = 0.2 * self.states[0];
        let z = C64::new(det::cos(phase), det::sin(phase)).scale(0.6);
        let v = C64::new(z.im, -z.re);
        let re = 3 * (2 * n - 1);
        seed[re + 1] = self.states[1] * z.re;
        seed[re + 2] = self.states[2] * v.re;
        seed[re + 4] = self.states[1] * z.im;
        seed[re + 5] = self.states[2] * v.im;
        *seed.last_mut().unwrap() = 1.15 * OMEGA / n as f64;
        seed
    }

    fn derivative(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        let u = x[1] / self.states[1];
        let v = x[2] / self.states[2];
        let physical = [
            [-1.0, 0.0, 0.0],
            [
                0.0,
                GROWTH * (1.0 - 3.0 * u * u - v * v),
                -OMEGA - 2.0 * GROWTH * u * v,
            ],
            [
                0.0,
                OMEGA - 2.0 * GROWTH * u * v,
                GROWTH * (1.0 - u * u - 3.0 * v * v),
            ],
        ];
        for i in 0..3 {
            for j in 0..3 {
                out[3 * i + j] = (self.equations[i] / self.states[j]) * physical[i][j];
            }
        }
        true
    }
}

impl OrbitProblem for Cycle {
    fn dim(&self) -> usize {
        3
    }

    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        let u = x[1] / self.states[1];
        let v = x[2] / self.states[2];
        let growth = GROWTH * (1.0 - u * u - v * v);
        out[0] = self.equations[0] * (BIAS - x[0] / self.states[0]);
        out[1] = self.equations[1] * (growth * u - OMEGA * v);
        out[2] = self.equations[2] * (OMEGA * u + growth * v);
    }

    fn port(&self, s: C64) -> Vec<C64> {
        let mut port = vec![C64::ZERO; 9];
        for k in 0..3 {
            port[4 * k] = s.scale(self.equations[k] / self.states[k]);
        }
        port
    }

    fn autonomous(&self) -> bool {
        self.autonomous
    }

    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.exact.then_some(Self::derivative)
    }
}

fn budget(harmonics: usize, deflate_equilibrium: bool) -> HbBudget {
    HbBudget {
        harmonics,
        tolerance: 1.0e-11,
        deflate_equilibrium,
        ..HbBudget::default()
    }
}

fn check_cycle(problem: &Cycle, orbit: &fs_orbit::HbOrbit, n: usize, phase: f64) {
    assert!((orbit.omega * n as f64 - OMEGA).abs() < 1.0e-9);
    assert!((orbit.coeffs[0][0].re / problem.states[0] - BIAS).abs() < 1.0e-10);
    let z = C64::new(det::cos(phase), det::sin(phase)).scale(0.5);
    let expected = [C64::ZERO, z, C64::new(z.im, -z.re)];
    for k in 0..3 {
        for h in 1..orbit.coeffs[k].len() {
            let target = if h == n { expected[k] } else { C64::ZERO };
            assert!((orbit.coeffs[k][h].scale(1.0 / problem.states[k]) - target).abs() < 1.0e-9);
        }
    }
}

#[test]
fn static_first_state_keeps_all_balances_phase_and_independent_physical_units() {
    for deflate in [false, true] {
        for exact in [false, true] {
            let problem = Cycle::new(exact);
            let phase = 0.71;
            let result = solve_hb_autonomous_scaled_seeded(
                &problem,
                problem.seed(1, phase),
                &budget(1, deflate),
                &problem.scales(),
            )
            .unwrap();
            check_cycle(&problem, &result.orbit, 1, phase);
            assert!(result.scaled_residual < 1.0e-11);
            assert!(result.scaled_sampling_defect < 1.0e-11);
            assert_eq!(problem.tangent_calls.get() > 0, exact);
        }
    }
}

#[test]
fn a_higher_harmonic_alone_can_define_phase_and_survive_order_refinement() {
    // This is the same cycle represented twice per period; primitive-period
    // selection is deliberately separate from phase admission and balance.
    for phase in [0.23, 1.17] {
        let problem = Cycle::new(true);
        let refinement = HbRefinementBudget {
            max_harmonics: 5,
            residual_tolerances: problem.equations.iter().map(|v| v * 1.0e-8).collect(),
            state_tolerances: problem.states.iter().map(|v| v * 1.0e-8).collect(),
            frequency_tolerance: 1.0e-9,
        };
        let result = solve_hb_autonomous_adaptive_scaled_seeded(
            &problem,
            problem.seed(2, phase),
            &budget(2, true),
            &refinement,
            &problem.scales(),
        )
        .unwrap();
        check_cycle(&problem, &result.orbit, 2, phase);
        assert_eq!(result.refinement_trace.len(), 2);
        assert_eq!(result.refinement_trace[0].harmonics, 2);
        assert_eq!(result.refinement_trace[1].harmonics, 5);
        assert!(result.refinement_trace[1].state_change_ratio.unwrap() < 1.0);
    }
}

#[test]
fn invalid_phase_inputs_refuse_before_original_source_or_tangent_evaluation() {
    let mut problem = Cycle::new(true);
    let hb = budget(1, true);
    let scales = problem.scales();
    let mut dc_only = vec![0.0; 10];
    dc_only[0] = BIAS * problem.states[0];
    dc_only[9] = OMEGA;
    assert_eq!(
        solve_hb_autonomous_scaled_seeded(&problem, dc_only, &hb, &scales).unwrap_err(),
        OrbitError::TrivialCollapse,
    );
    let mut bad_scales = scales.clone();
    bad_scales.states[1] = 0.0;
    let mut negative_frequency = problem.seed(1, 0.0);
    negative_frequency[9] = -1.0;
    for (seed, units) in [
        (problem.seed(1, 0.0), &bad_scales),
        (vec![], &scales),
        (negative_frequency, &scales),
    ] {
        assert!(matches!(
            solve_hb_autonomous_scaled_seeded(&problem, seed, &hb, units),
            Err(OrbitError::BadParameter { .. }),
        ));
    }
    problem.autonomous = false;
    assert!(matches!(
        solve_hb_autonomous_scaled_seeded(&problem, problem.seed(1, 0.0), &hb, &scales),
        Err(OrbitError::BadParameter { .. }),
    ));
    assert_eq!(problem.calls.get(), 0);
    assert_eq!(problem.tangent_calls.get(), 0);
}
