//! Exact physical tangents on a nonlinear, memory-coupled branch of cycles.

use fs_math::{c64::C64, det};
use fs_orbit::autonomous::{
    AutonomousBranch, AutonomousContinuationBudget, AutonomousContinuationError, AutonomousScales,
    continue_autonomous_branch,
};
use fs_orbit::{ContinuableProblem, HbOrbit, IslandStateJacobian, OrbitError, OrbitProblem};
use std::cell::Cell;
use std::rc::Rc;

#[derive(Clone, Copy)]
enum Tangent {
    Absent,
    Exact,
    Conditional,
    Refused,
    Incomplete,
}

// A saddle-node of circular cycles with two genuine relaxation states.
// a = x^2+y^2, b = zx^2+zy^2, lambda = (a-1)^2,
// Omega = 1+.2*a+.4*b, and z' = .7 * (x-z).
// The memories feed back into oscillator phase. On the exact cycle,
// b = a*.7^2/(.7^2+Omega^2). A parameter-dependent common mass/time
// factor tau=1+.15*lambda makes physical frequency Omega/tau.
struct MemoryCycle {
    parameter: f64,
    parameter_unit: f64,
    states: [f64; 5],
    equations: [f64; 5],
    tangent: Tangent,
    island_calls: Cell<usize>,
    tangent_calls: Rc<Cell<usize>>,
    tangent_parameter: Cell<f64>,
}

fn lambda(radius: f64) -> f64 {
    (radius * radius - 1.0).powi(2)
}

fn internal_frequency(radius: f64) -> f64 {
    let a = radius * radius;
    let (mut lo, mut hi) = (1.0 + 0.2 * a, 1.0 + 0.6 * a);
    // The scalar balance is strictly increasing; this independent closed-cycle
    // oracle contains no harmonic-balance or continuation implementation.
    for _ in 0..64 {
        let w = 0.5 * (lo + hi);
        if w - 1.0 - 0.2 * a - 0.4 * a * 0.49 / (0.49 + w * w) > 0.0 {
            hi = w;
        } else {
            lo = w;
        }
    }
    0.5 * (lo + hi)
}

impl MemoryCycle {
    fn new(tangent: Tangent, converted: bool) -> Self {
        Self {
            parameter: 0.0,
            parameter_unit: if converted { 1e-4 } else { 1.0 },
            states: if converted {
                [3e4, 1e-9, 2e3, 4e-7, 1e2]
            } else {
                [1.0; 5]
            },
            equations: if converted {
                [7.0, 1e2, 5e-6, 3e-8, 2e4]
            } else {
                [1.0; 5]
            },
            tangent,
            island_calls: Cell::new(0),
            tangent_calls: Rc::new(Cell::new(0)),
            tangent_parameter: Cell::new(f64::NAN),
        }
    }

    fn seed(&self, radius: f64, phase: f64) -> HbOrbit {
        let w = internal_frequency(radius);
        let x = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * radius);
        let y = x * C64::new(0.0, 1.0);
        let memory = C64::from_re(0.7) / C64::new(0.7, w);
        let mut coeffs = vec![vec![C64::ZERO; 4]; 5];
        coeffs[0][0] = C64::from_re(0.3 * self.states[0]);
        for (i, value) in [x, y, memory * x, memory * y].into_iter().enumerate() {
            coeffs[i + 1][1] = value.scale(self.states[i + 1]);
        }
        HbOrbit {
            omega: w / (1.0 + 0.15 * lambda(radius)),
            coeffs,
            residual_trace: Vec::new(),
            residual: f64::NAN,
            aft_samples: 0,
            sampling_defect: f64::NAN,
        }
    }

    fn budget(&self) -> AutonomousContinuationBudget {
        AutonomousContinuationBudget {
            harmonics: 3,
            steps: 24,
            max_newton: 40,
            max_attempts_per_step: 12,
            step_sizes: [0.045, 1e-5, 0.08],
            tolerance: 2e-10,
            minimum_oscillation: 1e-5,
            maximum_unknowns: 64,
            maximum_payload_bytes: 16 << 20,
            maximum_residual_evaluations: 20_000,
            scales: AutonomousScales {
                states: self.states.to_vec(),
                equations: self.equations.to_vec(),
                angular_frequency: 2.0,
                parameter: self.parameter_unit,
            },
        }
    }

    fn exact_tangent(&self, _: f64, physical: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        self.tangent_parameter.set(self.parameter);
        match self.tangent {
            Tangent::Refused => return false,
            Tangent::Conditional if self.parameter / self.parameter_unit >= 0.1 => return false,
            Tangent::Incomplete => {
                out[0] = -self.equations[0] / self.states[0];
                return true;
            }
            _ => {}
        }
        let x: [f64; 5] = std::array::from_fn(|i| physical[i] / self.states[i]);
        let a = x[1] * x[1] + x[2] * x[2];
        let b = x[3] * x[3] + x[4] * x[4];
        let radial = self.parameter / self.parameter_unit - (a - 1.0).powi(2);
        let w = 1.0 + 0.2 * a + 0.4 * b;
        let drdx = -4.0 * (a - 1.0) * x[1];
        let drdy = -4.0 * (a - 1.0) * x[2];
        out.fill(0.0);
        out[0] = -1.0;
        out[6] = 0.4 * x[1] * x[2] + radial + x[1] * drdx;
        out[7] = w + 0.4 * x[2] * x[2] + x[1] * drdy;
        out[11] = -w - 0.4 * x[1] * x[1] + x[2] * drdx;
        out[12] = -0.4 * x[1] * x[2] + radial + x[2] * drdy;
        for column in [3, 4] {
            out[5 + column] = 0.8 * x[column] * x[2];
            out[10 + column] = -0.8 * x[column] * x[1];
        }
        out[16] = 0.7;
        out[18] = -0.7;
        out[22] = 0.7;
        out[24] = -0.7;
        for row in 0..5 {
            for column in 0..5 {
                out[5 * row + column] *= self.equations[row] / self.states[column];
            }
        }
        true
    }
}

impl OrbitProblem for MemoryCycle {
    fn dim(&self) -> usize {
        5
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, physical: &[f64], out: &mut [f64]) {
        self.island_calls.set(self.island_calls.get() + 1);
        let x: [f64; 5] = std::array::from_fn(|i| physical[i] / self.states[i]);
        let a = x[1] * x[1] + x[2] * x[2];
        let b = x[3] * x[3] + x[4] * x[4];
        let radial = self.parameter / self.parameter_unit - (a - 1.0).powi(2);
        let w = 1.0 + 0.2 * a + 0.4 * b;
        out.copy_from_slice(&[
            -(x[0] - 0.3),
            w * x[2] + radial * x[1],
            -w * x[1] + radial * x[2],
            0.7 * (x[1] - x[3]),
            0.7 * (x[2] - x[4]),
        ]);
        for (value, unit) in out.iter_mut().zip(&self.equations) {
            *value *= unit;
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let time_factor = 1.0 + 0.15 * self.parameter / self.parameter_unit;
        let mut matrix = vec![C64::ZERO; 25];
        for i in 0..5 {
            matrix[i * 5 + i] = s.scale(time_factor * self.equations[i] / self.states[i]);
        }
        matrix
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        let available = match self.tangent {
            Tangent::Absent => false,
            Tangent::Conditional => self.parameter / self.parameter_unit < 0.1,
            _ => true,
        };
        available.then_some(Self::exact_tangent)
    }
}

impl ContinuableProblem for MemoryCycle {
    fn set_parameter(&mut self, value: f64) {
        self.parameter = value;
    }
}

fn branch(problem: &mut MemoryCycle) -> AutonomousBranch {
    let seeds = [problem.seed(0.65, 0.43), problem.seed(0.70, 0.43)];
    let unit = problem.parameter_unit;
    let budget = problem.budget();
    continue_autonomous_branch(
        problem,
        [
            (unit * lambda(0.65), &seeds[0]),
            (unit * lambda(0.70), &seeds[1]),
        ],
        &budget,
        || false,
    )
    .unwrap()
}

#[test]
fn exact_tangents_preserve_coupled_memory_cycles_and_units_with_fewer_source_calls() {
    let mut reference: Option<AutonomousBranch> = None;
    for converted in [false, true] {
        let mut exact = MemoryCycle::new(Tangent::Exact, converted);
        let mut probes = MemoryCycle::new(Tangent::Absent, converted);
        let a = branch(&mut exact);
        let b = branch(&mut probes);
        assert_eq!(a.points.len(), 26);
        assert!(a.tangent_preparations > 0);
        assert!(a.tangent_preparations < a.residual_evaluations);
        assert_eq!(b.tangent_preparations, 0);
        assert_eq!(probes.tangent_calls.get(), 0);
        assert!(4 * a.residual_evaluations < b.residual_evaluations);
        assert!(4 * exact.island_calls.get() < probes.island_calls.get());
        assert!(exact.island_calls.get() + exact.tangent_calls.get() < probes.island_calls.get());
        for (p, q) in a.points.iter().zip(&b.points) {
            let r = p.orbit.first_harmonic_amplitude(1) / exact.states[1];
            let w = internal_frequency(r);
            assert!((p.parameter / exact.parameter_unit - lambda(r)).abs() < 3e-8);
            assert!((p.orbit.omega - w / (1.0 + 0.15 * lambda(r))).abs() < 3e-8);
            let memory = C64::from_re(0.7) / C64::new(0.7, w);
            for (source, retained) in [(1, 3), (2, 4)] {
                let expected = memory * p.orbit.coeffs[source][1].scale(1.0 / exact.states[source]);
                let actual = p.orbit.coeffs[retained][1].scale(1.0 / exact.states[retained]);
                assert!((actual - expected).abs() < 3e-8);
            }
            assert!(((p.parameter - q.parameter) / exact.parameter_unit).abs() < 2e-7);
            assert!((p.orbit.omega - q.orbit.omega).abs() < 2e-7);
            for state in 0..5 {
                for order in 0..=3 {
                    assert!(
                        (p.orbit.coeffs[state][order] - q.orbit.coeffs[state][order]).abs()
                            / exact.states[state]
                            < 2e-7
                    );
                }
            }
            assert!(p.scaled_residual <= exact.budget().tolerance);
            assert!(p.scaled_sampling_defect <= exact.budget().tolerance);
            assert!(p.phase_residual.abs() <= exact.budget().tolerance);
            assert!(p.arclength_residual.abs() <= exact.budget().tolerance);
        }
        assert!(a.points.windows(2).any(|p| p[1].parameter > p[0].parameter));
        assert!(a.points.windows(2).any(|p| p[1].parameter < p[0].parameter));
        assert_eq!(exact.parameter, a.points.last().unwrap().parameter);
        if let Some(reference) = &reference {
            for (p, q) in a.points.iter().zip(&reference.points) {
                assert!((p.parameter / exact.parameter_unit - q.parameter).abs() < 2e-7);
                assert!((p.orbit.omega - q.orbit.omega).abs() < 2e-7);
            }
        } else {
            reference = Some(a);
        }
    }
}

#[test]
fn declared_refused_or_incomplete_tangents_restore_the_last_accepted_seed_without_fallback() {
    for tangent in [Tangent::Refused, Tangent::Incomplete] {
        let mut problem = MemoryCycle::new(tangent, true);
        let seeds = [problem.seed(0.65, 0.0), problem.seed(0.70, 0.0)];
        let unit = problem.parameter_unit;
        let budget = problem.budget();
        let result = continue_autonomous_branch(
            &mut problem,
            [
                (unit * lambda(0.65), &seeds[0]),
                (unit * lambda(0.70), &seeds[1]),
            ],
            &budget,
            || false,
        );
        assert!(matches!(
            result,
            Err(AutonomousContinuationError::Numerical(
                OrbitError::TangentRefused { sample: 0 }
            ))
        ));
        assert_eq!(problem.tangent_calls.get(), 1);
        // The physical second seed has a different relaxation phase lag. Its
        // integral phase alignment needs a correction before it is accepted.
        assert!((problem.tangent_parameter.get() / unit - lambda(0.70)).abs() < 1e-14);
        assert_eq!(problem.parameter, unit * lambda(0.65));
    }
}

#[test]
fn cached_workspace_rechecks_the_current_parameter_before_using_an_exact_tangent() {
    let mut exact = MemoryCycle::new(Tangent::Exact, false);
    let mut conditional = MemoryCycle::new(Tangent::Conditional, false);
    let a = branch(&mut exact);
    let b = branch(&mut conditional);
    assert!(b.tangent_preparations > 0);
    assert!(b.tangent_preparations < a.tangent_preparations);
    assert!(b.residual_evaluations > a.residual_evaluations);
    // The declared tangent is absent at both ends and present around the fold.
    // Its callback refuses outside that domain, so retaining a stale callback
    // across a parameter change would fail instead of following this branch.
    assert!(b.points.first().unwrap().parameter > 0.1);
    assert!(b.points.last().unwrap().parameter > 0.1);
    for (p, q) in a.points.iter().zip(&b.points) {
        assert!((p.parameter - q.parameter).abs() < 2e-7);
        assert!((p.orbit.omega - q.orbit.omega).abs() < 2e-7);
        for (p, q) in p.orbit.coeffs.iter().zip(&q.orbit.coeffs) {
            for (&p, &q) in p.iter().zip(q) {
                assert!((p - q).abs() < 2e-7);
            }
        }
    }
}

#[test]
fn cancellation_after_the_analytic_pass_restores_the_last_accepted_seed_without_a_branch() {
    let mut problem = MemoryCycle::new(Tangent::Exact, false);
    let calls = Rc::clone(&problem.tangent_calls);
    let seeds = [problem.seed(0.65, 0.0), problem.seed(0.70, 0.0)];
    let budget = problem.budget();
    let result = continue_autonomous_branch(
        &mut problem,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])],
        &budget,
        || calls.get() > 0,
    );
    assert!(matches!(
        result,
        Err(AutonomousContinuationError::Cancelled)
    ));
    assert!(calls.get() > 0);
    assert!((problem.tangent_parameter.get() - lambda(0.70)).abs() < 1e-14);
    assert_eq!(problem.parameter, lambda(0.65));
}

#[test]
fn analytic_workspace_is_admitted_before_tangent_callbacks_and_residual_limits_still_apply() {
    for tangent in [Tangent::Absent, Tangent::Exact] {
        let mut problem = MemoryCycle::new(tangent, false);
        let seeds = [problem.seed(0.65, 0.0), problem.seed(0.70, 0.0)];
        let mut budget = problem.budget();
        budget.steps = 1;
        // Enough for this problem's existing residual/dense solve path, but
        // not for the additional retained Fourier tangent and phase tables.
        budget.maximum_payload_bytes = 160_000;
        let result = continue_autonomous_branch(
            &mut problem,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])],
            &budget,
            || false,
        );
        if matches!(tangent, Tangent::Exact) {
            assert!(matches!(
                result,
                Err(AutonomousContinuationError::Invalid(_))
            ));
            assert!(problem.island_calls.get() > 0);
            assert_eq!(problem.tangent_calls.get(), 0);
            assert_eq!(problem.parameter, lambda(0.65));
        } else {
            assert_eq!(result.unwrap().points.len(), 3);
        }
    }
    let mut problem = MemoryCycle::new(Tangent::Exact, false);
    let seeds = [problem.seed(0.65, 0.0), problem.seed(0.70, 0.0)];
    let mut budget = problem.budget();
    budget.maximum_residual_evaluations = 5;
    assert!(matches!(
        continue_autonomous_branch(
            &mut problem,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])],
            &budget,
            || false,
        ),
        Err(AutonomousContinuationError::EvaluationLimit { evaluations: 5 })
    ));
    // The cap interrupts the second seed's phase correction after the first
    // seed's two independent residual evaluations accepted that first point.
    assert_eq!(problem.parameter, lambda(0.65));
}
