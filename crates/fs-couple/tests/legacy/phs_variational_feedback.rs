//! Original active port work and exact autonomous variational response.
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PhsOrbitError, audit_autonomous_feedback_orbit,
    solve_autonomous_feedback_shooting_variational,
    solve_autonomous_feedback_shooting_variational_analytic,
};
use fs_math::c64::C64;
use fs_orbit::{HbOrbit, ShootBudget, VariationalShootControls, VariationalShootError};
use fs_phs::{PortHamiltonian, QuadraticStorage, Storage};
use std::{cell::Cell, rc::Rc};

fn budget(steps: usize) -> ShootBudget {
    ShootBudget {
        steps_per_period: steps,
        max_newton: 20,
        tolerance: 1e-11,
    }
}
fn controls(scales: &[f64], period: f64) -> VariationalShootControls {
    VariationalShootControls {
        state_scales: scales.to_vec(),
        period_scale: period,
        minimum_autonomous_excursion: 0.01,
        maximum_dimension: 8,
        maximum_payload_bytes: 2_000_000,
        maximum_work: 10_000_000_000,
        maximum_period_maps: 64,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} within {tolerance:e}"
    );
}
fn spectrum(actual: &[(f64, f64)], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for &value in expected {
        let error = actual
            .iter()
            .map(|&(re, im)| (C64::new(re, im) - C64::from_re(value)).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(error < 2e-7, "missing multiplier {value:e} in {actual:?}");
    }
}
fn waveform(omega: f64, scales: &[f64], first: usize, phase: f64) -> HbOrbit {
    let mut coeffs = vec![vec![C64::ZERO; 2]; scales.len()];
    coeffs[first][1] = C64::new(phase.cos(), phase.sin()).scale(scales[first] / 2.0);
    coeffs[first + 1][1] = C64::new(phase.sin(), -phase.cos()).scale(scales[first + 1] / 2.0);
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

struct RadialFeedback {
    c: f64,
    mu: f64,
    calls: Cell<usize>,
}
impl AutonomousPortFeedback for RadialFeedback {
    fn port_count(&self) -> usize {
        2
    }
    fn input(&self, _: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        self.calls.set(self.calls.get() + 1);
        let gain = self.c + self.mu * (1.0 - flows[0] * flows[0] - flows[1] * flows[1]);
        for p in 0..2 {
            out[p] = gain * flows[p];
        }
        true
    }
}
fn radial(scales: &[f64; 3], omega: f64, c: f64, decay: f64) -> PortHamiltonian {
    let mut q = vec![0.0; 9];
    let mut j = vec![0.0; 9];
    let mut r = vec![0.0; 9];
    for i in 0..3 {
        q[i * 3 + i] = 1.0 / (scales[i] * scales[i]);
        r[i * 3 + i] = if i == 0 { decay } else { c } * scales[i] * scales[i];
    }
    j[5] = -omega * scales[1] * scales[2];
    j[7] = -j[5];
    PortHamiltonian::new(
        3,
        2,
        j,
        r,
        vec![0.0, 0.0, scales[1], 0.0, 0.0, scales[2]],
        Box::new(QuadraticStorage::new(q, 3).unwrap()),
    )
    .unwrap()
}
fn radial_jacobian(
    law: &RadialFeedback,
    scales: &[f64; 3],
    flows: &[f64],
    out: &mut [f64],
) -> bool {
    out.fill(0.0);
    let gain = law.c + law.mu * (1.0 - flows[0] * flows[0] - flows[1] * flows[1]);
    for p in 0..2 {
        for j in 0..2 {
            // This TOTAL derivative already contains d(G^T gradH)/dx=1/scale.
            out[p * 3 + j + 1] = ((if p == j { gain } else { 0.0 })
                - 2.0 * law.mu * flows[p] * flows[j])
                / scales[j + 1];
        }
    }
    true
}

#[test]
fn full_state_phase_static_first_coordinate_and_physical_units_keep_active_cycle_and_work() {
    let (omega, c, mu, decay) = (2.4_f64, 0.3, 0.4, 0.7);
    let period = core::f64::consts::TAU / omega;
    let mut reference: Option<Vec<f64>> = None;
    for unit in [1e-30, 1.0, 1e30] {
        let scales = [unit, 2.0 * unit, 0.5 * unit];
        let system = radial(&scales, omega, c, decay);
        let feedback = RadialFeedback {
            c,
            mu,
            calls: Cell::new(0),
        };
        let calls = Cell::new(0);
        let jacobian = |_: &[f64], flows: &[f64], out: &mut [f64]| {
            calls.set(calls.get() + 1);
            radial_jacobian(&feedback, &scales, flows, out)
        };
        // The first coordinate is identically static. Phase must use the full
        // original state; replacing a coordinate-zero closure cannot work here.
        let seed = [0.0, 0.9 * scales[1], 0.1 * scales[2]];
        let result = solve_autonomous_feedback_shooting_variational(
            &system,
            &feedback,
            &jacobian,
            &seed,
            0.95 * period,
            &budget(512),
            &controls(&scales, period),
            || false,
        )
        .unwrap();
        assert_eq!(feedback.calls.get(), result.rhs_evaluations);
        assert_eq!(calls.get(), result.jacobian_evaluations);
        assert_eq!(result.orbit.x0[0], 0.0);
        assert_eq!(result.monodromy.len(), 9);
        assert_eq!(result.period_derivative.as_ref().unwrap().len(), 3);
        let x = &result.orbit.x0;
        close((x[1] / scales[1]).hypot(x[2] / scales[2]), 1.0, 2e-8);
        close(result.orbit.period, period, 2e-8);
        spectrum(
            &result.orbit.multipliers,
            &[1.0, (-2.0 * mu * period).exp(), (-decay * period).exp()],
        );
        assert!(result.normalized_residual < 1e-11 && result.phase_residual < 1e-11);
        if let Some(previous) = &reference {
            for (&a, &b) in result.monodromy.iter().zip(previous) {
                close(a, b, 3e-9);
            }
        } else {
            reference = Some(result.monodromy.clone());
        }
        let phase = (x[2] / scales[2]).atan2(x[1] / scales[1]);
        // This independent analytic continuous waveform is an original H/R/G
        // work oracle, not quadrature of a second numerically integrated orbit.
        let exact = waveform(omega, &scales, 1, phase);
        let audit = audit_autonomous_feedback_orbit(&system, &feedback, &exact, 128).unwrap();
        close(audit.mean_stored_energy_j, 0.5, 1e-13);
        close(audit.mean_resistance_loss_w, c, 1e-13);
        close(audit.mean_supplied_power_w, c, 1e-13);
        assert!(audit.relative_power_residual < 1e-13);
    }
}

struct NonlinearStorage {
    beta: f64,
    gradients: Rc<Cell<usize>>,
}
impl Storage for NonlinearStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let r = x[0] * x[0] + x[1] * x[1];
        0.5 * r + 0.25 * self.beta * r * r
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.gradients.set(self.gradients.get() + 1);
        let g = 1.0 + self.beta * (x[0] * x[0] + x[1] * x[1]);
        out[0] = g * x[0];
        out[1] = g * x[1];
    }
    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        Some(fs_phs::linear_exact::QuadraticPart {
            matrix: vec![1.0, 0.0, 0.0, 1.0],
            exact: false,
        })
    }
}
struct NonlinearFeedback {
    beta: f64,
    c: f64,
    mu: f64,
}
impl AutonomousPortFeedback for NonlinearFeedback {
    fn port_count(&self) -> usize {
        2
    }
    fn input(&self, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        let r = x[0] * x[0] + x[1] * x[1];
        let gain = self.c + self.mu * (1.0 - r) / (1.0 + self.beta * r);
        out[0] = gain * flows[0];
        out[1] = gain * flows[1];
        true
    }
}

#[test]
fn nonlinear_storage_and_total_feedback_chain_keep_frequency_pulling_and_separate_heat() {
    let (omega, beta, c, mu) = (1.5_f64, 0.4, 0.3, 0.25);
    let period = core::f64::consts::TAU / (omega * (1.0 + beta));
    let gradients = Rc::new(Cell::new(0));
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, -omega, omega, 0.0],
        vec![c, 0.0, 0.0, c],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(NonlinearStorage {
            beta,
            gradients: gradients.clone(),
        }),
    )
    .unwrap();
    let feedback = NonlinearFeedback { beta, c, mu };
    let tangent_calls = Cell::new(0);
    let jacobian = |x: &[f64], flow: &[f64], out: &mut [f64]| {
        tangent_calls.set(tangent_calls.get() + 1);
        let r = x[0] * x[0] + x[1] * x[1];
        let g = 1.0 + beta * r;
        let gain = c + mu * (1.0 - r) / g;
        for i in 0..2 {
            for j in 0..2 {
                let hessian = (if i == j { g } else { 0.0 }) + 2.0 * beta * x[i] * x[j];
                let gain_derivative = -2.0 * mu * (1.0 + beta) * x[j] / (g * g);
                out[i * 2 + j] = gain * hessian + flow[i] * gain_derivative;
            }
        }
        true
    };
    let seed = [0.9, 0.1];
    let limits = controls(&[1.0, 1.0], period);
    assert!(matches!(
        solve_autonomous_feedback_shooting_variational(
            &system,
            &feedback,
            &jacobian,
            &seed,
            period,
            &budget(512),
            &limits,
            || false,
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::TangentUnavailable
        ))
    ));
    assert_eq!(gradients.get(), 0);
    assert_eq!(tangent_calls.get(), 0);
    let hessian_calls = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        hessian_calls.set(hessian_calls.get() + 1);
        let g = 1.0 + beta * (x[0] * x[0] + x[1] * x[1]);
        let dot = x[0] * v[0] + x[1] * v[1];
        for i in 0..2 {
            out[i] = g * v[i] + 2.0 * beta * x[i] * dot;
        }
        true
    };
    let result = solve_autonomous_feedback_shooting_variational_analytic(
        &system,
        &feedback,
        &jacobian,
        &seed,
        0.97 * period,
        &budget(512),
        &limits,
        &hessian,
        || false,
    )
    .unwrap();
    assert_eq!(hessian_calls.get(), 2 * result.jacobian_evaluations);
    assert_eq!(tangent_calls.get(), result.jacobian_evaluations);
    assert_eq!(
        gradients.get(),
        result.rhs_evaluations + result.jacobian_evaluations
    );
    let x = &result.orbit.x0;
    close(x[0].hypot(x[1]), 1.0, 3e-8);
    close(result.orbit.period, period, 3e-8);
    spectrum(
        &result.orbit.multipliers,
        &[1.0, (-2.0 * mu * period).exp()],
    );
    let exact = waveform(omega * (1.0 + beta), &[1.0, 1.0], 0, x[1].atan2(x[0]));
    let audit = audit_autonomous_feedback_orbit(&system, &feedback, &exact, 128).unwrap();
    close(audit.mean_stored_energy_j, 0.5 + 0.25 * beta, 1e-13);
    close(
        audit.mean_resistance_loss_w,
        c * (1.0 + beta).powi(2),
        1e-13,
    );
    close(audit.mean_supplied_power_w, c * (1.0 + beta).powi(2), 1e-13);
    assert!(audit.relative_power_residual < 1e-13);
    let partial_hessian = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        solve_autonomous_feedback_shooting_variational_analytic(
            &system,
            &feedback,
            &jacobian,
            &seed,
            period,
            &budget(64),
            &limits,
            &partial_hessian,
            || false,
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::SourceRefused { tangent: true, .. }
        ))
    ));
}

#[test]
fn complete_feedback_tangent_admission_budgets_cancellation_and_replay_are_transactional() {
    let scales = [1.0, 2.0, 0.5];
    let period = core::f64::consts::PI;
    let system = radial(&scales, 2.0, 0.3, 0.7);
    let feedback = RadialFeedback {
        c: 0.3,
        mu: 0.4,
        calls: Cell::new(0),
    };
    let jacobian =
        |_: &[f64], flow: &[f64], out: &mut [f64]| radial_jacobian(&feedback, &scales, flow, out);
    let seed = [0.0, 1.8, 0.05];
    let limits = controls(&scales, period);
    let grid = budget(128);
    let mut bad_limits = limits.clone();
    bad_limits.state_scales.clear();
    assert!(
        solve_autonomous_feedback_shooting_variational(
            &system,
            &feedback,
            &jacobian,
            &seed,
            period,
            &grid,
            &bad_limits,
            || false
        )
        .is_err()
    );
    for bad_scale in [0.0, f64::NAN] {
        let mut bad_limits = limits.clone();
        bad_limits.state_scales[1] = bad_scale;
        assert!(
            solve_autonomous_feedback_shooting_variational(
                &system,
                &feedback,
                &jacobian,
                &seed,
                period,
                &grid,
                &bad_limits,
                || false
            )
            .is_err()
        );
    }
    for memory in [false, true] {
        let mut bad_limits = limits.clone();
        if memory {
            bad_limits.maximum_payload_bytes = 1;
        } else {
            bad_limits.maximum_work = 1;
        }
        assert!(matches!(
            solve_autonomous_feedback_shooting_variational(
                &system,
                &feedback,
                &jacobian,
                &seed,
                period,
                &grid,
                &bad_limits,
                || false
            ),
            Err(PhsOrbitError::Variational(
                VariationalShootError::Budget { .. }
            ))
        ));
    }
    assert_eq!(feedback.calls.get(), 0);
    assert!(matches!(
        solve_autonomous_feedback_shooting_variational(
            &system,
            &feedback,
            &jacobian,
            &seed,
            period,
            &grid,
            &limits,
            || true
        ),
        Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
    ));
    assert_eq!(feedback.calls.get(), 0);
    for kind in 0..3 {
        let bad_tangent = move |_: &[f64], _: &[f64], out: &mut [f64]| {
            if kind == 0 {
                return false;
            }
            out[0] = 0.0;
            if kind == 2 {
                out.fill(0.0);
                out[1] = f64::INFINITY;
            }
            true
        };
        assert!(matches!(
            solve_autonomous_feedback_shooting_variational(
                &system,
                &feedback,
                &bad_tangent,
                &seed,
                period,
                &grid,
                &limits,
                || false
            ),
            Err(PhsOrbitError::Variational(
                VariationalShootError::SourceRefused { tangent: true, .. }
            ))
        ));
    }
    feedback.calls.set(0);
    assert!(matches!(
        solve_autonomous_feedback_shooting_variational(
            &system,
            &feedback,
            &jacobian,
            &seed,
            period,
            &grid,
            &limits,
            || feedback.calls.get() >= 20
        ),
        Err(PhsOrbitError::Variational(VariationalShootError::Cancelled))
    ));
    let first = solve_autonomous_feedback_shooting_variational(
        &system,
        &feedback,
        &jacobian,
        &seed,
        period,
        &grid,
        &limits,
        || false,
    )
    .unwrap();
    let replay = solve_autonomous_feedback_shooting_variational(
        &system,
        &feedback,
        &jacobian,
        &seed,
        period,
        &grid,
        &limits,
        || false,
    )
    .unwrap();
    assert_eq!(first.orbit.x0, replay.orbit.x0);
    assert_eq!(first.monodromy, replay.monodromy);
    assert_eq!(first.orbit.period.to_bits(), replay.orbit.period.to_bits());
    assert!(
        solve_autonomous_feedback_shooting_variational(
            &system,
            &feedback,
            &jacobian,
            &[0.0; 3],
            period,
            &grid,
            &limits,
            || false
        )
        .is_err()
    );
    struct Refused;
    impl AutonomousPortFeedback for Refused {
        fn port_count(&self) -> usize {
            2
        }
        fn input(&self, _: &[f64], _: &[f64], _: &mut [f64]) -> bool {
            false
        }
    }
    assert!(matches!(
        solve_autonomous_feedback_shooting_variational(
            &system,
            &Refused,
            &jacobian,
            &seed,
            period,
            &grid,
            &limits,
            || false
        ),
        Err(PhsOrbitError::SourceRefused)
    ));
}
