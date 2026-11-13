//! G1/G3: original active work, retained relaxation harmonics and physical units.
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PhsOrbitError, audit_autonomous_feedback_orbit,
    solve_autonomous_feedback_harmonic_balance_adaptive_scaled_seeded,
    solve_autonomous_feedback_harmonic_balance_analytic_adaptive_scaled_seeded,
    solve_autonomous_feedback_harmonic_balance_analytic_scaled_seeded,
    solve_autonomous_feedback_harmonic_balance_scaled_seeded,
};
use fs_math::c64::C64;
use fs_orbit::{HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, OrbitError};
use fs_phs::{PortHamiltonian, QuadraticStorage, Storage};
use std::{cell::Cell, rc::Rc};

fn budget(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 30,
        tolerance: 1e-11,
        deflate_equilibrium: true,
    }
}
fn scales(states: &[f64], omega: f64) -> HbScales {
    HbScales {
        states: states.to_vec(),
        equations: states.iter().map(|s| omega * s).collect(),
        angular_frequency: omega,
    }
}
fn refinement(states: &[f64], omega: f64) -> HbRefinementBudget {
    HbRefinementBudget {
        max_harmonics: 7,
        residual_tolerances: states.iter().map(|s| 1e-8 * omega * s).collect(),
        state_tolerances: states.iter().map(|s| 1e-8 * s).collect(),
        frequency_tolerance: 1e-8,
    }
}
fn seed(states: &[f64], first: usize, harmonics: usize, omega: f64) -> Vec<f64> {
    let phase = 0.27_f64;
    let mut coeffs = vec![vec![C64::ZERO; harmonics + 1]; states.len()];
    coeffs[first][1] = C64::new(phase.cos(), phase.sin()).scale(0.43 * states[first]);
    coeffs[first + 1][1] = C64::new(phase.sin(), -phase.cos()).scale(0.43 * states[first + 1]);
    HbOrbit {
        omega: 0.93 * omega,
        coeffs,
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
    .packed_seed(&HbAnchor::Autonomous { omega_guess: omega })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:e} != {b:e} within {tolerance:e}"
    );
}

struct MemoryDrive {
    c: f64,
    mu: f64,
    lambda: f64,
    amplitude: f64,
    calls: Cell<usize>,
    refuse: Cell<bool>,
}
impl AutonomousPortFeedback for MemoryDrive {
    fn port_count(&self) -> usize {
        3
    }
    fn input(&self, _: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        self.calls.set(self.calls.get() + 1);
        if self.refuse.get() {
            return false;
        }
        let gain = self.c + self.mu * (1.0 - flow[0] * flow[0] - flow[1] * flow[1]);
        out[0] = gain * flow[0];
        out[1] = gain * flow[1];
        out[2] = self.lambda * self.amplitude * flow[0].powi(3);
        true
    }
}
fn memory_system(states: &[f64; 4], omega: f64, drive: &MemoryDrive) -> PortHamiltonian {
    let mut q = vec![0.0; 16];
    let mut j = vec![0.0; 16];
    let mut r = vec![0.0; 16];
    let mut g = vec![0.0; 12];
    for i in 0..4 {
        q[5 * i] = 1.0 / states[i].powi(2);
        r[5 * i] = [0.7, drive.c, drive.c, drive.lambda][i] * states[i].powi(2);
        if i > 0 {
            g[3 * i + i - 1] = states[i];
        }
    }
    j[6] = -omega * states[1] * states[2];
    j[9] = -j[6];
    PortHamiltonian::new(
        4,
        3,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
}
fn memory_jacobian(drive: &MemoryDrive, states: &[f64; 4], y: &[f64], out: &mut [f64]) -> bool {
    out.fill(0.0);
    let gain = drive.c + drive.mu * (1.0 - y[0] * y[0] - y[1] * y[1]);
    for i in 0..2 {
        for k in 0..2 {
            out[4 * i + k + 1] =
                ((if i == k { gain } else { 0.0 }) - 2.0 * drive.mu * y[i] * y[k]) / states[k + 1];
        }
    }
    out[9] = 3.0 * drive.lambda * drive.amplitude * y[0] * y[0] / states[1];
    true
}
fn drive() -> MemoryDrive {
    MemoryDrive {
        c: 0.3,
        mu: 0.4,
        lambda: 0.9,
        amplitude: 0.6,
        calls: Cell::new(0),
        refuse: Cell::new(false),
    }
}
fn check_memory(orbit: &HbOrbit, states: &[f64; 4], omega: f64, drive: &MemoryDrive) {
    close(orbit.omega, omega, 1e-9);
    close(2.0 * orbit.coeffs[1][1].abs() / states[1], 1.0, 1e-9);
    close(2.0 * orbit.coeffs[2][1].abs() / states[2], 1.0, 1e-9);
    assert!(orbit.coeffs[0].iter().all(|z| z.abs() / states[0] < 1e-12));
    let phase = orbit.coeffs[1][1].im.atan2(orbit.coeffs[1][1].re);
    close(phase, 0.27, 1e-8);
    for (h, peak) in [(1, 0.75 * drive.amplitude), (3, 0.25 * drive.amplitude)] {
        // Original relaxation: m' + lambda*m = lambda*a*cos(theta)^3.
        // Its physical phase lag differs at omega and 3*omega.
        let w = h as f64 * omega;
        let gain = 0.5 * drive.lambda * peak / (drive.lambda.powi(2) + w * w);
        let shift = h as f64 * phase;
        let expected = C64::new(drive.lambda, -w).scale(gain) * C64::new(shift.cos(), shift.sin());
        let actual = orbit.coeffs[3][h].scale(1.0 / states[3]);
        assert!(
            (actual - expected).abs() < 1e-9,
            "h={h}: {actual:?} != {expected:?}"
        );
    }
}

#[test]
fn fixed_and_adaptive_full_state_feedback_keep_relaxation_harmonics_and_port_work_across_units() {
    let omega = 2.4;
    for unit in [1e-30, 1.0, 1e30] {
        let states = [unit, 2.0 * unit, 0.5 * unit, 3.0 * unit];
        let feedback = drive();
        let system = memory_system(&states, omega, &feedback);
        let calls = Cell::new(0);
        let jacobian = |_: &[f64], flow: &[f64], out: &mut [f64]| {
            calls.set(calls.get() + 1);
            memory_jacobian(&feedback, &states, flow, out)
        };
        let fixed = solve_autonomous_feedback_harmonic_balance_scaled_seeded(
            &system,
            &feedback,
            &jacobian,
            seed(&states, 1, 3, omega),
            &budget(3),
            &scales(&states, omega),
        )
        .unwrap();
        assert!(calls.get() > 0 && fixed.scaled_residual < 1e-11);
        check_memory(&fixed.orbit, &states, omega, &feedback);
        let adaptive = solve_autonomous_feedback_harmonic_balance_adaptive_scaled_seeded(
            &system,
            &feedback,
            &jacobian,
            seed(&states, 1, 1, omega),
            &budget(1),
            &refinement(&states, omega),
            &scales(&states, omega),
        )
        .unwrap();
        check_memory(&adaptive.orbit, &states, omega, &feedback);
        assert!(adaptive.refinement_trace[0].omitted_residual_ratio > 1.0);
        let final_check = adaptive.refinement_trace.last().unwrap();
        assert!(final_check.harmonics >= 3 && final_check.total_residual_ratio <= 1.0);
        assert!(final_check.state_change_ratio.unwrap() <= 1.0);
        let memory_energy: f64 = [1, 3]
            .iter()
            .map(|&h| {
                let peak = if h == 1 { 0.75 } else { 0.25 } * feedback.amplitude;
                0.25 * (feedback.lambda * peak).powi(2)
                    / (feedback.lambda.powi(2) + (h as f64 * omega).powi(2))
            })
            .sum();
        let audit =
            audit_autonomous_feedback_orbit(&system, &feedback, &adaptive.orbit, 256).unwrap();
        close(audit.mean_stored_energy_j, 0.5 + memory_energy, 1e-9);
        let loss = feedback.c + 2.0 * feedback.lambda * memory_energy;
        close(audit.mean_resistance_loss_w, loss, 1e-9);
        close(audit.mean_supplied_power_w, loss, 1e-9);
        assert!(audit.relative_power_residual < 1e-9);
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
        for i in 0..2 {
            out[i] = g * x[i];
        }
    }
    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        Some(fs_phs::linear_exact::QuadraticPart {
            matrix: vec![1.0, 0.0, 0.0, 1.0],
            exact: false,
        })
    }
}
struct NonlinearDrive {
    beta: f64,
    c: f64,
    mu: f64,
}
impl AutonomousPortFeedback for NonlinearDrive {
    fn port_count(&self) -> usize {
        2
    }
    fn input(&self, x: &[f64], y: &[f64], out: &mut [f64]) -> bool {
        let r = x[0] * x[0] + x[1] * x[1];
        let gain = self.c + self.mu * (1.0 - r) / (1.0 + self.beta * r);
        for i in 0..2 {
            out[i] = gain * y[i];
        }
        true
    }
}

#[test]
fn authored_full_hessian_and_total_feedback_chain_retain_nonlinear_frequency_pulling() {
    let (omega0, beta, c, mu) = (1.5, 0.4, 0.3, 0.25);
    let omega = omega0 * (1.0 + beta);
    let gradients = Rc::new(Cell::new(0));
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, -omega0, omega0, 0.0],
        vec![c, 0.0, 0.0, c],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(NonlinearStorage {
            beta,
            gradients: gradients.clone(),
        }),
    )
    .unwrap();
    let feedback = NonlinearDrive { beta, c, mu };
    let calls = Cell::new(0);
    let jacobian = |x: &[f64], y: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        let r = x[0] * x[0] + x[1] * x[1];
        let g = 1.0 + beta * r;
        let gain = c + mu * (1.0 - r) / g;
        for i in 0..2 {
            for j in 0..2 {
                let h = (if i == j { g } else { 0.0 }) + 2.0 * beta * x[i] * x[j];
                out[2 * i + j] = gain * h - y[i] * 2.0 * mu * (1.0 + beta) * x[j] / (g * g);
            }
        }
        true
    };
    let physical_scales = scales(&[1.0, 1.0], omega);
    assert!(matches!(
        solve_autonomous_feedback_harmonic_balance_scaled_seeded(
            &system,
            &feedback,
            &jacobian,
            seed(&[1.0, 1.0], 0, 1, omega),
            &budget(1),
            &physical_scales,
        ),
        Err(PhsOrbitError::Invalid { .. })
    ));
    assert_eq!(gradients.get(), 0);
    assert_eq!(calls.get(), 0);
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
    let fixed = solve_autonomous_feedback_harmonic_balance_analytic_scaled_seeded(
        &system,
        &feedback,
        &jacobian,
        seed(&[1.0, 1.0], 0, 1, omega),
        &budget(1),
        &physical_scales,
        &hessian,
    )
    .unwrap();
    close(fixed.orbit.omega, omega, 1e-9);
    assert!(calls.get() > 0);
    assert_eq!(hessian_calls.get(), 2 * calls.get());
    let adaptive = solve_autonomous_feedback_harmonic_balance_analytic_adaptive_scaled_seeded(
        &system,
        &feedback,
        &jacobian,
        seed(&[1.0, 1.0], 0, 1, omega),
        &budget(1),
        &refinement(&[1.0, 1.0], omega),
        &physical_scales,
        &hessian,
    )
    .unwrap();
    close(adaptive.orbit.omega, omega, 1e-9);
    close(2.0 * adaptive.orbit.coeffs[0][1].abs(), 1.0, 1e-9);
    let audit = audit_autonomous_feedback_orbit(&system, &feedback, &adaptive.orbit, 128).unwrap();
    close(audit.mean_stored_energy_j, 0.5 + 0.25 * beta, 1e-9);
    close(audit.mean_supplied_power_w, c * (1.0 + beta).powi(2), 1e-9);
    close(audit.mean_resistance_loss_w, c * (1.0 + beta).powi(2), 1e-9);
    let incomplete_hessian = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 1.0;
        true
    };
    assert!(matches!(
        solve_autonomous_feedback_harmonic_balance_analytic_scaled_seeded(
            &system,
            &feedback,
            &jacobian,
            seed(&[1.0, 1.0], 0, 1, omega),
            &budget(1),
            &physical_scales,
            &incomplete_hessian,
        ),
        Err(PhsOrbitError::Orbit(OrbitError::TangentRefused { .. }))
    ));
}

#[test]
fn malformed_input_and_source_or_incomplete_tangent_refuse_without_silent_fallback() {
    let states = [1.0, 2.0, 0.5, 3.0];
    let omega = 2.4;
    let feedback = drive();
    let system = memory_system(&states, omega, &feedback);
    let jacobian =
        |_: &[f64], y: &[f64], out: &mut [f64]| memory_jacobian(&feedback, &states, y, out);
    let chart = scales(&states, omega);
    assert!(
        solve_autonomous_feedback_harmonic_balance_scaled_seeded(
            &system,
            &feedback,
            &jacobian,
            vec![0.0],
            &budget(3),
            &chart,
        )
        .is_err()
    );
    assert_eq!(feedback.calls.get(), 0);
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        solve_autonomous_feedback_harmonic_balance_scaled_seeded(
            &system,
            &feedback,
            &incomplete,
            seed(&states, 1, 3, omega),
            &budget(3),
            &chart,
        ),
        Err(PhsOrbitError::Orbit(OrbitError::TangentRefused { .. }))
    ));
    feedback.refuse.set(true);
    assert!(matches!(
        solve_autonomous_feedback_harmonic_balance_scaled_seeded(
            &system,
            &feedback,
            &jacobian,
            seed(&states, 1, 3, omega),
            &budget(3),
            &chart,
        ),
        Err(PhsOrbitError::SourceRefused)
    ));
    feedback.refuse.set(false);
    let accepted = solve_autonomous_feedback_harmonic_balance_scaled_seeded(
        &system,
        &feedback,
        &jacobian,
        seed(&states, 1, 3, omega),
        &budget(3),
        &chart,
    )
    .unwrap();
    check_memory(&accepted.orbit, &states, omega, &feedback);
}
