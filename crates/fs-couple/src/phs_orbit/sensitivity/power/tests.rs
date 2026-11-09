//! G1/G3: original power, full Fourier covectors and nonlinear material chains.

use super::*;
use crate::phs_orbit::{
    PeriodicPortInput, solve_forced_harmonic_balance,
    solve_forced_harmonic_balance_analytic_scaled_seeded,
    solve_forced_harmonic_balance_scaled_seeded,
};
use fs_math::c64::C64;
use fs_orbit::sensitivity::{ForcedSensitivityControls, ForcedSensitivityDirection};
use fs_orbit::{HbAnchor, HbBudget, HbOrbit, HbScales};
use fs_phs::{ContinuousStorageTangent, PortHamiltonian, Storage};

fn controls() -> ForcedSensitivityControls {
    ForcedSensitivityControls {
        samples_per_period: 128,
        maximum_unknowns: 256,
        maximum_payload_bytes: 16 * 1024 * 1024,
        maximum_work: 500_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_condition_estimate: 1e9,
        maximum_relative_residual: 1e-12,
        maximum_derivative_residual: 1e-9,
    }
}
fn objectives() -> ForcedMeanObjectiveControls {
    ForcedMeanObjectiveControls {
        objective_scale: 1.0,
        maximum_value_sampling_defect: 1e-10,
        maximum_gradient_sampling_defect: 1e-10,
    }
}
fn scales(n: usize) -> HbScales {
    HbScales {
        states: vec![1.0; n],
        equations: vec![1.0; n],
        angular_frequency: 1.0,
    }
}
fn hb(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 40,
        tolerance: 1e-12,
        deflate_equilibrium: false,
    }
}
fn direction() -> ForcedSensitivityDirection {
    ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 1.0,
        angular_frequency_derivative: 0.0,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tol * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn phased_multiharmonic_supply_and_loss_have_complete_distinct_covectors() {
    let (mass, damping, omega) = (2.0, 1.3, 1.1);
    let system = fs_phs::mass_spring_damper(mass, 10.0, damping).unwrap();
    let drive = vec![
        C64::from_re(0.4),
        C64::new(0.5, -0.2),
        C64::ZERO,
        C64::new(0.08, 0.02),
    ];
    let input = PeriodicPortInput::new(omega, vec![drive.clone()]).unwrap();
    let orbit = solve_forced_harmonic_balance(&system, &input, 0.1, &hb(3)).unwrap();
    let mut prepared =
        PreparedPhsForcedSensitivity::new(&system, &input, &orbit, controls(), &scales(2), || {
            false
        })
        .unwrap();
    let supplied = prepared
        .mean_supplied_power(objectives(), || false)
        .unwrap();
    let loss = prepared
        .mean_resistance_power(objectives(), || false)
        .unwrap();
    let expected = 2.0 * damping / (mass * mass)
        * orbit.coeffs[1][1..]
            .iter()
            .map(|p| p.norm_sq())
            .sum::<f64>();
    close(supplied.value, expected, 3e-13);
    close(loss.value, expected, 3e-13);
    // Equal accepted powers do NOT have equal fixed-waveform covectors.
    let mut supply_gradient = vec![0.0; 14];
    let mut loss_gradient = vec![0.0; 14];
    supply_gradient[1] = drive[0].re / mass;
    loss_gradient[1] = 2.0 * damping * orbit.coeffs[1][0].re / (mass * mass);
    for n in 1..=3 {
        let re = 2 + 4 * (n - 1) + 1;
        let im = re + 2;
        supply_gradient[re] = 2.0 * drive[n].re / mass;
        supply_gradient[im] = 2.0 * drive[n].im / mass;
        loss_gradient[re] = 4.0 * damping * orbit.coeffs[1][n].re / (mass * mass);
        loss_gradient[im] = 4.0 * damping * orbit.coeffs[1][n].im / (mass * mass);
    }
    for (a, b) in supplied.gradient.iter().zip(supply_gradient) {
        close(*a, b, 3e-13);
    }
    for (a, b) in loss.gradient.iter().zip(loss_gradient) {
        close(*a, b, 3e-13);
    }
    for (objective, explicit) in [(&supplied, supplied.value), (&loss, 0.0)] {
        let mut adjoint = prepared
            .prepare_adjoint(&objective.gradient, 1.0, || false)
            .unwrap();
        let derivative = adjoint
            .gradient_port_direction(direction(), explicit, &input, || false)
            .unwrap();
        close(derivative.derivative, 2.0 * expected, 3e-12);
    }
}

struct Memory {
    beta: f64,
    native: bool,
}
impl Storage for Memory {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (2.0 * x[0] * x[0] + x[1] * x[1] + 0.5 * (x[2] - x[0]).powi(2))
            + 0.25 * self.beta * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[
            2.5 * x[0] - 0.5 * x[2] + self.beta * x[0].powi(3),
            x[1],
            0.5 * (x[2] - x[0]),
        ]);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        self.native.then_some(self)
    }
}
impl ContinuousStorageTangent for Memory {
    fn apply(&self, x: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        out.copy_from_slice(&[
            (2.5 + 3.0 * self.beta * x[0] * x[0]) * v[0] - 0.5 * v[2],
            v[1],
            0.5 * (v[2] - v[0]),
        ]);
        true
    }
}
fn memory(beta: f64, native: bool) -> PortHamiltonian {
    PortHamiltonian::new(
        3,
        1,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.15, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.7],
        vec![0.2, 1.0, 0.3],
        Box::new(Memory { beta, native }),
    )
    .unwrap()
}
fn waveform(orbit: &HbOrbit, phase: f64) -> Vec<f64> {
    orbit
        .coeffs
        .iter()
        .map(|row| {
            row[0].re
                + row
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(n, c)| {
                        let (s, z) = (n as f64 * phase).sin_cos();
                        2.0 * (c.re * z - c.im * s)
                    })
                    .sum::<f64>()
        })
        .collect()
}
fn physical_means(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    orbit: &HbOrbit,
) -> [f64; 4] {
    let mut result = [0.0; 4];
    for j in 0..2048 {
        let phase = std::f64::consts::TAU * j as f64 / 2048.0;
        let x = waveform(orbit, phase);
        let mut e = [0.0; 3];
        let mut u = [0.0];
        system.storage().gradient(&x, &mut e);
        assert!(input.sample_into(phase / orbit.omega, &mut u));
        // Independent original diagonal resistance and source map.
        result[0] += (0.2 * e[0] + e[1] + 0.3 * e[2]) * u[0] / 2048.0;
        result[1] += (0.15 * e[0] * e[0] + 0.4 * e[1] * e[1] + 0.7 * e[2] * e[2]) / 2048.0;
        result[2] += 0.2 * x[0].powi(3) * u[0] / 2048.0;
        result[3] += 0.3 * e[0] * x[0].powi(3) / 2048.0;
    }
    result
}

#[test]
fn nonlinear_storage_memory_power_matches_independent_neighboring_physical_cycles() {
    let beta = 0.7;
    let system = memory(beta, true);
    let input = PeriodicPortInput::new(
        1.1,
        vec![vec![
            C64::from_re(0.1),
            C64::new(0.18, -0.06),
            C64::new(0.03, 0.01),
        ]],
    )
    .unwrap();
    let orbit = solve_forced_harmonic_balance(&system, &input, 0.1, &hb(5)).unwrap();
    let mut prepared =
        PreparedPhsForcedSensitivity::new(&system, &input, &orbit, controls(), &scales(3), || {
            false
        })
        .unwrap();
    let supply = prepared
        .mean_supplied_power(objectives(), || false)
        .unwrap();
    let loss = prepared
        .mean_resistance_power(objectives(), || false)
        .unwrap();
    let baseline = physical_means(&system, &input, &orbit);
    close(supply.value, baseline[0], 3e-12);
    close(loss.value, baseline[1], 3e-12);
    let epsilon = 1e-4;
    let mut neighbors = [[0.0; 4]; 2];
    for (i, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let model = memory(beta + sign * epsilon, true);
        let solved = solve_forced_harmonic_balance_scaled_seeded(
            &model,
            &input,
            orbit.packed_seed(&HbAnchor::Forced { omega: orbit.omega }),
            &hb(5),
            &scales(3),
        )
        .unwrap();
        neighbors[i] = physical_means(&model, &input, &solved.orbit);
    }
    for (i, objective) in [&supply, &loss].into_iter().enumerate() {
        let mut adjoint = prepared
            .prepare_adjoint(&objective.gradient, 1.0, || false)
            .unwrap();
        let result = adjoint
            .gradient_storage_direction(
                direction(),
                baseline[2 + i],
                |x, out| {
                    out.copy_from_slice(&[x[0].powi(3), 0.0, 0.0]);
                    true
                },
                || false,
            )
            .unwrap();
        close(
            result.derivative,
            (neighbors[1][i] - neighbors[0][i]) / (2.0 * epsilon),
            3e-8,
        );
    }
}

#[test]
fn source_owned_resistance_gradient_and_explicit_hessian_refusals_remain_reusable() {
    let system = memory(0.7, false);
    let material = Memory {
        beta: 0.7,
        native: false,
    };
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| material.apply(x, v, out);
    let input =
        PeriodicPortInput::new(1.1, vec![vec![C64::from_re(0.1), C64::new(0.18, -0.06)]]).unwrap();
    let solved = solve_forced_harmonic_balance_analytic_scaled_seeded(
        &system,
        &input,
        vec![0.0; 33],
        &hb(5),
        &scales(3),
        &hessian,
    )
    .unwrap();
    let mut prepared = PreparedPhsForcedSensitivity::new_analytic(
        &system,
        &input,
        &solved.orbit,
        controls(),
        &scales(3),
        &hessian,
        || false,
    )
    .unwrap();
    assert!(
        prepared
            .mean_resistance_power(objectives(), || false)
            .is_err()
    );
    assert!(
        prepared
            .mean_supplied_power(objectives(), || false)
            .is_err()
    );
    let baseline = prepared
        .mean_resistance_power_analytic(objectives(), &hessian, || false)
        .unwrap();
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(
        prepared
            .mean_resistance_power_analytic(objectives(), &incomplete, || false)
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        prepared.mean_supplied_power_analytic(objectives(), &hessian, || {
            calls += 1;
            calls > 9
        }),
        Err(PhsForcedSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    assert_eq!(
        baseline,
        prepared
            .mean_resistance_power_analytic(objectives(), &hessian, || false)
            .unwrap()
    );
    let supply = prepared
        .mean_supplied_power_analytic(objectives(), &hessian, || false)
        .unwrap();
    let expected = physical_means(&system, &input, &solved.orbit);
    close(supply.value, expected[0], 3e-12);
    let mut gradient = [0.0; 3];
    system
        .resistance_power_effort_gradient_into(&[2.0, -1.0, 0.3], &mut gradient)
        .unwrap();
    for (a, b) in gradient.into_iter().zip([0.6, -0.8, 0.42]) {
        close(a, b, 2e-15);
    }
    assert!(
        system
            .resistance_power_effort_gradient_into(&[1.0], &mut gradient)
            .is_err()
    );
    let raw = PortHamiltonian::from_raw_parts(
        2,
        0,
        vec![0.0; 4],
        vec![1.0, 4.0, -2.0, 3.0],
        vec![],
        Box::new(fs_phs::QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
    );
    let mut raw_gradient = [0.0; 2];
    raw.resistance_power_effort_gradient_into(&[2.0, -1.0], &mut raw_gradient)
        .unwrap();
    assert_eq!(raw_gradient, [2.0, -2.0]);
    let strong = 2.0_f64.powi(54);
    let raw = PortHamiltonian::from_raw_parts(
        3,
        0,
        vec![0.0; 9],
        vec![strong, 1.0, -strong, 1.0, 0.0, 0.0, -strong, 0.0, 0.0],
        vec![],
        Box::new(Memory {
            beta: 0.0,
            native: true,
        }),
    );
    raw.resistance_power_effort_gradient_into(&[1.0, 1.0, 1.0], &mut gradient)
        .unwrap();
    assert_eq!(gradient[0], 2.0);
}
