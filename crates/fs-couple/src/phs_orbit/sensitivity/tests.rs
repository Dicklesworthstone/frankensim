use super::*;
use crate::phs_orbit::{
    solve_forced_harmonic_balance, solve_forced_harmonic_balance_scaled_seeded,
};
use fs_math::c64::C64;
use fs_orbit::HbBudget;
use fs_phs::{ContinuousStorageTangent, Storage};

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
        maximum_derivative_residual: 1e-10,
    }
}
fn scales(dimension: usize) -> HbScales {
    HbScales {
        states: vec![1.0; dimension],
        equations: vec![1.0; dimension],
        angular_frequency: 1.0,
    }
}
fn direction() -> ForcedSensitivityDirection {
    ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 1.0,
        angular_frequency_derivative: 0.0,
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

#[test]
fn original_port_and_storage_maps_match_dc_and_phased_multiharmonic_mechanics() {
    let (mass, stiffness, damping, omega) = (2.0, 10.0, 1.3, 1.1);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let input = PeriodicPortInput::new(
        omega,
        vec![vec![
            C64::from_re(0.4),
            C64::new(0.5, -0.2),
            C64::ZERO,
            C64::new(0.08, 0.02),
        ]],
    )
    .unwrap();
    let orbit = solve_forced_harmonic_balance(&system, &input, 0.1, &hb(3)).unwrap();
    let mut prepared =
        PreparedPhsForcedSensitivity::new(&system, &input, &orbit, controls(), &scales(2), || {
            false
        })
        .unwrap();
    let coefficients = vec![
        C64::from_re(0.3),
        C64::new(0.1, 0.2),
        C64::new(0.0, -0.02),
        C64::ZERO,
    ];
    let port = PeriodicPortInput::new(omega, vec![coefficients.clone()]).unwrap();
    let force = prepared
        .solve_port_direction(direction(), &port, || false)
        .unwrap();
    let material = prepared
        .solve_storage_direction(
            direction(),
            |x, out| {
                out.copy_from_slice(&[x[0], 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    for n in 0..=3 {
        let w = n as f64 * omega;
        let impedance = C64::new(stiffness - mass * w * w, damping * w);
        let dq = coefficients[n] / impedance;
        let dk = orbit.coeffs[0][n].scale(-1.0) / impedance;
        for (response, expected) in [(&force, dq), (&material, dk)] {
            assert!((response.coeffs[0][n] - expected).abs() < 2e-12);
            assert!((response.coeffs[1][n] - expected * C64::new(0.0, mass * w)).abs() < 2e-12);
        }
    }
    let mut unresolved = vec![C64::ZERO; 5];
    unresolved[4] = C64::from_re(f64::from_bits(1));
    let unresolved = PeriodicPortInput::new(omega, vec![unresolved]).unwrap();
    assert!(matches!(
        prepared.solve_port_direction(direction(), &unresolved, || false),
        Err(PhsForcedSensitivityError::Source(
            PhsOrbitError::Invalid { .. }
        ))
    ));
    assert!(matches!(
        prepared.solve_storage_direction(
            direction(),
            |_, out| {
                out[0] = 1.0;
                true
            },
            || false
        ),
        Err(PhsForcedSensitivityError::Sensitivity(
            ForcedSensitivityError::IslandDirectionRefused { .. }
        ))
    ));
}

// One mechanical Duffing coordinate and one explicit relaxing strain. The
// material parameter changes the ORIGINAL quartic potential, not a fitted
// amplitude-dependent modal frequency. p=(unit mass)*q_dot.
struct NonlinearMemory {
    beta: f64,
}
impl Storage for NonlinearMemory {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (2.0 * x[0] * x[0] + x[1] * x[1] + 0.5 * (x[2] - x[0]).powi(2))
            + 0.25 * self.beta * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[
            2.0 * x[0] + self.beta * x[0].powi(3) + 0.5 * (x[0] - x[2]),
            x[1],
            0.5 * (x[2] - x[0]),
        ]);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for NonlinearMemory {
    fn apply(&self, x: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        if x.len() != 3 || v.len() != 3 || out.len() != 3 {
            return false;
        }
        out.copy_from_slice(&[
            (2.5 + 3.0 * self.beta * x[0] * x[0]) * v[0] - 0.5 * v[2],
            v[1],
            0.5 * (v[2] - v[0]),
        ]);
        out.iter().all(|v| v.is_finite())
    }
}
fn memory_system(beta: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        3,
        1,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 2.0],
        vec![0.0, 1.0, 0.0],
        Box::new(NonlinearMemory { beta }),
    )
    .unwrap()
}
fn memory_orbit(system: &PortHamiltonian, input: &PeriodicPortInput) -> HbOrbit {
    solve_forced_harmonic_balance_scaled_seeded(
        system,
        input,
        vec![0.0; 3 * 19],
        &hb(9),
        &scales(3),
    )
    .unwrap()
    .orbit
}

#[test]
fn native_storage_tangent_retains_nonlinear_harmonics_and_explicit_memory() {
    let system = memory_system(1.4);
    let input = PeriodicPortInput::new(0.85, vec![vec![C64::ZERO, C64::new(0.7, -0.2)]]).unwrap();
    let orbit = memory_orbit(&system, &input);
    let mut prepared =
        PreparedPhsForcedSensitivity::new(&system, &input, &orbit, controls(), &scales(3), || {
            false
        })
        .unwrap();
    let response = prepared
        .solve_storage_direction(
            direction(),
            |x, out| {
                out.copy_from_slice(&[x[0].powi(3), 0.0, 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    assert!(response.coeffs[0][3].abs() > 1e-4);
    assert!(response.coeffs[2][3].abs() > 1e-5);
    let step = 1e-4;
    let lower = memory_orbit(&memory_system(1.4 - step), &input);
    let upper = memory_orbit(&memory_system(1.4 + step), &input);
    for state in 0..3 {
        for n in 0..=9 {
            let reference = (upper.coeffs[state][n] - lower.coeffs[state][n]).scale(0.5 / step);
            assert!((response.coeffs[state][n] - reference).abs() < 2e-9);
        }
    }
    // Authored override and native storage action differentiate the same source.
    let storage = NonlinearMemory { beta: 1.4 };
    let action = |x: &[f64], v: &[f64], out: &mut [f64]| storage.apply(x, v, out);
    let mut explicit = PreparedPhsForcedSensitivity::new_analytic(
        &system,
        &input,
        &orbit,
        controls(),
        &scales(3),
        &action,
        || false,
    )
    .unwrap();
    let repeated = explicit
        .solve_storage_direction(
            direction(),
            |x, out| {
                out.copy_from_slice(&[x[0].powi(3), 0.0, 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    assert_eq!(response.coeffs, repeated.coeffs);
}
