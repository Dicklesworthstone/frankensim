//! Original forced-response frequency tangents, including nonlinear harmonics.

use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput,
    sideband::{PhsPortSidebandControls, PhsPortSidebandError, PhsPortSidebandSource},
    solve_forced_harmonic_balance_scaled_seeded,
};
use fs_math::{c64::C64, det};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbScales,
    sideband::{SidebandControls, SidebandError, SidebandScales},
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn controls(n: usize, m: usize, h: usize, output: usize) -> PhsPortSidebandControls {
    PhsPortSidebandControls {
        sideband: SidebandControls {
            harmonics: h,
            samples_per_period: 128,
            maximum_unknowns: n * (2 * h + 1),
            maximum_payload_bytes: 32 << 20,
            maximum_work: 1_000_000_000,
            maximum_base_residual: 2e-10,
            maximum_tangent_sampling_defect: 2e-10,
            maximum_relative_residual: 1e-12,
            maximum_omitted_residual: 2e-8,
        },
        output_harmonics: output,
        initial_samples: 64,
        max_samples: 256,
        max_hessian_actions: 1024,
        max_workspace_bytes: 8 << 20,
        max_observation_work: 100_000_000,
        envelope_tolerances: vec![2e-10; m],
    }
}

fn scales(n: usize) -> SidebandScales {
    SidebandScales {
        states: (0..n).map(|j| 0.5 + j as f64).collect(),
        equations: (0..n).map(|j| 0.7 + 0.4 * j as f64).collect(),
    }
}

fn orbit(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs,
        // Prior solver metadata is deliberately untrusted; the source must
        // re-evaluate the complete base equations before its tangent is used.
        residual_trace: vec![0.0],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:?}, expected {expected:?}"
    );
}

#[test]
fn logarithmic_frequency_tangent_keeps_inertia_phase_and_material_memory() {
    let (mass, stiffness, damping, arm, tau) = (2.0, 4.0, 0.3, 3.0, 0.4);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let forcing = [
        C64::from_re(0.02),
        C64::new(0.04, -0.02),
        C64::new(-0.01, 0.03),
        C64::new(0.007, 0.002),
    ];
    for omega in [0.9, 1.4, 2.3] {
        // Independent original dynamic stiffness, retaining the Maxwell arm.
        let exact = |w: f64, f: C64| {
            let s = C64::new(0.0, w);
            let relaxation = C64::ONE + s.scale(tau);
            let dynamic = C64::new(stiffness - mass * w * w, damping * w)
                + (s.scale(tau) / relaxation).scale(arm);
            let log_dynamic = C64::new(-2.0 * mass * w * w, damping * w)
                + (s.scale(tau) / (relaxation * relaxation)).scale(arm);
            let q = f / dynamic;
            let log_q = -(q * log_dynamic / dynamic);
            let state = [q, s * q.scale(mass), q.scale(det::sqrt(arm)) / relaxation];
            let tangent = [
                log_q,
                s * (q + log_q).scale(mass),
                (log_q / relaxation - q * s.scale(tau) / (relaxation * relaxation))
                    .scale(det::sqrt(arm)),
            ];
            (state, tangent)
        };
        let mut coefficients = vec![vec![C64::ZERO; forcing.len()]; 3];
        for (k, &force) in forcing.iter().enumerate() {
            let state = exact(k as f64 * omega, force).0;
            for j in 0..3 {
                coefficients[j][k] = state[j];
            }
        }
        let base = orbit(omega, coefficients);
        let input = PeriodicPortInput::new(omega, vec![forcing.to_vec()]).unwrap();
        let mut source =
            PhsPortSidebandSource::forced(&system, &input, &base, &controls(3, 1, 3, 3)).unwrap();
        let mut prepared = source.prepare(&scales(3), || false).unwrap();
        let result = prepared.forced_frequency_tangent(|| false).unwrap();
        assert_eq!(result.state_dimension, 3);
        for (k, &force) in forcing.iter().enumerate() {
            let expected = exact(k as f64 * omega, force).1;
            for j in 0..3 {
                near(
                    result.state_log_frequency_derivative.coefficients[(3 + k) * 3 + j],
                    expected[j],
                    8e-12,
                );
                near(
                    result.state_log_frequency_derivative.coefficients[(3 - k) * 3 + j],
                    expected[j].conj(),
                    8e-12,
                );
            }
            near(
                result.port_log_frequency_derivative[3 + k],
                expected[1].scale(1.0 / mass),
                8e-12,
            );
        }
        assert_eq!(result.observation_trace.len(), 3);
        assert!(result.state_log_frequency_derivative.relative_residual < 1e-12);
        assert!(result.state_log_frequency_derivative.omitted_residual < 1e-12);
        if omega == 1.4 {
            drop(prepared);
            drop(source);
            let mut truncated =
                PhsPortSidebandSource::forced(&system, &input, &base, &controls(3, 1, 2, 3))
                    .unwrap();
            let mut prepared = truncated.prepare(&scales(3), || false).unwrap();
            assert!(matches!(
                prepared.forced_frequency_tangent(|| false),
                Err(PhsPortSidebandError::Invalid(_))
            ));
        }
    }
}

fn nonlinear_fixture() -> (PortHamiltonian, PeriodicPortInput, HbOrbit) {
    let (stiffness, cubic, damping, omega, amplitude, phase) = (4.0, 2.0, 0.3, 2.05, 0.3, 0.37);
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(SosModalStorage {
            omegas: vec![det::sqrt(stiffness)],
            channels: vec![StressChannel {
                coefficient: cubic,
                coupling: vec![1.0],
            }],
        }),
    )
    .unwrap();
    let q = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let mut states = vec![vec![C64::ZERO; 10]; 2];
    states[0][1] = q;
    states[1][1] = C64::new(0.0, omega) * q;
    // The original cubic law makes this an exact sinusoidal base under its
    // fundamental AND third harmonic. These phasors then stay fixed as omega
    // changes; the perturbed orbit is not forced to remain sinusoidal.
    let mut drive = vec![vec![C64::ZERO; 10]; 2];
    drive[1][1] = q * C64::new(
        stiffness - omega * omega + 0.75 * cubic * amplitude * amplitude,
        damping * omega,
    );
    drive[1][3] = (q * q * q).scale(cubic);
    (
        system,
        PeriodicPortInput::new(omega, drive).unwrap(),
        orbit(omega, states),
    )
}

#[test]
fn finite_amplitude_frequency_tangent_agrees_with_original_nonlinear_orbits() {
    let (system, input, base) = nonlinear_fixture();
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 9, 11)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let tangent = prepared.forced_frequency_tangent(|| false).unwrap();
    let state = |k: isize, j: usize| {
        tangent.state_log_frequency_derivative.coefficients[(k + 9) as usize * 2 + j]
    };
    // Changing frequency generates additional harmonics, even though the exact
    // authored base contains only a fundamental. A scalar fitted resonance or
    // replacing Hxx(theta) by its mean would discard this response.
    assert!(state(3, 0).abs() > 1e-4);
    let delta = 1e-5;
    let budget = HbBudget {
        harmonics: 9,
        max_newton: 40,
        tolerance: 2e-12,
        deflate_equilibrium: false,
    };
    let physical = HbScales {
        states: vec![0.3, 0.6],
        equations: vec![0.6, 1.2],
        angular_frequency: base.omega,
    };
    let solve = |sign: f64| {
        let omega = base.omega * det::exp(sign * delta);
        solve_forced_harmonic_balance_scaled_seeded(
            &system,
            &input.at_frequency(omega).unwrap(),
            base.packed_seed(&HbAnchor::Forced { omega }),
            &budget,
            &physical,
        )
        .unwrap()
        .orbit
    };
    let (minus, plus) = (solve(-1.0), solve(1.0));
    for j in 0..2 {
        for k in 0..=9 {
            let difference = (plus.coeffs[j][k] - minus.coeffs[j][k]).scale(0.5 / delta);
            near(state(k as isize, j), difference, 8e-7);
            near(state(-(k as isize), j), difference.conj(), 8e-7);
        }
    }
    // Check the nonlinear ORIGINAL physical outputs independently at common
    // forcing phases. State coefficients are not the restoring-effort output.
    for sample in 0..31 {
        let theta = 0.19 + core::f64::consts::TAU * sample as f64 / 31.0;
        let mut outputs = [[0.0; 2]; 2];
        for (i, candidate) in [&minus, &plus].into_iter().enumerate() {
            let state = [candidate.sample(0, theta), candidate.sample(1, theta)];
            let mut effort = [0.0; 2];
            system.storage().gradient(&state, &mut effort);
            system
                .output_from_effort_into(&effort, &mut outputs[i])
                .unwrap();
        }
        for j in 0..2 {
            let mut predicted = C64::ZERO;
            for k in -11..=11 {
                let (sine, cosine) = det::sin_cos(k as f64 * theta);
                predicted = predicted
                    + tangent.port_log_frequency_derivative[(k + 11) as usize * 2 + j]
                        * C64::new(cosine, sine);
            }
            near(
                predicted,
                C64::from_re((outputs[1][j] - outputs[0][j]) / (2.0 * delta)),
                2e-6,
            );
        }
    }
}

#[test]
fn frequency_tangent_preserves_source_truncation_and_cancellation_refusals() {
    let (system, input, mut base) = nonlinear_fixture();
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 9, 11)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    {
        let mut factor = prepared.factor(0.4, || false).unwrap();
        assert!(matches!(
            factor.forced_frequency_tangent(|| false),
            Err(PhsPortSidebandError::Invalid(_))
        ));
    }
    assert!(matches!(
        prepared.forced_frequency_tangent(|| true),
        Err(PhsPortSidebandError::Cancelled)
    ));
    drop(prepared);
    drop(source);
    // The base's zero padding is harmless, but the original tangent generates
    // omitted harmonics. Keeping only the base's active order must still fail
    // the existing differentiated-equation resolution gate.
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 1, 3)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    assert!(matches!(
        prepared.forced_frequency_tangent(|| false),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::OmittedSidebands { .. }
        ))
    ));
    drop(prepared);
    drop(source);
    // Producer metadata cannot admit a changed base; its original nonlinear
    // balance is checked before a derivative is constructed.
    base.coeffs[0][3] = C64::from_re(0.02);
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 9, 11)).unwrap();
    assert!(matches!(
        source.prepare(&scales(2), || false),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::BaseUnresolved { .. }
        ))
    ));
}

struct ZeroFeedback;
impl AutonomousPortFeedback for ZeroFeedback {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, _: &[f64], _: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        true
    }
}

#[test]
fn autonomous_clock_is_not_silently_treated_as_an_external_frequency_parameter() {
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![1.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let base = orbit(1.0, vec![vec![C64::ZERO; 2]]);
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let mut source = PhsPortSidebandSource::autonomous_feedback(
        &system,
        &ZeroFeedback,
        &derivative,
        &base,
        &controls(1, 1, 1, 1),
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(1), || false).unwrap();
    assert!(matches!(
        prepared.forced_frequency_tangent(|| false),
        Err(PhsPortSidebandError::Invalid(_))
    ));
}
