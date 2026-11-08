//! Physical material sensitivities of an original forced periodic response.

use core::cell::Cell;
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput,
    sideband::{
        PhsPortSidebandControls, PhsPortSidebandError, PhsPortSidebandSource,
        PhsStorageParameterControls,
    },
    solve_forced_harmonic_balance_scaled_seeded,
    spectrum::{PortSpectrumControls, observe_port_spectrum},
};
use fs_math::c64::C64;
use fs_nlmodal::{
    SosModalParameterDirection, SosModalStorage, StressChannel, StressChannelParameterDirection,
};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbScales,
    sideband::{SidebandControls, SidebandError, SidebandScales},
};
use fs_phs::{PortHamiltonian, RelaxationBranch};

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

fn parameter_controls(n: usize) -> PhsStorageParameterControls {
    PhsStorageParameterControls {
        forcing_rms_tolerances: vec![2e-10; n],
        maximum_parameter_actions: 1024,
    }
}

fn orbit(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

fn model(storage: SosModalStorage) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.3],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(storage),
    )
    .unwrap()
}

fn sinusoid(cubic: f64, amplitude: f64, h: usize) -> (PeriodicPortInput, HbOrbit) {
    let (omega, phase): (f64, f64) = (2.05, 0.37);
    let q = C64::new(phase.cos(), phase.sin()).scale(0.5 * amplitude);
    let mut state = vec![vec![C64::ZERO; h + 1]; 2];
    state[0][1] = q;
    state[1][1] = C64::new(0.0, omega) * q;
    let mut drive = vec![vec![C64::ZERO; h + 1]; 2];
    drive[1][1] = q * C64::new(
        4.0 - omega * omega + 0.75 * cubic * amplitude * amplitude,
        0.3 * omega,
    );
    if h >= 3 {
        drive[1][3] = (q * q * q).scale(cubic);
    }
    (
        PeriodicPortInput::new(omega, drive).unwrap(),
        orbit(omega, state),
    )
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} vs {expected:?}"
    );
}

#[test]
fn zero_cubic_baseline_keeps_the_exact_third_harmonic_and_direct_force_derivative() {
    let storage = SosModalStorage {
        omegas: vec![2.0],
        channels: vec![StressChannel {
            coefficient: 0.0,
            coupling: vec![1.0],
        }],
    };
    let system = model(SosModalStorage {
        omegas: vec![2.0],
        channels: vec![StressChannel {
            coefficient: 0.0,
            coupling: vec![1.0],
        }],
    });
    let channels = [StressChannelParameterDirection {
        coefficient: 1.0,
        coupling: None,
    }];
    let direction = SosModalParameterDirection {
        omegas: None,
        channels: Some(&channels),
    };
    let action = |x: &[f64], out: &mut [f64]| {
        storage
            .gradient_parameter_derivative(x, direction, out)
            .is_ok()
    };
    let mut measured = Vec::new();
    for amplitude in [0.4, 0.2] {
        let (input, base) = sinusoid(0.0, amplitude, 3);
        let mut source =
            PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 3, 3)).unwrap();
        let mut prepared = source.prepare(&scales(2), || false).unwrap();
        let mut factor = prepared.factor(0.0, || false).unwrap();
        let result = factor
            .forced_storage_parameter_tangent(&action, &parameter_controls(2), || false)
            .unwrap();
        let q = base.coeffs[0][1];
        for harmonic in 0..=3 {
            // Independently differentiate the original cubic force at beta=0:
            // q^3 has coefficients 3|Q|² Q at k=1 and Q³ at k=3.
            let direct = match harmonic {
                1 => q.scale(3.0 * q.abs().powi(2)),
                3 => q * q * q,
                _ => C64::ZERO,
            };
            let w = harmonic as f64 * base.omega;
            let dq = -direct / C64::new(4.0 - w * w, 0.3 * w);
            let dp = C64::new(0.0, w) * dq;
            near(
                result.state_parameter_derivative.coefficients[(3 + harmonic) * 2],
                dq,
                2e-11,
            );
            near(
                result.state_parameter_derivative.coefficients[(3 + harmonic) * 2 + 1],
                dp,
                2e-11,
            );
            near(
                result.port_parameter_derivative[(3 + harmonic) * 2],
                dq.scale(4.0) + direct,
                2e-11,
            );
            near(
                result.port_parameter_derivative[(3 + harmonic) * 2 + 1],
                dp,
                2e-11,
            );
            near(
                result.port_parameter_derivative[(3 - harmonic) * 2],
                (dq.scale(4.0) + direct).conj(),
                2e-11,
            );
        }
        assert!(result.state_parameter_derivative.coefficients[6 * 2].abs() > 1e-5);
        assert_eq!(result.forcing_projection_trace.len(), 3);
        assert_eq!(result.observation_trace.len(), 3);
        assert_eq!(result.parameter_actions, 2 * (64 + 128 + 256));
        assert_eq!(result.hessian_actions, 2 * (64 + 128 + 256));
        assert!(result.combined_omitted_equation_defect < 2e-8);
        measured.push(result.state_parameter_derivative.coefficients[4 * 2]);

        // The same factor serves a distinct material direction. A zero action
        // remains exactly zero; no previous response or action budget leaks.
        let zero = |_: &[f64], out: &mut [f64]| {
            out.fill(0.0);
            true
        };
        let zero = factor
            .forced_storage_parameter_tangent(&zero, &parameter_controls(2), || false)
            .unwrap();
        assert!(
            zero.state_parameter_derivative
                .coefficients
                .iter()
                .all(|&z| z == C64::ZERO)
        );
        assert!(
            zero.port_parameter_derivative
                .iter()
                .all(|&z| z == C64::ZERO)
        );
    }
    // Cubic material sensitivity vanishes as amplitude cubed in the linear limit.
    near(measured[0], measured[1].scale(8.0), 2e-11);
}

fn material(parameter: f64) -> SosModalStorage {
    SosModalStorage {
        omegas: vec![2.0 + 0.15 * parameter],
        channels: vec![StressChannel {
            coefficient: 2.0 + 0.5 * parameter,
            coupling: vec![1.0 + 0.1 * parameter],
        }],
    }
}

#[test]
fn nonlinear_material_sensitivity_matches_independently_resolved_original_orbits() {
    let storage = material(0.0);
    let system = model(material(0.0));
    let (input, base) = sinusoid(2.0, 0.3, 9);
    let channels = [StressChannelParameterDirection {
        coefficient: 0.5,
        coupling: Some(&[0.1]),
    }];
    let direction = SosModalParameterDirection {
        omegas: Some(&[0.15]),
        channels: Some(&channels),
    };
    let action = |x: &[f64], out: &mut [f64]| {
        storage
            .gradient_parameter_derivative(x, direction, out)
            .is_ok()
    };
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 9, 11)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let derivative = prepared
        .forced_storage_parameter_tangent(&action, &parameter_controls(2), || false)
        .unwrap();
    assert!(derivative.state_parameter_derivative.coefficients[(9 + 3) * 2].abs() > 1e-5);
    let budget = HbBudget {
        harmonics: 9,
        max_newton: 40,
        tolerance: 2e-12,
        deflate_equilibrium: false,
    };
    let hb_scales = HbScales {
        states: vec![0.3, 0.6],
        equations: vec![0.6, 1.2],
        angular_frequency: base.omega,
    };
    let observation = PortSpectrumControls {
        harmonics: 11,
        initial_samples: 32,
        max_samples: 128,
        max_source_evaluations: 256,
        max_workspace_bytes: 1 << 20,
        rms_tolerances: vec![2e-12; 2],
    };
    let solve = |parameter: f64| {
        let system = model(material(parameter)); // J/R/G and authored drive remain fixed.
        let result = solve_forced_harmonic_balance_scaled_seeded(
            &system,
            &input,
            base.packed_seed(&HbAnchor::Forced { omega: base.omega }),
            &budget,
            &hb_scales,
        )
        .unwrap();
        let outputs =
            observe_port_spectrum(&system, &result.orbit, &observation, || false).unwrap();
        (result.orbit, outputs)
    };
    let mut errors = Vec::new();
    for step in [0.004, 0.002, 0.001] {
        let (plus, yp) = solve(step);
        let (minus, ym) = solve(-step);
        let mut error = 0.0_f64;
        for harmonic in 0..=9 {
            for state in 0..2 {
                let finite = (plus.coeffs[state][harmonic] - minus.coeffs[state][harmonic])
                    .scale(0.5 / step);
                let exact =
                    derivative.state_parameter_derivative.coefficients[(9 + harmonic) * 2 + state];
                error = error.max((finite - exact).abs());
            }
        }
        for harmonic in 0..=11 {
            for port in 0..2 {
                let finite = if harmonic == 0 {
                    C64::from_re((yp.dc[port] - ym.dc[port]) / (2.0 * step))
                } else {
                    (yp.peak_phasors[port][harmonic - 1] - ym.peak_phasors[port][harmonic - 1])
                        .conj()
                        .scale(0.25 / step)
                };
                let exact = derivative.port_parameter_derivative[(11 + harmonic) * 2 + port];
                error = error.max((finite - exact).abs());
            }
        }
        errors.push(error);
    }
    assert!(errors[0] > 3.5 * errors[1], "{errors:?}");
    assert!(errors[1] > 3.5 * errors[2], "{errors:?}");
    assert!(errors[2] < 3e-6, "{errors:?}");
}

#[test]
fn material_direction_retains_the_original_maxwell_memory_response() {
    let (mass, stiffness, damping, arm, tau) = (2.0, 4.0, 0.3, 3.0, 0.4);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let omega = 1.4;
    let forces = [
        C64::from_re(0.02),
        C64::new(0.04, -0.02),
        C64::new(-0.01, 0.03),
    ];
    let exact = |harmonic: usize| {
        let w = harmonic as f64 * omega;
        let s = C64::new(0.0, w);
        let relaxation = C64::ONE + s.scale(tau);
        let dynamic = C64::new(stiffness - mass * w * w, damping * w)
            + (s.scale(tau) / relaxation).scale(arm);
        let q = forces[harmonic] / dynamic;
        let dq = -q / dynamic;
        (
            [q, s * q.scale(mass), q.scale(arm.sqrt()) / relaxation],
            [dq, s * dq.scale(mass), dq.scale(arm.sqrt()) / relaxation],
        )
    };
    let mut states = vec![vec![C64::ZERO; 3]; 3];
    for k in 0..=2 {
        for j in 0..3 {
            states[j][k] = exact(k).0[j];
        }
    }
    let base = orbit(omega, states);
    let input = PeriodicPortInput::new(omega, vec![forces.to_vec()]).unwrap();
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(3, 1, 2, 2)).unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    // Change only the relaxed spring stiffness; the memory energy/law is fixed.
    let action = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[x[0], 0.0, 0.0]);
        true
    };
    let result = prepared
        .forced_storage_parameter_tangent(&action, &parameter_controls(3), || false)
        .unwrap();
    assert_eq!(result.state_dimension, 3);
    for k in 0..=2 {
        for j in 0..3 {
            near(
                result.state_parameter_derivative.coefficients[(2 + k) * 3 + j],
                exact(k).1[j],
                2e-11,
            );
        }
        near(
            result.port_parameter_derivative[2 + k],
            exact(k).1[1].scale(1.0 / mass),
            2e-11,
        );
    }
}

#[test]
fn source_resolution_refusal_cancellation_and_retry_remain_explicit() {
    let system = model(material(0.0));
    let (input, base) = sinusoid(2.0, 0.3, 9);
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 9, 11)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let action = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[x[0].powi(3), 0.0]);
        true
    };
    {
        let mut nonzero = prepared.factor(0.2, || false).unwrap();
        assert!(matches!(
            nonzero.forced_storage_parameter_tangent(&action, &parameter_controls(2), || false),
            Err(PhsPortSidebandError::Invalid(_))
        ));
    }
    let mut factor = prepared.factor(0.0, || false).unwrap();
    let calls = Cell::new(0);
    let counted = |x: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        action(x, out)
    };
    let mut invalid = parameter_controls(2);
    invalid.forcing_rms_tolerances.pop();
    assert!(matches!(
        factor.forced_storage_parameter_tangent(&counted, &invalid, || false),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    invalid = parameter_controls(2);
    invalid.maximum_parameter_actions = 14 * 64 - 1;
    assert!(matches!(
        factor.forced_storage_parameter_tangent(&counted, &invalid, || false),
        Err(PhsPortSidebandError::Budget(_))
    ));
    assert_eq!(calls.get(), 0);
    for mode in 0..3 {
        let bad = |_: &[f64], out: &mut [f64]| {
            if mode == 0 {
                false
            } else if mode == 1 {
                out[0] = 0.0;
                true
            } else {
                out.fill(f64::INFINITY);
                true
            }
        };
        assert_eq!(
            factor
                .forced_storage_parameter_tangent(&bad, &parameter_controls(2), || false)
                .unwrap_err(),
            PhsPortSidebandError::SourceRefused
        );
    }
    assert_eq!(
        factor
            .forced_storage_parameter_tangent(&counted, &parameter_controls(2), || calls.get() > 10)
            .unwrap_err(),
        PhsPortSidebandError::Cancelled
    );
    let first = factor
        .forced_storage_parameter_tangent(&action, &parameter_controls(2), || false)
        .unwrap();
    let repeat = factor
        .forced_storage_parameter_tangent(&action, &parameter_controls(2), || false)
        .unwrap();
    assert_eq!(
        first.state_parameter_derivative.coefficients,
        repeat.state_parameter_derivative.coefficients
    );
    assert_eq!(
        first.port_parameter_derivative,
        repeat.port_parameter_derivative
    );

    // A source can require more harmonics than its linear base. Neither a
    // tight physical forcing gate nor a deliberately loose one may erase it.
    let linear = model(SosModalStorage {
        omegas: vec![2.0],
        channels: vec![],
    });
    let (input, base) = sinusoid(0.0, 0.4, 3);
    let mut source =
        PhsPortSidebandSource::forced(&linear, &input, &base, &controls(2, 2, 1, 3)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    assert!(matches!(
        prepared.forced_storage_parameter_tangent(&action, &parameter_controls(2), || false),
        Err(PhsPortSidebandError::ObservationUnresolved { .. })
    ));
    let mut loose = parameter_controls(2);
    loose.forcing_rms_tolerances.fill(1.0);
    assert!(matches!(
        prepared.forced_storage_parameter_tangent(&action, &loose, || false),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::OmittedSidebands { .. }
        ))
    ));
}

struct NoFeedback;
impl AutonomousPortFeedback for NoFeedback {
    fn port_count(&self) -> usize {
        2
    }
    fn input(&self, _: &[f64], _: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        true
    }
}

#[test]
fn feedback_sources_cannot_omit_their_direct_material_chain_rule() {
    let system = model(material(0.0));
    let base = orbit(2.0, vec![vec![C64::ZERO; 2]; 2]);
    let feedback = NoFeedback;
    let jacobian = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let mut source = PhsPortSidebandSource::autonomous_feedback(
        &system,
        &feedback,
        &jacobian,
        &base,
        &controls(2, 2, 1, 1),
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let calls = Cell::new(0);
    let action = |_: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out.fill(0.0);
        true
    };
    assert!(matches!(
        prepared.forced_storage_parameter_tangent(&action, &parameter_controls(2), || false),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    assert_eq!(calls.get(), 0);
}
