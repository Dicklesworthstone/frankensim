//! Physical material sensitivities of an original forced periodic response.

use core::cell::Cell;
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput,
    sideband::{
        PhsForcedParameterDirection, PhsPortSidebandControls, PhsPortSidebandError,
        PhsPortSidebandSource, PhsStorageParameterControls,
    },
    solve_forced_feedback_harmonic_balance_scaled_seeded,
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
use fs_phs::{PhsStructureParameterDirection, PortHamiltonian, QuadraticStorage, RelaxationBranch};

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

    // Change relaxation TIME and physical viscous damping together. The
    // normalized memory coordinate is fixed; its original R entry is 1/tau.
    let (dtau, dc) = (0.08, -0.1);
    let dr = [0.0, 0.0, 0.0, 0.0, dc, 0.0, 0.0, 0.0, -dtau / (tau * tau)];
    let result = prepared
        .forced_parameter_tangent(
            PhsForcedParameterDirection {
                structure: PhsStructureParameterDirection {
                    dissipation: Some(&dr),
                    ..Default::default()
                },
                ..Default::default()
            },
            &parameter_controls(3),
            || false,
        )
        .unwrap();
    for k in 0..=2 {
        let s = C64::new(0.0, k as f64 * omega);
        let denominator = C64::ONE + s.scale(tau);
        let dynamic = C64::new(
            stiffness - mass * (k as f64 * omega).powi(2),
            damping * k as f64 * omega,
        ) + (s.scale(tau) / denominator).scale(arm);
        let q = exact(k).0[0];
        let dq = -q * (s.scale(dc) + s.scale(arm * dtau) / (denominator * denominator)) / dynamic;
        let dz =
            (dq / denominator - q * s.scale(dtau) / (denominator * denominator)).scale(arm.sqrt());
        for (j, expected) in [dq, s * dq.scale(mass), dz].into_iter().enumerate() {
            near(
                result.state_parameter_derivative.coefficients[(2 + k) * 3 + j],
                expected,
                2e-11,
            );
        }
        near(result.port_parameter_derivative[2 + k], s * dq, 2e-11);
    }
}

#[test]
fn resonant_damping_reciprocal_gain_and_authored_drive_slopes_match_exact_response() {
    let (omega, stiffness, damping, gain): (f64, f64, f64, f64) = (2.0, 4.0, 0.02, 1.3);
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, gain],
        Box::new(QuadraticStorage::new(vec![stiffness, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    let input_coefficients = [
        C64::from_re(0.02),
        C64::new(0.03, -0.02),
        C64::new(0.01, 0.01),
        C64::ZERO,
    ];
    let input_slopes = [
        C64::from_re(-0.01),
        C64::new(-0.02, 0.03),
        C64::new(0.005, -0.02),
        C64::new(0.01, 0.02),
    ];
    let input = PeriodicPortInput::new(omega, vec![input_coefficients.to_vec()]).unwrap();
    let input_derivative = PeriodicPortInput::new(omega, vec![input_slopes.to_vec()]).unwrap();
    let mut coefficients = vec![vec![C64::ZERO; 4]; 2];
    for k in 0..=3 {
        let w = k as f64 * omega;
        coefficients[0][k] =
            input_coefficients[k].scale(gain) / C64::new(stiffness - w * w, damping * w);
        coefficients[1][k] = C64::new(0.0, w) * coefficients[0][k];
    }
    let base = orbit(omega, coefficients);
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 1, 3, 3)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let mut factor = prepared.factor(0.0, || false).unwrap();
    let (dc, dg) = (-0.03, 0.2);
    let dr = [0.0, 0.0, 0.0, dc];
    let port_direction = [0.0, dg];
    let result = factor
        .forced_parameter_tangent(
            PhsForcedParameterDirection {
                structure: PhsStructureParameterDirection {
                    dissipation: Some(&dr),
                    ports: Some(&port_direction),
                    ..Default::default()
                },
                input: Some(&input_derivative),
                ..Default::default()
            },
            &parameter_controls(2),
            || false,
        )
        .unwrap();
    for k in 0..=3 {
        let w = k as f64 * omega;
        let dq = (input_coefficients[k].scale(dg) + input_slopes[k].scale(gain)
            - C64::new(0.0, dc * w) * base.coeffs[0][k])
            / C64::new(stiffness - w * w, damping * w);
        let dp = C64::new(0.0, w) * dq;
        near(
            result.state_parameter_derivative.coefficients[(3 + k) * 2],
            dq,
            2e-10,
        );
        near(
            result.state_parameter_derivative.coefficients[(3 + k) * 2 + 1],
            dp,
            2e-10,
        );
        near(
            result.port_parameter_derivative[3 + k],
            dp.scale(gain) + base.coeffs[1][k].scale(dg),
            2e-10,
        );
        near(
            result.port_parameter_derivative[3 - k],
            (dp.scale(gain) + base.coeffs[1][k].scale(dg)).conj(),
            2e-10,
        );
    }
    // A changing reciprocal gain appears twice in y=G^T*effort, once through
    // the driven motion and once directly. The SAME factor serves this probe.
    let only_gain = factor
        .forced_parameter_tangent(
            PhsForcedParameterDirection {
                structure: PhsStructureParameterDirection {
                    ports: Some(&port_direction),
                    ..Default::default()
                },
                ..Default::default()
            },
            &parameter_controls(2),
            || false,
        )
        .unwrap();
    near(
        only_gain.state_parameter_derivative.coefficients[4 * 2],
        base.coeffs[0][1].scale(dg / gain),
        2e-11,
    );
    near(
        only_gain.port_parameter_derivative[4],
        base.coeffs[1][1].scale(2.0 * dg),
        2e-11,
    );
    assert!(result.state_parameter_derivative.coefficients[6 * 2].abs() > 1e-4); // newly driven third harmonic
}

fn mixed_model(parameter: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0 + 0.04 * parameter, -1.0 - 0.04 * parameter, 0.0],
        vec![0.0, 0.0, 0.0, 0.3 - 0.04 * parameter],
        vec![
            1.0 + 0.1 * parameter,
            -0.03 * parameter,
            0.07 * parameter,
            1.0 + 0.2 * parameter,
        ],
        Box::new(material(parameter)),
    )
    .unwrap()
}

fn mixed_input(base: &HbOrbit, parameter: f64) -> (PeriodicPortInput, PeriodicPortInput) {
    let q = base.coeffs[0][1];
    let h = base.coeffs[0].len() - 1;
    let mut input = vec![vec![C64::ZERO; h + 1]; 2];
    input[1][1] = q * C64::new(
        4.0 - base.omega * base.omega + 6.0 * q.abs().powi(2),
        0.3 * base.omega,
    );
    input[1][3] = (q * q * q).scale(2.0);
    let mut derivative = vec![vec![C64::ZERO; h + 1]; 2];
    derivative[0][1] = C64::new(0.02, -0.01);
    for k in [1, 3] {
        derivative[1][k] = input[1][k] * C64::new(0.1, -0.15);
    }
    for p in 0..2 {
        for k in 0..=h {
            input[p][k] = input[p][k] + derivative[p][k].scale(parameter);
        }
    }
    (
        PeriodicPortInput::new(base.omega, input).unwrap(),
        PeriodicPortInput::new(base.omega, derivative).unwrap(),
    )
}

#[test]
fn complete_nonlinear_physical_direction_matches_independently_resolved_orbits() {
    let storage = material(0.0);
    let system = mixed_model(0.0);
    let (_, base) = sinusoid(2.0, 0.3, 9);
    let (input, input_derivative) = mixed_input(&base, 0.0);
    let channels = [StressChannelParameterDirection {
        coefficient: 0.5,
        coupling: Some(&[0.1]),
    }];
    let storage_direction = SosModalParameterDirection {
        omegas: Some(&[0.15]),
        channels: Some(&channels),
    };
    let action = |state: &[f64], out: &mut [f64]| {
        storage
            .gradient_parameter_derivative(state, storage_direction, out)
            .is_ok()
    };
    let mut source =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 9, 11)).unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let derivative = prepared
        .forced_parameter_tangent(
            PhsForcedParameterDirection {
                structure: PhsStructureParameterDirection {
                    interconnection: Some(&[0.0, 0.04, -0.04, 0.0]),
                    dissipation: Some(&[0.0, 0.0, 0.0, -0.04]),
                    ports: Some(&[0.1, -0.03, 0.07, 0.2]),
                },
                storage: Some(&action),
                input: Some(&input_derivative),
                feedback: None,
            },
            &parameter_controls(2),
            || false,
        )
        .unwrap();
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
        let system = mixed_model(parameter);
        let input = mixed_input(&base, parameter).0;
        let result = solve_forced_harmonic_balance_scaled_seeded(
            &system,
            &input,
            base.packed_seed(&HbAnchor::Forced { omega: base.omega }),
            &budget,
            &hb_scales,
        )
        .unwrap();
        let output = observe_port_spectrum(&system, &result.orbit, &observation, || false).unwrap();
        (result.orbit, output)
    };
    let mut errors = Vec::new();
    for step in [0.002, 0.001, 0.0005] {
        let (plus, yp) = solve(step);
        let (minus, ym) = solve(-step);
        let mut error = 0.0_f64;
        for k in 0..=9 {
            for j in 0..2 {
                let finite = (plus.coeffs[j][k] - minus.coeffs[j][k]).scale(0.5 / step);
                error = error.max(
                    (finite - derivative.state_parameter_derivative.coefficients[(9 + k) * 2 + j])
                        .abs(),
                );
            }
        }
        for k in 0..=11 {
            for p in 0..2 {
                let finite = if k == 0 {
                    C64::from_re((yp.dc[p] - ym.dc[p]) / (2.0 * step))
                } else {
                    (yp.peak_phasors[p][k - 1] - ym.peak_phasors[p][k - 1])
                        .conj()
                        .scale(0.25 / step)
                };
                error = error
                    .max((finite - derivative.port_parameter_derivative[(11 + k) * 2 + p]).abs());
            }
        }
        errors.push(error);
    }
    assert!(errors[0] > 3.5 * errors[1], "{errors:?}");
    assert!(errors[1] > 3.5 * errors[2], "{errors:?}");
    assert!(errors[2] < 3e-6, "{errors:?}");
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
    let wrong_clock = PeriodicPortInput::new(base.omega * 1.01, vec![vec![C64::ONE]; 2]).unwrap();
    let mut wide = vec![vec![C64::ZERO; 33]; 2];
    wide[1][32] = C64::ONE;
    let wide = PeriodicPortInput::new(base.omega, wide).unwrap();
    for direction in [
        PhsForcedParameterDirection {
            storage: Some(&counted),
            input: Some(&wrong_clock),
            ..Default::default()
        },
        PhsForcedParameterDirection {
            storage: Some(&counted),
            input: Some(&wide),
            ..Default::default()
        },
        PhsForcedParameterDirection {
            structure: PhsStructureParameterDirection {
                interconnection: Some(&[0.0, 1.0, 1.0, 0.0]),
                ..Default::default()
            },
            storage: Some(&counted),
            ..Default::default()
        },
    ] {
        assert!(matches!(
            factor.forced_parameter_tangent(direction, &parameter_controls(2), || false),
            Err(PhsPortSidebandError::Invalid(_))
        ));
    }
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
            .forced_parameter_tangent(
                PhsForcedParameterDirection {
                    structure: PhsStructureParameterDirection {
                        dissipation: Some(&[0.0, 0.0, 0.0, -0.1]),
                        ..Default::default()
                    },
                    storage: Some(&counted),
                    ..Default::default()
                },
                &parameter_controls(2),
                || calls.get() > 10
            )
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
    let new_harmonic = PeriodicPortInput::new(
        base.omega,
        vec![
            vec![C64::ZERO; 4],
            vec![C64::ZERO, C64::ZERO, C64::ZERO, C64::new(0.02, -0.01)],
        ],
    )
    .unwrap();
    let direction = PhsForcedParameterDirection {
        input: Some(&new_harmonic),
        ..Default::default()
    };
    assert!(matches!(
        prepared.forced_parameter_tangent(direction, &parameter_controls(2), || false),
        Err(PhsPortSidebandError::ObservationUnresolved { .. })
    ));
    assert!(matches!(
        prepared.forced_parameter_tangent(direction, &loose, || false),
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

#[test]
fn forced_feedback_chain_is_explicit_and_its_unresolved_harmonics_cannot_disappear() {
    let system = model(SosModalStorage {
        omegas: vec![2.0],
        channels: vec![],
    });
    let (input, base) = sinusoid(0.0, 0.4, 3);
    let feedback = NoFeedback;
    let jacobian = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let calls = Cell::new(0);
    let fixed = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out.fill(0.0);
        true
    };
    let mut source = PhsPortSidebandSource::forced_feedback(
        &system,
        &input,
        &feedback,
        &jacobian,
        &base,
        &controls(2, 2, 1, 3),
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    assert!(matches!(
        prepared.forced_parameter_tangent(
            PhsForcedParameterDirection::default(),
            &parameter_controls(2),
            || false
        ),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    assert_eq!(calls.get(), 0);
    let result = prepared
        .forced_parameter_tangent(
            PhsForcedParameterDirection {
                feedback: Some(&fixed),
                ..Default::default()
            },
            &parameter_controls(2),
            || false,
        )
        .unwrap();
    assert!(
        result
            .state_parameter_derivative
            .coefficients
            .iter()
            .all(|z| *z == C64::ZERO)
    );
    assert_eq!(calls.get(), 64 + 128 + 256);
    let cubic = |_: &[f64], y: &[f64], _: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[0.0, y[0].powi(3)]);
        true
    };
    let direction = PhsForcedParameterDirection {
        feedback: Some(&cubic),
        ..Default::default()
    };
    assert!(matches!(
        prepared.forced_parameter_tangent(direction, &parameter_controls(2), || false),
        Err(PhsPortSidebandError::ObservationUnresolved { .. })
    ));
    let mut loose = parameter_controls(2);
    loose.forcing_rms_tolerances.fill(10.0);
    assert!(matches!(
        prepared.forced_parameter_tangent(direction, &loose, || false),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::OmittedSidebands { .. }
        ))
    ));
    let mut ordinary =
        PhsPortSidebandSource::forced(&system, &input, &base, &controls(2, 2, 1, 3)).unwrap();
    let mut ordinary = ordinary.prepare(&scales(2), || false).unwrap();
    assert!(matches!(
        ordinary.forced_parameter_tangent(
            PhsForcedParameterDirection {
                feedback: Some(&fixed),
                ..Default::default()
            },
            &parameter_controls(2),
            || false,
        ),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    assert_eq!(calls.get(), 64 + 128 + 256);
}

struct FlowFeedback {
    linear: f64,
    cubic: f64,
    position: f64,
}

impl AutonomousPortFeedback for FlowFeedback {
    fn port_count(&self) -> usize {
        1
    }

    fn input(&self, state: &[f64], flow: &[f64], out: &mut [f64]) -> bool {
        out[0] = self.position * state[0] - self.linear * flow[0] - self.cubic * flow[0].powi(3);
        true
    }
}

#[test]
fn feedback_parameter_chain_retains_mass_gain_and_maxwell_history() {
    let (mass, stiffness, damping, arm, tau, gain): (f64, f64, f64, f64, f64, f64) =
        (1.2, 4.2, 0.07, 1.4, 0.23, 1.1);
    let (linear, position) = (0.12, 0.06);
    let (dm, dk, dc, dtau, dg, dl, da) = (0.2, -0.3, 0.02, -0.04, 0.13, -0.03, 0.04);
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, gain],
        Box::new(QuadraticStorage::new(vec![stiffness, 0.0, 0.0, 1.0 / mass], 2).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0],
        stiffness: arm,
        relaxation_time_s: tau,
    }])
    .unwrap();
    let omega = 1.8;
    let forces = [
        C64::from_re(0.02),
        C64::new(0.04, -0.02),
        C64::new(-0.01, 0.03),
        C64::ZERO,
    ];
    let force_slopes = [
        C64::from_re(-0.01),
        C64::new(-0.02, 0.01),
        C64::new(0.01, 0.01),
        C64::new(0.02, -0.03),
    ];
    let exact = |k: usize| {
        let w = k as f64 * omega;
        let s = C64::new(0.0, w);
        let memory = C64::ONE + s.scale(tau);
        let dynamic = C64::new(
            stiffness - gain * position - mass * w * w,
            w * (damping + linear * gain * gain),
        ) + (s.scale(tau) / memory).scale(arm);
        let dynamic_slope = C64::new(
            dk - dg * position - gain * da - dm * w * w,
            w * (dc + dl * gain * gain + 2.0 * linear * gain * dg),
        ) + s.scale(arm * dtau) / (memory * memory);
        let q = forces[k].scale(gain) / dynamic;
        let dq = (forces[k].scale(dg) + force_slopes[k].scale(gain) - dynamic_slope * q) / dynamic;
        (
            [q, s * q.scale(mass), q.scale(arm.sqrt()) / memory],
            [
                dq,
                s * (dq.scale(mass) + q.scale(dm)),
                (dq / memory - q * s.scale(dtau) / (memory * memory)).scale(arm.sqrt()),
            ],
            s * (dq.scale(gain) + q.scale(dg)),
        )
    };
    let mut states = vec![vec![C64::ZERO; 4]; 3];
    for k in 0..=3 {
        for j in 0..3 {
            states[j][k] = exact(k).0[j];
        }
    }
    let base = orbit(omega, states);
    let input = PeriodicPortInput::new(omega, vec![forces.to_vec()]).unwrap();
    let input_slope = PeriodicPortInput::new(omega, vec![force_slopes.to_vec()]).unwrap();
    let feedback = FlowFeedback {
        linear,
        cubic: 0.0,
        position,
    };
    let jacobian = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[position, -linear * gain / mass, 0.0]);
        true
    };
    let storage = |x: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[dk * x[0], -dm * x[1] / (mass * mass), 0.0]);
        true
    };
    let feedback_calls = Cell::new(0);
    let parameter = |x: &[f64], y: &[f64], yp: &[f64], out: &mut [f64]| {
        feedback_calls.set(feedback_calls.get() + 1);
        assert!((y[0] - gain * x[1] / mass).abs() < 2e-14);
        assert!((yp[0] - (dg / mass - gain * dm / (mass * mass)) * x[1]).abs() < 2e-14);
        out[0] = da * x[0] - dl * y[0] - linear * yp[0];
        true
    };
    let dr = [0.0, 0.0, 0.0, 0.0, dc, 0.0, 0.0, 0.0, -dtau / (tau * tau)];
    let dg_map = [0.0, dg, 0.0];
    let mut source = PhsPortSidebandSource::forced_feedback(
        &system,
        &input,
        &feedback,
        &jacobian,
        &base,
        &controls(3, 1, 3, 3),
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    let mut factor = prepared.factor(0.0, || false).unwrap();
    let direction = PhsForcedParameterDirection {
        structure: PhsStructureParameterDirection {
            dissipation: Some(&dr),
            ports: Some(&dg_map),
            ..Default::default()
        },
        storage: Some(&storage),
        input: Some(&input_slope),
        feedback: Some(&parameter),
    };
    let result = factor
        .forced_parameter_tangent(direction, &parameter_controls(3), || false)
        .unwrap();
    for k in 0..=3 {
        for j in 0..3 {
            near(
                result.state_parameter_derivative.coefficients[(3 + k) * 3 + j],
                exact(k).1[j],
                3e-11,
            );
            near(
                result.state_parameter_derivative.coefficients[(3 - k) * 3 + j],
                exact(k).1[j].conj(),
                3e-11,
            );
        }
        near(result.port_parameter_derivative[3 + k], exact(k).2, 3e-11);
    }
    assert_eq!(feedback_calls.get(), 64 + 128 + 256);
    assert_eq!(result.parameter_actions, 2 * (64 + 128 + 256));
    assert!(result.state_parameter_derivative.coefficients[6 * 3 + 2].abs() > 1e-5);

    // Missing/failed/incomplete direct chains cannot borrow the state Jacobian
    // or a previously evaluated callback's scratch. A refusal preserves the LU.
    assert!(matches!(
        factor.forced_parameter_tangent(
            PhsForcedParameterDirection {
                feedback: None,
                ..direction
            },
            &parameter_controls(3),
            || false
        ),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    for mode in 0..3 {
        let bad = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
            if mode == 2 {
                out[0] = f64::INFINITY;
            }
            mode != 0
        };
        assert_eq!(
            factor
                .forced_parameter_tangent(
                    PhsForcedParameterDirection {
                        feedback: Some(&bad),
                        ..direction
                    },
                    &parameter_controls(3),
                    || false
                )
                .unwrap_err(),
            PhsPortSidebandError::SourceRefused
        );
    }
    let before = feedback_calls.get();
    assert_eq!(
        factor
            .forced_parameter_tangent(direction, &parameter_controls(3), || feedback_calls.get()
                > before + 10)
            .unwrap_err(),
        PhsPortSidebandError::Cancelled
    );
    let repeat = factor
        .forced_parameter_tangent(direction, &parameter_controls(3), || false)
        .unwrap();
    assert_eq!(
        result.state_parameter_derivative.coefficients,
        repeat.state_parameter_derivative.coefficients
    );
    assert_eq!(
        result.port_parameter_derivative,
        repeat.port_parameter_derivative
    );
}

fn nonlinear_feedback_model(parameter: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.3 + 0.02 * parameter],
        vec![0.0, 1.1 + 0.13 * parameter],
        Box::new(material(parameter)),
    )
    .unwrap()
}

#[test]
fn nonlinear_feedback_parameter_slope_matches_independently_resolved_original_orbits() {
    let (omega, phase, amplitude): (f64, f64, f64) = (2.05, 0.37, 0.3);
    let q = C64::new(phase.cos(), phase.sin()).scale(0.5 * amplitude);
    let p = C64::new(0.0, omega) * q;
    let mut states = vec![vec![C64::ZERO; 10]; 2];
    states[0][1] = q;
    states[1][1] = p;
    let base = orbit(omega, states);
    let law = |parameter: f64| FlowFeedback {
        linear: 0.12 - 0.03 * parameter,
        cubic: 0.2 + 0.06 * parameter,
        position: 0.07 + 0.04 * parameter,
    };
    let (gain, linear, cubic, position): (f64, f64, f64, f64) = (1.1, 0.12, 0.2, 0.07);
    let mut forces = vec![C64::ZERO; 10];
    forces[1] = (q * C64::new(
        4.0 - gain * position - omega * omega,
        omega * (0.3 + linear * gain * gain),
    ) + q.scale(6.0 * q.abs().powi(2))
        + p.scale(3.0 * cubic * gain.powi(4) * p.abs().powi(2)))
    .scale(1.0 / gain);
    forces[3] =
        ((q * q * q).scale(2.0) + (p * p * p).scale(cubic * gain.powi(4))).scale(1.0 / gain);
    let mut slopes = vec![C64::ZERO; 10];
    slopes[1] = C64::new(0.02, -0.01);
    slopes[3] = C64::new(-0.003, 0.004);
    let drive = |parameter: f64| {
        PeriodicPortInput::new(
            omega,
            vec![
                forces
                    .iter()
                    .zip(&slopes)
                    .map(|(&u, &up)| u + up.scale(parameter))
                    .collect(),
            ],
        )
        .unwrap()
    };
    let input = drive(0.0);
    let input_slope = PeriodicPortInput::new(omega, vec![slopes.clone()]).unwrap();
    let system = nonlinear_feedback_model(0.0);
    let feedback = law(0.0);
    let jacobian = |_: &[f64], y: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[position, -gain * (linear + 3.0 * cubic * y[0] * y[0])]);
        true
    };
    let storage = material(0.0);
    let channels = [StressChannelParameterDirection {
        coefficient: 0.5,
        coupling: Some(&[0.1]),
    }];
    let storage_direction = SosModalParameterDirection {
        omegas: Some(&[0.15]),
        channels: Some(&channels),
    };
    let action = |x: &[f64], out: &mut [f64]| {
        storage
            .gradient_parameter_derivative(x, storage_direction, out)
            .is_ok()
    };
    let parameter = |x: &[f64], y: &[f64], yp: &[f64], out: &mut [f64]| {
        out[0] = 0.04 * x[0] + 0.03 * y[0]
            - 0.06 * y[0].powi(3)
            - (linear + 3.0 * cubic * y[0] * y[0]) * yp[0];
        true
    };
    let mut source = PhsPortSidebandSource::forced_feedback(
        &system,
        &input,
        &feedback,
        &jacobian,
        &base,
        &controls(2, 1, 9, 11),
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let derivative = prepared
        .forced_parameter_tangent(
            PhsForcedParameterDirection {
                structure: PhsStructureParameterDirection {
                    dissipation: Some(&[0.0, 0.0, 0.0, 0.02]),
                    ports: Some(&[0.0, 0.13]),
                    ..Default::default()
                },
                storage: Some(&action),
                input: Some(&input_slope),
                feedback: Some(&parameter),
            },
            &parameter_controls(2),
            || false,
        )
        .unwrap();
    let hb_budget = HbBudget {
        harmonics: 9,
        max_newton: 40,
        tolerance: 2e-12,
        deflate_equilibrium: false,
    };
    let hb_scales = HbScales {
        states: vec![0.3, 0.6],
        equations: vec![0.6, 1.2],
        angular_frequency: omega,
    };
    let observation = PortSpectrumControls {
        harmonics: 11,
        initial_samples: 32,
        max_samples: 128,
        max_source_evaluations: 256,
        max_workspace_bytes: 1 << 20,
        rms_tolerances: vec![2e-12],
    };
    let solve = |value: f64| {
        let system = nonlinear_feedback_model(value);
        let feedback = law(value);
        let jacobian = |_: &[f64], y: &[f64], out: &mut [f64]| {
            out.copy_from_slice(&[
                feedback.position,
                -(1.1 + 0.13 * value) * (feedback.linear + 3.0 * feedback.cubic * y[0] * y[0]),
            ]);
            true
        };
        let result = solve_forced_feedback_harmonic_balance_scaled_seeded(
            &system,
            &drive(value),
            &feedback,
            &jacobian,
            base.packed_seed(&HbAnchor::Forced { omega }),
            &hb_budget,
            &hb_scales,
        )
        .unwrap();
        let output = observe_port_spectrum(&system, &result.orbit, &observation, || false).unwrap();
        (result.orbit, output)
    };
    let mut errors = Vec::new();
    for step in [0.002, 0.001, 0.0005] {
        let (plus, yp) = solve(step);
        let (minus, ym) = solve(-step);
        let mut error = 0.0_f64;
        for k in 0..=9 {
            for j in 0..2 {
                let finite = (plus.coeffs[j][k] - minus.coeffs[j][k]).scale(0.5 / step);
                error = error.max(
                    (finite - derivative.state_parameter_derivative.coefficients[(9 + k) * 2 + j])
                        .abs(),
                );
            }
        }
        for k in 1..=11 {
            let finite = (yp.peak_phasors[0][k - 1] - ym.peak_phasors[0][k - 1])
                .conj()
                .scale(0.25 / step);
            error = error.max((finite - derivative.port_parameter_derivative[11 + k]).abs());
        }
        errors.push(error);
    }
    assert!(errors[0] > 3.5 * errors[1], "{errors:?}");
    assert!(errors[1] > 3.5 * errors[2], "{errors:?}");
    assert!(errors[2] < 3e-6, "{errors:?}");
}
