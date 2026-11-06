use core::cell::Cell;
use fs_couple::phs_orbit::{
    PeriodicPortInput,
    sideband::{
        PhsPortSidebandControls, PhsPortSidebandError, solve_forced_port_sideband_response,
        solve_forced_port_sideband_response_analytic,
    },
};
use fs_math::{c64::C64, det};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_orbit::{
    HbOrbit,
    sideband::{SidebandControls, SidebandError, SidebandScales},
};
use fs_phs::{
    ContinuousStorageTangent, PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage,
};

#[path = "phs_sideband/autonomous.rs"]
mod autonomous;
#[path = "phs_sideband/driven_feedback.rs"]
mod driven_feedback;
#[path = "phs_sideband/resonance.rs"]
mod resonance;
#[path = "phs_sideband/stability.rs"]
mod stability;
#[path = "phs_sideband/transient.rs"]
mod transient;

fn controls(
    states: usize,
    ports: usize,
    harmonics: usize,
    output: usize,
) -> PhsPortSidebandControls {
    PhsPortSidebandControls {
        sideband: SidebandControls {
            harmonics,
            samples_per_period: 128,
            maximum_unknowns: states * (2 * harmonics + 1),
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
        envelope_tolerances: vec![2e-10; ports],
    }
}
fn scales(n: usize) -> SidebandScales {
    SidebandScales {
        states: (0..n).map(|k| 0.5 + k as f64).collect(),
        equations: (0..n).map(|k| 0.7 + 0.4 * k as f64).collect(),
    }
}

fn derivative_controls() -> fs_orbit::sideband::sensitivity::SidebandCarrierDerivativeControls {
    fs_orbit::sideband::sensitivity::SidebandCarrierDerivativeControls {
        maximum_payload_bytes: 32 << 20,
        maximum_work: 1_000_000_000,
        maximum_relative_residual: 1e-12,
        maximum_omitted_residual: 2e-8,
    }
}

fn base_stability_controls(
    scales: &[f64],
) -> fs_couple::phs_orbit::sideband::PhsBaseStabilityControls {
    fs_couple::phs_orbit::sideband::PhsBaseStabilityControls {
        initial_steps_per_period: 128,
        monodromy_change_tolerance: 2e-5,
        spectral_radius_change_tolerance: 1e-5,
        variational: fs_orbit::VariationalShootControls {
            state_scales: scales.to_vec(),
            period_scale: 1.0,
            minimum_autonomous_excursion: 1e-3,
            maximum_dimension: 16,
            maximum_payload_bytes: 2_000_000,
            maximum_work: 1_000_000_000,
            maximum_period_maps: 1,
        },
        floquet: fs_spectral::floquet::FloquetControls {
            spectrum: fs_spectral::dense::DenseSpectrumOptions {
                maximum_dimension: 16,
                maximum_dense_bytes: 1 << 20,
                maximum_work: 100_000_000,
                // The physical map already computed QR. The supplied-root
                // path must not need a second QR allowance.
                qr_max_sweeps_per_window: 0,
                svd_max_sweeps: 100,
                svd_tolerance: 2e-14,
                maximum_relative_residual: 1e-10,
                minimum_eigenbasis_reciprocal_condition: 1e-8,
            },
            unit_circle_guard: 1e-3,
            phase_multiplier_tolerance: 2e-6,
            maximum_phase_residual: 2e-6,
        },
        maximum_workspace_bytes: 8 << 20,
        maximum_work: 1_000_000_000,
    }
}
fn orbit(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    // Deliberately untrusted metadata: the adapter must re-evaluate the source.
    HbOrbit {
        omega,
        coeffs,
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
fn signed_multiforce_resonance_preserves_physical_inputs_velocity_and_phase() {
    let (mass, stiffness, damping, omega) = (2.0, 4.0, 0.3, 1.3);
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 0.0, 1.0, -2.0],
        Box::new(QuadraticStorage::new(vec![stiffness, 0.0, 0.0, 1.0 / mass], 2).unwrap()),
    )
    .unwrap();
    let (a, b) = (C64::new(0.03, -0.01), C64::new(-0.02, 0.005));
    let q = (a - b.scale(2.0)) / C64::new(stiffness - mass * omega * omega, damping * omega);
    let base = orbit(
        omega,
        vec![
            vec![C64::ZERO, q],
            vec![C64::ZERO, C64::new(0.0, mass * omega) * q],
        ],
    );
    let input =
        PeriodicPortInput::new(omega, vec![vec![C64::ZERO, a], vec![C64::ZERO, b]]).unwrap();
    let carrier = 1.8;
    let c = controls(2, 2, 2, 2);
    let mut forcing = vec![C64::ZERO; 10];
    forcing[0] = C64::new(0.4, -0.2); // negative physical frequency
    forcing[1] = C64::new(-0.2, 0.15);
    forcing[4] = C64::new(0.13, 0.09);
    forcing[5] = C64::new(-0.04, -0.02);
    forcing[8] = C64::new(-0.1, 0.3);
    let result = solve_forced_port_sideband_response(
        &system,
        &input,
        &base,
        carrier,
        &forcing,
        &c,
        &scales(2),
        || false,
    )
    .unwrap();
    for band in 0..5 {
        let w = carrier + (band as f64 - 2.0) * omega;
        let f = forcing[2 * band] - forcing[2 * band + 1].scale(2.0);
        let q = f / C64::new(stiffness - mass * w * w, damping * w);
        let velocity = C64::new(0.0, w) * q;
        near(result.state_response.coefficients[2 * band], q, 3e-12);
        near(
            result.state_response.coefficients[2 * band + 1],
            velocity.scale(mass),
            3e-12,
        );
        near(result.port_coefficients[2 * band], velocity, 3e-12);
        near(
            result.port_coefficients[2 * band + 1],
            velocity.scale(-2.0),
            3e-12,
        );
    }
    assert_eq!(result.observation_trace.len(), 3);
    assert_eq!(result.hessian_actions, 2 * (64 + 128 + 256));
    assert!(result.base_residual < 1e-12);
    assert!(result.state_response.omitted_residual < 1e-13);
}

#[test]
fn resonant_probe_keeps_complete_maxwell_memory_and_frequency_dependent_loss() {
    let (mass, stiffness, damping, arm, tau, omega) = (2.0, 4.0, 0.3, 3.0, 0.4, 1.3);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let exact = |w: f64, f: C64| {
        let s = C64::new(0.0, w);
        let relaxation = C64::ONE + s.scale(tau);
        let q = f
            / (C64::new(stiffness - mass * w * w, damping * w)
                + (s.scale(tau) / relaxation).scale(arm));
        [q, s * q.scale(mass), q.scale(det::sqrt(arm)) / relaxation]
    };
    let base_force = C64::new(0.04, -0.02);
    let base = orbit(
        omega,
        exact(omega, base_force)
            .into_iter()
            .map(|z| vec![C64::ZERO, z])
            .collect(),
    );
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, base_force]]).unwrap();
    let forcing = [
        C64::new(-0.2, 0.1),
        C64::new(0.1, -0.07),
        C64::new(0.03, 0.09),
    ];
    let carrier = 1.5;
    let result = solve_forced_port_sideband_response(
        &system,
        &input,
        &base,
        carrier,
        &forcing,
        &controls(3, 1, 1, 1),
        &scales(3),
        || false,
    )
    .unwrap();
    assert_eq!(result.state_dimension, 3);
    for (band, &force) in forcing.iter().enumerate() {
        let w = carrier + (band as f64 - 1.0) * omega;
        let expected = exact(w, force);
        for (j, &z) in expected.iter().enumerate() {
            near(result.state_response.coefficients[3 * band + j], z, 3e-12);
        }
        near(
            result.port_coefficients[band],
            expected[1].scale(1.0 / mass),
            3e-12,
        );
    }
    let collapsed = forcing[1] / C64::new(stiffness - mass * carrier * carrier, damping * carrier);
    assert!((result.state_response.coefficients[3] - collapsed).abs() > 0.05);
}

#[test]
fn native_nonlinear_orbit_converts_frequencies_and_matches_original_physical_law() {
    let (stiffness, cubic, damping, omega, amplitude, phase) = (4.0, 2.0, 0.3, 1.1, 0.3, 0.37);
    let storage = SosModalStorage {
        omegas: vec![det::sqrt(stiffness)],
        channels: vec![StressChannel {
            coefficient: cubic,
            coupling: vec![1.0],
        }],
    };
    // Both reciprocal physical ports are retained: generalized rate and force.
    // The first output is the nonlinear restoring effort, the second velocity.
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(storage),
    )
    .unwrap();
    let q1 = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let p1 = C64::new(0.0, omega) * q1;
    let f1 = q1
        * C64::new(
            stiffness - omega * omega + 0.75 * cubic * amplitude * amplitude,
            damping * omega,
        );
    let f3 = (q1 * q1 * q1).scale(cubic);
    let base = orbit(
        omega,
        vec![
            vec![C64::ZERO, q1, C64::ZERO, C64::ZERO],
            vec![C64::ZERO, p1, C64::ZERO, C64::ZERO],
        ],
    );
    let input = PeriodicPortInput::new(
        omega,
        vec![vec![C64::ZERO; 4], vec![C64::ZERO, f1, C64::ZERO, f3]],
    )
    .unwrap();
    let c = controls(2, 2, 8, 10);
    let (carrier, probe) = (1.9, C64::new(0.1, -0.04));
    let mut forcing = vec![C64::ZERO; 34];
    forcing[17] = probe;
    let result = solve_forced_port_sideband_response(
        &system,
        &input,
        &base,
        carrier,
        &forcing,
        &c,
        &scales(2),
        || false,
    )
    .unwrap();
    let state = |k: isize, j: usize| {
        if k.unsigned_abs() > 8 {
            C64::ZERO
        } else {
            result.state_response.coefficients[(k + 8) as usize * 2 + j]
        }
    };
    assert!(state(-2, 0).abs() > 1e-3);
    assert!(state(2, 0).abs() > 1e-3);
    let mean_stiffness = stiffness + 1.5 * cubic * amplitude * amplitude;
    let modulation = (q1 * q1).scale(3.0 * cubic);
    for k in -10..=10 {
        let effort = state(k, 0).scale(mean_stiffness)
            + modulation * state(k - 2, 0)
            + modulation.conj() * state(k + 2, 0);
        near(
            result.port_coefficients[(k + 10) as usize * 2],
            effort,
            3e-12,
        );
        near(
            result.port_coefficients[(k + 10) as usize * 2 + 1],
            state(k, 1),
            3e-12,
        );
    }
    // Independently differentiate the ORIGINAL nonlinear flow and physical
    // outputs with symmetric finite probes, without using the declared Hessian.
    let epsilon = 1e-4;
    for sample in 0..97 {
        let time = 0.19 + 0.137 * sample as f64;
        let theta = omega * time;
        let x = [base.sample(0, theta), base.sample(1, theta)];
        let mut dx = [0.0; 2];
        let mut derivative = [0.0; 2];
        let mut dy = [0.0; 2];
        for k in -8..=8 {
            let w = carrier + k as f64 * omega;
            let phase = C64::new(det::cos(w * time), det::sin(w * time));
            for j in 0..2 {
                dx[j] += (state(k, j) * phase).re;
                derivative[j] += (C64::new(0.0, w) * state(k, j) * phase).re;
            }
        }
        for k in -10..=10 {
            let angle = (carrier + k as f64 * omega) * time;
            let phase = C64::new(det::cos(angle), det::sin(angle));
            for (j, y) in dy.iter_mut().enumerate() {
                *y += (result.port_coefficients[(k + 10) as usize * 2 + j] * phase).re;
            }
        }
        let du = (probe * C64::new(det::cos(carrier * time), det::sin(carrier * time))).re;
        let mut flows = [[0.0; 2]; 2];
        let mut outputs = [[0.0; 2]; 2];
        for (index, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let perturbed = [x[0] + sign * epsilon * dx[0], x[1] + sign * epsilon * dx[1]];
            let mut effort = [0.0; 2];
            system.storage().gradient(&perturbed, &mut effort);
            system
                .vector_field_from_effort_into(
                    &effort,
                    &[0.0, sign * epsilon * du],
                    &mut flows[index],
                )
                .unwrap();
            system
                .output_from_effort_into(&effort, &mut outputs[index])
                .unwrap();
        }
        for j in 0..2 {
            assert!((derivative[j] - (flows[1][j] - flows[0][j]) / (2.0 * epsilon)).abs() < 3e-8);
            assert!((dy[j] - (outputs[1][j] - outputs[0][j]) / (2.0 * epsilon)).abs() < 2e-9);
        }
    }
}

struct QuarticOnly;
impl Storage for QuarticOnly {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        x[0] * x[0] + 0.75 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 2.0 * x[0] + 3.0 * x[0].powi(3);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for QuarticOnly {
    fn apply(&self, x: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        out[0] = (2.0 + 9.0 * x[0] * x[0]) * v[0];
        true
    }
}

#[test]
fn effort_modulation_hidden_from_flow_tangent_requires_independent_output_resolution() {
    // J-R=0 hides the ENTIRE nonlinear storage tangent from the state ODE.
    // It remains physical at the work-conjugate output y=gradient H.
    let system =
        PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(QuarticOnly)).unwrap();
    let (omega, amplitude, phase, carrier) = (1.1, 0.4, 0.2, 1.7);
    let q1 = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let base = orbit(omega, vec![vec![C64::ZERO, q1]]);
    let input =
        PeriodicPortInput::new(omega, vec![vec![C64::ZERO, C64::new(0.0, omega) * q1]]).unwrap();
    let forcing = [C64::new(0.3, -0.1)];
    let mut c = controls(1, 1, 0, 0);
    let error = solve_forced_port_sideband_response(
        &system,
        &input,
        &base,
        carrier,
        &forcing,
        &c,
        &scales(1),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        PhsPortSidebandError::ObservationUnresolved { .. }
    ));
    c.output_harmonics = 2;
    let result = solve_forced_port_sideband_response(
        &system,
        &input,
        &base,
        carrier,
        &forcing,
        &c,
        &scales(1),
        || false,
    )
    .unwrap();
    let dq = forcing[0] / C64::new(0.0, carrier);
    let modulation = (q1 * q1).scale(9.0);
    near(result.state_response.coefficients[0], dq, 2e-12);
    for (actual, expected) in result.port_coefficients.iter().zip([
        modulation.conj() * dq,
        C64::ZERO,
        dq.scale(2.0 + 4.5 * amplitude * amplitude),
        C64::ZERO,
        modulation * dq,
    ]) {
        near(*actual, expected, 2e-12);
    }
}

#[test]
fn source_refusal_cancellation_and_budget_exhaustion_publish_no_partial_response() {
    let system = fs_phs::mass_spring_damper(2.0, 4.0, 0.3).unwrap();
    let input = PeriodicPortInput::new(1.3, vec![vec![C64::ZERO]]).unwrap();
    let base = orbit(1.3, vec![vec![C64::ZERO], vec![C64::ZERO]]);
    let original = base.coeffs.clone();
    let c = controls(2, 1, 0, 0);
    let forcing = [C64::new(0.2, -0.1)];
    for incomplete in [false, true] {
        let action = |_: &[f64], _: &[f64], out: &mut [f64]| {
            if incomplete {
                out[0] = 1.0;
                true
            } else {
                false
            }
        };
        assert!(matches!(
            solve_forced_port_sideband_response_analytic(
                &system,
                &input,
                &base,
                1.4,
                &forcing,
                &c,
                &scales(2),
                &action,
                || false,
            ),
            Err(PhsPortSidebandError::Sideband(
                SidebandError::TangentRefused { .. }
            ))
        ));
    }
    // A valid tangent on basis directions can still refuse an observation
    // direction. The declared source must not be replaced after preparation.
    let action = |_: &[f64], v: &[f64], out: &mut [f64]| {
        if v != [1.0, 0.0] && v != [0.0, 1.0] {
            return false;
        }
        out[0] = 4.0 * v[0];
        out[1] = 0.5 * v[1];
        true
    };
    assert!(matches!(
        solve_forced_port_sideband_response_analytic(
            &system,
            &input,
            &base,
            1.4,
            &forcing,
            &c,
            &scales(2),
            &action,
            || false,
        ),
        Err(PhsPortSidebandError::SourceRefused)
    ));
    assert!(matches!(
        solve_forced_port_sideband_response(
            &system,
            &input,
            &base,
            1.4,
            &forcing,
            &c,
            &scales(2),
            || true,
        ),
        Err(PhsPortSidebandError::Cancelled)
    ));
    let observing = Cell::new(false);
    let action = |_: &[f64], v: &[f64], out: &mut [f64]| {
        if v != [1.0, 0.0] && v != [0.0, 1.0] {
            observing.set(true);
        }
        out[0] = 4.0 * v[0];
        out[1] = 0.5 * v[1];
        true
    };
    assert!(matches!(
        solve_forced_port_sideband_response_analytic(
            &system,
            &input,
            &base,
            1.4,
            &forcing,
            &c,
            &scales(2),
            &action,
            || observing.get(),
        ),
        Err(PhsPortSidebandError::Cancelled)
    ));
    let mut too_small = c;
    too_small.max_hessian_actions = 1;
    let calls = Cell::new(0);
    let unused = |_: &[f64], _: &[f64], _: &mut [f64]| {
        calls.set(calls.get() + 1);
        false
    };
    assert!(matches!(
        solve_forced_port_sideband_response_analytic(
            &system,
            &input,
            &base,
            1.4,
            &forcing,
            &too_small,
            &scales(2),
            &unused,
            || false,
        ),
        Err(PhsPortSidebandError::Budget(_))
    ));
    assert_eq!(calls.get(), 0);
    assert_eq!(base.coeffs, original);
}
