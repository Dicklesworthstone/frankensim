use super::{QuarticOnly, controls, derivative_controls, near, orbit, scales};
use core::cell::Cell;
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput,
    sideband::{
        PhsPortSidebandError, PhsPortSidebandSource, solve_forced_feedback_port_sideband_response,
        solve_forced_feedback_port_sideband_response_analytic,
    },
};
use fs_math::{c64::C64, det};
use fs_orbit::sideband::{SidebandError, refinement::SidebandRefinementControls};
use fs_phs::{PortHamiltonian, RelaxationBranch};

const MASS: f64 = 2.0;
const STIFFNESS: f64 = 4.0;
const DAMPING: f64 = 0.3;
const ARM: f64 = 3.0;
const TAU: f64 = 0.4;
const OMEGA: f64 = 1.3;

struct VelocityFeedback(f64);
impl AutonomousPortFeedback for VelocityFeedback {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, _: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        out[0] = -self.0 * flows[0];
        true
    }
}

// Independent physical transfer: the ORIGINAL Maxwell branch remains dynamic,
// while velocity feedback adds g to the closed-loop damping, not to material R.
fn exact(frequency: f64, force: C64, gain: f64) -> [C64; 3] {
    let s = C64::new(0.0, frequency);
    let relaxation = C64::ONE + s.scale(TAU);
    let q = force
        / (C64::new(
            STIFFNESS - MASS * frequency * frequency,
            (DAMPING + gain) * frequency,
        ) + (s.scale(TAU) / relaxation).scale(ARM));
    [q, s * q.scale(MASS), q.scale(det::sqrt(ARM)) / relaxation]
}

fn full_hessian(_: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
    out[0] = (STIFFNESS + ARM) * direction[0] - det::sqrt(ARM) * direction[2];
    out[1] = direction[1] / MASS;
    out[2] = direction[2] - det::sqrt(ARM) * direction[0];
    true
}

fn exact_derivative(frequency: f64, force: C64, force_derivative: C64, gain: f64) -> [C64; 3] {
    let state = exact(frequency, force, gain);
    let s = C64::new(0.0, frequency);
    let relaxation = C64::ONE + s.scale(TAU);
    let denominator = C64::new(
        STIFFNESS - MASS * frequency * frequency,
        (DAMPING + gain) * frequency,
    ) + (s.scale(TAU) / relaxation).scale(ARM);
    let denominator_derivative = C64::new(-2.0 * MASS * frequency, DAMPING + gain)
        + C64::new(0.0, ARM * TAU) / (relaxation * relaxation);
    let q = (force_derivative - denominator_derivative * state[0]) / denominator;
    [
        q,
        C64::new(0.0, MASS) * (state[0] + q.scale(frequency)),
        q.scale(det::sqrt(ARM)) / relaxation - C64::new(0.0, TAU) * state[2] / relaxation,
    ]
}

#[test]
fn driven_feedback_sweep_matches_original_maxwell_transfer_and_reuses_preparation() {
    let gain = 0.7;
    let feedback = VelocityFeedback(gain);
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: ARM,
            relaxation_time_s: TAU,
        }])
        .unwrap();
    let base_force = C64::new(0.04, -0.02);
    let base = orbit(
        OMEGA,
        exact(OMEGA, base_force, gain)
            .into_iter()
            .map(|z| vec![C64::ZERO, z])
            .collect(),
    );
    let original = base.coeffs.clone();
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO, base_force]]).unwrap();
    let calls = Cell::new(0_usize);
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out.fill(0.0);
        out[1] = -gain / MASS;
        true
    };
    let c = controls(3, 1, 1, 1);
    let forcing = [
        C64::new(-0.2, 0.1),
        C64::new(0.1, -0.07),
        C64::new(0.03, 0.09),
    ];
    let rotations = [C64::ONE, C64::new(-0.3, 0.7)];
    let mut source =
        PhsPortSidebandSource::forced_feedback(&system, &input, &feedback, &derivative, &base, &c)
            .unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    let preparations = 3 * c.sideband.samples_per_period + 1;
    assert_eq!(calls.get(), preparations);
    for carrier in [0.11, 1.2, -2.4] {
        let mut factor = prepared.factor(carrier, || false).unwrap();
        for rotation in rotations {
            let forcing = forcing.map(|z| z * rotation);
            let result = factor.solve(&forcing, || false).unwrap();
            for (band, force) in forcing.into_iter().enumerate() {
                let w = carrier + (band as f64 - 1.0) * OMEGA;
                let expected = exact(w, force, gain);
                for (j, expected) in expected.into_iter().enumerate() {
                    near(
                        result.state_response.coefficients[3 * band + j],
                        expected,
                        4e-12,
                    );
                }
                near(
                    result.port_coefficients[band],
                    expected[1].scale(1.0 / MASS),
                    4e-12,
                );
            }
            assert_eq!(calls.get(), preparations);
            assert!(result.base_residual < 1e-12);
            assert_eq!(result.hessian_actions, 896);
        }
    }
    // The two public convenience entries use the same physical source, complete
    // memory and output conventions, including the explicitly declared Hessian.
    let native = solve_forced_feedback_port_sideband_response(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        1.2,
        &forcing,
        &c,
        &scales(3),
        || false,
    )
    .unwrap();
    let analytic = solve_forced_feedback_port_sideband_response_analytic(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        1.2,
        &forcing,
        &c,
        &scales(3),
        &full_hessian,
        || false,
    )
    .unwrap();
    for (a, b) in native
        .state_response
        .coefficients
        .iter()
        .zip(&analytic.state_response.coefficients)
    {
        near(*a, *b, 4e-12);
    }
    for (a, b) in native
        .port_coefficients
        .iter()
        .zip(&analytic.port_coefficients)
    {
        near(*a, *b, 4e-12);
    }
    let physical = exact(1.2, forcing[1], gain);
    assert!((physical[0] - exact(1.2, forcing[1], 0.0)[0]).abs() > 1e-3);
    assert_eq!(base.coeffs, original);
}

#[test]
fn prepared_physical_observation_retains_hessian_modulation_hidden_by_the_flow_tangent() {
    let system =
        PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(QuarticOnly)).unwrap();
    let (omega, amplitude, phase) = (1.1, 0.4, 0.2);
    let q1 = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let base = orbit(omega, vec![vec![C64::ZERO, q1]]);
    let input =
        PeriodicPortInput::new(omega, vec![vec![C64::ZERO, C64::new(0.0, omega) * q1]]).unwrap();
    let c = controls(1, 1, 0, 2);
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    let mut prepared = source.prepare(&scales(1), || false).unwrap();
    for carrier in [1.7, -0.8] {
        let mut factor = prepared.factor(carrier, || false).unwrap();
        for force in [C64::new(0.3, -0.1), C64::new(-0.07, 0.11)] {
            let result = factor.solve(&[force], || false).unwrap();
            let dq = force / C64::new(0.0, carrier);
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
    }
}

#[test]
fn prepared_driven_feedback_refuses_invalid_sources_and_keeps_factor_after_cancelled_observation() {
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING).unwrap();
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO]]).unwrap();
    let base = orbit(OMEGA, vec![vec![C64::ZERO], vec![C64::ZERO]]);
    let original = base.coeffs.clone();
    let feedback = VelocityFeedback(0.7);
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        out[1] = -feedback.0 / MASS;
        true
    };
    let c = controls(2, 1, 0, 0);
    let bad_clock = PeriodicPortInput::new(OMEGA + 1e-6, vec![vec![C64::ZERO]]).unwrap();
    assert!(matches!(
        PhsPortSidebandSource::forced_feedback(
            &system,
            &bad_clock,
            &feedback,
            &derivative,
            &base,
            &c,
        ),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    let mut unresolved_drive = vec![C64::ZERO; 65];
    unresolved_drive[64] = C64::from_re(1e-100);
    let unresolved_drive = PeriodicPortInput::new(OMEGA, vec![unresolved_drive]).unwrap();
    assert!(matches!(
        PhsPortSidebandSource::forced_feedback(
            &system,
            &unresolved_drive,
            &feedback,
            &derivative,
            &base,
            &c,
        ),
        Err(PhsPortSidebandError::Invalid(_))
    ));
    let hessian_calls = Cell::new(0_usize);
    let hessian = |_: &[f64], v: &[f64], out: &mut [f64]| {
        hessian_calls.set(hessian_calls.get() + 1);
        out[0] = STIFFNESS * v[0];
        out[1] = v[1] / MASS;
        true
    };
    let mut insufficient = c.clone();
    insufficient.max_workspace_bytes = 0;
    assert!(matches!(
        PhsPortSidebandSource::forced_feedback_analytic(
            &system,
            &input,
            &feedback,
            &derivative,
            &base,
            &insufficient,
            &hessian,
        ),
        Err(PhsPortSidebandError::Budget(_))
    ));
    assert_eq!(hessian_calls.get(), 0);
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    let mut source =
        PhsPortSidebandSource::forced_feedback(&system, &input, &feedback, &incomplete, &base, &c)
            .unwrap();
    assert!(matches!(
        source.prepare(&scales(2), || false),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::TangentRefused { .. }
        ))
    ));
    let refused = |_: &[f64], _: &[f64], _: &mut [f64]| false;
    let mut source = PhsPortSidebandSource::forced_feedback_analytic(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        &c,
        &refused,
    )
    .unwrap();
    assert!(matches!(
        source.prepare(&scales(2), || false),
        Err(PhsPortSidebandError::Sideband(
            SidebandError::TangentRefused { .. }
        ))
    ));
    let mut source = PhsPortSidebandSource::forced_feedback_analytic(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        &c,
        &hessian,
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(2), || false).unwrap();
    let prepared_actions = hessian_calls.get();
    assert_eq!(
        prepared_actions,
        2 * (3 * c.sideband.samples_per_period + 1)
    );
    let mut factor = prepared.factor(1.4, || false).unwrap();
    let force = C64::new(0.2, -0.1);
    assert!(matches!(
        factor.solve(&[force], || hessian_calls.get() > prepared_actions),
        Err(PhsPortSidebandError::Cancelled)
    ));
    assert_eq!(hessian_calls.get(), prepared_actions + 1);
    let result = factor.solve(&[force], || false).unwrap();
    let expected = force / C64::new(STIFFNESS - MASS * 1.4 * 1.4, (DAMPING + feedback.0) * 1.4);
    near(result.state_response.coefficients[0], expected, 4e-12);
    near(
        result.port_coefficients[0],
        C64::new(0.0, 1.4) * expected,
        4e-12,
    );
    assert_eq!(base.coeffs, original);
}

#[test]
fn adaptive_driven_feedback_keeps_original_input_support_and_maxwell_loss() {
    let gain = 0.7;
    let feedback = VelocityFeedback(gain);
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: ARM,
            relaxation_time_s: TAU,
        }])
        .unwrap();
    let base_force = C64::new(0.04, -0.02);
    let base = orbit(
        OMEGA,
        exact(OMEGA, base_force, gain)
            .into_iter()
            .map(|z| vec![C64::ZERO, z])
            .collect(),
    );
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO, base_force]]).unwrap();
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        out[1] = -gain / MASS;
        true
    };
    let mut c = controls(3, 1, 1, 1);
    c.sideband.maximum_unknowns = 3 * 7;
    let refinement = SidebandRefinementControls {
        maximum_harmonics: 3,
        equation_tolerances: vec![2e-10; 3],
        state_tolerances: vec![2e-10; 3],
    };
    let forcing = [
        C64::new(-0.2, 0.1),
        C64::new(0.1, -0.07),
        C64::new(0.03, 0.09),
    ];
    let mut source = PhsPortSidebandSource::forced_feedback_analytic(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        &c,
        &full_hessian,
    )
    .unwrap();
    for carrier in [1.2, -0.4] {
        let adaptive = source
            .solve_adaptive(carrier, &forcing, &scales(3), &refinement, || false)
            .unwrap();
        assert_eq!(
            adaptive
                .refinement_trace
                .iter()
                .map(|check| check.harmonics)
                .collect::<Vec<_>>(),
            [1, 3]
        );
        let result = adaptive.response;
        assert_eq!(result.state_harmonics, 3);
        assert_eq!(result.output_harmonics, 1);
        for k in -3_isize..=3 {
            let force = if k.unsigned_abs() <= 1 {
                forcing[(k + 1) as usize]
            } else {
                C64::ZERO
            };
            let expected = exact(carrier + k as f64 * OMEGA, force, gain);
            for (j, expected) in expected.into_iter().enumerate() {
                near(
                    result.state_response.coefficients[(k + 3) as usize * 3 + j],
                    expected,
                    5e-12,
                );
            }
            if k.unsigned_abs() <= 1 {
                near(
                    result.port_coefficients[(k + 1) as usize],
                    expected[1].scale(1.0 / MASS),
                    5e-12,
                );
            }
        }
    }
}

#[test]
fn adaptive_original_output_refinement_refusal_grid_budgets_and_cancellation_remain_independent() {
    let system =
        PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(QuarticOnly)).unwrap();
    let (omega, amplitude, phase, carrier) = (1.1, 0.4, 0.2, 1.7);
    let q1 = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let base = orbit(omega, vec![vec![C64::ZERO, q1]]);
    let original = base.coeffs.clone();
    let input =
        PeriodicPortInput::new(omega, vec![vec![C64::ZERO, C64::new(0.0, omega) * q1]]).unwrap();
    let force = C64::new(0.3, -0.1);
    let refinement = SidebandRefinementControls {
        maximum_harmonics: 1,
        equation_tolerances: vec![2e-10],
        state_tolerances: vec![2e-10],
    };
    let mut c = controls(1, 1, 0, 0);
    c.sideband.maximum_unknowns = 3;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    // A(t)=0 makes the adaptive state solve exact at every order. Its success
    // must not conceal the nonzero +/-2 OUTPUT generated by Hxx(t).
    assert!(matches!(
        source.solve_adaptive(carrier, &[force], &scales(1), &refinement, || false),
        Err(PhsPortSidebandError::ObservationUnresolved { .. })
    ));
    c.output_harmonics = 2;
    let observing = Cell::new(false);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        if v != [1.0] {
            observing.set(true);
        }
        out[0] = (2.0 + 9.0 * x[0] * x[0]) * v[0];
        true
    };
    let mut source =
        PhsPortSidebandSource::forced_analytic(&system, &input, &base, &c, &hessian).unwrap();
    assert!(matches!(
        source.solve_adaptive(carrier, &[force], &scales(1), &refinement, || observing
            .get()),
        Err(PhsPortSidebandError::Cancelled)
    ));
    assert!(observing.get());
    observing.set(false);
    let result = source
        .solve_adaptive(carrier, &[force], &scales(1), &refinement, || false)
        .unwrap()
        .response;
    let dq = force / C64::new(0.0, carrier);
    let modulation = (q1 * q1).scale(9.0);
    for (actual, expected) in result.port_coefficients.iter().zip([
        modulation.conj() * dq,
        C64::ZERO,
        dq.scale(2.0 + 4.5 * amplitude * amplitude),
        C64::ZERO,
        modulation * dq,
    ]) {
        near(*actual, expected, 2e-12);
    }

    // The accepted state order can outgrow the authored initial observer grid.
    // Full accepted coefficients remain present; grow only the minimum grid,
    // never the caller's maximum or its independent output order.
    let mut c = controls(1, 1, 15, 2);
    c.sideband.maximum_unknowns = 63;
    c.initial_samples = 32;
    c.max_samples = 128;
    let refinement = SidebandRefinementControls {
        maximum_harmonics: 31,
        equation_tolerances: vec![2e-10],
        state_tolerances: vec![2e-10],
    };
    let mut forcing = vec![C64::ZERO; 31];
    forcing[15] = force;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    assert!(matches!(
        source.solve_adaptive(carrier, &forcing, &scales(1), &refinement, || false),
        Err(PhsPortSidebandError::Budget(_))
    ));
    c.max_samples = 256;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    let result = source
        .solve_adaptive(carrier, &forcing, &scales(1), &refinement, || false)
        .unwrap()
        .response;
    assert_eq!(result.state_harmonics, 31);
    assert_eq!(result.state_response.coefficients.len(), 63);
    assert_eq!(result.output_harmonics, 2);
    assert_eq!(
        result
            .observation_trace
            .iter()
            .map(|check| check.samples)
            .collect::<Vec<_>>(),
        [64, 128, 256]
    );
    assert_eq!(result.hessian_actions, 896);
    near(result.state_response.coefficients[31], dq, 2e-12);
    near(result.port_coefficients[0], modulation.conj() * dq, 2e-12);
    near(result.port_coefficients[4], modulation * dq, 2e-12);
    assert_eq!(base.coeffs, original);
}

#[test]
fn carrier_derivative_preserves_dc_kinematics_source_phase_and_full_maxwell_dispersion() {
    let gain = 0.7;
    let feedback = VelocityFeedback(gain);
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: ARM,
            relaxation_time_s: TAU,
        }])
        .unwrap();
    let base_force = C64::new(0.04, -0.02);
    let base = orbit(
        OMEGA,
        exact(OMEGA, base_force, gain)
            .into_iter()
            .map(|z| vec![C64::ZERO, z])
            .collect(),
    );
    let original = base.coeffs.clone();
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO, base_force]]).unwrap();
    let calls = Cell::new(0_usize);
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out.fill(0.0);
        out[1] = -gain / MASS;
        true
    };
    let mut c = controls(3, 1, 1, 1);
    c.max_hessian_actions = 1792;
    let forcing = [
        C64::from_re(0.2),
        C64::new(0.1, -0.07),
        C64::new(0.03, 0.09),
    ];
    let slope = [
        C64::new(-0.03, 0.04),
        C64::new(0.02, -0.05),
        C64::new(0.01, 0.03),
    ];
    let mut source = PhsPortSidebandSource::forced_feedback_analytic(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        &c,
        &full_hessian,
    )
    .unwrap();
    let mut prepared = source.prepare(&scales(3), || false).unwrap();
    let preparations = calls.get();
    assert_eq!(preparations, 385);
    for carrier in [OMEGA, 1.2, -0.4] {
        let mut factor = prepared.factor(carrier, || false).unwrap();
        for authored_slope in [None, Some(slope.as_slice())] {
            let pair = factor
                .solve_with_carrier_derivative(
                    &forcing,
                    authored_slope,
                    derivative_controls(),
                    &[2e-10],
                    || false,
                )
                .unwrap();
            assert_eq!(pair.total_hessian_actions, 1792);
            assert_eq!(pair.response.hessian_actions, 896);
            assert_eq!(pair.derivative_observation_trace.len(), 3);
            assert_eq!(calls.get(), preparations);
            for (band, force) in forcing.into_iter().enumerate() {
                let w = carrier + (band as f64 - 1.0) * OMEGA;
                let expected = exact(w, force, gain);
                let expected_derivative = exact_derivative(
                    w,
                    force,
                    authored_slope.map_or(C64::ZERO, |s| s[band]),
                    gain,
                );
                for j in 0..3 {
                    near(
                        pair.response.state_response.coefficients[3 * band + j],
                        expected[j],
                        5e-12,
                    );
                    near(
                        pair.state_angular_derivative.coefficients[3 * band + j],
                        expected_derivative[j],
                        6e-12,
                    );
                }
                near(
                    pair.response.port_coefficients[band],
                    expected[1].scale(1.0 / MASS),
                    5e-12,
                );
                near(
                    pair.port_angular_derivative[band],
                    expected_derivative[1].scale(1.0 / MASS),
                    6e-12,
                );
            }
            if carrier == OMEGA {
                // The DC band has zero velocity, but its continuous carrier
                // derivative contains i*q even with a frequency-independent force.
                assert!(pair.response.port_coefficients[0].abs() < 1e-13);
                near(
                    pair.port_angular_derivative[0],
                    C64::new(0.0, forcing[0].re / STIFFNESS),
                    5e-12,
                );
            }
        }
    }
    assert_eq!(base.coeffs, original);
}

#[test]
fn carrier_derivative_has_independent_output_checks_and_atomic_combined_observation_budgets() {
    let system =
        PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(QuarticOnly)).unwrap();
    let (omega, amplitude, phase, carrier) = (1.1, 0.4, 0.2, 0.01);
    let q1 = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let base = orbit(omega, vec![vec![C64::ZERO, q1]]);
    let original = base.coeffs.clone();
    let input =
        PeriodicPortInput::new(omega, vec![vec![C64::ZERO, C64::new(0.0, omega) * q1]]).unwrap();
    let force = C64::from_re(1e-9);
    let actions = Cell::new(0_usize);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        actions.set(actions.get() + 1);
        out[0] = (2.0 + 9.0 * x[0] * x[0]) * v[0];
        true
    };
    let mut c = controls(1, 1, 0, 0);
    c.max_hessian_actions = 1792;
    c.envelope_tolerances[0] = 1e-6;
    {
        let mut source =
            PhsPortSidebandSource::forced_analytic(&system, &input, &base, &c, &hessian).unwrap();
        let mut prepared = source.prepare(&scales(1), || false).unwrap();
        let mut factor = prepared.factor(carrier, || false).unwrap();
        // Its state derivative is exact. The primal output tail fits its own
        // units/tolerance, but the omitted derivative output is 100 times larger.
        assert!(matches!(
            factor.solve_with_carrier_derivative(
                &[force],
                None,
                derivative_controls(),
                &[1e-8],
                || false,
            ),
            Err(PhsPortSidebandError::DerivativeObservationUnresolved { .. })
        ));
        factor.solve(&[force], || false).unwrap();
    }
    c.output_harmonics = 2;
    for limitation in 0..3 {
        let mut limited = c.clone();
        match limitation {
            0 => limited.max_hessian_actions = 1791,
            1 => limited.max_observation_work = 1_000_000,
            // Enough for one 256-point observation, but not the first returned
            // output/history retained alongside the derivative's final grid.
            _ => limited.max_workspace_bytes = 43_592,
        }
        let mut source =
            PhsPortSidebandSource::forced_analytic(&system, &input, &base, &limited, &hessian)
                .unwrap();
        let mut prepared = source.prepare(&scales(1), || false).unwrap();
        let mut factor = prepared.factor(carrier, || false).unwrap();
        let before = actions.get();
        assert!(matches!(
            factor.solve_with_carrier_derivative(
                &[force],
                None,
                derivative_controls(),
                &[1e-8],
                || false,
            ),
            Err(PhsPortSidebandError::Budget(_))
        ));
        if limitation < 2 {
            assert_eq!(actions.get(), before);
        }
        factor.solve(&[force], || false).unwrap();
    }
    let mut source =
        PhsPortSidebandSource::forced_analytic(&system, &input, &base, &c, &hessian).unwrap();
    let mut prepared = source.prepare(&scales(1), || false).unwrap();
    let mut factor = prepared.factor(carrier, || false).unwrap();
    let before = actions.get();
    assert!(matches!(
        factor.solve_with_carrier_derivative(
            &[force],
            None,
            derivative_controls(),
            &[1e-8],
            || actions.get() > before + 896,
        ),
        Err(PhsPortSidebandError::Cancelled)
    ));
    assert_eq!(actions.get(), before + 897);
    let pair = factor
        .solve_with_carrier_derivative(&[force], None, derivative_controls(), &[1e-8], || false)
        .unwrap();
    let q = force / C64::new(0.0, carrier);
    let dq = -q.scale(1.0 / carrier);
    let modulation = (q1 * q1).scale(9.0);
    near(pair.state_angular_derivative.coefficients[0], dq, 2e-12);
    for (actual, expected) in pair.port_angular_derivative.iter().zip([
        modulation.conj() * dq,
        C64::ZERO,
        dq.scale(2.0 + 4.5 * amplitude * amplitude),
        C64::ZERO,
        modulation * dq,
    ]) {
        near(*actual, expected, 2e-12);
    }
    assert_eq!(pair.total_hessian_actions, 1792);
    assert_eq!(base.coeffs, original);
}

#[test]
fn joint_adaptive_driven_feedback_preserves_original_maxwell_source_derivative_and_dc_velocity() {
    use fs_orbit::sideband::sensitivity::refinement::{
        SidebandCarrierRefinementControls, SidebandCarrierRefinementTolerances,
    };
    let gain = 0.7;
    let feedback = VelocityFeedback(gain);
    let system = fs_phs::mass_spring_damper(MASS, STIFFNESS, DAMPING)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: ARM,
            relaxation_time_s: TAU,
        }])
        .unwrap();
    let base_force = C64::new(0.04, -0.02);
    let base = orbit(
        OMEGA,
        exact(OMEGA, base_force, gain)
            .into_iter()
            .map(|z| vec![C64::ZERO, z])
            .collect(),
    );
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO, base_force]]).unwrap();
    let derivative = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        out[1] = -gain / MASS;
        true
    };
    let mut c = controls(3, 1, 1, 1);
    c.sideband.maximum_unknowns = 21;
    c.max_hessian_actions = 1792;
    let refinement = SidebandCarrierRefinementControls {
        maximum_harmonics: 3,
        tolerances: SidebandCarrierRefinementTolerances {
            equation_tolerances: vec![2e-10; 3],
            state_tolerances: vec![2e-10; 3],
            derivative_equation_tolerances: vec![2e-10; 3],
            derivative_state_tolerances: vec![2e-10; 3],
        },
    };
    let forcing = [
        C64::from_re(0.2),
        C64::new(0.1, -0.07),
        C64::new(0.03, 0.09),
    ];
    let slope = [
        C64::new(-0.03, 0.04),
        C64::new(0.02, -0.05),
        C64::new(0.01, 0.03),
    ];
    let mut source = PhsPortSidebandSource::forced_feedback_analytic(
        &system,
        &input,
        &feedback,
        &derivative,
        &base,
        &c,
        &full_hessian,
    )
    .unwrap();
    assert!(matches!(
        source.solve_adaptive_with_carrier_derivative(
            OMEGA,
            &forcing,
            Some(&slope),
            &scales(3),
            derivative_controls(),
            &[2e-10],
            &refinement,
            || true,
        ),
        Err(PhsPortSidebandError::Cancelled)
    ));
    for carrier in [OMEGA, -0.4] {
        let resolved = source
            .solve_adaptive_with_carrier_derivative(
                carrier,
                &forcing,
                Some(&slope),
                &scales(3),
                derivative_controls(),
                &[2e-10],
                &refinement,
                || false,
            )
            .unwrap();
        assert_eq!(
            resolved
                .refinement_trace
                .iter()
                .map(|c| c.response.harmonics)
                .collect::<Vec<_>>(),
            [1, 3]
        );
        assert!(resolved.refinement_trace.last().unwrap().is_resolved());
        let pair = resolved.response;
        assert_eq!(pair.response.state_harmonics, 3);
        assert_eq!(pair.response.output_harmonics, 1);
        assert_eq!(pair.total_hessian_actions, 1792);
        for k in -3_isize..=3 {
            let (force, slope) = if k.unsigned_abs() <= 1 {
                (forcing[(k + 1) as usize], slope[(k + 1) as usize])
            } else {
                (C64::ZERO, C64::ZERO)
            };
            let w = carrier + k as f64 * OMEGA;
            let expected = exact(w, force, gain);
            let expected_slope = exact_derivative(w, force, slope, gain);
            for j in 0..3 {
                near(
                    pair.response.state_response.coefficients[(k + 3) as usize * 3 + j],
                    expected[j],
                    6e-12,
                );
                near(
                    pair.state_angular_derivative.coefficients[(k + 3) as usize * 3 + j],
                    expected_slope[j],
                    8e-12,
                );
            }
            if k.unsigned_abs() <= 1 {
                near(
                    pair.response.port_coefficients[(k + 1) as usize],
                    expected[1].scale(1.0 / MASS),
                    6e-12,
                );
                near(
                    pair.port_angular_derivative[(k + 1) as usize],
                    expected_slope[1].scale(1.0 / MASS),
                    8e-12,
                );
            }
        }
    }
}

#[test]
fn joint_adaptive_refuses_hidden_output_derivative_and_combined_physical_actions() {
    use fs_orbit::sideband::sensitivity::refinement::{
        SidebandCarrierRefinementControls, SidebandCarrierRefinementTolerances,
    };
    let system =
        PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(QuarticOnly)).unwrap();
    let (omega, amplitude, phase, carrier) = (1.1, 0.4, 0.2, 0.01);
    let q1 = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * amplitude);
    let base = orbit(omega, vec![vec![C64::ZERO, q1]]);
    let input =
        PeriodicPortInput::new(omega, vec![vec![C64::ZERO, C64::new(0.0, omega) * q1]]).unwrap();
    let force = C64::from_re(1e-9);
    let refinement = SidebandCarrierRefinementControls {
        maximum_harmonics: 1,
        tolerances: SidebandCarrierRefinementTolerances {
            equation_tolerances: vec![2e-10],
            state_tolerances: vec![2e-10],
            derivative_equation_tolerances: vec![2e-10],
            derivative_state_tolerances: vec![2e-10],
        },
    };
    let mut c = controls(1, 1, 0, 0);
    c.sideband.maximum_unknowns = 3;
    c.max_hessian_actions = 1792;
    c.envelope_tolerances = vec![1e-6];
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    assert!(matches!(
        source.solve_adaptive_with_carrier_derivative(
            carrier,
            &[force],
            None,
            &scales(1),
            derivative_controls(),
            &[1e-8],
            &refinement,
            || false,
        ),
        Err(PhsPortSidebandError::DerivativeObservationUnresolved { .. })
    ));
    c.output_harmonics = 2;
    c.max_hessian_actions = 1791;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    assert!(matches!(
        source.solve_adaptive_with_carrier_derivative(
            carrier,
            &[force],
            None,
            &scales(1),
            derivative_controls(),
            &[1e-8],
            &refinement,
            || false,
        ),
        Err(PhsPortSidebandError::Budget(_))
    ));
    c.max_hessian_actions = 1792;
    let mut source = PhsPortSidebandSource::forced(&system, &input, &base, &c).unwrap();
    let pair = source
        .solve_adaptive_with_carrier_derivative(
            carrier,
            &[force],
            None,
            &scales(1),
            derivative_controls(),
            &[1e-8],
            &refinement,
            || false,
        )
        .unwrap()
        .response;
    let state = force / C64::new(0.0, carrier);
    let slope = state.scale(-1.0 / carrier);
    near(pair.state_angular_derivative.coefficients[1], slope, 2e-12);
    near(
        pair.port_angular_derivative[0],
        (q1.conj() * q1.conj()).scale(9.0) * slope,
        2e-12,
    );
    near(
        pair.port_angular_derivative[4],
        (q1 * q1).scale(9.0) * slope,
        2e-12,
    );
}
