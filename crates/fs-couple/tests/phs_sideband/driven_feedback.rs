use super::{QuarticOnly, controls, near, orbit, scales};
use core::cell::Cell;
use fs_couple::phs_orbit::{
    AutonomousPortFeedback, PeriodicPortInput,
    sideband::{
        PhsPortSidebandError, PhsPortSidebandSource, solve_forced_feedback_port_sideband_response,
        solve_forced_feedback_port_sideband_response_analytic,
    },
};
use fs_math::{c64::C64, det};
use fs_orbit::sideband::SidebandError;
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
