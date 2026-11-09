//! G1/G3: exact native cycle power, active coefficients and parameter work.
use fs_mbd::rotordynamics::periodic::power::{
    RotorCyclePowerError, prepare_forced_power_cycle,
};
use fs_mbd::rotordynamics::periodic::{RotorFloquetBudget, RotorPeriodicSegment};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::f64::consts::PI;

fn model(m: f64, c: f64, k: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[m], &[c], &[k], &[0.0], || false).unwrap()
}
fn state(q: f64, v: f64) -> RotorTransientState {
    RotorTransientState {
        displacement: vec![q],
        velocity: vec![v],
    }
}
fn near(value: f64, expected: f64) {
    assert!(
        (value - expected).abs() < 2e-10 * (1.0 + expected.abs()),
        "{value:e} != {expected:e}"
    );
}

#[test]
fn shared_signed_quadratic_integral_matches_decay_and_refuses_caps() {
    use fs_phs::linear_exact::{DenseExponentialBudget, quadratic_flow_integral_with_budget};
    let budget = DenseExponentialBudget {
        max_dimension: 1,
        max_dense_bytes: 24 * 8,
        max_multiply_terms: 256,
    };
    // x'=-2x, observed signed power -3x². This is active work, not a
    // nonnegative damping law; the shared kernel must retain its sign.
    let result = quadratic_flow_integral_with_budget(&[-2.0], &[-3.0], 1, 0.7, budget).unwrap();
    near(result.transition[0], (-1.4_f64).exp());
    near(result.integral[0], 0.75 * (-2.8_f64).exp_m1());
    for limited in [
        DenseExponentialBudget { max_dimension: 0, ..budget },
        DenseExponentialBudget { max_dense_bytes: 24 * 8 - 1, ..budget },
        DenseExponentialBudget { max_multiply_terms: 255, ..budget },
    ] {
        assert!(quadratic_flow_integral_with_budget(&[-2.0], &[-3.0], 1, 0.7, limited).is_err());
    }
    assert!(quadratic_flow_integral_with_budget(&[-2.0], &[f64::NAN], 1, 0.7, budget).is_err());
    assert!(quadratic_flow_integral_with_budget(&[-2.0], &[-3.0], 1, -0.7, budget).is_err());
    assert!(quadratic_flow_integral_with_budget(&[-2.0], &[], 1, 0.7, budget).is_err());
}

#[test]
fn exact_velocity_reversal_retains_force_work_and_damper_heat() {
    let m = model(2.0, 3.0, 0.0);
    let h = 0.7_f64;
    let plan = prepare_forced_power_cycle(
        0.0,
        &[RotorPeriodicSegment { model: &m, duration_seconds: h }],
        &[&[5.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    ).unwrap();
    let (motion, power) = plan.advance(&state(0.0, -2.0), 1e-11, 1e-10, || false).unwrap();
    let equilibrium = 5.0 / 3.0;
    let difference = -2.0 - equilibrium;
    let rate = 1.5;
    let e = (-rate * h).exp();
    let travel = equilibrium * h + difference * (1.0 - e) / rate;
    let final_v = equilibrium + difference * e;
    let heat = 3.0 * (
        equilibrium * equilibrium * h
            + 2.0 * equilibrium * difference * (1.0 - e) / rate
            + difference * difference * (1.0 - e * e) / (2.0 * rate)
    );
    assert!(motion.boundary_states[1].velocity[0] > 0.0);
    near(motion.boundary_states[1].displacement[0], travel);
    near(power.applied_work_j, 5.0 * travel);
    near(power.velocity_force_work_j, heat);
    near(power.kinetic_energy_change_j, final_v * final_v - 4.0);
    near(power.average_applied_power_w, 5.0 * travel / h);
    assert!(power.balance_residual_j.abs() < 1e-12);
    // Neither midpoint velocity nor endpoint-power trapezoids integrate this
    // exponentially changing trajectory exactly.
    assert!((5.0 * h * (-2.0 + final_v) / 2.0 - power.applied_work_j).abs() > 0.1);
}

#[test]
fn parametric_resonance_and_period_seam_have_real_signed_work() {
    let soft = model(1.0, 0.0, 1.0);
    let stiff = model(1.0, 0.0, 4.0);
    let plan = prepare_forced_power_cycle(
        7.0,
        &[
            RotorPeriodicSegment { model: &soft, duration_seconds: PI / 2.0 },
            RotorPeriodicSegment { model: &stiff, duration_seconds: PI / 4.0 },
        ],
        &[&[0.0], &[0.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    ).unwrap();
    let (_, pumped) = plan.advance(&state(0.0, 1.0), 1e-11, 1e-10, || false).unwrap();
    near(pumped.stiffness_switch_work_j, 1.5);
    near(pumped.kinetic_energy_change_j, 1.5);
    near(pumped.applied_work_j, 0.0);
    near(pumped.velocity_force_work_j, 0.0);
    assert!(plan.motion_cycle().floquet_report().spectral_radius > 1.9);
    let (_, extracted) = plan.advance(&state(1.0, 0.0), 1e-11, 1e-10, || false).unwrap();
    near(extracted.segments[0].stiffness_switch_work_j, 0.0);
    near(extracted.segments[1].stiffness_switch_work_j, -0.375);
    near(extracted.elastic_energy_change_j, -0.375);
    near(extracted.stiffness_switch_work_j, -0.375);
}

#[test]
fn active_damping_and_circulatory_stiffness_remain_signed_work() {
    let active = model(1.0, -2.0, 0.0);
    let h = 0.2_f64;
    let plan = prepare_forced_power_cycle(
        0.0,
        &[RotorPeriodicSegment { model: &active, duration_seconds: h }],
        &[&[0.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    ).unwrap();
    let (_, power) = plan.advance(&state(0.0, 1.0), 1e-11, 1e-10, || false).unwrap();
    near(power.velocity_force_work_j, -(4.0 * h).exp_m1() / 2.0);
    assert!(power.velocity_force_work_j < 0.0);
    near(power.velocity_force_work_j + power.kinetic_energy_change_j, 0.0);

    // K[0,1]=2 and q_y=1 are constant, hence q_x=-t² and v_x=-2t.
    // Symmetric-K storage falls by t²; skew-K work is -t², while K grows 2t².
    let circulatory = RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[0.0, 2.0, 0.0, 0.0],
        &[0.0; 4],
        || false,
    ).unwrap();
    let plan = prepare_forced_power_cycle(
        0.0,
        &[RotorPeriodicSegment { model: &circulatory, duration_seconds: h }],
        &[&[0.0, 0.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    ).unwrap();
    let initial = RotorTransientState {
        displacement: vec![0.0, 1.0],
        velocity: vec![0.0, 0.0],
    };
    let (_, power) = plan.advance(&initial, 1e-11, 1e-10, || false).unwrap();
    near(power.kinetic_energy_change_j, 2.0 * h * h);
    near(power.elastic_energy_change_j, -h * h);
    near(power.nonconservative_stiffness_work_j, -h * h);
    near(power.velocity_force_work_j, 0.0);
}

#[test]
fn periodic_square_wave_work_closes_and_true_segment_refinement_is_invariant() {
    let m = model(1.0, 0.2, 1.0);
    let coarse = prepare_forced_power_cycle(
        0.0,
        &[RotorPeriodicSegment { model: &m, duration_seconds: PI }; 2],
        &[&[1.0], &[-1.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    ).unwrap();
    let (orbit, power) = coarse.periodic_orbit(1e-10, 1e8, 1e-8, || false).unwrap();
    let initial = &orbit.trajectory.boundary_states[0];
    // Held-force integration has an independent endpoint primitive.
    let expected_work = 2.0 * (
        orbit.trajectory.boundary_states[1].displacement[0] - initial.displacement[0]
    );
    near(power.applied_work_j, expected_work);
    near(power.velocity_force_work_j, expected_work);
    assert!(power.average_applied_power_w > 1.0);
    let refined = prepare_forced_power_cycle(
        0.0,
        &[RotorPeriodicSegment { model: &m, duration_seconds: PI / 2.0 }; 4],
        &[&[1.0], &[1.0], &[-1.0], &[-1.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    ).unwrap();
    let (_, subdivided) = refined.advance(initial, 1e-10, 1e-8, || false).unwrap();
    near(subdivided.applied_work_j, power.applied_work_j);
    near(subdivided.velocity_force_work_j, power.velocity_force_work_j);
}

#[test]
fn budget_energy_and_cancellation_refuse_without_mutating_state() {
    let m = model(1.0, 0.3, 1.0);
    let segments = [RotorPeriodicSegment { model: &m, duration_seconds: 0.4 }];
    assert!(prepare_forced_power_cycle(
        0.0,
        &segments,
        &[&[1.0]],
        1e-12,
        RotorFloquetBudget { max_dense_bytes: 0, ..RotorFloquetBudget::default() },
        || false,
    ).is_err());
    assert!(matches!(
        prepare_forced_power_cycle(
            0.0, &segments, &[&[1.0]], 1e-12, RotorFloquetBudget::default(), || true,
        ),
        Err(RotorError::Cancelled)
    ));
    let plan = prepare_forced_power_cycle(
        0.0, &segments, &[&[1.0]], 1e-12, RotorFloquetBudget::default(), || false,
    ).unwrap();
    let initial = state(0.7, -0.2);
    let saved = initial.clone();
    assert!(plan.advance(&initial, 1e-10, f64::NAN, || false).is_err());
    assert!(matches!(
        plan.advance(&initial, 1e-10, 1e-8, || true),
        Err(RotorCyclePowerError::Dynamics(RotorError::Cancelled))
    ));
    assert_eq!(initial, saved);
    assert_eq!(
        plan.advance(&initial, 1e-10, 1e-8, || false).unwrap(),
        plan.advance(&initial, 1e-10, 1e-8, || false).unwrap()
    );
}
