//! Original source laws, independent oscillator solutions and real event work.
#[path = "adhesive_dmt_coupled/common.rs"]
mod common;
use common::{Z, close, sources};
use fs_contact::normal_patch::dmt::elastic_events::*;
use fs_contact::normal_patch::dmt::events::{
    DmtSlipEventBudget, DmtSlipEventKind, DmtSlipRegime, MAX_SLIP_EVENTS,
};
use fs_exec::CancelGate;
use fs_math::det;

const M: f64 = 1e-9;
const K: f64 = 0.1;
const W: f64 = 1e4;
const A: f64 = 1e-7;
fn budget() -> DmtElasticEventBudget {
    DmtElasticEventBudget {
        motion: DmtSlipEventBudget {
            maximum_step_s: 0.1,
            minimum_indentation_m: -5.0 * Z,
            maximum_normal_increment_m: Z,
            maximum_accumulated_slip_m: 1e-3,
            maximum_tangent_speed_m_s: 1.0,
            maximum_force_n: 1e-2,
            maximum_events: MAX_SLIP_EVENTS,
            impulse_absolute_tolerance_n_s: 1e-24,
            impulse_relative_tolerance: 2e-11,
            energy_absolute_tolerance_j: 1e-26,
            energy_relative_tolerance: 2e-11,
        },
        maximum_spring_extension_m: 1e-3,
        maximum_stage_speed_m_s: 1.0,
        maximum_segments: 1024,
        maximum_root_iterations: 64,
        event_time_absolute_tolerance_s: 1e-16,
        maximum_numerical_work_uncertainty_j: 1e-22,
    }
}
fn check(a: f64, b: f64) {
    close(a, b, 2e-23, 3e-10);
}

#[test]
fn original_sourced_zero_shear_conservative_oscillator_has_no_spurious_loss() {
    let (n, s) = sources(1e-5, 0.0, 0.0);
    let owner = DmtElasticSlipEvents::new(&n, &s, M, K, 0.0, budget()).unwrap();
    let start = owner.initialize(0.5 * Z, A, 0.0, 0.0).unwrap();
    let t = 2.7 / W;
    let step = owner
        .advance(start, 0.0, 0.0, t, &CancelGate::new_clock_free())
        .unwrap();
    check(step.state().tangent_displacement_m(), A * det::cos(W * t));
    check(
        step.state().tangent_velocity_m_s(),
        -A * W * det::sin(W * t),
    );
    check(
        step.spring_energy_j + step.tangent_kinetic_energy_j,
        0.5 * K * A * A,
    );
    assert_eq!(step.friction_dissipation_j, 0.0);
    assert_eq!(step.damping_dissipation_j, 0.0);
    assert!(step.numerical_work_uncertainty_j > 0.0);
    assert_eq!(step.prescribed_normal_work_j, 0.0);
    check(step.reaction_impulse_n_s, M * A * W * det::sin(W * t));
    check(
        step.reaction_time_moment_n_s2,
        M * A * (W * t * det::sin(W * t) + det::cos(W * t) - 1.0),
    );
    assert!(std::ptr::eq(owner.normal(), &n) && std::ptr::eq(owner.shear(), &s));
}

#[test]
fn sliding_loading_spring_matches_independent_shifted_equilibrium_solution() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let f = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    for stiffness in [0.025, K, 0.4] {
        let owner = DmtElasticSlipEvents::new(&n, &s, M, stiffness, 0.0, budget()).unwrap();
        let omega = det::sqrt(stiffness / M);
        let v0 = 0.01;
        let t = 0.4 / omega;
        let initial = owner.initialize(0.5 * Z, 0.0, v0, 0.0).unwrap();
        let step = owner
            .advance(initial, 0.0, 0.0, t, &CancelGate::new_clock_free())
            .unwrap();
        let expected_x =
            v0 / omega * det::sin(omega * t) - f / stiffness * (1.0 - det::cos(omega * t));
        let expected_v = v0 * det::cos(omega * t) - f / (M * omega) * det::sin(omega * t);
        check(step.state().tangent_displacement_m(), expected_x);
        check(step.state().tangent_velocity_m_s(), expected_v);
        check(step.friction_dissipation_j, f * expected_x);
        assert_eq!(step.event_count, 0);
    }
}

#[test]
fn actual_stop_holds_then_moving_stage_restarts_at_original_static_threshold() {
    let (n, s) = sources(1e-5, 2e7, 1e7);
    let f = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let owner = DmtElasticSlipEvents::new(&n, &s, M, K, 0.0, budget()).unwrap();
    let v0 = 0.0005;
    let initial = owner.initialize(0.5 * Z, 0.0, v0, 0.0).unwrap();
    let stop_time = det::atan(M * W * v0 / f) / W;
    let stopped = owner
        .advance(
            initial,
            0.0,
            0.0,
            3.0 * stop_time,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(stopped.state().tangent_velocity_m_s(), 0.0);
    assert_eq!(stopped.state().regime(), DmtSlipRegime::Sticking);
    assert_eq!(stopped.event_count, 1);
    check(stopped.events[0].unwrap().offset_s, stop_time);
    let x = stopped.state().tangent_displacement_m();
    let u = 0.004;
    let release = (2.0 * f + K * x) / (K * u);
    let duration = 1.6 * release;
    let restart = owner
        .advance(
            stopped.state(),
            u * duration,
            0.0,
            duration,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        restart.events[0].unwrap().kind,
        DmtSlipEventKind::SlipStarted
    );
    check(restart.events[0].unwrap().offset_s, release);
    let h = duration - release;
    let a0 = f / M;
    let dx = u * (h - det::sin(W * h) / W) + a0 / (W * W) * (1.0 - det::cos(W * h));
    check(restart.tangent_increment_m, dx);
    check(
        restart.state().tangent_velocity_m_s(),
        u * (1.0 - det::cos(W * h)) + a0 / W * det::sin(W * h),
    );
    assert!(restart.stage_work_j > 0.0 && restart.friction_dissipation_j > 0.0);
    check(
        restart.spring_energy_j + restart.tangent_kinetic_energy_j + restart.friction_dissipation_j
            - restart.stage_work_j,
        stopped.spring_energy_j,
    );
}

#[test]
fn true_dashpot_integral_matches_independent_damped_solutions_in_all_regimes() {
    let (n, s) = sources(1e-5, 0.0, 0.0);
    for zeta in [0.2_f64, 1.0, 2.0, 10.0] {
        let c = 2.0 * M * W * zeta;
        let owner = DmtElasticSlipEvents::new(&n, &s, M, K, c, budget()).unwrap();
        let initial = owner.initialize(0.5 * Z, A, 0.0, 0.0).unwrap();
        let t = 0.7 / W;
        let step = owner
            .advance(initial, 0.0, 0.0, t, &CancelGate::new_clock_free())
            .unwrap();
        let (x, v, loss) = if zeta < 1.0 {
            let alpha = zeta * W;
            let beta = W * det::sqrt(1.0 - zeta * zeta);
            let decay = det::exp(-alpha * t);
            let x = A * decay * (det::cos(beta * t) + alpha / beta * det::sin(beta * t));
            let v = -A * decay * W * W / beta * det::sin(beta * t);
            let p = 2.0 * alpha;
            let q = 2.0 * beta;
            let cosine_integral = (det::exp(-p * t) * (-p * det::cos(q * t) + q * det::sin(q * t))
                + p)
                / (p * p + q * q);
            let integral = 0.5 * (-det::expm1(-p * t) / p - cosine_integral);
            (x, v, c * (A * W * W / beta).powi(2) * integral)
        } else if zeta == 1.0 {
            let decay = det::exp(-W * t);
            let integral = (2.0
                - det::exp(-2.0 * W * t) * ((2.0 * W * t).powi(2) + 4.0 * W * t + 2.0))
                / (2.0 * W).powi(3);
            (
                A * decay * (1.0 + W * t),
                -A * W * W * t * decay,
                c * A * A * W.powi(4) * integral,
            )
        } else {
            let a = -W * (zeta - det::sqrt(zeta * zeta - 1.0));
            let b = -W * (zeta + det::sqrt(zeta * zeta - 1.0));
            let coefficient = A * a * b / (b - a);
            let integral = det::expm1(2.0 * a * t) / (2.0 * a)
                - 2.0 * det::expm1((a + b) * t) / (a + b)
                + det::expm1(2.0 * b * t) / (2.0 * b);
            (
                A * (b * det::exp(a * t) - a * det::exp(b * t)) / (b - a),
                coefficient * (det::exp(a * t) - det::exp(b * t)),
                c * coefficient * coefficient * integral,
            )
        };
        check(step.state().tangent_displacement_m(), x);
        check(step.state().tangent_velocity_m_s(), v);
        check(step.damping_dissipation_j, loss);
        assert!(step.damping_dissipation_j > 0.0);
        assert_eq!(step.friction_dissipation_j, 0.0);
        check(
            step.spring_energy_j + step.tangent_kinetic_energy_j + step.damping_dissipation_j,
            0.5 * K * A * A,
        );
    }
}

#[test]
fn zero_net_slip_reversal_keeps_positive_actual_junction_loss() {
    let (n, s) = sources(1e-5, 1.5e7, 1e7);
    let f = 1e7 * n.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let owner = DmtElasticSlipEvents::new(&n, &s, M, K, 0.0, budget()).unwrap();
    let v0 = 0.01;
    let turn = det::atan(M * W * v0 / f) / W;
    let peak = det::sqrt((v0 / W).powi(2) + (f / K).powi(2)) - f / K;
    let returning = det::acos(-f / (K * peak - f)) / W;
    let initial = owner.initialize(0.5 * Z, 0.0, v0, 0.0).unwrap();
    let step = owner
        .advance(
            initial,
            0.0,
            0.0,
            turn + returning,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(step.tangent_increment_m, 0.0, 3e-18, 0.0);
    check(step.accumulated_slip_m, 2.0 * peak);
    check(step.friction_dissipation_j, 2.0 * f * peak);
    assert!(step.state().tangent_velocity_m_s() < 0.0);
    assert_eq!(step.event_count, 2);
    check(step.events[0].unwrap().offset_s, turn);
    check(
        step.tangent_kinetic_energy_j + step.friction_dissipation_j,
        0.5 * M * v0 * v0,
    );
}

#[test]
fn phase_resolved_underdamped_reversals_and_work_are_partition_invariant() {
    let (n, s) = sources(1e-5, 1.5e5, 1e5);
    let owner = DmtElasticSlipEvents::new(&n, &s, M, K, 0.08 * M * W, budget()).unwrap();
    let initial = owner.initialize(0.5 * Z, A, 0.0, 0.0).unwrap();
    let duration = 9.0 / W;
    let u = 0.0001;
    let applied = 1e-10;
    let gate = CancelGate::new_clock_free();
    let full = owner
        .advance(initial, u * duration, applied, duration, &gate)
        .unwrap();
    assert!(full.event_count >= 5);
    for count in [3, 17, 43] {
        let mut state = initial;
        let (mut loss, mut damping, mut work, mut impulse, mut moment, mut slip) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for i in 0..count {
            let step = owner
                .advance(
                    state,
                    u * duration * (i + 1) as f64 / count as f64,
                    applied,
                    duration / count as f64,
                    &gate,
                )
                .unwrap();
            loss += step.friction_dissipation_j;
            damping += step.damping_dissipation_j;
            work += step.stage_work_j;
            moment += step.reaction_time_moment_n_s2 + state.time_s() * step.reaction_impulse_n_s;
            impulse += step.reaction_impulse_n_s;
            slip += step.accumulated_slip_m;
            state = step.state();
        }
        check(
            state.tangent_displacement_m(),
            full.state().tangent_displacement_m(),
        );
        check(
            state.tangent_velocity_m_s(),
            full.state().tangent_velocity_m_s(),
        );
        check(loss, full.friction_dissipation_j);
        check(damping, full.damping_dissipation_j);
        check(work, full.stage_work_j);
        check(impulse, full.reaction_impulse_n_s);
        check(moment, full.reaction_time_moment_n_s2);
        check(slip, full.accumulated_slip_m);
    }
}

#[test]
fn source_cancel_phase_event_and_uncertainty_refusals_preserve_checkpoint_and_retry() {
    let (n, s) = sources(1e-5, 1.5e5, 1e5);
    let owner = DmtElasticSlipEvents::new(&n, &s, M, K, 0.1 * M * W, budget()).unwrap();
    let initial = owner.initialize(0.5 * Z, A, 0.0, 0.0).unwrap();
    let gate = CancelGate::new_clock_free();
    let expected = owner.advance(initial, 0.0, 0.0, 7.0 / W, &gate).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        owner.advance(initial, 0.0, 0.0, 7.0 / W, &cancelled),
        Err(DmtElasticEventError::Cancelled)
    ));
    let alien = DmtElasticSlipEvents::new(&n, &s, 2.0 * M, K, 0.0, budget()).unwrap();
    assert!(matches!(
        alien.advance(initial, 0.0, 0.0, 7.0 / W, &gate),
        Err(DmtElasticEventError::StateMismatch)
    ));
    for variant in 0..4 {
        let mut b = budget();
        match variant {
            0 => b.maximum_segments = 1,
            1 => b.motion.maximum_events = 1,
            2 => b.maximum_root_iterations = 1,
            _ => b.maximum_numerical_work_uncertainty_j = 1e-50,
        }
        let limited = DmtElasticSlipEvents::new(&n, &s, M, K, 0.1 * M * W, b).unwrap();
        assert!(
            limited.advance(initial, 0.0, 0.0, 7.0 / W, &gate).is_err(),
            "variant={variant}"
        );
    }
    assert_eq!(
        owner.advance(initial, 0.0, 0.0, 7.0 / W, &gate).unwrap(),
        expected
    );
}
