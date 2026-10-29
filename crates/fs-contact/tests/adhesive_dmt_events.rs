//! G1/G3/G4: actual stopping/restart paths, independent closed forms and work.
#[path = "adhesive_dmt_coupled/common.rs"]
mod common;

use common::{Z, close, sources};
use fs_contact::normal_patch::dmt::events::*;
use fs_exec::CancelGate;

const MASS: f64 = 1e-9;
const T: f64 = 1e-4;

fn budget() -> DmtSlipEventBudget {
    DmtSlipEventBudget {
        maximum_step_s: 0.01,
        minimum_indentation_m: -5.0 * Z,
        maximum_normal_increment_m: 3.0 * Z,
        maximum_accumulated_slip_m: 1e-3,
        maximum_tangent_speed_m_s: 1.0,
        maximum_force_n: 1e-2,
        maximum_events: MAX_SLIP_EVENTS,
        impulse_absolute_tolerance_n_s: 1e-25,
        impulse_relative_tolerance: 2e-12,
        energy_absolute_tolerance_j: 1e-27,
        energy_relative_tolerance: 2e-12,
    }
}
fn check(a: f64, b: f64) {
    close(a, b, 1e-25, 3e-12);
}

#[test]
fn g1_exact_stop_uses_actual_slip_and_supports_the_remaining_time() {
    let (normal, shear) = sources(1e-5, 2e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let area = normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let force = 1e7 * area;
    let stop = 0.25 * T;
    let initial_velocity = force * stop / MASS;
    let initial = owner.initialize(0.5 * Z, 0.0, initial_velocity).unwrap();
    let step = owner
        .advance(initial, 0.5 * Z, 0.0, T, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(step.state().tangent_velocity_m_s(), 0.0);
    assert_eq!(step.state().regime(), DmtSlipRegime::Sticking);
    assert_eq!(step.event_count, 1);
    let event = step.events[0].unwrap();
    assert_eq!(event.kind, DmtSlipEventKind::Stopped);
    check(event.offset_s, stop);
    check(step.tangent_increment_m, 0.5 * initial_velocity * stop);
    check(step.accumulated_slip_m, step.tangent_increment_m);
    check(
        step.friction_dissipation_j,
        0.5 * MASS * initial_velocity.powi(2),
    );
    check(step.friction_impulse_n_s, MASS * initial_velocity);
    check(step.friction_time_moment_n_s2, 0.5 * force * stop.powi(2));
    assert_eq!(step.prescribed_normal_work_j, 0.0);
    assert!(core::ptr::eq(owner.normal(), &normal));
    assert!(core::ptr::eq(owner.shear(), &shear));
}

#[test]
fn g1_zero_net_displacement_reversal_has_positive_slip_and_dissipation() {
    let (normal, shear) = sources(1e-5, 1.5e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let k = 1e7 * normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let stop = 1e-5;
    let end = (1.0 + 3.0_f64.sqrt()) * stop;
    let v0 = 3.0 * k * stop / MASS;
    let initial = owner.initialize(0.5 * Z, 0.0, v0).unwrap();
    let step = owner
        .advance(
            initial,
            0.5 * Z,
            -2.0 * k,
            end,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(step.tangent_increment_m, 0.0, 2e-24, 0.0);
    check(step.accumulated_slip_m, v0 * stop);
    check(step.friction_dissipation_j, k * v0 * stop);
    check(
        step.state().tangent_velocity_m_s(),
        -k * (end - stop) / MASS,
    );
    assert_eq!(step.event_count, 2);
    assert_eq!(step.events[0].unwrap().kind, DmtSlipEventKind::Stopped);
    assert_eq!(step.events[1].unwrap().kind, DmtSlipEventKind::SlipStarted);
    check(step.events[0].unwrap().offset_s, stop);
    check(step.events[1].unwrap().offset_s, stop);
    assert!(step.friction_dissipation_j > 0.0);
}

#[test]
fn g1_declining_real_area_releases_static_junction_at_the_actual_threshold() {
    let (normal, shear) = sources(1e-5, 2e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let x0 = 0.8 * Z;
    let k0 = 1e7 * normal.evaluate(x0).unwrap().contact_area_m2;
    let applied = 1.5 * k0;
    let state = owner.initialize(x0, 0.0, 0.0).unwrap();
    let step = owner
        .advance(state, x0 * 0.25, applied, T, &CancelGate::new_clock_free())
        .unwrap();
    let release = T / 3.0;
    assert_eq!(step.event_count, 1);
    assert_eq!(step.events[0].unwrap().kind, DmtSlipEventKind::SlipStarted);
    check(step.events[0].unwrap().offset_s, release);
    let h = T - release;
    let a0 = 0.75 * k0 / MASS;
    let a1 = 0.75 * k0 / (MASS * T);
    check(
        step.state().tangent_velocity_m_s(),
        a0 * h + 0.5 * a1 * h.powi(2),
    );
    check(
        step.tangent_increment_m,
        0.5 * a0 * h.powi(2) + a1 * h.powi(3) / 6.0,
    );
    check(
        step.prescribed_normal_work_j,
        normal.reversible_work_between(x0, 0.25 * x0).unwrap(),
    );
    assert!(step.prescribed_normal_work_j != 0.0);
    check(
        step.applied_tangent_work_j,
        step.tangent_kinetic_energy_j + step.friction_dissipation_j,
    );
}

#[test]
fn g1_contact_entry_accumulates_no_gap_friction_and_stops_on_growing_area() {
    let (normal, shear) = sources(1e-5, 2e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let final_k = 1e7 * normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let rate = final_k / (0.5 * T);
    let stop_after_entry = T * 0.25;
    let v0 = 0.5 * rate * stop_after_entry.powi(2) / MASS;
    let initial = owner.initialize(-0.5 * Z, 0.0, v0).unwrap();
    let step = owner
        .advance(initial, 0.5 * Z, 0.0, T, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(step.initial_contact_area_m2, 0.0);
    assert_eq!(step.event_count, 2);
    assert_eq!(
        step.events[0].unwrap().kind,
        DmtSlipEventKind::ContactEntered
    );
    check(step.events[0].unwrap().offset_s, 0.5 * T);
    assert_eq!(step.events[1].unwrap().kind, DmtSlipEventKind::Stopped);
    check(step.events[1].unwrap().offset_s, 0.75 * T);
    check(
        step.tangent_increment_m,
        v0 * (0.5 * T + (2.0 / 3.0) * stop_after_entry),
    );
    check(step.friction_dissipation_j, 0.5 * MASS * v0.powi(2));
    assert_eq!(step.state().tangent_velocity_m_s(), 0.0);
}

#[test]
fn g3_event_partition_preserves_physical_state_total_slip_and_time_moment() {
    let (normal, shear) = sources(1e-5, 1.4e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let start = 0.8 * Z;
    let end = -0.4 * Z;
    let k = 1e7 * normal.evaluate(start).unwrap().contact_area_m2;
    let force = -1.7 * k;
    let initial = owner.initialize(start, 0.0, 0.4 * k * T / MASS).unwrap();
    let gate = CancelGate::new_clock_free();
    let whole = owner.advance(initial, end, force, T, &gate).unwrap();
    for n in [2, 7, 31] {
        let mut state = initial;
        let mut impulse = 0.0;
        let mut moment = 0.0;
        let mut slip = 0.0;
        let mut dissipation = 0.0;
        for i in 0..n {
            let step = owner
                .advance(
                    state,
                    start + (end - start) * f64::from(i + 1) / f64::from(n),
                    force,
                    T / f64::from(n),
                    &gate,
                )
                .unwrap();
            moment += step.friction_time_moment_n_s2 + state.time_s() * step.friction_impulse_n_s;
            impulse += step.friction_impulse_n_s;
            slip += step.accumulated_slip_m;
            dissipation += step.friction_dissipation_j;
            state = step.state();
        }
        check(
            state.tangent_velocity_m_s(),
            whole.state().tangent_velocity_m_s(),
        );
        check(state.tangent_displacement_m(), whole.tangent_increment_m);
        check(impulse, whole.friction_impulse_n_s);
        check(moment, whole.friction_time_moment_n_s2);
        check(slip, whole.accumulated_slip_m);
        check(dissipation, whole.friction_dissipation_j);
    }
}

#[test]
fn g1_two_finite_bodies_receive_reciprocal_impulse_and_position_moment() {
    let (normal, shear) = sources(1e-5, 1.4e7, 1e7);
    let (ma, mb) = (2e-9, 3e-9);
    let mass = ma * mb / (ma + mb);
    let owner = DmtSlipEvents::new(&normal, &shear, mass, budget()).unwrap();
    let k = 1e7 * normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let (fa, fb) = (1.1 * k, -2.2 * k);
    let relative_force = mass * (fb / mb - fa / ma);
    let relative_velocity = 0.25 * k * T / mass;
    let common_velocity = 0.01;
    let va = common_velocity - mb / (ma + mb) * relative_velocity;
    let vb = common_velocity + ma / (ma + mb) * relative_velocity;
    let initial = owner.initialize(0.5 * Z, 0.0, relative_velocity).unwrap();
    let step = owner
        .advance(
            initial,
            0.5 * Z,
            relative_force,
            T,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let impulse = step.friction_impulse_n_s;
    let position_moment = T * impulse - step.friction_time_moment_n_s2;
    let next_va = va + (fa * T + impulse) / ma;
    let next_vb = vb + (fb * T - impulse) / mb;
    let xa = va * T + (0.5 * fa * T.powi(2) + position_moment) / ma;
    let xb = vb * T + (0.5 * fb * T.powi(2) - position_moment) / mb;
    check(next_vb - next_va, step.state().tangent_velocity_m_s());
    check(xb - xa, step.tangent_increment_m);
    check(ma * (next_va - va) + mb * (next_vb - vb), (fa + fb) * T);
    let dk =
        0.5 * ma * (next_va - va) * (next_va + va) + 0.5 * mb * (next_vb - vb) * (next_vb + vb);
    check(dk + step.friction_dissipation_j, fa * xa + fb * xb);
}

#[test]
fn g0_g4_source_event_work_speed_and_cancellation_refuse_without_losing_retry() {
    let (normal, shear) = sources(1e-5, 1.5e7, 1e7);
    let (other_normal, other_shear) = sources(2e-5, 1.5e7, 1e7);
    assert!(matches!(
        DmtSlipEvents::new(&normal, &other_shear, MASS, budget()),
        Err(DmtSlipEventError::StateMismatch)
    ));
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let initial = owner.initialize(0.5 * Z, 0.0, 1e-4).unwrap();
    let gate = CancelGate::new_clock_free();
    let k = 1e7 * normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let expected = owner.advance(initial, 0.5 * Z, -2.0 * k, T, &gate).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        owner.advance(initial, 0.5 * Z, -2.0 * k, T, &cancelled),
        Err(DmtSlipEventError::Cancelled)
    );
    let alien = DmtSlipEvents::new(&other_normal, &other_shear, MASS, budget())
        .unwrap()
        .initialize(0.5 * Z, 0.0, 1e-4)
        .unwrap();
    assert_eq!(
        owner.advance(alien, 0.5 * Z, -2.0 * k, T, &gate),
        Err(DmtSlipEventError::StateMismatch)
    );
    let mut low = budget();
    low.maximum_events = 1;
    let bounded = DmtSlipEvents::new(&normal, &shear, MASS, low).unwrap();
    assert_eq!(
        bounded.advance(initial, 0.5 * Z, -2.0 * k, T, &gate),
        Err(DmtSlipEventError::EventBudget)
    );
    low = budget();
    low.maximum_accumulated_slip_m = expected.accumulated_slip_m * 0.5;
    let bounded = DmtSlipEvents::new(&normal, &shear, MASS, low).unwrap();
    assert!(matches!(
        bounded.advance(initial, 0.5 * Z, -2.0 * k, T, &gate),
        Err(DmtSlipEventError::Range("accumulated slip"))
    ));
    // End speeds are equal; the true maximum lies inside this accelerating-
    // then-decelerating interval. Endpoint-only checking would admit it.
    low = budget();
    low.maximum_tangent_speed_m_s = 0.0011;
    let bounded = DmtSlipEvents::new(&normal, &shear, MASS, low).unwrap();
    let state = bounded.initialize(0.25 * Z, 0.0, 0.001).unwrap();
    let force = 1e7 * normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    assert_eq!(
        bounded.advance(state, 0.75 * Z, force, T, &gate),
        Err(DmtSlipEventError::Range("tangent speed"))
    );
    assert_eq!(
        owner.advance(initial, 0.5 * Z, -2.0 * k, T, &gate).unwrap(),
        expected
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.steps(), 0);
}

#[test]
fn g1_equal_strength_threshold_has_no_artificial_capture_deadband() {
    let (normal, shear) = sources(1e-5, 1e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let k = 1e7 * normal.evaluate(0.5 * Z).unwrap().contact_area_m2;
    let state = owner.initialize(0.5 * Z, 0.0, 0.0).unwrap();
    let gate = CancelGate::new_clock_free();
    let held = owner.advance(state, 0.5 * Z, k, T, &gate).unwrap();
    assert_eq!(held.state().tangent_velocity_m_s(), 0.0);
    assert_eq!(held.event_count, 0);
    let released = owner.advance(state, 0.25 * Z, k, T, &gate).unwrap();
    assert_eq!(
        released.events[0].unwrap().kind,
        DmtSlipEventKind::SlipStarted
    );
    assert_eq!(released.events[0].unwrap().offset_s, 0.0);
    check(released.state().tangent_velocity_m_s(), 0.25 * k * T / MASS);
    check(released.tangent_increment_m, k * T.powi(2) / (12.0 * MASS));
}

#[test]
fn g0_original_area_rounding_cannot_trigger_false_static_release() {
    let (normal, shear) = sources(1e-5, 2e7, 1e7);
    let owner = DmtSlipEvents::new(&normal, &shear, MASS, budget()).unwrap();
    let gate = CancelGate::new_clock_free();
    for i in 1..32 {
        let indentation = Z * f64::from(i) / 32.0;
        let area = normal.evaluate(indentation).unwrap().contact_area_m2;
        let force = shear.law().static_shear_strength_pa() * area;
        let initial = owner.initialize(indentation, 0.0, 0.0).unwrap();
        let step = owner
            .advance(initial, indentation, force, T, &gate)
            .unwrap();
        assert_eq!(step.event_count, 0);
        assert_eq!(step.state().tangent_velocity_m_s(), 0.0);
        assert_eq!(step.tangent_increment_m, 0.0);
        assert_eq!(step.friction_dissipation_j, 0.0);
        check(step.friction_impulse_n_s, force * T);
    }
}
