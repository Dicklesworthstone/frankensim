//! Prescribed force ramps retain their actual velocity path, work and heat.
use super::*;

fn interpolate(first: Vec3, last: Vec3, fraction: f64) -> Vec3 {
    first.scale(1.0 - fraction).add(last.scale(fraction))
}

// Independent standard-lib closed ODE solution and Simpson power quadrature.
// Used at moderate exposure, where terminal-coefficient subtraction is safe.
fn analytic_velocity(
    initial: Vec3,
    ambient: Vec3,
    first: Vec3,
    last: Vec3,
    gravity: Vec3,
    mass: f64,
    resistance: f64,
    duration: f64,
    time: f64,
) -> Vec3 {
    let rate = resistance / mass;
    let slope = last.sub(first).scale(1.0 / (resistance * duration));
    let offset = first.add(gravity.scale(mass)).scale(1.0 / resistance)
        .sub(slope.scale(1.0 / rate));
    let amplitude = initial.sub(ambient).sub(offset);
    ambient.add(offset).add(slope.scale(time)).add(amplitude.scale((-rate * time).exp()))
}

fn simpson(mut f: impl FnMut(f64) -> f64, duration: f64) -> f64 {
    let intervals = 8192;
    let dt = duration / f64::from(intervals);
    let mut sum = f(0.0) + f(duration);
    let mut correction = 0.0;
    for i in 1..intervals {
        let term = f(f64::from(i) * dt) * if i % 2 == 0 { 2.0 } else { 4.0 };
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    sum * dt / 3.0
}

#[test]
fn affine_force_velocity_and_independent_power_integrals_cover_both_bases() {
    let p = properties();
    let ambient = Vec3::new(0.5, -1.0, 0.25);
    let gravity = Gravity::new(Vec3::new(0.2, 0.0, -0.3)).unwrap();
    let initial_v = Vec3::new(1.25, -0.5, 0.75);
    let initial = state(p, initial_v, Vec3::new(0.2, -0.3, 0.1));
    let first = Vec3::new(0.3, -1.2, 2.0);
    let last = Vec3::new(-1.0, 0.7, -0.2);
    let c = 4.0;
    let flow = IsotropicLinearDrag::new(p, gravity, ambient, c, 6.0).unwrap();
    for exposure in [0.01, 0.25, 1.999_999, 2.0, 2.000_001, 8.0] {
        let duration = exposure * p.mass() / c;
        let step = flow.step_affine_force(initial, force(first), last, duration).unwrap();
        let velocity = |time| analytic_velocity(
            initial_v, ambient, first, last, gravity.acceleration_world(), p.mass(), c,
            duration, time,
        );
        close_vector(
            step.state_after.center_of_mass_velocity_world(p).unwrap(),
            velocity(duration),
            3e-12,
        );
        let displacement = Vec3::new(
            simpson(|t| velocity(t).x, duration),
            simpson(|t| velocity(t).y, duration),
            simpson(|t| velocity(t).z, duration),
        );
        close_vector(step.displacement_world, displacement, 3e-12);
        let heat = simpson(|t| c * velocity(t).sub(ambient).norm_squared(), duration);
        let work = simpson(|t| interpolate(first, last, t / duration).dot(velocity(t)), duration);
        close(step.translation_heat_j, heat, 3e-12);
        close(step.applied_work_j, work, 3e-12);
        close(step.background_work_j, step.drag_impulse_world.dot(ambient), 2e-14);
        close(step.work_energy_residual_j, 0.0, 3e-13);
        close_vector(step.linear_momentum_residual_world, Vec3::ZERO, 2e-13);
        let free_spin = flow.step(initial, Wrench::ZERO, duration).unwrap();
        assert_eq!(step.state_after.pose().orientation(), free_spin.state_after.pose().orientation());
        assert_eq!(step.state_after.angular_momentum_body(), free_spin.state_after.angular_momentum_body());
        assert_eq!(step.rotation_heat_j, free_spin.rotation_heat_j);
    }
}

#[test]
fn zero_applied_impulse_still_transfers_work_through_the_ramped_velocity_path() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let initial = state(p, Vec3::ZERO, Vec3::ZERO);
    let flow = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 2.0, 0.0).unwrap();
    let step = flow.step_affine_force(
        initial, force(Vec3::new(1.0, 0.0, 0.0)), Vec3::new(-1.0, 0.0, 0.0), 1.0,
    ).unwrap();
    let e = (-2.0_f64).exp();
    close(step.state_after.linear_momentum_world().x, -e, 2e-14);
    close(step.displacement_world.x, 0.5 * e, 2e-14);
    close(step.applied_work_j, 1.0 / 6.0 - e, 2e-14);
    close(step.translation_heat_j, 1.0 / 6.0 - e - 0.5 * e * e, 2e-14);
    assert!(step.applied_work_j > 0.0);
    assert!(step.translation_heat_j > 0.0);
    close(step.work_energy_residual_j, 0.0, 2e-14);
}

#[test]
fn affine_path_and_accounts_compose_at_the_actual_shared_force_endpoints() {
    let p = properties();
    let flow = IsotropicLinearDrag::new(
        p, Gravity::new(Vec3::new(0.1, 0.0, -0.2)).unwrap(),
        Vec3::new(0.2, -0.4, 0.3), 4.0, 6.0,
    ).unwrap();
    let initial = state(p, Vec3::new(-0.5, 0.3, 1.0), Vec3::new(0.2, 0.4, -0.3));
    let first = Vec3::new(0.5, -0.7, 1.4);
    let last = Vec3::new(-0.4, 0.8, -0.1);
    let duration = 20.0;
    let whole = flow.step_affine_force(initial, force(first), last, duration).unwrap();
    let mut current = initial;
    let mut time = 0.0;
    let mut heat = 0.0;
    let mut work = 0.0;
    let mut bath = 0.0;
    let mut impulse = Vec3::ZERO;
    for dt in [0.08, 0.17, 1.0, 18.75] {
        let part = flow.step_affine_force(
            current, force(interpolate(first, last, time / duration)),
            interpolate(first, last, (time + dt) / duration), dt,
        ).unwrap();
        current = part.state_after;
        time += dt;
        heat += part.dissipated_energy_j;
        work += part.applied_work_j;
        bath += part.background_work_j;
        impulse = impulse.add(part.drag_impulse_world);
    }
    close_vector(current.linear_momentum_world(), whole.state_after.linear_momentum_world(), 2e-13);
    close_vector(current.pose().position_world(), whole.state_after.pose().position_world(), 2e-13);
    close_vector(impulse, whole.drag_impulse_world, 2e-13);
    close(heat, whole.dissipated_energy_j, 2e-13);
    close(work, whole.applied_work_j, 2e-13);
    close(bath, whole.background_work_j, 2e-13);
}

#[test]
fn tiny_drag_resolves_positive_reversal_heat_after_mean_and_linear_modes_cancel() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let c = 1e-20;
    let flow = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, c, 0.0).unwrap();
    let step = flow.step_affine_force(
        state(p, Vec3::new(1.0 / 12.0, 0.0, 0.0), Vec3::ZERO),
        force(Vec3::new(-0.5, 0.0, 0.0)), Vec3::new(0.5, 0.0, 0.0), 1.0,
    ).unwrap();
    // At zero exposure v(s)=(s²-s)/2+1/12, orthogonal to both 1 and s.
    // Its mean-square is 1/720; the first drag correction is O(c).
    assert!(step.translation_heat_j > 0.0);
    close(step.translation_heat_j / (c / 720.0), 1.0, 2e-14);
    close(step.displacement_world.x, 0.0, 2e-17);
}

#[test]
fn zero_drag_retains_affine_ballistics_work_and_independence_from_bath_velocity() {
    let p = properties();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -0.5)).unwrap();
    let initial_v = Vec3::new(0.5, -1.0, 0.25);
    let initial = state(p, initial_v, Vec3::new(0.2, -0.1, 0.4));
    let first = Vec3::new(2.0, -1.0, 3.0);
    let last = Vec3::new(-3.0, 2.0, -1.0);
    let duration = 1.5;
    let baseline = IsotropicLinearDrag::new(p, gravity, Vec3::ZERO, 0.0, 0.0).unwrap()
        .step_affine_force(initial, force(first), last, duration).unwrap();
    let shifted = IsotropicLinearDrag::new(
        p, gravity, Vec3::new(f64::MAX, -f64::MAX, f64::MAX), 0.0, 0.0,
    ).unwrap().step_affine_force(initial, force(first), last, duration).unwrap();
    assert_eq!(baseline, shifted);
    let acceleration = first.scale(1.0 / p.mass()).add(gravity.acceleration_world());
    let change = last.sub(first).scale(1.0 / p.mass());
    let expected_v = initial_v.add(acceleration.scale(duration)).add(change.scale(0.5 * duration));
    let expected_dx = initial_v.scale(duration)
        .add(acceleration.scale(0.5 * duration * duration))
        .add(change.scale(duration * duration / 6.0));
    close_vector(baseline.state_after.center_of_mass_velocity_world(p).unwrap(), expected_v, 2e-14);
    close_vector(baseline.displacement_world, expected_dx, 2e-14);
    let work = simpson(|t| {
        let v = initial_v.add(acceleration.scale(t)).add(change.scale(t * t / (2.0 * duration)));
        interpolate(first, last, t / duration).dot(v)
    }, duration);
    close(baseline.applied_work_j, work, 2e-14);
    assert_eq!(baseline.translation_heat_j, 0.0);
    assert_eq!(baseline.background_work_j, 0.0);
    close(baseline.work_energy_residual_j, 0.0, 2e-14);
}

#[test]
fn ramped_drive_obeys_galilean_work_and_heat_transformation() {
    let p = properties();
    let ambient = Vec3::new(0.5, -0.25, 0.125);
    let velocity = Vec3::new(-0.25, 0.75, 0.5);
    let boost = Vec3::new(32.0, -16.0, 8.0);
    let first = Vec3::new(0.3, -0.5, 0.9);
    let last = Vec3::new(-0.7, 0.2, 0.1);
    let duration = 3.0;
    let step = |shift: Vec3| {
        IsotropicLinearDrag::new(p, Gravity::ZERO, ambient.add(shift), 4.0, 0.0).unwrap()
            .step_affine_force(state(p, velocity.add(shift), Vec3::ZERO), force(first), last, duration).unwrap()
    };
    let base = step(Vec3::ZERO);
    let shifted = step(boost);
    close(shifted.translation_heat_j, base.translation_heat_j, 2e-14);
    close_vector(shifted.drag_impulse_world, base.drag_impulse_world, 2e-14);
    close_vector(shifted.displacement_world.sub(base.displacement_world), boost.scale(duration), 2e-14);
    close(shifted.applied_work_j - base.applied_work_j,
        boost.dot(first.add(last).scale(0.5 * duration)), 2e-14);
    close(shifted.background_work_j - base.background_work_j,
        boost.dot(base.drag_impulse_world), 2e-14);
}

#[test]
fn held_force_parity_and_invalid_affine_inputs() {
    let p = properties();
    let initial = state(p, Vec3::new(0.5, 0.25, -0.125), Vec3::new(0.2, -0.3, 0.4));
    let flow = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 4.0, 6.0).unwrap();
    let first = Vec3::new(0.5, -0.1, 0.3);
    assert_eq!(flow.step(initial, force(first), 1.0).unwrap(),
        flow.step_affine_force(initial, force(first), first, 1.0).unwrap());
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(flow.step_affine_force(initial, force(first), Vec3::new(invalid, 0.0, 0.0), 1.0).is_err());
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(flow.step_affine_force(initial, force(first), Vec3::ZERO, dt).is_err());
    }
    assert!(flow.step_affine_force(initial,
        Wrench { force_world: first, torque_body: Vec3::new(0.0, 0.0, 1.0) },
        Vec3::ZERO, 1.0).is_err());
}


#[test]
fn whole_interval_speed_bound_covers_an_interior_peak_with_zero_endpoint_slip() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let initial = state(p, Vec3::ZERO, Vec3::ZERO);
    let flow = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 1.0, 0.0).unwrap();
    let first = Vec3::new(1.0, 0.0, 0.0);
    let last = Vec3::new(2.0 - core::f64::consts::E, 0.0, 0.0);
    let step = flow.step_affine_force(initial, force(first), last, 1.0).unwrap();
    let bound = flow
        .affine_force_relative_speed_bound_m_s(initial, force(first), last, 1.0)
        .unwrap();
    // v(t)=e+(1-e)t-e*exp(-t): both endpoint slips vanish, but the
    // derivative vanishes at an interior positive velocity maximum.
    let peak_time = -(1.0 - (-1.0_f64).exp()).ln();
    let peak = 1.0 + (1.0 - core::f64::consts::E) * peak_time;
    close(step.state_after.linear_momentum_world().x, 0.0, 2e-14);
    close(
        analytic_velocity(
            Vec3::ZERO, Vec3::ZERO, first, last, Vec3::ZERO,
            1.0, 1.0, 1.0, peak_time,
        ).x,
        peak,
        2e-14,
    );
    assert!(peak > 0.2);
    assert!(peak < bound);
    close(bound, -(-1.0_f64).exp_m1(), 2e-14);
}

#[test]
fn speed_bound_retains_ballistic_limit_gravity_and_large_common_translation() {
    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(3.0, 3.0, 3.0)).unwrap();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -0.25)).unwrap();
    // Weight cancels the common z component: both acceleration norms are one.
    let first = Vec3::new(2.0, 0.0, 0.5);
    let last = Vec3::new(-2.0, 0.0, 0.5);
    for boost in [Vec3::ZERO, Vec3::new(2.0_f64.powi(40), -2.0_f64.powi(39), 2.0_f64.powi(38))] {
        let initial = state(p, boost.add(Vec3::new(0.0, 1.0, 0.0)), Vec3::ZERO);
        for resistance in [0.0, 1e-20] {
            let flow = IsotropicLinearDrag::new(p, gravity, boost, resistance, 0.0).unwrap();
            let bound = flow
                .affine_force_relative_speed_bound_m_s(initial, force(first), last, 0.01)
                .unwrap();
            close(bound, 1.01, 2e-14);
        }
        // A comparison trajectory that decreases retains the initial maximum.
        let flow = IsotropicLinearDrag::new(p, gravity, boost, 4.0, 0.0).unwrap();
        assert_eq!(
            flow.affine_force_relative_speed_bound_m_s(initial, force(first), last, 10.0).unwrap(),
            1.0,
        );
        assert!(flow.affine_force_relative_speed_bound_m_s(
            initial,
            Wrench { force_world: first, torque_body: Vec3::new(0.0, 0.0, 1.0) },
            last,
            1.0,
        ).is_err());
        assert!(flow.affine_force_relative_speed_bound_m_s(
            initial, force(first), Vec3::new(f64::NAN, 0.0, 0.0), 1.0,
        ).is_err());
        assert!(flow.affine_force_relative_speed_bound_m_s(
            initial, force(first), last, 0.0,
        ).is_err());
    }
}
