//! Moving baths retain their actual slip path and independent source work.
use super::*;

fn interpolate(first: Vec3, last: Vec3, fraction: f64) -> Vec3 {
    first.scale(1.0 - fraction).add(last.scale(fraction))
}

// Independent absolute-velocity ODE solution, used only at moderate exposure.
fn velocity(
    initial: Vec3,
    ambient: Vec3,
    ambient_after: Vec3,
    first: Vec3,
    last: Vec3,
    gravity: Vec3,
    mass: f64,
    resistance: f64,
    duration: f64,
    time: f64,
) -> Vec3 {
    let rate = resistance / mass;
    let slope = last
        .sub(first)
        .scale(1.0 / (resistance * duration))
        .add(ambient_after.sub(ambient).scale(1.0 / duration));
    let offset = first
        .add(gravity.scale(mass))
        .scale(1.0 / resistance)
        .add(ambient)
        .sub(slope.scale(1.0 / rate));
    offset
        .add(slope.scale(time))
        .add(initial.sub(offset).scale((-rate * time).exp()))
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
fn changing_force_and_bath_match_independent_ode_and_all_power_integrals() {
    let p = properties();
    let ambient = Vec3::new(0.5, -1.0, 0.25);
    let ambient_after = Vec3::new(-0.3, 0.6, -0.2);
    let gravity = Gravity::new(Vec3::new(0.2, 0.0, -0.3)).unwrap();
    let initial_v = Vec3::new(1.25, -0.5, 0.75);
    let initial = state(p, initial_v, Vec3::new(0.2, -0.3, 0.1));
    let first = Vec3::new(0.3, -1.2, 2.0);
    let last = Vec3::new(-1.0, 0.7, -0.2);
    let c = 4.0;
    let flow = IsotropicLinearDrag::new(p, gravity, ambient, c, 6.0).unwrap();
    for exposure in [0.125, 1.999_999, 2.0, 8.0] {
        let duration = exposure * p.mass() / c;
        let step = flow
            .step_affine_force_and_flow(initial, force(first), last, ambient_after, duration)
            .unwrap();
        let v = |t| {
            velocity(
                initial_v,
                ambient,
                ambient_after,
                first,
                last,
                gravity.acceleration_world(),
                p.mass(),
                c,
                duration,
                t,
            )
        };
        let u = |t| interpolate(ambient, ambient_after, t / duration);
        close_vector(
            step.state_after.center_of_mass_velocity_world(p).unwrap(),
            v(duration),
            3e-12,
        );
        close_vector(
            step.displacement_world,
            Vec3::new(
                simpson(|t| v(t).x, duration),
                simpson(|t| v(t).y, duration),
                simpson(|t| v(t).z, duration),
            ),
            3e-12,
        );
        close(
            step.translation_heat_j,
            simpson(|t| c * v(t).sub(u(t)).norm_squared(), duration),
            3e-12,
        );
        close(
            step.applied_work_j,
            simpson(
                |t| interpolate(first, last, t / duration).dot(v(t)),
                duration,
            ),
            3e-12,
        );
        close(
            step.background_work_j,
            simpson(|t| -c * v(t).sub(u(t)).dot(u(t)), duration),
            3e-12,
        );
        assert!(
            (step.background_work_j
                - step
                    .drag_impulse_world
                    .dot(interpolate(ambient, ambient_after, 0.5)))
            .abs()
                > 1e-4
        );
        close(step.work_energy_residual_j, 0.0, 5e-13);
        close_vector(step.linear_momentum_residual_world, Vec3::ZERO, 3e-13);
        let spin = flow.step(initial, Wrench::ZERO, duration).unwrap();
        assert_eq!(
            step.state_after.pose().orientation(),
            spin.state_after.pose().orientation()
        );
        assert_eq!(step.rotation_heat_j, spin.rotation_heat_j);
    }
}

#[test]
fn motion_and_work_compose_when_both_actual_source_endpoints_are_carried_forward() {
    let p = properties();
    let gravity = Gravity::new(Vec3::new(0.1, 0.0, -0.2)).unwrap();
    let ambient = Vec3::new(0.2, -0.4, 0.3);
    let ambient_after = Vec3::new(-0.5, 1.0, -0.2);
    let initial = state(p, Vec3::new(-0.5, 0.3, 1.0), Vec3::new(0.2, 0.4, -0.3));
    let first = Vec3::new(0.5, -0.7, 1.4);
    let last = Vec3::new(-0.4, 0.8, -0.1);
    let duration = 20.0;
    let make = |u| IsotropicLinearDrag::new(p, gravity, u, 4.0, 6.0).unwrap();
    let whole = make(ambient)
        .step_affine_force_and_flow(initial, force(first), last, ambient_after, duration)
        .unwrap();
    let mut current = initial;
    let mut time = 0.0;
    let (mut heat, mut work, mut bath) = (0.0, 0.0, 0.0);
    let mut impulse = Vec3::ZERO;
    for dt in [0.08, 0.17, 1.0, 18.75] {
        let part = make(interpolate(ambient, ambient_after, time / duration))
            .step_affine_force_and_flow(
                current,
                force(interpolate(first, last, time / duration)),
                interpolate(first, last, (time + dt) / duration),
                interpolate(ambient, ambient_after, (time + dt) / duration),
                dt,
            )
            .unwrap();
        current = part.state_after;
        time += dt;
        heat += part.dissipated_energy_j;
        work += part.applied_work_j;
        bath += part.background_work_j;
        impulse = impulse.add(part.drag_impulse_world);
    }
    close_vector(
        current.linear_momentum_world(),
        whole.state_after.linear_momentum_world(),
        3e-13,
    );
    close_vector(
        current.pose().position_world(),
        whole.state_after.pose().position_world(),
        3e-13,
    );
    close_vector(impulse, whole.drag_impulse_world, 3e-13);
    close(heat, whole.dissipated_energy_j, 3e-13);
    close(work, whole.applied_work_j, 3e-13);
    close(bath, whole.background_work_j, 3e-13);
}

#[test]
fn moving_bath_obeys_galilean_heat_impulse_and_both_work_transformations() {
    let p = properties();
    let ambient = Vec3::new(0.5, -0.25, 0.125);
    let ambient_after = Vec3::new(-0.25, 0.75, -0.5);
    let initial_v = Vec3::new(-0.25, 0.75, 0.5);
    let boost = Vec3::new(32.0, -16.0, 8.0);
    let first = Vec3::new(0.3, -0.5, 0.9);
    let last = Vec3::new(-0.7, 0.2, 0.1);
    let duration = 3.0;
    let step = |shift: Vec3| {
        IsotropicLinearDrag::new(p, Gravity::ZERO, ambient.add(shift), 4.0, 0.0)
            .unwrap()
            .step_affine_force_and_flow(
                state(p, initial_v.add(shift), Vec3::ZERO),
                force(first),
                last,
                ambient_after.add(shift),
                duration,
            )
            .unwrap()
    };
    let base = step(Vec3::ZERO);
    let shifted = step(boost);
    close(shifted.translation_heat_j, base.translation_heat_j, 3e-14);
    close_vector(shifted.drag_impulse_world, base.drag_impulse_world, 3e-14);
    close_vector(
        shifted.displacement_world.sub(base.displacement_world),
        boost.scale(duration),
        3e-14,
    );
    close(
        shifted.applied_work_j - base.applied_work_j,
        boost.dot(first.add(last).scale(0.5 * duration)),
        3e-14,
    );
    close(
        shifted.background_work_j - base.background_work_j,
        boost.dot(base.drag_impulse_world),
        3e-14,
    );
}

#[test]
fn zero_and_tiny_drag_keep_their_physical_limits_with_a_changing_bath() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let initial = state(p, Vec3::ZERO, Vec3::ZERO);
    let first = Vec3::new(0.5, -0.25, 0.125);
    let last = first.scale(-1.0);
    let bath = Vec3::new(f64::MAX, -f64::MAX, f64::MAX);
    let ballistic = IsotropicLinearDrag::new(p, Gravity::ZERO, bath, 0.0, 0.0).unwrap();
    assert_eq!(
        ballistic
            .step_affine_force(initial, force(first), last, 1.0)
            .unwrap(),
        ballistic
            .step_affine_force_and_flow(initial, force(first), last, bath.scale(-1.0), 1.0)
            .unwrap()
    );
    let c = 1e-20;
    let tiny = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, c, 0.0)
        .unwrap()
        .step_affine_force_and_flow(
            initial,
            Wrench::ZERO,
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
        )
        .unwrap();
    close(
        tiny.state_after.linear_momentum_world().x / (c * 0.5),
        1.0,
        2e-14,
    );
    assert!(tiny.translation_heat_j > 0.0);
    close(tiny.translation_heat_j / (c / 3.0), 1.0, 2e-14);
    close(tiny.background_work_j / (c / 3.0), 1.0, 2e-14);
    assert_eq!(tiny.applied_work_j, 0.0);
}

#[test]
fn matching_applied_acceleration_keeps_zero_slip_in_a_ramping_bath() {
    let p = properties();
    let ambient = Vec3::new(0.25, -0.5, 0.125);
    let acceleration = Vec3::new(0.5, -0.25, 0.125);
    let duration = 2.0;
    let ambient_after = ambient.add(acceleration.scale(duration));
    let applied = acceleration.scale(p.mass());
    for c in [0.125, 2.0, 32.0] {
        let step = IsotropicLinearDrag::new(p, Gravity::ZERO, ambient, c, 0.0)
            .unwrap()
            .step_affine_force_and_flow(
                state(p, ambient, Vec3::ZERO),
                force(applied),
                applied,
                ambient_after,
                duration,
            )
            .unwrap();
        assert_eq!(step.translation_heat_j, 0.0);
        assert_eq!(step.background_work_j, 0.0);
        assert_eq!(step.drag_impulse_world, Vec3::ZERO);
        close_vector(
            step.state_after.center_of_mass_velocity_world(p).unwrap(),
            ambient_after,
            2e-14,
        );
        close_vector(
            step.displacement_world,
            ambient.add(ambient_after).scale(0.5 * duration),
            2e-14,
        );
        close(step.work_energy_residual_j, 0.0, 2e-14);
    }
}

#[test]
fn held_flow_parity_and_relative_path_bound_cover_interior_slip() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let initial = state(p, Vec3::ZERO, Vec3::ZERO);
    let flow = IsotropicLinearDrag::new(p, Gravity::ZERO, Vec3::ZERO, 1.0, 0.0).unwrap();
    let first = Vec3::new(2.0, 0.0, 0.0);
    let last = Vec3::new(3.0 - core::f64::consts::E, 0.0, 0.0);
    let ambient_after = Vec3::new(1.0, 0.0, 0.0);
    assert_eq!(
        flow.step_affine_force(initial, force(first), last, 1.0)
            .unwrap(),
        flow.step_affine_force_and_flow(initial, force(first), last, Vec3::ZERO, 1.0)
            .unwrap()
    );
    let step = flow
        .step_affine_force_and_flow(initial, force(first), last, ambient_after, 1.0)
        .unwrap();
    close(step.state_after.linear_momentum_world().x, 1.0, 2e-14);
    let peak_time = -(1.0 - (-1.0_f64).exp()).ln();
    let peak = velocity(
        Vec3::ZERO,
        Vec3::ZERO,
        ambient_after,
        first,
        last,
        Vec3::ZERO,
        1.0,
        1.0,
        1.0,
        peak_time,
    )
    .x - peak_time;
    assert!(peak > 0.21);
    let bound = flow
        .affine_force_and_flow_relative_speed_bound_m_s(
            initial,
            force(first),
            last,
            ambient_after,
            1.0,
        )
        .unwrap();
    close(bound, -(-1.0_f64).exp_m1(), 2e-14);
    assert!(bound >= peak);
    for invalid in [f64::NAN, f64::INFINITY] {
        let end = Vec3::new(invalid, 0.0, 0.0);
        assert!(
            flow.step_affine_force_and_flow(initial, force(first), last, end, 1.0)
                .is_err()
        );
        assert!(
            flow.affine_force_and_flow_relative_speed_bound_m_s(
                initial,
                force(first),
                last,
                end,
                1.0
            )
            .is_err()
        );
    }
}
