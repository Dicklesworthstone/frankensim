//! G1/G3/G4: nonlinear drag trajectories, reservoir work and regime refusal.
use super::*;
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion, Vec3 as NativeVec3};
use fs_mbd_material::creeping_flow::finite_reynolds::{FiniteReynoldsLimits, FiniteReynoldsSphere};

fn finite_limits() -> FiniteReynoldsLimits {
    FiniteReynoldsLimits {
        maximum_reynolds: 1000.0,
        maximum_diffusion_ratio: 0.05,
    }
}

fn state_at(
    body: &fs_mbd_material::principal::PrincipalBody<'_>,
    velocity: NativeVec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(NativeVec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(body.properties().mass()),
        NativeVec3::ZERO,
    )
    .unwrap()
}

#[test]
fn finite_reynolds_force_is_regular_at_rest_and_recovers_stokes_at_small_slip() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 0.001);
    let flow = FiniteReynoldsSphere::new(
        &material,
        RADIUS,
        Gravity::ZERO,
        NativeVec3::ZERO,
        finite_limits(),
    )
    .unwrap();
    assert!(core::ptr::eq(flow.source(), &material));
    assert_eq!(flow.transport(), material.transport_state().unwrap());
    let c = 6.0 * core::f64::consts::PI * 0.001 * RADIUS;
    let rest = flow
        .load(&body, state_at(&body, NativeVec3::ZERO), || false)
        .unwrap();
    assert_eq!(rest.drag_force_world_n, NativeVec3::ZERO);
    assert_eq!(rest.dissipated_power_w, 0.0);
    close(rest.tangent_resistance_kg_s, c);
    let slow = flow
        .load(
            &body,
            state_at(&body, NativeVec3::new(5e-11, 0.0, 0.0)),
            || false,
        )
        .unwrap();
    close(slow.reynolds, 1e-9);
    assert!((-slow.drag_force_world_n.x / (c * 5e-11) - 1.0).abs() < 2e-7);
    let finite = flow
        .load(
            &body,
            state_at(&body, NativeVec3::new(5.0, 0.0, 0.0)),
            || false,
        )
        .unwrap();
    close(finite.reynolds, 100.0);
    // The independent coefficient form F=-rho*A*Cd*v²/2 uses diameter Re.
    let cd = 24.0 / 100.0 * (1.0 + 0.15 * 100.0_f64.powf(0.687));
    close(
        finite.drag_force_world_n.x,
        -0.5 * core::f64::consts::PI * RADIUS.powi(2) * cd * 25.0,
    );
    assert!(finite.drag_force_world_n.x < -4.0 * c * 5.0);
    close(
        finite.dissipated_power_w,
        -finite.drag_force_world_n.x * 5.0,
    );
    let force_at = |speed| {
        flow.load(
            &body,
            state_at(&body, NativeVec3::new(speed, 0.0, 0.0)),
            || false,
        )
        .unwrap()
        .drag_force_world_n
        .x
    };
    let derivative = (force_at(5.0 - 1e-4) - force_at(5.0 + 1e-4)) / 2e-4;
    assert!((derivative / finite.tangent_resistance_kg_s - 1.0).abs() < 1e-9);
}

// A closed-form nonlinear decay oracle, independent of midpoint/Newton:
// v^-p(t) = (v0^-p + beta)*exp(p*alpha*t) - beta.
fn exact_decay(t: f64, body_density: f64, fluid_density: f64, viscosity: f64) -> f64 {
    let p = 0.687;
    let alpha = 9.0 * viscosity / (2.0 * body_density * RADIUS.powi(2));
    let initial_correction = 0.15 * (2.0 * RADIUS * 5.0 * fluid_density / viscosity).powf(p);
    5.0 / (1.0 + (1.0 + initial_correction) * (p * alpha * t).exp_m1()).powf(1.0 / p)
}

#[test]
fn nonlinear_native_motion_converges_at_second_order_to_closed_form_decay() {
    for (body_density, fluid_density, viscosity) in [(20000.0, 1.0, 0.001), (10000.0, 2.0, 0.002)] {
        let source = sphere(body_density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let material = fluid(fluid_density, viscosity);
        let flow = FiniteReynoldsSphere::new(
            &material,
            RADIUS,
            Gravity::ZERO,
            NativeVec3::ZERO,
            finite_limits(),
        )
        .unwrap();
        let exact_velocity = exact_decay(40.0, body_density, fluid_density, viscosity);
        // Composite Simpson quadrature of the independent exact velocity gives
        // a position oracle much tighter than the three native step sizes.
        let n = 8192;
        let h = 40.0 / f64::from(n);
        let mut exact_position = 0.0;
        for i in 0..=n {
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            exact_position +=
                weight * exact_decay(f64::from(i) * h, body_density, fluid_density, viscosity);
        }
        exact_position *= h / 3.0;
        let mut errors = Vec::new();
        for steps in [20, 40, 80] {
            let dt = 40.0 / f64::from(steps);
            let mut state = state_at(&body, NativeVec3::new(5.0, 0.0, 0.0));
            let initial_energy = state.kinetic_energy(body.properties()).unwrap();
            let mut heat = 0.0;
            for _ in 0..steps {
                let result = flow.step(&body, state, controls(dt), || false).unwrap();
                let native = RigidBodyIntegrator::new(Gravity::ZERO)
                    .step(state, body.properties(), result.load_midpoint.wrench, dt)
                    .unwrap();
                assert_eq!(result.dynamics, native);
                assert!(result.midpoint_velocity_residual_m_s < 2e-14);
                assert!(result.viscous_heat_j > 0.0);
                let previous = state
                    .center_of_mass_velocity_world(body.properties())
                    .unwrap()
                    .x;
                state = result.dynamics.state_after;
                let current = state
                    .center_of_mass_velocity_world(body.properties())
                    .unwrap()
                    .x;
                assert!(current > 0.0 && current < previous);
                heat += result.viscous_heat_j;
            }
            close(
                state.kinetic_energy(body.properties()).unwrap() + heat,
                initial_energy,
            );
            errors.push((
                (state
                    .center_of_mass_velocity_world(body.properties())
                    .unwrap()
                    .x
                    - exact_velocity)
                    .abs(),
                (state.pose().position_world().x - exact_position).abs(),
            ));
        }
        for pair in errors.windows(2) {
            assert!(pair[0].0 / pair[1].0 > 3.8, "velocity errors {errors:?}");
            assert!(pair[0].1 / pair[1].1 > 3.8, "position errors {errors:?}");
        }
    }
}

#[test]
fn finite_reynolds_buoyancy_and_drag_preserve_the_terminal_settling_solution() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 0.001);
    let c = 6.0 * core::f64::consts::PI * 0.001 * RADIUS * (1.0 + 0.15 * 100.0_f64.powf(0.687));
    let displaced_mass = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let acceleration = c * 5.0 / (body.properties().mass() - displaced_mass);
    let gravity = Gravity::new(NativeVec3::new(0.0, 0.0, -acceleration)).unwrap();
    let flow = FiniteReynoldsSphere::new(
        &material,
        RADIUS,
        gravity,
        NativeVec3::ZERO,
        finite_limits(),
    )
    .unwrap();
    let state = state_at(&body, NativeVec3::new(0.0, 0.0, -5.0));
    let result = flow.step(&body, state, controls(2.0), || false).unwrap();
    close(
        result
            .dynamics
            .state_after
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z,
        -5.0,
    );
    close(result.dynamics.state_after.pose().position_world().z, -10.0);
    close(result.stored_energy_change_j + result.viscous_heat_j, 0.0);
}

#[test]
fn maintained_flow_supplies_work_and_common_boost_preserves_relative_motion() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 0.001);
    let wind = NativeVec3::new(3.0, 4.0, 0.0);
    let boost = NativeVec3::new(-2.0, 1.0, 3.0);
    let flow =
        FiniteReynoldsSphere::new(&material, RADIUS, Gravity::ZERO, wind, finite_limits()).unwrap();
    let shifted = FiniteReynoldsSphere::new(
        &material,
        RADIUS,
        Gravity::ZERO,
        wind.add(boost),
        finite_limits(),
    )
    .unwrap();
    let result = flow
        .step(
            &body,
            state_at(&body, NativeVec3::ZERO),
            controls(2.0),
            || false,
        )
        .unwrap();
    let moved = shifted
        .step(&body, state_at(&body, boost), controls(2.0), || false)
        .unwrap();
    assert!(result.stored_energy_change_j > 0.0);
    assert!(result.background_work_j > result.viscous_heat_j);
    close(
        result.background_work_j,
        result.stored_energy_change_j + result.viscous_heat_j,
    );
    close(moved.viscous_heat_j, result.viscous_heat_j);
    let dv = moved
        .dynamics
        .state_after
        .center_of_mass_velocity_world(body.properties())
        .unwrap()
        .sub(
            result
                .dynamics
                .state_after
                .center_of_mass_velocity_world(body.properties())
                .unwrap(),
        )
        .sub(boost);
    assert!(
        dv.norm_squared() < 1e-24,
        "boost changed relative velocity: {dv:?}"
    );
    let dx = moved
        .dynamics
        .state_after
        .pose()
        .position_world()
        .sub(result.dynamics.state_after.pose().position_world())
        .sub(boost.scale(2.0));
    assert!(
        dx.norm_squared() < 1e-24,
        "boost changed relative displacement: {dx:?}"
    );
    close(
        moved.background_work_j - result.background_work_j,
        result.load_midpoint.drag_force_world_n.dot(boost) * 2.0,
    );
}

#[test]
fn nonlinear_regime_spin_large_steps_and_late_cancellation_refuse() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 0.001);
    let flow = FiniteReynoldsSphere::new(
        &material,
        RADIUS,
        Gravity::ZERO,
        NativeVec3::ZERO,
        finite_limits(),
    )
    .unwrap();
    let state = state_at(&body, NativeVec3::new(5.0, 0.0, 0.0));
    assert!(
        flow.load(
            &body,
            state_at(&body, NativeVec3::new(50.0, 0.0, 0.0)),
            || false
        )
        .is_ok()
    );
    assert!(matches!(
        flow.load(
            &body,
            state_at(&body, NativeVec3::new(50.01, 0.0, 0.0)),
            || false
        ),
        Err(CreepingFlowError::Reynolds { .. })
    ));
    let spinning = RigidBodyState::new(
        state.pose(),
        state.linear_momentum_world(),
        NativeVec3::new(0.0, 0.0, 1e-8),
    )
    .unwrap();
    assert!(flow.load(&body, spinning, || false).is_err());
    assert!(
        flow.step(&body, state, controls(10000.0), || false)
            .is_err()
    );
    let tangent = flow
        .load(&body, state, || false)
        .unwrap()
        .tangent_resistance_kg_s;
    let strong_step = flow
        .step(
            &body,
            state,
            CreepingStepControls {
                maximum_relaxation_fraction: 1.0,
                ..controls(0.99 * body.properties().mass() / tangent)
            },
            || false,
        )
        .unwrap();
    let strong_velocity = strong_step
        .dynamics
        .state_after
        .center_of_mass_velocity_world(body.properties())
        .unwrap()
        .x;
    assert!(strong_velocity > 0.0 && strong_velocity < 5.0);
    let light = sphere(1.0);
    let light_body = light.principal_body(1e-12, || false).unwrap();
    assert!(matches!(
        flow.step(
            &light_body,
            state_at(&light_body, NativeVec3::ZERO),
            controls(0.001),
            || false
        ),
        Err(CreepingFlowError::Unsteady { .. })
    ));
    let mut calls = 0;
    flow.step(&body, state, controls(2.0), || {
        calls += 1;
        false
    })
    .unwrap();
    let total = calls;
    calls = 0;
    assert!(matches!(
        flow.step(&body, state, controls(2.0), || {
            calls += 1;
            calls == total
        }),
        Err(CreepingFlowError::Cancelled)
    ));
    assert_eq!(
        state
            .center_of_mass_velocity_world(body.properties())
            .unwrap(),
        NativeVec3::new(5.0, 0.0, 0.0)
    );
}
