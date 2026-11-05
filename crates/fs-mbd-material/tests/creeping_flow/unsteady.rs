//! G1/G3/G4: analytic Basset impedance/settling, source changes and atomic memory.
use super::*;
use fs_mbd::{Vec3 as NativeVec3, Wrench};
use fs_mbd_material::creeping_flow::unsteady::{BassetResolution, UnsteadyStepControls};

fn resolution(lo: f64, hi: f64, intervals: usize) -> BassetResolution {
    BassetResolution {
        minimum_wavenumber_radius: lo,
        maximum_wavenumber_radius: hi,
        intervals,
    }
}

fn step_controls(dt: f64) -> UnsteadyStepControls {
    UnsteadyStepControls {
        duration_s: dt,
        maximum_energy_defect_j: 1e-18,
        maximum_momentum_defect_ns: 1e-16,
    }
}

#[test]
fn positive_memory_ladder_converges_to_the_analytic_unsteady_sphere_impedance() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let c = 6.0 * core::f64::consts::PI * 1.0 * RADIUS;
    let diffusion_time = RADIUS.powi(2) / 0.001;
    let added_mass = 0.5 * 1000.0 * 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    for omega in [1.0, 10.0, 100.0] {
        let history = c * (0.5_f64 * omega * diffusion_time).sqrt();
        let exact = [c + history, added_mass * omega + history];
        let mut previous = f64::INFINITY;
        for (lo, hi, n) in [
            (0.02, 40.0, 32),
            (0.001, 1000.0, 128),
            (0.0001, 10000.0, 256),
        ] {
            let motion = medium
                .unsteady_motion(&body, rest, resolution(lo, hi, n), || false)
                .unwrap();
            close(motion.added_mass_kg(), added_mass);
            assert!(
                motion
                    .branches()
                    .iter()
                    .all(|b| b.mass_kg > 0.0 && b.resistance_kg_s > 0.0)
            );
            let actual = motion.impedance_kg_s(omega).unwrap();
            let error = (actual[0] - exact[0]).hypot(actual[1] - exact[1]);
            assert!(
                error < previous * 0.15,
                "spectral ladder failed: {error} >= {previous}"
            );
            // This fixture's fine quadrature error is negligible compared with
            // the independently bounded two omitted continuum tails.
            assert!(error <= 1.001 * motion.omitted_tail_bound_kg_s(omega).unwrap());
            previous = error;
        }
        assert!(previous / exact[0].hypot(exact[1]) < 5e-4);
    }
}

#[test]
fn source_density_sets_buoyant_direction_and_actual_added_mass_acceleration() {
    for (solid_density, fluid_density) in [
        (500.0, 1000.0),
        (1000.0, 1000.0),
        (2000.0, 1000.0),
        (2000.0, 500.0),
    ] {
        let material = fluid(fluid_density, 1.0);
        let source = sphere(solid_density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let medium =
            CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits())
                .unwrap();
        let rest = body
            .state_from_kinematics(Se3::identity(), Twist::zero())
            .unwrap();
        // The old quasi-steady model correctly refuses these density ratios.
        assert!(matches!(
            medium.step(&body, rest, controls(1e-10), || false),
            Err(CreepingFlowError::Unsteady { .. })
        ));
        let motion = medium
            .unsteady_motion(&body, rest, resolution(1e-4, 1000.0, 128), || false)
            .unwrap();
        let dt = 1e-10;
        let (next, audit) = motion.step(step_controls(dt), || false).unwrap();
        let velocity = next
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z;
        let acceleration =
            -0.01 * (solid_density - fluid_density) / (solid_density + 0.5 * fluid_density);
        if acceleration == 0.0 {
            assert_eq!(next.state(), rest);
            assert!(
                next.fluid_velocities_world_m_s()
                    .iter()
                    .all(|u| *u == fs_mbd::Vec3::ZERO)
            );
            assert_eq!(audit.viscous_heat_j, 0.0);
        } else {
            assert!((velocity / (dt * acceleration) - 1.0).abs() < 3e-6);
            assert_eq!(velocity.is_sign_positive(), solid_density < fluid_density);
            assert!(audit.added_mass_energy_change_j > 0.0);
            assert!(audit.memory_energy_change_j > 0.0);
            assert!(audit.viscous_heat_j > 0.0);
        }
        let native = RigidBodyIntegrator::new(gravity())
            .step(rest, body.properties(), audit.fluid_wrench, dt)
            .unwrap();
        assert_eq!(native, audit.dynamics);
        assert_eq!(next.state(), native.state_after);
    }
}

/// Branch-cut inversion of V(s)=F/[s*(M*s+b*sqrt(s)+c)] from a quiescent start.
/// With alpha=b/sqrt(M*c), the normalized velocity is 1 minus this positive
/// decaying integral. Simpson quadrature is independent of the production
/// rational branch approximation and native midpoint stepping.
fn settling_fraction(time_s: f64, moving_mass: f64, c: f64, b: f64) -> f64 {
    let alpha = b / (moving_mass * c).sqrt();
    let scaled_time = time_s * c / moving_mass;
    let n = 8192;
    let dx = 12.0 / n as f64;
    let integrand =
        |x: f64| (-scaled_time * x * x).exp() / ((1.0 - x * x).powi(2) + alpha * alpha * x * x);
    let mut integral = integrand(0.0) + integrand(12.0);
    for i in 1..n {
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(i as f64 * dx);
    }
    1.0 - (2.0 * alpha / core::f64::consts::PI) * integral * dx / 3.0
}

#[test]
fn native_settling_converges_to_basset_history_with_closed_storage_and_heat() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let volume = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let moving_mass = 2500.0 * volume;
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let b = c * (RADIUS.powi(2) / 0.001).sqrt();
    let terminal = -1000.0 * volume * 0.01 / c;
    let target = 0.1;
    let expected = terminal * settling_fraction(target, moving_mass, c, b);
    let mut previous = f64::INFINITY;
    for (lo, hi, branches, steps) in [
        (0.01, 100.0, 32, 50),
        (0.001, 1000.0, 128, 200),
        (0.0001, 10000.0, 256, 800),
    ] {
        let mut motion = medium
            .unsteady_motion(&body, rest, resolution(lo, hi, branches), || false)
            .unwrap();
        let mut heat = 0.0;
        let mut stored = 0.0;
        let mut fluid_storage = 0.0;
        for _ in 0..steps {
            let (next, audit) = motion
                .step(step_controls(target / steps as f64), || false)
                .unwrap();
            heat += audit.viscous_heat_j;
            stored += audit.stored_energy_change_j;
            fluid_storage += audit.memory_energy_change_j + audit.added_mass_energy_change_j;
            assert!(audit.momentum_residual_ns.norm_squared().sqrt() < 1e-18);
            motion = next;
        }
        assert!(fluid_storage > 0.0 && heat > 0.0);
        assert!((stored + heat).abs() < 1e-11 * heat);
        let velocity = motion
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z;
        let relative_error = (velocity / expected - 1.0).abs();
        assert!(
            relative_error < previous * 0.2,
            "history ladder: {relative_error} >= {previous}"
        );
        previous = relative_error;
    }
    assert!(previous < 3e-4);
    // Fluid history changes the trajectory, rather than just its diagnostics.
    let no_history = terminal * (1.0 - (-c * target / moving_mass).exp());
    assert!((expected / no_history - 1.0).abs() > 0.1);
}

#[test]
fn restart_and_late_cancellation_preserve_the_same_body_and_fluid_history() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, gravity(), fs_mbd::Vec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let original = medium
        .unsteady_motion(&body, rest, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let (first, _) = original.step(step_controls(0.001), || false).unwrap();
    let (whole, audit) = first.step(step_controls(0.001), || false).unwrap();
    let (resumed, repeated) = first.clone().step(step_controls(0.001), || false).unwrap();
    assert_eq!(whole.state(), resumed.state());
    assert_eq!(
        whole.fluid_velocities_world_m_s(),
        resumed.fluid_velocities_world_m_s()
    );
    assert_eq!(audit, repeated);
    let mut polls = 0;
    first
        .step(step_controls(0.001), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        first.step(step_controls(0.001), || {
            visited += 1;
            visited == polls
        }),
        Err(CreepingFlowError::Cancelled)
    ));
    let (retried, retried_audit) = first.step(step_controls(0.001), || false).unwrap();
    assert_eq!(retried.state(), whole.state());
    assert_eq!(
        retried.fluid_velocities_world_m_s(),
        whole.fluid_velocities_world_m_s()
    );
    assert_eq!(retried_audit, audit);
    assert_eq!(original.state(), rest);
    assert_eq!(original.time_s(), 0.0);
    assert!(
        medium
            .unsteady_motion(&body, whole.state(), resolution(0.001, 1000.0, 64), || {
                false
            })
            .is_err()
    );
}

#[test]
fn imposed_force_converges_to_independent_basset_step_response_in_three_axes() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let applied = Wrench {
        force_world: NativeVec3::new(1e-4, -2e-4, 3e-4),
        torque_body: NativeVec3::ZERO,
    };
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let moving_mass = 2500.0 * 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let b = c * (RADIUS.powi(2) / 0.001).sqrt();
    let target = 0.1;
    // Independent branch-cut quadrature of the continuum transfer function.
    let exact = applied
        .force_world
        .scale(settling_fraction(target, moving_mass, c, b) / c);
    let mut previous = f64::INFINITY;
    for (lo, hi, branches, steps) in [
        (0.01, 100.0, 32, 50),
        (0.001, 1000.0, 128, 200),
        (0.0001, 10000.0, 256, 800),
    ] {
        let mut motion = medium
            .unsteady_motion(&body, rest, resolution(lo, hi, branches), || false)
            .unwrap();
        let dt = target / f64::from(steps);
        let mut heat = 0.0;
        let mut work = 0.0;
        let mut stored = 0.0;
        for _ in 0..steps {
            let before = motion.state();
            let (next, audit) = motion
                .step_with_wrench(applied, step_controls(dt), || false)
                .unwrap();
            let displacement = next
                .state()
                .pose()
                .position_world()
                .sub(before.pose().position_world());
            assert!((audit.applied_work_j - applied.force_world.dot(displacement)).abs() < 1e-18);
            assert!(audit.viscous_heat_j >= 0.0);
            assert!(audit.momentum_residual_ns.norm_squared().sqrt() < 1e-16);
            heat += audit.viscous_heat_j;
            work += audit.applied_work_j;
            stored += audit.stored_energy_change_j;
            motion = next;
        }
        let actual = motion
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap();
        let relative_error = actual.sub(exact).norm_squared().sqrt() / exact.norm_squared().sqrt();
        assert!(
            relative_error < previous * 0.2,
            "driven Basset ladder: {relative_error} >= {previous}"
        );
        previous = relative_error;
        assert!(work > heat && heat > 0.0);
        assert!((stored + heat - work).abs() < 1e-11 * work);
    }
    assert!(previous < 3e-4);
}

#[test]
fn changing_drives_retain_fluid_memory_and_close_native_work_and_impulse() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let mut motion = medium
        .unsteady_motion(&body, rest, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let dt = 0.001;
    let mut applied_work = 0.0;
    let mut heat = 0.0;
    let mut imposed_impulse = NativeVec3::ZERO;
    let mut exterior_impulse = NativeVec3::ZERO;
    for step in 0..400 {
        let t = (f64::from(step) + 0.5) * dt;
        let applied = Wrench {
            force_world: NativeVec3::new(
                1e-4 * (25.0 * t).cos(),
                2e-4 * (17.0 * t).sin(),
                if step < 200 { 1e-4 } else { -1e-4 },
            ),
            torque_body: NativeVec3::ZERO,
        };
        let before = motion.state();
        let (next, audit) = motion
            .step_with_wrench(applied, step_controls(dt), || false)
            .unwrap();
        let total = Wrench {
            force_world: audit.fluid_wrench.force_world.add(applied.force_world),
            torque_body: NativeVec3::ZERO,
        };
        // The original native owner supplies the accepted body trajectory.
        let native = RigidBodyIntegrator::new(Gravity::ZERO)
            .step(before, body.properties(), total, dt)
            .unwrap();
        assert_eq!(native, audit.dynamics);
        applied_work += applied.force_world.dot(
            next.state()
                .pose()
                .position_world()
                .sub(before.pose().position_world()),
        );
        heat += audit.viscous_heat_j;
        imposed_impulse = imposed_impulse.add(applied.force_world.scale(dt));
        exterior_impulse = exterior_impulse.add(audit.far_field_impulse_ns);
        motion = next;
    }
    let velocity = motion
        .state()
        .center_of_mass_velocity_world(body.properties())
        .unwrap();
    let mut energy =
        0.5 * (body.properties().mass() + motion.added_mass_kg()) * velocity.norm_squared();
    let mut momentum = velocity.scale(body.properties().mass() + motion.added_mass_kg());
    for (branch, fluid_velocity) in motion
        .branches()
        .iter()
        .zip(motion.fluid_velocities_world_m_s())
    {
        energy += 0.5 * branch.mass_kg * fluid_velocity.norm_squared();
        momentum = momentum.add(fluid_velocity.scale(branch.mass_kg));
    }
    assert!(
        motion
            .fluid_velocities_world_m_s()
            .iter()
            .any(|v| v.norm_squared() > 0.0)
    );
    assert!((energy + heat - applied_work).abs() < 1e-11 * applied_work.abs());
    assert!(
        momentum
            .add(exterior_impulse)
            .sub(imposed_impulse)
            .norm_squared()
            .sqrt()
            < 1e-15
    );
    let (ordinary, ordinary_audit) = motion.step(step_controls(dt), || false).unwrap();
    let (explicit_zero, zero_audit) = motion
        .step_with_wrench(Wrench::ZERO, step_controls(dt), || false)
        .unwrap();
    assert_eq!(ordinary.state(), explicit_zero.state());
    assert_eq!(
        ordinary.fluid_velocities_world_m_s(),
        explicit_zero.fluid_velocities_world_m_s()
    );
    assert_eq!(ordinary_audit, zero_audit);
}

#[test]
fn driven_unsteady_refusal_and_late_cancellation_leave_all_history_unchanged() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let motion = medium
        .unsteady_motion(&body, rest, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let drive = Wrench {
        force_world: NativeVec3::new(1e-4, 0.0, 0.0),
        torque_body: NativeVec3::ZERO,
    };
    let (motion, _) = motion
        .step_with_wrench(drive, step_controls(0.001), || false)
        .unwrap();
    let state = motion.state();
    let time = motion.time_s();
    let history = motion.fluid_velocities_world_m_s().to_vec();
    let mut polls = 0;
    let (expected, expected_audit) = motion
        .step_with_wrench(drive, step_controls(0.001), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        motion.step_with_wrench(drive, step_controls(0.001), || {
            visited += 1;
            visited == polls
        }),
        Err(CreepingFlowError::Cancelled)
    ));
    for invalid in [
        Wrench {
            force_world: NativeVec3::new(f64::NAN, 0.0, 0.0),
            ..Wrench::ZERO
        },
        Wrench {
            torque_body: NativeVec3::new(0.0, 1e-12, 0.0),
            ..Wrench::ZERO
        },
        Wrench {
            force_world: NativeVec3::new(1e20, 0.0, 0.0),
            ..Wrench::ZERO
        },
    ] {
        assert!(
            motion
                .step_with_wrench(invalid, step_controls(0.001), || false)
                .is_err()
        );
    }
    assert_eq!(motion.state(), state);
    assert_eq!(motion.time_s(), time);
    assert_eq!(motion.fluid_velocities_world_m_s(), history);
    let (retried, retried_audit) = motion
        .step_with_wrench(drive, step_controls(0.001), || false)
        .unwrap();
    assert_eq!(retried.state(), expected.state());
    assert_eq!(
        retried.fluid_velocities_world_m_s(),
        expected.fluid_velocities_world_m_s()
    );
    assert_eq!(retried_audit, expected_audit);
}

fn vector_close(actual: NativeVec3, expected: NativeVec3, tolerance: f64) {
    let error = actual.sub(expected);
    assert!(
        error.x.hypot(error.y).hypot(error.z) <= tolerance,
        "{actual:?} != {expected:?} within {tolerance:e}"
    );
}

#[test]
fn neutrally_buoyant_sphere_follows_accelerating_fluid_without_slip_heat() {
    let source = sphere(1000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let start = body
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.13, -0.07, 0.2)),
        )
        .unwrap();
    let mut velocity = NativeVec3::new(0.13, -0.07, 0.2);
    // The large common velocity is deliberately outside the old absolute-velocity
    // Reynolds gate. Only slip is a creeping-flow criterion.
    let medium = CreepingSphere::new(&material, RADIUS, gravity(), velocity, limits()).unwrap();
    let mut motion = medium
        .unsteady_motion(&body, start, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let mut position = start.pose().position_world();
    let dt = 0.002;
    for i in 0..100 {
        let acceleration = if i < 50 {
            NativeVec3::new(0.02, -0.01, 0.03)
        } else {
            NativeVec3::new(-0.04, 0.03, -0.01)
        };
        position = position.add(velocity.add(acceleration.scale(0.5 * dt)).scale(dt));
        velocity = velocity.add(acceleration.scale(dt));
        let (next, audit) = motion
            .step_with_flow_acceleration(Wrench::ZERO, acceleration, step_controls(dt), || false)
            .unwrap();
        vector_close(next.ambient_velocity_world_m_s(), velocity, 1e-15);
        vector_close(
            next.state()
                .center_of_mass_velocity_world(body.properties())
                .unwrap(),
            velocity,
            2e-15,
        );
        vector_close(next.state().pose().position_world(), position, 2e-15);
        vector_close(
            audit.ambient_pressure_force_world_n,
            acceleration.scale(body.properties().mass()),
            1e-18,
        );
        assert!(
            audit.viscous_heat_j < 1e-28,
            "neutral heat {:?}",
            audit.viscous_heat_j
        );
        assert!((audit.stored_energy_change_j - audit.background_work_j).abs() < 1e-18);
        assert_eq!(audit.applied_work_j, 0.0);
        motion = next;
    }
}

#[test]
fn ambient_pressure_and_added_mass_give_the_correct_density_dependent_startup() {
    let material = fluid(1000.0, 1.0);
    let acceleration = NativeVec3::new(0.02, -0.03, 0.01);
    let dt = 1e-10;
    for density in [500.0, 1000.0, 2000.0] {
        let source = sphere(density);
        let body = source.principal_body(1e-12, || false).unwrap();
        let medium =
            CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits())
                .unwrap();
        let rest = body
            .state_from_kinematics(Se3::identity(), Twist::zero())
            .unwrap();
        let motion = medium
            .unsteady_motion(&body, rest, resolution(1e-4, 1000.0, 128), || false)
            .unwrap();
        let (next, audit) = motion
            .step_with_flow_acceleration(Wrench::ZERO, acceleration, step_controls(dt), || false)
            .unwrap();
        let actual = next
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .scale(1.0 / dt);
        let expected = acceleration.scale(1500.0 / (density + 500.0));
        vector_close(actual, expected, 2e-7);
        let displaced_mass = 1000.0 * 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
        vector_close(
            audit.ambient_pressure_force_world_n,
            acceleration.scale(displaced_mass),
            1e-18,
        );
        let native = RigidBodyIntegrator::new(Gravity::ZERO)
            .step(rest, body.properties(), audit.fluid_wrench, dt)
            .unwrap();
        assert_eq!(native, audit.dynamics);
    }
}

#[test]
fn accelerating_current_converges_to_independent_continuum_basset_response() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let acceleration = NativeVec3::new(0.01, -0.02, 0.015);
    let volume = 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let moving_mass = 2500.0 * volume;
    let b = c * (RADIUS.powi(2) / 0.001).sqrt();
    let target = 0.1;
    let slip =
        acceleration.scale(-1000.0 * volume * settling_fraction(target, moving_mass, c, b) / c);
    let expected = acceleration.scale(target).add(slip);
    let mut previous = f64::INFINITY;
    for (lo, hi, branches, steps) in [
        (0.01, 100.0, 32, 50),
        (0.001, 1000.0, 128, 200),
        (0.0001, 10000.0, 256, 800),
    ] {
        let mut motion = medium
            .unsteady_motion(&body, rest, resolution(lo, hi, branches), || false)
            .unwrap();
        let mut accounts = [0.0; 3];
        for _ in 0..steps {
            let (next, audit) = motion
                .step_with_flow_acceleration(
                    Wrench::ZERO,
                    acceleration,
                    step_controls(target / f64::from(steps)),
                    || false,
                )
                .unwrap();
            accounts[0] += audit.stored_energy_change_j;
            accounts[1] += audit.viscous_heat_j;
            accounts[2] += audit.background_work_j;
            motion = next;
        }
        let actual = motion
            .state()
            .center_of_mass_velocity_world(body.properties())
            .unwrap();
        let error = actual.sub(expected).norm_squared().sqrt() / slip.norm_squared().sqrt();
        assert!(
            error < 0.2 * previous,
            "accelerating-flow ladder {error:e} >= {previous:e}"
        );
        assert!(accounts[1] > 0.0 && accounts[2] > 0.0);
        assert!((accounts[0] + accounts[1] - accounts[2]).abs() < 1e-11 * accounts[2]);
        previous = error;
    }
    assert!(previous < 3e-4);
}

#[test]
fn ambient_translation_preserves_slip_memory_and_heat_under_a_galilean_boost() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let boosted = body
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.125, -0.25, 0.0625)),
        )
        .unwrap();
    let boost = boosted
        .center_of_mass_velocity_world(body.properties())
        .unwrap();
    let still_medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let moving_medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, boost, limits()).unwrap();
    let mut still = still_medium
        .unsteady_motion(&body, rest, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let mut moving = moving_medium
        .unsteady_motion(&body, boosted, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let dt = 0.001;
    for i in 0..100 {
        let acceleration = if i < 50 {
            NativeVec3::new(0.01, -0.02, 0.015)
        } else {
            NativeVec3::ZERO
        };
        let applied = Wrench {
            force_world: NativeVec3::new(1e-5, -2e-5, 0.0),
            ..Wrench::ZERO
        };
        let (next_still, a) = still
            .step_with_flow_acceleration(applied, acceleration, step_controls(dt), || false)
            .unwrap();
        let (next_moving, b) = moving
            .step_with_flow_acceleration(applied, acceleration, step_controls(dt), || false)
            .unwrap();
        vector_close(
            next_moving
                .state()
                .center_of_mass_velocity_world(body.properties())
                .unwrap(),
            next_still
                .state()
                .center_of_mass_velocity_world(body.properties())
                .unwrap()
                .add(boost),
            2e-14,
        );
        vector_close(
            next_moving.state().pose().position_world(),
            next_still
                .state()
                .pose()
                .position_world()
                .add(boost.scale(next_still.time_s())),
            2e-15,
        );
        for (a, b) in next_still
            .fluid_velocities_world_m_s()
            .iter()
            .zip(next_moving.fluid_velocities_world_m_s())
        {
            vector_close(*a, *b, 2e-14);
        }
        assert!((a.viscous_heat_j - b.viscous_heat_j).abs() < 1e-20);
        assert!((a.added_mass_energy_change_j - b.added_mass_energy_change_j).abs() < 1e-20);
        still = next_still;
        moving = next_moving;
    }
    // The ordinary no-acceleration API holds the accepted current, not the
    // immutable reservoir's initial velocity and not a new quiescent history.
    let (held, held_audit) = moving
        .step_with_wrench(Wrench::ZERO, step_controls(dt), || false)
        .unwrap();
    let (explicit, explicit_audit) = moving
        .step_with_flow_acceleration(Wrench::ZERO, NativeVec3::ZERO, step_controls(dt), || false)
        .unwrap();
    assert_eq!(
        held.ambient_velocity_world_m_s(),
        moving.ambient_velocity_world_m_s()
    );
    assert_eq!(held.state(), explicit.state());
    assert_eq!(held_audit, explicit_audit);
}

#[test]
fn accelerating_flow_refusal_preserves_ambient_clock_body_and_history() {
    let source = sphere(2000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1000.0, 1.0);
    let medium =
        CreepingSphere::new(&material, RADIUS, Gravity::ZERO, NativeVec3::ZERO, limits()).unwrap();
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let motion = medium
        .unsteady_motion(&body, rest, resolution(0.001, 1000.0, 64), || false)
        .unwrap();
    let acceleration = NativeVec3::new(0.01, -0.02, 0.015);
    let (motion, _) = motion
        .step_with_flow_acceleration(Wrench::ZERO, acceleration, step_controls(0.001), || false)
        .unwrap();
    let mut polls = 0;
    let (expected, expected_audit) = motion
        .step_with_flow_acceleration(Wrench::ZERO, acceleration, step_controls(0.001), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        motion.step_with_flow_acceleration(
            Wrench::ZERO,
            acceleration,
            step_controls(0.001),
            || {
                visited += 1;
                visited == polls
            }
        ),
        Err(CreepingFlowError::Cancelled)
    ));
    for invalid in [
        NativeVec3::new(f64::NAN, 0.0, 0.0),
        NativeVec3::new(1e20, 0.0, 0.0),
        NativeVec3::new(1e-100, 0.0, 0.0),
    ] {
        assert!(
            motion
                .step_with_flow_acceleration(Wrench::ZERO, invalid, step_controls(0.001), || false)
                .is_err()
        );
    }
    let (retried, retried_audit) = motion
        .clone()
        .step_with_flow_acceleration(Wrench::ZERO, acceleration, step_controls(0.001), || false)
        .unwrap();
    assert_eq!(retried.state(), expected.state());
    assert_eq!(retried.time_s(), expected.time_s());
    assert_eq!(
        retried.ambient_velocity_world_m_s(),
        expected.ambient_velocity_world_m_s()
    );
    assert_eq!(
        retried.fluid_velocities_world_m_s(),
        expected.fluid_velocities_world_m_s()
    );
    assert_eq!(retried_audit, expected_audit);
    assert!(
        medium
            .unsteady_motion(&body, motion.state(), resolution(0.001, 1000.0, 64), || {
                false
            })
            .is_err()
    );
}
