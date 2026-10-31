//! End-to-end dynamics, not just plausible hydrostatic coefficients.
use super::*;
use fs_mbd_material::hydrostatic::{HydrostaticStepControls, HydrostaticStepError};

fn controls(dt: f64) -> HydrostaticStepControls {
    HydrostaticStepControls {
        duration_s: dt,
        maximum_energy_defect_j: 0.1,
    }
}

#[test]
fn step_reuses_the_existing_integrator_and_recomputes_the_actual_endpoint_load() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(0.3), Twist::zero())
        .unwrap();
    let mesh = cube();
    let result = medium
        .step(&mesh, &body, state, controls(0.001), || false)
        .unwrap();
    let (direct, endpoint) = RigidBodyIntegrator::new(medium.gravity())
        .step_verlet_with(
            state,
            body.properties(),
            result.load_before.wrench(),
            0.001,
            |drifted| {
                let load = medium.load(&mesh, &body, drifted, || false)?;
                Ok::<_, HydrostaticError>((load.wrench(), load))
            },
        )
        .unwrap();
    assert_eq!(result.dynamics, direct);
    assert_eq!(result.load_after, endpoint);
    assert_eq!(
        result.load_after,
        medium
            .load(&mesh, &body, direct.state_after, || false)
            .unwrap()
    );
    assert!(result.load_after.displaced().volume_m3() > result.load_before.displaced().volume_m3());
    assert!(
        ((result.total_energy_after_j - result.total_energy_before_j) - result.energy_defect_j)
            .abs()
            < 1e-8
    );
    assert_eq!(state, result.dynamics.state_before);
}

#[test]
fn floating_draft_depends_on_body_and_fluid_data_without_mass_override() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    for density in [250.0, 500.0, 750.0] {
        let source = body(density, sides);
        let body = source.principal_body(1e-12, || false).unwrap();
        let draft = density / 1000.0;
        let initial = body
            .state_from_kinematics(pose(0.5 - draft), Twist::zero())
            .unwrap();
        let mut state = initial;
        for _ in 0..100 {
            let result = medium
                .step(&mesh, &body, state, controls(0.005), || false)
                .unwrap();
            close(
                result.load_after.displaced().volume_m3() / mesh.volume_m3(),
                draft,
            );
            state = result.dynamics.state_after;
        }
        assert!(
            (state.pose().position_world().z - initial.pose().position_world().z).abs() < 1e-10
        );
        assert!(
            state
                .center_of_mass_velocity_world(body.properties())
                .unwrap()
                .norm_squared()
                < 1e-18
        );
        close(body.properties().mass(), density * 4.0);
    }
}

fn heave(dt: f64, steps: usize) -> (f64, f64, f64) {
    // Squat geometry is roll-stable, unlike an upright half-submerged cube.
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let mut state = body
        .state_from_kinematics(pose(0.02), Twist::zero())
        .unwrap();
    let mut cumulative = 0.0;
    let mut maximum = 0.0_f64;
    for _ in 0..steps {
        let result = medium
            .step(&mesh, &body, state, controls(dt), || false)
            .unwrap();
        cumulative += result.energy_defect_j;
        maximum = maximum.max(cumulative.abs());
        state = result.dynamics.state_after;
        // Independent upright-cuboid law, not a re-reading of load_after.
        let height = body.world_from_body(state).unwrap().translation().z;
        assert!((result.load_after.displaced().volume_m3() - 4.0 * (0.5 - height)).abs() < 1e-8);
    }
    let position = body.world_from_body(state).unwrap().translation().z;
    (
        position,
        state
            .center_of_mass_velocity_world(body.properties())
            .unwrap()
            .z,
        maximum,
    )
}

#[test]
fn hydrostatic_heave_converges_to_the_independent_continuum_oscillator() {
    // m=2000 kg; rho*g*waterplane_area=40000 N/m; no added mass in this model.
    let exact = 0.02 * (20.0_f64.sqrt() * 1.3).cos();
    let mut errors = Vec::new();
    for (dt, n) in [(0.02, 65), (0.01, 130), (0.005, 260)] {
        let (position, _, energy) = heave(dt, n);
        errors.push((position - exact).abs());
        assert!(energy < 0.4);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.8..4.2).contains(&ratio), "{errors:?}");
    }
}

#[test]
fn repeated_geometry_feedback_does_not_invent_hydrodynamic_damping() {
    let (position, velocity, maximum) = heave(0.002, 10_000);
    // Initial oscillation energy is 8 J above equilibrium. These bounds check
    // this trajectory and discretization, not a general stability certificate.
    let final_energy = 1000.0 * velocity * velocity + 20000.0 * position * position;
    assert!((final_energy - 8.0).abs() < 0.04);
    assert!(maximum < 0.04);
}

#[test]
fn floating_roll_has_second_order_convergence_to_the_metacentric_oscillator() {
    // Half-submerged 3 x 2 x 1 m cuboid: GM=B²/(6H)-H/4=5/12 m,
    // I_x/m=(B²+H²)/12=5/12 m². Thus omega²=g*GM/(I_x/m)=10 s^-2.
    // The 1e-4-rad amplitude leaves O(amplitude³) nonlinear error below the
    // integration errors here. The oracle does not use the clipping routine.
    let amplitude = 1e-4;
    let exact = amplitude * (10.0_f64.sqrt() * 1.3).cos();
    let sides = Vec3::new(3.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let initial_pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(amplitude, 0.0, 0.0))).unwrap(),
        origin(),
    )
    .unwrap();
    let initial = body
        .state_from_kinematics(initial_pose, Twist::zero())
        .unwrap();
    let mut errors = Vec::new();
    for (dt, steps) in [(0.02, 65), (0.01, 130), (0.005, 260)] {
        let mut current = initial;
        for _ in 0..steps {
            current = medium
                .step(&mesh, &body, current, controls(dt), || false)
                .unwrap()
                .dynamics
                .state_after;
        }
        let y = body
            .world_from_body(current)
            .unwrap()
            .rotation()
            .rotate(Vec3::new(0.0, 1.0, 0.0))
            .unwrap();
        errors.push((y.z.atan2(y.y) - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn tumbling_hull_retains_yaw_momentum_and_both_endpoint_pressure_moments() {
    let sides = Vec3::new(1.7, 2.1, 0.8);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.17, 0.23, -0.11))).unwrap(),
        Vec3::new(0.4, -0.2, 0.02),
    )
    .unwrap();
    let mut current = body
        .state_from_kinematics(
            pose,
            Twist::new(Vec3::new(0.3, -0.2, 0.4), Vec3::new(0.07, -0.05, 0.02)),
        )
        .unwrap();
    let integrator = RigidBodyIntegrator::new(gravity());
    let initial_yaw = integrator
        .diagnostics(current, body.properties())
        .unwrap()
        .angular_momentum_world
        .z;
    let weight = gravity()
        .acceleration_world()
        .scale(body.properties().mass());
    let dt = 0.005;
    for _ in 0..200 {
        let step = medium
            .step(&mesh, &body, current, controls(dt), || false)
            .unwrap();
        let after = step.dynamics.state_after;
        let moment = |state: RigidBodyState,
                      load: fs_mbd_material::hydrostatic::HydrostaticLoad| {
            let center = load.center_of_buoyancy_world().unwrap();
            let center = fs_mbd::Vec3::new(center.x, center.y, center.z);
            center
                .cross(load.wrench().force_world)
                .add(state.pose().position_world().cross(weight))
        };
        let expected = moment(current, step.load_before)
            .add(moment(after, step.load_after))
            .scale(0.5 * dt);
        let actual = step
            .dynamics
            .diagnostics_after
            .angular_momentum_world
            .sub(step.dynamics.diagnostics_before.angular_momentum_world);
        assert!(actual.sub(expected).norm_squared().sqrt() < 3e-10);
        let impulse = step
            .load_before
            .wrench()
            .force_world
            .add(step.load_after.wrench().force_world)
            .scale(0.5 * dt)
            .add(weight.scale(dt));
        assert!(
            after
                .linear_momentum_world()
                .sub(current.linear_momentum_world())
                .sub(impulse)
                .norm_squared()
                .sqrt()
                < 3e-10
        );
        assert!(
            (step.dynamics.diagnostics_after.angular_momentum_world.z - initial_yaw).abs() < 3e-10
        );
        current = after;
    }
}

#[test]
fn overlarge_steps_and_late_cancellation_publish_no_candidate() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(0.3), Twist::zero())
        .unwrap();
    let saved = state;
    let mesh = cube();
    assert!(matches!(
        medium.step(&mesh, &body, state, controls(0.2), || false),
        Err(HydrostaticStepError::EnergyDefect { .. })
    ));
    let mut count = 0;
    let accepted = medium
        .step(&mesh, &body, state, controls(0.001), || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        medium.step(&mesh, &body, state, controls(0.001), || {
            calls += 1;
            calls == count
        }),
        Err(HydrostaticStepError::Physics(HydrostaticError::Cancelled))
    ));
    assert_eq!(state, saved);
    assert_eq!(
        medium
            .step(&mesh, &body, state, controls(0.001), || false)
            .unwrap(),
        accepted
    );
    for c in [
        HydrostaticStepControls {
            duration_s: 0.0,
            ..controls(0.001)
        },
        HydrostaticStepControls {
            duration_s: f64::NAN,
            ..controls(0.001)
        },
        HydrostaticStepControls {
            maximum_energy_defect_j: f64::INFINITY,
            ..controls(0.001)
        },
        HydrostaticStepControls {
            maximum_energy_defect_j: 0.0,
            ..controls(0.001)
        },
    ] {
        assert!(medium.step(&mesh, &body, state, c, || false).is_err());
    }
}

#[test]
fn free_translation_parallel_to_the_surface_does_not_change_hydrostatic_energy() {
    let source = body(500.0, Vec3::new(2.0, 2.0, 1.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let mesh = DisplacementMesh::cuboid(Vec3::new(2.0, 2.0, 1.0), || false).unwrap();
    let state = body
        .state_from_kinematics(
            Se3::identity(),
            Twist::new(origin(), Vec3::new(2.0, -1.0, 0.0)),
        )
        .unwrap();
    let result = medium
        .step(&mesh, &body, state, controls(0.1), || false)
        .unwrap();
    close(result.dynamics.state_after.pose().position_world().x, 0.2);
    close(result.dynamics.state_after.pose().position_world().y, -0.1);
    close(result.energy_defect_j, 0.0);
}
