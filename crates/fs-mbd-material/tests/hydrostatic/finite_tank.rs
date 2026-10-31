//! G1/G3/G4 finite-volume buoyancy through the actual material/dynamics APIs.
use super::*;
use fs_mbd_material::hydrostatic::{
    FiniteHydrostaticTank, FiniteTankControls, HydrostaticStepControls, HydrostaticTankGeometry,
};

fn volume_controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-12,
        maximum_volume_evaluations: 64,
    }
}

fn tank_geometry(sides: Vec3) -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: sides,
    }
}

#[test]
fn finite_inventory_changes_waterline_load_and_actual_liquid_potential() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let current_pose = pose(1.5);
    let state = body
        .state_from_kinematics(current_pose, Twist::zero())
        .unwrap();
    let fill = tank.fill_at_pose(&mesh, current_pose, || false).unwrap();
    // 16*h - 4*(h-0.5) = 16: h=7/6, displaced volume=8/3.
    close(fill.surface_height_m, 7.0 / 6.0);
    close(fill.displaced.volume_m3(), 8.0 / 3.0);
    close(tank.liquid_mass_kg(), 16000.0);
    close(tank.liquid_volume_m3(), 16.0);
    assert!(fill.liquid_volume_defect_m3.abs() <= 1e-12);
    assert!(fill.volume_evaluations <= 3);
    let load = tank.load(&mesh, &body, state, || false).unwrap();
    close(load.wrench().force_world.z, 80000.0 / 3.0);
    close(load.center_of_buoyancy_world().unwrap().z, 5.0 / 6.0);
    // Independent tank-minus-submerged-body first moment about the floor.
    let direct_liquid_energy_change =
        10000.0 * (8.0 * fill.surface_height_m.powi(2) - (8.0 / 3.0) * (5.0 / 6.0) - 8.0);
    close(load.pressure_potential_j(), direct_liquid_energy_change);
    close(load.pressure_potential_j(), 20000.0 / 3.0);
    close(
        fill.liquid_gravitational_energy_change_j,
        load.pressure_potential_j(),
    );
    // The old prescribed-height model is still available and is physically
    // different: it displaces only 2 m3 at the unchanged h=1 surface.
    let infinite = HydrostaticMedium::new(&water, gravity(), Vec3::new(0.0, 0.0, 1.0)).unwrap();
    close(
        infinite
            .load(&mesh, &body, state, || false)
            .unwrap()
            .wrench()
            .force_world
            .z,
        20000.0,
    );
    assert!(std::ptr::eq(tank.source(), &water));
}

#[test]
fn dry_and_fully_submerged_bodies_use_the_same_conserved_inventory() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 5.0)),
        32.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let source = body(700.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    for z in [1.0, 1.2] {
        let fill = tank.fill_at_pose(&mesh, pose(z), || false).unwrap();
        close(fill.surface_height_m, 2.25);
        close(fill.displaced.volume_m3(), 4.0);
        close(
            fill.liquid_gravitational_energy_change_j,
            10000.0 * (4.0 * (2.25 - z) - 0.5),
        );
        let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
        close(
            tank.load(&mesh, &body, state, || false)
                .unwrap()
                .wrench()
                .force_world
                .z,
            40000.0,
        );
    }
    let dry = tank.fill_at_pose(&mesh, pose(4.0), || false).unwrap();
    close(dry.surface_height_m, 2.0);
    close(dry.displaced.volume_m3(), 0.0);
    close(dry.liquid_gravitational_energy_change_j, 0.0);
    assert_eq!(dry.volume_evaluations, 1);
}

#[test]
fn finite_bath_potential_gradient_matches_force_and_righting_torque() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(5.0, 5.0, 4.0)),
        25.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.2, 1.8, 0.8);
    let source = body(800.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.15, 0.3, 0.1))).unwrap();
    let pose = Se3::from_parts(rotation, Vec3::new(0.1, -0.2, 1.0)).unwrap();
    let state = body.state_from_kinematics(pose, Twist::zero()).unwrap();
    let load = tank.load(&mesh, &body, state, || false).unwrap();
    let h = 1e-5;
    for axis in 0..3 {
        let unit = match axis {
            0 => fs_mbd::Vec3::new(1.0, 0.0, 0.0),
            1 => fs_mbd::Vec3::new(0.0, 1.0, 0.0),
            _ => fs_mbd::Vec3::new(0.0, 0.0, 1.0),
        };
        for rotational in [false, true] {
            let energy = |sign: f64| {
                let old = state.pose();
                let position = if rotational {
                    old.position_world()
                } else {
                    old.position_world().add(unit.scale(sign * h))
                };
                let orientation = if rotational {
                    old.orientation().right_exp(unit.scale(sign * h)).unwrap()
                } else {
                    old.orientation()
                };
                let changed = RigidBodyState::new(
                    Pose::new(position, orientation).unwrap(),
                    state.linear_momentum_world(),
                    state.angular_momentum_body(),
                )
                .unwrap();
                tank.load(&mesh, &body, changed, || false)
                    .unwrap()
                    .pressure_potential_j()
            };
            let wrench = if rotational {
                load.wrench().torque_body
            } else {
                load.wrench().force_world
            };
            let derivative = -(energy(1.0) - energy(-1.0)) / (2.0 * h);
            assert!(
                (derivative - wrench.dot(unit)).abs() < 2e-4,
                "axis={axis}, rotational={rotational}, derivative={derivative}, load={wrench:?}"
            );
        }
    }
}

fn finite_heave(dt: f64, steps: usize) -> f64 {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let mut state = body
        .state_from_kinematics(pose(1.125 + 0.02), Twist::zero())
        .unwrap();
    let mut defect = 0.0;
    for _ in 0..steps {
        let next = tank
            .step(
                &mesh,
                &body,
                state,
                HydrostaticStepControls {
                    duration_s: dt,
                    maximum_energy_defect_j: 0.1,
                },
                || false,
            )
            .unwrap();
        defect += next.energy_defect_j;
        state = next.dynamics.state_after;
        let z = body.world_from_body(state).unwrap().translation().z;
        let fill = tank
            .fill_at_pose(&mesh, body.world_from_body(state).unwrap(), || false)
            .unwrap();
        close(fill.surface_height_m, 1.5 - z / 3.0);
        close(
            next.load_after.displaced().volume_m3(),
            8.0 - 16.0 * z / 3.0,
        );
    }
    assert!(defect.abs() < 0.5);
    body.world_from_body(state).unwrap().translation().z
}

#[test]
fn actual_native_motion_resolves_the_finite_bath_heave_stiffness() {
    // Bath A=16, body waterplane a=4, m=2000, rho*g=10000:
    // k=rho*g*a*A/(A-a), omega^2=80/3. The fixed-height answer is 20.
    let exact = 1.125 + 0.02 * ((80.0_f64 / 3.0).sqrt() * 1.3).cos();
    let mut errors = Vec::new();
    for (dt, steps) in [(0.02, 65), (0.01, 130), (0.005, 260)] {
        errors.push((finite_heave(dt, steps) - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn arbitrary_world_frame_preserves_inventory_and_native_body_torque() {
    let water = fluid(900.0);
    let sides = Vec3::new(1.2, 1.8, 0.8);
    let source = body(600.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let local_pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.2, 0.1, -0.15))).unwrap(),
        Vec3::new(0.1, -0.1, 1.0),
    )
    .unwrap();
    let first = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 5.0, 4.0)),
        20.0,
        volume_controls(),
    )
    .unwrap();
    let global = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(-0.3, 0.4, 0.1))).unwrap(),
        Vec3::new(3.0, -7.0, 11.0),
    )
    .unwrap();
    let rotated_g = global
        .rotation()
        .rotate(Vec3::new(0.0, 0.0, -10.0))
        .unwrap();
    let second = FiniteHydrostaticTank::new(
        &water,
        Gravity::new(fs_mbd::Vec3::new(rotated_g.x, rotated_g.y, rotated_g.z)).unwrap(),
        HydrostaticTankGeometry {
            world_from_tank: global,
            inside_dimensions_m: Vec3::new(4.0, 5.0, 4.0),
        },
        20.0,
        volume_controls(),
    )
    .unwrap();
    let state = body
        .state_from_kinematics(local_pose, Twist::zero())
        .unwrap();
    let moved = body
        .state_from_kinematics(global.compose(local_pose).unwrap(), Twist::zero())
        .unwrap();
    let a = first.load(&mesh, &body, state, || false).unwrap();
    let b = second.load(&mesh, &body, moved, || false).unwrap();
    close(a.displaced().volume_m3(), b.displaced().volume_m3());
    close(a.pressure_potential_j(), b.pressure_potential_j());
    let force = global
        .rotation()
        .rotate(Vec3::new(
            a.wrench().force_world.x,
            a.wrench().force_world.y,
            a.wrench().force_world.z,
        ))
        .unwrap();
    vec_close(
        force,
        Vec3::new(
            b.wrench().force_world.x,
            b.wrench().force_world.y,
            b.wrench().force_world.z,
        ),
    );
    close(a.wrench().torque_body.x, b.wrench().torque_body.x);
    close(a.wrench().torque_body.y, b.wrench().torque_body.y);
    close(a.wrench().torque_body.z, b.wrench().torque_body.z);
}

#[test]
fn overflow_contact_work_budget_and_cancellation_never_return_partial_loads() {
    let water = fluid(1000.0);
    let mesh = cube();
    let narrow = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 1.1)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    match narrow.fill_at_pose(&mesh, pose(1.5), || false) {
        Err(HydrostaticError::TankOverflow { excess_volume_m3 }) => close(excess_volume_m3, 0.8),
        other => panic!("expected overflow, got {other:?}"),
    }
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    assert!(tank.fill_at_pose(&mesh, pose(0.5), || false).is_err());
    let across_side = Se3::from_parts(So3::identity(), Vec3::new(1.1, 0.0, 1.5)).unwrap();
    assert!(tank.fill_at_pose(&mesh, across_side, || false).is_err());
    let budget = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        FiniteTankControls {
            maximum_volume_evaluations: 2,
            ..volume_controls()
        },
    )
    .unwrap();
    assert!(matches!(
        budget.fill_at_pose(&mesh, pose(1.5), || false),
        Err(HydrostaticError::TankVolumeNotConverged { evaluations: 2, .. })
    ));
    let mut calls = 0;
    let expected = tank
        .fill_at_pose(&mesh, pose(1.5), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in 1..=calls {
        let mut count = 0;
        assert_eq!(
            tank.fill_at_pose(&mesh, pose(1.5), || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticError::Cancelled)
        );
    }
    assert_eq!(
        tank.fill_at_pose(&mesh, pose(1.5), || false).unwrap(),
        expected
    );
    let tilted = HydrostaticTankGeometry {
        world_from_tank: Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.1, 0.0, 0.0))).unwrap(),
            origin(),
        )
        .unwrap(),
        ..tank_geometry(Vec3::new(4.0, 4.0, 4.0))
    };
    assert!(
        FiniteHydrostaticTank::new(&water, gravity(), tilted, 16.0, volume_controls()).is_err()
    );
}
