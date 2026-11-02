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

#[test]
fn changing_source_density_preserves_mass_and_changes_native_waterline_and_motion() {
    let initial = fluid(1000.0);
    let updated = fluid(800.0);
    let tank = FiniteHydrostaticTank::new(
        &initial,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let held = pose(1.5);
    let (next, audit) = tank
        .rebind_source_at_pose(&updated, &mesh, held, || false)
        .unwrap();
    assert_eq!(
        next.liquid_mass_kg().to_bits(),
        tank.liquid_mass_kg().to_bits()
    );
    close(next.liquid_volume_m3(), 20.0);
    assert!(std::ptr::eq(tank.source(), &initial));
    assert!(std::ptr::eq(next.source(), &updated));
    close(audit.fill_before.surface_height_m, 7.0 / 6.0);
    close(audit.fill_after.surface_height_m, 1.5);
    close(audit.fill_after.displaced.volume_m3(), 4.0);
    assert_eq!(audit.liquid_volume_m3, [16.0, 20.0]);
    assert_eq!(audit.density_kg_m3, [1000.0, 800.0]);

    // Independently integrate liquid height over the rectangular bath minus
    // the prism: G = rho*g*(A*h²/2 - a*(h²-bottom²)/2).
    let potential = |rho: f64, height: f64| {
        rho * 10.0 * (8.0 * height * height - 2.0 * (height * height - 0.25))
    };
    let expected = potential(800.0, 1.5) - potential(1000.0, 7.0 / 6.0);
    close(audit.empty_tank_gravitational_energy_change_j, 20000.0);
    close(audit.displacement_potential_change_j, 16000.0 / 3.0);
    close(audit.liquid_gravitational_energy_change_j, expected);
    close(audit.liquid_gravitational_energy_change_j, 76000.0 / 3.0);

    let material = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = material.principal_body(1e-12, || false).unwrap();
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    close(
        next.load(&mesh, &body, state, || false)
            .unwrap()
            .wrench()
            .force_world
            .z,
        32000.0,
    );
    let step = next
        .step(
            &mesh,
            &body,
            state,
            HydrostaticStepControls {
                duration_s: 0.001,
                maximum_energy_defect_j: 1e-6,
            },
            || false,
        )
        .unwrap();
    // Native first kick/drift: (32000 N - 4000 kg*10 m/s²)/4000 kg = -2 m/s².
    close(
        body.world_from_body(step.dynamics.state_after)
            .unwrap()
            .translation()
            .z,
        1.5 - 0.001_f64.powi(2),
    );

    let (restored, reverse) = next
        .rebind_source_at_pose(&initial, &mesh, held, || false)
        .unwrap();
    assert_eq!(restored.liquid_mass_kg(), tank.liquid_mass_kg());
    assert_eq!(restored.liquid_volume_m3(), tank.liquid_volume_m3());
    assert_eq!(reverse.fill_after, audit.fill_before);
    close(reverse.liquid_gravitational_energy_change_j, -expected);
}

#[test]
fn density_update_includes_empty_bath_energy_even_when_the_body_is_dry() {
    let initial = fluid(1000.0);
    let updated = fluid(800.0);
    let mesh = cube();
    for floor_z in [0.0, 10000.0] {
        let tank = FiniteHydrostaticTank::new(
            &initial,
            gravity(),
            HydrostaticTankGeometry {
                world_from_tank: pose(floor_z),
                inside_dimensions_m: Vec3::new(4.0, 4.0, 6.0),
            },
            16.0,
            volume_controls(),
        )
        .unwrap();
        let (_, audit) = tank
            .rebind_source_at_pose(&updated, &mesh, pose(floor_z + 4.0), || false)
            .unwrap();
        assert_eq!(audit.fill_before.displaced.volume_m3(), 0.0);
        assert_eq!(audit.fill_after.displaced.volume_m3(), 0.0);
        assert_eq!(audit.displacement_potential_change_j, 0.0);
        close(audit.liquid_gravitational_energy_change_j, 20000.0);
        close(audit.fill_after.surface_height_m, 1.25);
    }
}

#[test]
fn unchanged_density_retains_inventory_and_fill_bits_and_freezes_the_new_source_once() {
    struct Source {
        calls: std::cell::Cell<usize>,
        viscosity: f64,
    }
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.calls.set(self.calls.get() + 1);
            FluidTransportState::try_new(997.3, 4000.0, 0.5, self.viscosity)
        }
    }
    let original = Source {
        calls: std::cell::Cell::new(0),
        viscosity: 0.001,
    };
    let source = Source {
        calls: std::cell::Cell::new(0),
        viscosity: 0.002,
    };
    let tank = FiniteHydrostaticTank::new(
        &original,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        15.123,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let expected = tank.fill_at_pose(&mesh, pose(1.5), || false).unwrap();
    let (next, audit) = tank
        .rebind_source_at_pose(&source, &mesh, pose(1.5), || false)
        .unwrap();
    assert_eq!(original.calls.get(), 1);
    assert_eq!(source.calls.get(), 1);
    assert_eq!(
        next.liquid_mass_kg().to_bits(),
        tank.liquid_mass_kg().to_bits()
    );
    assert_eq!(
        next.liquid_volume_m3().to_bits(),
        tank.liquid_volume_m3().to_bits()
    );
    assert_eq!(audit.fill_before, expected);
    assert_eq!(audit.fill_after, expected);
    assert_eq!(
        audit.empty_tank_gravitational_energy_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        audit.displacement_potential_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        audit.liquid_gravitational_energy_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(next.transport().dynamic_viscosity_pa_s(), 0.002);
    next.fill_at_pose(&mesh, pose(1.5), || false).unwrap();
    assert_eq!(source.calls.get(), 1);

    // Zero change remains representable even if the absolute M*g scale is not.
    let large = FiniteHydrostaticTank::new(
        &original,
        Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -1e10)).unwrap(),
        tank_geometry(Vec3::new(1e150, 1e150, 6.0)),
        1e300,
        volume_controls(),
    )
    .unwrap();
    let (_, unchanged) = large
        .rebind_source_at_pose(&source, &mesh, pose(4.0), || false)
        .unwrap();
    assert_eq!(unchanged.liquid_gravitational_energy_change_j, 0.0);
}

#[test]
fn density_update_source_overflow_and_late_cancellation_leave_the_original_tank_intact() {
    struct Source {
        density: f64,
    }
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            FluidTransportState::try_new(self.density, 4000.0, 0.5, 0.001)
        }
    }
    let initial = Source { density: 1000.0 };
    let changed = Source { density: 800.0 };
    let tank = FiniteHydrostaticTank::new(
        &initial,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 1.7)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let mesh = cube();
    let held = pose(1.5);
    let expected = tank.fill_at_pose(&mesh, held, || false).unwrap();
    // Empty inventory fits, but its displacement at the held pose overflows.
    assert!(matches!(
        tank.rebind_source_at_pose(&Source { density: 600.0 }, &mesh, held, || false),
        Err(HydrostaticError::TankOverflow { .. })
    ));
    // Even the empty inventory cannot fit at this lower density.
    assert!(matches!(
        tank.rebind_source_at_pose(&Source { density: 400.0 }, &mesh, held, || false),
        Err(HydrostaticError::TankOverflow { .. })
    ));
    assert!(matches!(
        tank.rebind_source_at_pose(&Source { density: -1.0 }, &mesh, held, || false),
        Err(HydrostaticError::Material(_))
    ));
    let mut polls = 0;
    let (_, success) = tank
        .rebind_source_at_pose(&changed, &mesh, held, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            tank.rebind_source_at_pose(&changed, &mesh, held, || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticError::Cancelled)
        ));
    }
    assert_eq!(tank.fill_at_pose(&mesh, held, || false).unwrap(), expected);
    assert_eq!(tank.liquid_mass_kg(), 16000.0);
    assert_eq!(tank.liquid_volume_m3(), 16.0);
    let (_, replay) = tank
        .rebind_source_at_pose(&changed, &mesh, held, || false)
        .unwrap();
    assert_eq!(replay, success);
}
