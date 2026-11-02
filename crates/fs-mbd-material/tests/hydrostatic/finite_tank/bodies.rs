//! Shared liquid inventory, cross-body forces and native coupled motion.
use super::*;
use fs_mbd_material::hydrostatic::{FiniteTankBody, HydrostaticStepError, MAX_TANK_BODIES};

fn at(x: f64, z: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0.0, z)).unwrap()
}

#[test]
fn two_bodies_share_one_waterline_and_the_full_fluid_energy_cross_term() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let z = [17.0 / 16.0 + 0.03, 17.0 / 16.0 - 0.01];
    let states = [
        body.state_from_kinematics(at(-1.0, z[0]), Twist::zero())
            .unwrap(),
        body.state_from_kinematics(at(1.0, z[1]), Twist::zero())
            .unwrap(),
    ];
    let load = tank.load_bodies(&bindings, &states, || false).unwrap();
    // A=16, each body waterplane a=1: (16-2)*h=16+1-z0-z1.
    let h = (17.0 - z[0] - z[1]) / 14.0;
    let d = [h - z[0] + 0.5, h - z[1] + 0.5];
    close(load.surface.surface_height_m, h);
    close(load.surface.displaced_volume_m3, d[0] + d[1]);
    close(16.0 * h - d[0] - d[1], tank.liquid_volume_m3());
    close(tank.liquid_mass_kg(), 16000.0);
    for i in 0..2 {
        close(load.bodies[i].displaced_mass_kg(), 1000.0 * d[i]);
        close(load.bodies[i].wrench().force_world.z, 10000.0 * d[i]);
        close(
            load.bodies[i].center_of_buoyancy_world().unwrap().z,
            (h + z[i] - 0.5) / 2.0,
        );
        let isolated = tank.load(&mesh, &body, states[i], || false).unwrap();
        assert!(load.bodies[i].wrench().force_world.z > isolated.wrench().force_world.z);
    }
    // Independently integrate the bath first moment minus each submerged slab,
    // relative to the same empty liquid at h0=1. No per-body fluid copies.
    let expected = 10000.0
        * (8.0 * h * h
            - 8.0
            - 0.5 * (h * h - (z[0] - 0.5).powi(2))
            - 0.5 * (h * h - (z[1] - 0.5).powi(2)));
    close(load.surface.liquid_gravitational_energy_change_j, expected);
    let moved = [
        states[0],
        body.state_from_kinematics(at(1.0, z[1] + 0.02), Twist::zero())
            .unwrap(),
    ];
    let after = tank.load_bodies(&bindings, &moved, || false).unwrap();
    close(
        after.bodies[0].wrench().force_world.z - load.bodies[0].wrench().force_world.z,
        -10000.0 * 0.02 / 14.0,
    );
}

#[test]
fn shared_potential_gradient_matches_each_original_body_force_and_torque() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sizes = [Vec3::new(0.8, 1.2, 1.0), Vec3::new(1.0, 0.8, 1.3)];
    let sources = [body(500.0, sizes[0]), body(700.0, sizes[1])];
    let bodies = sources
        .each_ref()
        .map(|source| source.principal_body(1e-12, || false).unwrap());
    let meshes = sizes.map(|size| DisplacementMesh::cuboid(size, || false).unwrap());
    let bindings = [
        FiniteTankBody {
            envelope: &meshes[0],
            frame: bodies[0].frame(),
        },
        FiniteTankBody {
            envelope: &meshes[1],
            frame: bodies[1].frame(),
        },
    ];
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.15, 0.2, -0.1))).unwrap();
    let states = [0, 1].map(|i| {
        bodies[i]
            .state_from_kinematics(
                Se3::from_parts(
                    rotation,
                    Vec3::new(if i == 0 { -1.0 } else { 1.0 }, 0.0, 1.05),
                )
                .unwrap(),
                Twist::zero(),
            )
            .unwrap()
    });
    let load = tank.load_bodies(&bindings, &states, || false).unwrap();
    let h = 1e-5;
    for i in 0..2 {
        for unit in [
            fs_mbd::Vec3::new(1.0, 0.0, 0.0),
            fs_mbd::Vec3::new(0.0, 1.0, 0.0),
            fs_mbd::Vec3::new(0.0, 0.0, 1.0),
        ] {
            for rotational in [false, true] {
                let energy = |sign: f64| {
                    let old = states[i].pose();
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
                    let mut changed = states;
                    changed[i] = RigidBodyState::new(
                        Pose::new(position, orientation).unwrap(),
                        states[i].linear_momentum_world(),
                        states[i].angular_momentum_body(),
                    )
                    .unwrap();
                    tank.load_bodies(&bindings, &changed, || false)
                        .unwrap()
                        .surface
                        .liquid_gravitational_energy_change_j
                };
                let wrench = if rotational {
                    load.bodies[i].wrench().torque_body
                } else {
                    load.bodies[i].wrench().force_world
                };
                let derivative = -(energy(1.0) - energy(-1.0)) / (2.0 * h);
                assert!(
                    (derivative - wrench.dot(unit)).abs() < 4e-4,
                    "body={i}, rotational={rotational}, derivative={derivative}, wrench={wrench:?}"
                );
            }
        }
    }
}

fn coupled_heave(dt: f64, steps: usize) -> [f64; 2] {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let equilibrium = 17.0 / 16.0;
    let mut states = [
        body.state_from_kinematics(at(-1.0, equilibrium + 0.02), Twist::zero())
            .unwrap(),
        body.state_from_kinematics(at(1.0, equilibrium), Twist::zero())
            .unwrap(),
    ];
    let mut defect = 0.0;
    for _ in 0..steps {
        let next = tank
            .step_bodies(
                &bindings,
                &states,
                HydrostaticStepControls {
                    duration_s: dt,
                    maximum_energy_defect_j: 0.05,
                },
                || false,
            )
            .unwrap();
        defect += next.energy_defect_j;
        states = [next.dynamics[0].state_after, next.dynamics[1].state_after];
        let z = states.map(|state| body.world_from_body(state).unwrap().translation().z);
        close(
            next.load_after.surface.surface_height_m,
            (17.0 - z[0] - z[1]) / 14.0,
        );
        assert!(next.load_after.surface.liquid_volume_defect_m3.abs() <= 1e-12);
    }
    assert!(defect.abs() < 0.05);
    states.map(|state| body.world_from_body(state).unwrap().translation().z)
}

#[test]
fn synchronous_native_motion_converges_to_both_collective_heave_modes() {
    // K_ij=rho*g*a_i*(delta_ij+a_j/(A-sum(a))), M=500 each.
    // In-phase omega²=160/7; anti-phase omega²=20. Initially only body 0 moves.
    let in_phase = ((160.0_f64 / 7.0).sqrt() * 1.2).cos();
    let anti_phase = (20.0_f64.sqrt() * 1.2).cos();
    let exact = [
        17.0 / 16.0 + 0.01 * (in_phase + anti_phase),
        17.0 / 16.0 + 0.01 * (in_phase - anti_phase),
    ];
    let errors: Vec<_> = [(0.04, 30), (0.02, 60), (0.01, 120)]
        .into_iter()
        .map(|(dt, steps)| {
            let actual = coupled_heave(dt, steps);
            ((actual[0] - exact[0]).powi(2) + (actual[1] - exact[1]).powi(2)).sqrt()
        })
        .collect();
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn singleton_batches_preserve_original_native_load_step_and_force_admission() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let binding = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }];
    for z in [3.0, 1.04, 0.55] {
        let state = body
            .state_from_kinematics(at(0.0, z), Twist::zero())
            .unwrap();
        let controls = HydrostaticStepControls {
            duration_s: 0.001,
            maximum_energy_defect_j: 0.1,
        };
        let single = tank.step(&mesh, &body, state, controls, || false).unwrap();
        let shared = tank
            .step_bodies(&binding, &[state], controls, || false)
            .unwrap();
        assert_eq!(shared.dynamics[0], single.dynamics);
        assert_eq!(shared.load_before.bodies[0], single.load_before.buoyancy());
        assert_eq!(shared.load_after.bodies[0], single.load_after.buoyancy());
        assert_eq!(shared.total_energy_before_j, single.total_energy_before_j);
        assert_eq!(shared.total_energy_after_j, single.total_energy_after_j);
        assert_eq!(shared.energy_defect_j, single.energy_defect_j);
        let stage = RigidBodyIntegrator::new(gravity())
            .begin_verlet(
                state,
                body.frame().properties(),
                single.load_before.wrench(),
                0.001,
            )
            .unwrap();
        assert!(
            stage
                .finish(fs_mbd::Wrench {
                    force_world: fs_mbd::Vec3::new(f64::NAN, 0.0, 0.0),
                    torque_body: fs_mbd::Vec3::ZERO,
                })
                .is_err()
        );
    }
}

#[test]
fn collective_overflow_and_late_step_failure_publish_no_body_or_fluid_history() {
    struct Source(std::cell::Cell<usize>);
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get() + 1);
            FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001)
        }
    }
    let water = Source(std::cell::Cell::new(0));
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let binding = FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    };
    let bindings = [binding; 2];
    let states = [
        body.state_from_kinematics(at(-1.0, 17.0 / 16.0 + 0.02), Twist::zero())
            .unwrap(),
        body.state_from_kinematics(at(1.0, 17.0 / 16.0), Twist::zero())
            .unwrap(),
    ];
    let controls = HydrostaticStepControls {
        duration_s: 0.01,
        maximum_energy_defect_j: 0.05,
    };
    let mut polls = 0;
    let expected = tank
        .step_bodies(&bindings, &states, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            tank.step_bodies(&bindings, &states, controls, || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticStepError::Physics(HydrostaticError::Cancelled))
        ));
    }
    assert!(matches!(
        tank.step_bodies(
            &bindings,
            &states,
            HydrostaticStepControls {
                duration_s: 0.2,
                maximum_energy_defect_j: 1e-12,
            },
            || false
        ),
        Err(HydrostaticStepError::EnergyDefect { .. })
    ));
    assert_eq!(
        tank.step_bodies(&bindings, &states, controls, || false)
            .unwrap(),
        expected
    );
    assert_eq!(water.0.get(), 1);
    assert!(tank.load_bodies(&bindings, &states[..1], || false).is_err());
    assert!(
        tank.load_bodies(
            &vec![binding; MAX_TANK_BODIES + 1],
            &vec![states[0]; MAX_TANK_BODIES + 1],
            || false
        )
        .is_err()
    );

    let shallow = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 1.04)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    for state in states {
        shallow.load(&mesh, &body, state, || false).unwrap();
    }
    assert!(matches!(
        shallow.load_bodies(&bindings, &states, || false),
        Err(HydrostaticError::TankOverflow { .. })
    ));
}

#[test]
fn shared_density_update_matches_independent_liquid_energy_and_continues_native_motion() {
    let water = fluid(1000.0);
    let changed = fluid(800.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let z = [1.1_f64, 1.05];
    let poses = [at(-1.0, z[0]), at(1.0, z[1])];
    let states = poses.map(|pose| body.state_from_kinematics(pose, Twist::zero()).unwrap());
    let old_load = tank.load_bodies(&bindings, &states, || false).unwrap();
    let (next, audit) = tank
        .rebind_source_at_bodies(&changed, &bindings, &states, || false)
        .unwrap();
    assert_eq!(audit.world_from_body, poses);
    assert_eq!(audit.liquid_volume_m3, [16.0, 20.0]);
    assert_eq!(
        next.liquid_mass_kg().to_bits(),
        tank.liquid_mass_kg().to_bits()
    );
    let h = [16.0, 20.0].map(|v| (v + 1.0 - z[0] - z[1]) / 14.0);
    close(audit.fill_before.surface_height_m, h[0]);
    close(audit.fill_after.surface_height_m, h[1]);
    // Direct integrate rho*g*z over the actual bath minus its two body slabs.
    let energy = |rho: f64, h: f64| {
        10.0 * rho * (8.0 * h * h - 0.5 * z.iter().map(|z| h * h - (z - 0.5).powi(2)).sum::<f64>())
    };
    close(
        audit.liquid_gravitational_energy_change_j,
        energy(800.0, h[1]) - energy(1000.0, h[0]),
    );
    close(audit.empty_tank_gravitational_energy_change_j, 20000.0);
    let load = next.load_bodies(&bindings, &states, || false).unwrap();
    assert_eq!(load.surface, audit.fill_after);
    assert_eq!(old_load.surface, audit.fill_before);
    let dt = 0.001_f64;
    let step = next
        .step_bodies(
            &bindings,
            &states,
            HydrostaticStepControls {
                duration_s: dt,
                maximum_energy_defect_j: 1e-5,
            },
            || false,
        )
        .unwrap();
    for i in 0..2 {
        let mass = 800.0 * (h[1] - z[i] + 0.5);
        close(load.bodies[i].displaced_mass_kg(), mass);
        close(
            body.world_from_body(step.dynamics[i].state_after)
                .unwrap()
                .translation()
                .z,
            z[i] + 0.5 * dt * dt * 10.0 * (mass / 500.0 - 1.0),
        );
    }
    let (_, reverse) = next
        .rebind_source_at_bodies(&water, &bindings, &states, || false)
        .unwrap();
    assert_eq!(reverse.fill_after, audit.fill_before);
    close(
        reverse.liquid_gravitational_energy_change_j,
        -audit.liquid_gravitational_energy_change_j,
    );
    assert_eq!(
        tank.load_bodies(&bindings, &states, || false).unwrap(),
        old_load
    );
}

#[test]
fn shared_source_refresh_preserves_singleton_bits_and_refuses_collective_overflow_atomically() {
    let water = fluid(1000.0);
    let changed = fluid(800.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 4.0)),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let poses = [at(-1.0, 1.0), at(1.0, 1.0)];
    let states = poses.map(|pose| body.state_from_kinematics(pose, Twist::zero()).unwrap());
    let (_, a) = tank
        .rebind_source_at_bodies(&changed, &bindings[..1], &states[..1], || false)
        .unwrap();
    let (_, b) = tank
        .rebind_source_at_pose(&changed, &mesh, poses[0], || false)
        .unwrap();
    assert_eq!(
        a.liquid_gravitational_energy_change_j,
        b.liquid_gravitational_energy_change_j
    );
    assert_eq!(
        a.displacement_potential_change_j,
        b.displacement_potential_change_j
    );
    assert_eq!(a.fill_after.surface_height_m, b.fill_after.surface_height_m);
    assert_eq!(
        a.fill_after.displaced_volume_m3,
        b.fill_after.displaced.volume_m3()
    );
    assert_eq!(
        a.fill_after.volume_evaluations,
        b.fill_after.volume_evaluations
    );
    let (_, unchanged) = tank
        .rebind_source_at_bodies(&water, &bindings, &states, || false)
        .unwrap();
    assert_eq!(unchanged.fill_before, unchanged.fill_after);
    assert_eq!(
        unchanged.liquid_gravitational_energy_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    let original = tank.load_bodies(&bindings, &states, || false).unwrap();
    let mut polls = 0;
    tank.rebind_source_at_bodies(&changed, &bindings, &states, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            tank.rebind_source_at_bodies(&changed, &bindings, &states, || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticError::Cancelled)
        ));
    }
    assert!(
        tank.rebind_source_at_bodies(&changed, &bindings, &states[..1], || false)
            .is_err()
    );
    assert!(
        tank.rebind_source_at_bodies(&changed, &[], &[], || false)
            .is_err()
    );
    assert_eq!(
        tank.load_bodies(&bindings, &states, || false).unwrap(),
        original
    );
    let shallow = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(Vec3::new(4.0, 4.0, 2.0)),
        28.0,
        volume_controls(),
    )
    .unwrap();
    let expanded = fluid(28000.0 / 30.5);
    for pose in poses {
        shallow
            .rebind_source_at_pose(&expanded, &mesh, pose, || false)
            .unwrap();
    }
    assert!(matches!(
        shallow.rebind_source_at_bodies(&expanded, &bindings, &states, || false),
        Err(HydrostaticError::TankOverflow { .. })
    ));
    assert_eq!(shallow.liquid_mass_kg(), 28000.0);
    assert_eq!(shallow.liquid_volume_m3(), 28.0);
}
