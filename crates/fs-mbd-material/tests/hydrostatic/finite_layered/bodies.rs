//! G1/G3/G4: density contrast couples actual bodies through conserved layers.
use super::*;
use fs_mbd_material::hydrostatic::{FiniteTankBody, MAX_TANK_BODIES};

fn at(x: f64, z: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0.0, z)).unwrap()
}
fn stack<'a>(
    bottom: &'a dyn FluidTransport,
    top: &'a dyn FluidTransport,
    volumes: [f64; 2],
    height: f64,
) -> FiniteLayeredHydrostaticTank<'a> {
    FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: bottom,
                liquid_volume_m3: volumes[0],
            },
            FiniteHydrostaticLayer {
                source: top,
                liquid_volume_m3: volumes[1],
            },
        ],
        gravity(),
        geometry(height),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn two_liquids_and_bodies_keep_each_inventory_and_match_independent_mass_moments_and_energy() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let tank = stack(&bottom, &top, [16.0; 2], 4.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let z = [17.0 / 16.0 + 0.03, 17.0 / 16.0 - 0.01];
    let states = [0, 1].map(|i| {
        body.state_from_kinematics(at(if i == 0 { -1.0 } else { 1.0 }, z[i]), Twist::zero())
            .unwrap()
    });
    let interfaces = tank
        .interfaces_at_bodies(&bindings, &states, || false)
        .unwrap();
    let load = tank.load_bodies(&bindings, &states, || false).unwrap();
    let hb = (17.0 - z[0] - z[1]) / 14.0;
    let ht = 34.0 / 16.0;
    let d = z.map(|z| hb - z + 0.5);
    close(interfaces[0].surface_height_m, hb);
    close(interfaces[1].surface_height_m, ht);
    close(load.surface.surface_height_m, ht);
    close(load.surface.displaced_volume_m3, 2.0);
    close(tank.liquid_mass_kg(), 28800.0);
    assert_eq!(
        tank.layer_masses_kg().collect::<Vec<_>>(),
        vec![16000.0, 12800.0]
    );
    for (i, interface) in interfaces.iter().enumerate() {
        let displaced = interface
            .displaced_below
            .iter()
            .map(|cut| cut.volume_m3())
            .sum::<f64>();
        close(
            16.0 * interface.surface_height_m - displaced,
            16.0 * (i + 1) as f64,
        );
        assert!(interface.cumulative_volume_defect_m3.abs() <= 1e-12);
    }
    for i in 0..2 {
        close(interfaces[0].displaced_below[i].volume_m3(), d[i]);
        close(interfaces[1].displaced_below[i].volume_m3(), 1.0);
        let mass = 1000.0 * d[i] + 800.0 * (1.0 - d[i]);
        close(load.bodies[i].displaced_mass_kg(), mass);
        close(load.bodies[i].wrench().force_world.z, 10.0 * mass);
        let centroid = (1000.0 * d[i] * (z[i] - 0.5 + hb) / 2.0
            + 800.0 * (1.0 - d[i]) * (hb + z[i] + 0.5) / 2.0)
            / mass;
        close(
            load.bodies[i].center_of_buoyancy_world().unwrap().z,
            centroid,
        );
    }
    // Direct integrate rho*g*z over each liquid slab minus the body portions.
    // Empty reference heights are 1 and 2; no density-jump formula is used.
    let bottom_first_moment =
        8.0 * hb * hb - 0.5 * z.iter().map(|z| hb * hb - (z - 0.5).powi(2)).sum::<f64>();
    let top_first_moment = 8.0 * (ht * ht - hb * hb)
        - 0.5 * z.iter().map(|z| (z + 0.5).powi(2) - hb * hb).sum::<f64>();
    let independent_energy =
        10000.0 * (bottom_first_moment - 8.0) + 8000.0 * (top_first_moment - 24.0);
    close(
        load.surface.liquid_gravitational_energy_change_j,
        independent_energy,
    );

    let moved = [
        states[0],
        body.state_from_kinematics(at(1.0, z[1] + 0.02), Twist::zero())
            .unwrap(),
    ];
    let after = tank.load_bodies(&bindings, &moved, || false).unwrap();
    close(
        after.bodies[0].wrench().force_world.z - load.bodies[0].wrench().force_world.z,
        -2000.0 * 0.02 / 14.0,
    );
}

#[test]
fn stratified_shared_potential_derivatives_match_all_body_forces_and_torques() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let tank = stack(&bottom, &top, [16.0; 2], 4.0);
    let sizes = [Vec3::new(0.8, 1.2, 1.0), Vec3::new(1.0, 0.8, 1.3)];
    let sources = [body(900.0, sizes[0]), body(850.0, sizes[1])];
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
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let tank = stack(&bottom, &top, [16.0; 2], 4.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(900.0, sides);
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
        let interfaces = tank
            .interfaces_at_bodies(&bindings, &states, || false)
            .unwrap();
        close(interfaces[0].surface_height_m, (17.0 - z[0] - z[1]) / 14.0);
        close(interfaces[1].surface_height_m, 34.0 / 16.0);
    }
    assert!(defect.abs() < 0.05);
    states.map(|state| body.world_from_body(state).unwrap().translation().z)
}

#[test]
fn native_collective_motion_converges_to_density_contrast_modes() {
    // Entire bodies lie below the free surface and cross the lower interface.
    // m=900, a=1, A=16: restoring stiffness uses delta_rho=200, not rho=1000.
    // Anti-phase omega²=20/9; in-phase omega²=160/63. One displaced body excites both.
    let in_phase = ((160.0_f64 / 63.0).sqrt() * 1.6).cos();
    let anti_phase = ((20.0_f64 / 9.0).sqrt() * 1.6).cos();
    let exact = [
        17.0 / 16.0 + 0.01 * (in_phase + anti_phase),
        17.0 / 16.0 + 0.01 * (in_phase - anti_phase),
    ];
    let errors: Vec<_> = [(0.08, 20), (0.04, 40), (0.02, 80)]
        .into_iter()
        .map(|(dt, steps)| {
            let actual = coupled_heave(dt, steps);
            assert!((actual[1] - 17.0 / 16.0).abs() > 1e-4);
            ((actual[0] - exact[0]).powi(2) + (actual[1] - exact[1]).powi(2)).sqrt()
        })
        .collect();
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn uniform_stacks_and_single_body_batches_retain_original_native_mechanics_exactly() {
    let water = fluid(1000.0);
    let oil = fluid(800.0);
    let uniform =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(4.0), 32.0, controls()).unwrap();
    let layered_uniform = stack(&water, &water, [16.0; 2], 4.0);
    let layered = stack(&water, &oil, [16.0; 2], 4.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let states = [
        body.state_from_kinematics(at(-1.0, 2.15), Twist::zero())
            .unwrap(),
        body.state_from_kinematics(at(1.0, 2.03), Twist::zero())
            .unwrap(),
    ];
    let controls = HydrostaticStepControls {
        duration_s: 0.001,
        maximum_energy_defect_j: 0.1,
    };
    assert_eq!(
        layered_uniform
            .load_bodies(&bindings, &states, || false)
            .unwrap(),
        uniform.load_bodies(&bindings, &states, || false).unwrap()
    );
    assert_eq!(
        layered_uniform
            .step_bodies(&bindings, &states, controls, || false)
            .unwrap(),
        uniform
            .step_bodies(&bindings, &states, controls, || false)
            .unwrap()
    );
    for z in [1.05, 2.04, 3.0] {
        let held = at(0.0, z);
        let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
        let old = layered
            .step(&mesh, &body, state, controls, || false)
            .unwrap();
        let new = layered
            .step_bodies(&bindings[..1], &[state], controls, || false)
            .unwrap();
        assert_eq!(new.dynamics[0], old.dynamics);
        assert_eq!(new.load_before.bodies[0], old.load_before.buoyancy());
        assert_eq!(new.load_after.bodies[0], old.load_after.buoyancy());
        assert_eq!(new.total_energy_before_j, old.total_energy_before_j);
        assert_eq!(new.total_energy_after_j, old.total_energy_after_j);
        assert_eq!(new.energy_defect_j, old.energy_defect_j);
        let old_interfaces = layered.interfaces_at_pose(&mesh, held, || false).unwrap();
        let new_interfaces = layered
            .interfaces_at_bodies(&bindings[..1], &[state], || false)
            .unwrap();
        for (a, b) in old_interfaces.iter().zip(new_interfaces) {
            assert_eq!(a.surface_height_m, b.surface_height_m);
            assert_eq!(a.displaced_below, b.displaced_below[0]);
            assert_eq!(a.cumulative_volume_defect_m3, b.cumulative_volume_defect_m3);
            assert_eq!(a.volume_evaluations, b.volume_evaluations);
        }
    }
}

#[test]
fn sources_stay_frozen_and_collective_overflow_or_late_refusal_publishes_no_batch() {
    struct Source {
        rho: f64,
        calls: std::cell::Cell<usize>,
    }
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.calls.set(self.calls.get() + 1);
            FluidTransportState::try_new(self.rho, 4000.0, 0.5, 0.001)
        }
    }
    let bottom = Source {
        rho: 1000.0,
        calls: std::cell::Cell::new(0),
    };
    let top = Source {
        rho: 800.0,
        calls: std::cell::Cell::new(0),
    };
    let tank = stack(&bottom, &top, [16.0; 2], 4.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let states = [
        body.state_from_kinematics(at(-1.0, 1.09), Twist::zero())
            .unwrap(),
        body.state_from_kinematics(at(1.0, 1.05), Twist::zero())
            .unwrap(),
    ];
    let originals = states;
    let baseline = tank.load_bodies(&bindings, &states, || false).unwrap();
    assert!(tank.load_bodies(&[], &[], || false).is_err());
    assert!(tank.load_bodies(&bindings, &states[..1], || false).is_err());
    assert!(
        tank.load_bodies(
            &vec![bindings[0]; MAX_TANK_BODIES + 1],
            &vec![states[0]; MAX_TANK_BODIES + 1],
            || false
        )
        .is_err()
    );
    let controls = HydrostaticStepControls {
        duration_s: 0.02,
        maximum_energy_defect_j: 0.1,
    };
    let mut polls = 0;
    tank.step_bodies(&bindings, &states, controls, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [0, polls / 2, polls - 1] {
        let mut count = 0;
        assert!(
            tank.step_bodies(&bindings, &states, controls, || {
                let cancel = count == stop;
                count += 1;
                cancel
            })
            .is_err()
        );
    }
    assert!(
        tank.step_bodies(
            &bindings,
            &states,
            HydrostaticStepControls {
                maximum_energy_defect_j: 1e-30,
                ..controls
            },
            || false
        )
        .is_err()
    );
    assert_eq!((bottom.calls.get(), top.calls.get()), (1, 1));
    assert_eq!(states, originals);
    assert_eq!(
        tank.load_bodies(&bindings, &states, || false).unwrap(),
        baseline
    );

    let overflow = stack(&bottom, &top, [16.0, 14.5], 2.0);
    let submerged = [
        body.state_from_kinematics(at(-1.0, 1.0), Twist::zero())
            .unwrap(),
        body.state_from_kinematics(at(1.0, 1.0), Twist::zero())
            .unwrap(),
    ];
    for state in submerged {
        assert!(overflow.load(&mesh, &body, state, || false).is_ok());
    }
    assert!(matches!(
        overflow.load_bodies(&bindings, &submerged, || false),
        Err(HydrostaticError::TankOverflow { .. })
    ));
    assert!(
        overflow
            .interfaces_at_bodies(&bindings, &submerged, || false)
            .is_err()
    );
    assert!(
        overflow
            .step_bodies(&bindings, &submerged, controls, || false)
            .is_err()
    );
}
