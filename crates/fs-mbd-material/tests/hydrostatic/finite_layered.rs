//! Conserved immiscible inventories reaching actual native rigid-body motion.
#[path = "finite_layered/density_change.rs"]
mod density_change;
use super::*;
use fs_mbd_material::hydrostatic::{
    FiniteHydrostaticLayer, FiniteHydrostaticTank, FiniteLayeredHydrostaticTank,
    FiniteTankControls, HydrostaticStepControls, HydrostaticTankGeometry,
};

fn controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-12,
        maximum_volume_evaluations: 64,
    }
}

fn geometry(height: f64) -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: Vec3::new(4.0, 4.0, height),
    }
}

fn tank<'a>(
    bottom: &'a dyn FluidTransport,
    top: &'a dyn FluidTransport,
    top_volume: f64,
) -> FiniteLayeredHydrostaticTank<'a> {
    FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: bottom,
                liquid_volume_m3: 12.0,
            },
            FiniteHydrostaticLayer {
                source: top,
                liquid_volume_m3: top_volume,
            },
        ],
        gravity(),
        geometry(6.0),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn two_liquids_conserve_individual_mass_and_match_independent_fluid_energy() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let tank = tank(&bottom, &top, 8.0);
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    close(tank.liquid_mass_kg(), 18400.0);
    close(tank.liquid_volume_m3(), 20.0);
    let inventories: Vec<_> = tank
        .layers()
        .map(|(input, state)| input.liquid_volume_m3 * state.density_kg_m3())
        .collect();
    assert_eq!(inventories, vec![12000.0, 6400.0]);

    for z in [1.2_f64, 1.5, 4.0] {
        let held = pose(z);
        let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
        let interfaces = tank.interfaces_at_pose(&mesh, held, || false).unwrap();
        let load = tank.load(&mesh, &body, state, || false).unwrap();
        let lower = z - 1.0;
        let upper = z + 1.0;
        let height = |volume: f64| {
            if volume / 16.0 <= lower {
                volume / 16.0
            } else if (volume + 8.0) / 16.0 >= upper {
                (volume + 8.0) / 16.0
            } else {
                (volume - 4.0 * lower) / 12.0
            }
        };
        let h = [height(12.0), height(20.0)];
        let cut_volume = |height: f64| 4.0 * (height.clamp(lower, upper) - lower);
        let cut_first = |height: f64| 2.0 * (height.clamp(lower, upper).powi(2) - lower * lower);
        let d = [cut_volume(h[0]), cut_volume(h[1])];
        let first = [cut_first(h[0]), cut_first(h[1])];
        for i in 0..2 {
            close(interfaces[i].surface_height_m, h[i]);
            close(interfaces[i].displaced_below.volume_m3(), d[i]);
            assert!(interfaces[i].cumulative_volume_defect_m3.abs() <= 1e-12);
        }
        close(16.0 * h[0] - d[0], 12.0);
        close(16.0 * (h[1] - h[0]) - (d[1] - d[0]), 8.0);
        let mass = 1000.0 * d[0] + 800.0 * (d[1] - d[0]);
        close(load.displaced_mass_kg(), mass);
        close(load.displaced().volume_m3(), d[1]);
        close(load.wrench().force_world.z, 10.0 * mass);
        if mass > 0.0 {
            close(
                load.center_of_buoyancy_world().unwrap().z,
                (1000.0 * first[0] + 800.0 * (first[1] - first[0])) / mass,
            );
        } else {
            assert_eq!(load.center_of_buoyancy_world(), None);
        }
        // Independent first moments of the actual two liquid slabs, removing
        // each submerged body slice, relative to the same empty-tank liquids.
        let actual = 10.0
            * (1000.0 * (8.0 * h[0].powi(2) - first[0])
                + 800.0 * (8.0 * (h[1].powi(2) - h[0].powi(2)) - (first[1] - first[0])));
        let empty = 10.0
            * (1000.0 * 8.0 * 0.75_f64.powi(2)
                + 800.0 * 8.0 * (1.25_f64.powi(2) - 0.75_f64.powi(2)));
        close(load.pressure_potential_j(), actual - empty);
    }
}

#[test]
fn moving_interface_potential_gradient_is_the_native_force_and_torque() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let tank = tank(&bottom, &top, 8.0);
    let sides = Vec3::new(1.2, 1.8, 0.8);
    let source = body(800.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let held = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.15, 0.3, 0.1))).unwrap(),
        Vec3::new(0.1, -0.2, 1.0),
    )
    .unwrap();
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    let load = tank.load(&mesh, &body, state, || false).unwrap();
    let h = 1e-5;
    for unit in [
        fs_mbd::Vec3::new(1.0, 0.0, 0.0),
        fs_mbd::Vec3::new(0.0, 1.0, 0.0),
        fs_mbd::Vec3::new(0.0, 0.0, 1.0),
    ] {
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
                "rotation={rotational}, derivative={derivative}, wrench={wrench:?}"
            );
        }
    }
}

fn interfacial_heave(dt: f64, steps: usize) -> f64 {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let tank = tank(&bottom, &top, 28.0);
    let source = body(2500.0 / 3.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    let mut state = body
        .state_from_kinematics(pose(1.52), Twist::zero())
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
                    maximum_energy_defect_j: 0.01,
                },
                || false,
            )
            .unwrap();
        defect += next.energy_defect_j;
        state = next.dynamics.state_after;
        let held = body.world_from_body(state).unwrap();
        let z = held.translation().z;
        let interfaces = tank.interfaces_at_pose(&mesh, held, || false).unwrap();
        close(interfaces[0].surface_height_m, 4.0 / 3.0 - z / 3.0);
        close(interfaces[1].surface_height_m, 3.0);
        close(
            next.load_after.displaced_mass_kg(),
            (24800.0 - 3200.0 * z) / 3.0,
        );
    }
    assert!(defect.abs() < 0.01);
    body.world_from_body(state).unwrap().translation().z
}

#[test]
fn actual_native_motion_has_density_contrast_and_finite_area_restoring_stiffness() {
    // Body remains fully below the free surface, crossing the lower interface.
    // M=20000/3, k=(1000-800)*10*4*16/(16-4)=32000/3, so omega²=1.6.
    // Holding the internal interface fixed would instead give omega²=1.2.
    let exact = 1.5 + 0.02 * (1.6_f64.sqrt() * 1.6).cos();
    let errors: Vec<_> = [(0.04, 40), (0.02, 80), (0.01, 160)]
        .into_iter()
        .map(|(dt, steps)| (interfacial_heave(dt, steps) - exact).abs())
        .collect();
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn uniform_layers_are_the_exact_original_total_inventory_load_and_step() {
    let water = fluid(1000.0);
    let uniform =
        FiniteHydrostaticTank::new(&water, gravity(), geometry(6.0), 20.0, controls()).unwrap();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    for volumes in [&[20.0][..], &[12.0, 8.0][..]] {
        let inputs: Vec<_> = volumes
            .iter()
            .map(|&volume| FiniteHydrostaticLayer {
                source: &water,
                liquid_volume_m3: volume,
            })
            .collect();
        let layered = FiniteLayeredHydrostaticTank::new(
            &inputs,
            gravity(),
            geometry(6.0),
            controls(),
            || false,
        )
        .unwrap();
        for z in [1.2, 1.5, 4.0] {
            let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
            assert_eq!(
                layered.load(&mesh, &body, state, || false).unwrap(),
                uniform.load(&mesh, &body, state, || false).unwrap()
            );
            let step_controls = HydrostaticStepControls {
                duration_s: 0.001,
                maximum_energy_defect_j: 0.1,
            };
            assert_eq!(
                layered
                    .step(&mesh, &body, state, step_controls, || false)
                    .unwrap(),
                uniform
                    .step(&mesh, &body, state, step_controls, || false)
                    .unwrap()
            );
        }
    }
}

#[test]
fn sources_are_frozen_and_invalid_layers_overflow_and_late_cancellation_publish_nothing() {
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
    let tank = tank(&bottom, &top, 8.0);
    assert_eq!((bottom.calls.get(), top.calls.get()), (1, 1));
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    let held = pose(1.5);
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    let interfaces = tank.interfaces_at_pose(&mesh, held, || false).unwrap();
    let mut polls = 0;
    let expected = tank
        .load(&mesh, &body, state, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            tank.load(&mesh, &body, state, || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticError::Cancelled)
        ));
    }
    assert_eq!(tank.load(&mesh, &body, state, || false).unwrap(), expected);
    assert_eq!(
        tank.interfaces_at_pose(&mesh, held, || false).unwrap(),
        interfaces
    );
    assert_eq!((bottom.calls.get(), top.calls.get()), (1, 1));

    let inputs = [
        FiniteHydrostaticLayer {
            source: &bottom,
            liquid_volume_m3: 12.0,
        },
        FiniteHydrostaticLayer {
            source: &top,
            liquid_volume_m3: 8.0,
        },
    ];
    let shallow =
        FiniteLayeredHydrostaticTank::new(&inputs, gravity(), geometry(1.4), controls(), || false)
            .unwrap();
    assert!(matches!(
        shallow.load(&mesh, &body, state, || false),
        Err(HydrostaticError::TankOverflow { .. })
    ));
    let mut unstable = inputs;
    unstable.swap(0, 1);
    assert!(
        FiniteLayeredHydrostaticTank::new(&unstable, gravity(), geometry(6.0), controls(), || {
            false
        },)
        .is_err()
    );
    for volume in [0.0, f64::NAN, f64::MIN_POSITIVE] {
        let mut invalid = inputs;
        invalid[1].liquid_volume_m3 = volume;
        assert!(
            FiniteLayeredHydrostaticTank::new(
                &invalid,
                gravity(),
                geometry(6.0),
                controls(),
                || false,
            )
            .is_err()
        );
    }
    let refused = Source {
        rho: -1.0,
        calls: std::cell::Cell::new(0),
    };
    assert!(matches!(
        FiniteLayeredHydrostaticTank::new(
            &[FiniteHydrostaticLayer {
                source: &refused,
                liquid_volume_m3: 1.0
            }],
            gravity(),
            geometry(6.0),
            controls(),
            || false,
        ),
        Err(HydrostaticError::Material(_))
    ));
}
