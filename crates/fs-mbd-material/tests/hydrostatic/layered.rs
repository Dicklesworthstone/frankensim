//! G1/G3/G4: independently integrated layers, native motion and atomic refusal.
use super::*;
use fs_mbd_material::hydrostatic::{
    HydrostaticLayer, HydrostaticStepControls, HydrostaticStepError, LayeredHydrostaticMedium,
};

#[test]
fn gas_over_two_liquids_uses_density_weighted_moments_and_continuous_pressure() {
    let gas = fluid(2.0);
    let light = fluid(800.0);
    let dense = fluid(1000.0);
    let medium = LayeredHydrostaticMedium::new(
        Some(&gas),
        &[
            HydrostaticLayer {
                depth_m: -0.25,
                source: &light,
            },
            HydrostaticLayer {
                depth_m: 0.5,
                source: &dense,
            },
        ],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let load = medium.load(&cube(), &body, state, || false).unwrap();
    // Independent layer slabs: V=(3,3,2), z_centroid=(.625,-.125,-.75).
    let mass = 2.0 * 3.0 + 800.0 * 3.0 + 1000.0 * 2.0;
    close(load.displaced().volume_m3(), 8.0);
    vec_close(load.displaced().centroid_body_m().unwrap(), origin());
    close(load.displaced_mass_kg(), mass);
    close(load.wrench().force_world.z, 10.0 * mass);
    close(
        load.center_of_buoyancy_world().unwrap().z,
        (2.0 * 3.0 * 0.625 - 800.0 * 3.0 * 0.125 - 1000.0 * 2.0 * 0.75) / mass,
    );
    close(load.pressure_potential_j(), 25937.5);
    assert_eq!(medium.layers().len(), 2);
    assert_eq!(medium.above().unwrap().1, gas);

    // Move only the pressure datum while keeping the physical interfaces fixed;
    // inserting a zero-density jump must also leave every force moment intact.
    let shifted = LayeredHydrostaticMedium::new(
        Some(&gas),
        &[
            HydrostaticLayer {
                depth_m: 9.75,
                source: &light,
            },
            HydrostaticLayer {
                depth_m: 10.0,
                source: &light,
            },
            HydrostaticLayer {
                depth_m: 10.5,
                source: &dense,
            },
        ],
        gravity(),
        Vec3::new(0.0, 0.0, 10.0),
        || false,
    )
    .unwrap();
    let shifted_load = shifted.load(&cube(), &body, state, || false).unwrap();
    close(
        shifted_load.displaced().volume_m3(),
        load.displaced().volume_m3(),
    );
    close(shifted_load.displaced_mass_kg(), mass);
    close(
        shifted_load.wrench().force_world.z,
        load.wrench().force_world.z,
    );
    vec_close(
        shifted_load.center_of_buoyancy_world().unwrap(),
        load.center_of_buoyancy_world().unwrap(),
    );
    // Only the uniform ambient term shifts, by rho_top*g*datum_shift*V.
    close(
        shifted_load.pressure_potential_j() - load.pressure_potential_j(),
        1600.0,
    );
}

#[test]
fn one_layer_recovers_existing_free_surface_load_and_motion() {
    let water = fluid(1000.0);
    let single = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let layered = LayeredHydrostaticMedium::new(
        None,
        &[HydrostaticLayer {
            depth_m: 0.0,
            source: &water,
        }],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let mut left = body
        .state_from_kinematics(pose(0.02), Twist::zero())
        .unwrap();
    let mut right = left;
    let controls = HydrostaticStepControls {
        duration_s: 0.001,
        maximum_energy_defect_j: 0.1,
    };
    for _ in 0..100 {
        let a = single.step(&mesh, &body, left, controls, || false).unwrap();
        let b = layered
            .step(&mesh, &body, right, controls, || false)
            .unwrap();
        close(
            a.load_after.displaced_mass_kg(),
            b.load_after.displaced_mass_kg(),
        );
        close(
            a.load_after.pressure_potential_j(),
            b.load_after.pressure_potential_j(),
        );
        left = a.dynamics.state_after;
        right = b.dynamics.state_after;
    }
    close(
        left.pose().position_world().z,
        right.pose().position_world().z,
    );
    close(
        left.linear_momentum_world().z,
        right.linear_momentum_world().z,
    );
}

#[test]
fn gas_only_lift_uses_actual_source_pressure_and_temperature_and_native_inertia() {
    use fs_material::gas::{GasSpec, GasState};
    let first = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    let second = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200000.0).unwrap();
    let heated = GasState::try_new(&GasSpec::dry_air_ussa1976(), 600.0, 100000.0).unwrap();
    let rho = first.transport_state().unwrap().density_kg_m3();
    let source = body(rho / 2.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(3.0), Twist::zero())
        .unwrap();
    for (gas, acceleration) in [(&first, 10.0), (&second, 30.0), (&heated, 0.0)] {
        let medium =
            LayeredHydrostaticMedium::new(Some(gas), &[], gravity(), origin(), || false).unwrap();
        let step = medium
            .step(
                &cube(),
                &body,
                state,
                HydrostaticStepControls {
                    duration_s: 1e-4,
                    maximum_energy_defect_j: 0.01,
                },
                || false,
            )
            .unwrap();
        close(
            step.dynamics
                .state_after
                .center_of_mass_velocity_world(body.properties())
                .unwrap()
                .z
                / 1e-4,
            acceleration,
        );
        // A body above the arbitrary reference plane still displaces gas.
        close(step.load_after.displaced().volume_m3(), 8.0);
        assert!(step.load_after.pressure_potential_j() < 0.0);
    }
}

#[test]
fn layered_pressure_potential_differentiates_to_force_and_principal_torque() {
    let upper = fluid(400.0);
    let lower = fluid(1000.0);
    let medium = LayeredHydrostaticMedium::new(
        Some(&upper),
        &[HydrostaticLayer {
            depth_m: 0.0,
            source: &lower,
        }],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let sides = Vec3::new(1.2, 1.8, 0.8);
    let source = body(700.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.15, 0.3, 0.1))).unwrap();
    let state = body
        .state_from_kinematics(
            Se3::from_parts(rotation, Vec3::new(0.1, -0.2, 0.05)).unwrap(),
            Twist::zero(),
        )
        .unwrap();
    let load = medium.load(&mesh, &body, state, || false).unwrap();
    let h = 1e-6;
    for unit in [
        fs_mbd::Vec3::new(1.0, 0.0, 0.0),
        fs_mbd::Vec3::new(0.0, 1.0, 0.0),
        fs_mbd::Vec3::new(0.0, 0.0, 1.0),
    ] {
        for rotational in [false, true] {
            let energy = |sign: f64| {
                let old = state.pose();
                let pos = if rotational {
                    old.position_world()
                } else {
                    old.position_world().add(unit.scale(sign * h))
                };
                let orient = if rotational {
                    old.orientation().right_exp(unit.scale(sign * h)).unwrap()
                } else {
                    old.orientation()
                };
                let changed = RigidBodyState::new(
                    Pose::new(pos, orient).unwrap(),
                    state.linear_momentum_world(),
                    state.angular_momentum_body(),
                )
                .unwrap();
                medium
                    .load(&mesh, &body, changed, || false)
                    .unwrap()
                    .pressure_potential_j()
            };
            let force = if rotational {
                load.wrench().torque_body
            } else {
                load.wrench().force_world
            };
            let derivative = -(energy(1.0) - energy(-1.0)) / (2.0 * h);
            assert!(
                (derivative - force.dot(unit)).abs() < 1e-5,
                "rotation={rotational}: {derivative} != {force:?}"
            );
        }
    }
}

#[test]
fn density_contrast_controls_interfacial_equilibrium_and_heave_stiffness() {
    let light = fluid(800.0);
    let dense = fluid(1000.0);
    let medium = LayeredHydrostaticMedium::new(
        Some(&light),
        &[HydrostaticLayer {
            depth_m: 0.0,
            source: &dense,
        }],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let controls = HydrostaticStepControls {
        duration_s: 0.001,
        maximum_energy_defect_j: 0.01,
    };
    let rest = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let result = medium.step(&mesh, &body, rest, controls, || false).unwrap();
    close(result.dynamics.state_after.linear_momentum_world().z, 0.0);
    let mut state = body
        .state_from_kinematics(pose(0.02), Twist::zero())
        .unwrap();
    for _ in 0..1000 {
        state = medium
            .step(&mesh, &body, state, controls, || false)
            .unwrap()
            .dynamics
            .state_after;
    }
    // m=3600 kg, restoring stiffness=(1000-800)*g*A=8000 N/m.
    let exact = 0.02 * (8000.0_f64 / 3600.0).sqrt().cos();
    assert!((body.world_from_body(state).unwrap().translation().z - exact).abs() < 2e-5);
}

#[test]
fn inadmissible_layers_refuse_and_cancelled_steps_publish_no_candidate() {
    struct Counted {
        count: std::cell::Cell<usize>,
        state: FluidTransportState,
    }
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.count.set(self.count.get() + 1);
            Ok(self.state)
        }
    }
    let upper = Counted {
        count: std::cell::Cell::new(0),
        state: fluid(800.0),
    };
    let lower = Counted {
        count: std::cell::Cell::new(0),
        state: fluid(1000.0),
    };
    let medium = LayeredHydrostaticMedium::new(
        Some(&upper),
        &[HydrostaticLayer {
            depth_m: 0.0,
            source: &lower,
        }],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let source = body(900.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(0.02), Twist::zero())
        .unwrap();
    let mesh = cube();
    let controls = HydrostaticStepControls {
        duration_s: 0.001,
        maximum_energy_defect_j: 0.1,
    };
    let mut count = 0;
    let accepted = medium
        .step(&mesh, &body, state, controls, || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        medium.step(&mesh, &body, state, controls, || {
            calls += 1;
            calls == count
        }),
        Err(HydrostaticStepError::Physics(HydrostaticError::Cancelled))
    ));
    assert_eq!(
        medium
            .step(&mesh, &body, state, controls, || false)
            .unwrap(),
        accepted
    );
    assert_eq!(upper.count.get(), 1);
    assert_eq!(lower.count.get(), 1);
    assert!(
        LayeredHydrostaticMedium::new(
            Some(&lower),
            &[HydrostaticLayer {
                depth_m: 0.0,
                source: &upper
            }],
            gravity(),
            origin(),
            || false
        )
        .is_err()
    );
    for depths in [[0.0, 0.0], [0.0, -1.0], [0.0, f64::NAN]] {
        assert!(
            LayeredHydrostaticMedium::new(
                None,
                &depths.map(|depth_m| HydrostaticLayer {
                    depth_m,
                    source: &lower
                }),
                gravity(),
                origin(),
                || false
            )
            .is_err()
        );
    }
    assert!(LayeredHydrostaticMedium::new(None, &[], gravity(), origin(), || false).is_err());
    struct Refusing;
    impl FluidTransport for Refusing {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            FluidTransportState::try_new(0.0, 1.0, 1.0, 1.0)
        }
    }
    assert!(matches!(
        LayeredHydrostaticMedium::new(Some(&Refusing), &[], gravity(), origin(), || false),
        Err(HydrostaticError::Material(_))
    ));
}
