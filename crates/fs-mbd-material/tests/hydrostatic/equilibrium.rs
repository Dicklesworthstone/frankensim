//! G1/G3/G4: source-driven heave balance, finite inventories and atomic refusal.
use super::*;
use fs_mbd_material::hydrostatic::{
    DensityProfileHydrostaticMedium, FiniteHydrostaticLayer, FiniteHydrostaticTank,
    FiniteLayeredHydrostaticTank, FiniteTankControls, HeaveEquilibriumControls,
    HeaveEquilibriumError, HydrostaticDensitySample, HydrostaticLayer, HydrostaticLoad,
    HydrostaticReservoir, HydrostaticTankGeometry, LayeredHydrostaticMedium,
    LinearDensityHydrostaticMedium, MAX_HEAVE_EVALUATIONS, solve_heave_equilibrium,
};
use fs_mbd_material::principal::PrincipalFrame;
use std::cell::Cell;

fn controls(lower: f64, upper: f64) -> HeaveEquilibriumControls {
    HeaveEquilibriumControls {
        lower_offset_m: lower,
        upper_offset_m: upper,
        maximum_force_defect_n: 1e-8,
        maximum_evaluations: 64,
    }
}

fn tank_geometry(height: f64) -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: Vec3::new(4.0, 4.0, height),
    }
}

fn volume_controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-12,
        maximum_volume_evaluations: 64,
    }
}

#[test]
fn material_densities_set_analytic_draft_without_changing_attitude_or_lateral_position() {
    let mesh = cube();
    for (solid_density, liquid_density) in [(500.0, 1000.0), (750.0, 1000.0), (500.0, 1200.0)] {
        let liquid = fluid(liquid_density);
        let medium = HydrostaticMedium::new(&liquid, gravity(), origin()).unwrap();
        let source = body(solid_density, Vec3::new(2.0, 2.0, 2.0));
        let body = source.principal_body(1e-12, || false).unwrap();
        let held = Se3::from_parts(So3::identity(), Vec3::new(0.3, -0.4, 0.37)).unwrap();
        let pose = body
            .state_from_kinematics(held, Twist::zero())
            .unwrap()
            .pose();
        let reference = RigidBodyState::new(
            pose,
            fs_mbd::Vec3::new(1.0, 2.0, 3.0),
            fs_mbd::Vec3::new(4.0, 5.0, 6.0),
        )
        .unwrap();
        let result = solve_heave_equilibrium(
            &medium,
            &mesh,
            body.frame(),
            reference,
            controls(-2.0, 2.0),
            || false,
        )
        .unwrap();
        let expected_height = 1.0 - 2.0 * solid_density / liquid_density;
        vec_close(
            body.world_from_body(result.state).unwrap().translation(),
            Vec3::new(0.3, -0.4, expected_height),
        );
        close(result.offset_m, expected_height - 0.37);
        assert_eq!(
            result.state.pose().orientation(),
            reference.pose().orientation()
        );
        assert_eq!(result.state.linear_momentum_world(), fs_mbd::Vec3::ZERO);
        assert_eq!(result.state.angular_momentum_body(), fs_mbd::Vec3::ZERO);
        assert_eq!(reference.linear_momentum_world().z, 3.0);
        assert_eq!(reference.angular_momentum_body().z, 6.0);
        close(result.load.displaced_mass_kg(), body.properties().mass());
        assert!(result.heave_force_residual_n.abs() <= 1e-8);
        assert_eq!(
            result.load,
            medium.load(&mesh, &body, result.state, || false).unwrap()
        );
    }
}

#[test]
fn immersed_body_finds_the_density_interface_and_the_continuous_profile_balance() {
    let upper = fluid(400.0);
    let lower = fluid(1000.0);
    let layered = LayeredHydrostaticMedium::new(
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
    let source = body(850.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = body
        .state_from_kinematics(pose(0.4), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &layered,
        &cube(),
        body.frame(),
        reference,
        controls(-2.0, 1.0),
        || false,
    )
    .unwrap();
    // The lower fluid occupies 3/4 of the envelope: (850-400)/(1000-400).
    close(
        body.world_from_body(result.state).unwrap().translation().z,
        -0.5,
    );
    close(result.load.displaced().volume_m3(), 8.0);
    close(result.load.displaced_mass_kg(), 6800.0);

    let fluids = [fluid(1000.0), fluid(1400.0), fluid(1600.0)];
    let profile = DensityProfileHydrostaticMedium::new(
        &[
            HydrostaticDensitySample {
                depth_m: 0.0,
                source: &fluids[0],
            },
            HydrostaticDensitySample {
                depth_m: 1.0,
                source: &fluids[1],
            },
            HydrostaticDensitySample {
                depth_m: 3.0,
                source: &fluids[2],
            },
        ],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let source = super::body(1325.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &profile,
        &cube(),
        body.frame(),
        reference,
        controls(-1.6, -0.2),
        || false,
    )
    .unwrap();
    // At centre depth 1 m, integrate BOTH slopes over depths [0,2]:
    // 4*(mean(1000,1400)*1 + mean(1400,1500)*1) = 10600 kg.
    close(
        body.world_from_body(result.state).unwrap().translation().z,
        -1.0,
    );
    close(result.load.displaced_mass_kg(), 10600.0);
    close(profile.density_at_depth_kg_m3(1.0).unwrap(), 1400.0);
}

#[test]
fn continuous_linear_stratification_can_balance_a_fully_submerged_body() {
    let water = fluid(1000.0);
    let medium = LinearDensityHydrostaticMedium::new(&water, 250.0, gravity(), origin()).unwrap();
    let sides = Vec3::new(0.5, 0.5, 0.5);
    let source = body(1250.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let reference = body
        .state_from_kinematics(pose(-0.5), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &medium,
        &mesh,
        body.frame(),
        reference,
        controls(-1.0, 0.0),
        || false,
    )
    .unwrap();
    close(
        body.world_from_body(result.state).unwrap().translation().z,
        -1.0,
    );
    close(result.load.displaced().volume_m3(), 0.125);
    close(result.load.displaced_mass_kg(), 156.25);
    assert!(result.evaluations <= 3);
}

#[test]
fn conserved_finite_bath_waterline_moves_with_the_solved_body_draft() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(4.0),
        16.0,
        volume_controls(),
    )
    .unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let reference = body
        .state_from_kinematics(pose(1.0), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &tank,
        &mesh,
        body.frame(),
        reference,
        controls(-0.4, 0.8),
        || false,
    )
    .unwrap();
    let held = body.world_from_body(result.state).unwrap();
    let fill = tank.fill_at_pose(&mesh, held, || false).unwrap();
    // Mass=2000 kg displaces 2 m3; the conserved bath rises from 1 to 1.125 m.
    close(held.translation().z, 1.125);
    close(fill.surface_height_m, 1.125);
    close(fill.displaced.volume_m3(), 2.0);
    close(
        16.0 * fill.surface_height_m - fill.displaced.volume_m3(),
        16.0,
    );
    assert_eq!(
        result.load,
        tank.load(&mesh, &body, result.state, || false).unwrap()
    );
    close(tank.liquid_mass_kg(), 16000.0);
}

#[test]
fn finite_immiscible_layers_preserve_both_inventories_at_interfacial_equilibrium() {
    let dense = fluid(1000.0);
    let light = fluid(800.0);
    let tank = FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: &dense,
                liquid_volume_m3: 12.0,
            },
            FiniteHydrostaticLayer {
                source: &light,
                liquid_volume_m3: 28.0,
            },
        ],
        gravity(),
        tank_geometry(6.0),
        volume_controls(),
        || false,
    )
    .unwrap();
    let source = body(2500.0 / 3.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = body
        .state_from_kinematics(pose(1.5), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &tank,
        &cube(),
        body.frame(),
        reference,
        controls(-0.3, 0.3),
        || false,
    )
    .unwrap();
    let held = body.world_from_body(result.state).unwrap();
    close(held.translation().z, 1.5);
    let layers = tank.interfaces_at_pose(&cube(), held, || false).unwrap();
    close(layers[0].surface_height_m, 5.0 / 6.0);
    close(layers[1].surface_height_m, 3.0);
    close(layers[0].displaced_below.volume_m3(), 4.0 / 3.0);
    close(layers[1].displaced_below.volume_m3(), 8.0);
    close(
        16.0 * layers[0].surface_height_m - layers[0].displaced_below.volume_m3(),
        12.0,
    );
    close(
        16.0 * (layers[1].surface_height_m - layers[0].surface_height_m)
            - (layers[1].displaced_below.volume_m3() - layers[0].displaced_below.volume_m3()),
        28.0,
    );
    close(result.load.displaced_mass_kg(), body.properties().mass());
}

#[test]
fn heave_balance_reports_the_torque_required_to_hold_a_heterogeneous_body() {
    let left = density_card("light", 500.0, 500.0);
    let right = density_card("heavy", 1500.0, 1500.0);
    let half = VolumeMoments::cuboid(Vec3::new(1.0, 2.0, 2.0)).unwrap();
    let p = point(300.0);
    let source = resolve_body(
        &[
            RegionBinding {
                key: "left",
                material: &left,
                point: &p,
                selection: MaterialPropertySelection::SingleClaimOnly,
                geometry: half
                    .transformed(
                        Se3::from_parts(So3::identity(), Vec3::new(-0.5, 0.0, 0.0)).unwrap(),
                    )
                    .unwrap(),
            },
            RegionBinding {
                key: "right",
                material: &right,
                point: &p,
                selection: MaterialPropertySelection::SingleClaimOnly,
                geometry: half
                    .transformed(
                        Se3::from_parts(So3::identity(), Vec3::new(0.5, 0.0, 0.0)).unwrap(),
                    )
                    .unwrap(),
            },
        ],
        2,
        || false,
    )
    .unwrap();
    let body = source.principal_body(1e-12, || false).unwrap();
    let liquid = fluid(2000.0);
    let medium = HydrostaticMedium::new(&liquid, gravity(), origin()).unwrap();
    let reference = body
        .state_from_kinematics(pose(0.3), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &medium,
        &cube(),
        body.frame(),
        reference,
        controls(-1.0, 1.0),
        || false,
    )
    .unwrap();
    close(
        body.world_from_body(result.state).unwrap().translation().z,
        0.0,
    );
    close(result.state.pose().position_world().x, 0.25);
    let holding = result
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(result.holding_torque_body_nm);
    close(holding.x, 0.0);
    close(holding.y, -20000.0);
    close(holding.z, 0.0);
    assert_eq!(
        result.holding_torque_body_nm,
        result.load.wrench().torque_body.scale(-1.0)
    );
    assert!(result.heave_force_residual_n.abs() <= 1e-8);
}

#[test]
fn world_rotation_and_translation_preserve_heave_solution_and_held_attitude() {
    let liquid = fluid(1000.0);
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    let attitude = So3::exp(So3Tangent::new(Vec3::new(0.2, 0.3, 0.4))).unwrap();
    let initial = Se3::from_parts(attitude, Vec3::new(0.3, -0.4, 0.35)).unwrap();
    let transform = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.4, -0.6, 0.3))).unwrap(),
        Vec3::new(3.0, -2.0, 1.0),
    )
    .unwrap();
    for world in [Se3::identity(), transform] {
        let g = world.rotation().rotate(Vec3::new(0.0, 0.0, -10.0)).unwrap();
        let gravity = Gravity::new(fs_mbd::Vec3::new(g.x, g.y, g.z)).unwrap();
        let medium = HydrostaticMedium::new(&liquid, gravity, world.translation()).unwrap();
        let reference = body
            .state_from_kinematics(world.compose(initial).unwrap(), Twist::zero())
            .unwrap();
        let result = solve_heave_equilibrium(
            &medium,
            &mesh,
            body.frame(),
            reference,
            controls(-1.0, 1.0),
            || false,
        )
        .unwrap();
        // A plane through a centrally symmetric solid divides its volume equally,
        // independently of attitude. Only this translation gives half immersion.
        let expected = world.transform_point(Vec3::new(0.3, -0.4, 0.0)).unwrap();
        vec_close(
            body.world_from_body(result.state).unwrap().translation(),
            expected,
        );
        close(result.offset_m, -0.35);
        assert_eq!(
            result.state.pose().orientation(),
            reference.pose().orientation()
        );
        assert!(result.net_force_world_n.norm_squared().sqrt() <= 2e-8);
    }
}

#[test]
fn endpoint_roots_no_root_brackets_and_evaluation_caps_have_distinct_outcomes() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let mesh = cube();
    let source = body(1000.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let result = solve_heave_equilibrium(
        &medium,
        &mesh,
        body.frame(),
        reference,
        controls(-3.0, -2.0),
        || false,
    )
    .unwrap();
    assert_eq!(result.evaluations, 1);
    close(result.offset_m, -3.0); // Neutral plateau: no unique draft is claimed.

    let source = super::body(1200.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    match solve_heave_equilibrium(
        &medium,
        &mesh,
        body.frame(),
        reference,
        controls(-2.0, 2.0),
        || false,
    ) {
        Err(HeaveEquilibriumError::NotBracketed {
            lower_residual_n,
            upper_residual_n,
            ..
        }) => {
            close(lower_residual_n, -16000.0);
            close(upper_residual_n, -96000.0);
        }
        other => panic!("expected sinking-body refusal, got {other:?}"),
    }
    let source = super::body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let cap = HeaveEquilibriumControls {
        maximum_evaluations: 2,
        ..controls(-2.0, 2.0)
    };
    assert!(matches!(
        solve_heave_equilibrium(&medium, &mesh, body.frame(), reference, cap, || false),
        Err(HeaveEquilibriumError::NotConverged { evaluations: 2, .. })
    ));
}

#[test]
fn source_is_frozen_and_every_cancellation_boundary_refuses_without_publishing_state() {
    struct CountedFluid(Cell<usize>);
    impl FluidTransport for CountedFluid {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get() + 1);
            Ok(fluid(1000.0))
        }
    }
    struct Counted<'a, R>(&'a R, Cell<usize>);
    impl<R: HydrostaticReservoir> HydrostaticReservoir for Counted<'_, R> {
        fn gravity(&self) -> Gravity {
            self.0.gravity()
        }
        fn load_frame(
            &self,
            mesh: &DisplacementMesh,
            body: &PrincipalFrame,
            state: RigidBodyState,
            cancelled: impl FnMut() -> bool,
        ) -> Result<HydrostaticLoad, HydrostaticError> {
            self.1.set(self.1.get() + 1);
            self.0.load_frame(mesh, body, state, cancelled)
        }
    }
    let liquid = CountedFluid(Cell::new(0));
    let medium = HydrostaticMedium::new(&liquid, gravity(), origin()).unwrap();
    let counted = Counted(&medium, Cell::new(0));
    let source = body(730.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = body
        .state_from_kinematics(pose(0.25), Twist::zero())
        .unwrap();
    let mesh = cube();
    let controls = controls(-2.0, 1.0);
    let mut polls = 0;
    let expected =
        solve_heave_equilibrium(&counted, &mesh, body.frame(), reference, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(counted.1.get(), expected.evaluations);
    assert!(expected.evaluations <= controls.maximum_evaluations);
    for stop in 1..=polls {
        let mut seen = 0;
        assert_eq!(
            solve_heave_equilibrium(&counted, &mesh, body.frame(), reference, controls, || {
                seen += 1;
                seen == stop
            }),
            Err(HeaveEquilibriumError::Physics(HydrostaticError::Cancelled)),
        );
    }
    assert_eq!(
        solve_heave_equilibrium(&medium, &mesh, body.frame(), reference, controls, || false)
            .unwrap(),
        expected
    );
    assert_eq!(liquid.0.get(), 1);

    let invalid = [
        HeaveEquilibriumControls {
            lower_offset_m: f64::NAN,
            ..controls
        },
        HeaveEquilibriumControls {
            upper_offset_m: f64::INFINITY,
            ..controls
        },
        HeaveEquilibriumControls {
            upper_offset_m: controls.lower_offset_m,
            ..controls
        },
        HeaveEquilibriumControls {
            lower_offset_m: -f64::MAX,
            upper_offset_m: f64::MAX,
            ..controls
        },
        HeaveEquilibriumControls {
            maximum_force_defect_n: 0.0,
            ..controls
        },
        HeaveEquilibriumControls {
            maximum_force_defect_n: f64::NAN,
            ..controls
        },
        HeaveEquilibriumControls {
            maximum_evaluations: 1,
            ..controls
        },
        HeaveEquilibriumControls {
            maximum_evaluations: MAX_HEAVE_EVALUATIONS + 1,
            ..controls
        },
    ];
    counted.1.set(0);
    for control in invalid {
        assert!(matches!(
            solve_heave_equilibrium(&counted, &mesh, body.frame(), reference, control, || false),
            Err(HeaveEquilibriumError::Physics(HydrostaticError::Invalid(_)))
        ));
    }
    assert_eq!(counted.1.get(), 0);
}

#[test]
fn underlying_tank_overflow_contact_and_inner_work_refusals_are_preserved() {
    let water = fluid(1000.0);
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = body
        .state_from_kinematics(pose(1.5), Twist::zero())
        .unwrap();
    let mesh = cube();
    let narrow = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(1.1),
        16.0,
        volume_controls(),
    )
    .unwrap();
    assert!(matches!(
        solve_heave_equilibrium(
            &narrow,
            &mesh,
            body.frame(),
            reference,
            controls(0.0, 1.0),
            || false
        ),
        Err(HeaveEquilibriumError::Physics(
            HydrostaticError::TankOverflow { .. }
        ))
    ));
    let budget = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        tank_geometry(4.0),
        16.0,
        FiniteTankControls {
            maximum_volume_evaluations: 2,
            ..volume_controls()
        },
    )
    .unwrap();
    assert!(matches!(
        solve_heave_equilibrium(
            &budget,
            &mesh,
            body.frame(),
            reference,
            controls(0.0, 1.0),
            || false
        ),
        Err(HeaveEquilibriumError::Physics(
            HydrostaticError::TankVolumeNotConverged { evaluations: 2, .. }
        ))
    ));
    assert!(matches!(
        solve_heave_equilibrium(
            &budget,
            &mesh,
            body.frame(),
            reference,
            controls(-1.0, 1.0),
            || false
        ),
        Err(HeaveEquilibriumError::Physics(HydrostaticError::Invalid(_)))
    ));
    close(
        body.world_from_body(reference).unwrap().translation().z,
        1.5,
    );
    close(narrow.liquid_mass_kg(), 16000.0);
}
