//! G1/G3/G4: material-derived floating poses, actual torque balance and bounded refusal.
use super::*;
use fs_mbd_material::hydrostatic::{
    DensityProfileHydrostaticMedium, FiniteHydrostaticLayer, FiniteHydrostaticTank,
    FiniteLayeredHydrostaticTank, FiniteTankControls, FloatingEquilibrium,
    FloatingEquilibriumControls, FloatingEquilibriumError, HydrostaticDensitySample,
    HydrostaticLayer, HydrostaticLoad, HydrostaticReservoir, HydrostaticTankGeometry,
    LayeredHydrostaticMedium, LinearDensityHydrostaticMedium, MAX_FLOATING_EVALUATIONS,
    MAX_FLOATING_NEWTON_ITERATIONS, solve_floating_equilibrium,
};
use fs_mbd_material::principal::PrincipalFrame;
use std::cell::Cell;

fn controls() -> FloatingEquilibriumControls {
    FloatingEquilibriumControls {
        roll_direction_world: fs_mbd::Vec3::new(1.0, 0.0, 0.0),
        maximum_force_defect_n: 1e-6,
        maximum_torque_defect_nm: 1e-6,
        length_scale_m: 1.0,
        difference_step: 1e-6,
        maximum_heave_offset_m: 1.0,
        maximum_tilt_radians: 1.2,
        maximum_newton_iterations: 32,
        maximum_evaluations: 2048,
    }
}

fn norm(v: fs_mbd::Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn initial(frame: &PrincipalFrame, z: f64, tilt: Vec3) -> RigidBodyState {
    frame
        .state_from_kinematics(
            Se3::from_parts(
                So3::exp(So3Tangent::new(tilt)).unwrap(),
                Vec3::new(0.0, 0.0, z),
            )
            .unwrap(),
            Twist::zero(),
        )
        .unwrap()
}

fn balanced<R: HydrostaticReservoir>(
    result: FloatingEquilibrium,
    reservoir: &R,
    mesh: &DisplacementMesh,
    frame: &PrincipalFrame,
) {
    assert!(norm(result.net_force_world_n) <= controls().maximum_force_defect_n);
    assert!(norm(result.torque_world_nm) <= controls().maximum_torque_defect_nm);
    assert_eq!(result.state.linear_momentum_world(), fs_mbd::Vec3::ZERO);
    assert_eq!(result.state.angular_momentum_body(), fs_mbd::Vec3::ZERO);
    assert_eq!(
        result.load,
        reservoir
            .load_frame(mesh, frame, result.state, || false)
            .unwrap()
    );
    assert!(result.evaluations <= controls().maximum_evaluations);
    assert!(result.iterations <= controls().maximum_newton_iterations);
}

#[test]
fn heel_trim_and_draft_balance_the_original_body_and_fluid_densities() {
    let sides = Vec3::new(3.0, 2.0, 1.0);
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    for (body_density, fluid_density) in [(500.0, 1000.0), (700.0, 1000.0), (500.0, 1250.0)] {
        let water = fluid(fluid_density);
        let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
        let source = body(body_density, sides);
        let body = source.principal_body(1e-12, || false).unwrap();
        let pose = initial(body.frame(), 0.18, Vec3::new(0.13, -0.09, 0.0)).pose();
        let reference = RigidBodyState::new(
            pose,
            fs_mbd::Vec3::new(1.0, 2.0, 3.0),
            fs_mbd::Vec3::new(4.0, 5.0, 6.0),
        )
        .unwrap();
        let result =
            solve_floating_equilibrium(&medium, &mesh, body.frame(), reference, controls(), || {
                false
            })
            .unwrap();
        balanced(result, &medium, &mesh, body.frame());
        let resolved = body.world_from_body(result.state).unwrap();
        near(
            resolved.translation().z,
            0.5 - body_density / fluid_density,
            2e-9,
        );
        let vertical = resolved
            .rotation()
            .rotate(Vec3::new(0.0, 0.0, 1.0))
            .unwrap();
        near(vertical.x, 0.0, 2e-9);
        near(vertical.y, 0.0, 2e-9);
        near(vertical.z, 1.0, 2e-9);
        assert_eq!(
            result.state.pose().position_world().x,
            reference.pose().position_world().x
        );
        assert_eq!(
            result.state.pose().position_world().y,
            reference.pose().position_world().y
        );
        assert!(norm(result.tilt_rotation_vector_world_radians) > 0.1);
        assert!(result.iterations > 0);
        assert_eq!(reference.linear_momentum_world().z, 3.0);
        assert_eq!(reference.angular_momentum_body().z, 6.0);
    }
}

#[test]
fn heterogeneous_mass_resolves_the_independently_integrated_nonzero_pitch() {
    let left = density_card("light-half", 400.0, 400.0);
    let right = density_card("heavy-half", 600.0, 600.0);
    let half = VolumeMoments::cuboid(Vec3::new(1.0, 2.0, 1.0)).unwrap();
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
    let mesh = DisplacementMesh::cuboid(Vec3::new(2.0, 2.0, 1.0), || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let reference = initial(body.frame(), 0.2, origin());
    let result =
        solve_floating_equilibrium(&medium, &mesh, body.frame(), reference, controls(), || {
            false
        })
        .unwrap();
    balanced(result, &medium, &mesh, body.frame());

    // Independent Cartesian integration for a half-wet 2 x 2 x 1 box:
    // Bx=(2/3)t, Bz=-1/4+t²/3, Gx=1/10, t=tan(pitch).
    // (R(B-G))_x=0 gives 5*t+4*t³=1.2. Cardano's real root does
    // not call the clipped-moment or equilibrium implementation.
    let discriminant = 0.15_f64.powi(2) + (5.0_f64 / 12.0).powi(3);
    let t = (0.15 + discriminant.sqrt()).cbrt() + (0.15 - discriminant.sqrt()).cbrt();
    let angle = t.atan();
    assert!(t.abs() < 0.5); // The free surface crosses only the vertical sides.
    near(5.0 * t + 4.0 * t.powi(3), 1.2, 1e-14);
    let resolved = body.world_from_body(result.state).unwrap();
    let x = resolved
        .rotation()
        .rotate(Vec3::new(1.0, 0.0, 0.0))
        .unwrap();
    near((-x.z).atan2(x.x), angle, 2e-9);
    near(x.y, 0.0, 2e-9);
    near(resolved.translation().z, 0.0, 2e-9);
    near(result.state.pose().position_world().x, 0.1, 1e-12);
    near(
        result.state.pose().position_world().z,
        -0.1 * angle.sin(),
        2e-9,
    );
    near(result.load.displaced_mass_kg(), 2000.0, 1e-7);
    assert!(result.load.wrench().torque_body.norm_squared() < 1e-12);
}

#[test]
fn layered_and_continuous_profiles_select_draft_and_attitude_together() {
    let sides = Vec3::new(3.0, 2.0, 1.0);
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
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
    let source = body(700.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = initial(body.frame(), 0.13, Vec3::new(0.08, -0.06, 0.0));
    let result =
        solve_floating_equilibrium(&layered, &mesh, body.frame(), reference, controls(), || {
            false
        })
        .unwrap();
    balanced(result, &layered, &mesh, body.frame());
    let resolved = body.world_from_body(result.state).unwrap();
    near(resolved.translation().z, 0.0, 2e-9);
    near(
        resolved
            .rotation()
            .rotate(Vec3::new(0.0, 0.0, 1.0))
            .unwrap()
            .z,
        1.0,
        2e-9,
    );
    near(result.load.displaced_mass_kg(), 4200.0, 1e-7);

    let deep = fluid(1800.0);
    let profile = DensityProfileHydrostaticMedium::new(
        &[
            HydrostaticDensitySample {
                depth_m: 0.0,
                source: &lower,
            },
            HydrostaticDensitySample {
                depth_m: 4.0,
                source: &deep,
            },
        ],
        gravity(),
        origin(),
        || false,
    )
    .unwrap();
    let linear = LinearDensityHydrostaticMedium::new(&lower, 200.0, gravity(), origin()).unwrap();
    let source = super::body(1400.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = initial(body.frame(), -1.8, Vec3::new(0.08, -0.06, 0.0));
    let controls = FloatingEquilibriumControls {
        maximum_heave_offset_m: 0.4,
        maximum_tilt_radians: 0.4,
        ..controls()
    };
    let result =
        solve_floating_equilibrium(&profile, &mesh, body.frame(), reference, controls, || false)
            .unwrap();
    let linear_result =
        solve_floating_equilibrium(&linear, &mesh, body.frame(), reference, controls, || false)
            .unwrap();
    balanced(result, &profile, &mesh, body.frame());
    balanced(linear_result, &linear, &mesh, body.frame());
    // For a fully submerged centrally symmetric box, mean rho is rho at its
    // geometric centre. 1000 + 200*depth = 1400 puts that centre 2 m down.
    for candidate in [result, linear_result] {
        let resolved = body.world_from_body(candidate.state).unwrap();
        near(resolved.translation().z, -2.0, 2e-9);
        near(
            resolved
                .rotation()
                .rotate(Vec3::new(0.0, 0.0, 1.0))
                .unwrap()
                .z,
            1.0,
            2e-9,
        );
        near(candidate.load.displaced().volume_m3(), 6.0, 1e-12);
    }
}

fn tank_geometry() -> HydrostaticTankGeometry {
    HydrostaticTankGeometry {
        world_from_tank: Se3::identity(),
        inside_dimensions_m: Vec3::new(4.0, 4.0, 4.0),
    }
}

fn volume_controls() -> FiniteTankControls {
    FiniteTankControls {
        maximum_volume_defect_m3: 1e-12,
        maximum_volume_evaluations: 64,
    }
}

#[test]
fn finite_tank_roots_retain_each_inventory_as_the_body_translates_and_rotates() {
    let water = fluid(1000.0);
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), tank_geometry(), 16.0, volume_controls())
            .unwrap();
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let control = FloatingEquilibriumControls {
        maximum_heave_offset_m: 0.25,
        maximum_tilt_radians: 0.25,
        ..controls()
    };
    let result = solve_floating_equilibrium(
        &tank,
        &mesh,
        body.frame(),
        initial(body.frame(), 1.3, Vec3::new(0.07, -0.05, 0.0)),
        control,
        || false,
    )
    .unwrap();
    balanced(result, &tank, &mesh, body.frame());
    let resolved = body.world_from_body(result.state).unwrap();
    let fill = tank.fill_at_pose(&mesh, resolved, || false).unwrap();
    // Displacing 2 m3 raises the same 16 m3 inventory in a 16 m2 tank.
    near(resolved.translation().z, 1.125, 2e-9);
    near(fill.surface_height_m, 1.125, 2e-9);
    near(
        16.0 * fill.surface_height_m - fill.displaced.volume_m3(),
        16.0,
        2e-11,
    );
    assert_eq!(tank.liquid_mass_kg(), 16000.0);

    let light = fluid(800.0);
    let layered = FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: &water,
                liquid_volume_m3: 14.0,
            },
            FiniteHydrostaticLayer {
                source: &light,
                liquid_volume_m3: 30.0,
            },
        ],
        gravity(),
        tank_geometry(),
        volume_controls(),
        || false,
    )
    .unwrap();
    let source = super::body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let result = solve_floating_equilibrium(
        &layered,
        &mesh,
        body.frame(),
        initial(body.frame(), 1.1, Vec3::new(0.05, -0.04, 0.0)),
        control,
        || false,
    )
    .unwrap();
    balanced(result, &layered, &mesh, body.frame());
    let resolved = body.world_from_body(result.state).unwrap();
    let interfaces = layered
        .interfaces_at_pose(&mesh, resolved, || false)
        .unwrap();
    near(resolved.translation().z, 1.0, 2e-9);
    near(interfaces[0].surface_height_m, 1.0, 2e-9);
    near(interfaces[1].surface_height_m, 3.0, 2e-9);
    near(
        16.0 * interfaces[0].surface_height_m - interfaces[0].displaced_below.volume_m3(),
        14.0,
        2e-11,
    );
    near(
        16.0 * interfaces[1].surface_height_m - interfaces[1].displaced_below.volume_m3(),
        44.0,
        2e-11,
    );
    near(result.load.displaced_mass_kg(), 3600.0, 1e-7);
}

#[test]
fn common_world_transform_and_yaw_do_not_introduce_extra_equilibrium_constraints() {
    let water = fluid(1000.0);
    let sides = Vec3::new(3.0, 2.0, 1.0);
    let source = body(650.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let authored = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.08, -0.06, 0.4))).unwrap(),
        Vec3::new(0.3, -0.4, 0.12),
    )
    .unwrap();
    let world = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.4, -0.6, 0.3))).unwrap(),
        Vec3::new(3.0, -2.0, 1.0),
    )
    .unwrap();
    let mut solutions = Vec::new();
    for transform in [Se3::identity(), world] {
        let g = transform
            .rotation()
            .rotate(Vec3::new(0.0, 0.0, -10.0))
            .unwrap();
        let roll = transform
            .rotation()
            .rotate(Vec3::new(1.0, 0.0, 0.0))
            .unwrap();
        let medium = HydrostaticMedium::new(
            &water,
            Gravity::new(fs_mbd::Vec3::new(g.x, g.y, g.z)).unwrap(),
            transform.translation(),
        )
        .unwrap();
        let reference = body
            .state_from_kinematics(transform.compose(authored).unwrap(), Twist::zero())
            .unwrap();
        let result = solve_floating_equilibrium(
            &medium,
            &mesh,
            body.frame(),
            reference,
            FloatingEquilibriumControls {
                roll_direction_world: fs_mbd::Vec3::new(roll.x, roll.y, roll.z),
                ..controls()
            },
            || false,
        )
        .unwrap();
        balanced(result, &medium, &mesh, body.frame());
        solutions.push(
            transform
                .inverse()
                .unwrap()
                .compose(body.world_from_body(result.state).unwrap())
                .unwrap(),
        );
    }
    near(solutions[0].translation().z, -0.15, 2e-9);
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        let a = solutions[0].transform_point(axis).unwrap();
        let b = solutions[1].transform_point(axis).unwrap();
        near(a.x, b.x, 2e-8);
        near(a.y, b.y, 2e-8);
        near(a.z, b.z, 2e-8);
    }

    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let yaw = 0.73_f64;
    let reference = initial(body.frame(), 0.1, Vec3::new(0.0, 0.0, yaw));
    let result =
        solve_floating_equilibrium(&medium, &mesh, body.frame(), reference, controls(), || {
            false
        })
        .unwrap();
    let x = body
        .world_from_body(result.state)
        .unwrap()
        .rotation()
        .rotate(Vec3::new(1.0, 0.0, 0.0))
        .unwrap();
    near(x.x, yaw.cos(), 2e-9);
    near(x.y, yaw.sin(), 2e-9);
    near(x.z, 0.0, 2e-9);
    assert!(norm(result.tilt_rotation_vector_world_radians) < 2e-9);
}

#[test]
fn source_snapshots_cancellation_and_both_native_work_caps_remain_transactional() {
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
    let water = CountedFluid(Cell::new(0));
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let counted = Counted(&medium, Cell::new(0));
    let sides = Vec3::new(3.0, 2.0, 1.0);
    let source = body(650.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let reference = initial(body.frame(), 0.18, Vec3::new(0.1, -0.08, 0.0));
    let mut polls = 0;
    let expected =
        solve_floating_equilibrium(&counted, &mesh, body.frame(), reference, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(counted.1.get(), expected.evaluations);
    for stop in [1, 2, polls / 2, polls - 1, polls] {
        let mut seen = 0;
        assert_eq!(
            solve_floating_equilibrium(
                &counted,
                &mesh,
                body.frame(),
                reference,
                controls(),
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(FloatingEquilibriumError::Physics(
                HydrostaticError::Cancelled
            ))
        );
    }
    let cap = FloatingEquilibriumControls {
        maximum_evaluations: 1,
        ..controls()
    };
    assert!(matches!(
        solve_floating_equilibrium(&counted, &mesh, body.frame(), reference, cap, || false),
        Err(FloatingEquilibriumError::EvaluationBudget {
            evaluations: 1,
            maximum_evaluations: 1
        })
    ));
    let cap = FloatingEquilibriumControls {
        maximum_newton_iterations: 1,
        ..controls()
    };
    assert!(matches!(
        solve_floating_equilibrium(&counted, &mesh, body.frame(), reference, cap, || false),
        Err(FloatingEquilibriumError::NotConverged { iterations: 1, .. })
    ));
    assert_eq!(
        solve_floating_equilibrium(&medium, &mesh, body.frame(), reference, controls(), || {
            false
        })
        .unwrap(),
        expected
    );
    assert_eq!(water.0.get(), 1);

    counted.1.set(0);
    for invalid in [
        FloatingEquilibriumControls {
            maximum_force_defect_n: 0.0,
            ..controls()
        },
        FloatingEquilibriumControls {
            maximum_torque_defect_nm: f64::NAN,
            ..controls()
        },
        FloatingEquilibriumControls {
            length_scale_m: 0.0,
            ..controls()
        },
        FloatingEquilibriumControls {
            difference_step: 0.0,
            ..controls()
        },
        FloatingEquilibriumControls {
            maximum_heave_offset_m: 0.0,
            ..controls()
        },
        FloatingEquilibriumControls {
            maximum_tilt_radians: core::f64::consts::PI,
            ..controls()
        },
        FloatingEquilibriumControls {
            roll_direction_world: fs_mbd::Vec3::new(0.0, 0.0, 1.0),
            ..controls()
        },
        FloatingEquilibriumControls {
            maximum_newton_iterations: MAX_FLOATING_NEWTON_ITERATIONS + 1,
            ..controls()
        },
        FloatingEquilibriumControls {
            maximum_evaluations: MAX_FLOATING_EVALUATIONS + 1,
            ..controls()
        },
    ] {
        assert!(matches!(
            solve_floating_equilibrium(&counted, &mesh, body.frame(), reference, invalid, || false),
            Err(FloatingEquilibriumError::Physics(
                HydrostaticError::Invalid(_)
            ))
        ));
    }
    assert_eq!(counted.1.get(), 0);
}

#[test]
fn every_newton_trial_stays_in_the_pose_domain_and_an_existing_neutral_root_is_usable() {
    struct Bounded<'a> {
        medium: &'a HydrostaticMedium<'a, FluidTransportState>,
        reference: RigidBodyState,
        calls: Cell<usize>,
    }
    impl HydrostaticReservoir for Bounded<'_> {
        fn gravity(&self) -> Gravity {
            gravity()
        }
        fn load_frame(
            &self,
            mesh: &DisplacementMesh,
            body: &PrincipalFrame,
            state: RigidBodyState,
            cancelled: impl FnMut() -> bool,
        ) -> Result<HydrostaticLoad, HydrostaticError> {
            self.calls.set(self.calls.get() + 1);
            let a = self.reference.pose();
            let b = state.pose();
            assert!((a.position_world().z - b.position_world().z).abs() <= 0.010000000001);
            assert_eq!(a.position_world().x, b.position_world().x);
            assert_eq!(a.position_world().y, b.position_world().y);
            let dot: f64 = a
                .orientation()
                .components()
                .into_iter()
                .zip(b.orientation().components())
                .map(|(x, y)| x * y)
                .sum();
            assert!(2.0 * dot.abs().min(1.0).acos() <= 0.0100000001);
            HydrostaticReservoir::load_frame(self.medium, mesh, body, state, cancelled)
        }
    }
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let sides = Vec3::new(3.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let reference = initial(body.frame(), 0.2, origin());
    let bounded = Bounded {
        medium: &medium,
        reference,
        calls: Cell::new(0),
    };
    let tiny = FloatingEquilibriumControls {
        maximum_heave_offset_m: 0.01,
        maximum_tilt_radians: 0.01,
        maximum_evaluations: 256,
        ..controls()
    };
    assert!(matches!(
        solve_floating_equilibrium(&bounded, &mesh, body.frame(), reference, tiny, || false),
        Err(FloatingEquilibriumError::NotConverged { .. }
            | FloatingEquilibriumError::Difference(_)
            | FloatingEquilibriumError::EvaluationBudget { .. })
    ));
    assert!(bounded.calls.get() > 1 && bounded.calls.get() <= 256);

    let source = super::body(1000.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let reference = initial(body.frame(), -2.0, origin());
    let result = solve_floating_equilibrium(
        &medium,
        &mesh,
        body.frame(),
        reference,
        FloatingEquilibriumControls {
            maximum_evaluations: 1,
            ..controls()
        },
        || false,
    )
    .unwrap();
    assert_eq!(result.evaluations, 1);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.state, reference);
    // A neutral family is a valid force balance; the initializer does not claim stability.
}

#[test]
fn original_reservoir_refusals_and_unmodeled_force_components_cannot_be_hidden() {
    let water = fluid(1000.0);
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let tank =
        FiniteHydrostaticTank::new(&water, gravity(), tank_geometry(), 63.8, volume_controls())
            .unwrap();
    let reference = initial(body.frame(), 3.0, origin());
    assert!(matches!(
        solve_floating_equilibrium(&tank, &mesh, body.frame(), reference, controls(), || false),
        Err(FloatingEquilibriumError::Physics(
            HydrostaticError::TankOverflow { .. }
        ))
    ));
    assert_eq!(tank.liquid_volume_m3(), 63.8);

    // Deliberately inconsistent extension: report vertical gravity but return
    // another reservoir's oblique force. Three projected residuals could balance
    // while a lateral force remains; the final full-vector gate must refuse it.
    struct Inconsistent<'a>(&'a HydrostaticMedium<'a, FluidTransportState>);
    impl HydrostaticReservoir for Inconsistent<'_> {
        fn gravity(&self) -> Gravity {
            gravity()
        }
        fn load_frame(
            &self,
            mesh: &DisplacementMesh,
            body: &PrincipalFrame,
            state: RigidBodyState,
            cancelled: impl FnMut() -> bool,
        ) -> Result<HydrostaticLoad, HydrostaticError> {
            HydrostaticReservoir::load_frame(self.0, mesh, body, state, cancelled)
        }
    }
    let oblique = HydrostaticMedium::new(
        &water,
        Gravity::new(fs_mbd::Vec3::new(1.0, 0.0, -10.0)).unwrap(),
        origin(),
    )
    .unwrap();
    let reference = initial(body.frame(), 0.1, origin());
    assert!(
        matches!(solve_floating_equilibrium(&Inconsistent(&oblique), &mesh, body.frame(), reference, controls(), || false),
        Err(FloatingEquilibriumError::ResidualDefect { force_defect_n, .. }) if force_defect_n > 1900.0)
    );
}
