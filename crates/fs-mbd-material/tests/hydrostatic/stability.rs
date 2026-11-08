//! Independent waterplane, metacentric and stratified stiffness references.
use super::*;
use fs_mbd_material::hydrostatic::{
    FiniteHydrostaticTank, FiniteTankControls, HydrostaticReservoir, HydrostaticStability,
    HydrostaticStabilityControls, HydrostaticStabilityError, HydrostaticTankGeometry,
    LinearDensityHydrostaticMedium, analyze_hydrostatic_stability,
};

fn controls() -> HydrostaticStabilityControls {
    HydrostaticStabilityControls {
        roll_direction_world: fs_mbd::Vec3::new(1.0, 0.0, 0.0),
        length_scale_m: 1.0,
        heave_step_m: 1e-4,
        tilt_step_radians: 1e-4,
        maximum_force_defect_n: 1e-6,
        maximum_torque_defect_nm: 1e-6,
        maximum_relative_stiffness_defect: 1e-4,
        neutral_stiffness_tolerance_j: 1e-3,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}

#[test]
fn waterplane_and_metacentric_stiffness_have_correct_units_and_axes() {
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    for length_scale_m in [0.5, 1.0, 2.0] {
        let result = analyze_hydrostatic_stability(
            &medium,
            &mesh,
            body.frame(),
            state,
            HydrostaticStabilityControls {
                length_scale_m,
                ..controls()
            },
            || false,
        )
        .unwrap();
        // rho*g*A; roll/pitch = mass*g*(I_waterplane/V_displaced - BG).
        // V_displaced=3, BG=.25; I_roll=4.5, I_pitch=2 m4.
        let expected = [60000.0 * length_scale_m.powi(2), 37500.0, 12500.0];
        for (i, row) in result.restoring_stiffness_j.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                near(*value, if i == j { expected[i] } else { 0.0 }, 0.002);
            }
        }
        assert_eq!(result.classification, HydrostaticStability::Restoring);
        assert_eq!(result.evaluations, 13);
        assert!(result.estimated_resolution_j < 0.01);
    }
}

#[test]
fn force_balanced_tall_body_is_unstable_and_neutral_submarine_is_marginal() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    for (sides, density, z, expected) in [
        (
            Vec3::new(1.0, 1.0, 4.0),
            500.0,
            0.0,
            HydrostaticStability::Unstable,
        ),
        (
            Vec3::new(2.0, 2.0, 2.0),
            1000.0,
            -3.0,
            HydrostaticStability::Marginal,
        ),
    ] {
        let source = body(density, sides);
        let body = source.principal_body(1e-12, || false).unwrap();
        let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
        let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
        let result =
            analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, controls(), || {
                false
            })
            .unwrap();
        assert_eq!(result.classification, expected);
        if expected == HydrostaticStability::Unstable {
            near(result.eigenvalues_j[0], 20000.0 * (1.0 / 24.0 - 1.0), 0.002);
            near(result.eigenvalues_j[1], result.eigenvalues_j[0], 0.002);
        } else {
            assert!(
                result
                    .eigenvalues_j
                    .iter()
                    .all(|x| x.abs() <= result.estimated_resolution_j)
            );
        }
    }
}

#[test]
fn conserved_tank_waterline_increases_heave_stiffness() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        HydrostaticTankGeometry {
            world_from_tank: Se3::identity(),
            inside_dimensions_m: Vec3::new(5.0, 5.0, 5.0),
        },
        25.0,
        FiniteTankControls {
            maximum_volume_defect_m3: 1e-12,
            maximum_volume_evaluations: 64,
        },
    )
    .unwrap();
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(28.0 / 25.0), Twist::zero())
        .unwrap();
    let result =
        analyze_hydrostatic_stability(&tank, &mesh, body.frame(), state, controls(), || false)
            .unwrap();
    // Differentiate A_tank*h - A_body*(h-z+.5) = V_inventory.
    near(
        result.restoring_stiffness_j[0][0],
        10000.0 * 6.0 * 25.0 / 19.0,
        0.003,
    );
    near(result.restoring_stiffness_j[1][1], 37500.0, 0.003);
    near(result.restoring_stiffness_j[2][2], 12500.0, 0.003);
    assert_eq!(result.classification, HydrostaticStability::Restoring);
    close(tank.liquid_volume_m3(), 25.0);
}

#[test]
fn continuous_density_gradient_restores_fully_submerged_body() {
    let water = fluid(1000.0);
    let medium = LinearDensityHydrostaticMedium::new(&water, 400.0, gravity(), origin()).unwrap();
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(1800.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(-2.0), Twist::zero())
        .unwrap();
    let result =
        analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, controls(), || false)
            .unwrap();
    // Integrate rho(z)=1000-400*z over the FULL envelope. Heave = g*beta*V.
    // Tilt = g*beta*V*(horizontal_side^2 - vertical_side^2)/12.
    for (i, expected) in [24000.0, 16000.0, 6000.0].into_iter().enumerate() {
        near(result.restoring_stiffness_j[i][i], expected, 0.003);
    }
    assert_eq!(result.classification, HydrostaticStability::Restoring);
    close(result.load.displaced().volume_m3(), 6.0);
}

#[test]
fn rotating_world_gravity_and_heading_preserves_stiffness() {
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let baseline =
        analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, controls(), || false)
            .unwrap();
    let world = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.7, -0.4, 0.3))).unwrap(),
        Vec3::new(3.0, -4.0, 5.0),
    )
    .unwrap();
    let g = world.transform_vector(Vec3::new(0.0, 0.0, -10.0)).unwrap();
    let roll = world.transform_vector(Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let medium = HydrostaticMedium::new(
        &water,
        Gravity::new(fs_mbd::Vec3::new(g.x, g.y, g.z)).unwrap(),
        world.transform_point(origin()).unwrap(),
    )
    .unwrap();
    let state = body.state_from_kinematics(world, Twist::zero()).unwrap();
    let result = analyze_hydrostatic_stability(
        &medium,
        &mesh,
        body.frame(),
        state,
        HydrostaticStabilityControls {
            roll_direction_world: fs_mbd::Vec3::new(roll.x, roll.y, roll.z),
            ..controls()
        },
        || false,
    )
    .unwrap();
    assert_eq!(result.classification, baseline.classification);
    for i in 0..3 {
        for j in 0..3 {
            near(
                result.restoring_stiffness_j[i][j],
                baseline.restoring_stiffness_j[i][j],
                0.003,
            );
        }
    }
}

#[test]
fn original_source_is_frozen_and_every_cancellation_returns_no_report() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get() + 1);
            Ok(fluid(1000.0))
        }
    }
    let water = Counted(std::cell::Cell::new(0));
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let mut polls = 0;
    analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert!(matches!(
            analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, controls(), || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticStabilityError::Load(HydrostaticError::Cancelled))
        ));
    }
    assert_eq!(water.0.get(), 1);
    assert_eq!(state.linear_momentum_world(), fs_mbd::Vec3::ZERO);
}

#[test]
fn rejects_unbalanced_moving_unresolved_or_unrepresentable_inputs() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let off = body
        .state_from_kinematics(pose(0.1), Twist::zero())
        .unwrap();
    assert!(matches!(
        analyze_hydrostatic_stability(&medium, &mesh, body.frame(), off, controls(), || false),
        Err(HydrostaticStabilityError::NotEquilibrium { .. })
    ));
    let moving = RigidBodyState::new(
        state.pose(),
        fs_mbd::Vec3::new(0.0, 0.0, 1.0),
        fs_mbd::Vec3::ZERO,
    )
    .unwrap();
    assert!(matches!(
        analyze_hydrostatic_stability(&medium, &mesh, body.frame(), moving, controls(), || false),
        Err(HydrostaticStabilityError::Invalid(_))
    ));
    for c in [
        HydrostaticStabilityControls {
            heave_step_m: f64::NAN,
            ..controls()
        },
        HydrostaticStabilityControls {
            length_scale_m: 0.0,
            ..controls()
        },
        HydrostaticStabilityControls {
            roll_direction_world: fs_mbd::Vec3::new(0.0, 0.0, 1.0),
            ..controls()
        },
        HydrostaticStabilityControls {
            tilt_step_radians: 0.2,
            maximum_relative_stiffness_defect: 1e-8,
            ..controls()
        },
    ] {
        assert!(
            analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, c, || false)
                .is_err()
        );
    }
}

#[test]
fn finite_tank_failed_perturbation_is_original_error_without_state_change() {
    let water = fluid(1000.0);
    let tank = FiniteHydrostaticTank::new(
        &water,
        gravity(),
        HydrostaticTankGeometry {
            world_from_tank: Se3::identity(),
            inside_dimensions_m: Vec3::new(2.1, 3.1, 3.0),
        },
        6.51,
        FiniteTankControls {
            maximum_volume_defect_m3: 1e-12,
            maximum_volume_evaluations: 64,
        },
    )
    .unwrap();
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let state = body
        .state_from_kinematics(pose((6.51 + 3.0) / 6.51), Twist::zero())
        .unwrap();
    let original = tank
        .load_frame(&mesh, body.frame(), state, || false)
        .unwrap();
    let error = analyze_hydrostatic_stability(
        &tank,
        &mesh,
        body.frame(),
        state,
        HydrostaticStabilityControls {
            tilt_step_radians: 0.25,
            ..controls()
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        HydrostaticStabilityError::Load(HydrostaticError::Invalid(_))
    ));
    assert_eq!(
        tank.load_frame(&mesh, body.frame(), state, || false)
            .unwrap(),
        original
    );
}

#[test]
fn oblique_heave_cannot_silently_lose_a_world_coordinate() {
    // Keep the actual clipping calculation in a local frame. Force and COM
    // torque are translation invariant; this isolates pose quantization from
    // any loss of geometric precision at a large global origin.
    struct Translated<'a, R> {
        original: &'a R,
        offset: fs_mbd::Vec3,
    }
    impl<R: HydrostaticReservoir> HydrostaticReservoir for Translated<'_, R> {
        fn gravity(&self) -> Gravity {
            self.original.gravity()
        }
        fn load_frame(
            &self,
            mesh: &DisplacementMesh,
            body: &fs_mbd_material::principal::PrincipalFrame,
            state: RigidBodyState,
            cancelled: impl FnMut() -> bool,
        ) -> Result<fs_mbd_material::hydrostatic::HydrostaticLoad, HydrostaticError> {
            let local = RigidBodyState::new(
                Pose::new(
                    state.pose().position_world().sub(self.offset),
                    state.pose().orientation(),
                )?,
                state.linear_momentum_world(),
                state.angular_momentum_body(),
            )?;
            self.original.load_frame(mesh, body, local, cancelled)
        }
    }
    let water = fluid(1000.0);
    let gravity = Gravity::new(fs_mbd::Vec3::new(-6.0, 0.0, -8.0)).unwrap();
    let medium = HydrostaticMedium::new(&water, gravity, origin()).unwrap();
    let offset = fs_mbd::Vec3::new(2.0_f64.powi(52), 0.0, 0.0);
    let translated = Translated {
        original: &medium,
        offset,
    };
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let world = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.0, 0.6_f64.asin(), 0.0))).unwrap(),
        Vec3::new(offset.x, offset.y, offset.z),
    )
    .unwrap();
    let state = body.state_from_kinematics(world, Twist::zero()).unwrap();
    let error = analyze_hydrostatic_stability(
        &translated,
        &mesh,
        body.frame(),
        state,
        HydrostaticStabilityControls {
            roll_direction_world: fs_mbd::Vec3::new(0.8, 0.0, -0.6),
            ..controls()
        },
        || false,
    )
    .unwrap_err();
    // h and h/2 would both lose their x part while retaining z. Comparing only
    // stencil levels would accept 25600 N/m in place of the physical40000 N/m.
    assert_eq!(
        error,
        HydrostaticStabilityError::Invalid(
            "heave perturbation is not resolved in every world component"
        )
    );
}

#[test]
fn nearly_parallel_heading_still_returns_orthonormal_gravity_axes() {
    let water = fluid(1000.0);
    let g = Gravity::new(fs_mbd::Vec3::new(-6.0, 0.0, -8.0)).unwrap();
    let medium = HydrostaticMedium::new(&water, g, origin()).unwrap();
    let sides = Vec3::new(2.0, 2.0, 1.0);
    let source = body(500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let world = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.0, 0.6_f64.asin(), 0.0))).unwrap(),
        origin(),
    )
    .unwrap();
    let state = body.state_from_kinematics(world, Twist::zero()).unwrap();
    let report = analyze_hydrostatic_stability(
        &medium,
        &mesh,
        body.frame(),
        state,
        HydrostaticStabilityControls {
            roll_direction_world: fs_mbd::Vec3::new(0.6 + 0.8e-12, 0.0, 0.8 - 0.6e-12),
            ..controls()
        },
        || false,
    )
    .unwrap();
    for i in 0..3 {
        for j in 0..3 {
            near(
                report.axes_world[i].dot(report.axes_world[j]),
                if i == j { 1.0 } else { 0.0 },
                1e-14,
            );
        }
    }
    assert_eq!(report.classification, HydrostaticStability::Restoring);
}

#[test]
fn exact_waterline_contact_does_not_fake_smooth_restoring_stiffness() {
    let water = fluid(1000.0);
    let medium = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let sides = Vec3::new(2.0, 3.0, 1.0);
    let source = body(1000.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let state = body
        .state_from_kinematics(pose(-0.5), Twist::zero())
        .unwrap();
    // Downward heave stays fully immersed and neutral; upward heave loses
    // displaced mass with slope rho*g*A. Both central stencil sizes alone
    // report the same positive30000 N/m, concealing the derivative jump.
    match analyze_hydrostatic_stability(&medium, &mesh, body.frame(), state, controls(), || false) {
        Err(HydrostaticStabilityError::UnresolvedStiffness {
            one_sided_defect_j, ..
        }) => {
            assert!(one_sided_defect_j > 59000.0);
        }
        result => panic!("waterline contact must remain unresolved: {result:?}"),
    }
}
