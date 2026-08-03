use fs_flux::{
    AlternativeWrenchSet, ApplicabilityEnvelope, BodyKinematics, ClosedRange, ContributionFamily,
    CorrelationIdentity, CorrelationUncertainty, DiscGeometry, DiscPose, EdgeFlow, FormDrag,
    GasProperties, GasPropertyCard, OrientationRateDamping, ReducedAeroComponents,
    ReducedAeroError, ReducedAeroInput, ReducedAeroModel, RotationalSkinFriction, SurfaceRoughness,
    Vec3, WorkWindow,
};

fn range(minimum: f64, maximum: f64) -> ClosedRange {
    ClosedRange::try_new(minimum, maximum).expect("test range is valid")
}

fn envelope() -> ApplicabilityEnvelope {
    ApplicabilityEnvelope {
        translational_reynolds: range(0.0, 1.0e9),
        rotational_reynolds: range(0.0, 1.0e9),
        relative_roughness: range(0.0, 0.1),
        maximum_tip_mach: 0.8,
    }
}

fn identity(id: &str) -> CorrelationIdentity {
    CorrelationIdentity::try_new(id, "v1", "source.reduced-aero.fixture")
        .expect("test correlation identity")
}

fn uncertainty() -> CorrelationUncertainty {
    CorrelationUncertainty {
        source_id: "source.reduced-aero.uncertainty.fixture".to_owned(),
        coefficient_relative_half_width: 0.15,
    }
}

fn model(id: &str, form_coefficient: f64) -> ReducedAeroModel {
    let components = ReducedAeroComponents {
        form_drag: Some(FormDrag {
            coefficient: form_coefficient,
        }),
        rotational_skin_friction: Some(RotationalSkinFriction { coefficient: 0.01 }),
        edge_flow: Some(EdgeFlow { coefficient: 0.04 }),
        orientation_rate_damping: Some(OrientationRateDamping { coefficient: 0.02 }),
    };
    ReducedAeroModel::try_new(
        identity(id),
        envelope(),
        uncertainty(),
        components,
        &[
            ContributionFamily::TranslationalFormDrag,
            ContributionFamily::RotationalSkinFriction,
            ContributionFamily::EdgeFlow,
            ContributionFamily::OrientationRateDamping,
        ],
    )
    .expect("admitted fixture model")
}

fn gas(velocity_world_m_per_s: Vec3) -> GasProperties {
    GasProperties::try_from(GasPropertyCard {
        source_id: "gas.air.fixture".to_owned(),
        density_kg_per_m3: Some(1.2),
        dynamic_viscosity_pa_s: Some(1.8e-5),
        speed_of_sound_m_per_s: Some(343.0),
        velocity_world_m_per_s,
    })
    .expect("complete gas card")
}

fn input(gas_velocity: Vec3) -> ReducedAeroInput {
    ReducedAeroInput {
        world_frame_id: "world.inertial".to_owned(),
        geometry: DiscGeometry {
            radius_m: 0.12,
            exterior_thickness_m: 0.008,
        },
        pose: DiscPose::try_new(Vec3::new(0.0, 0.6, 0.8)).expect("unit normal"),
        kinematics: BodyKinematics {
            reference_point_world_m: Vec3::new(1.0, -2.0, 0.5),
            linear_velocity_world_m_per_s: Vec3::new(4.0, -1.0, 2.0),
            angular_velocity_world_rad_per_s: Vec3::new(5.0, 15.0, 20.0),
        },
        gas: gas(gas_velocity),
        roughness: SurfaceRoughness {
            source_id: "roughness.fixture".to_owned(),
            height_m: 1.0e-5,
        },
    }
}

fn assert_vec_close(left: Vec3, right: Vec3, tolerance: f64) {
    assert!(
        (left.x - right.x).abs() <= tolerance,
        "x: {left:?} vs {right:?}"
    );
    assert!(
        (left.y - right.y).abs() <= tolerance,
        "y: {left:?} vs {right:?}"
    );
    assert!(
        (left.z - right.z).abs() <= tolerance,
        "z: {left:?} vs {right:?}"
    );
}

#[test]
fn aero_001_zero_density_and_zero_motion_are_exact_zero_limits() {
    let model = model("fixture.zero", 1.1);
    let mut vacuum = input(Vec3::ZERO);
    vacuum.gas = GasProperties::try_from(GasPropertyCard {
        source_id: "gas.vacuum.fixture".to_owned(),
        density_kg_per_m3: Some(0.0),
        dynamic_viscosity_pa_s: Some(1.8e-5),
        speed_of_sound_m_per_s: Some(343.0),
        velocity_world_m_per_s: Vec3::ZERO,
    })
    .expect("vacuum limit is admitted");
    let vacuum_wrench = model.evaluate(&vacuum).expect("vacuum evaluates");
    assert_eq!(vacuum_wrench.force_world_n, Vec3::ZERO);
    assert_eq!(vacuum_wrench.torque_world_n_m, Vec3::ZERO);
    assert_eq!(vacuum_wrench.receipt.relative_power_w, 0.0);

    let mut stationary = input(Vec3::ZERO);
    stationary.kinematics = BodyKinematics {
        reference_point_world_m: stationary.kinematics.reference_point_world_m,
        linear_velocity_world_m_per_s: Vec3::ZERO,
        angular_velocity_world_rad_per_s: Vec3::ZERO,
    };
    let stationary_wrench = model.evaluate(&stationary).expect("static limit evaluates");
    assert_eq!(stationary_wrench.force_world_n, Vec3::ZERO);
    assert_eq!(stationary_wrench.torque_world_n_m, Vec3::ZERO);
    assert_eq!(stationary_wrench.receipt.dissipated_relative_power_w, 0.0);
}

#[test]
fn aero_002_reversal_and_passive_stationary_power_laws_hold() {
    let model = model("fixture.reversal", 1.1);
    let forward = model.evaluate(&input(Vec3::ZERO)).expect("forward result");
    let mut reverse_input = input(Vec3::ZERO);
    reverse_input.kinematics.linear_velocity_world_m_per_s = reverse_input
        .kinematics
        .linear_velocity_world_m_per_s
        .scaled(-1.0);
    reverse_input.kinematics.angular_velocity_world_rad_per_s = reverse_input
        .kinematics
        .angular_velocity_world_rad_per_s
        .scaled(-1.0);
    let reverse = model.evaluate(&reverse_input).expect("reverse result");

    assert_vec_close(
        reverse.force_world_n,
        forward.force_world_n.scaled(-1.0),
        1.0e-13,
    );
    assert_vec_close(
        reverse.torque_world_n_m,
        forward.torque_world_n_m.scaled(-1.0),
        1.0e-13,
    );
    assert!(forward.receipt.relative_power_w <= 1.0e-12);
    assert!(forward.receipt.dissipated_relative_power_w >= -1.0e-12);
    assert!((forward.receipt.body_power_w - forward.receipt.relative_power_w).abs() <= 1.0e-12);
}

#[test]
fn aero_003_scaling_is_quadratic_and_source_mutation_is_retained() {
    let model = model("fixture.scale", 1.1);
    let base = model.evaluate(&input(Vec3::ZERO)).expect("base result");
    let mut doubled = input(Vec3::ZERO);
    doubled.kinematics.linear_velocity_world_m_per_s =
        doubled.kinematics.linear_velocity_world_m_per_s.scaled(2.0);
    doubled.kinematics.angular_velocity_world_rad_per_s = doubled
        .kinematics
        .angular_velocity_world_rad_per_s
        .scaled(2.0);
    let doubled = model.evaluate(&doubled).expect("doubled result");
    assert_vec_close(
        doubled.force_world_n,
        base.force_world_n.scaled(4.0),
        1.0e-12,
    );
    assert_vec_close(
        doubled.torque_world_n_m,
        base.torque_world_n_m.scaled(4.0),
        1.0e-12,
    );

    let mut source_mutated = input(Vec3::ZERO);
    source_mutated.gas = gas(Vec3::ZERO);
    source_mutated.roughness.source_id = "roughness.mutated".to_owned();
    let source_mutated = model
        .evaluate(&source_mutated)
        .expect("mutated source evaluates");
    assert_ne!(
        base.receipt.roughness_source_id,
        source_mutated.receipt.roughness_source_id
    );
    assert_eq!(base.receipt.coefficient_relative_half_width, 0.15);
    assert_eq!(
        base.receipt.correlation_uncertainty_source_id,
        "source.reduced-aero.uncertainty.fixture"
    );
}

#[test]
fn aero_004_moving_ambient_uses_relative_motion_and_demands_total_accounting() {
    let model = model("fixture.ambient", 1.1);
    let still_air = model
        .evaluate(&input(Vec3::ZERO))
        .expect("still air result");
    let matching_ambient = model
        .evaluate(&input(Vec3::new(4.0, -1.0, 2.0)))
        .expect("matching ambient result");

    assert_eq!(matching_ambient.components.form_force_world_n, Vec3::ZERO);
    assert!(
        matching_ambient
            .receipt
            .moving_ambient_requires_total_energy_accounting
    );
    assert!(still_air.receipt.body_power_w <= 1.0e-12);
    assert!(matching_ambient.receipt.relative_power_w <= 1.0e-12);
    assert!(
        (matching_ambient.receipt.body_power_w
            - (matching_ambient.receipt.relative_power_w
                + matching_ambient.receipt.ambient_boundary_power_w))
            .abs()
            <= 1.0e-12
    );
}

#[test]
fn aero_005_frame_objectivity_under_a_rigid_quarter_turn() {
    let model = model("fixture.objective", 1.1);
    let original = model.evaluate(&input(Vec3::ZERO)).expect("original result");
    let rotate = |vector: Vec3| Vec3::new(-vector.y, vector.x, vector.z);
    let mut rotated = input(Vec3::ZERO);
    rotated.pose = DiscPose::try_new(rotate(rotated.pose.normal_world)).expect("rotated normal");
    rotated.kinematics.reference_point_world_m = rotate(rotated.kinematics.reference_point_world_m);
    rotated.kinematics.linear_velocity_world_m_per_s =
        rotate(rotated.kinematics.linear_velocity_world_m_per_s);
    rotated.kinematics.angular_velocity_world_rad_per_s =
        rotate(rotated.kinematics.angular_velocity_world_rad_per_s);
    let rotated = model.evaluate(&rotated).expect("rotated result");

    assert_vec_close(
        rotated.force_world_n,
        rotate(original.force_world_n),
        1.0e-13,
    );
    assert_vec_close(
        rotated.torque_world_n_m,
        rotate(original.torque_world_n_m),
        1.0e-13,
    );
    assert!(
        (rotated.receipt.relative_power_w - original.receipt.relative_power_w).abs() <= 1.0e-12
    );
}

#[test]
fn aero_006_domain_and_missing_property_refusals_are_typed() {
    let narrow = ApplicabilityEnvelope {
        translational_reynolds: range(0.0, 1.0),
        ..envelope()
    };
    let components = ReducedAeroComponents {
        form_drag: Some(FormDrag { coefficient: 1.0 }),
        ..ReducedAeroComponents::default()
    };
    let narrow_model = ReducedAeroModel::try_new(
        identity("fixture.narrow"),
        narrow,
        uncertainty(),
        components,
        &[ContributionFamily::TranslationalFormDrag],
    )
    .expect("narrow model is structurally valid");
    assert!(matches!(
        narrow_model.evaluate(&input(Vec3::ZERO)),
        Err(ReducedAeroError::OutsideCorrelationDomain {
            quantity: "translational_reynolds",
            ..
        })
    ));
    assert!(matches!(
        GasProperties::try_from(GasPropertyCard {
            source_id: "gas.incomplete".to_owned(),
            density_kg_per_m3: None,
            dynamic_viscosity_pa_s: Some(1.8e-5),
            speed_of_sound_m_per_s: Some(343.0),
            velocity_world_m_per_s: Vec3::ZERO
        }),
        Err(ReducedAeroError::MissingGasProperty("density_kg_per_m3"))
    ));
}

#[test]
fn aero_007_alternatives_remain_distinct_and_deterministically_sorted() {
    let alternatives = AlternativeWrenchSet::evaluate(
        &[model("zeta", 1.1), model("alpha", 0.6)],
        &input(Vec3::ZERO),
    )
    .expect("alternatives evaluate");
    assert_eq!(alternatives.candidates[0].correlation.id, "alpha");
    assert!(alternatives.has_force_or_torque_disagreement());
}

#[test]
fn aero_008_exactly_once_work_refuses_duplicate_application() {
    let wrench = model("fixture.work", 1.1)
        .evaluate(&input(Vec3::ZERO))
        .expect("candidate result");
    let mut window = WorkWindow::default();
    let receipt = window
        .record_once(7, 0.25, &wrench)
        .expect("first application succeeds");
    assert_eq!(receipt.body_work_j, wrench.receipt.body_power_w * 0.25);
    assert_eq!(
        window.relative_dissipation_j(),
        wrench.receipt.dissipated_relative_power_w * 0.25
    );
    assert!(matches!(
        window.record_once(7, 0.25, &wrench),
        Err(ReducedAeroError::DuplicateWorkExchange { key: 7 })
    ));
}

#[test]
fn aero_009_thin_gap_and_target_fits_are_hostile_admission_refusals() {
    let components = ReducedAeroComponents {
        form_drag: Some(FormDrag { coefficient: 1.0 }),
        ..ReducedAeroComponents::default()
    };
    for forbidden in [
        ContributionFamily::ThinGapPressure,
        ContributionFamily::TargetFitted,
    ] {
        assert!(
            matches!(ReducedAeroModel::try_new(identity("fixture.hostile"), envelope(), uncertainty(), components, &[ContributionFamily::TranslationalFormDrag, forbidden]), Err(ReducedAeroError::ForbiddenContribution { family }) if family == forbidden)
        );
    }
}

#[test]
fn aero_010_mutated_public_pose_is_revalidated_at_evaluation() {
    let model = model("fixture.hostile-pose", 1.0);
    let mut nonfinite = input(Vec3::ZERO);
    nonfinite.pose.normal_world = Vec3::new(f64::NAN, 0.0, 1.0);
    assert!(matches!(
        model.evaluate(&nonfinite),
        Err(ReducedAeroError::InvalidInput {
            field: "pose.normal_world"
        })
    ));

    let mut nonunit = input(Vec3::ZERO);
    nonunit.pose.normal_world = Vec3::new(0.0, 0.0, 2.0);
    assert!(matches!(
        model.evaluate(&nonunit),
        Err(ReducedAeroError::NonUnitDiscAxis { .. })
    ));
}

#[test]
fn aero_011_public_identity_and_range_literals_are_revalidated_at_admission() {
    let components = ReducedAeroComponents {
        form_drag: Some(FormDrag { coefficient: 1.0 }),
        ..ReducedAeroComponents::default()
    };
    let malformed_identity = CorrelationIdentity {
        id: "".to_owned(),
        version: "v1".to_owned(),
        source_id: "source.reduced-aero.fixture".to_owned(),
    };
    assert!(matches!(
        ReducedAeroModel::try_new(
            malformed_identity,
            envelope(),
            uncertainty(),
            components,
            &[ContributionFamily::TranslationalFormDrag],
        ),
        Err(ReducedAeroError::InvalidIdentity {
            field: "correlation.id"
        })
    ));

    for malformed_range in [
        ClosedRange {
            minimum: f64::NAN,
            maximum: 1.0,
        },
        ClosedRange {
            minimum: 2.0,
            maximum: 1.0,
        },
    ] {
        let malformed_envelope = ApplicabilityEnvelope {
            translational_reynolds: malformed_range,
            ..envelope()
        };
        assert!(matches!(
            ReducedAeroModel::try_new(
                identity("fixture.malformed-range"),
                malformed_envelope,
                uncertainty(),
                components,
                &[ContributionFamily::TranslationalFormDrag],
            ),
            Err(ReducedAeroError::InvalidInput {
                field: "correlation.range"
            })
        ));
    }
}

#[test]
fn aero_012_duplicate_correlation_identity_is_order_independent_refusal() {
    let first = model("fixture.duplicate", 0.2);
    let second = model("fixture.duplicate", 1.7);
    for models in [vec![first.clone(), second.clone()], vec![second, first]] {
        assert!(matches!(
            AlternativeWrenchSet::evaluate(&models, &input(Vec3::ZERO)),
            Err(ReducedAeroError::DuplicateCorrelationIdentity { correlation })
                if correlation.id == "fixture.duplicate"
        ));
    }
}

#[test]
fn aero_013_edge_on_form_drag_uses_rim_silhouette_not_wetted_rim_area() {
    let model = model("fixture.edge-on-silhouette", 1.1);
    let mut edge_on = input(Vec3::ZERO);
    edge_on.pose = DiscPose::try_new(Vec3::new(0.0, 0.0, 1.0)).expect("unit normal");
    edge_on.kinematics.linear_velocity_world_m_per_s = Vec3::new(4.0, 0.0, 0.0);
    edge_on.kinematics.angular_velocity_world_rad_per_s = Vec3::ZERO;
    let wrench = model.evaluate(&edge_on).expect("edge-on wrench");

    let rim_silhouette_m2 = 2.0 * edge_on.geometry.radius_m * edge_on.geometry.exterior_thickness_m;
    let roughness_factor = 1.0 + edge_on.roughness.height_m / edge_on.geometry.radius_m;
    let expected_force_n =
        -0.5 * 1.2 * 1.1 * roughness_factor * rim_silhouette_m2 * 4.0_f64.powi(2);
    assert!((wrench.force_world_n.x - expected_force_n).abs() <= 1.0e-14);
    assert_eq!(wrench.force_world_n.y, 0.0);
    assert_eq!(wrench.force_world_n.z, 0.0);
}

#[test]
fn aero_014_forged_wrenches_and_overflow_do_not_consume_work_keys_or_totals() {
    let wrench = model("fixture.transaction", 1.1)
        .evaluate(&input(Vec3::ZERO))
        .expect("candidate result");
    let mut window = WorkWindow::default();

    let mut forged_vector = wrench.clone();
    forged_vector.force_world_n.x = f64::NAN;
    assert!(matches!(
        window.record_once(29, 0.25, &forged_vector),
        Err(ReducedAeroError::InvalidCandidateWrench {
            field: "candidate.force_world_n"
        })
    ));
    assert_eq!(window.body_work_j(), 0.0);
    assert_eq!(window.relative_dissipation_j(), 0.0);

    let mut overflowing_work = wrench.clone();
    overflowing_work.receipt.relative_power_w = 0.0;
    overflowing_work.receipt.dissipated_relative_power_w = 0.0;
    overflowing_work.receipt.body_power_w = f64::MAX;
    overflowing_work.receipt.ambient_boundary_power_w = f64::MAX;
    assert!(matches!(
        window.record_once(29, 2.0, &overflowing_work),
        Err(ReducedAeroError::NonFiniteDerived {
            field: "work.body_work_j"
        })
    ));
    assert_eq!(window.body_work_j(), 0.0);
    assert_eq!(window.relative_dissipation_j(), 0.0);

    window
        .record_once(29, 0.25, &wrench)
        .expect("valid candidate must retain the unconsumed key");
}

#[test]
fn aero_015_derived_overflow_and_passivity_are_checked_in_all_build_modes() {
    let model = model("fixture.overflow", 1.1);
    let mut overflowing = input(Vec3::ZERO);
    overflowing.kinematics.linear_velocity_world_m_per_s = Vec3::new(f64::MAX, 0.0, 0.0);
    assert!(matches!(
        model.evaluate(&overflowing),
        Err(ReducedAeroError::NonFiniteDerived {
            field: "translational_reynolds"
        })
    ));

    let wrench = model
        .evaluate(&input(Vec3::ZERO))
        .expect("passive candidate");
    assert!(wrench.receipt.relative_power_w <= 0.0);
    let mut forged_passivity = wrench.clone();
    forged_passivity.receipt.relative_power_w = 1.0;
    forged_passivity.receipt.dissipated_relative_power_w = -1.0;
    forged_passivity.receipt.body_power_w = 1.0;
    forged_passivity.receipt.ambient_boundary_power_w = 0.0;
    let mut window = WorkWindow::default();
    assert!(matches!(
        window.record_once(31, 0.25, &forged_passivity),
        Err(ReducedAeroError::PassivePowerViolation {
            relative_power_w: 1.0
        })
    ));
    window
        .record_once(31, 0.25, &wrench)
        .expect("passivity refusal must not consume the key");
}
