#[path = "../src/external_air.rs"]
mod external_air;

use external_air::{
    EulerDiscBodyFrame, EulerDiscExteriorGeometry, EulerDiscExteriorState, EulerExternalAirInput,
    EulerExternalAirWorkWindow, ExteriorAirHeatDisposition, ExteriorAirPressure,
    ExteriorAirPressureScaling, ExternalAirDomain, ExternalAirError, ExternalAirIdentity,
    evaluate_euler_disc_external_air,
};
use fs_flux::{
    ApplicabilityEnvelope, ClosedRange, ContributionFamily, CorrelationIdentity,
    CorrelationUncertainty, DiscGeometry, EdgeFlow, FormDrag, GasPropertyCard,
    OrientationRateDamping, ReducedAeroComponents, ReducedAeroError, ReducedAeroModel,
    RotationalSkinFriction, SurfaceRoughness, Vec3,
};

fn range(minimum: f64, maximum: f64) -> ClosedRange {
    ClosedRange::try_new(minimum, maximum).expect("valid inclusive range")
}

fn model(id: &str, form: f64, spin: f64) -> ReducedAeroModel {
    let components = ReducedAeroComponents {
        form_drag: (form > 0.0).then_some(FormDrag { coefficient: form }),
        rotational_skin_friction: (spin > 0.0)
            .then_some(RotationalSkinFriction { coefficient: spin }),
        edge_flow: (spin > 0.0).then_some(EdgeFlow {
            coefficient: spin * 0.5,
        }),
        orientation_rate_damping: (spin > 0.0).then_some(OrientationRateDamping {
            coefficient: spin * 0.25,
        }),
    };
    let mut families = Vec::new();
    if form > 0.0 {
        families.push(ContributionFamily::TranslationalFormDrag);
    }
    if spin > 0.0 {
        families.extend([
            ContributionFamily::RotationalSkinFriction,
            ContributionFamily::EdgeFlow,
            ContributionFamily::OrientationRateDamping,
        ]);
    }
    ReducedAeroModel::try_new(
        CorrelationIdentity::try_new(id, "v1", format!("source:{id}")).expect("identity"),
        ApplicabilityEnvelope {
            translational_reynolds: range(0.0, 1.0e9),
            rotational_reynolds: range(0.0, 1.0e9),
            relative_roughness: range(0.0, 1.0),
            maximum_tip_mach: 1.0,
        },
        CorrelationUncertainty {
            source_id: format!("uncertainty:{id}"),
            coefficient_relative_half_width: 0.1,
        },
        components,
        &families,
    )
    .expect("admitted generic exterior correlation")
}

fn fixture() -> EulerExternalAirInput {
    EulerExternalAirInput {
        domain: ExternalAirDomain::ExteriorFreeGas,
        identity: ExternalAirIdentity {
            case_id: "euler-case".into(),
            world_frame_id: "world".into(),
            body_frame_id: "disc-body".into(),
            geometry_source_id: "geometry:disc".into(),
            state_source_id: "state:rigid".into(),
            domain_source_id: "domain:free-exterior-gas".into(),
        },
        geometry: EulerDiscExteriorGeometry {
            radius_m: 0.04,
            exterior_thickness_m: 0.003,
        },
        state: EulerDiscExteriorState {
            center_world_m: Vec3::ZERO,
            center_velocity_world_m_per_s: Vec3::new(4.0, 0.0, 0.0),
            angular_velocity_world_rad_per_s: Vec3::new(0.0, 0.0, 60.0),
            body_frame: EulerDiscBodyFrame {
                x_world: Vec3::new(1.0, 0.0, 0.0),
                z_world: Vec3::new(0.0, 0.0, 1.0),
            },
        },
        gas: GasPropertyCard {
            source_id: "gas:fixture".into(),
            density_kg_per_m3: Some(1.2),
            dynamic_viscosity_pa_s: Some(1.8e-5),
            speed_of_sound_m_per_s: Some(340.0),
            velocity_world_m_per_s: Vec3::ZERO,
        },
        pressure: ExteriorAirPressure {
            absolute_pressure_pa: 101_325.0,
            source_id: "pressure:fixture".into(),
        },
        exterior_roughness: SurfaceRoughness {
            source_id: "roughness:fixture".into(),
            height_m: 1.0e-5,
        },
        alternatives: vec![model("corr-a", 1.0, 0.02)],
    }
}

fn candidate(input: &EulerExternalAirInput) -> external_air::EulerExternalAirCandidate {
    evaluate_euler_disc_external_air(input)
        .expect("admitted exterior request")
        .candidates
        .into_iter()
        .next()
        .expect("one candidate")
}

fn close(left: f64, right: f64) {
    let scale = left.abs().max(right.abs()).max(1.0);
    assert!((left - right).abs() <= 1.0e-11 * scale, "{left} != {right}");
}

fn rotated_z(vector: Vec3) -> Vec3 {
    Vec3::new(-vector.y, vector.x, vector.z)
}

#[test]
fn g0_zero_density_and_zero_relative_speed_produce_no_wrench() {
    let mut vacuum = fixture();
    vacuum.gas.density_kg_per_m3 = Some(0.0);
    let result = candidate(&vacuum);
    assert_eq!(result.world_wrench.force_world_n, Vec3::ZERO);
    assert_eq!(result.world_wrench.torque_world_n_m, Vec3::ZERO);

    let mut co_moving = fixture();
    co_moving.state.angular_velocity_world_rad_per_s = Vec3::ZERO;
    co_moving.gas.velocity_world_m_per_s = co_moving.state.center_velocity_world_m_per_s;
    let result = candidate(&co_moving);
    assert_eq!(result.world_wrench.force_world_n, Vec3::ZERO);
    assert_eq!(result.world_wrench.torque_world_n_m, Vec3::ZERO);
}

#[test]
fn g0_reversal_reverses_translation_force_and_spin_torque() {
    let forward = candidate(&fixture());
    let mut reversed = fixture();
    reversed.state.center_velocity_world_m_per_s.x *= -1.0;
    reversed.state.angular_velocity_world_rad_per_s.z *= -1.0;
    let backward = candidate(&reversed);
    close(
        forward.world_wrench.force_world_n.x,
        -backward.world_wrench.force_world_n.x,
    );
    close(
        forward.world_wrench.torque_world_n_m.z,
        -backward.world_wrench.torque_world_n_m.z,
    );
}

#[test]
fn g0_correlation_specific_density_and_viscosity_behavior_is_not_overclaimed() {
    let baseline = candidate(&fixture());
    let mut doubled_density = fixture();
    doubled_density.gas.density_kg_per_m3 = Some(2.4);
    let dense = candidate(&doubled_density);
    close(
        dense.world_wrench.force_world_n.x,
        2.0 * baseline.world_wrench.force_world_n.x,
    );
    close(
        dense.world_wrench.torque_world_n_m.z,
        2.0 * baseline.world_wrench.torque_world_n_m.z,
    );

    let mut different_viscosity = fixture();
    different_viscosity.gas.dynamic_viscosity_pa_s = Some(3.6e-5);
    let viscous = candidate(&different_viscosity);
    close(
        viscous.world_wrench.force_world_n.x,
        baseline.world_wrench.force_world_n.x,
    );
    close(
        viscous.world_wrench.torque_world_n_m.z,
        baseline.world_wrench.torque_world_n_m.z,
    );

    let mut different_pressure = fixture();
    different_pressure.pressure.absolute_pressure_pa *= 2.0;
    different_pressure.pressure.source_id = "pressure:changed".into();
    let pressure_result = evaluate_euler_disc_external_air(&different_pressure)
        .expect("pressure state remains a retained input");
    assert_eq!(
        pressure_result.pressure_scaling,
        ExteriorAirPressureScaling::NoDirectScaling
    );
    let pressure_candidate = pressure_result
        .candidates
        .first()
        .expect("one retained correlation candidate");
    close(
        pressure_candidate.world_wrench.force_world_n.x,
        baseline.world_wrench.force_world_n.x,
    );
}

#[test]
fn g0_stationary_ambient_is_passive_and_work_is_exactly_once() {
    let result = candidate(&fixture());
    assert!(result.world_wrench.receipt.relative_power_w <= 0.0);
    assert!(result.world_wrench.receipt.dissipated_relative_power_w >= 0.0);
    assert!(
        !result
            .world_wrench
            .receipt
            .moving_ambient_requires_total_energy_accounting
    );
    assert_eq!(
        result.heat,
        ExteriorAirHeatDisposition::UnallocatedNoThermalModel
    );

    let mut work = EulerExternalAirWorkWindow::default();
    let receipt = work.record_once(7, 0.25, &result).expect("first exchange");
    assert!(receipt.relative_dissipation_j >= 0.0);
    assert!(matches!(
        work.record_once(7, 0.25, &result),
        Err(ExternalAirError::GenericRefusal {
            detail: ReducedAeroError::DuplicateWorkExchange { key: 7 }
        })
    ));
}

#[test]
fn g0_moving_ambient_can_do_positive_body_work_without_violating_relative_passivity() {
    let mut input = fixture();
    input.state.angular_velocity_world_rad_per_s = Vec3::ZERO;
    input.gas.velocity_world_m_per_s = Vec3::new(8.0, 0.0, 0.0);
    let result = candidate(&input);
    assert!(result.world_wrench.receipt.relative_power_w <= 0.0);
    assert!(result.world_wrench.receipt.body_power_w > 0.0);
    assert!(
        result
            .world_wrench
            .receipt
            .moving_ambient_requires_total_energy_accounting
    );
}

#[test]
fn g3_world_and_body_wrenches_are_rotation_covariant() {
    let baseline = candidate(&fixture());
    let mut rotated = fixture();
    rotated.identity.world_frame_id = "world-rotated".into();
    rotated.state.center_velocity_world_m_per_s =
        rotated_z(rotated.state.center_velocity_world_m_per_s);
    rotated.state.angular_velocity_world_rad_per_s =
        rotated_z(rotated.state.angular_velocity_world_rad_per_s);
    rotated.state.body_frame = EulerDiscBodyFrame {
        x_world: rotated_z(rotated.state.body_frame.x_world),
        z_world: rotated_z(rotated.state.body_frame.z_world),
    };
    let transformed = candidate(&rotated);
    assert_eq!(transformed.force_body_n, baseline.force_body_n);
    assert_eq!(transformed.torque_body_n_m, baseline.torque_body_n_m);
    assert_eq!(
        transformed.world_wrench.force_world_n,
        rotated_z(baseline.world_wrench.force_world_n)
    );
    assert_eq!(
        transformed.world_wrench.torque_world_n_m,
        rotated_z(baseline.world_wrench.torque_world_n_m)
    );
}

#[test]
fn g0_alternatives_remain_visible_and_forbidden_domains_refuse() {
    let mut alternatives = fixture();
    alternatives.alternatives = vec![model("corr-b", 1.4, 0.02), model("corr-a", 1.0, 0.02)];
    let result = evaluate_euler_disc_external_air(&alternatives).expect("two alternatives");
    assert_eq!(result.candidates.len(), 2);
    assert!(result.has_force_or_torque_disagreement);
    assert!(result.applicability.generic_correlation_domain_admitted);
    assert!(result.applicability.estimate_only);
    assert_eq!(result.candidates[0].world_wrench.correlation.id, "corr-a");

    alternatives.domain = ExternalAirDomain::ThinGap;
    assert_eq!(
        evaluate_euler_disc_external_air(&alternatives),
        Err(ExternalAirError::ThinGapDomainRejected)
    );
}

#[test]
fn g0_generic_correlation_domain_refusal_survives_the_adapter() {
    let constrained = ReducedAeroModel::try_new(
        CorrelationIdentity::try_new("narrow", "v1", "source:narrow").expect("identity"),
        ApplicabilityEnvelope {
            translational_reynolds: range(0.0, 1.0),
            rotational_reynolds: range(0.0, 1.0e9),
            relative_roughness: range(0.0, 1.0),
            maximum_tip_mach: 1.0,
        },
        CorrelationUncertainty {
            source_id: "uncertainty:narrow".into(),
            coefficient_relative_half_width: 0.0,
        },
        ReducedAeroComponents {
            form_drag: Some(FormDrag { coefficient: 1.0 }),
            ..ReducedAeroComponents::default()
        },
        &[ContributionFamily::TranslationalFormDrag],
    )
    .expect("model");
    let mut input = fixture();
    input.alternatives = vec![constrained];
    assert!(matches!(
        evaluate_euler_disc_external_air(&input),
        Err(ExternalAirError::GenericRefusal {
            detail: ReducedAeroError::OutsideCorrelationDomain {
                quantity: "translational_reynolds",
                ..
            }
        })
    ));
}

#[test]
fn g0_generic_model_refuses_thin_gap_and_target_fitted_families_before_adapter_use() {
    let identity = CorrelationIdentity::try_new("bad", "v1", "source:bad").expect("identity");
    let envelope = ApplicabilityEnvelope {
        translational_reynolds: range(0.0, 1.0),
        rotational_reynolds: range(0.0, 1.0),
        relative_roughness: range(0.0, 1.0),
        maximum_tip_mach: 1.0,
    };
    let uncertainty = CorrelationUncertainty {
        source_id: "uncertainty:bad".into(),
        coefficient_relative_half_width: 0.0,
    };
    let components = ReducedAeroComponents {
        form_drag: Some(FormDrag { coefficient: 1.0 }),
        ..ReducedAeroComponents::default()
    };
    for forbidden in [
        ContributionFamily::ThinGapPressure,
        ContributionFamily::TargetFitted,
    ] {
        assert!(matches!(
            ReducedAeroModel::try_new(
                identity.clone(),
                envelope,
                uncertainty.clone(),
                components,
                &[ContributionFamily::TranslationalFormDrag, forbidden],
            ),
            Err(ReducedAeroError::ForbiddenContribution { family }) if family == forbidden
        ));
    }
}

#[test]
fn g0_adapter_maps_only_admitted_exterior_disc_geometry() {
    let result = candidate(&fixture());
    assert_eq!(
        result.world_wrench.receipt.authority,
        fs_flux::EstimateAuthority::EstimateOnly
    );
    let _generic_geometry = DiscGeometry {
        radius_m: 0.04,
        exterior_thickness_m: 0.003,
    };
}
