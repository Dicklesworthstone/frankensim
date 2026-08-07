//! G0/G3 coverage for the Euler finite-patch normal-contact adapter.

#[path = "../src/contact_dynamics.rs"]
mod contact_dynamics;
#[path = "../src/normal_contact.rs"]
mod normal_contact;
#[path = "../src/patch_kinematics.rs"]
mod patch_kinematics;

use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, NormalPatchEmbedState,
    NormalPatchPort, NormalPatchReceipt,
};
use fs_couple::StableId;
use fs_mbd::{PointKinematics, Vec3};
use fs_rep_frep::AxisymmetricSupportAuthority;
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};
use normal_contact::{
    EulerNormalContactInput, EulerNormalContactOutcome, EulerNormalGeometry,
    NORMAL_CONTACT_ADAPTER_ID, NormalContactError, NormalContactIdentity, NormalDissipation,
    NormalElasticStorage, NormalMaterialInterface, evaluate_normal_contact,
};
use patch_kinematics::{
    Creepage, CurvatureMetadata, OrderedSurfacePair, PatchContactStatus, PatchGeometryMetadata,
    PatchKinematics, SurfaceOrder, TangentBasis, TangentComponents, TangentGaugeSource,
};

fn id(value: &str) -> StableId {
    StableId::new(value).expect("test identity")
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-12 * actual.abs().max(expected.abs()).max(1.0),
        "{actual} != {expected}"
    );
}

fn point(point_world: Vec3) -> PointKinematics {
    PointKinematics {
        arm_body: Vec3::new(1.0, 0.0, 0.0),
        arm_world: Vec3::new(1.0, 0.0, 0.0),
        point_world,
        center_of_mass_velocity_world: Vec3::ZERO,
        angular_velocity_body: Vec3::ZERO,
        angular_velocity_world: Vec3::ZERO,
        point_velocity_world: Vec3::ZERO,
    }
}

fn kinematics(
    status: PatchContactStatus,
    first_curvature_m_inverse: f64,
    second_curvature_m_inverse: f64,
    gap_m: f64,
    normal_velocity_m_per_s: f64,
) -> PatchKinematics {
    PatchKinematics {
        surfaces: OrderedSurfacePair::try_new(id("disc"), id("base"), SurfaceOrder::DiscThenBase)
            .expect("ordered surfaces"),
        patch: PatchGeometryMetadata {
            patch_identity: id("patch/rim-v1"),
            source_feature: 3,
            gap_uncertainty_m: 1.0e-7,
            curvature: CurvatureMetadata::Known {
                curvature_identity: id("curvature/rim-v1"),
                authority: InputAuthority::SyntheticFixture,
                first_principal_m_inverse: first_curvature_m_inverse,
                second_principal_m_inverse: second_curvature_m_inverse,
                uncertainty_m_inverse: 0.01,
            },
        },
        support_authority: AxisymmetricSupportAuthority::Estimate,
        disc_point: point(Vec3::new(1.0, 0.0, gap_m)),
        base_point: point(Vec3::ZERO),
        tangent_basis: TangentBasis {
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            first_world: Vec3::new(1.0, 0.0, 0.0),
            second_world: Vec3::new(0.0, 1.0, 0.0),
            source: TangentGaugeSource::CallerReference,
        },
        relative_velocity_world_m_per_s: Vec3::new(0.0, 0.0, normal_velocity_m_per_s),
        normal_relative_velocity_m_per_s: normal_velocity_m_per_s,
        tangential_relative_velocity_world_m_per_s: Vec3::ZERO,
        tangential_relative_velocity: TangentComponents {
            first: 0.0,
            second: 0.0,
        },
        rolling_entrainment_velocity_world_m_per_s: Vec3::ZERO,
        rolling_entrainment_tangent_world_m_per_s: Vec3::ZERO,
        reference_rolling_speed_m_per_s: 0.0,
        normal_spin_rad_per_s: 0.0,
        creepage: Creepage::Unavailable {
            reference_rolling_speed_m_per_s: 0.0,
            minimum_reference_rolling_speed_m_per_s: 1.0e-6,
        },
        tangential_power_w: None,
        status,
    }
}

fn material(dissipation: Option<f64>) -> NormalMaterialInterface {
    NormalMaterialInterface {
        material_card_id: "card/steel-rim-v1".into(),
        model_id: "law/hertz-half-space-v1".into(),
        source_id: "source/synthetic-coupon-v1".into(),
        interface: InterfaceSystemRef::new(
            "disc/rim->base/track",
            "history/normal-v1",
            "source/interface-v1",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .expect("interface"),
        reduced_modulus_pa: 2.0e9,
        hunt_crossley_dissipation_s_per_m: dissipation,
        applicability: ApplicabilityInput {
            half_space_depth_m: 1.0,
            layer_thickness_m: 1.0,
            yield_strength_pa: 1.0e12,
            characteristic_rate_m_per_s: 1.0,
            temperature_k: 293.15,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 0.2,
            max_strain: 0.1,
            max_patch_to_depth: 0.2,
            max_patch_to_layer: 0.2,
            max_pressure_to_yield: 0.2,
            max_rate_ratio: 1.0,
            min_temperature_k: 200.0,
            max_temperature_k: 400.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.001,
            modulus_relative: 0.02,
            load_relative: 0.03,
        },
    }
}

fn input(
    status: PatchContactStatus,
    curvature: (f64, f64),
    gap_m: f64,
    normal_velocity_m_per_s: f64,
    geometry: EulerNormalGeometry,
    dissipation: Option<f64>,
) -> EulerNormalContactInput {
    EulerNormalContactInput {
        identity: NormalContactIdentity {
            case_id: "case/normal-v1".into(),
            adapter_id: NORMAL_CONTACT_ADAPTER_ID.into(),
            solver_id: "solver/fixed-v1".into(),
            contact_id: "contact/rim-v1".into(),
            sample_id: "sample/1".into(),
        },
        kinematics: kinematics(
            status,
            curvature.0,
            curvature.1,
            gap_m,
            normal_velocity_m_per_s,
        ),
        material: material(dissipation),
        geometry,
        state: NormalPatchEmbedState::new(0.0, 1.0).expect("state"),
        time_s: 0.1,
        iteration: 1,
        step_s: 0.01,
        converged: true,
    }
}

#[test]
fn sphere_closing_maps_hertz_scale_force_and_application_point() {
    let result = evaluate_normal_contact(&input(
        PatchContactStatus::Approaching,
        (10.0, 10.0),
        -1.0e-4,
        -0.2,
        EulerNormalGeometry::SpherePlane,
        None,
    ))
    .expect("admitted sphere contact");
    let EulerNormalContactOutcome::Active(active) = result else {
        panic!("active contact expected");
    };
    close(active.curvature.reporting_radius_m, 0.1);
    assert_eq!(active.curvature.authority, InputAuthority::SyntheticFixture);
    close(active.application_point_world_m.x, 1.0);
    close(active.application_point_world_m.z, -1.0e-4);
    let expected_force = (4.0 / 3.0) * 2.0e9 * 0.1_f64.sqrt() * (1.0e-4_f64).powf(1.5);
    match &active.generic.receipt {
        NormalPatchReceipt::Point(receipt) => {
            close(receipt.normal_force_n, expected_force);
            close(receipt.approach_m, 1.0e-4);
        }
        NormalPatchReceipt::Line(_) => panic!("point receipt expected"),
    }
    match &active.generic.port {
        NormalPatchPort::Point(port) => {
            close(port.action_force_n[2], expected_force);
            close(port.action_moment_n_m[1], -expected_force);
            close(port.residual_force_n[0], 0.0);
        }
        NormalPatchPort::Line(_) => panic!("point port expected"),
    }
    assert!(
        matches!(active.elastic_storage, NormalElasticStorage::PointJoules(value) if value > 0.0)
    );
}

#[test]
fn separated_contact_is_inactive_without_consuming_work_state() {
    let mut input = input(
        PatchContactStatus::Separated,
        (10.0, 10.0),
        0.02,
        0.3,
        EulerNormalGeometry::SpherePlane,
        None,
    );
    // An open contact needs neither a material law nor an integration sample.
    input.material.material_card_id.clear();
    input.material.reduced_modulus_pa = f64::NAN;
    input.time_s = f64::NAN;
    input.step_s = 0.0;
    input.converged = false;
    let result = evaluate_normal_contact(&input).expect("separation is not a law call");
    match result {
        EulerNormalContactOutcome::InactiveSeparated {
            gap_m,
            normal_relative_velocity_m_per_s,
            state,
        } => {
            close(gap_m, 0.02);
            close(normal_relative_velocity_m_per_s, 0.3);
            assert_eq!(state, input.state);
        }
        EulerNormalContactOutcome::Active(_) => panic!("separated contact must stay inactive"),
    }
}

#[test]
fn initial_time_zero_is_a_valid_active_sample() {
    let mut initial = input(
        PatchContactStatus::Touching,
        (10.0, 10.0),
        -1.0e-4,
        0.0,
        EulerNormalGeometry::SpherePlane,
        None,
    );
    initial.time_s = 0.0;
    initial.iteration = 1;
    assert!(matches!(
        evaluate_normal_contact(&initial),
        Ok(EulerNormalContactOutcome::Active(_))
    ));
}

#[test]
fn unloading_hunt_crossley_keeps_dissipation_nonnegative() {
    let closing = evaluate_normal_contact(&input(
        PatchContactStatus::Approaching,
        (10.0, 10.0),
        -1.0e-4,
        -0.2,
        EulerNormalGeometry::SpherePlane,
        Some(1.0),
    ))
    .expect("closing step");
    let opening = evaluate_normal_contact(&input(
        PatchContactStatus::Grazing,
        (10.0, 10.0),
        -1.0e-4,
        0.1,
        EulerNormalGeometry::SpherePlane,
        Some(1.0),
    ))
    .expect("unloading step");
    let EulerNormalContactOutcome::Active(closing) = closing else {
        panic!("active closing")
    };
    let EulerNormalContactOutcome::Active(opening) = opening else {
        panic!("active opening")
    };
    let NormalPatchReceipt::Point(closing_receipt) = closing.generic.receipt else {
        panic!("point")
    };
    let NormalPatchReceipt::Point(opening_receipt) = opening.generic.receipt else {
        panic!("point")
    };
    assert!(opening_receipt.normal_force_n < closing_receipt.normal_force_n);
    assert!(
        matches!(opening.dissipation, NormalDissipation::Point { work_j, power_w } if work_j >= 0.0 && power_w >= 0.0)
    );
}

#[test]
fn line_contact_preserves_line_units_and_elastic_storage() {
    let result = evaluate_normal_contact(&input(
        PatchContactStatus::Touching,
        (10.0, 0.0),
        0.0,
        0.0,
        EulerNormalGeometry::CylinderPlane {
            line_load_n_per_m: 100.0,
        },
        None,
    ))
    .expect("admitted line contact");
    let EulerNormalContactOutcome::Active(active) = result else {
        panic!("active line")
    };
    assert!(matches!(
        &active.generic.receipt,
        NormalPatchReceipt::Line(_)
    ));
    assert!(matches!(&active.generic.port, NormalPatchPort::Line(_)));
    assert!(
        matches!(active.elastic_storage, NormalElasticStorage::LineJoulesPerMetre(value) if value > 0.0)
    );
    assert!(matches!(
        active.dissipation,
        NormalDissipation::Line {
            work_j_per_m: 0.0,
            power_w_per_m: 0.0
        }
    ));
}

#[test]
fn unequal_positive_curvatures_map_to_a_true_elliptic_patch() {
    let result = evaluate_normal_contact(&input(
        PatchContactStatus::Approaching,
        (20.0, 10.0),
        -1.0e-4,
        -0.1,
        EulerNormalGeometry::EllipticParaboloid,
        None,
    ))
    .expect("admitted elliptic contact");
    let EulerNormalContactOutcome::Active(active) = result else {
        panic!("active elliptic contact expected");
    };
    match &active.generic.receipt {
        NormalPatchReceipt::Point(receipt) => {
            let axes = receipt
                .elliptic_patch_axes
                .expect("elliptic semiaxes retained");
            assert!(axes.semi_major_axis_m > axes.semi_minor_axis_m);
            assert!(receipt.normal_force_n > 0.0);
        }
        NormalPatchReceipt::Line(_) => panic!("point receipt expected"),
    }
    close(active.curvature.reporting_radius_m, 1.0 / 200.0_f64.sqrt());

    let dissipative = evaluate_normal_contact(&input(
        PatchContactStatus::Approaching,
        (20.0, 10.0),
        -1.0e-4,
        -0.1,
        EulerNormalGeometry::EllipticParaboloid,
        Some(0.1),
    ));
    assert!(matches!(
        dissipative,
        Err(NormalContactError::DissipativeEllipticUnsupported)
    ));
}

#[test]
fn hostile_curvature_unknown_and_event_states_refuse() {
    let toroidal = evaluate_normal_contact(&input(
        PatchContactStatus::Touching,
        (10.0, 20.0),
        0.0,
        0.0,
        EulerNormalGeometry::SpherePlane,
        None,
    ));
    assert!(matches!(
        toroidal,
        Err(NormalContactError::SphereCurvatureMismatch { .. })
    ));
    let unknown = evaluate_normal_contact(&input(
        PatchContactStatus::Unknown,
        (10.0, 10.0),
        0.0,
        0.0,
        EulerNormalGeometry::SpherePlane,
        None,
    ));
    assert!(matches!(
        unknown,
        Err(NormalContactError::UnavailableKinematics {
            status: PatchContactStatus::Unknown
        })
    ));
    let event = evaluate_normal_contact(&input(
        PatchContactStatus::ImpactCandidate,
        (10.0, 10.0),
        0.0,
        -2.0,
        EulerNormalGeometry::SpherePlane,
        None,
    ));
    assert!(matches!(
        event,
        Err(NormalContactError::UnavailableKinematics {
            status: PatchContactStatus::ImpactCandidate
        })
    ));
}

#[test]
fn deterministic_receipt_replays_with_identical_source_identity() {
    let value = input(
        PatchContactStatus::Touching,
        (10.0, 10.0),
        -1.0e-4,
        0.0,
        EulerNormalGeometry::SpherePlane,
        None,
    );
    let left = evaluate_normal_contact(&value).expect("left");
    let right = evaluate_normal_contact(&value).expect("right");
    let (EulerNormalContactOutcome::Active(left), EulerNormalContactOutcome::Active(right)) =
        (left, right)
    else {
        panic!("active receipts")
    };
    assert_eq!(left.curvature.curvature_identity, "curvature/rim-v1");
    assert_eq!(left.generic.receipt_id, right.generic.receipt_id);
    assert_eq!(left.generic.embedding_id, right.generic.embedding_id);
    assert_eq!(left.generic.receipt, right.generic.receipt);
}
