use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::Pressure;

#[test]
fn original_ordered_source_survives_damage_and_frame_transport_without_reselection() {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured objective cohesive fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured fixture".into(),
            revision: 1,
        },
        texture_frame: name.into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit bond".into(),
            third_body: None,
            environment: "fixed source query".into(),
            history: "original bond".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let source = resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap();
    let damaged = source
        .law()
        .bonded_point()
        .trial([3e-4, 2e-4, -1e-4], 1e-8)
        .unwrap()
        .point;
    let numerical = ObjectiveCohesivePoint::new(
        damaged,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.003, -0.001, 0.002),
        0.002,
    )
    .unwrap();
    let receiving = [&card.surface_a().material, &card.surface_b().material];
    assert!(
        SourceBoundObjectiveCohesivePoint::new(numerical, &source, [receiving[1], receiving[0]])
            .is_err()
    );
    // Numerically equal data from a distinct original receipt must not silently
    // take custody of an existing point's history.
    let other_source = source.clone();
    assert!(matches!(
        SourceBoundObjectiveCohesivePoint::new(numerical, &other_source, receiving),
        Err(SourceBoundObjectiveCohesiveError::DifferentOriginalLaw)
    ));
    let bound = SourceBoundObjectiveCohesivePoint::new(numerical, &source, receiving).unwrap();
    let receipts = source.resolved().properties().as_ptr();
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 2.0).unwrap();
    let rotated = bound.transport_frame(frame).unwrap();
    assert_eq!(rotated.numerical().material(), damaged);
    let config = ObjectiveCohesiveConfig {
        maximum_cayley_increment: 8.0,
        kinematic_tolerance_m: 2e-14,
        material_energy_tolerance_j_m2: 1e-8,
        spatial_work_tolerance_j: 1e-10,
    };
    let accepted = rotated
        .trial_local_jump(
            frame.right_exp(Vec3::new(0.3, -0.2, 0.4)).unwrap(),
            [6e-4, 3e-4, -2e-4],
            0.02,
            config,
        )
        .unwrap()
        .point();
    assert!(accepted.numerical().material().damage() > damaged.damage());
    assert!(core::ptr::eq(accepted.source(), &source));
    assert!(core::ptr::eq(
        accepted.numerical().material().law(),
        source.law()
    ));
    assert_eq!(accepted.source().resolved().properties().as_ptr(), receipts);
}
