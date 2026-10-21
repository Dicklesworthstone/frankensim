//! Synthetic source support and complete-bundle propagation, not measured data.
#[path = "material_surface/path_physics.rs"]
mod path_physics;

#[path = "material_surface/caloric.rs"]
mod caloric;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};

#[test]
fn g0_surface_refresh_keeps_both_axes_and_does_not_freeze_at_a_vertex() {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    let selected = claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("young_modulus", Pressure::DIMS),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "pressure".into()],
                axis_dims: [Temperature::DIMS, Pressure::DIMS],
                samples: vec![
                    ([300.0, 100.0], 2e9),
                    ([400.0, 100.0], 3e9),
                    ([300.0, 200.0], 4e9),
                ],
                triangles: vec![[0, 1, 2]],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", absolute, 300.0, 400.0)
                .with_quantity(
                    "pressure",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    100.0,
                    200.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured two-variable elasticity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic two-axis source".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = |t, p| {
        QueryPoint::new()
            .with_quantity("T", absolute, t)
            .unwrap()
            .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), p)
            .unwrap()
    };
    let requirement = ScalarPropertyRequirement::try_new(
        "young_modulus",
        Pressure::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap();
    let vertex = resolve_material_state_point(
        &card,
        &query(300.0, 100.0),
        &[requirement],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let interior = vertex.at_query_point(&card, &query(325.0, 125.0)).unwrap();
    assert_eq!(
        interior.property("young_modulus").unwrap().value_si(),
        2.75e9
    );
    for state in [&vertex, &interior] {
        let property = state.property("young_modulus").unwrap();
        assert_eq!(
            property.source_axes(),
            &[
                ("T".into(), Temperature::DIMS),
                ("pressure".into(), Pressure::DIMS),
            ]
        );
        assert!(!property.is_constant_at_fixed_temperature());
        assert_eq!(property.answer().receipt.selected, selected);
        assert_eq!(state.card_identity(), card.content_hash());
    }
    // A rectangle around sparse observations is not source coverage.
    assert!(vertex.at_query_point(&card, &query(375.0, 175.0)).is_err());
    assert_eq!(vertex.query(), query(300.0, 100.0));
}
