//! Original source roles, units and state context survive spatial E/J refresh.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ElectricalPropertyBasis, ResolvedConductorState,
        response::{ConductorMagnitude, ScalarConductorResponse},
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::{Density, Frequency, QuantitySpec};

fn fixture(frequency_bound: bool) -> (MaterialCard, ResolvedConductorState) {
    let e = ConductorMagnitude::ElectricField;
    let j = ConductorMagnitude::CurrentDensity;
    let basis = ElectricalPropertyBasis::Conductivity;
    let electrical = PropertyKey::new("declared-transport", basis.dimensions());
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (
            electrical.clone(),
            PropertyValue::Surface2D {
                axes: ["field".into(), "current".into()],
                axis_dims: [e.quantity().dims(), j.quantity().dims()],
                samples: [[0.0, 0.0], [10.0, 0.0], [0.0, 20.0], [10.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, 2.0 + 0.5 * p[0] + 0.1 * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
        ),
        (
            PropertyKey::new("density", Density::DIMS),
            PropertyValue::Scalar {
                value: 9000.0,
                dims: Density::DIMS,
            },
        ),
    ] {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        let mut validity = ValidityDomain::unconstrained().with("T", 250.0, 500.0);
        if key == electrical {
            validity = validity
                .with_quantity("field", e.quantity(), 0.0, 10.0)
                .with_quantity("current", j.quantity(), 0.0, 20.0);
            if frequency_bound {
                validity = validity.with_quantity(
                    "omega",
                    QuantitySpec::dimensional(Frequency::DIMS),
                    0.0,
                    10.0,
                );
            }
        }
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured complete conductor source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut query = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity("field", e.quantity(), 0.0)
        .unwrap()
        .with_quantity("current", j.quantity(), 0.0)
        .unwrap();
    if frequency_bound {
        query = query
            .with_quantity("omega", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
            .unwrap();
    }
    let material = resolve_material_state_point(
        &card,
        &query,
        &[
            basis.requirement_for(&electrical).unwrap(),
            ScalarPropertyRequirement::try_new(
                "density",
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original =
        ResolvedConductorState::try_from_material_with_property(material, basis, &electrical)
            .unwrap();
    (card, original)
}

#[test]
fn actual_rotational_invariants_refresh_one_complete_pinned_material_bundle() {
    let (card, original) = fixture(false);
    let axes = [
        ("field", ConductorMagnitude::ElectricField),
        ("current", ConductorMagnitude::CurrentDensity),
    ];
    let law = ScalarConductorResponse::new(&card, original.clone(), &axes).unwrap();
    let seed = law
        .at_query_point(&original.material().query().with("T", 350.0).unwrap())
        .unwrap();
    let next = law.at_fields(&seed, [3.0, 4.0], [6.0, 8.0]).unwrap();
    assert!((next.resistivity_ohm_m() - 1.0 / 5.5).abs() < 1e-15);
    assert_eq!(next.source_key(), original.source_key());
    assert_eq!(next.material().query().axes()["field"], 5.0);
    assert_eq!(next.material().query().axes()["current"], 10.0);
    assert_eq!(next.material().query().axes()["T"], 350.0);
    assert_eq!(
        next.material().property("density").unwrap().value_si(),
        9000.0
    );
    for property in next.material().properties() {
        let before = original
            .material()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(property.requirement(), before.requirement());
        assert_eq!(
            property.answer().receipt.selected,
            before.answer().receipt.selected
        );
    }
    let reversed = law.at_fields(&seed, [-4.0, 3.0], [-8.0, 6.0]).unwrap();
    assert_eq!(next, reversed);
    assert_eq!(law.original(), &original);
}

#[test]
fn multiple_owned_axes_do_not_admit_unowned_or_retyped_conditions() {
    let (card, original) = fixture(false);
    let e = ConductorMagnitude::ElectricField;
    let j = ConductorMagnitude::CurrentDensity;
    let property = original.source_property();
    assert!(!property.has_only_owned_and_fixed_external_conditions("field", e.quantity()));
    assert!(
        property.has_only_owned_coordinates_and_fixed_external_conditions(&[
            ("field", e.quantity()),
            ("current", j.quantity())
        ])
    );
    for owned in [
        vec![("field", e.quantity()), ("field", e.quantity())],
        vec![("field", j.quantity()), ("current", j.quantity())],
        vec![("field", e.quantity()), ("missing", j.quantity())],
    ] {
        assert!(!property.has_only_owned_coordinates_and_fixed_external_conditions(&owned));
    }
    assert!(ScalarConductorResponse::new(&card, original.clone(), &[("field", e)]).is_err());
    assert!(
        ScalarConductorResponse::new(&card, original.clone(), &[("field", e), ("current", e)])
            .is_err()
    );
    let law =
        ScalarConductorResponse::new(&card, original.clone(), &[("field", e), ("current", j)])
            .unwrap();
    assert!(law.at_fields(&original, [11.0, 0.0], [1.0, 0.0]).is_err());
    assert!(
        law.at_fields(&original, [f64::NAN, 0.0], [0.0, 0.0])
            .is_err()
    );
    assert!(
        law.at_query_point(&original.material().query().with("T", 600.0).unwrap())
            .is_err()
    );
    let (card, original) = fixture(true);
    assert!(
        ScalarConductorResponse::new(&card, original, &[("field", e), ("current", j)]).is_err()
    );
}

#[test]
fn refreshing_context_cannot_drop_other_original_material_roles() {
    let (card, original) = fixture(false);
    let law = ScalarConductorResponse::new(
        &card,
        original.clone(),
        &[
            ("field", ConductorMagnitude::ElectricField),
            ("current", ConductorMagnitude::CurrentDensity),
        ],
    )
    .unwrap();
    let incomplete = resolve_material_state_point(
        &card,
        &original.material().query(),
        &[original
            .source_basis()
            .requirement_for(original.source_key())
            .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let incomplete = ResolvedConductorState::try_from_material_with_property(
        incomplete,
        original.source_basis(),
        original.source_key(),
    )
    .unwrap();
    assert!(law.admit(&incomplete).is_err());
    assert!(law.at_fields(&incomplete, [1.0, 0.0], [1.0, 0.0]).is_err());
}
