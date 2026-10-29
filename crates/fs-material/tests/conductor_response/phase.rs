//! G0/G3: one original source with independently owned phase and actual E/J.
use super::*;
use fs_material::conductor::planar::PlanarConductorMaterial;
use fs_qty::Dims;

const PHASE: &str = "liquid_mass_fraction";
const AXES: [(&str, ConductorMagnitude); 2] = [
    ("field", ConductorMagnitude::ElectricField),
    ("current", ConductorMagnitude::CurrentDensity),
];

fn phase_source(frequency_bound: bool) -> (MaterialCard, ResolvedConductorState) {
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = PropertyKey::new("declared_phase_transport", basis.dimensions());
    let scalar = QuantitySpec::dimensional(Dims([0; 6]));
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (
            key.clone(),
            PropertyValue::Surface2D {
                axes: [PHASE.into(), "field".into()],
                axis_dims: [scalar.dims(), AXES[0].1.quantity().dims()],
                samples: [[0., 0.], [1., 0.], [0., 10.], [1., 10.]]
                    .into_iter()
                    .map(|p| (p, 1. + 3. * p[0] + 0.1 * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
        ),
        (
            PropertyKey::new("density", Density::DIMS),
            PropertyValue::Scalar {
                value: 1.,
                dims: Density::DIMS,
            },
        ),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300., 350.);
        if key.name() == "declared_phase_transport" {
            validity = validity
                .with_quantity(PHASE, scalar, 0., 1.)
                .with_quantity("field", AXES[0].1.quantity(), 0., 10.)
                .with_quantity("current", AXES[1].1.quantity(), 0., 20.);
            if frequency_bound {
                validity = validity.with_quantity(
                    "frequency",
                    QuantitySpec::dimensional(Frequency::DIMS),
                    0.,
                    10.,
                );
            }
        }
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured phase and field conductor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid-liquid chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut point = QueryPoint::new()
        .with("T", 300.)
        .unwrap()
        .with_quantity(PHASE, scalar, 0.)
        .unwrap()
        .with_quantity("field", AXES[0].1.quantity(), 0.)
        .unwrap()
        .with_quantity("current", AXES[1].1.quantity(), 0.)
        .unwrap();
    if frequency_bound {
        point = point
            .with_quantity("frequency", QuantitySpec::dimensional(Frequency::DIMS), 0.)
            .unwrap();
    }
    let material = resolve_material_state_point(
        &card,
        &point,
        &[
            basis.requirement_for(&key).unwrap(),
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
        ResolvedConductorState::try_from_material_with_property(material, basis, &key).unwrap();
    (card, original)
}

#[test]
fn one_pinned_phase_field_source_refreshes_without_weakening_static_ingress() {
    let (card, original) = phase_source(false);
    assert!(
        original
            .clone()
            .with_liquid_mass_fraction_axis(PHASE)
            .is_err()
    );
    assert!(ScalarConductorResponse::new(&card, original.clone(), &AXES).is_err());
    let source = ScalarConductorResponse::new_with_liquid_mass_fraction_axis(
        &card,
        original.clone(),
        &AXES,
        PHASE,
    )
    .unwrap();
    assert_eq!(source.original().liquid_mass_fraction_axis(), Some(PHASE));
    assert!(!source.original().has_admitted_static_sources());
    assert!(
        source
            .original()
            .with_material(original.material().clone())
            .is_err()
    );
    let query = original
        .material()
        .query()
        .with("T", 310.)
        .unwrap()
        .with(PHASE, 0.4)
        .unwrap();
    let seed = source.at_query_point(&query).unwrap();
    let next = source.at_fields(&seed, [3., 4.], [6., 8.]).unwrap();
    assert!((next.resistivity_ohm_m() - 2.7).abs() < 1e-13);
    assert_eq!(next.material().query().axes()["field"], 5.);
    assert_eq!(next.material().query().axes()["current"], 10.);
    assert_eq!(next.material().query().axes()["T"], 310.);
    assert_eq!(next.material().query().axes()[PHASE], 0.4);
    assert_eq!(next.source_key(), original.source_key());
    assert_eq!(next.source_basis(), original.source_basis());
    for (before, after) in original
        .material()
        .properties()
        .iter()
        .zip(next.material().properties())
    {
        assert_eq!(before.requirement(), after.requirement());
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
    }
    let liquid = source
        .at_query_point(
            &next
                .material()
                .query()
                .with("T", 320.)
                .unwrap()
                .with(PHASE, 1.)
                .unwrap(),
        )
        .unwrap();
    assert!((liquid.resistivity_ohm_m() - 4.5).abs() < 1e-13);
    assert_eq!(liquid.liquid_mass_fraction_axis(), Some(PHASE));
    assert_eq!(original.material().query().axes()[PHASE], 0.);
}

#[test]
fn phase_field_binding_refuses_unowned_retyped_and_unsupported_coordinates() {
    let (card, original) = phase_source(false);
    for axis in ["", "missing", "field", "T"] {
        assert!(
            ScalarConductorResponse::new_with_liquid_mass_fraction_axis(
                &card,
                original.clone(),
                &AXES,
                axis
            )
            .is_err()
        );
    }
    assert!(
        ScalarConductorResponse::new_with_liquid_mass_fraction_axis(
            &card,
            original.clone(),
            &AXES[..1],
            PHASE
        )
        .is_err()
    );
    assert!(
        ScalarConductorResponse::new_with_liquid_mass_fraction_axis(
            &card,
            original.clone(),
            &[(PHASE, ConductorMagnitude::ElectricField)],
            PHASE
        )
        .is_err()
    );
    let source = ScalarConductorResponse::new_with_liquid_mass_fraction_axis(
        &card,
        original.clone(),
        &AXES,
        PHASE,
    )
    .unwrap();
    assert!(source.admit(&original).is_err());
    for (axis, value) in [
        (PHASE, 1.01),
        (PHASE, -0.1),
        ("T", 351.),
        ("field", 11.),
        ("current", 21.),
    ] {
        assert!(
            source
                .at_query_point(&original.material().query().with(axis, value).unwrap())
                .is_err()
        );
    }
    assert!(
        source
            .at_query_point(
                &original
                    .material()
                    .query()
                    .with_quantity(PHASE, QuantitySpec::dimensional(Frequency::DIMS), 0.5)
                    .unwrap()
            )
            .is_err()
    );
    let (card, original) = phase_source(true);
    assert!(
        ScalarConductorResponse::new_with_liquid_mass_fraction_axis(&card, original, &AXES, PHASE)
            .is_err()
    );
}
