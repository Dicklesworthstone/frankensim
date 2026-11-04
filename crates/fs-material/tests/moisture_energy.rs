//! Authored sorption caloric paths; no measured wood or sorbent validation.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};

const NAME: &str = MOISTURE_ISOSTERIC_HEAT_PROPERTY;
const ENERGY: Dims = MOISTURE_ISOSTERIC_HEAT_DIMS;
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(m: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "moisture_content",
            quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
            m,
        )
        .unwrap()
}
fn claim(value: PropertyValue) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(NAME, ENERGY),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                290.0,
                310.0,
            )
            .with_quantity(
                "moisture_content",
                quantity(QuantityKind::Moisture(MoistureBasis::DryMass)),
                0.0,
                1.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "authored differential isosteric heat, J/kg water".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn assemble(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored bound-water energy specimen".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn card(value: PropertyValue) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim(value)).unwrap();
    assemble(claims)
}
fn curve(knots: Vec<(f64, f64)>) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "moisture_content".into(),
        abscissa_dims: Dims::NONE,
        knots,
        dims: ENERGY,
    }
}
fn binding(c: &MaterialCard, p: &QueryPoint) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        c,
        p,
        &[ScalarPropertyRequirement::try_new(NAME, ENERGY, ScalarAdmissibility::Finite).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn chart(c: &MaterialCard, bounds: [f64; 2]) -> MoistureSorptionEnergy {
    MoistureSorptionEnergy::from_binding(c, &binding(c, &point(bounds[0])), NAME, bounds).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn all_loading_segments_contribute_and_small_increments_preserve_their_energy() {
    let c = card(curve(vec![(0.0, 3e6), (0.1, 4e6), (0.3, 2e6), (0.5, 3e6)]));
    let law = chart(&c, [0.0, 0.5]);
    // Independent trapezoids: .1*3.5e6 + .2*3e6 + .2*2.5e6.
    near(law.integral(0.0, 0.5).unwrap(), 1.45e6, 1e-9);
    assert_eq!(
        law.integral(0.5, 0.0).unwrap(),
        -law.integral(0.0, 0.5).unwrap()
    );
    near(
        law.integral(0.0, 0.17).unwrap() + law.integral(0.17, 0.5).unwrap(),
        1.45e6,
        1e-9,
    );
    assert_eq!(law.integration_segments(0.0, 0.5).unwrap(), 3);
    assert_eq!(law.integration_segments(0.1, 0.3).unwrap(), 1);
    assert_eq!(law.integration_segments(0.3, 0.3).unwrap(), 0);
    let m = 0.45_f64;
    let next = f64::from_bits(m.to_bits() + 1);
    near(law.integral(m, next).unwrap() / (next - m), 2.75e6, 2e-9);
    assert!(law.integral(m, next).unwrap() > 0.0);
    assert_eq!(law, chart(&c, [0.0, 0.5]));
}

#[test]
fn signed_heats_and_point_support_do_not_invent_an_exothermic_default() {
    let c = card(curve(vec![(0.0, -2e6), (1.0, 2e6)]));
    let law = chart(&c, [0.0, 1.0]);
    assert_eq!(law.integral(0.0, 1.0).unwrap(), 0.0);
    near(law.integral(0.0, 0.5).unwrap(), -5e5, 1e-10);
    near(law.integral(0.5, 1.0).unwrap(), 5e5, 1e-10);
    let single = chart(&c, [0.25, 0.25]);
    assert_eq!(single.knot_count(), 1);
    assert_eq!(single.isosteric_heat(0.25).unwrap(), -1e6);
    assert_eq!(single.integral(0.25, 0.25).unwrap(), 0.0);
    assert!(single.integral(0.25, 0.250001).is_err());
    for invalid in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(law.integral(0.5, invalid).is_err());
        assert!(law.integration_segments(invalid, invalid).is_err());
    }
}

#[test]
fn prepared_path_preserves_pinned_claim_fixed_context_and_original_card() {
    let mut claims = ClaimSet::new();
    let selected = claims
        .insert_claim(claim(curve(vec![(0.0, 2e6), (1.0, 3e6)])))
        .unwrap();
    let alternative = claims
        .insert_claim(claim(curve(vec![(0.0, 9e6), (1.0, 9e6)])))
        .unwrap();
    assert_ne!(selected, alternative);
    let c = assemble(claims);
    let source = resolve_material_state_point(
        &c,
        &point(0.1),
        &[ScalarPropertyRequirement::try_new(NAME, ENERGY, ScalarAdmissibility::Finite).unwrap()],
        MaterialPropertySelection::PinnedByProperty([(NAME.to_owned(), selected)].into()),
    )
    .unwrap();
    let law = MoistureSorptionEnergy::from_binding(&c, &source, NAME, [0.0, 1.0]).unwrap();
    assert_eq!(law.selected_claim(), selected);
    assert_eq!(law.card_identity(), c.content_hash());
    assert_eq!(law.source_point(), &point(0.1));
    near(law.integral(0.0, 1.0).unwrap(), 2.5e6, 1e-9);
    let other = card(curve(vec![(0.0, 2e6), (1.0, 3e6)]));
    assert!(MoistureSorptionEnergy::from_binding(&other, &source, NAME, [0.0, 1.0]).is_err());
    let wet_point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "moisture_content",
            quantity(QuantityKind::Moisture(MoistureBasis::WetMass)),
            0.1,
        )
        .unwrap();
    assert!(
        resolve_material_state_point(
            &other,
            &wet_point,
            &[
                ScalarPropertyRequirement::try_new(NAME, ENERGY, ScalarAdmissibility::Finite)
                    .unwrap()
            ],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn held_temperature_slice_integrates_its_actual_cells_and_refuses_a_source_hole() {
    let c = card(PropertyValue::Surface2D {
        axes: ["moisture_content".into(), "temperature".into()],
        axis_dims: [Dims::NONE, Temperature::DIMS],
        dims: ENERGY,
        samples: vec![
            ([0.0, 290.0], 1e6),
            ([1.0, 290.0], 3e6),
            ([0.0, 310.0], 2e6),
            ([1.0, 310.0], 6e6),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2]],
    });
    let law = chart(&c, [0.0, 1.0]);
    // The T=300 slice crosses two original triangles: q=1.5,2.5,4.5 MJ/kg.
    assert_eq!(law.knot_count(), 3);
    near(law.integral(0.0, 1.0).unwrap(), 2.75e6, 1e-9);
    let gap = card(PropertyValue::Surface2D {
        axes: ["moisture_content".into(), "temperature".into()],
        axis_dims: [Dims::NONE, Temperature::DIMS],
        dims: ENERGY,
        samples: vec![
            ([0.0, 290.0], 1e6),
            ([0.25, 290.0], 2e6),
            ([0.0, 310.0], 2e6),
            ([0.25, 310.0], 3e6),
            ([0.75, 290.0], 3e6),
            ([1.0, 290.0], 4e6),
            ([0.75, 310.0], 5e6),
            ([1.0, 310.0], 6e6),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
    });
    let begin = binding(&gap, &point(0.0));
    binding(&gap, &point(1.0));
    assert!(matches!(
        MoistureSorptionEnergy::from_binding(&gap, &begin, NAME, [0.0, 1.0]),
        Err(MaterialStatePointError::Support { .. })
    ));
}
