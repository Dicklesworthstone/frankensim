//! Fixed-context two-coordinate source data drives the existing enthalpy chart.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), t)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), p)
        .unwrap()
}

fn card(gap: bool, interpolation: InterpolationPolicy) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("specific_heat_capacity", CP, [2.0, 4.0, 4.0, 10.0]),
        ("density", DENSITY, [1000.0, 1000.0, 1100.0, 1100.0]),
    ] {
        let mut samples = [
            [300.0, 100.0],
            [400.0, 100.0],
            [300.0, 200.0],
            [400.0, 200.0],
        ]
        .into_iter()
        .zip(values)
        .collect();
        let mut triangles = vec![[0, 1, 2], [1, 3, 2]];
        if gap && name == "specific_heat_capacity" {
            let right = f64::from_bits(350_f64.to_bits() + 1);
            samples = [
                [300.0, 100.0],
                [350.0, 100.0],
                [300.0, 200.0],
                [350.0, 200.0],
                [right, 100.0],
                [400.0, 100.0],
                [right, 200.0],
                [400.0, 200.0],
            ]
            .map(|at| (at, 2.0))
            .to_vec();
            triangles = vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]];
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Surface2D {
                    axes: ["T".into(), "pressure".into()],
                    axis_dims: [Temperature::DIMS, Pressure::DIMS],
                    samples,
                    triangles,
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute(), 300.0, 400.0)
                    .with_quantity(
                        "pressure",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        100.0,
                        200.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic caloric surface; not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "analytic source fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(
    card: &MaterialCard,
    p: f64,
    reference: f64,
) -> Result<ResolvedHeatCapacityCurve, fs_material::phase::PhaseStateError> {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", CP),
            density: &PropertyKey::new("density", DENSITY),
            temperature_axis: "T",
            lower: &point(300.0, p),
            upper: &point(400.0, p),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        reference,
        0.1,
    )
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-10, "{a} != {b}");
}

#[test]
fn two_coordinate_capacity_retains_crossed_slopes_and_drives_conserved_enthalpy() {
    let card = card(false, InterpolationPolicy::LinearInside);
    for reference in [0.0, -1000.0] {
        let law = caloric(&card, 150.0, reference).unwrap();
        assert!(!law.curve().source_receipts().is_empty());
        assert_eq!(law.curve().source_receipts(), law.receipts());
        assert_eq!(law.curve().clone().source_receipts(), law.receipts());
        // Cp is 3→4 over [300,350], then 4→7 over [350,400].
        // Integrating only endpoint values would incorrectly give 500 J/kg.
        near(law.specific_enthalpy_change(300.0, 350.0).unwrap(), 175.0);
        near(law.specific_enthalpy_change(350.0, 400.0).unwrap(), 275.0);
        near(law.specific_enthalpy_change(300.0, 400.0).unwrap(), 450.0);
        near(law.specific_enthalpy_change(400.0, 300.0).unwrap(), -450.0);
        let initial = law.state_at_temperature(300.0).unwrap();
        let warm = law.curve().advance_specific_energy(initial, 175.0).unwrap();
        let hot = law.curve().advance_specific_energy(warm, 275.0).unwrap();
        near(warm.temperature_k(), 350.0);
        near(hot.temperature_k(), 400.0);
        near(
            hot.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg(),
            450.0,
        );
        near(hot.bulk_density_kg_m3(), 1050.0);
        assert_eq!(hot.material_card_identity(), card.content_hash());
        for receipt in law.receipts() {
            card.claims().verify_receipt(receipt).unwrap();
            assert!(receipt.query_point.contains(&("pressure".into(), 150.0)));
        }
    }
    // The same source data at another fixed context changes actual heat storage.
    let low = caloric(&card, 100.0, 0.0).unwrap();
    let high = caloric(&card, 200.0, 0.0).unwrap();
    near(low.specific_enthalpy_change(300.0, 400.0).unwrap(), 300.0);
    near(high.specific_enthalpy_change(300.0, 400.0).unwrap(), 700.0);
    assert_ne!(low.curve().identity(), high.curve().identity());
}

#[test]
fn capacity_paths_refuse_real_source_holes_and_discrete_only_samples() {
    assert!(caloric(&card(true, InterpolationPolicy::LinearInside), 150.0, 0.0).is_err());
    // Both endpoints are source vertices, but samples do not supply the path.
    assert!(caloric(&card(false, InterpolationPolicy::TabulatedOnly), 100.0, 0.0).is_err());
}

#[test]
fn untyped_temperature_cannot_erase_capacity_or_density_source_axis_dimensions() {
    let base = card(false, InterpolationPolicy::LinearInside);
    let keys = [
        PropertyKey::new("specific_heat_capacity", CP),
        PropertyKey::new("density", DENSITY),
    ];
    let query = |t| {
        QueryPoint::new()
            .with("T", t)
            .unwrap()
            .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 150.0)
            .unwrap()
    };
    for wrong_property in ["specific_heat_capacity", "density"] {
        let mut claims = ClaimSet::new();
        for key in &keys {
            let mut claim = base.claims().claims_for(key.name())[0].1.clone();
            claim.validity = ValidityDomain::unconstrained()
                .with("T", 300.0, 400.0)
                .with_quantity(
                    "pressure",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    100.0,
                    200.0,
                );
            if key.name() == wrong_property {
                let PropertyValue::Surface2D { axis_dims, .. } = &mut claim.value else {
                    unreachable!("fixture is a two-coordinate source");
                };
                axis_dims[0] = Pressure::DIMS;
            }
            claims.insert_claim(claim).unwrap();
        }
        let source = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
        for key in &keys {
            for t in [300.0, 350.0, 400.0] {
                source
                    .claims()
                    .query_typed(key, &query(t), SelectionPolicy::SingleClaimOnly)
                    .unwrap();
            }
        }
        let result = ResolvedHeatCapacityCurve::try_new(
            &source,
            HeatCapacitySourceQuery {
                heat_capacity: &keys[0],
                density: &keys[1],
                temperature_axis: "T",
                lower: &query(300.0),
                upper: &query(400.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Liquid,
            0.0,
            0.1,
        );
        assert!(
            matches!(result, Err(fs_material::phase::PhaseStateError::SourceQuery { property, reason })
            if property == wrong_property && reason == "source temperature coordinate must have temperature dimensions")
        );
    }
}
