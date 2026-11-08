use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, SolidLiquidPhase};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

const ROTATION: [[f64; 3]; 3] = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];

fn source() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, knots) in [
        ("kx", vec![(300.0, 2.0), (400.0, 10.0), (500.0, 4.0)]),
        ("ky", vec![(300.0, 3.0), (375.0, 9.0), (500.0, 5.0)]),
        ("kz", vec![(300.0, 6.0), (450.0, 3.0), (500.0, 12.0)]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, THERMAL_CONDUCTIVITY_DIMS),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute_temperature(),
                    300.0,
                    500.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured original directional conduction curves".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic anisotropic body".into(),
            phase: "declared solid and liquid branches".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn curve(card: &MaterialCard, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 300.0, 0.0),
            (400.0, 500.0, 0.0),
            (600.0, 500.0, 1.0),
            (800.0, 600.0, 1.0),
        ]
        .map(|(energy, temperature, fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: datum + energy,
            temperature_k: temperature,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: 1000.0,
        })
        .to_vec(),
    )
    .unwrap()
}

fn layer<'m>(
    card: &MaterialCard,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
    point: &QueryPoint,
    name: &str,
    range: [f64; 2],
) -> ResolvedConductionLayer<'m> {
    ResolvedConductionLayer::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new(name, THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        range,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn original_directional_integrals_rotate_without_copying_or_extrapolation() {
    let card = source();
    let caloric = curve(&card, 0.0);
    let point = QueryPoint::new()
        .with_quantity("T", absolute_temperature(), 300.0)
        .unwrap();
    let layers = [
        layer(&card, &caloric, &point, "kx", [300.0, 500.0]),
        layer(&card, &caloric, &point, "ky", [325.0, 475.0]),
        layer(&card, &caloric, &point, "kz", [350.0, 450.0]),
    ];
    let law = RotatedSourceConductivity::new(
        [&layers[0], &layers[1], &layers[2]],
        ROTATION,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(law.temperature_range_k(), [350.0, 450.0]);
    assert_eq!(law.specific_energy_range_j_kg(), [100.0, 300.0]);
    assert!(std::ptr::eq(law.layers()[1], &layers[1]));
    assert!(std::ptr::eq(law.caloric_curve(), &caloric));
    // Independent exact trapezoids over ORIGINAL differently located knots:
    // principal means = [8.25, 7.85, 4]. Apply the authored rotation by hand.
    let expected = [[8.106, 0.192, 0.0], [0.192, 7.994, 0.0], [0.0, 0.0, 4.0]];
    let mean = law.mean_conductivity(350.0, 450.0).unwrap().as_array();
    for i in 0..3 {
        for j in 0..3 {
            assert!((mean[i][j] - expected[i][j]).abs() < 5e-14);
        }
    }
    assert_eq!(
        mean,
        law.mean_conductivity(450.0, 350.0).unwrap().as_array()
    );
    assert_eq!(
        law.mean_conductivity(400.0, 400.0).unwrap(),
        law.conductivity_at_temperature(400.0).unwrap(),
    );
    assert_ne!(
        mean,
        law.conductivity_at_temperature(400.0).unwrap().as_array()
    );
    assert!(law.conductivity_at_temperature(349.0).is_err());
    assert!(law.mean_conductivity(350.0, 450.0_f64.next_up()).is_err());
}

#[test]
fn original_context_chart_support_frame_and_cancellation_are_admitted() {
    let card = source();
    let caloric = curve(&card, 0.0);
    let other_caloric = curve(&card, 10.0);
    let point = QueryPoint::new()
        .with_quantity("T", absolute_temperature(), 300.0)
        .unwrap();
    let original = layer(&card, &caloric, &point, "kx", [300.0, 500.0]);
    let different = layer(&card, &other_caloric, &point, "ky", [300.0, 500.0]);
    assert!(
        RotatedSourceConductivity::new(
            [&original, &different, &original],
            ROTATION,
            &CancelGate::new(),
        )
        .is_err()
    );
    let changed_context = point
        .clone()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            1e5,
        )
        .unwrap();
    let other_context = layer(&card, &caloric, &changed_context, "ky", [300.0, 500.0]);
    assert!(
        RotatedSourceConductivity::new(
            [&original, &other_context, &original],
            ROTATION,
            &CancelGate::new(),
        )
        .is_err()
    );
    let low = layer(&card, &caloric, &point, "kx", [300.0, 350.0]);
    let high = layer(&card, &caloric, &point, "ky", [400.0, 500.0]);
    assert!(
        RotatedSourceConductivity::new([&low, &high, &original], ROTATION, &CancelGate::new(),)
            .is_err()
    );
    let mut reflected = ROTATION;
    for row in &mut reflected {
        row[0] = -row[0];
    }
    assert!(
        RotatedSourceConductivity::new([&original; 3], reflected, &CancelGate::new(),).is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        RotatedSourceConductivity::new([&original; 3], ROTATION, &gate),
        Err(ChtError::Cancelled),
    ));
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn same_source_value_preserves_original_energy_support_and_solid_endpoint_tangent() {
    use super::super::super::super::TensorCaloricLaw;
    use fs_conduction::enthalpy::ThermalLaw;
    let card = source();
    let caloric = curve(&card, 0.0);
    let point = QueryPoint::new()
        .with_quantity("T", absolute_temperature(), 300.0)
        .unwrap();
    let layers =
        ["kx", "ky", "kz"].map(|name| layer(&card, &caloric, &point, name, [300.0, 500.0]));
    let law = RotatedSourceConductivity::new(
        [&layers[0], &layers[1], &layers[2]],
        ROTATION,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(law.specific_energy_bounds_j_kg(), [0.0, 400.0]);
    assert_eq!(
        law.temperature_and_slope_at_specific_energy(400.0).unwrap(),
        [500.0, 0.5],
    );
    assert_eq!(
        caloric
            .temperature_derivative_at_specific_enthalpy(400.0)
            .unwrap(),
        0.0
    );
    // Same temperature on the original latent branch does not extend the
    // pure-solid conductivity source's admitted energy domain.
    assert_eq!(caloric.temperature_k(450.0).unwrap(), 500.0);
    assert!(law.temperature_k(450.0).is_err());
    assert!(law.temperature_slope_kg_k_j(400.0_f64.next_up()).is_err());
    assert_eq!(law.bulk_density_kg_m3(100.0).unwrap(), Some(1000.0));
    let entropy = law
        .specific_entropy_change_j_kg_k(100.0, 200.0)
        .unwrap()
        .unwrap();
    assert!((entropy - 2.0 * (400.0_f64 / 350.0).ln()).abs() < 2e-13);
    assert!(law.specific_entropy_change_j_kg_k(100.0, 401.0).is_err());
}
