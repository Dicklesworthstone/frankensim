//! G1/G3: original property sources follow kinetic history, not equilibrium h.
//! Manufactured constitutive inputs verify the bridge, not a real polymer fit.

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::crystallization::properties::{
    CrystallizationMaterialPoint, CrystallizationPropertyError, CrystallizationPropertySource,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::crystallization::pressure::exchange::{
    IsobaricExchangeControls, IsobaricExchangeSideControls, exchange_heat_with_conductance,
};
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;
use fs_thermal::crystallization::{CrystallizationError, CrystallizationSolve};
use fs_thermal::enthalpy::series_contact_conductance_w_k;

const RANGE: [f64; 2] = [280.0, 420.0];

fn claim(name: &str, dims: Dims, axis: &str, knots: &[(f64, f64)]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: Dims::NONE,
            knots: knots.to_vec(),
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                RANGE[0],
                RANGE[1],
            )
            .with_quantity(
                axis,
                QuantitySpec::dimensional(Dims::NONE),
                knots[0].0,
                knots[knots.len() - 1].0,
            ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured kinetic-property fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured two-phase material".into(),
            phase: "amorphous solid and crystal".into(),
            process: "explicit kinetic history".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}

fn branches(card: &MaterialCard) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    [(20_000.0, 1000.0), (0.0, 2000.0)].map(|(reference, density)| {
        EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            card.content_hash(),
            // Both branches are SOLIDS: 1-X is not a liquid fraction.
            SolidLiquidPhase::Solid,
            reference,
            &RANGE.map(|temperature_k| HeatCapacityKnot {
                temperature_k,
                specific_heat_capacity_j_kg_k: 1000.0,
                bulk_density_kg_m3: density,
            }),
            1e-8,
        )
        .unwrap()
    })
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("X", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity("V", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap()
}

fn point<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> CrystallizationMaterialPoint<'a> {
    CrystallizationMaterialPoint {
        material: card,
        context,
        temperature_axis: "T",
        crystalline_mass_fraction_axis: Some("X"),
        crystalline_volume_fraction_axis: Some("V"),
    }
}

fn requirements() -> [ScalarPropertyRequirement; 2] {
    [
        ScalarPropertyRequirement::try_new(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "volume_observation",
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
    ]
}

fn material() -> MaterialCard {
    card(vec![
        claim(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            "X",
            &[(0.0, 0.2), (0.8, 1.0)],
        ),
        claim(
            "volume_observation",
            Dims::NONE,
            "V",
            &[(0.0, 0.0), (1.0, 1.0)],
        ),
    ])
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn kinetic_mass_and_volume_coordinates_share_original_claims_and_fixed_conditions() {
    let card = material();
    let [parent, crystal] = branches(&card);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, RANGE).unwrap();
    let rate = |_: f64| Ok(0.0);
    let law = pair
        .law(&rate, ContentHash([1; 32]), 1.0, 0.8, 420.0)
        .unwrap();
    let initial = law.state_at_temperature(350.0, 0.0).unwrap();
    let context = context();
    let source = CrystallizationPropertySource::new(
        pair,
        law,
        initial,
        point(&card, &context),
        &requirements(),
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
    )
    .unwrap();
    let partial = law
        .state_at_temperature(350.0, std::f64::consts::LN_2)
        .unwrap();
    let response = source.evaluate(partial).unwrap();
    assert_eq!(response.phase(), partial);
    close(response.volume().crystalline_mass_fraction, 0.4, 1e-15);
    close(response.volume().crystalline_volume_fraction, 0.25, 1e-15);
    close(response.volume().bulk_density_kg_m3, 1250.0, 1e-12);
    let bundle = response.material();
    close(
        bundle
            .property(THERMAL_CONDUCTIVITY_PROPERTY)
            .unwrap()
            .value_si(),
        0.6,
        1e-15,
    );
    close(
        bundle.property("volume_observation").unwrap().value_si(),
        0.25,
        1e-15,
    );
    assert_eq!(bundle.state_coordinate("pressure"), Some(100_000.0));
    assert_eq!(bundle.state_coordinate("T"), Some(350.0));
    close(bundle.state_coordinate("X").unwrap(), 0.4, 1e-15);
    close(bundle.state_coordinate("V").unwrap(), 0.25, 1e-15);
    for property in bundle.properties() {
        let original = source
            .original_source()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
    let amorphous = source.evaluate(initial).unwrap();
    close(
        amorphous
            .material()
            .property(THERMAL_CONDUCTIVITY_PROPERTY)
            .unwrap()
            .value_si(),
        0.2,
        1e-15,
    );
    // Equal enthalpy with a different history has a different T and conductivity.
    let same_h = law
        .state_at_enthalpy(
            initial.specific_enthalpy_j_kg(),
            partial.exposure(),
            CrystallizationSolve::default(),
            || false,
        )
        .unwrap();
    let same_h_response = source.evaluate(same_h).unwrap();
    close(same_h_response.phase().temperature_k(), 358.0, 1e-9);
    close(
        same_h_response
            .material()
            .property(THERMAL_CONDUCTIVITY_PROPERTY)
            .unwrap()
            .value_si(),
        0.6,
        1e-15,
    );
}

#[test]
fn foreign_histories_bad_axes_and_source_switches_refuse_without_mutation() {
    let card = material();
    let [parent, crystal] = branches(&card);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, RANGE).unwrap();
    let rate = |_: f64| Ok(0.0);
    let law = pair
        .law(&rate, ContentHash([1; 32]), 1.0, 0.8, 420.0)
        .unwrap();
    let initial = law.state_at_temperature(350.0, 0.0).unwrap();
    let context = context();
    let prepare = |point| {
        CrystallizationPropertySource::new(
            pair,
            law,
            initial,
            point,
            &requirements(),
            MaterialPropertySelection::SingleClaimOnly,
            CrystallizationSolve::default(),
        )
    };
    let source = prepare(point(&card, &context)).unwrap();
    let foreign = pair
        .law(&rate, ContentHash([2; 32]), 1.0, 0.8, 420.0)
        .unwrap();
    assert!(matches!(
        source.evaluate(foreign.state_at_temperature(350.0, 0.0).unwrap()),
        Err(CrystallizationPropertyError::Kinetic(
            CrystallizationError::StateMismatch
        ))
    ));
    let mut bad = point(&card, &context);
    bad.crystalline_volume_fraction_axis = Some("X");
    assert!(prepare(bad).is_err());
    bad = point(&card, &context);
    bad.crystalline_mass_fraction_axis = None;
    assert!(prepare(bad).is_err());
    bad = point(&card, &context);
    bad.temperature_axis = "missing";
    assert!(prepare(bad).is_err());
    assert_eq!(source.evaluate(initial).unwrap().phase(), initial);

    let split = card_for_switch();
    let [parent, crystal] = branches(&split);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, RANGE).unwrap();
    let law = pair
        .law(&rate, ContentHash([3; 32]), 1.0, 0.8, 420.0)
        .unwrap();
    let initial = law.state_at_temperature(350.0, 0.0).unwrap();
    let source = CrystallizationPropertySource::new(
        pair,
        law,
        initial,
        point(&split, &context),
        &requirements()[..1],
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
    )
    .unwrap();
    let original = source.original_source().clone();
    let later = law.state_at_temperature(350.0, 2.0).unwrap();
    assert!(later.crystalline_mass_fraction() > 0.5);
    // A second valid claim exists at this fraction, but refreshing may not switch.
    assert!(matches!(
        source.evaluate(later),
        Err(CrystallizationPropertyError::Material(_))
    ));
    assert_eq!(source.original_source(), &original);
    assert!(source.evaluate(initial).is_ok());
}

fn card_for_switch() -> MaterialCard {
    card(vec![
        claim(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            "X",
            &[(0.0, 0.2), (0.45, 0.65)],
        ),
        claim(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            "X",
            &[(0.5, 1.0), (0.8, 1.6)],
        ),
    ])
}

#[test]
fn finite_contact_uses_actual_candidate_crystallinity_and_conserves_heat() {
    let card = material();
    let [parent, crystal] = branches(&card);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, RANGE).unwrap();
    let rate = |t: f64| Ok(if t < 420.0 { 0.1 } else { 0.0 });
    let law = pair
        .law(&rate, ContentHash([4; 32]), 1.0, 0.8, 420.0)
        .unwrap();
    let cold = law.state_at_temperature(300.0, 0.0).unwrap();
    let hot = law.state_at_temperature(350.0, 0.0).unwrap();
    let context = context();
    let source = CrystallizationPropertySource::new(
        pair,
        law,
        hot,
        point(&card, &context),
        &requirements(),
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
    )
    .unwrap();
    let side = IsobaricExchangeSideControls::default();
    let left =
        IsobaricCrystallizationParcel::new(&law, hot, 1.0, 0.001, 100_000.0, side.caloric, || {
            false
        })
        .unwrap();
    let right =
        IsobaricCrystallizationParcel::new(&law, cold, 1.0, 0.001, 100_000.0, side.caloric, || {
            false
        })
        .unwrap();
    let mut conductance = |a: &IsobaricCrystallizationParcel<'_, '_>,
                           b: &IsobaricCrystallizationParcel<'_, '_>| {
        let values = [a, b].map(|parcel| {
            source
                .evaluate(parcel.state())
                .unwrap()
                .material()
                .property(THERMAL_CONDUCTIVITY_PROPERTY)
                .unwrap()
                .value_si()
        });
        series_contact_conductance_w_k(0.01, [0.005; 2], values, 0.0)
            .map_err(|_| CrystallizationError::Invalid("manufactured conduction geometry"))
    };
    let (left, right, report) = exchange_heat_with_conductance(
        &left,
        &right,
        &mut conductance,
        2.0,
        side,
        side,
        IsobaricExchangeControls {
            heat_bracket_j: [0.0, 100.0],
            residual_tolerance_j: 1e-6,
            energy_tolerance_j: 1e-5,
            maximum_evaluations: 96,
        },
        || false,
    )
    .unwrap();
    let x = 0.8 * (1.0 - (-0.2_f64).exp());
    let k = 0.2 + x;
    // Equal masses and Cp have deltaT_end = 50 - 2 Q/1000. Both
    // histories gain the same latent enthalpy, so it cancels in deltaT.
    let expected_heat = 2.0 * k * 50.0 / (1.0 + 4.0 * k / 1000.0);
    close(left.state().crystalline_mass_fraction(), x, 1e-12);
    close(right.state().crystalline_mass_fraction(), x, 1e-12);
    close(report.conductance_w_k, k, 1e-12);
    close(report.heat_left_to_right_j, expected_heat, 2e-6);
    close(report.combined_enthalpy_change_j, 0.0, 1e-8);
    close(
        report.combined_internal_energy_change_j + report.combined_boundary_work_j,
        0.0,
        1e-8,
    );
    close(
        conductance(&left, &right).unwrap(),
        report.conductance_w_k,
        1e-12,
    );
}
