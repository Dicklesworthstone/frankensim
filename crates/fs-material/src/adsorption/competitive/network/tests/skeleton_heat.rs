//! Manufactured caloric/flux references, not empirical porous-material data.
use super::*;
use crate::adsorption::competitive::network::conduction::{
    AdsorbentConductionLink, AdsorbentThermalContact,
};
use crate::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use crate::state_point::THERMAL_CONDUCTIVITY_DIMS;
use crate::thermal::{
    conduction::{LayerConductionConfig, LayerGeometry, source::ResolvedConductionLayer},
    interface::{
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        ThermalInterfacePoint,
    },
    source::ThermalMaterialPoint,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", temp_quantity(), 200.0, 800.0)
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                1000.0,
                1e7,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured finite skeleton conduction".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(name: &str, claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "analytic heat/site fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
pub(super) fn fixture_with_heat(
    k: f64,
    rate: f64,
) -> (
    MaterialCard,
    [CompetitiveSolid; 2],
    BTreeMap<u64, ClosedGasParcel>,
    ThermochemicalGasModel,
) {
    let (original, gas, catalog) = fixture([0.0, 0.0]);
    let cp = Dims([2, 0, -2, -1, 0, 0]);
    let rho = Dims([-3, 1, 0, 0, 0, 0]);
    let source = card(
        "shared finite dry skeleton",
        vec![
            claim("specific_heat_capacity", cp, 500.0),
            claim("density", rho, 1000.0),
            claim("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, k),
        ],
    );
    let dry = ResolvedHeatCapacityCurve::try_new(
        &source,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", cp),
            density: &PropertyKey::new("density", rho),
            temperature_axis: "temperature",
            lower: &point(200.0, 1e5),
            upper: &point(800.0, 1e5),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    let original = &original.model.channels[0];
    let law = LangmuirAdsorbent::try_new(
        &catalog,
        original.species.clone(),
        AdsorptionSource {
            material: card(
                "explicit single A site source",
                vec![
                    claim(LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_CAPACITY_DIMS, 0.5),
                    claim(
                        LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                        LANGMUIR_AFFINITY_DIMS,
                        2e-5,
                    ),
                    claim(LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, rate),
                ],
            ),
            ..original.source.clone()
        },
        AdsorptionCaloricBinding {
            dry_skeleton: dry,
            ..original.calorics.clone()
        },
    )
    .unwrap();
    let model = Arc::new(CompetitiveAdsorbent::try_new(vec![law]).unwrap());
    let solids = [301.0, 299.0].map(|temperature| {
        CompetitiveSolid::try_new(
            Arc::clone(&model),
            0.1,
            &[(SpeciesId::new("A").unwrap(), 0.0)],
            temperature,
        )
        .unwrap()
    });
    (
        source,
        solids,
        BTreeMap::from([(0, gas.clone()), (1, gas)]),
        catalog,
    )
}
pub(super) fn layer<'a>(
    solid: &'a CompetitiveSolid,
    card: &'a MaterialCard,
    context: &'a QueryPoint,
    range: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        solid.model().dry_skeleton_calorics().curve(),
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "temperature",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        range,
        || Ok(()),
    )
    .unwrap()
}
pub(super) fn heat_link<'a, 'm>(
    layers: [&'a ResolvedConductionLayer<'m>; 2],
    contact: Option<&'a AdsorbentThermalContact<'m>>,
) -> AdsorbentConductionLink<'a, 'm> {
    AdsorbentConductionLink {
        id: 41,
        solids: [0, 1],
        layers,
        geometry: LayerGeometry::planar([0.005, 0.005], 0.01).unwrap(),
        contact,
        solve: LayerConductionConfig {
            max_iterations: 128,
            absolute_heat_tolerance_w: 1e-10,
            relative_heat_tolerance: 1e-12,
        },
    }
}
pub(super) fn contact_card(
    material: &MaterialCard,
    resistance: f64,
    max_b: f64,
) -> InterfaceSystemCard {
    let mut claim = claim(
        AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        resistance,
    );
    claim.validity = ValidityDomain::unconstrained()
        .with_quantity("temperature_a", temp_quantity(), 200.0, 800.0)
        .with_quantity("temperature_b", temp_quantity(), 200.0, max_b);
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    let surface = || SurfaceSpec {
        material: material.id().clone(),
        texture_frame: "explicit ordered solid face".into(),
    };
    InterfaceSystemCard::assemble(
        surface(),
        surface(),
        SystemContext {
            medium: "explicit finite thermal contact".into(),
            third_body: None,
            environment: "manufactured caloric fixture".into(),
            history: "fixed original contact".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub(super) fn contact_point() -> QueryPoint {
    point(300.0, 1e5)
        .with_quantity("temperature_a", temp_quantity(), 301.0)
        .unwrap()
        .with_quantity("temperature_b", temp_quantity(), 299.0)
        .unwrap()
}
pub(super) fn contact<'a>(
    solids: &[CompetitiveSolid; 2],
    source: &'a InterfaceSystemCard,
    material: &'a MaterialCard,
    context: &'a QueryPoint,
) -> AdsorbentThermalContact<'a> {
    AdsorbentThermalContact::try_new(
        ThermalInterfacePoint {
            interface: source,
            materials: [material, material],
            context,
            temperature_axes: ["temperature_a", "temperature_b"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        solids.each_ref().map(|solid| {
            solid
                .model()
                .dry_skeleton_calorics()
                .state_at_temperature(solid.temperature_k())
                .unwrap()
        }),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn bindings(
    solids: &[CompetitiveSolid; 2],
    solve: CompetitiveSolve,
) -> [AdsorptionNetworkBinding<'_>; 2] {
    std::array::from_fn(|i| AdsorptionNetworkBinding {
        id: i as u64,
        parcel: i as u64,
        solid: &solids[i],
        solve,
    })
}

#[test]
fn original_layer_and_contact_heat_use_one_finite_caloric_inversion() {
    let (card, solids, gases, catalog) = fixture_with_heat(20.0, 0.0);
    let saved = (solids.clone(), gases.clone());
    let context = point(300.0, 1e5);
    let layers = solids
        .each_ref()
        .map(|solid| layer(solid, &card, &context, [200.0, 800.0]));
    let source = contact_card(&card, 0.002, 800.0);
    let contact_context = contact_point();
    let contact = contact(&solids, &source, &card, &contact_context);
    let bindings = bindings(&solids, pair_controls());
    let dt = 0.001;
    let run = |links: &[AdsorbentConductionLink<'_, '_>]| {
        step_candidate_with_conduction(
            &gases,
            GasTransportPaths::default(),
            &bindings,
            &catalog,
            dt,
            &endpoints(&gases),
            controls(),
            links,
            || false,
        )
        .unwrap()
    };
    assert_eq!(
        run(&[]),
        step_candidate(
            &gases,
            GasTransportPaths::default(),
            &bindings,
            &catalog,
            dt,
            &endpoints(&gases),
            controls(),
            || false
        )
        .unwrap()
    );
    for (contact, resistance) in [(None, 0.0), (Some(&contact), 0.002)] {
        let link = heat_link(layers.each_ref(), contact);
        let result = run(&[link]);
        let heat = dt * 0.01 * 2.0 / (0.005 / 20.0 + 0.005 / 20.0 + resistance);
        near(
            result.solids()[&0].temperature_k(),
            301.0 - heat / 50.0,
            2e-10,
        );
        near(
            result.solids()[&1].temperature_k(),
            299.0 + heat / 50.0,
            2e-10,
        );
        near(result.step().solids[0].conductive_heat_j, -heat, 1e-12);
        near(result.step().solids[1].conductive_heat_j, heat, 1e-12);
        near(
            result
                .step()
                .solids
                .iter()
                .map(|s| s.solid_energy_change_j)
                .sum(),
            0.0,
            2e-9,
        );
        near(
            result.step().total_energy_residual_j,
            0.0,
            controls().exchange.energy_tolerance_j,
        );
        assert!(result.step().entropy_change_j_k > 0.0);
        for (&id, before) in &gases {
            near(result.parcels()[&id].mass_kg(), before.mass_kg(), 1e-13);
            near(
                result.parcels()[&id].state().state().temperature,
                300.0,
                2e-7,
            );
            near(
                result.parcels()[&id].state().state().pressure,
                100000.0,
                1e-6,
            );
        }
        let sample = &result.step().conduction[0];
        near(sample.energy_a_to_b_j, heat, 1e-12);
        assert_eq!(sample.temperatures_k, [301.0, 299.0]);
        assert!(Arc::ptr_eq(
            &sample.source_lines[0],
            layers[0].source_line()
        ));
        assert_eq!(
            sample.caloric_curve_identities[0],
            layers[0].caloric_curve().identity()
        );
        assert_eq!(sample.interface.is_some(), contact.is_some());
        if let Some(contact) = contact {
            assert_eq!(
                sample.interface.as_ref().unwrap().properties()[0]
                    .answer()
                    .receipt
                    .selected,
                contact.source().properties()[0].answer().receipt.selected
            );
        }
    }
    assert_eq!((solids, gases), saved);
}

#[test]
fn data_only_conductivity_swap_changes_actual_finite_temperature() {
    let run = |conductivity| {
        let (card, solids, gases, catalog) = fixture_with_heat(conductivity, 0.0);
        let context = point(300.0, 1e5);
        let layers = solids
            .each_ref()
            .map(|solid| layer(solid, &card, &context, [200.0, 800.0]));
        step_candidate_with_conduction(
            &gases,
            GasTransportPaths::default(),
            &bindings(&solids, pair_controls()),
            &catalog,
            0.001,
            &endpoints(&gases),
            controls(),
            &[heat_link(layers.each_ref(), None)],
            || false,
        )
        .unwrap()
    };
    let low = run(20.0);
    let high = run(80.0);
    near(
        high.step().conduction[0].energy_a_to_b_j,
        4.0 * low.step().conduction[0].energy_a_to_b_j,
        1e-12,
    );
    near(
        301.0 - high.solids()[&0].temperature_k(),
        4.0 * (301.0 - low.solids()[&0].temperature_k()),
        1e-9,
    );
    assert_ne!(
        low.step().conduction[0].material_card_identities,
        high.step().conduction[0].material_card_identities
    );
}

#[test]
fn original_skeleton_identity_and_late_source_domains_are_mandatory() {
    let (card, solids, gases, catalog) = fixture_with_heat(20.0, 0.4);
    let saved = (solids.clone(), gases.clone());
    let context = point(300.0, 1e5);
    let layers = solids
        .each_ref()
        .map(|solid| layer(solid, &card, &context, [299.0, 301.0001]));
    let bindings = bindings(&solids, pair_controls());
    let link = heat_link(layers.each_ref(), None);
    let result = step_candidate_with_conduction(
        &gases,
        GasTransportPaths::default(),
        &bindings,
        &catalog,
        0.001,
        &endpoints(&gases),
        controls(),
        &[link],
        || false,
    );
    assert!(matches!(result, Err(Error::Conduction { id: 41, .. })));
    // An equal-valued clone is not the retained finite storage object.
    let foreign = solids[0].model().dry_skeleton_calorics().clone();
    let foreign_layer = ResolvedConductionLayer::new(
        foreign.curve(),
        ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "temperature",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        [200.0, 800.0],
        || Ok(()),
    )
    .unwrap();
    let link = heat_link([&foreign_layer, &layers[1]], None);
    assert!(matches!(
        step_candidate_with_conduction(
            &gases,
            GasTransportPaths::default(),
            &bindings,
            &catalog,
            0.001,
            &endpoints(&gases),
            controls(),
            &[link],
            || false
        ),
        Err(Error::Conduction { id: 41, .. })
    ));
    assert_eq!((solids, gases), saved);
}

#[test]
fn excessive_positive_heat_refuses_instead_of_clipping_the_transfer() {
    let (card, solids, gases, catalog) = fixture_with_heat(20.0, 0.0);
    let saved = (solids.clone(), gases.clone());
    let context = point(300.0, 1e5);
    let layers = solids
        .each_ref()
        .map(|solid| layer(solid, &card, &context, [200.0, 800.0]));
    let mut solve = pair_controls();
    solve.maximum_temperature_change_k = 10.0;
    let link = heat_link(layers.each_ref(), None);
    // Q=120 J moves both 50 J/K stores by 2.4 K. The unclipped explicit
    // candidate increases the temperature separation and fails total entropy.
    assert!(matches!(
        step_candidate_with_conduction(
            &gases,
            GasTransportPaths::default(),
            &bindings(&solids, solve),
            &catalog,
            3.0,
            &endpoints(&gases),
            controls(),
            &[link],
            || false
        ),
        Err(Error::Thermodynamic(
            AdsorptionError::EntropyDecreased { .. }
        ))
    ));
    assert_eq!((solids, gases), saved);
}
