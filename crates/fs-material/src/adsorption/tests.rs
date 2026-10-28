//! G0/G1/G3 synthetic Langmuir kinetics and conservative finite-gas coupling.
//! These fixtures prove the equations and accounting, not a named adsorbent.
use super::*;
use crate::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::tests::{model, point, species_with_energy_reference, transport};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const T0: f64 = 300.0;
const PRESSURE: f64 = 100_000.0;
const MASS: f64 = 0.028;
const HEAT: f64 = 6000.0;

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
            source: "synthetic reversible Langmuir and caloric fixture".into(),
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
            process: "analytic adsorption fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn chart(name: &str, cp: f64, datum: f64) -> ResolvedHeatCapacityCurve {
    let cp_dims = Dims([2, 0, -2, -1, 0, 0]);
    let rho_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let source = card(
        name,
        vec![
            claim("specific_heat_capacity", cp_dims, cp),
            claim("density", rho_dims, 1000.0),
        ],
    );
    ResolvedHeatCapacityCurve::try_new(
        &source,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", cp_dims),
            density: &PropertyKey::new("density", rho_dims),
            temperature_axis: "temperature",
            lower: &point(200.0, PRESSURE),
            upper: &point(800.0, PRESSURE),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        datum,
        1e-6,
    )
    .unwrap()
}
fn rate_source(affinity: f64, rate: f64) -> AdsorptionSource {
    AdsorptionSource {
        material: card(
            "synthetic adsorbent",
            vec![
                claim(LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_CAPACITY_DIMS, 0.5),
                claim(
                    LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                    LANGMUIR_AFFINITY_DIMS,
                    affinity,
                ),
                claim(LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, rate),
            ],
        ),
        point: point(T0, PRESSURE),
        selection: MaterialPropertySelection::SingleClaimOnly,
        reference_temperature_k: T0,
    }
}
fn calorics(
    catalog: &ThermochemicalGasModel,
    bound_datum: f64,
    heat: f64,
) -> AdsorptionCaloricBinding {
    let bound = chart("bound synthetic A", 3.5 * R / MASS, bound_datum);
    let hg = catalog
        .species_specific_enthalpy_j_kg(&SpeciesId::new("A").unwrap(), T0)
        .unwrap();
    let hb = bound
        .state_at_temperature(T0)
        .unwrap()
        .specific_enthalpy_j_kg();
    AdsorptionCaloricBinding {
        dry_skeleton: chart("dry skeleton", 500.0, 0.0),
        bound_adsorbate: bound,
        // Explicit synthetic heat calibration, not an inferred production heat.
        bound_enthalpy_offset_j_kg: hg - hb - heat / MASS,
    }
}
fn adsorbent(catalog: &ThermochemicalGasModel, affinity: f64, rate: f64) -> Arc<LangmuirAdsorbent> {
    Arc::new(
        LangmuirAdsorbent::try_new(
            catalog,
            SpeciesId::new("A").unwrap(),
            rate_source(affinity, rate),
            calorics(catalog, 0.0, HEAT),
        )
        .unwrap(),
    )
}
fn gas(catalog: &ThermochemicalGasModel) -> ClosedGasParcel {
    ClosedGasParcel::new(catalog.clone(), point(T0, PRESSURE), 0.001).unwrap()
}
fn controls() -> AdsorptionSolve {
    AdsorptionSolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 800.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        maximum_gas_gross_fraction: 0.1,
        maximum_site_turnover: 0.25,
        maximum_gas_surface_temperature_difference_k: 5.0,
        maximum_temperature_change_k: 2.0,
        equilibrium_coverage_tolerance: 1e-12,
        entropy_tolerance_j_k: 1e-10,
        maximum_iterations: 128,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn exact_frozen_langmuir_relaxation_retains_both_gross_streams_and_empty_site_turnover() {
    let capacity = 2.0;
    let initial = 0.4;
    let dt = 0.3;
    let result = integrate_sites(initial, capacity, 0.25, 0.1, 4.0, dt).unwrap();
    let target = 0.5;
    let lambda = 0.2;
    let expected = capacity * (target + (initial / capacity - target) * (-lambda * dt).exp());
    let integrated =
        target * dt + (initial / capacity - target) * (1.0 - (-lambda * dt).exp()) / lambda;
    near(result.final_moles, expected, 2e-15);
    near(result.desorption, capacity * 0.1 * integrated, 2e-15);
    near(result.adsorption, capacity * 0.1 * (dt - integrated), 2e-15);
    near(
        result.adsorption - result.desorption,
        result.final_moles - initial,
        2e-15,
    );
    let empty = integrate_sites(0.0, capacity, 0.25, 0.1, 4.0, dt).unwrap();
    assert!(empty.final_moles > 0.0 && empty.desorption > 0.0);
    assert!(empty.adsorption > empty.desorption);
    near(
        empty.adsorption - empty.desorption,
        empty.final_moles,
        2e-15,
    );
    let tiny = integrate_sites(0.0, capacity, 0.25, 0.1, 4.0, 1e-8).unwrap();
    assert!(tiny.desorption > 0.0);
    near(tiny.final_moles, capacity * 0.1 * 1e-8, 1e-17);
}

#[test]
fn thermodynamic_affinity_matches_constant_heat_vant_hoff_and_source_units() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let law = adsorbent(&catalog, 2e-5, 0.4);
    for temperature in [250.0, 300.0, 350.0, 700.0] {
        let properties = law.properties(temperature, PRESSURE).unwrap();
        let expected = 2e-5 * (HEAT / R * (1.0 / temperature - 1.0 / T0)).exp();
        near(properties.affinity_pa_inverse, expected, 2e-12 * expected);
        near(properties.adsorption_heat_j_mol, HEAT, 2e-9);
        near(properties.desorption_rate_s_inverse, 0.4, 0.0);
    }
    assert_ne!(LANGMUIR_CAPACITY_DIMS, Dims([0, 0, 0, 0, 0, 0]));
    assert_ne!(LANGMUIR_AFFINITY_DIMS, Dims([-3, 0, 0, 0, 0, 1]));
}

#[test]
fn finite_uptake_depletes_only_the_adsorbate_and_heats_the_real_solid_without_a_thermostat() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let original_gas = gas(&catalog);
    let original = AdsorptionSolid::try_new(adsorbent(&catalog, 2e-5, 0.4), 0.1, 0.0, T0).unwrap();
    let saved = (original.clone(), original_gas.clone());
    let result = advance_candidate(
        &original,
        &original_gas,
        &catalog,
        0.005,
        controls(),
        || false,
    )
    .unwrap();
    let n = result.solid().adsorbed_moles();
    assert!(n > 0.0);
    let a = SpeciesId::new("A").unwrap();
    let b = SpeciesId::new("B").unwrap();
    near(
        original_gas.species_mass_kg(&a).unwrap() - result.gas().species_mass_kg(&a).unwrap(),
        n * MASS,
        1e-13,
    );
    near(
        result.gas().species_mass_kg(&b).unwrap(),
        original_gas.species_mass_kg(&b).unwrap(),
        1e-13,
    );
    assert_eq!(result.solid().dry_mass_kg(), original.dry_mass_kg());
    assert!(result.solid().temperature_k() > T0);
    assert!(
        result.step().iterations <= 3,
        "constant-Cp solid should use its exact caloric slope"
    );
    assert!(result.gas().state().state().temperature < T0);
    assert!(result.step().entropy_change_j_k > 0.0);
    assert!(result.step().energy_residual_j.abs() < controls().gas.energy_tolerance_j);
    assert!(result.step().mass_residual_kg.abs() < controls().gas.mass_tolerance_kg);
    near(
        result.step().solid_energy_change_j + result.step().gas_energy_change_j,
        0.0,
        2e-8,
    );
    near(result.step().energy_into_gas_j, -n * 3.5 * R * T0, 2e-10);
    near(
        result.step().dry_sensible_energy_change_j + result.step().bound_sensible_energy_change_j,
        n * HEAT,
        1e-8,
    );
    assert_eq!((&original, &original_gas), (&saved.0, &saved.1));
}

#[test]
fn desorption_reintroduces_an_absent_gas_species_and_cools_the_finite_solid() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let carrier = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap())
        .unwrap();
    let original_gas = gas(&carrier);
    let original =
        AdsorptionSolid::try_new(adsorbent(&catalog, 2e-5, 0.4), 0.1, 0.025, T0).unwrap();
    let result = advance_candidate(
        &original,
        &original_gas,
        &catalog,
        0.002,
        controls(),
        || false,
    )
    .unwrap();
    let a = SpeciesId::new("A").unwrap();
    assert_eq!(original_gas.species_mass_kg(&a).unwrap(), 0.0);
    assert!(result.gas().species_mass_kg(&a).unwrap() > 0.0);
    assert!(result.solid().adsorbed_moles() < original.adsorbed_moles());
    assert!(result.solid().temperature_k() < T0);
    assert!(result.gas().state().state().temperature > T0);
    assert_eq!(result.step().gross_adsorption_mol, 0.0);
    assert!(result.step().gross_desorption_mol > 0.0);
    assert!(result.step().entropy_change_j_k > 0.0);
    near(
        result.gas().species_mass_kg(&a).unwrap(),
        (original.adsorbed_moles() - result.solid().adsorbed_moles()) * MASS,
        1e-13,
    );
    let continued = advance_candidate(
        result.solid(),
        result.gas(),
        &catalog,
        0.002,
        controls(),
        || false,
    )
    .unwrap();
    assert!(continued.solid().adsorbed_moles() < result.solid().adsorbed_moles());
    assert!(continued.step().energy_residual_j.abs() < controls().gas.energy_tolerance_j);
}

#[test]
fn equilibrium_has_gross_reversible_turnover_but_zero_net_mass_heat_and_entropy() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let original_gas = gas(&catalog);
    let original =
        AdsorptionSolid::try_new(adsorbent(&catalog, 2e-5, 0.4), 0.1, 0.025, T0).unwrap();
    let result = advance_candidate(
        &original,
        &original_gas,
        &catalog,
        0.005,
        controls(),
        || false,
    )
    .unwrap();
    assert!(result.step().gross_adsorption_mol > 0.0);
    near(
        result.step().gross_adsorption_mol,
        result.step().gross_desorption_mol,
        0.0,
    );
    assert_eq!(result.solid(), &original);
    assert_eq!(result.gas(), &original_gas);
    assert_eq!(result.step().energy_into_gas_j, 0.0);
    assert_eq!(result.step().entropy_change_j_k, 0.0);
    let zero = advance_candidate(&original, &original_gas, &catalog, 0.0, controls(), || {
        false
    })
    .unwrap();
    assert_eq!(zero.solid(), &original);
    assert_eq!(zero.gas(), &original_gas);
    assert_eq!(zero.step().gross_adsorption_mol, 0.0);
}

#[test]
fn bound_caloric_datum_shift_is_canceled_by_the_explicit_alignment() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let original_gas = gas(&catalog);
    let make = |datum| {
        Arc::new(
            LangmuirAdsorbent::try_new(
                &catalog,
                SpeciesId::new("A").unwrap(),
                rate_source(2e-5, 0.4),
                calorics(&catalog, datum, HEAT),
            )
            .unwrap(),
        )
    };
    let a = AdsorptionSolid::try_new(make(0.0), 0.1, 0.0, T0).unwrap();
    let b = AdsorptionSolid::try_new(make(1e6), 0.1, 0.0, T0).unwrap();
    let pa = advance_candidate(&a, &original_gas, &catalog, 0.005, controls(), || false).unwrap();
    let pb = advance_candidate(&b, &original_gas, &catalog, 0.005, controls(), || false).unwrap();
    near(pa.solid().temperature_k(), pb.solid().temperature_k(), 2e-8);
    near(
        pa.gas().state().state().temperature,
        pb.gas().state().state().temperature,
        2e-8,
    );
    near(
        pa.step().entropy_change_j_k,
        pb.step().entropy_change_j_k,
        1e-9,
    );
    near(
        pa.solid().adsorbed_moles(),
        pb.solid().adsorbed_moles(),
        0.0,
    );
}

#[test]
fn finite_bath_inventory_thermal_support_and_changed_sources_refuse_without_partial_updates() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let original_gas = gas(&catalog);
    let original = AdsorptionSolid::try_new(adsorbent(&catalog, 2e-5, 0.4), 0.1, 0.0, T0).unwrap();
    let saved = (original.clone(), original_gas.clone());
    assert!(matches!(
        advance_candidate(
            &original,
            &original_gas,
            &catalog,
            100.0,
            controls(),
            || false
        ),
        Err(AdsorptionError::InventoryLimit { .. })
    ));
    let strict = AdsorptionSolve {
        maximum_temperature_change_k: 0.0,
        ..controls()
    };
    assert!(
        advance_candidate(&original, &original_gas, &catalog, 0.005, strict, || false).is_err()
    );
    let hot = ClosedGasParcel::new(catalog.clone(), point(350.0, PRESSURE), 0.001).unwrap();
    assert!(advance_candidate(&original, &hot, &catalog, 0.005, controls(), || false).is_err());
    let changed = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", MASS, 3.5, 0.0, 1.0),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap();
    assert!(
        advance_candidate(
            &original,
            &original_gas,
            &changed,
            0.005,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!((&original, &original_gas), (&saved.0, &saved.1));
}

#[test]
fn source_units_aliases_and_out_of_domain_states_are_not_silently_accepted() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let bad = AdsorptionSource {
        material: card(
            "wrong capacity basis",
            vec![
                claim(LANGMUIR_CAPACITY_PROPERTY, Dims([0, 0, 0, 0, 0, 0]), 0.5),
                claim(
                    LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                    LANGMUIR_AFFINITY_DIMS,
                    2e-5,
                ),
                claim(LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, 0.4),
            ],
        ),
        ..rate_source(2e-5, 0.4)
    };
    assert!(
        LangmuirAdsorbent::try_new(
            &catalog,
            SpeciesId::new("A").unwrap(),
            bad,
            calorics(&catalog, 0.0, HEAT)
        )
        .is_err()
    );
    let mut legacy = claim(LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, 0.4);
    legacy.validity = legacy.validity.with("T", 290.0, 310.0);
    let bad = AdsorptionSource {
        material: card(
            "stale T alias",
            vec![
                claim(LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_CAPACITY_DIMS, 0.5),
                claim(
                    LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                    LANGMUIR_AFFINITY_DIMS,
                    2e-5,
                ),
                legacy,
            ],
        ),
        point: point(T0, PRESSURE).with("T", 300.0).unwrap(),
        ..rate_source(2e-5, 0.4)
    };
    assert!(
        LangmuirAdsorbent::try_new(
            &catalog,
            SpeciesId::new("A").unwrap(),
            bad,
            calorics(&catalog, 0.0, HEAT)
        )
        .is_err()
    );
    // Semantic aliases must not freeze another owned temperature/pressure
    // while the canonical gas/surface state advances.
    for (axis, quantity, value, range) in [
        ("sample_temperature", temp_quantity(), T0, [200.0, 800.0]),
        (
            "sample_pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            PRESSURE,
            [1000.0, 1e7],
        ),
    ] {
        let mut rate = claim(LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, 0.4);
        rate.validity = rate
            .validity
            .with_quantity(axis, quantity, range[0], range[1]);
        let aliased = AdsorptionSource {
            material: card(
                "stale semantic source alias",
                vec![
                    claim(LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_CAPACITY_DIMS, 0.5),
                    claim(
                        LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                        LANGMUIR_AFFINITY_DIMS,
                        2e-5,
                    ),
                    rate,
                ],
            ),
            point: point(T0, PRESSURE)
                .with_quantity(axis, quantity, value)
                .unwrap(),
            ..rate_source(2e-5, 0.4)
        };
        assert!(
            LangmuirAdsorbent::try_new(
                &catalog,
                SpeciesId::new("A").unwrap(),
                aliased,
                calorics(&catalog, 0.0, HEAT)
            )
            .is_err()
        );
    }
    let law = adsorbent(&catalog, 2e-5, 0.4);
    assert!(law.properties(199.0, PRESSURE).is_err());
    assert!(law.properties(T0, 100.0).is_err());
    assert!(AdsorptionSolid::try_new(Arc::clone(&law), 0.1, 0.051, T0).is_err());
    // A finite positive loading must not disappear when normalized by capacity.
    assert!(AdsorptionSolid::try_new(Arc::clone(&law), 2e300, 1e-30, T0).is_err());
    assert!(AdsorptionSolid::try_new(law, 0.1, -0.001, T0).is_err());
    // Nonzero forward kinetics cannot underflow into a fictitious zero stream.
    assert!(integrate_sites(0.0, 1.0, 1e-300, 1e-300, 1.0, 1.0).is_err());
    assert!(
        AdsorptionSolve {
            maximum_gas_gross_fraction: 1.0,
            ..controls()
        }
        .validate(0.1)
        .is_err()
    );
    assert!(
        AdsorptionSolve {
            entropy_tolerance_j_k: -1.0,
            ..controls()
        }
        .validate(0.1)
        .is_err()
    );
}

#[test]
fn cancellation_at_kinetics_thermal_and_final_audit_boundaries_preserves_both_owners() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let original_gas = gas(&catalog);
    let original = AdsorptionSolid::try_new(adsorbent(&catalog, 2e-5, 0.4), 0.1, 0.0, T0).unwrap();
    let saved = (original.clone(), original_gas.clone());
    let mut calls = 0;
    advance_candidate(
        &original,
        &original_gas,
        &catalog,
        0.005,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert!(matches!(
            advance_candidate(
                &original,
                &original_gas,
                &catalog,
                0.005,
                controls(),
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(AdsorptionError::Cancelled)
        ));
    }
    assert_eq!((&original, &original_gas), (&saved.0, &saved.1));
}

#[test]
fn unequal_caloric_heat_capacities_preserve_vant_hoff_and_the_langmuir_entropy_force() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let mut binding = calorics(&catalog, 0.0, HEAT);
    binding.bound_adsorbate = chart("bound A with distinct Cp", 2.0 * R / MASS, 0.0);
    let id = SpeciesId::new("A").unwrap();
    binding.bound_enthalpy_offset_j_kg = catalog.species_specific_enthalpy_j_kg(&id, T0).unwrap()
        - binding
            .bound_adsorbate
            .state_at_temperature(T0)
            .unwrap()
            .specific_enthalpy_j_kg()
        - HEAT / MASS;
    let law = Arc::new(
        LangmuirAdsorbent::try_new(&catalog, id, rate_source(2e-5, 0.4), binding).unwrap(),
    );
    let t = 400.0;
    let h = 0.01;
    let properties = law.properties(t, PRESSURE).unwrap();
    near(
        properties.adsorption_heat_j_mol,
        HEAT + 1.5 * R * (t - T0),
        2e-9,
    );
    let derivative = (law
        .properties(t + h, PRESSURE)
        .unwrap()
        .affinity_pa_inverse
        .ln()
        - law
            .properties(t - h, PRESSURE)
            .unwrap()
            .affinity_pa_inverse
            .ln())
        / (2.0 * h);
    near(
        derivative,
        -properties.adsorption_heat_j_mol / (R * t * t),
        1e-10,
    );

    // Independent chemical-potential identity at equal initial temperatures:
    // dS_total/dn_bound = R ln[b p_i (1-theta)/theta]. The full candidate
    // evaluates gas pressure entropy, both calorics and site mixing separately.
    let original_gas = gas(&catalog);
    let original = AdsorptionSolid::try_new(law, 0.1, 0.01, T0).unwrap();
    let mut solve = controls();
    solve.gas.energy_tolerance_j = 1e-10;
    solve.gas.gas.energy_tolerance_j_kg = 1e-8;
    let result =
        advance_candidate(&original, &original_gas, &catalog, 1e-5, solve, || false).unwrap();
    let dn = result.solid().adsorbed_moles() - original.adsorbed_moles();
    let expected_force =
        R * (2e-5 * PRESSURE * 0.5 * (1.0 - original.coverage()) / original.coverage()).ln();
    assert!(dn > 0.0 && result.step().entropy_change_j_k > 0.0);
    near(result.step().entropy_change_j_k / dn, expected_force, 2e-4);
}

#[test]
fn safeguarded_caloric_inverse_crosses_source_chart_knots_in_both_directions() {
    let catalog = model(vec![0.5, 0.5], 0.0);
    let cp_dims = Dims([2, 0, -2, -1, 0, 0]);
    let rho_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let mut cp = claim("specific_heat_capacity", cp_dims, 500.0);
    cp.value = PropertyValue::Curve {
        abscissa: "temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![
            (200.0, 500.0),
            (300.0, 500.0),
            (301.0, 1000.0),
            (800.0, 1000.0),
        ],
        dims: cp_dims,
    };
    cp.interpolation = InterpolationPolicy::LinearInside;
    let material = card(
        "piecewise caloric inverse fixture",
        vec![cp, claim("density", rho_dims, 1000.0)],
    );
    let dry = ResolvedHeatCapacityCurve::try_new(
        &material,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", cp_dims),
            density: &PropertyKey::new("density", rho_dims),
            temperature_axis: "temperature",
            lower: &point(200.0, PRESSURE),
            upper: &point(800.0, PRESSURE),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        100.0,
    )
    .unwrap();
    // This declared interpolation allowance retains the four original knots:
    // exact chart chord slopes are 500, 750 and 1000 J/(kg K).
    assert_eq!(dry.curve().knots().len(), 4);
    let mut binding = calorics(&catalog, 0.0, HEAT);
    binding.dry_skeleton = dry;
    let law = Arc::new(
        LangmuirAdsorbent::try_new(
            &catalog,
            SpeciesId::new("A").unwrap(),
            rate_source(2e-5, 0.4),
            binding,
        )
        .unwrap(),
    );
    let solid = AdsorptionSolid::try_new(law, 0.1, 0.01, T0).unwrap();
    let final_bound_moles = 0.02;
    for target in [250.0, 300.01, 600.0] {
        let dry_change = if target <= 300.0 {
            500.0 * (target - 300.0)
        } else if target <= 301.0 {
            750.0 * (target - 300.0)
        } else {
            750.0 + 1000.0 * (target - 301.0)
        };
        let heat = 0.1 * dry_change + final_bound_moles * 3.5 * R * (target - T0);
        let (temperature, change, iterations) = thermal_endpoint(
            &solid,
            final_bound_moles * MASS,
            heat,
            controls(),
            &mut || false,
        )
        .unwrap();
        near(temperature, target, 2e-10);
        near(
            change[0] + change[1],
            heat,
            controls().gas.energy_tolerance_j,
        );
        assert!(
            iterations <= 5,
            "safeguarded chart solve took {iterations} iterations"
        );
    }
    assert!(matches!(
        thermal_endpoint(
            &solid,
            final_bound_moles * MASS,
            1e8,
            controls(),
            &mut || false
        ),
        Err(AdsorptionError::ThermalRange)
    ));
}
