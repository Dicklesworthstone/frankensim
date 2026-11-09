//! Original molecular sources and finite-gas droplet integration.
//! This file is included under the existing evaporation fixture module so the
//! NASA/material owners are reused rather than replaced with a new gas model.
#[path = "forced/tests.rs"]
mod forced;

use super::*;
use crate::graph::thermal::medium::mass_transfer::evaporation::droplet::{
    SourceSphericalDropletFilm, SphericalDropletFilmError,
};
use fs_material::evaporation::gas::GasFilmHeatSolve;
use fs_material::evaporation::gas::droplet::{
    SphericalDropletGeometry, SphericalDropletState, advance_droplet_with_gas,
    advance_droplet_with_gas_and_heat,
};
use fs_thermal::evaporation::{InterfaceKinetics, KineticEvaporation};

fn geometry() -> SphericalDropletGeometry {
    SphericalDropletGeometry::from_mass(1e-6, 1000.0).unwrap()
}
fn sphere(caloric: &LiquidCaloric, t: f64, mass_kg: f64) -> SphericalDropletState {
    SphericalDropletState {
        mass_kg,
        specific_enthalpy_j_kg: caloric.h(t),
    }
}
fn inverse(caloric: &LiquidCaloric, h: f64) -> Result<f64, EvaporationError> {
    caloric
        .temperature_k(h)
        .map_err(|_| EvaporationError::Invalid("original caloric fixture refused"))
}
fn heat_controls(caloric: &LiquidCaloric) -> GasFilmHeatSolve {
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [caloric.h(250.0), caloric.h(500.0)],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn samples_refresh_actual_temperature_total_pressure_and_size_with_the_same_source() {
    let catalog = catalog();
    let first_gas = gas(&catalog, 400.0, 100_000.0, 0.0);
    // A spherical Stefan film does not inherit a dilute flat-plate cutoff.
    let changed_gas = gas(&catalog, 450.0, 130_000.0, 0.4);
    let context = point(250.0, 50_000.0);
    for (axis, knots, expected) in [
        ("temperature", vec![(300.0, 1e-5), (500.0, 3e-5)], 2.5e-5),
        (
            "pressure",
            vec![(50_000.0, 3e-5), (150_000.0, 1e-5)],
            1.4e-5,
        ),
    ] {
        let card = material(vec![scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            axis,
            knots,
            [200.0, 900.0],
        )]);
        let binding = source(&card, &context).source;
        let vapor = SpeciesId::new("B").unwrap();
        let prepared =
            SourceSphericalDropletFilm::new(&binding, &vapor, &first_gas, || false).unwrap();
        let pinned = prepared.selected_diffusivity().clone();
        let sample = prepared
            .resolve(geometry(), &changed_gas, || false)
            .unwrap();
        assert!(std::ptr::eq(prepared.diffusion_source().material, &card));
        assert!(std::ptr::eq(prepared.diffusion_source().context, &context));
        assert_eq!(prepared.vapor_species(), &vapor);
        assert_eq!(prepared.carrier_species().as_str(), "A");
        assert_eq!(
            sample.diffusion.species,
            [SpeciesId::new("A").unwrap(), vapor]
        );
        assert_eq!(sample.temperature_k, 450.0);
        assert_eq!(
            sample.total_pressure_pa,
            changed_gas.state().state().pressure
        );
        near(sample.vapor_mole_fraction, 0.4, 1e-15, 0.0);
        assert_eq!(
            sample.diffusion.material.state_coordinate("temperature"),
            Some(450.0)
        );
        assert_eq!(
            sample.diffusion.material.state_coordinate("pressure"),
            Some(sample.total_pressure_pa)
        );
        near(sample.diffusion.diffusivity_m2_s, expected, 1e-18, 1e-13);
        near(
            sample.mass_transfer_coefficient_m_s,
            expected / geometry().radius_m(),
            1e-15,
            1e-13,
        );
        let receipt = &sample
            .diffusion
            .material
            .property(BINARY_DIFFUSIVITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt;
        assert!(!receipt.source_hashes.is_empty());
        assert_eq!(
            receipt.selected,
            card.claims_for(BINARY_DIFFUSIVITY_PROPERTY)[0].0
        );
        assert_eq!(prepared.selected_diffusivity(), &pinned);
        assert_eq!(context, point(250.0, 50_000.0));
        let larger = SphericalDropletGeometry::from_mass(8e-6, 1000.0).unwrap();
        let larger_sample = prepared.resolve(larger, &changed_gas, || false).unwrap();
        near(
            larger_sample.mass_transfer_coefficient_m_s,
            0.5 * sample.mass_transfer_coefficient_m_s,
            1e-15,
            3e-13,
        );
    }
    assert_eq!(
        first_gas
            .species_mass_kg(&SpeciesId::new("B").unwrap())
            .unwrap(),
        0.0
    );
}

#[test]
fn active_pair_pinned_volatile_and_dynamic_source_coordinates_are_admitted_explicitly() {
    let catalog = catalog();
    let original = gas(&catalog, 400.0, 100_000.0, 0.003);
    let saved = original.clone();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let pair = source(&card, &context).source;
    let vapor = SpeciesId::new("B").unwrap();
    assert!(matches!(
        SourceSphericalDropletFilm::new(&pair, &vapor, &original, || true),
        Err(SphericalDropletFilmError::Film(GasFilmError::Cancelled))
    ));
    let mut wrong_pair = pair.clone();
    wrong_pair.species[0] = SpeciesId::new("C").unwrap();
    assert!(matches!(
        SourceSphericalDropletFilm::new(&wrong_pair, &vapor, &original, || false),
        Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(_)))
    ));
    wrong_pair.species[0] = vapor.clone();
    assert!(SourceSphericalDropletFilm::new(&wrong_pair, &vapor, &original, || false).is_err());
    assert!(
        SourceSphericalDropletFilm::new(&pair, &SpeciesId::new("C").unwrap(), &original, || false)
            .is_err()
    );
    let pure_vapor = gas(&catalog, 400.0, 100_000.0, 1.0);
    assert!(SourceSphericalDropletFilm::new(&pair, &vapor, &pure_vapor, || false).is_err());
    for invalid_context in [
        context.clone().with("humidity", 0.003).unwrap(),
        context.clone().with("T", 400.0).unwrap(),
        QueryPoint::new()
            .with_quantity("temperature", typed(QuantityKind::Pressure), 100_000.0)
            .unwrap(),
    ] {
        let invalid = BinaryDiffusionMaterial {
            context: &invalid_context,
            ..pair.clone()
        };
        assert!(SourceSphericalDropletFilm::new(&invalid, &vapor, &original, || false).is_err());
    }
    for value in [0.0, -1e-5] {
        let invalid_card = material(vec![scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(200.0, value), (900.0, value)],
            [200.0, 900.0],
        )]);
        assert!(
            SourceSphericalDropletFilm::new(
                &source(&invalid_card, &context).source,
                &vapor,
                &original,
                || false
            )
            .is_err()
        );
    }
    let prepared = SourceSphericalDropletFilm::new(&pair, &vapor, &original, || false).unwrap();
    let caloric = LiquidCaloric::new(&catalog);
    let wrong_binding = GasLiquidBinding {
        species: SpeciesId::new("A").unwrap(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let state = sphere(&caloric, 400.0, 1e-6);
    let no_caloric_call = |_| -> Result<f64, EvaporationError> {
        panic!("wrong volatile identity must refuse before the original caloric owner")
    };
    assert!(matches!(
        prepared.advance_with_gas(
            &LiquidPhase,
            &wrong_binding,
            state,
            &original,
            &catalog,
            no_caloric_call,
            0.01,
            phase_controls().phase,
            || false
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(
            "spherical film cannot change its pinned volatile species"
        )))
    ));
    assert!(matches!(
        prepared.advance_with_gas_and_heat(
            &LiquidPhase,
            &wrong_binding,
            state,
            &original,
            &catalog,
            no_caloric_call,
            80.0,
            0.01,
            phase_controls().phase,
            heat_controls(&caloric),
            || false
        ),
        Err(SphericalDropletFilmError::Film(GasFilmError::Invalid(
            "spherical film cannot change its pinned volatile species"
        )))
    ));
    assert_eq!(original, saved);
}

#[test]
fn sourced_phase_matches_original_owner_in_both_directions_and_refreshes_the_next_interval() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let card = diffusion();
    let context = point(250.0, 50_000.0);
    let selected = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 1e-5,
        },
    )
    .unwrap();
    for (t, vapor) in [(400.0, 0.0), (430.0, 0.02)] {
        let original = gas(&catalog, t, 100_000.0, vapor);
        let saved = original.clone();
        let prepared = SourceSphericalDropletFilm::new(
            &source(&card, &context).source,
            &binding.species,
            &original,
            || false,
        )
        .unwrap();
        let state = sphere(&caloric, 400.0, 1e-6);
        let next = prepared
            .advance_with_gas(
                &selected,
                &binding,
                state,
                &original,
                &catalog,
                |h| inverse(&caloric, h),
                0.25,
                phase_controls().phase,
                || false,
            )
            .unwrap();
        let d = 1e-5 + (t - 200.0) * (2e-5 / 700.0);
        let direct = advance_droplet_with_gas(
            &selected,
            &binding,
            state,
            &original,
            &catalog,
            |h| inverse(&caloric, h),
            next.films()[0].diffusion.diffusivity_m2_s,
            0.25,
            phase_controls().phase,
            || false,
        )
        .unwrap();
        assert_eq!(next.original(), &direct);
        near(next.films()[0].diffusion.diffusivity_m2_s, d, 1e-18, 1e-13);
        assert_eq!(next.films()[0].geometry, direct.geometry()[0]);
        assert_eq!(next.films()[1].geometry, direct.geometry()[1]);
        assert_eq!(
            next.films()[1].temperature_k,
            direct.gas().state().state().temperature
        );
        assert_eq!(
            next.films()[1].total_pressure_pa,
            direct.gas().state().state().pressure
        );
        assert_eq!(direct.step().transferred_mass_kg > 0.0, vapor == 0.0);
        near(
            direct.gas().species_mass_kg(&binding.species).unwrap() + direct.droplet().mass_kg,
            original.species_mass_kg(&binding.species).unwrap() + state.mass_kg,
            1e-12,
            0.0,
        );
        let carrier = SpeciesId::new("A").unwrap();
        near(
            direct.gas().species_mass_kg(&carrier).unwrap(),
            original.species_mass_kg(&carrier).unwrap(),
            1e-12,
            0.0,
        );
        assert!(
            direct.step().energy_residual_j.abs() <= phase_controls().phase.gas.energy_tolerance_j
        );
        assert_eq!(direct.gas().volume_m3(), original.volume_m3());
        assert_eq!(direct.gas().point().axes().get("local-lot"), Some(&9.0));
        assert!(std::ptr::eq(direct.catalog(), &catalog));
        let second = prepared
            .advance_with_gas(
                &selected,
                &binding,
                direct.droplet(),
                direct.gas(),
                &catalog,
                |h| inverse(&caloric, h),
                0.25,
                phase_controls().phase,
                || false,
            )
            .unwrap();
        assert_eq!(second.films()[0], next.films()[1]);
        assert_eq!(original, saved);
    }
}

#[test]
fn original_reciprocal_heat_uses_one_frozen_source_and_reports_the_actual_final_sample() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let original = gas(&catalog, 430.0, 100_000.0, 0.003);
    let saved = original.clone();
    let card = diffusion();
    let context = point(400.0, 100_000.0);
    let prepared = SourceSphericalDropletFilm::new(
        &source(&card, &context).source,
        &binding.species,
        &original,
        || false,
    )
    .unwrap();
    let state = sphere(&caloric, 400.0, 1e-6);
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let next = prepared
        .advance_with_gas_and_heat(
            &sealed,
            &binding,
            state,
            &original,
            &catalog,
            |h| inverse(&caloric, h),
            80.0,
            0.1,
            phase_controls().phase,
            heat_controls(&caloric),
            || false,
        )
        .unwrap();
    let direct = advance_droplet_with_gas_and_heat(
        &sealed,
        &binding,
        state,
        &original,
        &catalog,
        |h| inverse(&caloric, h),
        next.films()[0].diffusion.diffusivity_m2_s,
        80.0,
        0.1,
        phase_controls().phase,
        heat_controls(&caloric),
        || false,
    )
    .unwrap();
    assert_eq!(next.original(), &direct);
    assert!(direct.heat().unwrap().heat_into_liquid_j > 0.0);
    assert_eq!(direct.droplet().mass_kg, state.mass_kg);
    assert_eq!(direct.gas().mass_kg(), original.mass_kg());
    assert_eq!(direct.step().transferred_mass_kg, 0.0);
    assert_ne!(
        next.films()[0].diffusion.diffusivity_m2_s,
        next.films()[1].diffusion.diffusivity_m2_s
    );
    assert_eq!(
        direct.report().binary_diffusivity_m2_s,
        next.films()[0].diffusion.diffusivity_m2_s
    );
    near(
        direct.step().liquid_enthalpy_change_j + direct.step().gas_internal_energy_change_j,
        0.0,
        phase_controls().phase.gas.energy_tolerance_j,
        0.0,
    );
    assert_eq!(original, saved);
}

#[test]
fn a_pinned_claim_cannot_switch_and_late_source_exhaustion_discards_the_complete_candidate() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let binding = liquid_binding();
    let original = gas(&catalog, 400.0, 100_000.0, 0.003);
    let saved = original.clone();
    let card = material(vec![
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(200.0, 1e-5), (400.001, 1e-5)],
            [200.0, 400.001],
        ),
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(400.002, 2e-5), (900.0, 2e-5)],
            [400.002, 900.0],
        ),
    ]);
    let context = point(400.0, 100_000.0);
    let pair = source(&card, &context).source;
    let prepared =
        SourceSphericalDropletFilm::new(&pair, &binding.species, &original, || false).unwrap();
    let pinned = prepared.selected_diffusivity().clone();
    let changed = gas(&catalog, 450.0, 100_000.0, 0.003);
    assert!(matches!(
        prepared.resolve(geometry(), &changed, || false),
        Err(SphericalDropletFilmError::Film(GasFilmError::Source(_)))
    ));
    assert!(SourceSphericalDropletFilm::new(&pair, &binding.species, &changed, || false).is_ok());
    let state = sphere(&caloric, 450.0, 1e-4);
    let sealed = KineticEvaporation::try_new(
        &LiquidPhase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 0.0,
        },
    )
    .unwrap();
    let direct = advance_droplet_with_gas_and_heat(
        &sealed,
        &binding,
        state,
        &original,
        &catalog,
        |h| inverse(&caloric, h),
        1e-5,
        80.0,
        0.1,
        phase_controls().phase,
        heat_controls(&caloric),
        || false,
    )
    .unwrap();
    assert!(direct.gas().state().state().temperature > 400.002);
    let result = prepared.advance_with_gas_and_heat(
        &sealed,
        &binding,
        state,
        &original,
        &catalog,
        |h| inverse(&caloric, h),
        80.0,
        0.1,
        phase_controls().phase,
        heat_controls(&caloric),
        || false,
    );
    assert!(matches!(
        result,
        Err(SphericalDropletFilmError::Film(GasFilmError::Source(_)))
    ));
    assert_eq!(prepared.selected_diffusivity(), &pinned);
    assert_eq!(original, saved);
    assert_eq!(state, sphere(&caloric, 450.0, 1e-4));
}

fn nasa_species(name: &str, molar_mass: f64, text: &str) -> Nasa9StandardStateModelV1 {
    let seed = species(name, molar_mass, 3.5);
    let mut card = seed.regions()[0].card().clone();
    let mut count = 0;
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.first() == Some(&"coefficient") && fields[2] == "low" {
            assert_eq!(fields[1], name);
            card.parameters.get_mut(fields[3]).unwrap().value = fields[4].parse().unwrap();
            count += 1;
        }
    }
    assert_eq!(count, 9);
    let hash = fs_blake3::hash_bytes(text.as_bytes());
    card.sources = vec![hash];
    card.provenance = Provenance {
        source: "NASA/TP-2002-211556, Appendix D, unchanged low-temperature coefficients".into(),
        license: "Work-of-the-US-Government-Public-Use-Permitted".into(),
        artifact: Some(hash),
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(100_000.0),
        ElementalReferenceIdV1::new("NASA-TP-2002-211556-reference-elements-298.15K-1bar").unwrap(),
    )
    .unwrap();
    Nasa9StandardStateModelV1::new(
        seed.species().clone(),
        seed.molar_mass(),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

#[test]
fn original_if97_water_and_nasa_gas_keep_species_energy_and_source_identity() {
    use fs_material::evaporation::WaterEvaporation;
    use fs_material::phase::water::{IsobaricWater, WaterChartControls};
    use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

    let nitrogen = nasa_species(
        "N2",
        0.028_013_4,
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/nasa9.tsv"
        )),
    );
    let vapor = nasa_species(
        "H2O",
        0.018_015_28,
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/matdb/seed-v1/water-vapor-nasa9-mcbride2002/nasa9.tsv"
        )),
    );
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![nitrogen, vapor],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("N2", 1.8e-5, 0.026),
            transport("H2O", 1.0e-5, 0.02),
        ],
        Some(&SpeciesId::new("H2O").unwrap()),
    )
    .unwrap();
    let chart = IsobaricWater::new(
        Pressure::new(100_000.0),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 450.0,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65535,
            },
        },
        || false,
    )
    .unwrap();
    let water = WaterEvaporation::try_new(chart, [280.0, 340.0]).unwrap();
    let vapor = SpeciesId::new("H2O").unwrap();
    let carrier = SpeciesId::new("N2").unwrap();
    let binding = GasLiquidBinding {
        species: vapor.clone(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&vapor, 300.0)
            .unwrap()
            - water.vapor_specific_enthalpy_j_kg(300.0).unwrap(),
    };
    let state = SphericalDropletState {
        mass_kg: 1e-6,
        specific_enthalpy_j_kg: water.liquid_state(300.0).unwrap().specific_enthalpy_j_kg(),
    };
    // Explicit manufactured D(T) and accommodation exercise source plumbing;
    // these are not experimental water-evaporation validation data.
    let card = material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY,
        BINARY_DIFFUSIVITY_DIMS,
        "temperature",
        vec![(280.0, 1.8e-5), (340.0, 2.5e-5)],
        [280.0, 340.0],
    )]);
    let context = point(300.0, 100_000.0);
    let pair = BinaryDiffusionMaterial {
        species: [carrier.clone(), vapor.clone()],
        material: &card,
        context: &context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let selected = KineticEvaporation::try_new(
        &water,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 1e-5,
        },
    )
    .unwrap();
    let mut solve = phase_controls().phase;
    solve.gas.gas.minimum_temperature_k = 280.0;
    solve.gas.gas.maximum_temperature_k = 340.0;
    solve.gas.gas.energy_tolerance_j_kg = 1e-6;
    solve.gas.energy_tolerance_j = 1e-7;
    solve.maximum_temperature_increment_k = 1.0;
    solve.latent_compatibility_tolerance_j_kg = 5.0;
    for (t, fraction) in [(300.0, 0.01), (310.0, 0.05)] {
        // Composition coordinates follow the original catalog's canonical
        // species order (H2O sorts before N2), never a guessed vapor index.
        let fractions: Vec<_> = source.thermodynamics().components().iter()
            .map(|component| if component.species() == &vapor { fraction } else { 1.0 - fraction })
            .collect();
        let model = source.with_composition(
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        ).unwrap();
        let original = ClosedGasParcel::new(
            model, point(t, 100_000.0).with("local-lot", 9.0).unwrap(), 0.01,
        ).unwrap();
        let saved = original.clone();
        let prepared = SourceSphericalDropletFilm::new(&pair, &vapor, &original, || false).unwrap();
        let next = prepared
            .advance_with_gas(
                &selected,
                &binding,
                state,
                &original,
                &source,
                |h| water.temperature_at_specific_enthalpy_k(h),
                0.01,
                solve,
                || false,
            )
            .unwrap();
        let candidate = next.original();
        assert_eq!(candidate.step().transferred_mass_kg > 0.0, fraction < 0.02);
        assert_eq!(
            candidate.geometry()[1].radius_m() < candidate.geometry()[0].radius_m(),
            fraction < 0.02
        );
        near(
            candidate.gas().species_mass_kg(&vapor).unwrap() + candidate.droplet().mass_kg,
            original.species_mass_kg(&vapor).unwrap() + state.mass_kg,
            solve.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            candidate.gas().species_mass_kg(&carrier).unwrap(),
            original.species_mass_kg(&carrier).unwrap(),
            solve.gas.mass_tolerance_kg,
            0.0,
        );
        near(
            candidate.step().liquid_enthalpy_change_j
                + candidate.step().gas_internal_energy_change_j,
            0.0,
            solve.gas.energy_tolerance_j,
            0.0,
        );
        assert!(candidate.step().energy_residual_j.abs() <= solve.gas.energy_tolerance_j);
        assert_eq!(candidate.gas().volume_m3(), original.volume_m3());
        assert!(std::ptr::eq(candidate.catalog(), &source));
        assert!(std::ptr::eq(selected.source(), &water));
        assert_eq!(
            next.films()[1]
                .diffusion
                .material
                .state_coordinate("temperature"),
            Some(candidate.gas().state().state().temperature)
        );
        assert_eq!(original, saved);
    }
}
