//! Real IF97/NASA-9 caloric interoperability, distinct from the manufactured
//! phase-law tests. NASA coefficients are read unchanged from the checked-in
//! source packs. The supplied gas transport cards and 0.005 m/s film coefficient
//! are explicit test choices; this is not an experimental drying-rate validation.

use super::*;
use crate::evaporation::WaterEvaporation;
use crate::phase::water::{IsobaricWater, WaterChartControls};
use crate::thermochemical_gas::tests::species_in_domain;
use fs_qty::Pressure;
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::{
    ElementalReferenceIdV1, Nasa9RegionV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1,
    StandardStateConventionV1, StandardStatePhaseV1,
};

fn original_species(name: &str, molar_mass: f64, text: &str) -> Nasa9StandardStateModelV1 {
    let seed = species_in_domain(name, molar_mass, 3.5, 0.0, 0.0, [200.0, 1000.0]);
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
    let source = fs_blake3::hash_bytes(text.as_bytes());
    card.provenance = fs_matdb::Provenance {
        source: "NASA/TP-2002-211556, Appendix D, original 200..1000 K source coefficients".into(),
        license: "Work-of-the-US-Government-Public-Use-Permitted".into(),
        artifact: Some(source),
    };
    card.sources = vec![source];
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

fn original_gas() -> ThermochemicalGasModel {
    let nitrogen = original_species(
        "N2",
        0.028_013_4,
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/nasa9.tsv"
        )),
    );
    let vapor = original_species(
        "H2O",
        0.018_015_28,
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/matdb/seed-v1/water-vapor-nasa9-mcbride2002/nasa9.tsv"
        )),
    );
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![vapor, nitrogen],
            Composition::new(CompositionBasis::MoleFraction, vec![0.01, 0.99]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("H2O", 1.0e-5, 0.02),
            transport("N2", 1.8e-5, 0.026),
        ],
        Some(&SpeciesId::new("H2O").unwrap()),
    )
    .unwrap()
}

fn original_water() -> WaterEvaporation {
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
    WaterEvaporation::try_new(chart, [280.0, 340.0]).unwrap()
}

#[test]
fn original_if97_water_evaporates_into_finite_nasa_gas_with_one_fixed_datum() {
    let source = original_gas();
    let water = original_water();
    let vapor = SpeciesId::new("H2O").unwrap();
    let nitrogen = SpeciesId::new("N2").unwrap();
    let initial_h = water.liquid_state(300.0).unwrap().specific_enthalpy_j_kg();
    let liquid = LiquidSurfaceState {
        mass_kg: 0.02,
        area_m2: 0.01,
        specific_enthalpy_j_kg: initial_h,
    };
    let gas = ClosedGasParcel::new(source.clone(), point(300.0, 100_000.0), 0.01).unwrap();
    // Align ONLY the caloric datum at the starting temperature. The original
    // temperature dependence of neither source is modified to force agreement.
    let binding = GasLiquidBinding {
        species: vapor.clone(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&vapor, 300.0)
            .unwrap()
            - water.vapor_specific_enthalpy_j_kg(300.0).unwrap(),
    };
    let mut solve = controls();
    solve.gas.gas.minimum_temperature_k = 280.0;
    solve.gas.gas.maximum_temperature_k = 340.0;
    solve.gas.gas.energy_tolerance_j_kg = 1e-6;
    solve.gas.energy_tolerance_j = 1e-7;
    solve.maximum_temperature_increment_k = 1.0;
    solve.latent_compatibility_tolerance_j_kg = 5.0;
    let candidate = advance_surface_with_gas(
        &water,
        &binding,
        liquid,
        &gas,
        &source,
        |h| water.temperature_at_specific_enthalpy_k(h),
        0.005,
        0.5,
        solve,
        || false,
    )
    .unwrap();
    let step = candidate.step();
    assert!(step.transferred_mass_kg > 0.0);
    assert!(step.surface_temperature_k[1] < 300.0);
    assert!(step.vapor_partial_pressure_pa[1] > step.vapor_partial_pressure_pa[0]);
    assert_ne!(step.gas_temperature_k[1], step.gas_temperature_k[0]);
    assert!(candidate.gas().state().state().pressure > 100_000.0);
    near(
        candidate.gas().species_mass_kg(&nitrogen).unwrap(),
        gas.species_mass_kg(&nitrogen).unwrap(),
        1e-14,
    );
    near(
        gas.species_mass_kg(&vapor).unwrap() + liquid.mass_kg,
        candidate.gas().species_mass_kg(&vapor).unwrap() + candidate.liquid().mass_kg,
        1e-14,
    );
    near(
        step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
        0.0,
        1e-7,
    );
    near(step.energy_residual_j, 0.0, 1e-7);
    let mismatch = source
        .species_specific_enthalpy_j_kg(&vapor, step.surface_temperature_k[1])
        .unwrap()
        - water
            .vapor_specific_enthalpy_j_kg(step.surface_temperature_k[1])
            .unwrap()
        - binding.liquid_enthalpy_offset_j_kg;
    assert!(
        mismatch.abs() > 1e-4 && mismatch.abs() < solve.latent_compatibility_tolerance_j_kg,
        "original source dependence must remain visible: {mismatch} J/kg"
    );
    let mut tight = solve;
    tight.latent_compatibility_tolerance_j_kg = mismatch.abs() / 2.0;
    assert!(matches!(
        advance_surface_with_gas(
            &water,
            &binding,
            liquid,
            &gas,
            &source,
            |h| water.temperature_at_specific_enthalpy_k(h),
            0.005,
            0.5,
            tight,
            || false
        ),
        Err(GasEvaporationError::CaloricMismatch { .. })
    ));
}

#[test]
fn original_if97_condensation_keeps_upstream_nasa_vapor_sensible_heat() {
    let source = original_gas();
    let water = original_water();
    let vapor = SpeciesId::new("H2O").unwrap();
    let initial_h = water.liquid_state(300.0).unwrap().specific_enthalpy_j_kg();
    let liquid = LiquidSurfaceState {
        mass_kg: 0.02,
        area_m2: 0.01,
        specific_enthalpy_j_kg: initial_h,
    };
    let wet = source
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.05, 0.95]).unwrap(),
        )
        .unwrap();
    let gas = ClosedGasParcel::new(wet, point(310.0, 100_000.0), 0.01).unwrap();
    let binding = GasLiquidBinding {
        species: vapor.clone(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&vapor, 300.0)
            .unwrap()
            - water.vapor_specific_enthalpy_j_kg(300.0).unwrap(),
    };
    let mut solve = controls();
    solve.gas.gas.minimum_temperature_k = 280.0;
    solve.gas.gas.maximum_temperature_k = 340.0;
    solve.gas.gas.energy_tolerance_j_kg = 1e-6;
    solve.gas.energy_tolerance_j = 1e-7;
    solve.maximum_temperature_increment_k = 1.0;
    solve.latent_compatibility_tolerance_j_kg = 5.0;
    let candidate = advance_surface_with_gas(
        &water,
        &binding,
        liquid,
        &gas,
        &source,
        |h| water.temperature_at_specific_enthalpy_k(h),
        0.005,
        0.5,
        solve,
        || false,
    )
    .unwrap();
    let step = candidate.step();
    assert!(step.transferred_mass_kg < 0.0);
    assert!(step.surface_temperature_k[1] > 300.0);
    assert!(step.vapor_partial_pressure_pa[1] < step.vapor_partial_pressure_pa[0]);
    near(
        step.energy_into_gas_j,
        step.transferred_mass_kg
            * source
                .species_specific_enthalpy_j_kg(&vapor, 310.0)
                .unwrap(),
        1e-7,
    );
    let sensible = source
        .species_specific_enthalpy_change_j_kg(&vapor, 300.0, 310.0)
        .unwrap();
    let actual_heat =
        candidate.liquid().mass_kg * (candidate.liquid().specific_enthalpy_j_kg - initial_h);
    assert!(actual_heat > -step.latent_energy_j);
    near(
        actual_heat + step.latent_energy_j,
        -step.transferred_mass_kg * sensible,
        1e-7,
    );
    near(
        step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
        0.0,
        1e-7,
    );
}

#[test]
fn original_water_films_condense_and_evaporate_against_one_nasa_gas_endpoint() {
    use crate::evaporation::gas::chamber::{
        GasLiquidChamberSolve, GasLiquidChamberSurface, advance_chamber_with_gas,
    };
    let source = original_gas();
    let water = original_water();
    let vapor = SpeciesId::new("H2O").unwrap();
    let carrier = SpeciesId::new("N2").unwrap();
    let initial = ClosedGasParcel::new(
        source
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![0.04, 0.96]).unwrap(),
            )
            .unwrap(),
        point(300.0, 100_000.0),
        0.01,
    )
    .unwrap();
    let binding = GasLiquidBinding {
        species: vapor.clone(),
        liquid_enthalpy_offset_j_kg: source
            .species_specific_enthalpy_j_kg(&vapor, 300.0)
            .unwrap()
            - water.vapor_specific_enthalpy_j_kg(300.0).unwrap(),
    };
    let liquid = |temperature| LiquidSurfaceState {
        mass_kg: 0.02,
        area_m2: 0.01,
        specific_enthalpy_j_kg: water
            .liquid_state(temperature)
            .unwrap()
            .specific_enthalpy_j_kg(),
    };
    let inverse = |h| water.temperature_at_specific_enthalpy_k(h);
    let cold = liquid(295.0);
    let warm = liquid(305.0);
    let surfaces = [
        GasLiquidChamberSurface {
            id: 1,
            law: &water,
            binding: &binding,
            liquid: cold,
            temperature_from_enthalpy: &inverse,
            mass_transfer_coefficient_m_s: 0.005,
        },
        GasLiquidChamberSurface {
            id: 2,
            law: &water,
            binding: &binding,
            liquid: warm,
            temperature_from_enthalpy: &inverse,
            mass_transfer_coefficient_m_s: 0.005,
        },
    ];
    let mut phase = controls();
    phase.gas.gas.minimum_temperature_k = 280.0;
    phase.gas.gas.maximum_temperature_k = 340.0;
    phase.gas.gas.energy_tolerance_j_kg = 1e-6;
    phase.gas.energy_tolerance_j = 1e-7;
    phase.maximum_temperature_increment_k = 1.0;
    // One fixed datum; retain the genuine NASA/IF97 vapor-source difference
    // away from 300 K as part of carried energy, without fitting either law.
    phase.latent_compatibility_tolerance_j_kg = 1000.0;
    let candidate = advance_chamber_with_gas(
        &surfaces,
        &initial,
        &source,
        0.1,
        GasLiquidChamberSolve {
            surface: phase,
            maximum_surfaces: 2,
        },
        || false,
    )
    .unwrap();
    let cold_step = &candidate.surfaces()[0];
    let warm_step = &candidate.surfaces()[1];
    assert!(cold_step.transferred_mass_kg < 0.0);
    assert!(warm_step.transferred_mass_kg > 0.0);
    assert!(cold_step.liquid.mass_kg > cold.mass_kg);
    assert!(warm_step.liquid.mass_kg < warm.mass_kg);
    assert!(cold_step.surface_temperature_k[1] > 295.0);
    assert!(warm_step.surface_temperature_k[1] < 305.0);
    near(
        initial.species_mass_kg(&vapor).unwrap() + cold.mass_kg + warm.mass_kg,
        candidate.gas().species_mass_kg(&vapor).unwrap()
            + cold_step.liquid.mass_kg
            + warm_step.liquid.mass_kg,
        phase.gas.mass_tolerance_kg,
    );
    near(
        initial.species_mass_kg(&carrier).unwrap(),
        candidate.gas().species_mass_kg(&carrier).unwrap(),
        phase.gas.mass_tolerance_kg,
    );
    near(
        candidate.step().liquid_enthalpy_change_j + candidate.step().gas_internal_energy_change_j,
        0.0,
        phase.gas.energy_tolerance_j,
    );
    assert!(candidate.step().sum_absolute_energy_residual_j <= phase.gas.energy_tolerance_j);
    assert!(std::ptr::eq(candidate.catalog(), &source));
}
