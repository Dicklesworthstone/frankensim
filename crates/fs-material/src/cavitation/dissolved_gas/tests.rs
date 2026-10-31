//! Synthetic dilute-solution fixtures exercise the ORIGINAL gas/caloric owners.
//! They verify equations and accounting, not empirical properties of a named pair.
use super::*;
use crate::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::tests::{model, point};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_thermochem::{Composition, CompositionBasis};

const T: f64 = 300.0;
const H: f64 = 2e-5;
const D: f64 = 2e-9;
const MW: f64 = 0.028;
const VL: f64 = 0.002;
const VG: f64 = 0.001;
const SOLUTION_HEAT: f64 = 6000.0;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn claim(name: &str, dims: Dims, value: f64, pressures: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                200.0,
                800.0,
            )
            .with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                pressures[0],
                pressures[1],
            )
            .with_quantity(
                "dissolved_concentration",
                typed(QuantityKind::AmountConcentration),
                0.0,
                100.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic explicit Henry/diffusion/caloric law".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(name: &str, claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for c in claims {
        set.insert_claim(c).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "liquid".into(),
            process: "authored analytic dilute solution".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn context(temperature: f64) -> QueryPoint {
    point(temperature, 100_000.0)
        .with_quantity(
            "dissolved_concentration",
            typed(QuantityKind::AmountConcentration),
            0.0,
        )
        .unwrap()
}
fn calorics(catalog: &ThermochemicalGasModel) -> DissolvedGasCalorics {
    let cp = Dims([2, 0, -2, -1, 0, 0]);
    let rho = Dims([-3, 1, 0, 0, 0, 0]);
    let source = card(
        "dissolved A partial caloric fixture",
        vec![
            claim("partial_heat_capacity", cp, 500.0, [1000.0, 1e7]),
            claim("chart_density", rho, 1000.0, [1000.0, 1e7]),
        ],
    );
    let curve = ResolvedHeatCapacityCurve::try_new(
        &source,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("partial_heat_capacity", cp),
            density: &PropertyKey::new("chart_density", rho),
            temperature_axis: "temperature",
            lower: &context(200.0),
            upper: &context(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-8,
    )
    .unwrap();
    let original_gas_h = catalog.species_specific_enthalpy_j_kg(&id("A"), T).unwrap();
    let hd = curve
        .state_at_temperature(T)
        .unwrap()
        .specific_enthalpy_j_kg();
    DissolvedGasCalorics {
        partial_enthalpy: curve,
        enthalpy_offset_j_kg: original_gas_h - hd - SOLUTION_HEAT / MW,
    }
}
fn source(diffusivity: f64, pressures: [f64; 2]) -> HenrySource {
    HenrySource {
        material: card(
            "A in declared solvent",
            vec![
                claim(
                    HENRY_SOLUBILITY_PROPERTY,
                    HENRY_SOLUBILITY_DIMS,
                    H,
                    pressures,
                ),
                claim(
                    DISSOLVED_DIFFUSIVITY_PROPERTY,
                    DIFFUSIVITY_DIMS,
                    diffusivity,
                    pressures,
                ),
            ],
        ),
        point: context(T),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn law(catalog: &ThermochemicalGasModel) -> HenryDissolution {
    HenryDissolution::new(
        catalog,
        id("A"),
        source(D, [1000.0, 1e7]),
        calorics(catalog),
    )
    .unwrap()
}
fn gas(catalog: &ThermochemicalGasModel) -> ClosedGasParcel {
    ClosedGasParcel::new(catalog.clone(), point(T, 100_000.0), VG).unwrap()
}
fn controls() -> DissolvedGasControls {
    DissolvedGasControls {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 800.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 96,
            },
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_step_s: 1e6,
        maximum_concentration_mol_m3: 20.0,
        maximum_film_relaxation_ratio: 0.01,
    }
}
fn film() -> DissolutionFilm {
    DissolutionFilm {
        area_m2: 0.01,
        thickness_m: 1e-5,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn finite_dissolution_and_exsolution_follow_analytic_partition_and_preserve_inert_gas() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let law = law(&catalog);
    for initial_dissolved in [0.0, 0.006] {
        let mut gas = gas(&catalog);
        let initial_gas = gas.clone();
        let mut liquid = DissolvedGasInventory::new(VL, initial_dissolved, 0.0).unwrap();
        let ng0 = gas.species_mass_kg(&id("A")).unwrap() / MW;
        let inert = gas.species_mass_kg(&id("B")).unwrap();
        let total = ng0 + initial_dissolved;
        // Independent finite-volume equilibrium partition and exponential oracle.
        let equilibrium = total * (H * R * T * VL) / (VG + H * R * T * VL);
        let lambda = D * film().area_m2 / film().thickness_m * (1.0 / VL + H * R * T / VG);
        let duration = 200.0;
        let expected = equilibrium + (initial_dissolved - equilibrium) * (-lambda * duration).exp();
        let step = law
            .advance(&mut gas, &mut liquid, film(), duration, controls(), || {
                false
            })
            .unwrap();
        near(liquid.dissolved_moles(), expected, 2e-15);
        near(
            gas.species_mass_kg(&id("A")).unwrap() / MW + liquid.dissolved_moles(),
            total,
            1e-14,
        );
        near(gas.species_mass_kg(&id("B")).unwrap(), inert, 1e-16);
        assert_eq!(gas.state().state().temperature, T);
        assert_eq!(gas.volume_m3(), VG);
        assert_eq!(
            gas.model().transport_materials(),
            initial_gas.model().transport_materials()
        );
        near(
            gas.species_partial_pressure_pa(&id("A")).unwrap(),
            (total - expected) * R * T / VG,
            1e-7,
        );
        assert_eq!(step.dissolved_moles_change > 0.0, initial_dissolved == 0.0);
        // Native gas heat maintains T while h-u pressure-flow work is counted once.
        near(
            step.thermostat_heat_into_gas_j,
            step.dissolved_moles_change * R * T,
            1e-12,
        );
        near(
            step.thermostat_heat_into_liquid_j,
            -step.dissolved_moles_change * SOLUTION_HEAT,
            1e-11,
        );
        near(
            step.thermostat_heat_into_system_j,
            step.dissolved_moles_change * (R * T - SOLUTION_HEAT),
            1e-11,
        );
        near(
            step.gas.energy_into_gas_j,
            -step.partial_enthalpy_into_liquid_j + step.thermostat_heat_into_gas_j,
            1e-11,
        );
        near(step.energy_residual_j, 0.0, 1e-8);
        near(step.mass_residual_kg, 0.0, 1e-14);
        // A very long admitted interval approaches the finite pair equilibrium,
        // not H times the old reservoir pressure and not negative inventory.
        law.advance(&mut gas, &mut liquid, film(), 1e6, controls(), || false)
            .unwrap();
        near(liquid.dissolved_moles(), equilibrium, 2e-14);
        assert!(gas.species_mass_kg(&id("A")).unwrap() > 0.0);
        assert!(liquid.dissolved_moles() > 0.0);
    }
}

#[test]
fn original_catalog_allows_exsolution_into_a_carrier_without_seed_species() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let carrier = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap())
        .unwrap();
    let mut gas = gas(&carrier);
    let inert = gas.mass_kg();
    let mut liquid = DissolvedGasInventory::new(VL, 0.001, 0.0).unwrap();
    assert_eq!(gas.species_mass_kg(&id("A")).unwrap(), 0.0);
    let step = law(&catalog)
        .advance(&mut gas, &mut liquid, film(), 200.0, controls(), || false)
        .unwrap();
    assert!(gas.species_mass_kg(&id("A")).unwrap() > 0.0);
    assert!(liquid.dissolved_moles() < 0.001);
    near(gas.species_mass_kg(&id("B")).unwrap(), inert, 1e-16);
    near(
        gas.species_mass_kg(&id("A")).unwrap() / MW + liquid.dissolved_moles(),
        0.001,
        1e-14,
    );
    assert_eq!(
        gas.model().thermodynamics().components()[0],
        catalog.thermodynamics().components()[0]
    );
    assert!(step.thermostat_heat_into_system_j > 0.0);
}

#[test]
fn exact_relaxation_keeps_small_steps_partition_consistency_and_repeated_history() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let law = law(&catalog);
    let run = |times: &[f64]| {
        let mut gas = gas(&catalog);
        let mut liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
        let mut reports = Vec::new();
        for &time in times {
            reports.push(
                law.advance(&mut gas, &mut liquid, film(), time, controls(), || false)
                    .unwrap(),
            );
        }
        (gas, liquid, reports)
    };
    let tiny = run(&[1e-8]);
    let rate0 = D * film().area_m2 / film().thickness_m * H * 60_000.0;
    near(tiny.1.dissolved_moles(), rate0 * 1e-8, 1e-23);
    assert!(tiny.1.dissolved_moles() > 0.0);
    let full = run(&[400.0]);
    let split = run(&[200.0, 400.0]);
    near(full.1.dissolved_moles(), split.1.dissolved_moles(), 1e-15);
    near(
        full.0.state().state().pressure,
        split.0.state().state().pressure,
        1e-8,
    );
    assert_eq!(split, run(&[200.0, 400.0]));
}

#[test]
fn late_cancellation_source_limits_and_infinite_dilution_guard_are_atomic() {
    let catalog = model(vec![0.6, 0.4], 0.0);
    let law = law(&catalog);
    let initial_gas = gas(&catalog);
    let initial_liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
    let mut accepted_gas = initial_gas.clone();
    let mut accepted_liquid = initial_liquid;
    let mut total_polls = 0;
    let expected = law
        .advance(
            &mut accepted_gas,
            &mut accepted_liquid,
            film(),
            200.0,
            controls(),
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    let mut gas = initial_gas.clone();
    let mut liquid = initial_liquid;
    let mut polls = 0;
    assert_eq!(
        law.advance(&mut gas, &mut liquid, film(), 200.0, controls(), || {
            polls += 1;
            polls == total_polls
        }),
        Err(DissolvedGasError::Cancelled)
    );
    assert_eq!(gas, initial_gas);
    assert_eq!(liquid, initial_liquid);
    // A dilute bulk cannot authorize an unsupported interfacial Henry concentration.
    let mut dilute = controls();
    dilute.maximum_concentration_mol_m3 = 1.0;
    assert!(matches!(
        law.advance(&mut gas, &mut liquid, film(), 200.0, dilute, || false),
        Err(DissolvedGasError::Concentration { .. })
    ));
    let thick = DissolutionFilm {
        area_m2: 0.01,
        thickness_m: 0.01,
    };
    assert!(matches!(
        law.advance(&mut gas, &mut liquid, thick, 200.0, controls(), || false),
        Err(DissolvedGasError::FilmRelaxation { .. })
    ));
    assert_eq!(gas, initial_gas);
    assert_eq!(liquid, initial_liquid);
    assert_eq!(
        law.advance(&mut gas, &mut liquid, film(), 200.0, controls(), || false)
            .unwrap(),
        expected
    );

    // Source support fails only AFTER native pressure/composition reconstruction.
    let narrow = HenryDissolution::new(
        &catalog,
        id("A"),
        source(D, [90_000.0, 110_000.0]),
        calorics(&catalog),
    )
    .unwrap();
    let mut gas = initial_gas.clone();
    let mut liquid = DissolvedGasInventory::new(VL, 0.02, 0.0).unwrap();
    let before = liquid;
    assert!(matches!(
        narrow.advance(&mut gas, &mut liquid, film(), 2000.0, controls(), || false),
        Err(DissolvedGasError::Material(_))
    ));
    assert_eq!(gas, initial_gas);
    assert_eq!(liquid, before);
}
