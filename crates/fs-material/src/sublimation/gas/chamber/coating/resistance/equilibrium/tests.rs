//! Independent phase/energy limits on the genuine original gas fixtures.
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};
use fs_thermal::sublimation::{SolidSurfaceState, SublimationLaw};

const SUPPORT_MASS: f64 = 0.02;
const SUPPORT_CP: f64 = 500.0;

fn controls(datum: f64) -> GasResistiveEquilibriumSolve {
    GasResistiveEquilibriumSolve {
        substrate_enthalpy_bracket_j_kg: [
            datum + SUPPORT_CP * (250.0 - 400.0),
            datum + SUPPORT_CP * (650.0 - 400.0),
        ],
        absolute_energy_tolerance_j: 1e-7,
        absolute_mass_tolerance_kg: 1e-14,
        absolute_temperature_tolerance_k: 1e-7,
        absolute_pressure_tolerance_pa: 1e-7,
        latent_compatibility_tolerance_j_kg: 1e-5,
        maximum_deposit_volume_fraction: 0.01,
        maximum_evaluations: 128,
    }
}

fn state(ts: f64, td: f64, mass: f64, chart: (f64, f64), datum: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: SUPPORT_MASS,
        substrate_specific_enthalpy_j_kg: datum + SUPPORT_CP * (ts - 400.0),
        deposit: SolidSurfaceState {
            mass_kg: mass,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (td - 400.0),
            area_m2: 1e-4,
        },
    }
}

fn surface<'a>(
    state: ResistiveSurfaceState,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    support: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    deposit: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasResistiveSurface<'a> {
    GasResistiveSurface {
        id: 17,
        law,
        binding,
        state,
        conductance_w_k: 1.0,
        substrate_temperature: support,
        deposit_temperature: deposit,
    }
}

// Independent absolute endpoint accounting is meaningful for the unshifted
// fixtures. Datum metamorphisms below compare physical states separately.
fn endpoint_balances(
    initial: &GasResistiveSurface<'_>,
    before: &ClosedGasParcel,
    accepted: &GasResistiveEquilibriumCandidate<'_>,
) {
    let end = accepted.state();
    let initial_solid = initial.state.substrate_mass_kg
        * initial.state.substrate_specific_enthalpy_j_kg
        + initial.state.deposit.mass_kg
            * (initial.state.deposit.specific_enthalpy_j_kg
                + initial.binding.solid_enthalpy_offset_j_kg);
    let final_solid = end.substrate_mass_kg * end.substrate_specific_enthalpy_j_kg
        + end.deposit.mass_kg
            * (end.deposit.specific_enthalpy_j_kg
                + initial.binding.solid_enthalpy_offset_j_kg);
    let gas_change = accepted.gas().mass_kg()
        * accepted.gas().state().specific_internal_energy_j_kg()
        - before.mass_kg() * before.state().specific_internal_energy_j_kg();
    fixture::near(final_solid - initial_solid + gas_change, 0.0, 2e-7);
    fixture::near(
        end.deposit.mass_kg - initial.state.deposit.mass_kg
            + accepted.gas().mass_kg() - before.mass_kg(),
        0.0,
        1e-14,
    );
    fixture::near(accepted.step().gas_energy_change_j, gas_change, 2e-7);
    fixture::near(
        accepted.step().substrate_energy_change_j + accepted.step().deposit_energy_change_j,
        final_solid - initial_solid,
        2e-7,
    );
    fixture::near(accepted.step().energy_residual_j, 0.0, 1e-7);
    assert!(accepted.step().caloric_encoding_residual_j.abs() <= 1e-7);
    assert_eq!(accepted.gas().volume_m3(), before.volume_m3());
    assert_eq!(end.substrate_mass_kg, initial.state.substrate_mass_kg);
    assert_eq!(end.deposit.area_m2, initial.state.deposit.area_m2);
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(
        accepted.gas().species_mass_kg(&carrier).unwrap(),
        before.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
}

#[test]
fn bare_supersaturated_gas_forms_a_positive_film_with_one_conserved_endpoint() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let saved = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let row = surface(state(350.0, 350.0, 0.0, chart, 0.0), &law, &binding, &support, &deposit);
    let result = equilibrate_resistive_chamber(&row, &gas, &catalog, controls(0.0), || false)
        .unwrap();
    let end = result.state();
    let t = result.step().temperature_k;
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / law.molar_mass;
    let psat = law.properties(t).unwrap().saturation_pressure_pa;
    assert!(end.deposit.mass_kg > 0.0);
    assert!(result.step().transferred_mass_kg < 0.0);
    assert_eq!(result.step().phase, GasResistiveEquilibriumPhase::SaturatedCoating);
    fixture::near(result.gas().species_partial_pressure_pa(&fixture::species()).unwrap(), psat, 1e-7);
    fixture::near(
        result.gas().species_mass_kg(&fixture::species()).unwrap(),
        psat * gas.volume_m3() / (r * t),
        1e-14,
    );
    fixture::near(support(end.substrate_specific_enthalpy_j_kg).unwrap(), t, 1e-7);
    fixture::near(deposit(end.deposit.specific_enthalpy_j_kg).unwrap(), t, 1e-7);
    fixture::near(result.gas().state().state().temperature, t, 1e-7);
    assert!(result.step().evaluations <= controls(0.0).maximum_evaluations);
    assert_eq!(result.gas().point().axes().get("local-lot"), Some(&9.0));
    assert_eq!(result.gas().model().water_species(), Some(&fixture::species()));
    endpoint_balances(&row, &gas, &result);
    assert_eq!(gas, saved);
    assert_eq!(row.state.deposit.mass_kg, 0.0);
}

#[test]
fn complete_evaporation_matches_the_independent_constant_cp_energy_solution() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 380.0, 0.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let mass = 1e-5;
    let row = surface(state(420.0, 350.0, mass, chart, 0.0), &law, &binding, &support, &deposit);
    let result = equilibrate_resistive_chamber(&row, &gas, &catalog, controls(0.0), || false)
        .unwrap();
    // The original fixture explicitly has carrier cp=3.5 R and vapor cp=4 R.
    // With all volatile mass in gas, its original u contributes 3 R*T.
    let ra = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let rv = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / law.molar_mass;
    let carrier_capacity = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap() * 2.5 * ra;
    let support_capacity = SUPPORT_MASS * SUPPORT_CP;
    let latent = law.properties(350.0).unwrap().latent_heat_j_kg;
    let expected_t = (support_capacity * 420.0 + carrier_capacity * 380.0
        + mass * (4.0 * rv * 350.0 - latent))
        / (support_capacity + carrier_capacity + 3.0 * rv * mass);
    assert_eq!(result.state().deposit.mass_kg, 0.0);
    assert_eq!(result.step().phase, GasResistiveEquilibriumPhase::AllVapor);
    fixture::near(result.step().temperature_k, expected_t, 2e-8);
    fixture::near(result.step().transferred_mass_kg, mass, 1e-18);
    fixture::near(result.gas().species_mass_kg(&fixture::species()).unwrap(), mass, 1e-14);
    assert!(result.step().saturation_residual_pa < 0.0);
    endpoint_balances(&row, &gas, &result);
}

#[test]
fn existing_film_reaches_saturation_and_exact_stationarity_preserves_inputs() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 500.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let row = surface(state(450.0, 300.0, 0.002, chart, 0.0), &law, &binding, &support, &deposit);
    let result = equilibrate_resistive_chamber(&row, &gas, &catalog, controls(0.0), || false)
        .unwrap();
    assert!(result.state().deposit.mass_kg > 0.0);
    fixture::near(result.step().saturation_residual_pa, 0.0, 1e-7);
    endpoint_balances(&row, &gas, &result);

    let stationary_gas = fixture::gas(&catalog, 400.0, 1000.0);
    let stationary = surface(
        state(400.0, 400.0, 0.002, chart, 0.0),
        &law, &binding, &support, &deposit,
    );
    let unchanged = equilibrate_resistive_chamber(
        &stationary, &stationary_gas, &catalog, controls(0.0), || false,
    ).unwrap();
    assert_eq!(unchanged.state(), stationary.state);
    assert_eq!(unchanged.gas(), &stationary_gas);
    assert_eq!(unchanged.step().energy_residual_j, 0.0);
    assert_eq!(unchanged.step().transferred_mass_kg, 0.0);
    assert_eq!(unchanged.step().evaluations, 0);
}

#[test]
fn independent_support_deposit_and_species_datums_preserve_the_physical_equilibrium() {
    let mut reference: Option<(f64, f64, f64)> = None;
    for (species_reference, deposit_offset, support_datum) in [
        (0.0, 0.0, 0.0),
        (1e7, -3e6, 1e7),
        (-1e7, 2e6, -1e7),
    ] {
        let catalog = fixture::catalog(species_reference);
        let gas = fixture::gas(&catalog, 350.0, 500.0);
        let law = fixture::phase();
        let binding = fixture::binding(deposit_offset);
        let chart = fixture::chart(&catalog, &binding);
        let support = |h| Ok(400.0 + (h - support_datum) / SUPPORT_CP);
        let deposit = |h| fixture::temperature(h, chart);
        let row = surface(
            state(450.0, 300.0, 0.002, chart, support_datum),
            &law, &binding, &support, &deposit,
        );
        let result = equilibrate_resistive_chamber(
            &row, &gas, &catalog, controls(support_datum), || false,
        ).unwrap();
        let actual = (
            result.step().temperature_k,
            result.state().deposit.mass_kg,
            result.gas().species_partial_pressure_pa(&fixture::species()).unwrap(),
        );
        if let Some(expected) = reference {
            fixture::near(actual.0, expected.0, 2e-8);
            fixture::near(actual.1, expected.1, 1e-13);
            fixture::near(actual.2, expected.2, 1e-7);
        } else {
            reference = Some(actual);
        }
        fixture::near(result.step().energy_residual_j, 0.0, 1e-7);
        fixture::near(result.step().mass_residual_kg, 0.0, 1e-14);
    }
}

#[derive(Clone, Copy)]
struct ZeroSaturation(fixture::Phase);
impl SublimationLaw for ZeroSaturation {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.0.temperature_range_k()
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        Ok(SublimationProperties {
            saturation_pressure_pa: 0.0,
            saturation_pressure_derivative_pa_k: 0.0,
            ..self.0.properties(t)?
        })
    }
}

#[test]
fn exact_species_depletion_retains_its_source_but_total_gas_vacuum_refuses() {
    let catalog = fixture::catalog(0.0);
    let law = ZeroSaturation(fixture::phase());
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let row = surface(state(400.0, 400.0, 0.0, chart, 0.0), &law, &binding, &support, &deposit);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let result = equilibrate_resistive_chamber(&row, &gas, &catalog, controls(0.0), || false)
        .unwrap();
    assert_eq!(result.gas().species_mass_kg(&fixture::species()).unwrap(), 0.0);
    assert_eq!(result.gas().species_partial_pressure_pa(&fixture::species()).unwrap(), 0.0);
    fixture::near(
        result.state().deposit.mass_kg,
        gas.species_mass_kg(&fixture::species()).unwrap(),
        1e-14,
    );
    assert!(result.catalog().thermodynamics().components().iter()
        .any(|source| source.species() == &fixture::species()));
    endpoint_balances(&row, &gas, &result);

    let pure = fixture::gas(&catalog, 400.0, 100000.0);
    let saved = pure.clone();
    assert!(matches!(
        equilibrate_resistive_chamber(&row, &pure, &catalog, controls(0.0), || false),
        Err(GasResistiveEquilibriumError::Source(
            GasSurfaceError::Gas(GasEnergyError::Invalid(_))
        ))
    ));
    assert_eq!(pure, saved);
}

#[test]
fn blocked_exchange_bad_sources_volume_and_unbracketed_states_refuse_atomically() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let saved = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = state(350.0, 350.0, 0.0, chart, 0.0);
    let row = surface(initial, &law, &binding, &support, &deposit);
    let insulated = GasResistiveSurface { conductance_w_k: 0.0, ..row };
    let blocked = GasSurfaceBinding { accommodation_coefficient: 0.0, ..binding.clone() };
    let blocked_row = GasResistiveSurface { binding: &blocked, ..row };
    let hidden = GasResistiveSurface {
        state: ResistiveSurfaceState {
            deposit: SolidSurfaceState { area_m2: 0.0, ..initial.deposit },
            ..initial
        },
        ..row
    };
    for refused in [insulated, blocked_row, hidden] {
        assert!(matches!(
            equilibrate_resistive_chamber(&refused, &gas, &catalog, controls(0.0), || false),
            Err(GasResistiveEquilibriumError::Invalid(_))
        ));
    }
    let wrong_inverse = |h| {
        let t = fixture::temperature(h, chart)?;
        Ok(t + 0.001 * (t - 350.0).powi(2))
    };
    let wrong = GasResistiveSurface { deposit_temperature: &wrong_inverse, ..row };
    assert!(matches!(
        equilibrate_resistive_chamber(&wrong, &gas, &catalog, controls(0.0), || false),
        Err(GasResistiveEquilibriumError::Invalid(_))
    ));
    let mut tiny_volume = controls(0.0);
    tiny_volume.maximum_deposit_volume_fraction = 1e-10;
    assert!(matches!(
        equilibrate_resistive_chamber(&row, &gas, &catalog, tiny_volume, || false),
        Err(GasResistiveEquilibriumError::Invalid(_))
    ));
    let mut outside = controls(0.0);
    outside.substrate_enthalpy_bracket_j_kg = [50000.0, 75000.0];
    assert!(matches!(
        equilibrate_resistive_chamber(&row, &gas, &catalog, outside, || false),
        Err(GasResistiveEquilibriumError::OutsideBracket { .. })
    ));
    let mut unsupported = controls(0.0);
    unsupported.substrate_enthalpy_bracket_j_kg[0] = SUPPORT_CP * (150.0 - 400.0);
    assert!(matches!(
        equilibrate_resistive_chamber(&row, &gas, &catalog, unsupported, || false),
        Err(GasResistiveEquilibriumError::Source(GasSurfaceError::Surface(
            SublimationError::TemperatureOutOfRange { .. }
        )))
    ));
    assert_eq!(row.state, initial);
    assert_eq!(gas, saved);
}

#[test]
fn finite_work_and_late_cancellation_never_publish_an_equilibrium_candidate() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let saved = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = state(350.0, 350.0, 0.0, chart, 0.0);
    let row = surface(initial, &law, &binding, &support, &deposit);
    let mut exhausted = controls(0.0);
    exhausted.maximum_evaluations = 2;
    assert!(matches!(
        equilibrate_resistive_chamber(&row, &gas, &catalog, exhausted, || false),
        Err(GasResistiveEquilibriumError::NotConverged { evaluations: 2, .. })
    ));
    let mut polls = 0;
    equilibrate_resistive_chamber(&row, &gas, &catalog, controls(0.0), || {
        polls += 1;
        false
    }).unwrap();
    assert!(polls > 16);
    let cancel_at = polls - 1;
    let mut seen = 0;
    assert!(matches!(
        equilibrate_resistive_chamber(&row, &gas, &catalog, controls(0.0), || {
            seen += 1;
            seen == cancel_at
        }),
        Err(GasResistiveEquilibriumError::Source(GasSurfaceError::Surface(
            SublimationError::Cancelled
        )))
    ));
    assert_eq!(seen, cancel_at);
    assert_eq!(row.state, initial);
    assert_eq!(gas, saved);
}
