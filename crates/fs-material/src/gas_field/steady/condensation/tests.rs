//! Manufactured gas sources and geometry; IAPWS supplies the water phase curve.
use super::*;
use crate::gas_field::AxialGasSampling;
use crate::gas_field::steady::{AxialGasCell, GasInjection};
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::tests::{species_with_energy_reference, transport};
use fs_qty::Temperature;
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, mixture::FrozenIdealGasMixtureModelV1,
    water::SaturatedWaterState,
};

const FLOW: f64 = 5e-5;
const PRESSURE: f64 = 101325.0;
const LENGTH: f64 = 0.4;

fn saturated_fraction(t: f64) -> f64 {
    SaturatedWaterState::at_temperature(Temperature::new(t))
        .unwrap()
        .pressure()
        .value()
        / PRESSURE
}
fn primitive(t: f64, water: f64) -> PrimitiveGasState {
    PrimitiveGasState {
        temperature_k: t,
        composition: Composition::new(CompositionBasis::MoleFraction, vec![1.0 - water, water])
            .unwrap(),
    }
}
fn source(shift: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("carrier", 0.028, 3.5, 0.0, 0.0),
                species_with_energy_reference("water", 0.01801528, 4.0, 0.001, shift),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.95, 0.05]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("carrier", 1.7e-5, 0.025),
            transport("water", 9e-6, 0.019),
        ],
        Some(&SpeciesId::new("water").unwrap()),
    )
    .unwrap()
}
fn fixture(count: usize) -> SteadyAxialGas {
    let dx = LENGTH / count as f64;
    SteadyAxialGas {
        source: source(0.0),
        source_conditions: QueryPoint::new(),
        static_pressure_pa: PRESSURE,
        phase: ThermochemicalPhase::LiquidWaterIapwsLimit,
        inlet: primitive(310.0, saturated_fraction(310.0)),
        inlet_mass_flow_kg_s: FLOW,
        cells: (0..count)
            .map(|index| AxialGasCell {
                end_m: (index + 1) as f64 * dx,
                volume_m3: 1e-4 * dx,
                wall_area_m2: 0.04 * dx,
                wall_conductance_w_k: 0.08 / count as f64,
                wall_temperature_k: 285.0,
                initial_temperature_k: 300.0,
            })
            .collect(),
        injections: vec![],
    }
}
fn model(problem: &SteadyAxialGas) -> WaterCondensationModel {
    WaterCondensationModel::new(
        problem.source.clone(),
        SpeciesId::new("water").unwrap(),
        300.0,
    )
    .unwrap()
}
fn controls() -> CondensingSteadyGasSolve {
    CondensingSteadyGasSolve {
        caloric: GasEnergySolve {
            minimum_temperature_k: 274.0,
            maximum_temperature_k: 323.0,
            energy_tolerance_j_kg: 1e-9,
            maximum_iterations: 100,
        },
        maximum_temperature_trials_per_cell: 100,
        steady_power_tolerance_w: 1e-11,
        mass_flow_tolerance_kg_s: 1e-18,
        maximum_mean_mach: 0.05,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}
fn water_mass_fraction(x: f64) -> f64 {
    x * 0.01801528 / ((1.0 - x) * 0.028 + x * 0.01801528)
}

#[test]
fn g0_condensation_sink_ledger_reaches_every_cell_and_all_field_samples() {
    let problem = fixture(4);
    let phase = model(&problem);
    let solution = problem
        .solve_condensing(&phase, controls(), || false)
        .unwrap();
    assert_eq!(
        solution.field().sampling(),
        AxialGasSampling::CellwiseConstant
    );
    let mut previous_temperature = problem.inlet.temperature_k;
    let mut previous_flow = FLOW;
    let mut wall_w = 0.0;
    let mut latent_w = 0.0;
    for (index, cell) in solution.cells.iter().enumerate() {
        assert!(cell.temperature_k < previous_temperature);
        assert!(cell.condensed_mass_flow_kg_s > 0.0);
        assert!(cell.gas_mass_flow_kg_s < previous_flow);
        near(
            cell.gas_mass_flow_kg_s + cell.condensed_mass_flow_kg_s,
            previous_flow,
            1e-18,
        );
        assert!(cell.power_residual_bound_w <= controls().steady_power_tolerance_w);
        assert!(cell.source_relative_power_residual_w.abs() <= controls().steady_power_tolerance_w);
        assert!(cell.maximum_species_mass_residual_kg_s <= controls().mass_flow_tolerance_kg_s);
        assert!(cell.temperature_trials <= controls().maximum_temperature_trials_per_cell);
        let state = solution
            .field()
            .state_at(problem.cells[index].end_m)
            .unwrap();
        near(
            state.state().water_mole_fraction,
            saturated_fraction(cell.temperature_k),
            1e-14,
        );
        state.state().validate_coherence().unwrap();
        previous_temperature = cell.temperature_k;
        previous_flow = cell.gas_mass_flow_kg_s;
        wall_w += cell.wall_heat_into_gas_w;
        latent_w += cell.latent_heat_w;
    }
    let outlet = solution.field().state_at(LENGTH).unwrap();
    let inlet = solution.field().state_at(0.0).unwrap();
    let original_water = FLOW * water_mass_fraction(inlet.state().water_mole_fraction);
    let outgoing_water = previous_flow * water_mass_fraction(outlet.state().water_mole_fraction);
    near(
        outgoing_water + solution.total_condensed_mass_flow_kg_s,
        original_water,
        1e-12 * original_water,
    );
    near(
        previous_flow + solution.total_condensed_mass_flow_kg_s,
        FLOW,
        1e-12 * FLOW,
    );
    let incoming_h_w = FLOW * inlet.specific_enthalpy_j_kg();
    let outgoing_h_w =
        previous_flow * outlet.specific_enthalpy_j_kg() + solution.total_liquid_enthalpy_flow_w;
    near(
        incoming_h_w + wall_w,
        outgoing_h_w,
        1e-12 * incoming_h_w.abs(),
    );
    assert!(
        latent_w > 0.25 * wall_w.abs(),
        "latent heat must be materially live"
    );
    assert!(solution.maximum_mean_mach <= controls().maximum_mean_mach);
    for sample in 0..=100 {
        let gas = solution
            .field()
            .state_at(LENGTH * sample as f64 / 100.0)
            .unwrap();
        let saturation = saturated_fraction(gas.state().temperature);
        assert!(gas.state().water_mole_fraction <= saturation * (1.0 + 1e-13));
    }
    assert_eq!(
        solution,
        problem
            .solve_condensing(&phase, controls(), || false)
            .unwrap()
    );
}

#[test]
fn g0_dry_limit_matches_the_original_upwind_constant_cp_balance() {
    let mut problem = fixture(4);
    problem.inlet = primitive(310.0, 0.0);
    let solution = problem
        .solve_condensing(&model(&problem), controls(), || false)
        .unwrap();
    let cp = 3.5 * R / 0.028;
    let mut expected = 310.0;
    for (index, cell) in solution.cells.iter().enumerate() {
        let conductance = problem.cells[index].wall_conductance_w_k;
        expected = (FLOW * cp * expected + conductance * 285.0) / (FLOW * cp + conductance);
        near(cell.temperature_k, expected, 1e-7);
        assert_eq!(cell.condensed_mass_flow_kg_s, 0.0);
        assert_eq!(cell.gas_mass_flow_kg_s, FLOW);
        assert_eq!(cell.latent_heat_w, 0.0);
    }
    assert_eq!(solution.total_condensed_mass_flow_kg_s, 0.0);
}

#[test]
fn g0_dry_side_injection_keeps_its_enthalpy_and_cannot_reuse_drained_water() {
    let mut problem = fixture(4);
    problem.injections.push(GasInjection {
        cell: 2,
        mass_flow_kg_s: FLOW,
        state: primitive(315.0, 0.0),
    });
    for cell in &mut problem.cells[2..] {
        cell.wall_conductance_w_k = 0.0;
    }
    let solution = problem
        .solve_condensing(&model(&problem), controls(), || false)
        .unwrap();
    assert!(solution.cells[0].condensed_mass_flow_kg_s > 0.0);
    assert!(solution.cells[1].condensed_mass_flow_kg_s > 0.0);
    assert_eq!(solution.cells[2].condensed_mass_flow_kg_s, 0.0);
    assert_eq!(solution.cells[3].condensed_mass_flow_kg_s, 0.0);
    assert!(solution.cells[2].temperature_k > solution.cells[1].temperature_k);
    near(
        solution.cells[2].gas_mass_flow_kg_s,
        solution.cells[1].gas_mass_flow_kg_s + FLOW,
        1e-18,
    );
    let outlet = solution.field().state_at(LENGTH).unwrap();
    assert!(outlet.state().water_mole_fraction < saturated_fraction(outlet.state().temperature));
    let inlet = solution.field().state_at(0.0).unwrap();
    let out_water = solution.cells[3].gas_mass_flow_kg_s
        * water_mass_fraction(outlet.state().water_mole_fraction);
    near(
        out_water + solution.total_condensed_mass_flow_kg_s,
        FLOW * water_mass_fraction(inlet.state().water_mole_fraction),
        1e-18,
    );
    near(
        solution.cells[3].gas_mass_flow_kg_s + solution.total_condensed_mass_flow_kg_s,
        2.0 * FLOW,
        1e-18,
    );
    let point = QueryPoint::new()
        .with_quantity("pressure", typed(QuantityKind::Pressure), PRESSURE)
        .unwrap();
    let injected = resolve_primitive(
        &problem.source,
        &point,
        problem.phase,
        &problem.injections[0].state,
    )
    .unwrap();
    let wall_w: f64 = solution
        .cells
        .iter()
        .map(|cell| cell.wall_heat_into_gas_w)
        .sum();
    near(
        FLOW * (inlet.specific_enthalpy_j_kg() + injected.specific_enthalpy_j_kg()) + wall_w,
        solution.cells[3].gas_mass_flow_kg_s * outlet.specific_enthalpy_j_kg()
            + solution.total_liquid_enthalpy_flow_w,
        1e-10,
    );
}

#[test]
fn g1_cell_refinement_reduces_the_remaining_spatial_change() {
    let mut temperatures = Vec::new();
    let mut liquid_rates = Vec::new();
    for count in [4, 8, 16] {
        let problem = fixture(count);
        let solved = problem
            .solve_condensing(&model(&problem), controls(), || false)
            .unwrap();
        temperatures.push(solved.cells.last().unwrap().temperature_k);
        liquid_rates.push(solved.total_condensed_mass_flow_kg_s);
        assert!(solved.maximum_power_residual_bound_w() <= controls().steady_power_tolerance_w);
    }
    let coarse = (temperatures[1] - temperatures[0]).abs();
    let fine = (temperatures[2] - temperatures[1]).abs();
    assert!(coarse > 1e-3, "the cell resolution axis must be live");
    assert!(fine < 0.7 * coarse, "T={temperatures:?}");
    assert!(
        (liquid_rates[2] - liquid_rates[1]).abs() < 0.7 * (liquid_rates[1] - liquid_rates[0]).abs()
    );
}

#[test]
fn g4_work_phase_bracket_cancellation_and_mean_flow_refuse_without_publication() {
    let original = fixture(2);
    let before = original.clone();
    let model = model(&original);
    assert!(matches!(
        original.solve_condensing(&model, controls(), || true),
        Err(CondensingSteadyGasError::Steady(SteadyGasError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    let mut checks = 0;
    assert!(
        original
            .solve_condensing(&model, controls(), || {
                checks += 1;
                checks > 15
            })
            .is_err()
    );
    assert!(matches!(
        original.solve_condensing(
            &model,
            CondensingSteadyGasSolve {
                maximum_temperature_trials_per_cell: 2,
                ..controls()
            },
            || false
        ),
        Err(CondensingSteadyGasError::NotConverged { .. })
    ));
    assert!(matches!(
        original.solve_condensing(
            &model,
            CondensingSteadyGasSolve {
                caloric: GasEnergySolve {
                    minimum_temperature_k: 310.0,
                    ..controls().caloric
                },
                ..controls()
            },
            || false
        ),
        Err(CondensingSteadyGasError::TemperatureBracket { .. })
    ));
    let mut wrong_phase = original.clone();
    wrong_phase.phase = ThermochemicalPhase::FrozenGas;
    assert!(
        wrong_phase
            .solve_condensing(&model, controls(), || false)
            .is_err()
    );
    let mut fast = original.clone();
    fast.inlet_mass_flow_kg_s = 1.0;
    assert!(matches!(
        fast.solve_condensing(&model, controls(), || false),
        Err(CondensingSteadyGasError::Steady(
            SteadyGasError::MeanFlow { .. }
        ))
    ));
    assert_eq!(original, before);
}

#[test]
fn g4_unresolved_formation_datum_cancellation_is_not_a_false_steady_solution() {
    let mut problem = fixture(1);
    problem.source = source(1e20);
    let before = problem.clone();
    // At this datum a source-relative watt budget cannot be inferred from the
    // rounded gross transport powers. The original owner or audit must refuse.
    assert!(
        problem
            .solve_condensing(&model(&problem), controls(), || false)
            .is_err()
    );
    assert_eq!(problem, before);
}
