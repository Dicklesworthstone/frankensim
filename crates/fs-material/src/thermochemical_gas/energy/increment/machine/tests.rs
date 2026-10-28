use super::*;
use crate::thermochemical_gas::{ThermochemicalGasModel, tests::{model, model_in_domains,
    point, close, species_with_energy_reference, transport}};
use crate::thermochemical_gas::energy::GasEnergySolve;
use fs_thermochem::{Composition, CompositionBasis};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

pub(super) fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
pub(super) fn shifted(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", 0.028, 3.5, 0.002, reference),
            species_with_energy_reference("B", 0.004, 2.5, 0.0, -reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    ).unwrap(), vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap()
}

#[test]
fn compressor_and_turbine_match_independent_constant_capacity_work_and_entropy() {
    let source = model(vec![0.75, 0.25], 0.0);
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    for (mode, t0, p1) in [(GasMachineMode::Compressor, 400.0, 400000.0),
        (GasMachineMode::Turbine, 1000.0, 25000.0)] {
        for eta in [0.5, 0.8, 1.0] {
            let context = point(t0, 100000.0).with("material-lot", 42.0).unwrap();
            let result = branch.adiabatic_machine(&context, Pressure::new(p1), mode, eta, 1e-8, || false).unwrap();
            let gas = result.inlet().state();
            let ts = t0 * (p1 / 100000.0_f64).powf((gas.gamma - 1.0) / gas.gamma);
            let factor = if mode == GasMachineMode::Compressor { 1.0 / eta } else { eta };
            let expected = t0 + factor * (ts - t0);
            close(result.isentropic_temperature_k, ts);
            close(result.outlet().state().temperature, expected);
            close(result.shaft_work_by_gas_j_kg, -gas.specific_heat_cp * factor * (ts - t0));
            assert_eq!(result.outlet().state().pressure, p1);
            assert_eq!(result.outlet_point().axes().get("material-lot"), Some(&42.0));
            assert!(result.energy_residual_j_kg.abs() <= controls().energy_tolerance_j_kg);
            assert!(result.entropy_change_j_kg_k >= -1e-8);
            if eta < 1.0 { assert!(result.entropy_change_j_kg_k > 0.0); }
            let rates = result.rates(2.0).unwrap();
            assert_eq!(rates.shaft_power_by_gas_w, 2.0 * result.shaft_work_by_gas_j_kg);
            assert_eq!(result.rates(0.0).unwrap().shaft_power_by_gas_w, 0.0);
            assert!(result.rates(-1.0).is_err());
        }
    }
}

#[test]
fn variable_capacity_stages_and_tangents_survive_large_formation_references() {
    let ordinary = shifted(0.0);
    let context = point(500.0, 100000.0);
    let expected = ordinary.prepare_caloric_branch(controls(), || false).unwrap()
        .adiabatic_machine(&context, Pressure::new(400000.0), GasMachineMode::Compressor, 0.8, 1e-8, || false).unwrap();
    for reference in [-1e24, 1e24] {
        let source = shifted(reference);
        let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
        let result = branch.adiabatic_machine(&context, Pressure::new(400000.0), GasMachineMode::Compressor, 0.8, 1e-8, || false).unwrap();
        assert_eq!(result.outlet().state(), expected.outlet().state());
        assert_eq!(result.shaft_work_by_gas_j_kg, expected.shaft_work_by_gas_j_kg);
        assert_eq!(result.derivatives, expected.derivatives);
        assert!(result.outlet().specific_enthalpy_j_kg().abs() > 1e20);
        let sample = |p, eta| branch.adiabatic_machine(&context, Pressure::new(p),
            GasMachineMode::Compressor, eta, 1e-8, || false).unwrap();
        let (dp, de) = (40.0, 0.0001);
        let (up, down) = (sample(400000.0 + dp, 0.8), sample(400000.0 - dp, 0.8));
        let dt = (up.outlet().state().temperature - down.outlet().state().temperature) / (2.0 * dp);
        let dw = (up.shaft_work_by_gas_j_kg - down.shaft_work_by_gas_j_kg) / (2.0 * dp);
        assert!((dt / result.derivatives.outlet_temperature_pressure_k_pa - 1.0).abs() < 1e-6);
        assert!((dw / result.derivatives.shaft_work_pressure_j_kg_pa - 1.0).abs() < 1e-6);
        let (up, down) = (sample(400000.0, 0.8 + de), sample(400000.0, 0.8 - de));
        let dt = (up.outlet().state().temperature - down.outlet().state().temperature) / (2.0 * de);
        let dw = (up.shaft_work_by_gas_j_kg - down.shaft_work_by_gas_j_kg) / (2.0 * de);
        assert!((dt / result.derivatives.outlet_temperature_efficiency_k - 1.0).abs() < 1e-6);
        assert!((dw / result.derivatives.shaft_work_efficiency_j_kg - 1.0).abs() < 1e-6);
    }
}

#[test]
fn hypothetical_isentropic_transport_does_not_block_a_supported_real_outlet() {
    let source = model_in_domains(vec![0.75, 0.25], 0.0,
        [200.0, 2000.0], [750.0, 1500.0], [1000.0, 1e7]);
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    let context = point(1000.0, 100000.0);
    assert!(branch.reconstruct_pressure_entropy_change(&context,
        Pressure::new(25000.0), 0.0, 1e-8, || false).is_err());
    let result = branch.adiabatic_machine(&context, Pressure::new(25000.0),
        GasMachineMode::Turbine, 0.5, 1e-8, || false).unwrap();
    assert!(result.isentropic_temperature_k < 750.0);
    assert!(result.outlet().state().temperature > 750.0);
    assert!(branch.adiabatic_machine(&context, Pressure::new(25000.0),
        GasMachineMode::Turbine, 1.0, 1e-8, || false).is_err());
}

#[test]
fn lossless_variable_cp_compression_and_expansion_close_the_same_enthalpy_cycle() {
    let source = shifted(1e24);
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    let compressor = branch.adiabatic_machine(&point(400.0, 100000.0), Pressure::new(400000.0),
        GasMachineMode::Compressor, 1.0, 1e-8, || false).unwrap();
    let turbine = branch.adiabatic_machine(compressor.outlet_point(), Pressure::new(100000.0),
        GasMachineMode::Turbine, 1.0, 1e-8, || false).unwrap();
    close(turbine.outlet().state().temperature, 400.0);
    assert!((compressor.shaft_work_by_gas_j_kg + turbine.shaft_work_by_gas_j_kg).abs()
        <= 2.0 * controls().energy_tolerance_j_kg);
}

#[test]
fn machine_refuses_wrong_directions_invalid_efficiency_and_all_cancelled_candidates() {
    let source = shifted(0.0);
    let saved = source.clone();
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    let context = point(400.0, 100000.0);
    for eta in [0.0, -1.0, 1.0001, f64::NAN, f64::INFINITY] {
        assert!(branch.adiabatic_machine(&context, Pressure::new(400000.0),
            GasMachineMode::Compressor, eta, 1e-8, || false).is_err());
    }
    assert!(branch.adiabatic_machine(&context, Pressure::new(50000.0),
        GasMachineMode::Compressor, 0.8, 1e-8, || false).is_err());
    assert!(branch.adiabatic_machine(&context, Pressure::new(200000.0),
        GasMachineMode::Turbine, 0.8, 1e-8, || false).is_err());
    let mut total = 0;
    branch.adiabatic_machine(&context, Pressure::new(400000.0), GasMachineMode::Compressor,
        0.8, 1e-8, || { total += 1; false }).unwrap();
    for stop in 1..=total {
        let mut polls = 0;
        assert!(branch.adiabatic_machine(&context, Pressure::new(400000.0), GasMachineMode::Compressor,
            0.8, 1e-8, || { polls += 1; polls == stop }).is_err());
    }
    assert_eq!(source, saved);
    let zero = branch.adiabatic_machine(&context, Pressure::new(100000.0), GasMachineMode::Compressor,
        0.8, 1e-8, || false).unwrap();
    assert_eq!(zero.inlet(), zero.outlet());
    assert_eq!(zero.enthalpy_change_j_kg, 0.0);
    assert_eq!(zero.shaft_work_by_gas_j_kg, 0.0);
    assert_eq!(zero.iterations, 0);
}
