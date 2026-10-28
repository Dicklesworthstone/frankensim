use super::*;
use super::tests::{controls, shifted};
use crate::thermochemical_gas::{ThermochemicalGasModel,
    tests::{model, point, close, species_in_domain, transport}};
use fs_thermochem::{Composition, CompositionBasis, MolarEnergyQuantityV1,
    MolarThermalQuantityV1, Nasa9StandardStateModelV1};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

#[test]
fn polytropic_constant_capacity_matches_local_efficiency_not_overall_efficiency() {
    let source = model(vec![0.75, 0.25], 0.0);
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    for (mode, t0, p1) in [(GasMachineMode::Compressor, 350.0, 400000.0),
        (GasMachineMode::Turbine, 1100.0, 25000.0)] {
        for eta in [0.6, 0.8, 1.0] {
            let context = point(t0, 100000.0).with("material-lot", 42.0).unwrap();
            let result = branch.polytropic_machine(&context, Pressure::new(p1), mode, eta, 1e-8, || false).unwrap();
            let gas = result.inlet().state();
            let exponent = (gas.gamma - 1.0) / gas.gamma;
            let factor = if mode == GasMachineMode::Compressor { 1.0 / eta } else { eta };
            let expected = t0 * (p1 / 100000.0_f64).powf(exponent * factor);
            close(result.outlet().state().temperature, expected);
            close(result.shaft_work_by_gas_j_kg, -gas.specific_heat_cp * (expected - t0));
            assert_eq!(result.declared_efficiency, GasMachineEfficiency::Polytropic(eta));
            assert_eq!(result.outlet_point().axes().get("material-lot"), Some(&42.0));
            assert!(result.polytropic_entropy_residual_j_kg_k.unwrap().abs() <= 1e-8);
            let overall = branch.adiabatic_machine(&context, Pressure::new(p1), mode, eta, 1e-8, || false).unwrap();
            if eta == 1.0 {
                assert_eq!(result.outlet(), overall.outlet());
                assert_eq!(result.isentropic_efficiency, 1.0);
            } else {
                assert!((result.outlet().state().temperature - overall.outlet().state().temperature).abs() > 1.0);
                if mode == GasMachineMode::Compressor {
                    assert!(result.isentropic_efficiency < eta);
                } else {
                    assert!(result.isentropic_efficiency > eta);
                }
                assert!(result.entropy_change_j_kg_k > 0.0);
            }
        }
    }
}

#[test]
fn constant_local_efficiency_is_invariant_to_pressure_subdivision_with_variable_cp() {
    let source = shifted(1e24);
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    for (mode, t0, p1) in [(GasMachineMode::Compressor, 400.0, 600000.0),
        (GasMachineMode::Turbine, 1200.0, 20000.0)] {
        let context = point(t0, 100000.0);
        let whole = branch.polytropic_machine(&context, Pressure::new(p1), mode, 0.8, 1e-8, || false).unwrap();
        for count in [2, 4, 8] {
            let mut current = context.clone();
            let mut work = 0.0;
            let mut entropy = 0.0;
            for index in 1..=count {
                let pressure = if index == count { p1 }
                    else { 100000.0 * (p1 / 100000.0_f64).powf(f64::from(index) / f64::from(count)) };
                let step = branch.polytropic_machine(&current, Pressure::new(pressure), mode,
                    0.8, 1e-8, || false).unwrap();
                work += step.shaft_work_by_gas_j_kg;
                entropy += step.entropy_change_j_kg_k;
                current = step.outlet_point().clone();
            }
            let final_state = source.evaluate(&current).unwrap();
            assert!((final_state.state().temperature - whole.outlet().state().temperature).abs() < 1e-8);
            assert!((work - whole.shaft_work_by_gas_j_kg).abs()
                <= f64::from(count + 1) * controls().energy_tolerance_j_kg);
            assert!((entropy - whole.entropy_change_j_kg_k).abs() <= f64::from(count + 1) * 1e-8);
        }
    }
}

#[test]
fn polytropic_variable_cp_derivatives_and_reference_shifts_use_the_same_state() {
    let ordinary = shifted(0.0);
    for (mode, t0, p1) in [(GasMachineMode::Compressor, 400.0, 400000.0),
        (GasMachineMode::Turbine, 1000.0, 25000.0)] {
        let context = point(t0, 100000.0);
        let expected = ordinary.prepare_caloric_branch(controls(), || false).unwrap()
            .polytropic_machine(&context, Pressure::new(p1), mode, 0.8, 1e-8, || false).unwrap();
        for reference in [-1e24, 1e24] {
            let source = shifted(reference);
            let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
            let evaluate = |p, eta| branch.polytropic_machine(&context, Pressure::new(p), mode,
                eta, 1e-8, || false).unwrap();
            let actual = evaluate(p1, 0.8);
            assert_eq!(actual.outlet().state(), expected.outlet().state());
            assert_eq!(actual.shaft_work_by_gas_j_kg, expected.shaft_work_by_gas_j_kg);
            assert_eq!(actual.derivatives, expected.derivatives);
            let (dp, de) = (p1 * 1e-4, 1e-4);
            let (plus, minus) = (evaluate(p1 + dp, 0.8), evaluate(p1 - dp, 0.8));
            let dt = (plus.outlet().state().temperature - minus.outlet().state().temperature) / (2.0 * dp);
            let dw = (plus.shaft_work_by_gas_j_kg - minus.shaft_work_by_gas_j_kg) / (2.0 * dp);
            assert!((dt / actual.derivatives.outlet_temperature_pressure_k_pa - 1.0).abs() < 1e-6);
            assert!((dw / actual.derivatives.shaft_work_pressure_j_kg_pa - 1.0).abs() < 1e-6);
            let (plus, minus) = (evaluate(p1, 0.8 + de), evaluate(p1, 0.8 - de));
            let dt = (plus.outlet().state().temperature - minus.outlet().state().temperature) / (2.0 * de);
            let dw = (plus.shaft_work_by_gas_j_kg - minus.shaft_work_by_gas_j_kg) / (2.0 * de);
            assert!((dt / actual.derivatives.outlet_temperature_efficiency_k - 1.0).abs() < 1e-6);
            assert!((dw / actual.derivatives.shaft_work_efficiency_j_kg - 1.0).abs() < 1e-6);
        }
    }
}

#[test]
fn polytropic_path_integrates_original_heat_capacity_through_an_admitted_seam() {
    let low = species_in_domain("A", 0.028, 3.5, 0.0, 0.0, [200.0, 600.0]);
    let high = species_in_domain("A", 0.028, 4.0, 0.0, 0.0, [600.0, 2000.0]);
    let species = Nasa9StandardStateModelV1::new(low.species().clone(), low.molar_mass(),
        low.convention().clone(), vec![low.regions()[0].clone(), high.regions()[0].clone()]).unwrap();
    let raw = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport("A", 1.7e-5, 0.025)], None).unwrap();
    assert!(raw.prepare_caloric_branch(controls(), || false).is_err());
    let source = raw.with_caloric_continuity(Temperature::new(300.0),
        MolarEnergyQuantityV1::new(1e5), MolarThermalQuantityV1::new(1e3)).unwrap();
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    let pressure = 100000.0 * (0.8 * (3.5 * 2.0_f64.ln() + 4.0 * 1.5_f64.ln())).exp();
    let step = branch.polytropic_machine(&point(300.0, 100000.0), Pressure::new(pressure),
        GasMachineMode::Compressor, 0.8, 1e-8, || false).unwrap();
    close(step.outlet().state().temperature, 900.0);
    let r = step.inlet().state().specific_gas_constant;
    close(step.shaft_work_by_gas_j_kg, -r * (3.5 * 300.0 + 4.0 * 300.0));
}

#[test]
fn polytropic_validation_cancellation_and_equal_pressure_are_atomic() {
    let source = shifted(1e24);
    let saved = source.clone();
    let branch = source.prepare_caloric_branch(controls(), || false).unwrap();
    let context = point(400.0, 100000.0);
    for eta in [0.0, -1.0, 1.001, f64::NAN, f64::INFINITY] {
        assert!(branch.polytropic_machine(&context, Pressure::new(400000.0),
            GasMachineMode::Compressor, eta, 1e-8, || false).is_err());
    }
    assert!(branch.polytropic_machine(&context, Pressure::new(50000.0),
        GasMachineMode::Compressor, 0.8, 1e-8, || false).is_err());
    assert!(branch.polytropic_machine(&context, Pressure::new(1e7),
        GasMachineMode::Compressor, 0.1, 1e-8, || false).is_err());
    let mut total = 0;
    branch.polytropic_machine(&context, Pressure::new(400000.0), GasMachineMode::Compressor,
        0.8, 1e-8, || { total += 1; false }).unwrap();
    for stop in 1..=total {
        let mut polls = 0;
        assert!(branch.polytropic_machine(&context, Pressure::new(400000.0), GasMachineMode::Compressor,
            0.8, 1e-8, || { polls += 1; polls == stop }).is_err());
    }
    for mode in [GasMachineMode::Compressor, GasMachineMode::Turbine] {
        let zero = branch.polytropic_machine(&context, Pressure::new(100000.0), mode,
            0.8, 1e-8, || false).unwrap();
        assert_eq!(zero.inlet(), zero.outlet());
        assert_eq!(zero.shaft_work_by_gas_j_kg, 0.0);
        assert_eq!(zero.polytropic_entropy_residual_j_kg_k, Some(0.0));
        assert_eq!(zero.isentropic_efficiency, 0.8);
    }
    assert_eq!(source, saved);
}
