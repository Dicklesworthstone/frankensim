use super::*;
use crate::thermochemical_gas::{ThermochemicalGasModel, tests::{close, model, point,
    species_in_domain, species_with_energy_reference, transport}};
use crate::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_thermochem::{Composition, CompositionBasis, MolarEnergyQuantityV1,
    MolarThermalQuantityV1, Nasa9StandardStateModelV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}

fn shifted(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", 0.028, 3.5, 0.001, reference),
            species_with_energy_reference("B", 0.004, 2.5, 0.0, -reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    ).unwrap(), vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap()
}

fn specific_r() -> f64 { R / (0.75 * 0.028 + 0.25 * 0.004) }
fn dh(a: f64, b: f64, slope: f64) -> f64 {
    specific_r() * (b - a) * (3.25 + 0.375 * slope * (b + a))
}

#[test]
fn g1_isobaric_heating_and_cooling_close_independent_work_and_storage() {
    for slope in [0.0, 0.001] {
        let mut gas = ClosedGasParcel::new(model(vec![0.75, 0.25], slope),
            point(300.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.01).unwrap();
        let mut fixed_volume = gas.clone();
        let mass = gas.mass_kg();
        let heat = mass * dh(300.0, 800.0, slope);
        let step = gas.add_heat_at_constant_pressure(heat, controls(), || false).unwrap();
        close(step.final_temperature_k, 800.0);
        assert_eq!(gas.state().state().pressure, 100000.0);
        assert_eq!(gas.mass_kg(), mass);
        assert_eq!(gas.point().axes().get("material-lot"), Some(&42.0));
        close(gas.volume_m3(), 0.01 * 800.0 / 300.0);
        close(step.work_by_gas_j, mass * specific_r() * 500.0);
        close(step.stored_energy_change_j, heat - mass * specific_r() * 500.0);
        close(step.enthalpy_change_j, heat);
        let ds = mass * specific_r()
            * (3.25 * (800.0_f64 / 300.0).ln() + 0.75 * slope * 500.0);
        close(step.entropy_change_j_k, ds);
        assert!(step.energy_residual_j.abs() <= mass * controls().energy_tolerance_j_kg);
        assert!(step.enthalpy_residual_j.abs() <= mass * controls().energy_tolerance_j_kg);
        assert_eq!(step.work_by_gas_j, step.pressure_pa * (step.final_volume_m3 - step.initial_volume_m3));
        gas.state().state().validate_coherence().unwrap();
        fixed_volume.add_heat(heat, controls(), || false).unwrap();
        assert!(fixed_volume.state().state().temperature > gas.state().state().temperature);
        assert!(fixed_volume.state().state().pressure > gas.state().state().pressure);
        assert_eq!(fixed_volume.volume_m3(), 0.01);
        let cool = gas.add_heat_at_constant_pressure(-heat, controls(), || false).unwrap();
        close(gas.state().state().temperature, 300.0);
        close(gas.volume_m3(), 0.01);
        assert!(cool.work_by_gas_j < 0.0 && cool.stored_energy_change_j < 0.0);
        assert!(cool.entropy_change_j_k < 0.0);
        assert!(cool.energy_residual_j.abs() <= mass * controls().energy_tolerance_j_kg);
        close(cool.work_by_gas_j, -step.work_by_gas_j);
    }
}

#[test]
fn g3_isobaric_heat_and_isentropic_work_ignore_common_formation_offsets() {
    for target in [300.0 + 1e-6, 800.0] {
        let mut baseline = ClosedGasParcel::new(shifted(0.0), point(300.0, 100000.0), 0.01).unwrap();
        let heat = baseline.mass_kg() * dh(300.0, target, 0.001);
        let expected = baseline.add_heat_at_constant_pressure(heat, controls(), || false).unwrap();
        for reference in [-1e24, 1e24] {
            let mut gas = ClosedGasParcel::new(shifted(reference), point(300.0, 100000.0), 0.01).unwrap();
            let step = gas.add_heat_at_constant_pressure(heat, controls(), || false).unwrap();
            assert_eq!(step, expected);
            assert_eq!(gas.state().state(), baseline.state().state());
            assert!(gas.state().specific_enthalpy_j_kg().abs() > 1e20);
        }
    }
    let entropy = GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        entropy_tolerance_j_kg_k: 1e-8, maximum_iterations: 128 };
    let mut baseline = ClosedGasParcel::new(shifted(0.0), point(400.0, 100000.0), 0.01).unwrap();
    let expected = baseline.set_volume_isentropic(0.005, entropy, || false).unwrap();
    for reference in [-1e24, 1e24] {
        let mut gas = ClosedGasParcel::new(shifted(reference), point(400.0, 100000.0), 0.01).unwrap();
        let step = gas.set_volume_isentropic(0.005, entropy, || false).unwrap();
        assert_eq!(step, expected);
        assert!(step.work_by_gas_j < 0.0);
        assert_eq!(gas.state().state(), baseline.state().state());
        let reverse = gas.set_volume_isentropic(0.01, entropy, || false).unwrap();
        let reversible_tolerance = 2.0 * gas.mass_kg() * entropy.maximum_temperature_k
            * entropy.entropy_tolerance_j_kg_k + 128.0 * f64::EPSILON * step.work_by_gas_j.abs();
        assert!((reverse.work_by_gas_j + step.work_by_gas_j).abs() <= reversible_tolerance);
    }
}

#[test]
fn g3_isobaric_heat_crosses_only_an_explicitly_reconciled_source_seam() {
    let low = species_in_domain("A", 0.028, 3.5, 0.0, 0.0, [200.0, 600.0]);
    let high = species_in_domain("A", 0.028, 4.0, 0.0, 0.0, [600.0, 2000.0]);
    let species = Nasa9StandardStateModelV1::new(low.species().clone(), low.molar_mass(),
        low.convention().clone(), vec![low.regions()[0].clone(), high.regions()[0].clone()]).unwrap();
    let raw = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport("A", 1.7e-5, 0.025)], None).unwrap();
    let mut gas = ClosedGasParcel::new(raw.clone(), point(300.0, 100000.0), 0.01).unwrap();
    let saved = gas.clone();
    let heat = gas.mass_kg() * (R / 0.028) * (3.5 * 300.0 + 4.0 * 200.0);
    assert!(gas.add_heat_at_constant_pressure(heat, controls(), || false).is_err());
    assert_eq!(gas, saved);
    let continuous = raw.with_caloric_continuity(Temperature::new(300.0),
        MolarEnergyQuantityV1::new(1e5), MolarThermalQuantityV1::new(1e3)).unwrap();
    let mut gas = ClosedGasParcel::new(continuous, point(300.0, 100000.0), 0.01).unwrap();
    let step = gas.add_heat_at_constant_pressure(heat, controls(), || false).unwrap();
    close(gas.state().state().temperature, 800.0);
    close(step.work_by_gas_j, gas.mass_kg() * (R / 0.028) * 500.0);
    assert!(step.energy_residual_j.abs() <= gas.mass_kg() * controls().energy_tolerance_j_kg);
}

#[test]
fn g4_isobaric_failures_and_late_cancellation_preserve_all_state() {
    let mut gas = ClosedGasParcel::new(shifted(1e24), point(300.0, 100000.0), 0.01).unwrap();
    let saved = gas.clone();
    for (heat, limits) in [
        (f64::NAN, controls()), (1e9, controls()),
        (100.0, GasEnergySolve { maximum_iterations: 1, ..controls() }),
        (100.0, GasEnergySolve { minimum_temperature_k: 400.0, ..controls() }),
    ] {
        assert!(gas.add_heat_at_constant_pressure(heat, limits, || false).is_err());
        assert_eq!(gas, saved);
    }
    let mut completed = 0;
    gas.clone().add_heat_at_constant_pressure(100.0, controls(), || { completed += 1; false }).unwrap();
    for stop in [1, completed / 2, completed] {
        let mut polls = 0;
        assert!(matches!(gas.add_heat_at_constant_pressure(100.0, controls(), || {
            polls += 1; polls == stop
        }), Err(GasEnergyError::Cancelled)));
        assert_eq!(gas, saved);
    }
    let source = crate::thermochemical_gas::tests::model_in_domains(vec![0.75, 0.25],
        0.001, [200.0, 2000.0], [200.0, 500.0], [1000.0, 1e7]);
    let mut gas = ClosedGasParcel::new(source, point(300.0, 100000.0), 0.01).unwrap();
    let saved = gas.clone();
    let heat = gas.mass_kg() * dh(300.0, 800.0, 0.001);
    assert!(matches!(gas.add_heat_at_constant_pressure(heat, controls(), || false),
        Err(GasEnergyError::Material(_))));
    assert_eq!(gas, saved);
}

#[test]
fn g0_isobaric_zero_heat_is_exact_but_does_not_skip_control_validation() {
    let mut gas = ClosedGasParcel::new(shifted(1e24), point(300.0, 100000.0), 0.01).unwrap();
    let saved = gas.clone();
    let result = gas.add_heat_at_constant_pressure(-0.0, controls(), || false).unwrap();
    assert_eq!(gas, saved);
    assert_eq!(result.work_by_gas_j, 0.0);
    assert_eq!(result.iterations, 0);
    assert!(gas.add_heat_at_constant_pressure(0.0,
        GasEnergySolve { maximum_iterations: 0, ..controls() }, || false).is_err());
    assert_eq!(gas, saved);
    assert!(matches!(gas.add_heat_at_constant_pressure(0.0, controls(), || true),
        Err(GasEnergyError::Cancelled)));
    assert_eq!(gas, saved);
    // A loose tolerance must not silently round a nonzero heat input to no motion.
    assert!(gas.add_heat_at_constant_pressure(f64::from_bits(1), controls(), || false).is_err());
    assert_eq!(gas, saved);
}
