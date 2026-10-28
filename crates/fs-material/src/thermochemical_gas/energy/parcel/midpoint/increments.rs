use super::*;
use crate::thermochemical_gas::{ThermochemicalGasModel, tests::{point, species_with_energy_reference, transport}};
use fs_thermochem::{Composition, CompositionBasis};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn parcel(reference: f64) -> ClosedGasParcel {
    let source = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", 0.028, 3.5, 0.002, reference),
            species_with_energy_reference("B", 0.004, 2.5, 0.0, -reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    ).unwrap(), vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap();
    ClosedGasParcel::new(source, point(400.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.001).unwrap()
}

#[test]
fn moving_volume_heating_and_discrete_tangents_are_reference_invariant() {
    let baseline = parcel(0.0);
    for volume in [0.0007, 0.001, 0.0015] {
        for heat in [-10.0, 0.0, 20.0] {
            let expected = baseline.midpoint_heat_volume(volume, heat, controls(), || false).unwrap();
            for reference in [-1e24, 1e24] {
                let gas = parcel(reference);
                let saved = gas.clone();
                let step = gas.midpoint_heat_volume(volume, heat, controls(), || false).unwrap();
                assert_eq!(step.step(), expected.step());
                assert_eq!(step.midpoint().state(), expected.midpoint().state());
                assert_eq!(step.endpoint().state().state(), expected.endpoint().state().state());
                assert_eq!(step.endpoint().point().axes().get("material-lot"), Some(&42.0));
                assert!(step.endpoint().state().specific_internal_energy_j_kg().abs() > 1e20);
                assert_eq!(gas, saved);
            }
        }
    }
}

#[test]
fn large_reference_midpoint_keeps_physical_derivatives_and_energy_midpoint() {
    let gas = parcel(1e24);
    let (volume, heat) = (0.0008, 25.0);
    let result = gas.midpoint_heat_volume(volume, heat, controls(), || false).unwrap();
    let evaluate = |v, q| gas.midpoint_heat_volume(v, q, controls(), || false).unwrap().midpoint().state().pressure;
    let dv = 1e-8;
    let dq = 1e-3;
    let fdv = (evaluate(volume + dv, heat) - evaluate(volume - dv, heat)) / (2.0 * dv);
    let fdq = (evaluate(volume, heat + dq) - evaluate(volume, heat - dq)) / (2.0 * dq);
    assert!((fdv / result.step().pressure_final_volume_derivative_pa_m3 - 1.0).abs() < 1e-7);
    assert!((fdq / result.step().pressure_heat_derivative_pa_j - 1.0).abs() < 1e-7);
    let delta = |t| gas.model().thermodynamics().specific_internal_energy_change(
        Temperature::new(400.0), Temperature::new(t)).unwrap().value();
    assert!((2.0 * delta(result.midpoint().state().temperature)
        - delta(result.endpoint().state().state().temperature)).abs() <= controls().energy_tolerance_j_kg);
    assert!(result.step().energy_residual_j.abs() <= gas.mass_kg() * controls().energy_tolerance_j_kg);
}

#[test]
fn finite_isothermal_motion_retains_pressure_work_at_large_reference() {
    let gas = parcel(1e24);
    let volume = 0.0014;
    let vm = 0.5 * (volume + gas.volume_m3());
    let heat = gas.mass_kg() * gas.state().state().specific_gas_constant * 400.0
        * ((volume - gas.volume_m3()) / vm);
    let result = gas.midpoint_heat_volume(volume, heat, controls(), || false).unwrap();
    assert!((result.endpoint().state().state().temperature - 400.0).abs() < 1e-9);
    assert!((result.step().work_by_gas_j - heat).abs() <= gas.mass_kg() * controls().energy_tolerance_j_kg);
    assert_eq!(result.endpoint().volume_m3(), volume);
    assert_eq!(result.endpoint().mass_kg(), gas.mass_kg());
}
