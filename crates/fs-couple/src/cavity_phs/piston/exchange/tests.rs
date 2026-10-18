//! Original material owners, physical inventory balances and accepted mechanics.
use super::*;
use super::super::tests::{controls, drive, gas};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_ports::{ModuleWindowAudit, PortedModule};

fn exchange() -> GasNozzleStepSolve {
    GasNozzleStepSolve { nozzle: controls().gas, maximum_donor_fraction: 0.1,
        transfer: GasMassTransferSolve { donor: controls().gas, receiver: controls().gas,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn nozzle() -> GasNozzle { GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 } }
fn hot_reservoir(x: f64) -> ClosedGasParcel {
    let mut reservoir = gas(x, 0.001);
    reservoir.add_heat(100.0, controls().gas, || false).unwrap();
    reservoir
}
fn step(piston: &mut GasPiston, reservoir: &mut ClosedGasParcel) -> PistonNozzleStep {
    piston.step_with_nozzle(reservoir, nozzle(), drive(0.0), exchange(), controls(),
        1e-7, || false).unwrap()
}

#[test]
fn influx_composition_drives_the_actual_pressure_sound_speed_and_wall() {
    let mut piston = GasPiston::new(gas(0.99, 0.001), 1.0, 0.01, 0.0).unwrap();
    let mut reservoir = hot_reservoir(0.01);
    let original = piston.clone();
    let original_reservoir = reservoir.clone();
    let report = step(&mut piston, &mut reservoir);
    assert_eq!(report.direction, Some(NozzleDirection::IntoPiston));
    assert!(report.transfer.transferred_mass_kg > 0.0);
    assert!(piston.gas().mass_kg() > original.gas().mass_kg());
    assert!(reservoir.mass_kg() < original_reservoir.mass_kg());
    assert!(piston.velocity_m_s() > 0.0);
    assert!(piston.gas().state().state().pressure > original.gas().state().state().pressure);
    assert_ne!(piston.gas().state().state().sound_speed, original.gas().state().state().sound_speed);
    assert_ne!(piston.gas().state().state().specific_gas_constant,
        original.gas().state().state().specific_gas_constant);
    assert!(report.transfer.species.iter().any(|species| species.mass_kg > 0.0));
    assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
    assert!(report.total_mass_residual_kg.abs() <= exchange().transfer.mass_tolerance_kg);
    assert!(report.total_energy_residual_j.abs() <= 1e-7);
    // Enthalpy is deposited once, including upstream flow work and jet energy.
    let enthalpy = report.transfer.transferred_enthalpy_j;
    assert!((report.reservoir_energy_change_j + enthalpy).abs() < 1e-8);
    assert!((report.chamber_energy_change_j + report.piston.gas.work_by_gas_j - enthalpy).abs() < 1e-7);
    assert_eq!(report.piston.gas.supplied_heat_j, 0.0);
    assert!(matches!(piston.audit_last_window().unwrap(), ModuleWindowAudit::Unavailable { .. }));
    assert_eq!(piston.last_window_transfers().unwrap().len(), 2);

    // The subsequent mechanics uses the same accepted mixture and source state.
    let composition = piston.gas().state().thermodynamics().receipt().mole_fraction_bits().to_vec();
    piston.step(drive(0.0), controls(), || false).unwrap();
    assert_eq!(piston.gas().state().thermodynamics().receipt().mole_fraction_bits(), composition.as_slice());
}

#[test]
fn pressure_reversal_withdraws_chamber_species_and_enthalpy_conservatively() {
    let mut piston = GasPiston::new(hot_reservoir(0.75), 1.0, 0.01, 0.0).unwrap();
    let mut reservoir = gas(0.25, 0.001);
    let old_chamber = piston.gas().clone();
    let old_reservoir = reservoir.clone();
    let report = step(&mut piston, &mut reservoir);
    assert_eq!(report.direction, Some(NozzleDirection::OutOfPiston));
    assert!(piston.gas().mass_kg() < old_chamber.mass_kg());
    assert!(reservoir.mass_kg() > old_reservoir.mass_kg());
    assert_eq!(report.flow.upstream_enthalpy_j_kg, old_chamber.state().specific_enthalpy_j_kg());
    assert_ne!(reservoir.state().state().specific_gas_constant, old_reservoir.state().state().specific_gas_constant);
    assert!((report.reservoir_energy_change_j - report.transfer.transferred_enthalpy_j).abs() < 1e-8);
    assert!((report.chamber_energy_change_j + report.piston.gas.work_by_gas_j
        + report.transfer.transferred_enthalpy_j).abs() < 1e-7);
    assert!(report.total_mass_residual_kg.abs() < 1e-14);
    assert!(report.total_energy_residual_j.abs() < 1e-7);
}

#[test]
fn closed_nozzle_is_the_identical_ordinary_step_including_port_history() {
    for closed in [GasNozzle { area_m2: 0.0, ..nozzle() },
        GasNozzle { discharge_coefficient: 0.0, ..nozzle() }] {
        let mut actual = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
        let mut expected = actual.clone();
        let mut reservoir = hot_reservoir(0.25);
        let old_reservoir = reservoir.clone();
        let reference = expected.step(drive(0.2), controls(), || false).unwrap();
        let report = actual.step_with_nozzle(&mut reservoir, closed, drive(0.2), exchange(), controls(),
            1e-7, || false).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(actual.ports, expected.ports);
        assert_eq!(reservoir, old_reservoir);
        assert_eq!(report.piston, reference);
        assert_eq!(report.direction, None);
        assert_eq!(report.transfer.transferred_mass_kg, 0.0);
        assert_eq!(report.reservoir_energy_change_j, 0.0);
    }
}

#[test]
fn refusal_after_transfer_and_final_cancellation_preserve_both_states_and_accounts() {
    let mut initial = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    initial.step(drive(0.0), controls(), || false).unwrap();
    let reservoir = hot_reservoir(0.25);
    let mut piston = initial.clone();
    let mut source = reservoir.clone();
    let wrong_bracket = PistonSolve { volume_bracket_m3: [0.0011, 0.0012], ..controls() };
    let mut checked = 0;
    let refused = piston.step_with_nozzle(&mut source, nozzle(), drive(0.0), exchange(), wrong_bracket,
        1e-7, || { checked += 1; false });
    assert!(matches!(refused, Err(PistonError::Invalid(_))));
    assert!(checked > 5); // the original material transfer ran before mechanical admission
    assert_eq!(piston, initial); assert_eq!(piston.ports, initial.ports); assert_eq!(source, reservoir);

    let mut complete = initial.clone(); let mut complete_source = reservoir.clone(); let mut calls = 0;
    complete.step_with_nozzle(&mut complete_source, nozzle(), drive(0.0), exchange(), controls(),
        1e-7, || { calls += 1; false }).unwrap();
    let mut tick = 0;
    let result = piston.step_with_nozzle(&mut source, nozzle(), drive(0.0), exchange(), controls(),
        1e-7, || { tick += 1; tick == calls });
    assert_eq!(result, Err(PistonError::Cancelled));
    assert_eq!(piston, initial); assert_eq!(piston.ports, initial.ports); assert_eq!(source, reservoir);
}

#[test]
fn transfer_then_motion_refines_as_a_first_order_physical_split() {
    let mut velocities = Vec::new();
    for count in [16, 32, 64] {
        let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
        let mut reservoir = hot_reservoir(0.25);
        let input = PistonInput { duration_s: 0.002 / f64::from(count), ..drive(0.0) };
        let solve = PistonSolve { momentum_tolerance_n_s: 2e-7 * input.duration_s, ..controls() };
        for _ in 0..count {
            let report = piston.step_with_nozzle(&mut reservoir, nozzle(), input, exchange(), solve,
                1e-7, || false).unwrap();
            assert!(report.total_mass_residual_kg.abs() < 1e-14);
            assert!(report.total_energy_residual_j.abs() < 1e-7);
        }
        velocities.push(piston.velocity_m_s());
    }
    let coarse = (velocities[0] - velocities[1]).abs();
    let fine = (velocities[1] - velocities[2]).abs();
    assert!(fine > 0.0 && coarse / fine > 1.6 && coarse / fine < 2.6,
        "velocity refinement {velocities:?}, ratio {}", coarse / fine);
}
