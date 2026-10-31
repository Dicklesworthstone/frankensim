use super::*;
use super::super::super::nonideal_tests::{controls, drive, parcel};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::nonideal::energy::parcel::transfer::RealGasMassTransferSolve;
use fs_ports::{ModuleWindowAudit, PortedModule};

fn exchange() -> RealGasNozzleStepSolve {
    RealGasNozzleStepSolve {
        entropy: GasEntropySolve {
            minimum_temperature_k: 240.0,
            maximum_temperature_k: 450.0,
            entropy_tolerance_j_kg_k: 1e-9,
            maximum_iterations: 128,
        },
        transfer: RealGasMassTransferSolve {
            donor: controls().gas,
            receiver: controls().gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        minimum_density_fraction: 0.2,
        pressure_tolerance_pa: 1e-4,
        sonic_tolerance_m2_s2: 1e-4,
        maximum_iterations: 128,
        maximum_donor_fraction: 0.05,
    }
}
fn nozzle() -> GasNozzle {
    GasNozzle { area_m2: 2e-7, discharge_coefficient: 0.8 }
}
fn hotter(mut gas: IsochoricRealGasParcel) -> IsochoricRealGasParcel {
    let heat = gas.mass_kg() * 30_000.0;
    gas.add_heat(heat, controls().gas, || false).unwrap();
    gas
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:.15e} != {b:.15e}, tol={tolerance:.3e}");
}

#[test]
fn dense_reservoir_inflow_drives_the_same_nonideal_piston_pressure_and_energy() {
    let chamber = parcel();
    let mut reservoir = hotter(parcel());
    assert!(reservoir.state().pressure_pa() > chamber.state().pressure_pa());
    let old_chamber_mass = chamber.mass_kg();
    let old_reservoir_mass = reservoir.mass_kg();
    let mut piston = GasPiston::new(chamber, 2.0, 0.01, 0.0).unwrap();
    let input = PistonInput { external_force_n: 0.0, heat_j: 10.0, ..drive(0.0) };
    let report = piston.step_with_real_gas_nozzle(
        &mut reservoir, nozzle(), input, exchange(), controls(), 2e-6, || false,
    ).unwrap();

    assert_eq!(report.direction, Some(NozzleDirection::IntoPiston));
    assert!(report.flow.mass_flow_kg_s > 0.0);
    assert!(report.transfer.transferred_mass_kg > 0.0);
    assert!(piston.gas().mass_kg() > old_chamber_mass);
    assert!(reservoir.mass_kg() < old_reservoir_mass);
    assert!(piston.gas().state().pressure_pa() > 0.0);
    assert!(report.piston.midpoint.pressure_pa() > 0.0);
    close(
        report.transfer.transferred_enthalpy_j,
        report.transfer.transferred_mass_kg * report.flow.upstream_enthalpy_j_kg,
        1e-9,
    );
    close(
        report.total_energy_residual_j,
        report.transfer.energy_residual_j + report.piston.total_energy_residual_j,
        1e-15,
    );
    assert!(report.total_mass_residual_kg.abs() <= exchange().transfer.mass_tolerance_kg);
    assert!(report.total_energy_residual_j.abs() <= 2e-6);
    assert!(matches!(
        piston.audit_last_window().unwrap(),
        ModuleWindowAudit::Unavailable { .. }
    ));
}

#[test]
fn hotter_chamber_exhausts_real_gas_and_retains_pair_balance() {
    let chamber = hotter(parcel());
    let mut reservoir = parcel();
    assert!(chamber.state().pressure_pa() > reservoir.state().pressure_pa());
    let old_chamber = chamber.mass_kg();
    let old_reservoir = reservoir.mass_kg();
    let mut piston = GasPiston::new(chamber, 2.0, 0.01, 0.0).unwrap();
    let report = piston.step_with_real_gas_nozzle(
        &mut reservoir, nozzle(),
        PistonInput { external_force_n: 0.0, ..drive(0.0) },
        exchange(), controls(), 2e-6, || false,
    ).unwrap();

    assert_eq!(report.direction, Some(NozzleDirection::OutOfPiston));
    assert!(piston.gas().mass_kg() < old_chamber);
    assert!(reservoir.mass_kg() > old_reservoir);
    assert_eq!(
        report.flow.upstream_enthalpy_j_kg,
        report.transfer.source_enthalpy_j_kg
    );
    assert!(report.total_mass_residual_kg.abs() <= exchange().transfer.mass_tolerance_kg);
    assert!(report.total_energy_residual_j.abs() <= 2e-6);
}

#[test]
fn shut_path_matches_the_native_closed_real_gas_step() {
    let chamber = parcel();
    let mut actual = GasPiston::new(chamber.clone(), 2.0, 0.01, 0.0).unwrap();
    let mut expected = actual.clone();
    let mut reservoir = hotter(parcel());
    let old_reservoir = reservoir.clone();
    let input = PistonInput { external_force_n: 0.0, heat_j: 5.0, ..drive(0.0) };
    let reference = expected.step(input, controls(), || false).unwrap();
    let report = actual.step_with_real_gas_nozzle(
        &mut reservoir,
        GasNozzle { area_m2: 0.0, discharge_coefficient: 0.8 },
        input, exchange(), controls(), 2e-6, || false,
    ).unwrap();

    assert_eq!(actual, expected);
    assert_eq!(reservoir, old_reservoir);
    assert_eq!(report.direction, None);
    assert_eq!(report.transfer.transferred_mass_kg, 0.0);
    assert_eq!(report.piston, reference);
    assert_eq!(report.reservoir_energy_change_j, 0.0);
}

#[test]
fn post_transfer_mechanical_refusal_and_final_cancellation_publish_nothing() {
    let chamber = parcel();
    let reservoir = hotter(parcel());
    let mut piston = GasPiston::new(chamber, 2.0, 0.01, 0.0).unwrap();
    piston.step(PistonInput { external_force_n: 0.0, ..drive(0.0) }, controls(), || false).unwrap();
    let initial = piston.clone();
    let old_reservoir = reservoir.clone();

    let mut source = reservoir.clone();
    let wrong = PistonSolve { volume_bracket_m3: [0.0011, 0.0012], ..controls() };
    assert!(piston.step_with_real_gas_nozzle(
        &mut source, nozzle(), PistonInput { external_force_n: 0.0, ..drive(0.0) },
        exchange(), wrong, 2e-6, || false,
    ).is_err());
    assert_eq!(piston, initial);
    assert_eq!(source, old_reservoir);

    let mut probe = initial.clone();
    let mut probe_source = old_reservoir.clone();
    let mut polls = 0;
    probe.step_with_real_gas_nozzle(
        &mut probe_source, nozzle(), PistonInput { external_force_n: 0.0, ..drive(0.0) },
        exchange(), controls(), 2e-6, || { polls += 1; false },
    ).unwrap();
    assert!(polls > 5);
    let mut seen = 0;
    let result = piston.step_with_real_gas_nozzle(
        &mut source, nozzle(), PistonInput { external_force_n: 0.0, ..drive(0.0) },
        exchange(), controls(), 2e-6, || {
            seen += 1;
            seen == polls
        },
    );
    assert_eq!(result, Err(PistonError::Cancelled));
    assert_eq!(piston, initial);
    assert_eq!(source, old_reservoir);
}
