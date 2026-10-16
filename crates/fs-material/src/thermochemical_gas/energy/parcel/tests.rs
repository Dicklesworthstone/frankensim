//! Synthetic thermodynamic identities and rollback, not combustion validation.
use super::*;
use crate::fluid_transport::FluidTransport;
use crate::thermochemical_gas::tests::{close, model, point};

fn entropy_controls() -> GasEntropySolve {
    GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
}
fn energy_controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn parcel(slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![0.75, 0.25], slope),
        point(300.0, 100000.0).with("material-lot", 12.0).unwrap(), 0.01).unwrap()
}

#[test]
fn compression_and_expansion_match_independent_pressure_work() {
    let mut parcel = parcel(0.0);
    let before = parcel.state().clone();
    let mass = parcel.mass_kg();
    let gamma = before.state().gamma;
    let step = parcel.set_volume_isentropic(0.0025, entropy_controls(), || false).unwrap();
    close(parcel.state().state().temperature, 300.0*4.0_f64.powf(gamma-1.0));
    close(parcel.state().state().pressure, 100000.0*4.0_f64.powf(gamma));
    close(parcel.state().state().density, before.state().density*4.0);
    // Independently integrated p(V)=p0*(V0/V)^gamma.
    let work = (100000.0*0.01-parcel.state().state().pressure*0.0025)/(gamma-1.0);
    close(step.work_by_gas_j, work);
    assert!(step.work_by_gas_j < 0.0);
    assert!(step.entropy_change_j_k.abs() <= mass*entropy_controls().entropy_tolerance_j_kg_k);
    assert!(parcel.state().state().sound_speed > before.state().sound_speed);
    assert!(parcel.state().state().dynamic_viscosity > before.state().dynamic_viscosity);
    assert!(parcel.state().state().thermal_conductivity > before.state().thermal_conductivity);
    close(parcel.state().transport_state().unwrap().density_kg_m3(), parcel.state().state().density);
    let reversed = parcel.set_volume_isentropic(0.01, entropy_controls(), || false).unwrap();
    close(reversed.work_by_gas_j, -step.work_by_gas_j);
    close(parcel.state().state().temperature, 300.0);
    close(parcel.state().state().pressure, 100000.0);
    assert_eq!(parcel.mass_kg(), mass);
    assert_eq!(parcel.point().axes().get("material-lot"), Some(&12.0));
}

#[test]
fn material_data_alone_changes_compression_thermal_and_acoustic_response() {
    let mut constant = parcel(0.0);
    let mut varying = parcel(0.001);
    constant.set_volume_isentropic(0.0025, entropy_controls(), || false).unwrap();
    varying.set_volume_isentropic(0.0025, entropy_controls(), || false).unwrap();
    assert_eq!(constant.mass_kg(), varying.mass_kg());
    assert_eq!(constant.volume_m3(), varying.volume_m3());
    assert!(varying.state().state().temperature < constant.state().state().temperature);
    assert!(varying.state().state().pressure < constant.state().state().pressure);
    assert!(varying.state().state().sound_speed < constant.state().state().sound_speed);
    assert!(varying.state().state().thermal_conductivity < constant.state().state().thermal_conductivity);
    // Same change in composition, with the same general compression operation.
    let mut different = ClosedGasParcel::new(model(vec![0.25, 0.75], 0.001),
        point(300.0, 100000.0), 0.01).unwrap();
    different.set_volume_isentropic(0.0025, entropy_controls(), || false).unwrap();
    assert!(different.state().state().sound_speed > varying.state().state().sound_speed);
    assert_ne!(different.mass_kg(), varying.mass_kg());
}

#[test]
fn variable_cp_boundary_work_matches_pressure_integral() {
    let mut parcel = parcel(0.001);
    let source = parcel.parcel.model.clone();
    let initial = parcel.state().clone();
    let point = parcel.point().clone();
    let mass = parcel.mass_kg();
    let entropy = initial.specific_entropy_j_kg_k().unwrap();
    let report = parcel.set_volume_isentropic(0.005, entropy_controls(), || false).unwrap();
    // Independent composite Simpson quadrature of source-resolved p dV.
    // This checks the energy-derived work, not a second copy of that formula.
    let intervals = 128;
    let h = (0.005-0.01)/f64::from(intervals);
    let mut sum = 0.0;
    for i in 0..=intervals {
        let volume = 0.01+f64::from(i)*h;
        let value = source.reconstruct_density_entropy(&point, mass/volume, entropy,
            entropy_controls(), || false).unwrap();
        let weight = if i == 0 || i == intervals { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += weight*value.state().state().pressure;
    }
    let integrated_work = sum*h/3.0;
    assert!((report.work_by_gas_j/integrated_work-1.0).abs() < 1e-8);
}

#[test]
fn heating_and_volume_work_compose_into_a_closed_idealized_cycle() {
    let mut parcel = parcel(0.0);
    let initial = parcel.state().clone();
    let mass = parcel.mass_kg();
    let cv = initial.thermodynamics().properties().mass_specific_cv().value();
    let compression = parcel.set_volume_isentropic(0.0025, entropy_controls(), || false).unwrap();
    let heat_in = mass*cv*(1000.0-parcel.state().state().temperature);
    let heated = parcel.add_heat(heat_in, energy_controls(), || false).unwrap();
    close(parcel.state().state().temperature, 1000.0);
    let expansion = parcel.set_volume_isentropic(0.01, entropy_controls(), || false).unwrap();
    let heat_out = mass*cv*(300.0-parcel.state().state().temperature);
    let cooled = parcel.add_heat(heat_out, energy_controls(), || false).unwrap();
    let net_work = compression.work_by_gas_j+expansion.work_by_gas_j;
    close(net_work, heat_in+heat_out);
    close(net_work/heat_in, 1.0-4.0_f64.powf(1.0-initial.state().gamma));
    close(parcel.state().state().temperature, initial.state().temperature);
    close(parcel.state().state().pressure, initial.state().pressure);
    assert_eq!(parcel.mass_kg(), mass);
    assert!(heated.residual_j.abs() <= mass*energy_controls().energy_tolerance_j_kg);
    assert!(cooled.residual_j.abs() <= mass*energy_controls().energy_tolerance_j_kg);
}

#[test]
fn failed_volume_operations_and_late_cancellation_preserve_all_state() {
    let mut parcel = parcel(0.001);
    let saved = parcel.clone();
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300, 1e300] {
        assert!(parcel.set_volume_isentropic(volume, entropy_controls(), || false).is_err());
        assert_eq!(parcel, saved);
    }
    let one = GasEntropySolve { maximum_iterations: 1, ..entropy_controls() };
    assert!(matches!(parcel.set_volume_isentropic(0.0025, one, || false),
        Err(GasEntropyError::NotConverged { .. })));
    assert_eq!(parcel, saved);
    let mut count = 0;
    let mut probe = saved.clone();
    probe.set_volume_isentropic(0.0025, entropy_controls(), || { count += 1; false }).unwrap();
    let mut actual = 0;
    assert!(matches!(parcel.set_volume_isentropic(0.0025, entropy_controls(), || {
        actual += 1; actual == count
    }), Err(GasEntropyError::Cancelled)));
    assert_eq!(parcel, saved);
    // Thermodynamics can reach the temperature but the final material pressure
    // is outside source support. Never publish that tentative volume.
    let mut dense = ClosedGasParcel::new(model(vec![0.75,0.25],0.001),
        point(300.0,8e6),0.01).unwrap();
    let saved_dense = dense.clone();
    assert!(dense.set_volume_isentropic(0.005, entropy_controls(), || false).is_err());
    assert_eq!(dense, saved_dense);
}

#[test]
fn unchanged_volume_is_noop_but_does_not_bypass_admission() {
    let mut parcel = parcel(0.001);
    let saved = parcel.clone();
    let report = parcel.set_volume_isentropic(0.01, entropy_controls(), || false).unwrap();
    assert_eq!(report.work_by_gas_j, 0.0);
    assert_eq!(report.iterations, 0);
    assert_eq!(parcel, saved);
    for controls in [GasEntropySolve { maximum_iterations: 0, ..entropy_controls() },
        GasEntropySolve { minimum_temperature_k: 400.0, ..entropy_controls() }] {
        assert!(parcel.set_volume_isentropic(0.01, controls, || false).is_err());
        assert_eq!(parcel, saved);
    }
    assert!(matches!(parcel.set_volume_isentropic(0.01, entropy_controls(), || true),
        Err(GasEntropyError::Cancelled)));
    assert_eq!(parcel, saved);
    let fixed = parcel.into_isochoric();
    assert_eq!(ClosedGasParcel::from_isochoric(fixed), saved);
}
