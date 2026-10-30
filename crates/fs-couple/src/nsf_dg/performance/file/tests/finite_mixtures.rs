//! G0/G3: original named mixture storage feeds the resolved source interval.
//! Frozen composition and small numerical chambers are not measured anatomy.
use super::*;
use crate::nsf_dg::NsfFiniteReservoirSource;

const MIXTURE: &str =
    "mixture 300 100000 mole 5 N2 0.74 O2 0.2 CO2 0.04 He 0.015 H2O 0.005";
const VOLUME: f64 = 1e-5;
const CONTROL: &str = "boundary finite-reservoir 0 1\nhold 32";

fn finite_row(volume: f64) -> String {
    format!("finite_reservoir 0 {volume:.17e} 0.1 0.1 1e-16 1e-9 1e-9")
}

fn input(volume: f64) -> String {
    let reference = read_gas_record(&mut MIXTURE.split_ascii_whitespace(), None, None).unwrap();
    let source = MIXTURE.replace("300 100000", "300 100002");
    let source = read_gas_record(&mut source.split_ascii_whitespace(), None, None).unwrap();
    let gas = *reference.state();
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    fixture(gas, MIXTURE)
        .replace("fluid_limits 108 250 600", "fluid_limits 108 290 320")
        .replace("coupling 64 6912", "coupling 128 13824")
        .replace("initial rest", &format!(
            "reservoir_aperture 1 low {:.17e} 0 0 0 300 1\nreservoir_cell 5 0 1\n{}\ninitial rest",
            source.state().density, finite_row(volume)))
        .replace("observation inlet-pressure", &format!(
            "upstream_dg dg_tube 0.03 {radius:.17e} matched 1 2 32 50 2000 0.01 1e-11 512\nobservation inlet-pressure"))
        .replace("boundary 1\nhold 32 -0.1", CONTROL)
}

#[test]
fn g3_native_named_inventory_conserves_mass_and_energy_and_changes_next_source_pressure() {
    let text = input(VOLUME);
    let file = read(&text, 8);
    assert!(file.renderer().system().finite_reservoir().is_none());
    assert_eq!(file.renderer().system().finite_reservoir_aperture(), Some(0));
    assert!(file.renderer().drive_history().iter().all(|drive|
        *drive == NsfDgDrive::FiniteReservoir { aperture: 0 }));
    let initial = file.renderer().system().finite_reservoir_inventory().unwrap();
    assert_eq!(initial.source_kind, NsfFiniteReservoirSource::NamedMixture);
    assert_eq!(initial.aperture, 0);
    assert_eq!(initial.volume_m3, VOLUME);
    assert_eq!(initial.energy_reference_j_kg, 0.0);
    // Prescribe the exact reconstructed initial state, including its source
    // roundoff, to isolate finite inventory feedback from initial conditions.
    let prescribed = text.replace(&format!("{}\n", finite_row(VOLUME)), "")
        .replace(CONTROL, &format!("boundary aperture 0 1\nhold 32 {:.17e} {:.17e}",
            initial.gas.pressure, initial.gas.temperature));
    let mut held = read(&prescribed, 8).into_renderer();
    let mut finite = file.into_renderer();
    let mut pressure = vec![0.0; SAMPLES];
    let mut held_pressure = vec![0.0; SAMPLES];
    for index in 0..SAMPLES {
        let before = finite.system().finite_reservoir_inventory().unwrap();
        finite.block(&mut pressure[index..index + 1]).unwrap();
        held.block(&mut held_pressure[index..index + 1]).unwrap();
        let after = finite.system().finite_reservoir_inventory().unwrap();
        let frame = finite.last_frame().unwrap();
        let account = frame.finite_reservoir.as_ref().unwrap();
        assert_eq!(account.source_exchange, frame.reservoir_exchanges[0]);
        assert_eq!(account.initial_mass_kg, before.mass_kg);
        assert_eq!(account.final_mass_kg, after.mass_kg);
        assert_eq!(account.initial_pressure_pa, before.gas.pressure);
        assert_eq!(account.final_pressure_pa, after.gas.pressure);
        assert!((after.mass_kg - before.mass_kg - account.source_exchange.mass_out_kg).abs() <= 1e-16);
        assert!((after.sensible_energy_j - before.sensible_energy_j
            - account.source_exchange.total_energy_out_j).abs() <= 1e-9);
        assert!(account.mass_residual_kg.abs() <= 1e-16);
        assert!(account.energy_residual_j.abs() <= 1e-9);
        let cv = after.gas.specific_heat_cp - after.gas.specific_gas_constant;
        let storage = after.mass_kg * cv * after.gas.temperature;
        assert!((after.sensible_energy_j - storage).abs()
            <= 64.0 * f64::EPSILON * storage.abs());
        let NsfGasMedium::NamedMixture(source) = finite.system().zone().medium() else {
            panic!("finite named chamber lost the original mixture owner");
        };
        assert_eq!(after.gas, source.at(after.gas.temperature, after.gas.pressure).unwrap());
        if index == 0 {
            assert!(after.mass_kg < before.mass_kg, "the pressurized chamber must supply mass");
            assert_eq!(finite.system().zone().state_hash(), held.system().zone().state_hash());
            assert_eq!(finite.system().load().state(), held.system().load().state());
            assert_eq!(finite.system().upstream_load().unwrap().state(), held.system().upstream_load().unwrap().state());
        }
    }
    assert_ne!(pressure, held_pressure);
    assert_ne!(finite.system().finite_reservoir_inventory().unwrap().gas.pressure, initial.gas.pressure);
    let (partitioned, replay) = render(read(&text, 3), 3);
    assert_eq!(pressure, partitioned);
    assert_eq!(finite.system().finite_reservoir_inventory(), replay.system().finite_reservoir_inventory());
    assert_eq!(finite.system().zone().state_hash(), replay.system().zone().state_hash());
    assert_eq!(finite.system().load().state(), replay.system().load().state());
    assert_eq!(finite.system().upstream_load().unwrap().state(), replay.system().upstream_load().unwrap().state());
}

#[test]
fn g0_native_named_inventory_refusal_and_cancel_preserve_every_accepted_owner() {
    let snapshot = |renderer: &NsfDgPerformance<D3q27>| {
        let system = renderer.system();
        (system.finite_reservoir_inventory(), system.accepted_steps(), renderer.samples_rendered(),
            system.zone().state_hash(), system.zone().cumulative_ledger(),
            system.zone().reservoir_aperture_at(0).unwrap().clone(),
            system.zone().last_reservoir_exchanges().to_vec(), system.load().state().clone(),
            system.upstream_load().unwrap().state().clone())
    };
    // This source is initially valid, but a complete original transfer cannot
    // fit its explicitly tiny inventory. No partial fluid step may publish.
    let mut refused = read(&input(1e-18), 1).into_renderer();
    let before = snapshot(&refused);
    assert!(refused.block(&mut [0.0]).is_err());
    assert_eq!(snapshot(&refused), before);

    let text = input(VOLUME);
    let gate = Arc::new(CancelGate::new_clock_free());
    let mut accepted = ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 1,
        Arc::clone(&gate)).unwrap().into_renderer();
    accepted.block(&mut [0.0]).unwrap();
    let prefix = snapshot(&accepted);
    gate.request();
    assert!(accepted.block(&mut [0.0]).is_err());
    assert_eq!(snapshot(&accepted), prefix);
    // Neither an authored pressure/temperature nor an out-of-domain material
    // declaration may override the finite named source's original authority.
    for invalid in [
        text.replace(CONTROL, "boundary aperture 0 1\nhold 32 100002 300"),
        text.replace("fluid_limits 108 290 320", "fluid_limits 108 250 600"),
        text.replace("N2 0.74", "N2 0.69").replace("H2O 0.005", "H2O 0.055"),
    ] {
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 1,
            Arc::new(CancelGate::new_clock_free())).is_err());
    }
}
