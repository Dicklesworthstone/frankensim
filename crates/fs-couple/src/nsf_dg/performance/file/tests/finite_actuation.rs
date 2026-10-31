//! G0/G3: native finite chamber work/heat retains the original material and clock.
//! These bounded synthetic chambers are numerical fixtures, not measured anatomy.
use super::*;
use crate::nsf_dg::NsfFiniteReservoirSource;

const MIXTURE: &str =
    "mixture 300 100000 mole 5 N2 0.74 O2 0.2 CO2 0.04 He 0.015 H2O 0.005";
const THERMOCHEMICAL: &str = "thermochemical 300 100000 mole 1 synthetic-alpha 1";
const FINITE: &str = "finite_reservoir 0 0.00001 0.1 0.1 1e-16 1e-9 1e-9";
const CONTROL: &str = "boundary finite-reservoir-actuation 0 0.01 1\nhold 32 -0.00001 0.25";
const PASSIVE: &str = "boundary finite-reservoir 0 1\nhold 32";
const VOLUME_RATE: f64 = -1e-5;
const HEAT_RATE: f64 = 0.25;

#[test]
fn g0_named_actuation_starter_uses_both_ducts_and_reciprocal_escape_without_source_packs() {
    let text = include_str!("../../../../../examples/resolved-dg-actuated-mixture.performance");
    let file = read(text, 16);
    let initial = file.renderer().system().finite_reservoir_inventory().unwrap();
    assert_eq!(initial.source_kind, NsfFiniteReservoirSource::NamedMixture);
    assert!(file.renderer().system().finite_reservoir().is_none());
    assert_eq!(file.renderer().system().finite_reservoir_aperture(), Some(0));
    assert_eq!(file.renderer().system().reservoir_radiation_aperture(), Some(1));
    let mut renderer = file.into_renderer();
    renderer.validate_sample_rate(48000).unwrap();
    let mut pressure = [0.0; 16];
    renderer.block(&mut pressure).unwrap();
    assert!(pressure.iter().all(|value| value.is_finite()));
    let frame = renderer.last_frame().unwrap();
    assert!(frame.upstream.is_some());
    assert_eq!(frame.reservoir_radiation.as_ref().unwrap().resolved_exchange.aperture, 1);
    let source = frame.finite_reservoir.as_ref().unwrap();
    assert_eq!(source.source_exchange.aperture, 0);
    assert!(source.final_volume_m3 < initial.volume_m3);
    assert!(source.actuation.as_ref().unwrap().step.work_by_gas_j < 0.0);
    assert!(source.energy_residual_j.abs() <= 1e-9);
}

fn build(thermochemical: bool, block: usize, edit: impl FnOnce(String) -> String)
    -> Result<ResolvedDgFilePerformance, ResolvedDgInputError>
{
    let packs = support::packs("synthetic-alpha", 0.028, 3.5, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas,
    };
    let declaration = if thermochemical { THERMOCHEMICAL } else { MIXTURE };
    let record = read_gas_record(&mut declaration.split_ascii_whitespace(), None,
        Some(&sources)).unwrap();
    let gas = *record.state();
    let density = 100002.0 / (gas.specific_gas_constant * 300.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    let text = fixture(gas, declaration)
        .replace("fluid_limits 108 250 600", "fluid_limits 108 290 320")
        .replace("coupling 64 6912", "coupling 128 13824")
        .replace("initial rest", &format!(
            "reservoir_aperture 1 low {density:.17e} 0 0 0 300 1\nreservoir_cell 5 0 1\n{FINITE}\ninitial rest"))
        .replace("observation inlet-pressure", &format!(
            "upstream_dg dg_tube 0.03 {radius:.17e} matched 1 2 32 50 2000 0.01 1e-11 512\nobservation inlet-pressure"))
        .replace("boundary 1\nhold 32 -0.1", CONTROL);
    ResolvedDgFilePerformance::from_bytes_with_sources(edit(text).as_bytes(), block,
        Arc::new(CancelGate::new_clock_free()), Some(&sources),
        |_| Err("this fixture has no material wall attachment".into()))
}

#[test]
fn g3_named_and_nasa_chambers_integrate_physical_rates_and_replay_conserved_work() {
    for thermochemical in [false, true] {
        let file = build(thermochemical, 8, |text| text).unwrap();
        let initial = file.renderer().system().finite_reservoir_inventory().unwrap();
        assert_eq!(initial.source_kind, if thermochemical {
            NsfFiniteReservoirSource::Thermochemical
        } else { NsfFiniteReservoirSource::NamedMixture });
        assert_eq!(file.renderer().system().finite_reservoir().is_some(), thermochemical);
        assert!(file.renderer().drive_history().iter().all(|drive| *drive ==
            NsfDgDrive::FiniteReservoirActuation {
                aperture: 0, volume_rate_m3_s: VOLUME_RATE, heat_power_w: HEAT_RATE,
                maximum_relative_volume_change: 0.01,
            }));
        let dt = file.renderer().system().load().propagator().dt();
        let mut renderer = file.into_renderer();
        let mut pressure = vec![0.0; SAMPLES];
        for sample in &mut pressure {
            let before = renderer.system().finite_reservoir_inventory().unwrap();
            renderer.block(core::slice::from_mut(sample)).unwrap();
            let after = renderer.system().finite_reservoir_inventory().unwrap();
            let frame = renderer.last_frame().unwrap();
            let account = frame.finite_reservoir.as_ref().unwrap();
            let motion = account.actuation.as_ref().unwrap();
            assert_eq!(account.source_exchange, frame.reservoir_exchanges[0]);
            assert_eq!(motion.step.initial_volume_m3, before.volume_m3);
            assert_eq!(motion.step.final_volume_m3, after.volume_m3);
            assert_eq!(account.volume_m3, before.volume_m3);
            assert_eq!(account.final_volume_m3, after.volume_m3);
            assert!((after.volume_m3 - before.volume_m3 - VOLUME_RATE * dt).abs()
                <= 4.0 * f64::EPSILON * before.volume_m3);
            assert_eq!(motion.step.supplied_heat_j, HEAT_RATE * dt);
            assert!(motion.step.work_by_gas_j < 0.0);
            assert!((motion.step.work_by_gas_j - motion.midpoint_pressure_pa
                * (after.volume_m3 - before.volume_m3)).abs() <= 1e-12);
            assert!((after.mass_kg - before.mass_kg
                - account.source_exchange.mass_out_kg).abs() <= 1e-16);
            assert!((account.sensible_energy_change_j - motion.step.stored_energy_change_j
                - account.source_exchange.total_energy_out_j).abs() <= 1e-9);
            let input = account.source_exchange.total_energy_out_j
                + motion.step.supplied_heat_j - motion.step.work_by_gas_j;
            assert!((after.sensible_energy_j - before.sensible_energy_j - input).abs() <= 1e-9);
            assert!((account.sensible_energy_change_j - input).abs() <= 1e-9);
            assert!((motion.step.stored_energy_change_j
                - motion.step.supplied_heat_j + motion.step.work_by_gas_j).abs() <= 1e-9);
            assert!(account.energy_residual_j.abs() <= 1e-9);
            assert!(motion.step.energy_residual_j.abs() <= 1e-9);
            assert!(account.mass_residual_kg.abs() <= 1e-16);
        }
        assert!(pressure.iter().all(|value| value.is_finite()));
        assert!(renderer.system().finite_reservoir_inventory().unwrap().volume_m3 < initial.volume_m3);
        let (passive, _) = render(build(thermochemical, 8,
            |text| text.replace(CONTROL, PASSIVE)).unwrap(), 8);
        assert_ne!(pressure, passive, "physical compression and heat must reach the resolved bore");
        let (partitioned, replay) = render(build(thermochemical, 3, |text| text).unwrap(), 3);
        assert_eq!(pressure, partitioned);
        assert_eq!(renderer.system().finite_reservoir_inventory(), replay.system().finite_reservoir_inventory());
        assert_eq!(renderer.system().zone().state_hash(), replay.system().zone().state_hash());
        assert_eq!(renderer.system().load().state(), replay.system().load().state());
        assert_eq!(renderer.system().upstream_load().unwrap().state(), replay.system().upstream_load().unwrap().state());
    }
}

#[test]
fn g0_actuation_ramps_use_original_interval_starts_and_physical_dt_after_resampling() {
    let file = build(false, 8, |text| text.replace("clock 96000 1", "clock 72000 1")
        .replace(CONTROL, "boundary finite-reservoir-actuation 0 0.01 2\nramp 16 -0.00001 0.25 -0.00002 0.5\nhold 32 0 0")).unwrap();
    assert_eq!(file.info().boundary_ramps, 1);
    let history = file.renderer().drive_history().to_vec();
    assert_eq!(history.len(), 43);
    let dt = file.renderer().system().load().propagator().dt();
    assert!(!same_clock(dt, 1.0 / 72000.0));
    let initial = file.renderer().system().finite_reservoir_inventory().unwrap();
    let mut expected_volume = initial.volume_m3;
    for (index, drive) in history.iter().enumerate() {
        let NsfDgDrive::FiniteReservoirActuation {
            aperture, volume_rate_m3_s, heat_power_w, maximum_relative_volume_change,
        } = drive else {
            panic!("native actuation lost its indexed physical rates");
        };
        assert_eq!(*aperture, 0);
        assert_eq!(*maximum_relative_volume_change, 0.01);
        let time = index as f64 * dt;
        let (volume_rate, heat_rate) = if time < 16.0 / 72000.0 {
            let fraction = time / (16.0 / 72000.0);
            (-1e-5 * (1.0 + fraction), 0.25 * (1.0 + fraction))
        } else { (0.0, 0.0) };
        assert!((*volume_rate_m3_s - volume_rate).abs() <= 64.0 * f64::EPSILON * 2e-5);
        assert!((*heat_power_w - heat_rate).abs() <= 64.0 * f64::EPSILON * 0.5);
        expected_volume += volume_rate_m3_s * dt;
    }
    let mut renderer = file.into_renderer();
    renderer.validate_sample_rate(72000).unwrap();
    let mut pressure = [0.0; SAMPLES];
    for chunk in pressure.chunks_mut(3) { renderer.block(chunk).unwrap(); }
    assert!(pressure.iter().all(|value| value.is_finite()));
    assert_eq!(renderer.system().accepted_steps(), history.len() as u64);
    assert!((renderer.system().finite_reservoir_inventory().unwrap().volume_m3 - expected_volume).abs()
        <= 64.0 * f64::EPSILON * initial.volume_m3);
}

#[test]
fn g0_named_actuation_requires_owned_bounded_controls_and_preserves_the_accepted_prefix() {
    for (from, to) in [
        ("finite_reservoir 0 0.00001 0.1 0.1 1e-16 1e-9 1e-9\n", ""),
        ("boundary finite-reservoir-actuation 0 0.01 1", "boundary finite-reservoir-actuation 0 0 1"),
        ("boundary finite-reservoir-actuation 0 0.01 1", "boundary finite-reservoir-actuation 1 0.01 1"),
        (CONTROL, "boundary aperture 0 1\nhold 32 100002 300"),
        ("hold 32 -0.00001 0.25", "hold 32 -0.00001 0.25 100002 300"),
        ("hold 32 -0.00001 0.25", "ramp 32 -0.00001 0.25 NaN 0.5"),
    ] {
        assert!(build(false, 1, |text| text.replace(from, to)).is_err(), "accepted {to}");
    }
    let snapshot = |renderer: &NsfDgPerformance<D3q27>| {
        let system = renderer.system();
        (system.finite_reservoir_inventory(), system.accepted_steps(), renderer.samples_rendered(),
            system.zone().state_hash(), system.zone().cumulative_ledger(),
            system.zone().reservoir_aperture_at(0).unwrap().clone(),
            system.zone().last_reservoir_exchanges().to_vec(), system.load().state().clone(),
            system.upstream_load().unwrap().state().clone())
    };
    // A bounded drive does not force chamber motion. Zero rates must preserve
    // the original passive source path, including both original duct states.
    let mut passive = build(false, 1, |text| text.replace(CONTROL, PASSIVE)).unwrap().into_renderer();
    let mut zero = build(false, 1, |text| text.replace("hold 32 -0.00001 0.25", "hold 32 0 0"))
        .unwrap().into_renderer();
    let mut a = [0.0];
    let mut b = [0.0];
    passive.block(&mut a).unwrap();
    zero.block(&mut b).unwrap();
    assert_eq!(a, b);
    assert_eq!(snapshot(&passive), snapshot(&zero));

    // Refuse excessive motion and a heat increment that leaves the explicitly
    // retained material window, after a valid compression/heating prefix.
    for next in ["-1 0.25", "0 1000"] {
        let mut renderer = build(false, 1, |text| text
            .replace("fluid_limits 108 290 320", "fluid_limits 108 299.9 300.1")
            .replace(CONTROL, &format!(
                "boundary finite-reservoir-actuation 0 0.01 2\nhold 1 -0.00001 0.25\nhold 32 {next}")))
            .unwrap().into_renderer();
        renderer.block(&mut [0.0]).unwrap();
        let prefix = snapshot(&renderer);
        assert!(renderer.block(&mut [0.0]).is_err());
        assert_eq!(snapshot(&renderer), prefix);
    }
}
