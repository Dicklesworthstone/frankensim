//! Native rates drive original finite gas work/heat on the physical clock.
use super::*;

const ACTUATED: &str =
    "boundary finite-reservoir-actuation 0 0.05 2\nramp 16 -2e-7 0.01 2e-7 -0.01\nhold 32 0 0";

fn ramp(text: String) -> String {
    text.replace("boundary finite-reservoir 0 1\nhold 32", ACTUATED)
}

#[test]
fn native_actuation_ramp_replays_across_output_block_partitions() {
    let file = build(false, 8, ramp).unwrap();
    assert_eq!(file.info().boundary_segments, 2);
    assert_eq!(file.info().boundary_ramps, 1);
    let history = file.renderer().drive_history();
    assert_eq!(history.len(), SAMPLES);
    for (i, control) in history.iter().enumerate() {
        let NsfDgDrive::FiniteReservoirActuation {
            aperture, volume_rate_m3_s, heat_power_w, maximum_relative_volume_change,
        } = control else { panic!("finite source rate units lost"); };
        assert_eq!(*aperture, 0);
        assert_eq!(*maximum_relative_volume_change, 0.05);
        let fraction = i as f64 / 16.0;
        let (rate, power) = if i < 16 {
            (-2e-7 + 4e-7 * fraction, 0.01 - 0.02 * fraction)
        } else { (0.0, 0.0) };
        assert!((volume_rate_m3_s - rate).abs() < 1e-21);
        assert!((heat_power_w - power).abs() < 1e-16);
    }
    let (whole, a) = render(file, 8);
    let (parts, b) = render(build(false, 3, ramp).unwrap(), 3);
    assert_eq!(whole, parts);
    assert_eq!(a.system().zone().state_hash(), b.system().zone().state_hash());
    assert_eq!(a.system().load().state(), b.system().load().state());
    assert_eq!(a.system().finite_reservoir(), b.system().finite_reservoir());
    assert!(a.last_frame().unwrap().finite_reservoir.as_ref().unwrap().actuation.is_some());
}

#[test]
fn native_actuation_ramp_samples_physical_interval_starts_with_an_independent_output_clock() {
    let file = build(false, 8, |text| ramp(text).replace("clock 96000 1", "clock 72000 1")).unwrap();
    let history = file.renderer().drive_history();
    assert_eq!(history.len(), 43);
    let dt = file.renderer().system().load().propagator().dt();
    assert!(!same_clock(dt, 1.0 / 72000.0));
    for (i, control) in history.iter().enumerate() {
        let NsfDgDrive::FiniteReservoirActuation {
            volume_rate_m3_s, heat_power_w, ..
        } = control else { panic!("finite source rate units lost"); };
        let time = i as f64 * dt;
        let (rate, power) = if time < 16.0 / 72000.0 {
            let fraction = time / (16.0 / 72000.0);
            (-2e-7 + 4e-7 * fraction, 0.01 - 0.02 * fraction)
        } else { (0.0, 0.0) };
        assert!((volume_rate_m3_s - rate).abs() < 1e-21);
        assert!((heat_power_w - power).abs() < 1e-16);
    }
}

#[test]
fn native_actuation_keeps_passive_escape_both_ducts_and_original_physical_rate_units() {
    let mut renderer = build(true, 4, |text| text.replace(
        "boundary finite-reservoir 0 1\nhold 64",
        "boundary finite-reservoir-actuation 0 0.05 1\nhold 64 -2e-7 0.02",
    )).unwrap().into_renderer();
    renderer.validate_sample_rate(8000).unwrap();
    let dt = renderer.system().load().propagator().dt();
    assert!(!same_clock(dt, 1.0 / 8000.0));
    renderer.block(&mut [0.0]).unwrap();
    let frame = renderer.last_frame().unwrap();
    let account = frame.finite_reservoir.as_ref().unwrap();
    let motion = account.actuation.as_ref().unwrap();
    assert_eq!(account.source_exchange.aperture, 0);
    assert_eq!(frame.reservoir_radiation.as_ref().unwrap().resolved_exchange.aperture, 1);
    assert!(frame.upstream.is_some());
    assert_eq!(motion.step.supplied_heat_j, 0.02 * dt);
    assert!((account.final_volume_m3 - account.volume_m3 + 2e-7 * dt).abs()
        <= 4.0 * f64::EPSILON * account.volume_m3);
    assert_eq!(account.final_volume_m3, renderer.system().finite_reservoir().unwrap().volume_m3());
    assert!(motion.step.work_by_gas_j < 0.0);
    assert!(account.energy_residual_j.abs() <= 1e-9);
}

#[test]
fn native_actuation_requires_valid_owned_inventory_even_at_zero_rates() {
    for (from, to) in [
        ("0 0.05 2", "0 0 2"),
        ("0 0.05 2", "0 1 2"),
        ("0 0.05 2", "0 NaN 2"),
        ("0 0.05 2", "1 0.05 2"),
        ("ramp 16 -2e-7 0.01 2e-7 -0.01", "ramp 16 -2e-7 0.01 2e-7"),
        ("ramp 16 -2e-7 0.01 2e-7 -0.01", "ramp 16 -2e-7 0.01 NaN -0.01"),
        ("hold 32 0 0", "hold 32 0"),
        ("hold 32 0 0", "hold 32 0 0 extra"),
    ] {
        assert!(build(false, 4, |text| ramp(text).replace(from, to)).is_err(), "{from} -> {to}");
    }
    assert!(build(false, 4, |text| ramp(text).replace(&format!("{FINITE}\n"), "")).is_err());
    let named = |text: String| text.replace(DECLARATION, "mixture 300 100000 mole 1 N2 1")
        .replace("energy_solve 1e-9 96\n", "")
        .replace("fluid_limits 108 250 600", "fluid_limits 108 290 320");
    let mut passive = build(false, 4, named).unwrap().into_renderer();
    assert_eq!(passive.system().finite_reservoir_aperture(), Some(0));
    passive.validate_sample_rate(RATE).unwrap();
    passive.block(&mut [0.0]).unwrap();
    assert!(passive.last_frame().unwrap().finite_reservoir.as_ref().unwrap().actuation.is_none());
    let mut zero = build(false, 4, |text| named(text).replace(
        "boundary finite-reservoir 0 1\nhold 32",
        "boundary finite-reservoir-actuation 0 0.05 1\nhold 32 0 0",
    )).unwrap().into_renderer();
    zero.validate_sample_rate(RATE).unwrap();
    zero.block(&mut [0.0]).unwrap();
    assert_eq!(zero.system().zone().state_hash(), passive.system().zone().state_hash());
    assert_eq!(zero.system().load().state(), passive.system().load().state());
    assert_eq!(zero.system().finite_reservoir_inventory(), passive.system().finite_reservoir_inventory());
    let receipt = zero.last_frame().unwrap().finite_reservoir.as_ref().unwrap()
        .actuation.as_ref().unwrap();
    assert_eq!(receipt.step.supplied_heat_j, 0.0);
    assert_eq!(receipt.step.work_by_gas_j, 0.0);
    assert_eq!(receipt.step.stored_energy_change_j, 0.0);
}
