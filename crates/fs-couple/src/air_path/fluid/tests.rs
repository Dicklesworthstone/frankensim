//! Independent attenuation formula, exact-clock convolution and refusal tests.
use super::*;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;

fn liquid(bulk: f64) -> ClassicalBulkAbsorption {
    // Explicit authored Newtonian-liquid fixture, not measured material data.
    let state = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1200.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.05).unwrap(),
        3900.0,
    )
    .unwrap();
    ClassicalBulkAbsorption::from_material(&state, bulk).unwrap()
}

fn magnitude(taps: &[f64], omega_dt: f64) -> f64 {
    let (re, im) = taps
        .iter()
        .enumerate()
        .fold((0.0, 0.0), |(re, im), (i, tap)| {
            (
                re + tap * (i as f64 * omega_dt).cos(),
                im - tap * (i as f64 * omega_dt).sin(),
            )
        });
    re.hypot(im)
}

#[test]
fn g1_liquid_fir_matches_independent_pressure_attenuation_and_physical_clock() {
    let dt = 1.0 / 48_000.25;
    let range = 500.0;
    let filter = FluidObserverFilter::new(liquid(0.1), range, dt, 1e-12).unwrap();
    let report = filter.receipt();
    assert_eq!(report.sample_rate_hz, 1.0 / dt);
    assert_eq!(report.absorption, liquid(0.1));
    assert!(report.truncated_energy <= report.maximum_truncated_energy);
    assert!(
        (report.receiver_delay_s + report.low_frequency_delay_s - range / 1200.0).abs() < 1e-16
    );
    for frequency in [0.0, 1000.0, 6000.0, 12_000.0, 20_000.0] {
        let omega = core::f64::consts::TAU * frequency;
        // Direct independent constitutive expression, not a call back into
        // the material law or the minimum-phase design being tested.
        let alpha = omega * omega / (2.0 * 1000.0 * 1200.0_f64.powi(3))
            * (4.0 * 0.05 / 3.0 + 0.1 + 0.6 * (1.0 / 3900.0 - 1.0 / 4000.0));
        let expected = (-alpha * range).exp();
        let actual = magnitude(filter.taps(), omega * dt);
        assert!(
            (actual - expected).abs() < 8e-5,
            "{frequency}: {actual} vs {expected}"
        );
    }
}

#[test]
fn g3_streaming_and_batch_share_the_exact_filter_without_circular_wrap() {
    let dt = 1.0 / 48_000.25;
    let mut filter = FluidObserverFilter::new(liquid(0.1), 500.0, dt, 1e-12).unwrap();
    let count = 3 * filter.taps().len() + 37;
    let original: Vec<f64> = (0..count)
        .map(|i| (0.7 * i as f64).sin() + if i == count - 1 { 2.0 } else { 0.0 })
        .collect();
    let streaming: Vec<f64> = original
        .iter()
        .map(|&p| filter.process(p).unwrap())
        .collect();
    let mut batch = original;
    let report = propagate_pressure_history(&mut batch, dt, 500.0, liquid(0.1), 1e-12).unwrap();
    assert_eq!(&report, filter.receipt());
    for (actual, expected) in batch.iter().zip(&streaming) {
        assert!((actual - expected).abs() < 4e-12);
    }
    let mut late = vec![0.0; count];
    late[count - 1] = 1.0;
    propagate_pressure_history(&mut late, dt, 500.0, liquid(0.1), 1e-12).unwrap();
    assert!(late[..count - 1].iter().all(|value| value.abs() < 2e-14));
    assert!((late[count - 1] - filter.taps()[0]).abs() < 2e-14);
}

#[test]
fn g3_bulk_viscosity_changes_absorption_without_changing_density_or_sound_speed() {
    let low = FluidObserverFilter::new(liquid(0.0), 500.0, 1.0 / 48_000.0, 1e-12).unwrap();
    let high = FluidObserverFilter::new(liquid(0.1), 500.0, 1.0 / 48_000.0, 1e-12).unwrap();
    assert_eq!(
        low.receipt().absorption.state(),
        high.receipt().absorption.state()
    );
    let omega = core::f64::consts::TAU * 20_000.0;
    let ratio = magnitude(high.taps(), omega / 48_000.0) / magnitude(low.taps(), omega / 48_000.0);
    let expected = (-omega * omega * 0.1 * 500.0 / (2.0 * 1000.0 * 1200.0_f64.powi(3))).exp();
    assert!((ratio - expected).abs() < 8e-5);
    assert!(ratio < 0.8);
}

#[test]
fn g0_invalid_or_strong_loss_refusal_is_atomic_and_zero_range_is_exact() {
    let original = vec![1.0, -2.0, 0.5];
    for (dt, range, budget) in [
        (0.0, 1.0, 1e-9),
        (f64::NAN, 1.0, 1e-9),
        (1e-15, 1.0, 1e-9),
        (1.0 / 48_000.0, -1.0, 1e-9),
        (1.0 / 48_000.0, 1.0, 0.0),
        (1.0 / 48_000.0, 1.0, 1.0),
    ] {
        let mut pressure = original.clone();
        assert!(propagate_pressure_history(&mut pressure, dt, range, liquid(0.1), budget).is_err());
        assert_eq!(pressure, original);
    }
    let mut pressure = original.clone();
    let report =
        propagate_pressure_history(&mut pressure, 1.0 / 48_000.0, 0.0, liquid(0.1), 1e-9).unwrap();
    assert_eq!(pressure, original);
    assert_eq!(report.receiver_delay_s, 0.0);
    assert_eq!(report.taps, 1);
    let mut filter = FluidObserverFilter::new(liquid(0.1), 500.0, 1.0 / 48_000.0, 1e-9).unwrap();
    filter.process(0.75).unwrap();
    let mut unchanged = filter.clone();
    assert!(filter.process(f64::NAN).is_err());
    for sample in [1.0, -1.0, 0.0, 0.2] {
        assert_eq!(
            filter.process(sample).unwrap(),
            unchanged.process(sample).unwrap()
        );
    }
}
