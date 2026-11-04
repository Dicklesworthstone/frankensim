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

fn viscous_liquid() -> ClassicalBulkAbsorption {
    // Cp=Cv decouples temperature exactly. This manufactured fluid has the
    // independent longitudinal solution k=(omega/c)/sqrt(1-i*omega*nu/c²).
    let state = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1200.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, 5.0).unwrap(),
        4000.0,
    )
    .unwrap();
    ClassicalBulkAbsorption::from_material(&state, 50.0).unwrap()
}

fn viscous_wave(omega: f64) -> [f64; 2] {
    let x = omega * (4.0 * 5.0 / 3.0 + 50.0) / (1000.0 * 1200.0_f64.powi(2));
    let norm = 1.0_f64.hypot(x);
    let re = ((norm + 1.0) / (2.0 * norm * norm)).sqrt();
    let im = x / (2.0 * norm * norm * re);
    [omega / 1200.0 * re, omega / 1200.0 * im]
}

fn complex_transfer(taps: &[f64], omega_dt: f64) -> [f64; 2] {
    taps.iter().enumerate().fold([0.0; 2], |[re, im], (n, h)| {
        [
            re + h * (n as f64 * omega_dt).cos(),
            im - h * (n as f64 * omega_dt).sin(),
        ]
    })
}

#[test]
fn g1_thermoviscous_filter_retains_independent_complex_viscous_phase() {
    let source = viscous_liquid();
    let (dt, range) = (1.0 / 48_000.25, 20.0);
    let filter = FluidObserverFilter::with_model(
        source, range, dt, 1e-12, FluidPropagationModel::ThermoviscousDispersion,
    )
    .unwrap();
    let receipt = filter.receipt();
    assert_eq!(receipt.model, FluidPropagationModel::ThermoviscousDispersion);
    assert_eq!(receipt.sample_rate_hz, 1.0 / dt);
    assert!(receipt.truncated_energy <= 1e-12);
    assert!((receipt.receiver_delay_s + receipt.low_frequency_delay_s - range / 1200.0).abs() < 1e-16);
    for frequency in [0.0, 1500.0, 6000.0, 12_000.0, 18_000.0] {
        let omega = core::f64::consts::TAU * frequency;
        let [wave_re, wave_im] = viscous_wave(omega);
        let magnitude = (-wave_im * range).exp();
        let phase = -(wave_re - omega / 1200.0) * range;
        let expected = [magnitude * phase.cos(), magnitude * phase.sin()];
        let direct = thermoviscous_transfer(&source, omega, range).unwrap();
        assert!((direct.re - expected[0]).hypot(direct.im - expected[1]) < 2e-12);
        let [re, im] = complex_transfer(filter.taps(), omega * dt);
        // Undo ONLY reported numerical alignment/padding, not the physical
        // excess phase. A magnitude-only or wrong-conjugation filter fails.
        let (sin, cos) = (omega * receipt.latency_s).sin_cos();
        let actual = [re * cos - im * sin, re * sin + im * cos];
        assert!(
            (actual[0] - expected[0]).hypot(actual[1] - expected[1]) < 8e-5,
            "{frequency}: {actual:?} != {expected:?}"
        );
        if frequency > 0.0 {
            assert!(expected[1] > 0.0, "this viscous fixture has positive excess FFT phase");
        }
    }
    let legacy = FluidObserverFilter::new(source, range, dt, 1e-12).unwrap();
    let explicit = FluidObserverFilter::with_model(
        source, range, dt, 1e-12, FluidPropagationModel::ClassicalMinimumPhase,
    )
    .unwrap();
    assert_eq!(legacy.taps(), explicit.taps());
    assert_eq!(legacy.receipt(), explicit.receipt());
}

#[test]
fn g3_selected_models_share_batch_and_streaming_without_circular_wrap() {
    for model in [
        FluidPropagationModel::ClassicalMinimumPhase,
        FluidPropagationModel::ThermoviscousDispersion,
    ] {
        let (dt, range) = (1.0 / 48_000.25, 20.0);
        let mut filter = FluidObserverFilter::with_model(
            viscous_liquid(), range, dt, 1e-12, model,
        )
        .unwrap();
        let count = 3 * filter.taps().len() + 37;
        let original: Vec<f64> = (0..count)
            .map(|n| (0.31 * n as f64).sin() + if n == count - 1 { 2.0 } else { 0.0 })
            .collect();
        let streaming: Vec<f64> = original.iter()
            .map(|&p| filter.process(p).unwrap())
            .collect();
        let mut batch = original;
        let receipt = propagate_pressure_history_with_model(
            &mut batch, dt, range, viscous_liquid(), 1e-12, model,
        )
        .unwrap();
        assert_eq!(&receipt, filter.receipt());
        for (actual, expected) in batch.iter().zip(&streaming) {
            assert!((actual - expected).abs() < 4e-12, "{model:?}");
        }
        let mut late = vec![0.0; count];
        late[count - 1] = 1.0;
        propagate_pressure_history_with_model(
            &mut late, dt, range, viscous_liquid(), 1e-12, model,
        )
        .unwrap();
        assert!(late[..count - 1].iter().all(|value| value.abs() < 2e-14));
        assert!((late[count - 1] - filter.taps()[0]).abs() < 2e-14);
    }
}

#[test]
fn g3_thermoviscous_physical_path_transfer_composes() {
    let source = viscous_liquid();
    for frequency in [0.0, 6000.0, 18_000.0] {
        let omega = core::f64::consts::TAU * frequency;
        let a = thermoviscous_transfer(&source, omega, 1.5).unwrap();
        let b = thermoviscous_transfer(&source, omega, 2.25).unwrap();
        let whole = thermoviscous_transfer(&source, omega, 3.75).unwrap();
        assert!((whole.re - (a.re * b.re - a.im * b.im))
            .hypot(whole.im - (a.re * b.im + a.im * b.re)) < 1e-14);
    }
}

#[test]
fn g0_thermoviscous_admission_is_atomic_and_zero_range_is_exact() {
    let model = FluidPropagationModel::ThermoviscousDispersion;
    let original = vec![1.0, -2.0, 0.5];
    for (dt, range, budget) in [
        (0.0, 20.0, 1e-12),
        (f64::NAN, 20.0, 1e-12),
        (1e-15, 20.0, 1e-12),
        (1e200, 20.0, 1e-12), // nonzero acoustic loss underflows: refuse
        (1.0 / 48_000.0, -1.0, 1e-12),
        (1.0 / 48_000.0, f64::INFINITY, 1e-12),
        (1.0 / 48_000.0, 20.0, 0.0),
        (1.0 / 48_000.0, 20.0, 1.0),
        (1.0 / 48_000.0, 20.0, f64::NAN),
    ] {
        let mut pressure = original.clone();
        assert!(propagate_pressure_history_with_model(
            &mut pressure, dt, range, viscous_liquid(), budget, model,
        ).is_err());
        assert_eq!(pressure, original);
    }
    let mut pressure = original.clone();
    let receipt = propagate_pressure_history_with_model(
        &mut pressure, 1.0 / 48_000.0, 0.0, viscous_liquid(), 1e-12, model,
    )
    .unwrap();
    assert_eq!(pressure, original);
    assert_eq!(receipt.model, model);
    assert_eq!(receipt.taps, 1);
    assert_eq!(receipt.receiver_delay_s, 0.0);
    let mut filter = FluidObserverFilter::with_model(
        viscous_liquid(), 20.0, 1.0 / 48_000.0, 1e-12, model,
    )
    .unwrap();
    filter.process(0.75).unwrap();
    let mut unchanged = filter.clone();
    assert!(filter.process(f64::NAN).is_err());
    for sample in [1.0, -1.0, 0.0, 0.2] {
        assert_eq!(filter.process(sample).unwrap(), unchanged.process(sample).unwrap());
    }
}
