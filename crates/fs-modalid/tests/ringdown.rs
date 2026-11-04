use core::f64::consts::TAU;
use fs_fft::harmonic::{
    DampedToneBounds, DampedTonePole, DampedToneRefinementBudget, DampedToneRefinementError,
    HarmonicFitBudget, HarmonicFitError, refine_damped_tone_poles_channels,
};
use fs_math::c64::C64;
use fs_modalid::{
    ModalIdError, Mode, ModeCompareOptions, RingdownError, compare_mode_tables, identify_ringdown, mac,
};

// Independent free motion of m*q'' + c*q' + k*q = 0 from physical parameters
// and initial conditions, rather than sampled coefficients from a fit plan.
#[derive(Clone, Copy)]
struct Oscillator {
    mass: f64,
    stiffness: f64,
    damping: f64,
    q0: f64,
    v0: f64,
}

impl Oscillator {
    fn new(mass: f64, natural_hz: f64, zeta: f64, q0: f64, v0: f64) -> Self {
        Self {
            mass,
            stiffness: mass * (TAU * natural_hz).powi(2),
            damping: 2.0 * mass * TAU * natural_hz * zeta,
            q0,
            v0,
        }
    }

    fn pole(self, window_rate: f64) -> DampedTonePole {
        let alpha = self.damping / (2.0 * self.mass);
        DampedTonePole {
            frequency_hz: (self.stiffness / self.mass - alpha * alpha).sqrt() / TAU,
            decay_rate_per_s: alpha + window_rate,
        }
    }

    fn response(self, time: f64) -> f64 {
        let alpha = self.damping / (2.0 * self.mass);
        let wd = (self.stiffness / self.mass - alpha * alpha).sqrt();
        let (sin, cos) = (wd * time).sin_cos();
        (-alpha * time).exp() * (self.q0 * cos + (self.v0 + alpha * self.q0) / wd * sin)
    }

    fn residue(self) -> C64 {
        let alpha = self.damping / (2.0 * self.mass);
        let wd = (self.stiffness / self.mass - alpha * alpha).sqrt();
        C64::new(0.5 * self.q0, -0.5 * (self.v0 + alpha * self.q0) / wd)
    }
}

fn budget() -> DampedToneRefinementBudget {
    DampedToneRefinementBudget {
        absolute_frequency_tolerance_hz: 1e-8,
        absolute_decay_tolerance_per_s: 1e-7,
        max_evaluations: 60,
        min_pole_reciprocal_condition: 1e-8,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 64 * 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        },
    }
}

fn bounds(poles: &[DampedTonePole]) -> Vec<DampedToneBounds> {
    poles
        .iter()
        .map(|p| DampedToneBounds {
            frequency_hz: [p.frequency_hz - 0.4, p.frequency_hz + 0.4],
            decay_rate_per_s: [(p.decay_rate_per_s - 2.0).max(0.0), p.decay_rate_per_s + 2.0],
        })
        .collect()
}

#[test]
fn oscillator_ringdowns_feed_existing_mac_and_calibration_with_nodal_signed_sensors() {
    let oscillators = [
        Oscillator::new(0.3, 75.3, 0.009, 0.8, 10.0),
        Oscillator::new(0.05, 137.2, 0.012, -0.35, 250.0),
    ];
    // Sensor 0 is at the first mode's node; sensors 1/2 have opposite polarity
    // and cancel that mode under averaging. Every channel must remain separate.
    let shapes = [[0.0, 1.0, -1.0, 0.3], [1.0, 0.0, 0.0, -0.8]];
    let offsets = [-0.7, 2.0, -0.4, 0.9];
    let dt = 1.0 / 2048.0;
    let records: Vec<Vec<_>> = (0..4)
        .map(|c| {
            (0..769)
                .map(|i| {
                    offsets[c]
                        + (0..2)
                            .map(|k| shapes[k][c] * oscillators[k].response(i as f64 * dt))
                            .sum::<f64>()
                })
                .collect()
        })
        .collect();
    let channels: Vec<_> = records.iter().map(Vec::as_slice).collect();
    let truth = [oscillators[1].pole(0.0), oscillators[0].pole(0.0)];
    let seeds = truth.map(|p| DampedTonePole {
        frequency_hz: p.frequency_hz + 0.15,
        decay_rate_per_s: p.decay_rate_per_s - 0.4,
    });
    let result = identify_ringdown(&channels, dt, &seeds, &bounds(&truth), budget(), None).unwrap();
    assert_eq!(result.source_pole_indices, [1, 0]);
    assert!(result.fit.evaluations > 1, "the seeds need actual refinement");
    let reference: Vec<_> = [75.3, 137.2]
        .into_iter()
        .zip([0.009, 0.012])
        .enumerate()
        .map(|(k, (frequency_hz, damping_ratio))| Mode {
            frequency_hz,
            damping_ratio,
            shape: shapes[k].map(|v| C64::new(v, 0.0)).to_vec(),
            frequency_ci_hz: f64::NAN,
            damping_ci: f64::NAN,
        })
        .collect();
    for (k, mode) in result.modes.iter().enumerate() {
        assert!((mode.frequency_hz - reference[k].frequency_hz).abs() < 2e-7);
        assert!((mode.damping_ratio - reference[k].damping_ratio).abs() < 2e-9);
        assert!(mode.frequency_ci_hz.is_nan() && mode.damping_ci.is_nan());
        assert!(mac(&mode.shape, &reference[k].shape) > 1.0 - 1e-12);
        for (c, residue) in mode.shape.iter().enumerate() {
            assert!((*residue - oscillators[k].residue().scale(shapes[k][c])).abs() < 2e-8);
        }
        let source = result.source_pole_indices[k];
        assert!((result.fit.poles[source].frequency_hz - truth[source].frequency_hz).abs() < 2e-7);
    }
    for (fit, offset) in result.fit.fits.iter().zip(offsets) {
        assert!((fit.dc - offset).abs() < 2e-9);
        assert!(fit.residual_rms < 2e-9);
    }
    let comparison = compare_mode_tables(
        &reference,
        &result.modes,
        &ModeCompareOptions { mac_floor: 0.999, max_cents: 0.001 },
    )
    .unwrap();
    assert!(comparison.unmatched_reference.is_empty() && comparison.unmatched_candidate.is_empty());
    assert_eq!(comparison.pairs.len(), 2);
    assert!(comparison.pairs.iter().all(|p| p.within_tolerance && p.mac.unwrap() > 0.999999));
}

#[test]
fn modal_frequency_is_natural_not_the_damped_oscillation_frequency() {
    let oscillator = Oscillator::new(0.2, 12.5, 0.6, 0.2, 35.0);
    let dt = 1.0 / 1024.0;
    let record: Vec<_> = (0..257).map(|i| oscillator.response(i as f64 * dt)).collect();
    let truth = oscillator.pole(0.0);
    let seed = DampedTonePole {
        frequency_hz: truth.frequency_hz + 0.12,
        decay_rate_per_s: truth.decay_rate_per_s - 0.4,
    };
    let result = identify_ringdown(&[&record], dt, &[seed], &bounds(&[truth]), budget(), None).unwrap();
    assert!((result.fit.poles[0].frequency_hz - 10.0).abs() < 2e-7);
    assert!((result.modes[0].frequency_hz - 12.5).abs() < 2e-7);
    assert!((result.modes[0].damping_ratio - 0.6).abs() < 2e-8);
    assert!((result.modes[0].shape[0] - oscillator.residue()).abs() < 2e-8);
}

#[test]
fn window_correction_reorders_natural_frequencies_without_losing_pole_or_shape_identity() {
    let oscillators = [
        Oscillator::new(0.3, 10.0, 0.8, 0.7, 0.0),
        Oscillator::new(0.2, 11.0, 0.02, -0.3, 3.0),
    ];
    let rate = 30.0;
    let dt = 1.0 / 2048.0;
    let records: Vec<Vec<_>> = oscillators
        .iter()
        .map(|oscillator| {
            (0..513)
                .map(|i| {
                    let t = i as f64 * dt;
                    oscillator.response(t) * (-rate * t).exp()
                })
                .collect()
        })
        .collect();
    let channels: Vec<_> = records.iter().map(Vec::as_slice).collect();
    let measured = oscillators.map(|o| o.pole(rate));
    let bounds = bounds(&measured);
    let raw = identify_ringdown(&channels, dt, &measured, &bounds, budget(), None).unwrap();
    let corrected =
        identify_ringdown(&channels, dt, &measured, &bounds, budget(), Some(1.0 / rate)).unwrap();
    assert_eq!(raw.source_pole_indices, [1, 0]);
    assert_eq!(corrected.source_pole_indices, [0, 1]);
    assert_eq!(raw.fit, corrected.fit, "window correction must retain the measured fit");
    for (k, mode) in corrected.modes.iter().enumerate() {
        assert!((mode.frequency_hz - [10.0, 11.0][k]).abs() < 2e-7);
        assert!((mode.damping_ratio - [0.8, 0.02][k]).abs() < 2e-8);
        assert!((mode.shape[k] - oscillators[k].residue()).abs() < 2e-8);
        assert!(mode.shape[1 - k].abs() < 2e-8);
        assert_eq!(mode.shape, raw.modes[1 - k].shape);
    }
}

#[test]
fn storage_and_fit_refusals_remain_typed_and_windows_are_admitted_first() {
    let oscillator = Oscillator::new(0.2, 25.0, 0.02, 0.5, 1.0);
    let dt = 1.0 / 512.0;
    let record: Vec<_> = (0..129).map(|i| oscillator.response(i as f64 * dt)).collect();
    let poles = [oscillator.pole(0.0)];
    let bounds = bounds(&poles);
    for tau in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE / 4.0] {
        assert!(matches!(
            identify_ringdown(&[], 0.0, &[], &[], budget(), Some(tau)),
            Err(RingdownError::Mapping(ModalIdError::Shape { .. }))
        ));
    }
    assert!(matches!(
        identify_ringdown(&[&record, &record[..128]], dt, &poles, &bounds, budget(), None),
        Err(RingdownError::Fit(DampedToneRefinementError::InvalidParameter { .. }))
    ));
    let silent = vec![0.0; record.len()];
    assert!(matches!(
        identify_ringdown(&[&silent], dt, &poles, &bounds, budget(), None),
        Err(RingdownError::Fit(DampedToneRefinementError::Unidentifiable { .. }))
    ));
    let mut small = budget();
    small.fit.max_workspace_bytes = 1;
    assert!(matches!(
        identify_ringdown(&[&record], dt, &poles, &bounds, small, None),
        Err(RingdownError::Fit(DampedToneRefinementError::Fit(HarmonicFitError::BudgetExceeded {
            required_bytes, max_bytes: 1,
        }))) if required_bytes > 1
    ));
    // The shared fit's own minimum cap is insufficient for the consumer's
    // extra table. Check the public refusal, without pinning a private formula.
    small.fit.max_workspace_bytes = 0;
    let required = match refine_damped_tone_poles_channels(&[&record], dt, &poles, &bounds, small) {
        Err(DampedToneRefinementError::Fit(HarmonicFitError::BudgetExceeded { required_bytes, .. })) => required_bytes,
        result => panic!("expected a numeric storage refusal: {result:?}"),
    };
    small.fit.max_workspace_bytes = required;
    assert!(matches!(
        identify_ringdown(&[&record], dt, &poles, &bounds, small, None),
        Err(RingdownError::Fit(DampedToneRefinementError::Fit(HarmonicFitError::BudgetExceeded {
            required_bytes, max_bytes,
        }))) if required_bytes == required && max_bytes < required
    ));
}
