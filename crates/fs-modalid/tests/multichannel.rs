use fs_math::{c64::C64, det};
use fs_modalid::{
    FrfData, IdentifyOptions, Mode, ModeCompareOptions, StabilizationRule, compare_mode_tables,
    identify, identify_multichannel, mac,
};
use fs_vfit::mimo::common_pole_vector_fit;
use fs_vfit::{FitOptions, PoleTerm};

const TAU: f64 = core::f64::consts::TAU;

fn options(order: usize, tau: Option<f64>) -> IdentifyOptions {
    IdentifyOptions {
        min_order: order,
        max_order: order,
        order_step: 2,
        snr_floor: 10.0,
        rule: StabilizationRule {
            freq_tol: 0.02,
            damp_tol: 0.05,
            required_runs: 1,
        },
        window_tau: tau,
    }
}

fn receptance(w: f64, f: f64, zeta: f64, rate: f64) -> C64 {
    let s = C64::new(rate, w);
    let wn = TAU * f;
    (s * s + s.scale(2.0 * zeta * wn) + C64::from_re(wn * wn)).recip()
}

fn data(modes: &[(f64, f64)], shapes: &[Vec<f64>], rate: f64, scale: f64) -> FrfData {
    let omega: Vec<_> = (0..401)
        .map(|i| TAU * (20.0 + 240.0 * i as f64 / 400.0))
        .collect();
    let channels = shapes
        .iter()
        .map(|shape| {
            let values: Vec<_> = omega
                .iter()
                .map(|&w| {
                    modes
                        .iter()
                        .zip(shape)
                        .fold(C64::ZERO, |sum, (&(f, z), &weight)| {
                            sum + receptance(w, f, z, rate).scale(scale * weight)
                        })
                })
                .collect();
            let coherence = if shape.iter().all(|&x| x == 0.0) {
                None
            } else {
                Some(vec![1.0; omega.len()])
            };
            (values, coherence)
        })
        .collect();
    FrfData::new_laplace(omega, channels).unwrap()
}

fn truth_modes(modes: &[(f64, f64)], shapes: &[Vec<f64>]) -> Vec<Mode> {
    modes
        .iter()
        .enumerate()
        .map(|(k, &(f, z))| {
            let wd = TAU * f * det::sqrt((1.0 - z) * (1.0 + z));
            Mode {
                frequency_hz: f,
                damping_ratio: z,
                shape: shapes
                    .iter()
                    .map(|row| C64::new(0.0, -row[k] / (2.0 * wd)))
                    .collect(),
                frequency_ci_hz: f64::NAN,
                damping_ci: f64::NAN,
            }
        })
        .collect()
}

#[test]
fn joint_identification_recovers_modes_hidden_from_the_reference_sensor() {
    let modes = [(100.0, 0.025), (190.0, 0.04)];
    // Channel 0 sits at both modal nodes; neither nonzero channel sees both modes.
    let shapes = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.7]];
    let input = data(&modes, &shapes, 0.0, 1.0);
    let opts = options(4, None);
    let legacy = identify(&input, &opts);
    assert!(
        legacy
            .as_ref()
            .map_or(true, |result| result.modes.len() != 2)
    );
    let result = identify_multichannel(&input, &opts).unwrap();
    assert_eq!(result.modes.len(), 2);
    let truth = truth_modes(&modes, &shapes);
    for (actual, expected) in result.modes.iter().zip(&truth) {
        assert!((actual.frequency_hz - expected.frequency_hz).abs() < 1e-5);
        assert!((actual.damping_ratio - expected.damping_ratio).abs() < 1e-7);
        assert!(mac(&actual.shape, &expected.shape) > 0.999999);
        assert!(actual.frequency_ci_hz.is_finite() && actual.damping_ci.is_finite());
    }
    // The actual calibration consumer can now pair both observable modes.
    let comparison = compare_mode_tables(
        &truth,
        &result.modes,
        &ModeCompareOptions {
            mac_floor: 0.9999,
            max_cents: 0.001,
        },
    )
    .unwrap();
    assert_eq!(comparison.pairs.len(), 2);
    assert!(comparison.unmatched_reference.is_empty());
    assert!(comparison.unmatched_candidate.is_empty());
    assert!(comparison.pairs.iter().all(|pair| pair.within_tolerance));
}

#[test]
fn shared_denominator_does_not_cancel_modes_when_sensor_responses_cancel() {
    let modes = [(100.0, 0.025), (190.0, 0.04)];
    let shapes = vec![vec![1.0, 1.0], vec![1.0, -1.0], vec![0.0, 0.0]];
    let input = data(&modes, &shapes, 0.0, 1.0);
    let channels: Vec<_> = input
        .channels
        .iter()
        .map(|channel| channel.h.as_slice())
        .collect();
    let opts = FitOptions {
        fit_e: false,
        iterations: 60,
        ..FitOptions::new(4)
    };
    let (models, report) = common_pole_vector_fit(&input.omega, &channels, &opts).unwrap();
    assert_eq!(models.len(), 3);
    assert!(report.weighted_rms < 1e-8, "{report:?}");
    assert_eq!(models[0].poles_expanded(), models[1].poles_expanded());
    assert_eq!(models[0].poles_expanded(), models[2].poles_expanded());
    let mut frequencies: Vec<_> = models[0]
        .terms
        .iter()
        .filter_map(|term| match term {
            PoleTerm::Pair { pole, .. } => Some(pole.abs() / TAU),
            _ => None,
        })
        .collect();
    frequencies.sort_by(f64::total_cmp);
    assert_eq!(frequencies.len(), 2);
    for (got, &(f, _)) in frequencies.iter().zip(&modes) {
        assert!((got - f).abs() < 1e-5);
    }
    for (model, channel) in models.iter().zip(&channels) {
        for (&w, &want) in input.omega.iter().zip(channel.iter()) {
            assert!((model.eval_iw(w) - want).abs() < 1e-11);
        }
    }
    // Averaging these channels would erase the second resonance exactly.
    for (i, &w) in input.omega.iter().enumerate() {
        let mean = (channels[0][i] + channels[1][i]).scale(0.5);
        assert!((mean - receptance(w, 100.0, 0.025, 0.0)).abs() < 1e-18);
    }
}

#[test]
fn multichannel_modes_and_shapes_follow_sensor_permutation_and_unit_scaling() {
    let modes = [(100.0, 0.025), (190.0, 0.04)];
    let shapes = vec![vec![0.0, 0.0], vec![1.0, 0.2], vec![-0.3, 1.7]];
    let reference =
        identify_multichannel(&data(&modes, &shapes, 0.0, 1.0), &options(4, None)).unwrap();
    let permuted = vec![shapes[2].clone(), shapes[0].clone(), shapes[1].clone()];
    for scale in [1e-9, 1e9] {
        let result =
            identify_multichannel(&data(&modes, &permuted, 0.0, scale), &options(4, None)).unwrap();
        assert_eq!(result.modes.len(), 2);
        for (actual, old) in result.modes.iter().zip(&reference.modes) {
            assert!((actual.frequency_hz - old.frequency_hz).abs() < 1e-5);
            assert!((actual.damping_ratio - old.damping_ratio).abs() < 1e-7);
            let expected = [old.shape[2], old.shape[0], old.shape[1]];
            assert!(mac(&actual.shape, &expected) > 0.999999);
        }
    }
}

#[test]
fn multichannel_half_grid_spread_uses_observable_channels_and_unwindows_each_pole() {
    let rate = 200.0;
    let omega: Vec<_> = (0..802)
        .map(|i| TAU * (20.0 + 240.0 * i as f64 / 801.0))
        .collect();
    let mut channels = vec![vec![C64::ZERO; omega.len()]; 3];
    for (i, &w) in omega.iter().enumerate() {
        let sign = if i % 2 == 0 { -1.0 } else { 1.0 };
        channels[1][i] = receptance(w, 100.0 + sign * 0.15, 0.025 + sign * 0.001, rate);
        channels[2][i] = receptance(w, 190.0 + sign * 0.25, 0.04 + sign * 0.002, rate);
    }
    let input = FrfData::new_laplace(
        omega.clone(),
        channels
            .into_iter()
            .map(|h| (h, Some(vec![1.0; omega.len()])))
            .collect(),
    )
    .unwrap();
    let result = identify_multichannel(&input, &options(4, Some(1.0 / rate))).unwrap();
    assert_eq!(result.modes.len(), 2);
    for (mode, (f_spread, z_spread)) in result.modes.iter().zip([(0.15, 0.001), (0.25, 0.002)]) {
        assert!((mode.frequency_ci_hz - f_spread).abs() < 1e-5, "{mode:?}");
        assert!((mode.damping_ci - z_spread).abs() < 1e-7, "{mode:?}");
    }
}

#[test]
fn coincident_poles_report_a_combined_shape_without_inventing_multiplicity() {
    let modes = [(120.0, 0.03), (120.0, 0.03)];
    let shapes = vec![vec![1.0, 0.0], vec![0.0, 2.0], vec![1.0, 1.0]];
    let input = data(&modes, &shapes, 0.0, 1.0);
    let result = identify_multichannel(&input, &options(2, None)).unwrap();
    assert_eq!(result.modes.len(), 1);
    assert!((result.modes[0].frequency_hz - 120.0).abs() < 1e-5);
    let combined = [
        C64::new(0.0, -1.0),
        C64::new(0.0, -2.0),
        C64::new(0.0, -2.0),
    ];
    assert!(mac(&result.modes[0].shape, &combined) > 0.999999);
}

#[test]
fn opposite_polarity_accelerance_keeps_signed_feedthrough_out_of_mode_shapes() {
    let modes = [(100.0, 0.025), (190.0, 0.04)];
    let weights = [[1.0, 0.0], [0.0, -1.0], [0.0, 0.0]];
    let omega: Vec<_> = (0..401)
        .map(|i| TAU * (20.0 + 240.0 * i as f64 / 400.0))
        .collect();
    let channels = weights
        .iter()
        .map(|row| {
            let h = omega
                .iter()
                .map(|&w| {
                    modes
                        .iter()
                        .zip(row)
                        .fold(C64::ZERO, |sum, (&(f, z), &weight)| {
                            sum + receptance(w, f, z, 0.0).scale(-w * w * weight)
                        })
                })
                .collect();
            (h, Some(vec![1.0; omega.len()]))
        })
        .collect();
    let input = FrfData::new_laplace(omega, channels).unwrap();
    let result = identify_multichannel(&input, &options(4, None)).unwrap();
    assert_eq!(result.modes.len(), 2);
    for (k, mode) in result.modes.iter().enumerate() {
        let (f, z) = modes[k];
        let wn = TAU * f;
        let wd = wn * det::sqrt((1.0 - z) * (1.0 + z));
        let pole = C64::new(-z * wn, wd);
        for (channel, &actual) in mode.shape.iter().enumerate() {
            let want = pole * pole * C64::new(0.0, -weights[channel][k] / (2.0 * wd));
            assert!((actual - want).abs() < 1e-6 * wn, "{actual:?} vs {want:?}");
        }
    }
}
