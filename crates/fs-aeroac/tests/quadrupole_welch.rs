//! Recorded physical volume stresses retain source and tensor coherence.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{
    PreparedQuadrupoleObserver, QuadrupolePhasor, resolved_reynolds_moment_n_m,
};
use fs_aeroac::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const N: usize = 64;
const BIN: usize = 7;
const DT: f64 = 1.0 / 8192.0;
const POSITIONS: [[f64; 3]; 2] = [[0.02, -0.01, 0.03], [0.07, 0.04, -0.02]];
const RECEIVER: [f64; 3] = [0.11, 0.09, 0.08];
const WORK: CoherenceBudget = CoherenceBudget {
    max_modes: 3,
    max_source_mode_pairs: 6,
};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.19, 344.0).unwrap()
}

fn estimate_record<const D: usize>(record: &[[f64; D]], window: WelchWindow) -> WelchCrossSpectrum {
    WelchCrossSpectrum::estimate_interleaved(
        record,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 12,
            maximum_samples_per_channel: 3 * N,
            maximum_transforms: 36,
            maximum_workspace_bytes: 1024 * 1024,
        },
        || false,
    )
    .unwrap()
}

fn estimate<const D: usize>(segments: &[[C64; D]], window: WelchWindow) -> WelchCrossSpectrum {
    let record: Vec<_> = segments
        .iter()
        .flat_map(|amplitudes| {
            (0..N).map(move |j| {
                let theta = TAU * BIN as f64 * j as f64 / N as f64;
                // Actual real samples on the common negative-time source clock.
                amplitudes.map(|a| a.re * theta.cos() + a.im * theta.sin())
            })
        })
        .collect();
    estimate_record(&record, window)
}

fn sources() -> [[QuadrupolePhasor; 2]; 3] {
    let a = [
        C64::new(0.002, -0.001),
        C64::new(-0.001, 0.0007),
        C64::new(0.0008, 0.0003),
        C64::new(0.0015, -0.0004),
        C64::new(-0.0006, 0.0009),
        C64::new(0.0005, 0.0002),
    ];
    let b = [
        C64::new(-0.001, 0.0015),
        C64::new(0.0009, -0.0003),
        C64::new(-0.0004, -0.0006),
        C64::new(0.0007, 0.0011),
        C64::new(0.0002, -0.0008),
        C64::new(-0.0012, 0.0004),
    ];
    [C64::ONE, C64::new(0.0, 1.0), C64::new(-0.5, 0.2)].map(|phase| {
        [
            QuadrupolePhasor {
                stress_moment_n_m: a,
            },
            QuadrupolePhasor {
                stress_moment_n_m: b.map(|q| q * phase),
            },
        ]
    })
}

fn channels(sources: [QuadrupolePhasor; 2]) -> [C64; 12] {
    core::array::from_fn(|i| sources[i / 6].stress_moment_n_m[i % 6])
}

fn observer(frequency_hz: f64) -> PreparedQuadrupoleObserver {
    PreparedQuadrupoleObserver::new(medium(), frequency_hz, RECEIVER, &POSITIONS, 2).unwrap()
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-18),
        "actual={actual:?}, expected={expected:?}"
    );
}

// Independent source-side reference: propagate known physical tensor phasors,
// then average their p/u products. No Welch coefficient enters this reference.
fn check_moments(actual: &CompactFieldStatistics, fields: &[CoherentCompactField], weight: f64) {
    let moment = |i: usize, j: usize| {
        let component = |field: &CoherentCompactField, axis| {
            if axis == 0 {
                field.pressure_pa
            } else {
                field.particle_velocity_m_s[axis - 1]
            }
        };
        fields
            .iter()
            .fold(C64::ZERO, |sum, field| {
                sum + component(field, i) * component(field, j).conj()
            })
            .scale(0.5 * weight / fields.len() as f64)
    };
    close(C64::from_re(actual.mean_square_pressure_pa2), moment(0, 0));
    for axis in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            moment(0, axis + 1),
        );
        for other in 0..3 {
            close(
                actual.velocity_covariance_m2_s2[axis][other],
                moment(axis + 1, other + 1),
            );
        }
    }
    assert_eq!(actual.scope, fields[0].scope);
}

#[test]
fn g1_recorded_tensor_and_cell_cross_terms_keep_physical_bin_power_and_near_fields() {
    let inputs = sources();
    for (window, bins) in [
        (WelchWindow::Rectangular, vec![(BIN, 1.0)]),
        (
            WelchWindow::PeriodicHann,
            vec![(BIN - 1, 1.0 / 6.0), (BIN, 2.0 / 3.0), (BIN + 1, 1.0 / 6.0)],
        ),
    ] {
        let spectrum = estimate(&inputs.map(channels), window);
        // Periodic Hann spreads a bin-centered sinusoid's power 1/6, 2/3, 1/6.
        // Each recorded bin still propagates at its own physical frequency.
        for (bin, fraction) in bins {
            let prepared = observer(spectrum.frequency_hz(bin).unwrap());
            let fields: Vec<_> = inputs
                .iter()
                .map(|s| prepared.evaluate(s).unwrap())
                .collect();
            let actual = prepared
                .evaluate_welch_bin(&spectrum, bin, WORK, || false)
                .unwrap();
            check_moments(&actual, &fields, fraction);
            assert!(actual.velocity_covariance_m2_s2[0][1].im.abs() > 1e-12);
            assert!(
                actual
                    .reactive_intensity_w_m2
                    .iter()
                    .any(|q| q.abs() > 1e-9)
            );
            assert_eq!(prepared.acoustic(), medium());
        }
    }
}

#[test]
fn g1_original_velocity_history_produces_dc_and_twice_frequency_stress_radiation() {
    let velocity_bin = 3;
    let [a, b, c] = [3.0, 2.0, -1.0];
    let rho = 1.3;
    let volume = 0.002;
    let mass = rho * volume;
    let record: Vec<_> = (0..N)
        .map(|j| {
            let theta = TAU * velocity_bin as f64 * j as f64 / N as f64;
            resolved_reynolds_moment_n_m(
                rho,
                volume,
                [a * theta.cos(), b * theta.sin(), c * theta.cos()],
            )
            .unwrap()
        })
        .collect();
    let spectrum = estimate_record(&record, WelchWindow::Rectangular);
    let bin = 2 * velocity_bin;
    let prepared = PreparedQuadrupoleObserver::new(
        medium(),
        spectrum.frequency_hz(bin).unwrap(),
        RECEIVER,
        &[[0.0; 3]],
        1,
    )
    .unwrap();
    // cos²=(1+cos(2t))/2, sin²=(1-cos(2t))/2, sin*cos=sin(2t)/2.
    // Shear entries are physical Qxy etc.; the transfer applies symmetry once.
    let expected_source = QuadrupolePhasor {
        stress_moment_n_m: [
            C64::from_re(0.5 * mass * a * a),
            C64::from_re(-0.5 * mass * b * b),
            C64::from_re(0.5 * mass * c * c),
            C64::new(0.0, 0.5 * mass * a * b),
            C64::from_re(0.5 * mass * a * c),
            C64::new(0.0, 0.5 * mass * b * c),
        ],
    };
    let expected = prepared.evaluate(&[expected_source]).unwrap();
    let actual = prepared
        .evaluate_welch_bin(&spectrum, bin, WORK, || false)
        .unwrap();
    check_moments(&actual, &[expected], 1.0);
    assert!(spectrum.cross_spectrum(0, 0, 0).unwrap().re > 0.0);
    assert!(matches!(
        prepared.evaluate_welch_bin(&spectrum, 0, WORK, || false),
        Err(AeroacError::Welch(_))
    ));
    let fundamental = PreparedQuadrupoleObserver::new(
        medium(),
        spectrum.frequency_hz(velocity_bin).unwrap(),
        RECEIVER,
        &[[0.0; 3]],
        1,
    )
    .unwrap()
    .evaluate_welch_bin(&spectrum, velocity_bin, WORK, || false)
    .unwrap();
    assert!(fundamental.mean_square_pressure_pa2 < actual.mean_square_pressure_pa2 * 1e-24);
}

#[test]
fn g3_recorded_cell_interference_and_pressure_null_transverse_velocity_survive() {
    let frequency = BIN as f64 / (N as f64 * DT);
    let a = QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.003, 0.001),
            C64::ZERO,
            C64::ZERO,
            C64::new(0.001, -0.002),
            C64::ZERO,
            C64::ZERO,
        ],
    };
    let minus_a = QuadrupolePhasor {
        stress_moment_n_m: a.stress_moment_n_m.map(|q| -q),
    };
    let receiver = [0.21, 0.0, 0.0];
    let prepared =
        PreparedQuadrupoleObserver::new(medium(), frequency, receiver, &[[0.0; 3]; 2], 2).unwrap();
    let one =
        PreparedQuadrupoleObserver::new(medium(), frequency, receiver, &[[0.0; 3]], 1).unwrap();
    let single = one.evaluate(&[a]).unwrap();
    let single_power = 0.5 * single.pressure_pa.norm_sq();
    let destructive = estimate(&[channels([a, minus_a]); 2], WelchWindow::Rectangular);
    let quiet = prepared
        .evaluate_welch_bin(&destructive, BIN, WORK, || false)
        .unwrap();
    assert!(quiet.mean_square_pressure_pa2 < single_power * 1e-24);
    let decorrelated = estimate(
        &[channels([a, a]), channels([a, minus_a])],
        WelchWindow::Rectangular,
    );
    let restored = prepared
        .evaluate_welch_bin(&decorrelated, BIN, WORK, || false)
        .unwrap();
    close(
        C64::from_re(restored.mean_square_pressure_pa2),
        C64::from_re(2.0 * single_power),
    );
    let constructive = estimate(&[channels([a, a]); 2], WelchWindow::Rectangular);
    let loud = prepared
        .evaluate_welch_bin(&constructive, BIN, WORK, || false)
        .unwrap();
    close(
        C64::from_re(loud.mean_square_pressure_pa2),
        C64::from_re(4.0 * single_power),
    );

    // Qxy has exactly zero pressure on the x axis, but nonzero transverse u_y.
    let shear = [
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ONE,
        C64::ZERO,
        C64::ZERO,
    ];
    let transverse = estimate(&[shear], WelchWindow::Rectangular);
    let field = one
        .evaluate_welch_bin(&transverse, BIN, WORK, || false)
        .unwrap();
    assert_eq!(field.mean_square_pressure_pa2, 0.0);
    assert!(field.velocity_covariance_m2_s2[1][1].re > 0.0);
    assert_eq!(field.velocity_covariance_m2_s2[0][0], C64::ZERO);
    assert_eq!(field.velocity_covariance_m2_s2[2][2], C64::ZERO);
}

#[test]
fn g4_recorded_shape_frequency_work_and_cancellation_refuse_without_changing_inputs() {
    let inputs = sources();
    let spectrum = estimate(&inputs.map(channels), WelchWindow::Rectangular);
    let frequency = spectrum.frequency_hz(BIN).unwrap();
    let prepared = observer(frequency);
    let original_cross = spectrum.cross_spectrum(BIN, 4, 10).unwrap();
    let original_field = prepared.evaluate(&inputs[0]).unwrap();
    for bin in [0, N / 2, BIN + 1, usize::MAX] {
        assert!(matches!(
            prepared.evaluate_welch_bin(&spectrum, bin, WORK, || false),
            Err(AeroacError::Welch(_))
        ));
    }
    assert!(matches!(
        observer(f64::from_bits(frequency.to_bits() + 1)).evaluate_welch_bin(
            &spectrum,
            BIN,
            WORK,
            || false
        ),
        Err(AeroacError::Welch(_))
    ));
    let wrong_shape = estimate(&[[C64::ONE; 11]; 3], WelchWindow::Rectangular);
    assert!(matches!(
        prepared.evaluate_welch_bin(&wrong_shape, BIN, WORK, || false),
        Err(AeroacError::InvalidParameter { .. })
    ));
    for budget in [
        CoherenceBudget {
            max_modes: 2,
            ..WORK
        },
        CoherenceBudget {
            max_source_mode_pairs: 5,
            ..WORK
        },
    ] {
        assert!(matches!(
            prepared.evaluate_welch_bin(&spectrum, BIN, budget, || false),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    for refusal_poll in [1, 5] {
        let mut polls = 0;
        assert!(matches!(
            prepared.evaluate_welch_bin(&spectrum, BIN, WORK, || {
                polls += 1;
                polls == refusal_poll
            }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let result = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    let fields: Vec<_> = inputs
        .iter()
        .map(|s| prepared.evaluate(s).unwrap())
        .collect();
    check_moments(&result, &fields, 1.0);
    assert_eq!(spectrum.cross_spectrum(BIN, 4, 10).unwrap(), original_cross);
    let repeated = prepared.evaluate(&inputs[0]).unwrap();
    assert_eq!(repeated.pressure_pa, original_field.pressure_pa);
    assert_eq!(
        repeated.particle_velocity_m_s,
        original_field.particle_velocity_m_s
    );
}
