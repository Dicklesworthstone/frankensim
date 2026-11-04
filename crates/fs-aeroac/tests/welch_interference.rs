//! Recorded-source interference: common-clock phase and physical bin power.
use fs_aeroac::curle2d::{
    CoherenceBudget, CoherentDipoleField, DipoleFieldStatistics, PreparedDipoleObserver,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const N: usize = 32;
const BIN: usize = 5;
const DT: f64 = 1.0 / 16_000.0;
const WORK: CoherenceBudget = CoherenceBudget {
    max_modes: 3,
    max_source_mode_pairs: 6,
};

fn config(window: WelchWindow) -> WelchConfig {
    WelchConfig {
        segment_samples: N,
        hop_samples: N,
        window,
        detrend: WelchDetrend::None,
    }
}

fn spectral_budget() -> WelchBudget {
    WelchBudget {
        maximum_channels: 4,
        maximum_samples_per_channel: N * 3,
        maximum_transforms: 12,
        maximum_workspace_bytes: 100_000,
    }
}

// Each segment is a known real force history, in the negative-time phasor
// convention. No estimator output is used to build the physical reference.
fn record(segments: &[[[C64; 2]; 2]]) -> Vec<[f64; 4]> {
    segments
        .iter()
        .flat_map(|sources| {
            (0..N).map(move |j| {
                let theta = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                core::array::from_fn(|channel| {
                    let a = sources[channel / 2][channel % 2];
                    a.re * theta.cos() + a.im * theta.sin()
                })
            })
        })
        .collect()
}

fn estimate(rows: &[[f64; 4]], window: WelchWindow) -> WelchCrossSpectrum {
    WelchCrossSpectrum::estimate_interleaved(rows, DT, config(window), spectral_budget(), || false)
        .unwrap()
}

fn observer(spectrum: &WelchCrossSpectrum, positions: &[[f64; 2]]) -> PreparedDipoleObserver {
    PreparedDipoleObserver::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        spectrum.frequency_hz(BIN).unwrap(),
        [0.38, 0.21],
        positions,
        2,
    )
    .unwrap()
}

fn close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 2.0e-11 * b.abs().max(1.0e-20),
        "actual {a:?}, expected {b:?}"
    );
}

// Average independently propagated, known physical force phasors. The taper
// weight is an analytic bin-energy fraction, not a measured PSD or factor.
fn check_moments(actual: &DipoleFieldStatistics, fields: &[CoherentDipoleField], weight: f64) {
    let moment = |left: fn(&CoherentDipoleField) -> C64, right: fn(&CoherentDipoleField) -> C64| {
        fields
            .iter()
            .fold(C64::ZERO, |sum, field| {
                sum + left(field) * right(field).conj()
            })
            .scale(0.5 * weight / fields.len() as f64)
    };
    close(
        C64::from_re(actual.mean_square_pressure_pa2),
        moment(|f| f.pressure_pa, |f| f.pressure_pa),
    );
    let velocity: [fn(&CoherentDipoleField) -> C64; 2] = [
        |f| f.particle_velocity_m_s[0],
        |f| f.particle_velocity_m_s[1],
    ];
    for axis in 0..2 {
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            moment(|f| f.pressure_pa, velocity[axis]),
        );
        for other in 0..2 {
            close(
                actual.velocity_covariance_m2_s2[axis][other],
                moment(velocity[axis], velocity[other]),
            );
        }
    }
}

#[test]
fn g1_recorded_cross_source_and_cross_axis_phase_reaches_all_field_moments() {
    let segments = [
        [
            [C64::new(1.0, 0.3), C64::new(-0.2, 0.7)],
            [C64::new(0.4, -0.8), C64::new(-0.3, 0.2)],
        ],
        [
            [C64::new(0.2, -0.5), C64::new(0.6, 0.1)],
            [C64::new(-0.7, 0.4), C64::new(0.1, -0.9)],
        ],
        [
            [C64::new(-0.4, 0.9), C64::new(0.2, 0.3)],
            [C64::new(0.7, 0.6), C64::new(-0.5, 0.4)],
        ],
    ];
    let rows = record(&segments);
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0),
        // A periodic Hann puts 2/3 of a resolved tone's energy in its central
        // bin. Coherent-gain correction here would invent 50% more power.
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = estimate(&rows, window);
        let prepared = observer(&spectrum, &[[-0.07, 0.025], [0.11, -0.04]]);
        let physical: Vec<_> = segments
            .iter()
            .map(|force| prepared.evaluate(force).unwrap())
            .collect();
        let actual = prepared
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .unwrap();
        check_moments(&actual, &physical, fraction);
        assert!(actual.velocity_covariance_m2_s2[0][1].im.abs() > 1.0e-8);
    }
}

#[test]
fn g1_measured_equal_opposite_and_partial_coherence_have_distinct_power() {
    let force = [C64::new(1.0, 0.3), C64::new(-0.2, 0.7)];
    let opposite = force.map(|a| -a);
    for (segments, multiplier) in [
        ([[force, force], [force, force]], 4.0),
        ([[force, opposite], [force, opposite]], 0.0),
        ([[force, force], [force, opposite]], 2.0),
    ] {
        let spectrum = estimate(&record(&segments), WelchWindow::Rectangular);
        let prepared = observer(&spectrum, &[[0.0; 2]; 2]);
        let single = prepared.evaluate(&[force, [C64::ZERO; 2]]).unwrap();
        let actual = prepared
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .unwrap();
        close(
            C64::from_re(actual.mean_square_pressure_pa2),
            C64::from_re(0.5 * single.pressure_pa.norm_sq() * multiplier),
        );
        if multiplier == 0.0 {
            assert_eq!(actual.mean_square_pressure_pa2, 0.0);
            assert_eq!(actual.velocity_covariance_m2_s2, [[C64::ZERO; 2]; 2]);
            assert_eq!(actual.active_intensity_w_m2, [0.0; 2]);
        }
    }
}

#[test]
fn g1_recorded_pressure_antinode_carries_no_net_energy_flux() {
    let force = [C64::ONE, C64::ZERO];
    let rows = record(&[[force, force.map(|a| -a)]]);
    let spectrum = estimate(&rows, WelchWindow::Rectangular);
    let prepared = PreparedDipoleObserver::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        spectrum.frequency_hz(BIN).unwrap(),
        [0.0; 2],
        &[[-0.3, 0.0], [0.3, 0.0]],
        2,
    )
    .unwrap();
    let actual = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    assert!(actual.mean_square_pressure_pa2 > 1.0);
    assert_eq!(actual.velocity_covariance_m2_s2, [[C64::ZERO; 2]; 2]);
    assert_eq!(actual.active_intensity_w_m2, [0.0; 2]);
}

#[test]
fn g0_recorded_source_admission_cancellation_and_replay_are_explicit() {
    let force = [C64::new(1.0, 0.2), C64::new(-0.1, 0.4)];
    let rows = record(&[[force, force]; 2]);
    let spectrum = estimate(&rows, WelchWindow::Rectangular);
    let prepared = observer(&spectrum, &[[0.0; 2]; 2]);
    for bin in [
        0,
        spectrum.bin_count() - 1,
        spectrum.bin_count(),
        usize::MAX,
    ] {
        assert!(
            prepared
                .evaluate_welch_bin(&spectrum, bin, WORK, || false)
                .is_err()
        );
    }
    let wrong_sources = observer(&spectrum, &[[0.0; 2]]);
    assert!(
        wrong_sources
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .is_err()
    );
    let wrong_frequency = PreparedDipoleObserver::new(
        prepared.acoustic(),
        1.01 * prepared.frequency_hz(),
        [0.38, 0.21],
        &[[0.0; 2]; 2],
        2,
    )
    .unwrap();
    assert!(
        wrong_frequency
            .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
            .is_err()
    );
    for budget in [
        CoherenceBudget {
            max_modes: 1,
            max_source_mode_pairs: 4,
        },
        CoherenceBudget {
            max_modes: 2,
            max_source_mode_pairs: 3,
        },
    ] {
        assert!(
            prepared
                .evaluate_welch_bin(&spectrum, BIN, budget, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(
        prepared
            .evaluate_welch_bin(&spectrum, BIN, WORK, || {
                polls += 1;
                polls == 4
            })
            .is_err()
    );
    assert_eq!(polls, 4);
    let first = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    let replay = prepared
        .evaluate_welch_bin(&spectrum, BIN, WORK, || false)
        .unwrap();
    assert_eq!(
        first.mean_square_pressure_pa2.to_bits(),
        replay.mean_square_pressure_pa2.to_bits()
    );
    assert_eq!(
        first.velocity_covariance_m2_s2,
        replay.velocity_covariance_m2_s2
    );
    assert_eq!(first.active_intensity_w_m2, replay.active_intensity_w_m2);
    assert_eq!(
        first.reactive_intensity_w_m2,
        replay.reactive_intensity_w_m2
    );
}
