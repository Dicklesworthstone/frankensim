//! Original recorded stress histories through passive-wall reflection.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::quadrupole::recorded::{
    QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
    RecordedQuadrupoleSpectrum,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    ImpedanceWallFieldKind, PlanarImpedanceWall, PreparedImpedanceWallQuadrupoleObserver,
};
use fs_aeroac::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const BIN: usize = 5;
const HZ: f64 = 640.0;
const POSITIONS: [[f64; 3]; 2] = [[0.0, 0.0, -0.21], [-0.07, 0.03, -0.14]];
const RECEIVER: [f64; 3] = [0.13, -0.05, -0.06];

fn limits() -> FullWaveQuadrupoleSourceLimits {
    FullWaveQuadrupoleSourceLimits {
        max_sources: 4,
        maximum_total_evaluations: 524_288,
        quadrupole: FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 1e-8,
            absolute_velocity_tolerance_m_s_per_n_m: 1e-11,
            relative_tolerance: 1e-10,
            maximum_evaluations: 262_144,
            maximum_depth: 40,
        },
    }
}
fn work() -> CoherenceBudget {
    CoherenceBudget { max_modes: 2, max_source_mode_pairs: 8 }
}
fn wall(normalized: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        normalized.map(|z| z.scale(1.2 * 343.0)),
    ).unwrap()
}
fn sources() -> [QuadrupolePhasor; 2] {
    core::array::from_fn(|cell| QuadrupolePhasor {
        stress_moment_n_m: core::array::from_fn(|entry| C64::new(
            1e-4 * (entry + 1) as f64 * if cell == 0 { 1.0 } else { -0.4 },
            2e-5 * (entry + cell + 1) as f64,
        )),
    })
}
fn record(
    positions: &[[f64; 3]],
    amplitudes: &[QuadrupolePhasor],
    window: WelchWindow,
) -> RecordedQuadrupoleSpectrum {
    let clock = QuadrupoleSampleClock::new(0.125, 1.0 / 8192.0, 0.0).unwrap();
    let config = WelchConfig {
        segment_samples: 64, hop_samples: 32, window, detrend: WelchDetrend::None,
    };
    let budget = QuadrupoleRecordingBudget {
        maximum_sources: positions.len(),
        welch: WelchBudget {
            maximum_channels: 6 * positions.len(), maximum_samples_per_channel: 96,
            maximum_transforms: 128, maximum_workspace_bytes: 1_000_000,
        },
    };
    let mut stream =
        QuadrupoleWelchStream::new(positions, clock, 96, config, budget, || false).unwrap();
    for index in 0..96 {
        let phase = core::f64::consts::TAU * BIN as f64 * index as f64 / 64.0;
        let (sin, cos) = phase.sin_cos();
        // Re(M exp(-i omega t)) preserves every relative component/cell phase.
        let moments: Vec<_> = amplitudes.iter().map(|source| {
            source.stress_moment_n_m.map(|m| m.re * cos + m.im * sin)
        }).collect();
        stream.push_sample(positions, clock.time_at(index).unwrap(), &moments, || false)
            .unwrap();
    }
    stream.finish().unwrap()
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 1e-21 + 3e-10 * expected.abs(),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn field_moments(field: &CoherentCompactField) -> CompactFieldStatistics {
    CompactFieldStatistics {
        mean_square_pressure_pa2: 0.5 * field.pressure_pa.abs().powi(2),
        velocity_covariance_m2_s2: core::array::from_fn(|i| core::array::from_fn(|j| {
            field.particle_velocity_m_s[i].scale(0.5)
                * field.particle_velocity_m_s[j].conj()
        })),
        active_intensity_w_m2: field.active_intensity_w_m2,
        reactive_intensity_w_m2: field.reactive_intensity_w_m2,
        scope: "Independent peak-phasor second moments",
    }
}
fn same_moments(actual: &CompactFieldStatistics, expected: &CompactFieldStatistics, scale: f64) {
    close(C64::from_re(actual.mean_square_pressure_pa2),
        C64::from_re(scale * expected.mean_square_pressure_pa2));
    for i in 0..3 {
        for j in 0..3 {
            close(actual.velocity_covariance_m2_s2[i][j],
                expected.velocity_covariance_m2_s2[i][j].scale(scale));
        }
        close(C64::new(actual.active_intensity_w_m2[i], actual.reactive_intensity_w_m2[i]),
            C64::new(expected.active_intensity_w_m2[i], expected.reactive_intensity_w_m2[i])
                .scale(scale));
    }
}

#[test]
fn original_tone_bins_keep_peak_scaling_full_moments_and_wall_absorption() {
    let sources = sources();
    let recording = record(&POSITIONS, &sources, WelchWindow::Rectangular);
    let boundary = wall(Some(C64::new(1.7, 0.4)));
    let receiver = [0.13, -0.05, 0.0];
    let observer = recording.prepare_impedance_wall_total_observer(
        boundary, receiver, BIN, limits(), || false,
    ).unwrap();
    assert_eq!(observer.frequency_hz(), HZ);
    assert_eq!(observer.field_kind(), ImpedanceWallFieldKind::Total);
    assert_eq!(recording.clock().epoch_s(), 0.125);
    assert_eq!(recording.clock().timestep_s(), 1.0 / 8192.0);
    let expected = field_moments(&observer.evaluate(&sources).unwrap());
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0), (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let recording = record(&POSITIONS, &sources, window);
        let actual = recording.evaluate_impedance_wall_bin(&observer, BIN, work(), || false)
            .unwrap();
        // A centered periodic-Hann tone has 2/3 of its power in the center bin.
        // This must remain bin power, without coherent-gain tone renormalization.
        same_moments(&actual, &expected, fraction);
        same_moments(&observer.evaluate_welch_bin(recording.spectrum(), BIN, work(), || false)
            .unwrap(), &actual, 1.0);
        let z = boundary.impedance_pa_s_m().unwrap();
        let normal_variance = actual.velocity_covariance_m2_s2[2][2].re;
        assert!(normal_variance > 0.0);
        let absorbed = z.re * normal_variance;
        assert!((actual.active_intensity_w_m2[2] - absorbed).abs() < 1e-8 * absorbed + 1e-18);
        assert!((actual.reactive_intensity_w_m2[2] - z.im * normal_variance).abs()
            < 1e-8 * absorbed + 1e-18);
        assert!((actual.mean_square_pressure_pa2 - z.abs().powi(2) * normal_variance).abs()
            < 1e-8 * actual.mean_square_pressure_pa2 + 1e-15);
    }
}

#[test]
fn ideal_wall_recordings_match_original_direct_and_tensor_image_histories() {
    let amplitudes = sources();
    let recording = record(&POSITIONS, &amplitudes, WelchWindow::PeriodicHann);
    let image_positions = POSITIONS.map(|p| [p[0], p[1], -p[2]]);
    for (impedance, parity) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        let boundary = wall(impedance);
        // Both stress indices reflect: S M S with S=diag(1,1,-1).
        let image_sources = amplitudes.map(|source| QuadrupolePhasor {
            stress_moment_n_m: core::array::from_fn(|entry| {
                source.stress_moment_n_m[entry].scale(if entry >= 4 { -parity } else { parity })
            }),
        });
        for total in [false, true] {
            let observer = if total {
                recording.prepare_impedance_wall_total_observer(
                    boundary, RECEIVER, BIN, limits(), || false)
            } else {
                recording.prepare_impedance_wall_reflected_observer(
                    boundary, RECEIVER, BIN, limits(), || false)
            }.unwrap();
            let mut positions = if total { POSITIONS.to_vec() } else { Vec::new() };
            let mut moments = if total { amplitudes.to_vec() } else { Vec::new() };
            positions.extend_from_slice(&image_positions);
            moments.extend_from_slice(&image_sources);
            let reference = record(&positions, &moments, WelchWindow::PeriodicHann);
            let free = reference.prepare_observer(boundary.fluid(), RECEIVER, BIN, || false)
                .unwrap();
            let expected = reference.evaluate_bin(&free, BIN, work(), || false).unwrap();
            let actual = recording.evaluate_impedance_wall_bin(&observer, BIN, work(), || false)
                .unwrap();
            assert_eq!(observer.evaluations(), 0);
            same_moments(&actual, &expected, 1.0);
        }
    }
}

#[test]
fn measured_cell_phase_and_transverse_motion_survive_pressure_nulls() {
    let positions = [[0.0, 0.0, -0.2]; 2];
    let mut shear = QuadrupolePhasor::default();
    shear.stress_moment_n_m[4] = C64::new(0.001, -0.0003);
    let opposite = QuadrupolePhasor {
        stress_moment_n_m: shear.stress_moment_n_m.map(|m| -m),
    };
    let recording = record(&positions, &[shear, opposite], WelchWindow::Rectangular);
    let observer = recording.prepare_impedance_wall_total_observer(
        wall(None), [0.0, 0.0, -0.07], BIN, limits(), || false,
    ).unwrap();
    let actual = recording.evaluate_impedance_wall_bin(&observer, BIN, work(), || false)
        .unwrap();
    let coherent = field_moments(&observer.evaluate(&[shear, opposite]).unwrap());
    same_moments(&actual, &coherent, 1.0);
    let independent = observer.evaluate_independent(&[shear, opposite], || false).unwrap();
    assert!(independent.velocity_covariance_m2_s2[0][0].re > 0.0);
    assert!(actual.velocity_covariance_m2_s2[0][0].re
        < 1e-20 * independent.velocity_covariance_m2_s2[0][0].re);
    let one = [shear, QuadrupolePhasor::default()];
    let recording = record(&positions, &one, WelchWindow::Rectangular);
    let actual = recording.evaluate_impedance_wall_bin(&observer, BIN, work(), || false)
        .unwrap();
    same_moments(&actual, &field_moments(&observer.evaluate(&one).unwrap()), 1.0);
    assert!(actual.mean_square_pressure_pa2 < 1e-28);
    assert!(actual.velocity_covariance_m2_s2[0][0].re > 0.0);
}

#[test]
fn original_layout_frequency_work_and_self_field_limits_bind_recorded_observation() {
    let recording = record(&POSITIONS, &sources(), WelchWindow::Rectangular);
    let boundary = wall(None);
    let observer = recording.prepare_impedance_wall_reflected_observer(
        boundary, RECEIVER, BIN, limits(), || false,
    ).unwrap();
    for positions in [
        vec![POSITIONS[1], POSITIONS[0]],
        vec![[0.01, 0.0, -0.21], POSITIONS[1]],
        vec![POSITIONS[0]],
    ] {
        let other = PreparedImpedanceWallQuadrupoleObserver::new_reflected(
            boundary, HZ, RECEIVER, &positions, limits(), &mut || false,
        ).unwrap();
        assert!(recording.evaluate_impedance_wall_bin(&other, BIN, work(), || false).is_err());
    }
    assert!(recording.evaluate_impedance_wall_bin(&observer, BIN + 1, work(), || false).is_err());
    for bin in [0, 32] {
        assert!(recording.prepare_impedance_wall_total_observer(
            boundary, RECEIVER, bin, limits(), || false,
        ).is_err());
    }
    for budget in [
        CoherenceBudget { max_modes: 1, max_source_mode_pairs: 4 },
        CoherenceBudget { max_modes: 2, max_source_mode_pairs: 3 },
    ] {
        assert!(recording.evaluate_impedance_wall_bin(&observer, BIN, budget, || false).is_err());
    }
    assert!(matches!(
        recording.prepare_impedance_wall_total_observer(
            boundary, RECEIVER, BIN, limits(), || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        recording.evaluate_impedance_wall_bin(&observer, BIN, work(), || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let returned = recording.prepare_impedance_wall_reflected_observer(
        boundary, POSITIONS[0], BIN, limits(), || false,
    ).unwrap();
    assert!(recording.evaluate_impedance_wall_bin(&returned, BIN, work(), || false)
        .unwrap().mean_square_pressure_pa2.is_finite());
    assert!(recording.prepare_impedance_wall_total_observer(
        boundary, POSITIONS[0], BIN, limits(), || false,
    ).is_err());
}
