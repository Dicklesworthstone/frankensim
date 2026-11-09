use super::*;
use crate::spectral::doppler::VACUUM_LIGHT_SPEED_M_S;
use crate::thermal_surface::{
    GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
};
use fs_blake3::ContentHash;

const C: f64 = VACUUM_LIGHT_SPEED_M_S;

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}

fn gray(temperature: f64) -> GrayDiffuseThermalSurface {
    GrayDiffuseThermalSurface::try_new(temperature, 0.7, ContentHash([11; 32])).unwrap()
}

fn motion(beta: f64, history: [f64; 2]) -> UniformVacuumWorldline {
    UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [0.0; 3],
        [beta * C, 0.0, 0.0],
        history.map(dd),
    )
    .unwrap()
}

fn receiver(position: [f64; 3]) -> VacuumReception {
    VacuumReception::try_new(Dd::ZERO, position, [0.0; 3]).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs(),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn g1_retarded_proper_history_drives_original_source_and_receiver_photon_rates() {
    let mut sampled = None;
    let mut source_identity = None;
    let event = RetardedThermalSurface::observe(
        &motion(0.6, [-10.0, 10.0]),
        receiver([C, 0.0, 0.0]),
        1e-6,
        |proper| {
            sampled = Some(proper);
            let surface = gray(4000.0 + 500.0 * proper.to_f64());
            source_identity = Some(surface.identity());
            Ok(ThermalEmissionState {
                // Current source state is 4000 K; the received photon was
                // emitted at proper time -2 s, when the material was 3000 K.
                surface: surface.into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
    )
    .unwrap();
    close(sampled.unwrap().to_f64(), -2.0);
    assert_eq!(sampled, Some(event.ray().emission_proper_time_s()));
    close(event.observed().source().temperature_k(), 3000.0);
    assert_eq!(Some(event.observed().source().identity()), source_identity);
    let sensor = QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (700.0, [1.0])]).unwrap();
    let direct = event
        .observed()
        .with_sensor(&sensor)
        .unwrap()
        .sample(0.5)
        .unwrap();
    let detector = event.with_sensor(&sensor).unwrap();
    let sample = detector.sample(0.5).unwrap();
    assert_eq!(sample, direct);
    assert_eq!(detector.ray(), event.ray());
    assert_eq!(detector.emission_cosine(), event.emission_cosine());
    close(sample.source_wavelength_nm, 1200.0);
    close(sample.receiver_wavelength_nm, 600.0);
    let lambda = 600e-9_f64;
    let h = 6.626_070_15e-34;
    let k = 1.380_649e-23;
    let expected_radiance =
        0.7 * 2.0 * h * C * C / (lambda.powi(5) * (h * C / (lambda * k * 6000.0)).exp_m1()) * 1e-9;
    close(sample.spectral_radiance_w_m2_sr_nm, expected_radiance);
    close(
        sample.electron_rate_per_m2_sr[0],
        expected_radiance * 200.0 * lambda / (h * C),
    );
}

#[test]
fn g1_emitting_side_uses_aberrated_material_direction_without_cosine_double_count() {
    let observe = |normal| {
        RetardedThermalSurface::observe(
            &motion(0.6, [-10.0, 10.0]),
            receiver([0.0, C, 0.0]),
            1e-6,
            |_| {
                Ok(ThermalEmissionState {
                    surface: gray(3000.0).into(),
                    outward_normal: normal,
                })
            },
        )
    };
    // This normal faces the WORLD ray (.6,.8,0) but faces away from the actual
    // SOURCE-frame ray (0,1,0). The emitting material determines its hemisphere.
    assert_eq!(
        observe([1.0, -0.1, 0.0]),
        Err(RetardedThermalError::NonEmittingSide)
    );
    let tilted = observe([-1.0, 0.1, 0.0]).unwrap();
    let normal = observe([0.0, 1.0, 0.0]).unwrap();
    close(tilted.emission_cosine(), 0.1 / 1.01_f64.sqrt());
    assert_eq!(
        tilted.observed().radiance_w_m2_sr_nm(600.0).unwrap(),
        normal.observed().radiance_w_m2_sr_nm(600.0).unwrap()
    );
    assert_eq!(observe([0.0; 3]), Err(RetardedThermalError::InvalidNormal));
}

#[test]
fn g3_unavailable_geometry_never_reads_history_and_original_source_errors_survive() {
    let mut called = false;
    let result = RetardedThermalSurface::observe(
        &motion(0.0, [0.0, 1.0]),
        receiver([C, 0.0, 0.0]),
        1e-6,
        |_| {
            called = true;
            Ok(ThermalEmissionState {
                surface: gray(3000.0).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
    );
    assert_eq!(
        result,
        Err(RetardedThermalError::Geometry(
            RetardedVacuumError::OutsideHistory
        ))
    );
    assert!(!called);
    let source_error =
        GrayDiffuseThermalSurface::try_new(-1.0, 0.7, ContentHash([11; 32])).unwrap_err();
    let result = RetardedThermalSurface::observe(
        &motion(0.0, [-10.0, 10.0]),
        receiver([C, 0.0, 0.0]),
        1e-6,
        |_| Err(source_error),
    );
    assert_eq!(result, Err(RetardedThermalError::Source(source_error)));
}

#[test]
fn g3_retarded_spectral_history_keeps_original_identity_and_full_sensor_band_admission() {
    let source = SpectralDiffuseThermalSurface::try_new(
        3000.0,
        SpectralEmissivity::constant_in_band(0.5, (800.0, 1600.0)).unwrap(),
        ContentHash([12; 32]),
    )
    .unwrap();
    let event = RetardedThermalSurface::observe(
        &motion(0.6, [-10.0, 10.0]),
        receiver([C, 0.0, 0.0]),
        1e-6,
        |_| {
            Ok(ThermalEmissionState {
                surface: (&source).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
    )
    .unwrap();
    assert_eq!(event.observed().source().identity(), source.identity());
    let valid = QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (700.0, [1.0])]).unwrap();
    assert!(event.with_sensor(&valid).unwrap().sample(0.5).is_ok());
    let unsupported =
        QuantumEfficiencySensor::try_new(vec![(700.0, [1.0]), (900.0, [1.0])]).unwrap();
    assert!(matches!(
        event.with_sensor(&unsupported),
        Err(ThermalDopplerError::SourceBand)
    ));
}

#[test]
fn g1_piecewise_motion_samples_accumulated_proper_history_and_original_qe() {
    use crate::spectral::doppler::retarded::VacuumTrajectoryKnot;

    let knots = [(0.0, 0.0), (2.0, 1.2), (4.0, 1.2)].map(|(time, x)| VacuumTrajectoryKnot {
        coordinate_time_s: dd(time),
        position_m: [x * C, 0.0, 0.0],
    });
    let motion = PiecewiseVacuumWorldline::try_new(&knots, dd(10.0), 2).unwrap();
    let observation = VacuumReception::try_new(dd(4.8), [3.0 * C, 0.0, 0.0], [0.0; 3]).unwrap();
    let mut calls = 0;
    let event = RetardedThermalSurface::observe_piecewise(&motion, observation, 1e-6, |proper| {
        calls += 1;
        close(proper.to_f64(), 12.6);
        Ok(ThermalEmissionState {
            surface: gray(1000.0 + 100.0 * proper.to_f64()).into(),
            outward_normal: [1.0, 0.0, 0.0],
        })
    })
    .unwrap();
    assert_eq!(calls, 1);
    close(event.observed().source().temperature_k(), 2260.0);
    assert_eq!(event.ray().proper_time_rate(), 1.0);
    let sensor = QuantumEfficiencySensor::try_new(vec![(500.0, [0.5]), (700.0, [0.5])]).unwrap();
    let sample = event.with_sensor(&sensor).unwrap().sample(0.5).unwrap();
    let original = gray(2260.0).radiance_w_m2_sr_nm(600.0).unwrap();
    close(sample.spectral_radiance_w_m2_sr_nm, original);
    close(
        sample.electron_rate_per_m2_sr[0],
        original * 100.0 * 600e-9 / (6.626_070_15e-34 * C),
    );

    let expired = VacuumReception::try_new(dd(7.0), [3.0 * C, 0.0, 0.0], [0.0; 3]).unwrap();
    let result = RetardedThermalSurface::observe_piecewise(&motion, expired, 1e-6, |_| {
        calls += 1;
        Ok(ThermalEmissionState {
            surface: gray(2260.0).into(),
            outward_normal: [1.0, 0.0, 0.0],
        })
    });
    assert_eq!(
        result,
        Err(RetardedThermalError::Geometry(
            RetardedVacuumError::OutsideHistory
        ))
    );
    assert_eq!(calls, 1);
}

#[test]
fn g1_moving_detector_proper_clock_drives_retarded_thermal_photon_rate() {
    let detector_motion = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        dd(10.0),
        [3.0 * C, 0.0, 0.0],
        [-0.6 * C, 0.0, 0.0],
        [dd(10.0), dd(11.0)],
    )
    .unwrap();
    let proper = dd(10.0) + dd(0.8);
    let mut calls = 0;
    let event = RetardedThermalSurface::observe_at_receiver_proper_time(
        &motion(0.0, [-10.0, 10.0]),
        &detector_motion,
        proper,
        1e-6,
        |emitted| {
            calls += 1;
            close(emitted.to_f64(), -1.4);
            Ok(ThermalEmissionState {
                surface: gray(4000.0 + 500.0 * emitted.to_f64()).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(event.ray().reception_proper_time_s(), Some(proper));
    close(event.ray().reception_coordinate_time_s().to_f64(), 1.0);
    close(event.observed().source().temperature_k(), 3300.0);
    let sensor = QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (700.0, [1.0])]).unwrap();
    let detector = event.with_sensor(&sensor).unwrap();
    assert_eq!(detector.ray().reception_proper_time_s(), Some(proper));
    let sample = detector.sample(0.5).unwrap();
    close(sample.source_wavelength_nm, 1200.0);
    let lambda = 600e-9_f64;
    let h = 6.626_070_15e-34;
    let k = 1.380_649e-23;
    let expected_radiance =
        0.7 * 2.0 * h * C * C / (lambda.powi(5) * (h * C / (lambda * k * 6600.0)).exp_m1()) * 1e-9;
    close(sample.spectral_radiance_w_m2_sr_nm, expected_radiance);
    // No additional coordinate/proper-clock factor belongs in this count rate.
    close(
        sample.electron_rate_per_m2_sr[0],
        expected_radiance * 200.0 * lambda / (h * C),
    );
}

#[test]
fn g3_piecewise_detector_clock_and_history_use_the_original_thermal_qe_route() {
    use crate::spectral::doppler::retarded::VacuumTrajectoryKnot;

    let knots = |values: &[(f64, f64)]| -> Vec<VacuumTrajectoryKnot> {
        values
            .iter()
            .map(|&(time, x)| VacuumTrajectoryKnot {
                coordinate_time_s: dd(time),
                position_m: [x * C, 0.0, 0.0],
            })
            .collect()
    };
    let source =
        PiecewiseVacuumWorldline::try_new(&knots(&[(-10.0, 0.0), (10.0, 0.0)]), dd(-10.0), 1)
            .unwrap();
    let receiver = PiecewiseVacuumWorldline::try_new(
        &knots(&[(0.0, 3.0), (2.0, 1.8), (4.0, 1.8)]),
        dd(10.0),
        2,
    )
    .unwrap();
    let mut calls = 0;
    let event = RetardedThermalSurface::observe_piecewise_at_receiver_proper_time(
        &source,
        &receiver,
        dd(12.6),
        1e-6,
        |emitted| {
            calls += 1;
            close(emitted.to_f64(), 1.2);
            Ok(ThermalEmissionState {
                surface: gray(4000.0 + 500.0 * emitted.to_f64()).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
    )
    .unwrap();
    assert_eq!(calls, 1);
    close(event.ray().reception_coordinate_time_s().to_f64(), 3.0);
    assert_eq!(event.ray().reception_proper_time_s(), Some(dd(12.6)));
    close(event.observed().source().temperature_k(), 4600.0);
    let sensor = QuantumEfficiencySensor::try_new(vec![(500.0, [0.5]), (700.0, [0.5])]).unwrap();
    let sample = event.with_sensor(&sensor).unwrap().sample(0.5).unwrap();
    let original = gray(4600.0).radiance_w_m2_sr_nm(600.0).unwrap();
    close(sample.spectral_radiance_w_m2_sr_nm, original);
    close(
        sample.electron_rate_per_m2_sr[0],
        original * 100.0 * 600e-9 / (6.626_070_15e-34 * C),
    );

    let result = RetardedThermalSurface::observe_piecewise_at_receiver_proper_time(
        &source,
        &receiver,
        dd(14.0),
        1e-6,
        |_| {
            calls += 1;
            Ok(ThermalEmissionState {
                surface: gray(4600.0).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
    );
    assert_eq!(
        result,
        Err(RetardedThermalError::Geometry(
            RetardedVacuumError::OutsideHistory
        ))
    );
    assert_eq!(calls, 1);
}

#[test]
fn g1_retarded_rough_thermal_uses_aberrated_material_angle_and_receiver_photons() {
    use crate::rough_diffuse::RoughDiffuseSurface;
    use crate::spectral::SpectralReflectance;
    use crate::thermal_surface::RoughDiffuseThermalSurface;

    let source = RoughDiffuseThermalSurface::try_new(
        3000.0,
        SpectralReflectance::table_in_band(&[(400.0, 0.2), (750.0, 0.8), (2000.0, 0.4)])
            .unwrap(),
        RoughDiffuseSurface::try_new(0.8).unwrap(),
        ContentHash([13; 32]),
    )
    .unwrap();
    let mut calls = 0;
    let event = RetardedThermalSurface::observe(
        &motion(0.6, [-10.0, 10.0]),
        receiver([0.0, C, 0.0]),
        1e-6,
        |proper| {
            calls += 1;
            close(proper.to_f64(), -1.0);
            Ok(ThermalEmissionState {
                surface: (&source).into(),
                outward_normal: [1.0, 1.0, 0.0],
            })
        },
    )
    .unwrap();
    assert_eq!(calls, 1);
    let cosine = 1.0 / 2.0_f64.sqrt();
    close(event.emission_cosine(), cosine);
    close(event.ray().proper_time_rate(), 1.25);
    assert_eq!(event.observed().source().identity(), source.identity());
    let sensor = QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (700.0, [1.0])]).unwrap();
    let sample = event.with_sensor(&sensor).unwrap().sample(0.5).unwrap();
    close(sample.source_wavelength_nm, 750.0);
    let native = source.directional_radiance_w_m2_sr_nm(750.0, cosine).unwrap();
    let expected = 1.25_f64.powi(5) * native;
    close(sample.spectral_radiance_w_m2_sr_nm, expected);
    close(
        sample.electron_rate_per_m2_sr[0],
        expected * 200.0 * 600e-9 / (6.626_070_15e-34 * C),
    );
    // The world ray is (.6,.8,0), while the material ray is (0,1,0).
    // Using the world angle or a hemispherical mean changes this physical law.
    let wrong_angle = source
        .directional_radiance_w_m2_sr_nm(750.0, 1.4 * cosine)
        .unwrap();
    let angular_mean = source.hemispherical_mean_radiance_w_m2_sr_nm(750.0).unwrap();
    assert!((native - wrong_angle).abs() > 1e-5 * native);
    assert!((native - angular_mean).abs() > 1e-5 * native);
}
