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
