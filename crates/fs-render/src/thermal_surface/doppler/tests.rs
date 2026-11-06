use super::*;
use crate::spectral::doppler::VACUUM_LIGHT_SPEED_M_S;
use crate::thermal_surface::{
    GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
};
use fs_blake3::ContentHash;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn approaching(beta: f64) -> VacuumDoppler {
    VacuumDoppler::try_new(
        [1.0, 0.0, 0.0],
        [beta * VACUUM_LIGHT_SPEED_M_S, 0.0, 0.0],
        [0.0; 3],
    )
    .unwrap()
}

fn gray(temperature: f64, emissivity: f64) -> GrayDiffuseThermalSurface {
    GrayDiffuseThermalSurface::try_new(temperature, emissivity, ContentHash([1; 32])).unwrap()
}

// Independent direct SI expression, used only where its exponent is moderate.
fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let lambda = wavelength_nm * 1.0e-9;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    2.0 * h * c * c / (lambda.powi(5) * (h * c / (lambda * k * temperature)).exp_m1()) * 1.0e-9
}

#[test]
fn g1_gray_color_temperature_and_d5_preserve_thermal_owner() {
    let source = gray(3000.0, 0.7);
    let observed = ObservedThermalSurface::try_new(source.into(), approaching(0.6)).unwrap();
    let radiance = observed.radiance_w_m2_sr_nm(600.0).unwrap();
    close(radiance, 0.7 * planck(6000.0, 600.0), 2.0e-13);
    close(
        radiance,
        32.0 * source.radiance_w_m2_sr_nm(1200.0).unwrap(),
        2.0e-13,
    );
    assert_eq!(observed.source().temperature_k(), 3000.0);
    assert_eq!(observed.source().identity(), source.identity());
    assert_eq!(observed.receiver_wavelength_bounds_nm(), None);
}

#[test]
fn g1_selective_emissivity_uses_source_wavelength_and_retains_support() {
    let emissivity =
        SpectralEmissivity::table_in_band(&[(800.0, 0.1), (1200.0, 0.9), (1600.0, 0.2)]).unwrap();
    let source =
        SpectralDiffuseThermalSurface::try_new(3000.0, emissivity, ContentHash([2; 32])).unwrap();
    let observed = ObservedThermalSurface::try_new((&source).into(), approaching(0.6)).unwrap();
    close(
        observed.radiance_w_m2_sr_nm(600.0).unwrap(),
        0.9 * planck(6000.0, 600.0),
        2.0e-13,
    );
    let (low, high) = observed.receiver_wavelength_bounds_nm().unwrap();
    close(low, 400.0, 2.0e-15);
    close(high, 800.0, 2.0e-15);
    assert!(observed.radiance_w_m2_sr_nm(low).is_ok());
    assert!(observed.radiance_w_m2_sr_nm(high).is_ok());
    assert_eq!(
        observed.radiance_w_m2_sr_nm(f64::from_bits(low.to_bits() - 1)),
        Err(ThermalDopplerError::SourceBand)
    );
    assert_eq!(
        observed.radiance_w_m2_sr_nm(f64::from_bits(high.to_bits() + 1)),
        Err(ThermalDopplerError::SourceBand)
    );
    assert_eq!(observed.source().identity(), source.identity());
}

#[test]
fn g1_detector_uses_receiver_proposal_and_receiver_photon_energy() {
    let source = gray(3000.0, 0.7);
    let observed = ObservedThermalSurface::try_new(source.into(), approaching(0.6)).unwrap();
    let sensor = QuantumEfficiencySensor::try_new(vec![
        (500.0, [1.0, 0.25, 0.0]),
        (700.0, [1.0, 0.25, 0.0]),
    ])
    .unwrap();
    let detector = observed.with_sensor(&sensor).unwrap();
    let sample = detector.sample(0.5).unwrap();
    close(sample.receiver_wavelength_nm, 600.0, 2.0e-15);
    close(sample.source_wavelength_nm, 1200.0, 2.0e-15);
    let radiance = 0.7 * planck(6000.0, 600.0);
    // Uniform RECEIVER-frame pdf=1/200 nm^-1, independent of source D.
    close(
        sample.qe_weighted_radiance_w_m2_sr[0],
        radiance * 200.0,
        2.0e-13,
    );
    close(
        sample.qe_weighted_radiance_w_m2_sr[1],
        radiance * 50.0,
        2.0e-13,
    );
    close(
        sample.electron_rate_per_m2_sr[0],
        radiance * 200.0 * 600.0e-9 / (6.626_070_15e-34 * VACUUM_LIGHT_SPEED_M_S),
        2.0e-13,
    );
    assert_eq!(sample.qe_weighted_radiance_w_m2_sr[2], 0.0);
    assert_eq!(sample.electron_rate_per_m2_sr[2], 0.0);
    assert_eq!(detector.sample(0.5).unwrap(), sample);
    assert_eq!(
        detector.sample(f64::NAN),
        Err(ThermalDopplerError::Sensor(
            SpectralSensorError::InvalidVariate
        ))
    );
}

#[test]
fn g1_nonuniform_detector_proposal_matches_direct_observed_band_integral() {
    let source = gray(2800.0, 0.8);
    let observed = ObservedThermalSurface::try_new(source.into(), approaching(0.6)).unwrap();
    let sensor = QuantumEfficiencySensor::try_new(vec![(450.0, [0.1]), (750.0, [0.9])]).unwrap();
    let detector = observed.with_sensor(&sensor).unwrap();
    let count = 8192;
    let mut sampled = 0.0;
    let mut direct = 0.0;
    for i in 0..count {
        let u = (i as f64 + 0.5) / count as f64;
        sampled += detector.sample(u).unwrap().qe_weighted_radiance_w_m2_sr[0];
        let wavelength = 450.0 + 300.0 * u;
        direct += 0.8 * planck(5600.0, wavelength) * (0.1 + 0.8 * u) * 300.0;
    }
    close(sampled / count as f64, direct / count as f64, 1.0e-7);
}

#[test]
fn g3_full_detector_band_is_admitted_before_sampling() {
    let source = SpectralDiffuseThermalSurface::try_new(
        3000.0,
        SpectralEmissivity::constant_in_band(0.5, (1000.0, 1400.0)).unwrap(),
        ContentHash([3; 32]),
    )
    .unwrap();
    let observed = ObservedThermalSurface::try_new((&source).into(), approaching(0.6)).unwrap();
    let valid = QuantumEfficiencySensor::try_new(vec![(510.0, [1.0]), (690.0, [1.0])]).unwrap();
    assert!(observed.with_sensor(&valid).is_ok());
    let invalid = QuantumEfficiencySensor::try_new(vec![(490.0, [1.0]), (600.0, [1.0])]).unwrap();
    assert!(matches!(
        observed.with_sensor(&invalid),
        Err(ThermalDopplerError::SourceBand)
    ));
}

#[test]
fn g3_zero_shift_preserves_existing_radiance_bits_and_dark_source() {
    let source = gray(3000.0, 0.7);
    let observed = ObservedThermalSurface::try_new(source.into(), approaching(0.0)).unwrap();
    assert_eq!(
        observed.radiance_w_m2_sr_nm(601.234).unwrap().to_bits(),
        source.radiance_w_m2_sr_nm(601.234).unwrap().to_bits()
    );
    let dark = ObservedThermalSurface::try_new(gray(3000.0, 0.0).into(), approaching(0.6)).unwrap();
    assert_eq!(dark.radiance_w_m2_sr_nm(600.0).unwrap(), 0.0);
}

#[test]
fn g3_amplified_wien_tail_is_evaluated_before_rest_frame_underflow() {
    let source = gray(1.0, 1.0);
    let shift = approaching(1.0 - 2.0e-14);
    let source_nm = 18000.0;
    assert_eq!(source.radiance_w_m2_sr_nm(source_nm).unwrap(), 0.0);
    let observed = ObservedThermalSurface::try_new(source.into(), shift).unwrap();
    let radiance = observed
        .radiance_w_m2_sr_nm(shift.receiver_wavelength_nm(source_nm).unwrap())
        .unwrap();
    assert!(radiance.is_finite() && radiance > 0.0);
}
