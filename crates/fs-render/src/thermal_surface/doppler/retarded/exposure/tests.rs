use super::*;
use crate::spectral::doppler::{VACUUM_LIGHT_SPEED_M_S as LIGHT, VacuumDoppler};
use crate::thermal_surface::{
    GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
};
use fs_blake3::ContentHash;

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}

fn motion(x: f64, beta: f64) -> UniformVacuumWorldline {
    UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [x * LIGHT, 0.0, 0.0],
        [beta * LIGHT, 0.0, 0.0],
        [dd(-20.0), dd(20.0)],
    )
    .unwrap()
}

fn settings(subdivisions: usize) -> ThermalExposureSettings {
    ThermalExposureSettings {
        time: RetardedExposureSettings {
            max_ray_queries: 256,
            max_residual_m: 1e-6,
        },
        receiver_etendue_m2_sr: 1e-8,
        wavelength_subdivisions: subdivisions,
        max_spectral_evaluations: 100_000,
    }
}

fn gray(temperature: f64) -> ThermalEmissionState<'static> {
    ThermalEmissionState {
        surface: GrayDiffuseThermalSurface::try_new(temperature, 0.7, ContentHash([31; 32]))
            .unwrap()
            .into(),
        outward_normal: [1.0, 0.0, 0.0],
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let k = 1.380_649e-23;
    let wavelength = wavelength_nm * 1e-9;
    2.0 * h * LIGHT * LIGHT
        / (wavelength.powi(5) * (h * LIGHT / (wavelength * k * temperature)).exp_m1())
        * 1e-9
}

// Independent uniform Simpson quadrature and direct SI Planck expression.
fn reference(temperature: f64, band: (f64, f64)) -> [f64; 2] {
    let n = 16_384;
    let step = (band.1 - band.0) / f64::from(n);
    let mut sum = [0.0; 2];
    for i in 0..=n {
        let wavelength = band.0 + f64::from(i) * step;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let value = 0.7 * planck(temperature, wavelength);
        sum[0] += weight * value;
        sum[1] += weight * value * wavelength * 1e-9 / (6.626_070_15e-34 * LIGHT);
    }
    sum.map(|x| x * step / 3.0)
}

#[test]
fn g1_moving_thermal_exposure_matches_absolute_si_energy_and_photon_integrals() {
    let sensor =
        QuantumEfficiencySensor::try_new(vec![(500.0, [1.0, 0.0]), (800.0, [1.0, 0.0])]).unwrap();
    let mut calls = 0;
    let result = RetardedThermalExposure::uniform(
        &motion(0.0, 0.6),
        &motion(3.0, 0.0),
        &sensor,
        &[dd(0.0), dd(0.25)],
        settings(8),
        |_| {
            calls += 1;
            Ok(gray(3000.0))
        },
        || false,
    )
    .unwrap();
    let expected = reference(6000.0, (500.0, 800.0)).map(|x| x * 0.25 * 1e-8);
    close(result.qe_weighted_energy_j[0], expected[0], 2e-10);
    close(result.expected_electrons[0], expected[1], 2e-10);
    assert_eq!(result.qe_weighted_energy_j[1], 0.0);
    assert_eq!(result.expected_electrons[1], 0.0);
    assert_eq!(result.ray_queries, 6);
    assert_eq!(result.spectral_evaluations, 128);
    assert_eq!(calls, 4);
}

#[test]
fn g1_shifted_narrow_emissivity_and_qe_features_survive_the_exposure_integral() {
    let source = SpectralDiffuseThermalSurface::try_new(
        3000.0,
        SpectralEmissivity::table_in_band(&[
            (800.0, 0.0),
            (1001.0, 0.0),
            (1001.002, 1.0),
            (1001.004, 0.0),
            (1600.0, 0.0),
        ])
        .unwrap(),
        ContentHash([32; 32]),
    )
    .unwrap();
    let d = VacuumDoppler::try_new([1.0, 0.0, 0.0], [0.6 * LIGHT, 0.0, 0.0], [0.0; 3])
        .unwrap()
        .frequency_ratio();
    let center = 1001.002 / d;
    let sensor = QuantumEfficiencySensor::try_new(vec![
        (450.0, [1.0, 0.0]),
        (1001.001 / d, [1.0, 0.0]),
        (center, [1.0, 1.0]),
        (1001.003 / d, [1.0, 0.0]),
        (750.0, [1.0, 0.0]),
    ])
    .unwrap();
    let result = RetardedThermalExposure::uniform(
        &motion(0.0, 0.6),
        &motion(3.0, 0.0),
        &sensor,
        &[dd(0.0), dd(0.2)],
        settings(1),
        |_| {
            Ok(ThermalEmissionState {
                surface: (&source).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
        || false,
    )
    .unwrap();
    // Across this picometre-scale half-width, Planck varies smoothly. The exact
    // triangular areas are .002/D and (5/12)*.002/D; curvature is negligible
    // at this tolerance. Omitting shifted source knots loses its outer wings.
    let area = 0.002 / d;
    let energy = planck(3000.0 * d, center) * area * 0.2 * 1e-8;
    close(result.qe_weighted_energy_j[0], energy, 3e-8);
    close(result.qe_weighted_energy_j[1], energy * (5.0 / 12.0), 3e-8);
    close(
        result.expected_electrons[0],
        energy * center * PHOTONS_PER_JOULE_PER_NM,
        3e-8,
    );
    assert!(result.spectral_evaluations > 16);
}

#[test]
fn g1_broad_band_refinement_converges_without_changing_physical_qe_knots() {
    let band = (200.0, 5000.0);
    let sensor = QuantumEfficiencySensor::try_new(vec![(band.0, [1.0]), (band.1, [1.0])]).unwrap();
    let evaluate = |subdivisions| {
        RetardedThermalExposure::uniform(
            &motion(0.0, 0.0),
            &motion(1.0, 0.0),
            &sensor,
            &[dd(1.0), dd(2.0)],
            settings(subdivisions),
            |_| Ok(gray(2000.0)),
            || false,
        )
        .unwrap()
    };
    let coarse = evaluate(1);
    let fine = evaluate(32);
    let expected = reference(2000.0, band).map(|x| x * 1e-8);
    let coarse_error = (coarse.expected_electrons[0] - expected[1]).abs();
    let fine_error = (fine.expected_electrons[0] - expected[1]).abs();
    assert!(fine_error < coarse_error * 1e-4);
    close(fine.expected_electrons[0], expected[1], 2e-8);
    close(fine.qe_weighted_energy_j[0], expected[0], 2e-8);
    assert_eq!(fine, evaluate(32));
}

#[test]
fn g4_spectral_work_dark_etendue_and_cancellation_keep_original_admission() {
    let sensor = QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (700.0, [1.0])]).unwrap();
    let mut limited = settings(1);
    limited.max_spectral_evaluations = 15;
    let evaluate = |settings| {
        RetardedThermalExposure::uniform(
            &motion(0.0, 0.0),
            &motion(1.0, 0.0),
            &sensor,
            &[dd(1.0), dd(2.0)],
            settings,
            |_| Ok(gray(3000.0)),
            || false,
        )
    };
    assert_eq!(
        evaluate(limited),
        Err(RetardedExposureError::Observable(
            ThermalExposureError::SpectralBudget
        ))
    );
    limited.max_spectral_evaluations = 16;
    assert_eq!(evaluate(limited).unwrap().spectral_evaluations, 16);
    limited.receiver_etendue_m2_sr = 0.0;
    let dark = evaluate(limited).unwrap();
    assert_eq!(dark.expected_electrons, [0.0]);
    assert_eq!(dark.qe_weighted_energy_j, [0.0]);
    for etendue in [-1.0, f64::NAN, f64::INFINITY] {
        limited.receiver_etendue_m2_sr = etendue;
        assert_eq!(
            evaluate(limited),
            Err(RetardedExposureError::Observable(
                ThermalExposureError::Etendue
            ))
        );
    }
    let mut gate = 0;
    let result = RetardedThermalExposure::uniform(
        &motion(0.0, 0.0),
        &motion(1.0, 0.0),
        &sensor,
        &[dd(1.0), dd(2.0)],
        settings(1),
        |_| Ok(gray(3000.0)),
        || {
            gate += 1;
            gate >= 5
        },
    );
    assert_eq!(
        result,
        Err(RetardedExposureError::Observable(
            ThermalExposureError::Cancelled
        ))
    );

    let finite = SpectralDiffuseThermalSurface::try_new(
        3000.0,
        SpectralEmissivity::constant_in_band(0.5, (800.0, 1200.0)).unwrap(),
        ContentHash([33; 32]),
    )
    .unwrap();
    let result = RetardedThermalExposure::uniform(
        &motion(0.0, 0.0),
        &motion(1.0, 0.0),
        &sensor,
        &[dd(1.0), dd(2.0)],
        ThermalExposureSettings {
            receiver_etendue_m2_sr: 0.0,
            ..settings(1)
        },
        |_| {
            Ok(ThermalEmissionState {
                surface: (&finite).into(),
                outward_normal: [1.0, 0.0, 0.0],
            })
        },
        || false,
    );
    assert_eq!(
        result,
        Err(RetardedExposureError::Observable(
            ThermalExposureError::Detector(ThermalDopplerError::SourceBand)
        ))
    );
}
