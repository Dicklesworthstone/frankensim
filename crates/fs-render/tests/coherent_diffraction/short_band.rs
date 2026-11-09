//! G1 fourth-order temporal nulls through the original scalar pupil;
//! G3 authored covariance identity and G0/G4 exact modal work admission.
use super::*;
use fs_math::dd::Dd;

const FREQUENCIES: [f64; 3] = [5.0e14, 5.0e14 + 1.0, 5.0e14 + 2.0];
const SENSOR_M: f64 = 0.03;

// The clear pupil has a chromatic on-axis transfer. Two authored waves retain
// the DD inverse transfer so a rounded input coefficient cannot mask the
// physical T^4 null. No detector or temporal kernel is used in this reference.
fn normalized_waves(camera: &PhysicalCamera) -> [[CoherentPupilWave; 2]; 3] {
    std::array::from_fn(|i| {
        let wavelength_nm = (Dd::from_f64(299_792_458.0) / Dd::from_f64(FREQUENCIES[i])
            * Dd::from_f64(1.0e9))
        .to_f64();
        let unit = [wave([0.0; 2], C64::ONE)];
        let transfer = camera
            .prepare_coherent_diffraction(
                wavelength_nm,
                CoherentPupilIllumination {
                    mean: &unit,
                    factors: &[],
                },
                LIMITS,
                || false,
            )
            .unwrap()
            .observe([0.0; 2], || false)
            .unwrap()
            .mean_field;
        assert_eq!(transfer.im, 0.0);
        let amplitude = Dd::from_f64([1.0, -2.0, 1.0][i]) / Dd::from_f64(transfer.re);
        [
            wave([0.0; 2], C64::new(amplitude.hi, 0.0)),
            wave([0.0; 2], C64::new(amplitude.lo, 0.0)),
        ]
    })
}

fn null_shutter(duration: f64) -> CoherentDiffractionShutter {
    CoherentDiffractionShutter::try_new(SENSOR_M / 299_792_458.0 - duration * 0.5, duration)
        .unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual > 0.0 && (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_three_tone_scalar_shutter_retains_fourth_order_dark_signal() {
    let camera = focused();
    let waves = normalized_waves(&camera);
    let tones = std::array::from_fn::<_, 3, _>(|i| tone(FREQUENCIES[i], &waves[i], &[]));
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let query = camera
            .prepare_coherent_diffraction_exposure(
                &tones,
                SENSOR_M,
                null_shutter(duration),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let sample = query.observe([0.0; 2], || false).unwrap();
        // |1-2*exp(-i*2*pi*t)+exp(-i*4*pi*t)|^2 = 16*sin(pi*t)^4.
        // Its symmetric uniform mean is pi^4*T^4/5 + O(T^6).
        let expected = core::f64::consts::PI.powi(4) * duration.powi(4) / 5.0;
        relative(sample.mean_power_density_w, expected, 2e-6);
        assert_eq!(sample.fluctuation_power_density_w, 0.0);
        assert_eq!(sample.power_density_w, sample.mean_power_density_w);
    }
}

#[test]
fn g3_short_band_keeps_shared_and_independent_covariance_columns_distinct() {
    let camera = focused();
    let waves = normalized_waves(&camera);
    let shared = waves
        .each_ref()
        .map(|waves| [PupilCoherenceFactor { waves }]);
    let shared_tones = std::array::from_fn::<_, 3, _>(|i| tone(FREQUENCIES[i], &[], &shared[i]));
    let independent = std::array::from_fn::<_, 3, _>(|i| {
        std::array::from_fn::<_, 3, _>(|column| PupilCoherenceFactor {
            waves: if i == column { &waves[i] } else { &[] },
        })
    });
    let independent_tones =
        std::array::from_fn::<_, 3, _>(|i| tone(FREQUENCIES[i], &[], &independent[i]));
    let shutter = null_shutter(1e-9);
    let shared = camera
        .prepare_coherent_diffraction_exposure(
            &shared_tones,
            SENSOR_M,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap();
    let independent = camera
        .prepare_coherent_diffraction_exposure(
            &independent_tones,
            SENSOR_M,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap();
    relative(
        shared.fluctuation_power_density_w,
        core::f64::consts::PI.powi(4) * 1e-36 / 5.0,
        2e-6,
    );
    relative(independent.fluctuation_power_density_w, 6.0, 2e-14);
    assert_eq!(shared.mean_power_density_w, 0.0);
    assert_eq!(independent.mean_power_density_w, 0.0);
}

#[test]
fn g0_g4_short_band_preflights_preparation_pixel_work_and_cancelled_reuse() {
    let camera = focused();
    let waves = normalized_waves(&camera);
    let tones = std::array::from_fn::<_, 3, _>(|i| tone(FREQUENCIES[i], &waves[i], &[]));
    // Three modes, three projection terms and one norm in each mode.
    let work = 3 * (3 + 1);
    let limits = CoherentDiffractionExposureLimits {
        maximum_pair_contractions: work,
        ..EXPOSURE_LIMITS
    };
    assert!(matches!(
        camera.prepare_coherent_diffraction_exposure(
            &tones,
            SENSOR_M,
            null_shutter(1e-9),
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: work - 1,
                ..limits
            },
            || false,
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "exposure pair contractions",
            requested: 12,
            ..
        })
    ));
    let query = camera
        .prepare_coherent_diffraction_exposure(&tones, SENSOR_M, null_shutter(1e-9), limits, || {
            false
        })
        .unwrap();
    let footprint = pixel([-0.1; 2], [0.1; 2]);
    let sample = query
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    assert!(sample.energy_j > 0.0);
    assert_eq!(
        sample.pair_contractions,
        work * sample.average_power.quadrature_samples
    );
    relative(sample.energy_j, sample.average_power.power_w * 1e-9, 2e-14);
    let exact = CoherentDiffractionExposurePixelLimits {
        maximum_pair_contractions: sample.pair_contractions,
        ..EXPOSURE_PIXEL_LIMITS
    };
    assert_eq!(
        query.integrate_pixel(footprint, exact, || false).unwrap(),
        sample
    );
    let mut polls = 0;
    assert!(matches!(
        query.integrate_pixel(
            footprint,
            CoherentDiffractionExposurePixelLimits {
                maximum_pair_contractions: sample.pair_contractions - 1,
                ..exact
            },
            || {
                polls += 1;
                false
            }
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pixel exposure pair contractions",
            ..
        })
    ));
    assert_eq!(polls, 1);
    let mut polls = 0;
    assert_eq!(
        query.observe([0.0; 2], || {
            polls += 1;
            polls == 8
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    relative(
        query.observe([0.0; 2], || false).unwrap().power_density_w,
        core::f64::consts::PI.powi(4) * 1e-36 / 5.0,
        2e-6,
    );
}
