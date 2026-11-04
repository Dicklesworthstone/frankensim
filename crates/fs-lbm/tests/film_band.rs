//! G1 continuous-band extrema against independent dimensional coefficients.

use fs_lbm::film::{FilmError, FilmMaterial, LongWaveLimits, NusseltFilm};

const LIQUID: FilmMaterial = FilmMaterial {
    density_kg_m3: 1000.0,
    kinematic_viscosity_m2_s: 1e-6,
    surface_tension_n_m: 0.072,
};
const LIMITS: LongWaveLimits = LongWaveLimits {
    max_kh: 0.1,
    max_re_kh: 0.1,
};

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 4e-12 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn continuous_band_finds_instability_missed_by_endpoint_samples() {
    let film = NusseltFilm::new(LIQUID, 1e-6, 4.0, 4.0).unwrap();
    assert_eq!(film.long_wave_mode(0.0, LIMITS).unwrap().growth_per_s, 0.0);
    assert!(film.long_wave_mode(200.0, LIMITS).unwrap().growth_per_s < 0.0);
    let band = film.long_wave_band([0.0, 200.0], LIMITS).unwrap();
    let h = film.thickness_m();
    let nu = LIQUID.kinematic_viscosity_m2_s;
    // Derive the maximum from independent dimensional h-power coefficients.
    let a = 2.0 * 4.0_f64.powi(2) * h.powi(6) / (15.0 * nu.powi(3)) - 4.0 * h.powi(3) / (3.0 * nu);
    let b = LIQUID.surface_tension_n_m * h.powi(3) / (3.0 * LIQUID.density_kg_m3 * nu);
    near(band.wavenumber_rad_m, (a / (2.0 * b)).sqrt());
    near(band.mode.growth_per_s, a * a / (4.0 * b));
    assert!(band.mode.growth_per_s > 0.0);
    near(band.max_kh, 200.0 * h);
    for i in 0..=1000 {
        let sampled = film.long_wave_mode(i as f64 * 0.2, LIMITS).unwrap();
        assert!(sampled.growth_per_s <= band.mode.growth_per_s + 1e-15);
    }
}

#[test]
fn endpoint_maxima_capillarity_omission_and_neutral_ties_are_physical() {
    let film = NusseltFilm::new(LIQUID, 1e-6, 4.0, 4.0).unwrap();
    let rising = film.long_wave_band([1.0, 20.0], LIMITS).unwrap();
    assert_eq!(rising.wavenumber_rad_m, 20.0);
    let falling = film.long_wave_band([120.0, 180.0], LIMITS).unwrap();
    assert_eq!(falling.wavenumber_rad_m, 120.0);
    let stable = NusseltFilm::new(LIQUID, 1e-7, 4.0, 4.0).unwrap();
    assert_eq!(
        stable
            .long_wave_band([0.0, 180.0], LIMITS)
            .unwrap()
            .wavenumber_rad_m,
        0.0
    );
    let no_capillary = NusseltFilm::new(
        FilmMaterial {
            surface_tension_n_m: 0.0,
            ..LIQUID
        },
        1e-6,
        4.0,
        4.0,
    )
    .unwrap();
    assert_eq!(
        no_capillary
            .long_wave_band([0.0, 180.0], LIMITS)
            .unwrap()
            .wavenumber_rad_m,
        180.0
    );
    let point = film.long_wave_band([10.0, 10.0], LIMITS).unwrap();
    assert_eq!(point.mode, film.long_wave_mode(10.0, LIMITS).unwrap());
    let zero = film.long_wave_band([0.0, 0.0], LIMITS).unwrap();
    assert_eq!(zero.mode.growth_per_s, 0.0);
}

#[test]
fn unsupported_band_tail_refuses_even_when_maximum_is_admitted() {
    let film = NusseltFilm::new(LIQUID, 1e-6, 4.0, 4.0).unwrap();
    for interval in [
        [2.0, 1.0],
        [-1.0, 20.0],
        [0.0, f64::INFINITY],
        [f64::NAN, 2.0],
    ] {
        assert_eq!(
            film.long_wave_band(interval, LIMITS),
            Err(FilmError::InvalidInput("wavenumber band"))
        );
    }
    let outside = 1.01 * LIMITS.max_kh / film.thickness_m();
    assert_eq!(
        film.long_wave_band([0.0, outside], LIMITS),
        Err(FilmError::OutsideLongWaveRegime)
    );
    assert!(film.long_wave_band([0.0, 200.0], LIMITS).is_ok());
}
