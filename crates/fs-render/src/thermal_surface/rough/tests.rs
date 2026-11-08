use super::*;
use crate::thermal_surface::{
    DiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
};
use core::f64::consts::PI;

fn source(rho: f64, roughness: f64) -> RoughDiffuseThermalSurface {
    RoughDiffuseThermalSurface::try_new(
        1800.0,
        SpectralReflectance::constant_in_band(rho, (300.0, 3000.0)).unwrap(),
        RoughDiffuseSurface::try_new(roughness).unwrap(),
        ContentHash([63; 32]),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g0_rough_thermal_retains_explicit_facets_band_identity_and_direction() {
    let facets =
        SpectralReflectance::table_in_band(&[(700.0, 0.2), (1200.0, 0.8), (2500.0, 0.4)]).unwrap();
    let law = RoughDiffuseSurface::try_new(0.7).unwrap();
    let surface =
        RoughDiffuseThermalSurface::try_new(1200.0, facets, law, ContentHash([61; 32])).unwrap();
    assert_eq!(surface.facet_albedo(), &facets);
    assert_eq!(surface.reflection(), law);
    assert_eq!(surface.wavelength_bounds_nm(), (700.0, 2500.0));
    assert_eq!(
        surface.directional_emissivity(699.0, 1.0),
        Err(ThermalSurfaceError::Wavelength)
    );
    assert_eq!(
        surface.directional_emissivity(1200.0, f64::NAN),
        Err(ThermalSurfaceError::Direction)
    );
    assert_eq!(
        RoughDiffuseThermalSurface::try_new(0.0, facets, law, ContentHash([61; 32])),
        Err(ThermalSurfaceError::Temperature)
    );
    assert_eq!(
        RoughDiffuseThermalSurface::try_new(1200.0, facets, law, ContentHash([0; 32])),
        Err(ThermalSurfaceError::SourceIdentity)
    );
    let changed =
        RoughDiffuseThermalSurface::try_new(1201.0, facets, law, ContentHash([61; 32])).unwrap();
    assert_ne!(surface.identity(), changed.identity());
    let view = DiffuseThermalSurface::Rough(&surface);
    assert_eq!(
        view.radiance_w_m2_sr_nm(1200.0),
        Err(ThermalSurfaceError::DirectionRequired)
    );
    assert_eq!(
        view.reflectance(1200.0),
        Err(ThermalSurfaceError::DirectionRequired)
    );
    assert_eq!(
        view.directional_radiance_w_m2_sr_nm(1200.0, 0.4).unwrap(),
        surface
            .directional_radiance_w_m2_sr_nm(1200.0, 0.4)
            .unwrap()
    );
}

#[test]
fn g1_rough_thermal_smooth_black_white_and_weak_absorption_limits() {
    for rho in [0.0, 0.2, 0.8, f64::from_bits(1.0_f64.to_bits() - 1), 1.0] {
        let smooth = source(rho, 0.0);
        let original = SpectralDiffuseThermalSurface::try_new(
            1800.0,
            SpectralEmissivity::constant_in_band(1.0 - rho, (300.0, 3000.0)).unwrap(),
            ContentHash([63; 32]),
        )
        .unwrap();
        for wavelength in [550.0, 1500.0, 2800.0] {
            for mu in [0.0, 0.02, 0.5, 1.0] {
                assert_eq!(
                    smooth
                        .directional_radiance_w_m2_sr_nm(wavelength, mu)
                        .unwrap(),
                    original.radiance_w_m2_sr_nm(wavelength).unwrap()
                );
            }
            for roughness in [0.2, 1.0] {
                let rough = source(rho, roughness);
                if rho == 1.0 {
                    assert!(!rough.is_emitting());
                    assert_eq!(rough.spectral_exitance_w_m2_nm(wavelength).unwrap(), 0.0);
                } else {
                    assert!(
                        rough
                            .directional_radiance_w_m2_sr_nm(wavelength, 0.01)
                            .unwrap()
                            > 0.0
                    );
                }
                if rho == 0.0 {
                    assert_eq!(
                        rough
                            .directional_radiance_w_m2_sr_nm(wavelength, 0.3)
                            .unwrap(),
                        planck_radiance_w_m2_sr_nm(1800.0, wavelength).unwrap()
                    );
                }
            }
        }
    }
}

const GAUSS: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];

fn integrate(f: impl Fn(f64) -> f64, lo: f64, hi: f64, panels: usize) -> f64 {
    let h = (hi - lo) / (2.0 * panels as f64);
    let mut sum = 0.0;
    for panel in 0..panels {
        let mid = lo + (2 * panel + 1) as f64 * h;
        for (x, w) in GAUSS {
            sum += w * (f(mid - h * x) + f(mid + h * x));
        }
    }
    sum * h
}

#[test]
fn g1_rough_thermal_planck_enclosure_and_exitance_match_angular_quadrature() {
    for roughness in [0.0, 0.6, 1.0] {
        for rho in [0.0, 0.3, 0.9, 1.0] {
            let surface = source(rho, roughness);
            for wavelength in [550.0, 1500.0] {
                let blackbody = planck_radiance_w_m2_sr_nm(1800.0, wavelength).unwrap();
                for mu_o in [0.02, 0.5, 1.0] {
                    let polar = |mu_i: f64| {
                        let radius = ((1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o)).sqrt();
                        let azimuth = |phi: f64| {
                            surface
                                .reflection()
                                .brdf(mu_i, mu_o, radius * phi.cos(), rho)
                                .unwrap()
                        };
                        2.0 * mu_i
                            * (integrate(azimuth, 0.0, PI / 2.0, 2)
                                + integrate(azimuth, PI / 2.0, PI, 2))
                    };
                    let reflected =
                        integrate(polar, 0.0, mu_o, 24) + integrate(polar, mu_o, 1.0, 24);
                    let emitted = surface
                        .directional_radiance_w_m2_sr_nm(wavelength, mu_o)
                        .unwrap();
                    close(reflected * blackbody + emitted, blackbody, 4e-7);
                }
                let actual = 2.0
                    * PI
                    * integrate(
                        |mu| {
                            mu * surface
                                .directional_radiance_w_m2_sr_nm(wavelength, mu)
                                .unwrap()
                        },
                        0.0,
                        1.0,
                        24,
                    );
                let expected = surface.spectral_exitance_w_m2_nm(wavelength).unwrap();
                if expected == 0.0 {
                    assert_eq!(actual, 0.0);
                } else {
                    close(actual, expected, 3e-8);
                }
            }
        }
    }
}
