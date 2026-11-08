use super::*;

#[test]
fn g0_rough_diffuse_admission_and_exact_limits() {
    for r in [-1.0, 1.01, f64::INFINITY, f64::NAN] {
        assert_eq!(
            RoughDiffuseSurface::try_new(r),
            Err(RoughDiffuseError::Roughness)
        );
    }
    assert_eq!(
        RoughDiffuseSurface::try_new(-0.0).unwrap(),
        RoughDiffuseSurface::SMOOTH
    );
    for mu in [0.001, 0.3, 1.0] {
        for rho in [0.0, 0.4, 1.0] {
            assert_eq!(
                RoughDiffuseSurface::SMOOTH.brdf(mu, mu, 0.0, rho).unwrap(),
                rho / PI
            );
            for r in [0.0, 1e-18, 0.3, 1.0] {
                let surface = RoughDiffuseSurface::try_new(r).unwrap();
                assert_eq!(surface.brdf(mu, mu, 0.0, 0.0).unwrap(), 0.0);
                assert_eq!(surface.directional_reflectance(mu, 1.0).unwrap(), 1.0);
                assert_eq!(surface.directional_reflectance(mu, 0.0).unwrap(), 0.0);
            }
        }
    }
    let surface = RoughDiffuseSurface::try_new(0.7).unwrap();
    assert_eq!(surface.brdf(0.0, 0.0, 1.0, 1.0).unwrap(), 0.0);
    assert_eq!(
        surface.brdf(1.0, 1.0, 0.5, 0.7),
        Err(RoughDiffuseError::Direction)
    );
    assert_eq!(
        surface.brdf(0.5, 0.5, 0.0, -0.1),
        Err(RoughDiffuseError::Albedo)
    );
    assert_eq!(
        surface.directional_reflectance(f64::NAN, 0.5),
        Err(RoughDiffuseError::Direction)
    );
}

#[test]
fn g3_rough_diffuse_is_reciprocal_and_enhances_backscatter() {
    for r in [0.02, 0.4, 1.0] {
        let surface = RoughDiffuseSurface::try_new(r).unwrap();
        for mu_i in [1e-12, 0.05, 0.4, 0.95, 1.0] {
            for mu_o in [1e-12, 0.2, 0.7, 1.0] {
                let bound = ((1.0_f64 - mu_i * mu_i) * (1.0 - mu_o * mu_o)).sqrt();
                for azimuth in [-1.0, -0.3, 0.0, 0.7, 1.0] {
                    let s = bound * azimuth;
                    let f = surface.brdf(mu_i, mu_o, s, 0.8).unwrap();
                    let reverse = surface.brdf(mu_o, mu_i, s, 0.8).unwrap();
                    assert_eq!(f, reverse);
                    assert!(f.is_finite() && f >= 0.0);
                }
            }
        }
        let back = surface.brdf(0.1, 0.1, 0.99, 0.6).unwrap();
        let forward = surface.brdf(0.1, 0.1, -0.99, 0.6).unwrap();
        assert!(back > 0.6 / PI && back > forward);
    }
    let tiny = angular_deficit(1e-100);
    assert!((tiny / (1e-100 / PI) - 1.0).abs() < 2e-15);
}

// Independent projected-solid-angle quadrature, split at the FON max-cosine
// and azimuth sign changes. It integrates the actual public BRDF; it does not
// reuse its analytic directional-albedo or multiple-scattering coefficients.
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
        for (node, weight) in GAUSS {
            sum += weight * (f(mid - h * node) + f(mid + h * node));
        }
    }
    sum * h
}

fn integrated_reflectance(surface: RoughDiffuseSurface, mu_o: f64, rho: f64) -> f64 {
    let polar = |mu_i: f64| {
        let radius = ((1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o)).sqrt();
        let azimuth = |phi: f64| surface.brdf(mu_i, mu_o, radius * phi.cos(), rho).unwrap();
        2.0 * mu_i * (integrate(azimuth, 0.0, PI / 2.0, 2) + integrate(azimuth, PI / 2.0, PI, 2))
    };
    integrate(polar, 0.0, mu_o, 24) + integrate(polar, mu_o, 1.0, 24)
}

#[test]
fn g1_rough_diffuse_white_furnace_and_colored_energy_match_projected_quadrature() {
    for r in [0.0, 0.02, 0.4, 1.0] {
        let surface = RoughDiffuseSurface::try_new(r).unwrap();
        for mu in [1e-6, 0.05, 0.5, 1.0] {
            for rho in [0.1, 0.7, 1.0] {
                let expected = surface.directional_reflectance(mu, rho).unwrap();
                let actual = integrated_reflectance(surface, mu, rho);
                assert!(
                    (actual - expected).abs() < 3e-7,
                    "r={r}, mu={mu}, rho={rho}: integrated={actual}, analytic={expected}"
                );
                assert!(expected >= 0.0 && expected <= rho);
                if rho == 1.0 {
                    assert!((actual - 1.0).abs() < 3e-7);
                }
            }
        }
    }
}

#[test]
fn g1_kirchhoff_absorptance_closes_directional_and_hemispherical_power() {
    for r in [0.0, 0.02, 0.4, 1.0] {
        let surface = RoughDiffuseSurface::try_new(r).unwrap();
        for rho in [0.0, 0.1, 0.7, 0.99, 1.0] {
            for mu in [0.0, 1e-12, 0.05, 0.5, 1.0] {
                let reflection = surface.directional_reflectance(mu, rho).unwrap();
                let absorption = surface.directional_absorptance(mu, rho).unwrap();
                assert!((reflection + absorption - 1.0).abs() <= 2.0 * f64::EPSILON);
                assert!((0.0..=1.0).contains(&absorption));
                if mu > 0.0 && rho > 0.0 && rho < 1.0 {
                    // Independent angular BRDF quadrature closes a unit
                    // isotropic bath when emission follows Kirchhoff's law.
                    let reflected_bath = integrated_reflectance(surface, mu, rho);
                    assert!((reflected_bath + absorption - 1.0).abs() < 3e-7);
                }
            }
            let mean = surface.hemispherical_absorptance(rho).unwrap();
            let angular = integrate(
                |mu| 2.0 * mu * surface.directional_absorptance(mu, rho).unwrap(),
                0.0,
                1.0,
                64,
            );
            assert!((angular - mean).abs() <= 2e-10 * mean.max(f64::MIN_POSITIVE));
        }
    }
}

#[test]
fn g0_absorptance_retains_weak_loss_and_exact_limits() {
    let almost_white = f64::from_bits(1.0_f64.to_bits() - 1);
    let facet_loss = 1.0 - almost_white;
    for r in [0.0, 1e-18, 0.4, 1.0] {
        let surface = RoughDiffuseSurface::try_new(r).unwrap();
        for mu in [0.0, 1e-100, 0.1, 0.7, 1.0] {
            assert_eq!(surface.directional_absorptance(mu, 0.0).unwrap(), 1.0);
            assert_eq!(surface.directional_absorptance(mu, 1.0).unwrap(), 0.0);
            let loss = surface.directional_absorptance(mu, almost_white).unwrap();
            assert!(loss.is_finite() && loss >= facet_loss);
            if r == 0.0 {
                assert_eq!(loss, facet_loss);
            }
        }
        assert_eq!(surface.hemispherical_absorptance(0.0).unwrap(), 1.0);
        assert_eq!(surface.hemispherical_absorptance(1.0).unwrap(), 0.0);
        let mean = surface.hemispherical_absorptance(almost_white).unwrap();
        assert!(mean >= facet_loss);
        if r >= 0.4 {
            // The absorbed energy from extra returns is representable even
            // when subtracting the rounded reflected fraction would lose it.
            assert!(mean / facet_loss > 1.05);
        }
    }
    for invalid in [-0.1, 1.1, f64::INFINITY, f64::NAN] {
        assert_eq!(
            RoughDiffuseSurface::SMOOTH.directional_absorptance(invalid, 0.4),
            Err(RoughDiffuseError::Direction)
        );
        assert_eq!(
            RoughDiffuseSurface::SMOOTH.directional_absorptance(0.4, invalid),
            Err(RoughDiffuseError::Albedo)
        );
        assert_eq!(
            RoughDiffuseSurface::SMOOTH.hemispherical_absorptance(invalid),
            Err(RoughDiffuseError::Albedo)
        );
    }
}
