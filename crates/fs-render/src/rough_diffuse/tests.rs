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
fn g1_reflection_proposal_has_full_support_and_normalized_directional_mass() {
    for roughness in [0.0, 0.25, 1.0] {
        let surface = RoughDiffuseSurface::try_new(roughness).unwrap();
        for mu_o in [1e-6, 0.1, 0.6, 1.0] {
            let angular = |mu_i: f64| {
                let radius = ((1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o)).sqrt();
                let pdf = |phi: f64| surface.pdf_admitted(mu_o, mu_i, radius * phi.cos());
                2.0 * (integrate(pdf, 0.0, PI / 2.0, 4) + integrate(pdf, PI / 2.0, PI, 4))
            };
            assert!((integrate(angular, 0.0, 1.0, 4) - 1.0).abs() < 3e-14);
            for mu_i in [1e-12, 0.2, 0.9] {
                let radius = ((1.0_f64 - mu_i * mu_i) * (1.0 - mu_o * mu_o)).sqrt();
                // Backward azimuths retain positive cosine mass even when
                // the added retroreflection proposal has no support there.
                assert!(surface.pdf_admitted(mu_o, mu_i, -radius) > 0.0);
                if roughness == 0.0 || mu_o == 1.0 {
                    assert_eq!(surface.pdf_admitted(mu_o, mu_i, radius), mu_i / PI);
                }
            }
        }
    }
    let surface = RoughDiffuseSurface::try_new(1.0).unwrap();
    assert!((surface.pdf_admitted(0.1, 0.8, 0.2)
        - surface.pdf_admitted(0.8, 0.1, 0.2)).abs() > 0.05);
}

#[test]
fn g1_reflection_samples_match_distribution_moments_and_colored_furnace() {
    const RADIAL_COUNT: usize = 1024;
    const AZIMUTH_COUNT: usize = 128;
    for (roughness, mu_o) in [(0.0, 0.2), (0.25, 0.1), (1.0, 0.01), (1.0, 1.0)] {
        let surface = RoughDiffuseSurface::try_new(roughness).unwrap();
        let probability = surface.grazing_proposal_probability(mu_o);
        let mut direction_mean = [0.0; 3];
        let mut power_mean = [0.0; 2];
        for row in 0..RADIAL_COUNT {
            for column in 0..AZIMUTH_COUNT {
                let wi = surface.sample_local(mu_o,
                    [(row as f64 + 0.5) / RADIAL_COUNT as f64,
                        (column as f64 + 0.5) / AZIMUTH_COUNT as f64]);
                assert!((wi.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 2e-14);
                assert!(wi[2] > 0.0);
                let projected = (1.0 - mu_o * mu_o).sqrt() * wi[0];
                let pdf = surface.pdf_admitted(mu_o, wi[2], projected);
                assert!(pdf.is_finite() && pdf > 0.0);
                for axis in 0..3 {
                    direction_mean[axis] += wi[axis] / (RADIAL_COUNT * AZIMUTH_COUNT) as f64;
                }
                for (lane, rho) in [0.4, 1.0].into_iter().enumerate() {
                    let weight = surface.brdf(wi[2], mu_o, projected, rho).unwrap() * wi[2] / pdf;
                    assert!(weight.is_finite() && weight >= 0.0);
                    power_mean[lane] += weight / (RADIAL_COUNT * AZIMUTH_COUNT) as f64;
                }
            }
        }
        // Independent analytic moments of cosine and uniform-elevation /
        // cosine-azimuth distributions, not reconstructed sampling formulas.
        assert!((direction_mean[0] - probability * PI * PI / 16.0).abs() < 0.0015);
        assert!(direction_mean[1].abs() < 2e-14);
        assert!((direction_mean[2] - ((1.0-probability)*2.0/3.0 + probability*0.5)).abs() < 0.0015);
        for (lane, rho) in [0.4, 1.0].into_iter().enumerate() {
            let reference = surface.directional_reflectance(mu_o, rho).unwrap();
            assert!((power_mean[lane] - reference).abs() < 0.003,
                "roughness={roughness}, cosine={mu_o}, albedo={rho}: {} vs {reference}", power_mean[lane]);
        }
    }
}

fn proposal_second_moment(surface: RoughDiffuseSurface, mu_o: f64, rho: f64,
    mixture: bool, panels: usize) -> f64 {
    let polar = |mu_i: f64| {
        let radius = ((1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o)).sqrt();
        let angular = |phi: f64| {
            let projected = radius * phi.cos();
            let flux = surface.brdf(mu_i, mu_o, projected, rho).unwrap() * mu_i;
            let pdf = if mixture { surface.pdf_admitted(mu_o, mu_i, projected) } else { mu_i / PI };
            flux * flux / pdf
        };
        2.0 * (integrate(angular, 0.0, PI/2.0, 8) + integrate(angular, PI/2.0, PI, 8))
    };
    integrate(polar, 0.0, mu_o, panels) + integrate(polar, mu_o, 1.0, panels)
}

#[test]
fn g2_equal_path_count_grazing_reflection_variance_is_reduced() {
    let surface = RoughDiffuseSurface::try_new(1.0).unwrap();
    for rho in [0.4, 1.0] {
        for mu_o in [0.1, 0.01] {
            let mean = integrated_reflectance(surface, mu_o, rho);
            let reference = surface.directional_reflectance(mu_o, rho).unwrap();
            assert!((mean - reference).abs() < 3e-7);
            let old = proposal_second_moment(surface, mu_o, rho, false, 64) - mean * mean;
            let new = proposal_second_moment(surface, mu_o, rho, true, 64) - mean * mean;
            let coarse = proposal_second_moment(surface, mu_o, rho, true, 32) - mean * mean;
            assert!((new - coarse).abs() < 2e-6);
            println!("cosine={mu_o}, albedo={rho}, mean={mean:.9}, old_variance={old:.9}, mixture_variance={new:.9}, variance_ratio={:.3}, refinement={:.3e}",
                old / new, (new - coarse).abs());
            assert!(new > 0.0 && old > 6.0 * new,
                "cosine={mu_o}, albedo={rho}, old variance={old}, mixture variance={new}");
        }
    }
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

