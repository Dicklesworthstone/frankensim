//! FetchAdjustedMassConsistent battery (bead `frankensim-wf-root-guzez.4.7`, E3.3c).

use fs_atmo::fetch_mass_consistent::{
    ExponentialDensityProfile, FetchAdjustedAnelastic, FetchAdjustedMassConsistent,
    FetchRoughnessProfile, MODEL_ID_FETCH_MASS_CONSISTENT,
};

#[test]
fn fetch_mass_consistent_divergence_cancels_to_machine_precision() {
    let profile = FetchRoughnessProfile::new(0.01, 0.0005, 500.0).expect("valid profile");
    let law = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 15.0).expect("valid law");

    assert_eq!(
        MODEL_ID_FETCH_MASS_CONSISTENT,
        "FetchAdjustedMassConsistent"
    );

    // Sample across grid of (x, h) points
    for ix in 0..20 {
        let x = ix as f64 * 25.0;
        for ih in 1..20 {
            let h = ih as f64 * 1.5;
            let div = law.divergence(x, h);
            assert!(
                div.abs() < 1e-12,
                "divergence at x={x}, h={h} must cancel to machine precision, got {div}"
            );
        }
    }
}

#[test]
fn fetch_mass_consistent_wall_impermeability() {
    let profile = FetchRoughnessProfile::new(0.02, 0.0002, 300.0).expect("valid profile");
    let law = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 12.0).expect("valid law");

    for ix in 0..10 {
        let x = ix as f64 * 30.0;
        let v_ground = law.sample_velocity(x, 0.0);
        assert_eq!(
            v_ground[2], 0.0,
            "vertical velocity at ground must be exactly zero"
        );
    }
}

#[test]
fn fetch_mass_consistent_refusals_and_boundaries() {
    // Non-finite roughness refuses
    assert!(FetchRoughnessProfile::new(f64::NAN, 0.001, 100.0).is_err());
    // z0 below minimum domain refuses
    assert!(FetchRoughnessProfile::new(1e-7, 0.001, 100.0).is_err());
    // z0 above maximum domain refuses
    assert!(FetchRoughnessProfile::new(2.5, 0.001, 100.0).is_err());
    // Non-positive max fetch refuses
    assert!(FetchRoughnessProfile::new(0.01, 0.001, 0.0).is_err());

    // Reference height below d + z0 refuses
    let profile = FetchRoughnessProfile::new(0.05, 0.0001, 100.0).expect("valid profile");
    assert!(FetchAdjustedMassConsistent::new(profile, 2.0, 2.02, 10.0).is_err());
}

#[test]
fn anelastic_mass_flux_conserves_mass_for_both_roughness_gradients() {
    let inverse_scale = 0.002;
    let density = ExponentialDensityProfile::new(1.2, 10.0, inverse_scale, [0.0, 1000.0]).unwrap();
    for (z0, gradient) in [(0.02, 0.0004), (0.25, -0.0004)] {
        let roughness = FetchRoughnessProfile::new(z0, gradient, 500.0).unwrap();
        let wind = FetchAdjustedMassConsistent::new(roughness, 1.5, 10.0, 15.0).unwrap();
        let law = FetchAdjustedAnelastic::new(wind, density).unwrap();
        for x in [40.0, 160.0, 400.0] {
            for h in [5.0, 40.0, 200.0, 600.0] {
                // Independently differentiate the published physical flux, not an
                // analytic divergence method that repeats the construction.
                let step = 0.01;
                let xp = law.sample(x + step, h).unwrap();
                let xm = law.sample(x - step, h).unwrap();
                let hp = law.sample(x, h + step).unwrap();
                let hm = law.sample(x, h - step).unwrap();
                let dx = (xp.mass_flux_kg_m2_s[0] - xm.mass_flux_kg_m2_s[0]) / (2.0 * step);
                let dh = (hp.mass_flux_kg_m2_s[2] - hm.mass_flux_kg_m2_s[2]) / (2.0 * step);
                assert!(
                    (dx + dh).abs() < 1e-7 * (dx.abs() + dh.abs()),
                    "x={x}, h={h}"
                );

                let du_dx = (xp.velocity_mps[0] - xm.velocity_mps[0]) / (2.0 * step);
                let dw_dh = (hp.velocity_mps[2] - hm.velocity_mps[2]) / (2.0 * step);
                let sample = law.sample(x, h).unwrap();
                let div_u = du_dx + dw_dh;
                assert!(
                    div_u.abs() > 1e-6,
                    "stratified mass continuity is not solenoidality"
                );
                assert!(
                    (div_u - inverse_scale * sample.velocity_mps[2]).abs()
                        < 1e-7 * (du_dx.abs() + dw_dh.abs())
                );
                let expected_density = 1.2 * (-inverse_scale * (h - 10.0)).exp();
                assert!((sample.density_kg_m3 / expected_density - 1.0).abs() < 2e-14);
            }
        }
    }
}

#[test]
fn anelastic_constant_density_is_exact_and_weak_stratification_is_resolved() {
    let profile = FetchRoughnessProfile::new(0.02, 0.0004, 500.0).unwrap();
    let wind = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 15.0).unwrap();
    for inverse_scale in [0.0, 1e-14, f64::from_bits(1)] {
        let density =
            ExponentialDensityProfile::new(1.2, 10.0, inverse_scale, [0.0, 1000.0]).unwrap();
        let law = FetchAdjustedAnelastic::new(wind, density).unwrap();
        for x in [-5.0, 0.0, 20.0, 250.0, 500.0, 510.0] {
            for h in [0.0, 0.05, 4.0, 100.0, 600.0] {
                let sample = law.sample(x, h).unwrap();
                let original = wind.sample_velocity(x, h);
                if inverse_scale == 0.0 {
                    assert_eq!(sample.velocity_mps, original);
                    assert_eq!(sample.density_kg_m3, 1.2);
                    assert_eq!(sample.mass_flux_kg_m2_s, original.map(|v| 1.2 * v));
                } else {
                    for (actual, expected) in sample.velocity_mps.into_iter().zip(original) {
                        assert!((actual - expected).abs() <= 1e-11 * expected.abs());
                    }
                }
            }
        }
    }
}

#[test]
fn anelastic_preserves_cutoff_and_clamped_fetch_behavior() {
    let profile = FetchRoughnessProfile::new(0.125, 1.0 / 2048.0, 512.0).unwrap();
    let wind = FetchAdjustedMassConsistent::new(profile, 2.0, 10.0, 12.0).unwrap();
    let density = ExponentialDensityProfile::new(1.2, 10.0, 0.002, [0.0, 1000.0]).unwrap();
    let law = FetchAdjustedAnelastic::new(wind, density).unwrap();
    for x in [-10.0, 0.0, 128.0, 512.0, 520.0] {
        let cutoff = 2.0 + profile.z0_at(x);
        for h in [0.0, cutoff - 0.125, cutoff] {
            let sample = law.sample(x, h).unwrap();
            assert_eq!(sample.velocity_mps, [0.0; 3]);
            assert_eq!(sample.mass_flux_kg_m2_s, [0.0; 3]);
            assert!(sample.density_kg_m3 > 0.0);
        }
    }
    for (outside, endpoint) in [(-10.0, 0.0), (520.0, 512.0)] {
        let sample = law.sample(outside, 100.0).unwrap();
        let edge = law.sample(endpoint, 100.0).unwrap();
        assert_eq!(sample.velocity_mps[0], edge.velocity_mps[0]);
        assert_eq!(sample.velocity_mps[2], 0.0);
        assert_eq!(sample.mass_flux_kg_m2_s[2], 0.0);
        assert!(edge.velocity_mps[2] > 0.0);
    }
}

#[test]
fn anelastic_density_scale_changes_mass_flux_without_changing_velocity() {
    let profile = FetchRoughnessProfile::new(0.02, 0.0004, 500.0).unwrap();
    let wind = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 15.0).unwrap();
    let low = FetchAdjustedAnelastic::new(
        wind,
        ExponentialDensityProfile::new(1.2, 10.0, 0.002, [0.0, 1000.0]).unwrap(),
    )
    .unwrap();
    let high = FetchAdjustedAnelastic::new(
        wind,
        ExponentialDensityProfile::new(12.0, 10.0, 0.002, [0.0, 1000.0]).unwrap(),
    )
    .unwrap();
    for h in [0.0, 10.0, 200.0, 600.0] {
        let a = low.sample(160.0, h).unwrap();
        let b = high.sample(160.0, h).unwrap();
        assert_eq!(a.velocity_mps, b.velocity_mps);
        assert!((b.density_kg_m3 / a.density_kg_m3 - 10.0).abs() < 1e-14);
        for (original, scaled) in a.mass_flux_kg_m2_s.into_iter().zip(b.mass_flux_kg_m2_s) {
            assert!((scaled - 10.0 * original).abs() <= 1e-13 * scaled.abs());
        }
    }
}

#[test]
fn anelastic_refuses_invalid_sources_queries_and_overflow() {
    assert!(ExponentialDensityProfile::new(-1.0, 10.0, 0.0, [0.0, 100.0]).is_err());
    assert!(ExponentialDensityProfile::new(f64::NAN, 10.0, 0.0, [0.0, 100.0]).is_err());
    assert!(ExponentialDensityProfile::new(1.2, 10.0, -0.1, [0.0, 100.0]).is_err());
    assert!(ExponentialDensityProfile::new(1.2, 10.0, 0.0, [100.0, 0.0]).is_err());
    assert!(ExponentialDensityProfile::new(1.2, 10.0, 0.0, [0.0, f64::INFINITY]).is_err());
    assert!(ExponentialDensityProfile::new(1.2, 10.0, 0.0, [20.0, 100.0]).is_err());
    assert!(ExponentialDensityProfile::new(1.2, 0.0, 1.0, [0.0, 1000.0]).is_err());
    assert!(ExponentialDensityProfile::new(1.2, 1000.0, 1.0, [0.0, 1000.0]).is_err());

    let profile = FetchRoughnessProfile::new(0.02, 0.0004, 500.0).unwrap();
    let wind = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 15.0).unwrap();
    let density = ExponentialDensityProfile::new(1.2, 10.0, 0.002, [0.0, 1000.0]).unwrap();
    let law = FetchAdjustedAnelastic::new(wind, density).unwrap();
    for (x, h) in [
        (f64::NAN, 10.0),
        (10.0, f64::INFINITY),
        (10.0, -1.0),
        (10.0, 1001.0),
    ] {
        assert!(law.sample(x, h).is_err());
    }
    let mut invalid_wind = wind;
    invalid_wind.fetch_profile.dz0_dx = f64::NAN;
    assert!(FetchAdjustedAnelastic::new(invalid_wind, density).is_err());
    let excluded_cutoff = ExponentialDensityProfile::new(1.2, 10.0, 0.0, [1.0, 100.0]).unwrap();
    assert!(FetchAdjustedAnelastic::new(wind, excluded_cutoff).is_err());

    let dense = ExponentialDensityProfile::new(f64::MAX, 10.0, 0.0, [0.0, 100.0]).unwrap();
    assert!(
        FetchAdjustedAnelastic::new(wind, dense)
            .unwrap()
            .sample(0.0, 10.0)
            .is_err()
    );
    // Both endpoint densities are representable, but the necessary vertical
    // velocity at the top is not. Clamped fetch requires no such correction.
    let contrast = ExponentialDensityProfile::new(1.0, 500.0, 1.0, [0.0, 1000.0]).unwrap();
    let steep = FetchAdjustedAnelastic::new(wind, contrast).unwrap();
    assert!(steep.sample(50.0, 1000.0).is_err());
    assert!(steep.sample(-1.0, 1000.0).is_ok());
}
