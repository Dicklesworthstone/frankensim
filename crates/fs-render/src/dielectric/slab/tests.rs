//! G1/G3: independent internal-bounce sums, power closure, source inversion.

use super::*;
use crate::dielectric::{CauchyIor, GlassProvenance};

fn glass(alpha: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(alpha).unwrap(),
        GlassProvenance::Custom,
    )
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn clear_plate_retains_both_interface_reflections_without_invented_heat() {
    let response = glass(0.0)
        .incoherent_slab_response(550.0, 0.01, 1.0, 1.0)
        .unwrap();
    let power = response.unpolarized();
    // n=1.5 gives single-interface R=1/25; the incoherent plate gives 1/13.
    close(power.reflectance, 1.0 / 13.0, 3.0e-15);
    close(power.transmittance, 12.0 / 13.0, 3.0e-15);
    assert_eq!(power.absorptance, 0.0);
    assert_eq!(response.optical_depth, 0.0);
    assert_eq!(response.perpendicular, response.parallel);
}

#[test]
fn oblique_polarized_slab_matches_independent_bounce_sum() {
    let alpha: f64 = 42.0;
    let thickness: f64 = 0.018;
    let cosine: f64 = 0.3;
    let ct = (1.0 - (1.0 - cosine * cosine) / 2.25).sqrt();
    let rs = ((cosine - 1.5 * ct) / (cosine + 1.5 * ct)).powi(2);
    let rp = ((1.5 * cosine - ct) / (1.5 * cosine + ct)).powi(2);
    let response = glass(alpha)
        .incoherent_slab_response(600.0, thickness, cosine, 1.0)
        .unwrap();
    close(response.optical_depth, alpha * thickness / ct, 3.0e-15);
    for (power, r) in [(response.perpendicular, rs), (response.parallel, rp)] {
        let a = (-alpha * thickness / ct).exp();
        let mut entering = 1.0 - r;
        let mut reflected = r;
        let mut transmitted = 0.0;
        let mut absorbed = 0.0;
        for pass in 0..200 {
            absorbed += entering * (1.0 - a);
            let leaving = entering * a;
            if pass % 2 == 0 {
                transmitted += leaving * (1.0 - r);
            } else {
                reflected += leaving * (1.0 - r);
            }
            entering = leaving * r;
        }
        close(power.reflectance, reflected, 8.0e-15);
        close(power.transmittance, transmitted, 8.0e-15);
        close(power.absorptance, absorbed, 8.0e-15);
        close(
            power.reflectance + power.transmittance + power.absorptance,
            1.0,
            4.0e-15,
        );
    }
    // Oblique bounces preserve s/p polarization. Summing with mean R is wrong.
    let r_mean = 0.5 * (rs + rp);
    let a = (-alpha * thickness / ct).exp();
    let premature_average = (1.0 - r_mean).powi(2) * a / (1.0 - (r_mean * a).powi(2));
    assert!((response.unpolarized().transmittance - premature_average).abs() > 1.0e-3);
}

#[test]
fn thickness_and_angle_change_opacity_with_nonnegative_closed_power() {
    for alpha in [0.0, 1.0e-18, 0.02, 20.0, 1.0e5] {
        for cosine in [1.0, 0.8, 0.2, 1.0e-10] {
            let mut last_t = 1.0;
            let mut last_a = 0.0;
            for thickness in [1.0e-6, 0.001, 0.01, 1.0, 100.0] {
                let power = glass(alpha)
                    .incoherent_slab_response(550.0, thickness, cosine, 1.0)
                    .unwrap()
                    .unpolarized();
                assert!(power.transmittance <= last_t);
                assert!(power.absorptance >= last_a);
                for value in [power.reflectance, power.transmittance, power.absorptance] {
                    assert!((0.0..=1.0).contains(&value));
                }
                close(
                    power.reflectance + power.transmittance + power.absorptance,
                    1.0,
                    6.0e-15,
                );
                last_t = power.transmittance;
                last_a = power.absorptance;
            }
        }
    }
}

#[test]
fn thin_absorption_is_retained_separately_from_rounded_transmission() {
    let power = glass(1.0e-20)
        .incoherent_slab_response(550.0, 0.2, 1.0, 1.0)
        .unwrap()
        .unpolarized();
    // dA/dtau at tau=0 is one for a symmetric normal-incidence plate.
    close(power.absorptance, 2.0e-21, 2.0e-15);
    let thick = glass(1.0e5)
        .incoherent_slab_response(550.0, 1.0, 1.0, 1.0)
        .unwrap()
        .unpolarized();
    close(thick.reflectance, 0.04, 3.0e-15);
    close(thick.absorptance, 0.96, 3.0e-15);
    assert_eq!(thick.transmittance, 0.0);
}

#[test]
fn index_matching_removes_interfaces_and_total_reflection_admits_no_heat() {
    let medium = glass(20.0);
    let matched = medium
        .incoherent_slab_response(550.0, 0.01, 0.25, 1.5)
        .unwrap()
        .unpolarized();
    assert_eq!(matched.reflectance, 0.0);
    close(matched.transmittance, (-0.8_f64).exp(), 4.0e-15);
    close(matched.absorptance, -(-0.8_f64).exp_m1(), 4.0e-15);
    let reflected = medium
        .incoherent_slab_response(550.0, 1.0, 0.2, 2.0)
        .unwrap()
        .unpolarized();
    assert_eq!(
        reflected,
        SlabPowerBalance {
            reflectance: 1.0,
            transmittance: 0.0,
            absorptance: 0.0
        }
    );
    assert_eq!(
        medium.incoherent_slab_response(900.0, 1.0, 0.2, 2.0),
        Err(DielectricError::InvalidWavelength)
    );
    let grazing = medium
        .incoherent_slab_response(550.0, 0.01, 1.0e-170, 1.5)
        .unwrap()
        .unpolarized();
    assert_eq!(
        grazing,
        SlabPowerBalance {
            reflectance: 0.0,
            transmittance: 0.0,
            absorptance: 1.0,
        }
    );
}

#[test]
fn internal_transmission_import_uses_spectral_thickness_scaling() {
    let knots = [(8000.0, 1.0), (10_000.0, 0.8), (14_000.0, 0.3)];
    let law = BeerLambertAbsorption::try_from_spectral_internal_transmittance_in_band(&knots, 0.01)
        .unwrap();
    assert_eq!(law.wavelength_bounds_nm(), Some((8000.0, 14_000.0)));
    for (wavelength, transmission) in knots {
        close(
            law.transmittance(wavelength, 0.02).unwrap(),
            transmission.powi(2),
            4.0e-15,
        );
    }
    // Linear alpha interpolation corresponds to geometric transmission.
    close(
        law.transmittance(9000.0, 0.01).unwrap(),
        0.8_f64.sqrt(),
        4.0e-15,
    );
    assert_eq!(
        law.transmittance(7999.0, 0.01),
        Err(DielectricError::InvalidWavelength)
    );
}

#[test]
fn normal_slab_measurements_recover_native_bulk_law_for_new_thicknesses() {
    let thickness = 0.008;
    for exterior in [1.0, 1.33, 1.5, 2.0] {
        for alpha in [0.0, 1.0e-6, 0.1, 200.0, 20_000.0] {
            let original = glass(alpha);
            let measured = original
                .incoherent_slab_response(550.0, thickness, 1.0, exterior)
                .unwrap()
                .unpolarized()
                .transmittance;
            let imported = BeerLambertAbsorption::try_from_normal_slab_transmittance_in_band(
                &[(380.0, measured), (780.0, measured)],
                thickness,
                original.ior(),
                exterior,
            )
            .unwrap();
            if alpha == 0.0 {
                assert_eq!(imported.extinction_per_m(550.0), Ok(0.0));
            } else {
                close(
                    imported.extinction_per_m(550.0).unwrap(),
                    alpha,
                    if alpha < 1.0e-5 { 8.0e-8 } else { 8.0e-13 },
                );
            }
            let reconstructed =
                DielectricGlass::with_ior(*original.ior(), imported, GlassProvenance::Custom);
            for length in [0.002, 0.012] {
                let expected = original
                    .incoherent_slab_response(550.0, length, 0.8, exterior)
                    .unwrap()
                    .unpolarized();
                let actual = reconstructed
                    .incoherent_slab_response(550.0, length, 0.8, exterior)
                    .unwrap()
                    .unpolarized();
                close(actual.transmittance, expected.transmittance, 2.0e-12);
                if alpha > 0.0 {
                    close(actual.absorptance, expected.absorptance, 8.0e-8);
                }
            }
        }
    }
}

#[test]
fn external_clear_loss_is_not_misinterpreted_as_bulk_absorption() {
    let original = glass(0.0);
    let clear = original
        .incoherent_slab_response(550.0, 0.01, 1.0, 1.0)
        .unwrap()
        .unpolarized()
        .transmittance;
    let imported = BeerLambertAbsorption::try_from_normal_slab_transmittance_in_band(
        &[(380.0, clear), (780.0, clear)],
        0.01,
        original.ior(),
        1.0,
    )
    .unwrap();
    assert_eq!(imported.extinction_per_m(550.0), Ok(0.0));
    assert!(
        -clear.ln() / 0.01 > 8.0,
        "naive external Beer inversion invents strong absorption"
    );
    let slightly_lossy = f64::from_bits(clear.to_bits() - 1);
    let weak = BeerLambertAbsorption::try_from_normal_slab_transmittance_in_band(
        &[(380.0, slightly_lossy), (780.0, slightly_lossy)],
        0.01,
        original.ior(),
        1.0,
    )
    .unwrap();
    assert!(weak.extinction_per_m(550.0).unwrap() > 0.0);
}

#[test]
fn invalid_measurements_and_geometries_refuse_without_clipping() {
    let original = glass(1.0);
    for value in [0.0, -1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            BeerLambertAbsorption::try_from_spectral_internal_transmittance_in_band(
                &[(380.0, value), (780.0, 0.5)],
                0.01,
            )
            .is_err()
        );
        assert!(
            BeerLambertAbsorption::try_from_normal_slab_transmittance_in_band(
                &[(380.0, value), (780.0, 0.5)],
                0.01,
                original.ior(),
                1.0,
            )
            .is_err()
        );
    }
    assert!(
        BeerLambertAbsorption::try_from_normal_slab_transmittance_in_band(
            &[(380.0, 0.99), (780.0, 0.99)],
            0.01,
            original.ior(),
            1.0,
        )
        .is_err()
    );
    for thickness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            original
                .incoherent_slab_response(550.0, thickness, 1.0, 1.0)
                .is_err()
        );
    }
    for cosine in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            original
                .incoherent_slab_response(550.0, 0.01, cosine, 1.0)
                .is_err()
        );
    }
}
