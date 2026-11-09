//! Original current alpha(T) combined with immutable reference-only phase.
use super::*;
use fs_material::state_point::{
    SpectralDielectricAbsorptionStatePoint, resolve_spectral_dielectric_absorption_state_point,
};

fn source_card(slope: f64, overflowing_reference: bool) -> MaterialCard {
    let cold_alpha = core::array::from_fn(|i| {
        if overflowing_reference { f64::MAX } else if i == 4 { 2000.0 } else { 1000.0 }
    });
    let base = glass_card_spectrum(cold_alpha, 0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut value = original.clone();
        let name = value.key.name();
        if VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES.contains(&name)
            || VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES.contains(&name)
            || SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES.contains(&name)
            || SPECTRAL_DIELECTRIC_IOR_PROPERTIES.contains(&name)
        {
            value.validity = ValidityDomain::unconstrained().with("T", 300.0, 300.0);
        } else if let Some(i) = SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES
            .iter().position(|&n| n == name)
        {
            // Phase fit: 800..2000 nm. Reference alpha: 1000..1800 nm.
            // Current alpha: 900..1900 nm. Neither alpha band owns phase support.
            value.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, (1000.0 + 100.0 * i as f64) * 1e-9),
                    (400.0, (900.0 + 125.0 * i as f64) * 1e-9),
                ],
                dims: value.key.dims(),
            };
            value.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
            value.interpolation = InterpolationPolicy::LinearInside;
        } else if let Some(i) = SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES
            .iter().position(|&n| n == name)
        {
            value.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, cold_alpha[i]),
                    (400.0, 100_000.0 + slope * ((900.0 + 125.0 * i as f64) * 1e-9 - 1e-6)),
                ],
                dims: value.key.dims(),
            };
            value.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
            value.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(value).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn current(card: &MaterialCard, t: f64) -> SpectralDielectricAbsorptionStatePoint {
    let alpha = resolve_spectral_dielectric_absorption_state_point(
        card, &point(t), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let complete = alpha.resolved().with_additional_requirements(
        card,
        &[ScalarPropertyRequirement::try_new(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    SpectralDielectricAbsorptionStatePoint::try_from_material(complete).unwrap()
}

#[test]
fn current_loss_has_independent_support_and_exact_complex_frequency_derivatives() {
    let lambda_nm = 1535.0;
    let lambda_m = lambda_nm * 1e-9;
    let omega = FREQUENCY_WAVELENGTH_NM / lambda_nm;
    for slope in [0.0, 2e10] {
        let card = source_card(slope, false);
        let original = SpectralSellmeierStatePoint::try_from_material(
            resolve(&card, &point(300.0)),
        ).unwrap();
        assert!(original.at_query_point(&card, &point(400.0)).is_err());
        let absorption = current(&card, 400.0);
        assert_eq!(absorption.resolved().properties().len(), 19);
        let law = law();
        let response = ThermoOpticOpticalResponse::try_with_current_absorption(
            reference(&original), &law, &absorption, "T", 160e-9,
        ).unwrap();
        assert!(!response.retains_reference_absorption());
        assert!(core::ptr::eq(response.current_absorption().unwrap(), &absorption));
        assert!(core::ptr::eq(response.material(), absorption.resolved()));
        let alpha_at = |lambda: f64| 100_000.0 + slope * (lambda - 1e-6);
        let index = response.complex_index(omega).unwrap();
        near(index.re, thermal_phase(lambda_nm, 400.0), 2e-15);
        near(index.im, -alpha_at(lambda_m) * lambda_m / (4.0 * core::f64::consts::PI), 2e-17);
        let bulk = DielectricResponse::ThermoOpticOptical(response)
            .bulk_dispersion(omega).unwrap();
        near(-2.0 * bulk.wavenumber_rad_m().im, alpha_at(lambda_m), 5e-10);
        let group = slope * lambda_m.powi(2) / (4.0 * core::f64::consts::PI);
        near(bulk.wavenumber_derivative_s_m().im * VACUUM_SPEED_M_S, group, 3e-16);
        near(
            bulk.wavenumber_second_derivative_s2_m().im * omega * VACUUM_SPEED_M_S,
            -2.0 * group, 6e-16,
        );
        if slope == 0.0 {
            assert_eq!(bulk.wavenumber_derivative_s_m().im, 0.0);
            assert_eq!(bulk.wavenumber_second_derivative_s2_m().im, 0.0);
        }
        let epsilon = |w: f64| {
            let lambda = FREQUENCY_WAVELENGTH_NM / w;
            let n = C64::new(
                thermal_phase(lambda, 400.0),
                -alpha_at(lambda * 1e-9) * lambda * 1e-9 / (4.0 * core::f64::consts::PI),
            );
            n * n
        };
        let h = 1e-3;
        let f = [-2.0, -1.0, 0.0, 1.0, 2.0].map(|i| epsilon(omega * (1.0 + i * h)));
        let first = (f[0] - f[1].scale(8.0) + f[3].scale(8.0) - f[4])
            .scale(1.0 / (12.0 * h));
        let second = (-f[4] + f[3].scale(16.0) - f[2].scale(30.0)
            + f[1].scale(16.0) - f[0]).scale(1.0 / (12.0 * h * h));
        let jet = response.frequency_jet(omega).unwrap();
        let scale = response.permeability_h_m() * VACUUM_SPEED_M_S.powi(2);
        complex_near(jet.epsilon_derivatives[0].scale(scale), first, 2e-9);
        complex_near(jet.epsilon_derivatives[1].scale(scale), second, 2e-8);
        for edge in response.frequency_range_rad_s() {
            assert!(response.complex_index(edge).is_ok());
            assert!(response.frequency_jet(edge).is_err());
        }
        // These wavelengths lie outside discarded reference alpha support.
        for nm in [950.0, 1850.0] {
            let w = FREQUENCY_WAVELENGTH_NM / nm;
            assert!(reference(&original).complex_index(w).is_err());
            near(response.complex_index(w).unwrap().re, thermal_phase(nm, 400.0), 3e-15);
            assert!(response.frequency_jet(w).is_ok());
        }
        assert!(response.complex_index(FREQUENCY_WAVELENGTH_NM / 850.0).is_err());
        assert!(response.complex_index(FREQUENCY_WAVELENGTH_NM / 1950.0).is_err());
        assert!(response.at_state(absorption.resolved(), "T").is_err());
        assert!(ThermoOpticOpticalResponse::try_with_current_absorption(
            reference(&original), &law, &absorption, "T", 120e-9,
        ).is_err());
        assert!(absorption.at_query_point(&card, &point(401.0)).is_err());
    }
}

#[test]
fn selected_phase_and_current_loss_own_cusps_without_discarded_reference_loss() {
    let card = source_card(0.0, false);
    let original = SpectralSellmeierStatePoint::try_from_material(
        resolve(&card, &point(300.0)),
    ).unwrap();
    let absorption = current(&card, 400.0);
    let law = law();
    let response = ThermoOpticOpticalResponse::try_with_current_absorption(
        reference(&original), &law, &absorption, "T", 160e-9,
    ).unwrap();
    let cusp = core::f64::consts::TAU * VACUUM_SPEED_M_S
        / original.absorption_samples()[4].wavelength_m;
    assert!(reference(&original).frequency_jet(cusp).is_err());
    assert!(response.frequency_jet(cusp).is_ok());
    let cooler = absorption.at_query_point(&card, &point(350.0)).unwrap();
    let refreshed = response.at_absorption_state(&cooler, "T").unwrap();
    let current_cusp = core::f64::consts::TAU * VACUUM_SPEED_M_S
        / cooler.absorption_samples()[4].wavelength_m;
    assert!(refreshed.complex_index(current_cusp).is_ok());
    assert!(refreshed.refractive_index_jet(current_cusp).is_ok());
    assert!(refreshed.frequency_jet(current_cusp).is_err());
    let tabulated = SpectralTabulatedDielectricStatePoint::try_from_material(
        original.resolved().clone(),
    ).unwrap();
    let tabulated_reference = SpectralOpticalResponse::try_new(
        SpectralOpticalSource::TabulatedDielectric(&tabulated), 160e-9,
    ).unwrap();
    let tabulated_thermal = ThermoOpticOpticalResponse::try_with_current_absorption(
        tabulated_reference, &law, &absorption, "T", 160e-9,
    ).unwrap();
    assert!(tabulated_thermal.complex_index(cusp).is_ok());
    assert!(tabulated_thermal.frequency_jet(cusp).is_err());
    assert!(!response.same_source(refreshed));
    let other = current(&source_card(1e10, false), 400.0);
    assert!(response.at_absorption_state(&other, "T").is_err());

    let without_mu = resolve_spectral_dielectric_absorption_state_point(
        &card, &point(400.0), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    assert!(ThermoOpticOpticalResponse::try_with_current_absorption(
        reference(&original), &law, &without_mu, "T", 160e-9,
    ).is_err());

    let huge_card = source_card(0.0, true);
    let huge_reference = SpectralSellmeierStatePoint::try_from_material(
        resolve(&huge_card, &point(300.0)),
    ).unwrap();
    let finite_current = current(&huge_card, 400.0);
    let composed = ThermoOpticOpticalResponse::try_with_current_absorption(
        reference(&huge_reference), &law, &finite_current, "T", 160e-9,
    ).unwrap();
    assert!(reference(&huge_reference).permittivity_f_m(cusp).is_err());
    assert!(composed.frequency_jet(cusp).is_ok());
}

#[test]
fn current_absorption_changes_native_refraction_fields_and_local_heat_once() {
    let card = source_card(0.0, false);
    let original = SpectralSellmeierStatePoint::try_from_material(
        resolve(&card, &point(300.0)),
    ).unwrap();
    let absorption = current(&card, 400.0);
    let law = law();
    let response = ThermoOpticOpticalResponse::try_with_current_absorption(
        reference(&original), &law, &absorption, "T", 160e-9,
    ).unwrap();
    let ambient = port();
    let omega = FREQUENCY_WAVELENGTH_NM / 1535.0;
    let n = C64::new(
        thermal_phase(1535.0, 400.0), -100_000.0 * 1535e-9 / (4.0 * core::f64::consts::PI),
    );
    let k0 = omega / VACUUM_SPEED_M_S;
    let cosine = 0.8_f64;
    let q = (n * n - C64::from_re(1.0 - cosine * cosine)).sqrt();
    let thickness = 2e-6;
    let angle = -k0 * q.re * thickness;
    let p = C64::new(angle.cos(), angle.sin()).scale((k0 * q.im * thickness).exp());
    let films = [DielectricLayer {
        response: DielectricResponse::ThermoOpticOptical(response), thickness_m: thickness,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &ambient, DielectricResponse::Lossless(&ambient), omega, &films,
    ).unwrap();
    let frozen = ThermoOpticOpticalResponse::try_with_reference_absorption(
        reference(&original), &law, absorption.resolved(), "T",
    ).unwrap();
    let frozen_films = [DielectricLayer {
        response: DielectricResponse::ThermoOpticOptical(frozen), thickness_m: thickness,
    }];
    let frozen_stack = PlanarDielectricLayers::try_new(
        &ambient, DielectricResponse::Lossless(&ambient), omega, &frozen_films,
    ).unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let yi = match polarization {
            Polarization::Te => cosine, Polarization::Tm => cosine.recip(),
        } / (MU * VACUUM_SPEED_M_S);
        let yt = match polarization {
            Polarization::Te => q, Polarization::Tm => n * n / q,
        }.scale(1.0 / (response.permeability_h_m() * VACUUM_SPEED_M_S));
        let r = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
        let denominator = C64::ONE - r * r * p * p;
        let reflected = r * (C64::ONE - p * p) / denominator;
        let transmitted = p * (C64::ONE - r * r) / denominator;
        let field = stack.prepare_layer_field(0, cosine, polarization).unwrap();
        let expected = 0.5 * yi * (1.0 - reflected.norm_sq() - transmitted.norm_sq());
        near(field.absorbed_power_w_m2(), expected, 5e-17);
        let observed = stack.response(cosine, polarization).unwrap();
        complex_near(observed.reflected_tangential_e, reflected, 3e-14);
        complex_near(observed.transmitted_tangential_e, transmitted, 3e-14);
        assert!(field.absorbed_power_w_m2() > 20.0 * frozen_stack
            .prepare_layer_field(0, cosine, polarization).unwrap().absorbed_power_w_m2());
    }
}
