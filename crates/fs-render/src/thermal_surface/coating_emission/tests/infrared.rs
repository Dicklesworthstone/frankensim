use super::*;
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::thermal_surface::{
    CoatingExchangeBudget, CoatingExchangeQuadrature, CoatingExchangeResponse,
};

fn infrared_film(bounds: (f64, f64), thickness: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                1.5,
                0.3,
            )
            .unwrap()
        }),
        thickness,
        (bounds.1 - bounds.0) / 8.0,
        source(20),
    )
    .unwrap()
}

fn infrared_substrate(bounds: (f64, f64), eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new_in_band(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                eta,
                k,
            )
            .unwrap()
        }),
        (bounds.1 - bounds.0) / 8.0,
        source(21),
    )
    .unwrap()
}

fn infrared_vacuum(gap: f64) -> DielectricIor {
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(8000.0 + 1000.0 * i as f64, 1.0).unwrap()
        }),
        gap,
    )
    .unwrap()
    .into()
}

#[test]
fn infrared_nonisothermal_emission_matches_independent_depth_and_si_planck_sources() {
    let films = [
        infrared_film((4000.0, 20_000.0), 2e-6),
        infrared_film((5000.0, 18_000.0), 3e-6),
    ];
    let temperatures = [280.0, 310.0];
    let emitter = CoatingThermalEmission::try_opaque(
        &[
            CoatingLayer::Absorbing(&films[0]),
            CoatingLayer::Absorbing(&films[1]),
        ],
        &temperatures,
        infrared_substrate((6000.0, 24_000.0), 1.5, 0.3),
        340.0,
        ContentHash([22; 32]),
    )
    .unwrap();
    assert_eq!(emitter.wavelength_bounds_nm(), (6000.0, 18_000.0));
    for lambda in [6000.0, 10_000.0, 18_000.0] {
        let alpha = 4.0 * PI * 0.3 / (lambda * 1e-9);
        let entrance = 4.0 * 1.5 / (2.5 * 2.5 + 0.3 * 0.3);
        let sample = emitter.evaluate(lambda, 1.0).unwrap();
        let absorbed = [
            entrance * -(-alpha * 2e-6).exp_m1(),
            entrance * (-alpha * 2e-6).exp() * -(-alpha * 3e-6).exp_m1(),
        ];
        for index in 0..2 {
            let expected = 0.5 * absorbed[index] * planck_si(temperatures[index], lambda);
            for value in sample.film_radiance_w_m2_sr_nm()[index] {
                close(value, expected);
            }
        }
        let expected = 0.5 * entrance * (-alpha * 5e-6).exp() * planck_si(340.0, lambda);
        for value in sample.substrate_radiance_w_m2_sr_nm().unwrap() {
            close(value, expected);
        }
        assert!(sample.stokes().i > 0.0);
    }
    for lambda in [5999.0, 18_001.0, 550.0] {
        assert!(emitter.evaluate(lambda, 1.0).is_err());
    }
    assert_eq!(
        CoatingThermalEmission::try_opaque(
            &[CoatingLayer::Absorbing(&films[0])],
            &[300.0],
            substrate(1.5, 0.3),
            300.0,
            ContentHash([22; 32]),
        )
        .unwrap_err(),
        CoatingEmissionError::WavelengthDomain
    );
}

#[test]
fn opposite_medium_and_zero_thickness_laws_bound_infrared_modes_and_identity() {
    let film = infrared_film((4000.0, 20_000.0), 1e-6);
    let make = |ior| {
        CoatingThermalEmission::try_transmitting(
            &[CoatingLayer::Absorbing(&film)],
            &[300.0],
            ior,
            true,
            ContentHash([23; 32]),
        )
    };
    let emitter = make(infrared_vacuum(1000.0)).unwrap();
    assert_eq!(emitter.wavelength_bounds_nm(), (8000.0, 16_000.0));
    let sample = emitter.evaluate(10_000.0, 0.6).unwrap();
    for pol in 0..2 {
        close(
            sample.film_radiance_w_m2_sr_nm()[0][pol],
            sample.reciprocal_absorbed_fraction()[pol] * 0.5 * planck_si(300.0, 10_000.0),
        );
    }
    assert!(emitter.evaluate(7000.0, 0.6).is_err());
    assert!(
        make(vacuum()).is_err(),
        "a visible constant-index law does not authorize IR"
    );
    assert_ne!(
        emitter.optical_identity(),
        make(infrared_vacuum(2000.0)).unwrap().optical_identity()
    );
    let dark = CoatingThermalEmission::try_transmitting(
        &[],
        &[],
        infrared_vacuum(1000.0),
        false,
        ContentHash([23; 32]),
    )
    .unwrap();
    assert_eq!(dark.wavelength_bounds_nm(), (8000.0, 16_000.0));
    assert_eq!(dark.evaluate(10_000.0, 0.6).unwrap().stokes().i, 0.0);
    assert!(dark.evaluate(550.0, 0.6).is_err());
    let zero_visible = super::film(1.5, 0.3, 0.0, 24);
    assert!(
        CoatingThermalEmission::try_transmitting(
            &[CoatingLayer::Absorbing(&zero_visible)],
            &[300.0],
            infrared_vacuum(1000.0),
            false,
            ContentHash([23; 32]),
        )
        .is_err(),
        "a zero-thickness material cannot silently lose its source domain"
    );
}

#[test]
fn infrared_exchange_converges_to_independent_finite_band_planck_integral() {
    let bounds = (2000.0, 80_000.0);
    let emitter = CoatingThermalEmission::try_opaque(
        &[],
        &[],
        infrared_substrate(bounds, 1.0, 1e-20),
        300.0,
        ContentHash([25; 32]),
    )
    .unwrap();
    let area = 0.01;
    let reference = {
        let n = 16_384;
        let width = (bounds.1 - bounds.0) / n as f64;
        let value = |lambda| PI * area * (planck_si(300.0, lambda) - planck_si(250.0, lambda));
        let mut sum = value(bounds.0) + value(bounds.1);
        for i in 1..n {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * value(bounds.0 + i as f64 * width);
        }
        sum * width / 3.0
    };
    let evaluate = |bins| {
        let response = CoatingExchangeResponse::from_emitter(
            &emitter,
            area,
            CoatingExchangeQuadrature {
                wavelength_bins: bins,
                cosine_bins: 1,
            },
            CoatingExchangeBudget {
                max_spectral_nodes: bins,
                max_optical_evaluations: bins,
            },
            || Ok(()),
        )
        .unwrap();
        assert_eq!(response.wavelength_bounds_nm(), bounds);
        let rate = response.evaluate(300.0, Some(250.0), || Ok(())).unwrap();
        close(rate.secant_w_k, rate.outward_power_w / 50.0);
        rate.outward_power_w
    };
    let coarse = evaluate(64);
    let fine = evaluate(1024);
    assert!(
        fine > 1.0,
        "room-temperature IR exchange must not collapse to visible Wien-tail power"
    );
    assert!((fine - reference).abs() < (coarse - reference).abs() / 100.0);
    assert!((fine - reference).abs() < 2e-6 * reference);
}
