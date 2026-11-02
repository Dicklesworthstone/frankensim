use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use crate::dielectric::CauchyIor;
use core::f64::consts::PI;

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}
fn metal(bounds: (f64, f64), eta: f64, k: f64, tag: u8) -> ConductorOptics {
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
        source(tag),
    )
    .unwrap()
}
fn film(eta: f64, k: f64, thickness: f64, tag: u8) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                4000.0 + 2000.0 * i as f64,
                eta + 0.02 * i as f64,
                k + 0.005 * i as f64,
            )
            .unwrap()
        }),
        thickness,
        2000.0,
        source(tag),
    )
    .unwrap()
}
fn emitter(metal: ConductorOptics, temperature: f64, tag: u8) -> CoatingThermalEmission {
    CoatingThermalEmission::try_opaque(&[], &[], metal, temperature, ContentHash([tag; 32]))
        .unwrap()
}
fn budget(q: CoatingExchangeQuadrature) -> CoatingExchangeBudget {
    CoatingExchangeBudget {
        max_spectral_nodes: q.wavelength_bins,
        max_optical_evaluations: 2 * q.wavelength_bins * q.cosine_bins,
    }
}
fn response(
    a: &CoatingThermalEmission,
    b: &CoatingThermalEmission,
    area: f64,
    n: usize,
    m: usize,
) -> CoatingGapResponse {
    let q = CoatingExchangeQuadrature {
        wavelength_bins: n,
        cosine_bins: m,
    };
    CoatingGapResponse::from_emitters(a, b, area, q, budget(q), || Ok(())).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + 3.0 * f64::from_bits(1),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn planck(t: f64, wavelength_nm: f64) -> f64 {
    let lambda = wavelength_nm * 1e-9;
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    2.0 * h * c * c / (lambda.powi(5) * (h * c / (lambda * k * t)).exp_m1()) * 1e-9
}
fn simpson(bounds: (f64, f64), n: usize, mut value: impl FnMut(f64) -> f64) -> f64 {
    let step = (bounds.1 - bounds.0) / n as f64;
    let mut sum = value(bounds.0) + value(bounds.1);
    for i in 1..n {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * value(bounds.0 + i as f64 * step);
    }
    sum * step / 3.0
}
// Independent intensity ray history: each emitted packet may return to the
// other plate repeatedly; accumulate only work actually absorbed there.
fn reflected_ray_work(first: f64, second: f64) -> f64 {
    let mut departing = first;
    let mut received = 0.0;
    for _ in 0..8192 {
        received += departing * second;
        departing *= 1.0 - second;
        departing *= 1.0 - first;
    }
    received
}

#[test]
fn g1_gray_multiple_reflections_match_independent_ray_work_and_keep_weak_exchange() {
    for (first, second) in [(1.0, 1.0), (1.0, 0.2), (0.2, 0.4), (0.01, 0.03)] {
        let transfer = multiple_reflection_exchange(first, second).unwrap();
        close(transfer, reflected_ray_work(first, second), 2e-13);
        assert_eq!(
            transfer,
            multiple_reflection_exchange(second, first).unwrap()
        );
        assert!(transfer <= first.min(second));
    }
    assert_eq!(multiple_reflection_exchange(0.0, 0.0).unwrap(), 0.0);
    assert_eq!(multiple_reflection_exchange(0.0, 0.8).unwrap(), 0.0);
    let tiny = 1e-200;
    assert_eq!(tiny * tiny, 0.0);
    close(
        multiple_reflection_exchange(tiny, tiny).unwrap(),
        0.5 * tiny,
        2e-15,
    );

    // Both original absorptances lie below a rounded 1-R deficit. Even their
    // product with this area vanishes, while watts AND W/K remain representable.
    let body = emitter(metal((380.0, 780.0), 1e-120, 100.0, 1), 300.0, 2);
    let pair = response(&body, &body, 1e-202, 1, 1);
    let a = body.reciprocal_absorbed_fraction(580.0, 0.5).unwrap();
    assert!(a.into_iter().all(|value| value > 0.0 && 1.0 - value == 1.0));
    assert_eq!(1e-202 * (a[0] + a[1]), 0.0);
    let rate = pair.evaluate(1e300, 300.0, || Ok(())).unwrap();
    let slope = 2.0 * 299_792_458.0 * 1.380_649e-23 * 1e27 / 580.0_f64.powi(4);
    let expected = (slope * 1e300) * (a[0] + a[1]) * (0.5 * 1e-202) * (PI * 0.5 * 400.0);
    close(rate.outward_power_w, expected, 2e-12);
    assert!(rate.outward_power_w > 0.0 && rate.secant_w_k > 0.0);
    close(rate.secant_w_k, rate.outward_power_w / 1e300, 0.0);
}

#[test]
fn g1_black_infrared_gap_recovers_si_planck_integral_and_shared_source_band() {
    let first = emitter(metal((2000.0, 80_000.0), 1.0, 1e-20, 3), 310.0, 4);
    let second_metal = metal((4000.0, 100_000.0), 1.0, 1e-20, 5);
    let second = emitter(second_metal, 250.0, 6);
    let bounds = (4000.0, 80_000.0);
    let area = 0.01;
    let expected = area
        * PI
        * simpson(bounds, 16_384, |lambda| {
            planck(310.0, lambda) - planck(250.0, lambda)
        });
    let coarse = response(&first, &second, area, 64, 4);
    let fine = response(&first, &second, area, 1024, 4);
    assert_eq!(fine.wavelength_bounds_nm(), bounds);
    assert_eq!(fine.optical_evaluations(), 8192);
    assert_eq!(fine.area_m2(), area);
    let coarse_rate = coarse.evaluate(310.0, 250.0, || Ok(())).unwrap();
    let rate = fine.evaluate(310.0, 250.0, || Ok(())).unwrap();
    close(rate.outward_power_w, expected, 2e-5);
    assert!(
        (rate.outward_power_w - expected).abs()
            < (coarse_rate.outward_power_w - expected).abs() / 100.0
    );
    close(rate.secant_w_k, rate.outward_power_w / 60.0, 5e-13);
    let reverse = fine.evaluate(250.0, 310.0, || Ok(())).unwrap();
    assert_eq!(reverse.outward_power_w, -rate.outward_power_w);
    assert_eq!(reverse.secant_w_k, rate.secant_w_k);
    let equilibrium = fine.evaluate(300.0, 300.0, || Ok(())).unwrap();
    assert_eq!(equilibrium.outward_power_w, 0.0);
    let tangent = area
        * PI
        * simpson(bounds, 16_384, |lambda| {
            let x = 6.626_070_15e-34 * 299_792_458.0 / (lambda * 1e-9 * 1.380_649e-23 * 300.0);
            planck(300.0, lambda) * x / (300.0 * -(-x).exp_m1())
        });
    close(equilibrium.secant_w_k, tangent, 2e-5);
    let other_temperature = emitter(second_metal, 1000.0, 6);
    let other = response(&first, &other_temperature, area, 1024, 4);
    assert_eq!(fine.identity(), other.identity());
    assert_eq!(
        fine.optical_identities(),
        [first.optical_identity(), second.optical_identity()]
    );
}

#[test]
fn g1_coherent_spectral_polarized_gap_matches_native_fields_and_independent_bounces() {
    let first_films = [film(2.2, 0.18, 0.7e-6, 10), film(1.4, 0.04, 0.9e-6, 11)];
    let second_film = film(2.0, 0.05, 3e-6, 12);
    let first_layers = [
        CoatingLayer::Absorbing(&first_films[0]),
        CoatingLayer::Absorbing(&first_films[1]),
    ];
    let second_layers = [CoatingLayer::Absorbing(&second_film)];
    let substrate = [
        metal((4000.0, 20_000.0), 0.4, 3.0, 13),
        metal((4000.0, 20_000.0), 0.9, 2.1, 14),
    ];
    let first = CoatingThermalEmission::try_opaque(
        &first_layers,
        &[380.0; 2],
        substrate[0],
        380.0,
        ContentHash([15; 32]),
    )
    .unwrap();
    let second = CoatingThermalEmission::try_opaque(
        &second_layers,
        &[280.0],
        substrate[1],
        280.0,
        ContentHash([16; 32]),
    )
    .unwrap();
    let pair = response(&first, &second, 0.02, 12, 16);
    let stacks = [
        CoherentCoatingStack::try_new(&first_layers).unwrap(),
        CoherentCoatingStack::try_new(&second_layers).unwrap(),
    ];
    let mut expected = 0.0;
    let mut wrong_polarization_average = 0.0;
    for w in 0..12 {
        let width = 16_000.0 / 12.0;
        let lambda = 4000.0 + (w as f64 + 0.5) * width;
        for angle in 0..16 {
            let mu = (angle as f64 + 0.5) / 16.0;
            let a = stacks[0]
                .layer_power_balance(&substrate[0], lambda, 1.0, mu)
                .unwrap();
            let b = stacks[1]
                .layer_power_balance(&substrate[1], lambda, 1.0, mu)
                .unwrap();
            let work = |balance: &crate::conductor::ConductorLayerPowerBalance, pol| {
                balance
                    .layer_absorbed()
                    .iter()
                    .map(|entry| entry[pol])
                    .sum::<f64>()
                    + balance.substrate_absorbed[pol]
            };
            let aa = [work(&a, 0), work(&a, 1)];
            let bb = [work(&b, 0), work(&b, 1)];
            let common =
                0.02 * PI * mu * width / 16.0 * (planck(380.0, lambda) - planck(280.0, lambda));
            expected +=
                common * (reflected_ray_work(aa[0], bb[0]) + reflected_ray_work(aa[1], bb[1]));
            wrong_polarization_average +=
                2.0 * common * reflected_ray_work(0.5 * (aa[0] + aa[1]), 0.5 * (bb[0] + bb[1]));
        }
    }
    let rate = pair.evaluate(380.0, 280.0, || Ok(())).unwrap();
    close(rate.outward_power_w, expected, 5e-12);
    assert!((expected - wrong_polarization_average).abs() > 1e-5 * expected);
    let reversed = response(&second, &first, 0.02, 12, 16);
    assert_eq!(reversed.evaluate(380.0, 280.0, || Ok(())).unwrap(), rate);
    assert_ne!(reversed.optical_identity(), pair.optical_identity());
}

#[test]
fn g4_nonopaque_nonisothermal_domain_budget_and_late_cancellation_refuse_atomically() {
    let first = emitter(metal((4000.0, 20_000.0), 0.4, 3.0, 20), 300.0, 21);
    let second = emitter(metal((8000.0, 16_000.0), 0.8, 2.0, 22), 250.0, 23);
    let ids = [first.identity(), second.identity()];
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 4,
        cosine_bins: 3,
    };
    let make =
        |a, b, area, limits| CoatingGapResponse::from_emitters(a, b, area, q, limits, || Ok(()));
    let nonopaque = CoatingThermalEmission::try_transmitting(
        &[],
        &[],
        CauchyIor::try_constant(1.0).unwrap().into(),
        false,
        ContentHash([24; 32]),
    )
    .unwrap();
    assert_eq!(
        make(&nonopaque, &second, 1.0, budget(q)).unwrap_err(),
        CoatingExchangeError::OpaqueRequired
    );
    assert_eq!(
        make(&first, &nonopaque, 1.0, budget(q)).unwrap_err(),
        CoatingExchangeError::OpaqueRequired
    );
    let layer = film(1.5, 0.1, 1e-6, 25);
    let nonisothermal = CoatingThermalEmission::try_opaque(
        &[CoatingLayer::Absorbing(&layer)],
        &[350.0],
        metal((4000.0, 20_000.0), 0.4, 3.0, 20),
        300.0,
        ContentHash([21; 32]),
    )
    .unwrap();
    assert_eq!(
        make(&nonisothermal, &second, 1.0, budget(q)).unwrap_err(),
        CoatingExchangeError::Nonisothermal
    );
    let separate = emitter(metal((20_000.0, 40_000.0), 0.4, 3.0, 20), 300.0, 21);
    assert_eq!(
        make(&first, &separate, 1.0, budget(q)).unwrap_err(),
        CoatingExchangeError::WavelengthDomain
    );
    let mut insufficient = budget(q);
    insufficient.max_optical_evaluations -= 1;
    assert_eq!(
        make(&first, &second, 1.0, insufficient).unwrap_err(),
        CoatingExchangeError::Budget
    );
    assert_eq!(
        make(&first, &second, 0.0, budget(q)).unwrap_err(),
        CoatingExchangeError::Area
    );
    let mut calls = 0;
    assert_eq!(
        CoatingGapResponse::from_emitters(&first, &second, 1.0, q, budget(q), || {
            calls += 1;
            if calls == 2 * q.wavelength_bins * q.cosine_bins + 2 {
                Err(CoatingExchangeError::Cancelled)
            } else {
                Ok(())
            }
        })
        .unwrap_err(),
        CoatingExchangeError::Cancelled
    );
    assert_eq!([first.identity(), second.identity()], ids);
    let valid = make(&first, &second, 1.0, budget(q)).unwrap();
    assert_eq!(
        valid.evaluate(0.0, 300.0, || Ok(())).unwrap_err(),
        CoatingExchangeError::Thermal(ThermalSurfaceError::Temperature)
    );
    let mut calls = 0;
    assert_eq!(
        valid
            .evaluate(300.0, 250.0, || {
                calls += 1;
                if calls == q.wavelength_bins + 2 {
                    Err(CoatingExchangeError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .unwrap_err(),
        CoatingExchangeError::Cancelled
    );
    assert!(
        valid
            .evaluate(300.0, 250.0, || Ok(()))
            .unwrap()
            .outward_power_w
            > 0.0
    );
}
