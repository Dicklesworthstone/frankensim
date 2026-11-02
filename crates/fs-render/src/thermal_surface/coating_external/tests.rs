use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus,
    ConductorIorSample, ConductorOptics, ConductorSource,
};
use crate::dielectric::{DielectricIor, DielectricIorSample, TabulatedIor};
use crate::thermal_surface::CoatingExchangeResponse;

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}
fn film(bounds: (f64, f64), eta: f64, k: f64, thickness: f64, tag: u8) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                eta,
                k,
            )
            .unwrap()
        }),
        thickness,
        (bounds.1 - bounds.0) / 8.0,
        source(tag),
    )
    .unwrap()
}
fn metal(bounds: (f64, f64), eta: f64, k: f64) -> ConductorOptics {
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
        source(3),
    )
    .unwrap()
}
fn vacuum(bounds: (f64, f64)) -> DielectricIor {
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                1.0,
            )
            .unwrap()
        }),
        (bounds.1 - bounds.0) / 8.0,
    )
    .unwrap()
    .into()
}
fn budget(
    emitter: &CoatingThermalEmission,
    q: CoatingExchangeQuadrature,
) -> CoatingExternalBathBudget {
    CoatingExternalBathBudget {
        max_channel_spectral_nodes: (emitter.film_count()
            + usize::from(emitter.substrate_temperature_k().is_some()))
            * q.wavelength_bins,
        max_optical_evaluations: q.wavelength_bins * q.cosine_bins,
    }
}
fn response(
    emitter: &CoatingThermalEmission,
    area: f64,
    q: CoatingExchangeQuadrature,
) -> CoatingExternalBathResponse {
    CoatingExternalBathResponse::from_emitter(emitter, area, q, budget(emitter, q), || Ok(()))
        .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + 3.0 * f64::from_bits(1),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn planck(t: f64, lambda_nm: f64) -> f64 {
    let lambda = lambda_nm * 1e-9;
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    2.0 * h * c * c / (lambda.powi(5) * (h * c / (lambda * k * t)).exp_m1()) * 1e-9
}
fn derivative(t: f64, lambda_nm: f64) -> f64 {
    let x = 6.626_070_15e-34 * 299_792_458.0 / (lambda_nm * 1e-9 * 1.380_649e-23 * t);
    planck(t, lambda_nm) * x / (t * -(-x).exp_m1())
}

// Independent half-space Fresnel and Beer depth reference. All finite films
// have the substrate's index, so no internal reflection remains. q is the
// positive complex normal index sqrt((eta+ik)^2 - sin(theta)^2).
fn matched_depth_work(eta: f64, k: f64, lambda: f64, mu: f64, depths: [f64; 2]) -> [f64; 3] {
    let nr = eta * eta - k * k;
    let ni = 2.0 * eta * k;
    let z = nr - (1.0 - mu * mu);
    let qr = ((z.hypot(ni) + z) / 2.0).sqrt();
    let qi = ni / (2.0 * qr);
    let s = 4.0 * mu * qr / ((mu + qr).powi(2) + qi * qi);
    let p = 4.0 * mu * (nr * qr + ni * qi) / ((nr * mu + qr).powi(2) + (ni * mu + qi).powi(2));
    let alpha = 4.0 * PI * qi / (lambda * 1e-9);
    let first = -(-alpha * depths[0]).exp_m1();
    let reach_second = (-alpha * depths[0]).exp();
    let second = reach_second * -(-alpha * depths[1]).exp_m1();
    let reach_substrate = (-alpha * (depths[0] + depths[1])).exp();
    [first, second, reach_substrate].map(|fraction| (s + p) * fraction)
}
fn matched_reference(q: CoatingExchangeQuadrature, temperatures: [f64; 3], bath: f64) -> [f64; 3] {
    let mut expected = [0.0; 3];
    let width = 12_000.0 / q.wavelength_bins as f64;
    for i in 0..q.wavelength_bins {
        let lambda = 8000.0 + (i as f64 + 0.5) * width;
        for j in 0..q.cosine_bins {
            let mu = (j as f64 + 0.5) / q.cosine_bins as f64;
            let fractions = matched_depth_work(1.5, 0.3, lambda, mu, [2e-6, 3e-6]);
            for channel in 0..3 {
                expected[channel] += 0.02 * PI * mu * fractions[channel] * width
                    / q.cosine_bins as f64
                    * (planck(temperatures[channel], lambda) - planck(bath, lambda));
            }
        }
    }
    expected
}

#[test]
fn g1_nonisothermal_exterior_channels_match_independent_oblique_depth_and_planck_sources() {
    let films = [
        film((4000.0, 20_000.0), 1.5, 0.3, 2e-6, 1),
        film((6000.0, 22_000.0), 1.5, 0.3, 3e-6, 2),
    ];
    let layers = [
        CoatingLayer::Absorbing(&films[0]),
        CoatingLayer::Absorbing(&films[1]),
    ];
    let make = |t: [f64; 3]| {
        CoatingThermalEmission::try_opaque(
            &layers,
            &t[..2],
            metal((8000.0, 24_000.0), 1.5, 0.3),
            t[2],
            ContentHash([9; 32]),
        )
        .unwrap()
    };
    let temperatures = [280.0, 370.0, 330.0];
    let emitter = make(temperatures);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 48,
        cosine_bins: 32,
    };
    let response = response(&emitter, 0.02, q);
    assert_eq!(response.wavelength_bounds_nm(), (8000.0, 20_000.0));
    assert_eq!(response.film_count(), 2);
    assert!(response.has_substrate());
    assert_eq!(response.channel_spectral_nodes(), 3 * q.wavelength_bins);
    assert_eq!(
        response.optical_evaluations(),
        q.wavelength_bins * q.cosine_bins
    );
    let rates = response
        .evaluate(
            &temperatures[..2],
            Some(temperatures[2]),
            Some(310.0),
            || Ok(()),
        )
        .unwrap();
    let expected = matched_reference(q, temperatures, 310.0);
    assert!(rates.film_rates()[0].outward_power_w < 0.0);
    assert!(rates.film_rates()[1].outward_power_w > 0.0);
    for (channel, rate) in rates
        .film_rates()
        .iter()
        .copied()
        .chain(rates.substrate_rate())
        .enumerate()
    {
        close(rate.outward_power_w, expected[channel], 4e-11);
        close(
            rate.secant_w_k,
            expected[channel] / (temperatures[channel] - 310.0),
            4e-11,
        );
    }
    close(rates.total_outward_power_w(), expected.iter().sum(), 4e-11);
    let other = CoatingExternalBathResponse::from_emitter(
        &make([400.0, 200.0, 500.0]),
        0.02,
        q,
        budget(&emitter, q),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(other.identity(), response.identity());
    assert_eq!(response.optical_identity(), emitter.optical_identity());

    let reference = matched_reference(
        CoatingExchangeQuadrature {
            wavelength_bins: 512,
            cosine_bins: 256,
        },
        temperatures,
        310.0,
    );
    let mut errors = Vec::new();
    for count in [8, 64] {
        let prepared = CoatingExternalBathResponse::from_emitter(
            &emitter,
            0.02,
            CoatingExchangeQuadrature {
                wavelength_bins: count,
                cosine_bins: count,
            },
            CoatingExternalBathBudget {
                max_channel_spectral_nodes: 3 * count,
                max_optical_evaluations: count * count,
            },
            || Ok(()),
        )
        .unwrap();
        let r = prepared
            .evaluate(
                &temperatures[..2],
                Some(temperatures[2]),
                Some(310.0),
                || Ok(()),
            )
            .unwrap();
        errors.push(
            r.film_rates()
                .iter()
                .copied()
                .chain(r.substrate_rate())
                .enumerate()
                .map(|(i, rate)| (rate.outward_power_w - reference[i]).abs())
                .sum::<f64>(),
        );
    }
    assert!(errors[1] < errors[0] * 0.03, "{errors:?}");
}

#[test]
fn g1_isothermal_reduction_and_asymmetric_reverse_ports_preserve_authored_channels() {
    let bounds = (6000.0, 18_000.0);
    let films = [
        film(bounds, 2.2, 0.4, 2e-6, 1),
        film(bounds, 1.3, 0.1, 4e-6, 2),
    ];
    let layers = [
        CoatingLayer::Absorbing(&films[0]),
        CoatingLayer::Absorbing(&films[1]),
    ];
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 32,
        cosine_bins: 24,
    };
    let area = 0.01;
    for reverse in [false, true] {
        let emitter = CoatingThermalEmission::try_transmitting(
            &layers,
            &[400.0; 2],
            vacuum(bounds),
            reverse,
            ContentHash([10; 32]),
        )
        .unwrap();
        let prepared = response(&emitter, area, q);
        let scalar = CoatingExchangeResponse::from_emitter(
            &emitter,
            area,
            q,
            CoatingExchangeBudget {
                max_spectral_nodes: q.wavelength_bins,
                max_optical_evaluations: q.wavelength_bins * q.cosine_bins,
            },
            || Ok(()),
        )
        .unwrap();
        let rates = prepared
            .evaluate(&[400.0; 2], None, Some(300.0), || Ok(()))
            .unwrap();
        let total = scalar.evaluate(400.0, Some(300.0), || Ok(())).unwrap();
        close(rates.total_outward_power_w(), total.outward_power_w, 8e-13);
        close(
            rates.film_rates().iter().map(|r| r.secant_w_k).sum(),
            total.secant_w_k,
            8e-13,
        );
        let equilibrium = prepared
            .evaluate(&[300.0; 2], None, Some(300.0), || Ok(()))
            .unwrap();
        assert_eq!(equilibrium.total_outward_power_w(), 0.0);
        assert!(
            equilibrium
                .film_rates()
                .iter()
                .all(|r| r.outward_power_w == 0.0 && r.secant_w_k > 0.0)
        );
    }
    let opaque = CoatingThermalEmission::try_opaque(
        &layers,
        &[400.0; 2],
        metal(bounds, 0.7, 3.0),
        400.0,
        ContentHash([11; 32]),
    )
    .unwrap();
    let scalar = CoatingExchangeResponse::from_emitter(
        &opaque,
        area,
        q,
        CoatingExchangeBudget {
            max_spectral_nodes: q.wavelength_bins,
            max_optical_evaluations: q.wavelength_bins * q.cosine_bins,
        },
        || Ok(()),
    )
    .unwrap();
    let layered = response(&opaque, area, q)
        .evaluate(&[400.0; 2], Some(400.0), None, || Ok(()))
        .unwrap();
    close(
        layered.total_outward_power_w(),
        scalar
            .evaluate(400.0, None, || Ok(()))
            .unwrap()
            .outward_power_w,
        8e-13,
    );

    let reversed = CoatingThermalEmission::try_transmitting(
        &layers,
        &[280.0, 430.0],
        vacuum(bounds),
        true,
        ContentHash([10; 32]),
    )
    .unwrap();
    let mirrored_layers = [layers[1], layers[0]];
    let mirrored = CoatingThermalEmission::try_transmitting(
        &mirrored_layers,
        &[430.0, 280.0],
        vacuum(bounds),
        false,
        ContentHash([10; 32]),
    )
    .unwrap();
    let a = response(&reversed, area, q)
        .evaluate(&[280.0, 430.0], None, Some(310.0), || Ok(()))
        .unwrap();
    let b = response(&mirrored, area, q)
        .evaluate(&[430.0, 280.0], None, Some(310.0), || Ok(()))
        .unwrap();
    for i in 0..2 {
        close(
            a.film_rates()[i].outward_power_w,
            b.film_rates()[1 - i].outward_power_w,
            5e-13,
        );
    }
    assert_ne!(reversed.optical_identity(), mirrored.optical_identity());
}

// Independent bounded one-film amplitude formula in equal vacuum ports.
// It is the closed geometric amplitude sum from Fresnel boundary conditions,
// not the production multilayer transfer/work integrator (Byrnes §§3--4).
#[derive(Clone, Copy)]
struct C(f64, f64);
impl C {
    fn add(self, b: Self) -> Self {
        Self(self.0 + b.0, self.1 + b.1)
    }
    fn sub(self, b: Self) -> Self {
        Self(self.0 - b.0, self.1 - b.1)
    }
    fn mul(self, b: Self) -> Self {
        Self(self.0 * b.0 - self.1 * b.1, self.0 * b.1 + self.1 * b.0)
    }
    fn div(self, b: Self) -> Self {
        let norm = b.norm();
        Self(
            (self.0 * b.0 + self.1 * b.1) / norm,
            (self.1 * b.0 - self.0 * b.1) / norm,
        )
    }
    fn norm(self) -> f64 {
        self.0 * self.0 + self.1 * self.1
    }
}
fn single_film_absorptance(lambda: f64, mu: f64) -> f64 {
    let (eta, k, d) = (1.7, 0.2, 2e-6);
    let n2 = C(eta * eta - k * k, 2.0 * eta * k);
    let z = n2.0 - (1.0 - mu * mu);
    let qr = ((z.hypot(n2.1) + z) / 2.0).sqrt();
    let qi = n2.1 / (2.0 * qr);
    let q = C(qr, qi);
    let delta = 2.0 * PI * d / (lambda * 1e-9);
    let e = C(
        (-delta * qi).exp() * (delta * qr).cos(),
        (-delta * qi).exp() * (delta * qr).sin(),
    );
    let e2 = e.mul(e);
    let one = C(1.0, 0.0);
    let mut total = 0.0;
    for entrance in [C(mu, 0.0), C(n2.0 * mu, n2.1 * mu)] {
        let r = entrance.sub(q).div(entrance.add(q));
        let r2 = r.mul(r);
        let denominator = one.sub(r2.mul(e2));
        let reflected = r.mul(one.sub(e2)).div(denominator);
        let transmitted = one.sub(r2).mul(e).div(denominator);
        total += 1.0 - reflected.norm() - transmitted.norm();
    }
    total
}

#[test]
fn g1_single_absorbing_film_between_two_baths_closes_complete_radiative_boundary() {
    let bounds = (6000.0, 18_000.0);
    let film = film(bounds, 1.7, 0.2, 2e-6, 1);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 48,
        cosine_bins: 32,
    };
    let (temperature, area) = (310.0, 0.02);
    let baths = [270.0, 360.0];
    let mut total = 0.0;
    let mut independent = 0.0;
    for (reverse, bath) in [false, true].into_iter().zip(baths) {
        let emitter = CoatingThermalEmission::try_transmitting(
            &[CoatingLayer::Absorbing(&film)],
            &[temperature],
            vacuum(bounds),
            reverse,
            ContentHash([12; 32]),
        )
        .unwrap();
        let rates = response(&emitter, area, q)
            .evaluate(&[temperature], None, Some(bath), || Ok(()))
            .unwrap();
        assert!(rates.substrate_rate().is_none());
        assert_eq!(
            rates.film_rates()[0].outward_power_w.is_sign_positive(),
            temperature > bath
        );
        let width = (bounds.1 - bounds.0) / q.wavelength_bins as f64;
        let mut expected = 0.0;
        for i in 0..q.wavelength_bins {
            let lambda = bounds.0 + (i as f64 + 0.5) * width;
            for j in 0..q.cosine_bins {
                let mu = (j as f64 + 0.5) / q.cosine_bins as f64;
                expected += area * PI * mu * single_film_absorptance(lambda, mu) * width
                    / q.cosine_bins as f64
                    * (planck(temperature, lambda) - planck(bath, lambda));
            }
        }
        close(rates.film_rates()[0].outward_power_w, expected, 5e-12);
        close(
            rates.film_rates()[0].secant_w_k,
            expected / (temperature - bath),
            5e-12,
        );
        total += rates.total_outward_power_w();
        independent += expected;
    }
    close(total, independent, 5e-12);
}

#[test]
fn g1_weak_layer_work_and_scaled_exchange_survive_rounded_flux_and_area_products() {
    let bounds = (380.0, 780.0);
    let film = film(bounds, 1.0, 1e-20, 100e-9, 1);
    let emitter = CoatingThermalEmission::try_transmitting(
        &[CoatingLayer::Absorbing(&film)],
        &[1800.0],
        vacuum(bounds),
        false,
        ContentHash([13; 32]),
    )
    .unwrap();
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 1,
        cosine_bins: 1,
    };
    let work = emitter.reciprocal_layer_absorption(580.0, 0.5).unwrap();
    let a = work.film_absorbed()[0];
    assert!(a.into_iter().all(|v| v > 0.0 && 1.0 - v == 1.0));
    let prepared = response(&emitter, 0.02, q);
    let t = 1800.0_f64;
    let adjacent = f64::from_bits(t.to_bits() + 1);
    let rate = prepared
        .evaluate(&[adjacent], None, Some(t), || Ok(()))
        .unwrap()
        .film_rates()[0];
    let expected_slope = 0.02 * PI * 0.5 * (a[0] + a[1]) * 400.0 * derivative(t, 580.0);
    assert!(rate.outward_power_w > 0.0 && rate.secant_w_k > 0.0);
    close(rate.secant_w_k, expected_slope, 3e-12);
    close(rate.outward_power_w, expected_slope * (adjacent - t), 3e-12);

    let opaque = CoatingThermalEmission::try_opaque(
        &[],
        &[],
        metal(bounds, 1e-120, 100.0),
        300.0,
        ContentHash([14; 32]),
    )
    .unwrap();
    let absorbed = opaque.reciprocal_absorbed_fraction(580.0, 0.5).unwrap();
    assert_eq!(1e-202 * (absorbed[0] + absorbed[1]), 0.0);
    let rate = response(&opaque, 1e-202, q)
        .evaluate(&[], Some(1e300), Some(300.0), || Ok(()))
        .unwrap()
        .substrate_rate()
        .unwrap();
    let rayleigh_jeans_slope = 2.0 * 299_792_458.0 * 1.380_649e-23 * 1e27 / 580.0_f64.powi(4);
    let expected =
        (rayleigh_jeans_slope * 1e300) * (absorbed[0] + absorbed[1]) * 1e-202 * (PI * 0.5 * 400.0);
    close(rate.outward_power_w, expected, 2e-12);
    assert!(rate.secant_w_k > 0.0 && rate.secant_w_k < f64::MIN_POSITIVE);
    close(rate.secant_w_k, rate.outward_power_w / 1e300, 0.0);
}

#[test]
fn g4_empty_dark_channels_total_slot_budgets_shapes_and_late_cancellation_remain_binding() {
    let bounds = (6000.0, 18_000.0);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 3,
        cosine_bins: 2,
    };
    let empty = CoatingThermalEmission::try_transmitting(
        &[],
        &[],
        vacuum(bounds),
        false,
        ContentHash([15; 32]),
    )
    .unwrap();
    let dark = response(&empty, 0.02, q);
    assert_eq!(dark.channel_spectral_nodes(), 0);
    let rates = dark.evaluate(&[], None, Some(500.0), || Ok(())).unwrap();
    assert!(rates.film_rates().is_empty());
    assert_eq!(rates.total_outward_power_w(), 0.0);
    assert!(dark.evaluate(&[], None, Some(0.0), || Ok(())).is_err());
    let films = [
        film(bounds, 1.5, 0.0, 1e-6, 1),
        film(bounds, 2.0, 0.4, 0.0, 2),
    ];
    let emitter = CoatingThermalEmission::try_transmitting(
        &[
            CoatingLayer::Absorbing(&films[0]),
            CoatingLayer::Absorbing(&films[1]),
        ],
        &[300.0, 400.0],
        vacuum(bounds),
        true,
        ContentHash([16; 32]),
    )
    .unwrap();
    let original_identity = emitter.identity();
    let good = budget(&emitter, q);
    let prepared = response(&emitter, 0.02, q);
    assert_eq!(prepared.channel_spectral_nodes(), 6);
    assert_eq!(prepared.optical_evaluations(), 6);
    let rates = prepared
        .evaluate(&[300.0, 400.0], None, Some(500.0), || Ok(()))
        .unwrap();
    assert_eq!(rates.film_rates(), &[ZERO_RATE; 2]);
    for bad in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        assert!(
            prepared
                .evaluate(&[300.0, bad], None, None, || Ok(()))
                .is_err()
        );
    }
    assert_eq!(
        prepared
            .evaluate(&[300.0], None, None, || Ok(()))
            .unwrap_err(),
        CoatingExchangeError::TemperatureChannels
    );
    assert_eq!(
        prepared
            .evaluate(&[300.0; 2], Some(300.0), None, || Ok(()))
            .unwrap_err(),
        CoatingExchangeError::TemperatureChannels
    );
    for bad in [
        CoatingExternalBathBudget {
            max_channel_spectral_nodes: 5,
            ..good
        },
        CoatingExternalBathBudget {
            max_optical_evaluations: 5,
            ..good
        },
    ] {
        let mut calls = 0;
        assert_eq!(
            CoatingExternalBathResponse::from_emitter(&emitter, 0.02, q, bad, || {
                calls += 1;
                Ok(())
            })
            .unwrap_err(),
            CoatingExchangeError::Budget
        );
        assert_eq!(calls, 0);
    }
    assert_eq!(
        CoatingExternalBathResponse::from_emitter(&emitter, 0.0, q, good, || Ok(())).unwrap_err(),
        CoatingExchangeError::Area
    );
    assert_eq!(
        CoatingExternalBathResponse::from_emitter(
            &emitter,
            0.02,
            CoatingExchangeQuadrature {
                wavelength_bins: usize::MAX,
                cosine_bins: 1
            },
            good,
            || Ok(())
        )
        .unwrap_err(),
        CoatingExchangeError::Budget
    );
    let mut preparation_calls = 0;
    let completed = CoatingExternalBathResponse::from_emitter(&emitter, 0.02, q, good, || {
        preparation_calls += 1;
        Ok(())
    })
    .unwrap();
    let mut calls = 0;
    assert_eq!(
        CoatingExternalBathResponse::from_emitter(&emitter, 0.02, q, good, || {
            calls += 1;
            if calls == preparation_calls {
                Err(CoatingExchangeError::Cancelled)
            } else {
                Ok(())
            }
        })
        .unwrap_err(),
        CoatingExchangeError::Cancelled
    );
    let mut evaluation_calls = 0;
    completed
        .evaluate(&[300.0, 400.0], None, Some(500.0), || {
            evaluation_calls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    assert_eq!(
        completed
            .evaluate(&[300.0, 400.0], None, Some(500.0), || {
                calls += 1;
                if calls == evaluation_calls {
                    Err(CoatingExchangeError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .unwrap_err(),
        CoatingExchangeError::Cancelled
    );
    assert_eq!(emitter.identity(), original_identity);
    assert_eq!(completed.identity(), prepared.identity());
    assert_eq!(
        completed
            .evaluate(&[300.0, 400.0], None, Some(500.0), || Ok(()))
            .unwrap(),
        rates
    );
}
