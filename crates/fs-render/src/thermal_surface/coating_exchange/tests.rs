use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use crate::dielectric::{CauchyIor, DielectricIor};

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}
fn metal(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|index| {
            ConductorIorSample::try_new(380.0 + 50.0 * index as f64, eta, k).unwrap()
        }),
        source(1),
    )
    .unwrap()
}
fn film(eta: f64, k: f64, thickness_m: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|index| {
            CoatingIorSample::try_new(380.0 + 50.0 * index as f64, eta, k).unwrap()
        }),
        thickness_m,
        source(2),
    )
    .unwrap()
}
fn vacuum() -> DielectricIor {
    CauchyIor::try_constant(1.0).unwrap().into()
}
fn budget(q: CoatingExchangeQuadrature) -> CoatingExchangeBudget {
    CoatingExchangeBudget {
        max_spectral_nodes: q.wavelength_bins,
        max_optical_evaluations: q.wavelength_bins * q.cosine_bins,
    }
}
fn response(
    emitter: &CoatingThermalEmission,
    area: f64,
    bins: [usize; 2],
) -> CoatingExchangeResponse {
    let q = CoatingExchangeQuadrature {
        wavelength_bins: bins[0],
        cosine_bins: bins[1],
    };
    CoatingExchangeResponse::from_emitter(emitter, area, q, budget(q), || Ok(())).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + 2.0 * f64::from_bits(1),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn planck_si(t: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength = wavelength_nm * 1.0e-9;
    2.0 * h * c * c / (wavelength.powi(5) * (h * c / (wavelength * k * t)).exp_m1()) * 1.0e-9
}
fn planck_derivative_si(t: f64, wavelength_nm: f64) -> f64 {
    let x = 6.626_070_15e-34 * 299_792_458.0 / (wavelength_nm * 1.0e-9 * 1.380_649e-23 * t);
    planck_si(t, wavelength_nm) * x / (t * -(-x).exp_m1())
}
fn simpson_band(mut f: impl FnMut(f64) -> f64) -> f64 {
    let n = 8192;
    let width = 400.0 / n as f64;
    let mut sum = f(380.0) + f(780.0);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(380.0 + i as f64 * width);
    }
    sum * width / 3.0
}

#[test]
fn black_half_space_recovers_one_hemisphere_and_independent_finite_band_power() {
    // A matched, vanishing-loss but semi-infinite substrate absorbs all power
    // entering it. It approaches a black interface without inventing gray data.
    let emitter = CoatingThermalEmission::try_opaque(
        &[],
        &[],
        metal(1.0, 1.0e-20),
        1500.0,
        ContentHash([9; 32]),
    )
    .unwrap();
    let area = 0.37;
    let response = response(&emitter, area, [2048, 8]);
    assert_eq!(response.optical_evaluations(), 16_384);
    assert_eq!(response.area_m2(), area);
    let (t, bath) = (1800.0, 1000.0);
    let expected =
        PI * area * simpson_band(|lambda| planck_si(t, lambda) - planck_si(bath, lambda));
    let rates = response.evaluate(t, Some(bath), || Ok(())).unwrap();
    close(rates.outward_power_w, expected, 4.0e-7);
    close(rates.secant_w_k, expected / (t - bath), 4.0e-7);
    let vacuum = response.evaluate(t, None, || Ok(())).unwrap();
    close(
        vacuum.outward_power_w,
        PI * area * simpson_band(|lambda| planck_si(t, lambda)),
        4.0e-7,
    );
    close(vacuum.secant_w_k, vacuum.outward_power_w / t, 3.0e-13);
    let bolometric = area * 5.670_374_419e-8 * (t.powi(4) - bath.powi(4));
    assert!(rates.outward_power_w > 0.0 && rates.outward_power_w < bolometric);
}

#[test]
fn coherent_ports_match_integrated_native_emission_and_reuse_temperature_independent_optics() {
    let absorbing = [film(2.2, 0.15, 100.0e-9), film(1.4, 0.03, 270.0e-9)];
    let clear = CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.45, 0.006, 0.0).unwrap(),
        110.0e-9,
        source(3),
    )
    .unwrap();
    let layers = [
        CoatingLayer::Absorbing(&absorbing[0]),
        CoatingLayer::Dielectric(&clear),
        CoatingLayer::Absorbing(&absorbing[1]),
    ];
    let (t, bath) = (1800.0, 1300.0);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 48,
        cosine_bins: 24,
    };
    let area = 0.02;
    let mut powers = [0.0; 2];
    for (side, reverse) in [false, true].into_iter().enumerate() {
        let make = |temperature| {
            CoatingThermalEmission::try_transmitting(
                &layers,
                &[temperature; 3],
                vacuum(),
                reverse,
                ContentHash([9; 32]),
            )
            .unwrap()
        };
        let body = make(t);
        let environment_temperature_source = make(bath);
        assert_eq!(
            body.optical_identity(),
            environment_temperature_source.optical_identity()
        );
        assert_ne!(body.identity(), environment_temperature_source.identity());
        let compiled =
            CoatingExchangeResponse::from_emitter(&body, area, q, budget(q), || Ok(())).unwrap();
        let other_temperature = CoatingExchangeResponse::from_emitter(
            &environment_temperature_source,
            area,
            q,
            budget(q),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(compiled.identity(), other_temperature.identity());
        let mut independent = 0.0;
        let dlambda = 400.0 / q.wavelength_bins as f64;
        let dmu = 1.0 / q.cosine_bins as f64;
        for i in 0..q.wavelength_bins {
            let wavelength = 380.0 + (i as f64 + 0.5) * dlambda;
            for j in 0..q.cosine_bins {
                let mu = (j as f64 + 0.5) * dmu;
                // This uses the independently tested directional emission owner
                // and its per-polarization Planck sources, not prepared weights.
                let outgoing = body.evaluate(wavelength, mu).unwrap().stokes().i;
                let incoming = environment_temperature_source
                    .evaluate(wavelength, mu)
                    .unwrap()
                    .stokes()
                    .i;
                independent += area * 2.0 * PI * mu * dmu * dlambda * (outgoing - incoming);
            }
        }
        let rates = compiled.evaluate(t, Some(bath), || Ok(())).unwrap();
        close(rates.outward_power_w, independent, 4.0e-12);
        powers[side] = rates.outward_power_w;
        let reversed = compiled.evaluate(bath, Some(t), || Ok(())).unwrap();
        assert_eq!(
            rates.outward_power_w.to_bits(),
            (-reversed.outward_power_w).to_bits()
        );
        assert_eq!(rates.secant_w_k, reversed.secant_w_k);
    }
    assert!((powers[0] / powers[1] - 1.0).abs() > 0.01);
}

#[test]
fn weak_loss_and_adjacent_temperatures_retain_real_exchange_and_positive_secants() {
    let k = 1.0e-20;
    let thickness = 80.0e-9;
    let weak = film(1.0, k, thickness);
    let emitter = CoatingThermalEmission::try_transmitting(
        &[CoatingLayer::Absorbing(&weak)],
        &[1500.0],
        vacuum(),
        false,
        ContentHash([9; 32]),
    )
    .unwrap();
    let area = 0.3;
    let compiled = response(&emitter, area, [1, 1]);
    let cold: f64 = 1700.0;
    let hot = f64::from_bits(cold.to_bits() + 1);
    // To leading order in k, each oblique absorption is alpha*d/mu;
    // integrating the original midpoint mode cancels that path-length cosine.
    let extraction = 2.0 * PI * 4.0 * PI * k * thickness / 580.0e-9;
    let expected_secant = area * extraction * 400.0 * planck_derivative_si(cold, 580.0);
    let rates = compiled.evaluate(hot, Some(cold), || Ok(())).unwrap();
    assert!(rates.outward_power_w > 0.0 && rates.secant_w_k > 0.0);
    close(rates.secant_w_k, expected_secant, 2.0e-11);
    close(
        rates.outward_power_w,
        expected_secant * (hot - cold),
        2.0e-11,
    );
    let equal = compiled.evaluate(cold, Some(cold), || Ok(())).unwrap();
    assert_eq!(equal.outward_power_w, 0.0);
    close(equal.secant_w_k, expected_secant, 2.0e-11);

    // Area*absorption underflows, and gross B overflows. Their final physical
    // product is representable and must survive the two native log scales.
    let small_area = 1.0e-305;
    assert_eq!(small_area * extraction, 0.0);
    let scaled = response(&emitter, small_area, [1, 1]);
    let huge_temperature = 1.0e308;
    let rayleigh_jeans_slope = 2.0 * 299_792_458.0 * 1.380_649e-23 * 1.0e27 / 580.0_f64.powi(4);
    let scaled_rates = scaled.evaluate(huge_temperature, None, || Ok(())).unwrap();
    let expected = (small_area * huge_temperature) * extraction * 400.0 * rayleigh_jeans_slope;
    close(scaled_rates.outward_power_w, expected, 2.0e-11);
    assert!(scaled_rates.secant_w_k > 0.0);
}

#[test]
fn coherent_spectral_and_angular_midpoint_refinement_converges() {
    let absorber = film(2.1, 0.15, 320.0e-9);
    let clear = CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.45, 0.006, 0.0).unwrap(),
        120.0e-9,
        source(3),
    )
    .unwrap();
    let emitter = CoatingThermalEmission::try_opaque(
        &[
            CoatingLayer::Absorbing(&absorber),
            CoatingLayer::Dielectric(&clear),
        ],
        &[2200.0; 2],
        metal(0.8, 3.2),
        2200.0,
        ContentHash([9; 32]),
    )
    .unwrap();
    let power = |bins| {
        response(&emitter, 0.1, [bins, bins])
            .evaluate(2200.0, Some(1200.0), || Ok(()))
            .unwrap()
            .outward_power_w
    };
    let reference = power(256);
    let errors = [16, 32, 64].map(|bins| (power(bins) - reference).abs());
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "{errors:?}"
    );
    assert!(
        errors[2] < 1.0e-4 * reference,
        "{errors:?}, reference={reference}"
    );
}

#[test]
fn admission_and_cancellation_do_not_publish_incomplete_or_nonisothermal_cooling() {
    let optical_film = film(1.7, 0.2, 100.0e-9);
    let make = |film_t, substrate_t| {
        CoatingThermalEmission::try_opaque(
            &[CoatingLayer::Absorbing(&optical_film)],
            &[film_t],
            metal(0.8, 3.0),
            substrate_t,
            ContentHash([9; 32]),
        )
        .unwrap()
    };
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 4,
        cosine_bins: 3,
    };
    assert!(matches!(
        CoatingExchangeResponse::from_emitter(&make(1500.0, 1501.0), 1.0, q, budget(q), || Ok(()),),
        Err(CoatingExchangeError::Nonisothermal)
    ));
    let emitter = make(1500.0, 1500.0);
    let mut checkpoints = 0;
    assert!(matches!(
        CoatingExchangeResponse::from_emitter(
            &emitter,
            1.0,
            q,
            CoatingExchangeBudget {
                max_spectral_nodes: 4,
                max_optical_evaluations: 11
            },
            || {
                checkpoints += 1;
                Ok(())
            },
        ),
        Err(CoatingExchangeError::Budget)
    ));
    assert_eq!(checkpoints, 0);
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            CoatingExchangeResponse::from_emitter(&emitter, area, q, budget(q), || Ok(()),),
            Err(CoatingExchangeError::Area)
        ));
    }
    for invalid in [
        CoatingExchangeQuadrature {
            wavelength_bins: 0,
            cosine_bins: 3,
        },
        CoatingExchangeQuadrature {
            wavelength_bins: usize::MAX,
            cosine_bins: 2,
        },
    ] {
        assert!(matches!(
            CoatingExchangeResponse::from_emitter(&emitter, 1.0, invalid, budget(q), || Ok(()),),
            Err(CoatingExchangeError::Quadrature)
        ));
    }
    let mut checkpoints = 0;
    assert!(matches!(
        CoatingExchangeResponse::from_emitter(&emitter, 1.0, q, budget(q), || {
            checkpoints += 1;
            if checkpoints == 8 {
                Err(CoatingExchangeError::Cancelled)
            } else {
                Ok(())
            }
        },),
        Err(CoatingExchangeError::Cancelled)
    ));
    let compiled = response(&emitter, 1.0, [4, 3]);
    let expected = compiled.evaluate(1700.0, Some(1300.0), || Ok(())).unwrap();
    let mut checkpoints = 0;
    assert_eq!(
        compiled.evaluate(1700.0, Some(1300.0), || {
            checkpoints += 1;
            if checkpoints == 4 {
                Err(CoatingExchangeError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(CoatingExchangeError::Cancelled)
    );
    assert_eq!(
        compiled.evaluate(1700.0, Some(1300.0), || Ok(())).unwrap(),
        expected
    );
    let empty =
        CoatingThermalEmission::try_transmitting(&[], &[], vacuum(), false, ContentHash([9; 32]))
            .unwrap();
    let dark = response(&empty, 1.0, [4, 3]);
    assert_eq!(
        dark.evaluate(1500.0, None, || Ok(())).unwrap(),
        CoatingExchangeRates {
            outward_power_w: 0.0,
            secant_w_k: 0.0
        }
    );
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            dark.evaluate(temperature, None, || Ok(())),
            Err(CoatingExchangeError::Thermal(
                ThermalSurfaceError::Temperature
            ))
        );
        assert_eq!(
            dark.evaluate(1500.0, Some(temperature), || Ok(())),
            Err(CoatingExchangeError::Thermal(
                ThermalSurfaceError::Temperature
            ))
        );
    }
}
