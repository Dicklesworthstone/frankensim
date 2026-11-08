use super::*;
use crate::thermal_surface::SpectralDiffuseThermalSurface;

fn response(epsilon: SpectralEmissivity, subdivisions: usize) -> DiffuseExchangeResponse {
    DiffuseExchangeResponse::from_emissivity(
        epsilon,
        0.03,
        DiffuseExchangeQuadrature {
            subdivisions_per_interval: subdivisions,
        },
        DiffuseExchangeBudget {
            max_spectral_nodes: 4096,
        },
        || Ok(()),
    )
    .unwrap()
}
fn planck(t: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength = wavelength_nm * 1.0e-9;
    2.0 * h * c * c / wavelength.powi(5) / (h * c / (wavelength * k * t)).exp_m1() * 1.0e-9
}
fn integrate(epsilon: SpectralEmissivity, mut f: impl FnMut(f64) -> f64) -> f64 {
    let mut integral = 0.0;
    for pair in epsilon.knots().windows(2) {
        let width = (pair[1].0 - pair[0].0) / 4096.0;
        let mut sum = f(pair[0].0) + f(pair[1].0);
        for index in 1..4096 {
            sum += (if index % 2 == 0 { 2.0 } else { 4.0 }) * f(pair[0].0 + index as f64 * width);
        }
        integral += sum * width / 3.0;
    }
    integral
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= relative * expected.abs().max(1.0e-30),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn spectral_exchange_matches_original_light_and_independent_knot_split_integral() {
    let epsilon =
        SpectralEmissivity::table_in_band(&[(1000.0, 0.2), (4000.0, 0.9), (10_000.0, 0.1)])
            .unwrap();
    let (body, bath) = (1800.0, 900.0);
    let expected = 0.03
        * PI
        * integrate(epsilon, |lambda| {
            epsilon.at(lambda).unwrap() * (planck(body, lambda) - planck(bath, lambda))
        });
    let coarse = response(epsilon, 1)
        .evaluate(body, Some(bath), || Ok(()))
        .unwrap();
    let fine = response(epsilon, 8)
        .evaluate(body, Some(bath), || Ok(()))
        .unwrap();
    close(fine.outward_power_w, expected, 2.0e-9);
    close(fine.secant_w_k, expected / (body - bath), 2.0e-9);
    assert!((fine.outward_power_w - expected).abs() < (coarse.outward_power_w - expected).abs());
    let surface =
        SpectralDiffuseThermalSurface::try_new(body, epsilon, ContentHash([71; 32])).unwrap();
    let emitted = 0.03
        * PI
        * integrate(epsilon, |lambda| {
            surface.radiance_w_m2_sr_nm(lambda).unwrap()
        });
    close(
        response(epsilon, 8)
            .evaluate(body, None, || Ok(()))
            .unwrap()
            .outward_power_w,
        emitted,
        2.0e-9,
    );
}

#[test]
fn equilibrium_tangent_nearby_difference_and_bath_reversal_remain_physical() {
    let epsilon = SpectralEmissivity::constant_in_band(0.6, (500.0, 20_000.0)).unwrap();
    let response = response(epsilon, 16);
    let t: f64 = 2500.0;
    let equal = response.evaluate(t, Some(t), || Ok(())).unwrap();
    assert_eq!(equal.outward_power_w, 0.0);
    assert_eq!(equal.omitted_exchange_upper_bound_w, 0.0);
    let tangent = 0.03
        * PI
        * integrate(epsilon, |lambda| {
            let x = 6.626_070_15e-34 * 299_792_458.0 / (lambda * 1.0e-9 * 1.380_649e-23 * t);
            0.6 * planck(t, lambda) * x / (t * -(-x).exp_m1())
        });
    close(equal.secant_w_k, tangent, 3.0e-9);
    let nearby = response.evaluate(t.next_up(), Some(t), || Ok(())).unwrap();
    assert!(nearby.outward_power_w > 0.0);
    close(nearby.secant_w_k, equal.secant_w_k, 5.0e-13);
    close(
        nearby.outward_power_w,
        equal.secant_w_k * (t.next_up() - t),
        5.0e-13,
    );
    let forward = response.evaluate(t, Some(1800.0), || Ok(())).unwrap();
    let reverse = response.evaluate(1800.0, Some(t), || Ok(())).unwrap();
    assert_eq!(forward.outward_power_w, -reverse.outward_power_w);
    assert_eq!(forward.secant_w_k, reverse.secant_w_k);
    assert_eq!(
        forward.omitted_exchange_upper_bound_w,
        reverse.omitted_exchange_upper_bound_w
    );
    let vacuum = response.evaluate(t, None, || Ok(())).unwrap();
    close(vacuum.secant_w_k, vacuum.outward_power_w / t, 5.0e-13);
}

#[test]
fn source_updates_share_grid_and_finite_band_coverage_cannot_disappear() {
    let epsilon = SpectralEmissivity::constant(0.4).unwrap();
    let original = response(epsilon, 4);
    let dark = original
        .with_emissivity(SpectralEmissivity::constant(0.0).unwrap())
        .unwrap();
    assert!(Arc::ptr_eq(&original.nodes, &dark.nodes));
    assert_ne!(original.identity(), dark.identity());
    let rates = dark.evaluate(2500.0, None, || Ok(())).unwrap();
    assert_eq!(rates.outward_power_w, 0.0);
    assert!(
        rates.omitted_exchange_upper_bound_w > 60_000.0,
        "a dark visible band does not establish a dark infrared source"
    );
    assert!(
        original
            .with_emissivity(SpectralEmissivity::constant_in_band(0.4, (400.0, 780.0)).unwrap())
            .is_err()
    );
    assert!(dark.evaluate(f64::NAN, None, || Ok(())).is_err());
    assert_eq!(
        original.evaluate(2500.0, None, || Err(DiffuseExchangeError::Cancelled)),
        Err(DiffuseExchangeError::Cancelled)
    );
    assert!(
        DiffuseExchangeResponse::from_emissivity(
            epsilon,
            0.03,
            DiffuseExchangeQuadrature {
                subdivisions_per_interval: 2
            },
            DiffuseExchangeBudget {
                max_spectral_nodes: 7
            },
            || Ok(()),
        )
        .is_err()
    );
}
