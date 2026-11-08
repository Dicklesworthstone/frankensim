use super::*;
use crate::spectral::y_integral;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs() * 2e-12,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn observed_projection_uses_emission_and_transmitted_illumination() {
    let normalization = 1.0 / y_integral();
    let mut sum = XyzAccumulator::default();
    let wavelengths = [430.0, 550.0, 690.0];
    let weights = [20.0, 30.0, 40.0];
    let emission = [2.0, 3.0, 4.0];
    let transmission = [0.5, 0.25, 0.125];
    let background = [7.0, 11.0, 13.0];
    let mut expected = [0.0; 3];
    for lane in 0..3 {
        let observer = XyzObserverLane::new(wavelengths[lane], weights[lane], normalization)
            .expect("visible node");
        sum.add(
            &observer,
            emission[lane],
            transmission[lane],
            -det::ln(transmission[lane]),
            background[lane],
        )
        .expect("finite radiance");
        let radiance = emission[lane] + transmission[lane] * background[lane];
        for (channel, cmf) in [cie_x, cie_y, cie_z].into_iter().enumerate() {
            expected[channel] += radiance * cmf(wavelengths[lane]) * normalization * weights[lane];
        }
    }
    for (actual, expected) in sum.finish().unwrap().into_iter().zip(expected) {
        close(actual, expected);
    }
}

#[test]
fn projection_recovers_opaque_bright_background_and_tiny_bin_weights() {
    let normalization = 1.0 / y_integral();
    let observer = XyzObserverLane::new(550.0, 1.0, normalization).unwrap();
    let mut sum = XyzAccumulator::default();
    sum.add(&observer, 0.0, 0.0, 1000.0, 1e300).unwrap();
    let expected = det::exp(det::ln(1e300) - 1000.0);
    for (actual, cmf) in sum.finish().unwrap().into_iter().zip([cie_x, cie_y, cie_z]) {
        assert!(actual > 0.0);
        close(actual, expected * cmf(550.0) * normalization);
    }

    let weight = f64::from_bits(1);
    let observer = XyzObserverLane::new(780.0, weight, normalization).unwrap();
    assert_eq!(observer.weighted, [0.0; 3]);
    let mut sum = XyzAccumulator::default();
    sum.add(&observer, 1e300, 1.0, 0.0, 0.0).unwrap();
    for (actual, cmf) in sum.finish().unwrap().into_iter().zip([cie_x, cie_y, cie_z]) {
        assert!(actual > 0.0);
        close(actual, (1e300 * weight) * cmf(780.0) * normalization);
    }
}

#[test]
fn projected_terms_can_remain_finite_when_their_spectral_sum_overflows() {
    let normalization = 1.0 / y_integral();
    let observer = XyzObserverLane::new(550.0, 1.0, normalization).unwrap();
    let radiance = f64::MAX * 0.75;
    assert!((radiance + radiance).is_infinite());
    let mut sum = XyzAccumulator::default();
    sum.add(&observer, radiance, 1.0, 0.0, radiance).unwrap();
    for (actual, cmf) in sum.finish().unwrap().into_iter().zip([cie_x, cie_y, cie_z]) {
        close(actual, (radiance * (cmf(550.0) * normalization)) * 2.0);
    }
    let too_large = XyzObserverLane::new(550.0, f64::MAX, normalization).unwrap();
    assert_eq!(
        XyzAccumulator::default().add(&too_large, f64::MAX, 1.0, 0.0, 0.0),
        Err(XyzObserverError::Unrepresentable)
    );
}

#[test]
fn observer_refuses_unsupported_wavelengths_weights_and_invalid_radiance() {
    assert!(matches!(
        XyzObserverLane::new(900.0, 0.0, 1.0),
        Err(XyzObserverError::Wavelength)
    ));
    for weight in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            XyzObserverLane::new(550.0, weight, 1.0),
            Err(XyzObserverError::Weight)
        ));
    }
    for background in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            XyzAccumulator::default().add(&XyzObserverLane::ZERO, 0.0, 1.0, 0.0, background),
            Err(XyzObserverError::Radiance)
        );
    }
}
