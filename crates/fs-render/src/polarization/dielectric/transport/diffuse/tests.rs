use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::polarization::dielectric::transport::BeamMedium;
use fs_geom::Point3;

fn beam(direction: Vec3, intensity: f64) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(1.0, 2.0, 3.0),
        PolarizationFrame::try_new(direction, reference_axis(direction)).unwrap(),
        550.0,
        DielectricGlass::new(
            CauchyIor::try_constant(1.4).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        ),
        Stokes::try_new(intensity, intensity, 0.0, 0.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_cosine_proposals_conserve_hemispherical_power_and_depolarize_at_both_poles() {
    for normal in [
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, -1.0),
        unit(Vec3::new(1.0, -2.0, 3.0)).unwrap(),
    ] {
        let incoming = beam(normal.scale(-1.0), 2.0);
        let mut cosine_mean = 0.0;
        let mut reflected_mean = 0.0;
        let count = 256;
        for index in 0..count {
            let draw = [(index as f64 + 0.5) / count as f64, 0.371];
            let sample = incoming
                .sample_depolarizing_diffuse(normal, 0.3, draw)
                .unwrap();
            let flipped = incoming
                .sample_depolarizing_diffuse(normal.scale(-7.0), 0.3, draw)
                .unwrap();
            let cosine = normal.dot(sample.beam.ray().dir);
            assert!(cosine > 0.0);
            assert!((sample.beam.ray().dir.dot(sample.beam.ray().dir) - 1.0).abs() < 1e-14);
            assert!((sample.pdf_solid_angle * core::f64::consts::PI - cosine).abs() < 1e-14);
            let direction = sample.beam.ray().dir;
            let reverse = flipped.beam.ray().dir;
            let difference = Vec3::new(
                direction.x - reverse.x,
                direction.y - reverse.y,
                direction.z - reverse.z,
            );
            assert!(difference.dot(difference) < 1e-28);
            assert_eq!(sample.beam.power(), Stokes::unpolarized(0.6));
            assert_eq!(sample.beam.medium(), incoming.medium());
            assert_eq!(sample.beam.ray().origin, incoming.ray().origin);
            assert_eq!(sample.beam.wavelength_nm(), incoming.wavelength_nm());
            cosine_mean += cosine / count as f64;
            reflected_mean += sample.beam.power().i / count as f64;
        }
        // The cosine-weighted hemisphere has E[cos(theta)] = 2/3.
        assert!((cosine_mean - 2.0 / 3.0).abs() < 2e-5);
        assert!((reflected_mean - 0.6).abs() < 1e-14);
    }
}

#[test]
fn g0_diffuse_rejects_invalid_draws_and_preserves_positive_underflow_as_an_error() {
    let normal = Vec3::new(0.0, 0.0, -1.0);
    let incoming = beam(normal.scale(-1.0), 1.0);
    for reflectance in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            incoming.sample_depolarizing_diffuse(normal, reflectance, [0.25, 0.5]),
            Err(BeamError::InvalidDiffuseSample)
        );
    }
    for value in [-0.1, 1.0, f64::NAN, f64::INFINITY] {
        for draw in [[value, 0.5], [0.5, value]] {
            assert_eq!(
                incoming.sample_depolarizing_diffuse(normal, 0.5, draw),
                Err(BeamError::InvalidDiffuseSample)
            );
        }
    }
    for normal in [Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)] {
        assert_eq!(
            incoming.sample_depolarizing_diffuse(normal, 0.5, [0.5, 0.5]),
            Err(BeamError::InvalidFrame)
        );
    }
    let tiny = beam(Vec3::new(0.0, 0.0, 1.0), f64::from_bits(1));
    assert_eq!(
        tiny.sample_depolarizing_diffuse(normal, 0.25, [0.0, 0.0]),
        Err(BeamError::NumericalRange)
    );
    assert_eq!(
        tiny.sample_depolarizing_diffuse(normal, 0.0, [0.0, 0.0])
            .unwrap()
            .beam
            .power(),
        Stokes::unpolarized(0.0)
    );
    let grazing = incoming
        .sample_depolarizing_diffuse(normal, 1.0, [f64::from_bits(1.0_f64.to_bits() - 1), 0.0])
        .unwrap();
    assert!(grazing.pdf_solid_angle > 0.0);
    assert_eq!(grazing.beam.power(), Stokes::unpolarized(1.0));
}

struct Gradient;
impl GradientIndexField for Gradient {
    fn sample(&self, point: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        GradientIndexSample::try_new(
            1.2 + 0.01 * point.x,
            Vec3::new(0.01, 0.0, 0.0),
            0.2,
            Some(1.3 + 0.01 * point.x),
        )
    }
}

#[test]
fn g1_diffuse_retains_the_exact_gradient_medium_and_requires_native_curved_propagation() {
    let source = beam(Vec3::new(0.0, 0.0, 1.0), 2.0);
    let medium = BeamMedium::try_gradient_at(&Gradient, source.ray().origin, 550.0).unwrap();
    let source = PolarizedBeam::try_new(
        source.ray().origin,
        source.frame(),
        550.0,
        medium,
        source.power(),
    )
    .unwrap();
    let sample = source
        .sample_depolarizing_diffuse(Vec3::new(0.0, 0.0, -1.0), 0.25, [0.2, 0.3])
        .unwrap();
    assert_eq!(sample.beam.medium(), medium);
    assert_eq!(sample.beam.ray().origin, source.ray().origin);
    assert_eq!(sample.beam.power(), Stokes::unpolarized(0.5));
    assert_eq!(
        sample.beam.advance(1.0),
        Err(BeamError::GradientFieldRequired)
    );
}
