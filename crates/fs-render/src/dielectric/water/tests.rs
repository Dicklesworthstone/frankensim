use super::*;
use crate::charts::Ray;
use crate::dielectric::{
    BeerLambertAbsorption, DielectricEvent, DielectricGlass, DielectricSurface, GlassProvenance,
    fresnel_dielectric, sample_smooth_dielectric,
};
use crate::gradient_index::{GradientIndexOptions, HomogeneousIndexField, propagate};
use fs_geom::{Point3, Vec3};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:.3e}"
    );
}

#[test]
fn native_material_identity_retains_water_source_temperature_density_and_absorption() {
    let identity = |temperature, density, absorption| {
        crate::tracer::Material::Dielectric {
            glass: DielectricGlass::new(
                IapwsWaterIor::try_new(temperature, density).unwrap(),
                BeerLambertAbsorption::try_constant(absorption).unwrap(),
                GlassProvenance::Custom,
            ),
            surface: DielectricSurface::SMOOTH,
        }
        .content_identity()
    };
    let original = identity(300.0, 995.0, 0.1);
    assert_eq!(original, identity(300.0, 995.0, 0.1));
    for changed in [
        identity(301.0, 995.0, 0.1),
        identity(300.0, 996.0, 0.1),
        identity(300.0, 995.0, 0.2),
    ] {
        assert_ne!(original, changed);
    }
}

#[test]
fn original_water_phase_drives_native_snell_fresnel_and_critical_angle() {
    let law = IapwsWaterIor::try_new(298.15, 997.047_435).unwrap();
    let ior = DielectricIor::from(law);
    assert_eq!(ior.wavelength_bounds_nm(), (200.0, 1100.0));
    assert!(ior.is_dispersive());
    for wavelength in [226.5, 589.0, 1013.98] {
        let n = ior.eval(wavelength).unwrap();
        let group = ior.group_index(wavelength).unwrap();
        assert!(group > n);
        let sine = 0.6;
        let cosine = 0.8;
        let sample = sample_smooth_dielectric(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(sine, 0.0, cosine),
            1.0,
            n,
            0.9999,
        )
        .unwrap();
        assert_eq!(sample.event, DielectricEvent::Transmission);
        near(-sample.direction.x, sine / n, 2.0e-15);
        near(
            -sample.direction.z,
            (1.0 - (sine / n).powi(2)).sqrt(),
            2.0e-15,
        );
        assert!((-sample.direction.x - sine / group).abs() > 1.0e-3);
        let interface = fresnel_dielectric(cosine, 1.0, n).unwrap();
        near(
            interface.reflectance + interface.transmittance,
            1.0,
            3.0e-15,
        );
        let reciprocal = fresnel_dielectric(interface.transmitted_cosine, n, 1.0).unwrap();
        near(reciprocal.reflectance, interface.reflectance, 2.0e-15);
        let critical = (1.0 - 1.0 / (n * n)).sqrt();
        assert!(
            fresnel_dielectric(critical - 1.0e-5, n, 1.0)
                .unwrap()
                .total_internal_reflection
        );
        assert!(
            !fresnel_dielectric(critical + 1.0e-5, n, 1.0)
                .unwrap()
                .total_internal_reflection
        );
    }
}

#[test]
fn native_water_transport_retains_distinct_phase_group_and_explicit_absorption() {
    let law = IapwsWaterIor::try_new(298.15, 997.047_435).unwrap();
    let glass = DielectricGlass::new(
        law,
        BeerLambertAbsorption::try_constant(0.07).unwrap(),
        GlassProvenance::Custom,
    );
    let field = HomogeneousIndexField::new(&glass);
    let ray = Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    };
    for wavelength in [226.5, 589.0, 1013.98] {
        let result = propagate(
            &field,
            ray,
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            3.0,
            GradientIndexOptions::default(),
        )
        .unwrap();
        assert_eq!(result.ray.origin, Point3::new(0.0, 0.0, 3.0));
        assert_eq!(result.ray.dir, ray.dir);
        near(
            result.optical_path_length_m,
            3.0 * law.eval(wavelength).unwrap(),
            2.0e-15,
        );
        near(
            result.group_optical_path_length_m.unwrap(),
            3.0 * law.group_index(wavelength).unwrap(),
            2.0e-15,
        );
        near(result.power_transmittance(), (-0.21_f64).exp(), 3.0e-15);
    }
}

#[test]
fn water_state_and_spectral_domains_remain_original_and_vacuum_is_exact() {
    for (t, rho) in [
        (260.0, 998.0),
        (774.0, 998.0),
        (300.0, -1.0),
        (300.0, 1061.0),
        (f64::NAN, 998.0),
    ] {
        assert!(IapwsWaterIor::try_new(t, rho).is_err());
    }
    let vacuum = DielectricIor::from(IapwsWaterIor::try_new(300.0, 0.0).unwrap());
    assert!(!vacuum.is_dispersive());
    for lambda in [200.0, 589.0, 1100.0] {
        assert_eq!(vacuum.eval(lambda), Ok(1.0));
        assert_eq!(vacuum.group_index(lambda), Ok(1.0));
    }
    for lambda in [199.0, 1101.0, f64::NAN] {
        assert_eq!(vacuum.eval(lambda), Err(DielectricError::InvalidWavelength));
    }
    let cold = IapwsWaterIor::try_new(280.0, 999.0).unwrap();
    let warm = IapwsWaterIor::try_new(330.0, 980.0).unwrap();
    assert_eq!(cold.source_identity(), warm.source_identity());
    assert_ne!(cold, warm);
    assert!(cold.eval(589.0).unwrap() > warm.eval(589.0).unwrap());
}
