use super::hale_querry_water_25c;
use crate::charts::Ray;
use crate::dielectric::{
    BeerLambertAbsorption, BeerLambertParameters, DielectricError, DielectricEvent,
    DielectricGlass, DielectricIor, DielectricSurface, GlassProvenance, sample_smooth_dielectric,
};
use crate::gradient_index::{GradientIndexOptions, HomogeneousIndexField, propagate};
use crate::tracer::Material;
use fs_geom::{Point3, Vec3};
use fs_material::water_absorption::HaleQuerryWaterAbsorption;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn complete_water_source_uses_native_beer_with_original_edge_semantics() {
    let absorption = BeerLambertAbsorption::hale_querry_water_25c();
    assert_eq!(
        absorption.parameters(),
        BeerLambertParameters::HaleQuerryWater25C
    );
    assert_eq!(absorption.wavelength_bounds_nm(), Some((200.0, 200_000.0)));
    for wavelength in [200.0, 475.0, 512.5, 975.0, 3000.0, 200_000.0] {
        assert_eq!(
            absorption.extinction_per_m(wavelength).unwrap(),
            HaleQuerryWaterAbsorption::extinction_per_m(wavelength).unwrap()
        );
        assert_eq!(absorption.transmittance(wavelength, 0.0), Ok(1.0));
        assert_eq!(absorption.absorbed_fraction(wavelength, 0.0), Ok(0.0));
        assert_eq!(absorption.transmittance(wavelength, f64::MAX), Ok(0.0));
        assert_eq!(absorption.absorbed_fraction(wavelength, f64::MAX), Ok(1.0));
    }
    let alpha = 4.0 * core::f64::consts::PI * 1.0e-9 / 500.0e-9;
    near(
        absorption.transmittance(500.0, 3.0).unwrap(),
        (-3.0 * alpha).exp(),
        3.0e-16,
    );
    near(
        absorption.absorbed_fraction(500.0, 3.0).unwrap(),
        -(-3.0 * alpha).exp_m1(),
        3.0e-16,
    );
    let thin = absorption.absorbed_fraction(500.0, 1.0e-16).unwrap();
    assert!(thin > 0.0);
    near(thin, alpha * 1.0e-16, 2.0e-33);
    assert_eq!(absorption.transmittance(500.0, 1.0e-16), Ok(1.0));
    for distance in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            absorption.transmittance(500.0, distance),
            Err(DielectricError::InvalidAbsorption)
        );
        assert_eq!(
            absorption.absorbed_fraction(500.0, distance),
            Err(DielectricError::InvalidAbsorption)
        );
    }
    for wavelength in [0.0, 199.9, 200_000.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            absorption.transmittance(wavelength, 0.0),
            Err(DielectricError::InvalidWavelength)
        );
        assert_eq!(
            absorption.absorbed_fraction(wavelength, 0.0),
            Err(DielectricError::InvalidWavelength)
        );
    }
}

#[test]
fn reference_water_retains_if97_state_iapws_dispersion_and_support_intersection() {
    use fs_material::water_refraction::IapwsWaterRefraction;
    use fs_qty::{Pressure, Temperature};
    use fs_thermochem::water::if97::{If97Region, If97State};
    let glass = hale_querry_water_25c().unwrap();
    let DielectricIor::Water(law) = *glass.ior() else {
        panic!("the original IAPWS law must remain native");
    };
    let state = If97State::at_pressure_temperature(
        Pressure::new(101_325.0),
        Temperature::new(298.15),
        If97Region::Region1,
    )
    .unwrap();
    assert_eq!(law.temperature_k(), 298.15);
    assert_eq!(law.density_kg_m3(), state.density().value());
    assert_eq!(glass.wavelength_bounds_nm(), Ok((200.0, 1100.0)));
    for wavelength in [200.0, 512.5, 1100.0] {
        let original = IapwsWaterRefraction::prepare(wavelength)
            .unwrap()
            .sample(298.15, state.density().value())
            .unwrap();
        assert_eq!(law.eval(wavelength).unwrap(), original.phase_index);
        assert_eq!(law.group_index(wavelength).unwrap(), original.group_index);
    }
    let field = HomogeneousIndexField::new(&glass);
    for wavelength in [199.99, 1100.01, 3000.0] {
        assert!(glass.ior().eval(wavelength).is_err());
        assert!(
            propagate(
                &field,
                Ray {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    dir: Vec3::new(0.0, 0.0, 1.0)
                },
                Vec3::new(1.0, 0.0, 0.0),
                wavelength,
                0.01,
                GradientIndexOptions::default()
            )
            .is_err()
        );
    }
    // Raw source support does not extend the refraction law's applicability.
    assert!(glass.absorption_ref().extinction_per_m(3000.0).is_ok());
}

#[test]
fn refracted_water_path_uses_physical_length_for_spectral_loss_and_group_delay() {
    let glass = hale_querry_water_25c().unwrap();
    let field = HomogeneousIndexField::new(&glass);
    let depth = 0.5;
    let mut transmitted = [0.0; 3];
    for (i, (wavelength, k)) in [(450.0, 1.02e-9), (650.0, 1.64e-8), (975.0, 3.48e-6)]
        .into_iter()
        .enumerate()
    {
        let n = glass.ior().eval(wavelength).unwrap();
        let ng = glass.ior().group_index(wavelength).unwrap();
        let transmitted_cosine = (1.0 - (0.6 / n).powi(2)).sqrt();
        let entered = sample_smooth_dielectric(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.6, 0.0, 0.8),
            1.0,
            n,
            0.9999,
        )
        .unwrap();
        assert_eq!(entered.event, DielectricEvent::Transmission);
        near(entered.direction.x, -0.6 / n, 2.0e-15);
        let distance = depth / transmitted_cosine;
        let path = propagate(
            &field,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: entered.direction,
            },
            Vec3::new(0.0, 1.0, 0.0),
            wavelength,
            distance,
            GradientIndexOptions::default(),
        )
        .unwrap();
        near(path.ray.origin.z, -depth, 2.0e-15);
        near(path.ray.origin.x, -distance * 0.6 / n, 2.0e-15);
        near(path.optical_path_length_m, n * distance, 2.0e-15);
        near(
            path.group_optical_path_length_m.unwrap(),
            ng * distance,
            2.0e-15,
        );
        let alpha = 4.0 * core::f64::consts::PI * k / (wavelength * 1.0e-9);
        let expected = (-alpha * distance).exp();
        near(path.power_transmittance(), expected, 3.0e-14 * expected);
        assert!((path.power_transmittance() - (-alpha * depth).exp()).abs() > 1.0e-12);
        transmitted[i] = path.power_transmittance();
    }
    assert!(transmitted[0] > transmitted[1] && transmitted[1] > transmitted[2]);
    assert!(transmitted[2] < 1.0e-10);
}

#[test]
fn native_water_identity_binds_complete_absorption_source() {
    let water = hale_querry_water_25c().unwrap();
    let material = |glass| Material::Dielectric {
        glass,
        surface: DielectricSurface::SMOOTH,
    };
    let identity = material(water).content_identity();
    assert_eq!(
        identity,
        material(hale_querry_water_25c().unwrap()).content_identity()
    );
    for absorption in [
        BeerLambertAbsorption::CLEAR,
        BeerLambertAbsorption::try_constant(
            water.absorption_ref().extinction_per_m(550.0).unwrap(),
        )
        .unwrap(),
        BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
            (200.0, 1.1e-7),
            (200_000.0, 0.504),
        ])
        .unwrap(),
    ] {
        assert_ne!(
            identity,
            material(DielectricGlass::new(
                *water.ior(),
                absorption,
                GlassProvenance::Custom
            ))
            .content_identity()
        );
    }
}
