use super::*;
use crate::camera::{DispersiveThinLens, ThinLensElement};
use crate::charts::Ray;
use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorSource};
use crate::dielectric::{
    AirRelativeWavelength, BeerLambertAbsorption, CiddorAirIor, DielectricCoatingStack,
    DielectricEvent, DielectricGlass, DielectricIorSample, DielectricSurface, GlassProvenance,
    fresnel_dielectric, sample_smooth_dielectric,
};
use crate::gradient_index::{
    GradientIndexField, GradientIndexOptions, HomogeneousIndexField, propagate,
};
use crate::material_state::thermo_optic::schott_glass_at_temperature;
use fs_geom::{Point3, Vec3};
use fs_material::air_refraction::CiddorAir;
use fs_material::thermo_optic::ThermoOpticCoefficients;

const B: [f64; 3] = [1.1, 0.25, 0.8];
const C: [f64; 3] = [0.01, 0.03, 100.0];
const D: [f64; 3] = [2e-6, 1.5e-8, -5e-11];
const E: [f64; 2] = [2e-7, 3e-10];

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn reference() -> SellmeierIor {
    SellmeierIor::try_new_in_band(B, C, (400.0, 1700.0)).unwrap()
}

fn law() -> SchottThermoOpticLaw {
    // Manufactured coefficients exercise every term, not a measured specimen.
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: D,
            e: E,
            resonance_um: 0.24,
        },
        300.0,
        [230.0, 400.0],
        [350.0, 1600.0],
        ContentHash([71; 32]),
    )
    .unwrap()
}

// Independent direct differentiation of TIE-19 equation (3), with wavelength
// in micrometres. Native material preparation uses cubic coefficient arrays.
fn oracle(lambda_nm: f64, temperature: f64) -> (f64, f64) {
    let l = lambda_nm * 0.001;
    let l2 = l * l;
    let mut n2 = 1.0;
    let mut dn2 = 0.0;
    for (b, c) in B.into_iter().zip(C) {
        n2 += b * l2 / (l2 - c);
        dn2 -= 2.0 * b * c * l / (l2 - c).powi(2);
    }
    let n0 = n2.sqrt();
    let dn0 = dn2 / (2.0 * n0);
    let dt = temperature - 300.0;
    let q = l2 - 0.24 * 0.24;
    let pole = E[0] * dt + E[1] * dt * dt;
    let f = D[0] * dt + D[1] * dt * dt + D[2] * dt.powi(3) + pole / q;
    let a = (n2 - 1.0) / (2.0 * n0);
    let n = n0 + a * f;
    let dn = dn0 + 0.5 * (1.0 + 1.0 / n2) * dn0 * f - a * 2.0 * l * pole / q.powi(2);
    (n, n - l * dn)
}

#[test]
fn g1_native_thermal_dispersion_matches_phase_and_group_derivatives() {
    for temperature in [230.0, 300.0, 375.0, 400.0] {
        let heated = SchottThermoOpticIor::try_new(reference(), law(), temperature).unwrap();
        let native = DielectricIor::from(heated);
        assert_eq!(native.wavelength_bounds_nm(), (400.0, 1600.0));
        for wavelength in [400.0, 550.0, 1000.0, 1600.0] {
            let (n, ng) = oracle(wavelength, temperature);
            near(native.eval(wavelength).unwrap(), n, 3e-15);
            near(native.group_index(wavelength).unwrap(), ng, 5e-15);
            if temperature == 300.0 {
                assert_eq!(native.eval(wavelength), reference().eval(wavelength));
                assert_eq!(
                    native.group_index(wavelength),
                    reference().group_index(wavelength)
                );
            }
        }
    }
}

#[test]
fn g1_temperature_changes_native_snell_critical_angle_and_lens_power() {
    let cold = SchottThermoOpticIor::try_new(reference(), law(), 230.0).unwrap();
    let warm = cold.at_temperature(400.0).unwrap();
    for wavelength in [450.0, 550.0, 1500.0] {
        let nc = oracle(wavelength, 230.0).0;
        let nw = oracle(wavelength, 400.0).0;
        assert!(nw > nc);
        let sine = 0.5 * (1.0 / nc + 1.0 / nw);
        let cosine = (1.0 - sine * sine).sqrt();
        assert!(
            !fresnel_dielectric(cosine, cold.eval(wavelength).unwrap(), 1.0)
                .unwrap()
                .total_internal_reflection
        );
        assert!(
            fresnel_dielectric(cosine, warm.eval(wavelength).unwrap(), 1.0)
                .unwrap()
                .total_internal_reflection
        );
        for heated in [cold, warm] {
            let n = oracle(wavelength, heated.temperature_k()).0;
            let sample = sample_smooth_dielectric(
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.6, 0.0, 0.8),
                1.0,
                heated.eval(wavelength).unwrap(),
                0.999,
            )
            .unwrap();
            assert_eq!(sample.event, DielectricEvent::Transmission);
            near(-sample.direction.x, 0.6 / n, 2e-15);
            let lens = DispersiveThinLens::try_new(
                550.0,
                &[ThinLensElement::try_new(heated.into(), 20.0, -20.0).unwrap()],
            )
            .unwrap();
            near(
                lens.power_per_m(wavelength).unwrap(),
                40.0 * (n - 1.0),
                2e-13,
            );
        }
    }
}

#[test]
fn g1_native_bulk_path_preserves_thermal_group_delay_and_explicit_absorption() {
    let reference = DielectricGlass::new(
        reference(),
        BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
            (500.0, 0.05),
            (1500.0, 0.15),
        ])
        .unwrap(),
        GlassProvenance::Custom,
    );
    let heated = schott_glass_at_temperature(&reference, &law(), 375.0).unwrap();
    assert_eq!(heated.wavelength_bounds_nm().unwrap(), (500.0, 1500.0));
    assert_eq!(heated.absorption_ref(), reference.absorption_ref());
    for wavelength in [500.0, 750.0, 1500.0] {
        let path = propagate(
            &HomogeneousIndexField::new(&heated),
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            2.0,
            GradientIndexOptions::default(),
        )
        .unwrap();
        let (n, ng) = oracle(wavelength, 375.0);
        near(path.optical_path_length_m, 2.0 * n, 6e-15);
        near(path.group_optical_path_length_m.unwrap(), 2.0 * ng, 1e-14);
        near(
            path.power_transmittance(),
            (-2.0 * wavelength * 1e-4).exp(),
            3e-15,
        );
        assert_eq!(path.field_evaluations, 1);
    }
}

#[test]
fn g1_heated_coherent_film_uses_original_index_and_source_identity() {
    let source =
        ConductorSource::try_new(ContentHash([81; 32]), ConductorDataStatus::Representative)
            .unwrap();
    let cold = SchottThermoOpticIor::try_new(reference(), law(), 230.0).unwrap();
    let warm = cold.at_temperature(400.0).unwrap();
    let film = |ior| CoherentDielectricCoating::try_new(ior, 150e-9, source).unwrap();
    for heated in [cold, warm] {
        let layers = [film(DielectricIor::from(heated))];
        let stack = DielectricCoatingStack::try_new_in_band(&layers, (400.0, 1600.0)).unwrap();
        for wavelength in [450.0, 550.0, 1500.0] {
            let n = oracle(wavelength, heated.temperature_k()).0;
            // Air/film/air at normal incidence: the Airy denominator includes
            // the original thermally shifted optical thickness n*d.
            let r = (1.0 - n) / (1.0 + n);
            let r2 = r * r;
            let phase = 4.0 * core::f64::consts::PI * n * 150.0 / wavelength;
            let expected =
                2.0 * r2 * (1.0 - phase.cos()) / (1.0 + r2 * r2 - 2.0 * r2 * phase.cos());
            let response = stack.fresnel(wavelength, false, 1.0, 1.0, 1.0).unwrap();
            near(response.reflectance, expected, 2e-14);
            near(response.reflectance + response.transmittance, 1.0, 2e-14);
        }
    }
}

#[test]
fn g1_air_relative_reference_is_converted_before_glass_temperature_changes() {
    let air =
        CiddorAirIor::try_new(CiddorAir::try_new(0.0, 450.0).unwrap(), 293.15, 101_325.0).unwrap();
    let relative = AirRelativeSellmeierIor::try_new(
        B,
        C,
        (400.0, 1700.0),
        air,
        AirRelativeWavelength::Air(air),
    )
    .unwrap();
    let heated = SchottThermoOpticIor::try_new(relative, law(), 375.0).unwrap();
    for wavelength in [450.0, 550.0, 1500.0] {
        let original = law()
            .prepare(
                wavelength,
                relative.eval(wavelength).unwrap(),
                Some(relative.group_index(wavelength).unwrap()),
            )
            .unwrap()
            .sample(375.0)
            .unwrap();
        assert_eq!(heated.eval(wavelength).unwrap(), original.phase_index);
        assert_eq!(
            heated.group_index(wavelength).unwrap(),
            original.group_index.unwrap()
        );
    }
    assert_eq!(
        heated.reference_ior(),
        DielectricIor::AirRelativeSellmeier(relative)
    );
}

#[test]
fn g0_support_derivative_and_duplicate_correction_refusals_remain_explicit() {
    let heated = SchottThermoOpticIor::try_new(reference(), law(), 375.0).unwrap();
    for bad in [399.0, 1601.0, f64::NAN, f64::INFINITY] {
        assert_eq!(heated.eval(bad), Err(DielectricError::InvalidWavelength));
        assert_eq!(
            heated.group_index(bad),
            Err(DielectricError::InvalidWavelength)
        );
    }
    for bad in [229.0, 401.0, f64::NAN, f64::INFINITY] {
        assert_eq!(heated.at_temperature(bad), Err(DielectricError::InvalidIor));
    }
    assert_eq!(
        SchottThermoOpticIor::try_new(heated, law(), 375.0),
        Err(DielectricError::InvalidIor)
    );
    let table = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(
                800.0 + 100.0 * i as f64,
                1.5 + 0.001 * (i * i) as f64,
            )
            .unwrap()
        }),
        100.0,
    )
    .unwrap();
    let at_knot = SchottThermoOpticIor::try_new(table, law(), 375.0).unwrap();
    assert!(at_knot.eval(1200.0).is_ok());
    assert_eq!(
        at_knot.group_index(1200.0),
        Err(DielectricError::UndefinedGroupIndex)
    );
    assert!(at_knot.group_index(1250.0).is_ok());
    let disjoint = SellmeierIor::try_new_in_band(B, C, (1700.0, 1800.0)).unwrap();
    assert_eq!(
        SchottThermoOpticIor::try_new(disjoint, law(), 375.0),
        Err(DielectricError::InvalidWavelength)
    );
}

#[test]
fn g3_temperature_cycles_restore_native_identity_without_accumulated_corrections() {
    let cold = SchottThermoOpticIor::try_new(reference(), law(), 230.0).unwrap();
    let warm = cold.at_temperature(400.0).unwrap();
    assert_eq!(warm.at_temperature(230.0).unwrap(), cold);
    let identity = |ior| {
        crate::tracer::Material::Dielectric {
            glass: DielectricGlass::new(ior, BeerLambertAbsorption::CLEAR, GlassProvenance::Custom),
            surface: DielectricSurface::SMOOTH,
        }
        .content_identity()
    };
    assert_ne!(identity(cold), identity(warm));
    assert_eq!(
        identity(cold),
        identity(warm.at_temperature(230.0).unwrap())
    );
}

#[test]
fn g0_unavailable_group_delay_preserves_native_phase_transport() {
    let table = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(800.0 + 100.0 * i as f64, 1.1 + 0.2 * i as f64)
                .unwrap()
        }),
        100.0,
    )
    .unwrap();
    let heated = SchottThermoOpticIor::try_new(table, law(), 300.0).unwrap();
    assert_eq!(heated.eval(1250.0), table.eval(1250.0));
    assert_eq!(
        heated.group_index(1250.0),
        Err(DielectricError::UndefinedGroupIndex)
    );
    let glass = DielectricGlass::new(
        heated,
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    let sample = HomogeneousIndexField::new(&glass)
        .sample(Point3::new(0.0, 0.0, 0.0), 1250.0)
        .unwrap();
    assert_eq!(sample.phase_index(), table.eval(1250.0).unwrap());
    assert_eq!(sample.group_index(), None);

    // A manufactured strong negative spectral correction can invalidate group
    // delay while its phase law remains ordinary and positive at this point.
    let thermal = SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [0.0; 3],
            e: [-0.0049, 0.0],
            resonance_um: 0.0,
        },
        300.0,
        [300.0, 400.0],
        [650.0, 1600.0],
        ContentHash([91; 32]),
    )
    .unwrap();
    let heated = SchottThermoOpticIor::try_new(reference(), thermal, 400.0).unwrap();
    assert!(heated.eval(650.0).unwrap() >= 1.0);
    assert_eq!(
        heated.group_index(650.0),
        Err(DielectricError::UndefinedGroupIndex)
    );
}
