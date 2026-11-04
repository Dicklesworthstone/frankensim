//! Finite-band analytic dispersion, original-owner transport, and source limits.

use super::*;
use crate::charts::Ray;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use crate::gradient_index::{
    GradientIndexError, GradientIndexOptions, HomogeneousIndexField, propagate,
};
use crate::thermal_surface::CoatingThermalEmission;
use fs_blake3::ContentHash;
use fs_geom::Point3;

const STRENGTHS: [f64; 3] = [1.1, 0.2, 0.8];
const RESONANCES_UM2: [f64; 3] = [0.01, 0.04, 100.0];

fn law(bounds: (f64, f64)) -> SellmeierIor {
    SellmeierIor::try_new_in_band(STRENGTHS, RESONANCES_UM2, bounds).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_finite_band_sellmeier_retains_independent_phase_and_group_dispersion() {
    let model = law((1000.0, 2500.0));
    assert_eq!(model.wavelength_bounds_nm(), (1000.0, 2500.0));
    assert_eq!(
        DielectricIor::from(model).wavelength_bounds_nm(),
        (1000.0, 2500.0)
    );
    // Independent 90-digit Decimal evaluation of the original binary64
    // synthetic coefficients. Lambda is converted exactly from nm to um;
    // this is an algebra oracle, not a material-source or measurement claim.
    for (lambda, phase, group) in [
        (1000.0, 1.520_316_952_600_225_6, 1.538_777_795_563_982_3),
        (1064.0, 1.519_199_780_406_712, 1.536_807_392_945_431_3),
        (1310.0, 1.515_668_804_356_159, 1.532_548_843_377_148),
        (1550.0, 1.512_709_759_719_37, 1.531_377_184_970_966_5),
        (2000.0, 1.507_130_970_061_144_2, 1.533_357_378_278_836_3),
        (2500.0, 1.499_905_907_697_543_8, 1.539_875_799_652_972_6),
    ] {
        near(model.eval(lambda).unwrap(), phase, 7.0e-16);
        near(model.group_index(lambda).unwrap(), group, 7.0e-16);
        assert!(group > phase);
    }
    for lambda in [1064.0, 1310.0, 1550.0, 2000.0] {
        let h = 0.05;
        let derivative = (model.eval(lambda - 2.0 * h).unwrap()
            - 8.0 * model.eval(lambda - h).unwrap()
            + 8.0 * model.eval(lambda + h).unwrap()
            - model.eval(lambda + 2.0 * h).unwrap())
            / (12.0 * h);
        near(
            model.group_index(lambda).unwrap(),
            model.eval(lambda).unwrap() - lambda * derivative,
            2.0e-11,
        );
    }
}

#[test]
fn g0_finite_band_sellmeier_refuses_hidden_poles_and_extrapolation() {
    for bounds in [
        (0.0, 2000.0),
        (-1.0, 2000.0),
        (2000.0, 1000.0),
        (f64::NAN, 2000.0),
        (1000.0, f64::INFINITY),
        (f64::NEG_INFINITY, 2000.0),
        (f64::MIN_POSITIVE, 1000.0),
        (1000.0, f64::MAX),
    ] {
        assert_eq!(
            SellmeierIor::try_new_in_band(STRENGTHS, RESONANCES_UM2, bounds),
            Err(DielectricError::InvalidWavelength)
        );
    }
    // An arbitrarily weak active resonance is still forbidden, including
    // one hidden between otherwise ordinary endpoint phase indices.
    for resonance in [1.0, 2.25, 4.0] {
        assert_eq!(
            SellmeierIor::try_new_in_band(
                [1.25, 1.0e-30, 0.0],
                [0.0, resonance, 0.0],
                (1000.0, 2000.0),
            ),
            Err(DielectricError::InvalidIor)
        );
    }
    // Pole-free does not imply an admitted index: the monotone endpoints
    // must still keep the entire band's phase index within [1,3].
    for (strengths, resonances) in [
        ([9.0, 0.0, 0.0], [0.0; 3]),
        ([1.0, 0.0, 0.0], [100.0, 0.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.99, 0.0, 0.0]),
    ] {
        assert_eq!(
            SellmeierIor::try_new_in_band(strengths, resonances, (1000.0, 2000.0)),
            Err(DielectricError::InvalidIor)
        );
    }
    let absent =
        SellmeierIor::try_new_in_band([-0.0; 3], [1.0, 2.25, 4.0], (1000.0, 2000.0)).unwrap();
    assert_eq!(absent.coefficients(), ([0.0; 3], [0.0; 3]));
    assert_eq!(absent.eval(1550.0), Ok(1.0));
    assert_eq!(absent.group_index(1550.0), Ok(1.0));
    let single = law((1550.0, 1550.0));
    assert_eq!(single.eval(1550.0), law((1000.0, 2000.0)).eval(1550.0));
    for lambda in [
        550.0,
        f64::from_bits(1550.0_f64.to_bits() - 1),
        f64::from_bits(1550.0_f64.to_bits() + 1),
        f64::NAN,
        f64::INFINITY,
        -1.0,
    ] {
        assert_eq!(single.eval(lambda), Err(DielectricError::InvalidWavelength));
        assert_eq!(
            single.group_index(lambda),
            Err(DielectricError::InvalidWavelength)
        );
    }
}

#[test]
fn g3_explicit_visible_sellmeier_preserves_original_arithmetic_and_band() {
    let original = SellmeierIor::try_new(STRENGTHS, RESONANCES_UM2).unwrap();
    let explicit = law((LAMBDA_MIN, LAMBDA_MAX));
    let broader = law((300.0, 2500.0));
    assert_eq!(original, explicit);
    assert_ne!(original, broader);
    for lambda in [LAMBDA_MIN, 404.7, 550.0, 710.0, LAMBDA_MAX] {
        assert_eq!(
            original.eval(lambda).unwrap().to_bits(),
            broader.eval(lambda).unwrap().to_bits()
        );
        assert_eq!(
            original.group_index(lambda).unwrap().to_bits(),
            broader.group_index(lambda).unwrap().to_bits()
        );
    }
    for lambda in [300.0, 1000.0, 1550.0, 2500.0] {
        assert_eq!(
            original.eval(lambda),
            Err(DielectricError::InvalidWavelength)
        );
        assert_eq!(
            original.group_index(lambda),
            Err(DielectricError::InvalidWavelength)
        );
        assert!(broader.eval(lambda).is_ok());
    }
}

#[test]
fn g1_infrared_sellmeier_native_refraction_keeps_phase_delay_and_beer_loss() {
    let glass = DielectricGlass::new(
        law((1000.0, 2000.0)),
        BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
            (1100.0, 0.25),
            (1900.0, 0.65),
        ])
        .unwrap(),
        GlassProvenance::Custom,
    );
    assert_eq!(glass.wavelength_bounds_nm(), Ok((1100.0, 1900.0)));
    let field = HomogeneousIndexField::new(&glass);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.6, 0.0, 0.8);
    let length = 0.75;
    // Same independent high-precision synthetic law as the phase/group test.
    for (lambda, phase, group) in [
        (1200.0, 1.517_141_483_751_031_7, 1.533_938_027_661_466_3),
        (1550.0, 1.512_709_759_719_37, 1.531_377_184_970_966_5),
        (1800.0, 1.509_674_681_465_904_3, 1.531_952_518_573_226_3),
    ] {
        let interface = DielectricSurface::SMOOTH.interface(lambda, false).unwrap();
        let index = glass.ior().eval(lambda).unwrap();
        let transmitted = interface
            .sample_smooth(normal, wo, 1.0, index, 0.99)
            .unwrap();
        assert_eq!(transmitted.event, DielectricEvent::Transmission);
        let ct = (1.0_f64 - (0.6_f64 / phase).powi(2)).sqrt();
        let rs = (0.8 - phase * ct) / (0.8 + phase * ct);
        let rp = (phase * 0.8 - ct) / (phase * 0.8 + ct);
        near(transmitted.direction.x, -0.6 / phase, 5.0e-16);
        near(transmitted.direction.z, -ct, 5.0e-16);
        near(
            transmitted.probability,
            1.0 - 0.5 * (rs * rs + rp * rp),
            5.0e-16,
        );
        let path = propagate(
            &field,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: transmitted.direction,
            },
            Vec3::new(0.0, 1.0, 0.0),
            lambda,
            length,
            GradientIndexOptions::default(),
        )
        .unwrap();
        near(path.ray.origin.x, -length * 0.6 / phase, 5.0e-16);
        near(path.ray.origin.z, -length * ct, 5.0e-16);
        near(path.optical_path_length_m, length * phase, 1.0e-15);
        near(
            path.group_optical_path_length_m.unwrap(),
            length * group,
            1.0e-15,
        );
        let alpha = 0.25 + (lambda - 1100.0) * 0.0005;
        near(path.power_transmittance(), (-length * alpha).exp(), 5.0e-16);
        assert_eq!(path.field_evaluations, 1);
    }
    assert_eq!(
        propagate(
            &field,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: normal
            },
            Vec3::new(1.0, 0.0, 0.0),
            1050.0,
            length,
            GradientIndexOptions::default(),
        ),
        Err(GradientIndexError::Dielectric(
            DielectricError::InvalidWavelength
        ))
    );
}

#[test]
fn g1_infrared_sellmeier_film_keeps_analytic_optical_thickness_and_source_band() {
    let film = CoherentDielectricCoating::try_new(
        law((1000.0, 2000.0)),
        200.0e-9,
        ConductorSource::try_new(ContentHash([27; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap();
    let films = [film];
    let stack = DielectricCoatingStack::try_new_in_band(&films, (1100.0, 1900.0)).unwrap();
    for (lambda, phase) in [
        (1200.0, 1.517_141_483_751_031_7),
        (1550.0, 1.512_709_759_719_37),
        (1800.0, 1.509_674_681_465_904_3),
    ] {
        // Independent normal-incidence Airy intensity for a lossless film,
        // using the high-precision phase index and physical thickness.
        let r01 = (1.0 - phase) / (1.0 + phase);
        let r12 = (phase - 2.0) / (phase + 2.0);
        let interference = 2.0 * r01 * r12 * (4.0 * PI * phase * 200.0 / lambda).cos();
        let expected_r =
            (r01 * r01 + r12 * r12 + interference) / (1.0 + (r01 * r12).powi(2) + interference);
        let power = stack.fresnel(lambda, false, 1.0, 1.0, 2.0).unwrap();
        near(power.reflectance, expected_r, 2.0e-14);
        near(power.transmittance, 1.0 - expected_r, 2.0e-14);
    }
    assert!(matches!(
        DielectricCoatingStack::try_new_in_band(&films, (999.0, 1900.0)),
        Err(DielectricError::InvalidWavelength)
    ));
    assert_eq!(
        stack.fresnel(1050.0, false, 1.0, 1.0, 2.0),
        Err(DielectricError::InvalidWavelength)
    );
}

#[test]
fn g3_sellmeier_source_band_reaches_coating_and_thermal_optical_identity() {
    let source =
        ConductorSource::try_new(ContentHash([28; 32]), ConductorDataStatus::Representative)
            .unwrap();
    let coating = |model| CoherentDielectricCoating::try_new(model, 200.0e-9, source).unwrap();
    let visible_metal = ConductorOptics::representative_tungsten();
    let original = SellmeierIor::try_new(STRENGTHS, RESONANCES_UM2).unwrap();
    assert_eq!(
        visible_metal
            .with_coating(coating(original))
            .content_identity(),
        visible_metal
            .with_coating(coating(law((380.0, 780.0))))
            .content_identity()
    );
    let metal = ConductorOptics::try_new_in_band(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(1000.0 + 125.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        125.0,
        source,
    )
    .unwrap();
    let broad = law((900.0, 2100.0));
    let narrow = law((950.0, 2050.0));
    let broad_metal = metal.with_coating(coating(broad));
    let narrow_metal = metal.with_coating(coating(narrow));
    assert_ne!(
        broad_metal.content_identity(),
        narrow_metal.content_identity()
    );
    for lambda in [1000.0, 1550.0, 2000.0] {
        assert_eq!(
            broad_metal.fresnel(lambda, 1.0, 0.6).unwrap().to_bits(),
            narrow_metal.fresnel(lambda, 1.0, 0.6).unwrap().to_bits()
        );
    }
    let absorbing = CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(1000.0 + 125.0 * i as f64, 1.3, 0.03).unwrap()
        }),
        70.0e-9,
        125.0,
        source,
    )
    .unwrap();
    let emitter = |model: SellmeierIor| {
        CoatingThermalEmission::try_transmitting(
            &[CoatingLayer::Absorbing(&absorbing)],
            &[1200.0],
            model.into(),
            false,
            ContentHash([29; 32]),
        )
        .unwrap()
    };
    let broad_emitter = emitter(broad);
    let narrow_emitter = emitter(narrow);
    // The absorbing film clips both emitters to exactly the same common
    // support. Original opposite-medium source bands must still remain bound.
    assert_eq!(broad_emitter.wavelength_bounds_nm(), (1000.0, 2000.0));
    assert_eq!(narrow_emitter.wavelength_bounds_nm(), (1000.0, 2000.0));
    assert_ne!(
        broad_emitter.optical_identity(),
        narrow_emitter.optical_identity()
    );
    for lambda in [1200.0, 1550.0, 1800.0] {
        let a = broad_emitter.evaluate(lambda, 0.6).unwrap().stokes();
        let b = narrow_emitter.evaluate(lambda, 0.6).unwrap().stokes();
        assert!(a.i > 0.0);
        assert_eq!(a, b);
    }
}
