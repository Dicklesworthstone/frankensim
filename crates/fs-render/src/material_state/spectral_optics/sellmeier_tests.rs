use super::*;
use crate::charts::Ray;
use crate::dielectric::{DielectricEvent, fresnel_dielectric, sample_smooth_dielectric};
use crate::gradient_index::{
    GradientIndexField, GradientIndexOptions, HomogeneousIndexField, propagate,
};
use fs_evidence::ValidityDomain;
use fs_geom::{Point3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES,
    VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES, resolve_spectral_sellmeier_state_point,
};
use fs_qty::{Dims, QuantitySpec, Temperature};

const SURFACE: ContentHash = ContentHash([73; 32]);
const WAVELENGTHS_UM: [f64; 9] = [1.2, 1.25, 1.35, 1.4, 1.5, 1.6, 1.65, 1.7, 1.8];

fn point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature_k).unwrap()
}
fn claim(name: &str, quantity: QuantitySpec, cold: f64, hot: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, quantity),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, cold), (600.0, hot)],
            dims: quantity.dims(),
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("synthetic infrared Sellmeier fixture: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}
fn card(edit: impl Fn(PropertyClaim) -> PropertyClaim) -> MaterialCard {
    let mut sources = Vec::new();
    let length = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
    for i in 0..3 {
        sources.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
            QuantitySpec::dimensional(Dims::NONE),
            [1.1, 0.2, 0.8][i],
            [1.2, 0.2, 0.8][i],
        ));
        let c = [0.01, 0.04, 100.0][i] * 1e-12;
        sources.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
            c,
            c,
        ));
    }
    for (name, value) in SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES
        .into_iter()
        .zip([1e-6, 2e-6])
    {
        sources.push(claim(name, length, value, value));
    }
    for (i, wavelength) in WAVELENGTHS_UM.into_iter().enumerate() {
        sources.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            length,
            wavelength * 1e-6,
            wavelength * 1e-6,
        ));
        sources.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            0.1 + 0.02 * i as f64,
            0.4 + 0.02 * i as f64,
        ));
    }
    let mut claims = ClaimSet::new();
    for source in sources {
        claims.insert_claim(edit(source)).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic absolute infrared dielectric".into(),
            phase: "solid".into(),
            process: "equation verification, not measured glass accuracy".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}
fn constant(source: &mut PropertyClaim, value: f64) {
    source.value = PropertyValue::Scalar {
        value,
        dims: source.key.dims(),
    };
    source.interpolation = InterpolationPolicy::ConstantWithinValidity;
}
fn source(card: &MaterialCard, temperature: f64) -> SpectralSellmeierStatePoint {
    resolve_spectral_sellmeier_state_point(
        card,
        &point(temperature),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn binding(source: &SpectralSellmeierStatePoint) -> MaterialOpticalBinding {
    MaterialOpticalBinding::try_spectral_sellmeier(source.resolved(), source, None, SURFACE)
        .unwrap()
}
fn glass(binding: MaterialOpticalBinding) -> DielectricGlass {
    let OpticalAppearance::Dielectric { glass, .. } = binding.appearance() else {
        panic!("analytic dielectric required")
    };
    assert!(matches!(glass.ior(), DielectricIor::Sellmeier(_)));
    glass
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn reference(wavelength_nm: f64, temperature: f64) -> (f64, f64) {
    // Independent micrometre equation and analytic wavelength derivative.
    let b = [1.1 + (temperature - 300.0) / 3000.0, 0.2, 0.8];
    let c = [0.01, 0.04, 100.0];
    let x = (wavelength_nm / 1000.0).powi(2);
    let n = (1.0 + (0..3).map(|i| b[i] * x / (x - c[i])).sum::<f64>()).sqrt();
    let group = n
        + (0..3)
            .map(|i| b[i] * c[i] * x / (x - c[i]).powi(2))
            .sum::<f64>()
            / n;
    (n, group)
}

#[test]
fn g1_source_analytic_infrared_law_drives_native_snell_fresnel_group_and_beer_transport() {
    let card = card(|source| source);
    let source = source(&card, 450.0);
    let binding = binding(&source);
    let glass = glass(binding);
    let field = HomogeneousIndexField::new(&glass);
    for wavelength_nm in [1_250.0, 1_550.0, 1_700.0] {
        let (n, group) = reference(wavelength_nm, 450.0);
        let sample = field
            .sample(Point3::new(0.0, 0.0, 0.0), wavelength_nm)
            .unwrap();
        near(sample.phase_index(), n, 4e-15);
        near(sample.group_index().unwrap(), group, 4e-15);
        let refracted = sample_smooth_dielectric(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.6, 0.0, 0.8),
            1.0,
            n,
            0.9,
        )
        .unwrap();
        assert_eq!(refracted.event, DielectricEvent::Transmission);
        near(refracted.direction.x, -0.6 / n, 3e-15);
        let cos_t = (1.0 - (0.6 / n).powi(2)).sqrt();
        near(refracted.direction.z, -cos_t, 3e-15);
        let reflectance = 0.5
            * (((0.8 - n * cos_t) / (0.8 + n * cos_t)).powi(2)
                + ((n * 0.8 - cos_t) / (n * 0.8 + cos_t)).powi(2));
        near(
            fresnel_dielectric(0.8, 1.0, n).unwrap().reflectance,
            reflectance,
            3e-15,
        );
        let length = 0.37;
        let path = propagate(
            &field,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: refracted.direction,
            },
            Vec3::new(0.0, 1.0, 0.0),
            wavelength_nm,
            length,
            GradientIndexOptions::default(),
        )
        .unwrap();
        near(path.optical_path_length_m, n * length, 3e-15);
        near(
            path.group_optical_path_length_m.unwrap(),
            group * length,
            3e-15,
        );
        let alpha = source.absorption_per_m_at(wavelength_nm * 1e-9).unwrap();
        near(path.power_transmittance(), (-alpha * length).exp(), 3e-15);
    }
    // Absorption knots do not introduce corners into the original index law.
    for knot in source.absorption_samples() {
        assert!(
            field
                .sample(Point3::new(0.0, 0.0, 0.0), knot.wavelength_m * 1e9)
                .unwrap()
                .group_index()
                .is_some()
        );
    }
    assert_eq!(
        binding.accepted_state_identity(),
        source.resolved().identity()
    );
    for property in source.resolved().properties() {
        card.claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
}

#[test]
fn g0_fit_domain_absorption_domain_and_native_material_identity_remain_distinct() {
    let card = card(|source| source);
    let original = source(&card, 300.0);
    let binding = binding(&original);
    let glass = glass(binding);
    assert_eq!(glass.ior().wavelength_bounds_nm(), (1000.0, 2000.0));
    let (lo, hi) = glass.wavelength_bounds_nm().unwrap();
    near(lo, 1200.0, 3e-13);
    near(hi, 1800.0, 3e-13);
    let field = HomogeneousIndexField::new(&glass);
    for wavelength in [999.0, 1199.0, 1801.0, 2001.0] {
        assert!(
            field
                .sample(Point3::new(0.0, 0.0, 0.0), wavelength)
                .is_err()
        );
    }
    let narrower = self::card(|mut claim| {
        if claim.key.name() == SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[0] {
            constant(&mut claim, 1.1e-6);
        }
        claim
    });
    let narrower = self::binding(&source(&narrower, 300.0));
    assert_ne!(binding.identity(), narrower.identity());
    assert_ne!(
        binding.material().content_identity(),
        narrower.material().content_identity()
    );
    assert_eq!(
        self::glass(narrower).wavelength_bounds_nm(),
        glass.wavelength_bounds_nm()
    );
}

#[test]
fn g0_active_pole_outside_absorption_still_refuses_the_complete_original_fit() {
    let card = card(|mut source| {
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[0] {
            constant(&mut source, 1e-5);
        }
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[1] {
            constant(&mut source, 1.25);
        }
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[2] {
            constant(&mut source, 0.0);
        }
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[0] {
            constant(&mut source, 1.21e-12);
        }
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[1] {
            constant(&mut source, 0.0);
        }
        source
    });
    let source = source(&card, 300.0);
    assert!(
        SellmeierIor::try_new_in_band([1e-5, 1.25, 0.0], [1.21, 0.0, 0.0], (1200.0, 1800.0))
            .is_ok()
    );
    assert!(matches!(
        MaterialOpticalBinding::try_spectral_sellmeier(source.resolved(), &source, None, SURFACE),
        Err(MaterialOpticalError::Dielectric(
            DielectricError::InvalidIor
        ))
    ));
}

#[test]
fn g1_pinned_refresh_updates_native_dispersion_and_absorption_atomically() {
    let card = card(|source| source);
    let original = source(&card, 300.0);
    let binding = binding(&original);
    let accepted = original
        .resolved()
        .at_query_point(&card, &point(600.0))
        .unwrap();
    let source_enum = SpectralOpticalStatePoint::DielectricSellmeier(original.clone());
    let (hot, refreshed) = binding
        .at_spectral_material_state(&source_enum, &card, &accepted)
        .unwrap();
    let cold_glass = glass(binding);
    let hot_glass = glass(hot);
    near(
        hot_glass.ior().eval(1550.0).unwrap(),
        reference(1550.0, 600.0).0,
        4e-15,
    );
    assert!(hot_glass.ior().eval(1550.0).unwrap() > cold_glass.ior().eval(1550.0).unwrap());
    assert!(
        hot_glass.absorption_ref().extinction_per_m(1550.0).unwrap()
            > cold_glass
                .absorption_ref()
                .extinction_per_m(1550.0)
                .unwrap()
    );
    assert_ne!(hot.identity(), binding.identity());
    for (a, b) in source_enum
        .resolved()
        .properties()
        .iter()
        .zip(refreshed.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert_eq!(
        binding.accepted_state_identity(),
        original.resolved().identity()
    );
    assert!(
        MaterialOpticalBinding::try_spectral_sellmeier(&accepted, &original, None, SURFACE)
            .is_err()
    );
}

#[test]
fn g0_refresh_refuses_a_new_pole_without_mutating_the_bound_material() {
    let card = card(|mut source| {
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[0] {
            let PropertyValue::Curve { knots, .. } = &mut source.value else {
                unreachable!()
            };
            knots[1].1 = 1.21e-12;
        }
        source
    });
    let original = source(&card, 300.0);
    let binding = binding(&original);
    let before = binding;
    let accepted = original
        .resolved()
        .at_query_point(&card, &point(600.0))
        .unwrap();
    let source_enum = SpectralOpticalStatePoint::DielectricSellmeier(original.clone());
    assert!(source_enum.at_query_point(&card, &accepted.query()).is_ok());
    assert!(matches!(
        binding.at_spectral_material_state(&source_enum, &card, &accepted),
        Err(MaterialOpticalRefreshError::Binding(
            MaterialOpticalError::Dielectric(DielectricError::InvalidIor)
        ))
    ));
    assert_eq!(binding, before);
}

#[test]
fn g1_explicit_high_index_ir_source_reaches_native_fresnel_snell_and_group_transport() {
    let card = card(|mut source| {
        if let Some(index) = VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES
            .iter()
            .position(|name| *name == source.key.name())
        {
            constant(&mut source, if index == 0 { 15.0 } else { 0.0 });
        }
        if VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES.contains(&source.key.name()) {
            constant(&mut source, 0.0);
        }
        source
    });
    let original = source(&card, 300.0);
    let glass = glass(binding(&original));
    let field = HomogeneousIndexField::new(&glass);
    let sample = field.sample(Point3::new(0.0, 0.0, 0.0), 1550.0).unwrap();
    assert_eq!(sample.phase_index(), 4.0);
    assert_eq!(sample.group_index(), Some(4.0));
    near(
        fresnel_dielectric(1.0, 1.0, 4.0).unwrap().reflectance,
        0.36,
        4e-16,
    );
    let refracted = sample_smooth_dielectric(
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.6, 0.0, 0.8),
        1.0,
        4.0,
        0.9,
    )
    .unwrap();
    assert_eq!(refracted.event, DielectricEvent::Transmission);
    near(refracted.direction.x, -0.15, 4e-16);
    let path = propagate(
        &field,
        Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: refracted.direction,
        },
        Vec3::new(0.0, 1.0, 0.0),
        1550.0,
        0.5,
        GradientIndexOptions::default(),
    )
    .unwrap();
    assert_eq!(path.optical_path_length_m, 2.0);
    assert_eq!(path.group_optical_path_length_m, Some(2.0));
}
