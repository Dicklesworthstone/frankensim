//! Authored original absorption claims on the native curved thermal rays.

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
    SpectralDielectricAbsorptionStatePoint, SpectralDielectricAbsorptionTemperaturePath,
    resolve_spectral_dielectric_absorption_state_point,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn source() -> (MaterialCard, SpectralDielectricAbsorptionStatePoint) {
    let mut claims = ClaimSet::new();
    for i in 0..9 {
        let wavelength_nm = 1000.0 + 125.0 * i as f64;
        let alpha_300 = 0.6 + 0.0004 * (wavelength_nm - 1000.0);
        for (name, q, value, interpolation) in [
            (
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                PropertyValue::Scalar {
                    value: wavelength_nm * 1.0e-9,
                    dims: Dims([1, 0, 0, 0, 0, 0]),
                },
                InterpolationPolicy::ConstantWithinValidity,
            ),
            (
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(280.0, alpha_300 - 0.4), (360.0, alpha_300 + 1.2)],
                    dims: SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY.dims(),
                },
                InterpolationPolicy::LinearInside,
            ),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, q),
                    value,
                    validity: ValidityDomain::unconstrained().with_quantity(
                        "T",
                        temperature_quantity(),
                        280.0,
                        360.0,
                    ),
                    interpolation,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "authored thermal absorption ray oracle".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored independent absorption".into(),
            phase: "solid".into(),
            process: "numerical optical-depth oracle".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), 300.0)
        .unwrap();
    let source = resolve_spectral_dielectric_absorption_state_point(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, source)
}

pub(in crate::material_state::thermo_optic) fn absorption_path(
    range: [f64; 2],
) -> SpectralDielectricAbsorptionTemperaturePath {
    let (card, source) = source();
    SpectralDielectricAbsorptionTemperaturePath::try_new(
        &card,
        source,
        "T",
        range,
        128,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_curved_ray_and_bundle_integrate_original_local_absorption_with_analytic_optical_depth() {
    let glass = glass(BeerLambertAbsorption::try_constant(0.8).unwrap());
    let law = law([0.012, 0.0, 0.0], [0.0045, 0.0]);
    let temperature = temperature(Vec3::new(40.0, 0.0, 0.0), [0.0; 3]);
    let absorption = absorption_path([280.0, 360.0]);
    let field = MaterialThermoOpticField::try_with_absorption_source(
        &glass,
        &law,
        &temperature,
        &absorption,
        ContentHash([112; 32]),
    )
    .unwrap();
    let frozen = MaterialThermoOpticField::try_new(
        &glass,
        &law,
        &temperature,
        ContentHash([112; 32]),
    )
    .unwrap();
    assert!(core::ptr::eq(field.absorption_source().unwrap(), &absorption));
    assert_eq!(field.temperature_range_k(), [280.0, 360.0]);
    let prepared = field.prepare(1500.0).unwrap();
    let length = 0.8;
    let path = propagate(
        &prepared,
        ray(),
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        tight(),
    )
    .unwrap();
    let bundle = propagate_bundle(
        &prepared,
        ray(),
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        BundleOptions {
            integration: tight(),
            ..BundleOptions::default()
        },
    )
    .unwrap();
    // n = 1.5 + g*x, alpha = .8 + .02*(40*x). Eliminate x using the
    // independent exact affine-index OPL; no sampled path quadrature is used.
    let g = 7.0 / 30.0;
    let n_end = (2.25 + (g * length).powi(2)).sqrt();
    let z = (1.5 / g) * (g * length / 1.5).asinh();
    let phase = 0.5 * (length * n_end + 1.5 * z);
    let optical_depth = 0.8 * length + (0.02 * 40.0 / g) * (phase - 1.5 * length);
    near(path.ray.origin.x, (n_end - 1.5) / g, 2.0e-11);
    near(path.optical_path_length_m, phase, 2.0e-11);
    near(path.optical_depth, optical_depth, 1.0e-10);
    near(bundle.propagation().optical_depth, optical_depth, 1.0e-10);
    near(path.power_transmittance(), (-optical_depth).exp(), 6.0e-11);
    near(bundle.propagation().power_transmittance(), path.power_transmittance(), 8.0e-11);
    let held = propagate(
        &frozen.prepare(1500.0).unwrap(),
        ray(),
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        tight(),
    )
    .unwrap();
    near(held.optical_depth, 0.8 * length, 2.0e-12);
    assert!(path.optical_depth - held.optical_depth > 0.005);
    near(held.optical_path_length_m, path.optical_path_length_m, 4.0e-11);

    #[cfg(feature = "photometry")]
    {
        use crate::photometry::scene::{NativeMaterialFieldBinding, NativeOpticalMedium};
        let first = NativeMaterialFieldBinding::from_gradient_field(
            NativeOpticalMedium::Ambient,
            &field,
        );
        let second = NativeMaterialFieldBinding::from_differentiable_field(
            NativeOpticalMedium::Ambient,
            &field,
        );
        let held = NativeMaterialFieldBinding::from_gradient_field(
            NativeOpticalMedium::Ambient,
            &frozen,
        );
        assert_eq!(first.field_identity(), second.field_identity());
        assert_ne!(first.field_identity(), held.field_identity());
        assert_eq!(first.source_snapshot_identity(), held.source_snapshot_identity());
    }
}

#[test]
fn g0_replaced_reference_loss_does_not_narrow_phase_support_or_bypass_source_limits() {
    let glass = glass(
        BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
            (1400.0, 0.2),
            (1800.0, 0.6),
        ])
        .unwrap(),
    );
    let law = law([0.012, 0.0, 0.0], [0.0045, 0.0]);
    let temperature = temperature(Vec3::new(40.0, 0.0, 0.0), [0.0; 3]);
    let absorption = absorption_path([280.0, 360.0]);
    let field = MaterialThermoOpticField::try_with_absorption_source(
        &glass,
        &law,
        &temperature,
        &absorption,
        ContentHash([113; 32]),
    )
    .unwrap();
    let frozen = MaterialThermoOpticField::try_new(
        &glass,
        &law,
        &temperature,
        ContentHash([113; 32]),
    )
    .unwrap();
    assert!(frozen.prepare(1100.0).is_err());
    for wavelength in [1100.0, 1500.0, 1950.0] {
        let sample = field
            .prepare(wavelength)
            .unwrap()
            .optical_sample(Point3::new(0.5, 0.0, 0.0))
            .unwrap();
        near(
            sample.absorption_per_m(),
            0.6 + 0.0004 * (wavelength - 1000.0) + 0.4,
            2.0e-14,
        );
    }
    for wavelength in [999.0, 2001.0] {
        assert!(field.prepare(wavelength).is_err());
    }
    let prepared = field.prepare(1500.0).unwrap();
    for (x, alpha) in [(-0.5, 0.4), (1.5, 2.0)] {
        near(
            prepared.sample(Point3::new(x, 0.0, 0.0), 1500.0).unwrap().absorption_per_m(),
            alpha,
            2.0e-14,
        );
    }
    for x in [-0.5001, 1.5001] {
        assert_eq!(
            prepared.sample(Point3::new(x, 0.0, 0.0), 1500.0),
            Err(GradientIndexError::OutsideFieldDomain),
        );
    }
    let before = prepared.sample(origin(), 1500.0).unwrap();
    assert_eq!(
        propagate_with_control(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            1500.0,
            0.5,
            tight(),
            || false,
        ),
        Err(GradientIndexError::Cancelled),
    );
    assert_eq!(prepared.sample(origin(), 1500.0).unwrap(), before);
    let (card, original) = source();
    for (range, cancelled) in [([279.0, 360.0], false), ([280.0, 360.0], true)] {
        assert!(SpectralDielectricAbsorptionTemperaturePath::try_new(
            &card, original.clone(), "T", range, 128, || cancelled,
        ).is_err());
    }
}
