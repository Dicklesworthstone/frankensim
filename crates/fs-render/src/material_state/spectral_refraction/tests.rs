//! Original tabulated thermo-optics through native rays, bundles and scenes.

use super::*;
use crate::charts::Ray;
use crate::dielectric::{BeerLambertAbsorption, GlassProvenance};
use crate::gradient_index::bundle::{BundleOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate, propagate_with_control};
use crate::material_state::thermo_optic::TemperatureJet;
use fs_evidence::ValidityDomain;
use fs_geom::Vec3;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, resolve_spectral_tabulated_dielectric_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}

/// Manufactured original claims, shared with the native element-traversal test.
pub(in crate::material_state) fn source_path(
    moving_wavelengths: bool,
    temperature_cusp: bool,
    range: [f64; 2],
) -> (SpectralDielectricTemperaturePath, DielectricGlass) {
    let mut claims = ClaimSet::new();
    for i in 0..9 {
        let wavelength_nm = 1000.0 + 125.0 * i as f64;
        let phase = 1.5 - 0.0001 * (wavelength_nm - 1000.0);
        let alpha = 0.6 + 0.0004 * (wavelength_nm - 1000.0);
        let wavelength = |t: f64| wavelength_nm * 1e-9
            * if moving_wavelengths { 1.0 + 0.0001 * (t - 300.0) } else { 1.0 };
        let phase_knots = if temperature_cusp {
            vec![(280.0, phase - 0.1), (310.0, phase + 0.05), (360.0, phase + 0.4)]
        } else {
            vec![(280.0, phase - 0.1), (360.0, phase + 0.3)]
        };
        for (name, quantity, knots) in [
            (
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                vec![(280.0, wavelength(280.0)), (360.0, wavelength(360.0))],
            ),
            (
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i],
                QuantitySpec::dimensional(Dims::NONE),
                phase_knots,
            ),
            (
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                vec![(280.0, alpha - 0.4), (360.0, alpha + 1.2)],
            ),
        ] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots, dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute_temperature(), 280.0, 360.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("manufactured original thermal spectrum: {name}"),
                    license: "internal-test-use".into(), artifact: None,
                },
            }).unwrap();
        }
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "original table thermal ray control".into(), phase: "solid".into(),
        process: "authored affine claims".into(), revision: 0,
    }, claims, Vec::new()).unwrap();
    let point = QueryPoint::new().with_quantity("T", absolute_temperature(), 300.0).unwrap();
    let source = resolve_spectral_tabulated_dielectric_state_point(
        &card, &point, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let reference = crate::material_state::bulk_optics::spectral_tabulated(&source, 200.0).unwrap();
    let path = SpectralDielectricTemperaturePath::try_new(
        &card, source, "T", range, 128, || false,
    ).unwrap();
    (path, reference)
}

#[derive(Clone, Copy)]
struct Thermal {
    base: f64,
    g: [f64; 3],
    h: [f64; 6],
}
impl TemperatureField for Thermal {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        if [p.x, p.y, p.z].into_iter().any(|v| !v.is_finite() || v.abs() > 2.0) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let [xx, yy, zz, xy, xz, yz] = self.h;
        Ok(self.base + self.g[0] * p.x + self.g[1] * p.y + self.g[2] * p.z
            + 0.5 * (xx * p.x * p.x + yy * p.y * p.y + zz * p.z * p.z)
            + xy * p.x * p.y + xz * p.x * p.z + yz * p.y * p.z)
    }
}
impl TemperatureGradientField for Thermal {
    fn temperature_gradient(&self, p: Point3) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(self.temperature_jet(p)?.sample)
    }
}
impl TemperatureSecondOrderField for Thermal {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        let [xx, yy, zz, xy, xz, yz] = self.h;
        Ok(TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: self.temperature_k(p)?,
                gradient_k_per_m: Vec3::new(
                    self.g[0] + xx * p.x + xy * p.y + xz * p.z,
                    self.g[1] + xy * p.x + yy * p.y + yz * p.z,
                    self.g[2] + xz * p.x + yz * p.y + zz * p.z,
                ),
            },
            diagonal_k_per_m2: [xx, yy, zz],
            off_diagonal_k_per_m2: [xy, xz, yz],
        })
    }
}
fn ray() -> Ray {
    Ray { origin: Point3::new(0.0, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0) }
}
fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.01, max_step_m: 0.04,
        relative_tolerance: 2e-12, position_tolerance_m: 2e-12,
        direction_tolerance: 2e-12, optical_path_tolerance_m: 2e-12,
        optical_depth_tolerance: 2e-12,
        ..GradientIndexOptions::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:.3e}");
}

#[test]
fn original_phase_group_and_current_loss_follow_the_same_curved_ray_and_bundle() {
    let (source, reference) = source_path(false, false, [280.0, 360.0]);
    let thermal = Thermal { base: 300.0, g: [40.0, 0.0, 0.0], h: [0.0; 6] };
    let snapshot = ContentHash([91; 32]);
    let field = MaterialSpectralRefractiveField::try_new(&reference, &source, &thermal, snapshot).unwrap();
    let wavelength = 1437.5;
    let prepared = field.prepare(wavelength).unwrap();
    assert_eq!(prepared.snapshot_identity(), snapshot);
    assert_eq!(field.source().identity(), source.identity());
    assert!(!field.is_spatially_uniform() && !prepared.is_spatially_uniform());
    let length: f64 = 0.8;
    let n0: f64 = 1.45625;
    let g: f64 = 0.2;
    // Analytic affine-index ray, derived from n*t_z conservation and arc length.
    let n_end = n0.hypot(g * length);
    let z = n0 / g * (g * length / n0).asinh();
    let x = (n_end - n0) / g;
    let phase = 0.5 * (length * n_end + n0 * z);
    let tau = 0.775 * length + 4.0 * (phase - n0 * length);
    let group = phase + 0.0001 * wavelength * length;
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let path = propagate(&prepared, ray(), axis, wavelength, length, options()).unwrap();
    let bundle = propagate_bundle(&prepared, ray(), axis, wavelength, length, BundleOptions {
        integration: options(), position_derivative_tolerance_m: 2e-12,
        direction_derivative_tolerance: 2e-12,
    }).unwrap();
    for result in [&path, bundle.propagation()] {
        near(result.ray.origin.x, x, 4e-11);
        near(result.ray.origin.y, 0.0, 2e-14);
        near(result.ray.origin.z, z, 4e-11);
        near(result.ray.dir.x, g * length / n_end, 4e-11);
        near(result.ray.dir.z, n0 / n_end, 4e-11);
        near(result.optical_path_length_m, phase, 4e-11);
        near(result.group_optical_path_length_m.unwrap(), group, 4e-11);
        near(result.optical_depth, tau, 8e-11);
        near(result.power_transmittance(), (-tau).exp(), 8e-11);
        assert!(result.ray.origin.x > 0.04);
    }
    #[cfg(feature = "photometry")]
    {
        use crate::photometry::scene::{NativeMaterialFieldBinding, NativeOpticalMedium};
        let medium = NativeOpticalMedium::Ambient;
        let bindings = [
            NativeMaterialFieldBinding::from_spectral_field(medium, &field),
            NativeMaterialFieldBinding::from_differentiable_spectral_field(medium, &field),
            NativeMaterialFieldBinding::from_prepared_spectral_field(medium, &prepared),
            NativeMaterialFieldBinding::from_prepared_differentiable_spectral_field(medium, &prepared),
        ];
        for binding in &bindings {
            assert_eq!(binding.source_snapshot_identity(), snapshot);
            assert_eq!(binding.field_identity(), bindings[0].field_identity());
        }
        let later = MaterialSpectralRefractiveField::try_new(
            &reference, &source, &thermal, ContentHash([92; 32]),
        ).unwrap();
        assert_ne!(bindings[0].field_identity(),
            NativeMaterialFieldBinding::from_spectral_field(medium, &later).field_identity());
    }
}

#[test]
fn moving_original_wavelengths_retain_the_full_nonlinear_spatial_chain_rule() {
    let (source, reference) = source_path(true, false, [280.0, 360.0]);
    let thermal = Thermal {
        base: 300.0, g: [3.0, -2.0, 4.0], h: [0.4, -0.2, 0.6, 0.1, -0.15, 0.25],
    };
    let field = MaterialSpectralRefractiveField::try_new(
        &reference, &source, &thermal, ContentHash([93; 32]),
    ).unwrap();
    let wavelength = 1437.5;
    let prepared = field.prepare(wavelength).unwrap();
    let jet = prepared.optical_jet(ray().origin).unwrap();
    // Independent scalar law: n=1.6+.005*dT-.0001*lambda/(1+.0001*dT).
    let first = 0.005 + 1e-8 * wavelength;
    let second = -2e-12 * wavelength;
    let g = jet.optical().phase_gradient_per_m();
    for (actual, expected) in [g.x, g.y, g.z].into_iter().zip(thermal.g.map(|v| first * v)) {
        near(actual, expected, 3e-15);
    }
    for i in 0..3 {
        near(jet.hessian().diagonal()[i],
            first * thermal.h[i] + second * thermal.g[i] * thermal.g[i], 3e-15);
    }
    for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
        near(jet.hessian().off_diagonal()[k],
            first * thermal.h[3 + k] + second * thermal.g[i] * thermal.g[j], 3e-15);
    }
    near(jet.optical().phase_index(), 1.45625, 2e-15);
    near(jet.optical().group_index().unwrap(), 1.6, 5e-15);
    near(jet.optical().absorption_per_m(), 0.775, 3e-15);
    assert_eq!(jet.optical(), field.sample(ray().origin, wavelength).unwrap());
    let stationary = Thermal { g: [0.0; 3], ..thermal };
    let stationary_field = MaterialSpectralRefractiveField::try_new(
        &reference, &source, &stationary, ContentHash([94; 32]),
    ).unwrap();
    let stationary_jet = stationary_field.sample_second_order(ray().origin, wavelength).unwrap();
    assert_eq!(stationary_jet.optical().phase_gradient_per_m(), Vec3::new(0.0, 0.0, 0.0));
    assert_ne!(stationary_jet.hessian(), IndexHessian::ZERO);
    for i in 0..3 {
        near(stationary_jet.hessian().diagonal()[i], first * stationary.h[i], 3e-15);
    }
}

#[test]
fn source_cusps_reference_policy_support_and_cancellation_remain_explicit() {
    let thermal = Thermal { base: 300.0, g: [0.0; 3], h: [0.0; 6] };
    let snapshot = ContentHash([95; 32]);
    let (cusp, reference) = source_path(false, true, [280.0, 360.0]);
    let field = MaterialSpectralRefractiveField::try_new(&reference, &cusp, &thermal, snapshot).unwrap();
    assert!(matches!(field.prepare(1437.5),
        Err(SpectralRefractiveFieldError::Source(SpectralDielectricTemperatureError::UndefinedDerivative))));
    let (cut, reference) = source_path(false, true, [280.0, 305.0]);
    let field = MaterialSpectralRefractiveField::try_new(&reference, &cut, &thermal, snapshot).unwrap();
    let prepared = field.prepare(1437.5).unwrap();
    assert_eq!(prepared.sample(ray().origin, 1438.0),
        Err(DielectricError::InvalidWavelength.into()));
    assert_eq!(prepared.sample_second_order(ray().origin, 1438.0),
        Err(DielectricError::InvalidWavelength.into()));
    for wavelength in [999.0, 2001.0, f64::NAN, f64::INFINITY] {
        assert!(field.prepare(wavelength).is_err());
    }
    assert_eq!(prepared.sample(Point3::new(f64::NAN, 0.0, 0.0), 1437.5),
        Err(GradientIndexError::OutsideFieldDomain));
    let too_hot = Thermal { base: 306.0, ..thermal };
    let outside = MaterialSpectralRefractiveField::try_new(&reference, &cut, &too_hot, snapshot).unwrap();
    assert_eq!(outside.sample(ray().origin, 1437.5), Err(GradientIndexError::OutsideFieldDomain));
    let wrong_loss = DielectricGlass::with_ior(
        *reference.ior(), BeerLambertAbsorption::CLEAR, GlassProvenance::Custom,
    );
    assert!(matches!(MaterialSpectralRefractiveField::try_new(&wrong_loss, &cut, &thermal, snapshot),
        Err(SpectralRefractiveFieldError::ReferenceSpectrumMismatch)));
    assert!(matches!(MaterialSpectralRefractiveField::try_new(
        &reference, &cut, &thermal, ContentHash([0; 32])),
        Err(SpectralRefractiveFieldError::MissingSnapshotIdentity)));
    let before = cut.identity();
    assert_eq!(propagate_with_control(
        &prepared, ray(), Vec3::new(1.0, 0.0, 0.0), 1437.5, 0.8, options(), || false,
    ), Err(GradientIndexError::Cancelled));
    assert_eq!(cut.identity(), before);
    // The reference grid meets this ceiling; its original hot grid exceeds it.
    let (moving, _) = source_path(true, false, [280.0, 360.0]);
    let tight_reference = crate::material_state::bulk_optics::spectral_tabulated(
        moving.source(), 125.25,
    ).unwrap();
    assert!(matches!(MaterialSpectralRefractiveField::try_new(
        &tight_reference, &moving, &thermal, snapshot),
        Err(SpectralRefractiveFieldError::Optics(DielectricError::InvalidTable))));
}
