use super::*;
use crate::charts::Ray;
use crate::dielectric::{DielectricEvent, sample_smooth_dielectric};
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::TemperatureJet;
use fs_geom::Vec3;

#[derive(Clone, Copy)]
struct Thermal {
    base: f64,
    g: [f64; 3],
    h: [f64; 6],
}
impl Thermal {
    fn scalar(&self, [x, y, z]: [f64; 3]) -> f64 {
        self.base
            + self.g[0] * x
            + self.g[1] * y
            + self.g[2] * z
            + 0.5 * (self.h[0] * x * x + self.h[1] * y * y + self.h[2] * z * z)
            + self.h[3] * x * y
            + self.h[4] * x * z
            + self.h[5] * y * z
    }
}
impl TemperatureField for Thermal {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        if ![p.x, p.y, p.z].into_iter().all(f64::is_finite)
            || p.x.abs() > 0.05
            || p.y.abs() > 0.05
            || p.z.abs() > 1.0
        {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        Ok(self.scalar([p.x, p.y, p.z]))
    }
}
impl TemperatureGradientField for Thermal {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(self.temperature_jet(p)?.sample)
    }
    // Even a flat profile is bounded; retain the default nonuniform contract.
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
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn tight() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.01,
        max_step_m: 0.04,
        relative_tolerance: 2e-12,
        position_tolerance_m: 2e-12,
        direction_tolerance: 2e-12,
        optical_path_tolerance_m: 2e-12,
        optical_depth_tolerance: 2e-12,
        ..GradientIndexOptions::default()
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:.3e}"
    );
}

#[test]
fn g0_frozen_native_glass_retains_source_dispersion_snell_and_independent_beer_loss() {
    let absorption = BeerLambertAbsorption::try_constant(0.3).unwrap();
    for temperature in [30.0, 250.0, 300.0] {
        let source = Corning7980Silica::at_temperature(temperature).unwrap();
        let glass = corning7980_glass(temperature, absorption).unwrap();
        assert_eq!(glass.provenance(), GlassProvenance::CharmsCorning7980V1);
        assert_eq!(glass.wavelength_bounds_nm().unwrap(), (400.0, 2600.0));
        for wavelength in [400.0, 1550.0, 2600.0] {
            let sample = source.evaluate(wavelength).unwrap();
            assert_eq!(glass.ior().eval(wavelength).unwrap(), sample.phase_index);
            assert_eq!(
                glass.ior().group_index(wavelength).unwrap(),
                sample.group_index
            );
            let refracted = sample_smooth_dielectric(
                ray().dir,
                Vec3::new(0.6, 0.0, 0.8),
                1.0,
                sample.phase_index,
                0.9,
            )
            .unwrap();
            assert_eq!(refracted.event, DielectricEvent::Transmission);
            near(-refracted.direction.x * sample.phase_index, 0.6, 2e-15);
            near(
                refracted.radiance_weight,
                1.0 / sample.phase_index.powi(2),
                2e-16,
            );
            near(
                glass
                    .absorption_ref()
                    .transmittance(wavelength, 0.8)
                    .unwrap(),
                (-0.24_f64).exp(),
                2e-15,
            );
        }
    }
}

#[test]
fn g1_full_thermal_hessian_matches_independent_scalar_spatial_differences() {
    let thermal = Thermal {
        base: 220.0,
        g: [30.0, -20.0, 10.0],
        h: [300.0, -400.0, 200.0, 50.0, -75.0, 125.0],
    };
    let absorption = BeerLambertAbsorption::CLEAR;
    let field =
        MaterialSilicaRefractiveField::try_new(&thermal, &absorption, ContentHash([31; 32]))
            .unwrap();
    let p = [0.01, -0.02, 0.1];
    for wavelength in [450.0, 1550.0, 2500.0] {
        let prepared = field.prepare(wavelength).unwrap();
        let point = |p: [f64; 3]| Point3::new(p[0], p[1], p[2]);
        let jet = prepared.optical_jet(point(p)).unwrap();
        // Independent scalar composition does not consume either analytic jet.
        let phase = |p| {
            Corning7980Silica::at_temperature(thermal.scalar(p))
                .unwrap()
                .evaluate(wavelength)
                .unwrap()
                .phase_index
        };
        let g = jet.optical().phase_gradient_per_m();
        let h = 0.001;
        for (i, derivative) in [g.x, g.y, g.z].into_iter().enumerate() {
            let (mut a, mut b) = (p, p);
            a[i] += h;
            b[i] -= h;
            near(derivative, (phase(a) - phase(b)) / (2.0 * h), 3e-10);
            near(
                jet.hessian().diagonal()[i],
                (phase(a) - 2.0 * phase(p) + phase(b)) / (h * h),
                1.5e-9,
            );
        }
        for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            let (mut pp, mut pm, mut mp, mut mm) = (p, p, p, p);
            pp[i] += h;
            pp[j] += h;
            pm[i] += h;
            pm[j] -= h;
            mp[i] -= h;
            mp[j] += h;
            mm[i] -= h;
            mm[j] -= h;
            near(
                jet.hessian().off_diagonal()[k],
                (phase(pp) - phase(pm) - phase(mp) + phase(mm)) / (4.0 * h * h),
                1.5e-9,
            );
        }
        assert_eq!(jet.optical(), field.sample(point(p), wavelength).unwrap());
    }
}

#[test]
fn g1_original_transverse_temperature_bends_spectral_rays_and_preserves_optical_momentum() {
    let thermal = Thermal {
        base: 220.0,
        g: [400.0, 0.0, 0.0],
        h: [0.0; 6],
    };
    let absorption = BeerLambertAbsorption::try_constant(0.25).unwrap();
    let field =
        MaterialSilicaRefractiveField::try_new(&thermal, &absorption, ContentHash([32; 32]))
            .unwrap();
    let mut angles = [0.0; 2];
    for (i, wavelength) in [450.0, 1550.0].into_iter().enumerate() {
        let prepared = field.prepare(wavelength).unwrap();
        let path = propagate(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            0.8,
            tight(),
        )
        .unwrap();
        let initial = prepared.optical_sample(ray().origin).unwrap();
        let end = prepared.optical_sample(path.ray.origin).unwrap();
        // Translation symmetry in z conserves n*t_z along a stratified ray.
        near(
            end.phase_index() * path.ray.dir.z,
            initial.phase_index(),
            2e-11,
        );
        near(path.power_transmittance(), (-0.2_f64).exp(), 2e-13);
        assert!(path.ray.origin.x > 1e-4 && path.ray.dir.x > 0.0);
        near(path.ray.origin.y, 0.0, 2e-15);
        assert!(path.group_optical_path_length_m.unwrap() > path.optical_path_length_m);
        angles[i] = path.ray.dir.x;
    }
    assert!((angles[0] - angles[1]).abs() > 1e-6);
}

#[test]
fn g1_quadratic_thermal_lens_focuses_native_bundles_at_zero_axis_gradient() {
    let thermal = Thermal {
        base: 220.0,
        g: [0.0; 3],
        h: [-16000.0, -24000.0, 0.0, 0.0, 0.0, 0.0],
    };
    let absorption = BeerLambertAbsorption::try_constant(0.25).unwrap();
    let field =
        MaterialSilicaRefractiveField::try_new(&thermal, &absorption, ContentHash([33; 32]))
            .unwrap();
    let length = 0.8;
    for wavelength in [450.0, 1550.0] {
        let source = Corning7980Silica::at_temperature(220.0)
            .unwrap()
            .evaluate(wavelength)
            .unwrap();
        let prepared = field.prepare(wavelength).unwrap();
        let bundle = propagate_bundle(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            length,
            BundleOptions {
                integration: tight(),
                position_derivative_tolerance_m: 2e-12,
                direction_derivative_tolerance: 2e-12,
            },
        )
        .unwrap();
        let footprint = bundle
            .project_at_endpoint(ray().dir, ProjectionOptions::default())
            .unwrap();
        // On-axis Jacobi equation: X''=n_T*T_xx*X/n_axis, X(0)=0,X'(0)=1.
        let expected = [16000.0, 24000.0].map(|curvature| {
            let omega =
                (source.derivative_temperature_per_k * curvature / source.phase_index).sqrt();
            (omega * length).sin() / omega
        });
        near(
            footprint.area_per_solid_angle_m2(),
            expected[0] * expected[1],
            2e-10,
        );
        assert!(footprint.area_per_solid_angle_m2() < length * length - 0.001);
        near(
            bundle.propagation().optical_path_length_m,
            source.phase_index * length,
            2e-11,
        );
        near(
            bundle.propagation().group_optical_path_length_m.unwrap(),
            source.group_index * length,
            2e-11,
        );
        near(
            bundle.propagation().power_transmittance(),
            (-0.25 * length).exp(),
            2e-13,
        );
        assert_eq!(
            prepared
                .optical_sample(ray().origin)
                .unwrap()
                .phase_gradient_per_m(),
            Vec3::new(0.0, 0.0, 0.0)
        );
        assert_ne!(
            prepared.optical_jet(ray().origin).unwrap().hessian(),
            IndexHessian::ZERO
        );
    }
}

#[test]
fn g0_source_snapshot_domains_and_bounded_flat_fields_are_not_extended() {
    let thermal = Thermal {
        base: 220.0,
        g: [0.0; 3],
        h: [0.0; 6],
    };
    let absorption = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (1200.0, 0.1),
        (1800.0, 0.3),
    ])
    .unwrap();
    let snapshot = ContentHash([34; 32]);
    let field = MaterialSilicaRefractiveField::try_new(&thermal, &absorption, snapshot).unwrap();
    assert_eq!(field.law_identity(), Corning7980Silica::source_identity());
    assert_eq!(field.snapshot_identity(), snapshot);
    assert_eq!(field.wavelength_range_nm(), [1200.0, 1800.0]);
    assert_eq!(
        corning7980_glass(220.0, absorption)
            .unwrap()
            .wavelength_bounds_nm()
            .unwrap(),
        (1200.0, 1800.0)
    );
    let prepared = field.prepare(1500.0).unwrap();
    assert_eq!(prepared.snapshot_identity(), snapshot);
    assert!(!field.is_spatially_uniform() && !prepared.is_spatially_uniform());
    assert_eq!(
        prepared.sample(ray().origin, 1501.0),
        Err(DielectricError::InvalidWavelength.into())
    );
    assert_eq!(
        prepared.sample_second_order(ray().origin, 1501.0),
        Err(DielectricError::InvalidWavelength.into())
    );
    assert_eq!(
        propagate(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            1500.0,
            1.1,
            tight()
        ),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    for wavelength in [1199.0, 1801.0, f64::NAN, f64::INFINITY] {
        assert!(field.prepare(wavelength).is_err());
    }
    assert_eq!(
        prepared.sample(Point3::new(f64::NAN, 0.0, 0.0), 1500.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert!(matches!(
        MaterialSilicaRefractiveField::try_new(&thermal, &absorption, ContentHash([0; 32])),
        Err(SilicaFieldError::MissingSnapshotIdentity)
    ));
    let unsupported = Thermal {
        base: 310.0,
        ..thermal
    };
    let invalid =
        MaterialSilicaRefractiveField::try_new(&unsupported, &absorption, snapshot).unwrap();
    assert_eq!(
        invalid.sample(ray().origin, 1500.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert!(corning7980_glass(310.0, absorption).is_err());
    let disjoint = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (2800.0, 0.1),
        (3000.0, 0.2),
    ])
    .unwrap();
    assert!(matches!(
        MaterialSilicaRefractiveField::try_new(&thermal, &disjoint, snapshot),
        Err(SilicaFieldError::DisjointSpectralSupport)
    ));
    assert!(corning7980_glass(220.0, disjoint).is_err());
}
