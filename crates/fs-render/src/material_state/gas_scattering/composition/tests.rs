use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::{AVOGADRO_PER_MOL, MolecularRayleigh, MolecularRayleighMixture};
use fs_qty::semantic::{Composition, CompositionBasis};

mod spectral_band;

fn model() -> MolecularRefractionModel {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optics = [0.00015, 0.0006]
        .map(|r| MolecularRayleigh::try_new(r, 2.5e25, 550.0, 0.03, ContentHash([1; 32])).unwrap());
    let source = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optics,
        None,
    )
    .unwrap();
    MolecularRefractionModel::from_mixture(
        &source,
        source.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap()
}

fn dispersive_model() -> MolecularRefractionModel {
    use fs_material::gas_scattering::composition::MolecularCauchyDispersion;
    let reference = model();
    let shapes: Vec<_> = reference
        .components()
        .iter()
        .enumerate()
        .map(|(i, component)| {
            MolecularCauchyDispersion::try_new(
                component.optical(),
                [2.0e-6, 8.0e-6][i],
                [1.0e-7, 3.0e-7][i],
                [380.0, 780.0],
                ContentHash([i as u8 + 3; 32]),
            )
            .unwrap()
        })
        .collect();
    reference.with_dispersion(&shapes).unwrap()
}

// Two opposite partial-density quadratics: TOTAL mass density is constant.
struct CompositionField {
    identity: ContentHash,
    count: usize,
    written: usize,
    basis: PartialDensityBasis,
    scale: [f64; 2],
    total: f64,
    gradient: [f64; 3],
    diagonal: [f64; 3],
    off_diagonal: [f64; 3],
    refusal: Option<GradientIndexError>,
}

impl CompositionField {
    fn new(model: &MolecularRefractionModel, basis: PartialDensityBasis) -> Self {
        Self {
            identity: model.content_identity(),
            count: 2,
            written: 2,
            basis,
            scale: core::array::from_fn(|i| match basis {
                PartialDensityBasis::MassKgPerM3 => 1.0,
                PartialDensityBasis::MoleculesPerM3 => {
                    AVOGADRO_PER_MOL / model.components()[i].molar_mass_kg_mol()
                }
            }),
            total: 2.0,
            gradient: [0.0; 3],
            diagonal: [0.0; 3],
            off_diagonal: [0.0; 3],
            refusal: None,
        }
    }

    fn jet(&self, p: Point3, i: usize) -> GasPartialDensityJet {
        let [xx, yy, zz] = self.diagonal;
        let [xy, xz, yz] = self.off_diagonal;
        let [gx, gy, gz] = self.gradient;
        let a = 0.5 * self.total
            + gx * p.x
            + gy * p.y
            + gz * p.z
            + 0.5 * (xx * p.x * p.x + yy * p.y * p.y + zz * p.z * p.z)
            + xy * p.x * p.y
            + xz * p.x * p.z
            + yz * p.y * p.z;
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let scale = self.scale[i];
        GasPartialDensityJet {
            sample: GasPartialDensitySample {
                density: (if i == 0 { a } else { self.total - a }) * scale,
                gradient_per_m: Vec3::new(
                    gx + xx * p.x + xy * p.y + xz * p.z,
                    gy + xy * p.x + yy * p.y + yz * p.z,
                    gz + xz * p.x + yz * p.y + zz * p.z,
                )
                .scale(sign * scale),
            },
            hessian_diagonal_per_m2: self.diagonal.map(|h| h * sign * scale),
            hessian_off_diagonal_per_m2: self.off_diagonal.map(|h| h * sign * scale),
        }
    }
}

impl GasCompositionGradientField for CompositionField {
    fn model_identity(&self) -> ContentHash {
        self.identity
    }
    fn component_count(&self) -> usize {
        self.count
    }
    fn density_basis(&self) -> PartialDensityBasis {
        self.basis
    }
    fn sample_partial_densities(
        &self,
        p: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError> {
        if let Some(error) = &self.refusal {
            return Err(error.clone());
        }
        for (i, row) in output.iter_mut().take(self.written).enumerate() {
            *row = self.jet(p, i).sample;
        }
        Ok(())
    }
}
impl GasCompositionSecondOrderField for CompositionField {
    fn sample_partial_density_jets(
        &self,
        p: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError> {
        if let Some(error) = &self.refusal {
            return Err(error.clone());
        }
        for (i, row) in output.iter_mut().take(self.written).enumerate() {
            *row = self.jet(p, i);
        }
        Ok(())
    }
}

fn bind<'a>(
    model: &'a MolecularRefractionModel,
    provider: &'a CompositionField,
) -> MaterialGasCompositionRefractiveField<'a, CompositionField> {
    MaterialGasCompositionRefractiveField::try_new(model, provider, ContentHash([2; 32])).unwrap()
}
fn close(got: f64, expected: f64, tolerance: f64) {
    assert!(
        (got - expected).abs() <= tolerance,
        "{got:.17e} != {expected:.17e}"
    );
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(1.0, 0.0, 0.0),
    }
}
fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.1,
        max_step_m: 0.5,
        relative_tolerance: 1.0e-11,
        ..Default::default()
    }
}

#[test]
fn g1_partial_density_gradient_and_hessian_match_spatial_differences_in_both_bases() {
    let model = model();
    let p = Point3::new(0.1, 0.3, -0.2);
    for basis in [
        PartialDensityBasis::MassKgPerM3,
        PartialDensityBasis::MoleculesPerM3,
    ] {
        for (diagonal, off_diagonal) in [
            ([0.0; 3], [0.0; 3]),
            ([0.04, 0.03, -0.01], [0.02, -0.03, 0.015]),
        ] {
            let mut provider = CompositionField::new(&model, basis);
            provider.gradient = [0.2, -0.1, 0.05];
            provider.diagonal = diagonal;
            provider.off_diagonal = off_diagonal;
            let field = bind(&model, &provider);
            let jet = field.sample_second_order(p, 550.0).unwrap();
            assert_eq!(jet.optical(), field.sample(p, 550.0).unwrap());
            let [xx, yy, zz] = jet.hessian().diagonal();
            let [xy, xz, yz] = jet.hessian().off_diagonal();
            let gradient = jet.optical().phase_gradient_per_m();
            for (axis, column) in [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
                .into_iter()
                .enumerate()
            {
                let mut a = [p.x, p.y, p.z];
                let mut b = a;
                let h = 0.01;
                a[axis] -= h;
                b[axis] += h;
                let a = field.sample(Point3::new(a[0], a[1], a[2]), 550.0).unwrap();
                let b = field.sample(Point3::new(b[0], b[1], b[2]), 550.0).unwrap();
                close(
                    (b.phase_index() - a.phase_index()) / (2.0 * h),
                    [gradient.x, gradient.y, gradient.z][axis],
                    5.0e-12,
                );
                let ga = a.phase_gradient_per_m();
                let gb = b.phase_gradient_per_m();
                for ((a, b), expected) in [ga.x, ga.y, ga.z]
                    .into_iter()
                    .zip([gb.x, gb.y, gb.z])
                    .zip(column)
                {
                    close(
                        (b - a) / (2.0 * h),
                        expected,
                        2.0e-8 * expected.abs() + 2.0e-15,
                    );
                }
            }
            if diagonal == [0.0; 3] {
                assert!(
                    xx > 0.0 && yy > 0.0 && zz > 0.0,
                    "nonlinear phase curvature must survive affine composition"
                );
            }
        }
    }
}

#[test]
fn g1_composition_alone_bends_native_rays_at_constant_total_density() {
    let model = model();
    let mut provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    provider.gradient = [0.0, 0.5, 0.0];
    let field = bind(&model, &provider);
    let result = propagate(
        &field,
        ray(),
        Vec3::new(0.0, 0.0, 1.0),
        550.0,
        20.0,
        options(),
    )
    .unwrap();
    assert!(result.ray.dir.y < -1.0e-4);
    let start = field.sample(ray().origin, 550.0).unwrap();
    let end = field.sample(result.ray.origin, 550.0).unwrap();
    close(
        end.phase_index() * result.ray.dir.x,
        start.phase_index(),
        2.0e-11,
    );
    let a = provider.jet(result.ray.origin, 0).sample.density;
    let b = provider.jet(result.ray.origin, 1).sample.density;
    close(a + b, provider.total, 1.0e-15);
    assert_eq!(result.power_transmittance(), 1.0);
    assert_eq!(
        result.group_optical_path_length_m,
        Some(result.optical_path_length_m)
    );
}

#[test]
fn g1_composition_curvature_focuses_native_bundles_and_matches_neighboring_rays() {
    let model = model();
    let basis = PartialDensityBasis::MassKgPerM3;
    let index = model.refractive_index(&[1.0, 1.0], basis, 550.0).unwrap();
    let delta_k =
        model.components()[0].coefficient(basis) - model.components()[1].coefficient(basis);
    let length = 40.0;
    for h in [-0.5, 0.5] {
        let mut provider = CompositionField::new(&model, basis);
        provider.diagonal = [0.0, h, h];
        let field = bind(&model, &provider);
        let result = propagate_bundle(
            &field,
            ray(),
            Vec3::new(0.0, 1.0, 0.0),
            550.0,
            length,
            BundleOptions {
                integration: options(),
                ..Default::default()
            },
        )
        .unwrap();
        // On axis n is constant. Independent exact paraxial equation:
        // X'' = q X, X(0)=0, X'(0)=1, q = n_L (K_a-K_b) h / n.
        let q = index.derivative_lorentz_lorenz * delta_k * h / index.phase_index;
        let rate = q.abs().sqrt();
        let expected = if q < 0.0 {
            (rate * length).sin() / rate
        } else {
            (rate * length).sinh() / rate
        };
        let projection = result
            .project_at_endpoint(ray().dir, ProjectionOptions::default())
            .unwrap();
        close(
            projection.area_per_solid_angle_m2(),
            expected * expected,
            2.0e-8,
        );
        close(
            projection.point_source_irradiance(7.0).unwrap(),
            7.0 / (expected * expected),
            1.0e-13,
        );
        assert_eq!(result.propagation().ray.dir, ray().dir);
        let epsilon = 1.0e-5;
        let neighbor = propagate(
            &field,
            Ray {
                origin: ray().origin,
                dir: Vec3::new(1.0, epsilon, 0.0),
            },
            Vec3::new(0.0, 0.0, 1.0),
            550.0,
            length,
            options(),
        )
        .unwrap();
        close(neighbor.ray.origin.y / epsilon, expected, 2.0e-6);
    }
}

#[test]
fn g0_composition_binding_and_queries_refuse_missing_nonphysical_or_foreign_data() {
    let model = model();
    let mut provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    assert!(matches!(
        MaterialGasCompositionRefractiveField::try_new(&model, &provider, ContentHash([0; 32])),
        Err(GasRefractiveFieldError::MissingSnapshotIdentity)
    ));
    provider.identity = ContentHash([9; 32]);
    assert!(matches!(
        MaterialGasCompositionRefractiveField::try_new(&model, &provider, ContentHash([2; 32])),
        Err(GasRefractiveFieldError::Gas(
            GasScatteringError::StateMismatch
        ))
    ));
    provider.identity = model.content_identity();
    provider.count = 0;
    assert!(matches!(
        MaterialGasCompositionRefractiveField::try_new(&model, &provider, ContentHash([2; 32])),
        Err(GasRefractiveFieldError::Gas(GasScatteringError::Components))
    ));
    provider.count = 2;
    provider.written = 1;
    assert_eq!(
        bind(&model, &provider).sample(ray().origin, 550.0),
        Err(GradientIndexError::InvalidSample)
    );
    assert_eq!(
        bind(&model, &provider).sample_second_order(ray().origin, 550.0),
        Err(GradientIndexError::InvalidSample)
    );
    provider.written = 2;
    provider.total = -1.0;
    assert_eq!(
        bind(&model, &provider).sample(ray().origin, 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    provider.total = 2.0;
    provider.gradient[0] = f64::NAN;
    assert_eq!(
        bind(&model, &provider).sample(ray().origin, 550.0),
        Err(GradientIndexError::InvalidSample)
    );
    provider.gradient[0] = 0.0;
    provider.off_diagonal[0] = f64::NAN;
    assert!(
        bind(&model, &provider)
            .sample_second_order(ray().origin, 550.0)
            .is_err()
    );
    provider.off_diagonal[0] = 0.0;
    provider.refusal = Some(GradientIndexError::Cancelled);
    assert_eq!(
        bind(&model, &provider).sample_second_order(ray().origin, 550.0),
        Err(GradientIndexError::Cancelled)
    );
    provider.refusal = None;
    provider.total = 0.0;
    let field = bind(&model, &provider);
    let jet = field.sample_second_order(ray().origin, 550.0).unwrap();
    assert_eq!(jet.optical().phase_index(), 1.0);
    assert_eq!(
        jet.optical().phase_gradient_per_m(),
        Vec3::new(0.0, 0.0, 0.0)
    );
    assert_eq!(jet.hessian(), IndexHessian::ZERO);
    assert_eq!(
        field.sample(Point3::new(f64::NAN, 0.0, 0.0), 550.0),
        Err(GradientIndexError::InvalidSample)
    );
    let result = propagate(
        &field,
        ray(),
        Vec3::new(0.0, 1.0, 0.0),
        550.0,
        1.0,
        options(),
    )
    .unwrap();
    close(result.ray.origin.x, 1.0, 1.0e-14);
    assert_eq!(result.ray.dir, ray().dir);
}

#[test]
fn g1_molecular_dispersion_bends_chromatically_and_tracks_true_group_delay() {
    let model = dispersive_model();
    for basis in [
        PartialDensityBasis::MassKgPerM3,
        PartialDensityBasis::MoleculesPerM3,
    ] {
        let mut provider = CompositionField::new(&model, basis);
        provider.gradient = [0.0, 0.5, 0.0];
        let field = bind(&model, &provider);
        let colors = [380.0, 780.0].map(|lambda| {
            let start = field.sample(ray().origin, lambda).unwrap();
            assert!(start.group_index().unwrap() > start.phase_index());
            let prepared = field.at_wavelength(lambda).unwrap();
            let result = propagate(
                &prepared,
                ray(),
                Vec3::new(0.0, 0.0, 1.0),
                lambda,
                20.0,
                options(),
            )
            .unwrap();
            assert!(result.group_optical_path_length_m.unwrap() > result.optical_path_length_m);
            let end = field.sample(result.ray.origin, lambda).unwrap();
            close(
                end.phase_index() * result.ray.dir.x,
                start.phase_index(),
                2.0e-11,
            );
            result
        });
        assert!(
            colors[0].ray.dir.y.abs() > 1.03 * colors[1].ray.dir.y.abs(),
            "blue composition refraction must respond to the sourced dispersion"
        );
        provider.gradient = [0.0; 3];
        let field = bind(&model, &provider);
        let index = field.sample(ray().origin, 500.0).unwrap();
        let result = propagate(
            &field,
            ray(),
            Vec3::new(0.0, 1.0, 0.0),
            500.0,
            37.0,
            options(),
        )
        .unwrap();
        close(
            result.optical_path_length_m,
            37.0 * index.phase_index(),
            2.0e-12,
        );
        close(
            result.group_optical_path_length_m.unwrap(),
            37.0 * index.group_index().unwrap(),
            2.0e-12,
        );
    }
}

#[test]
fn g1_spectral_composition_hessian_and_bundle_focusing_use_the_same_phase_law() {
    let model = dispersive_model();
    let basis = PartialDensityBasis::MassKgPerM3;
    let mut provider = CompositionField::new(&model, basis);
    provider.diagonal = [0.0, 0.5, 0.5];
    let field = bind(&model, &provider);
    for lambda in [380.0, 780.0] {
        let p = Point3::new(0.0, 0.2, -0.1);
        let jet = field.sample_second_order(p, lambda).unwrap();
        assert_eq!(jet.optical(), field.sample(p, lambda).unwrap());
        let [xx, yy, zz] = jet.hessian().diagonal();
        let [xy, xz, yz] = jet.hessian().off_diagonal();
        for (axis, expected) in [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
            .into_iter()
            .enumerate()
        {
            let mut a = [p.x, p.y, p.z];
            let mut b = a;
            let h = 0.01;
            a[axis] -= h;
            b[axis] += h;
            let ga = field
                .sample(Point3::new(a[0], a[1], a[2]), lambda)
                .unwrap()
                .phase_gradient_per_m();
            let gb = field
                .sample(Point3::new(b[0], b[1], b[2]), lambda)
                .unwrap()
                .phase_gradient_per_m();
            for ((a, b), expected) in [ga.x, ga.y, ga.z]
                .into_iter()
                .zip([gb.x, gb.y, gb.z])
                .zip(expected)
            {
                close(
                    (b - a) / (2.0 * h),
                    expected,
                    2.0e-8 * expected.abs() + 2.0e-15,
                );
            }
        }
        let spectral = model.at_wavelength(lambda, basis).unwrap();
        let index = spectral.refractive_index(&[1.0, 1.0]).unwrap();
        let k = spectral.coefficients();
        let q = index.derivative_lorentz_lorenz * (k[0].coefficient - k[1].coefficient) * 0.5
            / index.phase_index;
        assert!(q < 0.0);
        let rate = (-q).sqrt();
        let expected = (40.0 * rate).sin() / rate;
        let prepared = field.at_wavelength(lambda).unwrap();
        let bundle = propagate_bundle(
            &prepared,
            ray(),
            Vec3::new(0.0, 1.0, 0.0),
            lambda,
            40.0,
            BundleOptions {
                integration: options(),
                ..Default::default()
            },
        )
        .unwrap();
        let projection = bundle
            .project_at_endpoint(ray().dir, ProjectionOptions::default())
            .unwrap();
        close(
            projection.area_per_solid_angle_m2(),
            expected * expected,
            2.0e-8,
        );
        close(
            bundle.propagation().group_optical_path_length_m.unwrap(),
            40.0 * index.group_index,
            2.0e-11,
        );
    }
}

#[test]
fn g3_prepared_wavelength_preserves_native_curved_rays_and_focusing_exactly() {
    for model in [model(), dispersive_model()] {
        for basis in [
            PartialDensityBasis::MassKgPerM3,
            PartialDensityBasis::MoleculesPerM3,
        ] {
            let mut provider = CompositionField::new(&model, basis);
            provider.gradient = [0.0, 0.3, -0.1];
            provider.diagonal = [0.0, 0.05, 0.02];
            provider.off_diagonal = [0.01, -0.01, 0.015];
            let field = bind(&model, &provider);
            for wavelength in [380.0, 550.0, 780.0] {
                let prepared = field.at_wavelength(wavelength).unwrap();
                assert_eq!(prepared.wavelength_nm(), wavelength);
                assert_eq!(prepared.content_identity(), field.content_identity());
                let point = Point3::new(0.1, -0.2, 0.3);
                assert_eq!(
                    prepared.sample_optics(point).unwrap(),
                    field.sample_optics(point, wavelength).unwrap()
                );
                assert_eq!(
                    prepared.sample_optical_jet(point).unwrap(),
                    field.sample_optical_jet(point, wavelength).unwrap()
                );
                let axis = Vec3::new(0.0, 1.0, 0.0);
                let original = propagate(&field, ray(), axis, wavelength, 5.0, options()).unwrap();
                let reused = propagate(&prepared, ray(), axis, wavelength, 5.0, options()).unwrap();
                assert!(reused.field_evaluations > 6);
                assert_eq!(reused, original);
                let bundle_options = BundleOptions {
                    integration: options(),
                    ..Default::default()
                };
                let original =
                    propagate_bundle(&field, ray(), axis, wavelength, 5.0, bundle_options).unwrap();
                let reused =
                    propagate_bundle(&prepared, ray(), axis, wavelength, 5.0, bundle_options)
                        .unwrap();
                assert_eq!(reused, original);
                assert_eq!(
                    reused
                        .project_at_endpoint(ray().dir, ProjectionOptions::default())
                        .unwrap(),
                    original
                        .project_at_endpoint(ray().dir, ProjectionOptions::default())
                        .unwrap()
                );
            }
        }
    }
}

#[test]
fn g0_prepared_wavelength_retains_spectral_and_provider_refusals() {
    let model = dispersive_model();
    let mut provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    let field = bind(&model, &provider);
    for bad in [379.0, 781.0, f64::NAN, f64::INFINITY] {
        assert!(field.at_wavelength(bad).is_err());
    }
    let prepared = field.at_wavelength(550.0).unwrap();
    for wrong in [380.0, 550.0 + f64::EPSILON * 550.0, 780.0, f64::NAN] {
        let expected =
            GradientIndexError::Dielectric(crate::dielectric::DielectricError::InvalidWavelength);
        assert_eq!(prepared.sample(ray().origin, wrong), Err(expected));
        assert_eq!(
            prepared.sample_second_order(ray().origin, wrong),
            Err(expected)
        );
    }
    assert_eq!(
        prepared.sample(Point3::new(f64::NAN, 0.0, 0.0), 550.0),
        Err(GradientIndexError::InvalidSample)
    );
    provider.refusal = Some(GradientIndexError::OutsideFieldDomain);
    let field = bind(&model, &provider);
    // Existing unprepared entry points retain provider-before-spectrum errors.
    assert_eq!(
        field.sample(ray().origin, f64::NAN),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert_eq!(
        field.sample_second_order(ray().origin, f64::NAN),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    let prepared = field.at_wavelength(550.0).unwrap();
    assert_eq!(
        prepared.sample(ray().origin, 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    assert_eq!(
        prepared.sample_second_order(ray().origin, 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
}
