//! Native propagation of explicitly supported broadband molecular sources.
use super::*;
use fs_material::gas_scattering::composition::MolecularCauchyDispersion;

fn broadband_model() -> MolecularRefractionModel {
    let reference = model();
    let sources: Vec<_> = reference
        .components()
        .iter()
        .enumerate()
        .map(|(i, component)| {
            MolecularCauchyDispersion::try_new_in_band(
                component.optical(),
                [2e-6, 8e-6][i],
                [1e-7, 3e-7][i],
                [300.0, 2500.0],
                ContentHash([i as u8 + 3; 32]),
            )
            .unwrap()
        })
        .collect();
    reference.with_dispersion(&sources).unwrap()
}

#[test]
fn infrared_phase_and_group_paths_use_the_original_molecular_source() {
    let model = broadband_model();
    let provider = CompositionField::new(&model, PartialDensityBasis::MassKgPerM3);
    let field = bind(&model, &provider);
    for lambda in [300.0, 1550.0, 2500.0] {
        let spectral = model
            .at_wavelength(lambda, PartialDensityBasis::MassKgPerM3)
            .unwrap();
        let expected = spectral.refractive_index(&[1.0, 1.0]).unwrap();
        let prepared = field.at_wavelength(lambda).unwrap();
        let result = propagate(
            &prepared,
            ray(),
            Vec3::new(0.0, 1.0, 0.0),
            lambda,
            37.0,
            options(),
        )
        .unwrap();
        close(
            result.optical_path_length_m,
            37.0 * expected.phase_index,
            2e-12,
        );
        close(
            result.group_optical_path_length_m.unwrap(),
            37.0 * expected.group_index,
            2e-12,
        );
        assert!(result.group_optical_path_length_m.unwrap() > result.optical_path_length_m);
        assert_eq!(
            prepared
                .sample(ray().origin, lambda)
                .unwrap()
                .absorption_per_m(),
            0.0
        );
        assert!(prepared.sample(ray().origin, 550.0).is_err());
    }
    assert!(field.at_wavelength(2500.0_f64.next_up()).is_err());
}

#[test]
fn broadband_composition_curvature_focuses_native_bundles_with_analytic_footprints() {
    let model = broadband_model();
    let basis = PartialDensityBasis::MassKgPerM3;
    let mut provider = CompositionField::new(&model, basis);
    provider.diagonal = [0.0, 0.5, 0.5];
    let field = bind(&model, &provider);
    let mut areas = Vec::new();
    for lambda in [300.0, 1550.0, 2500.0] {
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
        let area = projection.area_per_solid_angle_m2();
        close(area, expected * expected, 2e-8);
        close(
            bundle.propagation().group_optical_path_length_m.unwrap(),
            40.0 * index.group_index,
            2e-11,
        );
        areas.push(area);
    }
    assert!(areas[0] < areas[1] && areas[1] < areas[2]);
}
