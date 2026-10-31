use super::*;
use crate::charts::Ray;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::HomogeneousIndexField;
use crate::gradient_index::bundle::{
    BundleOptions, IndexHessian, QuadraticIndexField, propagate_bundle,
};
use crate::instances::RigidTransform;
use crate::photometry::source::{AngularDistribution, SourceSpectrum};
use fs_geom::Aabb;

fn fixture(position: [f64; 3]) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        1200.0,
        0.75,
        AngularDistribution::cosine_power(2.0, core::f64::consts::FRAC_PI_2).unwrap(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn glass(ior: CauchyIor, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-10),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_calibrated_bundles_recover_oblique_inverse_square_lux_and_absorption() {
    let source = fixture([0.0; 3]);
    let material = glass(CauchyIor::try_constant(1.5).unwrap(), 0.1);
    let field = HomogeneousIndexField::new(&material);
    let receiver = PlanarIrradianceReceiver::try_new(
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.6, 0.0, -0.8),
        1e-10,
    )
    .unwrap();
    let mut meter = BundleIrradianceMeter::new(&source, receiver, ProjectionOptions::default());
    for bin in 0..SPECTRAL_BINS {
        let bundle = propagate_bundle(
            &field,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            wavelength_nm(bin),
            2.0,
            BundleOptions::default(),
        )
        .unwrap();
        meter.record(bin, &bundle).unwrap();
    }
    let reading = meter.reading().unwrap();
    let intensity_cd = 900.0 * 3.0 / (2.0 * core::f64::consts::PI);
    let geometric_transmission = 0.8 / 4.0 * (-0.2_f64).exp();
    close(
        reading.illuminance_lux,
        intensity_cd * geometric_transmission,
        3e-13,
    );
    let radiant_per_lumen: f64 = source
        .spectrum()
        .bins_w_per_lumen_nm()
        .iter()
        .map(|v| v * BIN_WIDTH_NM)
        .sum();
    close(
        reading.radiant_irradiance_w_m2,
        intensity_cd * radiant_per_lumen * geometric_transmission,
        3e-13,
    );
}

#[test]
fn g1_calibrated_chromatic_grin_focusing_matches_independent_jacobi_solution() {
    let source = fixture([0.0; 3]);
    let material = glass(CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(), 0.15);
    let strength = 1.5_f64;
    let field = QuadraticIndexField::try_new(
        &material,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        IndexHessian::try_new([-strength, -strength, 0.0], [0.0; 3]).unwrap(),
        Aabb::new(Point3::new(-0.2, -0.2, -0.1), Point3::new(0.2, 0.2, 2.0)),
    )
    .unwrap();
    let distance = 1.4_f64;
    let receiver = PlanarIrradianceReceiver::try_new(
        Point3::new(0.0, 0.0, distance),
        Vec3::new(0.0, 0.0, -1.0),
        1e-9,
    )
    .unwrap();
    let mut meter = BundleIrradianceMeter::new(&source, receiver, ProjectionOptions::default());
    let mut expected_w_m2 = 0.0;
    let mut expected_lux = 0.0;
    for bin in 0..SPECTRAL_BINS {
        let wavelength = wavelength_nm(bin);
        let bundle = propagate_bundle(
            &field,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            distance,
            BundleOptions::default(),
        )
        .unwrap();
        meter.record(bin, &bundle).unwrap();
        // On the axis, the exact first variation satisfies X''=-(k/n0)X,
        // X(0)=0, X'(0)=1. Both transverse columns are sin(omega*s)/omega.
        // This reference uses the spectral law and oscillator closed form,
        // without the production propagator or its projection calculation.
        let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
        let omega = (strength / index).sqrt();
        let area = ((omega * distance).sin() / omega).powi(2);
        let spectral = source
            .spectral_intensity_w_sr_nm(ray().dir, wavelength)
            .unwrap()
            * (-0.15 * distance).exp()
            / area;
        close(
            meter.spectral_irradiance_w_m2_nm(bin).unwrap(),
            spectral,
            2e-7,
        );
        expected_w_m2 += spectral * BIN_WIDTH_NM;
        expected_lux += spectral * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
    }
    let reading = meter.reading().unwrap();
    close(reading.radiant_irradiance_w_m2, expected_w_m2, 2e-7);
    close(reading.illuminance_lux, expected_lux, 2e-7);
    let homogeneous_lux =
        source.intensity_cd(ray().dir).unwrap() * (-0.15 * distance).exp() / (distance * distance);
    assert!(reading.illuminance_lux > 1.8 * homogeneous_lux);
    assert!(
        meter.spectral_irradiance_w_m2_nm(0).unwrap()
            < meter
                .spectral_irradiance_w_m2_nm(SPECTRAL_BINS - 1)
                .unwrap()
    );
}

#[test]
fn g0_source_receiver_and_wavelength_binding_refuse_without_changing_meter() {
    let source = fixture([0.0; 3]);
    let material = glass(CauchyIor::try_constant(1.5).unwrap(), 0.0);
    let field = HomogeneousIndexField::new(&material);
    let receiver = PlanarIrradianceReceiver::try_new(
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 0.0, -1.0),
        1e-9,
    )
    .unwrap();
    let projection = ProjectionOptions::default();
    let bundle = propagate_bundle(
        &field,
        ray(),
        Vec3::new(1.0, 0.0, 0.0),
        wavelength_nm(0),
        2.0,
        BundleOptions::default(),
    )
    .unwrap();
    assert_eq!(
        fixture([1.0, 0.0, 0.0]).spectral_irradiance_from_bundle(&bundle, receiver, projection),
        Err(BundlePhotometryError::SourceOriginMismatch)
    );
    let displaced =
        PlanarIrradianceReceiver::try_new(Point3::new(0.1, 0.0, 2.0), receiver.normal(), 1e-9)
            .unwrap();
    assert_eq!(
        source.spectral_irradiance_from_bundle(&bundle, displaced, projection),
        Err(BundlePhotometryError::ReceiverPositionMismatch)
    );
    let reversed =
        PlanarIrradianceReceiver::try_new(receiver.position_m(), Vec3::new(0.0, 0.0, 1.0), 1e-9)
            .unwrap();
    assert_eq!(
        source.spectral_irradiance_from_bundle(&bundle, reversed, projection),
        Err(BundlePhotometryError::BackFacingReceiver)
    );
    let mut meter = BundleIrradianceMeter::new(&source, receiver, projection);
    assert_eq!(
        meter.record(1, &bundle),
        Err(BundlePhotometryError::WavelengthMismatch)
    );
    assert_eq!(
        meter.record(SPECTRAL_BINS, &bundle),
        Err(BundlePhotometryError::InvalidBin)
    );
    assert_eq!(meter.spectral_irradiance_w_m2_nm(0), None);
    assert_eq!(
        meter.reading(),
        Err(BundlePhotometryError::IncompleteSpectrum)
    );
    meter.record(0, &bundle).unwrap();
    let retained = meter.spectral_irradiance_w_m2_nm(0);
    assert_eq!(
        meter.record(0, &bundle),
        Err(BundlePhotometryError::DuplicateBin)
    );
    assert_eq!(meter.spectral_irradiance_w_m2_nm(0), retained);
    assert_eq!(
        meter.reading(),
        Err(BundlePhotometryError::IncompleteSpectrum)
    );
}
