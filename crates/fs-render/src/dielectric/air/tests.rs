use super::*;
use crate::camera::{DispersiveThinLens, ThinLensElement};
use crate::charts::Ray;
use crate::dielectric::{
    BeerLambertAbsorption, DielectricEvent, DielectricGlass, DielectricSurface, GlassProvenance,
    sample_smooth_dielectric,
};
use crate::gradient_index::{GradientIndexOptions, HomogeneousIndexField, propagate};
use fs_geom::{Point3, Vec3};

const STRENGTHS: [f64; 3] = [1.03961212, 0.231792344, 1.01046945];
const RESONANCES: [f64; 3] = [0.00600069867, 0.0200179144, 103.560653];

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.17e} != {expected:.17e}");
}

fn air(temperature: f64, pressure: f64) -> CiddorAirIor {
    CiddorAirIor::try_new(CiddorAir::try_new(0.0, 450.0).unwrap(), temperature, pressure).unwrap()
}

fn relative(reference: CiddorAirIor, coordinate: AirRelativeWavelength) -> AirRelativeSellmeierIor {
    AirRelativeSellmeierIor::try_new(
        STRENGTHS, RESONANCES, (320.0, 1800.0), reference, coordinate,
    ).unwrap()
}

// Analytic oscillator oracle evaluated in the source coordinate, with the
// derivative written independently from the native group-index implementation.
fn source_reference(wavelength_nm: f64) -> (f64, f64) {
    let lambda = wavelength_nm * 0.001;
    let square = lambda * lambda;
    let mut n2 = 1.0;
    let mut derivative_n2 = 0.0;
    for (b, c) in STRENGTHS.into_iter().zip(RESONANCES) {
        n2 += b * square / (square - c);
        derivative_n2 -= 2.0 * b * lambda * c / (square - c).powi(2);
    }
    let n = n2.sqrt();
    (n, n - lambda * derivative_n2 / (2.0 * n))
}

#[test]
fn g1_original_air_state_reaches_native_phase_and_group_without_fitting() {
    let composition = CiddorAir::try_new(0.012, 650.0).unwrap();
    let law = CiddorAirIor::try_new(composition, 302.15, 95_000.0).unwrap();
    let native = DielectricIor::from(law);
    assert_eq!(native.wavelength_bounds_nm(), (300.0, 1700.0));
    assert!(native.is_dispersive());
    for wavelength in [300.0, 450.0, 632.991, 1550.0, 1700.0] {
        let original = composition.prepare(wavelength).unwrap().sample(302.15, 95_000.0).unwrap();
        assert_eq!(native.eval(wavelength), Ok(original.phase_index));
        assert_eq!(native.group_index(wavelength), Ok(original.group_index));
        assert!(native.group_index(wavelength).unwrap() > native.eval(wavelength).unwrap());
    }
}

#[test]
fn g1_absolute_catalog_phase_and_group_use_each_declared_wavelength_convention() {
    let index_air = air(293.15, 101_330.0);
    let coordinate_air = air(288.15, 101_325.0);
    for coordinate in [
        AirRelativeWavelength::Vacuum,
        AirRelativeWavelength::Air(index_air),
        AirRelativeWavelength::Air(coordinate_air),
    ] {
        let law = relative(index_air, coordinate);
        for lambda in [350.0, 550.0, 1000.0, 1600.0] {
            let a = index_air.eval(lambda).unwrap();
            let ga = index_air.group_index(lambda).unwrap();
            let (b, gb) = match coordinate {
                AirRelativeWavelength::Vacuum => (1.0, 1.0),
                AirRelativeWavelength::Air(air) => (
                    air.eval(lambda).unwrap(), air.group_index(lambda).unwrap(),
                ),
            };
            let (r, gr) = source_reference(lambda / b);
            near(law.eval(lambda).unwrap(), a * r, 1e-14);
            near(law.group_index(lambda).unwrap(), ga * r + (a * gb / b) * (gr - r), 1e-14);
            assert!((law.eval(lambda).unwrap() - r).abs() > 3e-4);
            // Independent five-point derivative of the full converted phase
            // detects both the reference-air dispersion and coordinate chain rule.
            let h = 0.1;
            let derivative = (law.eval(lambda - 2.0 * h).unwrap()
                - 8.0 * law.eval(lambda - h).unwrap()
                + 8.0 * law.eval(lambda + h).unwrap()
                - law.eval(lambda + 2.0 * h).unwrap()) / (12.0 * h);
            near(law.group_index(lambda).unwrap(), law.eval(lambda).unwrap() - lambda * derivative, 2e-9);
        }
    }
    let vacuum = relative(index_air, AirRelativeWavelength::Vacuum);
    let reference_coordinate = relative(index_air, AirRelativeWavelength::Air(index_air));
    let standard_coordinate = relative(index_air, AirRelativeWavelength::Air(coordinate_air));
    assert!((vacuum.eval(450.0).unwrap() - reference_coordinate.eval(450.0).unwrap()).abs() > 1e-5);
    assert!((reference_coordinate.eval(450.0).unwrap() - standard_coordinate.eval(450.0).unwrap()).abs() > 1e-8);
}

#[test]
fn g1_source_band_is_mapped_inward_and_never_extrapolated() {
    let reference = air(293.15, 101_330.0);
    let coordinate_air = air(288.15, 101_325.0);
    let law = AirRelativeSellmeierIor::try_new(
        STRENGTHS, RESONANCES, (380.0, 780.0), reference,
        AirRelativeWavelength::Air(coordinate_air),
    ).unwrap();
    let (lo, hi) = law.wavelength_bounds_nm();
    assert!(lo > 380.0 && hi > 780.0);
    assert!(lo / coordinate_air.eval(lo).unwrap() >= 380.0);
    assert!(hi / coordinate_air.eval(hi).unwrap() <= 780.0);
    assert!(lo.next_down() / coordinate_air.eval(lo.next_down()).unwrap() < 380.0);
    assert!(hi.next_up() / coordinate_air.eval(hi.next_up()).unwrap() > 780.0);
    for edge in [lo, hi] {
        assert!(law.eval(edge).is_ok());
        assert!(law.group_index(edge).is_ok());
    }
    for outside in [380.0, lo.next_down(), hi.next_up(), f64::NAN] {
        assert_eq!(law.eval(outside), Err(DielectricError::InvalidWavelength));
        assert_eq!(law.group_index(outside), Err(DielectricError::InvalidWavelength));
    }
    assert_eq!(relative(reference, AirRelativeWavelength::Vacuum).wavelength_bounds_nm(), (320.0, 1700.0));
    assert!(AirRelativeSellmeierIor::try_new(
        STRENGTHS, RESONANCES, (1900.0, 2100.0), reference,
        AirRelativeWavelength::Air(coordinate_air),
    ).is_err());
    for (temperature, pressure) in [(200.0, 101_325.0), (300.0, 0.0), (f64::NAN, 100_000.0)] {
        assert!(CiddorAirIor::try_new(reference.composition(), temperature, pressure).is_err());
    }
}

#[test]
fn g1_native_snell_and_lens_power_keep_reference_air_separate_from_scene_air() {
    let reference = air(293.15, 101_330.0);
    let law = relative(reference, AirRelativeWavelength::Vacuum);
    let lens = DispersiveThinLens::try_new_in_medium(
        550.0, &[ThinLensElement::try_new(law.into(), 20.0, -20.0).unwrap()], reference.into(),
    ).unwrap();
    let warm_scene_air = air(333.15, 80_000.0);
    let warm_lens = DispersiveThinLens::try_new_in_medium(
        550.0, &[ThinLensElement::try_new(law.into(), 20.0, -20.0).unwrap()], warm_scene_air.into(),
    ).unwrap();
    for wavelength in [450.0, 550.0, 1000.0, 1600.0] {
        let relative_index = source_reference(wavelength).0;
        let absolute_index = law.eval(wavelength).unwrap();
        let incident_index = reference.eval(wavelength).unwrap();
        let sample = sample_smooth_dielectric(
            Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.6, 0.0, 0.8),
            incident_index, absolute_index, 0.999,
        ).unwrap();
        assert_eq!(sample.event, DielectricEvent::Transmission);
        near(-sample.direction.x, 0.6 / relative_index, 2e-15);
        near(lens.power_per_m(wavelength).unwrap(), 40.0 * (relative_index - 1.0), 2e-13);
        near(warm_lens.power_per_m(wavelength).unwrap(),
            40.0 * (absolute_index / warm_scene_air.eval(wavelength).unwrap() - 1.0), 2e-13);
        assert!((warm_lens.power_per_m(wavelength).unwrap() - lens.power_per_m(wavelength).unwrap()).abs() > 1e-3);
    }
}

#[test]
fn g1_native_bulk_paths_retain_absolute_phase_group_and_explicit_beer_loss() {
    let reference = air(293.15, 101_330.0);
    let law = relative(reference, AirRelativeWavelength::Air(air(288.15, 101_325.0)));
    let glass = DielectricGlass::new(
        law, BeerLambertAbsorption::try_constant(0.07).unwrap(), GlassProvenance::Custom,
    );
    let field = HomogeneousIndexField::new(&glass);
    for wavelength in [450.0, 633.0, 1550.0] {
        let path = propagate(&field, Ray {
            origin: Point3::new(0.0, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0),
        }, Vec3::new(1.0, 0.0, 0.0), wavelength, 3.0, GradientIndexOptions::default()).unwrap();
        near(path.optical_path_length_m, 3.0 * law.eval(wavelength).unwrap(), 2e-15);
        near(path.group_optical_path_length_m.unwrap(), 3.0 * law.group_index(wavelength).unwrap(), 2e-15);
        near(path.power_transmittance(), (-0.21_f64).exp(), 3e-15);
        assert_eq!(path.field_evaluations, 1);
    }
}

#[test]
fn g3_native_identity_binds_both_atmospheres_original_coefficients_and_source_support() {
    let reference = air(293.15, 101_330.0);
    let coordinate_air = air(288.15, 101_325.0);
    let base = relative(reference, AirRelativeWavelength::Air(coordinate_air));
    let identity = |law: AirRelativeSellmeierIor| {
        crate::tracer::Material::Dielectric {
            glass: DielectricGlass::new(law, BeerLambertAbsorption::CLEAR, GlassProvenance::Custom),
            surface: DielectricSurface::SMOOTH,
        }.content_identity()
    };
    let same = relative(reference, AirRelativeWavelength::Air(coordinate_air));
    assert_eq!(identity(base), identity(same));
    for changed in [
        relative(air(294.15, 101_330.0), AirRelativeWavelength::Air(coordinate_air)),
        relative(reference, AirRelativeWavelength::Air(air(289.15, 101_325.0))),
        relative(reference, AirRelativeWavelength::Vacuum),
        AirRelativeSellmeierIor::try_new(STRENGTHS, RESONANCES, (330.0, 1800.0), reference,
            AirRelativeWavelength::Air(coordinate_air)).unwrap(),
        AirRelativeSellmeierIor::try_new([1.04, STRENGTHS[1], STRENGTHS[2]], RESONANCES,
            (320.0, 1800.0), reference, AirRelativeWavelength::Air(coordinate_air)).unwrap(),
    ] {
        assert_ne!(identity(base), identity(changed));
    }
    for changed in [
        air(294.15, 101_330.0), air(293.15, 101_325.0),
        CiddorAirIor::try_new(CiddorAir::try_new(0.01, 450.0).unwrap(), 293.15, 101_330.0).unwrap(),
        CiddorAirIor::try_new(CiddorAir::try_new(0.0, 600.0).unwrap(), 293.15, 101_330.0).unwrap(),
    ] {
        assert_ne!(reference.content_identity(), changed.content_identity());
    }
}
