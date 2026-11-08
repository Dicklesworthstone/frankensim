//! A changed off-axis prescription cannot reuse a progressive camera binding.
use super::*;
use crate::camera::{Aperture, CameraProjection, LensSurface, SequentialLens};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};

fn camera(conic: f64, coefficients: [f64; 7], explicit_asphere: bool) -> PhysicalCamera {
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.5, 0.0, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.0).unwrap(),
        GlassProvenance::Custom,
    );
    let rear = if explicit_asphere {
        LensSurface::try_new_aspheric(0.08, 20.0, conic, coefficients, 0.012, Some(glass))
    } else {
        LensSurface::try_new(0.08, 20.0, 0.012, Some(glass))
    }
    .unwrap();
    let lens = SequentialLens::try_new(
        550.0,
        &[
            rear,
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

fn identity(camera: &PhysicalCamera) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.render.aspheric-camera-test.v1");
    update_physical_camera_continuity(&mut hasher, camera);
    hasher.finalize()
}

#[test]
fn every_aspheric_parameter_changes_continuity_without_changing_paraxial_focus() {
    let sphere = camera(0.0, [0.0; 7], false);
    let baseline = identity(&sphere);
    assert_eq!(baseline, identity(&camera(0.0, [0.0; 7], true)));
    let conic = camera(-0.5, [0.0; 7], true);
    assert_eq!(sphere.focus_distance_m(), conic.focus_distance_m());
    assert_ne!(baseline, identity(&conic));
    for term in 0..7 {
        let mut coefficients = [0.0; 7];
        coefficients[term] = 1.0;
        let asphere = camera(0.0, coefficients, true);
        assert_eq!(sphere.focus_distance_m(), asphere.focus_distance_m());
        assert_ne!(baseline, identity(&asphere), "unbound asphere term {term}");
        assert_eq!(identity(&asphere), identity(&asphere.clone()));
    }
}

#[test]
fn lens_coating_changes_cannot_reuse_progressive_camera_accumulation() {
    use crate::conductor::{
        CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
    };
    use crate::dielectric::DielectricCoating;

    let bare = camera(0.0, [0.0; 7], false);
    let baseline = identity(&bare);
    let coated = |surface_index: usize, thickness: f64, source_byte: u8| {
        let source = ConductorSource::try_new(
            ContentHash([source_byte; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let layer =
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.0, thickness)
                .unwrap();
        let coating = DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(&[layer], 1, source).unwrap(),
        )
        .unwrap();
        let mut surfaces = bare.sequential_lens().unwrap().surfaces().to_vec();
        surfaces[surface_index] = surfaces[surface_index].with_coating(coating);
        bare.clone()
            .with_sequential_lens(SequentialLens::try_new(550.0, &surfaces).unwrap())
            .unwrap()
    };
    let reference = coated(0, 100e-9, 1);
    assert_eq!(bare.focus_distance_m(), reference.focus_distance_m());
    let reference_identity = identity(&reference);
    assert_ne!(baseline, reference_identity);
    assert_eq!(reference_identity, identity(&coated(0, 100e-9, 1)));
    for variant in [
        coated(1, 100e-9, 1),
        coated(0, 110e-9, 1),
        coated(0, 100e-9, 2),
    ] {
        assert_ne!(reference_identity, identity(&variant));
    }
}
