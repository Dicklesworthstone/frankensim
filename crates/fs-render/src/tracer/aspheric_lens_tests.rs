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
