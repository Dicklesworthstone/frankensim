use super::*;
use crate::camera::{
    Aperture, CameraProjection, DefocusedDiffraction, DispersiveThinLens, ThinLensElement,
};
use crate::dielectric::{CauchyIor, DielectricIorSample, TabulatedIor};

fn camera(exterior: Option<DielectricIor>) -> PhysicalCamera {
    let element = ThinLensElement::try_new(
        CauchyIor::try_new(2.0, 0.0, 0.0).unwrap().into(),
        20.0,
        -20.0,
    )
    .unwrap();
    let lens = match exterior {
        Some(exterior) => DispersiveThinLens::try_new_in_medium(550.0, &[element], exterior),
        None => DispersiveThinLens::try_new(550.0, &[element]),
    }
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.002).unwrap(),
        2.0,
        Aperture::try_circular(0.0008).unwrap(),
    )
    .unwrap()
    .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(2.0), 0.01).unwrap())
    .unwrap()
    .with_dispersive_lens(lens)
    .unwrap()
}

#[test]
fn exterior_dispersion_rekeys_cinematic_camera_even_at_the_same_reference_focus() {
    let constant = CauchyIor::try_new(1.25, 0.0, 0.0).unwrap().into();
    let varying = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            let wave = 350.0 + i as f64 * 50.0;
            DielectricIorSample::try_new_in_band(wave, 1.25 + (wave - 550.0) * 1e-4).unwrap()
        }),
        50.0,
    )
    .unwrap()
    .into();
    let a = camera(Some(constant));
    let b = camera(Some(varying));
    assert_eq!(
        a.dispersive_lens().unwrap().reference_power_per_m(),
        b.dispersive_lens().unwrap().reference_power_per_m()
    );
    let identity = |camera: &PhysicalCamera| {
        let mut hasher = DomainHasher::new("frankensim.test.immersed-wave-camera");
        update_physical_camera_continuity(&mut hasher, camera);
        hasher.finalize()
    };
    assert_ne!(identity(&a), identity(&b));
    assert_ne!(identity(&a), identity(&camera(None)));
}
