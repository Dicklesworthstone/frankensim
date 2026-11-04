use super::*;
use crate::instances::{InstanceScene, RigidTransform};
use crate::polarization::dielectric::scene::anisotropic_tests::{
    binding, close, glass, input, plane, surface, with_cx,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};

#[test]
fn reciprocal_camera_receivers_and_recorded_pdfs_bind_the_same_body_axis() {
    with_cx(|cx| {
        let q = core::f64::consts::FRAC_1_SQRT_2;
        let placement = RigidTransform::try_new([0.0, 0.0, q, q], [0.0; 3]).unwrap();
        let geometry = InstanceScene::try_new(vec![plane(placement)]).unwrap();
        let media = [glass(1.0), glass(1.5)];
        let media = [
            SceneMedium::Homogeneous(&media[0]),
            SceneMedium::Homogeneous(&media[1]),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![binding()], 100_000).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let incoming = input(placement).advance(1.0).unwrap();
        let path = scene
            .trace_with_samples(
                cx,
                input(placement),
                0,
                CurvedTraceSettings {
                    max_surface_events: 1,
                    max_distance_m: 3.0,
                    ..CurvedTraceSettings::default()
                },
                |_| [0.04, 0.3, 0.0],
            )
            .unwrap();
        let contact = *path.last_contact().unwrap();
        let receiver = SurfaceReceiver::at_contact(&camera, contact, incoming, 0)
            .unwrap()
            .unwrap();
        let axis = placement.transform_vector(Vec3::new(1.0, 0.0, 0.0));
        let world_surface = surface().with_tangent_axis(axis).unwrap();
        for event in [0.0, 0.999] {
            let sample = incoming
                .sample_interface(
                    contact.witness.normal,
                    glass(1.5),
                    world_surface,
                    [0.04, 0.3, event],
                )
                .unwrap()
                .unwrap();
            let source_medium =
                usize::from(sample.event == crate::dielectric::DielectricEvent::Transmission);
            let frame = PolarizationFrame::try_new(
                sample.beam.frame().direction().scale(-1.0),
                Vec3::new(0.0, 1.0, 0.0),
            )
            .unwrap();
            for state in [
                Stokes::unpolarized(1.0),
                Stokes {
                    i: 1.0,
                    q: 0.2,
                    u: 0.4,
                    v: 0.6,
                },
            ] {
                let actual = receiver
                    .evaluate(
                        &camera,
                        contact,
                        incoming.ray().origin,
                        source_medium,
                        frame,
                        state,
                        0,
                    )
                    .unwrap();
                let expected = incoming
                    .evaluate_dielectric_connection(
                        contact.witness.normal,
                        glass(1.5),
                        world_surface,
                        frame,
                        state,
                    )
                    .unwrap();
                close(actual.factor, expected.factor);
                close(actual.pdf, expected.pdf);
                close(
                    receiver
                        .sampled_pdf(&camera, sample.beam.frame().direction(), source_medium, 0)
                        .unwrap(),
                    sample.pdf,
                );
            }
        }
    });
}
