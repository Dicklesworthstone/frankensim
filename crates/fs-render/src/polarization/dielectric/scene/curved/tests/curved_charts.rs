use super::*;
use fs_geom::fixtures::SphereChart;

fn sphere() -> GeometryInstance {
    GeometryInstance::try_new(
        1,
        hash_domain("curved-scene-test", b"native-unit-ball"),
        SharedGeometry::chart(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }),
        RigidTransform::identity(),
    )
    .unwrap()
}

#[derive(Clone, Copy)]
struct BallLensReference {
    detector_x: f64,
    direction: Vec3,
    outside_length: f64,
    inside_length: f64,
    s_transmission: f64,
    p_transmission: f64,
}

fn ball_lens_reference(index: f64) -> BallLensReference {
    // Exact meridional ball-lens geometry, R=1, launch x=.35,z=-2 along +z.
    // Entry incidence i=asin(.35), refraction t=asin(.35/n). The chord has
    // length 2*cos(t), exit normal angle 2t-i, outgoing angle -2(i-t).
    let i = 0.35_f64.asin();
    let t = (0.35 / index).asin();
    let exit_x = (2.0 * t - i).sin();
    let exit_z = (2.0 * t - i).cos();
    let direction = Vec3::new((-2.0 * (i - t)).sin(), 0.0, (2.0 * (i - t)).cos());
    let after = (2.0 - exit_z) / direction.z;
    let rs = ((i.cos() - index * t.cos()) / (i.cos() + index * t.cos())).powi(2);
    let rp = ((index * i.cos() - t.cos()) / (index * i.cos() + t.cos())).powi(2);
    BallLensReference {
        detector_x: exit_x + after * direction.x,
        direction,
        outside_length: 2.0 - i.cos() + after,
        inside_length: 2.0 * t.cos(),
        s_transmission: 1.0 - rs,
        p_transmission: 1.0 - rp,
    }
}

#[test]
fn g1_native_ball_lens_dispersion_reaches_finite_receiver_with_full_optical_lengths() {
    with_cx(|cx, _| {
        let wavelengths = [400.0_f64, 700.0_f64];
        let indices = wavelengths.map(|w| 1.45 + 0.012 / (w * 1e-3).powi(2));
        let reference = indices.map(ball_lens_reference);
        let receiver_x = 0.5 * (reference[0].detector_x + reference[1].detector_x);
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, 2.0, receiver_x, 0.1)]).unwrap();
        let lens = DielectricGlass::new(
            CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.17).unwrap(),
            GlassProvenance::Custom,
        );
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let mut limits = settings();
        limits.max_distance_m = 6.0;
        let mut positions = [0.0; 2];
        for (lane, wavelength) in wavelengths.into_iter().enumerate() {
            let input = scene
                .launch_beam(
                    Point3::new(0.35, 0.0, -2.0),
                    frame(),
                    wavelength,
                    1,
                    Stokes::unpolarized(1.0),
                )
                .unwrap();
            let mut interfaces = 0;
            let result = scene
                .trace_with_samples(cx, input, 1, limits, |_| {
                    interfaces += 1;
                    [0.2, 0.7, 0.999_999]
                })
                .unwrap();
            assert_eq!(interfaces, 2, "both entry and exit must refract");
            assert_eq!(result.termination(), PathTermination::Detected);
            assert_eq!(result.surface_events(), 3);
            assert_eq!(result.medium_index(), 1);
            assert_eq!(result.last_contact().unwrap().object_id, 9);
            assert_eq!(result.beam().medium(), BeamMedium::Homogeneous(air));
            let expected = reference[lane];
            positions[lane] = result.beam().ray().origin.x;
            assert!((positions[lane] - expected.detector_x).abs() < 3e-6);
            assert!((result.beam().ray().origin.z - 2.0).abs() < 2e-6);
            assert!((result.beam().ray().dir.x - expected.direction.x).abs() < 2e-6);
            assert!((result.beam().ray().dir.z - expected.direction.z).abs() < 2e-6);
            assert!(
                (result.traveled_m() - expected.outside_length - expected.inside_length).abs()
                    < 4e-6
            );
            assert!(
                (result.optical_path_length_m()
                    - expected.outside_length
                    - indices[lane] * expected.inside_length)
                    .abs()
                    < 5e-6
            );
            let group_index = 1.45 + 3.0 * 0.012 / (wavelength * 1e-3).powi(2);
            assert!(
                (result.group_optical_path_length_m().unwrap()
                    - expected.outside_length
                    - group_index * expected.inside_length)
                    .abs()
                    < 5e-6
            );
            let depth = 0.17 * expected.inside_length;
            assert!((result.optical_depth() - depth).abs() < 1e-6);
            // Each transmission uses the native unpolarized Fresnel proposal.
            // Retain the induced polarization across the second reciprocal face.
            let ts = expected.s_transmission;
            let tp = expected.p_transmission;
            let probability = 0.5 * (ts + tp);
            let scale = (-depth).exp() / (probability * probability);
            assert!((result.beam().power().i - 0.5 * (ts * ts + tp * tp) * scale).abs() < 2e-6);
            assert!((result.beam().power().q - 0.5 * (ts * ts - tp * tp) * scale).abs() < 2e-6);
        }
        assert!(
            positions[0] < positions[1] - 0.02,
            "blue light must refract more strongly"
        );
        // This receiver is much narrower than the unbent beam's displacement.
        assert!((0.35 - receiver_x).abs() > 0.4);
    });
}

#[test]
fn g1_native_sphere_tir_reaches_interior_receiver_with_phase_and_absorption() {
    with_cx(|cx, _| {
        // The first hit is (.8,0,.6): sin(i)=.8 exceeds 1/1.5.
        // Reflection is (-.96,0,.28); the detector at z=.8 is reached before
        // the next sphere boundary. There is no triangulated lens geometry.
        let after = 0.2 / 0.28;
        let receiver_x = 0.8 - 0.96 * after;
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, 0.8, receiver_x, 0.04)]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.8, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(),
            )
            .unwrap();
        let mut interfaces = 0;
        let result = scene
            .trace_with_samples(cx, input, 0, settings(), |_| {
                interfaces += 1;
                [0.2, 0.7, 0.999_999]
            })
            .unwrap();
        assert_eq!(interfaces, 1);
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(result.medium_index(), 0);
        assert_eq!(result.beam().medium(), BeamMedium::Homogeneous(lens));
        assert_eq!(result.surface_events(), 2);
        assert!((result.beam().ray().origin.x - receiver_x).abs() < 4e-6);
        assert!((result.beam().ray().dir.x + 0.96).abs() < 2e-6);
        assert!((result.beam().ray().dir.z - 0.28).abs() < 2e-6);
        let length = 0.6 + after;
        assert!((result.traveled_m() - length).abs() < 5e-6);
        assert!((result.optical_path_length_m() - 1.5 * length).abs() < 8e-6);
        assert!((result.optical_depth() - 0.2 * length).abs() < 2e-6);
        // TIR amplitudes have unit modulus; their relative phase is the
        // difference of two independently evaluated Fresnel phase angles.
        let gamma = (1.5_f64.powi(2) * 0.8_f64.powi(2) - 1.0).sqrt();
        let phase = 2.0 * ((1.5 * gamma).atan2(0.6) - gamma.atan2(1.5 * 0.6));
        let attenuation = (-0.2 * length).exp();
        assert!((result.beam().power().i - attenuation).abs() < 2e-6);
        assert!((result.beam().power().u - attenuation * phase.cos()).abs() < 3e-6);
        assert!((result.beam().power().v - attenuation * phase.sin()).abs() < 3e-6);
    });
}

#[test]
fn g1_affine_gradient_field_refracts_through_native_sphere_to_finite_receiver() {
    with_cx(|cx, _| {
        // Intersect the independent closed-form GRIN trajectory with x²+z²=1.
        // The scalar root is monotone on [0,1]; no production ODE or geometry
        // intersection is reused to compute this reference endpoint.
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        for _ in 0..60 {
            let z = 0.5 * (lo + hi);
            let x = 6.0 * ((z / 6.0).cosh() - 1.0);
            if x * x + z * z > 1.0 {
                hi = z;
            } else {
                lo = z;
            }
        }
        let z = 0.5 * (lo + hi);
        let x = 6.0 * ((z / 6.0).cosh() - 1.0);
        let inside_length = 6.0 * (z / 6.0).sinh();
        let index = 1.2 + 0.2 * x;
        let tangent = Vec3::new(0.2 * inside_length / index, 0.0, 1.2 / index);
        let normal = Vec3::new(x, 0.0, z);
        let cosine = normal.dot(tangent);
        let transmitted_cosine = (1.0 - index * index * (1.0 - cosine * cosine)).sqrt();
        // Tangential optical momentum n*t is conserved at the local surface.
        let outgoing = Vec3::new(
            index * (tangent.x - cosine * normal.x) + transmitted_cosine * normal.x,
            0.0,
            index * (tangent.z - cosine * normal.z) + transmitted_cosine * normal.z,
        );
        let after = (2.0 - z) / outgoing.z;
        let receiver_x = x + after * outgoing.x;
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, 2.0, receiver_x, 0.04)]).unwrap();
        let field = SlabField;
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([17; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(9), 100_000).unwrap();
        let input = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let result = scene
            .trace_with_samples(cx, input, 0, settings(), |_| [0.2, 0.7, 0.999_999])
            .unwrap();
        assert_eq!(result.termination(), PathTermination::Detected);
        assert_eq!(result.surface_events(), 2);
        assert_eq!(result.medium_index(), 1);
        assert_eq!(result.medium_snapshot_identity(), None);
        assert_eq!(result.last_contact().unwrap().object_id, 9);
        assert!((result.beam().ray().origin.x - receiver_x).abs() < 3e-6);
        assert!((result.beam().ray().dir.x - outgoing.x).abs() < 2e-6);
        assert!((result.traveled_m() - inside_length - after).abs() < 4e-6);
        let optical_length = 0.5 * (inside_length * index + 1.2 * z) + after;
        assert!((result.optical_path_length_m() - optical_length).abs() < 5e-6);
        assert!((result.optical_depth() - 0.3 * inside_length).abs() < 1e-6);
        assert!((result.beam().power().i - (-0.3 * inside_length).exp()).abs() < 2e-6);
    });
}

#[test]
fn g0_tangent_native_sphere_does_not_invent_a_refractive_crossing() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![sphere()]).unwrap();
        let lens = glass(1.5, 0.0);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![bindings(9)[0]], 100_000)
                .unwrap();
        let input = scene
            .launch_beam(
                Point3::new(1.0, 0.0, -2.0),
                frame(),
                550.0,
                1,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        assert!(matches!(
            scene.trace_with_samples(cx, input, 1, settings(), |_| {
                panic!("a grazing touch has no admitted entering interface")
            }),
            Err(CurvedSceneError::Geometry(ProximityError::AmbiguousContact))
        ));
        assert_eq!(input.ray().origin, Point3::new(1.0, 0.0, -2.0));
        assert_eq!(input.power(), Stokes::unpolarized(1.0));
    });
}
