//! Full reciprocal lens weights against physical scene-to-sensor Stokes transport.
#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, LensSample, LensSurface, PhysicalCamera,
    SequentialLens,
};
use fs_render::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use fs_render::polarization::Stokes;
use fs_render::polarization::dielectric::fresnel_dielectric_polarized;
use fs_render::polarization::dielectric::transport::PolarizationFrame;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 12,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}

fn camera(rotated: bool, sigma: f64, stopped: bool) -> PhysicalCamera {
    let glass = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(sigma).unwrap(),
        GlassProvenance::Custom,
    );
    let mut surfaces = vec![LensSurface::try_new(0.08, 20.0, 0.012, Some(glass)).unwrap()];
    if stopped {
        surfaces.push(LensSurface::try_new(0.085, 0.0, 0.0001, Some(glass)).unwrap());
    }
    surfaces.push(LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap());
    let (eye, target, up) = if rotated {
        (
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(2.0, 2.0, 3.0),
            Vec3::new(0.0, 0.0, 1.0),
        )
    } else {
        (
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
        )
    };
    PhysicalCamera::try_look_at(
        eye,
        target,
        up,
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(SequentialLens::try_new(550.0, &surfaces).unwrap())
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-12, "{actual} != {expected}");
}

fn components(s: Stokes) -> [f64; 4] {
    [s.i, s.q, s.u, s.v]
}
fn scale(s: Stokes, a: f64) -> Stokes {
    Stokes {
        i: a * s.i,
        q: a * s.q,
        u: a * s.u,
        v: a * s.v,
    }
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
fn world(camera: &PhysicalCamera, v: Vec3) -> Vec3 {
    let (r, u, f) = (camera.right(), camera.up(), camera.forward());
    Vec3::new(
        r.x * v.x + u.x * v.y + f.x * v.z,
        r.y * v.x + u.y * v.y + f.y * v.z,
        r.z * v.x + u.z * v.y + f.z * v.z,
    )
}
fn express(s: Stokes, from: PolarizationFrame, to: PolarizationFrame) -> Stokes {
    s.rotate_axes(
        from.vertical()
            .dot(to.horizontal())
            .atan2(from.horizontal().dot(to.horizontal())),
    )
}

// Forward physical transport in world coordinates, using independently rebuilt
// incidence frames. It neither reverses detector rows nor enters lens-local
// Stokes coordinates, so it detects both reciprocity and handedness mistakes.
fn forward_interface(
    state: Stokes,
    frame: PolarizationFrame,
    next_direction: Vec3,
    normal: Vec3,
    ni: f64,
    nt: f64,
) -> (Stokes, PolarizationFrame) {
    let s = unit(cross(frame.direction(), normal));
    let incident = PolarizationFrame::try_new(frame.direction(), s).unwrap();
    let response = fresnel_dielectric_polarized(-frame.direction().dot(normal), ni, nt).unwrap();
    assert!(!response.scalar.total_internal_reflection);
    (
        response.transmission.apply(express(state, frame, incident)),
        PolarizationFrame::try_new(next_direction, s).unwrap(),
    )
}

#[test]
fn skew_lens_matches_forward_stokes_with_sensor_axis_and_world_handedness() {
    with_cx(|_, cx| {
        for rotated in [false, true] {
            let camera = camera(rotated, 7.0, false);
            let axis = world(&camera, Vec3::new(0.3, 0.8, 0.2));
            let sample = LensSample::try_new(0.75, 0.5).unwrap();
            let row = Stokes {
                i: 0.5,
                q: 0.2,
                u: -0.1,
                v: 0.3,
            };
            let beam = camera
                .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                    cx,
                    0.08,
                    0.045,
                    sample,
                    &[550.0],
                    axis,
                    row,
                )
                .unwrap()[0]
                .unwrap();
            let ray = beam.ray();
            assert_eq!(ray.dir, beam.frame().direction());
            assert_eq!(beam.wavelength_nm(), 550.0);

            // This disk sample is exactly half the rear clear radius on +x.
            let x = 0.006_f64;
            let rear_z = 0.08 + 20.0 * x * x / (1.0 + (1.0 - (20.0 * x).powi(2)).sqrt());
            let first = camera
                .eye()
                .offset(world(&camera, Vec3::new(x, 0.0, rear_z)));
            let sensor = camera
                .eye()
                .offset(world(&camera, Vec3::new(-0.08 * 0.08, -0.045 * 0.08, 0.0)));
            let sensor_direction = unit(first.delta_from(sensor));
            let inside = unit(ray.origin.delta_from(first));
            let exit = ray.origin.delta_from(camera.eye());
            let exit = Vec3::new(
                exit.dot(camera.right()),
                exit.dot(camera.up()),
                exit.dot(camera.forward()),
            );
            let front_normal = world(
                &camera,
                unit(Vec3::new(
                    -20.0 * exit.x,
                    -20.0 * exit.y,
                    -20.0 * (exit.z - 0.09) - 1.0,
                )),
            )
            .scale(-1.0);
            let rear_normal = world(
                &camera,
                unit(Vec3::new(20.0 * x, 0.0, 20.0 * (rear_z - 0.08) - 1.0)),
            )
            .scale(-1.0);

            let source = Stokes::try_new(2.0, 0.7, -0.8, 0.9).unwrap();
            let source_frame =
                PolarizationFrame::try_new(ray.dir.scale(-1.0), beam.frame().horizontal()).unwrap();
            let (state, frame) = forward_interface(
                source,
                source_frame,
                inside.scale(-1.0),
                front_normal,
                1.0,
                1.5,
            );
            let state = scale(state, (-7.0 * ray.origin.delta_from(first).norm()).exp());
            let (state, frame) = forward_interface(
                state,
                frame,
                sensor_direction.scale(-1.0),
                rear_normal,
                1.5,
                1.0,
            );
            let detector_frame =
                PolarizationFrame::try_new(sensor_direction.scale(-1.0), axis).unwrap();
            let detected = express(state, frame, detector_frame);
            let forward = components(row)
                .into_iter()
                .zip(components(detected))
                .map(|(a, b)| a * b)
                .sum::<f64>();
            let b = beam.power();
            let reciprocal = source.i * b.i + source.q * b.q - source.u * b.u + source.v * b.v;
            close(reciprocal, forward);
            assert!(
                (reciprocal - source.i * b.i).abs() > 0.05,
                "scalar transmission loses the polarization measurement"
            );

            let reversed_axis = camera
                .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                    cx,
                    0.08,
                    0.045,
                    sample,
                    &[550.0],
                    axis.scale(-3.0),
                    row,
                )
                .unwrap()[0]
                .unwrap();
            for (a, b) in components(beam.power())
                .into_iter()
                .zip(components(reversed_axis.power()))
            {
                close(a, b);
            }
        }
    });
}

#[test]
fn normal_incidence_retains_circular_sign_and_complementary_detector_throughput() {
    with_cx(|_, cx| {
        for rotated in [false, true] {
            let camera = camera(rotated, 7.0, false);
            let row = Stokes {
                i: 0.5,
                q: 0.1,
                u: 0.2,
                v: 0.3,
            };
            let opposite = Stokes {
                i: 0.5,
                q: -row.q,
                u: -row.u,
                v: -row.v,
            };
            let generate = |row| {
                camera
                    .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                        cx,
                        0.0,
                        0.0,
                        LensSample::CENTER,
                        &[450.0, 650.0],
                        camera.right(),
                        row,
                    )
                    .unwrap()
            };
            let selected = generate(row);
            let rejected = generate(opposite);
            let identity = camera
                .generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[450.0, 650.0],
                )
                .unwrap();
            assert_eq!(generate(Stokes::unpolarized(1.0)), identity);
            let scalar = camera
                .generate_sequential_spectral_rays_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[450.0, 650.0],
                )
                .unwrap();
            let t = 0.96_f64.powi(2) * (-0.07_f64).exp();
            for i in 0..2 {
                let beam = selected[i].unwrap();
                let expected = scale(Stokes { u: -row.u, ..row }, t);
                for (a, b) in components(beam.power())
                    .into_iter()
                    .zip(components(expected))
                {
                    close(a, b);
                }
                close(beam.frame().horizontal().dot(camera.right()), 1.0);
                close(beam.geometry_weight(), 1.0);
                let sum = components(beam.power())
                    .into_iter()
                    .zip(components(rejected[i].unwrap().power()))
                    .map(|(a, b)| a + b);
                for (a, b) in sum.zip(components(identity[i].unwrap().power())) {
                    close(a, b);
                }
                close(
                    scalar[i].unwrap().weight,
                    identity[i].unwrap().power().i * identity[i].unwrap().geometry_weight(),
                );
                assert_eq!(beam.wavelength_nm(), [450.0, 650.0][i]);
            }
        }
    });
}

#[test]
fn clipped_paths_invalid_analyzers_and_attenuation_underflow_are_distinct() {
    with_cx(|gate, cx| {
        let stopped = camera(false, 0.0, true);
        let sample = LensSample::try_new(0.75, 0.5).unwrap();
        assert_eq!(
            stopped
                .generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    sample,
                    &[450.0, 650.0],
                )
                .unwrap(),
            [None, None]
        );
        assert!(
            stopped
                .generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    sample,
                    &[f64::NAN],
                )
                .is_err()
        );
        let dark = camera(false, 1e6, false);
        assert_eq!(
            dark.generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            ),
            Err(CameraError::SequentialLensOutsideModel)
        );
        assert_eq!(
            dark.generate_sequential_spectral_rays_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
                .unwrap()
                .weight,
            0.0
        );

        let camera = camera(false, 0.0, false);
        let query = |axis, row| {
            camera.generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
                axis,
                row,
            )
        };
        for row in [
            Stokes::unpolarized(0.0),
            Stokes::unpolarized(1.01),
            Stokes {
                i: 0.8,
                q: 0.3,
                u: 0.0,
                v: 0.0,
            },
            Stokes {
                i: 0.5,
                q: 0.0,
                u: f64::NAN,
                v: 0.0,
            },
        ] {
            assert!(query(camera.right(), row).is_err());
        }
        assert!(query(Vec3::new(0.0, 0.0, 0.0), Stokes::unpolarized(1.0)).is_err());
        let circular = Stokes {
            i: 0.5,
            q: 0.0,
            u: 0.0,
            v: 0.5,
        };
        assert!(query(camera.forward(), circular).is_err());
        assert!(query(camera.forward(), Stokes::unpolarized(1.0)).is_ok());
        gate.request();
        assert_eq!(query(camera.right(), circular), Err(CameraError::Cancelled));
    });
}
