use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, StreamKey};

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn rough() -> DielectricSurface {
    DielectricSurface::try_anisotropic_rough(0.14, 0.65, Vec3::new(1.0, 0.0, 0.0)).unwrap()
}

fn scene(placement: RigidTransform) -> Scene {
    Scene {
        primitives: vec![Primitive {
            shape: Shape::Instance(
                GeometryInstance::try_new(
                    83,
                    hash_domain("org.frankensim.test.anisotropic-glass", b"plane"),
                    SharedGeometry::mesh(TriMesh::new(
                        vec![
                            [-2.0, -2.0, 0.0],
                            [2.0, -2.0, 0.0],
                            [2.0, 2.0, 0.0],
                            [-2.0, 2.0, 0.0],
                        ],
                        vec![[0, 1, 2], [0, 2, 3]],
                    )),
                    placement,
                )
                .unwrap(),
            ),
            material: Material::Dielectric {
                glass: glass(),
                surface: rough(),
            },
            emission: None,
        }],
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.2,
        },
    }
}

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3.0e-10 * a.abs().max(b.abs()).max(1.0e-8),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn native_camera_and_bidirectional_paths_share_both_anisotropic_branches_and_reverse_pdfs() {
    let half = core::f64::consts::FRAC_1_SQRT_2;
    let placement = RigidTransform::try_new([0.0, 0.0, half, half], [0.0; 3]).unwrap();
    let scene = scene(placement);
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.2, 0.3, 1.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let intersection = intersect(&scene, cx, &ray, None).unwrap().unwrap();
        let frame = surface_frame(&intersection.hit, &ray).unwrap();
        let stack = MediumStack::new();
        let vertex =
            Vertex::surface(&scene, intersection, ray, frame, &stack, 1.0, 1.0, None).unwrap();
        let tangent = vertex.brushed_frame.unwrap().tangent;
        // A quaternion quarter-turn cancels unit-scale terms here; use an
        // absolute roundoff bound only for this theoretically zero component.
        assert!(
            tangent.x.abs() <= 8.0 * f64::EPSILON,
            "rotated tangent x={:.17e}",
            tangent.x
        );
        close(tangent.y, 1.0);
        let Material::Dielectric { glass, surface } = &scene.primitives[0].material else {
            unreachable!()
        };
        let bound = dielectric_surface_in_frame(*surface, vertex.brushed_frame.as_ref()).unwrap();
        let boundary = boundary_media(0, glass, true, &stack).unwrap();
        let wavelengths = [410.0, 500.0, 600.0, 730.0];
        let mut events = [false; 2];
        for branch in [0.0, 0.999] {
            let random = [0.25, 0.35, branch, 0.5];
            let camera = sample_dielectric_path(
                frame.oriented,
                ray.dir.scale(-1.0),
                bound,
                &boundary,
                &wavelengths,
                0,
                random[0],
                random[1],
                random[2],
            )
            .unwrap()
            .unwrap();
            let bdpt = sample_surface(&scene, &vertex, 410.0, TransportMode::Radiance, random)
                .unwrap()
                .unwrap();
            close(camera.direction.dot(bdpt.direction), 1.0);
            close(camera.weights[0], bdpt.weight);
            close(camera.pdf, bdpt.pdf_fwd_solid_angle);
            close(
                reverse_dielectric_pdf(
                    frame.oriented,
                    ray.dir.scale(-1.0),
                    camera.direction,
                    bound,
                    boundary,
                    410.0,
                )
                .unwrap(),
                bdpt.pdf_rev_solid_angle,
            );
            let importance =
                sample_surface(&scene, &vertex, 410.0, TransportMode::Importance, random)
                    .unwrap()
                    .unwrap();
            let transmission = camera.event == DielectricEvent::Transmission;
            events[usize::from(transmission)] = true;
            close(
                importance.weight,
                bdpt.weight * if transmission { 1.5_f64.powi(2) } else { 1.0 },
            );
            let value = surface_scattering(
                &scene,
                &vertex,
                ray.dir.scale(-1.0),
                camera.direction,
                410.0,
                TransportMode::Radiance,
            )
            .unwrap();
            close(
                value * frame.oriented.dot(camera.direction).abs() / camera.pdf,
                camera.weights[0],
            );
        }
        assert!(events.into_iter().all(|event| event));
    });
}

#[test]
fn dispersive_packet_lanes_and_material_identity_keep_the_authored_finish() {
    let glass = DielectricGlass::representative_borosilicate();
    let surface = rough();
    let material = Material::Dielectric { glass, surface };
    assert_ne!(
        material.content_identity(),
        Material::Dielectric {
            glass,
            surface: DielectricSurface::try_rough(0.14).unwrap(),
        }
        .content_identity()
    );
    assert_ne!(
        material.content_identity(),
        Material::Dielectric {
            glass,
            surface: surface.with_tangent_axis(Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        }
        .content_identity()
    );
    let stack = MediumStack::new();
    let boundary = boundary_media(0, &glass, true, &stack).unwrap();
    let n = Vec3::new(0.0, 0.0, 1.0);
    let wavelengths = [390.0, 470.0, 610.0, 770.0];
    let continuations = sample_dispersive_dielectric_lanes(
        n,
        n,
        surface,
        &boundary,
        &wavelengths,
        0.25,
        0.35,
        0.99,
    )
    .unwrap();
    for lane in 0..PACKET {
        let actual = continuations[lane].unwrap();
        let eta = glass.ior().eval(wavelengths[lane]).unwrap();
        let expected = surface
            .interface(wavelengths[lane], false)
            .unwrap()
            .sample_rough(n, n, 1.0, eta, 0.14, 0.25, 0.35, 0.99)
            .unwrap()
            .unwrap();
        close(actual.weight, expected.radiance_weight);
        close(actual.pdf, expected.pdf);
        close(actual.direction.dot(expected.direction), 1.0);
    }
    assert!(
        (continuations[0].unwrap().direction.x - continuations[3].unwrap().direction.x).abs()
            > 1.0e-6
    );
}

#[test]
fn polarized_queries_retain_the_same_anisotropic_scalar_proposal() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let frame = PolarizationFrame::try_new(n.scale(-1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let vacuum = DielectricGlass::new(
        CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    let beam = PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        frame,
        550.0,
        vacuum,
        Stokes::unpolarized(1.0),
    )
    .unwrap();
    let sample = beam
        .sample_interface(n, glass(), rough(), [0.2, 0.3, 0.4])
        .unwrap()
        .unwrap();
    let scalar = rough()
        .interface(550.0, false)
        .unwrap()
        .sample_rough(n, n, 1.0, 1.5, 0.14, 0.2, 0.3, 0.4)
        .unwrap()
        .unwrap();
    assert_eq!(sample.pdf, scalar.pdf);
    close(sample.beam.frame().direction().dot(scalar.direction), 1.0);
    let source =
        PolarizationFrame::try_new(scalar.direction.scale(-1.0), frame.horizontal()).unwrap();
    let connection = beam
        .evaluate_dielectric_connection(n, glass(), rough(), source, Stokes::unpolarized(1.0))
        .unwrap();
    let evaluation = rough()
        .interface(550.0, false)
        .unwrap()
        .evaluate_rough(n, n, scalar.direction, 1.0, 1.5, 0.14)
        .unwrap();
    close(connection.factor, evaluation.value);
    close(connection.pdf, scalar.pdf);
}
