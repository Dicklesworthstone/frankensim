use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::photometry::source::{
    AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, SourceSpectrum, photopic_weight_nm,
};
use crate::spectral::lift_rgb;
use crate::tracer::{Camera, Material, Primitive, Shape, TracerError};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

mod history;
mod reflection;

fn scope(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 71,
                    kernel_id: 52,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn source(position: [f64; 3], distribution: AngularDistribution) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        4.0 * core::f64::consts::PI,
        1.0,
        distribution,
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn scene(ior: CauchyIor, sigma: f64) -> Scene {
    let mesh = TriMesh::new(
        vec![
            [-6.0, -6.0, 0.0],
            [6.0, -6.0, 0.0],
            [6.0, 6.0, 0.0],
            [-6.0, 6.0, 0.0],
            [-6.0, -6.0, 0.2],
            [6.0, -6.0, 0.2],
            [6.0, 6.0, 0.2],
            [-6.0, 6.0, 0.2],
        ],
        vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    );
    Scene {
        primitives: vec![Primitive {
            shape: Shape::Instance(
                GeometryInstance::try_new(
                    1,
                    fs_blake3::hash_domain("fs-render:slab-photometry-test", b"slab"),
                    SharedGeometry::mesh(mesh),
                    RigidTransform::identity(),
                )
                .unwrap(),
            ),
            material: Material::Dielectric {
                glass: DielectricGlass::new(
                    ior,
                    BeerLambertAbsorption::try_constant(sigma).unwrap(),
                    GlassProvenance::Custom,
                ),
                surface: DielectricSurface::SMOOTH,
            },
            emission: None,
        }],
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, -2.0),
            forward: Vec3::new(0.0, 0.0, 1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.5,
        },
    }
}

fn probe(x: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x, 0.0, 1.2), Vec3::new(0.0, 0.0, -1.0)).unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-14),
        "{a:.16e} != {b:.16e}"
    );
}

fn panel(x: f64, z: f64, half: f64) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [x - half, -half, z],
                [x + half, -half, z],
                [x + half, half, z],
                [x - half, half, z],
            ],
            vec![[0, 2, 1], [0, 3, 2]],
        )),
        material: Material::Lambertian {
            reflectance: lift_rgb([0.0; 3]),
        },
        emission: None,
    }
}

#[test]
fn g1_normal_slab_image_has_exact_fresnel_beer_and_paraxial_spreading() {
    scope(|cx, _| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        for index in [1.0, 1.5, 2.0] {
            let scene = scene(CauchyIor::try_constant(index).unwrap(), 0.7);
            let image = source
                .transmitted_slab_irradiance(cx, &scene, 0, probe(0.0), TransportOptions::default())
                .unwrap();
            let reflection = ((index - 1.0) / (index + 1.0)).powi(2);
            let expected =
                (1.0 - reflection).powi(2) * (-0.7_f64 * 0.2).exp() / (2.0 + 0.2 / index).powi(2);
            close(image.reading().illuminance_lux, expected, 3e-13);
            assert!(image.unoccluded().iter().all(|v| *v));
            assert_eq!(image.launch_directions()[0], Vec3::new(0.0, 0.0, 1.0));
        }
    });
}

// Independent forward Snell ray and receiving-plane intersection. This does
// not use the native connector, its inverse solve or its analytic Jacobian.
fn endpoint(direction: Vec3, index: f64, receiver: IlluminanceProbe) -> Point3 {
    let inside = Vec3::new(
        direction.x / index,
        direction.y / index,
        (1.0 - (direction.x * direction.x + direction.y * direction.y) / (index * index)).sqrt(),
    );
    let exit = Point3::new(0.0, 0.0, -1.0)
        .offset(direction.scale(1.0 / direction.z))
        .offset(inside.scale(0.2 / inside.z));
    let parameter = receiver
        .position_m()
        .delta_from(exit)
        .dot(receiver.normal())
        / direction.dot(receiver.normal());
    exit.offset(direction.scale(parameter))
}

fn reference_area(direction: Vec3, index: f64, receiver: IlluminanceProbe) -> f64 {
    let epsilon = 2e-5;
    let u = Vec3::new(direction.z, 0.0, -direction.x);
    let v = Vec3::new(0.0, 1.0, 0.0);
    let derivative = |axis: Vec3| {
        let plus = Vec3::new(
            direction.x + epsilon * axis.x,
            direction.y + epsilon * axis.y,
            direction.z + epsilon * axis.z,
        );
        let minus = Vec3::new(
            direction.x - epsilon * axis.x,
            direction.y - epsilon * axis.y,
            direction.z - epsilon * axis.z,
        );
        endpoint(plus.scale(1.0 / plus.norm()), index, receiver)
            .delta_from(endpoint(minus.scale(1.0 / minus.norm()), index, receiver))
            .scale(0.5 / epsilon)
    };
    let a = derivative(u);
    let b = derivative(v);
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
    .norm()
}

#[test]
fn g1_dispersive_oblique_image_preserves_sp_history_focusing_and_actual_source_angle() {
    scope(|cx, _| {
        let source = source(
            [0.0, 0.0, -1.0],
            AngularDistribution::cosine_power(4.0, core::f64::consts::FRAC_PI_2).unwrap(),
        );
        let law = CauchyIor::try_new(1.42, 0.018, 0.0).unwrap();
        let scene = scene(law, 0.3);
        let receiver =
            IlluminanceProbe::try_new(Point3::new(3.6, 0.0, 1.2), Vec3::new(0.1, 0.2, -1.0))
                .unwrap();
        let image = source
            .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
            .unwrap();
        let mut expected_lux = 0.0;
        let mut expected_w = 0.0;
        let mut depolarized_lux = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let lambda = wavelength_nm(bin);
            let n = law.eval(lambda).unwrap();
            let d = image.launch_directions()[bin];
            assert!(
                endpoint(d, n, receiver)
                    .delta_from(receiver.position_m())
                    .norm()
                    < 1e-11
            );
            let cg = (1.0 - (d.x * d.x + d.y * d.y) / (n * n)).sqrt();
            let ts = 4.0 * n * d.z * cg / (d.z + n * cg).powi(2);
            let tp = 4.0 * n * d.z * cg / (n * d.z + cg).powi(2);
            let scale = source.spectral_intensity_w_sr_nm(d, lambda).unwrap()
                * (-0.3 * 0.2 / cg).exp()
                / reference_area(d, n, receiver);
            let expected = scale * 0.5 * (ts * ts + tp * tp);
            close(image.spectral_irradiance_w_m2_nm()[bin], expected, 5e-8);
            expected_lux += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
            expected_w += expected * BIN_WIDTH_NM;
            depolarized_lux += scale
                * (0.5 * (ts + tp)).powi(2)
                * photopic_weight_nm(bin)
                * LUMINOUS_EFFICACY_LM_W;
        }
        close(image.reading().illuminance_lux, expected_lux, 5e-8);
        close(image.reading().radiant_irradiance_w_m2, expected_w, 5e-8);
        assert!((expected_lux - depolarized_lux) / expected_lux > 0.005);
        assert!((image.launch_directions()[0].x - image.launch_directions()[79].x).abs() > 1e-4);
    });
}

#[test]
fn g3_visibility_uses_refracted_segments_and_preserves_receiver_endpoint() {
    scope(|cx, _| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let mut scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.0);
        let receiver = probe(3.6);
        let expected = source
            .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
            .unwrap();
        let direction = expected.launch_directions()[0];
        let chord_x = 3.6 * 0.5 / 2.2;
        assert!((chord_x - 0.5 * direction.x / direction.z).abs() > 0.01);
        // This occludes the straight chord but not the true source segment.
        scene.primitives.push(panel(chord_x, -0.5, 0.005));
        // The normal discovery probe is a geometric query, not an optical
        // segment. This embedded opaque panel must not veto the true image.
        scene.primitives.push(panel(3.6 / 2.2, 0.1, 0.005));
        scene.primitives.push(panel(3.6, 1.2, 0.2));
        let clear = source
            .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
            .unwrap();
        assert_eq!(clear, expected);
        let true_x = 3.6 - 0.5 * direction.x / direction.z;
        scene.primitives.push(panel(true_x, 0.7, 0.005));
        let blocked = source
            .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
            .unwrap();
        assert_eq!(blocked.reading().illuminance_lux, 0.0);
        assert!(blocked.unoccluded().iter().all(|v| !v));
    });
}

#[test]
fn g0_selected_image_refuses_unsupported_media_roughness_range_and_cancellation() {
    scope(|cx, gate| {
        let source = source([0.0, 0.0, -1.0], AngularDistribution::isotropic());
        let mut scene = scene(CauchyIor::try_constant(1.5).unwrap(), 0.0);
        let receiver = probe(3.6);
        let chord = receiver.position_m().delta_from(source.position_m()).norm();
        assert!(matches!(
            source.transmitted_slab_irradiance(
                cx,
                &scene,
                0,
                receiver,
                TransportOptions {
                    max_distance_m: chord,
                    ..TransportOptions::default()
                }
            ),
            Err(PhotometryError::DistanceLimitExceeded)
        ));
        let inside =
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.1), Vec3::new(0.0, 0.0, -1.0))
                .unwrap();
        assert!(
            source
                .transmitted_slab_irradiance(cx, &scene, 0, inside, TransportOptions::default())
                .is_err()
        );
        let back =
            IlluminanceProbe::try_new(receiver.position_m(), Vec3::new(0.0, 0.0, 1.0)).unwrap();
        assert!(
            source
                .transmitted_slab_irradiance(cx, &scene, 0, back, TransportOptions::default())
                .is_err()
        );
        let Shape::Instance(instance) = &scene.primitives[0].shape else {
            unreachable!()
        };
        let SharedGeometry::Mesh(mesh) = instance.geometry() else {
            unreachable!()
        };
        let open_shape = Shape::Instance(
            GeometryInstance::try_new(
                1,
                fs_blake3::hash_domain("fs-render:slab-photometry-test", b"open-sheets"),
                SharedGeometry::mesh(TriMesh::new(
                    mesh.vertices.clone(),
                    mesh.triangles[..4].to_vec(),
                )),
                instance.transform(),
            )
            .unwrap(),
        );
        let closed_shape = std::mem::replace(&mut scene.primitives[0].shape, open_shape);
        assert!(matches!(
            source.transmitted_slab_irradiance(
                cx,
                &scene,
                0,
                receiver,
                TransportOptions::default()
            ),
            Err(PhotometryError::Geometry(TracerError::InitialMedium {
                boundary_primitive: 0,
                ..
            }))
        ));
        scene.primitives[0].shape = closed_shape;
        if let Material::Dielectric { surface, .. } = &mut scene.primitives[0].material {
            *surface = DielectricSurface::try_rough(0.2).unwrap();
        }
        assert!(
            source
                .transmitted_slab_irradiance(cx, &scene, 0, receiver, TransportOptions::default())
                .is_err()
        );
        gate.request();
        assert!(matches!(
            source.transmitted_slab_irradiance(
                cx,
                &scene,
                0,
                receiver,
                TransportOptions::default()
            ),
            Err(PhotometryError::Cancelled)
        ));
    });
}
