//! G0/G3 camera tests: physical admission, ideal thin-lens geometry,
//! deterministic aperture sampling, animation, hard cuts, and cancellation.
#![cfg(feature = "chart-backends")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraError, CameraKeyframe, CameraProjection, CameraShot, CutSide,
    DefocusedDiffraction, FocusedDiffraction, LensSample, PhysicalCamera,
};
use fs_render::instances::RigidTransform;
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};

const TOLERANCE: f64 = 2.0e-12;

fn with_cx<R>(operation: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4341_4d45_5241,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&gate, &cx)
    })
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{context}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn assert_point_close(actual: Point3, expected: Point3, context: &str) {
    assert_close(actual.x, expected.x, TOLERANCE, context);
    assert_close(actual.y, expected.y, TOLERANCE, context);
    assert_close(actual.z, expected.z, TOLERANCE, context);
}

fn assert_vec_close(actual: Vec3, expected: Vec3, context: &str) {
    assert_close(actual.x, expected.x, TOLERANCE, context);
    assert_close(actual.y, expected.y, TOLERANCE, context);
    assert_close(actual.z, expected.z, TOLERANCE, context);
}

fn canonical_projection() -> CameraProjection {
    CameraProjection::try_half_tangent(0.5).unwrap()
}

fn pinhole() -> Aperture {
    Aperture::try_circular(0.0).unwrap()
}

fn look_at_camera(
    eye: Point3,
    target: Point3,
    focus_distance_m: f64,
    aperture: Aperture,
) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        eye,
        target,
        Vec3::new(0.0, 1.0, 0.0),
        canonical_projection(),
        focus_distance_m,
        aperture,
    )
    .unwrap()
}

fn axial_camera(aperture: Aperture, focus_distance_m: f64) -> PhysicalCamera {
    PhysicalCamera::try_legacy_compatible(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        0.5,
        focus_distance_m,
        aperture,
    )
    .unwrap()
}

fn lens_coordinates(camera: &PhysicalCamera, origin: Point3) -> [f64; 2] {
    let displacement = origin.delta_from(camera.eye());
    [
        displacement.dot(camera.right()),
        displacement.dot(camera.up()),
    ]
}

fn shutter(open_s: f64, close_s: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        open_s,
        close_s - open_s,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(0.0, 2.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn physical_projection_and_f_number_admission_are_fail_closed() {
    let physical = CameraProjection::try_focal_sensor(0.050, 0.024).unwrap();
    assert_eq!(physical.focal_length_m(), Some(0.050));
    assert_close(
        physical.vertical_half_tan(),
        0.24,
        f64::EPSILON,
        "50 mm lens on 24 mm sensor",
    );

    let fov = CameraProjection::try_vertical_fov(2.0 * 0.24_f64.atan()).unwrap();
    assert_close(
        fov.vertical_half_tan(),
        physical.vertical_half_tan(),
        4.0 * f64::EPSILON,
        "physical projection and equivalent FOV",
    );
    assert_eq!(fov.focal_length_m(), None);

    for invalid in [
        CameraProjection::try_focal_sensor(0.0, 0.024),
        CameraProjection::try_focal_sensor(0.050, -0.024),
        CameraProjection::try_focal_sensor(f64::NAN, 0.024),
        CameraProjection::try_vertical_fov(0.0),
        CameraProjection::try_vertical_fov(core::f64::consts::PI),
        CameraProjection::try_vertical_fov(f64::INFINITY),
        CameraProjection::try_half_tangent(-f64::EPSILON),
    ] {
        assert_eq!(invalid, Err(CameraError::InvalidProjection));
    }

    let aperture = Aperture::try_from_f_number(0.050, 2.0).unwrap();
    assert_close(aperture.radius_m(), 0.0125, f64::EPSILON, "f/2 radius");
    let physical_camera = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        physical,
        1.0,
        aperture.clone(),
    )
    .unwrap();
    assert_close(
        physical_camera.f_number().unwrap(),
        2.0,
        f64::EPSILON,
        "derived physical f-number",
    );
    let extreme_projection =
        CameraProjection::try_focal_sensor(f64::from_bits(1), f64::from_bits(1)).unwrap();
    let extreme_camera = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        extreme_projection,
        1.0,
        Aperture::try_circular(f64::MAX).unwrap(),
    )
    .unwrap();
    assert_eq!(
        extreme_camera.f_number(),
        None,
        "an underflowed f-number must not be reported as physical zero"
    );
    assert!(Aperture::try_circular(0.0).unwrap().is_pinhole());
    assert!(
        Aperture::try_regular_polygon(0.0, 7, 0.3)
            .unwrap()
            .is_pinhole()
    );
    for invalid in [
        Aperture::try_from_f_number(0.0, 2.0),
        Aperture::try_from_f_number(0.050, 0.0),
        Aperture::try_from_f_number(0.050, f64::NAN),
        Aperture::try_from_f_number(f64::from_bits(1), 2.0),
        Aperture::try_circular(-0.001),
        Aperture::try_regular_polygon(0.01, 2, 0.0),
        Aperture::try_regular_polygon(0.01, 65, 0.0),
        Aperture::try_regular_polygon(0.01, 6, f64::NAN),
    ] {
        assert_eq!(invalid, Err(CameraError::InvalidAperture));
    }
    assert_eq!((LensSample::CENTER.u(), LensSample::CENTER.v()), (0.5, 0.5));
    for invalid in [
        LensSample::try_new(-f64::EPSILON, 0.5),
        LensSample::try_new(1.0, 0.5),
        LensSample::try_new(0.5, f64::NAN),
    ] {
        assert_eq!(invalid, Err(CameraError::InvalidLensSample));
    }
}

#[test]
fn pinhole_raster_importance_round_trips_camera_ray_and_solid_angle_density() {
    let camera = axial_camera(pinhole(), 3.0);
    let width = 400_u32;
    let height = 200_u32;
    let aspect = f64::from(width) / f64::from(height);
    let half_tan = camera.projection().vertical_half_tan();
    let tangent_pixel_area =
        4.0 * aspect * half_tan * half_tan / (f64::from(width) * f64::from(height));

    with_cx(|_, cx| {
        for (pixel_x, pixel_y) in [(0_u32, 0_u32), (123, 77), (399, 199)] {
            let x_tan =
                (2.0 * (f64::from(pixel_x) + 0.5) / f64::from(width) - 1.0) * aspect * half_tan;
            let y_tan = (1.0 - 2.0 * (f64::from(pixel_y) + 0.5) / f64::from(height)) * half_tan;
            let ray = camera
                .generate_ray_from_tangent_offsets(cx, x_tan, y_tan, LensSample::CENTER)
                .unwrap();
            let point = ray.origin.offset(ray.dir.scale(7.0));
            let response = camera
                .pinhole_raster_sample(point, width, height)
                .unwrap()
                .expect("every pixel-centre ray projects inside its source pixel");
            assert_eq!(response.pixel, pixel_y * width + pixel_x);
            assert_vec_close(
                response.direction_from_camera,
                ray.dir,
                "light connection reverses the pixel-centre ray",
            );
            let optical_cosine = ray.dir.dot(camera.forward());
            let optical_cosine_cubed = optical_cosine * optical_cosine * optical_cosine;
            let expected_pdf = 1.0 / (tangent_pixel_area * optical_cosine_cubed);
            assert_close(
                response.pdf_solid_angle,
                expected_pdf,
                16.0 * f64::EPSILON * expected_pdf,
                "uniform tangent-pixel density converted to solid angle",
            );
            assert_close(
                response.depth_m,
                7.0 * optical_cosine,
                16.0 * f64::EPSILON * response.depth_m,
                "axial projection depth",
            );
        }
    });
}

#[test]
fn pinhole_raster_importance_refuses_wrong_camera_measure_and_excludes_boundaries() {
    let camera = axial_camera(pinhole(), 3.0);
    assert_eq!(
        camera.pinhole_raster_sample(Point3::new(0.0, 0.0, 1.0), 16, 9),
        Ok(None),
        "a point behind the pinhole has no sensor response",
    );
    assert_eq!(
        camera.pinhole_raster_sample(Point3::new(0.5, 0.0, -1.0), 16, 16),
        Ok(None),
        "the positive NDC boundary belongs to no half-open raster pixel",
    );
    assert_eq!(
        camera
            .pinhole_raster_sample(Point3::new(-0.5, 0.5, -1.0), 16, 16)
            .unwrap()
            .map(|sample| sample.pixel),
        Some(0),
        "the left and top boundaries belong to the first raster pixel",
    );
    assert_eq!(
        camera.pinhole_raster_sample(Point3::new(0.0, -0.5, -1.0), 16, 16),
        Ok(None),
        "the bottom NDC boundary belongs to no half-open raster pixel",
    );
    assert_eq!(
        camera.pinhole_raster_sample(Point3::new(0.0, 0.0, -1.0), 0, 16),
        Err(CameraError::InvalidProjection),
    );

    let finite_aperture = axial_camera(Aperture::try_circular(0.01).unwrap(), 3.0);
    assert_eq!(
        finite_aperture.pinhole_raster_sample(Point3::new(0.0, 0.0, -1.0), 16, 16),
        Err(CameraError::InvalidAperture),
        "an optical-centre splat cannot substitute for lens-area integration",
    );
}

#[test]
fn exact_pinhole_ray_matches_the_legacy_operation_order_and_known_bits() {
    let camera = axial_camera(pinhole(), 5.0);
    let x_tan = 0.25;
    let y_tan = -0.5;
    let legacy_raster = Vec3::new(x_tan, y_tan, -1.0);
    let legacy_dir = legacy_raster.scale(1.0 / legacy_raster.norm());

    with_cx(|_, cx| {
        let ray = camera
            .generate_ray_from_tangent_offsets(cx, x_tan, y_tan, LensSample::CENTER)
            .unwrap();
        let different_lens_draw = camera
            .generate_ray_from_tangent_offsets(
                cx,
                x_tan,
                y_tan,
                LensSample::try_new(0.03125, 0.96875).unwrap(),
            )
            .unwrap();

        assert_eq!(ray.origin, Point3::new(0.0, 0.0, 0.0));
        assert_eq!(ray, different_lens_draw, "pinhole must ignore lens draws");
        assert_eq!(ray.dir.x.to_bits(), legacy_dir.x.to_bits());
        assert_eq!(ray.dir.y.to_bits(), legacy_dir.y.to_bits());
        assert_eq!(ray.dir.z.to_bits(), legacy_dir.z.to_bits());
        assert_eq!(
            [
                ray.dir.x.to_bits(),
                ray.dir.y.to_bits(),
                ray.dir.z.to_bits()
            ],
            [
                0x3fcb_ee90_56fb_9c39,
                0xbfdb_ee90_56fb_9c39,
                0xbfeb_ee90_56fb_9c39,
            ]
        );
    });
}

#[test]
fn every_sampled_nonzero_lens_ray_converges_on_the_axial_focus_plane() {
    let focus_distance_m = 4.5;
    let apertures = [
        Aperture::try_circular(0.18).unwrap(),
        Aperture::try_regular_polygon(0.18, 7, 0.37).unwrap(),
    ];
    let lens_samples = [
        LensSample::try_new(0.07, 0.13).unwrap(),
        LensSample::try_new(0.23, 0.81).unwrap(),
        LensSample::try_new(0.61, 0.39).unwrap(),
        LensSample::try_new(0.93, 0.71).unwrap(),
    ];
    let raster_offsets = [(-0.3, -0.2), (0.0, 0.0), (0.4, -0.15)];

    with_cx(|_, cx| {
        for aperture in apertures {
            let camera = axial_camera(aperture, focus_distance_m);
            for (x_tan, y_tan) in raster_offsets {
                let raster = Vec3::new(
                    camera.forward().x + x_tan * camera.right().x + y_tan * camera.up().x,
                    camera.forward().y + x_tan * camera.right().y + y_tan * camera.up().y,
                    camera.forward().z + x_tan * camera.right().z + y_tan * camera.up().z,
                );
                let expected = camera.eye().offset(raster.scale(focus_distance_m));
                for lens_sample in lens_samples {
                    let ray = camera
                        .generate_ray_from_tangent_offsets(cx, x_tan, y_tan, lens_sample)
                        .unwrap();
                    assert!(ray.origin.delta_from(camera.eye()).norm() > 0.0);
                    let denominator = ray.dir.dot(camera.forward());
                    assert!(denominator > 0.0);
                    let at_focus = ray.at(focus_distance_m / denominator);
                    assert_point_close(at_focus, expected, "thin-lens focus convergence");
                }
            }
        }
    });
}

#[test]
fn camera_rays_are_equivariant_under_a_general_proper_rigid_transform() {
    let eye = Point3::new(0.7, -1.2, 2.4);
    let target = Point3::new(-0.3, 0.8, -1.1);
    let up_reference = Vec3::new(0.4, 1.0, 0.2);
    let aperture = Aperture::try_regular_polygon(0.12, 7, 0.31).unwrap();
    let camera = PhysicalCamera::try_look_at(
        eye,
        target,
        up_reference,
        canonical_projection(),
        3.7,
        aperture.clone(),
    )
    .unwrap();
    let half_angle = 0.37_f64;
    let axis_inverse_norm = 1.0 / 14.0_f64.sqrt();
    let sin_half = half_angle.sin();
    let transform = RigidTransform::try_new(
        [
            sin_half * axis_inverse_norm,
            2.0 * sin_half * axis_inverse_norm,
            3.0 * sin_half * axis_inverse_norm,
            half_angle.cos(),
        ],
        [1.1, -0.6, 0.9],
    )
    .unwrap();
    let transformed = PhysicalCamera::try_look_at(
        transform.transform_point(eye),
        transform.transform_point(target),
        transform.transform_vector(up_reference),
        canonical_projection(),
        3.7,
        aperture,
    )
    .unwrap();
    let lens_sample = LensSample::try_new(0.173, 0.827).unwrap();

    with_cx(|_, cx| {
        let ray = camera
            .generate_ray_from_tangent_offsets(cx, 0.23, -0.41, lens_sample)
            .unwrap();
        let transformed_ray = transformed
            .generate_ray_from_tangent_offsets(cx, 0.23, -0.41, lens_sample)
            .unwrap();
        assert_point_close(
            transformed_ray.origin,
            transform.transform_point(ray.origin),
            "rigidly transformed ray origin",
        );
        assert_vec_close(
            transformed_ray.dir,
            transform.transform_vector(ray.dir),
            "rigidly transformed ray direction",
        );
    });
}

#[test]
fn circular_aperture_is_contained_centered_and_has_uniform_area_moment() {
    let radius_m = 0.7;
    let camera = axial_camera(Aperture::try_circular(radius_m).unwrap(), 3.0);
    let grid = 128_u32;
    let mut sum = [0.0; 2];
    let mut sum_radius_squared = 0.0;

    with_cx(|_, cx| {
        for row in 0..grid {
            for column in 0..grid {
                let sample = LensSample::try_new(
                    (f64::from(column) + 0.5) / f64::from(grid),
                    (f64::from(row) + 0.5) / f64::from(grid),
                )
                .unwrap();
                let ray = camera
                    .generate_ray_from_tangent_offsets(cx, 0.0, 0.0, sample)
                    .unwrap();
                let displacement = ray.origin.delta_from(camera.eye());
                assert_close(
                    displacement.dot(camera.forward()),
                    0.0,
                    2.0e-15,
                    "circular aperture remains in lens plane",
                );
                let [x, y] = lens_coordinates(&camera, ray.origin);
                let radius_squared = x * x + y * y;
                assert!(
                    radius_squared <= radius_m * radius_m * (1.0 + 8.0 * f64::EPSILON),
                    "circular sample escaped: x={x}, y={y}"
                );
                sum[0] += x;
                sum[1] += y;
                sum_radius_squared += radius_squared;
            }
        }
    });

    let count = f64::from(grid * grid);
    assert_close(sum[0] / count, 0.0, 1.0e-14, "circular mean x");
    assert_close(sum[1] / count, 0.0, 1.0e-14, "circular mean y");
    assert_close(
        sum_radius_squared / count,
        0.5 * radius_m * radius_m,
        1.0e-4,
        "uniform disk E[r^2]",
    );
}

#[test]
fn polygon_aperture_is_uniform_contained_and_rotates_samplewise() {
    let radius_m = 0.8;
    let blades = 6_u8;
    let rotation = 0.19;
    let rotation_delta = 0.43;
    let camera = axial_camera(
        Aperture::try_regular_polygon(radius_m, blades, rotation).unwrap(),
        3.0,
    );
    let rotated_camera = axial_camera(
        Aperture::try_regular_polygon(radius_m, blades, rotation + rotation_delta).unwrap(),
        3.0,
    );
    let grid = 120_u32;
    let step = 2.0 * core::f64::consts::PI / f64::from(blades);
    let apothem = radius_m * (0.5 * step).cos();
    let (sin_delta, cos_delta) = rotation_delta.sin_cos();
    let mut sum = [0.0; 2];
    let mut sum_radius_squared = 0.0;

    with_cx(|_, cx| {
        for row in 0..grid {
            for column in 0..grid {
                let sample = LensSample::try_new(
                    (f64::from(column) + 0.5) / f64::from(grid),
                    (f64::from(row) + 0.5) / f64::from(grid),
                )
                .unwrap();
                let ray = camera
                    .generate_ray_from_tangent_offsets(cx, 0.0, 0.0, sample)
                    .unwrap();
                let rotated_ray = rotated_camera
                    .generate_ray_from_tangent_offsets(cx, 0.0, 0.0, sample)
                    .unwrap();
                let displacement = ray.origin.delta_from(camera.eye());
                assert_close(
                    displacement.dot(camera.forward()),
                    0.0,
                    2.0e-15,
                    "polygon aperture remains in lens plane",
                );
                let [x, y] = lens_coordinates(&camera, ray.origin);
                let [rotated_x, rotated_y] = lens_coordinates(&rotated_camera, rotated_ray.origin);

                assert_close(
                    rotated_x,
                    cos_delta * x - sin_delta * y,
                    3.0e-12,
                    "polygon sample rotation x",
                );
                assert_close(
                    rotated_y,
                    sin_delta * x + cos_delta * y,
                    3.0e-12,
                    "polygon sample rotation y",
                );
                for edge in 0..blades {
                    let normal_angle = rotation + (f64::from(edge) + 0.5) * step;
                    let support = x * normal_angle.cos() + y * normal_angle.sin();
                    assert!(
                        support <= apothem + 2.0e-12,
                        "polygon sample escaped edge {edge}: support={support}, apothem={apothem}"
                    );
                }
                sum[0] += x;
                sum[1] += y;
                sum_radius_squared += x * x + y * y;
            }
        }
    });

    let count = f64::from(grid * grid);
    let expected_moment = radius_m * radius_m * (2.0 + step.cos()) / 6.0;
    assert_close(sum[0] / count, 0.0, 2.0e-4, "polygon mean x");
    assert_close(sum[1] / count, 0.0, 2.0e-4, "polygon mean y");
    assert_close(
        sum_radius_squared / count,
        expected_moment,
        2.0e-3,
        "uniform regular-polygon E[r^2]",
    );
}

#[test]
fn nearly_collinear_look_at_is_refused_with_ranked_actionable_fixes() {
    let eye = Point3::new(0.0, 0.0, 0.0);
    let target = Point3::new(0.0, 0.0, -1.0);
    let result = PhysicalCamera::try_look_at(
        eye,
        target,
        Vec3::new(0.0, 1.0e-12, -1.0),
        canonical_projection(),
        1.0,
        pinhole(),
    );
    assert_eq!(result, Err(CameraError::NearlyCollinearUp));

    let fixes = CameraError::NearlyCollinearUp.ranked_fixes();
    assert_eq!(fixes.len(), 2);
    assert_eq!(
        (fixes[0].rank, fixes[0].code),
        (1, "choose_least_aligned_axis")
    );
    assert_eq!((fixes[1].rank, fixes[1].code), (2, "declare_roll"));
    assert!(fixes.iter().all(|fix| !fix.message.is_empty()));

    let suggested = PhysicalCamera::suggested_up_reference(target.delta_from(eye)).unwrap();
    assert_eq!(suggested, Vec3::new(1.0, 0.0, 0.0));
    PhysicalCamera::try_look_at(
        eye,
        target,
        suggested,
        canonical_projection(),
        1.0,
        pinhole(),
    )
    .expect("explicitly applying the ranked roll fix must admit the camera");

    PhysicalCamera::try_look_at(
        eye,
        target,
        Vec3::new(0.0, 1.0e-8, -1.0),
        canonical_projection(),
        1.0,
        pinhole(),
    )
    .expect("a non-collinear declared roll above the guard must be admitted");
    assert_eq!(
        PhysicalCamera::try_look_at(
            eye,
            eye,
            Vec3::new(0.0, 1.0, 0.0),
            canonical_projection(),
            1.0,
            pinhole(),
        ),
        Err(CameraError::DegenerateDirection)
    );
}

#[test]
fn keyframes_preserve_endpoints_and_interpolate_pose_and_focus_pull() {
    let left = look_at_camera(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        2.0,
        pinhole(),
    );
    let right = look_at_camera(
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        6.0,
        pinhole(),
    );
    let animation = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            7,
            0.0,
            2.0,
            vec![
                CameraKeyframe::try_new(0.0, left.clone()).unwrap(),
                CameraKeyframe::try_new(2.0, right.clone()).unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap();

    with_cx(|_, cx| {
        assert_eq!(animation.evaluate(cx, 0.0, CutSide::After).unwrap(), left);
        assert_eq!(animation.evaluate(cx, 2.0, CutSide::Before).unwrap(), right);
        let midpoint = animation.evaluate(cx, 1.0, CutSide::After).unwrap();
        assert_point_close(midpoint.eye(), Point3::new(1.0, 0.0, 0.0), "midpoint eye");
        assert_vec_close(
            midpoint.forward(),
            Vec3::new(-0.5_f64.sqrt(), 0.0, -0.5_f64.sqrt()),
            "shortest-path midpoint orientation",
        );
        assert_close(midpoint.focus_distance_m(), 4.0, TOLERANCE, "focus pull");
    });
}

#[test]
fn extreme_finite_keyframe_times_interpolate_without_overflow() {
    let left = axial_camera(pinhole(), 2.0);
    let right = PhysicalCamera::try_legacy_compatible(
        Point3::new(2.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        0.5,
        6.0,
        pinhole(),
    )
    .unwrap();
    let camera = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            71,
            -f64::MAX,
            f64::MAX,
            vec![
                CameraKeyframe::try_new(-f64::MAX, left).unwrap(),
                CameraKeyframe::try_new(f64::MAX, right).unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap();

    with_cx(|_, cx| {
        let midpoint = camera.evaluate(cx, 0.0, CutSide::After).unwrap();
        assert_point_close(
            midpoint.eye(),
            Point3::new(1.0, 0.0, 0.0),
            "overflow-safe midpoint eye",
        );
        assert_close(
            midpoint.focus_distance_m(),
            4.0,
            TOLERANCE,
            "overflow-safe focus pull",
        );
    });
}

#[test]
fn moving_world_focus_target_is_interpolated_then_projected_on_axis() {
    let left = axial_camera(pinhole(), 1.0);
    let right = PhysicalCamera::try_legacy_compatible(
        Point3::new(2.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        0.5,
        1.0,
        pinhole(),
    )
    .unwrap();
    let animation = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            8,
            0.0,
            2.0,
            vec![
                CameraKeyframe::try_with_world_focus(0.0, left, Point3::new(0.0, 0.0, -4.0))
                    .unwrap(),
                CameraKeyframe::try_with_world_focus(2.0, right, Point3::new(2.0, 0.0, -8.0))
                    .unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap();

    with_cx(|_, cx| {
        let left = animation.evaluate(cx, 0.0, CutSide::After).unwrap();
        let midpoint = animation.evaluate(cx, 1.0, CutSide::After).unwrap();
        let right = animation.evaluate(cx, 2.0, CutSide::Before).unwrap();
        assert_close(
            left.focus_distance_m(),
            4.0,
            TOLERANCE,
            "left tracked focus",
        );
        assert_close(
            midpoint.focus_distance_m(),
            6.0,
            TOLERANCE,
            "moving tracked focus",
        );
        assert_close(
            right.focus_distance_m(),
            8.0,
            TOLERANCE,
            "right tracked focus",
        );
    });
}

#[test]
fn hard_cuts_never_blend_and_crossing_or_out_of_range_queries_refuse() {
    let before = look_at_camera(
        Point3::new(-2.0, 0.0, 0.0),
        Point3::new(-2.0, 0.0, -1.0),
        3.0,
        pinhole(),
    );
    let after = look_at_camera(
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, -1.0),
        3.0,
        pinhole(),
    );
    let camera = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            1,
            0.0,
            1.0,
            vec![CameraKeyframe::try_new(0.0, before.clone()).unwrap()],
        )
        .unwrap(),
        CameraShot::try_new(
            2,
            1.0,
            2.0,
            vec![CameraKeyframe::try_new(1.0, after.clone()).unwrap()],
        )
        .unwrap(),
    ])
    .unwrap();

    with_cx(|_, cx| {
        assert_eq!(camera.evaluate(cx, 1.0, CutSide::Before).unwrap(), before);
        assert_eq!(camera.evaluate(cx, 1.0, CutSide::After).unwrap(), after);
        assert_eq!(
            camera.evaluate(cx, -f64::EPSILON, CutSide::After),
            Err(CameraError::Extrapolation)
        );
        assert_eq!(
            camera.evaluate(cx, 2.0 + 4.0 * f64::EPSILON, CutSide::Before),
            Err(CameraError::Extrapolation)
        );

        assert!(matches!(
            camera.admit_shutter(cx, shutter(0.75, 1.25), CutSide::After),
            Err(CameraError::ShutterCrossesCut)
        ));
        let before_exposure = camera
            .admit_shutter(cx, shutter(0.5, 1.0), CutSide::After)
            .unwrap();
        assert_eq!(before_exposure.shot_id(), 1);
        assert_eq!(
            camera.evaluate_exposure(cx, before_exposure, 1.0).unwrap(),
            before
        );
        assert_eq!(
            camera.evaluate_exposure(cx, before_exposure, 0.25),
            Err(CameraError::InvalidExposure)
        );

        let instant = shutter(1.0, 1.0);
        assert_eq!(
            camera
                .admit_shutter(cx, instant, CutSide::Before)
                .unwrap()
                .shot_id(),
            1
        );
        assert_eq!(
            camera
                .admit_shutter(cx, instant, CutSide::After)
                .unwrap()
                .shot_id(),
            2
        );
    });
}

#[test]
fn timeline_admission_rejects_three_shots_sharing_one_instant() {
    let physical = axial_camera(pinhole(), 3.0);
    let shot = |shot_id, start_s, end_s| {
        CameraShot::try_new(
            shot_id,
            start_s,
            end_s,
            vec![CameraKeyframe::try_new(start_s, physical.clone()).unwrap()],
        )
        .unwrap()
    };

    assert_eq!(
        AnimatedCamera::try_new(vec![
            shot(81, 0.0, 1.0),
            shot(82, 1.0, 1.0),
            shot(83, 1.0, 2.0),
        ]),
        Err(CameraError::InvalidShot)
    );

    AnimatedCamera::try_new(vec![
        shot(84, 0.0, 0.5),
        shot(85, 1.0, 1.0),
        shot(86, 1.5, 2.0),
    ])
    .expect("an isolated zero-duration shot is unambiguous and remains valid");
    AnimatedCamera::try_new(vec![
        shot(87, 0.0, 1.0),
        shot(88, 1.0, 2.0),
        shot(89, 2.0, 3.0),
    ])
    .expect("ordinary contiguous shots share only pairwise cut instants");

    let empty_fixes = CameraError::EmptyShot.ranked_fixes();
    assert_eq!(
        (empty_fixes[0].rank, empty_fixes[0].code),
        (1, "add_camera_keyframe")
    );
    let shot_fixes = CameraError::InvalidShot.ranked_fixes();
    assert_eq!(
        (shot_fixes[0].rank, shot_fixes[0].code),
        (1, "unique_shot_ids")
    );
    assert_eq!(
        (shot_fixes[1].rank, shot_fixes[1].code),
        (2, "partition_shot_timeline")
    );
}

#[test]
fn camera_ray_and_timeline_evaluation_observe_cancellation() {
    let physical = axial_camera(Aperture::try_circular(0.1).unwrap(), 3.0);
    with_cx(|gate, cx| {
        gate.request();
        assert_eq!(
            physical.generate_ray_from_tangent_offsets(cx, 0.1, -0.2, LensSample::CENTER),
            Err(CameraError::Cancelled)
        );
    });

    let animated = AnimatedCamera::try_static(9, 0.0, 2.0, physical).unwrap();
    with_cx(|gate, cx| {
        let exposure = animated
            .admit_shutter(cx, shutter(0.25, 0.75), CutSide::After)
            .unwrap();
        gate.request();
        assert_eq!(
            animated.evaluate(cx, 0.5, CutSide::After),
            Err(CameraError::Cancelled)
        );
        assert_eq!(
            animated.evaluate_exposure(cx, exposure, 0.5),
            Err(CameraError::Cancelled)
        );
    });
    with_cx(|gate, cx| {
        gate.request();
        assert!(matches!(
            animated.admit_shutter(cx, shutter(0.25, 0.75), CutSide::After),
            Err(CameraError::Cancelled)
        ));
    });
}

fn diffraction_camera(radius_m: f64, focus_m: f64) -> PhysicalCamera {
    axial_camera(Aperture::try_circular(radius_m).unwrap(), focus_m)
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.05).unwrap())
        .unwrap()
}

/// G2: the first circular-pupil null is at 1.21966989 lambda/D, and packet
/// wavelengths retain different ring positions even though they share a ray.
#[test]
fn diffraction_airy_zero_scale_and_packet_chromaticity() {
    let zero = 3.831_705_970_207_512_5_f64;
    let u = 1.0 - 1.0 / (1.0 + (zero / 2.0).powi(2)).sqrt();
    let sample = LensSample::try_new(u, 0.0).unwrap();
    let radius = 1.0e-3;
    with_cx(|_, cx| {
        let ray = diffraction_camera(radius, 2.0)
            .generate_spectral_ray_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                sample,
                &[400.0, 550.0, 700.0],
            )
            .unwrap();
        let tangent = ray.ray.dir.x / -ray.ray.dir.z;
        assert_close(
            tangent,
            zero * 550.0e-9 / (2.0 * core::f64::consts::PI * radius),
            1.0e-16,
            "Airy null angular radius",
        );
        assert!(
            ray.weights[1] < 1.0e-25,
            "central wavelength misses Airy null: {:?}",
            ray.weights
        );
        assert!(ray.weights[0] > 0.01 && ray.weights[2] > 0.01);
        assert_eq!(ray.ray.origin, Point3::new(0.0, 0.0, 0.0));
        let wider = diffraction_camera(2.0 * radius, 2.0)
            .generate_spectral_ray_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                sample,
                &[400.0, 550.0, 700.0],
            )
            .unwrap();
        assert_close(
            wider.ray.dir.x / -wider.ray.dir.z,
            0.5 * tangent,
            1.0e-16,
            "inverse pupil scale",
        );
        assert_eq!(wider.weights, ray.weights);
        let red = diffraction_camera(radius, 2.0)
            .generate_spectral_ray_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                sample,
                &[800.0, 1100.0, 1400.0],
            )
            .unwrap();
        assert_close(
            red.ray.dir.x / -red.ray.dir.z,
            2.0 * tangent,
            1.0e-16,
            "wavelength scale",
        );
        assert_eq!(red.weights, ray.weights);
    });
}

/// G2: integrate the actual wavelength-weighted sampling map. The Airy disk
/// contains 0.837784869... of the power; the rings retain the rest, so a
/// constant radiance field is preserved rather than clipped or peak-normalized.
#[test]
fn diffraction_packet_conserves_radiance_and_encircled_power() {
    let camera = diffraction_camera(1.0e-3, 2.0);
    let wavelengths = [400.0, 550.0, 700.0];
    let n = 8_192;
    with_cx(|_, cx| {
        let mut total = [0.0; 3];
        let mut inner = [0.0; 3];
        for i in 0..n {
            let u = (f64::from(i) + 0.5) / f64::from(n);
            let ray = camera
                .generate_spectral_ray_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    LensSample::try_new(u, 0.0).unwrap(),
                    &wavelengths,
                )
                .unwrap();
            let theta = ray.ray.dir.x / -ray.ray.dir.z;
            for lane in 0..3 {
                total[lane] += ray.weights[lane] / f64::from(n);
                let v = 2.0 * core::f64::consts::PI * 1.0e-3 * theta / (wavelengths[lane] * 1.0e-9);
                if v < 3.831_705_970_207_512_5 {
                    inner[lane] += ray.weights[lane] / f64::from(n);
                }
            }
        }
        for lane in 0..3 {
            assert_close(total[lane], 1.0, 2.0e-3, "flat radiance conservation");
            assert_close(
                inner[lane],
                0.837_784_869_173_314_3,
                1.0e-7,
                "Airy disk power",
            );
        }
    });
}

#[test]
fn diffraction_refuses_defocus_invalid_pupils_and_nonspectral_paths() {
    for value in [-1.0, 0.11, f64::NAN, f64::INFINITY] {
        assert!(FocusedDiffraction::try_new(value).is_err());
    }
    let policy = FocusedDiffraction::try_new(0.05).unwrap();
    let virtual_image = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        0.025,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap();
    assert!(virtual_image.with_in_focus_diffraction(policy).is_err());
    assert!(
        axial_camera(pinhole(), 2.0)
            .with_in_focus_diffraction(policy)
            .is_err()
    );
    assert!(
        axial_camera(Aperture::try_regular_polygon(0.01, 6, 0.0).unwrap(), 2.0)
            .with_in_focus_diffraction(policy)
            .is_ok()
    );
    let camera = diffraction_camera(1.0e-3, 2.0);
    assert!(
        camera
            .admit_diffraction_depth(Some(Point3::new(0.1, 0.0, -2.0)), &[400.0, 700.0])
            .is_ok()
    );
    assert_eq!(
        camera.admit_diffraction_depth(Some(Point3::new(0.0, 0.0, -1.0)), &[400.0]),
        Err(CameraError::DefocusOutsideDiffractionModel)
    );
    assert_eq!(
        camera.admit_diffraction_depth(None, &[400.0]),
        Err(CameraError::DefocusOutsideDiffractionModel)
    );
    assert!(
        diffraction_camera(1.0e-3, 1.0e8)
            .admit_diffraction_depth(None, &[400.0])
            .is_ok()
    );
    let depth = 1.0 / (0.5 - 0.04 * 700.0e-9 / (core::f64::consts::PI * 1.0e-6));
    let point = Some(Point3::new(0.0, 0.0, -depth));
    assert!(camera.admit_diffraction_depth(point, &[700.0]).is_ok());
    assert_eq!(
        camera.admit_diffraction_depth(point, &[700.0, 400.0]),
        Err(CameraError::DefocusOutsideDiffractionModel)
    );
    with_cx(|_, cx| {
        assert_eq!(
            camera.generate_ray_from_tangent_offsets(cx, 0.0, 0.0, LensSample::CENTER),
            Err(CameraError::DiffractionRequiresSpectralPath)
        );
        for wavelength in [0.0, -500.0, f64::NAN, f64::INFINITY, 1.0e9] {
            assert!(
                camera
                    .generate_spectral_ray_from_tangent_offsets(
                        cx,
                        0.0,
                        0.0,
                        LensSample::CENTER,
                        LensSample::CENTER,
                        &[wavelength],
                    )
                    .is_err()
            );
        }
        let geometric = axial_camera(Aperture::try_circular(0.01).unwrap(), 2.0);
        let old = geometric
            .generate_ray_from_tangent_offsets(cx, 0.1, -0.2, LensSample::CENTER)
            .unwrap();
        let spectral = geometric
            .generate_spectral_ray_from_tangent_offsets(
                cx,
                0.1,
                -0.2,
                LensSample::CENTER,
                LensSample::CENTER,
                &[400.0, 700.0],
            )
            .unwrap();
        assert_eq!(spectral.ray, old);
        assert_eq!(spectral.weights, [1.0; 2]);
    });
    with_cx(|gate, cx| {
        gate.request();
        assert!(matches!(
            camera.generate_spectral_ray_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                LensSample::CENTER,
                &[550.0],
            ),
            Err(CameraError::Cancelled)
        ));
    });
}

#[test]
fn diffraction_policy_survives_camera_animation() {
    let first = diffraction_camera(1.0e-3, 2.0);
    let last = diffraction_camera(1.0e-3, 3.0);
    let animated = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            74,
            0.0,
            1.0,
            vec![
                CameraKeyframe::try_new(0.0, first).unwrap(),
                CameraKeyframe::try_new(1.0, last).unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap();
    with_cx(|_, cx| {
        let camera = animated.evaluate(cx, 0.5, CutSide::After).unwrap();
        assert_eq!(
            camera.in_focus_diffraction(),
            Some(FocusedDiffraction::try_new(0.05).unwrap())
        );
        assert!(
            camera
                .admit_diffraction_depth(Some(Point3::new(0.0, 0.0, -2.5)), &[550.0])
                .is_ok()
        );
    });
}

/// A square pupil has a separable sinc² image, not circular Airy rings.
/// The central density is its actual area/lambda² and scales accordingly.
#[test]
fn bladed_diffraction_density_has_physical_area_and_spectral_scaling() {
    let radius = 0.001;
    let wavelength_nm = 550.0;
    let lambda = wavelength_nm * 1.0e-9;
    let camera = |radius, rotation| {
        axial_camera(
            Aperture::try_regular_polygon(radius, 4, rotation).unwrap(),
            2.0,
        )
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.05).unwrap())
        .unwrap()
    };
    let square = camera(radius, core::f64::consts::FRAC_PI_4);
    let peak = square
        .in_focus_diffraction_density([0.0; 2], wavelength_nm)
        .unwrap();
    assert_close(
        peak,
        2.0 * radius * radius / (lambda * lambda),
        1.0e-8,
        "physical square area",
    );
    let first_zero = lambda / (core::f64::consts::SQRT_2 * radius);
    assert!(
        square
            .in_focus_diffraction_density([first_zero, 0.0], wavelength_nm)
            .unwrap()
            < 1.0e-20 * peak
    );
    let offset = [0.37 * first_zero, -0.61 * first_zero];
    let density = square
        .in_focus_diffraction_density(offset, wavelength_nm)
        .unwrap();
    let red = square
        .in_focus_diffraction_density(offset.map(|v| 2.0 * v), 2.0 * wavelength_nm)
        .unwrap();
    assert_close(
        red,
        density / 4.0,
        1.0e-12 * peak,
        "inverse-square spectral density",
    );
    let wider = camera(2.0 * radius, core::f64::consts::FRAC_PI_4)
        .in_focus_diffraction_density(offset.map(|v| 0.5 * v), wavelength_nm)
        .unwrap();
    assert_close(
        wider,
        4.0 * density,
        1.0e-12 * peak,
        "pupil-area density scale",
    );
    let rotated = camera(radius, core::f64::consts::FRAC_PI_2);
    let angle = core::f64::consts::FRAC_PI_4;
    let rotated_offset = [
        offset[0] * angle.cos() - offset[1] * angle.sin(),
        offset[0] * angle.sin() + offset[1] * angle.cos(),
    ];
    assert_close(
        rotated
            .in_focus_diffraction_density(rotated_offset, wavelength_nm)
            .unwrap(),
        density,
        1.0e-12 * peak,
        "rotated diffraction spikes",
    );
    assert!(
        axial_camera(pinhole(), 2.0)
            .in_focus_diffraction_density([0.0; 2], wavelength_nm)
            .is_err()
    );
    assert!(
        square
            .in_focus_diffraction_density([f64::NAN, 0.0], wavelength_nm)
            .is_err()
    );
    let small_triangle = axial_camera(Aperture::try_regular_polygon(4.0e-5, 3, 0.0).unwrap(), 2.0)
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.05).unwrap())
        .unwrap();
    assert!(
        small_triangle
            .in_focus_diffraction_density([0.0; 2], wavelength_nm)
            .is_err(),
        "scalar admission must resolve the inscribed pupil diameter"
    );
}

#[test]
fn bladed_diffraction_sampling_rotates_with_the_actual_aperture() {
    let camera = |rotation| {
        axial_camera(
            Aperture::try_regular_polygon(0.001, 7, rotation).unwrap(),
            2.0,
        )
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.05).unwrap())
        .unwrap()
    };
    let (plain, rotated) = (camera(0.0), camera(0.4));
    with_cx(|_, cx| {
        for u in [0.0, 0.137, 0.37, 0.73, 0.999] {
            let sample = LensSample::try_new(u, 0.31).unwrap();
            let draw = |camera: &PhysicalCamera| {
                camera
                    .generate_spectral_ray_from_tangent_offsets(
                        cx,
                        0.0,
                        0.0,
                        LensSample::CENTER,
                        sample,
                        &[400.0, 550.0, 700.0],
                    )
                    .unwrap()
            };
            let a = draw(&plain);
            let b = draw(&rotated);
            let x = a.ray.dir.x;
            let y = a.ray.dir.y;
            assert_close(
                b.ray.dir.x,
                x * 0.4_f64.cos() - y * 0.4_f64.sin(),
                1.0e-12,
                "rotated sampled x",
            );
            assert_close(
                b.ray.dir.y,
                x * 0.4_f64.sin() + y * 0.4_f64.cos(),
                1.0e-12,
                "rotated sampled y",
            );
            assert_close(b.ray.dir.z, a.ray.dir.z, 1.0e-12, "rotated sampled z");
            if u > 0.0 {
                for (left, right) in a.weights.into_iter().zip(b.weights) {
                    assert_close(
                        left,
                        right,
                        1.0e-8 * left.max(1.0),
                        "rotation-invariant spectral weights",
                    );
                }
            }
        }
    });
}

#[test]
fn defocused_camera_preserves_physical_phase_and_spectral_psf() {
    let radius = 1.0e-3;
    let wavelength = 550.0;
    let lambda_m = wavelength * 1.0e-9;
    let focus = 0.25;
    let inverse_shift = 2.0 * lambda_m / (radius * radius);
    let object_depth = 1.0 / (1.0 / focus + inverse_shift);
    let policy = DefocusedDiffraction::try_new(Some(object_depth), 0.0).unwrap();
    let camera = axial_camera(Aperture::try_circular(radius).unwrap(), focus)
        .with_defocused_diffraction(policy)
        .unwrap();
    assert!(camera.has_diffraction());
    assert_eq!(camera.defocused_diffraction(), Some(policy));
    assert_eq!(camera.in_focus_diffraction(), None);
    let focused_peak = core::f64::consts::PI * radius * radius / (lambda_m * lambda_m);
    assert!(
        camera
            .defocused_diffraction_density([0.0; 2], wavelength)
            .unwrap()
            < 1.0e-25 * focused_peak,
        "two-pi rim phase has an axial null, not a broadened positive Airy core"
    );
    for wavelength_nm in [400.0, 700.0] {
        let phase = 2.0 * core::f64::consts::PI * wavelength / wavelength_nm;
        let expected = core::f64::consts::PI * radius * radius / (wavelength_nm * 1.0e-9).powi(2)
            * ((0.5 * phase).sin() / (0.5 * phase)).powi(2);
        assert_close(
            camera
                .defocused_diffraction_density([0.0; 2], wavelength_nm)
                .unwrap(),
            expected,
            1.0e-12 * focused_peak,
            "physical on-axis defocus and chromatic phase",
        );
    }
    let opposite = axial_camera(Aperture::try_circular(radius).unwrap(), focus)
        .with_defocused_diffraction(
            DefocusedDiffraction::try_new(Some(1.0 / (1.0 / focus - inverse_shift)), 0.0).unwrap(),
        )
        .unwrap();
    let offset = [2.0e-4, -1.0e-4];
    let density = camera
        .defocused_diffraction_density(offset, wavelength)
        .unwrap();
    assert_close(
        opposite
            .defocused_diffraction_density(offset, wavelength)
            .unwrap(),
        density,
        1.0e-12 * focused_peak,
        "opposite wavefront curvature has the same circular-pupil intensity",
    );
    let scaled = axial_camera(Aperture::try_circular(2.0 * radius).unwrap(), 2.0 * focus)
        .with_defocused_diffraction(
            DefocusedDiffraction::try_new(Some(2.0 * object_depth), 0.0).unwrap(),
        )
        .unwrap();
    assert_close(
        scaled
            .defocused_diffraction_density(offset, 2.0 * wavelength)
            .unwrap(),
        density,
        1.0e-12 * focused_peak,
        "common SI-length and wavelength scaling preserves angular imaging",
    );
    assert!(
        camera
            .admit_diffraction_depth(Some(Point3::new(0.1, 0.2, -object_depth)), &[400.0])
            .is_ok()
    );
    assert_eq!(
        camera.admit_diffraction_depth(Some(Point3::new(0.0, 0.0, -focus)), &[550.0]),
        Err(CameraError::DefocusOutsideDiffractionModel)
    );
    assert_eq!(
        camera.admit_diffraction_depth(None, &[550.0]),
        Err(CameraError::DefocusOutsideDiffractionModel)
    );
    with_cx(|_, cx| {
        let first = camera
            .generate_spectral_ray_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                LensSample::try_new(0.0, 0.0).unwrap(),
                &[400.0, 550.0, 700.0],
            )
            .unwrap();
        assert!(first.weights[1] < 1.0e-24);
        assert!(first.weights[0] > 0.01 && first.weights[2] > 0.01);
        for u in [0.37, 1.0 - 0.5 * f64::EPSILON] {
            let sample = LensSample::try_new(u, 0.23).unwrap();
            let draw = |lens_sample| {
                camera
                    .generate_spectral_ray_from_tangent_offsets(
                        cx,
                        0.0,
                        0.0,
                        lens_sample,
                        sample,
                        &[400.0, 550.0, 700.0],
                    )
                    .unwrap()
            };
            let a = draw(LensSample::CENTER);
            let b = draw(LensSample::try_new(0.23, 0.71).unwrap());
            assert_eq!(a.ray.origin, camera.eye());
            assert_eq!(a.ray.dir, b.ray.dir);
            assert_eq!(a.weights, b.weights, "no extra geometric pupil defocus");
            assert!(a.weights.iter().all(|v| v.is_finite() && *v >= 0.0));
        }
        assert_eq!(
            camera.generate_ray_from_tangent_offsets(cx, 0.0, 0.0, LensSample::CENTER),
            Err(CameraError::DiffractionRequiresSpectralPath)
        );
    });
}

#[test]
fn finite_defocus_admission_and_focus_rack_preserve_the_declared_object_plane() {
    for depth in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0e-320] {
        assert!(DefocusedDiffraction::try_new(Some(depth), 0.01).is_err());
    }
    for error in [-0.01, 0.10001, f64::NAN] {
        assert!(DefocusedDiffraction::try_new(None, error).is_err());
    }
    let policy = DefocusedDiffraction::try_new(Some(0.5), 0.01).unwrap();
    assert!(
        axial_camera(pinhole(), 0.5)
            .with_defocused_diffraction(policy)
            .is_err()
    );
    assert!(
        axial_camera(Aperture::try_regular_polygon(0.001, 6, 0.0).unwrap(), 0.5)
            .with_defocused_diffraction(policy)
            .is_err()
    );
    let camera = |focus| {
        axial_camera(Aperture::try_circular(0.001).unwrap(), focus)
            .with_defocused_diffraction(policy)
            .unwrap()
    };
    let animation = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            902,
            0.0,
            1.0,
            vec![
                CameraKeyframe::try_new(0.0, camera(0.5)).unwrap(),
                CameraKeyframe::try_new(1.0, camera(1.0)).unwrap(),
            ],
        )
        .unwrap(),
    ])
    .unwrap();
    with_cx(|_, cx| {
        let middle = animation.evaluate(cx, 0.5, CutSide::After).unwrap();
        assert_eq!(middle.defocused_diffraction(), Some(policy));
        assert_eq!(middle.focus_distance_m(), 0.75);
        assert!(
            middle
                .admit_diffraction_depth(Some(Point3::new(0.0, 0.0, -0.5)), &[400.0, 700.0])
                .is_ok()
        );
        let focused = middle
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.05).unwrap())
            .unwrap();
        assert_eq!(focused.defocused_diffraction(), None);
        assert!(focused.in_focus_diffraction().is_some());
    });
    let beyond = axial_camera(Aperture::try_circular(0.001).unwrap(), 2.0)
        .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(0.3), 0.01).unwrap())
        .unwrap();
    assert!(
        beyond
            .defocused_diffraction_density([0.0; 2], 700.0)
            .is_ok()
    );
    assert_eq!(
        beyond.defocused_diffraction_density([0.0; 2], 400.0),
        Err(CameraError::DefocusOutsideDiffractionModel)
    );
    let infinity = axial_camera(Aperture::try_circular(0.001).unwrap(), 1.0)
        .with_defocused_diffraction(DefocusedDiffraction::try_new(None, 0.0).unwrap())
        .unwrap();
    assert!(
        infinity
            .admit_diffraction_depth(None, &[400.0, 700.0])
            .is_ok()
    );
    assert!(
        infinity
            .defocused_diffraction_density([0.0; 2], 550.0)
            .unwrap()
            > 0.0
    );
}
