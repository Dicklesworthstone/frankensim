//! G1: exact refracting lens paths, thick-lens limits and sensor radiometry.
#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, LensSample, LensSurface, PhysicalCamera,
    SequentialLens, SequentialLensRay,
};
use fs_render::charts::Ray;
use fs_render::conductor::{
    CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use fs_render::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use fs_render::dielectric::{DielectricCoating, DielectricSurface};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 13,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

fn glass(ior: CauchyIor, extinction: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(extinction).unwrap(),
        GlassProvenance::Custom,
    )
}

fn surface(z: f64, curvature: f64, radius: f64, after: Option<DielectricGlass>) -> LensSurface {
    LensSurface::try_new(z, curvature, radius, after).unwrap()
}

fn lens(surfaces: &[LensSurface]) -> SequentialLens {
    SequentialLens::try_new(550.0, surfaces).unwrap()
}

fn biconvex(material: DielectricGlass) -> SequentialLens {
    lens(&[
        surface(0.08, 20.0, 0.012, Some(material)),
        surface(0.09, -20.0, 0.012, None),
    ])
}

fn trace(lens: &SequentialLens, cx: &Cx<'_>, ray: Ray, nm: f64) -> SequentialLensRay {
    lens.trace_ray(cx, ray, nm).unwrap().unwrap()
}

fn ray(sensor_x: f64, slope: f64) -> Ray {
    Ray {
        origin: Point3::new(sensor_x, 0.0, 0.0),
        dir: Vec3::new(slope, 0.0, 1.0),
    }
}

fn axial_intercept(ray: Ray) -> f64 {
    ray.origin.z - ray.origin.x * ray.dir.z / ray.dir.x
}

fn coating(layers: &[(f64, f64, f64)]) -> DielectricCoating {
    let layers: Vec<_> = layers
        .iter()
        .map(|&(n, k, thickness)| {
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(n).unwrap(), k, thickness).unwrap()
        })
        .collect();
    let source = ConductorSource::try_new(
        fs_blake3::hash_domain("org.frankensim.test.sequential-lens", b"authored coating"),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&layers, 1, source).unwrap())
        .unwrap()
}

#[test]
fn quarter_wave_lens_coatings_remove_fresnel_loss_without_changing_refraction() {
    let index = 1.5_f64;
    let film_index = index.sqrt();
    let ar = coating(&[(film_index, 0.0, 550e-9 / (4.0 * film_index))]);
    let uncoated = biconvex(glass(CauchyIor::try_constant(index).unwrap(), 7.0));
    let coated = lens(&[
        uncoated.surfaces()[0].with_coating(ar),
        uncoated.surfaces()[1].with_coating(ar),
    ]);
    assert_eq!(coated.surfaces()[0].coating(), Some(&ar));
    assert_eq!(coated.reference_focus(), uncoated.reference_focus());
    with_cx(|_, cx| {
        let bare = trace(&uncoated, cx, ray(0.0, 0.0), 550.0);
        let coated_axis = trace(&coated, cx, ray(0.0, 0.0), 550.0);
        // Quarter-wave n_f=sqrt(n_air*n_glass) cancels both normal reflections.
        let expected = (-7.0_f64 * 0.01).exp();
        close(coated_axis.transmittance, expected, 8e-14);
        close(bare.transmittance, 0.96_f64.powi(2) * expected, 8e-14);
        for slope in [0.0, 0.015, 0.075] {
            for wavelength in [450.0, 550.0, 650.0] {
                let bare = trace(&uncoated, cx, ray(0.0, slope), wavelength);
                let coated = trace(&coated, cx, ray(0.0, slope), wavelength);
                assert_eq!(coated.ray, bare.ray);
                assert_eq!(coated.distance_m, bare.distance_m);
                assert_eq!(coated.optical_path_m, bare.optical_path_m);
                assert!(coated.transmittance > bare.transmittance);
                let glass_path = (coated.optical_path_m - coated.distance_m) / (index - 1.0);
                assert!(coated.transmittance <= (-7.0 * glass_path).exp() + 1e-13);
            }
        }
        assert!(trace(&coated, cx, ray(0.0, 0.0), 450.0).transmittance < expected - 1e-3);
    });
}

#[test]
fn absorbing_lens_films_retain_authored_order_and_polarization_correlation() {
    let stack = coating(&[(1.3, 0.08, 120e-9), (2.2, 0.0, 80e-9)]);
    let material = glass(CauchyIor::try_constant(1.5).unwrap(), 0.7);
    let coated = lens(&[
        surface(0.02, 0.0, 0.1, Some(material)).with_coating(stack),
        surface(0.03, 0.0, 0.1, None).with_coating(stack),
    ]);
    with_cx(|_, cx| {
        let si = 0.8_f64;
        let ci = (1.0 - si * si).sqrt();
        let ct = (1.0 - (si / 1.5).powi(2)).sqrt();
        let film = DielectricSurface::SMOOTH.with_coating(stack);
        let response = film.interface(550.0, false).unwrap();
        let first = response.fresnel_polarized(ci, 1.0, 1.5).unwrap();
        let second = response.fresnel_polarized(ct, 1.5, 1.0).unwrap();
        let one = fs_render::polarization::Stokes::unpolarized(1.0);
        let expected = second
            .transmission
            .apply_diagonal_jones(first.transmission.apply_diagonal_jones(one).unwrap())
            .unwrap()
            .i
            * (-0.7_f64 * 0.01 / ct).exp();
        let traced = trace(&coated, cx, ray(0.0, si / ci), 550.0);
        close(traced.transmittance, expected, 5e-13);
        close(traced.ray.dir.x, si, 3e-14);
        assert!(first.scalar.absorption_residual() > 0.05);
        let independent_averages =
            first.scalar.transmittance * second.scalar.transmittance * (-0.7_f64 * 0.01 / ct).exp();
        assert!((expected - independent_averages).abs() > 1e-3);
        // The front face is authored sensor-first too, not implicitly reversed.
        let wrong_order = film
            .interface(550.0, true)
            .unwrap()
            .fresnel_polarized(ct, 1.5, 1.0)
            .unwrap()
            .transmission
            .apply_diagonal_jones(first.transmission.apply_diagonal_jones(one).unwrap())
            .unwrap()
            .i
            * (-0.7_f64 * 0.01 / ct).exp();
        assert!((expected - wrong_order).abs() > 1e-3);
    });
}

#[test]
fn lens_coatings_admit_source_band_and_refuse_unmodeled_group_delay() {
    let zero = coating(&[(1.3, 0.1, 0.0)]);
    let active = coating(&[(1.3, 0.0, 100e-9)]);
    let stop = surface(0.02, 0.0, 0.01, None).with_coating(zero);
    assert!(SequentialLens::try_new(1000.0, &[stop]).is_err());
    let stopped = lens(&[stop]);
    assert_eq!(stopped.wavelength_bounds_nm(), (380.0, 780.0));
    let bare = biconvex(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0));
    let camera = |coating| {
        PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
            1.0,
            Aperture::try_circular(0.001).unwrap(),
        )
        .unwrap()
        .with_sequential_lens(lens(&[
            bare.surfaces()[0].with_coating(coating),
            bare.surfaces()[1],
        ]))
        .unwrap()
    };
    with_cx(|_, cx| {
        assert_eq!(
            stopped.trace_ray(cx, ray(1.0, 0.0), 1000.0),
            Err(CameraError::SequentialLensOutsideModel)
        );
        let zero_lens = camera(zero);
        for slope in [0.0, 0.015, 0.075] {
            assert_eq!(
                trace(
                    zero_lens.sequential_lens().unwrap(),
                    cx,
                    ray(0.0, slope),
                    550.0
                ),
                trace(&bare, cx, ray(0.0, slope), 550.0)
            );
        }
        let timed = zero_lens
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert!(timed.group_delay_s().is_some());
        let active_lens = camera(active);
        assert!(
            active_lens
                .generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                )
                .unwrap()[0]
                .is_some()
        );
        for x in [0.0, 20.0] {
            assert_eq!(
                active_lens.generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    x,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                ),
                Err(CameraError::SequentialLensOutsideModel)
            );
        }
    });
}

// Thick-lens power P=P1+P2-d*P1*P2/n; reduced-angle ABCD imaging.
// Callegari, https://arxiv.org/abs/2104.00700, sections II–III.
fn paraxial(n: f64) -> (f64, f64) {
    let p = (n - 1.0) / 0.05;
    let power = 2.0 * p - 0.01 * p * p / n;
    let a = 1.0 - 0.01 * p / n;
    (
        1.0 / power,
        0.09 + (a * 0.08 + 0.01 / n) / (power * 0.08 - a),
    )
}

#[test]
fn slab_retains_polarization_and_uses_oblique_bulk_path() {
    let n = 1.5;
    let sigma = 7.0;
    let material = glass(CauchyIor::try_constant(n).unwrap(), sigma);
    let lens = lens(&[
        surface(0.02, 0.0, 0.05, Some(material)),
        surface(0.03, 0.0, 0.05, None),
    ]);
    with_cx(|_, cx| {
        for si in [0.0_f64, 0.8] {
            let ci = (1.0 - si * si).sqrt();
            let st = si / n;
            let ct = (1.0 - st * st).sqrt();
            let traced = trace(&lens, cx, ray(0.0, si / ci), 550.0);
            let rs = ((ci - n * ct) / (ci + n * ct)).powi(2);
            let rp = ((n * ci - ct) / (n * ci + ct)).powi(2);
            let beer = (-sigma * 0.01 / ct).exp();
            let expected = 0.5 * ((1.0 - rs).powi(2) + (1.0 - rp).powi(2)) * beer;
            close(traced.transmittance, expected, 3e-13);
            close(traced.distance_m, 0.02 / ci + 0.01 / ct, 2e-14);
            close(traced.optical_path_m, 0.02 / ci + n * 0.01 / ct, 2e-14);
            close(traced.ray.origin.x, 0.02 * si / ci + 0.01 * st / ct, 2e-14);
            close(traced.ray.origin.z, 0.03, 2e-14);
            close(traced.ray.dir.x, si, 2e-14);
            close(traced.ray.dir.z, ci, 2e-14);
            if si != 0.0 {
                let depolarized_product = (1.0 - 0.5 * (rs + rp)).powi(2) * beer;
                assert!(expected - depolarized_product > 1e-3);
            }
        }
    });
}

#[test]
fn thick_lens_focus_matches_lensmaker_and_near_axis_exact_ray() {
    let lens = biconvex(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0));
    let expected = paraxial(1.5);
    let focus = lens.reference_focus().unwrap();
    close(focus.0, expected.0, 2e-14);
    close(focus.1, expected.1, 2e-14);
    assert!(
        (focus.0 - 0.05).abs() > 1e-3,
        "thickness changes focal power"
    );
    with_cx(|_, cx| {
        let exact = trace(&lens, cx, ray(0.0, 1e-6), 550.0);
        close(axial_intercept(exact.ray), expected.1, 3e-10);
    });
}

#[test]
fn exact_lens_resolves_chromatic_focus_and_spherical_aberration() {
    let ior = CauchyIor::try_new(1.48, 0.006, 0.0).unwrap();
    let lens = biconvex(glass(ior, 0.0));
    with_cx(|_, cx| {
        let focus = [450.0, 650.0].map(|wavelength| {
            let near = trace(&lens, cx, ray(0.0, 1e-6), wavelength);
            let near_focus = axial_intercept(near.ray);
            close(near_focus, paraxial(ior.eval(wavelength).unwrap()).1, 3e-10);
            let marginal = trace(&lens, cx, ray(0.0, 0.1), wavelength);
            assert!(axial_intercept(marginal.ray) < near_focus - 0.005);
            near_focus
        });
        assert!(
            focus[0] < focus[1] - 0.005,
            "blue converges sooner than red"
        );
    });
}

#[test]
fn internal_stop_preserves_passed_rays_and_curved_exit_can_totally_reflect() {
    let material = glass(CauchyIor::try_constant(1.5).unwrap(), 0.0);
    let bare = biconvex(material);
    let stopped = lens(&[
        bare.surfaces()[0],
        surface(0.085, 0.0, 0.001, Some(material)),
        bare.surfaces()[1],
    ]);
    let tir = lens(&[
        surface(0.02, 0.0, 0.009, Some(material)),
        surface(0.04, -100.0, 0.009, None),
    ]);
    with_cx(|_, cx| {
        let plain = trace(&bare, cx, ray(0.0, 0.005), 550.0);
        let passed = trace(&stopped, cx, ray(0.0, 0.005), 550.0);
        close(passed.transmittance, plain.transmittance, 3e-13);
        close(passed.distance_m, plain.distance_m, 2e-14);
        close(passed.optical_path_m, plain.optical_path_m, 2e-14);
        close(passed.ray.origin.x, plain.ray.origin.x, 2e-14);
        close(passed.ray.dir.x, plain.ray.dir.x, 2e-14);
        assert!(
            stopped
                .trace_ray(cx, ray(0.0, 0.05), 550.0)
                .unwrap()
                .is_none()
        );
        // At this exit sphere, sin(theta_i)=|curvature*x|; n*sin(theta_i)>1.
        assert!(tir.trace_ray(cx, ray(0.005, 0.0), 550.0).unwrap().is_some());
        assert!(
            tir.trace_ray(cx, ray(0.0075, 0.0), 550.0)
                .unwrap()
                .is_none()
        );
    });
}

#[test]
fn clipped_rays_still_validate_spectral_support_and_cancellation() {
    let lens = biconvex(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0));
    with_cx(|gate, cx| {
        assert!(lens.trace_ray(cx, ray(1.0, 0.0), 550.0).unwrap().is_none());
        assert_eq!(
            lens.trace_ray(cx, ray(1.0, 0.0), 1000.0),
            Err(CameraError::SequentialLensOutsideModel)
        );
        gate.request();
        assert_eq!(
            lens.trace_ray(cx, ray(1.0, 0.0), 550.0),
            Err(CameraError::Cancelled)
        );
    });
}

#[test]
fn weighted_camera_uses_projected_cap_jacobian_and_world_pose() {
    let base = biconvex(glass(CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(), 2.0));
    with_cx(|_, cx| {
        for planar_rear in [false, true] {
            let lens = if planar_rear {
                lens(&[
                    surface(0.02, 0.0, 0.002, None),
                    base.surfaces()[0],
                    base.surfaces()[1],
                ])
            } else {
                base.clone()
            };
            let camera = PhysicalCamera::try_look_at(
                Point3::new(1.0, 2.0, 3.0),
                Point3::new(2.0, 2.0, 3.0),
                Vec3::new(0.0, 0.0, 1.0),
                CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
                1.0,
                Aperture::try_circular(0.001).unwrap(),
            )
            .unwrap()
            .with_sequential_lens(lens.clone())
            .unwrap();
            let rear = lens.surfaces()[0];
            close(
                camera.raster_vertical_half_tan(),
                0.012 / rear.vertex_z_m(),
                2e-14,
            );
            assert_eq!(
                camera.generate_ray_from_tangent_offsets(cx, 0.0, 0.0, LensSample::CENTER),
                Err(CameraError::SequentialLensRequiresWeightedRays)
            );
            // u=.75,v=.5 is the positive x-axis at half the projected disk radius.
            let x = 0.5 * rear.clear_radius_m();
            let c = rear.curvature_per_m();
            let root = (1.0 - (c * x).powi(2)).sqrt();
            let z = rear.vertex_z_m() + c * x * x / (1.0 + root);
            let sensor = Point3::new(-0.025 * rear.vertex_z_m(), 0.012 * rear.vertex_z_m(), 0.0);
            let displacement = Point3::new(x, 0.0, z).delta_from(sensor);
            let distance = displacement.dot(displacement).sqrt();
            let direction = displacement.scale(1.0 / distance);
            let incident = Ray {
                origin: sensor,
                dir: direction,
            };
            // cos(sensor)*cos(cap)*dA_cap/(r²*dA_projected), normalized by z_rear².
            let geometry = (rear.vertex_z_m() / distance).powi(2)
                * direction.z
                * (direction.z - direction.x * c * x / root);
            if planar_rear {
                close(geometry, direction.z.powi(4), 2e-14);
            } else {
                assert!((geometry - direction.z.powi(4)).abs() > 1e-4);
            }
            let wavelengths = [450.0, 650.0];
            let packet = camera
                .generate_sequential_spectral_rays_from_tangent_offsets(
                    cx,
                    0.025,
                    -0.012,
                    LensSample::try_new(0.75, 0.5).unwrap(),
                    &wavelengths,
                )
                .unwrap();
            for (output, wavelength) in packet.into_iter().zip(wavelengths) {
                let direct = trace(&lens, cx, incident, wavelength);
                let output = output.unwrap();
                close(output.weight, geometry * direct.transmittance, 3e-13);
                let p = direct.ray.origin;
                let d = direct.ray.dir;
                for (actual, expected) in [
                    (output.ray.origin.x, 1.0 + p.z),
                    (output.ray.origin.y, 2.0 - p.x),
                    (output.ray.origin.z, 3.0 + p.y),
                    (output.ray.dir.x, d.z),
                    (output.ray.dir.y, -d.x),
                    (output.ray.dir.z, d.y),
                ] {
                    close(actual, expected, 3e-13);
                }
            }
        }
    });
}

#[test]
fn native_film_preserves_spectral_paths_camera_weights_aovs_and_replay() {
    use fs_rand::philox::philox4x32_10;
    use fs_render::aov::{
        CinematicAovConfig, CinematicAovLimits, CinematicAovProfile, CinematicAovProvenance,
        validity,
    };
    use fs_render::camera::{AnimatedCamera, CutSide};
    use fs_render::charts::TriMesh;
    use fs_render::instances::{GeometryInstance, RigidTransform, SharedGeometry};
    use fs_render::lighting::EnvironmentMap;
    use fs_render::motion::{
        ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval,
    };
    use fs_render::spectral::{LAMBDA_MAX, LAMBDA_MIN, cie_x, cie_y, cie_z, lift_rgb, y_integral};
    use fs_render::tracer::{
        Camera, DirectStrategy, Film, Material, PACKET, Primitive, Sampler, Scene, Settings, Shape,
        render_cinematic, render_cinematic_range, render_cinematic_with_aovs,
        trace_cinematic_pixel_sample,
    };
    let mut scene = Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        // Angular stripes make accidentally sharing the hero direction observable.
        environment: Some(
            EnvironmentMap::try_from_linear_srgb(
                2048,
                1,
                (0..2048)
                    .map(|i| [if i % 2 == 0 { 1.0 } else { 0.05 }; 3])
                    .collect(),
                0.0,
            )
            .unwrap(),
        ),
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 0.0),
            forward: Vec3::new(0.0, 0.0, 1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.02,
        },
    };
    let physical = PhysicalCamera::try_look_at(
        scene.camera.eye,
        Point3::new(0.0, 0.0, 1.0),
        scene.camera.up,
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(biconvex(glass(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        2.0,
    )))
    .unwrap();
    let camera = AnimatedCamera::try_static(901, 0.0, 1.0, physical.clone()).unwrap();
    let shutter = ShutterInterval::resolve(
        0.5,
        0.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 16 },
        ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
    )
    .unwrap();
    let settings = Settings {
        width: 1,
        height: 1,
        spp: 16,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::BsdfOnly,
        seed: 13,
    };
    with_cx(|_, cx| {
        let film =
            render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter).unwrap();
        let mut native_sum = [0.0; 3];
        let mut direction_changes_source = false;
        let half_tan = physical.raster_vertical_half_tan();
        let range = LAMBDA_MAX - LAMBDA_MIN;
        let observer_normalization = 1.0 / y_integral();
        let radiance = |direction, wavelength| {
            let (spectrum, scale) = scene
                .environment
                .as_ref()
                .unwrap()
                .evaluate(direction)
                .unwrap()
                .emission;
            spectrum.eval(wavelength) * scale
        };
        for sample in 0..settings.spp {
            // Original IID film/hero and lens stream domains, not a new sampler.
            let key = [settings.seed as u32, (settings.seed >> 32) as u32];
            let unit = |x: u32| f64::from(x) / 4_294_967_296.0;
            let draws = philox4x32_10([0, sample, 0xdead_0001, 0], key).map(unit);
            let pupil = [1, 2].map(|d| unit(philox4x32_10([0, sample, 0x6c65_6e73, d], key)[0]));
            let wavelengths = fs_render::hero_wavelength_packet::<PACKET>(
                LAMBDA_MIN + draws[2] * range,
                LAMBDA_MIN,
                LAMBDA_MAX,
            );
            let rays = physical
                .generate_sequential_spectral_rays_from_tangent_offsets(
                    cx,
                    (2.0 * draws[0] - 1.0) * half_tan,
                    (1.0 - 2.0 * draws[1]) * half_tan,
                    LensSample::try_new(pupil[0], pupil[1]).unwrap(),
                    &wavelengths,
                )
                .unwrap();
            let shared_direction = rays.iter().flatten().next().map(|r| r.ray.dir);
            let mut expected = [0.0; 3];
            for (ray, wavelength) in rays.into_iter().zip(wavelengths) {
                let Some(ray) = ray else {
                    continue;
                };
                let source = radiance(ray.ray.dir, wavelength);
                direction_changes_source |=
                    (source - radiance(shared_direction.unwrap(), wavelength)).abs() > 0.1;
                let spectral_weight =
                    ray.weight * source * range / PACKET as f64 * observer_normalization;
                for (value, observer) in expected.iter_mut().zip([
                    cie_x(wavelength),
                    cie_y(wavelength),
                    cie_z(wavelength),
                ]) {
                    *value += spectral_weight * observer;
                }
            }
            let actual = trace_cinematic_pixel_sample(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter,
                0,
                sample,
            )
            .unwrap();
            for channel in 0..3 {
                close(actual.xyz[channel], expected[channel], 3e-12);
                native_sum[channel] += actual.xyz[channel];
            }
            assert_eq!(actual.emission_xyz, actual.xyz);
            assert_eq!(actual.direct_xyz, [0.0; 3]);
            assert_eq!(actual.indirect_xyz, [0.0; 3]);
            assert!(actual.primary.is_none());
        }
        assert!(
            direction_changes_source,
            "fixture must detect a reused spectral ray"
        );
        assert!(film.xyz[0].iter().all(|v| v.is_finite() && *v > 0.0));
        assert_eq!(film.xyz[0].map(f64::to_bits), native_sum.map(f64::to_bits));
        let mut progressive = Film::try_new(1, 1).unwrap();
        for (from, to) in [(0, 7), (7, settings.spp)] {
            render_cinematic_range(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                &mut progressive,
                from,
                to,
                shutter,
            )
            .unwrap();
        }
        assert_eq!(progressive.spp_done, settings.spp);
        assert_eq!(
            film.xyz[0].map(f64::to_bits),
            progressive.xyz[0].map(f64::to_bits)
        );
        // Axial depth remains defined even though this lens has no pinhole
        // inverse projection. Since b0e26775bf the film exports previous motion
        // through the lens's native inverse; for this static camera and identity
        // instance it is exactly zero. Clipped seeds and camera cuts that leave
        // motion unavailable are covered by sequential_motion.rs's
        // native_clipping_missing_witness_local_inverse_and_camera_cut_have_distinct_outcomes.
        let identity = fs_blake3::hash_domain("org.frankensim.test.sequential-lens", b"plane");
        let plane = TriMesh::new(
            vec![
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        );
        scene.primitives.push(Primitive {
            shape: Shape::Instance(
                GeometryInstance::try_new(
                    1,
                    identity,
                    SharedGeometry::mesh(plane),
                    RigidTransform::identity(),
                )
                .unwrap(),
            ),
            material: Material::Lambertian {
                reflectance: lift_rgb([1.0; 3]),
            },
            emission: None,
        });
        let config = CinematicAovConfig::new(
            CinematicAovProfile::FinalDiagnostic,
            CinematicAovProvenance::try_new(0, 0.5, 0.0, 1.0, identity, identity, identity)
                .unwrap(),
            CinematicAovLimits::default(),
        );
        let hit = render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &Settings { spp: 1, ..settings },
            shutter,
            config,
        )
        .unwrap();
        assert_eq!(hit.primary_count(0), Some(1));
        // Albedo and validity are public through the original EXR AOV output.
        let decoded = fs_img::read_exr(&hit.to_exr().unwrap()).unwrap();
        let channel = |name: &str| {
            f64::from(
                decoded
                    .channels
                    .iter()
                    .find(|c| c.name == name)
                    .unwrap()
                    .data[0],
            )
        };
        close(channel("depth.Z"), 1.0, 1e-6);
        close(channel("normal.Z"), -1.0, 1e-6);
        for name in ["albedo.R", "albedo.G", "albedo.B"] {
            close(channel(name), 1.0, 0.02);
        }
        let mask = channel("diagnostic.validity") as u32;
        assert_eq!(
            mask & (validity::PRIMARY | validity::ALBEDO),
            validity::PRIMARY | validity::ALBEDO
        );
        assert_eq!(mask & validity::PREVIOUS_MOTION, validity::PREVIOUS_MOTION);
        let guides = hit.denoise_guides().unwrap();
        assert_eq!(guides.motion_prev_x()[0].to_bits(), 0.0_f32.to_bits());
        assert_eq!(guides.motion_prev_y()[0].to_bits(), 0.0_f32.to_bits());
    });
}

#[path = "sequential_lens/projection.rs"]
mod projection;
