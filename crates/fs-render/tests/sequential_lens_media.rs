//! Authored immersion media, radiance invariance and explicit lens group transit.
#![cfg(feature = "tracer")]

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_math::dd::Dd;
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, LensSample, LensSurface, PhysicalCamera,
    SequentialLens,
};
use fs_render::charts::Ray;
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricIorSample, GlassProvenance,
    TabulatedIor,
};
use fs_render::polarization::Stokes;

const C0: f64 = 299_792_458.0;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 14,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn glass(n: f64, sigma: f64) -> DielectricGlass {
    cauchy(n, 0.0, sigma)
}
fn cauchy(a: f64, b: f64, sigma: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(a, b, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(sigma).unwrap(),
        GlassProvenance::Custom,
    )
}
fn surface(z: f64, c: f64, medium: Option<DielectricGlass>) -> LensSurface {
    LensSurface::try_new(z, c, 0.012, medium).unwrap()
}
fn lens(
    sensor: Option<DielectricGlass>,
    interior: DielectricGlass,
    scene: Option<DielectricGlass>,
) -> SequentialLens {
    SequentialLens::try_new_with_sensor_medium(
        550.0,
        sensor,
        &[
            surface(0.08, 20.0, Some(interior)),
            surface(0.09, -20.0, scene),
        ],
    )
    .unwrap()
}
fn camera(lens: SequentialLens) -> PhysicalCamera {
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
fn normal_t(ni: f64, nt: f64) -> f64 {
    4.0 * ni * nt / (ni + nt).powi(2)
}

#[test]
fn immersed_slab_retains_snell_power_fresnel_and_sensor_bulk_path() {
    with_cx(|cx, _| {
        let (ns, ng, nc) = (1.6, 2.0, 1.33);
        let surfaces = [
            LensSurface::try_new(0.02, 0.0, 0.1, Some(glass(ng, 0.7))).unwrap(),
            LensSurface::try_new(0.05, 0.0, 0.1, Some(glass(nc, 0.9))).unwrap(),
        ];
        assert_eq!(
            SequentialLens::try_new(550.0, &surfaces),
            Err(CameraError::InvalidSequentialLens)
        );
        let lens =
            SequentialLens::try_new_with_sensor_medium(550.0, Some(glass(ns, 0.3)), &surfaces)
                .unwrap();
        let si = 0.6_f64;
        let ci = (1.0 - si * si).sqrt();
        let sg = ns * si / ng;
        let cg = (1.0 - sg * sg).sqrt();
        let sc = ns * si / nc;
        let cc = (1.0 - sc * sc).sqrt();
        let result = lens
            .trace_ray(
                cx,
                Ray {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    dir: Vec3::new(si, 0.0, ci),
                },
                550.0,
            )
            .unwrap()
            .unwrap();
        let transmission = |a: f64, b: f64, ca: f64, cb: f64| {
            [
                4.0 * a * b * ca * cb / (a * ca + b * cb).powi(2),
                4.0 * a * b * ca * cb / (b * ca + a * cb).powi(2),
            ]
        };
        let first = transmission(ns, ng, ci, cg);
        let second = transmission(ng, nc, cg, cc);
        let beer = (-0.3 * 0.02 / ci - 0.7 * 0.03 / cg).exp();
        let power = 0.5 * (first[0] * second[0] + first[1] * second[1]) * beer;
        near(result.transmittance, power, 3e-14);
        near(result.ray.dir.x, sc, 3e-14);
        near(result.ray.dir.z, cc, 3e-14);
        near(result.distance_m, 0.02 / ci + 0.03 / cg, 3e-14);
        near(
            result.optical_path_m,
            ns * 0.02 / ci + ng * 0.03 / cg,
            3e-14,
        );
        assert!((result.transmittance - power * (ns / nc).powi(2)).abs() > 0.1);
    });
}

#[test]
fn immersed_camera_uses_sensor_focal_length_and_one_endpoint_radiance_factor() {
    with_cx(|cx, _| {
        let (ng, nc, z, t) = (2.0, 1.33, 0.08, 0.01);
        for ns in [1.0, 1.2, 1.6] {
            let lens = lens(Some(glass(ns, 0.3)), glass(ng, 0.7), Some(glass(nc, 0.9)));
            assert_eq!(lens.sensor_medium(), Some(&glass(ns, 0.3)));
            assert_eq!(lens.scene_medium(), Some(&glass(nc, 0.9)));
            // Thick-lens power and the reduced-angle imaging equation.
            let p1 = (ng - ns) * 20.0;
            let p2 = (ng - nc) * 20.0;
            let power = p1 + p2 - t * p1 * p2 / ng;
            let b = (z / ns) * (1.0 - t * p1 / ng) + t / ng;
            let d = 1.0 - p1 * z / ns - p2 * b;
            let focus = lens.reference_focus().unwrap();
            near(focus.0, ns / power, 3e-14);
            near(focus.1, 0.09 - nc * b / d, 3e-14);
            let paraxial = lens
                .trace_ray(
                    cx,
                    Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new(1e-7, 0.0, 1.0),
                    },
                    550.0,
                )
                .unwrap()
                .unwrap()
                .ray;
            near(
                paraxial.origin.z - paraxial.origin.x * paraxial.dir.z / paraxial.dir.x,
                focus.1,
                3e-9,
            );
            let camera = camera(lens);
            let beam = camera
                .generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                )
                .unwrap()[0]
                .unwrap();
            let scalar = camera
                .generate_sequential_spectral_rays_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                )
                .unwrap()[0]
                .unwrap();
            let expected = normal_t(ns, ng)
                * normal_t(ng, nc)
                * (-0.3 * z - 0.7 * t).exp()
                * (ns / nc).powi(2);
            near(beam.power().i, expected, 3e-14);
            near(scalar.weight, expected, 3e-14);
            assert_eq!(beam.geometry_weight(), 1.0);
            assert_eq!(beam.scene_phase_index(), nc);
            assert_eq!(scalar.scene_phase_index, nc);
            assert_eq!(beam.group_delay_s(), None);
            if ns > nc {
                assert!(beam.power().i > 1.0);
            }
            let row = Stokes {
                i: 0.5,
                q: 0.2,
                u: 0.1,
                v: 0.3,
            };
            let analyzed = camera
                .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                    camera.right(),
                    row,
                )
                .unwrap()[0]
                .unwrap()
                .power();
            near(analyzed.i, row.i * expected, 3e-14);
            near(analyzed.q, row.q * expected, 3e-14);
            near(analyzed.u, -row.u * expected, 3e-14);
            near(analyzed.v, row.v * expected, 3e-14);
        }
        let caps = [
            surface(0.08, 20.0, Some(glass(1.5, 0.0))),
            surface(0.09, -20.0, None),
        ];
        assert_eq!(
            SequentialLens::try_new(550.0, &caps).unwrap(),
            SequentialLens::try_new_with_sensor_medium(550.0, None, &caps).unwrap()
        );
    });
}

#[test]
fn timed_beams_integrate_original_group_laws_on_actual_segments_in_dd_seconds() {
    with_cx(|cx, _| {
        let camera = camera(lens(
            Some(cauchy(1.1, 0.003, 0.2)),
            cauchy(1.8, 0.01, 0.4),
            Some(glass(1.33, 0.8)),
        ));
        let sample = LensSample::try_new(0.75, 0.5).unwrap();
        let wavelengths = [450.0, 650.0];
        let steady = camera
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
            )
            .unwrap();
        let timed = camera
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
            )
            .unwrap();
        let analyzed = camera
            .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                cx,
                0.025,
                -0.013,
                sample,
                &wavelengths,
                camera.right(),
                Stokes {
                    i: 0.5,
                    q: 0.0,
                    u: 0.0,
                    v: 0.5,
                },
            )
            .unwrap();
        let x = 0.006_f64;
        let sag = 20.0 * x * x / (1.0 + (1.0 - (20.0 * x).powi(2)).sqrt());
        // This pose maps local +x to world -x. The sensor ray is fixed by the
        // projected rear-cap sample; the glass path bends with wavelength.
        let rear = Point3::new(-x, 0.0, 0.08 + sag);
        let sensor = Point3::new(0.025 * 0.08, 0.013 * 0.08, 0.0);
        let sensor_length = rear.delta_from(sensor).norm();
        let mut delays = [0.0; 2];
        for (i, nm) in wavelengths.into_iter().enumerate() {
            let beam = timed[i].unwrap();
            let original = steady[i].unwrap();
            assert_eq!(beam.ray(), original.ray());
            assert_eq!(beam.frame(), original.frame());
            assert_eq!(beam.power(), original.power());
            assert_eq!(beam.geometry_weight(), original.geometry_weight());
            assert_eq!(original.group_delay_s(), None);
            let delay = beam.group_delay_s().unwrap();
            assert_eq!(analyzed[i].unwrap().group_delay_s(), Some(delay));
            near(analyzed[i].unwrap().power().i, 0.5 * beam.power().i, 2e-14);
            let glass_length = beam.ray().origin.delta_from(rear).norm();
            let inv_um2 = 1.0 / (nm * 1e-3).powi(2);
            let gs = 1.1 + 3.0 * 0.003 * inv_um2;
            let gg = 1.8 + 3.0 * 0.01 * inv_um2;
            let expected = (gs * sensor_length + gg * glass_length) / C0;
            near(delay.to_f64(), expected, 3e-24);
            let phase_delay = ((1.1 + 0.003 * inv_um2) * sensor_length
                + (1.8 + 0.01 * inv_um2) * glass_length)
                / C0;
            assert!(delay.to_f64() - phase_delay > 1e-12);
            // The short delay disappears in a binary64 epoch but survives DD.
            let epoch = Dd::from_f64(1e9);
            assert_eq!(1e9_f64 - delay.to_f64(), 1e9);
            near((epoch - (epoch - delay) - delay).to_f64(), 0.0, 3e-23);
            delays[i] = delay.to_f64();
        }
        assert!(delays[0] - delays[1] > 1e-12);
    });
}

#[test]
fn group_derivatives_are_opt_in_and_final_exterior_transit_belongs_to_the_scene() {
    with_cx(|cx, gate| {
        let table = |negative_group: bool| {
            let ior = TabulatedIor::try_new(core::array::from_fn(|i| {
                let nm = 380.0 + 50.0 * i as f64;
                let index = if negative_group {
                    1.1 + 0.004 * (nm - 380.0)
                } else {
                    1.7 + 0.002 * i as f64 + if i == 3 { 0.01 } else { 0.0 }
                };
                DielectricIorSample::try_new(nm, index).unwrap()
            }))
            .unwrap();
            DielectricGlass::new(ior, BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
        };
        let negative = table(true);
        assert!(negative.ior().group_index(535.0).unwrap() < 0.0);
        for (material, wavelength) in [(table(false), 530.0), (negative, 535.0)] {
            let camera = camera(lens(None, material, None));
            assert!(
                camera
                    .generate_sequential_spectral_beams_from_tangent_offsets(
                        cx,
                        0.0,
                        0.0,
                        LensSample::CENTER,
                        &[wavelength],
                    )
                    .unwrap()[0]
                    .is_some()
            );
            assert_eq!(
                camera.generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[wavelength],
                ),
                Err(CameraError::SequentialLensOutsideModel)
            );
        }
        // The same cusp in the FINAL medium has no finite lens segment. The
        // later native scene walk must supply its own group-delay admission.
        let camera = camera(lens(None, glass(2.0, 0.0), Some(table(false))));
        let timed = camera
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[530.0],
            )
            .unwrap()[0]
            .unwrap();
        near(
            timed.group_delay_s().unwrap().to_f64(),
            (0.08 + 2.0 * (0.09 - 0.08)) / C0,
            3e-24,
        );
        gate.request();
        assert_eq!(
            camera.generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[530.0],
            ),
            Err(CameraError::Cancelled)
        );
    });
}
