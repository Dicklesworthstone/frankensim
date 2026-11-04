use super::super::tests::{settings, with_cx};
use super::*;
use crate::camera::{
    Aperture, CameraProjection, DefocusedDiffraction, DispersiveThinLens, ThinLensElement,
};
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricIor, GlassProvenance, SellmeierIor,
};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform};
use crate::polarization::dielectric::scene::SurfaceBinding;
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceEndpoint, RadianceSpectrum, SurfaceEmitter,
};
use fs_blake3::ContentHash;
use fs_geom::Vec3;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn constant(n: f64) -> DielectricIor {
    SellmeierIor::try_new_in_band([n * n - 1.0, 0.0, 0.0], [0.0; 3], (200.0, 1400.0))
        .unwrap()
        .into()
}

fn glass(ior: DielectricIor) -> DielectricGlass {
    DielectricGlass::new(ior, BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
}

fn lens(exterior: Option<DielectricIor>) -> DispersiveThinLens {
    let elements = [ThinLensElement::try_new(constant(2.0), 20.0, -20.0).unwrap()];
    exterior
        .map_or_else(
            || DispersiveThinLens::try_new(550.0, &elements),
            |ior| DispersiveThinLens::try_new_in_medium(550.0, &elements, ior),
        )
        .unwrap()
}

fn camera(lens: DispersiveThinLens, radius: f64, depths: &[f64]) -> PhysicalCamera {
    let camera = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.002).unwrap(),
        2.0,
        Aperture::try_circular(radius).unwrap(),
    )
    .unwrap();
    let camera = if depths.len() == 1 {
        camera.with_defocused_diffraction(
            DefocusedDiffraction::try_new(Some(depths[0]), 1e-5).unwrap(),
        )
    } else {
        camera.with_layered_diffraction(LayeredDiffraction::try_new(depths, 1e-5).unwrap())
    }
    .unwrap();
    camera.with_dispersive_lens(lens).unwrap()
}

fn plane(id: u64, depth: f64, x: [f64; 2], y: [f64; 2]) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        ContentHash([id as u8; 32]),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [x[0], y[0], depth],
                [x[1], y[0], depth],
                [x[1], y[1], depth],
                [x[0], y[1], depth],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn opaque(ids: &[u64]) -> Vec<SurfaceBinding> {
    ids.iter()
        .map(|&object_id| SurfaceBinding {
            object_id,
            interaction: Interaction::Opaque,
        })
        .collect()
}

fn emitter(id: u64, band: [f64; 2], radiance: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        id,
        RadianceSpectrum::try_new_in_band(vec![(band[0], radiance), (band[1], radiance)]).unwrap(),
        EmissionSides::Both,
    )
}

// Independently recover the raster-only tangent, separating the optical kick
// from source geometry and avoiding a second call to the production sampler.
fn optical_kick(
    camera: &PhysicalCamera,
    width: u32,
    height: u32,
    seed: u64,
    pixel: u32,
    sample: u32,
    ray: Ray,
) -> [f64; 2] {
    let block = fs_rand::philox::philox4x32_10(
        [pixel, sample, 0, 0x6372_7331],
        [seed as u32, (seed >> 32) as u32],
    );
    let u = |i: usize| (f64::from(block[i]) + 0.5) / 4_294_967_296.0;
    let h = camera.raster_vertical_half_tan();
    let x = (2.0 * (f64::from(pixel % width) + u(0)) / f64::from(width) - 1.0) * f64::from(width)
        / f64::from(height)
        * h;
    let y = (1.0 - 2.0 * (f64::from(pixel / width) + u(1)) / f64::from(height)) * h;
    let z = ray.dir.dot(camera.forward());
    [
        ray.dir.dot(camera.right()) / z - x,
        ray.dir.dot(camera.up()) / z - y,
    ]
}

#[cfg(feature = "material-state")]
#[test]
fn water_wave_sensor_retains_geometric_depth_optical_lengths_and_vacuum_photons() {
    use crate::dielectric::IapwsWaterIor;
    with_cx(|cx, _| {
        let water = IapwsWaterIor::try_new(293.15, 998.2).unwrap();
        let ior = DielectricIor::Water(water);
        let water_glass = glass(ior);
        let media = [SceneMedium::Homogeneous(&water_glass)];
        let geometry =
            InstanceScene::try_new(vec![plane(1, 2.0, [-50.0, 50.0], [-50.0, 50.0])]).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let emitters = [emitter(1, [800.0, 1000.0], 3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(lens(Some(ior)), 0.0003, &[2.0]);
        assert!(CurvedCameraRaster::try_new(&camera, 2, 2).is_err());
        let raster = CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 2, 2).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(800.0, [0.25, 0.75]), (1000.0, [0.75, 0.25])])
                .unwrap();
        let result = raster
            .sample_pixel_sensor(cx, &sensor, settings(), 71, 1, 7)
            .unwrap();
        let kick = optical_kick(&camera, 2, 2, 71, 1, 7, result.camera_ray());
        let mut power = [0.0; 2];
        let mut electrons = [0.0; 2];
        let mut densities = [0.0; CAMERA_WAVELENGTHS];
        for (lane, ((draw, path), &weight)) in result
            .detector_wavelengths()
            .iter()
            .zip(result.spectral_samples())
            .zip(result.diffraction_weights())
            .enumerate()
        {
            let wavelength = draw.wavelength_nm();
            let n = water.eval(wavelength).unwrap();
            let ng = water.group_index(wavelength).unwrap();
            assert!(n > 1.3 && ng > n);
            near(path.launch_phase_index(), n, 0.0);
            near(path.radiance_w_m2_sr_nm(), 3.0, 1e-14);
            assert_eq!(path.path().surface_events(), 1);
            let length = path.path().traveled_m();
            near(length, 2.0 / result.camera_rays()[lane].dir.z, 2e-7);
            near(path.path().optical_path_length_m(), n * length, 2e-10);
            near(
                path.path().group_optical_path_length_m().unwrap(),
                ng * length,
                2e-10,
            );
            let point = path.path().last_contact().unwrap().witness.surface_point;
            near(point.z, 2.0, 1e-14);
            camera
                .admit_diffraction_depth(Some(point), &[wavelength])
                .unwrap();
            assert!(
                camera
                    .admit_diffraction_depth(
                        Some(Point3::new(
                            point.x,
                            point.y,
                            path.path().optical_path_length_m()
                        )),
                        &[wavelength],
                    )
                    .is_err(),
                "optical length is not paraxial source depth"
            );
            densities[lane] = camera
                .defocused_diffraction_density(kick, wavelength)
                .unwrap();
            let eta0 = 0.25 + 0.5 * (wavelength - 800.0) / 200.0;
            for (channel, eta) in [eta0, 1.0 - eta0].into_iter().enumerate() {
                let contribution = 3.0 * weight * eta * 200.0 / 4.0;
                power[channel] += contribution;
                electrons[channel] +=
                    contribution * (wavelength * 1e-9) / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        for channel in 0..2 {
            near(
                result.qe_weighted_radiance_w_m2_sr()[channel],
                power[channel],
                1e-11 * power[channel],
            );
            near(
                result.photoelectron_rate_per_m2_sr()[channel],
                electrons[channel],
                1e-11 * electrons[channel],
            );
        }
        for lane in 1..CAMERA_WAVELENGTHS {
            near(
                result.diffraction_weights()[lane] / result.diffraction_weights()[0],
                densities[lane] / densities[0],
                2e-9 * (densities[lane] / densities[0]).max(1.0),
            );
        }
        assert_eq!(
            result,
            raster
                .sample_pixel_sensor(cx, &sensor, settings(), 71, 1, 7)
                .unwrap()
        );
    });
}

#[test]
fn explicit_vacuum_limit_and_black_background_preserve_native_wave_weights() {
    with_cx(|cx, _| {
        let vacuum = glass(constant(1.0));
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let geometry =
            InstanceScene::try_new(vec![plane(1, 2.0, [-50.0, 50.0], [-50.0, 50.0])]).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let emitters = [emitter(1, [380.0, 780.0], 2.5)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let implicit = camera(lens(None), 0.0004, &[2.0]);
        let explicit = camera(lens(Some(constant(1.0))), 0.0004, &[2.0]);
        let a = CurvedWaveCameraRaster::try_new(cx, &implicit, &scene, 0, 2, 2).unwrap();
        let b = CurvedWaveCameraRaster::try_new(cx, &explicit, &scene, 0, 2, 2).unwrap();
        for sample in 0..4 {
            assert_eq!(
                a.sample_pixel(cx, settings(), 19, 1, sample).unwrap(),
                b.sample_pixel(cx, settings(), 19, 1, sample).unwrap()
            );
        }
        let geometry = InstanceScene::try_new(vec![plane(1, 2.0, [5.0, 6.0], [5.0, 6.0])]).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        for camera in [implicit, camera(lens(None), 0.0004, &[1.6, 2.0])] {
            let raster = CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 2, 2).unwrap();
            let black = raster.sample_pixel(cx, settings(), 19, 1, 1).unwrap();
            assert_eq!(black.xyz(), [0.0; 3]);
            assert_eq!(black.diffraction_weights(), &[0.0; CAMERA_WAVELENGTHS]);
            assert!(
                black
                    .spectral_samples()
                    .iter()
                    .all(|s| s.endpoint() == RadianceEndpoint::Environment)
            );
        }
    });
}

#[test]
fn homogeneous_immersed_uniform_source_preserves_pupil_normalization() {
    with_cx(|cx, _| {
        let ior = constant(1.25);
        let carrier = glass(ior);
        let media = [SceneMedium::Homogeneous(&carrier)];
        let geometry =
            InstanceScene::try_new(vec![plane(1, 2.0, [-50.0, 50.0], [-50.0, 50.0])]).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let emitters = [emitter(1, [500.0, 600.0], 3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(lens(Some(ior)), 0.0005, &[2.0]);
        let raster = CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (600.0, [1.0])]).unwrap();
        let mut sum = 0.0;
        let mut square = 0.0;
        const COUNT: u32 = 512;
        for sample in 0..COUNT {
            let pixel = raster
                .sample_pixel_sensor(cx, &sensor, settings(), 109, 0, sample)
                .unwrap();
            let normalized = pixel.qe_weighted_radiance_w_m2_sr()[0] / 300.0;
            sum += normalized;
            square += normalized * normalized;
            for path in pixel.spectral_samples() {
                near(path.radiance_w_m2_sr_nm(), 3.0, 1e-13);
            }
        }
        let mean = sum / f64::from(COUNT);
        let variance = (square / f64::from(COUNT) - mean * mean).max(0.0);
        let tolerance = (6.0 * (variance / f64::from(COUNT)).sqrt()).min(0.2);
        near(mean, 1.0, tolerance);
    });
}

#[test]
fn separated_layers_resolve_original_first_hits_and_their_own_depth_kernels() {
    with_cx(|cx, _| {
        let ior: DielectricIor = CauchyIor::try_new(1.2, 0.01, 0.0).unwrap().into();
        let carrier = glass(ior);
        let media = [SceneMedium::Homogeneous(&carrier)];
        let geometry = InstanceScene::try_new(vec![
            plane(1, 1.6, [-10.0, -0.002], [-10.0, 10.0]),
            plane(2, 2.4, [0.002, 10.0], [-10.0, 10.0]),
        ])
        .unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1, 2]), 100_000).unwrap();
        let emitters = [
            emitter(1, [380.0, 780.0], 2.5),
            emitter(2, [380.0, 780.0], 6.0),
        ];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(lens(Some(ior)), 0.0006, &[1.6, 2.4]);
        let raster = CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 2, 1).unwrap();
        let mut seen = [false; 2];
        for pixel in 0..2 {
            let result = raster.sample_pixel(cx, settings(), 37, pixel, 3).unwrap();
            let kick = optical_kick(&camera, 2, 1, 37, pixel, 3, result.camera_ray());
            let mut densities = [0.0; CAMERA_WAVELENGTHS];
            for (lane, (&wavelength, native)) in result
                .wavelengths_nm()
                .iter()
                .zip(result.spectral_samples())
                .enumerate()
            {
                let contact = native.path().last_contact().unwrap();
                let index = (contact.object_id - 1) as usize;
                seen[index] = true;
                let depth = [1.6, 2.4][index];
                near(contact.witness.surface_point.z, depth, 1e-14);
                near(native.radiance_w_m2_sr_nm(), [2.5, 6.0][index], 1e-13);
                densities[lane] = camera
                    .layered_diffraction_density(kick, depth, wavelength)
                    .unwrap();
            }
            for lane in 1..CAMERA_WAVELENGTHS {
                near(
                    result.diffraction_weights()[lane] / result.diffraction_weights()[0],
                    densities[lane] / densities[0],
                    5e-9 * (densities[lane] / densities[0]).max(1.0),
                );
            }
        }
        assert_eq!(seen, [true, true]);
    });
}

#[test]
fn layered_token_resolves_each_native_source_and_black_lane_independently() {
    with_cx(|cx, _| {
        let camera = camera(lens(Some(constant(1.25))), 0.0006, &[1.6, 2.4]);
        let wavelengths = [400.0, 500.0, 600.0, 700.0];
        let draw = LensSample::try_new(0.31, 0.43).unwrap();
        let token = camera
            .sample_layered_diffraction_ray(cx, 0.0, 0.0, draw, &wavelengths)
            .unwrap();
        let ray = token.ray();
        let points = [
            Some(Point3::new(0.0, 0.0, 1.6)),
            Some(Point3::new(0.0, 0.0, 2.4)),
            None,
            Some(Point3::new(0.0, 0.0, 2.4)),
        ];
        let resolved = token.resolve_spectral_sources(cx, points).unwrap();
        let kick = [ray.dir.x / ray.dir.z, ray.dir.y / ray.dir.z];
        let reference = camera
            .layered_diffraction_density(kick, 1.6, 400.0)
            .unwrap();
        assert!(resolved.weights[0] > 0.0);
        assert_eq!(resolved.weights[2], 0.0);
        for lane in [1, 3] {
            let ratio = camera
                .layered_diffraction_density(kick, 2.4, wavelengths[lane])
                .unwrap()
                / reference;
            near(
                resolved.weights[lane] / resolved.weights[0],
                ratio,
                2e-10 * ratio.max(1.0),
            );
        }
        let mut invalid = points;
        invalid[3] = Some(Point3::new(0.0, 0.0, 1.9));
        assert!(
            camera
                .sample_layered_diffraction_ray(cx, 0.0, 0.0, draw, &wavelengths)
                .unwrap()
                .resolve_spectral_sources(cx, invalid)
                .is_err()
        );
    });
}

#[test]
fn full_pupil_admission_refuses_hidden_edge_occlusion_and_nearby_undeclared_depths() {
    with_cx(|cx, _| {
        let ior = constant(1.25);
        let carrier = glass(ior);
        let media = [SceneMedium::Homogeneous(&carrier)];
        let camera = camera(lens(Some(ior)), 0.001, &[1.0, 2.0]);
        // Every central ray to this far rectangle misses the foreground;
        // a ray from x=-radius reaches x=0.0016 on the foreground plane.
        assert!(0.0042 / 2.0 > 0.002);
        assert!((-0.001 + 0.0042) / 2.0 < 0.002);
        let geometry = InstanceScene::try_new(vec![
            plane(1, 1.0, [-0.002, 0.002], [-0.002, 0.002]),
            plane(2, 2.0, [0.0042, 0.005], [-0.002, 0.002]),
        ])
        .unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1, 2]), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        assert!(matches!(
            CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 1, 1),
            Err(CurvedRasterError::WaveScene(
                WaveSceneError::PupilOcclusion {
                    nearer: 0,
                    farther: 1
                }
            ))
        ));
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                1,
                ContentHash([1; 32]),
                SharedGeometry::mesh(TriMesh::new(
                    vec![
                        [-0.01, -0.01, 1.0],
                        [0.01, -0.01, 1.0],
                        [0.0, 0.01, 1.000001],
                    ],
                    vec![[0, 1, 2]],
                )),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        assert!(matches!(
            CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 1, 1),
            Err(CurvedRasterError::WaveScene(WaveSceneError::SourcePlane {
                object_id: 1
            }))
        ));
    });
}

struct ConstantField;
impl GradientIndexField for ConstantField {
    fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        GradientIndexSample::try_new(1.25, Vec3::new(0.0, 0.0, 0.0), 0.0, Some(1.25))
    }
}

#[test]
fn wave_admission_binds_original_medium_and_refuses_absorption_fields_and_interfaces() {
    with_cx(|cx, _| {
        let ior = constant(1.25);
        let clear = glass(ior);
        let absorbing = DielectricGlass::new(
            ior,
            BeerLambertAbsorption::try_constant(0.001).unwrap(),
            GlassProvenance::Custom,
        );
        let wrong = glass(constant(1.3));
        let same_reference = glass(
            CauchyIor::try_new(1.25 - 0.01 / (0.55 * 0.55), 0.01, 0.0)
                .unwrap()
                .into(),
        );
        near(same_reference.ior().eval(550.0).unwrap(), 1.25, 1e-15);
        let field = ConstantField;
        let media = [
            SceneMedium::Homogeneous(&clear),
            SceneMedium::Homogeneous(&absorbing),
            SceneMedium::Homogeneous(&wrong),
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([3; 32]),
            },
            SceneMedium::Homogeneous(&same_reference),
        ];
        let geometry =
            InstanceScene::try_new(vec![plane(1, 2.0, [-1.0, 1.0], [-1.0, 1.0])]).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let camera = camera(lens(Some(ior)), 0.0003, &[2.0]);
        CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 1, 1).unwrap();
        for (medium, expected) in [
            (1, WaveSceneError::AbsorbingMedium),
            (2, WaveSceneError::ExteriorMismatch),
            (3, WaveSceneError::InvalidMedium),
            (4, WaveSceneError::ExteriorMismatch),
            (5, WaveSceneError::InvalidMedium),
        ] {
            assert!(
                matches!(CurvedWaveCameraRaster::try_new(cx, &camera, &scene, medium, 1, 1),
                Err(CurvedRasterError::WaveScene(error)) if error == expected)
            );
        }
        let bindings = vec![SurfaceBinding {
            object_id: 1,
            interaction: Interaction::Dielectric {
                outside: 0,
                inside: 2,
                surface: crate::dielectric::DielectricSurface::SMOOTH,
            },
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings, 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        assert!(matches!(
            CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 1, 1),
            Err(CurvedRasterError::WaveScene(
                WaveSceneError::UnsupportedSurface { object_id: 1 }
            ))
        ));
    });
}

#[test]
fn wave_sensor_preserves_source_refusals_native_limits_and_cancellation() {
    with_cx(|cx, gate| {
        let ior = constant(1.25);
        let carrier = glass(ior);
        let media = [SceneMedium::Homogeneous(&carrier)];
        let geometry =
            InstanceScene::try_new(vec![plane(1, 2.0, [-50.0, 50.0], [-50.0, 50.0])]).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, opaque(&[1]), 100_000).unwrap();
        let emitters = [emitter(1, [500.0, 600.0], 3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = camera(lens(Some(ior)), 0.0003, &[2.0]);
        let raster = CurvedWaveCameraRaster::try_new(cx, &camera, &scene, 0, 1, 1).unwrap();
        let valid = QuantumEfficiencySensor::try_new(vec![(500.0, [1.0]), (600.0, [1.0])]).unwrap();
        let outside_source =
            QuantumEfficiencySensor::try_new(vec![(700.0, [1.0]), (800.0, [1.0])]).unwrap();
        assert!(matches!(
            raster.sample_pixel_sensor(cx, &outside_source, settings(), 3, 0, 0),
            Err(CurvedRasterError::Transport(_))
        ));
        let outside_material =
            QuantumEfficiencySensor::try_new(vec![(1500.0, [1.0]), (1600.0, [1.0])]).unwrap();
        assert!(matches!(
            raster.sample_pixel_sensor(cx, &outside_material, settings(), 3, 0, 0),
            Err(CurvedRasterError::Camera(_))
        ));
        let mut limited = settings();
        limited.max_distance_m = 0.5;
        assert!(matches!(
            raster.sample_pixel_sensor(cx, &valid, limited, 3, 0, 0),
            Err(CurvedRasterError::Transport(
                CameraRadianceError::Unresolved(PathTermination::DistanceLimit)
            ))
        ));
        assert!(matches!(
            raster.sample_pixel_sensor(cx, &valid, settings(), 3, 1, 0),
            Err(CurvedRasterError::InvalidRaster)
        ));
        gate.request();
        assert!(
            raster
                .sample_pixel_sensor(cx, &valid, settings(), 3, 0, 0)
                .is_err()
        );
    });
}
