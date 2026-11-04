use super::*;
use crate::camera::{
    Aperture, CameraError, CameraProjection, DispersiveThinLens, PhysicalCamera, ThinLensElement,
};
use crate::charts::TriMesh;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};

fn lens_camera() -> PhysicalCamera {
    let law = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            let wavelength = 8000.0 + 750.0 * i as f64;
            DielectricIorSample::try_new_in_band(wavelength, 1.5 + (11000.0 - wavelength) * 2e-5)
                .unwrap()
        }),
        750.0,
    )
    .unwrap();
    let lens = DispersiveThinLens::try_new(
        11000.0,
        &[ThinLensElement::try_new(law.into(), 20.0, -20.0).unwrap()],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        2.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_geometric_dispersive_lens(lens)
    .unwrap()
}

fn infrared_vacuum() -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.0).unwrap()
            }),
            750.0,
        )
        .unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

struct GradientFromVacuum;
impl GradientIndexField for GradientFromVacuum {
    fn sample(
        &self,
        p: Point3,
        wavelength: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if !(8000.0..=14000.0).contains(&wavelength) {
            return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
        }
        let color = (11000.0 / wavelength).powi(2);
        let slope = 0.02 + 0.01 * color;
        let n = 1.0 + slope * p.z;
        GradientIndexSample::try_new(
            n,
            Vec3::new(0.0, 0.0, slope),
            0.2,
            Some(n + 0.02 * p.z * color),
        )
    }
}

#[test]
fn g1_native_infrared_lens_and_grin_keep_each_camera_ray_in_energy_and_electron_response() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let field = GradientFromVacuum;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([81; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = lens_camera();
        let raster = CurvedCameraRaster::try_new(&camera, 2, 2).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [0.8]), (14000.0, [0.8])]).unwrap();
        let sample = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 79, 0, 6)
            .unwrap();
        let mut energy = 0.0;
        let mut electrons = 0.0;
        for ((initial, draw), path) in sample
            .camera_rays()
            .iter()
            .zip(sample.detector_wavelengths())
            .zip(sample.spectral_samples())
        {
            assert_eq!(initial.origin, sample.camera_rays()[0].origin);
            assert_eq!(path.launch_phase_index(), 1.0);
            // Plane-stratified n(z,lambda) conserves transverse n*sin(theta).
            // Its closed-form arc length uses THIS wavelength's camera ray.
            let n1 = 1.04 + 0.02 * (11000.0 / draw.wavelength_nm()).powi(2);
            let transverse = initial.dir.x.powi(2) + initial.dir.y.powi(2);
            let length =
                2.0 * (n1 + 1.0) / ((n1 * n1 - transverse).sqrt() + (1.0 - transverse).sqrt());
            let radiance = 3.0 / (n1 * n1) * (-0.2 * length).exp();
            near(path.path().traveled_m(), length, 3e-7);
            near(path.path().beam().ray().dir.x, initial.dir.x / n1, 3e-8);
            near(path.radiance_w_m2_sr_nm(), radiance, 3e-7);
            let power = radiance * 6000.0 * 0.8 / 4.0;
            energy += power;
            electrons += power * draw.wavelength_nm() * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
        }
        assert_ne!(sample.camera_rays()[0].dir, sample.camera_rays()[1].dir);
        near(sample.qe_weighted_radiance_w_m2_sr()[0], energy, 2e-3);
        near(
            sample.photoelectron_rate_per_m2_sr()[0],
            electrons,
            electrons * 2e-7,
        );
        assert_eq!(
            sample,
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 79, 0, 6)
                .unwrap()
        );
    });
}

#[test]
fn g1_chromatic_pupil_rays_change_actual_native_surface_visibility_and_preserve_zero_scattering() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let vacuum = infrared_vacuum();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = lens_camera();
        let raster = CurvedCameraRaster::try_new(&camera, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [1.0]), (14000.0, [1.0])]).unwrap();
        let full = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 0, 4)
            .unwrap();
        let intercepts = full
            .camera_rays()
            .map(|ray| ray.origin.x + 2.0 * ray.dir.x / ray.dir.z);
        let minimum = intercepts.into_iter().fold(f64::INFINITY, f64::min);
        let maximum = intercepts.into_iter().fold(f64::NEG_INFINITY, f64::max);
        assert!(maximum - minimum > 1e-6);
        let edge = 0.5 * minimum + 0.5 * maximum;
        let clipped = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                7,
                ContentHash([82; 32]),
                SharedGeometry::mesh(TriMesh::new(
                    vec![
                        [-50.0, -50.0, 2.0],
                        [edge, -50.0, 2.0],
                        [edge, 50.0, 2.0],
                        [-50.0, 50.0, 2.0],
                    ],
                    vec![[0, 1, 2], [0, 2, 3]],
                )),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &clipped, &media, bindings(), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let sample = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 0, 4)
            .unwrap();
        let mut lit = 0;
        for (&intercept, path) in intercepts.iter().zip(sample.spectral_samples()) {
            let expected = if intercept < edge {
                lit += 1;
                3.0
            } else {
                0.0
            };
            near(path.radiance_w_m2_sr_nm(), expected, 1e-12);
        }
        assert!((1..4).contains(&lit));
        near(
            sample.qe_weighted_radiance_w_m2_sr()[0],
            3.0 * 6000.0 * f64::from(lit) / 4.0,
            1e-9,
        );
        #[cfg(feature = "volumes")]
        {
            let spectrum = MassScatteringSpectrum::source_table(
                &[(8000.0, 0.0), (14000.0, 0.0)],
                ContentHash([83; 32]),
            )
            .unwrap();
            let phase = PhaseFunction::Rayleigh;
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let participating =
                CurvedRayleighCameraScene::try_new(cx, &scene, &scattering).unwrap();
            let replay = raster
                .sample_pixel_sensor_rayleigh(
                    cx,
                    &participating,
                    &sensor,
                    0,
                    CurvedRayleighTraceSettings {
                        surfaces: settings(),
                        max_volume_events: 1,
                    },
                    73,
                    0,
                    4,
                )
                .unwrap();
            assert_eq!(*replay.pixel_sample(), sample);
            assert_eq!(*replay.volume_events(), [0; CAMERA_WAVELENGTHS]);
        }
    });
}

#[test]
fn g0_lens_source_domain_and_actual_exterior_index_are_not_silently_changed() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let glass = infrared_glass(); // n=1.2 cannot use an index-one prescription.
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let camera = lens_camera();
        let raster = CurvedCameraRaster::try_new(&camera, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [1.0]), (14000.0, [1.0])]).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 2, 0, 0)
                .unwrap_err(),
            CurvedRasterError::Camera(CameraError::DispersiveLensOutsideModel)
        );
        let outside =
            QuantumEfficiencySensor::try_new(vec![(7000.0, [1.0]), (7500.0, [1.0])]).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor(cx, &scene, &outside, 0, settings(), 2, 0, 0)
                .unwrap_err(),
            CurvedRasterError::Camera(CameraError::DispersiveLensOutsideModel)
        );
        gate.request();
        assert_eq!(
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 2, 0, 0)
                .unwrap_err(),
            CurvedRasterError::Camera(CameraError::Cancelled)
        );
    });
}

#[test]
fn g1_immersed_dispersive_exterior_reaches_native_rays_beer_and_detector_electrons() {
    with_cx(|cx, _| {
        let exterior = TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                let wavelength = 8000.0 + 750.0 * i as f64;
                DielectricIorSample::try_new_in_band(
                    wavelength,
                    1.2 + (11000.0 - wavelength) * 1e-5,
                )
                .unwrap()
            }),
            750.0,
        )
        .unwrap();
        let glass = DielectricGlass::new(
            exterior,
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let lens_material = crate::dielectric::SellmeierIor::try_new_in_band(
            [3.0, 0.0, 0.0],
            [0.0; 3],
            (8000.0, 14000.0),
        )
        .unwrap();
        let lens = DispersiveThinLens::try_new_in_medium(
            11000.0,
            &[ThinLensElement::try_new(lens_material.into(), 20.0, -20.0).unwrap()],
            exterior.into(),
        )
        .unwrap();
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 2.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_focal_sensor(lens.reference_focal_length_m(), 0.002).unwrap(),
            2.0,
            Aperture::try_circular(0.001).unwrap(),
        )
        .unwrap()
        .with_geometric_dispersive_lens(lens)
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&camera, 2, 2).unwrap();
        let geometry = geometry();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [infrared_emitter(3.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [0.8]), (14000.0, [0.8])]).unwrap();
        let sample = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 1, 6)
            .unwrap();
        assert_ne!(sample.camera_rays()[0].dir, sample.camera_rays()[1].dir);
        let inverse_sensor = (2.0 / 1.2 - 1.0) * 40.0 - 0.5;
        let mut raster_tangent = None;
        let (mut energy, mut electrons) = (0.0, 0.0);
        for ((ray, draw), path) in sample
            .camera_rays()
            .iter()
            .zip(sample.detector_wavelengths())
            .zip(sample.spectral_samples())
        {
            let wavelength = draw.wavelength_nm();
            let ambient = 1.2 + (11000.0 - wavelength) * 1e-5;
            near(path.launch_phase_index(), ambient, 3e-16);
            let power = (2.0 / ambient - 1.0) * 40.0;
            let pupil = ray.origin.delta_from(camera.eye());
            let axial = ray.dir.dot(camera.forward());
            let recovered_raster = [
                ray.dir.dot(camera.right()) / axial
                    - pupil.dot(camera.right()) * (inverse_sensor - power),
                ray.dir.dot(camera.up()) / axial
                    - pupil.dot(camera.up()) * (inverse_sensor - power),
            ];
            if let Some(reference) = raster_tangent {
                for (a, b) in recovered_raster.into_iter().zip(reference) {
                    near(a, b, 2e-15);
                }
            } else {
                raster_tangent = Some(recovered_raster);
            }
            let length = 2.0 / axial;
            near(
                path.path().traveled_m(),
                length,
                settings().event.surface_tolerance_m / axial,
            );
            // Surface isolation retains its declared spatial tolerance;
            // attenuation uses the path length actually transported.
            let radiance = 3.0 * (-0.2 * path.path().traveled_m()).exp();
            // Camera and emitter occupy the SAME dispersive material. The
            // original n_camera²/n_source² factor is one, not an extra n².
            near(path.radiance_w_m2_sr_nm(), radiance, 2e-12);
            let power = radiance * 6000.0 * 0.8 / 4.0;
            energy += power;
            // Photon energy retains VACUUM wavelength although the lens's
            // opening admission uses its shorter medium wavelength.
            electrons += power * wavelength * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
        }
        near(sample.qe_weighted_radiance_w_m2_sr()[0], energy, 2e-8);
        near(
            sample.photoelectron_rate_per_m2_sr()[0],
            electrons,
            electrons * 2e-12,
        );
        assert_eq!(
            sample,
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 1, 6)
                .unwrap()
        );

        let different = infrared_glass();
        let media = [SceneMedium::Homogeneous(&different)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 1, 6)
                .unwrap_err(),
            CurvedRasterError::Camera(CameraError::DispersiveLensOutsideModel),
        );
    });
}
