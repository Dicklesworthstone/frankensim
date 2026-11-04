use super::*;
use crate::dielectric::{DielectricIorSample, TabulatedIor};

#[test]
fn g0_radiance_source_retains_its_band_and_legacy_visible_arithmetic() {
    let knots = vec![(1200.0, 2.0), (1500.0, 5.0), (1800.0, 3.0)];
    let source = RadianceSpectrum::try_new_in_band(knots.clone()).unwrap();
    assert_eq!(source.knots(), knots);
    assert_eq!(source.wavelength_bounds_nm(), (1200.0, 1800.0));
    for (wavelength, radiance) in [
        (1200.0, 2.0),
        (1350.0, 3.5),
        (1500.0, 5.0),
        (1650.0, 4.0),
        (1800.0, 3.0),
    ] {
        assert_eq!(source.evaluate(wavelength).unwrap(), radiance);
    }
    for outside in [550.0, 1199.0, 1801.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            source.evaluate(outside),
            Err(CameraRadianceError::Optics(
                DielectricError::InvalidWavelength
            ))
        ));
    }
    assert!(matches!(
        RadianceSpectrum::try_new(knots),
        Err(CameraRadianceError::InvalidSpectrum)
    ));
    for invalid in [
        vec![],
        vec![(1550.0, 1.0)],
        vec![(0.0, 1.0), (1800.0, 1.0)],
        vec![(-1200.0, 1.0), (1800.0, 1.0)],
        vec![(1200.0, 1.0), (f64::INFINITY, 1.0)],
        vec![(1200.0, 1.0), (1200.0, 2.0)],
        vec![(1800.0, 1.0), (1200.0, 2.0)],
        vec![(1200.0, -1.0), (1800.0, 2.0)],
        vec![(1200.0, 1.0), (1800.0, f64::NAN)],
        vec![(1200.0, 1.0); MAX_RADIANCE_KNOTS + 1],
    ] {
        assert!(matches!(
            RadianceSpectrum::try_new_in_band(invalid),
            Err(CameraRadianceError::InvalidSpectrum)
        ));
    }
    let visible = vec![(380.0, 1.5), (550.0, 0.0), (780.0, 8.0)];
    let legacy = RadianceSpectrum::try_new(visible.clone()).unwrap();
    let declared = RadianceSpectrum::try_new_in_band(visible).unwrap();
    assert_eq!(legacy, declared);
    for wavelength in [380.0, 430.0, 550.0, 610.0, 780.0] {
        assert_eq!(
            legacy.evaluate(wavelength).unwrap().to_bits(),
            declared.evaluate(wavelength).unwrap().to_bits()
        );
    }
}

#[test]
fn g0_rgb_environment_does_not_invent_nonvisible_source_radiance() {
    let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
    for scale in [0.0, 1.0] {
        let environment = EnvironmentRadiance::try_new(&map, scale).unwrap();
        let direction = Vec3::new(0.0, 0.0, 1.0);
        assert!(environment.evaluate(direction, 550.0).is_ok());
        for outside in [250.0, 1550.0, f64::NAN] {
            assert!(matches!(
                environment.evaluate(direction, outside),
                Err(CameraRadianceError::Optics(
                    DielectricError::InvalidWavelength
                ))
            ));
        }
    }
}

fn infrared_glass(
    index: impl Fn(f64) -> f64,
    absorption: BeerLambertAbsorption,
) -> DielectricGlass {
    let samples = core::array::from_fn(|i| {
        let wavelength = 1000.0 + 125.0 * i as f64;
        DielectricIorSample::try_new_in_band(wavelength, index(wavelength)).unwrap()
    });
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(samples, 125.0).unwrap(),
        absorption,
        GlassProvenance::Custom,
    )
}

#[test]
fn g1_infrared_slab_camera_keeps_snell_geometry_fresnel_and_beer_flux() {
    with_cx(|cx, _| {
        let host = infrared_glass(|_| 1.0, BeerLambertAbsorption::try_constant(0.1).unwrap());
        let slab = infrared_glass(
            |wavelength| 1.5 - 0.0001 * (wavelength - 1500.0),
            BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
                (1000.0, 0.1),
                (1500.0, 0.4),
                (2000.0, 0.7),
            ])
            .unwrap(),
        );
        let media = [
            SceneMedium::Homogeneous(&host),
            SceneMedium::Homogeneous(&slab),
        ];
        let spectrum =
            RadianceSpectrum::try_new_in_band(vec![(1200.0, 2.0), (1500.0, 5.0), (1800.0, 3.0)])
                .unwrap();
        for wavelength in [1200.0, 1500.0, 1800.0] {
            // Both interfaces use +z normals. The second reverses medium IDs.
            let n = 1.5 - 0.0001 * (wavelength - 1500.0);
            let sine = 0.6_f64;
            let cosine = 0.8_f64;
            let sine_t = sine / n;
            let cosine_t = (1.0 - sine_t * sine_t).sqrt();
            let slab_distance = 0.4 / cosine_t;
            let source_x = 0.8 * sine / cosine + slab_distance * sine_t;
            let geometry = InstanceScene::try_new(vec![
                plane(1, 0.4, 0.0, 10.0),
                plane(2, 0.8, 0.0, 10.0),
                plane(9, 1.2, source_x, 0.01),
            ])
            .unwrap();
            let native = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::Dielectric {
                            inside: 0,
                            outside: 1,
                            surface: DielectricSurface::SMOOTH,
                        },
                    },
                    SurfaceBinding {
                        object_id: 2,
                        interaction: Interaction::Dielectric {
                            inside: 1,
                            outside: 0,
                            surface: DielectricSurface::SMOOTH,
                        },
                    },
                    opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let emitters = [SurfaceEmitter::new(
                9,
                spectrum.clone(),
                EmissionSides::Both,
            )];
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let result = camera
                .sample_ray_with_samples(
                    cx,
                    Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new(sine, 0.0, cosine),
                    },
                    wavelength,
                    0,
                    settings(),
                    |_| [0.0, 0.0, 0.99],
                )
                .unwrap();
            // Independent lossless s/p Fresnel powers. The reciprocal second
            // interface has the same powers. Scalar proposal probabilities
            // are divided once; polarization survives both transmissions.
            let rs = (cosine - n * cosine_t) / (cosine + n * cosine_t);
            let rp = (n * cosine - cosine_t) / (n * cosine + cosine_t);
            let ts = 1.0 - rs * rs;
            let tp = 1.0 - rp * rp;
            let proposal = 0.5 * (ts + tp);
            let alpha = 0.1 + 0.0006 * (wavelength - 1000.0);
            let beer = (-0.1 * 0.8 / cosine - alpha * slab_distance).exp();
            let conditional = beer * 0.5 * (ts * ts + tp * tp) / (proposal * proposal);
            near(result.transport_weight(), conditional, 3e-6);
            near(
                result.radiance_w_m2_sr_nm(),
                spectrum.evaluate(wavelength).unwrap() * conditional,
                1e-5,
            );
            near(result.path().beam().ray().origin.x, source_x, 2e-6);
            near(result.path().beam().ray().dir.x, sine, 2e-12);
            near(result.path().beam().ray().dir.z, cosine, 2e-12);
            assert_eq!(result.path().medium_index(), 0);
            assert_eq!(result.path().surface_events(), 3);
            assert_eq!(
                result.endpoint(),
                RadianceEndpoint::Surface { object_id: 9 }
            );
            // The original emitted spectrum stays narrower than the media.
            assert!(matches!(
                camera.sample_ray_with_samples(
                    cx,
                    Ray {
                        origin: Point3::new(source_x, 0.0, 0.9),
                        dir: Vec3::new(0.0, 0.0, 1.0)
                    },
                    1100.0,
                    0,
                    settings(),
                    |_| [0.0; 3],
                ),
                Err(CameraRadianceError::Optics(
                    DielectricError::InvalidWavelength
                ))
            ));
        }
    });
}
