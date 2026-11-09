use super::*;

mod dielectric;
use crate::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use core::f64::consts::PI;

fn optics() -> ConductorOptics {
    let source =
        ConductorSource::try_new(ContentHash([37; 32]), ConductorDataStatus::Representative)
            .unwrap();
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        source,
    )
    .unwrap()
}

#[test]
fn g1_smith_point_camera_keeps_absolute_irradiance_beer_and_analytic_first_order() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, 0.0, 20.0), plane(9, 2.0, 0.0, 20.0)])
                .unwrap();
        let n = 1.33;
        let absorption = 0.12;
        let material = glass(n, absorption);
        let media = [SceneMedium::Homogeneous(&material)];
        let conductors = [optics()];
        let alpha = 1.0;
        let scene = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &conductors,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::RoughConductor {
                        medium: 0,
                        optics: 0,
                        surface: ConductorSurface::try_rough(alpha).unwrap(),
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(Point3::new(0.0, 0.0, 1.0));
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, -1.0))];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let smith = SmithPointLightCameraScene::new(&lit, SmithWalkSettings::try_new(2).unwrap());
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.5),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let one = smith
            .sample_ray_with_samples(
                cx,
                ray,
                550.0,
                0,
                settings(),
                |_| [0.0; 3],
                |_, _| [0.5, 0.0, 0.0],
                |_, _, _| [0.5, 0.0, 0.0],
            )
            .unwrap();
        let r = ((0.2 - n) * (0.2 - n) + 9.0) / ((0.2 + n) * (0.2 + n) + 9.0);
        let expected =
            intensity(&lamp) * r / (4.0 * PI * alpha * alpha) * (-1.5 * absorption).exp();
        near(one.fixture_radiance_w_m2_sr_nm(), expected);
        near(one.radiance_sample().radiance_w_m2_sr_nm(), expected);
        assert_eq!(one.connections(), 1);
        let ordinary = lit
            .sample_ray_with_samples(cx, ray, 550.0, 0, settings(), |_| [0.0; 3])
            .unwrap();
        near(
            one.fixture_radiance_w_m2_sr_nm(),
            ordinary.fixture_radiance_w_m2_sr_nm(),
        );
        const TWO: [[f64; 3]; 3] = [[0.5, 0.8, 0.125], [0.5, 0.1, 0.375], [0.5, 0.2, 0.3]];
        let mut evaluation_ids = Vec::new();
        let multiple = smith
            .sample_ray_with_samples(
                cx,
                ray,
                550.0,
                0,
                settings(),
                |_| [0.0; 3],
                |_, _| [0.5, 0.0, 0.0],
                |event, binding, facet| {
                    evaluation_ids.push((event, binding, facet));
                    TWO[facet as usize]
                },
            )
            .unwrap();
        assert!(multiple.fixture_radiance_w_m2_sr_nm() > expected);
        assert_eq!(
            multiple.radiance_sample().path(),
            one.radiance_sample().path()
        );
        assert_eq!(evaluation_ids.len(), 3);
        assert_eq!(
            evaluation_ids.iter().map(|id| id.2).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert!(
            evaluation_ids
                .iter()
                .all(|id| id.0 == evaluation_ids[0].0 && id.1 == 0)
        );
        let replay = smith
            .sample_ray_with_samples(
                cx,
                ray,
                550.0,
                0,
                settings(),
                |_| [0.0; 3],
                |_, _| [0.5, 0.0, 0.0],
                |_, _, _| [0.5, 0.0, 0.0],
            )
            .unwrap();
        assert_eq!(one, replay);
    });
}

#[test]
fn g0_smith_point_camera_propagates_evaluation_budget_and_cancellation() {
    with_cx(|cx, gate| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, 0.0, 20.0), plane(9, 2.0, 0.0, 20.0)])
                .unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let conductors = [optics()];
        let scene = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &conductors,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::RoughConductor {
                        medium: 0,
                        optics: 0,
                        surface: ConductorSurface::try_rough(1.0).unwrap(),
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(Point3::new(0.0, 0.0, 1.0));
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, -1.0))];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let smith = SmithPointLightCameraScene::new(&lit, SmithWalkSettings::try_new(1).unwrap());
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.5),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let budget = smith.sample_ray_with_samples(
            cx,
            ray,
            550.0,
            0,
            settings(),
            |_| [0.0; 3],
            |_, _| [0.5, 0.0, 0.0],
            |_, _, i| {
                if i == 0 {
                    [0.5, 0.8, 0.125]
                } else {
                    [0.5, 0.1, 0.375]
                }
            },
        );
        assert!(matches!(budget,Err(PointLightCameraError::Camera(
            CameraRadianceError::Scene(CurvedSceneError::Scene(
                crate::polarization::dielectric::scene::SceneError::Smith(
                    crate::polarization::dielectric::transport::microsurface::SmithBeamError::FacetLimit{..}
                )
            ))
        ))));
        assert!(
            smith
                .sample_ray_with_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    settings(),
                    |_| [0.0; 3],
                    |_, _| [0.5, 0.0, 0.0],
                    |_, _, _| {
                        gate.request();
                        [0.5, 0.0, 0.0]
                    }
                )
                .is_err()
        );
    });
}
