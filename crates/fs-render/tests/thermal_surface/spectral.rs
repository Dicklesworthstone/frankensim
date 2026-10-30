//! Wavelength-dependent Kirchhoff surfaces retain physical light transport.

use super::*;
use fs_exec::RunId;
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::instances::{RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::spectral::{cie_x, cie_y, cie_z, y_integral};
use fs_render::thermal_surface::{SpectralDiffuseThermalSurface, SpectralEmissivity};
use fs_render::tracer::{
    Film, RenderExecutionConfig, SceneExposure, SceneExposureBudget, SceneExposureFilm, SceneFrame,
    render_motion_range, render_range, render_range_with_execution, render_scene_exposure_range,
    render_with_execution,
};

fn spectral_surface(temperature: f64, curve: SpectralEmissivity) -> Material {
    Material::SpectralDiffuseThermal {
        surface: SpectralDiffuseThermalSurface::try_new(temperature, curve, ContentHash([57; 32]))
            .unwrap(),
    }
}

fn selective_curve() -> SpectralEmissivity {
    SpectralEmissivity::table(&[(380.0, 0.2), (580.0, 0.8), (780.0, 0.2)]).unwrap()
}

fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let wavelength_m = wavelength_nm * 1.0e-9;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    2.0 * h * c * c / wavelength_m.powi(5) / (h * c / (wavelength_m * k * temperature)).exp_m1()
        * 1.0e-9
}

fn quad(x: f64, z: f64, half_width: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [x - half_width, -half_width, z],
                [x + half_width, -half_width, z],
                [x + half_width, half_width, z],
                [x - half_width, half_width, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material,
        emission: None,
    }
}

fn illuminated(source: Material) -> Scene {
    let mut s = scene(Material::Lambertian {
        reflectance: lift_rgb([0.4; 3]),
    });
    s.camera.eye.z = 0.5;
    s.primitives.push(quad(0.0, 2.0, 1.0, source));
    s
}

// Independent area quadrature for a square [-1,1]^2 at z=2 and a diffuse
// receiver at the origin: (1/pi) integral cos(theta_i)cos(theta_l)/r^2 dA.
fn geometric_factor() -> f64 {
    let n = 128;
    let dx = 2.0 / f64::from(n);
    let mut integral = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = -1.0 + (f64::from(i) + 0.5) * dx;
            let y = -1.0 + (f64::from(j) + 0.5) * dx;
            let r2 = x * x + y * y + 4.0;
            integral += 4.0 / (r2 * r2);
        }
    }
    integral * dx * dx / core::f64::consts::PI
}

fn assert_mean(film: &Film, expected: [f64; 3], tolerance: f64) {
    for (actual, expected) in film.xyz[0].iter().zip(expected) {
        let actual = actual / f64::from(film.spp_done);
        assert!(
            expected > 0.0 && (actual / expected - 1.0).abs() < tolerance,
            "spectral radiance {actual} vs {expected}"
        );
    }
}

#[test]
fn g1_spectral_kirchhoff_emission_and_reflection_match_an_independent_spectrum() {
    with_cx(|cx| {
        let curve = SpectralEmissivity::table(&[(380.0, 0.1), (780.0, 0.9)]).unwrap();
        let mut s = scene(spectral_surface(2300.0, curve));
        s.environment =
            Some(EnvironmentMap::try_from_linear_srgb(1, 1, vec![[20.0; 3]], 0.0).unwrap());
        let exterior = s
            .environment
            .as_ref()
            .unwrap()
            .evaluate(Vec3::new(0.0, 0.0, 1.0))
            .unwrap()
            .emission;
        let expected = xyz_of_spectrum(|lambda| {
            let epsilon = 0.1 + 0.8 * (lambda - 380.0) / 400.0;
            epsilon * planck(2300.0, lambda)
                + (1.0 - epsilon) * exterior.0.eval(lambda) * exterior.1
        });
        let settings = Settings {
            spp: 8192,
            ..settings(DirectStrategy::BsdfOnly)
        };
        assert_mean(&render(&s, cx, &settings).unwrap(), expected, 0.015);
    });
}

#[test]
fn g3_constant_spectral_emissivity_reproduces_gray_emission_and_reflection() {
    with_cx(|cx| {
        let settings = Settings {
            spp: 256,
            sampler: Sampler::OwenSobolFullPath,
            ..settings(DirectStrategy::BsdfOnly)
        };
        for epsilon in [0.0, 0.35, 1.0] {
            let gray = Material::GrayDiffuseThermal {
                surface: surface(2400.0, epsilon),
            };
            let spectral = spectral_surface(2400.0, SpectralEmissivity::constant(epsilon).unwrap());
            for indirect in [false, true] {
                let make = |material| {
                    let mut s = if indirect {
                        illuminated(material)
                    } else {
                        scene(material)
                    };
                    s.environment = Some(
                        EnvironmentMap::try_from_linear_srgb(1, 1, vec![[2.0, 1.0, 0.5]], 0.0)
                            .unwrap(),
                    );
                    s
                };
                let reference = render(&make(gray), cx, &settings).unwrap();
                let actual = render(&make(spectral), cx, &settings).unwrap();
                for (actual, expected) in actual.xyz[0].iter().zip(reference.xyz[0]) {
                    assert!(
                        (actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0),
                        "epsilon {epsilon}, indirect {indirect}: {actual} vs {expected}"
                    );
                }
            }
        }
    });
}

#[test]
fn g1_spectral_mesh_lighting_preserves_energy_across_direct_estimators() {
    with_cx(|cx| {
        let s = illuminated(spectral_surface(2600.0, selective_curve()));
        let geometry = geometric_factor();
        let reflectance = lift_rgb([0.4; 3]);
        let expected = xyz_of_spectrum(|lambda| {
            let epsilon = 0.2 + 0.6 * (1.0 - (lambda - 580.0).abs() / 200.0);
            reflectance.eval(lambda) * epsilon * planck(2600.0, lambda) * geometry
        });
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let settings = Settings {
                spp: 16_384,
                sampler: Sampler::OwenSobolFullPath,
                ..settings(strategy)
            };
            assert_mean(&render(&s, cx, &settings).unwrap(), expected, 0.025);
        }
    });
}

fn moving_scene(source: Material) -> Scene {
    let mut s = illuminated(source);
    let Shape::Mesh(mesh) = quad(0.0, 0.0, 1.0, source).shape else {
        unreachable!()
    };
    let trajectory = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(
            0.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-0.6, 0.0, 2.0]).unwrap(),
            [1.2, 0.0, 0.4],
        )
        .unwrap(),
        TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.6, 0.0, 2.4]).unwrap(),
            [1.2, 0.0, 0.4],
        )
        .unwrap(),
    ])
    .unwrap();
    s.primitives[1].shape = Shape::AnimatedInstance(
        AnimatedGeometryInstance::try_new(
            570,
            ContentHash([57; 32]),
            SharedGeometry::mesh(mesh),
            trajectory,
        )
        .unwrap(),
    );
    s
}

#[test]
fn g1_moving_spectral_light_and_accepted_material_changes_use_one_physical_time() {
    with_cx(|cx| {
        let first = moving_scene(spectral_surface(2400.0, selective_curve()));
        let second = moving_scene(spectral_surface(
            2800.0,
            SpectralEmissivity::table(&[(380.0, 0.9), (780.0, 0.1)]).unwrap(),
        ));
        let states = [
            SceneFrame {
                start_s: 0.0,
                end_s: 0.25,
                scene: &first,
            },
            SceneFrame {
                start_s: 0.25,
                end_s: 1.0,
                scene: &second,
            },
        ];
        let shutter = ShutterInterval::resolve(
            0.0,
            1.0,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::StratifiedCounterV1 { strata: 16 },
            ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        let exposure = SceneExposure::try_new(
            cx,
            &states,
            shutter,
            SceneExposureBudget {
                max_frames: 2,
                max_total_primitives: 4,
                max_total_light_triangles: 4,
            },
        )
        .unwrap();
        let settings = Settings {
            spp: 64,
            max_depth: 1,
            sampler: Sampler::OwenSobolFullPath,
            ..settings(DirectStrategy::NeeOnly)
        };
        let mut actual = SceneExposureFilm::try_new(&exposure, &settings).unwrap();
        let mut frozen = [Film::new(1, 1), Film::new(1, 1)];
        let mut expected = [0.0; 3];
        let mut counts = [0; 2];
        for sample in 0..settings.spp {
            let before = [frozen[0].xyz[0], frozen[1].xyz[0]];
            for (index, scene) in [&first, &second].into_iter().enumerate() {
                render_motion_range(
                    scene,
                    cx,
                    &settings,
                    &mut frozen[index],
                    sample,
                    sample + 1,
                    shutter,
                )
                .unwrap();
            }
            let time =
                shutter.time_at(shutter.sample_for_stream(settings.seed, 0, u64::from(sample)));
            let state = usize::from(time >= 0.25);
            counts[state] += 1;
            for (channel, sum) in expected.iter_mut().enumerate() {
                *sum += frozen[state].xyz[0][channel] - before[state][channel];
            }
            render_scene_exposure_range(&exposure, cx, &settings, &mut actual, sample, sample + 1)
                .unwrap();
            for (actual, expected) in actual.film().xyz[0].iter().zip(expected) {
                assert!((actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0));
            }
        }
        assert_eq!(counts, [16, 48]);
        assert!(expected[1] > 0.0);
        assert_ne!(frozen[0].xyz, frozen[1].xyz);
    });
}

#[test]
fn g1_a_sub_nanometre_emission_feature_keeps_a_direct_light_proposal() {
    with_cx(|cx| {
        // All legacy five-nanometre midpoint samples miss this entire feature.
        let curve = SpectralEmissivity::table(&[
            (380.0, 0.0),
            (552.6, 0.0),
            (553.0, 1.0),
            (553.4, 0.0),
            (780.0, 0.0),
        ])
        .unwrap();
        let s = illuminated(spectral_surface(3000.0, curve));
        let geometry = geometric_factor();
        let reflectance = lift_rgb([0.4; 3]);
        let mut expected = [0.0; 3];
        let dl = 0.8 / 256.0;
        let normalization = 1.0 / y_integral();
        for index in 0..256 {
            let lambda = 552.6 + (f64::from(index) + 0.5) * dl;
            let epsilon = 1.0 - (lambda - 553.0).abs() / 0.4;
            let value = reflectance.eval(lambda)
                * epsilon
                * planck(3000.0, lambda)
                * geometry
                * dl
                * normalization;
            for (sum, observer) in
                expected
                    .iter_mut()
                    .zip([cie_x(lambda), cie_y(lambda), cie_z(lambda)])
            {
                *sum += value * observer;
            }
        }
        let settings = Settings {
            spp: 16_384,
            max_depth: 1,
            sampler: Sampler::OwenSobolFullPath,
            ..settings(DirectStrategy::NeeOnly)
        };
        let image = render(&s, cx, &settings).unwrap();
        assert_mean(&image, expected, 0.2);
        assert_eq!(
            render(
                &s,
                cx,
                &Settings {
                    strategy: DirectStrategy::BsdfOnly,
                    ..settings
                }
            )
            .unwrap()
            .xyz,
            vec![[0.0; 3]],
            "the final vertex must receive the narrow source through NEE"
        );
    });
}

#[test]
fn g5_spectral_light_and_receiver_replay_across_tiles_and_progressive_ranges() {
    with_cx(|cx| {
        let mut s = illuminated(spectral_surface(2600.0, selective_curve()));
        s.primitives[0].material = spectral_surface(1.0, selective_curve());
        s.primitives.push(quad(
            3.0,
            2.5,
            0.6,
            Material::GrayDiffuseThermal {
                surface: surface(2900.0, 0.4),
            },
        ));
        s.camera.half_tan = 0.2;
        let execution = |width, height, workers| {
            RenderExecutionConfig::try_new(width, height, workers, 8 << 20, RunId(570)).unwrap()
        };
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let settings = Settings {
                width: 3,
                height: 2,
                spp: 32,
                sampler,
                ..settings(DirectStrategy::PowerMis)
            };
            let serial = render(&s, cx, &settings).unwrap();
            assert!(serial.xyz.iter().all(|pixel| pixel[1] > 0.0));
            for config in [execution(1, 1, 1), execution(2, 1, 3)] {
                assert_eq!(
                    render_with_execution(&s, cx, &settings, &config)
                        .unwrap()
                        .film,
                    serial
                );
            }
            let mut resumed = Film::new(settings.width, settings.height);
            render_range(&s, cx, &settings, &mut resumed, 0, 11).unwrap();
            render_range_with_execution(
                &s,
                cx,
                &settings,
                &mut resumed,
                11,
                settings.spp,
                &execution(2, 1, 3),
            )
            .unwrap();
            assert_eq!(resumed, serial);
        }
    });
}
