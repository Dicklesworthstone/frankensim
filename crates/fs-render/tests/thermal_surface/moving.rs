//! Moving thermal surfaces illuminate the same accepted-time scene they inhabit.

use super::*;
use fs_exec::RunId;
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::tracer::{
    Film, RenderExecutionConfig, render_motion, render_motion_range,
    render_motion_range_with_execution, render_motion_with_execution,
};

const START_S: f64 = 10.0;
const HALF_WIDTH: f64 = 0.8;
const TEMPERATURE_K: f64 = 2600.0;

fn shutter(offset_s: f64, duration_s: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        START_S + offset_s,
        duration_s,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 1024 },
        ShotTimeBounds::try_new(START_S, START_S + 1.0).unwrap(),
    )
    .unwrap()
}

fn square(half_width: f64) -> SharedGeometry {
    SharedGeometry::mesh(TriMesh::new(
        vec![
            [-half_width, -half_width, 0.0],
            [half_width, -half_width, 0.0],
            [half_width, half_width, 0.0],
            [-half_width, half_width, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    ))
}

fn moving_square(
    object_id: u64,
    half_width: f64,
    material: Material,
    start: [f64; 3],
    end: [f64; 3],
    end_angle: f64,
) -> Primitive {
    let velocity = core::array::from_fn(|axis| end[axis] - start[axis]);
    let trajectory = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(
            START_S,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], start).unwrap(),
            velocity,
        )
        .unwrap(),
        TransformKeyframe::try_new(
            START_S + 1.0,
            RigidTransform::try_new(
                [0.0, (0.5 * end_angle).sin(), 0.0, (0.5 * end_angle).cos()],
                end,
            )
            .unwrap(),
            velocity,
        )
        .unwrap(),
    ])
    .unwrap();
    Primitive {
        shape: Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(
                object_id,
                ContentHash([object_id as u8; 32]),
                square(half_width),
                trajectory,
            )
            .unwrap(),
        ),
        material,
        emission: None,
    }
}

fn moving_scene() -> Scene {
    let mut s = scene(Material::Lambertian {
        reflectance: lift_rgb([0.4; 3]),
    });
    s.camera.eye.z = 0.5;
    s.primitives.push(moving_square(
        501,
        HALF_WIDTH,
        Material::GrayDiffuseThermal {
            surface: surface(TEMPERATURE_K, 1.0),
        },
        [-0.75, 0.0, 2.0],
        [0.75, 0.0, 2.4],
        0.65,
    ));
    s
}

// Independent midpoint area integration. The world point is Ry(angle)*(u,v,0)
// plus the translation, and the source normal is (sin(angle),0,cos(angle)).
// Lambertian receiver/source cosine factors give z*abs(n·p)/|p|^4 dA.
fn square_factor(half_width: f64, translation: [f64; 3], angle: f64, sigma: f64) -> f64 {
    let n = 64;
    let step = 2.0 * half_width / f64::from(n);
    let (sin, cos) = angle.sin_cos();
    let normal_dot_position = sin * translation[0] + cos * translation[2];
    let mut result = 0.0;
    for iu in 0..n {
        let u = -half_width + (f64::from(iu) + 0.5) * step;
        for iv in 0..n {
            let v = -half_width + (f64::from(iv) + 0.5) * step;
            let x = translation[0] + cos * u;
            let y = translation[1] + v;
            let z = translation[2] - sin * u;
            let r2 = x * x + y * y + z * z;
            // The optional homogeneous grid fills both the 0.5m camera leg
            // and the complete receiver-to-source segment.
            let transmission = (-sigma * (0.5 + r2.sqrt())).exp();
            result += z * normal_dot_position.abs() / (r2 * r2) * transmission;
        }
    }
    result * step * step / core::f64::consts::PI
}

fn exposure_factor(sigma_at: impl Fn(f64) -> f64) -> f64 {
    let n = 64;
    let mut result = 0.0;
    for it in 0..n {
        let t = (f64::from(it) + 0.5) / f64::from(n);
        result += square_factor(
            HALF_WIDTH,
            [-0.75 + 1.5 * t, 0.0, 2.0 + 0.4 * t],
            0.65 * t,
            sigma_at(t),
        );
    }
    result / f64::from(n)
}

fn expected(factor: f64) -> [f64; 3] {
    let receiver = lift_rgb([0.4; 3]);
    let body = surface(TEMPERATURE_K, 1.0);
    xyz_of_spectrum(|lambda| {
        receiver.eval(lambda) * body.radiance_w_m2_sr_nm(lambda).unwrap() * factor
    })
}

fn assert_mean(image: &Film, expected: [f64; 3], tolerance: f64) {
    for (actual, expected) in image.xyz[0].iter().zip(expected) {
        let mean = actual / f64::from(image.spp_done);
        assert!(
            (mean / expected - 1.0).abs() < tolerance,
            "moving thermal radiance: {mean} vs {expected}"
        );
    }
}

#[test]
fn g1_translating_rotating_thermal_mesh_lights_the_final_diffuse_vertex() {
    with_cx(|cx| {
        let s = moving_scene();
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.sampler = Sampler::OwenSobolFullPath;
        settings.spp = 16_384;
        settings.max_depth = 1;
        let exposure = shutter(0.0, 1.0);
        let reference = render_motion(&s, cx, &settings, exposure).unwrap();
        assert_mean(&reference, expected(exposure_factor(|_| 0.0)), 0.025);
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            assert_eq!(
                render_motion(&s, cx, &settings, exposure).unwrap(),
                reference,
                "a final-vertex source sample has no continuation competitor"
            );
        }
        settings.strategy = DirectStrategy::BsdfOnly;
        assert_eq!(
            render_motion(&s, cx, &settings, exposure).unwrap().xyz,
            vec![[0.0; 3]]
        );
    });
}

#[test]
fn g1_moving_and_static_thermal_techniques_preserve_energy_across_strategies() {
    with_cx(|cx| {
        let mut s = moving_scene();
        let stationary = surface(2800.0, 1.0);
        // This source occupies a disjoint angular region from the moving one.
        // Both have emissivity one, so this two-vertex scene has only direct
        // emission, without further source-to-source reflected illumination.
        s.primitives.push(Primitive {
            shape: Shape::Instance(
                GeometryInstance::try_new(
                    502,
                    ContentHash([52; 32]),
                    square(1.0),
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [3.0, 0.0, 2.0]).unwrap(),
                )
                .unwrap(),
            ),
            material: Material::GrayDiffuseThermal {
                surface: stationary,
            },
            emission: None,
        });
        let moving_factor = exposure_factor(|_| 0.0);
        let static_factor = square_factor(1.0, [3.0, 0.0, 2.0], 0.0, 0.0);
        let reflectance = lift_rgb([0.4; 3]);
        let moving = surface(TEMPERATURE_K, 1.0);
        let reference = xyz_of_spectrum(|lambda| {
            reflectance.eval(lambda)
                * (moving.radiance_w_m2_sr_nm(lambda).unwrap() * moving_factor
                    + stationary.radiance_w_m2_sr_nm(lambda).unwrap() * static_factor)
        });
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.sampler = Sampler::OwenSobolFullPath;
            settings.spp = 16_384;
            settings.max_depth = 2;
            assert_mean(
                &render_motion(&s, cx, &settings, shutter(0.0, 1.0)).unwrap(),
                reference,
                0.04,
            );
        }
    });
}

#[test]
fn g3_moving_source_and_occluder_share_absolute_shutter_time() {
    with_cx(|cx| {
        let mut s = scene(Material::Lambertian {
            reflectance: lift_rgb([0.4; 3]),
        });
        s.camera.eye.z = 0.5;
        s.primitives.push(moving_square(
            503,
            HALF_WIDTH,
            Material::GrayDiffuseThermal {
                surface: surface(TEMPERATURE_K, 1.0),
            },
            [-1.0, 0.0, 2.0],
            [1.0, 0.0, 2.0],
            0.0,
        ));
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.spp = 1024;
        settings.max_depth = 1;
        let exposure = shutter(0.0, 1.0);
        let clear = render_motion(&s, cx, &settings, exposure).unwrap();
        assert!(clear.xyz[0].iter().all(|value| *value > 0.0));

        // At every time, a half-scale blocker halfway along the source rays
        // covers the complete emitter. Its motion must be evaluated at the
        // source's time; a frozen endpoint leaves a substantial moving gap.
        s.primitives.push(moving_square(
            504,
            HALF_WIDTH * 0.5 + 1.0e-5,
            Material::Lambertian {
                reflectance: lift_rgb([0.0; 3]),
            },
            [-0.5, 0.0, 1.0],
            [0.5, 0.0, 1.0],
            0.0,
        ));
        assert_eq!(
            render_motion(&s, cx, &settings, exposure).unwrap().xyz,
            vec![[0.0; 3]]
        );
        let Shape::AnimatedInstance(blocker) = &s.primitives[2].shape else {
            unreachable!()
        };
        s.primitives[2].shape = Shape::Instance(blocker.instance_at(cx, START_S).unwrap());
        let frozen_blocker = render_motion(&s, cx, &settings, exposure).unwrap();
        assert!(frozen_blocker.xyz[0].iter().all(|value| *value > 0.0));
    });
}

fn execution(tile_width: u32, tile_height: u32, workers: usize) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(tile_width, tile_height, workers, 8 << 20, RunId(504)).unwrap()
}

#[test]
fn g5_moving_thermal_sources_replay_across_tiles_workers_and_progressive_ranges() {
    with_cx(|cx| {
        let mut s = moving_scene();
        s.camera.half_tan = 0.2;
        let exposure = shutter(0.0, 1.0);
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let mut settings = settings(DirectStrategy::PowerMis);
            settings.width = 3;
            settings.height = 2;
            settings.spp = 32;
            settings.max_depth = 2;
            settings.sampler = sampler;
            let serial = render_motion(&s, cx, &settings, exposure).unwrap();
            assert!(serial.xyz.iter().all(|pixel| pixel[1] > 0.0));
            for config in [execution(1, 1, 1), execution(2, 1, 3)] {
                assert_eq!(
                    render_motion_with_execution(&s, cx, &settings, exposure, &config)
                        .unwrap()
                        .film,
                    serial
                );
            }
            let mut resumed = Film::new(settings.width, settings.height);
            render_motion_range(&s, cx, &settings, &mut resumed, 0, 11, exposure).unwrap();
            render_motion_range_with_execution(
                &s,
                cx,
                &settings,
                &mut resumed,
                11,
                settings.spp,
                exposure,
                &execution(2, 1, 3),
            )
            .unwrap();
            assert_eq!(resumed, serial);
        }
    });
}

#[cfg(feature = "volumes")]
#[test]
fn g1_moving_thermal_light_reads_the_matching_physical_extinction_snapshot() {
    use fs_render::tracer::{MediumExposure, render_medium_exposure};
    use fs_render::volumes::VolumeGrid;
    use fs_render::volumes::transport::timeline::{
        MediumFrame, MediumTimeline, MediumTimelineBudget, PhysicalMedium,
    };
    use fs_render::volumes::transport::{
        AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum,
    };

    with_cx(|cx| {
        let density = [1.0];
        let make_grid = |sigma| {
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-6.0, -4.0, 0.0], [12.0, 8.0, 4.0]),
                MassAbsorptionSpectrum::constant(sigma).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    max_steps_per_segment: 8,
                },
            )
            .unwrap()
        };
        let first = make_grid(0.1);
        let second = make_grid(1.0);
        let frames = [
            MediumFrame {
                start_s: START_S,
                end_s: START_S + 0.25,
                medium: PhysicalMedium::Absorbing(&first),
            },
            MediumFrame {
                start_s: START_S + 0.25,
                end_s: START_S + 1.0,
                medium: PhysicalMedium::Absorbing(&second),
            },
        ];
        let history =
            MediumTimeline::try_new(cx, &frames, MediumTimelineBudget { max_frames: 2 }).unwrap();
        let exposure = MediumExposure::try_new(cx, &history, shutter(0.0, 1.0)).unwrap();
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.sampler = Sampler::OwenSobolFullPath;
        settings.spp = 16_384;
        settings.max_depth = 1;
        let image = render_medium_exposure(&moving_scene(), &exposure, cx, &settings).unwrap();
        assert_mean(
            &image,
            expected(exposure_factor(|t| if t < 0.25 { 0.1 } else { 1.0 })),
            0.035,
        );
    });
}

#[cfg(feature = "volumes")]
#[test]
fn g1_moving_thermal_source_illuminates_single_scattering_vertices() {
    use fs_render::tracer::{MediumExposure, render_medium_exposure};
    use fs_render::volumes::VolumeGrid;
    use fs_render::volumes::transport::scattering::{
        MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
    };
    use fs_render::volumes::transport::timeline::{
        MediumFrame, MediumTimeline, MediumTimelineBudget, PhysicalMedium,
    };
    use fs_render::volumes::transport::{
        AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum,
    };

    // Integrate the camera's first collision independently of the production
    // free-flight sampler, followed by isotropic source-to-collision transport.
    // All source rays leave through the slab's z=1 face before its side faces.
    let n = 32;
    let dt = 1.0 / f64::from(n);
    let dz = 0.5 / f64::from(n);
    let du = 2.0 * HALF_WIDTH / f64::from(n);
    let sigma_s = 0.5;
    let sigma_t = 0.7;
    let mut factor = 0.0;
    for it in 0..n {
        let t = (f64::from(it) + 0.5) * dt;
        let (sin, cos) = (0.65 * t).sin_cos();
        let translation_x = -0.75 + 1.5 * t;
        let translation_z = 2.0 + 0.4 * t;
        for iz in 0..n {
            let z = (f64::from(iz) + 0.5) * dz;
            let camera_collision = sigma_s * (-sigma_t * (0.5 - z)).exp();
            let normal_dot_displacement = sin * translation_x + cos * (translation_z - z);
            for iu in 0..n {
                let u = -HALF_WIDTH + (f64::from(iu) + 0.5) * du;
                let x = translation_x + cos * u;
                let light_z = translation_z - sin * u;
                let delta_z = light_z - z;
                for iv in 0..n {
                    let y = -HALF_WIDTH + (f64::from(iv) + 0.5) * du;
                    let r2 = x * x + y * y + delta_z * delta_z;
                    let distance = r2.sqrt();
                    let shadow_length = (1.0 - z) * distance / delta_z;
                    let source_kernel = normal_dot_displacement.abs()
                        / (4.0 * core::f64::consts::PI * r2 * distance);
                    factor += camera_collision * source_kernel * (-sigma_t * shadow_length).exp();
                }
            }
        }
    }
    factor *= dt * dz * du * du;
    let body = surface(TEMPERATURE_K, 1.0);
    let reference = xyz_of_spectrum(|lambda| body.radiance_w_m2_sr_nm(lambda).unwrap() * factor);

    with_cx(|cx| {
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-6.0, -4.0, 0.0], [12.0, 8.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    max_steps_per_segment: 8,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::constant(sigma_s).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let frames = [MediumFrame {
            start_s: START_S,
            end_s: START_S + 1.0,
            medium: PhysicalMedium::Scattering(&medium),
        }];
        let history =
            MediumTimeline::try_new(cx, &frames, MediumTimelineBudget { max_frames: 1 }).unwrap();
        let exposure = MediumExposure::try_new(cx, &history, shutter(0.0, 1.0)).unwrap();
        let mut s = moving_scene();
        // Keep only the hot square: the camera sees darkness until a physical
        // volume collision scatters its radiation toward the camera.
        s.primitives = s.primitives.into_iter().skip(1).collect();
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.sampler = Sampler::OwenSobolFullPath;
        settings.spp = 16_384;
        settings.max_depth = 1;
        let nee = render_medium_exposure(&s, &exposure, cx, &settings).unwrap();
        assert_mean(&nee, reference, 0.05);
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            assert_eq!(
                render_medium_exposure(&s, &exposure, cx, &settings).unwrap(),
                nee,
                "terminal phase source connections have no continuation competitor"
            );
        }
    });
}
