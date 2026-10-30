//! Independent geometric-optics references for a finite projector pupil.
//! The oracle integrates the declared pupil/gobo ray-power law, without
//! calling the production projector evaluator, sampler or surface BSDF.

#[path = "projector/splitter.rs"]
mod splitter;

#[path = "projector/diffuser.rs"]
mod diffuser;

use super::*;
use fs_exec::RunId;
use fs_render::instances::{GeometryInstance, SharedGeometry};
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{
    EmissionSidedness, GoboTransmission, PhotometricAreaEmitter, PhotometricEmitters,
    ThinLensProjector,
};
use fs_render::photometry::path::{
    IlluminanceSampling, PathIlluminanceEstimate, PathIlluminanceJob,
};
use fs_render::tracer::RenderExecutionConfig;

struct ProjectorCase {
    half_pupil: f64,
    gobo_size: [f64; 2],
    gobo_shape: [usize; 2],
    transmission: Vec<f64>,
    distance: f64,
    focal_length: f64,
    clear_flux: f64,
}

impl ProjectorCase {
    fn new(focus: f64, shape: [usize; 2], transmission: Vec<f64>) -> Self {
        let distance = 0.1;
        Self {
            half_pupil: 0.05,
            gobo_size: [0.04, 0.02],
            gobo_shape: shape,
            transmission,
            distance,
            focal_length: distance * focus / (distance + focus),
            clear_flux: 240.0,
        }
    }

    fn emitter(&self, dimmer: f64, frame: RigidTransform) -> PhotometricAreaEmitter {
        let gobo = GoboTransmission::try_new(
            self.gobo_size[0],
            self.gobo_size[1],
            self.gobo_shape[0],
            self.gobo_shape[1],
            self.transmission.clone(),
        )
        .unwrap();
        let projector =
            ThinLensProjector::try_new(frame, self.distance, self.focal_length, gobo).unwrap();
        PhotometricAreaEmitter::try_new(
            0,
            self.clear_flux,
            dimmer,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()
        .with_projector(projector)
        .unwrap()
    }

    fn scene(&self) -> Scene {
        let mut scene = empty_scene([0.0, 0.0, 3.0], [0.0, 0.0, -1.0]);
        scene
            .primitives
            .push(floor(0.0, self.half_pupil, LiftedSpectrum { c: [0.0; 3] }));
        scene
    }

    fn lux(&self, point: [f64; 2], plane: f64) -> f64 {
        // A ray from pupil a to screen p maps to g = c*a - s*p/D.
        // Integrating uniform pupil power is therefore a rectangle-overlap
        // calculation for each constant-transmission cell. At focus c=0:
        // the complete pupil contributes to exactly the inverted gobo point.
        let c = 1.0 - self.distance / self.focal_length + self.distance / plane;
        let mut coverage = 0.0;
        for row in 0..self.gobo_shape[1] {
            for column in 0..self.gobo_shape[0] {
                let mut fraction = 1.0;
                for (axis, cell) in [column, row].into_iter().enumerate() {
                    let size = self.gobo_size[axis];
                    let count = self.gobo_shape[axis] as f64;
                    let low = -0.5 * size + size * cell as f64 / count;
                    let high = -0.5 * size + size * (cell + 1) as f64 / count;
                    let offset = -self.distance * point[axis] / plane;
                    fraction *= if c.abs() < 1.0e-14 {
                        if offset >= low && offset < high {
                            1.0
                        } else {
                            0.0
                        }
                    } else {
                        let left = (low - offset) / c;
                        let right = (high - offset) / c;
                        (left.max(right).min(self.half_pupil)
                            - left.min(right).max(-self.half_pupil))
                        .max(0.0)
                            / (2.0 * self.half_pupil)
                    };
                }
                coverage += self.transmission[row * self.gobo_shape[0] + column] * fraction;
            }
        }
        self.clear_flux * self.distance.powi(2)
            / (self.gobo_size[0] * self.gobo_size[1] * plane.powi(2))
            * coverage
    }
}

fn probe(
    scene: &Scene,
    sources: &[PhotometricAreaEmitter],
    position: Point3,
    normal: Vec3,
    samples: u32,
    cx: &Cx<'_>,
) -> PathIlluminanceEstimate {
    let probes = [IlluminanceProbe::try_new(position, normal).unwrap()];
    let mut job = PathIlluminanceJob::new_with_emitters(
        scene,
        PhotometricEmitters {
            points: &[],
            surfaces: sources,
        },
        PathMedium::Vacuum,
        &probes,
        TransportOptions::default(),
        IlluminanceSampling {
            samples_per_probe: samples,
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            seed: 173,
        },
        cx,
    )
    .unwrap();
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let result = job.probe_estimate(0).unwrap();
    assert_eq!(result.direct_lux, 0.0);
    assert!(result.total_error_bound_lux.is_none());
    assert!(result.estimated_standard_error_lux.is_none());
    result
}

fn pixel(
    scene: &Scene,
    sources: &[PhotometricAreaEmitter],
    medium: PathMedium<'_, '_>,
    policy: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceEstimate {
    let mut job = PathLuminanceJob::new_with_emitters(
        scene,
        PhotometricEmitters {
            points: &[],
            surfaces: sources,
        },
        medium,
        TransportOptions::default(),
        policy,
        cx,
    )
    .unwrap();
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let result = job.pixel_estimate(0).unwrap();
    assert!(result.total_error_bound_cd_m2.is_none());
    result
}

#[test]
fn g1_projector_focus_has_inverted_image_size_and_absolute_lux() {
    with_cx(|cx, _| {
        for focus in [3.0, 6.0] {
            let case = ProjectorCase::new(focus, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
            let scene = case.scene();
            let sources = [case.emitter(1.0, RigidTransform::identity())];
            // Fixed gobo-cell centres map to positions proportional to q/s.
            // The unequal cells expose either sign error or a swapped axis.
            let image_cell_area =
                case.gobo_size[0] * case.gobo_size[1] * (focus / case.distance).powi(2) / 4.0;
            let mut received_flux = 0.0;
            for gobo in [
                [-0.01, -0.005],
                [0.01, -0.005],
                [-0.01, 0.005],
                [0.01, 0.005],
                [-0.022, 0.0],
            ] {
                let point = gobo.map(|coordinate| -focus / case.distance * coordinate);
                let expected = case.lux(point, focus);
                let result = probe(
                    &scene,
                    &sources,
                    Point3::new(point[0], point[1], focus),
                    Vec3::new(0.0, 0.0, -1.0),
                    8192,
                    cx,
                );
                close(
                    result.illuminance_lux,
                    expected,
                    0.003,
                    "focused inverted gobo lux",
                );
                received_flux += result.illuminance_lux * image_cell_area;
            }
            // The fifth sample is outside the image and contributes zero.
            // The four constant image cells partition all outgoing ray power.
            close(
                received_flux,
                sources[0].flux_lm(),
                0.003,
                "projected screen flux",
            );
        }
    });
}

#[test]
fn g1_projector_defocus_matches_exact_pupil_mask_overlaps() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 1], vec![1.0, 0.0]);
        let scene = case.scene();
        let sources = [case.emitter(1.0, RigidTransform::identity())];
        for plane in [1.5, 6.0] {
            let c = 1.0 - case.distance / case.focal_length + case.distance / plane;
            let blur_half_width = plane / case.distance * c.abs() * case.half_pupil;
            for fraction in [-0.5, 0.0, 0.5] {
                let point = [fraction * blur_half_width, 0.0];
                let expected = case.lux(point, plane);
                let full = case.clear_flux * case.distance.powi(2)
                    / (case.gobo_size[0] * case.gobo_size[1] * plane.powi(2));
                close(
                    expected / full,
                    0.5 + 0.5 * fraction,
                    1.0e-12,
                    "independent overlap fraction",
                );
                let result = probe(
                    &scene,
                    &sources,
                    Point3::new(point[0], point[1], plane),
                    Vec3::new(0.0, 0.0, -1.0),
                    32_768,
                    cx,
                );
                close(
                    result.illuminance_lux,
                    expected,
                    0.035,
                    "defocused finite pupil lux",
                );
            }
        }
    });
}

#[test]
fn g1_projector_passive_gobo_reduces_flux_without_renormalizing_clear_rays() {
    with_cx(|cx, _| {
        let mut case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let mut scene = case.scene();
        let pupil = Point3::new(0.02, -0.01, 0.0);
        let gobo = [-0.008, -0.004];
        let slope = [
            (pupil.x * (1.0 - case.distance / case.focal_length) - gobo[0]) / case.distance,
            (pupil.y * (1.0 - case.distance / case.focal_length) - gobo[1]) / case.distance,
        ];
        let direction = Vec3::new(slope[0], slope[1], 1.0);
        let direction = direction.scale(1.0 / direction.norm());
        scene.camera.eye = pupil.offset(direction.scale(3.0));
        scene.camera.forward = direction.scale(-1.0);
        let up = Vec3::new(0.0, 1.0, -slope[1]);
        scene.camera.up = up.scale(1.0 / up.norm());
        let clear_luminance = case.clear_flux / (4.0 * case.half_pupil.powi(2))
            * case.distance.powi(2)
            / (case.gobo_size[0] * case.gobo_size[1] * direction.z.powi(4));
        for (mask, mean, ray_transmission) in [
            (vec![1.0, 0.25, 0.5, 0.0], 0.4375, 1.0),
            (vec![0.25; 4], 0.25, 0.25),
            (vec![0.0; 4], 0.0, 0.0),
        ] {
            case.transmission = mask;
            let source = case.emitter(0.4, RigidTransform::identity());
            close(
                source.clear_flux_lm(),
                0.4 * case.clear_flux,
                1.0e-15,
                "dimmed clear flux",
            );
            close(
                source.flux_lm(),
                0.4 * case.clear_flux * mean,
                1.0e-15,
                "passive masked flux",
            );
            let value = pixel(&scene, &[source], PathMedium::Vacuum, settings(1, 8192), cx);
            close(
                value.luminance_cd_m2,
                0.4 * clear_luminance * ray_transmission,
                0.002,
                "directional source-hit radiance before display transforms",
            );
        }
    });
}

#[test]
fn g1_projector_nee_and_source_hits_share_native_lambertian_transport() {
    with_cx(|cx, _| {
        let mut case = ProjectorCase::new(1.0, [1, 1], vec![1.0]);
        case.half_pupil = 0.5;
        case.gobo_size = [0.2, 0.2];
        case.distance = 0.4;
        case.focal_length = 0.4 / 1.4;
        let mut scene = case.scene();
        scene.camera.eye = Point3::new(0.0, 0.0, 0.75);
        scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
        scene
            .primitives
            .push(floor(1.0, 2.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources = [case.emitter(1.0, RigidTransform::identity())];
        let expected = 0.5 * case.lux([0.0, 0.0], 1.0) / PI;
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let policy = Settings {
                strategy,
                sampler: Sampler::Iid,
                ..settings(2, 65_536)
            };
            let value = pixel(&scene, &sources, PathMedium::Vacuum, policy, cx);
            let error = value.estimated_standard_error_cd_m2.unwrap();
            assert!(
                (value.luminance_cd_m2 - expected).abs() <= 6.0 * error + 0.003 * expected,
                "{strategy:?}: observed={} reference={expected} estimated_se={error}",
                value.luminance_cd_m2
            );
            close(
                value.luminance_cd_m2,
                expected,
                0.04,
                "independent projector surface estimators",
            );
        }
    });
}

fn single_scatter_reference(case: &ProjectorCase, cells: usize) -> f64 {
    let mut sum = 0.0;
    let da = 2.0 * case.half_pupil / cells as f64;
    let dz = 1.0 / (2 * cells) as f64;
    for iz in 0..2 * cells {
        let z = 0.5 + (iz as f64 + 0.5) * dz;
        for ix in 0..cells {
            let a = -case.half_pupil + (ix as f64 + 0.5) * da;
            for iy in 0..cells {
                let b = -case.half_pupil + (iy as f64 + 0.5) * da;
                let r = ((0.1 - a).powi(2) + b * b + z * z).sqrt();
                let medium_distance = (1.5 - z) + (z - 0.5) * r / z;
                // L*cos(theta)/r² = Phi*s²/(A*G) * r/z³.
                // The full clear gobo supports these pupil/event pairs.
                sum += r / z.powi(3) * (-0.6 * medium_distance).exp() * da * da * dz;
            }
        }
    }
    case.clear_flux * case.distance.powi(2)
        / (4.0 * case.half_pupil.powi(2) * case.gobo_size[0] * case.gobo_size[1])
        * 0.4
        / (4.0 * PI)
        * sum
}

#[test]
fn g1_projector_haze_matches_homogeneous_single_scattering_reference() {
    with_cx(|cx, _| {
        let mut case = ProjectorCase::new(3.0, [1, 1], vec![1.0]);
        case.gobo_size = [0.08, 0.08];
        let mut scene = case.scene();
        // The camera ray misses the pupil, so the one-vertex result contains
        // only source->haze->camera transport, with no source-hit radiance.
        scene.camera.eye = Point3::new(0.1, 0.0, 2.0);
        let sources = [case.emitter(1.0, RigidTransform::identity())];
        let density = [1.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 0.5], [4.0, 4.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        assert_eq!(absorption.bounds_m(), ([-2.0, -2.0, 0.5], [2.0, 2.0, 1.5]));
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(0.4).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let expected = single_scatter_reference(&case, 64);
        close(
            single_scatter_reference(&case, 32),
            expected,
            0.0003,
            "independent haze quadrature refinement",
        );
        let value = pixel(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            settings(1, 65_536),
            cx,
        );
        close(
            value.luminance_cd_m2,
            expected,
            0.035,
            "projector single scattering with both Beer legs",
        );
    });
}

#[test]
fn g1_projector_optical_frame_moves_with_original_rigid_mesh() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let frame = RigidTransform::try_new(
            [
                0.0,
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
                core::f64::consts::FRAC_1_SQRT_2,
            ],
            [0.3, -0.2, 0.4],
        )
        .unwrap();
        let placement = RigidTransform::try_new(
            [
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
            ],
            [7.0, -3.0, 4.0],
        )
        .unwrap();
        let h = case.half_pupil;
        let vertices = [[-h, -h], [h, -h], [h, h], [-h, h]]
            .map(|p| {
                let point = frame.transform_point(Point3::new(p[0], p[1], 0.0));
                [point.x, point.y, point.z]
            })
            .to_vec();
        let instance = GeometryInstance::try_new(
            71,
            fs_blake3::hash_domain("org.frankensim.test.projector", b"rigid-aperture"),
            SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2], [0, 2, 3]])),
            placement,
        )
        .unwrap();
        let mut scene = empty_scene([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        scene.primitives.push(Primitive {
            shape: Shape::Instance(instance),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let sources = [case.emitter(1.0, frame)];
        let local_point = Point3::new(0.3, 0.15, 3.0);
        let point = placement.transform_point(frame.transform_point(local_point));
        let normal = placement.transform_vector(frame.transform_vector(Vec3::new(0.0, 0.0, -1.0)));
        let result = probe(&scene, &sources, point, normal, 8192, cx);
        close(
            result.illuminance_lux,
            case.lux([0.3, 0.15], 3.0),
            0.003,
            "body-local optical frame under pan/tilt",
        );
    });
}

#[test]
fn g0_projector_jobs_refuse_nonplanar_reversed_and_unrepresentable_pupils() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [1, 1], vec![1.0]);
        let sources = [case.emitter(1.0, RigidTransform::identity())];
        for warped in [true, false] {
            let mut scene = case.scene();
            let Shape::Mesh(mesh) = &scene.primitives[0].shape else {
                unreachable!()
            };
            let mut vertices = mesh.vertices.clone();
            let triangles = if warped {
                vertices[2][2] = 0.001;
                vec![[0, 1, 2], [0, 2, 3]]
            } else {
                vec![[0, 2, 1], [0, 3, 2]]
            };
            scene.primitives[0].shape = Shape::Mesh(TriMesh::new(vertices, triangles));
            assert!(
                matches!(
                    PathLuminanceJob::new_with_emitters(
                        &scene,
                        PhotometricEmitters {
                            points: &[],
                            surfaces: &sources
                        },
                        PathMedium::Vacuum,
                        TransportOptions::default(),
                        settings(1, 8),
                        cx,
                    ),
                    Err(PhotometryError::Source(_))
                ),
                "warped={warped}"
            );
        }

        let both = PhotometricAreaEmitter::try_new(
            0,
            case.clear_flux,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Both,
        )
        .unwrap();
        assert!(
            both.with_projector(sources[0].projector().unwrap().clone())
                .is_err()
        );

        let mut extreme = case;
        extreme.clear_flux = 1.0e200;
        extreme.gobo_size = [1.0e-100; 2];
        let scene = extreme.scene();
        let sources = [extreme.emitter(1.0, RigidTransform::identity())];
        // Source flux, pupil area, base Lambertian luminance and the isolated
        // lens law are finite. Their focused radiance product is not.
        assert!((extreme.clear_flux / (4.0 * extreme.half_pupil.powi(2)) / PI).is_finite());
        assert!(matches!(
            PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &sources
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(1, 8),
                cx,
            ),
            Err(PhotometryError::Source(_))
        ));
    });
}

type ImageBits = Vec<Option<([u64; 3], u32, Option<u64>)>>;
fn image_bits(job: &PathLuminanceJob<'_, '_>) -> ImageBits {
    job.estimates()
        .map(|value| value.map(estimate_bits))
        .collect()
}

fn probe_bits(job: &PathIlluminanceJob<'_, '_>) -> ImageBits {
    job.estimates()
        .map(|value| {
            value.map(|v| {
                (
                    [v.direct_lux, v.continuous_lux, v.illuminance_lux].map(f64::to_bits),
                    v.samples,
                    v.estimated_standard_error_lux.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn screen_scene(case: &ProjectorCase) -> Scene {
    let mut scene = case.scene();
    scene
        .primitives
        .push(floor(3.0, 3.0, LiftedSpectrum { c: [0.0; 3] }));
    scene.camera.eye = Point3::new(0.0, 0.0, 2.0);
    scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
    scene.camera.half_tan = 0.4;
    scene
}

#[test]
fn g4_g5_projector_parallel_cancellation_and_mixed_resume_retain_bits() {
    with_cx(|cx, gate| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let scene = screen_scene(&case);
        let sources = [case.emitter(1.0, RigidTransform::identity())];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &sources,
        };
        let probes = [-0.3, 0.0, 0.3].map(|x| {
            IlluminanceProbe::try_new(Point3::new(x, 0.1, 3.0), Vec3::new(0.0, 0.0, -1.0)).unwrap()
        });
        let new_image = || {
            PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    width: 4,
                    height: 3,
                    sampler: Sampler::Iid,
                    ..settings(2, 32)
                },
                cx,
            )
            .unwrap()
        };
        let new_probes = || {
            PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 32,
                    max_depth: 1,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap()
        };
        let mut image_reference = new_image();
        let mut probe_reference = new_probes();
        image_reference.advance(cx, usize::MAX).unwrap();
        probe_reference.advance(cx, usize::MAX).unwrap();
        let mut image = new_image();
        let mut probes = new_probes();
        let execution =
            RenderExecutionConfig::try_new(2, 2, 2, 64 << 20, RunId(0x676f626f)).unwrap();
        image.advance(cx, 7).unwrap();
        probes.advance_with_execution(cx, 11, &execution).unwrap();
        let image_prefix = image_bits(&image);
        let probe_prefix = probe_bits(&probes);
        gate.request();
        assert!(matches!(
            image.advance_with_execution(cx, 17, &execution),
            Err(PhotometryError::Cancelled)
        ));
        assert!(matches!(
            probes.advance(cx, 17),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(image_bits(&image), image_prefix);
        assert_eq!(probe_bits(&probes), probe_prefix);
        with_cx(|resumed, _| {
            let execution =
                RenderExecutionConfig::try_new(3, 1, 4, 64 << 20, RunId(0x676f6270)).unwrap();
            image
                .advance_with_execution(resumed, 23, &execution)
                .unwrap();
            image.advance(resumed, usize::MAX).unwrap();
            probes.advance(resumed, 5).unwrap();
            probes
                .advance_with_execution(resumed, usize::MAX, &execution)
                .unwrap();
        });
        assert_eq!(image_bits(&image), image_bits(&image_reference));
        assert_eq!(probe_bits(&probes), probe_bits(&probe_reference));
    });
}

#[test]
#[ignore = "run explicitly to qualify finite-pupil projector CPU throughput"]
fn photometric_projector_parallel_throughput() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let scene = screen_scene(&case);
        let sources = [case.emitter(1.0, RigidTransform::identity())];
        let new_job = || {
            PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &sources,
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    width: 32,
                    height: 32,
                    sampler: Sampler::Iid,
                    ..settings(2, 128)
                },
                cx,
            )
            .unwrap()
        };
        let mut serial = new_job();
        let start = std::time::Instant::now();
        serial.advance(cx, usize::MAX).unwrap();
        let seconds = start.elapsed().as_secs_f64();
        let paths = serial.completed_path_samples();
        assert_eq!(paths, 131_072);
        assert!(
            serial
                .estimates()
                .flatten()
                .any(|pixel| pixel.luminance_cd_m2 > 0.0)
        );
        let expected = image_bits(&serial);
        eprintln!(
            "photometric_projector_throughput execution=serial paths={paths} seconds={seconds:.6} paths_per_second={:.1}",
            paths as f64 / seconds
        );
        for workers in [1, 4, 8] {
            let mut job = new_job();
            let execution = RenderExecutionConfig::try_new(
                8,
                8,
                workers,
                64 << 20,
                RunId(0x676f6200 + workers as u64),
            )
            .unwrap();
            let start = std::time::Instant::now();
            job.advance_with_execution(cx, usize::MAX, &execution)
                .unwrap();
            let seconds = start.elapsed().as_secs_f64();
            assert!(job.is_complete());
            assert_eq!(image_bits(&job), expected, "workers={workers}");
            eprintln!(
                "photometric_projector_throughput workers={workers} paths={paths} seconds={seconds:.6} paths_per_second={:.1}",
                paths as f64 / seconds
            );
        }
    });
}
