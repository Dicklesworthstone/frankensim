//! G1/G3: actual source histories, dispersive delays and original native paths.

mod analyzer;
mod exposure;

use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::camera::{RadianceSpectrum, SurfaceEmitter};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use asupersync::types::Budget;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;
use std::sync::Mutex;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 53,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn plane(id: u64, z: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        fs_blake3::hash_domain("retarded-native-camera", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-10.0, -10.0, z], [10.0, -10.0, z], [0.0, 10.0, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn material(ior: CauchyIor, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn opaque(object_id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::Opaque,
    }
}

fn source(object_id: u64, radiance: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        object_id,
        RadianceSpectrum::try_constant(radiance).unwrap(),
        EmissionSides::Both,
    )
}

fn interface() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::Dielectric {
            inside: 0,
            outside: 1,
            surface: DielectricSurface::SMOOTH,
        },
    }
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 10.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tol={tolerance:e}"
    );
}

struct History<F> {
    evaluate: F,
    queries: Mutex<Vec<SurfaceEmissionQuery>>,
}

impl<F> History<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            queries: Mutex::new(Vec::new()),
        }
    }
}

impl<F> SurfaceEmissionHistory for History<F>
where
    F: Fn(SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> + Sync,
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        (self.evaluate)(query)
    }
}

fn binding(id: u64, history: &dyn SurfaceEmissionHistory) -> SurfaceEmissionHistoryBinding<'_> {
    SurfaceEmissionHistoryBinding {
        object_id: id,
        source_identity: ContentHash([73; 32]),
        history,
    }
}

#[test]
fn g3_unit_histories_keep_complete_native_samples_and_random_event_addresses() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4), plane(9, 1.0)]).unwrap();
        let incident = material(CauchyIor::try_constant(1.0).unwrap(), 0.1);
        let transmitted = material(CauchyIor::try_constant(1.5).unwrap(), 0.4);
        let media = [
            SceneMedium::Homogeneous(&incident),
            SceneMedium::Homogeneous(&transmitted),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![interface(), opaque(9)],
            100_000,
        )
        .unwrap();
        let sources = [source(1, 2.0), source(9, 5.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let unit = History::new(|_| Ok(1.0));
        let timed = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &unit), binding(9, &unit)],
        )
        .unwrap();
        assert!(core::ptr::eq(timed.source_camera(), &camera));
        assert_eq!(timed.history_identity(9), Some(ContentHash([73; 32])));
        assert_eq!(timed.history_identity(404), None);
        for branch in [0.0, 0.99] {
            let mut steady_events = Vec::new();
            let steady = camera
                .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |event| {
                    steady_events.push(event);
                    [0.25, 0.5, branch]
                })
                .unwrap();
            let mut timed_events = Vec::new();
            let observed = timed
                .sample_ray_with_samples_at_time(
                    cx,
                    ray(),
                    550.0,
                    0,
                    Dd::from_f64(1.0),
                    settings(),
                    |event| {
                        timed_events.push(event);
                        [0.25, 0.5, branch]
                    },
                )
                .unwrap();
            assert_eq!(*observed.radiance_sample(), steady);
            assert_eq!(steady_events, timed_events);
        }
        for path in 0..8 {
            let steady = camera
                .sample_ray(cx, ray(), 550.0, 0, settings(), 29, path)
                .unwrap();
            let observed = timed
                .sample_ray_at_time(cx, ray(), 550.0, 0, Dd::from_f64(1.0), settings(), 29, path)
                .unwrap();
            assert_eq!(*observed.radiance_sample(), steady);
        }
    });
}

#[test]
fn g1_each_emitting_encounter_uses_its_own_delay_before_the_next_refraction() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4), plane(9, 1.0)]).unwrap();
        let incident = material(CauchyIor::try_constant(1.0).unwrap(), 0.1);
        let transmitted = material(CauchyIor::try_constant(1.5).unwrap(), 0.4);
        let media = [
            SceneMedium::Homogeneous(&incident),
            SceneMedium::Homogeneous(&transmitted),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![interface(), opaque(9)],
            100_000,
        )
        .unwrap();
        let sources = [source(1, 2.0), source(9, 5.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        // One actual source history is shared by two physically distinct emitters.
        let ramp =
            History::new(|query: SurfaceEmissionQuery| Ok(query.emission_time_s.to_f64() / 1e-8));
        let timed = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &ramp), binding(9, &ramp)],
        )
        .unwrap();
        let observation = 2e-8;
        let sample = timed
            .sample_ray_with_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                Dd::from_f64(observation),
                settings(),
                |_| [0.5, 0.5, 0.99],
            )
            .unwrap();
        let first_time = observation - 0.4 / LIGHT_M_S;
        let last_time = observation - (0.4 + 1.5 * 0.6) / LIGHT_M_S;
        let first = 2.0 * (-0.04_f64).exp();
        let last = 5.0 * (-0.28_f64).exp() / 1.5_f64.powi(2);
        near(
            sample.radiance_sample().radiance_w_m2_sr_nm(),
            first * first_time / 1e-8 + last * last_time / 1e-8,
            3e-7,
        );
        let queries = ramp.queries.lock().unwrap();
        assert_eq!(queries.len(), 2);
        assert_eq!([queries[0].object_id, queries[1].object_id], [1, 9]);
        near(queries[0].emission_time_s.to_f64(), first_time, 2e-16);
        near(queries[1].emission_time_s.to_f64(), last_time, 2e-16);
        assert!(queries[1].emission_time_s.lt(queries[0].emission_time_s));
        assert_eq!(sample.radiance_sample().path().surface_events(), 2);
        let incorrect_final_gate = (first + last) * last_time / 1e-8;
        assert!(
            (sample.radiance_sample().radiance_w_m2_sr_nm() - incorrect_final_gate).abs() > 0.1
        );
    });
}

#[test]
fn g1_chromatic_pulse_arrival_uses_group_index_instead_of_phase_or_vacuum_delay() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 2.0)]).unwrap();
        let glass = material(CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(), 0.0);
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let sources = [source(9, 3.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let pulse = History::new(|query: SurfaceEmissionQuery| {
            Ok((1.0 - (query.emission_time_s.to_f64() / 1e-10).abs()).max(0.0))
        });
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[binding(9, &pulse)]).unwrap();
        // Independently differentiate n=a+b/lambda_um²: ng=a+3b/lambda_um².
        let blue_group = 1.4 + 0.06 / 0.45_f64.powi(2);
        let time = Dd::from_f64(2.0 * blue_group / LIGHT_M_S);
        let blue = timed
            .sample_ray_at_time(cx, ray(), 450.0, 0, time, settings(), 1, 0)
            .unwrap();
        let red = timed
            .sample_ray_at_time(cx, ray(), 650.0, 0, time, settings(), 1, 0)
            .unwrap();
        near(blue.radiance_sample().radiance_w_m2_sr_nm(), 3.0, 3e-6);
        assert_eq!(red.radiance_sample().radiance_w_m2_sr_nm(), 0.0);
        let q = pulse.queries.lock().unwrap();
        near(
            q[0].group_delay_s.to_f64(),
            2.0 * blue_group / LIGHT_M_S,
            2e-16,
        );
        near(
            q[1].group_delay_s.to_f64(),
            2.0 * (1.4 + 0.06 / 0.65_f64.powi(2)) / LIGHT_M_S,
            2e-16,
        );
        let phase_time = 2.0 * (1.4 + 0.02 / 0.45_f64.powi(2)) / LIGHT_M_S;
        assert!((time.to_f64() - phase_time).abs() > 1e-10);
    });
}

#[test]
fn g1_original_curved_path_group_integral_controls_emission_clock_and_radiance() {
    struct LinearDispersive;
    impl GradientIndexField for LinearDispersive {
        fn sample(
            &self,
            p: Point3,
            wavelength: f64,
        ) -> Result<GradientIndexSample, GradientIndexError> {
            let dispersion = 0.02 / (wavelength * 1e-3).powi(2);
            let n = 1.2 + 0.1 * p.z + dispersion;
            GradientIndexSample::try_new(
                n,
                Vec3::new(0.0, 0.0, 0.1),
                0.2,
                Some(n + 2.0 * dispersion),
            )
        }
    }
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0)]).unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &LinearDispersive,
            identity: ContentHash([81; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let sources = [source(9, 7.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let gaussian = History::new(|query: SurfaceEmissionQuery| {
            Ok(-(query.emission_time_s.to_f64() / 1e-9).powi(2)).map(f64::exp)
        });
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[binding(9, &gaussian)]).unwrap();
        let (n0, n1, gradient) = (
            1.2 + 0.02 / 0.55_f64.powi(2),
            1.3 + 0.02 / 0.55_f64.powi(2),
            0.1,
        );
        let momentum = n0 * 0.4;
        let root = |n: f64| (n * n - momentum * momentum).sqrt();
        let delta = 0.04 / 0.55_f64.powi(2);
        // Exact conserved tangential optical momentum. ng=n+delta, so the
        // group integral has a different antiderivative from the phase path.
        let primitive = |n: f64| {
            0.5 * (n * root(n) + momentum * momentum * (n + root(n)).ln()) + delta * root(n)
        };
        let group_m = (primitive(n1) - primitive(n0)) / gradient;
        let distance = (root(n1) - root(n0)) / gradient;
        let x = momentum / gradient * ((n1 + root(n1)) / (n0 + root(n0))).ln();
        let launch = Ray {
            origin: ray().origin,
            dir: Vec3::new(0.4, 0.0, (1.0_f64 - 0.16).sqrt()),
        };
        let observed = timed
            .sample_ray_at_time(
                cx,
                launch,
                550.0,
                0,
                Dd::from_f64(group_m / LIGHT_M_S + 0.35e-9),
                settings(),
                3,
                0,
            )
            .unwrap();
        let sample = observed.radiance_sample();
        near(sample.path().beam().ray().origin.x, x, 3e-8);
        near(sample.path().beam().ray().dir.x, momentum / n1, 3e-9);
        near(sample.path().traveled_m(), distance, 3e-8);
        near(
            sample.path().group_optical_path_length_m().unwrap(),
            group_m,
            5e-8,
        );
        near(
            sample.radiance_w_m2_sr_nm(),
            7.0 * (n0 / n1).powi(2) * (-0.2 * distance - 0.35_f64.powi(2)).exp(),
            2e-6,
        );
        near(
            gaussian.queries.lock().unwrap()[0].emission_time_s.to_f64(),
            0.35e-9,
            2e-16,
        );
    });
}

#[test]
fn g3_common_large_epoch_preserves_short_optical_delays_and_source_waveforms() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 2.0)]).unwrap();
        let glass = material(CauchyIor::try_constant(1.5).unwrap(), 0.1);
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let sources = [source(9, 3.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let clock = Dd::from_f64(3.0 / LIGHT_M_S) + Dd::from_f64(0.2e-9);
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        assert_eq!((epoch + clock).to_f64(), epoch.to_f64());
        let original = History::new(|query: SurfaceEmissionQuery| {
            Ok(-(query.emission_time_s.to_f64() / 2e-9).powi(2)).map(f64::exp)
        });
        let shifted = History::new(|query: SurfaceEmissionQuery| {
            Ok(-((query.emission_time_s - epoch).to_f64() / 2e-9).powi(2)).map(f64::exp)
        });
        let first = TransientCurvedCameraScene::try_new(cx, &camera, &[binding(9, &original)])
            .unwrap()
            .sample_ray_at_time(cx, ray(), 550.0, 0, clock, settings(), 5, 0)
            .unwrap();
        let second = TransientCurvedCameraScene::try_new(cx, &camera, &[binding(9, &shifted)])
            .unwrap()
            .sample_ray_at_time(cx, ray(), 550.0, 0, epoch + clock, settings(), 5, 0)
            .unwrap();
        near(
            first.radiance_sample().radiance_w_m2_sr_nm(),
            second.radiance_sample().radiance_w_m2_sr_nm(),
            1e-13,
        );
        assert_eq!(
            first.radiance_sample().path(),
            second.radiance_sample().path()
        );
        assert_eq!(second.observation_time_s(), epoch + clock);
        let a = original.queries.lock().unwrap()[0];
        let b = shifted.queries.lock().unwrap()[0];
        near(
            (b.emission_time_s - epoch).to_f64(),
            a.emission_time_s.to_f64(),
            1e-22,
        );
    });
}

#[test]
fn g0_missing_group_index_invisible_sides_and_source_bindings_remain_explicit() {
    struct NoGroup;
    impl GradientIndexField for NoGroup {
        fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            GradientIndexSample::try_new(1.3, Vec3::new(0.0, 0.0, 0.0), 0.0, None)
        }
    }
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0)]).unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &NoGroup,
            identity: ContentHash([82; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let history = History::new(|_| Ok(1.0));
        for sides in [EmissionSides::Both, EmissionSides::Front] {
            let sources = [SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(2.0).unwrap(),
                sides,
            )];
            let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
            let empty = TransientCurvedCameraScene::try_new(cx, &camera, &[]).unwrap();
            let original = camera
                .sample_ray(cx, ray(), 550.0, 0, settings(), 3, 0)
                .unwrap();
            let timed = empty
                .sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, settings(), 3, 0)
                .unwrap();
            assert_eq!(*timed.radiance_sample(), original);
            let bound =
                TransientCurvedCameraScene::try_new(cx, &camera, &[binding(9, &history)]).unwrap();
            let result = bound.sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, settings(), 3, 0);
            if sides == EmissionSides::Both {
                assert!(matches!(
                    result,
                    Err(TransientCameraError::MissingGroupDelay { object_id: 9 })
                ));
            } else {
                assert_eq!(result.unwrap().radiance_sample().radiance_w_m2_sr_nm(), 0.0);
            }
            assert!(matches!(
                TransientCurvedCameraScene::try_new(cx, &camera, &[binding(10, &history)]),
                Err(TransientCameraError::InvalidHistoryBinding { object_id: 10 })
            ));
            assert!(
                TransientCurvedCameraScene::try_new(
                    cx,
                    &camera,
                    &[binding(9, &history), binding(9, &history)]
                )
                .is_err()
            );
            let mut unidentified = binding(9, &history);
            unidentified.source_identity = ContentHash([0; 32]);
            assert!(matches!(
                TransientCurvedCameraScene::try_new(cx, &camera, &[unidentified]),
                Err(TransientCameraError::MissingHistoryIdentity { object_id: 9 })
            ));
        }
        assert!(history.queries.lock().unwrap().is_empty());
        let unlit = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        assert!(matches!(
            TransientCurvedCameraScene::try_new(cx, &unlit, &[binding(9, &history)]),
            Err(TransientCameraError::InvalidHistoryBinding { object_id: 9 })
        ));
    });
}

#[test]
fn g4_bad_histories_unresolved_clocks_limits_and_late_cancellation_publish_no_sample() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4), plane(9, 1.0)]).unwrap();
        let glass = material(CauchyIor::try_constant(1.3).unwrap(), 0.0);
        let media = [
            SceneMedium::Homogeneous(&glass),
            SceneMedium::Homogeneous(&glass),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![interface(), opaque(9)],
            100_000,
        )
        .unwrap();
        let sources = [source(1, 2.0), source(9, 5.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let one = History::new(|_| Ok(1.0));
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[binding(1, &one), binding(9, &one)])
                .unwrap();
        for time in [
            Dd::from_f64(f64::NAN),
            Dd {
                hi: 1.0,
                lo: f64::INFINITY,
            },
        ] {
            assert!(matches!(
                timed.sample_ray_at_time(cx, ray(), 550.0, 0, time, settings(), 3, 0),
                Err(TransientCameraError::InvalidObservationTime)
            ));
        }
        // The clock guard (transient.rs checked_retarded_time) refuses an emission time
        // when recovering the group delay G from it misses by more than 64 eps G.
        // Closed form, through the double-double sum of fs_math::dd (two_sum, then
        // quick_two_sum), for observation (oh, lo) and G = (gh, gl), with gh and lo
        // below half an ulp of oh:
        // - the emission time is (oh, fl(lo - gh));
        // - the recovered delay is gh - d exactly (Sterbenz), where
        //   d = fl(lo - gh) - (lo - gh) = -e for (_, e) = two_sum(lo, -gh);
        // - so the miss is |gl - e|.
        // A zero lo word misses only |gl| <= eps/2 G, so 1e40 alone cannot fire the
        // guard (R-c: 4.68e-26 s against 2.46e-23 s for object 1). A lo word can:
        // half an ulp of lo passes 64 eps G once lo > 128 G.
        // G is the scene's own reported delay, read from a resolved probe. This test
        // checks the guard given a delay, not the delay.
        let probe = History::new(|_| Ok(1.0));
        TransientCurvedCameraScene::try_new(cx, &camera, &[binding(1, &probe), binding(9, &probe)])
            .unwrap()
            .sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, settings(), 3, 0)
            .unwrap();
        let delays: Vec<(u64, Dd)> = probe
            .queries
            .lock()
            .unwrap()
            .iter()
            .map(|query| (query.object_id, query.group_delay_s))
            .collect();
        assert_eq!(delays.iter().map(|(id, _)| *id).collect::<Vec<_>>(), [1, 9]);
        // The first contact in path order whose delay the closed form refuses, if any.
        let refused_by_closed_form = |lo: f64| {
            delays
                .iter()
                .find(|(_, delay)| {
                    let (_, e) = fs_math::eft::two_sum(lo, -delay.hi);
                    let miss = (Dd::from_f64(delay.lo) - Dd::from_f64(e)).abs();
                    (*delay * Dd::from_f64(64.0 * f64::EPSILON)).lt(miss)
                })
                .map(|(id, _)| *id)
        };
        let at = |lo: f64| {
            timed.sample_ray_at_time(cx, ray(), 550.0, 0, Dd { hi: 1e40, lo }, settings(), 3, 0)
        };
        assert_eq!(refused_by_closed_form(0.0), None);
        assert!(at(0.0).is_ok());
        // Plant one lo word on each side of the threshold. Start at a power of two
        // 2^k <= 127 G1: half an ulp of lo plus |gl| stays within the allowance there,
        // for G1 and for the larger G9. Then step up one binade at a time. Inside one
        // binade every lo shares one ulp and gh is fixed, so e is the same for all of
        // them, and the first binade the closed form refuses is "just past".
        let mut inside = 2f64.powi((127.0 * delays[0].1.hi).log2().floor() as i32);
        assert_eq!(refused_by_closed_form(inside), None);
        let (past, object_id) = loop {
            let next = 2.0 * inside;
            if let Some(object_id) = refused_by_closed_form(next) {
                break (next, object_id);
            }
            inside = next;
            assert!(inside < 1.0, "no lo word below 1 s reaches the clock guard");
        };
        assert!(
            matches!(at(past), Err(TransientCameraError::UnresolvedClock { object_id: id }) if id == object_id),
            "lo = {past:e} must refuse the clock of object {object_id}"
        );
        assert!(at(inside).is_ok(), "lo = {inside:e} must resolve");
        let mut limited = settings();
        limited.max_surface_events = 1;
        assert!(matches!(
            timed.sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, limited, 3, 0),
            Err(TransientCameraError::Camera(
                CameraRadianceError::Unresolved(PathTermination::SurfaceLimit)
            ))
        ));
        for value in [-1.0, f64::NAN, f64::INFINITY] {
            let invalid = History::new(|_| Ok(value));
            let bad = TransientCurvedCameraScene::try_new(
                cx,
                &camera,
                &[binding(1, &one), binding(9, &invalid)],
            )
            .unwrap();
            assert!(matches!(
                bad.sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, settings(), 3, 0),
                Err(TransientCameraError::InvalidHistoryScale { object_id: 9 })
            ));
        }
        let outside = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let bad = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &one), binding(9, &outside)],
        )
        .unwrap();
        assert!(matches!(
            bad.sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, settings(), 3, 0),
            Err(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::OutsideTimeDomain
            })
        ));
        assert_eq!(outside.queries.lock().unwrap().len(), 1);
        let cancel = History::new(|_| {
            gate.request();
            Ok(1.0)
        });
        let interrupted = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &one), binding(9, &cancel)],
        )
        .unwrap();
        assert!(
            interrupted
                .sample_ray_at_time(cx, ray(), 550.0, 0, Dd::ZERO, settings(), 3, 0)
                .is_err()
        );
        assert_eq!(cancel.queries.lock().unwrap().len(), 1);
        assert!(TransientCurvedCameraScene::try_new(cx, &camera, &[]).is_err());
    });
}
