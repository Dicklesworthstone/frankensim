//! Temporal support, state identity and independent exposure integrals.

use super::*;
use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use crate::volumes::transport::{AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 41,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn absorbing<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    absorption: f64,
    temperatures: Option<&'a [f64]>,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1, 1, density.len()], density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(absorption).unwrap(),
        temperatures,
        AbsorbingGridBudget {
            max_grid_cells: 16,
            max_steps_per_segment: 16,
        },
    )
    .unwrap()
}

fn timeline_budget() -> MediumTimelineBudget {
    MediumTimelineBudget { max_frames: 16 }
}

fn shutter(open: f64, duration: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        open,
        duration,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(-8.0, 8.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_held_snapshots_match_duration_weighted_beer_lambert_integral() {
    with_cx(|cx, _| {
        let density = [1.0];
        let first = absorbing(cx, &density, 0.3, None);
        let second = absorbing(cx, &density, 1.4, None);
        let frames = [
            MediumFrame {
                start_s: -2.0,
                end_s: -1.0,
                medium: PhysicalMedium::Absorbing(&first),
            },
            MediumFrame {
                start_s: -1.0,
                end_s: 2.0,
                medium: PhysicalMedium::Absorbing(&second),
            },
        ];
        let timeline = MediumTimeline::try_new(cx, &frames, timeline_budget()).unwrap();
        timeline.admit_shutter(cx, shutter(-2.0, 4.0)).unwrap();
        let mut transmitted = 0.0;
        const SAMPLES: usize = 4096;
        for sample in 0..SAMPLES {
            let time = -2.0 + 4.0 * (sample as f64 + 0.5) / SAMPLES as f64;
            let PhysicalMedium::Absorbing(medium) = timeline.at_time(time).unwrap() else {
                panic!("wrong optical model")
            };
            transmitted += medium
                .transmittance(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    Vec3::new(1.0, 0.0, 0.0),
                    f64::INFINITY,
                    &[550.0],
                )
                .unwrap()[0];
        }
        let expected = 0.25 * (-0.3_f64).exp() + 0.75 * (-1.4_f64).exp();
        assert!((transmitted / SAMPLES as f64 / expected - 1.0).abs() < 2e-13);
        // Interpolating opacity or averaging density would give a different
        // radiance, even though both endpoints themselves remain unchanged.
        let averaged_density_result = (-(0.25_f64 * 0.3 + 0.75 * 1.4)).exp();
        assert!((expected / averaged_density_result - 1.0).abs() > 0.1);
    });
}

#[test]
fn g0_discontinuities_choose_later_state_and_final_endpoint_is_inclusive() {
    with_cx(|cx, _| {
        let density = [1.0];
        let first = absorbing(cx, &density, 0.3, None);
        let second = absorbing(cx, &density, 0.7, None);
        let frames = [
            MediumFrame {
                start_s: -2.0,
                end_s: 0.0,
                medium: PhysicalMedium::Absorbing(&first),
            },
            MediumFrame {
                start_s: 0.0,
                end_s: 2.0,
                medium: PhysicalMedium::Absorbing(&second),
            },
        ];
        let timeline = MediumTimeline::try_new(cx, &frames, timeline_budget()).unwrap();
        assert_eq!(timeline.frame_count(), 2);
        for time in [-2.0, -1.0, -f64::MIN_POSITIVE] {
            assert_eq!(
                timeline.at_time(time).unwrap().content_identity(),
                first.content_identity()
            );
        }
        for time in [-0.0, 0.0, 1.0, 2.0] {
            assert_eq!(
                timeline.at_time(time).unwrap().content_identity(),
                second.content_identity()
            );
        }
        for time in [-2.0, 0.0, 2.0] {
            timeline.admit_shutter(cx, shutter(time, 0.0)).unwrap();
        }
        for time in [-2.1, 2.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                timeline.at_time(time),
                Err(TimelineError::TimeOutsideHistory)
            ));
        }
        assert_eq!(
            timeline.admit_shutter(cx, shutter(-3.0, 2.0)),
            Err(TimelineError::ShutterOutsideHistory)
        );
        assert_eq!(
            timeline.admit_shutter(cx, shutter(1.0, 2.0)),
            Err(TimelineError::ShutterOutsideHistory)
        );
    });
}

#[test]
fn g3_identity_binds_numeric_state_and_times_without_addresses_or_budget() {
    with_cx(|cx, _| {
        let first_density = [1.0];
        let copied_density = [1.0];
        let changed_density = [2.0];
        let first = absorbing(cx, &first_density, 0.4, None);
        let copy = absorbing(cx, &copied_density, 0.4, None);
        let changed = absorbing(cx, &changed_density, 0.4, None);
        let frame = |start_s, end_s, medium| MediumFrame {
            start_s,
            end_s,
            medium,
        };
        let first_frames = [frame(-0.0, 1.0, PhysicalMedium::Absorbing(&first))];
        let copied_frames = [frame(0.0, 1.0, PhysicalMedium::Absorbing(&copy))];
        let changed_frames = [frame(0.0, 1.0, PhysicalMedium::Absorbing(&changed))];
        let changed_times = [frame(0.0, 2.0, PhysicalMedium::Absorbing(&first))];
        let original = MediumTimeline::try_new(cx, &first_frames, timeline_budget()).unwrap();
        let duplicate =
            MediumTimeline::try_new(cx, &copied_frames, MediumTimelineBudget { max_frames: 1 })
                .unwrap();
        assert_eq!(original.content_identity(), duplicate.content_identity());
        assert_eq!(original.start_s().to_bits(), 0.0f64.to_bits());
        assert_ne!(
            original.content_identity(),
            MediumTimeline::try_new(cx, &changed_frames, timeline_budget())
                .unwrap()
                .content_identity()
        );
        assert_ne!(
            original.content_identity(),
            MediumTimeline::try_new(cx, &changed_times, timeline_budget())
                .unwrap()
                .content_identity()
        );
    });
}

#[test]
fn g0_model_and_grid_changes_are_explicit_snapshot_discontinuities() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0];
        let temperature = [1800.0, 1900.0];
        let vacuum = absorbing(cx, &density[..1], 0.0, None);
        let hot = absorbing(cx, &density, 0.3, Some(&temperature));
        let scattering = ScatteringGrid::try_new(
            cx,
            hot,
            MassScatteringSpectrum::constant(0.7).unwrap(),
            PhaseFunction::Rayleigh,
        )
        .unwrap();
        let frames = [
            MediumFrame {
                start_s: 0.0,
                end_s: 1.0,
                medium: PhysicalMedium::Absorbing(&vacuum),
            },
            MediumFrame {
                start_s: 1.0,
                end_s: 2.0,
                medium: PhysicalMedium::Scattering(&scattering),
            },
        ];
        let timeline = MediumTimeline::try_new(cx, &frames, timeline_budget()).unwrap();
        assert!(!timeline.is_vacuum());
        assert!(timeline.has_thermal_emission());
        let PhysicalMedium::Absorbing(first) = timeline.at_time(0.5).unwrap() else {
            panic!("interpolated model")
        };
        let PhysicalMedium::Scattering(second) = timeline.at_time(1.0).unwrap() else {
            panic!("interpolated model")
        };
        assert!(core::ptr::eq(first, &vacuum));
        assert!(core::ptr::eq(second, &scattering));
        assert!(PhysicalMedium::Absorbing(first).is_vacuum());
        assert!(PhysicalMedium::Scattering(second).has_thermal_emission());
    });
}

#[test]
fn g4_missing_ambiguous_or_nonfinite_intervals_refuse_at_admission() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = absorbing(cx, &density, 0.3, None);
        let frame = |start_s, end_s| MediumFrame {
            start_s,
            end_s,
            medium: PhysicalMedium::Absorbing(&medium),
        };
        let first = frame(0.0, 1.0);
        assert_eq!(
            MediumTimeline::try_new(cx, &[], timeline_budget()).unwrap_err(),
            TimelineError::Empty
        );
        for (start, end) in [
            (1.0, 1.0),
            (1.0, 0.0),
            (f64::NAN, 2.0),
            (0.0, f64::INFINITY),
        ] {
            assert_eq!(
                MediumTimeline::try_new(cx, &[frame(start, end)], timeline_budget()).unwrap_err(),
                TimelineError::InvalidInterval { index: 0 }
            );
        }
        assert_eq!(
            MediumTimeline::try_new(cx, &[first, frame(1.0001, 2.0)], timeline_budget())
                .unwrap_err(),
            TimelineError::Gap { index: 1 }
        );
        assert_eq!(
            MediumTimeline::try_new(cx, &[first, frame(0.9999, 2.0)], timeline_budget())
                .unwrap_err(),
            TimelineError::Overlap { index: 1 }
        );
        assert_eq!(
            MediumTimeline::try_new(cx, &[frame(1.0, 2.0), first], timeline_budget()).unwrap_err(),
            TimelineError::Overlap { index: 1 }
        );
    });
}

#[test]
fn g4_frame_budget_and_cancellation_refuse_without_partial_history() {
    with_cx(|cx, gate| {
        let density = [1.0];
        let medium = absorbing(cx, &density, 0.0, None);
        let frames = [MediumFrame {
            start_s: 0.0,
            end_s: 1.0,
            medium: PhysicalMedium::Absorbing(&medium),
        }];
        assert_eq!(
            MediumTimeline::try_new(cx, &frames, MediumTimelineBudget { max_frames: 0 })
                .unwrap_err(),
            TimelineError::Budget
        );
        let timeline = MediumTimeline::try_new(cx, &frames, timeline_budget()).unwrap();
        assert!(timeline.is_vacuum());
        assert!(!timeline.has_thermal_emission());
        gate.request();
        assert_eq!(
            timeline.admit_shutter(cx, shutter(0.0, 1.0)),
            Err(TimelineError::Cancelled)
        );
        assert_eq!(
            MediumTimeline::try_new(cx, &frames, timeline_budget()).unwrap_err(),
            TimelineError::Cancelled
        );
    });
}
