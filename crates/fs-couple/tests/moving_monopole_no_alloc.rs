//! Warm reflected full fields and ideal pickup use only retained callback storage.
#![allow(unsafe_code)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! {
    static TRACKING: Cell<bool> = const { Cell::new(false) };
    static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}
struct CountingAllocator;
fn record() {
    if TRACKING.try_with(Cell::get).unwrap_or(false) {
        let _ = ALLOCATIONS.try_with(|count| count.set(count.get() + 1));
    }
}
// SAFETY: each operation delegates directly to System with unchanged arguments.
// Const-initialized thread-local counters allocate nothing and leave allocation
// behavior unchanged. Other test/harness threads do not affect the measurement.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record();
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record();
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        record();
        unsafe { System.realloc(ptr, layout, size) }
    }
}
#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, PreparedReflectedSourceTrajectory, PreparedSourceTrajectory,
    ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::observation::moving::field::{CompactFieldSample, MovingMonopoleRenderer};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;

#[test]
fn reflected_field_and_audio_callbacks_allocate_nothing_after_preparation() {
    let observer = MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        [12.0, 8.0, 3.0],
        MovingCompactLimits {
            maximum_mach: 0.1,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap();
    let knots = [-0.1, 0.5].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: [2.0 + 4.0 * source_time_s, 2.0 + 3.0 * source_time_s, 1.0],
            velocity_m_s: [4.0, 3.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        },
    });
    let trajectory = PreparedSourceTrajectory::new(
        observer,
        &knots,
        TrajectoryBudget {
            maximum_segments: 2,
            maximum_payload_bytes: 16_384,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
            maximum_root_iterations: 64,
        },
        &mut || false,
    )
    .unwrap();
    let walls = IdealReflectionObserver::new(
        observer,
        &[
            IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
            IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        ],
        4,
    )
    .unwrap();
    let reflected = PreparedReflectedSourceTrajectory::new(
        &trajectory,
        &walls,
        ReflectedTrajectoryBudget {
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
        },
        &mut || false,
    )
    .unwrap();
    let samples = [-0.1, 0.5].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength {
            volume_flow_m3_s: 1e-5,
            ..Default::default()
        },
    });
    let sources = [CompactSourceHistory {
        motion: &reflected,
        samples: &samples,
        source_bandwidth_hz: 0.0,
    }; 2];
    let mut renderer = MovingMonopoleRenderer::new(
        &walls,
        &sources,
        FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap(),
        MovingPressureConfig {
            first_receiver_time_s: 0.1,
            sample_rate_hz: 1000,
            samples: 128,
            max_block: 32,
            receiver_bandwidth_hz: 20.0,
            maximum_time_residual_s: 1e-12,
        },
        MovingPressureBudget {
            maximum_sources: 2,
            maximum_source_samples: 4,
            maximum_output_samples: 128,
            maximum_payload_bytes: 65_536,
            maximum_source_receiver_pairs: 8 * 128,
        },
        &mut || false,
    )
    .unwrap();
    let mut fields = [CompactFieldSample::default(); 32];
    let mut audio = [0.0; 32];
    // Warm all original path searches and the first-order pickup before counting.
    renderer.field_block(&mut fields).unwrap();
    ALLOCATIONS.with(|count| count.set(0));
    TRACKING.with(|enabled| enabled.set(true));
    let first = renderer.field_block(&mut fields);
    let second = PressureRenderer::block(&mut renderer, &mut audio);
    let third = renderer.field_block(&mut fields);
    TRACKING.with(|enabled| enabled.set(false));
    let count = ALLOCATIONS.with(Cell::get);
    first.unwrap();
    second.unwrap();
    third.unwrap();
    assert_eq!(count, 0, "warm physical field/audio callbacks allocated");
    assert_eq!(renderer.samples_rendered(), 128);
    assert_eq!(renderer.source_path_count(), 8);
    assert!(audio.iter().any(|pressure| pressure.abs() > 1e-8));
}
