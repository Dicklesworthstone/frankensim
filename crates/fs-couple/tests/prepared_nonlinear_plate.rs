//! The prepared body preserves physical dynamics while reusing fixed-clock work.
// The allocator forwards directly to System and is confined to this test binary.
#![allow(unsafe_code)]
use fs_couple::thin_plate::VkBody;
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_phs::RelaxationBranch;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! {
    static TRACK: Cell<bool> = const { Cell::new(false) };
    static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}
struct Counter;
fn count() {
    let _ = TRACK.try_with(|on| {
        if on.get() {
            let _ = ALLOCATIONS.try_with(|n| n.set(n.get() + 1));
        }
    });
}
// SAFETY: all allocation operations are forwarded unchanged to System.
unsafe impl GlobalAlloc for Counter {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        count();
        unsafe { System.realloc(ptr, layout, size) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
}
#[global_allocator]
static COUNTER: Counter = Counter;

fn body(nonlinear: bool) -> VkBody {
    VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![2.0, 3.0],
            channels: if nonlinear {
                vec![StressChannel {
                    coefficient: 30.0,
                    coupling: vec![1.0, 0.2, 0.2, 0.6],
                }]
            } else {
                vec![]
            },
        },
        &[0.4, -0.1, -0.1, 0.2],
        &[1.0, 0.3],
        vec![0.02, -0.01],
    )
    .unwrap()
    .with_relaxation_branches(
        vec![RelaxationBranch {
            projection: vec![1.0, 0.0, -0.5, 0.0],
            stiffness: 3.5,
            relaxation_time_s: 0.15,
        }],
        1,
    )
    .unwrap()
}

#[test]
fn prepared_body_preserves_motion_history_and_independent_physical_work() {
    for nonlinear in [false, true] {
        let dt = 0.001;
        let mut reference = body(nonlinear);
        let mut original = body(nonlinear);
        for _ in 0..7 {
            reference.drive_and_radiate(2.0, dt, 1.2, 1.5).unwrap();
            original.drive_and_radiate(2.0, dt, 1.2, 1.5).unwrap();
        }
        let initial = original.state().to_vec();
        let initial_energy = original.stored_energy_j();
        let mut prepared = original.into_prepared(dt).unwrap();
        assert_eq!(prepared.state(), initial);
        assert_eq!(prepared.stored_energy_j(), initial_energy);
        assert_eq!(prepared.time_step_s(), dt);
        for i in 0..150 {
            let force = if i < 50 {
                0.8
            } else if i < 100 {
                -0.3
            } else {
                0.0
            };
            let before = prepared.state().to_vec();
            let energy = prepared.stored_energy_j();
            let frame = prepared.drive_and_radiate(force, 1.2, 1.5).unwrap();
            let (pressure, flow) = reference.drive_and_radiate(force, dt, 1.2, 1.5).unwrap();
            for (actual, expected) in prepared.state().iter().zip(reference.state()) {
                assert!(
                    (actual - expected).abs() < 2e-11,
                    "nonlinear={nonlinear}: {actual} {expected}"
                );
            }
            assert!((frame.pressure_pa - pressure).abs() < 2e-11);
            assert!((frame.volume_velocity_m3_s - flow).abs() < 2e-12);
            assert_eq!(frame.volume_velocity_m3_s, prepared.volume_velocity());
            let source_work = force
                * ((prepared.state()[0] - before[0]) + 0.3 * (prepared.state()[2] - before[2]));
            assert!((frame.work.supplied - source_work).abs() < 2e-12);
            assert!((frame.work.delta_h - (prepared.stored_energy_j() - energy)).abs() < 2e-12);
            assert!(frame.work.dissipated >= -2e-14);
            assert!(frame.work.balance_residual().abs() < 2e-12);
        }
        assert!(prepared.state()[4].abs() > 1e-5);
    }
}

#[test]
fn warm_complete_linear_and_nonlinear_material_steps_allocate_nothing() {
    for nonlinear in [false, true] {
        let mut prepared = body(nonlinear).into_prepared(0.001).unwrap();
        prepared.drive_and_radiate(1.0, 1.2, 1.0).unwrap();
        ALLOCATIONS.with(|n| n.set(0));
        TRACK.with(|on| on.set(true));
        for _ in 0..128 {
            let frame = prepared.drive_and_radiate(1.0, 1.2, 1.0).unwrap();
            std::hint::black_box(frame);
        }
        TRACK.with(|on| on.set(false));
        let allocations = ALLOCATIONS.with(Cell::get);
        assert_eq!(allocations, 0, "nonlinear={nonlinear}");
    }
}

#[test]
fn fixed_clock_and_cancelled_or_invalid_observation_never_reset_or_publish_state() {
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(body(false).into_prepared(dt).is_err());
    }
    assert!(body(true).into_prepared(0.1).is_err());
    // The complete linear exponential does not inherit the nonlinear phase cap.
    assert!(body(false).into_prepared(0.1).is_ok());
    for nonlinear in [false, true] {
        let mut prepared = body(nonlinear).into_prepared(0.001).unwrap();
        prepared.drive_and_radiate(1.0, 1.2, 1.0).unwrap();
        let before = prepared.state().to_vec();
        assert!(prepared.drive_and_radiate(1.0, -1.0, 1.0).is_err());
        assert_eq!(prepared.state(), before);
        assert!(prepared.drive_and_radiate(f64::NAN, 1.2, 1.0).is_err());
        assert_eq!(prepared.state(), before);
        for stop in [1, 2] {
            let mut calls = 0;
            let result = prepared.drive_and_radiate_controlled(1.0, 1.2, 1.0, || {
                calls += 1;
                calls == stop
            });
            assert!(result.is_err());
            assert_eq!(prepared.state(), before);
        }
        // Refused scratch cannot poison the next accepted material step.
        let mut same = body(nonlinear).into_prepared(0.001).unwrap();
        same.drive_and_radiate(1.0, 1.2, 1.0).unwrap();
        let expected = same.drive_and_radiate(0.7, 1.2, 1.0).unwrap();
        let actual = prepared.drive_and_radiate(0.7, 1.2, 1.0).unwrap();
        assert_eq!(prepared.state(), same.state());
        assert_eq!(actual, expected);
    }
}
