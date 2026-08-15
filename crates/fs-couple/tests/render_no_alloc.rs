//! No-allocation gate for the block render hot path (music bead
//! `frankensim-music-v8-root-3ez8g.2.1`).
//!
//! A counting global allocator wraps `System`; after construction and one
//! warm-up block, rendering through the ADMITTED no-alloc assembly (the
//! massless-reed voice with an empty plate bank — scalar reed math plus
//! the exact-FIR ring buffer) must perform ZERO heap allocations per
//! block. This lives in its own test binary because a `#[global_allocator]`
//! is process-wide.
//!
//! Disclosed allocating voices (NOT gated here, by design): the
//! massive-reed lay path (`dissipative_modal_forces` returns a `Vec` per
//! sample) and the exact-ZOH modal voice (`ModalAcousticTimeModel::step`
//! builds a per-sample energy frame). Both are fusion candidates (bead
//! 3ez8g.15); admitting them to the no-alloc set without fixing the
//! kernels would make this gate a lie.

// A counting global allocator cannot exist without `unsafe impl
// GlobalAlloc`; the two blocks below forward verbatim to `System` and
// only bump a relaxed counter. This is test-binary-only code (the
// unsafe-capsule scanner covers crates/*/src; production stays
// deny(unsafe_code)).
#![allow(unsafe_code)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicU64, Ordering};

struct CountingAllocator;

static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);

// SAFETY: delegates directly to `System`; the counter is a relaxed atomic
// with no side effects on allocation behavior.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: forwarded verbatim.
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: forwarded verbatim.
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: forwarded verbatim.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static COUNTER: CountingAllocator = CountingAllocator;

use fs_couple::render::{ReedBoreVoice, RenderContext, RenderVoice};
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;

#[test]
fn admitted_no_alloc_assembly_allocates_nothing_per_block() {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: 0.50,
        }],
    };
    let reed = BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.013,
        closing_pressure_pa: 6_000.0,
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0, // massless: the admitted no-alloc reed path
        stiffness_n_m: 0.0,
    };
    let voice = ReedBoreVoice::new(
        &duct,
        &air,
        reed,
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        48_000,
        4_800,
        None,
    )
    .expect("voice admits");
    let mut context = RenderContext::new(vec![RenderVoice::ReedBore(voice)], 512);
    let mut block = vec![0.0; 512];
    // Warm-up: any lazy one-time allocation happens here, outside the gate.
    context.block(&mut block).expect("warm-up block");

    let before = ALLOCATIONS.load(Ordering::Relaxed);
    for _ in 0..32 {
        context.block(&mut block).expect("gated block");
    }
    let after = ALLOCATIONS.load(Ordering::Relaxed);
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"render-no-alloc\",\"blocks\":32,\
         \"allocations\":{}}}",
        after - before
    );
    assert_eq!(
        after - before,
        0,
        "the admitted no-alloc assembly allocated {} times across 32 blocks; find the \
         allocation and either remove it or demote the voice to the disclosed-allocating set",
        after - before
    );
    // The gate must not be vacuous: the stream is real sound.
    let rms = (block.iter().map(|p| p * p).sum::<f64>() / block.len() as f64).sqrt();
    assert!(rms > 1.0, "non-vacuous stream (rms {rms})");
}
