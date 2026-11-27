//! Zero-allocation proof for `set_input`, `step`, `read_probes`, and `render_f32`,
//! plus `from_system` binding validation and `render_f32` non-finite overflow silencing.

use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};

use fs_ports::PortKind;
use fs_runtime::{AssemblyRuntime, InputBinding, InputDescriptor, ProbeBinding, ProbeDescriptor};
use fs_system::elements::{ControlledSource, Mass, Motion, SourceAccount};
use fs_system::{NodeId, System, SystemBuilder};

struct Counting;
static TRACK: AtomicBool = AtomicBool::new(false);
static ALLOCS: AtomicUsize = AtomicUsize::new(0);
fn note() {
    if TRACK.load(Ordering::Relaxed) {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
    }
}
#[allow(unsafe_code)]
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        note();
        unsafe { SystemAllocator.alloc(l) }
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        note();
        unsafe { SystemAllocator.alloc_zeroed(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { SystemAllocator.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        note();
        unsafe { SystemAllocator.realloc(p, l, n) }
    }
}
#[global_allocator]
static ALLOCATOR: Counting = Counting;

fn osc_sys(cell: Arc<AtomicU64>) -> (System, NodeId) {
    let mut b = SystemBuilder::new();
    let m0 = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let m = b
        .add_with_tangent("m", Mass::new(1.0).unwrap(), m0)
        .unwrap();
    let src = ControlledSource::effort(PortKind::MechanicalForceVelocity);
    let f = b
        .add_with_tangent("f", src, SourceAccount::default())
        .unwrap();
    let (mp, fp, fc) = (
        b.port(m, "p").unwrap(),
        b.port(f, "p").unwrap(),
        b.port(f, "command").unwrap(),
    );
    b.common_flow(mp, &[fp]).unwrap();
    b.boundary(
        fc,
        Box::new(move |_, _, o| o[0] = f64::from_bits(cell.load(Ordering::Relaxed))),
    )
    .unwrap();
    (b.seal().unwrap(), m)
}

#[test]
fn hot_path_zero_alloc_and_binding_and_overflow_guards() {
    let cell = Arc::new(AtomicU64::new(999.0_f64.to_bits()));
    let (sys, m) = osc_sys(Arc::clone(&cell));
    let in0 = InputDescriptor {
        index: 0,
        name: "force".into(),
        target: "m".into(),
        port: "p".into(),
        unit: "N".into(),
        dims: [1, 1, -2, 0, 0, 0],
        width: 1,
        default_value: 0.0,
    };
    let pr0 = ProbeDescriptor {
        index: 0,
        name: "pos".into(),
        kind: "field-quantity".into(),
        target: "m".into(),
        quantity: "position".into(),
        unit: "m".into(),
        dims: [1, 0, 0, 0, 0, 0],
        offset: 0,
        len: 1,
        shape: Vec::new(),
    };
    let ins = vec![InputBinding {
        descriptor: in0.clone(),
        cell: Arc::clone(&cell),
    }];
    let prs = vec![ProbeBinding {
        descriptor: pr0.clone(),
        reader: Box::new(move |s, d| d[0] = s.state::<Motion>(m).unwrap().position_m),
    }];
    let mut rt = AssemblyRuntime::from_system("no_alloc_osc", sys, ins, prs, 7).unwrap();
    assert_eq!(cell.load(Ordering::Relaxed), 0.0_f64.to_bits());
    rt.step(1.0e-3).unwrap();
    let (mut probe_buf, mut audio_buf) = ([0.0_f64; 1], [0.0_f32; 128]);
    ALLOCS.store(0, Ordering::SeqCst);
    TRACK.store(true, Ordering::SeqCst);
    rt.set_input(0, 5.0).unwrap();
    for _ in 0..32 {
        rt.step(1.0e-3).unwrap();
        rt.read_probes(&mut probe_buf).unwrap();
    }
    assert_eq!(rt.render_f32(0, 1.0 / 48_000.0, 2, &mut audio_buf), 128);
    TRACK.store(false, Ordering::SeqCst);
    assert_eq!(ALLOCS.load(Ordering::SeqCst), 0);
    assert!(audio_buf[127].abs() > 0.0);

    let p3 = ProbeBinding {
        descriptor: pr0,
        reader: Box::new(|_, d| d[0] = 1.0e40),
    };
    let mut ovf = AssemblyRuntime::from_system("a", osc_sys(cell).0, vec![], vec![p3], 1).unwrap();
    let mut out = [1.0_f32; 4];
    assert_eq!(ovf.render_f32(0, 1.0e-3, 1, &mut out), 0);
    assert!(out.iter().all(|s| s.to_bits() == 0) && ovf.is_halted());
    assert_eq!(ovf.take_refusal().unwrap().code, "non-finite-probe");
}
