//! G1/G3: owned exact runtime retains full damping/memory with zero warm allocations.
// Test-only allocator forwards to System, as in fs-couple's render_no_alloc test.
#![allow(unsafe_code)]
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh, split_quadratic, step_render};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage, duffing_oscillator};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! { static TRACK: Cell<(bool, usize)> = const { Cell::new((false, 0)) }; }
struct Counter;
fn count() {
    let _ = TRACK.try_with(|cell| {
        let (enabled, n) = cell.get();
        if enabled {
            cell.set((true, n + 1));
        }
    });
}
// SAFETY: all allocations and layouts are forwarded unchanged to System. The
// thread-local counter affects only test observation, never allocation behavior.
unsafe impl GlobalAlloc for Counter {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        count();
        unsafe { System.realloc(ptr, layout, size) }
    }
}
#[global_allocator]
static ALLOCATOR: Counter = Counter;

struct Datum {
    storage: QuadraticStorage,
    datum: f64,
}
impl Storage for Datum {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.datum + self.storage.hamiltonian(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.storage.gradient(x, out)
    }
    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        self.storage.quadratic_part()
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], m: &[f64]) -> Option<f64> {
        self.storage.midpoint_energy_remainder(a, b, m)
    }
}
fn system(datum: f64) -> PortHamiltonian {
    let q = vec![
        4.0, 0.0, 0.6, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    r[5] = 0.4;
    r[7] = -0.1;
    r[13] = -0.1;
    r[15] = 0.2;
    PortHamiltonian::new(
        4,
        2,
        j,
        r,
        vec![0.0, 0.0, 1.0, 0.2, 0.0, 0.0, -0.3, 1.0],
        Box::new(Datum {
            storage: QuadraticStorage::new(q, 4).unwrap(),
            datum,
        }),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 3.5,
        relaxation_time_s: 0.15,
    }])
    .unwrap()
}
fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 5,
        max_omega_dt: 0.04,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-12 * (1.0 + a.abs().max(b.abs())),
        "{a} versus {b}"
    );
}

#[test]
fn owning_and_borrowed_flows_preserve_complete_work_and_warm_steps_allocate_nothing() {
    let original = system(0.0);
    let borrowed = PreparedExactZoh::new(&original, 0.003, budget()).unwrap();
    // Owning constructor permits returning/moving the complete physical owner.
    let make = || PreparedExactZoh::new_owned(system(0.0), 0.003, budget()).unwrap();
    let mut owned = make();
    assert_eq!(owned.system().state_dim(), 5);
    assert_eq!(owned.system().port_dim(), 2);
    assert_eq!(owned.duration(), 0.003);
    let mut x = [0.3, -0.1, -0.2, 0.4, 0.07];
    let mut next = [0.0; 5];
    let mut output = [0.0; 2];
    for sample in 0..40 {
        let input = [0.2 * (sample as f64 * 0.17).sin(), -0.3];
        let reference = borrowed.step(&x, &input).unwrap();
        let actual = owned.step_into(&x, &input, &mut next, &mut output).unwrap();
        assert_eq!(reference.x, next);
        assert_eq!(reference.y, output);
        assert_eq!(reference.delta_h, actual.delta_h);
        assert_eq!(reference.dissipated, actual.dissipated);
        assert_eq!(reference.supplied, actual.supplied);
        assert!(actual.balance_residual().abs() < 3e-12);
        x = next;
    }
    let before = owned.system().hamiltonian(&x);
    let mut heat = 0.0;
    let mut supplied = 0.0;
    TRACK.with(|cell| cell.set((true, 0)));
    for _ in 0..1000 {
        let record = owned
            .step_into(&x, &[0.17, -0.23], &mut next, &mut output)
            .unwrap();
        heat += record.dissipated;
        supplied += record.supplied;
        x = next;
    }
    let allocations = TRACK.with(|cell| {
        let (_, n) = cell.get();
        cell.set((false, n));
        n
    });
    assert_eq!(
        allocations, 0,
        "successful warm propagation must reuse every numerical buffer"
    );
    assert!(heat > 0.01);
    close(owned.system().hamiltonian(&x) - before + heat, supplied);
}

#[test]
fn cached_full_damping_and_memory_match_direct_exact_flow_and_energy_datum() {
    let original = system(0.0);
    let split = split_quadratic(&original, budget()).unwrap();
    assert!(split.is_exact_quadratic());
    let mut cached = PreparedExactZoh::new_owned(system(0.0), 0.1, budget()).unwrap();
    let mut shifted = PreparedExactZoh::new_owned(system(1e24), 0.1, budget()).unwrap();
    let x = [0.3, -0.1, -0.2, 0.4, 0.07];
    let input = [0.2, -0.3];
    let direct = step_render(&original, &split, &x, &input, 0.1).unwrap();
    let mut end = [0.0; 5];
    let mut shifted_end = [0.0; 5];
    let mut out = [0.0; 2];
    let mut shifted_out = [0.0; 2];
    let a = cached.step_into(&x, &input, &mut end, &mut out).unwrap();
    let b = shifted
        .step_into(&x, &input, &mut shifted_end, &mut shifted_out)
        .unwrap();
    assert_eq!(end, shifted_end);
    assert_eq!(out, shifted_out);
    assert_eq!(a, b);
    for (actual, truth) in end.iter().zip(&direct.x) {
        close(*actual, *truth);
    }
    for (actual, truth) in out.iter().zip(&direct.y) {
        close(*actual, *truth);
    }
    close(a.delta_h, direct.delta_h);
    close(a.dissipated, direct.dissipated);
    close(a.supplied, direct.supplied);
}

#[test]
fn complete_admission_and_failed_candidates_preserve_both_caller_outputs() {
    let nonlinear = duffing_oscillator(1.0, 1.0, 2.0, 0.1).unwrap();
    assert!(
        !split_quadratic(&nonlinear, budget())
            .unwrap()
            .is_exact_quadratic()
    );
    assert!(PreparedExactZoh::new_owned(nonlinear, 0.001, budget()).is_err());
    let too_small = LinearSplitBudget {
        max_modes: 4,
        ..budget()
    };
    assert!(PreparedExactZoh::new_owned(system(0.0), 0.01, too_small).is_err());
    let mut cached = PreparedExactZoh::new_owned(system(0.0), 0.1, budget()).unwrap();
    let x = [0.3, -0.1, -0.2, 0.4, 0.07];
    let mut end = [71.0; 5];
    let mut out = [83.0; 2];
    assert!(
        cached
            .step_into(&x, &[f64::NAN, 0.0], &mut end, &mut out)
            .is_err()
    );
    assert_eq!(end, [71.0; 5]);
    assert_eq!(out, [83.0; 2]);
    assert!(
        cached
            .step_into(&x, &[0.0, 0.0], &mut end[..4], &mut out)
            .is_err()
    );
    assert_eq!(end, [71.0; 5]);
    assert_eq!(out, [83.0; 2]);
    assert!(
        cached
            .step_into(&[f64::MAX; 5], &[f64::MAX; 2], &mut end, &mut out)
            .is_err()
    );
    assert_eq!(end, [71.0; 5]);
    assert_eq!(out, [83.0; 2]);
    let original = system(0.0);
    let expected = original.output(&x);
    let mut zero = PreparedExactZoh::new_owned(original, 0.0, budget()).unwrap();
    let record = zero.step_into(&x, &[0.0, 0.0], &mut end, &mut out).unwrap();
    assert_eq!(end, x);
    assert_eq!(out.as_slice(), expected.as_slice());
    assert_eq!(
        (record.delta_h, record.dissipated, record.supplied),
        (0.0, 0.0, 0.0)
    );
}
