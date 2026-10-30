//! G1/G3: selected-port wave closure retains original motion, memory and work.
// Test-only allocation observation follows prepared_exact_owned.rs.
#![allow(unsafe_code)]
use fs_phs::linear_exact::{LinearSplitBudget, PreparedScatteringZoh};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, duffing_oscillator};
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
// SAFETY: layouts and pointers pass unchanged to the system allocator. Only a
// thread-local test counter is added; production allocation behavior is untouched.
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

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.04,
    }
}

fn body() -> PortHamiltonian {
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    let mut r = vec![0.0; 16];
    // Positive definite, non-diagonally-dominant physical material damping.
    r[5] = 0.4;
    r[7] = 0.5;
    r[13] = 0.5;
    r[15] = 0.8;
    let q = vec![
        4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    PortHamiltonian::new(
        4,
        3,
        j,
        r,
        vec![
            0.0, 0.0, 0.0, 1.0, 0.7, -0.2, 0.0, 0.0, 0.0, -0.3, -0.4, 1.1,
        ],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, 0.0, -0.5, 0.0],
        stiffness: 4.0,
        relaxation_time_s: 0.15,
    }])
    .unwrap()
}

fn rk4<const N: usize>(
    mut x: [f64; N],
    time: f64,
    steps: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = time / steps as f64;
    for _ in 0..steps {
        let a = rhs(&x);
        let b = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(&core::array::from_fn(|i| x[i] + dt * c[i]));
        for i in 0..N {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn selected_wave_preserves_full_material_history_and_all_external_port_work() {
    let initial = [0.02, 0.1, -0.03, 0.04, 0.03];
    let inputs = [0.3, 0.07, -0.2];
    let z = 2.5;
    let dt = 0.037;
    let original = body();
    let structure = original.structure();
    let expected = (
        structure.0.to_vec(),
        structure.1.to_vec(),
        structure.2.to_vec(),
    );
    let mut runtime = PreparedScatteringZoh::new_owned(original, 1, z, dt, budget(5)).unwrap();
    assert_eq!(
        runtime.system().structure(),
        (
            expected.0.as_slice(),
            expected.1.as_slice(),
            expected.2.as_slice()
        )
    );
    assert_eq!(runtime.selected_port(), 1);
    let (free, beta) = runtime.average_relation(&initial, &inputs).unwrap();
    let (unforced, same_beta) = runtime
        .average_relation(&initial, &[0.0, 10.0, 0.0])
        .unwrap();
    assert_eq!(beta.to_bits(), same_beta.to_bits());
    assert!((free - unforced).abs() > 1e-5);
    let step = runtime.step(&initial, &inputs).unwrap();
    assert!((step.y[1] - (free + beta * inputs[1])).abs() < 3e-13);
    // Independent original physical equations: pressure/flow interpretation is
    // not assumed by the shared owner, and this mechanical port uses y/u units.
    let reference = rk4(
        [
            initial[0], initial[1], initial[2], initial[3], initial[4], 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0,
        ],
        dt,
        2048,
        |x| {
            let arm = 2.0 * x[0] - x[2] - x[4];
            let y = 0.7 * x[1] - 0.4 * x[3];
            let u = (2.0 * inputs[1] - y) / z;
            let heat = 0.4 * x[1] * x[1] + x[1] * x[3] + 0.8 * x[3] * x[3] + arm * arm / 0.15;
            let other = inputs[0] * (x[1] - 0.3 * x[3]) + inputs[2] * (-0.2 * x[1] + 1.1 * x[3]);
            [
                x[1],
                -4.0 * x[0] - 0.4 * x[1] - 0.5 * x[3] - 2.0 * arm + inputs[0] + 0.7 * u
                    - 0.2 * inputs[2],
                x[3],
                -9.0 * x[2] - 0.5 * x[1] - 0.8 * x[3] + arm - 0.3 * inputs[0] - 0.4 * u
                    + 1.1 * inputs[2],
                arm / 0.15,
                heat,
                other,
                u * y,
                inputs[1] * inputs[1] / z,
                (y - inputs[1]).powi(2) / z,
                y,
            ]
        },
    );
    for i in 0..5 {
        assert!((step.x[i] - reference[i]).abs() < 3e-12);
    }
    let record = step.account;
    for (actual, expected) in [
        (record.physical_heat, reference[5]),
        (record.external_work, reference[6]),
        (record.boundary_work, reference[7]),
        (record.incident_energy, reference[8]),
        (record.outgoing_energy, reference[9]),
        (step.y[1] * dt, reference[10]),
    ] {
        assert!(
            (actual - expected).abs() < 5e-13,
            "{actual:e} vs {expected:e}"
        );
    }
    assert!(
        (record.delta_h + record.physical_heat - record.external_work - record.boundary_work).abs()
            < 3e-13
    );
    let mut half = PreparedScatteringZoh::new_owned(body(), 1, z, dt / 2.0, budget(5)).unwrap();
    let a = half.step(&initial, &inputs).unwrap();
    let b = half.step(&a.x, &inputs).unwrap();
    for i in 0..5 {
        assert!((step.x[i] - b.x[i]).abs() < 3e-12);
    }
    assert!(
        (record.physical_heat - a.account.physical_heat - b.account.physical_heat).abs() < 5e-13
    );
    assert!(
        (record.outgoing_energy - a.account.outgoing_energy - b.account.outgoing_energy).abs()
            < 5e-13
    );
}

fn matched_port() -> PortHamiltonian {
    PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.25],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap()
}

#[test]
fn nearly_matched_outgoing_energy_is_positive_without_subtracting_incident_energy() {
    let mut runtime =
        PreparedScatteringZoh::new_owned(matched_port(), 0, 4.0, 2.0, budget(1)).unwrap();
    let matched = runtime.step(&[1.0], &[1.0]).unwrap();
    assert_eq!(matched.account.outgoing_energy, 0.0);
    assert!((matched.account.physical_heat - 0.5).abs() < 2e-13);
    let x = 1.0 + 1e-10;
    let step = runtime.step(&[x], &[1.0]).unwrap();
    let expected = (x - 1.0).powi(2) * -(-2.0_f64).exp_m1() / 4.0;
    assert!(step.account.outgoing_energy > 0.0);
    assert!((step.account.outgoing_energy / expected - 1.0).abs() < 2e-10);
}

#[test]
fn nonlinear_storage_invalid_ports_and_failed_steps_do_not_publish_outputs() {
    assert!(
        PreparedScatteringZoh::new_owned(
            duffing_oscillator(1.0, 1.0, 2.0, 0.1).unwrap(),
            0,
            1.0,
            0.01,
            budget(2)
        )
        .is_err()
    );
    for (port, z, dt, modes) in [
        (3, 1.0, 0.01, 5),
        (1, 0.0, 0.01, 5),
        (1, f64::NAN, 0.01, 5),
        (1, 1.0, 0.0, 5),
        (1, 1.0, 0.01, 4),
    ] {
        assert!(PreparedScatteringZoh::new_owned(body(), port, z, dt, budget(modes)).is_err());
    }
    let mut runtime = PreparedScatteringZoh::new_owned(body(), 1, 2.5, 0.001, budget(5)).unwrap();
    let mut next = vec![123.0; 5];
    let mut output = vec![456.0; 3];
    assert!(
        runtime
            .step_into(&[0.0; 5], &[0.0, f64::NAN, 0.0], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, vec![123.0; 5]);
    assert_eq!(output, vec![456.0; 3]);
    assert!(
        runtime
            .step_into(&[0.0; 4], &[0.0; 3], &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, vec![123.0; 5]);
    let inactive = PortHamiltonian::new(
        1,
        2,
        vec![0.0],
        vec![0.25],
        vec![1.0, 0.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let mut inactive = PreparedScatteringZoh::new_owned(inactive, 1, 2.0, 0.01, budget(1)).unwrap();
    let step = inactive.step(&[0.1], &[0.3, 4.0]).unwrap();
    assert_eq!(step.y[1], 0.0);
    assert_eq!(step.account.boundary_work, 0.0);
    assert_eq!(step.account.incident_energy, step.account.outgoing_energy);
}

#[test]
fn prepared_wave_steps_and_current_window_prediction_reuse_all_scratch() {
    let mut runtime = PreparedScatteringZoh::new_owned(body(), 1, 2.5, 0.001, budget(5)).unwrap();
    let mut state = vec![0.02, 0.1, -0.03, 0.04, 0.03];
    let mut next = vec![0.0; 5];
    let mut output = vec![0.0; 3];
    let inputs = [0.3, 0.07, -0.2];
    TRACK.with(|cell| cell.set((true, 0)));
    for _ in 0..1000 {
        let (free, beta) = runtime.average_relation(&state, &inputs).unwrap();
        let record = runtime
            .step_into(&state, &inputs, &mut next, &mut output)
            .unwrap();
        assert!((output[1] - free - beta * inputs[1]).abs() < 3e-13);
        assert!(record.physical_heat >= 0.0 && record.outgoing_energy >= 0.0);
        core::mem::swap(&mut state, &mut next);
    }
    let allocations = TRACK.with(|cell| {
        let (_, count) = cell.replace((false, 0));
        count
    });
    assert_eq!(allocations, 0);
}
