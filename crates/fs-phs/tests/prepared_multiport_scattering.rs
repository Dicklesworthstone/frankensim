//! Joint reflection retains cross-port dynamics, physical memory and wave power.
// Test-only allocation observation follows prepared_exact_owned.rs.
#![allow(unsafe_code)]
use fs_phs::linear_exact::{
    LinearSplitBudget, PreparedMultiportScatteringZoh, PreparedScatteringZoh,
    ScatteringPortRecord,
};
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


fn run(
    runtime: &mut PreparedMultiportScatteringZoh,
    state: &[f64],
    inputs: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<ScatteringPortRecord>, fs_phs::linear_exact::MultiportScatteringStepRecord) {
    let mut next = vec![0.0; runtime.system().state_dim()];
    let mut mean = vec![0.0; runtime.system().port_dim()];
    let mut ports = vec![ScatteringPortRecord::default(); runtime.ports().len()];
    let account = runtime.step_into(state, inputs, &mut next, &mut mean, &mut ports).unwrap();
    (next, mean, ports, account)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:e} vs {expected:e}");
}

#[test]
fn jointly_closed_ports_match_independent_physical_trajectory_and_each_wave_integral() {
    let initial = [0.02, 0.1, -0.03, 0.04, 0.03];
    let input = [0.3, 0.07, -0.2];
    let (z0, z2, dt) = (1.7, 3.2, 0.037);
    let original = body();
    let structure = original.structure();
    let source = (structure.0.to_vec(), structure.1.to_vec(), structure.2.to_vec());
    let mut runtime = PreparedMultiportScatteringZoh::new_owned(
        original, &[(0, z0), (2, z2)], dt, budget(5),
    ).unwrap();
    assert_eq!(runtime.system().structure(), (source.0.as_slice(), source.1.as_slice(), source.2.as_slice()));
    let (mut free, mut direct) = ([0.0; 2], [0.0; 4]);
    runtime.average_relation_into(&initial, &input, &mut free, &mut direct).unwrap();
    assert!(direct[1].abs() > 1e-5 && direct[2].abs() > 1e-5);
    let (next, mean, ports, account) = run(&mut runtime, &initial, &input);
    for row in 0..2 {
        close(mean[[0, 2][row]], free[row] + direct[2 * row] * input[0]
            + direct[2 * row + 1] * input[2], 3e-13);
    }
    let reference = rk4(
        [initial[0], initial[1], initial[2], initial[3], initial[4],
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        dt, 2048, |x| {
            let arm = 2.0 * x[0] - x[2] - x[4];
            let y0 = x[1] - 0.3 * x[3];
            let y1 = 0.7 * x[1] - 0.4 * x[3];
            let y2 = -0.2 * x[1] + 1.1 * x[3];
            let u0 = (2.0 * input[0] - y0) / z0;
            let u2 = (2.0 * input[2] - y2) / z2;
            [
                x[1],
                -4.0 * x[0] - 0.4 * x[1] - 0.5 * x[3] - 2.0 * arm
                    + u0 + 0.7 * input[1] - 0.2 * u2,
                x[3],
                -9.0 * x[2] - 0.5 * x[1] - 0.8 * x[3] + arm
                    - 0.3 * u0 - 0.4 * input[1] + 1.1 * u2,
                arm / 0.15,
                0.4 * x[1] * x[1] + x[1] * x[3] + 0.8 * x[3] * x[3]
                    + arm * arm / 0.15,
                input[1] * y1,
                u0 * y0, u2 * y2,
                input[0] * input[0] / z0, input[2] * input[2] / z2,
                (y0 - input[0]).powi(2) / z0, (y2 - input[2]).powi(2) / z2,
                y0, y2, y1,
            ]
        },
    );
    for i in 0..5 { close(next[i], reference[i], 3e-12); }
    close(account.physical_heat, reference[5], 5e-13);
    close(account.external_work, reference[6], 5e-13);
    for i in 0..2 {
        close(ports[i].boundary_work, reference[7 + i], 5e-13);
        close(ports[i].incident_energy, reference[9 + i], 5e-13);
        close(ports[i].outgoing_energy, reference[11 + i], 5e-13);
        close(ports[i].outgoing_wave, mean[[0, 2][i]] - input[[0, 2][i]], 0.0);
        close(mean[[0, 2][i]] * dt, reference[13 + i], 5e-13);
        // Time-varying reflected pressure carries variance as well as its mean.
        assert!(ports[i].outgoing_energy > dt * ports[i].outgoing_wave.powi(2) / [z0, z2][i]);
    }
    close(mean[1] * dt, reference[15], 5e-13);
    close(account.boundary_work, ports[0].boundary_work + ports[1].boundary_work, 0.0);
    close(account.outgoing_energy, ports[0].outgoing_energy + ports[1].outgoing_energy, 0.0);
    close(account.delta_h + account.physical_heat, account.external_work + account.boundary_work, 5e-13);
}

#[test]
fn complete_interval_and_two_half_intervals_have_the_same_port_energies() {
    let state = [0.02, 0.1, -0.03, 0.04, 0.03];
    let input = [0.3, 0.07, -0.2];
    let selection = [(0, 1.7), (2, 3.2)];
    let mut full = PreparedMultiportScatteringZoh::new_owned(body(), &selection, 0.037, budget(5)).unwrap();
    let mut half = PreparedMultiportScatteringZoh::new_owned(body(), &selection, 0.0185, budget(5)).unwrap();
    let (x, y, ports, all) = run(&mut full, &state, &input);
    let (a, ya, pa, aa) = run(&mut half, &state, &input);
    let (b, yb, pb, ab) = run(&mut half, &a, &input);
    for i in 0..5 { close(x[i], b[i], 3e-12); }
    for i in 0..3 { close(y[i], 0.5 * (ya[i] + yb[i]), 3e-12); }
    close(all.physical_heat, aa.physical_heat + ab.physical_heat, 5e-13);
    close(all.external_work, aa.external_work + ab.external_work, 5e-13);
    for i in 0..2 {
        close(ports[i].boundary_work, pa[i].boundary_work + pb[i].boundary_work, 5e-13);
        close(ports[i].outgoing_energy, pa[i].outgoing_energy + pb[i].outgoing_energy, 5e-13);
    }
}

#[test]
fn selection_order_changes_only_the_returned_selected_port_layout() {
    let state = [0.02, 0.1, -0.03, 0.04, 0.03];
    let input = [0.3, 0.07, -0.2];
    let mut a = PreparedMultiportScatteringZoh::new_owned(body(), &[(0, 1.7), (2, 3.2)], 0.037, budget(5)).unwrap();
    let mut b = PreparedMultiportScatteringZoh::new_owned(body(), &[(2, 3.2), (0, 1.7)], 0.037, budget(5)).unwrap();
    let (xa, ya, pa, _) = run(&mut a, &state, &input);
    let (xb, yb, pb, _) = run(&mut b, &state, &input);
    for i in 0..5 { close(xa[i], xb[i], 3e-13); }
    for i in 0..3 { close(ya[i], yb[i], 3e-13); }
    for i in 0..2 {
        close(pa[i].outgoing_energy, pb[1 - i].outgoing_energy, 3e-13);
        close(pa[i].boundary_work, pb[1 - i].boundary_work, 3e-13);
    }
    let (mut fa, mut fb, mut da, mut db) = ([0.0; 2], [0.0; 2], [0.0; 4], [0.0; 4]);
    a.average_relation_into(&state, &input, &mut fa, &mut da).unwrap();
    b.average_relation_into(&state, &[99.0, input[1], -99.0], &mut fb, &mut db).unwrap();
    for i in 0..2 {
        close(fa[i], fb[1 - i], 3e-13);
        for j in 0..2 { close(da[2 * i + j], db[2 * (1 - i) + 1 - j], 3e-13); }
    }
}

#[test]
fn more_selected_ports_than_states_preserve_inactive_reflection_and_analytic_decay() {
    let system = PortHamiltonian::new(
        1, 3, vec![0.0], vec![0.25], vec![1.0, 2.0, 0.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    ).unwrap();
    let mut runtime = PreparedMultiportScatteringZoh::new_owned(
        system, &[(0, 4.0), (1, 2.0), (2, 3.0)], 0.2, budget(1),
    ).unwrap();
    let (next, mean, ports, account) = run(&mut runtime, &[0.4], &[0.8, -0.3, 4.0]);
    let equilibrium = -0.2 / 2.5;
    close(next[0], equilibrium + (0.4 - equilibrium) * (-0.5_f64).exp(), 3e-13);
    let average = equilibrium + (0.4 - equilibrium) * -(-0.5_f64).exp_m1() / 0.5;
    close(mean[0], average, 3e-13);
    close(mean[1], 2.0 * average, 3e-13);
    assert_eq!(mean[2], 0.0);
    assert_eq!(ports[2].outgoing_wave, -4.0);
    assert_eq!(ports[2].boundary_work, 0.0);
    assert_eq!(ports[2].outgoing_energy, ports[2].incident_energy);
    close(account.delta_h + account.physical_heat, account.boundary_work, 3e-13);
}

fn matched_pair() -> PortHamiltonian {
    PortHamiltonian::new(
        2, 2, vec![0.0; 4], vec![0.25, 0.0, 0.0, 0.0625],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
    ).unwrap()
}

#[test]
fn each_nearly_matched_channel_retains_its_tiny_positive_reflected_energy() {
    let mut runtime = PreparedMultiportScatteringZoh::new_owned(
        matched_pair(), &[(0, 4.0), (1, 16.0)], 2.0, budget(2),
    ).unwrap();
    let (_, _, matched, _) = run(&mut runtime, &[1.0, 2.0], &[1.0, 2.0]);
    assert_eq!(matched[0].outgoing_energy, 0.0);
    assert_eq!(matched[1].outgoing_energy, 0.0);
    let state: [f64; 2] = [1.0 + 1e-10, 2.0 - 2e-10];
    let (_, _, ports, _) = run(&mut runtime, &state, &[1.0, 2.0]);
    for (i, z, rate, a) in [(0, 4.0, 0.5, 1.0), (1, 16.0, 0.125, 2.0)] {
        let expected = (state[i] - a).powi(2) * -(-4.0_f64 * rate).exp_m1() / (2.0 * rate * z);
        assert!(ports[i].outgoing_energy > 0.0);
        close(ports[i].outgoing_energy / expected, 1.0, 5e-10);
    }
}

#[test]
fn one_port_adapter_uses_the_same_complete_flow_and_accounts() {
    let state = [0.02, 0.1, -0.03, 0.04, 0.03];
    let input = [0.3, 0.07, -0.2];
    let mut scalar = PreparedScatteringZoh::new_owned(body(), 1, 2.5, 0.037, budget(5)).unwrap();
    let mut multi = PreparedMultiportScatteringZoh::new_owned(body(), &[(1, 2.5)], 0.037, budget(5)).unwrap();
    let scalar = scalar.step(&state, &input).unwrap();
    let (x, y, ports, all) = run(&mut multi, &state, &input);
    assert_eq!(scalar.x, x);
    assert_eq!(scalar.y, y);
    assert_eq!(scalar.account.physical_heat, all.physical_heat);
    assert_eq!(scalar.account.outgoing_energy, ports[0].outgoing_energy);
    assert_eq!(scalar.account.outgoing_wave, ports[0].outgoing_wave);
}

#[test]
fn invalid_selections_and_refused_candidates_publish_no_partial_results() {
    for selection in [
        vec![], vec![(0, 1.0), (0, 2.0)], vec![(3, 1.0)],
        vec![(0, 0.0)], vec![(0, f64::NAN)], vec![(0, f64::INFINITY)],
    ] {
        assert!(PreparedMultiportScatteringZoh::new_owned(body(), &selection, 0.01, budget(5)).is_err());
    }
    assert!(PreparedMultiportScatteringZoh::new_owned(
        duffing_oscillator(1.0, 1.0, 2.0, 0.1).unwrap(), &[(0, 1.0)], 0.01, budget(2),
    ).is_err());
    let mut runtime = PreparedMultiportScatteringZoh::new_owned(body(), &[(0, 1.7), (2, 3.2)], 0.01, budget(5)).unwrap();
    let sentinel = ScatteringPortRecord {
        boundary_work: 71.0, incident_energy: 72.0, outgoing_energy: 73.0, outgoing_wave: 74.0,
    };
    let (mut next, mut mean, mut ports) = ([123.0; 5], [456.0; 3], [sentinel; 2]);
    for input in [[0.0, f64::NAN, 0.0], [f64::MAX, f64::MAX, f64::MAX]] {
        assert!(runtime.step_into(&[0.0; 5], &input, &mut next, &mut mean, &mut ports).is_err());
        assert_eq!(next, [123.0; 5]);
        assert_eq!(mean, [456.0; 3]);
        assert_eq!(ports, [sentinel; 2]);
    }
    assert!(runtime.step_into(&[0.0; 5], &[0.0; 3], &mut next, &mut mean, &mut ports[..1]).is_err());
    let (mut free, mut direct) = ([123.0; 2], [456.0; 4]);
    assert!(runtime.average_relation_into(&[0.0; 4], &[0.0; 3], &mut free, &mut direct).is_err());
    assert_eq!(free, [123.0; 2]);
    assert_eq!(direct, [456.0; 4]);
}

#[test]
fn repeated_joint_steps_and_prediction_allocate_no_owner_memory() {
    let mut runtime = PreparedMultiportScatteringZoh::new_owned(body(), &[(0, 1.7), (2, 3.2)], 0.001, budget(5)).unwrap();
    let mut state = [0.02, 0.1, -0.03, 0.04, 0.03];
    let (mut next, mut mean, mut ports) = ([0.0; 5], [0.0; 3], [ScatteringPortRecord::default(); 2]);
    let (mut free, mut direct) = ([0.0; 2], [0.0; 4]);
    let input = [0.3, 0.07, -0.2];
    TRACK.with(|cell| cell.set((true, 0)));
    for _ in 0..1000 {
        runtime.average_relation_into(&state, &input, &mut free, &mut direct).unwrap();
        let record = runtime.step_into(&state, &input, &mut next, &mut mean, &mut ports).unwrap();
        assert!(record.physical_heat >= 0.0 && ports.iter().all(|p| p.outgoing_energy >= 0.0));
        core::mem::swap(&mut state, &mut next);
    }
    let allocations = TRACK.with(|cell| {
        let (_, count) = cell.replace((false, 0));
        count
    });
    assert_eq!(allocations, 0);
}
