//! A density wave at rest that varies along all three axes grows on D3Q19's reduced
//! cumulant and decays on D3Q27's cumulant (pexfp.138).
//!
//! At rest the relaxed rows are Hermite polynomials. On D3Q19 the fourth-order row H220
//! and the second-order row H002 have weighted inner product -1/27, the eight corner
//! velocities D3Q19 lacks, and they relax at different rates. D3Q27 has no such overlap,
//! so its linear step can't grow. The one-axis twin can't excite the overlap's mode on
//! either lattice, so it shows the test can tell the two apart.
//!
//! Periodic 4 x 4 x 4, no wall, no force, nu = 1e-3, every other rate 1. The test reads
//! only whether the departure from rest is larger or smaller after a fixed number of
//! steps; no rate enters an assertion.

use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, CollisionModel3, Q3, W3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, Grid27, Rates27, W27};
use fs_lbm::lattice3::GridBudget3;

const N: usize = 4;
const STEPS: usize = 300;
const AMPLITUDE: f64 = 1e-8;
/// Shear (and, on D3Q27, bulk) rate for nu = 1e-3.
const SECOND_ORDER_RATE: f64 = 1.0 / (3.0 * 1e-3 + 0.5);
/// Wavenumber 2 pi m / 4 on each axis: (pi, pi, pi/2).
const THREE_AXES: [usize; 3] = [2, 2, 1];
/// The same x wavenumber with y and z uniform: (pi, 0, 0).
const ONE_AXIS: [usize; 3] = [2, 0, 0];

/// A cos(k.x) for k = 2 pi m / 4, from quarter turns, so the start is exact.
fn wave(p: [usize; 3], m: [usize; 3]) -> f64 {
    let quarter = (m[0] * p[0] + m[1] * p[1] + m[2] * p[2]) % 4;
    AMPLITUDE * [1.0, 0.0, -1.0, 0.0][quarter]
}

/// Root mean square over cells of |f - w|.
fn departure<const Q: usize>(cells: &[[f64; Q]], weights: &[f64; Q]) -> f64 {
    let sum: f64 = cells
        .iter()
        .flat_map(|f| f.iter().zip(weights))
        .map(|(value, weight)| (value - weight) * (value - weight))
        .sum();
    (sum / cells.len() as f64).sqrt()
}

fn d3q19_cells(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut cells = Vec::with_capacity(N * N * N);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                cells.push(grid.populations(x, y, z));
            }
        }
    }
    cells
}

/// Departure after `STEPS` against the start, on D3Q19 `ReducedCumulant`.
fn d3q19_ratio(m: [usize; 3]) -> f64 {
    let mut grid = BoundaryGrid3::with_collision_model(
        N,
        N,
        N,
        CollisionModel3::ReducedCumulant {
            second_order_rate: SECOND_ORDER_RATE,
            third_order_rate: 1.0,
            fourth_order_rate: 1.0,
        },
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let density = 1.0 + wave([x, y, z], m);
                let f: [f64; Q3] = core::array::from_fn(|q| W3[q] * density);
                grid.set_populations(x, y, z, &f);
            }
        }
    }
    let start = departure(&d3q19_cells(&grid), &W3);
    for _ in 0..STEPS {
        grid.step();
    }
    departure(&d3q19_cells(&grid), &W3) / start
}

fn with_cx<R>(body: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 138,
                kernel_id: 138,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        body(&cx)
    })
}

/// Departure after `STEPS` against the start, on D3Q27 `Cumulant`, omega [r, r, 1 x 8].
fn d3q27_ratio(m: [usize; 3]) -> f64 {
    let mut omega = [1.0; 10];
    omega[0] = SECOND_ORDER_RATE;
    omega[1] = SECOND_ORDER_RATE;
    let mut grid = Grid27::try_new(
        [N; 3],
        CollisionConfig27 {
            model: CollisionModel27::Cumulant {
                rates: Rates27 { omega },
            },
            maximum_mach: 0.2,
        },
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: N * N * N,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let density = 1.0 + wave(p, m);
            core::array::from_fn(|q| W27[q] * density)
        },
    )
    .expect("periodic D3Q27 rest grid");
    let start = departure(grid.populations(), &W27);
    with_cx(|cx| {
        for step in 1..=STEPS {
            grid.step_serial(cx)
                .unwrap_or_else(|err| panic!("D3Q27 step {step}: {err:?}"));
        }
    });
    departure(grid.populations(), &W27) / start
}

#[test]
fn three_axis_wave_grows_on_reduced_cumulant() {
    let ratio = d3q19_ratio(THREE_AXES);
    println!("D3Q19 reduced cumulant, three axes: {ratio:e} of the start after {STEPS} steps");
    assert!(ratio > 1.0, "three-axis wave did not grow: {ratio:e}");
}

#[test]
fn three_axis_wave_decays_on_d3q27_cumulant() {
    let ratio = d3q27_ratio(THREE_AXES);
    println!("D3Q27 cumulant, three axes: {ratio:e} of the start after {STEPS} steps");
    assert!(ratio < 1.0, "three-axis wave did not decay: {ratio:e}");
}

#[test]
fn one_axis_twin_decays_on_both_lattices() {
    let d3q19 = d3q19_ratio(ONE_AXIS);
    let d3q27 = d3q27_ratio(ONE_AXIS);
    println!("one axis: D3Q19 {d3q19:e}, D3Q27 {d3q27:e} of the start after {STEPS} steps");
    assert!(d3q19 < 1.0, "one-axis wave grew on D3Q19: {d3q19:e}");
    assert!(d3q27 < 1.0, "one-axis wave grew on D3Q27: {d3q27:e}");
}
