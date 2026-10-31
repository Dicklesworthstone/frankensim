//! A grid-scale wave at the walls of a body-forced channel (pexfp.138).
//!
//! Four columns, ten fluid rows, every face periodic, reduced cumulant. The flow settles
//! on the scheme's own parabola, then a wave at wavenumber pi in x is added on the two
//! rows next to each wall. With no open face and every column the same, that is one
//! wavenumber, so the wave has one rate.
//!
//! The bound on the rate is zero. This flow is far below the linear limit for plane
//! Poiseuille flow (Orszag, J. Fluid Mech. 50 (1971): 5772 on centreline speed and
//! half-width, about 250 here at the higher Re below), so any growth is the scheme's.
//! No threshold comes from a measured rate.
//!
//! This test sees k = pi only. A fix can pass it while nearby wavenumbers still grow, so
//! it can't stand alone as the green-after test for pexfp.138.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::{BoundaryGrid3, BoundarySpec3, FaceBoundary3, equilibrium3};

const NX: usize = 4;
const NY: usize = 12;
const NZ: usize = 4;
const H: f64 = 10.0;
const AXIS: f64 = 6.0;
const U_LAT: f64 = 0.05;
const SEED: f64 = 1e-6;
/// The wave is read on every `DUMP`-th step and the step after it.
const DUMP: usize = 250;
const LAST_DUMP: usize = 2_500;
const FIT_FROM: usize = 750;
/// Three e-folds of h^2 / (pi^2 nu), rounded down to even, so a start-up loss has gone.
const SETTLE_RE_251: usize = 15_198;
const SETTLE_RE_502: usize = 30_396;

/// Adds `seed` to u on the rows next to the walls, (-1)^x, + on rows 1 and 10, - on rows
/// 2 and 9. Each population takes the change in equilibrium, so the settled stress stays.
fn add_wave(grid: &mut BoundaryGrid3, seed: f64) {
    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..NX {
                let row = match y.min(NY - 1 - y) {
                    1 => 1.0,
                    2 => -1.0,
                    _ => 0.0,
                };
                let sign = if x % 2 == 0 { 1.0 } else { -1.0 };
                let du = sign * row * seed;
                if du != 0.0 && !grid.is_solid(x, y, z) {
                    let rho = grid.density(x, y, z);
                    let u = grid.velocity(x, y, z);
                    let from = equilibrium3(rho, u);
                    let to = equilibrium3(rho, [u[0] + du, u[1], u[2]]);
                    let mut populations = grid.populations(x, y, z);
                    for (q, value) in populations.iter_mut().enumerate() {
                        *value += to[q] - from[q];
                    }
                    grid.set_populations(x, y, z, &populations);
                }
            }
        }
    }
}

/// The wave's size at each dump t: the mean over rows 1 and 10 of |1/4 sum (-1)^(x+s) u|
/// over steps t and t + 1 and one pair of columns, for the larger pair. `re` is the key
/// that sets nu; Re read as 1.5 bulk h / nu is about 1.0047 times it.
fn wave_sizes(re: f64, settle: usize, seed: f64) -> Vec<(usize, f64)> {
    let nu = U_LAT * H / re;
    let rate = 1.0 / (3.0 * nu + 0.5);
    // The scheme's own wall position: H_eff^2 = H^2 + (16 L - 3) / 3.
    let lambda = (1.0 / rate - 0.5) * 0.5;
    let h2 = H * H + (16.0 * lambda - 3.0) / 3.0;
    let mean_y2 = (0..10)
        .map(|j| (f64::from(j) + 0.5 - H / 2.0).powi(2))
        .sum::<f64>()
        / H;
    let amplitude = 0.67 * U_LAT / (h2 / 4.0 - mean_y2);

    let periodic = FaceBoundary3::Periodic;
    let mut grid = BoundaryGrid3::with_collision_model(
        NX,
        NY,
        NZ,
        CollisionModel3::ReducedCumulant {
            second_order_rate: rate,
            third_order_rate: 1.0,
            fourth_order_rate: 1.0,
        },
        [2.0 * nu * amplitude, 0.0, 0.0],
        BoundarySpec3::new([periodic; 6]),
    );
    grid.voxelize_sdf(|c| H / 2.0 - (c[1] - AXIS).abs());
    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..NX {
                if !grid.is_solid(x, y, z) {
                    let yc = y as f64 + 0.5 - AXIS;
                    let u = amplitude * (h2 / 4.0 - yc * yc);
                    grid.set_populations(x, y, z, &equilibrium3(1.0, [u, 0.0, 0.0]));
                }
            }
        }
    }

    let gate = CancelGate::new();
    let pool = TilePool::for_host(1, 0);
    pool.with_parked_crew_local(|parked| {
        for step in 1..=settle {
            grid.step_pooled(parked, &gate)
                .unwrap_or_else(|err| panic!("settle step {step}: {err}"));
        }
        add_wave(&mut grid, seed);

        let mut sizes = Vec::new();
        let mut signed = [[0.0; 2]; NX / 2];
        for step in 1..=LAST_DUMP + 1 {
            grid.step_pooled(parked, &gate)
                .unwrap_or_else(|err| panic!("step {step}: {err}"));
            let phase = step % DUMP;
            if step > 1 && phase <= 1 {
                for (pair, rows) in signed.iter_mut().enumerate() {
                    for (row, y) in [1, NY - 2].into_iter().enumerate() {
                        for x in [2 * pair, 2 * pair + 1] {
                            let sign = if (x + step) % 2 == 0 { 1.0 } else { -1.0 };
                            rows[row] += sign * grid.velocity(x, y, 0)[0];
                        }
                    }
                }
            }
            if step > 1 && phase == 1 {
                let size = signed
                    .iter()
                    .map(|rows| (rows[0].abs() / 4.0 + rows[1].abs() / 4.0) / 2.0)
                    .fold(0.0, f64::max);
                sizes.push((step - 1, size));
                signed = [[0.0; 2]; NX / 2];
            }
        }
        sizes
    })
}

/// Least-squares slope of ln(size) against step, dumps `FIT_FROM` to `LAST_DUMP`.
fn fitted_rate(sizes: &[(usize, f64)]) -> f64 {
    let points: Vec<(f64, f64)> = sizes
        .iter()
        .filter(|(step, _)| *step >= FIT_FROM)
        .map(|&(step, size)| (step as f64, size.ln()))
        .collect();
    let count = points.len() as f64;
    let mean_step = points.iter().map(|p| p.0).sum::<f64>() / count;
    let mean_log = points.iter().map(|p| p.1).sum::<f64>() / count;
    let rise: f64 = points
        .iter()
        .map(|p| (p.0 - mean_step) * (p.1 - mean_log))
        .sum();
    let run: f64 = points.iter().map(|p| (p.0 - mean_step).powi(2)).sum();
    rise / run
}

/// The seeded rate, after checking the measure can't pass by seeing nothing: seed off it
/// reads exactly 0, and seeded its first dump reads the seed to within a factor of 10.
fn seeded_rate(re: f64, settle: usize) -> f64 {
    let quiet = wave_sizes(re, settle, 0.0);
    assert!(
        quiet.iter().all(|&(_, size)| size == 0.0),
        "seed-off wave at re key {re}: {quiet:?}"
    );
    let sizes = wave_sizes(re, settle, SEED);
    let first = sizes[0].1;
    assert!(
        first >= SEED / 10.0,
        "first dump at re key {re} reads {first:e} for a seed of {SEED:e}"
    );
    let rate = fitted_rate(&sizes);
    println!("wave rate {rate:+.4e} per step at re key {re}");
    rate
}

/// This reads a mixture, not the wall wave alone. Below the crossing near Re 319 two real
/// modes trade character: at Re 251 the slowest mode is mostly a viscous mode of the
/// resting fluid, about -8.3e-4 per step, and the wall wave is second at -2.69e-3. The
/// fitted rate is a blend of the two. Both decay, so the assertion still holds.
#[test]
fn wall_wave_decays_at_re_251() {
    let rate = seeded_rate(250.0, SETTLE_RE_251);
    assert!(rate < 0.0, "wave rate {rate:e} per step at Re 251");
}

#[test]
#[ignore = "pexfp.138: fails, the wave grows at +2.94e-3 per step at Re 502, an e-fold of 340 steps"]
fn wall_wave_must_not_grow_at_re_502() {
    let rate = seeded_rate(500.0, SETTLE_RE_502);
    assert!(rate <= 0.0, "wave rate {rate:e} per step at Re 502");
}
