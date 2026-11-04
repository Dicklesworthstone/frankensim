//! Sponge (absorbing layer) conformance: MEASURED acoustic
//! reflection coefficients with a bounce-back-wall control (the
//! disabled-sponge mutation), steady-state non-interference, and the
//! constructor's typed panics.
//!
//! Method: D2Q9 is weakly compressible — a small density bump splits
//! into two acoustic pulses traveling at `c_s = 1/sqrt(3)` lattice
//! units per step. A probe between the source and the layer records
//! `max |rho - rho0|` in an INCIDENT time window (outbound pulse
//! passing the probe) and a REFLECTED window (whatever the layer
//! sends back). The ratio is the reflection coefficient R at that
//! pulse's spectral content; the same geometry with a bounce-back
//! wall in place of the layer measures the control. Both windows sit
//! at the SAME probe, so bulk viscous attenuation biases only the
//! extra sponge-to-probe leg (disclosed; the control is measured
//! identically, so the sponge-vs-wall comparison is fair).

use fs_lbm::core2::{Cell, Grid};
use fs_lbm::sponge::{Sponge2, SpongeSide};

const NX: usize = 1600;
const NY: usize = 3;
const TAU: f64 = 0.55;
const RHO0: f64 = 1.0;
const PULSE_X: usize = 400;
const PROBE_X: usize = 700;
/// Layer inner edge (width 200 reaching to the domain edge at 1600).
const SPONGE_W: usize = 200;

fn pulse_grid(pulse_half_width: f64) -> Grid {
    let mut grid = Grid::uniform(NX, NY, TAU);
    // Gaussian density bump, amplitude 1e-4 (linear acoustics).
    for x in 0..NX {
        #[allow(clippy::cast_precision_loss)]
        let dx = (x as f64 - PULSE_X as f64) / pulse_half_width;
        let rho = RHO0 + 1.0e-4 * (-dx * dx).exp();
        for y in 0..NY {
            let i = grid.idx(x, y);
            grid.f[i] = fs_lbm::equilibrium(rho, 0.0, 0.0);
        }
    }
    grid
}

fn probe_amplitude(grid: &Grid) -> f64 {
    let i = grid.idx(PROBE_X, 1);
    (grid.moments(i).rho - RHO0).abs()
}

/// Run `steps` steps; return (max amplitude in the incident window,
/// max amplitude in the reflected window). Window bounds are derived
/// from c_s = 1/sqrt(3): the rightward pulse passes the probe around
/// t = 300 * sqrt(3) = 520 and any reflection off the layer
/// (inner edge at x = 1400) returns to the probe around
/// t = (300 + 2 * 700) * sqrt(3) = 2944.
fn measure(mut grid: Grid, apply: impl Fn(&mut Grid)) -> (f64, f64) {
    let mut scratch = Vec::new();
    let mut incident = 0.0f64;
    let mut reflected = 0.0f64;
    for t in 0..3400 {
        grid.step(&mut scratch);
        apply(&mut grid);
        let a = probe_amplitude(&grid);
        if (300..900).contains(&t) {
            incident = incident.max(a);
        }
        if (2500..3400).contains(&t) {
            reflected = reflected.max(a);
        }
    }
    (incident, reflected)
}

/// The layer's measured reflection coefficient stays under the
/// authored ceiling for two pulse widths (different spectral
/// content), and the bounce-back-wall CONTROL measured identically is
/// orders of magnitude worse — the disabled-sponge mutation cannot
/// pass this test.
#[test]
fn sponge_reflection_coefficient_under_ceiling_with_wall_control() {
    let sponge = Sponge2::new(SpongeSide::RightX, SPONGE_W, 0.3, RHO0, [0.0, 0.0]);
    let mut results = Vec::new();
    for half_width in [12.0, 30.0] {
        let (incident, reflected) = measure(pulse_grid(half_width), |g| sponge.apply(g));
        assert!(
            incident > 1.0e-5,
            "incident pulse must actually pass the probe: {incident:e}"
        );
        let r = reflected / incident;
        // Measured floor: 1.1e-4 (half-width 12) / 2.3e-4 (30);
        // ceiling authored ~20x above the worse one (the
        // measured-tolerance doctrine: re-measure in the same commit
        // if the fixture is re-dimensioned).
        assert!(
            r < 5.0e-3,
            "sponge reflection coefficient {r:.4} over ceiling (half-width {half_width})"
        );
        results.push((half_width, r));
    }
    // CONTROL: bounce-back wall where the layer's inner edge sits —
    // the same measurement reads a large reflection (the
    // disabled-sponge mutation is caught by orders of magnitude).
    let mut walled = pulse_grid(12.0);
    for y in 0..NY {
        let i = walled.idx(NX - SPONGE_W, y);
        walled.flags[i] = Cell::Wall;
    }
    let (wi, wr) = measure(walled, |_| {});
    let r_wall = wr / wi;
    assert!(
        r_wall > 0.5,
        "wall control must reflect strongly: {r_wall:.4}"
    );
    assert!(
        results.iter().all(|&(_, r)| r_wall > 20.0 * r),
        "sponge must beat the wall control by >20x: wall {r_wall:.4} vs {results:?}"
    );
    println!(
        "{{\"suite\":\"fs-lbm\",\"case\":\"sponge-reflection\",\"r_pulse12\":{:.6},\"r_pulse30\":{:.6},\"r_wall_control\":{:.4},\"verdict\":\"pass\"}}",
        results[0].1, results[1].1, r_wall
    );
}

/// A sponge whose target equals the ambient state is a NO-OP on that
/// state: uniform rest fluid passes through unchanged to fp accuracy
/// (the layer must not inject anything into a matched far field).
#[test]
fn sponge_is_noop_on_matched_state() {
    let mut grid = Grid::uniform(600, 3, TAU);
    let sponge = Sponge2::new(SpongeSide::RightX, 100, 0.5, 1.0, [0.0, 0.0]);
    let before = grid.f.clone();
    let mut scratch = Vec::new();
    for _ in 0..50 {
        grid.step(&mut scratch);
        sponge.apply(&mut grid);
    }
    let worst = grid
        .f
        .iter()
        .zip(&before)
        .flat_map(|(a, b)| a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()))
        .fold(0.0f64, f64::max);
    assert!(worst < 1.0e-14, "matched sponge must be a no-op: {worst:e}");
}

/// Left-side layer absorbs a leftward pulse the same way (side
/// symmetry of the ramp indexing).
#[test]
fn sponge_left_side_absorbs() {
    let sponge = Sponge2::new(SpongeSide::LeftX, SPONGE_W, 0.3, RHO0, [0.0, 0.0]);
    // Mirror geometry: pulse at NX - PULSE_X, probe at NX - PROBE_X.
    let mut grid = Grid::uniform(NX, NY, TAU);
    for x in 0..NX {
        #[allow(clippy::cast_precision_loss)]
        let dx = (x as f64 - (NX - PULSE_X) as f64) / 12.0;
        let rho = RHO0 + 1.0e-4 * (-dx * dx).exp();
        for y in 0..NY {
            let i = grid.idx(x, y);
            grid.f[i] = fs_lbm::equilibrium(rho, 0.0, 0.0);
        }
    }
    let mut scratch = Vec::new();
    let mut incident = 0.0f64;
    let mut reflected = 0.0f64;
    let probe = grid.idx(NX - PROBE_X, 1);
    for t in 0..3400 {
        grid.step(&mut scratch);
        sponge.apply(&mut grid);
        let a = (grid.moments(probe).rho - RHO0).abs();
        if (300..900).contains(&t) {
            incident = incident.max(a);
        }
        if (2500..3400).contains(&t) {
            reflected = reflected.max(a);
        }
    }
    let r = reflected / incident;
    assert!(r < 5.0e-3, "left-side sponge reflection {r:.4}");
}

/// Constructor refusals (crate boundary convention: checked asserts,
/// matching `VelocityPressureX2`).
#[test]
#[should_panic(expected = "sigma_max must lie in (0, 1]")]
fn sponge_refuses_bad_sigma() {
    let _ = Sponge2::new(SpongeSide::RightX, 100, 1.5, 1.0, [0.0, 0.0]);
}

#[test]
#[should_panic(expected = "target density must be positive")]
fn sponge_refuses_bad_density() {
    let _ = Sponge2::new(SpongeSide::RightX, 100, 0.3, -1.0, [0.0, 0.0]);
}

#[test]
#[should_panic(expected = "width must be positive")]
fn sponge_refuses_zero_width() {
    let _ = Sponge2::new(SpongeSide::RightX, 0, 0.3, 1.0, [0.0, 0.0]);
}

#[test]
#[should_panic(expected = "low-Mach boundary envelope")]
fn sponge_refuses_fast_target() {
    let _ = Sponge2::new(SpongeSide::RightX, 100, 0.3, 1.0, [0.2, 0.0]);
}

/// Determinism: bitwise-identical reruns.
#[test]
fn sponge_determinism_bitwise() {
    let run = || {
        let mut grid = pulse_grid(12.0);
        let sponge = Sponge2::new(SpongeSide::RightX, SPONGE_W, 0.3, RHO0, [0.0, 0.0]);
        let mut scratch = Vec::new();
        for _ in 0..400 {
            grid.step(&mut scratch);
            sponge.apply(&mut grid);
        }
        grid.f[grid.idx(PROBE_X, 1)]
    };
    let a = run();
    let b = run();
    for (x, y) in a.iter().zip(b.iter()) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
}

/// Independent physical moments, including the first two Guo source moments.
/// A fringe must blend these moments rather than the shifted stored ones.
fn physical_moments(grid: &Grid, cell: usize) -> (f64, [f64; 2], [[f64; 2]; 2]) {
    let directions = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [-1.0, 0.0],
        [0.0, -1.0],
        [1.0, 1.0],
        [-1.0, 1.0],
        [-1.0, -1.0],
        [1.0, -1.0],
    ];
    let rho = grid.f[cell].iter().sum::<f64>();
    let force: [f64; 2] = std::array::from_fn(|a| rho * grid.g[a] + grid.fext[cell][a]);
    let momentum: [f64; 2] = std::array::from_fn(|a| {
        grid.f[cell]
            .iter()
            .zip(directions)
            .map(|(f, e)| f * e[a])
            .sum::<f64>()
            + 0.5 * force[a]
    });
    let velocity = momentum.map(|p| p / rho);
    let second = std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            grid.f[cell]
                .iter()
                .zip(directions)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
        })
    });
    (rho, momentum, second)
}

#[test]
fn sponge_g0_preserves_matched_physical_state_under_forcing() {
    for side in [SpongeSide::LeftX, SpongeSide::RightX] {
        for velocity in [[0.0, 0.0], [0.035, -0.02]] {
            let mut grid = Grid::uniform(5, 3, 0.8);
            grid.f
                .fill(fs_lbm::equilibrium(1.4, velocity[0], velocity[1]));
            let external: Vec<_> = (0..grid.f.len())
                .map(|i| [0.001 * (i % 5) as f64, -0.002 * (i / 5) as f64])
                .collect();
            grid.rebind_forcing([0.003, -0.004], &external, &mut Vec::new())
                .unwrap();
            let before = grid.f.clone();
            let sponge = Sponge2::new(side, 5, 1.0, 1.4, velocity);
            for _ in 0..20 {
                sponge.apply(&mut grid);
            }
            for (actual, expected) in grid.f.iter().flatten().zip(before.iter().flatten()) {
                assert!((actual - expected).abs() < 2.0e-14);
            }
            for cell in 0..grid.f.len() {
                for (actual, expected) in grid.moments(cell).u.into_iter().zip(velocity) {
                    assert!((actual - expected).abs() < 2.0e-14);
                }
            }
            assert_eq!(grid.g, [0.003, -0.004]);
            assert_eq!(grid.fext, external);
        }
    }
}

#[test]
fn sponge_g0_blends_physical_mass_momentum_and_full_stress() {
    let profile = [
        (0.9, [0.03, -0.01]),
        (1.3, [-0.02, 0.025]),
        (1.1, [0.01, 0.0]),
    ];
    for side in [SpongeSide::LeftX, SpongeSide::RightX] {
        for sigma_max in [0.4, 1.0] {
            let mut grid = Grid::uniform(6, profile.len(), 0.8);
            for cell in 0..grid.f.len() {
                let rho = 1.2 + 0.01 * cell as f64;
                grid.f[cell] = fs_lbm::equilibrium(rho, -0.01, 0.015);
                // Non-equilibrium normal and shear stress must relax with the
                // authored strength, not disappear in an equilibrium reset.
                grid.f[cell][1] += 0.002;
                grid.f[cell][3] += 0.002;
                grid.f[cell][5] += 0.001;
                grid.f[cell][7] += 0.001;
                grid.f[cell][0] -= 0.006;
            }
            for (x, flag) in [(1, Cell::Wall), (2, Cell::Interface), (3, Cell::Gas)] {
                let cell = grid.idx(x, 1);
                grid.flags[cell] = flag;
            }
            let external: Vec<_> = (0..grid.f.len())
                .map(|i| [-0.0004 * (i % 6) as f64, 0.0007 * (i / 6) as f64])
                .collect();
            grid.rebind_forcing([0.003, -0.004], &external, &mut Vec::new())
                .unwrap();
            let before = grid.clone();
            Sponge2::with_profile(side, 4, sigma_max, &profile).apply(&mut grid);
            for (y, &(rho_t, u_t)) in profile.iter().enumerate() {
                for x in 0..grid.nx {
                    let cell = grid.idx(x, y);
                    let depth = match side {
                        SpongeSide::LeftX => 3usize.checked_sub(x),
                        SpongeSide::RightX => x.checked_sub(grid.nx - 4),
                    };
                    if depth.is_none() || grid.flags[cell] != Cell::Fluid {
                        assert_eq!(grid.f[cell], before.f[cell]);
                        continue;
                    }
                    let t = (depth.unwrap() as f64 + 1.0) / 4.0;
                    let sigma = sigma_max * t * t;
                    let (rho, momentum, second) = physical_moments(&before, cell);
                    let (new_rho, new_momentum, new_second) = physical_moments(&grid, cell);
                    let expected_rho = (1.0 - sigma) * rho + sigma * rho_t;
                    assert!((new_rho - expected_rho).abs() < 2.0e-14);
                    for a in 0..2 {
                        let expected_p = (1.0 - sigma) * momentum[a] + sigma * rho_t * u_t[a];
                        assert!((new_momentum[a] - expected_p).abs() < 2.0e-14);
                        for b in 0..2 {
                            let target_second =
                                rho_t * (u_t[a] * u_t[b] + if a == b { 1.0 / 3.0 } else { 0.0 });
                            let expected_second =
                                (1.0 - sigma) * second[a][b] + sigma * target_second;
                            assert!((new_second[a][b] - expected_second).abs() < 2.0e-14);
                        }
                    }
                }
            }
            assert_eq!(grid.g, before.g);
            assert_eq!(grid.fext, before.fext);
            assert_eq!(grid.tau, before.tau);
            assert_eq!(grid.flags, before.flags);
        }
    }
}

#[test]
fn sponge_g3_unforced_populations_retain_original_blending_bits() {
    let mut grid = Grid::uniform(4, 1, 0.8);
    let targets = fs_lbm::equilibrium(0.9, 0.02, -0.01);
    let expected: Vec<_> = grid
        .f
        .iter()
        .enumerate()
        .map(|(x, f)| {
            let t = (x as f64 + 1.0) / 4.0;
            let sigma = 0.7 * t * t;
            std::array::from_fn(|q| f[q] + sigma * (targets[q] - f[q]))
        })
        .collect();
    Sponge2::new(SpongeSide::RightX, 4, 0.7, 0.9, [0.02, -0.01]).apply(&mut grid);
    assert_eq!(grid.f, expected);
}
