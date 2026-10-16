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
