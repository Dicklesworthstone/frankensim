//! Stage 3: REFINE. The screen winner's section is rasterized into an
//! fs-lbm channel as a bounce-back obstacle and driven to a steady
//! state; aerodynamic forces come from fs-lbm's native stationary-wall
//! momentum exchange, averaged over the final 100 steps. The measured
//! boundary impulse includes the pressure and viscous loads represented by
//! the lattice without dropping unsteady storage or inlet reactions. The gate is
//! MODEL-FORM honest: the panel method and the low-Re lattice flow are
//! different physics, so agreement is gated on lift SIGN and
//! order-of-magnitude band, never point-matching; the honesty label
//! travels in the report. Cumulant LES on sparse VDB lattices with
//! DWR adaptivity is the recorded successor.

use crate::param::OrnithCandidate;
use fs_lbm::{Cell, Grid, Q, equilibrium};

const FLOW_STEPS: usize = 3200;
const FORCE_SAMPLES: usize = 100;

/// The refinement report — with its honesty label.
#[derive(Debug, Clone)]
pub struct RefineReport {
    /// Mean section lift over the final 100 steps (lattice units, +y).
    pub lift: f64,
    /// Mean section drag over the final 100 steps (lattice units, +x).
    pub drag: f64,
    /// Lift coefficient, normalized by freestream dynamic pressure and chord.
    pub cl: f64,
    /// Drag coefficient, normalized by freestream dynamic pressure and chord.
    pub cd: f64,
    /// Panel-method cl at the same trim (the screen's number).
    pub panel_cl: f64,
    /// Model-form honesty label (travels with every consumer).
    pub honesty: &'static str,
    /// Steady-state residual (velocity change per step at the end).
    pub steadiness: f64,
}

/// Rasterize the section into the grid as Wall cells.
fn rasterize(grid: &mut Grid, c: &OrnithCandidate, chord: f64, x0: f64, y0: f64) {
    let foil = c.section(64);
    let (nx, ny) = (grid.nx, grid.ny);
    let (ca, sa) = (fs_math::det::cos(c.alpha), fs_math::det::sin(c.alpha));
    for gy in 1..ny - 1 {
        for gx in 1..nx - 1 {
            // Cell center in section coordinates (rotate by −α).
            let px = (gx as f64 - x0) / chord;
            let py = (gy as f64 - y0) / chord;
            let sx = ca.mul_add(px, -(sa * py));
            let sy = sa.mul_add(px, ca * py);
            if inside_section(&foil, sx, sy) {
                let i = grid.idx(gx, gy);
                grid.flags[i] = Cell::Wall;
            }
        }
    }
}

/// Point-in-polygon (even-odd) against the closed section.
fn inside_section(foil: &fs_bem::panel2d::Airfoil2d, x: f64, y: f64) -> bool {
    let n = foil.nodes().len();
    let mut inside = false;
    for i in 0..n {
        let a = foil.nodes()[i];
        let b = foil.nodes()[(i + 1) % n];
        if (a[1] > y) != (b[1] > y) {
            let t = (y - a[1]) / (b[1] - a[1]);
            let xc = t.mul_add(b[0] - a[0], a[0]);
            if x < xc {
                inside = !inside;
            }
        }
    }
    inside
}

/// Run the LBM refinement of one candidate.
///
/// # Panics
/// Only on fs-lbm programmer contracts (fixture-scale).
#[must_use]
pub fn refine(c: &OrnithCandidate) -> RefineReport {
    let (nx, ny) = (96usize, 48usize);
    let chord = 24.0;
    let (x0, y0) = (30.0f64, 24.0f64);
    let u_in = 0.06;
    let mut grid = Grid::uniform(nx, ny, 0.56);
    grid.periodic_x = false;
    grid.periodic_y = true;
    rasterize(&mut grid, c, chord, x0, y0);
    // Only the section's walls contribute to the force receipt; inlet/outlet
    // momentum injection is not part of the aerodynamic force on the section.
    let section_walls: Vec<bool> = grid.flags.iter().map(|&cell| cell == Cell::Wall).collect();
    // Freestream init.
    for y in 0..ny {
        for x in 0..nx {
            let i = grid.idx(x, y);
            if grid.flags[i] == Cell::Fluid {
                grid.f[i] = equilibrium(1.0, u_in, 0.0);
            }
        }
    }
    let mut scratch: Vec<[f64; Q]> = Vec::new();
    let mut prev_u = vec![[0.0f64; 2]; nx * ny];
    let mut steadiness = f64::INFINITY;
    let mut force_sum = [0.0f64; 2];
    // Two full flow-throughs (transit = nx/u_in = 1600 steps): one
    // transit was MEASURED unsettled (steadiness 1.1e-3).
    for step in 0..FLOW_STEPS {
        // Inlet/outlet columns pinned to freestream equilibrium (the
        // smoke boundary treatment; documented).
        for y in 0..ny {
            let i0 = grid.idx(0, y);
            let i1 = grid.idx(nx - 1, y);
            if grid.flags[i0] == Cell::Fluid {
                grid.f[i0] = equilibrium(1.0, u_in, 0.0);
            }
            if grid.flags[i1] == Cell::Fluid {
                grid.f[i1] = equilibrium(1.0, u_in, 0.0);
            }
        }
        if step >= FLOW_STEPS - FORCE_SAMPLES {
            let exchange = grid.step_with_wall_momentum(&mut scratch, &section_walls);
            force_sum[0] += exchange.wall_impulse[0];
            force_sum[1] += exchange.wall_impulse[1];
        } else {
            grid.step(&mut scratch);
        }
        if step % 100 == 99 {
            let mut worst = 0.0f64;
            for y in 0..ny {
                for x in 0..nx {
                    let i = grid.idx(x, y);
                    if grid.flags[i] == Cell::Fluid {
                        let m = grid.moments(i);
                        let du = (m.u[0] - prev_u[i][0])
                            .abs()
                            .max((m.u[1] - prev_u[i][1]).abs());
                        worst = worst.max(du);
                        prev_u[i] = m.u;
                    }
                }
            }
            steadiness = worst;
        }
    }
    // Lattice dt = 1: each wall impulse is the mean force over one step.
    // Average the settled tail to reduce the rasterized obstacle's fluctuations.
    let drag = force_sum[0] / FORCE_SAMPLES as f64;
    let lift = force_sum[1] / FORCE_SAMPLES as f64;
    // Freestream rho = 1 and unit span: q_inf * c = 0.5 * U_inf^2 * c.
    let reference_force = 0.5 * u_in * u_in * chord;
    let panel_cl = fs_bem::panel2d::solve(&c.section(64), c.alpha)
        .expect("ornith candidate must satisfy the bounded panel-solve contract")
        .cl;
    RefineReport {
        lift,
        drag,
        cl: lift / reference_force,
        cd: drag / reference_force,
        panel_cl,
        honesty: "panel (inviscid, attached) vs D2Q9 BGK low-Re channel: native wall \
                  momentum exchange averaged over the final 100 steps; SIGN and \
                  order-band agreement only; LES model cards are the successor",
        steadiness,
    }
}
