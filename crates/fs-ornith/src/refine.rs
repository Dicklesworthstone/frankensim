//! Stage 3: REFINE. The screen winner's section is rasterized into an
//! fs-lbm channel as a bounce-back obstacle and driven to a steady
//! state; aerodynamic forces come from a CONTROL-VOLUME momentum
//! balance over the PUBLIC cell moments (∮ (ρuu + p·I)·n dA with
//! p = ρ·c_s² — no reliance on fs-lbm internals). The gate is
//! MODEL-FORM honest: the panel method and the low-Re lattice flow are
//! different physics, so agreement is gated on lift SIGN and
//! order-of-magnitude band, never point-matching; the honesty label
//! travels in the report. Cumulant LES on sparse VDB lattices with
//! DWR adaptivity is the recorded successor.

use crate::param::OrnithCandidate;
use fs_lbm::{CS2, Cell, Grid, Q, equilibrium};

/// The refinement report — with its honesty label.
#[derive(Debug, Clone)]
pub struct RefineReport {
    /// Control-volume lift (lattice units, +y).
    pub lift: f64,
    /// Control-volume drag (lattice units, +x).
    pub drag: f64,
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
    let n = foil.nodes.len();
    let mut inside = false;
    for i in 0..n {
        let a = foil.nodes[i];
        let b = foil.nodes[(i + 1) % n];
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
    // Two full flow-throughs (transit = nx/u_in = 1600 steps): one
    // transit was MEASURED unsettled (steadiness 1.1e-3).
    for step in 0..3200 {
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
        grid.step(&mut scratch);
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
    // Control-volume momentum balance on a box around the obstacle:
    // F = ∮ (ρ u (u·n) + p n) dA, force ON the body = −flux balance.
    let (bx0, bx1) = (12usize, 72usize);
    let (by0, by1) = (4usize, ny - 5);
    let mut fx = 0.0f64;
    let mut fy = 0.0f64;
    let mut add_face = |i: usize, n: [f64; 2], grid: &Grid| {
        let m = grid.moments(i);
        let p = m.rho * CS2;
        let un = m.u[0] * n[0] + m.u[1] * n[1];
        fx += m.rho * m.u[0] * un + p * n[0];
        fy += m.rho * m.u[1] * un + p * n[1];
    };
    for y in by0..=by1 {
        add_face(grid.idx(bx0, y), [-1.0, 0.0], &grid);
        add_face(grid.idx(bx1, y), [1.0, 0.0], &grid);
    }
    for x in bx0..=bx1 {
        add_face(grid.idx(x, by0), [0.0, -1.0], &grid);
        add_face(grid.idx(x, by1), [0.0, 1.0], &grid);
    }
    let panel_cl = fs_bem::panel2d::solve(&c.section(64), c.alpha).cl;
    RefineReport {
        lift: -fy,
        drag: -fx,
        panel_cl,
        honesty: "panel (inviscid, attached) vs D2Q9 BGK low-Re channel: SIGN and \
                  order-band agreement only; LES model cards are the successor",
        steadiness,
    }
}
