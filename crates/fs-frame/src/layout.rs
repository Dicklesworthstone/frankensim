//! Stages 1–2: ground-structure LAYOUT (PDHG LP with the duality gap
//! as CERTIFICATE) and code-checked SIZING (Euler floors, catalog
//! up-snap, mandatory post-prune equilibrium re-verification, member
//! code rows) — thin, honest composition over fs-truss, whose own
//! battery (truss-001..006) carries the per-stage evidence. The
//! Michell-continuum catalogue comparison stays LEDGERED PENDING
//! exactly as in the fs-truss contract.

use fs_truss::ground::{GroundRules, GroundStructure};
use fs_truss::lp::{LayoutLp, PdhgSettings};
use fs_truss::sizing::{CatalogAudit, size_and_snap};

/// The layout + sizing record: solution, certificate, audit.
pub struct LayoutReport {
    /// The ground structure.
    pub gs: GroundStructure,
    /// Split LP solution (q⁺, q⁻).
    pub x: Vec<f64>,
    /// Relative duality gap (the certificate).
    pub gap: f64,
    /// Equilibrium residual.
    pub residual: f64,
    /// Certified volume Σ l·|q|/σ_y.
    pub volume: f64,
    /// Sizing/catalog audit (Euler floors, snap, code rows).
    pub audit: CatalogAudit,
}

/// Run layout + sizing for a grid ground structure with a downward
/// unit load at the top-right node and pins on the left edge (the
/// smoke-tier cantilever fixture; the geometry knobs are arguments so
/// study programs can vary them).
///
/// # Panics
/// On degenerate grids (fs-truss programmer contracts).
#[must_use]
pub fn layout_and_size(
    nx: usize,
    ny: usize,
    w: f64,
    h: f64,
    sigma_y: f64,
    youngs: f64,
    catalog: &[f64],
) -> LayoutReport {
    let rules = GroundRules::default();
    let gs = GroundStructure::grid(nx, ny, w, h, &rules);
    // Pins: left edge; load: unit down at the top-right node.
    let n_nodes = gs.nodes.len();
    let mut pinned = vec![false; n_nodes];
    for (i, p) in gs.nodes.iter().enumerate() {
        if p[0] <= 1e-12 {
            pinned[i] = true;
        }
    }
    let mut load = vec![[0.0f64; 2]; n_nodes];
    let mut best = 0usize;
    let mut score = f64::NEG_INFINITY;
    for (i, p) in gs.nodes.iter().enumerate() {
        let s = p[0] + p[1];
        if s > score {
            score = s;
            best = i;
        }
    }
    load[best] = [0.0, -1.0];
    let pin_fn = |node: usize, _dof: usize| pinned[node];
    let load_fn = |node: usize| load[node];
    // Solve the LP at σ_y = 1: the yield stress only scales the
    // objective (c = l/σ_y), never the equilibrium constraints, so
    // the optimal forces are σ_y-independent — and a 250 MPa σ_y
    // makes c ~ 1e-9 against unit loads, stalling PDHG's primal-dual
    // scaling (measured: gap stuck at 1.0). Volume is rescaled to
    // physical units on report; sizing gets the TRUE σ_y.
    let lp = LayoutLp::assemble(&gs, &pin_fn, &load_fn, 1.0);
    let (x, y, _report) = lp.solve(None, None, PdhgSettings::default());
    let bnorm = 1.0;
    let (gap, residual, volume_unit) = lp.certificate(&x, &y, bnorm);
    let volume = volume_unit / sigma_y;
    let audit = size_and_snap(&gs, &lp, &x, sigma_y, youngs, catalog, 1e-3);
    LayoutReport {
        gs,
        x,
        gap,
        residual,
        volume,
        audit,
    }
}
