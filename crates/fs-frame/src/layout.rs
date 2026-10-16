//! Stages 1–2: ground-structure layout (PDHG LP diagnostics plus an outward
//! optimum certificate) and code-checked sizing (Euler floors, catalog
//! up-snap, mandatory post-prune equilibrium re-verification, member
//! code rows) — thin, honest composition over fs-truss, whose own
//! battery (truss-001..006) carries the per-stage evidence. The
//! Michell-continuum catalogue comparison stays LEDGERED PENDING
//! exactly as in the fs-truss contract.

use fs_evidence::Color;
use fs_exec::Cx;
use fs_truss::ground::{GroundLimits, GroundRules, GroundStructure, TrussConstructionError};
use fs_truss::lp::{LayoutCase, LayoutLimits, LayoutLp, PdhgError, PdhgSettings};
use fs_truss::sizing::{CatalogAudit, size_and_snap};
use fs_truss::{LayoutCertificateError, LayoutCertificateLimits, LayoutCertificateProblem};
use fs_truss_e2e::{optimality_color_from_certificate, rescale_optimality_color};

/// Structured failure from the frame layout and sizing composition.
#[derive(Debug)]
pub enum LayoutError {
    /// Ground-structure, load-case, or LP construction was refused.
    Construction(TrussConstructionError),
    /// The admitted PDHG solve or its returned-state diagnostics failed.
    Solver(PdhgError),
    /// The outward optimum-certificate stage encountered malformed state,
    /// allocation failure, or cancellation.
    Certificate(LayoutCertificateError),
}

impl core::fmt::Display for LayoutError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Construction(error) => {
                write!(formatter, "frame layout construction failed: {error}")
            }
            Self::Solver(error) => write!(formatter, "frame layout solve failed: {error}"),
            Self::Certificate(error) => {
                write!(formatter, "frame layout certificate failed: {error}")
            }
        }
    }
}

impl std::error::Error for LayoutError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Construction(error) => Some(error),
            Self::Solver(error) => Some(error),
            Self::Certificate(error) => Some(error),
        }
    }
}

impl From<TrussConstructionError> for LayoutError {
    fn from(error: TrussConstructionError) -> Self {
        Self::Construction(error)
    }
}

impl From<PdhgError> for LayoutError {
    fn from(error: PdhgError) -> Self {
        Self::Solver(error)
    }
}

impl From<LayoutCertificateError> for LayoutError {
    fn from(error: LayoutCertificateError) -> Self {
        Self::Certificate(error)
    }
}

/// The layout + sizing record: approximate solution, diagnostics, and audit.
pub struct LayoutReport {
    /// The ground structure.
    pub gs: GroundStructure,
    /// Split LP solution (q⁺, q⁻).
    pub x: Vec<f64>,
    /// Relative primal/dual objective separation diagnostic.
    pub gap: f64,
    /// Equilibrium residual.
    pub residual: f64,
    /// Approximate returned-iterate volume Σ l·|q|/σ_y.
    pub volume: f64,
    /// Outward-certified physical optimum-volume interval, or an honest
    /// diagnostic estimate when the proof is numerically unavailable.
    pub optimality_color: Color,
    /// Sizing/catalog audit (Euler floors, snap, code rows).
    pub audit: CatalogAudit,
}

/// Run layout + sizing for a grid ground structure with a downward
/// unit load at the top-right node and pins on the left edge (the
/// smoke-tier cantilever fixture; the geometry knobs are arguments so
/// study programs can vary them).
///
/// # Errors
/// Returns [`LayoutError::Construction`] when the geometry, physical
/// parameters, catalog, construction budgets, or cancellation state are not
/// admitted. Returns [`LayoutError::Solver`] when the PDHG controls or its
/// returned state are invalid. Returns [`LayoutError::Certificate`] when the
/// outward proof stage encounters malformed state, allocation failure, or
/// cancellation.
#[allow(clippy::too_many_arguments)] // Five physical inputs, catalog, and explicit execution context.
pub fn layout_and_size(
    nx: usize,
    ny: usize,
    w: f64,
    h: f64,
    sigma_y: f64,
    youngs: f64,
    catalog: &[f64],
    cx: &Cx<'_>,
) -> Result<LayoutReport, LayoutError> {
    if !sigma_y.is_finite() || sigma_y <= 0.0 {
        return Err(TrussConstructionError::InvalidInput {
            field: "sigma_y",
            requirement: "must be finite and positive",
        }
        .into());
    }
    if !youngs.is_finite() || youngs <= 0.0 {
        return Err(TrussConstructionError::InvalidInput {
            field: "youngs",
            requirement: "must be finite and positive",
        }
        .into());
    }
    if catalog.is_empty()
        || catalog.iter().any(|area| !area.is_finite() || *area <= 0.0)
        || catalog.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(TrussConstructionError::InvalidInput {
            field: "catalog",
            requirement: "must contain finite positive areas in strictly increasing order",
        }
        .into());
    }
    let rules = GroundRules::default();
    let gs = GroundStructure::try_grid(nx, ny, w, h, &rules, GroundLimits::default(), cx)?;
    // Pins: left edge; load: unit down at the top-right node.
    let n_nodes = gs.nodes().len();
    let mut supported = vec![[false; 2]; n_nodes];
    for (i, p) in gs.nodes().iter().enumerate() {
        if p[0] <= 1e-12 {
            supported[i] = [true; 2];
        }
    }
    let mut load = vec![[0.0f64; 2]; n_nodes];
    let mut best = 0usize;
    let mut score = f64::NEG_INFINITY;
    for (i, p) in gs.nodes().iter().enumerate() {
        let s = p[0] + p[1];
        if s > score {
            score = s;
            best = i;
        }
    }
    load[best] = [0.0, -1.0];
    let case = LayoutCase::try_new(supported, load, n_nodes)?;
    // Solve the LP at σ_y = 1: the yield stress only scales the
    // objective (c = l/σ_y), never the equilibrium constraints, so
    // the optimal forces are σ_y-independent — and a 250 MPa σ_y
    // makes c ~ 1e-9 against unit loads, stalling PDHG's primal-dual
    // scaling (measured: gap stuck at 1.0). Volume is rescaled to
    // physical units on report; sizing gets the TRUE σ_y.
    let lp = LayoutLp::try_assemble(&gs, &case, 1.0, LayoutLimits::default(), cx)?;
    let settings = PdhgSettings::default();
    let (x, y, mut pdhg_report) = lp.solve(None, None, settings)?;
    let bnorm = 1.0;
    let (gap, residual, volume_unit) = lp.diagnostics(&x, &y, bnorm)?;
    let certificate_status = lp.certify_optimum_for_report(
        &x,
        &y,
        settings,
        &mut pdhg_report,
        LayoutCertificateLimits::default(),
        cx,
    )?;
    let certificate_problem = LayoutCertificateProblem::try_new(lp.a(), lp.c(), lp.b())?;
    let normalized_optimality = optimality_color_from_certificate(
        &certificate_problem,
        &x,
        &y,
        settings,
        &certificate_status,
        gap,
        residual,
        cx,
    )?;
    let optimality_color = rescale_optimality_color(&normalized_optimality, sigma_y);
    let volume = volume_unit / sigma_y;
    let audit = size_and_snap(&gs, &lp, &x, sigma_y, youngs, catalog, 1e-3);
    Ok(LayoutReport {
        gs,
        x,
        gap,
        residual,
        volume,
        optimality_color,
        audit,
    })
}
