//! fs-spectral — spectral health monitoring (plan addendum, Proposal 5).
//! Layer: L1.
//!
//! The fixed-iteration gauge-fit/candidate-remainder triage used by MERGE
//! (Proposal 10) and sheaf repair is only as trustworthy as the sheaf
//! Laplacian's SPECTRAL GAP. The gap conditions numerical confidence; it does
//! not certify H¹ or turn a remainder into a structural obstruction. A system
//! that can say "I am becoming
//! numerically unreliable in this region, and here is the reroute" is a
//! different species from one that silently degrades. This crate:
//!
//! - computes a legacy, unauthoritative λ-gap diagnostic for a (small, dense,
//!   symmetric) operator via an in-house Jacobi eigensolver —
//!   [`symmetric_eigenvalues`], [`spectral_gap`];
//! - classifies gap HEALTH with HYSTERESIS thresholds so a marginal region does
//!   not flap between healthy and degraded — [`GapHealthMonitor`];
//! - PROPAGATES low confidence: a degraded gap DEMOTES any downstream color
//!   (verified/validated → estimated) and never promotes — [`propagate`]. This
//!   is the flag Proposal 10's merge outputs MUST surface;
//! - composes per-op AMPLIFICATION factors into an end-to-end CONDITIONING
//!   estimate and gives the representation router a conditioning term, so it
//!   prefers paths that keep the whole pipeline well-posed —
//!   [`compose_conditioning`], [`route`].
//! - validates versioned spectral problem descriptors, proposition-bound
//!   authority, and exact method-family prerequisites — [`admission`];
//! - validates one-way, content-addressed physical-domain crosswalks without
//!   importing higher-layer physical models or inventing reverse maps —
//!   [`adapter`];
//! - validates orthogonal set-valued result truth (authority, algebraic
//!   coverage, cluster/internal separation, projective accounting, and
//!   termination) against the complete admitted problem — [`truth`].
//!
//! Everything is deterministic. `#![deny(unsafe_code)]` via the workspace lint.

use fs_evidence::{Color, ColorRank};

pub mod adapter;
pub mod admission;
pub mod service;
pub mod truth;

/// A structured spectral failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpectralError {
    /// The matrix is not square.
    NotSquare,
    /// The matrix is not symmetric (within tolerance).
    NotSymmetric,
    /// The matrix is empty.
    Empty,
}

/// Eigenvalues of a small dense SYMMETRIC matrix, ascending, via the cyclic
/// Jacobi rotation algorithm (pure Rust; production uses the BEDROCK sparse
/// eigensolver with warm starts).
///
/// # Errors
/// [`SpectralError`] if the matrix is empty, non-square, or non-symmetric.
// A dense Jacobi kernel: the transpose reads (`m[j][i]`) and two-column plane
// rotations (`m[k][p]`/`m[k][q]`) are inherently 2D-indexed; index loops are the
// correct, readable form here.
#[allow(clippy::needless_range_loop)]
pub fn symmetric_eigenvalues(matrix: &[Vec<f64>]) -> Result<Vec<f64>, SpectralError> {
    let n = matrix.len();
    if n == 0 {
        return Err(SpectralError::Empty);
    }
    for row in matrix {
        if row.len() != n {
            return Err(SpectralError::NotSquare);
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if (matrix[i][j] - matrix[j][i]).abs() > 1e-9 {
                return Err(SpectralError::NotSymmetric);
            }
        }
    }
    let mut m: Vec<Vec<f64>> = matrix.to_vec();
    // cyclic Jacobi: zero out off-diagonals by plane rotations until converged.
    for _sweep in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += m[i][j] * m[i][j];
            }
        }
        if off <= 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if m[p][q].abs() <= 1e-18 {
                    continue;
                }
                let (app, aqq, apq) = (m[p][p], m[q][q], m[p][q]);
                // stable rotation: t = tan θ.
                let theta = (aqq - app) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // apply the plane rotation on both sides (M ← Jᵀ M J).
                for k in 0..n {
                    let (mkp, mkq) = (m[k][p], m[k][q]);
                    m[k][p] = c * mkp - s * mkq;
                    m[k][q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let (mpk, mqk) = (m[p][k], m[q][k]);
                    m[p][k] = c * mpk - s * mqk;
                    m[q][k] = s * mpk + c * mqk;
                }
            }
        }
    }
    let mut eig: Vec<f64> = (0..n).map(|i| m[i][i]).collect();
    eig.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Ok(eig)
}

/// A legacy unauthoritative spectral-gap reading: the gap above the smallest
/// eigenvalue (the Fiedler / algebraic-connectivity gap) relative to the
/// spectral spread.
///
/// This diagnostic cannot satisfy any [`admission`] witness or [`truth`]
/// proposition. In particular, callers must not treat its historical
/// zero-spread ratio as cluster-separation authority; RB.1b owns that numerical
/// correction while RB.1a keeps the typed truth surface isolated from it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralGap {
    /// The smallest eigenvalue.
    pub lambda_min: f64,
    /// The next eigenvalue.
    pub lambda_next: f64,
    /// `lambda_next - lambda_min`.
    pub gap: f64,
    /// The spectral spread `lambda_max - lambda_min`.
    pub spread: f64,
    /// `gap / spread` in `[0, 1]` — the dimensionless health ratio (1.0 if the
    /// spectrum is a single point).
    pub ratio: f64,
}

/// The spectral gap of an ascending eigenvalue list (needs ≥ 2 eigenvalues).
#[must_use]
pub fn spectral_gap(eigenvalues: &[f64]) -> Option<SpectralGap> {
    if eigenvalues.len() < 2 {
        return None;
    }
    let lambda_min = eigenvalues[0];
    let lambda_next = eigenvalues[1];
    let lambda_max = *eigenvalues.last().unwrap();
    let gap = lambda_next - lambda_min;
    let spread = lambda_max - lambda_min;
    let ratio = if spread <= f64::EPSILON {
        1.0 // a degenerate single-point spectrum is trivially "well separated"
    } else {
        (gap / spread).clamp(0.0, 1.0)
    };
    Some(SpectralGap {
        lambda_min,
        lambda_next,
        gap,
        spread,
        ratio,
    })
}

/// Gap health.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Health {
    /// The gap is well separated — triage is trustworthy.
    Healthy,
    /// The gap has collapsed — triage in this region is low-confidence.
    Degraded,
}

/// A hysteresis health monitor: a region degrades when its gap ratio falls to
/// `degrade_below` and only recovers above `restore_above` (`> degrade_below`),
/// so a marginal region cannot flap.
#[derive(Debug, Clone)]
pub struct GapHealthMonitor {
    degrade_below: f64,
    restore_above: f64,
    state: Health,
}

impl GapHealthMonitor {
    /// A monitor, starting Healthy.
    ///
    /// # Panics
    /// If `restore_above < degrade_below` (inverted hysteresis band).
    #[must_use]
    pub fn new(degrade_below: f64, restore_above: f64) -> GapHealthMonitor {
        assert!(
            restore_above >= degrade_below,
            "inverted hysteresis band: restore_above must be >= degrade_below"
        );
        GapHealthMonitor {
            degrade_below,
            restore_above,
            state: Health::Healthy,
        }
    }

    /// Feed a new gap ratio and return the (possibly updated) health.
    pub fn update(&mut self, ratio: f64) -> Health {
        self.state = match self.state {
            Health::Healthy if ratio <= self.degrade_below => Health::Degraded,
            Health::Degraded if ratio > self.restore_above => Health::Healthy,
            other => other,
        };
        self.state
    }

    /// The current health (no update).
    #[must_use]
    pub fn health(&self) -> Health {
        self.state
    }
}

/// Propagate spectral health into a downstream color: a degraded gap DEMOTES a
/// verified/validated color to estimated (low confidence) and never promotes;
/// a healthy gap leaves the color unchanged. Merge/triage outputs in a degraded
/// region MUST carry this demotion.
#[must_use]
pub fn propagate(color: Color, health: Health) -> Color {
    match health {
        Health::Healthy => color,
        Health::Degraded => {
            if color.rank() == ColorRank::Estimated {
                color // already lowest confidence
            } else {
                Color::Estimated {
                    estimator: "spectral-gap-degraded".to_string(),
                    dispersion: f64::INFINITY,
                }
            }
        }
    }
}

/// Compose per-op AMPLIFICATION factors (local condition estimates, each `>=
/// 0`) into an end-to-end CONDITIONING estimate — condition numbers multiply
/// along a pipeline. An empty pipeline is perfectly conditioned (`1.0`).
///
/// # Errors
/// [`SpectralError::NotSquare`] is reused to signal a negative/non-finite
/// amplification factor (an invalid conditioning input).
pub fn compose_conditioning(amplifications: &[f64]) -> Result<f64, SpectralError> {
    let mut product = 1.0;
    for &a in amplifications {
        if !(a.is_finite() && a >= 0.0) {
            return Err(SpectralError::NotSquare);
        }
        product *= a;
    }
    Ok(product)
}

/// A candidate representation-routing path with its cost and conditioning.
#[derive(Debug, Clone, PartialEq)]
pub struct RouterPath {
    /// A label.
    pub label: String,
    /// The base (cheapness) cost.
    pub base_cost: f64,
    /// The end-to-end conditioning estimate (`>= 1` well-posed; larger = worse).
    pub conditioning: f64,
}

/// Choose a routing path by fitness `base_cost + conditioning_weight ·
/// ln(max(conditioning, 1))`. With `conditioning_weight = 0` this is pure
/// cheapness; with a positive weight the router prefers better-conditioned
/// paths that keep the whole pipeline well-posed. Ties break on order.
#[must_use]
pub fn route(paths: &[RouterPath], conditioning_weight: f64) -> Option<&RouterPath> {
    paths.iter().min_by(|a, b| {
        let fa = a.base_cost + conditioning_weight * a.conditioning.max(1.0).ln();
        let fb = b.base_cost + conditioning_weight * b.conditioning.max(1.0).ln();
        fa.partial_cmp(&fb).unwrap_or(std::cmp::Ordering::Equal)
    })
}
