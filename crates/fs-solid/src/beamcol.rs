//! Force-based distributed-plasticity beam-columns (bead tfz.14,
//! Spacone-style state determination): the equilibrium interpolation
//! is EXACT — for a cantilever under tip shear, `N(x) = N_tip` and
//! `M(x) = V·(L − x)` hold identically — so plasticity spreads through
//! the fiber sections at the integration points without displacement
//! shape-function error. Element flexibility integrates the section
//! flexibilities along the length; displacement control drives the
//! pushover.

use crate::SolidError;
use crate::fiber::Section;

/// A force-based cantilever element with fiber sections at
/// Gauss–Lobatto points (Lobatto puts a point AT the base hinge).
#[derive(Debug, Clone)]
pub struct ForceBasedElement {
    /// Length.
    pub length: f64,
    /// Sections at the integration points (base first).
    pub sections: Vec<Section>,
    /// Lobatto abscissae in [0, 1] (from the base).
    pub xi: Vec<f64>,
    /// Lobatto weights (sum 1).
    pub wi: Vec<f64>,
    /// Committed section strain pairs (ε₀, κ) per point.
    pub committed: Vec<(f64, f64)>,
}

impl ForceBasedElement {
    /// Five-point Gauss–Lobatto cantilever element.
    #[must_use]
    pub fn new(length: f64, make_section: &dyn Fn() -> Section) -> ForceBasedElement {
        let xi = vec![0.0, 0.172_673_164_646_011, 0.5, 0.827_326_835_353_989, 1.0];
        let wi = vec![
            0.05,
            0.272_222_222_222_222,
            0.355_555_555_555_556,
            0.272_222_222_222_222,
            0.05,
        ];
        ForceBasedElement {
            length,
            sections: (0..xi.len()).map(|_| make_section()).collect(),
            committed: vec![(0.0, 0.0); xi.len()],
            xi,
            wi,
        }
    }

    /// Tip deflection under tip shear `v` (state determination:
    /// per-section Newton on (ε₀, κ) against the EXACT internal forces
    /// N = 0, M = v·(L − x)), then curvature integration by virtual
    /// work. Commits the section states.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] if a section refuses to converge.
    pub fn tip_deflection_under_shear(&mut self, v: f64) -> Result<f64, SolidError> {
        let l = self.length;
        let mut defl = 0.0;
        for p in 0..self.xi.len() {
            let x = self.xi[p] * l;
            let m_target = v * (l - x);
            // Section Newton: find (ε₀, κ) with N = 0, M = m_target.
            let (mut e0, mut k0) = self.committed[p];
            let mut converged = false;
            let mut history = Vec::new();
            for _ in 0..60 {
                let r = self.sections[p].respond(e0, k0);
                let rn = [r.n, r.m - m_target];
                let norm = rn[0].hypot(rn[1]);
                history.push(norm);
                let scale = 1.0 + m_target.abs();
                if norm < 1e-8 * scale {
                    converged = true;
                    break;
                }
                let t = r.tangent;
                let det = t[0][0] * t[1][1] - t[0][1] * t[1][0];
                if det.abs() < 1e-30 {
                    break;
                }
                let de = (t[1][1] * rn[0] - t[0][1] * rn[1]) / det;
                let dk = (-t[1][0] * rn[0] + t[0][0] * rn[1]) / det;
                // Backtracking on the residual norm (fiber laws kink
                // at reversals; fixed damping oscillates there).
                let mut alpha = 1.0f64;
                for _ in 0..25 {
                    let r2 = self.sections[p].respond(e0 - alpha * de, k0 - alpha * dk);
                    if (r2.n).hypot(r2.m - m_target) < norm {
                        break;
                    }
                    alpha *= 0.5;
                }
                e0 -= alpha * de;
                k0 -= alpha * dk;
            }
            if !converged {
                return Err(SolidError::NewtonStalled { history });
            }
            self.sections[p].commit(e0, k0);
            self.committed[p] = (e0, k0);
            // Virtual work: δ_tip = ∫ κ(x)·(L − x) dx.
            defl += self.wi[p] * l * k0 * (l - x);
        }
        Ok(defl)
    }

    /// Energy dissipated between two committed sweeps of a cyclic
    /// moment-curvature history at the base section (trapezoidal work
    /// integral of the base hinge).
    #[must_use]
    pub fn base_committed(&self) -> (f64, f64) {
        self.committed[0]
    }
}

/// One pushover ledger row.
#[derive(Debug, Clone, Copy)]
pub struct PushoverStep {
    /// Applied tip shear.
    pub shear: f64,
    /// Measured tip deflection.
    pub deflection: f64,
    /// Base-section curvature.
    pub base_curvature: f64,
}
