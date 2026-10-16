//! Select the existing radiation formulation from the actual closed solids.
//!
//! The conventional exterior equation loses uniqueness at interior Dirichlet
//! eigenvalues. For a component contained in a box of side lengths L_i,
//! domain monotonicity gives k_1 >= pi * sqrt(sum_i 1/L_i^2). Separate solids
//! have separate interior spectra: their separation and a distant lid do NOT
//! create a fictitious interior mode of the enclosing empty space.
//!
//! Use triangle CBIE below 80% of the smallest component bound; otherwise use
//! the existing Burton--Miller arm. This replaces a k*scene_radius heuristic,
//! not a failed-solve retry, positivity repair or added physical damping.
//! All original wavelength, solve and radiation-power checks still apply.
//!
//! The bound excludes continuum fictitious resonances only. It does NOT bound
//! collocation/quadrature error, certify passivity of a discrete matrix, or
//! validate intersecting/nested solids. Components must be disjoint, outward
//! and non-self-intersecting, as in the existing exterior surface contract.
//!
//! References: "An improved form of the hypersingular boundary integral
//! equation for exterior acoustic problems", EABE 34 (2010), 189--195,
//! doi:10.1016/j.enganabound.2009.10.005 (interior Dirichlet nonuniqueness);
//! V. Ivrii, PDE textbook, sections 13.2--13.3 (box spectrum and min--max).
use crate::{
    helmholtz::{self, Formulation, HelmholtzError, Medium, RadiationSolution},
    near_field, panel3d::SpherePanels,
};
use fs_math::c64::C64;

/// Keep a fixed 20% wavenumber margin below the geometric exclusion bound.
/// This is a formulation-selection margin, not a power/error tolerance.
pub const CBIE_SPECTRAL_MARGIN: f64 = 0.8;

/// Geometry-bound policy, reusable across a frequency grid. The borrow keeps
/// a selected spectrum and its source triangles from being cross-wired.
pub struct GeometryPolicy<'a> {
    surface: &'a SpherePanels,
    component_bounds: Vec<f64>,
    plain_limit: f64,
}
fn bad(what: &'static str) -> HelmholtzError { HelmholtzError::BadParameter { what } }
impl<'a> GeometryPolicy<'a> {
    /// Derive each component's conservative axis-aligned box eigenvalue bound.
    /// The exact-coordinate closure/orientation owner is shared with receiver
    /// admission; no tolerance weld, geometry inflation, or missing-face fill.
    /// Bounds are rounded down and box widths up at every positive operation.
    /// # Errors
    /// Missing retained triangles, open/inward/degenerate or unbounded geometry,
    /// exceeded panel cap, or an unrepresentable positive spectral bound.
    pub fn new(surface: &'a SpherePanels) -> Result<Self, HelmholtzError> {
        let triangles = surface.triangles().ok_or_else(|| bad("geometric radiation policy requires retained triangles"))?;
        if triangles.is_empty() || triangles.len() > helmholtz::MAX_DENSE_PANELS
            || triangles.iter().flatten().flatten().any(|v| !v.is_finite() || v.abs() > 1e9) {
            return Err(bad("invalid bounded triangle surface for radiation policy"));
        }
        let components = near_field::components(triangles)?;
        let mut component_bounds = Vec::with_capacity(components.len());
        for component in components {
            let mut lo = [f64::INFINITY; 3]; let mut hi = [f64::NEG_INFINITY; 3];
            for i in component { for p in triangles[i] { for c in 0..3 {
                lo[c] = lo[c].min(p[c]); hi[c] = hi[c].max(p[c]);
            }}}
            let mut widths: [f64; 3] = std::array::from_fn(|c| (hi[c] - lo[c]).next_up());
            if widths.iter().any(|w| !w.is_finite() || *w <= 0.) {
                return Err(bad("radiating solid has unresolved enclosing box"));
            }
            // Fixed summation order also gives axis-permutation parity.
            widths.sort_by(f64::total_cmp);
            let mut inverse_square_sum = 0.;
            for width in widths {
                let inverse = (1. / width).next_down();
                let term = (inverse * inverse).next_down();
                inverse_square_sum = (inverse_square_sum + term).next_down();
            }
            let lower = (std::f64::consts::PI.next_down()
                * inverse_square_sum.sqrt().next_down()).next_down();
            if !lower.is_finite() || lower <= 0. {
                return Err(bad("unrepresentable interior Dirichlet exclusion bound"));
            }
            component_bounds.push(lower);
        }
        let smallest = component_bounds.iter().copied().fold(f64::INFINITY, f64::min);
        let plain_limit = (CBIE_SPECTRAL_MARGIN.next_down() * smallest).next_down();
        if !plain_limit.is_finite() || plain_limit <= 0. {
            return Err(bad("empty or unresolved radiation-selection band"));
        }
        Ok(Self { surface, component_bounds, plain_limit })
    }
    /// Component-wise lower bounds on the FIRST interior Dirichlet wavenumber
    /// [rad/m], not measured resonances. Axis-aligned boxes can be conservative
    /// after a general rigid rotation; validity is unchanged, tightness is not.
    pub fn first_dirichlet_bounds(&self) -> &[f64] { &self.component_bounds }
    /// Maximum k [rad/m] selecting the existing triangle-CBIE image.
    pub fn plain_cbie_limit(&self) -> f64 { self.plain_limit }
    /// Select before solving. A failed selected solve is never retried in the
    /// other formulation and never repaired by altering its radiation power.
    /// # Errors
    /// Nonpositive or nonfinite wavenumber.
    pub fn formulation(&self, k: f64) -> Result<Formulation, HelmholtzError> {
        if !k.is_finite() || k <= 0. { return Err(bad("radiation policy requires positive finite k")); }
        Ok(if k <= self.plain_limit { Formulation::PlainCbie } else { Formulation::BurtonMiller })
    }
    /// Solve all fields with the original shared Helmholtz factorization.
    /// # Errors
    /// Original batch, medium, geometry-resolution and factorization refusals.
    pub fn solve_batch(&self, k: f64, medium: Medium, fields: &[&[C64]])
        -> Result<Vec<RadiationSolution>, HelmholtzError> {
        helmholtz::solve_radiation_batch(self.surface, k, medium, fields, self.formulation(k)?)
    }
}

#[cfg(test)]
#[path = "radiation_policy_tests.rs"]
mod tests;
