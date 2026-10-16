//! Nonlinear films in the same impact solve as the stick, supports and air.
//!
//! fs-plate owns the prestressed pencil and statically relaxed von Karman
//! energy. fs-la factors its cold in-plane problem once. This module admits
//! that potential to the existing impact host; it adds neither an integrator
//! nor an amplitude-to-pitch control curve. The caller supplies damping and
//! launch motion separately through ImpactBody.
use std::sync::Arc;
use fs_plate::{ModePair, PlateError, PlateMesh, PlateModel, PlateSection};
use fs_plate::shell::head::nonlinear::{MembraneReduction, MembraneReductionBudget};
use super::{ImpactError, MAX_IMPACT_MODES};

/// Geometry-owned stretching energy and an explicit small-slope validity limit.
/// Clones share only immutable reduction data; every dynamic state remains in
/// its ImpactSystem. Internal Newton trials may extrapolate the law, but an
/// accepted initial/endpoint state must pass the declared slope limit.
#[derive(Debug, Clone)]
pub struct MembranePotential {
    reduction: Arc<MembraneReduction>,
    maximum_slope: f64,
}

/// Observation of an accepted body state, not a separate audio signal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MembraneObservation {
    /// Largest transverse P1 slope anywhere on the supplied reference mesh.
    pub maximum_slope: f64,
    /// Additional stretching storage, excluding the original tension/bending.
    pub stretching_energy_j: f64,
}

impl MembranePotential {
    /// Adopt an already prepared fs-plate reduction without changing its basis.
    /// A slope cap is a caller-selected model restriction, not an error bound.
    ///
    /// # Errors
    /// Empty/oversized basis, or a nonfinite cap outside the open interval (0,1).
    pub fn new(reduction: MembraneReduction, maximum_slope: f64) -> Result<Self, ImpactError> {
        if reduction.mode_count() == 0 || reduction.mode_count() > MAX_IMPACT_MODES
            || !maximum_slope.is_finite() || maximum_slope <= 0.0 || maximum_slope >= 1.0
        {
            return Err(ImpactError::Invalid("membrane impact needs bounded modes and an explicit slope limit in (0,1)"));
        }
        Ok(Self { reduction: Arc::new(reduction), maximum_slope })
    }

    /// Reduce an arbitrary planar film, retaining its actual modes and rim.
    /// Both incremental in-plane translations at `fixed_nodes` are constrained;
    /// all other nodes relax through the existing static elasticity solve.
    /// The original pencil supplies installed tension and bending exactly once.
    ///
    /// # Errors
    /// Existing fs-plate admission, a singular in-plane stiffness, insufficient
    /// cold-work budget, or an invalid impact slope/basis limit.
    #[allow(clippy::too_many_arguments)] // one geometry/basis/material admission
    pub fn from_pencil(
        mesh: &PlateMesh, section: &PlateSection, model: &PlateModel,
        modes: &[ModePair], fixed_nodes: &[usize], budget: MembraneReductionBudget,
        maximum_slope: f64,
    ) -> Result<Self, ImpactError> {
        if modes.is_empty() || modes.len() > MAX_IMPACT_MODES
            || !maximum_slope.is_finite() || maximum_slope <= 0.0 || maximum_slope >= 1.0
        {
            return Err(ImpactError::Invalid("membrane impact basis/slope admission failed before reduction"));
        }
        let reduction = MembraneReduction::new(mesh, section, model, modes, fixed_nodes, budget,
            |dimension, stiffness, loads| {
                let factor = fs_la::factor::cholesky(stiffness, dimension).map_err(|_|
                    PlateError::BadSection { what: "in-plane film stiffness is not positive definite under the supplied supports" })?;
                let mut solutions = loads.to_vec();
                for displacement in &mut solutions { factor.solve(displacement); }
                Ok(solutions)
            }).map_err(|e| ImpactError::Owner(e.to_string()))?;
        Self::new(reduction, maximum_slope)
    }

    /// Unmodified geometric reduction, including its cold-solve residual.
    #[must_use]
    pub fn reduction(&self) -> &MembraneReduction { &self.reduction }

    /// Largest admitted physical endpoint slope (dimensionless).
    #[must_use]
    pub const fn slope_limit(&self) -> f64 { self.maximum_slope }

    /// Inspect a complete mass-normalized displacement vector without mutation.
    /// Refusing is preferable to silently clamping strain, force or pitch.
    pub fn observe(&self, q: &[f64]) -> Result<MembraneObservation, ImpactError> {
        let slope = self.reduction.maximum_slope(q);
        let energy = self.reduction.stretching_energy(q);
        if !slope.is_finite() || slope > self.maximum_slope
            || !energy.is_finite() || energy < 0.0
        {
            return Err(ImpactError::Invalid("membrane state exceeds its declared slope or finite stretching-energy validity"));
        }
        Ok(MembraneObservation { maximum_slope: slope, stretching_energy_j: energy })
    }

    pub(super) fn observe_interleaved(&self, state: &[f64], offset: usize)
        -> Result<MembraneObservation, ImpactError>
    {
        let count = self.reduction.mode_count();
        let end = offset.checked_add(count).and_then(|n| n.checked_mul(2))
            .ok_or(ImpactError::Invalid("membrane state address overflow"))?;
        if end > state.len() { return Err(ImpactError::Invalid("membrane state does not contain its complete body")); }
        let mut q = [0.0; MAX_IMPACT_MODES];
        for i in 0..count { q[i] = state[2*(offset+i)]; }
        self.observe(&q[..count])
    }
}
