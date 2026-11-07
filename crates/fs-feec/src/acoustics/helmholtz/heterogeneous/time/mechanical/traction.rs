//! Applied surface traction in the original mechanical wall trace space.
//!
//! A dimensionless shape multiplies a physical pressure amplitude [Pa]. Positive
//! traction acts OUT of the fluid, like positive interior pressure. Its weak
//! force is W*shape*amplitude and its conjugate velocity is shape^T W*v [m^3/s].
//! This excites the existing wall momentum equation; it does not prescribe wall
//! motion or inject an additional fluid volume source.

use super::{HeterogeneousError, PreparedMechanicalTime, invalid, poll};
use std::collections::BTreeSet;

/// One physical traction shape on one original mechanical patch.
#[derive(Clone, Debug)]
pub struct WallTractionProfile {
    /// Index in the original `MechanicalBoundary` preparation order.
    pub wall: usize,
    /// Dimensionless coefficients indexed by ORIGINAL GLOBAL hierarchical FE
    /// DOF, not patch-local index or a nodal sample. Omitted coefficients are
    /// zero. Every named coordinate must belong to this patch, even if its
    /// fluid pressure is constrained by an adjacent pressure-release face.
    pub coefficients: Vec<(usize, f64)>,
}

impl WallTractionProfile {
    pub(super) fn local_coefficients(
        &self,
        dofs: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, HeterogeneousError> {
        poll(cancelled)?;
        if self.coefficients.len() > dofs.len() {
            return Err(invalid("traction profile exceeds its original wall trace extent").into());
        }
        let mut result = vec![0.0; dofs.len()];
        let mut seen = BTreeSet::new();
        for &(global, coefficient) in &self.coefficients {
            poll(cancelled)?;
            if !coefficient.is_finite() || !seen.insert(global) {
                return Err(
                    invalid("finite unique original traction coefficients required").into(),
                );
            }
            let local = dofs
                .binary_search(&global)
                .map_err(|_| invalid("traction coefficient is outside its mechanical patch"))?;
            result[local] = coefficient;
        }
        poll(cancelled)?;
        Ok(result)
    }
}

/// A reusable shape bound to one prepared wall law and original geometry.
/// Preparation needs O(wall DOFs) storage. Warm application neither assembles
/// a surface matrix nor factors a new acoustic or mechanical operator.
pub struct PreparedMechanicalTraction<'p, 'm, 'f, 'c> {
    pub(super) owner: &'p PreparedMechanicalTime<'m, 'f, 'c>,
    pub(super) wall: usize,
    pub(super) coefficients: Vec<f64>,
}

impl PreparedMechanicalTraction<'_, '_, '_, '_> {
    #[must_use]
    pub const fn wall(&self) -> usize {
        self.wall
    }

    /// Shape coefficients in the prepared wall's `wall_dofs` order.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }
}

impl<'m, 'f, 'c> PreparedMechanicalTime<'m, 'f, 'c> {
    /// Bind an applied traction shape to the same complete physical wall that
    /// exchanges pressure work with the cavity. This changes no material state.
    pub fn prepare_wall_traction(
        &self,
        profile: &WallTractionProfile,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTraction<'_, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let wall = self
            .walls
            .get(profile.wall)
            .ok_or_else(|| invalid("traction mechanical wall index"))?;
        let coefficients = profile.local_coefficients(&wall.dofs, cancelled)?;
        Ok(PreparedMechanicalTraction {
            owner: self,
            wall: profile.wall,
            coefficients,
        })
    }
}
