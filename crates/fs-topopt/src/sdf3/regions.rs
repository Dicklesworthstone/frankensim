//! Prescribed PHYSICAL material regions in the existing 3-D density map.
//!
//! The map is `physical = region_override(project(filter(raw)))`. A solid cell
//! has density one, a void cell density zero; the override has zero local
//! derivative. Apply that derivative BEFORE the filter transpose, not afterward:
//! a raw filter control located in a protected cell can still affect nearby
//! free cells. Raw controls are deliberately NOT prescribed material fractions.
//!
//! The physical mask enters the same stiffness, volume and sensitivity map used
//! by the numerical evaluators. Void retains the declared SIMP ersatz stiffness
//! `e_min` and does not remove quadrature, loads or boundary supports; use an
//! implicit-domain difference for an actual geometric hole. Prescribed solids
//! consume the material budget. Refinement inherits whole-cell labels through
//! the checked geometric parent map; it never reclassifies them by a new center.

use super::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::{EvaluationStop, SolveControl};

/// Physical material behavior of one active cut cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PhysicalRegion3 {
    /// Use the filtered and projected design control.
    #[default]
    Design,
    /// Preserve full material, independently of every raw design control.
    Solid,
    /// Preserve zero material density with the declared ersatz stiffness.
    Void,
}
impl PhysicalRegion3 {
    /// Stable human/machine-readable label for retained design exports.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self { Self::Design => "design", Self::Solid => "solid", Self::Void => "void" }
    }
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Bind immutable physical material labels in this operator's active-cell
    /// order. Call before optimizing; as with a freshly constructed study, a
    /// complete `evaluate` is required before treating operator scales as a
    /// solution at the new declaration. This performs no PDE solve.
    ///
    /// The filter and its raw controls are unchanged. Only the physical output
    /// and its local Jacobian are overridden. Every new continuation model
    /// therefore preserves the exact zero/one material assignment. The adaptive
    /// compliance driver automatically inherits labels to child cells. Other
    /// callers rebuilding a study must explicitly preserve/rebind their policy.
    ///
    /// # Errors
    /// Refuses a mask not matching the full active-cell count.
    pub fn with_physical_regions(mut self, regions: Vec<PhysicalRegion3>) -> Result<Self, EvaluationStop> {
        if regions.len() != self.cells() { return Err(failure("sdf3-physical-region-shape")); }
        self.physical_regions = Some(regions);
        Ok(self)
    }

    /// Exact retained labels, or `None` for the original all-design map.
    #[must_use]
    pub fn physical_regions(&self) -> Option<&[PhysicalRegion3]> {
        self.physical_regions.as_deref()
    }

    /// Numerical fraction of the domain occupied by prescribed solid cells.
    /// This is an unavoidable material cost, not a continuum-volume certificate.
    #[must_use]
    pub fn prescribed_solid_fraction(&self) -> f64 {
        self.physical_regions.as_ref().map_or(0.0, |regions| {
            regions.iter().zip(&self.mass).filter_map(|(region, mass)| {
                (*region == PhysicalRegion3::Solid).then_some(*mass)
            }).sum()
        })
    }
}

/// Copy labels through an already established fine-to-coarse active-cell map.
/// Source labels are borrowed and unchanged on every stop. The caller owns
/// geometric identity: this function cannot validate geometry from indices.
/// Invalid parents and cancellation produce no partial mask.
pub fn inherit_physical_regions(
    source: &[PhysicalRegion3], parents: &[usize], control: &mut SolveControl<'_>,
) -> Result<Vec<PhysicalRegion3>, EvaluationStop> {
    control.checkpoint("sdf3-region-transfer")?;
    let mut result = Vec::with_capacity(parents.len());
    for &parent in parents {
        control.checkpoint("sdf3-region-transfer")?;
        result.push(*source.get(parent).ok_or_else(|| failure("sdf3-region-parent"))?);
    }
    control.checkpoint("sdf3-region-transfer-publish")?;
    Ok(result)
}

#[cfg(test)]
mod tests;
