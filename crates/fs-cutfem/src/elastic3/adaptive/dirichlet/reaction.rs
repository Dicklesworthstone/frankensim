//! Reconstruct the physical trace and apply the same constraint transpose.
use super::*;
use crate::elastic3::dirichlet::EmbeddedReaction3;

impl AdaptiveElasticity3 {
    /// Evaluate the same Nitsche numerical reaction on a locally refined grid.
    /// Reconstruct u with T, then reduce its state derivative by T^T. Direct
    /// cell-scale derivatives keep their original active-cell ordering.
    /// No solve/model mutation or transfer of an already evaluated reaction.
    pub fn embedded_reaction(&self, u: &[f64], prescribed: Option<&PrescribedMotion3<'_>>,
        mode: &PrescribedMotion3<'_>, mut checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<EmbeddedReaction3, ElasticityError3> {
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        let physical = self.physical_displacements(u)?;
        let mut result = self.raw.embedded_reaction(&physical, prescribed, mode, &mut checkpoint)?;
        let mut gradient = vec![0.0; self.n()];
        self.restrict(&result.displacement_gradient, &mut gradient);
        if gradient.iter().any(|v| !v.is_finite()) {
            return Err(ElasticityError3::Invalid("reduced reaction gradient overflow"));
        }
        result.displacement_gradient = gradient;
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        Ok(result)
    }
}
