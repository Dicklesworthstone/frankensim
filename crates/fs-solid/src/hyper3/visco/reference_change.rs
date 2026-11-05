//! Held-configuration material-reference changes with retained native memory.

use std::sync::Arc;

use fs_exec::Cx;

use super::{
    FiniteViscoTetError, FiniteViscoTetSolid, HyperTetError, ViscoelasticLawField, checkpoint,
    invalid,
};

/// Complete proposed material-reference change at unchanged displacement/time.
/// Dropping the result leaves the previous solid and every history intact.
#[derive(Debug)]
pub struct FiniteViscoTetReferenceChange<'a> {
    pub solid: FiniteViscoTetSolid<'a>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Signed external parameter work: positive adds physical stored energy.
    /// This is neither nodal mechanical work nor viscous dissipation.
    pub reference_work_j: f64,
}

impl<'a> FiniteViscoTetSolid<'a> {
    /// Change prescribed stress-free references while holding actual deformation.
    ///
    /// The material owner retains each accepted viscous metric in its declared
    /// natural coordinates. Only the prescribed reference and relaxation times
    /// may differ; changing equilibrium or branch energies is refused. The
    /// signed storage exchange is reported separately from mechanical work and
    /// dissipation. No time elapses and no history is reset or relaxed.
    ///
    /// Original geometry and displacement buffers are shared. Changed stresses
    /// may create a force imbalance; a subsequent native [`Self::advance`] owns
    /// the mechanical response. All elements must admit before a successor is
    /// returned, and its new accepted identity invalidates old prepared contact
    /// or split-motion responses. The previous solid remains usable on refusal.
    pub fn change_reference<'b>(
        &self,
        target_material: ViscoelasticLawField<'b>,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoTetReferenceChange<'b>, FiniteViscoTetError>
    where
        'a: 'b,
    {
        checkpoint(cx)?;
        if let ViscoelasticLawField::PerElement(laws) = target_material
            && laws.len() != self.elements.len()
        {
            return Err(HyperTetError::LawCountMismatch {
                expected: self.elements.len(),
                got: laws.len(),
            }
            .into());
        }
        let mut states = Vec::with_capacity(self.states.len());
        let mut first_piola_pa = Vec::with_capacity(self.states.len());
        let mut stored_after_j = 0.0;
        let mut reference_work_j = 0.0;
        for (element, geometry) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let change = target_material
                .at(element)
                .change_reference_from(self.material.at(element), &self.states[element])
                .map_err(|source| FiniteViscoTetError::Material { element, source })?;
            stored_after_j += geometry.volume * change.state.stored;
            reference_work_j += geometry.volume * change.reference_work_pa;
            first_piola_pa.push(change.first_piola_pa);
            states.push(change.state);
        }
        if !stored_after_j.is_finite() || !reference_work_j.is_finite() {
            return Err(invalid(
                "material-reference storage or work is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        Ok(FiniteViscoTetReferenceChange {
            solid: FiniteViscoTetSolid {
                nodes_m: self.nodes_m,
                tetrahedra: self.tetrahedra,
                material: target_material,
                elements: Arc::clone(&self.elements),
                budget: self.budget,
                displacement_m: Arc::clone(&self.displacement_m),
                states: states.into(),
                first_piola_pa: first_piola_pa.into(),
                stored_energy_j: stored_after_j,
                time_s: self.time_s,
                accepted_identity: Arc::new(()),
            },
            stored_before_j: self.stored_energy_j,
            stored_after_j,
            reference_work_j,
        })
    }
}
