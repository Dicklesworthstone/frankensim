//! Retained owned motion across prescribed natural-reference changes.
//!
//! The native membrane transaction remains the sole owner of original-source,
//! topology, inventory, constitutive compatibility and work admission. This
//! adapter owns its accepted target bindings and buffers so a subsequent
//! temperature endpoint can outlive the caller's temporary prepared laws.
//! Reference input and endpoint defects retain their native signed meanings;
//! neither is converted into heat or silently included in a physical timestep.

use std::sync::Arc;

use super::{AcceptedMembrane, OwnedFiniteViscoMembraneMotion, OwnedLaws};
use crate::membrane::MembranePatch;
use crate::membrane::visco::{
    FiniteViscoMembraneError, FiniteViscoMembraneReferenceEnergy, FiniteViscoStepControls,
    ViscoelasticLawField, poll,
};

/// Complete owned reference proposal with the unchanged native energy receipt.
/// Local transverse equilibrium is restored; subsequent native advancement owns
/// the global mechanical response to any changed material forces.
#[derive(Debug)]
pub struct OwnedFiniteViscoMembraneMotionReferenceChange {
    pub motion: OwnedFiniteViscoMembraneMotion,
    pub energy: FiniteViscoMembraneReferenceEnergy,
    /// Native reference events supply no impulse and retain every nodal momentum.
    pub impulse_kg_m_s: [f64; 3],
    pub kinetic_energy_change_j: f64,
}

impl OwnedFiniteViscoMembraneMotion {
    /// Apply the original native reference transaction and own its admitted
    /// successor without retaining any borrow of the supplied material field.
    ///
    /// The target patch must retain original coordinates, connectivity, grain,
    /// source/model, thickness, strain envelope, corner inventory and selected
    /// full/exact-volume treatment. Its natural stretches must match the target
    /// laws. The native owner admits reference/kinetic changes, restores local
    /// transverse plane stress and gates every original work account.
    ///
    /// The supplied patch Arc is retained without rebuilding its geometry or
    /// mass. Target laws are copied only after native admission; uniform fields
    /// remain uniform. Actual displacement and velocity keep their original
    /// shared buffers; native material histories preserve the existing viscous
    /// metrics, accumulated loss and clock. No virgin history is constructed.
    ///
    /// Signed reference input, exact-volume constraint actuation and numerical
    /// defects are returned unchanged. This event supplies no heat, impulse or
    /// elapsed time and does not close a temperature-dependent caloric balance.
    /// Any refusal, including the final cancellation check, leaves `self` intact.
    pub fn change_reference(
        &self,
        target_patch: Arc<MembranePatch>,
        target_material: ViscoelasticLawField<'_>,
        controls: FiniteViscoStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OwnedFiniteViscoMembraneMotionReferenceChange, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        let native = self.as_native();
        let (accepted, energy, impulse_kg_m_s, kinetic_energy_change_j) = {
            let change = native.change_reference(
                &target_patch,
                target_material,
                controls,
                &mut cancelled,
            )?;
            (
                AcceptedMembrane::from_native(change.motion.membrane()),
                change.energy,
                change.impulse_kg_m_s,
                change.kinetic_energy_change_j,
            )
        };
        let motion = Self {
            patch: target_patch,
            laws: Arc::new(OwnedLaws::from_native(target_material)),
            accepted,
            // Native admission retained the complete original reference mass
            // and velocity. Share that exact buffer instead of its native copy.
            velocity_mps: Arc::clone(&self.velocity_mps),
            mass_image: self.mass_image,
        };
        poll(&mut cancelled)?;
        Ok(OwnedFiniteViscoMembraneMotionReferenceChange {
            motion,
            energy,
            impulse_kg_m_s,
            kinetic_energy_change_j,
        })
    }
}
