//! Prescribed natural-strain changes at fixed body deformation and memory.
//!
//! With `A = F G^-1`, the stored energy is `Psi = det(G) W(A, Cv)` per
//! original reference volume. Changing prescribed `G` at fixed total `F` and
//! fixed natural-coordinate viscous metrics supplies the signed energy
//! `Psi_after - Psi_before`. No mechanical displacement occurs, no physical
//! time advances and no Maxwell relaxation heat is generated. Subsequent
//! positive-duration steps use the new reference and the SAME retained metrics.
//!
//! This is a physical eigenstrain change in fixed material axes, not a change
//! of coordinates. It does not transport metrics/fibres to a renamed frame,
//! evolve a growth law or define thermal entropy/caloric exchange. Source
//! owners must retain their original thermal datum and material frame. Changing
//! moduli, branch energies/order or the elastic domain is a separate model.

use super::{
    Hyperelastic, MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic, ViscoRemainderState, ViscoState,
    parameters,
};
use crate::MaterialError;

/// Immutable proposal for a prescribed-reference change of a full history.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoReferenceChange {
    /// Same total F, viscous metrics, mechanical work and accumulated loss,
    /// with the target reference binding, storage and reference-work account.
    pub state: ViscoState,
    /// Instantaneous Piola stress after the change, conjugate to total F.
    pub first_piola_pa: [f64; 9],
    /// Signed energy supplied by this reference change, per original volume.
    /// Negative values release stored energy; neither sign is viscous heat.
    pub reference_work_pa: f64,
}

/// Immutable proposal with an explicitly bulk-free reference-energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoRemainderReferenceChange {
    /// Retained metrics and total F with the new remainder storage/binding.
    pub state: ViscoRemainderState,
    /// Instantaneous remainder Piola stress; no pressure stress is included.
    pub first_piola_pa: [f64; 9],
    /// Signed remainder-energy exchange. A mixed owner separately accounts
    /// for pressure-reconstructed bulk storage and re-admits its volume constraint.
    pub reference_work_pa: f64,
}

impl ViscoHyperelastic {
    fn validate_reference_change_from(&self, previous: &Self) -> Result<(), MaterialError> {
        if self.equilibrium != previous.equilibrium
            || self.branches.len() != previous.branches.len()
            || self.branches.len() > MAX_FINITE_STRAIN_BRANCHES
            || self
                .branches
                .iter()
                .zip(&previous.branches)
                .any(|(new, old)| new.energy != old.energy)
        {
            return Err(parameters(
                "reference changes require identical equilibrium, elastic domain and ordered branch energies"
                    .into(),
            ));
        }
        // Public material parameters may have changed since construction.
        Hyperelastic::new(
            self.equilibrium.model.clone(),
            self.equilibrium.stretch_limit,
        )?;
        for (new, old) in self.branches.iter().zip(&previous.branches) {
            new.energy.validate()?;
            if ![new.tau_s, old.tau_s]
                .iter()
                .all(|tau| tau.is_finite() && *tau > 0.0)
            {
                return Err(parameters(
                    "reference changes require finite positive relaxation times".into(),
                ));
            }
        }
        Ok(())
    }

    /// Change a prescribed natural configuration at fixed total deformation
    /// and fixed viscous metrics. Call this on the TARGET law, passing the
    /// previous law that owns `accepted`. Neither input is mutated.
    /// `accepted` must be an unmodified output of that owner's initialization,
    /// step or reference-change methods. Its physical storage ledger is trusted,
    /// as in ordinary stepping; finite public-field edits or imported histories
    /// are not certified by this transaction.
    ///
    /// Equilibrium parameters, stretch domain, branch energies and their order
    /// must agree exactly. Relaxation times may change, as for a temperature
    /// shift: no time passes during this event. Ordinary stepping resumes on
    /// this target law with the returned history. Its ledger obeys
    /// `mechanical work + reference work = change in storage + relaxation loss`.
    /// An unchanged G retains the complete accepted history exactly, even when
    /// kinetic times or Newton controls change.
    ///
    /// # Errors
    /// Refuses mismatched physical laws, an invalid old binding or metrics,
    /// nonfinite ledgers or negative accumulated relaxation loss,
    /// either elastic endpoint outside its domain, or unrepresentable storage,
    /// stress or signed reference-work accumulation. No candidate is published.
    pub fn change_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoState,
    ) -> Result<ViscoReferenceChange, MaterialError> {
        self.validate_reference_change_from(previous)?;
        previous.freeze(accepted, 0.0)?;
        previous.equilibrium.check_state(
            &previous
                .reference_deformation
                .elastic_deformation(&accepted.deformation_gradient),
        )?;
        let mut rebound = accepted.clone();
        rebound.reference_gradient = *self.reference_deformation.deformation_gradient();
        // The zero-duration evaluator supplies an energy increment, not a
        // mechanical increment here. Avoid an irrelevant intermediate sum
        // with the mechanical accumulator, which could otherwise overflow.
        rebound.work = 0.0;
        let endpoint = self.step_response(&rebound, &accepted.deformation_gradient, 0.0)?;
        let (state, reference_work_pa) = if self.reference_deformation.deformation_gradient()
            == previous.reference_deformation.deformation_gradient()
        {
            (accepted.clone(), 0.0)
        } else {
            let mut state = endpoint.state;
            state.work = accepted.work;
            state.dissipated = accepted.dissipated;
            state.reference_work = accepted.reference_work + endpoint.loading_work_pa;
            state.check_ledger()?;
            (state, endpoint.loading_work_pa)
        };
        Ok(ViscoReferenceChange {
            state,
            first_piola_pa: endpoint.first_piola_pa,
            reference_work_pa,
        })
    }

    /// Bulk-free counterpart of [`Self::change_reference_from`]. The physical
    /// remainder is evaluated directly; no potentially large bulk penalty is
    /// evaluated or subtracted. Returned reference work excludes bulk storage.
    /// The same unmodified-owner-history precondition applies to `accepted`.
    ///
    /// A spatial mixed owner must separately re-admit pressure and its volume
    /// constraint. In particular, changing det(G) at fixed F generally violates
    /// exact elastic incompressibility. This material proposal does not claim
    /// that such a body configuration is an admissible mixed equilibrium.
    pub fn change_remainder_reference_from(
        &self,
        previous: &Self,
        accepted: &ViscoRemainderState,
    ) -> Result<ViscoRemainderReferenceChange, MaterialError> {
        self.validate_reference_change_from(previous)?;
        previous.freeze_remainder(accepted, 0.0)?;
        let mut rebound = accepted.clone();
        rebound.reference_gradient = *self.reference_deformation.deformation_gradient();
        rebound.work = 0.0;
        let endpoint =
            self.step_remainder_response(&rebound, &accepted.deformation_gradient, 0.0)?;
        let (state, reference_work_pa) = if self.reference_deformation.deformation_gradient()
            == previous.reference_deformation.deformation_gradient()
        {
            (accepted.clone(), 0.0)
        } else {
            let mut state = endpoint.state;
            state.work = accepted.work;
            state.dissipated = accepted.dissipated;
            state.reference_work = accepted.reference_work + endpoint.loading_work_pa;
            state.check_ledger()?;
            (state, endpoint.loading_work_pa)
        };
        Ok(ViscoRemainderReferenceChange {
            state,
            first_piola_pa: endpoint.first_piola_pa,
            reference_work_pa,
        })
    }
}
