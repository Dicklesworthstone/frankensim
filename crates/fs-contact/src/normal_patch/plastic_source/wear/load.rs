//! Archard recession while an external normal force is maintained.
//!
//! The original accepted force must already match the requested maintained
//! load. Each native trial combines the proposed recession with the original
//! plastic history and solves the approach directly. There is no fictitious
//! held-approach indentation step before restoring the load.

use fs_exec::Cx;

use crate::normal_patch::{
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt, FiniteGapPlasticReceipt,
};

use super::{
    SourcePlasticArchardWear, SourcePlasticLoadError, SourcePlasticWearAdvance,
    SourcePlasticWearError, SourcePlasticWearState, checkpoint,
};

/// Normal loading with both original wear profiles and inventories retained.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearNormalLoadAdvance {
    /// Complete successor contact/wear state.
    pub candidate: SourcePlasticWearState,
    /// Original pressure, indentation, elastic storage, and plastic work.
    pub receipt: FiniteGapPlasticReceipt,
    /// Native force-control residual and approach search diagnostics.
    pub load: FiniteGapPlasticLoadReceipt,
}

/// Material removal and the equilibrium needed to retain the applied load.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearLoadAdvance {
    /// Original per-body Archard update and native energy accounts.
    ///
    /// Its continuation separates normal work `F * delta(approach)` from
    /// recession parameter work. Neither account is assigned to friction heat.
    pub wear: SourcePlasticWearAdvance,
    /// Prescribed-load equilibrium residual and approach search diagnostics.
    pub load: FiniteGapPlasticLoadReceipt,
}

impl SourcePlasticArchardWear {
    /// Change prescribed load while retaining removal and indentation history.
    ///
    /// The response is an endpoint equilibrium, not an unresolved loading-path
    /// work integral. Candidate contact still sees the remaining homogeneous
    /// material depth after all previously accepted removal.
    pub fn advance_normal_at_load(
        &self,
        before: &SourcePlasticWearState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearNormalLoadAdvance, SourcePlasticWearError> {
        self.admit(before)?;
        checkpoint(cx)?;
        let next = self
            .contact
            .advance_at_load(&before.contact, target_force_n, controls, cx)?;
        self.admit_receded_profile(&next.candidate, &before.recession_m, cx)?;
        Ok(SourcePlasticWearNormalLoadAdvance {
            candidate: SourcePlasticWearState {
                owner: self.identity,
                contact: next.candidate,
                volumes: before.volumes,
                recession_m: before.recession_m.clone(),
            },
            receipt: next.receipt,
            load: next.load,
        })
    }

    /// Integrate one beginning-pressure wear step at maintained normal force.
    ///
    /// Local slip remains unsigned irreversible relative sliding on the original
    /// grid. Each original body keeps its own coefficient and hardness. The
    /// accepted predecessor must already match the requested load; use
    /// `advance_normal_at_load` explicitly when changing that load.
    ///
    /// Recession and rigid approach evolve together in the native load solve.
    /// Physical limits and all per-body inventory checks apply atomically.
    /// Refine sliding steps when the local pressure distribution changes.
    pub fn advance_wear_at_load(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearLoadAdvance, SourcePlasticWearError> {
        self.admit(before)?;
        checkpoint(cx)?;
        self.contact.admit_load_controls(controls)?;
        let (wear, load) =
            self.advance_wear_with_recession(before, local_unsigned_slip_m, cx, |increment| {
                let response = self
                    .contact
                    .model()
                    .recede_at_load(
                        &before.contact.contact,
                        increment,
                        target_force_n,
                        controls,
                        cx,
                    )
                    .map_err(SourcePlasticLoadError::from)?;
                Ok((response.continuation, response.load))
            })?;
        Ok(SourcePlasticWearLoadAdvance { wear, load })
    }
}
