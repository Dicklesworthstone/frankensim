//! Prescribed normal force under the original material and geometry envelope.
//!
//! The native solve changes rigid approach to meet the supplied force. Every
//! trial starts from the original accepted indentation history; trial dents are
//! never committed. This quasistatic boundary condition supplies no actuator
//! dynamics or work integral for an unresolved change of applied force.

use fs_blake3::ContentHash;
use fs_exec::Cx;

use crate::normal_patch::{
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError, FiniteGapPlasticLoadReceipt,
    FiniteGapPlasticReceipt,
};

use super::{SourcePlasticContact, SourcePlasticContactError, SourcePlasticContactState};

/// Original source/applicability failure or native prescribed-load refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SourcePlasticLoadError {
    /// The original source owner or physical envelope refused the request.
    Contact(SourcePlasticContactError),
    /// Load is unattainable, ambiguous, unresolved, or otherwise refused.
    Load(FiniteGapPlasticLoadError),
}

impl From<SourcePlasticContactError> for SourcePlasticLoadError {
    fn from(error: SourcePlasticContactError) -> Self {
        Self::Contact(error)
    }
}

impl From<FiniteGapPlasticLoadError> for SourcePlasticLoadError {
    fn from(error: FiniteGapPlasticLoadError) -> Self {
        Self::Load(error)
    }
}

impl core::fmt::Display for SourcePlasticLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source hardness-limited load control: {self:?}")
    }
}

impl std::error::Error for SourcePlasticLoadError {}

/// Original-source candidate and the actual force/approach solve diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticContactLoadAdvance {
    /// Complete original source owner, including selected hardness tests.
    pub source_identity: ContentHash,
    /// Proposed successor; the accepted predecessor remains unchanged.
    pub candidate: SourcePlasticContactState,
    /// Native pressure, elastic storage, and incremental plastic work.
    pub receipt: FiniteGapPlasticReceipt,
    /// Requested/achieved load, residual, approach bracket and work count.
    pub load: FiniteGapPlasticLoadReceipt,
}

impl SourcePlasticContact {
    /// Solve for approach at a prescribed compressive normal force.
    ///
    /// The entire declared approach search must fit the original material's
    /// approach envelope. The accepted candidate also passes the original
    /// permanent-depth and half-space-footprint checks. A nonunique perfectly
    /// plastic load cannot choose an arbitrary additional indentation.
    ///
    /// No finite-load-ramp work is inferred from its two endpoints. Zero load
    /// follows the native explicit unloading convention, not free-body dynamics.
    pub fn advance_at_load(
        &self,
        before: &SourcePlasticContactState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticContactLoadAdvance, SourcePlasticLoadError> {
        if before.owner != self.identity {
            return Err(SourcePlasticContactError::Invalid(
                "load control requires the original source-bound history",
            )
            .into());
        }
        self.admit_load_controls(controls)?;
        let next = self
            .model
            .advance_at_load(&before.contact, target_force_n, controls, cx)?;
        self.admit_candidate(&next.contact.candidate, cx)?;
        Ok(SourcePlasticContactLoadAdvance {
            source_identity: self.identity,
            candidate: SourcePlasticContactState {
                owner: self.identity,
                contact: next.contact.candidate,
            },
            receipt: next.contact.receipt,
            load: next.load,
        })
    }

    pub(super) fn admit_load_controls(
        &self,
        controls: FiniteGapPlasticLoadControls,
    ) -> Result<(), SourcePlasticContactError> {
        super::within(
            "load-control maximum approach [m]",
            controls.maximum_approach_m,
            self.envelope.maximum_approach_m,
        )
    }
}
