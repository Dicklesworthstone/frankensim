//! A finite-stiffness loading device acting on original-source plastic contact.
//!
//! The device supplies `F = K * (position - approach)` with finite positive K.
//! Its compliance determines the approach even on a perfectly plastic contact
//! plateau. Material hardness, original test context, physical applicability,
//! and irreversible history remain owned by the existing contact law.

use fs_blake3::ContentHash;
use fs_exec::Cx;

use crate::normal_patch::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveError, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticReceipt,
};

use super::{SourcePlasticContact, SourcePlasticContactError, SourcePlasticContactState};

/// Original source/applicability refusal or unresolved spring/contact balance.
#[derive(Clone, Debug, PartialEq)]
pub enum SourcePlasticDriveError {
    /// Original source owner or physical envelope refused the candidate.
    Contact(SourcePlasticContactError),
    /// Native finite-stiffness equilibrium refused the request.
    Drive(FiniteGapPlasticDriveError),
}
impl From<SourcePlasticContactError> for SourcePlasticDriveError {
    fn from(error: SourcePlasticContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<FiniteGapPlasticDriveError> for SourcePlasticDriveError {
    fn from(error: FiniteGapPlasticDriveError) -> Self {
        Self::Drive(error)
    }
}
impl core::fmt::Display for SourcePlasticDriveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source hardness-limited finite-stiffness drive: {self:?}"
        )
    }
}
impl std::error::Error for SourcePlasticDriveError {}

/// Original-source history and the resolved loading-device equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticContactDriveAdvance {
    /// Complete original material/interface/test and applicability owner.
    pub source_identity: ContentHash,
    /// Immutable candidate, accepted only by the enclosing physical owner.
    pub candidate: SourcePlasticContactState,
    /// Native pressure field, stored contact energy, and plastic work.
    pub receipt: FiniteGapPlasticReceipt,
    /// Spring force/storage, equilibrium residual, bracket, and work budget.
    pub drive: FiniteGapPlasticDriveReceipt,
}

impl SourcePlasticContact {
    /// Solve contact and an explicitly supplied finite-stiffness device together.
    ///
    /// Every trial uses the accepted original indentation and recession history.
    /// The complete approach search must fit the original physical envelope;
    /// the accepted footprint and permanent depth must also remain admissible.
    /// Moving the device specifies a new endpoint, not its intervening work path.
    pub fn advance_with_drive(
        &self,
        before: &SourcePlasticContactState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticContactDriveAdvance, SourcePlasticDriveError> {
        if before.owner != self.identity {
            return Err(SourcePlasticContactError::Invalid(
                "finite-stiffness drive requires the original source-bound history",
            )
            .into());
        }
        self.admit_load_controls(controls)?;
        let next = self
            .model
            .advance_with_drive(&before.contact, drive, controls, cx)?;
        self.admit_candidate(&next.contact.candidate, cx)?;
        Ok(SourcePlasticContactDriveAdvance {
            source_identity: self.identity,
            candidate: SourcePlasticContactState {
                owner: self.identity,
                contact: next.contact.candidate,
            },
            receipt: next.contact.receipt,
            drive: next.drive,
        })
    }
}
