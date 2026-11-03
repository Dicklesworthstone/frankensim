//! Simultaneous native normal impacts with original-coordinate support work.

use super::*;
use crate::dynamics::impact::{
    self, SolidContactImpulseReport, SolidNormalContact, SolidNormalImpactControls,
    SolidNormalImpactError,
};

/// Contact-law or original supported-motion refusal. Neither publishes state.
#[derive(Debug)]
pub enum SupportedSolidNormalImpactError {
    Contact(SolidNormalImpactError),
    Solid(SupportedSolidError),
}

impl core::fmt::Display for SupportedSolidNormalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "supported solid normal impact refused: {self:?}")
    }
}
impl core::error::Error for SupportedSolidNormalImpactError {}
impl From<SolidNormalImpactError> for SupportedSolidNormalImpactError {
    fn from(error: SolidNormalImpactError) -> Self {
        Self::Contact(error)
    }
}
impl From<SupportedSolidError> for SupportedSolidNormalImpactError {
    fn from(error: SupportedSolidError) -> Self {
        Self::Solid(error)
    }
}

/// One accepted contact event at the existing trajectory time. The contact
/// law's determinacy classification is retained: redundant contact rows need
/// not determine individual impulses or their applied/support work partition.
/// The original represented momentum and combined physical work still close.
#[derive(Debug)]
pub struct SupportedSolidNormalImpactReport {
    /// Shared native contact law and actual endpoint/relative-work readback.
    pub contact: SolidContactImpulseReport,
    /// Original full-mass momentum change, support reactions and separate work.
    pub impulse: SupportedSolidImpulse,
}

type ImpactResult<T> = core::result::Result<T, SupportedSolidNormalImpactError>;

impl SupportedSpectralSolid {
    /// Resolve a declared set of simultaneous touching normal contacts through
    /// the original `fs-contact` law, then continue this supported trajectory.
    /// Rows use `represented_dofs()` order, INCLUDING prescribed coordinates.
    /// Positive row velocity opens the gap; each obstacle normal velocity is
    /// supplied explicitly by the contact. Geometry/event discovery remains
    /// caller-owned and no friction or contact-loss-to-heat law is inferred.
    ///
    /// The supports keep their accepted instantaneous velocities. The contact
    /// compliance is J_f M_ff^-1 J_f^T using every original free spectral mode;
    /// constrained rows have zero incremental mobility but retain their actual
    /// pre/post velocity in contact kinematics and work. The ORIGINAL full mass
    /// returns M_cf*Delta v_f - J_c as support reaction, including cross mass.
    /// Thus a passive relative impact can receive real work from a moving
    /// support without inventing kinetic energy or assigning it to heat.
    ///
    /// One native supported impulse candidate borrows all prepared matrices.
    /// Contact speeds, obstacle work, relative work, full momentum, support work
    /// and cancellation are admitted before either supported/free state commits.
    /// Positions, source preparation, support acceleration and clocks are fixed.
    pub fn apply_normal_impacts(
        &mut self,
        contacts: &[SolidNormalContact<'_>],
        controls: SolidNormalImpactControls,
        cx: &Cx<'_>,
    ) -> ImpactResult<SupportedSolidNormalImpactReport> {
        poll(cx)?;
        let proposal = impact::propose(
            &self.velocity,
            contacts,
            controls,
            |impulse, cx| self.momentum_response(impulse, cx),
            cx,
        )?;
        let mut support_velocity = zeros(self.support_rows.len())?;
        for (index, &row) in self.support_rows.iter().enumerate() {
            if index % 256 == 0 {
                poll(cx)?;
            }
            support_velocity[index] = self.velocity[row];
        }
        let trial = self.trial_momentum_impulse(
            &proposal.impulse_kg_m_s,
            &support_velocity,
            controls.maximum_impulse_error_kg_m_s,
            cx,
        )?;
        let contact = impact::audit(
            proposal,
            contacts,
            &trial.velocity_mps,
            trial.energy.applied_work_j,
            controls,
            cx,
        )?;
        poll(cx)?;
        let impulse = self.accept_momentum_impulse(trial);
        Ok(SupportedSolidNormalImpactReport { contact, impulse })
    }
}
