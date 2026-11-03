//! Simultaneous Coulomb impulses with original full-mass support reactions.

use super::*;
use crate::dynamics::frictional_impact::{
    self, SolidFrictionalContact, SolidFrictionalContactImpulseReport,
    SolidFrictionalImpactControls, SolidFrictionalImpactError,
};

/// A failed law or full supported-state audit publishes neither trajectory.
#[derive(Debug)]
pub enum SupportedSolidFrictionalImpactError {
    Contact(SolidFrictionalImpactError),
    Solid(SupportedSolidError),
}

impl core::fmt::Display for SupportedSolidFrictionalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "supported solid frictional impact refused: {self:?}")
    }
}
impl core::error::Error for SupportedSolidFrictionalImpactError {}
impl From<SolidFrictionalImpactError> for SupportedSolidFrictionalImpactError {
    fn from(error: SolidFrictionalImpactError) -> Self {
        Self::Contact(error)
    }
}
impl From<SupportedSolidError> for SupportedSolidFrictionalImpactError {
    fn from(error: SupportedSolidError) -> Self {
        Self::Solid(error)
    }
}

/// One event on the original trajectory. A redundant contact law retains its
/// determinacy qualification: local load and work partitions can depend on the
/// chosen impulse representative even when the complete motion is determined.
#[derive(Debug)]
pub struct SupportedSolidFrictionalImpactReport {
    pub contact: SolidFrictionalContactImpulseReport,
    pub impulse: SupportedSolidImpulse,
}

impl SupportedSpectralSolid {
    /// Resolve declared normal and two-tangent contacts simultaneously through
    /// the original non-associated Signorini--Coulomb law. Contact covectors use
    /// `represented_dofs()` order, including prescribed coordinates, in a
    /// caller-owned orthonormal contact frame. Positive normal velocity opens
    /// the gap. Every contact's opposing surface velocity remains explicit.
    ///
    /// Supports keep their accepted velocities during this event. The complete
    /// free spectrum supplies J_f M_ff^-1 J_f^T, including normal/tangent and
    /// inter-contact cross terms. Prescribed coordinates have zero incremental
    /// mobility but retain their actual velocity in slip and work. The original
    /// full mass then returns R_c = M_cf Delta-v_f - J_c; an impulse at a fixed
    /// coordinate is not silently discarded.
    ///
    /// Actual endpoint slip, normal law, friction cone, relative passivity,
    /// prescribed-surface work and original momentum/energy gates all precede
    /// the single acceptance of both supported and free trajectories. Support
    /// work is separate from applied contact work. Positions, preparation,
    /// history and clocks remain unchanged. This is an instantaneous event,
    /// not contact detection, sustained friction or a loss-to-heat model.
    pub fn apply_frictional_impacts(
        &mut self,
        contacts: &[SolidFrictionalContact<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> core::result::Result<
        SupportedSolidFrictionalImpactReport,
        SupportedSolidFrictionalImpactError,
    > {
        poll(cx)?;
        let proposal = frictional_impact::propose(
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
        let contact = frictional_impact::audit(
            proposal,
            contacts,
            &trial.velocity_mps,
            trial.energy.applied_work_j,
            controls,
            cx,
        )?;
        poll(cx)?;
        let impulse = self.accept_momentum_impulse(trial);
        Ok(SupportedSolidFrictionalImpactReport { contact, impulse })
    }
}
