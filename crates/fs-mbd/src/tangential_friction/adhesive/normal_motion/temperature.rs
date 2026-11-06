//! Original DMT thermal continuation at the captured actual native geometry.
//!
//! The prepared pair keeps its original native bodies, poses, mass properties,
//! geometry markers, reaction plane, time step and complete normal/shear
//! mobility. Only the original normal/shear sources and their authenticated
//! history change. The source supplies signed parameter work at fixed position;
//! there is no native impulse, slip work or inferred heat during this update.
//! Source receipts and successor shear remain caller-owned immutable borrows.

use super::{DmtCoupledState, PreparedRigidDmtPair, RigidDmtContact, RigidDmtError, finite, poll};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::SourceDmtThermalRefresh;
use fs_contact::normal_patch::dmt::{DmtResponse, coupled::DmtCoupledRegime};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Signed original-source change at unchanged body momenta and geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDmtSourceRefreshReport {
    pub normal_before: DmtResponse,
    pub normal_after: DmtResponse,
    /// Accepted normal/shear coordinates, unchanged by source continuation.
    pub position_m: [f64; 3],
    /// Original occupied branch, including Released in an attractive gap.
    pub regime: DmtCoupledRegime,
    /// Accepted source approach minus actual captured marker approach [m].
    pub geometry_residual_m: f64,
    /// Normal potential change [J]; the native kinetic change is exactly zero.
    pub normal_storage_change_j: f64,
    /// Signed supply of the explicit source-temperature driver [J].
    pub parameter_work_j: f64,
    /// Normal potential change minus its original parameter work [J].
    pub energy_residual_j: f64,
    /// Original contact-sized absolute gate, independent of body kinetic energy.
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum RigidDmtSourceRefreshError {
    Mechanical(RigidDmtError),
    Energy { residual_j: f64, tolerance_j: f64 },
}
impl From<RigidDmtError> for RigidDmtSourceRefreshError {
    fn from(error: RigidDmtError) -> Self {
        Self::Mechanical(error)
    }
}
impl core::fmt::Display for RigidDmtSourceRefreshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native DMT source-temperature change refused: {self:?}")
    }
}
impl std::error::Error for RigidDmtSourceRefreshError {}

impl<'a> PreparedRigidDmtPair<'a> {
    /// Stage the original joint at prescribed temperatures and fixed actual
    /// normal/shear coordinates. The source receipt must use held indentation
    /// at this captured geometry; a held-load change cannot move native bodies.
    /// `next_shear` must be the original selected junction claims' successor.
    ///
    /// The original source transfer retains branch occupancy, rechecks the
    /// complete mechanical/source envelope and preserves every compliance entry.
    /// Native mass/inertia and body geometry do not change with these contact
    /// source temperatures. The caller owns thermal expansion and heat transport.
    /// No accepted owner changes if any source, geometry, work or cancel gate
    /// refuses. After a mechanical step, consume its reported pose travel and
    /// prepare the updated geometry before refreshing that advanced history.
    pub fn with_temperature_refresh(
        &self,
        previous: DmtCoupledState,
        refreshed: &'a SourceDmtThermalRefresh,
        next_shear: &'a SourceAdhesiveShear,
        gate: &CancelGate,
    ) -> Result<(Self, DmtCoupledState, RigidDmtSourceRefreshReport), RigidDmtSourceRefreshError>
    {
        poll(gate).map_err(RigidDmtError::from)?;
        let position = previous.position_m();
        let geometry_residual = finite(position[0] - self.approach_m)?;
        if geometry_residual.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidDmtError::GeometryMismatch {
                residual_m: geometry_residual,
            }
            .into());
        }
        let (load_line, state) = self
            .load_line
            .with_temperature_refresh(previous, refreshed, next_shear, gate)
            .map_err(RigidDmtError::from)?;
        let before = *refreshed.before();
        let after = *refreshed.response();
        let change =
            Dd::from_f64(after.reversible_energy_j) - Dd::from_f64(before.reversible_energy_j);
        let work = finite(refreshed.parameter_work_j())?;
        let residual = finite((change - Dd::from_f64(work)).to_f64())?;
        let tolerance = self.controls.energy_tolerance_j;
        // Check the small source account directly. Subtracting complete body
        // energies would erase this work under a large common native motion.
        for error in [
            residual,
            finite((change - Dd::from_f64(refreshed.potential_change_j())).to_f64())?,
        ] {
            if error.abs() > tolerance {
                return Err(RigidDmtSourceRefreshError::Energy {
                    residual_j: error,
                    tolerance_j: tolerance,
                });
            }
        }
        let candidate = Self {
            contact: RigidDmtContact {
                normal: refreshed.source(),
                shear: next_shear,
                ..self.contact
            },
            before: self.before,
            properties: self.properties,
            reaction_arms: self.reaction_arms,
            axes: self.axes,
            markers_before: self.markers_before,
            reaction_before: self.reaction_before,
            approach_m: self.approach_m,
            mobility: self.mobility,
            load_line,
            dt_s: self.dt_s,
            controls: self.controls,
        };
        let report = RigidDmtSourceRefreshReport {
            normal_before: before,
            normal_after: after,
            position_m: position,
            regime: state.regime(),
            geometry_residual_m: geometry_residual,
            normal_storage_change_j: finite(change.to_f64())?,
            parameter_work_j: work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        };
        poll(gate).map_err(RigidDmtError::from)?;
        Ok((candidate, state, report))
    }

    /// Publish the same captured native pair with its original thermal sources
    /// and transferred history. The surrounding body array is never mutated.
    pub fn apply_temperature_refresh(
        &mut self,
        history: &mut DmtCoupledState,
        refreshed: &'a SourceDmtThermalRefresh,
        next_shear: &'a SourceAdhesiveShear,
        gate: &CancelGate,
    ) -> Result<RigidDmtSourceRefreshReport, RigidDmtSourceRefreshError> {
        self.apply_temperature_refresh_with_admission(
            history,
            refreshed,
            next_shear,
            gate,
            |report, _, _| Ok(*report),
        )
    }

    /// Stage another receiving owner before publishing the pair and history.
    /// The callback returns an infallibly installable candidate and must not
    /// mutate accepted external state. Receiver refusal or late cancellation
    /// retains the previous sources, branch, native snapshot and solve workspace.
    pub fn apply_temperature_refresh_with_admission<R, E>(
        &mut self,
        history: &mut DmtCoupledState,
        refreshed: &'a SourceDmtThermalRefresh,
        next_shear: &'a SourceAdhesiveShear,
        gate: &CancelGate,
        admit: impl FnOnce(&RigidDmtSourceRefreshReport, &Self, DmtCoupledState) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<RigidDmtSourceRefreshError>,
    {
        let (candidate, state, report) = self
            .with_temperature_refresh(*history, refreshed, next_shear, gate)
            .map_err(E::from)?;
        let receiving = admit(&report, &candidate, state)?;
        poll(gate)
            .map_err(RigidDmtError::from)
            .map_err(RigidDmtSourceRefreshError::from)
            .map_err(E::from)?;
        *self = candidate;
        *history = state;
        Ok(receiving)
    }
}
