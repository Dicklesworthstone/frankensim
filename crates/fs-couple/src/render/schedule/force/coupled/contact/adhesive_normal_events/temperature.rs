//! Prescribed-temperature continuation at an accepted native mechanical state.
//!
//! The original source changes contact storage at fixed geometry. Its signed
//! parameter work belongs to the prescribed driver, not to dissipative heat.
//! No modal coordinate, momentum, support history or clock changes in this
//! instantaneous split. Later intervals use the successor's actual force law.

use super::{
    JkrNormalEventModalError, JkrNormalEventModalSystem, check_energy, compatible, finite, limit,
    poll,
};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrResponse,
    motion::temperature::{JkrNormalTemperatureContinuation, JkrNormalTemperatureReport},
    temperature::SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;

/// Original coefficient supply and its independent native-geometry check.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalTemperatureModalReport {
    pub sample: u64,
    pub source: JkrNormalTemperatureReport,
    pub native_before: JkrResponse,
    pub native_after: JkrResponse,
    pub native_storage_change_j: f64,
    /// Actual native contact-storage change minus the original parameter supply.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

impl<'a> JkrNormalEventModalSystem<'a> {
    /// Continue the original source at held approach, without resetting motion.
    ///
    /// The caller retains the immutable original refresh for the lifetime of
    /// the continued owner. Its signed `parameter_work_j` must be supplied or
    /// absorbed by the prescribed temperature driver. This method neither
    /// derives an interface temperature nor withdraws heat from caloric cells.
    /// The last frame continues to describe its already accepted interval.
    pub fn continue_at_temperature(
        &mut self,
        refresh: &'a SourceJkrThermalRefresh,
        gate: &CancelGate,
    ) -> Result<JkrNormalTemperatureModalReport, JkrNormalEventModalError> {
        self.continue_at_temperature_with_admission(refresh, gate, |_, report| {
            Ok::<_, JkrNormalEventModalError>(*report)
        })
    }

    /// Stage receiving history before publishing the new mechanical source.
    /// A refusal or cancellation leaves every accepted mechanical state intact.
    /// The callback must return an uncommitted candidate, as in step admission.
    pub fn continue_at_temperature_with_admission<R, E>(
        &mut self,
        refresh: &'a SourceJkrThermalRefresh,
        gate: &CancelGate,
        admit: impl FnOnce(
            &JkrNormalTemperatureContinuation<'a>,
            &JkrNormalTemperatureModalReport,
        ) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<JkrNormalEventModalError>,
    {
        let (continuation, report) = self.stage_temperature(refresh, gate).map_err(E::from)?;
        let receiving = admit(&continuation, &report)?;
        poll(Some(gate))
            .map_err(JkrNormalEventModalError::from)
            .map_err(E::from)?;
        let (motion, state, _) = continuation.into_parts();
        self.motion = motion;
        self.state = state;
        Ok(receiving)
    }

    fn stage_temperature(
        &self,
        refresh: &'a SourceJkrThermalRefresh,
        gate: &CancelGate,
    ) -> Result<
        (
            JkrNormalTemperatureContinuation<'a>,
            JkrNormalTemperatureModalReport,
        ),
        JkrNormalEventModalError,
    > {
        poll(Some(gate))?;
        let continuation = self
            .motion
            .continue_at_temperature(self.state, refresh, gate)?;
        let source = *continuation.report();
        let before = self.native_contact(self.components(), self.state.contact_state().branch())?;
        let delta = before.indentation_m;
        compatible(
            "temperature refresh retains native approach",
            delta,
            self.state.contact_state().indentation_m(),
            self.config.coordinate_absolute_tolerance_m,
            self.config.coordinate_relative_tolerance,
        )?;
        let successor = continuation.motion().normal();
        let after = match before.state().branch() {
            JkrBranch::Attached => successor.evaluate_attached(delta)?,
            JkrBranch::Detached => {
                let state = successor.detached_state(delta)?;
                continuation
                    .motion()
                    .response(continuation.motion().initialize(state, 0.0)?)?
            }
        };
        limit(
            "native refreshed JKR contact force magnitude",
            after.normal_force_n.abs(),
            self.config.motion.maximum_force_n,
        )?;
        // Network storage is exactly unchanged. Gate the small contact supply
        // independently instead of subtracting two large network-energy datums.
        let change = finite(after.reversible_energy_j - before.reversible_energy_j)?;
        let residual = finite(change - source.parameter_work_j)?;
        let tolerance = finite(
            self.config.motion.energy_absolute_tolerance_j
                + self.config.motion.energy_relative_tolerance
                    * before
                        .reversible_energy_j
                        .abs()
                        .max(after.reversible_energy_j.abs())
                        .max(source.parameter_work_j.abs()),
        )?;
        check_energy(residual, tolerance)?;
        limit(
            "native refreshed JKR total energy magnitude",
            finite(self.network.total_energy_j()? + after.reversible_energy_j)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok((
            continuation,
            JkrNormalTemperatureModalReport {
                sample: self.samples_rendered(),
                source,
                native_before: before,
                native_after: after,
                native_storage_change_j: change,
                energy_residual_j: residual,
                energy_tolerance_j: tolerance,
            },
        ))
    }
}
