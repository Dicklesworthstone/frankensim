//! Joint native JKR mechanics and event-timed caloric transport.
//!
//! The existing mechanical owner resolves normal motion and branch events;
//! the existing thermal timeline transports their explicitly declared heat.
//! Original support damping has its own heat/outside partition and finite-time
//! source history. It never inherits the adhesive snap partition implicitly.
//! Both accepted owners publish together. Material/contact coefficients stay
//! frozen over an interval. Prescribed temperature continuation at its boundary
//! retains all three histories and reports signed driver work; it does not
//! infer autonomous contact entropy or temperature feedback from caloric heat.

use super::temperature::JkrNormalTemperatureModalReport;
use super::{JkrNormalEventModalError, JkrNormalEventModalFrame, JkrNormalEventModalSystem};
use crate::friction_heat::adhesive::snap::{
    JkrSnapHeat, JkrSnapHeatError,
    timed::{JkrTimedSnapHeatControls, JkrTimedSnapHeatReport},
    viscous::{JkrViscousHeatReport, JkrViscousHeatRoute},
};
use fs_contact::normal_patch::adhesive::temperature::SourceJkrThermalRefresh;
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::HeatPartition;

#[derive(Debug)]
pub enum JkrThermalModalError {
    Mechanical(JkrNormalEventModalError),
    Thermal(JkrSnapHeatError),
    Invalid(&'static str),
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}
impl From<JkrNormalEventModalError> for JkrThermalModalError {
    fn from(value: JkrNormalEventModalError) -> Self {
        Self::Mechanical(value)
    }
}
impl From<JkrSnapHeatError> for JkrThermalModalError {
    fn from(value: JkrSnapHeatError) -> Self {
        Self::Thermal(value)
    }
}
impl core::fmt::Display for JkrThermalModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint native JKR thermal motion refused: {self:?}")
    }
}
impl std::error::Error for JkrThermalModalError {}

#[derive(Clone, Debug, PartialEq)]
pub struct JkrThermalModalFrame {
    pub mechanical: JkrNormalEventModalFrame,
    pub heat: JkrTimedSnapHeatReport,
    /// Separate support-loss destinations; adhesive snap work stays in `heat`.
    pub support_heat: Option<JkrViscousHeatReport>,
    /// Native boundary/source heat excluding both internal mechanical transfers.
    pub external_transport_heat_j: f64,
    /// Difference between source-regenerated heat debit and mechanical release.
    pub snap_transfer_residual_j: f64,
    pub support_transfer_residual_j: f64,
    /// Composed native mechanical and thermal residuals after cancelling the
    /// same original release. Native birth-relative thermal audits are retained
    /// rather than subtracting large absolute formation-enthalpy datums.
    pub combined_energy_residual_j: f64,
    pub combined_energy_tolerance_j: f64,
    /// Actual endpoint mechanical plus native thermal storage, excluding the
    /// declared outside energy destination. It is an observation, not the
    /// small-work conservation gate above.
    pub stored_energy_j: f64,
}

/// Own the actual body network, source-bound heat history and thermal timeline.
/// Read-only access prevents advancing a participating owner independently.
pub struct JkrThermalModalSystem<'s, 'm> {
    mechanical: JkrNormalEventModalSystem<'s>,
    heat: JkrSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    controls: JkrTimedSnapHeatControls,
    support_heat: Option<JkrViscousHeatRoute>,
    thermal_time_at_motion_zero_s: f64,
    last: Option<JkrThermalModalFrame>,
}

impl<'s, 'm> JkrThermalModalSystem<'s, 'm> {
    /// Couple admitted free or single-spring native mechanics to two original
    /// caloric cells and their declared heat/outside partition. The explicit
    /// clock origin maps mechanical time into the original thermal timeline;
    /// the supplied thermal state must already be at the mapped initial time.
    ///
    /// Initial snapshots define this owner's accounting window. Every later
    /// mechanical interval advances transport at its actual snap timestamps,
    /// with native scheduled thermal events applied before a coincident snap.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mechanical: JkrNormalEventModalSystem<'s>,
        thermal: DepositionSimulation<'m>,
        cells: [CellId; 2],
        caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
        partition: HeatPartition,
        outside_destination: Option<String>,
        thermal_time_at_motion_zero_s: f64,
        controls: JkrTimedSnapHeatControls,
        gate: &CancelGate,
    ) -> Result<Self, JkrThermalModalError> {
        Self::new_inner(
            mechanical,
            thermal,
            cells,
            caloric,
            partition,
            outside_destination,
            thermal_time_at_motion_zero_s,
            controls,
            None,
            gate,
        )
    }

    /// Couple original support damping with its explicitly declared heat shares
    /// and outside destination. The two pinned caloric cells receive smooth
    /// support heat through native source powers over sealed mechanical segments;
    /// adhesive releases retain their independent partition and snap timestamps.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_support_heat(
        mechanical: JkrNormalEventModalSystem<'s>,
        thermal: DepositionSimulation<'m>,
        cells: [CellId; 2],
        caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
        partition: HeatPartition,
        outside_destination: Option<String>,
        thermal_time_at_motion_zero_s: f64,
        controls: JkrTimedSnapHeatControls,
        support_heat: JkrViscousHeatRoute,
        gate: &CancelGate,
    ) -> Result<Self, JkrThermalModalError> {
        Self::new_inner(
            mechanical,
            thermal,
            cells,
            caloric,
            partition,
            outside_destination,
            thermal_time_at_motion_zero_s,
            controls,
            Some(support_heat),
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_inner(
        mechanical: JkrNormalEventModalSystem<'s>,
        thermal: DepositionSimulation<'m>,
        cells: [CellId; 2],
        caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
        partition: HeatPartition,
        outside_destination: Option<String>,
        thermal_time_at_motion_zero_s: f64,
        controls: JkrTimedSnapHeatControls,
        support_heat: Option<JkrViscousHeatRoute>,
        gate: &CancelGate,
    ) -> Result<Self, JkrThermalModalError> {
        checkpoint(gate).map_err(JkrSnapHeatError::from)?;
        if mechanical.support_damping_n_s_per_m() > 0.0 && support_heat.is_none() {
            return Err(JkrThermalModalError::Invalid(
                "damped support requires an explicit heat/outside destination",
            ));
        }
        let time = finite(thermal_time_at_motion_zero_s + mechanical.state().time_s())?;
        if !thermal_time_at_motion_zero_s.is_finite()
            || time != thermal.time_s()
            || controls.maximum_transport_substeps == 0
            || controls.heat.maximum_segments == 0
            || !controls.transport.max_dt_s.is_finite()
            || controls.transport.max_dt_s <= 0.0
            || [
                controls.heat.cell_energy_tolerance_j,
                controls.heat.thermal_energy_tolerance_j,
                controls.heat.transfer_energy_tolerance_j,
                controls.transport.energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(JkrThermalModalError::Invalid(
                "matching finite clocks and explicit heat/transport budgets required",
            ));
        }
        let heat = JkrSnapHeat::new(
            mechanical.normal(),
            mechanical.state().contact_state(),
            &thermal,
            cells,
            caloric,
            partition,
            outside_destination,
        )?;
        finite(mechanical.total_energy_j()? + thermal.energy_j())?;
        checkpoint(gate).map_err(JkrSnapHeatError::from)?;
        Ok(Self {
            mechanical,
            heat,
            thermal,
            controls,
            support_heat,
            thermal_time_at_motion_zero_s,
            last: None,
        })
    }

    pub const fn mechanical(&self) -> &JkrNormalEventModalSystem<'s> {
        &self.mechanical
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn heat_history(&self) -> &JkrSnapHeat<'s> {
        &self.heat
    }
    pub const fn controls(&self) -> JkrTimedSnapHeatControls {
        self.controls
    }
    pub const fn support_heat_route(&self) -> Option<&JkrViscousHeatRoute> {
        self.support_heat.as_ref()
    }
    pub const fn thermal_time_at_motion_zero_s(&self) -> f64 {
        self.thermal_time_at_motion_zero_s
    }
    pub const fn last_frame(&self) -> Option<&JkrThermalModalFrame> {
        self.last.as_ref()
    }

    pub fn stored_energy_j(&self) -> Result<f64, JkrThermalModalError> {
        finite(self.mechanical.total_energy_j()? + self.thermal.energy_j())
    }

    /// Atomically continue mechanics and snap history under prescribed source
    /// temperatures at the accepted boundary. The native thermal state, clock,
    /// caloric storage and accumulated heat ledger remain exactly unchanged.
    /// The returned signed parameter work is supplied by the prescribed driver
    /// and must not be passed through either dissipative heat partition.
    ///
    /// Independent surface and interface temperatures come from the original
    /// source refresh. Inferring them from finite cells requires an additional
    /// thermodynamic interface law; this operation supplies no such closure.
    pub fn continue_at_temperature(
        &mut self,
        refresh: &'s SourceJkrThermalRefresh,
        gate: &CancelGate,
    ) -> Result<JkrNormalTemperatureModalReport, JkrThermalModalError> {
        let thermal = &self.thermal;
        let heat = &self.heat;
        let (history, report) = self.mechanical.continue_at_temperature_with_admission(
            refresh,
            gate,
            |continuation, report| {
                let history = heat.propose_temperature_continuation(
                    thermal,
                    heat.stamp(),
                    continuation,
                    || checkpoint(gate),
                )?;
                // All mechanical and thermal changes are now staged. No cell
                // operation is needed: only original contact storage changes.
                Ok::<_, JkrThermalModalError>((history, *report))
            },
        )?;
        self.heat = history;
        Ok(report)
    }

    pub fn step(
        &mut self,
        external_modal_forces: &[f64],
    ) -> Result<&JkrThermalModalFrame, JkrThermalModalError> {
        self.step_under_gate(external_modal_forces, &CancelGate::new_clock_free())
    }

    /// Publish all three accepted owners only after original mechanics,
    /// event-timed thermal transport and their shared-energy account succeed.
    /// A source/thermal/energy refusal or cancellation leaves body states,
    /// clocks, spring storage, contact/heat history and observations unchanged.
    pub fn step_under_gate(
        &mut self,
        external_modal_forces: &[f64],
        gate: &CancelGate,
    ) -> Result<&JkrThermalModalFrame, JkrThermalModalError> {
        let thermal = &self.thermal;
        let history = &self.heat;
        let controls = self.controls;
        let origin = self.thermal_time_at_motion_zero_s;
        let support_route = self.support_heat.as_ref();
        let (candidate, mechanical, transfer, support_transfer, residual, tolerance, stored) = self
            .mechanical
            .step_with_admission(external_modal_forces, gate, |mechanical, _| {
                let candidate = stage_heat(
                    history,
                    thermal,
                    mechanical,
                    support_route,
                    origin,
                    controls,
                    gate,
                )?;
                let transfer = finite(
                    (Dd::from_f64(candidate.report.snap.work.dissipated_work_j())
                        - Dd::from_f64(mechanical.released_energy_j))
                    .to_f64(),
                )?;
                energy_gate(transfer, controls.heat.transfer_energy_tolerance_j)?;
                let support_transfer = finite(
                    (Dd::from_f64(
                        candidate
                            .support
                            .as_ref()
                            .map_or(0.0, |s| s.work.dissipated_work_j()),
                    ) - Dd::from_f64(mechanical.damping_dissipation_j))
                    .to_f64(),
                )?;
                energy_gate(support_transfer, controls.heat.transfer_energy_tolerance_j)?;
                let residual = finite(
                    (Dd::from_f64(mechanical.energy_residual_j)
                        + Dd::from_f64(candidate.residual_j)
                        + Dd::from_f64(transfer)
                        + Dd::from_f64(support_transfer))
                    .to_f64(),
                )?;
                let tolerance = finite(
                    mechanical.energy_tolerance_j
                        + controls.total_energy_tolerance_j
                        + controls.heat.transfer_energy_tolerance_j
                            * if candidate.support.is_some() {
                                2.0
                            } else {
                                1.0
                            },
                )?;
                energy_gate(residual, tolerance)?;
                let stored = finite(
                    mechanical.network_energy_j
                        + mechanical.contact.reversible_energy_j
                        + candidate.thermal.energy_j(),
                )?;
                Ok::<_, JkrThermalModalError>((
                    candidate,
                    mechanical.clone(),
                    transfer,
                    support_transfer,
                    residual,
                    tolerance,
                    stored,
                ))
            })?;
        // Native mechanical publication has succeeded. Everything remaining
        // is an infallible move of the already-admitted receiving candidates.
        self.heat = candidate.history;
        self.thermal = candidate.thermal;
        self.last = Some(JkrThermalModalFrame {
            mechanical,
            heat: candidate.report,
            support_heat: candidate.support,
            external_transport_heat_j: candidate.external_heat_j,
            snap_transfer_residual_j: transfer,
            support_transfer_residual_j: support_transfer,
            combined_energy_residual_j: residual,
            combined_energy_tolerance_j: tolerance,
            stored_energy_j: stored,
        });
        Ok(self
            .last
            .as_ref()
            .expect("published joint JKR thermal frame"))
    }
}

struct StagedHeat<'s, 'm> {
    history: JkrSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrTimedSnapHeatReport,
    support: Option<JkrViscousHeatReport>,
    external_heat_j: f64,
    residual_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn stage_heat<'s, 'm>(
    history: &JkrSnapHeat<'s>,
    thermal: &DepositionSimulation<'m>,
    mechanical: &JkrNormalEventModalFrame,
    support_route: Option<&JkrViscousHeatRoute>,
    origin: f64,
    controls: JkrTimedSnapHeatControls,
    gate: &CancelGate,
) -> Result<StagedHeat<'s, 'm>, JkrThermalModalError> {
    if let Some(trace) = mechanical.dissipation_trace() {
        let route = support_route.ok_or(JkrThermalModalError::Invalid(
            "sealed support loss requires an explicit heat/outside destination",
        ))?;
        let candidate = history.propose_timed_damped_normal_motion(
            thermal,
            history.stamp(),
            trace,
            origin,
            route,
            controls,
            || checkpoint(gate),
        )?;
        let (history, thermal, report) = candidate.into_parts();
        Ok(StagedHeat {
            history,
            thermal,
            report: report.timed,
            support: Some(report.support),
            external_heat_j: report.external_transport_heat_j,
            residual_j: report.combined_energy_residual_j,
        })
    } else {
        let candidate = history.propose_timed_normal_motion(
            thermal,
            history.stamp(),
            &mechanical.path,
            origin,
            controls,
            || checkpoint(gate),
        )?;
        let (history, thermal, report) = candidate.into_parts();
        Ok(StagedHeat {
            history,
            thermal,
            external_heat_j: report.transport.external_heat_j,
            residual_j: report.combined_energy_residual_j,
            report,
            support: None,
        })
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), DepositionError> {
    if gate.is_requested() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, JkrThermalModalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrThermalModalError::Invalid(
            "unrepresentable joint JKR energy or clock",
        ))
    }
}
fn energy_gate(residual_j: f64, tolerance_j: f64) -> Result<(), JkrThermalModalError> {
    if residual_j.abs() > tolerance_j {
        Err(JkrThermalModalError::EnergyBalance {
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}
