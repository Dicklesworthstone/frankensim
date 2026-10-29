//! Prepared, power-conjugate mean-port connection to the original DG owner.
//!
//! This is a held-pressure/held-flow coupling clock, not a joint continuous
//! exponential. Each participant retains its own exact propagator and memory.
//! The scalar mobility is prepared once; no response fitting or root iteration
//! is performed during a step. The owning physical adapter keeps the fixed
//! source, column, clock and duct together and exposes no mutable participant.

use super::{
    CoupledModalFrame, CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::dg_port::{DgPortLoad, sum_terminal_work};
use fs_duct::dg::DgStepOutput;
use fs_exec::CancelGate;
pub use super::exact::dg::ExactDgBudget;
pub(crate) use super::exact::dg::PreparedJointDg;

/// A refusal from either original owner or their shared pressure/work admission.
#[derive(Debug)]
pub enum ModalDgError {
    /// Original mechanics, cancellation, or joint physical budget refusal.
    Modal(ModalCouplingError),
    /// Original duct, wall, branch, or terminal propagation refusal.
    Acoustic(AcousticRealizeError),
}
impl From<ModalCouplingError> for ModalDgError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Modal(error)
    }
}
impl From<AcousticRealizeError> for ModalDgError {
    fn from(error: AcousticRealizeError) -> Self {
        Self::Acoustic(error)
    }
}
impl core::fmt::Display for ModalDgError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Modal(e) => e.fmt(f),
            Self::Acoustic(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for ModalDgError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Modal(e) => Some(e),
            Self::Acoustic(e) => Some(e),
        }
    }
}

/// One accepted coupled step. Internal inlet work cancels in the joint account.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalDgFrame {
    /// True when one continuous modal/DG exponential advances all coordinates.
    pub exact_joint_evolution: bool,
    /// Actual integrated inlet work [J]. For exact joint evolution this is
    /// integral p(t) Q(t) dt, generally different from the product of means.
    pub inlet_work_j: f64,
    /// Original modal frame. Its external work includes the opposing DG effort;
    /// use `authored_external_work_j` for the external work of the joint system.
    pub modal: CoupledModalFrame,
    /// Original aggregate duct/branch/gap ledger and physical port observations.
    pub duct: DgStepOutput,
    /// Swept volume divided by the shared duration, INTO the downstream duct.
    pub mean_flow_m3_s: f64,
    /// Mean inlet pressure opposing this volume flow [Pa], not microphone sound.
    pub mean_pressure_pa: f64,
    /// Accepted mechanical, material, duct, wall, branch and gap storage [J].
    /// Rational terminal storage is not certified by its owner and is excluded.
    pub stored_energy_j: f64,
    /// Change in that same joint storage over this step [J].
    pub storage_change_j: f64,
    /// Authored generalized force work only; internal duct reaction cancels [J].
    pub authored_external_work_j: f64,
    /// Original modal loss plus DG loss [J], independently integrated on the
    /// exact joint backend and inferred on the dense mean-port backend. When the DG ledger has
    /// `balance_checked == false`, this is not an independent wall heat integral.
    pub inferred_total_dissipation_j: f64,
    /// Signed work into retained rational terminals, excluding stored gap work.
    /// This is boundary exchange, not another loss or a rational energy claim.
    pub terminal_work_j: f64,
    /// Original integrated DG inlet work minus opposing modal work [J].
    pub port_work_residual_j: f64,
    /// Joint storage change + loss + terminal exchange - authored work [J].
    pub balance_residual_j: f64,
    /// Original absolute/relative energy admission applied to the joint account.
    pub balance_tolerance_j: f64,
}

/// Private preparation retained by an owning physical connection. Reusing it
/// with an unrelated network is deliberately not a public operation.
pub(crate) struct PreparedModalDg {
    flow_column: Vec<f64>,
    applied_forces: Vec<f64>,
    mobility_m3_pa_s: f64,
    impedance_pa_s_m3: f64,
    denominator: f64,
}
impl PreparedModalDg {
    pub(crate) fn new(
        network: &CoupledModalSystem,
        load: &DgPortLoad,
        flow_column: Vec<f64>,
        gate: &CancelGate,
    ) -> Result<Self, ModalDgError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0
            || load.state().step_count() != 0
            || network.sample_period_s().to_bits() != load.propagator().dt().to_bits()
            || load.has_driven_mouth()
            || flow_column.len() != network.mode_count()
            || flow_column.iter().any(|b| !b.is_finite())
        {
            return Err(invalid("DG mean coupling requires cold participants, one clock, a complete finite volume column and an ordinary inlet").into());
        }
        let response = network.linear_response(&flow_column, gate)?;
        let mobility_m3_pa_s = finite(dot(&flow_column, &response)? / network.dt)?;
        let impedance_pa_s_m3 = load.port_response()?.impedance_pa_s_m3;
        if mobility_m3_pa_s <= 0.0 || !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0 {
            return Err(invalid("DG mean coupling requires positive finite displacement mobility and inlet impedance").into());
        }
        let denominator = finite(1.0 + mobility_m3_pa_s * impedance_pa_s_m3)?;
        let initial = finite(network.total_energy_j()? + load.stored_energy_j())?;
        limit(
            "modal and DG total energy",
            initial,
            network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(Self {
            applied_forces: vec![0.0; flow_column.len()],
            flow_column,
            mobility_m3_pa_s,
            impedance_pa_s_m3,
            denominator,
        })
    }

    pub(crate) fn applied_forces(&self) -> &[f64] {
        &self.applied_forces
    }

    pub(crate) fn step(
        &mut self,
        network: &mut CoupledModalSystem,
        load: &mut DgPortLoad,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<ModalDgFrame, ModalDgError> {
        poll(Some(gate))?;
        let relation = load.port_response()?;
        if network.samples_rendered() != load.state().step_count()
            || relation.impedance_pa_s_m3.to_bits() != self.impedance_pa_s_m3.to_bits()
        {
            return Err(
                invalid("prepared DG source clock or fixed inlet impedance changed").into(),
            );
        }
        let modal_before = network.total_energy_j()?;
        network.prepare_prediction(external, Some(gate))?;
        let free_flow = finite(dot(&self.flow_column, &network.free_delta)? / network.dt)?;
        let flow = finite(
            (free_flow - self.mobility_m3_pa_s * relation.history_pressure_pa) / self.denominator,
        )?;
        let pressure = finite(relation.history_pressure_pa + self.impedance_pa_s_m3 * flow)?;
        for ((applied, &authored), &b) in self
            .applied_forces
            .iter_mut()
            .zip(external)
            .zip(&self.flow_column)
        {
            *applied = finite(authored - b * pressure)?;
        }
        network.stage_inner(&self.applied_forces, Some(gate))?;
        let volume = dot(&self.flow_column, &network.free_delta)?;
        let actual_flow = finite(volume / network.dt)?;
        let scale = free_flow.abs() + (self.mobility_m3_pa_s * pressure).abs() + actual_flow.abs();
        if (actual_flow - flow).abs() > network.config.solve_relative_tolerance * scale.max(1e-300)
        {
            return Err(
                invalid("DG and modal candidates disagree on the same swept volume").into(),
            );
        }
        let authored_external_work_j = dot(external, &network.free_delta)?;
        let inlet = load.inlet();
        let acoustic = load.preview_step(actual_flow)?;
        let duct = acoustic.output();
        let port = duct
            .port_outputs
            .iter()
            .find(|port| port.id == inlet)
            .ok_or_else(|| invalid("DG candidate omitted the connected inlet"))?;
        let pressure_scale = relation.history_pressure_pa.abs()
            + (self.impedance_pa_s_m3 * actual_flow).abs()
            + pressure.abs();
        if !port.conjugate_rate.is_finite()
            || (port.conjugate_rate - pressure).abs() > 1e-10 * pressure_scale.max(1e-300)
        {
            return Err(invalid(
                "DG and modal candidates disagree on the same mean inlet pressure",
            )
            .into());
        }
        limit(
            "DG mean inlet pressure",
            pressure.abs(),
            network.config.maximum_abs_pressure_pa,
        )?;
        let modal = &network.staged_frame;
        let stored_energy_j = finite(modal.total_energy_j() + duct.ledger.stored_energy_j)?;
        let before = finite(modal_before + acoustic.previous_stored_energy_j())?;
        let storage_change_j = finite(stored_energy_j - before)?;
        let inferred_total_dissipation_j = finite(
            modal.total_dissipation_j()
                + (duct.ledger.cumulative_dissipated_j
                    - acoustic.previous_cumulative_dissipated_j()),
        )?;
        let terminal_work_j = finite(sum_terminal_work(acoustic.terminal_frames()))?;
        let inlet_work_j = acoustic.inlet_work_j();
        let port_work_residual_j = finite(inlet_work_j - pressure * volume)?;
        let balance_residual_j = finite(
            storage_change_j + inferred_total_dissipation_j + terminal_work_j
                - authored_external_work_j,
        )?;
        let scale = finite(
            before.abs()
                + stored_energy_j.abs()
                + authored_external_work_j.abs()
                + inferred_total_dissipation_j.abs()
                + terminal_work_j.abs(),
        )?;
        let balance_tolerance_j = finite(
            network.config.energy_absolute_tolerance_j
                + network.config.energy_relative_tolerance * scale,
        )?;
        limit(
            "modal and DG total energy",
            stored_energy_j,
            network.config.maximum_total_energy_j,
        )?;
        for residual in [port_work_residual_j, balance_residual_j] {
            if residual.abs() > balance_tolerance_j {
                return Err(ModalCouplingError::EnergyBalance {
                    residual_j: residual,
                    tolerance_j: balance_tolerance_j,
                }
                .into());
            }
        }
        let modal = modal.clone();
        network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        // Both accepts are infallible. No cancellation or error may separate them.
        let duct = acoustic.accept();
        network.publish_staged();
        Ok(ModalDgFrame {
            exact_joint_evolution: false,
            inlet_work_j,
            modal,
            duct,
            mean_flow_m3_s: actual_flow,
            mean_pressure_pa: pressure,
            stored_energy_j,
            storage_change_j,
            authored_external_work_j,
            inferred_total_dissipation_j,
            terminal_work_j,
            port_work_residual_j,
            balance_residual_j,
            balance_tolerance_j,
        })
    }
}
