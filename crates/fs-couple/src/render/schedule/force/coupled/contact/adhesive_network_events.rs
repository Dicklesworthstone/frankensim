//! Sourced adhesive events on the original complete linear elastic network.
//!
//! Every modal spring, modal resistance and Kelvin--Voigt connection enters
//! the same matrix-free trajectory. During sticking, the mass-orthogonal
//! constraint acts only on the declared tangent; all other motion continues
//! and can change the required reaction until release. Original component
//! budgets, connection ports and pressure observations remain authoritative.
//!
//! Normal indentation, source conditions, coefficients and tangent geometry
//! are held. This owner admits one scalar junction, not simultaneous/vector
//! friction, nonlinear material memory or acoustic radiation back-reaction.

use super::adhesive_events::linear_network::{
    LinearNetworkEventConfig, LinearNetworkEventError, LinearNetworkEventFlow,
    LinearNetworkEventStep, LinearNetworkPhase,
};
use super::adhesive_events::native::linear_network::{
    NativeLinearEventError, NativeLinearEventTrajectory, NativeLinearNetworkWorkspace,
    network_parameters, require_linear_network,
};
use super::adhesive_events::native::{projected_position, projected_velocity};
use super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::{DmtError, DmtResponse, SourceDmtSphereContact};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveNetworkEventModalConfig {
    pub events: LinearNetworkEventConfig,
    pub maximum_contact_force_n: f64,
    pub maximum_tangent_speed_m_s: f64,
    pub maximum_accumulated_slip_m: f64,
}

#[derive(Debug)]
pub enum AdhesiveNetworkEventModalError {
    Network(ModalCouplingError),
    Normal(DmtError),
    Trajectory(LinearNetworkEventError),
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for AdhesiveNetworkEventModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native elastic-network adhesion refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveNetworkEventModalError {}
impl From<ModalCouplingError> for AdhesiveNetworkEventModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<DmtError> for AdhesiveNetworkEventModalError {
    fn from(error: DmtError) -> Self {
        Self::Normal(error)
    }
}
impl From<LinearNetworkEventError> for AdhesiveNetworkEventModalError {
    fn from(error: LinearNetworkEventError) -> Self {
        Self::Trajectory(error)
    }
}
impl From<NativeLinearEventError> for AdhesiveNetworkEventModalError {
    fn from(error: NativeLinearEventError) -> Self {
        match error {
            NativeLinearEventError::Network(error) => Self::Network(error),
            NativeLinearEventError::EnergyAccuracy {
                residual_j,
                estimated_error_j,
                tolerance_j,
            } => Self::EnergyAccuracy {
                residual_j,
                estimated_error_j,
                tolerance_j,
            },
        }
    }
}

/// Simultaneous sourced junctions sharing the same original native mechanics.
pub mod multiple;

/// Joint two-dimensional Coulomb disks with explicit time-discretization loss.
pub mod vector;

/// Complete accepted motion, with source capacities and distinct loss channels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveNetworkEventModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub friction_dissipation_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub tangent_position_m: f64,
    pub tangent_velocity_m_s: f64,
    pub static_capacity_n: f64,
    pub kinetic_capacity_n: f64,
    pub path: LinearNetworkEventStep,
}

pub struct AdhesiveNetworkEventModalSystem<'a> {
    network: CoupledModalSystem,
    normal: &'a SourceDmtSphereContact,
    shear: &'a SourceAdhesiveShear,
    contact: DmtResponse,
    indentation_m: f64,
    column: Vec<f64>,
    flow: LinearNetworkEventFlow,
    phase: LinearNetworkPhase,
    workspace: NativeLinearNetworkWorkspace,
    config: AdhesiveNetworkEventModalConfig,
    last: Option<AdhesiveNetworkEventModalFrame>,
}

impl<'a> AdhesiveNetworkEventModalSystem<'a> {
    /// Attach the original ordered DMT/shear source to a complete fixed tangent
    /// column at sample zero. Existing modal frequencies, resistances, link
    /// columns and rest offsets are retained without modal condensation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        surface_materials: [&MaterialStateId; 2],
        indentation_m: f64,
        column: Vec<f64>,
        config: AdhesiveNetworkEventModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, AdhesiveNetworkEventModalError> {
        poll(Some(gate))?;
        let n = network.mode_count();
        require_linear_network(
            &network,
            config.events.maximum_modes,
            config.events.maximum_connections,
        )?;
        if column.len() != n || column.iter().any(|b| !b.is_finite()) {
            return Err(invalid("adhesive work column must cover every original mode").into());
        }
        let materials = normal
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        if surface_materials[0] != &materials[0]
            || surface_materials[1] != &materials[1]
            || shear.normal_state_identity() != normal.normal_state().identity()
        {
            return Err(
                invalid("adhesive network source and receiving ordered surfaces differ").into(),
            );
        }
        for value in [
            config.maximum_contact_force_n,
            config.maximum_tangent_speed_m_s,
            config.maximum_accumulated_slip_m,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid("positive finite network junction limits required").into());
            }
        }
        let contact = normal.evaluate(indentation_m)?;
        limit(
            "adhesive static capacity",
            finite(contact.contact_area_m2 * shear.law().static_shear_strength_pa())?,
            config.maximum_contact_force_n,
        )?;
        let (stiffness, damping, links) = network_parameters(&network)?;
        let flow =
            LinearNetworkEventFlow::new(stiffness, damping, links, column.clone(), config.events)?;
        let workspace = NativeLinearNetworkWorkspace::new(&network);
        let owner = Self {
            network,
            normal,
            shear,
            contact,
            indentation_m,
            column,
            flow,
            phase: LinearNetworkPhase::Fresh,
            workspace,
            config,
            last: None,
        };
        limit(
            "initial adhesive network energy magnitude",
            owner.total_energy_j()?.abs(),
            owner.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(owner)
    }

    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.normal
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.shear
    }
    pub const fn indentation_m(&self) -> f64 {
        self.indentation_m
    }
    pub const fn phase(&self) -> LinearNetworkPhase {
        self.phase
    }
    pub fn work_column(&self) -> &[f64] {
        &self.column
    }
    pub const fn last_frame(&self) -> Option<&AdhesiveNetworkEventModalFrame> {
        self.last.as_ref()
    }
    pub const fn config(&self) -> AdhesiveNetworkEventModalConfig {
        self.config
    }
    pub fn total_energy_j(&self) -> Result<f64, AdhesiveNetworkEventModalError> {
        Ok(finite(
            self.network.total_energy_j()? + self.contact.reversible_energy_j,
        )?)
    }

    pub fn step(
        &mut self,
        external: &[f64],
    ) -> Result<&AdhesiveNetworkEventModalFrame, AdhesiveNetworkEventModalError> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&AdhesiveNetworkEventModalFrame, AdhesiveNetworkEventModalError> {
        poll(Some(gate))?;
        if external.len() != self.mode_count() || external.iter().any(|f| !f.is_finite()) {
            return Err(invalid("external event forces must cover every original mode").into());
        }
        self.workspace.capture(&self.network);
        let static_capacity =
            finite(self.contact.contact_area_m2 * self.shear.law().static_shear_strength_pa())?;
        let kinetic_capacity =
            finite(self.contact.contact_area_m2 * self.shear.law().kinetic_shear_strength_pa())?;
        let path = self.flow.advance(
            self.workspace.positions(),
            self.workspace.velocities(),
            external,
            self.phase,
            self.network.dt,
            static_capacity,
            kinetic_capacity,
            gate,
        )?;
        limit(
            "junction accumulated physical slip",
            path.accumulated_slip_m,
            self.config.maximum_accumulated_slip_m,
        )?;
        let native = self.workspace.stage(
            &mut self.network,
            external,
            self.contact.reversible_energy_j,
            NativeLinearEventTrajectory {
                positions: self.flow.positions(),
                velocities: self.flow.velocities(),
                connection_dissipation_j: self.flow.connection_dissipation_j(),
                initial_energy_j: path.initial_energy_j,
                modal_energy_j: path.modal_energy_j,
                connection_energy_j: path.connection_energy_j,
                external_work_j: path.external_work_j,
                modal_damping_dissipation_j: path.modal_damping_dissipation_j,
                connection_damping_dissipation_j: path.connection_damping_dissipation_j,
                friction_dissipation_j: path.friction_dissipation_j,
                numerical_work_uncertainty_j: path.numerical_work_uncertainty_j,
            },
            gate,
        )?;
        limit(
            "junction speed over complete trajectory",
            finite(det::sqrt(dot(&self.column, &self.column)?) * native.root_energy_bound)?,
            self.config.maximum_tangent_speed_m_s,
        )?;
        let accepted = AdhesiveNetworkEventModalFrame {
            sample: native.sample,
            observer_pressure_pa: native.observer_pressure_pa,
            network_energy_j: native.network_energy_j,
            surface_energy_j: self.contact.reversible_energy_j,
            external_work_j: native.external_work_j,
            friction_dissipation_j: path.friction_dissipation_j,
            modal_damping_dissipation_j: path.modal_damping_dissipation_j,
            connection_damping_dissipation_j: native.connection_damping_dissipation_j,
            energy_residual_j: native.energy_residual_j,
            energy_tolerance_j: native.energy_tolerance_j,
            numerical_work_uncertainty_j: path.numerical_work_uncertainty_j,
            tangent_position_m: projected_position(&self.network.candidates, &self.column)?,
            tangent_velocity_m_s: projected_velocity(&self.network.candidates, &self.column)?,
            static_capacity_n: static_capacity,
            kinetic_capacity_n: kinetic_capacity,
            path,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        self.phase = path.phase;
        Ok(self.last.insert(accepted))
    }
}
