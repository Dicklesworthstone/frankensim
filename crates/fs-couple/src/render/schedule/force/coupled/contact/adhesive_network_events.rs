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

use super::super::{component_energy, extension};
use super::adhesive_events::linear_network::{
    LinearNetworkEventConfig, LinearNetworkEventError, LinearNetworkEventFlow,
    LinearNetworkEventStep, LinearNetworkLink, LinearNetworkPhase,
};
use super::adhesive_events::native::{
    EventNativeAccount, projected_position, projected_velocity, stage_account,
};
use super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel};
use crate::port_adoption::PortTransfer;
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
    positions: Vec<f64>,
    velocities: Vec<f64>,
    candidate_states: Vec<ModalAcousticState>,
    connection_work: Vec<PortTransfer>,
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
        if network.samples_rendered() != 0
            || n > config.events.maximum_modes
            || network.connections.len() > config.events.maximum_connections
            || column.len() != n
            || column.iter().any(|b| !b.is_finite())
            || !network.relaxation_branches.is_empty()
            || !network.plane_stress_material_points.is_empty()
            || network.radiation.is_some()
            || network
                .connections
                .iter()
                .any(|link| link.constant_loss_reference_omega_rad_s.is_some())
        {
            return Err(invalid("adhesive network events require a complete sample-zero linear modal network with Kelvin-Voigt connections and no constitutive/acoustic memory").into());
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
        let stiffness = network
            .models
            .iter()
            .flat_map(|model| model.modes())
            .map(|mode| finite(mode.angular_frequency_rad_s * mode.angular_frequency_rad_s))
            .collect::<Result<Vec<_>, _>>()?;
        let damping = network
            .models
            .iter()
            .flat_map(|model| model.modes())
            .map(|mode| finite(2.0 * mode.damping_ratio * mode.angular_frequency_rad_s))
            .collect::<Result<Vec<_>, _>>()?;
        let links = network
            .columns
            .iter()
            .zip(&network.connections)
            .map(|(b, link)| LinearNetworkLink {
                column: b.clone(),
                stiffness_n_per_m: link.stiffness_n_m,
                damping_n_s_per_m: link.damping_n_s_m,
                rest_extension_m: link.rest_extension_m,
            })
            .collect();
        let flow =
            LinearNetworkEventFlow::new(stiffness, damping, links, column.clone(), config.events)?;
        let connection_work = vec![PortTransfer::work(0.0, false); network.connections.len()];
        let candidate_states = network
            .models
            .iter()
            .flat_map(|model| model.states())
            .copied()
            .collect();
        let owner = Self {
            network,
            normal,
            shear,
            contact,
            indentation_m,
            column,
            flow,
            phase: LinearNetworkPhase::Fresh,
            positions: vec![0.0; n],
            velocities: vec![0.0; n],
            candidate_states,
            connection_work,
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
        let sample = self
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| invalid("sample clock overflow"))?;
        let before = self.network.total_energy_j()?;
        for (k, old) in self
            .network
            .models
            .iter()
            .flat_map(|model| model.states())
            .enumerate()
        {
            self.positions[k] = old.displacement_m_sqrt_kg;
            self.velocities[k] = old.velocity_m_sqrt_kg_per_s;
        }
        let static_capacity =
            finite(self.contact.contact_area_m2 * self.shear.law().static_shear_strength_pa())?;
        let kinetic_capacity =
            finite(self.contact.contact_area_m2 * self.shear.law().kinetic_shear_strength_pa())?;
        let path = self.flow.advance(
            &self.positions,
            &self.velocities,
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
        // Passivity supplies a bound over the entire trajectory, including
        // between resolved events. Original rest offsets remain in E itself.
        let root_energy = finite(
            det::sqrt(2.0 * before) + self.network.dt * det::sqrt(dot(external, external)?),
        )?;
        limit(
            "junction speed over complete trajectory",
            finite(det::sqrt(dot(&self.column, &self.column)?) * root_energy)?,
            self.config.maximum_tangent_speed_m_s,
        )?;

        for k in 0..self.mode_count() {
            self.candidate_states[k] = ModalAcousticState {
                displacement_m_sqrt_kg: self.flow.positions()[k],
                velocity_m_sqrt_kg_per_s: self.flow.velocities()[k],
            };
            self.network.old_q[k] = self.positions[k];
            self.network.free_delta[k] = finite(self.flow.positions()[k] - self.positions[k])?;
        }
        let mut pressure = 0.0;
        for i in 0..self.network.models.len() {
            poll(Some(gate))?;
            let candidate = &mut self.network.candidates[i];
            candidate
                .restore_states(
                    &self.candidate_states[self.network.offsets[i]..self.network.offsets[i + 1]],
                )
                .map_err(|source| ModalCouplingError::Component {
                    component: i,
                    source,
                })?;
            let mut component_pressure = 0.0;
            for (mode, state) in candidate.modes().iter().zip(candidate.states()) {
                component_pressure = finite(
                    component_pressure
                        + mode.pressure_per_modal_velocity.re * state.velocity_m_sqrt_kg_per_s
                        + mode.pressure_per_modal_velocity.im
                            * mode.angular_frequency_rad_s
                            * state.displacement_m_sqrt_kg,
                )?;
            }
            limit(
                "original component observer pressure",
                component_pressure.abs(),
                candidate.budget().maximum_abs_pressure_pa,
            )?;
            pressure = finite(pressure + component_pressure)?;
        }
        limit(
            "original combined observer pressure",
            pressure.abs(),
            self.network.config.maximum_abs_pressure_pa,
        )?;
        let modal_energy = component_energy(&self.network.candidates)?;
        let mut spring_energy = 0.0;
        let mut connection_loss = 0.0;
        for (i, (column, link)) in self
            .network
            .columns
            .iter()
            .zip(&self.network.connections)
            .enumerate()
        {
            poll(Some(gate))?;
            let z0 = extension(&self.network.models, column, link.rest_extension_m)?;
            let z1 = extension(&self.network.candidates, column, link.rest_extension_m)?;
            let speed = projected_velocity(&self.network.candidates, column)?;
            let loss = finite(self.flow.connection_dissipation_j()[i])?;
            if loss < 0.0 {
                return Err(invalid("negative original connection loss").into());
            }
            let change = finite(0.5 * link.stiffness_n_m * (z1 - z0) * (z1 + z0))?;
            self.connection_work[i] = PortTransfer::work(finite(-change - loss)?, false);
            connection_loss = finite(connection_loss + loss)?;
            spring_energy = finite(spring_energy + 0.5 * link.stiffness_n_m * z1 * z1)?;
            self.network.reactions[i] =
                finite(-link.stiffness_n_m * z1 - link.damping_n_s_m * speed)?;
            let bound = finite(
                (det::sqrt(link.stiffness_n_m)
                    + link.damping_n_s_m * det::sqrt(dot(column, column)?))
                    * root_energy,
            )?;
            limit(
                "original connection force over complete event trajectory",
                bound,
                self.network.config.maximum_abs_connection_force_n,
            )?;
        }
        let external_work = dot(external, &self.network.free_delta)?;
        let after = finite(modal_energy + spring_energy)?;
        let loss = finite(
            path.modal_damping_dissipation_j + connection_loss + path.friction_dissipation_j,
        )?;
        let uncertainty = finite(path.numerical_work_uncertainty_j)?;
        if path.modal_damping_dissipation_j < 0.0
            || path.friction_dissipation_j < 0.0
            || uncertainty < 0.0
        {
            return Err(invalid(
                "nonnegative modal/junction loss and numerical uncertainty required",
            )
            .into());
        }
        limit(
            "original network event energy",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        limit(
            "combined adhesive event energy magnitude",
            finite(after + self.contact.reversible_energy_j)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let residual = finite(after - before + loss - external_work)?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(external_work.abs())
            .max(loss)
            .max(self.contact.reversible_energy_j.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        for residual_j in [
            residual,
            before - path.initial_energy_j,
            modal_energy - path.modal_energy_j,
            spring_energy - path.connection_energy_j,
            external_work - path.external_work_j,
            connection_loss - path.connection_damping_dissipation_j,
        ] {
            if finite(residual_j.abs() + uncertainty)? > tolerance {
                return Err(AdhesiveNetworkEventModalError::EnergyAccuracy {
                    residual_j,
                    estimated_error_j: uncertainty,
                    tolerance_j: tolerance,
                });
            }
        }
        stage_account(
            &mut self.network,
            EventNativeAccount {
                sample,
                observer_pressure_pa: pressure,
                energy_before_j: before,
                modal_energy_j: modal_energy,
                component_dissipation_j: path.modal_damping_dissipation_j,
                connection_energy_j: spring_energy,
                connection_dissipation_j: connection_loss,
                work_j: finite(external_work - path.friction_dissipation_j)?,
                residual_j: residual,
                tolerance_j: tolerance,
            },
            self.connection_work.iter().copied(),
        );
        let accepted = AdhesiveNetworkEventModalFrame {
            sample,
            observer_pressure_pa: pressure,
            network_energy_j: after,
            surface_energy_j: self.contact.reversible_energy_j,
            external_work_j: external_work,
            friction_dissipation_j: path.friction_dissipation_j,
            modal_damping_dissipation_j: path.modal_damping_dissipation_j,
            connection_damping_dissipation_j: connection_loss,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            numerical_work_uncertainty_j: uncertainty,
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
