//! Event-resolved adhesive sliding through an original native spring/damper.
//!
//! One original connection and its complete work column reduce the free body
//! network exactly to one relative oscillator plus unaffected common motion.
//! Its actual stiffness, dashpot and rest extension enter the sourced event
//! owner. No equivalent body or replacement connection is constructed.
//!
//! Normal indentation and the original spring anchor are held. External modal
//! loads may change each sample. The constitutive owner resolves all intrastep
//! stops and releases; complete reaction impulse and first moment propagate
//! actual bodies. Spring storage, physical damper loss and junction loss remain
//! separate, with original connection-port work and joint state publication.

use super::adhesive_events::native::{
    EventNativeAccount, projected_position, projected_velocity, stage_account, stage_motion,
};
use super::adhesive_events::{AdhesiveEventModalError, compatible};
use super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::elastic_events::{
    DmtElasticEventBudget, DmtElasticEventError, DmtElasticEventState, DmtElasticEventStep,
    DmtElasticSlipEvents,
};
use fs_contact::normal_patch::dmt::{DmtError, SourceDmtSphereContact};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_math::{c64::C64, det};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveElasticEventModalConfig {
    pub events: DmtElasticEventBudget,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub velocity_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum AdhesiveElasticEventModalError {
    Network(ModalCouplingError),
    Contact(DmtElasticEventError),
    Motion(AdhesiveEventModalError),
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for AdhesiveElasticEventModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native elastic adhesive event refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveElasticEventModalError {}
impl From<ModalCouplingError> for AdhesiveElasticEventModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<DmtElasticEventError> for AdhesiveElasticEventModalError {
    fn from(error: DmtElasticEventError) -> Self {
        Self::Contact(error)
    }
}
impl From<DmtError> for AdhesiveElasticEventModalError {
    fn from(error: DmtError) -> Self {
        Self::Contact(error.into())
    }
}
impl From<AdhesiveEventModalError> for AdhesiveElasticEventModalError {
    fn from(error: AdhesiveEventModalError) -> Self {
        Self::Motion(error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveElasticEventModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub spring_energy_j: f64,
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub friction_dissipation_j: f64,
    pub damping_dissipation_j: f64,
    pub reaction_work_residual_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub coordinate_residual_m: f64,
    pub velocity_residual_m_s: f64,
    pub connection_force_n: f64,
    /// Passive energy bound for the entire original connection trajectory.
    pub connection_force_bound_n: f64,
    pub path: DmtElasticEventStep,
}

pub struct AdhesiveElasticEventModalSystem<'a> {
    network: CoupledModalSystem,
    events: DmtElasticSlipEvents<'a>,
    state: DmtElasticEventState,
    column: Vec<f64>,
    candidate_states: Vec<ModalAcousticState>,
    config: AdhesiveElasticEventModalConfig,
    last: Option<AdhesiveElasticEventModalFrame>,
}

impl<'a> AdhesiveElasticEventModalSystem<'a> {
    /// Bind the contact to the ONE original connection's complete signed column.
    /// Free modes and zero acoustic transfers are required; a modal spring,
    /// extra connection, material memory or non-collinear loss law is refused.
    /// The supplied tangent coordinate attaches to actual initial geometry;
    /// original rest extension supplies the immutable spring-stage offset.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        surface_materials: [&MaterialStateId; 2],
        indentation_m: f64,
        initial_tangent_coordinate_m: f64,
        config: AdhesiveElasticEventModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, AdhesiveElasticEventModalError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0
            || network.connections.len() != 1
            || !network.relaxation_branches.is_empty()
            || !network.plane_stress_material_points.is_empty()
            || network.radiation.is_some()
            || network.models.iter().flat_map(|m| m.modes()).any(|mode| {
                mode.angular_frequency_rad_s != 0.0
                    || mode.damping_ratio != 0.0
                    || mode.pressure_per_modal_velocity != C64::ZERO
            })
        {
            return Err(invalid(
                "elastic adhesive events require free modes and exactly one original connection",
            )
            .into());
        }
        let link = &network.connections[0];
        if link.stiffness_n_m <= 0.0 || link.constant_loss_reference_omega_rad_s.is_some() {
            return Err(invalid(
                "positive original spring and collinear Kelvin-Voigt dashpot required",
            )
            .into());
        }
        let materials = normal
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        if surface_materials[0] != &materials[0] || surface_materials[1] != &materials[1] {
            return Err(invalid(
                "elastic adhesive receiving surfaces differ from the ordered source",
            )
            .into());
        }
        if [
            config.coordinate_absolute_tolerance_m,
            config.velocity_absolute_tolerance_m_s,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
            || [
                config.coordinate_relative_tolerance,
                config.velocity_relative_tolerance,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(
                invalid("positive finite elastic event compatibility budgets required").into(),
            );
        }
        let column = network.columns[0].clone();
        let events = DmtElasticSlipEvents::new(
            normal,
            shear,
            finite(1.0 / dot(&column, &column)?)?,
            link.stiffness_n_m,
            link.damping_n_s_m,
            config.events,
        )?;
        let extension =
            finite(projected_position(network.components(), &column)? - link.rest_extension_m)?;
        let state = events.initialize(
            indentation_m,
            initial_tangent_coordinate_m,
            projected_velocity(network.components(), &column)?,
            finite(initial_tangent_coordinate_m - extension)?,
        )?;
        let candidate_states = network
            .models
            .iter()
            .flat_map(|m| m.states().iter().copied())
            .collect();
        let system = Self {
            network,
            events,
            state,
            column,
            candidate_states,
            config,
            last: None,
        };
        limit(
            "native elastic adhesive initial energy magnitude",
            system.total_energy_j()?.abs(),
            system.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(system)
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
    pub const fn state(&self) -> DmtElasticEventState {
        self.state
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.events.normal()
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.events.shear()
    }
    pub const fn last_frame(&self) -> Option<&AdhesiveElasticEventModalFrame> {
        self.last.as_ref()
    }
    pub fn work_column(&self) -> &[f64] {
        &self.column
    }
    pub const fn effective_mass_kg(&self) -> f64 {
        self.events.mass_kg()
    }
    pub const fn config(&self) -> AdhesiveElasticEventModalConfig {
        self.config
    }
    pub fn total_energy_j(&self) -> Result<f64, AdhesiveElasticEventModalError> {
        Ok(finite(
            self.network.total_energy_j()?
                + self
                    .normal()
                    .evaluate(self.state.indentation_m())?
                    .reversible_energy_j,
        )?)
    }

    pub fn step(
        &mut self,
        external: &[f64],
    ) -> Result<&AdhesiveElasticEventModalFrame, AdhesiveElasticEventModalError> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&AdhesiveElasticEventModalFrame, AdhesiveElasticEventModalError> {
        poll(Some(gate))?;
        if external.len() != self.mode_count() || external.iter().any(|f| !f.is_finite()) {
            return Err(
                invalid("external forces must match every original elastic-event mode").into(),
            );
        }
        let sample = self
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| invalid("sample clock overflow"))?;
        let dt = self.network.dt;
        let before = self.network.total_energy_j()?;
        let link = &self.network.connections[0];
        let (k, c, rest) = (
            link.stiffness_n_m,
            link.damping_n_s_m,
            link.rest_extension_m,
        );
        let extension_before = finite(projected_position(self.components(), &self.column)? - rest)?;
        self.check_coordinate(
            "initial original spring extension",
            extension_before,
            self.state.tangent_displacement_m() - self.state.stage_position_m(),
        )?;
        compatible(
            "initial elastic event velocity",
            projected_velocity(self.components(), &self.column)?,
            self.state.tangent_velocity_m_s(),
            self.config.velocity_absolute_tolerance_m_s,
            self.config.velocity_relative_tolerance,
        )?;
        let path = self.events.advance(
            self.state,
            self.state.stage_position_m(),
            finite(self.effective_mass_kg() * dot(&self.column, external)?)?,
            dt,
            gate,
        )?;
        let motion = stage_motion(
            &mut self.network,
            &mut self.candidate_states,
            &self.column,
            external,
            path.reaction_impulse_n_s,
            path.reaction_time_moment_n_s2,
            gate,
        )?;
        let coordinate_residual = self.check_coordinate(
            "actual elastic event travel",
            motion.tangent_increment_m,
            path.tangent_increment_m,
        )?;
        let velocity_residual = compatible(
            "actual elastic event velocity",
            motion.tangent_velocity_m_s,
            path.state().tangent_velocity_m_s(),
            self.config.velocity_absolute_tolerance_m_s,
            self.config.velocity_relative_tolerance,
        )?;
        let extension_after =
            finite(projected_position(&self.network.candidates, &self.column)? - rest)?;
        self.check_coordinate(
            "actual original spring extension",
            extension_after,
            path.state().tangent_displacement_m() - path.state().stage_position_m(),
        )?;
        let spring_after = finite(0.5 * k * extension_after * extension_after)?;
        let spring_change = finite(
            0.5 * k * (extension_after - extension_before) * (extension_after + extension_before),
        )?;
        let after = finite(motion.modal_energy_j + spring_after)?;
        let surface_energy = path.contact.reversible_energy_j;
        limit(
            "original elastic event network energy",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        limit(
            "elastic event signed energy magnitude",
            finite(after + surface_energy)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let friction_loss = finite(path.friction_dissipation_j)?;
        let damping_loss = finite(path.damping_dissipation_j)?;
        let uncertainty = finite(path.numerical_work_uncertainty_j)?;
        if friction_loss < 0.0 || damping_loss < 0.0 || uncertainty < 0.0 {
            return Err(
                invalid("nonnegative physical losses and numerical uncertainty required").into(),
            );
        }
        let kinetic_change = finite(motion.modal_energy_j - motion.modal_energy_before_j)?;
        let body_residual =
            finite(kinetic_change - motion.external_work_j - motion.reaction_work_j)?;
        let reaction_residual =
            finite(motion.reaction_work_j + spring_change + damping_loss + friction_loss)?;
        let residual = finite(
            kinetic_change + spring_change + damping_loss + friction_loss - motion.external_work_j,
        )?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(surface_energy.abs())
            .max(motion.external_work_j.abs())
            .max(damping_loss + friction_loss);
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        for residual_j in [body_residual, reaction_residual, residual] {
            if finite(residual_j.abs() + uncertainty)? > tolerance {
                return Err(AdhesiveElasticEventModalError::EnergyAccuracy {
                    residual_j,
                    estimated_error_j: uncertainty,
                    tolerance_j: tolerance,
                });
            }
        }
        // For passive spring/damper/friction dynamics, sqrt(2E(t)) is bounded
        // by sqrt(2E0)+dt*||external||. This bounds the connection effort over
        // the full sample, rather than checking only its endpoint.
        let force_bound = finite(
            (det::sqrt(k) + c * det::sqrt(dot(&self.column, &self.column)?))
                * (det::sqrt(2.0 * before) + dt * det::sqrt(dot(external, external)?)),
        )?;
        limit(
            "original connection effort over event trajectory",
            force_bound,
            self.network.config.maximum_abs_connection_force_n,
        )?;
        let connection_force = finite(-k * extension_after - c * motion.tangent_velocity_m_s)?;
        self.network.reactions[0] = connection_force;
        // Junction work crosses the native subsystem boundary; its physical
        // dashpot loss and spring storage remain inside that original owner.
        let native_work = finite(motion.external_work_j - friction_loss)?;
        stage_account(
            &mut self.network,
            EventNativeAccount {
                sample,
                energy_before_j: before,
                modal_energy_j: motion.modal_energy_j,
                connection_energy_j: spring_after,
                connection_dissipation_j: damping_loss,
                work_j: native_work,
                residual_j: residual,
                tolerance_j: tolerance,
            },
            [crate::port_adoption::PortTransfer::work(
                finite(-spring_change - damping_loss)?,
                false,
            )],
        );
        let accepted = AdhesiveElasticEventModalFrame {
            sample,
            observer_pressure_pa: 0.0,
            network_energy_j: after,
            spring_energy_j: spring_after,
            surface_energy_j: surface_energy,
            external_work_j: motion.external_work_j,
            friction_dissipation_j: friction_loss,
            damping_dissipation_j: damping_loss,
            reaction_work_residual_j: reaction_residual,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            numerical_work_uncertainty_j: uncertainty,
            coordinate_residual_m: coordinate_residual,
            velocity_residual_m_s: velocity_residual,
            connection_force_n: connection_force,
            connection_force_bound_n: force_bound,
            path,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        self.state = path.state();
        Ok(self.last.insert(accepted))
    }

    fn check_coordinate(
        &self,
        what: &'static str,
        actual: f64,
        expected: f64,
    ) -> Result<f64, AdhesiveElasticEventModalError> {
        Ok(compatible(
            what,
            actual,
            expected,
            self.config.coordinate_absolute_tolerance_m,
            self.config.coordinate_relative_tolerance,
        )?)
    }
}
