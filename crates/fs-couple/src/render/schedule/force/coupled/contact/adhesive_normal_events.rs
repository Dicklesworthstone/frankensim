//! Source JKR attachment and pull-off on the original free modal coordinates.
//!
//! The scalar owner condenses the kinetic/contact discrete gradient, then its
//! impulse and first time moment restore every original native coordinate.
//! Approach evolves with the bodies. A branch jump changes contact storage but
//! supplies no impulse: its released energy remains an explicit outward account.
//! One fixed normal and a frozen quasistatic JKR source are required. Elastic
//! networks, rotating normals, resolved contact waves and thermal feedback need
//! their own propagation; they are not represented by free coordinates here.

use super::adhesive_events::native::{
    EventNativeAccount, projected_velocity, stage_account, stage_motion,
};
use super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrError, JkrResponse, JkrState, SourceJkrSphereContact,
    motion::{
        JkrNormalMotion, JkrNormalMotionBudget, JkrNormalMotionError, JkrNormalMotionState,
        JkrNormalMotionStep,
    },
};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;

#[derive(Debug)]
pub enum JkrNormalEventModalError {
    Network(ModalCouplingError),
    Motion(JkrNormalMotionError),
    Contact(JkrError),
    Compatibility {
        what: &'static str,
        residual: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for JkrNormalEventModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native JKR normal motion refused: {self:?}")
    }
}
impl std::error::Error for JkrNormalEventModalError {}
impl From<ModalCouplingError> for JkrNormalEventModalError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Network(value)
    }
}
impl From<JkrNormalMotionError> for JkrNormalEventModalError {
    fn from(value: JkrNormalMotionError) -> Self {
        Self::Motion(value)
    }
}
impl From<JkrError> for JkrNormalEventModalError {
    fn from(value: JkrError) -> Self {
        Self::Contact(value)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalEventModalConfig {
    pub motion: JkrNormalMotionBudget,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub velocity_relative_tolerance: f64,
}

/// One accepted native interval. Release is not silently classified as heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalEventModalFrame {
    pub sample: u64,
    pub network_energy_j: f64,
    /// Original source evaluated at the actual restored native approach.
    pub contact: JkrResponse,
    pub external_work_j: f64,
    pub normal_work_on_network_j: f64,
    pub released_energy_j: f64,
    pub normal_work_residual_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub coordinate_residual_m: f64,
    pub velocity_residual_m_s: f64,
    pub path: JkrNormalMotionStep,
}

pub struct JkrNormalEventModalSystem<'a> {
    network: CoupledModalSystem,
    motion: JkrNormalMotion<'a>,
    state: JkrNormalMotionState,
    column: Vec<f64>,
    reference_positions: Vec<f64>,
    initial_indentation_m: f64,
    candidate_states: Vec<ModalAcousticState>,
    config: JkrNormalEventModalConfig,
    last: Option<JkrNormalEventModalFrame>,
}

impl<'a> JkrNormalEventModalSystem<'a> {
    /// Attach one source to a complete mass-normalized native work column B.
    /// Initial history defines approach at the current coordinates; later
    /// approach is delta_initial + B·(q-q_initial). The relative mass is
    /// 1/(B·B), while common and uncontacted motion stay in the native owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        normal: &'a SourceJkrSphereContact,
        surface_materials: [&MaterialStateId; 2],
        initial_contact: JkrState,
        column: Vec<f64>,
        config: JkrNormalEventModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, JkrNormalEventModalError> {
        poll(Some(gate))?;
        network.require_time_domain()?;
        if network.samples_rendered() != 0
            || !network.connections.is_empty()
            || !network.relaxation_branches.is_empty()
            || !network.plane_stress_material_points.is_empty()
            || network.radiation.is_some()
            || network
                .models
                .iter()
                .flat_map(|m| m.damping_rates_per_s())
                .any(|&rate| rate != 0.0)
            || network.models.iter().flat_map(|m| m.modes()).any(|m| {
                m.angular_frequency_rad_s != 0.0
                    || m.pressure_per_modal_velocity != fs_math::c64::C64::ZERO
            })
        {
            return Err(invalid("JKR normal events require a sample-zero free network without connections, memory or radiation").into());
        }
        let materials = normal
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        if surface_materials[0] != &materials[0] || surface_materials[1] != &materials[1] {
            return Err(
                invalid("JKR receiving surfaces differ from the ordered normal source").into(),
            );
        }
        if column.len() != network.mode_count()
            || column.iter().any(|x| !x.is_finite())
            || [
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
            return Err(invalid("complete finite JKR work column and positive native compatibility budgets required").into());
        }
        let motion =
            JkrNormalMotion::new(normal, finite(1.0 / dot(&column, &column)?)?, config.motion)?;
        let state = motion.initialize(
            initial_contact,
            projected_velocity(network.components(), &column)?,
        )?;
        let reference_positions = network
            .models
            .iter()
            .flat_map(|m| m.states())
            .map(|s| s.displacement_m_sqrt_kg)
            .collect();
        let candidate_states = network
            .models
            .iter()
            .flat_map(|m| m.states().iter().copied())
            .collect();
        let owner = Self {
            network,
            motion,
            state,
            column,
            reference_positions,
            initial_indentation_m: initial_contact.indentation_m(),
            candidate_states,
            config,
            last: None,
        };
        limit(
            "native JKR initial energy magnitude",
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
    pub const fn state(&self) -> JkrNormalMotionState {
        self.state
    }
    pub const fn normal(&self) -> &'a SourceJkrSphereContact {
        self.motion.normal()
    }
    pub const fn effective_mass_kg(&self) -> f64 {
        self.motion.effective_mass_kg()
    }
    pub fn work_column(&self) -> &[f64] {
        &self.column
    }
    pub const fn config(&self) -> JkrNormalEventModalConfig {
        self.config
    }
    pub const fn last_frame(&self) -> Option<&JkrNormalEventModalFrame> {
        self.last.as_ref()
    }

    pub fn total_energy_j(&self) -> Result<f64, JkrNormalEventModalError> {
        let contact =
            self.native_contact(self.components(), self.state.contact_state().branch())?;
        Ok(finite(
            self.network.total_energy_j()? + contact.reversible_energy_j,
        )?)
    }

    /// Original bodies, clock and contact history publish together. Released
    /// jump energy is an explicit outward transfer; thermal deposition requires
    /// a separately declared recipient and fraction.
    pub fn step(
        &mut self,
        external: &[f64],
    ) -> Result<&JkrNormalEventModalFrame, JkrNormalEventModalError> {
        self.step_under_gate(external, &CancelGate::new_clock_free())
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&JkrNormalEventModalFrame, JkrNormalEventModalError> {
        let accepted = self.stage_step(external, gate)?;
        self.publish(accepted);
        Ok(self.last.as_ref().expect("published JKR native frame"))
    }

    /// Stage receiving thermal/other state before any native body or contact
    /// history publishes. The callback must return a candidate, not mutate an
    /// accepted external owner. A final cancellation check precedes publication.
    /// Its returned candidate can then be installed by an infallible assignment.
    pub fn step_with_admission<R, E>(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
        admit: impl FnOnce(&JkrNormalEventModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<JkrNormalEventModalError>,
    {
        let accepted = self.stage_step(external, gate).map_err(E::from)?;
        let candidate = admit(&accepted, &self.network.candidates)?;
        poll(Some(gate))
            .map_err(JkrNormalEventModalError::from)
            .map_err(E::from)?;
        self.publish(accepted);
        Ok(candidate)
    }

    fn publish(&mut self, accepted: JkrNormalEventModalFrame) {
        self.network.publish_staged();
        self.state = accepted.path.state();
        self.last = Some(accepted);
    }

    fn stage_step(
        &mut self,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<JkrNormalEventModalFrame, JkrNormalEventModalError> {
        poll(Some(gate))?;
        if external.len() != self.mode_count() || external.iter().any(|x| !x.is_finite()) {
            return Err(
                invalid("held external forces must cover every original finite mode").into(),
            );
        }
        let sample = self
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| invalid("JKR native clock overflow"))?;
        let before = self.network.total_energy_j()?;
        let old = self.native_contact(self.components(), self.state.contact_state().branch())?;
        compatible(
            "accepted native approach",
            old.indentation_m,
            self.state.contact_state().indentation_m(),
            self.config.coordinate_absolute_tolerance_m,
            self.config.coordinate_relative_tolerance,
        )?;
        compatible(
            "accepted native normal velocity",
            projected_velocity(self.components(), &self.column)?,
            self.state.velocity_m_s(),
            self.config.velocity_absolute_tolerance_m_s,
            self.config.velocity_relative_tolerance,
        )?;
        let force = finite(self.effective_mass_kg() * dot(&self.column, external)?)?;
        let path = self
            .motion
            .advance(self.state, force, self.network.dt, gate)?;
        let motion = stage_motion(
            &mut self.network,
            &mut self.candidate_states,
            &self.column,
            external,
            path.normal_impulse_n_s,
            path.normal_time_moment_n_s2,
            gate,
        )?;
        let contact = self.native_contact(
            &self.network.candidates,
            path.state().contact_state().branch(),
        )?;
        let coordinate_residual = compatible(
            "restored native approach",
            contact.indentation_m,
            path.state().contact_state().indentation_m(),
            self.config.coordinate_absolute_tolerance_m,
            self.config.coordinate_relative_tolerance,
        )?;
        compatible(
            "restored native normal increment",
            motion.tangent_increment_m,
            path.normal_increment_m,
            self.config.coordinate_absolute_tolerance_m,
            self.config.coordinate_relative_tolerance,
        )?;
        let velocity_residual = compatible(
            "restored native normal velocity",
            motion.tangent_velocity_m_s,
            path.state().velocity_m_s(),
            self.config.velocity_absolute_tolerance_m_s,
            self.config.velocity_relative_tolerance,
        )?;
        let after = motion.modal_energy_j;
        limit(
            "native JKR kinetic energy",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        limit(
            "native JKR energy magnitude",
            finite(after + contact.reversible_energy_j)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let work = finite(motion.external_work_j + motion.reaction_work_j)?;
        let native_residual = finite(after - before - work)?;
        let normal_residual = finite(
            motion.reaction_work_j + contact.reversible_energy_j - old.reversible_energy_j
                + path.released_energy_j,
        )?;
        let residual = finite(
            after - before + contact.reversible_energy_j - old.reversible_energy_j
                + path.released_energy_j
                - motion.external_work_j,
        )?;
        let contact_scale = old
            .reversible_energy_j
            .abs()
            .max(contact.reversible_energy_j.abs())
            .max(motion.reaction_work_j.abs())
            .max(path.released_energy_j);
        let contact_tolerance = finite(
            self.config.motion.energy_absolute_tolerance_j
                + self.config.motion.energy_relative_tolerance * contact_scale,
        )?;
        check_energy(normal_residual, contact_tolerance)?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(contact_scale)
            .max(work.abs())
            .max(motion.external_work_j.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        check_energy(native_residual, tolerance)?;
        check_energy(residual, tolerance)?;
        stage_account(
            &mut self.network,
            EventNativeAccount {
                sample,
                observer_pressure_pa: 0.0,
                energy_before_j: before,
                modal_energy_j: after,
                component_dissipation_j: 0.0,
                connection_energy_j: 0.0,
                connection_dissipation_j: 0.0,
                work_j: work,
                residual_j: native_residual,
                tolerance_j: tolerance,
            },
            [],
        );
        let accepted = JkrNormalEventModalFrame {
            sample,
            network_energy_j: after,
            contact,
            external_work_j: motion.external_work_j,
            normal_work_on_network_j: motion.reaction_work_j,
            released_energy_j: path.released_energy_j,
            normal_work_residual_j: normal_residual,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            coordinate_residual_m: coordinate_residual,
            velocity_residual_m_s: velocity_residual,
            path,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok(accepted)
    }

    fn native_contact(
        &self,
        models: &[ModalAcousticTimeModel],
        phase: JkrBranch,
    ) -> Result<JkrResponse, JkrNormalEventModalError> {
        let mut delta = self.initial_indentation_m;
        for ((state, reference), weight) in models
            .iter()
            .flat_map(|m| m.states())
            .zip(&self.reference_positions)
            .zip(&self.column)
        {
            delta = finite(delta + weight * (state.displacement_m_sqrt_kg - reference))?;
        }
        match phase {
            JkrBranch::Attached => Ok(self.normal().evaluate_attached(delta)?),
            JkrBranch::Detached => {
                let state = self.normal().detached_state(delta)?;
                // Stay on the already selected detached phase. The scalar
                // event owner alone decides attachment, including exact touch.
                Ok(self.motion.response(self.motion.initialize(state, 0.0)?)?)
            }
        }
    }
}

fn check_energy(residual_j: f64, tolerance_j: f64) -> Result<(), JkrNormalEventModalError> {
    if residual_j.abs() > tolerance_j {
        return Err(ModalCouplingError::EnergyBalance {
            residual_j,
            tolerance_j,
        }
        .into());
    }
    Ok(())
}
fn compatible(
    what: &'static str,
    actual: f64,
    expected: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, JkrNormalEventModalError> {
    let residual = finite(actual - expected)?;
    let tolerance = finite(absolute + relative * actual.abs().max(expected.abs()))?;
    if residual.abs() > tolerance {
        return Err(JkrNormalEventModalError::Compatibility {
            what,
            residual,
            tolerance,
        });
    }
    Ok(residual)
}
