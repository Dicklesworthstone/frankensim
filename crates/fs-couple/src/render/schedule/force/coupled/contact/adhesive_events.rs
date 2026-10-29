//! Exact stop/restart junction impulses on the original free modal bodies.
//!
//! This lane admits only mass-normalized free coordinates, no connections,
//! material memory or radiation, and one fixed tangential work column `B`.
//! Normal approach is an explicit prescribed actuator, linear during a sample.
//! The original sourced DMT/constant-junction event owner determines every
//! stop, static hold and restart. Its impulse and first time moment propagate
//! all original body coordinates, including common and uncontacted motion.
//! Elastic modes, rotating contact frames and two-dimensional slip require a
//! different trajectory solver and are refused here. No velocity is clamped.

use super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::events::{
    DmtSlipEventBudget, DmtSlipEventError, DmtSlipEventState, DmtSlipEventStep, DmtSlipEvents,
};
use fs_contact::normal_patch::dmt::{DmtError, SourceDmtSphereContact};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;

#[derive(Debug)]
pub enum AdhesiveEventModalError {
    Network(ModalCouplingError),
    Contact(DmtSlipEventError),
    Compatibility {
        what: &'static str,
        residual: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for AdhesiveEventModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native adhesive event trial refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveEventModalError {}
impl From<ModalCouplingError> for AdhesiveEventModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}
impl From<DmtSlipEventError> for AdhesiveEventModalError {
    fn from(error: DmtSlipEventError) -> Self {
        Self::Contact(error)
    }
}
impl From<DmtError> for AdhesiveEventModalError {
    fn from(error: DmtError) -> Self {
        Self::Contact(error.into())
    }
}

/// Numerical agreement of original native motion with the analytic event path.
/// These tolerances never select sticking or set an endpoint velocity to zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveEventModalConfig {
    pub events: DmtSlipEventBudget,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub velocity_relative_tolerance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveEventModalFrame {
    pub sample: u64,
    /// Original free coordinates have zero acoustic pressure transfers.
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    pub prescribed_normal_work_j: f64,
    /// Signed work at the boundary of the original mechanical subsystem.
    pub friction_work_on_network_j: f64,
    /// Analytic sourced friction loss; no thermal destination is inferred.
    pub friction_dissipation_j: f64,
    pub friction_work_residual_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub coordinate_residual_m: f64,
    pub velocity_residual_m_s: f64,
    /// Event times, areas, total absolute slip, impulse and first time moment.
    pub path: DmtSlipEventStep,
}

pub struct AdhesiveEventModalSystem<'a> {
    network: CoupledModalSystem,
    events: DmtSlipEvents<'a>,
    state: DmtSlipEventState,
    column: Vec<f64>,
    candidate_states: Vec<ModalAcousticState>,
    config: AdhesiveEventModalConfig,
    last: Option<AdhesiveEventModalFrame>,
}

impl<'a> AdhesiveEventModalSystem<'a> {
    /// Attach supplied physical coordinates to the current sample-zero geometry.
    /// The column contains every component/mode in original network order:
    /// `dx = B·dq`, `v = B·qdot`, and the reciprocal resisting force is `-B F`.
    /// Consequently the exact relative mass is `1/(B·B)`. Normal motion belongs
    /// to the prescribed actuator and is not an unmodeled finite normal mass.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        normal: &'a SourceDmtSphereContact,
        shear: &'a SourceAdhesiveShear,
        surface_materials: [&MaterialStateId; 2],
        initial_indentation_m: f64,
        initial_tangent_coordinate_m: f64,
        column: Vec<f64>,
        config: AdhesiveEventModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, AdhesiveEventModalError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0
            || !network.connections.is_empty()
            || !network.relaxation_branches.is_empty()
            || !network.plane_stress_material_points.is_empty()
            || network.radiation.is_some()
            || network.models.iter().flat_map(|m| m.modes()).any(|mode| {
                mode.angular_frequency_rad_s != 0.0
                    || mode.damping_ratio != 0.0
                    || mode.pressure_per_modal_velocity != fs_math::c64::C64::ZERO
            })
        {
            return Err(invalid("adhesive events require a sample-zero free network without connections, memory or radiation").into());
        }
        let materials = normal
            .normal_state()
            .elastic()
            .interface()
            .surface_materials();
        if surface_materials[0] != &materials[0] || surface_materials[1] != &materials[1] {
            return Err(invalid(
                "adhesive event receiving surfaces differ from the ordered source",
            )
            .into());
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
            return Err(invalid(
                "complete finite adhesive event work column and positive compatibility budgets",
            )
            .into());
        }
        let events = DmtSlipEvents::new(
            normal,
            shear,
            finite(1.0 / dot(&column, &column)?)?,
            config.events,
        )?;
        let state = events.initialize(
            initial_indentation_m,
            initial_tangent_coordinate_m,
            projected_velocity(network.components(), &column)?,
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
            "native/event initial energy magnitude",
            system.total_energy_j()?.abs(),
            system.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(system)
    }

    /// Original accepted native owner, including observer, clock and port APIs.
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
    pub const fn state(&self) -> DmtSlipEventState {
        self.state
    }
    pub const fn normal(&self) -> &'a SourceDmtSphereContact {
        self.events.normal()
    }
    pub const fn shear(&self) -> &'a SourceAdhesiveShear {
        self.events.shear()
    }
    pub fn work_column(&self) -> &[f64] {
        &self.column
    }
    pub const fn effective_mass_kg(&self) -> f64 {
        self.events.mass_kg()
    }
    pub const fn config(&self) -> AdhesiveEventModalConfig {
        self.config
    }
    pub const fn last_frame(&self) -> Option<&AdhesiveEventModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, AdhesiveEventModalError> {
        Ok(finite(
            self.network.total_energy_j()?
                + self
                    .normal()
                    .evaluate(self.state.indentation_m())?
                    .reversible_energy_j,
        )?)
    }

    /// One sample of held external generalized forces and linear normal approach.
    /// The preallocated original candidate bodies and all event/port history
    /// publish together only after source, state, work and cancellation gates.
    pub fn step(
        &mut self,
        external: &[f64],
        next_indentation_m: f64,
    ) -> Result<&AdhesiveEventModalFrame, AdhesiveEventModalError> {
        self.step_under_gate(external, next_indentation_m, &CancelGate::new_clock_free())
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        next_indentation_m: f64,
        gate: &CancelGate,
    ) -> Result<&AdhesiveEventModalFrame, AdhesiveEventModalError> {
        poll(Some(gate))?;
        if external.len() != self.mode_count() || external.iter().any(|x| !x.is_finite()) {
            return Err(
                invalid("external forces must be finite and match every original mode").into(),
            );
        }
        let sample = self
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| invalid("coupled sample clock overflow"))?;
        let dt = self.network.dt;
        let before = self.network.total_energy_j()?;
        let surface_before = self
            .normal()
            .evaluate(self.state.indentation_m())?
            .reversible_energy_j;
        let v0 = projected_velocity(self.components(), &self.column)?;
        compatible(
            "initial relative velocity",
            v0,
            self.state.tangent_velocity_m_s(),
            self.config.velocity_absolute_tolerance_m_s,
            self.config.velocity_relative_tolerance,
        )?;
        let acceleration = dot(&self.column, external)?;
        let path = self.events.advance(
            self.state,
            next_indentation_m,
            finite(self.effective_mass_kg() * acceleration)?,
            dt,
            gate,
        )?;
        let impulse = path.friction_impulse_n_s;
        let moment = path.friction_time_moment_n_s2;
        let position_impulse = finite(dt * impulse - moment)?;
        let mut after = 0.0;
        for (k, old) in self
            .network
            .models
            .iter()
            .flat_map(|m| m.states())
            .enumerate()
        {
            let increment = finite(
                dt * old.velocity_m_sqrt_kg_per_s + 0.5 * dt * dt * external[k]
                    - self.column[k] * position_impulse,
            )?;
            let next = ModalAcousticState {
                displacement_m_sqrt_kg: finite(old.displacement_m_sqrt_kg + increment)?,
                velocity_m_sqrt_kg_per_s: finite(
                    old.velocity_m_sqrt_kg_per_s + dt * external[k] - self.column[k] * impulse,
                )?,
            };
            // Recompute actual represented displacement, rather than trusting
            // an increment that may be lost when added to a large coordinate.
            self.network.old_q[k] = old.displacement_m_sqrt_kg;
            self.network.free_delta[k] =
                finite(next.displacement_m_sqrt_kg - old.displacement_m_sqrt_kg)?;
            after = finite(after + 0.5 * next.velocity_m_sqrt_kg_per_s.powi(2))?;
            self.candidate_states[k] = next;
        }
        for i in 0..self.network.models.len() {
            poll(Some(gate))?;
            self.network.candidates[i]
                .restore_states(
                    &self.candidate_states[self.network.offsets[i]..self.network.offsets[i + 1]],
                )
                .map_err(|source| ModalCouplingError::Component {
                    component: i,
                    source,
                })?;
        }
        let coordinate_residual = compatible(
            "actual tangent increment",
            dot(&self.column, &self.network.free_delta)?,
            path.tangent_increment_m,
            self.config.coordinate_absolute_tolerance_m,
            self.config.coordinate_relative_tolerance,
        )?;
        let velocity_residual = compatible(
            "actual endpoint relative velocity",
            projected_velocity(&self.network.candidates, &self.column)?,
            path.state().tangent_velocity_m_s(),
            self.config.velocity_absolute_tolerance_m_s,
            self.config.velocity_relative_tolerance,
        )?;
        limit(
            "native event kinetic energy",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let surface_after = path.contact.reversible_energy_j;
        limit(
            "native/event energy magnitude",
            finite(after + surface_after)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let external_work = dot(external, &self.network.free_delta)?;
        // Independent exact work on the actual native trajectory:
        // v_rel(t)=v0+a*t-(B·B)*I(t), ∫F I dt=J²/2.
        // This retains common motion in the body energy and never infers work
        // from a zero net slip or manufactures a residual-energy correction.
        let friction_work = finite(
            -v0 * impulse - acceleration * moment
                + 0.5 * dot(&self.column, &self.column)? * impulse * impulse,
        )?;
        let mechanical_work = finite(external_work + friction_work)?;
        let native_residual = finite(after - before - mechanical_work)?;
        let friction_residual = finite(friction_work + path.friction_dissipation_j)?;
        let residual = finite(
            (after - before) + (surface_after - surface_before) + path.friction_dissipation_j
                - external_work
                - path.prescribed_normal_work_j,
        )?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(surface_before.abs())
            .max(surface_after.abs())
            .max(external_work.abs())
            .max(friction_work.abs())
            .max(path.friction_dissipation_j)
            .max(path.prescribed_normal_work_j.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        for residual_j in [native_residual, friction_residual, residual] {
            if residual_j.abs() > tolerance {
                return Err(ModalCouplingError::EnergyBalance {
                    residual_j,
                    tolerance_j: tolerance,
                }
                .into());
            }
        }
        let frame = &mut self.network.staged_frame;
        frame.sample = sample;
        frame.observer_pressure_pa = 0.0;
        frame.modal_energy_j = after;
        frame.connection_energy_j = 0.0;
        frame.material_energy_j = 0.0;
        frame.radiation_energy_j = 0.0;
        frame.external_work_j = mechanical_work;
        frame.component_dissipation_j = 0.0;
        frame.connection_dissipation_j = 0.0;
        frame.material_dissipation_j = 0.0;
        frame.radiation_dissipation_j = 0.0;
        frame.energy_residual_j = native_residual;
        frame.energy_tolerance_j = tolerance;
        frame.solve_relative_residual = 0.0;
        self.network.staged_ports.accept(
            crate::port_adoption::StepAccount {
                tick: sample,
                duration_s: dt,
                storage_change_j: Some(after - before),
                loss_j: Some(0.0),
                included_work_j: mechanical_work,
                entropy_change_j_k: 0.0,
                tolerance_j: Some(tolerance),
                unavailable: None,
            },
            [],
        );
        let accepted = AdhesiveEventModalFrame {
            sample,
            observer_pressure_pa: 0.0,
            network_energy_j: after,
            surface_energy_j: surface_after,
            external_work_j: external_work,
            prescribed_normal_work_j: path.prescribed_normal_work_j,
            friction_work_on_network_j: friction_work,
            friction_dissipation_j: path.friction_dissipation_j,
            friction_work_residual_j: friction_residual,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            coordinate_residual_m: coordinate_residual,
            velocity_residual_m_s: velocity_residual,
            path,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        self.state = path.state();
        Ok(self.last.insert(accepted))
    }
}

fn projected_velocity(
    models: &[ModalAcousticTimeModel],
    column: &[f64],
) -> Result<f64, ModalCouplingError> {
    models
        .iter()
        .flat_map(|m| m.states())
        .zip(column)
        .try_fold(0.0, |v, (state, b)| {
            finite(v + b * state.velocity_m_sqrt_kg_per_s)
        })
}

fn compatible(
    what: &'static str,
    actual: f64,
    expected: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, AdhesiveEventModalError> {
    let residual = finite(actual - expected)?;
    let tolerance = finite(absolute + relative * actual.abs().max(expected.abs()))?;
    if residual.abs() > tolerance {
        return Err(AdhesiveEventModalError::Compatibility {
            what,
            residual,
            tolerance,
        });
    }
    Ok(residual)
}
