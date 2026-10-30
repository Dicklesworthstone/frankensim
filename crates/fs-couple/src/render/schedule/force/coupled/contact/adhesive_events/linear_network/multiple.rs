//! Joint event resolution for independent scalar junctions on one elastic body network.
//!
//! All active static reactions solve the same Gram system, and their mass-orthogonal
//! projection constrains only the corresponding contact coordinates. Every original
//! mode, spring and dashpot remains in the shared matrix-free trajectory owner.
//! Sliding efforts are held kinetic capacities; simultaneous stop/release events
//! are resolved together within the declared event-time accuracy. No independent
//! contact solve, velocity smoothing or artificial contact stiffness is introduced.
//!
//! At zero-velocity transitions the owner retains the largest feasible sticking
//! set. Equally sized admissible alternatives refuse rather than selecting an
//! arbitrary unloading order. Full independent columns are required; redundant
//! contact force distributions, vector Coulomb cones, evolving capacities and
//! general nonlinear mechanics require other owners. Numerical work uncertainty
//! is an explicit estimate, not an interval certificate of every event.

use super::{
    Error, LinearNetworkEventConfig, LinearNetworkEventFlow, LinearNetworkEventKind,
    LinearNetworkLink, LinearNetworkPhase, MAX_LINEAR_NETWORK_EVENTS, ORDER, apply, dot, finite,
    poll, roots, sign,
};
use fs_exec::CancelGate;
use fs_math::det;

#[path = "multiple/contact.rs"]
mod contact;
pub(super) use contact::Contacts;

pub const MAX_LINEAR_NETWORK_CONTACTS: usize = 8;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleLinearNetworkEventConfig {
    pub network: LinearNetworkEventConfig,
    pub maximum_contacts: usize,
    /// Minimum Schur pivot of the Gram matrix of UNIT work columns, in (0,1).
    pub minimum_column_pivot: f64,
    /// Aggregate admissible-phase trials over one complete advance, not per event.
    pub maximum_phase_trials: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleLinearNetworkEvent {
    pub contact_index: usize,
    pub time_s: f64,
    pub kind: LinearNetworkEventKind,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleLinearNetworkEventStep {
    pub contact_count: usize,
    pub phases: [LinearNetworkPhase; MAX_LINEAR_NETWORK_CONTACTS],
    pub initial_energy_j: f64,
    pub modal_energy_j: f64,
    pub connection_energy_j: f64,
    pub external_work_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub contact_friction_dissipation_j: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub accumulated_slip_m: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub contact_impulse_n_s: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub contact_time_moment_n_s2: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub events: [Option<MultipleLinearNetworkEvent>; MAX_LINEAR_NETWORK_EVENTS],
    pub event_count: usize,
    pub cells: usize,
    pub operator_evaluations: usize,
    pub root_evaluations: usize,
    pub phase_trials: usize,
}

/// Shared physical cell account; both scalar and multiple drivers use it.
pub(super) struct ContactAccount {
    pub modal_loss: f64,
    pub friction_loss: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub slip: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub impulse: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub moment: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub uncertainty: f64,
}
impl Default for ContactAccount {
    fn default() -> Self {
        Self {
            modal_loss: 0.0,
            friction_loss: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            slip: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            impulse: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            moment: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            uncertainty: 0.0,
        }
    }
}

/// Preallocated joint-contact workspace. Accepted state/phase belong to the caller.
pub struct MultipleLinearNetworkEventFlow {
    flow: LinearNetworkEventFlow,
    config: MultipleLinearNetworkEventConfig,
}

impl MultipleLinearNetworkEventFlow {
    pub fn new(
        modal_stiffness: Vec<f64>,
        modal_damping: Vec<f64>,
        links: Vec<LinearNetworkLink>,
        columns: Vec<Vec<f64>>,
        config: MultipleLinearNetworkEventConfig,
    ) -> Result<Self, Error> {
        let count = columns.len();
        if config.maximum_contacts == 0
            || config.maximum_contacts > MAX_LINEAR_NETWORK_CONTACTS
            || count == 0
            || count > config.maximum_contacts
            || config.maximum_phase_trials == 0
        {
            return Err(Error::Budget("joint contact extent and phase trials"));
        }
        let n = modal_stiffness.len();
        let entries = n
            .checked_mul(2 * (ORDER + 1) + 8 + count + links.len())
            .and_then(|v| v.checked_add(3 * links.len()))
            .ok_or(Error::Budget("joint contact workspace extent"))?;
        let bytes = entries
            .checked_mul(8)
            .and_then(|v| v.checked_add(core::mem::size_of::<LinearNetworkEventFlow>()))
            .ok_or(Error::Budget("joint contact workspace bytes"))?;
        if bytes > config.network.maximum_workspace_bytes {
            return Err(Error::Budget("joint contact workspace bytes"));
        }
        let contacts = Contacts::new(columns, config.minimum_column_pivot)?;
        let mut flow = LinearNetworkEventFlow::new(
            modal_stiffness,
            modal_damping,
            links,
            contacts.columns[0].clone(),
            config.network,
        )?;
        flow.contacts = contacts;
        Ok(Self { flow, config })
    }

    pub fn contact_count(&self) -> usize {
        self.flow.contacts.count()
    }
    pub fn positions(&self) -> &[f64] {
        self.flow.positions()
    }
    pub fn velocities(&self) -> &[f64] {
        self.flow.velocities()
    }
    pub fn connection_dissipation_j(&self) -> &[f64] {
        self.flow.connection_dissipation_j()
    }
    pub fn connection_forces_n(&self) -> &[f64] {
        self.flow.connection_forces_n()
    }
    pub fn connection_force_bounds_n(&self) -> &[f64] {
        self.flow.connection_force_bounds_n()
    }
    pub const fn config(&self) -> MultipleLinearNetworkEventConfig {
        self.config
    }

    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        q0: &[f64],
        v0: &[f64],
        external: &[f64],
        previous_phases: &[LinearNetworkPhase],
        dt_s: f64,
        static_capacities_n: &[f64],
        kinetic_capacities_n: &[f64],
        gate: &CancelGate,
    ) -> Result<MultipleLinearNetworkEventStep, Error> {
        poll(gate)?;
        let n = self.flow.q.len();
        let count = self.contact_count();
        let cfg = self.config.network;
        if q0.len() != n
            || v0.len() != n
            || external.len() != n
            || previous_phases.len() != count
            || static_capacities_n.len() != count
            || kinetic_capacities_n.len() != count
            || q0.iter().chain(v0).chain(external).any(|x| !x.is_finite())
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || dt_s > cfg.maximum_step_s
            || static_capacities_n
                .iter()
                .zip(kinetic_capacities_n)
                .any(|(s, k)| !s.is_finite() || !k.is_finite() || *k < 0.0 || *s < *k)
        {
            return Err(Error::Invalid(
                "complete joint initial state, capacities and held load",
            ));
        }
        self.flow.q.copy_from_slice(q0);
        self.flow.v.copy_from_slice(v0);
        self.flow.link_losses.fill(0.0);
        self.flow.link_forces.fill(0.0);
        self.flow.link_force_bounds.fill(0.0);
        self.flow.operator_evaluations = 0;
        self.flow.roots.reset();
        self.flow.check_state()?;
        let (em0, el0) = self.flow.energy(q0, v0)?;
        let mut step = MultipleLinearNetworkEventStep {
            contact_count: count,
            phases: [LinearNetworkPhase::Fresh; MAX_LINEAR_NETWORK_CONTACTS],
            initial_energy_j: finite(em0 + el0)?,
            modal_energy_j: em0,
            connection_energy_j: el0,
            external_work_j: 0.0,
            modal_damping_dissipation_j: 0.0,
            connection_damping_dissipation_j: 0.0,
            friction_dissipation_j: 0.0,
            contact_friction_dissipation_j: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            accumulated_slip_m: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            contact_impulse_n_s: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            contact_time_moment_n_s2: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            energy_residual_j: 0.0,
            energy_tolerance_j: 0.0,
            numerical_work_uncertainty_j: 0.0,
            events: [None; MAX_LINEAR_NETWORK_EVENTS],
            event_count: 0,
            cells: 0,
            operator_evaluations: 0,
            root_evaluations: 0,
            phase_trials: 0,
        };
        let mut account = ContactAccount::default();
        for i in 0..count {
            self.flow.contacts.static_capacity[i] = static_capacities_n[i];
            self.flow.contacts.kinetic_capacity[i] = kinetic_capacities_n[i];
            let speed = dot(&self.flow.contacts.columns[i], v0)?;
            self.flow.contacts.phases[i] = match previous_phases[i] {
                LinearNetworkPhase::Fresh => {
                    if speed == 0.0 {
                        LinearNetworkPhase::Sticking
                    } else {
                        LinearNetworkPhase::Sliding {
                            direction: sign(speed),
                        }
                    }
                }
                LinearNetworkPhase::Sticking => {
                    if speed.abs() > cfg.sticking_velocity_tolerance_m_s {
                        return Err(Error::PhaseMismatch);
                    }
                    LinearNetworkPhase::Sticking
                }
                LinearNetworkPhase::Sliding { direction } => {
                    if direction.abs() != 1 || speed != 0.0 && sign(speed) != direction {
                        return Err(Error::PhaseMismatch);
                    }
                    // A sample boundary may change the load exactly at release.
                    if speed == 0.0 {
                        LinearNetworkPhase::Sticking
                    } else {
                        previous_phases[i]
                    }
                }
            };
        }
        self.flow.contacts.factor()?;
        self.project_velocity(&mut account)?;
        self.select_phase(external, 0.0, &mut step, gate)?;
        let mut time = 0.0;
        while time < dt_s {
            poll(gate)?;
            if step.cells >= cfg.maximum_cells {
                return Err(Error::Budget("joint trajectory cells"));
            }
            let remainder = finite(dt_s - time)?;
            self.flow.prepare_contacts_cell(external, remainder, gate)?;
            step.cells += 1;
            let h = self.flow.cell_h;
            let tolerance = (cfg.event_time_tolerance_s / h)
                .min(1e-3)
                .max(8.0 * f64::EPSILON);
            let mut events: [Option<roots::EventRoot>; MAX_LINEAR_NETWORK_CONTACTS] =
                [None; MAX_LINEAR_NETWORK_CONTACTS];
            let mut first: f64 = 1.0;
            let mut found = false;
            for i in 0..count {
                let root = match self.flow.contacts.phases[i] {
                    LinearNetworkPhase::Sliding { direction } => roots::first_zero(
                        &self.flow.contacts.velocity_coeff[i],
                        direction,
                        tolerance,
                        &mut self.flow.roots,
                        gate,
                    )?,
                    LinearNetworkPhase::Sticking => roots::first_exit(
                        &self.flow.contacts.reaction_coeff[i],
                        static_capacities_n[i],
                        tolerance,
                        &mut self.flow.roots,
                        gate,
                    )?,
                    LinearNetworkPhase::Fresh => return Err(Error::PhaseMismatch),
                };
                if let Some(root) = root {
                    first = first.min(root.fraction);
                    found = true;
                }
                events[i] = root;
            }
            self.flow
                .integrate_contacts_cell(first, time, external, &mut account, gate)?;
            let next_time = if !found && h == remainder {
                dt_s
            } else {
                finite(time + h * first)?
            };
            if found {
                let mut timing_impulse = 0.0;
                for i in 0..count {
                    let Some(event) = events[i] else {
                        continue;
                    };
                    if (event.fraction - first) * h > cfg.event_time_tolerance_s {
                        continue;
                    }
                    timing_impulse += (static_capacities_n[i] + kinetic_capacities_n[i])
                        * cfg.event_time_tolerance_s
                        * self.flow.contacts.columns[i]
                            .iter()
                            .map(|b| b.abs())
                            .fold(0.0, f64::max);
                    match self.flow.contacts.phases[i] {
                        LinearNetworkPhase::Sliding { .. } => {
                            self.flow.contacts.phases[i] = LinearNetworkPhase::Sticking;
                            self.event(&mut step, i, next_time, LinearNetworkEventKind::Stopped)?;
                        }
                        LinearNetworkPhase::Sticking => {
                            self.flow.contacts.phases[i] = LinearNetworkPhase::Sliding {
                                direction: event.direction,
                            };
                            self.event(
                                &mut step,
                                i,
                                next_time,
                                LinearNetworkEventKind::SlipStarted {
                                    direction: event.direction,
                                },
                            )?;
                        }
                        _ => return Err(Error::PhaseMismatch),
                    }
                }
                self.timing_uncertainty(external, timing_impulse, &mut account)?;
                self.flow.contacts.factor()?;
                self.project_velocity(&mut account)?;
                self.select_phase(external, next_time, &mut step, gate)?;
            } else if next_time == time {
                return Err(Error::EventUnresolved);
            }
            time = next_time;
        }
        let (em, el) = self.flow.energy(&self.flow.q, &self.flow.v)?;
        step.phases = self.flow.contacts.phases;
        step.modal_energy_j = em;
        step.connection_energy_j = el;
        step.external_work_j = finite(
            external
                .iter()
                .zip(self.flow.q.iter().zip(q0))
                .map(|(f, (q, q0))| f * (q - q0))
                .sum(),
        )?;
        step.modal_damping_dissipation_j = account.modal_loss;
        step.connection_damping_dissipation_j = finite(self.flow.link_losses.iter().sum())?;
        step.contact_friction_dissipation_j = account.friction_loss;
        step.accumulated_slip_m = account.slip;
        step.contact_impulse_n_s = account.impulse;
        step.contact_time_moment_n_s2 = account.moment;
        step.friction_dissipation_j = finite(account.friction_loss[..count].iter().sum())?;
        step.numerical_work_uncertainty_j = account.uncertainty;
        for (i, link) in self.flow.links.iter().enumerate() {
            self.flow.link_forces[i] = finite(
                -link.stiffness_n_per_m
                    * (dot(&link.column, &self.flow.q)? - link.rest_extension_m)
                    - link.damping_n_s_per_m * dot(&link.column, &self.flow.v)?,
            )?;
        }
        step.energy_residual_j = finite(
            em + el - step.initial_energy_j - step.external_work_j
                + step.modal_damping_dissipation_j
                + step.connection_damping_dissipation_j
                + step.friction_dissipation_j,
        )?;
        let scale = step
            .initial_energy_j
            .abs()
            .max((em + el).abs())
            .max(step.external_work_j.abs())
            .max(
                step.modal_damping_dissipation_j
                    + step.connection_damping_dissipation_j
                    + step.friction_dissipation_j,
            );
        step.energy_tolerance_j =
            finite(cfg.energy_absolute_tolerance_j + cfg.energy_relative_tolerance * scale)?;
        step.operator_evaluations = self.flow.operator_evaluations;
        step.root_evaluations = self.flow.roots.used();
        if step.numerical_work_uncertainty_j > cfg.maximum_numerical_work_uncertainty_j
            || step.energy_residual_j.abs() + step.numerical_work_uncertainty_j
                > step.energy_tolerance_j
        {
            return Err(Error::EnergyAccuracy {
                residual_j: step.energy_residual_j,
                estimated_error_j: step.numerical_work_uncertainty_j,
                tolerance_j: step
                    .energy_tolerance_j
                    .min(cfg.maximum_numerical_work_uncertainty_j),
            });
        }
        poll(gate)?;
        Ok(step)
    }

    fn event(
        &self,
        step: &mut MultipleLinearNetworkEventStep,
        contact_index: usize,
        time_s: f64,
        kind: LinearNetworkEventKind,
    ) -> Result<(), Error> {
        if step.event_count >= self.config.network.maximum_events {
            return Err(Error::Budget("joint physical events"));
        }
        step.events[step.event_count] = Some(MultipleLinearNetworkEvent {
            contact_index,
            time_s,
            kind,
        });
        step.event_count += 1;
        Ok(())
    }

    fn project_velocity(&mut self, account: &mut ContactAccount) -> Result<(), Error> {
        let norm = det::sqrt(dot(&self.flow.v, &self.flow.v)?);
        let (_, correction) = self.flow.contacts.project(&mut self.flow.v, None)?;
        account.uncertainty =
            finite(account.uncertainty + norm * correction + 0.5 * correction * correction)?;
        Ok(())
    }

    fn timing_uncertainty(
        &self,
        external: &[f64],
        velocity_error: f64,
        account: &mut ContactAccount,
    ) -> Result<(), Error> {
        let n = self.flow.q.len() as f64;
        let qmax = self.flow.q.iter().map(|q| q.abs()).fold(0.0, f64::max);
        let vmax = self.flow.v.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let rest = self
            .flow
            .links
            .iter()
            .map(|l| {
                l.stiffness_n_per_m
                    * l.rest_extension_m.abs()
                    * l.column.iter().map(|b| b.abs()).sum::<f64>()
            })
            .sum::<f64>();
        let friction = (0..self.contact_count())
            .map(|i| self.flow.contacts.kinetic_capacity[i] * self.flow.contacts.norm1[i])
            .sum::<f64>();
        let force = external.iter().map(|f| f.abs()).sum::<f64>()
            + n * self.flow.stiffness_norm * qmax
            + rest
            + friction;
        account.uncertainty = finite(
            account.uncertainty
                + n * (vmax * velocity_error + 0.5 * velocity_error * velocity_error)
                + velocity_error * self.config.network.event_time_tolerance_s * force,
        )?;
        Ok(())
    }

    fn select_phase(
        &mut self,
        external: &[f64],
        time: f64,
        step: &mut MultipleLinearNetworkEventStep,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let before = self.flow.contacts.phases;
        let count = self.contact_count();
        let mut candidates = [0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut candidate_count = 0;
        for i in 0..count {
            if before[i] == LinearNetworkPhase::Sticking {
                candidates[candidate_count] = i;
                candidate_count += 1;
            }
        }
        apply(
            &self.flow.stiffness,
            &self.flow.damping,
            &self.flow.links,
            &self.flow.q,
            &self.flow.v,
            Some(external),
            &mut self.flow.accel,
        )?;
        self.flow.operator_evaluations += 1;
        if self.flow.operator_evaluations > self.config.network.maximum_operator_evaluations {
            return Err(Error::Budget("joint operator evaluations"));
        }
        let choices = 3usize.pow(candidate_count as u32);
        let mut best = None;
        let mut best_sticks = 0;
        let mut ambiguous = false;
        for code in 0..choices {
            poll(gate)?;
            let mut digits = code;
            let mut phases = before;
            let mut sticks = 0;
            for &i in &candidates[..candidate_count] {
                phases[i] = match digits % 3 {
                    0 => {
                        sticks += 1;
                        LinearNetworkPhase::Sticking
                    }
                    1 => LinearNetworkPhase::Sliding { direction: 1 },
                    _ => LinearNetworkPhase::Sliding { direction: -1 },
                };
                digits /= 3;
            }
            if best.is_some() && sticks < best_sticks {
                continue;
            }
            if step.phase_trials >= self.config.maximum_phase_trials {
                return Err(Error::Budget("joint admissible phase trials"));
            }
            step.phase_trials += 1;
            self.flow.contacts.phases = phases;
            self.flow.contacts.factor()?;
            self.flow.node_v.copy_from_slice(&self.flow.accel);
            let reaction = self.flow.contacts.reactions(&mut self.flow.node_v, true)?;
            let mut admissible = true;
            for &i in &candidates[..candidate_count] {
                match phases[i] {
                    LinearNetworkPhase::Sticking => {
                        if reaction[i].abs() > self.flow.contacts.static_capacity[i] {
                            admissible = false;
                            break;
                        }
                    }
                    LinearNetworkPhase::Sliding { direction } => {
                        let acceleration = dot(&self.flow.contacts.columns[i], &self.flow.node_v)?;
                        // A newly chosen kinetic branch must actually depart in
                        // that direction. An unresolved grazing choice refuses.
                        if direction as f64 * acceleration <= 0.0 {
                            admissible = false;
                            break;
                        }
                    }
                    _ => return Err(Error::PhaseMismatch),
                }
            }
            if !admissible {
                continue;
            }
            if code == 0 {
                best = Some(phases);
                ambiguous = false;
                break;
            }
            if best.is_none() || sticks > best_sticks {
                best = Some(phases);
                best_sticks = sticks;
                ambiguous = false;
            } else {
                ambiguous = true;
            }
        }
        if ambiguous {
            return Err(Error::EventUnresolved);
        }
        self.flow.contacts.phases = best.ok_or(Error::EventUnresolved)?;
        self.flow.contacts.factor()?;
        for i in 0..count {
            if before[i] == LinearNetworkPhase::Sticking {
                if let LinearNetworkPhase::Sliding { direction } = self.flow.contacts.phases[i] {
                    self.event(
                        step,
                        i,
                        time,
                        LinearNetworkEventKind::SlipStarted { direction },
                    )?;
                }
            }
        }
        Ok(())
    }
}
