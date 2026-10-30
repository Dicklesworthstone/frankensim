//! Resolved scalar friction on the complete passive linear mechanical network.
//!
//! Coordinates remain the owner's mass-normalized coordinates. The operator is
//! diagonal modal stiffness/resistance plus every original rank-one spring and
//! Kelvin--Voigt dashpot. Sliding has constant signed friction; sticking uses
//! the exact mass-orthogonal constraint projection. No mode or link is removed.
//!
//! Nonstationary elastic Taylor cells have scaled operator norm at most 1/4;
//! exactly stationary and nilpotent free-body cells terminate analytically.
//! Degree-24 derivative
//! root isolation partitions the dense contact velocity/reaction polynomial
//! into monotone intervals before selecting its first event. Positive Gaussian
//! integration measures each actual resistance independently. Roundoff, Taylor
//! remainder and stop projection enter an explicit numerical-work estimate;
//! that estimate is not an interval certificate. An unresolved event refuses.
//!
//! All state and physical phase publication belongs to the caller. A successful
//! trial may still be rejected by original component budgets without changing
//! the next trial. Coefficients, contact frame, capacities and external force are
//! held over one call. Changing normal geometry, material memory, nonlinear
//! stiffness, vector slip, and radiative feedback require other owners.

use fs_exec::CancelGate;
use fs_math::det;

#[path = "linear_network/flow.rs"]
mod flow;
#[path = "linear_network/roots.rs"]
mod roots;

/// Joint stopping and release at several independent scalar junctions.
#[path = "linear_network/multiple.rs"]
pub mod multiple;

/// Joint two-dimensional Coulomb disks with separately measured time-discretization loss.
#[path = "linear_network/vector.rs"]
pub mod vector;

const ORDER: usize = 24;
pub const MAX_LINEAR_NETWORK_EVENTS: usize = 64;

#[derive(Clone, Debug, PartialEq)]
pub struct LinearNetworkLink {
    pub column: Vec<f64>,
    pub stiffness_n_per_m: f64,
    pub damping_n_s_per_m: f64,
    pub rest_extension_m: f64,
}

/// Explicit physical phase, retained by the outer physical transaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinearNetworkPhase {
    /// Initial admission only: every nonzero actual velocity is sliding.
    Fresh,
    Sticking,
    Sliding {
        direction: i8,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearNetworkEventConfig {
    pub maximum_modes: usize,
    pub maximum_connections: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_cells: usize,
    pub maximum_events: usize,
    pub maximum_operator_evaluations: usize,
    pub maximum_root_evaluations: usize,
    pub maximum_root_iterations: usize,
    pub maximum_step_s: f64,
    pub event_time_tolerance_s: f64,
    /// Only validates an explicitly retained Sticking phase; never captures a
    /// Fresh or Sliding velocity merely because it is numerically small.
    pub sticking_velocity_tolerance_m_s: f64,
    pub maximum_abs_modal_coordinate: f64,
    pub maximum_abs_modal_velocity: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
    pub maximum_numerical_work_uncertainty_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinearNetworkEventKind {
    Stopped,
    SlipStarted { direction: i8 },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearNetworkEvent {
    pub time_s: f64,
    pub kind: LinearNetworkEventKind,
}

#[derive(Clone, Debug, PartialEq)]
pub enum LinearNetworkEventError {
    Invalid(&'static str),
    Nonfinite,
    Cancelled,
    Budget(&'static str),
    PhaseMismatch,
    EventUnresolved,
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for LinearNetworkEventError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "linear network adhesive event refused: {self:?}")
    }
}
impl std::error::Error for LinearNetworkEventError {}
type Error = LinearNetworkEventError;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearNetworkEventStep {
    pub phase: LinearNetworkPhase,
    pub initial_energy_j: f64,
    pub modal_energy_j: f64,
    pub connection_energy_j: f64,
    pub external_work_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub accumulated_slip_m: f64,
    pub contact_impulse_n_s: f64,
    pub contact_time_moment_n_s2: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub numerical_work_uncertainty_j: f64,
    pub events: [Option<LinearNetworkEvent>; MAX_LINEAR_NETWORK_EVENTS],
    pub event_count: usize,
    pub cells: usize,
    pub operator_evaluations: usize,
    pub root_evaluations: usize,
}

/// Reusable matrix-free trajectory workspace. Allocation occurs at admission.
pub struct LinearNetworkEventFlow {
    stiffness: Vec<f64>,
    damping: Vec<f64>,
    links: Vec<LinearNetworkLink>,
    tangent: Vec<f64>,
    norm_b2: f64,
    norm_b1: f64,
    stiffness_norm: f64,
    damping_norm: f64,
    config: LinearNetworkEventConfig,
    contacts: multiple::Contacts,
    q: Vec<f64>,
    v: Vec<f64>,
    q_coeff: Vec<f64>,
    v_coeff: Vec<f64>,
    accel: Vec<f64>,
    node_q: Vec<f64>,
    node_v: Vec<f64>,
    link_losses: Vec<f64>,
    link_forces: Vec<f64>,
    link_force_bounds: Vec<f64>,
    velocity_coeff: [f64; ORDER + 1],
    reaction_coeff: [f64; ORDER + 1],
    cell_state_error: f64,
    cell_scale: f64,
    cell_rate: f64,
    cell_h: f64,
    operator_evaluations: usize,
    roots: roots::RootBudget,
}

impl LinearNetworkEventFlow {
    pub fn new(
        modal_stiffness: Vec<f64>,
        modal_damping: Vec<f64>,
        links: Vec<LinearNetworkLink>,
        tangent_column: Vec<f64>,
        config: LinearNetworkEventConfig,
    ) -> Result<Self, Error> {
        let n = modal_stiffness.len();
        validate_config(config)?;
        if n == 0 || n > config.maximum_modes || links.len() > config.maximum_connections {
            return Err(Error::Budget("complete original network extent"));
        }
        if modal_damping.len() != n
            || tangent_column.len() != n
            || modal_stiffness
                .iter()
                .chain(&modal_damping)
                .any(|x| !x.is_finite() || *x < 0.0)
            || tangent_column.iter().any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "passive modal diagonal and complete contact column",
            ));
        }
        for link in &links {
            if link.column.len() != n
                || link.column.iter().any(|x| !x.is_finite())
                || !link.stiffness_n_per_m.is_finite()
                || link.stiffness_n_per_m < 0.0
                || !link.damping_n_s_per_m.is_finite()
                || link.damping_n_s_per_m < 0.0
                || !link.rest_extension_m.is_finite()
            {
                return Err(Error::Invalid("passive original connection"));
            }
        }
        let entries = n
            .checked_mul(2 * (ORDER + 1) + 5)
            .and_then(|v| v.checked_add(3 * links.len()))
            .and_then(|v| v.checked_add(4 * (ORDER + 1)))
            .ok_or(Error::Budget("workspace extent"))?;
        if entries
            .checked_add(n)
            .and_then(|entries| entries.checked_mul(8))
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<multiple::Contacts>()))
            .ok_or(Error::Budget("workspace bytes"))?
            > config.maximum_workspace_bytes
        {
            return Err(Error::Budget("workspace bytes"));
        }
        let norm_b2 = dot(&tangent_column, &tangent_column)?;
        if norm_b2 <= 0.0 {
            return Err(Error::Invalid("nonzero contact column"));
        }
        let norm_b1 = finite(tangent_column.iter().map(|b| b.abs()).sum())?;
        let (stiffness_norm, damping_norm) = {
            let mut stiffness_rows = modal_stiffness.clone();
            let mut damping_rows = modal_damping.clone();
            for link in &links {
                let column_norm = finite(link.column.iter().map(|b| b.abs()).sum())?;
                for i in 0..n {
                    let weight = link.column[i].abs() * column_norm;
                    stiffness_rows[i] =
                        finite(stiffness_rows[i] + link.stiffness_n_per_m * weight)?;
                    damping_rows[i] = finite(damping_rows[i] + link.damping_n_s_per_m * weight)?;
                }
            }
            (
                stiffness_rows.into_iter().fold(0.0, f64::max),
                damping_rows.into_iter().fold(0.0, f64::max),
            )
        };
        let l = links.len();
        let contacts = multiple::Contacts::new(vec![tangent_column.clone()], 128.0 * f64::EPSILON)?;
        Ok(Self {
            stiffness: modal_stiffness,
            damping: modal_damping,
            links,
            tangent: tangent_column,
            norm_b2,
            norm_b1,
            stiffness_norm,
            damping_norm,
            config,
            contacts,
            q: vec![0.0; n],
            v: vec![0.0; n],
            q_coeff: vec![0.0; n * (ORDER + 1)],
            v_coeff: vec![0.0; n * (ORDER + 1)],
            accel: vec![0.0; n],
            node_q: vec![0.0; n],
            node_v: vec![0.0; n],
            link_losses: vec![0.0; l],
            link_forces: vec![0.0; l],
            link_force_bounds: vec![0.0; l],
            velocity_coeff: [0.0; ORDER + 1],
            reaction_coeff: [0.0; ORDER + 1],
            cell_state_error: 0.0,
            cell_scale: 0.0,
            cell_rate: 0.0,
            cell_h: 0.0,
            operator_evaluations: 0,
            roots: roots::RootBudget::new(
                config.maximum_root_evaluations,
                config.maximum_root_iterations,
            ),
        })
    }

    pub fn positions(&self) -> &[f64] {
        &self.q
    }
    pub fn velocities(&self) -> &[f64] {
        &self.v
    }
    pub fn connection_dissipation_j(&self) -> &[f64] {
        &self.link_losses
    }
    pub fn connection_forces_n(&self) -> &[f64] {
        &self.link_forces
    }
    pub fn connection_force_bounds_n(&self) -> &[f64] {
        &self.link_force_bounds
    }
    pub const fn config(&self) -> LinearNetworkEventConfig {
        self.config
    }

    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        q0: &[f64],
        v0: &[f64],
        external: &[f64],
        previous_phase: LinearNetworkPhase,
        dt_s: f64,
        static_capacity_n: f64,
        kinetic_capacity_n: f64,
        gate: &CancelGate,
    ) -> Result<LinearNetworkEventStep, Error> {
        poll(gate)?;
        let n = self.q.len();
        if q0.len() != n
            || v0.len() != n
            || external.len() != n
            || q0.iter().chain(v0).chain(external).any(|x| !x.is_finite())
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || dt_s > self.config.maximum_step_s
            || !static_capacity_n.is_finite()
            || static_capacity_n < 0.0
            || !kinetic_capacity_n.is_finite()
            || kinetic_capacity_n < 0.0
            || kinetic_capacity_n > static_capacity_n
        {
            return Err(Error::Invalid(
                "complete initial state, held load, duration and capacities",
            ));
        }
        self.q.copy_from_slice(q0);
        self.v.copy_from_slice(v0);
        self.link_losses.fill(0.0);
        self.link_forces.fill(0.0);
        self.link_force_bounds.fill(0.0);
        self.operator_evaluations = 0;
        self.roots.reset();
        self.check_state()?;
        let (em0, el0) = self.energy(&self.q, &self.v)?;
        let mut step = LinearNetworkEventStep {
            phase: previous_phase,
            initial_energy_j: finite(em0 + el0)?,
            modal_energy_j: em0,
            connection_energy_j: el0,
            external_work_j: 0.0,
            modal_damping_dissipation_j: 0.0,
            connection_damping_dissipation_j: 0.0,
            friction_dissipation_j: 0.0,
            accumulated_slip_m: 0.0,
            contact_impulse_n_s: 0.0,
            contact_time_moment_n_s2: 0.0,
            energy_residual_j: 0.0,
            energy_tolerance_j: 0.0,
            numerical_work_uncertainty_j: 0.0,
            events: [None; MAX_LINEAR_NETWORK_EVENTS],
            event_count: 0,
            cells: 0,
            operator_evaluations: 0,
            root_evaluations: 0,
        };
        let initial_speed = dot(&self.tangent, &self.v)?;
        let mut phase = match previous_phase {
            LinearNetworkPhase::Fresh => {
                if initial_speed == 0.0 {
                    LinearNetworkPhase::Sticking
                } else {
                    LinearNetworkPhase::Sliding {
                        direction: sign(initial_speed),
                    }
                }
            }
            LinearNetworkPhase::Sticking => {
                if initial_speed.abs() > self.config.sticking_velocity_tolerance_m_s {
                    return Err(Error::PhaseMismatch);
                }
                self.project_velocity(&mut step)?;
                LinearNetworkPhase::Sticking
            }
            LinearNetworkPhase::Sliding { direction } => {
                if direction.abs() != 1 || initial_speed != 0.0 && sign(initial_speed) != direction
                {
                    return Err(Error::PhaseMismatch);
                }
                // A caller may change its held force at a sample boundary
                // that coincides with release. Exact zero actual velocity can
                // support the new static load; a merely small velocity cannot.
                if initial_speed == 0.0 {
                    LinearNetworkPhase::Sticking
                } else {
                    previous_phase
                }
            }
        };
        let mut time = 0.0;
        while time < dt_s {
            poll(gate)?;
            if step.cells >= self.config.maximum_cells {
                return Err(Error::Budget("trajectory cells"));
            }
            if phase == LinearNetworkPhase::Sticking {
                let reaction = self.instant_reaction(external)?;
                if reaction.abs() > static_capacity_n {
                    phase = LinearNetworkPhase::Sliding {
                        direction: sign(reaction),
                    };
                    self.event(
                        &mut step,
                        time,
                        LinearNetworkEventKind::SlipStarted {
                            direction: sign(reaction),
                        },
                    )?;
                }
            }
            let remainder = finite(dt_s - time)?;
            self.prepare_cell(external, phase, kinetic_capacity_n, remainder, gate)?;
            step.cells += 1;
            let h = self.cell_h;
            let tolerance = (self.config.event_time_tolerance_s / h)
                .min(1e-3)
                .max(8.0 * f64::EPSILON);
            let found = match phase {
                LinearNetworkPhase::Sliding { direction } => roots::first_zero(
                    &self.velocity_coeff,
                    direction,
                    tolerance,
                    &mut self.roots,
                    gate,
                )?,
                LinearNetworkPhase::Sticking => roots::first_exit(
                    &self.reaction_coeff,
                    static_capacity_n,
                    tolerance,
                    &mut self.roots,
                    gate,
                )?,
                _ => None,
            };
            let z = found.map_or(1.0, |event| event.fraction);
            self.integrate_cell(
                z,
                time,
                phase,
                kinetic_capacity_n,
                external,
                &mut step,
                gate,
            )?;
            let next_time = if found.is_none() && h == remainder {
                dt_s
            } else {
                finite(time + h * z)?
            };
            if let Some(found) = found {
                // A bounded event-time bracket changes the applied generalized
                // effort over that bracket. Keep its state/work sensitivity
                // even when the trajectory happens to be conservative.
                let timing = self.config.event_time_tolerance_s;
                let impulse_error = finite((static_capacity_n + kinetic_capacity_n) * timing)?;
                let velocity_error =
                    impulse_error * self.tangent.iter().map(|b| b.abs()).fold(0.0, f64::max);
                let position_error = velocity_error * timing;
                let vmax = self.v.iter().map(|v| v.abs()).fold(0.0, f64::max);
                let qmax = self.q.iter().map(|q| q.abs()).fold(0.0, f64::max);
                let rest_force_bound = self
                    .links
                    .iter()
                    .map(|link| {
                        link.stiffness_n_per_m
                            * link.rest_extension_m.abs()
                            * link.column.iter().map(|b| b.abs()).sum::<f64>()
                    })
                    .sum::<f64>();
                let force_scale = external.iter().map(|f| f.abs()).sum::<f64>()
                    + self.stiffness_norm * qmax * n as f64
                    + rest_force_bound;
                step.numerical_work_uncertainty_j = finite(
                    step.numerical_work_uncertainty_j
                        + n as f64
                            * (vmax * velocity_error + 0.5 * velocity_error * velocity_error)
                        + position_error * (force_scale + kinetic_capacity_n * self.norm_b1),
                )?;
                if matches!(phase, LinearNetworkPhase::Sliding { .. }) {
                    self.project_velocity(&mut step)?;
                    self.event(&mut step, next_time, LinearNetworkEventKind::Stopped)?;
                    phase = LinearNetworkPhase::Sticking;
                } else {
                    phase = LinearNetworkPhase::Sliding {
                        direction: found.direction,
                    };
                    self.event(
                        &mut step,
                        next_time,
                        LinearNetworkEventKind::SlipStarted {
                            direction: found.direction,
                        },
                    )?;
                }
                // A zero-duration phase transition is allowed. Its next phase
                // has different dynamics, and the finite event budget bounds it.
            } else if next_time == time {
                return Err(Error::EventUnresolved);
            }
            time = next_time;
        }
        let (em, el) = self.energy(&self.q, &self.v)?;
        if static_capacity_n == 0.0 && kinetic_capacity_n == 0.0 {
            let speed = dot(&self.tangent, &self.v)?;
            phase = if speed == 0.0 {
                LinearNetworkPhase::Sticking
            } else {
                LinearNetworkPhase::Sliding {
                    direction: sign(speed),
                }
            };
        }
        step.phase = phase;
        step.modal_energy_j = em;
        step.connection_energy_j = el;
        step.external_work_j = finite(
            external
                .iter()
                .zip(self.q.iter().zip(q0))
                .map(|(f, (q, q0))| f * (q - q0))
                .sum(),
        )?;
        step.connection_damping_dissipation_j = finite(self.link_losses.iter().sum())?;
        for (i, link) in self.links.iter().enumerate() {
            self.link_forces[i] = finite(
                -link.stiffness_n_per_m * (dot(&link.column, &self.q)? - link.rest_extension_m)
                    - link.damping_n_s_per_m * dot(&link.column, &self.v)?,
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
        step.energy_tolerance_j = finite(
            self.config.energy_absolute_tolerance_j + self.config.energy_relative_tolerance * scale,
        )?;
        step.operator_evaluations = self.operator_evaluations;
        step.root_evaluations = self.roots.used();
        if step.numerical_work_uncertainty_j > self.config.maximum_numerical_work_uncertainty_j
            || step.energy_residual_j.abs() + step.numerical_work_uncertainty_j
                > step.energy_tolerance_j
        {
            return Err(Error::EnergyAccuracy {
                residual_j: step.energy_residual_j,
                estimated_error_j: step.numerical_work_uncertainty_j,
                tolerance_j: step
                    .energy_tolerance_j
                    .min(self.config.maximum_numerical_work_uncertainty_j),
            });
        }
        poll(gate)?;
        Ok(step)
    }

    fn event(
        &self,
        step: &mut LinearNetworkEventStep,
        time: f64,
        kind: LinearNetworkEventKind,
    ) -> Result<(), Error> {
        if step.event_count >= self.config.maximum_events {
            return Err(Error::Budget("physical events"));
        }
        step.events[step.event_count] = Some(LinearNetworkEvent { time_s: time, kind });
        step.event_count += 1;
        Ok(())
    }
    fn check_state(&self) -> Result<(), Error> {
        if self
            .q
            .iter()
            .any(|x| x.abs() > self.config.maximum_abs_modal_coordinate)
            || self
                .v
                .iter()
                .any(|v| v.abs() > self.config.maximum_abs_modal_velocity)
        {
            return Err(Error::Budget("physical modal state"));
        }
        Ok(())
    }
    fn energy(&self, q: &[f64], v: &[f64]) -> Result<(f64, f64), Error> {
        let modal = finite(
            q.iter()
                .zip(v)
                .zip(&self.stiffness)
                .map(|((q, v), k)| 0.5 * (v * v + k * q * q))
                .sum(),
        )?;
        let mut links = 0.0;
        for l in &self.links {
            let z = finite(dot(&l.column, q)? - l.rest_extension_m)?;
            links = finite(links + 0.5 * l.stiffness_n_per_m * z * z)?;
        }
        Ok((modal, links))
    }
    fn instant_reaction(&mut self, external: &[f64]) -> Result<f64, Error> {
        apply(
            &self.stiffness,
            &self.damping,
            &self.links,
            &self.q,
            &self.v,
            Some(external),
            &mut self.accel,
        )?;
        self.operator_evaluations += 1;
        if self.operator_evaluations > self.config.maximum_operator_evaluations {
            return Err(Error::Budget("operator evaluations"));
        }
        finite(dot(&self.tangent, &self.accel)? / self.norm_b2)
    }
    fn project_velocity(&mut self, step: &mut LinearNetworkEventStep) -> Result<(), Error> {
        let drift = dot(&self.tangent, &self.v)?;
        let correction = finite(drift / self.norm_b2)?;
        let old_norm = det::sqrt(dot(&self.v, &self.v)?);
        let correction_norm = correction.abs() * det::sqrt(self.norm_b2);
        step.numerical_work_uncertainty_j = finite(
            step.numerical_work_uncertainty_j
                + old_norm * correction_norm
                + 0.5 * correction_norm * correction_norm,
        )?;
        for (v, b) in self.v.iter_mut().zip(&self.tangent) {
            *v = finite(*v - b * correction)?;
        }
        Ok(())
    }
}

fn validate_config(c: LinearNetworkEventConfig) -> Result<(), Error> {
    if c.maximum_modes == 0
        || c.maximum_connections == usize::MAX
        || c.maximum_workspace_bytes == 0
        || c.maximum_cells == 0
        || c.maximum_events == 0
        || c.maximum_events > MAX_LINEAR_NETWORK_EVENTS
        || c.maximum_operator_evaluations == 0
        || c.maximum_root_evaluations == 0
        || c.maximum_root_iterations == 0
        || [
            c.maximum_step_s,
            c.event_time_tolerance_s,
            c.sticking_velocity_tolerance_m_s,
            c.maximum_abs_modal_coordinate,
            c.maximum_abs_modal_velocity,
            c.energy_absolute_tolerance_j,
            c.maximum_numerical_work_uncertainty_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || !c.energy_relative_tolerance.is_finite()
        || c.energy_relative_tolerance < 0.0
        || c.energy_relative_tolerance >= 1.0
    {
        return Err(Error::Invalid(
            "positive finite event accuracy and work budgets",
        ));
    }
    Ok(())
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Nonfinite)
    }
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, Error> {
    finite(a.iter().zip(b).map(|(a, b)| a * b).sum())
}
fn sign(x: f64) -> i8 {
    if x < 0.0 { -1 } else { 1 }
}
fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn apply(
    k: &[f64],
    c: &[f64],
    links: &[LinearNetworkLink],
    q: &[f64],
    v: &[f64],
    external: Option<&[f64]>,
    out: &mut [f64],
) -> Result<(), Error> {
    for i in 0..out.len() {
        out[i] = finite(external.map_or(0.0, |f| f[i]) - k[i] * q[i] - c[i] * v[i])?;
    }
    for link in links {
        let z = finite(
            dot(&link.column, q)?
                - if external.is_some() {
                    link.rest_extension_m
                } else {
                    0.0
                },
        )?;
        let reaction =
            finite(link.stiffness_n_per_m * z + link.damping_n_s_per_m * dot(&link.column, v)?)?;
        for (out, b) in out.iter_mut().zip(&link.column) {
            *out = finite(*out - b * reaction)?;
        }
    }
    Ok(())
}
