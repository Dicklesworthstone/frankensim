//! Finite normal motion on the original quasistatic JKR contact branches.
//!
//! One fixed normal coordinate has caller-supplied relative mass and held
//! external force. Each attached substep is the scalar condensation of the
//! complete kinetic/contact discrete gradient: its original path mean supplies
//! both impulse and work. Optional Kelvin--Voigt support participates in the
//! same equation with separate spring and viscous work/impulse accounts.
//! Viscous loss is the passive midpoint quadrature c*h*v_mid^2; the quadratic
//! geometry envelope is not an exact continuous dashpot trajectory. Detached
//! motion is midpoint motion, or exact flight without spring or damping.
//! Touch and the displacement fold split the trajectory at their original
//! thresholds. A snap changes only contact phase/storage, never momentum.
//!
//! The contact radius is still quasistatic and its source remains frozen.
//! Released snap energy is explicit outgoing energy, not invented heat or a
//! restitution impulse. Fast elastic waves, finite-range precontact attraction,
//! evolving geometry/material, and other contact coordinates require their own
//! owners. Smooth motion is second-order time discrete; root/energy gates are
//! floating-point checks, not interval certificates or time-error estimates.

use super::{JkrBranch, JkrError, JkrJump, JkrResponse, JkrState, SourceJkrSphereContact};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;

mod advance;
mod damping;
mod elastic;
mod solve;
pub use damping::{JkrNormalDissipationSegment, JkrNormalMotionTrace};

/// Fixed event storage; the caller may impose a smaller limit.
pub const MAX_JKR_NORMAL_EVENTS: usize = 16;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalMotionBudget {
    pub maximum_step_s: f64,
    /// Physical time-resolution control; no automatic accuracy claim is made.
    pub maximum_substep_s: f64,
    /// Bounds all trial substeps, including failed attempts before halving.
    pub maximum_substeps: usize,
    pub maximum_root_iterations: usize,
    pub maximum_events: usize,
    /// Lowest admitted detached approach; must be below the original fold.
    pub minimum_indentation_m: f64,
    /// Actual travel, including a within-branch reversal, over one call.
    pub maximum_accumulated_travel_m: f64,
    pub maximum_normal_speed_m_s: f64,
    /// Bounds held load, contact, spring, viscous and combined support forces.
    pub maximum_force_n: f64,
    pub coordinate_absolute_tolerance_m: f64,
    pub coordinate_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
    pub impulse_absolute_tolerance_n_s: f64,
    pub impulse_relative_tolerance: f64,
}

/// Source, mass and numerical policy are retained by this immutable token.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalMotionState {
    identity: ContentHash,
    contact: JkrState,
    velocity_m_s: f64,
    time_s: f64,
    steps: u64,
}
impl JkrNormalMotionState {
    pub const fn contact_state(self) -> JkrState {
        self.contact
    }
    pub const fn velocity_m_s(self) -> f64 {
        self.velocity_m_s
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    pub const fn steps(self) -> u64 {
        self.steps
    }
}

/// An original source transition at fixed approach and continuous velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalMotionEvent {
    pub offset_s: f64,
    pub jump: JkrJump,
    pub velocity_m_s: f64,
    before: JkrState,
    after: JkrState,
}
impl JkrNormalMotionEvent {
    pub const fn before(self) -> JkrState {
        self.before
    }
    pub const fn after(self) -> JkrState {
        self.after
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalMotionStep {
    initial: JkrNormalMotionState,
    state: JkrNormalMotionState,
    pub contact: JkrResponse,
    pub normal_increment_m: f64,
    pub accumulated_travel_m: f64,
    /// Signed contact reaction; its force on the relative mass is negative.
    pub normal_impulse_n_s: f64,
    /// Integral t F dt, with time measured from this call's start.
    pub normal_time_moment_n_s2: f64,
    /// Original smooth contact work, integral F d(delta), excluding snaps.
    pub normal_work_j: f64,
    /// Linear-support reaction; the original contact observations stay separate.
    pub support_impulse_n_s: f64,
    pub support_time_moment_n_s2: f64,
    pub support_work_j: f64,
    pub support_energy_j: f64,
    /// Maximum absolute support force over all accepted interpolation extrema.
    pub support_force_bound_n: f64,
    /// Endpoint constitutive viscous reaction c*v, distinct from its path mean.
    pub damping_force_n: f64,
    pub damping_impulse_n_s: f64,
    pub damping_time_moment_n_s2: f64,
    /// Nonnegative discrete-gradient viscous work, excluding source snaps.
    pub damping_dissipation_j: f64,
    pub damping_force_bound_n: f64,
    /// Maximum |spring + c*v| over the accepted quadratic geometry envelope.
    pub support_reaction_force_bound_n: f64,
    pub applied_work_j: f64,
    /// Energy leaving quasistatic storage at source-owned branch changes.
    pub released_energy_j: f64,
    pub kinetic_energy_j: f64,
    /// Endpoint kinetic, contact and support storage; released energy is separate.
    pub total_energy_j: f64,
    pub energy_residual_j: f64,
    pub impulse_residual_n_s: f64,
    pub coordinate_residual_m: f64,
    events: [Option<JkrNormalMotionEvent>; MAX_JKR_NORMAL_EVENTS],
    event_count: usize,
    pub substeps: usize,
    pub trials: usize,
}
impl JkrNormalMotionStep {
    pub fn reaction_impulse_n_s(&self) -> f64 {
        self.normal_impulse_n_s + self.support_impulse_n_s + self.damping_impulse_n_s
    }
    pub fn reaction_time_moment_n_s2(&self) -> f64 {
        self.normal_time_moment_n_s2 + self.support_time_moment_n_s2 + self.damping_time_moment_n_s2
    }
    pub const fn initial_state(&self) -> JkrNormalMotionState {
        self.initial
    }
    pub const fn state(&self) -> JkrNormalMotionState {
        self.state
    }
    pub const fn event_count(&self) -> usize {
        self.event_count
    }
    /// Original ordered transitions, retained privately with their endpoints.
    pub fn events(&self) -> impl Iterator<Item = &JkrNormalMotionEvent> {
        self.events[..self.event_count].iter().flatten()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrNormalMotionError {
    Normal(JkrError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    TrialBudget,
    EventBudget,
    Range(&'static str),
    Unresolved(&'static str),
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    ImpulseBalance {
        residual_n_s: f64,
        tolerance_n_s: f64,
    },
    CoordinateBalance {
        residual_m: f64,
        tolerance_m: f64,
    },
}
impl From<JkrError> for JkrNormalMotionError {
    fn from(value: JkrError) -> Self {
        Self::Normal(value)
    }
}
impl core::fmt::Display for JkrNormalMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sourced JKR normal motion refused: {self:?}")
    }
}
impl std::error::Error for JkrNormalMotionError {}

/// A source-owned normal trajectory, independent of the receiving body basis.
/// Two unconstrained translating bodies use their reduced mass. A fixed work
/// column in mass-normalized native coordinates uses `1/(B dot B)`.
pub struct JkrNormalMotion<'a> {
    normal: &'a SourceJkrSphereContact,
    mass_kg: f64,
    budget: JkrNormalMotionBudget,
    identity: ContentHash,
    support: elastic::Support,
    damping_n_s_per_m: f64,
}
impl<'a> JkrNormalMotion<'a> {
    pub fn new(
        normal: &'a SourceJkrSphereContact,
        mass_kg: f64,
        budget: JkrNormalMotionBudget,
    ) -> Result<Self, JkrNormalMotionError> {
        Self::new_with_support(normal, mass_kg, budget, 0.0, 0.0)
    }

    /// Add one undamped spring with reaction `k * (indentation - anchor)`.
    /// This is source-independent mechanical support, not contact stiffness or
    /// modified adhesion. The anchor is frozen with the support identity.
    pub fn new_with_support(
        normal: &'a SourceJkrSphereContact,
        mass_kg: f64,
        budget: JkrNormalMotionBudget,
        stiffness_n_per_m: f64,
        zero_force_indentation_m: f64,
    ) -> Result<Self, JkrNormalMotionError> {
        Self::new_with_damped_support(
            normal,
            mass_kg,
            budget,
            stiffness_n_per_m,
            zero_force_indentation_m,
            0.0,
        )
    }

    /// Add a frozen Kelvin--Voigt support. Viscous work is c*h*v_mid^2;
    /// it is outgoing mechanical energy with no implicit heat partition.
    pub fn new_with_damped_support(
        normal: &'a SourceJkrSphereContact,
        mass_kg: f64,
        budget: JkrNormalMotionBudget,
        stiffness_n_per_m: f64,
        zero_force_indentation_m: f64,
        damping_n_s_per_m: f64,
    ) -> Result<Self, JkrNormalMotionError> {
        let positive = [
            mass_kg,
            budget.maximum_step_s,
            budget.maximum_substep_s,
            budget.maximum_accumulated_travel_m,
            budget.maximum_normal_speed_m_s,
            budget.maximum_force_n,
            budget.coordinate_absolute_tolerance_m,
            budget.energy_absolute_tolerance_j,
            budget.impulse_absolute_tolerance_n_s,
        ];
        let relative = [
            budget.coordinate_relative_tolerance,
            budget.energy_relative_tolerance,
            budget.impulse_relative_tolerance,
        ];
        if positive.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || relative.iter().any(|x| !x.is_finite() || *x < 0.0)
            || !budget.minimum_indentation_m.is_finite()
            || budget.minimum_indentation_m >= normal.minimum_attached_indentation_m()
            || budget.maximum_substep_s > budget.maximum_step_s
            || !(1..=4096).contains(&budget.maximum_substeps)
            || !(1..=64).contains(&budget.maximum_root_iterations)
            || !(1..=MAX_JKR_NORMAL_EVENTS).contains(&budget.maximum_events)
            || !damping_n_s_per_m.is_finite()
            || damping_n_s_per_m < 0.0
        {
            return Err(JkrNormalMotionError::Invalid(
                "finite mass and bounded explicit motion controls",
            ));
        }
        let support = elastic::Support::new(
            normal,
            stiffness_n_per_m,
            zero_force_indentation_m,
            budget.maximum_root_iterations,
        )?;
        // Preserve the existing zero-damping identity and arithmetic path.
        let mut hash = DomainHasher::new(if damping_n_s_per_m == 0.0 {
            "fs-contact/jkr-normal-motion/v2"
        } else {
            "fs-contact/jkr-normal-motion/v3"
        });
        hash.update(normal.identity().as_bytes());
        for value in positive.into_iter().chain(relative).chain([
            budget.minimum_indentation_m,
            stiffness_n_per_m,
            zero_force_indentation_m,
        ]) {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for count in [
            budget.maximum_substeps,
            budget.maximum_root_iterations,
            budget.maximum_events,
        ] {
            hash.update(&(count as u64).to_le_bytes());
        }
        if damping_n_s_per_m > 0.0 {
            hash.update(&damping_n_s_per_m.to_bits().to_le_bytes());
        }
        Ok(Self {
            normal,
            mass_kg,
            budget,
            identity: hash.finalize(),
            support,
            damping_n_s_per_m,
        })
    }

    pub const fn normal(&self) -> &'a SourceJkrSphereContact {
        self.normal
    }
    pub const fn effective_mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn budget(&self) -> JkrNormalMotionBudget {
        self.budget
    }
    pub const fn support_stiffness_n_per_m(&self) -> f64 {
        self.support.stiffness
    }
    pub const fn support_zero_force_indentation_m(&self) -> f64 {
        self.support.anchor
    }
    pub const fn support_damping_n_s_per_m(&self) -> f64 {
        self.damping_n_s_per_m
    }
    pub fn support_energy_j(
        &self,
        state: JkrNormalMotionState,
    ) -> Result<f64, JkrNormalMotionError> {
        self.response(state)?;
        self.support.energy(state.contact.indentation_m())
    }

    /// The caller explicitly supplies a detached or stable attached source
    /// token. A marginal attached initial fold is refused; initialize with its
    /// source transition's detached successor instead. Detached touch at zero
    /// is retained and attaches at offset zero on the first accepted advance.
    pub fn initialize(
        &self,
        contact: JkrState,
        velocity_m_s: f64,
    ) -> Result<JkrNormalMotionState, JkrNormalMotionError> {
        if contact.model_identity != self.normal.identity() {
            return Err(JkrNormalMotionError::StateMismatch);
        }
        if contact.branch() == JkrBranch::Attached
            && contact.indentation_m() <= self.normal.minimum_attached_indentation_m()
        {
            return Err(JkrError::SingularTangent.into());
        }
        let state = JkrNormalMotionState {
            identity: self.identity,
            contact,
            velocity_m_s,
            time_s: 0.0,
            steps: 0,
        };
        let response = self.response(state)?;
        let kinetic = finite(
            0.5 * self.mass_kg * velocity_m_s * velocity_m_s,
            "initial normal kinetic energy",
        )?;
        finite(
            kinetic
                + response.reversible_energy_j
                + self.support.energy(contact.indentation_m())?,
            "initial normal total energy",
        )?;
        Ok(state)
    }

    /// Evaluate the retained phase without advancing history or causing a snap.
    pub fn response(
        &self,
        state: JkrNormalMotionState,
    ) -> Result<JkrResponse, JkrNormalMotionError> {
        if state.identity != self.identity || state.contact.model_identity != self.normal.identity()
        {
            return Err(JkrNormalMotionError::StateMismatch);
        }
        if !state.velocity_m_s.is_finite()
            || state.velocity_m_s.abs() > self.budget.maximum_normal_speed_m_s
            || !state.time_s.is_finite()
            || state.time_s < 0.0
            || state.contact.indentation_m() < self.budget.minimum_indentation_m
        {
            return Err(JkrNormalMotionError::Range(
                "accepted speed, time or detached approach",
            ));
        }
        let response = match state.contact.branch() {
            JkrBranch::Attached => self
                .normal
                .evaluate_attached(state.contact.indentation_m())?,
            JkrBranch::Detached => {
                self.normal.detached_state(state.contact.indentation_m())?;
                self.normal.detached_response(state.contact.indentation_m())
            }
        };
        if response.normal_force_n.abs() > self.budget.maximum_force_n {
            return Err(JkrNormalMotionError::Range("normal force"));
        }
        if self.support.force(state.contact.indentation_m())?.abs() > self.budget.maximum_force_n {
            return Err(JkrNormalMotionError::Range("support force"));
        }
        let viscous = damping::force(self, state.velocity_m_s)?;
        let combined = finite(
            self.support.force(state.contact.indentation_m())? + viscous,
            "endpoint support reaction",
        )?;
        if viscous.abs().max(combined.abs()) > self.budget.maximum_force_n {
            return Err(JkrNormalMotionError::Range("viscous support reaction"));
        }
        Ok(response)
    }

    /// Advance held external load over the complete requested interval.
    /// All trials and branch transitions remain local until the returned token
    /// is accepted; cancellation or a budget/domain refusal publishes nothing.
    pub fn advance(
        &self,
        previous: JkrNormalMotionState,
        external_force_n: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<JkrNormalMotionStep, JkrNormalMotionError> {
        advance::run(self, previous, external_force_n, dt_s, gate, None)
    }

    fn coordinate_tolerance(&self, scale: f64) -> f64 {
        self.budget.coordinate_absolute_tolerance_m
            + self.budget.coordinate_relative_tolerance
                * scale.max(self.normal.indentation_scale_m())
    }
}

fn poll(gate: &CancelGate) -> Result<(), JkrNormalMotionError> {
    if gate.is_requested() {
        Err(JkrNormalMotionError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64, field: &'static str) -> Result<f64, JkrNormalMotionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(JkrNormalMotionError::Range(field))
    }
}
