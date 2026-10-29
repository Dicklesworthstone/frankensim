//! Mixed-mode pull-off with finite normal AND tangential loading stiffness.
//!
//! A fixed vector shear stage D_t and a moving normal stage D_n load the
//! ORIGINAL conservative full-stick contact through positive springs K_t,K_n.
//! At radius a, u=K_t D_t/(K_t+k a), k=8G*, and the tensile Griffith branch is
//! delta=a^2/R-sqrt((2 pi W a-k|u|^2/2)/E*). The actual normal stage is
//! D_n=delta+N/K_n. Both spring energies belong to the total potential.
//! This is the lambda=1 specialization of Papangelo, Cricri and Ciavarella,
//! JMPS 144 (2020), 104099, equations 18--26,
//! <https://doi.org/10.1016/j.jmps.2020.104099>, using the original two-body
//! Mindlin modulus instead of assuming an incompressible rigid/soft pair.
//!
//! Normal withdrawal reaches the FIRST sourced gross-slip limit or apparatus
//! radius fold. A fold can produce a detached, relaxed successor only if its
//! normal stage is nonpositive and its stored-energy release is nonnegative.
//! That release is retained as an explicit cumulative transfer to an unresolved
//! fast transient. It is not silently assigned to heat or erased. Other folds
//! and gross slip return the actual boundary and untraversed stage travel.
//! Detached states never spontaneously rebond; positive detached approach
//! refuses. This is a quasistatic failure handoff, not an inertial snap history,
//! post-slip JKR law, partial-slip model or empirical mixed-mode toughness fit.

use super::loading::JkrShearBoundaryKind;
use super::{JkrMixedModeControl, JkrMixedModeError, JkrMixedModeResponse, SourceJkrMixedMode};
use crate::interface_binding::adhesive_shear::{AdhesiveShearError, SourceAdhesiveShear};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;

/// The maximum-radius attached state is the preloaded initial condition.
/// Its preparation work is not inferred. The admitted range must include the
/// first boundary; each construction/trial has one complete evaluation budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrPullOffLimits {
    pub maximum_contact_radius_m: f64,
    pub maximum_absolute_normal_stage_m: f64,
    pub maximum_iterations: usize,
    pub maximum_evaluations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrPullOffBranch {
    Attached,
    Detached,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrPullOffState {
    identity: ContentHash,
    branch: JkrPullOffBranch,
    normal_stage_m: f64,
    released_energy_j: f64,
}
impl JkrPullOffState {
    pub const fn branch(self) -> JkrPullOffBranch {
        self.branch
    }
    pub const fn normal_stage_m(self) -> f64 {
        self.normal_stage_m
    }
    /// Energy already transferred out of the quasistatic mechanical owner.
    pub const fn accumulated_released_energy_j(self) -> f64 {
        self.released_energy_j
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrPullOffResponse {
    state: JkrPullOffState,
    pub contact: Option<JkrMixedModeResponse>,
    /// Actual contact-tip approach and both elastic shear shifts.
    pub tip_coordinates_m: [f64; 3],
    pub stage_force_n: [f64; 3],
    pub normal_spring_energy_j: f64,
    pub tangential_spring_energy_j: f64,
    pub total_stored_energy_j: f64,
    pub static_capacity_n: f64,
    pub static_margin_n: f64,
    /// Radius Hessian after relaxing all three tip coordinates against springs.
    pub apparatus_radius_curvature_n_per_m: f64,
    /// Complete symmetric derivative of stage reactions with respect to stages.
    /// Absent at the radius fold. The zero matrix applies after detachment.
    pub stage_stiffness_n_per_m: Option<[[f64; 3]; 3]>,
    pub spring_force_residual_n: [f64; 3],
}
impl JkrPullOffResponse {
    pub const fn state(&self) -> JkrPullOffState {
        self.state
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrPullOffBoundary {
    pub kind: JkrShearBoundaryKind,
    pub response: JkrPullOffResponse,
    /// Numerical radius bracket; this is not an interval-arithmetic enclosure.
    pub radius_bracket_m: [f64; 2],
    pub relaxed_separation_admissible: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrPullOffEventKind {
    GrossSlip,
    RadiusInstability,
    Coincident,
    Separation,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrPullOffEvent {
    pub kind: JkrPullOffEventKind,
    pub normal_stage_m: f64,
    pub before_tip_coordinates_m: [f64; 3],
    pub after_tip_coordinates_m: [f64; 3],
    pub released_energy_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrPullOffTransition {
    pub response: JkrPullOffResponse,
    pub event: Option<JkrPullOffEvent>,
    pub requested_normal_stage_m: f64,
    pub untraversed_normal_stage_m: f64,
    /// Integral N dD_n. The held tangential stages perform zero external work.
    pub prescribed_work_j: f64,
    pub released_energy_j: f64,
    pub evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrPullOffError {
    Contact(JkrMixedModeError),
    Shear(AdhesiveShearError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    Budget,
    NoConvergence,
    OutsideStageEnvelope,
    NoAdmittedAttachedBranch,
    BeyondFirstBoundary,
    DetachedOverlap,
}
impl core::fmt::Display for JkrPullOffError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "mixed-mode finite-apparatus pull-off refused: {self:?}")
    }
}
impl std::error::Error for JkrPullOffError {}
impl From<JkrMixedModeError> for JkrPullOffError {
    fn from(error: JkrMixedModeError) -> Self {
        Self::Contact(error)
    }
}
impl From<AdhesiveShearError> for JkrPullOffError {
    fn from(error: AdhesiveShearError) -> Self {
        Self::Shear(error)
    }
}

/// Source-bound, immutable, allocation-free runtime continuation.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrPullOffFixture {
    contact: SourceJkrMixedMode,
    shear: SourceAdhesiveShear,
    normal_stiffness: f64,
    tangential_stiffness: f64,
    tangential_stage: [f64; 2],
    limits: JkrPullOffLimits,
    identity: ContentHash,
    initial: JkrPullOffResponse,
    boundary: JkrPullOffBoundary,
}

impl JkrPullOffFixture {
    pub fn new(
        contact: SourceJkrMixedMode,
        shear: SourceAdhesiveShear,
        normal_stiffness_n_per_m: f64,
        tangential_stiffness_n_per_m: f64,
        tangential_stage_m: [f64; 2],
        limits: JkrPullOffLimits,
        gate: &CancelGate,
    ) -> Result<Self, JkrPullOffError> {
        for value in [
            normal_stiffness_n_per_m,
            tangential_stiffness_n_per_m,
            limits.maximum_contact_radius_m,
            limits.maximum_absolute_normal_stage_m,
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(JkrPullOffError::Invalid("positive finite stiffness/range"));
            }
        }
        if tangential_stage_m.iter().any(|v| !v.is_finite())
            || limits.maximum_iterations == 0
            || limits.maximum_iterations > 256
            || limits.maximum_evaluations == 0
        {
            return Err(JkrPullOffError::Invalid("stage or solve budget"));
        }
        if shear.normal_state_identity() != contact.normal_source().normal_state().identity() {
            return Err(JkrPullOffError::StateMismatch);
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-pull-off.v1");
        hash.update(contact.identity().as_bytes());
        hash.update(shear.resolved().card_identity().as_bytes());
        for value in [
            normal_stiffness_n_per_m,
            tangential_stiffness_n_per_m,
            tangential_stage_m[0],
            tangential_stage_m[1],
            limits.maximum_contact_radius_m,
            limits.maximum_absolute_normal_stage_m,
            shear.law().static_shear_strength_pa(),
            shear.law().kinetic_shear_strength_pa(),
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&(limits.maximum_iterations as u64).to_le_bytes());
        hash.update(&(limits.maximum_evaluations as u64).to_le_bytes());
        let identity = hash.finalize();
        let parameters = Parameters {
            contact: &contact,
            shear: &shear,
            kn: normal_stiffness_n_per_m,
            kt: tangential_stiffness_n_per_m,
            stage: tangential_stage_m,
            identity,
        };
        let mut budget = Work::new(limits.maximum_evaluations, gate);
        let maximum = limits.maximum_contact_radius_m;
        let at_maximum = parameters.kinematics(maximum, &mut budget)?;
        if at_maximum.ratio <= 1.0 {
            return Err(JkrPullOffError::NoAdmittedAttachedBranch);
        }
        let initial = parameters.response(maximum, at_maximum, false, &mut budget)?;
        if initial.static_margin_n < 0.0 {
            return Err(JkrPullOffError::NoAdmittedAttachedBranch);
        }
        // The fold ratio is strictly increasing on the tensile radius branch.
        // Its extension is zero below the unique positive-opening radius.
        let mut lo = 0.0;
        let mut hi = maximum;
        let mut complete = false;
        for _ in 0..limits.maximum_iterations {
            let a = f64::midpoint(lo, hi);
            if a == lo || a == hi || hi - lo <= 16.0 * f64::EPSILON * hi {
                complete = true;
                break;
            }
            let ratio = parameters.fold_ratio(a, &mut budget)?;
            if ratio < 1.0 {
                lo = a;
            } else {
                hi = a;
            }
        }
        if !complete {
            return Err(JkrPullOffError::NoConvergence);
        }
        let fold_bracket = [lo, hi];
        let shear_magnitude = tangential_stage_m[0].hypot(tangential_stage_m[1]);
        let strength = shear.law().static_shear_strength_pa();
        // T/(tau*pi*a^2) decreases strictly with a. Its crossing is a positive
        // quadratic root, so it competes with the fold without scanning paths.
        let slip_radius = if shear_magnitude == 0.0 {
            0.0
        } else {
            if strength == 0.0 {
                return Err(JkrPullOffError::NoAdmittedAttachedBranch);
            }
            let b = tangential_stiffness_n_per_m / contact.tangential_stiffness_per_radius;
            let c = (tangential_stiffness_n_per_m * shear_magnitude) / (PI * strength);
            finite(b)?;
            finite(c)?;
            let root = c / (0.5 * b + 0.5 * b.hypot(2.0 * c.sqrt()));
            finite(root)?;
            root
        };
        let slip_bracket = [
            slip_radius * (1.0 - 8.0 * f64::EPSILON),
            slip_radius * (1.0 + 8.0 * f64::EPSILON),
        ];
        let (kind, bracket) = if slip_bracket[0] > hi {
            (JkrShearBoundaryKind::GrossSlip, slip_bracket)
        } else if slip_bracket[1] < lo {
            (JkrShearBoundaryKind::RadiusInstability, fold_bracket)
        } else {
            (
                JkrShearBoundaryKind::Coincident,
                [lo.min(slip_bracket[0]), hi.max(slip_bracket[1])],
            )
        };
        let boundary_radius = bracket[1];
        if boundary_radius > maximum {
            return Err(JkrPullOffError::NoAdmittedAttachedBranch);
        }
        let kin = parameters.kinematics(boundary_radius, &mut budget)?;
        let boundary_response = parameters.response(
            boundary_radius,
            kin,
            kind != JkrShearBoundaryKind::GrossSlip,
            &mut budget,
        )?;
        if boundary_response.static_margin_n < 0.0 {
            return Err(JkrPullOffError::NoConvergence);
        }
        for stage in [
            initial.state.normal_stage_m,
            boundary_response.state.normal_stage_m,
        ] {
            if stage.abs() > limits.maximum_absolute_normal_stage_m {
                return Err(JkrPullOffError::OutsideStageEnvelope);
            }
        }
        let boundary = JkrPullOffBoundary {
            kind,
            response: boundary_response,
            radius_bracket_m: bracket,
            relaxed_separation_admissible: kind == JkrShearBoundaryKind::RadiusInstability
                && boundary_response.state.normal_stage_m <= 0.0
                && boundary_response.total_stored_energy_j >= 0.0,
        };
        budget.poll()?;
        Ok(Self {
            contact,
            shear,
            normal_stiffness: normal_stiffness_n_per_m,
            tangential_stiffness: tangential_stiffness_n_per_m,
            tangential_stage: tangential_stage_m,
            limits,
            identity,
            initial,
            boundary,
        })
    }

    pub const fn contact(&self) -> &SourceJkrMixedMode {
        &self.contact
    }
    pub const fn shear_source(&self) -> &SourceAdhesiveShear {
        &self.shear
    }
    pub const fn normal_stiffness_n_per_m(&self) -> f64 {
        self.normal_stiffness
    }
    pub const fn tangential_stiffness_n_per_m(&self) -> f64 {
        self.tangential_stiffness
    }
    pub const fn tangential_stage_m(&self) -> [f64; 2] {
        self.tangential_stage
    }
    pub const fn initial_response(&self) -> JkrPullOffResponse {
        self.initial
    }
    pub const fn first_boundary(&self) -> JkrPullOffBoundary {
        self.boundary
    }

    pub fn evaluate_attached(
        &self,
        normal_stage_m: f64,
        gate: &CancelGate,
    ) -> Result<JkrPullOffResponse, JkrPullOffError> {
        let mut budget = Work::new(self.limits.maximum_evaluations, gate);
        let response = self.attached(normal_stage_m, &mut budget)?;
        budget.poll()?;
        Ok(response)
    }

    /// Follow an elastic segment or its first failure event. On separation,
    /// unused withdrawal continues with zero force on the detached branch.
    /// Gross-slip/unsupported-fold handoffs retain all actual stored energy.
    pub fn advance(
        &self,
        state: JkrPullOffState,
        requested_normal_stage_m: f64,
        gate: &CancelGate,
    ) -> Result<JkrPullOffTransition, JkrPullOffError> {
        if state.identity != self.identity {
            return Err(JkrPullOffError::StateMismatch);
        }
        self.check_stage(requested_normal_stage_m)?;
        let mut budget = Work::new(self.limits.maximum_evaluations, gate);
        budget.poll()?;
        let mut event = None;
        let mut release = 0.0;
        let before_energy;
        let mut response;
        if state.branch == JkrPullOffBranch::Detached {
            before_energy = 0.0;
            response = self.detached(requested_normal_stage_m, state.released_energy_j)?;
        } else {
            let before = self.attached(state.normal_stage_m, &mut budget)?;
            before_energy = before.total_stored_energy_j;
            let boundary = self.boundary.response;
            if requested_normal_stage_m < boundary.state.normal_stage_m {
                response = boundary;
                let mut kind = match self.boundary.kind {
                    JkrShearBoundaryKind::GrossSlip => JkrPullOffEventKind::GrossSlip,
                    JkrShearBoundaryKind::RadiusInstability => {
                        JkrPullOffEventKind::RadiusInstability
                    }
                    JkrShearBoundaryKind::Coincident => JkrPullOffEventKind::Coincident,
                };
                let mut after_tip = boundary.tip_coordinates_m;
                if self.boundary.relaxed_separation_admissible {
                    release = boundary.total_stored_energy_j;
                    let accumulated = state.released_energy_j + release;
                    finite(accumulated)?;
                    response = self.detached(requested_normal_stage_m, accumulated)?;
                    kind = JkrPullOffEventKind::Separation;
                    after_tip = [
                        boundary.state.normal_stage_m,
                        self.tangential_stage[0],
                        self.tangential_stage[1],
                    ];
                }
                event = Some(JkrPullOffEvent {
                    kind,
                    normal_stage_m: boundary.state.normal_stage_m,
                    before_tip_coordinates_m: boundary.tip_coordinates_m,
                    after_tip_coordinates_m: after_tip,
                    released_energy_j: release,
                });
            } else {
                response = self.attached(requested_normal_stage_m, &mut budget)?;
            }
            if response.state.branch == JkrPullOffBranch::Attached {
                response.state.released_energy_j = state.released_energy_j;
            }
        }
        let prescribed_work = response.total_stored_energy_j - before_energy + release;
        finite(prescribed_work)?;
        budget.poll()?;
        Ok(JkrPullOffTransition {
            response,
            event,
            requested_normal_stage_m,
            untraversed_normal_stage_m: requested_normal_stage_m - response.state.normal_stage_m,
            prescribed_work_j: prescribed_work,
            released_energy_j: release,
            evaluations: budget.evaluations,
        })
    }

    fn attached(
        &self,
        stage: f64,
        budget: &mut Work<'_>,
    ) -> Result<JkrPullOffResponse, JkrPullOffError> {
        self.check_stage(stage)?;
        if stage < self.boundary.response.state.normal_stage_m {
            return Err(JkrPullOffError::BeyondFirstBoundary);
        }
        if stage > self.initial.state.normal_stage_m {
            return Err(JkrPullOffError::OutsideStageEnvelope);
        }
        if stage == self.initial.state.normal_stage_m {
            budget.poll()?;
            return Ok(self.initial);
        }
        if stage == self.boundary.response.state.normal_stage_m {
            budget.poll()?;
            return Ok(self.boundary.response);
        }
        let parameters = self.parameters();
        let mut lo = self.boundary.radius_bracket_m[1];
        let mut hi = self.limits.maximum_contact_radius_m;
        for _ in 0..self.limits.maximum_iterations {
            let a = f64::midpoint(lo, hi);
            let kin = parameters.kinematics(a, budget)?;
            let residual = kin.stage - stage;
            let scale = stage.abs().max(kin.stage.abs()).max(kin.opening);
            if residual.abs() <= 32.0 * f64::EPSILON * scale {
                let mut response = parameters.response(a, kin, false, budget)?;
                // Store the actual requested stage and expose its force-balance
                // residual instead of accumulating the radius root's stage error.
                response.state.normal_stage_m = stage;
                let extension = stage - response.tip_coordinates_m[0];
                let spring_energy = 0.5 * self.normal_stiffness * extension * extension;
                response.total_stored_energy_j += spring_energy - response.normal_spring_energy_j;
                response.normal_spring_energy_j = spring_energy;
                response.spring_force_residual_n[0] =
                    response.stage_force_n[0] - self.normal_stiffness * extension;
                finite(response.total_stored_energy_j)?;
                return Ok(response);
            }
            if residual > 0.0 {
                hi = a;
            } else {
                lo = a;
            }
            if a == f64::midpoint(lo, hi) {
                break;
            }
        }
        Err(JkrPullOffError::NoConvergence)
    }

    fn detached(&self, stage: f64, released: f64) -> Result<JkrPullOffResponse, JkrPullOffError> {
        if stage > 0.0 {
            return Err(JkrPullOffError::DetachedOverlap);
        }
        Ok(JkrPullOffResponse {
            state: JkrPullOffState {
                identity: self.identity,
                branch: JkrPullOffBranch::Detached,
                normal_stage_m: stage,
                released_energy_j: released,
            },
            contact: None,
            tip_coordinates_m: [stage, self.tangential_stage[0], self.tangential_stage[1]],
            stage_force_n: [0.0; 3],
            normal_spring_energy_j: 0.0,
            tangential_spring_energy_j: 0.0,
            total_stored_energy_j: 0.0,
            static_capacity_n: 0.0,
            static_margin_n: 0.0,
            apparatus_radius_curvature_n_per_m: 0.0,
            stage_stiffness_n_per_m: Some([[0.0; 3]; 3]),
            spring_force_residual_n: [0.0; 3],
        })
    }
    fn check_stage(&self, stage: f64) -> Result<(), JkrPullOffError> {
        if !stage.is_finite() || stage.abs() > self.limits.maximum_absolute_normal_stage_m {
            Err(JkrPullOffError::OutsideStageEnvelope)
        } else {
            Ok(())
        }
    }
    fn parameters(&self) -> Parameters<'_> {
        Parameters {
            contact: &self.contact,
            shear: &self.shear,
            kn: self.normal_stiffness,
            kt: self.tangential_stiffness,
            stage: self.tangential_stage,
            identity: self.identity,
        }
    }
}

struct Work<'a> {
    gate: &'a CancelGate,
    maximum: usize,
    evaluations: usize,
}
impl<'a> Work<'a> {
    fn new(maximum: usize, gate: &'a CancelGate) -> Self {
        Self {
            gate,
            maximum,
            evaluations: 0,
        }
    }
    fn poll(&self) -> Result<(), JkrPullOffError> {
        if self.gate.is_requested() {
            Err(JkrPullOffError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn probe(&mut self) -> Result<(), JkrPullOffError> {
        self.poll()?;
        if self.evaluations >= self.maximum {
            return Err(JkrPullOffError::Budget);
        }
        self.evaluations += 1;
        Ok(())
    }
}

struct Parameters<'a> {
    contact: &'a SourceJkrMixedMode,
    shear: &'a SourceAdhesiveShear,
    kn: f64,
    kt: f64,
    stage: [f64; 2],
    identity: ContentHash,
}
#[derive(Clone, Copy)]
struct Kinematics {
    shift: [f64; 2],
    opening: f64,
    delta: f64,
    stage: f64,
    ratio: f64,
    derivative: f64,
}
impl Parameters<'_> {
    fn fold_ratio(&self, a: f64, work: &mut Work<'_>) -> Result<f64, JkrPullOffError> {
        work.probe()?;
        let k = self.contact.tangential_stiffness_per_radius;
        let w = self
            .contact
            .normal_source()
            .normal_state()
            .adhesion_energy_j_per_m2();
        let u = self.stage[0].hypot(self.stage[1]) * (self.kt / (self.kt + k * a));
        if a == 0.0 || 2.0 * PI * w * a <= 0.5 * k * u * u {
            return Ok(0.0);
        }
        Ok(self.kinematics_unmetered(a)?.ratio)
    }
    fn kinematics(&self, a: f64, work: &mut Work<'_>) -> Result<Kinematics, JkrPullOffError> {
        work.probe()?;
        self.kinematics_unmetered(a)
    }
    fn kinematics_unmetered(&self, a: f64) -> Result<Kinematics, JkrPullOffError> {
        let normal = self.contact.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let w = normal.normal_state().adhesion_energy_j_per_m2();
        let r = normal.effective_radius_m();
        let k = self.contact.tangential_stiffness_per_radius;
        let shear_denominator = self.kt + k * a;
        finite(shear_denominator)?;
        let shift = self.stage.map(|s| s * (self.kt / shear_denominator));
        let u = shift[0].hypot(shift[1]);
        let radicand = (2.0 * PI * w * a - 0.5 * k * u * u) / e;
        if !(radicand.is_finite() && radicand > 0.0) {
            return Err(JkrPullOffError::NoAdmittedAttachedBranch);
        }
        let opening = radicand.sqrt();
        let delta = a * (a / r) - opening;
        let normal_force = 2.0 * e * a * (delta - a * (a / r) / 3.0);
        let stage = delta + normal_force / self.kn;
        let alpha = 2.0 * e / self.kn;
        let q = 1.0 + alpha * a;
        let radicand_derivative = (2.0 * PI * w + k * k * u * u / (self.kt + k * a)) / e;
        let denominator = alpha * opening + q * radicand_derivative / (2.0 * opening);
        let ratio = (2.0 * a * q / r) / denominator;
        let derivative = denominator * (ratio - 1.0);
        for value in [delta, stage, ratio, derivative] {
            finite(value)?;
        }
        Ok(Kinematics {
            shift,
            opening,
            delta,
            stage,
            ratio,
            derivative,
        })
    }

    fn response(
        &self,
        a: f64,
        kin: Kinematics,
        marginal: bool,
        work: &mut Work<'_>,
    ) -> Result<JkrPullOffResponse, JkrPullOffError> {
        work.probe()?;
        if kin.ratio < 1.0 {
            return Err(JkrPullOffError::NoAdmittedAttachedBranch);
        }
        let normal = self.contact.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let w = normal.normal_state().adhesion_energy_j_per_m2();
        let r = normal.effective_radius_m();
        let k = self.contact.tangential_stiffness_per_radius;
        let h = 4.0 * e * a * kin.opening / r - 2.0 * PI * w;
        let contact = self.contact.response(
            a,
            kin.delta,
            kin.shift,
            JkrMixedModeControl::Displacement,
            h,
            h - 4.0 * PI * w,
        )?;
        let capacity = self
            .shear
            .law()
            .static_capacity_n(
                contact
                    .friction_contact_state()
                    .map_err(AdhesiveShearError::from)?,
            )
            .map_err(AdhesiveShearError::from)?;
        let normal_energy = 0.5 * contact.normal_force_n * (contact.normal_force_n / self.kn);
        let shear_force = contact.tangential_force_n[0].hypot(contact.tangential_force_n[1]);
        let shear_energy = 0.5 * shear_force * (shear_force / self.kt);
        let total = contact.reversible_energy_j + normal_energy + shear_energy;
        // Factored stage derivative retains the admissible side at the fold.
        let curvature = (2.0 * e * kin.opening / (1.0 + 2.0 * e * a / self.kn)) * kin.derivative;
        let diagonal = [2.0 * e * a, k * a, k * a];
        let springs = [self.kn, self.kt, self.kt];
        let cross = [-2.0 * e * kin.opening, k * kin.shift[0], k * kin.shift[1]];
        let b: [f64; 3] =
            core::array::from_fn(|i| (springs[i] / (diagonal[i] + springs[i])) * cross[i]);
        let tangent = if marginal {
            None
        } else {
            let matrix = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    (if i == j {
                        (springs[i] / (diagonal[i] + springs[i])) * diagonal[i]
                    } else {
                        0.0
                    }) - b[i] * b[j] / curvature
                })
            });
            Some(matrix)
        };
        let forces = [
            contact.normal_force_n,
            contact.tangential_force_n[0],
            contact.tangential_force_n[1],
        ];
        let tip = [kin.delta, kin.shift[0], kin.shift[1]];
        let stages = [kin.stage, self.stage[0], self.stage[1]];
        let residual = core::array::from_fn(|i| forces[i] - springs[i] * (stages[i] - tip[i]));
        for value in [capacity, normal_energy, shear_energy, total, curvature]
            .into_iter()
            .chain(residual)
            .chain(tangent.into_iter().flatten().flatten())
        {
            finite(value)?;
        }
        Ok(JkrPullOffResponse {
            state: JkrPullOffState {
                identity: self.identity,
                branch: JkrPullOffBranch::Attached,
                normal_stage_m: kin.stage,
                released_energy_j: 0.0,
            },
            contact: Some(contact),
            tip_coordinates_m: tip,
            stage_force_n: forces,
            normal_spring_energy_j: normal_energy,
            tangential_spring_energy_j: shear_energy,
            total_stored_energy_j: total,
            static_capacity_n: capacity,
            static_margin_n: capacity - shear_force,
            apparatus_radius_curvature_n_per_m: curvature,
            stage_stiffness_n_per_m: tangent,
            spring_force_residual_n: residual,
        })
    }
}

// Put L=s^2=C a-D/(a+b)^2, q=1+alpha a. The fold ratio is
// 4 a q sqrt(L)/(R M), M=2 alpha L+q L'. Here L'>0, L''<=0.
// With r=a L'/L>0 and y=alpha a/q in [0,1), its logarithmic
// derivative times a exceeds 1+y+r/2-3 r y/(r+2y). The numerator
// after multiplying by 2(r+2y) is r^2+2r(1-y)+4y+4y^2>0.
// Thus this scalar fold is unique; no numerical branch scan is inferred.
fn finite(value: f64) -> Result<(), JkrPullOffError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(JkrPullOffError::Invalid("finite arithmetic"))
    }
}
