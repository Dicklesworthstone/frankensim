//! Absolute elastic/surface energy and the noncontact Dugdale disk.
//!
//! Uses the SAME W, z0, E*, radius, source and envelope as the contacting owner.
//! The surface potential is sigma0*(gap-z0) for 0 <= gap <= z0 and zero beyond.
//! Noncontact deformation is the Boussinesq response of a uniformly loaded disk,
//! not an undeformed-sphere DMT tail. Contact energy follows the equivalent Abel
//! line-load energy and reciprocal work with that disk. No energy quadrature is
//! needed. These are nominal half-space estimates, not stability certificates
//! against nonaxisymmetric perturbations or a dynamic solution of snap events.
use super::{MaugisError, MaugisResponse, PointNormalPort, SourceMaugisSphereContact, kernel};
use core::f64::consts::PI;
use fs_blake3::ContentHash;
use fs_tribo::{InputAuthority, TriboError, adhesive::AdhesiveContactState};

mod history {
    //! Displacement-controlled metastable continuation of the original Maugis law.
    //! The transition is a proposal: the caller publishes its successor only after
    //! the complete mechanical transaction accepts. Jump release is not caloric heat.
    use super::*;
    use super::super::{MaugisWorkBudget, MaugisWorkEstimate};

    /// Which loss of stability caused a fixed-displacement transition.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum MaugisSnapKind {
        /// The separated state reaches the finite interaction range.
        AttractionOnset,
        /// A contacting or noncontacting branch reaches its displacement fold.
        Separation,
    }

    /// Localized fixed-displacement event. Its potential decrease is an explicit
    /// energy output, not automatically assigned to heat, damage or radiation.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct MaugisSnap {
        kind: MaugisSnapKind,
        indentation_m: f64,
        before: MaugisEquilibrium,
        after: MaugisEquilibrium,
        released_energy_j: f64,
    }
    impl MaugisSnap {
        #[must_use] pub const fn kind(&self) -> MaugisSnapKind { self.kind }
        #[must_use] pub const fn indentation_m(&self) -> f64 { self.indentation_m }
        #[must_use] pub const fn before(&self) -> &MaugisEquilibrium { &self.before }
        #[must_use] pub const fn after(&self) -> &MaugisEquilibrium { &self.after }
        #[must_use] pub const fn released_energy_j(&self) -> f64 { self.released_energy_j }
    }

    /// Immutable source-bound accepted history. The requested control coordinate is
    /// distinct from the contact root's reconstructed indentation (roundoff only).
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct MaugisDisplacementHistory {
        model_identity: ContentHash,
        control_indentation_m: f64,
        equilibrium: MaugisEquilibrium,
    }
    impl MaugisDisplacementHistory {
        #[must_use] pub const fn indentation_m(&self) -> f64 { self.control_indentation_m }
        #[must_use] pub const fn equilibrium(&self) -> &MaugisEquilibrium { &self.equilibrium }
        #[must_use] pub const fn model_identity(&self) -> ContentHash { self.model_identity }
    }

    /// Separate force-work integration and independent full-energy closure budgets.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct MaugisContinuationBudget {
        pub work: MaugisWorkBudget,
        pub energy_absolute_tolerance_j: f64,
        pub energy_relative_tolerance: f64,
    }
    impl MaugisContinuationBudget {
        fn validate(self) -> Result<(), MaugisContinuationError> {
            if !self.work.absolute_tolerance_j.is_finite() || self.work.absolute_tolerance_j <= 0.0
                || !self.work.relative_tolerance.is_finite() || !(0.0..1.0).contains(&self.work.relative_tolerance)
                || self.work.maximum_evaluations < 5
                || !self.energy_absolute_tolerance_j.is_finite() || self.energy_absolute_tolerance_j < 0.0
                || !self.energy_relative_tolerance.is_finite() || !(0.0..1.0).contains(&self.energy_relative_tolerance)
            {
                return Err(MaugisError::InvalidInput { field: "Maugis continuation work/energy budget" }.into());
            }
            Ok(())
        }
    }

    /// An accepted numerical proposal; no state inside the owner was advanced.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct MaugisDisplacementStep {
        pub next: MaugisDisplacementHistory,
        pub snap: Option<MaugisSnap>,
        /// Integral of original force over continuous portions only.
        pub work: MaugisWorkEstimate,
        pub potential_change_j: f64,
        pub released_energy_j: f64,
        /// work - potential change - released energy, independently reconstructed.
        pub energy_residual_j: f64,
        pub energy_tolerance_j: f64,
    }

    /// Prepared source-preserving continuation. Fold locations, jump states and
    /// their energies are computed once, not rediscovered in every mechanical step.
    /// This is massless quasistatic DISPLACEMENT control, not force control or a
    /// trajectory through the jump. Finite fixture compliance has a different fold.
    #[derive(Clone, Debug, PartialEq)]
    pub struct MaugisDisplacementContinuation {
        law: SourceMaugisSphereContact,
        noncontact: Option<MaugisNoncontactLimits>,
        onset: MaugisSnap,
        separation: MaugisSnap,
    }
    impl MaugisDisplacementContinuation {
        pub fn new(law: SourceMaugisSphereContact) -> Result<Self, MaugisContinuationError> {
            let noncontact = law.noncontact_limits()?;
            let onset_d = -law.interaction_range_m();
            let before = law.evaluate_detached(onset_d)?;
            let after = match noncontact {
                Some(limits) if onset_d < limits.first_contact_indentation_m => law.evaluate_noncontact(onset_d)?,
                _ => law.evaluate_contact_energy(onset_d)?,
            };
            let onset = snap(MaugisSnapKind::AttractionOnset, onset_d, before, after)?;
            let (separation_d, before) = match noncontact {
                Some(limits) => (limits.minimum_indentation_m, law.evaluate_noncontact(limits.minimum_indentation_m)?),
                None => {
                    let d = law.minimum_attached_indentation_m();
                    let before = law.evaluate_contact_energy(d)?;
                    if before.contact_response().is_some_and(|r| r.contact_radius_m == 0.0) {
                        return Err(MaugisError::Unrepresentable { field: "unresolved contact/noncontact stability crossover" }.into());
                    }
                    (d, before)
                }
            };
            if separation_d >= onset_d {
                return Err(MaugisError::Unrepresentable { field: "ordered Maugis stability limits" }.into());
            }
            let after = law.evaluate_detached(separation_d)?;
            let separation = snap(MaugisSnapKind::Separation, separation_d, before, after)?;
            Ok(Self { law, noncontact, onset, separation })
        }
        /// Original full source owner, not a fitted replacement potential.
        #[must_use] pub const fn law(&self) -> &SourceMaugisSphereContact { &self.law }
        #[must_use] pub const fn onset(&self) -> &MaugisSnap { &self.onset }
        #[must_use] pub const fn separation(&self) -> &MaugisSnap { &self.separation }

        /// Choose the initial physical branch explicitly. No unreported jump or
        /// previous loading history is inferred from an initial displacement.
        pub fn initialize(&self, phase: MaugisPhase, indentation_m: f64)
            -> Result<MaugisDisplacementHistory, MaugisContinuationError>
        {
            let equilibrium = match phase {
                MaugisPhase::Detached => self.law.evaluate_detached(indentation_m)?,
                MaugisPhase::Noncontact => self.law.evaluate_noncontact(indentation_m)?,
                MaugisPhase::Contact => self.law.evaluate_contact_energy(indentation_m)?,
            };
            Ok(MaugisDisplacementHistory { model_identity: self.law.identity(), control_indentation_m: indentation_m, equilibrium })
        }

        /// Follow one monotone segment, preserving metastability in a coexistence
        /// interval. Reversals require separate calls. A crossed endpoint is consumed
        /// at equality; a zero-length proposal never emits a new event. Retrying from
        /// the same token is deterministic; adopting `next` prevents repeated release.
        /// The caller remains responsible for publishing each proposal only once.
        pub fn advance(&self, old: &MaugisDisplacementHistory, target_m: f64, budget: MaugisContinuationBudget)
            -> Result<MaugisDisplacementStep, MaugisContinuationError>
        {
            budget.validate()?;
            if old.model_identity != self.law.identity() || old.equilibrium.model_identity() != self.law.identity() {
                return Err(MaugisContinuationError::HistoryMismatch);
            }
            if !target_m.is_finite() {
                return Err(MaugisError::InvalidInput { field: "continuation target indentation" }.into());
            }
            let from = old.control_indentation_m;
            let detached = old.equilibrium.phase() == MaugisPhase::Detached;
            let event = if target_m > from && detached && target_m >= self.onset.indentation_m {
                Some(self.onset)
            } else if target_m < from && !detached && target_m <= self.separation.indentation_m {
                Some(self.separation)
            } else { None };
            let (equilibrium, work) = if target_m == from {
                (old.equilibrium, zero_work())
            } else if let Some(event) = event {
                match event.kind {
                    MaugisSnapKind::AttractionOnset => (
                        self.connected_equilibrium(target_m)?,
                        self.connected_work(event.indentation_m, target_m, budget.work)?,
                    ),
                    MaugisSnapKind::Separation => (
                        self.law.evaluate_detached(target_m)?,
                        self.connected_work(from, event.indentation_m, budget.work)?,
                    ),
                }
            } else if detached {
                (self.law.evaluate_detached(target_m)?, zero_work())
            } else {
                (self.connected_equilibrium(target_m)?, self.connected_work(from, target_m, budget.work)?)
            };
            let released = event.map_or(0.0, |e| e.released_energy_j);
            let change = equilibrium.potential_j()-old.equilibrium.potential_j();
            let residual = (work.work_j-change)-released;
            let scale = old.equilibrium.potential_j().abs().max(equilibrium.potential_j().abs())
                .max(work.work_j.abs()).max(released);
            let tolerance = budget.energy_absolute_tolerance_j+budget.energy_relative_tolerance*scale;
            if [change, residual, tolerance].iter().any(|v| !v.is_finite()) || residual.abs() > tolerance {
                return Err(MaugisContinuationError::EnergyClosure { residual_j: residual, tolerance_j: tolerance });
            }
            Ok(MaugisDisplacementStep {
                next: MaugisDisplacementHistory { model_identity: self.law.identity(), control_indentation_m: target_m, equilibrium },
                snap: event, work, potential_change_j: change, released_energy_j: released,
                energy_residual_j: residual, energy_tolerance_j: tolerance,
            })
        }
        fn connected_equilibrium(&self, d: f64) -> Result<MaugisEquilibrium, MaugisContinuationError> {
            match self.noncontact {
                Some(limits) if d < limits.first_contact_indentation_m => self.law.evaluate_noncontact(d),
                _ => self.law.evaluate_contact_energy(d),
            }
        }
        fn connected_work(&self, from: f64, to: f64, budget: MaugisWorkBudget)
            -> Result<MaugisWorkEstimate, MaugisContinuationError>
        {
            if from == to { return Ok(zero_work()); }
            let Some(limits) = self.noncontact else {
                return self.law.reversible_work_between(from, to, budget).map_err(Into::into);
            };
            let sign = if to > from { 1.0 } else { -1.0 };
            let lo = from.min(to);
            let hi = from.max(to);
            let touch = limits.first_contact_indentation_m;
            let mut work = zero_work();
            if lo < touch {
                work.work_j = self.noncontact_work(lo, hi.min(touch))?;
            }
            if hi > touch {
                let contact = self.law.reversible_work_between(lo.max(touch), hi, budget)?;
                work.work_j += contact.work_j;
                work.estimated_absolute_error_j = contact.estimated_absolute_error_j;
                work.evaluations = contact.evaluations;
            }
            work.work_j *= sign;
            if !work.work_j.is_finite() { return Err(MaugisError::Unrepresentable { field: "continuation path work" }.into()); }
            Ok(work)
        }
        fn noncontact_work(&self, from: f64, to: f64) -> Result<f64, MaugisContinuationError> {
            let a = self.law.evaluate_noncontact(from)?;
            let b = self.law.evaluate_noncontact(to)?;
            if from == to { return Ok(0.0); }
            let c0 = a.cohesive_radius_m()/self.law.radius_scale_m;
            let c1 = b.cohesive_radius_m()/self.law.radius_scale_m;
            let fold = (8.0/(3.0*PI))*self.law.lambda;
            let x0 = c0-fold;
            let x1 = c1-fold;
            // Exact integral P(c)*(c-fold) dc. Translate to x=c-fold so the
            // polynomial terms have one sign and remain well conditioned at fold.
            // Rationalize x1-x0 from the displacement difference, rather than
            // subtracting nearly equal radii (or absolute endpoint energies).
            let dx = (2.0*((to-from)/self.law.indentation_scale_m))/(x0+x1);
            let i1 = 0.5*(x0+x1);
            let i2 = (x0*x0+x0*x1+x1*x1)/3.0;
            let i3 = 0.25*(x0+x1)*(x0*x0+x1*x1);
            let value = -(PI/2.0)*self.law.lambda*dx*(fold*fold*i1+2.0*fold*i2+i3)
                * self.law.continuation_energy_scale()?;
            if !dx.is_finite() || dx <= 0.0 || !value.is_finite() || value == 0.0 {
                return Err(MaugisError::Unrepresentable { field: "noncontact work increment" }.into());
            }
            Ok(value)
        }
    }
    fn zero_work() -> MaugisWorkEstimate {
        MaugisWorkEstimate { work_j: 0.0, estimated_absolute_error_j: 0.0, evaluations: 0 }
    }
    fn snap(kind: MaugisSnapKind, indentation_m: f64, before: MaugisEquilibrium, after: MaugisEquilibrium)
        -> Result<MaugisSnap, MaugisContinuationError>
    {
        let released = before.potential_j()-after.potential_j();
        if !released.is_finite() || released < 0.0 {
            return Err(MaugisContinuationError::NonpassiveSnap { released_energy_j: released });
        }
        Ok(MaugisSnap { kind, indentation_m, before, after, released_energy_j: released })
    }
}
pub use history::{MaugisContinuationBudget, MaugisDisplacementContinuation,
    MaugisDisplacementHistory, MaugisDisplacementStep, MaugisSnap, MaugisSnapKind};

/// Physical branch; attraction does not imply real solid contact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaugisPhase {
    Detached,
    Noncontact,
    Contact,
}

/// Original numerical/source refusal or an unavailable continuation branch.
#[derive(Clone, Debug, PartialEq)]
pub enum MaugisContinuationError {
    Contact(MaugisError),
    HistoryMismatch,
    NonpassiveSnap { released_energy_j: f64 },
    EnergyClosure { residual_j: f64, tolerance_j: f64 },
    NoNoncontactEquilibrium { indentation_m: f64 },
    NoDetachedEquilibrium { indentation_m: f64, cutoff_m: f64 },
}
impl From<MaugisError> for MaugisContinuationError {
    fn from(error: MaugisError) -> Self { Self::Contact(error) }
}
impl core::fmt::Display for MaugisContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Maugis continuation: {self:?}")
    }
}
impl std::error::Error for MaugisContinuationError {}

/// Absolute energy uses fully separated, undeformed surfaces as zero.
/// Fields are read-only: a continuation token cannot import edited energy data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisEquilibrium {
    model_identity: ContentHash,
    phase: MaugisPhase,
    indentation_m: f64,
    normal_force_n: f64,
    cohesive_radius_m: f64,
    central_gap_m: f64,
    elastic_energy_j: f64,
    surface_energy_j: f64,
    potential_j: f64,
    tangent_n_per_m: Option<f64>,
    authority: InputAuthority,
    contact: Option<MaugisResponse>,
}
impl MaugisEquilibrium {
    #[must_use] pub const fn model_identity(&self) -> ContentHash { self.model_identity }
    #[must_use] pub const fn phase(&self) -> MaugisPhase { self.phase }
    #[must_use] pub const fn indentation_m(&self) -> f64 { self.indentation_m }
    #[must_use] pub const fn normal_force_n(&self) -> f64 { self.normal_force_n }
    #[must_use] pub const fn cohesive_radius_m(&self) -> f64 { self.cohesive_radius_m }
    #[must_use] pub const fn central_gap_m(&self) -> f64 { self.central_gap_m }
    #[must_use] pub const fn elastic_energy_j(&self) -> f64 { self.elastic_energy_j }
    #[must_use] pub const fn surface_energy_j(&self) -> f64 { self.surface_energy_j }
    #[must_use] pub const fn potential_j(&self) -> f64 { self.potential_j }
    #[must_use] pub const fn tangent_n_per_m(&self) -> Option<f64> { self.tangent_n_per_m }
    #[must_use] pub const fn authority(&self) -> InputAuthority { self.authority }
    /// Unmodified contacting response for the existing source-bound shear adapter.
    #[must_use] pub const fn contact_response(&self) -> Option<&MaugisResponse> { self.contact.as_ref() }
    /// The entire attractive disk carries ZERO junction shear while noncontacting.
    pub fn friction_contact_state(&self) -> Result<AdhesiveContactState, TriboError> {
        match &self.contact {
            Some(contact) => contact.friction_contact_state(),
            None => AdhesiveContactState::new(0.0, 0.0),
        }
    }
    /// Signed reciprocal wrench; folds and branch junctions refuse a smooth port.
    pub fn point_port(&self, normal: [f64; 3], arm_m: [f64; 3]) -> Result<PointNormalPort, MaugisError> {
        let tangent = self.tangent_n_per_m.ok_or(MaugisError::EndpointTangent)?;
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !length.is_finite() || (length - 1.0).abs() > 1e-12 || arm_m.iter().any(|x| !x.is_finite()) {
            return Err(MaugisError::InvalidInput { field: "Maugis continuation port frame" });
        }
        let f = normal.map(|x| x * self.normal_force_n);
        let m = [arm_m[1]*f[2]-arm_m[2]*f[1], arm_m[2]*f[0]-arm_m[0]*f[2], arm_m[0]*f[1]-arm_m[1]*f[0]];
        if f.iter().chain(m.iter()).any(|x| !x.is_finite()) {
            return Err(MaugisError::Unrepresentable { field: "Maugis continuation wrench" });
        }
        Ok(PointNormalPort { action_force_n: f, action_moment_n_m: m,
            reaction_force_n: f.map(|x| -x), reaction_moment_n_m: m.map(|x| -x),
            residual_force_n: [0.0; 3], tangent_n_per_m: tangent,
            dissipated_power_w: 0.0, irreversible_work_j: 0.0 })
    }
}

/// The finite stable noncontact interval, including its marginal lower endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisNoncontactLimits {
    pub minimum_indentation_m: f64,
    pub first_contact_indentation_m: f64,
    pub minimum_cohesive_radius_m: f64,
    pub first_contact_cohesive_radius_m: f64,
}

// theta - sin(k*theta)/k. Fixed work, stable even for a very thin annulus.
fn angular_defect(theta: f64, k: f64) -> f64 {
    if theta >= 0.125 { return theta - (k * theta).sin() / k; }
    let x2 = (k * theta) * (k * theta);
    let mut term = theta * x2 / 6.0;
    let mut sum = term;
    for n in 2..=10 {
        term *= -x2 / f64::from((2*n)*(2*n+1));
        sum += term;
    }
    sum
}

// Energies divided by (pi W R)*(a_scale^2/R). Abel line-load q(x) is
// E*(delta-x^2/R) on [0,a], -2*sigma0*sqrt(c^2-x^2) on [a,c].
// Elastic energy = integral(q^2/E*) dx. The surface part follows reciprocal
// disk work; J = integral_a^c(c^2-x^2)dx is evaluated without c-a cancellation.
fn contact_energies(lambda: f64, s: kernel::State) -> (f64, f64, f64) {
    let a2 = s.a * s.a;
    let c2 = s.c * s.c;
    let width = s.l * (s.l / (s.c + s.a));
    let j = width * width * ((2.0*s.c+s.a)/3.0);
    let shift = s.d - a2/3.0;
    let inside = 0.75*s.a*(shift*shift + (4.0/45.0)*a2*a2);
    let outside = (4.0/3.0)*lambda*lambda*j;
    let theta = 1.0_f64.atan2(s.t);
    let surface_base = -c2 - lambda*s.d*c2*angular_defect(theta, 2.0)
        + 0.25*lambda*c2*c2*angular_defect(theta, 4.0);
    (inside+outside, surface_base-2.0*outside, inside-outside+surface_base)
}

impl SourceMaugisSphereContact {
    /// Contact branch with its absolute elastic-plus-surface potential.
    /// dU/d(delta)=F, with no numerical integration or reference-state fitting.
    pub fn evaluate_contact_energy(&self, indentation_m: f64) -> Result<MaugisEquilibrium, MaugisContinuationError> {
        let state = self.state_at_indentation(indentation_m)?;
        let contact = self.response(state)?;
        let (elastic, surface, potential) = contact_energies(self.lambda, state);
        let scale = self.continuation_energy_scale()?;
        let result = MaugisEquilibrium { model_identity: self.identity(), phase: MaugisPhase::Contact,
            indentation_m: contact.indentation_m, normal_force_n: contact.normal_force_n,
            cohesive_radius_m: contact.cohesive_radius_m, central_gap_m: 0.0,
            elastic_energy_j: elastic*scale, surface_energy_j: surface*scale, potential_j: potential*scale,
            tangent_n_per_m: contact.tangent_n_per_m, authority: contact.authority, contact: Some(contact) };
        check_energy(&result)?;
        Ok(result)
    }

    /// Stable noncontact states exist only while the disk's center gap is >=0
    /// and c >= 4*sigma0*R/(pi E*). None means contact intervenes before that fold.
    pub fn noncontact_limits(&self) -> Result<Option<MaugisNoncontactLimits>, MaugisContinuationError> {
        let zero = kernel::state(self.lambda, 0.0)?;
        let fold = (8.0/(3.0*PI))*self.lambda;
        if zero.c <= fold { return Ok(None); }
        let minimum = -self.interaction_range_m - (0.5*fold*self.indentation_scale_m)*fold;
        let first = zero.d*self.indentation_scale_m;
        if !minimum.is_finite() || minimum >= -self.interaction_range_m || minimum >= first {
            return Err(MaugisError::Unrepresentable { field: "distinct noncontact fold and cutoff" }.into());
        }
        Ok(Some(MaugisNoncontactLimits { minimum_indentation_m: minimum, first_contact_indentation_m: first,
            minimum_cohesive_radius_m: fold*self.radius_scale_m,
            first_contact_cohesive_radius_m: zero.c*self.radius_scale_m }))
    }

    /// Exact uniform-traction disk on the locally stable axisymmetric branch.
    /// delta = c^2/(2R)-4*sigma0*c/(pi E*)-z0, F=-pi*sigma0*c^2.
    /// Choosing a branch is explicit: this does not silently attach a far body.
    pub fn evaluate_noncontact(&self, indentation_m: f64) -> Result<MaugisEquilibrium, MaugisContinuationError> {
        if !indentation_m.is_finite() {
            return Err(MaugisError::InvalidInput { field: "noncontact indentation" }.into());
        }
        let limits = self.noncontact_limits()?.ok_or(MaugisContinuationError::NoNoncontactEquilibrium { indentation_m })?;
        if indentation_m < limits.minimum_indentation_m || indentation_m > limits.first_contact_indentation_m {
            return Err(MaugisContinuationError::NoNoncontactEquilibrium { indentation_m });
        }
        let zero = kernel::state(self.lambda, 0.0)?;
        let fold = (8.0/(3.0*PI))*self.lambda;
        let offset = (indentation_m-limits.minimum_indentation_m)/self.indentation_scale_m;
        let c = if indentation_m == limits.first_contact_indentation_m { zero.c }
            else { fold+(2.0*offset).sqrt() };
        if !c.is_finite() || c < fold || c > zero.c {
            return Err(MaugisError::Unrepresentable { field: "stable noncontact radius" }.into());
        }
        let radius = c*self.radius_scale_m;
        self.admit_noncontact_radius(radius)?;
        let c2 = c*c;
        let c3 = c2*c;
        let elastic = (8.0/9.0)*self.lambda*self.lambda*c3;
        let geometric = (PI/8.0)*self.lambda*c2*c2;
        let surface = -geometric-0.5*elastic;
        let potential = -geometric+0.5*elastic;
        let gap = (zero.c-c)*(0.5*(zero.c+c)+(4.0/3.0-8.0/(3.0*PI))*self.lambda)*self.indentation_scale_m;
        let endpoint = indentation_m == limits.minimum_indentation_m || indentation_m == limits.first_contact_indentation_m;
        let tangent = if endpoint { None } else {
            Some(-PI*self.lambda*(c/(c-fold))*(self.force_scale_n/self.indentation_scale_m))
        };
        let scale = self.continuation_energy_scale()?;
        let result = MaugisEquilibrium { model_identity: self.identity(), phase: MaugisPhase::Noncontact,
            indentation_m, normal_force_n: -(PI/2.0)*self.lambda*c2*self.force_scale_n,
            cohesive_radius_m: radius, central_gap_m: gap,
            elastic_energy_j: elastic*scale, surface_energy_j: surface*scale, potential_j: potential*scale,
            tangent_n_per_m: tangent, authority: self.continuation_authority(), contact: None };
        check_energy(&result)?;
        if radius <= 0.0 || !gap.is_finite() || gap < 0.0 || result.normal_force_n >= 0.0
            || tangent.is_some_and(|v| !v.is_finite()) {
            return Err(MaugisError::Unrepresentable { field: "dimensional noncontact response" }.into());
        }
        Ok(result)
    }

    /// Undeformed separated branch; its nearest point must lie beyond z0.
    /// At the exact cutoff the smooth tangent is withheld.
    pub fn evaluate_detached(&self, indentation_m: f64) -> Result<MaugisEquilibrium, MaugisContinuationError> {
        if !indentation_m.is_finite() {
            return Err(MaugisError::InvalidInput { field: "detached indentation" }.into());
        }
        let cutoff_m = -self.interaction_range_m;
        if indentation_m > cutoff_m {
            return Err(MaugisContinuationError::NoDetachedEquilibrium { indentation_m, cutoff_m });
        }
        Ok(MaugisEquilibrium { model_identity: self.identity(), phase: MaugisPhase::Detached,
            indentation_m, normal_force_n: 0.0, cohesive_radius_m: 0.0, central_gap_m: -indentation_m,
            elastic_energy_j: 0.0, surface_energy_j: 0.0, potential_j: 0.0,
            tangent_n_per_m: if indentation_m == cutoff_m { None } else { Some(0.0) },
            authority: self.continuation_authority(), contact: None })
    }

    fn continuation_energy_scale(&self) -> Result<f64, MaugisError> {
        let scale = self.force_scale_n*self.indentation_scale_m;
        if !scale.is_finite() || scale <= 0.0 {
            Err(MaugisError::Unrepresentable { field: "absolute Maugis energy scale" })
        } else { Ok(scale) }
    }
    fn continuation_authority(&self) -> InputAuthority {
        match self.normal.elastic().interface().interface().provenance().authority() {
            InputAuthority::SyntheticFixture => InputAuthority::SyntheticFixture,
            _ => InputAuthority::Estimated,
        }
    }
    fn admit_noncontact_radius(&self, radius: f64) -> Result<(), MaugisError> {
        for (ratio, value, limit) in [
            ("noncontact radius/radius", radius/self.radius_m, self.envelope.maximum_interaction_radius_to_radius),
            ("noncontact radius/depth", radius/self.envelope.half_space_depth_m, self.envelope.maximum_interaction_radius_to_depth),
            ("noncontact radius/layer", radius/self.envelope.layer_thickness_m, self.envelope.maximum_interaction_radius_to_layer),
        ] {
            if !value.is_finite() || value > limit { return Err(MaugisError::OutsideApplicability { ratio, value, limit }); }
        }
        Ok(())
    }
}
fn check_energy(result: &MaugisEquilibrium) -> Result<(), MaugisError> {
    if [result.elastic_energy_j, result.surface_energy_j, result.potential_j,
        result.normal_force_n, result.cohesive_radius_m, result.indentation_m, result.central_gap_m].iter().any(|v| !v.is_finite())
        || result.elastic_energy_j < 0.0 || result.surface_energy_j > 0.0 {
        Err(MaugisError::Unrepresentable { field: "absolute elastic/surface energy" })
    } else { Ok(()) }
}
