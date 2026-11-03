//! Finite loading compliance on the original Maugis--Dugdale sphere law.
//!
//! A massless tip obeys F(delta)=K(D-delta). In the original dimensionless
//! kernel, D/d0=d(t)+chi*p(t), chi=F0/(K*d0); stability ends where
//! d'(t)+chi*p'(t)=0. Thus a finite spring has its own fold, between the
//! displacement- and force-controlled limits. No JKR/DMT interpolation is used.
//!
//! For the deformed noncontact disk, D/d0=A*c²/2-beta*c-z0/d0, where
//! A=1-pi*lambda*chi and beta=8*lambda/(3*pi). A stable disk exists only
//! for A>0, c>=beta/A and a nonnegative center gap. All energies retain the
//! original separated reference, adding the actual spring storage F²/(2K).
//! Snap release is mechanical energy leaving this quasistatic model, not heat.
//! These equations follow directly by adding a linear loading spring to the
//! Maugis and noncontact equations cited by the parent and continuation modules.

use super::continuation::{
    MaugisContinuationBudget, MaugisContinuationError, MaugisEquilibrium, MaugisPhase,
    MaugisSnapKind,
};
use super::{MaugisError, MaugisWorkBudget, MaugisWorkEstimate, SourceMaugisSphereContact, kernel};
use crate::normal_patch::PointNormalPort;
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};

/// A refusal never changes the accepted fixture history.
#[derive(Clone, Debug, PartialEq)]
pub enum MaugisFixtureError {
    /// Original material/domain, contact-kernel or work-quadrature refusal.
    Contact(MaugisError),
    /// Original absolute-energy or deformed noncontact-branch refusal.
    Continuation(MaugisContinuationError),
    /// The requested branch has no stable equilibrium at this stage position.
    NoStableEquilibrium {
        /// Explicit physical phase requested by the caller.
        phase: MaugisPhase,
        /// Requested stage approach, m.
        stage_m: f64,
    },
    /// The history belongs to a different original law or loading stiffness.
    HistoryMismatch,
    /// Independent force-work and total-potential accounting disagree.
    EnergyClosure {
        /// Stage work minus potential change and jump release, J.
        residual_j: f64,
        /// Caller-specified absolute-plus-relative allowance, J.
        tolerance_j: f64,
    },
}
impl From<MaugisError> for MaugisFixtureError {
    fn from(error: MaugisError) -> Self {
        Self::Contact(error)
    }
}
impl From<MaugisContinuationError> for MaugisFixtureError {
    fn from(error: MaugisContinuationError) -> Self {
        Self::Continuation(error)
    }
}
impl core::fmt::Display for MaugisFixtureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Maugis loading fixture refused: {self:?}")
    }
}
impl std::error::Error for MaugisFixtureError {}

/// An immutable equilibrium of the original contact plus loading spring.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisFixtureResponse {
    equilibrium: MaugisEquilibrium,
    contact_coordinate: Option<f64>,
    stage_m: f64,
    spring_energy_j: f64,
    total_potential_j: f64,
    stage_tangent_n_per_m: Option<f64>,
    equilibrium_residual_n: f64,
}
impl MaugisFixtureResponse {
    /// Unmodified original source, real contact area and elastic/surface energy.
    #[must_use]
    pub const fn equilibrium(&self) -> &MaugisEquilibrium {
        &self.equilibrium
    }
    /// Prescribed loading-stage approach D, m; distinct from tip indentation.
    #[must_use]
    pub const fn stage_displacement_m(&self) -> f64 {
        self.stage_m
    }
    /// Nonnegative loading-spring elastic storage F²/(2K), J.
    #[must_use]
    pub const fn spring_energy_j(&self) -> f64 {
        self.spring_energy_j
    }
    /// Original contact potential plus loading-spring storage, J.
    #[must_use]
    pub const fn total_potential_j(&self) -> f64 {
        self.total_potential_j
    }
    /// dF/dD, withheld at fixture folds and geometric branch junctions.
    #[must_use]
    pub const fn stage_tangent_n_per_m(&self) -> Option<f64> {
        self.stage_tangent_n_per_m
    }
    /// Original F-K(D-delta), without projecting either physical quantity.
    #[must_use]
    pub const fn equilibrium_residual_n(&self) -> f64 {
        self.equilibrium_residual_n
    }
    /// Reciprocal resultant conjugate to STAGE motion; tip ports are distinct.
    pub fn stage_port(
        &self,
        normal: [f64; 3],
        arm_m: [f64; 3],
    ) -> Result<PointNormalPort, MaugisFixtureError> {
        let tangent = self
            .stage_tangent_n_per_m
            .ok_or(MaugisError::EndpointTangent)?;
        let mut port = self.equilibrium.point_port(normal, arm_m)?;
        port.tangent_n_per_m = tangent;
        Ok(port)
    }
}

/// Accepted source-and-stiffness-bound branch history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisFixtureState {
    fixture_identity: ContentHash,
    response: MaugisFixtureResponse,
}
impl MaugisFixtureState {
    /// Accepted equilibrium and spring state at the retained stage coordinate.
    #[must_use]
    pub const fn response(&self) -> &MaugisFixtureResponse {
        &self.response
    }
    /// Content identity of original sources, geometry, envelope and stiffness.
    #[must_use]
    pub const fn fixture_identity(&self) -> ContentHash {
        self.fixture_identity
    }
}

/// Fixed-stage instability, including finite tip motion and spring release.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisFixtureJump {
    kind: MaugisSnapKind,
    before: MaugisFixtureResponse,
    after: MaugisFixtureResponse,
    released_energy_j: f64,
}
impl MaugisFixtureJump {
    /// Attraction-onset snap or withdrawal separation.
    #[must_use]
    pub const fn kind(&self) -> MaugisSnapKind {
        self.kind
    }
    /// Fixed loading-stage position at the instability, m.
    #[must_use]
    pub const fn stage_displacement_m(&self) -> f64 {
        self.before.stage_m
    }
    /// Marginal pre-jump contact and spring equilibrium.
    #[must_use]
    pub const fn before(&self) -> &MaugisFixtureResponse {
        &self.before
    }
    /// Post-jump equilibrium at exactly the same prescribed stage.
    #[must_use]
    pub const fn after(&self) -> &MaugisFixtureResponse {
        &self.after
    }
    /// Positive total-potential decrease, J; no thermal destination is inferred.
    #[must_use]
    pub const fn released_energy_j(&self) -> f64 {
        self.released_energy_j
    }
}

/// One monotone stage proposal; publish `next` only when the outer step accepts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisFixtureTransition {
    /// Complete immutable successor, adopted only when the outer step accepts.
    pub next: MaugisFixtureState,
    /// At most one localized instability on a monotone stage segment.
    pub jump: Option<MaugisFixtureJump>,
    /// Integral F dD along continuous branches; never includes snap motion.
    pub work: MaugisWorkEstimate,
    /// Endpoint minus initial total contact-and-spring potential, J.
    pub potential_change_j: f64,
    /// Jump release, J; zero on a smooth segment.
    pub released_energy_j: f64,
    /// Continuous stage work minus potential change and jump release, J.
    pub energy_residual_j: f64,
    /// Supplied absolute-plus-relative energy allowance, J.
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Noncontact {
    quadratic: f64,
    fold_radius: f64,
    minimum_stage_m: f64,
    touch_stage_m: f64,
    touch_indentation_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
struct Branches {
    law: SourceMaugisSphereContact,
    stiffness: f64,
    compliance: f64,
    minimum_t: f64,
    minimum_contact_stage_m: f64,
    maximum_stage_m: f64,
    noncontact: Option<Noncontact>,
}

/// Prepared finite-stiffness continuation, with both instability states cached.
/// Root inversions have at most 128 iterations and allocate no memory. Work
/// uses the original bounded quadrature, with analytic noncontact integration.
#[derive(Clone, Debug, PartialEq)]
pub struct MaugisSpringFixture {
    branches: Branches,
    identity: ContentHash,
    onset: MaugisFixtureJump,
    separation: MaugisFixtureJump,
}
impl MaugisSpringFixture {
    /// Bind K>0 [N/m]. The original envelope must include the snap-in landing.
    /// Unresolved fold/range/junction separations refuse instead of erasing them.
    pub fn new(
        law: SourceMaugisSphereContact,
        stiffness_n_per_m: f64,
    ) -> Result<Self, MaugisFixtureError> {
        super::positive(stiffness_n_per_m, "loading spring stiffness")?;
        let compliance = (law.force_scale_n / stiffness_n_per_m) / law.indentation_scale_m;
        let eta = PI * law.lambda * compliance;
        if !compliance.is_finite() || compliance <= 0.0 || !eta.is_finite() {
            return Err(unrepresentable("dimensionless loading compliance"));
        }
        let slope = |t| -> Result<f64, MaugisError> {
            let s = kernel::state(law.lambda, t)?;
            Ok(s.ds + compliance * s.ps)
        };
        let minimum_t = if slope(0.0)? >= 0.0 {
            0.0
        } else {
            let mut lo = law.minimum_t;
            let mut hi = law.load_fold_t;
            if hi > law.maximum_t || hi <= lo || slope(hi)? <= 0.0 {
                return Err(unrepresentable("finite-stiffness contacting fold bracket"));
            }
            for _ in 0..128 {
                let mid = 0.5 * lo + 0.5 * hi;
                if mid == lo || mid == hi {
                    break;
                }
                if slope(mid)? < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            hi
        };
        if minimum_t > 0.0 {
            let fold = kernel::state(law.lambda, minimum_t)?;
            let residual = fold.ds + compliance * fold.ps;
            // Source derivatives themselves involve cancellation near either
            // ideal-control fold, before ps is multiplied by compliance.
            let scale = fold.ds.abs().max(fold.ps.abs()).max(1.0) * (1.0 + compliance);
            let tolerance = 4096.0 * f64::EPSILON * scale;
            if minimum_t <= law.minimum_t
                || minimum_t >= law.load_fold_t
                || fold.ds <= 0.0
                || fold.ps >= 0.0
                || !residual.is_finite()
                || !tolerance.is_finite()
                || residual.abs() > tolerance
            {
                return Err(unrepresentable("resolved finite-stiffness contacting fold"));
            }
        }
        let stage = |s: kernel::State| (s.d + compliance * s.p) * law.indentation_scale_m;
        let minimum_contact_stage_m = stage(kernel::state(law.lambda, minimum_t)?);
        let maximum_stage_m = stage(kernel::state(law.lambda, law.maximum_t)?);
        if minimum_t > 0.0 && law.minimum_t == 0.0 {
            let touch = stage(kernel::state(law.lambda, 0.0)?);
            let resolution = 64.0 * f64::EPSILON * touch.abs().max(minimum_contact_stage_m.abs());
            if touch - minimum_contact_stage_m <= resolution {
                return Err(unrepresentable(
                    "distinct contacting fold and zero-area junction",
                ));
            }
        }
        let noncontact = if minimum_t == 0.0 {
            let quadratic = 1.0 - eta;
            let beta = (8.0 / (3.0 * PI)) * law.lambda;
            let fold_radius = beta / quadratic;
            let zero = kernel::state(law.lambda, 0.0)?;
            let minimum_stage_m =
                -law.interaction_range_m - 0.5 * beta * fold_radius * law.indentation_scale_m;
            let span = 0.5 * quadratic * (zero.c - fold_radius).powi(2) * law.indentation_scale_m;
            let resolution =
                64.0 * f64::EPSILON * minimum_stage_m.abs().max(minimum_contact_stage_m.abs());
            if quadratic <= 0.0
                || !fold_radius.is_finite()
                || fold_radius >= zero.c
                || !minimum_stage_m.is_finite()
                || minimum_stage_m >= -law.interaction_range_m
                || minimum_stage_m >= minimum_contact_stage_m
                || !span.is_finite()
                || span <= resolution
            {
                return Err(unrepresentable("stable noncontact fixture branch"));
            }
            Some(Noncontact {
                quadratic,
                fold_radius,
                minimum_stage_m,
                touch_stage_m: minimum_contact_stage_m,
                touch_indentation_m: zero.d * law.indentation_scale_m,
            })
        } else {
            None
        };
        let minimum = noncontact.map_or(minimum_contact_stage_m, |n| n.minimum_stage_m);
        if !minimum.is_finite()
            || !maximum_stage_m.is_finite()
            || minimum >= -law.interaction_range_m
            || maximum_stage_m <= -law.interaction_range_m
        {
            return Err(unrepresentable("ordered finite fixture stage bounds"));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.maugis-spring-fixture.v1");
        hash.update(law.identity().as_bytes());
        hash.update(&stiffness_n_per_m.to_bits().to_le_bytes());
        let branches = Branches {
            law,
            stiffness: stiffness_n_per_m,
            compliance,
            minimum_t,
            minimum_contact_stage_m,
            maximum_stage_m,
            noncontact,
        };
        let onset_stage = -branches.law.interaction_range_m;
        let onset = jump(
            MaugisSnapKind::AttractionOnset,
            branches.detached(onset_stage)?,
            branches.connected(onset_stage)?,
        )?;
        let separation = jump(
            MaugisSnapKind::Separation,
            branches.connected(minimum)?,
            branches.detached(minimum)?,
        )?;
        branches.contact(maximum_stage_m)?;
        Ok(Self {
            branches,
            identity: hash.finalize(),
            onset,
            separation,
        })
    }
    /// Original immutable material/interface, geometry and operating envelope.
    #[must_use]
    pub const fn contact(&self) -> &SourceMaugisSphereContact {
        &self.branches.law
    }
    /// Explicit positive loading-machine stiffness, N/m.
    #[must_use]
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.branches.stiffness
    }
    /// Original law plus finite-stiffness content identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Prepared snap-in event at the finite interaction-range cutoff.
    #[must_use]
    pub const fn onset(&self) -> &MaugisFixtureJump {
        &self.onset
    }
    /// Prepared snap-off event at this fixture's connected-branch fold.
    #[must_use]
    pub const fn separation(&self) -> &MaugisFixtureJump {
        &self.separation
    }
    /// Smallest stage approach on the stable solid-contact branch, m.
    #[must_use]
    pub const fn minimum_contact_stage_m(&self) -> f64 {
        self.branches.minimum_contact_stage_m
    }
    /// Largest stage approach within the original indentation envelope, m.
    #[must_use]
    pub const fn maximum_stage_m(&self) -> f64 {
        self.branches.maximum_stage_m
    }
    /// Stable attractive disk interval, if it exists. Neither endpoint is detached.
    #[must_use]
    pub fn noncontact_stage_interval_m(&self) -> Option<[f64; 2]> {
        self.branches
            .noncontact
            .map(|n| [n.minimum_stage_m, n.touch_stage_m])
    }
    /// Evaluate the stable solid-contact branch, including its branch endpoint.
    pub fn evaluate_contact(
        &self,
        stage_m: f64,
    ) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        self.branches.contact(stage_m)
    }
    /// Evaluate the stable deformed attractive disk; solid contact area stays zero.
    pub fn evaluate_noncontact(
        &self,
        stage_m: f64,
    ) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        self.branches.noncontact(stage_m)
    }
    /// Evaluate the zero-force, relaxed-spring branch beyond the interaction range.
    pub fn evaluate_detached(
        &self,
        stage_m: f64,
    ) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        self.branches.detached(stage_m)
    }
    /// Explicit initial physical phase; no prehistory or initial jump is inferred.
    pub fn initialize(
        &self,
        phase: MaugisPhase,
        stage_m: f64,
    ) -> Result<MaugisFixtureState, MaugisFixtureError> {
        let response = match phase {
            MaugisPhase::Detached => self.evaluate_detached(stage_m)?,
            MaugisPhase::Noncontact => self.evaluate_noncontact(stage_m)?,
            MaugisPhase::Contact => self.evaluate_contact(stage_m)?,
        };
        Ok(MaugisFixtureState {
            fixture_identity: self.identity,
            response,
        })
    }
    /// Follow one monotone stage segment. Split reversals into separate calls.
    /// Equality consumes an approached event; zero motion never re-releases it.
    pub fn advance_monotone(
        &self,
        previous: &MaugisFixtureState,
        stage_m: f64,
        budget: MaugisContinuationBudget,
    ) -> Result<MaugisFixtureTransition, MaugisFixtureError> {
        validate_budget(budget)?;
        if previous.fixture_identity != self.identity {
            return Err(MaugisFixtureError::HistoryMismatch);
        }
        self.branches.validate_stage(stage_m)?;
        let old = previous.response;
        let detached = old.equilibrium.phase() == MaugisPhase::Detached;
        let event =
            if stage_m > old.stage_m && detached && stage_m >= self.onset.stage_displacement_m() {
                Some(self.onset)
            } else if stage_m < old.stage_m
                && !detached
                && stage_m <= self.separation.stage_displacement_m()
            {
                Some(self.separation)
            } else {
                None
            };
        let (response, work) = if stage_m == old.stage_m {
            (old, zero_work())
        } else if let Some(event) = event {
            match event.kind {
                MaugisSnapKind::AttractionOnset => {
                    let end = self.branches.connected(stage_m)?;
                    (end, self.branches.work(event.after, end, budget.work)?)
                }
                MaugisSnapKind::Separation => (
                    self.branches.detached(stage_m)?,
                    self.branches.work(old, event.before, budget.work)?,
                ),
            }
        } else if detached {
            (self.branches.detached(stage_m)?, zero_work())
        } else {
            let end = self.branches.connected(stage_m)?;
            (end, self.branches.work(old, end, budget.work)?)
        };
        let released = event.map_or(0.0, |j| j.released_energy_j);
        let change = response.total_potential_j - old.total_potential_j;
        let residual = (work.work_j - change) - released;
        let scale = old
            .total_potential_j
            .abs()
            .max(response.total_potential_j.abs())
            .max(work.work_j.abs())
            .max(released);
        let tolerance =
            budget.energy_absolute_tolerance_j + budget.energy_relative_tolerance * scale;
        if [change, residual, tolerance].iter().any(|x| !x.is_finite())
            || residual.abs() > tolerance
        {
            return Err(MaugisFixtureError::EnergyClosure {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        Ok(MaugisFixtureTransition {
            next: MaugisFixtureState {
                fixture_identity: self.identity,
                response,
            },
            jump: event,
            work,
            potential_change_j: change,
            released_energy_j: released,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
        })
    }
}

impl Branches {
    fn validate_stage(&self, stage: f64) -> Result<(), MaugisFixtureError> {
        if !stage.is_finite() {
            return Err(MaugisError::InvalidInput {
                field: "fixture stage",
            }
            .into());
        }
        if stage > self.maximum_stage_m {
            return Err(MaugisError::OutsideApplicability {
                ratio: "fixture stage",
                value: stage,
                limit: self.maximum_stage_m,
            }
            .into());
        }
        Ok(())
    }
    fn contact(&self, stage: f64) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        self.validate_stage(stage)?;
        if stage < self.minimum_contact_stage_m {
            return Err(MaugisFixtureError::NoStableEquilibrium {
                phase: MaugisPhase::Contact,
                stage_m: stage,
            });
        }
        let value = |s: kernel::State| (s.d + self.compliance * s.p) * self.law.indentation_scale_m;
        let mut lo = self.minimum_t;
        let mut hi = self.law.maximum_t;
        if stage == self.minimum_contact_stage_m {
            hi = lo;
        }
        if stage == self.maximum_stage_m {
            lo = hi;
        }
        for _ in 0..128 {
            let mid = 0.5 * lo + 0.5 * hi;
            if mid == lo || mid == hi {
                break;
            }
            if value(kernel::state(self.law.lambda, mid)?) < stage {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let left = kernel::state(self.law.lambda, lo)?;
        let right = kernel::state(self.law.lambda, hi)?;
        let s = if (value(left) - stage).abs() < (value(right) - stage).abs() {
            left
        } else {
            right
        };
        let equilibrium = self.law.contact_energy_at_state(s)?;
        let tangent = if stage == self.minimum_contact_stage_m || s.t == 0.0 {
            None
        } else {
            let slope = s.ds + self.compliance * s.ps;
            if slope <= 0.0 {
                return Err(unrepresentable("stable fixture tangent"));
            }
            Some((self.law.force_scale_n / self.law.indentation_scale_m) * (s.ps / slope))
        };
        self.response(equilibrium, stage, tangent, Some(s.t))
    }
    fn noncontact(&self, stage: f64) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        self.validate_stage(stage)?;
        let n = self
            .noncontact
            .filter(|n| stage >= n.minimum_stage_m && stage <= n.touch_stage_m)
            .ok_or(MaugisFixtureError::NoStableEquilibrium {
                phase: MaugisPhase::Noncontact,
                stage_m: stage,
            })?;
        let zero = kernel::state(self.law.lambda, 0.0)?;
        let beta = (8.0 / (3.0 * PI)) * self.law.lambda;
        let (c, delta) = if stage == n.touch_stage_m {
            (zero.c, n.touch_indentation_m)
        } else {
            let offset = ((stage - n.minimum_stage_m) / self.law.indentation_scale_m) / n.quadratic;
            let c = n.fold_radius + (2.0 * offset).sqrt();
            let delta = -self.law.interaction_range_m
                - 0.5 * beta * beta * self.law.indentation_scale_m
                + 0.5 * (c - beta).powi(2) * self.law.indentation_scale_m;
            (c, delta)
        };
        let equilibrium =
            self.law
                .noncontact_energy_at_radius(c, delta, c == beta || c == zero.c)?;
        let tangent = if stage == n.minimum_stage_m || stage == n.touch_stage_m {
            None
        } else {
            let slope = n.quadratic * (c - n.fold_radius);
            if slope <= 0.0 {
                return Err(unrepresentable("stable noncontact fixture tangent"));
            }
            Some(
                -PI * self.law.lambda
                    * (c / slope)
                    * (self.law.force_scale_n / self.law.indentation_scale_m),
            )
        };
        self.response(equilibrium, stage, tangent, None)
    }
    fn detached(&self, stage: f64) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        self.validate_stage(stage)?;
        let equilibrium = self.law.evaluate_detached(stage)?;
        self.response(equilibrium, stage, equilibrium.tangent_n_per_m(), None)
    }
    fn connected(&self, stage: f64) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        if self.noncontact.is_some_and(|n| stage < n.touch_stage_m) {
            self.noncontact(stage)
        } else {
            self.contact(stage)
        }
    }
    fn response(
        &self,
        equilibrium: MaugisEquilibrium,
        stage: f64,
        tangent: Option<f64>,
        contact_coordinate: Option<f64>,
    ) -> Result<MaugisFixtureResponse, MaugisFixtureError> {
        let force = equilibrium.normal_force_n();
        let extension = force / self.stiffness;
        let spring_energy_j = (0.5 * force) * extension;
        let total_potential_j = equilibrium.potential_j() + spring_energy_j;
        let residual = force - self.stiffness * (stage - equilibrium.indentation_m());
        let tolerance = 256.0
            * f64::EPSILON
            * self
                .law
                .force_scale_n
                .max(force.abs())
                .max((self.stiffness * stage).abs())
                .max((self.stiffness * equilibrium.indentation_m()).abs());
        if [
            extension,
            spring_energy_j,
            total_potential_j,
            residual,
            tolerance,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || tangent.is_some_and(|x| !x.is_finite())
            || residual.abs() > tolerance
        {
            return Err(unrepresentable("finite-stiffness equilibrium and storage"));
        }
        Ok(MaugisFixtureResponse {
            equilibrium,
            contact_coordinate,
            stage_m: stage,
            spring_energy_j,
            total_potential_j,
            stage_tangent_n_per_m: tangent,
            equilibrium_residual_n: residual,
        })
    }
    fn work(
        &self,
        from: MaugisFixtureResponse,
        to: MaugisFixtureResponse,
        budget: MaugisWorkBudget,
    ) -> Result<MaugisWorkEstimate, MaugisFixtureError> {
        if from.stage_m == to.stage_m {
            return Ok(zero_work());
        }
        let (lo, hi, sign) = if from.stage_m < to.stage_m {
            (from, to, 1.0)
        } else {
            (to, from, -1.0)
        };
        let mut work = zero_work();
        let mut contact_start = lo;
        if let Some(n) = self.noncontact.filter(|n| lo.stage_m < n.touch_stage_m) {
            let end = if hi.stage_m <= n.touch_stage_m {
                hi
            } else {
                self.noncontact(n.touch_stage_m)?
            };
            let c0 = lo.equilibrium.cohesive_radius_m() / self.law.radius_scale_m;
            let c1 = end.equilibrium.cohesive_radius_m() / self.law.radius_scale_m;
            let x0 = c0 - n.fold_radius;
            let x1 = c1 - n.fold_radius;
            let dx = (2.0 * ((end.stage_m - lo.stage_m) / self.law.indentation_scale_m))
                / (n.quadratic * (x0 + x1));
            let i1 = 0.5 * (x0 + x1);
            let i2 = (x0 * x0 + x0 * x1 + x1 * x1) / 3.0;
            let i3 = 0.25 * (x0 + x1) * (x0 * x0 + x1 * x1);
            work.work_j = -(PI / 2.0)
                * self.law.lambda
                * n.quadratic
                * dx
                * (n.fold_radius * n.fold_radius * i1 + 2.0 * n.fold_radius * i2 + i3)
                * self.law.force_scale_n
                * self.law.indentation_scale_m;
            if !dx.is_finite() || dx <= 0.0 {
                return Err(unrepresentable("noncontact stage work increment"));
            }
            contact_start = end;
        }
        if hi.stage_m > contact_start.stage_m {
            let a = contact_start.equilibrium;
            let b = hi.equilibrium;
            let from_t = match contact_start.contact_coordinate {
                Some(t) => t,
                None if self
                    .noncontact
                    .is_some_and(|n| contact_start.stage_m == n.touch_stage_m) =>
                {
                    0.0
                }
                None => return Err(unrepresentable("contact work start coordinate")),
            };
            let to_t = hi
                .contact_coordinate
                .ok_or_else(|| unrepresentable("contact work end coordinate"))?;
            let contact = self
                .law
                .reversible_work_between_coordinates(from_t, to_t, budget)?;
            let spring = (0.5 * a.normal_force_n() + 0.5 * b.normal_force_n())
                * ((b.normal_force_n() - a.normal_force_n()) / self.stiffness);
            work.work_j += contact.work_j + spring;
            work.estimated_absolute_error_j = contact.estimated_absolute_error_j;
            work.evaluations = contact.evaluations;
        }
        work.work_j *= sign;
        if !work.work_j.is_finite() {
            return Err(unrepresentable("finite-stiffness path work"));
        }
        Ok(work)
    }
}

fn jump(
    kind: MaugisSnapKind,
    before: MaugisFixtureResponse,
    after: MaugisFixtureResponse,
) -> Result<MaugisFixtureJump, MaugisFixtureError> {
    let released = before.total_potential_j - after.total_potential_j;
    if !released.is_finite() || released <= 0.0 {
        return Err(MaugisContinuationError::NonpassiveSnap {
            released_energy_j: released,
        }
        .into());
    }
    Ok(MaugisFixtureJump {
        kind,
        before,
        after,
        released_energy_j: released,
    })
}
fn zero_work() -> MaugisWorkEstimate {
    MaugisWorkEstimate {
        work_j: 0.0,
        estimated_absolute_error_j: 0.0,
        evaluations: 0,
    }
}
fn unrepresentable(field: &'static str) -> MaugisFixtureError {
    MaugisError::Unrepresentable { field }.into()
}
fn validate_budget(b: MaugisContinuationBudget) -> Result<(), MaugisFixtureError> {
    if !b.work.absolute_tolerance_j.is_finite()
        || b.work.absolute_tolerance_j <= 0.0
        || !b.work.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&b.work.relative_tolerance)
        || b.work.maximum_evaluations < 5
        || !b.energy_absolute_tolerance_j.is_finite()
        || b.energy_absolute_tolerance_j < 0.0
        || !b.energy_relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&b.energy_relative_tolerance)
    {
        return Err(MaugisError::InvalidInput {
            field: "fixture work/energy budget",
        }
        .into());
    }
    Ok(())
}
