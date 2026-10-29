//! Reversible shear loading through a spring, stopped at its first physical limit.
//!
//! At fixed normal approach delta, a scalar tangential stage D loads the first
//! contact tangent through a physical spring K>0. The ORIGINAL mixed-mode owner
//! supplies U(a,delta,u), T=k a u and the radius curvature h at fixed u, where
//! k=8G*. Equilibrium is D=u+T/K. Eliminating u from U+K(D-u)^2/2 gives radius
//! curvature h_K=h-(k u)^2/(k a+K). Its first zero is the apparatus-dependent
//! instability. The independent sourced static-capacity condition is checked
//! on the actual shrinking area and repulsive load, using the existing law.
//!
//! This follows the conservative, lambda=1 fixed-normal-approach specialization
//! of Papangelo, Cricri and Ciavarella, JMPS 144 (2020), 104099, equations 7-12
//! and 24-26, <https://doi.org/10.1016/j.jmps.2020.104099>. Their discussion
//! distinguishes the radius instability from the area-based gross-slip limit.
//!
//! This owner ends AT the first limit. It never carries singular full-stick
//! tractions into gross sliding, invents a partial-slip zone, selects a post-snap
//! landing point, or turns stored-energy release into heat. The boundary is a
//! handoff with actual contact/spring storage and untraversed stage travel.
//! On the admitted reversible branch, including unloading and sign reversal,
//! stage work is exactly the change in contact-plus-spring potential; the fixed
//! normal coordinate contributes no work although its reaction can change.

use super::{JkrMixedModeError, JkrMixedModeResponse, SourceJkrMixedMode};
use crate::interface_binding::adhesive_shear::{AdhesiveShearError, SourceAdhesiveShear};
use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::TriboError;

/// Explicit numerical work and prescribed-stage range budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrShearLoadingLimits {
    /// Each bracket/continuation solve uses at most this many iterations, 1..=256.
    /// The original mixed-mode owner retains its own bounded inner solves.
    pub maximum_iterations: usize,
    /// Every requested stage coordinate must lie in this symmetric range [m].
    /// The physical first boundary may be outside this requested operating range.
    pub maximum_absolute_stage_m: f64,
}

/// Which independent restriction first ends increasing full-stick loading.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrShearBoundaryKind {
    /// The source's static junction capacity is reached, before radius instability.
    /// This is a gross-slip admissibility limit, not a resolved microslip onset.
    GrossSlip,
    /// The contact/spring Hessian loses radius stability, before static capacity.
    /// The subsequent motion can be sliding or separation and is not prescribed.
    RadiusInstability,
    /// Both conditions coincide at the same represented equilibrium.
    Coincident,
}

/// Original source errors and explicit loading/budget refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum JkrShearLoadingError {
    Contact(JkrMixedModeError),
    Shear(AdhesiveShearError),
    Invalid {
        field: &'static str,
    },
    StateMismatch,
    StageOutsideEnvelope {
        stage_m: f64,
        maximum_absolute_stage_m: f64,
    },
    BeyondFirstBoundary {
        stage_m: f64,
        boundary_m: f64,
        kind: JkrShearBoundaryKind,
    },
    NoConvergence {
        operation: &'static str,
    },
}
impl From<JkrMixedModeError> for JkrShearLoadingError {
    fn from(error: JkrMixedModeError) -> Self {
        Self::Contact(error)
    }
}
impl From<AdhesiveShearError> for JkrShearLoadingError {
    fn from(error: AdhesiveShearError) -> Self {
        Self::Shear(error)
    }
}
impl From<TriboError> for JkrShearLoadingError {
    fn from(error: TriboError) -> Self {
        Self::Shear(error.into())
    }
}
impl core::fmt::Display for JkrShearLoadingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sourced JKR shear loading refusal: {self:?}")
    }
}
impl std::error::Error for JkrShearLoadingError {}

/// Immutable accepted stage checkpoint, bound to both original sources and K.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrShearFixtureState {
    identity: ContentHash,
    stage_m: f64,
}
impl JkrShearFixtureState {
    #[must_use]
    pub const fn stage_displacement_m(self) -> f64 {
        self.stage_m
    }
}

/// One actual full-stick contact plus its loading apparatus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrShearFixtureResponse {
    state: JkrShearFixtureState,
    pub contact: JkrMixedModeResponse,
    pub stage_displacement_m: f64,
    pub spring_extension_m: f64,
    pub spring_energy_j: f64,
    pub total_stored_energy_j: f64,
    /// Original source-law capacity on this actual area and repulsive load.
    pub static_capacity_n: f64,
    /// Capacity minus |T|; no clipping or replacement of the actual force.
    pub static_margin_n: f64,
    /// Radius curvature after relaxing the tangential tip against K.
    pub fixture_radius_curvature_n_per_m: f64,
    /// dT/dD along the retained branch. Absent at a radius-instability handoff.
    pub stage_tangent_n_per_m: Option<f64>,
    /// T-K(D-u), exposing the actual floating-point load-line residual.
    pub force_balance_residual_n: f64,
}
impl JkrShearFixtureResponse {
    #[must_use]
    pub const fn state(&self) -> JkrShearFixtureState {
        self.state
    }
}

/// The positive-loading boundary; the negative boundary is its reflected state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrShearBoundary {
    pub kind: JkrShearBoundaryKind,
    pub response: JkrShearFixtureResponse,
    /// Final numerical root bracket in actual contact shift [m]. The response
    /// uses the admissible endpoint; this is not an interval-arithmetic proof.
    pub contact_displacement_bracket_m: [f64; 2],
}

/// One admitted segment, possibly shortened at the first slip/peel boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrShearFixtureTransition {
    pub response: JkrShearFixtureResponse,
    pub boundary: Option<JkrShearBoundaryKind>,
    pub requested_stage_displacement_m: f64,
    /// Signed travel left for a separately selected slip/dynamic owner [m].
    pub untraversed_stage_displacement_m: f64,
    /// Integral T dD on this retained reversible branch, with no jump/heat term.
    pub prescribed_work_j: f64,
}

/// Finite tangential spring and constant normal approach on original sources.
/// Runtime continuation allocates nothing and retains no mutable trial state.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrShearSpringFixture {
    contact: SourceJkrMixedMode,
    shear: SourceAdhesiveShear,
    indentation_m: f64,
    stiffness_n_per_m: f64,
    limits: JkrShearLoadingLimits,
    identity: ContentHash,
    initial: JkrShearFixtureResponse,
    boundary: JkrShearBoundary,
}

impl JkrShearSpringFixture {
    /// Bind the exact original normal and shear source states. The original
    /// envelopes must cover the displacement-controlled branch used to locate
    /// the apparatus fold; a failure to resolve that physical boundary refuses.
    /// No source strength is inferred from W, E, yield, or a material name.
    pub fn new(
        contact: SourceJkrMixedMode,
        shear: SourceAdhesiveShear,
        indentation_m: f64,
        stiffness_n_per_m: f64,
        limits: JkrShearLoadingLimits,
    ) -> Result<Self, JkrShearLoadingError> {
        positive(stiffness_n_per_m, "tangential apparatus stiffness")?;
        positive(limits.maximum_absolute_stage_m, "maximum absolute stage")?;
        if limits.maximum_iterations == 0 || limits.maximum_iterations > 256 {
            return Err(JkrShearLoadingError::Invalid {
                field: "maximum_iterations",
            });
        }
        if shear.normal_state_identity() != contact.normal_source().normal_state().identity() {
            return Err(JkrShearLoadingError::StateMismatch);
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-shear-spring.v1");
        hash.update(contact.identity().as_bytes());
        hash.update(shear.resolved().card_identity().as_bytes());
        for value in [
            indentation_m,
            stiffness_n_per_m,
            limits.maximum_absolute_stage_m,
            shear.law().static_shear_strength_pa(),
            shear.law().kinetic_shear_strength_pa(),
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&(limits.maximum_iterations as u64).to_le_bytes());
        let identity = hash.finalize();
        let sample = |u| trial(&contact, &shear, indentation_m, stiffness_n_per_m, u);
        let start = sample(0.0)?;
        let initial = start.response(identity, stiffness_n_per_m, 0.0, false)?;
        let (kind, terminal, bracket) = if start.curvature <= 0.0 {
            (
                if start.capacity == 0.0 {
                    JkrShearBoundaryKind::Coincident
                } else {
                    JkrShearBoundaryKind::RadiusInstability
                },
                start,
                [0.0; 2],
            )
        } else if start.capacity == 0.0 {
            (JkrShearBoundaryKind::GrossSlip, start, [0.0; 2])
        } else {
            let rigid_limit = contact.critical_tangential_displacement_m(indentation_m)?;
            let end = sample(rigid_limit)?;
            let (peeling, peel_bracket) = last_admissible(
                start,
                end,
                limits.maximum_iterations,
                "apparatus radius instability",
                &sample,
                |r| r.curvature,
            )?;
            let margin = peeling.capacity - peeling.contact.tangential_force_n[0];
            if margin > 0.0 {
                (
                    JkrShearBoundaryKind::RadiusInstability,
                    peeling,
                    peel_bracket,
                )
            } else if margin == 0.0 {
                (JkrShearBoundaryKind::Coincident, peeling, peel_bracket)
            } else {
                // Capacity/(k a) decreases as area shrinks for the original
                // nonnegative junction-plus-pressure law. This analytic upper
                // bracket also resolves arbitrarily weak, nonzero junctions.
                let stiffness =
                    contact.tangential_stiffness_per_radius * start.contact.contact_radius_m;
                let upper = (start.capacity / stiffness).next_up().min(peeling.shift());
                let (slip, slip_bracket) = last_admissible(
                    start,
                    sample(upper)?,
                    limits.maximum_iterations,
                    "sourced gross-slip threshold",
                    &sample,
                    |r| r.capacity - r.contact.tangential_force_n[0],
                )?;
                (JkrShearBoundaryKind::GrossSlip, slip, slip_bracket)
            }
        };
        let response = terminal.response(
            identity,
            stiffness_n_per_m,
            terminal.stage_m,
            kind != JkrShearBoundaryKind::GrossSlip,
        )?;
        let boundary = JkrShearBoundary {
            kind,
            response,
            contact_displacement_bracket_m: bracket,
        };
        Ok(Self {
            contact,
            shear,
            indentation_m,
            stiffness_n_per_m,
            limits,
            identity,
            initial,
            boundary,
        })
    }
    #[must_use]
    pub const fn contact(&self) -> &SourceJkrMixedMode {
        &self.contact
    }
    #[must_use]
    pub const fn shear_source(&self) -> &SourceAdhesiveShear {
        &self.shear
    }
    #[must_use]
    pub const fn indentation_m(&self) -> f64 {
        self.indentation_m
    }
    #[must_use]
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_per_m
    }
    #[must_use]
    pub const fn limits(&self) -> JkrShearLoadingLimits {
        self.limits
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn initial_response(&self) -> JkrShearFixtureResponse {
        self.initial
    }
    #[must_use]
    pub const fn first_boundary(&self) -> JkrShearBoundary {
        self.boundary
    }

    /// Unique equilibrium on the connected admissible branch. A stage beyond
    /// the first source-capacity or apparatus-stability boundary never clamps.
    pub fn evaluate(&self, stage_m: f64) -> Result<JkrShearFixtureResponse, JkrShearLoadingError> {
        self.check_stage(stage_m)?;
        let magnitude = stage_m.abs();
        let maximum = self.boundary.response.stage_displacement_m;
        if magnitude > maximum {
            return Err(JkrShearLoadingError::BeyondFirstBoundary {
                stage_m,
                boundary_m: stage_m.signum() * maximum,
                kind: self.boundary.kind,
            });
        }
        if magnitude == 0.0 {
            return Ok(self.initial);
        }
        if magnitude == maximum {
            let mut response = self.boundary.response;
            if stage_m < 0.0 {
                reflect(&mut response);
            }
            return Ok(response);
        }
        let mut lo = 0.0;
        let mut hi = self.boundary.response.contact.tangential_displacement_m[0];
        let initial_stiffness =
            self.contact.tangential_stiffness_per_radius * self.initial.contact.contact_radius_m;
        let mut u = magnitude / (1.0 + initial_stiffness / self.stiffness_n_per_m);
        for _ in 0..self.limits.maximum_iterations {
            let current = trial(
                &self.contact,
                &self.shear,
                self.indentation_m,
                self.stiffness_n_per_m,
                u,
            )?;
            let residual = current.stage_m - magnitude;
            if residual.abs() <= 8.0 * f64::EPSILON * magnitude {
                let mut response =
                    current.response(self.identity, self.stiffness_n_per_m, magnitude, false)?;
                if stage_m < 0.0 {
                    reflect(&mut response);
                }
                return Ok(response);
            }
            if residual > 0.0 {
                hi = u;
            } else {
                lo = u;
            }
            let derivative = (1.0
                + current.contact.contact_radius_m * self.contact.tangential_stiffness_per_radius
                    / self.stiffness_n_per_m)
                * current.curvature
                / current.contact.displacement_radius_curvature_n_per_m;
            let newton = u - residual / derivative;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                f64::midpoint(lo, hi)
            };
            if next == u {
                break;
            }
            u = next;
        }
        Err(JkrShearLoadingError::NoConvergence {
            operation: "stage equilibrium",
        })
    }

    /// Follow a monotone stage segment, including elastic unloading and sign
    /// reversal. If the target crosses the first boundary, return that boundary
    /// and its unused travel. No invalid post-slip endpoint is manufactured.
    /// Errors leave the caller's accepted checkpoint unchanged.
    pub fn advance_until_event(
        &self,
        state: JkrShearFixtureState,
        target_stage_m: f64,
    ) -> Result<JkrShearFixtureTransition, JkrShearLoadingError> {
        if state.identity != self.identity {
            return Err(JkrShearLoadingError::StateMismatch);
        }
        self.check_stage(target_stage_m)?;
        let before = self.evaluate(state.stage_m)?;
        let maximum = self.boundary.response.stage_displacement_m;
        let hits = target_stage_m.abs() >= maximum && target_stage_m != state.stage_m;
        let accepted = if target_stage_m.abs() > maximum {
            target_stage_m.signum() * maximum
        } else {
            target_stage_m
        };
        let response = self.evaluate(accepted)?;
        let prescribed_work_j = response.total_stored_energy_j - before.total_stored_energy_j;
        finite(prescribed_work_j, "prescribed stage work")?;
        Ok(JkrShearFixtureTransition {
            response,
            boundary: if hits { Some(self.boundary.kind) } else { None },
            requested_stage_displacement_m: target_stage_m,
            untraversed_stage_displacement_m: target_stage_m - accepted,
            prescribed_work_j,
        })
    }
    fn check_stage(&self, stage_m: f64) -> Result<(), JkrShearLoadingError> {
        finite(stage_m, "stage displacement")?;
        if stage_m.abs() > self.limits.maximum_absolute_stage_m {
            return Err(JkrShearLoadingError::StageOutsideEnvelope {
                stage_m,
                maximum_absolute_stage_m: self.limits.maximum_absolute_stage_m,
            });
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Trial {
    contact: JkrMixedModeResponse,
    capacity: f64,
    stage_m: f64,
    curvature: f64,
    tangent: Option<f64>,
}
impl Trial {
    fn shift(self) -> f64 {
        self.contact.tangential_displacement_m[0]
    }
    fn response(
        self,
        identity: ContentHash,
        stiffness: f64,
        stage_m: f64,
        marginal: bool,
    ) -> Result<JkrShearFixtureResponse, JkrShearLoadingError> {
        let force = self.contact.tangential_force_n[0];
        let extension = force / stiffness;
        let spring_energy_j = (0.5 * force) * extension;
        let total = self.contact.reversible_energy_j + spring_energy_j;
        let residual = force - stiffness * (stage_m - self.shift());
        for value in [extension, spring_energy_j, total, residual] {
            finite(value, "shear fixture response")?;
        }
        Ok(JkrShearFixtureResponse {
            state: JkrShearFixtureState { identity, stage_m },
            contact: self.contact,
            stage_displacement_m: stage_m,
            spring_extension_m: extension,
            spring_energy_j,
            total_stored_energy_j: total,
            static_capacity_n: self.capacity,
            static_margin_n: self.capacity - force.abs(),
            fixture_radius_curvature_n_per_m: self.curvature,
            stage_tangent_n_per_m: if marginal { None } else { self.tangent },
            force_balance_residual_n: residual,
        })
    }
}

fn trial(
    contact: &SourceJkrMixedMode,
    shear: &SourceAdhesiveShear,
    delta: f64,
    spring: f64,
    shift: f64,
) -> Result<Trial, JkrShearLoadingError> {
    let response = contact.evaluate_displacement_controlled(delta, [shift, 0.0])?;
    let capacity = shear
        .law()
        .static_capacity_n(response.friction_contact_state()?)?;
    let k = contact.tangential_stiffness_per_radius;
    let stiffness = k * response.contact_radius_m;
    let combined = stiffness + spring;
    positive(combined, "combined tangential stiffness")?;
    let cross = k * shift;
    let correction = cross * (cross / combined);
    let h = response.displacement_radius_curvature_n_per_m;
    let curvature = h - correction;
    let stage_m = shift + response.tangential_force_n[0] / spring;
    let tangent = if curvature > 0.0 {
        Some(spring * ((stiffness / combined) * h - correction) / curvature)
    } else {
        None
    };
    for value in [capacity, stage_m, curvature] {
        finite(value, "shear fixture equilibrium")?;
    }
    if let Some(value) = tangent {
        finite(value, "shear fixture tangent")?;
    }
    Ok(Trial {
        contact: response,
        capacity,
        stage_m,
        curvature,
        tangent,
    })
}

/// Locate a monotone loss of admissibility, retaining its physical side. Root
/// width is reduced to 8*f64::EPSILON relative to its scale; budget failure refuses.
fn last_admissible(
    mut lo: Trial,
    mut hi: Trial,
    budget: usize,
    operation: &'static str,
    sample: &impl Fn(f64) -> Result<Trial, JkrShearLoadingError>,
    margin: impl Fn(Trial) -> f64,
) -> Result<(Trial, [f64; 2]), JkrShearLoadingError> {
    if margin(lo) < 0.0 || margin(hi) > 0.0 {
        return Err(JkrShearLoadingError::NoConvergence { operation });
    }
    if margin(hi) == 0.0 {
        return Ok((hi, [hi.shift(); 2]));
    }
    for _ in 0..budget {
        let left = lo.shift();
        let right = hi.shift();
        if right - left <= 8.0 * f64::EPSILON * right || left.next_up() >= right {
            return Ok((lo, [left, right]));
        }
        let mid = sample(f64::midpoint(left, right))?;
        let value = margin(mid);
        if value == 0.0 {
            return Ok((mid, [mid.shift(); 2]));
        }
        if value > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(JkrShearLoadingError::NoConvergence { operation })
}

fn reflect(response: &mut JkrShearFixtureResponse) {
    response.state.stage_m = -response.state.stage_m;
    response.stage_displacement_m = -response.stage_displacement_m;
    response.spring_extension_m = -response.spring_extension_m;
    response.force_balance_residual_n = -response.force_balance_residual_n;
    response.contact.tangential_displacement_m[0] = -response.contact.tangential_displacement_m[0];
    response.contact.tangential_force_n[0] = -response.contact.tangential_force_n[0];
    for matrix in [
        &mut response.contact.coupled_stiffness_n_per_m,
        &mut response.contact.coupled_compliance_m_per_n,
    ]
    .into_iter()
    .flatten()
    {
        for i in [0, 2] {
            matrix[i][1] = -matrix[i][1];
            matrix[1][i] = -matrix[1][i];
        }
    }
}
fn finite(value: f64, field: &'static str) -> Result<(), JkrShearLoadingError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(JkrShearLoadingError::Invalid { field })
    }
}
fn positive(value: f64, field: &'static str) -> Result<(), JkrShearLoadingError> {
    finite(value, field)?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(JkrShearLoadingError::Invalid { field })
    }
}
