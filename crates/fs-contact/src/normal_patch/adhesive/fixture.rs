//! Force control and a compliant, massless loading fixture for JKR contact.
//!
//! The local contact law remains the original JKR law. An explicitly supplied
//! linear spring connects the contact to a prescribed stage position D, so
//! F = K(D - delta). Stable equilibria require K + dF/d(delta) > 0 in addition
//! to local contact-radius stability. Their fold lies between the force- and
//! displacement-controlled JKR separation points.
//!
//! Total stored energy is U_contact + F^2/(2K), and its derivative with respect
//! to D is F. At fixed-stage attachment/separation the energy decrease is
//! reported separately for allocation to an unresolved fast transient, waves
//! or heat. No loss is inserted on a smooth elastic segment. The linear spring
//! and quasistatic, massless fixture are explicit assumptions: dynamics,
//! nonlinear fixtures, servo dynamics and arbitrary frame prestress are not
//! inferred. See Deng and Kesari, JMPS 131 (2019), 404-423,
//! <https://doi.org/10.1016/j.jmps.2019.07.009>.

use fs_blake3::{ContentHash, DomainHasher};

use super::{JkrBranch, JkrError, JkrJumpKind, JkrResponse, SourceJkrSphereContact, root_offset};
use crate::normal_patch::PointNormalPort;

/// Loading-specific refusals; the original contact applicability errors survive.
#[derive(Clone, Debug, PartialEq)]
pub enum JkrLoadingError {
    /// Original geometry, source, or constitutive numerical envelope refused.
    Contact(JkrError),
    /// No stable force-controlled equilibrium exists below this tensile load.
    NoForceEquilibrium {
        /// Requested signed load, N.
        force_n: f64,
        /// Force-controlled pull-off load, N.
        minimum_force_n: f64,
    },
    /// No attached fixture equilibrium exists below its stage-displacement fold.
    NoFixtureEquilibrium {
        /// Requested stage displacement, m.
        stage_displacement_m: f64,
        /// Lowest stable attached stage displacement, m.
        minimum_stage_displacement_m: f64,
    },
    /// The contact would exceed the original owner's maximum indentation.
    StageOutsideEnvelope {
        /// Requested stage displacement, m.
        stage_displacement_m: f64,
        /// Stage position at the original contact's upper operating boundary, m.
        maximum_stage_displacement_m: f64,
    },
    /// A history token belongs to a different contact or spring stiffness.
    StateMismatch,
    /// Finite-f64 arithmetic cannot resolve the declared loading regime.
    Unrepresentable {
        /// Quantity that cannot be resolved without widening the physical domain.
        field: &'static str,
    },
}

impl From<JkrError> for JkrLoadingError {
    fn from(error: JkrError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for JkrLoadingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Contact(error) => write!(f, "JKR loading: {error}"),
            Self::NoForceEquilibrium {
                force_n,
                minimum_force_n,
            } => write!(
                f,
                "JKR load {force_n} is below force-controlled pull-off {minimum_force_n}"
            ),
            Self::NoFixtureEquilibrium {
                stage_displacement_m,
                minimum_stage_displacement_m,
            } => write!(
                f,
                "JKR fixture stage {stage_displacement_m} is below its fold {minimum_stage_displacement_m}"
            ),
            Self::StageOutsideEnvelope {
                stage_displacement_m,
                maximum_stage_displacement_m,
            } => write!(
                f,
                "JKR fixture stage {stage_displacement_m} exceeds {maximum_stage_displacement_m}"
            ),
            Self::StateMismatch => write!(
                f,
                "JKR fixture history belongs to another source or stiffness"
            ),
            Self::Unrepresentable { field } => {
                write!(f, "unrepresentable JKR loading quantity: {field}")
            }
        }
    }
}
impl std::error::Error for JkrLoadingError {}

impl SourceJkrSphereContact {
    /// Stable force-controlled JKR branch, including the marginal pull-off point.
    /// The closed-form radius solves F/F0 = y(y-6), y=(a/a_s)^(3/2), on y>=3.
    /// The returned tangent is the CONTACT dF/d(delta), not a force-servo law.
    /// At the marginal point it is exactly zero; no finite load compliance is
    /// implied there. This query alone neither nucleates nor advances history.
    pub fn evaluate_force_controlled(&self, force_n: f64) -> Result<JkrResponse, JkrLoadingError> {
        if !force_n.is_finite() {
            return Err(JkrError::InvalidInput {
                field: "force-controlled normal load",
            }
            .into());
        }
        let minimum = self.force_controlled_pull_off_n();
        if force_n < minimum {
            return Err(JkrLoadingError::NoForceEquilibrium {
                force_n,
                minimum_force_n: minimum,
            });
        }
        let maximum = self.evaluate_attached(self.envelope.maximum_indentation_m)?;
        if force_n > maximum.normal_force_n {
            return Err(JkrError::OutsideApplicability {
                ratio: "force-controlled normal load",
                value: force_n,
                limit: maximum.normal_force_n,
            }
            .into());
        }
        if force_n == maximum.normal_force_n {
            return Ok(maximum);
        }
        // Subtract the dimensional pull-off first: this retains adjacent loads
        // instead of losing their difference in 9 + F/F0.
        let distance = (force_n - minimum) / self.force_scale_n;
        if !distance.is_finite() {
            return Err(JkrLoadingError::Unrepresentable {
                field: "relative applied load",
            });
        }
        let x = (3.0 + distance.sqrt()).cbrt();
        let t = x - 1.0;
        let indentation = indentation_at_offset(self, t);
        self.validate_indentation(indentation)?;
        let mut response = self.response_at_offset(t, indentation)?;
        if force_n == minimum {
            response.tangent_n_per_m = Some(0.0);
        }
        Ok(response)
    }
}

/// Immutable accepted branch history for one source and fixture stiffness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrFixtureState {
    identity: ContentHash,
    branch: JkrBranch,
    stage_displacement_m: f64,
}
impl JkrFixtureState {
    /// Attached or detached branch.
    #[must_use]
    pub const fn branch(self) -> JkrBranch {
        self.branch
    }
    /// Last accepted prescribed stage displacement, m.
    #[must_use]
    pub const fn stage_displacement_m(self) -> f64 {
        self.stage_displacement_m
    }
}

/// Contact and spring response at a prescribed stage displacement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrFixtureResponse {
    state: JkrFixtureState,
    /// Original source-bound contact response at the actual tip indentation.
    pub contact: JkrResponse,
    /// Prescribed stage displacement D, m; distinct from tip indentation.
    pub stage_displacement_m: f64,
    /// Signed spring extension F/K, m; negative on the tensile branch.
    pub spring_extension_m: f64,
    /// Nonnegative stored elastic spring energy, J.
    pub spring_energy_j: f64,
    /// Contact plus spring potential, J, relative to a detached relaxed fixture.
    pub total_stored_energy_j: f64,
    /// dF/dD, N/m; absent at the fixture fold where it diverges.
    pub stage_tangent_n_per_m: Option<f64>,
    /// Floating-point force-balance residual F - K(D-delta), N.
    pub equilibrium_residual_n: f64,
}
impl JkrFixtureResponse {
    /// Successor token, independent of mutable observation fields.
    #[must_use]
    pub const fn state(&self) -> JkrFixtureState {
        self.state
    }
    /// Signed port for the STAGE coordinate after statically condensing the
    /// massless contact/spring. `normal` and `moment_arm_m` describe that stage
    /// coordinate. For a tip-coordinate port use `contact.point_port` instead.
    pub fn stage_port(
        &self,
        normal: [f64; 3],
        moment_arm_m: [f64; 3],
    ) -> Result<PointNormalPort, JkrLoadingError> {
        let tangent = self
            .stage_tangent_n_per_m
            .ok_or(JkrError::SingularTangent)?;
        let mut port = self.contact.point_port(normal, moment_arm_m)?;
        port.tangent_n_per_m = tangent;
        Ok(port)
    }
}

/// A jump localized at a fixed stage position, with finite tip displacement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrFixtureJump {
    /// Attachment or separation.
    pub kind: JkrJumpKind,
    /// Prescribed stage displacement at the event, m.
    pub stage_displacement_m: f64,
    /// Actual tip indentation immediately before the event, m.
    pub before_indentation_m: f64,
    /// Actual tip indentation immediately after the event, m.
    pub after_indentation_m: f64,
    /// Energy released from contact plus spring storage at fixed stage, J.
    pub released_energy_j: f64,
}

/// One monotone prescribed-stage segment, with physical jump loss separated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrFixtureTransition {
    /// Endpoint with the successor history token.
    pub response: JkrFixtureResponse,
    /// At most one branch jump lies on a monotone stage segment.
    pub jump: Option<JkrFixtureJump>,
    /// Energy leaving the quasistatic model in the branch jump only, J.
    pub irreversible_work_j: f64,
    /// Integral F dD; equals delta(total stored energy) plus jump loss, J.
    pub prescribed_work_j: f64,
}

/// A JKR contact connected to its prescribed stage by an explicit linear spring.
/// Construction is allocation-free beyond ownership of the supplied source.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrSpringFixture {
    contact: SourceJkrSphereContact,
    stiffness_n_per_m: f64,
    identity: ContentHash,
    critical_offset: f64,
    maximum_offset: f64,
    critical_stage_m: f64,
    maximum_stage_m: f64,
    stage_scale_m: f64,
    excess_coefficients: [f64; 5],
    critical_contact: JkrResponse,
    attachment_contact: JkrResponse,
    attachment_release_j: f64,
    detachment_release_j: f64,
}

impl JkrSpringFixture {
    /// Supply the machine stiffness K>0 [N/m]. The original contact envelope
    /// must include the finite tip motion during attachment at D=0. A stiffness
    /// whose fold cannot be distinguished from rigid control in f64 refuses;
    /// use the original displacement-control owner for that exact limit.
    pub fn new(
        contact: SourceJkrSphereContact,
        stiffness_n_per_m: f64,
    ) -> Result<Self, JkrLoadingError> {
        super::positive(stiffness_n_per_m, "loading fixture stiffness")?;
        let spring_scale = contact.force_scale_n / stiffness_n_per_m;
        let stage_scale_m = contact.displacement_scale_m.max(spring_scale);
        let d = contact.displacement_scale_m / stage_scale_m;
        let f = spring_scale / stage_scale_m;
        if !spring_scale.is_finite() || spring_scale <= 0.0 || d <= 0.0 || f <= 0.0 {
            return Err(JkrLoadingError::Unrepresentable {
                field: "relative spring compliance",
            });
        }
        let mut lo = 0.0;
        let mut hi = (3.0_f64.cbrt() - 1.0).next_up();
        // g'(x)=4(x^3-1)+6 rho x^2(x^3-3). The offset form retains x^3-1
        // near the displacement fold; normalized coefficients cannot overflow.
        for _ in 0..96 {
            let t = 0.5 * lo + 0.5 * hi;
            if t == lo || t == hi {
                break;
            }
            let p = t * (t * t + 3.0 * t + 3.0);
            let x = 1.0 + t;
            let derivative = 4.0 * d * p + 6.0 * f * x * x * (p - 2.0);
            if derivative >= 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let critical_offset = hi;
        let delta_critical = indentation_at_offset(&contact, critical_offset);
        if delta_critical <= contact.minimum_attached_indentation_m() {
            return Err(JkrLoadingError::Unrepresentable {
                field: "finite-stiffness fold near rigid control",
            });
        }
        contact.validate_indentation(delta_critical)?;
        let critical_contact = contact.response_at_offset(critical_offset, delta_critical)?;
        let critical_stage_m = delta_critical + critical_contact.normal_force_n / stiffness_n_per_m;
        let maximum = contact.evaluate_attached(contact.envelope.maximum_indentation_m)?;
        let maximum_stage_m = maximum.indentation_m + maximum.normal_force_n / stiffness_n_per_m;
        if !critical_stage_m.is_finite()
            || critical_stage_m >= 0.0
            || !maximum_stage_m.is_finite()
            || maximum_stage_m < 0.0
        {
            return Err(JkrLoadingError::Unrepresentable {
                field: "fixture operating interval including attachment",
            });
        }
        let maximum_offset = root_offset(
            (contact.envelope.maximum_indentation_m - contact.minimum_attached_indentation_m())
                / contact.displacement_scale_m,
        );
        let c = 1.0 + critical_offset;
        let c2 = c * c;
        let c3 = c2 * c;
        // g(c+s)-g(c) has no linear term at the fold. All higher coefficients
        // are positive there. This polynomial avoids subtracting nearly equal
        // stage coordinates in the radius solve immediately after pull-off.
        let excess_coefficients = [
            6.0 * d * c2 + f * (15.0 * c2 * c2 - 18.0 * c),
            4.0 * d * c + f * (20.0 * c3 - 6.0),
            d + 15.0 * f * c2,
            6.0 * f * c,
            f,
        ];
        if excess_coefficients
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || maximum_offset <= critical_offset
        {
            return Err(JkrLoadingError::Unrepresentable {
                field: "fixture fold curvature or radius interval",
            });
        }
        let detachment_release_j = total_energy(critical_contact, stiffness_n_per_m)?.1;
        if detachment_release_j <= 0.0 {
            return Err(JkrLoadingError::Unrepresentable {
                field: "fixture separation release",
            });
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-spring-fixture.v1");
        hash.update(contact.identity().as_bytes());
        hash.update(&stiffness_n_per_m.to_bits().to_le_bytes());
        let mut model = Self {
            contact,
            stiffness_n_per_m,
            identity: hash.finalize(),
            critical_offset,
            maximum_offset,
            critical_stage_m,
            maximum_stage_m,
            stage_scale_m,
            excess_coefficients,
            critical_contact,
            attachment_contact: critical_contact,
            attachment_release_j: 0.0,
            detachment_release_j,
        };
        let attached = model.evaluate_attached(0.0)?;
        if attached.total_stored_energy_j >= 0.0 {
            return Err(JkrLoadingError::Unrepresentable {
                field: "fixture attachment release",
            });
        }
        model.attachment_contact = attached.contact;
        model.attachment_release_j = -attached.total_stored_energy_j;
        Ok(model)
    }

    /// Original source-bound contact owner, with unchanged material receipts.
    #[must_use]
    pub const fn contact(&self) -> &SourceJkrSphereContact {
        &self.contact
    }
    /// Explicit loading-machine stiffness, N/m.
    #[must_use]
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_per_m
    }
    /// Source, geometry, envelope and spring identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Compliance-dependent detachment stage coordinate, m.
    #[must_use]
    pub const fn minimum_attached_stage_m(&self) -> f64 {
        self.critical_stage_m
    }
    /// Upper stage coordinate admitted by the original contact envelope, m.
    #[must_use]
    pub const fn maximum_stage_m(&self) -> f64 {
        self.maximum_stage_m
    }
    /// Original contact response at the fixture fold, before separation.
    #[must_use]
    pub const fn critical_contact(&self) -> &JkrResponse {
        &self.critical_contact
    }

    /// Initialize a relaxed detached fixture at D<=0. It attaches at D=0 on
    /// the next monotone advance; a positive detached indentation is forbidden.
    pub fn detached_state(
        &self,
        stage_displacement_m: f64,
    ) -> Result<JkrFixtureState, JkrLoadingError> {
        self.validate_stage(stage_displacement_m)?;
        if stage_displacement_m > 0.0 {
            return Err(JkrError::InvalidInput {
                field: "positive detached fixture stage",
            }
            .into());
        }
        Ok(self.state(JkrBranch::Detached, stage_displacement_m))
    }

    /// Stable attached contact/spring equilibrium at the prescribed stage.
    /// The fold is included as a limiting response with no smooth stage tangent.
    /// Evaluation uses bounded scalar arithmetic and allocates nothing.
    pub fn evaluate_attached(
        &self,
        stage_displacement_m: f64,
    ) -> Result<JkrFixtureResponse, JkrLoadingError> {
        self.validate_stage(stage_displacement_m)?;
        if stage_displacement_m < self.critical_stage_m {
            return Err(JkrLoadingError::NoFixtureEquilibrium {
                stage_displacement_m,
                minimum_stage_displacement_m: self.critical_stage_m,
            });
        }
        let marginal = stage_displacement_m == self.critical_stage_m;
        let (contact, offset) = if marginal {
            (self.critical_contact, self.critical_offset)
        } else if stage_displacement_m == self.maximum_stage_m {
            (
                self.contact
                    .evaluate_attached(self.contact.envelope.maximum_indentation_m)?,
                self.maximum_offset,
            )
        } else {
            let q = (stage_displacement_m - self.critical_stage_m) / self.stage_scale_m;
            if !q.is_finite() {
                return Err(JkrLoadingError::Unrepresentable {
                    field: "stage displacement relative to fold",
                });
            }
            let excess = root_excess(
                q,
                self.maximum_offset - self.critical_offset,
                self.excess_coefficients,
            );
            let t = self.critical_offset + excess;
            // A scalar bracket endpoint may coincide with the admitted maximum
            // through roundoff; the original maximum response owns that endpoint.
            if t >= self.maximum_offset {
                (
                    self.contact
                        .evaluate_attached(self.contact.envelope.maximum_indentation_m)?,
                    self.maximum_offset,
                )
            } else {
                let delta = indentation_at_offset(&self.contact, t);
                self.contact.validate_indentation(delta)?;
                (self.contact.response_at_offset(t, delta)?, t)
            }
        };
        let spring_extension_m = contact.normal_force_n / self.stiffness_n_per_m;
        let (spring_energy_j, total_stored_energy_j) =
            total_energy(contact, self.stiffness_n_per_m)?;
        let stage_tangent_n_per_m = if marginal {
            None
        } else {
            let s = offset - self.critical_offset;
            let derivative = excess_derivative(s, self.excess_coefficients);
            let x = 1.0 + offset;
            let x2 = x * x;
            let value = (self.contact.force_scale_n / self.stage_scale_m)
                * (6.0 * x2 * (x2 * x - 3.0) / derivative);
            if !value.is_finite() {
                return Err(JkrLoadingError::Unrepresentable {
                    field: "condensed stage tangent",
                });
            }
            Some(value)
        };
        let equilibrium_residual_n = contact.normal_force_n
            - self.stiffness_n_per_m * (stage_displacement_m - contact.indentation_m);
        if !spring_extension_m.is_finite() || !equilibrium_residual_n.is_finite() {
            return Err(JkrLoadingError::Unrepresentable {
                field: "fixture extension or force residual",
            });
        }
        Ok(JkrFixtureResponse {
            state: self.state(JkrBranch::Attached, stage_displacement_m),
            contact,
            stage_displacement_m,
            spring_extension_m,
            spring_energy_j,
            total_stored_energy_j,
            stage_tangent_n_per_m,
            equilibrium_residual_n,
        })
    }

    /// Follow one monotone stage-displacement segment; split every reversal into
    /// a separate call. Only accepting the returned token advances history.
    pub fn advance_monotone(
        &self,
        previous: JkrFixtureState,
        stage_displacement_m: f64,
    ) -> Result<JkrFixtureTransition, JkrLoadingError> {
        if previous.identity != self.identity {
            return Err(JkrLoadingError::StateMismatch);
        }
        self.validate_stage(stage_displacement_m)?;
        let previous_energy = if previous.branch == JkrBranch::Attached {
            self.evaluate_attached(previous.stage_displacement_m)?
                .total_stored_energy_j
        } else {
            0.0
        };
        let (branch, jump) = match previous.branch {
            JkrBranch::Detached if stage_displacement_m >= 0.0 => (
                JkrBranch::Attached,
                Some(JkrFixtureJump {
                    kind: JkrJumpKind::Attach,
                    stage_displacement_m: 0.0,
                    before_indentation_m: 0.0,
                    after_indentation_m: self.attachment_contact.indentation_m,
                    released_energy_j: self.attachment_release_j,
                }),
            ),
            JkrBranch::Attached if stage_displacement_m <= self.critical_stage_m => (
                JkrBranch::Detached,
                Some(JkrFixtureJump {
                    kind: JkrJumpKind::Detach,
                    stage_displacement_m: self.critical_stage_m,
                    before_indentation_m: self.critical_contact.indentation_m,
                    after_indentation_m: self.critical_stage_m,
                    released_energy_j: self.detachment_release_j,
                }),
            ),
            _ => (previous.branch, None),
        };
        let response = if branch == JkrBranch::Attached {
            self.evaluate_attached(stage_displacement_m)?
        } else {
            self.detached_response(stage_displacement_m)
        };
        let irreversible_work_j = jump.map_or(0.0, |j| j.released_energy_j);
        let prescribed_work_j =
            response.total_stored_energy_j - previous_energy + irreversible_work_j;
        if !prescribed_work_j.is_finite() {
            return Err(JkrLoadingError::Unrepresentable {
                field: "fixture segment work",
            });
        }
        Ok(JkrFixtureTransition {
            response,
            jump,
            irreversible_work_j,
            prescribed_work_j,
        })
    }

    fn state(&self, branch: JkrBranch, stage_displacement_m: f64) -> JkrFixtureState {
        JkrFixtureState {
            identity: self.identity,
            branch,
            stage_displacement_m,
        }
    }
    fn detached_response(&self, stage_displacement_m: f64) -> JkrFixtureResponse {
        JkrFixtureResponse {
            state: self.state(JkrBranch::Detached, stage_displacement_m),
            contact: self.contact.detached_response(stage_displacement_m),
            stage_displacement_m,
            spring_extension_m: 0.0,
            spring_energy_j: 0.0,
            total_stored_energy_j: 0.0,
            stage_tangent_n_per_m: Some(0.0),
            equilibrium_residual_n: 0.0,
        }
    }
    fn validate_stage(&self, stage_displacement_m: f64) -> Result<(), JkrLoadingError> {
        if !stage_displacement_m.is_finite() {
            return Err(JkrError::InvalidInput {
                field: "loading stage displacement",
            }
            .into());
        }
        if stage_displacement_m > self.maximum_stage_m {
            return Err(JkrLoadingError::StageOutsideEnvelope {
                stage_displacement_m,
                maximum_stage_displacement_m: self.maximum_stage_m,
            });
        }
        Ok(())
    }
}

fn indentation_at_offset(contact: &SourceJkrSphereContact, t: f64) -> f64 {
    contact.minimum_attached_indentation_m()
        + contact.displacement_scale_m * t * t * (t * t + 4.0 * t + 6.0)
}
fn total_energy(contact: JkrResponse, stiffness: f64) -> Result<(f64, f64), JkrLoadingError> {
    let spring = (0.5 * contact.normal_force_n) * (contact.normal_force_n / stiffness);
    let total = contact.reversible_energy_j + spring;
    if !spring.is_finite() || spring < 0.0 || !total.is_finite() {
        return Err(JkrLoadingError::Unrepresentable {
            field: "contact plus fixture energy",
        });
    }
    Ok((spring, total))
}
fn excess_value(s: f64, c: [f64; 5]) -> f64 {
    s * s * (c[0] + s * (c[1] + s * (c[2] + s * (c[3] + s * c[4]))))
}
fn excess_derivative(s: f64, c: [f64; 5]) -> f64 {
    s * (2.0 * c[0] + s * (3.0 * c[1] + s * (4.0 * c[2] + s * (5.0 * c[3] + 6.0 * s * c[4]))))
}
fn root_excess(q: f64, maximum: f64, c: [f64; 5]) -> f64 {
    if q == 0.0 {
        return 0.0;
    }
    let mut lo = 0.0;
    let mut hi = (q / c[0]).sqrt().min(maximum);
    let mut s = hi;
    for _ in 0..80 {
        let value = excess_value(s, c);
        let residual = value - q;
        if residual.abs() <= 8.0 * f64::EPSILON * q {
            return s;
        }
        if residual > 0.0 {
            hi = s;
        } else {
            lo = s;
        }
        let next = s - residual / excess_derivative(s, c);
        let candidate = if next.is_finite() && next > lo && next < hi {
            next
        } else {
            0.5 * lo + 0.5 * hi
        };
        if candidate == s || candidate == lo || candidate == hi {
            return candidate;
        }
        s = candidate;
    }
    0.5 * lo + 0.5 * hi
}
