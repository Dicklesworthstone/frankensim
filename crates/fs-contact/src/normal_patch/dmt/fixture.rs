//! A massless linear loading spring on the original finite-range DMT law.
//!
//! At stage approach D and tip approach delta, F(delta)=K(D-delta).
//! Stability requires K+F'(delta)>0. The attractive tail has F'=-S, where
//! S=2*pi*R*W/z0. For K>S all three geometric branches connect reversibly.
//! For K<S the tail is unstable: approach snaps in at D=-z0 and withdrawal
//! snaps off at D=-2*pi*R*W/K. The released contact-plus-spring energy belongs
//! to the fast transient, not automatically to heat. This follows directly
//! from the source's piecewise potential; no JKR jump rule is substituted.
//!
//! The importance of machine stiffness for adhesive stability is discussed by
//! Deng & Kesari (2019), DOI 10.1016/j.jmps.2019.07.009. The contact and finite
//! range approximation, assumptions and source authority remain those of DMT.

use super::{DmtError, DmtRegime, DmtResponse, SourceDmtSphereContact};
use crate::normal_patch::PointNormalPort;
use fs_blake3::{ContentHash, DomainHasher};

/// A fixture/source refusal; an unstable attractive gap is never accepted.
#[derive(Clone, Debug, PartialEq)]
pub enum DmtFixtureError {
    /// Original contact/source/envelope refusal.
    Contact(DmtError),
    /// K=S leaves a continuum of neutral tail equilibria at D=-z0.
    NeutralAttraction,
    /// Stage is below the contacting branch's endpoint.
    NoContactEquilibrium,
    /// No stable noncontact branch at the requested stage.
    NoNoncontactEquilibrium,
    /// History belongs to a different normal source or spring stiffness.
    StateMismatch,
    /// A derived quantity cannot be represented with finite f64 precision.
    Unrepresentable { field: &'static str },
}
impl From<DmtError> for DmtFixtureError {
    fn from(value: DmtError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for DmtFixtureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "DMT fixture refusal: {self:?}")
    }
}
impl std::error::Error for DmtFixtureError {}

/// Accepted branch history, bound to the original source and spring stiffness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtFixtureState {
    identity: ContentHash,
    regime: DmtRegime,
    stage_displacement_m: f64,
}
impl DmtFixtureState {
    /// Physical contact/gap regime of this accepted equilibrium.
    #[must_use]
    pub const fn regime(self) -> DmtRegime {
        self.regime
    }
    /// Last accepted stage approach, m.
    #[must_use]
    pub const fn stage_displacement_m(self) -> f64 {
        self.stage_displacement_m
    }
}

/// A stable contact/spring equilibrium, or its limiting corner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtFixtureResponse {
    state: DmtFixtureState,
    /// Original contact response, including original source identity and area.
    pub contact: DmtResponse,
    /// Prescribed stage approach D, m; it is not the actual tip indentation.
    pub stage_displacement_m: f64,
    /// Signed spring extension F/K, m.
    pub spring_extension_m: f64,
    /// Nonnegative elastic spring storage, J.
    pub spring_energy_j: f64,
    /// Contact plus spring potential relative to a separated relaxed fixture, J.
    pub total_stored_energy_j: f64,
    /// Statically condensed dF/dD; absent at the geometric branch corners.
    pub stage_tangent_n_per_m: Option<f64>,
    /// F-K(D-delta), N, retained without silently correcting either force.
    pub equilibrium_residual_n: f64,
}
impl DmtFixtureResponse {
    /// Immutable successor history, independent of mutable observation fields.
    #[must_use]
    pub const fn state(&self) -> DmtFixtureState {
        self.state
    }
    /// Action/reaction at the prescribed STAGE coordinate. Tip ports remain
    /// available on the original contact response and have a different tangent.
    pub fn stage_port(
        &self,
        normal: [f64; 3],
        arm_m: [f64; 3],
    ) -> Result<PointNormalPort, DmtFixtureError> {
        let tangent = self.stage_tangent_n_per_m.ok_or(DmtError::CornerTangent)?;
        let mut port = self.contact.point_port(normal, arm_m)?;
        port.tangent_n_per_m = tangent;
        Ok(port)
    }
}

/// Instability of a soft fixture, localized at fixed stage displacement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DmtFixtureJumpKind {
    /// Motion from the range cutoff into finite-area solid contact.
    SnapIn,
    /// Motion from zero-area touch to separation beyond the attraction range.
    SnapOff,
}
/// Finite tip motion and energy release at one instability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtFixtureJump {
    /// Attachment or separation instability.
    pub kind: DmtFixtureJumpKind,
    /// Stage displacement at the instability, m.
    pub stage_displacement_m: f64,
    /// Tip indentation before the jump, m.
    pub before_indentation_m: f64,
    /// Tip indentation after the jump, m.
    pub after_indentation_m: f64,
    /// Nonnegative release of total contact and spring potential, J.
    pub released_energy_j: f64,
}
/// One monotone stage segment, separating smooth work from snap release.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtFixtureTransition {
    /// Endpoint with the successor history token.
    pub response: DmtFixtureResponse,
    /// At most one instability lies on a monotone segment.
    pub jump: Option<DmtFixtureJump>,
    /// Jump release only, J; its eventual heat/wave partition is external.
    pub irreversible_work_j: f64,
    /// Integral F dD along the actual stable path, J.
    pub prescribed_work_j: f64,
}

/// A source-preserving DMT contact with explicit finite loading compliance.
/// Evaluation is allocation-free and the contact root has at most 64 iterations.
#[derive(Clone, Debug, PartialEq)]
pub struct DmtSpringFixture {
    contact: SourceDmtSphereContact,
    stiffness_n_per_m: f64,
    identity: ContentHash,
    touch_stage_m: f64,
    maximum_stage_m: f64,
    soft: bool,
    // Positive coefficients of the contact equation in x=delta/delta_max.
    solve_scale_m: f64,
    linear: f64,
    hertz: f64,
}
impl DmtSpringFixture {
    /// Bind K>0 [N/m]. Exact neutral tail stiffness refuses; soft fixtures must
    /// admit the snap-in endpoint inside the original indentation envelope.
    pub fn new(
        contact: SourceDmtSphereContact,
        stiffness_n_per_m: f64,
    ) -> Result<Self, DmtFixtureError> {
        super::positive(stiffness_n_per_m, "loading spring stiffness")?;
        let s = contact.tail_stiffness_n_per_m;
        if stiffness_n_per_m == s {
            return Err(DmtFixtureError::NeutralAttraction);
        }
        let touch_stage_m = -contact.pull_off_magnitude_n / stiffness_n_per_m;
        let maximum = contact.evaluate(contact.envelope.maximum_indentation_m)?;
        let maximum_stage_m = maximum.indentation_m + maximum.normal_force_n / stiffness_n_per_m;
        let elastic_extension = maximum.repulsive_normal_force_n / stiffness_n_per_m;
        let solve_scale_m = maximum.indentation_m.max(elastic_extension);
        let linear = maximum.indentation_m / solve_scale_m;
        let hertz = elastic_extension / solve_scale_m;
        if !touch_stage_m.is_finite()
            || !maximum_stage_m.is_finite()
            || !solve_scale_m.is_finite()
            || linear <= 0.0
            || hertz <= 0.0
            || maximum_stage_m <= touch_stage_m
            || touch_stage_m == -contact.interaction_range_m
        {
            return Err(DmtFixtureError::Unrepresentable {
                field: "finite fixture compliance and stage bounds",
            });
        }
        let soft = stiffness_n_per_m < s;
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.dmt-spring-fixture.v1");
        hash.update(contact.identity().as_bytes());
        hash.update(&stiffness_n_per_m.to_bits().to_le_bytes());
        let model = Self {
            contact,
            stiffness_n_per_m,
            identity: hash.finalize(),
            touch_stage_m,
            maximum_stage_m,
            soft,
            solve_scale_m,
            linear,
            hertz,
        };
        if soft {
            model.snap_in()?;
            model.snap_off()?;
        }
        model.evaluate_contact(maximum_stage_m)?;
        Ok(model)
    }
    /// Retained original contact, material sources, geometry and envelope.
    #[must_use]
    pub const fn contact(&self) -> &SourceDmtSphereContact {
        &self.contact
    }
    /// Explicit loading stiffness, N/m.
    #[must_use]
    pub const fn stiffness_n_per_m(&self) -> f64 {
        self.stiffness_n_per_m
    }
    /// True only when the attraction gradient exceeds the loading stiffness.
    #[must_use]
    pub const fn has_snap_hysteresis(&self) -> bool {
        self.soft
    }
    /// Stage at zero-area touch, -2*pi*R*W/K [m]. For a soft fixture this is
    /// the snap-off threshold; for a stiff fixture it connects to a stable gap.
    #[must_use]
    pub const fn touch_stage_m(&self) -> f64 {
        self.touch_stage_m
    }
    /// Maximum stage admitted by the original positive-indentation envelope.
    #[must_use]
    pub const fn maximum_stage_m(&self) -> f64 {
        self.maximum_stage_m
    }
    /// Initialize a separated, relaxed fixture outside the interaction range.
    pub fn detached_state(&self, stage_m: f64) -> Result<DmtFixtureState, DmtFixtureError> {
        if !stage_m.is_finite() || stage_m > -self.contact.interaction_range_m {
            return Err(DmtError::InvalidInput {
                field: "detached fixture stage outside attraction range",
            }
            .into());
        }
        Ok(self
            .response(self.contact.evaluate(stage_m)?, stage_m)?
            .state())
    }
    /// Stable contacting branch, including its zero-area limiting endpoint.
    /// Soft fixtures can also have a separated equilibrium at the same D;
    /// callers must select physical branch history rather than minimum energy.
    pub fn evaluate_contact(&self, stage_m: f64) -> Result<DmtFixtureResponse, DmtFixtureError> {
        self.validate_stage(stage_m)?;
        if stage_m < self.touch_stage_m {
            return Err(DmtFixtureError::NoContactEquilibrium);
        }
        let delta = if stage_m == self.touch_stage_m {
            0.0
        } else if stage_m == self.maximum_stage_m {
            self.contact.envelope.maximum_indentation_m
        } else {
            let q = (stage_m - self.touch_stage_m) / self.solve_scale_m;
            if !q.is_finite() || q <= 0.0 {
                return Err(DmtFixtureError::Unrepresentable {
                    field: "contact stage excess",
                });
            }
            let mut lo = 0.0;
            let mut hi = (q / self.linear)
                .min((q / self.hertz).cbrt().powi(2))
                .min(1.0);
            let mut x = hi;
            for _ in 0..64 {
                let root = x.sqrt();
                let residual = x * (self.linear + self.hertz * root) - q;
                if residual.abs() <= 8.0 * f64::EPSILON * q {
                    break;
                }
                if residual > 0.0 {
                    hi = x;
                } else {
                    lo = x;
                }
                let next = x - residual / (self.linear + 1.5 * self.hertz * root);
                let next = if next.is_finite() && next > lo && next < hi {
                    next
                } else {
                    0.5 * lo + 0.5 * hi
                };
                if next == x {
                    break;
                }
                x = next;
            }
            let delta = x * self.contact.envelope.maximum_indentation_m;
            if delta <= 0.0 {
                return Err(DmtFixtureError::Unrepresentable {
                    field: "positive fixture indentation",
                });
            }
            delta
        };
        self.response(self.contact.evaluate(delta)?, stage_m)
    }
    /// Stable separated or attractive-gap branch. The soft fixture's unstable
    /// attractive tail is refused, even where its algebraic equilibrium exists.
    pub fn evaluate_noncontact(&self, stage_m: f64) -> Result<DmtFixtureResponse, DmtFixtureError> {
        self.validate_stage(stage_m)?;
        let cutoff = -self.contact.interaction_range_m;
        let delta = if stage_m <= cutoff {
            stage_m
        } else if self.soft || stage_m > self.touch_stage_m {
            return Err(DmtFixtureError::NoNoncontactEquilibrium);
        } else if stage_m == self.touch_stage_m {
            0.0
        } else {
            // Offset from touch avoids subtracting K*D and the pull-off force.
            let delta = (stage_m - self.touch_stage_m)
                / (1.0 - self.contact.tail_stiffness_n_per_m / self.stiffness_n_per_m);
            if !delta.is_finite() || delta <= cutoff || delta >= 0.0 {
                return Err(DmtFixtureError::Unrepresentable {
                    field: "stable attractive-gap displacement",
                });
            }
            delta
        };
        self.response(self.contact.evaluate(delta)?, stage_m)
    }
    /// Advance one monotone stage segment, splitting reversals into calls.
    /// Accepting the immutable endpoint token commits history; a refusal never
    /// mutates it. Snap releases are counted once even after endpoint overshoot.
    pub fn advance_monotone(
        &self,
        previous: DmtFixtureState,
        stage_m: f64,
    ) -> Result<DmtFixtureTransition, DmtFixtureError> {
        if previous.identity != self.identity {
            return Err(DmtFixtureError::StateMismatch);
        }
        self.validate_stage(stage_m)?;
        let start = if previous.regime == DmtRegime::Contact {
            self.evaluate_contact(previous.stage_displacement_m)?
        } else {
            self.evaluate_noncontact(previous.stage_displacement_m)?
        };
        let (response, jump, work) = if self.soft
            && previous.regime == DmtRegime::Separated
            && stage_m >= -self.contact.interaction_range_m
        {
            let (event, after) = self.snap_in()?;
            let end = self.evaluate_contact(stage_m)?;
            let work = self.smooth_work(after, end)?;
            (end, Some(event), work)
        } else if self.soft
            && previous.regime == DmtRegime::Contact
            && stage_m <= self.touch_stage_m
        {
            let (event, before) = self.snap_off()?;
            let end = self.evaluate_noncontact(stage_m)?;
            let work = self.smooth_work(start, before)?;
            (end, Some(event), work)
        } else {
            let end = if self.soft {
                if previous.regime == DmtRegime::Contact {
                    self.evaluate_contact(stage_m)?
                } else {
                    self.evaluate_noncontact(stage_m)?
                }
            } else if stage_m >= self.touch_stage_m {
                self.evaluate_contact(stage_m)?
            } else {
                self.evaluate_noncontact(stage_m)?
            };
            (end, None, self.smooth_work(start, end)?)
        };
        Ok(DmtFixtureTransition {
            response,
            jump,
            irreversible_work_j: jump.map_or(0.0, |event| event.released_energy_j),
            prescribed_work_j: work,
        })
    }
    fn snap_in(&self) -> Result<(DmtFixtureJump, DmtFixtureResponse), DmtFixtureError> {
        let stage = -self.contact.interaction_range_m;
        let after = self.evaluate_contact(stage)?;
        // At the event H+K*delta=(S-K)*z. Eliminate the large cancelling
        // adhesive/spring terms before evaluating the released energy.
        let delta = after.contact.indentation_m;
        let z = self.contact.interaction_range_m;
        let released =
            0.5 * z * (z * (self.contact.tail_stiffness_n_per_m - self.stiffness_n_per_m))
                + 0.5 * (self.stiffness_n_per_m * delta) * delta
                + 0.6 * after.contact.repulsive_normal_force_n * delta;
        if !released.is_finite() || released <= 0.0 {
            return Err(DmtFixtureError::Unrepresentable {
                field: "positive snap-in release",
            });
        }
        Ok((
            DmtFixtureJump {
                kind: DmtFixtureJumpKind::SnapIn,
                stage_displacement_m: stage,
                before_indentation_m: stage,
                after_indentation_m: after.contact.indentation_m,
                released_energy_j: released,
            },
            after,
        ))
    }
    fn snap_off(&self) -> Result<(DmtFixtureJump, DmtFixtureResponse), DmtFixtureError> {
        let before = self.evaluate_contact(self.touch_stage_m)?;
        // Factored exact release, avoiding subtraction of spring and adhesion
        // potentials as K approaches the critical attraction stiffness.
        let released = 0.5
            * self.contact.pull_off_magnitude_n
            * self.contact.interaction_range_m
            * ((self.contact.tail_stiffness_n_per_m - self.stiffness_n_per_m)
                / self.stiffness_n_per_m);
        if !released.is_finite() || released <= 0.0 {
            return Err(DmtFixtureError::Unrepresentable {
                field: "positive snap-off release",
            });
        }
        Ok((
            DmtFixtureJump {
                kind: DmtFixtureJumpKind::SnapOff,
                stage_displacement_m: self.touch_stage_m,
                before_indentation_m: 0.0,
                after_indentation_m: self.touch_stage_m,
                released_energy_j: released,
            },
            before,
        ))
    }
    fn response(
        &self,
        contact: DmtResponse,
        stage: f64,
    ) -> Result<DmtFixtureResponse, DmtFixtureError> {
        let extension = contact.normal_force_n / self.stiffness_n_per_m;
        let spring_energy_j = (0.5 * contact.normal_force_n) * extension;
        let total_stored_energy_j = match contact.regime {
            DmtRegime::Separated => 0.0,
            DmtRegime::AttractiveGap => {
                // The two quadratic storages nearly cancel as K approaches S.
                0.5 * contact.normal_force_n
                    * (contact.indentation_m + self.contact.interaction_range_m)
                    * ((self.stiffness_n_per_m - self.contact.tail_stiffness_n_per_m)
                        / self.stiffness_n_per_m)
            }
            DmtRegime::Contact
                if contact.repulsive_normal_force_n <= 0.5 * self.contact.pull_off_magnitude_n =>
            {
                let pull = self.contact.pull_off_magnitude_n;
                let hertz = contact.repulsive_normal_force_n;
                let touch_energy = 0.5
                    * pull
                    * self.contact.interaction_range_m
                    * ((self.contact.tail_stiffness_n_per_m - self.stiffness_n_per_m)
                        / self.stiffness_n_per_m);
                touch_energy + (0.4 * hertz - pull) * contact.indentation_m
                    - (hertz / self.stiffness_n_per_m) * (pull - 0.5 * hertz)
            }
            DmtRegime::Contact => contact.reversible_energy_j + spring_energy_j,
        };
        let residual =
            contact.normal_force_n - self.stiffness_n_per_m * (stage - contact.indentation_m);
        let tangent = contact.tangent_n_per_m.map(|k| {
            if k == 0.0 {
                0.0
            } else if k > self.stiffness_n_per_m {
                self.stiffness_n_per_m / (1.0 + self.stiffness_n_per_m / k)
            } else {
                k / (1.0 + k / self.stiffness_n_per_m)
            }
        });
        if !extension.is_finite()
            || !spring_energy_j.is_finite()
            || !total_stored_energy_j.is_finite()
            || !residual.is_finite()
            || tangent.is_some_and(|k| !k.is_finite())
        {
            return Err(DmtFixtureError::Unrepresentable {
                field: "fixture response",
            });
        }
        Ok(DmtFixtureResponse {
            state: DmtFixtureState {
                identity: self.identity,
                regime: contact.regime,
                stage_displacement_m: stage,
            },
            contact,
            stage_displacement_m: stage,
            spring_extension_m: extension,
            spring_energy_j,
            total_stored_energy_j,
            stage_tangent_n_per_m: tangent,
            equilibrium_residual_n: residual,
        })
    }
    fn smooth_work(
        &self,
        a: DmtFixtureResponse,
        b: DmtFixtureResponse,
    ) -> Result<f64, DmtFixtureError> {
        let da = a.contact.indentation_m;
        let db = b.contact.indentation_m;
        if a.stage_displacement_m != b.stage_displacement_m && da == db {
            return Err(DmtFixtureError::Unrepresentable {
                field: "distinct stage endpoints collapse to one tip coordinate",
            });
        }
        let contact = self.contact.reversible_work_between(da, db)?;
        // For nearby contacting endpoints use the factored Hertz difference;
        // the constant adhesion cancels before any floating-point subtraction.
        let force_change = if da >= 0.0 && db >= 0.0 && da.max(db) > 0.0 {
            let hi = da.max(db);
            let r = (da.min(db) / hi).sqrt();
            self.contact.hertz_coefficient * hi.sqrt() * ((1.0 + r + r * r) / (1.0 + r)) * (db - da)
        } else if da >= -self.contact.interaction_range_m
            && da <= 0.0
            && db >= -self.contact.interaction_range_m
            && db <= 0.0
        {
            -self.contact.tail_stiffness_n_per_m * (db - da)
        } else {
            b.contact.normal_force_n - a.contact.normal_force_n
        };
        let spring = (0.5 * a.contact.normal_force_n + 0.5 * b.contact.normal_force_n)
            * (force_change / self.stiffness_n_per_m);
        let work = contact + spring;
        if !work.is_finite() {
            return Err(DmtFixtureError::Unrepresentable {
                field: "stage work",
            });
        }
        Ok(work)
    }
    fn validate_stage(&self, stage: f64) -> Result<(), DmtFixtureError> {
        if !stage.is_finite() {
            return Err(DmtError::InvalidInput {
                field: "fixture stage",
            }
            .into());
        }
        if stage > self.maximum_stage_m {
            return Err(DmtError::OutsideApplicability {
                ratio: "fixture stage",
                value: stage,
                limit: self.maximum_stage_m,
            }
            .into());
        }
        Ok(())
    }
}
