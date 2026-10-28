//! Johnson--Kendall--Roberts adhesion of an elastic sphere and half space.
//!
//! The input `adhesion-energy` is the WORK OF ADHESION W [J/m2], not the surface
//! energy of one body. With c = sqrt(2 pi W / E*), the attached equilibrium is
//! delta = a^2/R - c sqrt(a), F = 4 E* a^3/(3R) - 2 E* c a^(3/2).
//! Positive delta denotes approach and positive F denotes compression. A finite
//! area and tensile force persist at negative approach. The stored potential
//! includes surface energy, so it can be negative; dU/d(delta) = F.
//!
//! This is the short-range, large-Tabor-parameter, frictionless JKR limit of
//! Johnson, Kendall and Roberts (1971), DOI: 10.1098/rspa.1971.0141. It is an
//! explicitly selected model estimate, not a material admission. Roughness,
//! viscoelasticity, finite-range attraction before touch, rate-dependent work
//! of adhesion and mixed-mode peeling are outside this law. A separate
//! [`fixture`] owner adds loading-machine stiffness without changing this law.
//! The local stable branch is DISPLACEMENT controlled; its
//! final separation point differs from the familiar force-controlled pull-off.

/// Stable applied-force and finite-loading-stiffness contact.
pub mod fixture;
/// Conservative full-stick Griffith coupling of normal adhesion and shear.
pub mod mixed_mode;

use core::f64::consts::PI;

use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::{InputAuthority, TriboError, adhesive::AdhesiveContactState};

use crate::interface_binding::{ADHESION_ENERGY_PROPERTY, BoundNormalInterfaceState};

use super::PointNormalPort;

/// Explicit geometry and operating limits of the selected JKR approximation.
/// Ratios are caller-authored limits, not propagated error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrEnvelope {
    /// Largest admitted positive approach, m.
    pub maximum_indentation_m: f64,
    /// Geometric depth available for the half-space approximation, m.
    pub half_space_depth_m: f64,
    /// Smallest relevant layer thickness, m.
    pub layer_thickness_m: f64,
    /// Upper bound on a/R, strictly between zero and one.
    pub maximum_patch_to_radius: f64,
    /// Upper bound on a/depth, strictly between zero and one.
    pub maximum_patch_to_depth: f64,
    /// Upper bound on a/thickness, strictly between zero and one.
    pub maximum_patch_to_layer: f64,
    /// Upper bound on a/R + sqrt(2 pi W/(E* a)), away from the edge singularity.
    pub maximum_nominal_strain: f64,
    /// Upper bound on the repulsive Hertz peak pressure / bulk yield stress.
    /// This is not a bound on the singular tensile edge stress.
    pub maximum_repulsive_pressure_to_yield: f64,
    /// Upper bound on the plane-stress Irwin estimate E*W/(pi yield^2 a).
    /// This small-scale-yielding screen does not resolve the cohesive zone.
    pub maximum_process_zone_to_patch: f64,
    /// Required Tabor parameter (R W^2 / (E*^2 z0^3))^(1/3), at least 5.
    /// Five is a model-selection threshold, not a specified accuracy guarantee.
    pub minimum_tabor_parameter: f64,
}

/// The two distinct local contact branches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrBranch {
    /// No attached area, force or surface-energy storage.
    Detached,
    /// The stable (or limiting marginal) displacement-controlled branch.
    Attached,
}

/// Immutable branch history bound to one material state and geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrState {
    model_identity: ContentHash,
    branch: JkrBranch,
    indentation_m: f64,
}

impl JkrState {
    /// Current branch. Returning a state does not decide global solver stability.
    #[must_use]
    pub const fn branch(self) -> JkrBranch {
        self.branch
    }
    /// Last accepted signed approach, m.
    #[must_use]
    pub const fn indentation_m(self) -> f64 {
        self.indentation_m
    }
}

/// One allocation-free local response; the original source owner is retained
/// by `SourceJkrSphereContact` rather than cloned at every force evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrResponse {
    state: JkrState,
    normal_state_identity: ContentHash,
    /// Current attached or detached branch.
    pub branch: JkrBranch,
    /// Signed approach, m; this must not be clamped to geometric penetration.
    pub indentation_m: f64,
    /// Actual circular contact radius, m.
    pub contact_radius_m: f64,
    /// Actual attached area, m2.
    pub contact_area_m2: f64,
    /// Hertz repulsive contribution before subtracting attraction, N.
    /// It is not the integral of the positive part of the JKR pressure field.
    pub repulsive_normal_force_n: f64,
    /// Signed resultant force, N; negative is adhesive tension.
    pub normal_force_n: f64,
    /// dF/d(delta), N/m. Absent at the displacement spinodal where it diverges.
    pub tangent_n_per_m: Option<f64>,
    /// Elastic plus surface energy, J, relative to the detached state.
    pub reversible_energy_j: f64,
    /// Model estimate, or `SyntheticFixture` when its original inputs are synthetic.
    pub authority: InputAuthority,
}

impl JkrResponse {
    /// State token for a subsequent monotone segment, independent of mutable
    /// observation fields in this public response.
    #[must_use]
    pub const fn state(&self) -> JkrState {
        self.state
    }
    /// Identity of the original bulk and ordered-interface state.
    #[must_use]
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_state_identity
    }
    /// Feed area and repulsive load to adhesive junction shear without treating
    /// a tensile or zero net force as absence of contact.
    pub fn friction_contact_state(&self) -> Result<AdhesiveContactState, TriboError> {
        AdhesiveContactState::new(self.contact_area_m2, self.repulsive_normal_force_n)
    }
    /// Signed action/reaction port at a common origin. The unit normal points
    /// from the reacting body to the acting body; gap = -indentation.
    /// Jump energy is reported only by `JkrTransition`, never as smooth power.
    pub fn point_port(
        &self,
        normal: [f64; 3],
        moment_arm_m: [f64; 3],
    ) -> Result<PointNormalPort, JkrError> {
        let tangent = self.tangent_n_per_m.ok_or(JkrError::SingularTangent)?;
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !length.is_finite()
            || (length - 1.0).abs() > 1e-12
            || moment_arm_m.iter().any(|x| !x.is_finite())
        {
            return Err(JkrError::InvalidInput {
                field: "JKR port frame",
            });
        }
        let f = normal.map(|n| n * self.normal_force_n);
        let r = moment_arm_m;
        let m = [
            r[1] * f[2] - r[2] * f[1],
            r[2] * f[0] - r[0] * f[2],
            r[0] * f[1] - r[1] * f[0],
        ];
        if f.iter().chain(m.iter()).any(|x| !x.is_finite()) {
            return Err(JkrError::Unrepresentable {
                field: "JKR port wrench",
            });
        }
        Ok(PointNormalPort {
            action_force_n: f,
            action_moment_n_m: m,
            reaction_force_n: f.map(|x| -x),
            reaction_moment_n_m: m.map(|x| -x),
            residual_force_n: [0.0; 3],
            tangent_n_per_m: tangent,
            dissipated_power_w: 0.0,
            irreversible_work_j: 0.0,
        })
    }
}

/// A discontinuity localized at its exact prescribed-approach threshold.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrJumpKind {
    /// Contact nucleates at zero approach in the zero-range JKR limit.
    Attach,
    /// The attached local minimum disappears at the displacement spinodal.
    Detach,
}

/// Energy released by a branch jump at fixed displacement, J.
/// The consumer owns its allocation to waves, heat or a resolved fast transient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrJump {
    /// Attachment or separation event.
    pub kind: JkrJumpKind,
    /// Approach at the event, even when the requested endpoint overshoots it.
    pub indentation_m: f64,
    /// Nonnegative energy leaving the quasistatic contact storage, J.
    pub released_energy_j: f64,
}

/// Exact work bookkeeping for one monotone prescribed-displacement segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTransition {
    /// Endpoint, carrying the successor history token.
    pub response: JkrResponse,
    /// At most one event lies on a monotone segment.
    pub jump: Option<JkrJump>,
    /// Released jump energy only, J; zero on a smooth elastic segment.
    pub irreversible_work_j: f64,
    /// Integral F d(delta) along the selected path, J, including its branches.
    pub prescribed_work_j: f64,
}

/// Errors retain applicability and history failures instead of clipping them.
#[derive(Clone, Debug, PartialEq)]
pub enum JkrError {
    /// Nonfinite, nonpositive or otherwise incompatible explicit data.
    InvalidInput { field: &'static str },
    /// A source property depends on an unowned physical state or rate.
    NonstaticSource,
    /// A dimensionless approximation or geometry bound is exceeded.
    OutsideApplicability {
        ratio: &'static str,
        value: f64,
        limit: f64,
    },
    /// The given approach is below the last attached equilibrium.
    NoAttachedEquilibrium {
        indentation_m: f64,
        minimum_indentation_m: f64,
    },
    /// A contact state was made by another material/geometry owner.
    StateMismatch,
    /// The tangent is unbounded at detachment; an event solver must take over.
    SingularTangent,
    /// The finite-f64 numerical envelope cannot represent this derived value.
    Unrepresentable { field: &'static str },
}
impl core::fmt::Display for JkrError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { field } => write!(f, "invalid JKR input: {field}"),
            Self::NonstaticSource => write!(
                f,
                "JKR requires static elastic and adhesion sources at fixed external conditions"
            ),
            Self::OutsideApplicability {
                ratio,
                value,
                limit,
            } => write!(f, "JKR {ratio}={value} exceeds limit {limit}"),
            Self::NoAttachedEquilibrium {
                indentation_m,
                minimum_indentation_m,
            } => write!(
                f,
                "JKR approach {indentation_m} is below spinodal {minimum_indentation_m}"
            ),
            Self::StateMismatch => write!(f, "JKR state belongs to another material or geometry"),
            Self::SingularTangent => write!(
                f,
                "JKR spinodal requires a branch event, not a smooth tangent"
            ),
            Self::Unrepresentable { field } => {
                write!(f, "JKR derived value is unrepresentable: {field}")
            }
        }
    }
}
impl std::error::Error for JkrError {}

/// Frozen source-bound elastic adhesion. Construction screens the whole local
/// operating box; evaluation uses a bounded scalar solve with no allocations.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrSphereContact {
    normal: BoundNormalInterfaceState,
    radius_m: f64,
    interaction_range_m: f64,
    envelope: JkrEnvelope,
    identity: ContentHash,
    spinodal_radius_m: f64,
    displacement_scale_m: f64,
    force_scale_n: f64,
    energy_scale_j: f64,
    tabor_parameter: f64,
}

impl SourceJkrSphereContact {
    /// Select JKR explicitly from a resolved dry interface; never infer W from
    /// bulk material names. `interaction_range_m` is an explicitly authored or
    /// measured attraction range used for the Tabor validity screen.
    pub fn new(
        normal: BoundNormalInterfaceState,
        effective_radius_m: f64,
        interaction_range_m: f64,
        envelope: JkrEnvelope,
    ) -> Result<Self, JkrError> {
        let w = normal.adhesion_energy_j_per_m2();
        for (field, value) in [
            ("work of adhesion", w),
            ("effective radius", effective_radius_m),
            ("interaction range", interaction_range_m),
            ("maximum indentation", envelope.maximum_indentation_m),
            ("half-space depth", envelope.half_space_depth_m),
            ("layer thickness", envelope.layer_thickness_m),
        ] {
            positive(value, field)?;
        }
        for (field, value) in [
            ("maximum patch/radius", envelope.maximum_patch_to_radius),
            ("maximum patch/depth", envelope.maximum_patch_to_depth),
            ("maximum patch/layer", envelope.maximum_patch_to_layer),
            ("maximum nominal strain", envelope.maximum_nominal_strain),
            (
                "maximum repulsive pressure/yield",
                envelope.maximum_repulsive_pressure_to_yield,
            ),
            (
                "maximum process zone/patch",
                envelope.maximum_process_zone_to_patch,
            ),
        ] {
            positive(value, field)?;
            if value >= 1.0 {
                return Err(JkrError::InvalidInput { field });
            }
        }
        if !envelope.minimum_tabor_parameter.is_finite() || envelope.minimum_tabor_parameter < 5.0 {
            return Err(JkrError::InvalidInput {
                field: "minimum Tabor parameter must be at least 5",
            });
        }
        let elastic = normal.elastic();
        if elastic
            .surface_a_state()
            .resolved()
            .properties()
            .iter()
            .chain(elastic.surface_b_state().resolved().properties())
            .any(|p| !p.is_static_at_fixed_external_conditions())
            || !normal
                .interface_state()
                .property(ADHESION_ENERGY_PROPERTY)
                .is_some_and(|p| p.is_static_at_fixed_external_conditions())
        {
            return Err(JkrError::NonstaticSource);
        }
        let e = elastic.reduced_modulus_pa();
        // a_s^3 = pi W R^2/(8 E*). Factored cube roots avoid squaring R first.
        let spinodal_radius_m =
            (PI / 8.0).cbrt() * (w / e).cbrt() * effective_radius_m.cbrt().powi(2);
        let displacement_scale_m = spinodal_radius_m * (spinodal_radius_m / effective_radius_m);
        let force_scale_n = (PI / 6.0) * w * effective_radius_m;
        let energy_scale_j = force_scale_n * displacement_scale_m;
        let tabor_parameter =
            (displacement_scale_m / interaction_range_m) / (PI / 8.0).cbrt().powi(2);
        for (field, value) in [
            ("spinodal radius", spinodal_radius_m),
            ("displacement scale", displacement_scale_m),
            ("force scale", force_scale_n),
            ("energy scale", energy_scale_j),
            ("Tabor parameter", tabor_parameter),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(JkrError::Unrepresentable { field });
            }
        }
        if tabor_parameter < envelope.minimum_tabor_parameter {
            return Err(JkrError::OutsideApplicability {
                ratio: "required Tabor/actual Tabor",
                value: envelope.minimum_tabor_parameter / tabor_parameter,
                limit: 1.0,
            });
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-sphere.v1");
        hash.update(normal.identity().as_bytes());
        for x in [
            effective_radius_m,
            interaction_range_m,
            envelope.maximum_indentation_m,
            envelope.half_space_depth_m,
            envelope.layer_thickness_m,
            envelope.maximum_patch_to_radius,
            envelope.maximum_patch_to_depth,
            envelope.maximum_patch_to_layer,
            envelope.maximum_nominal_strain,
            envelope.maximum_repulsive_pressure_to_yield,
            envelope.maximum_process_zone_to_patch,
            envelope.minimum_tabor_parameter,
        ] {
            hash.update(&x.to_bits().to_le_bytes());
        }
        let model = Self {
            normal,
            radius_m: effective_radius_m,
            interaction_range_m,
            envelope,
            identity: hash.finalize(),
            spinodal_radius_m,
            displacement_scale_m,
            force_scale_n,
            energy_scale_j,
            tabor_parameter,
        };
        // Radius, pressure and process-zone ratios are monotone; the nominal
        // strain has only an interior minimum. These two corners screen all a.
        for delta in [
            model.minimum_attached_indentation_m(),
            envelope.maximum_indentation_m,
        ] {
            model.check_applicability(model.evaluate_attached(delta)?.contact_radius_m)?;
        }
        Ok(model)
    }

    /// Original immutable bulk and interface receipts and declared authority.
    #[must_use]
    pub const fn normal_state(&self) -> &BoundNormalInterfaceState {
        &self.normal
    }
    /// Source, geometry and operating-box identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Explicit operating box.
    #[must_use]
    pub const fn envelope(&self) -> JkrEnvelope {
        self.envelope
    }
    /// Relative sphere/plane radius, m.
    #[must_use]
    pub const fn effective_radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Frozen attraction range used only for model applicability, m.
    #[must_use]
    pub const fn interaction_range_m(&self) -> f64 {
        self.interaction_range_m
    }
    /// Dimensionless Tabor parameter, with its convention stated in the envelope.
    #[must_use]
    pub const fn tabor_parameter(&self) -> f64 {
        self.tabor_parameter
    }
    /// Displacement-controlled spinodal, delta_s = -3 a_s^2/R, m.
    #[must_use]
    pub fn minimum_attached_indentation_m(&self) -> f64 {
        -3.0 * self.displacement_scale_m
    }
    /// Force-controlled pull-off is -3 pi W R/2, N. This is not the event
    /// threshold of the displacement-controlled state machine.
    #[must_use]
    pub fn force_controlled_pull_off_n(&self) -> f64 {
        -9.0 * self.force_scale_n
    }

    /// Initialize a detached contact at nonpositive approach. Attachment at
    /// exactly zero occurs on the next `advance_monotone` call.
    pub fn detached_state(&self, indentation_m: f64) -> Result<JkrState, JkrError> {
        self.validate_indentation(indentation_m)?;
        if indentation_m > 0.0 {
            return Err(JkrError::InvalidInput {
                field: "positive detached indentation",
            });
        }
        Ok(self.state(JkrBranch::Detached, indentation_m))
    }

    /// Evaluate the stable attached branch, including its marginal endpoint.
    /// The scalar equation is solved in t = sqrt(a/a_s)-1 using
    /// (delta-delta_s)/h = t^2 (t^2+4t+6), which avoids subtractive cancellation
    /// close to pull-off. The safeguarded solve is bounded by 64 iterations.
    pub fn evaluate_attached(&self, indentation_m: f64) -> Result<JkrResponse, JkrError> {
        self.validate_indentation(indentation_m)?;
        let minimum = self.minimum_attached_indentation_m();
        if indentation_m < minimum {
            return Err(JkrError::NoAttachedEquilibrium {
                indentation_m,
                minimum_indentation_m: minimum,
            });
        }
        let q = (indentation_m - minimum) / self.displacement_scale_m;
        if !q.is_finite() {
            return Err(JkrError::Unrepresentable {
                field: "dimensionless approach",
            });
        }
        self.response_at_offset(root_offset(q), indentation_m)
    }

    // Shared radius-coordinate evaluation for displacement, force and fixture
    // control. Callers admit the corresponding indentation before entering.
    fn response_at_offset(&self, t: f64, indentation_m: f64) -> Result<JkrResponse, JkrError> {
        let x = 1.0 + t;
        let x2 = x * x;
        let x3 = x2 * x;
        let a = self.spinodal_radius_m * x2;
        let area = PI * a * a;
        let repulsive = self.force_scale_n * x3 * x3;
        let force = self.force_scale_n * x3 * (x3 - 6.0);
        let energy = self.energy_scale_j * x2 * x2 * ((0.4 * x3 - 4.0) * x3 + 6.0);
        let tangent = if t == 0.0 {
            None
        } else {
            // x^3-1 = t(t^2+3t+3), evaluated without cancellation.
            Some(
                (1.5 * self.force_scale_n / self.displacement_scale_m)
                    * x2
                    * ((x3 - 3.0) / (t * (t * t + 3.0 * t + 3.0))),
            )
        };
        for (field, value) in [
            ("contact radius", a),
            ("contact area", area),
            ("repulsive force", repulsive),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(JkrError::Unrepresentable { field });
            }
        }
        if !force.is_finite() || !energy.is_finite() || tangent.is_some_and(|k| !k.is_finite()) {
            return Err(JkrError::Unrepresentable {
                field: "attached force, energy or tangent",
            });
        }
        Ok(JkrResponse {
            state: self.state(JkrBranch::Attached, indentation_m),
            normal_state_identity: self.normal.identity(),
            branch: JkrBranch::Attached,
            indentation_m,
            contact_radius_m: a,
            contact_area_m2: area,
            repulsive_normal_force_n: repulsive,
            normal_force_n: force,
            tangent_n_per_m: tangent,
            reversible_energy_j: energy,
            authority: self.authority(),
        })
    }

    /// Follow ONE monotone prescribed-approach segment. Split any reversal into
    /// separate calls. The returned token commits history only if the consumer
    /// accepts it; failed evaluations never mutate the previous state.
    pub fn advance_monotone(
        &self,
        previous: JkrState,
        indentation_m: f64,
    ) -> Result<JkrTransition, JkrError> {
        if previous.model_identity != self.identity {
            return Err(JkrError::StateMismatch);
        }
        self.validate_indentation(indentation_m)?;
        let start_energy = match previous.branch {
            JkrBranch::Detached => 0.0,
            JkrBranch::Attached => {
                self.evaluate_attached(previous.indentation_m)?
                    .reversible_energy_j
            }
        };
        let (branch, jump) = match previous.branch {
            JkrBranch::Detached if indentation_m >= 0.0 => {
                let released = -self.evaluate_attached(0.0)?.reversible_energy_j;
                (
                    JkrBranch::Attached,
                    Some(JkrJump {
                        kind: JkrJumpKind::Attach,
                        indentation_m: 0.0,
                        released_energy_j: released,
                    }),
                )
            }
            JkrBranch::Attached if indentation_m <= self.minimum_attached_indentation_m() => (
                JkrBranch::Detached,
                Some(JkrJump {
                    kind: JkrJumpKind::Detach,
                    indentation_m: self.minimum_attached_indentation_m(),
                    released_energy_j: 2.4 * self.energy_scale_j,
                }),
            ),
            _ => (previous.branch, None),
        };
        let response = match branch {
            JkrBranch::Attached => self.evaluate_attached(indentation_m)?,
            JkrBranch::Detached => self.detached_response(indentation_m),
        };
        let irreversible_work_j = jump.map_or(0.0, |event| event.released_energy_j);
        let prescribed_work_j = response.reversible_energy_j - start_energy + irreversible_work_j;
        if !irreversible_work_j.is_finite()
            || irreversible_work_j < 0.0
            || !prescribed_work_j.is_finite()
        {
            return Err(JkrError::Unrepresentable {
                field: "JKR transition work",
            });
        }
        Ok(JkrTransition {
            response,
            jump,
            irreversible_work_j,
            prescribed_work_j,
        })
    }

    fn state(&self, branch: JkrBranch, indentation_m: f64) -> JkrState {
        JkrState {
            model_identity: self.identity,
            branch,
            indentation_m,
        }
    }

    fn detached_response(&self, indentation_m: f64) -> JkrResponse {
        JkrResponse {
            state: self.state(JkrBranch::Detached, indentation_m),
            normal_state_identity: self.normal.identity(),
            branch: JkrBranch::Detached,
            indentation_m,
            contact_radius_m: 0.0,
            contact_area_m2: 0.0,
            repulsive_normal_force_n: 0.0,
            normal_force_n: 0.0,
            tangent_n_per_m: Some(0.0),
            reversible_energy_j: 0.0,
            authority: self.authority(),
        }
    }

    fn validate_indentation(&self, delta: f64) -> Result<(), JkrError> {
        if !delta.is_finite() {
            return Err(JkrError::InvalidInput {
                field: "signed indentation",
            });
        }
        if delta > self.envelope.maximum_indentation_m {
            return Err(JkrError::OutsideApplicability {
                ratio: "indentation",
                value: delta,
                limit: self.envelope.maximum_indentation_m,
            });
        }
        Ok(())
    }

    fn authority(&self) -> InputAuthority {
        match self
            .normal
            .elastic()
            .interface()
            .interface()
            .provenance()
            .authority()
        {
            InputAuthority::SyntheticFixture => InputAuthority::SyntheticFixture,
            InputAuthority::CallerDeclared | InputAuthority::Estimated => InputAuthority::Estimated,
        }
    }

    fn check_applicability(&self, a: f64) -> Result<(), JkrError> {
        let b = self.envelope;
        let elastic = self.normal.elastic();
        let e = elastic.reduced_modulus_pa();
        let yield_stress = elastic.limiting_yield_stress_pa();
        let w = self.normal.adhesion_energy_j_per_m2();
        for (ratio, value, limit) in [
            ("patch/radius", a / self.radius_m, b.maximum_patch_to_radius),
            (
                "patch/depth",
                a / b.half_space_depth_m,
                b.maximum_patch_to_depth,
            ),
            (
                "patch/layer",
                a / b.layer_thickness_m,
                b.maximum_patch_to_layer,
            ),
            (
                "nominal strain",
                a / self.radius_m + (2.0 * PI * (w / e) / a).sqrt(),
                b.maximum_nominal_strain,
            ),
            (
                "repulsive pressure/yield",
                (2.0 / PI) * (e / yield_stress) * (a / self.radius_m),
                b.maximum_repulsive_pressure_to_yield,
            ),
            (
                "process zone/patch",
                (e / yield_stress) * (w / yield_stress) / (PI * a),
                b.maximum_process_zone_to_patch,
            ),
        ] {
            if !value.is_finite() || value > limit {
                return Err(JkrError::OutsideApplicability {
                    ratio,
                    value,
                    limit,
                });
            }
        }
        Ok(())
    }
}

fn positive(value: f64, field: &'static str) -> Result<(), JkrError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(JkrError::InvalidInput { field })
    }
}

fn root_offset(q: f64) -> f64 {
    if q == 0.0 {
        return 0.0;
    }
    let mut lower = 0.0;
    let mut upper = (q.sqrt() / 6.0_f64.sqrt()).min(q.sqrt().sqrt()) * (1.0 + 8.0 * f64::EPSILON);
    let mut t = upper;
    for _ in 0..64 {
        let value = t * t * (t * t + 4.0 * t + 6.0);
        let residual = value - q;
        if residual.abs() <= 8.0 * f64::EPSILON * q {
            return t;
        }
        if residual > 0.0 {
            upper = t;
        } else {
            lower = t;
        }
        let next = t - residual / (4.0 * t * (t * t + 3.0 * t + 3.0));
        let candidate = if next.is_finite() && next > lower && next < upper {
            next
        } else {
            0.5 * lower + 0.5 * upper
        };
        if candidate == t || candidate == lower || candidate == upper {
            return candidate;
        }
        t = candidate;
    }
    0.5 * lower + 0.5 * upper
}
