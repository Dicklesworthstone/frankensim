//! Finite-zone Maugis--Dugdale adhesion for elastic sphere/plane contact.
//!
//! Solve the original cohesive-zone closure and load/indentation equations:
//! Maugis (1992), DOI 10.1016/0021-9797(92)90285-T; equations (2)--(8) in
//! Leroy & Willert (2021), DOI 10.1002/zamm.201900313. This is an explicit
//! finite traction sigma0=W/z0 model, not an interpolation between JKR and DMT.
//!
//! Only contacting equilibria and reversible work along their stable branch
//! are owned here. Zero-area contact can still attract. Neither that endpoint
//! nor a fold is a detached state: noncontact elastic deformation, nucleation,
//! snap energy and a detached-reference potential require a separate owner.
use super::PointNormalPort;
use crate::interface_binding::{ADHESION_ENERGY_PROPERTY, BoundNormalInterfaceState};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::{InputAuthority, TriboError, adhesive::AdhesiveContactState};

#[path = "maugis/kernel.rs"]
mod kernel;
#[path = "maugis/work.rs"]
mod work;
pub use work::{MaugisWorkBudget, MaugisWorkEstimate};

/// Explicit small-strain, smooth elastic half-space operating domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisEnvelope {
    /// Greatest positive indentation, m.
    pub maximum_indentation_m: f64,
    /// Available half-space depth, m.
    pub half_space_depth_m: f64,
    /// Smallest load-bearing layer thickness, m.
    pub layer_thickness_m: f64,
    /// Greatest outer cohesive radius c / sphere radius.
    pub maximum_interaction_radius_to_radius: f64,
    /// Greatest c / half-space depth.
    pub maximum_interaction_radius_to_depth: f64,
    /// Greatest c / layer thickness.
    pub maximum_interaction_radius_to_layer: f64,
    /// Greatest contact radius a / sphere radius.
    pub maximum_nominal_strain: f64,
    /// Greatest Hertz peak pressure / limiting bulk yield stress.
    pub maximum_repulsive_pressure_to_yield: f64,
    /// Greatest sigma0 / limiting bulk yield stress.
    pub maximum_adhesive_stress_to_yield: f64,
    /// Greatest z0 / sphere radius; cannot exceed 0.01.
    pub maximum_interaction_range_to_radius: f64,
}

/// A contacting-branch limit. It is not a detached equilibrium.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaugisEndpoint {
    /// Real contact area tends to zero while cohesive-zone attraction remains.
    ZeroArea,
    /// The stable displacement-controlled contact branch ends at a fold.
    DisplacementFold,
}

/// Immutable local response, retaining the original material/interface identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisResponse {
    normal_state_identity: ContentHash,
    /// Signed indentation, m.
    pub indentation_m: f64,
    /// Radius of true solid contact, m.
    pub contact_radius_m: f64,
    /// Outer radius of the finite cohesive zone, m.
    pub cohesive_radius_m: f64,
    /// True contact area, m2; excludes the adhesive annulus.
    pub contact_area_m2: f64,
    /// Hertz repulsive component, N, for the existing junction-shear owner.
    pub repulsive_normal_force_n: f64,
    /// Signed resultant normal force, N; attraction is negative.
    pub normal_force_n: f64,
    /// Analytic dF/d(delta), N/m; absent at a contacting-branch endpoint.
    pub tangent_n_per_m: Option<f64>,
    /// Contacting-branch endpoint, if any.
    pub endpoint: Option<MaugisEndpoint>,
    /// Explicit source/model authority, never promoted from synthetic input.
    pub authority: InputAuthority,
    adhesive_stress_pa: f64,
    hertz_peak_pressure_pa: f64,
    annulus_root_m: f64,
}
impl MaugisResponse {
    /// Identity of the unchanged normal interface source state.
    #[must_use]
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_state_identity
    }
    /// Supply only real solid area and repulsive load to adhesive shear.
    pub fn friction_contact_state(&self) -> Result<AdhesiveContactState, TriboError> {
        AdhesiveContactState::new(self.contact_area_m2, self.repulsive_normal_force_n)
    }
    /// Axisymmetric signed pressure, Pa: contact, constant-traction annulus,
    /// then zero. Pressure at the outer cohesive edge is the inner value.
    pub fn pressure_pa(&self, radius_m: f64) -> Result<f64, MaugisError> {
        if !radius_m.is_finite() || radius_m < 0.0 {
            return Err(MaugisError::InvalidInput {
                field: "radial pressure coordinate",
            });
        }
        let a = self.contact_radius_m;
        if radius_m > self.cohesive_radius_m {
            return Ok(0.0);
        }
        if radius_m >= a {
            return Ok(-self.adhesive_stress_pa);
        }
        let root = ((a - radius_m) * (a + radius_m)).sqrt();
        Ok(self.hertz_peak_pressure_pa * (root / a)
            - (2.0 / PI) * self.adhesive_stress_pa * self.annulus_root_m.atan2(root))
    }
    /// Signed action/reaction wrench at one origin on a smooth contact branch.
    pub fn point_port(
        &self,
        normal: [f64; 3],
        arm_m: [f64; 3],
    ) -> Result<PointNormalPort, MaugisError> {
        let tangent = self.tangent_n_per_m.ok_or(MaugisError::EndpointTangent)?;
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !length.is_finite()
            || (length - 1.0).abs() > 1e-12
            || arm_m.iter().any(|x| !x.is_finite())
        {
            return Err(MaugisError::InvalidInput {
                field: "Maugis port frame",
            });
        }
        let f = normal.map(|x| x * self.normal_force_n);
        let r = arm_m;
        let m = [
            r[1] * f[2] - r[2] * f[1],
            r[2] * f[0] - r[0] * f[2],
            r[0] * f[1] - r[1] * f[0],
        ];
        if f.iter().chain(m.iter()).any(|v| !v.is_finite()) {
            return Err(MaugisError::Unrepresentable {
                field: "Maugis port wrench",
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

/// Explicit refusal, preserving the distinction between a branch limit and separation.
#[derive(Clone, Debug, PartialEq)]
pub enum MaugisError {
    /// Invalid explicit data or kinematics.
    InvalidInput { field: &'static str },
    /// A frozen property depends on an unowned evolving response coordinate.
    NonstaticSource,
    /// An authored physical applicability limit is exceeded.
    OutsideApplicability {
        ratio: &'static str,
        value: f64,
        limit: f64,
    },
    /// No stable contacting equilibrium at this prescribed indentation.
    NoContactEquilibrium { indentation_m: f64, minimum_m: f64 },
    /// A load is below the stable force-controlled contact branch.
    BelowPullOff { load_n: f64, pull_off_n: f64 },
    /// A smooth tangent was requested at a contacting-branch endpoint.
    EndpointTangent,
    /// Finite precision cannot represent this response.
    Unrepresentable { field: &'static str },
    /// The explicit work-quadrature budget was exhausted.
    WorkBudgetExceeded,
}
impl core::fmt::Display for MaugisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Maugis contact refusal: {self:?}")
    }
}
impl std::error::Error for MaugisError {}

/// Original source-bound Maugis model. Construction precomputes stable-branch
/// bounds; evaluation is allocation-free and uses at most 128 scalar iterations.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceMaugisSphereContact {
    normal: BoundNormalInterfaceState,
    radius_m: f64,
    interaction_range_m: f64,
    envelope: MaugisEnvelope,
    identity: ContentHash,
    radius_scale_m: f64,
    indentation_scale_m: f64,
    force_scale_n: f64,
    lambda: f64,
    minimum_t: f64,
    load_fold_t: f64,
    maximum_t: f64,
}
impl SourceMaugisSphereContact {
    /// Explicitly select finite-range cohesive traction; W comes from the
    /// original ordered interface and z0 is authored/measured, never inferred.
    pub fn new(
        normal: BoundNormalInterfaceState,
        effective_radius_m: f64,
        interaction_range_m: f64,
        envelope: MaugisEnvelope,
    ) -> Result<Self, MaugisError> {
        let w = normal.adhesion_energy_j_per_m2();
        for (field, v) in [
            ("work of adhesion", w),
            ("effective radius", effective_radius_m),
            ("interaction range", interaction_range_m),
            ("maximum indentation", envelope.maximum_indentation_m),
            ("half-space depth", envelope.half_space_depth_m),
            ("layer thickness", envelope.layer_thickness_m),
        ] {
            positive(v, field)?;
        }
        for (field, v) in [
            (
                "cohesive radius/radius",
                envelope.maximum_interaction_radius_to_radius,
            ),
            (
                "cohesive radius/depth",
                envelope.maximum_interaction_radius_to_depth,
            ),
            (
                "cohesive radius/layer",
                envelope.maximum_interaction_radius_to_layer,
            ),
            ("nominal strain", envelope.maximum_nominal_strain),
            (
                "repulsive pressure/yield",
                envelope.maximum_repulsive_pressure_to_yield,
            ),
            (
                "adhesive stress/yield",
                envelope.maximum_adhesive_stress_to_yield,
            ),
        ] {
            positive(v, field)?;
            if v >= 1.0 {
                return Err(MaugisError::InvalidInput { field });
            }
        }
        positive(
            envelope.maximum_interaction_range_to_radius,
            "maximum range/radius",
        )?;
        if envelope.maximum_interaction_range_to_radius > 0.01 {
            return Err(MaugisError::InvalidInput {
                field: "maximum range/radius",
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
            return Err(MaugisError::NonstaticSource);
        }
        let e = elastic.reduced_modulus_pa();
        let radius_scale_m = (PI * w / (4.0 * e / 3.0)).cbrt() * effective_radius_m.cbrt().powi(2);
        let indentation_scale_m = radius_scale_m * (radius_scale_m / effective_radius_m);
        let force_scale_n = PI * w * effective_radius_m;
        let lambda = (2.0 / PI) * (indentation_scale_m / interaction_range_m);
        for (field, v) in [
            ("radius scale", radius_scale_m),
            ("indentation scale", indentation_scale_m),
            ("force scale", force_scale_n),
            ("Maugis parameter", lambda),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(MaugisError::Unrepresentable { field });
            }
        }
        let minimum_t = kernel::fold(lambda, false)?;
        let load_fold_t = kernel::fold(lambda, true)?;
        let target = envelope.maximum_indentation_m / indentation_scale_m;
        if !target.is_finite() {
            return Err(MaugisError::Unrepresentable {
                field: "scaled indentation",
            });
        }
        let mut hi = 1.0_f64.max(minimum_t * 2.0);
        for _ in 0..128 {
            if kernel::state(lambda, hi)?.d >= target {
                break;
            }
            hi *= 2.0;
        }
        if kernel::state(lambda, hi)?.d < target {
            return Err(MaugisError::Unrepresentable {
                field: "contact envelope bracket",
            });
        }
        let maximum_t = kernel::invert(lambda, minimum_t, hi, target, false)?.t;
        let upper = kernel::state(lambda, maximum_t)?;
        let a = upper.a * radius_scale_m;
        let c = upper.c * radius_scale_m;
        let y = elastic.limiting_yield_stress_pa();
        for (ratio, value, limit) in [
            (
                "range/radius",
                interaction_range_m / effective_radius_m,
                envelope.maximum_interaction_range_to_radius,
            ),
            (
                "cohesive radius/radius",
                c / effective_radius_m,
                envelope.maximum_interaction_radius_to_radius,
            ),
            (
                "cohesive radius/depth",
                c / envelope.half_space_depth_m,
                envelope.maximum_interaction_radius_to_depth,
            ),
            (
                "cohesive radius/layer",
                c / envelope.layer_thickness_m,
                envelope.maximum_interaction_radius_to_layer,
            ),
            (
                "nominal strain",
                a / effective_radius_m,
                envelope.maximum_nominal_strain,
            ),
            (
                "repulsive pressure/yield",
                (2.0 / PI) * (e / y) * (a / effective_radius_m),
                envelope.maximum_repulsive_pressure_to_yield,
            ),
            (
                "adhesive stress/yield",
                (w / interaction_range_m) / y,
                envelope.maximum_adhesive_stress_to_yield,
            ),
        ] {
            if !value.is_finite() || value > limit {
                return Err(MaugisError::OutsideApplicability {
                    ratio,
                    value,
                    limit,
                });
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.maugis-dugdale-sphere.v1");
        hash.update(normal.identity().as_bytes());
        for v in [
            effective_radius_m,
            interaction_range_m,
            envelope.maximum_indentation_m,
            envelope.half_space_depth_m,
            envelope.layer_thickness_m,
            envelope.maximum_interaction_radius_to_radius,
            envelope.maximum_interaction_radius_to_depth,
            envelope.maximum_interaction_radius_to_layer,
            envelope.maximum_nominal_strain,
            envelope.maximum_repulsive_pressure_to_yield,
            envelope.maximum_adhesive_stress_to_yield,
            envelope.maximum_interaction_range_to_radius,
        ] {
            hash.update(&v.to_bits().to_le_bytes());
        }
        let model = Self {
            normal,
            radius_m: effective_radius_m,
            interaction_range_m,
            envelope,
            identity: hash.finalize(),
            radius_scale_m,
            indentation_scale_m,
            force_scale_n,
            lambda,
            minimum_t,
            load_fold_t,
            maximum_t,
        };
        model.response(kernel::state(lambda, minimum_t)?)?;
        model.response(upper)?;
        Ok(model)
    }
    /// Retained original normal source with all source/query receipts.
    #[must_use]
    pub const fn normal_state(&self) -> &BoundNormalInterfaceState {
        &self.normal
    }
    /// Identity of this model, source, geometry and applicability envelope.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Effective sphere radius, m.
    #[must_use]
    pub const fn effective_radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Explicit interaction range z0, m.
    #[must_use]
    pub const fn interaction_range_m(&self) -> f64 {
        self.interaction_range_m
    }
    /// Dimensionless Maugis lambda (K=4E*/3 convention).
    #[must_use]
    pub const fn maugis_parameter(&self) -> f64 {
        self.lambda
    }
    /// Authored applicability envelope.
    #[must_use]
    pub const fn envelope(&self) -> MaugisEnvelope {
        self.envelope
    }
    /// Minimum indentation on the stable contacting branch, m.
    pub fn minimum_attached_indentation_m(&self) -> f64 {
        kernel::state(self.lambda, self.minimum_t)
            .expect("validated branch endpoint")
            .d
            * self.indentation_scale_m
    }
    /// Minimum load on the stable force-controlled contacting branch, N.
    pub fn force_controlled_pull_off_n(&self) -> f64 {
        kernel::state(self.lambda, self.load_fold_t)
            .expect("validated load endpoint")
            .p
            * self.force_scale_n
    }
    fn state_at_indentation(&self, indentation_m: f64) -> Result<kernel::State, MaugisError> {
        if !indentation_m.is_finite() {
            return Err(MaugisError::InvalidInput {
                field: "signed indentation",
            });
        }
        let min = self.minimum_attached_indentation_m();
        if indentation_m < min {
            return Err(MaugisError::NoContactEquilibrium {
                indentation_m,
                minimum_m: min,
            });
        }
        if indentation_m > self.envelope.maximum_indentation_m {
            return Err(MaugisError::OutsideApplicability {
                ratio: "indentation",
                value: indentation_m,
                limit: self.envelope.maximum_indentation_m,
            });
        }
        // Preserve the exact endpoint despite physical-unit division rounding.
        if indentation_m == min {
            return kernel::state(self.lambda, self.minimum_t);
        }
        kernel::invert(
            self.lambda,
            self.minimum_t,
            self.maximum_t,
            indentation_m / self.indentation_scale_m,
            false,
        )
    }
    /// Stable displacement-controlled contact, including negative indentation.
    pub fn evaluate_attached(&self, indentation_m: f64) -> Result<MaugisResponse, MaugisError> {
        self.response(self.state_at_indentation(indentation_m)?)
    }
    /// Stable force-controlled contact; a more tensile request refuses.
    pub fn evaluate_contact_load(&self, load_n: f64) -> Result<MaugisResponse, MaugisError> {
        if !load_n.is_finite() {
            return Err(MaugisError::InvalidInput {
                field: "applied load",
            });
        }
        let min = self.force_controlled_pull_off_n();
        if load_n < min {
            return Err(MaugisError::BelowPullOff {
                load_n,
                pull_off_n: min,
            });
        }
        let upper = kernel::state(self.lambda, self.maximum_t)?.p * self.force_scale_n;
        if load_n > upper {
            return Err(MaugisError::OutsideApplicability {
                ratio: "applied load",
                value: load_n,
                limit: upper,
            });
        }
        let s = if load_n == min {
            kernel::state(self.lambda, self.load_fold_t)?
        } else {
            kernel::invert(
                self.lambda,
                self.load_fold_t,
                self.maximum_t,
                load_n / self.force_scale_n,
                true,
            )?
        };
        self.response(s)
    }
    fn response(&self, s: kernel::State) -> Result<MaugisResponse, MaugisError> {
        let a = s.a * self.radius_scale_m;
        let c = s.c * self.radius_scale_m;
        let endpoint = if s.t == self.minimum_t {
            Some(if s.t == 0.0 {
                MaugisEndpoint::ZeroArea
            } else {
                MaugisEndpoint::DisplacementFold
            })
        } else {
            None
        };
        let tangent = if endpoint.is_some() {
            None
        } else {
            Some((s.ps / s.ds) * (self.force_scale_n / self.indentation_scale_m))
        };
        let area = PI * a * a;
        let force = s.p * self.force_scale_n;
        let repulsive = s.a * s.a * s.a * self.force_scale_n;
        let indentation = s.d * self.indentation_scale_m;
        let peak = (2.0 / PI) * self.normal.elastic().reduced_modulus_pa() * (a / self.radius_m);
        let sigma = self.normal.adhesion_energy_j_per_m2() / self.interaction_range_m;
        if [a, c, area, force, repulsive, indentation, peak, sigma]
            .iter()
            .any(|v| !v.is_finite())
            || tangent.is_some_and(|v| !v.is_finite())
            || (s.a > 0.0 && (a <= 0.0 || area <= 0.0 || repulsive <= 0.0))
        {
            return Err(MaugisError::Unrepresentable {
                field: "dimensional Maugis response",
            });
        }
        let authority = match self
            .normal
            .elastic()
            .interface()
            .interface()
            .provenance()
            .authority()
        {
            InputAuthority::SyntheticFixture => InputAuthority::SyntheticFixture,
            _ => InputAuthority::Estimated,
        };
        Ok(MaugisResponse {
            normal_state_identity: self.normal.identity(),
            indentation_m: indentation,
            contact_radius_m: a,
            cohesive_radius_m: c,
            contact_area_m2: area,
            repulsive_normal_force_n: repulsive,
            normal_force_n: force,
            tangent_n_per_m: tangent,
            endpoint,
            authority,
            adhesive_stress_pa: sigma,
            hertz_peak_pressure_pa: peak,
            annulus_root_m: s.l * self.radius_scale_m,
        })
    }
}
fn positive(v: f64, field: &'static str) -> Result<(), MaugisError> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(MaugisError::InvalidInput { field })
    }
}
