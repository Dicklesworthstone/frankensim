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
            Some(-PI*self.lambda*c*(c/(c-fold))*(self.force_scale_n/self.indentation_scale_m))
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
    if [result.elastic_energy_j, result.surface_energy_j, result.potential_j].iter().any(|v| !v.is_finite())
        || result.elastic_energy_j < 0.0 || result.surface_energy_j > 0.0 {
        Err(MaugisError::Unrepresentable { field: "absolute elastic/surface energy" })
    } else { Ok(()) }
}
