//! Source-pinned directional impulses sampled at actual preimpact slip speed.
//!
//! The caller explicitly chooses which original coefficient role defines the
//! fixed instantaneous ellipse. This is an impact constitutive assumption, not
//! an inferred compression history or a static/kinetic switch during impact.
//! Accumulated physical travel and native slip/temperature history do not advance.
//! Source receipts and the full static/kinetic law refresh together for later
//! sustained contact; impact loss remains unallocated mechanical loss.

use super::*;
use fs_tribo::tangential::elliptic::impulse::{
    EllipticImpulseControls, EllipticImpulseError, EllipticImpulseInput, EllipticImpulseReport,
};

/// Complete caller-owned thermal/external boundary for source contact phases.
///
/// Primary temperature/pressure must match the explicit physical variables.
/// Additional declared surface conditions come from the query. Motion owners
/// replace declared speed/path placeholders with their own physical values.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceFrictionConditions {
    pub query: QueryPoint,
    pub contact_temperature_k: Option<f64>,
    pub nominal_pressure_pa: Option<f64>,
}

/// Explicit choice of an existing source role for one fixed impact ellipse.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceImpulseCoefficientRole {
    /// Use each original kinetic-coefficient claim for the entire impact.
    Kinetic,
    /// Use each original static-coefficient claim for the entire impact.
    /// This is a caller-declared impact model, not automatic static capture.
    Static,
}

/// Original coefficient claims sampled before a simultaneous contact solve.
///
/// No normal impulse or mobility is needed to sample the declared impact law.
/// The consuming mechanical owner must solve and admit the full event before
/// publishing this candidate history. The retained texture is expressed in the
/// same tangent coordinates as `velocity_before_m_s` supplied to the sampler.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceImpulseSample {
    pub state: SourceThermalEllipticState,
    pub coefficient_role: SourceImpulseCoefficientRole,
    pub coefficients_mu: [f64; 2],
    pub queried_speed_m_s: f64,
}

/// Accepted native impulse and original source candidate; no borrowed mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceImpulseReport {
    pub impulse: EllipticImpulseReport,
    pub state: SourceThermalEllipticState,
    pub coefficient_role: SourceImpulseCoefficientRole,
    pub coefficients_mu: [f64; 2],
    pub queried_speed_m_s: f64,
}

#[derive(Debug)]
pub enum SourceImpulseError {
    Source(SourceThermalCoulombError),
    Native(EllipticImpulseError),
    Invalid(&'static str),
    Cancelled,
}
impl From<SourceThermalCoulombError> for SourceImpulseError {
    fn from(error: SourceThermalCoulombError) -> Self {
        Self::Source(error)
    }
}
impl From<TriboError> for SourceImpulseError {
    fn from(error: TriboError) -> Self {
        Self::Source(SourceThermalCoulombError::Tribo(error))
    }
}
impl From<EllipticImpulseError> for SourceImpulseError {
    fn from(error: EllipticImpulseError) -> Self {
        match error {
            EllipticImpulseError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl core::fmt::Display for SourceImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional preimpact impulse: {self:?}")
    }
}
impl std::error::Error for SourceImpulseError {}

impl SourceThermalEllipticCoulomb {
    /// Resolve one source-held impact ellipse without presupposing its impulse.
    ///
    /// Simultaneous normal/tangential contact owners can sample every ordered
    /// interface from its actual incoming velocity, then solve all contact
    /// impulses together. This performs the same original claim refresh as
    /// `solve_preimpact_impulse`; it neither invokes a mechanical solver nor
    /// synthesizes a normal load, pressure, time interval or travel distance.
    /// Every accepted path coordinate and native slip/texture history survive
    /// unchanged. The complete source state remains a private candidate until
    /// the caller admits the complete coupled mechanical response.
    pub fn sample_preimpact_coefficients(
        &self,
        state: &SourceThermalEllipticState,
        velocity_before_m_s: [f64; 2],
        conditions: &SourceFrictionConditions,
        coefficient_role: SourceImpulseCoefficientRole,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceImpulseSample, SourceImpulseError> {
        poll(cancelled)?;
        self.check(state)?;
        if velocity_before_m_s.into_iter().any(|v| !v.is_finite()) {
            return Err(SourceImpulseError::Invalid(
                "finite preimpact velocity required",
            ));
        }
        let speed = fs_math::det::hypot(velocity_before_m_s[0], velocity_before_m_s[1]);
        if !speed.is_finite() {
            return Err(SourceImpulseError::Invalid(
                "representable preimpact speed required",
            ));
        }
        let refresh = |axis: usize| -> Result<ResolvedInterfaceStatePoint, SourceImpulseError> {
            let previous = &state.sources[axis];
            let source = &self.principal[axis];
            let point = query_for_preimpact(
                &conditions.query,
                &previous.query(),
                source.query_axes(),
                speed,
            )?;
            Ok(source.source_at_query_conditions(
                previous,
                conditions.contact_temperature_k,
                conditions.nominal_pressure_pa,
                &point,
            )?)
        };
        let first = refresh(0)?;
        poll(cancelled)?;
        let second = refresh(1)?;
        poll(cancelled)?;
        let sources = [first, second];
        let (static_mu, kinetic_mu) =
            principal_coefficients(&self.principal, [&sources[0], &sources[1]])?;
        let (law, history) =
            state
                .law
                .continue_coefficients(&state.constitutive, static_mu, kinetic_mu)?;
        let coefficients_mu = match coefficient_role {
            SourceImpulseCoefficientRole::Kinetic => kinetic_mu,
            SourceImpulseCoefficientRole::Static => static_mu,
        };
        poll(cancelled)?;
        Ok(SourceImpulseSample {
            state: SourceThermalEllipticState {
                owner: self.identity,
                law,
                constitutive: history,
                sources,
            },
            coefficient_role,
            coefficients_mu,
            queried_speed_m_s: speed,
        })
    }

    /// Resolve original claims at preimpact speed and apply their chosen ellipse.
    ///
    /// Every declared SlidingSpeed coordinate uses the norm of the supplied
    /// actual preimpact tangent velocity. Every AccumulatedDistance coordinate
    /// retains its own accepted value exactly; template motion entries are only
    /// placeholders. The original query schema, fixed external coordinates,
    /// temperature/pressure agreement, source domains and claim pins are checked
    /// by the existing complete-query owner.
    ///
    /// Both principal source bundles are refreshed before the allocation-free
    /// native impulse solve. A temporary equal-coefficient native law expresses
    /// the explicitly selected impact role. The returned source state retains
    /// BOTH refreshed coefficient pairs for future sustained contact, along with
    /// unchanged native slip, texture and last finite-step temperature.
    ///
    /// No time interval, displacement, pressure from normal impulse, or impact
    /// heat allocation is inferred. Source and native failures or cancellation
    /// publish no candidate. Normal/tangent-coupled restitution remains the
    /// consuming mechanical owner's responsibility.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_preimpact_impulse(
        &self,
        state: &SourceThermalEllipticState,
        input: &EllipticImpulseInput,
        conditions: &SourceFrictionConditions,
        coefficient_role: SourceImpulseCoefficientRole,
        controls: EllipticImpulseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceImpulseReport, SourceImpulseError> {
        let sample = self.sample_preimpact_coefficients(
            state,
            input.velocity_before_m_s,
            conditions,
            coefficient_role,
            cancelled,
        )?;
        let (impact_law, impact_history) = sample.state.law.continue_coefficients(
            &sample.state.constitutive,
            sample.coefficients_mu,
            sample.coefficients_mu,
        )?;
        let impulse = impact_law.solve_impulse(&impact_history, input, controls, cancelled)?;
        poll(cancelled)?;
        Ok(SourceImpulseReport {
            impulse,
            state: sample.state,
            coefficient_role: sample.coefficient_role,
            coefficients_mu: sample.coefficients_mu,
            queried_speed_m_s: sample.queried_speed_m_s,
        })
    }
}

fn query_for_preimpact(
    template: &QueryPoint,
    accepted: &QueryPoint,
    axes: &[FrictionQueryAxis],
    speed: f64,
) -> Result<QueryPoint, SourceImpulseError> {
    // Check the caller's complete schema before replacing owned placeholders,
    // so filling a missing coordinate cannot conceal an incomplete boundary.
    if template.axis_quantities() != accepted.axis_quantities()
        || template.axes().keys().ne(accepted.axes().keys())
    {
        return Err(SourceImpulseError::Invalid(
            "original complete query schema required",
        ));
    }
    let mut point = template.clone();
    for axis in axes {
        let value = match axis.kind {
            FrictionQueryAxisKind::SlidingSpeed => speed,
            FrictionQueryAxisKind::AccumulatedDistance => accepted.axes()[&axis.name],
            _ => continue,
        };
        point = point
            .with(&axis.name, value)
            .map_err(SourceThermalCoulombError::from)?;
    }
    Ok(point)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourceImpulseError> {
    if cancelled() {
        Err(SourceImpulseError::Cancelled)
    } else {
        Ok(())
    }
}
