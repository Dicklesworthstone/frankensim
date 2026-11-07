//! Source-owned directional friction on the original Hertz pressure footprint.
//!
//! Normal geometry, load, pressure and normal loss come from `SourcePointContact`.
//! The original directional source owner refreshes both coefficient pairs, and
//! the native elliptical Hertz owner integrates sliding and spin. No circular
//! equivalent radius, new coefficient selector or independent pressure model is
//! inserted between those owners.
//!
//! Conditions are uniform over this instantaneous patch. Two explicitly named
//! temperature coordinates read the original ordered bulk temperatures. Local
//! slip speed is uniform only for pure translation: a speed-dependent source
//! refuses nonzero spin rather than sampling every point at the centre speed.
//! `resolved` evaluates speed-dependent spinning traction at actual local slip.
//! Spatially varying temperature, distributed path memory, static reactions and
//! physical-time advancement require separate owners. A uniform response's
//! returned checkpoint is only a source-refresh candidate.

/// Distributed local-speed source queries on the actual Hertz pressure field.
pub mod resolved;

use super::impulse::{SourceImpulseCoefficientRole, SourceImpulseError};
use super::*;
use crate::normal_patch::source::SourcePointContact;
use crate::normal_patch::{NormalPatchError, PointNormalPatchReceipt};
use fs_tribo::sliding_patch::HertzSlidingWrench;
use fs_tribo::sliding_patch::elliptic::{
    EllipticHertzError, EllipticHertzPatchParameters, EllipticHertzQuadrature,
    EllipticHertzSlidingPatch,
};

/// Actual normal sample and prescribed rigid slip in one tangent frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceEllipticHertzInput<'a> {
    pub normal_state_id: &'a str,
    pub indentation_m: f64,
    pub indentation_rate_m_s: f64,
    /// Original normal owner's sample duration. Its reported loss is still
    /// instantaneous power times this duration, not a path-integrated loss.
    pub sample_duration_s: f64,
    /// Major footprint direction for an ellipse; arbitrary nonzero direction
    /// for a circle. This is independent of the retained material texture.
    pub footprint_axis: [f64; 2],
    pub slip_velocity_m_s: [f64; 2],
    pub spin_rad_s: f64,
    /// Optional original source coordinates for bulk A and B temperatures.
    /// Every additional owned temperature coordinate must be assigned exactly
    /// once. Their template values are replaced by the admitted normal bulk
    /// states' actual T values; no surface mean is substituted.
    pub surface_temperature_axes: [Option<&'a str>; 2],
}

/// One complete trial: original normal response, original source candidate and
/// the native resultant acting on the side of the supplied relative velocity.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceEllipticHertzReport {
    pub normal_model_identity: ContentHash,
    pub normal: PointNormalPatchReceipt,
    pub state: SourceThermalEllipticState,
    /// None at zero indentation: an absent footprint is not assigned a radius.
    pub patch: Option<EllipticHertzSlidingPatch>,
    pub wrench: HertzSlidingWrench,
}

/// All original refusal surfaces, with no published source or mechanical state.
#[derive(Debug)]
pub enum SourceEllipticHertzError {
    Source(SourceThermalCoulombError),
    Sample(SourceImpulseError),
    Normal(NormalPatchError),
    Native(EllipticHertzError),
    Invalid(&'static str),
    Cancelled,
}
impl From<SourceThermalCoulombError> for SourceEllipticHertzError {
    fn from(error: SourceThermalCoulombError) -> Self {
        Self::Source(error)
    }
}
impl From<TriboError> for SourceEllipticHertzError {
    fn from(error: TriboError) -> Self {
        Self::Source(SourceThermalCoulombError::Tribo(error))
    }
}
impl From<SourceImpulseError> for SourceEllipticHertzError {
    fn from(error: SourceImpulseError) -> Self {
        match error {
            SourceImpulseError::Cancelled => Self::Cancelled,
            other => Self::Sample(other),
        }
    }
}
impl From<NormalPatchError> for SourceEllipticHertzError {
    fn from(error: NormalPatchError) -> Self {
        Self::Normal(error)
    }
}
impl From<EllipticHertzError> for SourceEllipticHertzError {
    fn from(error: EllipticHertzError) -> Self {
        match error {
            EllipticHertzError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl core::fmt::Display for SourceEllipticHertzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source elliptical Hertz sliding: {self:?}")
    }
}
impl std::error::Error for SourceEllipticHertzError {}

impl SourceThermalEllipticCoulomb {
    /// Evaluate original normal and orthotropic sliding laws on one actual patch.
    ///
    /// Both owners must bind the same complete ordered interface and authority.
    /// The primary friction temperature must equal the original normal interface
    /// temperature. Additional surface-temperature axes read the ordered bulk
    /// states through `surface_temperature_axes`; all other external coordinates
    /// retain their original source semantics, including supplied nominal pressure.
    /// That source pressure is not inferred from Hertz peak or mean pressure.
    ///
    /// The existing instantaneous coefficient sampler replaces every declared
    /// speed with actual translation speed and keeps every accepted path value.
    /// Consequently speed-dependent coefficients require exactly zero spin here.
    /// For speed-independent coefficients, translation and spin are integrated
    /// together by the original elliptical pressure law with its bounded angular
    /// quadrature. The current transported material axis, not its initial axis,
    /// is passed to that law. No slip, temperature history or physical time is
    /// advanced; the enclosing owner must accept this source candidate with its
    /// complete mechanical/thermal trial. Normal and tangential loss stay distinct.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_hertz_patch(
        &self,
        state: &SourceThermalEllipticState,
        normal: &SourcePointContact,
        input: SourceEllipticHertzInput<'_>,
        conditions: &SourceFrictionConditions,
        quadrature: EllipticHertzQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceEllipticHertzReport, SourceEllipticHertzError> {
        let actual = self.hertz_conditions(state, normal, input, conditions, cancelled)?;
        if !(1..=256).contains(&quadrature.angular_panels) {
            return Err(SourceEllipticHertzError::Invalid(
                "bounded angular quadrature required",
            ));
        }
        let axes = self.principal[0].query_axes();
        if input.spin_rad_s != 0.0
            && axes
                .iter()
                .any(|axis| axis.kind == FrictionQueryAxisKind::SlidingSpeed)
        {
            return Err(SourceEllipticHertzError::Invalid(
                "spinning speed-dependent friction requires resolved local speed queries",
            ));
        }
        poll(cancelled)?;
        let response = normal.evaluate(
            input.normal_state_id,
            input.indentation_m,
            input.indentation_rate_m_s,
            input.sample_duration_s,
        )?;
        let sample = self.sample_preimpact_coefficients(
            state,
            input.slip_velocity_m_s,
            &actual,
            SourceImpulseCoefficientRole::Kinetic,
            cancelled,
        )?;
        // The source checkpoint's owner retains its initial frame identity;
        // the instantaneous patch instead uses the currently transported axis.
        let law = EllipticCoulomb::new(
            sample.state.law.static_coefficients(),
            sample.state.law.kinetic_coefficients(),
            sample.state.constitutive.principal_axis(),
        )?;
        let semi_axes = response
            .elliptic_patch_axes
            .map_or([response.patch_radius_m; 2], |axes| {
                [axes.semi_major_axis_m, axes.semi_minor_axis_m]
            });
        let (patch, wrench) = if input.indentation_m == 0.0 {
            (
                None,
                HertzSlidingWrench {
                    force_n: [0.0; 2],
                    spin_moment_nm: 0.0,
                    dissipated_power_w: 0.0,
                    power_residual_w: 0.0,
                    at_rest: input.slip_velocity_m_s == [0.0; 2] && input.spin_rad_s == 0.0,
                },
            )
        } else {
            let patch = EllipticHertzSlidingPatch::new(
                EllipticHertzPatchParameters {
                    normal_force_n: response.normal_force_n,
                    semi_axes_m: semi_axes,
                    footprint_axis: input.footprint_axis,
                },
                law,
            )?;
            let wrench = patch.evaluate_cancellable(
                input.slip_velocity_m_s,
                input.spin_rad_s,
                quadrature,
                &mut *cancelled,
            )?;
            (Some(patch), wrench)
        };
        poll(cancelled)?;
        Ok(SourceEllipticHertzReport {
            normal_model_identity: normal.model().identity(),
            normal: response,
            state: sample.state,
            patch,
            wrench,
        })
    }

    fn hertz_conditions(
        &self,
        state: &SourceThermalEllipticState,
        normal: &SourcePointContact,
        input: SourceEllipticHertzInput<'_>,
        conditions: &SourceFrictionConditions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceFrictionConditions, SourceEllipticHertzError> {
        poll(cancelled)?;
        self.check(state)?;
        let elastic = normal.model().normal_state().elastic();
        if elastic.interface() != self.principal[0].interface()
            || conditions.contact_temperature_k != elastic.state_coordinate("T")
        {
            return Err(SourceEllipticHertzError::Invalid(
                "normal and friction require the same ordered interface, authority and interface temperature",
            ));
        }
        let accepted = state.sources[0].query();
        if conditions.query.axis_quantities() != accepted.axis_quantities()
            || conditions.query.axes().keys().ne(accepted.axes().keys())
        {
            return Err(SourceEllipticHertzError::Invalid(
                "original complete query schema required before replacing temperature placeholders",
            ));
        }
        let footprint_norm = fs_math::det::hypot(input.footprint_axis[0], input.footprint_axis[1]);
        if !input.spin_rad_s.is_finite()
            || input.slip_velocity_m_s.iter().any(|x| !x.is_finite())
            || input.footprint_axis.iter().any(|x| !x.is_finite())
            || !footprint_norm.is_finite()
            || footprint_norm == 0.0
        {
            return Err(SourceEllipticHertzError::Invalid(
                "finite rigid patch motion and nonzero footprint direction required",
            ));
        }
        let axes = self.principal[0].query_axes();
        let mut actual = conditions.clone();
        let bulk = [elastic.surface_a_state(), elastic.surface_b_state()];
        for (side, axis) in input.surface_temperature_axes.into_iter().enumerate() {
            let Some(axis) = axis else { continue };
            if !axes
                .iter()
                .any(|owned| owned.name == axis && owned.kind == FrictionQueryAxisKind::Temperature)
                || input.surface_temperature_axes[..side].contains(&Some(axis))
            {
                return Err(SourceEllipticHertzError::Invalid(
                    "each bulk temperature requires a distinct original owned temperature axis",
                ));
            }
            let temperature = bulk[side].resolved().state_coordinate("T").ok_or(
                SourceEllipticHertzError::Invalid("original bulk absolute temperature required"),
            )?;
            actual.query = actual
                .query
                .with(axis, temperature)
                .map_err(SourceThermalCoulombError::from)?;
        }
        if axes.iter().any(|axis| {
            axis.kind == FrictionQueryAxisKind::Temperature
                && !input
                    .surface_temperature_axes
                    .contains(&Some(axis.name.as_str()))
        }) {
            return Err(SourceEllipticHertzError::Invalid(
                "all additional temperature axes need an original ordered bulk owner",
            ));
        }
        Ok(actual)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourceEllipticHertzError> {
    if cancelled() {
        Err(SourceEllipticHertzError::Cancelled)
    } else {
        Ok(())
    }
}
