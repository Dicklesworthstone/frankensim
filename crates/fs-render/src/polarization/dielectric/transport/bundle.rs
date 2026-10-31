//! Differentiated specular paths with native polarized Fresnel power.
//!
//! The geometric bundle retains the original emission solid angle through
//! smooth fields and explicit interfaces. A unit-intensity Stokes state in its
//! transported frame and a separate logarithmic interface-power factor retain
//! partial/circular polarization across every face without confusing Fresnel
//! loss with Beer--Lambert bulk optical depth. No Monte Carlo branch probability
//! or radiance eta-squared factor is part of this deterministic branch weight.
//!
//! Every segment and boundary is supplied explicitly. This owner does not find
//! connecting rays, establish visibility, differentiate rough scattering or
//! enumerate multiple images. Its irradiance is the contribution of this one
//! branch; different admitted branches can be added by their caller.

use super::*;
use crate::conductor::CoherentDielectricCoating;
use crate::gradient_index::bundle::interface::{BundleInterfaceError, InterfaceOptions};
use crate::gradient_index::bundle::{
    BundleError, BundleOptions, ProjectionOptions, RayBundlePropagation, SecondOrderIndexField,
    propagate_seeded_bundle_with_control,
};
use crate::instances::proximity::InstanceNormalDifferential;
use fs_geom::NormalDifferential;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolarizedBundleError {
    Bundle(BundleError),
    Interface(BundleInterfaceError),
    Beam(BeamError),
    MissingInterfaceTransfer,
    InvalidSourcePolarization,
    ContactPositionMismatch,
    NumericalRange,
}
impl From<BundleError> for PolarizedBundleError {
    fn from(error: BundleError) -> Self {
        Self::Bundle(error)
    }
}
impl From<BundleInterfaceError> for PolarizedBundleError {
    fn from(error: BundleInterfaceError) -> Self {
        Self::Interface(error)
    }
}
impl From<BeamError> for PolarizedBundleError {
    fn from(error: BeamError) -> Self {
        Self::Beam(error)
    }
}
impl core::fmt::Display for PolarizedBundleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized refractive bundle: {self:?}")
    }
}
impl core::error::Error for PolarizedBundleError {}

/// One differentiated optical branch, including all its interface power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedRayBundle {
    bundle: RayBundlePropagation,
    interface_state: Stokes,
    interface_log_power: f64,
}

impl PolarizedRayBundle {
    /// Start from a bundle with no preceding discrete interface events.
    pub fn unpolarized(bundle: RayBundlePropagation) -> Result<Self, PolarizedBundleError> {
        Self::try_new(bundle, Stokes::unpolarized(1.0))
    }

    /// Source polarization has unit intensity in the bundle's transported
    /// transverse frame. Its absolute spectral intensity stays with the source.
    /// A geometric bundle that already crossed interfaces cannot reconstruct
    /// their omitted power transfer and is refused.
    pub fn try_new(
        bundle: RayBundlePropagation,
        source_polarization: Stokes,
    ) -> Result<Self, PolarizedBundleError> {
        if bundle.interface_events() != 0 {
            return Err(PolarizedBundleError::MissingInterfaceTransfer);
        }
        if source_polarization.i != 1.0
            || Stokes::try_new(
                1.0,
                source_polarization.q,
                source_polarization.u,
                source_polarization.v,
            )
            .is_none()
        {
            return Err(PolarizedBundleError::InvalidSourcePolarization);
        }
        let interface_state = physical_roundoff(source_polarization, 1.0)?;
        Ok(Self {
            bundle,
            interface_state,
            interface_log_power: 0.0,
        })
    }

    pub const fn propagation(&self) -> &RayBundlePropagation {
        &self.bundle
    }
    /// Unit-intensity polarization in the current bundle's transverse frame.
    /// For a completely dark branch this frame state is conventional; its log
    /// power is negative infinity and every irradiance contribution is zero.
    pub const fn interface_state(&self) -> Stokes {
        self.interface_state
    }
    /// Log of accumulated physical interface power, excluding bulk absorption.
    pub const fn interface_log_power(&self) -> f64 {
        self.interface_log_power
    }

    pub fn continue_in<F: SecondOrderIndexField + ?Sized>(
        &self,
        field: &F,
        distance_m: f64,
        options: BundleOptions,
    ) -> Result<Self, PolarizedBundleError> {
        self.continue_in_with_control(field, distance_m, options, || true)
    }

    /// One smooth segment. The original field sample must match the current
    /// endpoint; step budgets/cancellation belong to the native bundle owner.
    pub fn continue_in_with_control<F: SecondOrderIndexField + ?Sized>(
        &self,
        field: &F,
        distance_m: f64,
        options: BundleOptions,
        control: impl FnMut() -> bool,
    ) -> Result<Self, PolarizedBundleError> {
        let bundle = propagate_seeded_bundle_with_control(
            field,
            &self.bundle.continuation(),
            distance_m,
            options,
            control,
        )?;
        Ok(Self { bundle, ..*self })
    }

    /// Select one physical smooth-interface branch at the current endpoint.
    /// Curvature and spatial index gradients affect its angular derivatives;
    /// the native s/p Mueller response affects its power and polarization.
    pub fn smooth_interface(
        &self,
        surface: NormalDifferential,
        target: GradientIndexSample,
        event: DielectricEvent,
        options: InterfaceOptions,
    ) -> Result<Self, PolarizedBundleError> {
        self.interface(surface, target, None, event, options)
    }

    /// A native coherent transparent coating changes the physical Mueller
    /// response while the two bulk media retain the same geometric Snell map.
    /// Film thickness remains the native collapsed-interface optical model;
    /// no geometric film travel or group delay is invented.
    pub fn smooth_coated_interface(
        &self,
        surface: NormalDifferential,
        target: GradientIndexSample,
        coating: &CoherentDielectricCoating,
        event: DielectricEvent,
        options: InterfaceOptions,
    ) -> Result<Self, PolarizedBundleError> {
        self.interface(surface, target, Some(coating), event, options)
    }

    /// Use a differential admitted by the original prepared geometry owner.
    /// The query must refer to this exact endpoint; it is never moved to the
    /// surface witness or reused at another point.
    pub fn smooth_native_interface(
        &self,
        contact: &InstanceNormalDifferential,
        target: GradientIndexSample,
        event: DielectricEvent,
        options: InterfaceOptions,
    ) -> Result<Self, PolarizedBundleError> {
        if contact.point != self.bundle.propagation().ray.origin {
            return Err(PolarizedBundleError::ContactPositionMismatch);
        }
        self.smooth_interface(contact.differential, target, event, options)
    }

    fn interface(
        &self,
        surface: NormalDifferential,
        target: GradientIndexSample,
        coating: Option<&CoherentDielectricCoating>,
        event: DielectricEvent,
        options: InterfaceOptions,
    ) -> Result<Self, PolarizedBundleError> {
        let seed = self
            .bundle
            .smooth_interface(surface, target, event, options)?;
        let bundle = seed.at_endpoint();
        let previous = self.bundle.propagation();
        let origin = previous.ray.origin;
        let wavelength_nm = self.bundle.wavelength_nm();
        let medium = |sample| {
            BeamMedium::GradientIndex(LocalGradientMedium {
                position: origin,
                wavelength_nm,
                sample,
            })
        };
        let incident = PolarizedBeam::try_new(
            origin,
            PolarizationFrame::try_new(previous.ray.dir, previous.transverse_axis)?,
            wavelength_nm,
            medium(previous.final_sample),
            self.interface_state,
        )?;
        let normal = incident.incident_normal(surface.normal)?;
        let cosine = (-normal.dot(incident.frame().direction())).clamp(0.0, 1.0);
        let eta_i = previous.final_sample.phase_index();
        let eta_t = target.phase_index();
        let response = match coating {
            None => fresnel_dielectric_polarized(cosine, eta_i, eta_t),
            Some(coating) => super::super::fresnel_coated_dielectric_polarized(
                cosine,
                eta_i,
                eta_t,
                wavelength_nm,
                coating,
            ),
        }
        .map_err(BeamError::from)?;
        let matrix = match event {
            DielectricEvent::Reflection => response.reflection,
            DielectricEvent::Transmission => response.transmission,
        };
        if matrix.0[0][0] == 0.0 {
            return Err(BundleInterfaceError::UnavailableBranch.into());
        }
        // Reuse the actual native frame rotation, Mueller application and
        // physical-roundoff admission, with deterministic physical power.
        let output = incident
            .scatter(
                normal,
                bundle.propagation().ray.dir,
                matrix,
                event,
                medium(target),
            )?
            .with_reference_axis(bundle.propagation().transverse_axis)?;
        let power = output.power();
        if power.i > 1.0 + 512.0 * f64::EPSILON {
            return Err(PolarizedBundleError::NumericalRange);
        }
        let (interface_state, interface_log_power) = if power.i == 0.0 {
            (Stokes::unpolarized(1.0), f64::NEG_INFINITY)
        } else {
            let state = physical_roundoff(
                Stokes {
                    i: 1.0,
                    q: power.q / power.i,
                    u: power.u / power.i,
                    v: power.v / power.i,
                },
                1.0,
            )?;
            let log_power = self.interface_log_power + det::ln(power.i.min(1.0));
            if log_power.is_nan() || log_power > 0.0 {
                return Err(PolarizedBundleError::NumericalRange);
            }
            (state, log_power)
        };
        Ok(Self {
            bundle,
            interface_state,
            interface_log_power,
        })
    }

    /// Contribution of this optical branch to point-source irradiance. The
    /// source intensity [W/sr] (or W/(sr nm)) belongs to the retained launch.
    /// The area Jacobian includes receiver obliquity; physical Fresnel transfer
    /// and bulk absorption are applied once, including representable results
    /// whose separately evaluated attenuation would underflow.
    pub fn point_source_irradiance(
        &self,
        normal: Vec3,
        options: ProjectionOptions,
        intensity: f64,
    ) -> Result<f64, PolarizedBundleError> {
        if !intensity.is_finite() || intensity < 0.0 {
            return Err(BundleError::InvalidIntensity.into());
        }
        let projection = self.bundle.project_at_endpoint(normal, options)?;
        if intensity == 0.0 || self.interface_log_power == f64::NEG_INFINITY {
            return Ok(0.0);
        }
        let log_attenuation = self.interface_log_power - self.bundle.propagation().optical_depth;
        let area = projection.area_per_solid_angle_m2();
        let value = (intensity * det::exp(log_attenuation)) / area;
        if value.is_finite() && value > 0.0 {
            return Ok(value);
        }
        let value = det::exp(det::ln(intensity) + log_attenuation - det::ln(area));
        if !value.is_finite() {
            return Err(PolarizedBundleError::NumericalRange);
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
