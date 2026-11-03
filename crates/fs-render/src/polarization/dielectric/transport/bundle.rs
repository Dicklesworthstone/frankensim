//! Differentiated specular paths with native polarized Fresnel power.
//!
//! The geometric bundle retains the original emission solid angle through
//! smooth fields and explicit interfaces. A unit-intensity Stokes state in its
//! transported frame and a separate logarithmic interface-power factor retain
//! partial/circular polarization across every face without confusing Fresnel
//! loss with Beer--Lambert bulk optical depth. No Monte Carlo branch probability
//! or radiance eta-squared factor is part of this deterministic branch weight.
//!
//! Explicit continuations and the native scene's automatic boundary walk use
//! the same physical transfer. This owner does not find connecting rays,
//! differentiate rough scattering or enumerate multiple images. Its irradiance
//! is the contribution of this one branch; admitted branches add separately.

use super::super::mirrors::ConductorBeamError;
use super::*;
use crate::conductor::{CoherentDielectricCoating, ConductorOptics};
use crate::gradient_index::bundle::event::{
    BundleBoundaryPropagation, BundleEventOptions, propagate_seeded_bundle_to_boundary_with_control,
};
use crate::gradient_index::bundle::interface::{BundleInterfaceError, InterfaceOptions};
use crate::gradient_index::bundle::{
    BundleError, BundleOptions, ProjectionOptions, RayBundlePropagation, SecondOrderIndexField,
    propagate_seeded_bundle_with_control,
};
use crate::gradient_index::event::{BoundaryOracle, GradientDeparture};
use crate::instances::proximity::InstanceNormalDifferential;
use fs_geom::NormalDifferential;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolarizedBundleError {
    Bundle(BundleError),
    Interface(BundleInterfaceError),
    Beam(BeamError),
    Conductor(ConductorBeamError),
    MissingInterfaceTransfer,
    InvalidSourcePolarization,
    ContactPositionMismatch,
    UnsupportedSurface,
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
impl From<ConductorBeamError> for PolarizedBundleError {
    fn from(error: ConductorBeamError) -> Self {
        Self::Conductor(error)
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
    reverse_branch_log_probability: f64,
}

#[derive(Clone, Copy)]
enum InterfaceOptics<'a> {
    Bare,
    Film(&'a CoherentDielectricCoating),
    Surface(&'a DielectricSurface),
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
            reverse_branch_log_probability: 0.0,
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

    /// Log of the native scalar proposal for traversing the selected smooth
    /// interface chain backward. This is a product of discrete branch
    /// probabilities, not a solid-angle density or physical Mueller power.
    /// Bulk absorption, source departure and receiver response are excluded.
    /// Transmitted coatings use the reversed physical layer order and local
    /// indices; a smooth conductor has deterministic reflection probability one.
    pub const fn reverse_branch_log_probability(&self) -> f64 {
        self.reverse_branch_log_probability
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

    /// The scene consumes one native augmented boundary solve without
    /// reconstructing a second ray or dropping accumulated interface power.
    /// Geometric replacement stays private and is constructed from this seed.
    pub(crate) fn continue_to_boundary_with_control<
        F: SecondOrderIndexField + ?Sized,
        O: BoundaryOracle + ?Sized,
    >(
        &self,
        field: &F,
        oracle: &mut O,
        max_distance_m: f64,
        departure: Option<GradientDeparture>,
        options: BundleEventOptions,
        control: impl FnMut() -> bool,
    ) -> Result<(Self, BundleBoundaryPropagation), PolarizedBundleError> {
        let segment = propagate_seeded_bundle_to_boundary_with_control(
            field,
            oracle,
            &self.bundle.continuation(),
            max_distance_m,
            departure,
            options,
            control,
        )?;
        Ok((
            Self {
                bundle: segment.bundle,
                ..*self
            },
            segment,
        ))
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
        self.interface(surface, target, InterfaceOptics::Bare, event, options)
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
        self.interface(
            surface,
            target,
            InterfaceOptics::Film(coating),
            event,
            options,
        )
    }

    /// Preserve an authored smooth material's complete coherent multilayer.
    /// `surface.normal` points outward, so exiting rays use the native reversed
    /// layer order. No stochastic selection or camera eta factor is introduced.
    pub fn smooth_surface_interface(
        &self,
        surface: NormalDifferential,
        target: GradientIndexSample,
        material: &DielectricSurface,
        event: DielectricEvent,
        options: InterfaceOptions,
    ) -> Result<Self, PolarizedBundleError> {
        if !material.is_delta() {
            return Err(PolarizedBundleError::UnsupportedSurface);
        }
        self.interface(
            surface,
            target,
            InterfaceOptics::Surface(material),
            event,
            options,
        )
    }

    /// Native smooth-metal reflection with the same differentiated geometric
    /// reflection map. Complex optical constants, coating and Stokes phase
    /// remain owned by the existing polarized conductor response.
    pub fn smooth_conductor_interface(
        &self,
        surface: NormalDifferential,
        optics: &ConductorOptics,
        options: InterfaceOptions,
    ) -> Result<Self, PolarizedBundleError> {
        let bundle = self
            .bundle
            .smooth_interface(
                surface,
                self.bundle.propagation().final_sample,
                DielectricEvent::Reflection,
                options,
            )?
            .at_endpoint();
        let output = self
            .endpoint_beam()?
            .reflect_conductor(surface.normal, optics)?
            .beam
            .with_reference_axis(bundle.propagation().transverse_axis)?;
        self.retain_interface_power(bundle, output.power())
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
        optics: InterfaceOptics<'_>,
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
        let incident = self.endpoint_beam()?;
        let normal = incident.incident_normal(surface.normal)?;
        let cosine = (-normal.dot(incident.frame().direction())).clamp(0.0, 1.0);
        let eta_i = previous.final_sample.phase_index();
        let eta_t = target.phase_index();
        let response = match optics {
            InterfaceOptics::Bare => fresnel_dielectric_polarized(cosine, eta_i, eta_t),
            InterfaceOptics::Film(coating) => super::super::fresnel_coated_dielectric_polarized(
                cosine,
                eta_i,
                eta_t,
                wavelength_nm,
                coating,
            ),
            InterfaceOptics::Surface(material) => material
                .interface(wavelength_nm, surface.normal.dot(previous.ray.dir) > 0.0)
                .and_then(|interface| interface.fresnel_polarized(cosine, eta_i, eta_t)),
        }
        .map_err(BeamError::from)?;
        let matrix = match event {
            DielectricEvent::Reflection => response.reflection,
            DielectricEvent::Transmission => response.transmission,
        };
        if matrix.0[0][0] == 0.0 {
            return Err(BundleInterfaceError::UnavailableBranch.into());
        }
        let reverse_probability = reverse_probability(
            optics,
            response,
            normal,
            surface.normal.dot(previous.ray.dir) > 0.0,
            bundle.propagation().ray.dir,
            wavelength_nm,
            eta_i,
            eta_t,
            event,
        )?;
        let reverse_branch_log_probability =
            self.reverse_branch_log_probability + det::ln(reverse_probability);
        if !reverse_branch_log_probability.is_finite() || reverse_branch_log_probability > 0.0 {
            return Err(PolarizedBundleError::NumericalRange);
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
        let mut result = self.retain_interface_power(bundle, output.power())?;
        result.reverse_branch_log_probability = reverse_branch_log_probability;
        Ok(result)
    }

    fn endpoint_beam(&self) -> Result<PolarizedBeam, PolarizedBundleError> {
        let propagation = self.bundle.propagation();
        Ok(PolarizedBeam::try_new(
            propagation.ray.origin,
            PolarizationFrame::try_new(propagation.ray.dir, propagation.transverse_axis)?,
            self.bundle.wavelength_nm(),
            BeamMedium::GradientIndex(LocalGradientMedium {
                position: propagation.ray.origin,
                wavelength_nm: self.bundle.wavelength_nm(),
                sample: propagation.final_sample,
            }),
            self.interface_state,
        )?)
    }

    fn retain_interface_power(
        &self,
        bundle: RayBundlePropagation,
        power: Stokes,
    ) -> Result<Self, PolarizedBundleError> {
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
            reverse_branch_log_probability: self.reverse_branch_log_probability,
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

// The scalar sampler has full polarization support, but its proposal is not
// the polarization-dependent physical power retained above. In particular,
// an absorbing multilayer can have different survival-conditioned proposals
// in opposite directions even though its transmitted power is reciprocal.
#[allow(clippy::too_many_arguments)]
fn reverse_probability(
    optics: InterfaceOptics<'_>,
    forward: super::super::PolarizedFresnel,
    incident_normal: Vec3,
    forward_reverse_layers: bool,
    outgoing: Vec3,
    wavelength_nm: f64,
    eta_i: f64,
    eta_t: f64,
    event: DielectricEvent,
) -> Result<f64, PolarizedBundleError> {
    let transmission = event == DielectricEvent::Transmission;
    let reverse_layers = forward_reverse_layers ^ transmission;
    let cosine = incident_normal.dot(outgoing).abs().clamp(0.0, 1.0);
    let (incident_index, target_index) = if transmission { (eta_t, eta_i) } else { (eta_i, eta_t) };
    let response = if !transmission {
        forward
    } else {
        match optics {
            InterfaceOptics::Bare => {
                fresnel_dielectric_polarized(cosine, incident_index, target_index)
            }
            InterfaceOptics::Film(coating) => super::super::fresnel_coated_dielectric_polarized(
                cosine, incident_index, target_index, wavelength_nm, coating,
            ),
            InterfaceOptics::Surface(material) => material
                .interface(wavelength_nm, reverse_layers)
                .and_then(|interface| {
                    interface.fresnel_polarized(cosine, incident_index, target_index)
                }),
        }
        .map_err(BeamError::from)?
    };
    let probabilities = match optics {
        InterfaceOptics::Bare => [response.scalar.reflectance, response.scalar.transmittance],
        InterfaceOptics::Film(coating) if coating.thickness_m() == 0.0 => {
            [response.scalar.reflectance, response.scalar.transmittance]
        }
        InterfaceOptics::Film(_) => {
            let reflected = response.reflection.0[0][0];
            let transmitted = response.transmission.0[0][0];
            let total = reflected + transmitted;
            [reflected / total, transmitted / total]
        }
        InterfaceOptics::Surface(material) => material
            .interface(wavelength_nm, reverse_layers)
            .map_err(BeamError::from)?
            .event_probabilities(response.scalar),
    };
    let probability = probabilities[usize::from(transmission)];
    if !probability.is_finite() || probability <= 0.0 || probability > 1.0 {
        return Err(PolarizedBundleError::NumericalRange);
    }
    Ok(probability)
}

#[cfg(test)]
mod tests;
