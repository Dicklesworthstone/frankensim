//! Source-bound, lossless coherent coatings on transmitting dielectric bodies.
//!
//! The optical kernel is shared with conductor multilayers, including scaled
//! evanescent propagation and critical-angle limits. Unlike opaque substrates,
//! transmitted fields retain every normalization factor and are normalized by
//! the outgoing/incident normal Poynting flux before Mueller conversion.
//! Layer order reverses for an inside-to-outside traversal. Fields use
//! exp(-i omega t), +Q along s and the native right-handed ray frames.

use super::{
    DielectricError, FresnelDielectric, RoughDielectricEvaluation, RoughDielectricSample,
    SmoothDielectricSample,
};
use crate::conductor::{
    CoherentDielectricCoating, CoherentPeriodicCoating, CoherentTransmission,
    MAX_COHERENT_STACK_LAYERS,
};
use crate::polarization::Mueller;
use crate::polarization::dielectric::{
    PolarizedDielectricSample, PolarizedFresnel, fresnel_dielectric_polarized,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::Vec3;
use fs_math::det;

/// Compact transparent coating on a dielectric surface. Every layer must
/// declare exactly zero extinction: absorbing coatings are explicitly refused,
/// never treated as `T=1-R`. Bulk Beer-Lambert absorption remains independent.
/// The source-bound Cauchy unit cell and 64-layer work bound are unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricCoating(CoherentPeriodicCoating);

impl DielectricCoating {
    /// Admit a fully declared periodic coating only when every film is lossless.
    pub fn try_new(coating: CoherentPeriodicCoating) -> Result<Self, DielectricError> {
        if coating
            .unit_cell()
            .iter()
            .any(|layer| layer.extinction_k() != 0.0)
        {
            return Err(DielectricError::InvalidCoating);
        }
        Ok(Self(coating))
    }

    /// Original authored cell, repeat count, physical thicknesses and source.
    #[must_use]
    pub const fn periodic(&self) -> &CoherentPeriodicCoating {
        &self.0
    }

    /// Exact authored optical law and source identity, distinct from bulk glass.
    #[must_use]
    pub fn content_identity(&self) -> ContentHash {
        let mut hasher = DomainHasher::new("org.frankensim.dielectric-coating.v1");
        self.0.update_identity(&mut hasher);
        hasher.finalize()
    }

    fn is_active(&self) -> bool {
        self.0
            .unit_cell()
            .iter()
            .any(|layer| layer.thickness_m() != 0.0)
    }
}

/// Borrowed arbitrary ordered transparent spectral stack. This preserves each
/// full Cauchy, Sellmeier or tabulated law and source without fitting it to the
/// compact scene material's Cauchy cell. Thicknesses and wavelengths use the
/// existing coating contracts. The stack does not allocate or own scene data.
#[derive(Clone, Copy, Debug)]
pub struct DielectricCoatingStack<'a> {
    layers: &'a [CoherentDielectricCoating],
}

impl<'a> DielectricCoatingStack<'a> {
    /// Layers are authored outside first; at most 64 physical films are admitted.
    /// An empty or entirely zero-thickness stack preserves the bare interface.
    pub fn try_new(layers: &'a [CoherentDielectricCoating]) -> Result<Self, DielectricError> {
        if layers.len() > MAX_COHERENT_STACK_LAYERS {
            return Err(DielectricError::InvalidCoating);
        }
        Ok(Self { layers })
    }

    /// Reflected/transmitted ray-power matrices. `reverse` traverses the authored
    /// stack in the opposite order; callers also supply reversed boundary indices
    /// and the corresponding Snell angle. Exact grazing of a nonzero stack is
    /// refused because incoming normal flux is zero. No absorption is inferred.
    pub fn fresnel_polarized(
        self,
        wavelength_nm: f64,
        reverse: bool,
        cosine: f64,
        eta_i: f64,
        eta_t: f64,
    ) -> Result<PolarizedFresnel, DielectricError> {
        super::validate_wavelength(wavelength_nm)?;
        let bare = fresnel_dielectric_polarized(cosine, eta_i, eta_t)?;
        if self.layers.iter().all(|layer| layer.thickness_m() == 0.0) {
            return Ok(bare);
        }
        if cosine == 0.0 {
            return Err(DielectricError::InvalidInterface);
        }
        let response = crate::conductor::lossless_stack_transmission(
            self.layers,
            wavelength_nm,
            eta_i,
            eta_t,
            cosine,
            reverse,
            exit_normal_wave(bare.scalar, cosine, eta_i, eta_t),
        )
        .map_err(|_| DielectricError::InvalidInterface)?;
        polarized_response(response, bare, self.layers.len())
    }

    /// Scalar power and native refraction classification from the same solution.
    pub fn fresnel(
        self,
        wavelength_nm: f64,
        reverse: bool,
        cosine: f64,
        eta_i: f64,
        eta_t: f64,
    ) -> Result<FresnelDielectric, DielectricError> {
        Ok(self
            .fresnel_polarized(wavelength_nm, reverse, cosine, eta_i, eta_t)?
            .scalar)
    }
}

/// Borrowed native interface owner at one wavelength and physical layer order.
/// Its scalar and polarized methods share one coating solution; sampling keeps
/// the existing Snell geometry, GGX Jacobians and radiance eta-squared factors.
#[derive(Clone, Copy)]
pub struct DielectricInterface<'a> {
    coating: Option<&'a DielectricCoating>,
    wavelength_nm: f64,
    reverse: bool,
}

impl<'a> DielectricInterface<'a> {
    pub(super) const fn new(
        coating: Option<&'a DielectricCoating>,
        wavelength_nm: f64,
        reverse: bool,
    ) -> Self {
        Self {
            coating,
            wavelength_nm,
            reverse,
        }
    }

    pub(super) fn is_active(self) -> bool {
        self.coating.is_some_and(DielectricCoating::is_active)
    }

    pub(super) fn has_contrast(self, index: f64) -> bool {
        self.coating.is_some_and(|coating| {
            coating.0.unit_cell().iter().any(|layer| {
                layer.thickness_m() != 0.0
                    && layer.ior().eval_unchecked(self.wavelength_nm) != index
            })
        })
    }

    /// Scalar normal-power flux, with exact native bare-interface behavior when
    /// no nonzero coating is present. Coating transmission is evaluated directly.
    pub fn fresnel(
        self,
        cosine: f64,
        eta_i: f64,
        eta_t: f64,
    ) -> Result<FresnelDielectric, DielectricError> {
        if !self.is_active() {
            return super::fresnel_dielectric(cosine, eta_i, eta_t);
        }
        Ok(self.fresnel_polarized(cosine, eta_i, eta_t)?.scalar)
    }

    /// Reflected and transmitted signed Mueller matrices in the native s/p frames.
    pub fn fresnel_polarized(
        self,
        cosine: f64,
        eta_i: f64,
        eta_t: f64,
    ) -> Result<PolarizedFresnel, DielectricError> {
        let bare = fresnel_dielectric_polarized(cosine, eta_i, eta_t)?;
        let Some(coating) = self.coating.filter(|coating| coating.is_active()) else {
            return Ok(bare);
        };
        if cosine == 0.0 {
            return Err(DielectricError::InvalidInterface);
        }
        let response = coating
            .0
            .lossless_transmission(
                self.wavelength_nm,
                eta_i,
                eta_t,
                cosine,
                self.reverse,
                exit_normal_wave(bare.scalar, cosine, eta_i, eta_t),
            )
            .map_err(|_| DielectricError::InvalidInterface)?;
        polarized_response(response, bare, coating.0.layer_count())
    }

    /// Native smooth reflection/transmission sample using this coherent response.
    pub fn sample_smooth(
        self,
        normal: Vec3,
        wo: Vec3,
        eta_i: f64,
        eta_t: f64,
        event_sample: f64,
    ) -> Result<SmoothDielectricSample, DielectricError> {
        super::sample_smooth_with_interface(normal, wo, eta_i, eta_t, event_sample, Some(self))
    }

    /// Smooth polarized sample with the same scalar proposal and Snell ray.
    /// Mueller throughput is divided by that proposal; transmission retains
    /// the native radiance eta-squared factor without changing power flux.
    pub fn sample_smooth_polarized(
        self,
        normal: Vec3,
        wo: Vec3,
        eta_i: f64,
        eta_t: f64,
        event_sample: f64,
    ) -> Result<PolarizedDielectricSample, DielectricError> {
        let sample = self.sample_smooth(normal, wo, eta_i, eta_t, event_sample)?;
        let n = normal.scale(1.0 / normal.dot(normal).sqrt());
        let w = wo.scale(1.0 / wo.dot(wo).sqrt());
        let response = self.fresnel_polarized(n.dot(w).clamp(0.0, 1.0), eta_i, eta_t)?;
        let (matrix, transport) = match sample.event {
            super::DielectricEvent::Reflection => (response.reflection, 1.0),
            super::DielectricEvent::Transmission => {
                (response.transmission, (eta_i / eta_t).powi(2))
            }
        };
        Ok(PolarizedDielectricSample {
            direction: sample.direction,
            event: sample.event,
            probability: sample.probability,
            radiance_weight: Mueller(
                matrix
                    .0
                    .map(|row| row.map(|value| (value / sample.probability) * transport)),
            ),
        })
    }

    /// Native GGX value and proposal density with this microfacet's coherent films.
    pub fn evaluate_rough(
        self,
        normal: Vec3,
        wo: Vec3,
        wi: Vec3,
        eta_i: f64,
        eta_t: f64,
        alpha: f64,
    ) -> Result<RoughDielectricEvaluation, DielectricError> {
        super::evaluate_rough_with_interface(normal, wo, wi, eta_i, eta_t, alpha, Some(self))
    }

    /// Native GGX sample. A nonzero coating between exactly equal bulk indices
    /// explicitly refuses the unsupported continuous/delta mixture; a smooth
    /// coated equal-index interface and an uncoated rough null remain supported.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_rough(
        self,
        normal: Vec3,
        wo: Vec3,
        eta_i: f64,
        eta_t: f64,
        alpha: f64,
        microfacet_u: f64,
        microfacet_v: f64,
        event_sample: f64,
    ) -> Result<Option<RoughDielectricSample>, DielectricError> {
        super::sample_rough_with_interface(
            normal,
            wo,
            eta_i,
            eta_t,
            alpha,
            microfacet_u,
            microfacet_v,
            event_sample,
            Some(self),
        )
    }
}

fn polarized_response(
    mut response: CoherentTransmission,
    mut bare: PolarizedFresnel,
    layer_count: usize,
) -> Result<PolarizedFresnel, DielectricError> {
    if response.bare_equivalent {
        return Ok(bare);
    }
    let tolerance = 512.0 * f64::EPSILON * (layer_count + 1) as f64;
    for pol in 0..2 {
        let reflected = response.reflected[pol];
        let transmitted = response.transmitted[pol];
        let closure = reflected[0] * reflected[0]
            + reflected[1] * reflected[1]
            + transmitted[0] * transmitted[0]
            + transmitted[1] * transmitted[1];
        if !closure.is_finite() || (closure - 1.0).abs() > tolerance {
            return Err(DielectricError::InvalidInterface);
        }
        // Only lossless films are admitted. Normalize bounded arithmetic
        // roundoff in amplitude space, keeping both channels and their phases.
        let scale = 1.0 / det::sqrt(closure);
        response.reflected[pol] = reflected.map(|component| component * scale);
        response.transmitted[pol] = transmitted.map(|component| component * scale);
    }
    let matrix = |amplitudes: [[f64; 2]; 2]| {
        let [s, p] = amplitudes;
        let s_power = s[0] * s[0] + s[1] * s[1];
        let p_power = p[0] * p[0] + p[1] * p[1];
        let mean = 0.5 * s_power + 0.5 * p_power;
        let difference = 0.5 * s_power - 0.5 * p_power;
        let real = s[0] * p[0] + s[1] * p[1];
        let imaginary = s[1] * p[0] - s[0] * p[1];
        Mueller([
            [mean, difference, 0.0, 0.0],
            [difference, mean, 0.0, 0.0],
            [0.0, 0.0, real, -imaginary],
            [0.0, 0.0, imaginary, real],
        ])
    };
    bare.reflection = matrix(response.reflected);
    bare.transmission = matrix(response.transmitted);
    bare.scalar.reflectance = bare.reflection.0[0][0].clamp(0.0, 1.0);
    bare.scalar.transmittance = bare.transmission.0[0][0].clamp(0.0, 1.0);
    Ok(bare)
}

fn exit_normal_wave(native: FresnelDielectric, cosine: f64, eta_i: f64, eta_t: f64) -> [f64; 2] {
    if native.total_internal_reflection {
        let ratio = eta_i / eta_t;
        let imaginary = eta_t
            * (ratio * ratio * (1.0 - cosine * cosine) - 1.0)
                .max(0.0)
                .sqrt();
        [0.0, imaginary]
    } else {
        [eta_t * native.transmitted_cosine, 0.0]
    }
}

#[cfg(test)]
mod tests;
