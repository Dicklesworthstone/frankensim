//! Source-bound, passive coherent coatings on transmitting dielectric bodies.
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
    CoatingLayer, CoherentCoatingStack, CoherentDielectricCoating, CoherentPeriodicCoating, CoherentTransmission,
    MAX_COHERENT_STACK_LAYERS,
};
use crate::polarization::Mueller;
use crate::polarization::dielectric::{
    PolarizedDielectricSample, PolarizedFresnel, fresnel_dielectric_polarized,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::Vec3;
use fs_math::det;

/// Compact passive coating on a dielectric surface. Each original layer owns
/// its real phase index, nonnegative extinction and physical thickness.
/// Absorbing films retain independent reflected and transmitted fields;
/// survival is R+T, never an invented lossless partition. Bulk Beer-Lambert
/// absorption remains independent and is applied only along bulk ray segments.
/// The source-bound Cauchy unit cell and 64-layer work bound are unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricCoating(CoherentPeriodicCoating);

impl DielectricCoating {
    /// Admit a fully declared passive periodic coating. The original layer
    /// constructors enforce finite indices, nonnegative loss and bounded work.
    pub fn try_new(coating: CoherentPeriodicCoating) -> Result<Self, DielectricError> {
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

/// Borrowed arbitrary ordered passive spectral stack. This preserves each
/// full Cauchy, Sellmeier or tabulated law, spectral extinction and source
/// without fitting it to the compact scene material's Cauchy cell. Thicknesses
/// and wavelengths use the existing coating contracts. No allocation or
/// per-query copy of the original optical tables is needed.
#[derive(Clone, Copy, Debug)]
pub struct DielectricCoatingStack<'a> {
    layers: StackLayers<'a>,
}

#[derive(Clone, Copy, Debug)]
enum StackLayers<'a> {
    Transparent(&'a [CoherentDielectricCoating]),
    Mixed(CoherentCoatingStack<'a>),
}

/// Native coherent reflection/transmission and positive absorption in finite
/// dielectric coating films. All fractions refer to incident ray power, not
/// camera radiance. The real-index exit medium receives transmitted radiation;
/// its subsequent bulk absorption is a separate physical transport operation.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricLayerPowerBalance {
    /// Original signed, power-normalized Mueller fields and native Snell policy.
    pub fresnel: PolarizedFresnel,
    /// Total finite-film absorption, ordered s, p. This is a positive sum of
    /// physical dissipative work, independent of the rounded R+T residual.
    pub absorbed: [f64; 2],
    layers: [[f64; 2]; MAX_COHERENT_STACK_LAYERS],
    layer_count: usize,
}

impl DielectricLayerPowerBalance {
    /// Per-film s/p absorption in the original authored outside-to-inside order.
    /// Reverse incidence changes the coherent solution, never this association
    /// with physical layer identities. Zero-thickness/lossless layers stay zero.
    #[must_use]
    pub fn layer_absorbed(&self) -> &[[f64; 2]] {
        &self.layers[..self.layer_count]
    }
}

impl<'a> DielectricCoatingStack<'a> {
    /// Layers are authored outside first; at most 64 physical films are admitted.
    /// An empty or entirely zero-thickness stack preserves the bare interface.
    pub fn try_new(layers: &'a [CoherentDielectricCoating]) -> Result<Self, DielectricError> {
        if layers.len() > MAX_COHERENT_STACK_LAYERS {
            return Err(DielectricError::InvalidCoating);
        }
        Ok(Self { layers: StackLayers::Transparent(layers) })
    }

    /// Borrow arbitrarily ordered transparent and absorbing spectral films.
    /// Their authored order runs outside first. Negative extinction is already
    /// refused by the original film constructors. The same 64-layer bound
    /// applies, including empty and zero-thickness stacks.
    pub fn try_mixed(layers: &'a [CoatingLayer<'a>]) -> Result<Self, DielectricError> {
        let stack = CoherentCoatingStack::try_new(layers)
            .map_err(|_| DielectricError::InvalidCoating)?;
        Ok(Self { layers: StackLayers::Mixed(stack) })
    }

    /// Number of authored physical layers, including zero-thickness films.
    #[must_use]
    pub fn layer_count(self) -> usize {
        match self.layers {
            StackLayers::Transparent(layers) => layers.len(),
            StackLayers::Mixed(stack) => stack.layers().len(),
        }
    }

    fn is_active(self) -> bool {
        match self.layers {
            StackLayers::Transparent(layers) => layers.iter().any(|film| film.thickness_m() != 0.0),
            StackLayers::Mixed(stack) => stack.layers().iter().any(|film| match film {
                CoatingLayer::Dielectric(film) => film.thickness_m() != 0.0,
                CoatingLayer::Absorbing(film) => film.thickness_m() != 0.0,
            }),
        }
    }

    /// Reflected/transmitted ray-power matrices. `reverse` traverses the authored
    /// stack in the opposite order; callers also supply reversed boundary indices
    /// and the corresponding Snell angle. Exact grazing of a nonzero stack is
    /// refused because incoming normal flux is zero. Absorption reduces R+T;
    /// no lossless renormalization is applied to a lossy physical solution.
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
        if !self.is_active() {
            return Ok(bare);
        }
        if cosine == 0.0 {
            return Err(DielectricError::InvalidInterface);
        }
        let exit = exit_normal_wave(bare.scalar, cosine, eta_i, eta_t);
        let response = match self.layers {
            StackLayers::Transparent(layers) => crate::conductor::lossless_stack_transmission(
                layers, wavelength_nm, eta_i, eta_t, cosine, reverse, exit,
            ),
            StackLayers::Mixed(stack) => stack.transmission(
                wavelength_nm, eta_i, eta_t, cosine, reverse, exit,
            ),
        }.map_err(|_| DielectricError::InvalidInterface)?;
        polarized_response(response, bare, self.layer_count())
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

    /// Resolve physical coating heat without subtracting reflected and
    /// transmitted power. The existing factored-wave solve supplies both its
    /// native R/T response and bounded interface histories; analytical positive
    /// even/odd film work retains weak losses below the R+T rounding scale.
    ///
    /// `reverse` reverses optical traversal and expects the caller's incident/
    /// exit indices and corresponding angle. Returned layers always retain
    /// authored order. This allocation-free query keeps the original 64-layer,
    /// spectral-support, coherence and exact-grazing admission rules.
    pub fn layer_power_balance(
        self,
        wavelength_nm: f64,
        reverse: bool,
        cosine: f64,
        eta_i: f64,
        eta_t: f64,
    ) -> Result<DielectricLayerPowerBalance, DielectricError> {
        super::validate_wavelength(wavelength_nm)?;
        let bare = fresnel_dielectric_polarized(cosine, eta_i, eta_t)?;
        let mut result = DielectricLayerPowerBalance {
            fresnel: bare,
            absorbed: [0.0; 2],
            layers: [[0.0; 2]; MAX_COHERENT_STACK_LAYERS],
            layer_count: self.layer_count(),
        };
        if !self.is_active() {
            return Ok(result);
        }
        if cosine == 0.0 {
            return Err(DielectricError::InvalidInterface);
        }
        let exit = exit_normal_wave(bare.scalar, cosine, eta_i, eta_t);
        let response = match self.layers {
            StackLayers::Transparent(layers) => crate::conductor::lossless_stack_transmission(
                layers, wavelength_nm, eta_i, eta_t, cosine, reverse, exit,
            ).map_err(|_| DielectricError::InvalidInterface)?,
            StackLayers::Mixed(stack) => {
                let work = stack.transmission_with_absorption(
                    wavelength_nm, eta_i, eta_t, cosine, reverse, exit,
                ).map_err(|_| DielectricError::InvalidInterface)?;
                result.absorbed = work.absorbed;
                result.layers = work.layers;
                work.response
            }
        };
        result.fresnel = polarized_response(response, bare, self.layer_count())?;
        Ok(result)
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
                    && (layer.extinction_k() != 0.0
                        || layer.ior().eval_unchecked(self.wavelength_nm) != index)
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
            .transmission(
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

    /// Native smooth sample using this coherent response. The proposal is
    /// conditional on survival; its weight retains R+T so absorbing coatings
    /// remove power once without a separate absorption roulette event.
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
        if !closure.is_finite() || closure > 1.0 + tolerance
            || (response.lossless && (closure - 1.0).abs() > tolerance) {
            return Err(DielectricError::InvalidInterface);
        }
        // Preserve the established lossless arithmetic. A passive absorbing
        // stack retains its actual survival deficit; only an above-one excess
        // within the declared floating-point allowance may be removed.
        if response.lossless || closure > 1.0 {
            let scale = 1.0 / det::sqrt(closure);
            response.reflected[pol] = reflected.map(|component| component * scale);
            response.transmitted[pol] = transmitted.map(|component| component * scale);
        }
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

#[cfg(test)]
mod absorption_tests;
