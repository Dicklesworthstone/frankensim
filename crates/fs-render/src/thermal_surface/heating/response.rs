//! Compiled surface power partitions from the original optical laws.
//!
//! Smooth isotropic interfaces depend on the incident s/p power fractions.
//! Relative U/V phase does not change their absorbed power. No microfacet
//! rejection, detector attenuation or camera exposure is interpreted as heat.

use super::SurfaceHeatingError;
use crate::conductor::{
    CoatingLayer, CoherentCoatingStack, ConductorOptics, MAX_COHERENT_STACK_LAYERS,
};
use crate::dielectric::DielectricCoatingStack;
use crate::polarization::Mueller;
use crate::polarization::dielectric::fresnel_dielectric_polarized;
use crate::thermal_surface::{GrayDiffuseThermalSurface, SpectralEmissivity};
use fs_blake3::{ContentHash, DomainHasher};

// The original dielectric coating owner admits this per-layer power roundoff.
// This is a floating-point admission allowance, not a certified error bound.
pub(super) const RESPONSE_ROUNDOFF: f64 =
    512.0 * f64::EPSILON * (MAX_COHERENT_STACK_LAYERS + 1) as f64;

/// Immutable spectral response for one already resolved incidence geometry.
///
/// Fractions multiply physical incident watts already intercepted by the
/// surface. A response carries no area, visibility or cosine projection factor.
/// Constitutive parameters remain fixed over the declared thermal domain of
/// its [`super::SurfaceIlluminationSource`]. Source identity is retained without
/// upgrading representative/authored optical data to measured-material data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfacePowerResponse<const N: usize> {
    identity: ContentHash,
    wavelengths_nm: [f64; N],
    reflected: [f64; N],
    transmitted: [f64; N],
    absorbed: [f64; N],
}

impl<const N: usize> SurfacePowerResponse<N> {
    /// Smooth opaque half-space, including the original attached coherent films.
    /// The optical owner obtains absorption directly from boundary power flux,
    /// preserving weak heat when reflected power rounds to the incident power.
    /// All film/substrate absorption belongs to the one assigned receiver; no
    /// depth or individual-film partition is inferred.
    pub fn from_conductor(
        wavelengths_nm: [f64; N],
        optics: &ConductorOptics,
        incident_eta: [f64; N],
        incidence_cosine: f64,
        s_power_fraction: [f64; N],
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        let mut hash = response_hash(b"opaque-conductor", optics.content_identity())?;
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        let mut reflected = [0.0; N];
        let mut absorbed = [0.0; N];
        for lane in 0..N {
            let response =
                optics.power_balance(wavelengths_nm[lane], incident_eta[lane], incidence_cosine)?;
            reflected[lane] = mix(response.reflected, s_power_fraction[lane]);
            absorbed[lane] = mix(response.absorbed, s_power_fraction[lane]);
        }
        Self::finish(hash, wavelengths_nm, reflected, [0.0; N], absorbed)
    }

    /// Arbitrary ordered coherent films on the original opaque substrate owner.
    /// As in the native stack API, this explicit ordered stack REPLACES any
    /// coating attached to the substrate; it is not appended to that coating.
    pub fn from_coated_conductor(
        wavelengths_nm: [f64; N],
        optics: &ConductorOptics,
        stack: CoherentCoatingStack<'_>,
        incident_eta: [f64; N],
        incidence_cosine: f64,
        s_power_fraction: [f64; N],
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        let mut hash = response_hash(b"opaque-coating-stack", optics.content_identity())?;
        hash.update(stack.content_identity().as_bytes());
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        let mut reflected = [0.0; N];
        let mut absorbed = [0.0; N];
        for lane in 0..N {
            let response = stack.power_balance(
                optics,
                wavelengths_nm[lane],
                incident_eta[lane],
                incidence_cosine,
            )?;
            reflected[lane] = mix(response.reflected, s_power_fraction[lane]);
            absorbed[lane] = mix(response.absorbed, s_power_fraction[lane]);
        }
        Self::finish(hash, wavelengths_nm, reflected, [0.0; N], absorbed)
    }

    /// Lossless smooth dielectric interface. Transmitted normal power remains
    /// an outgoing channel, with exactly zero surface heat. Bulk absorption,
    /// when present, belongs to the subsequent volume transport owner.
    pub fn from_dielectric(
        wavelengths_nm: [f64; N],
        incident_eta: [f64; N],
        transmitted_eta: [f64; N],
        incidence_cosine: f64,
        s_power_fraction: [f64; N],
        source_identity: ContentHash,
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        let mut hash = response_hash(b"lossless-dielectric", source_identity)?;
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        hash_values(&mut hash, transmitted_eta);
        let mut reflected = [0.0; N];
        let mut transmitted = [0.0; N];
        for lane in 0..N {
            // Even though Fresnel depends only on the supplied indices, retain
            // the native renderer's explicitly supported wavelength interval.
            if !(crate::spectral::LAMBDA_MIN..=crate::spectral::LAMBDA_MAX)
                .contains(&wavelengths_nm[lane])
            {
                return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
            }
            let response = fresnel_dielectric_polarized(
                incidence_cosine,
                incident_eta[lane],
                transmitted_eta[lane],
            )?;
            reflected[lane] = polarized_power(response.reflection, s_power_fraction[lane])?;
            transmitted[lane] = polarized_power(response.transmission, s_power_fraction[lane])?;
        }
        Self::finish(hash, wavelengths_nm, reflected, transmitted, [0.0; N])
    }

    /// Passive coherent films between two real-index half-spaces. Only film
    /// absorption heats the assigned receiver. Reflection and transmitted flux
    /// come from the native power-normalized Mueller solution, without the
    /// camera transport's radiance eta-squared factor.
    ///
    /// Positive coherent film work retains weak absorption even when the
    /// original R+T rounds to one. Lossless films contribute exactly zero;
    /// neither transmitted scene power nor floating-point residuals become heat.
    #[allow(clippy::too_many_arguments)]
    pub fn from_dielectric_stack(
        wavelengths_nm: [f64; N],
        layers: &[CoatingLayer<'_>],
        incident_eta: [f64; N],
        transmitted_eta: [f64; N],
        incidence_cosine: f64,
        reverse: bool,
        s_power_fraction: [f64; N],
        source_identity: ContentHash,
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &s_power_fraction)?;
        // Surrounding real indices are explicitly supplied at these packet
        // wavelengths. Every original film still validates its own support.
        let stack = DielectricCoatingStack::try_mixed_in_band(
            layers,
            (wavelengths_nm[0], wavelengths_nm[N - 1]),
        )?;
        let mut hash = response_hash(b"transmitting-coating-stack", source_identity)?;
        hash.update(
            CoherentCoatingStack::try_new(layers)?
                .content_identity()
                .as_bytes(),
        );
        hash.update(&[u8::from(reverse)]);
        hash_geometry(
            &mut hash,
            incidence_cosine,
            &incident_eta,
            &s_power_fraction,
        );
        hash_values(&mut hash, transmitted_eta);
        let mut reflected = [0.0; N];
        let mut transmitted = [0.0; N];
        let mut absorbed = [0.0; N];
        for lane in 0..N {
            let wavelength = wavelengths_nm[lane];
            let response = stack.layer_power_balance(
                wavelength,
                reverse,
                incidence_cosine,
                incident_eta[lane],
                transmitted_eta[lane],
            )?;
            reflected[lane] = polarized_power(response.fresnel.reflection, s_power_fraction[lane])?;
            transmitted[lane] =
                polarized_power(response.fresnel.transmission, s_power_fraction[lane])?;
            absorbed[lane] = mix(response.absorbed, s_power_fraction[lane]);
        }
        Self::finish(hash, wavelengths_nm, reflected, transmitted, absorbed)
    }

    /// Original explicit gray/diffuse law, valid at every supplied wavelength.
    /// Its constant absorptance is the same declared emissivity; no gray law is
    /// inferred from a conductor table or a finite-band spectral measurement.
    pub fn from_gray_diffuse(
        wavelengths_nm: [f64; N],
        surface: &GrayDiffuseThermalSurface,
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &[0.5; N])?;
        let hash = response_hash(b"gray-diffuse", surface.identity())?;
        Self::finish(
            hash,
            wavelengths_nm,
            [surface.reflectance(); N],
            [0.0; N],
            [surface.emissivity(); N],
        )
    }

    /// Original visible-band diffuse emissivity law, with opaque Kirchhoff
    /// absorption. This does not extend that source's spectral support or infer
    /// a bolometric cooling model.
    pub fn from_spectral_diffuse(
        wavelengths_nm: [f64; N],
        emissivity: &SpectralEmissivity,
        source_identity: ContentHash,
    ) -> Result<Self, SurfaceHeatingError> {
        validate_packet(&wavelengths_nm, &[0.5; N])?;
        let mut hash = response_hash(b"spectral-diffuse", source_identity)?;
        hash.update(emissivity.identity().as_bytes());
        let mut absorbed = [0.0; N];
        for lane in 0..N {
            absorbed[lane] = emissivity.at(wavelengths_nm[lane])?;
        }
        Self::finish(
            hash,
            wavelengths_nm,
            absorbed.map(|a| 1.0 - a),
            [0.0; N],
            absorbed,
        )
    }

    pub fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn wavelengths_nm(&self) -> &[f64; N] {
        &self.wavelengths_nm
    }
    pub fn reflected_fraction(&self) -> &[f64; N] {
        &self.reflected
    }
    pub fn transmitted_fraction(&self) -> &[f64; N] {
        &self.transmitted
    }
    pub fn absorbed_fraction(&self) -> &[f64; N] {
        &self.absorbed
    }

    pub(super) fn finish(
        mut hash: DomainHasher,
        wavelengths_nm: [f64; N],
        reflected: [f64; N],
        transmitted: [f64; N],
        absorbed: [f64; N],
    ) -> Result<Self, SurfaceHeatingError> {
        for lane in 0..N {
            let values = [reflected[lane], transmitted[lane], absorbed[lane]];
            if !values
                .into_iter()
                .all(|v| v.is_finite() && (0.0..=1.0).contains(&v))
                || ((1.0 - reflected[lane]) - transmitted[lane] - absorbed[lane]).abs()
                    > RESPONSE_ROUNDOFF
            {
                return Err(invalid(
                    "finite passive surface power fractions must close to one",
                ));
            }
        }
        hash.update(&(N as u64).to_le_bytes());
        hash_values(&mut hash, wavelengths_nm);
        hash_values(&mut hash, reflected);
        hash_values(&mut hash, transmitted);
        hash_values(&mut hash, absorbed);
        Ok(Self {
            identity: hash.finalize(),
            wavelengths_nm,
            reflected,
            transmitted,
            absorbed,
        })
    }
}

pub(super) fn validate_packet<const N: usize>(
    wavelengths: &[f64; N],
    s_fraction: &[f64; N],
) -> Result<(), SurfaceHeatingError> {
    if N == 0
        || wavelengths.iter().any(|v| !v.is_finite() || *v <= 0.0)
        || wavelengths.windows(2).any(|pair| pair[0] >= pair[1])
        || s_fraction
            .iter()
            .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
    {
        return Err(invalid(
            "ordered positive wavelengths and physical s/p power fractions required",
        ));
    }
    Ok(())
}

pub(super) fn response_hash(
    kind: &[u8],
    source: ContentHash,
) -> Result<DomainHasher, SurfaceHeatingError> {
    if source == ContentHash([0; 32]) {
        return Err(invalid("explicit surface optical source identity required"));
    }
    let mut hash = DomainHasher::new("org.frankensim.render.surface-power-response.v1");
    hash.update(kind);
    hash.update(source.as_bytes());
    Ok(hash)
}

pub(super) fn hash_geometry<const N: usize>(
    hash: &mut DomainHasher,
    cosine: f64,
    incident_eta: &[f64; N],
    s_fraction: &[f64; N],
) {
    hash.update(&cosine.to_bits().to_le_bytes());
    hash_values(hash, *incident_eta);
    hash_values(hash, *s_fraction);
}

pub(super) fn hash_values<const N: usize>(hash: &mut DomainHasher, values: [f64; N]) {
    for value in values {
        hash.update(&value.to_bits().to_le_bytes());
    }
}

pub(super) fn mix(values: [f64; 2], s_fraction: f64) -> f64 {
    s_fraction * values[0] + (1.0 - s_fraction) * values[1]
}

pub(super) fn polarized_power(
    matrix: Mueller,
    s_fraction: f64,
) -> Result<f64, SurfaceHeatingError> {
    let value = matrix.0[0][0] + (2.0 * s_fraction - 1.0) * matrix.0[0][1];
    if !value.is_finite() || !(-RESPONSE_ROUNDOFF..=1.0 + RESPONSE_ROUNDOFF).contains(&value) {
        return Err(invalid("nonpassive polarized surface power"));
    }
    Ok(value.clamp(0.0, 1.0))
}

fn invalid(message: &'static str) -> SurfaceHeatingError {
    SurfaceHeatingError::Invalid(message)
}
