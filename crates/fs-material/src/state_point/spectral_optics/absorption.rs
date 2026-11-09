//! The existing canonical dielectric loss spectrum without a phase-index role.

use super::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    SPECTRAL_OPTICAL_SAMPLE_COUNT, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
    ScalarAdmissibility, ScalarPropertyRequirement, SpectralDielectricAbsorptionSample,
    WAVELENGTH_QUANTITY, resolve_material_state_point, source_value, source_wavelengths,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};

/// Independently supported intensity absorption at one complete material state.
///
/// This projects the SAME nine wavelength and nine alpha properties already
/// consumed by tabulated and Sellmeier dielectrics. It requires no phase-index
/// claims, so a thermal absorption source can be refreshed without extrapolating
/// a separately held reference-index fit. Every additional resolved requirement
/// and source receipt is retained and refreshed, including permeability when
/// needed by a Maxwell consumer. No absorption(T), scattering or phase law is
/// inferred; state dependence comes exclusively from the selected source claims.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralDielectricAbsorptionStatePoint {
    resolved: ResolvedMaterialStatePoint,
    samples: [SpectralDielectricAbsorptionSample; SPECTRAL_OPTICAL_SAMPLE_COUNT],
}

impl SpectralDielectricAbsorptionStatePoint {
    /// Adopt canonical scalar quantities and ordered positive SI wavelengths.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        let wavelengths = source_wavelengths(&resolved)?;
        let mut samples = [SpectralDielectricAbsorptionSample {
            wavelength_m: 0.0,
            absorption_per_m: 0.0,
        }; SPECTRAL_OPTICAL_SAMPLE_COUNT];
        for (index, sample) in samples.iter_mut().enumerate() {
            *sample = SpectralDielectricAbsorptionSample {
                wavelength_m: wavelengths[index],
                absorption_per_m: source_value(
                    &resolved,
                    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                    ScalarAdmissibility::NonNegative,
                )?,
            };
        }
        Ok(Self { resolved, samples })
    }

    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Identity includes all original selected properties and state coordinates.
    pub fn identity(&self) -> ContentHash {
        let mut hash =
            DomainHasher::new("org.frankensim.material.spectral-dielectric-absorption.v1");
        hash.update(self.resolved.identity().as_bytes());
        hash.finalize()
    }

    pub const fn absorption_samples(
        &self,
    ) -> &[SpectralDielectricAbsorptionSample; SPECTRAL_OPTICAL_SAMPLE_COUNT] {
        &self.samples
    }

    pub fn wavelength_bounds_m(&self) -> (f64, f64) {
        (self.samples[0].wavelength_m, self.samples[8].wavelength_m)
    }

    /// Original piecewise-linear alpha [power Np/m], without allocating or
    /// evaluating a refractive index. Transparent spectra retain their support.
    pub fn absorption_per_m_at(&self, wavelength_m: f64) -> Result<f64, MaterialStatePointError> {
        interpolate(wavelength_m, |index| self.samples[index])
    }

    /// Refresh the exact original selected claims, including source wavelengths.
    /// No reference phase properties or synthesized coefficients are inserted.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(self.resolved.at_query_point(card, point)?)
    }
}

/// Resolve the 18 original wavelength/alpha roles, independently of refraction.
pub fn resolve_spectral_dielectric_absorption_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<SpectralDielectricAbsorptionStatePoint, MaterialStatePointError> {
    let mut requirements = Vec::with_capacity(2 * SPECTRAL_OPTICAL_SAMPLE_COUNT);
    for (names, quantity, admissibility) in [
        (
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
            WAVELENGTH_QUANTITY,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            ScalarAdmissibility::NonNegative,
        ),
    ] {
        for name in names {
            requirements.push(ScalarPropertyRequirement::try_with_quantity(
                name,
                quantity,
                admissibility,
            )?);
        }
    }
    SpectralDielectricAbsorptionStatePoint::try_from_material(resolve_material_state_point(
        card,
        point,
        &requirements,
        selection,
    )?)
}

// One original interpolation law shared with BOTH existing dielectric sources.
// The closure borrows their native layouts; no extra table or resampling occurs.
pub(super) fn interpolate(
    wavelength_m: f64,
    sample: impl Fn(usize) -> SpectralDielectricAbsorptionSample,
) -> Result<f64, MaterialStatePointError> {
    if !wavelength_m.is_finite()
        || wavelength_m < sample(0).wavelength_m
        || wavelength_m > sample(SPECTRAL_OPTICAL_SAMPLE_COUNT - 1).wavelength_m
    {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "wavelength outside the original intensity-absorption spectrum",
        });
    }
    let (mut low, mut high) = (0, SPECTRAL_OPTICAL_SAMPLE_COUNT);
    while low < high {
        let middle = low + (high - low) / 2;
        if sample(middle).wavelength_m < wavelength_m {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    let right = sample(low);
    if right.wavelength_m == wavelength_m {
        return Ok(right.absorption_per_m);
    }
    let left = sample(low - 1);
    let t = (wavelength_m - left.wavelength_m) / (right.wavelength_m - left.wavelength_m);
    let [a, b] = [left.absorption_per_m, right.absorption_per_m];
    let value = if t <= 0.5 {
        a + t * (b - a)
    } else {
        b + (1.0 - t) * (a - b)
    };
    if !ScalarAdmissibility::NonNegative.admits(value) {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "unrepresentable source intensity-absorption interpolation",
        });
    }
    Ok(value)
}
