//! Spectral wave response from original metal and dielectric material sources.
//!
//! The source uses absolute vacuum-referenced n and extinction k. With this
//! owner's exp(+i omega t) convention N=n-i k and epsilon=N^2/(mu c0^2).
//! Permeability is explicitly sourced in the SAME complete material snapshot;
//! no nonmagnetic default, DC-conductivity fit or second loss term is inserted.
//! Metal n/k or dielectric phase index and intensity absorption retain their
//! original spectral parameterization. All interpolation stays in its band;
//! RGB transmittance cannot supply a monochromatic Maxwell constitutive law.
//!
//! Wave and absorption conventions: S. J. Byrnes, "Multilayer optical
//! calculations", sections 2--4, <https://arxiv.org/abs/1603.02720> (opposite
//! time convention). The general dispersion relation retains sourced mu.

use super::{ReflectionError, positive};
use fs_material::{
    field_storage::linear_magnetic_permeability,
    state_point::{
        ResolvedMaterialStatePoint, VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
        VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
        VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY, VISIBLE_DIELECTRIC_IOR_PROPERTIES,
        VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
        VISIBLE_DIELECTRIC_WAVELENGTHS_NM, VisibleConductorStatePoint, VisibleDielectricAbsorption,
        VisibleOpticalStatePoint,
    },
};
use fs_math::{c64::C64, det};

mod dispersion;

/// Exact SI vacuum light speed, used only to interpret vacuum wavelengths.
const VACUUM_SPEED_M_S: f64 = 299_792_458.0;
const FREQUENCY_WAVELENGTH_NM: f64 = core::f64::consts::TAU * VACUUM_SPEED_M_S * 1.0e9;

/// Original spectral complex-index law with its explicitly sourced permeability.
///
/// Adopt a complete material bundle with
/// [`VisibleConductorStatePoint::try_from_material`]; that bundle must contain
/// the optical knots and absolute scalar magnetic permeability. The original
/// spectrum supplies all electric loss. Ohmic loss must not be added again.
/// This homogeneous isotropic monochromatic law does not supply a causal
/// time-domain realization or extrapolate the source band.
#[derive(Clone, Copy, Debug)]
pub struct OpticalConductorResponse<'a> {
    source: &'a VisibleConductorStatePoint,
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> OpticalConductorResponse<'a> {
    /// Validate the original permeability and retain the complete source owner.
    pub fn try_new(source: &'a VisibleConductorStatePoint) -> Result<Self, ReflectionError> {
        let permeability_h_m = linear_magnetic_permeability(source.resolved())?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        Ok(Self {
            source,
            permeability_h_m,
            index_to_sqrt_permittivity,
        })
    }

    pub const fn source(&self) -> &'a VisibleConductorStatePoint {
        self.source
    }

    pub const fn material(&self) -> &'a ResolvedMaterialStatePoint {
        self.source.resolved()
    }

    pub const fn permeability_h_m(&self) -> f64 {
        self.permeability_h_m
    }

    /// Original n-i k at one angular frequency [rad/s], with no extrapolation.
    /// At most eight source intervals are visited and no allocation is made.
    pub fn complex_index(&self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let omega = positive(omega_rad_s)?;
        let samples = self.source.samples();
        let first = samples[0].wavelength_nm;
        let last = samples[samples.len() - 1].wavelength_nm;
        if omega < FREQUENCY_WAVELENGTH_NM / last || omega > FREQUENCY_WAVELENGTH_NM / first {
            return Err(ReflectionError::Invalid(
                "optical wave frequency is outside the original n/k source band",
            ));
        }
        // Frequency admission precedes this rounding guard: it cannot extend
        // the source band, even by one representable frequency step.
        let wavelength = (FREQUENCY_WAVELENGTH_NM / omega).clamp(first, last);
        for pair in samples.windows(2) {
            let [left, right] = [pair[0], pair[1]];
            if wavelength <= right.wavelength_nm {
                let t =
                    (wavelength - left.wavelength_nm) / (right.wavelength_nm - left.wavelength_nm);
                let n = (1.0 - t) * left.eta + t * right.eta;
                let k = (1.0 - t) * left.k + t * right.k;
                return Ok(C64::new(n, -k));
            }
        }
        Err(ReflectionError::Invalid(
            "optical source interval is missing",
        ))
    }

    /// Absolute complex epsilon [F/m], including the original source loss.
    /// Negative real epsilon is valid for a passive metal; passivity requires
    /// the imaginary part to be nonpositive under exp(+i omega t).
    pub fn permittivity_f_m(&self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let index = self.complex_index(omega_rad_s)?;
        dispersion::permittivity(index, self.index_to_sqrt_permittivity).ok_or(
            ReflectionError::Invalid(
                "optical permittivity or its original loss is unrepresentable",
            ),
        )
    }
}

/// Original visible optical family with explicit same-state permeability.
///
/// Metals retain their n/k table. Cauchy, tabulated and Sellmeier dielectrics
/// retain their phase-index law and typed spectral POWER absorption alpha.
/// For vacuum wavelength lambda, k=alpha*lambda/(4*pi); alpha is not a field
/// attenuation coefficient, a scattering loss, or external slab transmission.
/// RGB-only absorption refuses instead of inventing a wavelength spectrum.
///
/// Extend the original family with permeability through
/// [`VisibleOpticalStatePoint::with_additional_requirements`]. The source must
/// retain absolute scalar mu in that same complete snapshot. No second loss
/// term or nonmagnetic default is introduced. This homogeneous isotropic
/// monochromatic response does not imply a causal time-domain realization.
#[derive(Clone, Copy, Debug)]
pub struct OpticalMaterialResponse<'a> {
    source: &'a VisibleOpticalStatePoint,
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> OpticalMaterialResponse<'a> {
    pub fn try_new(source: &'a VisibleOpticalStatePoint) -> Result<Self, ReflectionError> {
        let permeability_h_m = linear_magnetic_permeability(source.resolved())?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        match source {
            VisibleOpticalStatePoint::Conductor(_) => {}
            VisibleOpticalStatePoint::Dielectric(value) => {
                admit_static_properties(
                    source.resolved(),
                    &[
                        VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY,
                        VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
                        VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
                    ],
                )?;
                admit_spectral_absorption(source.resolved(), value.absorption())?;
            }
            VisibleOpticalStatePoint::TabulatedDielectric(_) => {
                admit_static_properties(source.resolved(), &VISIBLE_DIELECTRIC_IOR_PROPERTIES)?;
                admit_static_properties(
                    source.resolved(),
                    &VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
                )?;
            }
            VisibleOpticalStatePoint::DielectricSellmeier(value) => {
                admit_static_properties(
                    source.resolved(),
                    &VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES,
                )?;
                admit_static_properties(
                    source.resolved(),
                    &VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
                )?;
                admit_spectral_absorption(source.resolved(), value.absorption())?;
                let (strengths, resonances) = value.sellmeier_coefficients_si();
                if !dispersion::admits_sellmeier_band(strengths, resonances) {
                    return Err(ReflectionError::Invalid(
                        "Sellmeier wave source must have positive finite index and no active pole in its complete visible band",
                    ));
                }
            }
        }
        let response = Self {
            source,
            permeability_h_m,
            index_to_sqrt_permittivity,
        };
        // Cauchy coefficients are nonnegative and tabulated knots are admitted
        // by their source. Endpoint checks additionally catch overflow over the
        // complete band; Sellmeier also excludes poles between those endpoints.
        for wavelength in dispersion::VISIBLE_BAND_NM {
            response.complex_index(FREQUENCY_WAVELENGTH_NM / wavelength)?;
        }
        Ok(response)
    }

    pub const fn source(&self) -> &'a VisibleOpticalStatePoint {
        self.source
    }
    pub const fn material(&self) -> &'a ResolvedMaterialStatePoint {
        self.source.resolved()
    }
    pub const fn permeability_h_m(&self) -> f64 {
        self.permeability_h_m
    }

    /// Absolute n-i*k at the requested angular frequency, without extrapolation.
    /// Dielectric n and alpha are interpolated before converting alpha to k;
    /// interpolating k instead would change the original absorption spectrum.
    pub fn complex_index(&self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let wavelength =
            dispersion::visible_wavelength_nm(omega_rad_s).ok_or(ReflectionError::Invalid(
                "optical wave frequency is outside the original visible source band",
            ))?;
        let (n, alpha) = match self.source {
            VisibleOpticalStatePoint::Conductor(source) => {
                return OpticalConductorResponse {
                    source,
                    permeability_h_m: self.permeability_h_m,
                    index_to_sqrt_permittivity: self.index_to_sqrt_permittivity,
                }
                .complex_index(omega_rad_s);
            }
            VisibleOpticalStatePoint::Dielectric(source) => (
                dispersion::cauchy_index(source.cauchy_coefficients_si(), wavelength).ok_or(
                    ReflectionError::Invalid("Cauchy phase index is unrepresentable"),
                )?,
                spectral_absorption(source.absorption(), wavelength)?,
            ),
            VisibleOpticalStatePoint::TabulatedDielectric(source) => {
                let samples = source.samples();
                (
                    interpolate(samples.map(|sample| sample.phase_index), wavelength)?,
                    interpolate(samples.map(|sample| sample.absorption_per_m), wavelength)?,
                )
            }
            VisibleOpticalStatePoint::DielectricSellmeier(source) => {
                let (strengths, resonances) = source.sellmeier_coefficients_si();
                (
                    dispersion::sellmeier_index(strengths, resonances, wavelength).ok_or(
                        ReflectionError::Invalid("Sellmeier phase index is unrepresentable"),
                    )?,
                    spectral_absorption(source.absorption(), wavelength)?,
                )
            }
        };
        dispersion::attenuating_index(n, alpha, wavelength).ok_or(ReflectionError::Invalid(
            "optical phase index or original intensity absorption is unrepresentable",
        ))
    }

    pub fn permittivity_f_m(&self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        dispersion::permittivity(
            self.complex_index(omega_rad_s)?,
            self.index_to_sqrt_permittivity,
        )
        .ok_or(ReflectionError::Invalid(
            "optical permittivity or original absorption is unrepresentable",
        ))
    }
}

fn admit_static_properties(
    material: &ResolvedMaterialStatePoint,
    names: &[&str],
) -> Result<(), ReflectionError> {
    for name in names {
        if !material
            .property(name)
            .is_some_and(|property| property.is_static_at_fixed_external_conditions())
        {
            return Err(ReflectionError::Invalid(
                "fixed-wavelength optical sources require fixed external conditions without unowned frequency, field or deformation coordinates",
            ));
        }
    }
    Ok(())
}

fn admit_spectral_absorption(
    material: &ResolvedMaterialStatePoint,
    source: VisibleDielectricAbsorption,
) -> Result<(), ReflectionError> {
    if !matches!(source, VisibleDielectricAbsorption::Spectral { .. }) {
        return Err(ReflectionError::Invalid(
            "monochromatic Maxwell absorption requires the original spectral alpha; RGB transmittance is insufficient",
        ));
    }
    admit_static_properties(material, &VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES)
}

fn spectral_absorption(
    source: VisibleDielectricAbsorption,
    wavelength_nm: f64,
) -> Result<f64, ReflectionError> {
    match source {
        VisibleDielectricAbsorption::Spectral { absorption_per_m } => {
            interpolate(absorption_per_m, wavelength_nm)
        }
        VisibleDielectricAbsorption::ReferenceRgb { .. } => Err(ReflectionError::Invalid(
            "RGB transmittance cannot define a unique monochromatic Maxwell absorption law",
        )),
    }
}

fn interpolate(values: [f64; 9], wavelength_nm: f64) -> Result<f64, ReflectionError> {
    for (i, pair) in VISIBLE_DIELECTRIC_WAVELENGTHS_NM.windows(2).enumerate() {
        if wavelength_nm <= pair[1] {
            let t = (wavelength_nm - pair[0]) / (pair[1] - pair[0]);
            return Ok((1.0 - t) * values[i] + t * values[i + 1]);
        }
    }
    Err(ReflectionError::Invalid(
        "optical source interval is missing",
    ))
}
