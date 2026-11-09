//! Original NASA CHARMS silica refraction in homogeneous Maxwell waves.
//!
//! The shared material owner retains the measured Corning 7980 Table-3 fit,
//! its absolute vacuum indices and its 30--300 K, 400--2600 nm support:
//! <https://ntrs.nasa.gov/citations/20070018851>.
//! The caller selects zero loss or a separate original current-state alpha
//! spectrum. CHARMS itself supplies no absorption, causal reservoir or stress law.

use super::{
    FREQUENCY_WAVELENGTH_NM, ReflectionError, VACUUM_SPEED_M_S,
    absorption::{self, SourcedAbsorption}, dispersion, positive,
};
use fs_material::{
    field_storage::linear_magnetic_permeability,
    silica_refraction::{Corning7980Silica, SilicaRefractionSample, WAVELENGTH_RANGE_NM},
    state_point::{ResolvedMaterialStatePoint, SpectralDielectricAbsorptionStatePoint},
};
use fs_math::{c64::C64, det};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};

/// One original temperature-resolved Corning 7980 law and its sourced mu.
///
/// The caller explicitly selects applicability of the named silica source;
/// no material-name match or generic fused-silica assumption is made. The
/// complete scalar source owns temperature and absolute SI permeability.
/// Pressure, strain, composition and frequency coordinates refuse because
/// this optical fit does not supply those constitutive dependences.
///
/// Preparation freezes the original six temperature polynomials once. All
/// wavelength queries allocate nothing and reuse that same material owner.
/// The complete n-i*k supplies epsilon=(n-i*k)²/(mu*c0²); phase/group refraction
/// always comes from the original analytic source. Optional current-state alpha
/// is independently sourced at that SAME temperature and permeability. This
/// homogeneous isotropic model does not add expansion, scattering or spatial GRIN.
#[derive(Clone, Copy, Debug)]
pub struct SilicaOpticalResponse<'a> {
    source: &'a ResolvedMaterialStatePoint,
    absorption: Option<SourcedAbsorption<'a>>,
    refractive_source: Corning7980Silica,
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> SilicaOpticalResponse<'a> {
    /// Explicitly choose zero absorption for the original CHARMS index law.
    ///
    /// `temperature_axis` must identify the sole static absolute-temperature
    /// coordinate. The existing untyped/dimensional `T` convention is kelvin;
    /// any other axis name requires absolute-temperature semantics. The
    /// source's temperature and wavelength support are never extrapolated.
    pub fn try_lossless(
        source: &'a ResolvedMaterialStatePoint,
        temperature_axis: &str,
    ) -> Result<Self, ReflectionError> {
        let mut temperature = None;
        for (axis, value) in source.query_point() {
            let quantity = source.axis_quantities().get(axis).copied();
            let is_absolute = quantity
                .and_then(QuantitySpec::semantic_type)
                .is_some_and(|ty| {
                    ty.kind() == QuantityKind::AbsoluteTemperature && ty.form() == ValueForm::Static
                })
                || axis == "T"
                    && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
            if !is_absolute || axis != temperature_axis || temperature.replace(*value).is_some() {
                return Err(ReflectionError::Invalid(
                    "CHARMS silica optics requires only one owned absolute-temperature coordinate",
                ));
            }
        }
        let temperature_k = temperature.ok_or(ReflectionError::Invalid(
            "CHARMS silica optics needs original absolute temperature in kelvin",
        ))?;
        let refractive_source = Corning7980Silica::at_temperature(temperature_k).map_err(|_| {
            ReflectionError::Invalid("silica temperature is outside the original CHARMS fit")
        })?;
        let permeability_h_m = linear_magnetic_permeability(source)?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        let response = Self {
            source,
            absorption: None,
            refractive_source,
            permeability_h_m,
            index_to_sqrt_permittivity,
        };
        for wavelength in WAVELENGTH_RANGE_NM {
            response.permittivity_f_m(FREQUENCY_WAVELENGTH_NM / wavelength)?;
        }
        Ok(response)
    }

    /// Compose original CHARMS phase with independently sourced current loss.
    ///
    /// The original canonical alpha snapshot also owns the actual temperature
    /// and absolute scalar permeability. Its wavelength support intersects the
    /// unchanged CHARMS band, including when every alpha sample is zero. The
    /// positive finite interpolation-gap admission is in vacuum metres.
    pub fn try_with_current_absorption(
        source: &'a SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        let absorption =
            SourcedAbsorption::try_new(source, maximum_knot_span_m, WAVELENGTH_RANGE_NM)?;
        let mut response = Self::try_lossless(source.resolved(), temperature_axis)?;
        response.absorption = Some(absorption);
        for omega in absorption.frequency_range_rad_s() {
            response.permittivity_f_m(omega)?;
        }
        Ok(response)
    }

    /// Original current-state intensity loss, or the explicit lossless policy.
    pub fn current_absorption(self) -> Option<&'a SpectralDielectricAbsorptionStatePoint> {
        self.absorption.map(SourcedAbsorption::source)
    }

    /// Original complete scalar snapshot, including its source-use identity.
    pub const fn source(self) -> &'a ResolvedMaterialStatePoint {
        self.source
    }
    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.source
    }
    /// The original source law at the retained temperature, without refitting.
    pub const fn refractive_source(self) -> Corning7980Silica {
        self.refractive_source
    }
    pub const fn temperature_k(self) -> f64 {
        self.refractive_source.temperature_k()
    }
    pub const fn permeability_h_m(self) -> f64 {
        self.permeability_h_m
    }

    /// Original absolute phase/group indices and thermal derivatives at the
    /// requested vacuum frequency. The shared source defines all derivatives.
    pub fn refractive_sample(
        self,
        omega_rad_s: f64,
    ) -> Result<SilicaRefractionSample, ReflectionError> {
        let wavelength = dispersion::wavelength_from_frequency(
            omega_rad_s,
            FREQUENCY_WAVELENGTH_NM,
            WAVELENGTH_RANGE_NM,
        )
        .ok_or(ReflectionError::Invalid(
            "silica optical frequency is outside the original CHARMS vacuum-wavelength band",
        ))?;
        self.refractive_source.evaluate(wavelength).map_err(|_| {
            ReflectionError::Invalid("original CHARMS silica refraction is unrepresentable")
        })
    }

    /// Absolute n-i*k; sourced intensity loss gives k=alpha*lambda/(4*pi).
    pub fn complex_index(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let phase = self.refractive_sample(omega_rad_s)?.phase_index;
        match self.absorption {
            Some(absorption) => absorption.complex_index(phase, omega_rad_s),
            None => Ok(C64::from_re(phase)),
        }
    }

    /// Absolute SI epsilon=(n-i*k)²/(mu*c0²), retaining the original source mu.
    pub fn permittivity_f_m(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        dispersion::permittivity(
            self.complex_index(omega_rad_s)?,
            self.index_to_sqrt_permittivity,
        )
        .ok_or(ReflectionError::Invalid(
            "silica optical permittivity is unrepresentable",
        ))
    }

    pub(in crate::plane_wave) fn frequency_jet(
        self,
        omega: f64,
    ) -> Result<super::frequency::OpticalFrequencyJet, ReflectionError> {
        super::frequency::silica(self, omega, self.index_to_sqrt_permittivity)
    }

    pub(in crate::plane_wave) fn same_source(self, other: Self) -> bool {
        // The immutable source owns the sole T coordinate and the same mu.
        core::ptr::eq(self.source, other.source)
            && absorption::same_source(self.absorption, other.absorption)
    }
}

#[cfg(test)]
mod tests;
