//! Ordinary-water phase dispersion in the native monochromatic Maxwell owner.
//!
//! The original IAPWS R9-97 source supplies absolute phase index at actual
//! temperature and mass density, for vacuum wavelengths 200--1100 nm. The
//! caller selects a lossless approximation or independently sourced current
//! absorption; the release supplies neither extinction nor a causal reservoir. Source:
//! <https://iapws.org/technical-guidance/release/Rindex.download>.

use super::{
    FREQUENCY_WAVELENGTH_NM, ReflectionError, VACUUM_SPEED_M_S,
    absorption::{self, SourcedAbsorption}, dispersion, positive,
};
use fs_material::{
    field_storage::linear_magnetic_permeability,
    state_point::{
        DENSITY_PROPERTY, ResolvedMaterialStatePoint, SpectralDielectricAbsorptionStatePoint,
    },
    water_refraction::{IapwsWaterRefraction, WAVELENGTH_RANGE_NM, WaterRefractiveSample},
};
use fs_math::{c64::C64, det};
use fs_qty::{
    Density, QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};

/// Original ordinary-water refraction with an explicit optical loss policy.
///
/// Temperature, bulk mass density and absolute scalar permeability all come
/// from the SAME complete resolved bundle. No pressure-to-density closure,
/// ambient index ratio, phase classification or nonmagnetic default is added.
/// Selecting [`Self::try_lossless`] asserts ordinary-water applicability and
/// chooses zero absorption; it does not infer purity from the material name or
/// claim that real water has zero absorption over this band. Mixtures, ice,
/// scattering and mixtures require their own owners. The separate absorption
/// constructor uses the canonical original alpha spectrum at the SAME T/rho/mu
/// state; its source band further restricts monochromatic Maxwell transport.
///
/// This finite-band, homogeneous, isotropic response retains IAPWS phase
/// dispersion for interfaces, films and monochromatic angular spectra. Its
/// group index is available separately; phase index alone forms Maxwell's
/// epsilon. It does not supply a causal transient storage or heat law.
#[derive(Clone, Copy, Debug)]
pub struct WaterOpticalResponse<'a> {
    source: &'a ResolvedMaterialStatePoint,
    absorption: Option<SourcedAbsorption<'a>>,
    temperature_k: f64,
    density_kg_m3: f64,
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> WaterOpticalResponse<'a> {
    /// Adopt a source state and explicitly choose zero optical absorption.
    ///
    /// `temperature_axis` names the unique static absolute-temperature query
    /// coordinate. The legacy untyped/dimension-only `T` axis is also kelvin;
    /// another name requires explicit absolute-temperature semantics. The
    /// canonical `density` property is bulk kg/m3, with its original static
    /// source admission. Source-dependent frequency/field laws refuse.
    /// The native IAPWS temperature/density domain is checked immediately.
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
            if is_absolute {
                if axis != temperature_axis || temperature.replace(*value).is_some() {
                    return Err(ReflectionError::Invalid(
                        "water optics requires one owned absolute-temperature coordinate",
                    ));
                }
            }
        }
        let temperature_k = temperature.ok_or(ReflectionError::Invalid(
            "water optics needs the original absolute temperature in kelvin",
        ))?;
        let density = source
            .property(DENSITY_PROPERTY)
            .ok_or(ReflectionError::Invalid(
                "water optics needs original bulk mass density in its complete source state",
            ))?;
        if density.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
            || !density.has_only_fixed_external_conditions(Some(temperature_axis), false)
        {
            return Err(ReflectionError::Invalid(
                "water optics needs original stationary bulk mass density in kg/m3",
            ));
        }
        let density_kg_m3 = density.value_si();
        let permeability_h_m = linear_magnetic_permeability(source)?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        let response = Self {
            source,
            absorption: None,
            temperature_k,
            density_kg_m3,
            permeability_h_m,
            index_to_sqrt_permittivity,
        };
        // Validate the state and representable SI epsilon before exposing the
        // response. Runtime evaluations retain the same exact source bounds.
        for wavelength_nm in WAVELENGTH_RANGE_NM {
            response.permittivity_f_m(FREQUENCY_WAVELENGTH_NM / wavelength_nm)?;
        }
        Ok(response)
    }

    /// Combine original IAPWS phase with source-supported intensity absorption.
    ///
    /// The complete absorption snapshot must also contain the original bulk
    /// density and scalar permeability. Its original temperature owns BOTH
    /// phase and loss; no separately supplied optical state can disagree.
    /// The positive finite gap admission is in vacuum metres. No scattering,
    /// phase correction or causal polarization reservoir is inferred from alpha.
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

    pub(in crate::plane_wave) fn frequency_jet(
        self,
        omega: f64,
    ) -> Result<super::frequency::OpticalFrequencyJet, ReflectionError> {
        super::frequency::water(self, omega, self.index_to_sqrt_permittivity)
    }

    /// Complete original T/rho/mu bundle, including its source-use identity.
    pub const fn source(self) -> &'a ResolvedMaterialStatePoint {
        self.source
    }

    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.source
    }

    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }

    pub const fn density_kg_m3(self) -> f64 {
        self.density_kg_m3
    }

    pub const fn permeability_h_m(self) -> f64 {
        self.permeability_h_m
    }

    /// Exact IAPWS phase/group indices and analytic partials at vacuum frequency.
    /// No source evaluation allocates. Derivatives retain the material owner's
    /// fixed-density/fixed-temperature meanings, not a guessed pressure path.
    pub fn refractive_sample(
        self,
        omega_rad_s: f64,
    ) -> Result<WaterRefractiveSample, ReflectionError> {
        let wavelength_nm = dispersion::wavelength_from_frequency(
            omega_rad_s,
            FREQUENCY_WAVELENGTH_NM,
            WAVELENGTH_RANGE_NM,
        )
        .ok_or(ReflectionError::Invalid(
            "water optical frequency is outside the original IAPWS vacuum-wavelength band",
        ))?;
        IapwsWaterRefraction::prepare(wavelength_nm)
            .and_then(|prepared| prepared.sample(self.temperature_k, self.density_kg_m3))
            .map_err(|_| {
                ReflectionError::Invalid(
                    "water optical state is outside the original IAPWS domain or is unrepresentable",
                )
            })
    }

    /// Absolute n-i*k, with k=alpha*lambda/(4*pi) for the sourced-loss policy.
    /// The explicit lossless approximation uses k=0 over the original phase band.
    pub fn complex_index(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let phase = self.refractive_sample(omega_rad_s)?.phase_index;
        match self.absorption {
            Some(absorption) => absorption.complex_index(phase, omega_rad_s),
            None => Ok(C64::from_re(phase)),
        }
    }

    /// Absolute SI epsilon = (n-i*k)^2/(mu*c0^2), in F/m.
    /// The permeability is the retained source value, not an implicit mu0.
    pub fn permittivity_f_m(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        dispersion::permittivity(
            self.complex_index(omega_rad_s)?,
            self.index_to_sqrt_permittivity,
        )
        .ok_or(ReflectionError::Invalid(
            "water optical permittivity is unrepresentable",
        ))
    }

    pub(in crate::plane_wave) fn same_source(self, other: Self) -> bool {
        // T/rho/mu are immutable projections of this exact complete snapshot.
        // Re-selecting a source at another temperature/density cannot match it.
        core::ptr::eq(self.source, other.source)
            && absorption::same_source(self.absorption, other.absorption)
    }
}

#[cfg(test)]
mod tests;
