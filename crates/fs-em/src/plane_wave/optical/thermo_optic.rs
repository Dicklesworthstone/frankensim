//! Source-owned SCHOTT temperature correction in native Maxwell transport.
//!
//! The original absolute spectral reference and the actual material temperature
//! remain distinct. This explicitly lossless composition never reheats the
//! catalog reference atmosphere or replaces the original dispersion by a fit.

use super::{
    FREQUENCY_WAVELENGTH_NM, ReflectionError, SpectralOpticalResponse, SpectralOpticalSource,
    VACUUM_SPEED_M_S, dispersion, frequency, positive,
};
use fs_material::{
    field_storage::linear_magnetic_permeability,
    state_point::ResolvedMaterialStatePoint,
    thermo_optic::{SchottThermoOpticLaw, ThermoOpticRefractiveIndexJet, ThermoOpticSample},
};
use fs_math::{c64::C64, det};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};

type Result<T> = core::result::Result<T, ReflectionError>;

fn invalid(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}

/// One original glass reference, a source-bounded thermal law and actual state.
///
/// The reference is a transparent dielectric table, analytic Sellmeier source,
/// or explicitly air-corrected catalog law. Material-card references must be
/// resolved at the thermal law's original reference temperature. Both magnetic
/// snapshots must belong to the same immutable material card, but Maxwell
/// transport always uses the actual state's permeability. A catalog's fixed
/// reference atmosphere does not supply the actual glass temperature.
///
/// Preparation refuses nonzero absorption instead of silently deleting it.
/// Thermal absorption, expansion, stress, spatial gradients and causal energy
/// reservoirs require additional source laws. Phase transport remains valid at
/// admitted table knots; group/GVD require the original two-sided derivative.
#[derive(Clone, Copy, Debug)]
pub struct ThermoOpticOpticalResponse<'a> {
    reference: SpectralOpticalResponse<'a>,
    law: &'a SchottThermoOpticLaw,
    source: &'a ResolvedMaterialStatePoint,
    temperature_k: f64,
    frequency_range_rad_s: [f64; 2],
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> ThermoOpticOpticalResponse<'a> {
    /// Select a lossless approximation without fitting or accumulating increments.
    ///
    /// Only one static absolute-temperature coordinate is supported. A named
    /// axis requires absolute-temperature semantics; the existing untyped or
    /// dimensional `T` convention is kelvin. Pressure/composition/stress axes
    /// cannot be silently carried through a temperature-only glass law.
    pub fn try_lossless(
        reference: SpectralOpticalResponse<'a>,
        law: &'a SchottThermoOpticLaw,
        source: &'a ResolvedMaterialStatePoint,
        temperature_axis: &str,
    ) -> Result<Self> {
        let temperature_k = temperature(source, temperature_axis)?;
        let [min_t, max_t] = law.temperature_range_k();
        if !temperature_k.is_finite() || !(min_t..=max_t).contains(&temperature_k) {
            return Err(invalid(
                "glass temperature is outside the original thermo-optic source",
            ));
        }
        if reference.material().card_identity() != source.card_identity()
            || reference.material().material() != source.material()
        {
            return Err(invalid(
                "thermal glass reference and actual state require the same original material card",
            ));
        }
        let constant_m = core::f64::consts::TAU * VACUUM_SPEED_M_S;
        let (constant, lower, upper) = match reference.source() {
            SpectralOpticalSource::DielectricSellmeier(value) => {
                if value
                    .absorption_samples()
                    .iter()
                    .any(|s| s.absorption_per_m != 0.0)
                {
                    return Err(invalid(
                        "lossless thermo-optic composition cannot discard original absorption",
                    ));
                }
                let (lo, hi) = value.wavelength_bounds_m();
                (constant_m, lo, hi)
            }
            SpectralOpticalSource::TabulatedDielectric(value) => {
                if value.samples().iter().any(|s| s.absorption_per_m != 0.0) {
                    return Err(invalid(
                        "lossless thermo-optic composition cannot discard original absorption",
                    ));
                }
                let (lo, hi) = value.wavelength_bounds_m();
                (constant_m, lo, hi)
            }
            SpectralOpticalSource::LosslessAirRelativeSellmeier { law, .. } => {
                let (lo, hi) = law.wavelength_bounds_nm();
                (FREQUENCY_WAVELENGTH_NM, lo, hi)
            }
            _ => {
                return Err(invalid(
                    "thermo-optic glass requires an original dielectric reference",
                ));
            }
        };
        if !matches!(
            reference.source(),
            SpectralOpticalSource::LosslessAirRelativeSellmeier { .. }
        ) && (temperature(reference.material(), temperature_axis)?
            != law.reference_temperature_k()
            || reference.material().axis_quantities() != source.axis_quantities())
        {
            return Err(invalid(
                "glass reference must retain the original thermal reference temperature and coordinates",
            ));
        }
        let [lo, hi] = law.wavelength_range_nm();
        // Intersect in frequency using EACH source's original units. A metres
        // to nanometres round trip must not expand or accidentally lose an edge.
        let frequency_range_rad_s = [
            (constant / upper).max(FREQUENCY_WAVELENGTH_NM / hi),
            (constant / lower).min(FREQUENCY_WAVELENGTH_NM / lo),
        ];
        if frequency_range_rad_s
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || frequency_range_rad_s[0] > frequency_range_rad_s[1]
        {
            return Err(invalid(
                "glass reference and thermo-optic law have no representable common band",
            ));
        }
        let permeability_h_m = linear_magnetic_permeability(source)?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        let response = Self {
            reference,
            law,
            source,
            temperature_k,
            frequency_range_rad_s,
            permeability_h_m,
            index_to_sqrt_permittivity,
        };
        for omega in frequency_range_rad_s {
            response.permittivity_f_m(omega)?;
        }
        Ok(response)
    }

    pub const fn reference(self) -> SpectralOpticalResponse<'a> {
        self.reference
    }
    pub const fn law(self) -> &'a SchottThermoOpticLaw {
        self.law
    }
    pub const fn source(self) -> &'a ResolvedMaterialStatePoint {
        self.source
    }
    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.source
    }
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }
    pub const fn frequency_range_rad_s(self) -> [f64; 2] {
        self.frequency_range_rad_s
    }
    pub const fn permeability_h_m(self) -> f64 {
        self.permeability_h_m
    }

    /// Rebind actual temperature and mu, always from the immutable reference.
    pub fn at_state(
        self,
        source: &'a ResolvedMaterialStatePoint,
        temperature_axis: &str,
    ) -> Result<Self> {
        Self::try_lossless(self.reference, self.law, source, temperature_axis)
    }

    /// Phase and thermal derivatives without requiring spectral differentiability.
    pub fn phase_sample(self, omega_rad_s: f64) -> Result<ThermoOpticSample> {
        let wavelength = self.wavelength(omega_rad_s)?;
        let n0 = self.reference.complex_index(omega_rad_s)?.re;
        self.law
            .prepare(wavelength, n0, None)
            .and_then(|value| value.sample(self.temperature_k))
            .map_err(|_| invalid("original thermo-optic phase response is unavailable"))
    }

    /// Analytic phase/group/curvature at fixed actual material temperature.
    pub fn refractive_index_jet(self, omega_rad_s: f64) -> Result<ThermoOpticRefractiveIndexJet> {
        let wavelength = self.wavelength(omega_rad_s)?;
        let [lo, hi] = self.frequency_range_rad_s;
        if omega_rad_s <= lo || omega_rad_s >= hi {
            return Err(invalid(
                "thermo-optic derivatives require a two-sided common source band",
            ));
        }
        let reference = self.reference.frequency_jet(omega_rad_s)?;
        let n0 = self.reference.complex_index(omega_rad_s)?.re;
        self.law
            .refractive_index_jet(
                wavelength,
                self.temperature_k,
                n0,
                reference.first_s_m.re * VACUUM_SPEED_M_S,
                reference.second_scaled_s_m.re * VACUUM_SPEED_M_S,
            )
            .map_err(|_| invalid("original thermo-optic frequency derivatives are unavailable"))
    }

    pub fn complex_index(self, omega_rad_s: f64) -> Result<C64> {
        Ok(C64::from_re(self.phase_sample(omega_rad_s)?.phase_index))
    }

    pub fn permittivity_f_m(self, omega_rad_s: f64) -> Result<C64> {
        dispersion::permittivity(
            self.complex_index(omega_rad_s)?,
            self.index_to_sqrt_permittivity,
        )
        .ok_or_else(|| invalid("thermo-optic permittivity is unrepresentable"))
    }

    fn wavelength(self, omega_rad_s: f64) -> Result<f64> {
        let [lo, hi] = self.frequency_range_rad_s;
        if !omega_rad_s.is_finite() || !(lo..=hi).contains(&omega_rad_s) {
            return Err(invalid(
                "frequency is outside the original common thermal glass band",
            ));
        }
        dispersion::wavelength_from_frequency(
            omega_rad_s,
            FREQUENCY_WAVELENGTH_NM,
            self.law.wavelength_range_nm(),
        )
        .ok_or_else(|| invalid("thermo-optic vacuum wavelength is unrepresentable"))
    }

    pub(in crate::plane_wave) fn frequency_jet(
        self,
        omega: f64,
    ) -> Result<frequency::OpticalFrequencyJet> {
        let sample = self.refractive_index_jet(omega)?;
        frequency::finish(
            C64::from_re(sample.phase_index),
            self.index_to_sqrt_permittivity,
            C64::from_re(sample.scaled_wavelength_derivatives[0]),
            C64::from_re(sample.scaled_wavelength_derivatives[1]),
            C64::from_re(sample.group_index),
        )
    }

    pub(in crate::plane_wave) fn same_source(self, other: Self) -> bool {
        self.reference.same_source(other.reference)
            && core::ptr::eq(self.law, other.law)
            && core::ptr::eq(self.source, other.source)
    }
}

fn temperature(source: &ResolvedMaterialStatePoint, temperature_axis: &str) -> Result<f64> {
    let [(axis, value)] = source.query_point() else {
        return Err(invalid(
            "thermo-optic glass needs one original absolute-temperature coordinate",
        ));
    };
    let quantity = source.axis_quantities().get(axis).copied();
    let absolute = quantity
        .and_then(QuantitySpec::semantic_type)
        .is_some_and(|ty| {
            ty.kind() == QuantityKind::AbsoluteTemperature && ty.form() == ValueForm::Static
        })
        || axis == "T"
            && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
    if axis != temperature_axis || !absolute {
        return Err(invalid(
            "thermo-optic glass requires the named static absolute temperature in kelvin",
        ));
    }
    Ok(*value)
}

#[cfg(test)]
mod tests;
