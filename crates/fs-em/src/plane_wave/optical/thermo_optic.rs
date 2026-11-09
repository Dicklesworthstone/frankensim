//! Source-owned SCHOTT temperature correction in native Maxwell transport.
//!
//! The original absolute spectral reference and the actual material temperature
//! remain distinct. The caller selects zero loss, explicitly holds the original
//! absorption fixed, or supplies independent current-state absorption. No mode
//! reheats the catalog reference atmosphere or replaces dispersion by a fit.

use super::{
    FREQUENCY_WAVELENGTH_NM, ReflectionError, SpectralOpticalResponse, SpectralOpticalSource,
    VACUUM_SPEED_M_S, dispersion, frequency, positive,
};
use fs_material::{
    field_storage::{FieldStorageError, linear_magnetic_permeability},
    state_point::{ResolvedMaterialStatePoint, SpectralDielectricAbsorptionStatePoint},
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

#[derive(Clone, Copy, Debug)]
enum Absorption<'a> {
    Lossless,
    Reference,
    Current {
        source: &'a SpectralDielectricAbsorptionStatePoint,
        maximum_knot_span_m: f64,
    },
}

/// One original glass reference, a source-bounded thermal law and actual state.
///
/// The reference is a dielectric table, analytic Sellmeier source,
/// or explicitly air-corrected catalog law. Material-card references must be
/// resolved at the thermal law's original reference temperature. Both magnetic
/// snapshots must belong to the same immutable material card, but Maxwell
/// transport always uses the actual state's permeability. A catalog's fixed
/// reference atmosphere does not supply the actual glass temperature.
///
/// The lossless constructor refuses nonzero absorption. The separate reference-
/// absorption constructor preserves the original intensity loss at every
/// wavelength while changing only real phase index. This is an explicit frozen-
/// loss approximation, not a supplied temperature-dependent absorption law.
/// The current-absorption constructor instead consumes original alpha claims at
/// the actual state, independently of reference phase claims. That source also
/// owns the actual permeability and every supplied current-state requirement.
/// Expansion, stress, spatial gradients and causal energy reservoirs require
/// additional source laws. Phase transport remains valid at
/// admitted table knots; group/GVD require the original two-sided derivative.
#[derive(Clone, Copy, Debug)]
pub struct ThermoOpticOpticalResponse<'a> {
    reference: SpectralOpticalResponse<'a>,
    law: &'a SchottThermoOpticLaw,
    source: &'a ResolvedMaterialStatePoint,
    absorption: Absorption<'a>,
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
        Self::try_new(reference, law, source, temperature_axis, Absorption::Lossless)
    }

    /// Apply the original real-index thermal law while retaining reference loss.
    ///
    /// The reference source's intensity absorption alpha(lambda) and complete
    /// complex-index derivatives remain unchanged at the same vacuum frequency.
    /// This explicitly selects temperature-independent loss within the original
    /// reference band; it does not infer an absorption law at the new temperature.
    /// Native field heating uses the resulting complete complex permittivity.
    pub fn try_with_reference_absorption(
        reference: SpectralOpticalResponse<'a>,
        law: &'a SchottThermoOpticLaw,
        source: &'a ResolvedMaterialStatePoint,
        temperature_axis: &str,
    ) -> Result<Self> {
        Self::try_new(reference, law, source, temperature_axis, Absorption::Reference)
    }

    /// Combine reference thermal phase with independently sourced CURRENT loss.
    ///
    /// The absorption projection retains same-state magnetic permeability.
    /// Its original support and interpolation admission restrict the PHASE
    /// reference/thermal band, independently of the discarded reference alpha.
    /// No second reference loss or fitted absorption(T) correction is added.
    pub fn try_with_current_absorption(
        reference: SpectralOpticalResponse<'a>,
        law: &'a SchottThermoOpticLaw,
        absorption: &'a SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self> {
        Self::try_new(
            reference,
            law,
            absorption.resolved(),
            temperature_axis,
            Absorption::Current { source: absorption, maximum_knot_span_m },
        )
    }

    fn try_new(
        reference: SpectralOpticalResponse<'a>,
        law: &'a SchottThermoOpticLaw,
        source: &'a ResolvedMaterialStatePoint,
        temperature_axis: &str,
        absorption: Absorption<'a>,
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
                if matches!(absorption, Absorption::Lossless)
                    && value
                        .absorption_samples()
                        .iter()
                        .any(|s| s.absorption_per_m != 0.0)
                {
                    return Err(invalid(
                        "lossless thermo-optic composition cannot discard original absorption",
                    ));
                }
                let (lo, hi) = if matches!(absorption, Absorption::Current { .. }) {
                    value.dispersion_wavelength_bounds_m()
                } else {
                    value.wavelength_bounds_m()
                };
                (constant_m, lo, hi)
            }
            SpectralOpticalSource::TabulatedDielectric(value) => {
                if matches!(absorption, Absorption::Lossless)
                    && value.samples().iter().any(|s| s.absorption_per_m != 0.0)
                {
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
        let mut frequency_range_rad_s = [
            (constant / upper).max(FREQUENCY_WAVELENGTH_NM / hi),
            (constant / lower).min(FREQUENCY_WAVELENGTH_NM / lo),
        ];
        if let Absorption::Current { source, maximum_knot_span_m } = absorption {
            if !maximum_knot_span_m.is_finite()
                || maximum_knot_span_m <= 0.0
                || source.absorption_samples().windows(2).any(|pair| {
                    let gap = pair[1].wavelength_m - pair[0].wavelength_m;
                    !gap.is_finite() || gap <= 0.0 || gap > maximum_knot_span_m
                })
            {
                return Err(invalid(
                    "current absorption exceeds its original interpolation-gap admission",
                ));
            }
            let (lo, hi) = source.wavelength_bounds_m();
            frequency_range_rad_s[0] = frequency_range_rad_s[0].max(constant_m / hi);
            frequency_range_rad_s[1] = frequency_range_rad_s[1].min(constant_m / lo);
        }
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
            absorption,
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
    /// Whether this response explicitly holds the original loss spectrum fixed.
    pub const fn retains_reference_absorption(self) -> bool {
        matches!(self.absorption, Absorption::Reference)
    }
    /// Original independently resolved loss and its complete actual state.
    pub const fn current_absorption(self) -> Option<&'a SpectralDielectricAbsorptionStatePoint> {
        match self.absorption {
            Absorption::Current { source, .. } => Some(source),
            _ => None,
        }
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
        if matches!(self.absorption, Absorption::Current { .. }) {
            return Err(invalid("current absorption must be refreshed with at_absorption_state"));
        }
        Self::try_new(
            self.reference,
            self.law,
            source,
            temperature_axis,
            self.absorption,
        )
    }

    /// Rebind current alpha/mu together while retaining the original phase and
    /// interpolation policies. A bare state cannot silently freeze current loss.
    pub fn at_absorption_state(
        self,
        absorption: &'a SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
    ) -> Result<Self> {
        let Absorption::Current { maximum_knot_span_m, .. } = self.absorption else {
            return Err(invalid("absorption refresh requires the original current-loss policy"));
        };
        Self::try_with_current_absorption(
            self.reference, self.law, absorption, temperature_axis, maximum_knot_span_m,
        )
    }

    fn reference_index(self, omega: f64) -> Result<C64> {
        if matches!(self.absorption, Absorption::Current { .. }) {
            frequency::dielectric_phase_index(self.reference, omega).map(C64::from_re)
        } else {
            self.reference.complex_index(omega)
        }
    }

    /// Phase and thermal derivatives without requiring spectral differentiability.
    pub fn phase_sample(self, omega_rad_s: f64) -> Result<ThermoOpticSample> {
        self.phase_with_reference(omega_rad_s)
            .map(|(sample, _)| sample)
    }

    fn phase_with_reference(self, omega_rad_s: f64) -> Result<(ThermoOpticSample, C64)> {
        let wavelength = self.wavelength(omega_rad_s)?;
        let reference = self.reference_index(omega_rad_s)?;
        self.law
            .prepare(wavelength, reference.re, None)
            .and_then(|value| value.sample(self.temperature_k))
            .map(|sample| (sample, reference))
            .map_err(|_| invalid("original thermo-optic phase response is unavailable"))
    }

    /// Analytic REAL phase/group/curvature at fixed actual material temperature.
    /// Native Maxwell frequency jets additionally retain the selected original
    /// extinction derivatives; real group index is not an energy velocity.
    pub fn refractive_index_jet(self, omega_rad_s: f64) -> Result<ThermoOpticRefractiveIndexJet> {
        self.spectral_with_reference(omega_rad_s)
            .map(|(sample, _, _)| sample)
    }

    fn spectral_with_reference(
        self,
        omega_rad_s: f64,
    ) -> Result<(
        ThermoOpticRefractiveIndexJet,
        C64,
        frequency::OpticalFrequencyJet,
    )> {
        let wavelength = self.wavelength(omega_rad_s)?;
        let [lo, hi] = self.frequency_range_rad_s;
        if omega_rad_s <= lo || omega_rad_s >= hi {
            return Err(invalid(
                "thermo-optic derivatives require a two-sided common source band",
            ));
        }
        let index = self.reference_index(omega_rad_s)?;
        let reference = if matches!(self.absorption, Absorption::Current { .. }) {
            frequency::dielectric_phase(
                self.reference, omega_rad_s, index.re, self.index_to_sqrt_permittivity,
            )?
        } else {
            self.reference.frequency_jet(omega_rad_s)?
        };
        self.law
            .refractive_index_jet(
                wavelength,
                self.temperature_k,
                index.re,
                reference.first_s_m.re * VACUUM_SPEED_M_S,
                reference.second_scaled_s_m.re * VACUUM_SPEED_M_S,
            )
            .map(|sample| (sample, index, reference))
            .map_err(|_| invalid("original thermo-optic frequency derivatives are unavailable"))
    }

    pub fn complex_index(self, omega_rad_s: f64) -> Result<C64> {
        let (sample, reference) = self.phase_with_reference(omega_rad_s)?;
        let extinction = match self.absorption {
            Absorption::Lossless => 0.0,
            Absorption::Reference => reference.im,
            Absorption::Current { source, .. } => {
                let (lo, hi) = source.wavelength_bounds_m();
                let wavelength = dispersion::wavelength_from_frequency(
                    omega_rad_s, core::f64::consts::TAU * VACUUM_SPEED_M_S, [lo, hi],
                ).ok_or_else(|| invalid("frequency is outside current absorption support"))?;
                let alpha = source.absorption_per_m_at(wavelength)
                    .map_err(FieldStorageError::from)?;
                dispersion::attenuating_index_m(sample.phase_index, alpha, wavelength)
                    .ok_or_else(|| invalid("current intensity absorption is unrepresentable"))?.im
            }
        };
        Ok(C64::new(sample.phase_index, extinction))
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
        let (sample, reference_index, reference) = self.spectral_with_reference(omega)?;
        let (extinction, group_loss, curvature_loss) = match self.absorption {
            Absorption::Reference => (
                reference_index.im,
                reference.first_s_m.im * VACUUM_SPEED_M_S,
                reference.second_scaled_s_m.im * VACUUM_SPEED_M_S,
            ),
            Absorption::Lossless => (0.0, 0.0, 0.0),
            Absorption::Current { source, .. } => frequency::dielectric_absorption(source, omega)?,
        };
        // N=n-i*k. Preserve the ORIGINAL imaginary wavelength jet, including
        // alpha slopes: lambda*N_im' = N_im - N_group,im. finish forms all
        // cross terms of epsilon=N²/(mu*c0²) with the heated real index.
        frequency::finish(
            C64::new(sample.phase_index, extinction),
            self.index_to_sqrt_permittivity,
            C64::new(
                sample.scaled_wavelength_derivatives[0],
                extinction - group_loss,
            ),
            C64::new(sample.scaled_wavelength_derivatives[1], curvature_loss),
            C64::new(sample.group_index, group_loss),
        )
    }

    pub(in crate::plane_wave) fn same_source(self, other: Self) -> bool {
        let same_absorption = match (self.absorption, other.absorption) {
            (Absorption::Lossless, Absorption::Lossless)
            | (Absorption::Reference, Absorption::Reference) => true,
            (
                Absorption::Current { source: a, maximum_knot_span_m: ag },
                Absorption::Current { source: b, maximum_knot_span_m: bg },
            ) => core::ptr::eq(a, b) && ag == bg,
            _ => false,
        };
        same_absorption && self.reference.same_source(other.reference)
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
