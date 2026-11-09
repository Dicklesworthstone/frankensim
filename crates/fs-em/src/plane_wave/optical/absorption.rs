//! Original canonical intensity loss composed with an independently owned phase law.

use super::{FREQUENCY_WAVELENGTH_NM, ReflectionError, VACUUM_SPEED_M_S, dispersion};
use fs_material::{
    field_storage::FieldStorageError,
    state_point::SpectralDielectricAbsorptionStatePoint,
};
use fs_math::c64::C64;

/// Prepared source/gap admission; runtime sampling borrows the original table.
#[derive(Clone, Copy, Debug)]
pub(super) struct SourcedAbsorption<'a> {
    source: &'a SpectralDielectricAbsorptionStatePoint,
    maximum_knot_span_m: f64,
    frequency_range_rad_s: [f64; 2],
}

impl<'a> SourcedAbsorption<'a> {
    pub(super) fn try_new(
        source: &'a SpectralDielectricAbsorptionStatePoint,
        maximum_knot_span_m: f64,
        phase_wavelength_range_nm: [f64; 2],
    ) -> Result<Self, ReflectionError> {
        if !maximum_knot_span_m.is_finite()
            || maximum_knot_span_m <= 0.0
            || source.absorption_samples().windows(2).any(|pair| {
                let gap = pair[1].wavelength_m - pair[0].wavelength_m;
                !gap.is_finite() || gap <= 0.0 || gap > maximum_knot_span_m
            })
        {
            return Err(ReflectionError::Invalid(
                "analytic optical absorption exceeds its original interpolation-gap admission",
            ));
        }
        let constant_m = core::f64::consts::TAU * VACUUM_SPEED_M_S;
        let (lo, hi) = source.wavelength_bounds_m();
        let frequency_range_rad_s = [
            (constant_m / hi).max(FREQUENCY_WAVELENGTH_NM / phase_wavelength_range_nm[1]),
            (constant_m / lo).min(FREQUENCY_WAVELENGTH_NM / phase_wavelength_range_nm[0]),
        ];
        if !frequency_range_rad_s.into_iter().all(f64::is_finite)
            || frequency_range_rad_s[0] <= 0.0
            || frequency_range_rad_s[0] >= frequency_range_rad_s[1]
        {
            return Err(ReflectionError::Invalid(
                "analytic phase and original absorption need a common finite wavelength band",
            ));
        }
        Ok(Self {
            source,
            maximum_knot_span_m,
            frequency_range_rad_s,
        })
    }

    pub(super) const fn source(self) -> &'a SpectralDielectricAbsorptionStatePoint {
        self.source
    }

    pub(super) const fn frequency_range_rad_s(self) -> [f64; 2] {
        self.frequency_range_rad_s
    }

    pub(super) fn complex_index(self, phase: f64, omega: f64) -> Result<C64, ReflectionError> {
        if !omega.is_finite()
            || omega < self.frequency_range_rad_s[0]
            || omega > self.frequency_range_rad_s[1]
        {
            return Err(ReflectionError::Invalid(
                "optical frequency is outside the original phase/absorption intersection",
            ));
        }
        let (lo, hi) = self.source.wavelength_bounds_m();
        let wavelength = dispersion::wavelength_from_frequency(
            omega,
            core::f64::consts::TAU * VACUUM_SPEED_M_S,
            [lo, hi],
        )
        .ok_or(ReflectionError::Invalid(
            "optical frequency is outside the original intensity-absorption band",
        ))?;
        let alpha = self
            .source
            .absorption_per_m_at(wavelength)
            .map_err(FieldStorageError::from)?;
        dispersion::attenuating_index_m(phase, alpha, wavelength).ok_or(ReflectionError::Invalid(
            "original intensity absorption gives an unrepresentable complex index",
        ))
    }

    pub(super) fn same_source(self, other: Self) -> bool {
        core::ptr::eq(self.source, other.source)
            && self.maximum_knot_span_m == other.maximum_knot_span_m
    }
}

pub(super) fn same_source(
    left: Option<SourcedAbsorption<'_>>,
    right: Option<SourcedAbsorption<'_>>,
) -> bool {
    match (left, right) {
        (None, None) => true,
        (Some(left), Some(right)) => left.same_source(right),
        _ => false,
    }
}

#[cfg(test)]
mod tests;
