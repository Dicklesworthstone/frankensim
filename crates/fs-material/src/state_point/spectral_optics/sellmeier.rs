//! Analytic Sellmeier source with independent fitted and absorption domains.
//!
//! The six oscillator coefficients retain their original SI convention. Two
//! additional source properties delimit the fit's vacuum-wavelength support;
//! the presence of absorption samples never establishes that support. Native
//! optical consumers admit the complete transparent fit band before using its
//! intersection with the original intensity-absorption spectrum.

use super::{
    INDEX_QUANTITY, MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    SPECTRAL_OPTICAL_SAMPLE_COUNT, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility,
    ScalarPropertyRequirement, WAVELENGTH_QUANTITY, resolve_material_state_point, source_value,
    source_wavelengths,
};
use crate::state_point::{
    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

/// Explicit lower and upper vacuum-wavelength limits of the original analytic
/// fit, in SI metres. These are separate from the absorption sample coordinates.
pub const SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES: [&str; 2] = [
    "optical_sellmeier_wavelength_min_m",
    "optical_sellmeier_wavelength_max_m",
];

/// One original internal intensity-absorption sample, without a fitted index.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDielectricAbsorptionSample {
    /// Vacuum wavelength in metres.
    pub wavelength_m: f64,
    /// Internal intensity absorption in power nepers per metre.
    pub absorption_per_m: f64,
}

/// Original analytic dispersion and independent spectral absorption at one
/// complete material state. No sampled phase-index approximation is introduced.
///
/// The source must express absolute indices relative to vacuum. This resolver
/// admits quantities, nonnegative coefficients and explicit finite domains;
/// consumers remain responsible for pole-free positive-index admission of the
/// complete fitted band. No RGB fallback, default transparency, relative-air
/// conversion, or material-name inference supplies missing source data.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralSellmeierStatePoint {
    resolved: ResolvedMaterialStatePoint,
    strengths: [f64; 3],
    resonances_m2: [f64; 3],
    dispersion_bounds_m: (f64, f64),
    absorption: [SpectralDielectricAbsorptionSample; SPECTRAL_OPTICAL_SAMPLE_COUNT],
}

impl SpectralSellmeierStatePoint {
    /// Adopt the source roles from one original complete snapshot, preserving
    /// every additional physical property and receipt. Requirements are checked
    /// again even if the original consumer admitted more permissive scalars.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        let wavelengths = source_wavelengths(&resolved)?;
        let mut strengths = [0.0; 3];
        let mut resonances_m2 = [0.0; 3];
        for index in 0..3 {
            strengths[index] = source_value(
                &resolved,
                VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[index],
                INDEX_QUANTITY,
                ScalarAdmissibility::NonNegative,
            )?;
            resonances_m2[index] = source_value(
                &resolved,
                VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[index],
                QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
                ScalarAdmissibility::NonNegative,
            )?;
        }
        let mut bounds = [0.0; 2];
        for index in 0..2 {
            bounds[index] = source_value(
                &resolved,
                SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[index],
                WAVELENGTH_QUANTITY,
                ScalarAdmissibility::StrictlyPositive,
            )?;
        }
        if bounds[0] > bounds[1] || bounds[0].max(wavelengths[0]) > bounds[1].min(wavelengths[8]) {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "sellmeier fit and absorption require ordered overlapping source bands",
            });
        }
        let mut absorption = [SpectralDielectricAbsorptionSample {
            wavelength_m: 0.0,
            absorption_per_m: 0.0,
        }; SPECTRAL_OPTICAL_SAMPLE_COUNT];
        for (index, sample) in absorption.iter_mut().enumerate() {
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
        Ok(Self {
            resolved,
            strengths,
            resonances_m2,
            dispersion_bounds_m: (bounds[0], bounds[1]),
            absorption,
        })
    }

    /// Original complete material state, including all selected source receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Analytic family and all original source values, domains and quantities.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.material.spectral-sellmeier.v1");
        hash.update(self.resolved.identity().as_bytes());
        hash.finalize()
    }

    /// Original dimensionless strengths and squared resonance wavelengths [m²].
    #[must_use]
    pub const fn sellmeier_coefficients_si(&self) -> ([f64; 3], [f64; 3]) {
        (self.strengths, self.resonances_m2)
    }

    /// Original inclusive fitted-index support in vacuum metres. Consumers must
    /// admit this entire band, even where absorption has narrower support.
    #[must_use]
    pub const fn dispersion_wavelength_bounds_m(&self) -> (f64, f64) {
        self.dispersion_bounds_m
    }

    /// Common inclusive support of the analytic fit and absorption spectrum.
    #[must_use]
    pub fn wavelength_bounds_m(&self) -> (f64, f64) {
        (
            self.dispersion_bounds_m
                .0
                .max(self.absorption[0].wavelength_m),
            self.dispersion_bounds_m
                .1
                .min(self.absorption[8].wavelength_m),
        )
    }

    /// Original independent intensity-loss samples, with their SI coordinates.
    #[must_use]
    pub const fn absorption_samples(
        &self,
    ) -> &[SpectralDielectricAbsorptionSample; SPECTRAL_OPTICAL_SAMPLE_COUNT] {
        &self.absorption
    }

    /// Evaluate the original piecewise-linear absorption spectrum without
    /// allocations. Its domain is independent of the fitted-index domain;
    /// consumers requiring both must also admit `wavelength_bounds_m()`.
    pub fn absorption_per_m_at(&self, wavelength_m: f64) -> Result<f64, MaterialStatePointError> {
        if !wavelength_m.is_finite()
            || wavelength_m < self.absorption[0].wavelength_m
            || wavelength_m > self.absorption[8].wavelength_m
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "wavelength outside the original Sellmeier absorption spectrum",
            });
        }
        let upper = self
            .absorption
            .partition_point(|sample| sample.wavelength_m < wavelength_m);
        if self.absorption[upper].wavelength_m == wavelength_m {
            return Ok(self.absorption[upper].absorption_per_m);
        }
        let [left, right] = [self.absorption[upper - 1], self.absorption[upper]];
        let t = (wavelength_m - left.wavelength_m) / (right.wavelength_m - left.wavelength_m);
        let [a, b] = [left.absorption_per_m, right.absorption_per_m];
        let value = if t <= 0.5 {
            a + t * (b - a)
        } else {
            b + (1.0 - t) * (a - b)
        };
        if !ScalarAdmissibility::NonNegative.admits(value) {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "unrepresentable Sellmeier-source absorption interpolation",
            });
        }
        Ok(value)
    }

    /// Re-evaluate the exact selected coefficient, bound and absorption claims
    /// at one new complete physical condition. No source switching, resampling
    /// or support along the intervening state path is inferred.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(self.resolved.at_query_point(card, point)?)
    }
}

/// Resolve 26 explicit source properties: six coefficients, two fit limits,
/// nine absorption coordinates and nine typed intensity-loss values. Legacy
/// visible family selection remains separate; every zero must be declared.
pub fn resolve_spectral_sellmeier_state_point(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<SpectralSellmeierStatePoint, MaterialStatePointError> {
    let mut requirements = Vec::with_capacity(26);
    for (names, quantity, admissibility) in [
        (
            &VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[..],
            INDEX_QUANTITY,
            ScalarAdmissibility::NonNegative,
        ),
        (
            &VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[..],
            QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
            ScalarAdmissibility::NonNegative,
        ),
        (
            &SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[..],
            WAVELENGTH_QUANTITY,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            &SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[..],
            WAVELENGTH_QUANTITY,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            &SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[..],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            ScalarAdmissibility::NonNegative,
        ),
    ] {
        for name in names {
            requirements.push(ScalarPropertyRequirement::try_with_quantity(
                *name,
                quantity,
                admissibility,
            )?);
        }
    }
    SpectralSellmeierStatePoint::try_from_material(resolve_material_state_point(
        card,
        point,
        &requirements,
        selection,
    )?)
}

#[cfg(test)]
mod tests;
