//! Explicit wavelength bindings for an original material's complex index.
//!
//! Each knot binds two complete source keys to a vacuum wavelength. This is
//! an authored interpretation of the selected source data, not a measurement
//! inferred from a material name or from visible-band availability. The whole
//! original material snapshot, including density and other physics, remains
//! intact. Between the declared knots the nominal constitutive model is linear
//! in n and k; neither unsampled spectral structure nor uncertainty is inferred.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

use super::{MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility};

/// Maximum number of independently bound wavelength samples in one spectrum.
pub const MAX_SPECTRAL_COMPLEX_INDEX_KNOTS: usize = 16;

/// One explicitly sourced pair of absolute optical constants.
///
/// The complete keys retain their quantity conventions. Real index requires a
/// dimension-only scalar; extinction additionally admits the existing static
/// `OpticalExtinctionCoefficient` kind. Tensor or hardness contexts cannot be
/// reinterpreted as scalar optical constants.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralComplexIndexKnot {
    /// Vacuum wavelength [nm], positive and increasing across the spectrum.
    pub wavelength_nm: f64,
    /// Exact source key for the absolute real refractive index.
    pub real_index_property: PropertyKey,
    /// Exact source key for the nonnegative extinction coefficient.
    pub extinction_index_property: PropertyKey,
}

/// Nominal complex index at one vacuum wavelength and held material state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralComplexIndexSample {
    /// Vacuum wavelength [nm].
    pub wavelength_nm: f64,
    /// Absolute positive real refractive index.
    pub eta: f64,
    /// Absolute nonnegative extinction coefficient.
    pub k: f64,
}

/// A bounded piecewise-linear optical projection of one complete material
/// snapshot. Its wavelength support is exactly the first through last declared
/// knots; there is no visible-to-infrared extrapolation or default spectrum.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralComplexIndexStatePoint {
    resolved: ResolvedMaterialStatePoint,
    bindings: Vec<SpectralComplexIndexKnot>,
    samples: Vec<SpectralComplexIndexSample>,
    identity: ContentHash,
}

impl SpectralComplexIndexStatePoint {
    /// Adopt 2–16 explicit wavelength/key pairs from an already resolved bundle.
    /// Every original requirement, answer and receipt is retained unchanged.
    /// Optical values may depend on held temperature, pressure, composition or
    /// moisture, but cannot hide a second unowned spectral or dynamic axis.
    ///
    /// # Errors
    /// Refuses missing exact keys, repeated source keys, unordered/nonpositive
    /// wavelengths, wrong scalar conventions, unsupported source coordinates,
    /// nonpositive n or negative k. Repeated equal values are allowed when
    /// separately declared source keys explicitly support those sample roles.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        bindings: &[SpectralComplexIndexKnot],
    ) -> Result<Self, MaterialStatePointError> {
        if bindings.len() < 2 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "complex-index spectrum requires at least two explicit wavelength knots",
            });
        }
        if bindings.len() > MAX_SPECTRAL_COMPLEX_INDEX_KNOTS {
            return Err(MaterialStatePointError::RequirementCount {
                observed: bindings.len(),
                maximum: MAX_SPECTRAL_COMPLEX_INDEX_KNOTS,
            });
        }
        let mut samples = Vec::with_capacity(bindings.len());
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.spectral-complex-index.v1");
        identity.update(resolved.identity().as_bytes());
        identity.update(&(bindings.len() as u64).to_le_bytes());
        for (index, binding) in bindings.iter().enumerate() {
            if !binding.wavelength_nm.is_finite()
                || binding.wavelength_nm <= 0.0
                || index > 0 && binding.wavelength_nm <= bindings[index - 1].wavelength_nm
            {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "complex-index wavelengths must be finite, positive and strictly increasing",
                });
            }
            if binding.real_index_property == binding.extinction_index_property {
                return Err(MaterialStatePointError::DuplicateRequirement {
                    property: binding.real_index_property.name().to_owned(),
                });
            }
            identity.update(&binding.wavelength_nm.to_bits().to_le_bytes());
            for key in [
                &binding.real_index_property,
                &binding.extinction_index_property,
            ] {
                if bindings[..index].iter().any(|previous| {
                    key == &previous.real_index_property
                        || key == &previous.extinction_index_property
                }) {
                    return Err(MaterialStatePointError::DuplicateRequirement {
                        property: key.name().to_owned(),
                    });
                }
                identity.update(&(key.name().len() as u64).to_le_bytes());
                identity.update(key.name().as_bytes());
                identity.update(&key.quantity().canonical_bytes());
            }
            samples.push(SpectralComplexIndexSample {
                wavelength_nm: binding.wavelength_nm,
                eta: source_value(&resolved, &binding.real_index_property, false)?,
                k: source_value(&resolved, &binding.extinction_index_property, true)?,
            });
        }
        Ok(Self {
            resolved,
            bindings: bindings.to_vec(),
            samples,
            identity: identity.finalize(),
        })
    }

    /// Entire original material bundle, including auxiliary properties.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    /// Exact wavelength/source-key interpretation retained by this projection.
    #[must_use]
    pub fn bindings(&self) -> &[SpectralComplexIndexKnot] {
        &self.bindings
    }

    /// Original source-knot values, without resampling or conversion through RGB.
    #[must_use]
    pub fn samples(&self) -> &[SpectralComplexIndexSample] {
        &self.samples
    }

    /// Closed vacuum-wavelength support [nm].
    #[must_use]
    pub fn wavelength_bounds_nm(&self) -> [f64; 2] {
        [
            self.samples[0].wavelength_nm,
            self.samples[self.samples.len() - 1].wavelength_nm,
        ]
    }

    /// Binds the complete original state and ordered wavelength/property roles.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Evaluate the declared piecewise-linear n/k model without material queries
    /// or allocations. Source knots return their exact stored values. This is a
    /// nominal interpolation, not a certificate of unresolved spectral features.
    ///
    /// # Errors
    /// Nonfinite or out-of-support wavelength, or unrepresentable passive index.
    pub fn at_wavelength(
        &self,
        wavelength_nm: f64,
    ) -> Result<SpectralComplexIndexSample, MaterialStatePointError> {
        let [low, high] = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || wavelength_nm < low || wavelength_nm > high {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "complex-index wavelength lies outside its explicit closed support",
            });
        }
        let index = self
            .samples
            .partition_point(|sample| sample.wavelength_nm < wavelength_nm);
        if self.samples[index].wavelength_nm == wavelength_nm {
            return Ok(self.samples[index]);
        }
        let a = self.samples[index - 1];
        let b = self.samples[index];
        let fraction = (wavelength_nm - a.wavelength_nm) / (b.wavelength_nm - a.wavelength_nm);
        // Evaluate from the nearer endpoint so a large positive index contrast
        // cannot cancel a small endpoint value to zero near that endpoint.
        let interpolate = |a: f64, b: f64| {
            if fraction <= 0.5 {
                a + fraction * (b - a)
            } else {
                b + (1.0 - fraction) * (a - b)
            }
        };
        let eta = interpolate(a.eta, b.eta);
        let k = interpolate(a.k, b.k);
        if !ScalarAdmissibility::StrictlyPositive.admits(eta)
            || !ScalarAdmissibility::NonNegative.admits(k)
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "unrepresentable passive complex-index interpolation",
            });
        }
        Ok(SpectralComplexIndexSample {
            wavelength_nm,
            eta,
            k,
        })
    }

    /// Refresh the complete material bundle through its original pinned claims,
    /// then re-admit every optical knot under the same bindings. A failed optical
    /// or auxiliary query returns no partial candidate. This is an endpoint
    /// observation, not temperature evolution or source-support proof between
    /// endpoints.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(self.resolved.at_query_point(card, point)?, &self.bindings)
    }
}

fn source_value(
    resolved: &ResolvedMaterialStatePoint,
    key: &PropertyKey,
    extinction: bool,
) -> Result<f64, MaterialStatePointError> {
    let quantity = key.quantity();
    let typed_extinction = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::OpticalExtinctionCoefficient,
        ValueForm::Static,
    ));
    if !(quantity == QuantitySpec::dimensional(Dims::NONE)
        || extinction && quantity == typed_extinction)
        || key.hardness_test().is_some()
        || key.elastic_component().is_some()
        || key.strain_component().is_some()
        || key.stress_component().is_some()
    {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: key.name().to_owned(),
            reason: "complex-index source role requires its exact static dimensionless scalar convention",
        });
    }
    let property = resolved.property_by_key(key).ok_or_else(|| {
        MaterialStatePointError::InvalidRequirement {
            property: key.name().to_owned(),
            reason: "complex-index source key is absent from the original resolved bundle",
        }
    })?;
    if !property.is_static_at_fixed_external_conditions() {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: key.name().to_owned(),
            reason: "complex-index knots require fixed external source conditions without an unowned dynamic axis",
        });
    }
    let admissibility = if extinction {
        ScalarAdmissibility::NonNegative
    } else {
        ScalarAdmissibility::StrictlyPositive
    };
    let value = property.value_si();
    if !admissibility.admits(value) {
        return Err(MaterialStatePointError::OutsideConsumerDomain {
            property: key.name().to_owned(),
            value,
            admissibility,
        });
    }
    Ok(value)
}

#[cfg(test)]
mod tests;
