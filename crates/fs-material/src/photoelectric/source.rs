//! Photoemission from an explicitly prepared, source-bound film.
//!
//! A bulk material name does not determine a surface work function. The
//! retained material card must describe the preparation, surface condition and
//! original calibration of this emitter. Its thickness is part of that
//! calibration: the supplied internal quantum yield includes transport and
//! escape losses through that entire film. It is electrons emitted into vacuum
//! per photon absorbed in that film, not external collected-electron QE.
//!
//! Each wavelength explicitly selects two original scalar claims: internal
//! yield and the mean fraction of the available excess energy carried away by
//! emitted electrons. Neither quantity is inferred from optical extinction or
//! from the maximum electron energy. Source knots may vary with held material
//! conditions; the declared wavelength interpolation is nominal, not evidence
//! for unsampled spectral structure. No extrapolation or default tail is used.

use super::{PhotoelectricError, PhotoelectricLaw};
use crate::state_point::{MaterialStatePointError, ResolvedMaterialStatePoint};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

/// Bounded source table; each wavelength consumes two original property roles.
pub const MAX_PHOTOELECTRIC_KNOTS: usize = 16;

/// An authored vacuum-wavelength interpretation of exact source property keys.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotoelectricKnot {
    /// Vacuum wavelength in nanometres, positive and strictly increasing.
    pub wavelength_nm: f64,
    /// Emitted electrons per photon absorbed in the calibrated film, in [0,1].
    pub absorbed_quantum_yield_property: PropertyKey,
    /// Mean emitted kinetic energy divided by available excess energy, in [0,1].
    pub mean_excess_fraction_property: PropertyKey,
}

/// Original source, finite-band or numerical constitutive refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum PhotoelectricSourceError {
    Material(MaterialStatePointError),
    Physics(PhotoelectricError),
    Invalid(&'static str),
}

impl core::fmt::Display for PhotoelectricSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "photoelectric material: {self:?}")
    }
}
impl std::error::Error for PhotoelectricSourceError {}
impl From<MaterialStatePointError> for PhotoelectricSourceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<PhotoelectricError> for PhotoelectricSourceError {
    fn from(error: PhotoelectricError) -> Self {
        Self::Physics(error)
    }
}

/// Original finite-film photoemission law at one complete physical state.
/// The identity binds source claims and their authored wavelength/role mapping;
/// it does not confer experimental accuracy on a synthetic or fitted source.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotoelectricMaterial {
    resolved: ResolvedMaterialStatePoint,
    work_function_property: PropertyKey,
    thickness_property: PropertyKey,
    knots: Vec<PhotoelectricKnot>,
    values: Vec<[f64; 2]>,
    maximum_gap_nm: f64,
    work_function_j: f64,
    thickness_m: f64,
    identity: ContentHash,
}

impl PhotoelectricMaterial {
    /// Bind exact source roles without selecting another material claim.
    ///
    /// Work function is joules per electron (not eV or J/mol); thickness is
    /// metres. Fractions require static dimensionless scalars. Source keys may
    /// not be reused for distinct physical roles. Optical, caloric and other
    /// properties in the original bundle remain intact for domain consumers.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        work_function_property: PropertyKey,
        thickness_property: PropertyKey,
        knots: &[PhotoelectricKnot],
        maximum_gap_nm: f64,
    ) -> Result<Self, PhotoelectricSourceError> {
        let invalid = PhotoelectricSourceError::Invalid;
        if !(2..=MAX_PHOTOELECTRIC_KNOTS).contains(&knots.len()) {
            return Err(invalid(
                "photoemission requires 2..=16 explicit wavelength knots",
            ));
        }
        if !maximum_gap_nm.is_finite() || maximum_gap_nm <= 0.0 {
            return Err(invalid(
                "photoemission interpolation gap must be finite and positive",
            ));
        }
        let work_function_j = source_value(
            &resolved,
            &work_function_property,
            Dims([2, 1, -2, 0, 0, 0]),
        )?;
        let thickness_m = source_value(&resolved, &thickness_property, fs_qty::Length::DIMS)?;
        if work_function_j <= 0.0 || thickness_m <= 0.0 {
            return Err(invalid(
                "source work function and calibrated film thickness must be positive",
            ));
        }
        let mut identity = DomainHasher::new("org.frankensim.material.photoelectric-film.v1");
        identity.update(resolved.identity().as_bytes());
        identity.update(&maximum_gap_nm.to_bits().to_le_bytes());
        identity.update(&(knots.len() as u64).to_le_bytes());
        let mut keys = vec![&work_function_property, &thickness_property];
        let mut values = Vec::with_capacity(knots.len());
        for (index, knot) in knots.iter().enumerate() {
            if !knot.wavelength_nm.is_finite()
                || knot.wavelength_nm <= 0.0
                || index > 0 && knot.wavelength_nm <= knots[index - 1].wavelength_nm
            {
                return Err(invalid(
                    "photoemission wavelengths must be positive and strictly increasing",
                ));
            }
            if index > 0 && knot.wavelength_nm - knots[index - 1].wavelength_nm > maximum_gap_nm {
                return Err(invalid(
                    "photoemission source gap exceeds the declared interpolation limit",
                ));
            }
            let yield_value =
                source_value(&resolved, &knot.absorbed_quantum_yield_property, Dims::NONE)?;
            let mean_fraction =
                source_value(&resolved, &knot.mean_excess_fraction_property, Dims::NONE)?;
            PhotoelectricLaw::try_new(work_function_j, yield_value, mean_fraction)?;
            values.push([yield_value, mean_fraction]);
            keys.push(&knot.absorbed_quantum_yield_property);
            keys.push(&knot.mean_excess_fraction_property);
            identity.update(&knot.wavelength_nm.to_bits().to_le_bytes());
        }
        for (index, key) in keys.iter().enumerate() {
            if keys[..index].contains(key) {
                return Err(invalid(
                    "photoemission source roles must use distinct exact keys",
                ));
            }
            identity.update(&(key.name().len() as u64).to_le_bytes());
            identity.update(key.name().as_bytes());
            identity.update(&key.quantity().canonical_bytes());
        }
        Ok(Self {
            resolved,
            work_function_property,
            thickness_property,
            knots: knots.to_vec(),
            values,
            maximum_gap_nm,
            work_function_j,
            thickness_m,
            identity: identity.finalize(),
        })
    }

    /// Entire original immutable source snapshot, including other physics.
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Energy required to remove one electron into vacuum, in joules.
    pub const fn work_function_j(&self) -> f64 {
        self.work_function_j
    }

    /// Source-calibrated emitting film thickness. Changing it requires a new source.
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    pub fn knots(&self) -> &[PhotoelectricKnot] {
        &self.knots
    }

    pub fn wavelength_bounds_nm(&self) -> [f64; 2] {
        [
            self.knots[0].wavelength_nm,
            self.knots[self.knots.len() - 1].wavelength_nm,
        ]
    }

    /// Allocation-free original-band interpolation followed by the shared law.
    /// Thresholding remains in `PhotoelectricLaw::evaluate` at this wavelength;
    /// broadband consumers must evaluate each physical spectral line separately.
    pub fn at_wavelength(
        &self,
        wavelength_nm: f64,
    ) -> Result<PhotoelectricLaw, PhotoelectricSourceError> {
        let [low, high] = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || wavelength_nm < low || wavelength_nm > high {
            return Err(PhotoelectricSourceError::Invalid(
                "photoemission wavelength is outside original closed support",
            ));
        }
        let index = self
            .knots
            .partition_point(|knot| knot.wavelength_nm < wavelength_nm);
        let values = if self.knots[index].wavelength_nm == wavelength_nm {
            self.values[index]
        } else {
            let fraction = (wavelength_nm - self.knots[index - 1].wavelength_nm)
                / (self.knots[index].wavelength_nm - self.knots[index - 1].wavelength_nm);
            core::array::from_fn(|role| {
                let a = self.values[index - 1][role];
                let b = self.values[index][role];
                if fraction <= 0.5 {
                    a + fraction * (b - a)
                } else {
                    b + (1.0 - fraction) * (a - b)
                }
            })
        };
        Ok(PhotoelectricLaw::try_new(
            self.work_function_j,
            values[0],
            values[1],
        )?)
    }

    /// Refresh the original claims atomically at a new supported physical point.
    /// An evolving consumer must separately account for changed removal energy,
    /// thickness or stored charge; this operation advances no physical history.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, PhotoelectricSourceError> {
        Self::try_from_resolved(
            self.resolved.at_query_point(card, point)?,
            self.work_function_property.clone(),
            self.thickness_property.clone(),
            &self.knots,
            self.maximum_gap_nm,
        )
    }
}

fn source_value(
    resolved: &ResolvedMaterialStatePoint,
    key: &PropertyKey,
    dims: Dims,
) -> Result<f64, PhotoelectricSourceError> {
    let invalid = PhotoelectricSourceError::Invalid;
    if key.quantity() != QuantitySpec::dimensional(dims)
        || key.hardness_test().is_some()
        || key.elastic_component().is_some()
        || key.strain_component().is_some()
        || key.stress_component().is_some()
    {
        return Err(invalid(
            "photoemission source role has incompatible quantity or tensor context",
        ));
    }
    let value = resolved.property_by_key(key).ok_or(invalid(
        "photoemission exact source key is absent from the original bundle",
    ))?;
    if !value.is_static_at_fixed_external_conditions() || !value.value_si().is_finite() {
        return Err(invalid(
            "photoemission source must be finite at fixed, explicitly owned external conditions",
        ));
    }
    Ok(value.value_si())
}

#[cfg(test)]
mod tests;
