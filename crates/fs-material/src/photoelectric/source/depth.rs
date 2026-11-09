//! Depth-conditioned emission from one explicitly prepared finite film.
//!
//! Each bin supplies the probability of escape through the selected vacuum face
//! **conditional on absorption in that bin**, and the mean excess-energy fraction
//! of those escaping electrons. These are not whole-film effective yields. The
//! optical consumer integrates its actual absorption over every bin before
//! applying the corresponding law, retaining changes in standing-wave profiles.
//!
//! Bins cover the full film in distance from the declared escape face. Fractions
//! are authored coordinates on the source-calibrated thickness, not an inferred
//! escape length. Piecewise constancy in depth and linear spectral interpolation
//! are explicit constitutive approximations; no unresolved depth or spectral
//! structure is certified. The source must support the surface preparation,
//! transport and escape conditions. No scattering-heat deposition kernel,
//! angular distribution, collector response or charging feedback is supplied.

use super::{MAX_PHOTOELECTRIC_KNOTS, PhotoelectricKnot, PhotoelectricSourceError, source_value};
use crate::photoelectric::PhotoelectricLaw;
use crate::state_point::ResolvedMaterialStatePoint;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::Dims;

/// Bound on piecewise-constant depth laws, independent of the native 64-claim
/// material-state bound. Shared source keys may explicitly repeat a role.
pub const MAX_PHOTOELECTRIC_DEPTH_BINS: usize = 16;

/// One bin, beginning at the previous bin's endpoint (or zero for the first).
#[derive(Clone, Debug, PartialEq)]
pub struct PhotoelectricDepthBin {
    /// Fraction of the calibrated thickness measured from the escape face.
    /// Endpoints are strictly increasing in (0,1], ending at exactly one.
    pub end_fraction: f64,
    /// Conditional escape probability and mean exit-energy fraction at each
    /// vacuum wavelength. All bins use the same explicitly authored grid.
    pub knots: Vec<PhotoelectricKnot>,
}

/// One source snapshot with independently calibrated local escape laws.
/// Unlike a whole-film effective yield, these probabilities are weighted by
/// the actual absorbed power in each bin, not the bin's geometric volume.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthPhotoelectricMaterial {
    resolved: ResolvedMaterialStatePoint,
    work_function_property: PropertyKey,
    thickness_property: PropertyKey,
    bins: Vec<PhotoelectricDepthBin>,
    values: Vec<Vec<[f64; 2]>>,
    maximum_gap_nm: f64,
    work_function_j: f64,
    thickness_m: f64,
    identity: ContentHash,
}

impl DepthPhotoelectricMaterial {
    /// Bind conditional escape claims to the original complete material state.
    ///
    /// Work function uses joules per electron, film thickness metres, and both
    /// local roles dimensionless scalars in [0,1]. Reusing a source key within
    /// the same role explicitly declares equal values in those bins/wavelengths;
    /// a key cannot serve both escape probability and mean-energy fraction.
    pub fn try_from_resolved(
        resolved: ResolvedMaterialStatePoint,
        work_function_property: PropertyKey,
        thickness_property: PropertyKey,
        bins: &[PhotoelectricDepthBin],
        maximum_gap_nm: f64,
    ) -> Result<Self, PhotoelectricSourceError> {
        let invalid = PhotoelectricSourceError::Invalid;
        if bins.is_empty() || bins.len() > MAX_PHOTOELECTRIC_DEPTH_BINS {
            return Err(invalid("depth photoemission requires 1..=16 bins"));
        }
        if !maximum_gap_nm.is_finite() || maximum_gap_nm <= 0.0 {
            return Err(invalid(
                "spectral interpolation gap must be finite and positive",
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
                "surface work and calibrated film thickness must be positive",
            ));
        }
        let mut identity = DomainHasher::new("org.frankensim.material.photoelectric-depth.v1");
        identity.update(resolved.identity().as_bytes());
        identity.update(&maximum_gap_nm.to_bits().to_le_bytes());
        identity.update(&(bins.len() as u64).to_le_bytes());
        hash_key(&mut identity, &work_function_property);
        hash_key(&mut identity, &thickness_property);
        let mut roles: Vec<(&PropertyKey, usize)> = Vec::new();
        let mut values = Vec::with_capacity(bins.len());
        let mut start = 0.0;
        for (bin_index, bin) in bins.iter().enumerate() {
            if !bin.end_fraction.is_finite()
                || bin.end_fraction <= start
                || bin.end_fraction > 1.0
                || !(2..=MAX_PHOTOELECTRIC_KNOTS).contains(&bin.knots.len())
            {
                return Err(invalid(
                    "depth bins need increasing endpoints and 2..=16 spectral knots",
                ));
            }
            if bin.knots.len() != bins[0].knots.len() {
                return Err(invalid(
                    "depth bins must share the original wavelength grid",
                ));
            }
            identity.update(&bin.end_fraction.to_bits().to_le_bytes());
            identity.update(&(bin.knots.len() as u64).to_le_bytes());
            let mut samples = Vec::with_capacity(bin.knots.len());
            for (index, knot) in bin.knots.iter().enumerate() {
                if !knot.wavelength_nm.is_finite()
                    || knot.wavelength_nm <= 0.0
                    || index > 0 && knot.wavelength_nm <= bin.knots[index - 1].wavelength_nm
                    || index > 0
                        && knot.wavelength_nm - bin.knots[index - 1].wavelength_nm > maximum_gap_nm
                    || bin_index > 0 && knot.wavelength_nm != bins[0].knots[index].wavelength_nm
                {
                    return Err(invalid(
                        "unsupported or inconsistent depth-profile wavelength grid",
                    ));
                }
                let mut sample = [0.0; 2];
                for (role, key) in [
                    &knot.absorbed_quantum_yield_property,
                    &knot.mean_excess_fraction_property,
                ]
                .into_iter()
                .enumerate()
                {
                    match roles.iter().find(|(previous, _)| *previous == key) {
                        Some((_, previous_role)) if *previous_role != role => {
                            return Err(invalid(
                                "one source key cannot describe both escape probability and energy",
                            ));
                        }
                        Some(_) => {}
                        None => roles.push((key, role)),
                    }
                    sample[role] = source_value(&resolved, key, Dims::NONE)?;
                    hash_key(&mut identity, key);
                }
                PhotoelectricLaw::try_new(work_function_j, sample[0], sample[1])?;
                samples.push(sample);
                identity.update(&knot.wavelength_nm.to_bits().to_le_bytes());
            }
            values.push(samples);
            start = bin.end_fraction;
        }
        if start != 1.0 {
            return Err(invalid(
                "depth profile must cover the entire calibrated film",
            ));
        }
        Ok(Self {
            resolved,
            work_function_property,
            thickness_property,
            bins: bins.to_vec(),
            values,
            maximum_gap_nm,
            work_function_j,
            thickness_m,
            identity: identity.finalize(),
        })
    }

    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn work_function_j(&self) -> f64 {
        self.work_function_j
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    pub fn bins(&self) -> &[PhotoelectricDepthBin] {
        &self.bins
    }

    /// Closed bin bounds measured from the selected escape face, not the
    /// optical entrance. Consumers reverse the coordinate for exit emission.
    pub fn bin_bounds_fraction(
        &self,
        bin_index: usize,
    ) -> Result<[f64; 2], PhotoelectricSourceError> {
        let bin = self
            .bins
            .get(bin_index)
            .ok_or(PhotoelectricSourceError::Invalid(
                "depth bin index is out of range",
            ))?;
        let start = if bin_index == 0 {
            0.0
        } else {
            self.bins[bin_index - 1].end_fraction
        };
        Ok([start, bin.end_fraction])
    }

    pub fn wavelength_bounds_nm(&self) -> [f64; 2] {
        let knots = &self.bins[0].knots;
        [knots[0].wavelength_nm, knots[knots.len() - 1].wavelength_nm]
    }

    /// Allocation-free projection of one conditional law, without depth
    /// averaging or extrapolation. Apply it only to absorption in this bin.
    pub fn at_wavelength(
        &self,
        bin_index: usize,
        wavelength_nm: f64,
    ) -> Result<PhotoelectricLaw, PhotoelectricSourceError> {
        let bin = self
            .bins
            .get(bin_index)
            .ok_or(PhotoelectricSourceError::Invalid(
                "depth bin index is out of range",
            ))?;
        let [low, high] = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || wavelength_nm < low || wavelength_nm > high {
            return Err(PhotoelectricSourceError::Invalid(
                "depth-profile wavelength is outside original closed support",
            ));
        }
        let index = bin
            .knots
            .partition_point(|knot| knot.wavelength_nm < wavelength_nm);
        let samples = &self.values[bin_index];
        let values = if bin.knots[index].wavelength_nm == wavelength_nm {
            samples[index]
        } else {
            let fraction = (wavelength_nm - bin.knots[index - 1].wavelength_nm)
                / (bin.knots[index].wavelength_nm - bin.knots[index - 1].wavelength_nm);
            core::array::from_fn(|role| {
                let a = samples[index - 1][role];
                let b = samples[index][role];
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

    /// Atomically refresh the same selected source claims. Fractional depths
    /// stay on the newly calibrated thickness; a consumer must rebind geometry.
    /// This operation advances neither surface charge nor thermal history.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, PhotoelectricSourceError> {
        Self::try_from_resolved(
            self.resolved.at_query_point(card, point)?,
            self.work_function_property.clone(),
            self.thickness_property.clone(),
            &self.bins,
            self.maximum_gap_nm,
        )
    }
}

fn hash_key(identity: &mut DomainHasher, key: &PropertyKey) {
    identity.update(&(key.name().len() as u64).to_le_bytes());
    identity.update(key.name().as_bytes());
    identity.update(&key.quantity().canonical_bytes());
}

#[cfg(test)]
mod tests;
