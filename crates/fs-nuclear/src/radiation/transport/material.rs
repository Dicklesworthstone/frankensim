//! Original-material mass coefficients for the planar photon transport owner.
//!
//! These are supplied, group-collapsed energy coefficients in m²/kg. Total
//! photon extinction is separate from local heat and outgoing photon energy.
//! In particular, neither attenuation nor a mass energy-absorption table alone
//! determines the redistribution kernel. Local heat assumes the declared
//! charged-particle deposition model; secondary electron transport is not
//! inferred. Coefficients are fixed for their declared material applicability.

use std::sync::Arc;

use super::{PhotonCollision, PhotonSlabCell, PhotonSlabLimits, validate_coefficients};
use crate::endf::decay::spectra::DecayRadiation;
use crate::radiation::{ROW_CLOSURE_TOLERANCE, RadiationError, Sum, finite, nonnegative, product};

/// Source-bound temperature interpolation of the complete photon energy kernel.
pub mod temperature;

/// One incident photon group's supplied mass energy coefficients [m²/kg].
#[derive(Clone, Debug, PartialEq)]
pub struct PhotonMassCollision {
    /// Removal of incident-group photon energy by all modeled collisions.
    pub total_m2_kg: f64,
    /// Energy transferred to locally thermalized charged particles/material.
    pub heat_m2_kg: f64,
    /// Outgoing photon ENERGY per incident energy, times the collision mass
    /// coefficient, in the original group order. These are not number yields.
    /// Same-group elastic scattering is allowed; upscatter is not supplied.
    pub redistribution_m2_kg: Vec<f64>,
}

/// Immutable coefficient table tied to an original material identity and grid.
/// Clones share the data. No material mixture, spectrum or density is invented.
#[derive(Clone, Debug)]
pub struct PhotonMassMaterial(Arc<MaterialData>);

/// One original constituent and its fraction of the mixture's total MASS.
/// These are neither volume fractions nor normalized photon-interaction weights.
#[derive(Clone, Debug)]
pub struct PhotonMassConstituent {
    pub mass_fraction: f64,
    pub material: PhotonMassMaterial,
}

/// Retained physical assumptions and source tables of a homogeneous mixture.
#[derive(Clone, Debug)]
pub struct PhotonMassMixture {
    /// Caller-supplied basis for independent photon interactions, compatible
    /// group collapsing and the applicability of every constituent's complete
    /// outgoing-photon/local-deposition kernel IN this target mixture.
    pub kernel_applicability: String,
    /// Original authored order, fractions and immutable source material tables.
    pub constituents: Vec<PhotonMassConstituent>,
}

#[derive(Debug)]
struct MaterialData {
    data_id: String,
    material_id: u64,
    radiation: DecayRadiation,
    edges_j: Vec<f64>,
    groups: Vec<PhotonMassCollision>,
    mixture: Option<PhotonMassMixture>,
}

impl PhotonMassMaterial {
    /// Admit a supplied constant-material photon kernel. All groups retain
    /// their original coefficients; closure admission never renormalizes them.
    pub fn new(
        data_id: String,
        material_id: u64,
        radiation: DecayRadiation,
        group_edges_j: Vec<f64>,
        groups: Vec<PhotonMassCollision>,
        limits: PhotonSlabLimits,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<Self, RadiationError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || !matches!(radiation, DecayRadiation::Gamma | DecayRadiation::XRay)
            || groups.is_empty()
            || group_edges_j.len().checked_sub(1) != Some(groups.len())
        {
            return Err(RadiationError::InvalidInput(
                "identified photon mass data and matching nonempty group grid required",
            ));
        }
        let count = groups.len();
        let values = count
            .checked_add(2)
            .and_then(|width| count.checked_mul(width))
            .and_then(|rows| rows.checked_add(group_edges_j.len()))
            .ok_or(RadiationError::Budget("photon mass coefficient storage"))?;
        if count > limits.max_groups || values > limits.max_values || values > limits.max_work {
            return Err(RadiationError::Budget(
                "photon mass coefficient storage/work",
            ));
        }
        for (index, &edge) in group_edges_j.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            if !nonnegative(edge) || (index > 0 && edge <= group_edges_j[index - 1]) {
                return Err(RadiationError::InvalidInput(
                    "finite nonnegative ascending photon energy edges required",
                ));
            }
        }
        for (incident, row) in groups.iter().enumerate() {
            checkpoint()?;
            // The closure and downscatter rules are dimensionless, so mass
            // and macroscopic tables share the transport owner's admission.
            validate_coefficients(
                row.total_m2_kg,
                row.heat_m2_kg,
                &row.redistribution_m2_kg,
                incident,
                count,
            )?;
        }
        checkpoint()?;
        Ok(Self(Arc::new(MaterialData {
            data_id,
            material_id,
            radiation,
            edges_j: group_edges_j,
            groups,
            mixture: None,
        })))
    }

    /// Prepare the complete mass kernel of an explicitly homogeneous mixture.
    ///
    /// Each coefficient is `sum_i(w_i * coefficient_i)`, using compensated
    /// summation; positive product underflow refuses. Fractions must close one
    /// within [`ROW_CLOSURE_TOLERANCE`]; neither composition nor kernels are normalized.
    /// All constituents, including zero-weight sources, must have the same
    /// photon kind and exact energy grid. The result uses the target material
    /// identity; its later density must be the actual MIXTURE density.
    ///
    /// Mass attenuation additivity follows NIST section 2, equation (6):
    /// <https://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html>.
    /// This does NOT establish additivity of mass energy-absorption coefficients
    /// (NIST section 3). `kernel_applicability` must justify each supplied
    /// constituent heat/redistribution kernel in the target mixture, including
    /// secondary-particle deposition and common within-group spectral weighting.
    /// Matching grids alone cannot establish those physical assumptions. No
    /// chemical absorption-edge structure or electron stopping law is inferred.
    /// Resolved layers must use slab geometry instead of this homogenization.
    ///
    /// Preparation bounds cover immediate source kernels, result, fractions and
    /// coefficient visits before allocation. Original source objects are shared;
    /// repeated transport/heating uses the prepared table without mixture work.
    pub fn from_mass_fractions(
        data_id: String,
        material_id: u64,
        constituents: Vec<PhotonMassConstituent>,
        kernel_applicability: String,
        limits: PhotonSlabLimits,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<Self, RadiationError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || kernel_applicability.trim().is_empty()
            || constituents.is_empty()
        {
            return Err(RadiationError::InvalidInput(
                "identified photon mixture, constituent mass fractions and kernel applicability required",
            ));
        }
        let first = &constituents[0].material;
        let count = first.groups().len();
        let values = count
            .checked_add(2)
            .and_then(|width| count.checked_mul(width))
            .and_then(|rows| rows.checked_add(first.group_edges_j().len()))
            .and_then(|kernel| {
                constituents
                    .len()
                    .checked_add(1)
                    .and_then(|tables| kernel.checked_mul(tables))
            })
            .and_then(|tables| tables.checked_add(constituents.len()))
            .and_then(|tables| tables.checked_add(2))
            .ok_or(RadiationError::Budget("photon mixture storage/work"))?;
        if count > limits.max_groups || values > limits.max_values || values > limits.max_work {
            return Err(RadiationError::Budget("photon mixture storage/work"));
        }
        let mut fractions = Sum::default();
        for constituent in &constituents {
            checkpoint()?;
            if !nonnegative(constituent.mass_fraction)
                || constituent.mass_fraction > 1.0
                || constituent.material.radiation() != first.radiation()
                || constituent.material.group_edges_j() != first.group_edges_j()
            {
                return Err(RadiationError::InvalidInput(
                    "photon mixture needs fractions in [0,1] and one photon kind/group grid",
                ));
            }
            fractions.add(constituent.mass_fraction)?;
        }
        fractions.add(-1.0)?;
        if fractions.hi.abs() + fractions.lo.abs() > ROW_CLOSURE_TOLERANCE {
            return Err(RadiationError::InvalidInput(
                "photon constituent mass fractions must close one without normalization",
            ));
        }
        let mut groups = Vec::with_capacity(count);
        for incident in 0..count {
            checkpoint()?;
            let total_m2_kg = weighted_coefficient(
                &constituents,
                incident,
                |row| row.total_m2_kg,
                &mut checkpoint,
            )?;
            let heat_m2_kg = weighted_coefficient(
                &constituents,
                incident,
                |row| row.heat_m2_kg,
                &mut checkpoint,
            )?;
            let mut redistribution_m2_kg = Vec::with_capacity(count);
            for outgoing in 0..count {
                redistribution_m2_kg.push(weighted_coefficient(
                    &constituents,
                    incident,
                    |row| row.redistribution_m2_kg[outgoing],
                    &mut checkpoint,
                )?);
            }
            groups.push(PhotonMassCollision {
                total_m2_kg,
                heat_m2_kg,
                redistribution_m2_kg,
            });
        }
        // Reuse the original table's energy closure and downscatter admission.
        let mut result = Self::new(
            data_id,
            material_id,
            first.radiation(),
            first.group_edges_j().to_vec(),
            groups,
            limits,
            &mut checkpoint,
        )?;
        Arc::get_mut(&mut result.0)
            .expect("newly admitted unshared photon mixture")
            .mixture = Some(PhotonMassMixture {
            kernel_applicability,
            constituents,
        });
        Ok(result)
    }

    pub fn data_id(&self) -> &str {
        &self.0.data_id
    }

    /// Identity of the physical material for which the supplied data applies.
    pub fn material_id(&self) -> u64 {
        self.0.material_id
    }

    pub fn radiation(&self) -> DecayRadiation {
        self.0.radiation
    }

    pub fn group_edges_j(&self) -> &[f64] {
        &self.0.edges_j
    }

    pub fn groups(&self) -> &[PhotonMassCollision] {
        &self.0.groups
    }

    /// Original composition and its physical applicability, when prepared by
    /// [`Self::from_mass_fractions`]. Directly supplied tables have no mixture.
    pub fn mixture(&self) -> Option<&PhotonMassMixture> {
        self.0.mixture.as_ref()
    }

    /// Scale by the explicitly supplied original reference density [kg/m³].
    /// The thermal adapter obtains it from stored mass / geometric volume.
    /// Positive underflow/overflow is refused, and no partial cell is returned.
    pub fn macroscopic_cell(
        &self,
        cell_id: u64,
        thickness_m: f64,
        density_kg_m3: f64,
        mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<PhotonSlabCell, RadiationError> {
        checkpoint()?;
        if !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
            || !thickness_m.is_finite()
            || thickness_m <= 0.0
        {
            return Err(RadiationError::InvalidInput(
                "positive finite reference density and photon layer thickness required",
            ));
        }
        let mut groups = Vec::with_capacity(self.0.groups.len());
        for row in &self.0.groups {
            checkpoint()?;
            let mut redistribution_per_m = Vec::with_capacity(self.0.groups.len());
            for (outgoing, &value) in row.redistribution_m2_kg.iter().enumerate() {
                if outgoing % 64 == 0 {
                    checkpoint()?;
                }
                redistribution_per_m.push(product(density_kg_m3, value)?.0);
            }
            groups.push(PhotonCollision {
                total_per_m: product(density_kg_m3, row.total_m2_kg)?.0,
                heat_per_m: product(density_kg_m3, row.heat_m2_kg)?.0,
                redistribution_per_m,
            });
        }
        checkpoint()?;
        Ok(PhotonSlabCell {
            cell_id,
            thickness_m,
            groups,
        })
    }
}

fn weighted_coefficient(
    constituents: &[PhotonMassConstituent],
    incident: usize,
    select: impl Fn(&PhotonMassCollision) -> f64,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<f64, RadiationError> {
    let mut sum = Sum::default();
    for (index, constituent) in constituents.iter().enumerate() {
        if index % 64 == 0 {
            checkpoint()?;
        }
        let value = select(&constituent.material.groups()[incident]);
        sum.add(product(constituent.mass_fraction, value)?.0)?;
    }
    finite(sum.hi + sum.lo)
}

#[cfg(test)]
mod tests;
