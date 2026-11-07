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
use crate::radiation::{RadiationError, nonnegative, product};

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

#[derive(Debug)]
struct MaterialData {
    data_id: String,
    material_id: u64,
    radiation: DecayRadiation,
    edges_j: Vec<f64>,
    groups: Vec<PhotonMassCollision>,
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
        })))
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

#[cfg(test)]
mod tests;
