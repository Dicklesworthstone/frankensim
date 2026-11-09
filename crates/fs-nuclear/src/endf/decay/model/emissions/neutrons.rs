//! Evaluated parent-wide neutron spectra on an authored transport energy grid.
//!
//! Each original spectral bin contributes its independently integrated particle
//! count and first energy moment. Its mean is projected onto two neighboring
//! representative neutron energies, preserving both moments without clipping,
//! extrapolation or another decay-branch multiplier. The original spectrum and
//! model remain attached. This is a declared spectral condensation; it does not
//! resolve within-bin reaction thresholds or supply angular transport.
//!
//! ENDF-102 (2023), §8.4 and Eq. (8.1), defines spectra per parent transformation,
//! summed over all decay modes, including beta-delayed neutron emission:
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.
//! The original heavy-particle radiation account must pay for the transported
//! kinetic energy. Deposited heat, escape and electromagnetic energy cannot be
//! repurposed as a neutron source. No source normalization or data reconciliation
//! changes the evaluated neutron count or energy.
//!
//! The resulting yields enter the original coupled decay/transport operator;
//! this preparation performs no time integration or interval-mean source fit.

use std::sync::Arc;

use super::{DecayEmissionCoverage, PreparedDecayEmission, visit};
use crate::diffusion::decay::{DecayNeutronEmission, DecayNeutronError};
use crate::endf::EndfError;
use crate::endf::decay::spectra::DecayRadiation;
use crate::endf::spectrum::PositiveIntegral;
use crate::fusion_neutrons::NeutronKineticModel;
use crate::irradiation::NuclideId;
use crate::{Accumulator, FissionError, nonnegative, sum};

/// Original moments and signed differences introduced by spectral condensation.
/// These are per parent transformation, not additional particles or energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayNeutronProjection {
    pub particles_per_decay: f64,
    pub kinetic_energy_j_per_decay: f64,
    pub particle_projection_residual: f64,
    pub energy_projection_residual_j: f64,
}

/// A neutron emission law retaining the complete evaluated source and the exact
/// kinetic-energy model used to project it. Clones share immutable source data.
#[derive(Clone, Debug)]
pub struct EvaluatedDecayNeutrons {
    source: Arc<PreparedDecayEmission>,
    kinetic_model: Arc<NeutronKineticModel>,
    emission: Arc<DecayNeutronEmission>,
    confirmed_nonemitting_parents: Arc<[usize]>,
    projection: Arc<[DecayNeutronProjection]>,
}

impl EvaluatedDecayNeutrons {
    pub fn source(&self) -> &Arc<PreparedDecayEmission> {
        &self.source
    }
    pub fn kinetic_model(&self) -> &Arc<NeutronKineticModel> {
        &self.kinetic_model
    }
    pub fn emission(&self) -> &Arc<DecayNeutronEmission> {
        &self.emission
    }
    /// Exact original nuclide indices, retained in caller declaration order.
    /// The caller supplies the independent evidence; a missing ENDF spectrum
    /// does not establish zero neutron emission.
    pub fn confirmed_nonemitting_parents(&self) -> &[usize] {
        &self.confirmed_nonemitting_parents
    }
    /// One original-order entry per nuclide, including stable and confirmed zeros.
    pub fn projection(&self) -> &[DecayNeutronProjection] {
        &self.projection
    }
}

/// Refusal before publishing a partially bound neutron source.
#[derive(Clone, Debug, PartialEq)]
pub enum DecayNeutronBindingError {
    Invalid(&'static str),
    /// No neutron spectrum and no independently supplied zero-emission declaration.
    Unreported {
        parent: usize,
        nuclide: NuclideId,
    },
    /// A spectral mean falls outside the authored representative energy range.
    OutsideRepresentativeEnergies {
        parent: usize,
        source_group: usize,
        mean_j: f64,
    },
    /// Original heavy-particle outward energy cannot pay for this neutron source.
    RadiationBudget {
        parent: usize,
        kinetic_energy_j: f64,
        available_energy_j: f64,
    },
    /// Condensation failed the fixed floating-point count or energy account.
    Projection {
        parent: usize,
        particles: f64,
        energy_j: f64,
    },
    Source(EndfError),
    Coupling(DecayNeutronError),
}

impl From<EndfError> for DecayNeutronBindingError {
    fn from(error: EndfError) -> Self {
        Self::Source(error)
    }
}
impl From<FissionError> for DecayNeutronBindingError {
    fn from(error: FissionError) -> Self {
        Self::Source(error.into())
    }
}
impl From<DecayNeutronError> for DecayNeutronBindingError {
    fn from(error: DecayNeutronError) -> Self {
        Self::Coupling(error)
    }
}
impl core::fmt::Display for DecayNeutronBindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "evaluated decay neutron binding refused: {self:?}")
    }
}
impl std::error::Error for DecayNeutronBindingError {}

impl PreparedDecayEmission {
    /// Condense evaluated neutron bin means onto an explicit kinetic model.
    ///
    /// `confirmed_nonemitting_parents` names every unstable nuclide with missing
    /// neutron spectra which the caller independently confirms emits no neutrons.
    /// Indices must be unique, in range, and refer only to `Unreported` entries.
    /// The supplied `data_id` identifies this projection and those declarations;
    /// neither a nuclide name nor an omitted record proves zero emission.
    /// Stable laws and explicitly reported zero spectra need no declaration.
    /// A positive-probability supplied neutron-emitting mode cannot be declared
    /// nonemitting merely because its spectrum is missing.
    ///
    /// Source and target energy grids may differ. The original preparation's
    /// allocation/work limits also bound this operation. Positive subnormal
    /// coefficients refuse here because the prepared physical source may later
    /// amplify them by arbitrarily large expected inventories.
    pub fn bind_neutron_group_means(
        self: &Arc<Self>,
        data_id: String,
        model: Arc<NeutronKineticModel>,
        confirmed_nonemitting_parents: &[usize],
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<EvaluatedDecayNeutrons, DecayNeutronBindingError> {
        checkpoint()?;
        if self.radiation != DecayRadiation::Neutron || data_id.trim().is_empty() {
            return Err(DecayNeutronBindingError::Invalid(
                "identified neutron spectrum and projection required",
            ));
        }
        let parents = self.coverage.len();
        let energies = model.group_energies_j();
        let groups = energies.len();
        if groups > self.limits.max_groups
            || parents
                .checked_mul(groups + 4)
                .is_none_or(|count| count > self.limits.max_group_values)
        {
            return Err(FissionError::Budget("evaluated neutron projection storage").into());
        }
        let mut confirmed = vec![false; parents];
        let mut work = 0;
        for &parent in confirmed_nonemitting_parents {
            visit(&mut work, self.limits, &mut checkpoint)?;
            if parent >= parents
                || confirmed[parent]
                || self.coverage[parent] != DecayEmissionCoverage::Unreported
            {
                return Err(DecayNeutronBindingError::Invalid(
                    "nonemitting declarations must uniquely identify unreported parents",
                ));
            }
            for mode in self.assembly.evaluations()[parent].modes() {
                visit(&mut work, self.limits, &mut checkpoint)?;
                if mode.probability > 0.0 && mode.rtyp.to_string().contains('5') {
                    return Err(DecayNeutronBindingError::Invalid(
                        "a supplied neutron-emitting mode contradicts the nonemitting declaration",
                    ));
                }
            }
            confirmed[parent] = true;
        }

        let mut yields = Vec::with_capacity(parents);
        let mut projection = Vec::with_capacity(parents);
        for (parent, coverage) in self.coverage.iter().enumerate() {
            visit(&mut work, self.limits, &mut checkpoint)?;
            if *coverage == DecayEmissionCoverage::Unreported && !confirmed[parent] {
                return Err(DecayNeutronBindingError::Unreported {
                    parent,
                    nuclide: self.assembly.nuclides()[parent],
                });
            }
            let mut row = vec![Accumulator::new(0.0); groups];
            let mut original_count = PositiveIntegral::new();
            let mut original_energy = PositiveIntegral::new();
            for (source_group, (&count, &energy)) in self
                .particles_per_decay(parent)
                .unwrap()
                .iter()
                .zip(self.energy_per_decay_j(parent).unwrap())
                .enumerate()
            {
                visit(&mut work, self.limits, &mut checkpoint)?;
                if count == 0.0 {
                    continue;
                }
                if count.is_subnormal() || energy.is_subnormal() {
                    return Err(FissionError::Arithmetic.into());
                }
                original_count.add_product(1.0, count, 1.0, 1.0)?;
                original_energy.add_product(1.0, energy, 1.0, 1.0)?;
                let mean = energy / count;
                if !mean.is_finite() || mean < energies[0] || mean > energies[groups - 1] {
                    return Err(DecayNeutronBindingError::OutsideRepresentativeEnergies {
                        parent,
                        source_group,
                        mean_j: mean,
                    });
                }
                let upper = energies.partition_point(|&value| value < mean);
                if energies[upper] == mean {
                    row[upper].add(count)?;
                } else {
                    let lower = upper - 1;
                    let span = energies[upper] - energies[lower];
                    let low = count * ((energies[upper] - mean) / span);
                    let high = count * ((mean - energies[lower]) / span);
                    if !(low >= f64::MIN_POSITIVE && high >= f64::MIN_POSITIVE)
                        || !low.is_finite()
                        || !high.is_finite()
                    {
                        return Err(FissionError::Arithmetic.into());
                    }
                    row[lower].add(low)?;
                    row[upper].add(high)?;
                }
            }
            let original_count = original_count.value()?;
            let original_energy = original_energy.value()?;
            let row: Vec<_> = row.into_iter().map(Accumulator::value).collect();
            let mut actual_count = PositiveIntegral::new();
            let mut actual_energy = PositiveIntegral::new();
            for (&count, &energy) in row.iter().zip(energies) {
                visit(&mut work, self.limits, &mut checkpoint)?;
                actual_count.add_product(1.0, count, 1.0, 1.0)?;
                actual_energy.add_product(1.0, count, energy, 1.0)?;
            }
            let actual_energy = actual_energy.value()?;
            let number_residual = sum(&[actual_count.value()?, -original_count]);
            let energy_residual = sum(&[actual_energy, -original_energy]);
            if !nonnegative(number_residual.abs())
                || !nonnegative(energy_residual.abs())
                || number_residual.abs() > 256.0 * f64::EPSILON * original_count
                || energy_residual.abs() > 256.0 * f64::EPSILON * original_energy
            {
                return Err(DecayNeutronBindingError::Projection {
                    parent,
                    particles: number_residual,
                    energy_j: energy_residual,
                });
            }
            // Full spectral moments are transferred. The original caller must
            // already have reserved that energy in the heavy-particle outward
            // channel; gamma radiation and local heat are not substitute credit.
            let evaluation = &self.assembly.evaluations()[parent];
            let heavy = evaluation.average_energies()[2].value_j.ok_or(
                DecayNeutronBindingError::Invalid("unknown original heavy-particle energy"),
            )?;
            let mut available = PositiveIntegral::new();
            available.add_product(
                1.0,
                heavy,
                self.assembly.routing()[parent].heavy_particles.radiation,
                1.0,
            )?;
            let available = available.value()?;
            if original_energy > available || actual_energy > available {
                return Err(DecayNeutronBindingError::RadiationBudget {
                    parent,
                    kinetic_energy_j: original_energy.max(actual_energy),
                    available_energy_j: available,
                });
            }
            yields.push(row);
            projection.push(DecayNeutronProjection {
                particles_per_decay: original_count,
                kinetic_energy_j_per_decay: original_energy,
                particle_projection_residual: number_residual,
                energy_projection_residual_j: energy_residual,
            });
        }
        checkpoint()?;
        let emission = Arc::new(DecayNeutronEmission::new(
            data_id,
            Arc::clone(self.assembly.model()),
            Arc::clone(model.diffusion().groups()),
            yields,
        )?);
        checkpoint()?;
        Ok(EvaluatedDecayNeutrons {
            source: Arc::clone(self),
            kinetic_model: model,
            emission,
            confirmed_nonemitting_parents: confirmed_nonemitting_parents.into(),
            projection: projection.into(),
        })
    }
}

#[cfg(test)]
mod tests;
