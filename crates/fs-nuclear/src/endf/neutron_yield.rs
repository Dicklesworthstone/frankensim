//! Evaluated fission-neutron multiplicities from ENDF-6 MF1.
//!
//! MT452, MT456 and MT455 remain total, prompt and delayed yields respectively.
//! Tabulated incident-energy laws retain their original nonnegative values;
//! multiplicities are not probabilities. MT455 with LDG=0 additionally retains
//! the ordered precursor decay constants. Family abundances and outgoing
//! neutron spectra are separate MF5 data and are never inferred from lifetimes.
//!
//! The original bounded reader and reaction-weighted collapse are reused.
//! MF1 yields can come from an unreconstructed, zero-temperature evaluation;
//! their MF3 partner must still satisfy its own processed-data admission.
//! Collapse checks the isotope and evaluation revision before combining them.
//! Polynomial LNU=1 and energy-dependent delayed-family LDG=1 representations
//! refuse explicitly. No processing, data validation or reactor qualification
//! is implied. Format: ENDF-102 (2023), sections 1.2–1.4,
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.

use std::sync::Arc;

use super::reaction_weighted::{
    ReactionWeightedCollapse, ReactionWeightedLimits, TabulatedYield, TabulatedYieldData,
    TabulatedYieldLimits, collapse_reaction_weighted,
};
use super::spectrum::FluxSpectrum;
use super::{
    EndfError, EndfMetadata, EndfReadLimits, EndfSelection, PointwiseCrossSection, Reader, Record,
    budget, read_metadata, read_nonnegative_tab1,
};
use crate::cross_sections::NeutronReaction;
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{FissionError, nonnegative, positive};

/// Distinct source quantities; a total multiplicity is not a prompt yield.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FissionNeutronKind {
    Total,
    Prompt,
    Delayed,
}

impl FissionNeutronKind {
    #[must_use]
    pub const fn mt(self) -> i32 {
        match self {
            Self::Total => 452,
            Self::Prompt => 456,
            Self::Delayed => 455,
        }
    }
}

/// Exact material, nuclide/isomer and MF1 temperature selection. Zero kelvin is
/// admitted for the original evaluation; it is not a broadened cross section.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EndfNeutronYieldSelection {
    pub material_number: u16,
    pub nuclide: NuclideId,
    pub temperature_k: f64,
    pub kind: FissionNeutronKind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfNeutronYieldReadLimits {
    pub records: EndfReadLimits,
    pub max_delayed_families: usize,
}

impl Default for EndfNeutronYieldReadLimits {
    fn default() -> Self {
        Self {
            records: EndfReadLimits::default(),
            max_delayed_families: 64,
        }
    }
}

/// Immutable identified yield law and, for MT455, ordered precursor lifetimes.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionNeutronYield {
    metadata: EndfMetadata,
    kind: FissionNeutronKind,
    curve: TabulatedYield,
    decay_constants_s_inv: Arc<[f64]>,
}

/// A condensation retains the yield's physical meaning and its original MF1
/// metadata, as well as the independently processed MF3 metadata in `collapse`.
#[derive(Clone, Debug, PartialEq)]
pub struct CollapsedFissionNeutronYield {
    pub kind: FissionNeutronKind,
    pub evaluation_metadata: EndfMetadata,
    pub decay_constants_s_inv: Arc<[f64]>,
    /// Exactly one yield row, in the original neutron-group order.
    pub collapse: ReactionWeightedCollapse,
}

impl FissionNeutronYield {
    pub fn parse(
        data_id: String,
        selection: EndfNeutronYieldSelection,
        text: &str,
        limits: EndfNeutronYieldReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || selection.material_number == 0
            || selection.material_number > 9999
            || !nonnegative(selection.temperature_k)
        {
            return Err(FissionError::InvalidInput(
                "identified neutron yield, valid MAT and nonnegative MF1 temperature required",
            )
            .into());
        }
        if text.len() > limits.records.max_bytes {
            return Err(budget("ENDF neutron-yield bytes"));
        }
        let mut reader = Reader {
            lines: text.lines(),
            count: 0,
            limits: limits.records,
        };
        let mut metadata = None;
        let mut candidate = None;
        let mat = i32::from(selection.material_number);
        while let Some(record) = reader.next(&mut checkpoint)? {
            if record.mat != mat {
                metadata = None;
                continue;
            }
            if record.mf == 1 && record.mt == 451 {
                metadata = read_metadata(
                    &mut reader,
                    record,
                    EndfSelection {
                        material_number: selection.material_number,
                        nuclide: selection.nuclide,
                        reaction: NeutronReaction::Fission,
                        temperature_k: selection.temperature_k,
                    },
                    false,
                    &mut checkpoint,
                )?;
            } else if record.mf == 1 && record.mt == selection.kind.mt() {
                if let Some(identity) = &metadata {
                    if candidate.is_some() {
                        return Err(EndfError::AmbiguousSelection);
                    }
                    candidate = Some(read_yield(
                        &mut reader,
                        record,
                        identity.clone(),
                        selection.kind,
                        limits.max_delayed_families,
                        &mut checkpoint,
                    )?);
                }
            }
        }
        let (metadata, curve, decay_constants_s_inv) =
            candidate.ok_or(EndfError::MissingSelection)?;
        let curve = TabulatedYield::new_with_checkpoint(
            TabulatedYieldData {
                data_id,
                energies_j: curve.energies_j,
                values: curve.values,
                interpolation_regions: curve.interpolation_regions,
            },
            TabulatedYieldLimits {
                max_points: limits.records.max_points,
                max_interpolation_regions: limits.records.max_interpolation_regions,
            },
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(Self {
            metadata,
            kind: selection.kind,
            curve,
            decay_constants_s_inv: decay_constants_s_inv.into(),
        })
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        self.curve.data_id()
    }

    #[must_use]
    pub const fn metadata(&self) -> &EndfMetadata {
        &self.metadata
    }

    #[must_use]
    pub const fn kind(&self) -> FissionNeutronKind {
        self.kind
    }

    /// The original law, usable by the existing batch reaction-weighted owner.
    #[must_use]
    pub const fn curve(&self) -> &TabulatedYield {
        &self.curve
    }

    /// Empty for prompt/total data. MT455 order is retained without sorting or
    /// an invented family abundance, equilibrium population or emission spectrum.
    #[must_use]
    pub fn decay_constants_s_inv(&self) -> &[f64] {
        &self.decay_constants_s_inv
    }

    /// Preserve integral(sigma_f * phi * nu) using the original condensation
    /// owner. Different processing temperatures/flags are allowed, but mixing
    /// isotopes, target masses or library/revision identities refuses.
    pub fn collapse(
        &self,
        fission: &PointwiseCrossSection,
        spectrum: &FluxSpectrum,
        groups: Arc<NeutronGroups>,
        limits: ReactionWeightedLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<CollapsedFissionNeutronYield, EndfError> {
        checkpoint()?;
        let a = self.metadata();
        let b = fission.metadata();
        if fission.reaction() != NeutronReaction::Fission
            || a.nuclide != b.nuclide
            || a.atomic_weight_ratio != b.atomic_weight_ratio
            || a.library_id != b.library_id
            || a.library_version != b.library_version
            || a.library_release != b.library_release
            || a.material_revision != b.material_revision
        {
            return Err(FissionError::InvalidInput(
                "neutron yield and fission cross section require the same isotope and evaluation revision",
            )
            .into());
        }
        let collapse = collapse_reaction_weighted(
            fission,
            spectrum,
            std::slice::from_ref(&self.curve),
            groups,
            limits,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(CollapsedFissionNeutronYield {
            kind: self.kind,
            evaluation_metadata: self.metadata.clone(),
            decay_constants_s_inv: Arc::clone(&self.decay_constants_s_inv),
            collapse,
        })
    }
}

fn read_yield(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    metadata: EndfMetadata,
    kind: FissionNeutronKind,
    max_families: usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(EndfMetadata, super::NonnegativeTabulation, Vec<f64>), EndfError> {
    let za = u32::from(metadata.nuclide.atomic_number()) * 1000
        + u32::from(metadata.nuclide.mass_number());
    if head.number(0)? != f64::from(za) || head.number(1)? != metadata.atomic_weight_ratio {
        return Err(head.error("neutron-yield ZA/AWR disagrees with MF1 identity"));
    }
    if head.integer(3)? != 2 {
        return Err(EndfError::Unsupported(
            "fission-neutron yield requires tabulated LNU=2",
        ));
    }
    head.zero_integer(4)?;
    head.zero_integer(5)?;
    let mut decay_constants = Vec::new();
    if kind == FissionNeutronKind::Delayed {
        if head.integer(2)? != 0 {
            return Err(EndfError::Unsupported(
                "energy-dependent delayed families (LDG=1) require a separate model",
            ));
        }
        let list = reader.section(head.mat, 1, head.mt, checkpoint)?;
        zero_header(list)?;
        list.zero_integer(5)?;
        let families = list.count(4)?;
        if families == 0 {
            return Err(list.error("delayed-neutron families must not be empty"));
        }
        if families > max_families
            || families.div_ceil(6) > reader.limits.max_records.saturating_sub(reader.count)
        {
            return Err(budget("ENDF delayed-neutron families"));
        }
        decay_constants
            .try_reserve_exact(families)
            .map_err(|_| budget("ENDF delayed-neutron memory"))?;
        for _ in 0..families.div_ceil(6) {
            let row = reader.section(head.mat, 1, head.mt, checkpoint)?;
            for field in 0..6 {
                let value = row.number(field)?;
                if decay_constants.len() == families {
                    if value != 0.0 {
                        return Err(row.error("nonzero delayed-neutron LIST padding"));
                    }
                } else if !positive(value) {
                    return Err(row.error("positive delayed-neutron decay constants required"));
                } else {
                    decay_constants.push(value);
                }
            }
        }
    } else {
        head.zero_integer(2)?;
    }
    let tab = reader.section(head.mat, 1, head.mt, checkpoint)?;
    zero_header(tab)?;
    let curve = read_nonnegative_tab1(reader, tab, 1.0, checkpoint)?;
    reader.send(head.mat, 1, checkpoint)?;
    Ok((metadata, curve, decay_constants))
}

fn zero_header(record: Record<'_>) -> Result<(), EndfError> {
    if record.number(0)? != 0.0 || record.number(1)? != 0.0 {
        return Err(record.error("nonzero reserved neutron-yield header field"));
    }
    record.zero_integer(2)?;
    record.zero_integer(3)?;
    Ok(())
}

#[cfg(test)]
mod tests;
