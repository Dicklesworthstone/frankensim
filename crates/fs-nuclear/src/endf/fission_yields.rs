//! Independent fission-product yields from supplied ENDF-6 MF8/MT454 data.
//!
//! Independent yields describe products before their subsequent radioactive
//! decays. Cumulative MT459 yields are deliberately not accepted here: using
//! them as births in a full decay network would count daughter production twice.
//! Multiplicities and their reported uncertainties are retained without
//! normalization. No missing product or nuclide evaluation is assumed stable.
//!
//! MF1/MT451 supplies the parent isomer and neutron-induced-yield sublibrary
//! identity; the MF8 header alone cannot establish either. Histogram and linear
//! interpolation are supported. Every energy row must explicitly contain the
//! same ordered product identities. ENDF's LE=0 declares energy-independent
//! yields; its single incident energy remains recorded as source metadata.
//! Format: ENDF-102, section 8.3,
//! <https://nds.iaea.org/public/endf/endf-manual.pdf>.

use std::sync::Arc;

use super::{
    EV_J, EndfError, EndfInterpolation, EndfReadLimits, Reader, Record, budget, convert_si,
};
use crate::irradiation::NuclideId;
use crate::{FissionError, nonnegative, positive};

/// Reaction-weighted condensation and disjoint-group assembly into burnup.
pub mod model;

/// Temperature-dependent removal, product production and prompt event energy.
pub mod temperature;

/// Both the selected MAT and the full parent identity must match the source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfFissionYieldSelection {
    pub material_number: u16,
    pub nuclide: NuclideId,
}

/// Source declarations, not a data-quality or applicability certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionYieldMetadata {
    pub material_number: u16,
    pub nuclide: NuclideId,
    pub atomic_weight_ratio: f64,
    pub library_id: i32,
    pub library_version: i32,
    pub library_release: i32,
    pub material_revision: i32,
}

/// Bounds are checked before source-controlled allocation. Product values count
/// all (yield, uncertainty) pairs across all incident-energy rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfFissionYieldReadLimits {
    pub max_bytes: usize,
    pub max_records: usize,
    pub max_energy_points: usize,
    pub max_products_per_energy: usize,
    pub max_product_values: usize,
}

impl Default for EndfFissionYieldReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 128 * 1024 * 1024,
            max_records: 1_600_000,
            max_energy_points: 4096,
            max_products_per_energy: 8192,
            max_product_values: 1_000_000,
        }
    }
}

/// Expected immediate products per fission and the source's one-sigma uncertainty.
/// Uncertainties are retained data; they are not independently propagated errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FissionYieldValue {
    pub yield_per_fission: f64,
    pub uncertainty: f64,
}

/// One incident-energy row, in the evaluation's common product order.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionYieldPoint {
    pub energy_j: f64,
    /// The source's law between the previous point and this point; none at E1.
    pub interpolation_from_previous: Option<EndfInterpolation>,
    pub values: Arc<[FissionYieldValue]>,
}

/// Immutable independent yields. Product identities are shared by every row.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionYieldEvaluation {
    data_id: Arc<str>,
    metadata: FissionYieldMetadata,
    products: Arc<[NuclideId]>,
    points: Arc<[FissionYieldPoint]>,
}

impl FissionYieldEvaluation {
    /// Read identified neutron-induced independent yields. This accepts neither
    /// spontaneous/photo-fission sublibraries nor a cumulative MT459 substitute.
    pub fn parse(
        data_id: String,
        selection: EndfFissionYieldSelection,
        text: &str,
        limits: EndfFissionYieldReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || selection.material_number == 0
            || selection.material_number > 9999
        {
            return Err(invalid(
                "identified fission-yield source and MAT in 1..=9999 required",
            ));
        }
        if text.len() > limits.max_bytes {
            return Err(budget("ENDF fission-yield bytes"));
        }
        let mut reader = Reader {
            lines: text.lines(),
            count: 0,
            limits: EndfReadLimits {
                max_bytes: limits.max_bytes,
                max_records: limits.max_records,
                max_points: 0,
                max_interpolation_regions: 0,
            },
        };
        let mut metadata = None;
        let mut candidate = None;
        while let Some(record) = reader.next(&mut checkpoint)? {
            if record.mat != i32::from(selection.material_number) {
                metadata = None;
            } else if (record.mf, record.mt) == (1, 451) {
                metadata = Some(read_metadata(
                    &mut reader,
                    record,
                    selection,
                    &mut checkpoint,
                )?);
            } else if (record.mf, record.mt) == (8, 454) {
                if candidate.is_some() {
                    return Err(EndfError::AmbiguousSelection);
                }
                let identity = metadata.as_ref().ok_or_else(|| {
                    record
                        .error("MF1/MT451 parent and neutron-yield identity must precede MF8/MT454")
                })?;
                candidate = Some(read_evaluation(
                    &mut reader,
                    record,
                    identity.clone(),
                    limits,
                    &mut checkpoint,
                )?);
            }
        }
        checkpoint()?;
        let mut evaluation = candidate.ok_or(EndfError::MissingSelection)?;
        evaluation.data_id = data_id.into();
        Ok(evaluation)
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }

    #[must_use]
    pub fn metadata(&self) -> &FissionYieldMetadata {
        &self.metadata
    }

    #[must_use]
    pub fn products(&self) -> &[NuclideId] {
        &self.products
    }

    #[must_use]
    pub fn points(&self) -> &[FissionYieldPoint] {
        &self.points
    }

    /// LE=0 explicitly declares no incident-energy dependence.
    #[must_use]
    pub fn is_energy_independent(&self) -> bool {
        self.points.len() == 1
    }
}

fn read_metadata(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    selection: EndfFissionYieldSelection,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<FissionYieldMetadata, EndfError> {
    let second = reader.section(head.mat, 1, 451, checkpoint)?;
    let third = reader.section(head.mat, 1, 451, checkpoint)?;
    let fourth = reader.section(head.mat, 1, 451, checkpoint)?;
    let texts = fourth.count(4)?;
    let directories = fourth.count(5)?;
    let remaining = texts
        .checked_add(directories)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| budget("ENDF fission-yield metadata records"))?;
    require_records(reader, remaining)?;
    for _ in 0..texts {
        reader.section(head.mat, 1, 451, checkpoint)?;
    }
    for _ in 0..directories {
        let record = reader.section(head.mat, 1, 451, checkpoint)?;
        for field in 2..6 {
            record.count(field)?;
        }
    }
    reader.send(head.mat, 1, checkpoint)?;
    if second.integer(5)? != 6 || third.integer(4)? != 11 {
        return Err(EndfError::Unsupported(
            "ENDF-6 neutron-induced fission-yield sublibrary NSUB=11 required",
        ));
    }
    let isomer = u16::try_from(second.integer(3)?)
        .map_err(|_| second.error("invalid parent isomer LISO"))?;
    let nuclide = nuclide(head.number(0)?, f64::from(isomer), head)?;
    if nuclide != selection.nuclide || second.integer(2)? < i32::from(isomer) {
        return Err(head.error("selected fission-yield isotope/isomer or level disagrees"));
    }
    let atomic_weight_ratio = head.number(1)?;
    if !positive(atomic_weight_ratio)
        || !nonnegative(second.number(0)?)
        || !nonnegative(fourth.number(0)?)
    {
        return Err(head.error("invalid parent mass, excitation or temperature"));
    }
    second.zero_integer(4)?;
    third.zero_integer(3)?;
    fourth.zero_integer(3)?;
    if fourth.number(1)? != 0.0 {
        return Err(fourth.error("nonzero reserved metadata field"));
    }
    Ok(FissionYieldMetadata {
        material_number: selection.material_number,
        nuclide,
        atomic_weight_ratio,
        library_id: head.integer(4)?,
        library_version: third.integer(5)?,
        library_release: third.integer(2)?,
        material_revision: head.integer(5)?,
    })
}

fn read_evaluation(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    metadata: FissionYieldMetadata,
    limits: EndfFissionYieldReadLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<FissionYieldEvaluation, EndfError> {
    let za = u32::from(metadata.nuclide.atomic_number()) * 1000
        + u32::from(metadata.nuclide.mass_number());
    if head.number(0)? != f64::from(za) || head.number(1)? != metadata.atomic_weight_ratio {
        return Err(head.error("fission-yield ZA/AWR disagrees with MF1 identity"));
    }
    for field in 3..6 {
        head.zero_integer(field)?;
    }
    let energies = head.count(2)?;
    if energies == 0 {
        return Err(head.error("at least one yield energy row required"));
    }
    if energies > limits.max_energy_points {
        return Err(budget("fission-yield energy points"));
    }
    require_records(reader, energies.saturating_mul(2).saturating_add(1))?;
    let mut points = reserved(energies, "fission-yield row memory")?;
    let mut products = Vec::new();
    let mut total_values = 0_usize;
    for row in 0..energies {
        let list = reader.section(head.mat, 8, 454, checkpoint)?;
        if list.number(1)? != 0.0 {
            return Err(list.error("nonzero reserved yield field"));
        }
        list.zero_integer(3)?;
        let energy_j = convert_si(list.number(0)?, EV_J)?;
        if !nonnegative(energy_j)
            || points
                .last()
                .is_some_and(|last: &FissionYieldPoint| last.energy_j >= energy_j)
        {
            return Err(list.error("finite nonnegative ascending yield energies required"));
        }
        let interpolation_from_previous = if row == 0 {
            if list.count(2)? != energies - 1 {
                return Err(list.error("LE disagrees with HEAD"));
            }
            None
        } else {
            Some(match list.integer(2)? {
                1 => EndfInterpolation::Histogram,
                2 => EndfInterpolation::Linear,
                _ => {
                    return Err(EndfError::Unsupported(
                        "histogram or linear fission yields required",
                    ));
                }
            })
        };
        let count = list.count(5)?;
        if count == 0 {
            return Err(list.error("explicit nonempty product list required"));
        }
        if count > limits.max_products_per_energy {
            return Err(budget("fission products per energy"));
        }
        total_values = total_values
            .checked_add(count)
            .ok_or_else(|| budget("fission-yield values"))?;
        if total_values > limits.max_product_values {
            return Err(budget("fission-yield values"));
        }
        let scalar_count = count
            .checked_mul(4)
            .ok_or_else(|| budget("fission-yield LIST size"))?;
        if list.count(4)? != scalar_count {
            return Err(list.error("yield NN must equal 4*NFP"));
        }
        require_records(reader, scalar_count.div_ceil(6))?;
        let mut scalars = reserved(scalar_count, "fission-yield LIST memory")?;
        while scalars.len() < scalar_count {
            let record = reader.section(head.mat, 8, 454, checkpoint)?;
            for field in 0..6 {
                let value = record.number(field)?;
                if scalars.len() < scalar_count {
                    scalars.push(value);
                } else if value != 0.0 {
                    return Err(record.error("nonzero yield LIST padding"));
                }
            }
        }
        if row == 0 {
            products = reserved(count, "fission-product identity memory")?;
        } else if count != products.len() {
            return Err(list.error("product set changes between energies"));
        }
        let mut values = reserved(count, "fission-yield value memory")?;
        for (index, entry) in scalars.chunks_exact(4).enumerate() {
            checkpoint()?;
            let id = nuclide(entry[0], entry[1], list)?;
            if row == 0 {
                if products.last().is_some_and(|last| *last >= id) {
                    return Err(list.error("products must be unique and ordered by ZA then isomer"));
                }
                products.push(id);
            } else if products[index] != id {
                return Err(list.error("product identities differ between yield energies"));
            }
            if !nonnegative(entry[2]) || !nonnegative(entry[3]) {
                return Err(list.error("finite nonnegative yields and uncertainties required"));
            }
            values.push(FissionYieldValue {
                yield_per_fission: entry[2],
                uncertainty: entry[3],
            });
        }
        points.push(FissionYieldPoint {
            energy_j,
            interpolation_from_previous,
            values: values.into(),
        });
    }
    reader.send(head.mat, 8, checkpoint)?;
    Ok(FissionYieldEvaluation {
        data_id: "".into(),
        metadata,
        products: products.into(),
        points: points.into(),
    })
}

fn nuclide(za: f64, state: f64, record: Record<'_>) -> Result<NuclideId, EndfError> {
    if za.fract() != 0.0
        || !(1001.0..=999_999.0).contains(&za)
        || state.fract() != 0.0
        || !(0.0..=f64::from(u16::MAX)).contains(&state)
    {
        return Err(record.error("integer isotope ZA and nonnegative isomer required"));
    }
    let za = za as u32;
    Ok(NuclideId::new(
        (za / 1000) as u16,
        (za % 1000) as u16,
        state as u16,
    )?)
}

fn require_records(reader: &Reader<'_>, required: usize) -> Result<(), EndfError> {
    if required > reader.limits.max_records.saturating_sub(reader.count) {
        return Err(budget("ENDF fission-yield records"));
    }
    Ok(())
}

fn reserved<T>(count: usize, detail: &'static str) -> Result<Vec<T>, EndfError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| budget(detail))?;
    Ok(values)
}

fn invalid(detail: &'static str) -> EndfError {
    FissionError::InvalidInput(detail).into()
}

#[cfg(test)]
mod tests;
