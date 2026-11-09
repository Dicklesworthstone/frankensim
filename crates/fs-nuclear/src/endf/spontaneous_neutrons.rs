//! Spontaneous-fission neutron data from the radioactive-decay sublibrary.
//!
//! ENDF-102 (2023), §§0.3.2.4, 1.2–1.4 and 5.4.1, places constant total,
//! prompt and delayed neutron multiplicities in NSUB=4 MF1/MT452,456,455;
//! prompt and delayed energy distributions are MF5/MT18,455 respectively:
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.
//! These are per fission, not per transformation of a mixed-mode parent.
//! Neither prompt/delayed subtraction nor a second decay-branch factor occurs.
//!
//! The opaque types share the original bounded ENDF reader but cannot enter
//! incident-neutron collapse. Source values, component probabilities, energy
//! grids and interpolation are retained without normalization repair. Missing
//! sections are not zero emission. Separate source binding must establish the
//! physical normalization and compatible decay inventory/energy routing.

use std::sync::Arc;

use super::neutron_yield::FissionNeutronKind;
use super::{EndfError, EndfMetadata, EndfReadLimits, Reader, Record, budget};
use crate::irradiation::NuclideId;
use crate::{FissionError, nonnegative, positive};

mod spectrum;
pub use spectrum::{
    SpontaneousFissionNeutronSpectrum, SpontaneousFissionShapePoint,
    SpontaneousFissionSpectrumComponent, SpontaneousFissionSpectrumLaw,
};

/// Both the MAT and full parent isotope/isomer identity must agree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfSpontaneousNeutronSelection {
    pub material_number: u16,
    pub nuclide: NuclideId,
}

/// Explicit record, allocation and source-work bounds. Spectrum point and
/// interpolation counts include repeated original incident-energy rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpontaneousNeutronReadLimits {
    pub records: EndfReadLimits,
    pub max_delayed_families: usize,
    pub max_components: usize,
    pub max_incident_points: usize,
    pub max_spectrum_points: usize,
    pub max_interpolation_regions: usize,
}

impl Default for SpontaneousNeutronReadLimits {
    fn default() -> Self {
        Self {
            records: EndfReadLimits::default(),
            max_delayed_families: 64,
            max_components: 64,
            max_incident_points: 4096,
            max_spectrum_points: 1_000_000,
            max_interpolation_regions: 8192,
        }
    }
}

/// One explicitly supplied LNU=1, NC=1 neutron multiplicity. Delayed lifetimes
/// retain original family order; they do not imply abundances or a spectrum.
#[derive(Clone, Debug, PartialEq)]
pub struct SpontaneousFissionNeutronYield {
    data_id: Arc<str>,
    metadata: EndfMetadata,
    kind: FissionNeutronKind,
    yield_per_fission: f64,
    decay_constants_s_inv: Arc<[f64]>,
}

impl SpontaneousFissionNeutronYield {
    pub fn parse(
        data_id: String,
        selection: EndfSpontaneousNeutronSelection,
        kind: FissionNeutronKind,
        text: &str,
        limits: SpontaneousNeutronReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        let (metadata, yield_per_fission, decay_constants_s_inv) = selected(
            &data_id,
            selection,
            text,
            limits,
            1,
            kind.mt(),
            &mut checkpoint,
            |reader, head, metadata, checkpoint| {
                match_identity(head, &metadata)?;
                if head.integer(3)? != 1 {
                    return Err(EndfError::Unsupported(
                        "spontaneous neutron yields require constant LNU=1, NC=1",
                    ));
                }
                head.zero_integer(4)?;
                head.zero_integer(5)?;
                let mut families = Vec::new();
                if kind == FissionNeutronKind::Delayed {
                    if head.integer(2)? != 0 {
                        return Err(EndfError::Unsupported(
                            "spontaneous delayed lifetimes require energy-independent LDG=0",
                        ));
                    }
                    let list = reader.section(head.mat, 1, head.mt, checkpoint)?;
                    zero_prefix(list)?;
                    list.zero_integer(5)?;
                    let count = list.count(4)?;
                    if count > limits.max_delayed_families {
                        return Err(budget("spontaneous delayed-neutron families"));
                    }
                    families = values(reader, list, count, checkpoint)?;
                    if !families.iter().copied().all(positive) {
                        return Err(list.error("positive delayed-family decay constants required"));
                    }
                } else {
                    head.zero_integer(2)?;
                }
                let list = reader.section(head.mat, 1, head.mt, checkpoint)?;
                zero_prefix(list)?;
                list.zero_integer(5)?;
                if list.count(4)? != 1 {
                    return Err(EndfError::Unsupported(
                        "spontaneous neutron yields require one constant polynomial coefficient",
                    ));
                }
                let yield_ = values(reader, list, 1, checkpoint)?[0];
                if !nonnegative(yield_)
                    || (kind == FissionNeutronKind::Delayed && yield_ > 0.0 && families.is_empty())
                {
                    return Err(list
                        .error("nonnegative yield and explicit positive-yield families required"));
                }
                reader.send(head.mat, 1, checkpoint)?;
                Ok((metadata, yield_, families))
            },
        )?;
        Ok(Self {
            data_id: data_id.into(),
            metadata,
            kind,
            yield_per_fission,
            decay_constants_s_inv: decay_constants_s_inv.into(),
        })
    }

    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    pub fn metadata(&self) -> &EndfMetadata {
        &self.metadata
    }
    pub fn kind(&self) -> FissionNeutronKind {
        self.kind
    }
    pub fn yield_per_fission(&self) -> f64 {
        self.yield_per_fission
    }
    pub fn decay_constants_s_inv(&self) -> &[f64] {
        &self.decay_constants_s_inv
    }
}

fn selected<'a, T, C>(
    data_id: &str,
    selection: EndfSpontaneousNeutronSelection,
    text: &'a str,
    limits: SpontaneousNeutronReadLimits,
    mf: i32,
    mt: i32,
    checkpoint: &mut C,
    mut read: impl FnMut(&mut Reader<'a>, Record<'a>, EndfMetadata, &mut C) -> Result<T, EndfError>,
) -> Result<T, EndfError>
where
    C: FnMut() -> Result<(), FissionError>,
{
    checkpoint()?;
    if data_id.trim().is_empty() || !(1..=9999).contains(&selection.material_number) {
        return Err(FissionError::InvalidInput(
            "identified spontaneous-neutron source and valid MAT required",
        )
        .into());
    }
    if text.len() > limits.records.max_bytes {
        return Err(budget("spontaneous-neutron source bytes"));
    }
    let mut reader = Reader {
        lines: text.lines(),
        count: 0,
        limits: limits.records,
    };
    let mut metadata = None;
    let mut candidate = None;
    while let Some(record) = reader.next(checkpoint)? {
        if record.mat != i32::from(selection.material_number) {
            metadata = None;
        } else if (record.mf, record.mt) == (1, 451) {
            metadata = Some(read_metadata(&mut reader, record, selection, checkpoint)?);
        } else if (record.mf, record.mt) == (mf, mt) {
            if candidate.is_some() {
                return Err(EndfError::AmbiguousSelection);
            }
            let identity = metadata.as_ref().ok_or_else(|| {
                record.error(
                "MF1/MT451 radioactive-decay identity must precede the selected neutron section",
            )
            })?;
            candidate = Some(read(&mut reader, record, identity.clone(), checkpoint)?);
        }
    }
    checkpoint()?;
    candidate.ok_or(EndfError::MissingSelection)
}

fn read_metadata(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    selection: EndfSpontaneousNeutronSelection,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<EndfMetadata, EndfError> {
    let second = reader.section(head.mat, 1, 451, checkpoint)?;
    let third = reader.section(head.mat, 1, 451, checkpoint)?;
    let fourth = reader.section(head.mat, 1, 451, checkpoint)?;
    let texts = fourth.count(4)?;
    let directories = fourth.count(5)?;
    require_records(
        reader,
        texts
            .checked_add(directories)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| budget("spontaneous-neutron metadata records"))?,
    )?;
    for _ in 0..texts {
        reader.section(head.mat, 1, 451, checkpoint)?;
    }
    for _ in 0..directories {
        let row = reader.section(head.mat, 1, 451, checkpoint)?;
        for field in 2..6 {
            row.count(field)?;
        }
    }
    reader.send(head.mat, 1, checkpoint)?;
    if second.integer(5)? != 6 || third.integer(4)? != 4 {
        return Err(EndfError::Unsupported(
            "ENDF-6 radioactive-decay sublibrary NSUB=4 required",
        ));
    }
    let za = head.number(0)?;
    if za.fract() != 0.0 || !(1001.0..=999_999.0).contains(&za) {
        return Err(head.error("integer isotope ZA required"));
    }
    let za = za as u32;
    let isomer =
        u16::try_from(second.integer(3)?).map_err(|_| second.error("invalid parent LISO"))?;
    let nuclide = NuclideId::new((za / 1000) as u16, (za % 1000) as u16, isomer)?;
    if nuclide != selection.nuclide || second.integer(2)? < i32::from(isomer) {
        return Err(head.error("selected spontaneous-neutron isotope/isomer disagrees"));
    }
    let atomic_weight_ratio = head.number(1)?;
    let temperature_k = fourth.number(0)?;
    if !positive(atomic_weight_ratio)
        || !nonnegative(second.number(0)?)
        || ![0.0, 1.0].contains(&second.number(1)?)
        || third.number(0)? != 0.0
        || !nonnegative(third.number(1)?)
        || !nonnegative(temperature_k)
        || !(-1..=2).contains(&head.integer(2)?)
        || !(0..=1).contains(&head.integer(3)?)
    {
        return Err(head.error("invalid spontaneous-neutron source identity fields"));
    }
    second.zero_integer(4)?;
    third.zero_integer(3)?;
    fourth.zero_integer(3)?;
    if fourth.number(1)? != 0.0 || fourth.integer(2)? < 0 {
        return Err(fourth.error("invalid metadata derivation/reserved field"));
    }
    Ok(EndfMetadata {
        material_number: selection.material_number,
        nuclide,
        temperature_k,
        atomic_weight_ratio,
        library_id: head.integer(4)?,
        library_version: third.integer(5)?,
        library_release: third.integer(2)?,
        material_revision: head.integer(5)?,
        resonance_flag: head.integer(2)?,
        derivation_flag: fourth.integer(2)?,
    })
}

fn match_identity(head: Record<'_>, metadata: &EndfMetadata) -> Result<(), EndfError> {
    let za = u32::from(metadata.nuclide.atomic_number()) * 1000
        + u32::from(metadata.nuclide.mass_number());
    if head.number(0)? != f64::from(za) || head.number(1)? != metadata.atomic_weight_ratio {
        return Err(head.error("neutron-section ZA/AWR disagrees with MF1 identity"));
    }
    Ok(())
}

fn zero_prefix(record: Record<'_>) -> Result<(), EndfError> {
    if record.number(0)? != 0.0 || record.number(1)? != 0.0 {
        return Err(record.error("nonzero reserved neutron header field"));
    }
    record.zero_integer(2)?;
    record.zero_integer(3)?;
    Ok(())
}

fn require_records(reader: &Reader<'_>, count: usize) -> Result<(), EndfError> {
    if count > reader.limits.max_records.saturating_sub(reader.count) {
        return Err(budget("spontaneous-neutron source records"));
    }
    Ok(())
}

fn reserved<T>(count: usize, name: &'static str) -> Result<Vec<T>, EndfError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| budget(name))?;
    Ok(values)
}

fn values(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    count: usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Vec<f64>, EndfError> {
    require_records(reader, count.div_ceil(6))?;
    let mut values = reserved(count, "spontaneous-neutron scalar memory")?;
    while values.len() < count {
        let record = reader.section(head.mat, head.mf, head.mt, checkpoint)?;
        for field in 0..6 {
            let value = record.number(field)?;
            if values.len() < count {
                values.push(value);
            } else if value != 0.0 {
                return Err(record.error("nonzero neutron LIST padding"));
            }
        }
    }
    Ok(values)
}

#[cfg(test)]
pub(in crate::endf) mod tests;
