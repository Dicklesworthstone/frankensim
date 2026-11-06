//! Bounded import of processed neutron ENDF-6/PENDF MF3 reaction curves.
//!
//! This deliberately narrow reader accepts total fission (MT18) and radiative
//! capture (MT102), with MF1/MT451 identification, a positive processed target
//! temperature, and LRP=0 or 2. LRP=1 is rejected: those MF3 values are resonance
//! backgrounds, not complete reconstructed cross sections. Natural-element
//! evaluations, unprocessed zero-temperature data and interpolation laws other
//! than histogram or linear-linear require another processor.
//!
//! Energies and Q values are converted from eV to J and cross sections from
//! barns to m². Q values are metadata only: neither is assigned to prompt heat.
//! Group collapse integrates the admitted pointwise law with explicitly flat
//! flux density per unit energy inside each group. It does not reconstruct
//! resonances, Doppler broaden, self-shield, infer a spectrum or certify data.
//! Format: ENDF-102 (2023), §§0.5–0.6, 1.1 and 3.2–3.3,
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.

use std::sync::Arc;

/// Explicit incident spectra and reaction-rate-preserving group collapse.
pub mod spectrum;

/// Evaluated prompt, delayed and total fission-neutron multiplicities.
pub mod neutron_yield;

/// Independent fission products connected to the original burnup inventory.
pub mod fission_yields;

/// Evaluated half-lives, decay modes and aggregate emitted energies.
pub mod decay;

/// Dimensionless reaction observables condensed with their production rates.
pub mod reaction_weighted;

use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{Accumulator, FissionError, nonnegative, positive};

const EV_J: f64 = 1.602_176_634e-19;
const BARN_M2: f64 = 1.0e-28;

/// All four identities must match. Repeated MAT values at different processed
/// temperatures are allowed; two sections matching the entire selection are not.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EndfSelection {
    pub material_number: u16,
    pub nuclide: NuclideId,
    pub reaction: NeutronReaction,
    pub temperature_k: f64,
}

/// Bounds include unrelated records scanned while selecting a section.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfReadLimits {
    pub max_bytes: usize,
    pub max_records: usize,
    pub max_points: usize,
    pub max_interpolation_regions: usize,
}

impl Default for EndfReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 128 * 1024 * 1024,
            max_records: 1_600_000,
            max_points: 1_000_000,
            max_interpolation_regions: 4096,
        }
    }
}

/// Segment visits bound integration work, including zero-width discontinuities.
/// A table assembly shares this budget across all supplied temperature curves.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfCollapseLimits {
    pub max_neutron_groups: usize,
    pub max_segment_visits: usize,
}

impl Default for EndfCollapseLimits {
    fn default() -> Self {
        Self {
            max_neutron_groups: 4096,
            max_segment_visits: 4_000_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EndfError {
    Format { line: usize, reason: &'static str },
    Unsupported(&'static str),
    MissingSelection,
    AmbiguousSelection,
    Physics(FissionError),
}

impl From<FissionError> for EndfError {
    fn from(value: FissionError) -> Self {
        Self::Physics(value)
    }
}

impl core::fmt::Display for EndfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Format { line, reason } => write!(f, "ENDF record {line}: {reason}"),
            Self::Unsupported(reason) => write!(f, "unsupported ENDF input: {reason}"),
            Self::MissingSelection => {
                f.write_str("identified ENDF reaction and temperature not found")
            }
            Self::AmbiguousSelection => f.write_str("multiple ENDF sections match the selection"),
            Self::Physics(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for EndfError {}

/// Retained MF1/MT451 identity. These are declarations from the source file,
/// not independent evidence that its evaluation or processing is accurate.
#[derive(Clone, Debug, PartialEq)]
pub struct EndfMetadata {
    pub material_number: u16,
    pub nuclide: NuclideId,
    pub temperature_k: f64,
    pub atomic_weight_ratio: f64,
    pub library_id: i32,
    pub library_version: i32,
    pub library_release: i32,
    pub material_revision: i32,
    pub resonance_flag: i32,
    pub derivation_flag: i32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndfInterpolation {
    Histogram,
    Linear,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfInterpolationRegion {
    /// ENDF NBT: one-based index of the region's final point.
    pub end_point: usize,
    pub interpolation: EndfInterpolation,
}

/// Immutable processed curve. Equal adjacent energies represent an ENDF
/// discontinuity: the first value is the left limit, the second the right.
/// Evaluation exactly at a discontinuity returns its right-hand value.
#[derive(Clone, Debug, PartialEq)]
pub struct PointwiseCrossSection {
    data_id: Arc<str>,
    metadata: EndfMetadata,
    reaction: NeutronReaction,
    energies_j: Arc<[f64]>,
    cross_sections_m2: Arc<[f64]>,
    interpolation_regions: Arc<[EndfInterpolationRegion]>,
    mass_difference_q_j: f64,
    reaction_q_j: f64,
}

impl PointwiseCrossSection {
    /// Read a supplied source without network access or built-in isotope data.
    /// MF1 temperature must equal the explicit selection; no snapping occurs.
    pub fn parse(
        data_id: String,
        selection: EndfSelection,
        text: &str,
        limits: EndfReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || selection.material_number == 0
            || selection.material_number > 9999
            || !positive(selection.temperature_k)
        {
            return Err(FissionError::InvalidInput(
                "source identity, MAT and positive target temperature required",
            )
            .into());
        }
        if text.len() > limits.max_bytes {
            return Err(budget("ENDF bytes"));
        }
        let mut reader = Reader {
            lines: text.lines(),
            count: 0,
            limits,
        };
        let mut metadata = None;
        let mut candidate = None;
        let mat = i32::from(selection.material_number);
        let mt = reaction_mt(selection.reaction);
        while let Some(record) = reader.next(&mut checkpoint)? {
            if record.mat != mat {
                metadata = None;
                continue;
            }
            if record.mf == 1 && record.mt == 451 {
                metadata = read_metadata(&mut reader, record, selection, true, &mut checkpoint)?;
            } else if record.mf == 3 && record.mt == mt {
                if let Some(identity) = &metadata {
                    if candidate.is_some() {
                        return Err(EndfError::AmbiguousSelection);
                    }
                    candidate = Some(read_curve(
                        &mut reader,
                        record,
                        identity.clone(),
                        selection.reaction,
                        &mut checkpoint,
                    )?);
                }
            }
        }
        checkpoint()?;
        let mut curve = candidate.ok_or(EndfError::MissingSelection)?;
        curve.data_id = data_id.into();
        Ok(curve)
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    #[must_use]
    pub const fn metadata(&self) -> &EndfMetadata {
        &self.metadata
    }
    #[must_use]
    pub const fn reaction(&self) -> NeutronReaction {
        self.reaction
    }
    #[must_use]
    pub fn energies_j(&self) -> &[f64] {
        &self.energies_j
    }
    #[must_use]
    pub fn cross_sections_m2(&self) -> &[f64] {
        &self.cross_sections_m2
    }
    #[must_use]
    pub fn interpolation_regions(&self) -> &[EndfInterpolationRegion] {
        &self.interpolation_regions
    }
    /// Mass-difference Q [J]; not a prompt deposited-energy coefficient.
    #[must_use]
    pub const fn mass_difference_q_j(&self) -> f64 {
        self.mass_difference_q_j
    }
    /// Reaction Q [J]; not a prompt deposited-energy coefficient.
    #[must_use]
    pub const fn reaction_q_j(&self) -> f64 {
        self.reaction_q_j
    }

    /// Binary-search evaluation of the supplied law, with no extrapolation.
    pub fn cross_section_at_energy(&self, energy_j: f64) -> Result<f64, EndfError> {
        self.check_domain(energy_j, energy_j)?;
        let left = self
            .energies_j
            .partition_point(|&energy| energy <= energy_j)
            - 1;
        if self.energies_j[left] == energy_j {
            return Ok(self.cross_sections_m2[left]);
        }
        let region = self
            .interpolation_regions
            .partition_point(|r| r.end_point <= left + 1);
        self.segment_value(left, region, energy_j)
    }

    /// Group mean σ = ∫σ(E)dE / ΔE. Each histogram/linear segment is integrated
    /// analytically, up to floating-point rounding; there is no midpoint sample.
    pub fn collapse_flat_in_energy(
        &self,
        groups: &NeutronGroups,
        limits: EndfCollapseLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<f64>, EndfError> {
        self.collapse(groups, limits, &mut 0, &mut checkpoint)
    }

    fn check_domain(&self, lower: f64, upper: f64) -> Result<(), EndfError> {
        if !nonnegative(lower)
            || !nonnegative(upper)
            || lower < self.energies_j[0]
            || upper > self.energies_j[self.energies_j.len() - 1]
        {
            return Err(
                FissionError::InvalidInput("energy outside the supplied ENDF curve").into(),
            );
        }
        Ok(())
    }

    fn segment_value(&self, segment: usize, region: usize, energy: f64) -> Result<f64, EndfError> {
        let lower = self.cross_sections_m2[segment];
        if self.interpolation_regions[region].interpolation == EndfInterpolation::Histogram {
            return Ok(lower);
        }
        let upper = self.cross_sections_m2[segment + 1];
        let width = self.energies_j[segment + 1] - self.energies_j[segment];
        let lower_weight = (self.energies_j[segment + 1] - energy) / width;
        let upper_weight = (energy - self.energies_j[segment]) / width;
        if (lower_weight == 0.0 && energy != self.energies_j[segment + 1])
            || (upper_weight == 0.0 && energy != self.energies_j[segment])
        {
            return Err(FissionError::Arithmetic.into());
        }
        // Compute both distances directly: 1 - weight can erase a physically
        // relevant tail when the opposite endpoint's coefficient is large.
        convex(lower, upper, lower_weight, upper_weight)
    }

    fn collapse(
        &self,
        groups: &NeutronGroups,
        limits: EndfCollapseLimits,
        visits: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<f64>, EndfError> {
        checkpoint()?;
        if groups.len() > limits.max_neutron_groups {
            return Err(budget("ENDF collapse groups"));
        }
        let edges = groups.edges_j();
        self.check_domain(edges[0], edges[edges.len() - 1])?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(groups.len())
            .map_err(|_| budget("ENDF collapse memory"))?;
        let mut segment = self.energies_j.partition_point(|&e| e <= edges[0]) - 1;
        let mut region = self
            .interpolation_regions
            .partition_point(|r| r.end_point <= segment + 1);
        for bin in edges.windows(2) {
            checkpoint()?;
            let width = bin[1] - bin[0];
            let mut left = bin[0];
            let mut area = Accumulator::new(0.0);
            while left < bin[1] {
                checkpoint()?;
                if *visits >= limits.max_segment_visits {
                    return Err(budget("ENDF collapse segment visits"));
                }
                *visits += 1;
                // A bin boundary can coincide with both sides of a jump.
                if self.energies_j[segment + 1] <= left {
                    segment += 1;
                    while self.interpolation_regions[region].end_point <= segment + 1 {
                        region += 1;
                    }
                    continue;
                }
                let right = bin[1].min(self.energies_j[segment + 1]);
                let mean = convex(
                    self.segment_value(segment, region, left)?,
                    self.segment_value(segment, region, right)?,
                    0.5,
                    0.5,
                )?;
                // Normalize before multiplying; σ * ΔE could overflow even
                // when the physically needed group mean is representable.
                let weight = (right - left) / width;
                let contribution = mean * weight;
                if !nonnegative(contribution)
                    || weight == 0.0
                    || (mean > 0.0 && contribution == 0.0)
                {
                    return Err(FissionError::Arithmetic.into());
                }
                area.add(contribution)?;
                left = right;
            }
            let mean = area.value();
            if !nonnegative(mean) {
                return Err(FissionError::Arithmetic.into());
            }
            result.push(mean);
        }
        checkpoint()?;
        Ok(result)
    }
}

/// Build the common reaction table from strictly ascending processed source
/// temperatures. Each curve retains its own source ID; `data_id` identifies
/// the caller's assembled table and declared flat-in-energy group collapse.
/// Nuclide, reaction and ENDF library/version/release must agree across rows.
pub fn collapse_temperature_table(
    data_id: String,
    curves: &[PointwiseCrossSection],
    groups: Arc<NeutronGroups>,
    collapse_limits: EndfCollapseLimits,
    table_limits: CrossSectionLimits,
    mut checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<CrossSectionTable, EndfError> {
    checkpoint()?;
    let first = curves.first().ok_or(FissionError::InvalidInput(
        "at least one processed ENDF curve required",
    ))?;
    if data_id.trim().is_empty() {
        return Err(FissionError::InvalidInput("assembled source identity required").into());
    }
    if curves.len() > table_limits.max_temperatures {
        return Err(budget("cross-section temperatures"));
    }
    if groups.len() > table_limits.max_neutron_groups {
        return Err(budget("neutron groups"));
    }
    let count = curves
        .len()
        .checked_mul(groups.len())
        .ok_or_else(|| budget("cross-section table size"))?;
    if count > table_limits.max_values {
        return Err(budget("cross-section table values"));
    }
    let mut temperatures = Vec::new();
    let mut rows = Vec::new();
    temperatures
        .try_reserve_exact(curves.len())
        .map_err(|_| budget("ENDF temperature memory"))?;
    rows.try_reserve_exact(curves.len())
        .map_err(|_| budget("ENDF row memory"))?;
    let mut visits = 0;
    let mut previous_temperature = 0.0;
    for curve in curves {
        checkpoint()?;
        let a = &first.metadata;
        let b = &curve.metadata;
        if b.nuclide != a.nuclide
            || curve.reaction != first.reaction
            || (
                b.library_id,
                b.library_version,
                b.library_release,
                b.material_revision,
            ) != (
                a.library_id,
                a.library_version,
                a.library_release,
                a.material_revision,
            )
            || b.temperature_k <= previous_temperature
        {
            return Err(FissionError::InvalidInput(
                "matching ENDF evaluation/reaction and ascending temperatures required",
            )
            .into());
        }
        previous_temperature = b.temperature_k;
        temperatures.push(b.temperature_k);
        rows.push(curve.collapse(&groups, collapse_limits, &mut visits, &mut checkpoint)?);
    }
    CrossSectionTable::new_with_checkpoint(
        CrossSectionData {
            data_id,
            nuclide: first.metadata.nuclide,
            reaction: first.reaction,
            groups,
            temperatures_k: temperatures,
            cross_sections_m2: rows,
        },
        table_limits,
        checkpoint,
    )
    .map_err(Into::into)
}

fn convex(lower: f64, upper: f64, lower_weight: f64, upper_weight: f64) -> Result<f64, EndfError> {
    if !nonnegative(lower_weight)
        || lower_weight > 1.0
        || !nonnegative(upper_weight)
        || upper_weight > 1.0
    {
        return Err(FissionError::Arithmetic.into());
    }
    let value = if lower <= upper {
        lower + (upper - lower) * upper_weight
    } else {
        upper + (lower - upper) * lower_weight
    };
    if !nonnegative(value)
        || (value == 0.0
            && ((lower > 0.0 && lower_weight > 0.0) || (upper > 0.0 && upper_weight > 0.0)))
    {
        return Err(FissionError::Arithmetic.into());
    }
    Ok(value)
}

fn reaction_mt(reaction: NeutronReaction) -> i32 {
    match reaction {
        NeutronReaction::Fission => 18,
        NeutronReaction::Capture => 102,
    }
}

fn budget(reason: &'static str) -> EndfError {
    FissionError::Budget(reason).into()
}

fn convert_si(value: f64, scale: f64) -> Result<f64, EndfError> {
    let converted = value * scale;
    if !converted.is_finite() || (value != 0.0 && converted == 0.0) {
        return Err(FissionError::Arithmetic.into());
    }
    Ok(converted)
}

#[derive(Clone, Copy)]
struct Record<'a> {
    text: &'a str,
    line: usize,
    mat: i32,
    mf: i32,
    mt: i32,
}

impl Record<'_> {
    fn error(self, reason: &'static str) -> EndfError {
        EndfError::Format {
            line: self.line,
            reason,
        }
    }
    fn integer(self, field: usize) -> Result<i32, EndfError> {
        integer(&self.text[field * 11..(field + 1) * 11])
            .ok_or_else(|| self.error("invalid integer field"))
    }
    fn number(self, field: usize) -> Result<f64, EndfError> {
        number(&self.text[field * 11..(field + 1) * 11])
            .ok_or_else(|| self.error("invalid or unrepresentable real field"))
    }
    fn count(self, field: usize) -> Result<usize, EndfError> {
        usize::try_from(self.integer(field)?)
            .map_err(|_| self.error("negative record or point count"))
    }
    fn zero_integer(self, field: usize) -> Result<(), EndfError> {
        if self.integer(field)? != 0 {
            return Err(self.error("nonzero reserved integer field"));
        }
        Ok(())
    }
}

struct Reader<'a> {
    lines: std::str::Lines<'a>,
    count: usize,
    limits: EndfReadLimits,
}

impl<'a> Reader<'a> {
    fn next(
        &mut self,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Option<Record<'a>>, EndfError> {
        checkpoint()?;
        let Some(text) = self.lines.next() else {
            return Ok(None);
        };
        if self.count >= self.limits.max_records {
            return Err(budget("ENDF records"));
        }
        self.count += 1;
        let invalid = |reason| EndfError::Format {
            line: self.count,
            reason,
        };
        if !(75..=80).contains(&text.len()) || !text.bytes().all(|b| (b' '..=b'~').contains(&b)) {
            return Err(invalid("75–80 printable ASCII columns required"));
        }
        let mat = integer(&text[66..70]).ok_or_else(|| invalid("invalid MAT"))?;
        let mf = integer(&text[70..72]).ok_or_else(|| invalid("invalid MF"))?;
        let mt = integer(&text[72..75]).ok_or_else(|| invalid("invalid MT"))?;
        if !(-1..=9999).contains(&mat)
            || !(0..=99).contains(&mf)
            || !(0..=999).contains(&mt)
            || integer(&text[75..]).is_none()
        {
            return Err(invalid("invalid record identifiers"));
        }
        Ok(Some(Record {
            text,
            line: self.count,
            mat,
            mf,
            mt,
        }))
    }

    fn section(
        &mut self,
        mat: i32,
        mf: i32,
        mt: i32,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Record<'a>, EndfError> {
        let record = self.next(checkpoint)?.ok_or(EndfError::Format {
            line: self.count + 1,
            reason: "truncated section",
        })?;
        if (record.mat, record.mf, record.mt) != (mat, mf, mt) {
            return Err(record.error("MAT/MF/MT changed inside a selected section"));
        }
        Ok(record)
    }

    fn send(
        &mut self,
        mat: i32,
        mf: i32,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), EndfError> {
        let record = self.section(mat, mf, 0, checkpoint)?;
        if record.number(0)? != 0.0 || record.number(1)? != 0.0 {
            return Err(record.error("nonzero SEND payload"));
        }
        for field in 2..6 {
            record.zero_integer(field)?;
        }
        Ok(())
    }
}

fn read_metadata(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    selection: EndfSelection,
    require_reconstructed: bool,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Option<EndfMetadata>, EndfError> {
    let second = reader.section(head.mat, 1, 451, checkpoint)?;
    let third = reader.section(head.mat, 1, 451, checkpoint)?;
    let fourth = reader.section(head.mat, 1, 451, checkpoint)?;
    let text_count = fourth.count(4)?;
    let directory_count = fourth.count(5)?;
    let records = text_count
        .checked_add(directory_count)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| budget("ENDF metadata records"))?;
    if records > reader.limits.max_records.saturating_sub(reader.count) {
        return Err(budget("ENDF metadata records"));
    }
    for _ in 0..text_count {
        reader.section(head.mat, 1, 451, checkpoint)?;
    }
    for _ in 0..directory_count {
        let entry = reader.section(head.mat, 1, 451, checkpoint)?;
        // Only the directory's fixed integer fields are interpreted here.
        for field in 2..6 {
            entry.count(field)?;
        }
    }
    reader.send(head.mat, 1, checkpoint)?;
    let temperature = fourth.number(0)?;
    if temperature != selection.temperature_k {
        return Ok(None);
    }
    if second.integer(5)? != 6 || third.integer(4)? != 10 {
        return Err(EndfError::Unsupported(
            "ENDF-6 incident-neutron sublibrary required",
        ));
    }
    let resonance_flag = head.integer(2)?;
    if !(0..=2).contains(&resonance_flag) {
        return Err(head.error("invalid resonance flag"));
    }
    if require_reconstructed && resonance_flag != 0 && resonance_flag != 2 {
        return Err(EndfError::Unsupported(
            "MF3 must contain complete reconstructed cross sections (LRP=0 or 2)",
        ));
    }
    let za = head.number(0)?;
    if za.fract() != 0.0 || za < 1001.0 || za > 999_999.0 {
        return Err(head.error("isotope ZA required"));
    }
    let za = za as u32;
    let isomer = u16::try_from(second.integer(3)?).map_err(|_| second.error("invalid LISO"))?;
    let nuclide = NuclideId::new((za / 1000) as u16, (za % 1000) as u16, isomer)?;
    if nuclide != selection.nuclide {
        return Err(head.error("selected MAT/temperature has a different isotope or isomer"));
    }
    if second.integer(2)? < i32::from(isomer) {
        return Err(second.error("LIS must be at least LISO"));
    }
    let atomic_weight_ratio = head.number(1)?;
    if !positive(atomic_weight_ratio)
        || !nonnegative(second.number(0)?)
        || !positive(third.number(0)?)
        || !positive(third.number(1)?)
    {
        return Err(head.error("invalid target/projectile mass or excitation"));
    }
    if ![0.0, 1.0].contains(&second.number(1)?) {
        return Err(second.error("invalid target stability flag"));
    }
    let fission_flag = head.integer(3)?;
    if !(0..=1).contains(&fission_flag)
        || (selection.reaction == NeutronReaction::Fission && fission_flag != 1)
    {
        return Err(head.error("fission flag disagrees with the selected reaction"));
    }
    second.zero_integer(4)?;
    third.zero_integer(3)?;
    fourth.zero_integer(3)?;
    if fourth.number(1)? != 0.0 {
        return Err(fourth.error("nonzero reserved real field"));
    }
    let derivation_flag = fourth.integer(2)?;
    if derivation_flag < 0 {
        return Err(fourth.error("invalid derivation flag"));
    }
    Ok(Some(EndfMetadata {
        material_number: selection.material_number,
        nuclide,
        temperature_k: temperature,
        atomic_weight_ratio,
        library_id: head.integer(4)?,
        library_version: third.integer(5)?,
        library_release: third.integer(2)?,
        material_revision: head.integer(5)?,
        resonance_flag,
        derivation_flag,
    }))
}

fn read_curve(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    metadata: EndfMetadata,
    reaction: NeutronReaction,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<PointwiseCrossSection, EndfError> {
    let za = u32::from(metadata.nuclide.atomic_number()) * 1000
        + u32::from(metadata.nuclide.mass_number());
    if head.number(0)? != f64::from(za) || head.number(1)? != metadata.atomic_weight_ratio {
        return Err(head.error("MF3 ZA/AWR disagrees with MF1 identity"));
    }
    for field in 2..6 {
        head.zero_integer(field)?;
    }
    let tab = reader.section(head.mat, 3, head.mt, checkpoint)?;
    tab.zero_integer(2)?;
    if tab.integer(3)? != 0 {
        return Err(EndfError::Unsupported("MF3 breakup flag LR must be zero"));
    }
    let curve = read_nonnegative_tab1(reader, tab, BARN_M2, checkpoint)?;
    reader.send(head.mat, 3, checkpoint)?;
    Ok(PointwiseCrossSection {
        data_id: Arc::from(""),
        metadata,
        reaction,
        energies_j: curve.energies_j.into(),
        cross_sections_m2: curve.values.into(),
        interpolation_regions: curve.interpolation_regions.into(),
        mass_difference_q_j: convert_si(tab.number(0)?, EV_J)?,
        reaction_q_j: convert_si(tab.number(1)?, EV_J)?,
    })
}

/// The same bounded TAB1 grammar serves cross sections and dimensionless
/// yields. Callers own the header's physical meaning and consume its SEND.
struct NonnegativeTabulation {
    energies_j: Vec<f64>,
    values: Vec<f64>,
    interpolation_regions: Vec<EndfInterpolationRegion>,
}

fn read_nonnegative_tab1(
    reader: &mut Reader<'_>,
    tab: Record<'_>,
    value_scale: f64,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<NonnegativeTabulation, EndfError> {
    let regions_count = tab.count(4)?;
    let points_count = tab.count(5)?;
    if regions_count == 0 || points_count < 2 || regions_count >= points_count {
        return Err(tab.error("TAB1 needs at least two points and valid interpolation regions"));
    }
    if points_count > reader.limits.max_points {
        return Err(budget("ENDF points"));
    }
    if regions_count > reader.limits.max_interpolation_regions {
        return Err(budget("ENDF interpolation regions"));
    }
    let needed = regions_count
        .div_ceil(3)
        .checked_add(points_count.div_ceil(3))
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| budget("ENDF TAB1 records"))?;
    if needed > reader.limits.max_records.saturating_sub(reader.count) {
        return Err(budget("ENDF TAB1 records"));
    }
    let mut regions = Vec::new();
    regions
        .try_reserve_exact(regions_count)
        .map_err(|_| budget("ENDF interpolation memory"))?;
    for _ in 0..regions_count.div_ceil(3) {
        let row = reader.section(tab.mat, tab.mf, tab.mt, checkpoint)?;
        for pair in 0..3 {
            if regions.len() == regions_count {
                row.zero_integer(pair * 2)?;
                row.zero_integer(pair * 2 + 1)?;
                continue;
            }
            let end_point = row.count(pair * 2)?;
            let previous = regions
                .last()
                .map_or(1, |r: &EndfInterpolationRegion| r.end_point);
            if end_point <= previous || end_point > points_count {
                return Err(row.error("ascending interpolation endpoints bounded by NP required"));
            }
            let interpolation = match row.integer(pair * 2 + 1)? {
                1 => EndfInterpolation::Histogram,
                2 => EndfInterpolation::Linear,
                _ => {
                    return Err(EndfError::Unsupported(
                        "only processed histogram and linear-linear TAB1 laws are supported",
                    ));
                }
            };
            regions.push(EndfInterpolationRegion {
                end_point,
                interpolation,
            });
        }
    }
    if regions.last().is_none_or(|r| r.end_point != points_count) {
        return Err(tab.error("last interpolation endpoint must equal NP"));
    }
    let mut energies = Vec::new();
    let mut values = Vec::new();
    energies
        .try_reserve_exact(points_count)
        .map_err(|_| budget("ENDF point memory"))?;
    values
        .try_reserve_exact(points_count)
        .map_err(|_| budget("ENDF point memory"))?;
    let mut previous_ev = None;
    for _ in 0..points_count.div_ceil(3) {
        let row = reader.section(tab.mat, tab.mf, tab.mt, checkpoint)?;
        for pair in 0..3 {
            let energy_ev = row.number(pair * 2)?;
            let source_value = row.number(pair * 2 + 1)?;
            if energies.len() == points_count {
                if energy_ev != 0.0 || source_value != 0.0 {
                    return Err(row.error("nonzero TAB1 padding"));
                }
                continue;
            }
            if !nonnegative(energy_ev) || !nonnegative(source_value) {
                return Err(row.error("nonnegative energy and tabulated value required"));
            }
            let energy_j = convert_si(energy_ev, EV_J)?;
            let value = convert_si(source_value, value_scale)?;
            if let Some(previous) = previous_ev {
                let last = energies[energies.len() - 1];
                if energy_ev < previous || (energy_ev > previous && energy_j <= last) {
                    return Err(row.error("ascending representable energies required"));
                }
                if energy_ev == previous
                    && energies.len() >= 2
                    && energies[energies.len() - 2] == energy_j
                {
                    return Err(
                        row.error("a discontinuity may contain only two equal energy points")
                    );
                }
            }
            previous_ev = Some(energy_ev);
            energies.push(energy_j);
            values.push(value);
        }
    }
    if energies[0] >= energies[energies.len() - 1] {
        return Err(tab.error("nonempty energy domain required"));
    }
    Ok(NonnegativeTabulation {
        energies_j: energies,
        values,
        interpolation_regions: regions,
    })
}

// ENDF allows embedded spaces, implicit signed exponents and blank numeric
// fields. Stack buffers bound temporary storage; D exponents are not ENDF-6.
fn integer(text: &str) -> Option<i32> {
    let mut buffer = [0_u8; 11];
    let mut length = 0;
    for byte in text.bytes().filter(|&b| b != b' ') {
        if length == buffer.len() {
            return None;
        }
        buffer[length] = byte;
        length += 1;
    }
    if length == 0 {
        return Some(0);
    }
    std::str::from_utf8(&buffer[..length]).ok()?.parse().ok()
}

fn number(text: &str) -> Option<f64> {
    let mut compact = [0_u8; 11];
    let mut length = 0;
    for byte in text.bytes().filter(|&b| b != b' ') {
        if length == compact.len()
            || !matches!(byte, b'0'..=b'9' | b'+' | b'-' | b'.' | b'e' | b'E')
        {
            return None;
        }
        compact[length] = byte;
        length += 1;
    }
    if length == 0 {
        return Some(0.0);
    }
    let explicit = compact[..length]
        .iter()
        .position(|&b| b == b'e' || b == b'E');
    let implicit = if explicit.is_none() {
        (1..length).find(|&i| compact[i] == b'+' || compact[i] == b'-')
    } else {
        None
    };
    let mut expanded = [0_u8; 12];
    let parsed = if let Some(position) = implicit {
        expanded[..position].copy_from_slice(&compact[..position]);
        expanded[position] = b'e';
        expanded[position + 1..length + 1].copy_from_slice(&compact[position..length]);
        std::str::from_utf8(&expanded[..length + 1])
            .ok()?
            .parse::<f64>()
            .ok()?
    } else {
        std::str::from_utf8(&compact[..length])
            .ok()?
            .parse::<f64>()
            .ok()?
    };
    let mantissa_end = explicit.or(implicit).unwrap_or(length);
    let nonzero_mantissa = compact[..mantissa_end]
        .iter()
        .any(|b| (b'1'..=b'9').contains(b));
    if !parsed.is_finite() || (parsed == 0.0 && nonzero_mantissa) {
        return None;
    }
    Some(parsed)
}

#[cfg(test)]
mod tests;
