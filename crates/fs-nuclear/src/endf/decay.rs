//! Evaluated spontaneous-decay data from ENDF-6 MF8/MT457.
//!
//! The original bounded ENDF reader selects an explicit MAT and nuclide/isomer.
//! Half-lives retain their supplied seconds and enter the shared `DecayLawV1`;
//! average emitted energies and decay-mode Q values are converted from eV to J.
//! Unknown average energies remain `None`. Source uncertainties are retained
//! declarations, not certified intervals or an assumed probability distribution.
//!
//! The first three average energies are electron-related, electromagnetic and
//! heavy-particle energy, including recoil. Their sum excludes neutrinos. The
//! optional seventeen-component form retains its neutrino component separately.
//! These are averages over ALL modes of the parent, not conditional branch
//! energies. Neither Q values nor emitted energies imply local deposition.
//!
//! Radiation spectra and their covariances are structurally traversed under the
//! same record budget, but are not interpreted, integrated or retained. Only
//! caller-supplied source identity and MF8 identity are claimed; this reader
//! neither requires nor certifies an MF1 library declaration. No isotope catalog,
//! branch normalization, daughter inference or radiation transport is supplied.
//! Format: ENDF-102 (2023), §§8.4–8.4.2,
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.

use std::sync::Arc;

use fs_qty::Time;
use fs_thermochem::decay::DecayLawV1;

use super::{EV_J, EndfError, EndfReadLimits, Reader, Record, budget, convert_si};
use crate::irradiation::NuclideId;
use crate::{Accumulator, FissionError, nonnegative, positive};

/// Both source identities must agree; no ground-state/isomer substitution occurs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfDecaySelection {
    pub material_number: u16,
    pub nuclide: NuclideId,
}

/// Bounds include all scanned and skipped records. Average-energy storage has
/// the format's hard maximum of seventeen value/uncertainty pairs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfDecayReadLimits {
    pub max_bytes: usize,
    pub max_records: usize,
    pub max_decay_modes: usize,
}

impl Default for EndfDecayReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 128 * 1024 * 1024,
            max_records: 1_600_000,
            max_decay_modes: 256,
        }
    }
}

/// Source-declared MF8 identity and original nuclear-state properties.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayMetadata {
    pub material_number: u16,
    pub nuclide: NuclideId,
    pub atomic_weight_ratio: f64,
    /// LIS, distinct from the nuclide's isomer number LISO.
    pub excited_state: u16,
    /// Number of structurally traversed, unretained radiation spectra.
    pub spectrum_count: usize,
    /// Source spin in units of hbar; `None` preserves ENDF's unknown sentinel.
    pub spin: Option<f64>,
    /// Source parity (+1 or -1); `None` denotes unknown.
    pub parity: Option<i8>,
}

/// One aggregate emitted-energy component in joules per parent decay.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayEnergy {
    /// ENDF's -1 eV sentinel is unknown, never an inferred zero contribution.
    pub value_j: Option<f64>,
    pub uncertainty_j: f64,
}

/// One source-declared decay mode. RTYP remains uninterpreted here, including
/// sequential emissions such as 1.5 and modes without a unique daughter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayMode {
    pub rtyp: f64,
    pub daughter_isomer: u16,
    /// Total available Q of this mode, not its deposited or emitted energy.
    pub q_j: f64,
    pub q_uncertainty_j: f64,
    pub probability: f64,
    pub probability_uncertainty: f64,
}

/// Immutable supplied evaluation. Parsing retains the source probabilities
/// without repairing their sum. Physical model construction admits its own
/// complete branch closure, topology and energy-destination assumptions.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayEvaluation {
    data_id: Arc<str>,
    metadata: DecayMetadata,
    law: DecayLawV1,
    half_life_uncertainty_s: f64,
    average_energies: Arc<[DecayEnergy]>,
    modes: Arc<[DecayMode]>,
    branch_fraction_sum: f64,
}

impl DecayEvaluation {
    /// Read caller-supplied text, accepting stable or radioactive MF8/MT457.
    /// Spectral/covariance payloads are skipped by their declared record counts,
    /// so truncation cannot silently become a complete decay evaluation.
    pub fn parse(
        data_id: String,
        selection: EndfDecaySelection,
        text: &str,
        limits: EndfDecayReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || selection.material_number == 0
            || selection.material_number > 9999
        {
            return Err(FissionError::InvalidInput(
                "decay source identity and MAT in 1..=9999 required",
            )
            .into());
        }
        if text.len() > limits.max_bytes {
            return Err(budget("ENDF decay bytes"));
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
        let mut candidate = None;
        while let Some(record) = reader.next(&mut checkpoint)? {
            if (record.mat, record.mf, record.mt) == (i32::from(selection.material_number), 8, 457)
            {
                if candidate.is_some() {
                    return Err(EndfError::AmbiguousSelection);
                }
                candidate = Some(read_evaluation(
                    &mut reader,
                    record,
                    selection,
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
    pub const fn metadata(&self) -> &DecayMetadata {
        &self.metadata
    }

    #[must_use]
    pub const fn law(&self) -> DecayLawV1 {
        self.law
    }

    #[must_use]
    pub const fn half_life_uncertainty_s(&self) -> f64 {
        self.half_life_uncertainty_s
    }

    /// Standard ENDF order: electron-related, electromagnetic, heavy-particle;
    /// then, when present, beta-minus, beta-plus, Auger, conversion electron,
    /// gamma, X-ray, internal bremsstrahlung, annihilation, alpha, recoil,
    /// spontaneous-fission, neutron, proton, neutrino. Detailed components are
    /// NOT additional heat on top of their first three aggregate totals.
    #[must_use]
    pub fn average_energies(&self) -> &[DecayEnergy] {
        &self.average_energies
    }

    #[must_use]
    pub fn modes(&self) -> &[DecayMode] {
        &self.modes
    }

    /// Compensated reduction of the supplied BR values; no normalization.
    #[must_use]
    pub const fn branch_fraction_sum(&self) -> f64 {
        self.branch_fraction_sum
    }
}

fn read_evaluation(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    selection: EndfDecaySelection,
    limits: EndfDecayReadLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<DecayEvaluation, EndfError> {
    let za = head.number(0)?;
    if za.fract() != 0.0 || !(1001.0..=999_999.0).contains(&za) {
        return Err(head.error("MF8 requires an isotope ZA"));
    }
    let za = za as u32;
    let excited_state =
        u16::try_from(head.integer(2)?).map_err(|_| head.error("invalid parent level LIS"))?;
    let isomer =
        u16::try_from(head.integer(3)?).map_err(|_| head.error("invalid parent isomer LISO"))?;
    let nuclide = NuclideId::new((za / 1000) as u16, (za % 1000) as u16, isomer)?;
    if nuclide != selection.nuclide || excited_state < isomer {
        return Err(head.error("selected MF8 nuclide/isomer or parent level disagrees"));
    }
    let atomic_weight_ratio = head.number(1)?;
    if !positive(atomic_weight_ratio) {
        return Err(head.error("positive nuclide-to-neutron mass ratio required"));
    }
    let stable = match head.integer(4)? {
        0 => false,
        1 => true,
        _ => return Err(head.error("NST must be zero or one")),
    };
    let spectrum_count = head.count(5)?;
    if stable && spectrum_count != 0 {
        return Err(head.error("stable nuclei have no decay spectra"));
    }
    require_records(reader, spectrum_count.saturating_mul(2))?;

    let half = reader.section(head.mat, 8, 457, checkpoint)?;
    for field in [2, 3, 5] {
        half.zero_integer(field)?;
    }
    let value_count = half.count(4)?;
    if (stable && value_count != 6) || (!stable && !matches!(value_count, 6 | 34)) {
        return Err(half.error("decay-energy LIST requires NC=3 or 17 (NC=3 when stable)"));
    }
    let half_life_s = half.number(0)?;
    let half_life_uncertainty_s = half.number(1)?;
    if !nonnegative(half_life_uncertainty_s)
        || (stable && (half_life_s != 0.0 || half_life_uncertainty_s != 0.0))
        || (!stable && !positive(half_life_s))
    {
        return Err(half.error("invalid half-life or half-life uncertainty"));
    }
    let law = if stable {
        DecayLawV1::stable()
    } else {
        DecayLawV1::from_half_life(Time::new(half_life_s))
            .map_err(|_| half.error("unrepresentable decay constant from half-life"))?
    };
    let mut raw_energies = [0.0; 34];
    read_values(
        reader,
        head.mat,
        &mut raw_energies[..value_count],
        checkpoint,
    )?;
    let mut energies = Vec::with_capacity(value_count / 2);
    for pair in raw_energies[..value_count].chunks_exact(2) {
        if !nonnegative(pair[1])
            || (pair[0] != -1.0 && !nonnegative(pair[0]))
            || (stable && (pair[0] != 0.0 || pair[1] != 0.0))
        {
            return Err(half.error("invalid average energy or its uncertainty"));
        }
        energies.push(DecayEnergy {
            value_j: if pair[0] == -1.0 {
                None
            } else {
                Some(convert_si(pair[0], EV_J)?)
            },
            uncertainty_j: convert_si(pair[1], EV_J)?,
        });
    }

    let mode_head = reader.section(head.mat, 8, 457, checkpoint)?;
    mode_head.zero_integer(2)?;
    mode_head.zero_integer(3)?;
    let spin = match mode_head.number(0)? {
        -77.777 => None,
        value if nonnegative(value) => Some(value),
        _ => return Err(mode_head.error("invalid spin or unknown-spin sentinel")),
    };
    let parity = match mode_head.number(1)? {
        -1.0 => Some(-1),
        0.0 => None,
        1.0 => Some(1),
        _ => return Err(mode_head.error("parity must be -1, 0 or +1")),
    };
    let mode_count = mode_head.count(5)?;
    if mode_count > limits.max_decay_modes {
        return Err(budget("ENDF decay modes"));
    }
    let mode_values = mode_head.count(4)?;
    let expected_values = if stable {
        6
    } else {
        mode_count
            .checked_mul(6)
            .ok_or_else(|| budget("ENDF decay-mode values"))?
    };
    if mode_values != expected_values || (stable && mode_count != 0) || (!stable && mode_count == 0)
    {
        return Err(mode_head.error("decay-mode LIST count disagrees with NST/NDK"));
    }
    require_records(reader, mode_values.div_ceil(6))?;
    let mut modes = Vec::with_capacity(mode_count);
    let mut probability_sum = Accumulator::new(0.0);
    if stable {
        let mut padding = [0.0; 6];
        read_values(reader, head.mat, &mut padding, checkpoint)?;
        if padding != [0.0; 6] {
            return Err(mode_head.error("stable decay-mode LIST must contain six zeros"));
        }
    } else {
        for _ in 0..mode_count {
            let record = reader.section(head.mat, 8, 457, checkpoint)?;
            let rtyp = record.number(0)?;
            let rfs = record.number(1)?;
            let q_j = convert_si(record.number(2)?, EV_J)?;
            let q_uncertainty = record.number(3)?;
            let probability = record.number(4)?;
            let probability_uncertainty = record.number(5)?;
            if !positive(rtyp)
                || rfs.fract() != 0.0
                || !(0.0..=f64::from(u16::MAX)).contains(&rfs)
                || !nonnegative(q_uncertainty)
                || !(0.0..=1.0).contains(&probability)
                || !nonnegative(probability_uncertainty)
            {
                return Err(
                    record.error("invalid decay-mode code, isomer, fraction or uncertainty")
                );
            }
            probability_sum.add(probability)?;
            modes.push(DecayMode {
                rtyp,
                daughter_isomer: rfs as u16,
                q_j,
                q_uncertainty_j: convert_si(q_uncertainty, EV_J)?,
                probability,
                probability_uncertainty,
            });
        }
    }

    skip_spectra(reader, head.mat, spectrum_count, checkpoint)?;
    reader.send(head.mat, 8, checkpoint)?;
    Ok(DecayEvaluation {
        data_id: Arc::from(""),
        metadata: DecayMetadata {
            material_number: selection.material_number,
            nuclide,
            atomic_weight_ratio,
            excited_state,
            spectrum_count,
            spin,
            parity,
        },
        law,
        half_life_uncertainty_s,
        average_energies: energies.into(),
        modes: modes.into(),
        branch_fraction_sum: probability_sum.value(),
    })
}

fn require_records(reader: &Reader<'_>, count: usize) -> Result<(), EndfError> {
    if count > reader.limits.max_records.saturating_sub(reader.count) {
        return Err(budget("ENDF decay payload records"));
    }
    Ok(())
}

fn read_values(
    reader: &mut Reader<'_>,
    mat: i32,
    values: &mut [f64],
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    require_records(reader, values.len().div_ceil(6))?;
    for chunk in values.chunks_mut(6) {
        let record = reader.section(mat, 8, 457, checkpoint)?;
        for field in 0..6 {
            let value = record.number(field)?;
            if field < chunk.len() {
                chunk[field] = value;
            } else if value != 0.0 {
                return Err(record.error("nonzero LIST padding"));
            }
        }
    }
    Ok(())
}

// Unretained data require no data-dependent allocation. Record counts, fixed
// numeric fields and section controls still detect malformed/truncated input.
fn skip_values(
    reader: &mut Reader<'_>,
    mat: i32,
    count: usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    require_records(reader, count.div_ceil(6))?;
    let mut remaining = count;
    while remaining != 0 {
        let record = reader.section(mat, 8, 457, checkpoint)?;
        for field in 0..6 {
            let value = record.number(field)?;
            if field >= remaining && value != 0.0 {
                return Err(record.error("nonzero skipped LIST/TAB1 padding"));
            }
        }
        remaining = remaining.saturating_sub(6);
    }
    Ok(())
}

fn skip_list(
    reader: &mut Reader<'_>,
    mat: i32,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    let head = reader.section(mat, 8, 457, checkpoint)?;
    head.number(0)?;
    head.number(1)?;
    head.integer(2)?;
    head.integer(3)?;
    head.count(5)?;
    skip_values(reader, mat, head.count(4)?, checkpoint)
}

fn skip_spectra(
    reader: &mut Reader<'_>,
    mat: i32,
    count: usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    let mut previous_type = -1.0;
    for _ in 0..count {
        let head = reader.section(mat, 8, 457, checkpoint)?;
        let styp = head.number(1)?;
        let lcon = head.integer(2)?;
        let lcov = head.integer(3)?;
        let discrete = head.count(5)?;
        if head.number(0)? != 0.0
            || styp.fract() != 0.0
            || styp <= previous_type
            || !(0..=2).contains(&lcon)
            || !(0..=3).contains(&lcov)
            || head.count(4)? != 6
            || (lcon == 1 && discrete != 0)
        {
            return Err(head.error("invalid radiation-spectrum header"));
        }
        previous_type = styp;
        skip_values(reader, mat, 6, checkpoint)?;
        if lcon != 1 {
            // Each discrete line has at least one LIST header. Remaining work
            // is checked again before consuming its actual value payload.
            require_records(reader, discrete)?;
            for _ in 0..discrete {
                skip_list(reader, mat, checkpoint)?;
            }
        }
        if lcon != 0 {
            skip_tab1(reader, mat, checkpoint)?;
            if lcov == 1 || lcov == 3 {
                skip_list(reader, mat, checkpoint)?;
            }
        }
        if lcov == 2 || lcov == 3 {
            skip_list(reader, mat, checkpoint)?;
        }
    }
    Ok(())
}

fn skip_tab1(
    reader: &mut Reader<'_>,
    mat: i32,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    let head = reader.section(mat, 8, 457, checkpoint)?;
    head.number(0)?;
    if head.number(1)? != 0.0 {
        return Err(head.error("nonzero reserved spectral TAB1 field"));
    }
    head.zero_integer(2)?;
    head.zero_integer(3)?;
    let regions = head.count(4)?;
    let points = head.count(5)?;
    if regions == 0 || points < 2 || regions >= points {
        return Err(head.error("spectral TAB1 needs points and interpolation ranges"));
    }
    let point_values = points
        .checked_mul(2)
        .ok_or_else(|| budget("ENDF spectral TAB1 values"))?;
    require_records(
        reader,
        regions
            .div_ceil(3)
            .checked_add(point_values.div_ceil(6))
            .ok_or_else(|| budget("ENDF spectral TAB1 records"))?,
    )?;
    let mut remaining = regions;
    let mut previous = 1;
    while remaining != 0 {
        let record = reader.section(mat, 8, 457, checkpoint)?;
        for pair in 0..3 {
            let endpoint = record.count(2 * pair)?;
            let law = record.integer(2 * pair + 1)?;
            if pair < remaining {
                if endpoint <= previous || endpoint > points || !(1..=6).contains(&law) {
                    return Err(record.error("invalid spectral interpolation range"));
                }
                previous = endpoint;
            } else if endpoint != 0 || law != 0 {
                return Err(record.error("nonzero spectral interpolation padding"));
            }
        }
        remaining = remaining.saturating_sub(3);
    }
    if previous != points {
        return Err(head.error("spectral interpolation endpoint disagrees with NP"));
    }
    skip_values(reader, mat, point_values, checkpoint)
}

#[cfg(test)]
mod tests;
