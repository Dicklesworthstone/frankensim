//! Parent-wide emitted radiation from ENDF-6 MF8/MT457, §8.4, Eq. (8.1).
//!
//! The absolute spectral intensity is `FC * RP(E) + FD * sum(RI * delta(E-ER))`.
//! These normalizations include decay branching already. The records are not
//! conditional on a chosen daughter and do not add heat to the independently
//! supplied aggregate energy categories. In particular, beta records contain
//! endpoint energies and transition information, not monoenergetic electrons.
//! Their shape needs a supplied continuum or a separately admitted beta model.
//!
//! Energies are J, continuum densities are J^-1, and normalization/intensity
//! factors are dimensionless. Values and declared uncertainties remain separate;
//! no distribution, covariance propagation or normalization repair is inferred.

use std::sync::Arc;

use super::{
    DecayEnergy, EV_J, EndfDecayReadLimits, EndfError, FissionError, Reader, Record, budget,
    convert_si, nonnegative, read_values, require_records, skip_list,
};

/// ENDF radiation identity (STYP), distinct from the originating decay mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DecayRadiation {
    Gamma,
    BetaMinus,
    ElectronCapturePositron,
    Alpha,
    Neutron,
    FissionFragment,
    Proton,
    Electron,
    XRay,
    ElectronAntineutrino,
    ElectronNeutrino,
    /// Retained without assigning a particle species or emission semantics.
    Other(u16),
}

impl DecayRadiation {
    fn from_code(code: u16) -> Self {
        match code {
            0 => Self::Gamma,
            1 => Self::BetaMinus,
            2 => Self::ElectronCapturePositron,
            4 => Self::Alpha,
            5 => Self::Neutron,
            6 => Self::FissionFragment,
            7 => Self::Proton,
            8 => Self::Electron,
            9 => Self::XRay,
            10 => Self::ElectronAntineutrino,
            11 => Self::ElectronNeutrino,
            other => Self::Other(other),
        }
    }
}

/// A supplied nonnegative value and its supplied nonnegative uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectrumValue {
    pub value: f64,
    pub uncertainty: f64,
}

/// One original discrete LIST. `energy_j` is an endpoint for beta records and
/// an electron-capture energy for STYP=2, not an emitted monoenergetic particle.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscreteDecayEmission {
    pub energy_j: SpectrumValue,
    /// RTYP of the originating decay, including sequential codes or unknown 10.
    pub decay_mode: f64,
    /// Original TYPE: beta allowed/forbidden-transition classification.
    pub transition_type: f64,
    /// RI in the source's relative units, before the spectrum's FD multiplier.
    pub intensity: SpectrumValue,
    /// RIS: pair coefficient for gamma, positron intensity for STYP=2.
    /// It must not be added to RI as an independent line.
    pub pair_or_positron: Option<SpectrumValue>,
    pub conversion_total: Option<SpectrumValue>,
    pub conversion_k: Option<SpectrumValue>,
    pub conversion_l: Option<SpectrumValue>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecaySpectrumPoint {
    pub energy_j: f64,
    /// RP(E), transformed from probability/eV to probability/J.
    pub density_per_j: f64,
}

/// Original ENDF interpolation range. Retention does not imply that every
/// consumer implements every law; unsupported numerical interpretations refuse.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DecayInterpolation {
    /// One-based NBT endpoint; the final endpoint equals the point count.
    pub end_point: usize,
    /// ENDF INT in 1..=6, without converting logarithmic laws to linear ones.
    pub law: u8,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DecayContinuum {
    pub decay_mode: f64,
    pub points: Arc<[DecaySpectrumPoint]>,
    pub interpolation: Arc<[DecayInterpolation]>,
}

/// Retained original spectrum. Clone shares the immutable point/line arrays.
#[derive(Clone, Debug, PartialEq)]
pub struct DecaySpectrum {
    pub radiation: DecayRadiation,
    pub discrete_normalization: SpectrumValue,
    /// ERAV, which can include recoil for alpha and is not inferred from lines.
    /// Unknown (-1 eV) stays unknown instead of becoming zero.
    pub average_energy_j: DecayEnergy,
    pub continuum_normalization: SpectrumValue,
    pub discrete: Arc<[DiscreteDecayEmission]>,
    pub continuum: Option<DecayContinuum>,
    /// LCOV is retained; covariance payloads are checked but not interpreted.
    pub covariance_flag: u8,
}

fn pair(record: &Record<'_>, values: [f64; 2], scale: f64) -> Result<SpectrumValue, EndfError> {
    if !values.into_iter().all(nonnegative) {
        return Err(record.error("negative spectral value or uncertainty"));
    }
    Ok(SpectrumValue {
        value: convert_si(values[0], scale)?,
        uncertainty: convert_si(values[1], scale)?,
    })
}

fn consume(
    count: &mut usize,
    added: usize,
    maximum: usize,
    name: &'static str,
) -> Result<(), EndfError> {
    *count = count.checked_add(added).ok_or_else(|| budget(name))?;
    if *count > maximum {
        return Err(budget(name));
    }
    Ok(())
}

pub(super) fn read_spectra(
    reader: &mut Reader<'_>,
    mat: i32,
    count: usize,
    limits: EndfDecayReadLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Vec<DecaySpectrum>, EndfError> {
    let mut spectra = Vec::with_capacity(count);
    let mut previous_type = -1.0;
    let mut line_count = 0;
    let mut point_count = 0;
    let mut region_count = 0;
    for _ in 0..count {
        let head = reader.section(mat, 8, 457, checkpoint)?;
        let styp = head.number(1)?;
        let lcon = head.integer(2)?;
        let lcov = head.integer(3)?;
        let discrete_count = head.count(5)?;
        if head.number(0)? != 0.0
            || styp.fract() != 0.0
            || !(0.0..=f64::from(u16::MAX)).contains(&styp)
            || styp <= previous_type
            || !(0..=2).contains(&lcon)
            || !(0..=3).contains(&lcov)
            || head.count(4)? != 6
            || (lcon == 1 && discrete_count != 0)
        {
            return Err(head.error("invalid radiation-spectrum header"));
        }
        previous_type = styp;
        consume(
            &mut line_count,
            discrete_count,
            limits.max_discrete_lines,
            "ENDF discrete decay lines",
        )?;
        // Each line needs a header and at least one data record. Admission must
        // precede capacity allocation, even for an enormous declared NER.
        require_records(reader, discrete_count.saturating_mul(2).saturating_add(1))?;
        let mut normalization = [0.0; 6];
        read_values(reader, mat, &mut normalization, checkpoint)?;
        let discrete_normalization = pair(&head, [normalization[0], normalization[1]], 1.0)?;
        let continuum_normalization = pair(&head, [normalization[4], normalization[5]], 1.0)?;
        let average_energy_j = if normalization[2] == -1.0 {
            DecayEnergy {
                value_j: None,
                uncertainty_j: pair(&head, [0.0, normalization[3]], EV_J)?.uncertainty,
            }
        } else {
            let energy = pair(&head, [normalization[2], normalization[3]], EV_J)?;
            DecayEnergy {
                value_j: Some(energy.value),
                uncertainty_j: energy.uncertainty,
            }
        };
        let mut discrete = Vec::with_capacity(discrete_count);
        for _ in 0..discrete_count {
            let line = reader.section(mat, 8, 457, checkpoint)?;
            for field in [2, 3, 5] {
                line.zero_integer(field)?;
            }
            let nt = line.count(4)?;
            if !matches!(nt, 6 | 8 | 12) {
                return Err(line.error("discrete radiation LIST requires NT=6, 8 or 12"));
            }
            let energy_j = pair(&line, [line.number(0)?, line.number(1)?], EV_J)?;
            let mut raw = [0.0; 12];
            read_values(reader, mat, &mut raw[..nt], checkpoint)?;
            if !nonnegative(raw[0]) || !nonnegative(raw[1]) || raw[1].fract() != 0.0 {
                return Err(line.error("invalid radiation origin or transition type"));
            }
            discrete.push(DiscreteDecayEmission {
                energy_j,
                decay_mode: raw[0],
                transition_type: raw[1],
                intensity: pair(&line, [raw[2], raw[3]], 1.0)?,
                pair_or_positron: Some(pair(&line, [raw[4], raw[5]], 1.0)?),
                conversion_total: if nt >= 8 {
                    Some(pair(&line, [raw[6], raw[7]], 1.0)?)
                } else {
                    None
                },
                conversion_k: if nt == 12 {
                    Some(pair(&line, [raw[8], raw[9]], 1.0)?)
                } else {
                    None
                },
                conversion_l: if nt == 12 {
                    Some(pair(&line, [raw[10], raw[11]], 1.0)?)
                } else {
                    None
                },
            });
        }
        let continuum = if lcon != 0 {
            Some(read_continuum(
                reader,
                mat,
                limits,
                &mut point_count,
                &mut region_count,
                checkpoint,
            )?)
        } else {
            None
        };
        if lcon != 0 && (lcov == 1 || lcov == 3) {
            skip_list(reader, mat, checkpoint)?;
        }
        if lcov == 2 || lcov == 3 {
            skip_list(reader, mat, checkpoint)?;
        }
        spectra.push(DecaySpectrum {
            radiation: DecayRadiation::from_code(styp as u16),
            discrete_normalization,
            average_energy_j,
            continuum_normalization,
            discrete: discrete.into(),
            continuum,
            covariance_flag: lcov as u8,
        });
    }
    Ok(spectra)
}

fn read_continuum(
    reader: &mut Reader<'_>,
    mat: i32,
    limits: EndfDecayReadLimits,
    point_count: &mut usize,
    region_count: &mut usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<DecayContinuum, EndfError> {
    let head = reader.section(mat, 8, 457, checkpoint)?;
    let decay_mode = head.number(0)?;
    if !nonnegative(decay_mode) || head.number(1)? != 0.0 {
        return Err(head.error("invalid spectral TAB1 origin or reserved field"));
    }
    head.zero_integer(2)?;
    head.zero_integer(3)?;
    let regions = head.count(4)?;
    let points = head.count(5)?;
    if regions == 0 || points < 2 || regions >= points {
        return Err(head.error("spectral TAB1 needs points and interpolation ranges"));
    }
    consume(
        point_count,
        points,
        limits.max_continuum_points,
        "ENDF decay continuum points",
    )?;
    consume(
        region_count,
        regions,
        limits.max_interpolation_regions,
        "ENDF decay interpolation regions",
    )?;
    require_records(
        reader,
        regions.div_ceil(3).saturating_add(points.div_ceil(3)),
    )?;
    let mut interpolation = Vec::with_capacity(regions);
    let mut previous = 1;
    while interpolation.len() < regions {
        let record = reader.section(mat, 8, 457, checkpoint)?;
        for slot in 0..3 {
            let end_point = record.count(2 * slot)?;
            let law = record.integer(2 * slot + 1)?;
            if interpolation.len() < regions {
                if end_point <= previous || end_point > points || !(1..=6).contains(&law) {
                    return Err(record.error("invalid spectral interpolation range"));
                }
                previous = end_point;
                interpolation.push(DecayInterpolation {
                    end_point,
                    law: law as u8,
                });
            } else if end_point != 0 || law != 0 {
                return Err(record.error("nonzero spectral interpolation padding"));
            }
        }
    }
    if previous != points {
        return Err(head.error("spectral interpolation endpoint disagrees with NP"));
    }
    let mut data: Vec<DecaySpectrumPoint> = Vec::with_capacity(points);
    while data.len() < points {
        let record = reader.section(mat, 8, 457, checkpoint)?;
        for slot in 0..3 {
            let energy = record.number(2 * slot)?;
            let density = record.number(2 * slot + 1)?;
            if data.len() < points {
                if !nonnegative(energy) || !nonnegative(density) {
                    return Err(record
                        .error("finite nonnegative continuum energies and densities required"));
                }
                let energy_j = convert_si(energy, EV_J)?;
                let density_per_j = convert_si(density, EV_J.recip())?;
                let n = data.len();
                if (n > 0 && energy_j < data[n - 1].energy_j)
                    || (n > 1
                        && energy_j == data[n - 1].energy_j
                        && energy_j == data[n - 2].energy_j)
                {
                    return Err(record.error(
                        "spectral energies must ascend with at most paired discontinuities",
                    ));
                }
                data.push(DecaySpectrumPoint {
                    energy_j,
                    density_per_j,
                });
            } else if energy != 0.0 || density != 0.0 {
                return Err(record.error("nonzero spectral TAB1 padding"));
            }
        }
    }
    if data[0].energy_j >= data[points - 1].energy_j {
        return Err(head.error("nonempty continuous spectrum energy domain required"));
    }
    Ok(DecayContinuum {
        decay_mode,
        points: data.into(),
        interpolation: interpolation.into(),
    })
}
