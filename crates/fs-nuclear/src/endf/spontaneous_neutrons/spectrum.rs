//! Energy-independent MF5 LF=1 partial distributions. The common continuum
//! point/interpolation records carry SI quantities without fabricating an MF8
//! decay spectrum or a normalization uncertainty that the source did not give.

use super::*;
use crate::endf::decay::spectra::{DecayInterpolation, DecaySpectrumPoint};
use crate::endf::{EV_J, convert_si};

/// One original partial distribution (or delayed family), in source order.
/// Original incident-grid fields are retained as metadata: every supplied row
/// must declare the same outgoing law and constant component probability.
#[derive(Clone, Debug, PartialEq)]
pub struct SpontaneousFissionSpectrumComponent {
    pub probability: f64,
    pub probability_energies_j: Arc<[f64]>,
    pub probability_interpolation: Arc<[DecayInterpolation]>,
    pub incident_energies_j: Arc<[f64]>,
    pub incident_interpolation: Arc<[DecayInterpolation]>,
    pub points: Arc<[DecaySpectrumPoint]>,
    pub interpolation: Arc<[DecayInterpolation]>,
}

/// Explicit prompt MT18 or delayed MT455 source. Neither this source type nor
/// its parser supplies family dynamics, a decay-branch factor, or renormalized
/// probabilities. A consumer must admit normalization using its own allowance.
#[derive(Clone, Debug, PartialEq)]
pub struct SpontaneousFissionNeutronSpectrum {
    data_id: Arc<str>,
    metadata: EndfMetadata,
    kind: FissionNeutronKind,
    components: Arc<[SpontaneousFissionSpectrumComponent]>,
}

impl SpontaneousFissionNeutronSpectrum {
    pub fn parse(
        data_id: String,
        selection: EndfSpontaneousNeutronSelection,
        kind: FissionNeutronKind,
        text: &str,
        limits: SpontaneousNeutronReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        let mt = match kind {
            FissionNeutronKind::Prompt => 18,
            FissionNeutronKind::Delayed => 455,
            FissionNeutronKind::Total => {
                return Err(EndfError::Unsupported(
                    "a total neutron yield does not identify a prompt or delayed MF5 spectrum",
                ));
            }
        };
        let (metadata, components) = selected(
            &data_id,
            selection,
            text,
            limits,
            5,
            mt,
            &mut checkpoint,
            |reader, head, metadata, checkpoint| {
                match_identity(head, &metadata)?;
                head.zero_integer(2)?;
                head.zero_integer(3)?;
                head.zero_integer(5)?;
                let count = head.count(4)?;
                if count == 0 {
                    return Err(head.error("nonempty MF5 component declaration required"));
                }
                if count > limits.max_components {
                    return Err(budget("spontaneous-neutron spectrum components"));
                }
                require_records(reader, count.saturating_mul(3).saturating_add(1))?;
                let mut components = reserved(count, "spontaneous-neutron component memory")?;
                let mut work = SourceCounts::default();
                for _ in 0..count {
                    checkpoint()?;
                    components.push(component(reader, head, limits, &mut work, checkpoint)?);
                }
                reader.send(head.mat, 5, checkpoint)?;
                Ok((metadata, components))
            },
        )?;
        Ok(Self {
            data_id: data_id.into(),
            metadata,
            kind,
            components: components.into(),
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
    pub fn components(&self) -> &[SpontaneousFissionSpectrumComponent] {
        &self.components
    }
}

#[derive(Default)]
struct SourceCounts {
    incident: usize,
    outgoing: usize,
    regions: usize,
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

fn component(
    reader: &mut Reader<'_>,
    section: Record<'_>,
    limits: SpontaneousNeutronReadLimits,
    counts: &mut SourceCounts,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<SpontaneousFissionSpectrumComponent, EndfError> {
    let probability = reader.section(section.mat, 5, section.mt, checkpoint)?;
    if probability.number(0)? != 0.0 || probability.number(1)? != 0.0 {
        return Err(probability.error("nonzero reserved LF1 probability field"));
    }
    probability.zero_integer(2)?;
    if probability.integer(3)? != 1 {
        return Err(EndfError::Unsupported(
            "spontaneous-neutron spectra currently require tabulated MF5 LF=1",
        ));
    }
    let probability_points = probability.count(5)?;
    if probability_points == 0 || probability_points > limits.records.max_points {
        return Err(budget("spontaneous-neutron probability points"));
    }
    consume(
        &mut counts.incident,
        probability_points,
        limits.max_incident_points,
        "spontaneous-neutron incident points",
    )?;
    let probability_interpolation = interpolation(
        reader,
        probability,
        probability_points,
        limits,
        counts,
        checkpoint,
    )?;
    require_records(reader, probability_points.div_ceil(3))?;
    let mut probability_energies =
        reserved(probability_points, "spontaneous-neutron probability memory")?;
    let mut constant = None;
    while probability_energies.len() < probability_points {
        let row = reader.section(section.mat, 5, section.mt, checkpoint)?;
        for slot in 0..3 {
            let energy = row.number(2 * slot)?;
            let value = row.number(2 * slot + 1)?;
            if probability_energies.len() == probability_points {
                if energy != 0.0 || value != 0.0 {
                    return Err(row.error("nonzero MF5 probability padding"));
                }
                continue;
            }
            if !nonnegative(energy) || !(0.0..=1.0).contains(&value) {
                return Err(row.error(
                    "finite nonnegative incident energies and probabilities in [0,1] required",
                ));
            }
            let energy = convert_si(energy, EV_J)?;
            if probability_energies
                .last()
                .is_some_and(|&previous| previous >= energy)
            {
                return Err(row.error("strictly ascending incident probability energies required"));
            }
            if constant.is_some_and(|original| original != value) {
                return Err(EndfError::Unsupported(
                    "spontaneous-neutron component probability must be energy-independent",
                ));
            }
            constant = Some(value);
            probability_energies.push(energy);
        }
    }
    log_domain(
        &probability_energies,
        &probability_interpolation,
        checkpoint,
    )?;

    let tab2 = reader.section(section.mat, 5, section.mt, checkpoint)?;
    zero_prefix(tab2)?;
    let rows = tab2.count(5)?;
    if rows == 0 || rows > limits.records.max_points {
        return Err(budget("spontaneous-neutron incident spectrum rows"));
    }
    consume(
        &mut counts.incident,
        rows,
        limits.max_incident_points,
        "spontaneous-neutron incident points",
    )?;
    let incident_interpolation = interpolation(reader, tab2, rows, limits, counts, checkpoint)?;
    require_records(reader, rows.saturating_mul(2))?;
    let mut incident_energies = reserved(rows, "spontaneous-neutron incident memory")?;
    let mut common: Option<(Vec<DecaySpectrumPoint>, Vec<DecayInterpolation>)> = None;
    for _ in 0..rows {
        let tab = reader.section(section.mat, 5, section.mt, checkpoint)?;
        if tab.number(0)? != 0.0 || !nonnegative(tab.number(1)?) {
            return Err(tab.error("finite nonnegative MF5 incident energy required"));
        }
        tab.zero_integer(2)?;
        tab.zero_integer(3)?;
        let incident = convert_si(tab.number(1)?, EV_J)?;
        if incident_energies
            .last()
            .is_some_and(|&previous| previous >= incident)
        {
            return Err(tab.error("strictly ascending incident spectrum energies required"));
        }
        incident_energies.push(incident);
        let points = tab.count(5)?;
        if points < 2 {
            return Err(tab.error("at least two outgoing spectrum points required"));
        }
        if points > limits.records.max_points {
            return Err(budget("spontaneous-neutron spectrum points"));
        }
        consume(
            &mut counts.outgoing,
            points,
            limits.max_spectrum_points,
            "spontaneous-neutron spectrum points",
        )?;
        let laws = interpolation(reader, tab, points, limits, counts, checkpoint)?;
        require_records(reader, points.div_ceil(3))?;
        let mut outgoing: Vec<DecaySpectrumPoint> =
            reserved(points, "spontaneous-neutron spectrum memory")?;
        while outgoing.len() < points {
            let row = reader.section(section.mat, 5, section.mt, checkpoint)?;
            for slot in 0..3 {
                let energy = row.number(2 * slot)?;
                let density = row.number(2 * slot + 1)?;
                if outgoing.len() == points {
                    if energy != 0.0 || density != 0.0 {
                        return Err(row.error("nonzero MF5 spectrum padding"));
                    }
                    continue;
                }
                if !nonnegative(energy) || !nonnegative(density) {
                    return Err(
                        row.error("finite nonnegative outgoing energy and density required")
                    );
                }
                let energy_j = convert_si(energy, EV_J)?;
                let density_per_j = convert_si(density, EV_J.recip())?;
                let n = outgoing.len();
                if (n > 0 && energy_j < outgoing[n - 1].energy_j)
                    || (n > 1
                        && energy_j == outgoing[n - 1].energy_j
                        && energy_j == outgoing[n - 2].energy_j)
                {
                    return Err(row.error(
                        "outgoing energies must ascend with at most paired discontinuities",
                    ));
                }
                outgoing.push(DecaySpectrumPoint {
                    energy_j,
                    density_per_j,
                });
            }
        }
        if outgoing[0].energy_j >= outgoing[points - 1].energy_j
            || outgoing[0].density_per_j != 0.0
            || outgoing[points - 1].density_per_j != 0.0
        {
            return Err(
                tab.error("LF1 requires a nonempty outgoing domain with zero density endpoints")
            );
        }
        let mut region = 0;
        for (index, pair) in outgoing.windows(2).enumerate() {
            checkpoint()?;
            while laws[region].end_point < index + 2 {
                region += 1;
            }
            if laws[region].law == 3
                && pair[0].energy_j < pair[1].energy_j
                && pair[0].energy_j == 0.0
            {
                return Err(
                    tab.error("log-energy interpolation requires positive segment energies")
                );
            }
        }
        if let Some((original_points, original_laws)) = &common {
            if original_points != &outgoing || original_laws != &laws {
                return Err(EndfError::Unsupported(
                    "spontaneous-neutron outgoing spectrum must be energy-independent",
                ));
            }
        } else {
            common = Some((outgoing, laws));
        }
    }
    log_domain(&incident_energies, &incident_interpolation, checkpoint)?;
    let (points, interpolation) = common.ok_or_else(|| tab2.error("missing outgoing spectrum"))?;
    Ok(SpontaneousFissionSpectrumComponent {
        probability: constant.ok_or_else(|| probability.error("missing component probability"))?,
        probability_energies_j: probability_energies.into(),
        probability_interpolation: probability_interpolation.into(),
        incident_energies_j: incident_energies.into(),
        incident_interpolation: incident_interpolation.into(),
        points: points.into(),
        interpolation: interpolation.into(),
    })
}

fn interpolation(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    points: usize,
    limits: SpontaneousNeutronReadLimits,
    counts: &mut SourceCounts,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Vec<DecayInterpolation>, EndfError> {
    let regions = head.count(4)?;
    if (points > 1 && regions == 0) || regions > points {
        return Err(head.error("explicit interpolation ranges bounded by point count required"));
    }
    if regions > limits.records.max_interpolation_regions {
        return Err(budget("spontaneous-neutron interpolation ranges"));
    }
    consume(
        &mut counts.regions,
        regions,
        limits.max_interpolation_regions,
        "spontaneous-neutron interpolation ranges",
    )?;
    require_records(reader, regions.div_ceil(3))?;
    let mut result = reserved(regions, "spontaneous-neutron interpolation memory")?;
    let mut previous = usize::from(points > 1);
    while result.len() < regions {
        let row = reader.section(head.mat, head.mf, head.mt, checkpoint)?;
        for slot in 0..3 {
            let endpoint = row.count(2 * slot)?;
            let law = row.integer(2 * slot + 1)?;
            if result.len() == regions {
                if endpoint != 0 || law != 0 {
                    return Err(row.error("nonzero MF5 interpolation padding"));
                }
                continue;
            }
            if endpoint <= previous || endpoint > points {
                return Err(
                    row.error("ascending interpolation endpoints bounded by point count required")
                );
            }
            if !(1..=3).contains(&law) {
                return Err(EndfError::Unsupported(
                    "MF5 LF1 interpolation requires ENDF INT=1,2,3",
                ));
            }
            result.push(DecayInterpolation {
                end_point: endpoint,
                law: law as u8,
            });
            previous = endpoint;
        }
    }
    if regions > 0 && previous != points {
        return Err(head.error("last interpolation endpoint must equal point count"));
    }
    Ok(result)
}

fn log_domain(
    energies: &[f64],
    interpolation: &[DecayInterpolation],
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    let mut region = 0;
    for (index, pair) in energies.windows(2).enumerate() {
        checkpoint()?;
        while interpolation[region].end_point < index + 2 {
            region += 1;
        }
        if interpolation[region].law == 3 && pair[0] == 0.0 {
            return Err(EndfError::Unsupported(
                "log-energy interpolation requires positive incident energies",
            ));
        }
    }
    Ok(())
}
