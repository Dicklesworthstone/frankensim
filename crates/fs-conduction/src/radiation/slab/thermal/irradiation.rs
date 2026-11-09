//! Prescribed directional power on the original radiative thermal cells.

use super::*;
use fs_thermal::enthalpy::conduction::implicit::source::HeatSourceSchedule;

/// Maximum independently directed sources on one gray or spectral slab.
pub const MAX_SLAB_THERMAL_BEAMS: usize = 64;

/// One mutually incoherent, fixed-optics beam and its original power history.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabThermalBeam {
    /// Unique physical source identifier; fixes the compensated reduction order.
    pub id: u64,
    /// Original finite-band index, or zero for a full-spectrum gray driver.
    /// The beam must share that band's exact layers and exterior index sources.
    pub band: usize,
    pub response: super::super::SlabBeamResponse,
    /// Power crossing the original [front, back] face per slab area [W/m²],
    /// before Fresnel entry, at schedule multiplier one. This is ALREADY
    /// projected irradiance, not beam-normal irradiance or radiance. Neither
    /// another incidence cosine nor an exterior n² factor is applied.
    /// The two directions use the response's conjugate Snell cosines.
    pub irradiation_w_m2: [f64; 2],
    pub schedule: HeatSourceSchedule,
}

/// Exact prescribed beam supply over one native interval, for frozen optics.
/// These diagnostics are not an additional thermal-energy credit. The native
/// external-heat account already contains `mean_sources_w * duration_s` once.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabBeamInterval {
    pub duration_s: f64,
    pub mean_sources_w: BTreeMap<CellId, f64>,
    pub incident_j: [f64; 2],
    /// Energy leaving the named face after reflection.
    pub reflected_j: [f64; 2],
    /// Energy leaving the named face after transmission from the opposite face.
    pub transmitted_j: [f64; 2],
    pub absorbed_j: f64,
    /// Signed incident minus reflected minus transmitted minus actual cell heat.
    pub balance_residual_j: f64,
}

#[derive(Clone, Debug)]
struct PreparedBeam {
    source: SlabThermalBeam,
    cell_power_w: Arc<[f64]>,
    // Incident, reflected and transmitted powers, each in front/back order.
    powers_w: [f64; 6],
}

/// Prepared original beam fields integrated over the existing thermal cells.
/// Preparation is O(beams*cells); no new optical interface, angular quadrature
/// or thermal solver is introduced. Storage and each hot reduction are bounded
/// by 64 beams times 64 cells. Original schedule storage/work is additionally
/// bounded by the caller's total-knot allowance at attachment.
#[derive(Clone, Debug)]
pub struct SlabBeamIllumination {
    beams: Arc<[PreparedBeam]>,
    cells: Arc<[CellId]>,
}

impl SlabBeamIllumination {
    pub(super) fn prepare(
        mut beams: Vec<SlabThermalBeam>,
        maximum_schedule_knots: usize,
        drivers: &[GraySlabThermal<'_>],
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if beams.is_empty() || beams.len() > MAX_SLAB_THERMAL_BEAMS {
            return Err(DepositionError::Budget("slab directional source count"));
        }
        let mut knots = 0usize;
        for beam in &beams {
            checkpoint()?;
            knots = knots
                .checked_add(beam.schedule.knots().len())
                .ok_or(DepositionError::Budget("slab beam schedule knots"))?;
            if knots > maximum_schedule_knots {
                return Err(DepositionError::Budget("slab beam schedule knots"));
            }
        }
        beams.sort_unstable_by_key(|beam| beam.id);
        if beams.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(DepositionError::InvalidInput(
                "slab beam IDs must be unique",
            ));
        }
        let mut prepared = Vec::with_capacity(beams.len());
        for source in beams {
            checkpoint()?;
            let driver = drivers
                .get(source.band)
                .ok_or(DepositionError::InvalidInput(
                    "slab beam must select an original optical band",
                ))?;
            let slab = &driver.slab;
            if source.response.layers() != slab.layers()
                || source.response.exterior_refractive_indices()
                    != slab.exterior_refractive_indices()
                || source.response.exterior_index_sources() != slab.exterior_index_sources()
            {
                return Err(DepositionError::InvalidInput(
                    "slab beam must retain original band layers and both exterior index sources",
                ));
            }
            source.schedule.multiplier_at(time_s)?;
            let mut whole = [0.0; MAX_SLAB_LAYERS];
            let exchange = source
                .response
                .exchange(source.irradiation_w_m2, &mut whole[..slab.layer_count()])
                .map_err(|_| {
                    DepositionError::InvalidInput("invalid original slab beam irradiation")
                })?;
            let mut cell_power_w = Vec::with_capacity(driver.layers.len());
            let mut absorbed = 0.0;
            let mut correction = 0.0;
            for (i, layer) in driver.layers.iter().enumerate() {
                checkpoint()?;
                let fractions = if let Some(response) = &driver.cell_response {
                    source
                        .response
                        .layer_interval_absorptance(
                            layer.optical_layer,
                            response.cells()[i].interval_m,
                        )
                        .map_err(|_| DepositionError::Arithmetic)?
                } else {
                    source.response.layer_absorptance()[layer.optical_layer]
                };
                let flux = product(fractions[0], source.irradiation_w_m2[0])?
                    + product(fractions[1], source.irradiation_w_m2[1])?;
                add(&mut absorbed, &mut correction, flux);
                cell_power_w.push(product(driver.area_m2, flux)?);
            }
            let scale = source.irradiation_w_m2[0].max(source.irradiation_w_m2[1]);
            if !absorbed.is_finite()
                || (absorbed - exchange.absorbed_w_m2).abs()
                    > 1024.0 * f64::EPSILON * driver.layers.len() as f64 * scale
            {
                return Err(DepositionError::Arithmetic);
            }
            let mut powers_w = [0.0; 6];
            for (out, flux) in powers_w.iter_mut().zip(
                source
                    .irradiation_w_m2
                    .into_iter()
                    .chain(exchange.reflected_w_m2)
                    .chain(exchange.transmitted_w_m2),
            ) {
                *out = product(driver.area_m2, flux)?;
            }
            prepared.push(PreparedBeam {
                source,
                cell_power_w: cell_power_w.into(),
                powers_w,
            });
        }
        let result = Self {
            beams: prepared.into(),
            cells: drivers[0].layers.iter().map(|layer| layer.id).collect(),
        };
        checkpoint()?;
        Ok(result)
    }

    /// Original sources in deterministic ID order, including their schedules.
    pub fn beams(&self) -> impl ExactSizeIterator<Item = &SlabThermalBeam> {
        self.beams.iter().map(|beam| &beam.source)
    }

    pub(super) fn admit_time(&self, time_s: f64) -> Result<(), DepositionError> {
        for beam in self.beams.iter() {
            // Dark or transparent beams still retain their source time domain.
            beam.source.schedule.multiplier_at(time_s)?;
        }
        Ok(())
    }

    pub(super) fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.beams
            .iter()
            .filter_map(|beam| beam.source.schedule.next_transition_time_s(time_s))
            .min_by(f64::total_cmp)
    }

    /// Analytically integrate the ORIGINAL continuous piecewise-linear power
    /// schedules, including every crossed knot. The actual native duration is
    /// retained if event-time subtraction differs by roundoff from `end-start`.
    /// Unsupported times, nonrepresentable positive powers/energies, or
    /// cancellation return no partial result. No source history is advanced.
    pub fn interval(
        &self,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SlabBeamInterval, DepositionError> {
        checkpoint()?;
        admit_interval(start_s, duration_s, end_s)?;
        let mut heat = [0.0; MAX_SLAB_LAYERS];
        let mut heat_correction = heat;
        let mut powers = [0.0; 6];
        let mut correction = powers;
        for beam in self.beams.iter() {
            let mean = schedule_mean(&beam.source.schedule, start_s, end_s, &mut checkpoint)?;
            for (i, &base) in beam.cell_power_w.iter().enumerate() {
                checkpoint()?;
                add(&mut heat[i], &mut heat_correction[i], product(base, mean)?);
            }
            for i in 0..6 {
                add(
                    &mut powers[i],
                    &mut correction[i],
                    product(beam.powers_w[i], mean)?,
                );
            }
        }
        let mut mean_sources_w = BTreeMap::new();
        let mut absorbed = 0.0;
        let mut correction = 0.0;
        for (i, &id) in self.cells.iter().enumerate() {
            if !heat[i].is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            add(&mut absorbed, &mut correction, heat[i]);
            if heat[i] != 0.0 {
                mean_sources_w.insert(id, heat[i]);
            }
        }
        let mut energies = [0.0; 6];
        for i in 0..6 {
            energies[i] = product(powers[i], duration_s)?;
        }
        let absorbed_j = product(absorbed, duration_s)?;
        let mut residual = 0.0;
        let mut correction = 0.0;
        for side in 0..2 {
            for value in [energies[side], -energies[2 + side], -energies[5 - side]] {
                add(&mut residual, &mut correction, value);
            }
        }
        add(&mut residual, &mut correction, -absorbed_j);
        if !residual.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok(SlabBeamInterval {
            duration_s,
            mean_sources_w,
            incident_j: [energies[0], energies[1]],
            reflected_j: [energies[2], energies[3]],
            transmitted_j: [energies[4], energies[5]],
            absorbed_j,
            balance_residual_j: residual,
        })
    }
}

fn admit_interval(start: f64, duration: f64, end: f64) -> Result<(), DepositionError> {
    if !start.is_finite()
        || start < 0.0
        || !end.is_finite()
        || end <= start
        || !duration.is_finite()
        || duration <= 0.0
        || ((end - start) - duration).abs() > 8.0 * f64::EPSILON * end.max(duration)
    {
        return Err(DepositionError::InvalidInput(
            "slab beams require the actual positive native interval",
        ));
    }
    Ok(())
}

fn schedule_mean(
    schedule: &HeatSourceSchedule,
    start: f64,
    end: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<f64, DepositionError> {
    let mut left = start;
    let mut a = schedule.multiplier_at(left)?;
    schedule.multiplier_at(end)?;
    let mut sum = 0.0;
    let mut correction = 0.0;
    let knots = schedule.knots();
    let mut upper = knots.partition_point(|knot| knot.time_s <= start);
    while left < end {
        checkpoint()?;
        let right = knots[upper].time_s.min(end);
        let b = schedule.multiplier_at(right)?;
        // Midpoint of nonnegative levels without an overflowing sum or losing
        // equal subnormal endpoints through two premature divisions by two.
        let average = a.min(b) + 0.5 * (a.max(b) - a.min(b));
        let weight = (right - left) / (end - start);
        if (average == 0.0 && a.max(b) > 0.0) || (weight == 0.0 && average > 0.0) {
            return Err(DepositionError::Arithmetic);
        }
        add(&mut sum, &mut correction, product(average, weight)?);
        left = right;
        a = b;
        upper += 1;
    }
    if sum.is_finite() {
        Ok(sum)
    } else {
        Err(DepositionError::Arithmetic)
    }
}

pub(super) fn add_beam_interval(
    beams: Option<&SlabBeamIllumination>,
    sources: &mut BTreeMap<CellId, f64>,
    start: f64,
    duration: f64,
    end: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    admit_interval(start, duration, end)?;
    if let Some(beams) = beams {
        for (id, power) in beams
            .interval(start, duration, end, checkpoint)?
            .mean_sources_w
        {
            let combined = sources.get(&id).copied().unwrap_or(0.0) + power;
            if !combined.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            if combined == 0.0 {
                sources.remove(&id);
            } else {
                sources.insert(id, combined);
            }
        }
    }
    Ok(())
}
