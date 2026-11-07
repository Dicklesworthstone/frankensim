//! Depth-resolved cycle-mean acoustic heat on the original thermal timeline.
//!
//! Every thermal cell is an explicitly bound depth bin inside an ORIGINAL
//! physical solid layer. Exact positive strain-work integrals supply its bulk
//! and shear heat; no acoustic interfaces or midpoint heat samples are added.
//! The native enthalpy owner retains mass, latent storage, conduction, boundary
//! events, timestep selection and atomic acceptance. Absorbed watts enter its
//! source hook once, and accepted native durations own every joule.
//!
//! The coherent carrier, density, geometry and Maxwell laws remain frozen.
//! A declared temperature band only limits that caller-supplied approximation;
//! it does not infer temperature/phase-dependent moduli or melting acoustics.
//! The optional power envelope scales both coherent inputs by its square root.
//! This is quasi-steady cycle-mean heating: carrier transients, propagation
//! delay, resonator ring-up, mechanical storage changes and acoustic streaming
//! are not resolved. Envelopes must be slow relative to those acoustic scales.
//! Thermal bins still need spatial/time refinement under the original native
//! cell-uniformity and supplied-conductance assumptions.

use super::{
    ViscoelasticSlabError,
    layers::{ViscoelasticLayers, ViscoelasticLayersPower},
};
use fs_blake3::ContentHash;
use fs_math::{c64::C64, dd::Dd};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalStorageBasis,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Hard storage bound for resolved acoustic depth bins.
pub const MAX_ACOUSTIC_THERMAL_BINS: usize = 16_384;
/// Hard storage bound for the optional prescribed power envelope.
pub const MAX_ACOUSTIC_POWER_SEGMENTS: usize = 4_096;

/// One original physical layer interval and its original thermal control volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticThermalBin {
    /// Physical layer in the prepared acoustic stack.
    pub layer: usize,
    /// Ordered local depths [m]; bins must exactly partition each positive layer.
    pub depths_m: [f64; 2],
    /// Unique active native thermal cell for this interval.
    pub cell: CellId,
    /// Closed positive temperature range [K] where the caller admits the SAME
    /// frozen acoustic law. This is not a calibrated temperature-response model.
    pub temperature_band_k: [f64; 2],
}

/// One affine segment of a prescribed dimensionless POWER fraction in [0,1].
/// Outside supplied segments, including gaps, the source is dark. Both carrier
/// phasors scale by sqrt(fraction), so their relative phase stays coherent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticPowerSegment {
    /// Absolute segment start [s], inclusive.
    pub start_s: f64,
    /// Absolute segment end [s], exclusive; must exceed start.
    pub end_s: f64,
    /// Power fraction at the start.
    pub start_fraction: f64,
    /// Left-limit power fraction at the end; jumps between segments are allowed.
    pub end_fraction: f64,
}

/// Explicit storage, work and numerical-account allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticHeatingBudget {
    pub maximum_bins: usize,
    pub maximum_power_segments: usize,
    /// Accepted interval preparations in one advance; reset on native begin.
    pub maximum_source_evaluations: usize,
    /// Whole-history bound [J] on accumulated absolute boundary/local-heat and
    /// local-heat/thermal-credit discrepancies, plus native clock discrepancy.
    /// This admits arithmetic only, not physical or discretization error.
    pub energy_absolute_tolerance_j: f64,
}

/// One cached physical bin contribution at the authored peak coherent drive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticBinPower {
    pub binding: AcousticThermalBin,
    pub source_identity: ContentHash,
    pub bulk_absorbed_power_w: f64,
    pub shear_absorbed_power_w: f64,
}

/// Cumulative accepted cycle-mean acoustic energy, distinct from the native
/// thermal audit's other boundaries, births and internal exchanges.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AcousticHeatingEnergy {
    pub duration_s: f64,
    /// Integral of the prescribed power fraction over accepted native time [s].
    pub exposure_s: f64,
    pub incident_j: [f64; 2],
    pub outgoing_j: [f64; 2],
    pub bulk_absorbed_j: f64,
    pub shear_absorbed_j: f64,
    /// Actual acoustic source watts times the duration credited by the native
    /// thermal owner. Already included in its external_heat_j; do not add twice.
    pub thermal_input_j: f64,
    /// Independent incident minus outgoing minus constitutive absorption [J].
    pub balance_residual_j: f64,
    /// Constitutive absorption minus actual cell source credits [J].
    pub credit_residual_j: f64,
    /// Sum of absolute per-interval balance and credit discrepancies [J].
    pub accumulated_absolute_defect_j: f64,
    /// Peak incident power times accumulated native endpoint/duration mismatch [J].
    pub timing_discrepancy_j: f64,
}

/// A refusal publishes neither new enthalpy nor an acoustic energy account.
#[derive(Debug)]
pub enum AcousticHeatingError {
    Invalid(&'static str),
    Acoustic(ViscoelasticSlabError),
    Thermal(DepositionError),
    Budget,
    ClockMismatch {
        thermal_time_s: f64,
        acoustic_time_s: f64,
    },
    TemperatureSupport {
        cell: CellId,
        temperature_k: f64,
        band_k: [f64; 2],
    },
}
impl core::fmt::Display for AcousticHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "acoustic heating refused: {self:?}")
    }
}
impl core::error::Error for AcousticHeatingError {}
impl From<ViscoelasticSlabError> for AcousticHeatingError {
    fn from(value: ViscoelasticSlabError) -> Self {
        Self::Acoustic(value)
    }
}
impl From<DepositionError> for AcousticHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}

type Result<T> = core::result::Result<T, AcousticHeatingError>;

#[derive(Clone, Debug)]
struct BoundBin {
    power: AcousticBinPower,
    material_id: MaterialId,
    segment_id: u64,
    body_id: u64,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric_identity: Option<[u8; 32]>,
    storage_basis: ThermalStorageBasis,
    heat_w: f64,
}

impl BoundBin {
    fn validate(&self, thermal: &DepositionSimulation<'_>) -> Result<()> {
        let cell =
            thermal
                .cells()
                .get(&self.power.binding.cell)
                .ok_or(AcousticHeatingError::Invalid(
                    "acoustic thermal cell is inactive",
                ))?;
        if cell.material_id != self.material_id
            || cell.segment_id != self.segment_id
            || cell.body_id != self.body_id
            || cell.mass_kg != self.mass_kg
            || cell.geometry != self.geometry
            || thermal.material_caloric_identity(cell.material_id)? != self.caloric_identity
            || thermal.material_storage_basis(cell.material_id)? != self.storage_basis
        {
            return Err(AcousticHeatingError::Invalid(
                "acoustic cell material, mass, geometry or original caloric binding changed",
            ));
        }
        let temperature_k = thermal.temperature_k(cell.id)?;
        let [lo, hi] = self.power.binding.temperature_band_k;
        if temperature_k < lo || temperature_k > hi {
            return Err(AcousticHeatingError::TemperatureSupport {
                cell: cell.id,
                temperature_k,
                band_k: [lo, hi],
            });
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct Account {
    duration: Dd,
    exposure: Dd,
    incident: [Dd; 2],
    outgoing: [Dd; 2],
    heat: [Dd; 2],
    credit: Dd,
    absolute_defect: Dd,
    timing: Dd,
}
impl Default for Account {
    fn default() -> Self {
        Self {
            duration: Dd::ZERO,
            exposure: Dd::ZERO,
            incident: [Dd::ZERO; 2],
            outgoing: [Dd::ZERO; 2],
            heat: [Dd::ZERO; 2],
            credit: Dd::ZERO,
            absolute_defect: Dd::ZERO,
            timing: Dd::ZERO,
        }
    }
}
impl Account {
    fn value(&self) -> Result<AcousticHeatingEnergy> {
        let heat = self.heat[0] + self.heat[1];
        let balance =
            self.incident[0] + self.incident[1] - self.outgoing[0] - self.outgoing[1] - heat;
        let values = [
            self.duration,
            self.exposure,
            self.incident[0],
            self.incident[1],
            self.outgoing[0],
            self.outgoing[1],
            self.heat[0],
            self.heat[1],
            self.credit,
            self.absolute_defect,
            self.timing,
        ];
        if values.iter().any(|v| !nonnegative(v.to_f64()))
            || !balance.to_f64().is_finite()
            || !(heat - self.credit).to_f64().is_finite()
        {
            return Err(AcousticHeatingError::Invalid(
                "acoustic energy account is unrepresentable",
            ));
        }
        Ok(AcousticHeatingEnergy {
            duration_s: self.duration.to_f64(),
            exposure_s: self.exposure.to_f64(),
            incident_j: self.incident.map(Dd::to_f64),
            outgoing_j: self.outgoing.map(Dd::to_f64),
            bulk_absorbed_j: self.heat[0].to_f64(),
            shear_absorbed_j: self.heat[1].to_f64(),
            thermal_input_j: self.credit.to_f64(),
            balance_residual_j: balance.to_f64(),
            credit_residual_j: (heat - self.credit).to_f64(),
            accumulated_absolute_defect_j: self.absolute_defect.to_f64(),
            timing_discrepancy_j: self.timing.to_f64(),
        })
    }
}

#[derive(Clone, Debug)]
struct Pending {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    fraction: f64,
    account: Account,
    energy: AcousticHeatingEnergy,
}

/// Prepared spatial acoustic heat source, borrowing the original solid stack.
/// Immutable constitutive powers are prepared once. Native steps make bounded
/// passes over the bins without repeating acoustic wave/material/exterior
/// queries. Native ordered-cell lookups and source-map insertions cost O(log n)
/// per bin; the original observer protocol allocates that per-step source map.
#[derive(Clone)]
pub struct ViscoelasticHeating<'s, 'a> {
    stack: &'s ViscoelasticLayers<'a>,
    incident_peak_pressure_pa: [C64; 2],
    area_m2: f64,
    bins: Arc<[BoundBin]>,
    peak_power: ViscoelasticLayersPower,
    envelope: Option<Arc<[AcousticPowerSegment]>>,
    budget: AcousticHeatingBudget,
    time_s: f64,
    last_fraction: f64,
    evaluations: usize,
    account: Account,
    energy: AcousticHeatingEnergy,
    pending: Option<Pending>,
}

impl<'s, 'a> ViscoelasticHeating<'s, 'a> {
    /// Bind a complete nonoverlapping depth partition to the active native cells.
    /// Cell volume must equal area*depth and mass must equal that volume times
    /// the ORIGINAL acoustic density, within 64 eps relative arithmetic.
    /// No mass, density, caloric law or acoustic material is replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        stack: &'s ViscoelasticLayers<'a>,
        incident_peak_pressure_pa: [C64; 2],
        area_m2: f64,
        bins: &[AcousticThermalBin],
        thermal: &DepositionSimulation<'_>,
        budget: AcousticHeatingBudget,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Self> {
        checkpoint()?;
        if budget.maximum_bins > MAX_ACOUSTIC_THERMAL_BINS
            || budget.maximum_power_segments > MAX_ACOUSTIC_POWER_SEGMENTS
            || bins.len() > budget.maximum_bins
            || !nonnegative(budget.energy_absolute_tolerance_j)
        {
            return Err(AcousticHeatingError::Budget);
        }
        if bins.is_empty() || !positive(area_m2) {
            return Err(AcousticHeatingError::Invalid(
                "acoustic heating needs positive area and thermal bins",
            ));
        }
        let mut ordered = bins.to_vec();
        ordered.sort_by(|a, b| {
            a.layer
                .cmp(&b.layer)
                .then(a.depths_m[0].total_cmp(&b.depths_m[0]))
        });
        let mut identities = BTreeSet::new();
        for bin in &ordered {
            checkpoint()?;
            let [lo, hi] = bin.temperature_band_k;
            if bin.layer >= stack.layer_count()
                || !positive(lo)
                || !positive(hi)
                || hi < lo
                || !identities.insert(bin.cell)
            {
                return Err(AcousticHeatingError::Invalid(
                    "acoustic bin index, unique cell or temperature support is invalid",
                ));
            }
        }
        let mut interruption = None;
        let drive = stack.drive(incident_peak_pressure_pa, area_m2, &mut || {
            cancelled(&mut checkpoint, &mut interruption)
        });
        let drive = acoustic_result(drive, interruption.take())?;
        let mut bound = Vec::with_capacity(ordered.len());
        let mut cursor = 0;
        for layer_index in 0..stack.layer_count() {
            checkpoint()?;
            let layer = stack.layer(layer_index).expect("original layer");
            let mut depth = 0.0;
            let mut heat = [Dd::ZERO; 2];
            let first = cursor;
            while cursor < ordered.len() && ordered[cursor].layer == layer_index {
                checkpoint()?;
                let bin = ordered[cursor];
                let [start, end] = bin.depths_m;
                if !nonnegative(start)
                    || !positive(end)
                    || start != depth
                    || end <= start
                    || end > layer.thickness_m
                {
                    return Err(AcousticHeatingError::Invalid(
                        "acoustic depth bins must partition each layer without gaps or overlap",
                    ));
                }
                let cell = thermal
                    .cells()
                    .get(&bin.cell)
                    .ok_or(AcousticHeatingError::Invalid(
                        "acoustic thermal bin is inactive",
                    ))?;
                let volume = product(area_m2, end - start)?;
                let mass = product(layer.source.density_kg_m3, volume)?;
                if !roundoff_equal(cell.geometry.volume_m3, volume)
                    || !roundoff_equal(cell.mass_kg, mass)
                {
                    return Err(AcousticHeatingError::Invalid(
                        "acoustic bin volume or mass disagrees with the original solid",
                    ));
                }
                let local = drive.absorption_between(layer_index, [start, end], &mut || {
                    cancelled(&mut checkpoint, &mut interruption)
                });
                let local = acoustic_result(local, interruption.take())?;
                let q = [local.bulk_absorbed_power_w, local.shear_absorbed_power_w];
                let heat_w = checked((Dd::from_f64(q[0]) + Dd::from_f64(q[1])).to_f64())?;
                for channel in 0..2 {
                    heat[channel] = heat[channel] + Dd::from_f64(q[channel]);
                }
                let binding = BoundBin {
                    power: AcousticBinPower {
                        binding: bin,
                        source_identity: local.source_identity,
                        bulk_absorbed_power_w: q[0],
                        shear_absorbed_power_w: q[1],
                    },
                    material_id: cell.material_id,
                    segment_id: cell.segment_id,
                    body_id: cell.body_id,
                    mass_kg: cell.mass_kg,
                    geometry: cell.geometry,
                    caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
                    storage_basis: thermal.material_storage_basis(cell.material_id)?,
                    heat_w,
                };
                binding.validate(thermal)?;
                bound.push(binding);
                depth = end;
                cursor += 1;
            }
            if depth != layer.thickness_m {
                return Err(AcousticHeatingError::Invalid(
                    "every positive acoustic layer needs a complete thermal depth partition",
                ));
            }
            let original = drive.layer_absorption(layer_index, &mut || {
                cancelled(&mut checkpoint, &mut interruption)
            });
            let original = acoustic_result(original, interruption.take())?;
            for (sum, total) in heat.into_iter().zip([
                original.bulk_absorbed_power_w,
                original.shear_absorbed_power_w,
            ]) {
                let discrepancy = (sum - Dd::from_f64(total)).to_f64().abs();
                let scale = sum.to_f64().max(total);
                let allowance = 4096.0 * f64::EPSILON * (cursor - first + 1) as f64 * scale;
                if !nonnegative(discrepancy) || discrepancy > allowance {
                    return Err(AcousticHeatingError::Invalid(
                        "positive depth-bin heat does not close against the original layer",
                    ));
                }
            }
        }
        if cursor != ordered.len() {
            return Err(AcousticHeatingError::Invalid("unbound acoustic depth bin"));
        }
        checkpoint()?;
        Ok(Self {
            stack,
            incident_peak_pressure_pa,
            area_m2,
            bins: bound.into(),
            peak_power: drive.power(),
            envelope: None,
            budget,
            time_s: thermal.time_s(),
            last_fraction: 1.0,
            evaluations: 0,
            account: Account::default(),
            energy: AcousticHeatingEnergy::default(),
            pending: None,
        })
    }

    /// Bind a bounded affine POWER envelope before this source's first advance.
    /// Constant operation remains the default. Gaps and times after the final
    /// segment are dark; native steps clip at every jump and ramp endpoint.
    /// This prescribes a quasi-steady intensity, not carrier/ring-up dynamics.
    pub fn with_power_envelope(
        mut self,
        thermal: &DepositionSimulation<'_>,
        segments: &[AcousticPowerSegment],
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Self> {
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        self.no_pending()?;
        if self.energy.duration_s != 0.0 {
            return Err(AcousticHeatingError::Invalid(
                "bind the acoustic power envelope before the first advance",
            ));
        }
        if segments.len() > self.budget.maximum_power_segments {
            return Err(AcousticHeatingError::Budget);
        }
        let mut previous = 0.0;
        for s in segments {
            checkpoint()?;
            if !nonnegative(s.start_s)
                || !positive(s.end_s)
                || s.end_s <= s.start_s
                || s.start_s < previous
                || !fraction(s.start_fraction)
                || !fraction(s.end_fraction)
            {
                return Err(AcousticHeatingError::Invalid(
                    "acoustic power envelope needs ordered bounded affine segments",
                ));
            }
            previous = s.end_s;
        }
        self.last_fraction = fraction_at(segments, self.time_s);
        self.envelope = Some(segments.into());
        checkpoint()?;
        Ok(self)
    }

    pub const fn original_stack(&self) -> &'s ViscoelasticLayers<'a> {
        self.stack
    }
    pub const fn incident_peak_pressure_pa(&self) -> [C64; 2] {
        self.incident_peak_pressure_pa
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn peak_power(&self) -> ViscoelasticLayersPower {
        self.peak_power
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn energy(&self) -> &AcousticHeatingEnergy {
        &self.energy
    }
    pub const fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    /// Mean multiplier applied on the last accepted interval; initial
    /// instantaneous multiplier before any advance.
    pub const fn last_power_fraction(&self) -> f64 {
        self.last_fraction
    }
    pub fn bin_count(&self) -> usize {
        self.bins.len()
    }
    pub fn bin_power(&self, index: usize) -> Option<&AcousticBinPower> {
        self.bins.get(index).map(|bin| &bin.power)
    }
    pub fn power_envelope(&self) -> Option<&[AcousticPowerSegment]> {
        self.envelope.as_deref()
    }

    /// Advance the original enthalpy history and this acoustic account together.
    /// Other native boundary sources and conductive contacts remain independent;
    /// the returned thermal external_heat_j already includes the acoustic credit.
    pub fn advance_to<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'thermal>, Self, EnergyAudit)> {
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        self.no_pending()?;
        let result = thermal.advance_observed(self, target_s, config, &mut checkpoint)?;
        result
            .1
            .validate(&result.0, result.1.time_s, &mut checkpoint)?;
        result.1.admit_energy(&result.1.energy)?;
        Ok(result)
    }

    fn no_pending(&self) -> Result<()> {
        if self.pending.is_some() {
            Err(AcousticHeatingError::Invalid(
                "uncommitted acoustic heating interval",
            ))
        } else {
            Ok(())
        }
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<()> {
        checkpoint()?;
        if thermal.time_s() != time_s {
            return Err(AcousticHeatingError::ClockMismatch {
                thermal_time_s: thermal.time_s(),
                acoustic_time_s: time_s,
            });
        }
        for bin in self.bins.iter() {
            checkpoint()?;
            bin.validate(thermal)?;
        }
        Ok(())
    }

    fn admit_energy(&self, energy: &AcousticHeatingEnergy) -> Result<()> {
        if energy.accumulated_absolute_defect_j + energy.timing_discrepancy_j
            > self.budget.energy_absolute_tolerance_j
        {
            Err(AcousticHeatingError::Invalid(
                "acoustic energy arithmetic exceeds the declared joule allowance",
            ))
        } else {
            Ok(())
        }
    }
}

impl DepositionObserver for ViscoelasticHeating<'_, '_> {
    type Error = AcousticHeatingError;

    fn begin(&mut self) -> Result<()> {
        self.no_pending()?;
        self.evaluations = 0;
        Ok(())
    }

    fn requires_step_interval(&self) -> bool {
        true
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Option<f64>> {
        checkpoint()?;
        if thermal.time_s() != self.time_s {
            return Err(AcousticHeatingError::ClockMismatch {
                thermal_time_s: thermal.time_s(),
                acoustic_time_s: self.time_s,
            });
        }
        self.no_pending()?;
        Ok(self
            .envelope
            .as_deref()
            .and_then(|s| next_transition(s, self.time_s)))
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>> {
        self.validate(thermal, self.time_s, checkpoint)?;
        self.no_pending()?;
        if self.evaluations >= self.budget.maximum_source_evaluations {
            return Err(AcousticHeatingError::Budget);
        }
        if start_s != self.time_s || !positive(duration_s) || !positive(end_s) || end_s <= start_s {
            return Err(AcousticHeatingError::Invalid(
                "invalid native acoustic heating interval",
            ));
        }
        let mean = match &self.envelope {
            Some(segments) => interval_fraction(segments, start_s, end_s)?,
            None => 1.0,
        };
        let mut sources = BTreeMap::new();
        let mut credit = Dd::ZERO;
        for bin in self.bins.iter() {
            checkpoint()?;
            let watts = product(bin.heat_w, mean)?;
            sources.insert(bin.power.binding.cell, watts);
            credit = credit + Dd::from_f64(product(watts, duration_s)?);
        }
        let energy_pair = |powers: [f64; 2]| -> Result<[f64; 2]> {
            Ok([
                product(product(powers[0], mean)?, duration_s)?,
                product(product(powers[1], mean)?, duration_s)?,
            ])
        };
        let incident = energy_pair(self.peak_power.incident_power_w)?;
        let outgoing = energy_pair(self.peak_power.outgoing_power_w)?;
        let heat = energy_pair([
            self.peak_power.bulk_absorbed_power_w,
            self.peak_power.shear_absorbed_power_w,
        ])?;
        let physical_heat = Dd::from_f64(heat[0]) + Dd::from_f64(heat[1]);
        let balance = Dd::from_f64(incident[0]) + Dd::from_f64(incident[1])
            - Dd::from_f64(outgoing[0])
            - Dd::from_f64(outgoing[1])
            - physical_heat;
        let mut account = self.account.clone();
        add(&mut account.duration, duration_s)?;
        add(&mut account.exposure, product(mean, duration_s)?)?;
        for i in 0..2 {
            add(&mut account.incident[i], incident[i])?;
            add(&mut account.outgoing[i], outgoing[i])?;
            add(&mut account.heat[i], heat[i])?;
        }
        account.credit = account.credit + credit;
        add(&mut account.absolute_defect, balance.to_f64().abs())?;
        add(
            &mut account.absolute_defect,
            (physical_heat - credit).to_f64().abs(),
        )?;
        let mismatch = (Dd::from_f64(end_s) - Dd::from_f64(start_s) - Dd::from_f64(duration_s))
            .to_f64()
            .abs();
        let peak_incident = checked(
            (Dd::from_f64(self.peak_power.incident_power_w[0])
                + Dd::from_f64(self.peak_power.incident_power_w[1]))
            .to_f64(),
        )?;
        add(&mut account.timing, checked(peak_incident * mismatch)?)?;
        let energy = account.value()?;
        self.admit_energy(&energy)?;
        checkpoint()?;
        self.pending = Some(Pending {
            start_s,
            duration_s,
            end_s,
            fraction: mean,
            account,
            energy,
        });
        Ok(sources)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<()> {
        self.no_pending()?;
        self.validate(thermal, self.time_s, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<()> {
        let pending = self.pending.as_ref().ok_or(AcousticHeatingError::Invalid(
            "native acoustic step has no prepared interval",
        ))?;
        if start_s != self.time_s || pending.start_s != start_s || pending.duration_s != duration_s
        {
            return Err(AcousticHeatingError::Invalid(
                "native acoustic interval changed before acceptance",
            ));
        }
        self.validate(thermal, pending.end_s, checkpoint)?;
        self.admit_energy(&pending.energy)?;
        checkpoint()?;
        let pending = self.pending.take().expect("validated acoustic interval");
        self.time_s = pending.end_s;
        self.last_fraction = pending.fraction;
        self.account = pending.account;
        self.energy = pending.energy;
        self.evaluations += 1;
        Ok(())
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn fraction(value: f64) -> bool {
    nonnegative(value) && value <= 1.0
}
fn checked(value: f64) -> Result<f64> {
    if nonnegative(value) {
        Ok(value)
    } else {
        Err(AcousticHeatingError::Invalid(
            "acoustic power or energy is unrepresentable",
        ))
    }
}
fn product(a: f64, b: f64) -> Result<f64> {
    let value = checked(a * b)?;
    if value == 0.0 && a > 0.0 && b > 0.0 {
        Err(AcousticHeatingError::Invalid(
            "positive acoustic power or energy underflows",
        ))
    } else {
        Ok(value)
    }
}
fn add(sum: &mut Dd, value: f64) -> Result<()> {
    checked(value)?;
    *sum = *sum + Dd::from_f64(value);
    checked(sum.to_f64())?;
    Ok(())
}
fn roundoff_equal(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}
fn cancelled(
    checkpoint: &mut dyn FnMut() -> core::result::Result<(), DepositionError>,
    interruption: &mut Option<DepositionError>,
) -> bool {
    match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            *interruption = Some(error);
            true
        }
    }
}
fn acoustic_result<T>(
    result: core::result::Result<T, ViscoelasticSlabError>,
    interruption: Option<DepositionError>,
) -> Result<T> {
    match interruption {
        Some(error) => Err(error.into()),
        None => result.map_err(Into::into),
    }
}
fn next_transition(segments: &[AcousticPowerSegment], time_s: f64) -> Option<f64> {
    let index = segments.partition_point(|s| s.end_s <= time_s);
    segments.get(index).map(|s| {
        if time_s < s.start_s {
            s.start_s
        } else {
            s.end_s
        }
    })
}
fn segment_fraction(segment: &AcousticPowerSegment, time_s: f64) -> f64 {
    let x = (time_s - segment.start_s) / (segment.end_s - segment.start_s);
    let a = segment.start_fraction;
    let b = segment.end_fraction;
    if x == 0.0 || a == b {
        a
    } else if x == 1.0 {
        b
    } else if a < b {
        a + (b - a) * x
    } else {
        b + (a - b) * (1.0 - x)
    }
}
fn fraction_at(segments: &[AcousticPowerSegment], time_s: f64) -> f64 {
    let index = segments.partition_point(|s| s.end_s <= time_s);
    segments
        .get(index)
        .filter(|s| s.start_s <= time_s)
        .map_or(0.0, |s| segment_fraction(s, time_s))
}
fn interval_fraction(segments: &[AcousticPowerSegment], start_s: f64, end_s: f64) -> Result<f64> {
    if next_transition(segments, start_s).is_some_and(|t| t < end_s) {
        return Err(AcousticHeatingError::Invalid(
            "native acoustic interval crosses a power transition",
        ));
    }
    let index = segments.partition_point(|s| s.end_s <= start_s);
    let Some(s) = segments.get(index).filter(|s| s.start_s <= start_s) else {
        return Ok(0.0);
    };
    let a = segment_fraction(s, start_s);
    let b = segment_fraction(s, end_s);
    let mean = a.min(b) + 0.5 * (a - b).abs();
    if mean == 0.0 && a.max(b) > 0.0 {
        return Err(AcousticHeatingError::Invalid(
            "positive acoustic envelope mean underflows",
        ));
    }
    Ok(mean)
}
