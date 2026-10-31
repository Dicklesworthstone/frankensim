//! Temperature-dependent beam absorption on the native enthalpy timeline.
//!
//! The original DDA is traversed once. At each accepted thermal substep the
//! original caloric laws supply each intersected cell's temperature; the
//! declared mass-absorption table supplies its spectrum. Serial Beer transport
//! therefore updates downstream shielding as upstream cells heat. Source rates
//! are frozen only over that one native explicit step, giving first-order time
//! accuracy. `advance_refined` compares complete native trajectories, including
//! latent enthalpy and outgoing radiation, before publishing a candidate.
//!
//! Geometry, density, source powers and spectral quadrature stay prescribed.
//! This is non-scattering incident-beam heating, without thermal emission,
//! deformation, chemistry, inferred material data or extrapolation. Source
//! identities bind caller-supplied data; they do not establish measurement
//! accuracy or prove that temperature alone is a complete material state.

use super::{
    CheckedSum, OpticalEnergyBudget, OpticalHeatingError, OpticalPowerBudget, finite, nonnegative,
};
use crate::volumes::transport::deposition::{BeamCellSegment, absorb_beam};
use crate::volumes::transport::{AbsorbingGrid, MassAbsorptionSpectrum, VolumeTransportError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;
use std::sync::Arc;

mod refinement;
pub use refinement::{FeedbackConvergence, FeedbackRefinement};

/// Explicitly source-bound samples of mass absorption versus temperature.
/// Interpolation is linear in coefficient, first in wavelength through the
/// original spectrum owner, then in kelvin. Neither coordinate extrapolates.
#[derive(Clone, Debug)]
pub struct TemperatureAbsorptionTable {
    source_identity: ContentHash,
    rows: Vec<(f64, MassAbsorptionSpectrum)>,
    wavelength_bounds_nm: (f64, f64),
}

impl TemperatureAbsorptionTable {
    pub fn try_new(
        source_identity: ContentHash,
        rows: &[(f64, MassAbsorptionSpectrum)],
        max_temperature_knots: usize,
    ) -> Result<Self, OpticalFeedbackError> {
        if rows.len() < 2
            || rows.len() > max_temperature_knots
            || source_identity.as_bytes().iter().all(|b| *b == 0)
            || rows.iter().any(|(t, _)| !t.is_finite() || *t <= 0.0)
            || rows.windows(2).any(|p| p[0].0 >= p[1].0)
        {
            return Err(invalid(
                "ordered positive temperatures, at least two bounded rows and an explicit source identity required",
            ));
        }
        let bounds = rows[0].1.wavelength_bounds_nm();
        if rows
            .iter()
            .any(|(_, spectrum)| spectrum.wavelength_bounds_nm() != bounds)
        {
            return Err(invalid(
                "temperature rows must share one declared wavelength domain",
            ));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(rows.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        retained.extend_from_slice(rows);
        Ok(Self {
            source_identity,
            rows: retained,
            wavelength_bounds_nm: bounds,
        })
    }

    pub fn temperature_bounds_k(&self) -> (f64, f64) {
        (self.rows[0].0, self.rows[self.rows.len() - 1].0)
    }

    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }

    pub fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
}

/// A finite collimated segment with already quadrature-weighted lane watts.
#[derive(Clone, Copy, Debug)]
pub struct IncidentBeam<const N: usize> {
    pub origin: Point3,
    pub unit_direction: Vec3,
    pub max_distance_m: f64,
    pub wavelengths_nm: [f64; N],
    pub incident_power_w: [f64; N],
}

/// Every non-vacuum traversed cell needs an explicit material and heat receiver,
/// including cells whose absorption happens to be zero at the current state.
#[derive(Clone, Copy, Debug)]
pub struct ThermalOpticalCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub absorption_table: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct OpticalFeedbackBudget {
    pub max_path_cells: usize,
    pub max_tables: usize,
    pub max_total_temperature_knots: usize,
    /// Whole native advance; counts each attempted source preparation.
    pub max_source_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum OpticalFeedbackError {
    Invalid(&'static str),
    Budget,
    TemperatureDomain {
        thermal_cell: CellId,
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    Transport(VolumeTransportError),
    Heating(OpticalHeatingError),
    Thermal(DepositionError),
    /// The complete native refinement sequence did not meet time accuracy.
    Convergence {
        normalized_error: f64,
        timestep_s: f64,
    },
}

impl From<DepositionError> for OpticalFeedbackError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<VolumeTransportError> for OpticalFeedbackError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Transport(error)
    }
}
impl From<OpticalHeatingError> for OpticalFeedbackError {
    fn from(error: OpticalHeatingError) -> Self {
        Self::Heating(error)
    }
}
impl core::fmt::Display for OpticalFeedbackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "temperature-dependent optical heating refused: {self:?}")
    }
}
impl std::error::Error for OpticalFeedbackError {}

#[derive(Debug)]
struct CompiledTable<const N: usize> {
    rows: Vec<(f64, [f64; N])>,
}

impl<const N: usize> CompiledTable<N> {
    fn check_temperature(
        &self,
        id: CellId,
        temperature_k: f64,
    ) -> Result<(), OpticalFeedbackError> {
        let minimum_k = self.rows[0].0;
        let maximum_k = self.rows[self.rows.len() - 1].0;
        if !temperature_k.is_finite() || !(minimum_k..=maximum_k).contains(&temperature_k) {
            return Err(OpticalFeedbackError::TemperatureDomain {
                thermal_cell: id,
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        Ok(())
    }

    fn at(&self, id: CellId, temperature_k: f64) -> Result<[f64; N], OpticalFeedbackError> {
        self.check_temperature(id, temperature_k)?;
        let upper = self
            .rows
            .partition_point(|(t, _)| *t < temperature_k)
            .clamp(1, self.rows.len() - 1);
        let (a, ka) = self.rows[upper - 1];
        let (b, kb) = self.rows[upper];
        let fraction = (temperature_k - a) / (b - a);
        let values = std::array::from_fn(|lane| {
            if ka[lane] == kb[lane] {
                ka[lane]
            } else {
                (1.0 - fraction) * ka[lane] + fraction * kb[lane]
            }
        });
        if !values.into_iter().all(nonnegative) {
            return Err(invalid("unrepresentable interpolated absorption"));
        }
        Ok(values)
    }
}

#[derive(Debug)]
struct BoundSegment {
    path: BeamCellSegment,
    thermal_cell: CellId,
    table: usize,
}

#[derive(Debug)]
struct FeedbackModel<const N: usize> {
    identity: ContentHash,
    beam: IncidentBeam<N>,
    wavelength_bounds_nm: (f64, f64),
    tables: Vec<CompiledTable<N>>,
    path: Vec<BoundSegment>,
    budget: OpticalFeedbackBudget,
    power_tolerance_w: f64,
}

/// Immutable optical data are shared by cheap candidate clones; accumulated
/// energy, clocks and prepared sources remain private to each candidate.
#[derive(Clone, Debug)]
pub struct TemperatureOpticalHeating<const N: usize> {
    model: Arc<FeedbackModel<N>>,
    time_s: f64,
    power: OpticalPowerBudget<N>,
    sources_w: BTreeMap<CellId, f64>,
    sums: EnergySums<N>,
    energy: OpticalEnergyBudget<N>,
    evaluations: usize,
}

impl<const N: usize> TemperatureOpticalHeating<N> {
    /// Cache the original physical DDA path and sample the original spectral
    /// laws once per temperature knot. Per-substep work is O(path cells × N)
    /// without spectral searches, geometry traversal or retained grid fields.
    /// The admitted grid supplies geometry/density; these explicit tables own
    /// absorption instead of its fixed spectrum. Any prescribed temperature
    /// field on that grid is not used: the receiving thermal cells supply the
    /// actual evolving state. This source accounts only for incident beam heat;
    /// a separate radiation owner must supply emission/cooling when required.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        grid: &AbsorbingGrid<'_>,
        beam: IncidentBeam<N>,
        tables: &[TemperatureAbsorptionTable],
        bindings: &[ThermalOpticalCell],
        thermal: &DepositionSimulation<'_>,
        budget: OpticalFeedbackBudget,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        if N == 0
            || thermal.cells().is_empty()
            || tables.is_empty()
            || tables.len() > budget.max_tables
            || bindings.len() > budget.max_path_cells
            || budget.max_source_evaluations == 0
            || !nonnegative(power_tolerance_w)
            || !beam.incident_power_w.into_iter().all(nonnegative)
        {
            return Err(invalid(
                "nonempty source tables, finite nonnegative powers/tolerance and explicit work limits required",
            ));
        }
        let path = grid.beam_cell_path(
            cx,
            beam.origin,
            beam.unit_direction,
            beam.max_distance_m,
            budget.max_path_cells,
        )?;
        let mut selected = BTreeMap::new();
        for binding in bindings {
            checkpoint()?;
            if selected.insert(binding.grid_cell, *binding).is_some() {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell).into());
            }
            if binding.absorption_table >= tables.len() {
                return Err(invalid("unknown temperature-dependent absorption table"));
            }
            if !thermal.cells().contains_key(&binding.thermal_cell) {
                return Err(OpticalHeatingError::InactiveCell(binding.thermal_cell).into());
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-beam.v1");
        hash.update(grid.content_identity().as_bytes());
        for value in [
            beam.origin.x,
            beam.origin.y,
            beam.origin.z,
            beam.unit_direction.x,
            beam.unit_direction.y,
            beam.unit_direction.z,
            beam.max_distance_m,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for lane in 0..N {
            hash.update(&beam.wavelengths_nm[lane].to_bits().to_le_bytes());
            hash.update(&beam.incident_power_w[lane].to_bits().to_le_bytes());
        }
        let mut compiled = Vec::new();
        compiled
            .try_reserve_exact(tables.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        let mut count = 0usize;
        let mut bounds = (0.0_f64, f64::MAX);
        for table in tables {
            checkpoint()?;
            count = count
                .checked_add(table.rows.len())
                .ok_or(OpticalFeedbackError::Budget)?;
            if count > budget.max_total_temperature_knots {
                return Err(OpticalFeedbackError::Budget);
            }
            bounds.0 = bounds.0.max(table.wavelength_bounds_nm.0);
            bounds.1 = bounds.1.min(table.wavelength_bounds_nm.1);
            hash.update(table.source_identity.as_bytes());
            hash.update(&table.wavelength_bounds_nm.0.to_bits().to_le_bytes());
            hash.update(&table.wavelength_bounds_nm.1.to_bits().to_le_bytes());
            hash.update(&(table.rows.len() as u64).to_le_bytes());
            let mut rows = Vec::new();
            rows.try_reserve_exact(table.rows.len())
                .map_err(|_| OpticalFeedbackError::Budget)?;
            for &(temperature, spectrum) in &table.rows {
                checkpoint()?;
                hash.update(&temperature.to_bits().to_le_bytes());
                hash.update(&[u8::from(spectrum.source_identity().is_some())]);
                if let Some(source) = spectrum.source_identity() {
                    hash.update(source.as_bytes());
                }
                let mut values = [0.0; N];
                for lane in 0..N {
                    values[lane] = spectrum.at(beam.wavelengths_nm[lane])?;
                    hash.update(&values[lane].to_bits().to_le_bytes());
                }
                rows.push((temperature, values));
            }
            compiled.push(CompiledTable { rows });
        }
        hash.update(&bounds.0.to_bits().to_le_bytes());
        hash.update(&bounds.1.to_bits().to_le_bytes());
        let mut bound_path = Vec::new();
        bound_path
            .try_reserve_exact(path.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        for segment in path {
            checkpoint()?;
            let binding = selected
                .remove(&segment.cell_index)
                .ok_or(OpticalHeatingError::UnmappedCell(segment.cell_index))?;
            hash.update(&(segment.cell_index as u64).to_le_bytes());
            hash.update(&binding.thermal_cell.0.to_le_bytes());
            hash.update(&(binding.absorption_table as u64).to_le_bytes());
            bound_path.push(BoundSegment {
                path: segment,
                thermal_cell: binding.thermal_cell,
                table: binding.absorption_table,
            });
        }
        if !selected.is_empty() {
            return Err(invalid(
                "binding names a cell outside the non-vacuum beam path",
            ));
        }
        let model = Arc::new(FeedbackModel {
            identity: hash.finalize(),
            beam,
            wavelength_bounds_nm: bounds,
            tables: compiled,
            path: bound_path,
            budget,
            power_tolerance_w,
        });
        let (power, sources_w) = model.evaluate(thermal, &mut checkpoint)?;
        let sums = EnergySums::default();
        let energy = sums.value()?;
        Ok(Self {
            model,
            time_s: thermal.time_s(),
            power,
            sources_w,
            sums,
            energy,
            evaluations: 0,
        })
    }

    /// Rates applied over the last accepted interval (initial rates before any
    /// advance). They are not a claim that the endpoint has the same opacity.
    pub fn power(&self) -> &OpticalPowerBudget<N> {
        &self.power
    }
    pub fn energy(&self) -> &OpticalEnergyBudget<N> {
        &self.energy
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }

    fn check_clock(&self, thermal: &DepositionSimulation<'_>) -> Result<(), OpticalFeedbackError> {
        if self.time_s != thermal.time_s() {
            return Err(invalid("optical and thermal accepted clocks differ"));
        }
        Ok(())
    }
}

impl<const N: usize> FeedbackModel<N> {
    fn evaluate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(OpticalPowerBudget<N>, BTreeMap<CellId, f64>), OpticalFeedbackError> {
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-beam-state.v1");
        hash.update(self.identity.as_bytes());
        let mut transmitted = self.beam.incident_power_w;
        let mut absorbed = [CheckedSum::default(); N];
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        for segment in &self.path {
            checkpoint()?;
            let temperature = thermal.temperature_k(segment.thermal_cell)?;
            let coefficients = self.tables[segment.table].at(segment.thermal_cell, temperature)?;
            hash.update(&temperature.to_bits().to_le_bytes());
            let mut heat = CheckedSum::default();
            for lane in 0..N {
                let (removed, remaining) = absorb_beam(
                    transmitted[lane],
                    segment.path.density_kg_m3,
                    segment.path.length_m,
                    coefficients[lane],
                );
                absorbed[lane].add(removed)?;
                heat.add(removed)?;
                transmitted[lane] = finite(remaining)?;
            }
            sources
                .entry(segment.thermal_cell)
                .or_default()
                .add(heat.value()?)?;
        }
        let mut power = OpticalPowerBudget {
            medium_identity: hash.finalize(),
            wavelength_bounds_nm: self.wavelength_bounds_nm,
            wavelengths_nm: self.beam.wavelengths_nm,
            incident_w: self.beam.incident_power_w,
            transmitted_w: transmitted,
            absorbed_w: [0.0; N],
            scattered_w: [0.0; N],
            residual_w: [0.0; N],
            thermal_input_w: 0.0,
        };
        let mut declared = CheckedSum::default();
        for lane in 0..N {
            power.absorbed_w[lane] = absorbed[lane].value()?;
            declared.add(power.absorbed_w[lane])?;
            let mut balance = CheckedSum::default();
            balance.add(power.incident_w[lane])?;
            balance.add(-transmitted[lane])?;
            balance.add(-power.absorbed_w[lane])?;
            power.residual_w[lane] = balance.value()?;
            if power.residual_w[lane].abs() > self.power_tolerance_w {
                return Err(invalid(
                    "temperature-dependent beam power is not conserved within roundoff allowance",
                ));
            }
        }
        let mut result = BTreeMap::new();
        let mut total = CheckedSum::default();
        for (id, heat) in sources {
            let watts = heat.value()?;
            total.add(watts)?;
            result.insert(id, watts);
        }
        power.thermal_input_w = total.value()?;
        if (power.thermal_input_w - declared.value()?).abs() > self.power_tolerance_w {
            return Err(invalid(
                "mapped heat disagrees with temperature-dependent absorption",
            ));
        }
        checkpoint()?;
        Ok((power, result))
    }
}

impl<const N: usize> DepositionObserver for TemperatureOpticalHeating<N> {
    type Error = OpticalFeedbackError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.evaluations = 0;
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        if self.evaluations >= self.model.budget.max_source_evaluations {
            return Err(OpticalFeedbackError::Budget);
        }
        self.evaluations += 1;
        let (power, sources) = self.model.evaluate(thermal, checkpoint)?;
        self.power = power;
        self.sources_w = sources;
        Ok(self.sources_w.clone())
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_clock(thermal)?;
        // Final-time events are observable even when no further source query
        // follows. Validate the same bindings/domains at that accepted state.
        for segment in &self.model.path {
            self.model.tables[segment.table].check_temperature(
                segment.thermal_cell,
                thermal.temperature_k(segment.thermal_cell)?,
            )?;
        }
        Ok(())
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || thermal.time_s() <= start_s
        {
            return Err(invalid(
                "matching accepted interval required for variable optical heat",
            ));
        }
        // Also admit the endpoint, including the final endpoint of an advance:
        // a coarse step must not escape a source's domain without a next query.
        for segment in &self.model.path {
            checkpoint()?;
            self.model.tables[segment.table].check_temperature(
                segment.thermal_cell,
                thermal.temperature_k(segment.thermal_cell)?,
            )?;
        }
        let mut sums = self.sums.clone();
        sums.add(&self.power, duration_s)?;
        let energy = sums.value()?;
        checkpoint()?;
        self.sums = sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct EnergySums<const N: usize> {
    duration: CheckedSum,
    incident: [CheckedSum; N],
    transmitted: [CheckedSum; N],
    absorbed: [CheckedSum; N],
    residual: [CheckedSum; N],
    thermal: CheckedSum,
}
impl<const N: usize> Default for EnergySums<N> {
    fn default() -> Self {
        Self {
            duration: CheckedSum::default(),
            incident: [CheckedSum::default(); N],
            transmitted: [CheckedSum::default(); N],
            absorbed: [CheckedSum::default(); N],
            residual: [CheckedSum::default(); N],
            thermal: CheckedSum::default(),
        }
    }
}
impl<const N: usize> EnergySums<N> {
    fn add(&mut self, power: &OpticalPowerBudget<N>, dt: f64) -> Result<(), OpticalFeedbackError> {
        self.duration.add(dt)?;
        self.thermal.add(finite(power.thermal_input_w * dt)?)?;
        for lane in 0..N {
            self.incident[lane].add(finite(power.incident_w[lane] * dt)?)?;
            self.transmitted[lane].add(finite(power.transmitted_w[lane] * dt)?)?;
            self.absorbed[lane].add(finite(power.absorbed_w[lane] * dt)?)?;
            self.residual[lane].add(finite(power.residual_w[lane] * dt)?)?;
        }
        Ok(())
    }
    fn value(&self) -> Result<OpticalEnergyBudget<N>, OpticalFeedbackError> {
        let mut energy = OpticalEnergyBudget {
            duration_s: self.duration.value()?,
            incident_j: [0.0; N],
            transmitted_j: [0.0; N],
            absorbed_j: [0.0; N],
            scattered_j: [0.0; N],
            residual_j: [0.0; N],
            thermal_input_j: self.thermal.value()?,
        };
        for lane in 0..N {
            energy.incident_j[lane] = self.incident[lane].value()?;
            energy.transmitted_j[lane] = self.transmitted[lane].value()?;
            energy.absorbed_j[lane] = self.absorbed[lane].value()?;
            energy.residual_j[lane] = self.residual[lane].value()?;
        }
        Ok(energy)
    }
}

fn invalid(reason: &'static str) -> OpticalFeedbackError {
    OpticalFeedbackError::Invalid(reason)
}
