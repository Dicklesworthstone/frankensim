//! Data-only assembly of the original NSF volume, DG bore and pressure observer.
//!
//! Geometry is an explicit Cartesian grid and stationary solid-cell mask. Gas
//! evolution retains either the original dry-air law or complete thermochemical
//! source models; a mixture snapshot is never converted into an invented law.
//! Initial physical states and finite pressure/characteristic controls are data. No jet
//! card, noise generator, desired pitch or instrument name enters this path.
//! See `docs/RESOLVED_NSF_DG_COUPLING.md` for the grammar and physical limits.

use super::{NsfDgPerformance, ReservoirObservation, source_intervals_for};
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::dg_port::input::{DgSourceError, DgSourceInput};
use crate::dg_port::input::wall::regions::{DgMaterialWallLayout, DgMaterialWallLayoutInput};
use crate::nsf_dg::{NsfDg, NsfDgDrive, NsfDgError, NsfDgLimits, same_clock};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::pcm_wav::baffled::gas::GasPropagation;
use crate::render::RenderError;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::CancelGate;
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfBoundaryDrive, NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfCharacteristicPatch,
    NsfError, NsfGasMedium, NsfGrid3, NsfOpenState, NsfReservoirAperture, NsfSgsConfig, NsfSpec,
    ThermochemicalNsfGas, WallThermal,
};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::input::{
    GasInputError, GasRecord, ThermochemicalSources, read_frozen_thermochemical_model,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use std::sync::Arc;

pub const RESOLVED_DG_PERFORMANCE_SCHEMA: &str = "frankensim-resolved-dg-performance-v1";
pub const MAX_RESOLVED_PERFORMANCE_BYTES: usize = 1024 * 1024;
/// Includes stationary solid cells and all cells retained by each trial.
pub const MAX_RESOLVED_CELLS: usize = 262_144;
pub const MAX_RESOLVED_SAMPLES: usize = 2_880_000;

/// Refusals retain the original physical owner; parsing does not repair inputs.
#[derive(Debug)]
pub enum ResolvedDgInputError {
    Input { line: usize, what: &'static str },
    Gas(GasInputError),
    Flow(NsfError),
    Duct(DgSourceError),
    Coupling(NsfDgError),
    Render(RenderError),
    Material { line: usize, source: String },
}
impl core::fmt::Display for ResolvedDgInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved flow input: {self:?}")
    }
}
impl std::error::Error for ResolvedDgInputError {}
impl From<DgSourceError> for ResolvedDgInputError {
    fn from(value: DgSourceError) -> Self { Self::Duct(value) }
}
fn bad(line: usize, what: &'static str) -> ResolvedDgInputError {
    ResolvedDgInputError::Input { line, what }
}
fn poll(gate: &CancelGate) -> Result<(), ResolvedDgInputError> {
    if gate.is_requested() { Err(bad(0, "resolved source construction cancelled")) } else { Ok(()) }
}

struct Rows<'a> { lines: Lines<'a>, line: usize }
struct Row<'a> { fields: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Rows<'a> {
    fn next_is(&self, key: &str) -> bool {
        self.lines.clone().next().is_some_and(|s| s.split_ascii_whitespace().next() == Some(key))
    }
    fn next(&mut self) -> Result<(usize, SplitAsciiWhitespace<'a>), ResolvedDgInputError> {
        self.line += 1;
        let text = self.lines.next().ok_or_else(|| bad(self.line, "missing source record"))?;
        Ok((self.line, text.split_ascii_whitespace()))
    }
    fn row(&mut self, key: &'static str) -> Result<Row<'a>, ResolvedDgInputError> {
        let (line, mut fields) = self.next()?;
        if fields.next() != Some(key) { return Err(bad(line, key)); }
        Ok(Row { fields, line })
    }
}
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str, ResolvedDgInputError> {
        self.fields.next().ok_or_else(|| bad(self.line, "missing source field"))
    }
    fn number<T: FromStr>(&mut self) -> Result<T, ResolvedDgInputError> {
        self.word()?.parse().map_err(|_| bad(self.line, "invalid numeric field"))
    }
    fn scalar(&mut self) -> Result<f64, ResolvedDgInputError> {
        let value: f64 = self.number()?;
        if !value.is_finite() { return Err(bad(self.line, "source scalars must be finite")); }
        Ok(value)
    }
    fn flag(&mut self) -> Result<bool, ResolvedDgInputError> {
        match self.word()? { "0" => Ok(false), "1" => Ok(true),
            _ => Err(bad(self.line, "periodicity is explicitly 0 or 1")) }
    }
    fn end(mut self) -> Result<(), ResolvedDgInputError> {
        if self.fields.next().is_some() { Err(bad(self.line, "unexpected trailing fields")) } else { Ok(()) }
    }
    fn index(&mut self, dims: [usize; 3]) -> Result<usize, ResolvedDgInputError> {
        let x: [usize; 3] = [self.number()?, self.number()?, self.number()?];
        if x.iter().zip(dims).any(|(&i, n)| i >= n) {
            return Err(bad(self.line, "cell coordinate outside the physical grid"));
        }
        Ok((x[2] * dims[1] + x[1]) * dims[0] + x[0])
    }
    fn state(&mut self) -> Result<[f64; 5], ResolvedDgInputError> {
        let value = [self.scalar()?, self.scalar()?, self.scalar()?, self.scalar()?, self.scalar()?];
        if value[0] <= 0.0 || value[4] <= 0.0 {
            return Err(bad(self.line, "physical density and absolute temperature must be positive"));
        }
        Ok(value)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ResolvedDgPerformanceInfo {
    pub input_hash: ContentHash,
    pub sample_rate_hz: u32,
    pub samples: u64,
    pub full_scale_pa: f64,
    pub dimensions: [usize; 3],
    pub solid_cells: usize,
    pub initial_state_cells: usize,
    pub boundary_segments: usize,
    /// Original physical NSF/DG intervals required for the observation window.
    pub physical_intervals: u64,
    /// Explicit exploratory transport declaration; zero coefficients retain
    /// the molecular trajectory but remain visible as a declared control.
    pub sgs: Option<NsfSgsConfig>,
}

/// The retained owners, not a baked recording or separate synthesis engine.
pub struct ResolvedDgFilePerformance {
    info: ResolvedDgPerformanceInfo,
    renderer: NsfDgPerformance<D3q27>,
    gas_records: Vec<GasRecord>,
    wall: DgMaterialWallLayout,
    upstream_wall: Option<DgMaterialWallLayout>,
}
impl ResolvedDgFilePerformance {
    pub fn from_bytes(bytes: &[u8], maximum_block: usize, gate: Arc<CancelGate>)
        -> Result<Self, ResolvedDgInputError>
    {
        Self::from_bytes_with_sources(bytes, maximum_block, gate, None,
            |_| Err("solid duct walls require an explicit material source loader".into()))
    }

    /// Source callbacks are supplied by the caller's existing attachment/store
    /// boundary. This L3 reader opens no files and retains the same cancellation
    /// gate around bounded construction kernels, every physical interval and
    /// observation. The original NSF cold constructor itself has no internal
    /// cancellation callback; its finite cell/transport budgets still apply.
    pub fn from_bytes_with_sources(
        bytes: &[u8], maximum_block: usize, gate: Arc<CancelGate>,
        sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, ResolvedDgInputError> {
        poll(&gate)?;
        if bytes.len() > MAX_RESOLVED_PERFORMANCE_BYTES || !(1..=65_536).contains(&maximum_block) {
            return Err(bad(0, "resolved source exceeds its byte or callback budget"));
        }
        let text = core::str::from_utf8(bytes).map_err(|_| bad(0, "resolved source must be UTF-8"))?;
        let (schema, body) = text.split_once('\n').ok_or_else(|| bad(1, "missing source schema"))?;
        if schema != RESOLVED_DG_PERFORMANCE_SCHEMA { return Err(bad(1, "unsupported resolved source schema")); }
        let mut rows = Rows { lines: body.lines(), line: 1 };
        let mut row = rows.row("clock")?;
        let sample_rate_hz: u32 = row.number()?;
        let substeps: usize = row.number()?;
        row.end()?;
        if !(8_000..=3_072_000).contains(&sample_rate_hz) || !(1..=65_536).contains(&substeps) {
            return Err(bad(rows.line, "clock needs 8..3072 kHz and 1..65536 original NSF substeps"));
        }
        let mut row = rows.row("output")?;
        let samples: usize = row.number()?;
        let full_scale_pa = row.scalar()?;
        row.end()?;
        if !(1..=MAX_RESOLVED_SAMPLES).contains(&samples) || full_scale_pa <= 0.0 {
            return Err(bad(rows.line, "finite source window or physical PCM scale is out of bounds"));
        }

        let mut row = rows.row("gas")?;
        let (gas, model, mut gas_records) = match row.fields.clone().next() {
            Some("dry-air") => {
                row.word()?;
                let temperature = row.scalar()?;
                let pressure = row.scalar()?;
                let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure)
                    .map_err(|e| ResolvedDgInputError::Flow(NsfError::Gas(e.to_string())))?;
                (gas, None, vec![GasRecord::Legacy(gas)])
            }
            Some("thermochemical") => {
                let registry = sources.ok_or_else(|| bad(row.line, "resolved mixture needs original thermochemical source packs"))?;
                let (model, point) = read_frozen_thermochemical_model(&mut row.fields, registry)
                    .map_err(ResolvedDgInputError::Gas)?;
                let resolved = model.evaluate(&point).map_err(|e|
                    ResolvedDgInputError::Gas(GasInputError::Material(e)))?;
                (*resolved.state(), Some((model, point)), vec![GasRecord::Thermochemical(resolved)])
            }
            _ => return Err(bad(row.line, "resolved gas requires dry-air T P or a complete thermochemical source; snapshots cannot evolve")),
        };
        row.end()?;
        let mut row = rows.row("grid")?;
        let dims: [usize; 3] = [row.number()?, row.number()?, row.number()?];
        let dx_m = row.scalar()?;
        let periodic = [row.flag()?, row.flag()?, row.flag()?];
        row.end()?;
        let cells = dims.iter().try_fold(1usize, |n, &d| n.checked_mul(d))
            .ok_or_else(|| bad(rows.line, "grid cell count overflow"))?;
        if dims.contains(&0) || cells > MAX_RESOLVED_CELLS || dx_m <= 0.0 {
            return Err(bad(rows.line, "positive physical grid exceeds the resolved cell budget"));
        }
        let mut row = rows.row("fluid_limits")?;
        let maximum_cells: usize = row.number()?;
        let temperature_window_k = [row.scalar()?, row.scalar()?];
        let maximum_mach = row.scalar()?;
        let maximum_temperature_ratio = row.scalar()?;
        let transport_points: usize = row.number()?;
        let workers: usize = row.number()?;
        row.end()?;
        if cells > maximum_cells || maximum_cells > MAX_RESOLVED_CELLS
            || !(2..=16_384).contains(&transport_points) || !(1..=256).contains(&workers)
        {
            return Err(bad(rows.line, "fluid dimensions, transport quadrature or workers exceed their explicit budget"));
        }
        let medium: NsfGasMedium = if let Some((model, point)) = model {
            let mut row = rows.row("energy_solve")?;
            let energy_tolerance_j_kg = row.scalar()?;
            let maximum_iterations = row.number()?;
            row.end()?;
            ThermochemicalNsfGas::new(model, point, GasEnergySolve {
                minimum_temperature_k: temperature_window_k[0],
                maximum_temperature_k: temperature_window_k[1],
                energy_tolerance_j_kg, maximum_iterations,
            }).map_err(ResolvedDgInputError::Flow)?.into()
        } else { GasSpec::dry_air_ussa1976().into() };
        let mut row = rows.row("bulk_viscosity")?;
        let bulk_to_shear_ratio = row.scalar()?;
        let bulk_ratio_source = row.word()?.to_string();
        row.end()?;
        let mut row = rows.row("wall")?;
        let wall_thermal = match row.word()? {
            "adiabatic" => WallThermal::Adiabatic,
            "isothermal" => WallThermal::Isothermal(row.scalar()?),
            _ => return Err(bad(row.line, "wall requires adiabatic or isothermal T")),
        };
        row.end()?;
        let sgs = if rows.next_is("sgs") {
            let mut row = rows.row("sgs")?;
            if row.word()? != "exploration" {
                return Err(bad(row.line, "SGS is an explicit exploration closure, not an acceptance render"));
            }
            let kind = row.word()?;
            let coefficient = row.scalar()?;
            let model = match kind {
                "smagorinsky" => SgsModel::Smagorinsky { cs: coefficient },
                "wale" => SgsModel::Wale { cw: coefficient },
                "vreman" => SgsModel::Vreman { cv: coefficient },
                _ => return Err(bad(row.line, "SGS requires smagorinsky, wale or vreman")),
            };
            let config = NsfSgsConfig {
                closure: SgsConfig {
                    model, filter_width_lu: row.scalar()?,
                    maximum_eddy_to_molecular: row.scalar()?, usage: SgsUse::Exploration,
                },
                turbulent_prandtl: row.scalar()?,
            };
            row.end()?;
            config.validate().map_err(ResolvedDgInputError::Flow)?;
            Some(config)
        } else { None };
        let mut row = rows.row("characteristic")?;
        let axis: usize = row.number()?;
        let high = match row.word()? { "low" => false, "high" => true,
            _ => return Err(bad(row.line, "coupled face requires low or high")) };
        let maximum_relative_perturbation = row.scalar()?;
        row.end()?;
        if axis >= 3 || periodic[axis] { return Err(bad(rows.line, "characteristic axis must be a nonperiodic physical axis")); }
        let coupled_patch = if rows.next_is("coupled_patch") {
            let mut row = rows.row("coupled_patch")?;
            let start0 = row.number()?;
            let end0 = row.number()?;
            let start1 = row.number()?;
            let end1 = row.number()?;
            let background = match row.word()? {
                "pressure" => NsfBoundaryDrive::PressureOffset(row.scalar()?),
                "pressure-temperature" => NsfBoundaryDrive::PressureTemperature {
                    pressure_offset_pa: row.scalar()?,
                    inflow_temperature_k: row.scalar()?,
                },
                "characteristic" => NsfBoundaryDrive::IncomingCharacteristic(row.scalar()?),
                _ => return Err(bad(row.line, "patch remainder requires pressure PA, pressure-temperature PA K, or characteristic M_S")),
            };
            if !background.is_finite() {
                return Err(bad(row.line, "patch controls require finite pressure and positive inflow temperature"));
            }
            row.end()?;
            Some((NsfCharacteristicPatch { start: [start0, start1], end: [end0, end1] }, background))
        } else { None };
        let reservoir_aperture = if rows.next_is("reservoir_aperture") {
            let mut row = rows.row("reservoir_aperture")?;
            let axis: usize = row.number()?;
            let high = match row.word()? {
                "low" => false, "high" => true,
                _ => return Err(bad(row.line, "reservoir face requires low or high")),
            };
            let state = NsfOpenState {
                density_kg_m3: row.scalar()?,
                velocity_m_s: [row.scalar()?, row.scalar()?, row.scalar()?],
                temperature_k: row.scalar()?,
            };
            let count: usize = row.number()?;
            row.end()?;
            if axis >= 3 || count == 0 || count > cells {
                return Err(bad(rows.line, "reservoir aperture requires a physical axis and bounded nonempty cell list"));
            }
            let mut selected = Vec::with_capacity(count);
            for n in 0..count {
                if n % 256 == 0 { poll(&gate)?; }
                let mut row = rows.row("reservoir_cell")?;
                let index = row.index(dims)?;
                row.end()?;
                selected.push([index % dims[0], (index / dims[0]) % dims[1], index / (dims[0] * dims[1])]);
            }
            Some(NsfReservoirAperture { axis, high, cells: selected, state })
        } else { None };
        let mut row = rows.row("initial")?;
        let initial = match row.word()? {
            "rest" => [gas.density, 0.0, 0.0, 0.0, gas.temperature],
            "uniform" => row.state()?,
            _ => return Err(bad(row.line, "initial requires rest or uniform RHO VX VY VZ T")),
        };
        row.end()?;
        let mut row = rows.row("solids")?;
        let solid_cells: usize = row.number()?;
        row.end()?;
        if solid_cells >= cells { return Err(bad(rows.line, "solid mask must leave actual fluid cells")); }
        let mut solid = vec![false; cells];
        for n in 0..solid_cells {
            if n % 256 == 0 { poll(&gate)?; }
            let mut row = rows.row("solid")?;
            let index = row.index(dims)?;
            if solid[index] { return Err(bad(row.line, "duplicate solid cell")); }
            solid[index] = true;
            row.end()?;
        }
        let mut row = rows.row("states")?;
        let initial_state_cells: usize = row.number()?;
        row.end()?;
        if initial_state_cells > cells - solid_cells { return Err(bad(rows.line, "state count exceeds actual fluid cells")); }
        let mut states = vec![None; cells];
        for n in 0..initial_state_cells {
            if n % 256 == 0 { poll(&gate)?; }
            let mut row = rows.row("state")?;
            let index = row.index(dims)?;
            if solid[index] || states[index].is_some() {
                return Err(bad(row.line, "initial state must name one unique fluid cell"));
            }
            states[index] = Some(row.state()?);
            row.end()?;
        }
        let mut row = rows.row("coupling")?;
        let velocity_unit = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt();
        let physical_dt = dx_m / velocity_unit * substeps as f64;
        let physical_intervals = source_intervals_for(physical_dt, sample_rate_hz, samples as u64)
            .map_err(ResolvedDgInputError::Render)?;
        if physical_intervals > MAX_RESOLVED_SAMPLES {
            return Err(bad(rows.line, "resolved observation exceeds the original physical-interval budget"));
        }
        let limits = NsfDgLimits {
            substeps, maximum_cells, maximum_steps: physical_intervals as u64,
            maximum_evaluations: row.number()?, maximum_cell_steps: row.number()?,
            incoming_bracket_m_s: [row.scalar()?, row.scalar()?],
            pressure_absolute_tolerance_pa: row.scalar()?, pressure_relative_tolerance: row.scalar()?,
            projection_absolute_tolerance_j: row.scalar()?, projection_relative_tolerance: row.scalar()?,
        };
        row.end()?;
        limits.validate_for_cells(cells).map_err(ResolvedDgInputError::Coupling)?;
        let (source, wall) = read_load(&mut rows, false, &gas, sources,
            &mut gas_records, &gate, &mut load_material)?;
        let upstream = if rows.next_is("upstream_dg") {
            Some(read_load(&mut rows, true, &gas, sources,
                &mut gas_records, &gate, &mut load_material)?)
        } else { None };
        let mut row = rows.row("observation")?;
        let observation = match row.word()? {
            "inlet-pressure" => ApertureObservation::Inlet,
            "baffled" => ApertureObservation::TubeBaffled(CircularOutletReceiver {
                position_m: [row.scalar()?, row.scalar()?, row.scalar()?],
                radial_rings: row.number()?, angular_points: row.number()?,
                maximum_frequency_hz: row.scalar()?,
            }),
            _ => return Err(bad(row.line, "observation requires inlet-pressure or baffled X Y Z RINGS ANGLES FMAX")),
        };
        row.end()?;
        let reservoir_observation = if rows.next_is("reservoir_observation") {
            let mut row = rows.row("reservoir_observation")?;
            if row.word()? != "baffled" {
                return Err(bad(row.line, "reservoir_observation requires baffled OX OY ANGLE_RAD"));
            }
            let input = ReservoirObservation { origin_m: [row.scalar()?, row.scalar()?], angle_rad: row.scalar()? };
            row.end()?;
            Some(input)
        } else { None };
        let observer_gas = if rows.next_is("observer_gas") {
            let mut row = rows.row("observer_gas")?;
            let exterior = source.exterior_gas_record()
                .ok_or_else(|| bad(row.line, "observer_gas requires the original DG exterior gas record"))?;
            let value = GasPropagation::read_record_for_gas(&mut row.fields, exterior)
                .map_err(|error| ResolvedDgInputError::Gas(GasInputError::Material(error)))?;
            row.end()?;
            Some(value)
        } else { None };
        let mut row = rows.row("boundary")?;
        // The legacy bare count denotes m/s, exactly as before. Explicit
        // pressure data is in Pa relative to the source face's gas reference.
        let boundary_kind = match row.fields.clone().next() {
            Some("pressure") => { row.word()?; 1 }
            Some("upstream-pressure") => { row.word()?; 2 }
            Some("reservoir") => { row.word()?; 3 }
            Some("pressure-temperature") => { row.word()?; 4 }
            Some("characteristic") => { row.word()?; 0 }
            _ => 0,
        };
        if boundary_kind == 3 && reservoir_aperture.is_none() {
            return Err(bad(row.line, "reservoir controls require the declared transverse aperture"));
        }
        if boundary_kind != 3 && upstream.is_some() != (boundary_kind == 2) {
            return Err(bad(row.line, "upstream-pressure requires its declared upstream DG load; face controls require an unloaded source face"));
        }
        let boundary_segments: usize = row.number()?;
        row.end()?;
        if boundary_segments == 0 || boundary_segments > 16_384 || boundary_segments > samples {
            return Err(bad(rows.line, "finite boundary history needs 1..16384 nonempty segments"));
        }
        let mut boundaries = Vec::with_capacity(boundary_segments);
        let mut previous_end = 0;
        for n in 0..boundary_segments {
            if n % 256 == 0 { poll(&gate)?; }
            let mut row = rows.row("hold")?;
            let end: usize = row.number()?;
            let value = row.scalar()?;
            let drive = match boundary_kind {
                4 => {
                    let temperature = row.scalar()?;
                    if !(temperature_window_k[0]..=temperature_window_k[1]).contains(&temperature) {
                        return Err(bad(row.line, "source inflow temperature lies outside the declared fluid window"));
                    }
                    NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
                        pressure_offset_pa: value, inflow_temperature_k: temperature,
                    })
                },
                3 => {
                    let temperature_k = row.scalar()?;
                    if value <= 0.0 || temperature_k <= 0.0 {
                        return Err(bad(row.line, "reservoir controls require positive absolute pressure and temperature"));
                    }
                    NsfDgDrive::Reservoir { pressure_pa: value, temperature_k }
                },
                2 => NsfDgDrive::UpstreamPressureRise(value),
                1 => NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(value)),
                _ => NsfDgDrive::Boundary(NsfBoundaryDrive::IncomingCharacteristic(value)),
            };
            row.end()?;
            if end <= previous_end || end > samples {
                return Err(bad(rows.line, "boundary endpoints must increase inside the declared sample window"));
            }
            boundaries.push((end, drive));
            previous_end = end;
        }
        if previous_end != samples { return Err(bad(rows.line, "boundary history must cover the complete declared window")); }
        rows.row("end")?.end()?;
        if rows.lines.next().is_some() { return Err(bad(rows.line + 1, "unexpected records after end")); }
        poll(&gate)?;

        // Control endpoints remain on the caller's physical observation clock
        // when the atmosphere changes the lattice clock. Each original source
        // interval holds the value at its start. Event timing error is bounded
        // by that interval and refines with the same physical discretization.
        let matching_clock = same_clock(physical_dt, 1.0 / f64::from(sample_rate_hz));
        let mut boundary = Vec::with_capacity(physical_intervals);
        let mut segment = 0;
        for index in 0..physical_intervals {
            if index % 256 == 0 { poll(&gate)?; }
            let start = if matching_clock { index as f64 }
                else { index as f64 * (physical_dt * f64::from(sample_rate_hz)) };
            while segment < boundaries.len() && start >= boundaries[segment].0 as f64 {
                segment += 1;
            }
            let value = boundaries.get(segment).ok_or_else(|| bad(0,
                "physical source interval starts beyond the declared boundary window"))?.1;
            boundary.push(value);
        }
        let coordinate = |x: [usize; 3]| (x[2] * dims[1] + x[1]) * dims[0] + x[0];
        let zone = NsfGrid3::<D3q27>::new(NsfSpec {
            dims, periodic, gas: medium, temperature_k: gas.temperature,
            pressure_pa: gas.pressure, dx_m, temperature_window_k,
            maximum_mach, maximum_temperature_ratio, bulk_to_shear_ratio, bulk_ratio_source,
            maximum_cells, wall: wall_thermal, workers, transport_points,
        }, |x| solid[coordinate(x)], |x| {
            let state = states[coordinate(x)].unwrap_or(initial);
            (state[0] / gas.density,
                [state[1] / velocity_unit, state[2] / velocity_unit, state[3] / velocity_unit], state[4])
        }).map_err(ResolvedDgInputError::Flow)?;
        let dt = zone.law().dt_s * substeps as f64;
        if !same_clock(dt, physical_dt) {
            return Err(bad(0, "original NSF clock disagrees with its source-derived physical units"));
        }
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: gas.pressure, reference_temperature_k: gas.temperature,
            reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
        };
        let mut zone = zone.with_characteristic_axis(NsfCharacteristicAxis {
            axis, faces: [face; 2], maximum_relative_perturbation,
        }).map_err(ResolvedDgInputError::Flow)?;
        if let Some((patch, background)) = coupled_patch {
            let mut drives = zone.boundary_drives().expect("constructed characteristic axis");
            drives[usize::from(high)] = background;
            zone.set_boundary_drives(drives).map_err(ResolvedDgInputError::Flow)?;
            zone = zone.with_characteristic_patch(high, patch,
                NsfBoundaryDrive::IncomingCharacteristic(0.0)).map_err(ResolvedDgInputError::Flow)?;
        }
        if let Some(aperture) = reservoir_aperture {
            zone = zone.with_reservoir_aperture(aperture).map_err(ResolvedDgInputError::Flow)?;
        }
        if let Some(config) = sgs {
            zone = zone.with_sgs(config).map_err(ResolvedDgInputError::Flow)?;
        }
        poll(&gate)?;
        let load = source.prepare_with_material_wall_layout(&gas, dt, &gate, &wall)?;
        let mut system = if coupled_patch.is_some() {
            NsfDg::new_on_patch(zone, load, high, limits, &gate)
        } else {
            NsfDg::new(zone, load, high, limits, &gate)
        }.map_err(ResolvedDgInputError::Coupling)?;
        let upstream_wall = if let Some((source, wall)) = upstream {
            let load = source.prepare_with_material_wall_layout(&gas, dt, &gate, &wall)?;
            system = system.with_upstream_load(load, &gate).map_err(ResolvedDgInputError::Coupling)?;
            Some(wall)
        } else { None };
        let mut renderer = NsfDgPerformance::new_resampled_with_drives(system, observation, sample_rate_hz,
            boundary, samples as u64, maximum_block, Arc::clone(&gate)).map_err(ResolvedDgInputError::Render)?;
        if let Some(input) = reservoir_observation {
            renderer = renderer.with_reservoir_observation(input).map_err(ResolvedDgInputError::Render)?;
        }
        if let Some(propagation) = observer_gas {
            renderer = renderer.with_observer_gas(&propagation).map_err(ResolvedDgInputError::Render)?;
        }
        poll(&gate)?;
        Ok(Self {
            info: ResolvedDgPerformanceInfo {
                input_hash: hash_domain("org.frankensim.fs-couple.resolved-dg-performance.v1", bytes),
                sample_rate_hz, samples: samples as u64, full_scale_pa, dimensions: dims,
                solid_cells, initial_state_cells, boundary_segments,
                physical_intervals: physical_intervals as u64,
                sgs,
            }, renderer, gas_records, wall, upstream_wall,
        })
    }
    pub const fn info(&self) -> ResolvedDgPerformanceInfo { self.info }
    pub const fn renderer(&self) -> &NsfDgPerformance<D3q27> { &self.renderer }
    pub fn gas_records(&self) -> &[GasRecord] { &self.gas_records }
    pub const fn wall_layout(&self) -> &DgMaterialWallLayout { &self.wall }
    pub fn upstream_wall_layout(&self) -> Option<&DgMaterialWallLayout> { self.upstream_wall.as_ref() }
    pub fn into_renderer(self) -> NsfDgPerformance<D3q27> { self.renderer }
}

// Both sides consume exactly the same geometry, gas, wall and molecular source
// readers. Prefixes only identify the physical participant; they add no law.
#[allow(clippy::too_many_arguments)]
fn read_load(
    rows: &mut Rows<'_>, upstream: bool, gas: &GasState,
    sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>, gate: &CancelGate,
    load_material: &mut impl FnMut(&str)->Result<NormalizedMaterialCardPack, String>,
) -> Result<(DgSourceInput, DgMaterialWallLayout), ResolvedDgInputError> {
    let (kind, mut row) = if upstream {
        let mut row = rows.row("upstream_dg")?;
        let kind = match row.word()? {
            "dg_tube" => "dg_tube", "dg_bore" => "dg_bore",
            _ => return Err(bad(row.line, "upstream_dg requires the original dg_tube or dg_bore grammar")),
        };
        (kind, row)
    } else {
        let kind = if rows.next_is("dg_tube") { "dg_tube" } else { "dg_bore" };
        (kind, rows.row(kind)?)
    };
    let source = DgSourceInput::read(kind, row.line, &mut row.fields,
        &mut || rows.next(), gas, sources, records, gate)?;
    row.end()?;
    let molecular_key = if upstream { "upstream_molecular" } else { "dg_molecular" };
    let source = if rows.next_is(molecular_key) {
        let mut row = rows.row(molecular_key)?;
        let law = fs_material::relaxation::input::read_relaxation(&mut row.fields)
            .map_err(|error| ResolvedDgInputError::Gas(GasInputError::Material(error)))?;
        row.end()?;
        source.with_molecular_relaxation(law, &records[0])?
    } else { source };
    let wall_key = if upstream { "upstream_wall" } else { "dg_wall" };
    let wall = if rows.next_is(wall_key) {
        let mut row = rows.row(wall_key)?;
        let declaration = DgMaterialWallLayoutInput::read(row.line, &mut row.fields,
            &mut || rows.next())?;
        row.end()?;
        declaration.resolve(gate, &mut |path, line| {
            let pack = load_material(path).map_err(|source|
                ResolvedDgInputError::Material { line, source })?;
            Ok::<_, ResolvedDgInputError>(Arc::new(pack))
        })?
    } else { DgMaterialWallLayout::default() };
    Ok((source, wall))
}

#[cfg(test)]
mod tests;
