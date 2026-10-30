//! One authored stationary gas builder for characteristic-line and DG sources.
//! Geometry, species transport and caloric inversion remain with their owners.

use fs_duct::{DuctError, Segment};
use fs_matdb::QueryPoint;
use fs_material::gas_field::input::PrimitiveGasColumns;
use fs_material::gas_field::steady::{
    AxialGasCell, GasInjection, MAX_STEADY_GAS_CELLS, MAX_STEADY_GAS_INJECTIONS,
    SteadyAxialGas, SteadyGasError, SteadyGasSolution, SteadyGasSolve,
};
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use std::str::SplitAsciiWhitespace;

/// Geometry-bound stationary gas transport through a bore and open side paths.
pub mod branched;

/// The actual thermal/species problem and accepted solution behind a performance.
/// These are owner results, not a separate reconstruction from acoustic output.
#[derive(Clone, Debug)]
pub struct SteadyGasPerformance {
    problem: SteadyAxialGas,
    controls: SteadyGasSolve,
    solution: SteadyGasSolution,
}
impl SteadyGasPerformance {
    /// Source-bound inlet, injections and wall cells, with geometry from the bore.
    pub fn problem(&self) -> &SteadyAxialGas { &self.problem }
    /// Explicit finite work, caloric bracket, balance and low-Mach admission limits.
    pub fn controls(&self) -> SteadyGasSolve { self.controls }
    /// Published mass/species/enthalpy balances and the SAME acoustic gas field.
    pub fn solution(&self) -> &SteadyGasSolution { &self.solution }
}

#[derive(Debug)]
pub(crate) enum SteadyGasInputError {
    Input { line: usize, what: &'static str },
    Gas { line: usize, source: fs_material::MaterialError },
    Geometry(DuctError),
    SteadyGas { line: usize, source: SteadyGasError },
}
pub(crate) enum SteadyGasReadError<E> {
    Input(SteadyGasInputError),
    Reader(E),
}
impl<E> From<SteadyGasInputError> for SteadyGasReadError<E> {
    fn from(error: SteadyGasInputError) -> Self { Self::Input(error) }
}
impl<E> SteadyGasReadError<E> {
    pub(crate) fn into_error(self, source: impl FnOnce(SteadyGasInputError) -> E) -> E {
        match self { Self::Input(error) => source(error), Self::Reader(error) => error }
    }
}
fn bad(line: usize, what: &'static str) -> SteadyGasInputError {
    SteadyGasInputError::Input { line, what }
}
fn column_error(line: usize, error: GasInputError) -> SteadyGasInputError {
    match error {
        GasInputError::Input(what) => bad(line, what),
        GasInputError::Material(source) => SteadyGasInputError::Gas { line, source },
    }
}
fn poll(line: usize, cancelled: &mut impl FnMut() -> bool) -> Result<(), SteadyGasInputError> {
    if cancelled() {
        Err(SteadyGasInputError::SteadyGas {
            line,
            source: SteadyGasError::Thermal(fs_thermal::enthalpy::DepositionError::Cancelled),
        })
    } else { Ok(()) }
}

struct Row<'a> { line: usize, words: SplitAsciiWhitespace<'a> }
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str, SteadyGasInputError> {
        self.words.next().ok_or_else(|| bad(self.line, "missing gas_steady field"))
    }
    fn scalar(&mut self) -> Result<f64, SteadyGasInputError> {
        let value: f64 = self.word()?.parse()
            .map_err(|_| bad(self.line, "malformed gas_steady scalar"))?;
        if !value.is_finite() { return Err(bad(self.line, "physical scalars must be finite")); }
        Ok(value)
    }
    fn positive(&mut self) -> Result<f64, SteadyGasInputError> {
        let value = self.scalar()?;
        if value <= 0.0 { return Err(bad(self.line, "physical dimension/scale must be positive")); }
        Ok(value)
    }
    fn count(&mut self, maximum: usize) -> Result<usize, SteadyGasInputError> {
        let value: usize = self.word()?.parse()
            .map_err(|_| bad(self.line, "malformed gas_steady count"))?;
        if value > maximum { return Err(bad(self.line, "count exceeds the input budget")); }
        Ok(value)
    }
    fn finish(&mut self) -> Result<(), SteadyGasInputError> {
        if self.words.next().is_some() { return Err(bad(self.line, "unexpected extra field")); }
        Ok(())
    }
}
fn row<'a, E>(
    next: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
    key: &str,
) -> Result<Row<'a>, SteadyGasReadError<E>> {
    let (line, words) = next().map_err(SteadyGasReadError::Reader)?;
    let mut row = Row { line, words };
    if row.word()? != key { return Err(bad(line, "unexpected gas_steady record").into()); }
    Ok(row)
}

/// Specific, bounded gas_steady header; acoustic discretization stays with its caller.
pub(crate) struct SteadyGasInput {
    line: usize,
    count: usize,
    injection_count: usize,
    columns: PrimitiveGasColumns,
}
impl SteadyGasInput {
    /// CELLS INJECTIONS mole SPECIES_COUNT ID... after the caller's acoustic axis.
    pub(crate) fn header(
        line: usize, fields: &mut SplitAsciiWhitespace<'_>,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, SteadyGasInputError> {
        let mut row = Row { line, words: fields.clone() };
        let count = row.count(MAX_STEADY_GAS_CELLS)?;
        let injection_count = row.count(MAX_STEADY_GAS_INJECTIONS)?;
        if count == 0 { return Err(bad(line, "gas_steady requires at least one axial cell")); }
        let columns = PrimitiveGasColumns::read(&mut row.words, sources, count + injection_count + 1)
            .map_err(|error| column_error(line, error))?;
        row.finish()?;
        *fields = row.words;
        Ok(Self { line, count, injection_count, columns })
    }

    /// Consume the existing inlet/cell/injection/solve rows and publish the
    /// original stationary solution only after complete source and work admission.
    pub(crate) fn read<'a, E>(
        self,
        next: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
        pressure_pa: f64,
        segments: &[Segment],
        records: &mut Vec<GasRecord>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SteadyGasPerformance, SteadyGasReadError<E>> {
        let Self { line, count, injection_count, columns } = self;
        poll(line, cancelled)?;
        let length = segments.iter().try_fold(0.0, |sum, segment| match segment {
            Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => Ok(sum + length),
            Segment::Profile { profile, .. } => Ok(sum + profile.length()),
            Segment::ToneHole { .. } => Err(bad(line,
                "gas_steady needs an unbranched bore; side-branch mean flow is not specified")),
        })?;
        if !length.is_finite() || length <= 0.0 {
            return Err(bad(line, "gas_steady requires finite positive complete bore length").into());
        }
        let mut inlet_row = row(next, "gas_inlet")?;
        let inlet_mass_flow_kg_s = inlet_row.positive()?;
        let inlet = columns.read_state(&mut inlet_row.words)
            .map_err(|error| column_error(inlet_row.line, error))?;
        inlet_row.finish()?;
        let mut cells = Vec::with_capacity(count);
        let mut start = 0.0;
        let mut geometry_steps = 0;
        let mut last_line = line;
        for _ in 0..count {
            poll(line, cancelled)?;
            let mut cell = row(next, "gas_cell")?;
            last_line = cell.line;
            let end_m = cell.positive()?;
            if end_m <= start || end_m > length {
                return Err(bad(cell.line, "gas cell ends must increase within the authored bore").into());
            }
            let (volume_m3, wall_area_m2) = geometry(
                segments, start, end_m, cell.line, &mut geometry_steps, cancelled,
            )?;
            cells.push(AxialGasCell {
                end_m, volume_m3, wall_area_m2,
                wall_conductance_w_k: cell.scalar()?,
                wall_temperature_k: cell.positive()?,
                initial_temperature_k: cell.positive()?,
            });
            cell.finish()?;
            start = end_m;
        }
        if start != length {
            return Err(bad(last_line, "gas cells must cover the exact complete authored bore length").into());
        }
        let mut injections = Vec::with_capacity(injection_count);
        for _ in 0..injection_count {
            poll(line, cancelled)?;
            let mut injection = row(next, "gas_injection")?;
            let cell = injection.count(count - 1)?;
            let mass_flow_kg_s = injection.positive()?;
            let state = columns.read_state(&mut injection.words)
                .map_err(|error| column_error(injection.line, error))?;
            injection.finish()?;
            injections.push(GasInjection { cell, mass_flow_kg_s, state });
        }
        let mut solve = row(next, "gas_solve")?;
        let controls = SteadyGasSolve {
            caloric: GasEnergySolve {
                minimum_temperature_k: solve.positive()?,
                maximum_temperature_k: solve.positive()?,
                energy_tolerance_j_kg: solve.positive()?,
                maximum_iterations: solve.count(128)?,
            },
            maximum_dt_s: solve.positive()?,
            maximum_relaxations: solve.count(1_000_000)?,
            maximum_substeps: solve.count(1_000_000)?,
            steady_power_tolerance_w: solve.positive()?,
            energy_tolerance_j: solve.scalar()?,
            maximum_mean_mach: solve.positive()?,
        };
        solve.finish()?;
        let (model, phase) = columns.into_parts();
        let problem = SteadyAxialGas {
            source: model, source_conditions: QueryPoint::new(), static_pressure_pa: pressure_pa,
            phase, inlet, inlet_mass_flow_kg_s, cells, injections,
        };
        let solution = problem.solve(controls, &mut *cancelled)
            .map_err(|source| SteadyGasInputError::SteadyGas { line, source })?;
        for &position in solution.field().stations_m() {
            poll(line, cancelled)?;
            let state = solution.field().state_at(position)
                .map_err(|source| SteadyGasInputError::Gas { line, source })?;
            records.push(GasRecord::Thermochemical(state));
        }
        Ok(SteadyGasPerformance { problem, controls, solution })
    }
}

/// Exact frusta and the original profile Gauss measures, sharing one panel cap.
pub(crate) fn geometry(
    segments: &[Segment], start: f64, end: f64, line: usize,
    geometry_steps: &mut usize, cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), SteadyGasInputError> {
    let (mut offset, mut volume, mut area) = (0.0_f64, 0.0_f64, 0.0_f64);
    for segment in segments {
        poll(line, cancelled)?;
        if matches!(segment, Segment::ToneHole { .. }) {
            return Err(bad(line, "gas cell geometry requires an unbranched bore"));
        }
        let length = segment.axial_length();
        let lo = start.max(offset);
        let hi = end.min(offset + length);
        if hi > lo {
            if let Segment::Profile { profile, steps } = segment {
                let local_start = if lo == offset { 0.0 } else { lo - offset };
                let local_end = if hi == offset + length { length } else { hi - offset };
                let mut exhausted = false;
                let measures = profile.geometry_controlled(local_start, local_end, *steps, &mut || {
                    if cancelled() { return true; }
                    if *geometry_steps >= fs_duct::profile::MAX_RICCATI_STEPS {
                        exhausted = true;
                        true
                    } else {
                        *geometry_steps += 1;
                        false
                    }
                });
                if exhausted {
                    return Err(bad(line, "thermal profile geometry exceeds MAX_RICCATI_STEPS"));
                }
                let measures = measures.map_err(SteadyGasInputError::Geometry)?;
                volume += measures.volume_m3;
                area += measures.lateral_area_m2;
            } else {
                let inlet = segment.inlet_radius();
                let outlet = segment.outlet_radius();
                let r0 = inlet + (outlet - inlet) * ((lo - offset) / length);
                let r1 = inlet + (outlet - inlet) * ((hi - offset) / length);
                let dx = hi - lo;
                volume += core::f64::consts::PI * dx * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
                area += core::f64::consts::PI * (r0 + r1) * dx.hypot(r1 - r0);
            }
        }
        offset += length;
    }
    if !volume.is_finite() || volume <= 0.0 || !area.is_finite() || area <= 0.0 {
        return Err(bad(line, "gas cell volume and lateral wall area must be representable"));
    }
    Ok((volume, area))
}
