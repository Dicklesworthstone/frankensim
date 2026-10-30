//! Performance data for the existing conservative stationary gas generator.
//! Geometry comes from the SAME axial segment records used by acoustics.

use super::{
    GasFieldResolution, GasRecord, Reader, ReedPerformanceError, SourceColumns,
    ThermochemicalSources, bad, column_error, read_resolution,
};
use fs_duct::Segment;
use crate::acoustic_realize::AcousticRealizeError;
use fs_matdb::QueryPoint;
use fs_material::gas_field::steady::{
    AxialGasCell, GasInjection, MAX_STEADY_GAS_CELLS, MAX_STEADY_GAS_INJECTIONS, SteadyAxialGas,
    SteadyGasSolution, SteadyGasSolve,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

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
    pub fn problem(&self) -> &SteadyAxialGas {
        &self.problem
    }
    /// Explicit finite work, caloric bracket, balance and low-Mach admission limits.
    pub fn controls(&self) -> SteadyGasSolve {
        self.controls
    }
    /// Published mass/species/enthalpy balances and the SAME acoustic gas field.
    pub fn solution(&self) -> &SteadyGasSolution {
        &self.solution
    }
}

/// gas_steady MODE RESOLUTION CELLS INJECTIONS mole SPECIES_COUNT ID...
/// gas_inlet MASS_FLOW_KG_S T_K FRACTION...
/// gas_cell END_M WALL_CONDUCTANCE_W_K WALL_T_K INITIAL_T_K (per cell)
/// gas_injection CELL MASS_FLOW_KG_S T_K FRACTION... (zero-based cell)
/// gas_solve MIN_T MAX_T H_TOL MAX_ITER DT MAX_RELAX MAX_SUBSTEPS P_TOL E_TOL MAX_MACH
///
/// The wall conductance is an explicit total W/K, never an inferred transport
/// coefficient. This one-dimensional generator has no side-branch mass balance;
/// a tone-hole network must supply its own prescribed field instead.
pub(in super::super) fn read(
    reader: &mut Reader<'_>,
    sources: Option<&ThermochemicalSources<'_>>,
    pressure_pa: f64,
    segments: &[Segment],
    records: &mut Vec<GasRecord>,
) -> Result<(SteadyGasPerformance, GasFieldResolution), ReedPerformanceError> {
    let mut row = reader.row("gas_steady")?;
    let line = row.line;
    let resolution = read_resolution(&mut row)?;
    let count = row.count(MAX_STEADY_GAS_CELLS)?;
    let injection_count = row.count(MAX_STEADY_GAS_INJECTIONS)?;
    if count == 0 {
        return Err(bad(line, "gas_steady requires at least one axial cell"));
    }
    let columns = SourceColumns::read(&mut row.words, sources, count + injection_count + 1)
        .map_err(|error| column_error(row.line, error))?;
    row.finish()?;
    let length = segments.iter().try_fold(0.0, |sum, segment| match segment {
        Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => Ok(sum + length),
        Segment::Profile { profile, .. } => Ok(sum + profile.length()),
        Segment::ToneHole { .. } => Err(bad(line,
            "gas_steady needs an unbranched bore; side-branch mean flow is not specified")),
    })?;
    if !length.is_finite() || length <= 0.0 {
        return Err(bad(
            line,
            "gas_steady requires finite positive complete bore length",
        ));
    }
    let mut row = reader.row("gas_inlet")?;
    let inlet_mass_flow_kg_s = row.positive()?;
    let inlet = columns.read_state(&mut row.words)
        .map_err(|error| column_error(row.line, error))?;
    row.finish()?;
    let mut cells = Vec::with_capacity(count);
    let mut start = 0.0;
    let mut geometry_steps = 0;
    for _ in 0..count {
        let mut row = reader.row("gas_cell")?;
        let end_m = row.positive()?;
        if end_m <= start || end_m > length {
            return Err(bad(
                row.line,
                "gas cell ends must increase within the authored bore",
            ));
        }
        let (volume_m3, wall_area_m2) = geometry(segments, start, end_m, row.line, &mut geometry_steps)?;
        let cell = AxialGasCell {
            end_m,
            volume_m3,
            wall_area_m2,
            wall_conductance_w_k: row.scalar()?,
            wall_temperature_k: row.positive()?,
            initial_temperature_k: row.positive()?,
        };
        row.finish()?;
        cells.push(cell);
        start = end_m;
    }
    if start != length {
        return Err(bad(
            reader.line,
            "gas cells must cover the exact complete authored bore length",
        ));
    }
    let mut injections = Vec::with_capacity(injection_count);
    for _ in 0..injection_count {
        let mut row = reader.row("gas_injection")?;
        let cell = row.count(count - 1)?;
        let mass_flow_kg_s = row.positive()?;
        let state = columns.read_state(&mut row.words)
            .map_err(|error| column_error(row.line, error))?;
        row.finish()?;
        injections.push(GasInjection {
            cell,
            mass_flow_kg_s,
            state,
        });
    }
    let mut row = reader.row("gas_solve")?;
    let controls = SteadyGasSolve {
        caloric: GasEnergySolve {
            minimum_temperature_k: row.positive()?,
            maximum_temperature_k: row.positive()?,
            energy_tolerance_j_kg: row.positive()?,
            maximum_iterations: row.count(128)?,
        },
        maximum_dt_s: row.positive()?,
        maximum_relaxations: row.count(1_000_000)?,
        maximum_substeps: row.count(1_000_000)?,
        steady_power_tolerance_w: row.positive()?,
        energy_tolerance_j: row.scalar()?,
        maximum_mean_mach: row.positive()?,
    };
    row.finish()?;
    let (model, phase) = columns.into_parts();
    let problem = SteadyAxialGas {
        source: model,
        source_conditions: QueryPoint::new(),
        static_pressure_pa: pressure_pa,
        phase,
        inlet,
        inlet_mass_flow_kg_s,
        cells,
        injections,
    };
    let solution = problem
        .solve(controls, || false)
        .map_err(|source| ReedPerformanceError::SteadyGas { line, source })?;
    // Retain the generated exact stations on the original observation surface.
    // The admitted problem separately preserves the supplied injection states.
    for &position in solution.field().stations_m() {
        let state = solution
            .field()
            .state_at(position)
            .map_err(|source| ReedPerformanceError::Gas { line, source })?;
        records.push(GasRecord::Thermochemical(state));
    }
    Ok((
        SteadyGasPerformance {
            problem,
            controls,
            solution,
        },
        resolution,
    ))
}

/// Exact circular frusta and the profile owner's refinable Gauss measures over
/// each cell/segment intersection. The complete thermal-cell binding is bounded
/// by the original profile work cap. No bore geometry is copied into the file.
fn geometry(
    segments: &[Segment],
    start: f64,
    end: f64,
    line: usize,
    geometry_steps: &mut usize,
) -> Result<(f64, f64), ReedPerformanceError> {
    let (mut offset, mut volume, mut area) = (0.0_f64, 0.0_f64, 0.0_f64);
    for segment in segments {
        if matches!(segment, Segment::ToneHole { .. }) {
            return Err(bad(line, "gas cell geometry requires an unbranched bore"));
        }
        let length = segment.axial_length();
        let lo = start.max(offset);
        let hi = end.min(offset + length);
        if hi > lo {
            if let Segment::Profile { profile, steps } = segment {
                // Preserve exact local endpoints under accumulated segment
                // offsets: subtracting them can round outside the support.
                let local_start = if lo == offset { 0.0 } else { lo - offset };
                let local_end = if hi == offset + length { length } else { hi - offset };
                let mut exhausted = false;
                let measures = profile.geometry_controlled(local_start, local_end, *steps, &mut || {
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
                let measures = measures
                    .map_err(|source| ReedPerformanceError::Model(AcousticRealizeError::Duct(source)))?;
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
        return Err(bad(
            line,
            "gas cell volume and lateral wall area must be representable",
        ));
    }
    Ok((volume, area))
}

#[cfg(test)]
mod tests;
