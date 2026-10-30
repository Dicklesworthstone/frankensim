//! One stationary thermal/species graph for a main bore and its open side paths.
//!
//! Authored positive outlet mass rates determine the decreasing main-bore flow.
//! Every side junction is an exact thermal-cell boundary, and the DG reader
//! additionally checks its original acoustic face. Volumes and wall areas come
//! from the actual bore/chimney geometry, never from an acoustic fit. The shared
//! material owner solves all cells together; this is a frozen played field,
//! not transient gas transport or a prediction of mean pressure drop.

use super::*;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_field::steady::network::{
    GasNetworkCell, GasOutlet, SteadyGasNetwork, SteadyGasNetworkSolution,
};
use fs_thermal::enthalpy::{CellGeometry, CellId};
use fs_thermal::enthalpy::advection::MassFlow;

/// Original open side geometry. A flush opening has no finite thermal volume.
pub(crate) struct GasSideGeometry {
    pub port: u32,
    pub cylinder: Option<(f64, f64)>,
}

/// The source-backed field in one original finite side branch.
#[derive(Clone, Debug)]
pub struct SteadyBranchField {
    port: u32,
    junction_m: f64,
    field: AxialGasField,
}
impl SteadyBranchField {
    pub fn port(&self) -> u32 { self.port }
    pub fn junction_m(&self) -> f64 { self.junction_m }
    pub fn field(&self) -> &AxialGasField { &self.field }
}

/// One accepted mass/species/enthalpy solution and its regional acoustic fields.
#[derive(Clone, Debug)]
pub struct BranchedGasPerformance {
    problem: SteadyGasNetwork,
    controls: SteadyGasSolve,
    solution: SteadyGasNetworkSolution,
    main: AxialGasField,
    branches: Vec<SteadyBranchField>,
    junctions: Vec<(u32, f64)>,
}
impl BranchedGasPerformance {
    pub fn problem(&self) -> &SteadyGasNetwork { &self.problem }
    pub fn controls(&self) -> SteadyGasSolve { self.controls }
    pub fn solution(&self) -> &SteadyGasNetworkSolution { &self.solution }
    pub fn main_field(&self) -> &AxialGasField { &self.main }
    pub fn branches(&self) -> &[SteadyBranchField] { &self.branches }
    /// Every open side path, including flush openings, in original port order.
    pub fn junctions(&self) -> &[(u32, f64)] { &self.junctions }
}

struct Region {
    port: u32,
    junction_m: f64,
    mass_flow_kg_s: f64,
    cells: Vec<AxialGasCell>,
    first: usize,
    donor: usize,
}

fn read_cells<'a, E>(
    next: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
    segments: &[Segment], count: usize, line: usize, geometry_steps: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<AxialGasCell>, SteadyGasReadError<E>> {
    let length: f64 = segments.iter().map(Segment::axial_length).sum();
    let mut cells = Vec::with_capacity(count);
    let mut start = 0.0;
    for _ in 0..count {
        poll(line, cancelled)?;
        let mut cell = row(next, "gas_cell")?;
        let end_m = cell.positive()?;
        if end_m <= start || end_m > length {
            return Err(bad(cell.line, "gas region cell ends must increase inside the original geometry").into());
        }
        let (volume_m3, wall_area_m2) = geometry(
            segments, start, end_m, cell.line, geometry_steps, cancelled,
        )?;
        cells.push(AxialGasCell { end_m, volume_m3, wall_area_m2,
            wall_conductance_w_k: cell.scalar()?, wall_temperature_k: cell.positive()?,
            initial_temperature_k: cell.positive()? });
        cell.finish()?;
        start = end_m;
    }
    if start != length {
        return Err(bad(line, "gas cells must cover the exact complete regional geometry").into());
    }
    Ok(cells)
}

fn append_cells(cells: &mut Vec<GasNetworkCell>, region: &[AxialGasCell]) {
    let mut start = 0.0;
    for cell in region {
        cells.push(GasNetworkCell {
            geometry: CellGeometry { start_m: [start, 0.0, 0.0],
                end_m: [cell.end_m, 0.0, 0.0], volume_m3: cell.volume_m3,
                surface_area_m2: cell.wall_area_m2 },
            flow_area_m2: cell.volume_m3 / (cell.end_m - start),
            wall_conductance_w_k: cell.wall_conductance_w_k,
            wall_temperature_k: cell.wall_temperature_k,
            initial_temperature_k: cell.initial_temperature_k,
        });
        start = cell.end_m;
    }
}

fn field(
    problem: &SteadyGasNetwork, solution: &SteadyGasNetworkSolution,
    cells: &[AxialGasCell], first: usize, inlet: PrimitiveGasState,
    line: usize, cancelled: &mut impl FnMut() -> bool,
) -> Result<AxialGasField, SteadyGasInputError> {
    let mut stations = vec![0.0];
    let mut states = vec![inlet];
    for (offset, cell) in cells.iter().enumerate() {
        poll(line, cancelled)?;
        stations.push(cell.end_m);
        states.push(PrimitiveGasState {
            temperature_k: solution.states()[first + offset].state().temperature,
            composition: solution.caloric_laws()[first + offset].mole_composition().clone(),
        });
    }
    AxialGasField::new(stations, problem.static_pressure_pa, states,
        problem.source.clone(), problem.source_conditions.clone(), problem.phase)
        .map_err(|source| SteadyGasInputError::Gas { line, source })
}

/// Read after all original side geometry and terminal declarations, before the
/// exterior radiation record. This selects the material source for the entire
/// open flow network; sealed finite branches retain their independently supplied
/// gas because a stagnant thermal region needs its own physical heat model.
#[allow(clippy::too_many_arguments)]
pub(crate) fn read<'a, E>(
    line: usize, fields: &mut SplitAsciiWhitespace<'a>,
    next: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
    pressure_pa: f64, segments: &[Segment], sides: &[GasSideGeometry],
    sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BranchedGasPerformance, SteadyGasReadError<E>> {
    poll(line, cancelled)?;
    let columns = PrimitiveGasColumns::read(fields, sources, MAX_STEADY_GAS_CELLS + sides.len() + 1)
        .map_err(|error| column_error(line, error))?;
    if fields.next().is_some() {
        return Err(bad(line, "unexpected branched gas header field").into());
    }
    let mut inlet_row = row(next, "gas_inlet")?;
    let inlet_mass_flow = inlet_row.positive()?;
    let inlet = columns.read_state(&mut inlet_row.words)
        .map_err(|error| column_error(inlet_row.line, error))?;
    inlet_row.finish()?;
    let mut bore = row(next, "gas_bore")?;
    let count = bore.count(MAX_STEADY_GAS_CELLS)?;
    let bore_outlet = bore.positive()?;
    bore.finish()?;
    if count == 0 { return Err(bad(line, "branched gas needs main-bore thermal cells").into()); }
    let mut geometry_steps = 0;
    let main = read_cells(next, segments, count, line, &mut geometry_steps, cancelled)?;
    let mut total_cells = count;
    let mut regions = Vec::with_capacity(sides.len());
    for side in sides {
        poll(line, cancelled)?;
        let mut branch = row(next, "gas_side")?;
        let port = branch.word()?.parse::<u32>()
            .map_err(|_| bad(branch.line, "gas side needs its original port ID"))?;
        let junction_m = branch.positive()?;
        let mass_flow_kg_s = branch.positive()?;
        let count = branch.count(MAX_STEADY_GAS_CELLS - total_cells)?;
        branch.finish()?;
        if port != side.port || side.cylinder.is_some() != (count > 0) {
            return Err(bad(branch.line, "gas sides must follow original open ports; only finite chimneys own thermal cells").into());
        }
        let donor = main.iter().position(|cell| cell.end_m == junction_m)
            .filter(|index| *index + 1 < main.len())
            .ok_or_else(|| bad(branch.line, "each gas side junction must be an interior main-bore cell boundary"))?;
        let cells = if let Some((radius, length)) = side.cylinder {
            read_cells(next, &[Segment::Cylinder { radius, length }], count,
                branch.line, &mut geometry_steps, cancelled)?
        } else { Vec::new() };
        regions.push(Region { port, junction_m, mass_flow_kg_s, cells,
            first: total_cells, donor });
        total_cells += count;
    }
    let total_outlet = bore_outlet + regions.iter().map(|region| region.mass_flow_kg_s).sum::<f64>();
    if !total_outlet.is_finite() || (inlet_mass_flow - total_outlet).abs()
        > 32.0 * f64::EPSILON * inlet_mass_flow.max(total_outlet)
    {
        return Err(bad(line, "inlet mass rate must equal all original bore and side discharges").into());
    }
    let mut solve = row(next, "gas_solve")?;
    let controls = SteadyGasSolve {
        caloric: GasEnergySolve { minimum_temperature_k: solve.positive()?,
            maximum_temperature_k: solve.positive()?, energy_tolerance_j_kg: solve.positive()?,
            maximum_iterations: solve.count(128)? },
        maximum_dt_s: solve.positive()?, maximum_relaxations: solve.count(1_000_000)?,
        maximum_substeps: solve.count(1_000_000)?, steady_power_tolerance_w: solve.positive()?,
        energy_tolerance_j: solve.scalar()?, maximum_mean_mach: solve.positive()?,
    };
    solve.finish()?;
    let mut cells = Vec::with_capacity(total_cells);
    append_cells(&mut cells, &main);
    for region in &regions { append_cells(&mut cells, &region.cells); }
    let mut flows = Vec::with_capacity(total_cells - 1);
    let mut connect = |from: usize, to: usize, kg_s: f64| {
        flows.push(MassFlow { id: flows.len() as u64,
            a: CellId(from as u64), b: CellId(to as u64), kg_s });
    };
    for index in 0..main.len() - 1 {
        let rate = bore_outlet + regions.iter().filter(|region| region.donor > index)
            .map(|region| region.mass_flow_kg_s).sum::<f64>();
        connect(index, index + 1, rate);
    }
    let mut outlets = vec![GasOutlet { cell: main.len() - 1, mass_flow_kg_s: bore_outlet }];
    for region in &regions {
        let mut donor = region.donor;
        for offset in 0..region.cells.len() {
            let cell = region.first + offset;
            connect(donor, cell, region.mass_flow_kg_s);
            donor = cell;
        }
        outlets.push(GasOutlet { cell: donor, mass_flow_kg_s: region.mass_flow_kg_s });
    }
    let (source, phase) = columns.into_parts();
    let problem = SteadyGasNetwork { source, source_conditions: QueryPoint::new(),
        static_pressure_pa: pressure_pa, phase, cells, flows,
        inlets: vec![GasInjection { cell: 0, mass_flow_kg_s: inlet_mass_flow, state: inlet.clone() }],
        outlets };
    let solution = problem.solve(controls, &mut *cancelled)
        .map_err(|source| SteadyGasInputError::SteadyGas { line, source })?;
    let main_field = field(&problem, &solution, &main, 0, inlet, line, cancelled)?;
    let mut branches = Vec::with_capacity(regions.len());
    for region in &regions {
        if region.cells.is_empty() { continue; }
        let donor = PrimitiveGasState {
            temperature_k: solution.states()[region.donor].state().temperature,
            composition: solution.caloric_laws()[region.donor].mole_composition().clone(),
        };
        branches.push(SteadyBranchField { port: region.port, junction_m: region.junction_m,
            field: field(&problem, &solution, &region.cells, region.first, donor, line, cancelled)? });
    }
    // No receipt is published before the complete network and regional fields pass.
    let mut accepted_records = Vec::new();
    for field in std::iter::once(&main_field).chain(branches.iter().map(|branch| &branch.field)) {
        for &x in field.stations_m() {
            poll(line, cancelled)?;
            accepted_records.push(GasRecord::Thermochemical(field.state_at(x)
                .map_err(|source| SteadyGasInputError::Gas { line, source })?));
        }
    }
    records.extend(accepted_records);
    Ok(BranchedGasPerformance { problem, controls, solution, main: main_field, branches,
        junctions: regions.iter().map(|region| (region.port, region.junction_m)).collect() })
}
