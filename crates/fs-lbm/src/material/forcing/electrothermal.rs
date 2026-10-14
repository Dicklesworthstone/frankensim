//! DC Joule heating coupled to the existing material-bound thermal/flow step.
//!
//! The electrical mesh is derived from the slab, not supplied independently.
//! Current balance, bulk/contact heat partition, transport, buoyancy and wall
//! exchange retain their existing owners. No material presets are inferred.
//! Electrical outer faces are insulated; lateral fluid boundaries are periodic.
//! Electrodes are ideal *electrical* cell constraints, not mechanical obstacles.
//! This is ohmic sensible heating, not electrochemistry, motional EMF, magnetic
//! force, phase change, compressible work or a total-energy model.

pub mod midpoint;

use fs_exec::CancelGate;
use fs_material::fluid_transport::FluidTransportState;

use super::{MaterialSlab, ThermalLbm, ThermalLbm3};
use crate::conjugate::{ChtError, ElectricSetup, ElectricSolution, ElectricalContact, VoxelDomain};
use crate::thermal::sources::{HeatBalance, HeatSourceLimits, Slab, ThermalSourceError};

/// Limits on the existing thermal step and the electrical-to-heat handoff.
#[derive(Debug, Clone, Copy)]
pub struct ElectricHeatLimits {
    /// Cell, absolute-temperature and thermal balance bounds.
    pub thermal: HeatSourceLimits,
    /// Absolute tolerance [J] for each electrical power/deposition and transfer
    /// check. Combined balance uses `thermal.balance_tolerance_j + 2*this`.
    pub coupling_tolerance_j: f64,
}

/// Electrical inputs for one step. Borrowing them does not persist a drive.
#[derive(Clone, Copy)]
pub struct ElectricHeatInput<'a> {
    /// Electrodes and their drives. The legacy solid preset table is ignored:
    /// the callback alone declares bulk conduction versus insulation.
    pub setup: &'a ElectricSetup,
    /// Explicit zero-thickness resistance and heat partition. Already included
    /// in the electrical solver's cell powers; never added a second time.
    pub contacts: &'a [ElectricalContact],
    /// Numerical bounds, including both thermal wall planes in the cell count.
    pub limits: ElectricHeatLimits,
}

/// Accepted numerical fields and independently checked energy transfers.
#[derive(Debug, Clone, PartialEq)]
pub struct ElectricHeatStep {
    /// The actual potential/current solve used to generate this step's heat.
    pub electric: ElectricSolution,
    /// Actual sensible-heat storage and independently measured wall exchange.
    pub heat: HeatBalance,
    /// `sum_e voltage * current * dt` [J].
    pub electrode_energy_j: f64,
    /// Sum of the electrical solver's cell heat powers times `dt` [J].
    pub joule_energy_j: f64,
    /// Electrode energy minus deposited electrical energy [J].
    pub electrical_defect_j: f64,
    /// Thermal source energy minus deposited electrical energy [J].
    pub transfer_defect_j: f64,
    /// Stored change minus wall inputs minus electrode energy [J].
    pub combined_defect_j: f64,
}

/// Refusal leaves the slab's populations, coefficients and material units intact.
/// A caller's externally captured constitutive receipts are not rolled back.
#[derive(Debug, Clone, PartialEq)]
pub enum ElectricHeatError {
    /// Bad geometry, limits, sources, or a failed energy handoff.
    Invalid(&'static str),
    /// Existing current solver or the caller's material law refused the state.
    Electric(ChtError),
    /// Existing material-bound thermal step refused the candidate.
    Thermal(ThermalSourceError),
    /// Cancelled before publishing the candidate.
    Cancelled,
}
impl core::fmt::Display for ElectricHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "electric/thermal flow: {s}"),
            Self::Electric(e) => write!(f, "electric/thermal flow: {e}"),
            Self::Thermal(e) => write!(f, "electric/thermal flow: {e}"),
            Self::Cancelled => write!(f, "electric/thermal flow cancelled"),
        }
    }
}
impl std::error::Error for ElectricHeatError {}
impl From<ChtError> for ElectricHeatError {
    fn from(error: ChtError) -> Self {
        if error == ChtError::Cancelled { Self::Cancelled } else { Self::Electric(error) }
    }
}
impl From<ThermalSourceError> for ElectricHeatError {
    fn from(error: ThermalSourceError) -> Self { Self::Thermal(error) }
}

fn poll(gate: &CancelGate) -> Result<(), ElectricHeatError> {
    if gate.is_requested() { Err(ElectricHeatError::Cancelled) } else { Ok(()) }
}
fn invalid(s: &'static str) -> ElectricHeatError { ElectricHeatError::Invalid(s) }
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }

fn temperatures<S: Slab<D, Q>, const D: usize, const Q: usize>(
    slab: &MaterialSlab<S>, limits: HeatSourceLimits, gate: &CancelGate,
) -> Result<Vec<f64>, ElectricHeatError> {
    let [lo, hi] = limits.temperature_range_k;
    slab.lattice.populations().iter().enumerate().map(|(cell, p)| {
        if cell % 256 == 0 { poll(gate)?; }
        if !p.iter().all(|x| x.is_finite()) { return Err(invalid("nonfinite temperature population")); }
        let t = slab.temperature_si(p.iter().sum())
            .map_err(ThermalSourceError::from)?;
        if t < lo || t > hi { return Err(invalid("temperature outside the admitted material interval")); }
        Ok(t)
    }).collect()
}

fn admit<S: Slab<D, Q>, const D: usize, const Q: usize>(
    slab: &MaterialSlab<S>, input: ElectricHeatInput<'_>, forces: &[[f64; D]],
    normal: usize, normal_cells: usize, gate: &CancelGate,
) -> Result<(usize, usize), ElectricHeatError> {
    poll(gate)?;
    let (cells, plane) = slab.lattice.layout()?;
    let limits = input.limits;
    let [lo, hi] = limits.thermal.temperature_range_k;
    let combined = limits.thermal.balance_tolerance_j + 2.0 * limits.coupling_tolerance_j;
    if cells > limits.thermal.max_cells || forces.len() != cells
        || !forces.iter().flatten().all(|v| v.is_finite())
        || !positive(lo) || !positive(hi) || lo > hi
        || !limits.thermal.balance_tolerance_j.is_finite() || limits.thermal.balance_tolerance_j < 0.0
        || !limits.coupling_tolerance_j.is_finite() || limits.coupling_tolerance_j < 0.0
        || !combined.is_finite()
    {
        return Err(invalid("invalid force field, temperature/energy limits or cell budget"));
    }
    if slab.lattice.walls() != (1.0, 0.0) || !slab.lattice.flow_valid() {
        return Err(invalid("invalid material slab or independently changed normalized walls"));
    }
    // Fixed-temperature wall cells have no evolved storage in this solver.
    // Do not silently discard electrical/contact heating into a reservoir.
    if input.setup.electrodes.iter().any(|e| e.lo[normal] == 0 || e.hi[normal] >= normal_cells) {
        return Err(invalid("electrodes must lie strictly between the thermal wall planes"));
    }
    Ok((cells, plane))
}

struct JouleSource {
    electric: ElectricSolution,
    power_w_m3: Vec<f64>,
    electrode_j: f64,
    joule_j: f64,
    electrical_defect_j: f64,
}

#[allow(clippy::too_many_arguments)] // geometry, field and constitutive law are independent inputs
fn solve_source(
    domain: &VoxelDomain, plane: usize, temperature: &[f64],
    transport: &FluidTransportState, dt: f64, input: ElectricHeatInput<'_>,
    gate: &CancelGate,
    resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
) -> Result<JouleSource, ElectricHeatError> {
    poll(gate)?;
    let cells = domain.cell_count();
    let electric = input.setup.solve_with_cell_contacts(
        domain, temperature, input.contacts,
        |cell, _, t| {
            if cell < plane || cell >= cells - plane { Ok(None) }
            else { resolve(cell, t, transport) }
        }, gate,
    )?;
    if electric.joule_w.len() != cells || !electric.total_joule_w.is_finite()
        || !electric.delivered_w.is_finite()
    { return Err(invalid("electrical output has invalid shape or total power")); }
    let mut power_w_m3 = Vec::with_capacity(cells);
    let mut joule_j = 0.0;
    for (cell, &watts) in electric.joule_w.iter().enumerate() {
        if cell % 256 == 0 { poll(gate)?; }
        let volume = domain.volume(cell);
        if !watts.is_finite() || watts < 0.0 || !positive(volume)
            || ((cell < plane || cell >= cells-plane) && watts != 0.0)
        { return Err(invalid("invalid electrical heat or heat deposited in a fixed-temperature wall")); }
        let density = watts / volume;
        let energy = watts * dt;
        if !density.is_finite() || !energy.is_finite()
            || (watts != 0.0 && (density == 0.0 || energy == 0.0))
        { return Err(invalid("electrical heating is not representable in physical thermal units")); }
        power_w_m3.push(density);
        joule_j += energy;
    }
    let electrode_j: f64 = electric.electrodes.iter().map(|e| (e.voltage_v * e.current_a) * dt).sum();
    let electrical_defect_j = electrode_j - joule_j;
    let tol = input.limits.coupling_tolerance_j;
    if ![joule_j, electrode_j, electrical_defect_j].iter().all(|x| x.is_finite())
        || electrical_defect_j.abs() > tol
        || (electric.total_joule_w * dt - joule_j).abs() > tol
        || (electric.delivered_w * dt - electrode_j).abs() > tol
    { return Err(invalid("electrode work and deposited Joule heat do not close")); }
    Ok(JouleSource { electric, power_w_m3, electrode_j, joule_j, electrical_defect_j })
}

fn balance(source: JouleSource, heat: HeatBalance, limits: ElectricHeatLimits) -> Result<ElectricHeatStep, ElectricHeatError> {
    let transfer_defect_j = heat.source_j - source.joule_j;
    let combined_defect_j = heat.stored_change_j - heat.bottom_wall_j - heat.top_wall_j - source.electrode_j;
    if !transfer_defect_j.is_finite() || transfer_defect_j.abs() > limits.coupling_tolerance_j
        || !combined_defect_j.is_finite()
        || combined_defect_j.abs() > limits.thermal.balance_tolerance_j + 2.0 * limits.coupling_tolerance_j
    { return Err(invalid("electrical-to-thermal transfer or combined energy balance failed")); }
    Ok(ElectricHeatStep {
        electric: source.electric, heat, electrode_energy_j: source.electrode_j,
        joule_energy_j: source.joule_j, electrical_defect_j: source.electrical_defect_j,
        transfer_defect_j, combined_defect_j,
    })
}

fn candidate<S: Slab<D, Q>, const D: usize, const Q: usize>(slab: &MaterialSlab<S>) -> MaterialSlab<S> {
    MaterialSlab {
        lattice: slab.lattice.snapshot().0,
        top_temperature_k: slab.top_temperature_k,
        temperature_difference_k: slab.temperature_difference_k,
        force_units: slab.force_units,
    }
}

#[allow(clippy::too_many_arguments)] // shared 2D/3D orchestration without duplicating solvers
fn advance<S: Slab<D, Q>, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, forces: &[[f64; D]], gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    step: impl FnOnce(&mut MaterialSlab<S>, &[f64], &[[f64; D]], HeatSourceLimits) -> Result<HeatBalance, ThermalSourceError>,
) -> Result<ElectricHeatStep, ElectricHeatError> {
    let temperature = temperatures::<S, D, Q>(slab, input.limits.thermal, gate)?;
    let source = solve_source(&domain, plane, &temperature, &slab.force_units.transport,
        slab.force_units.timestep_s, input, gate, &mut resolve)?;
    poll(gate)?;
    let mut next = candidate::<S, D, Q>(slab);
    let heat = step(&mut next, &source.power_w_m3, forces, input.limits.thermal)?;
    let report = balance(source, heat, input.limits)?;
    poll(gate)?;
    *slab = next;
    Ok(report)
}

impl MaterialSlab<ThermalLbm> {
    /// Solve DC current at the accepted temperatures, deposit its bulk AND
    /// contact Joule heat, then advance the existing 2-D heat/flow equations.
    /// `resolve(cell, kelvin, retained_transport)` supplies resistivity [ohm m],
    /// or explicit insulation (`None`); missing material data must be an error.
    /// Only non-electrode fluid cells invoke the callback. It owns source pins,
    /// pressure/composition coordinates, and any externally retained receipts.
    ///
    /// The electrical domain is `nx * ny * 1`, with the supplied depth [m].
    /// Electrode and contact cell indices are the existing slab's linear indices.
    /// Electrical outer faces are insulated, even where flow is periodic.
    /// Electrical conductivity can change each step; thermal transport remains
    /// the slab's retained state. This is explicit coupling, not an implicit
    /// temperature/resistivity solve or a second-order temporal accuracy claim.
    ///
    /// Admission, electrical, thermal, transfer and cancellation failures do not
    /// publish a partial step. Peak storage includes original, candidate and the
    /// thermal owner's tentative grid. Cancellation is polled around complete
    /// thermal steps and within the electrical solve, not inside every LBM cell.
    pub fn step_electrothermal(
        &mut self, input: ElectricHeatInput<'_>, forces_n_m3: &[[f64; 2]],
        depth_m: f64, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<ElectricHeatStep, ElectricHeatError> {
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces_n_m3, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) { return Err(invalid("invalid represented 2D cell volume")); }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ])?;
        advance::<ThermalLbm, 2, 5>(self, domain, plane, input, forces_n_m3, gate, resolve,
            |s, p, f, limits| s.step_with_power_density(p, f, depth_m, limits))
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart of [`MaterialSlab::step_electrothermal`].
    /// Electrical and thermal cells share x-fastest indexing and the retained
    /// `dx^3` volume. All source, rollback and boundary semantics are shared.
    pub fn step_electrothermal(
        &mut self, input: ElectricHeatInput<'_>, forces_n_m3: &[[f64; 3]],
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<Option<f64>, ChtError>,
    ) -> Result<ElectricHeatStep, ElectricHeatError> {
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces_n_m3, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume")); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, 3, 7>(self, domain, plane, input, forces_n_m3, gate, resolve,
            |s, p, f, limits| s.step_with_power_density(p, f, limits))
    }
}
