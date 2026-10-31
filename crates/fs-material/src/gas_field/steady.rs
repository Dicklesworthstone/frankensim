//! Steady axial gas heating, cooling and species injection.
//!
//! The existing Advected observer owns upwind enthalpy and stationary species
//! fluxes. The existing thermal timeline owns wall heat, CFL and energy updates.
//! Fixed reference cell masses are used ONLY as pseudo-time relaxation weights.
//! A converged state is reconstructed at common pressure through NASA-9/EOS and
//! becomes the same AxialGasField consumed by ducts, reeds and radiation.
//!
//! This model has prescribed positive through-flow and perfectly mixed axial
//! control volumes. Wall temperatures and total conductances are supplied
//! physical data. It does not predict pressure drop, wall heat-transfer
//! coefficients, turbulence, diffusion, reaction or condensation. No transient
//! density/mass history or finite-Mach acoustic accuracy is claimed.

/// Shared source-bound steady transport on branched and recirculating graphs.
pub mod network;

/// Equilibrium liquid-water removal through the same species/enthalpy observer.
pub mod condensation;

use super::{AxialGasField, PrimitiveGasState, resolve_primitive, typed};
use crate::{MaterialError, gas_species::input::ThermochemicalPhase};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::{GasEnergyError, GasEnergySolve}};
use fs_matdb::QueryPoint;
use fs_qty::semantic::QuantityKind;
use fs_thermal::enthalpy::{CellGeometry, CellId, DepositionError};
use fs_thermal::enthalpy::advection::MassFlow;
use fs_thermal::enthalpy::advection::species::StationarySpecies;
use network::{GasNetworkCell, GasOutlet, SteadyGasNetwork};

/// Bound the thermal/spatial cold path independently of audio sample count.
pub const MAX_STEADY_GAS_CELLS: usize = 64;
/// Extra gas supplies in addition to the upstream inlet.
pub const MAX_STEADY_GAS_INJECTIONS: usize = 64;

/// One finite volume, from the preceding endpoint (zero for the first) to end_m.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialGasCell {
    /// Downstream coordinate [m].
    pub end_m: f64,
    /// Physical volume [m3]; also sets the mean flow area volume/length.
    pub volume_m3: f64,
    /// Gas-facing wall area [m2].
    pub wall_area_m2: f64,
    /// Supplied TOTAL conductance to the wall [W/K], including zero insulation.
    /// A different material/flow requires the corresponding resolved conductance.
    pub wall_conductance_w_k: f64,
    /// Prescribed wall reservoir temperature [K].
    pub wall_temperature_k: f64,
    /// Initial temperature for pseudo-time relaxation, not an observed transient.
    pub initial_temperature_k: f64,
}

/// An explicit gas supply mixed into one control volume.
#[derive(Clone, Debug, PartialEq)]
pub struct GasInjection {
    /// Zero-based receiving cell. An axial path increases its downstream flow.
    pub cell: usize,
    /// Positive supplied mass rate [kg/s].
    pub mass_flow_kg_s: f64,
    /// Source temperature and full canonical mole composition.
    pub state: PrimitiveGasState,
}

/// Replaceable source, geometry, flow and wall data for one steady gas field.
#[derive(Clone, Debug, PartialEq)]
pub struct SteadyAxialGas {
    /// The single shared NASA-9/transport source set for every local mixture.
    pub source: ThermochemicalGasModel,
    /// Extra common material coordinates; temperature and pressure must be absent.
    pub source_conditions: QueryPoint,
    /// Common static pressure [Pa]; pressure drop is outside this model.
    pub static_pressure_pa: f64,
    /// Explicit gas-phase/saturation admission.
    pub phase: ThermochemicalPhase,
    /// Upstream source state.
    pub inlet: PrimitiveGasState,
    /// Positive upstream mass rate [kg/s].
    pub inlet_mass_flow_kg_s: f64,
    /// Ordered finite volumes.
    pub cells: Vec<AxialGasCell>,
    /// Optional explicitly supplied gases, preserving their own enthalpies.
    pub injections: Vec<GasInjection>,
}

/// Independent caloric, steady-residual and work admission controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SteadyGasSolve {
    /// Shared NASA-9 inversion bracket, absolute J/kg tolerance and iterations.
    pub caloric: GasEnergySolve,
    /// Maximum pseudo-time step [s]; the original combined CFL can tighten it.
    pub maximum_dt_s: f64,
    /// At most this many residual/relaxation requests.
    pub maximum_relaxations: usize,
    /// Total actual thermal substeps across all requests, never reset per request.
    pub maximum_substeps: usize,
    /// Maximum absolute net inward power at EVERY cell [W], including
    /// the conservative error bound from publishing caloric inverse states.
    pub steady_power_tolerance_w: f64,
    /// Per-request absolute energy-audit tolerance [J].
    pub energy_tolerance_j: f64,
    /// Caller-selected low-Mach admission, positive and at most 0.1.
    /// This is an explicit model limit, not a bound on omitted acoustic convection.
    pub maximum_mean_mach: f64,
}

/// No failed or unconverged field is published.
#[derive(Clone, Debug, PartialEq)]
pub enum SteadyGasError {
    /// Invalid or over-budget model inputs.
    Invalid(&'static str),
    /// Original material/source/phase refusal.
    Material(MaterialError),
    /// Original caloric inversion/source-bound refusal.
    Caloric(GasEnergyError),
    /// Original shared transport/integration refusal.
    Thermal(DepositionError),
    /// Requested steady tolerance was not reached within declared work.
    NotConverged { substeps: usize, maximum_power_residual_w: f64 },
    /// Caloric closure would publish a field outside the steady power tolerance.
    PublishedBalance { cell: usize, residual_bound_w: f64 },
    /// Prescribed flow leaves the declared low-Mach envelope.
    MeanFlow { cell: usize, mach: f64, limit: f64 },
}
impl core::fmt::Display for SteadyGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "steady axial gas: {self:?}")
    }
}
impl std::error::Error for SteadyGasError {}
impl From<MaterialError> for SteadyGasError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl From<GasEnergyError> for SteadyGasError {
    fn from(error: GasEnergyError) -> Self { Self::Caloric(error) }
}
impl From<DepositionError> for SteadyGasError {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}

/// A usable shared field and the actual conservative steady-solve observations.
#[derive(Clone, Debug, PartialEq)]
pub struct SteadyGasSolution {
    field: AxialGasField,
    species: StationarySpecies,
    mass_flow_kg_s: Vec<f64>,
    cell_power_residuals_w: Vec<f64>,
    maximum_power_residual_bound_w: f64,
    maximum_energy_residual_j: f64,
    substeps: usize,
    relaxation_time_s: f64,
    maximum_mean_mach: f64,
}
impl SteadyGasSolution {
    /// Feed directly to RegionalDuct/ReedBoreVoice/ReedPerformance field APIs.
    pub fn field(&self) -> &AxialGasField { &self.field }
    /// Move the same field into an existing consumer without reevaluation.
    pub fn into_field(self) -> AxialGasField { self.field }
    /// Full stationary species balances on the actual prescribed flow graph.
    pub fn species(&self) -> &StationarySpecies { &self.species }
    /// Downstream through-flow per axial cell [kg/s].
    pub fn mass_flow_kg_s(&self) -> &[f64] { &self.mass_flow_kg_s }
    /// Net inward watts at the emitted field states: the original owner balance
    /// plus the actual closure corrections on its upstream/outgoing enthalpy fluxes.
    pub fn cell_power_residuals_w(&self) -> &[f64] { &self.cell_power_residuals_w }
    /// Conservative maximum W residual including absolute closure-flux defects.
    /// This bounds BOTH the emitted-grid balance and the exact supplied-inlet
    /// balance; common inverse errors cannot cancel a source enthalpy mismatch.
    pub fn maximum_power_residual_bound_w(&self) -> f64 { self.maximum_power_residual_bound_w }
    /// Largest absolute per-request stored-minus-external energy defect [J].
    pub fn maximum_energy_residual_j(&self) -> f64 { self.maximum_energy_residual_j }
    /// Actual shared thermal substeps used to reach the residual tolerance.
    pub fn substeps(&self) -> usize { self.substeps }
    /// Pseudo-time relaxation coordinate; not physical elapsed breath time.
    pub fn relaxation_time_s(&self) -> f64 { self.relaxation_time_s }
    /// Largest final mean speed/sound-speed ratio over cell and inlet states.
    pub fn maximum_mean_mach(&self) -> f64 { self.maximum_mean_mach }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn bad(message: &'static str) -> SteadyGasError { SteadyGasError::Invalid(message) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DepositionError> {
    if cancelled() { Err(DepositionError::Cancelled) } else { Ok(()) }
}
fn point_at(context: &QueryPoint, t: f64) -> Result<QueryPoint, SteadyGasError> {
    context.clone().with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)
        .map_err(|_| bad("invalid temperature coordinate"))
}

impl SteadyAxialGas {
    /// Lower axial geometry, supplies and wall data into the same reusable gas
    /// network used by branched/recirculating systems. This validates geometry
    /// and flow bookkeeping; source and steady-state admission occurs in solve.
    pub fn to_network(&self) -> Result<SteadyGasNetwork, SteadyGasError> {
        if self.cells.is_empty() || self.cells.len() > MAX_STEADY_GAS_CELLS
            || self.injections.len() > MAX_STEADY_GAS_INJECTIONS
            || !positive(self.inlet_mass_flow_kg_s) || !positive(self.static_pressure_pa)
        {
            return Err(bad("bounded cells/injections and positive inlet flow/pressure required"));
        }
        if self.source_conditions.axes().contains_key("temperature")
            || self.source_conditions.axes().contains_key("pressure")
        {
            return Err(bad("steady gas owns temperature and common pressure"));
        }
        let mut previous = 0.0;
        for cell in &self.cells {
            if !cell.end_m.is_finite() || cell.end_m <= previous
                || !positive(cell.volume_m3) || !positive(cell.wall_area_m2)
                || !cell.wall_conductance_w_k.is_finite() || cell.wall_conductance_w_k < 0.0
                || !positive(cell.wall_temperature_k) || !positive(cell.initial_temperature_k)
            {
                return Err(bad("invalid axial cell geometry, wall data or initial temperature"));
            }
            previous = cell.end_m;
        }
        if self.injections.iter().any(|inlet|
            inlet.cell >= self.cells.len() || !positive(inlet.mass_flow_kg_s))
        {
            return Err(bad("gas injection needs an active axial cell and positive mass rate"));
        }
        let count = self.cells.len();
        let mut inlets = vec![GasInjection { cell: 0,
            mass_flow_kg_s: self.inlet_mass_flow_kg_s, state: self.inlet.clone() }];
        inlets.extend(self.injections.iter().cloned());
        let mut incoming_mass = vec![0.0; count];
        for inlet in &inlets {
            incoming_mass[inlet.cell] += inlet.mass_flow_kg_s;
        }
        let mut flows = Vec::with_capacity(count - 1);
        let mut rate = 0.0;
        for (index, &incoming) in incoming_mass.iter().enumerate() {
            rate += incoming;
            if !positive(rate) { return Err(bad("axial mass flow is not representable")); }
            if index + 1 < count {
                flows.push(MassFlow { id: index as u64,
                    a: CellId(index as u64), b: CellId((index + 1) as u64), kg_s: rate });
            }
        }
        let mut cells = Vec::with_capacity(count);
        previous = 0.0;
        for cell in &self.cells {
            cells.push(GasNetworkCell {
                geometry: CellGeometry { start_m: [previous, 0.0, 0.0], end_m: [cell.end_m, 0.0, 0.0],
                    volume_m3: cell.volume_m3, surface_area_m2: cell.wall_area_m2 },
                flow_area_m2: cell.volume_m3 / (cell.end_m - previous),
                wall_conductance_w_k: cell.wall_conductance_w_k,
                wall_temperature_k: cell.wall_temperature_k,
                initial_temperature_k: cell.initial_temperature_k,
            });
            previous = cell.end_m;
        }
        Ok(SteadyGasNetwork { source: self.source.clone(),
            source_conditions: self.source_conditions.clone(), static_pressure_pa: self.static_pressure_pa,
            phase: self.phase, cells, flows, inlets,
            outlets: vec![GasOutlet { cell: count - 1, mass_flow_kg_s: rate }] })
    }

    /// Solve through the shared gas network and publish a source-bound axial
    /// field for existing ducts, reeds and radiation. Each downstream face takes
    /// its cell's well-mixed state; spatial refinement remains independent of
    /// thermal convergence and of the consumer's duct discretization.
    pub fn solve(
        &self, controls: SteadyGasSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<SteadyGasSolution, SteadyGasError> {
        poll(&mut cancelled)?;
        let network = self.to_network()?;
        let mut solution = network.solve(controls, &mut cancelled)?;
        // The upstream boundary is an enthalpy-balanced mixture of all supplies
        // entering cell zero, before any of its wall heat is exchanged.
        let upstream_h = solution.inlet_enthalpy_w[0] / solution.inlet_mass_kg_s[0];
        let upstream = solution.laws[0].resolve_enthalpy(upstream_h, &mut cancelled)?;
        let mut stations = vec![0.0];
        let mut states = vec![PrimitiveGasState {
            temperature_k: upstream.state().state().temperature,
            composition: solution.laws[0].mole_composition().clone(),
        }];
        // The network checks the exact supplied enthalpies. The axial field
        // additionally publishes their reconstructed mixture at its x=0 face.
        let correction = solution.inlet_mass_kg_s[0] * upstream.residual_j_kg();
        solution.cell_power_residuals_w[0] += correction;
        solution.cell_power_residual_bounds_w[0] += correction.abs();
        let bound = solution.cell_power_residual_bounds_w[0];
        if !bound.is_finite() || bound > controls.steady_power_tolerance_w {
            return Err(SteadyGasError::PublishedBalance { cell: 0, residual_bound_w: bound });
        }
        let gas = upstream.state().state();
        let sonic_mass_rate = gas.density * network.cells[0].flow_area_m2 * gas.sound_speed;
        if !positive(sonic_mass_rate) { return Err(bad("inlet sonic mass rate is not representable")); }
        let mach = solution.inlet_mass_kg_s[0] / sonic_mass_rate;
        if !positive(mach) || mach > controls.maximum_mean_mach {
            return Err(SteadyGasError::MeanFlow { cell: 0, mach, limit: controls.maximum_mean_mach });
        }
        solution.maximum_mean_mach = solution.maximum_mean_mach.max(mach);
        for (index, cell) in self.cells.iter().enumerate() {
            poll(&mut cancelled)?;
            stations.push(cell.end_m);
            states.push(PrimitiveGasState { temperature_k: solution.states()[index].state().temperature,
                composition: solution.laws[index].mole_composition().clone() });
        }
        poll(&mut cancelled)?;
        let field = AxialGasField::new(stations, self.static_pressure_pa, states,
            self.source.clone(), self.source_conditions.clone(), self.phase)?;
        poll(&mut cancelled)?;
        let maximum_power_residual_bound_w = solution.maximum_power_residual_bound_w();
        Ok(SteadyGasSolution { field, species: solution.species, mass_flow_kg_s: solution.mass_flow_kg_s,
            cell_power_residuals_w: solution.cell_power_residuals_w, maximum_power_residual_bound_w,
            maximum_energy_residual_j: solution.maximum_energy_residual_j, substeps: solution.substeps,
            relaxation_time_s: solution.relaxation_time_s, maximum_mean_mach: solution.maximum_mean_mach })
    }
}

#[cfg(test)]
pub(crate) mod tests;
