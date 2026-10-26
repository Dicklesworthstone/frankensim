//! Source-bound steady gas transport on a prescribed branched flow network.
//!
//! Species balance and enthalpy transport use the same `fs-thermal` graph.
//! Each accepted mixture retains the original NASA-9 and transport sources;
//! density, sound speed and viscosity are consequences of that one state.
//! This is fixed-pressure, well-mixed, low-Mach steady transport. Prescribed
//! flows and wall conductances are inputs, not predictions of a pressure-drop
//! or turbulence model. Reference cell masses only set pseudo-time relaxation.

use super::{GasInjection, MAX_STEADY_GAS_CELLS, MAX_STEADY_GAS_INJECTIONS,
    PrimitiveGasState, SteadyGasError, SteadyGasSolve, bad, point_at, poll, positive,
    resolve_primitive, typed};
use crate::gas_species::input::ThermochemicalPhase;
use crate::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};
use crate::thermochemical_gas::energy::enthalpy::GasEnthalpyLaw;
use fs_matdb::QueryPoint;
use fs_qty::semantic::QuantityKind;
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositedCell,
    DepositionBudget, DepositionSimulation, Event, MaterialId, Reservoir,
    StepConfig, ThermalLaw, TimedEvent};
use fs_thermal::enthalpy::advection::{Advected, MassBoundary, MassFlow};
use fs_thermal::enthalpy::advection::species::{InletSpecies, StationarySpecies};
use fs_thermochem::{Composition, CompositionBasis};

/// Bound all interior and boundary connections in one cold-path solve.
pub const MAX_STEADY_GAS_CONNECTIONS: usize = 256;

/// One well-mixed control volume with an explicit bulk-flow area.
#[derive(Clone, Debug, PartialEq)]
pub struct GasNetworkCell {
    /// Existing physical volume, wall area and geometric support.
    pub geometry: CellGeometry,
    /// Effective cross section for the cell's total outgoing bulk flow [m2].
    /// Branch/nozzle throat speeds require their own geometry and model.
    pub flow_area_m2: f64,
    /// Supplied total wall conductance [W/K], including zero insulation.
    pub wall_conductance_w_k: f64,
    /// Prescribed wall reservoir temperature [K].
    pub wall_temperature_k: f64,
    /// Initial pseudo-time relaxation temperature [K].
    pub initial_temperature_k: f64,
}

/// Prescribed discharge; its enthalpy is always the current local gas enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasOutlet {
    /// Zero-based source cell.
    pub cell: usize,
    /// Positive outward mass rate [kg/s].
    pub mass_flow_kg_s: f64,
}

/// Data-only material, geometry, flow and wall description.
#[derive(Clone, Debug, PartialEq)]
pub struct SteadyGasNetwork {
    /// Common original species thermodynamics and material transport sources.
    pub source: ThermochemicalGasModel,
    /// Fixed external coordinates, excluding owned temperature and pressure.
    pub source_conditions: QueryPoint,
    /// Common static pressure [Pa]. No pressure-drop law is implied.
    pub static_pressure_pa: f64,
    /// Shared gas-phase and saturation policy.
    pub phase: ThermochemicalPhase,
    /// Cells numbered `CellId(0)..CellId(cells.len()-1)`.
    pub cells: Vec<GasNetworkCell>,
    /// The original thermal owner's signed flows and stable connection IDs.
    /// Boundary IDs are allocated after the largest interior ID.
    pub flows: Vec<MassFlow>,
    /// Explicit supplies with their own compositions and specific enthalpies.
    pub inlets: Vec<GasInjection>,
    /// Explicit discharges. Every active cell must balance total mass.
    pub outlets: Vec<GasOutlet>,
}

/// Complete material states and independently checked steady transport balances.
#[derive(Clone, Debug, PartialEq)]
pub struct SteadyGasNetworkSolution {
    states: Vec<ResolvedThermochemicalGasState>,
    pub(super) laws: Vec<GasEnthalpyLaw>,
    pub(super) species: StationarySpecies,
    pub(super) mass_flow_kg_s: Vec<f64>,
    pub(super) cell_power_residuals_w: Vec<f64>,
    pub(super) cell_power_residual_bounds_w: Vec<f64>,
    pub(super) maximum_energy_residual_j: f64,
    pub(super) substeps: usize,
    pub(super) relaxation_time_s: f64,
    pub(super) maximum_mean_mach: f64,
    pub(super) inlet_mass_kg_s: Vec<f64>,
    pub(super) inlet_enthalpy_w: Vec<f64>,
}

impl SteadyGasNetworkSolution {
    /// Source-resolved gas states in cell order, usable by existing physics consumers.
    pub fn states(&self) -> &[ResolvedThermochemicalGasState] { &self.states }
    /// Frozen source-backed caloric laws in cell order, with full species order.
    pub fn caloric_laws(&self) -> &[GasEnthalpyLaw] { &self.laws }
    /// Stationary species balance on the same mass-flow graph.
    pub fn species(&self) -> &StationarySpecies { &self.species }
    /// Total outgoing flow per cell, including external discharge [kg/s].
    pub fn mass_flow_kg_s(&self) -> &[f64] { &self.mass_flow_kg_s }
    /// Net inward power at each published state [W].
    pub fn cell_power_residuals_w(&self) -> &[f64] { &self.cell_power_residuals_w }
    /// Bounds include each caloric reconstruction's absolute flux correction.
    pub fn cell_power_residual_bounds_w(&self) -> &[f64] { &self.cell_power_residual_bounds_w }
    /// Largest per-cell conservative steady residual bound [W].
    pub fn maximum_power_residual_bound_w(&self) -> f64 {
        self.cell_power_residual_bounds_w.iter().copied().fold(0.0, f64::max)
    }
    /// Largest actual thermal stored-minus-external energy residual [J].
    pub fn maximum_energy_residual_j(&self) -> f64 { self.maximum_energy_residual_j }
    /// Actual thermal substeps, bounded across the complete solve.
    pub fn substeps(&self) -> usize { self.substeps }
    /// Pseudo-time, not a physical density or mass history [s].
    pub fn relaxation_time_s(&self) -> f64 { self.relaxation_time_s }
    /// Largest cell bulk-flow Mach number, using each cell's complete gas state.
    pub fn maximum_mean_mach(&self) -> f64 { self.maximum_mean_mach }
}

impl SteadyGasNetwork {
    /// Solve stationary species, then relax the source-specific enthalpy/wall
    /// balance through the existing conservative thermal integrator. Branched
    /// and inlet-fed recirculating graphs use the same laws as straight ducts.
    /// An unanchored closed circulation has no prescribed mixture and refuses.
    ///
    /// No partial, unconverged, supersaturated or out-of-source state is returned.
    /// The residual is a balance check, not a spatial discretization error bound.
    pub fn solve(
        &self, controls: SteadyGasSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<SteadyGasNetworkSolution, SteadyGasError> {
        poll(&mut cancelled)?;
        let count = self.cells.len();
        if count == 0 || count > MAX_STEADY_GAS_CELLS
            || self.inlets.is_empty() || self.inlets.len() > MAX_STEADY_GAS_INJECTIONS + 1
            || self.flows.len().saturating_add(self.inlets.len())
                .saturating_add(self.outlets.len()) > MAX_STEADY_GAS_CONNECTIONS
            || !positive(self.static_pressure_pa) || !positive(controls.maximum_dt_s)
            || !positive(controls.steady_power_tolerance_w)
            || !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j < 0.0
            || controls.maximum_relaxations == 0 || controls.maximum_relaxations > 1_000_000
            || controls.maximum_substeps == 0 || controls.maximum_substeps > 1_000_000
            || !positive(controls.maximum_mean_mach) || controls.maximum_mean_mach > 0.1
        {
            return Err(bad("bounded network/work, explicit supplies, positive pressure/tolerances and a low-Mach ceiling required"));
        }
        if self.source_conditions.axes().contains_key("temperature")
            || self.source_conditions.axes().contains_key("pressure")
        {
            return Err(bad("steady gas owns temperature and common pressure"));
        }
        for cell in &self.cells {
            if !positive(cell.flow_area_m2) || !positive(cell.geometry.volume_m3)
                || !positive(cell.geometry.surface_area_m2)
                || !cell.wall_conductance_w_k.is_finite() || cell.wall_conductance_w_k < 0.0
                || !positive(cell.wall_temperature_k) || !positive(cell.initial_temperature_k)
            {
                return Err(bad("invalid gas cell geometry, wall data or initial temperature"));
            }
        }
        if self.inlets.iter().any(|inlet| inlet.cell >= count || !positive(inlet.mass_flow_kg_s))
            || self.outlets.iter().any(|outlet| outlet.cell >= count || !positive(outlet.mass_flow_kg_s))
            || self.flows.iter().any(|flow| flow.a.0 >= count as u64 || flow.b.0 >= count as u64)
        {
            return Err(bad("gas connections require active cells and positive boundary mass rates"));
        }
        let pressure_point = self.source_conditions.clone()
            .with_quantity("pressure", typed(QuantityKind::Pressure), self.static_pressure_pa)
            .map_err(|_| bad("invalid common pressure coordinate"))?;
        let mut next_id = self.flows.iter().map(|flow| flow.id).max()
            .map_or(Some(0), |id| id.checked_add(1))
            .ok_or(bad("no remaining gas boundary connection identities"))?;
        let mut boundary_id = || {
            let id = next_id;
            next_id = next_id.checked_add(1)
                .ok_or(bad("gas boundary connection identity overflow"))?;
            Ok::<_, SteadyGasError>(id)
        };
        let mut inlet_laws = Vec::with_capacity(self.inlets.len());
        let mut boundaries = Vec::with_capacity(self.inlets.len() + self.outlets.len());
        let mut inlet_species = Vec::with_capacity(self.inlets.len());
        let mut inlet_mass = vec![0.0; count];
        let mut inlet_enthalpy = vec![0.0; count];
        for (index, inlet) in self.inlets.iter().enumerate() {
            poll(&mut cancelled)?;
            resolve_primitive(&self.source, &pressure_point, self.phase, &inlet.state)?;
            let law = GasEnthalpyLaw::from_mole(&self.source, inlet.state.composition.clone(),
                point_at(&pressure_point, inlet.state.temperature_k)?, controls.caloric)?;
            let h = law.enthalpy_at_temperature(inlet.state.temperature_k)?;
            let id = boundary_id()?;
            boundaries.push(MassBoundary::Inlet { id, cell: CellId(inlet.cell as u64),
                kg_s: inlet.mass_flow_kg_s, material_id: MaterialId((count + index) as u64),
                specific_h_j_kg: h });
            inlet_species.push(InletSpecies { boundary_id: id,
                mass_fractions: law.advected_composition().expect("gas composition").mass_fractions.to_vec() });
            inlet_mass[inlet.cell] += inlet.mass_flow_kg_s;
            inlet_enthalpy[inlet.cell] += inlet.mass_flow_kg_s * h;
            inlet_laws.push(law);
        }
        for outlet in &self.outlets {
            boundaries.push(MassBoundary::Outlet { id: boundary_id()?,
                cell: CellId(outlet.cell as u64), kg_s: outlet.mass_flow_kg_s });
        }
        if inlet_mass.iter().chain(&inlet_enthalpy).any(|value| !value.is_finite()) {
            return Err(bad("supplied gas mass or enthalpy rate is not representable"));
        }
        let species = self.source.thermodynamics().components().iter()
            .map(|component| component.species().to_string()).collect();
        let observer = Advected::new((), self.flows.clone(), MAX_STEADY_GAS_CONNECTIONS)?
            .with_boundaries(boundaries, MAX_STEADY_GAS_CONNECTIONS)?
            .with_stationary_species(species, inlet_species, count, 128)?;
        let stationary = observer.stationary_species().expect("just admitted").clone();
        if stationary.cells().len() != count {
            return Err(bad("every steady gas cell needs inlet-connected positive through-flow"));
        }
        let mut laws = Vec::with_capacity(count + inlet_laws.len());
        for (index, cell) in self.cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let fractions = &stationary.cells()[&CellId(index as u64)];
            // Uniform runs retain the authored mole coefficients bit for bit.
            let same = inlet_laws.iter().find(|law|
                law.advected_composition().expect("gas law").mass_fractions == fractions.as_slice());
            let law = if let Some(law) = same { law.clone() } else {
                GasEnthalpyLaw::new(&self.source,
                    Composition::new(CompositionBasis::MassFraction, fractions.clone())
                        .map_err(|_| bad("stationary composition refused"))?,
                    point_at(&pressure_point, cell.initial_temperature_k)?, controls.caloric)?
            };
            laws.push(law);
        }
        laws.extend(inlet_laws);
        let mut schedule = Vec::with_capacity(count);
        for (index, cell) in self.cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let gas = resolve_primitive(&self.source, &pressure_point, self.phase,
                &PrimitiveGasState { temperature_k: cell.initial_temperature_k,
                    composition: laws[index].mole_composition().clone() })?;
            let mass_kg = gas.state().density * cell.geometry.volume_m3;
            if !positive(mass_kg) { return Err(bad("initial EOS relaxation mass is not representable")); }
            schedule.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(index as u64), segment_id: index as u64, body_id: 0,
                material_id: MaterialId(index as u64), mass_kg,
                specific_h_j_kg: laws[index].enthalpy_at_temperature(cell.initial_temperature_k)?,
                geometry: cell.geometry,
                boundary: Boundary { ambient_k: cell.wall_temperature_k, htc_w_m2_k: 0.0,
                    reservoir: Some(Reservoir { area_m2: cell.geometry.surface_area_m2,
                        conductance_w_k: cell.wall_conductance_w_k,
                        temperature_k: cell.wall_temperature_k }), source_w: 0.0 },
            }) });
        }
        let references: Vec<_> = laws.iter().enumerate()
            .map(|(index, law)| (MaterialId(index as u64), law as &dyn ThermalLaw)).collect();
        let budget = DepositionBudget { max_cells: laws.len(), max_links: self.flows.len(),
            max_events: count, max_substeps: controls.maximum_substeps };
        let step = StepConfig { max_dt_s: controls.maximum_dt_s,
            energy_tolerance_j: controls.energy_tolerance_j };
        let (mut thermal, _) = DepositionSimulation::new(schedule, &references, budget)?
            .advance_to(0.0, step, || poll(&mut cancelled))?;
        let mut observer = observer;
        let mut substeps = 0usize;
        let mut maximum_energy_residual_j: f64 = 0.0;
        let mut balance = thermal.power_balance_observed(&observer, || poll(&mut cancelled))?;
        let mut maximum_power = balance.net_power_w.values().fold(0.0f64, |a, b| a.max(b.abs()));
        let mut relaxations = 0usize;
        while maximum_power > controls.steady_power_tolerance_w {
            if relaxations >= controls.maximum_relaxations || substeps >= controls.maximum_substeps {
                return Err(SteadyGasError::NotConverged { substeps, maximum_power_residual_w: maximum_power });
            }
            let target = thermal.time_s() + controls.maximum_dt_s;
            if !target.is_finite() || target <= thermal.time_s() {
                return Err(bad("pseudo-time increment is not representable"));
            }
            let (candidate, prepared, audit) = thermal.advance_observed_limited(
                &observer, target, step, controls.maximum_substeps - substeps,
                || poll(&mut cancelled))?;
            substeps += audit.substeps;
            maximum_energy_residual_j = maximum_energy_residual_j.max(audit.residual_j.abs());
            thermal = candidate;
            observer = prepared;
            relaxations += 1;
            balance = thermal.power_balance_observed(&observer, || poll(&mut cancelled))?;
            maximum_power = balance.net_power_w.values().fold(0.0f64, |a, b| a.max(b.abs()));
        }
        let mut states = Vec::with_capacity(count);
        let mut closure = Vec::with_capacity(count);
        for (index, law) in laws[..count].iter().enumerate() {
            poll(&mut cancelled)?;
            let h = thermal.cells()[&CellId(index as u64)].specific_h_j_kg;
            let resolved = law.resolve_enthalpy(h, &mut cancelled)?;
            let gas = resolve_primitive(&self.source, &pressure_point, self.phase,
                &PrimitiveGasState { temperature_k: resolved.state().state().temperature,
                    composition: law.mole_composition().clone() })?;
            closure.push(gas.specific_enthalpy_j_kg() - h);
            states.push(gas);
        }
        // Correct the original balance using every actual published donor state.
        // Internal energy transfers appear with opposite signs at their two ends.
        // Absolute corrections are checked separately so common caloric errors
        // cannot conceal a mismatch against exact supplied inlet enthalpies.
        let mut outgoing = vec![0.0; count];
        let mut corrections = vec![0.0; count];
        let mut closure_bounds = vec![0.0; count];
        for flow in observer.flows() {
            let (up, down, rate) = if flow.kg_s >= 0.0 {
                (flow.a.0 as usize, flow.b.0 as usize, flow.kg_s)
            } else { (flow.b.0 as usize, flow.a.0 as usize, -flow.kg_s) };
            outgoing[up] += rate;
            let flux = rate * closure[up];
            corrections[up] -= flux;
            corrections[down] += flux;
            closure_bounds[up] += flux.abs();
            closure_bounds[down] += flux.abs();
        }
        for outlet in &self.outlets {
            outgoing[outlet.cell] += outlet.mass_flow_kg_s;
            let flux = outlet.mass_flow_kg_s * closure[outlet.cell];
            corrections[outlet.cell] -= flux;
            closure_bounds[outlet.cell] += flux.abs();
        }
        let mut residuals = Vec::with_capacity(count);
        let mut residual_bounds = Vec::with_capacity(count);
        let mut maximum_mean_mach: f64 = 0.0;
        for (index, cell) in self.cells.iter().enumerate() {
            let owner_residual_w = balance.net_power_w[&CellId(index as u64)];
            let bound = owner_residual_w.abs() + closure_bounds[index];
            if !bound.is_finite() || bound > controls.steady_power_tolerance_w {
                return Err(SteadyGasError::PublishedBalance { cell: index, residual_bound_w: bound });
            }
            let gas = states[index].state();
            let sonic_mass_rate = gas.density * cell.flow_area_m2 * gas.sound_speed;
            if !positive(sonic_mass_rate) { return Err(bad("sonic mass rate is not representable")); }
            let mach = outgoing[index] / sonic_mass_rate;
            if !positive(mach) || mach > controls.maximum_mean_mach {
                return Err(SteadyGasError::MeanFlow { cell: index, mach, limit: controls.maximum_mean_mach });
            }
            maximum_mean_mach = maximum_mean_mach.max(mach);
            residuals.push(owner_residual_w + corrections[index]);
            residual_bounds.push(bound);
        }
        let relaxation_time_s = thermal.time_s();
        // End borrowed transport state before moving its source laws.
        drop(thermal);
        laws.truncate(count);
        poll(&mut cancelled)?;
        Ok(SteadyGasNetworkSolution { states, laws, species: stationary,
            mass_flow_kg_s: outgoing, cell_power_residuals_w: residuals,
            cell_power_residual_bounds_w: residual_bounds, maximum_energy_residual_j,
            substeps, relaxation_time_s, maximum_mean_mach,
            inlet_mass_kg_s: inlet_mass, inlet_enthalpy_w: inlet_enthalpy })
    }
}

#[cfg(test)]
mod tests;
