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

use super::{AxialGasField, PrimitiveGasState, resolve_primitive, typed};
use crate::{MaterialError, gas_species::input::ThermochemicalPhase};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::{GasEnergyError, GasEnergySolve}};
use crate::thermochemical_gas::energy::enthalpy::GasEnthalpyLaw;
use fs_matdb::QueryPoint;
use fs_qty::semantic::QuantityKind;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::enthalpy::advection::{Advected, MassBoundary, MassFlow};
use fs_thermal::enthalpy::advection::species::{InletSpecies, StationarySpecies};
use fs_thermochem::{Composition, CompositionBasis};

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

/// An additional gas supply mixed into one axial control volume.
#[derive(Clone, Debug, PartialEq)]
pub struct GasInjection {
    /// Zero-based receiving cell. Downstream through-flow increases by this rate.
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
    /// Solve conservative stationary composition first, then relax the coupled
    /// species-specific enthalpy/wall balance with the original thermal owner.
    /// Only a converged, source-admitted, fixed-pressure field is returned.
    ///
    /// The grid axis is first-order donor-cell transport: each downstream face
    /// takes its cell's mixed state. The field interpolates primitive values
    /// between those face samples. Refine this grid independently of the duct
    /// sampler and thermal pseudo-time step; a residual tolerance is not a
    /// spatial-accuracy certificate.
    pub fn solve(
        &self, controls: SteadyGasSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<SteadyGasSolution, SteadyGasError> {
        poll(&mut cancelled)?;
        if self.cells.is_empty() || self.cells.len() > MAX_STEADY_GAS_CELLS
            || self.injections.len() > MAX_STEADY_GAS_INJECTIONS
            || !positive(self.inlet_mass_flow_kg_s) || !positive(self.static_pressure_pa)
            || !positive(controls.maximum_dt_s) || !positive(controls.steady_power_tolerance_w)
            || !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j < 0.0
            || controls.maximum_relaxations == 0 || controls.maximum_relaxations > 1_000_000
            || controls.maximum_substeps == 0 || controls.maximum_substeps > 1_000_000
            || !positive(controls.maximum_mean_mach) || controls.maximum_mean_mach > 0.1
        {
            return Err(bad("bounded cells/injections/work, positive flow/pressure/tolerances and a low-Mach ceiling required"));
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
        let pressure_point = self.source_conditions.clone()
            .with_quantity("pressure", typed(QuantityKind::Pressure), self.static_pressure_pa)
            .map_err(|_| bad("invalid common pressure coordinate"))?;
        let count = self.cells.len();
        let inlet_count = self.injections.len() + 1;
        let mut inlet_data = vec![(0, self.inlet_mass_flow_kg_s, &self.inlet)];
        inlet_data.extend(self.injections.iter()
            .map(|inlet| (inlet.cell, inlet.mass_flow_kg_s, &inlet.state)));
        let mut inlet_laws = Vec::with_capacity(inlet_count);
        let mut boundaries = Vec::with_capacity(inlet_count + 1);
        let mut inlet_species = Vec::with_capacity(inlet_count);
        let mut incoming_mass = vec![0.0; count];
        let mut upstream_enthalpy_rate = 0.0;
        for (index, &(cell, kg_s, primitive)) in inlet_data.iter().enumerate() {
            poll(&mut cancelled)?;
            // Source/phase admission precedes any mixing or energy conversion.
            resolve_primitive(&self.source, &pressure_point, self.phase, primitive)?;
            let law = GasEnthalpyLaw::from_mole(&self.source, primitive.composition.clone(),
                point_at(&pressure_point, primitive.temperature_k)?, controls.caloric)?;
            let h = law.enthalpy_at_temperature(primitive.temperature_k)?;
            let id = index as u64;
            boundaries.push(MassBoundary::Inlet { id, cell: CellId(cell as u64), kg_s,
                material_id: MaterialId((count + index) as u64), specific_h_j_kg: h });
            inlet_species.push(InletSpecies { boundary_id: id,
                mass_fractions: law.advected_composition().expect("gas composition").mass_fractions.to_vec() });
            incoming_mass[cell] += kg_s;
            if cell == 0 { upstream_enthalpy_rate += kg_s * h; }
            inlet_laws.push(law);
        }
        let mut flows = Vec::with_capacity(count - 1);
        let mut mass_flow = Vec::with_capacity(count);
        let mut rate = 0.0;
        for (index, &incoming) in incoming_mass.iter().enumerate() {
            rate += incoming;
            if !positive(rate) { return Err(bad("axial mass flow is not representable")); }
            mass_flow.push(rate);
            if index + 1 < count {
                flows.push(MassFlow { id: (inlet_count + index) as u64,
                    a: CellId(index as u64), b: CellId((index + 1) as u64), kg_s: rate });
            }
        }
        if !upstream_enthalpy_rate.is_finite() {
            return Err(bad("upstream enthalpy rate is not representable"));
        }
        boundaries.push(MassBoundary::Outlet { id: (inlet_count + count) as u64,
            cell: CellId((count - 1) as u64), kg_s: rate });
        let species = self.source.thermodynamics().components().iter()
            .map(|component| component.species().to_string()).collect();
        let observer = Advected::new((), flows, MAX_STEADY_GAS_CELLS)?
            .with_boundaries(boundaries, MAX_STEADY_GAS_CELLS + MAX_STEADY_GAS_INJECTIONS + 2)?
            .with_stationary_species(species, inlet_species, count, 128)?;
        let stationary = observer.stationary_species().expect("just admitted").clone();
        let mut laws = Vec::with_capacity(count + inlet_count);
        for index in 0..count {
            poll(&mut cancelled)?;
            let fractions = &stationary.cells()[&CellId(index as u64)];
            // Preserve the authored mole coefficients on uniform composition runs.
            let same = inlet_laws.iter().find(|law|
                law.advected_composition().expect("gas law").mass_fractions == fractions.as_slice());
            let law = if let Some(law) = same { law.clone() } else {
                GasEnthalpyLaw::new(&self.source,
                    Composition::new(CompositionBasis::MassFraction, fractions.clone())
                        .map_err(|_| bad("stationary composition refused"))?,
                    point_at(&pressure_point, self.cells[index].initial_temperature_k)?, controls.caloric)?
            };
            laws.push(law);
        }
        laws.extend(inlet_laws);
        let mut schedule = Vec::with_capacity(count);
        previous = 0.0;
        for (index, cell) in self.cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let gas = laws[index].model().evaluate(&point_at(&pressure_point, cell.initial_temperature_k)?)?;
            let mass_kg = gas.state().density * cell.volume_m3;
            if !positive(mass_kg) { return Err(bad("initial EOS relaxation mass is not representable")); }
            schedule.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(index as u64), segment_id: index as u64, body_id: 0,
                material_id: MaterialId(index as u64), mass_kg,
                specific_h_j_kg: laws[index].enthalpy_at_temperature(cell.initial_temperature_k)?,
                geometry: CellGeometry { start_m: [previous, 0.0, 0.0], end_m: [cell.end_m, 0.0, 0.0],
                    volume_m3: cell.volume_m3, surface_area_m2: cell.wall_area_m2 },
                boundary: Boundary { ambient_k: cell.wall_temperature_k, htc_w_m2_k: 0.0,
                    reservoir: Some(Reservoir { area_m2: cell.wall_area_m2,
                        conductance_w_k: cell.wall_conductance_w_k,
                        temperature_k: cell.wall_temperature_k }), source_w: 0.0 },
            }) });
            previous = cell.end_m;
        }
        let references: Vec<_> = laws.iter().enumerate()
            .map(|(index, law)| (MaterialId(index as u64), law as &dyn ThermalLaw)).collect();
        let budget = DepositionBudget { max_cells: laws.len(), max_links: count,
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
        // The upstream boundary is an enthalpy-balanced mixture of all supplies
        // entering cell zero, before any of its wall heat is exchanged.
        let upstream_h = upstream_enthalpy_rate / incoming_mass[0];
        let upstream = laws[0].resolve_enthalpy(upstream_h, &mut cancelled)?;
        let mut stations = vec![0.0];
        let mut states = vec![PrimitiveGasState {
            temperature_k: upstream.state().state().temperature,
            composition: laws[0].mole_composition().clone(),
        }];
        let mut residuals = Vec::with_capacity(count);
        let mut maximum_power_residual_bound_w: f64 = 0.0;
        let mut incoming_closure_power_w = incoming_mass[0] * upstream.residual_j_kg();
        if !incoming_closure_power_w.is_finite() {
            return Err(bad("upstream caloric flux correction is not representable"));
        }
        let mut maximum_mean_mach: f64 = 0.0;
        previous = 0.0;
        for (index, cell) in self.cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let h = thermal.cells()[&CellId(index as u64)].specific_h_j_kg;
            let resolved = laws[index].resolve_enthalpy(h, &mut cancelled)?;
            let gas = resolved.state().state();
            // A loose J/kg inverse tolerance must never bypass the independent
            // W steady tolerance. Bound each actual closure flux separately:
            // equal inlet/outlet errors must not hide loss of supplied enthalpy.
            let outgoing_closure_power_w = mass_flow[index] * resolved.residual_j_kg();
            let owner_residual_w = balance.net_power_w[&CellId(index as u64)];
            let residual_bound_w = owner_residual_w.abs()
                + incoming_closure_power_w.abs() + outgoing_closure_power_w.abs();
            if !residual_bound_w.is_finite() {
                return Err(bad("published caloric power balance is not representable"));
            }
            if residual_bound_w > controls.steady_power_tolerance_w {
                return Err(SteadyGasError::PublishedBalance { cell: index, residual_bound_w });
            }
            maximum_power_residual_bound_w = maximum_power_residual_bound_w.max(residual_bound_w);
            let area = cell.volume_m3 / (cell.end_m - previous);
            let sonic_mass_rate = gas.density * area * gas.sound_speed;
            if !positive(area) || !positive(sonic_mass_rate) {
                return Err(bad("mean-flow area or sonic mass rate is not representable"));
            }
            let mach = mass_flow[index] / sonic_mass_rate;
            if !positive(mach) || mach > controls.maximum_mean_mach {
                return Err(SteadyGasError::MeanFlow { cell: index, mach, limit: controls.maximum_mean_mach });
            }
            maximum_mean_mach = maximum_mean_mach.max(mach);
            if index == 0 {
                let gas = upstream.state().state();
                let sonic_mass_rate = gas.density * area * gas.sound_speed;
                if !positive(sonic_mass_rate) {
                    return Err(bad("inlet sonic mass rate is not representable"));
                }
                let mach = incoming_mass[0] / sonic_mass_rate;
                if !positive(mach) || mach > controls.maximum_mean_mach {
                    return Err(SteadyGasError::MeanFlow { cell: index, mach, limit: controls.maximum_mean_mach });
                }
                maximum_mean_mach = maximum_mean_mach.max(mach);
            }
            stations.push(cell.end_m);
            states.push(PrimitiveGasState { temperature_k: gas.temperature,
                composition: laws[index].mole_composition().clone() });
            residuals.push(owner_residual_w + incoming_closure_power_w - outgoing_closure_power_w);
            incoming_closure_power_w = outgoing_closure_power_w;
            previous = cell.end_m;
        }
        poll(&mut cancelled)?;
        let field = AxialGasField::new(stations, self.static_pressure_pa, states,
            self.source.clone(), self.source_conditions.clone(), self.phase)?;
        poll(&mut cancelled)?;
        Ok(SteadyGasSolution { field, species: stationary, mass_flow_kg_s: mass_flow,
            cell_power_residuals_w: residuals, maximum_power_residual_bound_w,
            maximum_energy_residual_j, substeps,
            relaxation_time_s: thermal.time_s(), maximum_mean_mach })
    }
}

#[cfg(test)]
pub(crate) mod tests;
