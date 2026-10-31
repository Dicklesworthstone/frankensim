//! Steady axial equilibrium condensation through the original transport owner.
//!
//! Each well-mixed cell receives the previous cell's actual gas discharge and
//! any explicit injections. The shared constant-pressure water partition sets
//! its gas and liquid outlets at a trial temperature. `Advected` owns BOTH
//! their species balance and their enthalpy fluxes, and the original thermal
//! timeline supplies the wall power. A bounded scalar solve closes that same
//! owner residual; no second advection or caloric kernel is introduced.
//!
//! Liquid leaves each cell at its equilibrium temperature with aligned liquid
//! enthalpy. The returned liquid rates describe a drained wall-film reservoir,
//! not an accumulating steady liquid store or liquid advected downstream.
//! A one-second reference packet converts rates to masses for the homogeneous
//! phase partition; it is neither a physical time step nor stored gas volume.
//! There is no finite-rate condensation, film momentum, pressure drop or wall
//! thermal evolution. Wall conductances and temperatures are supplied data.
//! The original near-atmospheric liquid-water model limits remain active.

use super::{SteadyAxialGas, SteadyGasError, point_at, poll, positive};
use crate::MaterialError;
use crate::gas_field::{AxialGasField, PrimitiveGasState, resolve_primitive, typed};
use crate::gas_species::input::ThermochemicalPhase;
use crate::phase::condensation::{
    MAXIMUM_TEMPERATURE_K, MINIMUM_TEMPERATURE_K, WaterCondensationError, WaterCondensationModel,
};
use crate::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve, enthalpy::GasEnthalpyLaw};
use fs_matdb::QueryPoint;
use fs_qty::semantic::QuantityKind;
use fs_thermal::enthalpy::advection::species::{InletSpecies, StationarySpecies};
use fs_thermal::enthalpy::advection::{Advected, MassBoundary, SelectiveOutlet};
use fs_thermal::enthalpy::{
    Boundary, CellId, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation,
    Event, MaterialId, Reservoir, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};
use std::collections::BTreeMap;

/// Independent caloric, nonlinear-work, species and steady-power budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensingSteadyGasSolve {
    /// Original NASA inversion controls; its temperature interval brackets
    /// every cell root and must remain within the condensation model domain.
    pub caloric: GasEnergySolve,
    /// Complete temperature proposals per cell, including bracket endpoints.
    pub maximum_temperature_trials_per_cell: usize,
    /// Per-cell W budget, including caloric reconstruction and datum defects.
    pub steady_power_tolerance_w: f64,
    /// Absolute species/total-flow budget [kg/s], also applied to the full path.
    pub mass_flow_tolerance_kg_s: f64,
    /// Explicit low-Mach admission, positive and at most 0.1.
    pub maximum_mean_mach: f64,
}

/// Failed proposals publish no field, drained liquid or partial solution.
#[derive(Clone, Debug, PartialEq)]
pub enum CondensingSteadyGasError {
    Steady(SteadyGasError),
    Phase(WaterCondensationError),
    /// Net inward W at the low/high ends do not bracket a steady state.
    TemperatureBracket {
        cell: usize,
        residuals_w: [f64; 2],
    },
    NotConverged {
        cell: usize,
        trials: usize,
        residual_bound_w: f64,
    },
    SpeciesBalance {
        cell: usize,
        residual_kg_s: f64,
    },
    /// Original transport and independent source-relative caloric audit differ.
    CaloricBalance {
        cell: usize,
        residual_w: f64,
    },
}
impl core::fmt::Display for CondensingSteadyGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "steady condensing gas: {self:?}")
    }
}
impl std::error::Error for CondensingSteadyGasError {}
impl From<SteadyGasError> for CondensingSteadyGasError {
    fn from(error: SteadyGasError) -> Self {
        Self::Steady(error)
    }
}
impl From<MaterialError> for CondensingSteadyGasError {
    fn from(error: MaterialError) -> Self {
        Self::Steady(error.into())
    }
}
impl From<GasEnergyError> for CondensingSteadyGasError {
    fn from(error: GasEnergyError) -> Self {
        Self::Steady(error.into())
    }
}
impl From<DepositionError> for CondensingSteadyGasError {
    fn from(error: DepositionError) -> Self {
        Self::Steady(error.into())
    }
}
impl From<WaterCondensationError> for CondensingSteadyGasError {
    fn from(error: WaterCondensationError) -> Self {
        Self::Phase(error)
    }
}
fn invalid(message: &'static str) -> CondensingSteadyGasError {
    SteadyGasError::Invalid(message).into()
}
fn finite(value: f64) -> Result<f64, CondensingSteadyGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("steady condensation arithmetic is unrepresentable"))
    }
}

/// One accepted cell's original transport ledger, with drained liquid explicit.
#[derive(Clone, Debug, PartialEq)]
pub struct CondensingGasCell {
    pub temperature_k: f64,
    pub gas_mass_flow_kg_s: f64,
    pub condensed_mass_flow_kg_s: f64,
    /// Carried liquid enthalpy, including the declared common source datum [W].
    pub liquid_enthalpy_flow_w: f64,
    /// Diagnostic latent release, already included by the liquid outlet [W].
    pub latent_heat_w: f64,
    /// Wall reservoir heat into the gas/phase control volume [W].
    pub wall_heat_into_gas_w: f64,
    /// Original owner balance corrected to the emitted primitive temperature.
    pub power_residual_w: f64,
    /// Independent NASA sensible increments plus shared phase latent gap.
    pub source_relative_power_residual_w: f64,
    /// Bounds original and published balances without cancelling inverse errors.
    pub power_residual_bound_w: f64,
    pub maximum_species_mass_residual_kg_s: f64,
    pub temperature_trials: usize,
    /// Original stationary species result; this map contains this cell's ID.
    pub species: StationarySpecies,
}

/// A shared axial field plus all gas, liquid and energy fluxes that produced it.
#[derive(Clone, Debug, PartialEq)]
pub struct CondensingSteadyGasSolution {
    field: AxialGasField,
    pub cells: Vec<CondensingGasCell>,
    pub total_condensed_mass_flow_kg_s: f64,
    pub total_liquid_enthalpy_flow_w: f64,
    pub maximum_mean_mach: f64,
    pub maximum_species_mass_residual_kg_s: f64,
    /// Sum of independently checked cell balances, with internal streams cancelled.
    pub source_relative_power_residual_w: f64,
}
impl CondensingSteadyGasSolution {
    pub fn field(&self) -> &AxialGasField {
        &self.field
    }
    pub fn into_field(self) -> AxialGasField {
        self.field
    }
    pub fn maximum_power_residual_bound_w(&self) -> f64 {
        self.cells
            .iter()
            .map(|cell| cell.power_residual_bound_w)
            .fold(0.0, f64::max)
    }
}

#[derive(Clone)]
struct Stream {
    kg_s: f64,
    primitive: PrimitiveGasState,
    law: GasEnthalpyLaw,
}
impl Stream {
    fn new(
        problem: &SteadyAxialGas,
        point: &QueryPoint,
        kg_s: f64,
        primitive: PrimitiveGasState,
        controls: CondensingSteadyGasSolve,
    ) -> Result<Self, CondensingSteadyGasError> {
        resolve_primitive(
            &problem.source,
            point,
            ThermochemicalPhase::LiquidWaterIapwsLimit,
            &primitive,
        )?;
        let law = GasEnthalpyLaw::from_mole(
            &problem.source,
            primitive.composition.clone(),
            point_at(point, primitive.temperature_k)?,
            controls.caloric,
        )?;
        Ok(Self {
            kg_s,
            primitive,
            law,
        })
    }
    fn fractions(&self) -> &[f64] {
        self.law
            .advected_composition()
            .expect("gas law")
            .mass_fractions
    }
}

struct Candidate {
    stream: Stream,
    report: CondensingGasCell,
    mach: f64,
}

fn add_species(
    target: &mut BTreeMap<SpeciesId, f64>,
    source: &WaterCondensationModel,
    stream: &Stream,
) -> Result<(), CondensingSteadyGasError> {
    for (species, &fraction) in source
        .source()
        .thermodynamics()
        .components()
        .iter()
        .zip(stream.fractions())
    {
        let rate = finite(stream.kg_s * fraction)?;
        if fraction > 0.0 && rate == 0.0 {
            return Err(invalid("positive supplied species rate was lost"));
        }
        let stored = target.entry(species.species().clone()).or_default();
        let next = finite(*stored + rate)?;
        if rate > 0.0 && next == *stored {
            return Err(invalid("positive injected species rate was lost"));
        }
        *stored = next;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    problem: &SteadyAxialGas,
    model: &WaterCondensationModel,
    cell_index: usize,
    streams: &[Stream],
    incoming: &BTreeMap<SpeciesId, f64>,
    pressure_point: &QueryPoint,
    temperature_k: f64,
    controls: CondensingSteadyGasSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Candidate, CondensingSteadyGasError> {
    poll(cancelled)?;
    // Homogeneity permits one second's flux packet; no physical time is marched.
    let phase = model.partition_at_pressure_temperature(
        incoming,
        0.0,
        temperature_k,
        problem.static_pressure_pa,
        &problem.source_conditions,
        &mut *cancelled,
    )?;
    let id = CellId(cell_index as u64);
    let cell = &problem.cells[cell_index];
    let mut boundaries = Vec::with_capacity(streams.len() + 2);
    let mut inlet_species = Vec::with_capacity(streams.len());
    for (index, stream) in streams.iter().enumerate() {
        boundaries.push(MassBoundary::Inlet {
            id: index as u64,
            cell: id,
            kg_s: stream.kg_s,
            material_id: MaterialId((index + 1) as u64),
            specific_h_j_kg: stream
                .law
                .enthalpy_at_temperature(stream.primitive.temperature_k)?,
        });
        inlet_species.push(InletSpecies {
            boundary_id: index as u64,
            mass_fractions: stream.fractions().to_vec(),
        });
    }
    boundaries.push(MassBoundary::Outlet {
        id: streams.len() as u64,
        cell: id,
        kg_s: phase.gas_mass_kg,
    });
    let mut selective = Vec::new();
    if phase.liquid_mass_kg > 0.0 {
        let boundary_id = (streams.len() + 1) as u64;
        boundaries.push(MassBoundary::Outlet {
            id: boundary_id,
            cell: id,
            kg_s: phase.liquid_mass_kg,
        });
        let fractions = problem
            .source
            .thermodynamics()
            .components()
            .iter()
            .map(|species| {
                if species.species() == model.water_species() {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        selective.push(SelectiveOutlet {
            boundary_id,
            mass_fractions: fractions,
            specific_h_j_kg: phase.liquid_specific_enthalpy_j_kg,
            enthalpy_reference: streams[0]
                .law
                .advected_composition()
                .expect("gas law")
                .enthalpy_reference
                .to_owned(),
        });
    }
    let names = problem
        .source
        .thermodynamics()
        .components()
        .iter()
        .map(|species| species.species().to_string())
        .collect();
    let observer = Advected::new((), vec![], super::network::MAX_STEADY_GAS_CONNECTIONS)?
        .with_boundaries(boundaries, super::network::MAX_STEADY_GAS_CONNECTIONS)?
        .with_selective_outlets(selective, 128)?
        .with_stationary_species(names, inlet_species, 1, 128)?;
    let species = observer
        .stationary_species()
        .expect("just admitted")
        .clone();
    let fractions = &species.cells()[&id];
    let law = GasEnthalpyLaw::new(
        &problem.source,
        Composition::new(CompositionBasis::MassFraction, fractions.clone())
            .map_err(|_| invalid("original stationary composition refused"))?,
        point_at(pressure_point, temperature_k)?,
        controls.caloric,
    )?;
    let primitive = PrimitiveGasState {
        temperature_k,
        composition: law.mole_composition().clone(),
    };
    let gas = resolve_primitive(
        &problem.source,
        pressure_point,
        ThermochemicalPhase::LiquidWaterIapwsLimit,
        &primitive,
    )?;
    let mut mass_residual = phase
        .mass_residual_kg
        .abs()
        .max(phase.maximum_species_mass_residual_kg)
        .max(species.maximum_residual_kg_s());
    for (source, &fraction) in problem
        .source
        .thermodynamics()
        .components()
        .iter()
        .zip(fractions)
    {
        mass_residual = mass_residual.max(
            (phase.gas_mass_kg * fraction - phase.gas_species_masses_kg[source.species()]).abs(),
        );
    }
    if !mass_residual.is_finite() || mass_residual > controls.mass_flow_tolerance_kg_s {
        return Err(CondensingSteadyGasError::SpeciesBalance {
            cell: cell_index,
            residual_kg_s: mass_residual,
        });
    }
    let start_m = if cell_index == 0 {
        0.0
    } else {
        problem.cells[cell_index - 1].end_m
    };
    let geometry = fs_thermal::enthalpy::CellGeometry {
        start_m: [start_m, 0.0, 0.0],
        end_m: [cell.end_m, 0.0, 0.0],
        volume_m3: cell.volume_m3,
        surface_area_m2: cell.wall_area_m2,
    };
    let mass_kg = gas.state().density * cell.volume_m3;
    if !positive(mass_kg) {
        return Err(invalid("cell reference gas mass is unrepresentable"));
    }
    let birth = TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id,
            segment_id: cell_index as u64,
            body_id: 0,
            material_id: MaterialId(0),
            mass_kg,
            specific_h_j_kg: law.enthalpy_at_temperature(temperature_k)?,
            geometry,
            boundary: Boundary {
                ambient_k: cell.wall_temperature_k,
                htc_w_m2_k: 0.0,
                reservoir: Some(Reservoir {
                    area_m2: cell.wall_area_m2,
                    conductance_w_k: cell.wall_conductance_w_k,
                    temperature_k: cell.wall_temperature_k,
                }),
                source_w: 0.0,
            },
        }),
    };
    let mut laws: Vec<(MaterialId, &dyn ThermalLaw)> = vec![(MaterialId(0), &law)];
    for (index, stream) in streams.iter().enumerate() {
        laws.push((MaterialId((index + 1) as u64), &stream.law));
    }
    let budget = DepositionBudget {
        max_cells: laws.len(),
        max_links: 0,
        max_events: 1,
        max_substeps: 1,
    };
    let (thermal, _) = DepositionSimulation::new(vec![birth], &laws, budget)?.advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 0.0,
        },
        || poll(cancelled),
    )?;
    let balance = thermal.power_balance_observed(&observer, || poll(cancelled))?;
    let owner_residual = balance.net_power_w[&id];
    let correction = cell.wall_conductance_w_k * (balance.temperatures_k[&id] - temperature_k);
    let power_residual = finite(owner_residual + correction)?;
    let wall_power = finite(cell.wall_conductance_w_k * (cell.wall_temperature_k - temperature_k))?;
    // Audit only: the original observer above owns every actual flux. Stable
    // species sensible increments cancel formation datums before subtraction.
    let mut sensible_w = 0.0;
    for stream in streams {
        for (species, &fraction) in problem
            .source
            .thermodynamics()
            .components()
            .iter()
            .zip(stream.fractions())
        {
            if fraction != 0.0 {
                sensible_w = finite(
                    sensible_w
                        + stream.kg_s
                            * fraction
                            * problem.source.species_specific_enthalpy_change_j_kg(
                                species.species(),
                                stream.primitive.temperature_k,
                                temperature_k,
                            )?,
                )?;
            }
        }
    }
    let latent_w = finite(phase.condensed_mass_kg * phase.latent_enthalpy_j_kg)?;
    let source_relative = finite(wall_power - sensible_w + latent_w)?;
    let datum_defect = finite(power_residual - source_relative)?;
    if datum_defect.abs() > controls.steady_power_tolerance_w {
        return Err(CondensingSteadyGasError::CaloricBalance {
            cell: cell_index,
            residual_w: datum_defect,
        });
    }
    let residual_bound = finite(owner_residual.abs() + correction.abs() + datum_defect.abs())?;
    let flow_area = cell.volume_m3 / (cell.end_m - start_m);
    let sonic_mass_rate = gas.state().density * flow_area * gas.state().sound_speed;
    if !positive(sonic_mass_rate) {
        return Err(invalid("sonic mass rate is unrepresentable"));
    }
    let mach = phase.gas_mass_kg / sonic_mass_rate;
    if !positive(mach) {
        return Err(invalid("mean gas flow is unrepresentable"));
    }
    // End original borrowed thermal state before moving the original gas law.
    drop(thermal);
    drop(laws);
    poll(cancelled)?;
    Ok(Candidate {
        stream: Stream {
            kg_s: phase.gas_mass_kg,
            primitive,
            law,
        },
        mach,
        report: CondensingGasCell {
            temperature_k,
            gas_mass_flow_kg_s: phase.gas_mass_kg,
            condensed_mass_flow_kg_s: phase.condensed_mass_kg,
            liquid_enthalpy_flow_w: finite(
                phase.liquid_mass_kg * phase.liquid_specific_enthalpy_j_kg,
            )?,
            latent_heat_w: latent_w,
            wall_heat_into_gas_w: wall_power,
            power_residual_w: power_residual,
            source_relative_power_residual_w: source_relative,
            power_residual_bound_w: residual_bound,
            maximum_species_mass_residual_kg_s: mass_residual,
            temperature_trials: 0,
            species,
        },
    })
}

fn solve_cell(
    cell: usize,
    controls: CondensingSteadyGasSolve,
    mut trial: impl FnMut(f64) -> Result<Candidate, CondensingSteadyGasError>,
) -> Result<Candidate, CondensingSteadyGasError> {
    let mut lo = controls.caloric.minimum_temperature_k;
    let mut hi = controls.caloric.maximum_temperature_k;
    let mut lower = trial(lo)?;
    lower.report.temperature_trials = 1;
    if lower.report.power_residual_bound_w <= controls.steady_power_tolerance_w {
        return Ok(lower);
    }
    let mut upper = trial(hi)?;
    upper.report.temperature_trials = 2;
    if upper.report.power_residual_bound_w <= controls.steady_power_tolerance_w {
        return Ok(upper);
    }
    // Net inward power decreases with temperature on the admitted sensible/
    // latent branch; solve its negative without inventing a thermal tangent.
    let mut f_lo = -lower.report.power_residual_w;
    let mut f_hi = -upper.report.power_residual_w;
    if f_lo > 0.0 || f_hi < 0.0 {
        return Err(CondensingSteadyGasError::TemperatureBracket {
            cell,
            residuals_w: [-f_lo, -f_hi],
        });
    }
    let mut last_bound = lower
        .report
        .power_residual_bound_w
        .min(upper.report.power_residual_bound_w);
    let mut trials = 2;
    while trials < controls.maximum_temperature_trials_per_cell {
        let fraction = -f_lo / (f_hi - f_lo);
        let weight = if (0.1..=0.9).contains(&fraction) {
            fraction
        } else {
            0.5
        };
        let t = lo + weight * (hi - lo);
        if t <= lo || t >= hi {
            break;
        }
        let mut candidate = trial(t)?;
        trials += 1;
        candidate.report.temperature_trials = trials;
        last_bound = candidate.report.power_residual_bound_w;
        if last_bound <= controls.steady_power_tolerance_w {
            return Ok(candidate);
        }
        let residual = -candidate.report.power_residual_w;
        if residual < 0.0 {
            lo = t;
            f_lo = residual;
        } else {
            hi = t;
            f_hi = residual;
        }
    }
    Err(CondensingSteadyGasError::NotConverged {
        cell,
        trials,
        residual_bound_w: last_bound,
    })
}

impl SteadyAxialGas {
    /// Generate a condensation-capable field through original thermal/species
    /// owners. Requires the explicit IAPWS phase policy and exactly this model's
    /// source catalog. The same common-pressure, positive-flow and geometry
    /// admission as `to_network` applies. Existing dry solves are unchanged.
    ///
    /// Field x=0 is the authored upstream inlet; side supplies enter their
    /// receiving finite volume. Each emitted interval is its well-mixed cell,
    /// cached in `AxialGasField` with explicit piecewise-constant sampling.
    /// No partial field or film flux is returned on any refusal/cancellation.
    pub fn solve_condensing(
        &self,
        model: &WaterCondensationModel,
        controls: CondensingSteadyGasSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CondensingSteadyGasSolution, CondensingSteadyGasError> {
        poll(&mut cancelled)?;
        self.to_network()?;
        controls.caloric.validate()?;
        if self.source != *model.source()
            || self.phase != ThermochemicalPhase::LiquidWaterIapwsLimit
            || !(MINIMUM_TEMPERATURE_K..=MAXIMUM_TEMPERATURE_K)
                .contains(&controls.caloric.minimum_temperature_k)
            || !(MINIMUM_TEMPERATURE_K..=MAXIMUM_TEMPERATURE_K)
                .contains(&controls.caloric.maximum_temperature_k)
            || !(2..=128).contains(&controls.maximum_temperature_trials_per_cell)
            || !positive(controls.steady_power_tolerance_w)
            || !positive(controls.mass_flow_tolerance_kg_s)
            || !positive(controls.maximum_mean_mach)
            || controls.maximum_mean_mach > 0.1
        {
            return Err(invalid(
                "condensing solve needs the same source, IAPWS phase, bounded temperature/work and positive balance budgets",
            ));
        }
        let pressure_point = self
            .source_conditions
            .clone()
            .with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                self.static_pressure_pa,
            )
            .map_err(|_| invalid("invalid common pressure coordinate"))?;
        let inlet = Stream::new(
            self,
            &pressure_point,
            self.inlet_mass_flow_kg_s,
            self.inlet.clone(),
            controls,
        )?;
        let mut supplies = vec![Vec::<Stream>::new(); self.cells.len()];
        let mut total_supplied_species = BTreeMap::new();
        add_species(&mut total_supplied_species, model, &inlet)?;
        for injection in &self.injections {
            poll(&mut cancelled)?;
            let stream = Stream::new(
                self,
                &pressure_point,
                injection.mass_flow_kg_s,
                injection.state.clone(),
                controls,
            )?;
            add_species(&mut total_supplied_species, model, &stream)?;
            supplies[injection.cell].push(stream);
        }
        let inlet_gas = resolve_primitive(&self.source, &pressure_point, self.phase, &self.inlet)?;
        let inlet_area = self.cells[0].volume_m3 / self.cells[0].end_m;
        let mut maximum_mach = self.inlet_mass_flow_kg_s
            / (inlet_gas.state().density * inlet_area * inlet_gas.state().sound_speed);
        if !positive(maximum_mach) || maximum_mach > controls.maximum_mean_mach {
            return Err(SteadyGasError::MeanFlow {
                cell: 0,
                mach: maximum_mach,
                limit: controls.maximum_mean_mach,
            }
            .into());
        }
        let mut upstream = inlet;
        let mut reports = Vec::with_capacity(self.cells.len());
        let mut stations = vec![0.0];
        let mut primitives = vec![self.inlet.clone()];
        let mut liquid_rate = 0.0;
        let mut liquid_power = 0.0;
        let mut maximum_species_residual: f64 = 0.0;
        let mut power_residual = 0.0;
        for (index, extra) in supplies.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let mut streams = Vec::with_capacity(extra.len() + 1);
            streams.push(upstream);
            streams.extend(extra);
            let mut incoming = BTreeMap::new();
            for stream in &streams {
                add_species(&mut incoming, model, stream)?;
            }
            let candidate = solve_cell(index, controls, |temperature_k| {
                evaluate(
                    self,
                    model,
                    index,
                    &streams,
                    &incoming,
                    &pressure_point,
                    temperature_k,
                    controls,
                    &mut cancelled,
                )
            })?;
            if candidate.mach > controls.maximum_mean_mach {
                return Err(SteadyGasError::MeanFlow {
                    cell: index,
                    mach: candidate.mach,
                    limit: controls.maximum_mean_mach,
                }
                .into());
            }
            maximum_mach = maximum_mach.max(candidate.mach);
            maximum_species_residual =
                maximum_species_residual.max(candidate.report.maximum_species_mass_residual_kg_s);
            liquid_rate = finite(liquid_rate + candidate.report.condensed_mass_flow_kg_s)?;
            liquid_power = finite(liquid_power + candidate.report.liquid_enthalpy_flow_w)?;
            power_residual =
                finite(power_residual + candidate.report.source_relative_power_residual_w)?;
            stations.push(self.cells[index].end_m);
            primitives.push(candidate.stream.primitive.clone());
            reports.push(candidate.report);
            upstream = candidate.stream;
        }
        let mut discharged = BTreeMap::new();
        add_species(&mut discharged, model, &upstream)?;
        *discharged.entry(model.water_species().clone()).or_default() += liquid_rate;
        for (species, incoming) in &total_supplied_species {
            maximum_species_residual =
                maximum_species_residual.max((discharged[species] - incoming).abs());
        }
        if !maximum_species_residual.is_finite()
            || maximum_species_residual > controls.mass_flow_tolerance_kg_s
        {
            return Err(CondensingSteadyGasError::SpeciesBalance {
                cell: self.cells.len() - 1,
                residual_kg_s: maximum_species_residual,
            });
        }
        let field = AxialGasField::new_piecewise_constant(
            stations,
            self.static_pressure_pa,
            primitives,
            self.source.clone(),
            self.source_conditions.clone(),
            self.phase,
        )?;
        poll(&mut cancelled)?;
        Ok(CondensingSteadyGasSolution {
            field,
            cells: reports,
            total_condensed_mass_flow_kg_s: liquid_rate,
            total_liquid_enthalpy_flow_w: liquid_power,
            maximum_mean_mach: maximum_mach,
            maximum_species_mass_residual_kg_s: maximum_species_residual,
            source_relative_power_residual_w: power_residual,
        })
    }
}

#[cfg(test)]
mod tests;
