//! Pressure-selected gas branches between common imposed static-pressure ports.
//!
//! Each branch reuses the original thermal/momentum/gravity solve. The ideal
//! zero-storage outlet mixes ORIGINAL source enthalpy at the common composition
//! and EOS pressure. Unequal densities change the flow split and local volume
//! rates; neither temperature nor volume flow is averaged to conserve energy.
//! This is not a general pressure graph or a manifold momentum/junction model.

/// Native density-aware fan selection of this same thermal flow split.
pub mod fan;

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_material::gas_field::steady::{MAX_STEADY_GAS_CELLS, SteadyAxialGas, SteadyGasSolve};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::{GasEnergyError, enthalpy::GasEnthalpyLaw};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec};

use super::super::{GasCellTransport, GasTransportFeedback, positive};
use super::{
    GasGravityBudget, GasPressureError, PressureDrivenGasSolution, SeriesGravity,
    solve_pressure_driven_material_gas, solve_pressure_gravity_driven_material_gas,
};

/// Bound the complete retained set of original branch solutions.
pub const MAX_PARALLEL_GAS_BRANCHES: usize = 64;

/// One existing axial source/geometry description and its supported flow controls.
#[derive(Clone, Copy, Debug)]
pub struct ParallelGasBranch<'a> {
    pub name: &'a str,
    pub problem: &'a SteadyAxialGas,
    pub transport: &'a [GasCellTransport],
    /// Original pressure/profile controls and branch-local potential allowance.
    /// Without gravity the potential allowance must be zero.
    pub budget: GasGravityBudget,
    /// Either every branch supplies a path or none does. Paths share uniform
    /// acceleration and exact inlet/outlet vertices; their interiors may differ.
    pub gravity: Option<&'a SeriesGravity>,
}

/// Aggregate bounds, shared across every branch and all its native trials.
#[derive(Clone, Copy, Debug)]
pub struct ParallelGasBudget {
    pub maximum_flow_evaluations: usize,
    pub maximum_substeps: usize,
    pub maximum_neglected_hydraulic_power_w: f64,
    /// Branch regional changes PLUS inlet/outlet kinetic-energy flux magnitudes.
    /// The latter explicitly bound energy ignored by the ideal split/merge.
    /// They are not a hydraulic junction-loss coefficient or pressure closure.
    pub maximum_neglected_kinetic_power_w: f64,
    pub maximum_neglected_potential_power_w: f64,
}

impl ParallelGasBudget {
    // Shared admission precedes any outer fan work clipping.
    fn admit(self) -> Result<(), ParallelGasError> {
        let nonnegative = |x: f64| x.is_finite() && x >= 0.0;
        if !(2..=128 * MAX_PARALLEL_GAS_BRANCHES).contains(&self.maximum_flow_evaluations)
            || self.maximum_substeps == 0
            || !nonnegative(self.maximum_neglected_hydraulic_power_w)
            || !nonnegative(self.maximum_neglected_kinetic_power_w)
            || !nonnegative(self.maximum_neglected_potential_power_w)
        {
            return Err(invalid(
                "positive shared work and finite omission allowances required",
            ));
        }
        Ok(())
    }
}

/// A refusal publishes no branch collection and never mutates a problem.
#[derive(Debug)]
pub enum ParallelGasError {
    Invalid(&'static str),
    Branch {
        index: usize,
        name: String,
        source: GasPressureError,
    },
    Caloric(GasEnergyError),
    WorkBudget {
        flow_evaluations: usize,
        substeps: usize,
    },
    OmittedPower {
        hydraulic_w: f64,
        kinetic_w: f64,
        potential_w: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ParallelGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "parallel source gas: {self:?}")
    }
}
impl std::error::Error for ParallelGasError {}
impl From<GasEnergyError> for ParallelGasError {
    fn from(error: GasEnergyError) -> Self {
        Self::Caloric(error)
    }
}

/// Original named branch and its complete accepted native result.
#[derive(Debug)]
pub struct ParallelGasBranchState {
    pub name: String,
    pub solution: PressureDrivenGasSolution,
}

/// Density-responsive flow split and source-caloric adiabatic outlet mixture.
#[derive(Debug)]
pub struct ParallelGasSolution {
    pub branches: Vec<ParallelGasBranchState>,
    pub driving_pressure_pa: f64,
    pub inlet: ResolvedThermochemicalGasState,
    pub mixed_outlet: ResolvedThermochemicalGasState,
    pub total_mass_flow_kg_per_s: f64,
    pub inlet_volume_flow_m3_per_s: f64,
    /// Total volume flow AFTER mixing, from actual mixed source density.
    pub outlet_volume_flow_m3_per_s: f64,
    pub total_wall_heat_w: f64,
    /// Sum mdot * original-law Delta h relative to the common inlet [W].
    pub total_enthalpy_gain_w: f64,
    /// Mixed original enthalpy flux minus the sum of branch outlet fluxes [W].
    pub mixing_enthalpy_residual_w: f64,
    pub neglected_hydraulic_power_w: f64,
    pub neglected_kinetic_power_bound_w: f64,
    pub neglected_potential_power_bound_w: f64,
    /// Included once in the aggregate kinetic bound, never deposited as heat.
    pub manifold_kinetic_power_bound_w: f64,
    pub flow_evaluations: usize,
    pub substeps: usize,
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ParallelGasError> {
    cx.checkpoint().map_err(|_| ParallelGasError::Cancelled)
}
fn invalid(reason: &'static str) -> ParallelGasError {
    ParallelGasError::Invalid(reason)
}

/// Select positive branch mass flows at one imposed static-pressure difference.
///
/// Zero branches refuse; one branch is supported. All branches must have the
/// identical source, external conditions, phase, common pressure and inlet.
/// Side injections are excluded by the native pressure owner. Each area stays
/// constant along its branch; different branches may have different areas.
/// Supplied gravity paths must connect the same two physical port locations.
///
/// The original branch limits, source/correlation support and pressure brackets
/// remain authoritative. No monotonicity, uniqueness or stability is inferred.
/// Each call receives only remaining shared work; a failing branch stops the
/// whole calculation, without retry or a second publication solve.
///
/// An ideal adiabatic zero-storage manifold mixes specific enthalpy. Kinetic
/// energy entering/leaving the two manifolds is explicitly neglected and bounded;
/// it is NOT silently converted into heat or a junction pressure loss. Authored
/// branch resistances own any such losses. This is not a compressible pressure
/// field, recirculating inventory, or complete manifold momentum/energy model.
pub fn solve_parallel_pressure_driven_material_gas(
    cx: &Cx<'_>,
    branches: &[ParallelGasBranch<'_>],
    drive: Pressure,
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
    budget: ParallelGasBudget,
) -> Result<ParallelGasSolution, ParallelGasError> {
    checkpoint(cx)?;
    budget.admit()?;
    if branches.is_empty() || branches.len() > MAX_PARALLEL_GAS_BRANCHES {
        return Err(invalid("bounded nonempty parallel branches required"));
    }
    let first = branches[0].problem;
    let mut names = BTreeSet::new();
    let mut count = 0usize;
    for branch in branches {
        checkpoint(cx)?;
        let p = branch.problem;
        if branch.name.trim().is_empty()
            || !names.insert(branch.name)
            || p.source != first.source
            || p.source_conditions != first.source_conditions
            || p.static_pressure_pa != first.static_pressure_pa
            || p.phase != first.phase
            || p.inlet != first.inlet
        {
            return Err(invalid(
                "parallel ports require unique names and identical source, pressure, phase and inlet",
            ));
        }
        if !(2..=128).contains(&branch.budget.flow.pressure.max_flow_evaluations) {
            return Err(invalid(
                "original branch pressure solve budget must be in 2..=128",
            ));
        }
        count = count.saturating_add(p.cells.len());
        if count > MAX_STEADY_GAS_CELLS
            || (branch.gravity.is_none()
                && branch.budget.maximum_neglected_potential_power_w != 0.0)
        {
            return Err(invalid(
                "whole-network cell cap or gravity omission ownership violated",
            ));
        }
        match (branches[0].gravity, branch.gravity) {
            (None, None) => {}
            (Some(a), Some(b))
                if a.gravity_m_s2() == b.gravity_m_s2()
                    && a.vertices_m().first() == b.vertices_m().first()
                    && a.vertices_m().last() == b.vertices_m().last() => {}
            _ => {
                return Err(invalid(
                    "gravity branches require the same acceleration and physical inlet/outlet",
                ));
            }
        }
    }
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let context = first
        .source_conditions
        .clone()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            first.inlet.temperature_k,
        )
        .and_then(|p| {
            p.with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                first.static_pressure_pa,
            )
        })
        .map_err(|_| invalid("invalid common source coordinates"))?;
    let law = GasEnthalpyLaw::from_mole(
        &first.source,
        first.inlet.composition.clone(),
        context,
        controls.caloric,
    )?;
    let mut results = Vec::with_capacity(branches.len());
    let (mut evaluations, mut substeps) = (0, 0);
    let (mut mass, mut heat, mut gain) = (0.0, 0.0, 0.0);
    let (mut hydraulic, mut kinetic, mut potential, mut manifold) = (0.0, 0.0, 0.0, 0.0);
    let mut inlet: Option<ResolvedThermochemicalGasState> = None;
    for (index, branch) in branches.iter().enumerate() {
        checkpoint(cx)?;
        let remaining = budget.maximum_flow_evaluations.saturating_sub(evaluations);
        if remaining < 2 || substeps >= budget.maximum_substeps {
            return Err(ParallelGasError::WorkBudget {
                flow_evaluations: evaluations,
                substeps,
            });
        }
        let mut local = branch.budget;
        local.flow.pressure.max_flow_evaluations =
            local.flow.pressure.max_flow_evaluations.min(remaining);
        local.flow.maximum_substeps = local
            .flow
            .maximum_substeps
            .min(budget.maximum_substeps - substeps);
        let solved = if let Some(path) = branch.gravity {
            solve_pressure_gravity_driven_material_gas(
                cx,
                branch.problem,
                branch.transport,
                drive,
                path,
                controls,
                feedback,
                local,
            )
        } else {
            solve_pressure_driven_material_gas(
                cx,
                branch.problem,
                branch.transport,
                drive,
                controls,
                feedback,
                local.flow,
            )
        }
        .map_err(|source| ParallelGasError::Branch {
            index,
            name: branch.name.to_owned(),
            source,
        })?;
        evaluations += solved.flow_evaluations;
        substeps += solved.substeps;
        if inlet.is_none() {
            inlet = Some(
                solved
                    .gas
                    .thermal
                    .field()
                    .state_at(0.0)
                    .map_err(GasEnergyError::from)?,
            );
        }
        let inlet_rho = inlet
            .as_ref()
            .expect("first accepted branch")
            .state()
            .density;
        let outlet = solved
            .gas
            .cells
            .last()
            .ok_or(invalid("native branch lost its outlet"))?;
        let m = solved.mass_flow_kg_per_s;
        let area = branch.transport[0].resistance.flow_area().value();
        let inlet_speed = (m / inlet_rho) / area;
        let outlet_speed = outlet.velocity_m_s;
        let endpoint_kinetic = (0.5 * solved.kinetic_energy_flux_factor * m)
            * (inlet_speed * inlet_speed + outlet_speed * outlet_speed);
        mass += m;
        heat += solved
            .gas
            .cells
            .iter()
            .map(|cell| cell.wall_heat_w)
            .sum::<f64>();
        gain += m * law
            .relative_enthalpy_and_heat_capacity(
                first.inlet.temperature_k,
                outlet.gas.state().temperature,
            )?
            .0;
        hydraulic += solved.gas.hydraulic_dissipation_w;
        manifold += endpoint_kinetic;
        kinetic += solved.neglected_kinetic_power_bound_w + endpoint_kinetic;
        potential += solved
            .gravity
            .as_ref()
            .map_or(0.0, |state| state.neglected_potential_power_bound_w);
        if !positive(mass)
            || ![heat, gain, hydraulic, kinetic, potential, manifold]
                .iter()
                .all(|x| x.is_finite())
        {
            return Err(invalid(
                "aggregate parallel source mass, enthalpy or mechanical power overflows",
            ));
        }
        if hydraulic > budget.maximum_neglected_hydraulic_power_w
            || kinetic > budget.maximum_neglected_kinetic_power_w
            || potential > budget.maximum_neglected_potential_power_w
        {
            return Err(ParallelGasError::OmittedPower {
                hydraulic_w: hydraulic,
                kinetic_w: kinetic,
                potential_w: potential,
            });
        }
        results.push(ParallelGasBranchState {
            name: branch.name.to_owned(),
            solution: solved,
        });
    }
    let mixed = law.model().reconstruct_pressure_enthalpy_increment(
        law.context(),
        gain / mass,
        controls.caloric,
        || cx.checkpoint().is_err(),
    )?;
    // Retain the original gas/water phase admission at the new mixed temperature.
    let point = law
        .context()
        .clone()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            mixed.state().state().temperature,
        )
        .map_err(|_| invalid("mixed temperature coordinate refused"))?;
    let mixed_outlet = match first.phase {
        ThermochemicalPhase::FrozenGas => mixed.state().clone(),
        ThermochemicalPhase::LiquidWaterLimit | ThermochemicalPhase::LiquidWaterIapwsLimit => {
            if first.source.water_species().is_none() {
                return Err(invalid(
                    "liquid-water phase requires an original source water species",
                ));
            }
            if law.model().water_species().is_none() {
                mixed.state().clone()
            } else if first.phase == ThermochemicalPhase::LiquidWaterIapwsLimit {
                law.model()
                    .evaluate_with_iapws_liquid_water_limit(&point)
                    .map_err(GasEnergyError::from)?
            } else {
                law.model()
                    .evaluate_with_liquid_water_limit(&point)
                    .map_err(GasEnergyError::from)?
            }
        }
    };
    let inlet = inlet.expect("nonempty accepted branches");
    let inlet_volume = mass / inlet.state().density;
    let outlet_volume = mass / mixed_outlet.state().density;
    let mixing_residual = mass * mixed.residual_j_kg();
    if !positive(inlet_volume) || !positive(outlet_volume) || !mixing_residual.is_finite() {
        return Err(invalid(
            "mixed volume flow or caloric residual is not representable",
        ));
    }
    checkpoint(cx)?;
    Ok(ParallelGasSolution {
        branches: results,
        driving_pressure_pa: drive.value(),
        inlet,
        mixed_outlet,
        total_mass_flow_kg_per_s: mass,
        inlet_volume_flow_m3_per_s: inlet_volume,
        outlet_volume_flow_m3_per_s: outlet_volume,
        total_wall_heat_w: heat,
        total_enthalpy_gain_w: gain,
        mixing_enthalpy_residual_w: mixing_residual,
        neglected_hydraulic_power_w: hydraulic,
        neglected_kinetic_power_bound_w: kinetic,
        neglected_potential_power_bound_w: potential,
        manifold_kinetic_power_bound_w: manifold,
        flow_evaluations: evaluations,
        substeps,
    })
}
