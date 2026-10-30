//! Material-driven parallel coolant branches between ideal pressure manifolds.
//!
//! Every branch solves its own existing pressure/temperature operating point.
//! They share one liquid source, inlet state and differential pressure, not a
//! frozen split fraction. Their outlets mix by conserved specific enthalpy on
//! the SAME source chart. No second material, caloric or resistance law lives
//! here. Junction losses must be included in the authored branch resistances;
//! no junction geometry, recirculation or transient storage is inferred.

/// Common-solid feedback over all pressure-driven material branches.
pub mod solid;
/// Equipment-source operating points over the same thermal network.
pub mod pump;

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::phase::EquilibriumPhaseState;
use fs_qty::{Pressure, Temperature};

use super::porous::{LiquidResistanceModel, path_outlet_pressure};
use super::MaterialFlowResistance;
use super::{LiquidEnthalpyError, LiquidFeedbackBudget, LiquidFeedbackError,
    LiquidHeatMarch, LiquidHeatPath, LiquidResistanceSection, PressureDrivenBudget,
    PressureDrivenError, PressureDrivenSolution, checkpoint, positive};

/// Hard bound on simultaneously retained hydraulic branches.
pub const MAX_PARALLEL_LIQUID_BRANCHES: usize = 64;
/// Bound on the complete wall interface, not a separate allowance per branch.
pub const MAX_PARALLEL_LIQUID_REGIONS: usize = 1024;

/// A complete existing thermal path and its exact pressure-loss assignments.
#[derive(Clone, Copy, Debug)]
pub struct ParallelLiquidBranch<'p, 'a, 'card, H = MaterialFlowResistance> {
    /// Unique branch label; thermal region names must also be globally unique.
    pub name: &'p str,
    /// Existing source-backed thermal description, including its common inlet.
    pub path: &'p LiquidHeatPath<'a, 'card>,
    /// Pressure loss assignments in this path's exact thermal-region order.
    pub resistances: &'p [LiquidResistanceSection<'p, H>],
    /// Explicit supported mass-flow bracket and local pressure tolerances.
    pub pressure_budget: PressureDrivenBudget,
}

/// A work and omitted-power allowance shared by the ENTIRE network calculation.
#[derive(Clone, Copy, Debug)]
pub struct ParallelPressureBudget {
    /// Total full branch thermal evaluations across every pressure trial.
    pub max_flow_evaluations: usize,
    /// Total omitted hydraulic dissipation, W, over all branches together.
    /// Branch-local omission allowances remain active as well.
    pub maximum_neglected_hydraulic_power_w: f64,
}
impl ParallelPressureBudget {
    fn admit(self) -> Result<(), ParallelPressureError> {
        if !(2..=1_048_576).contains(&self.max_flow_evaluations)
            || !self.maximum_neglected_hydraulic_power_w.is_finite()
            || self.maximum_neglected_hydraulic_power_w < 0.0
        {
            return Err(ParallelPressureError::Invalid("invalid whole-network work or omitted-power allowance"));
        }
        Ok(())
    }
}

/// A network refuses atomically; successful earlier branches are not an answer.
#[derive(Debug)]
pub enum ParallelPressureError {
    /// Incompatible source/inlet, region ownership or numerical configuration.
    Invalid(&'static str),
    /// The original named branch refused its material, correlation or flow solve.
    Branch { branch: String, source: PressureDrivenError },
    /// Shared source, caloric, cancellation or solid-solver refusal.
    Physics(PressureDrivenError),
    /// The shared thermal-evaluation allowance was exhausted.
    FlowBudget { performed: usize, maximum: usize },
    /// The common-solid calculation exhausted its joint convergence budget.
    CouplingNotConverged { iterations: usize, residuals: solid::ParallelLiquidResiduals },
    /// Branch-local allowances cannot hide excessive aggregate omitted power.
    NeglectedDissipation { power_w: f64, allowance_w: f64 },
}
impl core::fmt::Display for ParallelPressureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "parallel material flow: {self:?}")
    }
}
impl std::error::Error for ParallelPressureError {}
impl From<PressureDrivenError> for ParallelPressureError {
    fn from(error: PressureDrivenError) -> Self { Self::Physics(error) }
}
impl From<LiquidEnthalpyError> for ParallelPressureError {
    fn from(error: LiquidEnthalpyError) -> Self { Self::Physics(error.into()) }
}
impl From<LiquidFeedbackError> for ParallelPressureError {
    fn from(error: LiquidFeedbackError) -> Self { Self::Physics(error.into()) }
}

/// Named branch result at the common manifold pressure.
#[derive(Debug)]
pub struct ParallelLiquidBranchState<'a> {
    /// Original branch label, never inferred from result ordering.
    pub name: String,
    /// Original pressure, material and caloric result for this branch.
    pub solution: PressureDrivenSolution<LiquidHeatMarch<'a>>,
}

/// Complete material-consistent split, transport and adiabatic outlet mixing.
#[derive(Debug)]
pub struct ParallelLiquidSolution<'a> {
    /// Branch-major results; each branch retains stream-wise region order.
    pub branches: Vec<ParallelLiquidBranchState<'a>>,
    /// Actual common inlet material, including its source receipts.
    pub inlet_material: ResolvedLiquidState,
    /// Pressure imposed across every branch, Pa.
    pub driving_pressure_pa: f64,
    /// Sum of branch mass flows, kg/s.
    pub total_mass_flow_kg_per_s: f64,
    /// Total supply volume flow using density at the COMMON INLET, m3/s.
    pub inlet_volume_flow_m3_per_s: f64,
    /// Source-chart state from enthalpy mixing, not temperature averaging.
    pub mixed_outlet: EquilibriumPhaseState,
    /// Sum of signed coolant enthalpy transfers, W.
    pub total_heat_rate_w: f64,
    /// Sum of all local hydraulic power omissions, W.
    pub total_dissipated_power_w: f64,
    /// Full branch thermal evaluations consumed by this network march.
    pub flow_evaluations: usize,
}

/// Parallel topology with one common source and inlet, and no shared wall rows.
#[derive(Debug)]
pub struct ParallelLiquidNetwork<'p, 'a, 'card, H = MaterialFlowResistance> {
    branches: Vec<ParallelLiquidBranch<'p, 'a, 'card, H>>,
    region_count: usize,
}
impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    /// Admit material/topology ownership before running a hydraulic or solid solve.
    /// All paths must retain identical liquid-law reference states and the same
    /// caloric chart identity and inlet. Equivalent independently resolved source
    /// objects are accepted; pointer identity is not required. Mixing unrelated
    /// material cards, phases or enthalpy reference conventions is not supported.
    pub fn try_new(
        cx: &Cx<'_>, branches: &[ParallelLiquidBranch<'p, 'a, 'card, H>],
    ) -> Result<Self, ParallelPressureError> {
        checkpoint(cx)?;
        if branches.is_empty() || branches.len() > MAX_PARALLEL_LIQUID_BRANCHES {
            return Err(ParallelPressureError::Invalid("parallel branch count is outside its bounded nonempty domain"));
        }
        let first = branches[0].path;
        let mut names = BTreeSet::new();
        let mut regions = BTreeSet::new();
        let mut outlet_pressure = None;
        for branch in branches {
            checkpoint(cx)?;
            if branch.name.trim().is_empty() || !names.insert(branch.name) {
                return Err(ParallelPressureError::Invalid("parallel branch names must be nonempty and unique"));
            }
            if let Some(value) = path_outlet_pressure(branch.resistances)? {
                if outlet_pressure.is_some_and(|previous| previous != value) {
                    return Err(ParallelPressureError::Invalid("parallel porous branches require one absolute outlet manifold pressure"));
                }
                outlet_pressure = Some(value);
            }
            let path = branch.path;
            if path.inlet != first.inlet
                || path.source.curve().identity() != first.source.curve().identity()
                || path.law.reference_state() != first.law.reference_state()
                || path.law.reference_point() != first.law.reference_point()
                || path.law.temperature_axis() != first.law.temperature_axis()
                || path.law.viscosity_basis() != first.law.viscosity_basis()
            {
                return Err(ParallelPressureError::Invalid("parallel manifolds require the same liquid source, caloric chart and inlet"));
            }
            path.admit_resistances(branch.resistances).map_err(|source| ParallelPressureError::Branch {
                branch: branch.name.to_owned(), source,
            })?;
            for section in &path.sections {
                if !regions.insert(section.name.as_str()) {
                    return Err(ParallelPressureError::Invalid("a solid wall region cannot belong to two parallel streams"));
                }
                if regions.len() > MAX_PARALLEL_LIQUID_REGIONS {
                    return Err(ParallelPressureError::Invalid("whole-network region count exceeds its bound"));
                }
            }
        }
        let region_count = regions.len();
        checkpoint(cx)?;
        Ok(Self { branches: branches.to_vec(), region_count })
    }

    /// Number of wall rows in the complete branch-major interface.
    #[must_use]
    pub const fn region_count(&self) -> usize { self.region_count }

    /// Original descriptions with their independent, supported flow brackets.
    #[must_use]
    pub fn branches(&self) -> &[ParallelLiquidBranch<'p, 'a, 'card, H>] { &self.branches }

    /// Solve branch flow redistribution at common pressure and supplied walls.
    /// Every branch reuses the existing pressure/material/enthalpy solver; there
    /// is no frozen hydraulic scaling. Outlet mixing assumes an ideal adiabatic
    /// zero-storage manifold at the common source pressure. Viscous heating,
    /// kinetic-energy recovery, junction geometry and recirculation are excluded.
    pub fn solve_pressure_driven(
        &self, cx: &Cx<'_>, walls: &[f64], drive: Pressure,
        local: LiquidFeedbackBudget, budget: ParallelPressureBudget,
    ) -> Result<ParallelLiquidSolution<'a>, ParallelPressureError> {
        let mut performed = 0;
        self.march_metered(cx, walls, drive, local, budget, &mut performed)
    }

    // A coupled owner shares this exact meter across all wall/source trials.
    fn march_metered(
        &self, cx: &Cx<'_>, walls: &[f64], drive: Pressure,
        local: LiquidFeedbackBudget, budget: ParallelPressureBudget, performed: &mut usize,
    ) -> Result<ParallelLiquidSolution<'a>, ParallelPressureError> {
        checkpoint(cx)?;
        budget.admit()?;
        local.admit()?;
        if walls.len() != self.region_count {
            return Err(ParallelPressureError::Invalid("one wall temperature per branch-major region is required"));
        }
        // Admit all descriptions/walls before executing even the first branch.
        let first = self.branches[0].path;
        for branch in &self.branches { branch.pressure_budget.admit(drive)?; }
        for &wall in walls { super::super::state_at_temperature(first.source.curve(), wall)?; }
        let inlet_material = first.law.at_temperature(Temperature::new(first.inlet.temperature_k()))
            .map_err(LiquidFeedbackError::from)?;
        let start = *performed;
        let mut results = Vec::with_capacity(self.branches.len());
        let mut offset = 0;
        let mut mass = 0.0;
        let mut heat = 0.0;
        let mut power = 0.0;
        for branch in &self.branches {
            checkpoint(cx)?;
            let remaining = budget.max_flow_evaluations.saturating_sub(*performed);
            if remaining < 2 {
                return Err(ParallelPressureError::FlowBudget { performed: *performed, maximum: budget.max_flow_evaluations });
            }
            let count = branch.path.sections.len();
            let branch_budget = PressureDrivenBudget {
                max_flow_evaluations: branch.pressure_budget.max_flow_evaluations.min(remaining),
                ..branch.pressure_budget
            };
            let solution = branch.path.solve_pressure_driven(cx, &walls[offset..offset+count],
                branch.resistances, drive, local, branch_budget)
                .map_err(|source| ParallelPressureError::Branch { branch: branch.name.to_owned(), source })?;
            *performed += solution.flow_evaluations;
            offset += count;
            mass += solution.mass_flow_kg_per_s;
            heat += solution.thermal.total_heat_rate_w;
            power += solution.hydraulics.total_dissipated_power_w;
            if !positive(mass) || !heat.is_finite() || !power.is_finite() {
                return Err(ParallelPressureError::Invalid("aggregate parallel mass, heat or hydraulic power is unrepresentable"));
            }
            if power > budget.maximum_neglected_hydraulic_power_w {
                return Err(ParallelPressureError::NeglectedDissipation {
                    power_w: power, allowance_w: budget.maximum_neglected_hydraulic_power_w,
                });
            }
            results.push(ParallelLiquidBranchState { name: branch.name.to_owned(), solution });
        }
        // Subtract the common enthalpy reference BEFORE weighting. Thus a large
        // arbitrary enthalpy origin cannot overflow individual mdot*h products.
        let mixed_h = first.inlet.specific_enthalpy_j_kg() + heat/mass;
        let mixed_outlet = first.source.curve().state_at_specific_enthalpy(mixed_h)
            .map_err(LiquidEnthalpyError::from)?;
        let inlet_volume = mass/inlet_material.density_kg_m3();
        if !positive(inlet_volume) {
            return Err(ParallelPressureError::Invalid("common inlet volume flow is unrepresentable"));
        }
        checkpoint(cx)?;
        Ok(ParallelLiquidSolution { branches: results, inlet_material, driving_pressure_pa: drive.value(),
            total_mass_flow_kg_per_s: mass, inlet_volume_flow_m3_per_s: inlet_volume,
            mixed_outlet, total_heat_rate_w: heat, total_dissipated_power_w: power,
            flow_evaluations: *performed-start })
    }
}
