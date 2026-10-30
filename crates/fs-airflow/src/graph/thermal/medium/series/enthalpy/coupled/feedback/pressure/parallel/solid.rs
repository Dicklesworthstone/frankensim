//! One solid solve coupled to all material-dependent pressure branches.
//!
//! The unknown is the complete wall-temperature vector. Each map sample first
//! solves every branch operating point, sends their actual Robin rows to ONE
//! common solid solve, then recomputes all branch flows at its returned walls.
//! The solid callback must always evaluate from the same accepted history.

use fs_convection::CorrelationRobinBoundary;
use fs_couple::AitkenRelaxation;
use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig};
use fs_exec::Cx;
use fs_qty::{Pressure, Temperature, VolumetricFlowRate};

use crate::AirflowError;
use crate::conjugate::{ConjugateConfig, Relaxation};
use crate::graph::thermal::{admit, checked, iqn_error};
use crate::graph::thermal::medium::{FlowHeatSection, FrozenFlowMedium, MediumFlowError};
use crate::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use super::LiquidResistanceModel;
use super::{ParallelLiquidNetwork, ParallelLiquidSolution, ParallelPressureBudget,
    ParallelPressureError, LiquidFeedbackBudget, LiquidFeedbackError, checkpoint, positive};

impl From<AirflowError> for ParallelPressureError {
    fn from(error: AirflowError) -> Self { Self::from(LiquidFeedbackError::from(error)) }
}
impl From<MediumFlowError> for ParallelPressureError {
    fn from(error: MediumFlowError) -> Self { Self::from(LiquidFeedbackError::from(error)) }
}

/// Explicit numerical controls; none modifies source or correlation validity.
#[derive(Clone, Copy, Debug)]
pub struct ParallelLiquidCouplingControls {
    /// Wall/reference and per-region power tolerances; maximum solid solves.
    pub coupling: ConjugateConfig,
    /// Local material/coefficient/outlet convergence controls.
    pub local: LiquidFeedbackBudget,
    /// Existing vector IQN-ILS history/rank policy, shared across ALL walls.
    pub acceleration: IqnIlsConfig,
    /// Require each recomputed branch mass flow to agree to this relative error.
    pub mass_flow_relative_tolerance: f64,
    /// One aggregate allowance across every wall iteration and pressure trial.
    pub work: ParallelPressureBudget,
}

/// Unrelaxed changes and local power mismatch in the last joint sample.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParallelLiquidResiduals {
    /// Maximum absolute wall-temperature change, K.
    pub wall_temperature_k: f64,
    /// Maximum absolute effective coolant-reference change, K.
    pub reference_temperature_k: f64,
    /// Maximum relative change of an actual convection coefficient.
    pub coefficient_relative: f64,
    /// Maximum relative change of an actual branch mass flow.
    pub mass_flow_relative: f64,
    /// Maximum absolute regional solid/coolant power mismatch, W.
    pub heat_imbalance_w: f64,
}

/// Original solid result and its exact sent/recomputed material-flow states.
#[derive(Debug)]
pub struct CoupledParallelLiquidSolution<'a, R> {
    /// Actual common-solid field/report, evaluated only once per iteration.
    pub solid: LiquidSolidResponse<R>,
    /// Exact original branch flows and Robin rows supplied to that solid solve.
    pub boundary_state: ParallelLiquidSolution<'a>,
    /// Recomputed branch flows at the returned solid's actual wall temperatures.
    pub fluid: ParallelLiquidSolution<'a>,
    /// Number of common-solid callbacks, not a sum of independent branch solves.
    pub iterations: usize,
    /// Full branch thermal evaluations across every nonlinear/pressure trial.
    pub total_flow_evaluations: usize,
    /// Final unrelaxed physical residuals, separate from relaxation size.
    pub residuals: ParallelLiquidResiduals,
}

impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    /// Close pressure-driven flow redistribution against one common solid.
    ///
    /// The callback receives every actual Robin row in branch-major/stream-wise
    /// order. It must report the sent reference for every matching solid region.
    /// All candidates must start from the same accepted solid history. Original
    /// fields and material/pressure states are published only after wall, fluid
    /// reference, h, branch flow and EACH regional power gate pass. Opposing
    /// branch/region errors cannot cancel into an accepted global total.
    ///
    /// Uses the existing IQN-ILS accelerator and fixed/Aitken fallback. Numerical
    /// extrapolation outside the caloric span discards secants; unsupported
    /// actual source/correlation/solid states remain terminal. The existing
    /// low-compressibility, ideal-manifold and omitted-heating limits remain.
    pub fn solve_coupled_iqn<R, F>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], drive: Pressure,
        controls: &ParallelLiquidCouplingControls, mut solid: F,
    ) -> Result<CoupledParallelLiquidSolution<'a, R>, ParallelPressureError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        checkpoint(cx)?;
        controls.local.admit()?;
        controls.work.admit()?;
        if !positive(controls.mass_flow_relative_tolerance) || controls.mass_flow_relative_tolerance >= 1.0 {
            return Err(ParallelPressureError::Invalid("branch flow tolerance must be finite and strictly between zero and one"));
        }
        let mut templates = Vec::with_capacity(self.branches.len());
        for branch in &self.branches {
            checkpoint(cx)?;
            branch.pressure_budget.admit(drive)?;
            let path = branch.path;
            let descriptions: Vec<_> = path.sections.iter().map(|s| FlowHeatSection {
                region: &s.name, geometry: s.geometry, correlation: s.correlation,
                direction: fs_convection::ThermalDirection::HeatingFluid,
            }).collect();
            let material = path.law.reference_state();
            templates.push(FrozenFlowMedium::Liquid(material).channel_series(
                VolumetricFlowRate::new(path.mass_flow/material.density_kg_m3()),
                Temperature::new(path.inlet.temperature_k()), &descriptions,
            )?);
        }
        // Reuse existing configuration and globally unique region admission.
        let paths: Vec<_> = templates.iter().map(|p| p.path().clone()).collect();
        admit(&paths, &controls.coupling, initial_walls_k)?;
        let mut accelerator = IqnIls::new(self.region_count, controls.acceleration).map_err(iqn_error)?;
        let mut relaxer = match controls.coupling.relaxation {
            Relaxation::Fixed { .. } => None,
            Relaxation::Aitken { omega_init, omega_max } => Some(AitkenRelaxation::new(omega_init,omega_max)),
        };
        let areas: Vec<_> = paths.iter().flat_map(|p| p.segments().iter().map(|s| s.area_m2())).collect();
        let total_area = checked("parallel total exchange area", areas.iter().sum())?;
        let mut walls = initial_walls_k.to_vec();
        let mut performed = 0;
        let mut residuals = ParallelLiquidResiduals::default();
        for iteration in 0..controls.coupling.max_iterations {
            let boundary_state = self.march_metered(cx,&walls,drive,controls.local,controls.work,&mut performed)?;
            let before: Vec<_> = boundary_state.branches.iter()
                .flat_map(|b| b.solution.thermal.regions.iter()).collect();
            let rows: Vec<_> = before.iter().map(|r| r.boundary.clone()).collect();
            let references: Vec<_> = rows.iter().map(|r| r.reference_temperature().value()).collect();
            checkpoint(cx)?;
            let response = solid(cx,&rows)?;
            checkpoint(cx)?;
            if response.regions.len() != self.region_count {
                return Err(AirflowError::SolidResponseArity { expected: self.region_count,found: response.regions.len() }.into());
            }
            let mut offset = 0;
            for template in &templates {
                let end = offset+template.path().segments().len();
                template.check_enthalpy_solid_response(&references[offset..end],&response.regions[offset..end])?;
                offset = end;
            }
            let next: Vec<_> = response.regions.iter().map(|r| r.mean_wall_temperature_k).collect();
            let fluid = self.march_metered(cx,&next,drive,controls.local,controls.work,&mut performed)?;
            let after = fluid.branches.iter().flat_map(|b| b.solution.thermal.regions.iter());
            residuals = ParallelLiquidResiduals::default();
            let mut weighted = 0.0;
            let mut failed_balance = None;
            for (index,((old,new),body)) in before.iter().zip(after).zip(&response.regions).enumerate() {
                let delta = checked("parallel wall residual",next[index]-walls[index])?;
                weighted = checked("parallel weighted wall residual",weighted+(areas[index]/total_area)*delta)?;
                residuals.wall_temperature_k = residuals.wall_temperature_k.max(delta.abs());
                let ref_delta = checked("parallel coolant reference residual",new.boundary.reference_temperature().value()-references[index])?;
                residuals.reference_temperature_k = residuals.reference_temperature_k.max(ref_delta.abs());
                let (old_h,new_h) = (old.boundary.coefficient().value.value(),new.boundary.coefficient().value.value());
                residuals.coefficient_relative = residuals.coefficient_relative.max((new_h-old_h).abs()/new_h.max(old_h));
                let imbalance = checked("parallel regional heat imbalance",body.heat_rate_w-new.heat_rate_w)?.abs();
                residuals.heat_imbalance_w = residuals.heat_imbalance_w.max(imbalance);
                let threshold = controls.coupling.balance_tolerance_w.max(controls.coupling.balance_relative_tolerance
                    * body.heat_rate_w.abs().max(new.heat_rate_w.abs()));
                if imbalance > threshold && failed_balance.is_none_or(|(worst,_)| imbalance > worst) {
                    failed_balance = Some((imbalance,threshold));
                }
            }
            for (before,after) in boundary_state.branches.iter().zip(&fluid.branches) {
                let (old,new) = (before.solution.mass_flow_kg_per_s,after.solution.mass_flow_kg_per_s);
                residuals.mass_flow_relative = residuals.mass_flow_relative.max((new-old).abs()/new.max(old));
            }
            let states_closed = residuals.wall_temperature_k <= controls.coupling.temperature_tolerance_k
                && residuals.reference_temperature_k <= controls.coupling.temperature_tolerance_k
                && residuals.coefficient_relative <= controls.local.coefficient_relative_tolerance
                && residuals.mass_flow_relative <= controls.mass_flow_relative_tolerance;
            if states_closed && failed_balance.is_none() {
                checkpoint(cx)?;
                return Ok(CoupledParallelLiquidSolution { solid: response,boundary_state,fluid,
                    iterations: iteration+1,total_flow_evaluations: performed,residuals });
            }
            if states_closed && (residuals.wall_temperature_k == 0.0 || iteration+1 == controls.coupling.max_iterations)
                && let Some((imbalance,threshold)) = failed_balance
            {
                return Err(AirflowError::ConjugateBalanceUnclosed { iterations: iteration+1,
                    max_region_imbalance_bits: imbalance.to_bits(),tolerance_bits: threshold.to_bits() }.into());
            }
            if iteration+1 == controls.coupling.max_iterations { break; }
            let omega = match (&controls.coupling.relaxation,relaxer.as_mut()) {
                (Relaxation::Fixed { omega },_) => *omega,
                (Relaxation::Aitken { .. },Some(r)) => r.next_omega(weighted),
                (Relaxation::Aitken { omega_init,.. },None) => *omega_init,
            };
            let mut proposed = accelerator.step(&walls,&next,omega).map_err(iqn_error)?.values;
            let knots = self.branches[0].path.source.curve().knots();
            let inside = |t: &f64| t.is_finite() && *t >= knots[0].temperature_k && *t <= knots[knots.len()-1].temperature_k;
            if !proposed.iter().all(inside) {
                accelerator.reset();
                proposed = walls.iter().zip(&next).map(|(&old,&new)| old+omega*(new-old)).collect();
                if !proposed.iter().all(inside) {
                    return Err(ParallelPressureError::Invalid("configured relaxation leaves the caloric wall span"));
                }
            }
            checkpoint(cx)?;
            walls = proposed;
        }
        Err(ParallelPressureError::CouplingNotConverged { iterations: controls.coupling.max_iterations,residuals })
    }
}
