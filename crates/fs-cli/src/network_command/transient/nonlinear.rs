//! Product adapter for nonlinear solid endpoints; no duplicate material model.
use super::*;
use fs_conduction::{ConductionError, ConductionProblem};
use fs_conduction::solve::LineSearch;
use fs_conduction::transient::backward_euler::{NonlinearStepConfig, NonlinearStepSolution};

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    pub(super) policy: NonlinearStepConfig,
}

impl Config {
    pub(super) fn parse(value: &J) -> Result<Self> {
        object(value, &["max_iterations", "residual_rtol", "residual_atol_j",
            "armijo_c", "shrink", "max_backtracks"], "transient.nonlinear")?;
        let policy = NonlinearStepConfig {
            max_iterations: count(get(value,"max_iterations")?,"nonlinear.max_iterations",10_000)?,
            residual_rtol: positive(get(value,"residual_rtol")?,"nonlinear.residual_rtol")?,
            residual_atol_j: number(get(value,"residual_atol_j")?,"nonlinear.residual_atol_j")?,
            line_search: LineSearch {
                armijo_c: number(get(value,"armijo_c")?,"nonlinear.armijo_c")?,
                shrink: number(get(value,"shrink")?,"nonlinear.shrink")?,
                max_backtracks: integer(get(value,"max_backtracks")?,"nonlinear.max_backtracks",64)?,
            },
        };
        policy.validate().map_err(producer)?;
        Ok(Self { policy })
    }
}

/// Check actual assigned materials, not the inactive scalar fallback. This is
/// also called when a schedule is supplied directly by an internal consumer.
pub(super) fn admit(request: &Request, cx: &Cx<'_>, config: Option<Config>) -> Result<()> {
    if let Some(config) = config { config.policy.validate().map_err(producer)?; }
    if let Some(materials) = &request.solid_data.element_materials {
        for element in 0..request.mesh.element_count() {
            if element % 512 == 0 { poll(cx)?; }
            if materials.model_for(element).map_err(producer)?.is_temperature_dependent()
                && config.is_none()
            {
                return Err(bad("transient k(T) requires an explicit transient.nonlinear policy; conductivity is never frozen at the old temperature"));
            }
        }
    }
    poll(cx)
}

/// Work from every successful solid endpoint evaluation, including discarded
/// coupling and adaptive trials. Their heat never enters physical history.
#[derive(Debug, Default)]
pub(super) struct Stats {
    solves: usize,
    updates: usize,
    krylov: usize,
    backtracks: usize,
    worst_residual_ratio: f64,
}

impl Stats {
    fn observe(&mut self, result: &NonlinearStepSolution) -> Result<()> {
        fn add(a: usize, b: usize) -> Result<usize> {
            a.checked_add(b).ok_or_else(|| budget("nonlinear transient work counter overflow"))
        }
        self.solves = add(self.solves,1)?;
        self.updates = add(self.updates,result.nonlinear_iterations)?;
        self.krylov = add(self.krylov,result.step.krylov_iterations)?;
        self.backtracks = add(self.backtracks,result.backtracks)?;
        let ratio = if result.threshold_j == 0.0 { 0.0 }
            else { finite(result.residual_j/result.threshold_j)? };
        self.worst_residual_ratio = self.worst_residual_ratio.max(ratio);
        Ok(())
    }

    pub(super) fn render(&self, config: Option<Config>) -> Result<String> {
        let Some(config) = config else { return Ok("null".into()); };
        let p = config.policy;
        Ok(format!(
            "{{\"method\":\"endpoint-newton-fgmres\",\"max_iterations_per_solid_solve\":{},\"residual_rtol\":{},\"residual_atol_j\":{},\"armijo_c\":{},\"shrink\":{},\"max_backtracks_per_update\":{},\"solid_solves\":{},\"newton_updates\":{},\"krylov_iterations\":{},\"backtracks\":{},\"worst_accepted_residual_ratio\":{},\"scope\":\"endpoint k(T), constant heat capacity and contact resistance; includes successful discarded coupling/adaptive trial work; no history advance before physical acceptance; linear_iterations caps all inner iterations of each solid endpoint solve\"}}",
            p.max_iterations,num(p.residual_rtol)?,num(p.residual_atol_j)?,
            num(p.line_search.armijo_c)?,num(p.line_search.shrink)?,p.line_search.max_backtracks,
            self.solves,self.updates,self.krylov,self.backtracks,num(self.worst_residual_ratio)?,
        ))
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    engine: &BackwardEuler<'_>, cx: &Cx<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&fs_conduction::ThermalInterfaces>, old: &[f64], dt: f64,
    linear: StepConfig, config: Option<Config>, stats: &mut Stats,
) -> Result<StepSolution> {
    match config {
        None => engine.advance(cx,problem,interfaces,old,dt,linear).map_err(producer),
        Some(config) => {
            let result = engine.advance_nonlinear(cx,problem,interfaces,old,dt,linear,config.policy)
                .map_err(|error| match error {
                    ConductionError::NotConverged { .. } => Failure {
                        code:"cooling-network-transient-budget", message:error.to_string(),
                    },
                    other => producer(other),
                })?;
            stats.observe(&result)?;
            Ok(result.step)
        }
    }
}
