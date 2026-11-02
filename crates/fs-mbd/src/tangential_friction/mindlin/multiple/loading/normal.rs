//! Original Hertz residual and analytic Jacobian for the shared normal load line.
//! Globalization, Krylov iteration and numerical root selection stay in fs-solver.
use super::*;
use fs_contact::normal_patch::source::path::SourceElasticPathAverage;
use fs_solver::{NewtonKrylovState, NonlinearProblem, SolverRunProgress};
use std::cell::Cell;

pub(super) struct NormalSolution {
    pub indentations: Vec<f64>,
    pub paths: Vec<SourceElasticPathAverage>,
    pub residual_m: f64,
    pub iterations: usize,
    pub evaluations: usize,
}
struct Cache {
    x: Vec<f64>,
    paths: Vec<SourceElasticPathAverage>,
}
struct Problem<'a, 'b> {
    island: &'a PreparedRigidHertzMindlinIsland<'b>,
    cache: RefCell<Cache>,
    evaluations: Cell<usize>,
    exhausted: Cell<bool>,
    last_refusal: RefCell<Option<NormalPatchError>>,
}
impl Problem<'_, '_> {
    fn evaluate(&self, x: &[f64]) -> bool {
        if self.exhausted.get() {
            return false;
        }
        if self.cache.borrow().x == x {
            return true;
        }
        // Invalid trial depths have no constitutive authority. The ORIGINAL
        // Newton line search may shorten that numerical step; accepted points
        // still require the exact original source box/rate admission below.
        if x.iter().any(|d| !d.is_finite() || *d <= 0.0) {
            return false;
        }
        let mut paths = Vec::with_capacity(x.len());
        for (patch, &indentation) in self.island.patches.iter().zip(x) {
            if self.evaluations.get() == self.island.controls.maximum_normal_evaluations {
                self.exhausted.set(true);
                return false;
            }
            self.evaluations.set(self.evaluations.get() + 1);
            let source = patch.law.fixed_source().normal_source();
            match source.elastic_path_average(
                patch.history.indentation_m(),
                indentation,
                self.island.dt_s,
            ) {
                Ok(path) => paths.push(path),
                Err(error) => {
                    *self.last_refusal.borrow_mut() = Some(error);
                    return false;
                }
            }
        }
        *self.cache.borrow_mut() = Cache {
            x: x.to_vec(),
            paths,
        };
        true
    }
}
impl NonlinearProblem for Problem<'_, '_> {
    fn dimension(&self) -> usize {
        self.island.patches.len()
    }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.evaluate(x) {
            return;
        }
        let cache = self.cache.borrow();
        let n = self.dimension();
        for (i, patch) in self.island.patches.iter().enumerate() {
            let normal_recoil: f64 = (0..n)
                .map(|j| self.island.normal_compliance[i * n + j] * cache.paths[j].mean_force_n)
                .sum();
            output[i] = x[i] - patch.history.indentation_m() - patch.released[0] + normal_recoil;
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.evaluate(x) {
            return;
        }
        let cache = self.cache.borrow();
        let n = self.dimension();
        for i in 0..n {
            output[i] = direction[i]
                + (0..n)
                    .map(|j| {
                        self.island.normal_compliance[i * n + j]
                            * cache.paths[j].endpoint_derivative_n_per_m
                            * direction[j]
                    })
                    .sum::<f64>();
        }
    }
    fn preconditioner_apply(
        &self,
        x: &[f64],
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.fill(f64::NAN);
        if !self.evaluate(x) {
            return;
        }
        let cache = self.cache.borrow();
        let n = self.dimension();
        for i in 0..n {
            output[i] = residual[i]
                / (1.0
                    + self.island.normal_compliance[i * n + i]
                        * cache.paths[i].endpoint_derivative_n_per_m);
        }
    }
}

pub(super) fn solve(
    island: &PreparedRigidHertzMindlinIsland<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NormalSolution, RigidHertzMindlinIslandError> {
    poll(cancelled)?;
    let problem = Problem {
        island,
        cache: RefCell::new(Cache {
            x: vec![],
            paths: vec![],
        }),
        evaluations: Cell::new(0),
        exhausted: Cell::new(false),
        last_refusal: RefCell::new(None),
    };
    let initial: Vec<_> = island
        .patches
        .iter()
        .map(|p| p.history.indentation_m())
        .collect();
    let construction = NewtonKrylovState::new(&problem, initial, island.controls.normal_solver);
    if problem.exhausted.get() {
        return Err(RigidHertzMindlinIslandError::NormalEvaluationBudget {
            evaluations: problem.evaluations.get(),
        });
    }
    let mut solver = construction?;
    poll(cancelled)?;
    let result = solver.run_controlled(
        &problem,
        island.controls.maximum_normal_iterations,
        &mut || cancelled() || problem.exhausted.get(),
    );
    if problem.exhausted.get() {
        return Err(RigidHertzMindlinIslandError::NormalEvaluationBudget {
            evaluations: problem.evaluations.get(),
        });
    }
    if result.progress == SolverRunProgress::Paused {
        return Err(RigidHertzMindlinIslandError::Cancelled);
    }
    poll(cancelled)?;
    if result.progress != SolverRunProgress::Complete || !result.report.converged {
        return Err(RigidHertzMindlinIslandError::NormalNotConverged {
            report: result.report,
            last_source_refusal: problem.last_refusal.into_inner(),
        });
    }
    let mut residual = vec![0.0; island.patches.len()];
    problem.residual(&solver.x, &mut residual);
    if problem.exhausted.get() {
        return Err(RigidHertzMindlinIslandError::NormalEvaluationBudget {
            evaluations: problem.evaluations.get(),
        });
    }
    let mut maximum: f64 = 0.0;
    for value in residual {
        audit(
            value,
            island
                .controls
                .native
                .multiple
                .local
                .displacement_tolerance_m,
            "original full normal load line",
        )?;
        maximum = maximum.max(value.abs());
    }
    poll(cancelled)?;
    Ok(NormalSolution {
        indentations: solver.x,
        paths: problem.cache.into_inner().paths,
        residual_m: maximum,
        iterations: result.report.iterations,
        evaluations: problem.evaluations.get(),
    })
}
