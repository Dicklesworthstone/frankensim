//! Original material path tangents composed with fixed mass and elasticity.
use super::*;
use fs_solver::NonlinearProblem;

/// Fallible native Jacobian with the existing source-domain-aware numerical
/// adapter retained for material switches and friction laws without a tangent.
pub(super) struct StepJacobian<'s, 'p, 'm, 'a, F> {
    problem: &'s StepProblem<'p, 'm, 'a>,
    numerical: FallibleFiniteDifference<F, Error>,
    failure: RefCell<Option<Error>>,
    failed: Cell<bool>,
    cache: RefCell<Cache>,
}

struct Cache {
    point: Vec<f64>,
    blocks: Vec<[[f64; 3]; 3]>,
    elastic: Vec<f64>,
    valid: bool,
    supported: bool,
}

impl<'s, 'p, 'm, 'a, F> StepJacobian<'s, 'p, 'm, 'a, F>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, Error>,
{
    pub(super) fn new(
        problem: &'s StepProblem<'p, 'm, 'a>,
        dimension: usize,
        difference_step: f64,
        evaluate: F,
    ) -> Result<Self, Error> {
        let n = problem.motion.state.system.n;
        Ok(Self {
            problem,
            numerical: FallibleFiniteDifference::new(dimension, difference_step, evaluate)
                .map_err(difference_error)?,
            failure: RefCell::new(None),
            failed: Cell::new(false),
            cache: RefCell::new(Cache {
                point: vec![0.0; n],
                blocks: vec![[[0.0; 3]; 3]; problem.motion.state.points.len()],
                elastic: vec![0.0; n],
                valid: false,
                supported: false,
            }),
        })
    }

    pub(super) fn has_failed(&self) -> bool {
        self.failed.get() || self.numerical.has_failed()
    }

    pub(super) fn take_failure(&self) -> Option<Error> {
        self.failure
            .borrow_mut()
            .take()
            .or_else(|| self.numerical.take_failure().map(difference_error))
    }

    fn fail(&self, error: Error) {
        if !self.failed.replace(true) {
            *self.failure.borrow_mut() = Some(error);
        }
    }

    fn original_action(
        &self,
        x: &[f64],
        direction: &[f64],
        output: &mut [f64],
    ) -> Result<bool, Error> {
        let problem = self.problem;
        problem.poll()?;
        // These owners add their own damage-weighted pressure and projection
        // derivatives; the cohesive tangent alone does not represent them.
        if problem.motion.coulomb.is_some()
            || problem
                .motion
                .friction
                .as_ref()
                .is_some_and(|models| models.iter().flatten().copied().any(DamagedAreaFriction::active))
        {
            return Ok(false);
        }
        let state = &problem.motion.state;
        let mut cache = self.cache.borrow_mut();
        if cache.valid && cache.point.as_slice() == x && !cache.supported {
            return Ok(false);
        }
        problem.charge(true)?;
        if !cache.valid || cache.point.as_slice() != x {
            cache.valid = false;
            cache.supported = false;
            // Admit the actual Newton point, never a displaced constitutive
            // probe. Tangent directions themselves are not physical motion.
            for (i, ((&u1, &u0), &v0)) in x
                .iter()
                .zip(&state.displacements)
                .zip(&problem.motion.velocity_m_s)
                .enumerate()
            {
                if i % 256 == 0 {
                    problem.poll()?;
                }
                let velocity = finite(2.0 * ((u1 - u0) / problem.dt_s) - v0)?;
                if u1.abs() > problem.motion.config.solve.maximum_displacement_m
                    || velocity.abs() > problem.motion.config.maximum_velocity_m_s
                {
                    return Err(Error::Invalid("trial motion applicability limit"));
                }
            }
            for (index, patch) in state.system.patches.iter().enumerate() {
                problem.poll()?;
                let jump = patch.jump(x, problem.offsets[index])?;
                let (_, _, tangent) = state.points[index]
                    .trial_over_time_with_mean_traction_tangent(
                        jump,
                        problem.dt_s,
                        problem.motion.config.solve.local_energy_tolerance_j_m2,
                    )
                    .map_err(|error| Error::Material {
                        patch: patch.id,
                        error,
                    })?;
                let Some(tangent) = tangent else {
                    cache.point.copy_from_slice(x);
                    cache.valid = true;
                    return Ok(false);
                };
                cache.blocks[index] = tangent;
            }
            cache.point.copy_from_slice(x);
            cache.valid = true;
            cache.supported = true;
        }

        // J v = (2 M / dt^2 + K / 2) v + sum A B^T D(t_mean) B v.
        // Reuse each 3x3 material block throughout FGMRES at this exact iterate.
        // Native operators and sparse jump maps retain all cross-coordinate
        // terms; neither a dense global tangent nor SPD assumption is needed.
        cache.elastic.fill(f64::NAN);
        state.system.elastic.apply(direction, &mut cache.elastic);
        problem.poll()?;
        problem.motion.mass.apply(direction, output);
        problem.poll()?;
        for (i, (value, &elastic)) in output.iter_mut().zip(&cache.elastic).enumerate() {
            if i % 256 == 0 {
                problem.poll()?;
            }
            finite(*value)?;
            finite(elastic)?;
            *value = finite(2.0 * (*value / problem.dt_s) / problem.dt_s + 0.5 * elastic)?;
        }
        for (index, patch) in state.system.patches.iter().enumerate() {
            problem.poll()?;
            let jump_direction = patch.jump(direction, [0.0; 3])?;
            let mut traction = [0.0; 3];
            for (i, value) in traction.iter_mut().enumerate() {
                for (j, &increment) in jump_direction.iter().enumerate() {
                    *value = finite(*value + cache.blocks[index][i][j] * increment)?;
                }
            }
            patch.scatter(traction, output)?;
        }
        problem.poll()?;
        Ok(true)
    }
}

impl<F> NonlinearProblem for StepJacobian<'_, '_, '_, '_, F>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, Error>,
{
    fn dimension(&self) -> usize {
        self.numerical.dimension()
    }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.has_failed() {
            self.numerical.residual(x, output);
        }
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.has_failed() {
            return;
        }
        let n = self.dimension();
        if x.len() != n
            || direction.len() != n
            || output.len() != n
            || x.iter().chain(direction).any(|value| !value.is_finite())
        {
            self.fail(Error::Invalid("finite complete Jacobian coordinates required"));
            return;
        }
        if direction.iter().all(|&value| value == 0.0) {
            output.fill(0.0);
            return;
        }
        match self.original_action(x, direction, output) {
            Ok(true) => {}
            Ok(false) => self.numerical.jacobian_apply(x, direction, output),
            Err(error) => {
                output.fill(f64::NAN);
                self.fail(error);
            }
        }
    }
}
