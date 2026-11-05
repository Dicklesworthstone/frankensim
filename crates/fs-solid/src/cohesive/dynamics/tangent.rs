//! Original material path tangents composed with fixed mass and elasticity.
use super::*;
use fs_solver::NonlinearProblem;

/// Fallible native Jacobian with the existing source-domain-aware numerical
/// adapter retained for material switches and friction projection boundaries.
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
    coulomb_blocks: Vec<Option<CoulombBlock>>,
    elastic: Vec<f64>,
    valid: bool,
    supported: bool,
}

// Derivatives of the original traction-valued natural map with respect to
// endpoint jump and the two metre-valued reaction coordinates, respectively.
type CoulombBlock = ([[f64; 3]; 2], [[f64; 2]; 2]);

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
        let extra = dimension
            .checked_sub(n)
            .filter(|extra| extra % 2 == 0)
            .ok_or(Error::Invalid("complete Coulomb Jacobian dimension"))?;
        Ok(Self {
            problem,
            numerical: FallibleFiniteDifference::new(dimension, difference_step, evaluate)
                .map_err(difference_error)?,
            failure: RefCell::new(None),
            failed: Cell::new(false),
            cache: RefCell::new(Cache {
                point: vec![0.0; dimension],
                blocks: vec![[[0.0; 3]; 3]; problem.motion.state.points.len()],
                coulomb_blocks: vec![None; extra / 2],
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
        let state = &problem.motion.state;
        let n = state.system.n;
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
                let jump = patch.jump(&x[..n], problem.offsets[index])?;
                let old = state.points[index];
                let friction = problem
                    .motion
                    .friction
                    .as_ref()
                    .and_then(|models| models[index])
                    .filter(|model| model.active());
                let coulomb = problem
                    .motion
                    .coulomb
                    .as_ref()
                    .filter(|binding| binding.needs_pressure(index));
                let tolerance = problem.motion.config.solve.local_energy_tolerance_j_m2;
                let (tangent, pressure, pressure_gradient) = if friction.is_some()
                    || coulomb.is_some()
                {
                    old.trial_over_time_with_damage_weighted_compression_tangent(
                        jump,
                        problem.dt_s,
                        tolerance,
                    )
                    .map(|(_, _, pressure, tangent, gradient)| (tangent, pressure, gradient))
                } else {
                    old.trial_over_time_with_mean_traction_tangent(jump, problem.dt_s, tolerance)
                        .map(|(_, _, tangent)| (tangent, 0.0, Some([0.0; 3])))
                }
                .map_err(|error| Error::Material {
                    patch: patch.id,
                    error,
                })?;
                let Some(mut tangent) = tangent else {
                    cache.point.copy_from_slice(x);
                    cache.valid = true;
                    return Ok(false);
                };
                if let Some(model) = friction {
                    let Some(pressure_gradient) = pressure_gradient else {
                        cache.point.copy_from_slice(x);
                        cache.valid = true;
                        return Ok(false);
                    };
                    let previous = old.jump_m();
                    let slip = [
                        finite(jump[1] - previous[1])?,
                        finite(jump[2] - previous[2])?,
                    ];
                    let Some(shear) =
                        model.response_tangent(pressure, pressure_gradient, slip, problem.dt_s)?
                    else {
                        cache.point.copy_from_slice(x);
                        cache.valid = true;
                        return Ok(false);
                    };
                    for i in 0..2 {
                        for j in 0..3 {
                            tangent[i + 1][j] = finite(tangent[i + 1][j] + shear[i][j])?;
                        }
                    }
                }
                if let Some(binding) = coulomb {
                    let Some(pressure_gradient) = pressure_gradient else {
                        cache.point.copy_from_slice(x);
                        cache.valid = true;
                        return Ok(false);
                    };
                    let previous = old.jump_m();
                    let slip = [
                        finite(jump[1] - previous[1])?,
                        finite(jump[2] - previous[2])?,
                    ];
                    let force = binding.traction(index, &x[n..])?;
                    let Some(block) = binding.natural_map_tangent(
                        patch.id,
                        force,
                        slip,
                        pressure,
                        pressure_gradient,
                        problem.coulomb_mu.borrow()[index],
                    )?
                    else {
                        cache.point.copy_from_slice(x);
                        cache.valid = true;
                        return Ok(false);
                    };
                    let coordinate = binding
                        .coordinate(index)
                        .ok_or(Error::Invalid("missing Coulomb Jacobian coordinates"))?;
                    cache.coulomb_blocks[coordinate / 2] = Some(block);
                }
                cache.blocks[index] = tangent;
            }
            cache.point.copy_from_slice(x);
            cache.valid = true;
            cache.supported = true;
        }

        // Mechanical rows: (2 M / dt^2 + K / 2) du + A B^T(D du + d_lambda).
        // Coulomb rows: A D(lambda - projection(lambda + k slip, mu P)).
        // Reuse the full coupled blocks throughout FGMRES at this exact iterate.
        // Native operators and sparse jump maps retain all cross-coordinate
        // terms; neither a dense global tangent nor SPD assumption is needed.
        let (native_output, constraint_output) = output.split_at_mut(n);
        constraint_output.fill(0.0);
        cache.elastic.fill(f64::NAN);
        state
            .system
            .elastic
            .apply(&direction[..n], &mut cache.elastic);
        problem.poll()?;
        problem.motion.mass.apply(&direction[..n], native_output);
        problem.poll()?;
        for (i, (value, &elastic)) in native_output.iter_mut().zip(&cache.elastic).enumerate() {
            if i % 256 == 0 {
                problem.poll()?;
            }
            finite(*value)?;
            finite(elastic)?;
            *value = finite(2.0 * (*value / problem.dt_s) / problem.dt_s + 0.5 * elastic)?;
        }
        for (index, patch) in state.system.patches.iter().enumerate() {
            problem.poll()?;
            let jump_direction = patch.jump(&direction[..n], [0.0; 3])?;
            let mut traction = [0.0; 3];
            for (i, value) in traction.iter_mut().enumerate() {
                for (j, &increment) in jump_direction.iter().enumerate() {
                    *value = finite(*value + cache.blocks[index][i][j] * increment)?;
                }
            }
            if let Some(binding) = &problem.motion.coulomb
                && let Some(coordinate) = binding.coordinate(index)
            {
                let reaction_direction = binding.traction(index, &direction[n..])?;
                let (jump_block, reaction_block) = cache.coulomb_blocks[coordinate / 2]
                    .ok_or(Error::Invalid("missing Coulomb Jacobian block"))?;
                for i in 0..2 {
                    traction[i + 1] = finite(traction[i + 1] + reaction_direction[i])?;
                    let mut constraint = 0.0;
                    for (j, &increment) in jump_direction.iter().enumerate() {
                        constraint = finite(constraint + jump_block[i][j] * increment)?;
                    }
                    for j in 0..2 {
                        constraint = finite(
                            constraint + reaction_block[i][j] * direction[n + coordinate + j],
                        )?;
                    }
                    constraint_output[coordinate + i] = finite(patch.area_m2 * constraint)?;
                }
            }
            patch.scatter(traction, native_output)?;
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
            self.fail(Error::Invalid(
                "finite complete Jacobian coordinates required",
            ));
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
