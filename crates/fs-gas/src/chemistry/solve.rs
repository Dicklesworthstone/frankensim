use super::*;

#[derive(Default)]
pub(super) struct Statistics {
    pub nonlinear_iterations: usize,
    pub temperature_evaluations: usize,
}

struct StateWorkspace {
    concentrations: Vec<f64>,
    kinetics: KineticWorkspaceV1,
}

/// One allocation per bounded vector/matrix at reactor construction. Stages,
/// rejected trials, temperature roots and Newton iterations reuse that storage.
pub(super) struct Reactor<'a> {
    mechanism: &'a IdealGasKineticsV1,
    controls: IsochoricKineticsSolve,
    volume: f64,
    target_density: f64,
    pub initial: Vec<f64>,
    nu: Vec<f64>,
    state: StateWorkspace,
    trial: Vec<f64>,
    matrix: Vec<f64>,
    residual: Vec<f64>,
    delta: Vec<f64>,
    row_scale: Vec<f64>,
    temperature_extent_derivatives: Vec<f64>,
    pub stats: Statistics,
}
impl<'a> Reactor<'a> {
    pub fn new(
        mechanism: &'a IdealGasKineticsV1,
        controls: IsochoricKineticsSolve,
        volume: Volume,
        target: Energy,
        initial: Vec<f64>,
        kinetics: KineticWorkspaceV1,
    ) -> Result<Self, IsochoricKineticsError> {
        let ns = initial.len();
        let nr = mechanism.reactions().len();
        let mut nu = zeros(ns * nr)?;
        for i in 0..ns {
            for j in 0..nr {
                nu[i * nr + j] = mechanism
                    .stoichiometry()
                    .get(i, j)
                    .expect("admitted stoichiometry") as f64;
            }
        }
        let target_density = target.value() / volume.value();
        if !target_density.is_finite() {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable internal-energy density",
            ));
        }
        Ok(Self {
            mechanism,
            controls,
            volume: volume.value(),
            target_density,
            initial,
            nu,
            state: StateWorkspace {
                concentrations: zeros(ns)?,
                kinetics,
            },
            trial: zeros(nr)?,
            matrix: zeros(nr * nr)?,
            residual: zeros(nr)?,
            delta: zeros(nr)?,
            row_scale: zeros(nr)?,
            temperature_extent_derivatives: zeros(nr)?,
            stats: Statistics::default(),
        })
    }
    pub fn concentrations(&self) -> &[f64] {
        &self.state.concentrations
    }
    pub fn resolve(
        &mut self,
        extents: &[f64],
        hint: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
        resolve(
            self.mechanism,
            self.controls,
            self.volume,
            self.target_density,
            &self.initial,
            &self.nu,
            extents,
            hint,
            &mut self.state,
            &mut self.stats,
            cancelled,
        )
    }

    pub fn implicit_step(
        &mut self,
        base: &[f64],
        step: f64,
        out: &mut [f64],
        hint: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
        if !positive(step) {
            return Err(failed(step, "unrepresentable implicit stage duration"));
        }
        out.copy_from_slice(base);
        let nr = base.len();
        let ns = self.initial.len();
        let mut state = self.resolve(out, hint, cancelled)?;
        for iteration in 0..self.controls.maximum_newton_iterations {
            poll(cancelled)?;
            self.stats.nonlinear_iterations += 1;
            self.mechanism.evaluate(
                state.temperature,
                &self.state.concentrations,
                &mut self.state.kinetics,
                &mut *cancelled,
            )?;
            for j in 0..nr {
                self.residual[j] =
                    (out[j] - base[j]) - step * self.state.kinetics.reaction_rates()[j];
            }
            let norm = residual_norm(
                self.controls,
                self.volume,
                &self.nu,
                &self.state.concentrations,
                &self.residual,
            );
            // A nonzero source must receive at least one Newton correction.
            // Otherwise many small admitted steps could silently freeze slow
            // chemistry merely because each individual increment is small.
            if norm <= 0.02 && (iteration > 0 || self.residual.iter().all(|r| *r == 0.0)) {
                return Ok(state);
            }
            if !norm.is_finite() {
                return Err(failed(step, "non-finite implicit residual"));
            }
            let u = self.state.kinetics.species_internal_energies();
            for k in 0..nr {
                let mut du = Sum::default();
                for i in 0..ns {
                    du.add(u[i] * self.nu[i * nr + k]);
                }
                self.temperature_extent_derivatives[k] =
                    -du.value() / state.heat_capacity_v_density_j_m3_k;
            }
            let jac = self.state.kinetics.rate_jacobian();
            let dt = self.state.kinetics.temperature_derivatives();
            for j in 0..nr {
                for k in 0..nr {
                    let mut derivative = Sum::default();
                    for i in 0..ns {
                        derivative.add(jac[j * ns + i] * self.nu[i * nr + k]);
                    }
                    derivative.add(dt[j] * self.temperature_extent_derivatives[k]);
                    self.matrix[j * nr + k] =
                        if j == k { 1.0 } else { 0.0 } - step * derivative.value();
                }
                self.delta[j] = -self.residual[j];
            }
            solve_linear(
                &mut self.matrix,
                &mut self.delta,
                &mut self.row_scale,
                nr,
                cancelled,
            )?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..self.controls.maximum_backtracks {
                poll(cancelled)?;
                for j in 0..nr {
                    self.trial[j] = out[j] + fraction * self.delta[j];
                }
                let candidate = resolve(
                    self.mechanism,
                    self.controls,
                    self.volume,
                    self.target_density,
                    &self.initial,
                    &self.nu,
                    &self.trial,
                    state.temperature.value(),
                    &mut self.state,
                    &mut self.stats,
                    cancelled,
                );
                match candidate {
                    Err(IsochoricKineticsError::Cancelled) => {
                        return Err(IsochoricKineticsError::Cancelled);
                    }
                    Ok(candidate) => {
                        let rates = self.mechanism.evaluate(
                            candidate.temperature,
                            &self.state.concentrations,
                            &mut self.state.kinetics,
                            &mut *cancelled,
                        );
                        match rates {
                            Err(KineticsErrorV1::Cancelled) => {
                                return Err(IsochoricKineticsError::Cancelled);
                            }
                            Ok(_) => {
                                for j in 0..nr {
                                    self.residual[j] = (self.trial[j] - base[j])
                                        - step * self.state.kinetics.reaction_rates()[j];
                                }
                                let trial_norm = residual_norm(
                                    self.controls,
                                    self.volume,
                                    &self.nu,
                                    &self.state.concentrations,
                                    &self.residual,
                                );
                                if trial_norm.is_finite()
                                    && (trial_norm <= 0.02
                                        || trial_norm < norm * (1.0 - 1e-4 * fraction))
                                {
                                    accepted = Some((candidate, trial_norm));
                                    break;
                                }
                            }
                            Err(_) => {}
                        }
                    }
                    Err(_) => {}
                }
                fraction *= 0.5;
            }
            let Some((candidate, trial_norm)) = accepted else {
                return Err(failed(
                    step,
                    "positivity/energy Newton line search exhausted",
                ));
            };
            out.copy_from_slice(&self.trial);
            state = candidate;
            if trial_norm <= 0.02 {
                return Ok(state);
            }
        }
        Err(failed(step, "implicit Newton iteration budget exhausted"))
    }
}

fn failed(step: f64, reason: &'static str) -> IsochoricKineticsError {
    IsochoricKineticsError::StepFailure {
        elapsed_time_s: 0.0,
        step_s: step,
        reason,
    }
}

fn residual_norm(
    controls: IsochoricKineticsSolve,
    volume: f64,
    nu: &[f64],
    c: &[f64],
    residual: &[f64],
) -> f64 {
    let nr = residual.len();
    let absolute = controls.absolute_amount_tolerance_mol / volume;
    let mut norm = 0.0_f64;
    for j in 0..nr {
        if !residual[j].is_finite() {
            return f64::INFINITY;
        }
        let mut scale = f64::INFINITY;
        for i in 0..c.len() {
            let coefficient = nu[i * nr + j].abs();
            if coefficient > 0.0 {
                scale = scale.min((absolute + controls.relative_tolerance * c[i]) / coefficient);
            }
        }
        let ratio = residual[j].abs() / (scale / nr as f64);
        if !positive(scale) || !ratio.is_finite() {
            return f64::INFINITY;
        }
        norm = norm.max(ratio);
    }
    norm
}

#[allow(clippy::too_many_arguments)]
fn resolve(
    mechanism: &IdealGasKineticsV1,
    controls: IsochoricKineticsSolve,
    volume: f64,
    target_density: f64,
    initial: &[f64],
    nu: &[f64],
    extents: &[f64],
    hint: f64,
    work: &mut StateWorkspace,
    statistics: &mut Statistics,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
    poll(cancelled)?;
    let nr = extents.len();
    for (i, &value) in initial.iter().enumerate() {
        let mut c = Sum::default();
        c.add(value);
        for j in 0..nr {
            c.add(nu[i * nr + j] * extents[j]);
        }
        let c = c.value();
        if !c.is_finite() || c < 0.0 {
            return Err(failed(
                0.0,
                "reaction extent leaves the nonnegative species domain",
            ));
        }
        work.concentrations[i] = c;
    }
    let tolerance = controls.energy_tolerance_j / volume;
    let evaluate = |temperature: f64,
                    work: &mut StateWorkspace,
                    statistics: &mut Statistics,
                    cancelled: &mut _| {
        statistics.temperature_evaluations += 1;
        mechanism
            .evaluate_thermodynamics(
                Temperature::new(temperature),
                &work.concentrations,
                &mut work.kinetics,
                cancelled,
            )
            .map_err(IsochoricKineticsError::from)
    };
    let mut lo = controls.minimum_temperature_k;
    let mut hi = controls.maximum_temperature_k;
    let mut t = hint.clamp(lo, hi);
    let first = evaluate(t, work, statistics, &mut *cancelled)?;
    let mut residual = first.internal_energy_density_j_m3 - target_density;
    if residual.abs() <= tolerance {
        return Ok(first);
    }
    let lower = evaluate(lo, work, statistics, &mut *cancelled)?;
    let upper = evaluate(hi, work, statistics, &mut *cancelled)?;
    if (lower.internal_energy_density_j_m3 - target_density).abs() <= tolerance {
        return evaluate(lo, work, statistics, &mut *cancelled);
    }
    if (upper.internal_energy_density_j_m3 - target_density).abs() <= tolerance {
        return Ok(upper);
    }
    if !(lower.internal_energy_density_j_m3 < target_density
        && target_density < upper.internal_energy_density_j_m3)
    {
        return Err(IsochoricKineticsError::EnergyNotBracketed {
            lower_energy_j: lower.internal_energy_density_j_m3 * volume,
            upper_energy_j: upper.internal_energy_density_j_m3 * volume,
        });
    }
    let mut state = first;
    for iteration in 0..controls.maximum_temperature_iterations {
        poll(cancelled)?;
        if residual < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        let next = t - residual / state.heat_capacity_v_density_j_m3_k;
        t = if next > lo && next < hi {
            next
        } else {
            lo + 0.5 * (hi - lo)
        };
        state = evaluate(t, work, statistics, &mut *cancelled)?;
        residual = state.internal_energy_density_j_m3 - target_density;
        if residual.abs() <= tolerance {
            return Ok(state);
        }
        if t == lo || t == hi {
            return Err(IsochoricKineticsError::TemperatureNoConvergence {
                iterations: iteration + 1,
                residual_j: residual * volume,
            });
        }
    }
    Err(IsochoricKineticsError::TemperatureNoConvergence {
        iterations: controls.maximum_temperature_iterations,
        residual_j: residual * volume,
    })
}

// Purpose-specific bounded dense Newton solve (at most128reaction extents).
// Scaled partial pivoting is deterministic; it refuses unresolved pivots and
// nonfinite arithmetic. The nonlinear residual remains the acceptance owner.
fn solve_linear(
    a: &mut [f64],
    b: &mut [f64],
    scale: &mut [f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), IsochoricKineticsError> {
    for i in 0..n {
        scale[i] = a[i * n..(i + 1) * n]
            .iter()
            .fold(0.0_f64, |s, x| s.max(x.abs()));
        if !positive(scale[i])
            || !b[i].is_finite()
            || a[i * n..(i + 1) * n].iter().any(|x| !x.is_finite())
        {
            return Err(failed(0.0, "non-finite Newton linear system"));
        }
    }
    for k in 0..n {
        poll(cancelled)?;
        let mut pivot = k;
        let mut best = a[k * n + k].abs() / scale[k];
        for i in k + 1..n {
            let value = a[i * n + k].abs() / scale[i];
            if value > best {
                best = value;
                pivot = i;
            }
        }
        if !best.is_finite() || best <= f64::EPSILON * n as f64 {
            return Err(failed(0.0, "numerically singular Newton matrix"));
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            b.swap(k, pivot);
            scale.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    for i in (0..n).rev() {
        let mut value = Sum::default();
        value.add(b[i]);
        for j in i + 1..n {
            value.add(-a[i * n + j] * b[j]);
        }
        b[i] = value.value() / a[i * n + i];
        if !b[i].is_finite() {
            return Err(failed(0.0, "non-finite Newton correction"));
        }
    }
    Ok(())
}
