//! Original-energy exponential discrete gradients for a declared quadratic part.
//!
//! Write H(x)=xᵀQx/2+U(x), D=J−R. Over one interval hold f=DG U(x0,x1)
//! and evolve x'=D(Qx+f)+Gu through the existing complete exact linear owner.
//! The endpoint solves that SAME trajectory's implicit equation. Consequently
//! ΔH=∫(Qx+f)·dx and the original R dissipates ∫(Qx+f)ᵀR(Qx+f)dt.
//! This is the exponential discrete-gradient construction of Li and Wu (2016),
//! DOI 10.1137/15M1023257. Q, physical coefficients and state order never change.
//!
//! Only the quadratic propagation is exact. The nonlinear effort is discrete;
//! its heat/output are not the pointwise continuous nonlinear trajectory's
//! heat/output. Smooth convergence and nonlinear/contact temporal resolution
//! still require time refinement. A converged energy account is not an accuracy
//! certificate, uniqueness proof, or permission to ignore an excitation band.

use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite,
    positive_heat::PositiveHeat, prepared::PreparedQuadraticFlow,
    prepared_zoh::energetic_scale_into, refuse, split_quadratic,
};
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, QuadraticStorage, StepRecord,
    Storage, discrete_gradient_into_unchecked, quadratic_midpoint_remainder, storage_energy_change,
};
use fs_la::LuWorkspace;
use fs_math::det;

const MAX_NEWTON: usize = 50;
const MAX_BACKTRACK: usize = 24;
const RELATIVE_RESIDUAL: f64 = 1e-13;

type Feedback<'a> = (usize, &'a mut dyn FnMut(f64, f64) -> Result<f64, PhsError>);

struct Remainder<'a> {
    original: &'a dyn Storage,
    q: &'a [f64],
    exact: bool,
}
impl Storage for Remainder<'_> {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if self.exact {
            return 0.0;
        }
        self.original.hamiltonian(x) - 0.5 * super::quadratic(self.q, x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if self.exact {
            out.fill(0.0);
            return;
        }
        self.original.gradient(x, out);
        let n = x.len();
        for i in 0..n {
            out[i] -= (0..n).map(|j| self.q[i * n + j] * x[j]).sum::<f64>();
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        if self.exact {
            return Some(0.0);
        }
        // Keep the original constitutive cancellation-stable numerator. A
        // large energy datum must not enter a difference of residual energies.
        Some(
            self.original.midpoint_energy_remainder(a, b, mid)?
                - quadratic_midpoint_remainder(self.q, a, b, mid),
        )
    }
}

struct Model {
    system: PortHamiltonian,
    q: Vec<f64>,
    exact: bool,
    dt: f64,
    flow: PreparedQuadraticFlow,
    heat: PositiveHeat,
    forward: Vec<f64>,
    inverse: Vec<f64>,
    mean_state: Vec<f64>,
    mean_input: Vec<f64>,
}

struct Evaluation {
    midpoint: Vec<f64>,
    effort: Vec<f64>,
    input: Vec<f64>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: Vec<f64>,
}
struct Scratch {
    x: Vec<f64>,
    trial: Vec<f64>,
    residual: Vec<f64>,
    plus: Vec<f64>,
    minus: Vec<f64>,
    delta: Vec<f64>,
    jacobian: Vec<f64>,
    output: Vec<f64>,
    effort: Vec<f64>,
    derivative: Vec<f64>,
    heat_coordinates: Vec<f64>,
    lu: LuWorkspace,
    evaluation: Evaluation,
}

/// Prepared original-energy nonlinear evolution with an exact quadratic block.
///
/// All original coordinates, including zero-Q directions and material memory,
/// remain present. Cold preparation reuses the existing exponential and positive
/// heat owners. A bounded central-difference Newton solve uses invertible energy
/// coordinates and the existing allocation-free LU workspace. The nonlinear
/// residual, callback, original storage and independent positive loss account
/// must all pass before any caller output is published.
///
/// Warm calls allocate nothing in this owner; storage and feedback callbacks
/// must themselves be read-only and allocation-free for that guarantee. No
/// physical accepted history is stored here. The work cap bounds Newton updates,
/// not nonlinear temporal error. Dense admission covers 2*n+m+1 coordinates.
pub struct PreparedNonlinearZoh {
    model: Model,
    scratch: Scratch,
    max_iterations: usize,
}

impl PreparedNonlinearZoh {
    /// Prepare a finite positive clock and the original explicitly declared Q.
    /// Missing decompositions, unchanged dense limits and nonfinite caches refuse.
    pub fn new_owned(
        system: PortHamiltonian,
        dt: f64,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(refuse(
                "exponential discrete gradient requires a positive finite clock",
            ));
        }
        let (n, m) = (system.n, system.m);
        let held = m
            .checked_add(n)
            .ok_or_else(|| refuse("nonlinear held-force dimension"))?;
        let dimension = n
            .checked_add(held)
            .and_then(|d| d.checked_add(1))
            .ok_or_else(|| refuse("nonlinear complete dimension"))?;
        if dimension > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
            return Err(refuse(
                "nonlinear complete operator exceeds existing exponential capacity",
            ));
        }
        let split = split_quadratic(&system, budget)?;
        let q = split.quadratic_matrix().to_vec();
        let exact = split.is_exact_quadratic();
        // This is the declared quadratic block, not a replacement physical
        // storage. Extra inputs are the ORIGINAL D*f, held only for this solve.
        let mut g = vec![0.0; n * held];
        for i in 0..n {
            g[i * held..i * held + m].copy_from_slice(&system.g[i * m..(i + 1) * m]);
            for j in 0..n {
                g[i * held + m + j] = system.j[i * n + j] - system.r[i * n + j];
            }
        }
        if !finite(&g) {
            return Err(refuse("nonlinear held-force map is nonfinite"));
        }
        let shadow = PortHamiltonian {
            n,
            m: held,
            j: system.j.clone(),
            r: system.r.clone(),
            g,
            storage: Box::new(QuadraticStorage { n, q: q.clone() }),
            resistance_witness: system.resistance_witness.clone(),
        };
        let shadow_split = split_quadratic(&shadow, budget)?;
        let outputs: Vec<Vec<f64>> = (0..m)
            .map(|p| {
                (0..n)
                    .map(|j| (0..n).map(|i| system.g[i * m + p] * q[i * n + j]).sum())
                    .collect()
            })
            .collect();
        let flow = PreparedQuadraticFlow::new(
            &shadow,
            &shadow_split,
            dt,
            &[],
            &outputs,
            DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX,
                max_multiply_terms: u64::MAX,
            },
        )?;
        let mut mean_state = Vec::with_capacity(m * n);
        let mut mean_input = Vec::with_capacity(m * held);
        for p in 0..m {
            let (x, u) = flow.linear_integral_map(p)?;
            mean_state.extend(x.into_iter().map(|v| v / dt));
            mean_input.extend(u.into_iter().map(|v| v / dt));
        }
        if !finite(&mean_state) || !finite(&mean_input) {
            return Err(refuse("nonlinear mean-output relation is nonfinite"));
        }
        let heat = PositiveHeat::new(&system, &split, dt, &system.r)?;
        let mut forward = vec![0.0; n * n];
        let mut inverse = vec![0.0; n * n];
        for i in 0..n {
            let weight = if split.values[i] > 0.0 {
                det::sqrt(split.values[i])
            } else {
                1.0
            };
            for j in 0..n {
                forward[i * n + j] = weight * split.basis[j * n + i];
                inverse[j * n + i] = split.basis[j * n + i] / weight;
            }
        }
        Ok(Self {
            model: Model {
                system,
                q,
                exact,
                dt,
                flow,
                heat,
                forward,
                inverse,
                mean_state,
                mean_input,
            },
            scratch: Scratch {
                x: vec![0.0; n],
                trial: vec![0.0; n],
                residual: vec![0.0; n],
                plus: vec![0.0; n],
                minus: vec![0.0; n],
                delta: vec![0.0; n],
                jacobian: vec![0.0; n * n],
                output: vec![0.0; m],
                effort: vec![0.0; n],
                derivative: vec![0.0; n],
                heat_coordinates: vec![0.0; 2 * n],
                lu: LuWorkspace::new(n).map_err(|_| refuse("nonlinear LU capacity"))?,
                evaluation: Evaluation {
                    midpoint: vec![0.0; n],
                    effort: vec![0.0; n],
                    input: vec![0.0; held],
                    augmented: vec![0.0; dimension],
                    balanced: vec![0.0; n],
                    endpoint: vec![0.0; n],
                    integrals: vec![0.0; m],
                },
            },
            max_iterations: MAX_NEWTON,
        })
    }

    pub fn system(&self) -> &PortHamiltonian {
        &self.model.system
    }
    pub fn duration(&self) -> f64 {
        self.model.dt
    }

    /// Bound nonlinear work; zero admits only an already resolved predictor.
    /// This changes no convergence threshold or physical coefficient.
    pub fn set_max_iterations(&mut self, limit: usize) -> Result<(), PhsError> {
        if limit > MAX_NEWTON {
            return Err(refuse("nonlinear Newton budget exceeds reference cap"));
        }
        self.max_iterations = limit;
        Ok(())
    }

    pub fn step(&mut self, initial: &[f64], input: &[f64]) -> Result<StepRecord, PhsError> {
        let mut x = vec![0.0; self.model.system.n];
        let mut y = vec![0.0; self.model.system.m];
        let record = self.step_into(initial, input, &mut x, &mut y)?;
        Ok(StepRecord {
            x,
            y,
            delta_h: record.delta_h,
            dissipated: record.dissipated,
            supplied: record.supplied,
            newton_iters: record.newton_iters,
            solver_residual: record.solver_residual,
            solver_certificate: record.solver_certificate,
        })
    }

    pub fn step_into(
        &mut self,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        mean_output: &mut [f64],
    ) -> Result<PreparedStepRecord, PhsError> {
        self.step_into_controlled(initial, input, next, mean_output, || false)
            .map_err(uncancelled)
    }

    /// Poll before Newton updates, Jacobian columns, line-search trials and
    /// final publication. Refusal or cancellation leaves both outputs intact.
    pub fn step_into_controlled(
        &mut self,
        initial: &[f64],
        input: &[f64],
        next: &mut [f64],
        mean_output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        admit_outputs(&self.model, next, mean_output)?;
        let record = solve(
            &self.model,
            &mut self.scratch,
            self.max_iterations,
            initial,
            input,
            &mut None,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&self.scratch.x);
        mean_output.copy_from_slice(&self.scratch.output);
        Ok(record)
    }

    /// Close one original input against its SAME-interval mean output on every
    /// nonlinear endpoint trial. `feedback(free,beta)` returns a held input u
    /// for y_bar=free+beta*u; all other external inputs retain their meaning.
    /// The selected input slot is published with state/output only on success.
    /// The final successful feedback call belongs to the accepted candidate;
    /// the callback must not itself commit physical state during trial solves.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_port_into(
        &mut self,
        initial: &[f64],
        input: &mut [f64],
        port: usize,
        next: &mut [f64],
        mean_output: &mut [f64],
        feedback: impl FnMut(f64, f64) -> Result<f64, PhsError>,
    ) -> Result<PreparedStepRecord, PhsError> {
        self.step_with_port_into_controlled(
            initial,
            input,
            port,
            next,
            mean_output,
            feedback,
            || false,
        )
        .map_err(uncancelled)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_with_port_into_controlled(
        &mut self,
        initial: &[f64],
        input: &mut [f64],
        port: usize,
        next: &mut [f64],
        mean_output: &mut [f64],
        mut feedback: impl FnMut(f64, f64) -> Result<f64, PhsError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        if port >= self.model.system.m {
            return Err(PhsError::BadPortPairing.into());
        }
        admit_outputs(&self.model, next, mean_output)?;
        let record = solve(
            &self.model,
            &mut self.scratch,
            self.max_iterations,
            initial,
            input,
            &mut Some((port, &mut feedback)),
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&self.scratch.x);
        mean_output.copy_from_slice(&self.scratch.output);
        input[port] = self.scratch.evaluation.input[port];
        Ok(record)
    }
}

fn uncancelled(error: PreparedStepError) -> PhsError {
    match error {
        PreparedStepError::Solver(error) => error,
        PreparedStepError::Cancelled => refuse("uncancelled nonlinear propagation interrupted"),
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_outputs(model: &Model, next: &[f64], output: &[f64]) -> Result<(), PhsError> {
    if next.len() != model.system.n || output.len() != model.system.m {
        return Err(refuse(
            "nonlinear output buffers require the complete state and ports",
        ));
    }
    Ok(())
}

fn norm(values: &[f64]) -> Result<f64, PhsError> {
    if !finite(values) {
        return Err(refuse("nonlinear candidate is nonfinite"));
    }
    Ok(values.iter().fold(0.0_f64, |a, v| a.max(v.abs())))
}
fn balanced_scale(model: &Model, x: &[f64]) -> f64 {
    let n = model.system.n;
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| model.forward[i * n + j] * x[j])
                .sum::<f64>()
                .abs()
        })
        .fold(f64::MIN_POSITIVE, f64::max)
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    model: &Model,
    e: &mut Evaluation,
    initial: &[f64],
    candidate: &[f64],
    input: &[f64],
    feedback: &mut Option<Feedback<'_>>,
    residual: &mut [f64],
) -> Result<f64, PhsError> {
    let (n, m) = (model.system.n, model.system.m);
    let held = n + m;
    norm(candidate)?;
    discrete_gradient_into_unchecked(
        &Remainder {
            original: model.system.storage.as_ref(),
            q: &model.q,
            exact: model.exact,
        },
        initial,
        candidate,
        &mut e.midpoint,
        &mut e.effort,
    );
    norm(&e.effort)?;
    e.input[..m].copy_from_slice(input);
    e.input[m..].copy_from_slice(&e.effort);
    if let Some((port, callback)) = feedback.as_mut() {
        let state_row = &model.mean_state[*port * n..(*port + 1) * n];
        let input_row = &model.mean_input[*port * held..(*port + 1) * held];
        let free = state_row
            .iter()
            .zip(initial)
            .map(|(a, x)| a * x)
            .sum::<f64>()
            + input_row
                .iter()
                .zip(&e.input)
                .enumerate()
                .filter(|(p, _)| *p != *port)
                .map(|(_, (a, u))| a * u)
                .sum::<f64>()
            + (0..n)
                .map(|i| model.system.g[i * m + *port] * e.effort[i])
                .sum::<f64>();
        if !free.is_finite() {
            return Err(refuse("nonlinear free mean output is nonfinite"));
        }
        e.input[*port] = callback(free, input_row[*port])?;
        if !e.input[*port].is_finite() {
            return Err(refuse("nonlinear feedback returned nonfinite input"));
        }
    }
    model.flow.evaluate_into(
        initial,
        &e.input,
        &mut e.augmented,
        &mut e.balanced,
        &mut e.endpoint,
        &mut e.integrals,
    )?;
    for i in 0..n {
        residual[i] = (0..n)
            .map(|j| model.forward[i * n + j] * (candidate[j] - e.endpoint[j]))
            .sum();
    }
    norm(residual)
}

#[allow(clippy::too_many_arguments)]
fn solve(
    model: &Model,
    s: &mut Scratch,
    max_iterations: usize,
    initial: &[f64],
    input: &[f64],
    feedback: &mut Option<Feedback<'_>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedStepRecord, PreparedStepError> {
    poll(cancelled)?;
    let (n, m) = (model.system.n, model.system.m);
    let selected = feedback.as_ref().map(|(port, _)| *port);
    if initial.len() != n
        || input.len() != m
        || !finite(initial)
        || input
            .iter()
            .enumerate()
            .any(|(p, u)| Some(p) != selected && !u.is_finite())
        || !model.system.hamiltonian(initial).is_finite()
    {
        return Err(
            refuse("nonlinear step requires complete finite state and held sources").into(),
        );
    }
    // The held initial nonlinear effort is only a predictor. Even zero Newton
    // updates require re-evaluation of the implicit endpoint equation.
    s.x.copy_from_slice(initial);
    evaluate(
        model,
        &mut s.evaluation,
        initial,
        &s.x,
        input,
        feedback,
        &mut s.residual,
    )?;
    s.x.copy_from_slice(&s.evaluation.endpoint);
    let initial_scale = balanced_scale(model, initial);
    let mut iterations = 0;
    let solver_residual = loop {
        poll(cancelled)?;
        let rnorm = evaluate(
            model,
            &mut s.evaluation,
            initial,
            &s.x,
            input,
            feedback,
            &mut s.residual,
        )?;
        let scale = initial_scale
            .max(balanced_scale(model, &s.x))
            .max(balanced_scale(model, &s.evaluation.endpoint));
        if !scale.is_finite() {
            return Err(refuse("nonlinear residual scale is nonfinite").into());
        }
        let target = RELATIVE_RESIDUAL * scale;
        if rnorm <= target {
            break rnorm;
        }
        if iterations >= max_iterations {
            return Err(PhsError::NewtonStalled { residual: rnorm }.into());
        }
        for col in 0..n {
            poll(cancelled)?;
            let h = 1e-6 * scale.max(rnorm);
            if !h.is_finite() || h <= 0.0 {
                return Err(refuse("nonlinear Jacobian increment").into());
            }
            let mut increment = h;
            let mut accepted = false;
            for _ in 0..=12 {
                poll(cancelled)?;
                for i in 0..n {
                    s.trial[i] = s.x[i] + increment * model.inverse[i * n + col];
                }
                let plus = evaluate(
                    model,
                    &mut s.evaluation,
                    initial,
                    &s.trial,
                    input,
                    feedback,
                    &mut s.plus,
                );
                for i in 0..n {
                    s.trial[i] = s.x[i] - increment * model.inverse[i * n + col];
                }
                let minus = evaluate(
                    model,
                    &mut s.evaluation,
                    initial,
                    &s.trial,
                    input,
                    feedback,
                    &mut s.minus,
                );
                if plus.is_ok() && minus.is_ok() {
                    accepted = true;
                    break;
                }
                increment *= 0.5;
            }
            if !accepted {
                return Err(
                    refuse("nonlinear Jacobian has no finite two-sided neighborhood").into(),
                );
            }
            for row in 0..n {
                s.jacobian[row * n + col] = (s.plus[row] - s.minus[row]) / (2.0 * increment);
            }
        }
        s.lu.solve_into(&s.jacobian, &s.residual, &mut s.delta)
            .map_err(|_| PhsError::NewtonStalled { residual: rnorm })?;
        let mut fraction = 1.0;
        let mut accepted = false;
        for _ in 0..=MAX_BACKTRACK {
            poll(cancelled)?;
            for i in 0..n {
                let correction = (0..n)
                    .map(|j| model.inverse[i * n + j] * s.delta[j])
                    .sum::<f64>();
                s.trial[i] = s.x[i] - fraction * correction;
            }
            if let Ok(next_norm) = evaluate(
                model,
                &mut s.evaluation,
                initial,
                &s.trial,
                input,
                feedback,
                &mut s.plus,
            ) {
                if next_norm < rnorm || next_norm <= target {
                    s.x.copy_from_slice(&s.trial);
                    accepted = true;
                    break;
                }
            }
            fraction *= 0.5;
        }
        if !accepted {
            return Err(PhsError::NewtonStalled { residual: rnorm }.into());
        }
        iterations += 1;
    };
    // evaluate above is the FINAL callback at the accepted endpoint. Its fbar,
    // held source and exact trajectory integrals stay paired through accounting.
    for i in 0..n {
        s.effort[i] =
            (0..n).map(|j| model.q[i * n + j] * initial[j]).sum::<f64>() + s.evaluation.effort[i];
    }
    let dissipated = model.heat.evaluate_with_effort_into(
        &model.system,
        &s.effort,
        &s.evaluation.input[..m],
        &mut s.derivative,
        &mut s.heat_coordinates,
    )?;
    for p in 0..m {
        s.output[p] = s.evaluation.integrals[p] / model.dt
            + (0..n)
                .map(|i| model.system.g[i * m + p] * s.evaluation.effort[i])
                .sum::<f64>();
    }
    let supplied = model.dt
        * (0..m)
            .map(|p| s.evaluation.input[p] * s.output[p])
            .sum::<f64>();
    let delta_h = storage_energy_change(
        model.system.storage.as_ref(),
        initial,
        &s.x,
        &mut s.evaluation.midpoint,
        &mut s.effort,
    );
    let scale = energetic_scale_into(&model.system, initial, &mut s.effort)
        .max(energetic_scale_into(&model.system, &s.x, &mut s.effort))
        .max(delta_h.abs())
        .max(supplied.abs())
        .max(dissipated.abs());
    if !finite(&s.output)
        || !finite(&[delta_h, supplied, dissipated, scale])
        || !model.system.hamiltonian(&s.x).is_finite()
        || dissipated < 0.0
        || (delta_h + dissipated - supplied).abs() > 1e-12 * scale
    {
        return Err(
            refuse("nonlinear original-storage work and positive heat balance refused").into(),
        );
    }
    Ok(PreparedStepRecord {
        delta_h,
        dissipated,
        supplied,
        newton_iters: iterations,
        solver_residual,
        solver_certificate: None,
    })
}
