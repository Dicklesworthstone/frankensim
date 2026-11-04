//! Physical parameter tangents of the original driven descriptor equations.
//!
//! With L=(J-R)W-T, A=-i omega E-L and A x=G u, a fixed-chart tangent obeys
//! A dx=dG u+G du+dL x. The same factor solves that full coupled equation;
//! algebraic reactions and internal memory are never eliminated or refitted.

use super::{
    C64, HarmonicError, HarmonicFactor, HarmonicModel, HarmonicPortPower, HarmonicResponse, finite,
    invalid, poll, resource_bounds,
};

/// Derivatives per unit of ONE caller-declared real parameter, evaluated at the
/// prepared model. None means exactly zero. Matrices are row-major in original
/// state/port order. The differential/algebraic split and state chart stay fixed.
/// Storage and dissipation tangents are symmetric, interconnection is skew.
/// Signed tangents need NOT be positive semidefinite: decreasing mass, stiffness
/// or loss is a valid local direction. No finite parameter interval is admitted.
#[derive(Clone, Copy, Debug, Default)]
pub struct HarmonicParameterDirection<'a> {
    pub interconnection: Option<&'a [f64]>,
    pub dissipation: Option<&'a [f64]>,
    /// Full n*n dQ with exactly zero algebraic rows/columns.
    pub storage: Option<&'a [f64]>,
    /// n*ports dG, used in BOTH source and observation semantics. Free-mode
    /// eigenvalues do not depend on G, but their observed output does.
    pub ports: Option<&'a [f64]>,
    /// Full n*n TOTAL derivative of T=D_x+D_e W, not an extra R matrix.
    /// No symmetry/PSD assumption. Include every constitutive, Hessian and
    /// operating-point chain-rule term; None holds T fixed independently of Q.
    /// A changing bias is not solved or differentiated automatically here.
    pub resisting_flow: Option<&'a [f64]>,
}

impl HarmonicParameterDirection<'_> {
    pub(super) fn validate(
        &self,
        model: &HarmonicModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HarmonicError> {
        let n = model.n;
        for matrix in [
            self.interconnection,
            self.dissipation,
            self.storage,
            self.resisting_flow,
        ]
        .into_iter()
        .flatten()
        {
            poll(cancelled)?;
            if matrix.len() != n * n || matrix.iter().any(|v| !v.is_finite()) {
                return Err(invalid("complete finite n*n parameter matrices required"));
            }
        }
        if let Some(g) = self.ports {
            if g.len() != n * model.m || g.iter().any(|v| !v.is_finite()) {
                return Err(invalid("complete finite n*ports port derivative required"));
            }
        }
        // Reuse the source owner's symmetry admission, never repair a tangent.
        if let Some(j) = self.interconnection {
            crate::require_skew(j, n, "interconnection derivative")
                .map_err(|_| invalid("interconnection derivative must be skew"))?;
        }
        for matrix in [self.dissipation, self.storage].into_iter().flatten() {
            crate::require_symmetric(matrix, n, "material derivative")
                .map_err(|_| invalid("storage and dissipation derivatives must be symmetric"))?;
        }
        if let Some(q) = self.storage {
            for i in 0..n {
                poll(cancelled)?;
                for j in 0..n {
                    if (i >= model.n_diff || j >= model.n_diff) && q[i * n + j] != 0.0 {
                        return Err(invalid(
                            "algebraic efforts cannot acquire storage through a derivative",
                        ));
                    }
                }
            }
        }
        poll(cancelled)
    }
}

pub(super) fn at(matrix: Option<&[f64]>, index: usize) -> f64 {
    matrix.map_or(0.0, |a| a[index])
}

pub(super) fn accumulate(sum: &mut C64, correction: &mut C64, term: C64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

/// Shared original-coordinate action for forced and free-mode tangents.
pub(super) struct ParameterAction {
    /// dQ x, padded by exact zero on algebraic rows.
    pub(super) storage: Vec<C64>,
    /// dL x=(dJ-dR)e+(J-R)dQ x-dT x.
    pub(super) generator: Vec<C64>,
    /// Sum of original uncancelled term magnitudes, for physical residuals.
    pub(super) envelope: Vec<f64>,
}

/// Apply a previously validated direction to a full original state and effort.
/// No derivative matrix, storage inverse, frequency shift or solve is formed.
pub(super) fn generator_action(
    model: &HarmonicModel,
    direction: &HarmonicParameterDirection<'_>,
    state: &[C64],
    effort: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ParameterAction, HarmonicError> {
    let n = model.n;
    let mut storage = vec![C64::ZERO; n];
    let mut storage_envelope = vec![0.0; n];
    for i in 0..model.n_diff {
        poll(cancelled)?;
        let mut correction = C64::ZERO;
        for j in 0..model.n_diff {
            let q = at(direction.storage, i * n + j);
            accumulate(&mut storage[i], &mut correction, state[j].scale(q));
            storage_envelope[i] += q.abs() * state[j].abs();
        }
    }
    let mut generator = vec![C64::ZERO; n];
    let mut envelope = vec![0.0; n];
    for i in 0..n {
        poll(cancelled)?;
        let mut correction = C64::ZERO;
        for j in 0..n {
            let dj = at(direction.interconnection, i * n + j);
            let dr = at(direction.dissipation, i * n + j);
            let term = effort[j].scale(dj - dr) + storage[j].scale(model.structure[i * n + j]);
            accumulate(&mut generator[i], &mut correction, term);
            envelope[i] += (dj.abs() + dr.abs()) * effort[j].abs()
                + model.structure[i * n + j].abs() * storage_envelope[j];
            if let Some(dt) = direction.resisting_flow {
                accumulate(
                    &mut generator[i],
                    &mut correction,
                    -state[j].scale(dt[i * n + j]),
                );
                envelope[i] += dt[i * n + j].abs() * state[j].abs();
            }
        }
    }
    if generator.iter().chain(&storage).any(|&z| !finite(z))
        || envelope
            .iter()
            .chain(&storage_envelope)
            .any(|x| !x.is_finite())
    {
        return Err(invalid("physical parameter action is not representable"));
    }
    Ok(ParameterAction {
        storage,
        generator,
        envelope,
    })
}

/// Original driven response and its derivative per unit of a real parameter.
///
/// The frequency, state chart and differential/algebraic split stay fixed. The
/// caller supplies complete physical coefficient and optional input tangents;
/// no finite difference, pole approximation or automatic bias derivative is
/// used. Physical derivatives are signed slopes, not additional stored energy
/// or heat. A backward residual is not a near-pole forward-error certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicParameterResponse {
    /// Identical physical response to [`HarmonicFactor::response`].
    pub response: HarmonicResponse,
    /// Differential coordinates and algebraic reactions in original order.
    pub state_derivative: Vec<C64>,
    /// W dx+dQ x; algebraic efforts have no storage derivative.
    pub effort_derivative: Vec<C64>,
    /// G^T de+dG^T e, retaining the changing physical output map.
    pub output_derivative: Vec<C64>,
    /// Derivative of cycle-mean x*Q*x/4, excluding algebraic storage.
    pub mean_stored_energy_derivative: f64,
    /// Supplied mean-power derivatives in original port order.
    pub port_power_derivative: Vec<f64>,
    pub supplied_power_derivative: f64,
    /// Independent derivative of Re(e*R*e+e*T*x)/2. For nonlinear local
    /// models this remains incremental work, not total biased-source heat.
    pub dissipated_power_derivative: f64,
    /// Dissipation derivative minus supply derivative.
    pub power_residual_derivative: f64,
    /// Original-coordinate backward residual of A dx=dG u+G du+dL x.
    pub relative_derivative_residual: f64,
}

impl HarmonicParameterResponse {
    /// (dy/dtheta)/y at one original port: log-amplitude and phase derivatives
    /// per parameter unit. Includes any supplied source derivative. Exactly
    /// zero output has no defined phase/log amplitude and returns None.
    pub fn output_logarithmic_derivative(&self, port: usize) -> Result<Option<C64>, HarmonicError> {
        let output = *self
            .response
            .output
            .get(port)
            .ok_or_else(|| invalid("parameter-slope port is out of range"))?;
        let derivative = *self
            .output_derivative
            .get(port)
            .ok_or_else(|| invalid("parameter-slope derivative port is out of range"))?;
        if !finite(output) || !finite(derivative) {
            return Err(invalid("finite output and parameter derivative required"));
        }
        if output == C64::ZERO {
            return Ok(None);
        }
        let slope = derivative / output;
        if !finite(slope) {
            return Err(invalid("output parameter slope is not representable"));
        }
        Ok(Some(slope))
    }
}

/// Shared admission for one complete original physical parameter tangent.
pub(super) fn admit_response(model: &HarmonicModel) -> Result<(), HarmonicError> {
    let (n, m) = (model.n, model.m);
    let (base_bytes, base_work) = resource_bounds(n, m)?;
    let bytes = n
        .checked_mul(192)
        .and_then(|v| m.checked_mul(64).and_then(|p| v.checked_add(p)))
        .and_then(|v| base_bytes.checked_add(v))
        .ok_or_else(|| invalid("parameter derivative payload bound overflow"))?;
    let additional_work = (n as u64)
        .checked_mul(n as u64 + m as u64)
        .and_then(|v| v.checked_mul(8))
        .ok_or_else(|| invalid("parameter derivative work bound overflow"))?;
    let work = base_work
        .checked_mul(2)
        .and_then(|v| v.checked_add(additional_work))
        .ok_or_else(|| invalid("parameter derivative work bound overflow"))?;
    if bytes > model.budget.max_dense_bytes || work > model.budget.max_dense_work {
        return Err(invalid(
            "insufficient parameter derivative byte or work budget",
        ));
    }
    Ok(())
}

/// Physical state/material observations before the actual source-work account.
/// Loaded boundaries supply their restored input tangent only after dy is known.
pub(super) struct ParameterStateDerivative {
    pub(super) state: Vec<C64>,
    pub(super) effort: Vec<C64>,
    pub(super) output: Vec<C64>,
    pub(super) energy: f64,
    pub(super) dissipation: f64,
    pub(super) relative_residual: f64,
}

impl HarmonicFactor<'_> {
    /// Differentiate the complete driven physical response at fixed frequency.
    ///
    /// `direction` gives dJ,dR,dQ,dG and optionally total dT per parameter unit;
    /// `input_derivative` gives complex du in the original source-port order.
    /// None holds that input fixed. Both dG*u and dG^T*e are retained. The full
    /// tangent reuses this factor's LU, scaling and bounded residual refinement.
    ///
    /// The fixed state chart includes every Maxwell history and algebraic
    /// reaction. Changing material laws/bias must supply their complete matrix
    /// chain rule, including dT; no state-coordinate transformation, topology
    /// change, bias solve or finite parameter interval is inferred.
    ///
    /// Admits twice the ordinary response work score plus 8*(n*n+n*ports), and
    /// the ordinary payload plus 192*n+64*ports bytes before allocation. Failed
    /// admission, cancellation or arithmetic/residual checks publish no partial
    /// result and leave this factor reusable.
    pub fn response_with_parameter_derivative(
        &self,
        input: &[C64],
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicParameterResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let model = self.model;
        let (n, m) = (model.n, model.m);
        if input_derivative
            .is_some_and(|values| values.len() != m || values.iter().any(|&value| !finite(value)))
        {
            return Err(invalid(
                "one finite input derivative per source port required",
            ));
        }
        admit_response(model)?;
        direction.validate(model, &mut cancelled)?;
        let response = self.response(input, &mut cancelled)?;
        let action = generator_action(
            model,
            &direction,
            &response.state,
            &response.effort,
            &mut cancelled,
        )?;
        let mut rhs = action.generator.clone();
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut correction = C64::ZERO;
            for p in 0..m {
                accumulate(
                    &mut rhs[i],
                    &mut correction,
                    input[p].scale(at(direction.ports, i * m + p)),
                );
                if let Some(du) = input_derivative {
                    accumulate(
                        &mut rhs[i],
                        &mut correction,
                        du[p].scale(model.g[i * m + p]),
                    );
                }
            }
        }
        let (state_derivative, relative_derivative_residual) =
            self.solve_rhs(&rhs, &mut cancelled)?;
        let observed = observe_state(
            model,
            &response,
            direction,
            state_derivative,
            &action.storage,
            relative_derivative_residual,
            &mut cancelled,
        )?;
        finish_response(response, observed, input, input_derivative, &mut cancelled)
    }
}

/// The original storage, output and independent material loss chain rules.
/// No boundary feedback is included in R or counted as material heat.
#[allow(clippy::too_many_arguments)]
pub(super) fn observe_state(
    model: &HarmonicModel,
    response: &HarmonicResponse,
    direction: HarmonicParameterDirection<'_>,
    state_derivative: Vec<C64>,
    storage_derivative: &[C64],
    relative_derivative_residual: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ParameterStateDerivative, HarmonicError> {
    let (n, m) = (model.n, model.m);
    let mut effort_derivative = state_derivative.clone();
    let mut energy_derivative = 0.0;
    for i in 0..model.n_diff {
        poll(&mut cancelled)?;
        effort_derivative[i] = storage_derivative[i];
        let mut correction = C64::ZERO;
        for j in 0..model.n_diff {
            accumulate(
                &mut effort_derivative[i],
                &mut correction,
                state_derivative[j].scale(model.q[i * n + j]),
            );
        }
        energy_derivative += 0.5 * half_pair(state_derivative[i], response.effort[i])?
            + 0.5 * half_pair(response.state[i], effort_derivative[i])?;
    }
    let mut output_derivative = vec![C64::ZERO; m];
    for p in 0..m {
        poll(&mut cancelled)?;
        let mut correction = C64::ZERO;
        for i in 0..n {
            accumulate(
                &mut output_derivative[p],
                &mut correction,
                effort_derivative[i].scale(model.g[i * m + p])
                    + response.effort[i].scale(at(direction.ports, i * m + p)),
            );
        }
    }
    let mut loss_derivative = 0.0;
    for i in 0..n {
        poll(&mut cancelled)?;
        let mut resistance = C64::ZERO;
        let mut d_resistance = C64::ZERO;
        let mut correction = C64::ZERO;
        let mut d_correction = C64::ZERO;
        for j in 0..n {
            accumulate(
                &mut resistance,
                &mut correction,
                response.effort[j].scale(model.r[i * n + j]),
            );
            accumulate(
                &mut d_resistance,
                &mut d_correction,
                effort_derivative[j].scale(model.r[i * n + j])
                    + response.effort[j].scale(at(direction.dissipation, i * n + j)),
            );
            if let Some(t) = &model.resisting_tangent {
                accumulate(
                    &mut resistance,
                    &mut correction,
                    response.state[j].scale(t[i * n + j]),
                );
                accumulate(
                    &mut d_resistance,
                    &mut d_correction,
                    state_derivative[j].scale(t[i * n + j]),
                );
            }
            if let Some(dt) = direction.resisting_flow {
                accumulate(
                    &mut d_resistance,
                    &mut d_correction,
                    response.state[j].scale(dt[i * n + j]),
                );
            }
        }
        loss_derivative += half_pair(effort_derivative[i], resistance)?
            + half_pair(response.effort[i], d_resistance)?;
    }
    if effort_derivative
        .iter()
        .chain(&output_derivative)
        .any(|&value| !finite(value))
        || !energy_derivative.is_finite()
        || !loss_derivative.is_finite()
    {
        return Err(invalid("parameter derivative observation overflow"));
    }
    poll(&mut cancelled)?;
    Ok(ParameterStateDerivative {
        state: state_derivative,
        effort: effort_derivative,
        output: output_derivative,
        energy: energy_derivative,
        dissipation: loss_derivative,
        relative_residual: relative_derivative_residual,
    })
}

/// Close the original physical account using actual inputs and their tangents.
pub(super) fn finish_response(
    response: HarmonicResponse,
    observed: ParameterStateDerivative,
    input: &[C64],
    input_derivative: Option<&[C64]>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<HarmonicParameterResponse, HarmonicError> {
    let m = response.output.len();
    let ParameterStateDerivative {
        state: state_derivative,
        effort: effort_derivative,
        output: output_derivative,
        energy: energy_derivative,
        dissipation: loss_derivative,
        relative_residual: relative_derivative_residual,
    } = observed;
    let mut port_power_derivative = Vec::with_capacity(m);
    for p in 0..m {
        poll(&mut cancelled)?;
        let derivative = if let Some(du) = input_derivative {
            half_pair(du[p], response.output[p])?
        } else {
            0.0
        } + half_pair(input[p], output_derivative[p])?;
        port_power_derivative.push(derivative);
    }
    let supply_derivative: f64 = port_power_derivative.iter().sum();
    let balance_derivative = loss_derivative - supply_derivative;
    if effort_derivative
        .iter()
        .chain(&output_derivative)
        .any(|&v| !finite(v))
        || port_power_derivative.iter().any(|v| !v.is_finite())
        || [
            energy_derivative,
            loss_derivative,
            supply_derivative,
            balance_derivative,
        ]
        .iter()
        .any(|v| !v.is_finite())
    {
        return Err(invalid("parameter derivative observation overflow"));
    }
    poll(&mut cancelled)?;
    Ok(HarmonicParameterResponse {
        response,
        state_derivative,
        effort_derivative,
        output_derivative,
        mean_stored_energy_derivative: energy_derivative,
        port_power_derivative,
        supplied_power_derivative: supply_derivative,
        dissipated_power_derivative: loss_derivative,
        power_residual_derivative: balance_derivative,
        relative_derivative_residual,
    })
}

fn half_pair(left: C64, right: C64) -> Result<f64, HarmonicError> {
    Ok(HarmonicPortPower::from_peak(left, right)?.active_w)
}
