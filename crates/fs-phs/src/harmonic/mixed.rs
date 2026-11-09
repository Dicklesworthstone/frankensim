//! Full original-coordinate response jets through one physical factor.
//!
//! E is the fixed diag(I_n_diff,0) chart of this owner. Material and source
//! directions act through J/R/Q/G/T; changing descriptor rank or coordinates
//! is not an authored E tangent. All implicit solves reuse the same LU.

use super::{
    C64, HarmonicError, HarmonicFactor, HarmonicFrequencyResponse, HarmonicModel,
    HarmonicParameterDirection, HarmonicParameterResponse, Sum, finite, invalid, parameter, poll,
    resource_bounds,
};

/// Original response, first frequency tangent, and exact frequency curvature.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicFrequencyCurvatureResponse {
    /// Original physical response and first angular-frequency derivative.
    pub frequency_response: HarmonicFrequencyResponse,
    /// d²x/domega² in original differential/algebraic state order.
    pub state_second_derivative: Vec<C64>,
    /// d²e/domega² with original storage and algebraic effort semantics.
    pub effort_second_derivative: Vec<C64>,
    /// d²y/domega² in original power-conjugate port order.
    pub output_second_derivative: Vec<C64>,
    /// Curvature of original cycle-mean storage [J per (rad/s)²].
    pub mean_stored_energy_second_derivative: f64,
    /// Per-port physical supplied-power curvature [W per (rad/s)²].
    pub port_power_second_derivative: Vec<f64>,
    /// Curvature of aggregate original source work [W per (rad/s)²].
    pub supplied_power_second_derivative: f64,
    /// Independent original R/T work curvature [W per (rad/s)²].
    pub dissipated_power_second_derivative: f64,
    /// Dissipation curvature minus supply curvature.
    pub power_residual_second_derivative: f64,
    /// Original-equation backward residual of the second-frequency solve.
    pub relative_second_derivative_residual: f64,
}

/// Full frequency/parameter jet of the original physical response.
///
/// Parameter derivatives are local signed slopes per one authored real
/// parameter. All source-law terms remain; neither a bias derivative nor a
/// coordinate/topology change is inferred. Backward residuals are not forward
/// error bounds near a pole.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicFrequencyParameterResponse {
    /// Original response and its first two frequency derivatives.
    pub frequency_response: HarmonicFrequencyCurvatureResponse,
    /// The same original response and its physical parameter derivative.
    pub parameter_response: HarmonicParameterResponse,
    /// d²x/(domega dtheta) in original state order.
    pub state_mixed_derivative: Vec<C64>,
    /// W*x_omega_theta + Q_theta*x_omega in original effort order.
    pub effort_mixed_derivative: Vec<C64>,
    /// G^T*e_omega_theta + G_theta^T*e_omega in original port order.
    pub output_mixed_derivative: Vec<C64>,
    /// Mixed derivative of cycle-mean storage [J per rad/s per parameter].
    pub mean_stored_energy_mixed_derivative: f64,
    /// Per-port physical supplied-power mixed derivatives.
    pub port_power_mixed_derivative: Vec<f64>,
    /// Mixed derivative of aggregate original source work.
    pub supplied_power_mixed_derivative: f64,
    /// Independent original R/T work mixed derivative.
    pub dissipated_power_mixed_derivative: f64,
    /// Dissipation mixed derivative minus supply mixed derivative.
    pub power_residual_mixed_derivative: f64,
    /// Original-equation backward residual of the mixed solve.
    pub relative_mixed_derivative_residual: f64,
}

#[derive(Clone, Copy)]
pub(in crate::harmonic) struct SourceJet<'a> {
    pub input: &'a [C64],
    pub frequency: Option<&'a [C64]>,
    pub second_frequency: Option<&'a [C64]>,
    pub parameter: Option<&'a [C64]>,
    pub mixed: Option<&'a [C64]>,
}

pub(in crate::harmonic) struct TangentState {
    pub state: Vec<C64>,
    pub effort: Vec<C64>,
    pub output: Vec<C64>,
    pub relative_residual: f64,
}

pub(in crate::harmonic) fn admit_curvature(model: &HarmonicModel) -> Result<(), HarmonicError> {
    admit(model, 3, 256, 128)
}

pub(in crate::harmonic) fn admit_mixed(model: &HarmonicModel) -> Result<(), HarmonicError> {
    admit(model, 6, 512, 256)
}

fn admit(
    model: &HarmonicModel,
    multiplier: u128,
    state_bytes: u128,
    port_bytes: u128,
) -> Result<(), HarmonicError> {
    let (bytes, work) = resource_bounds(model.n, model.m)?;
    if bytes as u128 + state_bytes * model.n as u128 + port_bytes * model.m as u128
        > model.budget.max_dense_bytes as u128
        || multiplier * work as u128 > model.budget.max_dense_work as u128
    {
        return Err(invalid(
            "insufficient simultaneous harmonic jet payload or work budget",
        ));
    }
    Ok(())
}

fn validate_source(model: &HarmonicModel, source: SourceJet<'_>) -> Result<(), HarmonicError> {
    for values in [
        Some(source.input),
        source.frequency,
        source.second_frequency,
        source.parameter,
        source.mixed,
    ]
    .into_iter()
    .flatten()
    {
        if values.len() != model.m || values.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "each supplied source jet must contain one finite value per original port",
            ));
        }
    }
    Ok(())
}

impl HarmonicFactor<'_> {
    /// Evaluate the original response and its first two frequency derivatives.
    ///
    /// Optional source slices supply du/domega and d²u/domega² in original
    /// port order; None means exactly zero. The shared factor solves
    /// A*x_omegaomega = G*u_omegaomega + 2*i*E*x_omega. E remains the owner's
    /// fixed differential/algebraic mask; every original material-memory state
    /// and algebraic reaction remains in the solve and physical observations.
    ///
    /// Three ordinary work scores and the ordinary payload plus
    /// 256*n+128*ports bytes are admitted before allocation. No perturbed solve,
    /// fit, extra factorization or state mutation occurs. Failed admission,
    /// cancellation or arithmetic/residual checks publish no partial jet.
    pub fn response_with_frequency_curvature(
        &self,
        input: &[C64],
        input_angular_derivative: Option<&[C64]>,
        input_angular_second_derivative: Option<&[C64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicFrequencyCurvatureResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let source = SourceJet {
            input,
            frequency: input_angular_derivative,
            second_frequency: input_angular_second_derivative,
            parameter: None,
            mixed: None,
        };
        validate_source(self.model, source)?;
        admit_curvature(self.model)?;
        let frequency = self.response_with_frequency_derivative(
            input,
            input_angular_derivative,
            &mut cancelled,
        )?;
        let mut rhs = vec![C64::ZERO; self.model.n];
        for (i, value) in rhs.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            let row = &self.model.g[i * self.model.m..(i + 1) * self.model.m];
            let mut terms = super::input_row_sums(row, |p| at(source.second_frequency, p))?;
            if i < self.model.n_diff {
                add_i_state(&mut terms, frequency.state_derivative[i], 2.0)?;
            }
            *value = super::input_value(terms, C64::ZERO)?;
        }
        let (state, residual) = self.solve_rhs(&rhs, &mut cancelled)?;
        let second = recover_tangent(
            self.model,
            state,
            HarmonicParameterDirection::default(),
            &frequency.response.state,
            &frequency.response.effort,
            residual,
            &mut cancelled,
        )?;
        observe_curvature(self.model, source, frequency, second, &mut cancelled)
    }

    /// Evaluate the full original frequency/parameter response jet.
    ///
    /// The optional source slices supply u_omega, u_omegaomega, u_theta and
    /// u_omega_theta, all independently authored in original port units.
    /// J/R/Q/G/T directions include every constitutive or operating-point term.
    /// The fixed descriptor E=diag(I_n_diff,0) has zero parameter derivative.
    /// Five state solves share this factor: original, first/second frequency,
    /// parameter, and mixed. The original response is reused by the parameter
    /// solve, rather than solved a second time.
    ///
    /// A*x_omega_theta = L_theta*x_omega + i*E*x_theta
    ///                    + G_theta*u_omega + G*u_omega_theta.
    /// Changing storage/output maps and all four physical work products are
    /// retained. Signed nonlinear incremental T work is not total biased heat.
    /// This admits six ordinary work scores and ordinary payload plus
    /// 512*n+256*ports bytes, leaving the reusable factor unchanged on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn response_with_frequency_parameter_derivative(
        &self,
        input: &[C64],
        direction: HarmonicParameterDirection<'_>,
        input_angular_derivative: Option<&[C64]>,
        input_angular_second_derivative: Option<&[C64]>,
        input_parameter_derivative: Option<&[C64]>,
        input_angular_parameter_derivative: Option<&[C64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicFrequencyParameterResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let source = SourceJet {
            input,
            frequency: input_angular_derivative,
            second_frequency: input_angular_second_derivative,
            parameter: input_parameter_derivative,
            mixed: input_angular_parameter_derivative,
        };
        validate_source(self.model, source)?;
        direction.validate(self.model, &mut cancelled)?;
        admit_mixed(self.model)?;
        let frequency = self.response_with_frequency_curvature(
            input,
            source.frequency,
            source.second_frequency,
            &mut cancelled,
        )?;
        let parameter = self.response_parameter_from_primal(
            input,
            direction,
            source.parameter,
            frequency.frequency_response.response.clone(),
            &mut cancelled,
        )?;
        let first = &frequency.frequency_response;
        let action = parameter::generator_action(
            self.model,
            &direction,
            &first.state_derivative,
            &first.effort_derivative,
            &mut cancelled,
        )?;
        let mut rhs = action.generator;
        for (i, value) in rhs.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            let range = i * self.model.m..(i + 1) * self.model.m;
            let mut terms = source_mixed_row(self.model, range, source, direction, |_| true)?;
            if i < self.model.n_diff {
                add_i_state(&mut terms, parameter.state_derivative[i], 1.0)?;
            }
            *value = super::input_value(terms, *value)?;
        }
        let (state, residual) = self.solve_rhs(&rhs, &mut cancelled)?;
        let mixed = recover_tangent(
            self.model,
            state,
            direction,
            &first.state_derivative,
            &first.effort_derivative,
            residual,
            &mut cancelled,
        )?;
        observe_mixed(
            self.model,
            source,
            direction,
            frequency,
            parameter,
            mixed,
            &mut cancelled,
        )
    }
}

fn at(values: Option<&[C64]>, index: usize) -> C64 {
    values.map_or(C64::ZERO, |v| v[index])
}

fn add_i_state(terms: &mut [Sum; 2], state: C64, factor: f64) -> Result<(), HarmonicError> {
    terms[0]
        .product(-state.im, factor)
        .and_then(|()| terms[1].product(state.re, factor))
        .ok_or_else(|| invalid("harmonic jet right-hand side is not representable"))
}

fn source_mixed_row(
    model: &HarmonicModel,
    range: core::ops::Range<usize>,
    source: SourceJet<'_>,
    direction: HarmonicParameterDirection<'_>,
    keep: impl Fn(usize) -> bool,
) -> Result<[Sum; 2], HarmonicError> {
    // If u_omega is absent, its G_theta contraction is exactly zero.
    let dg = if source.frequency.is_some() {
        direction.ports.map(|g| &g[range.clone()])
    } else {
        None
    };
    super::input_tangent_sums(
        &model.g[range],
        dg,
        source.frequency.unwrap_or(source.input),
        source.mixed,
        keep,
    )
}

/// Recover original efforts/outputs once, before a loaded owner restores its
/// actual source derivatives. The same maps serve ordinary and closed solves.
pub(in crate::harmonic) fn recover_tangent(
    model: &HarmonicModel,
    state: Vec<C64>,
    direction: HarmonicParameterDirection<'_>,
    base_state: &[C64],
    base_effort: &[C64],
    residual: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TangentState, HarmonicError> {
    let mut effort = state.clone();
    for (i, value) in effort[..model.n_diff].iter_mut().enumerate() {
        poll(cancelled)?;
        *value = super::input_value(
            super::storage_tangent_row(model, i, base_state, &state, direction.storage)?,
            C64::ZERO,
        )?;
    }
    let mut output = vec![C64::ZERO; model.m];
    for (p, value) in output.iter_mut().enumerate() {
        poll(cancelled)?;
        *value = super::port_output_tangent(
            &model.g,
            model.m,
            p,
            &effort,
            direction.ports,
            base_effort,
        )?;
    }
    if effort.iter().chain(&output).any(|&z| !finite(z)) {
        return Err(invalid("harmonic jet field is not representable"));
    }
    poll(cancelled)?;
    Ok(TangentState {
        state,
        effort,
        output,
        relative_residual: residual,
    })
}

/// Original physical curvature from already solved/recovered field jets.
pub(in crate::harmonic) fn observe_curvature(
    model: &HarmonicModel,
    source: SourceJet<'_>,
    frequency: HarmonicFrequencyResponse,
    second: TangentState,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicFrequencyCurvatureResponse, HarmonicError> {
    let original = &frequency.response;
    let mut energy = Sum::default();
    let mut loss = Sum::default();
    for i in 0..model.n {
        poll(cancelled)?;
        if i < model.n_diff {
            let value = super::storage_row(model, i, &original.state)?;
            let first = super::storage_row(model, i, &frequency.state_derivative)?;
            let second_row = super::storage_row(model, i, &second.state)?;
            super::add_stored_work(&mut energy, value, second.state[i])?;
            super::add_stored_work(&mut energy, first, frequency.state_derivative[i])?;
            super::add_stored_work(&mut energy, first, frequency.state_derivative[i])?;
            super::add_stored_work(&mut energy, second_row, original.state[i])?;
        }
        let value = super::resisting_row(model, i, &original.state, &original.effort)?;
        let first = super::resisting_row(
            model,
            i,
            &frequency.state_derivative,
            &frequency.effort_derivative,
        )?;
        let second_row = super::resisting_row(model, i, &second.state, &second.effort)?;
        super::add_dissipative_work(&mut loss, value, second.effort[i])?;
        super::add_dissipative_work(&mut loss, first, frequency.effort_derivative[i])?;
        super::add_dissipative_work(&mut loss, first, frequency.effort_derivative[i])?;
        super::add_dissipative_work(&mut loss, second_row, original.effort[i])?;
    }
    let mut port_power = Vec::with_capacity(model.m);
    for p in 0..model.m {
        poll(cancelled)?;
        let mut power = Sum::default();
        add_peak(&mut power, source.input[p], second.output[p])?;
        add_peak(
            &mut power,
            at(source.frequency, p),
            frequency.output_derivative[p],
        )?;
        add_peak(
            &mut power,
            at(source.frequency, p),
            frequency.output_derivative[p],
        )?;
        add_peak(
            &mut power,
            at(source.second_frequency, p),
            original.output[p],
        )?;
        port_power.push(power_value(power)?);
    }
    let energy = power_value(energy)?;
    let loss = power_value(loss)?;
    let supply = power_value(physical_power_curvature_sum(
        model,
        source,
        &original.effort,
        &frequency.effort_derivative,
        &second.effort,
        |_| true,
        cancelled,
    )?)?;
    let balance = loss - supply;
    if !balance.is_finite() {
        return Err(invalid("harmonic curvature balance is not representable"));
    }
    poll(cancelled)?;
    Ok(HarmonicFrequencyCurvatureResponse {
        frequency_response: frequency,
        state_second_derivative: second.state,
        effort_second_derivative: second.effort,
        output_second_derivative: second.output,
        mean_stored_energy_second_derivative: energy,
        port_power_second_derivative: port_power,
        supplied_power_second_derivative: supply,
        dissipated_power_second_derivative: loss,
        power_residual_second_derivative: balance,
        relative_second_derivative_residual: second.relative_residual,
    })
}

/// Original physical mixed work from already solved/recovered field jets.
pub(in crate::harmonic) fn observe_mixed(
    model: &HarmonicModel,
    source: SourceJet<'_>,
    direction: HarmonicParameterDirection<'_>,
    frequency: HarmonicFrequencyCurvatureResponse,
    parameter: HarmonicParameterResponse,
    mixed: TangentState,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicFrequencyParameterResponse, HarmonicError> {
    let first = &frequency.frequency_response;
    let original = &first.response;
    let mut energy = Sum::default();
    let mut loss = Sum::default();
    for i in 0..model.n {
        poll(cancelled)?;
        if i < model.n_diff {
            let value = super::storage_row(model, i, &original.state)?;
            let angular = super::storage_row(model, i, &first.state_derivative)?;
            let parameter_row = super::storage_tangent_row(
                model,
                i,
                &original.state,
                &parameter.state_derivative,
                direction.storage,
            )?;
            let mixed_row = super::storage_tangent_row(
                model,
                i,
                &first.state_derivative,
                &mixed.state,
                direction.storage,
            )?;
            super::add_stored_work(&mut energy, value, mixed.state[i])?;
            super::add_stored_work(&mut energy, angular, parameter.state_derivative[i])?;
            super::add_stored_work(&mut energy, parameter_row, first.state_derivative[i])?;
            super::add_stored_work(&mut energy, mixed_row, original.state[i])?;
        }
        let value = super::resisting_row(model, i, &original.state, &original.effort)?;
        let angular =
            super::resisting_row(model, i, &first.state_derivative, &first.effort_derivative)?;
        let parameter_row = super::resisting_tangent_row(
            model,
            i,
            &original.state,
            &original.effort,
            &parameter.state_derivative,
            &parameter.effort_derivative,
            direction,
        )?;
        let mixed_row = super::resisting_tangent_row(
            model,
            i,
            &first.state_derivative,
            &first.effort_derivative,
            &mixed.state,
            &mixed.effort,
            direction,
        )?;
        super::add_dissipative_work(&mut loss, value, mixed.effort[i])?;
        super::add_dissipative_work(&mut loss, angular, parameter.effort_derivative[i])?;
        super::add_dissipative_work(&mut loss, parameter_row, first.effort_derivative[i])?;
        super::add_dissipative_work(&mut loss, mixed_row, original.effort[i])?;
    }
    let mut port_power = Vec::with_capacity(model.m);
    for p in 0..model.m {
        poll(cancelled)?;
        let mut power = Sum::default();
        add_peak(&mut power, source.input[p], mixed.output[p])?;
        add_peak(
            &mut power,
            at(source.frequency, p),
            parameter.output_derivative[p],
        )?;
        add_peak(
            &mut power,
            at(source.parameter, p),
            first.output_derivative[p],
        )?;
        add_peak(&mut power, at(source.mixed, p), original.output[p])?;
        port_power.push(power_value(power)?);
    }
    let energy = power_value(energy)?;
    let loss = power_value(loss)?;
    let supply = power_value(physical_power_mixed_sum(
        model,
        source,
        direction,
        &original.effort,
        &first.effort_derivative,
        &parameter.effort_derivative,
        &mixed.effort,
        |_| true,
        cancelled,
    )?)?;
    let balance = loss - supply;
    if !balance.is_finite() {
        return Err(invalid("harmonic mixed balance is not representable"));
    }
    poll(cancelled)?;
    Ok(HarmonicFrequencyParameterResponse {
        frequency_response: frequency,
        parameter_response: parameter,
        state_mixed_derivative: mixed.state,
        effort_mixed_derivative: mixed.effort,
        output_mixed_derivative: mixed.output,
        mean_stored_energy_mixed_derivative: energy,
        port_power_mixed_derivative: port_power,
        supplied_power_mixed_derivative: supply,
        dissipated_power_mixed_derivative: loss,
        power_residual_mixed_derivative: balance,
        relative_mixed_derivative_residual: mixed.relative_residual,
    })
}

/// Unrounded selected-port source curvature. Loaded owners add independent
/// wave fluxes to this sum before comparing them with original material work.
#[allow(clippy::too_many_arguments)]
pub(in crate::harmonic) fn physical_power_curvature_sum(
    model: &HarmonicModel,
    source: SourceJet<'_>,
    effort: &[C64],
    frequency_effort: &[C64],
    second_effort: &[C64],
    keep: impl Fn(usize) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Sum, HarmonicError> {
    let mut power = Sum::default();
    for i in 0..model.n {
        poll(cancelled)?;
        let row = &model.g[i * model.m..(i + 1) * model.m];
        let value =
            super::input_row_sums(row, |p| if keep(p) { source.input[p] } else { C64::ZERO })?;
        let first = super::input_row_sums(row, |p| {
            if keep(p) {
                at(source.frequency, p)
            } else {
                C64::ZERO
            }
        })?;
        let second = super::input_row_sums(row, |p| {
            if keep(p) {
                at(source.second_frequency, p)
            } else {
                C64::ZERO
            }
        })?;
        super::add_dissipative_work(&mut power, value, second_effort[i])?;
        super::add_dissipative_work(&mut power, first, frequency_effort[i])?;
        super::add_dissipative_work(&mut power, first, frequency_effort[i])?;
        super::add_dissipative_work(&mut power, second, effort[i])?;
    }
    Ok(power)
}

/// Unrounded selected-port mixed source work, including both original G maps.
#[allow(clippy::too_many_arguments)]
pub(in crate::harmonic) fn physical_power_mixed_sum(
    model: &HarmonicModel,
    source: SourceJet<'_>,
    direction: HarmonicParameterDirection<'_>,
    effort: &[C64],
    frequency_effort: &[C64],
    parameter_effort: &[C64],
    mixed_effort: &[C64],
    keep: impl Fn(usize) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Sum, HarmonicError> {
    let mut power = Sum::default();
    for i in 0..model.n {
        poll(cancelled)?;
        let range = i * model.m..(i + 1) * model.m;
        let row = &model.g[range.clone()];
        let value =
            super::input_row_sums(row, |p| if keep(p) { source.input[p] } else { C64::ZERO })?;
        let angular = super::input_row_sums(row, |p| {
            if keep(p) {
                at(source.frequency, p)
            } else {
                C64::ZERO
            }
        })?;
        let parameter = super::input_tangent_sums(
            row,
            direction.ports.map(|dg| &dg[range.clone()]),
            source.input,
            source.parameter,
            &keep,
        )?;
        let mixed = source_mixed_row(model, range, source, direction, &keep)?;
        super::add_dissipative_work(&mut power, value, mixed_effort[i])?;
        super::add_dissipative_work(&mut power, angular, parameter_effort[i])?;
        super::add_dissipative_work(&mut power, parameter, frequency_effort[i])?;
        super::add_dissipative_work(&mut power, mixed, effort[i])?;
    }
    Ok(power)
}

fn add_peak(power: &mut Sum, left: C64, right: C64) -> Result<(), HarmonicError> {
    power
        .half_product(left.re, right.re)
        .and_then(|()| power.half_product(left.im, right.im))
        .ok_or_else(|| invalid("harmonic jet port work is not representable"))
}

fn power_value(power: Sum) -> Result<f64, HarmonicError> {
    power
        .value()
        .ok_or_else(|| invalid("harmonic jet physical observation is not representable"))
}
