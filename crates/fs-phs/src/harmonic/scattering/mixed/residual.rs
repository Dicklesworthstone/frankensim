//! Independent original equations for the two higher loaded derivatives.

use super::{
    C64, HarmonicError, HarmonicParameterDirection, ScatteringModel, Sum, TangentState, finite,
    invalid, parameter, poll,
};

/// Audit sE*x_h -(J-R)*e_h +T*x_h-G*u_h plus authored coefficient
/// derivatives on the omega state, and -iE times the appropriate lower state.
/// The closed generator and the imposed wave balance never certify themselves.
#[allow(clippy::too_many_arguments)]
pub(in crate::harmonic::scattering) fn original_higher_residual(
    owner: &ScatteringModel,
    omega: f64,
    higher: &TangentState,
    input: &[C64],
    base_state: &[C64],
    base_effort: &[C64],
    base_input: &[C64],
    forcing_state: &[C64],
    forcing_multiplier: f64,
    direction: HarmonicParameterDirection<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, HarmonicError> {
    let model = &owner.closed;
    let mut relative = 0.0_f64;
    for i in 0..model.n {
        poll(cancelled)?;
        let mut terms = [Sum::default(); 2];
        let mut scale = Sum::default();
        let mut add = |value: C64, coefficient: f64| -> Result<(), HarmonicError> {
            if !finite(value) || !coefficient.is_finite() {
                return Err(invalid("original loaded jet row is not representable"));
            }
            terms[0]
                .product(value.re, coefficient)
                .and_then(|()| terms[1].product(value.im, coefficient))
                .and_then(|()| scale.product(value.abs(), coefficient.abs()))
                .ok_or_else(|| invalid("original loaded jet row is not representable"))
        };
        if i < model.n_diff {
            let x = higher.state[i];
            add(C64::new(x.im, -x.re), omega)?;
            let lower = forcing_state[i];
            add(C64::new(lower.im, -lower.re), forcing_multiplier)?;
        }
        for j in 0..model.n {
            let index = i * model.n + j;
            add(higher.effort[j], -model.structure[index])?;
            add(
                base_effort[j],
                -parameter::at(direction.interconnection, index),
            )?;
            add(base_effort[j], parameter::at(direction.dissipation, index))?;
            add(
                higher.state[j],
                model.resisting_tangent.as_ref().map_or(0.0, |t| t[index]),
            )?;
            add(
                base_state[j],
                parameter::at(direction.resisting_flow, index),
            )?;
        }
        for p in 0..model.m {
            let index = i * model.m + p;
            add(input[p], -model.g[index])?;
            add(base_input[p], -parameter::at(direction.ports, index))?;
        }
        let row = crate::harmonic::input_value(terms, C64::ZERO)?;
        let scale = scale
            .value()
            .ok_or_else(|| invalid("original loaded jet residual scale overflows"))?;
        relative = relative.max(if scale == 0.0 { 0.0 } else { row.abs() / scale });
    }
    if relative > model.budget.relative_residual_tolerance {
        return Err(HarmonicError::Residual {
            relative,
            tolerance: model.budget.relative_residual_tolerance,
        });
    }
    Ok(relative)
}
