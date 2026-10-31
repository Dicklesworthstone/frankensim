//! Shared fixed-stack Dormand--Prince integration of physical radial models.

use super::{
    BubbleNumerics, BubbleState, BubbleStep, CavitationError, RadialDynamics, finite, poll,
};

pub(super) fn advance<Model, Pressure, Cancel>(
    model: &Model,
    state: &mut BubbleState,
    end_time_s: f64,
    numerics: &BubbleNumerics,
    pressure: &Pressure,
    cancelled: &mut Cancel,
) -> Result<BubbleStep, CavitationError>
where
    Model: RadialDynamics,
    Pressure: Fn(f64) -> Model::Pressure,
    Cancel: FnMut() -> bool,
{
    poll(cancelled)?;
    numerics.validate()?;
    if !end_time_s.is_finite()
        || end_time_s < state.time_s
        || !(end_time_s - state.time_s).is_finite()
    {
        return Err(CavitationError::InvalidInput("forward time interval"));
    }
    let mut current = *state;
    let mut rhs_evaluations = 0;
    let mut derivative = rhs(model, current, pressure, cancelled, &mut rhs_evaluations)?;
    let mut step = numerics.initial_step_s;
    let mut accepted_steps = 0;
    let mut rejected_steps = 0;
    let mut maximum_normalized_local_error: f64 = 0.0;
    while current.time_s < end_time_s {
        if accepted_steps + rejected_steps == numerics.maximum_attempts {
            return Err(CavitationError::StepBudgetExceeded {
                attempts: numerics.maximum_attempts,
                accepted_steps,
            });
        }
        poll(cancelled)?;
        let remaining = end_time_s - current.time_s;
        step = step.min(remaining);
        let next_time = if step == remaining {
            end_time_s
        } else {
            current.time_s + step
        };
        if next_time <= current.time_s || next_time > end_time_s {
            return Err(CavitationError::StepUnderflow);
        }
        step = next_time - current.time_s;
        // Comparing represented endpoint times avoids repeatedly retrying an
        // already minimum-size step whose subtraction rounded slightly up.
        let can_refine = next_time > current.time_s + numerics.minimum_step_s;
        let trial = trial(
            model,
            current,
            derivative,
            next_time,
            step,
            numerics,
            pressure,
            cancelled,
            &mut rhs_evaluations,
        );
        let (next, endpoint_derivative, error) = match trial {
            Ok(trial) => trial,
            Err(error) if refinable(&error) && can_refine => {
                rejected_steps += 1;
                step = (0.5 * step).max(numerics.minimum_step_s);
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected_steps += 1;
            if !can_refine {
                return Err(CavitationError::ToleranceNotMet {
                    normalized_error: error,
                    step_s: step,
                });
            }
        } else {
            current = next;
            derivative = endpoint_derivative;
            accepted_steps += 1;
            maximum_normalized_local_error = maximum_normalized_local_error.max(error);
        }
        let factor = if error == 0.0 {
            5.0
        } else {
            (0.9 * error.powf(-0.2)).clamp(0.2, 5.0)
        };
        step = (step * factor).clamp(numerics.minimum_step_s, numerics.maximum_step_s);
    }
    poll(cancelled)?;
    let observation = model.observe(&current, pressure(current.time_s))?;
    poll(cancelled)?;
    *state = current;
    Ok(BubbleStep {
        observation,
        accepted_steps,
        rejected_steps,
        rhs_evaluations,
        maximum_normalized_local_error,
    })
}

fn rhs<Model, Pressure, Cancel>(
    model: &Model,
    state: BubbleState,
    pressure: &Pressure,
    cancelled: &mut Cancel,
    evaluations: &mut usize,
) -> Result<[f64; 2], CavitationError>
where
    Model: RadialDynamics,
    Pressure: Fn(f64) -> Model::Pressure,
    Cancel: FnMut() -> bool,
{
    poll(cancelled)?;
    *evaluations += 1;
    let (_, acceleration) = model.dynamics(&state, pressure(state.time_s))?;
    Ok([state.wall_velocity_m_s, acceleration])
}

#[allow(clippy::too_many_arguments)]
fn trial<Model, Pressure, Cancel>(
    model: &Model,
    state: BubbleState,
    first_derivative: [f64; 2],
    next_time: f64,
    step: f64,
    numerics: &BubbleNumerics,
    pressure: &Pressure,
    cancelled: &mut Cancel,
    evaluations: &mut usize,
) -> Result<(BubbleState, [f64; 2], f64), CavitationError>
where
    Model: RadialDynamics,
    Pressure: Fn(f64) -> Model::Pressure,
    Cancel: FnMut() -> bool,
{
    let mut stages = [[0.0; 2]; 7];
    stages[0] = first_derivative;
    let initial = [state.radius_m, state.wall_velocity_m_s];
    let mut values = initial;
    for i in 0..6 {
        for component in 0..2 {
            let mut sum = 0.0;
            for (j, stage) in stages.iter().enumerate().take(i + 1) {
                sum += A[i][j] * stage[component];
            }
            values[component] = initial[component] + step * sum;
        }
        let stage_time = if C[i] == 1.0 {
            next_time
        } else {
            state.time_s + C[i] * step
        };
        stages[i + 1] = rhs(
            model,
            BubbleState {
                time_s: stage_time,
                radius_m: values[0],
                wall_velocity_m_s: values[1],
            },
            pressure,
            cancelled,
            evaluations,
        )?;
    }
    let absolute = [
        numerics.radius_absolute_tolerance_m,
        numerics.velocity_absolute_tolerance_m_s,
    ];
    let mut error: f64 = 0.0;
    for component in 0..2 {
        let mut difference = 0.0;
        for (stage, weight) in stages.iter().zip(ERROR_WEIGHT) {
            difference += weight * stage[component];
        }
        let scale = finite(
            absolute[component]
                + numerics.relative_tolerance
                    * initial[component].abs().max(values[component].abs()),
            "local error scale",
        )?;
        error = error.max(finite(
            (step * difference).abs() / scale,
            "embedded local error",
        )?);
    }
    Ok((
        BubbleState {
            time_s: next_time,
            radius_m: values[0],
            wall_velocity_m_s: values[1],
        },
        stages[6],
        error,
    ))
}

// Standard Dormand--Prince coefficients. The final A row is the fifth-order
// endpoint; ERROR_WEIGHT is its difference from the embedded fourth-order rule.
const C: [f64; 6] = [0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
const ERROR_WEIGHT: [f64; 7] = [
    35.0 / 384.0 - 5179.0 / 57600.0,
    0.0,
    500.0 / 1113.0 - 7571.0 / 16695.0,
    125.0 / 192.0 - 393.0 / 640.0,
    -2187.0 / 6784.0 + 92097.0 / 339200.0,
    11.0 / 84.0 - 187.0 / 2100.0,
    -1.0 / 40.0,
];

fn refinable(error: &CavitationError) -> bool {
    matches!(
        error,
        CavitationError::RadiusOutsideDomain { .. }
            | CavitationError::WallMachExceeded { .. }
            | CavitationError::PressureRatioExceeded { .. }
            | CavitationError::NonFiniteDerived(_)
    )
}
