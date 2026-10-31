//! Shared fixed-stack Dormand--Prince integration of physical radial models.
//!
//! The two-state RP/Keller--Miksis adapters and interacting clusters use the
//! same tableau, componentwise error rule, rejection policy and work budget.

use super::{
    BubbleNumerics, BubbleState, BubbleStep, CavitationError, RadialDynamics, finite, poll,
    positive,
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
    let mut rhs = |time_s, values: &[f64; 2]| -> Result<[f64; 2], CavitationError> {
        let current = BubbleState {
            time_s,
            radius_m: values[0],
            wall_velocity_m_s: values[1],
        };
        let (_, acceleration) = model.dynamics(&current, pressure(time_s))?;
        Ok([values[1], acceleration])
    };
    let result = advance_array(
        state.time_s,
        [state.radius_m, state.wall_velocity_m_s],
        end_time_s,
        numerics,
        [
            numerics.radius_absolute_tolerance_m,
            numerics.velocity_absolute_tolerance_m_s,
        ],
        &mut rhs,
        cancelled,
        &refinable,
    )?;
    let current = BubbleState {
        time_s: result.time_s,
        radius_m: result.values[0],
        wall_velocity_m_s: result.values[1],
    };
    poll(cancelled)?;
    let observation = model.observe(&current, pressure(current.time_s))?;
    poll(cancelled)?;
    *state = current;
    Ok(BubbleStep {
        observation,
        accepted_steps: result.accepted_steps,
        rejected_steps: result.rejected_steps,
        rhs_evaluations: result.rhs_evaluations,
        maximum_normalized_local_error: result.maximum_normalized_local_error,
    })
}

/// Complete candidate only; each physical owner controls final publication.
pub(super) struct ArrayStep<const N: usize> {
    pub time_s: f64,
    pub values: [f64; N],
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub rhs_evaluations: usize,
    pub maximum_normalized_local_error: f64,
}

/// One fixed-stack numerical owner for individual and jointly coupled states.
/// Extra coordinates supply their own physical absolute error scales.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance_array<const N: usize, Error, Rhs, Cancel, Refine>(
    initial_time_s: f64,
    initial: [f64; N],
    end_time_s: f64,
    numerics: &BubbleNumerics,
    absolute: [f64; N],
    rhs: &mut Rhs,
    cancelled: &mut Cancel,
    is_refinable: &Refine,
) -> Result<ArrayStep<N>, Error>
where
    Error: From<CavitationError>,
    Rhs: FnMut(f64, &[f64; N]) -> Result<[f64; N], Error>,
    Cancel: FnMut() -> bool,
    Refine: Fn(&Error) -> bool,
{
    poll(cancelled).map_err(Error::from)?;
    numerics.validate().map_err(Error::from)?;
    if !end_time_s.is_finite()
        || end_time_s < initial_time_s
        || !(end_time_s - initial_time_s).is_finite()
    {
        return Err(CavitationError::InvalidInput("forward time interval").into());
    }
    if N == 0 || absolute.into_iter().any(|value| !positive(value)) {
        return Err(CavitationError::InvalidInput("component integration error scales").into());
    }
    let mut current_time = initial_time_s;
    let mut current = initial;
    let mut rhs_evaluations = 0;
    let mut derivative = evaluate(rhs, current_time, &current, cancelled, &mut rhs_evaluations)?;
    let mut step = numerics.initial_step_s;
    let mut accepted_steps = 0;
    let mut rejected_steps = 0;
    let mut maximum_normalized_local_error: f64 = 0.0;
    while current_time < end_time_s {
        if accepted_steps + rejected_steps == numerics.maximum_attempts {
            return Err(CavitationError::StepBudgetExceeded {
                attempts: numerics.maximum_attempts,
                accepted_steps,
            }
            .into());
        }
        poll(cancelled).map_err(Error::from)?;
        let remaining = end_time_s - current_time;
        step = step.min(remaining);
        let next_time = if step == remaining {
            end_time_s
        } else {
            current_time + step
        };
        if next_time <= current_time || next_time > end_time_s {
            return Err(CavitationError::StepUnderflow.into());
        }
        step = next_time - current_time;
        let can_refine = next_time > current_time + numerics.minimum_step_s;
        let trial = trial(
            rhs,
            current_time,
            &current,
            derivative,
            next_time,
            step,
            numerics.relative_tolerance,
            &absolute,
            cancelled,
            &mut rhs_evaluations,
        );
        let (next, endpoint_derivative, error) = match trial {
            Ok(trial) => trial,
            Err(error) if is_refinable(&error) && can_refine => {
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
                }
                .into());
            }
        } else {
            current = next;
            current_time = next_time;
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
    Ok(ArrayStep {
        time_s: current_time,
        values: current,
        accepted_steps,
        rejected_steps,
        rhs_evaluations,
        maximum_normalized_local_error,
    })
}

fn evaluate<const N: usize, Error, Rhs, Cancel>(
    rhs: &mut Rhs,
    time_s: f64,
    values: &[f64; N],
    cancelled: &mut Cancel,
    evaluations: &mut usize,
) -> Result<[f64; N], Error>
where
    Error: From<CavitationError>,
    Rhs: FnMut(f64, &[f64; N]) -> Result<[f64; N], Error>,
    Cancel: FnMut() -> bool,
{
    poll(cancelled).map_err(Error::from)?;
    *evaluations += 1;
    rhs(time_s, values)
}

#[allow(clippy::too_many_arguments)]
fn trial<const N: usize, Error, Rhs, Cancel>(
    rhs: &mut Rhs,
    time_s: f64,
    initial: &[f64; N],
    first_derivative: [f64; N],
    next_time: f64,
    step: f64,
    relative_tolerance: f64,
    absolute: &[f64; N],
    cancelled: &mut Cancel,
    evaluations: &mut usize,
) -> Result<([f64; N], [f64; N], f64), Error>
where
    Error: From<CavitationError>,
    Rhs: FnMut(f64, &[f64; N]) -> Result<[f64; N], Error>,
    Cancel: FnMut() -> bool,
{
    let mut stages = [[0.0; N]; 7];
    stages[0] = first_derivative;
    let mut values = *initial;
    for i in 0..6 {
        for component in 0..N {
            let mut sum = 0.0;
            for (j, stage) in stages.iter().enumerate().take(i + 1) {
                sum += A[i][j] * stage[component];
            }
            values[component] = initial[component] + step * sum;
        }
        let stage_time = if C[i] == 1.0 {
            next_time
        } else {
            time_s + C[i] * step
        };
        stages[i + 1] = evaluate(rhs, stage_time, &values, cancelled, evaluations)?;
    }
    let mut error: f64 = 0.0;
    for component in 0..N {
        let mut difference = 0.0;
        for (stage, weight) in stages.iter().zip(ERROR_WEIGHT) {
            difference += weight * stage[component];
        }
        let scale = finite(
            absolute[component]
                + relative_tolerance * initial[component].abs().max(values[component].abs()),
            "local error scale",
        )
        .map_err(Error::from)?;
        error = error.max(
            finite((step * difference).abs() / scale, "embedded local error")
                .map_err(Error::from)?,
        );
    }
    Ok((values, stages[6], error))
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

pub(super) fn refinable(error: &CavitationError) -> bool {
    matches!(
        error,
        CavitationError::RadiusOutsideDomain { .. }
            | CavitationError::WallMachExceeded { .. }
            | CavitationError::PressureRatioExceeded { .. }
            | CavitationError::NonFiniteDerived(_)
    )
}
