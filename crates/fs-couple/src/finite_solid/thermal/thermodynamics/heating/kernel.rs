//! One bounded source-temperature solve for the actual constrained inventory.
//! Both callers retain their native mechanical potential and immutable sources.

use super::{
    FiniteThermoelasticHeatError, FiniteThermoelasticHeatSolve, Result, finite, invalid, poll,
};
use fs_thermochem::solution::caloric::temperature::{self, TemperatureSample};

pub(super) trait Endpoint: Clone {
    fn temperature_k(&self) -> f64;
    fn temperature_range_k(&self) -> [f64; 2];
    fn heat_capacity_j_per_k(&self) -> f64;
    fn at_sample_temperature(&self, temperature: f64) -> Result<Self>;
    fn internal_energy_change_to(&self, other: &Self) -> Result<f64>;
    fn entropy_change_to(&self, other: &Self) -> Result<f64>;
    fn mechanical_entropy_scale(&self, other: &Self) -> f64;
}

#[derive(Clone, Debug)]
pub(super) struct Balance<S> {
    pub(super) final_state: S,
    pub(super) heat_j: f64,
    pub(super) energy_change_j: f64,
    pub(super) entropy_change_j_per_k: f64,
    pub(super) temperature_bracket_k: [f64; 2],
    pub(super) iterations: usize,
    pub(super) state_evaluations: usize,
}

pub(super) fn entropy_allowance<S: Endpoint>(
    initial: &S,
    final_state: &S,
    controls: FiniteThermoelasticHeatSolve,
    thermal_entropy_scale_j_per_k: f64,
) -> f64 {
    controls.entropy_tolerance_j_per_k
        + controls.thermal.energy_tolerance_j
            / initial.temperature_k().min(final_state.temperature_k())
        + 64.0
            * f64::EPSILON
            * (thermal_entropy_scale_j_per_k + initial.mechanical_entropy_scale(final_state))
}

pub(super) fn solve<S: Endpoint>(
    initial: &S,
    controls: FiniteThermoelasticHeatSolve,
    heat_at_temperature: impl Fn(f64) -> Result<f64>,
    heat_temperature_derivative_j_per_k: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Balance<S>> {
    poll(&mut cancelled)?;
    controls.thermal.validate()?;
    let heat_j = heat_at_temperature(initial.temperature_k())?;
    let [lo, hi] = controls.thermal.bracket;
    let range = initial.temperature_range_k();
    if !heat_j.is_finite()
        || !heat_temperature_derivative_j_per_k.is_finite()
        || heat_temperature_derivative_j_per_k > 0.0
        || lo < range[0]
        || hi > range[1]
        || initial.temperature_k() < lo
        || initial.temperature_k() > hi
        || !controls.maximum_temperature_step_k.is_finite()
        || controls.maximum_temperature_step_k <= 0.0
        || !(1..=4096).contains(&controls.maximum_state_evaluations)
        || !controls.entropy_tolerance_j_per_k.is_finite()
        || controls.entropy_tolerance_j_per_k < 0.0
    {
        return Err(invalid(
            "finite heat needs an original bracket, temperature spacing and finite work/tolerance budgets",
        ));
    }
    if heat_j == 0.0 {
        return Ok(Balance {
            final_state: initial.clone(),
            heat_j,
            energy_change_j: 0.0,
            entropy_change_j_per_k: 0.0,
            temperature_bracket_k: [initial.temperature_k(); 2],
            iterations: 0,
            state_evaluations: 0,
        });
    }
    let divisions = ((hi - lo) / controls.maximum_temperature_step_k)
        .ceil()
        .max(1.0);
    if divisions + 1.0 > controls.maximum_state_evaluations as f64 {
        return Err(FiniteThermoelasticHeatError::StateBudget { completed: 0 });
    }
    let divisions = divisions as usize;
    let mut states = vec![initial.clone()];
    let mut evaluations = 0;
    for index in 0..=divisions {
        let t = if index == divisions {
            hi
        } else {
            (hi - lo).mul_add(index as f64 / divisions as f64, lo)
        };
        sample(
            initial,
            t,
            &mut states,
            &mut evaluations,
            controls.maximum_state_evaluations,
            &mut cancelled,
        )?;
    }
    let root = temperature::solve(
        controls.thermal,
        |t| -> Result<_> {
            let state = sample(
                initial,
                t,
                &mut states,
                &mut evaluations,
                controls.maximum_state_evaluations,
                &mut cancelled,
            )?;
            let change_j = initial.internal_energy_change_to(&state)?;
            Ok(TemperatureSample {
                state,
                change_j,
                residual_j: finite(change_j - heat_at_temperature(t)?)?,
            })
        },
        |state| Some(state.heat_capacity_j_per_k() - heat_temperature_derivative_j_per_k),
    )?;
    let final_state = root.sample.state;
    let heat_j = heat_at_temperature(final_state.temperature_k())?;
    let change = initial.internal_energy_change_to(&final_state)?;
    let entropy = initial.entropy_change_to(&final_state)?;
    if final_state.temperature_k() == initial.temperature_k()
        || (change - heat_j).abs() > controls.thermal.energy_tolerance_j
    {
        return Err(invalid(
            "nonzero heat must resolve a changed temperature and close actual internal energy",
        ));
    }
    let a = change / initial.temperature_k();
    let b = change / final_state.temperature_k();
    let bounds = [a.min(b), a.max(b)];
    let allowance = entropy_allowance(
        initial,
        &final_state,
        controls,
        bounds[0].abs().max(bounds[1].abs()),
    );
    if !allowance.is_finite() || entropy < bounds[0] - allowance || entropy > bounds[1] + allowance
    {
        return Err(FiniteThermoelasticHeatError::ClausiusBounds {
            entropy_change_j_per_k: entropy,
            bounds_j_per_k: bounds,
        });
    }
    poll(&mut cancelled)?;
    Ok(Balance {
        final_state,
        heat_j,
        energy_change_j: change,
        entropy_change_j_per_k: entropy,
        temperature_bracket_k: root.bracket,
        iterations: root.iterations,
        state_evaluations: evaluations,
    })
}

fn sample<S: Endpoint>(
    initial: &S,
    temperature: f64,
    states: &mut Vec<S>,
    evaluations: &mut usize,
    maximum: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<S> {
    poll(cancelled)?;
    let index = match states.binary_search_by(|state| state.temperature_k().total_cmp(&temperature))
    {
        Ok(index) => return Ok(states[index].clone()),
        Err(index) => index,
    };
    if *evaluations >= maximum {
        return Err(FiniteThermoelasticHeatError::StateBudget {
            completed: *evaluations,
        });
    }
    let state = initial.at_sample_temperature(temperature)?;
    *evaluations += 1;
    for (left, right) in index
        .checked_sub(1)
        .map(|i| (&states[i], &state))
        .into_iter()
        .chain(states.get(index).map(|next| (&state, next)))
    {
        let change = left.internal_energy_change_to(right)?;
        if change <= 0.0 {
            return Err(FiniteThermoelasticHeatError::NonMonotone {
                temperatures_k: [left.temperature_k(), right.temperature_k()],
                energy_change_j: change,
            });
        }
    }
    poll(cancelled)?;
    states.insert(index, state.clone());
    Ok(state)
}
