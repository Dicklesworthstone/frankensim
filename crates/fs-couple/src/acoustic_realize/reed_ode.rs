//! Complete vibrating surfaces on the existing reed and duct owners.
//!
//! Quasistatic and massive reeds close against the complete body's interval
//! pressure under held volume flow. Quadratic bodies use exact propagation;
//! nonlinear bodies retain their original storage through an implicit discrete
//! gradient while the complete stiff linear part evolves exponentially.
use super::*;

mod nonlinear;
use nonlinear::realize_nonlinear;

#[allow(clippy::too_many_arguments)]
pub(super) fn realize_reed_ode(
    sections: &[AcousticSection],
    outlet: ChainOutlet,
    taps: &[AcousticTap],
    reed: BeatingReed,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<LocallyReactingWall>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    if !crate::reed_bore::reed_parameters_valid(reed) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "reed parameters must be physical and finite",
        });
    }
    let inlets = if plates.is_empty() { 1 } else { 2 };
    let line = ode_bore_chain(sections, outlet, inlets, taps, gas, wall)?;
    let (bore, state, observation) = compose_plate_chain(plates, line, 1, gas, sample_rate_hz)?;
    // A bank keeps its independent point-force port at zero. Only the fluid
    // inlet supplies volume flow; surface reaction uses the original area map.
    let source = usize::from(observation.is_some());
    let dt = 1.0 / f64::from(sample_rate_hz);
    let linear = prepare_render_linear(&bore)?;
    if reed.mass_kg == 0.0 && linear.is_exact_quadratic() {
        return realize_quasistatic(
            bore,
            state,
            observation,
            source,
            reed,
            gas,
            listener_m,
            dt,
            n,
        );
    }
    if reed.mass_kg > 0.0 && linear.is_exact_quadratic() {
        return realize_massive(
            bore,
            state,
            observation,
            source,
            reed,
            gas,
            listener_m,
            dt,
            n,
        );
    }
    realize_nonlinear(
        bore,
        state,
        observation,
        source,
        reed,
        gas,
        listener_m,
        dt,
        n,
    )
}

fn prepare_reed_load(
    reed: BeatingReed,
    density: f64,
    dt: f64,
) -> Result<crate::reed_bore::loaded::PreparedReedLoad, AcousticRealizeError> {
    let (_, damping) = reed_structural(reed);
    // Retain the existing reduced lay and its loss calibration for both body
    // propagators. These coefficients do not become measured specimen data.
    let stiffness = 1.0e7 * reed.width_m;
    let chi = damping / (stiffness * reed.rest_opening_m * reed.rest_opening_m).max(1.0e-18);
    let lay = slit_lay(stiffness, 2.0)
        .and_then(|obstacle| obstacle.with_internal_loss(chi))
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
    crate::reed_bore::loaded::PreparedReedLoad::new(reed, density, dt, lay)
}

/// The exact body remains a held-flow owner: its port work is exactly
/// `dt * mean_pressure * flow`. Using a held incoming wave here instead would
/// leave the within-window pressure variance unmatched by the reed's work.
/// Nonlinear pressure/contact forces are held on the original exact reed
/// mechanics; the aperture is evaluated at the common midpoint configuration.
#[allow(clippy::too_many_arguments)]
fn realize_massive(
    system: fs_phs::PortHamiltonian,
    mut state: Vec<f64>,
    observation: Option<PlateBankObservation>,
    source: usize,
    reed: BeatingReed,
    gas: &GasState,
    listener_m: f64,
    dt: f64,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let budget = fs_phs::linear_exact::LinearSplitBudget {
        max_modes: system.state_dim(),
        max_omega_dt: core::f64::consts::PI,
    };
    let mut propagator = fs_phs::linear_exact::PreparedExactZoh::new_owned(system, dt, budget)
        .map_err(render_phs_error)?;
    let load = prepare_reed_load(reed, gas.density, dt)?;
    let mut input = vec![0.0; propagator.system().port_dim()];
    let mut next_state = vec![0.0; state.len()];
    let mut average_output = vec![0.0; input.len()];
    let (mut opening, mut velocity) = (reed.rest_opening_m, 0.0);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (free, impedance) = propagator
            .average_relation(&state, &input, source)
            .map_err(render_phs_error)?;
        let next_reed = load.step(
            opening,
            velocity,
            blowing_envelope(reed, i as f64 * dt),
            free,
            impedance,
        )?;
        input[source] = next_reed.volume_flow_m3_s;
        let body = propagator
            .step_into(&state, &input, &mut next_state, &mut average_output)
            .map_err(render_phs_error)?;
        // Check the actual accepted port, not only the condensed prediction.
        // Every other external input is zero; the body work must equal the
        // reed's independently evaluated bore work on this same interval.
        let pressure_scale = free.abs()
            + (impedance * input[source]).abs()
            + average_output[source].abs()
            + next_reed.pressure_pa.abs();
        let work_scale = body.supplied.abs()
            + next_reed.account.bore_work_j.abs()
            + dt * pressure_scale * input[source].abs();
        if (average_output[source] - next_reed.pressure_pa).abs()
            > 2.0e-11 * pressure_scale.max(f64::MIN_POSITIVE)
            || (body.supplied - next_reed.account.bore_work_j).abs()
                > 2.0e-11 * work_scale.max(f64::MIN_POSITIVE)
        {
            return Err(AcousticRealizeError::Reed {
                what: "coupled reed and complete body did not agree on pressure work",
            });
        }
        let pressure = average_output[source]
            + plate_chain_radiation(
                observation.as_ref(),
                &next_state,
                &state,
                gas,
                listener_m,
                dt,
            )?;
        if !pressure.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "coupled reed and complete body pressure left the finite set",
            });
        }
        // The candidate reed and body states are published together, after
        // their physical accounts and the existing observer have succeeded.
        core::mem::swap(&mut state, &mut next_state);
        opening = next_reed.opening_m;
        velocity = next_reed.velocity_m_s;
        out.push(pressure);
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn realize_quasistatic(
    system: fs_phs::PortHamiltonian,
    mut state: Vec<f64>,
    observation: Option<PlateBankObservation>,
    source: usize,
    reed: BeatingReed,
    gas: &GasState,
    listener_m: f64,
    dt: f64,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let budget = fs_phs::linear_exact::LinearSplitBudget {
        max_modes: system.state_dim(),
        max_omega_dt: core::f64::consts::PI,
    };
    let mut propagator = fs_phs::linear_exact::PreparedExactZoh::new_owned(system, dt, budget)
        .map_err(render_phs_error)?;
    let mut input = vec![0.0; propagator.system().port_dim()];
    let mut next_state = vec![0.0; state.len()];
    let mut average_output = vec![0.0; input.len()];
    let mut previous_flow = 0.0;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        // The same held-flow boundary as the massive reed makes body work
        // exactly dt*p_bar*q, without an unmatched Robin pressure variance.
        let (free, impedance) = propagator
            .average_relation(&state, &input, source)
            .map_err(render_phs_error)?;
        let valve = crate::reed_bore::loaded::solve_quasistatic_load(
            reed,
            gas.density,
            dt,
            blowing_envelope(reed, i as f64 * dt),
            free,
            impedance,
            previous_flow,
        )?;
        input[source] = valve.volume_flow_m3_s;
        let body = propagator
            .step_into(&state, &input, &mut next_state, &mut average_output)
            .map_err(render_phs_error)?;
        let pressure_scale = free.abs()
            + (impedance * input[source]).abs()
            + average_output[source].abs()
            + valve.pressure_pa.abs();
        let work_scale = body.supplied.abs()
            + valve.bore_work_j.abs()
            + dt * pressure_scale * input[source].abs();
        if (average_output[source] - valve.pressure_pa).abs()
            > 2.0e-11 * pressure_scale.max(f64::MIN_POSITIVE)
            || (body.supplied - valve.bore_work_j).abs()
                > 2.0e-11 * work_scale.max(f64::MIN_POSITIVE)
        {
            return Err(AcousticRealizeError::Reed {
                what: "quasistatic valve and complete body did not agree on pressure work",
            });
        }
        let pressure = average_output[source]
            + plate_chain_radiation(
                observation.as_ref(),
                &next_state,
                &state,
                gas,
                listener_m,
                dt,
            )?;
        if !pressure.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "ODE reed pressure left the finite set",
            });
        }
        previous_flow = valve.volume_flow_m3_s;
        core::mem::swap(&mut state, &mut next_state);
        out.push(pressure);
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
