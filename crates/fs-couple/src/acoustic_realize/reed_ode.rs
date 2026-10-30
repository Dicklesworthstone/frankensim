//! Complete vibrating surfaces on the existing reed and duct owners.
//!
//! The quasistatic valve closes against the exact current-window scattering
//! response. Massive reeds retain their existing held-force/contact scheme;
//! this module does not claim a monolithic reed/contact/bore energy solve.
use super::*;

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
    let (bore, mut state, observation) = compose_plate_chain(plates, line, 1, gas, sample_rate_hz)?;
    // A bank keeps its independent point-force port at zero. Only the fluid
    // inlet supplies volume flow; surface reaction uses the original area map.
    let source = usize::from(observation.is_some());
    let dt = 1.0 / f64::from(sample_rate_hz);
    let linear = prepare_render_linear(&bore)?;
    if reed.mass_kg == 0.0 && linear.is_exact_quadratic() {
        let radius = sections.first().map_or(0.0, |section| section.radius);
        let area = core::f64::consts::PI * radius * radius;
        let impedance = gas.density * gas.sound_speed / area;
        return realize_quasistatic(
            bore,
            state,
            observation,
            source,
            impedance,
            reed,
            gas,
            listener_m,
            dt,
            n,
        );
    }
    let exact = prepare_chain_exact(&bore, &linear, dt)?;
    let massive = if reed.mass_kg > 0.0 {
        let (stiffness, damping) = reed_structural(reed);
        let system =
            mass_spring_damper(reed.mass_kg, stiffness, damping).map_err(render_phs_error)?;
        let face = reed_pressure_face(reed);
        let k_lay = 1.0e7 * reed.width_m;
        let chi = damping / (k_lay * reed.rest_opening_m * reed.rest_opening_m).max(1.0e-18);
        let lay = slit_lay(k_lay, 2.0)
            .and_then(|obstacle| obstacle.with_internal_loss(chi))
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let linear = prepare_render_linear(&system)?;
        Some((system, linear, face, lay))
    } else {
        None
    };
    let reed_exact = massive
        .as_ref()
        .map(|(system, linear, _, _)| prepare_chain_exact(system, linear, dt))
        .transpose()?
        .flatten();
    let mut reed_state = vec![0.0; 2];
    let mut input = vec![0.0; bore.port_dim()];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mouth = blowing_envelope(reed, i as f64 * dt);
        let pressure = bore.output(&state)[source];
        let drop = mouth - pressure;
        input[source] = if let Some((system, linear, face, lay)) = &massive {
            let opening = reed.rest_opening_m + reed_state[0];
            let velocity = system.output(&reed_state)[0];
            let jet = bernoulli_volume_flow(reed.width_m, opening.max(0.0), drop, gas.density);
            let swept = reed_swept_flow(*face, velocity);
            let contact = slit_contact_force(lay, opening)
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
            let loss = lay
                .dissipative_modal_forces(1, &[opening, velocity], &[velocity])
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?[0];
            let force = -face * drop + contact + loss;
            let next = step_chain_flow(
                system,
                linear,
                reed_exact.as_ref(),
                &reed_state,
                &[force],
                dt,
            )?;
            reed_state = next.x;
            jet + swept
        } else {
            // A nonlinear structural body retains its complete Gonzalez flow
            // and original clock admission. It is never probed into a linear
            // scattering model. This branch keeps the held-flow approximation.
            let opening =
                quasistatic_aperture_opening(reed.rest_opening_m, reed.closing_pressure_pa, drop);
            bernoulli_volume_flow(reed.width_m, opening, drop, gas.density)
        };
        let next = step_chain_flow(&bore, &linear, exact.as_ref(), &state, &input, dt)?;
        let pressure = next.y[source]
            + plate_chain_radiation(observation.as_ref(), &next.x, &state, gas, listener_m, dt)?;
        if !pressure.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "ODE reed pressure left the finite set",
            });
        }
        state = next.x;
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
    impedance: f64,
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
    let mut propagator = fs_phs::linear_exact::PreparedScatteringZoh::new_owned(
        system, source, impedance, dt, budget,
    )
    .map_err(render_phs_error)?;
    let mut input = vec![0.0; propagator.system().port_dim()];
    let mut next_state = vec![0.0; state.len()];
    let mut average_output = vec![0.0; input.len()];
    let mut previous_incoming = 0.0;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        // The shared cache gives the SAME interval's mean-pressure response:
        // p = free + beta*a, hence b = free + (beta-1)*a. The existing valve
        // root therefore enforces its pressure/flow law without a sample lag.
        let (free, beta) = propagator
            .average_relation(&state, &input)
            .map_err(render_phs_error)?;
        let incoming = crate::reed_bore::solve_reed_wave(
            reed,
            gas.density,
            impedance,
            beta - 1.0,
            free,
            blowing_envelope(reed, i as f64 * dt),
            previous_incoming,
            0.0,
        )?;
        input[source] = incoming;
        propagator
            .step_into(&state, &input, &mut next_state, &mut average_output)
            .map_err(render_phs_error)?;
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
        previous_incoming = incoming;
        core::mem::swap(&mut state, &mut next_state);
        out.push(pressure);
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
