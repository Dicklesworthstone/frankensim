//! One nonlinear body/reed work interval over the complete retained state.
//!
//! At each body endpoint trial, the shared exponential discrete-gradient owner
//! supplies the affine mean pressure for that trial's held nonlinear effort.
//! The original reed/contact solver closes against that pressure before the
//! body residual is evaluated. No previous pressure or frozen body correction
//! enters the accepted candidate. Time refinement remains necessary: resistance
//! heat and port work use the scheme's discrete effort, not an exact continuous
//! nonlinear trajectory.

use super::*;
use crate::reed_bore::loaded::{LoadedReedStep, QuasistaticLoadStep, solve_quasistatic_load};

enum SourceCandidate {
    Moving(LoadedReedStep),
    Quasistatic(QuasistaticLoadStep),
}

#[allow(clippy::too_many_arguments)]
pub(super) fn realize_nonlinear(
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
    let mut propagator = fs_phs::linear_exact::PreparedNonlinearZoh::new_owned(system, dt, budget)
        .map_err(render_phs_error)?;
    let load = if reed.mass_kg > 0.0 {
        Some(prepare_reed_load(reed, gas.density, dt)?)
    } else {
        None
    };
    let mut input = vec![0.0; propagator.system().port_dim()];
    let mut next_state = vec![0.0; state.len()];
    let mut average_output = vec![0.0; input.len()];
    let (mut opening, mut velocity, mut previous_flow) = (reed.rest_opening_m, 0.0, 0.0);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mouth = blowing_envelope(reed, i as f64 * dt);
        let mut candidate = None;
        let mut source_error = None;
        let mut pressure_scale = 0.0;
        let body = propagator.step_with_port_into(
            &state,
            &mut input,
            source,
            &mut next_state,
            &mut average_output,
            |free, impedance| {
                let trial = if let Some(load) = &load {
                    load.step(opening, velocity, mouth, free, impedance)
                        .map(SourceCandidate::Moving)
                } else {
                    solve_quasistatic_load(
                        reed,
                        gas.density,
                        dt,
                        mouth,
                        free,
                        impedance,
                        previous_flow,
                    )
                    .map(SourceCandidate::Quasistatic)
                };
                match trial {
                    Ok(trial) => {
                        source_error = None;
                        let (pressure, flow) = match &trial {
                            SourceCandidate::Moving(step) => {
                                (step.pressure_pa, step.volume_flow_m3_s)
                            }
                            SourceCandidate::Quasistatic(step) => {
                                (step.pressure_pa, step.volume_flow_m3_s)
                            }
                        };
                        pressure_scale = free.abs() + (impedance * flow).abs() + pressure.abs();
                        candidate = Some(trial);
                        Ok(flow)
                    }
                    Err(error) => {
                        source_error = Some(error);
                        Err(fs_phs::PhsError::Dimension {
                            what: "nonlinear body reed/contact trial refused",
                        })
                    }
                }
            },
        );
        let body = match body {
            Ok(body) => body,
            Err(error) => return Err(source_error.unwrap_or_else(|| render_phs_error(error))),
        };
        // The shared owner repeats the callback at the accepted endpoint before
        // publishing any buffers; this is that exact candidate, not a Jacobian
        // probe. Its original source work must equal the accepted body work.
        let candidate = candidate.ok_or(AcousticRealizeError::Reed {
            what: "nonlinear body did not produce an accepted reed candidate",
        })?;
        let (pressure, work, flow) = match &candidate {
            SourceCandidate::Moving(step) => (
                step.pressure_pa,
                step.account.bore_work_j,
                step.volume_flow_m3_s,
            ),
            SourceCandidate::Quasistatic(step) => {
                (step.pressure_pa, step.bore_work_j, step.volume_flow_m3_s)
            }
        };
        pressure_scale += average_output[source].abs();
        let work_scale = body.supplied.abs() + work.abs() + dt * pressure_scale * flow.abs();
        if (average_output[source] - pressure).abs() > 2e-11 * pressure_scale.max(f64::MIN_POSITIVE)
            || (body.supplied - work).abs() > 2e-11 * work_scale.max(f64::MIN_POSITIVE)
            || input[source] != flow
        {
            return Err(AcousticRealizeError::Reed {
                what: "nonlinear body and reed did not agree on the accepted pressure work",
            });
        }
        let observed = average_output[source]
            + plate_chain_radiation(
                observation.as_ref(),
                &next_state,
                &state,
                gas,
                listener_m,
                dt,
            )?;
        if !observed.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "nonlinear body and reed pressure left the finite set",
            });
        }
        if let SourceCandidate::Moving(step) = candidate {
            opening = step.opening_m;
            velocity = step.velocity_m_s;
        }
        previous_flow = flow;
        core::mem::swap(&mut state, &mut next_state);
        out.push(observed);
    }
    Ok(out)
}
