//! Shared original-body propagation for analytically integrated scalar efforts.
//!
//! The event owners supply the physical reaction impulse and first time moment.
//! This module has no contact law or event selection: it restores the original
//! candidate bodies, retaining their budgets and actual represented motion.

use super::super::{CoupledModalSystem, ModalCouplingError, dot, finite, invalid, poll};
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeModel};
use crate::port_adoption::{PortTransfer, StepAccount};
use fs_exec::CancelGate;

#[derive(Clone, Copy)]
pub(crate) struct EventNativeMotion {
    pub modal_energy_before_j: f64,
    pub modal_energy_j: f64,
    pub tangent_increment_m: f64,
    pub tangent_velocity_m_s: f64,
    pub external_work_j: f64,
    /// Independent integral of the complete scalar reaction on the bodies.
    pub reaction_work_j: f64,
}

pub(crate) fn projected_velocity(
    models: &[ModalAcousticTimeModel],
    column: &[f64],
) -> Result<f64, ModalCouplingError> {
    models
        .iter()
        .flat_map(|m| m.states())
        .zip(column)
        .try_fold(0.0, |v, (state, b)| {
            finite(v + b * state.velocity_m_sqrt_kg_per_s)
        })
}

pub(crate) fn projected_position(
    models: &[ModalAcousticTimeModel],
    column: &[f64],
) -> Result<f64, ModalCouplingError> {
    models
        .iter()
        .flat_map(|m| m.states())
        .zip(column)
        .try_fold(0.0, |x, (state, b)| {
            finite(x + b * state.displacement_m_sqrt_kg)
        })
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn stage_motion(
    network: &mut CoupledModalSystem,
    candidate_states: &mut [ModalAcousticState],
    column: &[f64],
    external: &[f64],
    impulse_n_s: f64,
    time_moment_n_s2: f64,
    gate: &CancelGate,
) -> Result<EventNativeMotion, ModalCouplingError> {
    poll(Some(gate))?;
    let count = network.mode_count();
    if column.len() != count || candidate_states.len() != count || external.len() != count {
        return Err(invalid("event trajectory must cover every original mode"));
    }
    let dt = network.dt;
    let v0 = projected_velocity(&network.models, column)?;
    let acceleration = dot(column, external)?;
    let position_impulse = finite(dt * impulse_n_s - time_moment_n_s2)?;
    let mut modal_energy_before_j = 0.0;
    let mut modal_energy_j = 0.0;
    for (k, old) in network.models.iter().flat_map(|m| m.states()).enumerate() {
        let increment = finite(
            dt * old.velocity_m_sqrt_kg_per_s + 0.5 * dt * dt * external[k]
                - column[k] * position_impulse,
        )?;
        let next = ModalAcousticState {
            displacement_m_sqrt_kg: finite(old.displacement_m_sqrt_kg + increment)?,
            velocity_m_sqrt_kg_per_s: finite(
                old.velocity_m_sqrt_kg_per_s + dt * external[k] - column[k] * impulse_n_s,
            )?,
        };
        network.old_q[k] = old.displacement_m_sqrt_kg;
        // Use actual represented travel, including any loss of significance.
        network.free_delta[k] = finite(next.displacement_m_sqrt_kg - old.displacement_m_sqrt_kg)?;
        modal_energy_before_j =
            finite(modal_energy_before_j + 0.5 * old.velocity_m_sqrt_kg_per_s.powi(2))?;
        modal_energy_j = finite(modal_energy_j + 0.5 * next.velocity_m_sqrt_kg_per_s.powi(2))?;
        candidate_states[k] = next;
    }
    for i in 0..network.models.len() {
        poll(Some(gate))?;
        network.candidates[i]
            .restore_states(&candidate_states[network.offsets[i]..network.offsets[i + 1]])
            .map_err(|source| ModalCouplingError::Component {
                component: i,
                source,
            })?;
    }
    // v_rel(t)=v0+a*t-(B·B)*I(t), and integral F I dt=J²/2.
    // This is an independent work integral, not an energy-residual correction.
    let reaction_work_j = finite(
        -v0 * impulse_n_s - acceleration * time_moment_n_s2
            + 0.5 * dot(column, column)? * impulse_n_s * impulse_n_s,
    )?;
    Ok(EventNativeMotion {
        modal_energy_before_j,
        modal_energy_j,
        tangent_increment_m: dot(column, &network.free_delta)?,
        tangent_velocity_m_s: projected_velocity(&network.candidates, column)?,
        external_work_j: dot(external, &network.free_delta)?,
        reaction_work_j,
    })
}

pub(crate) struct EventNativeAccount {
    pub sample: u64,
    pub energy_before_j: f64,
    pub modal_energy_j: f64,
    pub connection_energy_j: f64,
    pub connection_dissipation_j: f64,
    pub work_j: f64,
    pub residual_j: f64,
    pub tolerance_j: f64,
}

/// Stage the original native diagnostics and work ports without publication.
/// Admitted callers have free modes, zero acoustic transfers and no memories.
pub(crate) fn stage_account(
    network: &mut CoupledModalSystem,
    account: EventNativeAccount,
    connection_transfers: impl IntoIterator<Item = PortTransfer>,
) {
    let frame = &mut network.staged_frame;
    frame.sample = account.sample;
    frame.observer_pressure_pa = 0.0;
    frame.modal_energy_j = account.modal_energy_j;
    frame.connection_energy_j = account.connection_energy_j;
    frame.material_energy_j = 0.0;
    frame.radiation_energy_j = 0.0;
    frame.external_work_j = account.work_j;
    frame.component_dissipation_j = 0.0;
    frame.connection_dissipation_j = account.connection_dissipation_j;
    frame.material_dissipation_j = 0.0;
    frame.radiation_dissipation_j = 0.0;
    frame.energy_residual_j = account.residual_j;
    frame.energy_tolerance_j = account.tolerance_j;
    frame.solve_relative_residual = 0.0;
    frame
        .connection_forces_n
        .copy_from_slice(&network.reactions);
    network.staged_ports.accept(
        StepAccount {
            tick: account.sample,
            duration_s: network.dt,
            storage_change_j: Some(
                account.modal_energy_j + account.connection_energy_j - account.energy_before_j,
            ),
            loss_j: Some(account.connection_dissipation_j),
            included_work_j: account.work_j,
            entropy_change_j_k: 0.0,
            tolerance_j: Some(account.tolerance_j),
            unavailable: None,
        },
        connection_transfers,
    );
}
