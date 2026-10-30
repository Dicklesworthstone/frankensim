//! One original-body transaction for single and simultaneous junction events.
//! The event propagator supplies a trajectory and its error estimate. This
//! module independently restores original component states and computes their
//! physical pressure, spring storage, reactions and work ports. It has no
//! contact law, phase selection or authority to publish accepted mechanics.

use super::super::super::super::{component_energy, extension};
use super::super::super::{
    CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, poll,
};
use super::super::linear_network::LinearNetworkLink;
use super::{EventNativeAccount, projected_velocity, stage_account_with_numerical_loss};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::port_adoption::PortTransfer;
use fs_exec::CancelGate;
use fs_math::det;

pub(crate) enum NativeLinearEventError {
    Network(ModalCouplingError),
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
}
impl From<ModalCouplingError> for NativeLinearEventError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Network(error)
    }
}

/// Borrowed event result; all these quantities are checked against the
/// supplied native bodies rather than installed as native stored energy.
pub(crate) struct NativeLinearEventTrajectory<'a> {
    pub positions: &'a [f64],
    pub velocities: &'a [f64],
    pub connection_dissipation_j: &'a [f64],
    pub initial_energy_j: f64,
    pub modal_energy_j: f64,
    pub connection_energy_j: f64,
    pub external_work_j: f64,
    pub modal_damping_dissipation_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub numerical_work_uncertainty_j: f64,
}

#[derive(Clone, Copy)]
pub(crate) struct NativeLinearEventFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub external_work_j: f64,
    pub connection_damping_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Bounds sqrt(2 E(t)) over the entire passive trajectory under the held
    /// original external force. Contact owners use it for interval speed gates.
    pub root_energy_bound: f64,
}

pub(crate) fn require_linear_network(
    network: &CoupledModalSystem,
    maximum_modes: usize,
    maximum_connections: usize,
) -> Result<(), ModalCouplingError> {
    if network.samples_rendered() != 0
        || network.mode_count() > maximum_modes
        || network.connections.len() > maximum_connections
        || !network.relaxation_branches.is_empty()
        || !network.plane_stress_material_points.is_empty()
        || network.radiation.is_some()
        || network
            .connections
            .iter()
            .any(|link| link.constant_loss_reference_omega_rad_s.is_some())
    {
        return Err(invalid(
            "adhesive network events require a complete sample-zero linear modal network with Kelvin-Voigt connections and no constitutive/acoustic memory",
        ));
    }
    Ok(())
}

pub(crate) fn network_parameters(
    network: &CoupledModalSystem,
) -> Result<(Vec<f64>, Vec<f64>, Vec<LinearNetworkLink>), ModalCouplingError> {
    let stiffness = network
        .models
        .iter()
        .flat_map(|model| model.modes())
        .map(|mode| finite(mode.angular_frequency_rad_s * mode.angular_frequency_rad_s))
        .collect::<Result<Vec<_>, _>>()?;
    let damping = network
        .models
        .iter()
        .flat_map(|model| model.modes())
        .map(|mode| finite(2.0 * mode.damping_ratio * mode.angular_frequency_rad_s))
        .collect::<Result<Vec<_>, _>>()?;
    let links = network
        .columns
        .iter()
        .zip(&network.connections)
        .map(|(column, link)| LinearNetworkLink {
            column: column.clone(),
            stiffness_n_per_m: link.stiffness_n_m,
            damping_n_s_per_m: link.damping_n_s_m,
            rest_extension_m: link.rest_extension_m,
        })
        .collect();
    Ok((stiffness, damping, links))
}

pub(crate) struct NativeLinearNetworkWorkspace {
    positions: Vec<f64>,
    velocities: Vec<f64>,
    candidate_states: Vec<ModalAcousticState>,
    connection_work: Vec<PortTransfer>,
}
impl NativeLinearNetworkWorkspace {
    pub fn new(network: &CoupledModalSystem) -> Self {
        Self {
            positions: vec![0.0; network.mode_count()],
            velocities: vec![0.0; network.mode_count()],
            candidate_states: network
                .models
                .iter()
                .flat_map(|model| model.states())
                .copied()
                .collect(),
            connection_work: vec![PortTransfer::work(0.0, false); network.connections.len()],
        }
    }
    pub fn capture(&mut self, network: &CoupledModalSystem) {
        for (k, state) in network
            .models
            .iter()
            .flat_map(|model| model.states())
            .enumerate()
        {
            self.positions[k] = state.displacement_m_sqrt_kg;
            self.velocities[k] = state.velocity_m_sqrt_kg_per_s;
        }
    }
    pub fn positions(&self) -> &[f64] {
        &self.positions
    }
    pub fn velocities(&self) -> &[f64] {
        &self.velocities
    }

    /// The caller retains accepted junction phases/history and decides whether
    /// to publish after its additional source, speed and slip checks succeed.
    pub fn stage(
        &mut self,
        network: &mut CoupledModalSystem,
        external: &[f64],
        surface_energy_j: f64,
        path: NativeLinearEventTrajectory<'_>,
        gate: &CancelGate,
    ) -> Result<NativeLinearEventFrame, NativeLinearEventError> {
        self.stage_with_numerical_loss(network, external, surface_energy_j, path, 0.0, None, gate)
    }

    /// Original mechanics and work ports for a dissipative time discretization.
    /// The numerical loss is independently derived by that integrator, never a
    /// residual correction, and is retained separately from physical losses.
    #[allow(clippy::too_many_arguments)]
    pub fn stage_with_numerical_loss(
        &mut self,
        network: &mut CoupledModalSystem,
        external: &[f64],
        surface_energy_j: f64,
        path: NativeLinearEventTrajectory<'_>,
        numerical_loss_j: f64,
        connection_numerical_loss_j: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<NativeLinearEventFrame, NativeLinearEventError> {
        poll(Some(gate))?;
        if !numerical_loss_j.is_finite() || numerical_loss_j < 0.0 {
            return Err(invalid("nonnegative finite numerical dissipation required").into());
        }
        if let Some(losses) = connection_numerical_loss_j {
            if losses.len() != network.connections.len()
                || losses.iter().any(|loss| !loss.is_finite() || *loss < 0.0)
            {
                return Err(invalid("each original connection requires a nonnegative discrete spring loss").into());
            }
            let sum = losses.iter().try_fold(0.0, |sum, loss| finite(sum + loss))?;
            if sum > numerical_loss_j * (1.0 + 64.0 * f64::EPSILON) {
                return Err(invalid("connection numerical loss exceeds the original total").into());
            }
        }
        let n = network.mode_count();
        if external.len() != n
            || external.iter().any(|f| !f.is_finite())
            || path.positions.len() != n
            || path.velocities.len() != n
            || path.connection_dissipation_j.len() != network.connections.len()
        {
            return Err(invalid(
                "event trajectory and force must cover every original mode and connection",
            )
            .into());
        }
        let sample = network
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| invalid("sample clock overflow"))?;
        let before = network.total_energy_j()?;
        let root_energy =
            finite(det::sqrt(2.0 * before) + network.dt * det::sqrt(dot(external, external)?))?;
        for k in 0..n {
            self.candidate_states[k] = ModalAcousticState {
                displacement_m_sqrt_kg: path.positions[k],
                velocity_m_sqrt_kg_per_s: path.velocities[k],
            };
            network.old_q[k] = self.positions[k];
            network.free_delta[k] = finite(path.positions[k] - self.positions[k])?;
        }
        let mut pressure = 0.0;
        for i in 0..network.models.len() {
            poll(Some(gate))?;
            let candidate = &mut network.candidates[i];
            candidate
                .restore_states(&self.candidate_states[network.offsets[i]..network.offsets[i + 1]])
                .map_err(|source| ModalCouplingError::Component {
                    component: i,
                    source,
                })?;
            let mut component_pressure = 0.0;
            for (mode, state) in candidate.modes().iter().zip(candidate.states()) {
                component_pressure = finite(
                    component_pressure
                        + mode.pressure_per_modal_velocity.re * state.velocity_m_sqrt_kg_per_s
                        + mode.pressure_per_modal_velocity.im
                            * mode.angular_frequency_rad_s
                            * state.displacement_m_sqrt_kg,
                )?;
            }
            limit(
                "original component observer pressure",
                component_pressure.abs(),
                candidate.budget().maximum_abs_pressure_pa,
            )?;
            pressure = finite(pressure + component_pressure)?;
        }
        limit(
            "original combined observer pressure",
            pressure.abs(),
            network.config.maximum_abs_pressure_pa,
        )?;
        let modal_energy = component_energy(&network.candidates)?;
        let mut spring_energy = 0.0;
        let mut connection_loss = 0.0;
        for (i, (column, link)) in network.columns.iter().zip(&network.connections).enumerate() {
            poll(Some(gate))?;
            let z0 = extension(&network.models, column, link.rest_extension_m)?;
            let z1 = extension(&network.candidates, column, link.rest_extension_m)?;
            let speed = projected_velocity(&network.candidates, column)?;
            let loss = finite(path.connection_dissipation_j[i])?;
            if loss < 0.0 {
                return Err(invalid("negative original connection loss").into());
            }
            let change = finite(0.5 * link.stiffness_n_m * (z1 - z0) * (z1 + z0))?;
            let numerical_spring_loss = connection_numerical_loss_j.map_or(0.0, |losses| losses[i]);
            self.connection_work[i] =
                PortTransfer::work(finite(-change - loss - numerical_spring_loss)?, false);
            connection_loss = finite(connection_loss + loss)?;
            spring_energy = finite(spring_energy + 0.5 * link.stiffness_n_m * z1 * z1)?;
            network.reactions[i] = finite(-link.stiffness_n_m * z1 - link.damping_n_s_m * speed)?;
            let bound = finite(
                (det::sqrt(link.stiffness_n_m)
                    + link.damping_n_s_m * det::sqrt(dot(column, column)?))
                    * root_energy,
            )?;
            limit(
                "original connection force over complete event trajectory",
                bound,
                network.config.maximum_abs_connection_force_n,
            )?;
        }
        let external_work = dot(external, &network.free_delta)?;
        let after = finite(modal_energy + spring_energy)?;
        let loss = finite(
            path.modal_damping_dissipation_j + connection_loss + path.friction_dissipation_j,
        )?;
        let uncertainty = finite(path.numerical_work_uncertainty_j)?;
        if path.modal_damping_dissipation_j < 0.0
            || path.friction_dissipation_j < 0.0
            || uncertainty < 0.0
        {
            return Err(invalid(
                "nonnegative modal/junction loss and numerical uncertainty required",
            )
            .into());
        }
        limit(
            "original network event energy",
            after,
            network.config.maximum_total_energy_j,
        )?;
        limit(
            "combined adhesive event energy magnitude",
            finite(after + surface_energy_j)?.abs(),
            network.config.maximum_total_energy_j,
        )?;
        let residual = finite(after - before + loss + numerical_loss_j - external_work)?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(external_work.abs())
            .max(loss)
            .max(numerical_loss_j)
            .max(surface_energy_j.abs());
        let tolerance = finite(
            network.config.energy_absolute_tolerance_j
                + network.config.energy_relative_tolerance * scale,
        )?;
        for residual_j in [
            residual,
            before - path.initial_energy_j,
            modal_energy - path.modal_energy_j,
            spring_energy - path.connection_energy_j,
            external_work - path.external_work_j,
            connection_loss - path.connection_damping_dissipation_j,
        ] {
            if finite(residual_j.abs() + uncertainty)? > tolerance {
                return Err(NativeLinearEventError::EnergyAccuracy {
                    residual_j,
                    estimated_error_j: uncertainty,
                    tolerance_j: tolerance,
                });
            }
        }
        stage_account_with_numerical_loss(
            network,
            EventNativeAccount {
                sample,
                observer_pressure_pa: pressure,
                energy_before_j: before,
                modal_energy_j: modal_energy,
                component_dissipation_j: path.modal_damping_dissipation_j,
                connection_energy_j: spring_energy,
                connection_dissipation_j: connection_loss,
                work_j: finite(external_work - path.friction_dissipation_j)?,
                residual_j: residual,
                tolerance_j: tolerance,
            },
            numerical_loss_j,
            self.connection_work.iter().copied(),
        );
        Ok(NativeLinearEventFrame {
            sample,
            observer_pressure_pa: pressure,
            network_energy_j: after,
            external_work_j: external_work,
            connection_damping_dissipation_j: connection_loss,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            root_energy_bound: root_energy,
        })
    }
}
