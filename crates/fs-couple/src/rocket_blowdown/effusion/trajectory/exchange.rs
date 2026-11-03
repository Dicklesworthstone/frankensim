//! Adaptive finite-vessel supply and discharge through the shared controller.
//!
//! Prescribed sources remain instantaneously comoving and thermochemically held
//! throughout the requested interval. Each native half step resamples the actual
//! finite gas. Its original source admission, gross withdrawal and complete
//! momentum/energy budgets remain binding; rejected trials contribute no receipts.

use super::super::{
    EffusionOutlet, EffusionRocketEnvironment, EffusionRocketSolve, GasRocketState,
    exchange::{EffusionExchangeStep, EffusionExchangeStream, EffusionReservoirPort},
};
use super::{EffusionAdaptiveSolve, EffusionTrajectoryError, integrate, poll};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_TRANSFERS;
use fs_mbd::Vec3;

/// Accepted native history and compensated whole-interval physical receipts.
/// Temporal tolerances estimate local error; physical budgets apply to each
/// accepted half step. Neither constitutes a certified global accuracy bound.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionExchangeTrajectoryStep {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_intervals: usize,
    /// Two original half steps per accepted interval, including each port's receipts.
    pub accepted_steps: Vec<EffusionExchangeStep>,
    pub maximum_accepted_error: f64,
    pub total_incoming: EffusionExchangeStream,
    pub total_outgoing: EffusionExchangeStream,
    /// Authoritative signed source-coordinate totals, accumulated separately
    /// from gross streams so counterflow retains its small net mass and heat.
    pub total_net_mass_into_gas_kg: f64,
    pub total_net_relative_energy_into_gas_j: f64,
    pub total_net_momentum_into_vessel_world_kg_m_s: Vec3,
    pub total_net_energy_into_vessel_world_j: f64,
    pub total_external_impulse_n_s: Vec3,
    pub total_external_work_j: f64,
    pub total_mass_residual_kg: f64,
    pub total_momentum_residual_kg_m_s: Vec3,
    pub total_energy_residual_j: f64,
    pub total_gravity_impulse_uncertainty_n_s: Vec3,
    pub total_gravity_work_uncertainty_j: f64,
}

/// No partial endpoint or receipt is published when any trial ultimately refuses.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionExchangeTrajectoryCandidate {
    endpoint: GasRocketState,
    step: EffusionExchangeTrajectoryStep,
}
impl EffusionExchangeTrajectoryCandidate {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionExchangeTrajectoryStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(self) -> (GasRocketState, EffusionExchangeTrajectoryStep) {
        (self.endpoint, self.step)
    }
}

fn stream_values(stream: &EffusionExchangeStream) -> [f64; 9] {
    let relative = stream.relative_momentum_kg_m_s;
    let world = stream.momentum_world_kg_m_s;
    [
        stream.mass_kg,
        relative.x,
        relative.y,
        relative.z,
        stream.relative_energy_j,
        world.x,
        world.y,
        world.z,
        stream.energy_world_j,
    ]
}
fn values(step: &EffusionExchangeStep) -> [f64; 37] {
    let mut values = [0.0; 37];
    values[..9].copy_from_slice(&stream_values(&step.incoming));
    values[9..18].copy_from_slice(&stream_values(&step.outgoing));
    let net = step.net_momentum_into_vessel_world_kg_m_s;
    let external = step.mechanics.external_impulse_n_s;
    let residual = step.momentum_residual_kg_m_s;
    let reserve = step.gravity_impulse_uncertainty_n_s;
    values[18..].copy_from_slice(&[
        step.net_mass_into_gas_kg,
        step.net_relative_energy_into_gas_j,
        net.x,
        net.y,
        net.z,
        step.net_energy_into_vessel_world_j,
        external.x,
        external.y,
        external.z,
        step.mechanics.external_work_j,
        step.gas.exchange.total_mass_residual_kg,
        residual.x,
        residual.y,
        residual.z,
        step.energy_residual_j,
        reserve.x,
        reserve.y,
        reserve.z,
        step.gravity_work_uncertainty_j,
    ]);
    values
}
fn stream(values: &[f64]) -> EffusionExchangeStream {
    EffusionExchangeStream {
        mass_kg: values[0],
        relative_momentum_kg_m_s: Vec3::new(values[1], values[2], values[3]),
        relative_energy_j: values[4],
        momentum_world_kg_m_s: Vec3::new(values[5], values[6], values[7]),
        energy_world_j: values[8],
    }
}

impl GasRocketState {
    /// Resolve one complete interval of molecular filling and discharge.
    ///
    /// This uses the same step-doubling controller, exact time-end handling and
    /// finite work limits as vacuum recoil. The error axis is the union of the
    /// original vessel and active prescribed-source species: absent species
    /// have zero mass, and every endpoint retains its original caloric identity.
    /// Thus species entering a previously pure gas participate immediately in
    /// refinement. Geometrically closed ports do not inspect inactive sources.
    ///
    /// Only original typed gross-withdrawal failures can trigger physical-step
    /// retries. Invalid sources, physical balance failures and cancellation
    /// remain atomic refusals. Loads, source states and aperture normals stay
    /// fixed; this is not a stationary ambient or a rotating-vessel model.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange_adaptive(
        &self,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        solve: EffusionRocketSolve,
        adaptive: EffusionAdaptiveSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionExchangeTrajectoryCandidate, EffusionTrajectoryError> {
        poll(&mut cancelled)?;
        if vacuum_outlets
            .len()
            .checked_add(reservoir_ports.len())
            .is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(EffusionTrajectoryError::Invalid(
                "adaptive molecular port work bound",
            ));
        }
        let sources = reservoir_ports
            .iter()
            .filter(|p| p.area_m2 > 0.0)
            .map(|p| p.source);
        let window = integrate(
            self,
            sources,
            duration_s,
            adaptive,
            |state, dt, cancel| {
                state
                    .advance_effusion_exchange(
                        environment,
                        vacuum_outlets,
                        reservoir_ports,
                        dt,
                        endpoint,
                        solve,
                        cancel,
                    )
                    .map(|candidate| candidate.into_parts())
            },
            values,
            &mut cancelled,
        )?;
        let t = window.totals;
        Ok(EffusionExchangeTrajectoryCandidate {
            endpoint: window.endpoint,
            step: EffusionExchangeTrajectoryStep {
                duration_s: window.duration_s,
                attempted_intervals: window.attempted_intervals,
                rejected_intervals: window.rejected_intervals,
                accepted_intervals: window.accepted_intervals,
                accepted_steps: window.accepted_steps,
                maximum_accepted_error: window.maximum_accepted_error,
                total_incoming: stream(&t[..9]),
                total_outgoing: stream(&t[9..18]),
                total_net_mass_into_gas_kg: t[18],
                total_net_relative_energy_into_gas_j: t[19],
                total_net_momentum_into_vessel_world_kg_m_s: Vec3::new(t[20], t[21], t[22]),
                total_net_energy_into_vessel_world_j: t[23],
                total_external_impulse_n_s: Vec3::new(t[24], t[25], t[26]),
                total_external_work_j: t[27],
                total_mass_residual_kg: t[28],
                total_momentum_residual_kg_m_s: Vec3::new(t[29], t[30], t[31]),
                total_energy_residual_j: t[32],
                total_gravity_impulse_uncertainty_n_s: Vec3::new(t[33], t[34], t[35]),
                total_gravity_work_uncertainty_j: t[36],
            },
        })
    }
}
