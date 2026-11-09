//! Adaptive finite-wall heating, molecular supply/discharge and actual recoil.
//!
//! Prescribed sources remain external reservoirs. Every trial samples the
//! original thermal/particle owner at its actual composition, wall state and
//! velocity; inertial ambient velocity never follows the accelerating vessel.
//! Accepted-only receipts and fixed whole-window budgets share the vacuum
//! controller. This estimates local splitting error, not constitutive accuracy.

use super::*;
use crate::rocket_blowdown::effusion::{
    exchange::{EffusionAmbientPort, EffusionExchangeStream, EffusionReservoirPort},
    thermal::exchange::EffusionWallExchangeStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_TRANSFERS;
use fs_thermochem::{SpeciesId, mixture::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1};
use std::collections::BTreeMap;

/// One complete physical window. External sources accept the opposite of the
/// retained original port receipts; their finite storage is not advanced here.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallExchangeTrajectoryStep<S = ()> {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_intervals: usize,
    pub accepted_steps: Vec<EffusionWallExchangeStep<S>>,
    pub maximum_accepted_error: f64,
    pub incoming: EffusionExchangeStream,
    pub outgoing: EffusionExchangeStream,
    /// Independent gross energies in the same inertial frame with each original
    /// species' fixed internal-energy baseline removed BEFORE multiplication.
    pub incoming_reference_energy_world_j: f64,
    pub outgoing_reference_energy_world_j: f64,
    pub net_reference_energy_into_vessel_world_j: f64,
    /// Signed native accounts, never reconstructed by subtracting gross sums.
    pub net_mass_into_gas_kg: f64,
    pub net_relative_energy_into_gas_j: f64,
    pub net_momentum_into_vessel_world_kg_m_s: Vec3,
    pub net_energy_into_vessel_world_j: f64,
    pub external_impulse_n_s: Vec3,
    pub external_work_j: f64,
    /// Prescribed wall/gas heat plus radiation; molecular energy is separate.
    pub external_thermal_heat_j: f64,
    pub wall_storage_change_j: f64,
    pub prescribed_heat_residual_j: [f64; 2],
    pub residual: EffusionWallWindowBalance,
    pub debit: EffusionWallWindowBalance,
    pub maximum_species_mass_residual_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallExchangeTrajectoryCandidate<S = ()> {
    endpoint: GasRocketState,
    wall: EquilibriumPhaseState,
    step: EffusionWallExchangeTrajectoryStep<S>,
}
impl<S> EffusionWallExchangeTrajectoryCandidate<S> {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn wall(&self) -> &EquilibriumPhaseState {
        &self.wall
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionWallExchangeTrajectoryStep<S> {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasRocketState,
        EquilibriumPhaseState,
        EffusionWallExchangeTrajectoryStep<S>,
    ) {
        (self.endpoint, self.wall, self.step)
    }
}

impl<S> HistoryStep for EffusionWallExchangeStep<S> {
    fn exchange(&self) -> &SpeciesNetworkStep {
        &self.exchange.gas.exchange
    }
    fn momenta(&self) -> [Vec3; 3] {
        [
            self.exchange.incoming.momentum_world_kg_m_s,
            self.exchange.outgoing.momentum_world_kg_m_s,
            self.exchange.net_momentum_into_vessel_world_kg_m_s,
        ]
    }
}

#[derive(Default)]
struct StreamTotals {
    mass: Sum,
    relative_momentum: [Sum; 3],
    relative_energy: Sum,
    momentum: [Sum; 3],
    energy: Sum,
    reference_energy: Sum,
}
impl StreamTotals {
    fn add(
        &mut self,
        stream: EffusionExchangeStream,
        reference: f64,
    ) -> Result<(), EffusionRocketError> {
        self.mass.add(stream.mass_kg)?;
        add_vector(&mut self.relative_momentum, stream.relative_momentum_kg_m_s)?;
        self.relative_energy.add(stream.relative_energy_j)?;
        add_vector(&mut self.momentum, stream.momentum_world_kg_m_s)?;
        self.energy.add(stream.energy_world_j)?;
        self.reference_energy.add(reference)
    }
    fn stream(&self) -> EffusionExchangeStream {
        EffusionExchangeStream {
            mass_kg: self.mass.total(),
            relative_momentum_kg_m_s: vector(&self.relative_momentum),
            relative_energy_j: self.relative_energy.total(),
            momentum_world_kg_m_s: vector(&self.momentum),
            energy_world_j: self.energy.total(),
        }
    }
}

// Reference energies use original caloric integrals and the original incident
// beam correction. Reconstructing a tiny incoming tail from net+outgoing, or
// subtracting a formation datum from a gross reported energy, loses information.
fn reference_energies<S, E>(
    cx: &Cx<'_>,
    step: &EffusionWallExchangeStep<S>,
    species: &SpeciesAxis<'_>,
    references: &BTreeMap<SpeciesId, f64>,
    sources: &BTreeMap<u64, &ClosedGasParcel>,
) -> Result<[f64; 3], EffusionWallTrajectoryError<E>> {
    let mut energy = [Sum::default(), Sum::default(), Sum::default()];
    let exchange = &step.exchange;
    let gas_t = step.thermal.report.next_gas.state().state().temperature;
    for row in &exchange.gas.exchange.boundaries {
        check(cx)?;
        let transfer = &row.transfer;
        let index = species
            .sources
            .binary_search_by(|s| s.species().cmp(&transfer.species))
            .map_err(|_| {
                EffusionTrajectoryError::Invalid("molecular receipt lost its original species")
            })?;
        let source = species.sources[index];
        let reference =
            references
                .get(&transfer.species)
                .ok_or(EffusionTrajectoryError::Invalid(
                    "molecular receipt lost its fixed caloric baseline",
                ))?;
        let du = |t| {
            source
                .molar_internal_energy_change(Temperature::new(*reference), Temperature::new(t))
                .map(|value| value.value() / source.molar_mass().value())
                .map_err(EffusionWallTrajectoryError::Caloric)
        };
        let gas_u = du(gas_t)?;
        let gas_rt = R * gas_t / source.molar_mass().value();
        energy[1].add(transfer.outgoing_mass_kg * (gas_u + 0.5 * gas_rt))?;
        energy[2].add(transfer.net_mass_into_gas_kg * (gas_u + gas_rt))?;
        energy[2].add(transfer.relative_energy_into_gas_j)?;
        if transfer.incoming_mass_kg > 0.0 {
            let incoming_t = sources
                .get(&row.id)
                .ok_or(EffusionTrajectoryError::Invalid(
                    "incoming receipt lost its external source",
                ))?
                .state()
                .state()
                .temperature;
            let rt = R * incoming_t / source.molar_mass().value();
            let correction = match exchange
                .gas
                .drifting_reservoirs
                .binary_search_by_key(&row.id, |p| p.id)
            {
                Ok(i) => {
                    let beams = &exchange.gas.drifting_reservoirs[i].incoming.species;
                    let j = beams
                        .binary_search_by(|b| b.species.cmp(&transfer.species))
                        .map_err(|_| {
                            EffusionTrajectoryError::Invalid(
                                "incoming receipt lost its original drift beam",
                            )
                        })?;
                    beams[j].incoming.additional_specific_energy_j_kg.ok_or(
                        EffusionTrajectoryError::Invalid(
                            "active drift beam lost its crossing-energy correction",
                        ),
                    )?
                }
                Err(_) => -0.5 * rt,
            };
            energy[0].add(transfer.incoming_mass_kg * (du(incoming_t)? + rt + correction))?;
        }
    }
    let mechanics = &exchange.mechanics;
    let mean = mechanics.mean_velocity_world_m_s;
    let mean_squared = mechanics.mean_speed_squared_m2_s2;
    for (total, (stream, native)) in energy[..2].iter_mut().zip([
        (exchange.incoming, mechanics.incoming),
        (exchange.outgoing, mechanics.outgoing),
    ]) {
        total.add(native.kinetic_transport_energy_j)?;
        total.add(
            stream
                .relative_momentum_kg_m_s
                .sub(native.relative_momentum_n_s)
                .dot(mean),
        )?;
        total.add(0.5 * (stream.mass_kg - native.mass_kg) * mean_squared)?;
    }
    energy[2].add(
        exchange
            .incoming
            .relative_momentum_kg_m_s
            .sub(exchange.outgoing.relative_momentum_kg_m_s)
            .dot(mean),
    )?;
    energy[2].add(0.5 * exchange.net_mass_into_gas_kg * mean_squared)?;
    Ok(energy.map(|s| s.total()))
}

impl GasRocketState {
    /// Advance one finite-wall window with vacuum, comoving and inertial ports.
    /// Whole-window heat totals are held powers. Actual wall/gas heat is internal;
    /// incoming molecular energy, including drift, enters exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange_with_wall_adaptive(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        ambient_ports: &[EffusionAmbientPort<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        solve: EffusionWallSolve,
        adaptive: EffusionWallAdaptiveSolve,
    ) -> Result<EffusionWallExchangeTrajectoryCandidate, EffusionWallTrajectoryError> {
        self.advance_effusion_exchange_with_wall_convection_adaptive(
            cx,
            environment,
            vacuum_outlets,
            reservoir_ports,
            ambient_ports,
            wall,
            body,
            gas_heat_j,
            solve,
            adaptive,
            |_| {
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: body.convection_w_per_m2_k,
                    source: (),
                })
            },
        )
    }

    /// Use the shared bounded full/two-half controller with a repeatable original
    /// convection source. Only typed gross withdrawal can cause physical retry;
    /// source limits, balance failures, cancellation and work exhaustion refuse
    /// atomically. Callback side effects from rejected trials are not undone.
    ///
    /// Every active external species joins the error axis, initially at zero
    /// mass when absent. Its baseline is the original vessel temperature if
    /// present, otherwise its first source in canonical port-ID order. All
    /// accepted thermal, incoming/outgoing and signed net receipts use the same
    /// species baselines and declared inertial frame. Export tolerances are
    /// frame-specific; changing observer does not promise identical step sizes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange_with_wall_convection_adaptive<S, E>(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        ambient_ports: &[EffusionAmbientPort<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        solve: EffusionWallSolve,
        adaptive: EffusionWallAdaptiveSolve,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<EffusionWallExchangeTrajectoryCandidate<S>, EffusionWallTrajectoryError<E>> {
        check(cx)?;
        if vacuum_outlets
            .len()
            .checked_add(reservoir_ports.len())
            .and_then(|n| n.checked_add(ambient_ports.len()))
            .is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(
                EffusionTrajectoryError::Invalid("adaptive molecular port work bound").into(),
            );
        }
        let mut sources = BTreeMap::new();
        for (id, source) in reservoir_ports
            .iter()
            .filter(|p| p.area_m2 > 0.0)
            .map(|p| (p.id, p.source))
            .chain(
                ambient_ports
                    .iter()
                    .filter(|p| p.area_m2 > 0.0)
                    .map(|p| (p.id, p.source)),
            )
        {
            check(cx)?;
            if sources.insert(id, source).is_some() {
                return Err(EffusionTrajectoryError::Invalid(
                    "external molecular source IDs must be unique",
                )
                .into());
            }
        }
        let mut references = BTreeMap::new();
        for gas in std::iter::once(self.reservoir()).chain(sources.values().copied()) {
            check(cx)?;
            for source in gas.model().thermodynamics().components() {
                check(cx)?;
                references
                    .entry(source.species().clone())
                    .or_insert(gas.state().state().temperature);
                if references.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
                    return Err(EffusionTrajectoryError::Invalid(
                        "adaptive species union exceeds the original gas work bound",
                    )
                    .into());
                }
            }
        }
        let (mut incoming, mut outgoing) = (StreamTotals::default(), StreamTotals::default());
        let mut net_relative_energy = Sum::default();
        let history = integrate(
            self,
            sources.values().copied(),
            cx,
            wall,
            gas_heat_j,
            adaptive,
            |state, input, heat| {
                state
                    .advance_effusion_exchange_with_wall_convection(
                        cx,
                        environment,
                        vacuum_outlets,
                        reservoir_ports,
                        ambient_ports,
                        input,
                        body,
                        heat,
                        solve,
                        &mut resolve,
                    )
                    .map(|candidate| candidate.into_parts())
            },
            |step, species| reference_energies(cx, step, species, &references, &sources),
            |totals, step, energies, heat| {
                let exchange = &step.exchange;
                incoming.add(exchange.incoming, energies[0])?;
                outgoing.add(exchange.outgoing, energies[1])?;
                net_relative_energy.add(exchange.net_relative_energy_into_gas_j)?;
                // The shared audit stores net OUTWARD mass/energy/momentum.
                totals.mass.add(-exchange.net_mass_into_gas_kg)?;
                add_vector(
                    &mut totals.momentum,
                    exchange.net_momentum_into_vessel_world_kg_m_s.scale(-1.0),
                )?;
                totals
                    .energy
                    .add(-exchange.net_energy_into_vessel_world_j)?;
                totals.reference_energy.add(-energies[2])?;
                add_vector(&mut totals.impulse, exchange.mechanics.external_impulse_n_s)?;
                totals.work.add(exchange.mechanics.external_work_j)?;
                totals.radiation.add(step.thermal.report.radiation_heat_j)?;
                for (sum, value) in totals.prescribed.iter_mut().zip(heat) {
                    sum.add(value)?;
                }
                totals.add_balance(BalanceReceipt {
                    exchange: &exchange.gas.exchange,
                    initial_mass_residual_kg: exchange.initial_mass_residual_kg,
                    momentum_residual: exchange.momentum_residual_kg_m_s,
                    initial_momentum_correction: exchange.initial_inventory_momentum_correction_n_s,
                    gravity_impulse_uncertainty: exchange.gravity_impulse_uncertainty_n_s,
                    energy_residual_j: step.combined_energy_residual_j,
                    absolute_energy_residual_sum_j: step.absolute_energy_residual_sum_j,
                    energy_uncertainty_bound_j: step.energy_uncertainty_bound_j,
                })
            },
        )?;
        let totals = history.totals;
        Ok(EffusionWallExchangeTrajectoryCandidate {
            endpoint: history.endpoint,
            wall: history.wall,
            step: EffusionWallExchangeTrajectoryStep {
                duration_s: history.duration_s,
                attempted_intervals: history.attempted_intervals,
                rejected_intervals: history.rejected_intervals,
                accepted_intervals: history.accepted_intervals,
                accepted_steps: history.accepted_steps,
                maximum_accepted_error: history.maximum_accepted_error,
                incoming: incoming.stream(),
                outgoing: outgoing.stream(),
                incoming_reference_energy_world_j: incoming.reference_energy.total(),
                outgoing_reference_energy_world_j: outgoing.reference_energy.total(),
                net_reference_energy_into_vessel_world_j: -totals.reference_energy.total(),
                net_mass_into_gas_kg: -totals.mass.total(),
                net_relative_energy_into_gas_j: net_relative_energy.total(),
                net_momentum_into_vessel_world_kg_m_s: vector(&totals.momentum).scale(-1.0),
                net_energy_into_vessel_world_j: -totals.energy.total(),
                external_impulse_n_s: vector(&totals.impulse),
                external_work_j: totals.work.total(),
                external_thermal_heat_j: history.external_thermal_heat_j,
                wall_storage_change_j: history.wall_storage_change_j,
                prescribed_heat_residual_j: history.prescribed_heat_residual_j,
                residual: history.residual,
                debit: history.debit,
                maximum_species_mass_residual_kg: history.maximum_species_mass_residual_kg,
            },
        })
    }
}
