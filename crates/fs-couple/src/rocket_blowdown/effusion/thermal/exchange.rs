//! Finite material heat storage with molecular supply, discharge and recoil.
//!
//! The SAME implicit wall/gas owner as vacuum blowdown heats the original gas,
//! then the original bidirectional particle owner samples that heated state.
//! Vacuum, comoving reservoirs and fixed-inertial ambient Maxwellians all share
//! one original species withdrawal budget and one mechanical update. Both
//! stages represent the same physical interval, a first-order Lie split.
//!
//! Material wall heat is an internal transfer. Incoming molecular energy and
//! tangential drift arrive through their actual particles; neither is added as
//! a second heat source. Prescribed source inventories are not depleted. The
//! existing fixed-volume, constant-density, nonrotating and dilute-aperture
//! assumptions, source qualification and thermal Biot gates remain binding.

use super::{
    ConvectionCoefficient, Cx, EffusionOutlet, EffusionRocketEnvironment, EffusionWallError,
    EffusionWallSolve, EquilibriumPhaseState, FiniteGasThermalReport, GasRocketState,
    LumpedThermalEnvironment, ResolvedThermochemicalGasState, UniformEnthalpyStep,
    UniformEnthalpyStepInput, admit_original_gas, admit_wall, energy_account, exchange_heat,
};
use crate::rocket_blowdown::effusion::exchange::{
    EffusionAmbientPort, EffusionExchangeStep, EffusionReservoirPort,
};

/// Original thermal and particle reports for one physical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallExchangeStep<S = ()> {
    pub duration_s: f64,
    /// Accepted finite wall and heated gas BEFORE molecular supply/discharge.
    pub thermal: UniformEnthalpyStep<FiniteGasThermalReport<S>>,
    /// Original source streams, actual finite inventory and mechanical update.
    pub exchange: EffusionExchangeStep,
    /// Body/gas prescribed heat plus enclosing radiation. Molecular supply is
    /// kept in `exchange`, and reciprocal wall/gas heat cancels internally.
    pub external_thermal_heat_j: f64,
    /// Original source-relative thermal plus gas/mechanical balance defects.
    /// No gross formation energies are subtracted to form this residual.
    pub combined_energy_residual_j: f64,
    pub absolute_energy_residual_sum_j: f64,
    /// Initial native inventory remapping and gravity mass-history reserves.
    pub energy_uncertainty_bound_j: f64,
}

/// Wall, finite gas and moving vessel publish together after all original gates.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallExchangeCandidate<S = ()> {
    endpoint: GasRocketState,
    step: EffusionWallExchangeStep<S>,
}
impl<S> EffusionWallExchangeCandidate<S> {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn wall(&self) -> &EquilibriumPhaseState {
        &self.step.thermal.state
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionWallExchangeStep<S> {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasRocketState,
        EquilibriumPhaseState,
        EffusionWallExchangeStep<S>,
    ) {
        (self.endpoint, self.step.thermal.state, self.step)
    }
}

impl GasRocketState {
    /// Exchange finite wall heat and molecular particles with all original ports.
    ///
    /// Reservoir sources comove with the aperture; ambient sources retain their
    /// declared inertial velocity. Both thermochemical source states are held
    /// during this interval and reported as external particle exchanges. Normal
    /// stresses and tangential momentum use the original crossing distributions.
    /// The wall is already included in vessel dry mass; its heat storage is
    /// counted once without adding material mass or kinetic energy again.
    ///
    /// `body` supplies the independently qualified whole-surface heat coefficient.
    /// Prescribed body/gas heat excludes both reciprocal wall heat and molecular
    /// energy. Time refinement must resolve the first-order heat/exchange split.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange_with_wall(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        ambient_ports: &[EffusionAmbientPort<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        controls: EffusionWallSolve,
    ) -> Result<EffusionWallExchangeCandidate, EffusionWallError> {
        self.advance_effusion_exchange_with_wall_convection(
            cx,
            environment,
            vacuum_outlets,
            reservoir_ports,
            ambient_ports,
            wall,
            body,
            gas_heat_j,
            controls,
            |_| {
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: body.convection_w_per_m2_k,
                    source: (),
                })
            },
        )
    }

    /// Resolve the original convection callback on each thermal trial, then
    /// sample ALL particle paths from the same accepted heated gas. Subsequent
    /// calls feed actual new composition, gas temperature and wall storage back
    /// to that owner; actual vessel velocity selects ambient incidence.
    ///
    /// Every active aperture admits the unheated gas first, so heating cannot
    /// hide an originally unsupported source. The original particle owner then
    /// admits external sources, full heated/final states, IDs, normals, frames,
    /// drift and pooled gross withdrawals. Zero-area sources remain inactive.
    /// Any late refusal leaves wall, gas, vehicle and prescribed sources intact.
    /// Callback side effects are not rolled back and callbacks must be repeatable.
    /// No additional exterior wall-pressure load or photon impulse is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange_with_wall_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        ambient_ports: &[EffusionAmbientPort<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        controls: EffusionWallSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<EffusionWallExchangeCandidate<S>, EffusionWallError<E>> {
        cx.checkpoint().map_err(|_| EffusionWallError::Cancelled)?;
        let count = vacuum_outlets
            .len()
            .checked_add(reservoir_ports.len())
            .and_then(|count| count.checked_add(ambient_ports.len()))
            .ok_or(EffusionWallError::Invalid(
                "molecular port work bound overflow",
            ))?;
        admit_wall(self, cx, environment, wall, controls, count)?;
        admit_original_gas(
            self,
            cx,
            vacuum_outlets
                .iter()
                .map(|p| (p.id, p.area_m2, p.interface.interface()))
                .chain(
                    reservoir_ports
                        .iter()
                        .map(|p| (p.id, p.area_m2, p.interface.interface())),
                )
                .chain(ambient_ports.iter().map(|p| {
                    (
                        p.id,
                        p.area_m2,
                        p.interface.momentum_interface().interface(),
                    )
                })),
        )?;
        let (heated, thermal) = exchange_heat(self, cx, wall, body, gas_heat_j, controls, resolve)?;
        let (endpoint, exchange) = heated
            .advance_effusion_exchange_with_ambient(
                environment,
                vacuum_outlets,
                reservoir_ports,
                ambient_ports,
                wall.duration_s,
                controls.gas,
                controls.recoil,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let account = energy_account(
            &thermal,
            wall,
            gas_heat_j,
            exchange.energy_residual_j,
            exchange.initial_inventory_kinetic_correction_j.abs()
                + exchange.gravity_work_uncertainty_j,
            controls.combined_energy_tolerance_j,
        )?;
        cx.checkpoint().map_err(|_| EffusionWallError::Cancelled)?;
        Ok(EffusionWallExchangeCandidate {
            endpoint,
            step: EffusionWallExchangeStep {
                duration_s: wall.duration_s,
                thermal,
                exchange,
                external_thermal_heat_j: account.external_thermal_heat_j,
                combined_energy_residual_j: account.combined_energy_residual_j,
                absolute_energy_residual_sum_j: account.absolute_energy_residual_sum_j,
                energy_uncertainty_bound_j: account.energy_uncertainty_bound_j,
            },
        })
    }
}
