//! Reaction products, reciprocal material-wall heat and motion on one interval.
//!
//! The original reaction owner changes species at fixed total internal energy.
//! The original thermal owner then solves gas pressure/work, wall momentum and
//! finite-body heat exchange together. Both stages approximate the same dt;
//! chemistry is never repeated by a heat or momentum residual iteration.

/// Bounded refinement of chemistry, wall heat and piston motion together.
pub mod subcycling;
/// Same-endpoint chemistry, motion and reciprocal finite-wall heat exchange.
pub mod coupled;

use super::{
    GasPiston, GasReaction, GasReactionError, GasReactionSolve, GasReactionStep, PistonInput,
    PistonSolve,
};
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError, ThermalPistonStep,
};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Original chemical and thermal failures; neither owner publishes on failure.
#[derive(Clone, Debug, PartialEq)]
pub enum ReactiveThermalPistonError<E> {
    /// Invalid common clock or a failed whole-system balance.
    Invalid(&'static str),
    /// The original kinetic, species or gas-energy owner refused the reaction.
    Reaction(GasReactionError),
    /// The original heat, source-domain or mechanical owner refused a proposal.
    Thermal(ThermalPistonError<E>),
    /// No error-admitted whole-system trajectory fits the explicit work limits.
    SubcyclingLimit {
        evaluations: usize,
        substeps: usize,
        last_error_ratio: Option<f64>,
    },
    /// Cancellation before returning the joint candidate.
    Cancelled,
}
impl<E> From<GasReactionError> for ReactiveThermalPistonError<E> {
    fn from(error: GasReactionError) -> Self {
        if error == GasReactionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Reaction(error)
        }
    }
}
impl<E> From<ThermalPistonError<E>> for ReactiveThermalPistonError<E> {
    fn from(error: ThermalPistonError<E>) -> Self {
        if matches!(error, ThermalPistonError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E: core::fmt::Display> core::fmt::Display for ReactiveThermalPistonError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "reactive gas/body coupling: {message}"),
            Self::Reaction(error) => write!(f, "reactive gas/body coupling: {error}"),
            Self::Thermal(error) => write!(f, "reactive gas/body coupling: {error}"),
            Self::SubcyclingLimit {
                evaluations,
                substeps,
                last_error_ratio,
            } => write!(
                f,
                "reactive gas/body refinement exhausted after {evaluations} steps at mesh {substeps}, error {last_error_ratio:?}"
            ),
            Self::Cancelled => write!(f, "reactive gas/body coupling cancelled"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ReactiveThermalPistonError<E> {}

/// Original stage reports plus an independently measured enclosing balance.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveThermalPistonStep {
    /// One species update at the accepted initial state and fixed total U.
    pub reaction: GasReactionStep,
    /// Existing simultaneous thermal/motion solve starting from reacted gas.
    /// Its local balance covers that stage, not chemical reconstruction roundoff.
    pub thermal: ThermalPistonStep,
    /// Actual gas internal-energy change across BOTH stages [J].
    pub gas_energy_change_j: f64,
    /// dU_gas+dH_body+dK+outside_work-force_work-prescribed_gas_heat
    /// minus prescribed_body_heat [J]. Chemical energy and exchanged Q are
    /// internal transfers; neither is added as another heat source.
    pub total_energy_residual_j: f64,
}

/// All accepted endpoints together, retaining the thermal chart and sources.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveThermalPistonCandidate<S, P = GasPiston> {
    piston: P,
    thermal: UniformEnthalpyStep<S>,
    step: ReactiveThermalPistonStep,
}
impl<S, P> ReactiveThermalPistonCandidate<S, P> {
    /// Joint gas, volume and momentum endpoint with its original body binding.
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    /// Finite material endpoint and actual midpoint interface/source evaluation.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> {
        &self.thermal
    }
    /// Stage accounting and full-system energy residual.
    pub const fn step(&self) -> &ReactiveThermalPistonStep {
        &self.step
    }
    /// Publish both endpoints together in the caller's existing state owner.
    pub fn into_parts(self) -> (P, UniformEnthalpyStep<S>, ReactiveThermalPistonStep) {
        (self.piston, self.thermal, self.step)
    }
}
impl<S> ReactiveThermalPistonCandidate<S> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston) -> P,
    ) -> ReactiveThermalPistonCandidate<S, P> {
        ReactiveThermalPistonCandidate {
            piston: bind(self.piston),
            thermal: self.thermal,
            step: self.step,
        }
    }
}

impl GasPiston {
    /// React once, then solve reciprocal material heat, gas energy and momentum.
    ///
    /// Every source callback sees the reacted composition's actual gas-energy
    /// midpoint and the finite body's enthalpy/phase midpoint. Species sources
    /// therefore determine Cp/Cv, pressure, conductivity, viscosity and acoustic
    /// state together. No independently authored product gas or heat release is
    /// accepted. The callback must be repeatable and cannot debit live state.
    ///
    /// `input.heat_j` and `body.internal_heat_j` are external sources ONLY.
    /// Species formation references already live inside total gas internal
    /// energy; reciprocal wall heat is solved once and cancels in the enclosing
    /// balance. `exchange.energy_tolerance_j` also bounds the complete step.
    /// The chamber heat-port includes exchanged heat and the whole-interval
    /// chemical entropy change, but advances its original clock by only one dt.
    ///
    /// This is first-order reaction/thermal-mechanical splitting. The reaction's
    /// donor-consumption and entropy gates are necessary admission, not temporal
    /// accuracy or stiff chemistry qualification. Refine dt for the reaction and
    /// splitting errors. The body keeps its fixed-volume, negligible-pV-work
    /// enthalpy model; it is not automatically the rigid translating piston.
    /// No flame front, shock, engine cycle, wall deformation or material-law
    /// extrapolation is inferred. All inputs remain unchanged on any refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_candidate<S, E>(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        kinetics: GasReactionSolve,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveThermalPistonCandidate<S>, ReactiveThermalPistonError<E>> {
        super::poll(&mut cancelled)?;
        if body.duration_s != input.duration_s {
            return Err(ReactiveThermalPistonError::Invalid(
                "reaction, body heat and motion need the same physical interval",
            ));
        }
        let (gas, reaction) = reaction
            .parcel_candidate(&self.gas, input.duration_s, kinetics, &mut cancelled)?
            .into_parts();
        let mut reacted = self.clone();
        reacted.gas = gas;
        let (mut piston, thermal, thermal_step) = reacted
            .thermal_candidate(input, body, controls, exchange, resolve, &mut cancelled)?
            .into_parts();
        let gas_energy = self.gas.mass_kg()
            * (piston.gas.state().specific_internal_energy_j_kg()
                - self.gas.state().specific_internal_energy_j_kg());
        let body_energy = body.mass_kg
            * (thermal.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        let motion = &thermal_step.piston;
        let total = ((gas_energy + body_energy) + motion.kinetic_energy_change_j)
            + motion.exterior_pressure_work_j
            - motion.external_force_work_j
            - input.heat_j
            - body.internal_heat_j;
        if !gas_energy.is_finite()
            || !body_energy.is_finite()
            || !total.is_finite()
            || total.abs() > exchange.energy_tolerance_j
        {
            return Err(ReactiveThermalPistonError::Invalid(
                "chemical, finite-body and mechanical endpoints do not close total energy",
            ));
        }
        // The chamber alone is an open thermal subsystem. Its one port retains
        // external heat PLUS Q; Q cancels only when the finite body is enclosed.
        let gas_input = PistonInput {
            heat_j: motion.gas.supplied_heat_j,
            ..input
        };
        self.record_reactive_account(
            &mut piston,
            gas_input,
            motion,
            &reaction,
            gas_energy,
            exchange.energy_tolerance_j,
        )?;
        super::poll(&mut cancelled)?;
        Ok(ReactiveThermalPistonCandidate {
            piston,
            thermal,
            step: ReactiveThermalPistonStep {
                reaction,
                thermal: thermal_step,
                gas_energy_change_j: gas_energy,
                total_energy_residual_j: total,
            },
        })
    }

    /// Reuse the existing material-bound series conduction geometry and law.
    /// Wall conductivity and enthalpy must belong to the same source card;
    /// gas conductivity comes only from the evolving reacted gas midpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_layers_candidate(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        kinetics: GasReactionSolve,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        layers: MaterialConductionLayers<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ReactiveThermalPistonCandidate<ResolvedMaterialStatePoint>,
        ReactiveThermalPistonError<LayerConductionError>,
    > {
        self.reactive_thermal_candidate(
            reaction,
            input,
            body,
            kinetics,
            controls,
            exchange,
            |gas, thermal| layers.resolve(gas, thermal),
            cancelled,
        )
    }
}
