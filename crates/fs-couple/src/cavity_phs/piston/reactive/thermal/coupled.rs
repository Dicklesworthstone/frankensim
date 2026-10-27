//! Same-endpoint reaction, piston work and reciprocal finite-wall heat.
//!
//! A single signed Q is supplied to the original coupled reacting piston and
//! removed from its finite counterpart. The shared scalar residual evaluates
//! conductance at the actual endpoint species, gas T/p and wall enthalpy/phase.
//! Every trial starts at the accepted histories; neither chemical formation
//! energy nor an earlier trial's reaction progress is applied as extra heat.
//! This is a first-order endpoint heat/kinetics scheme with midpoint mechanical
//! kinematics, not an exact reacting-flow integrator or resolved wall profile.

use super::super::coupled::{
    CoupledReactivePistonError, CoupledReactivePistonSolve, CoupledReactivePistonStep,
};
use super::super::{GasReaction, GasReactionSolve};
use crate::cavity_phs::piston::thermal::boundary::heat::conductive_heat;
use crate::cavity_phs::piston::thermal::boundary::ports::{ThermalPort, ThermalPortSample};
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError, feedback,
};
use crate::cavity_phs::piston::{GasPiston, PistonInput};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use std::collections::BTreeMap;

/// No source, chemical, mechanical, heat or cancellation refusal publishes state.
#[derive(Clone, Debug, PartialEq)]
pub enum CoupledReactiveThermalError<E> {
    Invalid(&'static str),
    Motion(CoupledReactivePistonError),
    Thermal(ThermalPistonError<E>),
    Cancelled,
}
impl<E> From<CoupledReactivePistonError> for CoupledReactiveThermalError<E> {
    fn from(error: CoupledReactivePistonError) -> Self {
        match error {
            CoupledReactivePistonError::Cancelled => Self::Cancelled,
            other => Self::Motion(other),
        }
    }
}
impl<E> From<ThermalPistonError<E>> for CoupledReactiveThermalError<E> {
    fn from(error: ThermalPistonError<E>) -> Self {
        match error {
            ThermalPistonError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for CoupledReactiveThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled reacting gas and wall: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for CoupledReactiveThermalError<E> {}

/// Original gas chemistry and motion plus the independent boundary heat law.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactiveBoundaryStep {
    pub piston: CoupledReactivePistonStep,
    /// Signed heat supplied to gas, to be debited exactly once from counterparts.
    pub exchanged_heat_j: f64,
    /// Sum of actual endpoint contact heat laws, before any residual admission.
    pub evaluated_heat_j: f64,
    pub source_residual_j: f64,
    /// Gas/kinetic change and external work minus prescribed and evaluated heat.
    pub total_energy_residual_j: f64,
    /// Complete pressure/chemistry solves, including discarded heat trials.
    pub exchange_trials: usize,
    /// All pressure proposals nested inside those solves.
    pub pressure_trials: usize,
    /// All kinetic endpoint/Jacobian evaluations nested inside those proposals.
    pub chemical_evaluations: usize,
}

/// Gas and all contact source records from one proposed interval.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactiveBoundaryCandidate<S, P = GasPiston> {
    piston: P,
    source: S,
    step: CoupledReactiveBoundaryStep,
}
impl<S, P> CoupledReactiveBoundaryCandidate<S, P> {
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    pub const fn source(&self) -> &S {
        &self.source
    }
    pub const fn step(&self) -> &CoupledReactiveBoundaryStep {
        &self.step
    }
    pub fn into_parts(self) -> (P, S, CoupledReactiveBoundaryStep) {
        (self.piston, self.source, self.step)
    }
}

/// Enclosing gas/kinetic/body energy account. No independent body entropy model.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactiveThermalStep {
    pub boundary: CoupledReactiveBoundaryStep,
    /// Actual source-resolved endpoint interface conductance [W/K].
    pub conductance_w_k: f64,
    pub thermal_energy_change_j: f64,
    pub thermal_energy_residual_j: f64,
    /// Internal Q cancels; only the two prescribed heats and mechanical work remain.
    pub total_energy_residual_j: f64,
}

/// Both continuing histories; consuming the report must not deposit Q again.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactiveThermalCandidate<S, P = GasPiston> {
    piston: P,
    thermal: UniformEnthalpyStep<S>,
    step: CoupledReactiveThermalStep,
}
impl<S, P> CoupledReactiveThermalCandidate<S, P> {
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> {
        &self.thermal
    }
    pub const fn step(&self) -> &CoupledReactiveThermalStep {
        &self.step
    }
    pub fn into_parts(self) -> (P, UniformEnthalpyStep<S>, CoupledReactiveThermalStep) {
        (self.piston, self.thermal, self.step)
    }
}
impl<S> CoupledReactiveThermalCandidate<S> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston) -> P,
    ) -> CoupledReactiveThermalCandidate<S, P> {
        CoupledReactiveThermalCandidate {
            piston: bind(self.piston),
            thermal: self.thermal,
            step: self.step,
        }
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), CoupledReactiveThermalError<E>> {
    if cancelled() {
        Err(CoupledReactiveThermalError::Cancelled)
    } else {
        Ok(())
    }
}

impl GasPiston {
    /// Close reaction progress, pressure work and EACH thermal contact together.
    ///
    /// Port temperatures belong to the same external thermal trial. Each source
    /// sees the actual reacted gas endpoint. Hot and cold contacts retain their
    /// signed transfers even when their sum is zero. A finite receiving network
    /// must debit every returned heat and publish all endpoints together.
    /// Source callbacks must be repeatable and leave accepted state untouched.
    ///
    /// Requires the native coupled implicit chemistry policy. The pressure,
    /// heat and reaction budgets bound every nested solve. No source selection,
    /// temperature, mixture transport or reaction rate is silently frozen.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_ports_candidate<S, E>(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        chemistry: GasReactionSolve,
        controls: CoupledReactivePistonSolve,
        exchange: HeatExchangeSolve,
        ports: &[ThermalPort],
        maximum_ports: usize,
        mut resolve: impl FnMut(
            ThermalPort,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        CoupledReactiveBoundaryCandidate<Vec<ThermalPortSample<S>>>,
        CoupledReactiveThermalError<E>,
    > {
        use CoupledReactiveThermalError::Invalid;
        poll(&mut cancelled)?;
        if ports.len() > maximum_ports {
            return Err(Invalid("reacting gas thermal port budget exceeded"));
        }
        let mut ordered = BTreeMap::new();
        for port in ports {
            poll(&mut cancelled)?;
            if !positive(port.temperature_k) || ordered.insert(port.id, *port).is_some() {
                return Err(Invalid(
                    "unique thermal port IDs and positive endpoint temperatures required",
                ));
            }
        }
        solve(
            self,
            reaction,
            input,
            chemistry,
            controls,
            exchange,
            |_, gas, cancel| {
                let mut total = 0.0;
                let mut samples = Vec::with_capacity(ordered.len());
                for &port in ordered.values() {
                    poll(cancel)?;
                    let sample = resolve(port, gas).map_err(ThermalPistonError::Interface)?;
                    let heat_j = conductive_heat(
                        input.duration_s,
                        sample.conductance_w_k,
                        port.temperature_k,
                        gas.state().temperature,
                    )?;
                    total += heat_j;
                    if !total.is_finite() {
                        return Err(Invalid("summed reactive contact heat is unrepresentable"));
                    }
                    samples.push(ThermalPortSample {
                        port,
                        conductance_w_k: sample.conductance_w_k,
                        heat_j,
                        source: sample.source,
                    });
                }
                Ok((total, samples))
            },
            &mut cancelled,
        )
    }

    /// Evolve a finite material enthalpy store and reacting moving gas together.
    ///
    /// The actual ENDPOINT gas and body states determine G and Q=dt G(Tb-Tg).
    /// Formation energy is already in the gas internal energy. Body latent heat
    /// is already in its enthalpy chart. Neither is an extra heat source.
    /// The material body retains fixed geometry and negligible pV-work scope;
    /// it need not be the translating piston. Material phase change does not
    /// silently change that piston's inertia, shape or mechanical support.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_candidate<S, E>(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        chemistry: GasReactionSolve,
        controls: CoupledReactivePistonSolve,
        exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledReactiveThermalCandidate<S>, CoupledReactiveThermalError<E>> {
        use CoupledReactiveThermalError::Invalid;
        poll(&mut cancelled)?;
        if !positive(body.mass_kg)
            || !positive(body.volume_m3)
            || !positive(body.surface_area_m2)
            || !positive(body.duration_s)
            || body.duration_s != input.duration_s
            || !body.internal_heat_j.is_finite()
        {
            return Err(Invalid(
                "matching positive clocks, body geometry/mass and finite body heat required",
            ));
        }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body
            .curve
            .state_at_specific_enthalpy(h0)
            .map_err(ThermalPistonError::Thermal)?
            != body.initial
        {
            return Err(Invalid(
                "accepted body state must belong to its original enthalpy chart",
            ));
        }
        let candidate = solve(
            self,
            reaction,
            input,
            chemistry,
            controls,
            exchange,
            |q, gas, cancel| {
                poll(cancel)?;
                let body_heat = body.internal_heat_j - q;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                if !body_heat.is_finite()
                    || !dh.is_finite()
                    || !h1.is_finite()
                    || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0 && body_heat == body.internal_heat_j)
                {
                    return Err(Invalid(
                        "reactive wall exchange vanishes or overflows in enthalpy",
                    ));
                }
                let state = body
                    .curve
                    .state_at_specific_enthalpy(h1)
                    .map_err(ThermalPistonError::Thermal)?;
                let sample = resolve(gas, state).map_err(ThermalPistonError::Interface)?;
                let evaluated = conductive_heat(
                    input.duration_s,
                    sample.conductance_w_k,
                    state.temperature_k(),
                    gas.state().temperature,
                )?;
                let stored = body.mass_kg * (state.specific_enthalpy_j_kg() - h0);
                let residual = stored + q - body.internal_heat_j;
                if !stored.is_finite()
                    || !residual.is_finite()
                    || residual.abs() > exchange.energy_tolerance_j
                {
                    return Err(Invalid(
                        "reacting wall enthalpy and its actual heat debit do not close",
                    ));
                }
                Ok((
                    evaluated,
                    (
                        UniformEnthalpyStep {
                            state,
                            external_heat_j: -q,
                            energy_residual_tolerance_j: exchange.energy_tolerance_j,
                            report: sample.source,
                        },
                        stored,
                        residual,
                        sample.conductance_w_k,
                    ),
                ))
            },
            &mut cancelled,
        )?;
        let (piston, (thermal, stored, residual, conductance_w_k), boundary) =
            candidate.into_parts();
        let motion = &boundary.piston;
        let total = motion.reaction.reaction.stored_energy_change_j
            + motion.kinetic_energy_change_j
            + stored
            + motion.exterior_pressure_work_j
            - motion.external_force_work_j
            - motion.support_work_j
            - input.heat_j
            - body.internal_heat_j;
        if !total.is_finite() || total.abs() > exchange.energy_tolerance_j {
            return Err(Invalid(
                "reacting gas, piston and finite-wall energy do not close",
            ));
        }
        poll(&mut cancelled)?;
        Ok(CoupledReactiveThermalCandidate {
            piston,
            thermal,
            step: CoupledReactiveThermalStep {
                boundary,
                conductance_w_k,
                thermal_energy_change_j: stored,
                thermal_energy_residual_j: residual,
                total_energy_residual_j: total,
            },
        })
    }

    /// Same-card source conductivity and enthalpy, with gas transport resolved
    /// from the SAME reacted endpoint used for pressure and chemical rates.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_layers_candidate(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        chemistry: GasReactionSolve,
        controls: CoupledReactivePistonSolve,
        exchange: HeatExchangeSolve,
        layers: MaterialConductionLayers<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        CoupledReactiveThermalCandidate<ResolvedMaterialStatePoint>,
        CoupledReactiveThermalError<LayerConductionError>,
    > {
        if layers.area_m2 > body.surface_area_m2 {
            return Err(CoupledReactiveThermalError::Invalid(
                "conducting patch exceeds the finite body's exposed surface",
            ));
        }
        self.coupled_reactive_thermal_candidate(
            reaction,
            input,
            body,
            chemistry,
            controls,
            exchange,
            |gas, state| layers.resolve(gas, state),
            cancelled,
        )
    }
}

// One heat-root owner for finite bodies and arbitrary many-contact networks.
#[allow(clippy::too_many_arguments)]
fn solve<S, E, C: FnMut() -> bool>(
    piston: &GasPiston,
    reaction: &GasReaction,
    input: PistonInput,
    chemistry: GasReactionSolve,
    controls: CoupledReactivePistonSolve,
    exchange: HeatExchangeSolve,
    mut resolve: impl FnMut(
        f64,
        &ResolvedThermochemicalGasState,
        &mut C,
    ) -> Result<(f64, S), CoupledReactiveThermalError<E>>,
    cancelled: &mut C,
) -> Result<CoupledReactiveBoundaryCandidate<S>, CoupledReactiveThermalError<E>> {
    use CoupledReactiveThermalError::{Invalid, Thermal};
    feedback::admit(exchange)?;
    if !input.heat_j.is_finite() {
        return Err(Invalid("finite prescribed gas heat required"));
    }
    let mut pressure_trials = 0usize;
    let mut chemical_evaluations = 0usize;
    let candidate = feedback::solve_scalar(
        exchange.heat_bracket_j,
        exchange.source_tolerance_j,
        exchange.maximum_iterations,
        |q, trials, cancel| {
            let heat_j = input.heat_j + q;
            if !heat_j.is_finite() || (q != 0.0 && heat_j == input.heat_j) {
                return Err(Invalid(
                    "reacting boundary heat vanishes or overflows in gas energy",
                ));
            }
            let proposal = piston.coupled_reactive_candidate(
                reaction,
                PistonInput { heat_j, ..input },
                chemistry,
                controls,
                &mut *cancel,
            )?;
            pressure_trials = pressure_trials
                .checked_add(proposal.step().trials)
                .ok_or(Invalid("reacting pressure work count overflowed"))?;
            chemical_evaluations = chemical_evaluations
                .checked_add(proposal.step().chemical_evaluations)
                .ok_or(Invalid("chemical work count overflowed"))?;
            poll(cancel)?;
            let (evaluated, source) = resolve(q, proposal.piston().gas().state(), cancel)?;
            let residual = q - evaluated;
            let (piston, motion) = proposal.into_parts();
            let total = motion.reaction.reaction.stored_energy_change_j
                + motion.kinetic_energy_change_j
                + motion.exterior_pressure_work_j
                - motion.external_force_work_j
                - motion.support_work_j
                - input.heat_j
                - evaluated;
            if !evaluated.is_finite() || !residual.is_finite() || !total.is_finite() {
                return Err(Invalid(
                    "reacting boundary source or work is unrepresentable",
                ));
            }
            Ok((
                CoupledReactiveBoundaryCandidate {
                    piston,
                    source,
                    step: CoupledReactiveBoundaryStep {
                        piston: motion,
                        exchanged_heat_j: q,
                        evaluated_heat_j: evaluated,
                        source_residual_j: residual,
                        total_energy_residual_j: total,
                        exchange_trials: trials,
                        pressure_trials,
                        chemical_evaluations,
                    },
                },
                residual,
            ))
        },
        cancelled,
    )
    .map_err(|error| match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Cancelled => CoupledReactiveThermalError::Cancelled,
        feedback::ScalarError::Invalid(message) => Invalid(message),
        feedback::ScalarError::OutsideBracket { lower, upper } => {
            Thermal(ThermalPistonError::OutsideBracket {
                lower_residual_j: lower,
                upper_residual_j: upper,
            })
        }
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => Thermal(ThermalPistonError::NotConverged {
            iterations,
            residual_j: residual,
        }),
    })?;
    if candidate.step.total_energy_residual_j.abs() > exchange.energy_tolerance_j {
        return Err(Invalid(
            "evaluated reactive boundary heat does not close gas and mechanical work",
        ));
    }
    poll(cancelled)?;
    Ok(candidate)
}
