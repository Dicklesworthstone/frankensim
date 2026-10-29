//! Reciprocal heat exchange on the existing finite-gas string transaction.
//! Both endpoint temperatures enter a backward-Euler interface residual. The
//! thermal string, source gas and interface keep their existing physical laws.

use super::{
    FiniteGasStringCandidate, FiniteGasStringError, FiniteGasStringFrame, FiniteGasStringHeat,
    FiniteGasStringSolve,
};
use crate::acoustic_realize::{
    AcousticRealizeError, ThermalMaterialStringRuntime, ThermalStringError,
};
use crate::cavity_phs::piston::{
    PistonError,
    thermal::{
        ConductanceSample, HeatExchangeSolve, ThermalPistonError, feedback,
        layers::{LayerConductionError, MaterialConductionLayers},
    },
};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasHeatStep, IsochoricGasParcel};

/// String admission or the original bounded thermal-feedback owner's refusal.
/// Reusing that owner's numerical error does not introduce piston dynamics.
#[derive(Debug)]
pub enum FiniteGasStringExchangeError<E> {
    /// Original string, gas, epoch or publication refusal.
    String(FiniteGasStringError),
    /// Original bounded exchange, enthalpy or source-law refusal.
    Exchange(ThermalPistonError<E>),
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteGasStringExchangeError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas string exchange: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FiniteGasStringExchangeError<E> {}
impl<E> From<FiniteGasStringError> for FiniteGasStringExchangeError<E> {
    fn from(error: FiniteGasStringError) -> Self {
        Self::String(error)
    }
}
impl<E> From<ThermalStringError> for FiniteGasStringExchangeError<E> {
    fn from(error: ThermalStringError) -> Self {
        Self::String(error.into())
    }
}
impl<E> From<AcousticRealizeError> for FiniteGasStringExchangeError<E> {
    fn from(error: AcousticRealizeError) -> Self {
        Self::from(ThermalStringError::Acoustic(error))
    }
}
impl<E> From<ThermalPistonError<E>> for FiniteGasStringExchangeError<E> {
    fn from(error: ThermalPistonError<E>) -> Self {
        Self::Exchange(error)
    }
}

/// One accepted coupled frame and the actual interface law at its endpoint.
#[derive(Debug)]
pub struct FiniteGasStringExchangeFrame<R> {
    /// Combined account. Internal exchange cancels from this energy residual.
    /// Its solid subaccount includes -Q in `string.external_heat_j`, and its
    /// gas subaccount includes +Q in `gas.supplied_heat_j`.
    pub coupled: FiniteGasStringFrame,
    /// Caller-prescribed solid heat only [J], excluding the internal -Q.
    pub external_string_heat_j: f64,
    /// Signed heat removed from the solid and added to the gas [J].
    pub solid_to_gas_heat_j: f64,
    /// Actual nonnegative endpoint interface conductance [W/K].
    pub conductance_w_k: f64,
    /// Q - dt G(T_solid,end - T_gas,end) [J].
    pub source_residual_j: f64,
    /// Caloric solid + gas storage minus their non-exchange heat sources [J],
    /// including the signed thermoelastic heat supplied to the solid.
    pub thermal_energy_residual_j: f64,
    /// Concrete source evaluation retained from the accepted interface trial.
    pub source: R,
    /// Complete interface/source proposals over all thermoelastic trials.
    /// Each thermal trial uses at most exchange.maximum_iterations + 2;
    /// the outer reversible-heat solve has its own explicit work ceiling.
    pub trials: usize,
}

struct ExchangeCandidate<R, G> {
    gas: IsochoricGasParcel,
    gas_step: GasHeatStep,
    gas_report: G,
    solid: EquilibriumPhaseState,
    heat_j: f64,
    sample: ConductanceSample<R>,
    source_residual_j: f64,
    energy_residual_j: f64,
    trials: usize,
}
fn gas_error<E>(error: GasEnergyError) -> ThermalPistonError<E> {
    PistonError::from(error).into()
}

fn feedback_error<E>(
    error: feedback::ScalarError<FiniteGasStringExchangeError<E>>,
) -> FiniteGasStringExchangeError<E> {
    use feedback::ScalarError;
    let error = match error {
        ScalarError::Evaluation(error) => return error,
        ScalarError::Invalid(_) => ThermalPistonError::Invalid("nonfinite physical energy feedback residual"),
        ScalarError::Cancelled => ThermalPistonError::Cancelled,
        ScalarError::OutsideBracket { lower, upper } => ThermalPistonError::OutsideBracket {
            lower_residual_j: lower, upper_residual_j: upper,
        },
        ScalarError::NotConverged { iterations, residual } => ThermalPistonError::NotConverged {
            iterations, residual_j: residual,
        },
    };
    error.into()
}

pub(super) fn resolve_conduction(
    layers: &MaterialConductionLayers<'_>,
    volume: f64,
    gas: &ResolvedThermochemicalGasState,
    solid: EquilibriumPhaseState,
) -> Result<ConductanceSample<ResolvedMaterialStatePoint>, LayerConductionError> {
    let sample = layers.resolve(gas, solid)?;
    let conductivity = sample.source.property_by_key(layers.conductivity)
        .ok_or(LayerConductionError::Invalid("resolved layer conductivity is absent"))?
        .value_si();
    let biot = if sample.conductance_w_k == 0.0 { 0.0 }
        else { sample.conductance_w_k * (volume / layers.area_m2) / layers.area_m2 / conductivity };
    if !biot.is_finite() || biot > fs_conduction::lumped::LUMPED_BIOT_CEILING {
        return Err(LayerConductionError::Invalid(
            "source conduction exceeds the lumped solid Biot limit; use spatial transport",
        ));
    }
    Ok(sample)
}

impl ThermalMaterialStringRuntime {
    /// Couple the same string and finite gas through a source-bearing heat law.
    ///
    /// Mechanics freezes the initial material viscosity/gas drag for one exact
    /// sample. Material loss, reversible thermoelastic exchange and prescribed
    /// solid heat enter h_s; aerodynamic loss and prescribed gas heat enter u_g.
    /// One implicit transfer Q then
    /// satisfies Q = dt G(T_s,end,T_g,end)(T_s,end-T_g,end), using the original
    /// solid enthalpy chart and gas internal-energy inversion at every trial.
    /// No midpoint temperatures or constant heat capacities are substituted.
    ///
    /// `resolve` receives actual ENDPOINT states and returns the existing
    /// `ConductanceSample`, retaining its concrete source. It must be repeatable
    /// and proposal-only. Constant, contact, conduction or convection laws are
    /// valid only when their own source, geometry and lumped assumptions admit
    /// them. The concrete `step_with_isochoric_gas_conduction` path supplies
    /// the existing source-bound series-layer law and a pointwise solid Biot gate.
    ///
    /// The bracket includes zero, and every evaluated trial must be supported
    /// by both original sources. A source-only state or bracket probe outside
    /// a chart refuses; there is no extrapolation or unadvertised substep.
    /// The existing bounded feedback owner controls iterations and absolute
    /// source/thermal-energy residuals. Nonnegative G and a final downhill
    /// transfer prevent endpoint temperature crossing by the accepted exchange.
    /// Nonmonotone state-dependent G is not given a uniqueness guarantee.
    ///
    /// The same final gas observer, epoch, independent combined energy check
    /// and atomic publication are used as in `step_with_isochoric_gas`.
    /// The mechanical/thermal split remains first order; this is not a spatial
    /// boundary layer, chamber resonance or gas pressure-work model.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_exchange<R, E>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<R>, E>,
    ) -> Result<FiniteGasStringExchangeFrame<R>, FiniteGasStringExchangeError<E>> {
        self.step_with_isochoric_gas_exchange_endpoint(
            cx, expected_epoch, forces, gas, heat, controls, exchange, resolve,
            |initial, heat_j, _, cancelled| {
                let mut next = initial.clone();
                let step = next.add_heat(heat_j, controls.gas, cancelled).map_err(gas_error)?;
                Ok((next, step, ()))
            },
        ).map(|(frame, ())| frame)
    }

    // One thermal proposal/feedback/publication owner for both frozen and
    // reacting gas. Every endpoint is proposed from the accepted inventory;
    // rejected heat trials cannot consume species or advance mechanical state.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_with_isochoric_gas_exchange_endpoint<R, E, G>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<R>, E>,
        mut endpoint: impl FnMut(
            &IsochoricGasParcel, f64, f64, &mut dyn FnMut() -> bool,
        ) -> Result<(IsochoricGasParcel, GasHeatStep, G), FiniteGasStringExchangeError<E>>,
    ) -> Result<(FiniteGasStringExchangeFrame<R>, G), FiniteGasStringExchangeError<E>> {
        self.admit_finite_gas_step(cx, expected_epoch, gas, heat, controls)?;
        feedback::admit::<E>(exchange)?;
        let mut mechanical = self.mechanical.clone();
        let vibration = mechanical.step(cx, forces)?;
        let gas_source = heat.gas_heat_j + vibration.dissipation.air_loss_j;
        if !gas_source.is_finite()
            || (vibration.dissipation.air_loss_j != 0.0 && gas_source == heat.gas_heat_j)
            || (heat.gas_heat_j != 0.0 && gas_source == vibration.dissipation.air_loss_j)
        {
            return Err(FiniteGasStringError::Admission(
                "gas thermal sources overflow or disappear when combined",
            )
            .into());
        }
        let mut total_trials = 0usize;
        // Every reversible-heat proposal restarts this SAME finite-gas exchange
        // from the accepted reservoirs. Neither the parcel nor the string is
        // published until both nested physical balances have been accepted.
        let proposed = self.propose_thermal_step(
            cx,
            &mechanical,
            &vibration,
            |input| -> Result<_, FiniteGasStringExchangeError<E>> {
                let solid_source = input.internal_heat_j + heat.string_heat_j;
                if !solid_source.is_finite() {
                    return Err(FiniteGasStringError::Admission(
                        "solid thermal sources overflow when combined",
                    ).into());
                }
                let solved = feedback::solve_scalar(
                    exchange.heat_bracket_j,
                    exchange.source_tolerance_j,
                    exchange.maximum_iterations,
                    |q, trials, cancelled| -> Result<_, FiniteGasStringExchangeError<E>> {
                        let solid_heat = solid_source - q;
                        let gas_heat = gas_source + q;
                        if !solid_heat.is_finite()
                            || !gas_heat.is_finite()
                            || (q != 0.0 && (solid_heat == solid_source || gas_heat == gas_source))
                        {
                            return Err(ThermalPistonError::Invalid(
                                "exchange heat overflows or disappears in either reservoir",
                            ).into());
                        }
                        let solid = input.curve
                            .advance_specific_energy(input.initial, solid_heat / input.mass_kg)
                            .map_err(ThermalPistonError::from)?;
                        let (next, gas_step, gas_report) = endpoint(
                            gas, gas_heat, input.duration_s, cancelled,
                        )?;
                        let sample = resolve(next.state(), solid).map_err(ThermalPistonError::Interface)?;
                        let g = sample.conductance_w_k;
                        let delta = solid.temperature_k() - next.state().state().temperature;
                        let power = g * delta;
                        let source_heat = input.duration_s * power;
                        let residual = q - source_heat;
                        let solid_storage = input.mass_kg
                            * (solid.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
                        let solid_residual = solid_storage - solid_source + q;
                        let energy_residual =
                            (solid_storage + gas_step.stored_energy_change_j) - solid_source - gas_source;
                        if !g.is_finite()
                            || g < 0.0
                            || !power.is_finite()
                            || !source_heat.is_finite()
                            || !residual.is_finite()
                            || (g != 0.0 && delta != 0.0 && power == 0.0)
                            || (power != 0.0 && source_heat == 0.0)
                        {
                            return Err(ThermalPistonError::Invalid(
                                "endpoint interface heat law is negative, nonfinite or unrepresentable",
                            ).into());
                        }
                        if !solid_residual.is_finite()
                            || !energy_residual.is_finite()
                            || solid_residual.abs() > exchange.energy_tolerance_j
                            || energy_residual.abs() > exchange.energy_tolerance_j
                        {
                            return Err(ThermalPistonError::Invalid(
                                "finite solid/gas thermal storage does not close its energy allowance",
                            ).into());
                        }
                        Ok((
                            ExchangeCandidate {
                                gas: next,
                                gas_step,
                                gas_report,
                                solid,
                                heat_j: q,
                                sample,
                                source_residual_j: residual,
                                energy_residual_j: energy_residual,
                                trials,
                            },
                            residual,
                        ))
                    },
                    &mut || cx.checkpoint().is_err(),
                ).map_err(feedback_error)?;
                let difference = solved.solid.temperature_k() - solved.gas.state().state().temperature;
                if (solved.heat_j > 0.0 && difference < 0.0)
                    || (solved.heat_j < 0.0 && difference > 0.0)
                {
                    return Err(FiniteGasStringError::Admission(
                        "accepted heat exchange would cross endpoint temperatures; tighten the source tolerance",
                    ).into());
                }
                total_trials += solved.trials;
                Ok(UniformEnthalpyStep {
                    state: solved.solid,
                    external_heat_j: heat.string_heat_j - solved.heat_j,
                    energy_residual_tolerance_j: exchange.energy_tolerance_j,
                    report: solved,
                })
            },
        )?;
        let UniformEnthalpyStep {
            state: solid,
            external_heat_j,
            energy_residual_tolerance_j,
            report: solved,
        } = proposed;
        let ExchangeCandidate {
            gas: next_gas,
            gas_step,
            gas_report,
            heat_j: q,
            sample,
            source_residual_j,
            energy_residual_j,
            ..
        } = solved;
        let mut candidate = self.clone();
        let string = candidate
            .finish_thermal_step(
                cx,
                forces,
                mechanical,
                vibration,
                UniformEnthalpyStep {
                    state: solid,
                    external_heat_j,
                    energy_residual_tolerance_j,
                    report: (),
                },
            )?
            .coupled;
        let coupled = self.publish_finite_gas(
            cx,
            forces,
            gas,
            heat,
            controls,
            FiniteGasStringCandidate {
                runtime: candidate,
                gas: next_gas,
                string,
                gas_step,
            },
        )?;
        Ok((FiniteGasStringExchangeFrame {
            coupled,
            external_string_heat_j: heat.string_heat_j,
            solid_to_gas_heat_j: q,
            conductance_w_k: sample.conductance_w_k,
            source_residual_j,
            thermal_energy_residual_j: energy_residual_j,
            source: sample.source,
            trials: total_trials,
        }, gas_report))
    }

    /// Source-bound series conduction using the SAME solid card/conditions and
    /// the gas endpoint's actual conductivity. Layer distances are explicit
    /// effective lengths; the transfer area cannot exceed the solid's surface.
    ///
    /// The existing layer evaluator is supplied endpoint states on this
    /// backward-Euler path. Its source and property-key checks remain unchanged.
    /// Every trial uses the existing lumped convention Bi = G(V/A)/(A k_s)
    /// with A equal to the declared transfer area and the existing 0.1 ceiling.
    /// This pointwise gate is not a bound on
    /// unresolved internal gradients, local spreading or a changing phase.
    /// Zero area or explicitly zero solid conductivity is insulation, while
    /// source binding is still checked. No missing coefficient is filled in.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
    ) -> Result<
        FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        let volume = self.admit_finite_gas_conduction(cx, layers)?;
        self.step_with_isochoric_gas_exchange(
            cx, expected_epoch, forces, gas, heat, controls, exchange,
            |gas_end, solid_end| resolve_conduction(layers, volume, gas_end, solid_end),
        )
    }

    pub(super) fn admit_finite_gas_conduction(
        &self,
        cx: &fs_exec::Cx<'_>,
        layers: &MaterialConductionLayers<'_>,
    ) -> Result<f64, FiniteGasStringExchangeError<LayerConductionError>> {
        super::check(cx)?;
        let specimen = self.mechanical.specimen();
        let radius = specimen.radius_m();
        let length = specimen.string().length_m;
        let area = 2.0 * core::f64::consts::PI * radius * (length + radius);
        let volume = specimen.area_m2() * length;
        let same_conditions = layers
            .context
            .axes()
            .iter()
            .filter(|(axis, _)| axis.as_str() != "T")
            .map(|(axis, value)| (axis, value.to_bits()))
            .eq(self
                .point
                .axes()
                .iter()
                .filter(|(axis, _)| axis.as_str() != "T")
                .map(|(axis, value)| (axis, value.to_bits())));
        if layers.material.content_hash() != self.card.content_hash()
            || layers.temperature_axis != "T"
            || layers.liquid_fraction_axis.is_some()
            || !same_conditions
            || layers.context.axis_quantities() != self.point.axis_quantities()
            || !layers.area_m2.is_finite()
            || layers.area_m2 < 0.0
            || layers.area_m2 > area
        {
            return Err(FiniteGasStringError::Admission(
                "conduction layers must retain the solid card/conditions, T axis and an actual surface area",
            ).into());
        }
        Ok(volume)
    }
}
