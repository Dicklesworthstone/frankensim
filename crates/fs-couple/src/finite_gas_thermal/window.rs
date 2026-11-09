//! Adaptive finite-gas cooling/heating with original sensible and latent storage.
//!
//! Compare one native backward-Euler interval with two native half intervals;
//! accept only the fine states. Rejected heat, entropy and elapsed time never
//! enter the physical account. This is local time-error control, not a global
//! accuracy bound or a new material, phase, convection or gas model.

use super::*;
use fs_material::phase::EquilibriumPhaseState;
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

/// Whole-window work and conservation limits plus local physical-time accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasWindowSolve {
    /// Initial and maximum full comparison span. Accepted halves are shorter.
    pub maximum_step_s: f64,
    /// Minimum retried comparison span; a final remainder may be shorter.
    pub minimum_step_s: f64,
    /// All native interval calls, including coarse, fine and refused proposals.
    /// Each call retains the original exchange and nested material work limits.
    pub maximum_attempts: usize,
    /// Applies separately to body enthalpy and fixed-density gas internal energy.
    pub energy_absolute_tolerance_j: f64,
    /// Applies separately to body and gas temperatures.
    pub temperature_absolute_tolerance_k: f64,
    pub phase_fraction_absolute_tolerance: f64,
    /// Scales each coordinate's CHANGE, never an absolute formation datum or T.
    pub relative_tolerance: f64,
    /// Whole-window storage, source-apportionment and absolute local defects [J].
    pub energy_balance_tolerance_j: f64,
}

/// Compact accepted history, retaining convection sources without gas catalogs.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasWindowSample<S> {
    pub time_s: f64,
    pub body_state: EquilibriumPhaseState,
    pub gas_temperature_k: f64,
    pub convection: ConvectionCoefficient<S>,
    pub solid_to_gas_heat_j: f64,
    pub radiation_heat_j: f64,
    pub entropy: FiniteGasEntropyReport,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasWindowReport<S, G> {
    pub next_gas: G,
    /// Every accepted half-step endpoint, in physical-time order.
    pub samples: Vec<FiniteGasWindowSample<S>>,
    pub attempted_intervals: usize,
    pub rejected_comparisons: usize,
    pub maximum_accepted_normalized_error: f64,
    pub solid_to_gas_heat_j: f64,
    pub radiation_heat_j: f64,
    pub internal_body_heat_j: f64,
    pub external_gas_heat_j: f64,
    pub body_heat_apportionment_residual_j: f64,
    pub gas_heat_apportionment_residual_j: f64,
    pub body_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub combined_energy_residual_j: f64,
    pub sum_absolute_energy_defect_j: f64,
    /// Sum of the original accepted entropy accounts, including temporal defects.
    pub entropy: FiniteGasEntropyReport,
}

#[derive(Debug)]
pub enum FiniteGasWindowError<E> {
    Invalid(&'static str),
    Step(ConvectionError<E>),
    Cancelled,
    WorkLimit { attempted_intervals: usize, accepted_time_s: f64 },
    TimeAccuracy { accepted_time_s: f64, span_s: f64, normalized_error: f64 },
    EnergyBalance { residual_j: f64, absolute_defect_j: f64, tolerance_j: f64 },
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteGasWindowError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite gas/body thermal window: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FiniteGasWindowError<E> {}
impl<E> From<ConvectionError<E>> for FiniteGasWindowError<E> {
    fn from(error: ConvectionError<E>) -> Self {
        match error {
            ConvectionError::Thermal(Error::Cancelled) => Self::Cancelled,
            error => Self::Step(error),
        }
    }
}
impl<E> From<Error> for FiniteGasWindowError<E> {
    fn from(error: Error) -> Self { ConvectionError::Thermal(error).into() }
}

#[derive(Clone, Copy, Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), FiniteGasWindowError<E>> {
        let (next, residual) = fs_math::eft::two_sum(self.value, value);
        self.value = next;
        self.correction += residual;
        finite(self.get()).map(|_| ())
    }
    fn get(self) -> f64 { self.value + self.correction }
}
fn finite<E>(value: f64) -> Result<f64, FiniteGasWindowError<E>> {
    if value.is_finite() { Ok(value) }
    else { Err(FiniteGasWindowError::Invalid("unrepresentable thermal window arithmetic")) }
}
fn poll<E>(cx: &Cx<'_>) -> Result<(), FiniteGasWindowError<E>> {
    cx.checkpoint().map_err(|_| FiniteGasWindowError::Cancelled)
}
fn source<E>(total: f64, start: f64, end: f64, duration: f64)
    -> Result<f64, FiniteGasWindowError<E>> {
    let at = |time| if time == duration { total } else { total * (time / duration) };
    let value = finite(at(end) - at(start))?;
    if total != 0.0 && value == 0.0 {
        return Err(FiniteGasWindowError::Invalid("nonzero interval heat is unrepresentable"));
    }
    Ok(value)
}

impl<G: IsochoricGasStorage> FiniteGasThermalEnvironment<'_, G> {
    /// Advance a complete physical window with the declared constant coefficient.
    pub fn advance_adaptive(
        &self, cx: &Cx<'_>, input: UniformEnthalpyStepInput<'_>, controls: FiniteGasWindowSolve,
    ) -> Result<UniformEnthalpyStep<FiniteGasWindowReport<(), G>>,
        FiniteGasWindowError<core::convert::Infallible>> {
        self.advance_adaptive_with_convection(cx, input, controls, |_| {
            Ok(ConvectionCoefficient {
                coefficient_w_per_m2_k: self.body.convection_w_per_m2_k, source: (),
            })
        })
    }

    /// Adapt physical time using native body h/phase and fixed-density gas Δu.
    ///
    /// `gas_heat_j` and `input.internal_heat_j` are TOTAL window energies,
    /// supplied at constant power. Cumulative source endpoints apportion them
    /// once, including the final remainder. A new authored source/environment
    /// window is required at a physical discontinuity. Geometry and mass stay
    /// fixed; the original body phase curve and gas source laws never change.
    ///
    /// Only exchange-bracket/increment admission failures may shorten a trial.
    /// Material, convection, nonlinear, conservation and cancellation errors
    /// propagate. Callbacks must be repeatable; coarse/fine/rejected evaluations
    /// cannot roll back their external side effects. Neither accepted physical
    /// input changes, even if work or accuracy fails late in the window.
    /// Native algebraic tolerances should be smaller than the time-error targets.
    pub fn advance_adaptive_with_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: FiniteGasWindowSolve,
        mut resolve: impl FnMut(&G::State) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<UniformEnthalpyStep<FiniteGasWindowReport<S, G>>, FiniteGasWindowError<E>> {
        self.advance_adaptive_with_interval(cx, input, controls, |environment, interval| {
            environment.advance_with_convection(cx, interval, &mut resolve)
        })
    }

    /// Adapt physical time while resolving the actual wall and gas endpoint
    /// together in EVERY native coarse, fine and retried interval.
    ///
    /// The original wall transaction retains both its queried and endpoint
    /// source samples in each accepted half step. Its temperature and boundary
    /// work tolerances remain independent of the local time-error targets.
    /// Source energies are apportioned once by the SAME window controller.
    /// Original wall/source failures propagate; only the existing outer heat
    /// bracket and gas-increment admissions may cause an interval retry.
    pub fn advance_adaptive_with_wall_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: FiniteGasWindowSolve,
        wall: WallConvectionSolve,
        mut resolve: impl FnMut(&G::State, f64) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<
        UniformEnthalpyStep<FiniteGasWindowReport<WallConvectionSample<S>, G>>,
        FiniteGasWindowError<WallConvectionError<E>>,
    > {
        self.advance_adaptive_with_wall_convection_checked(
            cx,
            input,
            controls,
            wall,
            |gas, temperature, _| resolve(gas, temperature),
            |_, _| Ok(()),
        )
    }

    /// Source-adapter variant retaining interval-specific work and state gates.
    ///
    /// The resolver also receives the actual native interval input, including
    /// its duration and apportioned body heat. The admission callback sees the
    /// original initial and actual final caloric states of EVERY proposal,
    /// including coarse and subsequently rejected intervals. It can refuse a
    /// source/phase/domain exit, but cannot replace a native thermal candidate.
    ///
    /// All body/gas heat, nonlinear wall feedback, endpoint reevaluation and
    /// energy ownership remain with advance_with_wall_convection. The adaptive
    /// controller is private, so a source callback cannot mint a physical report.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_with_wall_convection_checked<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: FiniteGasWindowSolve,
        wall: WallConvectionSolve,
        mut resolve: impl FnMut(
            &G::State,
            f64,
            UniformEnthalpyStepInput<'_>,
        ) -> Result<ConvectionCoefficient<S>, E>,
        mut admit: impl FnMut(UniformEnthalpyStepInput<'_>, EquilibriumPhaseState) -> Result<(), E>,
    ) -> Result<
        UniformEnthalpyStep<FiniteGasWindowReport<WallConvectionSample<S>, G>>,
        FiniteGasWindowError<WallConvectionError<E>>,
    > {
        self.advance_adaptive_with_interval(cx, input, controls, |environment, interval| {
            let source_error =
                |error| ConvectionError::Convection(WallConvectionError::Source(error));
            cx.checkpoint()
                .map_err(|_| ConvectionError::Thermal(Error::Cancelled))?;
            admit(interval, interval.initial).map_err(source_error)?;
            let next = environment.advance_with_wall_convection(
                cx,
                interval,
                wall,
                |gas, temperature| resolve(gas, temperature, interval),
            )?;
            admit(interval, next.state).map_err(source_error)?;
            cx.checkpoint()
                .map_err(|_| ConvectionError::Thermal(Error::Cancelled))?;
            Ok(next)
        })
    }

    // The only interval producers are this module's original native gas-only
    // and wall-feedback transactions. Source adapters can add refusal gates,
    // never a substitute state, heat account or integration method.
    fn advance_adaptive_with_interval<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: FiniteGasWindowSolve,
        mut advance: impl FnMut(
            &FiniteGasThermalEnvironment<'_, G>,
            UniformEnthalpyStepInput<'_>,
        ) -> Result<
            UniformEnthalpyStep<FiniteGasThermalReport<S, G>>,
            ConvectionError<E>,
        >,
    ) -> Result<UniformEnthalpyStep<FiniteGasWindowReport<S, G>>, FiniteGasWindowError<E>> {
        poll(cx)?;
        let c = controls;
        if [input.duration_s, c.maximum_step_s, c.minimum_step_s,
            c.energy_absolute_tolerance_j, c.temperature_absolute_tolerance_k,
            c.phase_fraction_absolute_tolerance, c.energy_balance_tolerance_j]
            .into_iter().any(|v| !v.is_finite() || v <= 0.0)
            || c.minimum_step_s > c.maximum_step_s || c.maximum_attempts == 0
            || c.maximum_attempts > 65_536 || c.phase_fraction_absolute_tolerance > 1.0
            || !c.relative_tolerance.is_finite() || !(0.0..1.0).contains(&c.relative_tolerance)
        { return Err(FiniteGasWindowError::Invalid("invalid finite gas/body window controls")); }
        let order = ConductionTimeOrder::First;
        let mut gas = G::clone(self.gas);
        let mut body = input.initial;
        let (mut time, mut span) = (0.0, c.maximum_step_s.min(input.duration_s));
        let (mut attempts, mut rejected) = (0, 0);
        let mut maximum_error = 0.0_f64;
        let mut samples = Vec::new();
        // Q, radiation, body source, gas source, absolute local energy defects.
        let mut totals = [Sum::default(); 5];
        let mut entropy = [Sum::default(); 13];
        while time < input.duration_s {
            poll(cx)?;
            let end = if span >= input.duration_s - time { input.duration_s } else { time + span };
            let dt = end - time;
            let middle = f64::midpoint(time, end);
            if !(time < middle && middle < end) {
                return Err(FiniteGasWindowError::Invalid("no representable thermal half interval"));
            }
            let mut interval = |gas: &G, body: EquilibriumPhaseState, start: f64, end: f64|
                -> Result<UniformEnthalpyStep<FiniteGasThermalReport<S, G>>, FiniteGasWindowError<E>> {
                if attempts == c.maximum_attempts {
                    return Err(FiniteGasWindowError::WorkLimit {
                        attempted_intervals: attempts, accepted_time_s: time,
                    });
                }
                attempts += 1;
                let body_heat = source(input.internal_heat_j, start, end, input.duration_s)?;
                let gas_heat = source(self.gas_heat_j, start, end, input.duration_s)?;
                let environment = FiniteGasThermalEnvironment {
                    gas, gas_heat_j: gas_heat, body: self.body, gas_solve: self.gas_solve,
                    exchange: self.exchange,
                    maximum_gas_temperature_increment_k: self.maximum_gas_temperature_increment_k,
                };
                advance(&environment, UniformEnthalpyStepInput {
                    initial: body, internal_heat_j: body_heat, duration_s: end - start, ..input
                }).map_err(FiniteGasWindowError::from)
            };
            let proposals = (|| {
                let coarse = interval(&gas, body, time, end)?;
                let first = interval(&gas, body, time, middle)?;
                let fine = interval(&first.report.next_gas, first.state, middle, end)?;
                Ok::<_, FiniteGasWindowError<E>>((coarse, first, fine))
            })();
            let (cause, factor) = match proposals {
                Ok((coarse, first, fine)) => {
                    let correction = order.fine_error_correction(middle - time, end - middle)
                        .ok_or(FiniteGasWindowError::Invalid("invalid time error correction"))?;
                    let normalized = compare(&gas, body, input.mass_kg, &coarse, &fine, c, correction)?;
                    if normalized <= 1.0 {
                        for (start, end, step) in [(time, middle, first), (middle, end, fine)] {
                            let r = step.report;
                            let values = [r.solid_to_gas_heat_j, r.radiation_heat_j,
                                source(input.internal_heat_j, start, end, input.duration_s)?,
                                r.external_gas_heat_j, r.solid_energy_residual_j.abs()
                                    + r.gas.residual_j.abs() + r.source_residual_j.abs()];
                            for (sum, value) in totals.iter_mut().zip(values) { sum.add(value)?; }
                            for (sum, value) in entropy.iter_mut().zip(entropy_terms(r.entropy)) {
                                sum.add(value)?;
                            }
                            body = step.state;
                            gas = r.next_gas;
                            samples.push(FiniteGasWindowSample {
                                time_s: end, body_state: body, gas_temperature_k: gas.temperature_k(),
                                convection: r.convection, solid_to_gas_heat_j: r.solid_to_gas_heat_j,
                                radiation_heat_j: r.radiation_heat_j, entropy: r.entropy,
                            });
                        }
                        time = end;
                        maximum_error = maximum_error.max(normalized);
                        let growth = order.accepted_growth_factor(normalized)
                            .ok_or(FiniteGasWindowError::Invalid("invalid accepted time error"))?;
                        span = (dt.min(c.maximum_step_s / growth) * growth)
                            .clamp(c.minimum_step_s, c.maximum_step_s);
                        continue;
                    }
                    (FiniteGasWindowError::TimeAccuracy { accepted_time_s: time,
                        span_s: end - time, normalized_error: normalized },
                        order.rejected_reduction_factor(normalized).unwrap_or(0.5))
                }
                Err(error @ FiniteGasWindowError::Step(ConvectionError::Thermal(
                    Error::OutsideBracket { .. }))) => (error, 0.5),
                Err(error @ FiniteGasWindowError::Step(ConvectionError::Thermal(
                    Error::Invalid(message)))) if message == EMPTY_HEAT_INTERVAL => (error, 0.5),
                Err(error) => return Err(error),
            };
            rejected += 1;
            let smaller = ((end - time) * factor).max(c.minimum_step_s);
            if end - time <= c.minimum_step_s || !(time < time + smaller && time + smaller < end) {
                return Err(cause);
            }
            span = smaller;
        }
        let [q, radiation, body_heat, gas_heat, absolute] = totals.map(Sum::get);
        let body_storage = finite(input.mass_kg
            * (body.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg()))?;
        let gas_storage = finite(self.gas.mass_kg()
            * self.gas.specific_internal_energy_change_to_temperature_j_kg(gas.temperature_k())
                .map_err(gas_error)?)?;
        let external_heat_j = finite(radiation - q)?;
        let body_apportionment = finite(body_heat - input.internal_heat_j)?;
        let gas_apportionment = finite(gas_heat - self.gas_heat_j)?;
        let residual = finite(body_storage + gas_storage - input.internal_heat_j - self.gas_heat_j - radiation)?;
        let body_residual = finite(body_storage - input.internal_heat_j - external_heat_j)?;
        let defect = finite(absolute + body_apportionment.abs() + gas_apportionment.abs())?;
        if residual.abs() > c.energy_balance_tolerance_j || body_residual.abs() > c.energy_balance_tolerance_j
            || defect > c.energy_balance_tolerance_j {
            return Err(FiniteGasWindowError::EnergyBalance { residual_j: residual,
                absolute_defect_j: defect, tolerance_j: c.energy_balance_tolerance_j });
        }
        poll(cx)?;
        Ok(UniformEnthalpyStep { state: body, external_heat_j,
            energy_residual_tolerance_j: c.energy_balance_tolerance_j,
            report: FiniteGasWindowReport { next_gas: gas, samples, attempted_intervals: attempts,
                rejected_comparisons: rejected, maximum_accepted_normalized_error: maximum_error,
                solid_to_gas_heat_j: q, radiation_heat_j: radiation,
                internal_body_heat_j: body_heat, external_gas_heat_j: gas_heat,
                body_heat_apportionment_residual_j: body_apportionment,
                gas_heat_apportionment_residual_j: gas_apportionment,
                body_enthalpy_change_j: body_storage, gas_internal_energy_change_j: gas_storage,
                combined_energy_residual_j: residual, sum_absolute_energy_defect_j: defect,
                entropy: entropy_report(entropy.map(Sum::get)),
            },
        })
    }
}

fn compare<S, G: IsochoricGasStorage, E>(gas: &G, body: EquilibriumPhaseState, mass: f64,
    coarse: &UniformEnthalpyStep<FiniteGasThermalReport<S, G>>,
    fine: &UniformEnthalpyStep<FiniteGasThermalReport<S, G>>,
    c: FiniteGasWindowSolve, correction: f64) -> Result<f64, FiniteGasWindowError<E>> {
    let a = coarse.state;
    let b = fine.state;
    let cg = &coarse.report.next_gas;
    let fg = &fine.report.next_gas;
    let energy = |t| gas.specific_internal_energy_change_to_temperature_j_kg(t)
        .map(|u| gas.mass_kg() * u).map_err(gas_error);
    let rows = [
        (mass * (a.specific_enthalpy_j_kg() - body.specific_enthalpy_j_kg()),
         mass * (b.specific_enthalpy_j_kg() - body.specific_enthalpy_j_kg()),
         mass * (b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg()), c.energy_absolute_tolerance_j),
        (energy(cg.temperature_k())?, energy(fg.temperature_k())?,
         cg.mass_kg() * cg.specific_internal_energy_change_to_temperature_j_kg(fg.temperature_k())
            .map_err(gas_error)?, c.energy_absolute_tolerance_j),
        (a.temperature_k() - body.temperature_k(), b.temperature_k() - body.temperature_k(),
         b.temperature_k() - a.temperature_k(), c.temperature_absolute_tolerance_k),
        (cg.temperature_k() - gas.temperature_k(), fg.temperature_k() - gas.temperature_k(),
         fg.temperature_k() - cg.temperature_k(), c.temperature_absolute_tolerance_k),
        (a.liquid_mass_fraction() - body.liquid_mass_fraction(), b.liquid_mass_fraction() - body.liquid_mass_fraction(),
         b.liquid_mass_fraction() - a.liquid_mass_fraction(), c.phase_fraction_absolute_tolerance),
    ];
    let mut normalized = 0.0_f64;
    for (a, b, difference, tolerance) in rows {
        let scale = finite(tolerance + c.relative_tolerance * a.abs().max(b.abs()))?;
        normalized = normalized.max(finite(correction * difference.abs() / scale)?);
    }
    Ok(normalized)
}

fn entropy_terms(e: FiniteGasEntropyReport) -> [f64; 13] {
    [e.storage_change_j_k[0], e.storage_change_j_k[1], e.convection_production_j_k,
     e.radiation_reservoir_input_j_k, e.radiation_transfer_production_j_k,
     e.unclassified_heat_at_endpoint_j_k[0], e.unclassified_heat_at_endpoint_j_k[1],
     e.temporal_defect_j_k[0], e.temporal_defect_j_k[1], e.source_residual_j_k,
     e.algebraic_residual_j_k[0], e.algebraic_residual_j_k[1], e.balance_residual_j_k]
}
fn entropy_report([sb, sg, convection, radiation, radiation_production, hb, hg, tb, tg, source, ab, ag, balance]: [f64; 13])
    -> FiniteGasEntropyReport {
    FiniteGasEntropyReport { storage_change_j_k: [sb, sg], convection_production_j_k: convection,
        radiation_reservoir_input_j_k: radiation, radiation_transfer_production_j_k: radiation_production,
        unclassified_heat_at_endpoint_j_k: [hb, hg], temporal_defect_j_k: [tb, tg],
        source_residual_j_k: source, algebraic_residual_j_k: [ab, ag], balance_residual_j_k: balance }
}
