//! Adaptive physical time for distinct bodies sharing one finite gas.
//!
//! Every coarse/fine proposal is the original simultaneous native transaction.
//! A component cannot hide its error behind an opposing heat transfer: enthalpy,
//! temperature and phase coordinates are checked independently for every body.
//! Gas storage is checked and counted once. Only accepted fine intervals enter
//! the retained history and whole-window energy/entropy account.

use super::*;
use crate::finite_gas_thermal::window::FiniteGasWindowSolve;
use fs_material::phase::EquilibriumPhaseState;
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

/// Local time accuracy, total native work, and retained component-history cap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasBodiesWindowSolve {
    pub time: FiniteGasWindowSolve,
    /// Counts each body at each accepted half-step, not just shared gas steps.
    /// Source payloads keep their original inner work/storage restrictions.
    pub maximum_retained_body_steps: usize,
}

/// The complete original simultaneous report at one accepted half endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodiesWindowSample<S> {
    pub time_s: f64,
    pub gas_temperature_k: f64,
    /// Original native body reports in canonical physical-ID order.
    pub bodies: Vec<UniformEnthalpyStep<FiniteGasBodyReport<S>>>,
    pub shared: FiniteGasBodiesReport,
}

/// One body's independent whole-window storage, sources and interface account.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodyWindowReport {
    pub id: u64,
    pub solid_to_gas_heat_j: f64,
    pub radiation_heat_j: f64,
    pub internal_heat_j: f64,
    pub heat_apportionment_residual_j: f64,
    pub enthalpy_change_j: f64,
    pub energy_residual_j: f64,
    pub entropy: FiniteGasBodyEntropyReport,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodiesWindowReport<S> {
    pub samples: Vec<FiniteGasBodiesWindowSample<S>>,
    pub attempted_intervals: usize,
    pub rejected_comparisons: usize,
    pub maximum_accepted_normalized_error: f64,
    pub heat_into_gas_j: f64,
    pub heat_from_bodies_j: f64,
    pub source_residual_j: f64,
    pub radiation_heat_j: f64,
    pub internal_body_heat_j: f64,
    pub external_gas_heat_j: f64,
    pub gas_heat_apportionment_residual_j: f64,
    pub body_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub combined_energy_residual_j: f64,
    pub sum_absolute_energy_defect_j: f64,
    /// Sum of independently recomputed final body, gas and matching defects.
    pub endpoint_absolute_energy_defect_j: f64,
    pub entropy: FiniteGasEntropyReport,
}

/// All physical endpoints remain a single immutable proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodiesWindowCandidate<S, G> {
    gas: G,
    bodies: Vec<UniformEnthalpyStep<FiniteGasBodyWindowReport>>,
    report: FiniteGasBodiesWindowReport<S>,
}
impl<S, G> FiniteGasBodiesWindowCandidate<S, G> {
    pub const fn gas(&self) -> &G {
        &self.gas
    }
    pub fn bodies(&self) -> &[UniformEnthalpyStep<FiniteGasBodyWindowReport>] {
        &self.bodies
    }
    pub const fn report(&self) -> &FiniteGasBodiesWindowReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        G,
        Vec<UniformEnthalpyStep<FiniteGasBodyWindowReport>>,
        FiniteGasBodiesWindowReport<S>,
    ) {
        (self.gas, self.bodies, self.report)
    }
}

#[derive(Debug)]
pub enum FiniteGasBodiesWindowError<E> {
    Invalid(&'static str),
    Step(FiniteGasBodiesError<E>),
    Cancelled,
    WorkLimit {
        attempted_intervals: usize,
        accepted_time_s: f64,
    },
    HistoryLimit {
        retained_body_steps: usize,
        maximum_body_steps: usize,
    },
    TimeAccuracy {
        accepted_time_s: f64,
        span_s: f64,
        normalized_error: f64,
    },
    EnergyBalance {
        residual_j: f64,
        absolute_defect_j: f64,
        tolerance_j: f64,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteGasBodiesWindowError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shared gas/body thermal window: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FiniteGasBodiesWindowError<E> {}
impl<E> From<FiniteGasBodiesError<E>> for FiniteGasBodiesWindowError<E> {
    fn from(error: FiniteGasBodiesError<E>) -> Self {
        match error {
            FiniteGasBodiesError::Cancelled => Self::Cancelled,
            error => Self::Step(error),
        }
    }
}
impl<E> From<Error> for FiniteGasBodiesWindowError<E> {
    fn from(error: Error) -> Self {
        FiniteGasBodiesError::from(error).into()
    }
}

fn finite<E>(value: f64) -> Result<f64, FiniteGasBodiesWindowError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FiniteGasBodiesWindowError::Invalid(
            "unrepresentable shared thermal window arithmetic",
        ))
    }
}
fn poll<E>(cx: &Cx<'_>) -> Result<(), FiniteGasBodiesWindowError<E>> {
    cx.checkpoint()
        .map_err(|_| FiniteGasBodiesWindowError::Cancelled)
}
fn source<E>(
    total: f64,
    start: f64,
    end: f64,
    duration: f64,
) -> Result<f64, FiniteGasBodiesWindowError<E>> {
    let at = |time| {
        if time == duration {
            total
        } else {
            total * (time / duration)
        }
    };
    let value = finite(at(end) - at(start))?;
    if total != 0.0 && value == 0.0 {
        return Err(FiniteGasBodiesWindowError::Invalid(
            "nonzero interval heat is unrepresentable",
        ));
    }
    Ok(value)
}
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), FiniteGasBodiesWindowError<E>> {
        let (next, residual) = fs_math::eft::two_sum(self.value, value);
        self.value = next;
        self.correction += residual;
        finite(self.get()).map(|_| ())
    }
    fn get(self) -> f64 {
        self.value + self.correction
    }
}
#[derive(Clone, Default)]
struct BodyAccount {
    heat: [Sum; 3],
    entropy: [Sum; 7],
}

impl<G: IsochoricGasStorage> FiniteGasBodiesEnvironment<'_, G> {
    pub fn advance_adaptive(
        &self,
        cx: &Cx<'_>,
        bodies: &[FiniteGasBody<'_>],
        controls: FiniteGasBodiesWindowSolve,
    ) -> Result<
        FiniteGasBodiesWindowCandidate<(), G>,
        FiniteGasBodiesWindowError<core::convert::Infallible>,
    > {
        self.advance_adaptive_with_convection(cx, bodies, controls, |body, _| {
            Ok(ConvectionCoefficient {
                coefficient_w_per_m2_k: body.environment.convection_w_per_m2_k,
                source: (),
            })
        })
    }

    /// Advance all original bodies and the ONE original gas together.
    ///
    /// Authored body and gas heats are constant-power total window energies.
    /// The actual interval descriptor retains each body's curve, geometry,
    /// mass, original environment, current state and apportioned source heat.
    /// No component is updated sequentially through the shared gas.
    ///
    /// Only original outer heat-bracket/increment refusals can shorten a
    /// proposal. Body/material/source/nonlinear errors retain their identity
    /// and propagate. Callbacks must be repeatable: external effects cannot
    /// be rolled back across coarse/fine/rejected proposals.
    pub fn advance_adaptive_with_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        bodies: &[FiniteGasBody<'_>],
        controls: FiniteGasBodiesWindowSolve,
        mut resolve: impl FnMut(&FiniteGasBody<'_>, &G::State) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<FiniteGasBodiesWindowCandidate<S, G>, FiniteGasBodiesWindowError<E>> {
        self.advance_adaptive_native(cx, bodies, controls, |environment, interval| {
            environment.advance_with_convection(cx, interval, &mut resolve)
        })
    }

    pub fn advance_adaptive_with_wall_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        bodies: &[FiniteGasBody<'_>],
        controls: FiniteGasBodiesWindowSolve,
        wall: WallConvectionSolve,
        mut resolve: impl FnMut(
            &FiniteGasBody<'_>,
            &G::State,
            f64,
        ) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<
        FiniteGasBodiesWindowCandidate<WallConvectionSample<S>, G>,
        FiniteGasBodiesWindowError<WallConvectionError<E>>,
    > {
        self.advance_adaptive_with_wall_convection_checked(
            cx,
            bodies,
            controls,
            wall,
            |body, gas, temperature| resolve(body, gas, temperature),
            |_, _| Ok(()),
        )
    }

    /// Source-conjugate integration with each original body's actual interval.
    ///
    /// Initial and final admission visits EVERY body of every coarse, fine
    /// and subsequently rejected native proposal. The caller can refuse an
    /// unsupported caloric/phase/conductivity state but cannot supply a new
    /// thermal candidate. Wall feedback retains both original source samples.
    /// Native per-body wall residuals share ONE original source-work allowance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_with_wall_convection_checked<S, E>(
        &self,
        cx: &Cx<'_>,
        bodies: &[FiniteGasBody<'_>],
        controls: FiniteGasBodiesWindowSolve,
        wall: WallConvectionSolve,
        mut resolve: impl FnMut(
            &FiniteGasBody<'_>,
            &G::State,
            f64,
        ) -> Result<ConvectionCoefficient<S>, E>,
        mut admit: impl FnMut(&FiniteGasBody<'_>, EquilibriumPhaseState) -> Result<(), E>,
    ) -> Result<
        FiniteGasBodiesWindowCandidate<WallConvectionSample<S>, G>,
        FiniteGasBodiesWindowError<WallConvectionError<E>>,
    > {
        self.advance_adaptive_native(cx, bodies, controls, |environment, interval| {
            let source_error = |id, error| FiniteGasBodiesError::Body {
                id,
                error: ConvectionError::Convection(WallConvectionError::Source(error)),
            };
            for body in interval {
                cx.checkpoint()
                    .map_err(|_| FiniteGasBodiesError::Cancelled)?;
                admit(body, body.input.initial).map_err(|error| source_error(body.id, error))?;
            }
            let next =
                environment.advance_with_wall_convection(cx, interval, wall, &mut resolve)?;
            for (body, endpoint) in interval.iter().zip(next.bodies()) {
                cx.checkpoint()
                    .map_err(|_| FiniteGasBodiesError::Cancelled)?;
                admit(body, endpoint.state).map_err(|error| source_error(body.id, error))?;
            }
            cx.checkpoint()
                .map_err(|_| FiniteGasBodiesError::Cancelled)?;
            Ok(next)
        })
    }

    // Both interval producers are this module's original simultaneous owners.
    // This private composition cannot accept an externally minted candidate.
    fn advance_adaptive_native<S, E>(
        &self,
        cx: &Cx<'_>,
        bodies: &[FiniteGasBody<'_>],
        controls: FiniteGasBodiesWindowSolve,
        mut advance: impl FnMut(
            &FiniteGasBodiesEnvironment<'_, G>,
            &[FiniteGasBody<'_>],
        )
            -> Result<FiniteGasBodiesCandidate<S, G>, FiniteGasBodiesError<E>>,
    ) -> Result<FiniteGasBodiesWindowCandidate<S, G>, FiniteGasBodiesWindowError<E>> {
        poll(cx)?;
        let c = controls.time;
        if self.maximum_bodies == 0
            || self.maximum_bodies > MAX_BODIES
            || bodies.is_empty()
            || bodies.len() > self.maximum_bodies
        {
            return Err(FiniteGasBodiesWindowError::Invalid(
                "bounded nonempty shared body set required",
            ));
        }
        let duration = bodies[0].input.duration_s;
        if [
            duration,
            c.maximum_step_s,
            c.minimum_step_s,
            c.energy_absolute_tolerance_j,
            c.temperature_absolute_tolerance_k,
            c.phase_fraction_absolute_tolerance,
            c.energy_balance_tolerance_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
            || c.minimum_step_s > c.maximum_step_s
            || c.maximum_attempts == 0
            || c.maximum_attempts > 65_536
            || c.phase_fraction_absolute_tolerance > 1.0
            || !c.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&c.relative_tolerance)
            || controls.maximum_retained_body_steps < 2 * bodies.len()
            || bodies.iter().any(|body| {
                body.input.duration_s != duration || !body.input.internal_heat_j.is_finite()
            })
            || !self.gas_heat_j.is_finite()
        {
            return Err(FiniteGasBodiesWindowError::Invalid(
                "invalid shared gas/body window controls or sources",
            ));
        }
        let mut ordered = bodies.to_vec();
        ordered.sort_unstable_by_key(|body| body.id);
        if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(FiniteGasBodiesWindowError::Invalid(
                "unique physical body IDs required",
            ));
        }
        let order = ConductionTimeOrder::First;
        let mut gas = G::clone(self.gas);
        let mut states: Vec<_> = ordered.iter().map(|body| body.input.initial).collect();
        let mut accounts = vec![BodyAccount::default(); ordered.len()];
        let (mut time, mut span) = (0.0, c.maximum_step_s.min(duration));
        let (mut attempts, mut rejected, mut retained) = (0, 0, 0);
        let mut maximum_error = 0.0_f64;
        let mut samples = Vec::new();
        // Native total-Q, measured body Q, gas heat, absolute native defects.
        let mut totals = [Sum::default(); 4];
        let mut entropy = [Sum::default(); 13];
        while time < duration {
            poll(cx)?;
            let added = 2 * ordered.len();
            if retained > controls.maximum_retained_body_steps - added {
                return Err(FiniteGasBodiesWindowError::HistoryLimit {
                    retained_body_steps: retained,
                    maximum_body_steps: controls.maximum_retained_body_steps,
                });
            }
            let end = if span >= duration - time {
                duration
            } else {
                time + span
            };
            let dt = end - time;
            let middle = f64::midpoint(time, end);
            if !(time < middle && middle < end) {
                return Err(FiniteGasBodiesWindowError::Invalid(
                    "no representable shared thermal half interval",
                ));
            }
            let mut interval = |gas: &G,
                                states: &[EquilibriumPhaseState],
                                start: f64,
                                end: f64|
             -> Result<
                FiniteGasBodiesCandidate<S, G>,
                FiniteGasBodiesWindowError<E>,
            > {
                if attempts == c.maximum_attempts {
                    return Err(FiniteGasBodiesWindowError::WorkLimit {
                        attempted_intervals: attempts,
                        accepted_time_s: time,
                    });
                }
                attempts += 1;
                let mut inputs = Vec::with_capacity(ordered.len());
                for (body, initial) in ordered.iter().zip(states) {
                    inputs.push(FiniteGasBody {
                        id: body.id,
                        environment: body.environment,
                        input: UniformEnthalpyStepInput {
                            initial: *initial,
                            duration_s: end - start,
                            internal_heat_j: source(
                                body.input.internal_heat_j,
                                start,
                                end,
                                duration,
                            )?,
                            ..body.input
                        },
                    });
                }
                let environment = FiniteGasBodiesEnvironment {
                    gas,
                    gas_heat_j: source(self.gas_heat_j, start, end, duration)?,
                    gas_solve: self.gas_solve,
                    exchange: self.exchange,
                    maximum_gas_temperature_increment_k: self.maximum_gas_temperature_increment_k,
                    maximum_bodies: self.maximum_bodies,
                };
                advance(&environment, &inputs).map_err(Into::into)
            };
            let proposals = (|| {
                let coarse = interval(&gas, &states, time, end)?;
                let first = interval(&gas, &states, time, middle)?;
                let middle_states: Vec<_> = first.bodies().iter().map(|body| body.state).collect();
                let fine = interval(first.gas(), &middle_states, middle, end)?;
                Ok::<_, FiniteGasBodiesWindowError<E>>((coarse, first, fine))
            })();
            let (cause, factor) = match proposals {
                Ok((coarse, first, fine)) => {
                    let correction = order
                        .fine_error_correction(middle - time, end - middle)
                        .ok_or(FiniteGasBodiesWindowError::Invalid(
                            "invalid shared time error correction",
                        ))?;
                    let error = compare(&gas, &ordered, &states, &coarse, &fine, c, correction)?;
                    if error <= 1.0 {
                        for (start, stop, step) in [(time, middle, first), (middle, end, fine)] {
                            let FiniteGasBodiesCandidate {
                                gas: next_gas,
                                bodies: next_bodies,
                                report,
                            } = step;
                            for (((original, state), account), next) in ordered
                                .iter()
                                .zip(&mut states)
                                .zip(&mut accounts)
                                .zip(&next_bodies)
                            {
                                let values = [
                                    next.report.solid_to_gas_heat_j,
                                    next.report.radiation_heat_j,
                                    source(original.input.internal_heat_j, start, stop, duration)?,
                                ];
                                for (sum, value) in account.heat.iter_mut().zip(values) {
                                    sum.add(value)?;
                                }
                                for (sum, value) in account
                                    .entropy
                                    .iter_mut()
                                    .zip(body_entropy_terms(next.report.entropy))
                                {
                                    sum.add(value)?;
                                }
                                *state = next.state;
                            }
                            for (sum, value) in totals.iter_mut().zip([
                                report.heat_into_gas_j,
                                report.heat_from_bodies_j,
                                report.external_gas_heat_j,
                                report.absolute_energy_defect_j,
                            ]) {
                                sum.add(value)?;
                            }
                            for (sum, value) in
                                entropy.iter_mut().zip(entropy_terms(report.entropy))
                            {
                                sum.add(value)?;
                            }
                            gas = next_gas;
                            samples.push(FiniteGasBodiesWindowSample {
                                time_s: stop,
                                gas_temperature_k: gas.temperature_k(),
                                bodies: next_bodies,
                                shared: report,
                            });
                        }
                        retained += added;
                        time = end;
                        maximum_error = maximum_error.max(error);
                        let growth = order.accepted_growth_factor(error).ok_or(
                            FiniteGasBodiesWindowError::Invalid(
                                "invalid accepted shared time error",
                            ),
                        )?;
                        span = (dt.min(c.maximum_step_s / growth) * growth)
                            .clamp(c.minimum_step_s, c.maximum_step_s);
                        continue;
                    }
                    (
                        FiniteGasBodiesWindowError::TimeAccuracy {
                            accepted_time_s: time,
                            span_s: dt,
                            normalized_error: error,
                        },
                        order.rejected_reduction_factor(error).unwrap_or(0.5),
                    )
                }
                Err(
                    error @ FiniteGasBodiesWindowError::Step(FiniteGasBodiesError::Thermal(
                        Error::OutsideBracket { .. },
                    )),
                ) => (error, 0.5),
                Err(
                    error @ FiniteGasBodiesWindowError::Step(FiniteGasBodiesError::Thermal(
                        Error::Invalid(message),
                    )),
                ) if message == EMPTY_HEAT_INTERVAL => (error, 0.5),
                Err(error) => return Err(error),
            };
            rejected += 1;
            let smaller = (dt * factor).max(c.minimum_step_s);
            if dt <= c.minimum_step_s || !(time < time + smaller && time + smaller < end) {
                return Err(cause);
            }
            span = smaller;
        }
        let [q, measured_q, gas_heat, absolute] = totals.map(Sum::get);
        let mut body_storage = Sum::default();
        let mut internal_heat = Sum::default();
        let mut authored_heat = Sum::default();
        let mut radiation = Sum::default();
        let mut apportionment_defect = Sum::default();
        let mut endpoint_defect = Sum::default();
        let mut endpoints = Vec::with_capacity(ordered.len());
        for ((original, state), account) in ordered.iter().zip(states).zip(accounts) {
            let [body_q, body_radiation, body_heat] = account.heat.map(Sum::get);
            let storage = finite(
                original.input.mass_kg
                    * (state.specific_enthalpy_j_kg()
                        - original.input.initial.specific_enthalpy_j_kg()),
            )?;
            let external = finite(body_radiation - body_q)?;
            let apportionment = finite(body_heat - original.input.internal_heat_j)?;
            let residual = finite(storage - original.input.internal_heat_j - external)?;
            if residual.abs() > c.energy_balance_tolerance_j {
                return Err(FiniteGasBodiesWindowError::EnergyBalance {
                    residual_j: residual,
                    absolute_defect_j: finite(absolute + apportionment.abs())?,
                    tolerance_j: c.energy_balance_tolerance_j,
                });
            }
            endpoint_defect.add(residual.abs())?;
            body_storage.add(storage)?;
            internal_heat.add(body_heat)?;
            authored_heat.add(original.input.internal_heat_j)?;
            radiation.add(body_radiation)?;
            apportionment_defect.add(apportionment.abs())?;
            endpoints.push(UniformEnthalpyStep {
                state,
                external_heat_j: external,
                energy_residual_tolerance_j: c.energy_balance_tolerance_j,
                report: FiniteGasBodyWindowReport {
                    id: original.id,
                    solid_to_gas_heat_j: body_q,
                    radiation_heat_j: body_radiation,
                    internal_heat_j: body_heat,
                    heat_apportionment_residual_j: apportionment,
                    enthalpy_change_j: storage,
                    energy_residual_j: residual,
                    entropy: body_entropy_report(account.entropy.map(Sum::get)),
                },
            });
        }
        let gas_storage = finite(
            self.gas.mass_kg()
                * self
                    .gas
                    .specific_internal_energy_change_to_temperature_j_kg(gas.temperature_k())
                    .map_err(gas_error)?,
        )?;
        let gas_apportionment = finite(gas_heat - self.gas_heat_j)?;
        let source_residual = finite(q - measured_q)?;
        endpoint_defect.add(finite(gas_storage - self.gas_heat_j - q)?.abs())?;
        endpoint_defect.add(source_residual.abs())?;
        let residual = finite(
            body_storage.get() + gas_storage
                - authored_heat.get()
                - self.gas_heat_j
                - radiation.get(),
        )?;
        let defect = finite(absolute + apportionment_defect.get() + gas_apportionment.abs())?;
        if residual.abs() > c.energy_balance_tolerance_j
            || defect > c.energy_balance_tolerance_j
            || endpoint_defect.get() > c.energy_balance_tolerance_j
        {
            return Err(FiniteGasBodiesWindowError::EnergyBalance {
                residual_j: residual,
                absolute_defect_j: defect.max(endpoint_defect.get()),
                tolerance_j: c.energy_balance_tolerance_j,
            });
        }
        poll(cx)?;
        Ok(FiniteGasBodiesWindowCandidate {
            gas,
            bodies: endpoints,
            report: FiniteGasBodiesWindowReport {
                samples,
                attempted_intervals: attempts,
                rejected_comparisons: rejected,
                maximum_accepted_normalized_error: maximum_error,
                heat_into_gas_j: q,
                heat_from_bodies_j: measured_q,
                source_residual_j: source_residual,
                radiation_heat_j: radiation.get(),
                internal_body_heat_j: internal_heat.get(),
                external_gas_heat_j: gas_heat,
                gas_heat_apportionment_residual_j: gas_apportionment,
                body_enthalpy_change_j: body_storage.get(),
                gas_internal_energy_change_j: gas_storage,
                combined_energy_residual_j: residual,
                sum_absolute_energy_defect_j: defect,
                endpoint_absolute_energy_defect_j: endpoint_defect.get(),
                entropy: entropy_report(entropy.map(Sum::get)),
            },
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn compare<S, G: IsochoricGasStorage, E>(
    gas: &G,
    bodies: &[FiniteGasBody<'_>],
    states: &[EquilibriumPhaseState],
    coarse: &FiniteGasBodiesCandidate<S, G>,
    fine: &FiniteGasBodiesCandidate<S, G>,
    c: FiniteGasWindowSolve,
    correction: f64,
) -> Result<f64, FiniteGasBodiesWindowError<E>> {
    let mut maximum = 0.0_f64;
    let mut coordinate = |a: f64,
                          b: f64,
                          difference: f64,
                          absolute: f64|
     -> Result<(), FiniteGasBodiesWindowError<E>> {
        let scale = finite(absolute + c.relative_tolerance * a.abs().max(b.abs()))?;
        maximum = maximum.max(finite(correction * difference.abs() / scale)?);
        Ok(())
    };
    for (((original, initial), a), b) in bodies
        .iter()
        .zip(states)
        .zip(coarse.bodies())
        .zip(fine.bodies())
    {
        let mass = original.input.mass_kg;
        coordinate(
            mass * (a.state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg()),
            mass * (b.state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg()),
            mass * (b.state.specific_enthalpy_j_kg() - a.state.specific_enthalpy_j_kg()),
            c.energy_absolute_tolerance_j,
        )?;
        coordinate(
            a.state.temperature_k() - initial.temperature_k(),
            b.state.temperature_k() - initial.temperature_k(),
            b.state.temperature_k() - a.state.temperature_k(),
            c.temperature_absolute_tolerance_k,
        )?;
        coordinate(
            a.state.liquid_mass_fraction() - initial.liquid_mass_fraction(),
            b.state.liquid_mass_fraction() - initial.liquid_mass_fraction(),
            b.state.liquid_mass_fraction() - a.state.liquid_mass_fraction(),
            c.phase_fraction_absolute_tolerance,
        )?;
    }
    let cg = coarse.gas();
    let fg = fine.gas();
    let energy = |t| {
        gas.specific_internal_energy_change_to_temperature_j_kg(t)
            .map(|u| gas.mass_kg() * u)
            .map_err(gas_error)
    };
    coordinate(
        energy(cg.temperature_k())?,
        energy(fg.temperature_k())?,
        cg.mass_kg()
            * cg.specific_internal_energy_change_to_temperature_j_kg(fg.temperature_k())
                .map_err(gas_error)?,
        c.energy_absolute_tolerance_j,
    )?;
    coordinate(
        cg.temperature_k() - gas.temperature_k(),
        fg.temperature_k() - gas.temperature_k(),
        fg.temperature_k() - cg.temperature_k(),
        c.temperature_absolute_tolerance_k,
    )?;
    Ok(maximum)
}
fn body_entropy_terms(e: FiniteGasBodyEntropyReport) -> [f64; 7] {
    [
        e.storage_change_j_k,
        e.convection_production_j_k,
        e.radiation_reservoir_input_j_k,
        e.radiation_transfer_production_j_k,
        e.unclassified_heat_at_endpoint_j_k,
        e.temporal_defect_j_k,
        e.algebraic_residual_j_k,
    ]
}
fn body_entropy_report(
    [
        storage,
        convection,
        radiation,
        radiation_production,
        heat,
        time,
        algebraic,
    ]: [f64; 7],
) -> FiniteGasBodyEntropyReport {
    FiniteGasBodyEntropyReport {
        storage_change_j_k: storage,
        convection_production_j_k: convection,
        radiation_reservoir_input_j_k: radiation,
        radiation_transfer_production_j_k: radiation_production,
        unclassified_heat_at_endpoint_j_k: heat,
        temporal_defect_j_k: time,
        algebraic_residual_j_k: algebraic,
    }
}
fn entropy_terms(e: FiniteGasEntropyReport) -> [f64; 13] {
    [
        e.storage_change_j_k[0],
        e.storage_change_j_k[1],
        e.convection_production_j_k,
        e.radiation_reservoir_input_j_k,
        e.radiation_transfer_production_j_k,
        e.unclassified_heat_at_endpoint_j_k[0],
        e.unclassified_heat_at_endpoint_j_k[1],
        e.temporal_defect_j_k[0],
        e.temporal_defect_j_k[1],
        e.source_residual_j_k,
        e.algebraic_residual_j_k[0],
        e.algebraic_residual_j_k[1],
        e.balance_residual_j_k,
    ]
}
fn entropy_report(
    [
        sb,
        sg,
        convection,
        radiation,
        radiation_production,
        hb,
        hg,
        tb,
        tg,
        source,
        ab,
        ag,
        balance,
    ]: [f64; 13],
) -> FiniteGasEntropyReport {
    FiniteGasEntropyReport {
        storage_change_j_k: [sb, sg],
        convection_production_j_k: convection,
        radiation_reservoir_input_j_k: radiation,
        radiation_transfer_production_j_k: radiation_production,
        unclassified_heat_at_endpoint_j_k: [hb, hg],
        temporal_defect_j_k: [tb, tg],
        source_residual_j_k: source,
        algebraic_residual_j_k: [ab, ag],
        balance_residual_j_k: balance,
    }
}
