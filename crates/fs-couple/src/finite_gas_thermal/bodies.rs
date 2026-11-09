//! One shared gas endpoint for independent native enthalpy/phase bodies.
//!
//! Each total-Q trial heats the original gas once and solves every original
//! body against that same endpoint. Q minus the SUM of actual body heat flows
//! is the only outer unknown. Individual transfers are never redistributed,
//! even when opposite transfers cancel in the shared gas. This fixed-volume
//! backward-Euler source operator introduces no circuit, piston, spatial gas
//! field, direct body/body contact or alternate material law.

use super::*;

/// Explicit work bound on independent native body solves per scalar trial.
pub const MAX_BODIES: usize = 128;

/// One distinct physical body with its own caloric curve, geometry and sources.
/// IDs designate bodies, not material or boundary objects; distinct specimens
/// may share the same immutable source curve/environment. Inputs share one dt.
#[derive(Clone, Copy, Debug)]
pub struct FiniteGasBody<'a> {
    pub id: u64,
    pub environment: &'a LumpedThermalEnvironment,
    pub input: UniformEnthalpyStepInput<'a>,
}

#[derive(Clone, Debug)]
pub struct FiniteGasBodiesEnvironment<'a, G = IsochoricGasParcel> {
    pub gas: &'a G,
    /// Total prescribed gas heat over the shared interval, excluding all contacts.
    pub gas_heat_j: f64,
    pub gas_solve: GasEnergySolve,
    /// TOTAL signed gas exchange bracket and all-owner energy allowance.
    pub exchange: HeatExchangeSolve,
    pub maximum_gas_temperature_increment_k: f64,
    pub maximum_bodies: usize,
}

/// Body storage and this body's interface contributions, in J/K.
/// Gas storage and its sources are intentionally absent; they are counted once
/// in the shared report. Transfer production belongs to the body/gas interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasBodyEntropyReport {
    pub storage_change_j_k: f64,
    pub convection_production_j_k: f64,
    pub radiation_reservoir_input_j_k: f64,
    pub radiation_transfer_production_j_k: f64,
    pub unclassified_heat_at_endpoint_j_k: f64,
    pub temporal_defect_j_k: f64,
    pub algebraic_residual_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodyReport<S> {
    pub id: u64,
    pub body: LumpedEnthalpyMarch,
    pub convection: ConvectionCoefficient<S>,
    /// Minus dt times the original body's endpoint convective power [J].
    pub solid_to_gas_heat_j: f64,
    pub radiation_heat_j: f64,
    pub enthalpy_change_j: f64,
    pub energy_residual_j: f64,
    pub entropy: FiniteGasBodyEntropyReport,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodiesReport {
    pub gas: GasHeatStep,
    /// Total Q used by the original gas inversion [J].
    pub heat_into_gas_j: f64,
    /// Sum of individual measured body-to-gas transfers [J].
    pub heat_from_bodies_j: f64,
    /// Q minus the measured sum; never reassigned to a physical body or source.
    pub source_residual_j: f64,
    pub external_gas_heat_j: f64,
    pub internal_body_heat_j: f64,
    pub radiation_heat_j: f64,
    pub body_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub combined_energy_residual_j: f64,
    /// Sum |body defects| + |gas defect| + |source mismatch|. One common
    /// allowance bounds all bodies together, not one renewed allowance per body.
    pub absolute_energy_defect_j: f64,
    /// Body slot aggregates distinct bodies at their own temperatures; gas once.
    pub entropy: FiniteGasEntropyReport,
    pub exchange_trials: usize,
    /// Final native body proposals for all gas trials. Additional wall-feedback
    /// proposals, when selected, are retained in each body's convection sample.
    pub body_trials: usize,
}

/// Every endpoint is proposed together. Neither original gas nor body mutates.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasBodiesCandidate<S, G> {
    gas: G,
    bodies: Vec<UniformEnthalpyStep<FiniteGasBodyReport<S>>>,
    report: FiniteGasBodiesReport,
}
impl<S, G> FiniteGasBodiesCandidate<S, G> {
    pub const fn gas(&self) -> &G { &self.gas }
    /// Canonical physical body-ID order, independent of input slice order.
    pub fn bodies(&self) -> &[UniformEnthalpyStep<FiniteGasBodyReport<S>>] { &self.bodies }
    pub const fn report(&self) -> &FiniteGasBodiesReport { &self.report }
    pub fn into_parts(self) -> (G, Vec<UniformEnthalpyStep<FiniteGasBodyReport<S>>>, FiniteGasBodiesReport) {
        (self.gas, self.bodies, self.report)
    }
}

#[derive(Debug)]
pub enum FiniteGasBodiesError<E> {
    Thermal(Error),
    Body { id: u64, error: ConvectionError<E> },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteGasBodiesError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shared gas/enthalpy bodies: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FiniteGasBodiesError<E> {}
impl<E> From<Error> for FiniteGasBodiesError<E> {
    fn from(error: Error) -> Self {
        match error { Error::Cancelled => Self::Cancelled, error => Self::Thermal(error) }
    }
}
fn body_error<E>(id: u64, error: Error) -> FiniteGasBodiesError<E> {
    match error {
        Error::Cancelled | Error::Interface(ConductionError::Cancelled { .. }) => FiniteGasBodiesError::Cancelled,
        error => FiniteGasBodiesError::Body { id, error: ConvectionError::Thermal(error) },
    }
}
fn scalar_error<E>(error: feedback::ScalarError<FiniteGasBodiesError<E>>) -> FiniteGasBodiesError<E> {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Cancelled => FiniteGasBodiesError::Cancelled,
        feedback::ScalarError::Invalid(message) => Error::Invalid(message).into(),
        feedback::ScalarError::OutsideBracket { lower, upper } => Error::OutsideBracket {
            lower_residual_j: lower, upper_residual_j: upper,
        }.into(),
        feedback::ScalarError::NotConverged { iterations, residual } => Error::NotConverged {
            iterations, residual_j: residual,
        }.into(),
    }
}

struct BodyTrial<S> {
    body: UniformEnthalpyStep<LumpedEnthalpyMarch>,
    convection: ConvectionCoefficient<S>,
    heat_j: f64,
}
struct Trial<S, G> {
    next_gas: G,
    gas: GasHeatStep,
    bodies: Vec<BodyTrial<S>>,
    heat_j: f64,
    residual_j: f64,
    trials: usize,
}

impl<G: IsochoricGasStorage> FiniteGasBodiesEnvironment<'_, G> {
    /// Use each original boundary's declared whole-surface convection coefficient.
    pub fn advance(&self, cx: &Cx<'_>, bodies: &[FiniteGasBody<'_>])
        -> Result<FiniteGasBodiesCandidate<(), G>, FiniteGasBodiesError<core::convert::Infallible>> {
        self.advance_with_convection(cx, bodies, |body, _| Ok(ConvectionCoefficient {
            coefficient_w_per_m2_k: body.environment.convection_w_per_m2_k, source: (),
        }))
    }

    /// Resolve every body's convection at its OWN wall endpoint and the SAME
    /// finite gas endpoint. Film-temperature and natural-convection laws use
    /// the original single-body wall feedback, without sequential gas updates.
    ///
    /// `resolve` receives the original body descriptor, the common candidate
    /// gas and that body's proposed wall temperature [K]. It may query the
    /// original gas source at a declared film temperature without changing the
    /// physical bulk inventory. Source identities and both final wall-law
    /// evaluations survive in each body's report. All callbacks must be
    /// repeatable and leave external histories unchanged.
    ///
    /// Wall temperature controls apply independently to every original body.
    /// The shared source-work tolerance is divided equally between bodies so
    /// the sum of absolute wall-law mismatches stays within ONE allowance.
    /// Only this numerical allowance is partitioned; no physical heat or
    /// conductance is redistributed. Existing gas storage, body phase/Biot,
    /// external radiation, energy/entropy and atomicity checks remain active.
    ///
    /// Each wall root takes at most `wall.maximum_iterations + 2` native body
    /// proposals per gas trial, plus the original final body proposal. Arbitrary
    /// material laws do not acquire a convergence or uniqueness guarantee.
    pub fn advance_with_wall_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        bodies: &[FiniteGasBody<'_>],
        wall: WallConvectionSolve,
        mut resolve: impl FnMut(
            &FiniteGasBody<'_>,
            &G::State,
            f64,
        ) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<
        FiniteGasBodiesCandidate<WallConvectionSample<S>, G>,
        FiniteGasBodiesError<WallConvectionError<E>>,
    > {
        // The shared owner admits the bounded body set before invoking this
        // callback. Allocate a counter only for each admitted physical body.
        let mut wall_trials = std::collections::BTreeMap::new();
        let source_tolerance_j = self.exchange.source_tolerance_j / bodies.len().max(1) as f64;
        let result = self.advance_with_convection(cx, bodies, |body, gas| {
            super::wall::resolve_wall_convection(
                cx,
                body.environment,
                body.input,
                gas,
                wall,
                source_tolerance_j,
                wall_trials.entry(body.id).or_insert(0),
                |gas, temperature_k| resolve(body, gas, temperature_k),
            )
        });
        match result {
            Err(FiniteGasBodiesError::Body {
                error: ConvectionError::Convection(WallConvectionError::Cancelled),
                ..
            }) => Err(FiniteGasBodiesError::Cancelled),
            other => other,
        }
    }

    /// Resolve distinct whole-body convection laws at one shared gas endpoint.
    ///
    /// Every trial starts from the same accepted inputs. Each native body owns
    /// its h/phase inverse, internal heat, radiation enclosure, transport and
    /// Biot admission. Opposite body transfers remain physical even at zero net
    /// gas heat. The gas caloric/increment bound is shared, including real-gas
    /// storage through IsochoricGasStorage. No individual qi is clipped or
    /// corrected to force the gas sum to close.
    ///
    /// The callback sees a body's original descriptor and the common proposed
    /// gas state, not an independently prescribed wall endpoint. It must be
    /// repeatable; external callback effects cannot be rolled back. Native
    /// source refusals and cancellation propagate with body identity. One
    /// complete trial performs one gas inversion and exactly N body solves.
    pub fn advance_with_convection<S, E>(
        &self, cx: &Cx<'_>, bodies: &[FiniteGasBody<'_>],
        mut resolve: impl FnMut(&FiniteGasBody<'_>, &G::State) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<FiniteGasBodiesCandidate<S, G>, FiniteGasBodiesError<E>> {
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        self.gas_solve.validate().map_err(gas_error)?;
        feedback::admit::<ConductionError>(self.exchange)?;
        if self.maximum_bodies == 0 || self.maximum_bodies > MAX_BODIES
            || bodies.is_empty() || bodies.len() > self.maximum_bodies
            || !self.gas_heat_j.is_finite() || !self.maximum_gas_temperature_increment_k.is_finite()
            || self.maximum_gas_temperature_increment_k <= 0.0 {
            return Err(Error::Invalid("bounded nonempty body set and finite gas sources/increment required").into());
        }
        let mut ordered: Vec<_> = bodies.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        let duration = ordered[0].input.duration_s;
        if !duration.is_finite() || duration <= 0.0
            || ordered.iter().any(|body| body.input.duration_s != duration)
            || ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(Error::Invalid("unique physical body IDs and one finite positive interval required").into());
        }
        // Reuse the established GAS-ONLY bracket helper; this descriptor's
        // body is not advanced separately or given a privileged thermal role.
        let bounds = FiniteGasThermalEnvironment { gas: self.gas, body: ordered[0].environment,
            gas_heat_j: self.gas_heat_j, gas_solve: self.gas_solve, exchange: self.exchange,
            maximum_gas_temperature_increment_k: self.maximum_gas_temperature_increment_k };
        let mut evaluate = |q: f64, trials: usize|
            -> Result<(Trial<S, G>, f64), FiniteGasBodiesError<E>> {
            let heat = self.gas_heat_j + q;
            if !heat.is_finite() || (q != 0.0 && heat == self.gas_heat_j)
                || (self.gas_heat_j != 0.0 && heat == q) {
                return Err(Error::Invalid("shared gas source/exchange heat is unrepresentable").into());
            }
            let mut next_gas = G::clone(self.gas);
            let gas = next_gas.add_heat(heat, self.gas_solve, || cx.checkpoint().is_err())
                .map_err(gas_error)?;
            let mut candidates = Vec::with_capacity(ordered.len());
            for body in &ordered {
                cx.checkpoint().map_err(|_| Error::Cancelled)?;
                let convection = resolve(body, next_gas.state())
                    .map_err(|error| FiniteGasBodiesError::Body { id: body.id, error: ConvectionError::Convection(error) })?;
                if !convection.coefficient_w_per_m2_k.is_finite() || convection.coefficient_w_per_m2_k < 0.0 {
                    return Err(body_error(body.id, Error::Invalid("finite nonnegative convection coefficient required")));
                }
                cx.checkpoint().map_err(|_| Error::Cancelled)?;
                let environment = LumpedThermalEnvironment { temperature_k: next_gas.temperature_k(),
                    convection_w_per_m2_k: convection.coefficient_w_per_m2_k, ..body.environment.clone() };
                let candidate = environment.advance(cx, body.input)
                    .map_err(|error| body_error(body.id, Error::Interface(error)))?;
                let endpoint = candidate.report.samples().last().ok_or_else(|| body_error(body.id,
                    Error::Invalid("native body returned no endpoint")))?;
                let heat_j = -duration * endpoint.convection_into_body_w;
                if !heat_j.is_finite() || (heat_j == 0.0 && endpoint.convection_into_body_w != 0.0) {
                    return Err(body_error(body.id, Error::Invalid("body convective heat is unrepresentable")));
                }
                candidates.push(BodyTrial { body: candidate, convection, heat_j });
            }
            let residual = sum(std::iter::once(q).chain(candidates.iter().map(|body| -body.heat_j)))?;
            Ok((Trial { next_gas, gas, bodies: candidates, heat_j: q, residual_j: residual, trials }, residual))
        };
        // Exact gas idle can still carry simultaneous opposite body transfers.
        // Evaluate ALL bodies and reuse that proposal if a scalar solve follows.
        let zero = if self.gas_heat_j == 0.0 { Some(evaluate(0.0, 1)?) } else { None };
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        let solved = match zero {
            Some((candidate, residual)) if residual.abs() <= self.exchange.source_tolerance_j => candidate,
            mut zero => {
                let bracket = bounds.supported_heat_bracket()?;
                if zero.is_some() && (bracket[0] > 0.0 || bracket[1] < 0.0) {
                    return Err(Error::Invalid("gas caloric storage domain excludes zero heat").into());
                }
                feedback::solve_scalar_interval(bracket, self.exchange.source_tolerance_j,
                    self.exchange.maximum_iterations, |q, trials, _| {
                        if let Some(first) = zero.take() { Ok(first) } else { evaluate(q, trials) }
                    }, &mut || cx.checkpoint().is_err()).map_err(scalar_error)?
            }
        };
        let gas_k = solved.next_gas.temperature_k();
        if (gas_k - self.gas.temperature_k()).abs() > self.maximum_gas_temperature_increment_k {
            return Err(Error::Invalid("shared gas temperature increment exceeds its limit").into());
        }
        let gas_storage = self.gas.mass_kg() * self.gas
            .specific_internal_energy_change_to_temperature_j_kg(gas_k).map_err(gas_error)?;
        let gas_residual = sum([gas_storage, -self.gas_heat_j, -solved.heat_j])?;
        let mut accepted = Vec::with_capacity(ordered.len());
        for (original, trial) in ordered.iter().zip(solved.bodies) {
            cx.checkpoint().map_err(|_| Error::Cancelled)?;
            let input = original.input;
            let endpoint = trial.body.report.samples().last().ok_or_else(|| body_error(original.id,
                Error::Invalid("native body returned no endpoint")))?;
            let radiation = duration * endpoint.radiation_into_body_w;
            let storage = input.mass_kg * (trial.body.state.specific_enthalpy_j_kg()
                - input.initial.specific_enthalpy_j_kg());
            let residual = sum([storage, -input.internal_heat_j, -radiation, trial.heat_j])?;
            let difference = trial.body.state.temperature_k() - gas_k;
            if residual.abs() > self.exchange.energy_tolerance_j
                || (trial.heat_j > 0.0 && difference < 0.0) || (trial.heat_j < 0.0 && difference > 0.0) {
                return Err(body_error(original.id, Error::Invalid("native body heat fails energy or direction admission")));
            }
            let terms = entropy::body_terms(input, trial.body.state, gas_k,
                original.environment.radiation_temperature_k, storage, [trial.heat_j, radiation], [0.0, residual])
                .map_err(|error| body_error(original.id, error))?;
            accepted.push(UniformEnthalpyStep { state: trial.body.state,
                external_heat_j: sum([radiation, -trial.heat_j])?,
                energy_residual_tolerance_j: self.exchange.energy_tolerance_j,
                report: FiniteGasBodyReport { id: original.id, body: trial.body.report,
                    convection: trial.convection, solid_to_gas_heat_j: trial.heat_j,
                    radiation_heat_j: radiation, enthalpy_change_j: storage, energy_residual_j: residual,
                    entropy: FiniteGasBodyEntropyReport {
                        storage_change_j_k: terms.storage_change_j_k[0],
                        convection_production_j_k: terms.convection_production_j_k,
                        radiation_reservoir_input_j_k: terms.radiation_reservoir_input_j_k,
                        radiation_transfer_production_j_k: terms.radiation_transfer_production_j_k,
                        unclassified_heat_at_endpoint_j_k: terms.unclassified_heat_at_endpoint_j_k[0],
                        temporal_defect_j_k: terms.temporal_defect_j_k[0],
                        algebraic_residual_j_k: terms.algebraic_residual_j_k[0],
                    },
                },
            });
        }
        let heat_from_bodies = sum(accepted.iter().map(|b| b.report.solid_to_gas_heat_j))?;
        let radiation = sum(accepted.iter().map(|b| b.report.radiation_heat_j))?;
        let body_storage = sum(accepted.iter().map(|b| b.report.enthalpy_change_j))?;
        let internal_heat = sum(ordered.iter().map(|b| b.input.internal_heat_j))?;
        let residual = sum([body_storage, gas_storage, -internal_heat, -self.gas_heat_j, -radiation])?;
        let absolute = sum([gas_residual.abs(), solved.residual_j.abs()].into_iter()
            .chain(accepted.iter().map(|b| b.report.energy_residual_j.abs())))?;
        if residual.abs() > self.exchange.energy_tolerance_j || absolute > self.exchange.energy_tolerance_j {
            return Err(Error::Invalid("shared gas/body storage exceeds its whole-account energy allowance").into());
        }
        let gas_entropy = self.gas.mass_kg() * self.gas.specific_entropy_change_to_temperature_j_kg_k(gas_k)
            .map_err(GasEnergyError::from).map_err(gas_error)?;
        let mut entropy = FiniteGasEntropyReport {
            storage_change_j_k: [sum(accepted.iter().map(|b| b.report.entropy.storage_change_j_k))?, gas_entropy],
            convection_production_j_k: sum(accepted.iter().map(|b| b.report.entropy.convection_production_j_k))?,
            radiation_reservoir_input_j_k: sum(accepted.iter().map(|b| b.report.entropy.radiation_reservoir_input_j_k))?,
            radiation_transfer_production_j_k: sum(accepted.iter().map(|b| b.report.entropy.radiation_transfer_production_j_k))?,
            unclassified_heat_at_endpoint_j_k: [sum(accepted.iter().map(|b| b.report.entropy.unclassified_heat_at_endpoint_j_k))?, self.gas_heat_j / gas_k],
            temporal_defect_j_k: [sum(accepted.iter().map(|b| b.report.entropy.temporal_defect_j_k))?, gas_entropy - gas_storage / gas_k],
            source_residual_j_k: solved.residual_j / gas_k,
            algebraic_residual_j_k: [sum(accepted.iter().map(|b| b.report.entropy.algebraic_residual_j_k))?, gas_residual / gas_k],
            balance_residual_j_k: 0.0,
        };
        entropy::finish_balance(&mut entropy)?;
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        Ok(FiniteGasBodiesCandidate { gas: solved.next_gas, bodies: accepted,
            report: FiniteGasBodiesReport { gas: solved.gas, heat_into_gas_j: solved.heat_j,
                heat_from_bodies_j: heat_from_bodies, source_residual_j: solved.residual_j,
                external_gas_heat_j: self.gas_heat_j, internal_body_heat_j: internal_heat,
                radiation_heat_j: radiation, body_enthalpy_change_j: body_storage,
                gas_internal_energy_change_j: gas_storage, combined_energy_residual_j: residual,
                absolute_energy_defect_j: absolute, entropy, exchange_trials: solved.trials,
                body_trials: solved.trials * ordered.len(),
            },
        })
    }
}

// Fixed canonical-ID order and compensation preserve opposing real transfers.
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, Error> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let (next, residual) = fs_math::eft::two_sum(total, value);
        total = next;
        correction += residual;
        if !total.is_finite() || !correction.is_finite() {
            return Err(Error::Invalid("shared gas/body account is unrepresentable"));
        }
    }
    let total = total + correction;
    if total.is_finite() { Ok(total) }
    else { Err(Error::Invalid("shared gas/body account is unrepresentable")) }
}
