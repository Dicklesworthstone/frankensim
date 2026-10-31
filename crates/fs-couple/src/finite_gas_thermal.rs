//! Finite source gas on the existing uniform-body thermal callback.
//!
//! The original lumped enthalpy owner supplies conductivity, radiation, Biot
//! admission and the body's implicit endpoint. The original gas owner supplies
//! internal energy, EOS and transport. The shared bounded feedback solver finds
//! one reciprocal convective transfer Q, positive from the body into the gas.
//! Maxwell, plasticity and other callers may repeat this entire proposal during
//! their own thermal feedback. Neither accepted body nor gas changes here.

mod entropy;
pub use entropy::FiniteGasEntropyReport;

/// Finite water condensation/evaporation in the same reciprocal thermal solve.
pub mod condensation;

/// Dry onset and complete evaporation of a retained film on a finite wall.
pub mod wall_film;

mod wall;
pub use wall::{WallConvectionError, WallConvectionSample, WallConvectionSolve};

use crate::cavity_phs::piston::{
    PistonError,
    thermal::{HeatExchangeSolve, ThermalPistonError, feedback},
};
use fs_conduction::{
    ConductionError,
    lumped::{LumpedEnthalpyMarch, LumpedThermalEnvironment},
};
use fs_exec::Cx;
use fs_material::{
    phase::{UniformEnthalpyStep, UniformEnthalpyStepInput},
    thermochemical_gas::energy::{
        GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel,
        storage::IsochoricGasStorage,
    },
};

/// Original thermal, gas or bounded-feedback refusal. No piston is introduced.
pub type Error = ThermalPistonError<ConductionError>;

/// Original coupled-thermal refusal or the caller's concrete convection law.
#[derive(Debug)]
pub enum ConvectionError<E> {
    Thermal(Error),
    Convection(E),
}
impl<E> From<Error> for ConvectionError<E> {
    fn from(error: Error) -> Self { Self::Thermal(error) }
}
impl<E: core::fmt::Display> core::fmt::Display for ConvectionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "{error}"),
            Self::Convection(error) => write!(f, "finite-gas convection: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ConvectionError<E> {}

/// Whole-surface coefficient and its retained constitutive evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvectionCoefficient<S> {
    /// Nonnegative h [W/(m² K)]; zero explicitly insulates the body.
    pub coefficient_w_per_m2_k: f64,
    pub source: S,
}

/// One finite gas inventory and the existing body boundary law for a proposal.
///
/// `body.temperature_k` is replaced by the gas endpoint at each trial.
/// The radiation temperature remains an independent external enclosure.
/// The declared film coefficient is not inferred from a gas name: its physical
/// applicability, mixing and flow regime remain the caller's responsibility.
/// `advance_with_convection` instead resolves h from the same candidate gas.
/// Both ideal and nonideal parcels use their original caloric/EOS owners;
/// the default type preserves the ideal-gas call sites.
#[derive(Clone, Debug)]
pub struct FiniteGasThermalEnvironment<'a, G = IsochoricGasParcel> {
    pub gas: &'a G,
    pub body: &'a LumpedThermalEnvironment,
    /// Prescribed heat into the gas over this interval [J], excluding Q.
    pub gas_heat_j: f64,
    pub gas_solve: GasEnergySolve,
    pub exchange: HeatExchangeSolve,
    /// Accepted gas temperature change allowed over one interval [K].
    /// This is an admission limit, not a temporal-error estimate.
    pub maximum_gas_temperature_increment_k: f64,
}

/// Both original owner reports and the candidate finite gas at one endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasThermalReport<S = (), G = IsochoricGasParcel> {
    /// Publish only after the enclosing mechanical/thermal proposal accepts.
    pub next_gas: G,
    /// Original gas storage and inversion report.
    pub gas: GasHeatStep,
    /// Original body trajectory, including actual endpoint boundary powers.
    pub body: LumpedEnthalpyMarch,
    /// Exact thermal storage, transfer and numerical entropy account for this
    /// accepted proposal. No additional body/gas feedback solve is performed.
    pub entropy: FiniteGasEntropyReport,
    /// Coefficient actually used at the accepted gas endpoint and its sources.
    pub convection: ConvectionCoefficient<S>,
    /// One signed interface transfer, removed from body and added to gas [J].
    pub solid_to_gas_heat_j: f64,
    /// Caller-prescribed gas heat, excluding the reciprocal transfer [J].
    pub external_gas_heat_j: f64,
    /// Independent enclosing radiation field's heat into the body [J].
    pub radiation_heat_j: f64,
    /// Q + dt times the original owner's convective power into the body [J].
    pub source_residual_j: f64,
    /// Measured body storage - internal heat - radiation heat + Q [J].
    pub solid_energy_residual_j: f64,
    /// Measured body + gas storage - internal body/gas sources - radiation [J].
    /// Internal Q cancels exactly from this accounting boundary.
    pub combined_energy_residual_j: f64,
    /// Complete body/gas proposals, at most maximum_iterations + 2.
    pub trials: usize,
}

struct Candidate<S, G> {
    next_gas: G,
    gas: GasHeatStep,
    body: UniformEnthalpyStep<LumpedEnthalpyMarch>,
    convection: ConvectionCoefficient<S>,
    heat_j: f64,
    source_residual_j: f64,
    trials: usize,
}

fn gas_error(error: GasEnergyError) -> Error {
    PistonError::from(error).into()
}

fn feedback_error<E>(error: feedback::ScalarError<ConvectionError<E>>) -> ConvectionError<E> {
    let thermal = match error {
        feedback::ScalarError::Evaluation(error) => return error,
        feedback::ScalarError::Invalid(_) => Error::Invalid("nonfinite physical energy feedback residual"),
        feedback::ScalarError::Cancelled => Error::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => Error::OutsideBracket {
            lower_residual_j: lower, upper_residual_j: upper,
        },
        feedback::ScalarError::NotConverged { iterations, residual } => Error::NotConverged {
            iterations, residual_j: residual,
        },
    };
    thermal.into()
}

impl<G: IsochoricGasStorage> FiniteGasThermalEnvironment<'_, G> {
    /// Propose reciprocal convection while retaining the shared thermal API.
    ///
    /// Each trial starts from the SAME accepted gas and body state. The gas
    /// receives prescribed gas heat + Q; the existing body owner solves its
    /// backward-Euler enthalpy equation against that trial gas temperature.
    /// The interface residual is Q - dt h A (T_body,end - T_gas,end).
    /// No constant body heat capacity, averaged temperature or replacement gas
    /// law is used. Source-domain and original Biot refusals propagate.
    ///
    /// The returned boundary heat is radiation_heat - Q, so exactly the same
    /// transfer appears on both sides. The original body's evaluated convection
    /// remains inspectable in `report.body`; its discrepancy from Q is the
    /// explicit source residual. Measured body and combined energy residuals
    /// must separately meet `exchange.energy_tolerance_j` after convergence.
    /// Numerical mismatch is never relabeled as heat or corrected by rescaling.
    ///
    /// Radiation is exchanged with its independently prescribed enclosure,
    /// never deposited into this gas. Spatial gas gradients, pressure work,
    /// gas motion, reaction and acoustic loading are outside this fixed-volume,
    /// well-mixed thermal coupling. Phase and geometry transitions remain with
    /// the caller that accepts the `UniformEnthalpyStep`.
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
    ) -> Result<UniformEnthalpyStep<FiniteGasThermalReport<(), G>>, Error> {
        self.advance_with_convection(cx, input, |_| {
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                coefficient_w_per_m2_k: self.body.convection_w_per_m2_k,
                source: (),
            })
        }).map_err(|error| match error {
            ConvectionError::Thermal(error) => error,
            ConvectionError::Convection(never) => match never {},
        })
    }

    /// Use a state-dependent whole-surface convection coefficient in the SAME
    /// reciprocal heat solve. `resolve` sees the actual candidate gas endpoint,
    /// including its source-selected Cp, viscosity, conductivity and composition.
    /// Its result replaces only `body.convection_w_per_m2_k`; the original body
    /// conductivity, radiation, enthalpy and Biot checks remain active.
    ///
    /// The repeatable callback owns correlation/geometry validity and returns
    /// its concrete source evaluation. This seam admits laws of the gas endpoint
    /// and prescribed geometry/flow, not laws requiring an unknown wall endpoint
    /// or a film temperature. Arbitrary state-dependent h has no monotonicity or
    /// unconditional convergence guarantee; bracket/work/refusal limits apply.
    /// No callback effects are rolled back. Neither physical input is mutated.
    pub fn advance_with_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        mut resolve: impl FnMut(&G::State) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<UniformEnthalpyStep<FiniteGasThermalReport<S, G>>, ConvectionError<E>> {
        self.advance_with_convection_feedback(cx, input, |gas, _| resolve(gas))
    }

    /// Resolve the same convection law with its proposed signed interface heat.
    ///
    /// `resolve(gas_endpoint, heat_into_gas_j)` may select heating/cooling
    /// branches from the actual proposed exchange. At acceptance that sign must
    /// agree with the original owner's endpoint temperature ordering. At zero
    /// exchange a direction-dependent law must supply an admissible neutral
    /// seed; either branch has zero flux at exact thermal equilibrium.
    ///
    /// All source, energy, Biot, radiation and transactional guarantees of
    /// `advance_with_convection` apply. The wall endpoint remains with the
    /// original body solver; this does not admit a film-temperature law or
    /// arbitrary dependence on an unresolved wall state.
    pub fn advance_with_convection_feedback<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        mut resolve: impl FnMut(&G::State, f64)
            -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<UniformEnthalpyStep<FiniteGasThermalReport<S, G>>, ConvectionError<E>> {
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        self.gas_solve.validate().map_err(gas_error)?;
        feedback::admit::<ConductionError>(self.exchange)?;
        if !self.gas_heat_j.is_finite()
            || !self.maximum_gas_temperature_increment_k.is_finite()
            || self.maximum_gas_temperature_increment_k <= 0.0
        {
            return Err(Error::Invalid(
                "finite prescribed gas heat and positive gas temperature-increment limit required",
            ).into());
        }
        let solved = feedback::solve_scalar(
            self.exchange.heat_bracket_j,
            self.exchange.source_tolerance_j,
            self.exchange.maximum_iterations,
            |q, trials, cancelled| {
                let gas_heat = self.gas_heat_j + q;
                if !gas_heat.is_finite()
                    || (q != 0.0 && gas_heat == self.gas_heat_j)
                    || (self.gas_heat_j != 0.0 && gas_heat == q)
                {
                    return Err(Error::Invalid(
                        "gas source or exchange heat overflows or disappears when combined",
                    ).into());
                }
                let mut next_gas = G::clone(self.gas);
                let gas = next_gas
                    .add_heat(gas_heat, self.gas_solve, &mut *cancelled)
                    .map_err(gas_error)?;
                let mut environment = self.body.clone();
                environment.temperature_k = next_gas.temperature_k();
                let convection = resolve(next_gas.state(), q).map_err(ConvectionError::Convection)?;
                if !convection.coefficient_w_per_m2_k.is_finite()
                    || convection.coefficient_w_per_m2_k < 0.0
                {
                    return Err(Error::Invalid("convection coefficient must be finite and nonnegative").into());
                }
                cx.checkpoint().map_err(|_| Error::Cancelled)?;
                environment.convection_w_per_m2_k = convection.coefficient_w_per_m2_k;
                let body = environment.advance(cx, input).map_err(Error::Interface)?;
                let endpoint = body.report.samples().last().ok_or(Error::Invalid(
                    "thermal owner returned no accepted endpoint",
                ))?;
                let convective_heat = input.duration_s * endpoint.convection_into_body_w;
                let residual = q + convective_heat;
                if !convective_heat.is_finite()
                    || !residual.is_finite()
                    || (endpoint.convection_into_body_w != 0.0 && convective_heat == 0.0)
                {
                    return Err(Error::Invalid("convective heat is unrepresentable").into());
                }
                Ok((
                    Candidate {
                        next_gas,
                        gas,
                        body,
                        convection,
                        heat_j: q,
                        source_residual_j: residual,
                        trials,
                    },
                    residual,
                ))
            },
            &mut || cx.checkpoint().is_err(),
        ).map_err(feedback_error)?;
        let Candidate {
            next_gas,
            gas,
            body,
            convection,
            heat_j: q,
            source_residual_j,
            trials,
        } = solved;
        let temperature_difference =
            body.state.temperature_k() - next_gas.temperature_k();
        if (q > 0.0 && temperature_difference < 0.0)
            || (q < 0.0 && temperature_difference > 0.0)
            || (next_gas.temperature_k() - self.gas.temperature_k()).abs()
                > self.maximum_gas_temperature_increment_k
        {
            return Err(Error::Invalid(
                "accepted exchange crosses endpoint temperatures or exceeds the gas increment limit",
            ).into());
        }
        let endpoint = body.report.samples().last().ok_or(Error::Invalid(
            "thermal owner returned no accepted endpoint",
        ))?;
        let radiation_heat_j = input.duration_s * endpoint.radiation_into_body_w;
        let external_heat_j = radiation_heat_j - q;
        let body_storage = input.mass_kg
            * (body.state.specific_enthalpy_j_kg()
                - input.initial.specific_enthalpy_j_kg());
        // Replay the original fixed-density caloric increment. Subtracting
        // absolute state energies can erase small heat beside formation energy.
        let gas_storage = self.gas.mass_kg()
            * self.gas.specific_internal_energy_change_to_temperature_j_kg(
                next_gas.temperature_k(),
            ).map_err(gas_error)?;
        let solid_energy_residual_j =
            body_storage - input.internal_heat_j - external_heat_j;
        let combined_energy_residual_j =
            body_storage + gas_storage - input.internal_heat_j - self.gas_heat_j
                - radiation_heat_j;
        if !radiation_heat_j.is_finite()
            || !external_heat_j.is_finite()
            || !solid_energy_residual_j.is_finite()
            || !combined_energy_residual_j.is_finite()
            || solid_energy_residual_j.abs() > self.exchange.energy_tolerance_j
            || combined_energy_residual_j.abs() > self.exchange.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "finite body/gas storage does not close the declared energy allowance",
            ).into());
        }
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        let entropy = entropy::evaluate(
            self,
            input,
            &next_gas,
            body.state,
            [body_storage, gas_storage],
            [q, radiation_heat_j],
            [source_residual_j, solid_energy_residual_j],
        )?;
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        Ok(UniformEnthalpyStep {
            state: body.state,
            external_heat_j,
            energy_residual_tolerance_j: self.exchange.energy_tolerance_j,
            report: FiniteGasThermalReport {
                next_gas,
                gas,
                body: body.report,
                entropy,
                convection,
                solid_to_gas_heat_j: q,
                external_gas_heat_j: self.gas_heat_j,
                radiation_heat_j,
                source_residual_j,
                solid_energy_residual_j,
                combined_energy_residual_j,
                trials,
            },
        })
    }
}
