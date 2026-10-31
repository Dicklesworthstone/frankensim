//! Finite wet-gas storage coupled to the original enthalpy body.
//!
//! Condensate stays in the material owner's finite liquid inventory. Gas and
//! liquid share a temperature and a consistent caloric datum; cooling through
//! the dew point changes both composition and storage. This is homogeneous
//! equilibrium condensation, not a wall-film transport or nucleation model.
//! The fixed volume is the gas volume: liquid displacement, pressure work and
//! mechanical/acoustic feedback remain outside this thermal accounting window.

use super::ConvectionCoefficient;
use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, feedback};
use fs_conduction::{
    ConductionError,
    lumped::{LumpedEnthalpyMarch, LumpedThermalEnvironment},
};
use fs_exec::Cx;
use fs_material::{
    phase::condensation::{WaterCondensationError, WaterCondensationStep, WetGasParcel},
    phase::{UniformEnthalpyStep, UniformEnthalpyStepInput},
    thermochemical_gas::ResolvedThermochemicalGasState,
};

/// Absolute energy admission for the complete gas plus retained liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WetGasHeatSolve {
    /// Admitted absolute temperature interval [K], containing the initial T.
    pub temperature_bracket_k: [f64; 2],
    /// Maximum actual gas-plus-liquid storage defect [J].
    pub energy_tolerance_j: f64,
    /// At most 128 interior safeguarded scalar iterations.
    pub maximum_iterations: usize,
}

/// Original phase, thermal or caller-supplied convection refusal.
#[derive(Debug)]
pub enum Error<E = core::convert::Infallible> {
    /// The source-bound mixed-gas/phase owner refused a candidate.
    Phase(WaterCondensationError),
    /// Original body or bounded energy-feedback failure.
    Thermal(super::Error),
    /// The supplied physical convection law refused its candidate state.
    Convection(E),
}
impl<E: core::fmt::Display> core::fmt::Display for Error<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Phase(error) => write!(f, "wet-gas phase: {error}"),
            Self::Thermal(error) => write!(f, "wet-gas thermal: {error}"),
            Self::Convection(error) => write!(f, "wet-gas convection: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for Error<E> {}
impl<E> Error<E> {
    fn from_heat(error: Error) -> Self {
        match error {
            Error::Phase(error) => Self::Phase(error),
            Error::Thermal(error) => Self::Thermal(error),
            Error::Convection(never) => match never {},
        }
    }
}
fn invalid<E>(message: &'static str) -> Error<E> {
    Error::Thermal(super::Error::Invalid(message))
}
fn scalar_error<E>(error: feedback::ScalarError<Error<E>>) -> Error<E> {
    Error::Thermal(match error {
        feedback::ScalarError::Evaluation(error) => return error,
        feedback::ScalarError::Invalid(message) => super::Error::Invalid(message),
        feedback::ScalarError::Cancelled => super::Error::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => super::Error::OutsideBracket {
            lower_residual_j: lower,
            upper_residual_j: upper,
        },
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => super::Error::NotConverged {
            iterations,
            residual_j: residual,
        },
    })
}

/// Accepted energy inversion, retaining the actual phase/species balance.
#[derive(Clone, Debug, PartialEq)]
pub struct WetGasHeatStep {
    /// Publish the gas and retained liquid together from this one endpoint.
    pub phase: WaterCondensationStep,
    /// Prescribed energy into gas plus liquid [J].
    pub supplied_heat_j: f64,
    /// Actual gas-plus-liquid internal-energy change minus supplied heat [J].
    pub energy_residual_j: f64,
    /// Material endpoint evaluations in this inversion.
    pub trials: usize,
}

/// Resolve heat into a finite wet parcel without mutating the accepted state.
///
/// Every trial equilibrates the ORIGINAL total water and carrier inventories.
/// The energy coordinate includes the real change in retained liquid storage;
/// latent heat is never added as a second external source. Zero supplied heat
/// may change an initially nonequilibrium parcel's temperature and phase split.
pub fn heat_candidate(
    parcel: &WetGasParcel,
    heat_j: f64,
    controls: WetGasHeatSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WetGasHeatStep, Error> {
    let t0 = parcel.temperature_k();
    let [lo, hi] = controls.temperature_bracket_k;
    if !heat_j.is_finite()
        || !lo.is_finite()
        || !hi.is_finite()
        || lo <= 0.0
        || lo > t0
        || hi < t0
        || lo >= hi
    {
        return Err(invalid(
            "wet-gas heat needs finite energy and an explicit temperature bracket containing the accepted state",
        ));
    }
    feedback::solve_scalar(
        [lo - t0, hi - t0],
        controls.energy_tolerance_j,
        controls.maximum_iterations,
        |delta_t, trials, cancelled| {
            let phase = parcel
                .at_temperature(t0 + delta_t, &mut *cancelled)
                .map_err(Error::Phase)?;
            let residual = phase.total_internal_energy_change_j - heat_j;
            Ok((
                WetGasHeatStep {
                    phase,
                    supplied_heat_j: heat_j,
                    energy_residual_j: residual,
                    trials,
                },
                residual,
            ))
        },
        &mut cancelled,
    )
    .map_err(scalar_error)
}

/// The existing uniform-body boundary, now with finite wet-gas storage.
#[derive(Clone, Debug)]
pub struct CondensingGasThermalEnvironment<'a> {
    /// Accepted gas and retained liquid, both at their common temperature.
    pub parcel: &'a WetGasParcel,
    /// Original conductivity, radiation and Biot admission. Its ambient T is
    /// replaced by each actual wet-gas endpoint, as in the dry-gas path.
    pub body: &'a LumpedThermalEnvironment,
    /// Independent prescribed heat into gas plus liquid [J].
    pub gas_heat_j: f64,
    /// Inner gas-plus-liquid energy inversion and source-domain limits.
    pub gas_solve: WetGasHeatSolve,
    /// Reciprocal body-to-wet-gas heat and enclosing energy limits.
    pub exchange: HeatExchangeSolve,
    /// Largest admitted wet-gas temperature change per interval [K].
    pub maximum_gas_temperature_increment_k: f64,
}

/// One reciprocal body/gas/liquid endpoint and its measured balances.
#[derive(Clone, Debug, PartialEq)]
pub struct CondensingGasThermalReport<S = ()> {
    /// Gas AND liquid must publish from this same accepted proposal.
    pub gas: WetGasHeatStep,
    /// Original body trajectory, with the actual endpoint boundary powers.
    pub body: LumpedEnthalpyMarch,
    /// Convection evaluated at the accepted gas composition and temperature.
    pub convection: ConvectionCoefficient<S>,
    /// One signed Q, removed from the body and added to gas plus liquid [J].
    pub solid_to_gas_heat_j: f64,
    /// Independent prescribed heat into gas plus liquid [J].
    pub external_gas_heat_j: f64,
    /// Heat from the separate radiative enclosure into the body [J].
    pub radiation_heat_j: f64,
    /// Q minus the actual endpoint convection law [J].
    pub source_residual_j: f64,
    /// Actual body storage minus body sources and boundary heat [J].
    pub solid_energy_residual_j: f64,
    /// Body + gas + retained liquid storage minus external heat [J].
    pub combined_energy_residual_j: f64,
    /// Number of whole body/phase proposals in the reciprocal heat solve.
    pub trials: usize,
}

impl CondensingGasThermalEnvironment<'_> {
    /// Advance with the explicitly supplied whole-surface heat coefficient.
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
    ) -> Result<UniformEnthalpyStep<CondensingGasThermalReport>, Error> {
        self.advance_with_convection(cx, input, |_| {
            Ok(ConvectionCoefficient {
                coefficient_w_per_m2_k: self.body.convection_w_per_m2_k,
                source: (),
            })
        })
    }

    /// Reevaluate transport on the actual changing vapor/carrier composition.
    ///
    /// The original body owner solves its backward-Euler enthalpy equation.
    /// The same bounded scalar feedback as dry-gas coupling enforces
    /// Q = dt h A (T_body,end - T_wet,end). Each trial starts from the same
    /// accepted body and wet parcel. Source, energy, work or cancellation errors
    /// publish neither input. Callback side effects cannot be rolled back.
    ///
    /// No entropy theorem is inherited from the frozen-composition dry-gas
    /// report: this lane's guarantee is the explicitly measured energy balance.
    pub fn advance_with_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<UniformEnthalpyStep<CondensingGasThermalReport<S>>, Error<E>> {
        cx.checkpoint()
            .map_err(|_| Error::Thermal(super::Error::Cancelled))?;
        feedback::admit::<ConductionError>(self.exchange).map_err(Error::Thermal)?;
        if !self.gas_heat_j.is_finite()
            || !self.maximum_gas_temperature_increment_k.is_finite()
            || self.maximum_gas_temperature_increment_k <= 0.0
        {
            return Err(invalid(
                "finite gas heat and positive wet-gas temperature increment limit required",
            ));
        }
        let (gas, body, convection, q, source_residual_j, trials) = feedback::solve_scalar(
            self.exchange.heat_bracket_j,
            self.exchange.source_tolerance_j,
            self.exchange.maximum_iterations,
            |q, trials, cancelled| {
                let gas_heat = self.gas_heat_j + q;
                if !gas_heat.is_finite()
                    || (q != 0.0 && gas_heat == self.gas_heat_j)
                    || (self.gas_heat_j != 0.0 && gas_heat == q)
                {
                    return Err(invalid(
                        "wet-gas source or exchange heat is unrepresentable",
                    ));
                }
                let gas = heat_candidate(self.parcel, gas_heat, self.gas_solve, &mut *cancelled)
                    .map_err(Error::from_heat)?;
                let endpoint = gas.phase.parcel.gas().state();
                let convection = resolve(endpoint).map_err(Error::Convection)?;
                if !convection.coefficient_w_per_m2_k.is_finite()
                    || convection.coefficient_w_per_m2_k < 0.0
                {
                    return Err(invalid(
                        "wet-gas convection coefficient must be finite and nonnegative",
                    ));
                }
                let mut environment = self.body.clone();
                environment.temperature_k = endpoint.state().temperature;
                environment.convection_w_per_m2_k = convection.coefficient_w_per_m2_k;
                let body = environment
                    .advance(cx, input)
                    .map_err(|error| Error::Thermal(super::Error::Interface(error)))?;
                let sample = body
                    .report
                    .samples()
                    .last()
                    .ok_or_else(|| invalid("thermal owner returned no accepted endpoint"))?;
                let convective_heat = input.duration_s * sample.convection_into_body_w;
                let residual = q + convective_heat;
                if !convective_heat.is_finite()
                    || (sample.convection_into_body_w != 0.0 && convective_heat == 0.0)
                {
                    return Err(invalid("wet-gas convective heat is unrepresentable"));
                }
                Ok(((gas, body, convection, q, residual, trials), residual))
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(scalar_error)?;
        let end_t = gas.phase.parcel.temperature_k();
        let difference = body.state.temperature_k() - end_t;
        if (q > 0.0 && difference < 0.0)
            || (q < 0.0 && difference > 0.0)
            || (end_t - self.parcel.temperature_k()).abs()
                > self.maximum_gas_temperature_increment_k
        {
            return Err(invalid(
                "wet-gas exchange crosses endpoint temperatures or exceeds the increment limit",
            ));
        }
        let sample = body
            .report
            .samples()
            .last()
            .ok_or_else(|| invalid("thermal owner returned no accepted endpoint"))?;
        let radiation_heat_j = input.duration_s * sample.radiation_into_body_w;
        let external_heat_j = radiation_heat_j - q;
        let body_storage = input.mass_kg
            * (body.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
        let solid_energy_residual_j = body_storage - input.internal_heat_j - external_heat_j;
        let combined_energy_residual_j = body_storage + gas.phase.total_internal_energy_change_j
            - input.internal_heat_j
            - self.gas_heat_j
            - radiation_heat_j;
        if !external_heat_j.is_finite()
            || !solid_energy_residual_j.is_finite()
            || !combined_energy_residual_j.is_finite()
            || solid_energy_residual_j.abs() > self.exchange.energy_tolerance_j
            || combined_energy_residual_j.abs() > self.exchange.energy_tolerance_j
        {
            return Err(invalid(
                "body, gas and condensate storage do not close the energy allowance",
            ));
        }
        cx.checkpoint()
            .map_err(|_| Error::Thermal(super::Error::Cancelled))?;
        Ok(UniformEnthalpyStep {
            state: body.state,
            external_heat_j,
            energy_residual_tolerance_j: self.exchange.energy_tolerance_j,
            report: CondensingGasThermalReport {
                gas,
                body: body.report,
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
