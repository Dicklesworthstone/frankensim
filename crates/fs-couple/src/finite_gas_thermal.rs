//! Finite source gas on the existing uniform-body thermal callback.
//!
//! The original lumped enthalpy owner supplies conductivity, radiation, Biot
//! admission and the body's implicit endpoint. The original gas owner supplies
//! internal energy, EOS and transport. The shared bounded feedback solver finds
//! one reciprocal convective transfer Q, positive from the body into the gas.
//! Maxwell, plasticity and other callers may repeat this entire proposal during
//! their own thermal feedback. Neither accepted body nor gas changes here.

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
    },
};

/// Original thermal, gas or bounded-feedback refusal. No piston is introduced.
pub type Error = ThermalPistonError<ConductionError>;

/// One finite gas inventory and the existing body boundary law for a proposal.
///
/// Only `body.temperature_k` is replaced, by the gas endpoint at each trial.
/// The radiation temperature remains an independent external enclosure.
/// The declared film coefficient is not inferred from a gas name: its physical
/// applicability, mixing and flow regime remain the caller's responsibility.
#[derive(Clone, Debug)]
pub struct FiniteGasThermalEnvironment<'a> {
    pub gas: &'a IsochoricGasParcel,
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
pub struct FiniteGasThermalReport {
    /// Publish only after the enclosing mechanical/thermal proposal accepts.
    pub next_gas: IsochoricGasParcel,
    /// Original gas storage and inversion report.
    pub gas: GasHeatStep,
    /// Original body trajectory, including actual endpoint boundary powers.
    pub body: LumpedEnthalpyMarch,
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

struct Candidate {
    next_gas: IsochoricGasParcel,
    gas: GasHeatStep,
    body: UniformEnthalpyStep<LumpedEnthalpyMarch>,
    heat_j: f64,
    source_residual_j: f64,
    trials: usize,
}

fn gas_error(error: GasEnergyError) -> Error {
    PistonError::from(error).into()
}

impl FiniteGasThermalEnvironment<'_> {
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
    ) -> Result<UniformEnthalpyStep<FiniteGasThermalReport>, Error> {
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        self.gas_solve.validate().map_err(gas_error)?;
        feedback::admit::<ConductionError>(self.exchange)?;
        if !self.gas_heat_j.is_finite()
            || !self.maximum_gas_temperature_increment_k.is_finite()
            || self.maximum_gas_temperature_increment_k <= 0.0
        {
            return Err(Error::Invalid(
                "finite prescribed gas heat and positive gas temperature-increment limit required",
            ));
        }
        let solved = feedback::solve(
            self.exchange,
            |q, trials, cancelled| {
                let gas_heat = self.gas_heat_j + q;
                if !gas_heat.is_finite()
                    || (q != 0.0 && gas_heat == self.gas_heat_j)
                    || (self.gas_heat_j != 0.0 && gas_heat == q)
                {
                    return Err(Error::Invalid(
                        "gas source or exchange heat overflows or disappears when combined",
                    ));
                }
                let mut next_gas = self.gas.clone();
                let gas = next_gas
                    .add_heat(gas_heat, self.gas_solve, &mut *cancelled)
                    .map_err(gas_error)?;
                let mut environment = self.body.clone();
                environment.temperature_k = next_gas.state().state().temperature;
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
                    return Err(Error::Invalid("convective heat is unrepresentable"));
                }
                Ok((
                    Candidate {
                        next_gas,
                        gas,
                        body,
                        heat_j: q,
                        source_residual_j: residual,
                        trials,
                    },
                    residual,
                ))
            },
            &mut || cx.checkpoint().is_err(),
        )?;
        let Candidate {
            next_gas,
            gas,
            body,
            heat_j: q,
            source_residual_j,
            trials,
        } = solved;
        let temperature_difference =
            body.state.temperature_k() - next_gas.state().state().temperature;
        if (q > 0.0 && temperature_difference < 0.0)
            || (q < 0.0 && temperature_difference > 0.0)
            || (next_gas.state().state().temperature - self.gas.state().state().temperature).abs()
                > self.maximum_gas_temperature_increment_k
        {
            return Err(Error::Invalid(
                "accepted exchange crosses endpoint temperatures or exceeds the gas increment limit",
            ));
        }
        let endpoint = body.report.samples().last().ok_or(Error::Invalid(
            "thermal owner returned no accepted endpoint",
        ))?;
        let radiation_heat_j = input.duration_s * endpoint.radiation_into_body_w;
        let external_heat_j = radiation_heat_j - q;
        let body_storage = input.mass_kg
            * (body.state.specific_enthalpy_j_kg()
                - input.initial.specific_enthalpy_j_kg());
        let gas_storage = self.gas.mass_kg()
            * (next_gas.state().specific_internal_energy_j_kg()
                - self.gas.state().specific_internal_energy_j_kg());
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
            ));
        }
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        Ok(UniformEnthalpyStep {
            state: body.state,
            external_heat_j,
            energy_residual_tolerance_j: self.exchange.energy_tolerance_j,
            report: FiniteGasThermalReport {
                next_gas,
                gas,
                body: body.report,
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
