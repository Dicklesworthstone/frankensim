//! Entropy-compatible contact between the original gas/piston and a finite body.
//!
//! One signed Q enters the gas and leaves the body. Gas work uses the material
//! discrete Gibbs rule; the body uses ds = dh/T on its original enthalpy chart.
//! The contact uses Q = dt G (T_bar_body - T_bar_gas), with T_bar_body = dh/ds.
//! Thus its evaluated production is dt G (T_bar_body - T_bar_gas)^2
//! / (T_bar_body T_bar_gas), including sensible heat and latent plateaus.
//! Neither an unspecified heat reservoir nor spatial equilibration is inferred.

use super::super::{PistonSupport, discrete_gradient::DiscreteGradientPistonStep};
use super::{
    ConductanceSample, GasPiston, HeatExchangeSolve, PistonInput, PistonSolve,
    ResolvedThermochemicalGasState, ThermalPistonError, feedback, poll, positive,
};
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// Independently reconstructed energy and entropy of the accepted contact.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscreteGradientThermalStep {
    /// Original gas, momentum, pressure work and discrete Gibbs account.
    pub piston: DiscreteGradientPistonStep,
    /// Original material state at midpoint body ENTHALPY. The interface source
    /// uses this unique phase coordinate, including fraction on a latent plateau.
    pub thermal_midpoint: EquilibriumPhaseState,
    /// Body entropy-conjugate temperature dh/ds [K], or T0 for zero body energy.
    /// This, rather than the enthalpy midpoint temperature, drives heat transfer.
    pub thermal_mean_temperature_k: f64,
    /// Actual nonnegative interface conductance [W/K].
    pub conductance_w_k: f64,
    /// Signed heat transferred from the body to the gas [J].
    pub exchanged_heat_j: f64,
    /// Duration times the accepted conductance and conjugate temperature gap [J].
    pub evaluated_heat_j: f64,
    /// Applied minus evaluated interface heat [J].
    pub source_residual_j: f64,
    /// Represented body enthalpy change m_body*(h1-h0) [J].
    pub thermal_energy_change_j: f64,
    /// Body storage plus exchanged heat minus its prescribed source [J].
    pub thermal_energy_residual_j: f64,
    /// Gas internal + wall kinetic + body enthalpy + exterior pressure work
    /// - applied force/support work - both prescribed heat inputs [J].
    /// The single exchanged Q cancels internally.
    pub total_energy_residual_j: f64,
    /// Actual m_body*integral(dh/T) on the original admitted phase chart [J/K].
    pub thermal_entropy_change_j_k: f64,
    /// Actual gas plus body entropy change [J/K]; mechanical storage adds none.
    pub total_entropy_change_j_k: f64,
    /// Q_gas_prescribed/T_bar_gas + Q_body_prescribed/T_bar_body [J/K].
    /// This is the two inputs' contribution to the material entropy account.
    /// It does not identify reservoir entropy or distinguish external heat from
    /// internally dissipated work; the source owner must supply that physics.
    pub prescribed_source_entropy_j_k: f64,
    /// Evaluated nonnegative conduction production [J/K], using the actual
    /// conductance law. Applied exchange differs by the source solve residual.
    pub contact_entropy_production_j_k: f64,
    /// Matter entropy change - prescribed contribution - contact production.
    /// Includes the gas/body state residuals and the interface source residual;
    /// it is not an extra physical entropy source or an integration-error bound.
    pub entropy_residual_j_k: f64,
    /// Numerical entropy allowance from the gas/body energy controls, interface
    /// source tolerance and roundoff in the actual terms [J/K]. This bounds
    /// accepted algebraic closure, not physical data or time-integration error.
    pub entropy_tolerance_j_k: f64,
    /// Complete gas/piston/interface proposals, including rejected heat trials.
    pub exchange_trials: usize,
    /// Material proposals summed across every nested mechanical trial.
    pub material_trials: usize,
}

/// A single proposal containing both endpoints and the concrete interface source.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscreteGradientThermalCandidate<S> {
    piston: GasPiston<ClosedGasParcel>,
    thermal: UniformEnthalpyStep<S>,
    step: DiscreteGradientThermalStep,
}

impl<S> DiscreteGradientThermalCandidate<S> {
    /// Proposed gas, volume and momentum; the input owner remains unchanged.
    pub const fn piston(&self) -> &GasPiston<ClosedGasParcel> {
        &self.piston
    }
    /// Original thermal carrier with the accepted phase state and source.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> {
        &self.thermal
    }
    /// Actual combined energy, entropy and contact balances.
    pub const fn step(&self) -> &DiscreteGradientThermalStep {
        &self.step
    }
    /// Publish both physical states together. Do not apply exchanged heat again.
    pub fn into_parts(
        self,
    ) -> (
        GasPiston<ClosedGasParcel>,
        UniformEnthalpyStep<S>,
        DiscreteGradientThermalStep,
    ) {
        (self.piston, self.thermal, self.step)
    }
}

impl GasPiston<ClosedGasParcel> {
    /// Solve finite-body conduction and discrete-gradient gas/piston motion.
    ///
    /// The two prescribed inputs `input.heat_j` and `body.internal_heat_j`
    /// exclude the reciprocal exchange Q. Their separate entropy contribution
    /// is reported at each material's conjugate temperature; assigning source
    /// or reservoir entropy remains the caller's responsibility. With both
    /// inputs zero, gas/body entropy change equals the nonnegative conduction
    /// production plus the reported numerical residual, even during motion.
    ///
    /// `resolve(gas_thermodynamic_mean, body_enthalpy_midpoint)` preserves the
    /// existing source/phase interface. Conductance is sampled at these actual
    /// states; the thermal FORCE uses the entropy-conjugate body temperature.
    /// An enthalpy midpoint retains a unique phase fraction where an isothermal
    /// melting plateau makes temperature insufficient. No fictitious phase state
    /// is constructed at the mean temperature. The callback must be repeatable;
    /// its irreversible side effects cannot be rolled back.
    ///
    /// The body's original fixed-mass, negligible-pV-work enthalpy storage model
    /// is retained. Geometry/density consistency and uniform-body validity are
    /// caller assumptions; phase change does not silently move the piston.
    /// The gas remains the original frozen-composition ideal mixture.
    ///
    /// Work is bounded by the existing exchange, mechanical and material solve
    /// controls. Every queried heat must lie in all original source domains.
    /// A positive conductance ensures nonnegative evaluated contact production,
    /// not uniqueness or temporal accuracy. A closed fixed-wall pair refuses an
    /// accepted endpoint that reverses its initial temperature ordering: refine
    /// that timestep. Moving or independently driven pairs have no such gate,
    /// because mechanical work or prescribed heat can physically reverse it.
    /// The returned entropy residual includes finite solve/roundoff error.
    /// Failure or cancellation publishes neither state nor a partial port window.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_discrete_gradient_candidate<S, E>(
        &self,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiscreteGradientThermalCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        feedback::admit(exchange)?;
        if !positive(body.mass_kg)
            || !positive(body.volume_m3)
            || !positive(body.surface_area_m2)
            || !positive(body.duration_s)
            || body.duration_s != input.duration_s
            || !body.internal_heat_j.is_finite()
            || !input.heat_j.is_finite()
        {
            return Err(ThermalPistonError::Invalid(
                "invalid discrete-gradient body geometry, clocks or prescribed heat",
            ));
        }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(ThermalPistonError::Invalid(
                "initial thermal state belongs to a different source chart",
            ));
        }
        let mut material_trials = 0usize;
        let candidate = feedback::solve(
            exchange,
            |q, trials, cancelled| {
                let body_heat = body.internal_heat_j - q;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                let gas_heat = input.heat_j + q;
                if !body_heat.is_finite()
                    || !dh.is_finite()
                    || !h1.is_finite()
                    || !gas_heat.is_finite()
                    || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0 && (body_heat == body.internal_heat_j || gas_heat == input.heat_j))
                {
                    return Err(ThermalPistonError::Invalid(
                        "heat exchange disappears or overflows in conserved energy coordinates",
                    ));
                }
                let thermal_midpoint = body
                    .curve
                    .state_at_specific_enthalpy(f64::midpoint(h0, h1))?;
                let endpoint = body.curve.state_at_specific_enthalpy(h1)?;
                poll(cancelled)?;
                let thermal_mean = body
                    .curve
                    .thermodynamic_mean_temperature(
                        h0,
                        h1,
                        body.curve.knots().len() - 1,
                        &mut *cancelled,
                    )
                    .map_err(|error| match error {
                        PhaseStateError::EntropyEvaluationCancelled => {
                            ThermalPistonError::Cancelled
                        }
                        other => ThermalPistonError::Thermal(other),
                    })?;
                let ds = thermal_mean.specific_entropy_change_j_kg_k;
                let body_temperature = thermal_mean.temperature_k;
                let mut next = self.clone();
                let piston = next.step_discrete_gradient(
                    PistonInput {
                        heat_j: gas_heat,
                        ..input
                    },
                    controls,
                    &mut *cancelled,
                )?;
                material_trials += piston.material_trials;
                poll(cancelled)?;
                let sample = resolve(&piston.thermodynamic_mean, thermal_midpoint)
                    .map_err(ThermalPistonError::Interface)?;
                let g = sample.conductance_w_k;
                let gas_temperature = piston.thermodynamic_mean.state().temperature;
                let delta = body_temperature - gas_temperature;
                let power = g * delta;
                let heat = input.duration_s * power;
                let residual = q - heat;
                if !g.is_finite()
                    || g < 0.0
                    || !power.is_finite()
                    || !heat.is_finite()
                    || !residual.is_finite()
                    || (g != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && heat == 0.0)
                {
                    return Err(ThermalPistonError::Invalid(
                        "interface heat law is negative, nonfinite or unrepresentable",
                    ));
                }
                let stored = body.mass_kg * (endpoint.specific_enthalpy_j_kg() - h0);
                let body_residual = (stored + q) - body.internal_heat_j;
                let total = ((piston.gas.stored_energy_change_j + piston.kinetic_energy_change_j)
                    + stored)
                    + piston.exterior_pressure_work_j
                    - piston.external_force_work_j
                    - piston.support_work_j
                    - input.heat_j
                    - body.internal_heat_j;
                let body_entropy = body.mass_kg * ds;
                let matter_entropy = piston.gas.entropy_change_j_k + body_entropy;
                let gas_prescribed_entropy = input.heat_j / gas_temperature;
                let body_prescribed_entropy = body.internal_heat_j / body_temperature;
                let prescribed_entropy = gas_prescribed_entropy + body_prescribed_entropy;
                // This stable form avoids subtracting nearly equal inverse temperatures.
                let production = (heat / gas_temperature) * (delta / body_temperature);
                let entropy_residual = (matter_entropy - prescribed_entropy) - production;
                // Delta S_g = (Q_g+q)/T_g + r_g, while Delta S_b =
                // (Q_b-q+r_H)/T_b up to its original-chart roundoff. The
                // remaining contact defect is (q-Q_evaluated)*(1/T_g-1/T_b).
                // Derive this allowance from those independent controls.
                let entropy_tolerance = self.gas().mass_kg() * controls.gas.energy_tolerance_j_kg
                    / gas_temperature
                    + exchange.energy_tolerance_j / body_temperature
                    + exchange.source_tolerance_j
                        * ((delta / gas_temperature) / body_temperature).abs()
                    + 128.0
                        * f64::EPSILON
                        * (piston.gas.entropy_change_j_k.abs()
                            + body_entropy.abs()
                            + gas_prescribed_entropy.abs()
                            + body_prescribed_entropy.abs()
                            + production.abs());
                if ![
                    stored,
                    body_residual,
                    total,
                    body_entropy,
                    matter_entropy,
                    prescribed_entropy,
                    production,
                    entropy_residual,
                    entropy_tolerance,
                ]
                .iter()
                .all(|value| value.is_finite())
                    || production < 0.0
                    || entropy_tolerance < 0.0
                    || body_residual.abs() > exchange.energy_tolerance_j
                    || total.abs() > exchange.energy_tolerance_j
                {
                    return Err(ThermalPistonError::Invalid(
                        "finite-body energy or entropy account is nonfinite or does not close",
                    ));
                }
                poll(cancelled)?;
                Ok((
                    DiscreteGradientThermalCandidate {
                        piston: next,
                        thermal: UniformEnthalpyStep {
                            state: endpoint,
                            external_heat_j: -q,
                            energy_residual_tolerance_j: exchange.energy_tolerance_j,
                            report: sample.source,
                        },
                        step: DiscreteGradientThermalStep {
                            piston,
                            thermal_midpoint,
                            thermal_mean_temperature_k: body_temperature,
                            conductance_w_k: g,
                            exchanged_heat_j: q,
                            evaluated_heat_j: heat,
                            source_residual_j: residual,
                            thermal_energy_change_j: stored,
                            thermal_energy_residual_j: body_residual,
                            total_energy_residual_j: total,
                            thermal_entropy_change_j_k: body_entropy,
                            total_entropy_change_j_k: matter_entropy,
                            prescribed_source_entropy_j_k: prescribed_entropy,
                            contact_entropy_production_j_k: production,
                            entropy_residual_j_k: entropy_residual,
                            entropy_tolerance_j_k: entropy_tolerance,
                            exchange_trials: trials,
                            material_trials,
                        },
                    },
                    residual,
                ))
            },
            &mut cancelled,
        )?;
        // Bracket trials may have a large source residual. Admit the actual
        // entropy balance only after the existing heat-law solve has accepted.
        let report = &candidate.step;
        if report.entropy_residual_j_k.abs() > report.entropy_tolerance_j_k
            || (input.heat_j == 0.0
                && body.internal_heat_j == 0.0
                && report.total_entropy_change_j_k < -report.entropy_tolerance_j_k)
        {
            return Err(ThermalPistonError::Invalid(
                "finite-body contact entropy balance exceeds the energy-derived numerical allowance",
            ));
        }
        if self.support == PistonSupport::Fixed
            && input.heat_j == 0.0
            && body.internal_heat_j == 0.0
        {
            let initial_gap = body.initial.temperature_k() - self.gas().state().state().temperature;
            let final_gap = candidate.thermal.state.temperature_k()
                - candidate.piston.gas().state().state().temperature;
            if (initial_gap > 0.0 && final_gap < 0.0) || (initial_gap < 0.0 && final_gap > 0.0) {
                return Err(ThermalPistonError::Invalid(
                    "closed fixed-wall heat contact crossed thermal equilibrium; refine the timestep",
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(candidate)
    }
}

#[cfg(test)]
mod tests;
