//! Reciprocal finite-body storage on the same isobaric chemical endpoint.
//!
//! Q enters gas enthalpy and leaves the original body enthalpy chart exactly
//! once. Body temperature/phase and gas composition/transport feed the SAME
//! endpoint interface law. No split burn, guessed wall capacity or latent-heat
//! source is introduced. The body is uniform, fixed-geometry and negligible-pV
//! work; its material chart, not this adapter, supplies caloric fidelity.

use super::*;
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Shared thermal/chemical work limits and a separate total entropy allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricBodySolve {
    pub heat: IsobaricHeatSolve,
    /// Absolute allowed roundoff in body/contact/combined entropy residuals [J/K].
    /// This is not a temporal error estimate or an external-reservoir model.
    pub entropy_tolerance_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBodyStep {
    pub boundary: IsobaricBoundaryStep,
    /// Source-resolved conductance at the accepted gas/body endpoint [W/K].
    pub conductance_w_k: f64,
    pub thermal_energy_change_j: f64,
    pub thermal_energy_residual_j: f64,
    /// Delta H_gas + Delta H_body - both external heats [J].
    pub total_enthalpy_residual_j: f64,
    /// Delta U_gas + p Delta V + Delta H_body - both external heats [J].
    pub total_energy_residual_j: f64,
    /// Integral of dh/T on the actual original body chart [J/K].
    pub thermal_entropy_change_j_k: f64,
    /// Delta S_body - (Q_body_external - Q)/T_body_end [J/K].
    pub thermal_entropy_residual_j_k: f64,
    /// Applied Q*(1/T_gas_end - 1/T_body_end), distinct from evaluated-law heat.
    pub applied_contact_entropy_j_k: f64,
    /// Delta S_gas+Delta S_body-Q_gas_external/T_gas_end-Q_body_external/T_body_end.
    /// Includes backward-Euler numerical entropy; not isolated chemical
    /// production. External reservoirs need their own physical entropy account.
    pub total_entropy_residual_j_k: f64,
}

/// Both continuing histories from one accepted interval, with no input mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBodyCandidate<S> {
    gas: ClosedGasParcel,
    thermal: UniformEnthalpyStep<S>,
    step: IsobaricBodyStep,
}
impl<S> IsobaricBodyCandidate<S> {
    pub const fn gas(&self) -> &ClosedGasParcel { &self.gas }
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> { &self.thermal }
    pub const fn step(&self) -> &IsobaricBodyStep { &self.step }
    /// Commit both histories together. The returned body already contains -Q.
    pub fn into_parts(self) -> (ClosedGasParcel, UniformEnthalpyStep<S>, IsobaricBodyStep) {
        (self.gas, self.thermal, self.step)
    }
}

impl GasReaction {
    /// Solve constant-pressure chemistry and reciprocal finite-wall heat together.
    ///
    /// Each heat trial starts at the same original gas inventory and original
    /// wall enthalpy. The existing coupled implicit reaction owner solves all
    /// channels at H_gas1=H_gas0+external_heat+Q. The original phase chart supplies
    /// T_body at h_body1=h_body0+(body.internal_heat-Q)/mass. A shared scalar root
    /// closes Q=dt G(T_body1-T_gas1), with actual endpoint source selection.
    ///
    /// Latent plateaus evolve fraction without inventing a temperature change.
    /// Source callbacks must be repeatable and must not mutate live counterparts.
    /// Zero conductance still queries source validity. Every queried proposal
    /// must fit both caloric domains and all work limits. Original gas, body and
    /// sources remain unchanged on any refusal, including late entropy failure.
    ///
    /// Uniform body, fixed geometry and negligible body pV work are explicit
    /// assumptions. This is first-order heat/chemistry feedback, not a resolved
    /// wall profile, flame, open-flow reactor or arbitrary-step stability claim.
    #[allow(clippy::too_many_arguments)]
    pub fn isobaric_thermal_body_candidate<S, E>(
        &self, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: IsobaricBodySolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricBodyCandidate<S>, IsobaricThermalError<E>> {
        self.isobaric_thermal_body_counted(gas, duration_s, external_heat_j, body,
            chemistry, controls, resolve, cancelled, &Cell::new(0))
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::cavity_phs::piston::reactive::implicit::coupled::isobaric)
    fn isobaric_thermal_body_counted<S, E>(
        &self, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: IsobaricBodySolve,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool, chemical_work: &Cell<usize>,
    ) -> Result<IsobaricBodyCandidate<S>, IsobaricThermalError<E>> {
        use IsobaricThermalError::Invalid;
        poll(&mut cancelled)?;
        if !positive(body.mass_kg) || !positive(body.volume_m3) || !positive(body.surface_area_m2)
            || !positive(duration_s) || body.duration_s != duration_s
            || !body.internal_heat_j.is_finite() || !positive(controls.entropy_tolerance_j_k)
        {
            return Err(Invalid("matching positive clocks, body geometry/mass and finite thermal controls required"));
        }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0).map_err(ThermalPistonError::Thermal)? != body.initial {
            return Err(Invalid("body state does not belong to its original enthalpy chart"));
        }
        let tolerance = controls.heat.exchange.energy_tolerance_j;
        let candidate = solve_counted(self, gas, duration_s, external_heat_j, chemistry, controls.heat,
            |q, gas, cancel| {
                poll(cancel)?;
                let heat = body.internal_heat_j - q;
                let dh = heat/body.mass_kg;
                let h1 = h0+dh;
                if !heat.is_finite() || !dh.is_finite() || !h1.is_finite()
                    || (heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0 && heat == body.internal_heat_j)
                    || (body.internal_heat_j != 0.0 && heat == -q)
                {
                    return Err(Invalid("body heat exchange vanished or overflowed in enthalpy"));
                }
                let state = body.curve.state_at_specific_enthalpy(h1).map_err(ThermalPistonError::Thermal)?;
                let source = resolve(gas, state).map_err(ThermalPistonError::Interface)?;
                let evaluated = conductive_heat(duration_s, source.conductance_w_k,
                    state.temperature_k(), gas.state().temperature)?;
                let entropy = contact_entropy(evaluated, state.temperature_k(), gas.state().temperature)?;
                let stored = body.mass_kg*(state.specific_enthalpy_j_kg()-h0);
                let residual = (stored+q)-body.internal_heat_j;
                if !stored.is_finite() || !residual.is_finite() || residual.abs() > tolerance {
                    return Err(Invalid("finite-wall storage did not close its reciprocal heat debit"));
                }
                Ok(HeatSample { heat_j: evaluated, entropy_j_k: entropy,
                    source: (UniformEnthalpyStep { state, external_heat_j: -q,
                        energy_residual_tolerance_j: tolerance, report: source.source },
                        stored, residual, source.conductance_w_k) })
            }, &mut cancelled, chemical_work)?;
        let (next, (thermal, stored, residual, conductance), boundary) = candidate.into_parts();
        let total_h = (boundary.chemistry.enthalpy_change_j+stored)-external_heat_j-body.internal_heat_j;
        let total_u = (boundary.chemistry.boundary.reaction.stored_energy_change_j+stored)
            +boundary.chemistry.boundary.boundary.work_by_gas_j-external_heat_j-body.internal_heat_j;
        poll(&mut cancelled)?;
        let body_entropy = body.mass_kg*body.curve.specific_entropy_change_j_kg_k(
            h0, thermal.state.specific_enthalpy_j_kg()).map_err(ThermalPistonError::Thermal)?;
        let tg = next.state().state().temperature;
        let tb = thermal.state.temperature_k();
        let q = boundary.exchanged_heat_j;
        let body_residual = body_entropy-(body.internal_heat_j-q)/tb;
        let applied_entropy = (q/tg)*((tb-tg)/tb);
        let total_entropy = (boundary.chemistry.boundary.reaction.entropy_change_j_k+body_entropy)
            -external_heat_j/tg-body.internal_heat_j/tb;
        if [total_h, total_u, body_entropy, body_residual, applied_entropy, total_entropy]
            .iter().any(|value| !value.is_finite())
            || total_h.abs() > tolerance || total_u.abs() > tolerance
            || body_residual < -controls.entropy_tolerance_j_k
            || applied_entropy < -controls.entropy_tolerance_j_k
            || total_entropy < -controls.entropy_tolerance_j_k
        {
            return Err(Invalid("gas and finite-wall enthalpy, pressure work or entropy did not close"));
        }
        poll(&mut cancelled)?;
        Ok(IsobaricBodyCandidate { gas: next, thermal,
            step: IsobaricBodyStep { boundary, conductance_w_k: conductance,
                thermal_energy_change_j: stored, thermal_energy_residual_j: residual,
                total_enthalpy_residual_j: total_h, total_energy_residual_j: total_u,
                thermal_entropy_change_j_k: body_entropy, thermal_entropy_residual_j_k: body_residual,
                applied_contact_entropy_j_k: applied_entropy, total_entropy_residual_j_k: total_entropy } })
    }

    /// Reuse the existing source-bound gas/body series conduction adapter.
    /// The body's actual phase state and original material binding are retained;
    /// gas conductivity belongs to the same reacting endpoint as its rates.
    #[allow(clippy::too_many_arguments)]
    pub fn isobaric_thermal_layers_candidate(
        &self, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: IsobaricBodySolve, layers: MaterialConductionLayers<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricBodyCandidate<ResolvedMaterialStatePoint>, IsobaricThermalError<LayerConductionError>> {
        if layers.area_m2 > body.surface_area_m2 {
            return Err(IsobaricThermalError::Invalid("conducting patch exceeds body exposed area"));
        }
        self.isobaric_thermal_body_candidate(gas, duration_s, external_heat_j, body,
            chemistry, controls, |gas, state| layers.resolve(gas, state), cancelled)
    }
}
