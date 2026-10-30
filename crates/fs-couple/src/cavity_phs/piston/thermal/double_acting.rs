//! One finite enthalpy/phase store exchanging heat with both moving gas volumes.
//!
//! Both heat laws see the SAME gas/momentum proposal and body energy midpoint.
//! Neither face can spend the original reservoir twice. The original bounded
//! heat-feedback and shared-piston solvers retain all constitutive ownership.

use super::{ConductanceSample, EquilibriumPhaseState, HeatExchangeSolve,
    ResolvedThermochemicalGasState, ThermalPistonError, UniformEnthalpyStep,
    UniformEnthalpyStepInput, feedback, poll, positive};
use super::super::double_acting::{DoubleActingGasPiston, DoublePistonInput,
    DoublePistonSolve, DoublePistonStep};

/// Separate material/interface brackets plus a total bound on nested proposals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleHeatExchangeSolve {
    /// [Primary, opposing] face controls, including separate source tolerances.
    /// The smaller energy tolerance bounds the common thermal and total balance.
    pub faces: [HeatExchangeSolve; 2],
    /// Maximum complete paired gas/momentum/body/interface evaluations across
    /// BOTH nested searches, including every rejected and bracket proposal.
    pub maximum_joint_trials: usize,
}

/// One accepted two-gas, one-body and one-kinetic thermal/mechanical update.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleThermalPistonStep {
    /// Original two-gas/one-momentum result; contains both material midpoints.
    pub piston: DoublePistonStep,
    /// Body phase state at midpoint ENTHALPY, not averaged temperature.
    pub thermal_midpoint: EquilibriumPhaseState,
    /// Actual material/interface conductances [W/K], [primary, opposing].
    pub conductance_w_k: [f64; 2],
    /// Signed heat from the ONE finite body into each gas [J].
    pub exchanged_heat_j: [f64; 2],
    /// Each Q minus dt*G*(T_body_mid-T_gas_mid) [J].
    pub source_residual_j: [f64; 2],
    /// Independently measured body mass times enthalpy change [J].
    pub thermal_energy_change_j: f64,
    /// Body storage + Q_primary + Q_opposing - body internal heat [J].
    pub thermal_energy_residual_j: f64,
    /// Both gas internal energies + ONE kinetic + ONE body enthalpy change,
    /// minus the three prescribed heat sources and external force work [J].
    pub total_energy_residual_j: f64,
    /// All complete coupled proposals, not just the final inner search.
    pub joint_trials: usize,
    /// All gas material proposals across both searches and both gases.
    pub material_trials: usize,
}

/// Publish both gases, their common momentum and the thermal body together.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleThermalPistonCandidate<S> {
    piston: DoubleActingGasPiston,
    thermal: UniformEnthalpyStep<[S; 2]>,
    step: DoubleThermalPistonStep,
}
impl<S> DoubleThermalPistonCandidate<S> {
    /// Tentative shared mechanics, geometry and both source-resolved gases.
    pub const fn piston(&self) -> &DoubleActingGasPiston { &self.piston }
    /// Shared thermal carrier; its report retains BOTH original interface sources.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<[S; 2]> { &self.thermal }
    /// Independent source and energy accounting on the accepted proposal.
    pub const fn step(&self) -> &DoubleThermalPistonStep { &self.step }
    /// Return all endpoints together; do not apply the boundary heat a second time.
    pub fn into_parts(self) -> (DoubleActingGasPiston, UniformEnthalpyStep<[S; 2]>, DoubleThermalPistonStep) {
        (self.piston, self.thermal, self.step)
    }
}

impl DoubleActingGasPiston {
    /// Couple a single finite thermal store to BOTH faces during one physical step.
    ///
    /// The body supplies h->(T,phase), each gas supplies u/EOS/transport/acoustics,
    /// and the existing shared piston supplies reciprocal pressure work. For
    /// each trial [Qa,Qb], body enthalpy changes by (H-Qa-Qb)/mass and both gas
    /// heats include their own Q. Both conductances are resolved only after the
    /// SAME coupled proposal, from its actual gas and body energy midpoints.
    ///
    /// `resolve` must be repeatable and must not debit live external state. It
    /// may reuse the existing layered-conductance/material interface providers.
    /// The returned source evaluations are retained without coefficient caches.
    ///
    /// Nested bounded scalar solves enforce BOTH physical heat residuals. They
    /// do not imply uniqueness for arbitrary state-dependent laws or unconditional
    /// stability. Every inner search starts from the original accepted states.
    /// The total trial cap prevents hidden multiplicative work. Source, geometry,
    /// work-budget or cancellation errors return no partially advanced owner.
    ///
    /// The finite body has fixed geometry and negligible pV/deformation work,
    /// following the existing enthalpy carrier. It can be a separate wall/store;
    /// phase change does not silently change piston rigidity, mass or geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_candidate<S, E>(
        &self, input: DoublePistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve, exchange: DoubleHeatExchangeSolve,
        resolve: impl FnMut(&[ResolvedThermochemicalGasState; 2], EquilibriumPhaseState)
            -> Result<[ConductanceSample<S>; 2], E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleThermalPistonCandidate<S>, ThermalPistonError<E>> {
        self.thermal_candidate_with_midpoint_resistance(
            input, body, controls, exchange, 0.0, resolve, cancelled,
        )
    }

    // The exact mean acoustic law is affine in the SAME displacement velocity.
    // Every heat trial retains this load in the original momentum solve; no
    // acoustic candidate advances inside the nested thermal search.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::cavity_phs::piston) fn thermal_candidate_with_midpoint_resistance<S, E>(
        &self, input: DoublePistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve, exchange: DoubleHeatExchangeSolve,
        resistance_n_s_m: f64,
        mut resolve: impl FnMut(&[ResolvedThermochemicalGasState; 2], EquilibriumPhaseState)
            -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleThermalPistonCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        if !resistance_n_s_m.is_finite() || resistance_n_s_m < 0.0 {
            return Err(ThermalPistonError::Invalid("invalid shared midpoint force resistance"));
        }
        for face in exchange.faces { feedback::admit::<E>(face)?; }
        if exchange.maximum_joint_trials == 0 || !positive(body.mass_kg)
            || !positive(body.volume_m3) || !positive(body.surface_area_m2)
            || !positive(body.duration_s) || body.duration_s != input.duration_s
            || !body.internal_heat_j.is_finite() || !input.heat_j.iter().all(|q| q.is_finite())
        { return Err(ThermalPistonError::Invalid("invalid shared thermal body, clock or total work bound")); }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(ThermalPistonError::Invalid("initial body state belongs to another source chart"));
        }
        let energy_tolerance = exchange.faces[0].energy_tolerance_j.min(exchange.faces[1].energy_tolerance_j);
        let mut joint_trials = 0usize;
        let mut material_trials = 0usize;
        feedback::solve(exchange.faces[0], |qa, _, cancelled| {
            let candidate = feedback::solve(exchange.faces[1], |qb, _, cancelled| {
                if joint_trials >= exchange.maximum_joint_trials {
                    return Err(ThermalPistonError::Invalid("shared thermal joint trial budget exhausted"));
                }
                joint_trials += 1;
                let q = [qa, qb];
                let net = qa + qb;
                let body_heat = body.internal_heat_j - net;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                let hm = h0 + 0.5 * dh;
                let gas_heat = [input.heat_j[0] + qa, input.heat_j[1] + qb];
                if ![net, body_heat, dh, h1, hm, gas_heat[0], gas_heat[1]].iter().all(|x| x.is_finite())
                    || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (net != 0.0 && body_heat == body.internal_heat_j)
                    || (0..2).any(|i| q[i] != 0.0 && gas_heat[i] == input.heat_j[i])
                { return Err(ThermalPistonError::Invalid("shared thermal exchange is unrepresentable in energy coordinates")); }
                let midpoint = body.curve.state_at_specific_enthalpy(hm)?;
                let endpoint = body.curve.state_at_specific_enthalpy(h1)?;
                let mut next = self.clone();
                let piston = next.step_with_midpoint_resistance(
                    DoublePistonInput { heat_j: gas_heat, ..input }, controls,
                    resistance_n_s_m, &mut *cancelled,
                )?;
                material_trials = material_trials.checked_add(piston.material_trials)
                    .ok_or(ThermalPistonError::Invalid("material proposal count overflow"))?;
                poll(cancelled)?;
                let samples = resolve(&piston.midpoint, midpoint).map_err(ThermalPistonError::Interface)?;
                let mut conductance = [0.0; 2];
                let mut residual = [0.0; 2];
                for i in 0..2 {
                    let g = samples[i].conductance_w_k;
                    let delta = midpoint.temperature_k() - piston.midpoint[i].state().temperature;
                    let power = g * delta;
                    let heat = input.duration_s * power;
                    residual[i] = q[i] - heat;
                    conductance[i] = g;
                    if !g.is_finite() || g < 0.0 || !power.is_finite() || !heat.is_finite()
                        || !residual[i].is_finite() || (g != 0.0 && delta != 0.0 && power == 0.0)
                        || (power != 0.0 && heat == 0.0)
                    { return Err(ThermalPistonError::Invalid("invalid or unrepresentable paired interface heat law")); }
                }
                let stored = body.mass_kg * (endpoint.specific_enthalpy_j_kg() - h0);
                let body_residual = (stored + net) - body.internal_heat_j;
                let total = ((piston.gas[0].stored_energy_change_j + piston.gas[1].stored_energy_change_j)
                    + piston.kinetic_energy_change_j) + stored - input.heat_j[0] - input.heat_j[1]
                    - body.internal_heat_j - piston.external_force_work_j;
                if ![stored, body_residual, total].iter().all(|x| x.is_finite())
                    || body_residual.abs() > energy_tolerance || total.abs() > energy_tolerance
                { return Err(ThermalPistonError::Invalid("shared body and paired gas/momentum energy do not close")); }
                Ok((DoubleThermalPistonCandidate {
                    piston: next,
                    thermal: UniformEnthalpyStep { state: endpoint, external_heat_j: -net,
                        energy_residual_tolerance_j: energy_tolerance, report: samples.map(|s| s.source) },
                    step: DoubleThermalPistonStep {
                        piston, thermal_midpoint: midpoint, conductance_w_k: conductance,
                        exchanged_heat_j: q, source_residual_j: residual,
                        thermal_energy_change_j: stored, thermal_energy_residual_j: body_residual,
                        total_energy_residual_j: total, joint_trials, material_trials,
                    },
                }, residual[1]))
            }, cancelled)?;
            let residual = candidate.step.source_residual_j[0];
            Ok((candidate, residual))
        }, &mut cancelled)
    }
}

impl DoubleActingGasPiston {
    /// Source-bound series conduction on both faces, with no copied gas or body
    /// coefficients. Each descriptor must bind the same card as the thermal
    /// chart; geometry may differ between the faces. This reuses the existing
    /// layer owner, including missing-data, source-identity and phase checks.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_layers_candidate(
        &self, input: DoublePistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve, exchange: DoubleHeatExchangeSolve,
        layers: [super::layers::MaterialConductionLayers<'_>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleThermalPistonCandidate<fs_material::state_point::ResolvedMaterialStatePoint>,
        ThermalPistonError<super::layers::LayerConductionError>>
    {
        self.thermal_candidate(input, body, controls, exchange,
            |gases, thermal| Ok([
                layers[0].resolve(&gases[0], thermal)?,
                layers[1].resolve(&gases[1], thermal)?,
            ]), cancelled)
    }
}

#[cfg(test)]
mod tests;
