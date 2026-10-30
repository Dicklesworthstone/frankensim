//! Two thermal ports on ONE two-gas/mechanical proposal, without extra stores.
//! Finite counterparts must debit the signed heat and publish all endpoints.

/// Existing finite circuit/enthalpy network on both gas ports.
pub mod network;

use super::{ConductanceSample, ResolvedThermochemicalGasState, ThermalPistonError,
    feedback, heat, poll, positive};
use super::super::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::double_acting::{DoubleActingGasPiston,
    DoublePistonInput, DoublePistonSolve, DoublePistonStep};

/// Thermal boundary and shared mechanical accounts for one accepted proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleBoundaryStep {
    /// Both source-resolved gases, ONE displacement and ONE kinetic account.
    pub piston: DoublePistonStep,
    /// Counterpart temperatures for this trial [K], [primary, opposing].
    pub boundary_temperature_k: [f64; 2],
    /// Actual midpoint conductances [W/K], in the same face order.
    pub conductance_w_k: [f64; 2],
    /// Signed applied heats INTO the gases [J]. Debit each counterpart once.
    pub exchanged_heat_j: [f64; 2],
    /// Independently evaluated dt*G*(T_boundary-T_gas_mid) [J].
    pub evaluated_heat_j: [f64; 2],
    /// Applied minus evaluated heat [J], checked SEPARATELY for each gas.
    pub source_residual_j: [f64; 2],
    /// Both gas internal changes + ONE kinetic change - external force work
    /// - prescribed gas heats - evaluated boundary heats [J].
    pub total_energy_residual_j: f64,
    /// All complete paired proposals across both bounded source searches.
    pub joint_trials: usize,
    /// All gas-material proposals, including discarded source/impulse trials.
    pub material_trials: usize,
}

/// Both gases, common momentum, and both interface evaluations accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleBoundaryCandidate<S> {
    piston: DoubleActingGasPiston,
    sources: [S; 2],
    step: DoubleBoundaryStep,
}
impl<S> DoubleBoundaryCandidate<S> {
    /// Tentative shared mechanism. Neither accepted gas was advanced in trials.
    pub const fn piston(&self) -> &DoubleActingGasPiston { &self.piston }
    /// Original interface evaluations at the SAME coupled proposal.
    pub const fn sources(&self) -> &[S; 2] { &self.sources }
    /// Independent source and combined-energy admission.
    pub const fn step(&self) -> &DoubleBoundaryStep { &self.step }
    /// Publish with every finite counterpart; the heat is already in the gases.
    pub fn into_parts(self) -> (DoubleActingGasPiston, [S; 2], DoubleBoundaryStep) {
        (self.piston, self.sources, self.step)
    }
}

impl DoubleActingGasPiston {
    /// Solve both thermal boundaries and the shared gas/momentum equations.
    ///
    /// Every [Qa,Qb] trial advances the ORIGINAL paired piston, never two
    /// independent walls. The existing caloric/EOS/transport owners resolve both
    /// gas midpoints; both conductance laws then see that SAME physical proposal.
    /// The existing bounded scalar feedback owner solves both heat residuals.
    /// A total proposal limit includes all rejected inner and outer trials.
    ///
    /// Temperatures are held for reservoirs. A finite network must supply its
    /// current nonlinear-trial temperatures, debit the returned applied heats
    /// from the respective cells, and admit all physical endpoints together.
    /// This API creates no fictitious finite body or infinite heat-capacity law.
    /// The repeatable callback must not modify external accepted physical state.
    ///
    /// Source residuals are admitted separately even when Qa+Qb is zero.
    /// Isolation still calls both source providers: missing data is not zero G.
    /// Failure/cancellation returns no partial candidate. Arbitrary conductance
    /// laws do not imply uniqueness, monotonicity or unconditional stability;
    /// the original uniform, closed-gas and rigid-wall assumptions still apply.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_boundaries_candidate<S, E>(
        &self, input: DoublePistonInput, controls: DoublePistonSolve,
        exchange: DoubleHeatExchangeSolve, boundary_temperature_k: [f64; 2],
        resolve: impl FnMut(&[ResolvedThermochemicalGasState; 2])
            -> Result<[ConductanceSample<S>; 2], E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleBoundaryCandidate<S>, ThermalPistonError<E>> {
        self.thermal_boundaries_candidate_with_midpoint_resistance(input, controls,
            exchange, boundary_temperature_k, 0.0, resolve, cancelled)
    }

    /// Shared composition seam: every thermal trial sees the same affine
    /// mechanical load, while the original impulse owner supplies its tangent.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::cavity_phs::piston) fn thermal_boundaries_candidate_with_midpoint_resistance<S, E>(
        &self, input: DoublePistonInput, controls: DoublePistonSolve,
        exchange: DoubleHeatExchangeSolve, boundary_temperature_k: [f64; 2],
        resistance_n_s_m: f64,
        mut resolve: impl FnMut(&[ResolvedThermochemicalGasState; 2])
            -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleBoundaryCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        for face in exchange.faces { feedback::admit::<E>(face)?; }
        if !positive(input.duration_s) || !boundary_temperature_k.iter().all(|t| positive(*t))
            || !input.heat_j.iter().all(|q| q.is_finite()) || exchange.maximum_joint_trials == 0
            || !resistance_n_s_m.is_finite() || resistance_n_s_m < 0.0 {
            return Err(ThermalPistonError::Invalid("finite paired boundary temperatures, heat, clock and work bound required"));
        }
        let tolerance = exchange.faces[0].energy_tolerance_j.min(exchange.faces[1].energy_tolerance_j);
        let mut joint_trials = 0usize;
        let mut material_trials = 0usize;
        let candidate = feedback::solve(exchange.faces[0], |qa, _, cancelled| {
            let candidate = feedback::solve(exchange.faces[1], |qb, _, cancelled| {
                if joint_trials >= exchange.maximum_joint_trials {
                    return Err(ThermalPistonError::Invalid("paired boundary joint trial budget exhausted"));
                }
                joint_trials += 1;
                let q = [qa, qb];
                let heat_j = [input.heat_j[0] + qa, input.heat_j[1] + qb];
                if !heat_j.iter().all(|q| q.is_finite())
                    || (0..2).any(|i| q[i] != 0.0 && heat_j[i] == input.heat_j[i]) {
                    return Err(ThermalPistonError::Invalid("boundary heat vanishes or overflows in paired gas input"));
                }
                let mut next = self.clone();
                let piston = next.step_with_midpoint_resistance(DoublePistonInput { heat_j, ..input },
                    controls, resistance_n_s_m, &mut *cancelled)?;
                material_trials = material_trials.checked_add(piston.material_trials)
                    .ok_or(ThermalPistonError::Invalid("paired material proposal count overflow"))?;
                poll(cancelled)?;
                let samples = resolve(&piston.midpoint).map_err(ThermalPistonError::Interface)?;
                let conductance = [samples[0].conductance_w_k, samples[1].conductance_w_k];
                let mut evaluated = [0.0; 2];
                for i in 0..2 {
                    evaluated[i] = heat::conductive_heat(input.duration_s, conductance[i],
                        boundary_temperature_k[i], piston.midpoint[i].state().temperature)?;
                }
                let residual = [qa-evaluated[0], qb-evaluated[1]];
                let total = (piston.gas[0].stored_energy_change_j + piston.gas[1].stored_energy_change_j)
                    + piston.kinetic_energy_change_j - piston.external_force_work_j
                    - input.heat_j[0] - input.heat_j[1] - evaluated[0] - evaluated[1];
                if !residual.iter().all(|r| r.is_finite()) || !total.is_finite() {
                    return Err(ThermalPistonError::Invalid("unrepresentable paired thermal boundary balance"));
                }
                Ok((DoubleBoundaryCandidate { piston: next, sources: samples.map(|s| s.source),
                    step: DoubleBoundaryStep { piston, boundary_temperature_k, conductance_w_k: conductance,
                        exchanged_heat_j: q, evaluated_heat_j: evaluated, source_residual_j: residual,
                        total_energy_residual_j: total, joint_trials, material_trials } }, residual[1]))
            }, cancelled)?;
            let residual = candidate.step.source_residual_j[0];
            Ok((candidate, residual))
        }, &mut cancelled)?;
        if candidate.step.total_energy_residual_j.abs() > tolerance {
            return Err(ThermalPistonError::Invalid("both evaluated boundary heats must close shared gas/mechanical energy"));
        }
        poll(&mut cancelled)?;
        Ok(candidate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cavity_phs::piston::{GasPiston, tests::{controls, gas}};
    use crate::cavity_phs::piston::thermal::HeatExchangeSolve;
    use std::convert::Infallible;

    fn piston() -> DoubleActingGasPiston {
        DoubleActingGasPiston::new(GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap(),
            gas(0.75,0.001),0.01).unwrap()
    }
    fn drive() -> DoublePistonInput {
        // About 118 rad/s: 0.00005 s keeps omega*dt inside the 0.04-rad band.
        DoublePistonInput { duration_s: 0.00005, heat_j: [0.0;2], external_force_n: 0.0 }
    }
    fn solve() -> DoublePistonSolve {
        DoublePistonSolve { piston: controls(), opposing_gas: controls().gas, volume_tolerance_m3: 1e-15 }
    }
    fn exchange() -> DoubleHeatExchangeSolve {
        DoubleHeatExchangeSolve { faces: [HeatExchangeSolve { heat_bracket_j: [-1.0,1.0],
            source_tolerance_j: 1e-10, energy_tolerance_j: 1e-7, maximum_iterations: 128 };2],
            maximum_joint_trials: 20000 }
    }
    fn links(g: [f64;2]) -> impl FnMut(&[ResolvedThermochemicalGasState;2])
        -> Result<[ConductanceSample<f64>;2],Infallible> {
        move |gas| Ok([0,1].map(|i| ConductanceSample { conductance_w_k: g[i],
            source: gas[i].state().temperature }))
    }
    #[test]
    fn equal_reservoirs_heat_both_gases_without_inventing_relative_motion() {
        let old = piston();
        let c = old.thermal_boundaries_candidate(drive(),solve(),exchange(),[500.0;2],links([1.0;2]),||false).unwrap();
        for i in 0..2 {
            assert!(c.step().exchanged_heat_j[i] > 0.0);
            assert!(c.step().source_residual_j[i].abs() <= exchange().faces[i].source_tolerance_j);
            assert_eq!(c.sources()[i],c.step().piston.midpoint[i].state().temperature);
        }
        assert!(c.piston().primary().velocity_m_s().abs() < 1e-8);
        assert!(c.step().total_energy_residual_j.abs() < 1e-7);
        assert_eq!(old,piston());
    }
    #[test]
    fn unequal_boundary_heat_changes_both_volumes_and_keeps_one_energy_balance() {
        let old = piston();
        let c = old.thermal_boundaries_candidate(drive(),solve(),exchange(),[600.0,300.0],links([1.0,2.0]),||false).unwrap();
        assert!(c.step().exchanged_heat_j[0] > 0.0 && c.step().exchanged_heat_j[1] < 0.0);
        assert!(c.piston().primary().velocity_m_s() > 0.0);
        assert!(c.piston().primary().gas().volume_m3() > old.primary().gas().volume_m3());
        assert!(c.piston().opposing().volume_m3() < old.opposing().volume_m3());
        let applied = c.step().exchanged_heat_j;
        let mut independent = old.clone();
        let expected = independent.step(DoublePistonInput { heat_j: applied, ..drive() },solve(),||false).unwrap();
        assert_eq!(c.piston(),&independent);
        assert_eq!(c.step().piston,expected);
        assert!(c.step().total_energy_residual_j.abs() < 1e-7);
    }
    #[test]
    fn insulated_boundaries_reproduce_original_mechanics_bitwise() {
        let old = piston();
        let input = DoublePistonInput { heat_j: [0.2,0.1], external_force_n: -0.1, ..drive() };
        let c = old.thermal_boundaries_candidate(input,solve(),exchange(),[500.0,300.0],links([0.0;2]),||false).unwrap();
        let mut direct = old.clone();
        let expected = direct.step(input,solve(),||false).unwrap();
        assert_eq!(c.piston(),&direct); assert_eq!(c.step().piston,expected);
        assert_eq!(c.step().joint_trials,1);
        assert_eq!(c.step().exchanged_heat_j,[0.0;2]);
    }
    #[test]
    fn source_refusal_budget_and_final_cancellation_return_neither_gas() {
        let old = piston(); let before = old.clone();
        let refused = |_: &[ResolvedThermochemicalGasState;2]| -> Result<[ConductanceSample<()>;2],&str> { Err("missing interface") };
        assert!(matches!(old.thermal_boundaries_candidate(drive(),solve(),exchange(),[500.0;2],
            refused,||false),Err(ThermalPistonError::Interface("missing interface"))));
        assert!(old.thermal_boundaries_candidate(drive(),solve(),DoubleHeatExchangeSolve {
            maximum_joint_trials: 1, ..exchange() },[500.0;2],links([1.0;2]),||false).is_err());
        assert!(old.thermal_boundaries_candidate(drive(),solve(),exchange(),[0.0,500.0],links([1.0;2]),||false).is_err());
        assert!(old.thermal_boundaries_candidate(drive(),solve(),exchange(),[500.0;2],links([1.0,-1.0]),||false).is_err());
        let mut polls = 0;
        old.thermal_boundaries_candidate(drive(),solve(),exchange(),[500.0;2],links([1.0;2]),
            || { polls+=1; false }).unwrap();
        let mut count = 0;
        assert!(matches!(old.thermal_boundaries_candidate(drive(),solve(),exchange(),[500.0;2],links([1.0;2]),
            || { count+=1; count==polls }),Err(ThermalPistonError::Cancelled)));
        assert_eq!(old,before);
    }
}
