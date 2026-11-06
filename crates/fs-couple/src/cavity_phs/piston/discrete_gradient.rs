//! Entropy-compatible ideal-gas work on the original inertial piston.
//!
//! Gas calorics and the discrete Gibbs identity belong to the material parcel.
//! This adapter uses that SAME pressure in the original wall kinematics,
//! safeguarded impulse solve and port account. It retains both the material
//! entropy residual and the independently measured mechanical energy residual.
//! The implicit midpoint momentum rule remains second order. Entropy balance
//! does not establish spatial accuracy or resolve a physical heat reservoir.

use super::{
    ClosedGasParcel, GasPiston, MechanicalTrial, PistonError, PistonInput,
    PistonSolve, PistonSupport, ResolvedThermochemicalGasState,
    impulse, invalid, poll, positive,
};
use fs_material::thermochemical_gas::energy::parcel::discrete_gradient::{
    GasDiscreteGradientCandidate, GasDiscreteGradientStep,
};

/// Accepted material and wall state for the discrete-gradient gas work rule.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscreteGradientPistonStep {
    /// Original material energy, entropy, Gibbs balance and pressure derivatives.
    pub gas: GasDiscreteGradientStep,
    /// Source-bearing thermodynamic mean used in pressure work and the heat port.
    /// It is not an arithmetic midpoint of temperatures or volumes.
    pub thermodynamic_mean: ResolvedThermochemicalGasState,
    /// Signed wall displacement [m], positive towards gas expansion.
    pub displacement_m: f64,
    /// Initial and final wall velocities [m/s].
    pub velocity_m_s: [f64; 2],
    /// Actual change of wall kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Additional prescribed force times actual displacement [J].
    pub external_force_work_j: f64,
    /// Work against prescribed exterior pressure [J].
    pub exterior_pressure_work_j: f64,
    /// Impulse exerted on the wall by its fixed support [N s].
    pub support_impulse_n_s: f64,
    /// Support impulse times the actual midpoint wall velocity [J].
    pub support_work_j: f64,
    /// Independently reconstructed wall impulse residual [N s].
    pub momentum_residual_n_s: f64,
    /// Wall kinetic change minus all mechanical input work [J].
    pub mechanical_energy_residual_j: f64,
    /// Gas plus wall energy change minus supplied heat and external work [J].
    pub total_energy_residual_j: f64,
    /// Interior mechanical iterations.
    pub iterations: usize,
    /// Original source-material proposals, including bracket endpoints.
    pub material_trials: usize,
}

struct Trial {
    gas: GasDiscreteGradientCandidate,
    motion: MechanicalTrial,
}

impl GasPiston<ClosedGasParcel> {
    fn discrete_gradient_trial(
        &self,
        volume_m3: f64,
        input: PistonInput,
        controls: PistonSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<impulse::Trial<Trial>, PistonError> {
        poll(cancelled)?;
        let gas = self.gas.discrete_gradient_heat_volume(
            volume_m3, input.heat_j, controls.gas, &mut *cancelled,
        )?;
        let motion = self.mechanical_trial(
            volume_m3, gas.thermodynamic_mean().state().pressure, input,
        )?;
        let tangent = (2.0 * (self.mass_kg / self.area_m2)) / input.duration_s
            - (input.duration_s * self.area_m2)
                * gas.step().pressure_final_volume_derivative_pa_m3;
        if self.support == PistonSupport::Free && !positive(tangent) {
            return Err(invalid("discrete-gradient piston impulse tangent is unrepresentable"));
        }
        Ok(impulse::Trial {
            residual: motion.residual,
            tangent,
            value: Trial { gas, motion },
        })
    }

    /// Advance the original wall using source-law discrete-gradient gas work.
    ///
    /// The material owner solves dU + p_bar*dV = Q and
    /// dS = Q/T_bar within its declared residual tolerance. This adapter solves
    /// dV = A*dt*(v0+v1)/2 and M*dv = dt*(A*(p_bar-p_ext)+F)
    /// with that SAME p_bar and its implicit final-volume derivative.
    /// Adiabatic steps therefore avoid midpoint gas entropy drift while closing
    /// gas and wall energy together. The heat port uses the conjugate T_bar.
    ///
    /// Frozen-composition ideal mixtures only. This does not extend the rule
    /// to a real-gas EOS or assign entropy to an unspecified external reservoir.
    /// Existing midpoint, reacting, nozzle and finite-wall thermal entry points
    /// retain their own integration rules; this choice applies to this call.
    /// Fixed support preserves the original volume and zero wall velocity.
    /// Failure or cancellation leaves gas, motion and the accepted port window
    /// unchanged. The original acoustic clock and nonlinear budgets still apply.
    pub fn step_discrete_gradient(
        &mut self,
        input: PistonInput,
        controls: PistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiscreteGradientPistonStep, PistonError> {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.volume_bracket_m3;
        if !positive(input.duration_s)
            || !input.heat_j.is_finite()
            || !input.exterior_pressure_pa.is_finite()
            || input.exterior_pressure_pa < 0.0
            || !input.external_force_n.is_finite()
            || !positive(lo)
            || !positive(hi)
            || lo > hi
            || (lo == hi && self.support == PistonSupport::Free)
            || self.gas.volume_m3() < lo
            || self.gas.volume_m3() > hi
            || !positive(controls.momentum_tolerance_n_s)
            || !positive(controls.energy_tolerance_j)
            || !(1..=128).contains(&controls.maximum_iterations)
        {
            return Err(invalid("invalid discrete-gradient piston input, bracket or residual controls"));
        }
        self.admit_linear_clock(self.gas.state(), self.gas.volume_m3(), input.duration_s)?;
        let (trial, iterations, material_trials) = if self.support == PistonSupport::Fixed {
            let trial = self.discrete_gradient_trial(
                self.gas.volume_m3(), input, controls, &mut cancelled,
            )?;
            (trial.value, 0, 1)
        } else {
            let solved = impulse::solve(self.gas.volume_m3(), controls, |volume| {
                self.discrete_gradient_trial(volume, input, controls, &mut cancelled)
            })?;
            (solved.value, solved.iterations, solved.proposals)
        };

        let mean = trial.gas.thermodynamic_mean();
        let mean_volume = self.gas.mass_kg() / mean.state().density;
        self.admit_linear_clock(mean, mean_volume, input.duration_s)?;
        self.admit_linear_clock(
            trial.gas.endpoint().state(), trial.gas.endpoint().volume_m3(), input.duration_s,
        )?;
        let before = self.velocity_m_s;
        let after = trial.motion.velocity;
        let kinetic = (0.5 * self.mass_kg * (after - before)) * (after + before);
        let exterior = input.exterior_pressure_pa
            * (trial.gas.endpoint().volume_m3() - self.gas.volume_m3());
        let external = input.external_force_n * trial.motion.displacement;
        let support = trial.motion.support_impulse * f64::midpoint(before, after);
        let force = self.area_m2 * (mean.state().pressure - input.exterior_pressure_pa)
            + input.external_force_n;
        let momentum = self.mass_kg * (after - before)
            - input.duration_s * force - trial.motion.support_impulse;
        let gas_step = trial.gas.step();
        let mechanical =
            (kinetic - gas_step.work_by_gas_j) + exterior - external - support;
        let total = (gas_step.stored_energy_change_j + kinetic)
            + exterior - external - support - input.heat_j;
        if ![kinetic, exterior, external, support, momentum, mechanical, total]
            .iter().all(|value| value.is_finite())
            || momentum.abs() > controls.momentum_tolerance_n_s
            || mechanical.abs() > controls.energy_tolerance_j
            || total.abs() > controls.energy_tolerance_j
        {
            return Err(invalid("discrete-gradient piston impulse or total energy did not close"));
        }
        let thermodynamic_mean = mean.clone();
        let (gas, gas_step) = trial.gas.into_parts();
        let report = DiscreteGradientPistonStep {
            gas: gas_step,
            thermodynamic_mean,
            displacement_m: trial.motion.displacement,
            velocity_m_s: [before, after],
            kinetic_energy_change_j: kinetic,
            external_force_work_j: external,
            exterior_pressure_work_j: exterior,
            support_impulse_n_s: trial.motion.support_impulse,
            support_work_j: support,
            momentum_residual_n_s: momentum,
            mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total,
            iterations,
            material_trials,
        };
        poll(&mut cancelled)?;
        self.ports.accept(
            crate::port_adoption::StepAccount {
                tick: self.ports.next_tick(),
                duration_s: input.duration_s,
                storage_change_j: Some(report.gas.stored_energy_change_j + kinetic),
                loss_j: Some(0.0),
                included_work_j: 0.0,
                entropy_change_j_k: report.gas.entropy_change_j_k,
                tolerance_j: Some(controls.energy_tolerance_j),
                unavailable: None,
            },
            [
                crate::port_adoption::PortTransfer::work(external - exterior + support, true),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j, report.thermodynamic_mean.state().temperature,
                ),
            ],
        );
        self.gas = gas;
        self.velocity_m_s = after;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
