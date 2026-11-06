//! One inertial wall between two entropy-compatible source-gas stores.
//!
//! Each original material supplies its discrete Gibbs pressure and temperature.
//! Both pressures act on the SAME displacement, so their work cancels the
//! corresponding gas energy debit. Each chamber retains its own heat/entropy
//! account; an adiabatic closed pair preserves both original source entropies.
//! This specializes the existing shared-wall geometry and impulse solve, without
//! inventing a second wall inertia, an exterior reservoir or copied caloric law.

use super::super::{
    ClosedGasParcel, PistonError, PistonGas, PistonGasState, PistonInput,
    ResolvedThermochemicalGasState, discrete_gradient::Trial as PrimaryTrial, impulse, invalid,
    poll, positive,
};
use super::{DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve};
use fs_material::thermochemical_gas::energy::parcel::discrete_gradient::{
    GasDiscreteGradientCandidate, GasDiscreteGradientStep,
};

/// Energy, entropy and momentum of the two original gases and their shared wall.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscreteGradientDoublePistonStep<S = ResolvedThermochemicalGasState> {
    /// Original gas heat/work/entropy results, [primary, opposing].
    pub gas: [GasDiscreteGradientStep; 2],
    /// Each original material at its entropy-conjugate temperature, used for
    /// source sampling and acoustics. The two work pressures are carried in
    /// `gas`; a nonideal material's sampling pressure need not equal its effort.
    pub thermodynamic_mean: [S; 2],
    /// Signed common wall displacement [m], expanding the primary chamber.
    pub displacement_m: f64,
    /// Initial and final velocity of the ONE wall [m/s].
    pub velocity_m_s: [f64; 2],
    /// Actual change in the ONE wall's kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Prescribed external force times actual displacement [J].
    pub external_force_work_j: f64,
    /// M*dv - dt*(A_primary*p_primary - A_opposing*p_opposing + F) [N s].
    pub momentum_residual_n_s: f64,
    /// dK - W_primary - W_opposing - W_external [J].
    pub mechanical_energy_residual_j: f64,
    /// dU_primary + dU_opposing + dK - Q_primary - Q_opposing - W_external [J].
    pub total_energy_residual_j: f64,
    /// dV_opposing + A_opposing*dx [m3].
    pub volume_exchange_residual_m3: f64,
    /// Interior shared momentum iterations.
    pub iterations: usize,
    /// Original material proposals, TWO per impulse trial.
    pub material_trials: usize,
}

struct PairTrial<G: PistonGas = ClosedGasParcel> {
    primary: PrimaryTrial<G>,
    opposing: GasDiscreteGradientCandidate<G, G::State>,
    primary_input: PistonInput,
    volume_residual_m3: f64,
}

impl<G: PistonGas> DoubleActingGasPiston<G> {
    fn discrete_gradient_trial(
        &self,
        volume_m3: f64,
        input: DoublePistonInput,
        controls: DoublePistonSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<impulse::Trial<PairTrial<G>>, PistonError> {
        poll(cancelled)?;
        let opposing_volume = self.opposing_volume(volume_m3)?;
        let opposing = self.opposing.discrete_gradient_heat_volume(
            opposing_volume,
            input.heat_j[1],
            controls.opposing_gas,
            &mut *cancelled,
        )?;
        let ratio = self.opposing_area_m2 / self.primary.area_m2();
        let exterior = opposing.step().work_pressure_pa * ratio;
        if !positive(exterior) {
            return Err(invalid(
                "opposing discrete-gradient pressure load is unrepresentable",
            ));
        }
        let primary_input = PistonInput {
            duration_s: input.duration_s,
            heat_j: input.heat_j[0],
            exterior_pressure_pa: exterior,
            external_force_n: input.external_force_n,
        };
        let primary = self.primary.discrete_gradient_trial(
            volume_m3,
            primary_input,
            controls.piston,
            cancelled,
        )?;
        // The opposing volume shrinks by A_opposing/A_primary times the
        // primary increment. Differentiate BOTH original pressure-work solves.
        let tangent = primary.tangent
            - input.duration_s
                * self.opposing_area_m2
                * ratio
                * opposing.step().pressure_final_volume_derivative_pa_m3;
        let volume_residual_m3 = (opposing.endpoint().volume_m3() - self.opposing.volume_m3())
            + self.opposing_area_m2 * primary.value.motion.displacement;
        if !positive(tangent)
            || !volume_residual_m3.is_finite()
            || volume_residual_m3.abs() > controls.volume_tolerance_m3
        {
            return Err(invalid(
                "two-gas discrete-gradient tangent or volume exchange failed admission",
            ));
        }
        Ok(impulse::Trial {
            residual: primary.residual,
            tangent,
            value: PairTrial {
                primary: primary.value,
                opposing,
                primary_input,
                volume_residual_m3,
            },
        })
    }

    /// Jointly evolve both original gas energies and the shared wall momentum.
    ///
    /// Each gas solves dU_i + p_bar_i*dV_i = Q_i and dS_i = Q_i/T_bar_i.
    /// The wall uses dV_primary=A_primary*dx, dV_opposing=-A_opposing*dx,
    /// dx=dt*(v0+v1)/2 and M*dv=dt*(A_primary*p_bar_primary
    /// - A_opposing*p_bar_opposing + F). Thus internal pressure work cancels
    /// exactly at the discrete equation level; accepted residuals remain explicit.
    ///
    /// Each chamber keeps its own species, caloric/transport data, mass,
    /// temperature support and energy budget. The coupled acoustic clock uses
    /// BOTH gas stiffnesses at the initial, thermodynamic-mean and final states.
    /// At most 2*(maximum_iterations+3) source-material proposals are made.
    ///
    /// Both chambers may use the original ideal owner or the original
    /// homogeneous Peng--Robinson owner, retaining its stable-branch admission.
    /// This does not change reacting, valve, circuit or acoustic-feedback rules.
    /// Supplied heat has a gas-side entropy contribution, but no unspecified
    /// reservoir entropy is inferred. Failure or cancellation preserves both
    /// physical states and every previously accepted port window.
    pub fn step_discrete_gradient(
        &mut self,
        input: DoublePistonInput,
        controls: DoublePistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiscreteGradientDoublePistonStep<G::State>, PistonError> {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.piston.volume_bracket_m3;
        if !positive(input.duration_s)
            || !input.heat_j.iter().all(|q| q.is_finite())
            || !input.external_force_n.is_finite()
            || !positive(controls.volume_tolerance_m3)
            || !positive(controls.piston.energy_tolerance_j)
            || !positive(lo)
            || !positive(hi)
            || lo >= hi
            || self.primary.gas().volume_m3() < lo
            || self.primary.gas().volume_m3() > hi
        {
            return Err(invalid(
                "invalid two-gas discrete-gradient input, bracket or residual controls",
            ));
        }
        self.opposing_volume(lo)?;
        self.opposing_volume(hi)?;
        self.admit_linear_clock(
            [self.primary.gas().state(), self.opposing.state()],
            [self.primary.gas().volume_m3(), self.opposing.volume_m3()],
            input.duration_s,
        )?;
        let solved = impulse::solve(self.primary.gas().volume_m3(), controls.piston, |volume| {
            self.discrete_gradient_trial(volume, input, controls, &mut cancelled)
        })?;
        let PairTrial {
            primary,
            opposing,
            primary_input,
            volume_residual_m3,
        } = solved.value;
        self.admit_linear_clock(
            [
                primary.gas.thermodynamic_mean(),
                opposing.thermodynamic_mean(),
            ],
            [
                self.primary.gas().mass_kg() / primary.gas.thermodynamic_mean().density_kg_m3(),
                self.opposing.mass_kg() / opposing.thermodynamic_mean().density_kg_m3(),
            ],
            input.duration_s,
        )?;
        self.admit_linear_clock(
            [primary.gas.endpoint().state(), opposing.endpoint().state()],
            [
                primary.gas.endpoint().volume_m3(),
                opposing.endpoint().volume_m3(),
            ],
            input.duration_s,
        )?;

        // Finalize the original one-sided owner on a PRIVATE candidate. Its
        // exterior-pressure work is the opposing gas's internal transfer.
        let opposing_mean = opposing.thermodynamic_mean().clone();
        let mut next_primary = self.primary.clone();
        let primary_step = next_primary.finish_discrete_gradient(
            primary,
            primary_input,
            controls.piston,
            solved.iterations,
            solved.proposals,
            &mut cancelled,
        )?;
        let (next_opposing, opposing_step) = opposing.into_parts();
        let force = self.primary.area_m2() * primary_step.gas.work_pressure_pa
            - self.opposing_area_m2 * opposing_step.work_pressure_pa
            + input.external_force_n;
        let momentum = self.primary.mass_kg()
            * (primary_step.velocity_m_s[1] - primary_step.velocity_m_s[0])
            - input.duration_s * force;
        let mechanical = (primary_step.kinetic_energy_change_j - primary_step.gas.work_by_gas_j)
            - opposing_step.work_by_gas_j
            - primary_step.external_force_work_j;
        let total = ((primary_step.gas.stored_energy_change_j
            + opposing_step.stored_energy_change_j)
            + primary_step.kinetic_energy_change_j)
            - input.heat_j[0]
            - input.heat_j[1]
            - primary_step.external_force_work_j;
        if ![force, momentum, mechanical, total]
            .iter()
            .all(|x| x.is_finite())
            || momentum.abs() > controls.piston.momentum_tolerance_n_s
            || mechanical.abs() > controls.piston.energy_tolerance_j
            || total.abs() > controls.piston.energy_tolerance_j
        {
            return Err(invalid(
                "two-gas discrete-gradient momentum or total energy did not close",
            ));
        }
        let report = DiscreteGradientDoublePistonStep {
            gas: [primary_step.gas, opposing_step],
            thermodynamic_mean: [primary_step.thermodynamic_mean, opposing_mean],
            displacement_m: primary_step.displacement_m,
            velocity_m_s: primary_step.velocity_m_s,
            kinetic_energy_change_j: primary_step.kinetic_energy_change_j,
            external_force_work_j: primary_step.external_force_work_j,
            momentum_residual_n_s: momentum,
            mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total,
            volume_exchange_residual_m3: volume_residual_m3,
            iterations: solved.iterations,
            material_trials: 2 * solved.proposals,
        };
        poll(&mut cancelled)?;
        let forces = [
            self.primary.area_m2() * report.gas[0].work_pressure_pa,
            -self.opposing_area_m2 * report.gas[1].work_pressure_pa,
        ];
        self.ports.accept(
            crate::port_adoption::StepAccount {
                tick: self.ports.next_tick(),
                duration_s: input.duration_s,
                storage_change_j: Some(
                    report.gas[0].stored_energy_change_j
                        + report.gas[1].stored_energy_change_j
                        + report.kinetic_energy_change_j,
                ),
                loss_j: Some(0.0),
                included_work_j: report.external_force_work_j,
                entropy_change_j_k: report.gas[0].entropy_change_j_k
                    + report.gas[1].entropy_change_j_k,
                tolerance_j: Some(controls.piston.energy_tolerance_j),
                unavailable: None,
            },
            [
                crate::port_adoption::PortTransfer::work(forces[0] * report.displacement_m, false),
                crate::port_adoption::PortTransfer::work(forces[1] * report.displacement_m, false),
                crate::port_adoption::PortTransfer::work(
                    -(forces[0] + forces[1]) * report.displacement_m,
                    false,
                ),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j[0],
                    report.thermodynamic_mean[0].temperature_k(),
                ),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j[1],
                    report.thermodynamic_mean[1].temperature_k(),
                ),
            ],
        );
        self.primary = next_primary;
        self.opposing = next_opposing;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
