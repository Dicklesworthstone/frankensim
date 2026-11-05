//! Original suspended mechanics with an explicit physical heat destination.

use super::*;

impl MixtureEuler1d {
    /// Couple both chamber pressures, a suspended piston and one finite wall.
    ///
    /// The original spring and damper participate in every joint midpoint-
    /// velocity proposal. After accepting that mechanical interval, the
    /// authored `damping_heat_fraction_to_wall` in [0, 1] deposits its share of
    /// the actual nonnegative `c*w_mid^2*dt` into the ONE original wall store.
    /// The remainder leaves this represented thermal system. No material or
    /// heat partition is inferred from a damping coefficient.
    ///
    /// The composition is H(dt/2), A_with_suspension(dt), deposit, H(dt/2).
    /// The deposit uses the original finite wall's calorics, never its frozen
    /// capacity tangent. The second H uses the heated wall and moved gas-cell
    /// volumes, so dissipation changes physical gas pressure and subsequent
    /// motion. Heat generated in this A does not retrospectively change its
    /// already accepted force; this is a split thermomechanical step, not a
    /// simultaneous temperature-dependent suspension solve.
    ///
    /// The enclosing suspended total accounts for actual gas and wall storage,
    /// piston kinetic and spring energy, EXTERNAL damping remainder, and drive
    /// work. Its `damping_energy_j` still reports the full mechanical loss;
    /// `piston.gas.heat` separates the returned heat and external remainder.
    /// The nested piston account keeps spring/damper force work and subtracts
    /// the returned heat once. Its pressure impulses/work never include heat.
    ///
    /// All source, pressure, energy and cancellation admissions precede joint
    /// publication. Shortened intervals and failed complete attempts restore
    /// both gases, the piston and the wall before retrying; no trial loss is
    /// deposited twice. Native pressure and caloric workspaces are reused.
    /// The spring/anchor stay fixed during the step. No chemistry, molecular
    /// transport, wall gradients, temperature-dependent damping, entropy
    /// inequality or general second-order/FSI-stability claim is added.
    ///
    /// # Errors
    /// Invalid fraction, original model/source/geometry refusal, exhausted work
    /// or retries, failed balances and cancellation preserve accepted state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_suspended_piston_and_heat(
        &mut self,
        right: &mut MixtureEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        suspension: PistonSuspension,
        damping_heat_fraction_to_wall: f64,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SuspendedPistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
        thermal.advance(
            self,
            right,
            piston,
            maximum_duration,
            gas_controls,
            heat_controls,
            external_force_n,
            Some(suspension),
            damping_heat_fraction_to_wall,
            piston_controls,
            None,
            cancelled,
        )
    }

    /// Preserve independently admitted material contacts during suspended heat.
    ///
    /// This uses the same physical composition and accounts as
    /// [`Self::step_with_shared_suspended_piston_and_heat`]. The two contact
    /// workspaces bind the actual fields after H1, and each enforces its own
    /// numerical-defect budget during every pressure proposal. Physical heat
    /// and resolved pressure gradients are not equilibrium corrections. The
    /// full energy account subtracts only the separately reported signed
    /// contact defects, never spring storage or returned damper heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_suspended_piston_heat_and_contact_equilibrium(
        &mut self,
        right: &mut MixtureEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        suspension: PistonSuspension,
        damping_heat_fraction_to_wall: f64,
        piston_controls: PistonControl,
        equilibrium: [ContactEquilibriumControl; 2],
        workspaces: [&mut ContactWorkspace; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SuspendedPistonStep<
            SharedPistonThermalStepReport<ContactEquilibriumStepReport<StepReport>>,
        >,
        SharedPistonThermalError,
    > {
        let [a, b] = workspaces;
        let step = thermal.advance(
            self,
            right,
            piston,
            maximum_duration,
            gas_controls,
            heat_controls,
            external_force_n,
            Some(suspension),
            damping_heat_fraction_to_wall,
            piston_controls,
            Some((equilibrium, [&mut *a, &mut *b])),
            cancelled,
        )?;
        Ok(step.map_gas(|gas| SharedPistonThermalStepReport {
            transport: contact_reports(gas.transport, [a.receipt(), b.receipt()]),
            before: gas.before,
            after: gas.after,
            heat: gas.heat,
            restarts: gas.restarts,
        }))
    }
}
