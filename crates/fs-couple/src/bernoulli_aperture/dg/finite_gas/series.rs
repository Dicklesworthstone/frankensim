//! One explicitly shared, well-mixed finite background for both original tracts.
//!
//! Both inlet gases already belong to the same volume-flow valve model. One
//! original series sample advances its mechanics and acoustic owners; one
//! scalar thermochemical solve then pays their SUM of native parameter work.
//! This does not replace two independently evolving inlet inventories.

use super::super::tract::prepare_upstream_step;
use super::super::{ApertureDg, DgApertureFrame};
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::bernoulli_aperture::wall::material::{
    admit_finite_gas_controls, finite_gas_field, heated_finite_gas,
};
use crate::cavity_phs::piston::thermal::feedback;
use crate::dg_port::finite_gas::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgSolve, acoustic, admit_background,
    admitted_volume, check, feedback_error, refresh_kind, refresh_output,
};
use crate::dg_port::{DgPortGasUpdate, DgPortLoad, sum_terminal_work};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

/// One accepted gas, source and two-tract interval. Gas storage and heat appear
/// once; the two refresh reports are ordered upstream, downstream.
#[derive(Debug)]
pub struct FiniteGasSeriesApertureDgStep<S> {
    /// Original interval forces/flows/work, with refreshed endpoint duct reports.
    pub coupled: DgApertureFrame,
    pub refresh: [DgPortGasUpdate; 2],
    /// Original caloric/transport result and any still-private companion states.
    pub background: S,
    pub supplied_heat_j: f64,
    pub transported_enthalpy_j: f64,
    pub gas_energy_change_j: f64,
    /// One debit OUT of the shared gas, paying both complete native transfers.
    pub work_by_gas_j: f64,
    pub work_residual_j: f64,
    pub gas_energy_residual_j: f64,
    /// Gas + original mechanics + both complete loads + losses/exterior work,
    /// minus heat, transported enthalpy and actual external drives [J].
    pub combined_energy_residual_j: f64,
    pub trials: usize,
}

struct Candidate<S> {
    loads: [DgPortLoad; 2],
    background: FiniteGasDgBackground<S>,
    refresh: [DgPortGasUpdate; 2],
    gas_change: f64,
    gas_residual: f64,
    work: f64,
    work_residual: f64,
    trials: usize,
}

impl ApertureDg<LoadedAperture> {
    /// Advance the original two-sided source with ONE shared finite gas parcel.
    ///
    /// The parcel must match both original complete uniform source states and
    /// its volume must equal BOTH main bores plus their main-gas finite terminal
    /// passages. Independently sourced branches retain their own inventories.
    /// This is an explicitly well-mixed shared background, not a mixing model
    /// or a substitute for two independent finite gas reservoirs.
    ///
    /// The callback receives the unchanged original parcel, dt, TOTAL work OUT
    /// of gas and cancellation. Every work trial begins with that same parcel
    /// and the same already-marched acoustic histories. Each original native
    /// transfer retains momentum and memory; one allowance bounds the sum of
    /// ABSOLUTE bore/channel work across both sides. `controls.loss` selects
    /// the native gas-memory transfer on each tract without dropping histories.
    ///
    /// Held-force/flow integrals retain their actual frozen-medium interval;
    /// endpoint pressures, storage, density, viscosity and summed exact source
    /// impedance use the accepted gas. The whole gas/mechanical/acoustic account
    /// is checked before publishing any participant. Return companion candidates
    /// through `source` and publish them only after success.
    ///
    /// Fixed geometry and the first-order partition remain explicit. No mean
    /// pressure prestress, thermodynamic mass flow inferred from acoustic flow,
    /// or thermal deposition of mechanical/acoustic loss is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_finite_gas_controlled<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        controls: FiniteGasDgSolve,
        mut evolve: impl FnMut(
            &ClosedGasParcel,
            f64,
            f64,
            &mut C,
        ) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasSeriesApertureDgStep<S>, FiniteGasDgError<E>> {
        check(cancelled)?;
        if self.supply_channel.is_some() {
            return Err(FiniteGasDgError::Invalid(
                "shared finite gas requires the fixed supply passage's explicit inventory and parameter work",
            ));
        }
        admit_finite_gas_controls(controls)?;
        let upstream = self.upstream.as_ref().ok_or(FiniteGasDgError::Invalid(
            "shared finite gas requires the original upstream and downstream tracts",
        ))?;
        let original = [upstream, &self.load];
        if original[0].state().step_count() != original[1].state().step_count()
            || original[0].state().step_count() != self.aperture.accepted_steps()
            || original[0].propagator().dt().to_bits() != original[1].propagator().dt().to_bits()
        {
            return Err(FiniteGasDgError::Invalid(
                "shared finite gas requires both tracts and original mechanics on one clock",
            ));
        }
        let volume = admitted_volume(original[0], gas, cancelled)?
            + admitted_volume(original[1], gas, cancelled)?;
        if !volume.is_finite()
            || volume <= 0.0
            || (gas.volume_m3() / volume - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(FiniteGasDgError::Invalid(
                "shared finite inventory volume must equal both bores and their main-gas terminal passages",
            ));
        }
        let initial_energy = self.stored_energy_j();
        let dt = original[0].propagator().dt();
        let fields = original.map(|load| load.duct().gas_field().expect("source admitted"));
        let mut marched = [original[0].clone(), original[1].clone()];
        let [up, down] = &mut marched;
        let (mechanical, mut frame) =
            prepare_upstream_step(&self.aperture, up, down, drive, forces_n, pressures_pa)
                .map_err(acoustic)?;
        check(cancelled)?;
        let [lo, hi] = controls.work_bracket_j;
        let maximum_work = lo.abs().max(hi.abs());
        let refresh = refresh_kind(controls.loss);
        let mut candidate = feedback::solve_scalar(
            controls.work_bracket_j,
            controls.work_tolerance_j,
            controls.maximum_iterations,
            |work, trials, cancel| {
                let background =
                    evolve(gas, dt, work, cancel).map_err(FiniteGasDgError::Background)?;
                check(cancel)?;
                let (gas_change, gas_residual) =
                    admit_background(gas, &background, work, controls)?;
                let mut loads = marched.clone();
                let up_field = finite_gas_field(fields[0], &background.gas)?;
                let up = loads[0]
                    .rebind_gas_field_controlled(&up_field, refresh, maximum_work, cancel)
                    .map_err(acoustic)?;
                let down_field = finite_gas_field(fields[1], &background.gas)?;
                let down = loads[1]
                    .rebind_gas_field_controlled(&down_field, refresh, maximum_work, cancel)
                    .map_err(acoustic)?;
                let actual_work = up.parameter_work_j + down.parameter_work_j;
                let absolute_work = up.absolute_parameter_work_j + down.absolute_parameter_work_j;
                if !actual_work.is_finite()
                    || !absolute_work.is_finite()
                    || absolute_work > maximum_work
                {
                    return Err(FiniteGasDgError::Invalid(
                        "shared gas intervention exceeds total absolute two-tract work allowance",
                    ));
                }
                let work_residual = work - actual_work;
                Ok((
                    Candidate {
                        loads,
                        background,
                        refresh: [up, down],
                        gas_change,
                        gas_residual,
                        work,
                        work_residual,
                        trials,
                    },
                    work_residual,
                ))
            },
            cancelled,
        )
        .map_err(feedback_error)?;

        let tract = frame.upstream.as_mut().expect("original series frame");
        refresh_output(&mut candidate.loads[0], &mut tract.duct, cancelled)?;
        refresh_output(&mut candidate.loads[1], &mut frame.duct, cancelled)?;
        let terminal_work =
            sum_terminal_work(&tract.terminals) + sum_terminal_work(&frame.terminals);
        frame.parameter_work_j =
            candidate.refresh[0].parameter_work_j + candidate.refresh[1].parameter_work_j;
        frame.stored_energy_j = frame.aperture.stored_energy_j
            + candidate.loads[0].stored_energy_j()
            + candidate.loads[1].stored_energy_j();
        frame.storage_change_j = frame.stored_energy_j - initial_energy;
        let combined = ((candidate.gas_change
            + frame.storage_change_j
            + frame.inferred_total_dissipation_j
            + terminal_work)
            - frame.upstream_work_j
            - frame.body_work_j
            - frame.aperture.mechanical_work_j
            - candidate.background.supplied_heat_j)
            - candidate.background.transported_enthalpy_j;
        let impedance =
            candidate.loads[0].impedance_pa_s_m3() + candidate.loads[1].impedance_pa_s_m3();
        if !frame.stored_energy_j.is_finite()
            || !frame.storage_change_j.is_finite()
            || !frame.balance_residual_j().is_finite()
            || frame.balance_residual_j().abs() > controls.energy_tolerance_j
            || !combined.is_finite()
            || combined.abs() > controls.energy_tolerance_j
            || candidate.loads[0].duct().inlet_gas() != candidate.loads[1].duct().inlet_gas()
            || !impedance.is_finite()
            || impedance <= 0.0
        {
            return Err(FiniteGasDgError::Invalid(
                "shared finite gas, original mechanics and both complete tracts do not close their common endpoint",
            ));
        }
        check(cancelled)?;
        // Source validation precedes its first mutation; publication is then
        // infallible for mechanics, both complete loads and the one shared gas.
        self.aperture
            .rebind_exact_gas_load(candidate.loads[1].duct().inlet_gas(), impedance)
            .map_err(acoustic)?;
        self.aperture.accept(mechanical);
        let [up, down] = candidate.loads;
        self.upstream = Some(up);
        self.load = down;
        *gas = candidate.background.gas;
        Ok(FiniteGasSeriesApertureDgStep {
            coupled: frame,
            refresh: candidate.refresh,
            background: candidate.background.source,
            supplied_heat_j: candidate.background.supplied_heat_j,
            transported_enthalpy_j: candidate.background.transported_enthalpy_j,
            gas_energy_change_j: candidate.gas_change,
            work_by_gas_j: candidate.work,
            work_residual_j: candidate.work_residual,
            gas_energy_residual_j: candidate.gas_residual,
            combined_energy_residual_j: combined,
            trials: candidate.trials,
        })
    }

    /// Heat one original shared finite inventory, paying both native acoustic
    /// work accounts once through the existing caloric owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_finite_gas_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        cancelled: &mut C,
    ) -> Result<FiniteGasSeriesApertureDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>>
    {
        self.step_with_shared_finite_gas_controlled(
            gas,
            drive,
            forces_n,
            pressures_pa,
            controls,
            |old, _, work, cancel| heated_finite_gas(old, heat_j, work, gas_solve, cancel),
            cancelled,
        )
    }
}
