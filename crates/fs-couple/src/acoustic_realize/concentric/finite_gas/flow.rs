//! Species injection/withdrawal on the continuing layered thermal string.
//!
//! The original layer/heat/contact owner advances once. The native nozzle
//! transfers upstream species and enthalpy between its finite atmosphere and a
//! second parcel. Final gas loading and pressure use that conserved mixture,
//! without replacing the material layers or restarting their vibration.

use super::{
    ConcentricFiniteGasFrame, ConcentricFiniteGasSolve, ConcentricGasContact,
    ConcentricStringError, ConductanceSample, DepositionSimulation, EquilibriumPhaseState,
    OhmicDrive, ResolvedThermochemicalGasState, StepCoefficients, ThermalConcentricStringRuntime,
    checkpoint, finite, invalid,
};
use crate::acoustic_realize::finite_gas::flow::FiniteGasStringFlowDirection;
use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
};

/// One accepted motion/thermal/transport sample, with all four owners retained.
#[derive(Debug)]
pub struct ConcentricFiniteGasFlowFrame<S, R> {
    /// Original thermal substep accounts. Its gas source snapshot is the
    /// pre-transfer endpoint used by the heat law. The vibration's pressure
    /// observation uses the FINAL mixture and previous accepted acceleration.
    pub thermal: ConcentricFiniteGasFrame<S, R>,
    /// Upstream/downstream ownership of the original transfer report.
    pub direction: FiniteGasStringFlowDirection,
    pub flow: GasNozzleFlow,
    pub transfer: GasMassTransferStep,
    /// Actual surrounding gas change of m*u, including heat and mass flow [J].
    pub surrounding_energy_change_j: f64,
    /// Actual second parcel change of m*u [J].
    pub reservoir_energy_change_j: f64,
    /// Independent supplies minus ALL mechanical, layer and gas storage
    /// changes, authored loss and explicit damping-partition roundoff [J].
    /// Transferred enthalpy and interface/drag heat are internal transfers.
    pub energy_balance_residual_j: f64,
    pub energy_tolerance_j: f64,
}

impl ThermalConcentricStringRuntime<'_> {
    /// Advance unequal-temperature material layers, optional electrical drive,
    /// reciprocal surface heat and pressure-driven atmospheric species exchange.
    ///
    /// Every physical owner advances once. First the original layered heat
    /// transaction solves its material/thermal endpoint at fixed gas inventory;
    /// then the higher-pressure finite parcel supplies the existing nozzle for
    /// the same interval. Its source-law join conserves each species and total
    /// enthalpy, retaining formation-energy references and original transport
    /// data. Conflicting source models, over-withdrawal and pressure reversal
    /// retain the native refusal. No second jet-kinetic or flow-work heat is added.
    ///
    /// The final gas changes drag, sound speed and the pressure observation on
    /// the continuing coordinates. The old accepted acceleration is used even
    /// when the discarded thermal-stage observer has already been evaluated.
    /// One clock tick is consumed. All four owners publish together only after
    /// final source, temperature, observer, energy and cancellation admission.
    /// Callbacks retain the original repeatable-proposal contract.
    ///
    /// This is a first-order partitioned model: layer heat precedes gas flow.
    /// Both gas volumes are well mixed and fixed, with a quasi-steady converging
    /// nozzle. Refine physical time for the split and independent string modes
    /// for spatial accuracy. No resolved jet force, gas chamber/radiation load,
    /// propagation through an evolving medium or acoustic equivalence is claimed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_nozzle<'t, S, R, E>(
        &mut self,
        cx: &Cx<'_>,
        thermal: &mut DepositionSimulation<'t>,
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        controls: ConcentricFiniteGasSolve,
        energy_guess_j: &[f64],
        forces: &[f64],
        gas_heat_j: f64,
        drive: Option<OhmicDrive>,
        contact: ConcentricGasContact,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
        resolve_thermal: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        resolve_exchange: impl Fn(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<R>, E>,
    ) -> Result<ConcentricFiniteGasFlowFrame<S, R>, ConcentricStringError<E>> {
        checkpoint(cx)?;
        let mut candidate = self.clone();
        let mut next_thermal = thermal.clone();
        let mut heated_gas = gas.clone();
        let mut heated = candidate.step_with_isochoric_gas_exchange(
            cx,
            &mut next_thermal,
            &mut heated_gas,
            controls,
            energy_guess_j,
            forces,
            gas_heat_j,
            drive,
            contact,
            resolve_thermal,
            resolve_exchange,
        )?;
        let local = ClosedGasParcel::from_isochoric(heated_gas);
        let into = reservoir.state().state().pressure >= local.state().state().pressure;
        let (donor, receiver) = if into {
            (&*reservoir, &local)
        } else {
            (&local, &*reservoir)
        };
        let (proposal, flow) = donor.through_nozzle_to(
            receiver,
            nozzle,
            self.member.model.sample_period_s(),
            flow_controls,
            || cx.checkpoint().is_err(),
        )?;
        let (next_donor, next_receiver, transfer) = proposal.into_parts();
        let (next_gas, next_reservoir, direction) = if into {
            (
                next_receiver.into_isochoric(),
                next_donor,
                FiniteGasStringFlowDirection::IntoSurroundings,
            )
        } else {
            (
                next_donor.into_isochoric(),
                next_receiver,
                FiniteGasStringFlowDirection::OutOfSurroundings,
            )
        };
        if finite(next_gas.state().state().temperature - gas.state().state().temperature)?.abs()
            > controls.gas.maximum_gas_temperature_increment_k
        {
            return Err(invalid(
                "gas heat and species transport exceeded the per-sample temperature increment",
            ));
        }
        self.refresh_finite_gas_medium(
            cx,
            forces,
            &mut candidate,
            &mut heated.string.vibration,
            next_gas.state().state(),
        )?;

        // Extensive storage is essential: composition changes BOTH m and u.
        // Difference specific energies before multiplying to retain small heat
        // increments against a large valid formation-energy reference.
        let old_u = gas.state().specific_internal_energy_j_kg();
        let new_u = next_gas.state().specific_internal_energy_j_kg();
        let reservoir_u = reservoir.state().specific_internal_energy_j_kg();
        let next_reservoir_u = next_reservoir.state().specific_internal_energy_j_kg();
        let surrounding_energy_change_j = finite(
            next_gas.mass_kg() * (new_u - old_u) + (next_gas.mass_kg() - gas.mass_kg()) * old_u,
        )?;
        let reservoir_energy_change_j = finite(
            next_reservoir.mass_kg() * (next_reservoir_u - reservoir_u)
                + (next_reservoir.mass_kg() - reservoir.mass_kg()) * reservoir_u,
        )?;
        let mechanical_storage = finite(
            heated.string.mechanical_internal_energy_after_j
                - heated.string.mechanical_internal_energy_before_j,
        )?;
        let solid_storage = finite(next_thermal.energy_j() - thermal.energy_j())?;
        let input = finite(
            heated.string.vibration.acoustic.input_work_j
                + heated.external_solid_heat_j
                + gas_heat_j,
        )?;
        let loss = &heated.string.vibration.dissipation;
        let residual = finite(
            input
                - (mechanical_storage
                    + solid_storage
                    + surrounding_energy_change_j
                    + reservoir_energy_change_j
                    + loss.authored_loss_j
                    + loss.roundoff_residual_j),
        )?;
        let scale = finite(
            input.abs()
                + heated.string.mechanical_internal_energy_before_j.abs()
                + heated.string.mechanical_internal_energy_after_j.abs()
                + thermal.energy_j().abs()
                + next_thermal.energy_j().abs()
                + surrounding_energy_change_j.abs()
                + reservoir_energy_change_j.abs()
                + gas.mass_kg() * old_u.abs()
                + next_gas.mass_kg() * new_u.abs()
                + reservoir.mass_kg() * reservoir_u.abs()
                + next_reservoir.mass_kg() * next_reservoir_u.abs(),
        )?;
        let tolerance = finite(
            heated.energy_tolerance_j
                + flow_controls.transfer.energy_tolerance_j
                + 32.0 * f64::EPSILON * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(invalid(
                "layered mechanical, caloric and two-volume gas energy do not close",
            ));
        }
        checkpoint(cx)?;
        // Retain old owners until every endpoint has been published.
        let _old_runtime = core::mem::replace(self, candidate);
        let _old_thermal = core::mem::replace(thermal, next_thermal);
        let _old_gas = core::mem::replace(gas, next_gas);
        let _old_reservoir = core::mem::replace(reservoir, next_reservoir);
        Ok(ConcentricFiniteGasFlowFrame {
            thermal: heated,
            direction,
            flow,
            transfer,
            surrounding_energy_change_j,
            reservoir_energy_change_j,
            energy_balance_residual_j: residual,
            energy_tolerance_j: tolerance,
        })
    }
}
