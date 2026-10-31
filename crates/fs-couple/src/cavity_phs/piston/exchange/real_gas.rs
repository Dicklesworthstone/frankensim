//! Finite-reservoir Peng--Robinson nozzle exchange on the native piston step.
//!
//! The source-backed real-gas nozzle first transfers mass and upstream enthalpy
//! between fixed well-mixed volumes. The resulting chamber then advances through
//! the existing nonideal midpoint piston mechanics. Stream enthalpy and boundary
//! work therefore have distinct owners and are each counted exactly once.

use super::super::{
    GasPiston, PistonError, PistonInput, PistonSolve, PistonStep, invalid, poll, positive,
};
use super::NozzleDirection;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;
use fs_material::thermochemical_gas::nonideal::energy::parcel::transfer::{
    RealGasMassTransferStep,
    nozzle::{RealGasNozzleFlow, RealGasNozzleStepSolve},
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;

/// One accepted real-gas stream transfer followed by native pressure/work motion.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasPistonNozzleStep {
    /// None when the pressure path is shut and no mass moved.
    pub direction: Option<NozzleDirection>,
    /// Source-backed real-gas throat, choking and mass rate.
    pub flow: RealGasNozzleFlow,
    /// Fixed-volume conservative mass/enthalpy transfer.
    pub transfer: RealGasMassTransferStep,
    /// Existing nonideal midpoint piston mechanics from the transferred chamber.
    pub piston: PistonStep<ResolvedThermochemicalRealGasState>,
    /// Complete chamber internal-energy change over transfer plus piston stages [J].
    pub chamber_energy_change_j: f64,
    /// Complete finite-reservoir internal-energy change [J].
    pub reservoir_energy_change_j: f64,
    /// Net mass defect across the two gas volumes [kg].
    pub total_mass_residual_kg: f64,
    /// Transfer residual plus native piston total-energy residual [J].
    pub total_energy_residual_j: f64,
}

struct Transfer {
    piston: GasPiston<IsochoricRealGasParcel>,
    reservoir: IsochoricRealGasParcel,
    direction: Option<NozzleDirection>,
    flow: RealGasNozzleFlow,
    transfer: RealGasMassTransferStep,
}

impl GasPiston<IsochoricRealGasParcel> {
    fn real_nozzle_transfer(
        &self,
        reservoir: &IsochoricRealGasParcel,
        nozzle: GasNozzle,
        duration_s: f64,
        exchange: RealGasNozzleStepSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Transfer, PistonError> {
        let inflow = reservoir.state().pressure_pa() >= self.gas.state().pressure_pa();
        let (candidate, flow) = if inflow {
            reservoir.through_nozzle_to(
                &self.gas, nozzle, duration_s, exchange, &mut *cancelled,
            )?
        } else {
            self.gas.through_nozzle_to(
                reservoir, nozzle, duration_s, exchange, &mut *cancelled,
            )?
        };
        let (donor, receiver, transfer) = candidate.into_parts();
        let (reservoir, chamber) = if inflow {
            (donor, receiver)
        } else {
            (receiver, donor)
        };
        poll(cancelled)?;
        let mut piston = self.clone();
        piston.gas = chamber;
        let direction = (transfer.transferred_mass_kg != 0.0).then_some(
            if inflow { NozzleDirection::IntoPiston } else { NozzleDirection::OutOfPiston },
        );
        Ok(Transfer { piston, reservoir, direction, flow, transfer })
    }

    /// Exchange with a finite real-gas reservoir, then advance the same chamber
    /// through the native nonideal piston pressure/work solve.
    ///
    /// The flow stage uses the original PR isentrope and donor enthalpy. The
    /// mechanics stage starts from that transferred chamber and owns only
    /// subsequent prescribed heat and moving-boundary work. The pair residual is
    /// the sum of those independently reference-safe ledgers; no absolute energy
    /// subtraction or ideal-gas correction is used.
    ///
    /// This remains a first-order transfer-then-motion split. The nozzle is a
    /// quasi-steady converging restriction and both gas volumes are homogeneous.
    /// No jet thrust, moving reservoir wall, phase equilibrium or composition
    /// change is inferred. Refine physical time for split accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_real_gas_nozzle(
        &mut self,
        reservoir: &mut IsochoricRealGasParcel,
        nozzle: GasNozzle,
        input: PistonInput,
        exchange: RealGasNozzleStepSolve,
        controls: PistonSolve,
        total_energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasPistonNozzleStep, PistonError> {
        poll(&mut cancelled)?;
        if !positive(total_energy_tolerance_j) {
            return Err(invalid(
                "positive finite real-gas nozzle/piston total energy tolerance required",
            ));
        }
        let mut next = self.real_nozzle_transfer(
            reservoir, nozzle, input.duration_s, exchange, &mut cancelled,
        )?;
        let piston = next.piston.step(input, controls, &mut cancelled)?;
        let (transfer_chamber, transfer_reservoir) = match next.direction {
            Some(NozzleDirection::IntoPiston) => (
                next.transfer.receiver_energy_change_j,
                next.transfer.donor_energy_change_j,
            ),
            Some(NozzleDirection::OutOfPiston) => (
                next.transfer.donor_energy_change_j,
                next.transfer.receiver_energy_change_j,
            ),
            None => (0.0, 0.0),
        };
        let chamber_energy_change_j =
            transfer_chamber + piston.gas.stored_energy_change_j;
        let reservoir_energy_change_j = transfer_reservoir;
        let total_mass_residual_kg = next.transfer.mass_residual_kg;
        let total_energy_residual_j =
            next.transfer.energy_residual_j + piston.total_energy_residual_j;
        if ![
            chamber_energy_change_j,
            reservoir_energy_change_j,
            total_mass_residual_kg,
            total_energy_residual_j,
        ].into_iter().all(f64::is_finite)
            || total_mass_residual_kg.abs() > exchange.transfer.mass_tolerance_kg
            || total_energy_residual_j.abs() > total_energy_tolerance_j
        {
            return Err(invalid(
                "real-gas nozzle/piston endpoints do not close whole-pair mass and energy",
            ));
        }

        // The existing wall/heat schema has no stream port. Preserve its actual
        // mechanical and heat transfers but explicitly mark the open-mass audit
        // unavailable rather than presenting a closed-chamber certificate.
        if next.direction.is_some() {
            next.piston.ports.accept(
                crate::port_adoption::StepAccount {
                    tick: self.ports.next_tick(),
                    duration_s: input.duration_s,
                    storage_change_j: Some(
                        chamber_energy_change_j + piston.kinetic_energy_change_j,
                    ),
                    loss_j: Some(0.0),
                    included_work_j: 0.0,
                    entropy_change_j_k: piston.gas.entropy_change_j_k,
                    tolerance_j: Some(total_energy_tolerance_j),
                    unavailable: Some(
                        "real-gas nozzle mass exchange needs a stream port; use the paired nozzle balance",
                    ),
                },
                [
                    crate::port_adoption::PortTransfer::work(
                        piston.external_force_work_j
                            - piston.exterior_pressure_work_j
                            + piston.support_work_j,
                        true,
                    ),
                    crate::port_adoption::PortTransfer::heat(
                        input.heat_j,
                        piston.midpoint.temperature_k(),
                    ),
                ],
            );
        }

        let report = RealGasPistonNozzleStep {
            direction: next.direction,
            flow: next.flow,
            transfer: next.transfer,
            piston,
            chamber_energy_change_j,
            reservoir_energy_change_j,
            total_mass_residual_kg,
            total_energy_residual_j,
        };
        poll(&mut cancelled)?;
        *reservoir = next.reservoir;
        *self = next.piston;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
