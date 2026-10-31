//! Finite-reservoir gas exchange coupled to the existing inertial wall.
//!
//! This is an explicit first-order transfer-then-motion split over ONE physical
//! interval. The nozzle owner moves species and upstream total enthalpy at fixed
//! volumes; the original midpoint owner then advances the wall and gas energy.
//! It is not a simultaneous compressible-flow solve or a jet-momentum model.

/// Nozzle exchange and reciprocal material-wall heat transfer on one accepted step.
pub mod thermal;
/// Species exchange, electrical storage, finite thermal cells and wall motion.
pub mod network;
/// Shared-state multiport flow driving independent moving gas chambers.
pub mod graph;
/// Source-backed Peng--Robinson reservoir flow driving the same piston mechanics.
pub mod real_gas;

use super::{ClosedGasParcel, GasEnergyError, GasPiston, PistonError, PistonInput,
    PistonSolve, PistonStep, invalid, poll, positive};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
};

/// Upstream side selected from the two actual accepted pressures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NozzleDirection {
    /// The finite reservoir supplies the piston cavity.
    IntoPiston,
    /// The piston cavity supplies the finite reservoir.
    OutOfPiston,
}

/// Separate stage accounts plus independent balances enclosing BOTH gases.
/// Internal transferred enthalpy cancels within this whole-system boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct PistonNozzleStep {
    /// None when the admitted path transfers no mass.
    pub direction: Option<NozzleDirection>,
    /// Original source-resolved upstream/throat flow, with choking information.
    pub flow: GasNozzleFlow,
    /// Fixed-volume species/enthalpy transfer in its original donor order.
    pub transfer: GasMassTransferStep,
    /// Original midpoint mechanics from the post-transfer chamber state.
    pub piston: PistonStep,
    /// Full chamber internal-energy change, including transferred enthalpy [J].
    pub chamber_energy_change_j: f64,
    /// Full finite-reservoir internal-energy change [J].
    pub reservoir_energy_change_j: f64,
    /// Independently measured total mass discrepancy across both endpoints [kg].
    pub total_mass_residual_kg: f64,
    /// dU_chamber+dU_reservoir+dK+outside_work-force_work-prescribed_heat [J].
    pub total_energy_residual_j: f64,
}

// The two consumers share transfer, endpoint accounting and open-boundary
// observation. Each retains its original mechanical/thermal solve owner.
struct NozzleTransfer {
    piston: GasPiston,
    reservoir: ClosedGasParcel,
    direction: Option<NozzleDirection>,
    flow: GasNozzleFlow,
    transfer: GasMassTransferStep,
}

struct NozzleBalance {
    chamber_energy_j: f64,
    reservoir_energy_j: f64,
    mass_residual_kg: f64,
    kinetic_energy_j: f64,
    exterior_work_j: f64,
    force_work_j: f64,
    residual_j: f64,
}

fn stored_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    // Stable extensive difference, retaining each source's energy reference.
    next.mass_kg() * (next.state().specific_internal_energy_j_kg()
        - old.state().specific_internal_energy_j_kg())
        + (next.mass_kg() - old.mass_kg()) * old.state().specific_internal_energy_j_kg()
}

impl GasPiston {
    fn nozzle_transfer(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle, duration_s: f64,
        exchange: GasNozzleStepSolve, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NozzleTransfer, PistonError> {
        let inflow = reservoir.state().state().pressure >= self.gas.state().state().pressure;
        let (candidate, flow) = if inflow {
            reservoir.through_nozzle_to(&self.gas, nozzle, duration_s, exchange, &mut *cancelled)?
        } else {
            self.gas.through_nozzle_to(reservoir, nozzle, duration_s, exchange, &mut *cancelled)?
        };
        let (donor, receiver, transfer) = candidate.into_parts();
        let (reservoir, chamber) = if inflow { (donor, receiver) } else { (receiver, donor) };
        poll(cancelled)?;
        let mut piston = self.clone();
        piston.gas = chamber;
        let direction = (transfer.transferred_mass_kg != 0.0).then_some(
            if inflow { NozzleDirection::IntoPiston } else { NozzleDirection::OutOfPiston });
        Ok(NozzleTransfer { piston, reservoir, direction, flow, transfer })
    }

    fn nozzle_balance(
        &self, reservoir: &ClosedGasParcel, next: &NozzleTransfer, input: PistonInput,
        exchange: GasNozzleStepSolve, tolerance_j: f64,
    ) -> Result<NozzleBalance, PistonError> {
        let chamber_energy_j = stored_change(&self.gas, &next.piston.gas);
        let reservoir_energy_j = stored_change(reservoir, &next.reservoir);
        let mass_residual_kg = (next.piston.gas.mass_kg() - self.gas.mass_kg())
            + (next.reservoir.mass_kg() - reservoir.mass_kg());
        let kinetic_energy_j = (0.5 * self.mass_kg * (next.piston.velocity_m_s - self.velocity_m_s))
            * (next.piston.velocity_m_s + self.velocity_m_s);
        let volume_change = next.piston.gas.volume_m3() - self.gas.volume_m3();
        let exterior_work_j = input.exterior_pressure_pa * volume_change;
        let force_work_j = input.external_force_n * (volume_change / self.area_m2);
        let residual_j = ((chamber_energy_j + reservoir_energy_j) + kinetic_energy_j)
            + exterior_work_j - force_work_j - input.heat_j;
        if ![chamber_energy_j, reservoir_energy_j, mass_residual_kg, kinetic_energy_j,
            exterior_work_j, force_work_j, residual_j].iter().all(|value| value.is_finite())
            || mass_residual_kg.abs() > exchange.transfer.mass_tolerance_kg
            || residual_j.abs() > tolerance_j {
            return Err(invalid("nozzle/piston accepted endpoints do not close whole-pair mass and energy"));
        }
        Ok(NozzleBalance { chamber_energy_j, reservoir_energy_j, mass_residual_kg,
            kinetic_energy_j, exterior_work_j, force_work_j, residual_j })
    }

    fn record_nozzle_account(
        &self, next: &mut NozzleTransfer, input: PistonInput, midpoint_temperature_k: f64,
        balance: &NozzleBalance, tolerance_j: f64,
    ) -> Result<(), PistonError> {
        if next.direction.is_some() {
            let initial_entropy = self.gas.state().specific_entropy_j_kg_k().map_err(GasEnergyError::from)?;
            let final_entropy = next.piston.gas.state().specific_entropy_j_kg_k().map_err(GasEnergyError::from)?;
            let entropy_change = next.piston.gas.mass_kg() * (final_entropy - initial_entropy)
                + (next.piston.gas.mass_kg() - self.gas.mass_kg()) * initial_entropy;
            if !entropy_change.is_finite() {
                return Err(invalid("nozzle/piston chamber entropy change is unrepresentable"));
            }
            next.piston.ports.accept(crate::port_adoption::StepAccount {
                tick: self.ports.next_tick(), duration_s: input.duration_s,
                storage_change_j: Some(balance.chamber_energy_j + balance.kinetic_energy_j), loss_j: Some(0.0),
                included_work_j: 0.0, entropy_change_j_k: entropy_change,
                tolerance_j: Some(tolerance_j),
                unavailable: Some("nozzle mass exchange needs a stream port; use the nozzle step's paired balance"),
            }, [crate::port_adoption::PortTransfer::work(balance.force_work_j - balance.exterior_work_j, true),
                crate::port_adoption::PortTransfer::heat(input.heat_j, midpoint_temperature_k)]);
        }
        Ok(())
    }

    /// Exchange through a nozzle and advance the existing gas/wall equations.
    ///
    /// Flow direction uses the accepted initial pressures. Transfer uses the
    /// original source-backed species, caloric, EOS and transport owners; no
    /// nominal gas constant, averaged composition or incoming-heat copy is used.
    /// `input.heat_j` is external heat ONLY: incoming total enthalpy already enters
    /// the transfer stage, including thermalized jet energy and flow work.
    ///
    /// The two stages each approximate their part of the same `duration_s`.
    /// This first-order splitting needs timestep refinement. The nozzle's donor
    /// fraction, source and pressure-crossing bounds remain active; piston clock,
    /// bracket and energy checks also remain active. Geometry has constant area,
    /// the reservoir is a well-mixed fixed volume, and no jet thrust is invented.
    ///
    /// Both parcels, wall state and prior port history publish together. With a
    /// nonzero mass transfer the existing two-port wall/heat audit is marked
    /// unavailable: it lacks the stream port needed for the chamber's open mass
    /// boundary. Its real wall/heat transfers remain visible. Use this report's
    /// separate stage accounts and independent whole-pair balance instead.
    /// Zero transfer preserves the ordinary closed-step port account.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_nozzle(
        &mut self, reservoir: &mut ClosedGasParcel, nozzle: GasNozzle,
        input: PistonInput, exchange: GasNozzleStepSolve, controls: PistonSolve,
        total_energy_tolerance_j: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonNozzleStep, PistonError> {
        poll(&mut cancelled)?;
        if !positive(total_energy_tolerance_j) {
            return Err(invalid("positive finite nozzle/piston total energy tolerance required"));
        }
        let mut next = self.nozzle_transfer(reservoir, nozzle, input.duration_s, exchange, &mut cancelled)?;
        let piston = next.piston.step(input, controls, &mut cancelled)?;
        let balance = self.nozzle_balance(reservoir, &next, input, exchange, total_energy_tolerance_j)?;
        self.record_nozzle_account(&mut next, input, piston.midpoint.state().temperature,
            &balance, total_energy_tolerance_j)?;
        let report = PistonNozzleStep { direction: next.direction, flow: next.flow,
            transfer: next.transfer, piston,
            chamber_energy_change_j: balance.chamber_energy_j,
            reservoir_energy_change_j: balance.reservoir_energy_j,
            total_mass_residual_kg: balance.mass_residual_kg, total_energy_residual_j: balance.residual_j };
        poll(&mut cancelled)?;
        *reservoir = next.reservoir;
        *self = next.piston;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;

