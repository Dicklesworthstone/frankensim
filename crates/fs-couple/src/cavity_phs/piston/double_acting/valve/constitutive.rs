//! Original material ownership for metered transfer and subsequent closed motion.
//!
//! The sealed piston gas supplies its existing nozzle, transfer report and
//! caloric increment. This seam adds no flow law, EOS or source interpretation.

use crate::cavity_phs::piston::PistonGas;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
};
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;
use fs_material::thermochemical_gas::nonideal::energy::parcel::transfer::RealGasMassTransferStep;
use fs_material::thermochemical_gas::nonideal::energy::parcel::transfer::nozzle::{
    RealGasNozzleFlow, RealGasNozzleStepSolve,
};
use fs_qty::Temperature;

/// Existing nozzle and conservation accounts for the sealed piston materials.
///
/// Implementations remain limited by [`PistonGas`] to the original ideal and
/// Peng--Robinson parcel owners. Their source, phase, inventory, work-budget and
/// cancellation admissions remain authoritative for every transfer proposal.
pub trait MeteredPistonGas: PistonGas {
    /// Original nozzle, endpoint-energy and inventory controls.
    type NozzleSolve: Copy + core::fmt::Debug + PartialEq;
    /// Original source-bearing throat state and flow account.
    type Flow: Clone + core::fmt::Debug + PartialEq;
    /// Original donor/receiver species and extensive-energy account.
    type Transfer: Clone + core::fmt::Debug + PartialEq;

    /// Accepted conserved inventory mass [kg].
    fn mass_kg(&self) -> f64;

    /// Propose both finite endpoints through the original material nozzle.
    ///
    /// Returns donor, receiver, transfer and throat flow in that order. Inputs
    /// remain unchanged; the enclosing valve publishes both endpoints together
    /// only after its subsequent motion, heat and conservation gates succeed.
    fn through_metering_nozzle(
        &self,
        receiver: &Self,
        nozzle: GasNozzle,
        duration_s: f64,
        controls: Self::NozzleSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, Self, Self::Transfer, Self::Flow), GasEnergyError>;

    /// Nonnegative mass moved from donor to receiver [kg].
    fn transferred_mass_kg(transfer: &Self::Transfer) -> f64;

    /// Original extensive storage changes [J], ordered [donor, receiver].
    /// These displayed quantities retain the material's formation-energy datum.
    fn transfer_energy_change_j(transfer: &Self::Transfer) -> [f64; 2];

    /// Original whole-pair energy defect [J], with transported datums cancelled
    /// by the material owner before rounding. Do not reconstruct from totals.
    fn transfer_energy_residual_j(transfer: &Self::Transfer) -> f64;

    /// Original absolute inventory-conservation tolerance [kg].
    fn mass_tolerance_kg(controls: Self::NozzleSolve) -> f64;

    /// Closed-motion storage change [J] from an explicit post-transfer state.
    ///
    /// The accepted final mass multiplies the original model's specific-energy
    /// increment. Composition and mass must remain unchanged over this closed
    /// motion. Ideal-gas energy depends only on temperature; real-gas energy
    /// also retains the supplied reference density and actual final density.
    fn closed_energy_change_j(
        &self,
        reference_density_kg_m3: f64,
        reference_temperature_k: f64,
    ) -> Result<f64, GasEnergyError>;
}

impl MeteredPistonGas for ClosedGasParcel {
    type NozzleSolve = GasNozzleStepSolve;
    type Flow = GasNozzleFlow;
    type Transfer = GasMassTransferStep;

    fn mass_kg(&self) -> f64 {
        ClosedGasParcel::mass_kg(self)
    }

    fn through_metering_nozzle(
        &self,
        receiver: &Self,
        nozzle: GasNozzle,
        duration_s: f64,
        controls: Self::NozzleSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, Self, Self::Transfer, Self::Flow), GasEnergyError> {
        let (candidate, flow) =
            self.through_nozzle_to(receiver, nozzle, duration_s, controls, cancelled)?;
        let (donor, receiver, transfer) = candidate.into_parts();
        Ok((donor, receiver, transfer, flow))
    }

    fn transferred_mass_kg(transfer: &Self::Transfer) -> f64 {
        transfer.transferred_mass_kg
    }

    fn transfer_energy_change_j(transfer: &Self::Transfer) -> [f64; 2] {
        [
            transfer.donor_energy_change_j,
            transfer.receiver_energy_change_j,
        ]
    }

    fn transfer_energy_residual_j(transfer: &Self::Transfer) -> f64 {
        transfer.energy_residual_j
    }

    fn mass_tolerance_kg(controls: Self::NozzleSolve) -> f64 {
        controls.transfer.mass_tolerance_kg
    }

    fn closed_energy_change_j(
        &self,
        _reference_density_kg_m3: f64,
        reference_temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        let change = self
            .model()
            .thermodynamics()
            .specific_internal_energy_change(
                Temperature::new(reference_temperature_k),
                Temperature::new(self.state().state().temperature),
            )
            .map_err(|error| {
                GasEnergyError::Material(fs_material::MaterialError::State {
                    what: format!("metering-valve endpoint caloric increment: {error}"),
                })
            })?
            .value();
        Ok(self.mass_kg() * change)
    }
}

impl MeteredPistonGas for IsochoricRealGasParcel {
    type NozzleSolve = RealGasNozzleStepSolve;
    type Flow = RealGasNozzleFlow;
    type Transfer = RealGasMassTransferStep;

    fn mass_kg(&self) -> f64 {
        IsochoricRealGasParcel::mass_kg(self)
    }

    fn through_metering_nozzle(
        &self,
        receiver: &Self,
        nozzle: GasNozzle,
        duration_s: f64,
        controls: Self::NozzleSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, Self, Self::Transfer, Self::Flow), GasEnergyError> {
        let (candidate, flow) =
            self.through_nozzle_to(receiver, nozzle, duration_s, controls, cancelled)?;
        let (donor, receiver, transfer) = candidate.into_parts();
        Ok((donor, receiver, transfer, flow))
    }

    fn transferred_mass_kg(transfer: &Self::Transfer) -> f64 {
        transfer.transferred_mass_kg
    }

    fn transfer_energy_change_j(transfer: &Self::Transfer) -> [f64; 2] {
        [
            transfer.donor_energy_change_j,
            transfer.receiver_energy_change_j,
        ]
    }

    fn transfer_energy_residual_j(transfer: &Self::Transfer) -> f64 {
        transfer.energy_residual_j
    }

    fn mass_tolerance_kg(controls: Self::NozzleSolve) -> f64 {
        controls.transfer.mass_tolerance_kg
    }

    fn closed_energy_change_j(
        &self,
        reference_density_kg_m3: f64,
        reference_temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        Ok(self.mass_kg()
            * self
                .model()
                .specific_internal_energy_change_between_densities(
                    reference_density_kg_m3,
                    reference_temperature_k,
                    self.density_kg_m3(),
                    self.state().temperature_k(),
                )?)
    }
}
