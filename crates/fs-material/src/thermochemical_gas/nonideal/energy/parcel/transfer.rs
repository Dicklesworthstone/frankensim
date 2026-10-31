//! Conservative upstream-enthalpy transfer between two fixed real-gas volumes.
//!
//! Both parcels must retain the exact same frozen composition, NASA calorics,
//! Peng--Robinson interaction model and transport sources. A prescribed positive
//! mass leaves the donor carrying its upstream total enthalpy and enters the
//! receiver once. Each fixed vessel then reconstructs temperature and pressure
//! from its new density and conserved internal energy through the original PR
//! owner. No ideal-gas pressure, gamma or caloric surrogate is introduced.
//!
//! This is a first-order well-mixed finite-volume transfer primitive. It does
//! not choose a flow rate, resolve a nozzle/jet, evolve momentum, perform phase
//! equilibrium or change composition. A higher-level pressure/port owner must
//! choose and refine the transferred mass.

/// Quasi-steady pressure-driven flow using the same source-backed PR branch.
pub mod nozzle;

use super::IsochoricRealGasParcel;
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, poll, updated_point,
};
use crate::thermochemical_gas::positive;

/// Bounded endpoint solves and absolute extensive conservation tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasMassTransferSolve {
    /// Donor caloric bracket and specific-energy tolerance.
    pub donor: GasEnergySolve,
    /// Receiver caloric bracket and specific-energy tolerance.
    pub receiver: GasEnergySolve,
    /// Absolute mass-balance/representability tolerance [kg].
    pub mass_tolerance_kg: f64,
    /// Absolute pair internal-energy tolerance [J].
    pub energy_tolerance_j: f64,
}
impl RealGasMassTransferSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.donor.validate()?;
        self.receiver.validate()?;
        if !positive(self.mass_tolerance_kg) || !positive(self.energy_tolerance_j) {
            return Err(GasEnergyError::Invalid(
                "positive finite real-gas transfer mass and energy tolerances required",
            ));
        }
        Ok(())
    }
}

/// Conserved two-volume accounting for one prescribed upstream transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasMassTransferStep {
    /// Positive mass moved from donor to receiver [kg].
    pub transferred_mass_kg: f64,
    /// Donor upstream specific enthalpy on the original NASA/PR datum [J/kg].
    pub source_enthalpy_j_kg: f64,
    /// Signed enthalpy carried into the receiver [J].
    pub transferred_enthalpy_j: f64,
    /// Donor extensive internal-energy change [J].
    pub donor_energy_change_j: f64,
    /// Receiver extensive internal-energy change [J].
    pub receiver_energy_change_j: f64,
    /// Pair energy defect after the transported enthalpy cancels [J].
    pub energy_residual_j: f64,
    /// Net mass defect across both endpoints [kg].
    pub mass_residual_kg: f64,
    /// Sum of interior iterations in both original energy inversions.
    pub iterations: usize,
}

/// Immutable tentative endpoints; callers publish both or neither.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasMassTransferCandidate {
    donor: IsochoricRealGasParcel,
    receiver: IsochoricRealGasParcel,
    step: RealGasMassTransferStep,
}
impl RealGasMassTransferCandidate {
    #[must_use]
    pub const fn donor(&self) -> &IsochoricRealGasParcel { &self.donor }
    #[must_use]
    pub const fn receiver(&self) -> &IsochoricRealGasParcel { &self.receiver }
    #[must_use]
    pub const fn step(&self) -> &RealGasMassTransferStep { &self.step }
    #[must_use]
    pub fn into_parts(self) -> (
        IsochoricRealGasParcel, IsochoricRealGasParcel, RealGasMassTransferStep,
    ) {
        (self.donor, self.receiver, self.step)
    }
}

fn endpoint(
    old: &IsochoricRealGasParcel,
    mass_kg: f64,
    change_j_kg: f64,
    controls: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(IsochoricRealGasParcel, f64, usize), GasEnergyError> {
    let density = mass_kg / old.volume_m3;
    if !positive(density) || (mass_kg != old.mass_kg && density == old.density_kg_m3) {
        return Err(GasEnergyError::Invalid(
            "real-gas transferred endpoint density is unrepresentable",
        ));
    }
    let solved = old.model.reconstruct_density_internal_energy_increment_between_densities(
        &old.point,
        old.density_kg_m3,
        density,
        change_j_kg,
        controls,
        &mut *cancelled,
    )?;
    let residual = solved.residual_j_kg();
    let iterations = solved.iterations();
    let state = solved.into_state();
    let point = updated_point(&old.point, state.temperature_k(), state.pressure_pa())?;
    let mut next = old.clone();
    next.point = point;
    next.state = state;
    next.density_kg_m3 = density;
    next.mass_kg = mass_kg;
    Ok((next, residual, iterations))
}

impl IsochoricRealGasParcel {
    /// Transfer a prescribed positive mass to another fixed-volume real-gas parcel.
    ///
    /// The two parcels must use identical frozen composition and PR/material
    /// source models. Therefore mass transfer changes density and energy but not
    /// species fractions; composition-changing real-gas mixing requires a
    /// separate EOS/source owner and is deliberately refused here.
    ///
    /// The donor loses `dm*h_d`; the receiver gains the same quantity. The
    /// energy solves are written in reference-safe differences:
    /// donor `Delta u = -(dm/m1) p_d/rho_d`, while the receiver uses
    /// `(dm/m1)[u_d-u_r+p_d/rho_d]`, with `u_d-u_r` evaluated directly by
    /// the original NASA/PR increment law. Large formation-energy constants
    /// therefore do not contaminate the nonlinear target.
    ///
    /// Zero mass is exact identity. Nonzero transfer never mutates either input;
    /// source refusal, an unsupported homogeneous branch, cancellation or failed
    /// conservation returns no physical endpoint.
    pub fn transfer_mass_to(
        &self,
        receiver: &Self,
        mass_kg: f64,
        controls: RealGasMassTransferSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasMassTransferCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !mass_kg.is_finite() || mass_kg < 0.0 || mass_kg >= self.mass_kg {
            return Err(GasEnergyError::Invalid(
                "finite nonnegative real-gas transfer smaller than donor inventory required",
            ));
        }
        if self.model != receiver.model {
            return Err(GasEnergyError::Invalid(
                "real-gas mass transfer requires identical frozen composition and PR/material sources",
            ));
        }
        let source_h = self.state.specific_enthalpy_j_kg();
        if mass_kg == 0.0 {
            return Ok(RealGasMassTransferCandidate {
                donor: self.clone(),
                receiver: receiver.clone(),
                step: RealGasMassTransferStep {
                    transferred_mass_kg: 0.0,
                    source_enthalpy_j_kg: source_h,
                    transferred_enthalpy_j: 0.0,
                    donor_energy_change_j: 0.0,
                    receiver_energy_change_j: 0.0,
                    energy_residual_j: 0.0,
                    mass_residual_kg: 0.0,
                    iterations: 0,
                },
            });
        }
        if core::ptr::eq(self, receiver) {
            return Err(GasEnergyError::Invalid(
                "nonzero real-gas transfer requires two distinct physical owners",
            ));
        }
        let donor_mass = self.mass_kg - mass_kg;
        let receiver_mass = receiver.mass_kg + mass_kg;
        if mass_kg <= controls.mass_tolerance_kg
            || !positive(donor_mass)
            || !positive(receiver_mass)
            || donor_mass == self.mass_kg
            || receiver_mass == receiver.mass_kg
            || ((self.mass_kg - donor_mass) - mass_kg).abs() > controls.mass_tolerance_kg
            || ((receiver_mass - receiver.mass_kg) - mass_kg).abs()
                > controls.mass_tolerance_kg
        {
            return Err(GasEnergyError::Invalid(
                "real-gas transfer is not representable within the absolute mass tolerance",
            ));
        }

        let flow_work_j_kg = self.state.pressure_pa() / self.density_kg_m3;
        let donor_change_j_kg = -(mass_kg / donor_mass) * flow_work_j_kg;
        let donor_to_receiver_u = self.model.specific_internal_energy_change_between_densities(
            receiver.density_kg_m3,
            receiver.state.temperature_k(),
            self.density_kg_m3,
            self.state.temperature_k(),
        )?;
        let receiver_change_j_kg =
            (mass_kg / receiver_mass) * (donor_to_receiver_u + flow_work_j_kg);
        if ![flow_work_j_kg, donor_change_j_kg, donor_to_receiver_u, receiver_change_j_kg]
            .into_iter().all(f64::is_finite)
        {
            return Err(GasEnergyError::Invalid(
                "real-gas transfer enthalpy target is unrepresentable",
            ));
        }

        let (donor, donor_residual, donor_iterations) = endpoint(
            self, donor_mass, donor_change_j_kg, controls.donor, &mut cancelled,
        )?;
        let (receiver, receiver_residual, receiver_iterations) = endpoint(
            receiver, receiver_mass, receiver_change_j_kg, controls.receiver, &mut cancelled,
        )?;
        let transferred_enthalpy_j = mass_kg * source_h;
        let donor_energy_change_j =
            -transferred_enthalpy_j + donor_mass * donor_residual;
        let receiver_energy_change_j =
            transferred_enthalpy_j + receiver_mass * receiver_residual;
        let energy_residual_j =
            donor_mass * donor_residual + receiver_mass * receiver_residual;
        let mass_residual_kg =
            (donor.mass_kg - self.mass_kg) + (receiver.mass_kg - receiver_mass + mass_kg);
        if ![
            transferred_enthalpy_j, donor_energy_change_j, receiver_energy_change_j,
            energy_residual_j, mass_residual_kg,
        ].into_iter().all(f64::is_finite)
            || energy_residual_j.abs() > controls.energy_tolerance_j
            || mass_residual_kg.abs() > controls.mass_tolerance_kg
        {
            return Err(GasEnergyError::Invalid(
                "real-gas paired transfer did not close mass and energy",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RealGasMassTransferCandidate {
            donor,
            receiver,
            step: RealGasMassTransferStep {
                transferred_mass_kg: mass_kg,
                source_enthalpy_j_kg: source_h,
                transferred_enthalpy_j,
                donor_energy_change_j,
                receiver_energy_change_j,
                energy_residual_j,
                mass_residual_kg,
                iterations: donor_iterations + receiver_iterations,
            },
        })
    }
}

#[cfg(test)]
mod tests;
