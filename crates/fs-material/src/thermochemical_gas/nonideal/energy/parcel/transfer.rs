//! Conservative upstream-enthalpy transfer between two fixed real-gas volumes.
//!
//! A prescribed positive mass carries its upstream species and total enthalpy
//! into the receiver once. Compatible source models may have different species
//! fractions or active species; every new PR pair interaction must already be
//! supplied by an original source. Each fixed vessel reconstructs its new
//! density, composition, temperature and pressure through the original owners.
//! No ideal-gas pressure, gamma or caloric surrogate is introduced.
//!
//! This is a first-order well-mixed finite-volume transfer primitive. It does
//! not choose a flow rate, resolve a nozzle/jet, evolve momentum, perform phase
//! equilibrium or infer missing cross-species interaction data. A higher-level
//! pressure/port owner must choose and refine the transferred mass.

/// Quasi-steady pressure-driven flow using the same source-backed PR branch.
pub mod nozzle;
mod mixing;

use super::IsochoricRealGasParcel;
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, poll, updated_point,
};
use crate::thermochemical_gas::positive;

/// Original canonical species-flux and extensive-balance report, now also used
/// by PR storage. Its energy entries include the original nonideal departures;
/// conserved species datums cancel before evaluating the reported residual.
pub use crate::thermochemical_gas::energy::parcel::transfer::{
    GasMassTransferStep as RealGasMassTransferStep, GasSpeciesTransfer,
};

/// Bounded endpoint solves and absolute extensive conservation tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasMassTransferSolve {
    /// Donor caloric bracket and specific-energy tolerance.
    pub donor: GasEnergySolve,
    /// Receiver caloric bracket and specific-energy tolerance.
    pub receiver: GasEnergySolve,
    /// Absolute total AND individual-species mass/representability tolerance [kg].
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
    /// Every transported species keeps its original NASA, transport, critical
    /// properties and pack bindings. Shared source laws and known pair entries
    /// must agree exactly; a pair absent from BOTH catalogs is refused, including
    /// when each separate catalog uses zero interactions. No missing cross-pair
    /// coefficient is guessed. The receiver combines extensive species amounts.
    ///
    /// The donor loses `dm*h_d`; the receiver gains the same quantity. The
    /// energy solves are written in reference-safe differences:
    /// donor `Delta u = -(dm/m1) p_d/rho_d`, while the receiver uses
    /// `(dm/m1)[u_d-u_r+p_d/rho_d]` for identical models. Changing composition
    /// uses the new mixture at receiver T0 and final density as the energy anchor.
    /// Its ideal species datums cancel by inventory before arithmetic; original
    /// NASA temperature increments and all three PR departure energies remain.
    /// Incoming species and the new mixture must support this common anchor and
    /// the caloric bracket. Large formation constants do not enter the root.
    /// Actual species-inventory roundoff has a separate absolute mass budget;
    /// the pair residual does not subtract enormous displayed energy reports.
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
                    maximum_species_mass_residual_kg: 0.0,
                    species: Vec::new(),
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

        let mixed = if self.model == receiver.model { None } else {
            Some(mixing::receiver_model(self, receiver, mass_kg, receiver_mass)?)
        };
        let flow_work_j_kg = self.state.pressure_pa() / self.density_kg_m3;
        let donor_change_j_kg = -(mass_kg / donor_mass) * flow_work_j_kg;
        if ![flow_work_j_kg, donor_change_j_kg]
            .into_iter().all(f64::is_finite)
        {
            return Err(GasEnergyError::Invalid(
                "real-gas transfer enthalpy target is unrepresentable",
            ));
        }

        let (donor, donor_residual, donor_iterations) = endpoint(
            self, donor_mass, donor_change_j_kg, controls.donor, &mut cancelled,
        )?;
        let old_receiver = receiver;
        let (receiver, receiver_residual, receiver_iterations) = if let Some(model) = mixed {
            mixing::receiver_endpoint(self, old_receiver, model, receiver_mass, mass_kg,
                controls.receiver, &mut cancelled)?
        } else {
            // Preserve the original identical-model arithmetic and solve.
            let donor_to_receiver_u = self.model.specific_internal_energy_change_between_densities(
                receiver.density_kg_m3, receiver.state.temperature_k(),
                self.density_kg_m3, self.state.temperature_k(),
            )?;
            let receiver_change_j_kg =
                (mass_kg / receiver_mass) * (donor_to_receiver_u + flow_work_j_kg);
            if !donor_to_receiver_u.is_finite() || !receiver_change_j_kg.is_finite() {
                return Err(GasEnergyError::Invalid("real-gas receiving enthalpy target is unrepresentable"));
            }
            endpoint(receiver, receiver_mass, receiver_change_j_kg, controls.receiver, &mut cancelled)?
        };
        let (species, maximum_species_mass_residual_kg) = mixing::species_balance(
            self, old_receiver, &donor, &receiver, mass_kg,
        )?;
        let transferred_enthalpy_j = mass_kg * source_h;
        let donor_energy_change_j =
            -transferred_enthalpy_j + donor_mass * donor_residual;
        let receiver_energy_change_j =
            transferred_enthalpy_j + receiver_mass * receiver_residual;
        let energy_residual_j =
            donor_mass * donor_residual + receiver_mass * receiver_residual;
        let mass_residual_kg =
            (donor.mass_kg - self.mass_kg) + (receiver.mass_kg - old_receiver.mass_kg);
        if ![
            transferred_enthalpy_j, donor_energy_change_j, receiver_energy_change_j,
            energy_residual_j, mass_residual_kg,
        ].into_iter().all(f64::is_finite)
            || energy_residual_j.abs() > controls.energy_tolerance_j
            || mass_residual_kg.abs() > controls.mass_tolerance_kg
            || maximum_species_mass_residual_kg > controls.mass_tolerance_kg
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
                maximum_species_mass_residual_kg,
                species,
                iterations: donor_iterations + receiver_iterations,
            },
        })
    }
}

#[cfg(test)]
mod tests;
