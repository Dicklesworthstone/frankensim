//! Retain a finite donor's original molecular beam in an evacuated rigid volume.
//!
//! The native vacuum transaction debits the donor once. Its exact accepted
//! species masses and molecular correction initialize the receiver through the
//! original caloric/EOS owner. No seed mass, pressure floor or vacuum temperature
//! is introduced. Once filled, both ordinary parcels can use reciprocal effusion.
//!
//! This is a first-order initial interval with zero incident receiver beam on
//! the original snapshot. All emitted molecules are retained and their complete
//! beam energy thermalizes into a well-mixed endpoint. Fixed rigid walls absorb
//! momentum without doing work; their inertia and kinetic transit are not
//! resolved. The model does not describe the nonequilibrium vacuum transient,
//! delayed molecular flight, adsorption, pumping or wall heat accommodation.

/// Adaptive initialization followed by the original finite-gas trajectory.
pub mod trajectory;

use super::{
    face_error,
    vacuum::{self, VacuumEffusionLink, VacuumEffusionStep},
};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        evacuated::{EvacuatedGasFillStep, EvacuatedGasVolume},
        transfer::species::{
            diffusion::network::DiffusionNetworkError,
            effusion::EffusiveInterface,
            entropy,
            network::transport::{GasTransportError, GasTransportPaths, GasTransportSolve},
        },
    },
};
use fs_math::dd::Dd;
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

/// Original aperture identity and the two distinct finite volume identities.
#[derive(Clone, Copy, Debug)]
pub struct EffusionFillLink<'a> {
    pub id: u64,
    pub donor: u64,
    pub receiver: u64,
    pub interface: &'a EffusiveInterface,
    /// Exposed area of the already qualified thin apertures [m2].
    pub area_m2: f64,
}

/// An inactive aperture leaves a true evacuated volume, without a gas state.
#[derive(Clone, Debug, PartialEq)]
pub enum EffusionReceiver {
    Evacuated(EvacuatedGasVolume),
    Filled(ClosedGasParcel),
}
impl EffusionReceiver {
    pub fn gas(&self) -> Option<&ClosedGasParcel> {
        match self {
            Self::Evacuated(_) => None,
            Self::Filled(gas) => Some(gas),
        }
    }
}

/// Independently measured initial-to-final conservation of the complete pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionFillBalance {
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Sum of actual final m_i * Delta u_i(T_donor, T_final) [J]. The initial
    /// donor is at that reference and the initial receiver has zero inventory.
    /// Original species references cancel before arithmetic; mass is checked
    /// separately. No absolute formation-energy difference scales this budget.
    pub energy_residual_j: f64,
    /// Complete closed-pair entropy change, including expansion into the
    /// initially empty volume through the original source-qualified helper.
    pub entropy_change_j_k: f64,
}

/// Original donor/receiver reports and their common closed-pair balance.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionFillStep {
    pub id: u64,
    pub donor: u64,
    pub receiver: u64,
    pub duration_s: f64,
    /// Native outgoing-beam transaction. Its external boundary rows are
    /// captured by THIS receiver, so they are internal to the complete pair.
    pub transport: VacuumEffusionStep,
    pub fill: Option<EvacuatedGasFillStep>,
    /// Exact accepted boundary masses, never rates integrated a second time.
    pub captured_species_kg: BTreeMap<SpeciesId, f64>,
    pub captured_mass_kg: f64,
    /// Numerical caloric chart, not a temperature assigned to the empty volume.
    pub reference_temperature_k: f64,
    /// U_receiver - sum_i m_i*u_i(reference_temperature_k) [J]. This comes
    /// directly from the original accepted molecular corrections.
    pub reference_relative_energy_j: f64,
    pub balance: EffusionFillBalance,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EffusionFillError {
    Invalid(&'static str),
    Transport(GasTransportError),
    Receiver(GasEnergyError),
    Conservation(EffusionFillBalance),
    /// A large explicit zero-return-beam interval may overfill the receiver.
    /// This physical refusal can be retried with a smaller first interval.
    EntropyDecrease {
        change_j_k: f64,
        tolerance_j_k: f64,
    },
    Cancelled,
}
impl core::fmt::Display for EffusionFillError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite molecular filling: {self:?}")
    }
}
impl std::error::Error for EffusionFillError {}
use EffusionFillError as Error;

/// Both complete endpoints, or neither. No caller-owned state is mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionFillCandidate {
    donor: ClosedGasParcel,
    receiver: EffusionReceiver,
    step: EffusionFillStep,
}
impl EffusionFillCandidate {
    pub const fn donor(&self) -> &ClosedGasParcel {
        &self.donor
    }
    pub const fn receiver(&self) -> &EffusionReceiver {
        &self.receiver
    }
    pub const fn step(&self) -> &EffusionFillStep {
        &self.step
    }
    pub fn into_parts(self) -> (ClosedGasParcel, EffusionReceiver, EffusionFillStep) {
        (self.donor, self.receiver, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn transport_error(error: GasTransportError) -> Error {
    if error == GasTransportError::Cancelled {
        Error::Cancelled
    } else {
        Error::Transport(error)
    }
}
fn receiver_error(error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Receiver(error)
    }
}
fn add(sum: &mut Dd, value: Dd) -> Result<(), Error> {
    *sum = *sum + value;
    if !value.hi.is_finite() || !value.lo.is_finite() || !sum.hi.is_finite() || !sum.lo.is_finite()
    {
        return Err(Error::Invalid(
            "captured molecular account is unrepresentable",
        ));
    }
    Ok(())
}

/// Audit a complete closed interval that began with one donor and an empty volume.
///
/// This also admits endpoints after subsequent ordinary reciprocal transport.
/// The original source owner checks actual total/species inventories, unchanged
/// donor volume, source-relative energy and complete entropy change. Every
/// physical budget applies once to the full interval, without step-count scaling.
///
/// The caller must retain the receiver's originally authored rigid volume.
/// External heat/work, reactions and moving boundaries are outside this closed
/// pair account. `None` represents a receiver that remains exactly empty; no
/// temperature or gas state is fabricated for that phase.
pub fn audit_pair(
    original: &ClosedGasParcel,
    donor: (u64, &ClosedGasParcel),
    receiver: Option<(u64, &ClosedGasParcel)>,
    controls: GasTransportSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EffusionFillBalance, Error> {
    poll(cancelled)?;
    if [
        controls.exchange.mass_tolerance_kg,
        controls.exchange.energy_tolerance_j,
        controls.entropy_tolerance_j_k,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value <= 0.0)
    {
        return Err(Error::Invalid(
            "positive finite pair conservation budgets required",
        ));
    }
    let change = match receiver {
        Some(receiver) => entropy::newly_filled_pair_change(
            original,
            donor,
            receiver,
            controls.exchange.mass_tolerance_kg,
            &mut *cancelled,
        ),
        None => entropy::closed_inventory_change(
            [(donor.0, original, donor.1)],
            controls.exchange.mass_tolerance_kg,
            &mut *cancelled,
        ),
    }
    .map_err(|error| match error {
        entropy::Error::Invalid(reason) => Error::Invalid(reason),
        entropy::Error::Parcel { id, error } => transport_error(GasTransportError::Diffusion(
            DiffusionNetworkError::Parcel { id, error },
        )),
        entropy::Error::Cancelled => Error::Cancelled,
    })?;
    let balance = EffusionFillBalance {
        total_mass_residual_kg: change.total_mass_residual_kg,
        maximum_species_mass_residual_kg: change.maximum_species_mass_residual_kg,
        energy_residual_j: change.energy_change_j,
        entropy_change_j_k: change.entropy_change_j_k,
    };
    if balance.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || balance.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
        || balance.energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(Error::Conservation(balance));
    }
    if balance.entropy_change_j_k < -controls.entropy_tolerance_j_k {
        return Err(Error::EntropyDecrease {
            change_j_k: balance.entropy_change_j_k,
            tolerance_j_k: controls.entropy_tolerance_j_k,
        });
    }
    poll(cancelled)?;
    Ok(balance)
}

/// Fill a true evacuated volume from one original finite molecular donor.
///
/// The original vacuum owner samples the initial zero-return-beam state and
/// enforces every gross-species withdrawal limit. Its accepted outgoing masses
/// enter the receiver unchanged. Its positive relative-energy correction is
/// exactly the captured `sum dm_i R T_donor/(2 M_i)` relative to the incoming
/// species' ORIGINAL internal energies at T_donor. Thus the empty receiver's
/// temperature follows its captured beam energy and actual composition; it is
/// neither prescribed nor copied from the donor.
///
/// Both actual endpoint gas states must satisfy the same aperture qualification.
/// A tiny positive fill below a source pressure domain refuses without applying
/// a floor. Zero area or zero transmitted inventory preserves the evacuated
/// owner. The original solver controls and ONE fixed pair mass/energy budget
/// apply. A refusal, including late cancellation, publishes neither endpoint.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    link: EffusionFillLink<'_>,
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EffusionFillCandidate, Error> {
    poll(&mut cancelled)?;
    if link.donor == link.receiver {
        return Err(Error::Invalid(
            "donor and evacuated receiver need distinct identities",
        ));
    }
    receiver_solve.validate().map_err(receiver_error)?;
    let original = BTreeMap::from([(link.donor, donor.clone())]);
    let endpoints = BTreeMap::from([(link.donor, donor_solve)]);
    let outlet = VacuumEffusionLink {
        id: link.id,
        parcel: link.donor,
        interface: link.interface,
        area_m2: link.area_m2,
    };
    let (mut next, transport) = vacuum::step_candidate(
        &original,
        GasTransportPaths::default(),
        &[outlet],
        duration_s,
        &endpoints,
        controls,
        &mut cancelled,
    )
    .map_err(transport_error)?
    .into_parts();
    let next_donor = next.remove(&link.donor).ok_or(Error::Invalid(
        "native vacuum transaction lost its original donor",
    ))?;
    let mut captured = BTreeMap::new();
    let mut mass = Dd::ZERO;
    let mut relative = Dd::ZERO;
    for row in &transport.exchange.boundaries {
        poll(&mut cancelled)?;
        let beam = &row.transfer;
        if row.id != link.id
            || row.parcel != link.donor
            || beam.incoming_mass_kg != 0.0
            || beam.outgoing_mass_kg <= 0.0
            || beam.net_mass_into_gas_kg != -beam.outgoing_mass_kg
            || captured
                .insert(beam.species.clone(), beam.outgoing_mass_kg)
                .is_some()
        {
            return Err(Error::Invalid(
                "native captured beam changed identity or direction",
            ));
        }
        add(&mut mass, Dd::from_f64(beam.outgoing_mass_kg))?;
        add(&mut relative, Dd::from_f64(beam.relative_energy_into_gas_j))?;
    }
    let reference_temperature_k = donor.state().state().temperature;
    let (next_receiver, fill) = if captured.is_empty() {
        (EffusionReceiver::Evacuated(receiver.clone()), None)
    } else {
        let (gas, step) = receiver
            .fill_candidate(
                donor.model(),
                &captured,
                reference_temperature_k,
                relative.to_f64(),
                receiver_solve,
                controls.exchange.mass_tolerance_kg,
                controls.exchange.energy_tolerance_j,
                &mut cancelled,
            )
            .map_err(receiver_error)?
            .into_parts();
        for endpoint in [&next_donor, &gas] {
            link.interface
                .admit_parcel(endpoint, &mut cancelled)
                .map_err(|error| transport_error(face_error(link.id, error).into()))?;
        }
        (EffusionReceiver::Filled(gas), Some(step))
    };
    let balance = audit_pair(
        donor,
        (link.donor, &next_donor),
        next_receiver.gas().map(|gas| (link.receiver, gas)),
        controls,
        &mut cancelled,
    )?;
    let step = EffusionFillStep {
        id: link.id,
        donor: link.donor,
        receiver: link.receiver,
        duration_s,
        transport,
        fill,
        captured_species_kg: captured,
        captured_mass_kg: mass.to_f64(),
        reference_temperature_k,
        reference_relative_energy_j: relative.to_f64(),
        balance,
    };
    poll(&mut cancelled)?;
    Ok(EffusionFillCandidate {
        donor: next_donor,
        receiver: next_receiver,
        step,
    })
}

#[cfg(test)]
mod tests;
