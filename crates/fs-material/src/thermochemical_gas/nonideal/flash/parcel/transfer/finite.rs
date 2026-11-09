//! Finite-tank transfer with integrated original flow enthalpy and bulk composition.
//!
//! The original donor S,V owner supplies the well-mixed adiabatic withdrawal.
//! Its integrated outgoing H, including cooling during withdrawal, drives the
//! original receiving U,V/species owner. A finite donor is never replaced by a
//! held reservoir, and no averaged outlet temperature or Cp is constructed.

use super::super::discharge::{
    MixtureEquilibriumDischargeControl, MixtureEquilibriumDischargeStep,
};
use super::{
    IsochoricMixtureEquilibriumParcel, MixtureEquilibriumInflow, MixtureEquilibriumInflowControl,
    MixtureEquilibriumInflowError, ThermochemicalFlashError, ThermochemicalPengRobinsonModel,
    inventory, positive, products, source_specific_volume, sum,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

/// Native donor/receiver work caps and independent whole-pair balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureEquilibriumTransferControl {
    /// Original S,V withdrawal; mass, phase and source domains remain authoritative.
    pub donor: MixtureEquilibriumDischargeControl,
    /// Original receiving composition join and U,V closure. Its PT work cap
    /// covers only NEW reference/trial work; the actual donor anchor is reused.
    pub receiver: MixtureEquilibriumInflowControl,
    /// Absolute paired internal-energy discrepancy allowance [J].
    pub energy_tolerance_j: f64,
    /// Allowed numerical shortfall from nonnegative total entropy production [J/K].
    pub entropy_tolerance_j_k: f64,
}
impl MixtureEquilibriumTransferControl {
    pub fn new(
        donor: MixtureEquilibriumDischargeControl,
        receiver: MixtureEquilibriumInflowControl,
    ) -> Self {
        Self {
            donor,
            receiver,
            energy_tolerance_j: 1e-3,
            entropy_tolerance_j_k: 1e-6,
        }
    }
    pub fn validate(self) -> Result<(), MixtureEquilibriumTransferError> {
        self.donor.validate().map_err(donor_error)?;
        self.receiver.validate().map_err(receiver_error)?;
        if !positive(self.energy_tolerance_j)
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k < 0.0
        {
            return Err(MixtureEquilibriumTransferError::Invalid(
                "finite paired energy/entropy allowances required",
            ));
        }
        Ok(())
    }
}

/// Every refusal leaves both original physical inventories unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureEquilibriumTransferError {
    Invalid(&'static str),
    Donor(ThermochemicalFlashError),
    Receiver(MixtureEquilibriumInflowError),
    MassBalance {
        species: Option<SpeciesId>,
        residual_kg: f64,
        tolerance_kg: f64,
    },
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Entropy {
        production_j_k: f64,
        tolerance_j_k: f64,
    },
    PressureHead {
        donor_pressure_pa: f64,
        receiver_pressure_pa: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MixtureEquilibriumTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite equilibrium tank transfer: {self:?}")
    }
}
impl std::error::Error for MixtureEquilibriumTransferError {}
type Error = MixtureEquilibriumTransferError;

/// Actual canonical species masses, ordered [donor, receiver].
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureEquilibriumSpeciesTransfer {
    pub species: SpeciesId,
    pub initial_mass_kg: [f64; 2],
    pub final_mass_kg: [f64; 2],
    /// Original bulk-composition portion delivered to the receiving source join.
    pub transferred_mass_kg: f64,
    pub residual_kg: f64,
}

/// Accepted paired storage and transport, with no invented mean stream state.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureEquilibriumTransferStep {
    /// Native finite withdrawal, including its datum-free integrated H correction.
    pub discharge: MixtureEquilibriumDischargeStep,
    pub species: Vec<MixtureEquilibriumSpeciesTransfer>,
    /// Actual total pair mass discrepancy, distinct from prescribed mass rounding.
    pub mass_residual_kg: f64,
    /// Exactly the same integrated H leaving the donor and entering the receiver.
    /// This original-datum absolute report may be negative; it is not added heat.
    pub transported_enthalpy_j: f64,
    /// Original receiving inventory U change; use the relative balances below
    /// when absolute source datums would erase a small difference.
    pub receiver_internal_energy_change_j: f64,
    pub receiver_energy_residual_j: f64,
    pub receiver_volume_residual_m3: f64,
    /// Native donor storage plus integrated H, with species references cancelled.
    pub donor_energy_residual_j: f64,
    /// Original paired U change at unchanged vessel volumes [J].
    pub combined_energy_residual_j: f64,
    /// Receiver S change minus transported donor S, without datum subtraction.
    pub receiver_entropy_production_j_k: f64,
    /// Receiver production plus native donor S,V entropy residual [J/K].
    pub entropy_production_j_k: f64,
    /// Actual pressure/PH/PS inversions across both native owners.
    pub pressure_evaluations: usize,
    /// New native PT selections across both solves; no donor anchor reflash.
    pub pt_evaluations: usize,
    /// Sum of original EOS work ceilings, not an observed cubic call count.
    pub eos_evaluation_upper_bound: usize,
}

/// Transfer one prescribed finite bulk portion between two rigid equilibrium tanks.
///
/// Both phases are entrained in their instantaneous equilibrium proportions,
/// without slip. The donor retains its overall species composition and specific
/// entropy while its T,p and phase amounts evolve. Its original S,V endpoint
/// gives `H_out = dm*u_initial - m_remaining*Delta_u`, integrated over withdrawal.
/// The receiving source joins those original species and solves
/// `U_after = U_before + H_out` at its unchanged physical vessel volume.
/// Pressure-flow work is already inside H_out; no second pV or latent debit occurs.
///
/// Initial donor pressure must exceed initial receiver pressure; final donor
/// pressure must also cover the final receiver. These necessary endpoint checks
/// do not prove the intervening pressure path. No duration, valve/nozzle rate,
/// selective vapor/liquid venting, resolved jet energy, heat or shaft work is
/// inferred. Subdivide a portion when its actual path requires a finer model.
/// Both tanks must remain nonempty. Catalog/source/phase compatibility is owned
/// by the original receiver; an optional original catalog supplies known cross
/// interactions, never fabricated pair data. Zero requested mass is exact idle.
///
/// One native S,V and one U,V closure have separately declared finite work caps.
/// Candidates remain private until all paired balances and final cancellation
/// pass. Neither owner changes on a late receiving or pressure-head refusal.
pub fn transfer_adiabatic_mixture_equilibrium(
    donor: &mut IsochoricMixtureEquilibriumParcel,
    receiver: &mut IsochoricMixtureEquilibriumParcel,
    requested_mass_kg: f64,
    interaction_catalog: Option<&ThermochemicalPengRobinsonModel>,
    controls: MixtureEquilibriumTransferControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureEquilibriumTransferStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if !requested_mass_kg.is_finite() || requested_mass_kg < 0.0 {
        return Err(Error::Invalid(
            "finite nonnegative transferred mass required",
        ));
    }
    if requested_mass_kg > 0.0 {
        if donor.state().pressure_pa() <= receiver.state().pressure_pa() {
            return Err(Error::PressureHead {
                donor_pressure_pa: donor.state().pressure_pa(),
                receiver_pressure_pa: receiver.state().pressure_pa(),
            });
        }
        if donor.state().phase().options() != controls.receiver.equilibrium.phase.phase_options()
            || receiver.state().phase().options()
                != controls.receiver.equilibrium.phase.phase_options()
        {
            return Err(Error::Invalid(
                "finite transfer requires a common original phase policy",
            ));
        }
    }
    let mut next_donor = donor.clone();
    let mut next_receiver = receiver.clone();
    let discharge = next_donor
        .discharge_adiabatic(
            requested_mass_kg,
            Pressure::new(receiver.state().pressure_pa()),
            controls.donor,
            &mut cancelled,
        )
        .map_err(donor_error)?;
    // Subtraction rounding belongs to the donor's measured debit. Credit that
    // SAME actual amount, not a second independently rounded requested portion.
    let amount = discharge.removed_mass_kg;
    let inflow = next_receiver
        .receive_accounted_stream(
            MixtureEquilibriumInflow {
                model: donor.model(),
                temperature: Temperature::new(donor.state().temperature_k()),
                pressure: Pressure::new(donor.state().pressure_pa()),
                mass_kg: amount,
            },
            interaction_catalog,
            controls.receiver,
            Some(donor.state()),
            discharge.exported_enthalpy_adjustment_j,
            &mut cancelled,
        )
        .map_err(receiver_error)?;
    if amount > 0.0 && next_donor.state().pressure_pa() < next_receiver.state().pressure_pa() {
        return Err(Error::PressureHead {
            donor_pressure_pa: next_donor.state().pressure_pa(),
            receiver_pressure_pa: next_receiver.state().pressure_pa(),
        });
    }
    let mass_residual = sum([
        next_donor.mass_kg(),
        -donor.mass_kg(),
        next_receiver.mass_kg(),
        -receiver.mass_kg(),
    ]);
    let mass_scale = donor
        .mass_kg()
        .max(receiver.mass_kg())
        .max(next_receiver.mass_kg());
    admit_mass(None, mass_residual, mass_scale, controls.receiver)?;
    let species = audit_species(
        donor,
        receiver,
        &next_donor,
        &next_receiver,
        amount,
        controls.receiver,
        &mut cancelled,
    )?;
    // Independently replay donor Delta-U + delivered H, after its common u0
    // reference cancels. The native correction already integrates its cooling.
    let donor_residual = products([
        (
            next_donor.mass_kg(),
            discharge.remaining_specific_internal_energy_change_j_kg,
        ),
        (
            amount,
            donor.state().pressure_pa() * source_specific_volume(donor.state()),
        ),
        (1.0, discharge.exported_enthalpy_adjustment_j),
    ]);
    let energy = sum([donor_residual, inflow.energy_residual_j]);
    let entropy = sum([
        discharge.entropy_balance_residual_j_k,
        inflow.entropy_production_j_k,
    ]);
    if !donor_residual.is_finite()
        || !energy.is_finite()
        || energy.abs() > controls.energy_tolerance_j
        || donor_residual.abs() + inflow.energy_residual_j.abs() > controls.energy_tolerance_j
        || discharge.exported_enthalpy_j != inflow.incoming_enthalpy_j
    {
        return Err(Error::Energy {
            residual_j: energy,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    if !entropy.is_finite() || entropy < -controls.entropy_tolerance_j_k {
        return Err(Error::Entropy {
            production_j_k: entropy,
            tolerance_j_k: controls.entropy_tolerance_j_k,
        });
    }
    let report = MixtureEquilibriumTransferStep {
        discharge,
        species,
        mass_residual_kg: mass_residual,
        transported_enthalpy_j: discharge.exported_enthalpy_j,
        receiver_internal_energy_change_j: inflow.stored_internal_energy_change_j,
        receiver_energy_residual_j: inflow.energy_residual_j,
        receiver_volume_residual_m3: inflow.volume_residual_m3,
        donor_energy_residual_j: donor_residual,
        combined_energy_residual_j: energy,
        receiver_entropy_production_j_k: inflow.entropy_production_j_k,
        entropy_production_j_k: entropy,
        pressure_evaluations: add_work(
            discharge.pressure_evaluations,
            inflow.pressure_evaluations,
        )?,
        pt_evaluations: add_work(discharge.pt_evaluations, inflow.pt_evaluations)?,
        eos_evaluation_upper_bound: add_work(
            discharge.eos_evaluation_upper_bound,
            inflow.eos_evaluation_upper_bound,
        )?,
    };
    poll(&mut cancelled)?;
    *donor = next_donor;
    *receiver = next_receiver;
    Ok(report)
}

fn audit_species(
    donor: &IsochoricMixtureEquilibriumParcel,
    receiver: &IsochoricMixtureEquilibriumParcel,
    next_donor: &IsochoricMixtureEquilibriumParcel,
    next_receiver: &IsochoricMixtureEquilibriumParcel,
    amount: f64,
    controls: MixtureEquilibriumInflowControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<MixtureEquilibriumSpeciesTransfer>, Error> {
    let mut amounts: BTreeMap<SpeciesId, [f64; 5]> = BTreeMap::new();
    for (index, (model, mass)) in [
        (donor.model(), donor.mass_kg()),
        (receiver.model(), receiver.mass_kg()),
        (next_donor.model(), next_donor.mass_kg()),
        (next_receiver.model(), next_receiver.mass_kg()),
        (donor.model(), amount),
    ]
    .into_iter()
    .enumerate()
    {
        poll(cancelled)?;
        // Native inventory accepts only positive portions; idle has no transfer.
        if mass == 0.0 {
            continue;
        }
        for (species, value) in inventory(model.source(), mass)
            .map_err(|e| receiver_error(MixtureEquilibriumInflowError::Catalog(e)))?
        {
            amounts.entry(species).or_insert([0.0; 5])[index] = value;
        }
    }
    let mut report = Vec::with_capacity(amounts.len());
    for (species, [d0, r0, d1, r1, transferred]) in amounts {
        poll(cancelled)?;
        let residual = sum([d1, -d0, r1, -r0]);
        admit_mass(
            Some(species.clone()),
            residual,
            d0.max(r0).max(d1).max(r1),
            controls,
        )?;
        report.push(MixtureEquilibriumSpeciesTransfer {
            species,
            initial_mass_kg: [d0, r0],
            final_mass_kg: [d1, r1],
            transferred_mass_kg: transferred,
            residual_kg: residual,
        });
    }
    Ok(report)
}
fn admit_mass(
    species: Option<SpeciesId>,
    residual: f64,
    scale: f64,
    controls: MixtureEquilibriumInflowControl,
) -> Result<(), Error> {
    let tolerance = controls.mass_tolerance(scale);
    if !positive(tolerance) || !residual.is_finite() || residual.abs() > tolerance {
        return Err(Error::MassBalance {
            species,
            residual_kg: residual,
            tolerance_kg: tolerance,
        });
    }
    Ok(())
}
fn donor_error(error: ThermochemicalFlashError) -> Error {
    if error == ThermochemicalFlashError::Cancelled {
        Error::Cancelled
    } else {
        Error::Donor(error)
    }
}
fn add_work(left: usize, right: usize) -> Result<usize, Error> {
    left.checked_add(right)
        .ok_or(Error::Invalid("paired native work counter overflow"))
}
fn receiver_error(error: MixtureEquilibriumInflowError) -> Error {
    if error == MixtureEquilibriumInflowError::Cancelled {
        Error::Cancelled
    } else {
        Error::Receiver(error)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
