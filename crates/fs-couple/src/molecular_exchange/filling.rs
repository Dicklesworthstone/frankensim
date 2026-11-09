//! Molecular startup of an exactly evacuated finite moving-gas receiver.
//!
//! One original drifting Maxwellian supplies the SAME accepted mass, momentum
//! and energy to both finite endpoints. The receiver retains the captured bulk
//! momentum; only energy beyond its actual bulk K enters its original caloric
//! state. Once filled, the original reciprocal pair owns every subsequent step.
//! The homogeneous endpoint is an admitted mean-drift closure, not a resolved
//! collisionless filling front, container motion or remaining wall-force model.

use super::{
    EffusiveDriftInterface, GasEnergyError, MAX_EFFUSIVE_SPECIES, MolecularExchangeControls,
    MolecularExchangeError, MolecularExchangeStep, MolecularExchangeSurface, MovingGasBoundaryStep,
    MovingGasParcel, SpeciesBoundaryRelativeTransfer, SpeciesId, admit, allowance, balance, finite,
    moving_error, poll, positive,
};
use crate::liquid_gas::evaporation::momentum::boundary::evacuated::MovingGasFillStep;
use fs_material::thermochemical_gas::energy::parcel::{
    evacuated::EvacuatedGasVolume, transfer::species::effusion::drift::DriftingEffusiveFlux,
};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

/// A true empty geometry or the original positive-inventory moving gas.
#[derive(Clone, Debug, PartialEq)]
pub enum MolecularReceiver {
    Evacuated(EvacuatedGasVolume),
    Filled(MovingGasParcel),
}
impl MolecularReceiver {
    pub const fn gas(&self) -> Option<&MovingGasParcel> {
        match self {
            Self::Evacuated(_) => None,
            Self::Filled(gas) => Some(gas),
        }
    }
}

#[derive(Debug)]
pub enum MolecularFillError {
    Exchange(MolecularExchangeError),
    Caloric {
        species: SpeciesId,
        error: ThermochemErrorV1,
    },
}
impl From<MolecularExchangeError> for MolecularFillError {
    fn from(error: MolecularExchangeError) -> Self {
        Self::Exchange(error)
    }
}
impl From<GasEnergyError> for MolecularFillError {
    fn from(error: GasEnergyError) -> Self {
        Self::Exchange(error.into())
    }
}
impl core::fmt::Display for MolecularFillError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving molecular receiver filling refused: {self:?}")
    }
}
impl std::error::Error for MolecularFillError {}
type Error = MolecularFillError;

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularCapturedSpecies {
    pub species: SpeciesId,
    pub mass_kg: f64,
    pub momentum_aperture_n_s: [f64; 3],
    /// Complete beam energy beyond mass times original u_i(T_donor).
    pub reference_energy_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularFirstFillStep {
    pub flux: DriftingEffusiveFlux,
    pub captured: Vec<MolecularCapturedSpecies>,
    pub donor: MovingGasBoundaryStep,
    pub receiver: MovingGasFillStep,
    pub captured_mass_kg: f64,
    pub captured_momentum_n_s: [f64; 3],
    /// Caloric coordinate ONLY; the initial vacuum has no temperature.
    pub reference_temperature_k: f64,
    pub captured_reference_energy_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub momentum_residual_n_s: [f64; 3],
    /// Independently reconstructed original-to-final U plus bulk K.
    pub energy_residual_j: f64,
    pub laboratory_energy_residual_j: f64,
    /// Native donor/receiver defects and reconstruction differences cannot cancel.
    pub absolute_energy_defect_j: f64,
    pub mass_allowance_kg: f64,
    pub momentum_allowance_n_s: [f64; 3],
    pub energy_allowance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MolecularFillStage {
    /// No received inventory; keep the original donor and exact vacuum.
    Closed(DriftingEffusiveFlux),
    FirstFill(MolecularFirstFillStep),
    Reciprocal(MolecularExchangeStep),
}
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularFillStep {
    pub duration_s: f64,
    pub surface: MolecularExchangeSurface,
    pub stage: MolecularFillStage,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularFillCandidate {
    donor: MovingGasParcel,
    receiver: MolecularReceiver,
    step: MolecularFillStep,
}
impl MolecularFillCandidate {
    pub const fn donor(&self) -> &MovingGasParcel {
        &self.donor
    }
    pub const fn receiver(&self) -> &MolecularReceiver {
        &self.receiver
    }
    pub const fn step(&self) -> &MolecularFillStep {
        &self.step
    }
    pub fn into_parts(self) -> (MovingGasParcel, MolecularReceiver, MolecularFillStep) {
        (self.donor, self.receiver, self.step)
    }
}

fn represented(value: Dd) -> Result<f64, Error> {
    let result = finite(value)?;
    if value != Dd::ZERO && result == 0.0 {
        return Err(
            MolecularExchangeError::Invalid("captured molecular coordinate underflow").into(),
        );
    }
    Ok(result)
}

/// First capture from the donor at the left of the authored aperture. An
/// already filled receiver delegates exactly to the native reciprocal pair.
///
/// All rates are sampled at the interval start. The same original aperture,
/// drift and caloric laws must admit both actual positive endpoints. No return
/// beam is fabricated for vacuum, nor is initial directional beam energy
/// thermalized twice. Only the complete candidate is returned on success.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    donor: &MovingGasParcel,
    receiver: &MolecularReceiver,
    interface: &EffusiveDriftInterface,
    surface: MolecularExchangeSurface,
    duration_s: f64,
    controls: MolecularExchangeControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MolecularFillCandidate, Error> {
    poll(&mut cancelled)?;
    let empty = match receiver {
        MolecularReceiver::Filled(gas) => {
            let result = super::step_candidate(
                [donor, gas],
                interface,
                surface,
                duration_s,
                controls,
                &mut cancelled,
            )?;
            let ([donor, receiver], report) = result.into_parts();
            poll(&mut cancelled)?;
            return Ok(MolecularFillCandidate {
                donor,
                receiver: MolecularReceiver::Filled(receiver),
                step: MolecularFillStep {
                    duration_s,
                    surface,
                    stage: MolecularFillStage::Reciprocal(report),
                },
            });
        }
        MolecularReceiver::Evacuated(volume) => volume,
    };
    let frame = admit(surface, duration_s, controls)?;
    let active = surface.area_m2 > 0.0;
    let velocity = donor.velocity_m_s();
    let relative = if active {
        core::array::from_fn(|d| velocity[d] - frame[d])
    } else {
        [0.0; 3]
    };
    let flux = donor.parcel().drifting_effusive_flux(
        interface,
        surface.area_m2,
        relative,
        surface.normal_left_to_right.map(|x| -x),
        &mut cancelled,
    )?;
    if flux.species.len() > MAX_EFFUSIVE_SPECIES {
        return Err(MolecularExchangeError::Invalid(
            "captured source exceeds the original species bound",
        )
        .into());
    }
    let catalog = interface.momentum_interface().interface().gas_catalog();
    let inventory = if active {
        donor.parcel().species_inventory_kg()?
    } else {
        BTreeMap::new()
    };
    let reference = donor.parcel().state().state().temperature;
    let mut captured = Vec::with_capacity(flux.species.len());
    let mut received = BTreeMap::new();
    let mut rows = Vec::with_capacity(flux.species.len());
    let mut mass = Dd::ZERO;
    let mut impulse = [Dd::ZERO; 3];
    let mut relative_energy = Dd::ZERO;
    let mut energy_scale = Dd::ZERO;
    for beam in &flux.species {
        poll(&mut cancelled)?;
        if beam.incoming.mass_rate_kg_s == 0.0 {
            continue;
        }
        let dm =
            represented(Dd::from_f64(beam.incoming.mass_rate_kg_s) * Dd::from_f64(duration_s))?;
        let available = inventory.get(&beam.species).copied().unwrap_or(0.0);
        if !positive(available) {
            return Err(
                MolecularExchangeError::Invalid("capture lost its original donor species").into(),
            );
        }
        if (Dd::from_f64(available) * Dd::from_f64(controls.maximum_species_fraction))
            .lt(Dd::from_f64(dm))
        {
            return Err(MolecularExchangeError::Withdrawal {
                side: 0,
                species: beam.species.clone(),
                fraction: dm / available,
                maximum: controls.maximum_species_fraction,
            }
            .into());
        }
        let correction = beam.incoming.additional_specific_energy_j_kg.ok_or(
            MolecularExchangeError::Invalid(
                "active captured beam lost its crossing-energy correction",
            ),
        )?;
        let index = catalog
            .thermodynamics()
            .components()
            .binary_search_by(|source| source.species().cmp(&beam.species))
            .map_err(|_| {
                MolecularExchangeError::Invalid("captured beam lost its original caloric source")
            })?;
        let molar_mass = catalog.thermodynamics().components()[index]
            .molar_mass()
            .value();
        let relative_h = represented(Dd::from_f64(dm) * Dd::from_f64(correction))?;
        let relative_u = represented(
            Dd::from_f64(dm)
                * (Dd::from_f64(R) * Dd::from_f64(reference) / Dd::from_f64(molar_mass)
                    + Dd::from_f64(correction)),
        )?;
        let momentum = beam.impulse_for_mass_kg(dm)?;
        mass = mass + Dd::from_f64(dm);
        relative_energy = relative_energy + Dd::from_f64(relative_u);
        energy_scale =
            energy_scale + Dd::from_f64(relative_h).abs() + Dd::from_f64(relative_u).abs();
        for d in 0..3 {
            impulse[d] = impulse[d] + Dd::from_f64(momentum[d]);
        }
        received.insert(beam.species.clone(), dm);
        rows.push(SpeciesBoundaryRelativeTransfer {
            species: beam.species.clone(),
            mass_into_gas_kg: -dm,
            relative_energy_into_gas_j: -relative_h,
        });
        captured.push(MolecularCapturedSpecies {
            species: beam.species.clone(),
            mass_kg: dm,
            momentum_aperture_n_s: momentum,
            reference_energy_j: relative_u,
        });
    }
    if received.is_empty() {
        poll(&mut cancelled)?;
        return Ok(MolecularFillCandidate {
            donor: donor.clone(),
            receiver: receiver.clone(),
            step: MolecularFillStep {
                duration_s,
                surface,
                stage: MolecularFillStage::Closed(flux),
            },
        });
    }
    let captured_mass_kg = represented(mass)?;
    for d in 0..3 {
        impulse[d] = impulse[d] + Dd::from_f64(frame[d]) * mass;
    }
    let captured_momentum_n_s = [
        represented(impulse[0])?,
        represented(impulse[1])?,
        represented(impulse[2])?,
    ];
    let captured_reference_energy_j = represented(relative_energy)?;
    let left = donor
        .exchange_species_momentum_relative(
            catalog,
            &rows,
            captured_momentum_n_s.map(|p| -p),
            frame,
            controls.gas[0],
            &mut cancelled,
        )
        .map_err(|error| moving_error(0, error))?;
    let right = MovingGasParcel::fill_evacuated_candidate(
        empty,
        catalog,
        &received,
        reference,
        captured_reference_energy_j,
        captured_momentum_n_s,
        frame,
        controls.gas[1],
        &mut cancelled,
    )
    .map_err(|error| moving_error(1, error))?;
    let (next_donor, _, donor_step) = left.into_parts();
    let (next_receiver, receiver_step) = right.into_parts();
    for gas in [&next_donor, &next_receiver] {
        poll(&mut cancelled)?;
        interface
            .admit_relative_velocity(core::array::from_fn(|d| gas.velocity_m_s()[d] - frame[d]))?;
        interface
            .momentum_interface()
            .interface()
            .admit_parcel(gas.parcel(), &mut cancelled)?;
    }
    let after = [
        next_donor.parcel().species_inventory_kg()?,
        next_receiver.parcel().species_inventory_kg()?,
    ];
    let mass_allowance_kg =
        allowance(controls.mass_absolute_kg, controls.relative_tolerance, mass)?;
    let mass_residual = Dd::from_f64(next_donor.parcel().mass_kg())
        + Dd::from_f64(next_receiver.parcel().mass_kg())
        - Dd::from_f64(donor.parcel().mass_kg());
    let mass_residual_kg = balance(
        "actual first-fill gas mass",
        mass_residual,
        mass_allowance_kg,
    )?;
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    let mut internal = Dd::ZERO;
    for source in donor.parcel().model().thermodynamics().components() {
        poll(&mut cancelled)?;
        let species = source.species();
        let before = inventory.get(species).copied().unwrap_or(0.0);
        let final_mass = after
            .each_ref()
            .map(|map| map.get(species).copied().unwrap_or(0.0));
        let limit = allowance(
            controls.mass_absolute_kg,
            controls.relative_tolerance,
            Dd::from_f64(received.get(species).copied().unwrap_or(0.0)),
        )?;
        let residual =
            Dd::from_f64(final_mass[0]) + Dd::from_f64(final_mass[1]) - Dd::from_f64(before);
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg
            .max(balance("actual first-fill species mass", residual, limit)?.abs());
        for (gas, amount) in [&next_donor, &next_receiver].into_iter().zip(final_mass) {
            if amount == 0.0 {
                continue;
            }
            let du = source
                .molar_internal_energy_change(
                    Temperature::new(reference),
                    Temperature::new(gas.parcel().state().state().temperature),
                )
                .map_err(|error| Error::Caloric {
                    species: species.clone(),
                    error,
                })?
                .value();
            internal = internal
                + Dd::from_f64(amount) * Dd::from_f64(du)
                    / Dd::from_f64(source.molar_mass().value());
        }
    }
    let mut momentum_residual_n_s = [0.0; 3];
    let mut momentum_allowance_n_s = [0.0; 3];
    let mut conversion = Dd::ZERO;
    for d in 0..3 {
        let residual = Dd::from_f64(next_donor.momentum_kg_m_s()[d])
            + Dd::from_f64(next_receiver.momentum_kg_m_s()[d])
            - Dd::from_f64(donor.momentum_kg_m_s()[d]);
        momentum_allowance_n_s[d] = allowance(
            controls.momentum_absolute_n_s,
            controls.relative_tolerance,
            impulse[d],
        )?;
        momentum_residual_n_s[d] = balance(
            "actual first-fill momentum",
            residual,
            momentum_allowance_n_s[d],
        )?;
        let v = Dd::from_f64(frame[d]);
        conversion = conversion + v * residual - (mass_residual * Dd::from_f64(0.5) * v) * v;
    }
    let residual = internal
        + Dd::from_f64(donor_step.kinetic_energy_change_j)
        + Dd::from_f64(receiver_step.kinetic_energy_j);
    let native = Dd::from_f64(donor_step.total_energy_residual_j)
        + Dd::from_f64(receiver_step.total_energy_residual_j);
    let defect = Dd::from_f64(donor_step.total_energy_residual_j).abs()
        + Dd::from_f64(receiver_step.absolute_energy_defect_j)
        + (residual - native).abs();
    energy_scale = energy_scale
        + Dd::from_f64(donor_step.kinetic_energy_change_j).abs()
        + Dd::from_f64(receiver_step.kinetic_energy_j).abs();
    let energy_allowance_j = allowance(
        controls.energy_absolute_j,
        controls.relative_tolerance,
        energy_scale,
    )?;
    let energy_residual_j = balance(
        "actual first-fill aperture-frame energy",
        residual,
        energy_allowance_j,
    )?;
    let laboratory_energy_residual_j = balance(
        "actual first-fill laboratory energy",
        residual + conversion,
        energy_allowance_j,
    )?;
    let absolute_energy_defect_j = balance(
        "absolute first-fill energy defects",
        defect + conversion.abs(),
        energy_allowance_j,
    )?;
    let step = MolecularFirstFillStep {
        flux,
        captured,
        donor: donor_step,
        receiver: receiver_step,
        captured_mass_kg,
        captured_momentum_n_s,
        reference_temperature_k: reference,
        captured_reference_energy_j,
        mass_residual_kg,
        maximum_species_mass_residual_kg,
        momentum_residual_n_s,
        energy_residual_j,
        laboratory_energy_residual_j,
        absolute_energy_defect_j,
        mass_allowance_kg,
        momentum_allowance_n_s,
        energy_allowance_j,
    };
    poll(&mut cancelled)?;
    Ok(MolecularFillCandidate {
        donor: next_donor,
        receiver: MolecularReceiver::Filled(next_receiver),
        step: MolecularFillStep {
            duration_s,
            surface,
            stage: MolecularFillStage::FirstFill(step),
        },
    })
}
