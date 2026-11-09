//! Prescribed quasineutral fueling and selective thermal-particle exhaust.
//!
//! These are local inventory sources, not hydrodynamic inlet/outlet fluxes.
//! Injected fully ionized particles bring their original sensible internal
//! energy and mix into the two Maxwellian pools. Departing thermal particles
//! carry the current mean internal energy of their pool. Companion electrons
//! enter/leave with each ion's charge; no electrostatic work is inferred.
//! Independent energy-confinement losses must exclude this carried energy.

use super::solve::Sum;
use super::{DtFusionError, R, poll, positive};
use fs_math::decay::{attenuate, decayed_fraction};
use fs_qty::{Frequency, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::fusion::SPEED_OF_LIGHT_M_S;
use fs_thermochem::fusion::plasma::{DtFusionPlasmaErrorV1, DtFusionPlasmaLimitsV1};

const CHARGE: [f64; 5] = [1.0, 1.0, 1.0, 2.0, 2.0];
pub(super) const BARYONS: [f64; 5] = [1.0, 2.0, 3.0, 3.0, 4.0];

/// Constant fully-ionized supply and first-order thermal-particle removal.
/// All ion arrays use physical order p, D, T, He3, He4, independently of the
/// source bank's canonical order. Zero removal means retained particles.
///
/// Rates and incoming temperatures are explicitly authored. No pellet ablation,
/// neutral ionization, recycling, momentum, ambipolar potential or magnetic
/// confinement law is supplied. In a kinetic-alpha solve, He4 removal acts only
/// on the Maxwellian ash; the separate alpha distribution is untouched.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionParticleExchange {
    injection_mol_s: [f64; 5],
    removal_s_inv: [f64; 5],
    incoming_temperatures: [Temperature; 2],
}

impl FusionParticleExchange {
    pub fn new(
        injection_mol_s: [f64; 5],
        removal_hazards: [Frequency; 5],
        incoming_heavy_temperature: Temperature,
        incoming_electron_temperature: Temperature,
    ) -> Result<Self, DtFusionError> {
        let removal_s_inv = removal_hazards.map(|rate| rate.value());
        let incoming_temperatures = [incoming_heavy_temperature, incoming_electron_temperature];
        if injection_mol_s
            .iter()
            .chain(removal_s_inv.iter())
            .any(|rate| !rate.is_finite() || *rate < 0.0)
            || incoming_temperatures.iter().any(|t| !positive(t.value()))
        {
            return Err(DtFusionError::Invalid(
                "finite nonnegative particle rates and positive incoming temperatures required",
            ));
        }
        Ok(Self {
            injection_mol_s,
            removal_s_inv,
            incoming_temperatures,
        })
    }
    pub const fn injection_rates_mol_s(&self) -> [f64; 5] {
        self.injection_mol_s
    }
    pub fn removal_hazards(&self) -> [Frequency; 5] {
        self.removal_s_inv.map(Frequency::new)
    }
    pub const fn incoming_temperatures(&self) -> [Temperature; 2] {
        self.incoming_temperatures
    }
}

/// Actual particle transfers and their thermal energy over the accepted solve.
/// Ion arrays use p, D, T, He3, He4 order. Energy arrays are heavy, electron.
/// The reservoir itself is prescribed; these exports are available for an
/// explicitly coupled receiver, without inventing its state or caloric law.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FusionParticleBalance {
    pub injected_ions_mol: [f64; 5],
    pub exhausted_ions_mol: [f64; 5],
    pub injected_electrons_mol: f64,
    pub exhausted_electrons_mol: f64,
    pub injected_thermal_energy_j: [f64; 2],
    pub exhausted_thermal_energy_j: [f64; 2],
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Accounts {
    pub(super) injected: [Sum; 5],
    pub(super) exhausted: [Sum; 5],
    pub(super) electrons_in: Sum,
    pub(super) electrons_out: Sum,
    pub(super) energy_in: [Sum; 2],
    pub(super) energy_out: [Sum; 2],
}
impl Accounts {
    pub(super) fn balance(self) -> FusionParticleBalance {
        FusionParticleBalance {
            injected_ions_mol: self.injected.map(|n| n.value),
            exhausted_ions_mol: self.exhausted.map(|n| n.value),
            injected_electrons_mol: self.electrons_in.value,
            exhausted_electrons_mol: self.electrons_out.value,
            injected_thermal_energy_j: self.energy_in.map(|u| u.value),
            exhausted_thermal_energy_j: self.energy_out.map(|u| u.value),
        }
    }
    pub(super) fn is_zero(self) -> bool {
        self.injected
            .iter()
            .chain(&self.exhausted)
            .all(|n| n.value == 0.0)
    }
    pub(super) fn finite(self) -> bool {
        self.injected
            .iter()
            .chain(&self.exhausted)
            .chain(&self.energy_in)
            .chain(&self.energy_out)
            .chain([&self.electrons_in, &self.electrons_out])
            .all(|x| x.value.is_finite() && x.value >= 0.0)
    }
    pub(super) fn nuclei_net_in(self) -> [f64; 2] {
        let mut baryons = Sum::default();
        let mut charge = Sum::default();
        for i in 0..5 {
            baryons.add(BARYONS[i] * self.injected[i].value);
            baryons.add(-BARYONS[i] * self.exhausted[i].value);
            charge.add(CHARGE[i] * self.injected[i].value);
            charge.add(-CHARGE[i] * self.exhausted[i].value);
        }
        [baryons.value, charge.value]
    }
    pub(super) fn nuclei_activity(self) -> [f64; 2] {
        let mut result = [Sum::default(); 2];
        for i in 0..5 {
            let amount = self.injected[i].value + self.exhausted[i].value;
            result[0].add(BARYONS[i] * amount);
            result[1].add(CHARGE[i] * amount);
        }
        result.map(|n| n.value)
    }
    pub(super) fn compare(
        self,
        full: Self,
        fine: Self,
        amount_tolerance: f64,
        energy_tolerance: f64,
        mut add: impl FnMut(f64, f64, f64),
    ) {
        for i in 0..5 {
            add(
                full.injected[i].increment_since(self.injected[i]),
                fine.injected[i].increment_since(self.injected[i]),
                amount_tolerance,
            );
            add(
                full.exhausted[i].increment_since(self.exhausted[i]),
                fine.exhausted[i].increment_since(self.exhausted[i]),
                amount_tolerance,
            );
        }
        for i in 0..2 {
            add(
                full.energy_in[i].increment_since(self.energy_in[i]),
                fine.energy_in[i].increment_since(self.energy_in[i]),
                energy_tolerance,
            );
            add(
                full.energy_out[i].increment_since(self.energy_out[i]),
                fine.energy_out[i].increment_since(self.energy_out[i]),
                energy_tolerance,
            );
        }
    }
}

/// Source-resolved immutable preparation, once per burn rather than per stage.
pub(super) struct Prepared {
    law: FusionParticleExchange,
    indices: [Option<usize>; 5],
    electron: usize,
    incoming_u_j_mol: [f64; 6], // physical ions, then electron
}
impl Prepared {
    pub(super) fn new(
        law: FusionParticleExchange,
        bank: &IdealGasEquilibriumV1,
        indices: [Option<usize>; 5],
        electron: usize,
        limits: DtFusionPlasmaLimitsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        let mut incoming_u_j_mol = [0.0; 6];
        for slot in 0..5 {
            poll(cancelled)?;
            if indices[slot].is_none()
                && (law.injection_mol_s[slot] != 0.0 || law.removal_s_inv[slot] != 0.0)
            {
                return Err(DtFusionError::Invalid(
                    "particle exchange requires the original supported ion species",
                ));
            }
            if law.injection_mol_s[slot] != 0.0 {
                let source = &bank.models()[indices[slot].expect("admitted ion")];
                if R * law.incoming_temperatures[0].value()
                    / (source.molar_mass().value() * SPEED_OF_LIGHT_M_S.powi(2))
                    > limits.maximum_thermal_energy_rest_energy_ratio
                {
                    return Err(DtFusionError::OutsideDomain);
                }
                incoming_u_j_mol[slot] = source
                    .evaluate_properties(law.incoming_temperatures[0])
                    .map_err(|error| {
                        DtFusionError::Plasma(DtFusionPlasmaErrorV1::Source {
                            species: source.species().clone(),
                            source: error,
                        })
                    })?
                    .u()
                    .value();
            }
        }
        if law.injection_mol_s.iter().any(|rate| *rate != 0.0) {
            let source = &bank.models()[electron];
            if R * law.incoming_temperatures[1].value()
                / (source.molar_mass().value() * SPEED_OF_LIGHT_M_S.powi(2))
                > limits.maximum_thermal_energy_rest_energy_ratio
            {
                return Err(DtFusionError::OutsideDomain);
            }
            incoming_u_j_mol[5] = source
                .evaluate_properties(law.incoming_temperatures[1])
                .map_err(|error| {
                    DtFusionError::Plasma(DtFusionPlasmaErrorV1::Source {
                        species: source.species().clone(),
                        source: error,
                    })
                })?
                .u()
                .value();
        }
        if incoming_u_j_mol.iter().any(|u| !u.is_finite() || *u < 0.0) {
            return Err(DtFusionError::OutsideDomain);
        }
        Ok(Self {
            law,
            indices,
            electron,
            incoming_u_j_mol,
        })
    }

    /// Symmetric exact exhaust/injection maps. Their composition is second
    /// order for simultaneous unlike species hazards and thermal mixing; the
    /// original burn's step doubling resolves that coupling error as well.
    /// Driven exhaust must first resolve its relaxation time: saturated coarse
    /// and fine maps can agree while erasing the physical supplied stock S/lambda.
    /// Pure exhaust retains its unrestricted exact map.
    pub(super) fn apply<const N: usize>(
        &self,
        amounts: &mut [Sum; N],
        thermal_he4: Option<&mut Sum>,
        kinetic_he4: f64,
        energy: &mut [Sum; 2],
        accounts: &mut Accounts,
        dt: f64,
        other_sources: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DtFusionError> {
        poll(cancelled)?;
        if dt == 0.0 {
            return Ok(());
        }
        if other_sources || self.law.injection_mol_s.iter().any(|rate| *rate > 0.0) {
            for slot in 0..5 {
                let Some(index) = self.indices[slot] else {
                    continue;
                };
                let available = if slot == 4 {
                    thermal_he4
                        .as_ref()
                        .map_or(amounts[index].value, |n| n.value)
                } else {
                    amounts[index].value
                };
                if (available > 0.0 || self.law.injection_mol_s[slot] > 0.0)
                    && self.law.removal_s_inv[slot] * dt > 1.0
                {
                    return Err(DtFusionError::StepBudget);
                }
            }
        }
        let mut ash = thermal_he4.as_ref().map(|n| **n);
        self.exhaust(amounts, &mut ash, kinetic_he4, energy, accounts, 0.5 * dt)?;
        let mut supplied_electrons = Sum::default();
        let mut supplied_energy = [Sum::default(); 2];
        for slot in 0..5 {
            poll(cancelled)?;
            let supplied = self.law.injection_mol_s[slot] * dt;
            if !supplied.is_finite() || (self.law.injection_mol_s[slot] > 0.0 && supplied == 0.0) {
                return Err(DtFusionError::Arithmetic("particle injection amount"));
            }
            if supplied == 0.0 {
                continue;
            }
            amounts[self.indices[slot].expect("admitted ion")].add(supplied);
            if slot == 4 {
                if let Some(thermal) = &mut ash {
                    thermal.add(supplied);
                }
            }
            accounts.injected[slot].add(supplied);
            supplied_electrons.add(CHARGE[slot] * supplied);
            supplied_energy[0].add(self.incoming_u_j_mol[slot] * supplied);
        }
        amounts[self.electron].add(supplied_electrons.value);
        accounts.electrons_in.add(supplied_electrons.value);
        supplied_energy[1].add(self.incoming_u_j_mol[5] * supplied_electrons.value);
        for pool in 0..2 {
            energy[pool].add(supplied_energy[pool].value);
            accounts.energy_in[pool].add(supplied_energy[pool].value);
        }
        self.exhaust(amounts, &mut ash, kinetic_he4, energy, accounts, 0.5 * dt)?;
        if let (Some(target), Some(value)) = (thermal_he4, ash) {
            *target = value;
        }
        if amounts
            .iter()
            .chain(energy.iter())
            .any(|x| !x.value.is_finite() || x.value < 0.0)
            || !accounts.finite()
        {
            return Err(DtFusionError::Arithmetic(
                "finite particle inventory and thermal accounts",
            ));
        }
        poll(cancelled)
    }

    fn exhaust<const N: usize>(
        &self,
        amounts: &mut [Sum; N],
        thermal_he4: &mut Option<Sum>,
        kinetic_he4: f64,
        energy: &mut [Sum; 2],
        accounts: &mut Accounts,
        dt: f64,
    ) -> Result<(), DtFusionError> {
        let mut heavy = [Sum::default(); 2]; // original and surviving thermal ions
        let mut charge_remaining = Sum::new(2.0 * kinetic_he4);
        let mut charge_error = Sum::new(amounts[self.electron].value);
        let mut removed = [Sum::default(); 2];
        for slot in 0..5 {
            let Some(index) = self.indices[slot] else {
                continue;
            };
            charge_error.add(-CHARGE[slot] * amounts[index].value);
            let old = if slot == 4 {
                thermal_he4.unwrap_or(amounts[index])
            } else {
                amounts[index]
            };
            let exposure = self.law.removal_s_inv[slot] * dt;
            let loss = old.value * decayed_fraction(exposure);
            let survivor = attenuate(old.value, exposure);
            if !loss.is_finite()
                || !survivor.is_finite()
                || (old.value > 0.0 && self.law.removal_s_inv[slot] > 0.0 && exposure == 0.0)
            {
                return Err(DtFusionError::Arithmetic(
                    "thermal-particle removal exposure",
                ));
            }
            let remaining = withdraw(old, loss, survivor);
            if slot == 4 && thermal_he4.is_some() {
                *thermal_he4 = Some(remaining);
                amounts[index] = withdraw(amounts[index], loss, survivor + kinetic_he4);
            } else {
                amounts[index] = remaining;
            }
            heavy[0].add(old.value);
            heavy[1].add(remaining.value);
            charge_remaining.add(CHARGE[slot] * remaining.value);
            removed[0].add(loss);
            removed[1].add(CHARGE[slot] * loss);
            accounts.exhausted[slot].add(loss);
        }
        charge_remaining.add(charge_error.value); // retain the admitted initial charge discrepancy
        let electrons = amounts[self.electron].value;
        amounts[self.electron] = withdraw(
            amounts[self.electron],
            removed[1].value,
            charge_remaining.value,
        );
        accounts.electrons_out.add(removed[1].value);
        for (pool, original, remaining) in [
            (0, heavy[0].value, heavy[1].value),
            (1, electrons, amounts[self.electron].value),
        ] {
            if removed[pool].value == 0.0 {
                continue;
            }
            if !positive(original) || remaining < 0.0 {
                return Err(DtFusionError::OutsideDomain);
            }
            let per_particle = energy[pool].value / original;
            let outgoing = per_particle * removed[pool].value;
            energy[pool] = withdraw(energy[pool], outgoing, per_particle * remaining);
            accounts.energy_out[pool].add(outgoing);
        }
        Ok(())
    }
}

pub(super) fn withdraw(old: Sum, removed: f64, survivor: f64) -> Sum {
    if removed <= 0.5 * old.value {
        let mut next = old;
        next.add(-removed);
        next
    } else {
        Sum::new(survivor)
    }
}
