//! Species-selective streams between two original homogeneous real-gas states.
//!
//! Each species carries its upstream partial enthalpy, including the original
//! mixture-dependent PR departure. Both endpoint inventories and energies are
//! advanced from one snapshot. Prescribed streams supply no rate law, membrane
//! work or entropy-production guarantee; those belong to the transport owner.

/// Simultaneous selective streams and additional heat on real-gas graphs.
pub mod network;

use super::{
    GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, poll, positive, updated_point,
};
use crate::thermochemical_gas::energy::parcel::transfer::inventory;
use crate::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_qty::Temperature;
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::MAX_PENG_ROBINSON_COMPONENTS_V1;
use std::collections::BTreeMap;

/// The same signed species-stream controls and extensive accounts as ideal gas
/// exchange; source enthalpies here include their actual PR partial departures.
pub use crate::thermochemical_gas::energy::parcel::transfer::species::{
    GasSpeciesExchangeSolve, GasSpeciesExchangeStep, SpeciesEnergyTransfer, SpeciesMassTransfer,
};

/// Both reconstructed real-gas endpoints, or neither. Inputs remain immutable.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasSpeciesExchangeCandidate {
    left: IsochoricRealGasParcel,
    right: IsochoricRealGasParcel,
    step: GasSpeciesExchangeStep,
}
impl RealGasSpeciesExchangeCandidate {
    #[must_use]
    pub const fn left(&self) -> &IsochoricRealGasParcel {
        &self.left
    }
    #[must_use]
    pub const fn right(&self) -> &IsochoricRealGasParcel {
        &self.right
    }
    #[must_use]
    pub const fn step(&self) -> &GasSpeciesExchangeStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        IsochoricRealGasParcel,
        IsochoricRealGasParcel,
        GasSpeciesExchangeStep,
    ) {
        (self.left, self.right, self.step)
    }
}

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), GasEnergyError> {
        let adjusted = term - self.correction;
        let value = self.value + adjusted;
        let correction = (value - self.value) - adjusted;
        if ![term, adjusted, value, correction]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "real-gas selective exchange arithmetic is unrepresentable",
            ));
        }
        self.value = value;
        self.correction = correction;
        Ok(())
    }
}

// A balanced counterstream changes composition but not the accepted density.
// Avoid an unnecessary rho*V/V round trip when the extensive mass is unchanged.
fn endpoint_density(old: &IsochoricRealGasParcel, mass: f64) -> f64 {
    if mass == old.mass_kg() {
        old.density_kg_m3()
    } else {
        mass / old.volume_m3()
    }
}

// The new inventory at its final density and old temperature owns the anchor.
// NASA formation datums have already cancelled species by species before this
// original fixed-composition finite caloric solve is entered.
fn endpoint(
    old: &IsochoricRealGasParcel,
    model: ThermochemicalPengRobinsonModel,
    mass: f64,
    target_j: f64,
    controls: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(IsochoricRealGasParcel, f64, usize), GasEnergyError> {
    let density = endpoint_density(old, mass);
    let temperature = old.state().temperature_k();
    let change = target_j / mass;
    if !change.is_finite() || (target_j != 0.0 && change == 0.0) {
        return Err(invalid(
            "real-gas selective caloric increment is unrepresentable",
        ));
    }
    // As in original ideal-gas network exchange, admit the closest temperature
    // only below half of one original-law temperature ULP and inside the user's
    // energy budget. The unmodified target still determines the reported defect.
    let mut represented_change = change;
    if change != 0.0 && change.abs() <= controls.energy_tolerance_j_kg {
        let neighbor = if change > 0.0 {
            temperature.next_up()
        } else {
            temperature.next_down()
        };
        if (controls.minimum_temperature_k..=controls.maximum_temperature_k).contains(&neighbor) {
            let adjacent =
                model.specific_internal_energy_change_at_density(density, temperature, neighbor)?;
            if !adjacent.is_finite()
                || adjacent == 0.0
                || adjacent.is_sign_positive() != change.is_sign_positive()
            {
                return Err(invalid(
                    "real-gas selective neighboring caloric increment is unrepresentable",
                ));
            }
            if change.abs() < 0.5 * adjacent.abs() {
                represented_change = 0.0;
            }
        }
    }
    let solved = model.reconstruct_density_internal_energy_increment(
        old.point(),
        density,
        represented_change,
        controls,
        &mut *cancelled,
    )?;
    let iterations = solved.iterations();
    let state = solved.into_state();
    let stored = mass
        * model.specific_internal_energy_change_at_density(
            density,
            temperature,
            state.temperature_k(),
        )?;
    let residual = stored - target_j;
    if !residual.is_finite() || (residual / mass).abs() > controls.energy_tolerance_j_kg {
        return Err(GasEnergyError::NotConverged {
            iterations,
            residual_j_kg: residual / mass,
        });
    }
    let point = updated_point(old.point(), state.temperature_k(), state.pressure_pa())?;
    let mut next = old.clone();
    next.model = model;
    next.state = state;
    next.point = point;
    next.mass_kg = mass;
    next.density_kg_m3 = density;
    Ok((next, residual, iterations))
}

impl IsochoricRealGasParcel {
    /// Apply signed species streams to the SAME pair of fixed-volume states.
    /// Positive mass moves left to right; negative mass moves right to left.
    /// Each stream carries its original upstream partial enthalpy, including
    /// flow work once. The simultaneous network owner supplies the same source
    /// joins, reference-safe caloric target, endpoint solves and balance audits.
    ///
    /// The pair interface retains unique SpeciesId rows and its strict local
    /// increment representability checks. No accepted input is mutated. These
    /// prescribed streams imply neither rates, separation work nor entropy
    /// production; the physical transport owner must supply them separately.
    pub fn exchange_species_with(
        &self,
        right: &Self,
        transfers: &[SpeciesMassTransfer],
        controls: GasSpeciesExchangeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasSpeciesExchangeCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if transfers.len() > MAX_PENG_ROBINSON_COMPONENTS_V1 {
            return Err(invalid(
                "real-gas selective exchange exceeds the original species bound",
            ));
        }
        let old = [self, right];
        let mut ordered: BTreeMap<SpeciesId, f64> = BTreeMap::new();
        for row in transfers {
            poll(&mut cancelled)?;
            if !row.mass_left_to_right_kg.is_finite()
                || ordered
                    .insert(row.species.clone(), row.mass_left_to_right_kg)
                    .is_some()
                || !old.iter().any(|gas| {
                    gas.model()
                        .source()
                        .thermodynamics()
                        .components()
                        .iter()
                        .any(|species| species.species() == &row.species)
                })
            {
                return Err(invalid(
                    "real-gas selective streams need finite unique masses and original species",
                ));
            }
        }
        if ordered.values().all(|dm| *dm == 0.0) {
            poll(&mut cancelled)?;
            return Ok(RealGasSpeciesExchangeCandidate {
                left: self.clone(),
                right: right.clone(),
                step: GasSpeciesExchangeStep {
                    species: Vec::new(),
                    net_mass_left_to_right_kg: 0.0,
                    enthalpy_left_to_right_j: 0.0,
                    energy_change_j: [0.0; 2],
                    mass_residual_kg: 0.0,
                    maximum_species_mass_residual_kg: 0.0,
                    energy_residual_j: 0.0,
                    iterations: 0,
                },
            });
        }
        if core::ptr::eq(self, right) {
            return Err(invalid(
                "nonzero real-gas selective exchange needs two distinct physical owners",
            ));
        }
        let before = [
            inventory(self.model().source(), self.mass_kg())?,
            inventory(right.model().source(), right.mass_kg())?,
        ];
        let mut mass = Sum::default();
        let mut rows = Vec::with_capacity(ordered.len());
        for (species, &dm) in &ordered {
            poll(&mut cancelled)?;
            if dm == 0.0 {
                continue;
            }
            let donor = usize::from(dm < 0.0);
            let available = before[donor].get(species).copied().unwrap_or(0.0);
            if !positive(available) || dm.abs() / available > controls.maximum_species_fraction {
                return Err(invalid(
                    "real-gas selective stream exceeds original species withdrawal allowance",
                ));
            }
            for (side, sign) in [(0, -1.0), (1, 1.0)] {
                let previous = before[side].get(species).copied().unwrap_or(0.0);
                let next = previous + sign * dm;
                if !positive(next)
                    || next == previous
                    || ((next - previous) - sign * dm).abs() > controls.mass_tolerance_kg
                {
                    return Err(invalid(
                        "real-gas selective species increment vanished or exhausted its inventory",
                    ));
                }
            }
            mass.add(dm)?;
            rows.push(network::SpeciesNetworkTransfer {
                donor: donor as u64,
                receiver: (1 - donor) as u64,
                species: species.clone(),
                mass_kg: dm.abs(),
            });
        }
        for (side, sign) in [(0, -1.0), (1, 1.0)] {
            let next = old[side].mass_kg() + sign * mass.value;
            let density = endpoint_density(old[side], next);
            if !positive(next)
                || !positive(density)
                || (mass.value != 0.0
                    && (next == old[side].mass_kg() || density == old[side].density_kg_m3()))
                || ((next - old[side].mass_kg()) - sign * mass.value).abs()
                    > controls.mass_tolerance_kg
            {
                return Err(invalid(
                    "real-gas selective total mass or density is unrepresentable",
                ));
            }
        }
        let candidate = network::exchange_borrowed(
            &BTreeMap::from([(0, self), (1, right)]),
            &rows,
            &[],
            &[],
            &BTreeMap::from([(0, controls.gas[0]), (1, controls.gas[1])]),
            network::SpeciesNetworkSolve {
                maximum_species_fraction: controls.maximum_species_fraction,
                mass_tolerance_kg: controls.mass_tolerance_kg,
                energy_tolerance_j: controls.energy_tolerance_j,
            },
            &mut cancelled,
        )
        .map_err(|error| match error {
            network::SpeciesNetworkError::Cancelled => GasEnergyError::Cancelled,
            network::SpeciesNetworkError::Invalid(message) => invalid(message),
            network::SpeciesNetworkError::Parcel { error, .. } => error,
            network::SpeciesNetworkError::Withdrawal { .. } => {
                invalid("real-gas selective stream exceeds original species withdrawal allowance")
            }
        })?;
        let (mut parcels, step) = candidate.into_parts();
        let mut species: Vec<_> = step
            .fluxes
            .iter()
            .map(|flux| {
                let sign = if flux.transfer.donor == 0 { 1.0 } else { -1.0 };
                SpeciesEnergyTransfer {
                    transfer: SpeciesMassTransfer {
                        species: flux.transfer.species.clone(),
                        mass_left_to_right_kg: sign * flux.transfer.mass_kg,
                    },
                    source_enthalpy_j_kg: flux.source_enthalpy_j_kg,
                    enthalpy_left_to_right_j: sign * flux.carried_enthalpy_j,
                }
            })
            .collect();
        // Preserve the pair's original SpeciesId order, regardless of direction.
        species.sort_by(|a, b| a.transfer.species.cmp(&b.transfer.species));
        let mut carried = Sum::default();
        for row in &species {
            carried.add(row.enthalpy_left_to_right_j)?;
        }
        let residual = [
            step.parcels[0].energy_residual_j,
            step.parcels[1].energy_residual_j,
        ];
        poll(&mut cancelled)?;
        Ok(RealGasSpeciesExchangeCandidate {
            left: parcels.remove(&0).expect("admitted pair endpoint"),
            right: parcels.remove(&1).expect("admitted pair endpoint"),
            step: GasSpeciesExchangeStep {
                species,
                net_mass_left_to_right_kg: mass.value,
                enthalpy_left_to_right_j: carried.value,
                energy_change_j: [-carried.value + residual[0], carried.value + residual[1]],
                mass_residual_kg: step.total_mass_residual_kg,
                maximum_species_mass_residual_kg: step.maximum_species_mass_residual_kg,
                energy_residual_j: step.total_energy_residual_j,
                iterations: step.parcels.iter().map(|p| p.iterations).sum(),
            },
        })
    }
}

#[cfg(test)]
mod tests;
