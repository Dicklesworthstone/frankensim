//! Species-selective streams between two original homogeneous real-gas states.
//!
//! Each species carries its upstream partial enthalpy, including the original
//! mixture-dependent PR departure. Both endpoint inventories and energies are
//! advanced from one snapshot. Prescribed streams supply no rate law, membrane
//! work or entropy-production guarantee; those belong to the transport owner.

use super::{
    GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, mixing, poll, positive, updated_point,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::model_for_inventory;
use crate::thermochemical_gas::energy::parcel::transfer::{inventory, mixed_model_from_portions};
use crate::thermochemical_gas::{failure, nonideal::ThermochemicalPengRobinsonModel};
use fs_qty::Temperature;
use fs_thermochem::eos::peng_robinson::MAX_PENG_ROBINSON_COMPONENTS_V1;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
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
    /// Apply all signed species streams to the SAME pair of fixed-volume states.
    ///
    /// Positive mass moves left to right; negative mass moves right to left.
    /// Counterflow with zero net mass can still exchange substantial enthalpy.
    /// Rows are canonicalized by SpeciesId. Each nonzero stream must remain
    /// representable and below its original donor's species withdrawal bound.
    /// NASA/transport/pack/critical sources and every PR pair entry must agree;
    /// a missing cross-species interaction is never inferred as zero.
    ///
    /// A species carries its actual upstream partial h_i, with flow work once.
    /// Each new composition at its final density and old temperature supplies a
    /// caloric anchor. Conserved ideal species datums cancel analytically before
    /// evaluating NASA Delta h, RT and PR departures. The original finite PR
    /// caloric owner reconstructs both final temperatures and complete states.
    /// Actual inventory roundoff retains its separate absolute mass budget.
    ///
    /// This first-order stream primitive does not select rates, add separation
    /// work, infer phase equilibrium or assert entropy production for arbitrary
    /// prescribed fluxes. Volumes and extra material axes remain unchanged.
    /// Cancellation or any source/conservation refusal returns neither endpoint.
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
        let mut ordered = BTreeMap::new();
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
        let mut target = before.clone();
        let mut mass = Sum::default();
        for (id, &dm) in &ordered {
            poll(&mut cancelled)?;
            if dm == 0.0 {
                continue;
            }
            let donor = usize::from(dm < 0.0);
            let available = before[donor].get(id).copied().unwrap_or(0.0);
            if !positive(available) || dm.abs() / available > controls.maximum_species_fraction {
                return Err(invalid(
                    "real-gas selective stream exceeds original species withdrawal allowance",
                ));
            }
            for (side, sign) in [(0, -1.0), (1, 1.0)] {
                let previous = before[side].get(id).copied().unwrap_or(0.0);
                let next = previous + sign * dm;
                if !positive(next)
                    || next == previous
                    || ((next - previous) - sign * dm).abs() > controls.mass_tolerance_kg
                {
                    return Err(invalid(
                        "real-gas selective species increment vanished or exhausted its inventory",
                    ));
                }
                target[side].insert(id.clone(), next);
            }
            mass.add(dm)?;
        }
        let masses = [self.mass_kg() - mass.value, right.mass_kg() + mass.value];
        for (side, sign) in [(0, -1.0), (1, 1.0)] {
            let density = endpoint_density(old[side], masses[side]);
            if !positive(masses[side])
                || !positive(density)
                || (mass.value != 0.0
                    && (masses[side] == old[side].mass_kg()
                        || density == old[side].density_kg_m3()))
                || ((masses[side] - old[side].mass_kg()) - sign * mass.value).abs()
                    > controls.mass_tolerance_kg
            {
                return Err(invalid(
                    "real-gas selective total mass or density is unrepresentable",
                ));
            }
        }
        let catalog = mixed_model_from_portions(
            &[
                (self.model().source(), self.mass_kg()),
                (right.model().source(), right.mass_kg()),
            ],
            self.mass_kg() + right.mass_kg(),
        )?;
        let originals = [self.model(), right.model()];
        let models = [
            mixing::model_for_source(&originals, model_for_inventory(&catalog, &target[0])?)?,
            mixing::model_for_source(&originals, model_for_inventory(&catalog, &target[1])?)?,
        ];
        let partial = [
            self.model()
                .partial_calorics_at_density(self.density_kg_m3(), self.state().temperature_k())?,
            right.model().partial_calorics_at_density(
                right.density_kg_m3(),
                right.state().temperature_k(),
            )?,
        ];
        let mut caloric = [Sum::default(), Sum::default()];
        let mut molar_flow = [Sum::default(), Sum::default()];
        for (side, sign) in [(0, -1.0), (1, 1.0)] {
            let anchor = models[side]
                .caloric_at_density(
                    endpoint_density(old[side], masses[side]),
                    old[side].state().temperature_k(),
                )?
                .residual_energy_j_kg;
            let departure = old[side].state().root().residual_internal_energy_j_mol()
                / old[side].model().mechanics().molar_mass().value();
            // Separate departure change from extensive mass change before adding
            // streams; absolute ideal energy references never enter this sum.
            caloric[side].add(old[side].mass_kg() * (departure - anchor))?;
            caloric[side].add(-sign * mass.value * anchor)?;
        }
        let mut carried = Sum::default();
        let mut species = Vec::with_capacity(ordered.len());
        for (id, &dm) in &ordered {
            poll(&mut cancelled)?;
            if dm == 0.0 {
                continue;
            }
            let donor = usize::from(dm < 0.0);
            let (index, source) = old[donor]
                .model()
                .source()
                .thermodynamics()
                .components()
                .iter()
                .enumerate()
                .find(|(_, species)| species.species() == id)
                .ok_or(invalid(
                    "real-gas selective species absent from its upstream owner",
                ))?;
            let molar_mass = source.molar_mass().value();
            let residual_h =
                partial[donor].departure().components()[index].residual_enthalpy_j_mol();
            let h = partial[donor].partial_specific_enthalpies_j_kg()[index];
            let stream_h = dm * h;
            carried.add(stream_h)?;
            for (side, sign) in [(0, -1.0), (1, 1.0)] {
                let sensible_h = source
                    .molar_enthalpy_change(
                        Temperature::new(old[side].state().temperature_k()),
                        Temperature::new(old[donor].state().temperature_k()),
                    )
                    .map_err(|error| failure("real-gas selective caloric increment", error))?
                    .value();
                let amount = sign * dm / molar_mass;
                if amount == 0.0 {
                    return Err(invalid(
                        "real-gas selective molar stream is unrepresentable",
                    ));
                }
                molar_flow[side].add(amount)?;
                caloric[side].add(amount * sensible_h)?;
                caloric[side].add(amount * residual_h)?;
            }
            species.push(SpeciesEnergyTransfer {
                transfer: SpeciesMassTransfer {
                    species: id.clone(),
                    mass_left_to_right_kg: dm,
                },
                source_enthalpy_j_kg: h,
                enthalpy_left_to_right_j: stream_h,
            });
        }
        for side in 0..2 {
            // h_i(Tup)-u_i(T0) = Delta h_i(T0,Tup)+R*T0. Summing molar
            // counterflow before multiplying RT retains balanced-stream heat.
            caloric[side].add(R * old[side].state().temperature_k() * molar_flow[side].value)?;
        }
        let [left_model, right_model] = models;
        let (left, left_residual, left_iterations) = endpoint(
            self,
            left_model,
            masses[0],
            caloric[0].value,
            controls.gas[0],
            &mut cancelled,
        )?;
        let (right_next, right_residual, right_iterations) = endpoint(
            right,
            right_model,
            masses[1],
            caloric[1].value,
            controls.gas[1],
            &mut cancelled,
        )?;
        let after = [
            inventory(left.model().source(), left.mass_kg())?,
            inventory(right_next.model().source(), right_next.mass_kg())?,
        ];
        let mut maximum_species_mass_residual_kg = 0.0_f64;
        for source in catalog.thermodynamics().components() {
            let id: &SpeciesId = source.species();
            let dm = ordered.get(id).copied().unwrap_or(0.0);
            let changes = core::array::from_fn::<_, 2, _>(|side| {
                after[side].get(id).copied().unwrap_or(0.0)
                    - before[side].get(id).copied().unwrap_or(0.0)
            });
            for defect in [changes[0] + dm, changes[1] - dm, changes[0] + changes[1]] {
                if !defect.is_finite() {
                    return Err(invalid(
                        "real-gas selective species balance is unrepresentable",
                    ));
                }
                maximum_species_mass_residual_kg =
                    maximum_species_mass_residual_kg.max(defect.abs());
            }
        }
        let energy_residual_j = left_residual + right_residual;
        let mass_residual_kg =
            (left.mass_kg() - self.mass_kg()) + (right_next.mass_kg() - right.mass_kg());
        let energy_change_j = [
            -carried.value + left_residual,
            carried.value + right_residual,
        ];
        if ![
            energy_residual_j,
            mass_residual_kg,
            energy_change_j[0],
            energy_change_j[1],
        ]
        .into_iter()
        .all(f64::is_finite)
            || [left_residual, right_residual, energy_residual_j]
                .into_iter()
                .any(|error| error.abs() > controls.energy_tolerance_j)
            || mass_residual_kg.abs() > controls.mass_tolerance_kg
            || maximum_species_mass_residual_kg > controls.mass_tolerance_kg
        {
            return Err(invalid(
                "real-gas selective exchange did not close local and paired mass/energy",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RealGasSpeciesExchangeCandidate {
            left,
            right: right_next,
            step: GasSpeciesExchangeStep {
                species,
                net_mass_left_to_right_kg: mass.value,
                enthalpy_left_to_right_j: carried.value,
                energy_change_j,
                mass_residual_kg,
                maximum_species_mass_residual_kg,
                energy_residual_j,
                iterations: left_iterations + right_iterations,
            },
        })
    }
}

#[cfg(test)]
mod tests;
