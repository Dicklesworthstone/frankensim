//! Simultaneous molecular apertures on ORIGINAL finite moving-gas inventories.
//!
//! All links sample one snapshot. Gross withdrawals pool before incoming mass
//! can be credited, and each gas receives one original caloric/momentum update.
//! Apertures may move at different prescribed velocities; crossing energies are
//! transformed into one explicitly declared inertial frame before pooling.
//! This is a first-order fixed-volume interface operator. Container motion,
//! flight time, remaining wall forces and closed-system passivity are not inferred.

/// Adaptive histories of the same simultaneous finite-inventory transaction.
pub mod trajectory;

/// Finite material-body heat feeding the same moving-gas transaction.
pub mod thermal;

use super::{
    DriftingEffusiveExchange, EffusiveDriftInterface, GasEnergyError, MAX_EFFUSIVE_SPECIES,
    MolecularExchangeError, MolecularExchangeSurface, MovingGasBoundaryError,
    MovingGasBoundarySolve, MovingGasBoundaryStep, MovingGasParcel,
    SpeciesBoundaryRelativeTransfer, SpeciesId, allowance, balance, exchange, finite, positive,
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
};
use fs_math::dd::Dd;
use std::collections::{BTreeMap, BTreeSet};

/// Separately exposed aperture, oriented from `gases[0]` to `gases[1]`.
/// Active links incident on one gas must retain the SAME original catalog.
#[derive(Clone, Copy, Debug)]
pub struct MolecularNetworkLink<'a> {
    pub id: u64,
    pub gases: [u64; 2],
    pub interface: &'a EffusiveDriftInterface,
    pub surface: MolecularExchangeSurface,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularNetworkControls {
    /// Common energy frame; velocities/momenta still use the original inertial axes.
    pub energy_frame_velocity_m_s: [f64; 3],
    pub maximum_duration_s: f64,
    pub maximum_species_fraction: f64,
    pub mass_absolute_kg: f64,
    pub momentum_absolute_n_s: f64,
    pub energy_absolute_j: f64,
    /// Relative accounts scale exchanged quantities, never original inventories.
    pub relative_tolerance: f64,
}

#[derive(Debug)]
pub enum MolecularNetworkError {
    Invalid(&'static str),
    Link {
        id: u64,
        error: GasEnergyError,
    },
    Gas {
        id: u64,
        error: MovingGasBoundaryError,
    },
    Withdrawal {
        gas: u64,
        species: SpeciesId,
        fraction: f64,
        maximum: f64,
    },
    Accounting(MolecularExchangeError),
    Cancelled,
}
impl core::fmt::Display for MolecularNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous molecular network refused: {self:?}")
    }
}
impl std::error::Error for MolecularNetworkError {}
impl From<MolecularExchangeError> for MolecularNetworkError {
    fn from(error: MolecularExchangeError) -> Self {
        match error {
            MolecularExchangeError::Cancelled => Self::Cancelled,
            error => Self::Accounting(error),
        }
    }
}
type Error = MolecularNetworkError;

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkLinkStep {
    pub id: u64,
    pub gases: [u64; 2],
    pub surface: MolecularExchangeSurface,
    pub exchange: DriftingEffusiveExchange,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkStep {
    pub duration_s: f64,
    pub energy_frame_velocity_m_s: [f64; 3],
    /// Canonical link order, independent of input ordering.
    pub links: Vec<MolecularNetworkLinkStep>,
    /// One original accepted reconstruction per gas, including isolated gases.
    pub gas: BTreeMap<u64, MovingGasBoundaryStep>,
    pub maximum_stream_mass_rounding_kg: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub momentum_residual_n_s: [f64; 3],
    pub stream_energy_residual_j: f64,
    /// Actual rounded pooled rows in a common original-species enthalpy datum.
    pub pooled_stream_energy_residual_j: f64,
    pub energy_residual_j: f64,
    pub laboratory_energy_residual_j: f64,
    pub mass_allowance_kg: f64,
    pub momentum_allowance_n_s: [f64; 3],
    pub energy_allowance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkCandidate {
    gases: BTreeMap<u64, MovingGasParcel>,
    step: MolecularNetworkStep,
}
impl MolecularNetworkCandidate {
    pub fn gases(&self) -> &BTreeMap<u64, MovingGasParcel> {
        &self.gases
    }
    pub const fn step(&self) -> &MolecularNetworkStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, MovingGasParcel>, MolecularNetworkStep) {
        (self.gases, self.step)
    }
}

#[derive(Default)]
struct SpeciesPool {
    mass: Dd,
    energy: Dd,
    outgoing: Dd,
    gross: Dd,
    rounding: Dd,
}
#[derive(Default)]
struct Pool<'a> {
    catalog: Option<&'a ThermochemicalGasModel>,
    species: BTreeMap<SpeciesId, SpeciesPool>,
    impulse: [Dd; 3],
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn source(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Link { id, error }
    }
}
fn gas_error(id: u64, error: MovingGasBoundaryError) -> Error {
    match error {
        MovingGasBoundaryError::Cancelled
        | MovingGasBoundaryError::Gas(GasEnergyError::Cancelled) => Error::Cancelled,
        error => Error::Gas { id, error },
    }
}

/// Sample every aperture before reconstructing any finite gas.
///
/// Separate link IDs must identify disjoint exposed physical areas. Reverse a
/// link by swapping its gas IDs and negating its normal, keeping its velocity.
/// Closed links do not query source/drift data. Gas solve IDs must exactly match
/// the input gas IDs. An empty link set preserves all original gases exactly.
/// Common-frame conversions use net stream increments, so opposite gross
/// formation energies and common inertial drift are never subtracted.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    controls: MolecularNetworkControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MolecularNetworkCandidate, Error> {
    poll(&mut cancelled)?;
    if gases.is_empty()
        || gases.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || !gases.keys().eq(endpoints.keys())
        || !positive(duration_s)
        || !positive(controls.maximum_duration_s)
        || !positive(controls.maximum_species_fraction)
        || controls.maximum_species_fraction >= 1.0
        || !positive(controls.mass_absolute_kg)
        || !positive(controls.momentum_absolute_n_s)
        || !positive(controls.energy_absolute_j)
        || !controls.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_tolerance)
        || controls
            .energy_frame_velocity_m_s
            .iter()
            .any(|v| !v.is_finite())
    {
        return Err(Error::Invalid(
            "bounded gas/link counts, exact solve IDs and positive physical budgets required",
        ));
    }
    if duration_s > controls.maximum_duration_s {
        return Err(MolecularExchangeError::Duration {
            duration_s,
            maximum_duration_s: controls.maximum_duration_s,
        }
        .into());
    }
    for (&id, solve) in endpoints {
        poll(&mut cancelled)?;
        solve
            .gas
            .gas
            .validate()
            .map_err(|error| gas_error(id, MovingGasBoundaryError::Gas(error)))?;
        if !positive(solve.gas.mass_tolerance_kg)
            || !positive(solve.gas.energy_tolerance_j)
            || !positive(solve.momentum_absolute_n_s)
            || !solve.momentum_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&solve.momentum_relative_tolerance)
        {
            return Err(Error::Invalid("positive original per-gas budgets required"));
        }
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        poll(&mut cancelled)?;
        if link.gases[0] == link.gases[1]
            || link.gases.iter().any(|id| !gases.contains_key(id))
            || !link.surface.area_m2.is_finite()
            || link.surface.area_m2 < 0.0
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "unique aperture IDs, distinct existing endpoints and nonnegative area required",
            ));
        }
    }
    let frame = controls.energy_frame_velocity_m_s;
    let mut pools: BTreeMap<_, Pool<'_>> = gases.keys().map(|&id| (id, Pool::default())).collect();
    let mut samples = Vec::with_capacity(ordered.len());
    let mut stream_energy = Dd::ZERO;
    let mut energy_scale = Dd::ZERO;
    let mut gross_mass = Dd::ZERO;
    let mut maximum_rounding = 0.0_f64;
    let mut energy_references = BTreeMap::new();
    for (id, link) in ordered {
        poll(&mut cancelled)?;
        let active = link.surface.area_m2 > 0.0;
        let aperture = if active {
            link.surface.velocity_world_m_s
        } else {
            frame
        };
        if aperture.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("active aperture velocity must be finite"));
        }
        let original = link.gases.map(|gas| &gases[&gas]);
        let relative = if active {
            original.map(|gas| core::array::from_fn(|d| gas.velocity_m_s()[d] - aperture[d]))
        } else {
            [[0.0; 3]; 2]
        };
        let sampled = exchange::sample(
            original[0].parcel(),
            original[1].parcel(),
            link.interface,
            link.surface.area_m2,
            relative,
            link.surface.normal_left_to_right,
            duration_s,
            &mut cancelled,
        )
        .map_err(|error| source(id, error))?;
        let catalog = link
            .interface
            .momentum_interface()
            .interface()
            .gas_catalog();
        if !sampled.species.is_empty() {
            for gas in link.gases {
                let pool = pools.get_mut(&gas).expect("admitted gas ID");
                if pool.catalog.is_some_and(|old| old != catalog) {
                    return Err(Error::Invalid(
                        "incident active apertures must retain the same original gas catalog",
                    ));
                }
                pool.catalog = Some(catalog);
            }
        }
        for row in &sampled.species {
            poll(&mut cancelled)?;
            let caloric = catalog
                .thermodynamics()
                .components()
                .iter()
                .find(|species| species.species() == &row.species)
                .ok_or(Error::Invalid("original molecular source lost its species"))?;
            let reference = energy_references
                .entry(row.species.clone())
                .or_insert((caloric, original[0].parcel().state().state().temperature));
            if reference.0 != caloric {
                return Err(Error::Invalid(
                    "one network species cannot substitute its original caloric source",
                ));
            }
            let dm = Dd::from_f64(row.net_mass_left_to_right_kg);
            let gross = Dd::from_f64(row.mass_kg[0]) + Dd::from_f64(row.mass_kg[1]);
            let rounding = Dd::from_f64(row.mass_kg[0]) - Dd::from_f64(row.mass_kg[1]) - dm;
            let limit = endpoints[&link.gases[0]]
                .gas
                .mass_tolerance_kg
                .min(endpoints[&link.gases[1]].gas.mass_tolerance_kg);
            maximum_rounding = maximum_rounding
                .max(balance("individual network gross/net mass", rounding, limit)?.abs());
            gross_mass = gross_mass + gross;
            let impulse: [Dd; 3] = core::array::from_fn(|d| {
                Dd::from_f64(row.momentum_relative_kg_m_s[0][d])
                    - Dd::from_f64(row.momentum_relative_kg_m_s[1][d])
            });
            // E_common = E_aperture + (V_aperture-V_common).P_aperture
            //            + |V_aperture-V_common|² dm/2.
            let mut conversion = Dd::ZERO;
            let mut world_impulse = [Dd::ZERO; 3];
            for d in 0..3 {
                let dv = Dd::from_f64(aperture[d]) - Dd::from_f64(frame[d]);
                conversion = conversion + dv * impulse[d] + (dm * Dd::from_f64(0.5) * dv) * dv;
                world_impulse[d] = impulse[d] + Dd::from_f64(aperture[d]) * dm;
            }
            let dh = if row.net_mass_left_to_right_kg == 0.0 {
                0.0
            } else {
                catalog
                    .species_specific_enthalpy_change_j_kg(
                        &row.species,
                        original[0].parcel().state().state().temperature,
                        original[1].parcel().state().state().temperature,
                    )
                    .map_err(|error| source(id, error))?
            };
            stream_energy = stream_energy
                + Dd::from_f64(row.relative_energy_into_gas_j[0])
                + Dd::from_f64(row.relative_energy_into_gas_j[1])
                + dm * Dd::from_f64(dh);
            for side in 0..2 {
                let sign = Dd::from_f64(if side == 0 { -1.0 } else { 1.0 });
                let pool = pools.get_mut(&link.gases[side]).expect("admitted gas ID");
                let species = pool.species.entry(row.species.clone()).or_default();
                species.mass = species.mass + sign * dm;
                let energy = Dd::from_f64(row.relative_energy_into_gas_j[side]) + sign * conversion;
                species.energy = species.energy + energy;
                species.outgoing = species.outgoing + Dd::from_f64(row.mass_kg[side]);
                species.gross = species.gross + gross;
                // Absolute error cannot cancel across separately funded apertures.
                species.rounding = species.rounding + rounding.abs();
                energy_scale = energy_scale + energy.abs();
                for (sum, increment) in pool.impulse.iter_mut().zip(world_impulse) {
                    *sum = *sum + sign * increment;
                }
                if pool.species.len() > MAX_EFFUSIVE_SPECIES {
                    return Err(Error::Invalid(
                        "incident aperture species exceed the original work bound",
                    ));
                }
            }
        }
        samples.push(MolecularNetworkLinkStep {
            id,
            gases: link.gases,
            surface: link.surface,
            exchange: sampled,
        });
    }
    // Admit ALL gross withdrawals from ORIGINAL inventories before any endpoint.
    let mut inventories = BTreeMap::new();
    for (&id, pool) in &pools {
        poll(&mut cancelled)?;
        let inventory = if pool.species.is_empty() {
            BTreeMap::new()
        } else {
            gases[&id]
                .parcel()
                .species_inventory_kg()
                .map_err(|error| gas_error(id, MovingGasBoundaryError::Gas(error)))?
        };
        let mut total_rounding = Dd::ZERO;
        for (species, row) in &pool.species {
            poll(&mut cancelled)?;
            let available = inventory.get(species).copied().unwrap_or(0.0);
            if Dd::ZERO.lt(row.outgoing) && !positive(available) {
                return Err(Error::Invalid(
                    "outgoing molecular species lost its original inventory",
                ));
            }
            if (Dd::from_f64(available) * Dd::from_f64(controls.maximum_species_fraction))
                .lt(row.outgoing)
            {
                return Err(Error::Withdrawal {
                    gas: id,
                    species: species.clone(),
                    fraction: finite(row.outgoing)? / available,
                    maximum: controls.maximum_species_fraction,
                });
            }
            total_rounding = total_rounding + row.rounding;
        }
        maximum_rounding = maximum_rounding.max(
            balance(
                "pooled network gross/net mass uncertainty",
                total_rounding,
                endpoints[&id].gas.mass_tolerance_kg,
            )?
            .abs(),
        );
        inventories.insert(id, inventory);
    }
    let mut next = BTreeMap::new();
    let mut gas_steps = BTreeMap::new();
    let mut mass_residual = Dd::ZERO;
    let mut momentum_residual = [Dd::ZERO; 3];
    let mut momentum_scale = [Dd::ZERO; 3];
    let mut caloric_residual = Dd::ZERO;
    let mut pooled_stream_energy = Dd::ZERO;
    let mut species_accounts: BTreeMap<SpeciesId, (Dd, Dd)> = BTreeMap::new();
    for (&id, pool) in &pools {
        poll(&mut cancelled)?;
        let original = &gases[&id];
        let rows: Vec<_> = pool
            .species
            .iter()
            .map(|(species, row)| {
                Ok(SpeciesBoundaryRelativeTransfer {
                    species: species.clone(),
                    mass_into_gas_kg: finite(row.mass)?,
                    relative_energy_into_gas_j: finite(row.energy)?,
                })
            })
            .collect::<Result<_, MolecularExchangeError>>()?;
        let catalog = pool.catalog.unwrap_or(original.parcel().model());
        for row in &rows {
            poll(&mut cancelled)?;
            let reference_temperature = energy_references[&row.species].1;
            let dh = if row.mass_into_gas_kg == 0.0 {
                0.0
            } else {
                catalog
                    .species_specific_enthalpy_change_j_kg(
                        &row.species,
                        reference_temperature,
                        original.parcel().state().state().temperature,
                    )
                    .map_err(|error| gas_error(id, MovingGasBoundaryError::Gas(error)))?
            };
            // Audit what the caloric owner ACTUALLY receives, after pooling and
            // f64 conversion. A shared species datum cancels formation constants.
            pooled_stream_energy = pooled_stream_energy
                + Dd::from_f64(row.relative_energy_into_gas_j)
                + Dd::from_f64(row.mass_into_gas_kg) * Dd::from_f64(dh);
        }
        let impulse = [
            finite(pool.impulse[0])?,
            finite(pool.impulse[1])?,
            finite(pool.impulse[2])?,
        ];
        let accepted = original
            .exchange_species_momentum_relative(
                catalog,
                &rows,
                impulse,
                frame,
                endpoints[&id],
                &mut cancelled,
            )
            .map_err(|error| gas_error(id, error))?;
        let (gas, _, step) = accepted.into_parts();
        mass_residual = mass_residual + Dd::from_f64(gas.parcel().mass_kg())
            - Dd::from_f64(original.parcel().mass_kg());
        for d in 0..3 {
            momentum_residual[d] = momentum_residual[d] + Dd::from_f64(gas.momentum_kg_m_s()[d])
                - Dd::from_f64(original.momentum_kg_m_s()[d]);
            momentum_scale[d] = momentum_scale[d] + pool.impulse[d].abs() * Dd::from_f64(0.5);
        }
        caloric_residual = caloric_residual + Dd::from_f64(step.total_energy_residual_j);
        energy_scale = energy_scale + Dd::from_f64(step.kinetic_energy_change_j).abs();
        if !pool.species.is_empty() {
            let inventory = gas
                .parcel()
                .species_inventory_kg()
                .map_err(|error| gas_error(id, MovingGasBoundaryError::Gas(error)))?;
            let union: BTreeSet<_> = inventories[&id].keys().chain(inventory.keys()).collect();
            for species in union {
                poll(&mut cancelled)?;
                let account = species_accounts
                    .entry(species.clone())
                    .or_insert((Dd::ZERO, Dd::ZERO));
                account.0 = account.0
                    + Dd::from_f64(inventory.get(species).copied().unwrap_or(0.0))
                    - Dd::from_f64(inventories[&id].get(species).copied().unwrap_or(0.0));
                if let Some(row) = pool.species.get(species) {
                    account.1 = account.1 + row.gross * Dd::from_f64(0.5);
                }
            }
        }
        next.insert(id, gas);
        gas_steps.insert(id, step);
    }
    for row in &samples {
        poll(&mut cancelled)?;
        let relative = row.gases.map(|id| {
            core::array::from_fn(|d| {
                next[&id].velocity_m_s()[d] - row.surface.velocity_world_m_s[d]
            })
        });
        row.exchange
            .admit_endpoints(
                [next[&row.gases[0]].parcel(), next[&row.gases[1]].parcel()],
                relative,
                &mut cancelled,
            )
            .map_err(|error| source(row.id, error))?;
    }
    let mass_allowance_kg = allowance(
        controls.mass_absolute_kg,
        controls.relative_tolerance,
        gross_mass,
    )?;
    let mass_residual_kg = balance("actual network gas mass", mass_residual, mass_allowance_kg)?;
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for (residual, gross) in species_accounts.values() {
        poll(&mut cancelled)?;
        let limit = allowance(
            controls.mass_absolute_kg,
            controls.relative_tolerance,
            *gross,
        )?;
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg
            .max(balance("actual network species mass", *residual, limit)?.abs());
    }
    let mut momentum_residual_n_s = [0.0; 3];
    let mut momentum_allowance_n_s = [0.0; 3];
    let mut conversion = Dd::ZERO;
    for d in 0..3 {
        momentum_allowance_n_s[d] = allowance(
            controls.momentum_absolute_n_s,
            controls.relative_tolerance,
            momentum_scale[d],
        )?;
        momentum_residual_n_s[d] = balance(
            "actual network momentum",
            momentum_residual[d],
            momentum_allowance_n_s[d],
        )?;
        let velocity = Dd::from_f64(frame[d]);
        conversion = conversion + velocity * momentum_residual[d]
            - (mass_residual * Dd::from_f64(0.5) * velocity) * velocity;
    }
    let energy_allowance_j = allowance(
        controls.energy_absolute_j,
        controls.relative_tolerance,
        energy_scale,
    )?;
    let stream_energy_residual_j = balance(
        "network reciprocal stream energy",
        stream_energy,
        energy_allowance_j,
    )?;
    let pooled_stream_energy_residual_j = balance(
        "rounded pooled molecular energy",
        pooled_stream_energy,
        energy_allowance_j,
    )?;
    let energy_residual = pooled_stream_energy + caloric_residual;
    let energy_residual_j = balance(
        "network common-frame energy",
        energy_residual,
        energy_allowance_j,
    )?;
    let laboratory_energy_residual_j = balance(
        "network inertial-frame energy",
        energy_residual + conversion,
        energy_allowance_j,
    )?;
    poll(&mut cancelled)?;
    Ok(MolecularNetworkCandidate {
        gases: next,
        step: MolecularNetworkStep {
            duration_s,
            energy_frame_velocity_m_s: frame,
            links: samples,
            gas: gas_steps,
            maximum_stream_mass_rounding_kg: maximum_rounding,
            mass_residual_kg,
            maximum_species_mass_residual_kg,
            momentum_residual_n_s,
            stream_energy_residual_j,
            pooled_stream_energy_residual_j,
            energy_residual_j,
            laboratory_energy_residual_j,
            mass_allowance_kg,
            momentum_allowance_n_s,
            energy_allowance_j,
        },
    })
}
