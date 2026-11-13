//! Fixed original caloric references and actual accepted molecular accounts.

use super::super::MolecularNetworkLinkStep;
use super::{
    Error, MAX_EFFUSIVE_SPECIES, MolecularNetworkAdaptiveControls, MolecularNetworkLink,
    MolecularNetworkLinkTotals, MolecularNetworkStep, MolecularNetworkStreamMoments,
    MolecularNetworkWindowBalance, MovingGasParcel, SpeciesId, poll, ratio, unresolved,
};
use crate::liquid_gas::evaporation::momentum::boundary::kinetic_change;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::{BTreeMap, BTreeSet};

struct Reference<'a> {
    source: &'a Nasa9StandardStateModelV1,
    temperature_k: f64,
}

pub(in crate::molecular_exchange::network) struct Sources<'a> {
    pub axes: BTreeMap<u64, &'a [Nasa9StandardStateModelV1]>,
    references: BTreeMap<SpeciesId, Reference<'a>>,
}

fn inventory(
    gas: u64,
    state: &MovingGasParcel,
    at: f64,
) -> Result<BTreeMap<SpeciesId, f64>, Error> {
    state
        .parcel()
        .species_inventory_kg()
        .map_err(|error| Error::Inventory {
            gas,
            elapsed_s: at,
            error,
        })
}

impl<'a> Sources<'a> {
    pub fn new(
        gases: &'a BTreeMap<u64, MovingGasParcel>,
        links: &'a [MolecularNetworkLink<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut ordered = BTreeMap::new();
        for link in links {
            poll(cancelled, 0.0)?;
            if link.gases[0] == link.gases[1]
                || link.gases.iter().any(|id| !gases.contains_key(id))
                || !link.surface.area_m2.is_finite()
                || link.surface.area_m2 < 0.0
                || ordered.insert(link.id, link).is_some()
            {
                return Err(Error::Invalid(
                    "unique molecular links, distinct existing gases and nonnegative areas required",
                ));
            }
        }
        let mut catalogs: BTreeMap<u64, &ThermochemicalGasModel> = BTreeMap::new();
        for link in ordered.values().filter(|link| link.surface.area_m2 > 0.0) {
            poll(cancelled, 0.0)?;
            let catalog = link
                .interface
                .momentum_interface()
                .interface()
                .gas_catalog();
            if catalog.thermodynamics().components().len() > MAX_EFFUSIVE_SPECIES {
                return Err(Error::Invalid(
                    "original molecular species work bound exceeded",
                ));
            }
            for id in link.gases {
                if catalogs
                    .insert(id, catalog)
                    .is_some_and(|old| old != catalog)
                {
                    return Err(Error::Invalid(
                        "incident active apertures must retain the same original gas catalog",
                    ));
                }
            }
        }
        let axes: BTreeMap<_, _> = catalogs
            .iter()
            .map(|(&id, catalog)| (id, catalog.thermodynamics().components()))
            .collect();
        let mut references: BTreeMap<SpeciesId, Reference<'a>> = BTreeMap::new();
        // A real positive initial inventory, never an absent receiver's T,
        // establishes each fixed reference. Newly arriving species need no seed.
        for (&id, axis) in &axes {
            poll(cancelled, 0.0)?;
            let original = &gases[&id];
            let masses = inventory(id, original, 0.0)?;
            for source in original.parcel().model().thermodynamics().components() {
                poll(cancelled, 0.0)?;
                if masses.get(source.species()).copied().unwrap_or(0.0) == 0.0 {
                    continue;
                }
                let index = axis
                    .binary_search_by(|s| s.species().cmp(source.species()))
                    .map_err(|_| {
                        Error::Invalid("original gas species absent from incident catalog")
                    })?;
                if &axis[index] != source {
                    return Err(Error::Invalid(
                        "original gas caloric source differs from its incident catalog",
                    ));
                }
                if let Some(reference) = references.get(source.species()) {
                    if reference.source != source {
                        return Err(Error::Invalid(
                            "active network species must retain one original caloric source",
                        ));
                    }
                } else {
                    references.insert(
                        source.species().clone(),
                        Reference {
                            source,
                            temperature_k: original.parcel().state().state().temperature,
                        },
                    );
                }
            }
        }
        Ok(Self { axes, references })
    }

    /// Register only gases explicitly heated by a composing thermal owner.
    /// The original molecular-only path still ignores isolated source catalogs.
    pub fn with_thermal_gases(
        gases: &'a BTreeMap<u64, MovingGasParcel>,
        links: &'a [MolecularNetworkLink<'_>],
        thermal_gases: impl IntoIterator<Item = u64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut result = Self::new(gases, links, cancelled)?;
        let mut ordered = BTreeSet::new();
        for id in thermal_gases {
            poll(cancelled, 0.0)?;
            ordered.insert(id);
        }
        for id in ordered {
            poll(cancelled, 0.0)?;
            let gas = gases.get(&id).ok_or(Error::Invalid(
                "thermal source references require an existing original gas",
            ))?;
            if result.axes.contains_key(&id) {
                continue;
            }
            let axis = gas.parcel().model().thermodynamics().components();
            if axis.len() > MAX_EFFUSIVE_SPECIES {
                return Err(Error::Invalid(
                    "original thermal species work bound exceeded",
                ));
            }
            let masses = inventory(id, gas, 0.0)?;
            for source in axis {
                poll(cancelled, 0.0)?;
                if masses.get(source.species()).copied().unwrap_or(0.0) == 0.0 {
                    continue;
                }
                if let Some(reference) = result.references.get(source.species()) {
                    if reference.source != source {
                        return Err(Error::Invalid(
                            "thermal network species must retain one original caloric source",
                        ));
                    }
                } else {
                    result.references.insert(
                        source.species().clone(),
                        Reference {
                            source,
                            temperature_k: gas.parcel().state().state().temperature,
                        },
                    );
                }
            }
            result.axes.insert(id, axis);
        }
        Ok(result)
    }

    pub fn masses(
        &self,
        id: u64,
        gas: &MovingGasParcel,
        at: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; MAX_EFFUSIVE_SPECIES], Error> {
        let axis = self.axes[&id];
        let mut result = [0.0; MAX_EFFUSIVE_SPECIES];
        for (species, mass) in inventory(id, gas, at)? {
            poll(cancelled, at)?;
            let index = axis
                .binary_search_by(|source| source.species().cmp(&species))
                .map_err(|_| {
                    Error::Invalid("accepted species absent from original incident catalog")
                })?;
            result[index] = mass;
        }
        Ok(result)
    }
}

fn value(x: Dd, at: f64) -> Result<f64, Error> {
    let value = x.to_f64();
    if x.hi.is_finite() && x.lo.is_finite() && value.is_finite() {
        Ok(value)
    } else {
        Err(unresolved(at, "molecular trajectory account"))
    }
}
fn maximum(a: Dd, b: Dd) -> Dd {
    if a.lt(b) { b } else { a }
}
fn du(source: &Nasa9StandardStateModelV1, from: f64, to: f64, at: f64) -> Result<Dd, Error> {
    if from == to {
        return Ok(Dd::ZERO);
    }
    let change = source
        .molar_internal_energy_change(Temperature::new(from), Temperature::new(to))
        .map_err(|error| Error::Caloric {
            species: source.species().clone(),
            elapsed_s: at,
            error,
        })?;
    Ok(Dd::from_f64(change.value()) / Dd::from_f64(source.molar_mass().value()))
}
fn reference_enthalpy(reference: &Reference<'_>, t: f64, at: f64) -> Result<Dd, Error> {
    Ok(du(reference.source, reference.temperature_k, t, at)?
        + Dd::from_f64(R) * Dd::from_f64(t) / Dd::from_f64(reference.source.molar_mass().value()))
}

#[derive(Clone, Copy)]
struct Moments {
    mass: Dd,
    momentum: [Dd; 3],
    energy: Dd,
}
impl Default for Moments {
    fn default() -> Self {
        Self {
            mass: Dd::ZERO,
            momentum: [Dd::ZERO; 3],
            energy: Dd::ZERO,
        }
    }
}
impl Moments {
    fn add(&mut self, other: Self) {
        self.mass = self.mass + other.mass;
        self.energy = self.energy + other.energy;
        for d in 0..3 {
            self.momentum[d] = self.momentum[d] + other.momentum[d];
        }
    }
    fn report(self, at: f64) -> Result<MolecularNetworkStreamMoments, Error> {
        Ok(MolecularNetworkStreamMoments {
            mass_kg: value(self.mass, at)?,
            momentum_n_s: [
                value(self.momentum[0], at)?,
                value(self.momentum[1], at)?,
                value(self.momentum[2], at)?,
            ],
            reference_energy_j: value(self.energy, at)?,
        })
    }
    // Original aperture P and E become common-frame P and E for THESE particles.
    fn translate(&mut self, delta_velocity: [Dd; 3]) {
        for (d, v) in delta_velocity.into_iter().enumerate() {
            self.energy =
                self.energy + v * self.momentum[d] + (self.mass * Dd::from_f64(0.5) * v) * v;
            self.momentum[d] = self.momentum[d] + v * self.mass;
        }
    }
}

fn streams(
    before: &BTreeMap<u64, MovingGasParcel>,
    step: &MolecularNetworkLinkStep,
    sources: &Sources<'_>,
    frame: [f64; 3],
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BTreeMap<SpeciesId, [Moments; 3]>, Error> {
    let mut result = BTreeMap::new();
    // Closed apertures do not query normals, velocity or inactive sources.
    if step.exchange.species.is_empty() {
        return Ok(result);
    }
    let velocities = core::array::from_fn(|d| {
        Dd::from_f64(step.surface.velocity_world_m_s[d]) - Dd::from_f64(frame[d])
    });
    let temperature = step
        .gases
        .map(|id| before[&id].parcel().state().state().temperature);
    let beams = [&step.exchange.left_to_right, &step.exchange.right_to_left];
    for row in &step.exchange.species {
        poll(cancelled, at)?;
        let reference = sources.references.get(&row.species).ok_or(Error::Invalid(
            "exchanged species has no original positive source inventory",
        ))?;
        let mut moments = [Moments::default(); 3];
        let mut enthalpy = [Dd::ZERO; 2];
        for side in 0..2 {
            if row.mass_kg[side] != 0.0 || (side == 0 && row.net_mass_left_to_right_kg != 0.0) {
                enthalpy[side] = reference_enthalpy(reference, temperature[side], at)?;
            }
            if row.mass_kg[side] == 0.0 {
                continue;
            }
            let beam = beams[side]
                .species
                .binary_search_by(|b| b.species.cmp(&row.species))
                .ok()
                .map(|i| &beams[side].species[i])
                .ok_or(Error::Invalid(
                    "accepted molecular mass lost its original beam",
                ))?;
            let correction =
                beam.incoming
                    .additional_specific_energy_j_kg
                    .ok_or(Error::Invalid(
                        "accepted molecular beam lost its source-relative energy",
                    ))?;
            moments[side] = Moments {
                mass: Dd::from_f64(row.mass_kg[side]),
                momentum: row.momentum_relative_kg_m_s[side].map(Dd::from_f64),
                energy: Dd::from_f64(row.mass_kg[side])
                    * (enthalpy[side] + Dd::from_f64(correction)),
            };
            moments[side].translate(velocities);
        }
        // Preserve the sampler's small-drift net coordinate. Subtracting two
        // rounded gross masses or absolute formation energies would lose it.
        moments[2] = Moments {
            mass: Dd::from_f64(row.net_mass_left_to_right_kg),
            momentum: core::array::from_fn(|d| {
                Dd::from_f64(row.momentum_relative_kg_m_s[0][d])
                    - Dd::from_f64(row.momentum_relative_kg_m_s[1][d])
            }),
            energy: Dd::from_f64(row.net_mass_left_to_right_kg) * enthalpy[0]
                - Dd::from_f64(row.relative_energy_into_gas_j[0]),
        };
        moments[2].translate(velocities);
        result.insert(row.species.clone(), moments);
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn stream_error(
    before: &BTreeMap<u64, MovingGasParcel>,
    middle: &BTreeMap<u64, MovingGasParcel>,
    coarse: &MolecularNetworkStep,
    first: &MolecularNetworkStep,
    second: &MolecularNetworkStep,
    sources: &Sources<'_>,
    controls: MolecularNetworkAdaptiveControls,
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    stream_error_from_states(
        [before, before, middle],
        [coarse, first, second],
        sources,
        controls,
        at,
        cancelled,
    )
}

/// Distinct physical source snapshots for a composed coarse/half/half trial.
pub(in crate::molecular_exchange::network) fn stream_error_from_states(
    states: [&BTreeMap<u64, MovingGasParcel>; 3],
    steps: [&MolecularNetworkStep; 3],
    sources: &Sources<'_>,
    controls: MolecularNetworkAdaptiveControls,
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let [before, first_before, middle] = states;
    let [coarse, first, second] = steps;
    if coarse.links.len() != first.links.len() || coarse.links.len() != second.links.len() {
        return Err(Error::Invalid(
            "native molecular link set changed during a trial",
        ));
    }
    let mut error = 0.0_f64;
    for ((a, b), c) in coarse.links.iter().zip(&first.links).zip(&second.links) {
        poll(cancelled, at)?;
        if a.id != b.id || a.id != c.id || a.gases != b.gases || a.gases != c.gases {
            return Err(Error::Invalid(
                "native molecular link identities changed during a trial",
            ));
        }
        let values = [
            streams(
                before,
                a,
                sources,
                coarse.energy_frame_velocity_m_s,
                at,
                cancelled,
            )?,
            streams(
                first_before,
                b,
                sources,
                first.energy_frame_velocity_m_s,
                at,
                cancelled,
            )?,
            streams(
                middle,
                c,
                sources,
                second.energy_frame_velocity_m_s,
                at,
                cancelled,
            )?,
        ];
        let species: BTreeSet<_> = values.iter().flat_map(|v| v.keys()).collect();
        for species in species {
            poll(cancelled, at)?;
            let [full, one, two] = values
                .each_ref()
                .map(|v| v.get(species).copied().unwrap_or([Moments::default(); 3]));
            for side in 0..3 {
                let mut fine = one[side];
                fine.add(two[side]);
                let (a, b) = (full[side].report(at)?, fine.report(at)?);
                for (a, b, absolute) in [
                    (
                        a.mass_kg,
                        b.mass_kg,
                        controls.states.species_mass_tolerance_kg,
                    ),
                    (
                        a.reference_energy_j,
                        b.reference_energy_j,
                        controls.stream_energy_tolerance_j,
                    ),
                    (
                        a.momentum_n_s[0],
                        b.momentum_n_s[0],
                        controls.stream_momentum_tolerance_n_s,
                    ),
                    (
                        a.momentum_n_s[1],
                        b.momentum_n_s[1],
                        controls.stream_momentum_tolerance_n_s,
                    ),
                    (
                        a.momentum_n_s[2],
                        b.momentum_n_s[2],
                        controls.stream_momentum_tolerance_n_s,
                    ),
                ] {
                    error = error.max(ratio(
                        a,
                        b,
                        a.abs().max(b.abs()),
                        absolute,
                        controls.states.relative_tolerance,
                        at,
                    )?);
                }
            }
        }
    }
    Ok(error)
}

/// Original-to-final physical U + K, with one fixed original u_i datum removed
/// per conserved species. For species absent initially, only the actual final
/// temperature is queried. K uses the original changing-mass momentum identity.
pub(in crate::molecular_exchange::network) fn window_balance(
    before: &BTreeMap<u64, MovingGasParcel>,
    after: &BTreeMap<u64, MovingGasParcel>,
    sources: &Sources<'_>,
    frame: [f64; 3],
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MolecularNetworkWindowBalance, Error> {
    let (mut mass, mut energy, mut laboratory) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    let mut momentum = [Dd::ZERO; 3];
    let mut species: BTreeMap<SpeciesId, Dd> = BTreeMap::new();
    for (&id, old) in before {
        poll(cancelled, at)?;
        let new = &after[&id];
        if new == old {
            continue;
        }
        let axis = sources.axes.get(&id).ok_or(Error::Invalid(
            "isolated molecular gas changed during exchange",
        ))?;
        mass = mass + Dd::from_f64(new.parcel().mass_kg()) - Dd::from_f64(old.parcel().mass_kg());
        for d in 0..3 {
            momentum[d] = momentum[d] + Dd::from_f64(new.momentum_kg_m_s()[d])
                - Dd::from_f64(old.momentum_kg_m_s()[d]);
        }
        let m0 = inventory(id, old, at)?;
        let m1 = inventory(id, new, at)?;
        let t0 = old.parcel().state().state().temperature;
        let t1 = new.parcel().state().state().temperature;
        let mut internal = Dd::ZERO;
        for source in *axis {
            poll(cancelled, at)?;
            let a = m0.get(source.species()).copied().unwrap_or(0.0);
            let b = m1.get(source.species()).copied().unwrap_or(0.0);
            if a == 0.0 && b == 0.0 {
                continue;
            }
            let reference = sources
                .references
                .get(source.species())
                .ok_or(Error::Invalid(
                    "accepted species lost its original caloric reference",
                ))?;
            if reference.source != source {
                return Err(Error::Invalid(
                    "accepted species changed its original caloric source",
                ));
            }
            let dm = Dd::from_f64(b) - Dd::from_f64(a);
            let total = species.entry(source.species().clone()).or_insert(Dd::ZERO);
            *total = *total + dm;
            internal = internal
                + if a == 0.0 {
                    Dd::from_f64(b) * du(source, reference.temperature_k, t1, at)?
                } else {
                    let baseline = du(source, reference.temperature_k, t0, at)?;
                    if b == 0.0 {
                        dm * baseline
                    } else {
                        Dd::from_f64(b) * du(source, t0, t1, at)? + dm * baseline
                    }
                };
        }
        let kinetic = |frame| {
            kinetic_change(
                old.parcel().mass_kg(),
                new.parcel().mass_kg(),
                old.momentum_kg_m_s(),
                new.momentum_kg_m_s(),
                frame,
            )
            .map(Dd::from_f64)
            .map_err(|error| Error::Inventory {
                gas: id,
                elapsed_s: at,
                error,
            })
        };
        energy = energy + internal + kinetic(frame)?;
        laboratory = laboratory + internal + kinetic([0.0; 3])?;
    }
    let maximum_species = species
        .values()
        .fold(Dd::ZERO, |max, x| maximum(max, x.abs()));
    Ok(MolecularNetworkWindowBalance {
        mass_kg: value(mass, at)?,
        maximum_species_mass_kg: value(maximum_species, at)?,
        momentum_n_s: [
            value(momentum[0], at)?,
            value(momentum[1], at)?,
            value(momentum[2], at)?,
        ],
        energy_j: value(energy, at)?,
        laboratory_energy_j: value(laboratory, at)?,
    })
}

struct LinkTotal {
    gases: [u64; 2],
    moments: [Moments; 3],
}
pub(in crate::molecular_exchange::network) struct Totals {
    links: BTreeMap<u64, LinkTotal>,
    // mass, maximum-species bound, three momenta, common energy, lab energy.
    debit: [Dd; 7],
}
impl Default for Totals {
    fn default() -> Self {
        Self {
            links: BTreeMap::new(),
            debit: [Dd::ZERO; 7],
        }
    }
}
impl Totals {
    pub fn accept(
        &mut self,
        before: &BTreeMap<u64, MovingGasParcel>,
        after: &BTreeMap<u64, MovingGasParcel>,
        step: &MolecularNetworkStep,
        sources: &Sources<'_>,
        at: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let mut rounding = Dd::ZERO;
        for link in &step.links {
            poll(cancelled, at)?;
            let total = self.links.entry(link.id).or_insert(LinkTotal {
                gases: link.gases,
                moments: [Moments::default(); 3],
            });
            if total.gases != link.gases {
                return Err(Error::Invalid("accepted native molecular topology changed"));
            }
            for row in streams(
                before,
                link,
                sources,
                step.energy_frame_velocity_m_s,
                at,
                cancelled,
            )?
            .values()
            {
                for side in 0..3 {
                    total.moments[side].add(row[side]);
                }
            }
            for row in &link.exchange.species {
                poll(cancelled, at)?;
                rounding = rounding
                    + (Dd::from_f64(row.mass_kg[0])
                        - Dd::from_f64(row.mass_kg[1])
                        - Dd::from_f64(row.net_mass_left_to_right_kg))
                    .abs();
            }
        }
        let actual = window_balance(
            before,
            after,
            sources,
            step.energy_frame_velocity_m_s,
            at,
            cancelled,
        )?;
        let (mut mass, mut species, mut energy) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        let (mut momentum, mut requested) = ([Dd::ZERO; 3], [Dd::ZERO; 3]);
        for gas in step.gas.values() {
            poll(cancelled, at)?;
            mass = mass + Dd::from_f64(gas.caloric.mass_residual_kg).abs();
            species = species + Dd::from_f64(gas.caloric.maximum_species_mass_residual_kg).abs();
            energy = energy + Dd::from_f64(gas.total_energy_residual_j).abs();
            for d in 0..3 {
                momentum[d] = momentum[d] + Dd::from_f64(gas.momentum_residual_n_s[d]).abs();
                requested[d] = requested[d] + Dd::from_f64(gas.requested_impulse_n_s[d]);
            }
        }
        // Each rounded gross/net uncertainty belongs to BOTH finite endpoints.
        mass = mass + rounding * Dd::from_f64(2.0);
        species = species + rounding * Dd::from_f64(2.0);
        energy = energy + Dd::from_f64(step.pooled_stream_energy_residual_j).abs();
        self.debit[0] = self.debit[0] + maximum(mass, Dd::from_f64(actual.mass_kg).abs());
        self.debit[1] =
            self.debit[1] + maximum(species, Dd::from_f64(actual.maximum_species_mass_kg));
        let mut conversion = Dd::ZERO;
        for d in 0..3 {
            self.debit[2 + d] = self.debit[2 + d]
                + maximum(
                    momentum[d] + requested[d].abs(),
                    Dd::from_f64(actual.momentum_n_s[d]).abs(),
                );
            let v = Dd::from_f64(step.energy_frame_velocity_m_s[d]);
            conversion = conversion + v * Dd::from_f64(actual.momentum_n_s[d])
                - (Dd::from_f64(actual.mass_kg) * Dd::from_f64(0.5) * v) * v;
        }
        self.debit[5] = self.debit[5] + maximum(energy, Dd::from_f64(actual.energy_j).abs());
        self.debit[6] = self.debit[6]
            + maximum(
                energy + conversion.abs(),
                Dd::from_f64(actual.laboratory_energy_j).abs(),
            );
        Ok(())
    }
    pub fn debit(&self, at: f64) -> Result<MolecularNetworkWindowBalance, Error> {
        Ok(MolecularNetworkWindowBalance {
            mass_kg: value(self.debit[0], at)?,
            maximum_species_mass_kg: value(self.debit[1], at)?,
            momentum_n_s: [
                value(self.debit[2], at)?,
                value(self.debit[3], at)?,
                value(self.debit[4], at)?,
            ],
            energy_j: value(self.debit[5], at)?,
            laboratory_energy_j: value(self.debit[6], at)?,
        })
    }
    pub fn links(
        self,
        at: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<MolecularNetworkLinkTotals>, Error> {
        self.links
            .into_iter()
            .map(|(id, total)| {
                poll(cancelled, at)?;
                Ok(MolecularNetworkLinkTotals {
                    id,
                    gases: total.gases,
                    gross: [total.moments[0].report(at)?, total.moments[1].report(at)?],
                    net_left_to_right: total.moments[2].report(at)?,
                })
            })
            .collect()
    }
}
