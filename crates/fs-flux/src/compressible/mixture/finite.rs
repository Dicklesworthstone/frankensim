//! Shared finite, well-mixed mixture storage for Cartesian pressure ports.
//!
//! Dimensions only supply the native side-ordered face impulses. The original
//! mixture plenum owns the extensive species/energy update and caloric recovery.

use fs_gas::multicomponent::network::Plenum;
use fs_gas::multicomponent::{Cell, Conserved, Error, MAX_SPECIES, MixtureGas, Primitive};

/// Extensive, stationary tank inventory in the original canonical species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirInventory {
    pub species_masses: [f64; MAX_SPECIES],
    pub species_count: usize,
    /// Rigid physical volume [m³].
    pub volume_m3: f64,
    /// Species-referenced internal energy [J].
    pub energy: f64,
    pub pressure_pa: f64,
    pub temperature_k: f64,
}
impl FiniteReservoirInventory {
    pub fn mass(&self) -> f64 {
        sum(self.species_masses[..self.species_count].iter().copied())
    }

    pub fn absolute_energy(&self, gas: &MixtureGas) -> Result<f64, Error> {
        absolute_energy(gas, self.species_count, &self.species_masses, self.energy)
    }

    fn from_tank(tank: &Plenum, primitive: Primitive) -> Self {
        let q = tank.state();
        let mut species_masses = [0.0; MAX_SPECIES];
        for (mass, density) in species_masses.iter_mut().zip(q.partial_densities) {
            *mass = density * tank.volume_m3();
        }
        Self {
            species_masses,
            species_count: q.partial_densities.len(),
            volume_m3: tank.volume_m3(),
            energy: q.energy * tank.volume_m3(),
            pressure_pa: primitive.pressure,
            temperature_k: primitive.temperature,
        }
    }
}

/// Outward field transfer received by one tank and its stationary support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirTransfer<const D: usize> {
    /// Species masses entering the tank [kg]; negative values drain it.
    pub species_masses: [f64; MAX_SPECIES],
    pub species_count: usize,
    /// Complete Cartesian impulse delivered to the fixed support [N s].
    ///
    /// The tank has no resolved bulk momentum. The fixed support performs no
    /// work; incoming kinetic energy remains in the tank's internal energy.
    pub momentum: [f64; D],
    /// Complete transported total energy, including flow pressure work [J].
    pub energy: f64,
}
impl<const D: usize> FiniteReservoirTransfer<D> {
    pub fn mass(&self) -> f64 {
        sum(self.species_masses[..self.species_count].iter().copied())
    }

    pub fn absolute_energy(&self, gas: &MixtureGas) -> Result<f64, Error> {
        absolute_energy(gas, self.species_count, &self.species_masses, self.energy)
    }
}

/// Last accepted tank exchange over the same interval as the native field step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteReservoirReport<const D: usize> {
    pub before: FiniteReservoirInventory,
    pub after: FiniteReservoirInventory,
    pub transfer: FiniteReservoirTransfer<D>,
    /// Per-species inventory change less the native signed face transfer [kg].
    pub species_residual_kg: [f64; MAX_SPECIES],
    pub mass_residual_kg: f64,
    /// Referenced internal-energy change less native total-energy transfer [J].
    pub energy_residual_j: f64,
    /// The same residual restored to the original chemical-energy gauge [J].
    pub absolute_energy_residual_j: f64,
}

/// Prepared accepted and private stage/candidate storage. Side slots are
/// 2*axis + usize::from(positive_side); native rows are species, D momenta, E.
pub(in crate::compressible) struct FiniteReservoirs<const D: usize, const S: usize> {
    accepted: [Option<Plenum>; S],
    trial: [Option<Plenum>; S],
    primitive: [Option<Primitive>; S],
    trial_primitive: [Option<Primitive>; S],
    reports: [Option<FiniteReservoirReport<D>>; S],
    trial_reports: [Option<FiniteReservoirReport<D>>; S],
    allocated_bytes: usize,
}
impl<const D: usize, const S: usize> FiniteReservoirs<D, S> {
    pub(in crate::compressible) fn empty() -> Self {
        Self {
            accepted: std::array::from_fn(|_| None),
            trial: std::array::from_fn(|_| None),
            primitive: [None; S],
            trial_primitive: [None; S],
            reports: [None; S],
            trial_reports: [None; S],
            allocated_bytes: 0,
        }
    }

    pub(in crate::compressible) fn required_bytes(
        gas: &MixtureGas,
        count: usize,
    ) -> Result<usize, Error> {
        if !(D == 2 || D == 3) || S != 2 * D || count > S {
            return Err(Error::Invalid(
                "Cartesian finite mixture reservoir dimensions",
            ));
        }
        count
            .checked_mul(2 * gas.species_count() * size_of::<f64>())
            .and_then(|bytes| bytes.checked_add(size_of::<Self>()))
            .ok_or(Error::Invalid(
                "Cartesian finite mixture workspace overflow",
            ))
    }

    pub(in crate::compressible) fn workspace_bytes(&self) -> usize {
        size_of::<Self>() + self.allocated_bytes
    }

    pub(in crate::compressible) fn has_any(&self) -> bool {
        self.accepted.iter().any(Option::is_some)
    }

    pub(in crate::compressible) fn has_side(&self, slot: usize) -> bool {
        self.accepted[slot].is_some()
    }

    pub(in crate::compressible) fn bind(
        &mut self,
        gas: &MixtureGas,
        slot: usize,
        volume_m3: f64,
        initial: Cell<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if slot >= S
            || self.accepted[slot].is_some()
            || initial.partial_densities.len() != gas.species_count()
            || initial.momentum != 0.0
        {
            return Err(Error::Invalid(
                "unique stationary finite mixture reservoir required",
            ));
        }
        let primitive = gas.primitive(initial, &mut *cancelled)?;
        let make = || -> Result<(Plenum, usize), Error> {
            let mut partial_densities = Vec::new();
            partial_densities
                .try_reserve_exact(initial.partial_densities.len())
                .map_err(|_| Error::Allocation)?;
            partial_densities.extend_from_slice(initial.partial_densities);
            let bytes = partial_densities.capacity() * size_of::<f64>();
            Ok((
                Plenum::new(
                    volume_m3,
                    Conserved {
                        partial_densities,
                        momentum: 0.0,
                        energy: initial.energy,
                    },
                )?,
                bytes,
            ))
        };
        let (accepted, accepted_bytes) = make()?;
        let (trial, trial_bytes) = make()?;
        let bytes = self
            .allocated_bytes
            .checked_add(accepted_bytes)
            .and_then(|bytes| bytes.checked_add(trial_bytes))
            .ok_or(Error::Invalid(
                "Cartesian finite mixture workspace overflow",
            ))?;
        poll(cancelled)?;
        self.accepted[slot] = Some(accepted);
        self.trial[slot] = Some(trial);
        self.primitive[slot] = Some(primitive);
        self.trial_primitive[slot] = Some(primitive);
        self.allocated_bytes = bytes;
        Ok(())
    }

    pub(in crate::compressible) fn inventory(&self, slot: usize) -> Option<&Plenum> {
        self.accepted[slot].as_ref()
    }

    pub(in crate::compressible) fn state(&self, slot: usize, trial: bool) -> Option<Cell<'_>> {
        if trial {
            &self.trial[slot]
        } else {
            &self.accepted[slot]
        }
        .as_ref()
        .map(Plenum::state)
    }

    pub(in crate::compressible) fn primitive(&self, slot: usize, trial: bool) -> Option<Primitive> {
        if trial {
            self.trial_primitive[slot]
        } else {
            self.primitive[slot]
        }
    }

    pub(in crate::compressible) fn volume(&self, slot: usize) -> Option<f64> {
        self.accepted[slot].as_ref().map(Plenum::volume_m3)
    }

    pub(in crate::compressible) fn report(&self, slot: usize) -> Option<&FiniteReservoirReport<D>> {
        self.reports[slot].as_ref()
    }

    /// First capture holds half the forward-Euler impulse: weight two makes the
    /// full first tank stage. Weight one uses both captured SSP half impulses
    /// to make the final candidate, always relative to accepted inventory.
    pub(in crate::compressible) fn stage(
        &mut self,
        gas: &MixtureGas,
        boundary: &[f64],
        face_counts: [usize; S],
        weight: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if !self.has_any() {
            return Ok(());
        }
        poll(cancelled)?;
        let n = gas.species_count();
        let width = n + D + 1;
        let faces = face_counts
            .into_iter()
            .try_fold(0usize, |sum, count| sum.checked_add(count));
        if faces.and_then(|count| count.checked_mul(width)) != Some(boundary.len())
            || !(weight == 1.0 || weight == 2.0)
        {
            return Err(Error::Invalid(
                "Cartesian finite mixture boundary transfer shape",
            ));
        }
        let mut offset = 0;
        for (slot, count) in face_counts.into_iter().enumerate() {
            let end = offset + count * width;
            let rows = &boundary[offset..end];
            offset = end;
            let Some(before) = self.accepted[slot].as_ref() else {
                continue;
            };
            poll(cancelled)?;
            let mut totals = [0.0; MAX_SPECIES + 4];
            let mut corrections = [0.0; MAX_SPECIES + 4];
            for row in rows.chunks_exact(width) {
                poll(cancelled)?;
                for d in 0..width {
                    let y = row[d] - corrections[d];
                    let next = totals[d] + y;
                    corrections[d] = (next - totals[d]) - y;
                    totals[d] = next;
                }
            }
            for value in &mut totals[..width] {
                *value *= weight;
                if !value.is_finite() {
                    return Err(Error::Invalid(
                        "unrepresentable finite mixture port transfer",
                    ));
                }
            }
            let after = self.trial[slot]
                .as_mut()
                .expect("prepared finite mixture output");
            before.exchange_candidate_into(
                gas,
                &totals[..n],
                totals[width - 1],
                after,
                &mut *cancelled,
            )?;
            let primitive = gas.primitive(after.state(), &mut *cancelled)?;
            let before_inventory = FiniteReservoirInventory::from_tank(
                before,
                self.primitive[slot].expect("accepted finite mixture primitive"),
            );
            let after_inventory = FiniteReservoirInventory::from_tank(after, primitive);
            let mut species_masses = [0.0; MAX_SPECIES];
            species_masses[..n].copy_from_slice(&totals[..n]);
            let transfer = FiniteReservoirTransfer {
                species_masses,
                species_count: n,
                momentum: std::array::from_fn(|d| totals[n + d]),
                energy: totals[width - 1],
            };
            let mut species_residual_kg = [0.0; MAX_SPECIES];
            for d in 0..n {
                species_residual_kg[d] = after_inventory.species_masses[d]
                    - before_inventory.species_masses[d]
                    - transfer.species_masses[d];
            }
            let energy_residual_j =
                after_inventory.energy - before_inventory.energy - transfer.energy;
            let mass_residual_kg = sum(species_residual_kg[..n].iter().copied());
            let absolute_energy_residual_j =
                absolute_energy(gas, n, &species_residual_kg, energy_residual_j)?;
            if !mass_residual_kg.is_finite()
                || !energy_residual_j.is_finite()
                || species_residual_kg[..n]
                    .iter()
                    .any(|value| !value.is_finite())
            {
                return Err(Error::Invalid(
                    "unrepresentable finite mixture port balance",
                ));
            }
            self.trial_primitive[slot] = Some(primitive);
            self.trial_reports[slot] = Some(FiniteReservoirReport {
                before: before_inventory,
                after: after_inventory,
                transfer,
                species_residual_kg,
                mass_residual_kg,
                energy_residual_j,
                absolute_energy_residual_j,
            });
        }
        Ok(())
    }

    /// Infallible publication after the composing owner's last cancellation poll.
    pub(in crate::compressible) fn publish(&mut self) {
        for slot in 0..S {
            if self.accepted[slot].is_some() {
                std::mem::swap(&mut self.accepted[slot], &mut self.trial[slot]);
                self.primitive[slot] = self.trial_primitive[slot];
                self.reports[slot] = self.trial_reports[slot];
            }
        }
    }

    pub(in crate::compressible) fn try_clone(
        &self,
        gas: &MixtureGas,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let mut result = Self::empty();
        for slot in 0..S {
            if let Some(tank) = &self.accepted[slot] {
                result.bind(gas, slot, tank.volume_m3(), tank.state(), cancelled)?;
            }
        }
        result.reports = self.reports;
        result.trial_reports = self.reports;
        poll(cancelled)?;
        Ok(result)
    }

    /// Restore a prepared split/retry owner without allocating. The source's
    /// accepted state and report replace all destinations together after a poll.
    pub(in crate::compressible) fn copy_accepted_from(
        &mut self,
        gas: &MixtureGas,
        source: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        for slot in 0..S {
            if self.has_side(slot) != source.has_side(slot) {
                return Err(Error::Invalid(
                    "prepared finite mixture reservoir bindings differ",
                ));
            }
        }
        let zero = [0.0; MAX_SPECIES];
        for slot in 0..S {
            if let Some(tank) = source.accepted[slot].as_ref() {
                tank.exchange_candidate_into(
                    gas,
                    &zero[..gas.species_count()],
                    0.0,
                    self.trial[slot]
                        .as_mut()
                        .expect("prepared finite mixture output"),
                    &mut *cancelled,
                )?;
                self.trial_primitive[slot] = source.primitive[slot];
                self.trial_reports[slot] = source.reports[slot];
            }
        }
        poll(cancelled)?;
        self.publish();
        Ok(())
    }
}

fn absolute_energy(
    gas: &MixtureGas,
    species_count: usize,
    masses: &[f64; MAX_SPECIES],
    energy: f64,
) -> Result<f64, Error> {
    if species_count != gas.species_count() {
        return Err(Error::Invalid(
            "finite mixture inventory/source species axis",
        ));
    }
    let value = sum(std::iter::once(energy).chain(
        masses[..species_count]
            .iter()
            .zip(gas.energy_references())
            .map(|(mass, reference)| mass * reference),
    ));
    if !value.is_finite() {
        return Err(Error::Invalid(
            "unrepresentable finite mixture absolute energy",
        ));
    }
    Ok(value)
}

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let y = value - correction;
        let next = total + y;
        correction = (next - total) - y;
        total = next;
    }
    total
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
