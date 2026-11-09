//! Shared finite, well-mixed mixture storage for Cartesian pressure ports.
//!
//! Dimensions only supply the native side-ordered face impulses. The original
//! mixture plenum owns the extensive species/energy update and caloric recovery.

use fs_gas::multicomponent::network::Plenum;
use fs_gas::multicomponent::{Cell, Conserved, Error, MAX_SPECIES, MixtureGas, Primitive};
use fs_gas::multicomponent::wall_exchange::{WallCellReport, WallControl, WallThermal};

mod heat;
pub use heat::FiniteReservoirHeatReport;
pub(in crate::compressible) use heat::FiniteStage;

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
    /// Original thermal-source receipt when this tank has an authored wall/bath.
    pub heat: Option<FiniteReservoirHeatReport>,
    /// Per-species inventory change less the native signed face transfer [kg].
    pub species_residual_kg: [f64; MAX_SPECIES],
    pub mass_residual_kg: f64,
    /// Referenced internal-energy change less port energy and authored wall heat [J].
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
    // Each bound side prepares one slot in each bank. A report contains four
    // maximum-species arrays, so keeping even empty Options inline would make
    // every accepted/trial flow owner tens of kilobytes larger on the stack.
    // Length zero means no accepted report yet; length one stores the report.
    reports: [Vec<FiniteReservoirReport<D>>; S],
    trial_reports: [Vec<FiniteReservoirReport<D>>; S],
    allocated_bytes: usize,
    // Only heated sides allocate a third inventory. It is the unpublished
    // transport start, then is reused for the final thermal endpoint.
    heat_base: [Option<Plenum>; S],
    heat_state: [Option<WallThermal>; S],
    trial_heat_state: [Option<WallThermal>; S],
    heat_controls: [Option<WallControl>; S],
    trial_heat_reports: [Option<WallCellReport>; S],
}
impl<const D: usize, const S: usize> FiniteReservoirs<D, S> {
    pub(in crate::compressible) fn empty() -> Self {
        Self {
            accepted: std::array::from_fn(|_| None),
            trial: std::array::from_fn(|_| None),
            primitive: [None; S],
            trial_primitive: [None; S],
            reports: std::array::from_fn(|_| Vec::new()),
            trial_reports: std::array::from_fn(|_| Vec::new()),
            allocated_bytes: 0,
            heat_base: std::array::from_fn(|_| None),
            heat_state: [None; S],
            trial_heat_state: [None; S],
            heat_controls: [None; S],
            trial_heat_reports: [None; S],
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
            .checked_mul(
                2 * gas.species_count() * size_of::<f64>()
                    + 2 * size_of::<FiniteReservoirReport<D>>(),
            )
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
        // Prepare both report buffers before publishing any binding. Every
        // later stage, retry and accepted-state copy reuses this capacity.
        let mut reports = Vec::<FiniteReservoirReport<D>>::new();
        reports
            .try_reserve_exact(1)
            .map_err(|_| Error::Allocation)?;
        let mut trial_reports = Vec::<FiniteReservoirReport<D>>::new();
        trial_reports
            .try_reserve_exact(1)
            .map_err(|_| Error::Allocation)?;
        let report_bytes = reports
            .capacity()
            .checked_add(trial_reports.capacity())
            .and_then(|count| count.checked_mul(size_of::<FiniteReservoirReport<D>>()))
            .ok_or(Error::Invalid(
                "Cartesian finite mixture workspace overflow",
            ))?;
        let bytes = self
            .allocated_bytes
            .checked_add(accepted_bytes)
            .and_then(|bytes| bytes.checked_add(trial_bytes))
            .and_then(|bytes| bytes.checked_add(report_bytes))
            .ok_or(Error::Invalid(
                "Cartesian finite mixture workspace overflow",
            ))?;
        poll(cancelled)?;
        self.accepted[slot] = Some(accepted);
        self.trial[slot] = Some(trial);
        self.primitive[slot] = Some(primitive);
        self.trial_primitive[slot] = Some(primitive);
        self.reports[slot] = reports;
        self.trial_reports[slot] = trial_reports;
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
        self.reports[slot].first()
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
            let start = self.heat_base[slot].as_ref().unwrap_or(before);
            if self.heat_state[slot].is_some() && self.trial_heat_reports[slot].is_none() {
                return Err(Error::Invalid("finite tank heat stage was not prepared"));
            }
            start.exchange_candidate_into(
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
            let heat_to_gas_j = self.trial_heat_reports[slot].map_or(0.0, |q| q.heat_to_gas_j);
            let energy_residual_j =
                after_inventory.energy - before_inventory.energy - transfer.energy - heat_to_gas_j;
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
            let heat = self.heat_receipt(slot, energy_residual_j);
            self.trial_primitive[slot] = Some(primitive);
            self.trial_reports[slot].clear();
            self.trial_reports[slot].push(FiniteReservoirReport {
                before: before_inventory,
                after: after_inventory,
                transfer,
                heat,
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
                self.heat_state[slot] = self.trial_heat_state[slot];
                std::mem::swap(&mut self.reports[slot], &mut self.trial_reports[slot]);
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
                if let Some(thermal) = self.heat_state[slot] {
                    result.bind_heat(gas, slot, thermal, self.heat_controls[slot].expect("bound heat control"), cancelled)?;
                }
                result.reports[slot].extend_from_slice(&self.reports[slot]);
                result.trial_reports[slot].extend_from_slice(&self.reports[slot]);
            }
        }
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
            if self.has_side(slot) != source.has_side(slot)
                || self.heat_state[slot].is_some() != source.heat_state[slot].is_some()
            {
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
                self.trial_heat_state[slot] = source.heat_state[slot];
                self.trial_reports[slot].clear();
                self.trial_reports[slot].extend_from_slice(&source.reports[slot]);
            }
        }
        poll(cancelled)?;
        self.heat_controls = source.heat_controls;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_empty_finite_storage_keeps_large_species_reports_off_the_owner_stack() {
        // An unused finite reservoir is embedded in every native field and
        // both molecular split frames. Optional storage must stay small even
        // when the public report retains the complete maximum species axis.
        assert!(size_of::<FiniteReservoirs<3, 6>>() < 4096);
        assert!(size_of::<FiniteReservoirs<2, 4>>() < 4096);
        let storage = FiniteReservoirs::<3, 6>::empty();
        assert!(!storage.has_any());
        assert_eq!(
            storage.workspace_bytes(),
            size_of::<FiniteReservoirs<3, 6>>()
        );
        for slot in 0..6 {
            assert_eq!(storage.reports[slot].capacity(), 0);
            assert_eq!(storage.trial_reports[slot].capacity(), 0);
            assert!(storage.report(slot).is_none());
        }
    }

    #[test]
    fn g2_finite_report_clone_and_cancelled_restore_reuse_prepared_storage() {
        let gas = crate::compressible::mixture::tests::standard(0.2);
        let initial = gas.conserved(&[0.6, 0.4, 0.0], 0.0, 1.0).unwrap();
        let slot = 3;
        let mut owner = FiniteReservoirs::<3, 6>::empty();
        owner
            .bind(&gas, slot, 0.7, initial.as_cell(), &mut || false)
            .unwrap();
        let bytes = owner.workspace_bytes();
        assert_eq!(
            bytes,
            FiniteReservoirs::<3, 6>::required_bytes(&gas, 1).unwrap()
        );
        let faces = [0, 0, 0, 1, 0, 0];
        let transfer = [0.001, -0.0003, 0.0, 0.005, -0.003, 0.004, 0.002];
        owner
            .stage(&gas, &transfer, faces, 1.0, &mut || false)
            .unwrap();
        assert!(owner.report(slot).is_none());
        owner.publish();
        let mut target = owner.try_clone(&gas, &mut || false).unwrap();
        assert_eq!(target.report(slot), owner.report(slot));
        assert_eq!(target.workspace_bytes(), bytes);
        let buffers = |storage: &FiniteReservoirs<3, 6>| {
            let mut values = [
                storage.reports[slot].as_ptr(),
                storage.trial_reports[slot].as_ptr(),
            ];
            values.sort_unstable();
            values
        };
        let prepared = buffers(&target);
        target
            .stage(&gas, &transfer, faces, 2.0, &mut || false)
            .unwrap();
        target.publish();
        let previous = *target.report(slot).unwrap();
        let prior_state = target.inventory(slot).unwrap().state();
        let prior_densities = prior_state.partial_densities.to_vec();
        let prior_energy = prior_state.energy;
        let mut probe = target.try_clone(&gas, &mut || false).unwrap();
        let mut calls = 0;
        probe
            .copy_accepted_from(&gas, &owner, &mut || {
                calls += 1;
                false
            })
            .unwrap();
        let mut seen = 0;
        assert_eq!(
            target.copy_accepted_from(&gas, &owner, &mut || {
                seen += 1;
                seen == calls
            }),
            Err(Error::Cancelled),
        );
        assert_eq!(target.report(slot), Some(&previous));
        assert_eq!(
            target.inventory(slot).unwrap().state().partial_densities,
            prior_densities
        );
        assert_eq!(target.inventory(slot).unwrap().state().energy, prior_energy);
        assert_eq!(buffers(&target), prepared);
        assert_eq!(target.workspace_bytes(), bytes);
        target
            .copy_accepted_from(&gas, &owner, &mut || false)
            .unwrap();
        assert_eq!(target.report(slot), owner.report(slot));
        assert_eq!(
            target.inventory(slot).unwrap().state(),
            owner.inventory(slot).unwrap().state()
        );
        assert_eq!(buffers(&target), prepared);
        assert_eq!(target.workspace_bytes(), bytes);
        for side in 0..6 {
            if side != slot {
                assert_eq!(target.reports[side].capacity(), 0);
                assert_eq!(target.trial_reports[side].capacity(), 0);
            }
        }
    }
}
