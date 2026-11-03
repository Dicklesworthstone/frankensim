//! Shared finite inventories on the existing mixture-network transaction.

use super::*;
use crate::multicomponent::wall_exchange::shared::{Group, Work, mixture};
pub use crate::multicomponent::wall_exchange::shared::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallReport,
    SharedWallStorage,
};
use std::cell::Cell as Counter;

mod links;
use links::WallLinks;
pub use links::{SharedWallLink, SharedWallLinkReport};

struct BoundGroup<'law> {
    storage: Group<'law>,
    indices: Vec<usize>,
}

pub(super) struct SharedWalls<'law> {
    groups: Vec<BoundGroup<'law>>,
    work: Work,
    links: Option<WallLinks>,
    candidate_report: Option<SharedWallReport>,
    last_report: Option<SharedWallReport>,
}

impl<'law> NetworkWallOperator<'law> {
    /// Attach distinct gas cells to one original finite storage per shared wall.
    ///
    /// Every port supplies its complete conductance [W/K]; the storage's own
    /// conductance must be zero. Selected local contacts must be adiabatic and
    /// have no competing finite storage. Their original drag remains active.
    /// Sharing heat preserves every gas cell's mass and species separately.
    ///
    /// Each source half solves all ports of a group simultaneously against
    /// one original wall h/u increment. Gas inversion and heat capacity use
    /// that half's CURRENT partial densities, including chemistry and transport
    /// between halves. The candidate wall carries through those same halves
    /// and publishes only with the complete network composition. Source work
    /// is bounded across both halves; a retry starts from accepted storage.
    ///
    /// These are spatially uniform walls with authored, constant conductances.
    /// Existing source/transport splitting remains first order. This binding
    /// does not infer through-wall gradients or temperature-dependent films.
    pub fn with_shared_walls(
        mut self,
        network: &MixtureNetwork1d,
        walls: Vec<SharedWall<'law>>,
        controls: SharedWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        self.admit_plenum_volumes(network, &mut cancelled)?;
        if self.last_report.is_some()
            || self.trial_active
            || self.shared.is_some()
            || walls.is_empty()
            || controls.maximum_source_evaluations == 0
            || self.gas != *network.gas()
            || self.branches.len() != network.branches.len()
            || self
                .branches
                .iter()
                .zip(&network.branches)
                .any(|(shape, branch)| shape.mesh != branch.mesh || branch.gas != self.gas)
        {
            return Err(WallExchangeError::Invalid(
                "shared wall binding or source ownership",
            ));
        }
        let mut used = allocated(self.contacts.len(), false)?;
        let mut groups = Vec::new();
        groups
            .try_reserve_exact(walls.len())
            .map_err(|_| WallExchangeError::Invalid("shared wall allocation"))?;
        for mut wall in walls {
            poll(&mut cancelled)?;
            let mut contact = WallContact {
                drag: WallDrag::None,
                thermal: wall.storage.thermal,
            };
            contact.validate()?;
            let temperature = if let Some(source) = &wall.storage.heat_storage {
                contact = source.bind_contact(contact)?;
                source.temperature_k(contact.thermal)?
            } else if let WallThermal::Lumped { temperature_k, .. } = contact.thermal {
                temperature_k
            } else {
                return Err(WallExchangeError::Invalid(
                    "shared wall needs original finite calorics",
                ));
            };
            let conductance = match contact.thermal {
                WallThermal::Lumped {
                    conductance_w_per_k,
                    ..
                }
                | WallThermal::Caloric {
                    conductance_w_per_k,
                    ..
                } => conductance_w_per_k,
                _ => {
                    return Err(WallExchangeError::Invalid(
                        "shared wall needs one finite inventory",
                    ));
                }
            };
            if conductance != 0.0
                || wall.ports.is_empty()
                || !(self.gas.domain().minimum_temperature..=self.gas.domain().maximum_temperature)
                    .contains(&temperature)
            {
                return Err(WallExchangeError::Invalid(
                    "shared wall films, ports or source interval",
                ));
            }
            wall.storage.thermal = contact.thermal;
            let mut indices = allocated(wall.ports.len(), 0)?;
            let mut volumes = allocated(wall.ports.len(), 0.0)?;
            for (i, port) in wall.ports.iter().enumerate() {
                poll(&mut cancelled)?;
                let (index, volume) = match port.gas {
                    SharedWallGas::Plenum(p) if p < self.plenum_volumes.len() => {
                        (self.plenum_start + p, self.plenum_volumes[p])
                    }
                    SharedWallGas::Branch { branch, cell } => {
                        let shape = self
                            .branches
                            .get(branch)
                            .ok_or(WallExchangeError::Invalid("shared wall branch"))?;
                        let volume = *shape
                            .mesh
                            .volumes()
                            .get(cell)
                            .ok_or(WallExchangeError::Invalid("shared wall cell"))?;
                        (shape.start + cell, volume)
                    }
                    _ => return Err(WallExchangeError::Invalid("shared wall plenum")),
                };
                if used[index]
                    || self.contacts[index].thermal != WallThermal::Adiabatic
                    || self.caloric[index].is_some()
                    || !port.conductance_w_per_k.is_finite()
                    || port.conductance_w_per_k < 0.0
                {
                    return Err(WallExchangeError::Invalid(
                        "shared wall duplicate, overlapping source or conductance",
                    ));
                }
                used[index] = true;
                indices[i] = index;
                volumes[i] = volume;
            }
            groups.push(BoundGroup {
                storage: Group::new(wall, volumes)?,
                indices,
            });
        }
        self.shared = Some(SharedWalls {
            groups,
            links: None,
            work: Work {
                used: Counter::new(0),
                contacts: Counter::new(0),
                limit: controls.maximum_source_evaluations,
            },
            candidate_report: None,
            last_report: None,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted original wall inventory and its ports; never a trial half state.
    pub fn shared_wall(&self, index: usize) -> Option<&SharedWall<'law>> {
        self.shared
            .as_ref()?
            .groups
            .get(index)
            .map(|group| &group.storage.accepted)
    }

    /// Accepted shared-source balances, including both composed wall halves.
    pub fn last_shared_wall_report(&self) -> Option<SharedWallReport> {
        self.shared.as_ref()?.last_report
    }
}

impl SharedWalls<'_> {
    pub(super) fn begin_trial(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        self.candidate_report = None;
        self.work.used.set(0);
        self.work.contacts.set(0);
        if let Some(links) = &mut self.links {
            links.candidate.clear();
        }
        for group in &mut self.groups {
            poll(cancelled)?;
            group.storage.candidate = group.storage.accepted.storage.thermal;
        }
        Ok(())
    }

    pub(super) fn prepare(
        &mut self,
        network: &MixtureNetwork1d,
        candidate_states: &mut [[f64; 2]],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, WallExchangeError> {
        let mut total = Balance::new();
        let mut maximum_wall = 0.0_f64;
        let mut iterations = 0_usize;
        for bound in &mut self.groups {
            poll(cancelled)?;
            let group = &mut bound.storage;
            for (i, (port, &index)) in group.accepted.ports.iter().zip(&bound.indices).enumerate() {
                poll(cancelled)?;
                group.initial[i] = candidate_states[index];
                // Only the composing volume-motion owner may change a trial
                // plenum volume; its current volume enters this same heat solve.
                if let SharedWallGas::Plenum(p) = port.gas {
                    group.volumes[i] = network.plenums[p].volume_m3;
                }
            }
            let original = |location, _temperature| {
                let densities = match location {
                    SharedWallGas::Plenum(p) => network.plenums[p].state().partial_densities,
                    SharedWallGas::Branch { branch, cell } => {
                        network.branches[branch]
                            .cell(cell)
                            .expect("bound original cell")
                            .partial_densities
                    }
                };
                mixture(network.gas(), densities)
            };
            let initial = group.candidate;
            let report = group.solve(
                &original, duration, controls, &self.work, initial, cancelled,
            )?;
            total.add(report)?;
            maximum_wall = maximum_wall.max(report.energy_residual_j.abs());
            iterations = iterations
                .checked_add(report.iterations)
                .ok_or(WallExchangeError::Invalid("shared wall iteration count"))?;
            for (&index, &energy) in bound.indices.iter().zip(&group.energies) {
                candidate_states[index][1] = energy;
            }
        }
        if let Some(links) = &mut self.links {
            // Prepared storage retains at most the two halves of one complete
            // network transaction. A retry clears it before any new solve.
            if links.candidate.len() > links.ports.len() {
                return Err(WallExchangeError::Invalid("wall link source-stage count"));
            }
            for &link in &links.ports {
                poll(cancelled)?;
                let receipt = crate::multicomponent::wall_exchange::shared::links::exchange(
                    link.walls.map(|i| &self.groups[i].storage),
                    link,
                    duration,
                    controls,
                    &self.work,
                    cancelled,
                )?;
                for i in 0..2 {
                    self.groups[link.walls[i]].storage.candidate = receipt.endpoint_state[i];
                }
                maximum_wall = maximum_wall.max(receipt.energy_residual_j.abs());
                iterations = iterations
                    .checked_add(receipt.iterations)
                    .ok_or(WallExchangeError::Invalid("wall link iteration count"))?;
                total.add(WallCellReport {
                    elapsed: duration,
                    wall_impulse_ns: 0.0,
                    drag_heating_j: 0.0,
                    heat_to_gas_j: 0.0,
                    wall_energy_change_j: receipt.energy_residual_j,
                    external_heat_to_gas_j: 0.0,
                    energy_residual_j: receipt.energy_residual_j,
                    thermal_residual_j: receipt.flux_residual_j,
                    iterations: receipt.iterations,
                })?;
                links.candidate.push(receipt);
            }
        }
        let report = total.report(duration);
        let mut shared_report = SharedWallReport {
            elapsed: duration,
            energy_residual_j: report.energy_residual_j,
            maximum_wall_residual_j: maximum_wall,
            maximum_contact_residual_j: report.maximum_thermal_residual_j,
            source_evaluations: self.work.used.get(),
            wall_iterations: iterations,
            contact_solves: self.work.contacts.get(),
        };
        if let Some(previous) = self.candidate_report {
            shared_report.elapsed = Time::new(previous.elapsed.value() + duration.value());
            shared_report.energy_residual_j += previous.energy_residual_j;
            shared_report.maximum_wall_residual_j =
                maximum_wall.max(previous.maximum_wall_residual_j);
            shared_report.maximum_contact_residual_j = shared_report
                .maximum_contact_residual_j
                .max(previous.maximum_contact_residual_j);
            shared_report.wall_iterations =
                iterations.checked_add(previous.wall_iterations).ok_or(
                    WallExchangeError::Invalid("composed shared wall iterations"),
                )?;
            if !shared_report.elapsed.value().is_finite()
                || !shared_report.energy_residual_j.is_finite()
            {
                return Err(WallExchangeError::Invalid("composed shared wall balance"));
            }
        }
        poll(cancelled)?;
        self.candidate_report = Some(shared_report);
        Ok(WallCellReport {
            elapsed: duration,
            wall_impulse_ns: 0.0,
            drag_heating_j: 0.0,
            heat_to_gas_j: report.heat_to_gas_j,
            wall_energy_change_j: report.wall_energy_change_j,
            external_heat_to_gas_j: 0.0,
            energy_residual_j: report.energy_residual_j,
            thermal_residual_j: report.maximum_thermal_residual_j,
            iterations: report.maximum_iterations,
        })
    }

    pub(super) fn commit(&mut self) {
        for group in &mut self.groups {
            group.storage.accepted.storage.thermal = group.storage.candidate;
        }
        self.last_report = self.candidate_report;
        if let Some(links) = &mut self.links {
            core::mem::swap(&mut links.accepted, &mut links.candidate);
        }
    }
}
