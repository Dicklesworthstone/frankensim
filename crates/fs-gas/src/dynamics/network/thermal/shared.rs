//! One original finite wall inventory shared by several gas control volumes.

use super::*;
use crate::multicomponent::wall_exchange::shared::{Group, Work};
pub use crate::multicomponent::wall_exchange::shared::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallReport,
};
use crate::multicomponent::wall_exchange::solve_scalar;
use std::cell::Cell;

mod conduction;
use conduction::WallLinks;
pub use conduction::{SharedWallLink, SharedWallLinkReport};

pub(super) struct SharedWalls<'law> {
    groups: Vec<Group<'law>>,
    controls: SharedWallControl,
    links: Option<WallLinks>,
    branch_correction: Vec<f64>,
    candidate_report: Option<SharedWallReport>,
    last_report: Option<SharedWallReport>,
}

impl<'law> NetworkWallOperator<'law> {
    /// Bind finite shared walls before the first accepted source interval.
    ///
    /// Referenced branch cells must already have original branch bindings.
    /// Each selected local contact must be adiabatic with no finite storage;
    /// its original drag remains active and thermalizes kinetic energy FIRST.
    /// Duplicate targets, overlapping thermal owners, and reservoir storage
    /// refuse. Initial shared-wall temperatures must lie in the gas source's
    /// declared interval. No temperature or source support is extrapolated.
    ///
    /// Every group closes its simultaneous backward-Euler endpoint equations,
    /// using one wall-energy increment unknown and the original gas solver.
    /// Transport/source composition remains first order. Optional finite wall
    /// links compose afterward through [`Self::with_shared_wall_links`].
    /// Geometric spatial conduction, changing geometry and multiple shared
    /// walls on one gas cell are not inferred by this uniform-cell reduction.
    pub fn with_shared_walls(
        mut self,
        network: &GasNetwork1d,
        walls: Vec<SharedWall<'law>>,
        controls: SharedWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        self.validate_binding(network)?;
        if self.last_report.is_some()
            || self.shared.is_some()
            || walls.is_empty()
            || controls.maximum_source_evaluations == 0
        {
            return Err(Error::Invalid("shared wall binding state or work budget"));
        }
        let mut used = allocated(self.contacts.len(), false)?;
        let mut groups = Vec::new();
        groups
            .try_reserve_exact(walls.len())
            .map_err(|_| Error::Allocation)?;
        for mut wall in walls {
            poll(&mut cancelled)?;
            let contact = admit_contact(
                WallContact {
                    drag: WallDrag::None,
                    thermal: wall.storage.thermal,
                },
                wall.storage.heat_storage.as_ref(),
            )?;
            let conductance = match contact.thermal {
                WallThermal::Lumped {
                    conductance_w_per_k,
                    ..
                }
                | WallThermal::Caloric {
                    conductance_w_per_k,
                    ..
                } => conductance_w_per_k,
                _ => return Err(Error::Invalid("shared wall needs one finite inventory")),
            };
            if conductance != 0.0 || wall.ports.is_empty() {
                return Err(Error::Invalid("shared wall films belong to nonempty ports"));
            }
            wall.storage.thermal = contact.thermal;
            let temperature = storage_temperature(&wall.storage)?;
            if !(self.gas.domain().minimum_temperature..=self.gas.domain().maximum_temperature)
                .contains(&temperature)
            {
                return Err(Error::Invalid(
                    "shared wall outside original gas temperature interval",
                ));
            }
            let n = wall.ports.len();
            let mut volumes = allocated(n, 0.0)?;
            for (i, port) in wall.ports.iter().enumerate() {
                poll(&mut cancelled)?;
                let index = match port.gas {
                    SharedWallGas::Plenum(p) if p < self.plenum_count => p,
                    SharedWallGas::Branch { branch, cell } => {
                        let binding =
                            self.branches.get(branch).and_then(Option::as_ref).ok_or(
                                Error::Invalid("shared wall needs original branch binding"),
                            )?;
                        if cell >= binding.mesh.cell_count() {
                            return Err(Error::Invalid("shared wall branch cell"));
                        }
                        binding.start + cell
                    }
                    _ => return Err(Error::Invalid("shared wall plenum")),
                };
                if used[index]
                    || self.contacts[index].thermal != WallThermal::Adiabatic
                    || self.heat_storage[index].is_some()
                    || !port.conductance_w_per_k.is_finite()
                    || port.conductance_w_per_k < 0.0
                {
                    return Err(Error::Invalid(
                        "shared wall duplicate, overlapping source or conductance",
                    ));
                }
                used[index] = true;
                volumes[i] = self.volumes[index];
            }
            groups.push(Group::new(wall, volumes)?);
        }
        self.shared = Some(SharedWalls {
            groups,
            controls,
            links: None,
            branch_correction: allocated(network.branches.len(), 0.0)?,
            candidate_report: None,
            last_report: None,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted single wall inventory, original source and its physical ports.
    pub fn shared_wall(&self, index: usize) -> Option<&SharedWall<'law>> {
        self.shared
            .as_ref()?
            .groups
            .get(index)
            .map(|group| &group.accepted)
    }

    /// Accepted coupled equations; failed whole-network intervals retain it.
    pub fn last_shared_wall_report(&self) -> Option<SharedWallReport> {
        self.shared.as_ref()?.last_report
    }
}

fn storage_temperature(storage: &PlenumThermalContact<'_>) -> Result<f64, Error> {
    match &storage.heat_storage {
        Some(source) => Ok(source.temperature_k(storage.thermal)?),
        None => match storage.thermal {
            WallThermal::Lumped { temperature_k, .. } => Ok(temperature_k),
            _ => Err(Error::Invalid("shared wall needs original finite calorics")),
        },
    }
}

fn gas_state<'a>(
    location: SharedWallGas,
    branches: &'a mut [Euler1d],
    plenums: &'a mut [Conserved],
) -> &'a mut Conserved {
    match location {
        SharedWallGas::Plenum(p) => &mut plenums[p],
        SharedWallGas::Branch { branch, cell } => &mut branches[branch].work.candidate[cell],
    }
}

impl SharedWalls<'_> {
    pub(super) fn prepare(
        &mut self,
        gas: &FrozenGas,
        branches: &mut [Euler1d],
        plenums: &mut [Conserved],
        branch_sources: &mut [Totals],
        duration: Time,
        control: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, Error> {
        control.validate()?;
        let work = Work {
            used: Cell::new(0),
            contacts: Cell::new(0),
            limit: self.controls.maximum_source_evaluations,
        };
        let mut total = HeatSum::default();
        let mut maximum_wall = 0.0_f64;
        let mut iterations = 0_usize;
        self.branch_correction.fill(0.0);
        for group in &mut self.groups {
            poll(cancelled)?;
            for (i, port) in group.accepted.ports.iter().enumerate() {
                let state = gas_state(port.gas, branches, plenums);
                group.initial[i] = [state.momentum, state.energy];
            }
            let original = |location, temperature| {
                let density = match location {
                    SharedWallGas::Plenum(p) => plenums[p].density,
                    SharedWallGas::Branch { branch, cell } => {
                        branches[branch].work.candidate[cell].density
                    }
                };
                FrozenCaloric {
                    gas,
                    density,
                    hint: Some(temperature),
                }
            };
            let initial = group.accepted.storage.thermal;
            let report = group.solve(&original, duration, control, &work, initial, cancelled)?;
            maximum_wall = maximum_wall.max(report.energy_residual_j.abs());
            iterations = iterations
                .checked_add(report.iterations)
                .ok_or(Error::Allocation)?;
            total.add(report)?;
            for (i, port) in group.accepted.ports.iter().enumerate() {
                poll(cancelled)?;
                gas_state(port.gas, branches, plenums).energy = group.energies[i];
                if let SharedWallGas::Branch { branch, .. } = port.gas {
                    let increment = group.heat[i] - self.branch_correction[branch];
                    let next = branch_sources[branch].energy + increment;
                    if !next.is_finite() {
                        return Err(Error::Invalid("shared wall branch heat"));
                    }
                    self.branch_correction[branch] =
                        (next - branch_sources[branch].energy) - increment;
                    branch_sources[branch].energy = next;
                }
            }
        }
        if let Some(links) = &mut self.links {
            for (link, receipt) in links.ports.iter().zip(&mut links.candidate) {
                poll(cancelled)?;
                *receipt =
                    conduction::exchange(&self.groups, *link, duration, control, &work, cancelled)?;
                for i in 0..2 {
                    self.groups[link.walls[i]].candidate = receipt.endpoint_state[i];
                }
                maximum_wall = maximum_wall.max(receipt.energy_residual_j.abs());
                iterations = iterations
                    .checked_add(receipt.iterations)
                    .ok_or(Error::Allocation)?;
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
            }
        }
        poll(cancelled)?;
        let report = total.report(duration);
        self.candidate_report = Some(SharedWallReport {
            elapsed: duration,
            energy_residual_j: report.energy_residual_j,
            maximum_wall_residual_j: maximum_wall,
            maximum_contact_residual_j: report.maximum_thermal_residual_j,
            source_evaluations: work.used.get(),
            wall_iterations: iterations,
            contact_solves: work.contacts.get(),
        });
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
            group.accepted.storage.thermal = group.candidate;
        }
        self.last_report = self.candidate_report;
        if let Some(links) = &mut self.links {
            links.accepted.clone_from(&links.candidate);
        }
    }
}
