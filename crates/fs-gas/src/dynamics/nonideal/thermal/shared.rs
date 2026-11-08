//! One original finite wall inventory coupled to multiple real-gas cells.

use super::*;
use crate::multicomponent::wall_exchange::WallDrag;
use crate::multicomponent::wall_exchange::shared::{Group, Work};
pub use crate::multicomponent::wall_exchange::shared::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallReport,
};
use std::cell::Cell;

pub(super) struct SharedWalls<'law> {
    groups: Vec<Group<'law>>,
    controls: SharedWallControl,
    candidate_report: Option<SharedWallReport>,
    last_report: Option<SharedWallReport>,
}

impl<'law> RealGasWallOperator<'law> {
    /// Bind one finite wall amount to several cells on this original mesh.
    /// Ports use `SharedWallGas::Branch { branch: 0, cell }`; plenums and other
    /// branches are not part of this single-field owner. Every gas cell belongs
    /// to at most one shared wall and must have an adiabatic local contact with
    /// no local heat storage. Its original drag still acts before heat exchange.
    ///
    /// Each group solves all backward-Euler contact equations against ONE wall
    /// energy increment. Every gas retains its current candidate density and
    /// original NASA/PR Delta u and Cv. Authored constant wall capacity or the
    /// original nonlinear/latent h/u law supplies wall storage exactly once.
    /// Wall source temperature must lie in the declared gas caloric interval.
    ///
    /// Bind before the first accepted wall interval. The outer flow owner
    /// publishes gas, wall, reports and clock only after all candidate cells
    /// pass complete source, transport and acoustic admission. This is a
    /// uniform-wall reduction with authored conductances, not a spatial wall
    /// conduction model or an inferred dense-gas heat-transfer correlation.
    pub fn with_shared_walls(
        mut self,
        walls: Vec<SharedWall<'law>>,
        controls: SharedWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        poll(&mut cancelled)?;
        if self.last_report.is_some()
            || self.shared.is_some()
            || walls.is_empty()
            || controls.maximum_source_evaluations == 0
        {
            return Err(RealGasFlowError::Invalid("shared real-gas wall binding state or work budget"));
        }
        let mut used = allocated(self.contacts.len(), false)?;
        let mut groups = Vec::new();
        groups.try_reserve_exact(walls.len()).map_err(|_| crate::dynamics::Error::Allocation)?;
        for mut wall in walls {
            poll(&mut cancelled)?;
            let mut contact = WallContact { drag: WallDrag::None, thermal: wall.storage.thermal };
            contact.validate()?;
            let temperature = match &wall.storage.heat_storage {
                Some(source) => {
                    contact = source.bind_contact(contact)?;
                    source.temperature_k(contact.thermal)?
                }
                None => match contact.thermal {
                    WallThermal::Lumped { temperature_k, .. } => temperature_k,
                    _ => return Err(RealGasFlowError::Invalid("shared real-gas wall needs finite storage")),
                },
            };
            let conductance = match contact.thermal {
                WallThermal::Lumped { conductance_w_per_k, .. }
                | WallThermal::Caloric { conductance_w_per_k, .. } => conductance_w_per_k,
                _ => return Err(RealGasFlowError::Invalid("shared real-gas wall needs finite storage")),
            };
            if conductance != 0.0 || wall.ports.is_empty() {
                return Err(RealGasFlowError::Invalid("shared real-gas wall films belong to nonempty ports"));
            }
            if !self.gas.domain.admits_temperature(temperature) {
                return Err(RealGasFlowError::Invalid("shared wall outside original real-gas temperature interval"));
            }
            wall.storage.thermal = contact.thermal;
            let mut volumes = allocated(wall.ports.len(), 0.0)?;
            for (i, port) in wall.ports.iter().enumerate() {
                poll(&mut cancelled)?;
                let cell = match port.gas {
                    SharedWallGas::Branch { branch: 0, cell } if cell < self.contacts.len() => cell,
                    _ => return Err(RealGasFlowError::Invalid("shared real-gas wall mesh cell")),
                };
                if used[cell]
                    || self.contacts[cell].thermal != WallThermal::Adiabatic
                    || self.heat_storage[cell].is_some()
                    || !port.conductance_w_per_k.is_finite()
                    || port.conductance_w_per_k < 0.0
                {
                    return Err(RealGasFlowError::Invalid("shared real-gas wall duplicate, overlapping source or conductance"));
                }
                used[cell] = true;
                volumes[i] = self.mesh.volumes()[cell];
            }
            groups.push(Group::new(wall, volumes)?);
        }
        self.shared = Some(SharedWalls {
            groups, controls, candidate_report: None, last_report: None,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted single wall inventory, source and its unchanged physical ports.
    pub fn shared_wall(&self, index: usize) -> Option<&SharedWall<'law>> {
        self.shared.as_ref()?.groups.get(index).map(|group| &group.accepted)
    }

    /// Last accepted shared storage/contact equations and actual bounded work.
    pub fn last_shared_wall_report(&self) -> Option<SharedWallReport> {
        self.shared.as_ref()?.last_report
    }
}

fn cell(location: SharedWallGas) -> usize {
    match location {
        SharedWallGas::Branch { branch: 0, cell } => cell,
        _ => unreachable!("original real-gas mesh binding admitted"),
    }
}

impl SharedWalls<'_> {
    pub(super) fn prepare(
        &mut self,
        gas: &RealGas,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, RealGasFlowError> {
        controls.validate()?;
        let work = Work {
            used: Cell::new(0), contacts: Cell::new(0),
            limit: self.controls.maximum_source_evaluations,
        };
        let mut total = HeatSum::default();
        let mut maximum_wall = 0.0_f64;
        let mut iterations = 0_usize;
        for group in &mut self.groups {
            poll(cancelled)?;
            for (i, port) in group.accepted.ports.iter().enumerate() {
                let state = states[cell(port.gas)];
                group.initial[i] = [state.momentum, state.energy];
            }
            let original = |location, temperature| RealGasCaloric {
                gas, density: states[cell(location)].density, hint: Some(temperature),
            };
            let initial = group.accepted.storage.thermal;
            let report = group.solve(&original, duration, controls, &work, initial, cancelled)?;
            maximum_wall = maximum_wall.max(report.energy_residual_j.abs());
            iterations = iterations.checked_add(report.iterations)
                .ok_or(RealGasFlowError::Invalid("shared real-gas wall iteration count overflow"))?;
            total.add(report)?;
            for (i, port) in group.accepted.ports.iter().enumerate() {
                poll(cancelled)?;
                let index = cell(port.gas);
                let candidate = Conserved { energy: group.energies[i], ..states[index] };
                // The scalar heat owner can reuse a temperature hint. Complete
                // material/transport/acoustic admission is unconditional here.
                work.charge()?;
                primitives[index] = gas.decode(candidate, cancelled)?;
                states[index] = candidate;
            }
        }
        poll(cancelled)?;
        let report = total.report(duration);
        self.candidate_report = Some(SharedWallReport {
            elapsed: duration,
            energy_residual_j: report.source_storage_residual_j,
            maximum_wall_residual_j: maximum_wall,
            maximum_contact_residual_j: report.maximum_thermal_residual_j,
            source_evaluations: work.used.get(),
            wall_iterations: iterations,
            contact_solves: work.contacts.get(),
        });
        Ok(WallCellReport {
            elapsed: duration, wall_impulse_ns: 0.0, mechanical_work_to_gas_j: 0.0,
            drag_heating_j: 0.0, heat_to_gas_j: report.heat_to_gas_j,
            wall_energy_change_j: report.wall_storage_change_j,
            external_heat_to_gas_j: 0.0,
            energy_residual_j: report.source_storage_residual_j,
            thermal_residual_j: report.maximum_thermal_residual_j,
            iterations: report.maximum_iterations,
        })
    }

    pub(super) fn commit(&mut self) {
        for group in &mut self.groups {
            group.accepted.storage.thermal = group.candidate;
        }
        self.last_report = self.candidate_report;
    }
}
