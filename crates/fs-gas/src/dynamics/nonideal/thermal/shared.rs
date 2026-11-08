//! One original finite wall inventory coupled to multiple real-gas cells.

use super::*;
use crate::multicomponent::network::wall_exchange::MaterialWallConduction;
use crate::multicomponent::wall_exchange::WallDrag;
use crate::multicomponent::wall_exchange::shared::{Group, Work};
pub use crate::multicomponent::wall_exchange::shared::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallReport,
};
use std::cell::Cell;

mod conduction;
use conduction::WallLinks;
pub use conduction::{SharedWallLink, SharedWallLinkReport};

pub(super) struct SharedWalls<'law> {
    groups: Vec<Group<'law>>,
    controls: SharedWallControl,
    links: Option<WallLinks>,
    candidate_report: Option<SharedWallReport>,
    last_report: Option<SharedWallReport>,
    midpoint: Vec<WallThermal>,
    first_report: Option<SharedWallReport>,
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
    /// Gas and wall retain independent caloric intervals. Constant or lagged
    /// conductance admits an initially hotter/colder wall when the coupled gas
    /// endpoint is supported. Implicit conductivity and integrated profiles
    /// require initial wall temperature inside every active gas interval.
    ///
    /// Bind before the first accepted wall interval. The outer flow owner
    /// publishes gas, wall, reports and clock only after all candidate cells
    /// pass complete source, transport and acoustic admission. This is a
    /// uniform storage reduction within each wall region. Material resistance
    /// and quasi-steady profiles bind through `with_shared_wall_conduction`;
    /// distinct regions exchange heat through `with_shared_wall_links`.
    /// No transient spatial PDE or dense-gas film correlation is inferred.
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
            return Err(RealGasFlowError::Invalid(
                "shared real-gas wall binding state or work budget",
            ));
        }
        let mut used = allocated(self.contacts.len(), false)?;
        let mut groups = Vec::new();
        groups
            .try_reserve_exact(walls.len())
            .map_err(|_| crate::dynamics::Error::Allocation)?;
        for mut wall in walls {
            poll(&mut cancelled)?;
            let mut contact = WallContact {
                drag: WallDrag::None,
                thermal: wall.storage.thermal,
            };
            contact.validate()?;
            match &wall.storage.heat_storage {
                Some(source) => {
                    contact = source.bind_contact(contact)?;
                    source.temperature_k(contact.thermal)?
                }
                None => match contact.thermal {
                    WallThermal::Lumped { temperature_k, .. } => temperature_k,
                    _ => {
                        return Err(RealGasFlowError::Invalid(
                            "shared real-gas wall needs finite storage",
                        ));
                    }
                },
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
                    return Err(RealGasFlowError::Invalid(
                        "shared real-gas wall needs finite storage",
                    ));
                }
            };
            if conductance != 0.0 || wall.ports.is_empty() {
                return Err(RealGasFlowError::Invalid(
                    "shared real-gas wall films belong to nonempty ports",
                ));
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
                    return Err(RealGasFlowError::Invalid(
                        "shared real-gas wall duplicate, overlapping source or conductance",
                    ));
                }
                used[cell] = true;
                volumes[i] = self.mesh.volumes()[cell];
            }
            groups.push(Group::new(wall, volumes)?);
        }
        let midpoint = allocated(groups.len(), WallThermal::Adiabatic)?;
        self.shared = Some(SharedWalls {
            groups,
            controls,
            links: None,
            candidate_report: None,
            last_report: None,
            midpoint,
            first_report: None,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted single wall inventory, source and its unchanged physical ports.
    pub fn shared_wall(&self, index: usize) -> Option<&SharedWall<'law>> {
        self.shared
            .as_ref()?
            .groups
            .get(index)
            .map(|group| &group.accepted)
    }

    /// Last accepted shared storage/contact equations and actual bounded work.
    pub fn last_shared_wall_report(&self) -> Option<SharedWallReport> {
        self.shared.as_ref()?.last_report
    }

    /// Bind original solid conductivity and resistance geometry to gas ports.
    /// Rows follow shared-wall order; entries follow each wall's port order.
    /// A bound port's authored G remains its gas-side film. Its material layer
    /// adds solid resistance once, without another mass, capacity or latent
    /// inventory. `None` retains the port's complete authored conductance.
    /// Bind after shared walls and before the first accepted source interval.
    ///
    /// Initial-state layers sample the original wall state at each source
    /// stage; implicit layers follow its proposed nonlinear or latent storage.
    /// An integrated profile instead uses the original conductivity integral
    /// between the wall trace and the solved gas-facing surface, with exact
    /// admitted planar or radial resistance geometry. The existing caloric
    /// owner still holds all storage; transient internal gradients require
    /// spatial refinement. Phase-dependent k requires its matching phase chart.
    ///
    /// Every runtime source query spends the shared attempt's work allowance,
    /// across both thermal halves of a piston step. A later wall link must
    /// leave every retained conductivity source and material face admissible.
    pub fn with_shared_wall_conduction(
        mut self,
        layers: Vec<Vec<Option<MaterialWallConduction>>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        poll(&mut cancelled)?;
        let shared = self.shared.as_mut().ok_or(RealGasFlowError::Invalid(
            "real-gas material conduction requires shared walls",
        ))?;
        if self.last_report.is_some()
            || layers.len() != shared.groups.len()
            || shared
                .groups
                .iter()
                .any(|group| group.conduction_layers().is_some())
            || !layers.iter().flatten().any(Option::is_some)
        {
            return Err(RealGasFlowError::Invalid(
                "shared real-gas wall conduction state or topology",
            ));
        }
        for (group, layers) in shared.groups.iter_mut().zip(layers) {
            group.bind_conduction(layers, &mut cancelled)?;
        }
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original material, resistance geometry and time evaluation per gas port.
    pub fn shared_wall_conduction(&self, wall: usize) -> Option<&[Option<MaterialWallConduction>]> {
        self.shared.as_ref()?.groups.get(wall)?.conduction_layers()
    }
}

fn cell(location: SharedWallGas) -> usize {
    match location {
        SharedWallGas::Branch { branch: 0, cell } => cell,
        _ => unreachable!("original real-gas mesh binding admitted"),
    }
}

impl SharedWalls<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        &mut self,
        gas: &RealGas,
        mesh: &Mesh1d,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, RealGasFlowError> {
        self.prepare_stage(
            gas, mesh, states, primitives, duration, controls, false, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_piston_first(
        &mut self,
        gas: &RealGas,
        mesh: &Mesh1d,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, RealGasFlowError> {
        let report = self.prepare_stage(
            gas, mesh, states, primitives, duration, controls, false, cancelled,
        )?;
        for (state, group) in self.midpoint.iter_mut().zip(&self.groups) {
            poll(cancelled)?;
            *state = group.candidate;
        }
        self.first_report = self.candidate_report;
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_piston_second(
        &mut self,
        gas: &RealGas,
        mesh: &Mesh1d,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, RealGasFlowError> {
        self.prepare_stage(
            gas, mesh, states, primitives, duration, controls, true, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_stage(
        &mut self,
        gas: &RealGas,
        mesh: &Mesh1d,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        continuation: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, RealGasFlowError> {
        controls.validate()?;
        let first = if continuation {
            Some(self.first_report.ok_or(RealGasFlowError::Invalid(
                "shared wall needs its first thermal half",
            ))?)
        } else {
            None
        };
        let work = Work {
            // Both source halves spend ONE attempt's original work allowance.
            used: Cell::new(first.map_or(0, |r| r.source_evaluations)),
            contacts: Cell::new(first.map_or(0, |r| r.contact_solves)),
            limit: self.controls.maximum_source_evaluations,
        };
        let mut total = HeatSum::default();
        let mut maximum_wall = first.map_or(0.0_f64, |r| r.maximum_wall_residual_j);
        let mut iterations = first.map_or(0_usize, |r| r.wall_iterations);
        for (index, group) in self.groups.iter_mut().enumerate() {
            poll(cancelled)?;
            for (i, port) in group.accepted.ports.iter().enumerate() {
                poll(cancelled)?;
                let state = states[cell(port.gas)];
                group.initial[i] = [state.momentum, state.energy];
                // These are private source scratch, refreshed on EVERY stage,
                // including a stationary call following a refused piston step.
                group.volumes[i] = mesh.volumes()[cell(port.gas)];
            }
            let original = |location, temperature| RealGasCaloric {
                gas,
                density: states[cell(location)].density,
                hint: Some(temperature),
            };
            let initial = if continuation {
                self.midpoint[index]
            } else {
                group.accepted.storage.thermal
            };
            let report = group.solve(&original, duration, controls, &work, initial, cancelled)?;
            maximum_wall = maximum_wall.max(report.energy_residual_j.abs());
            iterations =
                iterations
                    .checked_add(report.iterations)
                    .ok_or(RealGasFlowError::Invalid(
                        "shared real-gas wall iteration count overflow",
                    ))?;
            total.add(report)?;
            for (i, port) in group.accepted.ports.iter().enumerate() {
                poll(cancelled)?;
                let index = cell(port.gas);
                let candidate = Conserved {
                    energy: group.energies[i],
                    ..states[index]
                };
                // The scalar heat owner can reuse a temperature hint. Complete
                // material/transport/acoustic admission is unconditional here.
                work.charge()?;
                primitives[index] = gas.decode(candidate, cancelled)?;
                states[index] = candidate;
            }
        }
        if let Some(links) = &mut self.links {
            conduction::prepare(
                &mut self.groups,
                links,
                duration,
                controls,
                &work,
                &mut total,
                &mut maximum_wall,
                &mut iterations,
                cancelled,
            )?;
            // Later links can move a wall beyond an earlier gas-port or link
            // source's support. Admit the complete stage without transferring
            // heat again. This runs in both stationary and piston source halves.
            for group in &self.groups {
                poll(cancelled)?;
                group.admit_conduction_endpoint(group.candidate, true, &work, cancelled)?;
                let original = |location, temperature| RealGasCaloric {
                    gas,
                    density: states[cell(location)].density,
                    hint: Some(temperature),
                };
                group.admit_integrated_profile_endpoints(
                    &original,
                    group.candidate,
                    duration,
                    controls,
                    &work,
                    cancelled,
                )?;
            }
            links.admit_conduction_endpoints(|i| &self.groups[i], &work, cancelled)?;
        }
        poll(cancelled)?;
        let report = total.report(duration);
        self.candidate_report = Some(SharedWallReport {
            elapsed: Time::new(duration.value() + first.map_or(0.0, |r| r.elapsed.value())),
            energy_residual_j: report.source_storage_residual_j
                + first.map_or(0.0, |r| r.energy_residual_j),
            maximum_wall_residual_j: maximum_wall,
            maximum_contact_residual_j: report
                .maximum_thermal_residual_j
                .max(first.map_or(0.0, |r| r.maximum_contact_residual_j)),
            source_evaluations: work.used.get(),
            wall_iterations: iterations,
            contact_solves: work.contacts.get(),
        });
        Ok(WallCellReport {
            elapsed: duration,
            wall_impulse_ns: 0.0,
            mechanical_work_to_gas_j: 0.0,
            drag_heating_j: 0.0,
            heat_to_gas_j: report.heat_to_gas_j,
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
        if let Some(links) = &mut self.links {
            links.accepted.clone_from(&links.candidate);
        }
    }
}
