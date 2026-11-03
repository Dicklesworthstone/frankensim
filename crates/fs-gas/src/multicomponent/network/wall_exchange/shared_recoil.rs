//! One axial mechanical support shared by every cell of a selected branch.
//!
//! Mechanical inertia and thermal storage have separate identities. The shared
//! friction solve runs before the original caloric/thermal owner and deposits
//! slip loss exactly once through total energy minus the new kinetic energy.

use super::{
    Cell, MixtureGas, MixtureNetwork1d, NetworkWallOperator, SupportExchange, WallCellReport,
    WallContact, WallControl, WallDrag, WallExchangeError, WallThermal, accumulate, allocated,
    poll, positive, recoil, sum,
};
use crate::multicomponent::wall_exchange::{
    RecoilControl, SharedRecoilCell, SharedRecoilControl, SharedRecoilWorkspace, SlidingSupport,
    WallUpdate,
};
use fs_qty::Time;

pub(super) struct SharedBranchSupports {
    branches: Vec<Option<BranchSupport>>,
    control: SharedRecoilControl,
    momentum_tolerance_ns: f64,
}

struct BranchSupport {
    accepted: SlidingSupport,
    trial: SlidingSupport,
    exchange: Option<SupportExchange>,
    candidate_exchange: SupportExchange,
    workspace: SharedRecoilWorkspace,
    inputs: Vec<SharedRecoilCell>,
    motions: Vec<WallUpdate>,
}

impl NetworkWallOperator<'_> {
    /// Bind one finite axial support to all cells of each selected branch.
    ///
    /// Rows follow the original branch order. `None` retains the existing
    /// prescribed or independent-support path. A selected branch has ONE mass:
    /// its cells exchange momentum with that same support simultaneously to the
    /// shared solver's declared local time-error estimate. The original exact
    /// pair flows are composed with bounded, rate-resolved time refinement.
    /// An exhausted work budget refuses the whole transaction.
    ///
    /// Coefficients are frozen at each source half's entering state. Authored
    /// friction or original circular transport supplies those coefficients;
    /// actual support-relative endpoints retain the original material checks.
    /// All frictional slip heat stays in gas, then the original local, caloric,
    /// conduction or shared thermal owner performs heat exchange.
    ///
    /// Supports follow both original source halves and publish with gas and
    /// walls. Allocation occurs only at binding. Different branch axes never
    /// share or sum inertia. This does not move the mesh, represent a piston,
    /// add a flow clock, or infer a 3D body/torque. Independent finite supports
    /// may coexist only on other branches. Bind before evolution.
    pub fn with_shared_branch_supports(
        mut self,
        branches: Vec<Option<SlidingSupport>>,
        control: SharedRecoilControl,
        momentum_tolerance_ns: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        control.validate()?;
        if self.shared_supports.is_some()
            || self.last_report.is_some()
            || self.trial_active
            || branches.len() != self.branches.len()
            || !positive(momentum_tolerance_ns)
        {
            return Err(WallExchangeError::Invalid(
                "shared branch supports must bind once before evolution",
            ));
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(branches.len())
            .map_err(|_| WallExchangeError::Invalid("shared support allocation"))?;
        for (shape, support) in self.branches.iter().zip(branches) {
            poll(&mut cancelled)?;
            let Some(support) = support else {
                prepared.push(None);
                continue;
            };
            let count = shape.mesh.cell_count();
            if !positive(support.mass_kg) || !support.velocity_m_s.is_finite() {
                return Err(WallExchangeError::Invalid(
                    "shared support mass or velocity",
                ));
            }
            if self
                .supports
                .as_ref()
                .is_some_and(|supports| supports.any_selected(shape.start, count))
            {
                return Err(WallExchangeError::Invalid(
                    "independent and shared supports overlap",
                ));
            }
            let inputs = allocated(
                count,
                SharedRecoilCell {
                    density_kg_m3: 1.0,
                    momentum_kg_m2_s: 0.0,
                    volume_m3: 1.0,
                    linear_rate_per_s: 0.0,
                    quadratic_per_m: 0.0,
                },
            )?;
            let motions = allocated(count, empty_motion())?;
            let workspace = SharedRecoilWorkspace::new(count, &mut cancelled)?;
            for contact in &mut self.contacts[shape.start..shape.start + count] {
                poll(&mut cancelled)?;
                contact.drag = recoil::at_velocity(contact.drag, support.velocity_m_s);
            }
            prepared.push(Some(BranchSupport {
                accepted: support,
                trial: support,
                exchange: None,
                candidate_exchange: SupportExchange::default(),
                workspace,
                inputs,
                motions,
            }));
        }
        poll(&mut cancelled)?;
        self.shared_supports = Some(SharedBranchSupports {
            branches: prepared,
            control,
            momentum_tolerance_ns,
        });
        Ok(self)
    }

    /// Accepted common support in this branch's original axial coordinate.
    pub fn branch_shared_support(&self, branch: usize) -> Option<SlidingSupport> {
        Some(
            self.shared_supports
                .as_ref()?
                .branches
                .get(branch)?
                .as_ref()?
                .accepted,
        )
    }

    /// Complete accepted friction-only exchange with this branch's ONE mass.
    ///
    /// Gas energy and impulse are summed over its cells and all source halves;
    /// support momentum and kinetic energy use the actual complete endpoints.
    /// Original thermal exchange remains in `NetworkWallReport`.
    pub fn accepted_branch_shared_support_exchange(
        &self,
        branch: usize,
    ) -> Option<SupportExchange> {
        self.last_report?;
        self.shared_supports
            .as_ref()?
            .branches
            .get(branch)?
            .as_ref()?
            .exchange
    }

    pub(super) fn prepare_shared_supports(
        &mut self,
        network: &MixtureNetwork1d,
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        if self.shared_supports.is_none() {
            return Ok(());
        }
        let n = self.gas.species_count();
        let width = n + 2;
        for (branch_index, (branch, shape)) in
            network.branches.iter().zip(&self.branches).enumerate()
        {
            let Some(support) = self
                .shared_supports
                .as_ref()
                .and_then(|supports| supports.branches[branch_index].as_ref())
                .map(|support| support.trial)
            else {
                continue;
            };
            poll(cancelled)?;
            if branch.gas != self.gas || branch.mesh != shape.mesh {
                return Err(WallExchangeError::Invalid(
                    "shared support branch geometry changed",
                ));
            }
            let count = shape.mesh.cell_count();
            for (i, values) in branch.states.chunks_exact(width).enumerate() {
                let index = shape.start + i;
                let state = Cell {
                    partial_densities: &values[..n],
                    momentum: values[n],
                    energy: values[n + 1],
                };
                let mut contact = self.candidate_contacts[index];
                contact.drag = recoil::at_velocity(contact.drag, support.velocity_m_s);
                #[cfg(feature = "material")]
                let contact = self.film_contact(index, contact);
                let contact = match self.transport.as_mut() {
                    Some(transport) => transport.resolve(
                        &self.gas,
                        branch_index,
                        state,
                        shape.mesh.volumes()[i],
                        contact,
                        cancelled,
                    )?,
                    None => contact,
                };
                self.candidate_contacts[index] = contact;
            }
            let shared = self
                .shared_supports
                .as_mut()
                .ok_or(WallExchangeError::Invalid(
                    "shared support owner disappeared",
                ))?;
            let prepared =
                shared.branches[branch_index]
                    .as_mut()
                    .ok_or(WallExchangeError::Invalid(
                        "shared support branch disappeared",
                    ))?;
            prepared.prepare(
                &self.gas,
                &branch.states,
                shape.mesh.volumes(),
                &self.candidate_contacts[shape.start..shape.start + count],
                duration,
                shared.control,
                RecoilControl {
                    wall: controls,
                    momentum_tolerance_ns: shared.momentum_tolerance_ns,
                },
                cancelled,
            )?;
        }
        Ok(())
    }
}

impl SharedBranchSupports {
    pub(super) fn has_branch(&self, branch: usize) -> bool {
        self.branches.get(branch).is_some_and(Option::is_some)
    }

    pub(super) fn begin_trial(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        for branch in &mut self.branches {
            poll(cancelled)?;
            if let Some(branch) = branch {
                branch.trial = branch.accepted;
                branch.candidate_exchange = SupportExchange::default();
            }
        }
        Ok(())
    }

    pub(super) fn motion(&self, branch: usize, cell: usize) -> Option<WallUpdate> {
        self.branches
            .get(branch)?
            .as_ref()?
            .motions
            .get(cell)
            .copied()
    }

    pub(super) fn commit(&mut self) {
        for branch in self.branches.iter_mut().flatten() {
            branch.accepted = branch.trial;
            branch.exchange = Some(branch.candidate_exchange);
        }
    }
}

impl BranchSupport {
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        &mut self,
        gas: &MixtureGas,
        states: &[f64],
        volumes: &[f64],
        contacts: &[WallContact],
        duration: Time,
        control: SharedRecoilControl,
        balance: RecoilControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        let n = gas.species_count();
        let width = n + 2;
        for (i, values) in states.chunks_exact(width).enumerate() {
            poll(cancelled)?;
            let state = Cell {
                partial_densities: &values[..n],
                momentum: values[n],
                energy: values[n + 1],
            };
            let primitive = gas.primitive(state, &mut *cancelled)?;
            let law = recoil::coefficients(contacts[i].drag);
            self.inputs[i] = SharedRecoilCell {
                density_kg_m3: primitive.density,
                momentum_kg_m2_s: state.momentum,
                volume_m3: volumes[i],
                linear_rate_per_s: law.linear_rate_per_s,
                quadratic_per_m: law.quadratic_per_m,
            };
        }
        let motion = self.workspace.advance(
            &self.inputs,
            duration,
            self.trial,
            control,
            balance,
            &mut *cancelled,
        )?;
        let mut totals = [0.0; 4];
        let mut corrections = [0.0; 4];
        for (i, (values, cell)) in states.chunks_exact(width).zip(motion.cells).enumerate() {
            poll(cancelled)?;
            let energy = values[n + 1] + cell.work_to_gas_j / volumes[i];
            let gas_energy = (energy - values[n + 1]) * volumes[i];
            let impulse = (cell.momentum - values[n]) * volumes[i];
            let residual = gas_energy - cell.work_to_gas_j;
            if !energy.is_finite()
                || !residual.is_finite()
                || residual.abs()
                    > balance
                        .wall
                        .tolerance(cell.work_to_gas_j.abs().max(cell.slip_dissipation_j))
            {
                return Err(WallExchangeError::Balance {
                    energy_residual_j: residual,
                    thermal_residual_j: 0.0,
                });
            }
            let mut contact = contacts[i];
            contact.drag = recoil::at_velocity(contact.drag, motion.support.velocity_m_s);
            self.motions[i] = WallUpdate {
                momentum: cell.momentum,
                energy,
                contact,
                report: WallCellReport {
                    elapsed: duration,
                    wall_impulse_ns: -impulse,
                    mechanical_work_to_gas_j: cell.work_to_gas_j,
                    drag_heating_j: cell.slip_dissipation_j,
                    energy_residual_j: residual,
                    ..empty_motion().report
                },
            };
            for (j, value) in [
                impulse,
                cell.work_to_gas_j,
                gas_energy,
                cell.slip_dissipation_j,
            ]
            .into_iter()
            .enumerate()
            {
                accumulate(&mut totals[j], &mut corrections[j], value)?;
            }
        }
        let previous = self.candidate_exchange;
        let gas_impulse = sum([previous.gas_impulse_ns, totals[0]]);
        let work = sum([previous.work_to_gas_j, totals[1]]);
        let gas_energy = sum([previous.gas_energy_change_j, totals[2]]);
        let heating = sum([previous.drag_heating_j, totals[3]]);
        let impulse =
            self.accepted.mass_kg * (motion.support.velocity_m_s - self.accepted.velocity_m_s);
        let kinetic =
            impulse * (0.5 * motion.support.velocity_m_s + 0.5 * self.accepted.velocity_m_s);
        let report = SupportExchange {
            gas_impulse_ns: gas_impulse,
            support_impulse_ns: impulse,
            work_to_gas_j: work,
            gas_energy_change_j: gas_energy,
            kinetic_energy_change_j: kinetic,
            drag_heating_j: heating,
            momentum_residual_ns: gas_impulse + impulse,
            energy_residual_j: gas_energy + kinetic,
        };
        if ![
            gas_impulse,
            work,
            gas_energy,
            heating,
            impulse,
            kinetic,
            report.momentum_residual_ns,
            report.energy_residual_j,
        ]
        .into_iter()
        .all(f64::is_finite)
            || heating < 0.0
            || report.momentum_residual_ns.abs() > balance.momentum_tolerance_ns
        {
            return Err(WallExchangeError::Invalid(
                "complete shared support momentum balance",
            ));
        }
        if report.energy_residual_j.abs() > balance.wall.tolerance(work.abs().max(heating)) {
            return Err(WallExchangeError::Balance {
                energy_residual_j: report.energy_residual_j,
                thermal_residual_j: 0.0,
            });
        }
        poll(cancelled)?;
        self.trial = motion.support;
        self.candidate_exchange = report;
        Ok(())
    }
}

fn empty_motion() -> WallUpdate {
    WallUpdate {
        momentum: 0.0,
        energy: 0.0,
        contact: WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        },
        report: WallCellReport {
            elapsed: Time::new(0.0),
            wall_impulse_ns: 0.0,
            mechanical_work_to_gas_j: 0.0,
            drag_heating_j: 0.0,
            heat_to_gas_j: 0.0,
            wall_energy_change_j: 0.0,
            external_heat_to_gas_j: 0.0,
            energy_residual_j: 0.0,
            thermal_residual_j: 0.0,
            iterations: 0,
        },
    }
}
