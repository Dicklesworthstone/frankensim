//! One finite thermal inventory on the two faces of an inertial partition.

use super::*;
use crate::multicomponent::wall_exchange::shared::{
    Group, SharedWall, SharedWallGas, SharedWallPort, Work, mixture,
};
use crate::multicomponent::wall_exchange::{
    SharedWallControl, SharedWallStorage, WallCaloricConstraint, WallCellReport, WallContact,
    WallControl, WallDrag, WallExchangeError, WallThermal,
};
use std::cell::Cell as Counter;

/// Actual heat exchanged with one original, finite piston inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedPistonHeatReport {
    /// Before and after the mechanical step, each over half its duration.
    /// Both reports refer to the SAME wall, with the second starting at the
    /// first endpoint. No thermal energy or mechanical work is duplicated.
    pub stages: [WallCellReport; 2],
    /// Heat into [left, right] gas across both half steps [J].
    pub heat_to_gas_j: [f64; 2],
    /// Original wall's represented final minus initial stored energy [J].
    pub wall_energy_change_j: f64,
    /// Sum of both gas heats plus the one wall-storage change [J].
    pub energy_residual_j: f64,
    /// Original caloric/state calls across both accepted heat stages.
    pub source_evaluations: usize,
    /// Implicit gas-contact solves, including the wall root's trial endpoints.
    pub contact_solves: usize,
}

/// Full thermal/mechanical state change with the original transport reports.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedPistonThermalStepReport<R = StepReport> {
    /// Actual pressure/work stages, between the two heat stages. These native
    /// before/after fields deliberately exclude the surrounding heat exchange.
    pub transport: SharedPistonGasReport<R>,
    /// Both chambers before any heat or pressure work.
    pub before: SharedPistonGasReport<Totals>,
    /// Both chambers after the complete heat/pressure/heat composition.
    pub after: SharedPistonGasReport<Totals>,
    /// Separate physical heat and finite-wall storage accounts.
    pub heat: SharedPistonHeatReport,
    /// Complete thermal/mechanical restarts, including acoustic shortening.
    pub restarts: usize,
}

/// A refused coupled increment leaves both gases, piston and wall unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum SharedPistonThermalError {
    /// Original pressure, geometry, source or inertial balance refusal.
    Piston(PistonError),
    /// Original shared caloric/heat solve, retaining its concrete cause.
    Heat(WallExchangeError),
    /// Cancellation before coupled publication.
    Cancelled,
    /// Every bounded complete composition was refused.
    StepRetries {
        /// Initial attempt plus complete restarts.
        attempts: usize,
        /// Last concrete heat or pressure refusal.
        last: Box<SharedPistonThermalError>,
    },
}
impl core::fmt::Display for SharedPistonThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shared thermal piston: {self:?}")
    }
}
impl std::error::Error for SharedPistonThermalError {}
impl From<PistonError> for SharedPistonThermalError {
    fn from(error: PistonError) -> Self {
        match error {
            PistonError::Cancelled => Self::Cancelled,
            other => Self::Piston(other),
        }
    }
}
impl From<Error> for SharedPistonThermalError {
    fn from(error: Error) -> Self {
        PistonError::from(error).into()
    }
}
impl From<WallExchangeError> for SharedPistonThermalError {
    fn from(error: WallExchangeError) -> Self {
        match error {
            WallExchangeError::Cancelled => Self::Cancelled,
            other => Self::Heat(other),
        }
    }
}

/// One spatially uniform, finite thermal store shared by both piston faces.
///
/// Each gas retains its original calorics, species axis and end-cell volume.
/// Conductances are complete authored interface G=hA [W/K], not per-area h;
/// no second multiplication by the mechanical cross section is performed.
/// Thermal mass belongs to this inventory once, and does not add to the
/// separately authored total piston inertia. Translation adds no wall p dV.
/// Source enthalpy is refused unless its existing caloric binding explicitly
/// declares negligible solid expansion; rigid translation cannot pay hidden
/// pressure/deformation work.
///
/// Two complete private flow workspaces are allocated once. Every duration
/// retry resets them from accepted fields; velocity proposals reuse the
/// original shared-piston scratch. The initial heat solve is independent of
/// midpoint velocity and runs once per duration, not once per root proposal.
/// The caloric work budget covers both heat halves of ONE complete attempt;
/// each restart resets it, under `StepControl::maximum_halvings`. Accepted
/// source-work counts describe the accepted attempt, not discarded retries.
pub struct SharedPistonThermal<'law> {
    group: Group<'law>,
    trial: [MixtureEuler1d; 2],
    work: Work,
    last_report: Option<SharedPistonHeatReport>,
}

impl<'law> SharedPistonThermal<'law> {
    /// Bind one finite wall to the last left cell and first right cell.
    /// The storage's own conductance must be zero; the two supplied values
    /// represent its distinct faces. A prescribed reservoir is not a finite
    /// partition and is rejected. Wall/source temperatures must fit both gas
    /// caloric paths. The fixed thermal amount does not follow gas density.
    pub fn new(
        left: &MixtureEuler1d,
        right: &MixtureEuler1d,
        mut storage: SharedWallStorage<'law>,
        conductance_w_per_k: [f64; 2],
        controls: SharedWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SharedPistonThermalError> {
        poll(&mut cancelled)?;
        admit_geometry(left, right)?;
        let mut contact = WallContact {
            drag: WallDrag::None,
            thermal: storage.thermal,
        };
        contact.validate()?;
        let temperature = if let Some(source) = &storage.heat_storage {
            if source.constraint() == WallCaloricConstraint::SourceEnthalpy {
                return Err(WallExchangeError::Invalid(
                    "translating shared piston requires internal energy or explicit negligible expansion",
                ).into());
            }
            contact = source.bind_contact(contact)?;
            source.temperature_k(contact.thermal)?
        } else if let WallThermal::Lumped { temperature_k, .. } = contact.thermal {
            temperature_k
        } else {
            return Err(WallExchangeError::Invalid(
                "shared thermal piston requires original finite wall storage",
            )
            .into());
        };
        let wall_conductance = match contact.thermal {
            WallThermal::Lumped {
                conductance_w_per_k,
                ..
            }
            | WallThermal::Caloric {
                conductance_w_per_k,
                ..
            } => conductance_w_per_k,
            _ => {
                return Err(
                    WallExchangeError::Invalid("shared piston needs one finite inventory").into(),
                );
            }
        };
        if wall_conductance != 0.0
            || controls.maximum_source_evaluations == 0
            || conductance_w_per_k
                .iter()
                .any(|g| !g.is_finite() || *g < 0.0)
            || [left, right].iter().any(|flow| {
                let domain = flow.gas.domain();
                !(domain.minimum_temperature..=domain.maximum_temperature).contains(&temperature)
            })
        {
            return Err(WallExchangeError::Invalid(
                "shared piston heat conductances, work budget or source interval",
            )
            .into());
        }
        storage.thermal = contact.thermal;
        let cells = [left.cell_count() - 1, 0];
        let mut ports = crate::multicomponent::allocated(
            2,
            SharedWallPort {
                gas: SharedWallGas::Branch {
                    branch: 0,
                    cell: cells[0],
                },
                conductance_w_per_k: conductance_w_per_k[0],
            },
        )?;
        ports[1] = SharedWallPort {
            gas: SharedWallGas::Branch {
                branch: 1,
                cell: cells[1],
            },
            conductance_w_per_k: conductance_w_per_k[1],
        };
        let mut volumes = crate::multicomponent::allocated(2, left.mesh.volumes()[cells[0]])?;
        volumes[1] = right.mesh.volumes()[cells[1]];
        let result = Self {
            group: Group::new(SharedWall { storage, ports }, volumes)?,
            trial: [
                left.duplicate(&mut cancelled)?,
                right.duplicate(&mut cancelled)?,
            ],
            work: Work {
                used: Counter::new(0),
                contacts: Counter::new(0),
                limit: controls.maximum_source_evaluations,
            },
            last_report: None,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Accepted original storage; no rejected or half-step temperature leaks.
    pub fn storage(&self) -> &SharedWallStorage<'law> {
        &self.group.accepted.storage
    }

    /// Authored complete conductances of [left, right] piston faces [W/K].
    pub fn conductances_w_per_k(&self) -> [f64; 2] {
        [
            self.group.accepted.ports[0].conductance_w_per_k,
            self.group.accepted.ports[1].conductance_w_per_k,
        ]
    }

    /// Heat account from the last complete accepted gas/wall/piston step.
    pub const fn last_report(&self) -> Option<SharedPistonHeatReport> {
        self.last_report
    }

    fn heat_stage(
        &mut self,
        duration: Time,
        control: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(WallCellReport, [f64; 2]), SharedPistonThermalError> {
        for (side, flow) in self.trial.iter().enumerate() {
            let cell = if side == 0 { flow.cell_count() - 1 } else { 0 };
            let state = flow.cell(cell).expect("bound piston end cell");
            self.group.initial[side] = [state.momentum, state.energy];
            self.group.volumes[side] = flow.mesh.volumes()[cell];
        }
        let fields = &self.trial;
        let report = self.group.solve(
            &|port, _wall_temperature| {
                let SharedWallGas::Branch { branch, cell } = port else {
                    unreachable!("only the two original piston-end ports are bound")
                };
                let flow = &fields[branch];
                mixture(
                    &flow.gas,
                    flow.cell(cell).expect("bound end cell").partial_densities,
                )
            },
            duration,
            control,
            &self.work,
            self.group.candidate,
            cancelled,
        )?;
        for (side, flow) in self.trial.iter_mut().enumerate() {
            let cell = if side == 0 { flow.cell_count() - 1 } else { 0 };
            let stride = flow.gas.species_count() + 2;
            flow.states[cell * stride + stride - 1] = self.group.energies[side];
        }
        Ok((report, [self.group.heat[0], self.group.heat[1]]))
    }

    fn complete_heat(
        &self,
        stages: [WallCellReport; 2],
        heat: [[f64; 2]; 2],
        controls: WallControl,
    ) -> Result<SharedPistonHeatReport, SharedPistonThermalError> {
        let initial = self.group.accepted.storage.thermal;
        let final_state = self.group.candidate;
        let wall_change = if let Some(source) = &self.group.accepted.storage.heat_storage {
            source.energy_change_j(
                WallContact {
                    drag: WallDrag::None,
                    thermal: initial,
                },
                WallContact {
                    drag: WallDrag::None,
                    thermal: final_state,
                },
            )?
        } else {
            let (
                WallThermal::Lumped {
                    temperature_k: before,
                    heat_capacity_j_per_k,
                    ..
                },
                WallThermal::Lumped {
                    temperature_k: after,
                    ..
                },
            ) = (initial, final_state)
            else {
                unreachable!("one admitted original finite inventory")
            };
            heat_capacity_j_per_k * (after - before)
        };
        let heat_to_gas_j = [sum([heat[0][0], heat[1][0]]), sum([heat[0][1], heat[1][1]])];
        let residual = sum([heat_to_gas_j[0], heat_to_gas_j[1], wall_change]);
        let scale = sum(heat.into_iter().flatten().map(f64::abs)).max(wall_change.abs());
        if !residual.is_finite() || residual.abs() > controls.tolerance(scale) {
            return Err(WallExchangeError::Balance {
                energy_residual_j: residual,
                thermal_residual_j: stages[0]
                    .thermal_residual_j
                    .abs()
                    .max(stages[1].thermal_residual_j.abs()),
            }
            .into());
        }
        Ok(SharedPistonHeatReport {
            stages,
            heat_to_gas_j,
            wall_energy_change_j: wall_change,
            energy_residual_j: residual,
            source_evaluations: self.work.used.get(),
            contact_solves: self.work.contacts.get(),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn advance(
        &mut self,
        left: &mut MixtureEuler1d,
        right: &mut MixtureEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        piston_controls: PistonControl,
        mut equilibrium: Option<([ContactEquilibriumControl; 2], [&mut ContactWorkspace; 2])>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
        poll(&mut cancelled)?;
        admit_geometry(left, right)?;
        validate(
            *piston,
            maximum_duration,
            gas_controls,
            PistonLoad {
                exterior_pressure_pa: 0.0,
                external_force_n,
            },
            piston_controls,
        )?;
        heat_controls.validate()?;
        if let Some((controls, workspaces)) = equilibrium.as_ref() {
            workspaces[0].admit(left.cell_count(), controls[0])?;
            workspaces[1].admit(right.cell_count(), controls[1])?;
        }
        // Both fields participate in the final workspace swap. Admit their
        // original geometry once, before any motion: reconstructing scratch
        // from rounded moved face differences can lose the original metrics.
        // As in native pressure stepping, a refusal may retain allocated scratch
        // but never publishes physical fields, geometry, clock or reports.
        for flow in [&mut *left, &mut *right] {
            if flow.moving.is_none() {
                flow.moving = Some(MovingWorkspace::new(&flow.mesh).map_err(geometry_error)?);
            }
        }
        self.trial[0].copy_physical_with_motion_from(left, &mut cancelled)?;
        self.trial[1].copy_physical_with_motion_from(right, &mut cancelled)?;
        let bracket = piston_controls.midpoint_velocity_bracket_m_s;
        let (left_before, left_rate) = prepare(&mut self.trial[0], bracket, &mut cancelled)?;
        let (right_before, right_rate) = prepare(&mut self.trial[1], bracket, &mut cancelled)?;
        let before = SharedPistonGasReport {
            left: left_before,
            right: right_before,
        };
        let mut dt = maximum_duration
            .value()
            .min(0.9 * gas_controls.cfl / left_rate.max(right_rate));
        let interface = left.mesh.positions()[left.cell_count()];
        if bracket[0] < 0.0 {
            dt = dt.min(
                0.9 * gas_controls.cfl * ((interface - left.mesh.positions()[0]) / -bracket[0]),
            );
        }
        if bracket[1] > 0.0 {
            dt = dt.min(
                0.9 * gas_controls.cfl
                    * ((right.mesh.positions()[right.cell_count()] - interface) / bracket[1]),
            );
        }
        let original_piston = *piston;
        let mut last = PistonError::Invalid("no complete thermal shared-piston increment").into();
        for restarts in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt) || !positive(0.5 * dt) || !positive(left.time + dt - left.time) {
                return Err(PistonError::Invalid(
                    "thermal shared-piston duration cannot advance time",
                )
                .into());
            }
            self.trial[0].copy_physical_with_motion_from(left, &mut cancelled)?;
            self.trial[1].copy_physical_with_motion_from(right, &mut cancelled)?;
            self.group.candidate = self.group.accepted.storage.thermal;
            self.work.used.set(0);
            self.work.contacts.set(0);
            let mut next_dt = 0.5 * dt;
            let attempt = (|| {
                let (first, first_heat) =
                    self.heat_stage(Time::new(0.5 * dt), heat_controls, &mut cancelled)?;
                let mut trial_piston = original_piston;
                let [trial_left, trial_right] = &mut self.trial;
                let mut step = trial_left
                    .advance_shared_piston(
                        trial_right,
                        &mut trial_piston,
                        Time::new(dt),
                        gas_controls,
                        external_force_n,
                        None,
                        piston_controls,
                        equilibrium
                            .as_mut()
                            .map(|(controls, [a, b])| (*controls, [&mut **a, &mut **b])),
                        &mut cancelled,
                    )?
                    .piston;
                let accepted_dt = step.gas.left.elapsed.value();
                if accepted_dt != dt {
                    next_dt = accepted_dt;
                    return Ok(None);
                }
                let (second, second_heat) =
                    self.heat_stage(Time::new(0.5 * dt), heat_controls, &mut cancelled)?;
                let heat =
                    self.complete_heat([first, second], [first_heat, second_heat], heat_controls)?;
                let after = SharedPistonGasReport {
                    left: totals(
                        &self.trial[0].states,
                        self.trial[0].mesh.volumes(),
                        self.trial[0].gas.species_count(),
                        &mut cancelled,
                    )?,
                    right: totals(
                        &self.trial[1].states,
                        self.trial[1].mesh.volumes(),
                        self.trial[1].gas.species_count(),
                        &mut cancelled,
                    )?,
                };
                let defect = equilibrium.as_ref().map_or(0.0, |(_, workspaces)| {
                    sum([
                        workspaces[0].receipt().signed_energy_j,
                        workspaces[1].receipt().signed_energy_j,
                    ])
                });
                step.total_energy_residual_j = sum([
                    after.left.energy - before.left.energy,
                    after.right.energy - before.right.energy,
                    heat.wall_energy_change_j,
                    step.kinetic_energy_change_j,
                    -step.external_force_work_j,
                    -defect,
                ]);
                if !step.total_energy_residual_j.is_finite() || !step.admitted(piston_controls) {
                    return Err(PistonError::Balance {
                        momentum_residual_n_s: step
                            .momentum_residual_n_s
                            .abs()
                            .max(step.combined_momentum_residual_n_s.abs()),
                        energy_residual_j: step
                            .mechanical_energy_residual_j
                            .abs()
                            .max(step.total_energy_residual_j.abs()),
                    }
                    .into());
                }
                Ok(Some(step.map_gas(|transport| {
                    SharedPistonThermalStepReport {
                        transport,
                        before,
                        after,
                        heat,
                        restarts,
                    }
                })))
            })();
            match attempt {
                Ok(Some(step)) => {
                    poll(&mut cancelled)?;
                    core::mem::swap(left, &mut self.trial[0]);
                    core::mem::swap(right, &mut self.trial[1]);
                    self.group.accepted.storage.thermal = self.group.candidate;
                    self.last_report = Some(step.gas.heat);
                    piston.accept(&step);
                    return Ok(step);
                }
                Ok(None) => {
                    last = PistonError::Invalid(
                        "heated shared-piston acoustics shortened the complete interval",
                    )
                    .into();
                }
                Err(SharedPistonThermalError::Cancelled) => {
                    return Err(SharedPistonThermalError::Cancelled);
                }
                Err(
                    error @ SharedPistonThermalError::Piston(PistonError::MixtureGas(
                        Error::Allocation,
                    )),
                )
                | Err(
                    error @ SharedPistonThermalError::Heat(WallExchangeError::Flow(
                        Error::Allocation,
                    )),
                ) => return Err(error),
                Err(error) => last = error,
            }
            dt = next_dt;
        }
        Err(SharedPistonThermalError::StepRetries {
            attempts: gas_controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}

impl MixtureEuler1d {
    /// Couple resolved pressures and piston inertia to heat through ONE finite wall.
    ///
    /// The stages are H(dt/2), A(dt), H(dt/2): both gas faces share the same
    /// implicit wall endpoint in each H, while A is the existing two-field
    /// pressure/work solve. The first H uses original cell volumes; the second
    /// uses actual moved volumes and current compositions. Its two backward-
    /// Euler heat stages do not establish second order or fully simultaneous
    /// thermomechanical coupling. A shortened A restarts the entire composition.
    ///
    /// Accepted gas energy, one finite-wall storage change, piston kinetic
    /// change and applied-force work close before publication. Heat never enters
    /// the piston pressure-work arrays, and both outer gas walls remain adiabatic.
    /// No chemistry, leakage, molecular transport, expansion or through-wall
    /// temperature gradient is inferred. Repeated steps reuse all trial buffers.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston_and_heat(
        &mut self,
        right: &mut MixtureEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
        thermal.advance(
            self,
            right,
            piston,
            maximum_duration,
            gas_controls,
            heat_controls,
            external_force_n,
            piston_controls,
            None,
            cancelled,
        )
    }

    /// Retain independently bounded material contacts during shared-wall heating.
    ///
    /// Contact targets are taken AFTER the first physical heat stage. Resolved
    /// heating gradients and actual two-sided pressure differences remain in the
    /// original mechanical solve. The full physical balance subtracts only the
    /// separately reported signed numerical defects; neither gas heat nor wall
    /// energy is relabeled as a contact correction. Both chambers independently
    /// enforce their existing absolute defect budget.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston_heat_and_contact_equilibrium(
        &mut self,
        right: &mut MixtureEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        piston_controls: PistonControl,
        equilibrium: [ContactEquilibriumControl; 2],
        workspaces: [&mut ContactWorkspace; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        PistonStep<SharedPistonThermalStepReport<ContactEquilibriumStepReport<StepReport>>>,
        SharedPistonThermalError,
    > {
        let [a, b] = workspaces;
        let step = thermal.advance(
            self,
            right,
            piston,
            maximum_duration,
            gas_controls,
            heat_controls,
            external_force_n,
            piston_controls,
            Some((equilibrium, [&mut *a, &mut *b])),
            cancelled,
        )?;
        Ok(step.map_gas(|gas| SharedPistonThermalStepReport {
            transport: contact_reports(gas.transport, [a.receipt(), b.receipt()]),
            before: gas.before,
            after: gas.after,
            heat: gas.heat,
            restarts: gas.restarts,
        }))
    }
}

fn admit_geometry(left: &MixtureEuler1d, right: &MixtureEuler1d) -> Result<(), PistonError> {
    let closed = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    if left.boundaries != closed
        || right.boundaries != closed
        || left.mesh.uniform_cross_section().is_none()
        || left.mesh.uniform_cross_section() != right.mesh.uniform_cross_section()
        || left.mesh.positions()[left.cell_count()] != right.mesh.positions()[0]
        || left.time != right.time
    {
        return Err(PistonError::Invalid(
            "shared thermal piston requires closed chambers with equal area, interface and clocks",
        ));
    }
    Ok(())
}
