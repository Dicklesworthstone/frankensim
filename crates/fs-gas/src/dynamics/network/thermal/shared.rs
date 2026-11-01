//! One original finite wall inventory shared by several gas control volumes.

use super::*;
use crate::multicomponent::wall_exchange::solve_scalar;
use std::cell::Cell;

mod conduction;
use conduction::WallLinks;
pub use conduction::{SharedWallLink, SharedWallLinkReport};

/// Original gas control volume attached to a shared, spatially uniform wall.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SharedWallGas {
    /// A retained compact volume, in network order.
    Plenum(usize),
    /// An original resolved branch cell.
    Branch { branch: usize, cell: usize },
}

/// One heat-transfer interface. This contains no second wall heat inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedWallPort {
    /// Original gas amount and volume.
    pub gas: SharedWallGas,
    /// Complete, nonnegative film/interface conductance G [W/K].
    pub conductance_w_per_k: f64,
}

/// One finite solid amount and all gas interfaces sharing its temperature.
/// Storage must be Lumped or original Caloric, with its own conductance zero:
/// all actual interface conductances belong to ports, not to another film.
/// Gas amounts are separate; sharing heat never mixes mass or species.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedWall<'law> {
    /// Original finite storage, counted exactly once.
    pub storage: PlenumThermalContact<'law>,
    /// Distinct original gas cells. A gas cell may belong to only one group.
    pub ports: Vec<SharedWallPort>,
}

/// Total nested work in one complete shared-source attempt. Each original gas
/// primitive/caloric query and wall-state query spends one evaluation; original
/// source inversions retain their own bounded controls. Network retries remain
/// bounded by StepControl.maximum_halvings. No numerical work is unbounded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SharedWallControl {
    /// Across every group, including brackets, failed iterates and final checks.
    pub maximum_source_evaluations: usize,
}
impl Default for SharedWallControl {
    fn default() -> Self {
        Self {
            maximum_source_evaluations: 65_536,
        }
    }
}

/// Accepted gas-contact equations and actual nested work, separate from the
/// complete network's original heat, gravity and transport accounts. With
/// wall links, these contacts precede the ordered link stages; they do not
/// claim a simultaneous endpoint for the entire thermal graph.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedWallReport {
    /// Actual accepted source interval.
    pub elapsed: Time,
    /// Sum of represented shared-wall plus gas storage residuals [J].
    pub energy_residual_j: f64,
    /// Largest absolute gas-group or wall-link storage balance [J].
    pub maximum_wall_residual_j: f64,
    /// Largest absolute endpoint heat-equation residual among gas and wall
    /// interfaces, each at its own accepted composition stage [J].
    pub maximum_contact_residual_j: f64,
    /// Actual calls to the original caloric/state owners.
    pub source_evaluations: usize,
    /// Total safeguarded outer iterations, excluding bracket queries.
    pub wall_iterations: usize,
    /// Number of original implicit gas-interface solves, including trial solves.
    pub contact_solves: usize,
}

struct Work {
    used: Cell<usize>,
    contacts: Cell<usize>,
    limit: usize,
}
impl Work {
    fn charge(&self) -> Result<(), WallExchangeError> {
        let next = self
            .used
            .get()
            .checked_add(1)
            .filter(|&next| next <= self.limit)
            .ok_or(WallExchangeError::Invalid(
                "shared-wall nested source work budget",
            ))?;
        self.used.set(next);
        Ok(())
    }
}

struct MeteredGas<'a> {
    original: FrozenCaloric<'a>,
    work: &'a Work,
}
impl WallGasCaloric for MeteredGas<'_> {
    fn domain(&self) -> CaloricDomain {
        self.original.domain()
    }
    fn primitive(
        &self,
        momentum: f64,
        energy: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallGasState, WallExchangeError> {
        self.work.charge()?;
        self.original.primitive(momentum, energy, cancelled)
    }
    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        self.work.charge()?;
        self.original.caloric(temperature)
    }
}

struct Group<'law> {
    accepted: SharedWall<'law>,
    candidate: WallThermal,
    volumes: Vec<f64>,
    initial: Vec<Conserved>,
    energies: Vec<f64>,
    heat: Vec<f64>,
    residuals: Vec<f64>,
}

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
            groups.push(Group {
                candidate: wall.storage.thermal,
                accepted: wall,
                volumes,
                initial: allocated(
                    n,
                    Conserved {
                        density: 1.0,
                        momentum: 0.0,
                        energy: 0.0,
                    },
                )?,
                energies: allocated(n, 0.0)?,
                heat: allocated(n, 0.0)?,
                residuals: allocated(n, 0.0)?,
            });
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

struct WallPoint {
    thermal: WallThermal,
    temperature: f64,
    change: f64,
    slope: Option<f64>,
}

impl Group<'_> {
    fn after(&self, change: f64, work: &Work) -> Result<WallPoint, WallExchangeError> {
        self.after_state(self.accepted.storage.thermal, change, work)
    }

    // The same original storage, now about a candidate composition-stage state.
    // Neither gas contact nor wall links store a second thermal inventory.
    fn after_state(
        &self,
        thermal: WallThermal,
        change: f64,
        work: &Work,
    ) -> Result<WallPoint, WallExchangeError> {
        work.charge()?;
        let storage = &self.accepted.storage;
        if let Some(source) = &storage.heat_storage {
            let point = source.after(
                WallContact {
                    drag: WallDrag::None,
                    thermal,
                },
                change,
            )?;
            return Ok(WallPoint {
                thermal: point.contact.thermal,
                temperature: point.temperature_k,
                change: point.energy_change_j,
                slope: point.temperature_slope_k_per_j,
            });
        }
        let WallThermal::Lumped {
            temperature_k,
            heat_capacity_j_per_k,
            ..
        } = thermal
        else {
            return Err(WallExchangeError::Invalid(
                "shared wall missing caloric source",
            ));
        };
        let temperature = temperature_k + change / heat_capacity_j_per_k;
        let actual = heat_capacity_j_per_k * (temperature - temperature_k);
        if !positive(temperature) || !actual.is_finite() {
            return Err(WallExchangeError::Invalid(
                "shared wall represented storage",
            ));
        }
        Ok(WallPoint {
            thermal: WallThermal::Lumped {
                temperature_k: temperature,
                heat_capacity_j_per_k,
                conductance_w_per_k: 0.0,
            },
            temperature,
            change: actual,
            slope: Some(1.0 / heat_capacity_j_per_k),
        })
    }

    // Restrict the ORIGINAL storage's energy domain to the passive initial
    // temperature hull. Bisection never approximates or replaces its calorics.
    fn bound(
        &self,
        bound: f64,
        temperature: f64,
        lower: bool,
        control: WallControl,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, WallExchangeError> {
        let inside = |t| {
            if lower {
                t >= temperature
            } else {
                t <= temperature
            }
        };
        if inside(self.after(bound, work)?.temperature) {
            return Ok(bound);
        }
        let (mut outer, mut inner) = (bound, 0.0);
        for _ in 0..128 {
            poll(cancelled)?;
            let mid = 0.5 * outer + 0.5 * inner;
            if mid == outer
                || mid == inner
                || (outer - inner).abs() <= 0.0625 * control.tolerance(inner.abs())
            {
                return Ok(inner);
            }
            if inside(self.after(mid, work)?.temperature) {
                inner = mid;
            } else {
                outer = mid;
            }
        }
        Err(WallExchangeError::TemperatureIterations)
    }

    fn evaluate(
        &mut self,
        gas: &FrozenGas,
        change: f64,
        duration: Time,
        control: WallControl,
        work: &Work,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(f64, Option<f64>, f64), WallExchangeError> {
        let wall = self.after(change, work)?;
        let mut sums = [0.0; 3];
        let mut corrections = [0.0; 3];
        for (i, port) in self.accepted.ports.iter().enumerate() {
            if cancelled() {
                return Err(WallExchangeError::Cancelled);
            }
            let state = self.initial[i];
            let caloric = MeteredGas {
                original: FrozenCaloric {
                    gas,
                    density: state.density,
                    hint: Some(wall.temperature),
                },
                work,
            };
            work.contacts.set(work.contacts.get().checked_add(1).ok_or(
                WallExchangeError::Invalid("shared wall contact count overflow"),
            )?);
            let update = exchange_caloric(
                &caloric,
                state.momentum,
                state.energy,
                self.volumes[i],
                duration,
                WallContact {
                    drag: WallDrag::None,
                    thermal: WallThermal::Reservoir {
                        temperature_k: wall.temperature,
                        conductance_w_per_k: port.conductance_w_per_k,
                    },
                },
                control,
                &mut *cancelled,
            )?;
            self.energies[i] = update.energy;
            self.heat[i] = update.report.heat_to_gas_j;
            self.residuals[i] = update.report.thermal_residual_j;
            let tangent = if wall.slope.is_some() {
                let endpoint = caloric.primitive(state.momentum, update.energy, cancelled)?;
                let capacity = caloric.caloric(endpoint.temperature)?.1 * self.volumes[i];
                let integrated = port.conductance_w_per_k * duration.value();
                if !positive(capacity) {
                    return Err(WallExchangeError::Invalid("shared gas capacity"));
                }
                if integrated <= capacity {
                    integrated / (1.0 + integrated / capacity)
                } else {
                    capacity / (1.0 + capacity / integrated)
                }
            } else {
                0.0
            };
            accumulate(
                &mut sums,
                &mut corrections,
                [self.heat[i], self.heat[i].abs(), tangent],
            )?;
        }
        let residual = wall.change + sums[0];
        let slope = wall.slope.map(|s| 1.0 + s * sums[2]);
        if !residual.is_finite() || slope.is_some_and(|s| !positive(s)) {
            return Err(WallExchangeError::Invalid("shared wall coupled equation"));
        }
        Ok((residual, slope, wall.change.abs().max(sums[1])))
    }

    fn solve(
        &mut self,
        gas: &FrozenGas,
        duration: Time,
        control: WallControl,
        work: &Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, WallExchangeError> {
        let wall_initial = self.after(0.0, work)?.temperature;
        let initial = self.evaluate(gas, 0.0, duration, control, work, cancelled)?;
        let mut iterations = 0;
        let mut accepted_change = 0.0;
        if initial.0.abs() > control.tolerance(initial.2) {
            let mut minimum = wall_initial;
            let mut maximum = wall_initial;
            for state in &self.initial {
                let original = MeteredGas {
                    original: FrozenCaloric {
                        gas,
                        density: state.density,
                        hint: Some(wall_initial),
                    },
                    work,
                };
                let t = original
                    .primitive(state.momentum, state.energy, cancelled)?
                    .temperature;
                minimum = minimum.min(t);
                maximum = maximum.max(t);
            }
            let (mut low, mut high) = if let Some(source) = &self.accepted.storage.heat_storage {
                let contact = WallContact {
                    drag: WallDrag::None,
                    thermal: self.accepted.storage.thermal,
                };
                work.charge()?;
                let cold = source.heat_limits(contact, minimum)?;
                work.charge()?;
                let hot = source.heat_limits(contact, maximum)?;
                (-cold[1].max(hot[1]), -cold[0].min(hot[0]))
            } else {
                let WallThermal::Lumped {
                    heat_capacity_j_per_k,
                    ..
                } = self.accepted.storage.thermal
                else {
                    unreachable!("finite storage admitted");
                };
                (
                    heat_capacity_j_per_k * (minimum - wall_initial),
                    heat_capacity_j_per_k * (maximum - wall_initial),
                )
            };
            if low < 0.0 {
                low = low.next_up();
            }
            if high > 0.0 {
                high = high.next_down();
            }
            low = self.bound(low, minimum, true, control, work, cancelled)?;
            high = self.bound(high, maximum, false, control, work, cancelled)?;
            if !low.is_finite() || !high.is_finite() || low >= high {
                return Err(WallExchangeError::Invalid(
                    "shared wall supported heat interval",
                ));
            }
            let result = solve_scalar(
                low,
                high,
                0.0,
                control,
                |change, checkpoint| {
                    self.evaluate(gas, change, duration, control, work, checkpoint)
                },
                cancelled,
            )?;
            iterations = result.1;
            self.evaluate(gas, result.0, duration, control, work, cancelled)?;
            accepted_change = result.0;
        }
        let mut sums = [0.0; 2];
        let mut corrections = [0.0; 2];
        for heat in &self.heat {
            accumulate(&mut sums, &mut corrections, [*heat, heat.abs()])?;
        }
        // Retain the ONE endpoint inventory whose ORIGINAL represented change
        // closed the coupled equation. Reprojecting it after the root would
        // change every interface flux, especially at very stiff conductance.
        let wall = self.after(accepted_change, work)?;
        let residual = sums[0] + wall.change;
        let mut maximum_flux = 0.0_f64;
        for i in 0..self.accepted.ports.len() {
            poll(cancelled)?;
            let flux_residual = self.residuals[i];
            if !flux_residual.is_finite()
                || flux_residual.abs() > control.tolerance(self.heat[i].abs())
            {
                return Err(WallExchangeError::Balance {
                    energy_residual_j: residual,
                    thermal_residual_j: flux_residual,
                });
            }
            maximum_flux = maximum_flux.max(flux_residual.abs());
        }
        if !residual.is_finite()
            || residual.abs() > control.tolerance(sums[1].max(wall.change.abs()))
        {
            return Err(WallExchangeError::Balance {
                energy_residual_j: residual,
                thermal_residual_j: maximum_flux,
            });
        }
        self.candidate = wall.thermal;
        Ok(WallCellReport {
            elapsed: duration,
            wall_impulse_ns: 0.0,
            drag_heating_j: 0.0,
            heat_to_gas_j: sums[0],
            wall_energy_change_j: wall.change,
            external_heat_to_gas_j: 0.0,
            energy_residual_j: residual,
            thermal_residual_j: maximum_flux,
            iterations,
        })
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
                group.initial[i] = *gas_state(port.gas, branches, plenums);
            }
            let report = group.solve(gas, duration, control, &work, cancelled)?;
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
