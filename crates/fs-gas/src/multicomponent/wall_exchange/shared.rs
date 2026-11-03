//! Common implicit energy balance for one original wall and several gas amounts.
//!
//! Frozen and evolving-mixture networks share this scalar owner. Each caller
//! supplies original calorics at its current composition and owns publication.

use super::*;
use crate::dynamics::CaloricDomain;
use std::cell::Cell;

pub(crate) mod links;

/// One compact volume's thermal boundary and optional original finite storage.
/// Distinct contacts represent distinct physical wall masses. A shared wall
/// requires its own thermal owner, rather than duplicating its stored energy.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedWallStorage<'law> {
    /// Prescribed reservoir, authored constant capacity, or original h/u state.
    pub thermal: WallThermal,
    /// Original material calorics, required for `WallThermal::Caloric` and
    /// optional for a temperature-coordinate finite wall. None retains the
    /// explicitly authored constant capacity; no Cp-to-Cv conversion is inferred.
    pub heat_storage: Option<MaterialWallCaloric<'law>>,
}

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
    pub storage: SharedWallStorage<'law>,
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

pub(crate) struct Work {
    pub(crate) used: Cell<usize>,
    pub(crate) contacts: Cell<usize>,
    pub(crate) limit: usize,
}
impl Work {
    pub(crate) fn charge(&self) -> Result<(), WallExchangeError> {
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

struct MeteredGas<'a, G> {
    original: G,
    work: &'a Work,
}
impl<G: WallGasCaloric> WallGasCaloric for MeteredGas<'_, G> {
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

pub(crate) struct Group<'law> {
    pub(crate) accepted: SharedWall<'law>,
    pub(crate) candidate: WallThermal,
    initial_wall: WallThermal,
    pub(crate) volumes: Vec<f64>,
    pub(crate) initial: Vec<[f64; 2]>,
    pub(crate) energies: Vec<f64>,
    pub(crate) heat: Vec<f64>,
    residuals: Vec<f64>,
}

impl<'law> Group<'law> {
    pub(crate) fn new(
        wall: SharedWall<'law>,
        volumes: Vec<f64>,
    ) -> Result<Self, WallExchangeError> {
        let n = wall.ports.len();
        Ok(Self {
            candidate: wall.storage.thermal,
            initial_wall: wall.storage.thermal,
            accepted: wall,
            volumes,
            initial: super::super::allocated(n, [0.0; 2])?,
            energies: super::super::allocated(n, 0.0)?,
            heat: super::super::allocated(n, 0.0)?,
            residuals: super::super::allocated(n, 0.0)?,
        })
    }
}

// Each source interval borrows the CURRENT constituent densities. No frozen
// mixture, copied caloric table, or transformed formation reference is created.
pub(crate) fn mixture<'a>(
    gas: &'a MixtureGas,
    partial_densities: &'a [f64],
) -> impl WallGasCaloric + 'a {
    MixtureCaloric {
        gas,
        partial_densities,
    }
}

pub(crate) struct WallPoint {
    pub(crate) thermal: WallThermal,
    pub(crate) temperature: f64,
    pub(crate) change: f64,
    pub(crate) slope: Option<f64>,
}

impl Group<'_> {
    fn after(&self, change: f64, work: &Work) -> Result<WallPoint, WallExchangeError> {
        self.after_state(self.initial_wall, change, work)
    }

    // The same original storage, now about a candidate composition-stage state.
    // Neither gas contact nor wall links store a second thermal inventory.
    pub(crate) fn after_state(
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

    fn evaluate<G: WallGasCaloric>(
        &mut self,
        gas: &impl Fn(SharedWallGas, f64) -> G,
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
                original: gas(port.gas, wall.temperature),
                work,
            };
            work.contacts.set(work.contacts.get().checked_add(1).ok_or(
                WallExchangeError::Invalid("shared wall contact count overflow"),
            )?);
            let update = exchange_caloric(
                &caloric,
                state[0],
                state[1],
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
                let endpoint = caloric.primitive(state[0], update.energy, cancelled)?;
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

    pub(crate) fn solve<G: WallGasCaloric>(
        &mut self,
        gas: &impl Fn(SharedWallGas, f64) -> G,
        duration: Time,
        control: WallControl,
        work: &Work,
        initial_wall: WallThermal,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallCellReport, WallExchangeError> {
        self.initial_wall = initial_wall;
        let wall_initial = self.after(0.0, work)?.temperature;
        let initial = self.evaluate(gas, 0.0, duration, control, work, cancelled)?;
        let mut iterations = 0;
        let mut accepted_change = 0.0;
        if initial.0.abs() > control.tolerance(initial.2) {
            let mut minimum = wall_initial;
            let mut maximum = wall_initial;
            for (port, state) in self.accepted.ports.iter().zip(&self.initial) {
                let original = MeteredGas {
                    original: gas(port.gas, wall_initial),
                    work,
                };
                let t = original
                    .primitive(state[0], state[1], cancelled)?
                    .temperature;
                minimum = minimum.min(t);
                maximum = maximum.max(t);
            }
            let (mut low, mut high) = if let Some(source) = &self.accepted.storage.heat_storage {
                let contact = WallContact {
                    drag: WallDrag::None,
                    thermal: self.initial_wall,
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
                } = self.initial_wall
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

fn accumulate<const N: usize>(
    sum: &mut [f64; N],
    correction: &mut [f64; N],
    values: [f64; N],
) -> Result<(), WallExchangeError> {
    for (k, value) in values.into_iter().enumerate() {
        let increment = value - correction[k];
        let next = sum[k] + increment;
        if !next.is_finite() {
            return Err(WallExchangeError::Invalid(
                "unrepresentable combined wall source",
            ));
        }
        correction[k] = (next - sum[k]) - increment;
        sum[k] = next;
    }
    Ok(())
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WallExchangeError> {
    if cancelled() {
        Err(WallExchangeError::Cancelled)
    } else {
        Ok(())
    }
}
