//! Prepared finite-rate ion Newton solve in the original charge-nullspace.
//!
//! Each free species has its own transfer. One charged pivot enforces charge;
//! absent pools and closed channels are fixed exactly. The native LU operates
//! on scaled original backward-Euler rows, using the membrane's analytic
//! response after its floating-potential solve. Every original row, including
//! the eliminated pivot and both compartments, still owns acceptance.

use super::*;
use fs_la::LuWorkspace;
use fs_material::osmosis::electrodiffusion::ElectrodiffusionDirection;

mod seed;

#[derive(Clone, Debug)]
pub(super) struct Workspace {
    active: [bool; CAPACITY],
    scales: [f64; CAPACITY],
    dimension: usize,
    jacobian: Box<[f64]>,
    factor: LuWorkspace,
}

struct Chart {
    lower: IonicSpeciesValues,
    upper: IonicSpeciesValues,
    pivot: Option<usize>,
    free: [usize; CAPACITY],
    count: usize,
}

impl Workspace {
    pub(super) fn new(
        membrane: &ElectrodiffusiveMembrane,
        totals: &IonicSpeciesValues,
        base_bytes: usize,
        maximum_bytes: usize,
    ) -> Result<Self> {
        let mut active = [false; CAPACITY];
        let mut scales = [0.0; CAPACITY];
        let mut count = 0;
        let mut charged = false;
        for i in 0..totals.len() {
            active[i] = totals[i] > 0.0 && membrane.permeances()[i].value() > 0.0;
            scales[i] = totals[i];
            if active[i] {
                count += 1;
                charged |= membrane.valences()[i] != 0;
            }
        }
        let dimension = count - usize::from(charged);
        let storage = storage_bytes(dimension);
        let required = base_bytes
            .checked_add(storage)
            .ok_or_else(|| invalid("ionic workspace byte overflow"))?;
        if required > maximum_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_bytes,
            }
            .into());
        }
        let mut jacobian = Vec::new();
        jacobian
            .try_reserve_exact(dimension * dimension)
            .map_err(|_| invalid("could not prepare the bounded ionic Jacobian"))?;
        jacobian.resize(dimension * dimension, 0.0);
        let factor = LuWorkspace::new(dimension)
            .map_err(|_| invalid("could not prepare the native ionic LU"))?;
        Ok(Self {
            active,
            scales,
            dimension,
            jacobian: jacobian.into_boxed_slice(),
            factor,
        })
    }

    pub(super) fn storage_bytes(&self) -> usize {
        storage_bytes(self.dimension)
    }

    fn chart(&self, pair: &IonicOsmoticPair, volumes: [f64; 2]) -> Result<Chart> {
        let n = pair.species_count();
        let z = pair.membrane.valences();
        let mut lower = IonicSpeciesValues::zeroed(n);
        let mut upper = lower;
        for i in 0..n {
            let (left, right) = (
                pair.compartments[0].solution().solute_moles()[i],
                pair.compartments[1].solution().solute_moles()[i],
            );
            let (mut lo, mut hi) = (-right, left);
            if let Some(model) = pair.compartments[0].solution().nonideal_model() {
                let cap = model.maximum_concentrations_mol_m3()[i];
                [lo, hi] = supported_cap_bounds(left, right, volumes, cap)?;
            }
            if !self.active[i] {
                if lo > 0.0 || hi < 0.0 {
                    return Err(OsmoticNetworkError::Unbracketed.into());
                }
                (lo, hi) = (0.0, 0.0);
            }
            if !lo.is_finite() || !hi.is_finite() || lo > hi {
                return Err(OsmoticNetworkError::Unbracketed.into());
            }
            lower.values[i] = lo;
            upper.values[i] = hi;
        }
        // At a one-sided conducting support, charge continuity can freeze the
        // entire charged pool. Recognize the exact zero extremum without
        // underflow-prone weighted sums or a fabricated concentration floor.
        let minimum_zero = (0..n).all(|i| {
            z[i] == 0
                || if z[i] > 0 {
                    lower[i] == 0.0
                } else {
                    upper[i] == 0.0
                }
        });
        let maximum_zero = (0..n).all(|i| {
            z[i] == 0
                || if z[i] > 0 {
                    upper[i] == 0.0
                } else {
                    lower[i] == 0.0
                }
        });
        if minimum_zero || maximum_zero {
            for i in 0..n {
                if z[i] != 0 {
                    lower.values[i] = 0.0;
                    upper.values[i] = 0.0;
                }
            }
        }
        let select = |varying: bool| {
            // Neutral or frozen backgrounds cannot set the charge scale: a
            // huge neutral pool would otherwise erase charged trace weights.
            let eligible =
                |i: usize| self.active[i] && z[i] != 0 && (!varying || lower[i] < upper[i]);
            let maximum = (0..n)
                .filter(|&i| eligible(i))
                .map(|i| self.scales[i])
                .fold(0.0_f64, f64::max);
            let mut selected = None;
            let mut weight = -1.0;
            for i in 0..n {
                if eligible(i) {
                    let candidate = (self.scales[i] / maximum) * f64::from(z[i]).abs();
                    if candidate > weight {
                        weight = candidate;
                        selected = Some(i);
                    }
                }
            }
            selected
        };
        let pivot = select(true).or_else(|| select(false));
        let mut free = [0; CAPACITY];
        let mut count = 0;
        for i in 0..n {
            if self.active[i] && lower[i] < upper[i] && Some(i) != pivot {
                free[count] = i;
                count += 1;
            }
        }
        if count > self.dimension {
            return Err(invalid("ionic chart exceeds prepared LU dimension"));
        }
        Ok(Chart {
            lower,
            upper,
            pivot,
            free,
            count,
        })
    }

    pub(super) fn solve(
        &mut self,
        pair: &mut IonicOsmoticPair,
        probe: Probe,
        controls: IonicStepControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        let start = work.salt;
        let chart = self.chart(pair, probe.volumes)?;
        let n = pair.species_count();
        let mut transfers = IonicSpeciesValues::zeroed(n);
        seed::feasible_seed(
            &chart.lower,
            &chart.upper,
            pair.membrane.valences(),
            chart.pivot,
            transfers.as_mut_slice(),
        )?;
        let mut trial =
            pair.sample_transfers(probe, &transfers, controls, work, start, cancelled)?;
        let mut restarted = false;
        loop {
            poll(cancelled)?;
            if trial.ion_error <= 1.0 {
                return Ok(trial);
            }
            if chart.count == 0 || work.salt - start >= controls.maximum_root_evaluations {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }

            // A single fresh, counted original response supplies all analytic
            // columns. No perturbed voltages or finite-difference force model.
            let mut rhs = [0.0; CAPACITY];
            {
                spend(work, start, controls, cancelled)?;
                let remaining = controls
                    .maximum_potential_evaluations
                    .checked_sub(work.potential)
                    .filter(|n| *n > 0)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                let mut electrical = controls.electrical;
                electrical.maximum_potential_evaluations =
                    electrical.maximum_potential_evaluations.min(remaining);
                let mut flux = IonicSpeciesValues::zeroed(n);
                let response = pair.membrane.flux_response_into_counted(
                    ElectrodiffusionProblem {
                        left: &pair.candidates[0],
                        right: &pair.candidates[1],
                        pressure_drop_pa: probe.pressure_drop,
                        trapped_particles_mol: pair.trapped_particles_mol,
                    },
                    electrical,
                    flux.as_mut_slice(),
                    &mut work.potential,
                    &mut *cancelled,
                );
                let response = match response {
                    Err(ElectrodiffusionError::Cancelled) => {
                        return Err(OsmoticNetworkError::Cancelled.into());
                    }
                    Err(_) if work.potential >= controls.maximum_potential_evaluations => {
                        return Err(OsmoticNetworkError::WorkBudget.into());
                    }
                    other => other?,
                };
                self.jacobian.fill(0.0);
                for row in 0..self.dimension {
                    self.jacobian[row * self.dimension + row] = 1.0;
                }
                for row in 0..chart.count {
                    let i = chart.free[row];
                    rhs[row] = finite(
                        (integrated(probe.dt, pair.area_m2, flux[i]) - transfers[i])
                            / self.scales[i],
                    )?;
                }
                for column in 0..chart.count {
                    poll(cancelled)?;
                    let i = chart.free[column];
                    let scale = self.scales[i];
                    let mut left = IonicSpeciesValues::zeroed(n);
                    let mut right = left;
                    left.values[i] = -scale;
                    right.values[i] = scale;
                    if let Some(pivot) = chart.pivot {
                        let value = finite(
                            (f64::from(pair.membrane.valences()[i])
                                / f64::from(pair.membrane.valences()[pivot]))
                                * scale,
                        )?;
                        left.values[pivot] = value;
                        right.values[pivot] = -value;
                    }
                    let mut derivative = IonicSpeciesValues::zeroed(n);
                    response.directional_derivative_into(
                        ElectrodiffusionDirection {
                            left_amounts_mol: &left,
                            right_amounts_mol: &right,
                            volumes_m3: [0.0; 2],
                            pressure_drop_pa: 0.0,
                        },
                        derivative.as_mut_slice(),
                        &mut *cancelled,
                    )?;
                    for row in 0..chart.count {
                        let j = chart.free[row];
                        let original = if row == column { scale } else { 0.0 };
                        self.jacobian[row * self.dimension + column] = finite(
                            (original - integrated(probe.dt, pair.area_m2, derivative[j]))
                                / self.scales[j],
                        )?;
                    }
                }
            }
            // The immutable original-state response ends before probing new
            // candidate inventories. Every trial below is an independent replay.
            let mut update = [0.0; CAPACITY];
            self.factor
                .solve_into(
                    &self.jacobian,
                    &rhs[..self.dimension],
                    &mut update[..self.dimension],
                )
                .map_err(|_| {
                    IonicOsmoticError::Network(OsmoticNetworkError::NonlinearConvergence)
                })?;
            let mut direction = IonicSpeciesValues::zeroed(n);
            for column in 0..chart.count {
                let i = chart.free[column];
                direction.values[i] = finite(self.scales[i] * update[column])?;
            }
            recover_pivot(
                direction.as_mut_slice(),
                pair.membrane.valences(),
                chart.pivot,
            )?;
            let mut limit = 1.0_f64;
            for i in 0..n {
                if direction[i] > 0.0 {
                    limit = limit.min((chart.upper[i] - transfers[i]) / direction[i]);
                }
                if direction[i] < 0.0 {
                    limit = limit.min((chart.lower[i] - transfers[i]) / direction[i]);
                }
            }
            let mut alpha = if limit < 1.0 { 0.95 * limit } else { 1.0 };
            let mut accepted = None;
            while positive(alpha) {
                let mut candidate = transfers;
                for &i in &chart.free[..chart.count] {
                    candidate.values[i] = transfers[i] + alpha * direction[i];
                }
                recover_pivot(
                    candidate.as_mut_slice(),
                    pair.membrane.valences(),
                    chart.pivot,
                )?;
                if candidate == transfers {
                    break;
                }
                if (0..n).any(|i| candidate[i] < chart.lower[i] || candidate[i] > chart.upper[i]) {
                    spend(work, start, controls, cancelled)?;
                    alpha *= 0.5;
                    continue;
                }
                match pair.sample_transfers(probe, &candidate, controls, work, start, cancelled) {
                    Ok(next)
                        if next.ion_error <= 1.0
                            || next.ion_error < trial.ion_error * (1.0 - 1e-4 * alpha) =>
                    {
                        accepted = Some((candidate, next));
                        break;
                    }
                    Err(error)
                        if fatal(&error)
                            || work.salt - start >= controls.maximum_root_evaluations =>
                    {
                        return Err(error);
                    }
                    _ => {}
                }
                alpha *= 0.5;
            }
            if let Some((candidate, next)) = accepted {
                transfers = candidate;
                trial = next;
                continue;
            }
            // A boundary Newton direction can have no admissible step. One
            // interior feasibility restart is numerical only; the original
            // finite-rate equations must still accept every species afterward.
            if restarted {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            restarted = true;
            let mut center = IonicSpeciesValues::from_fn(n, |i| {
                let (lo, hi) = (chart.lower[i], chart.upper[i]);
                if lo.signum() == hi.signum() {
                    lo + 0.5 * (hi - lo)
                } else {
                    0.5 * lo + 0.5 * hi
                }
            });
            seed::feasible_seed(
                &chart.lower,
                &chart.upper,
                pair.membrane.valences(),
                chart.pivot,
                center.as_mut_slice(),
            )?;
            if center == transfers {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            transfers = center;
            trial = pair.sample_transfers(probe, &transfers, controls, work, start, cancelled)?;
        }
    }
}

fn storage_bytes(n: usize) -> usize {
    core::mem::size_of::<Workspace>()
        + (2 * n * n + n) * core::mem::size_of::<f64>()
        + n * core::mem::size_of::<usize>()
}

/// Only refused rounded endpoints move inward. The original amount/volume
/// division owns support, including a valid zero or singleton transfer. The
/// correction uses inventory differences rather than walking tiny transfer
/// ulps beside a large background amount.
fn supported_cap_bounds(left: f64, right: f64, volumes: [f64; 2], cap: f64) -> Result<[f64; 2]> {
    let mut bounds = [
        (-right).max(left - cap * volumes[0]),
        left.min(cap * volumes[1] - right),
    ];
    if left / volumes[0] <= cap {
        bounds[0] = bounds[0].min(0.0);
    }
    if right / volumes[1] <= cap {
        bounds[1] = bounds[1].max(0.0);
    }
    for side in 0..2 {
        let amount = if side == 0 { left } else { right };
        let sign = if side == 0 { -1.0 } else { 1.0 };
        let mut admitted = cap * volumes[side];
        for _ in 0..4 {
            if admitted / volumes[side] <= cap {
                break;
            }
            admitted = admitted.next_down();
        }
        if admitted.is_finite() && admitted / volumes[side] > cap {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
        for _ in 0..4 {
            let actual = amount + sign * bounds[side];
            if actual / volumes[side] <= cap {
                break;
            }
            let corrected = bounds[side] - sign * (actual - admitted);
            bounds[side] = if side == 0 {
                corrected.max(bounds[side].next_up())
            } else {
                corrected.min(bounds[side].next_down())
            };
        }
        let actual = amount + sign * bounds[side];
        if !nonnegative(actual) || actual / volumes[side] > cap {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
    }
    Ok(bounds)
}

fn recover_pivot(transfers: &mut [f64], valences: &[i32], pivot: Option<usize>) -> Result<()> {
    if let Some(pivot) = pivot {
        let (mut sum, mut correction) = (0.0, 0.0);
        let z = f64::from(valences[pivot]);
        for (i, &value) in transfers.iter().enumerate() {
            if i != pivot {
                compensated_add(
                    &mut sum,
                    &mut correction,
                    -(f64::from(valences[i]) / z) * value,
                );
            }
        }
        transfers[pivot] = finite(sum)?;
    }
    Ok(())
}
fn spend(
    work: &mut Work,
    start: usize,
    controls: IonicStepControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    if work.salt - start >= controls.maximum_root_evaluations {
        return Err(OsmoticNetworkError::NonlinearConvergence.into());
    }
    bump(&mut work.salt)
}
