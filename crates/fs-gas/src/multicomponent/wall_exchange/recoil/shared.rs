//! One finite axial support coupled to several fixed gas amounts.
//!
//! Common linear resistance uses exact collective and deviation modes. Other
//! laws compose exact pair maps symmetrically with bounded time refinement;
//! a resolved-rate cap prevents falsely converged saturated pair projections.

use super::{RecoilControl, RecoilLaw, SlidingSupport, Time, WallExchangeError, mechanics};
use super::{nonnegative, poll, positive};

mod collective_linear;

/// One fixed gas amount and its frozen acceleration law against the support.
/// Coefficients are active even when this cell initially has zero slip.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedRecoilCell {
    /// Positive gas density [kg/m³].
    pub density_kg_m3: f64,
    /// Original axial momentum density [kg/(m² s)].
    pub momentum_kg_m2_s: f64,
    /// Positive fixed gas volume [m³].
    pub volume_m3: f64,
    /// Nonnegative linear coefficient a in du/dt = -a*s - b*abs(s)*s [1/s].
    pub linear_rate_per_s: f64,
    /// Nonnegative quadratic coefficient b [1/m].
    pub quadratic_per_m: f64,
}

/// Per-cell exchanges over the complete common interval. All slip heat is
/// retained in that gas; the caller updates total energy by work/volume only.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SharedRecoilCellUpdate {
    /// Returned axial momentum density [kg/(m² s)].
    pub momentum: f64,
    /// Actual complete endpoint gas momentum change [N s].
    pub gas_impulse_ns: f64,
    /// Integrated support work entering this gas amount [J].
    pub work_to_gas_j: f64,
    /// Nonnegative physical relative-slip loss [J], integrated analytically
    /// for common linear rates or through the refined exact pair subflows.
    pub slip_dissipation_j: f64,
}

/// Temporal splitting accuracy and bounded work, separate from the original
/// RecoilControl dimensional momentum/energy conservation budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedRecoilControl {
    /// Absolute whole-interval velocity error estimate allowance [m/s].
    pub absolute_velocity_tolerance_m_s: f64,
    /// Absolute whole-interval per-cell heat/work error estimate allowance [J].
    pub absolute_energy_tolerance_j: f64,
    /// Relative allowance, at most 1e-3. Velocity uses the initial velocity
    /// hull width; each cell's energy uses mass*width²/2, never a bulk boost.
    pub relative_tolerance: f64,
    /// Maximum attempted adaptive intervals, including rejected attempts,
    /// in 1..=1_048_576. Each attempt uses six bounded pair sweeps.
    pub maximum_steps: usize,
}

impl Default for SharedRecoilControl {
    fn default() -> Self {
        Self {
            absolute_velocity_tolerance_m_s: 1e-7,
            absolute_energy_tolerance_j: 1e-7,
            relative_tolerance: 1e-6,
            maximum_steps: 16_384,
        }
    }
}

impl SharedRecoilControl {
    pub(crate) fn validate(self) -> Result<(), WallExchangeError> {
        if !positive(self.absolute_velocity_tolerance_m_s)
            || !positive(self.absolute_energy_tolerance_j)
            || !positive(self.relative_tolerance)
            || self.relative_tolerance > 1e-3
            || !(1..=1_048_576).contains(&self.maximum_steps)
        {
            return Err(WallExchangeError::Invalid(
                "shared support integration controls",
            ));
        }
        Ok(())
    }
}

/// Complete candidate result. Borrowing this result prevents workspace reuse;
/// no caller-owned gas or support state is changed on success or refusal.
#[derive(Debug)]
pub struct SharedRecoilUpdate<'a> {
    /// Exchanges in exactly the input cell order.
    pub cells: &'a [SharedRecoilCellUpdate],
    /// One support after the complete simultaneous interval.
    pub support: SlidingSupport,
    /// Actual complete support momentum change [N s].
    pub support_impulse_ns: f64,
    /// Actual complete support kinetic-energy change [J].
    pub support_kinetic_energy_change_j: f64,
    /// Actual gas plus support endpoint momentum change [N s].
    pub momentum_residual_ns: f64,
    /// Actual gas/support kinetic change plus all slip heat [J].
    pub energy_residual_j: f64,
    /// Sum of accepted local maximum velocity estimates [m/s], not a bound.
    /// Exact pair/collective modes have no splitting error and return zero;
    /// this is not a floating-point error certificate.
    pub velocity_error_estimate_m_s: f64,
    /// Sum of accepted local maximum per-cell heat/work estimates [J], not a bound.
    pub energy_error_estimate_j: f64,
    /// Accepted intervals: one exact pair/collective interval, or adaptive
    /// intervals that each use the two refined Strang steps.
    pub accepted_steps: usize,
    /// Attempted intervals, including all discarded refinements.
    pub attempted_steps: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    high: f64,
    low: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.high + value;
        self.low += if self.high.abs() >= value.abs() {
            (self.high - next) + value
        } else {
            (value - next) + self.high
        };
        self.high = next;
    }
    fn value(self) -> f64 {
        self.high + self.low
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Stage {
    momentum: f64,
    // Work is accumulated in the initial support frame. Restoring U0*Δp at
    // the end retains weak transfers beside a large laboratory-frame boost.
    work: Sum,
    heat: Sum,
}

/// Reusable O(number of cells) storage for one shared axial support. Allocate
/// at binding, then reuse without allocation on every original source half.
pub struct SharedRecoilWorkspace {
    accepted: Vec<Stage>,
    coarse: Vec<Stage>,
    fine: Vec<Stage>,
    updates: Vec<SharedRecoilCellUpdate>,
}

impl SharedRecoilWorkspace {
    /// Allocate all candidate and refinement storage once. Cell count is fixed.
    pub fn new(
        cell_count: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        Ok(Self {
            accepted: allocated(cell_count, Stage::default(), &mut cancelled)?,
            coarse: allocated(cell_count, Stage::default(), &mut cancelled)?,
            fine: allocated(cell_count, Stage::default(), &mut cancelled)?,
            updates: allocated(
                cell_count,
                SharedRecoilCellUpdate::default(),
                &mut cancelled,
            )?,
        })
    }

    /// Advance du_i/dt = -a_i*s_i - b_i*|s_i|*s_i and
    /// M*dU/dt = sum(m_i*(a_i*s_i + b_i*|s_i|*s_i)), with s_i=u_i-U.
    ///
    /// Positive exact pair flows preserve the initial velocity hull and
    /// nonnegative physical slip heat. Palindromic Strang steps are compared
    /// with two half steps; the refined solution is accepted without
    /// extrapolation. Estimates cover every gas velocity, the support velocity,
    /// per-cell heat and work in the initial support frame. Local allowances
    /// consume their duration fraction of the whole-interval allowance.
    ///
    /// For multiple active cells, h*(max(k_i)+sum(m_i*k_i/M)) <= 1/2 with
    /// k_i=a_i+2*b_i*initial_hull_width resolves every pair before its error
    /// estimate is trusted. Very stiff multi-cell problems can exhaust the
    /// explicit work budget; they never silently become sequential equilibrium
    /// projections. A single active cell uses the original exact pair once.
    /// A common positive linear rate with no quadratic terms instead uses the
    /// exact collective/deviation modes in one O(cells) pass, including each
    /// cell's physical heat and work, without a stiffness step restriction.
    ///
    /// These are temporal error estimates, not certified error enclosures.
    /// No moving volume, thermal law, external drive or branch-axis mapping is
    /// inferred. Cancellation/refusal only changes private workspace scratch.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'a>(
        &'a mut self,
        cells: &[SharedRecoilCell],
        duration: Time,
        support: SlidingSupport,
        control: SharedRecoilControl,
        balance: RecoilControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SharedRecoilUpdate<'a>, WallExchangeError> {
        poll(&mut cancelled)?;
        control.validate()?;
        balance.wall.validate()?;
        let dt = duration.value();
        if cells.len() != self.accepted.len()
            || !positive(support.mass_kg)
            || !support.velocity_m_s.is_finite()
            || !positive(balance.momentum_tolerance_ns)
            || !nonnegative(dt)
        {
            return Err(WallExchangeError::Invalid(
                "shared support shape, mass, clock or balance budget",
            ));
        }
        let mut minimum_slip = 0.0_f64;
        let mut maximum_slip = 0.0_f64;
        let mut active = 0;
        let mut only_active = 0;
        for (i, cell) in cells.iter().copied().enumerate() {
            poll(&mut cancelled)?;
            let velocity = cell.momentum_kg_m2_s / cell.density_kg_m3;
            let slip = velocity - support.velocity_m_s;
            if !positive(cell.density_kg_m3)
                || !positive(cell.volume_m3)
                || !positive(cell.density_kg_m3 * cell.volume_m3)
                || !cell.momentum_kg_m2_s.is_finite()
                || !slip.is_finite()
                || !nonnegative(cell.linear_rate_per_s)
                || !nonnegative(cell.quadratic_per_m)
            {
                return Err(WallExchangeError::Invalid(
                    "shared support gas mass, momentum or drag",
                ));
            }
            self.accepted[i] = Stage {
                momentum: cell.momentum_kg_m2_s,
                ..Stage::default()
            };
            if cell.linear_rate_per_s > 0.0 || cell.quadratic_per_m > 0.0 {
                minimum_slip = minimum_slip.min(slip);
                maximum_slip = maximum_slip.max(slip);
                active += 1;
                only_active = i;
            }
        }
        let width = maximum_slip - minimum_slip;
        if !nonnegative(width) {
            return Err(WallExchangeError::Invalid(
                "shared support initial velocity hull",
            ));
        }
        let mut next_support = support;
        let mut accepted_steps = 0;
        let mut attempted_steps = 0;
        let mut single_work = None;
        let mut velocity_estimate = Sum::default();
        let mut energy_estimate = Sum::default();
        if dt > 0.0 && width > 0.0 && active > 0 {
            if active == 1 {
                single_work = Some(pair(
                    cells[only_active],
                    &mut self.accepted[only_active],
                    &mut next_support,
                    dt,
                    support.velocity_m_s,
                )?);
                accepted_steps = 1;
                attempted_steps = 1;
            } else if let Some(updated) =
                collective_linear::advance(cells, &mut self.accepted, dt, support, &mut cancelled)?
            {
                next_support = updated;
                accepted_steps = 1;
                attempted_steps = 1;
            } else {
                let mut weighted_rate = Sum::default();
                let mut largest_rate = 0.0_f64;
                for cell in cells {
                    poll(&mut cancelled)?;
                    if cell.linear_rate_per_s == 0.0 && cell.quadratic_per_m == 0.0 {
                        continue;
                    }
                    let rate = cell.linear_rate_per_s + (cell.quadratic_per_m * width) * 2.0;
                    let weighted = (cell.density_kg_m3 * cell.volume_m3 / support.mass_kg) * rate;
                    if !nonnegative(rate) || !nonnegative(weighted) {
                        return Err(WallExchangeError::Invalid(
                            "shared support resolved relaxation rate",
                        ));
                    }
                    largest_rate = largest_rate.max(rate);
                    weighted_rate.add(weighted);
                }
                let rate = largest_rate + weighted_rate.value();
                let maximum_step = 0.5 / rate;
                if !positive(rate) || !positive(maximum_step) {
                    return Err(WallExchangeError::Invalid(
                        "shared support resolved time interval",
                    ));
                }
                let velocity_allowance =
                    control.absolute_velocity_tolerance_m_s + control.relative_tolerance * width;
                if !positive(velocity_allowance) {
                    return Err(WallExchangeError::Invalid(
                        "shared support velocity allowance",
                    ));
                }
                let mut elapsed = 0.0;
                let mut step = dt.min(maximum_step);
                while elapsed < dt {
                    poll(&mut cancelled)?;
                    if attempted_steps == control.maximum_steps {
                        return Err(WallExchangeError::Invalid(
                            "shared support integration work exhausted",
                        ));
                    }
                    let remaining = dt - elapsed;
                    step = step.min(remaining).min(maximum_step);
                    if !positive(step) || elapsed + step == elapsed || step * 0.25 == 0.0 {
                        return Err(WallExchangeError::Invalid("shared support time resolution"));
                    }
                    attempted_steps += 1;
                    for i in 0..cells.len() {
                        poll(&mut cancelled)?;
                        self.coarse[i] = Stage {
                            momentum: self.accepted[i].momentum,
                            ..Stage::default()
                        };
                        self.fine[i] = self.coarse[i];
                    }
                    let mut coarse_support = next_support;
                    let mut fine_support = next_support;
                    sweep(
                        cells,
                        &mut self.coarse,
                        &mut coarse_support,
                        step,
                        support.velocity_m_s,
                        &mut cancelled,
                    )?;
                    sweep(
                        cells,
                        &mut self.fine,
                        &mut fine_support,
                        step * 0.5,
                        support.velocity_m_s,
                        &mut cancelled,
                    )?;
                    sweep(
                        cells,
                        &mut self.fine,
                        &mut fine_support,
                        step * 0.5,
                        support.velocity_m_s,
                        &mut cancelled,
                    )?;
                    let fraction = step / dt;
                    let velocity_budget = velocity_allowance * fraction;
                    if !positive(fraction) || !positive(velocity_budget) {
                        return Err(WallExchangeError::Invalid(
                            "shared support unresolved velocity allowance",
                        ));
                    }
                    let mut velocity_error =
                        (coarse_support.velocity_m_s - fine_support.velocity_m_s).abs() / 3.0;
                    let mut energy_error = 0.0_f64;
                    let mut ratio = velocity_error / velocity_budget;
                    for (i, cell) in cells.iter().enumerate() {
                        poll(&mut cancelled)?;
                        if cell.linear_rate_per_s == 0.0 && cell.quadratic_per_m == 0.0 {
                            continue;
                        }
                        let velocity = ((self.coarse[i].momentum - self.fine[i].momentum)
                            / cell.density_kg_m3)
                            .abs()
                            / 3.0;
                        let work =
                            (self.coarse[i].work.value() - self.fine[i].work.value()).abs() / 3.0;
                        let heat =
                            (self.coarse[i].heat.value() - self.fine[i].heat.value()).abs() / 3.0;
                        let scale = ((0.5 * (cell.density_kg_m3 * cell.volume_m3)) * width) * width;
                        let energy_allowance = control.absolute_energy_tolerance_j
                            + control.relative_tolerance * scale;
                        let energy_budget = energy_allowance * fraction;
                        if !nonnegative(scale)
                            || !positive(energy_allowance)
                            || !positive(energy_budget)
                            || !velocity.is_finite()
                            || !work.is_finite()
                            || !heat.is_finite()
                        {
                            return Err(WallExchangeError::Invalid(
                                "shared support error estimate scale",
                            ));
                        }
                        velocity_error = velocity_error.max(velocity);
                        energy_error = energy_error.max(work.max(heat));
                        ratio = ratio
                            .max(velocity / velocity_budget)
                            .max(work.max(heat) / energy_budget);
                    }
                    if !ratio.is_finite() {
                        return Err(WallExchangeError::Invalid(
                            "shared support unresolved error allowance",
                        ));
                    }
                    if ratio <= 1.0 {
                        for i in 0..cells.len() {
                            poll(&mut cancelled)?;
                            self.accepted[i].momentum = self.fine[i].momentum;
                            self.accepted[i].work.add(self.fine[i].work.value());
                            self.accepted[i].heat.add(self.fine[i].heat.value());
                        }
                        next_support = fine_support;
                        accepted_steps += 1;
                        velocity_estimate.add(velocity_error);
                        energy_estimate.add(energy_error);
                        elapsed = if step == remaining {
                            dt
                        } else {
                            elapsed + step
                        };
                        // Dyadic refinement keeps the same deterministic arithmetic
                        // and avoids a platform-dependent fractional-power controller.
                        if ratio < 0.125 {
                            step = (step * 2.0).min(maximum_step);
                        }
                    } else {
                        step *= 0.5;
                    }
                }
            }
        }

        let support_impulse = support.mass_kg * (next_support.velocity_m_s - support.velocity_m_s);
        let support_kinetic =
            support_impulse * (0.5 * next_support.velocity_m_s + 0.5 * support.velocity_m_s);
        let mut momentum_residual = Sum::default();
        momentum_residual.add(support_impulse);
        let mut energy_residual = Sum::default();
        energy_residual.add(support_kinetic);
        let mut exchange_scale = Sum::default();
        for (i, cell) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let state = self.accepted[i];
            let impulse = (state.momentum - cell.momentum_kg_m2_s) * cell.volume_m3;
            let work = single_work
                .filter(|_| i == only_active)
                .unwrap_or_else(|| state.work.value() + support.velocity_m_s * impulse);
            let heat = state.heat.value();
            let kinetic = impulse
                * (0.5 * (state.momentum / cell.density_kg_m3)
                    + 0.5 * (cell.momentum_kg_m2_s / cell.density_kg_m3));
            let local_residual = (kinetic + heat) - work;
            if ![impulse, work, heat, kinetic, local_residual]
                .into_iter()
                .all(f64::is_finite)
                || heat < 0.0
            {
                return Err(WallExchangeError::Invalid(
                    "shared support cell work or heat",
                ));
            }
            if local_residual.abs() > balance.wall.tolerance(work.abs().max(heat)) {
                return Err(WallExchangeError::Balance {
                    energy_residual_j: local_residual,
                    thermal_residual_j: 0.0,
                });
            }
            self.updates[i] = SharedRecoilCellUpdate {
                momentum: state.momentum,
                gas_impulse_ns: impulse,
                work_to_gas_j: work,
                slip_dissipation_j: heat,
            };
            momentum_residual.add(impulse);
            energy_residual.add(kinetic);
            energy_residual.add(heat);
            exchange_scale.add(work.abs().max(heat));
        }
        let momentum_residual = momentum_residual.value();
        let energy_residual = energy_residual.value();
        let scale = exchange_scale.value().max(support_kinetic.abs());
        if ![
            support_impulse,
            support_kinetic,
            momentum_residual,
            energy_residual,
            scale,
        ]
        .into_iter()
        .all(f64::is_finite)
            || momentum_residual.abs() > balance.momentum_tolerance_ns
        {
            return Err(WallExchangeError::Invalid(
                "shared support complete momentum balance",
            ));
        }
        if energy_residual.abs() > balance.wall.tolerance(scale) {
            return Err(WallExchangeError::Balance {
                energy_residual_j: energy_residual,
                thermal_residual_j: 0.0,
            });
        }
        if !nonnegative(velocity_estimate.value()) || !nonnegative(energy_estimate.value()) {
            return Err(WallExchangeError::Invalid(
                "shared support complete error estimate",
            ));
        }
        poll(&mut cancelled)?;
        Ok(SharedRecoilUpdate {
            cells: &self.updates,
            support: next_support,
            support_impulse_ns: support_impulse,
            support_kinetic_energy_change_j: support_kinetic,
            momentum_residual_ns: momentum_residual,
            energy_residual_j: energy_residual,
            velocity_error_estimate_m_s: velocity_estimate.value(),
            energy_error_estimate_j: energy_estimate.value(),
            accepted_steps,
            attempted_steps,
        })
    }
}

fn allocated<T: Copy>(
    count: usize,
    value: T,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<T>, WallExchangeError> {
    poll(cancelled)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| WallExchangeError::Invalid("shared support workspace allocation"))?;
    for _ in 0..count {
        poll(cancelled)?;
        values.push(value);
    }
    Ok(values)
}

fn pair(
    cell: SharedRecoilCell,
    state: &mut Stage,
    support: &mut SlidingSupport,
    dt: f64,
    frame_velocity: f64,
) -> Result<f64, WallExchangeError> {
    if cell.linear_rate_per_s == 0.0 && cell.quadratic_per_m == 0.0 {
        return Ok(0.0);
    }
    let motion = mechanics::advance(
        state.momentum,
        cell.density_kg_m3,
        state.momentum / cell.density_kg_m3,
        cell.volume_m3,
        dt,
        *support,
        RecoilLaw {
            linear_rate_per_s: cell.linear_rate_per_s,
            quadratic_per_m: cell.quadratic_per_m,
            gas_heat_fraction: 1.0,
        },
    )?;
    let relative_support = support.velocity_m_s - frame_velocity;
    let support_delta = -motion.gas_impulse_ns / support.mass_kg;
    state
        .work
        .add(motion.gas_impulse_ns * (relative_support + 0.5 * support_delta));
    state.heat.add(motion.dissipation_j);
    state.momentum = motion.momentum;
    support.velocity_m_s = motion.support_velocity_m_s;
    Ok(motion.work_to_gas_j)
}

fn sweep(
    cells: &[SharedRecoilCell],
    states: &mut [Stage],
    support: &mut SlidingSupport,
    dt: f64,
    frame_velocity: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), WallExchangeError> {
    for (cell, state) in cells.iter().zip(states.iter_mut()) {
        poll(cancelled)?;
        pair(*cell, state, support, dt * 0.5, frame_velocity)?;
    }
    for (cell, state) in cells.iter().zip(states.iter_mut()).rev() {
        poll(cancelled)?;
        pair(*cell, state, support, dt * 0.5, frame_velocity)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
