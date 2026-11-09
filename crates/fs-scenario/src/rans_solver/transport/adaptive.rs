//! Bounded physical-window integration of the original implicit channel steps.
//!
//! One full step is compared with two half steps. Only the latter trajectory
//! is accepted, without Richardson extrapolation or positivity clipping. This
//! estimates local temporal error; it does not certify global or model error.

use super::{
    RansChannelGrid, RansFieldState, RansHeatLoad, RansSolver, RansThermalTransportReport,
    RansTransportError, RansTransportReport, RansTransportWorkspace, checkpoint,
};

/// Absolute local-error allowances, in the physical units of each field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansFieldTolerances {
    /// Mean velocity allowance [m/s].
    pub velocity_m_s: f64,
    /// Turbulent kinetic energy allowance [m²/s²].
    pub kinetic_energy_m2_s2: f64,
    /// Modified dissipation allowance [m²/s³].
    pub dissipation_m2_s3: f64,
    /// Eddy viscosity allowance [m²/s].
    pub eddy_viscosity_m2_s: f64,
    /// Temperature allowance [K], used only when heat is advanced.
    pub temperature_k: f64,
}

impl RansFieldTolerances {
    fn values(self) -> [f64; 5] {
        [
            self.velocity_m_s,
            self.kinetic_energy_m2_s2,
            self.dissipation_m2_s3,
            self.eddy_viscosity_m2_s,
            self.temperature_k,
        ]
    }
}

/// All work/time/error limits are explicit. Intervals denote the full trial;
/// each accepted interval consists of two native half steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansAdaptiveControls {
    /// First proposed trial interval [s].
    pub initial_interval_s: f64,
    /// Smallest retry interval [s]. A final endpoint clip may be smaller.
    pub minimum_interval_s: f64,
    /// Largest proposed trial interval [s].
    pub maximum_interval_s: f64,
    /// Includes every accepted, rejected and numerically refused trial.
    pub maximum_attempts: usize,
    /// Common relative local-error allowance; finite and nonnegative.
    pub relative_tolerance: f64,
    /// Positive finite absolute allowances. These never modify state values.
    pub absolute: RansFieldTolerances,
    /// Bound on each field's change relative to its magnitude plus absolute
    /// allowance, in (0, 0.5]. Checked at the half and full trial endpoints.
    /// This prevents a large, stiff step erasing a resolved mode while full
    /// and half-step endpoints both approach zero. Eddy viscosity is derived
    /// from k/epsilon and does not use a potentially stale input nu_t here.
    pub maximum_fractional_change: f64,
}

/// Every failure leaves the caller's whole physical window uncommitted.
#[derive(Clone, Debug, PartialEq)]
pub enum RansAdaptiveError {
    /// Invalid interval, tolerances, bounds or buffer shapes.
    InvalidInput(&'static str),
    /// Original solver refusal, including cancellation.
    Transport(RansTransportError),
    /// The bounded number of trials was exhausted; no partial state publishes.
    AttemptBudget {
        /// Number of trial intervals attempted.
        attempts: usize,
        /// Internally advanced time before the whole window was refused [s].
        elapsed_s: f64,
    },
    /// A rejected interval cannot be refined within the declared lower bound
    /// or represented physical time. Infinity denotes an arithmetic refusal
    /// or an unrepresentable interval.
    MinimumInterval {
        /// Internally advanced time before the whole window was refused [s].
        elapsed_s: f64,
        /// Trial interval that could not be refined [s].
        interval_s: f64,
        /// Normalized rejection indicator, or infinity for arithmetic limits.
        indicator: f64,
    },
}

impl From<RansTransportError> for RansAdaptiveError {
    fn from(error: RansTransportError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for RansAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "RANS adaptive channel: {self:?}")
    }
}
impl std::error::Error for RansAdaptiveError {}

/// Original numerical workspace plus three complete reusable state images.
/// The payload is 24 f64 fields (192 bytes/node); stepping allocates nothing.
pub struct RansAdaptiveWorkspace {
    kernel: RansTransportWorkspace,
    accepted: RansFieldState,
    coarse: RansFieldState,
    fine: RansFieldState,
}

impl RansAdaptiveWorkspace {
    /// Admit a node ceiling before reserving every numerical buffer.
    pub fn try_new(
        grid: &RansChannelGrid,
        maximum_nodes: usize,
    ) -> Result<Self, RansTransportError> {
        let n = grid
            .num_cells
            .checked_add(1)
            .ok_or(RansTransportError::InvalidInput("node count overflow"))?;
        if n > maximum_nodes {
            return Err(RansTransportError::NodeBudget {
                needed: n,
                maximum: maximum_nodes,
            });
        }
        n.checked_mul(24 * core::mem::size_of::<f64>())
            .ok_or(RansTransportError::Allocation)?;
        let field = || {
            let mut values = Vec::new();
            values
                .try_reserve_exact(n)
                .map_err(|_| RansTransportError::Allocation)?;
            values.resize(n, 0.0);
            Ok::<_, RansTransportError>(values)
        };
        let state = || {
            Ok::<_, RansTransportError>(RansFieldState {
                u: field()?,
                k: field()?,
                eps: field()?,
                nu_t: field()?,
                temp_k: field()?,
            })
        };
        Ok(Self {
            kernel: RansTransportWorkspace::try_new(grid, maximum_nodes)?,
            accepted: state()?,
            coarse: state()?,
            fine: state()?,
        })
    }

    /// Numerical payload, excluding fixed metadata and allocator bookkeeping.
    #[must_use]
    pub fn numerical_bytes(&self) -> usize {
        self.accepted.u.len() * 24 * core::mem::size_of::<f64>()
    }
}

/// Heat account over only the accepted fine trajectory, per unit wall area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansAdaptiveHeatBalance {
    /// Original sensible inventory [J/m²].
    pub energy_before_j_m2: f64,
    /// Final sensible inventory [J/m²].
    pub energy_after_j_m2: f64,
    /// Compensated sum of actual local accepted inventory changes [J/m²].
    pub energy_change_j_m2: f64,
    /// Integrated imposed wall input [J/m²].
    pub wall_heat_j_m2: f64,
    /// Integrated volumetric input [J/m²].
    pub volume_heat_j_m2: f64,
    /// Sum of accepted local heat balance defects [J/m²].
    pub balance_defect_j_m2: f64,
}

/// Completed window; refused/coarse trial accounts never contribute here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansAdaptiveReport {
    /// Full physical interval successfully advanced [s].
    pub duration_s: f64,
    /// Number of coarse-versus-two-half-step trials, including refusals.
    pub attempts: usize,
    /// Number of accepted intervals; each contains two native steps.
    pub accepted_intervals: usize,
    /// Error, change-limit or arithmetic rejections.
    pub rejected_intervals: usize,
    /// Actual native step calls, including the failed call in a refused trial.
    pub step_evaluations: usize,
    /// Largest accepted full-versus-half field error ratio, at most one.
    pub maximum_accepted_error_ratio: f64,
    /// Largest accepted field change fraction at either fine endpoint.
    pub maximum_accepted_fractional_change: f64,
    /// Smallest accepted full trial interval [s].
    pub minimum_accepted_interval_s: f64,
    /// Largest accepted full trial interval [s].
    pub maximum_accepted_interval_s: f64,
    /// Initial/final energies and time-averaged rates along the accepted fine
    /// trajectory. `time_step_s` is the full window. Balance defects are
    /// averaged local defects, not a newly applied single backward-Euler step.
    pub mechanical: RansTransportReport,
    /// Present only when passive temperature was advanced.
    pub heat: Option<RansAdaptiveHeatBalance>,
}

#[derive(Clone, Copy)]
struct Step {
    mechanical: RansTransportReport,
    heat: Option<RansThermalTransportReport>,
}

impl RansSolver {
    /// Advance a complete constant-forcing window with bounded local accuracy.
    ///
    /// The original positivity-preserving implicit solver is the sole stepper.
    /// One full trial and two half trials start from the same accepted state;
    /// their maximum componentwise difference is scaled by authored absolute
    /// and relative tolerances. Only the two-half-step trajectory is accepted.
    /// A fractional-change limit also guards unresolved stiff transients.
    /// Interval adaptation uses the first-order step-doubling estimate, safety
    /// 0.9 and factors in [0.2, 2]. Arithmetic refusals retry at one quarter.
    ///
    /// Pressure and optional wall/volume heat inputs are held for this window.
    /// `None` retains temperature exactly and imposes no thermal-material gate.
    /// Every accepted half step refreshes turbulence; optional heat consumes
    /// its own mechanical endpoint. No extra physical model, extrapolation,
    /// clipping or numerical buffer allocation is introduced.
    ///
    /// The entire window publishes together. Cancellation, minimum-interval
    /// refusal, exhausted work or a fatal solver error preserves all caller
    /// fields, even after many internally accepted intervals. Scratch may
    /// change. This is first-order local error control, not a global error,
    /// spatial-resolution, model-validation or intermediate-event guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_launder_sharma_adaptive(
        &self,
        grid: &RansChannelGrid,
        state: &mut RansFieldState,
        pressure_gradient_pa_m: f64,
        duration_s: f64,
        heat: Option<RansHeatLoad<'_>>,
        controls: RansAdaptiveControls,
        work: &mut RansAdaptiveWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RansAdaptiveReport, RansAdaptiveError> {
        checkpoint(0, &mut cancelled)?;
        validate(grid, state, duration_s, controls, work)?;
        copy_state(&mut work.accepted, state, &mut cancelled)?;
        let mut elapsed = 0.0;
        let mut proposal = controls.initial_interval_s;
        let mut attempts = 0;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut evaluations = 0;
        let mut max_error = 0.0_f64;
        let mut max_change = 0.0_f64;
        let mut min_interval = f64::INFINITY;
        let mut max_interval = 0.0_f64;
        let mut account = Account::default();
        while elapsed < duration_s {
            checkpoint(0, &mut cancelled)?;
            if attempts == controls.maximum_attempts {
                return Err(RansAdaptiveError::AttemptBudget {
                    attempts,
                    elapsed_s: elapsed,
                });
            }
            let remaining = duration_s - elapsed;
            let h = proposal.min(remaining);
            let end = if h == remaining {
                duration_s
            } else {
                elapsed + h
            };
            let half = 0.5 * h;
            if !end.is_finite() || end <= elapsed || half <= 0.0 {
                return Err(RansAdaptiveError::MinimumInterval {
                    elapsed_s: elapsed,
                    interval_s: h,
                    indicator: f64::INFINITY,
                });
            }
            attempts += 1;
            let trial = (|| {
                copy_state(&mut work.coarse, &work.accepted, &mut cancelled)?;
                native_step(
                    self,
                    grid,
                    &mut work.coarse,
                    pressure_gradient_pa_m,
                    h,
                    heat,
                    &mut work.kernel,
                    &mut evaluations,
                    &mut cancelled,
                )?;
                copy_state(&mut work.fine, &work.accepted, &mut cancelled)?;
                let first = native_step(
                    self,
                    grid,
                    &mut work.fine,
                    pressure_gradient_pa_m,
                    half,
                    heat,
                    &mut work.kernel,
                    &mut evaluations,
                    &mut cancelled,
                )?;
                let mid_change = change_ratio(
                    &work.accepted,
                    &work.fine,
                    controls,
                    heat.is_some(),
                    &mut cancelled,
                )?;
                let second = native_step(
                    self,
                    grid,
                    &mut work.fine,
                    pressure_gradient_pa_m,
                    half,
                    heat,
                    &mut work.kernel,
                    &mut evaluations,
                    &mut cancelled,
                )?;
                let change = mid_change.max(change_ratio(
                    &work.accepted,
                    &work.fine,
                    controls,
                    heat.is_some(),
                    &mut cancelled,
                )?);
                let error = error_ratio(
                    &work.accepted,
                    &work.coarse,
                    &work.fine,
                    controls,
                    heat.is_some(),
                    &mut cancelled,
                )?;
                Ok::<_, RansTransportError>((first, second, error, change))
            })();
            let (factor, indicator, success) = match trial {
                Ok((first, second, error, change)) => {
                    let error_factor = if error == 0.0 {
                        2.0
                    } else {
                        0.9 / error.sqrt()
                    };
                    let change_factor = if change == 0.0 {
                        2.0
                    } else {
                        0.9 * controls.maximum_fractional_change / change
                    };
                    let success = error <= 1.0 && change <= controls.maximum_fractional_change;
                    if success {
                        account.add(first, half)?;
                        account.add(second, half)?;
                        core::mem::swap(&mut work.accepted, &mut work.fine);
                        elapsed = end;
                        accepted += 1;
                        max_error = max_error.max(error);
                        max_change = max_change.max(change);
                        min_interval = min_interval.min(h);
                        max_interval = max_interval.max(h);
                    }
                    (
                        error_factor.min(change_factor).clamp(0.2, 2.0),
                        error.max(change / controls.maximum_fractional_change),
                        success,
                    )
                }
                Err(RansTransportError::Arithmetic { .. }) => (0.25, f64::INFINITY, false),
                Err(error) => return Err(error.into()),
            };
            if !success {
                rejected += 1;
            }
            if elapsed == duration_s {
                break;
            }
            let next = (h * factor).min(controls.maximum_interval_s);
            if !success
                && (next < controls.minimum_interval_s || next >= h || elapsed + next == elapsed)
            {
                return Err(RansAdaptiveError::MinimumInterval {
                    elapsed_s: elapsed,
                    interval_s: h,
                    indicator,
                });
            }
            proposal = next.max(controls.minimum_interval_s);
        }
        let (mechanical, heat) = account.finish(duration_s)?;
        let report = RansAdaptiveReport {
            duration_s,
            attempts,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            step_evaluations: evaluations,
            maximum_accepted_error_ratio: max_error,
            maximum_accepted_fractional_change: max_change,
            minimum_accepted_interval_s: min_interval,
            maximum_accepted_interval_s: max_interval,
            mechanical,
            heat,
        };
        checkpoint(0, &mut cancelled)?;
        // All fallible work precedes this indivisible publication boundary.
        state.u.copy_from_slice(&work.accepted.u);
        state.k.copy_from_slice(&work.accepted.k);
        state.eps.copy_from_slice(&work.accepted.eps);
        state.nu_t.copy_from_slice(&work.accepted.nu_t);
        state.temp_k.copy_from_slice(&work.accepted.temp_k);
        Ok(report)
    }
}

fn validate(
    grid: &RansChannelGrid,
    state: &RansFieldState,
    duration: f64,
    c: RansAdaptiveControls,
    work: &RansAdaptiveWorkspace,
) -> Result<(), RansAdaptiveError> {
    let n = grid
        .num_cells
        .checked_add(1)
        .ok_or(RansAdaptiveError::InvalidInput("node count overflow"))?;
    if n < 5 || work.accepted.u.len() != n || fields(state).iter().any(|s| s.len() != n) {
        return Err(RansAdaptiveError::InvalidInput(
            "grid/state/workspace shape mismatch",
        ));
    }
    if [
        duration,
        c.initial_interval_s,
        c.minimum_interval_s,
        c.maximum_interval_s,
    ]
    .iter()
    .any(|&x| !x.is_finite() || x <= 0.0)
        || c.minimum_interval_s > c.initial_interval_s
        || c.initial_interval_s > c.maximum_interval_s
        || c.absolute
            .values()
            .iter()
            .any(|&x| !x.is_finite() || x <= 0.0)
        || !c.relative_tolerance.is_finite()
        || c.relative_tolerance < 0.0
        || !c.maximum_fractional_change.is_finite()
        || c.maximum_fractional_change <= 0.0
        || c.maximum_fractional_change > 0.5
        || c.maximum_attempts > usize::MAX / 3
    {
        return Err(RansAdaptiveError::InvalidInput(
            "time, tolerance or work controls",
        ));
    }
    Ok(())
}

fn fields(state: &RansFieldState) -> [&[f64]; 5] {
    [&state.u, &state.k, &state.eps, &state.nu_t, &state.temp_k]
}

fn copy_state(
    to: &mut RansFieldState,
    from: &RansFieldState,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), RansTransportError> {
    for i in (0..from.u.len()).step_by(64) {
        checkpoint(i, cancelled)?;
        let end = (i + 64).min(from.u.len());
        to.u[i..end].copy_from_slice(&from.u[i..end]);
        to.k[i..end].copy_from_slice(&from.k[i..end]);
        to.eps[i..end].copy_from_slice(&from.eps[i..end]);
        to.nu_t[i..end].copy_from_slice(&from.nu_t[i..end]);
        to.temp_k[i..end].copy_from_slice(&from.temp_k[i..end]);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn native_step(
    solver: &RansSolver,
    grid: &RansChannelGrid,
    state: &mut RansFieldState,
    pressure: f64,
    dt: f64,
    heat: Option<RansHeatLoad<'_>>,
    work: &mut RansTransportWorkspace,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Step, RansTransportError> {
    *evaluations += 1;
    if let Some(load) = heat {
        let report = solver
            .advance_launder_sharma_with_heat(grid, state, pressure, dt, load, work, cancelled)?;
        Ok(Step {
            mechanical: report.mechanical,
            heat: Some(report),
        })
    } else {
        let report = solver.advance_launder_sharma(grid, state, pressure, dt, work, cancelled)?;
        Ok(Step {
            mechanical: report,
            heat: None,
        })
    }
}

fn error_ratio(
    before: &RansFieldState,
    coarse: &RansFieldState,
    fine: &RansFieldState,
    c: RansAdaptiveControls,
    heat: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, RansTransportError> {
    let mut result = 0.0_f64;
    let initial = fields(before);
    for (field, ((full, refined), absolute)) in fields(coarse)
        .into_iter()
        .zip(fields(fine))
        .zip(c.absolute.values())
        .enumerate()
    {
        if field == 4 && !heat {
            continue;
        }
        for i in 0..full.len() {
            checkpoint(i, cancelled)?;
            // nu_t is derived, and its incoming cached value need not have
            // been populated by the caller before the first physical step.
            let baseline = if field == 3 {
                full[i]
            } else {
                initial[field][i]
            };
            let scale = absolute + c.relative_tolerance * baseline.abs().max(refined[i].abs());
            if !scale.is_finite() || scale <= 0.0 {
                return Err(RansTransportError::InvalidInput(
                    "unrepresentable adaptive error scale",
                ));
            }
            result = result.max((full[i] - refined[i]).abs() / scale);
        }
    }
    Ok(result)
}

fn change_ratio(
    before: &RansFieldState,
    after: &RansFieldState,
    c: RansAdaptiveControls,
    heat: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, RansTransportError> {
    let mut result = 0.0_f64;
    for (field, ((old, new), absolute)) in fields(before)
        .into_iter()
        .zip(fields(after))
        .zip(c.absolute.values())
        .enumerate()
    {
        if field == 3 || (field == 4 && !heat) {
            continue;
        }
        for i in 0..old.len() {
            checkpoint(i, cancelled)?;
            let scale = absolute + old[i].abs().max(new[i].abs());
            if !scale.is_finite() || scale <= 0.0 {
                return Err(RansTransportError::InvalidInput(
                    "unrepresentable adaptive change scale",
                ));
            }
            result = result.max((new[i] - old[i]).abs() / scale);
        }
    }
    Ok(result)
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let part = value - self.correction;
        let next = self.value + part;
        self.correction = (next - self.value) - part;
        self.value = next;
    }
}

#[derive(Default)]
struct Account {
    first: Option<Step>,
    last: Option<Step>,
    // Eight mechanical rate integrals, then wall/volume heat, change, defect.
    integral: [Sum; 12],
}

impl Account {
    fn add(&mut self, step: Step, dt: f64) -> Result<(), RansTransportError> {
        let m = step.mechanical;
        let values = [
            m.pressure_power,
            m.molecular_dissipation,
            m.turbulence_production,
            m.turbulent_destruction,
            m.turbulent_wall_loss,
            m.temporal_dissipation,
            m.mean_balance_defect,
            m.turbulent_balance_defect,
        ];
        for (sum, rate) in self.integral[..8].iter_mut().zip(values) {
            sum.add(dt * rate);
        }
        if let Some(h) = step.heat {
            self.integral[8].add(dt * h.wall_heat_flux_w_m2);
            self.integral[9].add(dt * h.volumetric_heat_input_w_m2);
            self.integral[10].add(h.thermal_energy_change_j_m2);
            self.integral[11].add(dt * h.thermal_balance_defect_w_m2);
        }
        if self.integral.iter().any(|sum| !sum.value.is_finite()) {
            return Err(RansTransportError::Arithmetic {
                equation: "adaptive integrated balance",
                node: 0,
            });
        }
        self.first.get_or_insert(step);
        self.last = Some(step);
        Ok(())
    }

    fn finish(
        self,
        duration: f64,
    ) -> Result<(RansTransportReport, Option<RansAdaptiveHeatBalance>), RansTransportError> {
        let first = self
            .first
            .ok_or(RansTransportError::InvalidInput("no accepted interval"))?;
        let last = self
            .last
            .ok_or(RansTransportError::InvalidInput("no accepted interval"))?;
        let rates: [f64; 8] = core::array::from_fn(|i| self.integral[i].value / duration);
        if rates.iter().any(|value| !value.is_finite()) {
            return Err(RansTransportError::Arithmetic {
                equation: "adaptive average balance",
                node: 0,
            });
        }
        let mechanical = RansTransportReport {
            time_step_s: duration,
            mean_energy_before: first.mechanical.mean_energy_before,
            mean_energy_after: last.mechanical.mean_energy_after,
            turbulent_energy_before: first.mechanical.turbulent_energy_before,
            turbulent_energy_after: last.mechanical.turbulent_energy_after,
            pressure_power: rates[0],
            molecular_dissipation: rates[1],
            turbulence_production: rates[2],
            turbulent_destruction: rates[3],
            turbulent_wall_loss: rates[4],
            temporal_dissipation: rates[5],
            mean_balance_defect: rates[6],
            turbulent_balance_defect: rates[7],
        };
        let heat = first
            .heat
            .zip(last.heat)
            .map(|(before, after)| RansAdaptiveHeatBalance {
                energy_before_j_m2: before.thermal_energy_before_j_m2,
                energy_after_j_m2: after.thermal_energy_after_j_m2,
                energy_change_j_m2: self.integral[10].value,
                wall_heat_j_m2: self.integral[8].value,
                volume_heat_j_m2: self.integral[9].value,
                balance_defect_j_m2: self.integral[11].value,
            });
        Ok((mechanical, heat))
    }
}

#[cfg(test)]
mod tests;
