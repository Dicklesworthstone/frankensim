//! Actual-state constitutive drag on a fixed gas amount and prescribed wall.
//!
//! Each elementary map uses the existing exact frozen resistance owner. Its
//! coefficient comes from the original source at a forward predicted midpoint
//! on the exact mechanical energy path. Step doubling controls slip and sensible
//! heating errors; a midpoint ordering guard rejects stiff, collapsing-rate
//! predictors. The error estimates are numerical estimates, not certified bounds.

use super::{
    Cell, MixtureGas, Time, WallCellReport, WallContact, WallControl, WallDrag, WallExchangeError,
    WallUpdate, drag, nonnegative, poll, positive,
};

/// Accuracy and whole-interval work limits for actual-state wall resistance.
///
/// Relative velocity accuracy uses entering wall-relative slip, never a large
/// laboratory boost. Heating accuracy uses physical slip loss, never a species
/// formation-energy datum. These are separate from the physical balance budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceCoupledDragControl {
    /// Positive absolute allowance on accumulated slip-error estimates [m/s].
    pub absolute_slip_tolerance_m_s: f64,
    /// Positive relative allowance against entering slip, at most 1e-2.
    pub relative_slip_tolerance: f64,
    /// Positive absolute allowance on sensible-heating error estimates [J].
    pub absolute_heating_tolerance_j: f64,
    /// Positive relative allowance against actual slip heating, at most 1e-2.
    pub relative_heating_tolerance: f64,
    /// Whole-interval attempted steps, including rejected predictors, in 1..=1e6.
    pub maximum_attempted_steps: usize,
    /// Whole-interval source evaluations, including rejected trials, in 1..=1e7.
    pub maximum_source_evaluations: usize,
}

impl Default for SourceCoupledDragControl {
    fn default() -> Self {
        Self {
            absolute_slip_tolerance_m_s: 1e-9,
            relative_slip_tolerance: 1e-6,
            absolute_heating_tolerance_j: 1e-10,
            relative_heating_tolerance: 1e-6,
            maximum_attempted_steps: 16384,
            maximum_source_evaluations: 131072,
        }
    }
}

impl SourceCoupledDragControl {
    pub(crate) fn validate(self) -> Result<(), WallExchangeError> {
        if !positive(self.absolute_slip_tolerance_m_s)
            || !positive(self.relative_slip_tolerance)
            || self.relative_slip_tolerance > 1e-2
            || !positive(self.absolute_heating_tolerance_j)
            || !positive(self.relative_heating_tolerance)
            || self.relative_heating_tolerance > 1e-2
            || !(1..=1_000_000).contains(&self.maximum_attempted_steps)
            || !(1..=10_000_000).contains(&self.maximum_source_evaluations)
        {
            return Err(WallExchangeError::Invalid(
                "source-coupled wall drag controls",
            ));
        }
        Ok(())
    }
}

/// Work and accumulated embedded estimates for one complete mechanical interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SourceCoupledDragReport {
    /// Attempted full-versus-two-half intervals, including all rejections.
    pub attempted_steps: usize,
    /// Accepted two-half-step endpoints.
    pub accepted_steps: usize,
    /// Original source selections, including rejected trials and final admission.
    pub source_evaluations: usize,
    /// Sum of absolute full-versus-two-half slip defects divided by three [m/s].
    pub estimated_slip_error_m_s: f64,
    /// Sum of absolute physical-heating defects divided by three [J].
    pub estimated_heating_error_j: f64,
}

/// Source-aware mechanical endpoint; thermal exchange remains a separate owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceCoupledDragUpdate {
    /// Actual gas momentum, total energy, wall work, impulse and slip heating.
    pub gas: WallUpdate,
    /// Numerical work and accuracy estimates; these are never deposited as heat.
    pub source: SourceCoupledDragReport,
}

/// Integrate current-state resistance against a stationary or prescribed wall.
///
/// The resolver must select the original constitutive law at the supplied cell,
/// including its current temperature, pressure, composition and relative speed.
/// No coefficient extrapolation or fitted viscosity is supplied here. Its wall
/// speed must remain the prescribed speed in `contact`; a finite inertial support
/// requires the separate coupled momentum owner.
///
/// Species, density and volume stay fixed. Every source query uses the exact
/// invariant `E = E0 + U_wall * (momentum - momentum0)`. The final actual impulse
/// determines work and physical slip heating once. The thermal contact is copied
/// unchanged; callers can pass the mechanical endpoint to the original heat
/// owner with drag disabled. Source refusals and cancellation return no endpoint.
///
/// Exponential midpoint with accepted step doubling is second order for smooth
/// admitted resistance. The explicit ordering guard prevents a midpoint which
/// has already relaxed past its proposed endpoint from hiding a stiff rate drop.
/// It does not certify arbitrary discontinuous constitutive sources or promote
/// the subsequent thermal splitting to second order.
#[allow(clippy::too_many_arguments)]
pub fn exchange_cell_with_source_drag(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    control: WallControl,
    source_control: SourceCoupledDragControl,
    mut resolve: impl FnMut(Cell<'_>, &mut dyn FnMut() -> bool) -> Result<WallDrag, WallExchangeError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SourceCoupledDragUpdate, WallExchangeError> {
    poll(&mut cancelled)?;
    contact.validate()?;
    control.validate()?;
    source_control.validate()?;
    let dt = duration.value();
    if !positive(volume_m3) || !nonnegative(dt) {
        return Err(WallExchangeError::Invalid(
            "source-coupled wall volume or interval",
        ));
    }
    let primitive = gas.primitive(state, &mut cancelled)?;
    let wall = contact.drag.moving_velocity().unwrap_or(0.0);
    let initial_slip = primitive.velocity - wall;
    if !initial_slip.is_finite() {
        return Err(WallExchangeError::Invalid(
            "source-coupled wall relative velocity",
        ));
    }
    let path = Path {
        gas,
        state,
        density: primitive.density,
        wall,
        volume: volume_m3,
    };
    let mut report = SourceCoupledDragReport::default();
    let mut current = Stage {
        momentum: state.momentum,
        law: evaluate(
            &path,
            state.momentum,
            source_control,
            &mut report,
            &mut resolve,
            &mut cancelled,
        )?,
    };
    let slip_allowance = source_control.absolute_slip_tolerance_m_s
        + source_control.relative_slip_tolerance * initial_slip.abs();
    if !positive(slip_allowance) {
        return Err(WallExchangeError::Invalid(
            "unrepresentable wall slip allowance",
        ));
    }
    let mut elapsed = 0.0;
    let mut step = dt;
    while elapsed < dt && initial_slip != 0.0 {
        poll(&mut cancelled)?;
        if report.attempted_steps >= source_control.maximum_attempted_steps {
            return Err(work_error(report));
        }
        report.attempted_steps += 1;
        let remaining = dt - elapsed;
        step = step.min(remaining);
        if !positive(step) || elapsed + step == elapsed || 0.5 * step == 0.0 {
            return Err(WallExchangeError::Invalid(
                "unrepresentable source-coupled wall step",
            ));
        }
        let Some(coarse) = midpoint(
            &path,
            current,
            step,
            source_control,
            &mut report,
            &mut resolve,
            &mut cancelled,
        )?
        else {
            step *= 0.5;
            continue;
        };
        let Some(half) = midpoint(
            &path,
            current,
            0.5 * step,
            source_control,
            &mut report,
            &mut resolve,
            &mut cancelled,
        )?
        else {
            step *= 0.5;
            continue;
        };
        let Some(fine) = midpoint(
            &path,
            half,
            0.5 * step,
            source_control,
            &mut report,
            &mut resolve,
            &mut cancelled,
        )?
        else {
            step *= 0.5;
            continue;
        };
        let defect = fine.momentum - coarse.momentum;
        let slip_error = (defect / path.density).abs() / 3.0;
        // Difference of two slip-energy losses, factored so weak heat changes
        // are not erased by subtracting almost equal full kinetic energies.
        let mean_slip = path.slip(current.momentum)
            + 0.5
                * ((fine.momentum - current.momentum) / path.density
                    + (coarse.momentum - current.momentum) / path.density);
        let heat_error = (defect * mean_slip * volume_m3).abs() / 3.0;
        let heating = path.heating(current.momentum, fine.momentum);
        let fraction = step / dt;
        let slip_progress =
            ((fine.momentum - current.momentum) / path.density).abs() / initial_slip.abs();
        // The energy fraction is evaluated as a ratio of slip changes, avoiding
        // construction of a potentially overflowing initial extensive energy.
        let heat_progress = slip_progress
            * (path.slip(current.momentum).abs() / initial_slip.abs()
                + path.slip(fine.momentum).abs() / initial_slip.abs());
        // Half the allowance follows time and half follows monotone physical
        // relaxation. Both fractions sum to at most one over the whole interval,
        // so a rapid initial transient can be resolved without resetting budgets
        // or receiving a vanishing share merely because its time scale is short.
        // Reserve ten percent for accumulation and representation.
        let slip_budget = 0.45 * slip_allowance * (fraction + slip_progress);
        let heat_budget = 0.9
            * (0.5 * source_control.absolute_heating_tolerance_j * (fraction + heat_progress)
                + source_control.relative_heating_tolerance * heating);
        if !nonnegative(slip_error)
            || !nonnegative(heat_error)
            || !nonnegative(heating)
            || !positive(slip_budget)
            || !positive(heat_budget)
        {
            return Err(WallExchangeError::Invalid(
                "source-coupled wall accuracy scales",
            ));
        }
        let ratio = (slip_error / slip_budget).max(heat_error / heat_budget);
        if ratio > 1.0 {
            step *= 0.5;
            continue;
        }
        current = fine;
        report.accepted_steps += 1;
        report.estimated_slip_error_m_s += slip_error;
        report.estimated_heating_error_j += heat_error;
        elapsed = if step == remaining {
            dt
        } else {
            elapsed + step
        };
        // Dyadic adaptation is deterministic and needs no extra transcendental
        // evaluations. Cubic local defects permit doubling below one eighth.
        if ratio <= 0.125 {
            step = (2.0 * step).min(dt - elapsed);
        }
    }
    let endpoint = path.cell(current.momentum)?;
    gas.primitive(endpoint, &mut cancelled)?;
    // Admit the actual returned state again, including identity/zero-time cases.
    current.law = evaluate(
        &path,
        current.momentum,
        source_control,
        &mut report,
        &mut resolve,
        &mut cancelled,
    )?;
    let change = current.momentum - state.momentum;
    let work_density = wall * change;
    let work = work_density * volume_m3;
    let impulse = -change * volume_m3;
    let heating = path.heating(state.momentum, current.momentum);
    let residual = ((endpoint.energy - state.energy) - work_density) * volume_m3;
    let heat_allowance = source_control.absolute_heating_tolerance_j
        + source_control.relative_heating_tolerance * heating;
    if ![
        current.momentum,
        endpoint.energy,
        work,
        impulse,
        heating,
        residual,
        heat_allowance,
        report.estimated_slip_error_m_s,
        report.estimated_heating_error_j,
    ]
    .into_iter()
    .all(f64::is_finite)
        || heating < 0.0
    {
        return Err(WallExchangeError::Invalid(
            "source-coupled wall mechanical balance",
        ));
    }
    if report.estimated_slip_error_m_s > slip_allowance
        || report.estimated_heating_error_j > heat_allowance
    {
        return Err(work_error(report));
    }
    if residual.abs() > control.tolerance(work.abs()) {
        return Err(WallExchangeError::Balance {
            energy_residual_j: residual,
            thermal_residual_j: 0.0,
        });
    }
    poll(&mut cancelled)?;
    Ok(SourceCoupledDragUpdate {
        gas: WallUpdate {
            momentum: current.momentum,
            energy: endpoint.energy,
            contact: WallContact {
                drag: current.law,
                ..contact
            },
            report: WallCellReport {
                elapsed: duration,
                wall_impulse_ns: impulse,
                mechanical_work_to_gas_j: work,
                drag_heating_j: heating,
                heat_to_gas_j: 0.0,
                wall_energy_change_j: 0.0,
                external_heat_to_gas_j: 0.0,
                energy_residual_j: residual,
                thermal_residual_j: 0.0,
                iterations: 0,
            },
        },
        source: report,
    })
}

struct Path<'a> {
    gas: &'a MixtureGas,
    state: Cell<'a>,
    density: f64,
    wall: f64,
    volume: f64,
}

impl Path<'_> {
    fn cell(&self, momentum: f64) -> Result<Cell<'_>, WallExchangeError> {
        // Anchor every trial to the entering actual impulse; repeated rounded
        // work additions must not drift the source's temperature path.
        let energy = self.state.energy + self.wall * (momentum - self.state.momentum);
        if !momentum.is_finite() || !energy.is_finite() {
            return Err(WallExchangeError::Invalid("source-coupled wall state"));
        }
        Ok(Cell {
            momentum,
            energy,
            ..self.state
        })
    }

    fn slip(&self, momentum: f64) -> f64 {
        momentum / self.density - self.wall
    }

    fn heating(&self, initial: f64, final_momentum: f64) -> f64 {
        let change = final_momentum - initial;
        -change * (self.slip(initial) + 0.5 * (change / self.density)) * self.volume
    }

    fn frozen(&self, state: Stage, dt: f64) -> Result<f64, WallExchangeError> {
        let cell = self.cell(state.momentum)?;
        // Original exact exponential/Darcy/Riccati mechanics is the only
        // frozen resistance integrator, including weak and stiff limits.
        Ok(drag::advance(
            state.law,
            cell.momentum,
            cell.energy,
            self.density,
            cell.momentum / self.density,
            self.volume,
            dt,
        )?
        .momentum)
    }
}

#[derive(Clone, Copy)]
struct Stage {
    momentum: f64,
    law: WallDrag,
}

#[allow(clippy::too_many_arguments)]
fn midpoint(
    path: &Path<'_>,
    entering: Stage,
    dt: f64,
    control: SourceCoupledDragControl,
    report: &mut SourceCoupledDragReport,
    resolve: &mut impl FnMut(Cell<'_>, &mut dyn FnMut() -> bool) -> Result<WallDrag, WallExchangeError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Stage>, WallExchangeError> {
    let middle = path.frozen(entering, 0.5 * dt)?;
    let middle_law = evaluate(path, middle, control, report, resolve, cancelled)?;
    let momentum = path.frozen(
        Stage {
            law: middle_law,
            ..entering
        },
        dt,
    )?;
    let start_slip = path.slip(entering.momentum).abs();
    let middle_slip = path.slip(middle).abs();
    let end_slip = path.slip(momentum).abs();
    let roundoff = 16.0 * f64::EPSILON * start_slip;
    if middle_slip > start_slip + roundoff || end_slip > middle_slip + roundoff {
        return Ok(None);
    }
    Ok(Some(Stage {
        momentum,
        law: evaluate(path, momentum, control, report, resolve, cancelled)?,
    }))
}

fn evaluate(
    path: &Path<'_>,
    momentum: f64,
    control: SourceCoupledDragControl,
    report: &mut SourceCoupledDragReport,
    resolve: &mut impl FnMut(Cell<'_>, &mut dyn FnMut() -> bool) -> Result<WallDrag, WallExchangeError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WallDrag, WallExchangeError> {
    poll(cancelled)?;
    if report.source_evaluations >= control.maximum_source_evaluations {
        return Err(work_error(*report));
    }
    report.source_evaluations += 1;
    let state = path.cell(momentum)?;
    // The public callback can be a constant law; original caloric admission is
    // still required for every constitutive stage, not only the final endpoint.
    path.gas.primitive(state, &mut *cancelled)?;
    let law = resolve(state, cancelled)?;
    law.validate()?;
    if law.moving_velocity().unwrap_or(0.0) != path.wall {
        return Err(WallExchangeError::Invalid(
            "source changed prescribed wall velocity",
        ));
    }
    poll(cancelled)?;
    Ok(law)
}

fn work_error(report: SourceCoupledDragReport) -> WallExchangeError {
    WallExchangeError::SourceDragWork {
        attempted_steps: report.attempted_steps,
        source_evaluations: report.source_evaluations,
    }
}

#[cfg(test)]
mod tests;
