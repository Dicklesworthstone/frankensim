//! Nonautonomous homogeneous viscosity on the original conservative Strang map.

use super::*;
use crate::BudgetRefusal;
use crate::diffusion::viscosity_history::{
    PseViscosityError, PseViscosityHistory, ViscosityHistoryError,
};

/// Original coefficient or native evolution refusal; no partial cloud escapes.
#[derive(Debug, Clone, PartialEq)]
pub enum TimeDependentVpmError {
    /// Coefficient support or source failure.
    Coefficient(PseViscosityError),
    /// Original mechanical, resource or cancellation refusal.
    Evolution(VpmError),
}
impl core::fmt::Display for TimeDependentVpmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Coefficient(error) => error.fmt(f),
            Self::Evolution(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for TimeDependentVpmError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Coefficient(error) => Some(error),
            Self::Evolution(error) => Some(error),
        }
    }
}
impl From<PseViscosityError> for TimeDependentVpmError {
    fn from(error: PseViscosityError) -> Self {
        Self::Coefficient(error)
    }
}
impl From<ViscosityHistoryError> for TimeDependentVpmError {
    fn from(error: ViscosityHistoryError) -> Self {
        PseViscosityError::from(error).into()
    }
}
impl From<VpmError> for TimeDependentVpmError {
    fn from(error: VpmError) -> Self {
        Self::Evolution(error)
    }
}
impl From<BudgetRefusal> for TimeDependentVpmError {
    fn from(error: BudgetRefusal) -> Self {
        VpmError::from(error).into()
    }
}

/// Native particle result and the exact completed physical window.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeDependentVpmRun {
    /// Fully accepted native state and original resource receipt.
    pub run: VpmRun,
    /// Exact represented start/end times [s].
    pub time_s: [f64; 2],
}

/// Evolve a moving cloud under an explicitly interpolated homogeneous `nu(t)`.
///
/// Every step ends at the next coefficient knot, the maximum physical step or
/// the requested endpoint. The two original positive PSE half sweeps use their
/// respective exactly integrated mean viscosities. Thus a changing coefficient
/// is resolved in its physical order around conservative advection. Replacing
/// both halves by a full-step average is not this nonautonomous composition.
/// The spatial kernel remains fixed; interpolation is the declared coefficient
/// law, not an inferred material constitutive relation between source samples.
///
/// This is second-order split evolution on each smooth admitted segment, not
/// adaptive error control. Refine time and quadrature/kernel independently. A
/// constant history on the same represented step partition uses the identical
/// native trajectory and pair/work/memory accounting as the original constant
/// owner. The maximum step is a strict represented-time bound; a final tiny
/// remainder can be necessary when decimal endpoints do not share that grid.
///
/// One particle vector is allocated. Planning and stepping both traverse the
/// same borrowed history without a timeline allocation or source callbacks.
/// Each step charges `N*(N-1)` advection pairs and that many additional pairs
/// for each positive-viscosity half; total work is `N + steps + pairs`.
/// Planning is capped by `budget.max_steps` and checks cancellation every 256
/// steps before one original ambient work admission. Publication is atomic for
/// the whole window. The original finite-cloud impulse and boundary limitations
/// of PSE continue to apply; density evolution and spatial viscosity gradients
/// need a separate physical formulation.
///
/// # Errors
/// Refuses unsupported/unrepresentable time or coefficient windows, original
/// numerical inputs, local/ambient budgets, cancellation and allocation errors.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    start_s: f64,
    end_s: f64,
    maximum_step_s: f64,
    core: f64,
    history: PseViscosityHistory<'_>,
    budget: VpmBudget,
) -> Result<TimeDependentVpmRun, TimeDependentVpmError> {
    diffusion::validate_inputs(particles.len(), areas_m2, end_s - start_s)?;
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None).into());
    }
    admit_cap("particles", particles.len(), budget.max_particles)?;
    let profile = history.profile();
    let mut time = start_s;
    let mut steps = 0_usize;
    let mut pair_evaluations = 0_u64;
    while let Some(interval) = profile.next_interval(time, end_s, maximum_step_s)? {
        if steps % 256 == 0 {
            cx.checkpoint().map_err(|_| BudgetRefusal::Cancelled {
                phase: "fs-vpm.viscosity-history-plan",
            })?;
        }
        steps = steps.checked_add(1).ok_or(VpmError::PlanOverflow {
            resource: "viscosity history steps",
        })?;
        admit_cap("steps", steps, budget.max_steps)?;
        let sweeps = 1 + interval
            .half_mean_viscosities_m2_s
            .iter()
            .filter(|&&nu| nu > 0.0)
            .count() as u64;
        pair_evaluations = pair_evaluations
            .checked_add(diffusion::pair_work(particles.len(), 1, sweeps)?)
            .ok_or(VpmError::PlanOverflow {
                resource: "viscosity history pair work",
            })?;
        if pair_evaluations > budget.max_pair_evaluations {
            return Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: pair_evaluations.into(),
                maximum: budget.max_pair_evaluations.into(),
            }
            .into());
        }
        time = interval.time_s[1];
    }
    let work_units = pair_evaluations
        .checked_add(particles.len() as u64)
        .and_then(|work| work.checked_add(steps as u64))
        .ok_or(VpmError::PlanOverflow {
            resource: "viscosity history total work",
        })?;
    let peak_live_bytes = particles
        .len()
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .ok_or(VpmError::PlanOverflow {
            resource: "viscosity history live bytes",
        })?;
    admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
    accountant.checkpoint("fs-vpm.viscosity-history-admission", cx)?;
    let mut state = reserved_vec("viscosity-history-state", particles.len())?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(
            &mut accountant,
            cx,
            "fs-vpm.viscosity-history-input",
            || {
                validate_particle(particle, index)?;
                diffusion::validate_area(areas_m2[index], index)?;
                state.push(particle);
                Ok(())
            },
        )?;
    }
    meter.boundary(
        &mut accountant,
        cx,
        "fs-vpm.viscosity-history-input-complete",
    )?;
    let mut completed_pairs = 0;
    let mut backend = MovingDiffusion::dense();
    time = start_s;
    while let Some(interval) = profile.next_interval(time, end_s, maximum_step_s)? {
        let [left, right] = interval.time_s;
        super::weighted_step_with_laws_observed(
            &mut state,
            areas_m2,
            right - left,
            1.0,
            core,
            history.half_laws(interval)?,
            &mut backend,
            &mut completed_pairs,
            &mut meter,
            &mut accountant,
            cx,
            &mut |_, _, _, _| Ok(()),
        )?;
        time = right;
    }
    meter.flush(&mut accountant, "fs-vpm.viscosity-history-finalize")?;
    let consumption = accountant.consumption();
    if meter.completed != work_units
        || completed_pairs != pair_evaluations
        || consumption.cost_charged != work_units
    {
        return Err(VpmError::AccountingMismatch {
            planned: work_units,
            completed: meter.completed,
            charged: consumption.cost_charged,
        }
        .into());
    }
    accountant.checkpoint("fs-vpm.viscosity-history-publication", cx)?;
    Ok(TimeDependentVpmRun {
        run: VpmRun {
            particles: state,
            steps_completed: steps,
            pair_evaluations,
            work_units,
            peak_live_bytes,
            budget: accountant.consumption(),
        },
        time_s: [start_s, time],
    })
}
