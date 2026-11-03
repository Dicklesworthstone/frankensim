//! Conservative viscous evolution with explicitly scheduled regular remeshing.
//!
//! Original D/2-A-D/2 steps and the original M4' remap share one execution
//! account. Target grids, padding, resolution and cadence are authored data.

use crate::{
    AdmittedBudget, BudgetConsumption, BudgetRefusal, Cx, VortexParticle, VpmError, WorkMeter,
    admit_cap,
    diffusion::{self, PseDiffusion},
    invalid_scalar,
    remesh::{self, CartesianGrid2, RemeshBudget, RemeshError, RemeshedCloud2},
    reserved_vec,
};

/// Initial grid followed by each periodic/final endpoint grid.
///
/// There must be exactly `1 + ceil(steps/every_steps)` grids. The first remap
/// occurs before evolution; later ones follow each complete interval and the
/// final partial interval. Zero steps performs only the initial remap. Zero
/// time still executes the authored remaps. Each target must contain every
/// nonzero M4' stencil. Widening padding is usually needed once PSE fills the
/// old boundary; the driver does not discard those physical tails.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemeshSchedule<'a> {
    /// A positive number of original evolution steps between remaps.
    pub every_steps: usize,
    /// Explicit grids, including the initial and final grids.
    pub grids: &'a [CartesianGrid2],
}

/// Combined admission envelope for the entire authored trajectory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RemeshedVpmBudget {
    /// Maximum borrowed particles before the initial remap.
    pub max_input_particles: usize,
    /// Maximum nodes in any authored grid, including zero-strength support.
    pub max_grid_particles: usize,
    /// Maximum total original D/2-A-D/2 steps.
    pub max_steps: usize,
    /// Maximum total original pair visits over all changing grid sizes.
    pub max_pair_evaluations: u64,
    /// Maximum pair, step-boundary and remap visits for the complete run.
    pub max_work_units: u64,
    /// Maximum owned payload reserved for both alternating buffer banks.
    /// Borrowed inputs/grids and allocator metadata/rounding are excluded.
    pub max_live_bytes: usize,
}

/// Final regular cloud and the complete execution account.
#[derive(Debug, Clone, PartialEq)]
pub struct RemeshedVpmRun {
    /// Final row-major grid particles, including zero-strength nodes.
    pub particles: Vec<VortexParticle>,
    /// Actual final-grid quadrature areas, independent of the velocity core.
    pub areas_m2: Vec<f64>,
    /// The authored final target grid.
    pub grid: CartesianGrid2,
    /// Completed physical integration steps.
    pub steps_completed: usize,
    /// Completed remaps, including the initial and final ones.
    pub remeshes_completed: usize,
    /// All original advection and PSE pair visits.
    pub pair_evaluations: u64,
    /// All M4' stencil visits, including zero-weight visits.
    pub stencil_deposits: u64,
    /// Exact total successful kernel work charged to the one ambient account.
    pub work_units: u64,
    /// Peak logical reserved payload across all simultaneously owned banks.
    pub peak_live_bytes: usize,
    /// Original ambient budget and exact complete-run consumption.
    pub budget: BudgetConsumption,
}

/// Original evolution and remap refusals, with no partial trajectory published.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RemeshedVpmError {
    /// Original physical input, planning, allocation or execution refusal.
    Evolution(VpmError),
    /// Original stencil support, geometry, remap or schedule refusal.
    Remesh(RemeshError),
}
impl From<VpmError> for RemeshedVpmError {
    fn from(error: VpmError) -> Self {
        Self::Evolution(error)
    }
}
impl From<RemeshError> for RemeshedVpmError {
    fn from(error: RemeshError) -> Self {
        Self::Remesh(error)
    }
}
impl From<BudgetRefusal> for RemeshedVpmError {
    fn from(error: BudgetRefusal) -> Self {
        Self::Evolution(error.into())
    }
}
impl core::fmt::Display for RemeshedVpmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Evolution(error) => error.fmt(f),
            Self::Remesh(error) => error.fmt(f),
        }
    }
}
impl core::error::Error for RemeshedVpmError {}

/// Evolve one immutable circulation inventory through the authored remaps.
///
/// The initial remap supplies regular quadrature before any diffusion. Each
/// intervening step calls the unchanged conservative viscous step owner. Each
/// later remap preserves circulation and first/quadratic vortex moments while
/// restoring the declared grid; the velocity core and viscosity stay separate.
///
/// One budget admits all stages before allocating. At most two banks of
/// particles plus f64 scratch/areas are reserved to the largest target grid,
/// reused thereafter without per-step copies or later allocation. Remapping
/// and pair work poll the same execution account within 256 logical visits.
/// Cancellation or failure at any late stage returns no trajectory.
///
/// This inherits finite-cloud PSE spatial/padding error and the signed M4'
/// interpolation error. Remapping does not make diffusion impulse-exact or
/// positive. Cadence, target-grid/kernel resolution and timestep need separate
/// refinement; an arbitrary changing remesh schedule is not a second-order
/// temporal convergence claim. Walls, automatic adaptivity and 3-D stretching
/// are not inferred.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    law: PseDiffusion,
    schedule: RemeshSchedule<'_>,
    budget: RemeshedVpmBudget,
) -> Result<RemeshedVpmRun, RemeshedVpmError> {
    if !dt.is_finite() || dt < 0.0 {
        return Err(invalid_scalar("dt", dt, None).into());
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None).into());
    }
    let plan = Plan::new(cx, particles.len(), dt, steps, law, schedule, budget)?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, plan.work)?;
    accountant.checkpoint("fs-vpm.remeshed-admission", cx)?;
    let first = workspace(plan.max_grid)?;
    let mut spare = if schedule.grids.len() > 1 {
        Some(workspace(plan.max_grid)?)
    } else {
        None
    };
    let local = RemeshBudget {
        max_input_particles: particles.len().max(plan.max_grid),
        max_grid_particles: plan.max_grid,
        max_work_units: plan.remesh_work,
        max_live_bytes: plan.bank_bytes,
    };
    let mut cloud = remap(
        cx,
        particles,
        schedule.grids[0],
        local,
        first,
        &mut accountant,
    )?;
    let mut remesh_work = cloud.work_units;
    let mut deposits = cloud.stencil_deposits;
    let mut remeshes = 1;
    let mut meter = WorkMeter::default();
    let mut pairs = 0;
    for step in 0..steps {
        super::viscous::step(
            &mut cloud.particles,
            &cloud.areas_m2,
            dt,
            core,
            law,
            &mut pairs,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        if (step + 1) % schedule.every_steps == 0 || step + 1 == steps {
            // The original step flushes its meter before remap work enters the
            // same account. The old bank becomes the next remap's workspace.
            let next = remap(
                cx,
                &cloud.particles,
                schedule.grids[remeshes],
                local,
                spare.take().expect("admitted second bank"),
                &mut accountant,
            )?;
            remesh_work += next.work_units;
            deposits += next.stencil_deposits;
            let previous = core::mem::replace(&mut cloud, next);
            spare = Some((previous.particles, previous.areas_m2));
            remeshes += 1;
        }
    }
    meter.flush(&mut accountant, "fs-vpm.remeshed-finalize")?;
    let completed = meter.completed + remesh_work;
    let consumption = accountant.consumption();
    if completed != plan.work
        || pairs != plan.pairs
        || remesh_work != plan.remesh_work
        || deposits != plan.deposits
        || consumption.cost_charged != plan.work
    {
        return Err(VpmError::AccountingMismatch {
            planned: plan.work,
            completed,
            charged: consumption.cost_charged,
        }
        .into());
    }
    accountant.checkpoint("fs-vpm.remeshed-publication", cx)?;
    Ok(RemeshedVpmRun {
        particles: cloud.particles,
        areas_m2: cloud.areas_m2,
        grid: cloud.grid,
        steps_completed: steps,
        remeshes_completed: remeshes,
        pair_evaluations: pairs,
        stencil_deposits: deposits,
        work_units: completed,
        peak_live_bytes: plan.peak_bytes,
        budget: accountant.consumption(),
    })
}

type Workspace = (Vec<VortexParticle>, Vec<f64>);

fn workspace(nodes: usize) -> Result<Workspace, VpmError> {
    Ok((
        reserved_vec("remeshed-particles", nodes)?,
        reserved_vec("remeshed-areas", nodes)?,
    ))
}

fn remap(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    grid: CartesianGrid2,
    local: RemeshBudget,
    workspace: Workspace,
    accountant: &mut AdmittedBudget<'_>,
) -> Result<RemeshedCloud2, RemeshedVpmError> {
    let mut charged = 0;
    remesh::remesh_with_progress(particles, grid, local, Some(workspace), |completed| {
        accountant.charge_cost("fs-vpm.remeshed-remap", completed - charged)?;
        charged = completed;
        accountant.checkpoint("fs-vpm.remeshed-remap", cx)?;
        Ok(())
    })
}

struct Plan {
    max_grid: usize,
    pairs: u64,
    remesh_work: u64,
    deposits: u64,
    work: u64,
    bank_bytes: usize,
    peak_bytes: usize,
}

impl Plan {
    fn new(
        cx: &Cx<'_>,
        inputs: usize,
        dt: f64,
        steps: usize,
        law: PseDiffusion,
        schedule: RemeshSchedule<'_>,
        budget: RemeshedVpmBudget,
    ) -> Result<Self, RemeshedVpmError> {
        if schedule.every_steps == 0 {
            return Err(RemeshError::Invalid("positive remesh cadence required").into());
        }
        admit_cap("input particles", inputs, budget.max_input_particles)?;
        admit_cap("steps", steps, budget.max_steps)?;
        let count =
            steps
                .div_ceil(schedule.every_steps)
                .checked_add(1)
                .ok_or(VpmError::PlanOverflow {
                    resource: "scheduled remaps",
                })?;
        if schedule.grids.len() != count {
            return Err(VpmError::LengthMismatch {
                field: "remesh grids",
                expected: count,
                actual: schedule.grids.len(),
            }
            .into());
        }
        let sweeps = if law.kinematic_viscosity_m2_s() > 0.0 {
            3
        } else {
            1
        };
        let (mut previous, mut max_grid, mut pairs) = (inputs, 0, 0_u64);
        let (mut remesh_work, mut deposits) = (0_u128, 0_u128);
        let mut remaining = steps;
        for (index, grid) in schedule.grids.iter().enumerate() {
            // Bounded cancellation while inspecting the already-authored plan;
            // the complete ambient cost/deadline admission follows this pass.
            if index % 256 == 0 {
                cx.checkpoint().map_err(|_| BudgetRefusal::Cancelled {
                    phase: "fs-vpm.remeshed-plan",
                })?;
            }
            let nodes = grid.node_count();
            admit_cap("grid particles", nodes, budget.max_grid_particles)?;
            max_grid = max_grid.max(nodes);
            remesh_work += 17 * previous as u128 + 2 * nodes as u128;
            deposits += 16 * previous as u128;
            limit("remap work", remesh_work, budget.max_work_units as u128)?;
            if index + 1 < count {
                let interval = remaining.min(schedule.every_steps);
                if dt > 0.0 {
                    pairs = pairs
                        .checked_add(diffusion::pair_work(nodes, interval, sweeps)?)
                        .ok_or(VpmError::PlanOverflow {
                            resource: "scheduled pair visits",
                        })?;
                }
                remaining -= interval;
            }
            previous = nodes;
        }
        limit(
            "pair evaluations",
            pairs as u128,
            budget.max_pair_evaluations as u128,
        )?;
        let work = remesh_work + pairs as u128 + steps as u128;
        limit("work units", work, budget.max_work_units as u128)?;
        let bank_bytes = max_grid as u128
            * (core::mem::size_of::<VortexParticle>() + core::mem::size_of::<f64>()) as u128;
        let peak_bytes = bank_bytes * if count > 1 { 2 } else { 1 };
        limit(
            "logical live bytes",
            peak_bytes,
            budget.max_live_bytes as u128,
        )?;
        Ok(Self {
            max_grid,
            pairs,
            remesh_work: remesh_work as u64,
            deposits: deposits as u64,
            work: work as u64,
            bank_bytes: bank_bytes as usize,
            peak_bytes: peak_bytes as usize,
        })
    }
}

fn limit(resource: &'static str, required: u128, maximum: u128) -> Result<(), VpmError> {
    if required > maximum {
        Err(VpmError::CapExceeded {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
