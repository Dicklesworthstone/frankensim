//! Conservative evolution with explicitly scheduled regular remeshing.
//!
//! Original viscous steps or reversible inviscid compositions, prescribed
//! affine flows and the original M4' remap share one execution account. Target
//! grids, padding, resolution, cadence and held flows are authored data.

use super::ambient::{AffineFlow2, AmbientMomentChange, HalfFlow};
use super::{CompositionOrder, composition::Composition};
use fs_math::dd::Dd;

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
    /// Maximum total physical integration steps.
    pub max_steps: usize,
    /// Maximum total original pair visits over all changing grid sizes.
    pub max_pair_evaluations: u64,
    /// Maximum pair, step-boundary, source and remap visits for the complete run.
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

/// A remeshed trajectory with separately measured prescribed-flow effects.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRemeshedVpmRun {
    /// Original final cloud, quadrature and complete shared execution account.
    pub run: RemeshedVpmRun,
    /// Actual signed vortex-moment changes during the source substeps alone.
    /// PSE/remap errors remain separate; these are not fluid/body momentum or
    /// pressure-work receipts for an infinite prescribed background.
    pub ambient_moment_change: AmbientMomentChange,
    /// Two visits per current grid node and active-source evolution step,
    /// or four with the optional fourth-order inviscid composition.
    /// Included in `run.work_units`, independent of pair/remap visits.
    pub source_particle_updates: u64,
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
    run(
        cx,
        particles,
        dt,
        steps,
        core,
        Evolution::Viscous(law),
        schedule,
        None,
        budget,
    )
    .map(|result| result.run)
}

/// Evolve the same scheduled trajectory in piecewise held affine backgrounds.
///
/// `flows` contains exactly `ceil(steps/every_steps)` entries. Flow `k` acts
/// after remapping to `grids[k]` and before remapping to `grids[k+1]`, including
/// the final partial interval. The source changes at that declared boundary;
/// no interpolation or smooth time-dependence is inferred. Zero steps needs
/// an empty source list and performs the original initial remap only.
///
/// Each step is `F(dt/2) D(dt/2) A(dt) D(dt/2) F(dt/2)`, followed by remapping
/// only at the authored endpoints. All stages reuse their original owners.
/// The affine map is prepared once on entering each active interval, with
/// constant stack storage; both cloud banks remain reserved once to the
/// largest grid. There are no per-step copies or further heap allocations.
///
/// Active intervals add exactly `2 * grid_nodes * interval_steps` visits to
/// the original whole-run work plan; zero source or zero time adds none.
/// Source moment increments are accumulated across all grid/flow changes in
/// double-double arithmetic. The final grid has its own original quadrature,
/// while every source subflow preserves its incoming particle areas.
///
/// All target support, input, physical, execution and publication checks remain
/// authoritative. Late map/support failure or cancellation returns no result.
/// This combines the existing source, PSE and remapping error models; arbitrary
/// remap/flow schedules do not inherit a global second-order convergence claim.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_affine_flows_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    law: PseDiffusion,
    schedule: RemeshSchedule<'_>,
    flows: &[AffineFlow2],
    budget: RemeshedVpmBudget,
) -> Result<AmbientRemeshedVpmRun, RemeshedVpmError> {
    run(
        cx,
        particles,
        dt,
        steps,
        core,
        Evolution::Viscous(law),
        schedule,
        Some(flows),
        budget,
    )
}

/// Evolve an inviscid cloud with the selected order between authored remaps.
///
/// `Second` recovers the original zero-viscosity trajectory and execution
/// receipt. `Fourth` uses three signed symmetric exact-pair sweeps per positive
/// physical step, retaining the same vortex-impulse invariants without another
/// cloud buffer. The original M4' remap runs only after complete positive-time
/// intervals, never inside a signed composition. There is no viscosity
/// parameter: irreversible diffusion cannot enter a negative internal stage.
///
/// Time, support, memory, cancellation and atomic-publication rules are the
/// same as [`simulate_with_cx`]. The fourth-order claim concerns inviscid time
/// integration between fixed remap endpoints; changing the remap cadence or
/// grid introduces its separate spatial error and is not a global fourth-order
/// convergence claim. Zero time still executes the authored remaps.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_order_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    order: CompositionOrder,
    schedule: RemeshSchedule<'_>,
    budget: RemeshedVpmBudget,
) -> Result<RemeshedVpmRun, RemeshedVpmError> {
    run(
        cx,
        particles,
        dt,
        steps,
        core,
        Evolution::Inviscid(order),
        schedule,
        None,
        budget,
    )
    .map(|result| result.run)
}

/// Evolve selected-order inviscid steps in the scheduled held affine fields.
///
/// The flow list and remap boundaries follow [`simulate_with_affine_flows_cx`].
/// `Second` is its exact zero-viscosity path. `Fourth` composes the existing
/// reversible pair and affine maps with signed internal times, merging adjacent
/// maps of the same held field into four source passes per physical step.
/// Weighted maps are prepared once per active interval with constant stack
/// storage. All three pair sweeps finish before any remap or source change.
///
/// For each positive-time interval with N grid nodes and S steps, `Fourth`
/// charges `3*S*N*(N-1)` pair visits and, if its source is nonzero, `4*S*N`
/// source visits. The one whole-run account and original two-bank payload bound
/// include every interval and remap. Zero sources add no work or moment change.
/// Source moment changes remain separate from remap error. The same temporal
/// versus spatial convergence qualification as [`simulate_with_order_cx`]
/// applies, including at discontinuities between authored held fields.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_affine_flows_order_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    order: CompositionOrder,
    schedule: RemeshSchedule<'_>,
    flows: &[AffineFlow2],
    budget: RemeshedVpmBudget,
) -> Result<AmbientRemeshedVpmRun, RemeshedVpmError> {
    run(
        cx,
        particles,
        dt,
        steps,
        core,
        Evolution::Inviscid(order),
        schedule,
        Some(flows),
        budget,
    )
}

#[derive(Clone, Copy)]
enum Evolution {
    Viscous(PseDiffusion),
    Inviscid(CompositionOrder),
}

impl Evolution {
    fn order(self) -> CompositionOrder {
        match self {
            Self::Viscous(_) => CompositionOrder::Second,
            Self::Inviscid(order) => order,
        }
    }

    fn pair_sweeps(self) -> u64 {
        match self {
            Self::Viscous(law) if law.kinematic_viscosity_m2_s() > 0.0 => 3,
            _ => self.order().sweeps() as u64,
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    evolution: Evolution,
    schedule: RemeshSchedule<'_>,
    flows: Option<&[AffineFlow2]>,
    budget: RemeshedVpmBudget,
) -> Result<AmbientRemeshedVpmRun, RemeshedVpmError> {
    if !dt.is_finite() || dt < 0.0 {
        return Err(invalid_scalar("dt", dt, None).into());
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None).into());
    }
    let plan = Plan::new(
        cx,
        particles.len(),
        dt,
        steps,
        evolution,
        schedule,
        flows,
        budget,
    )?;
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
    let composition = Composition::new(evolution.order());
    let mut source = None;
    let mut source_updates = 0;
    let mut source_moments = [Dd::ZERO; 3];
    for step in 0..steps {
        if step % schedule.every_steps == 0 {
            source = match flows.map(|flows| flows[step / schedule.every_steps]) {
                Some(flow) if dt > 0.0 && !flow.is_zero() => Some((
                    match evolution {
                        Evolution::Viscous(_) => HalfFlow::new(flow, dt)?,
                        Evolution::Inviscid(_) => {
                            HalfFlow::weighted(flow, dt, composition.weights()[0])?
                        }
                    },
                    composition
                        .bridge_weight()
                        .map(|weight| HalfFlow::weighted(flow, dt, weight))
                        .transpose()?,
                )),
                _ => None,
            };
        }
        if let Some((half, _)) = &source {
            half.apply(
                &mut cloud.particles,
                &mut source_moments,
                &mut meter,
                &mut accountant,
                cx,
            )?;
            source_updates += cloud.particles.len() as u64;
        }
        match evolution {
            Evolution::Viscous(law) => super::viscous::step(
                &mut cloud.particles,
                &cloud.areas_m2,
                dt,
                core,
                law,
                &mut pairs,
                &mut meter,
                &mut accountant,
                cx,
            )?,
            Evolution::Inviscid(_) => {
                if dt > 0.0 {
                    for (index, &weight) in composition.weights().iter().enumerate() {
                        if index > 0
                            && let Some((_, Some(bridge))) = &source
                        {
                            bridge.apply(
                                &mut cloud.particles,
                                &mut source_moments,
                                &mut meter,
                                &mut accountant,
                                cx,
                            )?;
                            source_updates += cloud.particles.len() as u64;
                        }
                        super::weighted_symmetric_sweep(
                            &mut cloud.particles,
                            dt,
                            core,
                            weight,
                            &mut pairs,
                            &mut meter,
                            &mut accountant,
                            cx,
                        )?;
                    }
                }
                meter.unit(&mut accountant, cx, "fs-vpm.conservative-step", || Ok(()))?;
                meter.boundary(&mut accountant, cx, "fs-vpm.conservative-step-complete")?;
            }
        }
        if let Some((half, _)) = &source {
            half.apply(
                &mut cloud.particles,
                &mut source_moments,
                &mut meter,
                &mut accountant,
                cx,
            )?;
            source_updates += cloud.particles.len() as u64;
            // The second source half follows the original step's checkpoint.
            // Flush it before remap progress charges the same ambient account.
            meter.boundary(&mut accountant, cx, "fs-vpm.remeshed-source-complete")?;
        }
        if (step + 1) % schedule.every_steps == 0 || step + 1 == steps {
            // The old bank becomes the next remap's workspace.
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
        || source_updates != plan.source_updates
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
    let run = RemeshedVpmRun {
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
    };
    Ok(AmbientRemeshedVpmRun {
        run,
        ambient_moment_change: AmbientMomentChange {
            first_moment_m3_s: [source_moments[0].to_f64(), source_moments[1].to_f64()],
            quadratic_moment_m4_s: source_moments[2].to_f64(),
        },
        source_particle_updates: source_updates,
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
    source_updates: u64,
    remesh_work: u64,
    deposits: u64,
    work: u64,
    bank_bytes: usize,
    peak_bytes: usize,
}

impl Plan {
    #[allow(clippy::too_many_arguments)]
    fn new(
        cx: &Cx<'_>,
        inputs: usize,
        dt: f64,
        steps: usize,
        evolution: Evolution,
        schedule: RemeshSchedule<'_>,
        flows: Option<&[AffineFlow2]>,
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
        if let Some(flows) = flows
            && flows.len() != count - 1
        {
            return Err(VpmError::LengthMismatch {
                field: "interval affine flows",
                expected: count - 1,
                actual: flows.len(),
            }
            .into());
        }
        let sweeps = evolution.pair_sweeps();
        let (mut previous, mut max_grid, mut pairs) = (inputs, 0, 0_u64);
        let mut source_updates = 0_u64;
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
                    if flows.is_some_and(|flows| !flows[index].is_zero()) {
                        source_updates = (nodes as u64)
                            .checked_mul(interval as u64)
                            .and_then(|visits| {
                                visits.checked_mul(evolution.order().source_passes())
                            })
                            .and_then(|visits| source_updates.checked_add(visits))
                            .ok_or(VpmError::PlanOverflow {
                                resource: "scheduled affine source visits",
                            })?;
                    }
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
        let work = remesh_work + u128::from(pairs) + steps as u128 + u128::from(source_updates);
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
            source_updates,
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
