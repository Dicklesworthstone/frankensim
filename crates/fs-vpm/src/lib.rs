//! fs-vpm — vortex particle method (2-D core). Layer: L3.
//!
//! Vorticity, not velocity, is the natural state variable for a wake: it is
//! compact (nonzero only where the fluid actually rotates) and it advects with
//! the flow. This v0 is the inviscid 2-D core — point vortices carrying
//! circulation, inducing velocity on each other by a DESINGULARIZED
//! BIOT–SAVART kernel, advanced with RK4.
//!
//! The load-bearing checks are analytic: a single vortex of circulation `Γ`
//! induces a purely TANGENTIAL field of magnitude `Γ/(2πr)` (and does not move
//! itself); a counter-rotating PAIR separated by `d` self-propels in a straight
//! line at exactly `Γ/(2πd)` (the 2-D analog of a vortex ring's translation);
//! and the total circulation is CONSERVED. The checked production entry point
//! admits particle, step, pair-work, logical live-memory, and ambient `Cx` budgets
//! before allocation, then polls cancellation inside every all-pairs sweep.

use core::f64::consts::PI;
use fs_exec::{AdmittedBudget, BudgetConsumption, BudgetRefusal, Cx};

/// Maximum completed logical work units between production-path checkpoints.
///
/// Cancellation latency is stated in logical work, not wall time. Pair
/// contributions, particle updates, and combined input-validation/copy visits
/// are each one unit.
pub const VPM_WORK_CHECKPOINT_STRIDE: u64 = 256;

/// A point vortex: a position carrying a circulation `Γ`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VortexParticle {
    /// The position.
    pub pos: [f64; 2],
    /// The circulation `Γ` (signed strength).
    pub circulation: f64,
}

impl VortexParticle {
    /// A vortex particle.
    #[must_use]
    pub fn new(pos: [f64; 2], circulation: f64) -> VortexParticle {
        VortexParticle { pos, circulation }
    }
}

/// Explicit local admission envelope for one production VPM run.
///
/// The ambient [`Cx`] budget is enforced independently. Both envelopes must
/// admit the same checked plan before work begins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VpmBudget {
    /// Maximum input particles.
    pub max_particles: usize,
    /// Maximum RK4 steps, including the zero-particle case where pair work is
    /// zero but step overhead is not.
    pub max_steps: usize,
    /// Maximum source-target contributions (`4 * particles^2 * steps`).
    pub max_pair_evaluations: u64,
    /// Maximum logical peak bytes owned by the run's particle/velocity
    /// payload buffers. Allocator metadata and capacity rounding are not
    /// claimed.
    pub max_live_bytes: usize,
}

impl VpmBudget {
    /// Construct one explicit local envelope.
    #[must_use]
    pub const fn new(
        max_particles: usize,
        max_steps: usize,
        max_pair_evaluations: u64,
        max_live_bytes: usize,
    ) -> Self {
        Self {
            max_particles,
            max_steps,
            max_pair_evaluations,
            max_live_bytes,
        }
    }
}

/// Why a checked VPM run refused admission or publication.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VpmError {
    /// A scalar input was NaN, infinite, or outside its admitted sign domain.
    InvalidScalar {
        /// Stable field name, including an indexed particle component when
        /// applicable.
        field: &'static str,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
        /// Particle index for particle-local fields.
        particle: Option<usize>,
    },
    /// Checked admission arithmetic could not represent the requested plan.
    PlanOverflow {
        /// Formula that overflowed.
        resource: &'static str,
    },
    /// A local explicit cap refused the checked request.
    CapExceeded {
        /// Bounded resource.
        resource: &'static str,
        /// Checked requirement.
        required: u128,
        /// Caller-authorized maximum.
        maximum: u128,
    },
    /// The ambient `Cx` budget refused admission or a checkpoint.
    ExecutionBudget(BudgetRefusal),
    /// A fallible vector reservation failed before writing the corresponding
    /// phase output.
    AllocationFailed {
        /// Allocation phase.
        phase: &'static str,
        /// Requested elements.
        elements: usize,
    },
    /// Finite admitted inputs produced a non-finite intermediate, so no state
    /// was published.
    NonFiniteResult {
        /// RK4 phase that lost the finite domain.
        phase: &'static str,
        /// Particle whose result was non-finite.
        particle: usize,
    },
    /// Internal progress failed to match the checked plan, so publication was
    /// refused rather than emitting an under- or over-charged receipt.
    AccountingMismatch {
        /// Planned logical work.
        planned: u64,
        /// Logical work completed by the private transaction.
        completed: u64,
        /// Work charged to the ambient budget.
        charged: u64,
    },
}

impl core::fmt::Display for VpmError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidScalar {
                field,
                bits,
                particle,
            } => {
                if let Some(particle) = particle {
                    write!(
                        formatter,
                        "invalid VPM scalar {field} for particle {particle}: bits=0x{bits:016x}"
                    )
                } else {
                    write!(formatter, "invalid VPM scalar {field}: bits=0x{bits:016x}")
                }
            }
            Self::PlanOverflow { resource } => {
                write!(formatter, "VPM admission formula overflowed for {resource}")
            }
            Self::CapExceeded {
                resource,
                required,
                maximum,
            } => write!(
                formatter,
                "VPM {resource} requires {required}, exceeding the explicit maximum {maximum}"
            ),
            Self::ExecutionBudget(refusal) => write!(formatter, "VPM {refusal}"),
            Self::AllocationFailed { phase, elements } => write!(
                formatter,
                "VPM allocation failed in {phase} while reserving {elements} elements"
            ),
            Self::NonFiniteResult { phase, particle } => write!(
                formatter,
                "VPM {phase} produced a non-finite result for particle {particle}; state refused"
            ),
            Self::AccountingMismatch {
                planned,
                completed,
                charged,
            } => write!(
                formatter,
                "VPM accounting mismatch: planned {planned}, completed {completed}, charged {charged}"
            ),
        }
    }
}

impl core::error::Error for VpmError {}

impl From<BudgetRefusal> for VpmError {
    fn from(refusal: BudgetRefusal) -> Self {
        Self::ExecutionBudget(refusal)
    }
}

/// Complete production-run result and the exact resource receipt enforced for
/// it. An error returns no partial state.
#[derive(Debug, Clone, PartialEq)]
pub struct VpmRun {
    /// Final particle state, published only after the final checkpoint.
    pub particles: Vec<VortexParticle>,
    /// Completed RK4 steps.
    pub steps_completed: usize,
    /// Exact source-target contributions charged by the run.
    pub pair_evaluations: u64,
    /// Exact logical work charged to the ambient `Cx`: input validation and
    /// copy, pair contributions, RK4 displacement/combine work, and one
    /// boundary unit per completed step.
    pub work_units: u64,
    /// Checked logical peak payload bytes admitted before allocation.
    pub peak_live_bytes: usize,
    /// Ambient `Cx` budget contract and exact successful consumption.
    pub budget: BudgetConsumption,
}

#[derive(Debug, Clone, Copy)]
struct VpmPlan {
    pair_evaluations: u64,
    work_units: u64,
    peak_live_bytes: usize,
}

impl VpmPlan {
    fn checked(particles: usize, steps: usize) -> Result<Self, VpmError> {
        let particle_count = u64::try_from(particles).map_err(|_| VpmError::PlanOverflow {
            resource: "particle count",
        })?;
        let step_count = u64::try_from(steps).map_err(|_| VpmError::PlanOverflow {
            resource: "step count",
        })?;
        let pair_evaluations = particle_count
            .checked_mul(particle_count)
            .and_then(|value| value.checked_mul(4))
            .and_then(|value| value.checked_mul(step_count))
            .ok_or(VpmError::PlanOverflow {
                resource: "pair evaluations",
            })?;
        let linear_step_work = particle_count
            .checked_mul(4)
            .ok_or(VpmError::PlanOverflow {
                resource: "linear RK4 work",
            })?;
        let step_work = particle_count
            .checked_mul(particle_count)
            .and_then(|value| value.checked_mul(4))
            .and_then(|value| value.checked_add(linear_step_work))
            .and_then(|value| value.checked_add(1))
            .ok_or(VpmError::PlanOverflow {
                resource: "per-step work",
            })?;
        let work_units = step_work
            .checked_mul(step_count)
            .and_then(|value| value.checked_add(particle_count))
            .ok_or(VpmError::PlanOverflow {
                resource: "total work",
            })?;

        let particle_buffers = if steps == 0 { 1 } else { 5 };
        let velocity_buffers = if steps == 0 { 0 } else { 4 };
        let bytes_per_particle = core::mem::size_of::<VortexParticle>()
            .checked_mul(particle_buffers)
            .and_then(|value| {
                core::mem::size_of::<[f64; 2]>()
                    .checked_mul(velocity_buffers)
                    .and_then(|velocity_bytes| value.checked_add(velocity_bytes))
            })
            .ok_or(VpmError::PlanOverflow {
                resource: "scratch bytes per particle",
            })?;
        let peak_live_bytes =
            particles
                .checked_mul(bytes_per_particle)
                .ok_or(VpmError::PlanOverflow {
                    resource: "scratch bytes",
                })?;

        Ok(Self {
            pair_evaluations,
            work_units,
            peak_live_bytes,
        })
    }
}

#[derive(Debug, Default)]
struct WorkMeter {
    pending: u64,
    completed: u64,
}

impl WorkMeter {
    fn unit<T>(
        &mut self,
        accountant: &mut AdmittedBudget<'_>,
        cx: &Cx<'_>,
        phase: &'static str,
        operation: impl FnOnce() -> Result<T, VpmError>,
    ) -> Result<T, VpmError> {
        if self.pending == 0 {
            accountant.checkpoint(phase, cx)?;
        }
        let result = operation();
        self.pending = self.pending.checked_add(1).ok_or(VpmError::PlanOverflow {
            resource: "pending work accounting",
        })?;
        self.completed = self
            .completed
            .checked_add(1)
            .ok_or(VpmError::PlanOverflow {
                resource: "completed work accounting",
            })?;
        if self.pending == VPM_WORK_CHECKPOINT_STRIDE || result.is_err() {
            self.flush(accountant, phase)?;
        }
        result
    }

    fn flush(
        &mut self,
        accountant: &mut AdmittedBudget<'_>,
        phase: &'static str,
    ) -> Result<(), VpmError> {
        if self.pending != 0 {
            accountant.charge_cost(phase, self.pending)?;
            self.pending = 0;
        }
        Ok(())
    }

    fn boundary(
        &mut self,
        accountant: &mut AdmittedBudget<'_>,
        cx: &Cx<'_>,
        phase: &'static str,
    ) -> Result<(), VpmError> {
        self.flush(accountant, phase)?;
        accountant.checkpoint(phase, cx)?;
        Ok(())
    }
}

struct RkScratch {
    k1: Vec<[f64; 2]>,
    k2: Vec<[f64; 2]>,
    k3: Vec<[f64; 2]>,
    k4: Vec<[f64; 2]>,
    s2: Vec<VortexParticle>,
    s3: Vec<VortexParticle>,
    s4: Vec<VortexParticle>,
    next: Vec<VortexParticle>,
}

impl RkScratch {
    fn try_new(elements: usize) -> Result<Self, VpmError> {
        Ok(Self {
            k1: reserved_vec("rk4-k1", elements)?,
            k2: reserved_vec("rk4-k2", elements)?,
            k3: reserved_vec("rk4-k3", elements)?,
            k4: reserved_vec("rk4-k4", elements)?,
            s2: reserved_vec("rk4-s2", elements)?,
            s3: reserved_vec("rk4-s3", elements)?,
            s4: reserved_vec("rk4-s4", elements)?,
            next: reserved_vec("rk4-next", elements)?,
        })
    }
}

fn reserved_vec<T>(phase: &'static str, elements: usize) -> Result<Vec<T>, VpmError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(elements)
        .map_err(|_| VpmError::AllocationFailed { phase, elements })?;
    Ok(values)
}

/// The velocity induced at `point` by all `particles` via the desingularized
/// 2-D Biot–Savart kernel `u = Σ (Γⱼ/2π)·perp(r)/(|r|² + core²)`, with
/// `perp([a,b]) = [−b, a]`. A `core` of `0` skips coincident particles.
#[must_use]
pub fn induced_velocity(particles: &[VortexParticle], point: [f64; 2], core: f64) -> [f64; 2] {
    let mut u = [0.0, 0.0];
    let core2 = core * core;
    for p in particles {
        add_induced_contribution(&mut u, p, point, core2);
    }
    u
}

fn add_induced_contribution(
    velocity: &mut [f64; 2],
    source: &VortexParticle,
    point: [f64; 2],
    core2: f64,
) {
    let r = [point[0] - source.pos[0], point[1] - source.pos[1]];
    let r2 = r[0] * r[0] + r[1] * r[1] + core2;
    if r2 <= 1e-30 {
        return; // coincident (self) with no desingularization
    }
    let g = source.circulation / (2.0 * PI * r2);
    velocity[0] += -g * r[1];
    velocity[1] += g * r[0];
}

fn velocities(particles: &[VortexParticle], core: f64) -> Vec<[f64; 2]> {
    particles
        .iter()
        .map(|p| induced_velocity(particles, p.pos, core))
        .collect()
}

fn displaced(particles: &[VortexParticle], vel: &[[f64; 2]], dt: f64) -> Vec<VortexParticle> {
    particles
        .iter()
        .zip(vel)
        .map(|(p, v)| {
            VortexParticle::new([p.pos[0] + dt * v[0], p.pos[1] + dt * v[1]], p.circulation)
        })
        .collect()
}

/// Advance the vortex system by one unchecked compatibility RK4 step of size
/// `dt` (circulations are invariant under inviscid advection).
///
/// This infallible helper exists for analytic fixtures and existing callers;
/// it has no explicit work, memory, finite-domain, or cancellation admission.
/// Use [`simulate_with_cx`] for production work.
#[must_use]
pub fn advect(particles: &[VortexParticle], dt: f64, core: f64) -> Vec<VortexParticle> {
    let k1 = velocities(particles, core);
    let s2 = displaced(particles, &k1, 0.5 * dt);
    let k2 = velocities(&s2, core);
    let s3 = displaced(particles, &k2, 0.5 * dt);
    let k3 = velocities(&s3, core);
    let s4 = displaced(particles, &k3, dt);
    let k4 = velocities(&s4, core);
    particles
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let vx = k1[i][0] + 2.0 * k2[i][0] + 2.0 * k3[i][0] + k4[i][0];
            let vy = k1[i][1] + 2.0 * k2[i][1] + 2.0 * k3[i][1] + k4[i][1];
            VortexParticle::new(
                [p.pos[0] + dt / 6.0 * vx, p.pos[1] + dt / 6.0 * vy],
                p.circulation,
            )
        })
        .collect()
}

/// Advance the system by `steps` unchecked compatibility RK4 steps.
///
/// This preserves the original analytic-fixture API. It may allocate without
/// an explicit cap and cannot observe cancellation; production callers should
/// use [`simulate_with_cx`].
#[must_use]
pub fn simulate(
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
) -> Vec<VortexParticle> {
    let mut state = particles.to_vec();
    for _ in 0..steps {
        state = advect(&state, dt, core);
    }
    state
}

fn invalid_scalar(field: &'static str, value: f64, particle: Option<usize>) -> VpmError {
    VpmError::InvalidScalar {
        field,
        bits: value.to_bits(),
        particle,
    }
}

fn admit_cap(resource: &'static str, required: usize, maximum: usize) -> Result<(), VpmError> {
    if required > maximum {
        return Err(VpmError::CapExceeded {
            resource,
            required: required as u128,
            maximum: maximum as u128,
        });
    }
    Ok(())
}

fn validate_particle(particle: VortexParticle, index: usize) -> Result<(), VpmError> {
    if !particle.pos[0].is_finite() {
        return Err(invalid_scalar("position.x", particle.pos[0], Some(index)));
    }
    if !particle.pos[1].is_finite() {
        return Err(invalid_scalar("position.y", particle.pos[1], Some(index)));
    }
    if !particle.circulation.is_finite() {
        return Err(invalid_scalar(
            "circulation",
            particle.circulation,
            Some(index),
        ));
    }
    Ok(())
}

fn checked_induced_contribution(
    velocity: &mut [f64; 2],
    source: &VortexParticle,
    point: [f64; 2],
    core2: f64,
    phase: &'static str,
    target: usize,
) -> Result<(), VpmError> {
    let r = [point[0] - source.pos[0], point[1] - source.pos[1]];
    if !r[0].is_finite() || !r[1].is_finite() {
        return Err(VpmError::NonFiniteResult {
            phase,
            particle: target,
        });
    }
    let r2 = r[0] * r[0] + r[1] * r[1] + core2;
    if !r2.is_finite() {
        return Err(VpmError::NonFiniteResult {
            phase,
            particle: target,
        });
    }
    if r2 <= 1e-30 {
        return Ok(());
    }
    let denominator = 2.0 * PI * r2;
    let g = source.circulation / denominator;
    let increments = [-g * r[1], g * r[0]];
    if !denominator.is_finite()
        || !g.is_finite()
        || !increments[0].is_finite()
        || !increments[1].is_finite()
    {
        return Err(VpmError::NonFiniteResult {
            phase,
            particle: target,
        });
    }
    velocity[0] += increments[0];
    velocity[1] += increments[1];
    if !velocity[0].is_finite() || !velocity[1].is_finite() {
        return Err(VpmError::NonFiniteResult {
            phase,
            particle: target,
        });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn checked_velocities(
    particles: &[VortexParticle],
    core2: f64,
    output: &mut Vec<[f64; 2]>,
    pair_evaluations: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
    phase: &'static str,
) -> Result<(), VpmError> {
    output.clear();
    for (target, particle) in particles.iter().enumerate() {
        let mut velocity = [0.0, 0.0];
        for source in particles {
            meter.unit(accountant, cx, phase, || {
                *pair_evaluations =
                    (*pair_evaluations)
                        .checked_add(1)
                        .ok_or(VpmError::PlanOverflow {
                            resource: "executed pair evaluations",
                        })?;
                checked_induced_contribution(
                    &mut velocity,
                    source,
                    particle.pos,
                    core2,
                    phase,
                    target,
                )
            })?;
        }
        output.push(velocity);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn checked_displaced(
    particles: &[VortexParticle],
    velocities: &[[f64; 2]],
    dt: f64,
    output: &mut Vec<VortexParticle>,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
    phase: &'static str,
) -> Result<(), VpmError> {
    output.clear();
    for (index, (particle, velocity)) in particles.iter().zip(velocities).enumerate() {
        meter.unit(accountant, cx, phase, || {
            let pos = [
                particle.pos[0] + dt * velocity[0],
                particle.pos[1] + dt * velocity[1],
            ];
            if !pos[0].is_finite() || !pos[1].is_finite() {
                return Err(VpmError::NonFiniteResult {
                    phase,
                    particle: index,
                });
            }
            output.push(VortexParticle::new(pos, particle.circulation));
            Ok(())
        })?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn checked_combined(
    particles: &[VortexParticle],
    dt: f64,
    k1: &[[f64; 2]],
    k2: &[[f64; 2]],
    k3: &[[f64; 2]],
    k4: &[[f64; 2]],
    output: &mut Vec<VortexParticle>,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    const PHASE: &str = "fs-vpm.rk4-combine";
    output.clear();
    for (index, particle) in particles.iter().enumerate() {
        meter.unit(accountant, cx, PHASE, || {
            let vx = k1[index][0] + 2.0 * k2[index][0] + 2.0 * k3[index][0] + k4[index][0];
            let vy = k1[index][1] + 2.0 * k2[index][1] + 2.0 * k3[index][1] + k4[index][1];
            let pos = [
                particle.pos[0] + dt / 6.0 * vx,
                particle.pos[1] + dt / 6.0 * vy,
            ];
            if !vx.is_finite() || !vy.is_finite() || !pos[0].is_finite() || !pos[1].is_finite() {
                return Err(VpmError::NonFiniteResult {
                    phase: PHASE,
                    particle: index,
                });
            }
            output.push(VortexParticle::new(pos, particle.circulation));
            Ok(())
        })?;
    }
    Ok(())
}

/// Run the direct RK4 VPM transaction under explicit local and ambient
/// admission.
///
/// The input is borrowed and never modified. Every allocation is fallible and
/// occurs after the first ambient checkpoint; every output remains private
/// until exact work accounting and a final cancellation/deadline checkpoint
/// succeed. Consequently every error returns no partial particle state.
///
/// Logical work is charged as one combined validation/copy unit per input,
/// four source-target visits and four particle updates per RK4 step, plus one
/// completed-step boundary unit. Logical live bytes exclude the borrowed
/// input, allocator metadata, and capacity rounding.
///
/// # Errors
/// Refuses invalid numeric inputs, checked-plan overflow, either explicit
/// local cap, ambient budget/cancellation/deadline exhaustion, allocation
/// failure, or any non-finite numerical intermediate.
#[allow(clippy::too_many_lines)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    budget: VpmBudget,
) -> Result<VpmRun, VpmError> {
    if !dt.is_finite() {
        return Err(invalid_scalar("dt", dt, None));
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None));
    }
    let core2 = core * core;
    if !core2.is_finite() {
        return Err(invalid_scalar("core", core, None));
    }

    admit_cap("particles", particles.len(), budget.max_particles)?;
    admit_cap("steps", steps, budget.max_steps)?;
    let plan = VpmPlan::checked(particles.len(), steps)?;
    if plan.pair_evaluations > budget.max_pair_evaluations {
        return Err(VpmError::CapExceeded {
            resource: "pair evaluations",
            required: u128::from(plan.pair_evaluations),
            maximum: u128::from(budget.max_pair_evaluations),
        });
    }
    admit_cap(
        "logical live bytes",
        plan.peak_live_bytes,
        budget.max_live_bytes,
    )?;

    let mut accountant = AdmittedBudget::admit_ambient(cx, plan.work_units)?;
    accountant.checkpoint("fs-vpm.admission", cx)?;

    let mut state = reserved_vec("state", particles.len())?;
    let mut scratch = RkScratch::try_new(if steps == 0 { 0 } else { particles.len() })?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.input", || {
            validate_particle(particle, index)?;
            state.push(particle);
            Ok(())
        })?;
    }
    meter.boundary(&mut accountant, cx, "fs-vpm.input-complete")?;

    let mut pair_evaluations = 0_u64;
    for _ in 0..steps {
        checked_velocities(
            &state,
            core2,
            &mut scratch.k1,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-k1",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-k1-complete")?;

        checked_displaced(
            &state,
            &scratch.k1,
            0.5 * dt,
            &mut scratch.s2,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-s2",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-s2-complete")?;
        checked_velocities(
            &scratch.s2,
            core2,
            &mut scratch.k2,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-k2",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-k2-complete")?;

        checked_displaced(
            &state,
            &scratch.k2,
            0.5 * dt,
            &mut scratch.s3,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-s3",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-s3-complete")?;
        checked_velocities(
            &scratch.s3,
            core2,
            &mut scratch.k3,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-k3",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-k3-complete")?;

        checked_displaced(
            &state,
            &scratch.k3,
            dt,
            &mut scratch.s4,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-s4",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-s4-complete")?;
        checked_velocities(
            &scratch.s4,
            core2,
            &mut scratch.k4,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
            "fs-vpm.rk4-k4",
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-k4-complete")?;

        checked_combined(
            &state,
            dt,
            &scratch.k1,
            &scratch.k2,
            &scratch.k3,
            &scratch.k4,
            &mut scratch.next,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        meter.boundary(&mut accountant, cx, "fs-vpm.rk4-combine-complete")?;
        meter.unit(&mut accountant, cx, "fs-vpm.step-commit", || Ok(()))?;
        core::mem::swap(&mut state, &mut scratch.next);
        scratch.next.clear();
        meter.boundary(&mut accountant, cx, "fs-vpm.step-complete")?;
    }

    meter.flush(&mut accountant, "fs-vpm.finalize")?;
    let consumption = accountant.consumption();
    if meter.completed != plan.work_units
        || pair_evaluations != plan.pair_evaluations
        || consumption.cost_charged != plan.work_units
    {
        return Err(VpmError::AccountingMismatch {
            planned: plan.work_units,
            completed: meter.completed,
            charged: consumption.cost_charged,
        });
    }
    accountant.checkpoint("fs-vpm.publication", cx)?;
    let consumption = accountant.consumption();

    Ok(VpmRun {
        particles: state,
        steps_completed: steps,
        pair_evaluations,
        work_units: plan.work_units,
        peak_live_bytes: plan.peak_live_bytes,
        budget: consumption,
    })
}

/// The total circulation `Σ Γᵢ` (a conserved invariant).
#[must_use]
pub fn total_circulation(particles: &[VortexParticle]) -> f64 {
    particles.iter().map(|p| p.circulation).sum()
}

/// The centroid of vorticity `Σ Γᵢ xᵢ / Σ Γᵢ` (the "center of vorticity";
/// invariant for an isolated inviscid system). Returns `None` when the total
/// circulation is zero (e.g. a symmetric pair — use per-particle tracking).
#[must_use]
pub fn vorticity_centroid(particles: &[VortexParticle]) -> Option<[f64; 2]> {
    let total = total_circulation(particles);
    if total.abs() <= 1e-30 {
        return None;
    }
    let mut c = [0.0, 0.0];
    for p in particles {
        c[0] += p.circulation * p.pos[0];
        c[1] += p.circulation * p.pos[1];
    }
    Some([c[0] / total, c[1] / total])
}
