//! Conservative two-dimensional particle-strength exchange (PSE).
//!
//! Each particle carries circulation Gamma [m²/s] and an independently supplied,
//! positive quadrature area A [m²]. Its sampled vorticity is Gamma/A [1/s]. For
//! constant kinematic viscosity nu [m²/s] and kernel width epsilon [m],
//!
//! `dGamma_i/dt = sum_j 4 nu/(pi epsilon^4) exp(-r_ij²/epsilon²)
//!                         (A_i Gamma_j - A_j Gamma_i)`.
//!
//! The Gaussian has second moment 2 in each coordinate, giving nu times the
//! Laplacian in the resolved-kernel limit. This is the integral-operator PSE
//! construction of Degond and Mas-Gallic, Math. Comp. 53 (1989), 485–507,
//! DOI 10.1090/S0025-5718-1989-0983559-9. Area is quadrature support, not the
//! Biot–Savart regularization core; neither is inferred from circulation.
//!
//! [`PseDiffusion::compact_wendland`] instead selects the positive compact
//! shape `(1-q)^4 (1+4q)`, q=r/epsilon<1, exactly zero outside that radius.
//! Its coefficient is `1008 nu/(5 pi epsilon^4)`: the radial integral
//! `integral_0^1 q^3 (1-q)^4 (1+4q) dq = 5/504` gives the same second
//! moment 2 and resolved Laplacian limit. The shape is Wendland's C2 radial
//! function, Adv. Comput. Math. 4 (1995), 389–396, DOI 10.1007/BF02123482;
//! this diffusion normalization is derived from its second moment here.
//! Width denotes a support radius for this law. Resolve that support by
//! quadrature and refine radius and particle spacing together; changing
//! kernel family at unchanged width does not give the same regularization.
//!
//! Each pair is integrated exactly with one equal-and-opposite transfer. A
//! forward half sweep followed by the reversed half sweep is a second-order
//! symmetric splitting of the fixed-site diffusion operator. Every pair map
//! conserves circulation and contracts the vorticity range and area-weighted
//! enstrophy, without an explicit diffusion CFL or post-step clipping. Splitting
//! error still requires time refinement. The moving-particle entry composes
//! D(dt/2), the original checked RK4 advection, and D(dt/2).
//!
//! Direct pair sweeps use O(N²) work and O(N) live storage; there is no dense
//! diffusion matrix or per-step allocation. Gaussian tails are not cut off.
//! Gaussian, relaxation and logarithmic fallback evaluations use the strict
//! deterministic elementary functions from `fs-math`, with no platform libm.
//! Continuum consistency requires a well-resolved, overlapping quadrature
//! cloud and domain padding. A finite cloud cannot diffuse beyond its supplied
//! sites. No wall condition, remeshing, particle insertion, distorted-cloud
//! accuracy guarantee, variable-viscosity stress law, or 3-D stretching is
//! inferred. Exact circulation is not exact linear/angular impulse on an
//! arbitrary finite cloud. The algebraic velocity core remains the original
//! advection model, not a Gaussian heat-kernel solution.

use super::*;
use fs_math::det;

/// Positive-time adaptivity on an existing compact diffusion graph.
pub mod adaptive;
pub(crate) mod moving;
/// Reusable spatial graphs for fixed-site compact-support diffusion.
pub mod sparse;
pub use sparse::{PreparationBudget, PreparationReport, PreparedPse};

/// Spatial kernel with unit-coordinate second moment 2 in each direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PseKernel {
    /// Infinite-support Gaussian, `4 exp(-q²)/pi`.
    Gaussian,
    /// Positive C2 Wendland shape, `1008 (1-q)^4 (1+4q)/(5 pi)` for q<1.
    /// It vanishes exactly at and beyond the declared support radius.
    CompactWendland,
}

/// Immutable constant-viscosity PSE law with an explicitly selected kernel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PseDiffusion {
    viscosity: f64,
    width: f64,
    kernel: PseKernel,
}

impl PseDiffusion {
    /// Admit finite nonnegative nu [m²/s] and positive kernel width [m].
    /// Zero viscosity is exact identity, not an inferred material property.
    ///
    /// # Errors
    /// Refuses non-finite, negative or zero-width input.
    pub fn new(kinematic_viscosity_m2_s: f64, kernel_width_m: f64) -> Result<Self, VpmError> {
        Self::with_kernel(
            kinematic_viscosity_m2_s,
            kernel_width_m,
            PseKernel::Gaussian,
        )
    }

    /// Select a positive compact-support diffusion kernel with the same
    /// normalized Laplacian limit as the Gaussian law.
    ///
    /// `support_radius_m` is the physical radius beyond which interactions
    /// vanish exactly, independent of the velocity core. The spatial law is
    /// `k_ij = 1008 nu/(5 pi radius^4) * (1-q)^4 (1+4q)` for q<1.
    /// Positive pair weights retain the same circulation conservation and
    /// vorticity/enstrophy contraction as the original pair integrator.
    /// This kernel is C2; adequate support overlap and quadrature resolution
    /// remain required for consistency with continuum diffusion.
    ///
    /// # Errors
    /// Refuses non-finite or negative viscosity and nonpositive/non-finite radius.
    pub fn compact_wendland(
        kinematic_viscosity_m2_s: f64,
        support_radius_m: f64,
    ) -> Result<Self, VpmError> {
        Self::with_kernel(
            kinematic_viscosity_m2_s,
            support_radius_m,
            PseKernel::CompactWendland,
        )
    }

    fn with_kernel(
        kinematic_viscosity_m2_s: f64,
        kernel_width_m: f64,
        kernel: PseKernel,
    ) -> Result<Self, VpmError> {
        if !kinematic_viscosity_m2_s.is_finite() || kinematic_viscosity_m2_s < 0.0 {
            return Err(invalid_scalar(
                "kinematic viscosity",
                kinematic_viscosity_m2_s,
                None,
            ));
        }
        if !kernel_width_m.is_finite() || kernel_width_m <= 0.0 {
            return Err(invalid_scalar(
                "diffusion kernel width",
                kernel_width_m,
                None,
            ));
        }
        Ok(Self {
            viscosity: kinematic_viscosity_m2_s,
            width: kernel_width_m,
            kernel,
        })
    }

    /// Declared kinematic viscosity [m²/s].
    #[must_use]
    pub const fn kinematic_viscosity_m2_s(self) -> f64 {
        self.viscosity
    }

    /// Declared Gaussian width or compact support radius [m], independent
    /// of the velocity core. The two kernel families have different shapes.
    #[must_use]
    pub const fn kernel_width_m(self) -> f64 {
        self.width
    }

    /// Explicit spatial kernel; the original constructor selects Gaussian.
    #[must_use]
    pub const fn kernel(self) -> PseKernel {
        self.kernel
    }

    /// Finite interaction radius [m], or None for the infinite Gaussian tail.
    #[must_use]
    pub const fn support_radius_m(self) -> Option<f64> {
        match self.kernel {
            PseKernel::Gaussian => None,
            PseKernel::CompactWendland => Some(self.width),
        }
    }
}

pub(super) fn validate_inputs(n: usize, areas: &[f64], dt: f64) -> Result<(), VpmError> {
    if areas.len() != n {
        return Err(VpmError::LengthMismatch {
            field: "particle areas",
            expected: n,
            actual: areas.len(),
        });
    }
    if !dt.is_finite() || dt < 0.0 {
        return Err(invalid_scalar("diffusion dt", dt, None));
    }
    Ok(())
}

pub(super) fn validate_area(area: f64, particle: usize) -> Result<(), VpmError> {
    if !area.is_finite() || area <= 0.0 {
        return Err(invalid_scalar("particle area", area, Some(particle)));
    }
    Ok(())
}

// One symmetric sweep visits N(N-1) pairs. The split advection path has two.
pub(super) fn pair_work(n: usize, steps: usize, sweeps: u64) -> Result<u64, VpmError> {
    let overflow = || VpmError::PlanOverflow {
        resource: "diffusion pair evaluations",
    };
    let n = u64::try_from(n).map_err(|_| overflow())?;
    let steps = u64::try_from(steps).map_err(|_| overflow())?;
    n.checked_mul(n.saturating_sub(1))
        .and_then(|v| v.checked_mul(steps))
        .and_then(|v| v.checked_mul(sweeps))
        .ok_or_else(overflow)
}

/// Diffuse circulation on the supplied fixed quadrature sites.
///
/// `dt` is nonnegative seconds. All areas are validated even for zero steps or
/// zero viscosity. Successful nonzero diffusion charges `S N(N-1)` pair visits
/// and `N + S(N(N-1)+1)` logical work units. Identity runs skip pair work.
/// The live payload is exactly one N-particle vector; borrowed areas are not
/// copied. Each failed run leaves input untouched and publishes no candidate.
///
/// # Errors
/// Refuses invalid inputs, local/ambient work or memory budgets, cancellation,
/// allocation failure and nonrepresentable finite circulation updates.
#[allow(clippy::too_many_arguments)]
pub fn diffuse_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    dt: f64,
    steps: usize,
    law: PseDiffusion,
    budget: VpmBudget,
) -> Result<VpmRun, VpmError> {
    validate_inputs(particles.len(), areas_m2, dt)?;
    admit_cap("particles", particles.len(), budget.max_particles)?;
    admit_cap("steps", steps, budget.max_steps)?;
    let active = dt > 0.0 && law.viscosity > 0.0;
    let pair_evaluations = if active {
        pair_work(particles.len(), steps, 1)?
    } else {
        0
    };
    if pair_evaluations > budget.max_pair_evaluations {
        return Err(VpmError::CapExceeded {
            resource: "pair evaluations",
            required: pair_evaluations.into(),
            maximum: budget.max_pair_evaluations.into(),
        });
    }
    let work_units = pair_evaluations
        .checked_add(particles.len() as u64)
        .and_then(|n| n.checked_add(steps as u64))
        .ok_or(VpmError::PlanOverflow {
            resource: "diffusion total work",
        })?;
    let peak_live_bytes = particles
        .len()
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .ok_or(VpmError::PlanOverflow {
            resource: "diffusion live bytes",
        })?;
    admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
    accountant.checkpoint("fs-vpm.diffusion-admission", cx)?;
    let mut state = reserved_vec("diffusion-state", particles.len())?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.diffusion-input", || {
            validate_particle(particle, index)?;
            validate_area(areas_m2[index], index)?;
            state.push(particle);
            Ok(())
        })?;
    }
    meter.boundary(&mut accountant, cx, "fs-vpm.diffusion-input-complete")?;
    let mut completed_pairs = 0;
    for _ in 0..steps {
        if active {
            symmetric_sweep(
                &mut state,
                areas_m2,
                law,
                dt,
                0.5,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
            )?;
        }
        meter.unit(&mut accountant, cx, "fs-vpm.diffusion-step", || Ok(()))?;
        meter.boundary(&mut accountant, cx, "fs-vpm.diffusion-step-complete")?;
    }
    meter.flush(&mut accountant, "fs-vpm.diffusion-finalize")?;
    if meter.completed != work_units
        || completed_pairs != pair_evaluations
        || accountant.consumption().cost_charged != work_units
    {
        return Err(VpmError::AccountingMismatch {
            planned: work_units,
            completed: meter.completed,
            charged: accountant.consumption().cost_charged,
        });
    }
    accountant.checkpoint("fs-vpm.diffusion-publication", cx)?;
    Ok(VpmRun {
        particles: state,
        steps_completed: steps,
        pair_evaluations,
        work_units,
        peak_live_bytes,
        budget: accountant.consumption(),
    })
}

/// Advect and diffuse the same circulation inventory under one transaction.
///
/// Areas move with their original particle and remain constant, as appropriate
/// to the admitted incompressible 2-D advection. Each step uses symmetric
/// diffusion halves around the original RK4 owner. Positive viscosity adds
/// `2 S N(N-1)` pair visits to the inviscid plan without extra live buffers.
/// Zero viscosity recovers the inviscid trajectory and work counts exactly.
/// Time and cloud/kernel refinement are distinct caller responsibilities.
///
/// # Errors
/// The same typed refusals as [`diffuse_with_cx`] and [`simulate_with_cx`].
#[allow(clippy::too_many_arguments)]
pub fn simulate_viscous_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    dt: f64,
    steps: usize,
    core: f64,
    law: PseDiffusion,
    budget: VpmBudget,
) -> Result<VpmRun, VpmError> {
    simulate_impl(
        cx,
        particles,
        dt,
        steps,
        core,
        budget,
        Some((areas_m2, law)),
    )
}

// dt and fraction stay separate: a tiny representable dt must not disappear
// before multiplication by a large physical diffusion rate.
#[allow(clippy::too_many_arguments)]
pub(super) fn symmetric_sweep(
    state: &mut [VortexParticle],
    areas: &[f64],
    law: PseDiffusion,
    dt: f64,
    pair_fraction: f64,
    visits: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    symmetric_sweep_observed(
        state,
        areas,
        law,
        dt,
        pair_fraction,
        visits,
        meter,
        accountant,
        cx,
        &mut |_, _, _, _| Ok(()),
    )
}

// Observe the existing signed transfer and its actual represented update in
// the same charged pair visit. The no-op specialization retains the original
// arithmetic/order without extra cloud passes, allocations or dynamic dispatch.
#[allow(clippy::too_many_arguments)]
pub(super) fn symmetric_sweep_observed(
    state: &mut [VortexParticle],
    areas: &[f64],
    law: PseDiffusion,
    dt: f64,
    pair_fraction: f64,
    visits: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
    observer: &mut impl FnMut(
        usize,
        [VortexParticle; 2],
        [VortexParticle; 2],
        f64,
    ) -> Result<(), VpmError>,
) -> Result<(), VpmError> {
    let kernel = PairKernel::new(law, dt, pair_fraction);
    let n = state.len();
    let mut visit = |i, j| {
        meter.unit(accountant, cx, "fs-vpm.diffusion-pair", || {
            *visits = visits.checked_add(1).ok_or(VpmError::PlanOverflow {
                resource: "executed diffusion pairs",
            })?;
            let amount = kernel.fraction(state[i], state[j], [areas[i], areas[j]]);
            exchange_fraction_observed(state, areas, amount, i, j, observer)
        })
    };
    for i in 0..n {
        for j in i + 1..n {
            visit(i, j)?;
        }
    }
    // Reverse the complete PAIR order, not merely the particle order.
    for i in (0..n).rev() {
        for j in (i + 1..n).rev() {
            visit(i, j)?;
        }
    }
    meter.boundary(accountant, cx, "fs-vpm.diffusion-sweep-complete")
}

struct PairKernel {
    shape: PseKernel,
    width: f64,
    width2: f64,
    rate_time: f64,
    log_rate_time: f64,
}
impl PairKernel {
    fn new(law: PseDiffusion, dt: f64, fraction: f64) -> Self {
        let width2 = law.width * law.width;
        let normalization = match law.kernel {
            PseKernel::Gaussian => 4.0 / PI,
            PseKernel::CompactWendland => 1008.0 / (5.0 * PI),
        };
        Self {
            shape: law.kernel,
            width: law.width,
            width2,
            rate_time: normalization * (law.viscosity / width2) * dt * fraction,
            log_rate_time: det::ln(normalization)
                + det::ln(law.viscosity)
                + det::ln(dt)
                + det::ln(fraction)
                - 4.0 * det::ln(law.width),
        }
    }

    fn fraction(&self, a: VortexParticle, b: VortexParticle, areas: [f64; 2]) -> f64 {
        let r2 = scaled_distance2(a, b, self.width);
        if r2.is_infinite() || (self.shape == PseKernel::CompactWendland && r2 >= 1.0) {
            return 0.0;
        }
        let (weight, edge, radius) = match self.shape {
            PseKernel::Gaussian => (det::exp(-r2), 0.0, 0.0),
            PseKernel::CompactWendland => {
                let radius = det::sqrt(r2);
                // Rationalization retains 1-sqrt(r2) near support, including
                // when the rounded square root itself becomes exactly one.
                let edge = (1.0 - r2) / (1.0 + radius);
                let square = edge * edge;
                (square * square * (1.0 + 4.0 * radius), edge, radius)
            }
        };
        let mut exponent =
            self.rate_time * (areas[0] / self.width2 + areas[1] / self.width2) * weight;
        if !exponent.is_normal() {
            // Avoid inf*0, lost narrow-kernel factors, and underflowed Gaussian
            // tails times large rates. The usual physical-scale path uses no logs.
            let scale = areas[0].max(areas[1]);
            let log_area = det::ln(scale) + det::ln(areas[0] / scale + areas[1] / scale);
            let log_weight = match self.shape {
                PseKernel::Gaussian => -r2,
                PseKernel::CompactWendland => 4.0 * det::ln(edge) + det::ln(1.0 + 4.0 * radius),
            };
            exponent = det::exp(self.log_rate_time + log_area + log_weight);
        }
        // expm1 retains weak diffusion and tends to exact equilibration at
        // unrepresentably large rates; it never introduces a negative weight.
        -det::expm1(-exponent)
    }
}

// Shared support metric for the dense owner and prepared neighbor graph.
// Validated finite positions and a positive finite width give a nonnegative
// value or infinity; enormous separations carry no representable interaction.
fn scaled_distance2(a: VortexParticle, b: VortexParticle, width: f64) -> f64 {
    let mut r2 = 0.0;
    for axis in 0..2 {
        let delta = a.pos[axis] - b.pos[axis];
        let scaled = if delta.is_finite() {
            delta / width
        } else {
            a.pos[axis] / width - b.pos[axis] / width
        };
        r2 += scaled * scaled;
    }
    r2
}

// Prepared graphs cache amount once per run and reuse this exact transfer.
fn exchange_fraction(
    state: &mut [VortexParticle],
    areas: &[f64],
    amount: f64,
    i: usize,
    j: usize,
) -> Result<(), VpmError> {
    exchange_fraction_observed(state, areas, amount, i, j, &mut |_, _, _, _| Ok(()))
}

fn exchange_fraction_observed(
    state: &mut [VortexParticle],
    areas: &[f64],
    amount: f64,
    i: usize,
    j: usize,
    observer: &mut impl FnMut(
        usize,
        [VortexParticle; 2],
        [VortexParticle; 2],
        f64,
    ) -> Result<(), VpmError>,
) -> Result<(), VpmError> {
    // An absent or unrepresentably weak interaction is exact identity,
    // including signed zeros; disconnected nodes need no transfer arithmetic.
    if amount == 0.0 {
        return Ok(());
    }
    let scale = areas[i].max(areas[j]);
    let a = areas[i] / scale;
    let b = areas[j] / scale;
    let incoming = if a.is_normal() {
        (a / (a + b)) * state[j].circulation
    } else {
        small_area_share(areas[i], areas[j], state[j].circulation)
    };
    let outgoing = if b.is_normal() {
        (b / (a + b)) * state[i].circulation
    } else {
        small_area_share(areas[j], areas[i], state[i].circulation)
    };
    let transfer = amount * (incoming - outgoing);
    let next_i = state[i].circulation + transfer;
    let next_j = state[j].circulation - transfer;
    if !transfer.is_finite() || !next_i.is_finite() || !next_j.is_finite() {
        return Err(VpmError::NonFiniteResult {
            phase: "fs-vpm.diffusion-pair",
            particle: i,
        });
    }
    let before = [state[i], state[j]];
    state[i].circulation = next_i;
    state[j].circulation = next_j;
    observer(i, before, [state[i], state[j]], transfer)
}

// Only used when small/large is subnormal or has underflowed, so
// 1 + small/large rounds to exactly one. A normalized area fraction can lose
// the entire transfer even when multiplication by circulation would restore
// a representable result. Divide circulation by its large support first;
// logarithms are needed only if that intermediate also loses representation.
fn small_area_share(small: f64, large: f64, circulation: f64) -> f64 {
    if circulation == 0.0 {
        return circulation;
    }
    let vorticity = circulation / large;
    if vorticity.is_normal() {
        small * vorticity
    } else {
        det::exp(det::ln(circulation.abs()) + det::ln(small) - det::ln(large)).copysign(circulation)
    }
}
