//! Held incompressible affine flow composed with the native vortex owners.
//!
//! For `u(x) = U + A (x-o)`, `A = [[a,b],[c,-a]]`, the background map is exact
//! up to floating-point arithmetic. Symmetric `F(dt/2) N(dt) F(dt/2)` splitting
//! retains second-order time accuracy of the original conservative or viscous
//! owner; inviscid evolution also offers the fourth-order reversible composition.
//! This prescribed external flow can change vortex impulse; those
//! changes are measured, not projected away or called closed-system momentum.

use super::*;
use crate::diffusion::{self, PseDiffusion};
use fs_math::{dd::Dd, det};

/// A held divergence-free background in metres and seconds.
///
/// `u(x) = velocity + [[a,b],[c,-a]] * (x-origin)`. Three gradient entries
/// make zero divergence structural, without snapping an almost traceless
/// matrix. Translation, rigid rotation, shear and planar extension are covered.
/// This is an externally maintained flow, not a solved wall or body response.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AffineFlow2 {
    origin: [f64; 2],
    velocity: [f64; 2],
    gradient: [f64; 3],
}

impl AffineFlow2 {
    /// Declare a finite source. `gradient_s_inv = [a,b,c]` denotes
    /// `du/dx=a`, `du/dy=b`, `dv/dx=c` and `dv/dy=-a`.
    pub fn new(
        origin_m: [f64; 2],
        velocity_m_s: [f64; 2],
        gradient_s_inv: [f64; 3],
    ) -> Result<Self, VpmError> {
        for (name, values) in [
            ("ambient origin", origin_m.as_slice()),
            ("ambient velocity", velocity_m_s.as_slice()),
            ("ambient gradient", gradient_s_inv.as_slice()),
        ] {
            for &value in values {
                if !value.is_finite() {
                    return Err(invalid_scalar(name, value, None));
                }
            }
        }
        Ok(Self {
            origin: origin_m,
            velocity: velocity_m_s,
            gradient: gradient_s_inv,
        })
    }

    pub(super) fn is_zero(self) -> bool {
        self.velocity == [0.0; 2] && self.gradient == [0.0; 3]
    }
}

/// Actual vortex-moment increments during the background subflows alone.
///
/// These signed moments are about the world origin. They are not finite-volume
/// fluid momentum or mechanical work. PSE changes moments separately; pair
/// advection retains its original roundoff-level impulse defect.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AmbientMomentChange {
    /// Change in `sum Gamma * x` and `sum Gamma * y` [m³/s].
    pub first_moment_m3_s: [f64; 2],
    /// Change in `sum Gamma * |x|²` [m⁴/s].
    pub quadratic_moment_m4_s: f64,
}

/// Complete native trajectory plus the external source's moment account.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientVpmRun {
    /// Original resource receipt and the sole owned particle vector.
    pub run: VpmRun,
    /// Actual source-driven increments, accumulated without subtracting large
    /// whole-cloud inventories or dividing by net circulation.
    pub ambient_moment_change: AmbientMomentChange,
    /// Visits per particle and nonzero-time step for an active source: two at
    /// second order, four at fourth order. Included in `run.work_units`, not
    /// pair evaluations.
    pub source_particle_updates: u64,
}

/// Advance inviscid vortices in one held affine background. Negative time is
/// supported. A zero background, zero time, zero steps or empty cloud delegates
/// to the original owner, preserving its result and resource counts exactly.
///
/// The nontrivial source adds `2*N*S` work units and no particle buffers to the
/// original plan. One constant-size exact half-map is prepared per run; all
/// subsequent source updates are O(N), with no per-step allocation. Original
/// pair, local/ambient budget and cancellation owners remain authoritative.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    flow: AffineFlow2,
    budget: VpmBudget,
) -> Result<AmbientVpmRun, VpmError> {
    simulate_with_order_cx(
        cx,
        particles,
        dt,
        steps,
        core,
        flow,
        CompositionOrder::Second,
        budget,
    )
}

/// Choose second- or fourth-order inviscid evolution in one held affine field.
///
/// Fourth order composes three original reversible source/pair steps. Adjacent
/// exact maps of the same held source combine, so one physical step uses three
/// pair sweeps and four source passes. Both distinct signed source maps are
/// prepared before evolution. The sole particle vector, original pair owner,
/// compensated source-moment account and execution budget remain shared.
///
/// The negative internal interval is meaningful only for reversible inviscid
/// motion in a held field. No viscous, remeshing, varying-source or finite-domain
/// reaction claim follows. Source moment changes include every signed stage;
/// fourth order does not impose inviscid impulse invariants on a forced cloud.
///
/// # Errors
/// Original input/resource/numerical/cancellation refusals, including any
/// unrepresentable internal pair update or prepared signed affine map. All
/// failures publish no partial cloud. Empty/zero-source paths retain the chosen
/// conservative owner's complete result and resource counts.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_order_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    flow: AffineFlow2,
    order: CompositionOrder,
    budget: VpmBudget,
) -> Result<AmbientVpmRun, VpmError> {
    run(
        cx,
        particles,
        dt,
        steps,
        core,
        flow,
        Evolution::Inviscid(order),
        budget,
    )
}

/// Compose the held source around original `D(dt/2) A(dt) D(dt/2)` evolution.
/// Time is nonnegative and positive quadrature areas stay attached to their
/// particles, since every source map has unit area determinant. PSE retains
/// its finite-cloud spatial error; neither source nor diffusion is projected
/// to an inviscid constant-impulse state. Zero viscosity matches the inviscid
/// entry within the nonnegative-time domain, including its accounting.
#[allow(clippy::too_many_arguments)]
pub fn simulate_viscous_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    dt: f64,
    steps: usize,
    core: f64,
    flow: AffineFlow2,
    law: PseDiffusion,
    budget: VpmBudget,
) -> Result<AmbientVpmRun, VpmError> {
    run(
        cx,
        particles,
        dt,
        steps,
        core,
        flow,
        Evolution::Viscous(areas_m2, law),
        budget,
    )
}

// A fourth-order request cannot carry a diffusion law. This is the only
// dispatch into the signed composition; viscous callers retain the old step.
#[derive(Clone, Copy)]
enum Evolution<'a> {
    Inviscid(CompositionOrder),
    Viscous(&'a [f64], PseDiffusion),
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    flow: AffineFlow2,
    evolution: Evolution<'_>,
    budget: VpmBudget,
) -> Result<AmbientVpmRun, VpmError> {
    let (viscosity, order) = match evolution {
        Evolution::Inviscid(order) => (None, order),
        Evolution::Viscous(areas, law) => (Some((areas, law)), CompositionOrder::Second),
    };
    if flow.is_zero() || dt == 0.0 || steps == 0 || particles.is_empty() {
        let run = match viscosity {
            Some((areas, law)) => super::viscous::simulate_with_cx(
                cx, particles, areas, dt, steps, core, law, budget,
            )?,
            None => super::simulate_with_order_cx(cx, particles, dt, steps, core, order, budget)?,
        };
        return Ok(AmbientVpmRun {
            run,
            ambient_moment_change: AmbientMomentChange::default(),
            source_particle_updates: 0,
        });
    }
    if let Some((areas, _)) = viscosity {
        diffusion::validate_inputs(particles.len(), areas, dt)?;
    } else if !dt.is_finite() {
        return Err(invalid_scalar("dt", dt, None));
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None));
    }
    admit_cap("particles", particles.len(), budget.max_particles)?;
    admit_cap("steps", steps, budget.max_steps)?;
    let sweeps = if viscosity.is_some_and(|(_, law)| law.kinematic_viscosity_m2_s() > 0.0) {
        3
    } else {
        order.sweeps() as u64
    };
    let pair_evaluations = diffusion::pair_work(particles.len(), steps, sweeps)?;
    if pair_evaluations > budget.max_pair_evaluations {
        return Err(VpmError::CapExceeded {
            resource: "pair evaluations",
            required: pair_evaluations.into(),
            maximum: budget.max_pair_evaluations.into(),
        });
    }
    let source_particle_updates = (particles.len() as u64)
        .checked_mul(steps as u64)
        .and_then(|value| value.checked_mul(order.source_passes()))
        .ok_or(VpmError::PlanOverflow {
            resource: "ambient source visits",
        })?;
    let work_units = pair_evaluations
        .checked_add(source_particle_updates)
        .and_then(|value| value.checked_add(particles.len() as u64))
        .and_then(|value| value.checked_add(steps as u64))
        .ok_or(VpmError::PlanOverflow {
            resource: "ambient total work",
        })?;
    let peak_live_bytes = particles
        .len()
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .ok_or(VpmError::PlanOverflow {
            resource: "ambient live bytes",
        })?;
    admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
    accountant.checkpoint("fs-vpm.ambient-admission", cx)?;
    let composition = Composition::new(order);
    let half = HalfFlow::weighted(flow, dt, composition.weights()[0])?;
    let bridge = composition
        .bridge_weight()
        .map(|weight| HalfFlow::weighted(flow, dt, weight))
        .transpose()?;
    let mut state = reserved_vec("ambient-state", particles.len())?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.ambient-input", || {
            validate_particle(particle, index)?;
            if let Some((areas, _)) = viscosity {
                diffusion::validate_area(areas[index], index)?;
            }
            state.push(particle);
            Ok(())
        })?;
    }
    meter.boundary(&mut accountant, cx, "fs-vpm.ambient-input-complete")?;
    let mut completed_pairs = 0;
    let mut moments = [Dd::ZERO; 3];
    for _ in 0..steps {
        half.apply(&mut state, &mut moments, &mut meter, &mut accountant, cx)?;
        match viscosity {
            Some((areas, law)) if law.kinematic_viscosity_m2_s() > 0.0 => super::viscous::step(
                &mut state,
                areas,
                dt,
                core,
                law,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
            )?,
            _ => {
                for (index, &weight) in composition.weights().iter().enumerate() {
                    if index != 0 {
                        // Only a three-stage inviscid composition has bridges.
                        if let Some(bridge) = &bridge {
                            bridge.apply(
                                &mut state,
                                &mut moments,
                                &mut meter,
                                &mut accountant,
                                cx,
                            )?;
                        }
                    }
                    super::weighted_symmetric_sweep(
                        &mut state,
                        dt,
                        core,
                        weight,
                        &mut completed_pairs,
                        &mut meter,
                        &mut accountant,
                        cx,
                    )?;
                }
                meter.unit(&mut accountant, cx, "fs-vpm.conservative-step", || Ok(()))?;
            }
        }
        half.apply(&mut state, &mut moments, &mut meter, &mut accountant, cx)?;
        meter.boundary(&mut accountant, cx, "fs-vpm.ambient-step-complete")?;
    }
    meter.flush(&mut accountant, "fs-vpm.ambient-finalize")?;
    let consumption = accountant.consumption();
    if meter.completed != work_units
        || completed_pairs != pair_evaluations
        || consumption.cost_charged != work_units
    {
        return Err(VpmError::AccountingMismatch {
            planned: work_units,
            completed: meter.completed,
            charged: consumption.cost_charged,
        });
    }
    accountant.checkpoint("fs-vpm.ambient-publication", cx)?;
    Ok(AmbientVpmRun {
        run: VpmRun {
            particles: state,
            steps_completed: steps,
            pair_evaluations,
            work_units,
            peak_live_bytes,
            budget: accountant.consumption(),
        },
        ambient_moment_change: AmbientMomentChange {
            first_moment_m3_s: [moments[0].to_f64(), moments[1].to_f64()],
            quadratic_moment_m4_s: moments[2].to_f64(),
        },
        source_particle_updates,
    })
}

// Cayley–Hamilton gives B² = d I for B = dt*A/2. For small d (including
// exactly nilpotent shear), phi_1(B) = s I + t B is evaluated by its entire
// series. Elliptic flow uses stable sine/cosine functions. Hyperbolic flow uses
// spectral projectors: subtracting cosh and sinh would erase contraction.
pub(super) struct HalfFlow {
    origin: [f64; 2],
    generator: [[Dd; 2]; 2],
    offset: [Dd; 2],
    map: HalfMap,
}

enum HalfMap {
    Increment {
        matrix: [[Dd; 2]; 2],
        shift: [Dd; 2],
    },
    Hyperbolic {
        projectors: [[[Dd; 2]; 2]; 2],
        multipliers: [Dd; 2],
        equilibria: [[Dd; 2]; 2],
    },
}

impl HalfFlow {
    pub(super) fn new(flow: AffineFlow2, dt: f64) -> Result<Self, VpmError> {
        Self::weighted(flow, dt, 1.0)
    }

    // Form the original DD half-product before applying a composition weight:
    // rounding weight*dt first can erase resolvable subnormal-time transport.
    // Weight one takes the original arithmetic exactly.
    #[allow(clippy::float_cmp)]
    fn weighted(flow: AffineFlow2, dt: f64, weight: f64) -> Result<Self, VpmError> {
        let product = |value| {
            let half = half_product(value, dt);
            if weight == 1.0 {
                half
            } else {
                let scaled = half * dd(weight);
                if value != 0.0 && scaled == Dd::ZERO {
                    dd(f64::NAN)
                } else {
                    scaled
                }
            }
        };
        let [a, b, c] = flow.gradient.map(product);
        let offset = flow.velocity.map(product);
        let delta = a * a + b * c;
        let generator = [[a, b], [c, -a]];
        if ![a, b, c, offset[0], offset[1], delta]
            .into_iter()
            .all(finite)
        {
            return Err(range_error("fs-vpm.ambient-map", 0));
        }
        let map = if delta.to_f64() <= 0.0625 {
            let (s, t) = if delta.abs().to_f64() <= 0.0625 {
                let (mut s, mut t) = (Dd::ONE, Dd::from_f64(0.5));
                let (mut term_s, mut term_t) = (s, t);
                for k in 1..=6 {
                    term_s = term_s * delta / dd((2 * k * (2 * k + 1)) as f64);
                    term_t = term_t * delta / dd(((2 * k + 1) * (2 * k + 2)) as f64);
                    s = s + term_s;
                    t = t + term_t;
                }
                (s, t)
            } else {
                let angle = (-delta).sqrt();
                let (sin, _) = sin_cos_dd(angle);
                let (half_sin, _) = sin_cos_dd(angle * dd(0.5));
                (sin / angle, dd(2.0) * half_sin * half_sin / (-delta))
            };
            let diagonal = delta * t;
            let matrix = [[diagonal + s * a, s * b], [s * c, diagonal - s * a]];
            let product = matvec(generator, offset);
            let shift = [
                s * offset[0] + t * product[0],
                s * offset[1] + t * product[1],
            ];
            if !matrix.iter().flatten().copied().all(finite) || !shift.into_iter().all(finite) {
                return Err(range_error("fs-vpm.ambient-map", 0));
            }
            HalfMap::Increment { matrix, shift }
        } else {
            let rate = delta.sqrt();
            let plus = exp_dd(rate);
            let minus = exp_dd(-rate);
            let p = [
                [(Dd::ONE + a / rate) * dd(0.5), b / rate * dd(0.5)],
                [c / rate * dd(0.5), (Dd::ONE - a / rate) * dd(0.5)],
            ];
            let q = [[Dd::ONE - p[0][0], -p[0][1]], [-p[1][0], Dd::ONE - p[1][1]]];
            let projectors = [p, q];
            let multipliers = [plus, minus];
            let equilibria = [
                matvec(p, offset).map(|value| -value / rate),
                matvec(q, offset).map(|value| value / rate),
            ];
            if !projectors.iter().flatten().flatten().copied().all(finite)
                || !multipliers.into_iter().all(finite)
                || !equilibria.iter().flatten().copied().all(finite)
            {
                return Err(range_error("fs-vpm.ambient-map", 0));
            }
            HalfMap::Hyperbolic {
                projectors,
                multipliers,
                equilibria,
            }
        };
        Ok(Self {
            // A uniform source does not depend on its authored chart origin.
            // Avoid an irrelevant overflowing subtraction at distant origins.
            origin: if flow.gradient == [0.0; 3] {
                [0.0; 2]
            } else {
                flow.origin
            },
            generator,
            offset,
            map,
        })
    }

    pub(super) fn apply(
        &self,
        state: &mut [VortexParticle],
        moments: &mut [Dd; 3],
        meter: &mut WorkMeter,
        accountant: &mut AdmittedBudget<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), VpmError> {
        for (index, particle) in state.iter_mut().enumerate() {
            meter.unit(accountant, cx, "fs-vpm.ambient-source", || {
                let before = particle.pos;
                let relative = [
                    dd(before[0]) - dd(self.origin[0]),
                    dd(before[1]) - dd(self.origin[1]),
                ];
                let tangent = matvec(self.generator, relative);
                if tangent[0] + self.offset[0] == Dd::ZERO
                    && tangent[1] + self.offset[1] == Dd::ZERO
                {
                    return Ok(());
                }
                let endpoint: [Dd; 2] = match self.map {
                    HalfMap::Increment { matrix, shift } => {
                        let mapped = matvec(matrix, relative);
                        core::array::from_fn(|i| dd(before[i]) + mapped[i] + shift[i])
                    }
                    HalfMap::Hyperbolic {
                        projectors,
                        multipliers,
                        equilibria,
                    } => {
                        // Project BEFORE exponentiation: assembling a full
                        // matrix loses a contracting rotated eigenspace when
                        // expanding coefficients dominate. Center each source
                        // mode on its equilibrium for the same reason.
                        let modes = projectors.map(|projector| matvec(projector, relative));
                        core::array::from_fn(|i| {
                            dd(self.origin[i])
                                + equilibria[0][i]
                                + multipliers[0] * (modes[0][i] - equilibria[0][i])
                                + equilibria[1][i]
                                + multipliers[1] * (modes[1][i] - equilibria[1][i])
                        })
                    }
                };
                let after = endpoint.map(Dd::to_f64);
                if after.iter().any(|x| !x.is_finite()) {
                    return Err(range_error("fs-vpm.ambient-source", index));
                }
                if after == before && (endpoint[0] != dd(before[0]) || endpoint[1] != dd(before[1]))
                {
                    return Err(range_error("fs-vpm.ambient-position-resolution", index));
                }
                if particle.circulation != 0.0 {
                    for axis in 0..2 {
                        let change =
                            dd(particle.circulation) * (dd(after[axis]) - dd(before[axis]));
                        moments[axis] = moments[axis] + change;
                        moments[2] =
                            moments[2] + change * dd(before[axis]) + change * dd(after[axis]);
                    }
                    if !moments.iter().copied().all(finite) {
                        return Err(range_error("fs-vpm.ambient-moments", index));
                    }
                }
                particle.pos = after;
                Ok(())
            })?;
        }
        Ok(())
    }
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn range_error(phase: &'static str, particle: usize) -> VpmError {
    VpmError::NonFiniteResult { phase, particle }
}
fn matvec(matrix: [[Dd; 2]; 2], vector: [Dd; 2]) -> [Dd; 2] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1])
}

// Retain the product's low part. Reorder only if the full product overflows;
// a nonzero half coefficient that cannot be represented refuses, rather than
// silently suppressing a source whose action on a large coordinate may matter.
fn half_product(value: f64, dt: f64) -> Dd {
    if value == 0.0 {
        return Dd::ZERO;
    }
    let product = if (value * dt).is_finite() {
        dd(value) * dd(dt) * dd(0.5)
    } else {
        dd(value) * dd(dt * 0.5)
    };
    if product.hi == 0.0 && product.lo == 0.0 {
        dd(f64::NAN)
    } else {
        product
    }
}
fn sin_cos_dd(angle: Dd) -> (Dd, Dd) {
    let (sh, ch) = det::sin_cos(angle.hi);
    let (sl, cl) = det::sin_cos(angle.lo);
    (
        dd(sh) * dd(cl) + dd(ch) * dd(sl),
        dd(ch) * dd(cl) - dd(sh) * dd(sl),
    )
}
fn exp_dd(value: Dd) -> Dd {
    dd(det::exp(value.hi)) * (Dd::ONE + dd(det::expm1(value.lo)))
}
