//! Adaptive, positive-time evolution on an existing compact PSE graph.
//!
//! A trial compares one symmetric pair sweep with two half sweeps. The
//! accepted candidate is the two-half-step state, never an extrapolation:
//! every published step retains the original conservative positive pair maps.
//! The local error estimate is |fine-coarse|/3, normalized by explicit absolute
//! vorticity and relative circulation tolerances. This controls time splitting
//! only; it does not estimate spatial quadrature, support or boundary error.
//!
//! Rejected trials halve the step, well-resolved trials may double it. All
//! work, including rejected trials, is charged. Three particle buffers and one
//! pair of fractions per graph edge are reserved once; no trial allocates.

use super::{PairKernel, PreparedPse, exchange_fraction};
use crate::{
    VortexParticle, VpmError, WorkMeter, admit_cap, invalid_scalar, reserved_vec, validate_particle,
};
use fs_exec::{AdmittedBudget, BudgetConsumption, BudgetRefusal, Cx};
use fs_math::det;

/// Explicit time window and local error controls; all times are in seconds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveOptions {
    /// First trial step, between the minimum and maximum.
    pub initial_dt_s: f64,
    /// Smallest retry step. A final endpoint remainder may be shorter.
    pub min_dt_s: f64,
    /// Largest trial step.
    pub max_dt_s: f64,
    /// Strictly positive absolute tolerance in sampled vorticity [1/s].
    pub absolute_vorticity_tolerance_s_inv: f64,
    /// Nonnegative relative tolerance on circulation.
    pub relative_tolerance: f64,
}

/// Independent caps for the entire requested window, including rejected work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdaptiveBudget {
    /// Maximum particles.
    pub max_particles: usize,
    /// Maximum complete or attempted coarse/fine trials.
    pub max_attempts: usize,
    /// Maximum pair updates; every active trial uses exactly 6 E.
    pub max_pair_evaluations: u64,
    /// Logical work ceiling admitted against the ambient context up front.
    pub max_work_units: u64,
    /// Maximum owned vector payload, excluding the borrowed prepared graph.
    pub max_live_bytes: usize,
}

/// Complete requested window. No partial candidate is returned on failure.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveRun {
    /// Accepted final particle inventory.
    pub particles: Vec<VortexParticle>,
    /// Exactly the requested finite duration [s].
    pub elapsed_s: f64,
    /// Accepted trials; each consists of two positive half steps.
    pub accepted_steps: usize,
    /// Rejected trials, which never modify the accepted inventory.
    pub rejected_steps: usize,
    /// All accepted and rejected trials.
    pub attempts: usize,
    /// Actual pair visits, including rejected candidates.
    pub pair_evaluations: u64,
    /// Actual pair-kernel evaluations; two per edge when the trial dt changes.
    pub fraction_evaluations: u64,
    /// Actual charged logical work.
    pub work_units: u64,
    /// Peak owned vector payload, excluding allocator metadata/rounding.
    pub peak_live_bytes: usize,
    /// Largest normalized local estimate among accepted trials, at most one.
    pub max_accepted_error_ratio: f64,
    /// Admitted work ceiling and actual ambient consumption.
    pub budget: BudgetConsumption,
}

/// Input, resource or time-accuracy refusal; the caller's state stays untouched.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveError {
    /// The existing checked PSE owner refused its physical or execution domain.
    Engine(VpmError),
    /// The local estimate still exceeds tolerance at the minimum retry step.
    ToleranceUnreachable {
        /// Exact rejected step bits.
        dt_bits: u64,
        /// Exact normalized error estimate bits; may represent positive infinity.
        error_ratio_bits: u64,
    },
    /// Floating-point time can no longer advance or a positive retry is lost.
    TimeStagnation {
        /// Last accepted elapsed time bits.
        elapsed_bits: u64,
        /// Attempted positive step bits.
        dt_bits: u64,
    },
}

impl core::fmt::Display for AdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Engine(error) => core::fmt::Display::fmt(error, f),
            Self::ToleranceUnreachable {
                dt_bits,
                error_ratio_bits,
            } => write!(
                f,
                "PSE local time tolerance is unreachable at dt={} s (error ratio={})",
                f64::from_bits(*dt_bits),
                f64::from_bits(*error_ratio_bits),
            ),
            Self::TimeStagnation {
                elapsed_bits,
                dt_bits,
            } => write!(
                f,
                "PSE time cannot advance from {} s by {} s",
                f64::from_bits(*elapsed_bits),
                f64::from_bits(*dt_bits),
            ),
        }
    }
}

impl core::error::Error for AdaptiveError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Engine(error) => Some(error),
            _ => None,
        }
    }
}

impl From<VpmError> for AdaptiveError {
    fn from(error: VpmError) -> Self {
        Self::Engine(error)
    }
}

impl From<BudgetRefusal> for AdaptiveError {
    fn from(error: BudgetRefusal) -> Self {
        Self::Engine(error.into())
    }
}

impl PreparedPse {
    /// Complete a fixed-site diffusion window with automatic positive steps.
    ///
    /// For each particle, acceptance uses
    /// |Gamma_fine-Gamma_coarse|/[3(A*atol + rtol*max(|Gamma_old|,|Gamma_fine|))].
    /// The maximum over particles must be at most one. This is a second-order
    /// step-doubling estimate, not a rigorous global time-error certificate.
    /// No Richardson correction, clipping, negative diffusion time or
    /// postprocessing of the conserved circulation is performed.
    ///
    /// Each trial charges 6 E pair visits and 2 N + 6 E + 1 logical units:
    /// N copies, 6 E exchanges, N error estimates and one trial boundary.
    /// When the trial dt changes, E additional work units prepare both pair
    /// fractions (2 E kernel evaluations); unchanged dt reuses the cache.
    /// Initial validation costs N. Peak
    /// owned payload is 3 N particles + 2 E f64 values. Identity windows
    /// (zero duration, zero viscosity or no edges) own only N particles and
    /// charge N validation units, while still checking every input and cap.
    ///
    /// # Errors
    /// Refuses invalid options or changed prepared geometry, failed original
    /// pair arithmetic, cancellation, exhausted aggregate caps, a tolerance
    /// that needs a smaller-than-admitted step, or unrepresentable progress.
    /// Errors return no partial window and leave the graph/input unchanged.
    #[allow(clippy::too_many_lines)]
    pub fn advance_adaptive_with_cx(
        &self,
        cx: &Cx<'_>,
        particles: &[VortexParticle],
        duration_s: f64,
        options: AdaptiveOptions,
        budget: AdaptiveBudget,
    ) -> Result<AdaptiveRun, AdaptiveError> {
        validate_options(duration_s, options)?;
        let n = self.positions().len();
        if particles.len() != n {
            return Err(VpmError::LengthMismatch {
                field: "prepared particles",
                expected: n,
                actual: particles.len(),
            }
            .into());
        }
        admit_cap("particles", n, budget.max_particles)?;
        let active = duration_s > 0.0
            && self.law().kinematic_viscosity_m2_s() > 0.0
            && !self.edges().is_empty();
        let e = if active { self.edges().len() } else { 0 };
        let n64 = count(n, "adaptive particles")?;
        let e64 = count(e, "adaptive edges")?;
        let pairs_per_trial = e64.checked_mul(6).ok_or(overflow("adaptive pair work"))?;
        let base_trial_work = n64
            .checked_mul(2)
            .and_then(|v| e64.checked_mul(6).and_then(|edges| v.checked_add(edges)))
            .and_then(|v| v.checked_add(1))
            .ok_or(overflow("adaptive trial work"))?;
        admit_count("adaptive work", n64, budget.max_work_units)?;
        let particle_copies = if active { 3 } else { 1 };
        let peak_live_bytes = n
            .checked_mul(core::mem::size_of::<VortexParticle>())
            .and_then(|v| v.checked_mul(particle_copies))
            .and_then(|v| {
                e.checked_mul(core::mem::size_of::<[f64; 2]>())
                    .and_then(|edges| v.checked_add(edges))
            })
            .ok_or(overflow("adaptive live bytes"))?;
        admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
        let planned_work = if active { budget.max_work_units } else { n64 };
        let mut accountant = AdmittedBudget::admit_ambient(cx, planned_work)?;
        accountant.checkpoint("fs-vpm.pse-adaptive-admission", cx)?;
        let mut meter = WorkMeter::default();
        let mut state = reserved_vec("pse-adaptive-state", n)?;
        let mut coarse = reserved_vec("pse-adaptive-coarse", if active { n } else { 0 })?;
        let mut fine = reserved_vec("pse-adaptive-fine", if active { n } else { 0 })?;
        let mut fractions = reserved_vec("pse-adaptive-fractions", e)?;
        for (i, particle) in particles.iter().copied().enumerate() {
            meter.unit(&mut accountant, cx, "fs-vpm.pse-adaptive-input", || {
                validate_particle(particle, i)?;
                for axis in 0..2 {
                    if particle.pos[axis].to_bits() != self.positions()[i][axis].to_bits() {
                        return Err(invalid_scalar(
                            if axis == 0 {
                                "prepared position.x"
                            } else {
                                "prepared position.y"
                            },
                            particle.pos[axis],
                            Some(i),
                        ));
                    }
                }
                state.push(particle);
                Ok(())
            })?;
        }
        meter.boundary(&mut accountant, cx, "fs-vpm.pse-adaptive-input-complete")?;

        let mut elapsed = if active { 0.0 } else { duration_s };
        let mut dt = options.initial_dt_s;
        let mut attempts = 0_usize;
        let mut accepted = 0_usize;
        let mut rejected = 0_usize;
        let mut pair_evaluations = 0_u64;
        let mut fraction_evaluations = 0_u64;
        let mut cached_dt_bits = None;
        let mut expected_work = n64;
        let mut max_accepted_error_ratio = 0.0_f64;

        while elapsed < duration_s {
            let next_attempt = attempts
                .checked_add(1)
                .ok_or(overflow("adaptive attempts"))?;
            admit_cap("adaptive attempts", next_attempt, budget.max_attempts)?;
            let next_pairs = pair_evaluations
                .checked_add(pairs_per_trial)
                .ok_or(overflow("adaptive accumulated pair work"))?;
            admit_count("pair evaluations", next_pairs, budget.max_pair_evaluations)?;
            let remaining = duration_s - elapsed;
            let h = dt.min(remaining);
            let refresh_fractions = cached_dt_bits != Some(h.to_bits());
            let trial_work = base_trial_work
                .checked_add(if refresh_fractions { e64 } else { 0 })
                .ok_or(overflow("adaptive trial work"))?;
            let next_work = expected_work
                .checked_add(trial_work)
                .ok_or(overflow("adaptive accumulated work"))?;
            admit_count("adaptive work", next_work, budget.max_work_units)?;
            let endpoint = if h == remaining {
                duration_s
            } else {
                elapsed + h
            };
            if h <= 0.0 || endpoint <= elapsed {
                return Err(AdaptiveError::TimeStagnation {
                    elapsed_bits: elapsed.to_bits(),
                    dt_bits: h.to_bits(),
                });
            }
            // Buffers retain their reserved capacity between all trials.
            coarse.clear();
            fine.clear();
            for &particle in &state {
                meter.unit(&mut accountant, cx, "fs-vpm.pse-adaptive-copy", || {
                    coarse.push(particle);
                    fine.push(particle);
                    Ok(())
                })?;
            }
            if refresh_fractions {
                fractions.clear();
                let full = PairKernel::new(self.law(), h, 0.5);
                let half = PairKernel::new(self.law(), h, 0.25);
                for &[i, j] in self.edges() {
                    meter.unit(&mut accountant, cx, "fs-vpm.pse-adaptive-fraction", || {
                        let areas = [self.areas_m2()[i], self.areas_m2()[j]];
                        let pair = [
                            full.fraction(state[i], state[j], areas),
                            half.fraction(state[i], state[j], areas),
                        ];
                        if pair
                            .iter()
                            .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
                        {
                            return Err(VpmError::NonFiniteResult {
                                phase: "fs-vpm.pse-adaptive-fraction",
                                particle: i,
                            });
                        }
                        fractions.push(pair);
                        Ok(())
                    })?;
                }
                meter.boundary(
                    &mut accountant,
                    cx,
                    "fs-vpm.pse-adaptive-fractions-complete",
                )?;
                cached_dt_bits = Some(h.to_bits());
                fraction_evaluations = fraction_evaluations
                    .checked_add(e64 * 2)
                    .ok_or(overflow("adaptive fraction evaluations"))?;
            }
            sweep(
                self,
                &mut coarse,
                &fractions,
                0,
                &mut meter,
                &mut accountant,
                cx,
            )?;
            sweep(
                self,
                &mut fine,
                &fractions,
                1,
                &mut meter,
                &mut accountant,
                cx,
            )?;
            sweep(
                self,
                &mut fine,
                &fractions,
                1,
                &mut meter,
                &mut accountant,
                cx,
            )?;
            let mut error_ratio = 0.0_f64;
            for i in 0..n {
                meter.unit(&mut accountant, cx, "fs-vpm.pse-adaptive-error", || {
                    error_ratio = error_ratio.max(normalized_error(
                        state[i].circulation,
                        coarse[i].circulation,
                        fine[i].circulation,
                        self.areas_m2()[i],
                        options,
                    ));
                    Ok(())
                })?;
            }
            meter.unit(&mut accountant, cx, "fs-vpm.pse-adaptive-trial", || Ok(()))?;
            meter.boundary(&mut accountant, cx, "fs-vpm.pse-adaptive-trial-complete")?;
            attempts = next_attempt;
            pair_evaluations = next_pairs;
            expected_work = next_work;
            if error_ratio <= 1.0 {
                core::mem::swap(&mut state, &mut fine);
                elapsed = endpoint;
                accepted += 1;
                max_accepted_error_ratio = max_accepted_error_ratio.max(error_ratio);
                if error_ratio <= 0.125 {
                    dt = (h + h).min(options.max_dt_s).max(options.min_dt_s);
                } else {
                    dt = h.max(options.min_dt_s);
                }
            } else {
                rejected += 1;
                if h <= options.min_dt_s {
                    return Err(AdaptiveError::ToleranceUnreachable {
                        dt_bits: h.to_bits(),
                        error_ratio_bits: error_ratio.to_bits(),
                    });
                }
                dt = (h * 0.5).max(options.min_dt_s);
                if dt >= h {
                    return Err(AdaptiveError::TimeStagnation {
                        elapsed_bits: elapsed.to_bits(),
                        dt_bits: dt.to_bits(),
                    });
                }
            }
        }
        meter.flush(&mut accountant, "fs-vpm.pse-adaptive-finalize")?;
        if meter.completed != expected_work
            || accountant.consumption().cost_charged != expected_work
        {
            return Err(VpmError::AccountingMismatch {
                planned: expected_work,
                completed: meter.completed,
                charged: accountant.consumption().cost_charged,
            }
            .into());
        }
        accountant.checkpoint("fs-vpm.pse-adaptive-publication", cx)?;
        Ok(AdaptiveRun {
            particles: state,
            elapsed_s: duration_s,
            accepted_steps: accepted,
            rejected_steps: rejected,
            attempts,
            pair_evaluations,
            fraction_evaluations,
            work_units: expected_work,
            peak_live_bytes,
            max_accepted_error_ratio,
            budget: accountant.consumption(),
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn sweep(
    graph: &PreparedPse,
    particles: &mut [VortexParticle],
    fractions: &[[f64; 2]],
    coefficient: usize,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    let pairs = graph.edges().iter().zip(fractions);
    for (&[i, j], fraction) in pairs.clone().chain(pairs.rev()) {
        meter.unit(accountant, cx, "fs-vpm.pse-adaptive-pair", || {
            exchange_fraction(particles, graph.areas_m2(), fraction[coefficient], i, j)
        })?;
    }
    meter.boundary(accountant, cx, "fs-vpm.pse-adaptive-sweep-complete")
}

fn normalized_error(old: f64, coarse: f64, fine: f64, area: f64, options: AdaptiveOptions) -> f64 {
    if fine == coarse {
        return 0.0;
    }
    let scale = old.abs().max(coarse.abs()).max(fine.abs());
    let delta = fine - coarse;
    let numerator = if delta.is_finite() {
        (delta / scale).abs() / 3.0
    } else {
        (fine / scale - coarse / scale).abs() / 3.0
    };
    let absolute = area * options.absolute_vorticity_tolerance_s_inv;
    let normalized_absolute = if absolute.is_normal() {
        absolute / scale
    } else {
        det::exp(
            det::ln(area) + det::ln(options.absolute_vorticity_tolerance_s_inv) - det::ln(scale),
        )
    };
    let relative = options.relative_tolerance * (old.abs().max(fine.abs()) / scale);
    let denominator = normalized_absolute + relative;
    let ratio = numerator / denominator;
    if numerator > 0.0 && denominator.is_normal() && ratio.is_finite() {
        return ratio;
    }
    // An extreme old inventory or area can erase both normalized terms.
    // Recover their ratio in log space instead of accepting a NaN estimate.
    let difference = (fine - coarse).abs();
    let log_difference = if difference.is_finite() {
        det::ln(difference)
    } else {
        let field_scale = fine.abs().max(coarse.abs());
        det::ln(field_scale) + det::ln((fine / field_scale - coarse / field_scale).abs())
    };
    let absolute_log = det::ln(area) + det::ln(options.absolute_vorticity_tolerance_s_inv);
    let relative_log = if options.relative_tolerance == 0.0 {
        f64::NEG_INFINITY
    } else {
        det::ln(options.relative_tolerance) + det::ln(old.abs().max(fine.abs()))
    };
    let largest = absolute_log.max(relative_log);
    let denominator_log =
        largest + det::ln(det::exp(absolute_log - largest) + det::exp(relative_log - largest));
    det::exp(log_difference - det::ln(3.0) - denominator_log)
}

fn validate_options(duration: f64, options: AdaptiveOptions) -> Result<(), VpmError> {
    if !duration.is_finite() || duration < 0.0 {
        return Err(invalid_scalar("adaptive duration", duration, None));
    }
    for (field, value) in [
        ("adaptive initial dt", options.initial_dt_s),
        ("adaptive minimum dt", options.min_dt_s),
        ("adaptive maximum dt", options.max_dt_s),
        (
            "adaptive absolute vorticity tolerance",
            options.absolute_vorticity_tolerance_s_inv,
        ),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid_scalar(field, value, None));
        }
    }
    if options.min_dt_s > options.initial_dt_s || options.initial_dt_s > options.max_dt_s {
        return Err(invalid_scalar(
            "adaptive dt ordering",
            options.initial_dt_s,
            None,
        ));
    }
    if !options.relative_tolerance.is_finite() || options.relative_tolerance < 0.0 {
        return Err(invalid_scalar(
            "adaptive relative tolerance",
            options.relative_tolerance,
            None,
        ));
    }
    Ok(())
}

fn overflow(resource: &'static str) -> VpmError {
    VpmError::PlanOverflow { resource }
}

fn count(value: usize, resource: &'static str) -> Result<u64, VpmError> {
    u64::try_from(value).map_err(|_| overflow(resource))
}

fn admit_count(resource: &'static str, required: u64, maximum: u64) -> Result<(), VpmError> {
    if required > maximum {
        return Err(VpmError::CapExceeded {
            resource,
            required: required.into(),
            maximum: maximum.into(),
        });
    }
    Ok(())
}
