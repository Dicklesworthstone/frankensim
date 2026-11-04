//! Continuous native forcing on the existing held-interface coupling owner.
use super::{IntegrationMethod, MultirateError, MultiratePair, SourceInput, check, finite};
use crate::linear_exact::{DenseExponentialBudget, LinearSplitBudget, PreparedExactTones};
use crate::{PhsError, PortHamiltonian, PreparedStepRecord};
use fs_math::det;

/// Cold preparation of one original subsystem's external source family.
/// Frequencies and envelopes do not change the original physical operators.
#[derive(Clone, Copy, Debug)]
pub struct ContinuousSourcePlan<'a> {
    /// Angular frequencies [rad/s], in caller coefficient order; DC is allowed.
    pub angular_frequencies_rad_s: &'a [f64],
    /// Permit linearly changing cosine/sine envelopes across a common window.
    /// Constant preparation keeps the original smaller two-coordinate source.
    pub affine_envelopes: bool,
    /// Complete source preparation, including the additional held interface
    /// coordinate and this adapter's reusable buffers, must fit this budget.
    pub budget: DenseExponentialBudget,
}

pub(super) struct ContinuousClock<'a> {
    source: PreparedExactTones<'a>,
    tones: usize,
    substeps: usize,
    phases: Vec<(f64, f64)>,
    physical: [Vec<f64>; 4],
    external: [Vec<f64>; 4],
    work_input: Vec<f64>,
}

fn bad(what: &'static str) -> PhsError {
    PhsError::Dimension { what }
}

impl<'a> ContinuousClock<'a> {
    fn new(
        system: &'a PortHamiltonian,
        substeps: usize,
        dt: f64,
        plan: ContinuousSourcePlan<'_>,
    ) -> Result<Self, PhsError> {
        let tones = plan.angular_frequencies_rad_s.len();
        let count = tones
            .checked_add(1)
            .and_then(|n| n.checked_mul(system.port_dim()))
            .ok_or_else(|| bad("continuous multirate source extent"))?;
        let phase_count = tones
            .checked_mul(substeps + 1)
            .ok_or_else(|| bad("continuous multirate phase extent"))?;
        let stride: usize = if plan.affine_envelopes { 4 } else { 2 };
        let bytes =
            (count as u128 * (8 + stride) as u128 + 2 * phase_count as u128 + tones as u128 + 1)
                * 8;
        if tones == 0 || bytes > plan.budget.max_dense_bytes as u128 {
            return Err(bad(
                "continuous multirate source buffers exceed explicit budget",
            ));
        }
        let budget = DenseExponentialBudget {
            max_dense_bytes: plan.budget.max_dense_bytes - bytes as usize,
            ..plan.budget
        };
        // The first DC coefficient is the held interface input. It is never
        // included in the independently integrated external-source work.
        let frequencies: Vec<_> = std::iter::once(0.0)
            .chain(plan.angular_frequencies_rad_s.iter().copied())
            .collect();
        let split = LinearSplitBudget {
            max_modes: system.state_dim(),
            max_omega_dt: 0.04,
        };
        let source = if plan.affine_envelopes {
            PreparedExactTones::new_affine(system, &frequencies, dt, split, budget)?
        } else {
            PreparedExactTones::new(system, &frequencies, dt, split, budget)?
        };
        let mut phases = Vec::with_capacity(phase_count);
        for step in 0..=substeps {
            for &omega in plan.angular_frequencies_rad_s {
                let angle = omega * (dt * step as f64);
                if !angle.is_finite() {
                    return Err(bad("continuous multirate source phase overflow"));
                }
                phases.push(det::sin_cos(angle));
            }
        }
        Ok(Self {
            source,
            tones,
            substeps,
            phases,
            physical: std::array::from_fn(|_| vec![0.0; count]),
            external: std::array::from_fn(|_| vec![0.0; count]),
            work_input: vec![0.0; stride * count],
        })
    }

    pub(super) fn valid(&self, input: SourceInput<'_>) -> bool {
        let count = self.tones * self.source.port_dim();
        let valid = |values: &[f64]| values.len() == count && finite(values);
        match input {
            SourceInput::Tones { cosine, sine } => valid(cosine) && valid(sine),
            SourceInput::AffineTones {
                cosine_start,
                sine_start,
                cosine_end,
                sine_end,
            } => {
                self.source.is_affine()
                    && [cosine_start, sine_start, cosine_end, sine_end]
                        .into_iter()
                        .all(valid)
            }
            _ => false,
        }
    }

    fn rotate(&self, step: usize, tone: usize, cosine: f64, sine: f64) -> (f64, f64) {
        let (sn, cs) = self.phases[step * self.tones + tone];
        (cs * cosine + sn * sine, cs * sine - sn * cosine)
    }

    fn phase_after_window(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), PhsError> {
        if tone >= self.tones || !finite(&[cosine, sine]) {
            return Err(bad("continuous multirate phase index or coefficients"));
        }
        let next = self.rotate(self.substeps, tone, cosine, sine);
        if !finite(&[next.0, next.1]) {
            return Err(bad("continuous multirate phase overflow"));
        }
        Ok(next)
    }

    pub(super) fn step(
        &mut self,
        initial: &[f64],
        input: SourceInput<'_>,
        exchange: &[f64],
        step: usize,
        next: &mut [f64],
        output: &mut [f64],
    ) -> Result<(PreparedStepRecord, f64), PhsError> {
        let (c0, s0, c1, s1) = match input {
            SourceInput::Tones { cosine, sine } => (cosine, sine, cosine, sine),
            SourceInput::AffineTones {
                cosine_start,
                sine_start,
                cosine_end,
                sine_end,
            } => (cosine_start, sine_start, cosine_end, sine_end),
            _ => return Err(bad("continuous multirate source variant")),
        };
        let fraction0 = step as f64 / self.substeps as f64;
        let fraction1 = (step + 1) as f64 / self.substeps as f64;
        let width = self.tones + 1;
        for (port, &coupled) in exchange.iter().enumerate() {
            let base = port * width;
            self.physical[0][base] = coupled;
            self.physical[2][base] = coupled;
            for tone in 0..self.tones {
                let at = port * self.tones + tone;
                // Interpolate the physical envelope, then express BOTH native
                // endpoints in that interval's start carrier phase.
                let start = self.rotate(
                    step,
                    tone,
                    blend(c0[at], c1[at], fraction0),
                    blend(s0[at], s1[at], fraction0),
                );
                let end = self.rotate(
                    step,
                    tone,
                    blend(c0[at], c1[at], fraction1),
                    blend(s0[at], s1[at], fraction1),
                );
                for (field, value) in [start.0, start.1, end.0, end.1].into_iter().enumerate() {
                    self.physical[field][base + tone + 1] = value;
                    self.external[field][base + tone + 1] = value;
                }
            }
        }
        self.source.step_with_source_work_into(
            initial,
            &self.physical[0],
            &self.physical[1],
            self.source
                .is_affine()
                .then_some((&self.physical[2], &self.physical[3])),
            &self.external[0],
            &self.external[1],
            self.source
                .is_affine()
                .then_some((&self.external[2], &self.external[3])),
            &mut self.work_input,
            next,
            output,
        )
    }
}

fn blend(start: f64, end: f64, fraction: f64) -> f64 {
    if fraction == 0.0 || start == end {
        start
    } else if fraction == 1.0 {
        end
    } else {
        (1.0 - fraction) * start + fraction * end
    }
}

impl<'a> MultiratePair<'a> {
    /// Prepare continuous forcing on either or both original exact-quadratic
    /// owners. The existing interface LU is reused: additive forcing changes its
    /// right-hand side, not the held-connection response. A nonlinear side may
    /// remain unconfigured and use its original held/substep source path.
    ///
    /// This consumes only numerical preparation, never physical state. A side
    /// configured for discrete gradients refuses a continuous plan rather than
    /// silently replacing its requested integrator. `None` preserves that side.
    pub fn with_continuous_sources(
        self,
        left: Option<ContinuousSourcePlan<'_>>,
        right: Option<ContinuousSourcePlan<'_>>,
    ) -> Result<Self, MultirateError> {
        self.with_continuous_sources_controlled(left, right, || false)
    }

    /// Poll around the bounded original source preparations and before publishing
    /// either cache. Their synchronous dense kernels are not preemptible.
    pub fn with_continuous_sources_controlled(
        mut self,
        left: Option<ContinuousSourcePlan<'_>>,
        right: Option<ContinuousSourcePlan<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MultirateError> {
        let mut candidates = [None, None];
        for (index, (owner, plan)) in [(&self.left, left), (&self.right, right)]
            .into_iter()
            .enumerate()
        {
            check(&mut cancelled)?;
            if let Some(plan) = plan {
                if owner.clock.method != IntegrationMethod::ExactQuadratic {
                    return Err(MultirateError::Invalid(
                        "continuous sources require the declared exact quadratic integrator",
                    ));
                }
                candidates[index] = Some(
                    ContinuousClock::new(owner.system, owner.clock.substeps, owner.dt, plan)
                        .map_err(|error| MultirateError::Subsystem {
                            side: owner.side,
                            error,
                        })?,
                );
            }
        }
        check(&mut cancelled)?;
        let [left, right] = candidates;
        if left.is_some() {
            self.left.continuous = left;
        }
        if right.is_some() {
            self.right.continuous = right;
        }
        Ok(self)
    }

    /// Rotate a left source's final envelope to the next common window's start
    /// phase. Call only after accepting a step; no source state is stored here.
    pub fn advance_left_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), MultirateError> {
        self.left
            .continuous
            .as_ref()
            .ok_or(MultirateError::Invalid(
                "left continuous source is not prepared",
            ))?
            .phase_after_window(tone, cosine, sine)
            .map_err(|error| MultirateError::Subsystem {
                side: "left",
                error,
            })
    }

    /// Right-source counterpart of [`Self::advance_left_force_phase`].
    pub fn advance_right_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), MultirateError> {
        self.right
            .continuous
            .as_ref()
            .ok_or(MultirateError::Invalid(
                "right continuous source is not prepared",
            ))?
            .phase_after_window(tone, cosine, sine)
            .map_err(|error| MultirateError::Subsystem {
                side: "right",
                error,
            })
    }
}
