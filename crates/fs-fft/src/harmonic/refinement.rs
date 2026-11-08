//! Shared bounded variable projection for stationary tones, poles and chirps.

use super::{
    C64, ChirpedTone, ChirpedToneBounds, ChirpedToneFitPlan, DampedToneBounds, DampedToneFitPlan,
    DampedTonePole, FitDesign, HarmonicFit, HarmonicFitBudget, HarmonicFitError, ToneFitPlan,
    sample_sin_cos,
};
use fs_la::factor::{Qr, qr};
use fs_math::{dd::Dd, det};

pub(super) enum Domain<'a> {
    Stationary(&'a [f64], &'a [[f64; 2]]),
    Damped(&'a [DampedTonePole], &'a [DampedToneBounds]),
    Chirped(&'a [ChirpedTone], &'a [ChirpedToneBounds]),
}

impl Domain<'_> {
    fn tones(&self) -> usize {
        match self {
            Self::Stationary(seeds, _) => seeds.len(),
            Self::Damped(seeds, _) => seeds.len(),
            Self::Chirped(seeds, _) => seeds.len(),
        }
    }

    fn damped(&self) -> bool {
        matches!(self, Self::Damped(..))
    }

    fn paired(&self) -> bool {
        !matches!(self, Self::Stationary(..))
    }

    fn bound(&self, coordinate: usize) -> [f64; 2] {
        match self {
            Self::Stationary(_, bounds) => bounds[coordinate],
            Self::Damped(_, bounds) if coordinate % 2 == 0 => bounds[coordinate / 2].frequency_hz,
            Self::Damped(_, bounds) => bounds[coordinate / 2].decay_rate_per_s,
            Self::Chirped(_, bounds) if coordinate % 2 == 0 => {
                bounds[coordinate / 2].start_frequency_hz
            }
            Self::Chirped(_, bounds) => bounds[coordinate / 2].end_frequency_hz,
        }
    }

    fn initial_plan(&self, n: usize, dt: f64, budget: HarmonicFitBudget) -> Result<Plan, Error> {
        Ok(match self {
            Self::Stationary(seeds, _) => Plan::Stationary(ToneFitPlan::new(n, dt, seeds, budget)?),
            Self::Damped(seeds, _) => Plan::Damped(DampedToneFitPlan::new(n, dt, seeds, budget)?),
            Self::Chirped(seeds, _) => {
                Plan::Chirped(ChirpedToneFitPlan::new(n, dt, seeds, budget)?)
            }
        })
    }

    fn trial_plan(
        &self,
        n: usize,
        dt: f64,
        parameters: &[f64],
        budget: HarmonicFitBudget,
    ) -> Result<Plan, Error> {
        Ok(match self {
            Self::Stationary(..) => Plan::Stationary(ToneFitPlan::new(n, dt, parameters, budget)?),
            Self::Damped(..) => {
                let poles: Vec<_> = parameters
                    .chunks_exact(2)
                    .map(|p| DampedTonePole {
                        frequency_hz: p[0],
                        decay_rate_per_s: p[1],
                    })
                    .collect();
                Plan::Damped(DampedToneFitPlan::new(n, dt, &poles, budget)?)
            }
            Self::Chirped(..) => {
                let tones: Vec<_> = parameters
                    .chunks_exact(2)
                    .map(|p| ChirpedTone {
                        start_frequency_hz: p[0],
                        end_frequency_hz: p[1],
                    })
                    .collect();
                Plan::Chirped(ChirpedToneFitPlan::new(n, dt, &tones, budget)?)
            }
        })
    }
}

#[derive(Clone, Copy)]
pub(super) struct Controls {
    pub tolerances: [f64; 2],
    pub max_evaluations: usize,
    pub min_reciprocal_condition: f64,
    pub fit: HarmonicFitBudget,
}

pub(super) struct Report {
    pub parameters: Vec<f64>,
    pub fits: Vec<HarmonicFit>,
    pub evaluations: usize,
    pub corrections: Vec<f64>,
    pub reciprocal_condition: f64,
}

pub(super) enum Error {
    Fit(HarmonicFitError),
    InvalidParameter {
        what: &'static str,
    },
    Unidentifiable {
        reciprocal_condition: f64,
        minimum: f64,
    },
    NotConverged {
        parameters: Vec<f64>,
        residual_rms: f64,
        evaluations: usize,
        corrections: Vec<f64>,
    },
}

impl From<HarmonicFitError> for Error {
    fn from(error: HarmonicFitError) -> Self {
        Self::Fit(error)
    }
}

fn invalid(what: &'static str) -> Error {
    Error::InvalidParameter { what }
}

fn unidentifiable(condition: f64, controls: Controls) -> Error {
    Error::Unidentifiable {
        reciprocal_condition: condition,
        minimum: controls.min_reciprocal_condition,
    }
}

struct Candidate {
    parameters: Vec<f64>,
    fits: Vec<HarmonicFit>,
    residual_rms: f64,
    residual_squared: f64,
    corrections: Vec<f64>,
    reciprocal_condition: f64,
}

impl Candidate {
    fn not_converged(self, evaluations: usize) -> Error {
        Error::NotConverged {
            parameters: self.parameters,
            residual_rms: self.residual_rms,
            evaluations,
            corrections: self.corrections,
        }
    }
}

// The public owners validate their own units, seeds, intervals and tolerances.
// Clock, data, storage, QR folding and line search remain one numerical owner.
pub(super) fn refine_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    domain: Domain<'_>,
    controls: Controls,
) -> Result<Report, Error> {
    let Some(samples) = channels.first() else {
        return Err(invalid("at least one simultaneous channel is required"));
    };
    let n = samples.len();
    let tones = domain.tones();
    let parameters = tones
        .checked_mul(if domain.paired() { 2 } else { 1 })
        .ok_or_else(|| invalid("tone parameter count overflows"))?;
    let columns = tones
        .checked_mul(2)
        .and_then(|k| k.checked_add(1))
        .ok_or_else(|| invalid("tone design size overflows"))?;
    if n <= columns || channels.iter().any(|channel| channel.len() != n) {
        return Err(invalid(
            "equal-length channels need more than 1 + 2 * tone_count samples",
        ));
    }
    if (n - columns)
        .checked_mul(channels.len())
        .is_none_or(|rows| rows < parameters)
    {
        return Err(invalid(
            "insufficient joint residual dimensions to identify all tone parameters",
        ));
    }
    let duration = (n - 1) as f64 * timestep_s;
    if !(duration.is_finite() && duration > 0.0) {
        return Err(invalid("record duration must be finite and positive"));
    }
    if channels.iter().any(|c| c.iter().any(|x| !x.is_finite())) {
        return Err(HarmonicFitError::NonFinite {
            what: "input sample",
        }
        .into());
    }
    let width = parameters + 1;
    let rows = if domain.paired() { n } else { n - columns }.max(width);
    let extra_bytes = refinement_bytes(n, parameters, channels.len(), rows, domain.paired())?;
    let plan_extra = if domain.paired() {
        tones.checked_mul(core::mem::size_of::<f64>())
    } else {
        Some(0)
    }
    .ok_or_else(|| invalid("paired tone parameter storage size overflows"))?;
    let required_bytes = super::fit_workspace_bytes(n, tones)?
        .checked_add(plan_extra)
        .and_then(|v| v.checked_add(extra_bytes))
        .ok_or_else(|| invalid("tone refinement workspace size overflows"))?;
    if required_bytes > controls.fit.max_workspace_bytes {
        return Err(HarmonicFitError::BudgetExceeded {
            required_bytes,
            max_bytes: controls.fit.max_workspace_bytes,
        }
        .into());
    }
    let plan_budget = HarmonicFitBudget {
        max_workspace_bytes: controls.fit.max_workspace_bytes - extra_bytes,
        ..controls.fit
    };
    let plan = domain.initial_plan(n, timestep_s, plan_budget)?;
    let mut ratios: Vec<f64> = channels
        .iter()
        .map(|channel| {
            let lo = channel.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = channel.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let center = lo.midpoint(hi);
            channel
                .iter()
                .map(|x| (x - center).abs())
                .fold(0.0, f64::max)
        })
        .collect();
    let physical_scale = ratios.iter().copied().fold(0.0, f64::max);
    if physical_scale == 0.0 {
        return Err(unidentifiable(0.0, controls));
    }
    for ratio in &mut ratios {
        *ratio /= physical_scale;
    }
    let mut workspace = Workspace {
        scratch: vec![0.0; n],
        derivative: vec![0.0; n],
        residual: vec![0.0; if domain.paired() { n } else { 0 }],
        augmented: vec![0.0; rows * width],
        folded: vec![0.0; 2 * width * width],
        projected_folded: vec![
            0.0;
            if domain.paired() {
                2 * width * width
            } else {
                0
            }
        ],
    };
    let mut current = evaluate(
        channels,
        &ratios,
        physical_scale,
        plan,
        duration,
        controls,
        &mut workspace,
    )?;
    let mut evaluations = 1;
    loop {
        if current
            .corrections
            .iter()
            .enumerate()
            .all(|(i, d)| d.abs() <= controls.tolerances[usize::from(domain.paired()) * (i % 2)])
        {
            return Ok(Report {
                parameters: current.parameters,
                fits: current.fits,
                evaluations,
                corrections: current.corrections,
                reciprocal_condition: current.reciprocal_condition,
            });
        }
        // Preserve the stationary quarter-cycle limit. For chirps both
        // endpoint phase weights are nonnegative and sum to t/duration, so
        // this cap bounds their combined phase change too. Damped trials also
        // limit the amplitude envelope change to half an e-fold per record.
        let mut fraction = 1.0_f64;
        for (i, d) in current.corrections.iter().enumerate() {
            let cap = if domain.damped() && i % 2 == 1 {
                0.5
            } else {
                0.25
            };
            fraction = fraction.min(cap / duration / d.abs());
        }
        loop {
            if evaluations == controls.max_evaluations {
                return Err(current.not_converged(evaluations));
            }
            let trial_parameters: Vec<f64> = current
                .parameters
                .iter()
                .zip(&current.corrections)
                .enumerate()
                .map(|(i, (&x, &d))| {
                    let [lo, hi] = domain.bound(i);
                    (x + fraction * d).clamp(lo, hi)
                })
                .collect();
            if trial_parameters == current.parameters {
                return Err(current.not_converged(evaluations));
            }
            evaluations += 1;
            let trial = domain
                .trial_plan(n, timestep_s, &trial_parameters, plan_budget)
                .and_then(|plan| {
                    evaluate(
                        channels,
                        &ratios,
                        physical_scale,
                        plan,
                        duration,
                        controls,
                        &mut workspace,
                    )
                });
            match trial {
                Ok(candidate) if candidate.residual_squared < current.residual_squared => {
                    current = candidate;
                    break;
                }
                Ok(_)
                | Err(Error::Fit(HarmonicFitError::IllConditioned { .. }))
                | Err(Error::Unidentifiable { .. }) => fraction *= 0.5,
                Err(error) => return Err(error),
            }
        }
    }
}

// Includes both candidates, previous/current small QR, and full-Jacobian
// residual reconstruction. Using parameter count in the report allowance
// conservatively covers paired-parameter conversion temporaries before allocation.
fn refinement_bytes(
    n: usize,
    k: usize,
    channels: usize,
    rows: usize,
    paired: bool,
) -> Result<usize, Error> {
    let s = k + 1;
    let scalars = n
        .checked_mul(if paired { 3 } else { 2 })
        .and_then(|v| rows.checked_mul(s)?.checked_mul(2)?.checked_add(v))
        .and_then(|v| {
            s.checked_mul(s)?
                .checked_mul(if paired { 12 } else { 6 })?
                .checked_add(v)
        })
        .and_then(|v| s.checked_mul(3)?.checked_add(v))
        .and_then(|v| k.checked_mul(if paired { 12 } else { 8 })?.checked_add(v))
        .and_then(|v| {
            k.checked_mul(4)?
                .checked_add(1)?
                .checked_mul(channels)?
                .checked_add(v)
        })
        .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()));
    scalars
        .and_then(|v| {
            channels
                .checked_mul(2 * core::mem::size_of::<HarmonicFit>())?
                .checked_add(v)
        })
        .ok_or_else(|| invalid("tone refinement workspace size overflows"))
}

struct Workspace {
    scratch: Vec<f64>,
    derivative: Vec<f64>,
    residual: Vec<f64>,
    augmented: Vec<f64>,
    folded: Vec<f64>,
    projected_folded: Vec<f64>,
}

enum Plan {
    Stationary(ToneFitPlan),
    Damped(DampedToneFitPlan),
    Chirped(ChirpedToneFitPlan),
}

impl Plan {
    fn design(&self) -> &FitDesign {
        match self {
            Self::Stationary(plan) => &plan.design,
            Self::Damped(plan) => &plan.design,
            Self::Chirped(plan) => &plan.design,
        }
    }

    fn into_parameters(self) -> Vec<f64> {
        match self {
            Self::Stationary(plan) => plan.frequencies_hz,
            Self::Damped(plan) => plan
                .poles
                .iter()
                .flat_map(|p| [p.frequency_hz, p.decay_rate_per_s])
                .collect(),
            Self::Chirped(plan) => plan
                .tones
                .iter()
                .flat_map(|p| [p.start_frequency_hz, p.end_frequency_hz])
                .collect(),
        }
    }

    // Derivative of the fitted model with respect to parameter * duration,
    // in Q coordinates. Stationary calls retain the original Kaufman columns.
    // For damped and chirped calls the complete variable-projection derivative is
    // P_perp D c + (A^+)^T D^T r. A=Q1 R gives [R^-T D^T r; Q2^T D c],
    // so neither an explicit projector nor normal equations are required.
    fn derivative(
        &self,
        coordinate: usize,
        coefficients: &[f64],
        residual: &[f64],
        out: &mut [f64],
    ) -> f64 {
        let design = self.design();
        let paired = !matches!(self, Self::Stationary(_));
        let tone = if paired { coordinate / 2 } else { coordinate };
        let (frequency, decay) = match self {
            Self::Stationary(plan) => (plan.frequencies_hz[tone], 0.0),
            Self::Damped(plan) => {
                let p = plan.poles[tone];
                (p.frequency_hz, p.decay_rate_per_s)
            }
            Self::Chirped(plan) => (plan.tones[tone].start_frequency_hz, 0.0),
        };
        let a = coefficients[2 * tone + 1];
        let b = coefficients[2 * tone + 2];
        let inverse_intervals = if matches!(self, Self::Chirped(_)) {
            Dd::ONE / Dd::from_f64((design.samples - 1) as f64)
        } else {
            Dd::ZERO
        };
        let (mut cos_residual, mut sin_residual, mut raw_scale) = (0.0, 0.0, 0.0_f64);
        for (j, value) in out.iter_mut().enumerate() {
            // Shifting t by a constant adds only columns already in A. Both
            // P_perp A=0 and A^T r=0 make this the same profiled derivative.
            let time = j as f64 / (design.samples - 1) as f64 - 0.5;
            let (s, c) = if let Self::Chirped(plan) = self {
                super::chirp::sample_sin_cos(
                    plan.tones[tone],
                    design.timestep_s,
                    inverse_intervals,
                    j,
                )
            } else {
                sample_sin_cos(frequency, design.timestep_s, j)
            };
            if paired {
                let (dc, ds) = if matches!(self, Self::Chirped(_)) {
                    let u = j as f64 / (design.samples - 1) as f64;
                    // Center each integrated-frequency weight by its value at
                    // t=T/2. The removed constant is a fitted basis column.
                    let weight = if coordinate % 2 == 0 {
                        u - 0.5 * u * u - 0.375
                    } else {
                        0.5 * u * u - 0.125
                    };
                    let time = core::f64::consts::TAU * weight;
                    (-time * s, time * c)
                } else {
                    let envelope = if j == 0 || decay == 0.0 {
                        1.0
                    } else {
                        det::exp(-((decay * design.timestep_s) * j as f64))
                    };
                    let time = time * envelope;
                    if coordinate % 2 == 0 {
                        let time = core::f64::consts::TAU * time;
                        (-time * s, time * c)
                    } else {
                        (-time * c, -time * s)
                    }
                };
                *value = a.mul_add(dc, b * ds);
                cos_residual = dc.mul_add(residual[j], cos_residual);
                sin_residual = ds.mul_add(residual[j], sin_residual);
            } else {
                *value = core::f64::consts::TAU * time * (-a).mul_add(s, b * c);
            }
            raw_scale = raw_scale.max(value.abs());
        }
        design.factor.apply_qt(out);
        if paired {
            let columns = 2 * design.tones + 1;
            out[..columns].fill(0.0);
            out[2 * tone + 1] = cos_residual;
            out[2 * tone + 2] = sin_residual;
            for i in 0..columns {
                let mut value = out[i];
                for (j, &entry) in out[..i].iter().enumerate() {
                    value = (-design.factor.r(j, i)).mul_add(entry, value);
                }
                out[i] = value / design.factor.r(i, i);
            }
        }
        raw_scale
    }
}

// Each channel's QR([D, residual]) compresses the full least-squares geometry.
// QR([R_previous; R_channel]) preserves common-unit channel weights and avoids
// both normal equations and a C*N-row derivative allocation.
fn fold_channel(augmented: &[f64], rows: usize, width: usize, folded: &mut [f64]) -> Qr {
    let local = qr(augmented, rows, width);
    for i in 0..width {
        for j in 0..width {
            folded[(i + width) * width + j] = local.r(i, j);
        }
    }
    drop(local);
    let combined = qr(folded, 2 * width, width);
    for i in 0..width {
        for j in 0..width {
            folded[i * width + j] = combined.r(i, j);
        }
    }
    combined
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    channels: &[&[f64]],
    ratios: &[f64],
    physical_scale: f64,
    plan: Plan,
    duration: f64,
    controls: Controls,
    work: &mut Workspace,
) -> Result<Candidate, Error> {
    let design = plan.design();
    let tones = design.tones;
    let columns = 2 * tones + 1;
    let paired = !matches!(plan, Plan::Stationary(_));
    let parameters = if paired { 2 * tones } else { tones };
    let offset = if paired { 0 } else { columns };
    let width = parameters + 1;
    let rows = work.augmented.len() / width;
    work.folded.fill(0.0);
    work.projected_folded.fill(0.0);
    let mut factor: Option<Qr> = None;
    let mut projected_factor: Option<Qr> = None;
    let mut fits = Vec::with_capacity(channels.len());
    let (mut residual_squared, mut joint_scale, mut raw_scale) = (0.0, 0.0_f64, 0.0_f64);
    let mut projected_scale = 0.0_f64;
    for (&samples, &ratio) in channels.iter().zip(ratios) {
        let mut phasors = vec![C64::default(); tones];
        let summary = design.fit_into(samples, &mut work.scratch, &mut phasors)?;
        fits.push(HarmonicFit {
            dc: summary.dc,
            phasors,
            residual_rms: summary.residual_rms,
        });
        if ratio == 0.0 {
            continue;
        }
        work.augmented.fill(0.0);
        for (j, &residual) in work.scratch[columns..].iter().enumerate() {
            let weighted = residual * ratio;
            work.augmented[(j + columns - offset) * width + parameters] = weighted;
            residual_squared = weighted.mul_add(weighted, residual_squared);
        }
        if paired {
            work.residual[..columns].fill(0.0);
            work.residual[columns..].copy_from_slice(&work.scratch[columns..]);
            design.factor.apply_q(&mut work.residual);
        }
        let mut channel_scale = 0.0_f64;
        for parameter in 0..parameters {
            let raw = plan.derivative(
                parameter,
                &work.scratch[..columns],
                &work.residual,
                &mut work.derivative,
            );
            raw_scale = raw_scale.max(raw * ratio);
            for (j, &derivative) in work.derivative[offset..].iter().enumerate() {
                let weighted = derivative * ratio;
                if !weighted.is_finite() {
                    return Err(HarmonicFitError::NonFinite {
                        what: "profiled tone derivative",
                    }
                    .into());
                }
                work.augmented[j * width + parameter] = weighted;
                channel_scale = channel_scale.max(weighted.abs());
                if paired && j >= columns {
                    projected_scale = projected_scale.max(weighted.abs());
                }
            }
        }
        if channel_scale > joint_scale {
            let change = joint_scale / channel_scale;
            for row in work.folded[..width * width].chunks_exact_mut(width) {
                for value in &mut row[..parameters] {
                    *value *= change;
                }
            }
            if paired {
                for row in work.projected_folded[..width * width].chunks_exact_mut(width) {
                    for value in &mut row[..parameters] {
                        *value *= change;
                    }
                }
            }
            joint_scale = channel_scale;
        }
        if joint_scale == 0.0 {
            continue;
        }
        for row in work.augmented.chunks_exact_mut(width) {
            for value in &mut row[..parameters] {
                *value /= joint_scale;
            }
        }
        factor = Some(fold_channel(&work.augmented, rows, width, &mut work.folded));
        if paired {
            // The residual-dependent term can remain nonzero for an absent
            // pole. It must improve the Newton direction without inventing
            // observability: also gate the projected physical sensitivity
            // P_perp*D*c. Reuse the same buffer, zeroing only its Q1 rows.
            work.augmented[..columns * width].fill(0.0);
            projected_factor = Some(fold_channel(
                &work.augmented,
                rows,
                width,
                &mut work.projected_folded,
            ));
        }
    }
    if joint_scale <= 64.0 * f64::EPSILON * raw_scale {
        return Err(unidentifiable(0.0, controls));
    }
    // Condition is scale invariant: an entirely unexcited pole set can leave
    // a well-conditioned but uniformly roundoff-sized P_perp*D*c, while D^T*r
    // supplies a large full-Jacobian block. Demand resolved sensitivity in the
    // original globally normalized signal units as well. Unit changes cancel
    // in physical_scale/ratios, so physically weak but resolved records remain
    // admissible. The raw scale also protects cancellation-heavy projections.
    if paired && projected_scale <= 64.0 * f64::EPSILON * raw_scale.max(1.0) {
        return Err(unidentifiable(0.0, controls));
    }
    let factor = factor.ok_or_else(|| unidentifiable(0.0, controls))?;
    let mut condition = super::triangular_reciprocal_condition(&factor, parameters);
    if !condition.is_finite() {
        return Err(unidentifiable(condition, controls));
    }
    if let Some(projected) = projected_factor {
        let projected_condition = super::triangular_reciprocal_condition(&projected, parameters);
        if !projected_condition.is_finite() {
            return Err(unidentifiable(projected_condition, controls));
        }
        condition = condition.min(projected_condition);
    }
    if !condition.is_finite() || condition < controls.min_reciprocal_condition {
        return Err(unidentifiable(condition, controls));
    }
    let mut correction: Vec<f64> = (0..parameters).map(|i| factor.r(i, parameters)).collect();
    for i in (0..parameters).rev() {
        let mut value = correction[i];
        for (j, &x) in correction.iter().enumerate().skip(i + 1) {
            value = (-factor.r(i, j)).mul_add(x, value);
        }
        correction[i] = value / factor.r(i, i);
    }
    for value in &mut correction {
        *value = (*value / joint_scale) / duration;
    }
    let residual_rms =
        det::sqrt(residual_squared / work.scratch.len() as f64 / channels.len() as f64)
            * physical_scale;
    if correction.iter().any(|x| !x.is_finite())
        || !residual_squared.is_finite()
        || !residual_rms.is_finite()
    {
        return Err(HarmonicFitError::NonFinite {
            what: "joint tone parameter correction",
        }
        .into());
    }
    Ok(Candidate {
        parameters: plan.into_parameters(),
        fits,
        residual_rms,
        residual_squared,
        corrections: correction,
        reciprocal_condition: condition,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn signal(p: DampedTonePole, a: C64, t: f64) -> f64 {
        let (s, c) = (core::f64::consts::TAU * p.frequency_hz * t).sin_cos();
        (-p.decay_rate_per_s * t).exp() * (a.re * c + a.im * s)
    }

    #[test]
    fn full_profiled_ringdown_jacobian_matches_refits_away_from_a_zero_residual() {
        // G1: a detuned record with a substantial omitted component. The
        // residual-dependent term is nonzero, so an amplitude-held derivative
        // or the Kaufman approximation cannot pass this finite-difference test.
        let (n, dt) = (401, 1.0 / 2048.0);
        let duration = (n - 1) as f64 * dt;
        let samples: Vec<_> = (0..n)
            .map(|j| {
                let t = j as f64 * dt;
                0.8 + signal(
                    DampedTonePole {
                        frequency_hz: 174.0,
                        decay_rate_per_s: 5.3,
                    },
                    C64::new(1.0, 0.5),
                    t,
                ) + signal(
                    DampedTonePole {
                        frequency_hz: 177.5,
                        decay_rate_per_s: 0.3,
                    },
                    C64::new(-0.4, 0.2),
                    t,
                )
            })
            .collect();
        let seed = DampedTonePole {
            frequency_hz: 173.2,
            decay_rate_per_s: 4.1,
        };
        let budget = HarmonicFitBudget {
            max_workspace_bytes: 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        };
        let plan = Plan::Damped(DampedToneFitPlan::new(n, dt, &[seed], budget).unwrap());
        let mut scratch = vec![0.0; n];
        plan.design()
            .fit_into(&samples, &mut scratch, &mut [C64::default()])
            .unwrap();
        let mut residual = vec![0.0; n];
        residual[3..].copy_from_slice(&scratch[3..]);
        plan.design().factor.apply_q(&mut residual);
        let lo = samples.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let center = lo.midpoint(hi);
        let scale = samples
            .iter()
            .map(|x| (x - center).abs())
            .fold(0.0, f64::max);
        for coordinate in 0..2 {
            let mut derivative = vec![0.0; n];
            plan.derivative(coordinate, &scratch[..3], &residual, &mut derivative);
            let span_norm: f64 = derivative[..3].iter().map(|x| x * x).sum();
            let total_norm: f64 = derivative.iter().map(|x| x * x).sum();
            assert!(
                span_norm > 1e-4 * total_norm,
                "residual-dependent term was not exercised"
            );
            plan.design().factor.apply_q(&mut derivative);
            let epsilon = 1e-4;
            let mut sides = Vec::new();
            for sign in [-1.0, 1.0] {
                let mut p = seed;
                if coordinate == 0 {
                    p.frequency_hz += sign * epsilon;
                } else {
                    p.decay_rate_per_s += sign * epsilon;
                }
                let fit = DampedToneFitPlan::new(n, dt, &[p], budget)
                    .unwrap()
                    .fit(&samples)
                    .unwrap();
                sides.push(
                    (0..n)
                        .map(|j| fit.dc + signal(p, fit.phasors[0], j as f64 * dt))
                        .collect::<Vec<_>>(),
                );
            }
            let (mut error, mut reference) = (0.0, 0.0);
            for j in 0..n {
                let fd = (sides[1][j] - sides[0][j]) / (2.0 * epsilon);
                error += (derivative[j] * scale * duration - fd).powi(2);
                reference += fd * fd;
            }
            assert!(
                error.sqrt() <= 2e-7 * reference.sqrt(),
                "coordinate {coordinate}, relative error {}",
                (error / reference).sqrt()
            );
        }
    }

    #[test]
    fn full_chirp_jacobian_matches_refits_with_substantial_unmodeled_energy() {
        let (n, dt) = (401, 1.0 / 2048.0);
        let duration = (n - 1) as f64 * dt;
        let value = |p: ChirpedTone, a: C64, t: f64| {
            let rate = (p.end_frequency_hz - p.start_frequency_hz) / duration;
            let phase = core::f64::consts::TAU * (p.start_frequency_hz * t + 0.5 * rate * t * t);
            let (s, c) = phase.sin_cos();
            a.re * c + a.im * s
        };
        let truth = ChirpedTone {
            start_frequency_hz: 174.0,
            end_frequency_hz: 183.5,
        };
        let omitted = DampedTonePole {
            frequency_hz: 177.5,
            decay_rate_per_s: 0.0,
        };
        let samples: Vec<_> = (0..n)
            .map(|j| {
                let t = j as f64 * dt;
                0.8 + value(truth, C64::new(1.0, 0.5), t) + signal(omitted, C64::new(-0.4, 0.2), t)
            })
            .collect();
        let seed = ChirpedTone {
            start_frequency_hz: 173.2,
            end_frequency_hz: 182.8,
        };
        let budget = HarmonicFitBudget {
            max_workspace_bytes: 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        };
        let plan = Plan::Chirped(ChirpedToneFitPlan::new(n, dt, &[seed], budget).unwrap());
        let mut scratch = vec![0.0; n];
        plan.design()
            .fit_into(&samples, &mut scratch, &mut [C64::default()])
            .unwrap();
        let mut residual = vec![0.0; n];
        residual[3..].copy_from_slice(&scratch[3..]);
        plan.design().factor.apply_q(&mut residual);
        let lo = samples.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let center = lo.midpoint(hi);
        let scale = samples
            .iter()
            .map(|x| (x - center).abs())
            .fold(0.0, f64::max);
        for coordinate in 0..2 {
            let mut derivative = vec![0.0; n];
            plan.derivative(coordinate, &scratch[..3], &residual, &mut derivative);
            let span: f64 = derivative[..3].iter().map(|x| x * x).sum();
            let total: f64 = derivative.iter().map(|x| x * x).sum();
            assert!(
                span > 1e-3 * total,
                "the full residual-dependent term must be exercised"
            );
            plan.design().factor.apply_q(&mut derivative);
            let epsilon = 1e-4;
            let mut sides = Vec::new();
            for sign in [-1.0, 1.0] {
                let mut p = seed;
                if coordinate == 0 {
                    p.start_frequency_hz += sign * epsilon;
                } else {
                    p.end_frequency_hz += sign * epsilon;
                }
                let fit = ChirpedToneFitPlan::new(n, dt, &[p], budget)
                    .unwrap()
                    .fit(&samples)
                    .unwrap();
                sides.push(
                    (0..n)
                        .map(|j| fit.dc + value(p, fit.phasors[0], j as f64 * dt))
                        .collect::<Vec<_>>(),
                );
            }
            let (mut error, mut reference) = (0.0, 0.0);
            for j in 0..n {
                let fd = (sides[1][j] - sides[0][j]) / (2.0 * epsilon);
                error += (derivative[j] * scale * duration - fd).powi(2);
                reference += fd * fd;
            }
            assert!(
                error.sqrt() < 3e-7 * reference.sqrt(),
                "endpoint {coordinate}, relative Jacobian error {}",
                (error / reference).sqrt()
            );
        }
    }
}
