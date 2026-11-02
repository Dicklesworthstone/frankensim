//! Frequency knots for prescribed physical loads, with explicit reconstruction
//! error. These controls never change a resonator or material frequency.
use super::{Reader, tone_envelope::ToneEnvelopePlan, tone_record::PhysicalTone};
use fs_modal::time::harmonic::ToneModulationStep;

const MAX_FREQUENCY_POINTS: usize = 65_536;

struct FrequencyKnot {
    sample: u64,
    omega_rad_s: f64,
}

/// Immutable source history plus the cursor of the last accepted physical step.
/// The driver owns the unscaled quadrature phase; this owner supplies its exact
/// frequency integral and a bounded affine reconstruction for each interval.
pub(super) struct FrequencySweep {
    points: Vec<FrequencyKnot>,
    interval: usize,
    accepted_sample: u64,
    samples: u64,
    duration_s: f64,
    carrier_omega_rad_s: f64,
    force_amplitude_n: f64,
    maximum_force_error_n: f64,
    force_error_bound_n: f64,
}

/// Original tone-row indexing, before grouping voices or coalescing carriers.
pub(super) struct ToneFrequencyPlan {
    curves: Vec<Option<FrequencySweep>>,
    count: usize,
    maximum_force_error_bound_n: f64,
}

impl ToneFrequencyPlan {
    /// Optional `frequency_sweeps COUNT`, then `sweep TONE POINTS ERROR_N`
    /// and `omega SAMPLE RAD_S`. The error allowance applies separately to
    /// each original physical force row. It bounds interpolation of the complex
    /// amplitude around the fixed carrier, not the response or acoustic fit.
    ///
    /// Entire-curve extrema admit a conservative bound without scanning every
    /// sample. Bounds from disjoint gain/frequency segments may make admission
    /// conservative. This is an analytic interpolation bound evaluated in f64,
    /// not a validated interval bound including accumulated phase roundoff.
    pub(super) fn read_optional(
        reader: &mut Reader<'_>,
        records: &[PhysicalTone],
        samples: u64,
        sample_rate_hz: u32,
        nyquist_guard_fraction: f64,
        envelopes: Option<&ToneEnvelopePlan>,
        remaining_storage: &mut usize,
    ) -> Result<Option<Self>, String> {
        if reader
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next())
            != Some("frequency_sweeps")
        {
            return Ok(None);
        }
        let count: usize = reader.one("frequency_sweeps").map_err(|e| e.to_string())?;
        if count == 0 {
            return Ok(None);
        }
        if count > records.len()
            || samples == 0
            || sample_rate_hz == 0
            || !nyquist_guard_fraction.is_finite()
            || !(0.0 < nyquist_guard_fraction && nyquist_guard_fraction <= 1.0)
        {
            return Err(
                "frequency sweeps require existing tones and the original positive guarded clock"
                    .into(),
            );
        }
        let duration_s = 1.0 / f64::from(sample_rate_hz);
        let maximum_omega =
            core::f64::consts::PI * nyquist_guard_fraction * f64::from(sample_rate_hz);
        let metadata = records
            .len()
            .checked_mul(core::mem::size_of::<Option<FrequencySweep>>())
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<Self>()))
            .ok_or("frequency-sweep metadata size overflow")?;
        reserve(remaining_storage, metadata)?;
        let mut curves: Vec<Option<FrequencySweep>> = (0..records.len()).map(|_| None).collect();
        let mut total_points = 0usize;
        let mut maximum_force_error_bound_n = 0.0_f64;
        for _ in 0..count {
            let mut row = reader.row("sweep").map_err(|e| e.to_string())?;
            let tone: usize = row.parse().map_err(|e| e.to_string())?;
            let point_count = row
                .count(MAX_FREQUENCY_POINTS - total_points)
                .map_err(|e| e.to_string())?;
            let maximum_force_error_n = row.scalar().map_err(|e| e.to_string())?;
            row.finish().map_err(|e| e.to_string())?;
            if tone >= records.len()
                || curves[tone].is_some()
                || point_count == 0
                || maximum_force_error_n <= 0.0
            {
                return Err("each frequency sweep requires a unique tone, knots and positive force-error allowance".into());
            }
            let source = &records[tone];
            if !source.omega.is_finite()
                || source.omega < 0.0
                || source.omega >= maximum_omega
                || !source.cosine.is_finite()
                || !source.sine.is_finite()
            {
                return Err("frequency-sweep carrier must retain its original physical source and guarded clock".into());
            }
            reserve(
                remaining_storage,
                point_count
                    .checked_mul(core::mem::size_of::<FrequencyKnot>())
                    .ok_or("frequency-sweep knot size overflow")?,
            )?;
            total_points += point_count;
            let mut points: Vec<FrequencyKnot> = Vec::with_capacity(point_count);
            let mut maximum_detuning = 0.0_f64;
            let mut maximum_frequency_change = 0.0_f64;
            for point in 0..point_count {
                let mut row = reader.row("omega").map_err(|e| e.to_string())?;
                let sample: u64 = row.parse().map_err(|e| e.to_string())?;
                let omega_rad_s = row.scalar().map_err(|e| e.to_string())?;
                row.finish().map_err(|e| e.to_string())?;
                if sample > samples
                    || (point == 0 && sample != 0)
                    || points.last().is_some_and(|last| sample <= last.sample)
                    || omega_rad_s < 0.0
                    || omega_rad_s >= maximum_omega
                {
                    return Err("frequency knots must start at zero, increase within the horizon and remain below the original guarded Nyquist limit".into());
                }
                maximum_detuning =
                    maximum_detuning.max(((omega_rad_s - source.omega) * duration_s).abs());
                if let Some(last) = points.last() {
                    maximum_frequency_change = maximum_frequency_change.max(
                        (omega_rad_s - last.omega_rad_s).abs() / (sample - last.sample) as f64,
                    );
                }
                points.push(FrequencyKnot {
                    sample,
                    omega_rad_s,
                });
            }
            let [gain, gain_change] = envelopes.map_or(Ok([1.0, 0.0]), |plan| plan.bounds(tone))?;
            let force_amplitude_n = quadrature_norm(source.cosine, source.sine)?;
            // h² |z''| / 8, z=g exp(i integral(omega-carrier)):
            // D bounds |omega-carrier|h; gain_change and frequency_change
            // are changes per one sample. This also covers held final knots.
            let coefficient_error = 0.25 * gain_change * maximum_detuning
                + 0.125
                    * gain
                    * (maximum_frequency_change * duration_s + maximum_detuning * maximum_detuning);
            let force_error_bound_n = force_amplitude_n * coefficient_error;
            if !force_error_bound_n.is_finite() || force_error_bound_n > maximum_force_error_n {
                return Err("frequency sweep exceeds its physical force-reconstruction error allowance; refine the clock or choose a closer fixed carrier".into());
            }
            maximum_force_error_bound_n = maximum_force_error_bound_n.max(force_error_bound_n);
            curves[tone] = Some(FrequencySweep {
                points,
                interval: 0,
                accepted_sample: 0,
                samples,
                duration_s,
                carrier_omega_rad_s: source.omega,
                force_amplitude_n,
                maximum_force_error_n,
                force_error_bound_n,
            });
        }
        Ok(Some(Self {
            curves,
            count,
            maximum_force_error_bound_n,
        }))
    }

    pub(super) const fn count(&self) -> usize {
        self.count
    }

    pub(super) fn contains(&self, tone: usize) -> bool {
        self.curves.get(tone).is_some_and(Option::is_some)
    }

    /// Largest admitted per-original-tone force bound [N], not a sum of
    /// unrelated physical ports and not an acoustic pressure bound.
    pub(super) const fn maximum_force_error_bound_n(&self) -> f64 {
        self.maximum_force_error_bound_n
    }

    pub(super) fn take(&mut self, tone: usize) -> Option<FrequencySweep> {
        self.curves.get_mut(tone).and_then(Option::take)
    }

    pub(super) fn step(
        &self,
        tone: usize,
        sample: u64,
        gain: [f64; 2],
    ) -> Result<Option<ToneModulationStep>, String> {
        self.curves
            .get(tone)
            .ok_or("frequency sweep tone index is outside its original source")?
            .as_ref()
            .map(|curve| curve.step(sample, gain))
            .transpose()
    }

    pub(super) fn commit(&mut self, next_sample: u64) {
        for curve in self.curves.iter_mut().flatten() {
            curve.commit(next_sample);
        }
    }
}

impl FrequencySweep {
    pub(super) const fn force_error_bound_n(&self) -> f64 {
        self.force_error_bound_n
    }

    pub(super) fn step(&self, sample: u64, gain: [f64; 2]) -> Result<ToneModulationStep, String> {
        if sample != self.accepted_sample || sample >= self.samples {
            return Err("frequency sweep query does not match its accepted physical clock".into());
        }
        let next_sample = sample
            .checked_add(1)
            .ok_or("frequency sample clock overflow")?;
        let left = &self.points[self.interval];
        let omega = if let Some(right) = self.points.get(self.interval + 1) {
            if sample < left.sample || next_sample > right.sample {
                return Err("frequency sweep cursor does not cover this physical interval".into());
            }
            let value = |at| {
                if at == left.sample {
                    return left.omega_rad_s;
                }
                if at == right.sample {
                    return right.omega_rad_s;
                }
                let t = (at - left.sample) as f64 / (right.sample - left.sample) as f64;
                left.omega_rad_s * (1.0 - t) + right.omega_rad_s * t
            };
            [value(sample), value(next_sample)]
        } else {
            [left.omega_rad_s; 2]
        };
        let step =
            ToneModulationStep::try_new(self.duration_s, self.carrier_omega_rad_s, omega, gain)
                .map_err(|e| e.to_string())?;
        let error = self.force_amplitude_n * step.coefficient_error_bound();
        if !error.is_finite() || error > self.maximum_force_error_n {
            return Err(
                "frequency interval exceeds its physical force-reconstruction error allowance"
                    .into(),
            );
        }
        Ok(step)
    }

    /// Publish only after the mechanical and receiver step has accepted.
    pub(super) fn commit(&mut self, next_sample: u64) {
        self.accepted_sample = next_sample;
        while self
            .points
            .get(self.interval + 1)
            .is_some_and(|point| point.sample <= next_sample)
        {
            self.interval += 1;
        }
    }
}

fn quadrature_norm(cosine: f64, sine: f64) -> Result<f64, String> {
    let scale = cosine.abs().max(sine.abs());
    let norm = if scale == 0.0 {
        0.0
    } else {
        let c = cosine / scale;
        let s = sine / scale;
        scale * fs_math::det::sqrt(c * c + s * s)
    };
    if !norm.is_finite() {
        return Err("physical tone magnitude is not representable".into());
    }
    Ok(norm)
}

fn reserve(remaining: &mut usize, bytes: usize) -> Result<(), String> {
    *remaining = remaining
        .checked_sub(bytes)
        .ok_or("physical frequency sweeps exceed source storage")?;
    Ok(())
}

#[cfg(test)]
#[path = "tone_frequency/tests.rs"]
mod tests;
