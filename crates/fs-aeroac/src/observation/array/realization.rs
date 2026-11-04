//! Coherent physical waveform realizations of the original mean/factor fields.
//!
//! One supplied coefficient multiplies its named factor at every frequency and
//! receiver. It remains fixed throughout the realized history. This represents
//! one member of the authored harmonic ensemble, not fresh independent noise
//! per sample, receiver or source. No distribution or random phase is invented.

use super::*;

/// One caller-supplied coefficient for a named physical fluctuation factor.
///
/// The ensemble law concerns a collection of draws: real coefficients have
/// E[xi]=0, E[xi²]=1; proper coefficients have E[z]=E[z²]=0, E[|z|²]=1.
/// Individual draws need not have unit magnitude and are never normalized.
/// Realized peak fields multiply these values directly, without a sqrt(2),
/// half-power, conjugation or additional phase convention conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HarmonicFactorDraw {
    /// Finite signed real coefficient for a Real coordinate.
    Real(f64),
    /// Finite complex coefficient for a ProperComplex coordinate.
    ProperComplex(C64),
}

/// Preparation work over all physical receivers or one processed channel.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicRealizationBudget {
    /// Maximum receiver fields, or at least one for a processed readout.
    pub max_channels: usize,
    /// Maximum shared coefficient draws, including silent factor coordinates.
    pub max_factors: usize,
    /// Sum of supplied mean/factor lines over the admitted channels.
    pub max_total_lines: usize,
    /// Sum over channels of (mean lines + all factor lines)². Realization
    /// combines factors coherently, so this includes cross-factor grouping.
    pub max_line_pairs: usize,
}

/// One original receiver and its coherent deterministic realization.
#[derive(Clone, Debug)]
pub struct RealizedHarmonicReceiver<const D: usize> {
    /// Original physical receiver position [m], with no new propagation.
    pub position_m: [f64; D],
    /// Original material coefficients used by physical microphone processing.
    pub acoustic: FluidAcousticState,
    /// Realized peak pressure and velocity, retaining all physical frequencies.
    pub field: PreparedHarmonicField<D>,
}

/// One fixed physical draw shared by every receiver and processed readout.
///
/// All source and propagation phases come from the original ensembles. Exact
/// frequency sums retain DD coefficients across the mean and every realized
/// factor, including weak cancellation residuals. Zero draws and coherent
/// nulls retain their physical lines for strict subsequent Nyquist admission.
/// Preparation is bounded and atomic; sampling allocates nothing. No ensemble
/// statistics can be inferred from one deterministic realization.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicArrayRealization<'a, const D: usize> {
    source_clock: &'a str,
    coordinates: &'a [HarmonicFactorCoordinate<'a>],
    draws: Vec<HarmonicFactorDraw>,
    receivers: Vec<RealizedHarmonicReceiver<D>>,
}

impl<'a, const D: usize> PreparedHarmonicArray<'a, D> {
    /// Apply the SAME supplied factor draw to every receiver and frequency.
    /// A draw's variant must match the corresponding declared factor law.
    pub fn realize(
        &self,
        draws: &[HarmonicFactorDraw],
        budget: HarmonicRealizationBudget,
    ) -> Result<PreparedHarmonicArrayRealization<'a, D>, AeroacError> {
        self.realize_cancellable(draws, budget, || false)
    }

    /// Bounded realization with no partial result on cancellation or refusal.
    pub fn realize_cancellable(
        &self,
        draws: &[HarmonicFactorDraw],
        budget: HarmonicRealizationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicArrayRealization<'a, D>, AeroacError> {
        poll(&mut cancelled)?;
        if draws.len() != self.coordinates.len()
            || draws.len() > budget.max_factors
            || self.receivers.len() > budget.max_channels
        {
            return Err(invalid(
                "harmonic realization factor/channel count or budget",
            ));
        }
        for (&draw, coordinate) in draws.iter().zip(self.coordinates) {
            poll(&mut cancelled)?;
            if !draw.matches(coordinate.statistics) {
                return Err(invalid(
                    "harmonic realization draw law or finite coefficient",
                ));
            }
        }
        let mut total_lines = 0usize;
        let mut line_pairs = 0usize;
        for receiver in &self.receivers {
            poll(&mut cancelled)?;
            admit_lines(receiver.ensemble, &mut total_lines, &mut line_pairs, budget)?;
        }
        let mut retained_draws = Vec::new();
        retained_draws
            .try_reserve_exact(draws.len())
            .map_err(|_| invalid("harmonic realization draw allocation failed"))?;
        retained_draws.extend_from_slice(draws);
        let mut receivers = Vec::new();
        receivers
            .try_reserve_exact(self.receivers.len())
            .map_err(|_| invalid("harmonic realization receiver allocation failed"))?;
        for receiver in &self.receivers {
            poll(&mut cancelled)?;
            receivers.push(RealizedHarmonicReceiver {
                position_m: receiver.position_m,
                acoustic: receiver.acoustic,
                field: realize_field(receiver.ensemble, draws, &mut cancelled)?,
            });
        }
        poll(&mut cancelled)?;
        Ok(PreparedHarmonicArrayRealization {
            source_clock: self.source_clock,
            coordinates: self.coordinates,
            draws: retained_draws,
            receivers,
        })
    }
}

impl<'a, const D: usize> PreparedHarmonicArrayRealization<'a, D> {
    /// Original common source-clock identity.
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }

    /// Original ordered physical coefficient identities and ensemble laws.
    pub fn factor_coordinates(&self) -> &[HarmonicFactorCoordinate<'a>] {
        self.coordinates
    }

    /// Exact supplied coefficient draws; no normalization or RNG is hidden.
    pub fn factor_draws(&self) -> &[HarmonicFactorDraw] {
        &self.draws
    }

    /// Original receiver order, geometry and material with realized fields.
    pub fn receivers(&self) -> &[RealizedHarmonicReceiver<D>] {
        &self.receivers
    }

    /// Realize an original microphone/beamforming readout with this SAME draw.
    /// The readout must carry matching clock and ordered physical coordinates.
    /// Its original p/u projection, signed gains and causal delays are already
    /// present, so they are neither re-derived nor applied a second time.
    pub fn realize_readout(
        &self,
        readout: &PreparedHarmonicReadout<'_>,
        budget: HarmonicRealizationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicField<0>, AeroacError> {
        poll(&mut cancelled)?;
        if self.source_clock != readout.source_clock()
            || self.coordinates != readout.factor_coordinates()
            || budget.max_channels == 0
            || self.draws.len() > budget.max_factors
        {
            return Err(invalid(
                "harmonic readout realization clock, coordinates or budget",
            ));
        }
        admit_lines(readout.ensemble(), &mut 0, &mut 0, budget)?;
        let field = realize_field(readout.ensemble(), &self.draws, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(field)
    }

    /// Sample all physical receivers on one absolute clock in sample-major
    /// order: output[frame * receiver_count + receiver]. Pressure remains Pa
    /// and each actual velocity component remains m/s. `max_samples` bounds
    /// receiver samples (output length), while `max_tone_samples` bounds
    /// frame count times the total distinct line count over all receivers.
    ///
    /// Shape, clock, work and every receiver's original-line Nyquist admission
    /// are checked before writing output. Runtime non-finite results or
    /// cancellation may leave a prefix; discard the buffer on any error.
    /// Blocks use absolute sample indices and cannot reset phases or redraw
    /// factors. This is settled harmonic continuation, not a transient onset.
    pub fn sample_fields_into(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [HarmonicFieldSample<D>],
        budget: HarmonicSampleBudget,
    ) -> Result<(), AeroacError> {
        self.sample_fields_into_cancellable(clock, first_sample, output, budget, || false)
    }

    /// Allocation-free joint output with bounded cancellation polls.
    pub fn sample_fields_into_cancellable(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [HarmonicFieldSample<D>],
        budget: HarmonicSampleBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(&mut cancelled)?;
        let receiver_count = self.receivers.len();
        if output.len() % receiver_count != 0 || output.len() > budget.max_samples {
            return Err(invalid("harmonic array sample shape or output budget"));
        }
        let frames = output.len() / receiver_count;
        let mut lines = 0usize;
        for receiver in &self.receivers {
            poll(&mut cancelled)?;
            lines = lines
                .checked_add(receiver.field.lines.len())
                .ok_or_else(|| invalid("harmonic array sample line count overflow"))?;
        }
        if frames
            .checked_mul(lines)
            .is_none_or(|work| work > budget.max_tone_samples)
        {
            return Err(invalid("harmonic array sampling exceeds tone work budget"));
        }
        let count =
            u64::try_from(frames).map_err(|_| invalid("harmonic array frame count exceeds u64"))?;
        for receiver in &self.receivers {
            poll(&mut cancelled)?;
            receiver.field.validate_sample_range_cancellable(
                clock,
                first_sample,
                count,
                &mut cancelled,
            )?;
        }
        for (offset, frame) in output.chunks_exact_mut(receiver_count).enumerate() {
            let time = clock.time_at(first_sample + offset as u64)?;
            for (sample, receiver) in frame.iter_mut().zip(&self.receivers) {
                *sample = receiver
                    .field
                    .field_at_clock_dd_cancellable(time, &mut cancelled)?;
            }
        }
        poll(&mut cancelled)
    }
}

impl HarmonicFactorDraw {
    fn matches(self, law: HarmonicFactorStatistics) -> bool {
        match (self, law) {
            (Self::Real(value), HarmonicFactorStatistics::Real) => value.is_finite(),
            (Self::ProperComplex(value), HarmonicFactorStatistics::ProperComplex) => finite(value),
            _ => false,
        }
    }

    fn parts(self) -> [Dd; 2] {
        match self {
            Self::Real(value) => [Dd::from_f64(value), Dd::ZERO],
            Self::ProperComplex(value) => [Dd::from_f64(value.re), Dd::from_f64(value.im)],
        }
    }
}

fn line_count<const D: usize>(
    ensemble: &PreparedHarmonicEnsemble<D>,
) -> Result<usize, AeroacError> {
    ensemble
        .factors
        .iter()
        .try_fold(ensemble.mean.lines.len(), |sum, (_, field)| {
            sum.checked_add(field.lines.len())
                .ok_or_else(|| invalid("harmonic realization line count overflow"))
        })
}

fn admit_lines<const D: usize>(
    ensemble: &PreparedHarmonicEnsemble<D>,
    total: &mut usize,
    pairs: &mut usize,
    budget: HarmonicRealizationBudget,
) -> Result<(), AeroacError> {
    let count = line_count(ensemble)?;
    *total = total
        .checked_add(count)
        .ok_or_else(|| invalid("harmonic realization total lines overflow"))?;
    *pairs = count
        .checked_mul(count)
        .and_then(|n| pairs.checked_add(n))
        .ok_or_else(|| invalid("harmonic realization pair work overflow"))?;
    if *total > budget.max_total_lines || *pairs > budget.max_line_pairs {
        return Err(invalid(
            "harmonic realization exceeds line/pair work budget",
        ));
    }
    Ok(())
}

#[allow(clippy::float_cmp)] // Only exactly equal physical frequencies coalesce.
fn realize_field<const D: usize>(
    ensemble: &PreparedHarmonicEnsemble<D>,
    draws: &[HarmonicFactorDraw],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedHarmonicField<D>, AeroacError> {
    let mut lines: Vec<PreparedLine<D>> = Vec::new();
    lines
        .try_reserve_exact(line_count(ensemble)?)
        .map_err(|_| invalid("harmonic realized field allocation failed"))?;
    for bank in 0..=draws.len() {
        poll(cancelled)?;
        let [re, im] = if bank == 0 {
            [Dd::ONE, Dd::ZERO]
        } else {
            draws[bank - 1].parts()
        };
        for line in &field(ensemble, bank).lines {
            poll(cancelled)?;
            let scale = |value: ComplexSum| ComplexSum {
                re: value.re * re - value.im * im,
                im: value.re * im + value.im * re,
            };
            let scaled = PreparedLine {
                frequency_hz: line.frequency_hz,
                pressure: scale(line.pressure),
                velocity: line.velocity.map(scale),
            };
            if !scaled.pressure.is_finite() || scaled.velocity.iter().any(|v| !v.is_finite()) {
                return Err(nonfinite("harmonic realization factor multiplication"));
            }
            let mut group = None;
            for (index, previous) in lines.iter().enumerate() {
                poll(cancelled)?;
                if previous.frequency_hz == line.frequency_hz {
                    group = Some(index);
                    break;
                }
            }
            if let Some(index) = group {
                let append = |sum: &mut ComplexSum, value: ComplexSum| {
                    sum.re = sum.re + value.re;
                    sum.im = sum.im + value.im;
                };
                append(&mut lines[index].pressure, scaled.pressure);
                for (sum, value) in lines[index].velocity.iter_mut().zip(scaled.velocity) {
                    append(sum, value);
                }
            } else {
                lines.push(scaled);
            }
        }
    }
    if lines
        .iter()
        .any(|line| !line.pressure.is_finite() || line.velocity.iter().any(|v| !v.is_finite()))
    {
        return Err(nonfinite("coalesced harmonic realization"));
    }
    Ok(PreparedHarmonicField { lines })
}
