//! Physical first-order microphone projection and authored delay-and-sum signals.

use super::*;
use fs_fft::harmonic_window::phase_hz_from_clock_dd;

/// One signed, causal tap on an original physical receiver.
/// The ideal first-order microphone reads alpha*p - (1-alpha)*rho*c*n.u,
/// with n pointing toward its front. Actual near-field velocity is used;
/// no pressure-to-velocity or intensity approximation is introduced.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReadoutTap<const D: usize> {
    /// Physical receiver index in the original array.
    pub receiver: usize,
    /// Finite signed dimensionless channel weight.
    pub gain: f64,
    /// Finite nonnegative processing delay [s].
    pub delay_s: f64,
    /// In [0,1]: 1 is pressure, 0 is figure eight, 1/2 is cardioid.
    pub pressure_fraction: f64,
    /// Finite nonzero front direction, normalized during preparation.
    /// It may be zero for a pure pressure readout.
    pub front_axis: [f64; D],
}

impl<const D: usize> HarmonicReadoutTap<D> {
    /// A pressure-only tap without directional velocity response.
    pub const fn pressure(receiver: usize, gain: f64, delay_s: f64) -> Self {
        Self {
            receiver,
            gain,
            delay_s,
            pressure_fraction: 1.0,
            front_axis: [0.0; D],
        }
    }
}

/// Admission limits on readout preparation and subsequent observations.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicReadoutBudget {
    /// Maximum signed receiver taps, including repeated receivers.
    pub max_taps: usize,
    /// Sum of tap/line contributions over the mean and every shared factor.
    pub max_total_lines: usize,
    /// Sum of squared contribution counts per bank, bounding exact-frequency
    /// coalescing and subsequent finite-window observations.
    pub max_line_pairs: usize,
}

/// Processed pressure-equivalent signal; it is not a new spatial receiver,
/// particle-velocity field, energy-flux estimate or acoustic loading port.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicReadout<'a> {
    source_clock: &'a str,
    coordinates: &'a [HarmonicFactorCoordinate<'a>],
    ensemble: PreparedHarmonicEnsemble<0>,
}

/// Joint statistics of two processed pressure-equivalent channels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowedHarmonicReadoutPair {
    /// Common physical observer interval [s].
    pub window: HarmonicWindow,
    /// <E[y_left(t) y_right(t)]> [Pa²], including deterministic waveforms.
    pub mean_product_pa2: f64,
    /// Cov(<y_left>, <y_right>) [Pa²], from fluctuations only.
    pub averaged_covariance_pa2: f64,
}

struct PreparedTap<const D: usize> {
    receiver: usize,
    delay_s: Dd,
    gain: Dd,
    pressure: Dd,
    velocity: [Dd; D],
}

impl<'a, const D: usize> PreparedHarmonicArray<'a, D> {
    /// Prepare y(t)=sum_tap gain*mic(t-delay), retaining one common mean and
    /// the SAME real/proper coefficient coordinates across every tap.
    /// Pressure and velocity interfere before any moment is formed. Delays
    /// add exp(+i*2*pi*f*delay) once under exp(-i omega t); original propagation
    /// is already in the fields. This uses their authored harmonic continuation,
    /// not a finite pulse onset or sampled fractional-delay filter.
    ///
    /// Exact-frequency sums, projection and delay phase retain DD coefficients
    /// through cancellation. Preparation is bounded; subsequent observations
    /// delegate to the original harmonic ensemble and allocate nothing.
    pub fn prepare_readout(
        &self,
        taps: &[HarmonicReadoutTap<D>],
        budget: HarmonicReadoutBudget,
    ) -> Result<PreparedHarmonicReadout<'a>, AeroacError> {
        self.prepare_readout_cancellable(taps, budget, || false)
    }

    /// Prepare a readout atomically, polling through taps and line coalescing.
    #[allow(clippy::float_cmp)] // Exactly alpha=1 removes the velocity response.
    pub fn prepare_readout_cancellable(
        &self,
        taps: &[HarmonicReadoutTap<D>],
        budget: HarmonicReadoutBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicReadout<'a>, AeroacError> {
        poll(&mut cancelled)?;
        if taps.len() > budget.max_taps {
            return Err(invalid("acoustic readout exceeds tap budget"));
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(taps.len())
            .map_err(|_| invalid("acoustic readout tap allocation failed"))?;
        for tap in taps {
            poll(&mut cancelled)?;
            let receiver = self.receiver(tap.receiver)?;
            if !tap.gain.is_finite()
                || !tap.delay_s.is_finite()
                || tap.delay_s < 0.0
                || !tap.pressure_fraction.is_finite()
                || !(0.0..=1.0).contains(&tap.pressure_fraction)
                || tap.front_axis.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "acoustic readout gain, delay or microphone response",
                ));
            }
            let mut velocity = [Dd::ZERO; D];
            if tap.pressure_fraction != 1.0 {
                let scale = tap.front_axis.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
                if scale == 0.0 {
                    return Err(invalid("directional microphone needs a front axis"));
                }
                let normal = tap.front_axis.map(|v| v / scale);
                let length = normal.iter().fold(0.0, |r, &v| det::hypot(r, v));
                let impedance = Dd::from_f64(receiver.acoustic.characteristic_impedance_pa_s_m());
                velocity = normal.map(|v| {
                    -Dd::from_f64(v / length)
                        * (Dd::ONE - Dd::from_f64(tap.pressure_fraction))
                        * impedance
                });
            }
            prepared.push(PreparedTap {
                receiver: tap.receiver,
                delay_s: Dd::from_f64(tap.delay_s),
                gain: Dd::from_f64(tap.gain),
                pressure: Dd::from_f64(tap.pressure_fraction),
                velocity,
            });
        }
        let mut total = 0usize;
        let mut pairs = 0usize;
        for bank in 0..=self.coordinates.len() {
            let mut lines = 0usize;
            for tap in &prepared {
                poll(&mut cancelled)?;
                lines = lines
                    .checked_add(
                        field(self.receivers[tap.receiver].ensemble, bank)
                            .lines
                            .len(),
                    )
                    .ok_or_else(|| invalid("acoustic readout line count overflows"))?;
            }
            total = total
                .checked_add(lines)
                .ok_or_else(|| invalid("acoustic readout line count overflows"))?;
            pairs = lines
                .checked_mul(lines)
                .and_then(|n| pairs.checked_add(n))
                .ok_or_else(|| invalid("acoustic readout pair work overflows"))?;
        }
        if total > budget.max_total_lines || pairs > budget.max_line_pairs {
            return Err(invalid("acoustic readout exceeds line/pair budget"));
        }
        let mean = self.readout_bank(0, &prepared, &mut cancelled)?;
        let mut factors = Vec::new();
        factors
            .try_reserve_exact(self.coordinates.len())
            .map_err(|_| invalid("acoustic readout factor allocation failed"))?;
        for (index, coordinate) in self.coordinates.iter().enumerate() {
            factors.push((
                coordinate.statistics,
                self.readout_bank(index + 1, &prepared, &mut cancelled)?,
            ));
        }
        poll(&mut cancelled)?;
        Ok(PreparedHarmonicReadout {
            source_clock: self.source_clock,
            coordinates: self.coordinates,
            ensemble: PreparedHarmonicEnsemble { mean, factors },
        })
    }

    #[allow(clippy::float_cmp)] // Only exactly equal frequencies may coalesce.
    fn readout_bank(
        &self,
        bank: usize,
        taps: &[PreparedTap<D>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicField<0>, AeroacError> {
        let count = taps
            .iter()
            .map(|tap| {
                field(self.receivers[tap.receiver].ensemble, bank)
                    .lines
                    .len()
            })
            .sum();
        let mut lines: Vec<PreparedLine<0>> = Vec::new();
        lines
            .try_reserve_exact(count)
            .map_err(|_| invalid("acoustic readout line allocation failed"))?;
        for tap in taps {
            for line in &field(self.receivers[tap.receiver].ensemble, bank).lines {
                poll(cancelled)?;
                let mut value = ComplexSum {
                    re: line.pressure.re * tap.pressure,
                    im: line.pressure.im * tap.pressure,
                };
                for (u, coefficient) in line.velocity.iter().zip(tap.velocity) {
                    value.re = value.re + u.re * coefficient;
                    value.im = value.im + u.im * coefficient;
                }
                let [re, im] = phase_hz_from_clock_dd(line.frequency_hz, -tap.delay_s)
                    .map_err(window_error)?;
                let rotated = ComplexSum {
                    re: (value.re * re - value.im * im) * tap.gain,
                    im: (value.re * im + value.im * re) * tap.gain,
                };
                if !rotated.is_finite() {
                    return Err(nonfinite("acoustic readout projection or delay"));
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
                    lines[index].pressure.re = lines[index].pressure.re + rotated.re;
                    lines[index].pressure.im = lines[index].pressure.im + rotated.im;
                } else {
                    lines.push(PreparedLine {
                        frequency_hz: line.frequency_hz,
                        pressure: rotated,
                        velocity: [],
                    });
                }
            }
        }
        if lines.iter().any(|line| !line.pressure.is_finite()) {
            return Err(nonfinite("coalesced acoustic readout"));
        }
        Ok(PreparedHarmonicField { lines })
    }
}

impl PreparedHarmonicReadout<'_> {
    /// Original authored source-clock identity, retained after processing delay.
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }

    /// Ordered physical factor identities retained through microphone processing.
    pub fn factor_coordinates(&self) -> &[HarmonicFactorCoordinate<'_>] {
        self.coordinates
    }

    /// Original pressure-only ensemble. Its zero velocity dimension explicitly
    /// prevents a processed channel from claiming a spatial intensity vector.
    pub fn ensemble(&self) -> &PreparedHarmonicEnsemble<0> {
        &self.ensemble
    }

    /// Original ensemble moments of the fully interfered processed waveform.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<WindowedAcousticEnsemble<0>, AeroacError> {
        self.ensemble.observe_window(window)
    }

    /// Allocation-free observation with the original ensemble cancellation gate.
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedAcousticEnsemble<0>, AeroacError> {
        self.ensemble.observe_window_cancellable(window, cancelled)
    }

    /// Joint statistics of two authored processed channels. Shared factor
    /// identities/laws and clock must match exactly. `maximum_line_pairs`
    /// bounds sum_banks(left_lines * right_lines); a covariance query needs at
    /// most two additional single-line integrals per admitted pair.
    pub fn observe_pair(
        &self,
        right: &Self,
        window: HarmonicWindow,
        maximum_line_pairs: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WindowedHarmonicReadoutPair, AeroacError> {
        poll(&mut cancelled)?;
        if self.source_clock != right.source_clock || self.coordinates != right.coordinates {
            return Err(invalid(
                "acoustic readouts do not share clock and factor coordinates",
            ));
        }
        let mut pairs = 0usize;
        for bank in 0..=self.coordinates.len() {
            poll(&mut cancelled)?;
            pairs = field(&self.ensemble, bank)
                .lines
                .len()
                .checked_mul(field(&right.ensemble, bank).lines.len())
                .and_then(|n| pairs.checked_add(n))
                .ok_or_else(|| invalid("acoustic readout cross work overflows"))?;
        }
        if pairs > maximum_line_pairs {
            return Err(invalid("acoustic readout cross work exceeds budget"));
        }
        let (moments, covariance) = pair(
            &self.ensemble,
            &right.ensemble,
            window,
            core::ptr::eq(self, right),
            &mut cancelled,
        )?;
        Ok(WindowedHarmonicReadoutPair {
            window,
            mean_product_pa2: moments.pressure_pressure_pa2,
            averaged_covariance_pa2: covariance.pressure_pressure_pa2,
        })
    }
}
