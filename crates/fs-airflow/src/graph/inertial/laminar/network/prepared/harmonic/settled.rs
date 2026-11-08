//! Settled simultaneous tones in every original network history coordinate.

use super::*;
use crate::graph::inertial::laminar::tones::LaminarPressureTone;

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Initialize the fully settled response to simultaneous pressure tones.
    ///
    /// Sources use exactly the port-major layout of `prepare_tones`: reservoir
    /// pressures first, then one directed head per branch, with every frequency
    /// present at every port. At the returned clock, each source starts as
    /// `cosine*cos(omega*tau) + sine*sin(omega*tau)`, where `tau = 0` at
    /// `accepted_time_s`. The clock does not rotate the supplied phases.
    ///
    /// DC uses the cached original steady solution. Positive frequencies reuse
    /// one complete native harmonic model, with one factorization per frequency.
    /// Under its `exp(-i*omega*tau)` convention the input is `cosine + i*sine`;
    /// the real response at `tau = 0` includes every original viscous-memory
    /// coordinate. Frequencies may be incommensurate, so a settled response need
    /// not have a common period.
    ///
    /// Each positive frequency must lie in every original pipe's admitted band.
    /// Each branch must also admit `abs(Q_dc) + sum(abs(Q_tone))`, with a
    /// roundoff allowance, over the entire waveform. This phase-independent
    /// bound can conservatively refuse tones whose peaks cannot coincide.
    ///
    /// Frequencies must be nonempty, finite, nonnegative and strictly increasing.
    /// DC sine coefficients must be zero. The source extent is bounded by
    /// `MAX_LAMINAR_NETWORK_AUGMENTED`; the harmonic budget applies to each
    /// native preparation/factor/response and is unused for DC alone. All work
    /// is cold and staged locally. Cancellation or refusal publishes no state.
    pub fn settled_tone_state(
        &self,
        cx: &Cx<'_>,
        frequencies_rad_s: &[f64],
        pressures: &[LaminarPressureTone],
        accepted_time_s: f64,
        budget: HarmonicBudget,
    ) -> Result<LaminarNetworkState, HarmonicError> {
        poll(cx)?;
        let count = frequencies_rad_s.len();
        let ports = self.port_count();
        let extent = ports
            .checked_mul(count)
            .ok_or_else(|| invalid("settled network tone extent overflow"))?;
        if count == 0
            || extent > MAX_LAMINAR_NETWORK_AUGMENTED
            || pressures.len() != extent
            || !accepted_time_s.is_finite()
            || accepted_time_s < 0.0
            || frequencies_rad_s.iter().any(|w| !w.is_finite() || *w < 0.0)
            || frequencies_rad_s.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(invalid(
                "bounded complete settled tones, increasing nonnegative frequencies and finite nonnegative clock required",
            ));
        }
        let has_dc = frequencies_rad_s[0] == 0.0;
        for (index, pressure) in pressures.iter().enumerate() {
            if !pressure.cosine.value().is_finite()
                || !pressure.sine.value().is_finite()
                || (has_dc && index % count == 0 && pressure.sine.value() != 0.0)
            {
                return Err(invalid(
                    "finite settled pressure tones and zero DC sine coefficients required",
                ));
            }
        }
        for &omega in frequencies_rad_s {
            if omega > 0.0 {
                self.admit_harmonic_frequency(cx, omega)?;
            }
        }

        let boundary = self.topology.pressure_nodes.len();
        let dc: Vec<_> = (0..ports)
            .map(|port| {
                if has_dc {
                    pressures[port * count].cosine
                } else {
                    Pressure::new(0.0)
                }
            })
            .collect();
        let steady = self
            .steady_state_with_heads(&dc[..boundary], &dc[boundary..], accepted_time_s)
            .map_err(|_| invalid("settled DC source or original steady flow is inadmissible"))?;
        poll(cx)?;
        if has_dc && count == 1 {
            return Ok(steady);
        }

        let tolerance = harmonic_tolerance(self, ports, budget)?;
        let model = self.harmonic_model(cx, budget)?;
        let n = self.system.state_dim();
        let mut amplitudes = vec![Sum::default(); self.branches.len()];
        for (edge, amplitude) in amplitudes.iter_mut().enumerate() {
            amplitude.add(dot(&self.flow_rows[edge], &steady.coordinates).abs())?;
            let mut projection = Sum::default();
            for (&coefficient, &value) in self.flow_rows[edge].iter().zip(&steady.coordinates) {
                projection.add(coefficient.abs() * value.abs())?;
            }
            amplitude.add(ROUNDOFF * n as f64 * projection.total()?)?;
        }
        let mut coordinates: Vec<_> = steady
            .coordinates
            .into_iter()
            .map(|sum| Sum {
                sum,
                correction: 0.0,
            })
            .collect();
        let mut phasors = vec![C64::ZERO; ports];
        for (tone, &omega) in frequencies_rad_s.iter().enumerate() {
            if omega == 0.0 {
                continue;
            }
            poll(cx)?;
            for (port, phasor) in phasors.iter_mut().enumerate() {
                let source = pressures[port * count + tone];
                *phasor = C64::new(source.cosine.value(), source.sine.value());
            }
            let input =
                self.harmonic_input_phasors(&phasors[..boundary], Some(&phasors[boundary..]))?;
            let factor = model.factor(omega, || cx.checkpoint().is_err())?;
            let response = self.response_with_factor(cx, &factor, &input, tolerance)?;
            for (coordinate, value) in coordinates.iter_mut().zip(&response.state) {
                coordinate.add(value.re)?;
            }
            for (edge, (amplitude, flow)) in amplitudes
                .iter_mut()
                .zip(&response.output[boundary..])
                .enumerate()
            {
                amplitude.add(flow.abs())?;
                // A small net flow can hide large cancelling projections. Bound
                // their absolute contributions before applying summation slack.
                let mut projection = Sum::default();
                for (&coefficient, value) in self.flow_rows[edge].iter().zip(&response.state) {
                    projection.add(coefficient.abs() * value.abs())?;
                }
                amplitude.add(ROUNDOFF * n as f64 * projection.total()?)?;
            }
        }
        let coordinates = coordinates
            .into_iter()
            .map(Sum::total)
            .collect::<Result<Vec<_>, _>>()?;
        let state = self
            .state(coordinates, accepted_time_s)
            .map_err(|_| invalid("settled complete network history is not representable"))?;
        let allowance = 1.0 + ROUNDOFF * (n + count) as f64;
        for (edge, amplitude) in amplitudes.into_iter().enumerate() {
            poll(cx)?;
            let peak = amplitude.total()?;
            let bound = if peak == 0.0 {
                0.0
            } else {
                (peak * allowance).next_up()
            };
            if !bound.is_finite() {
                return Err(invalid(
                    "settled branch waveform bound is not representable",
                ));
            }
            let actual = dot(&self.flow_rows[edge], &state.coordinates).abs();
            if !actual.is_finite() {
                return Err(invalid("settled branch phase flow is not representable"));
            }
            self.branches[edge]
                .admit_flow(bound.max(actual))
                .map_err(|_| {
                    invalid("settled simultaneous branch peaks exceed original laminar admission")
                })?;
        }
        poll(cx)?;
        Ok(state)
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) -> Result<(), HarmonicError> {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        if !self.sum.is_finite() || !self.correction.is_finite() {
            return Err(invalid("settled tone accumulation is not representable"));
        }
        Ok(())
    }

    fn total(self) -> Result<f64, HarmonicError> {
        let value = self.sum + self.correction;
        if !value.is_finite() {
            return Err(invalid("settled tone accumulation is not representable"));
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
