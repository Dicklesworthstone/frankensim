//! Original-history sampling and complete, coherent causal field contraction.

use super::super::history::PreparedQuadrupoleHistory;
use super::super::{Projection, Transfer, TransientQuadrupoleSource};
use super::*;

fn finite_source(source: &TransientQuadrupoleSource) -> bool {
    source
        .stress_moment_n_m
        .iter()
        .chain(&source.stress_rate_n_m_s)
        .chain(&source.stress_acceleration_n_m_s2)
        .chain(&source.stress_integral_n_m_s)
        .all(|v| v.is_finite())
}

fn compatible_initial(source: &TransientQuadrupoleSource) -> bool {
    finite_source(source)
        && source
            .stress_moment_n_m
            .iter()
            .chain(&source.stress_rate_n_m_s)
            .chain(&source.stress_integral_n_m_s)
            .all(|&v| v == 0.0)
}

// Row m is acted on by the same F_m for ALL four physical field components.
// Cached scales already contain geometric spreading, rho, c_inf and image sign.
pub(super) fn inputs(
    path: &Transfer,
    source: TransientQuadrupoleSource,
) -> Result<[[f64; 4]; 4], AeroacError> {
    if !finite_source(&source) {
        return Err(nonfinite("molecular quadrupole source history"));
    }
    let n = path.direction.map(Dd::from_f64);
    let integral = Projection::new(source.stress_integral_n_m_s, n);
    let moment = Projection::new(source.stress_moment_n_m, n);
    let rate = Projection::new(source.stress_rate_n_m_s, n);
    let acceleration = Projection::new(source.stress_acceleration_n_m_s2, n);
    let p = path.pressure_scales.map(Dd::from_f64);
    let v = path.velocity_scales.map(Dd::from_f64);
    let mut rows = [[Dd::ZERO; 4]; 4];
    rows[0][0] = moment.pressure_near() * p[0];
    rows[1][0] = rate.pressure_near() * p[1];
    rows[2][0] = acceleration.radial * p[2];
    for axis in 0..3 {
        rows[0][axis + 1] = integral.velocity_near(n[axis], axis) * v[0];
        rows[1][axis + 1] = moment.velocity_near(n[axis], axis) * v[1];
        rows[2][axis + 1] = ((rate.radial * Dd::from_f64(6.0) - rate.trace) * n[axis]
            - rate.vector[axis] * Dd::from_f64(2.0))
            * v[2];
        rows[3][axis + 1] = acceleration.radial * n[axis] * v[3];
    }
    if let Some(inverse) = path.inverse_image {
        for row in &mut rows {
            let velocity = [row[1], row[2], row[3]];
            for axis in 0..3 {
                row[axis + 1] = (0..3).fold(Dd::ZERO, |sum, virtual_axis| {
                    sum + Dd::from_f64(inverse[virtual_axis][axis]) * velocity[virtual_axis]
                });
            }
        }
    }
    let rows = rows.map(|row| row.map(Dd::to_f64));
    if rows.iter().flatten().any(|v| !v.is_finite()) {
        return Err(nonfinite("molecular quadrupole projected source"));
    }
    Ok(rows)
}

impl PreparedMolecularQuadrupoleObserver {
    /// Observe an explicitly initially resting causal source on the original
    /// integer clock. By calling this method the caller declares Q(t)=0 for
    /// t<epoch and no initial exterior/material excitation. Q(0), Qdot(0) and
    /// its supplied integral MUST be zero; their missing distributional onset
    /// terms are not fabricated. A finite Qddot(0) is allowed, but a jump does
    /// not receive the smooth-onset second-order accuracy claim.
    ///
    /// The callback supplies one internally consistent physical history and
    /// may be queried out of order for different paths. All queries lie in the
    /// declared finite prefix. Before a path's exact r/c_inf wavefront the
    /// output is zero by the caller's explicit rest condition. No phase delay
    /// is added a second time, and no source or received tail is extrapolated.
    /// Cancellation publishes no partial record and does not mutate this cache.
    pub fn observe_initially_resting(
        &self,
        clock: QuadrupoleSampleClock,
        mut source_at: impl FnMut(usize, Dd) -> Result<TransientQuadrupoleSource, AeroacError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularQuadrupoleRecord, AeroacError> {
        poll(&mut cancelled)?;
        if clock.timestep_s() != self.config.timestep_s {
            return Err(invalid(
                "molecular transient source and receiver timestep mismatch",
            ));
        }
        clock.validate_last_index(self.config.samples - 1)?;
        let epoch = Dd::from_f64(clock.epoch_s());
        for source in 0..self.geometry.source_count() {
            poll(&mut cancelled)?;
            if !compatible_initial(&source_at(source, epoch)?) {
                return Err(invalid(
                    "molecular history needs original zero initial stress, rate and integral",
                ));
            }
        }
        let samples = self.config.samples;
        let nfft = self.config.transform_len;
        let h = Dd::from_f64(clock.timestep_s());
        let mut sums = vec![[Dd::ZERO; 4]; samples];
        let mut signals = vec![[[0.0; 4]; 4]; samples];
        let mut scratch = vec![FftComplex::default(); nfft];
        let mut input = vec![FftComplex::default(); nfft];
        let mut output = vec![FftComplex::default(); nfft];
        for (path_index, path) in self.geometry.transfers.iter().enumerate() {
            poll(&mut cancelled)?;
            let source_index = path_index % self.geometry.source_count();
            let delay = Dd::from_f64(path.delay_s);
            let mut first_active = samples;
            for (index, signal) in signals.iter_mut().enumerate() {
                if index.is_multiple_of(64) {
                    poll(&mut cancelled)?;
                }
                *signal = [[0.0; 4]; 4];
                // Subtract the exact frozen-speed delay in elapsed DD time,
                // avoiding a large-epoch add/subtract round-trip in the history.
                let elapsed = h * Dd::from_f64(index as f64);
                let retarded = elapsed - delay;
                if retarded.lt(Dd::ZERO) {
                    continue;
                }
                let receiver = epoch + elapsed;
                let emission = receiver - delay;
                if !emission.hi.is_finite() || !emission.lo.is_finite() || !emission.lt(receiver) {
                    return Err(invalid(
                        "molecular path retardation is unresolved on the DD clock",
                    ));
                }
                *signal = inputs(path, source_at(source_index, epoch + retarded)?)?;
                if first_active == samples && signal.iter().flatten().any(|&v| v != 0.0) {
                    first_active = index;
                }
            }
            if first_active == samples {
                continue;
            }
            if let Some(kernels) = &self.kernels {
                for field in 0..4 {
                    output.fill(FftComplex::default());
                    for (order, kernel) in kernels[path_index].iter().enumerate() {
                        poll(&mut cancelled)?;
                        input.fill(FftComplex::default());
                        let mut active = false;
                        for (index, signal) in signals.iter().enumerate() {
                            input[index].re = signal[order][field];
                            active |= input[index].re != 0.0;
                        }
                        if !active {
                            continue;
                        }
                        self.fft.forward(&mut input, &mut scratch);
                        for index in 0..nfft {
                            if index.is_multiple_of(256) {
                                poll(&mut cancelled)?;
                            }
                            let (a, b, acc) = (input[index], kernel[index], output[index]);
                            output[index] = FftComplex::new(
                                acc.re + a.re.mul_add(b.re, -(a.im * b.im)),
                                acc.im + a.re.mul_add(b.im, a.im * b.re),
                            );
                        }
                    }
                    poll(&mut cancelled)?;
                    self.fft.inverse(&mut output, &mut scratch);
                    // Exact physical/front silence also avoids displaying FFT
                    // roundoff from future samples before a known causal onset.
                    for index in first_active..samples {
                        if !output[index].re.is_finite() || !output[index].im.is_finite() {
                            return Err(nonfinite("molecular transient convolution result"));
                        }
                        sums[index][field] = sums[index][field] + Dd::from_f64(output[index].re);
                    }
                }
            } else {
                // Every F_m is exactly identity for no active reservoir. Reuse
                // the complete physical near/far field without any FFT error.
                for index in first_active..samples {
                    for row in signals[index] {
                        for field in 0..4 {
                            sums[index][field] = sums[index][field] + Dd::from_f64(row[field]);
                        }
                    }
                }
            }
        }
        let mut fields = Vec::with_capacity(samples);
        for (index, sum) in sums.into_iter().enumerate() {
            if index.is_multiple_of(128) {
                poll(&mut cancelled)?;
            }
            let [pressure_pa, ux, uy, uz] = sum.map(Dd::to_f64);
            let particle_velocity_m_s = [ux, uy, uz];
            let instantaneous_intensity_w_m2 = particle_velocity_m_s.map(|u| pressure_pa * u);
            if !pressure_pa.is_finite()
                || particle_velocity_m_s
                    .iter()
                    .chain(&instantaneous_intensity_w_m2)
                    .any(|v| !v.is_finite())
            {
                return Err(nonfinite(
                    "molecular transient accumulated field or intensity",
                ));
            }
            fields.push(RetardedQuadrupoleField {
                pressure_pa,
                particle_velocity_m_s,
                instantaneous_intensity_w_m2,
                scope: SCOPE,
            });
        }
        poll(&mut cancelled)?;
        Ok(MolecularQuadrupoleRecord { clock, fields })
    }

    /// Bind the original geometry, physical tensor order and clock of a finite
    /// native stress record. This explicitly declares rest before its first
    /// sample; incompatible authored initial conditions refuse. The recorded
    /// spline remains an offline finite-bandwidth reconstruction, not a causal
    /// predictor. The entire requested output prefix must fit inside the
    /// original record's sample count; later echoes/tails remain unspecified.
    pub fn observe_history(
        &self,
        history: &PreparedQuadrupoleHistory,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularQuadrupoleRecord, AeroacError> {
        if history.source_count() != self.geometry.source_count()
            || history.sample_count() < self.config.samples
            || history
                .source_positions_m()
                .iter()
                .enumerate()
                .any(|(i, &x)| self.geometry.source_position_m(i) != Some(x))
        {
            return Err(invalid(
                "molecular observer and original recorded source layout/horizon mismatch",
            ));
        }
        let epoch = Dd::from_f64(history.clock().epoch_s());
        self.observe_initially_resting(
            history.clock(),
            |index, time| {
                if time == epoch {
                    history.initial_source(index)
                } else {
                    history.source_at_clock_dd(index, time)
                }
            },
            cancelled,
        )
    }
}
