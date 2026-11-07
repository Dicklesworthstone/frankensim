//! Causal MIMO observation FIR data and replay. These filters never return a
//! reaction force or assert positive realness: observation is outside the loop.
use fs_math::{
    c64::C64,
    det,
    eft::{two_prod, two_sum},
};

/// Compensated dot product (Ogita--Rump--Oishi Dot2). Both product roundoff
/// and addition roundoff survive destructive interference. Keeping this pair
/// across kernels matters: rounding a source's partial sum first can erase
/// the quiet remainder of two large opposing sources.
#[derive(Clone, Copy, Default)]
struct Dot {
    sum: f64,
    correction: f64,
}

impl Dot {
    fn add_product(&mut self, a: f64, b: f64) {
        let (product, product_error) = two_prod(a, b);
        let (sum, sum_error) = two_sum(self.sum, product);
        self.sum = sum;
        self.correction += product_error + sum_error;
    }

    fn total(self) -> f64 {
        self.sum + self.correction
    }
}

/// Total retained coefficients and shared input history, each separately capped.
pub const MAX_OBSERVATION_ENTRIES: usize = 4_194_304;
/// Complete input basis cap; never silently truncate modes/apertures.
pub const MAX_OBSERVATION_INPUTS: usize = 4096;
/// Output components cap.
pub const MAX_OBSERVATION_OUTPUTS: usize = 256;

/// One source-to-receiver-component kernel. Fractional propagation belongs to
/// the taps, not to an unaccounted allpass. The integer prefix is exact.
#[derive(Debug, Clone, PartialEq)]
pub struct FirKernel {
    pub input: usize,
    pub output: usize,
    pub integer_delay_samples: usize,
    pub taps: Vec<f64>,
}
impl FirKernel {
    /// DFT convention exp(-i*omega*t), including the integer propagation delay.
    pub fn transfer_at(&self, omega_rad_per_sample: f64) -> Result<C64, FirError> {
        if !omega_rad_per_sample.is_finite() {
            return Err(FirError::Invalid("nonfinite frequency"));
        }
        let mut real = Dot::default();
        let mut imaginary = Dot::default();
        for (i, &tap) in self.taps.iter().enumerate() {
            let lag = self
                .integer_delay_samples
                .checked_add(i)
                .ok_or(FirError::Invalid("delay length overflow"))?;
            let phase = -omega_rad_per_sample * lag as f64;
            if !phase.is_finite() {
                return Err(FirError::Invalid("phase overflow"));
            }
            real.add_product(tap, det::cos(phase));
            imaginary.add_product(tap, det::sin(phase));
        }
        let h = C64::new(real.total(), imaginary.total());
        if !h.re.is_finite() || !h.im.is_finite() {
            return Err(FirError::Invalid("nonfinite FIR transfer"));
        }
        Ok(h)
    }
}

/// Typed cold admission or transactional replay refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FirError {
    Invalid(&'static str),
    Cap { requested: usize, cap: usize },
}
impl core::fmt::Display for FirError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "observation FIR: {self:?}")
    }
}
impl std::error::Error for FirError {}

/// Allocation-free complete-shape admission, also usable by offline producers.
pub fn admit_observation_work(
    inputs: usize,
    outputs: usize,
    coefficients: usize,
    history_entries: usize,
) -> Result<(), FirError> {
    if inputs == 0 || outputs == 0 || coefficients == 0 || history_entries == 0 {
        return Err(FirError::Invalid("empty observation operator"));
    }
    for (requested, cap) in [
        (inputs, MAX_OBSERVATION_INPUTS),
        (outputs, MAX_OBSERVATION_OUTPUTS),
        (coefficients, MAX_OBSERVATION_ENTRIES),
        (history_entries, MAX_OBSERVATION_ENTRIES),
    ] {
        if requested > cap {
            return Err(FirError::Cap { requested, cap });
        }
    }
    Ok(())
}

/// Real FIR convolution with compensated products/sums, zero prehistory and
/// fixed traversal order. Input histories are shared across all receiver
/// components. No source is omitted. Compensation preserves cancellation
/// residuals but is not an error certificate at arbitrary interference nulls.
#[derive(Clone)]
pub struct ObservationFir {
    inputs: usize,
    outputs: usize,
    kernels: Vec<FirKernel>,
    history: Vec<Vec<f64>>,
    sums: Vec<Dot>,
    cursor: usize,
    length: usize,
}
impl ObservationFir {
    /// Construct from explicit kernels. Duplicate input/output kernels add;
    /// missing pairs are zero only because the producer explicitly supplied
    /// that sparse operator. A BEM bake supplies every pair and component.
    pub fn new(inputs: usize, outputs: usize, kernels: Vec<FirKernel>) -> Result<Self, FirError> {
        let mut length = 0;
        let mut count = 0usize;
        for k in &kernels {
            if k.input >= inputs
                || k.output >= outputs
                || k.taps.is_empty()
                || k.taps.iter().any(|x| !x.is_finite())
            {
                return Err(FirError::Invalid("invalid kernel index or taps"));
            }
            let end = k
                .integer_delay_samples
                .checked_add(k.taps.len())
                .ok_or(FirError::Invalid("delay length overflow"))?;
            length = length.max(end);
            count = count
                .checked_add(k.taps.len())
                .ok_or(FirError::Invalid("tap count overflow"))?;
        }
        let entries = length
            .checked_mul(inputs)
            .ok_or(FirError::Invalid("history shape overflow"))?;
        admit_observation_work(inputs, outputs, count, entries)?;
        Ok(Self {
            inputs,
            outputs,
            kernels,
            history: vec![vec![0.0; length]; inputs],
            sums: vec![Dot::default(); outputs],
            cursor: 0,
            length,
        })
    }
    /// Immutable operator data, suitable for content-addressed artifacts.
    pub fn kernels(&self) -> &[FirKernel] {
        &self.kernels
    }
    /// Transactional one-sample observation. A rejected input/output changes no
    /// history, so retry and replay remain identical after an error.
    pub fn step(&mut self, input: &[f64]) -> Result<Vec<f64>, FirError> {
        let mut output = vec![0.0; self.outputs];
        self.step_into(input, &mut output)?;
        Ok(output)
    }
    /// Allocation-free one-sample observation into exactly the complete output
    /// shape. An error changes neither the caller's output nor any history.
    /// Cold-allocated scratch holds two scalars per admitted output; no
    /// sample-dependent allocation or source/tap truncation occurs.
    pub fn step_into(&mut self, input: &[f64], output: &mut [f64]) -> Result<(), FirError> {
        if input.len() != self.inputs || input.iter().any(|x| !x.is_finite()) {
            return Err(FirError::Invalid("input shape or nonfinite value"));
        }
        if output.len() != self.outputs {
            return Err(FirError::Invalid("output shape"));
        }
        self.sums.fill(Dot::default());
        for kernel in &self.kernels {
            let history = &self.history[kernel.input];
            // Lag l reads slot (cursor - l) mod length, so walk down the ring
            // instead of dividing per tap. Keep compensation across sources.
            let delay = kernel.integer_delay_samples;
            let mut slot = (self.cursor + self.length - delay) % self.length;
            let sum = &mut self.sums[kernel.output];
            for (i, &tap) in kernel.taps.iter().enumerate() {
                let value = if i + delay == 0 {
                    input[kernel.input]
                } else {
                    history[slot]
                };
                sum.add_product(tap, value);
                slot = if slot == 0 { self.length - 1 } else { slot - 1 };
            }
        }
        if self.sums.iter().any(|sum| !sum.total().is_finite()) {
            return Err(FirError::Invalid("output overflow"));
        }
        for (row, &value) in self.history.iter_mut().zip(input) {
            row[self.cursor] = value;
        }
        self.cursor = (self.cursor + 1) % self.length;
        for (value, sum) in output.iter_mut().zip(&self.sums) {
            *value = sum.total();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_causal_multi_input_replay_is_transactional() {
        let kernels = vec![
            FirKernel {
                input: 0,
                output: 0,
                integer_delay_samples: 2,
                taps: vec![1.0, -0.5],
            },
            FirKernel {
                input: 1,
                output: 0,
                integer_delay_samples: 0,
                taps: vec![2.0],
            },
        ];
        let mut a = ObservationFir::new(2, 1, kernels).unwrap();
        let mut b = a.clone();
        assert!(a.step(&[f64::NAN, 0.0]).is_err());
        for n in 0..8 {
            let x = [
                if n == 0 { 1.0 } else { 0.0 },
                if n == 1 { 3.0 } else { 0.0 },
            ];
            let expected = match n {
                1 => 6.0,
                2 => 1.0,
                3 => -0.5,
                _ => 0.0,
            };
            assert_eq!(a.step(&x).unwrap(), vec![expected]);
            assert_eq!(b.step(&x).unwrap(), vec![expected]);
        }
    }
    #[test]
    fn ring_walk_matches_independent_convolution_and_replays_bitwise() {
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        // Zero and long delays, duplicate pairs and every output, over many ring turns.
        let shapes = [
            (0, 1, 0, 5),
            (1, 0, 3, 17),
            (2, 2, 0, 1),
            (0, 2, 41, 9),
            (1, 1, 1, 64),
            (0, 1, 2, 3),
        ];
        let kernels: Vec<FirKernel> = shapes
            .iter()
            .map(|&(input, output, integer_delay_samples, taps)| FirKernel {
                input,
                output,
                integer_delay_samples,
                taps: (0..taps).map(|_| next()).collect(),
            })
            .collect();
        let mut fir = ObservationFir::new(3, 3, kernels.clone()).unwrap();
        let mut replay = fir.clone();
        let length = fir.length;
        let mut past: Vec<Vec<f64>> = Vec::new();
        for _ in 0..7 * length {
            let x: Vec<f64> = (0..3).map(|_| next()).collect();
            past.push(x.clone());
            let mut expected = [fs_math::dd::Dd::ZERO; 3];
            for k in &kernels {
                for (i, &tap) in k.taps.iter().enumerate() {
                    let lag = i + k.integer_delay_samples;
                    let value = past
                        .len()
                        .checked_sub(lag + 1)
                        .map_or(0.0, |n| past[n][k.input]);
                    expected[k.output] = expected[k.output]
                        + fs_math::dd::Dd::from_f64(tap) * fs_math::dd::Dd::from_f64(value);
                }
            }
            let actual = fir.step(&x).unwrap();
            let mut repeated = [0.0; 3];
            replay.step_into(&x, &mut repeated).unwrap();
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&actual), bits(&repeated));
            for (&actual, expected) in actual.iter().zip(expected) {
                let expected = expected.to_f64();
                assert!((actual - expected).abs() <= 2.0 * f64::EPSILON * expected.abs().max(1.0));
            }
        }
    }

    #[test]
    fn g1_destructive_interference_retains_quiet_source_across_kernel_boundaries() {
        let loud = 2.0_f64.powi(54);
        let kernels = vec![
            FirKernel {
                input: 0,
                output: 0,
                integer_delay_samples: 2,
                taps: vec![loud],
            },
            FirKernel {
                input: 1,
                output: 0,
                integer_delay_samples: 0,
                taps: vec![1.0],
            },
            FirKernel {
                input: 0,
                output: 0,
                integer_delay_samples: 2,
                taps: vec![-loud],
            },
        ];
        let mut observer = ObservationFir::new(2, 1, kernels).unwrap();
        // Two equal, oppositely signed delayed radiators cancel. The quiet
        // source must survive before/after the wavefront and every ring turn.
        for n in 0..24 {
            let quiet = if n % 2 == 0 { 1.0 } else { -0.25 };
            assert_eq!(observer.step(&[1.0, quiet]).unwrap(), vec![quiet]);
        }
        let kernel = FirKernel {
            input: 0,
            output: 0,
            integer_delay_samples: 3,
            taps: vec![loud, 1.0, -loud],
        };
        assert_eq!(kernel.transfer_at(0.0).unwrap(), C64::ONE);
        let mut observer = ObservationFir::new(1, 1, vec![kernel]).unwrap();
        for _ in 0..8 {
            let _ = observer.step(&[1.0]).unwrap();
        }
        assert_eq!(observer.step(&[1.0]).unwrap(), vec![1.0]);
    }

    #[test]
    fn g1_cancellation_retains_the_exact_product_remainder() {
        // (1+2^-27)(1-2^-27)-1 = -2^-54. Each rounded product
        // alone is +/-1, so even compensated addition alone gives false zero.
        let epsilon = 2.0_f64.powi(-27);
        let kernels = vec![
            FirKernel {
                input: 0,
                output: 0,
                integer_delay_samples: 0,
                taps: vec![1.0 + epsilon],
            },
            FirKernel {
                input: 1,
                output: 0,
                integer_delay_samples: 0,
                taps: vec![-1.0],
            },
        ];
        let mut observer = ObservationFir::new(2, 1, kernels).unwrap();
        assert_eq!(
            observer.step(&[1.0 - epsilon, 1.0]).unwrap(),
            vec![-epsilon * epsilon]
        );
    }

    #[test]
    fn g0_into_refusal_keeps_output_and_delayed_history_transactional() {
        let kernels = vec![FirKernel {
            input: 0,
            output: 0,
            integer_delay_samples: 0,
            taps: vec![2.0, -1.0],
        }];
        let mut observer = ObservationFir::new(1, 1, kernels).unwrap();
        observer.step(&[0.75]).unwrap();
        let mut reference = observer.clone();
        let mut output = [17.0];
        assert!(observer.step_into(&[f64::NAN], &mut output).is_err());
        assert_eq!(output, [17.0]);
        assert!(observer.step_into(&[f64::MAX], &mut output).is_err());
        assert_eq!(output, [17.0]);
        let mut wrong_shape = [19.0, 23.0];
        assert!(observer.step_into(&[1.0], &mut wrong_shape).is_err());
        assert_eq!(wrong_shape, [19.0, 23.0]);
        for input in [0.25, -0.5, 0.0] {
            observer.step_into(&[input], &mut output).unwrap();
            assert_eq!(output.as_slice(), reference.step(&[input]).unwrap());
        }
    }

    #[test]
    fn g0_observation_caps_and_plus_one() {
        assert!(
            admit_observation_work(4096, 256, MAX_OBSERVATION_ENTRIES, MAX_OBSERVATION_ENTRIES)
                .is_ok()
        );
        for request in [
            (4097, 256, 1, 1),
            (4096, 257, 1, 1),
            (1, 1, MAX_OBSERVATION_ENTRIES + 1, 1),
            (1, 1, 1, MAX_OBSERVATION_ENTRIES + 1),
        ] {
            assert!(matches!(
                admit_observation_work(request.0, request.1, request.2, request.3),
                Err(FirError::Cap { .. })
            ));
        }
        assert!(ObservationFir::new(1, 1, vec![]).is_err());
        assert!(
            ObservationFir::new(
                1,
                1,
                vec![FirKernel {
                    input: 1,
                    output: 0,
                    integer_delay_samples: 0,
                    taps: vec![1.0]
                }]
            )
            .is_err()
        );
    }
}
