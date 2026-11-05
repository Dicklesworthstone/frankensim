//! Bounded in-place transforms for short real-time convolution partitions.
// Owner-authorized backport of FrankenPiano receiver_fft::Plan (daf1eb708131).
use crate::C64;

/// Cold-planned radix-2 FFT for lengths 16..=256, without runtime allocation.
///
/// Uses ordinary multiply/add arithmetic, not fused multiply-add. Its rounding
/// differs from [`crate::Fft`]; it does not implement that plan's transform-bit
/// contract. Twiddles use the same deterministic fs-math trigonometry. The
/// receiver checks numerical equivalence against direct convolution separately.
#[derive(Clone, Debug)]
pub struct SmallFft {
    roots: Vec<C64>,
    reverse: Vec<usize>,
}

impl SmallFft {
    /// Plan a power-of-two size in 16..=256.
    ///
    /// # Panics
    /// If the requested size is outside this bounded transform's domain.
    #[must_use]
    pub fn new(n: usize) -> Self {
        assert!((16..=256).contains(&n) && n.is_power_of_two());
        let roots = (0..n / 2)
            .map(|i| {
                let angle = -std::f64::consts::TAU * i as f64 / n as f64;
                C64 {
                    re: fs_math::det::cos(angle),
                    im: fs_math::det::sin(angle),
                }
            })
            .collect();
        let reverse = (0..n)
            .map(|i| i.reverse_bits() >> (usize::BITS - n.trailing_zeros()))
            .collect();
        Self { roots, reverse }
    }

    /// Forward DFT, unnormalized, with the negative-exponent convention.
    ///
    /// # Panics
    /// If `data` does not have the planned length.
    pub fn forward(&self, data: &mut [C64]) {
        self.transform(data, false);
    }

    /// Inverse DFT, normalized by the planned length.
    ///
    /// # Panics
    /// If `data` does not have the planned length.
    pub fn inverse(&self, data: &mut [C64]) {
        self.transform(data, true);
    }

    fn transform(&self, data: &mut [C64], inverse: bool) {
        let n = self.reverse.len();
        assert_eq!(data.len(), n);
        for (i, &j) in self.reverse.iter().enumerate() {
            if i < j {
                data.swap(i, j);
            }
        }
        for power in 1..=n.trailing_zeros() {
            let size = 1usize << power;
            for chunk in data.chunks_exact_mut(size) {
                let (left, right) = chunk.split_at_mut(size / 2);
                for (k, (a, b)) in left.iter_mut().zip(right).enumerate() {
                    let mut w = self.roots[k * (n / size)];
                    if inverse {
                        w.im = -w.im;
                    }
                    let z = C64 {
                        re: b.re * w.re - b.im * w.im,
                        im: b.re * w.im + b.im * w.re,
                    };
                    let old = *a;
                    *a = C64 {
                        re: old.re + z.re,
                        im: old.im + z.im,
                    };
                    *b = C64 {
                        re: old.re - z.re,
                        im: old.im - z.im,
                    };
                }
            }
        }
        if inverse {
            for z in data {
                z.re /= n as f64;
                z.im /= n as f64;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn small_plan_matches_general_transform_and_round_trips_every_size() {
        for n in [16, 32, 64, 128, 256] {
            let plan = SmallFft::new(n);
            let reference = crate::Fft::new(n);
            let input: Vec<_> = (0..n)
                .map(|i| C64 {
                    re: (i as f64 * 0.17).sin(),
                    im: (i as f64 * 0.31).cos(),
                })
                .collect();
            let mut actual = input.clone();
            let mut expected = input.clone();
            let mut scratch = vec![C64::default(); n];
            plan.forward(&mut actual);
            reference.forward(&mut expected, &mut scratch);
            for (a, b) in actual.iter().zip(&expected) {
                assert!((a.re - b.re).abs() < 1e-11 && (a.im - b.im).abs() < 1e-11);
            }
            plan.inverse(&mut actual);
            for (a, b) in actual.iter().zip(&input) {
                assert!((a.re - b.re).abs() < 1e-12 && (a.im - b.im).abs() < 1e-12);
            }
        }
    }
}
