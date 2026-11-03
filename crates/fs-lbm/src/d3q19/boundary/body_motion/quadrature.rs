//! First spatial moments of the existing uniform midpoint occupancy samples.
//!
//! Integer offsets avoid cancellation in full/symmetric cells. These moments
//! integrate a piecewise-constant cell momentum or force over the same occupied
//! subvolume as the mass fraction; they do not alter the sampling rule.

/// Sum of twice-the-sample-count-scaled offsets from the cell center.
/// The caller supplies the admitted sample count in `1..=8` on each axis.
#[derive(Default)]
pub(super) struct FirstMoment3 {
    offset_sum: [i32; 3],
}

impl FirstMoment3 {
    /// Include one occupied midpoint, with each sample index in `0..samples`.
    pub(super) fn add(&mut self, sample: [usize; 3], samples: usize) {
        for (sum, index) in self.offset_sum.iter_mut().zip(sample) {
            *sum += 2 * index as i32 + 1 - samples as i32;
        }
    }

    /// `integral_inside (position - cell_center) dV` in lattice length^4.
    /// Each sample represents volume `1/samples^3` in the unit cell.
    pub(super) fn finish(self, samples: usize) -> [f64; 3] {
        let denominator = (2 * samples.pow(4)) as f64;
        self.offset_sum.map(|sum| f64::from(sum) / denominator)
    }
}

/// Angular momentum or force moment of a uniform cell density over its
/// occupied subvolume, about an arbitrary fixed origin. `arm` is from that
/// origin to the CELL center; `first_moment` is relative to the cell center.
#[inline]
pub(super) fn angular_integral(
    fraction: f64,
    first_moment: [f64; 3],
    arm: [f64; 3],
    density: [f64; 3],
) -> [f64; 3] {
    core::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        fraction * (arm[j] * density[k] - arm[k] * density[j]) + first_moment[j] * density[k]
            - first_moment[k] * density[j]
    })
}

/// Native tile/lane ordering for a row-major cell key. The grid owner admits
/// dimensions that are positive multiples of `tile_side`.
pub(super) fn cell_order(index: usize, dims: [usize; 3], tile_side: usize) -> (usize, usize) {
    let x = index % dims[0];
    let y = index / dims[0] % dims[1];
    let z = index / (dims[0] * dims[1]);
    let tile = (z / tile_side * (dims[1] / tile_side) + y / tile_side) * (dims[0] / tile_side)
        + x / tile_side;
    let lane = (z % tile_side * tile_side + y % tile_side) * tile_side + x % tile_side;
    (tile, lane)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn integrate(samples: usize, occupied: impl Fn([usize; 3]) -> bool) -> (f64, [f64; 3]) {
        let mut count = 0;
        let mut moment = FirstMoment3::default();
        for z in 0..samples {
            for y in 0..samples {
                for x in 0..samples {
                    let sample = [x, y, z];
                    if occupied(sample) {
                        count += 1;
                        moment.add(sample, samples);
                    }
                }
            }
        }
        (count as f64 / samples.pow(3) as f64, moment.finish(samples))
    }

    fn close(actual: [f64; 3], expected: [f64; 3]) {
        for i in 0..3 {
            assert!(
                (actual[i] - expected[i]).abs() < 2e-14,
                "{actual:?} != {expected:?}"
            );
        }
    }

    #[test]
    fn empty_full_and_symmetric_cells_have_exactly_zero_first_moment() {
        for samples in 1..=8 {
            assert_eq!(integrate(samples, |_| false), (0.0, [0.0; 3]));
            assert_eq!(integrate(samples, |_| true), (1.0, [0.0; 3]));
            let (_, moment) = integrate(samples, |[x, y, z]| {
                (x == 0 && y == 0 && z == 0)
                    || (x + 1 == samples && y + 1 == samples && z + 1 == samples)
            });
            assert_eq!(moment, [0.0; 3]);
        }
    }

    #[test]
    fn asymmetric_half_cell_has_its_real_moment_arm() {
        for samples in [2, 4, 6, 8] {
            let (fraction, first) = integrate(samples, |[x, _, _]| x < samples / 2);
            assert_eq!(fraction, 0.5);
            assert_eq!(first, [-0.125, 0.0, 0.0]);
            // Left half of [-.5,.5]^3 has centroid x=-.25. With a
            // y-directed density of 2, its angular momentum is -0.25 z.
            // The old fraction-times-cell-center expression yields zero.
            assert_eq!(
                angular_integral(fraction, first, [0.0; 3], [0.0, 2.0, 0.0]),
                [0.0, 0.0, -0.25]
            );
        }
    }

    #[test]
    fn angular_integral_matches_direct_sample_moments_and_origin_shift() {
        let samples = 5;
        let arm = [2.1, -0.7, 1.3];
        let density = [-0.3, 0.4, 0.2];
        let occupied = |[x, y, z]: [usize; 3]| 2 * x + y + 3 * z < 11;
        let (fraction, first) = integrate(samples, occupied);
        let mut direct = [0.0; 3];
        for z in 0..samples {
            for y in 0..samples {
                for x in 0..samples {
                    if occupied([x, y, z]) {
                        let r: [f64; 3] = core::array::from_fn(|i| {
                            arm[i] + ([x, y, z][i] as f64 + 0.5) / samples as f64 - 0.5
                        });
                        let weighted = density.map(|p| p / samples.pow(3) as f64);
                        direct[0] += r[1] * weighted[2] - r[2] * weighted[1];
                        direct[1] += r[2] * weighted[0] - r[0] * weighted[2];
                        direct[2] += r[0] * weighted[1] - r[1] * weighted[0];
                    }
                }
            }
        }
        let angular = angular_integral(fraction, first, arm, density);
        close(angular, direct);
        let shift = [0.4, -1.1, 0.2];
        let shifted = angular_integral(
            fraction,
            first,
            core::array::from_fn(|i| arm[i] - shift[i]),
            density,
        );
        let p = density.map(|p| fraction * p);
        close(
            shifted,
            [
                angular[0] - (shift[1] * p[2] - shift[2] * p[1]),
                angular[1] - (shift[2] * p[0] - shift[0] * p[2]),
                angular[2] - (shift[0] * p[1] - shift[1] * p[0]),
            ],
        );
    }

    #[test]
    fn complementary_subvolumes_recover_full_cell_momentum_and_torque() {
        for samples in 1..=8 {
            let selected = |[x, y, z]: [usize; 3]| x + 2 * y + z < samples;
            let (fraction_a, first_a) = integrate(samples, selected);
            let (fraction_b, first_b) = integrate(samples, |point| !selected(point));
            close(core::array::from_fn(|i| first_a[i] + first_b[i]), [0.0; 3]);
            let arm = [0.2, -0.9, 1.7];
            let density = [0.03, 0.01, -0.02];
            let a = angular_integral(fraction_a, first_a, arm, density);
            let b = angular_integral(fraction_b, first_b, arm, density);
            close(
                core::array::from_fn(|i| a[i] + b[i]),
                angular_integral(1.0, [0.0; 3], arm, density),
            );
            assert!((fraction_a + fraction_b - 1.0).abs() < 2e-16);
        }
    }

    #[test]
    fn sparse_moments_follow_native_population_order_across_tile_boundaries() {
        let mut cells = [1, 4, 8, 15, 64, 63, 56, 5, 91];
        cells.sort_unstable_by_key(|&index| cell_order(index, [8; 3], 4));
        assert_eq!(cells, [1, 8, 64, 91, 4, 5, 15, 56, 63]);
    }
}
