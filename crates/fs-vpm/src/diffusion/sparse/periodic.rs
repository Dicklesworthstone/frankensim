//! Exact minimum-image geometry for compact PSE on a rectangular torus.

use crate::{VpmError, invalid_scalar};

/// Rectangular periodic fundamental cell, in metres, with origin `[0, 0]`.
///
/// Coordinates belong to `[0,Lx) × [0,Ly)`; the caller explicitly wraps or
/// recentres external geometry before preparation. The cell is immutable and
/// does not imply a periodic velocity or pressure solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicBox2 {
    lengths: [f64; 2],
}

impl PeriodicBox2 {
    /// Admit two finite positive periods [m].
    ///
    /// # Errors
    /// Refuses zero, negative or nonfinite lengths. Compatibility with a
    /// compact support radius is checked when the graph is prepared.
    pub fn new(lengths_m: [f64; 2]) -> Result<Self, VpmError> {
        for length in lengths_m {
            if !length.is_finite() || length <= 0.0 {
                return Err(invalid_scalar("periodic length", length, None));
            }
        }
        Ok(Self { lengths: lengths_m })
    }

    /// Original x/y periods [m].
    #[must_use]
    pub const fn lengths_m(self) -> [f64; 2] {
        self.lengths
    }

    pub(super) fn grid(self, radius: f64) -> Result<PeriodicGrid, VpmError> {
        let mut counts = [0; 2];
        for (axis, count) in counts.iter_mut().enumerate() {
            // Radius <= L/2 makes the compact periodized kernel equal to its
            // minimum image. Equality is safe because the kernel vanishes
            // exactly at support, including both images on the half-period.
            if radius > self.lengths[axis] * 0.5 {
                return Err(invalid_scalar(
                    "periodic support exceeds half period",
                    radius,
                    None,
                ));
            }
            let ratio = self.lengths[axis] / radius;
            if !ratio.is_finite() || !(2.0..=1_125_899_906_842_624.0).contains(&ratio) {
                return Err(invalid_scalar("periodic cell count", ratio, None));
            }
            *count = ratio.floor() as i64;
        }
        Ok(PeriodicGrid {
            domain: self,
            counts,
        })
    }

    pub(super) fn scaled_distance2(self, a: [f64; 2], b: [f64; 2], width: f64) -> f64 {
        let mut distance2 = 0.0;
        for axis in 0..2 {
            let direct = (a[axis] - b[axis]).abs();
            let separation = if direct <= self.lengths[axis] * 0.5 {
                direct
            } else {
                // Sum distances to the seam instead of L - |a-b|, which
                // loses a small but represented distance near opposite ends.
                (self.lengths[axis] - a[axis].max(b[axis])) + a[axis].min(b[axis])
            };
            let scaled = separation / width;
            distance2 += scaled * scaled;
        }
        distance2
    }
}

#[derive(Clone, Copy)]
pub(super) struct PeriodicGrid {
    pub(super) domain: PeriodicBox2,
    counts: [i64; 2],
}

impl PeriodicGrid {
    pub(super) fn cell(self, position: [f64; 2], particle: usize) -> Result<[i64; 2], VpmError> {
        let mut cell = [0; 2];
        for axis in 0..2 {
            if !(0.0..self.domain.lengths[axis]).contains(&position[axis]) {
                return Err(invalid_scalar(
                    "position outside periodic box",
                    position[axis],
                    Some(particle),
                ));
            }
            // Cells tile the full period, and L/floor(L/radius) >= radius.
            // At the upper edge the quotient may round up to count. Clamp
            // that single rounding case to the last cell, not across the seam.
            cell[axis] = ((position[axis] / self.domain.lengths[axis] * self.counts[axis] as f64)
                .floor() as i64)
                .min(self.counts[axis] - 1);
        }
        Ok(cell)
    }

    pub(super) fn halo_end(self) -> [i64; 2] {
        // The original two-cell rounding halo remains, but each residue is
        // visited once when the period has fewer than five cells.
        self.counts.map(|count| -2 + count.min(5))
    }

    pub(super) fn wrap(self, cell: [i64; 2]) -> [i64; 2] {
        core::array::from_fn(|axis| cell[axis].rem_euclid(self.counts[axis]))
    }
}
