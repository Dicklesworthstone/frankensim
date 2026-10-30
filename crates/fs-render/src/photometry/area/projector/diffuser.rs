//! Passive, finite slope-space redistribution of the projector's ray power.
use super::{GoboTransmission, MAX_GOBO_TEXELS, SourceError};
use fs_math::det;
use std::sync::Arc;

/// An ideal achromatic exit-pupil diffuser with independent uniform X/Y slopes.
///
/// Scattered rays add slopes uniform in +/-tan(half_angle) in each axis.
/// The remaining fraction retains its original direction. Transmission is an
/// absolute surviving power fraction. This declared paraxial ray-power model
/// is not a reciprocal surface BSDF or a measured manufactured frost model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParaxialBoxDiffuser {
    half_angles_rad: [f64; 2],
    half_slope_widths: [f64; 2],
    scattered_fraction: f64,
    transmission: f64,
}

impl ParaxialBoxDiffuser {
    /// Admit positive half angles below pi/2 and passive fractions in [0,1].
    /// Widths are retained even when the scattered fraction is zero.
    pub fn try_new(
        half_angles_rad: [f64; 2],
        scattered_fraction: f64,
        transmission: f64,
    ) -> Result<Self, SourceError> {
        if half_angles_rad.iter().any(|&angle| {
            !angle.is_finite() || angle <= 0.0 || angle >= core::f64::consts::FRAC_PI_2
        }) {
            return Err(SourceError::InvalidDistribution);
        }
        if [scattered_fraction, transmission]
            .iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        {
            return Err(SourceError::InvalidFilter);
        }
        let half_slope_widths = half_angles_rad.map(det::tan);
        if half_slope_widths
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            half_angles_rad,
            half_slope_widths,
            scattered_fraction,
            transmission,
        })
    }
    /// Declared positive X/Y half angles in radians.
    #[must_use]
    pub const fn half_angles_rad(&self) -> [f64; 2] {
        self.half_angles_rad
    }
    /// Positive dimensionless half widths of outgoing slope deflections.
    #[must_use]
    pub const fn half_slope_widths(&self) -> [f64; 2] {
        self.half_slope_widths
    }
    /// Fraction of surviving power redistributed over the slope rectangle.
    #[must_use]
    pub const fn scattered_fraction(&self) -> f64 {
        self.scattered_fraction
    }
    /// Absolute surviving power fraction, independent of gobo and splitter loss.
    #[must_use]
    pub const fn transmission(&self) -> f64 {
        self.transmission
    }
}

#[derive(Debug, Clone)]
pub(super) struct PreparedDiffuser {
    pub(super) law: ParaxialBoxDiffuser,
    pub(super) half_width_m: [f64; 2],
    integral: Option<Arc<PositiveIntegral>>,
}

impl PreparedDiffuser {
    pub(super) fn try_new(
        law: ParaxialBoxDiffuser,
        gobo: &GoboTransmission,
        distance_m: f64,
    ) -> Result<Self, SourceError> {
        let active =
            law.scattered_fraction > 0.0 && law.transmission > 0.0 && gobo.data.maximum > 0.0;
        let (half_width_m, integral) = if active {
            let widths = law.half_slope_widths.map(|slope| distance_m * slope);
            for (width, size) in widths.iter().zip([gobo.width_m(), gobo.height_m()]) {
                if !width.is_finite()
                    || *width <= 0.0
                    || !(0.5 * size + *width).next_up().is_finite()
                {
                    return Err(SourceError::NumericalRange);
                }
            }
            (widths, Some(Arc::new(PositiveIntegral::try_new(gobo)?)))
        } else {
            ([0.0; 2], None)
        };
        Ok(Self {
            law,
            half_width_m,
            integral,
        })
    }
    pub(super) fn is_identity(&self) -> bool {
        self.law.scattered_fraction == 0.0 && self.law.transmission == 1.0
    }
    pub(super) fn mask(
        &self,
        gobo: &GoboTransmission,
        x_m: f64,
        y_m: f64,
    ) -> Result<f64, SourceError> {
        if self.law.transmission == 0.0 {
            return Ok(0.0);
        }
        let clear = if self.law.scattered_fraction < 1.0 {
            product(
                1.0 - self.law.scattered_fraction,
                gobo.transmission_at(x_m, y_m)?,
            )?
        } else {
            0.0
        };
        let scattered = match &self.integral {
            Some(integral) => product(
                self.law.scattered_fraction,
                integral.average(gobo, [x_m, y_m], self.half_width_m)?,
            )?,
            None => 0.0,
        };
        product(self.law.transmission, clear + scattered)
    }
}

pub(super) fn product(left: f64, right: f64) -> Result<f64, SourceError> {
    let value = left * right;
    if !value.is_finite() || (left > 0.0 && right > 0.0 && value == 0.0) {
        Err(SourceError::NumericalRange)
    } else {
        Ok(value)
    }
}

/// Tensor segment tree: canonical rectangles add only positive nodes.
/// No large-prefix subtraction can erase a small independently queried texel.
/// Power-of-two padding is charged against four million double storage slots.
#[derive(Debug)]
struct PositiveIntegral {
    row_leaf: usize,
    column_leaf: usize,
    stride: usize,
    sums: Vec<f64>,
}

impl PositiveIntegral {
    fn try_new(gobo: &GoboTransmission) -> Result<Self, SourceError> {
        let row_leaf = gobo.rows().next_power_of_two();
        let column_leaf = gobo.columns().next_power_of_two();
        let length = row_leaf
            .checked_mul(column_leaf)
            .and_then(|n| n.checked_mul(4))
            .ok_or(SourceError::DistributionTooLarge)?;
        if length > 4 * MAX_GOBO_TEXELS {
            return Err(SourceError::DistributionTooLarge);
        }
        if gobo.width_m() / gobo.columns() as f64 == 0.0
            || gobo.height_m() / gobo.rows() as f64 == 0.0
        {
            return Err(SourceError::NumericalRange);
        }
        let stride = 2 * column_leaf;
        let mut sums = Vec::new();
        sums.try_reserve_exact(length)
            .map_err(|_| SourceError::DistributionTooLarge)?;
        sums.resize(length, 0.0);
        for row in 0..gobo.rows() {
            let base = (row_leaf + row) * stride;
            for column in 0..gobo.columns() {
                sums[base + column_leaf + column] = gobo.data.values[row * gobo.columns() + column];
            }
            for column in (1..column_leaf).rev() {
                sums[base + column] = sums[base + 2 * column] + sums[base + 2 * column + 1];
            }
        }
        for row in (1..row_leaf).rev() {
            for column in 1..stride {
                sums[row * stride + column] =
                    sums[2 * row * stride + column] + sums[(2 * row + 1) * stride + column];
            }
        }
        Ok(Self {
            row_leaf,
            column_leaf,
            stride,
            sums,
        })
    }
    fn rectangle_sum(&self, x: [usize; 2], y: [usize; 2]) -> f64 {
        let (mut top, mut bottom) = (y[0] + self.row_leaf, y[1] + self.row_leaf);
        let (mut sum, mut correction) = (0.0, 0.0);
        while top < bottom {
            if top & 1 != 0 {
                self.add_row(top, x, &mut sum, &mut correction);
                top += 1;
            }
            if bottom & 1 != 0 {
                bottom -= 1;
                self.add_row(bottom, x, &mut sum, &mut correction);
            }
            top /= 2;
            bottom /= 2;
        }
        sum + correction
    }
    fn add_row(&self, row: usize, x: [usize; 2], sum: &mut f64, correction: &mut f64) {
        let (mut left, mut right) = (x[0] + self.column_leaf, x[1] + self.column_leaf);
        let base = row * self.stride;
        while left < right {
            if left & 1 != 0 {
                add(sum, correction, self.sums[base + left]);
                left += 1;
            }
            if right & 1 != 0 {
                right -= 1;
                add(sum, correction, self.sums[base + right]);
            }
            left /= 2;
            right /= 2;
        }
    }
    fn average(
        &self,
        gobo: &GoboTransmission,
        center: [f64; 2],
        half_width: [f64; 2],
    ) -> Result<f64, SourceError> {
        let x = axis_ranges(gobo.width_m(), gobo.columns(), center[0], half_width[0])?;
        let y = axis_ranges(gobo.height_m(), gobo.rows(), center[1], half_width[1])?;
        let (mut sum, mut correction) = (0.0, 0.0);
        for xr in x {
            if xr.first == xr.end {
                continue;
            }
            for yr in y {
                if yr.first == yr.end {
                    continue;
                }
                let mass = self.rectangle_sum([xr.first, xr.end], [yr.first, yr.end]);
                let value = product(product(mass, xr.weight)?, yr.weight)?;
                add(&mut sum, &mut correction, value);
            }
        }
        let value = sum + correction;
        if !value.is_finite() {
            return Err(SourceError::NumericalRange);
        }
        Ok(value)
    }
}

#[derive(Clone, Copy, Default)]
struct WeightedRange {
    first: usize,
    end: usize,
    // Fraction of the complete kernel interval assigned to each full cell.
    weight: f64,
}

fn axis_ranges(
    size: f64,
    count: usize,
    center: f64,
    half: f64,
) -> Result<[WeightedRange; 3], SourceError> {
    if !center.is_finite() || !half.is_finite() || half <= 0.0 {
        return Err(SourceError::NumericalRange);
    }
    let boundary = |index: usize| (index as f64 / count as f64 - 0.5) * size;
    // Comparing distances avoids collapsing narrow positive kernels by first
    // rounding centre +/- half at an otherwise ordinary world coordinate.
    let (mut lo, mut hi) = (0, count);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if boundary(mid + 1) - center <= -half {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let first = lo;
    let (mut lo, mut hi) = (first, count);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if boundary(mid) - center < half {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let end = lo;
    let mut result = [WeightedRange::default(); 3];
    if first == end {
        return Ok(result);
    }
    let fraction = |cell: usize| {
        let left = ((boundary(cell) - center) / half).clamp(-1.0, 1.0);
        let right = ((boundary(cell + 1) - center) / half).clamp(-1.0, 1.0);
        let value = 0.5 * (right - left);
        if !value.is_finite() || value <= 0.0 {
            Err(SourceError::NumericalRange)
        } else {
            Ok(value)
        }
    };
    result[0] = WeightedRange {
        first,
        end: first + 1,
        weight: fraction(first)?,
    };
    if end - first == 1 {
        return Ok(result);
    }
    result[2] = WeightedRange {
        first: end - 1,
        end,
        weight: fraction(end - 1)?,
    };
    if end - first > 2 {
        let weight = 0.5 * ((size / count as f64) / half);
        if !weight.is_finite() || weight <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        result[1] = WeightedRange {
            first: first + 1,
            end: end - 1,
            weight,
        };
    }
    Ok(result)
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let next = *sum + value;
    *correction += if *sum >= value {
        (*sum - next) + value
    } else {
        (value - next) + *sum
    };
    *sum = next;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::RigidTransform;
    use crate::photometry::area::{ParaxialBeamSplitter, ProjectorBeam, ThinLensProjector};
    use fs_geom::{Point3, Vec3};

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2.0e-14 * expected.abs().max(1.0),
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }
    fn gobo() -> GoboTransmission {
        GoboTransmission::try_new(
            4.0,
            2.0,
            4,
            2,
            vec![0.0, 0.25, 0.5, 1.0, 1.0, 0.5, 0.25, 0.0],
        )
        .unwrap()
    }

    #[test]
    fn g0_box_diffuser_admits_only_finite_positive_widths_and_passive_fractions() {
        let law = ParaxialBoxDiffuser::try_new([0.1, 0.2], 0.75, 0.8).unwrap();
        assert_eq!(law.half_angles_rad(), [0.1, 0.2]);
        close(law.half_slope_widths()[0], 0.1_f64.tan());
        close(law.half_slope_widths()[1], 0.2_f64.tan());
        assert_eq!(law.scattered_fraction(), 0.75);
        assert_eq!(law.transmission(), 0.8);
        for angle in [
            0.0,
            -0.1,
            core::f64::consts::FRAC_PI_2,
            f64::INFINITY,
            f64::NAN,
        ] {
            assert!(ParaxialBoxDiffuser::try_new([angle, 0.1], 1.0, 1.0).is_err());
            assert!(ParaxialBoxDiffuser::try_new([0.1, angle], 1.0, 1.0).is_err());
        }
        for fraction in [-0.1, 1.1, f64::INFINITY, f64::NAN] {
            assert!(ParaxialBoxDiffuser::try_new([0.1; 2], fraction, 1.0).is_err());
            assert!(ParaxialBoxDiffuser::try_new([0.1; 2], 1.0, fraction).is_err());
        }
        assert!(ParaxialBoxDiffuser::try_new([0.1; 2], 0.0, 0.0).is_ok());
    }

    #[test]
    fn g1_positive_rectangles_match_independent_texel_overlap_and_keep_dark_exterior() {
        let gobo = gobo();
        let integral = PositiveIntegral::try_new(&gobo).unwrap();
        close(
            integral.average(&gobo, [0.0; 2], [2.0, 1.0]).unwrap(),
            0.4375,
        );
        close(
            integral.average(&gobo, [0.0; 2], [4.0, 2.0]).unwrap(),
            0.4375 / 4.0,
        );
        close(
            integral.average(&gobo, [2.0, 0.0], [0.5, 1.0]).unwrap(),
            0.25,
        );
        assert_eq!(integral.average(&gobo, [3.0, 0.0], [0.5; 2]).unwrap(), 0.0);
        for row in 0..2 {
            for column in 0..4 {
                let center = [column as f64 - 1.5, row as f64 - 0.5];
                assert_eq!(
                    integral.average(&gobo, center, [0.5; 2]).unwrap(),
                    gobo.data.values[row * 4 + column]
                );
            }
        }
        let gobo = GoboTransmission::try_new(
            3.0,
            2.0,
            5,
            3,
            (0..15).map(|index| (index % 7) as f64 / 6.0).collect(),
        )
        .unwrap();
        let integral = PositiveIntegral::try_new(&gobo).unwrap();
        for index in 0..37 {
            let center = [index as f64 * 0.137 - 2.0, index as f64 * 0.071 - 1.3];
            let half = [
                0.23 + (index % 4) as f64 * 0.31,
                0.19 + (index % 3) as f64 * 0.47,
            ];
            let mut expected = 0.0;
            for row in 0..3 {
                for column in 0..5 {
                    let mut weight = 1.0;
                    for (axis, cell, count, size) in [(0, column, 5, 3.0), (1, row, 3, 2.0)] {
                        let low: f64 = -size / 2.0 + cell as f64 * size / count as f64;
                        let high: f64 = -size / 2.0 + (cell + 1) as f64 * size / count as f64;
                        let overlap = (high.min(center[axis] + half[axis])
                            - low.max(center[axis] - half[axis]))
                        .max(0.0);
                        weight *= overlap / (2.0 * half[axis]);
                    }
                    expected += gobo.data.values[row * 5 + column] * weight;
                }
            }
            close(integral.average(&gobo, center, half).unwrap(), expected);
        }
    }

    #[test]
    fn g3_positive_tree_preserves_tiny_cells_and_sub_ulp_kernel_widths() {
        let gobo = GoboTransmission::try_new(2.0, 1.0, 2, 1, vec![1.0, 1.0e-280]).unwrap();
        let integral = PositiveIntegral::try_new(&gobo).unwrap();
        assert_eq!(
            integral.average(&gobo, [0.5, 0.0], [0.4; 2]).unwrap(),
            1.0e-280
        );
        assert_eq!(
            integral.average(&gobo, [0.5, 0.1], [1.0e-300; 2]).unwrap(),
            1.0e-280
        );
        assert_eq!(
            integral.average(&gobo, [0.0, 0.1], [1.0e-300; 2]).unwrap(),
            0.5
        );
        close(
            integral.average(&gobo, [1.0, 0.0], [1.0e-300; 2]).unwrap() / 1.0e-280,
            0.5,
        );
        assert!(matches!(
            integral.average(&gobo, [0.5, 0.0], [1.0e300; 2]),
            Err(SourceError::NumericalRange)
        ));
    }

    #[test]
    fn g0_frost_power_composes_and_replacement_does_not_compound_losses() {
        let base =
            ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 0.5, gobo()).unwrap();
        let splitter = ParaxialBeamSplitter::try_new(vec![
            ProjectorBeam::try_new([-0.1, 0.0], 0.4).unwrap(),
            ProjectorBeam::try_new([0.1, 0.0], 0.3).unwrap(),
        ])
        .unwrap();
        let frost = ParaxialBoxDiffuser::try_new([0.2, 0.3], 0.8, 0.6).unwrap();
        let first = base
            .clone()
            .with_beam_splitter(splitter.clone())
            .unwrap()
            .with_diffuser(frost)
            .unwrap();
        let other_order = base
            .clone()
            .with_diffuser(frost)
            .unwrap()
            .with_beam_splitter(splitter)
            .unwrap();
        close(first.transmitted_fraction(), 0.4375 * 0.7 * 0.6);
        assert_eq!(
            first.transmitted_fraction().to_bits(),
            other_order.transmitted_fraction().to_bits()
        );
        assert_eq!(first.diffuser(), Some(&frost));
        let replaced = first
            .with_diffuser(ParaxialBoxDiffuser::try_new([0.1; 2], 0.2, 0.5).unwrap())
            .unwrap();
        close(replaced.transmitted_fraction(), 0.4375 * 0.7 * 0.5);
        let replaced = replaced
            .with_beam_splitter(
                ParaxialBeamSplitter::try_new(vec![ProjectorBeam::try_new([0.0; 2], 0.2).unwrap()])
                    .unwrap(),
            )
            .unwrap();
        close(replaced.transmitted_fraction(), 0.4375 * 0.2 * 0.5);
        let opaque = base
            .with_diffuser(ParaxialBoxDiffuser::try_new([0.1; 2], 1.0, 0.0).unwrap())
            .unwrap();
        assert_eq!(opaque.transmitted_fraction(), 0.0);
        assert_eq!(
            opaque
                .luminance_per_clear_flux_per_area(
                    Point3::new(0.0, 0.0, 0.0),
                    Vec3::new(0.0, 0.0, 1.0)
                )
                .unwrap(),
            0.0
        );
    }

    #[test]
    fn g1_identity_frost_keeps_bits_and_scattered_support_obeys_the_radiance_bound() {
        let base =
            ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 0.5, gobo()).unwrap();
        for base in [
            base.clone(),
            base.with_beam_splitter(
                ParaxialBeamSplitter::try_new(vec![
                    ProjectorBeam::try_new([0.5, -0.2], 0.6).unwrap(),
                    ProjectorBeam::try_new([-0.3, 0.1], 0.3).unwrap(),
                ])
                .unwrap(),
            )
            .unwrap(),
        ] {
            let identity = base
                .clone()
                .with_diffuser(ParaxialBoxDiffuser::try_new([0.2; 2], 0.0, 1.0).unwrap())
                .unwrap();
            assert_eq!(
                base.transmitted_fraction().to_bits(),
                identity.transmitted_fraction().to_bits()
            );
            assert_eq!(
                base.maximum_luminance_factor([0.2; 2]).unwrap().to_bits(),
                identity
                    .maximum_luminance_factor([0.2; 2])
                    .unwrap()
                    .to_bits()
            );
            let frosted = base
                .clone()
                .with_diffuser(ParaxialBoxDiffuser::try_new([0.2, 0.3], 0.7, 0.8).unwrap())
                .unwrap();
            let bound = frosted.maximum_luminance_factor([0.2; 2]).unwrap();
            for x in -16..=16 {
                for y in -12..=12 {
                    let point = Point3::new(0.2, -0.2, 0.0);
                    let direction = Vec3::new(x as f64 * 0.2, y as f64 * 0.2, 1.0);
                    assert_eq!(
                        base.luminance_per_clear_flux_per_area(point, direction)
                            .unwrap()
                            .to_bits(),
                        identity
                            .luminance_per_clear_flux_per_area(point, direction)
                            .unwrap()
                            .to_bits()
                    );
                    let value = frosted
                        .luminance_per_clear_flux_per_area(point, direction)
                        .unwrap();
                    assert!(value >= 0.0 && value <= bound * (1.0 + 1.0e-14));
                }
            }
        }
    }

    #[test]
    fn g3_frost_tree_padding_and_unrepresentable_kernel_sizes_refuse() {
        let too_padded =
            GoboTransmission::try_new(2.0, 1.0, 1025, 513, vec![1.0; 1025 * 513]).unwrap();
        assert!(matches!(
            PositiveIntegral::try_new(&too_padded),
            Err(SourceError::DistributionTooLarge)
        ));
        let gobo = GoboTransmission::try_new(1.0, 1.0, 1, 1, vec![1.0]).unwrap();
        let law = ParaxialBoxDiffuser::try_new([1.0e-200; 2], 1.0, 1.0).unwrap();
        assert!(matches!(
            PreparedDiffuser::try_new(law, &gobo, 1.0e-200),
            Err(SourceError::NumericalRange)
        ));
    }
}
