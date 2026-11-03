//! Receiving-disk averages of the existing straight finite-window wake.
//!
//! Parallel equal-radius disks use exact circular overlap and its centroid.
//! Canted disks use a deterministic 32 x 64 equal-area polar rule, prepared
//! once per geometry. This integrates the declared field; it is not a free
//! wake, an error-certified quadrature or a partial-overlap blade-load model.

use fs_math::det;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct DiskWash {
    /// Receiving-area average of the source's axial wash factor.
    pub mean: f64,
    /// Average factor times source-relative transverse position / radius.
    /// Dot with a source-frame harmonic amplitude vector to transfer its mean.
    pub transverse: [f64; 3],
}

pub(crate) fn axial_factor(dx: f64, radius: f64) -> f64 {
    if dx.abs() >= 2.0 * radius {
        return 0.0;
    }
    let root = (dx * dx + radius * radius).sqrt();
    if dx >= 0.0 {
        1.0 - dx / root
    } else {
        1.0 + (-dx) / root
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let xy = a[0] * b[0] + a[1] * b[1];
    xy + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn between(
    source: [f64; 3],
    source_axis: [f64; 3],
    receiver: [f64; 3],
    receiver_axis: [f64; 3],
    radius: f64,
) -> Result<DiskWash, &'static str> {
    let offset = core::array::from_fn(|j| (receiver[j] - source[j]) / radius);
    if !radius.is_finite() || radius <= 0.0 || !offset.iter().all(|x| x.is_finite()) {
        return Err("invalid disk radius or relative geometry");
    }
    let normal_cross = cross(source_axis, receiver_axis);
    let cross2 = dot(normal_cross, normal_cross);
    let axial = dot(offset, source_axis);
    let radial: [f64; 3] = core::array::from_fn(|j| offset[j] - axial * source_axis[j]);
    let distance = dot(radial, radial).sqrt();
    if !axial.is_finite() || !distance.is_finite() {
        return Err("disk-relative geometry overflow");
    }
    // The receiving disk projects to a subset of a unit-radius circle; its
    // axial half-width is |source_axis x receiver_axis|. Reject disjoint
    // support before preparing any quadrature.
    let axial_half_width = cross2.sqrt();
    if distance >= 2.0 || axial.abs() >= 2.0 + axial_half_width {
        return Ok(DiskWash::default());
    }
    // Exact common-axis geometry, including opposite shaft conventions.
    if cross2 == 0.0 {
        if distance >= 2.0 || axial.abs() >= 2.0 {
            return Ok(DiskWash::default());
        }
        let fraction = if distance == 0.0 {
            1.0
        } else {
            let theta = det::acos(0.5 * distance);
            let segment = if theta < 0.001 {
                // theta - sin(theta)cos(theta), without cancellation for a
                // vanishing lens. Retain positive weak momentum transfer.
                let t2 = theta * theta;
                theta * t2 * (2.0 / 3.0 + t2 * (-2.0 / 15.0 + t2 * 4.0 / 315.0))
            } else {
                theta - 0.5 * distance * (1.0 - 0.25 * distance * distance).sqrt()
            };
            2.0 * segment / core::f64::consts::PI
        };
        let mean = axial_factor(axial, 1.0) * fraction;
        return Ok(DiskWash {
            mean,
            // The equal-circle intersection centroid is halfway between
            // projected centers; axial position is not a first-harmonic mode.
            transverse: radial.map(|x| mean * (0.5 * x)),
        });
    }
    // Geometry-defined bases rotate with both disks; no arbitrary world-axis
    // choice enters this approximate receiving-area integration.
    let radial_basis = normal_cross.map(|x| x / axial_half_width);
    let azimuth_basis = cross(receiver_axis, radial_basis);
    Ok(quadrature(
        offset,
        source_axis,
        radial_basis,
        azimuth_basis,
        32,
        64,
    ))
}

fn quadrature(
    offset: [f64; 3],
    source_axis: [f64; 3],
    radial_basis: [f64; 3],
    azimuth_basis: [f64; 3],
    radial_cells: usize,
    azimuth_cells: usize,
) -> DiskWash {
    let mut result = DiskWash::default();
    let weight = 1.0 / (radial_cells * azimuth_cells) as f64;
    for j in 0..azimuth_cells {
        let angle = core::f64::consts::TAU * (j as f64 + 0.5) / azimuth_cells as f64;
        let c = det::cos(angle);
        let s = det::sin(angle);
        let direction: [f64; 3] =
            core::array::from_fn(|k| c * radial_basis[k] + s * azimuth_basis[k]);
        for i in 0..radial_cells {
            let radius = ((i as f64 + 0.5) / radial_cells as f64).sqrt();
            let point = core::array::from_fn(|k| offset[k] + radius * direction[k]);
            let axial = dot(point, source_axis);
            let radial = core::array::from_fn(|k| point[k] - axial * source_axis[k]);
            if dot(radial, radial) >= 1.0 {
                continue;
            }
            let factor = weight * axial_factor(axial, 1.0);
            result.mean += factor;
            for k in 0..3 {
                result.transverse[k] += factor * radial[k];
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parallel_partial_overlap_and_harmonic_centroid_are_analytic() {
        let axis = [1.0, 0.0, 0.0];
        let result = between([0.0; 3], axis, [-1.0, 1.5, 0.0], axis, 1.0).unwrap();
        let lens = (2.0 * 0.75_f64.acos() - 0.75 * 1.75_f64.sqrt()) / core::f64::consts::PI;
        let expected = (1.0 + 1.0 / 2.0_f64.sqrt()) * lens;
        assert!((result.mean - expected).abs() < 1e-15);
        assert!((result.transverse[1] - 0.75 * expected).abs() < 1e-15);
        assert_eq!(result.transverse[0], 0.0);
        assert!(result.mean > 0.2); // the receiving center itself is outside the wake
        assert_eq!(
            between([0.0; 3], axis, [-1.0, 2.0, 0.0], axis, 1.0).unwrap(),
            DiskWash::default()
        );
        assert_eq!(
            between([0.0; 3], axis, [-2.0, 0.0, 0.0], axis, 1.0).unwrap(),
            DiskWash::default()
        );
        let touching = between([0.0; 3], axis, [0.0, 2.0 - 1e-10, 0.0], axis, 1.0).unwrap();
        assert!(touching.mean > 0.0);
        let scaled = between([3.0, -4.0, 2.0], axis, [1.0, -1.0, 2.0], axis, 2.0).unwrap();
        assert_eq!(scaled, result);
    }

    #[test]
    fn canted_area_rule_resolves_partial_overlap_and_refines() {
        let axis = [1.0, 0.0, 0.0];
        let length = 1.13_f64.sqrt();
        let other = [1.0 / length, 0.3 / length, 0.2 / length];
        let normal = cross(axis, other);
        let radial = normal.map(|x| x / dot(normal, normal).sqrt());
        let azimuth = cross(other, radial);
        let offset = [-0.5, 1.4, 0.2];
        let coarse = quadrature(offset, axis, radial, azimuth, 8, 16);
        let normal = between([0.0; 3], axis, offset, other, 1.0).unwrap();
        let fine = quadrature(offset, axis, radial, azimuth, 128, 256);
        assert!(normal.mean > 0.1);
        assert!((normal.mean - fine.mean).abs() < 0.003);
        assert!((normal.mean - fine.mean).abs() < (coarse.mean - fine.mean).abs());
        for k in 0..3 {
            assert!((normal.transverse[k] - fine.transverse[k]).abs() < 0.003);
        }
        // A centered tilted receiver is wholly inside the source cylinder;
        // opposite points have axial factors summing to exactly two.
        let centered = between([0.0; 3], axis, [0.0; 3], other, 1.0).unwrap();
        assert!((centered.mean - 1.0).abs() < 1e-14);
    }
}
