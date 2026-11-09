//! P1 traces on fluid triangles nested in one original planar shell facet.

use super::TensorShellInteriorError;
use fs_math::dd::Dd;

type Vector = [Dd; 3];

fn dot(a: Vector, b: Vector) -> Dd {
    a.into_iter()
        .zip(b)
        .fold(Dd::ZERO, |sum, (a, b)| sum + a * b)
}
fn norm(a: Vector) -> Dd {
    dot(a, a).sqrt()
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn invalid(what: &'static str) -> TensorShellInteriorError {
    TensorShellInteriorError::Invalid(what)
}

/// Rows address fluid vertices; columns address original shell vertices.
/// Geometry is prepared once, before the native pressure operator is assembled.
pub(super) fn contained_p1(
    fluid: [[f64; 3]; 3],
    shell: [[f64; 3]; 3],
    tolerance_m: f64,
) -> Result<[[f64; 3]; 3], TensorShellInteriorError> {
    if fluid.iter().chain(&shell).flatten().any(|x| !x.is_finite()) {
        return Err(invalid("finite interface coordinates required"));
    }
    let edge = |point: [f64; 3]| -> Vector {
        core::array::from_fn(|d| Dd::from_f64(point[d]) - Dd::from_f64(shell[0][d]))
    };
    let a = edge(shell[1]);
    let b = edge(shell[2]);
    let scale = a
        .iter()
        .chain(&b)
        .map(|x| x.hi.abs())
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid("unresolved original shell facet geometry"));
    }
    // Binary scaling avoids adding division residuals to coplanar geometry.
    let scale = f64::from_bits(scale.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE);
    let scale = Dd::from_f64(scale);
    let a = a.map(|x| x / scale);
    let b = b.map(|x| x / scale);
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let area = norm(cross);
    let bound = Dd::from_f64(tolerance_m) / scale;
    // Signed edge distance is barycentric coordinate times that vertex's
    // altitude. Compare it in the same scaled geometry as the plane distance.
    let altitudes = [
        area / norm(core::array::from_fn(|d| b[d] - a[d])),
        area / norm(b),
        area / norm(a),
    ];
    if !finite(area)
        || area.hi <= 0.0
        || !finite(bound)
        || altitudes.iter().any(|&h| !finite(h) || h.hi <= 0.0)
    {
        return Err(invalid("unresolved original shell facet geometry"));
    }
    // Drop the largest normal component instead of solving nearly dependent
    // three-dimensional normal equations. DD retains translated edge differences.
    let dropped = (1..3).fold(0, |best, d| {
        if cross[d].abs().hi > cross[best].abs().hi {
            d
        } else {
            best
        }
    });
    let (j, k) = match dropped {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let determinant = a[j] * b[k] - a[k] * b[j];
    let mut shapes = [[0.0; 3]; 3];
    for (row, point) in fluid.into_iter().enumerate() {
        if let Some(vertex) = shell.iter().position(|&source| source == point) {
            shapes[row][vertex] = 1.0;
            continue;
        }
        let r = edge(point).map(|x| x / scale);
        let plane_gap = dot(r, cross).abs() / area;
        if !finite(plane_gap) || bound.lt(plane_gap) {
            return Err(invalid(
                "fluid face leaves original shell plane beyond geometry allowance",
            ));
        }
        let l1 = (r[j] * b[k] - r[k] * b[j]) / determinant;
        let l2 = (a[j] * r[k] - a[k] * r[j]) / determinant;
        let coordinates = [Dd::ONE - l1 - l2, l1, l2];
        for i in 0..3 {
            if !finite(coordinates[i])
                || (coordinates[i].lt(Dd::ZERO) && bound.lt(-coordinates[i] * altitudes[i]))
            {
                return Err(invalid(
                    "fluid face crosses original shell edge beyond geometry allowance",
                ));
            }
            shapes[row][i] = coordinates[i].to_f64();
        }
    }
    Ok(shapes)
}

#[cfg(test)]
mod tests {
    use super::contained_p1;

    #[test]
    fn translated_tilted_midpoints_keep_zero_allowance_for_all_windings() {
        let shell = [
            [1e12, -2e12, 3e12],
            [1e12 + 2.0, -2e12 + 1.0, 3e12 + 1.0],
            [1e12 - 1.0, -2e12 + 2.0, 3e12 + 3.0],
        ];
        let fluid = [
            [1e12 + 1.0, -2e12 + 0.5, 3e12 + 0.5],
            [1e12 - 0.5, -2e12 + 1.0, 3e12 + 1.5],
            [1e12 + 0.5, -2e12 + 1.5, 3e12 + 2.0],
        ];
        let expected = [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]];
        for order in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let shapes = contained_p1(fluid, order.map(|i| shell[i]), 0.0).unwrap();
            for (actual, original) in shapes.into_iter().zip(expected) {
                for (value, reference) in actual.into_iter().zip(order.map(|i| original[i])) {
                    assert!((value - reference).abs() < 2e-15);
                }
            }
        }
    }
}
