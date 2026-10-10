//! The explicit quadratic-complete edge-cubic transverse plate field.
//! One Bernstein owner supplies panel inertia, reciprocal point/volume maps,
//! and analytic physical rotations; the DKT stiffness field is unchanged.

/// Degree-three Bernstein displacement coefficients (three vertices, six
/// edge controls, one interior control) as linear forms in `(w, wx, wy)`.
/// On each edge, the curve is the cubic Hermite interpolant of its endpoints
/// and tangential slopes. The symmetric interior rule reproduces all
/// quadratic polynomials exactly, while making no claim about DKT's absent
/// interior displacement field.
fn edge_cubic_displacement(x: &[f64; 3], y: &[f64; 3]) -> [[f64; 9]; 10] {
    let mut b = [[0.0; 9]; 10];
    for node in 0..3 {
        b[node][3 * node] = 1.0;
    }
    for (edge, (start, end)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
        let dx = (x[end] - x[start]) / 3.0;
        let dy = (y[end] - y[start]) / 3.0;
        b[3 + 2 * edge][3 * start] = 1.0;
        b[3 + 2 * edge][3 * start + 1] = dx;
        b[3 + 2 * edge][3 * start + 2] = dy;
        b[4 + 2 * edge][3 * end] = 1.0;
        b[4 + 2 * edge][3 * end + 1] = -dx;
        b[4 + 2 * edge][3 * end + 2] = -dy;
    }
    for dof in 0..9 {
        b[9][dof] = 0.25 * (3..9).map(|edge| b[edge][dof]).sum::<f64>()
            - (b[0][dof] + b[1][dof] + b[2][dof]) / 6.0;
    }
    b
}

const EDGE_CUBIC_EXP: [[i32; 3]; 10] = [
    [3, 0, 0],
    [0, 3, 0],
    [0, 0, 3],
    [2, 1, 0],
    [1, 2, 0],
    [0, 2, 1],
    [0, 1, 2],
    [1, 0, 2],
    [2, 0, 1],
    [1, 1, 1],
];
const EDGE_CUBIC_MULT: [f64; 10] = [1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 6.0];

/// Transverse displacement shape for the opt-in edge-cubic plate mass law.
/// The caller supplies a validated triangle and barycentric coordinates on
/// it. Entries use the local `(w, wx, wy)` order for each vertex. Apply the
/// same shape to point effort and motion to preserve virtual work.
#[must_use]
pub fn edge_cubic_transverse_shape(x: &[f64; 3], y: &[f64; 3], barycentric: [f64; 3]) -> [f64; 9] {
    let b = edge_cubic_displacement(x, y);
    let mut shape = [0.0; 9];
    for i in 0..10 {
        let value = EDGE_CUBIC_MULT[i]
            * (0..3)
                .map(|axis| barycentric[axis].powi(EDGE_CUBIC_EXP[i][axis]))
                .product::<f64>();
        for dof in 0..9 {
            shape[dof] += value * b[i][dof];
        }
    }
    shape
}

/// Physical x/y derivatives of the same edge-cubic transverse field. Rows
/// are `dw/dx` and `dw/dy`; columns retain the vertex `(w, wx, wy)` order.
/// The caller supplies a validated nondegenerate triangle and barycentric
/// point. These analytic derivatives include the interior Bernstein control;
/// interpolating the nodal slopes instead would describe a different field.
#[must_use]
pub fn edge_cubic_transverse_gradient_shape(
    x: &[f64; 3], y: &[f64; 3], barycentric: [f64; 3],
) -> [[f64; 9]; 2] {
    let area2 = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
    let gradient = [
        [(y[1] - y[2]) / area2, (y[2] - y[0]) / area2, (y[0] - y[1]) / area2],
        [(x[2] - x[1]) / area2, (x[0] - x[2]) / area2, (x[1] - x[0]) / area2],
    ];
    let b = edge_cubic_displacement(x, y);
    let mut shape = [[0.0; 9]; 2];
    for i in 0..10 {
        for axis in 0..3 {
            let exponent = EDGE_CUBIC_EXP[i][axis];
            if exponent == 0 { continue; }
            let derivative = EDGE_CUBIC_MULT[i] * f64::from(exponent)
                * (0..3).map(|j| barycentric[j].powi(
                    EDGE_CUBIC_EXP[i][j] - if j == axis { 1 } else { 0 },
                )).product::<f64>();
            for c in 0..2 {
                for dof in 0..9 {
                    shape[c][dof] += derivative * gradient[c][axis] * b[i][dof];
                }
            }
        }
    }
    shape
}

/// Area-average transverse displacement shape for the same cubic field.
#[must_use]
pub fn edge_cubic_transverse_mean_shape(x: &[f64; 3], y: &[f64; 3]) -> [f64; 9] {
    let b = edge_cubic_displacement(x, y);
    let mut mean = [0.0; 9];
    for row in &b {
        for dof in 0..9 {
            mean[dof] += row[dof] / 10.0;
        }
    }
    mean
}

/// Exact `rho*h*∫w²` matrix for the declared cubic Bernstein field.
pub(crate) fn edge_cubic_mass(x: &[f64; 3], y: &[f64; 3], rho_h_area: f64) -> [f64; 81] {
    const FACT: [f64; 7] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0];
    let b = edge_cubic_displacement(x, y);
    let mut mass = [0.0; 81];
    for i in 0..10 {
        for j in 0..10 {
            let gram = 2.0
                * rho_h_area
                * EDGE_CUBIC_MULT[i]
                * EDGE_CUBIC_MULT[j]
                * (0..3)
                    .map(|axis| FACT[(EDGE_CUBIC_EXP[i][axis] + EDGE_CUBIC_EXP[j][axis]) as usize])
                    .product::<f64>()
                / 40320.0;
            for r in 0..9 {
                for c in 0..9 {
                    mass[9 * r + c] += gram * b[i][r] * b[j][c];
                }
            }
        }
    }
    mass
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edge_cubic_mass_reproduces_quadratic_fields_and_physical_integrals() {
        let x = [0.0, 1.0, 0.0];
        let y = [0.0, 0.0, 1.0];
        let b = edge_cubic_displacement(&x, &y);
        let exponents = [
            [3, 0, 0],
            [0, 3, 0],
            [0, 0, 3],
            [2, 1, 0],
            [1, 2, 0],
            [0, 2, 1],
            [0, 1, 2],
            [1, 0, 2],
            [2, 0, 1],
            [1, 1, 1],
        ];
        let multiplicity = [1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 6.0];
        let check = |field: fn(f64, f64) -> (f64, f64, f64)| {
            let mut dofs = [0.0; 9];
            for node in 0..3 {
                let (w, wx, wy) = field(x[node], y[node]);
                dofs[3 * node..3 * node + 3].copy_from_slice(&[w, wx, wy]);
            }
            for bary in [[0.2_f64, 0.3, 0.5], [0.6, 0.1, 0.3], [0.1, 0.8, 0.1]] {
                let px = bary[1];
                let py = bary[2];
                let actual = exponents
                    .iter()
                    .enumerate()
                    .map(|(i, e)| {
                        let bernstein = multiplicity[i]
                            * (0..3)
                                .map(|axis| bary[axis].powi(e[axis] as i32))
                                .product::<f64>();
                        bernstein * b[i].iter().zip(dofs).map(|(v, u)| v * u).sum::<f64>()
                    })
                    .sum::<f64>();
                assert!((actual - field(px, py).0).abs() < 1e-14);
                let point_shape = edge_cubic_transverse_shape(&x, &y, bary);
                let sampled = point_shape
                    .iter()
                    .zip(dofs)
                    .map(|(v, u)| v * u)
                    .sum::<f64>();
                assert!((sampled - field(px, py).0).abs() < 1e-14);
            }
        };
        check(|_, _| (1.0, 0.0, 0.0));
        check(|x, _y| (x, 1.0, 0.0));
        check(|_, y| (y, 0.0, 1.0));
        check(|x, _| (x * x, 2.0 * x, 0.0));
        check(|x, y| (x * y, y, x));
        check(|_, y| (y * y, 0.0, 2.0 * y));

        let mean = edge_cubic_transverse_mean_shape(&x, &y);
        let three_point = [
            [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
        ];
        for local in 0..9 {
            let sampled = three_point
                .iter()
                .map(|&barycentric| edge_cubic_transverse_shape(&x, &y, barycentric)[local])
                .sum::<f64>()
                / 3.0;
            assert!((sampled - mean[local]).abs() < 1e-14);
        }
        let mean_value = |dofs: [f64; 9]| mean.iter().zip(dofs).map(|(v, u)| v * u).sum::<f64>();
        assert!((mean_value([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]) - 1.0).abs() < 1e-14);
        assert!(
            (mean_value([0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0]) - 1.0 / 3.0).abs() < 1e-14
        );
        assert!(
            (mean_value([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]) - 1.0 / 6.0).abs() < 1e-14
        );

        let mass = edge_cubic_mass(&x, &y, 0.5);
        let energy = |dofs: [f64; 9]| -> f64 {
            (0..9)
                .map(|r| {
                    (0..9)
                        .map(|c| dofs[r] * mass[9 * r + c] * dofs[c])
                        .sum::<f64>()
                })
                .sum()
        };
        assert!((energy([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]) - 0.5).abs() < 1e-14);
        assert!((energy([0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0]) - 1.0 / 12.0).abs() < 1e-14);
        assert!((energy([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]) - 1.0 / 30.0).abs() < 1e-14);
    }

    #[test]
    fn edge_cubic_gradients_reproduce_physical_quadratic_slopes_on_oblique_triangles() {
        let fields: [fn(f64, f64) -> [f64; 3]; 6] = [
            |_, _| [1.0, 0.0, 0.0], |x, _| [x, 1.0, 0.0], |_, y| [y, 0.0, 1.0],
            |x, _| [x*x, 2.0*x, 0.0], |x, y| [x*y, y, x], |_, y| [y*y, 0.0, 2.0*y],
        ];
        for order in [[0, 1, 2], [2, 1, 0]] {
            let x = order.map(|i| [-0.3, 1.7, 0.4][i]);
            let y = order.map(|i| [-0.2, 0.6, 1.3][i]);
            for field in fields {
                let dofs: [f64; 9] = std::array::from_fn(|i| field(x[i/3], y[i/3])[i%3]);
                for bary in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
                    [0.2, 0.8, 0.0], [0.2, 0.3, 0.5]] {
                    let px = (0..3).map(|i| bary[i]*x[i]).sum();
                    let py = (0..3).map(|i| bary[i]*y[i]).sum();
                    let expected = field(px, py);
                    let shape = edge_cubic_transverse_gradient_shape(&x, &y, bary);
                    for c in 0..2 {
                        let actual: f64 = shape[c].iter().zip(dofs).map(|(a,b)| a*b).sum();
                        assert!((actual - expected[c+1]).abs() < 2e-14);
                    }
                }
            }
        }
    }

    #[test]
    fn edge_cubic_gradients_differentiate_the_actual_cubic_interior_and_edges() {
        // Nodal values/slopes of x^3 reconstruct x^3 - (1-x-y)*x*y,
        // not x^3: the existing symmetric interior rule is quadratic-complete.
        let x = [0.0, 1.0, 0.0]; let y = [0.0, 0.0, 1.0];
        let dofs = [0.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0];
        for [px, py] in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
            [0.37, 0.0], [0.23, 0.41], [0.6, 0.4]] {
            let bary = [1.0-px-py, px, py];
            let w: f64 = edge_cubic_transverse_shape(&x, &y, bary).iter()
                .zip(dofs).map(|(a,b)| a*b).sum();
            assert!((w - (px*px*px - (1.0-px-py)*px*py)).abs() < 1e-14);
            let gradient = edge_cubic_transverse_gradient_shape(&x, &y, bary);
            let expected = [3.0*px*px + 2.0*px*py + py*py - py,
                px*px + 2.0*px*py - px];
            for c in 0..2 {
                let actual: f64 = gradient[c].iter().zip(dofs).map(|(a,b)| a*b).sum();
                assert!((actual - expected[c]).abs() < 1e-14);
            }
        }
    }
}
