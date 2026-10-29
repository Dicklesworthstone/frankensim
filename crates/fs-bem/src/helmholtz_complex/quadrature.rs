//! Fixed geometric quadrature for analytic complex continuation.
//!
//! These are the original triangle-CBIE fan/Duffy and tensor-product node
//! locations and weights, independent of frequency. Keeping values and their
//! derivatives in the same traversal is essential for a holomorphic pencil.
//! The existing real-frequency assembly remains the bitwise regression oracle.
use super::{GL8_W, GL8_X, affine3, cross, dot, norm, sub};

pub(super) fn triangle_remainder_rule<const ORDER: usize>(
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    ny: [f64; 3],
    rule: (&[f64; ORDER], &[f64; ORDER]),
    mut visit: impl FnMut(f64, f64),
) {
    let (nodes, weights) = rule;
    let [a, b, c] = triangle;
    let ab = sub(b, a);
    let ac = sub(c, a);
    let double_area = norm(cross(ab, ac));
    let height = dot(sub(a, x), ny);
    let x_proj = [x[0] + height * ny[0], x[1] + height * ny[1], x[2] + height * ny[2]];
    let sub_areas: [f64; 3] = core::array::from_fn(|edge| {
        let e0 = sub(triangle[edge], x_proj);
        let e1 = sub(triangle[(edge + 1) % 3], x_proj);
        dot(cross(e0, e1), ny)
    });
    if sub_areas.iter().all(|&sa| sa > 0.0) {
        for (edge, &sub_double_area) in sub_areas.iter().enumerate() {
            let e0 = sub(triangle[edge], x_proj);
            let e1 = sub(triangle[(edge + 1) % 3], x_proj);
            for (&radial, &wr) in nodes.iter().zip(weights) {
                for (&tangent, &wt) in nodes.iter().zip(weights) {
                    let ray = [
                        (1.0 - tangent).mul_add(e0[0], tangent * e1[0]),
                        (1.0 - tangent).mul_add(e0[1], tangent * e1[1]),
                        (1.0 - tangent).mul_add(e0[2], tangent * e1[2]),
                    ];
                    let y = [
                        radial.mul_add(ray[0], x_proj[0]),
                        radial.mul_add(ray[1], x_proj[1]),
                        radial.mul_add(ray[2], x_proj[2]),
                    ];
                    visit(norm(sub(x, y)), wr * wt * radial * sub_double_area);
                }
            }
        }
        return;
    }
    for (&u, &wu) in nodes.iter().zip(weights) {
        for (&v, &wv) in nodes.iter().zip(weights) {
            let y = affine3(a, ab, u, ac, (1.0 - u) * v);
            visit(norm(sub(x, y)), wu * wv * (1.0 - u) * double_area);
        }
    }
}

pub(super) fn triangle_self_rule(
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    mut visit: impl FnMut(f64, f64),
) {
    for edge in 0..3 {
        let e0 = sub(triangle[edge], x);
        let e1 = sub(triangle[(edge + 1) % 3], x);
        let double_area = norm(cross(e0, e1));
        for (&radial, &wr) in GL8_X.iter().zip(&GL8_W) {
            for (&tangent, &wt) in GL8_X.iter().zip(&GL8_W) {
                let ray = [
                    (1.0 - tangent).mul_add(e0[0], tangent * e1[0]),
                    (1.0 - tangent).mul_add(e0[1], tangent * e1[1]),
                    (1.0 - tangent).mul_add(e0[2], tangent * e1[2]),
                ];
                visit(radial * norm(ray), wr * wt * radial * double_area);
            }
        }
    }
}
