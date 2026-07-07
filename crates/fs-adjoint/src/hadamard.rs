//! Hadamard shape gradients on boundary traces — the mathematically
//! clean boundary-integral form of the shape derivative. Two
//! classics, both VERIFIED against perturb-and-resolve finite
//! differences in the battery (the acceptance's sphere-volume /
//! compliance pattern, on the cube fixtures):
//!
//! - VOLUME: J(Ω) = |Ω| ⇒ dJ[V] = ∫_∂Ω V·n dA. On the discrete mesh
//!   this is EXACT for the discrete volume (both sides are polynomial
//!   in the perturbation parameter).
//! - COMPLIANCE (Dirichlet Poisson, J = ∫ f·u): the boundary form is
//!   dJ[V] = −∫_∂Ω (∂u/∂n)²·(V·n) dA (self-adjoint case λ = u). On P1
//!   discrete solutions this carries discretization error — the
//!   battery gates RELATIVE agreement with FD and reports the number
//!   instead of pretending exactness.

use fs_feec::ElementGeometry;
use fs_rep_mesh::TetComplex;

/// Boundary faces of a complex (faces incident to exactly one tet),
/// with their owning tet.
#[must_use]
pub fn boundary_faces(complex: &TetComplex) -> Vec<(usize, usize)> {
    let d2 = complex.d2();
    let mut face_use: Vec<Vec<usize>> = vec![Vec::new(); complex.faces.len()];
    for (t, row) in d2.rows.iter().enumerate() {
        for &(f, _) in row {
            face_use[f].push(t);
        }
    }
    face_use
        .iter()
        .enumerate()
        .filter(|(_, ts)| ts.len() == 1)
        .map(|(f, ts)| (f, ts[0]))
        .collect()
}

/// Outward unit normal and area of face `f` owned by tet `t`
/// (outward = away from the tet's fourth vertex).
fn face_normal_area(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    f: usize,
    t: usize,
) -> ([f64; 3], f64) {
    let tri = complex.faces[f];
    let (pa, pb, pc) = (
        positions[tri[0] as usize],
        positions[tri[1] as usize],
        positions[tri[2] as usize],
    );
    let e1 = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
    let e2 = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
    let mut n = [
        e1[1].mul_add(e2[2], -(e1[2] * e2[1])),
        e1[2].mul_add(e2[0], -(e1[0] * e2[2])),
        e1[0].mul_add(e2[1], -(e1[1] * e2[0])),
    ];
    let len = fs_math::det::sqrt(n[0].mul_add(n[0], n[1].mul_add(n[1], n[2] * n[2])));
    let area = 0.5 * len;
    for c in &mut n {
        *c /= len;
    }
    // Orient outward: opposite the vector to the tet's off-face vertex.
    let tet = complex.tets[t];
    let opp = tet
        .iter()
        .find(|v| !tri.contains(v))
        .expect("tet has a vertex off the face");
    let po = positions[*opp as usize];
    let centroid = [
        (pa[0] + pb[0] + pc[0]) / 3.0,
        (pa[1] + pb[1] + pc[1]) / 3.0,
        (pa[2] + pb[2] + pc[2]) / 3.0,
    ];
    let to_opp = [po[0] - centroid[0], po[1] - centroid[1], po[2] - centroid[2]];
    let dot = n[0].mul_add(to_opp[0], n[1].mul_add(to_opp[1], n[2] * to_opp[2]));
    if dot > 0.0 {
        for c in &mut n {
            *c = -*c;
        }
    }
    (n, area)
}

/// Hadamard VOLUME shape gradient: dJ[V] = ∫_∂Ω V·n dA with V given
/// nodally (P1 on the boundary; exact one-point-per-vertex
/// quadrature of the affine integrand: face average).
#[must_use]
pub fn volume_shape_gradient(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    velocity: &dyn Fn([f64; 3]) -> [f64; 3],
) -> f64 {
    let mut total = 0.0f64;
    for (f, t) in boundary_faces(complex) {
        let (n, area) = face_normal_area(complex, positions, f, t);
        let tri = complex.faces[f];
        // Affine V ⇒ vertex-average is the exact face integral.
        let mut vn = 0.0f64;
        for &v in &tri {
            let vel = velocity(positions[v as usize]);
            vn += n[0].mul_add(vel[0], n[1].mul_add(vel[1], n[2] * vel[2])) / 3.0;
        }
        total = vn.mul_add(area, total);
    }
    total
}

/// Hadamard COMPLIANCE shape gradient for the Dirichlet Poisson
/// problem: dJ[V] = −∫_∂Ω (∂u_h/∂n)²·(V·n) dA, with ∂u_h/∂n taken
/// from the owning tet's constant P1 gradient (u_h given at ALL
/// vertices, boundary values included).
#[must_use]
pub fn compliance_shape_gradient(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    u: &[f64],
    velocity: &dyn Fn([f64; 3]) -> [f64; 3],
) -> f64 {
    let mut total = 0.0f64;
    for (f, t) in boundary_faces(complex) {
        let (n, area) = face_normal_area(complex, positions, f, t);
        // Constant ∇u_h on tet t.
        let tet = complex.tets[t];
        let mut grad = [0.0f64; 3];
        for (a, &v) in tet.iter().enumerate() {
            for c in 0..3 {
                grad[c] = geo.grads[t][a][c].mul_add(u[v as usize], grad[c]);
            }
        }
        let dudn = n[0].mul_add(grad[0], n[1].mul_add(grad[1], n[2] * grad[2]));
        let tri = complex.faces[f];
        let mut vn = 0.0f64;
        for &v in &tri {
            let vel = velocity(positions[v as usize]);
            vn += n[0].mul_add(vel[0], n[1].mul_add(vel[1], n[2] * vel[2])) / 3.0;
        }
        total = (-dudn * dudn * vn).mul_add(area, total);
    }
    total
}
