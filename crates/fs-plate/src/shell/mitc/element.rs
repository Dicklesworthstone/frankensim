use super::{PlateError, ShearShellSection, ShellMass, ShellMesh, cross, dot};
use crate::shell::FacetGeometry;
use fs_math::det;

const N: usize = 20;
type Matrix = [[f64; N]; N];
pub(super) type Condensed = [[f64; 18]; 18];

fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = det::sqrt(dot(v, v));
    v.map(|x| x / n)
}

struct Kinematics {
    normals: [[f64; 3]; 3],
    average: [f64; 3],
    bubble_axes: [[f64; 3]; 2],
    arms: [[f64; 3]; 2],
    director_gradient: [[f64; 3]; 2],
}
impl Kinematics {
    fn new(mesh: &ShellMesh, tri: &[usize; 3], g: &FacetGeometry, normals: [[f64; 3]; 3]) -> Self {
        let average =
            core::array::from_fn(|k| (normals[0][k] + normals[1][k] + normals[2][k]) / 3.0);
        let normal = unit(average);
        let ex = unit(core::array::from_fn(|k| {
            g.frame[0][k] - dot(g.frame[0], normal) * normal[k]
        }));
        let ey = cross(normal, ex);
        let arms = [
            core::array::from_fn(|k| mesh.nodes[tri[1]][k] - mesh.nodes[tri[0]][k]),
            core::array::from_fn(|k| mesh.nodes[tri[2]][k] - mesh.nodes[tri[0]][k]),
        ];
        let director_gradient = core::array::from_fn(|a| {
            core::array::from_fn(|k| (0..3).map(|i| g.gradient[i][a] * normals[i][k]).sum())
        });
        Self {
            normals,
            average,
            bubble_axes: [ex, ey],
            arms,
            director_gradient,
        }
    }

    fn basis(
        &self,
        g: &FacetGeometry,
        r: f64,
        s: f64,
        col: usize,
    ) -> ([f64; 3], [f64; 3], [f64; 2]) {
        let h = [1.0 - r - s, r, s];
        let bubble = 27.0 * r * s * h[0];
        let db = [
            27.0 * s * (1.0 - 2.0 * r - s),
            27.0 * r * (1.0 - r - 2.0 * s),
        ];
        let db_xy = core::array::from_fn::<_, 2, _>(|k| {
            db[0] * g.gradient[1][k] + db[1] * g.gradient[2][k]
        });
        let d: [f64; 3] = core::array::from_fn(|k| (0..3).map(|i| h[i] * self.normals[i][k]).sum());
        let mut ux = [0.0; 3];
        let mut uy = [0.0; 3];
        let mut ur = [0.0; 3];
        let mut us = [0.0; 3];
        let mut v = [0.0; 3];
        let mut vx = [0.0; 3];
        let mut vy = [0.0; 3];
        if col < 18 {
            let node = col / 6;
            let component = col % 6;
            if component < 3 {
                ux[component] = g.gradient[node][0];
                uy[component] = g.gradient[node][1];
                ur[component] = [-1.0, 1.0, 0.0][node];
                us[component] = [-1.0, 0.0, 1.0][node];
            } else {
                let mut rotation = [0.0; 3];
                rotation[component - 3] = 1.0;
                let motion = cross(rotation, self.normals[node]);
                let f = h[node] - bubble / 3.0;
                let fx = g.gradient[node][0] - db_xy[0] / 3.0;
                let fy = g.gradient[node][1] - db_xy[1] / 3.0;
                v = motion.map(|x| f * x);
                vx = motion.map(|x| fx * x);
                vy = motion.map(|x| fy * x);
            }
        } else {
            let motion = cross(self.bubble_axes[col - 18], self.average);
            v = motion.map(|x| bubble * x);
            vx = motion.map(|x| db_xy[0] * x);
            vy = motion.map(|x| db_xy[1] * x);
        }
        // Coefficients of the linearized covariant strain at z=0 and its
        // first thickness moment, Eq.3 with Eq.6 geometry. The d_,a · u_,b
        // terms are essential for rigid rotation on varying directors.
        let membrane = [
            dot(g.frame[0], ux),
            dot(g.frame[1], uy),
            dot(g.frame[0], uy) + dot(g.frame[1], ux),
        ];
        let bending = [
            dot(g.frame[0], vx) + dot(self.director_gradient[0], ux),
            dot(g.frame[1], vy) + dot(self.director_gradient[1], uy),
            dot(g.frame[0], vy)
                + dot(g.frame[1], vx)
                + dot(self.director_gradient[0], uy)
                + dot(self.director_gradient[1], ux),
        ];
        let shear = [
            dot(self.arms[0], v) + dot(d, ur),
            dot(self.arms[1], v) + dot(d, us),
        ];
        (membrane, bending, shear)
    }

    fn shear_ties(&self, g: &FacetGeometry) -> ([[f64; N]; 2], [f64; N]) {
        let mut constant = [[0.0; N]; 2];
        let mut twisting = [0.0; N];
        let delta = 1.0 / 10_000.0;
        let third = 1.0 / 3.0;
        for c in 0..N {
            let a = self.basis(g, 1.0 / 6.0, 2.0 / 3.0, c).2;
            let b = self.basis(g, 2.0 / 3.0, 1.0 / 6.0, c).2;
            let cc = self.basis(g, 1.0 / 6.0, 1.0 / 6.0, c).2;
            let d = self.basis(g, third + delta, third - 2.0 * delta, c).2;
            let e = self.basis(g, third - 2.0 * delta, third + delta, c).2;
            let f = self.basis(g, third + delta, third + delta, c).2;
            constant[0][c] = (2.0 * b[0] - b[1] + cc[0] + cc[1]) / 3.0;
            constant[1][c] = (2.0 * a[1] - a[0] + cc[0] + cc[1]) / 3.0;
            twisting[c] = (f[0] - d[0]) - (f[1] - e[1]);
        }
        (constant, twisting)
    }

    fn inertia_basis(&self, r: f64, s: f64, col: usize) -> ([f64; 3], [f64; 3]) {
        let h = [1.0 - r - s, r, s];
        let bubble = 27.0 * r * s * h[0];
        let mut translation = [0.0; 3];
        let mut rotation = [0.0; 3];
        if col < 18 {
            let node = col / 6;
            let c = col % 6;
            if c < 3 {
                translation[c] = h[node];
            } else {
                rotation[c - 3] = h[node];
                // Enrich the two director rotations only; the drilling
                // coordinate remains a nodal gauge. Constant angular fields
                // are retained when the internal rotation is their projection.
                for axis in self.bubble_axes {
                    for k in 0..3 {
                        rotation[k] -= bubble / 3.0 * axis[k] * axis[c - 3];
                    }
                }
            }
        } else {
            rotation = self.bubble_axes[col - 18].map(|v| bubble * v);
        }
        (translation, rotation)
    }
}

fn stiffness_rule() -> Vec<([f64; 3], f64)> {
    // Positive degree-five symmetric triangle rule; weights sum to one.
    let mut q = vec![([1.0 / 3.0; 3], 0.225)];
    for (a, b, w) in [
        (0.059715871789770, 0.470142064105115, 0.132394152788506),
        (0.797426985353087, 0.101286507323456, 0.125939180544827),
    ] {
        q.extend([([a, b, b], w), ([b, a, b], w), ([b, b, a], w)]);
    }
    q
}

fn mass_rule() -> Vec<([f64; 3], f64)> {
    // Positive collapsed five-by-five Gauss rule, exact beyond the degree-six
    // enriched rotation products. Reference weights normalized to unit area.
    let nodes = [
        -0.9061798459386640,
        -0.5384693101056831,
        0.0,
        0.5384693101056831,
        0.9061798459386640,
    ];
    let weights = [
        0.2369268850561891,
        0.4786286704993665,
        0.5688888888888889,
        0.4786286704993665,
        0.2369268850561891,
    ];
    let mut q = Vec::with_capacity(25);
    for i in 0..5 {
        let r = 0.5 * (1.0 + nodes[i]);
        for j in 0..5 {
            let s = 0.5 * (1.0 + nodes[j]) * (1.0 - r);
            q.push((
                [1.0 - r - s, r, s],
                0.5 * weights[i] * weights[j] * (1.0 - r),
            ));
        }
    }
    q
}

pub(super) fn matrices(
    mesh: &ShellMesh,
    tri: &[usize; 3],
    g: &FacetGeometry,
    normals: [[f64; 3]; 3],
    section: &ShearShellSection,
    mass: ShellMass,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Condensed, Condensed), PlateError> {
    let kin = Kinematics::new(mesh, tri, g, normals);
    let (constant, twisting) = kin.shear_ties(g);
    let (a, d, ks) = section.matrices();
    let mut k: Matrix = [[0.0; N]; N];
    let mut m: Matrix = [[0.0; N]; N];
    for (lambda, weight) in stiffness_rule() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let r = lambda[1];
        let s = lambda[2];
        let mut bm = [[0.0; N]; 3];
        let mut bb = [[0.0; N]; 3];
        let mut bs = [[0.0; N]; 2];
        for c in 0..N {
            let (membrane, bending, _) = kin.basis(g, r, s, c);
            for i in 0..3 {
                bm[i][c] = membrane[i];
                bb[i][c] = bending[i];
            }
            let cov = [
                constant[0][c] + twisting[c] * (s - 1.0 / 3.0),
                constant[1][c] + twisting[c] * (1.0 / 3.0 - r),
            ];
            for i in 0..2 {
                bs[i][c] = g.gradient[1][i] * cov[0] + g.gradient[2][i] * cov[1];
            }
        }
        let mut am = [[0.0; N]; 3];
        let mut db = [[0.0; N]; 3];
        let mut ss = [[0.0; N]; 2];
        for p in 0..3 {
            for c in 0..N {
                for q in 0..3 {
                    am[p][c] += a[3 * p + q] * bm[q][c];
                    db[p][c] += d[3 * p + q] * bb[q][c];
                }
            }
        }
        for p in 0..2 {
            for c in 0..N {
                for q in 0..2 {
                    ss[p][c] += ks[p][q] * bs[q][c];
                }
            }
        }
        for i in 0..N {
            for j in i..N {
                let v = (0..3)
                    .map(|p| bm[p][i] * am[p][j] + bb[p][i] * db[p][j])
                    .sum::<f64>()
                    + (0..2).map(|p| bs[p][i] * ss[p][j]).sum::<f64>();
                k[i][j] += g.area_m2 * weight * v;
            }
        }
    }
    if mass == ShellMass::Consistent {
        for (lambda, weight) in mass_rule() {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let basis: Vec<_> = (0..N)
                .map(|c| kin.inertia_basis(lambda[1], lambda[2], c))
                .collect();
            for i in 0..N {
                for j in i..N {
                    m[i][j] += weight
                        * g.area_m2
                        * (section.areal_mass_kg_m2 * dot(basis[i].0, basis[j].0)
                            + section.rotary_inertia_kg * dot(basis[i].1, basis[j].1));
                }
            }
        }
    } else {
        for i in 0..18 {
            m[i][i] = g.area_m2 / 3.0
                * if i % 6 < 3 {
                    section.areal_mass_kg_m2
                } else {
                    section.rotary_inertia_kg
                };
        }
    }
    for i in 0..N {
        for j in i..N {
            k[j][i] = k[i][j];
            m[j][i] = m[i][j];
        }
    }
    let c00 = k[18][18];
    let c01 = k[18][19];
    let c11 = k[19][19];
    let scale = det::sqrt(c00) * det::sqrt(c11);
    let correlation = c01 / scale;
    if !scale.is_finite()
        || c00 <= 0.0
        || c11 <= 0.0
        || !correlation.is_finite()
        || correlation.abs() >= 1.0
    {
        return Err(PlateError::BadSection {
            what: "MITC bubble block is not positive definite",
        });
    }
    let determinant = c00 * c11 - c01 * c01;
    let mut transform = [[0.0; 18]; N];
    for i in 0..18 {
        transform[i][i] = 1.0;
    }
    for i in 0..18 {
        transform[18][i] = (-c11 * k[18][i] + c01 * k[19][i]) / determinant;
        transform[19][i] = (c01 * k[18][i] - c00 * k[19][i]) / determinant;
    }
    let mut kc = [[0.0; 18]; 18];
    let mut mc = [[0.0; 18]; 18];
    let mut mt = [[0.0; 18]; N];
    for a in 0..N {
        for j in 0..18 {
            for b in 0..N {
                mt[a][j] += m[a][b] * transform[b][j];
            }
        }
    }
    for i in 0..18 {
        for j in i..18 {
            kc[i][j] = k[i][j] + k[i][18] * transform[18][j] + k[i][19] * transform[19][j];
            for a in 0..N {
                mc[i][j] += transform[a][i] * mt[a][j];
            }
            if !kc[i][j].is_finite() || !mc[i][j].is_finite() {
                return Err(PlateError::BadSection {
                    what: "nonfinite condensed MITC pencil",
                });
            }
            kc[j][i] = kc[i][j];
            mc[j][i] = mc[i][j];
        }
    }
    Ok((kc, mc))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_transverse_shear_patch() {
        let mesh = ShellMesh::new(
            vec![[0.0; 3], [1.3, 0.2, 0.0], [0.1, 0.8, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let g = mesh.facet(0).unwrap();
        let kin = Kinematics::new(&mesh, &mesh.tris[0], &g, [[0.0, 0.0, 1.0]; 3]);
        let (constant, twisting) = kin.shear_ties(&g);
        let gamma = [0.17, -0.08];
        let mut u = [0.0; 20];
        for (i, p) in mesh.nodes.iter().enumerate() {
            u[6 * i + 2] = gamma[0] * p[0] + gamma[1] * p[1];
        }
        for (lambda, _) in stiffness_rule() {
            let cov = core::array::from_fn::<_, 2, _>(|a| {
                (0..20)
                    .map(|c| {
                        let linear = if a == 0 {
                            lambda[2] - 1.0 / 3.0
                        } else {
                            1.0 / 3.0 - lambda[1]
                        };
                        (constant[a][c] + twisting[c] * linear) * u[c]
                    })
                    .sum::<f64>()
            });
            for a in 0..2 {
                let observed = g.gradient[1][a] * cov[0] + g.gradient[2][a] * cov[1];
                let expected = gamma[0] * g.frame[a][0] + gamma[1] * g.frame[a][1];
                assert!((observed - expected).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn quadrature_integrates_enriched_rotation_products() {
        fn factorial(n: usize) -> f64 {
            (1..=n).map(|x| x as f64).product()
        }
        for (rule, degree) in [(stiffness_rule(), 5), (mass_rule(), 6)] {
            for a in 0..=degree {
                for b in 0..=degree - a {
                    let observed: f64 = rule
                        .iter()
                        .map(|(l, w)| w * l[1].powi(a as i32) * l[2].powi(b as i32))
                        .sum();
                    let exact = 2.0 * factorial(a) * factorial(b) / factorial(a + b + 2);
                    assert!((observed - exact).abs() < 2e-14);
                }
            }
        }
    }
}
