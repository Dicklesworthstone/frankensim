//! The TET10 integration/assembly path, shared with quadratic hexahedra.
//! Local upper triangles, Mandel arithmetic and COO insertion retain the
//! original TET10 operation order. Geometry and reference rules are inputs.
use super::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetMaterialField, free_dof_map,
};
use fs_exec::Cx;
use fs_sparse::Coo;

pub(super) struct Point<const N: usize> {
    pub shapes: [f64; N],
    pub gradients: [[f64; 3]; N],
    pub determinant: f64,
    pub scaled_jacobian: f64,
}

pub(super) fn assemble<const N: usize, const D: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    fixed_dofs: &[usize],
    budget: TetAssemblyBudget,
    stiffness_rule: &[(R, f64)],
    mass_rule: &[(R, f64)],
    evaluate: impl Fn([[f64; 3]; N], R, usize, usize) -> Result<Point<N>, TetElasticError>,
    cx: &Cx<'_>,
) -> Result<TetElasticAssembly, TetElasticError> {
    debug_assert_eq!(D, 3 * N);
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    let (free_dofs, reduced) = free_dof_map(nodes_m.len(), fixed_dofs, budget.maximum_free_dofs)?;
    let n = free_dofs.len();
    let mut stiffness = Coo::new(n, n);
    let mut mass = Coo::new(n, n);
    let mut total_mass_kg = 0.0;
    let mut minimum_scaled_jacobian = f64::INFINITY;
    let mut element_volumes_m3 = Vec::with_capacity(elements.len());
    for (element, connectivity) in elements.iter().enumerate() {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let nodes = connectivity.map(|i| nodes_m[i]);
        let material = materials.at(element);
        let mut ke = [[0.0; D]; D];
        let mut me = [[0.0; D]; D];
        let mut volume = 0.0;
        for (point, &(reference, weight)) in stiffness_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let p = evaluate(nodes, reference, element, point)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(p.scaled_jacobian);
            let b = mandel_b::<N, D>(p.gradients);
            let mut cb = [[0.0; D]; 6];
            for i in 0..6 {
                for j in 0..D {
                    for k in 0..6 {
                        cb[i][j] += material.stiffness_mandel_pa[i][k] * b[k][j];
                    }
                }
            }
            let dv = weight * p.determinant;
            for i in 0..D {
                for j in i..D {
                    let entry: f64 = (0..6).map(|k| b[k][i] * cb[k][j]).sum();
                    ke[i][j] += dv * entry;
                }
            }
        }
        for (point, &(reference, weight)) in mass_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let p = evaluate(nodes, reference, element, stiffness_rule.len() + point)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(p.scaled_jacobian);
            let dv = weight * p.determinant;
            volume += dv;
            for a in 0..N {
                for b in a..N {
                    let v = material.density_kg_m3 * dv * p.shapes[a] * p.shapes[b];
                    for c in 0..3 {
                        me[3 * a + c][3 * b + c] += v;
                    }
                }
            }
        }
        for i in 0..D {
            for j in i..D {
                if !ke[i][j].is_finite() || !me[i][j].is_finite() {
                    return Err(TetElasticError::InvalidMaterial {
                        what: "nonfinite quadratic element matrix",
                    });
                }
                ke[j][i] = ke[i][j];
                me[j][i] = me[i][j];
            }
        }
        for i in 0..D {
            let Some(row) = reduced[3 * connectivity[i / 3] + i % 3] else {
                continue;
            };
            for j in 0..D {
                let Some(col) = reduced[3 * connectivity[j / 3] + j % 3] else {
                    continue;
                };
                if ke[i][j] != 0.0 {
                    stiffness.push(row, col, ke[i][j]);
                }
                if me[i][j] != 0.0 {
                    mass.push(row, col, me[i][j]);
                }
            }
        }
        element_volumes_m3.push(volume);
        total_mass_kg += material.density_kg_m3 * volume;
    }
    if !(total_mass_kg.is_finite() && total_mass_kg > 0.0) {
        return Err(TetElasticError::InvalidMaterial {
            what: "nonfinite quadratic total mass",
        });
    }
    let stiffness = stiffness.assemble();
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    let mass = mass.assemble();
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    if (0..n).any(|row| {
        stiffness
            .row(row)
            .1
            .iter()
            .chain(mass.row(row).1)
            .any(|v| !v.is_finite())
    }) {
        return Err(TetElasticError::InvalidMaterial {
            what: "nonfinite assembled quadratic matrix",
        });
    }
    Ok(TetElasticAssembly {
        stiffness,
        mass,
        free_dofs,
        total_mass_kg,
        minimum_scaled_jacobian,
        element_volumes_m3,
    })
}
fn mandel_b<const N: usize, const D: usize>(gradients: [[f64; 3]; N]) -> [[f64; D]; 6] {
    let mut b = [[0.0; D]; 6];
    for (a, g) in gradients.iter().enumerate() {
        let c = 3 * a;
        let s = core::f64::consts::FRAC_1_SQRT_2;
        b[0][c] = g[0];
        b[1][c + 1] = g[1];
        b[2][c + 2] = g[2];
        b[3][c] = s * g[1];
        b[3][c + 1] = s * g[0];
        b[4][c + 1] = s * g[2];
        b[4][c + 2] = s * g[1];
        b[5][c] = s * g[2];
        b[5][c + 2] = s * g[0];
    }
    b
}

// Public P1 support-motion methods reuse the private linear3 assembly/solve owners.
mod prescribed;

// Full-row P1 dynamic response reuses the same canonical element kernels.
mod response;
