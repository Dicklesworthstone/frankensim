use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_solid::linear3::tet10::{
    QuadraticTetProblem, TET10_EDGES, TetQuadrature, promote_p1_to_p2, tet10_point,
};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticMaterial, TetMaterialField,
};
use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};
use std::collections::BTreeMap;

pub const Q: TetQuadrature = TetQuadrature {
    stiffness_points_per_dir: 3,
    mass_points_per_dir: 4,
};

pub fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 58,
                kernel_id: 58,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

pub fn material(e: f64, nu: f64) -> TetElasticMaterial {
    TetElasticMaterial::try_new(1.0, e, nu, [58; 32]).unwrap()
}

#[derive(Clone)]
pub struct Mesh {
    pub nodes: Vec<[f64; 3]>,
    pub tets: Vec<[usize; 4]>,
}

/// Freudenthal block fixture; meshing algorithm/refinement stays in fs-mesh.
pub fn block(n: [usize; 3], extent: [f64; 3]) -> Mesh {
    let [nx, ny, nz] = n;
    let index = |i, j, k| (i * (ny + 1) + j) * (nz + 1) + k;
    let mut nodes = Vec::new();
    for i in 0..=nx {
        for j in 0..=ny {
            for k in 0..=nz {
                nodes.push([
                    extent[0] * i as f64 / nx as f64,
                    extent[1] * j as f64 / ny as f64,
                    extent[2] * k as f64 / nz as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let v = [
                    index(i, j, k),
                    index(i + 1, j, k),
                    index(i, j + 1, k),
                    index(i + 1, j + 1, k),
                    index(i, j, k + 1),
                    index(i + 1, j, k + 1),
                    index(i, j + 1, k + 1),
                    index(i + 1, j + 1, k + 1),
                ];
                tets.extend([
                    [v[0], v[1], v[3], v[7]],
                    [v[0], v[3], v[2], v[7]],
                    [v[0], v[2], v[6], v[7]],
                    [v[0], v[6], v[4], v[7]],
                    [v[0], v[4], v[5], v[7]],
                    [v[0], v[5], v[1], v[7]],
                ]);
            }
        }
    }
    Mesh { nodes, tets }
}

pub fn jittered_cube(n: usize) -> Mesh {
    let mut mesh = block([n, n, n], [1.0; 3]);
    for (i, p) in mesh.nodes.iter_mut().enumerate() {
        if p.iter().all(|v| *v > 0.0 && *v < 1.0) {
            for (k, x) in p.iter_mut().enumerate() {
                *x += 0.15 / n as f64 * (((i + 11 * k) * 17 % 19) as f64 / 19.0 - 0.5);
            }
        }
    }
    mesh
}

pub fn promote(mesh: &Mesh, cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(&mesh.nodes, &mesh.tets, None, cx).unwrap()
}

pub fn assembly(
    nodes: &[[f64; 3]],
    elements: &[[usize; 10]],
    m: &TetElasticMaterial,
    fixed: &[usize],
    cx: &Cx<'_>,
) -> TetElasticAssembly {
    QuadraticTetProblem {
        nodes_m: nodes,
        elements,
        materials: TetMaterialField::Uniform(m),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: Q,
    }
    .assemble(cx)
    .unwrap()
}

pub fn quadratic_form(k: &Csr, u: &[f64]) -> f64 {
    let mut ku = vec![0.0; u.len()];
    k.spmv(u, &mut ku);
    u.iter().zip(ku).map(|(a, b)| a * b).sum()
}

/// Sparse Dirichlet lifting against the FULL assembled operator. The exact
/// boundary values and right-hand side are independent inputs to this solve.
pub fn solve(k: &Csr, rhs: &[f64], prescribed: &[Option<f64>]) -> Vec<f64> {
    assert_eq!(k.nrows(), prescribed.len());
    let free: Vec<_> = (0..prescribed.len())
        .filter(|&i| prescribed[i].is_none())
        .collect();
    let mut reduced = vec![None; prescribed.len()];
    for (i, &full) in free.iter().enumerate() {
        reduced[full] = Some(i);
    }
    let mut matrix = Coo::new(free.len(), free.len());
    let mut force = vec![0.0; free.len()];
    for (i, &full) in free.iter().enumerate() {
        force[i] = rhs[full];
        let (cols, vals) = k.row(full);
        for (&col, &v) in cols.iter().zip(vals) {
            if let Some(j) = reduced[col] {
                matrix.push(i, j, v);
            } else {
                force[i] -= v * prescribed[col].unwrap();
            }
        }
    }
    let values = if free.is_empty() {
        Vec::new()
    } else {
        let matrix = matrix.assemble();
        let sym = SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd).unwrap();
        let factor = sym.factor(&matrix, &LdltOptions::default()).unwrap();
        let mut values = factor.solve(&force);
        let mut check = vec![0.0; values.len()];
        let norm = force.iter().map(|x| x * x).sum::<f64>().sqrt();
        let mut residual = f64::INFINITY;
        let mut initial_residual = 0.0;
        let mut refinements = 0;
        // Slender bending fixtures amplify factorization roundoff. Reuse the
        // same factors for at most three correction solves, rather than relax
        // the independently evaluated original SpMV residual gate below.
        for step in 0..=3 {
            matrix.spmv(&values, &mut check);
            residual = check
                .iter()
                .zip(&force)
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            if step == 0 {
                initial_residual = residual;
            }
            if residual <= 1e-9 * norm.max(1e-20) || step == 3 {
                break;
            }
            // Compensated products and sums keep cancellation in b-Au from
            // setting the correction accuracy. Final acceptance still uses
            // the unchanged ordinary SpMV, not this more accurate residual.
            let correction_rhs: Vec<_> = (0..values.len())
                .map(|row| {
                    let (cols, vals) = matrix.row(row);
                    let mut sum = force[row];
                    let mut low = 0.0;
                    for (&col, &a) in cols.iter().zip(vals) {
                        let product = a * values[col];
                        let term = -product;
                        let next = sum + term;
                        low += if sum.abs() >= term.abs() {
                            (sum - next) + term
                        } else {
                            (term - next) + sum
                        };
                        low -= a.mul_add(values[col], -product);
                        sum = next;
                    }
                    sum + low
                })
                .collect();
            let correction = factor.solve(&correction_rhs);
            for (value, delta) in values.iter_mut().zip(correction) {
                *value += delta;
            }
            refinements += 1;
        }
        if refinements > 0 {
            println!(
                "{{\"suite\":\"irp-s1\",\"stage\":\"static-iterative-refinement\",\"free_dofs\":{},\"steps\":{refinements},\"initial_relative_residual\":{},\"final_relative_residual\":{},\"tolerance\":1e-9}}",
                free.len(),
                initial_residual / norm.max(1e-20),
                residual / norm.max(1e-20)
            );
        }
        assert!(
            residual <= 1e-9 * norm.max(1e-20),
            "linear residual {residual:e}/{norm:e}"
        );
        values
    };
    (0..prescribed.len())
        .map(|i| prescribed[i].unwrap_or_else(|| values[reduced[i].unwrap()]))
        .collect()
}

pub fn volume_load(
    nodes: &[[f64; 3]],
    elements: &[[usize; 10]],
    f: impl Fn([f64; 3]) -> [f64; 3],
) -> Vec<f64> {
    let mut rhs = vec![0.0; 3 * nodes.len()];
    for (e, t) in elements.iter().enumerate() {
        for (q, (lambda, w)) in duffy_quadrature(5).into_iter().enumerate() {
            let p = tet10_point(t.map(|i| nodes[i]), lambda, 1e-12, e, q).unwrap();
            let force = f(p.position_m);
            for a in 0..10 {
                for c in 0..3 {
                    rhs[3 * t[a] + c] += w * p.determinant_m3 * p.shapes[a] * force[c];
                }
            }
        }
    }
    rhs
}

/// Constant traction on straight triangular traces, using the exact TET10
/// trace integrals: corner weights zero; each midside weight A/3.
pub fn surface_load(
    nodes: &[[f64; 3]],
    elements: &[[usize; 10]],
    on: impl Fn([f64; 3]) -> bool,
    traction: [f64; 3],
) -> Vec<f64> {
    let mut faces = BTreeMap::<[usize; 3], (usize, [usize; 6])>::new();
    for tet in elements {
        for local in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let [a, b, c] = local;
            let mut trace = [tet[a], tet[b], tet[c], 0, 0, 0];
            for (k, [u, v]) in [[a, b], [b, c], [a, c]].into_iter().enumerate() {
                let edge = TET10_EDGES
                    .iter()
                    .position(|pair| pair.contains(&u) && pair.contains(&v))
                    .unwrap();
                trace[3 + k] = tet[4 + edge];
            }
            let mut key = [tet[a], tet[b], tet[c]];
            key.sort_unstable();
            let entry = faces.entry(key).or_insert((0, trace));
            entry.0 += 1;
        }
    }
    let mut out = vec![0.0; 3 * nodes.len()];
    for (_, (count, face)) in faces {
        if count != 1 || !face[..3].iter().all(|&i| on(nodes[i])) {
            continue;
        }
        let a = nodes[face[0]];
        let b = nodes[face[1]];
        let c = nodes[face[2]];
        let u: [f64; 3] = core::array::from_fn(|k| b[k] - a[k]);
        let v: [f64; 3] = core::array::from_fn(|k| c[k] - a[k]);
        let cross = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let area = 0.5 * cross.iter().map(|x| x * x).sum::<f64>().sqrt();
        for &node in &face[3..] {
            for k in 0..3 {
                out[3 * node + k] += area / 3.0 * traction[k];
            }
        }
    }
    out
}

pub fn log(
    test: &str,
    tier: &str,
    fixture: &str,
    values: &[f64],
    dofs: &[usize],
    observed: &[f64],
    tolerance: f64,
    provenance: &str,
    elapsed_s: f64,
) {
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.3-W1a\",\"test\":\"{test}\",\"tier\":\"{tier}\",\"inputs\":{{\"fixture\":\"{fixture}\",\"element\":\"tet10\",\"quadrature_n\":[3,4]}},\"provenance\":\"{provenance}\",\"rungs\":{{\"axis\":\"h\",\"values\":{values:?},\"dofs\":{dofs:?}}},\"observed\":{observed:?},\"tolerance\":{tolerance},\"verdict\":\"pass\",\"elapsed_s\":{elapsed_s}}}"
    );
}
