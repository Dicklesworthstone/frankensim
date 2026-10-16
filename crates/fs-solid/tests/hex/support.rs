pub use super::tet_support::{material, quadratic_form, with_cx};
use fs_exec::Cx;
pub use fs_mesh::structured_hex::{
    HexElements, HexMeshBudget, StructuredHexMesh, extruded_box_hex,
};
pub use fs_solid::linear3::hex::*;
pub use fs_solid::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial, TetMaterialField,
};

pub const FAMILIES: [HexFamily; 2] = [HexFamily::Serendipity20, HexFamily::Lagrange27];
pub fn family_name(f: HexFamily) -> &'static str {
    match f {
        HexFamily::Serendipity20 => "hex20",
        HexFamily::Lagrange27 => "hex27",
    }
}
pub fn mesh(n: [usize; 3], size: [f64; 3], f: HexFamily, cx: &Cx<'_>) -> StructuredHexMesh {
    extruded_box_hex(
        n,
        size,
        f,
        HexMeshBudget {
            maximum_nodes: 1_000_000,
            maximum_elements: 1_000_000,
        },
        cx,
    )
    .unwrap()
}
pub fn connectivity(elements: &HexElements) -> HexConnectivity<'_> {
    match elements {
        HexElements::Hex20(e) => HexConnectivity::Hex20(e),
        HexElements::Hex27(e) => HexConnectivity::Hex27(e),
    }
}
pub fn problem<'a>(
    m: &'a StructuredHexMesh,
    law: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> HexElasticProblem<'a> {
    HexElasticProblem {
        nodes_m: &m.nodes_m,
        family: m.elements.family(),
        elements: connectivity(&m.elements),
        materials: TetMaterialField::Uniform(law),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: HexQuadrature { points_per_dir: 3 },
    }
}
pub fn assembly(
    m: &StructuredHexMesh,
    law: &TetElasticMaterial,
    fixed: &[usize],
    cx: &Cx<'_>,
) -> TetElasticAssembly {
    problem(m, law, fixed).assemble(cx).unwrap()
}
pub fn modes(a: &TetElasticAssembly, window: (f64, f64)) -> Vec<fs_modal::ModePair> {
    fs_modal::slice_window(
        &a.stiffness,
        &a.mass,
        window,
        &fs_modal::SliceOptions {
            max_lanczos: 512,
            max_restarts: 32,
            ritz_tol: 1e-13,
            mass_diagonal_equilibration: true,
            ..fs_modal::SliceOptions::default()
        },
    )
    .unwrap()
    .modes
}
pub fn bending(p: [f64; 3], nu: f64) -> [f64; 3] {
    let [x, y, z] = p;
    [-x * y, (x * x + nu * (y * y - z * z)) / 2.0, nu * y * z]
}
pub fn bending_gradient(p: [f64; 3], nu: f64) -> [[f64; 3]; 3] {
    let [x, y, z] = p;
    [[-y, -x, 0.0], [x, nu * y, -nu * z], [0.0, nu * z, nu * y]]
}
pub fn distortion(p: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = p;
    [x + 0.12 * y * z, y + 0.09 * x * z, z + 0.08 * x * y]
}
pub fn point_data(
    mesh: &StructuredHexMesh,
    r: [f64; 3],
) -> ([f64; 3], Vec<f64>, Vec<[f64; 3]>, Vec<usize>) {
    match &mesh.elements {
        HexElements::Hex20(e) => {
            let p = hex20_point(e[0].map(|i| mesh.nodes_m[i]), r, 1e-12, 0, 0).unwrap();
            (
                p.position_m,
                p.shapes.to_vec(),
                p.gradients.to_vec(),
                e[0].to_vec(),
            )
        }
        HexElements::Hex27(e) => {
            let p = hex27_point(e[0].map(|i| mesh.nodes_m[i]), r, 1e-12, 0, 0).unwrap();
            (
                p.position_m,
                p.shapes.to_vec(),
                p.gradients.to_vec(),
                e[0].to_vec(),
            )
        }
    }
}

/// Test-only diagnostic seam. It cannot be selected by a production problem.
/// Independent Mandel strain contraction makes the hourglass rank observable.
pub fn diagnostic_hex20_stiffness(
    nodes: [[f64; 3]; 20],
    points_per_dir: usize,
    law: &TetElasticMaterial,
) -> Vec<f64> {
    let (x, w) = fs_feec::highorder::quad1d::gauss_legendre(points_per_dir);
    let mut matrix = vec![0.0; 60 * 60];
    for i in 0..x.len() {
        for j in 0..x.len() {
            for k in 0..x.len() {
                let p = hex20_point(nodes, [x[i], x[j], x[k]], 1e-12, 0, 0).unwrap();
                let mut strains = [[0.0; 6]; 60];
                for a in 0..20 {
                    let [gx, gy, gz] = p.gradients[a];
                    let s = std::f64::consts::FRAC_1_SQRT_2;
                    strains[3 * a] = [gx, 0.0, 0.0, s * gy, 0.0, s * gz];
                    strains[3 * a + 1] = [0.0, gy, 0.0, s * gx, s * gz, 0.0];
                    strains[3 * a + 2] = [0.0, 0.0, gz, 0.0, s * gy, s * gx];
                }
                let dv = w[i] * w[j] * w[k] * p.determinant_m3;
                for a in 0..60 {
                    for b in 0..60 {
                        for u in 0..6 {
                            for v in 0..6 {
                                matrix[60 * a + b] += dv
                                    * strains[a][u]
                                    * law.stiffness_mandel_pa()[u][v]
                                    * strains[b][v];
                            }
                        }
                    }
                }
            }
        }
    }
    matrix
}

/// Three-rung Richardson estimate. fs-mms owns the log-linear order fit;
/// its ladder is formed from distances to the fitted extrapolated limit.
/// Eigen residual boxes are propagated through the extrapolant separately.
pub fn richardson(values: [f64; 3], residuals: [f64; 3]) -> (f64, f64, f64) {
    let [x, y, z] = values;
    let d0 = x - y;
    let d1 = y - z;
    assert!(
        d0 * d1 > 0.0 && d0.abs() > d1.abs(),
        "contracting same-sign differences required: {values:?}"
    );
    let e = residuals;
    assert!(
        e[1] + e[2] <= 0.02 * d1.abs(),
        "eigensolver must resolve the finest difference: {e:?}, {d1}"
    );
    let s = d0 - d1;
    let es = e[0] + 2.0 * e[1] + e[2];
    assert!(s.abs() > es);
    let limit = z - d1 * d1 / s;
    let ladder = fs_mms::RefinementLadder::new(
        vec![1.0, 0.5, 0.25],
        values.map(|v| (v - limit).abs()).to_vec(),
    )
    .unwrap();
    let order = fs_mms::fit_order(&ladder).observed;
    let lo = (d1.abs() - e[1] - e[2]).powi(2) / (s.abs() + es);
    let hi = (d1.abs() + e[1] + e[2]).powi(2) / (s.abs() - es);
    let correction = d1 * d1 / s.abs();
    let uncertainty = e[2] + (hi - correction).max(correction - lo);
    (limit, (z - limit).abs() + uncertainty, order)
}
pub fn log(
    test: &str,
    tier: &str,
    family: &str,
    hs: &[f64],
    dofs: &[usize],
    observed: &[f64],
    eb: &[f64],
    order: &[f64],
    tolerance: f64,
    elapsed: f64,
) {
    let (geometry, provenance) = match tier {
        "G1" => (
            "\"extent_m\":[1.0,0.7,0.4],\"nu\":0.3",
            "exact roller-face dilatational eigenmode",
        ),
        "G3" => (
            "\"extent_m\":[1.0,0.1,0.07],\"nu\":0.25",
            "same solid operator, HEX27 and TET10 Richardson estimates",
        ),
        _ => (
            "\"extent_m\":[1.0,0.02,0.02],\"L_over_h\":50,\"nu\":0.3",
            "Timoshenko free-free, Cowper shear coefficient; Euler-Bernoulli diagnostic only",
        ),
    };
    let verdict = if hs.len() == 1 { "diagnostic" } else { "pass" };
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.3-W1b\",\"test\":\"{test}\",\"tier\":\"{tier}\",\"inputs\":{{\"family\":\"{family}\",\"mesh\":\"affine\",{geometry}}},\"provenance\":\"{provenance}\",\"rungs\":{{\"axis\":\"h\",\"values\":{hs:?},\"dofs\":{dofs:?}}},\"observed\":{{\"rel_err\":{observed:?},\"eb_rel_err\":{eb:?}}},\"order\":{order:?},\"tolerance\":{tolerance},\"verdict\":\"{verdict}\",\"elapsed_s\":{elapsed}}}"
    );
}
