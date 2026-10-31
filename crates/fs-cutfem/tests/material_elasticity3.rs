//! Native material energy, inventory and heterogeneous stabilization derivatives.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::op::LinearOp;

struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.0
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(self.0, self.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
fn state(scale: f64, density: f64, id: u8) -> ElasticAcousticState {
    let mut c = [[0.0; 6]; 6];
    for (i, row) in c.iter_mut().enumerate() {
        row[i] = scale * (15.0 + i as f64);
    }
    // Coupled normal/shear terms distinguish a full tensor from an E/nu proxy.
    for (i, j, v) in [(0, 3, 1.2), (1, 4, -0.7), (2, 5, 0.8), (3, 5, 0.4)] {
        c[i][j] = scale * v;
        c[j][i] = scale * v;
    }
    ElasticAcousticState::try_new_mandel(c, density, ContentHash([id; 32]), 0.1).unwrap()
}
fn build(field: ElasticMaterialField3<'_>, cut: f64, gamma: f64) -> CutElasticity3 {
    let mut checkpoint = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    CutElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [2; 3],
        &Slab(cut),
        field,
        &|p| p == [0.0; 3],
        ElasticityOptions3 {
            ghost_gamma: gamma,
            ..Default::default()
        },
        &mut control,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}
fn bilinear(op: &CutElasticity3, u: &[f64], z: &[f64]) -> f64 {
    let mut ku = vec![0.0; op.n()];
    op.apply(u, &mut ku);
    z.iter().zip(ku).map(|(z, ku)| z * ku).sum()
}

#[test]
fn coupled_tensor_matches_affine_continuum_energy_and_retains_mass_and_identity() {
    let native = state(200.0, 1700.0, 1);
    let op = build(ElasticMaterialField3::Uniform(&native), 1.1, 0.1);
    let a = [
        [0.012, 0.003, -0.002],
        [0.001, -0.004, 0.005],
        [0.004, -0.001, 0.007],
    ];
    let u: Vec<f64> = op
        .nodes()
        .iter()
        .flat_map(|p| a.map(|row| row.iter().zip(p).map(|(a, p)| a * p).sum::<f64>()))
        .collect();
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let e = [
        a[0][0],
        a[1][1],
        a[2][2],
        (a[0][1] + a[1][0]) * q,
        (a[1][2] + a[2][1]) * q,
        (a[0][2] + a[2][0]) * q,
    ];
    let expected: f64 = native
        .stiffness_mandel_pa()
        .iter()
        .enumerate()
        .map(|(i, row)| e[i] * row.iter().zip(e).map(|(c, e)| c * e).sum::<f64>())
        .sum();
    close(bilinear(&op, &u, &u), expected, 2e-13);
    close(
        op.reference_cell_masses_kg().unwrap().iter().sum(),
        1700.0,
        2e-13,
    );
    for i in 0..op.cells() {
        assert_eq!(op.cell_material_state(i), Some(&native));
    }
}

#[test]
fn heterogeneous_cut_interfaces_have_exact_side_specific_scale_derivatives() {
    let left = state(1.0, 2.0, 2);
    let right = state(4.0, 9.0, 3);
    let field = |bounds: HexCell| {
        Ok(if bounds.lo()[0] < 0.5 {
            left.clone()
        } else {
            right.clone()
        })
    };
    let mut op = build(ElasticMaterialField3::ByCell(&field), 0.73, 0.2);
    let base: Vec<_> = (0..op.cells()).map(|i| 0.41 + 0.03 * i as f64).collect();
    op.set_scales(&base).unwrap();
    let u: Vec<_> = (0..op.n())
        .map(|i| {
            if op.fixed()[i / 3] {
                0.0
            } else {
                (i % 7) as f64 / 17.0
            }
        })
        .collect();
    let z: Vec<_> = (0..op.n())
        .map(|i| {
            if op.fixed()[i / 3] {
                0.0
            } else {
                (i % 11) as f64 / 19.0
            }
        })
        .collect();
    let work = op
        .scale_bilinear_forms(&z, &u, || ControlFlow::Continue(()))
        .unwrap();
    let energy = op.scale_quadratic_forms(&u).unwrap();
    close(
        bilinear(&op, &u, &u),
        energy.iter().zip(&base).map(|(e, s)| e * s).sum(),
        2e-13,
    );
    close(
        bilinear(&op, &u, &z),
        work.iter().zip(&base).map(|(w, s)| w * s).sum(),
        2e-13,
    );
    for cell in 0..op.cells() {
        let h = 1e-5;
        let mut scales = base.clone();
        scales[cell] += h;
        op.set_scales(&scales).unwrap();
        let plus = bilinear(&op, &u, &z);
        scales[cell] -= 2.0 * h;
        op.set_scales(&scales).unwrap();
        let minus = bilinear(&op, &u, &z);
        close((plus - minus) / (2.0 * h), work[cell], 2e-8);
    }
    op.set_scales(&base).unwrap();
    let mut unstabilized = build(ElasticMaterialField3::ByCell(&field), 0.73, 0.0);
    unstabilized.set_scales(&base).unwrap();
    assert!(bilinear(&op, &u, &u) > bilinear(&unstabilized, &u, &u) + 1e-3);
    // Material density weights physical inventory; cut measure remains geometric.
    close(op.volumes().iter().sum(), 0.73, 2e-13);
    close(
        op.reference_cell_masses_kg().unwrap().iter().sum(),
        0.73 * (2.0 + 9.0) / 2.0,
        2e-13,
    );
}
