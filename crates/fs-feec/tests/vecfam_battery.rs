//! Simplicial vector-family battery (bead dcng): dimension counts and
//! the exact-sequence Euler identity, dof-Kronecker unisolvence + mass
//! SPD, tangential/normal trace CONFORMITY across shared entities,
//! dd = 0 through the discrete grad/curl/div chain, interpolation
//! ladders at order r, the two-tier G3 relabeling battery (edge dofs
//! transform with definite parity; interpolated fields are label-
//! invariant), r = 1 Whitney cross-checks, and a frozen golden.

use fs_feec::highorder::simplex::SimplexSpace;
use fs_feec::highorder::vecfam::{
    DgSpace, Family, VecSpace, build_element, curl_matrix, dg_cell_dofs, div_matrix, grad_matrix,
    nedelec_entity_dofs, rt_entity_dofs, tri_quad3d,
};
use fs_feec::{deram1, deram2, element_geometry, kuhn_cube, two_tets};
use fs_rep_mesh::TetComplex;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-feec-vecfam\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

const PI: f64 = std::f64::consts::PI;

fn smooth_field(p: [f64; 3]) -> [f64; 3] {
    [
        (PI * p[1]).sin() * (PI * p[2]).cos(),
        (PI * p[2]).sin() * (PI * p[0]).cos(),
        (PI * p[0]).sin() * (PI * p[1]).cos(),
    ]
}

/// vec-001: entity/global dimension counts r = 1..4 and the
/// exact-sequence Euler identity Σ(−1)ᵏ dim Vᵏ = 1 on contractible
/// complexes (both fixtures) — the global fingerprint of exactness.
#[test]
fn vec_001_dimensions() {
    for r in 1..=4usize {
        let (ne, nf, nc) = nedelec_entity_dofs(r);
        let (rf, rc) = rt_entity_dofs(r);
        // One-tet dims against the closed forms.
        let ned_dim = 6 * ne + 4 * nf + nc;
        let rt_dim = 4 * rf + rc;
        assert_eq!(ned_dim, r * (r + 2) * (r + 3) / 2, "N_r dim r={r}");
        assert_eq!(rt_dim, r * (r + 1) * (r + 3) / 2, "RT_r dim r={r}");
        for (complex, positions, name) in [
            {
                let (c, p) = two_tets();
                (c, p, "two_tets")
            },
            {
                let (c, p) = kuhn_cube(1);
                (c, p, "kuhn1")
            },
        ] {
            let h1 = SimplexSpace::new(&complex, r);
            let ned = VecSpace::new(&complex, &positions, r, Family::Nedelec);
            let rt = VecSpace::new(&complex, &positions, r, Family::Rt);
            let dg = DgSpace::new(&complex, r);
            let euler = h1.ndof as i64 - ned.ndof as i64 + rt.ndof as i64 - dg.ndof as i64;
            assert_eq!(
                euler, 1,
                "{name} r={r}: Euler sum {euler} (H1 {} N {} RT {} DG {})",
                h1.ndof, ned.ndof, rt.ndof, dg.ndof
            );
        }
    }
    assert_eq!(dg_cell_dofs(3), 10, "P_2 dim");
    log(
        "vec-001",
        "pass",
        "dims + Euler alternating sum = 1, r=1..4",
    );
}

/// vec-002: unisolvence — every element basis reproduces its defining
/// dofs (Kronecker via independent re-application of the functionals)
/// and the global mass matrices are SPD (Cholesky succeeds).
#[test]
fn vec_002_unisolvence() {
    let (complex, positions) = two_tets();
    let mut worst = 0.0f64;
    for r in 1..=4usize {
        for family in [Family::Nedelec, Family::Rt] {
            for t in 0..complex.tets.len() {
                let el = build_element(&complex, &positions, t, r, family);
                // Re-apply the functionals through the public
                // interpolation path: interpolate each basis function
                // of THIS element (zero-extended globally).
                let space = VecSpace::new(&complex, &positions, r, family);
                let dofs = space.element_dofs(t);
                for (j, f) in el.funcs.iter().enumerate() {
                    let field = |p: [f64; 3]| f.eval_local(&el.monos, el.chart.local(p));
                    let vals = space.interpolate(&positions, &field);
                    for (i, &g) in dofs.iter().enumerate() {
                        let want = if i == j { 1.0 } else { 0.0 };
                        worst = worst.max((vals[g] - want).abs());
                    }
                }
                if t == 0 && r == 4 {
                    // one full sweep is enough at the top order
                    break;
                }
            }
        }
    }
    assert!(worst < 1e-9, "dof Kronecker worst {worst:.2e}");
    // Mass SPD on the 6-tet cube at r = 2, 3.
    let (cube, cpos) = kuhn_cube(1);
    for r in [2usize, 3] {
        for family in [Family::Nedelec, Family::Rt] {
            let sp = VecSpace::new(&cube, &cpos, r, family);
            let m = sp.mass(&cpos);
            let dense = m.to_dense();
            let f = fs_la::factor::cholesky(&dense, sp.ndof);
            assert!(f.is_ok(), "mass SPD r={r}");
        }
    }
    log(
        "vec-002",
        "pass",
        &format!("Kronecker worst {worst:.2e}; mass SPD r=2,3"),
    );
}

/// vec-003: conformity across the shared face of two_tets — for every
/// shared-entity dof's global basis function, the TANGENTIAL trace
/// (Nédélec) / NORMAL trace (RT) agrees from both sides at face
/// quadrature points. The sorted-global frame convention under test.
#[test]
fn vec_003_conformity() {
    let (complex, positions) = two_tets();
    let shared: [u32; 3] = [1, 2, 3];
    let f_idx = complex.faces.binary_search(&shared).expect("shared face");
    let tri: [[f64; 3]; 3] = core::array::from_fn(|k| positions[shared[k] as usize]);
    let e1 = [
        tri[1][0] - tri[0][0],
        tri[1][1] - tri[0][1],
        tri[1][2] - tri[0][2],
    ];
    let e2 = [
        tri[2][0] - tri[0][0],
        tri[2][1] - tri[0][1],
        tri[2][2] - tri[0][2],
    ];
    let cross = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let nrm = cross(e1, e2);
    let mut worst_t = 0.0f64;
    let mut worst_n = 0.0f64;
    for r in 2..=4usize {
        for family in [Family::Nedelec, Family::Rt] {
            let sp = VecSpace::new(&complex, &positions, r, family);
            // Global dofs living on the shared face (and, for Nédélec,
            // its edges too).
            let mut shared_dofs: Vec<usize> = Vec::new();
            if sp.per_edge > 0 {
                for pq in [[1u32, 2], [1, 3], [2, 3]] {
                    let e = complex.edges.binary_search(&pq).expect("edge");
                    for k in 0..sp.per_edge {
                        shared_dofs.push(e * sp.per_edge + k);
                    }
                }
            }
            for k in 0..sp.per_face {
                shared_dofs.push(sp.face_off + f_idx * sp.per_face + k);
            }
            for &g in &shared_dofs {
                let mut u = vec![0.0f64; sp.ndof];
                u[g] = 1.0;
                for (p, _) in tri_quad3d(tri, 3) {
                    let v0 = sp.eval_in(0, &u, p);
                    let v1 = sp.eval_in(1, &u, p);
                    let d = [v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]];
                    match family {
                        Family::Nedelec => {
                            // Tangential jump: d × n must vanish.
                            let j = cross(d, nrm);
                            let m = (j[0] * j[0] + j[1] * j[1] + j[2] * j[2]).sqrt();
                            worst_t = worst_t.max(m);
                        }
                        Family::Rt => {
                            let jn = d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2];
                            worst_n = worst_n.max(jn.abs());
                        }
                    }
                }
            }
        }
    }
    assert!(
        worst_t < 1e-9 && worst_n < 1e-9,
        "tangential jump {worst_t:.2e}, normal jump {worst_n:.2e}"
    );
    log(
        "vec-003",
        "pass",
        &format!("tangential jump {worst_t:.2e}, normal jump {worst_n:.2e}, r=2..4"),
    );
}

/// vec-004: dd = 0 through the discrete chain — ‖curl∘grad‖ and
/// ‖div∘curl‖ at roundoff (grad P_r ⊂ N_r and curl N_r ⊂ RT_r make
/// the interpolation matrices exact chain maps).
#[test]
fn vec_004_dd_zero() {
    for r in [2usize, 3] {
        let (complex, positions) = two_tets();
        let h1 = SimplexSpace::new(&complex, r);
        let ned = VecSpace::new(&complex, &positions, r, Family::Nedelec);
        let rt = VecSpace::new(&complex, &positions, r, Family::Rt);
        let dg = DgSpace::new(&complex, r);
        let d0 = grad_matrix(&h1, &ned, &positions);
        let d1 = curl_matrix(&ned, &rt, &positions);
        let d2 = div_matrix(&rt, &dg, &positions);
        // Scales for relative gates.
        let scale = |m: &fs_sparse::Csr, rows: usize, cols: usize| -> f64 {
            let mut s = 0.0f64;
            let dense = m.to_dense();
            for v in &dense {
                s = s.max(v.abs());
            }
            let _ = (rows, cols);
            s
        };
        let s0 = scale(&d0, ned.ndof, h1.ndof);
        let s1 = scale(&d1, rt.ndof, ned.ndof);
        let s2 = scale(&d2, dg.ndof, rt.ndof);
        // curl∘grad: apply to every H1 basis vector.
        let mut worst_cg = 0.0f64;
        let mut x = vec![0.0f64; h1.ndof];
        let mut y = vec![0.0f64; ned.ndof];
        let mut z = vec![0.0f64; rt.ndof];
        for j in 0..h1.ndof {
            x.iter_mut().for_each(|v| *v = 0.0);
            x[j] = 1.0;
            d0.spmv(&x, &mut y);
            d1.spmv(&y, &mut z);
            for v in &z {
                worst_cg = worst_cg.max(v.abs());
            }
        }
        // div∘curl.
        let mut worst_dc = 0.0f64;
        let mut w = vec![0.0f64; dg.ndof];
        let mut yn = vec![0.0f64; ned.ndof];
        for j in 0..ned.ndof {
            yn.iter_mut().for_each(|v| *v = 0.0);
            yn[j] = 1.0;
            d1.spmv(&yn, &mut z);
            d2.spmv(&z, &mut w);
            for v in &w {
                worst_dc = worst_dc.max(v.abs());
            }
        }
        let gate = 1e-10 * (s0 * s1).max(s1 * s2).max(1.0);
        assert!(
            worst_cg < gate && worst_dc < gate,
            "r={r}: curl.grad {worst_cg:.2e}, div.curl {worst_dc:.2e}, gate {gate:.2e}"
        );
        log(
            "vec-004",
            "pass",
            &format!("r={r}: curl.grad {worst_cg:.2e}, div.curl {worst_dc:.2e}"),
        );
    }
}

/// vec-005: canonical-interpolation ladders — smooth-field
/// interpolation error decays at order r on the Kuhn refinement
/// ladder for both families at r = 2, 3.
#[test]
fn vec_005_interpolation_ladder() {
    for r in [2usize, 3] {
        for family in [Family::Nedelec, Family::Rt] {
            let mut errs = Vec::new();
            for n in [2usize, 3] {
                let (complex, positions) = kuhn_cube(n);
                let sp = VecSpace::new(&complex, &positions, r, family);
                let u = sp.interpolate(&positions, &smooth_field);
                errs.push(sp.l2_error(&positions, &u, &smooth_field));
            }
            let slope = (errs[0] / errs[1]).ln() / (3.0f64 / 2.0).ln();
            let fam = if matches!(family, Family::Nedelec) {
                "N"
            } else {
                "RT"
            };
            assert!(
                slope > r as f64 - 0.45 && errs[1] < errs[0],
                "{fam} r={r}: errs {errs:?} slope {slope:.2}"
            );
            log(
                "vec-005",
                "pass",
                &format!(
                    "{fam} r={r}: errs {:.3e}/{:.3e} slope {slope:.2}",
                    errs[0], errs[1]
                ),
            );
        }
    }
}

/// vec-006 (G3): the two-tier orientation battery. Relabeling global
/// vertices changes each entity's sorted order, so entity bases are
/// NOT naively label-invariant — instead: (tier 1, signed-permutation
/// level) edge dofs transform with DEFINITE PARITY (−1)^{k+1} when the
/// edge direction flips (tangent flips AND P_k(−s) = (−1)^k P_k(s));
/// (tier 2, physics level, where face/interior bases mix) the
/// canonically interpolated FIELD of a fixed analytic function is
/// pointwise label-invariant.
#[test]
fn vec_006_relabeling() {
    let (complex, positions) = two_tets();
    let perm: [u32; 5] = [3, 0, 4, 1, 2];
    let mut new_pos = vec![[0.0f64; 3]; 5];
    for (old, &np) in perm.iter().enumerate() {
        new_pos[np as usize] = positions[old];
    }
    let new_tets: Vec<[u32; 4]> = complex
        .tets
        .iter()
        .map(|t| {
            let mut nt = [0u32; 4];
            for (k, &v) in t.iter().enumerate() {
                nt[k] = perm[v as usize];
            }
            nt
        })
        .collect();
    let complex2 = TetComplex::from_tets(5, new_tets);
    let mut worst_sign = 0.0f64;
    let mut worst_field = 0.0f64;
    for r in [2usize, 3] {
        for family in [Family::Nedelec, Family::Rt] {
            let sp1 = VecSpace::new(&complex, &positions, r, family);
            let sp2 = VecSpace::new(&complex2, &new_pos, r, family);
            let u1 = sp1.interpolate(&positions, &smooth_field);
            let u2 = sp2.interpolate(&new_pos, &smooth_field);
            // Tier 1: edge dofs are a SIGNED permutation.
            if sp1.per_edge > 0 {
                for (e1, &[a, b]) in complex.edges.iter().enumerate() {
                    let (na, nb) = (perm[a as usize], perm[b as usize]);
                    let flipped = na > nb;
                    let key = if flipped { [nb, na] } else { [na, nb] };
                    let e2 = complex2.edges.binary_search(&key).expect("edge");
                    for k in 0..sp1.per_edge {
                        let sign = if flipped {
                            if k % 2 == 0 { -1.0 } else { 1.0 } // (−1)^{k+1}
                        } else {
                            1.0
                        };
                        let d =
                            (u2[e2 * sp2.per_edge + k] - sign * u1[e1 * sp1.per_edge + k]).abs();
                        worst_sign = worst_sign.max(d);
                    }
                }
            }
            // Tier 2: the interpolated FIELD is label-invariant.
            for t in 0..complex.tets.len() {
                let tet1 = complex.tets[t];
                for lam in [
                    [0.4f64, 0.3, 0.2, 0.1],
                    [0.1, 0.2, 0.3, 0.4],
                    [0.25, 0.25, 0.25, 0.25],
                ] {
                    let mut p = [0.0f64; 3];
                    for (a, &v) in tet1.iter().enumerate() {
                        for k in 0..3 {
                            p[k] += lam[a] * positions[v as usize][k];
                        }
                    }
                    let v1 = sp1.eval_in(t, &u1, p);
                    let v2 = sp2.eval_in(t, &u2, p);
                    for k in 0..3 {
                        worst_field = worst_field.max((v1[k] - v2[k]).abs());
                    }
                }
            }
        }
    }
    assert!(
        worst_sign < 1e-10 && worst_field < 1e-9,
        "signed-permutation dev {worst_sign:.2e}, field dev {worst_field:.2e}"
    );
    log(
        "vec-006",
        "pass",
        &format!(
            "edge signed-permutation dev {worst_sign:.2e}; physics-level field dev {worst_field:.2e}"
        ),
    );
}

/// vec-007: r = 1 members coincide with the Whitney forms — the
/// canonical interpolants match deram1/deram2 dof-for-dof, and a
/// frozen golden pins the r = 2 construction bits.
#[test]
fn vec_007_whitney_and_golden() {
    let (complex, positions) = two_tets();
    let ned = VecSpace::new(&complex, &positions, 1, Family::Nedelec);
    let rt = VecSpace::new(&complex, &positions, 1, Family::Rt);
    let a = |p: [f64; 3]| -> [f64; 3] { [p[1], p[2], p[0]] };
    let mine = ned.interpolate(&positions, &a);
    let whitney = deram1(&complex, &positions, &a);
    let mut worst = 0.0f64;
    for (m, w) in mine.iter().zip(&whitney) {
        worst = worst.max((m - w).abs());
    }
    // RT dofs are area-normalized moments; deram2 is the raw flux.
    let b = |p: [f64; 3]| -> [f64; 3] { [p[0], -p[1], 2.0 * p[2]] };
    let mine_rt = rt.interpolate(&positions, &b);
    let flux = deram2(&complex, &positions, &b);
    let geo = element_geometry(&complex, &positions);
    let _ = geo;
    let mut worst_rt = 0.0f64;
    for (f_idx, tri) in complex.faces.iter().enumerate() {
        let p: [[f64; 3]; 3] = core::array::from_fn(|k| positions[tri[k] as usize]);
        let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
        let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
        let c = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let area = 0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        // vecfam face dof = mean normal moment = flux / area.
        worst_rt = worst_rt.max((mine_rt[f_idx] - flux[f_idx] / area).abs());
    }
    assert!(
        worst < 1e-11 && worst_rt < 1e-11,
        "whitney cross-check: edge {worst:.2e}, face {worst_rt:.2e}"
    );
    // Golden: FNV over r=2 mass-matrix action bits.
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    for family in [Family::Nedelec, Family::Rt] {
        let sp = VecSpace::new(&complex, &positions, 2, family);
        let m = sp.mass(&positions);
        let x: Vec<f64> = (0..sp.ndof)
            .map(|i| ((i * 2_654_435_761 + 12_345) % 1000) as f64 / 1000.0 - 0.5)
            .collect();
        let mut y = vec![0.0f64; sp.ndof];
        m.spmv(&x, &mut y);
        for v in y.iter().step_by(7) {
            feed(*v);
        }
    }
    log("vec-007-golden", "info", &format!("{acc:#018x}"));
    log(
        "vec-007",
        "pass",
        &format!("whitney edge {worst:.2e} face {worst_rt:.2e}; golden {acc:#018x}"),
    );
    // Frozen on this platform (bead dcng); the cross-ISA confirmation
    // is LEDGERED PENDING, same policy as the simplex golden.
    const GOLDEN_HASH: u64 = 0x2737_401f_72b4_99f2;
    assert_eq!(
        acc, GOLDEN_HASH,
        "vecfam bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
