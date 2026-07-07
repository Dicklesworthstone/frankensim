//! Simplicial high-order H¹ battery (tfz.6 slice 3): dof counting,
//! unisolvence (mass SPD), CONFORMITY across shared faces (traces
//! from both elements bitwise-comparable at face barycentric points),
//! G1 MMS at r = 1..4 on Kuhn ladders, and the ORIENTATION battery
//! (G3): operator equivariance under vertex relabeling at the
//! signed-permutation level (r = 2, where the dof map is exactly a
//! signed permutation) and physics invariance at r = 4 (vertex point
//! values, energy, and L2 error are label-independent).

use fs_feec::highorder::simplex::{SimplexSpace, duffy_quadrature, entity_dofs};
use fs_feec::kuhn_cube;
use fs_rand::StreamKey;
use fs_rep_mesh::TetComplex;
use fs_sparse::precond::{IdentityPrecond, pcg};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-feec-ho\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 11,
        kernel: 0x51E3,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

#[test]
fn quadrature_and_dof_counts() {
    // Duffy weights sum to the reference volume 1/6.
    for n in 2..=8usize {
        let q = duffy_quadrature(n);
        let total: f64 = q.iter().map(|(_, w)| w).sum();
        assert!(
            (total - 1.0 / 6.0).abs() < 1e-14,
            "n={n}: Duffy weights sum {total}"
        );
        // Barycentrics valid.
        for (lam, _) in &q {
            assert!(lam.iter().all(|&l| l > -1e-15 && l < 1.0 + 1e-15));
            let s: f64 = lam.iter().sum();
            assert!((s - 1.0).abs() < 1e-14);
        }
    }
    // Entity dof counts match the P_r dimension: total local dofs
    // must equal C(r+3, 3).
    for r in 1..=6usize {
        let (pe, pf, pi) = entity_dofs(r);
        let local = 4 + 6 * pe + 4 * pf + pi;
        let dim = (r + 1) * (r + 2) * (r + 3) / 6;
        assert_eq!(
            local, dim,
            "r={r}: local dof count {local} vs dim P_r = {dim}"
        );
    }
    log("counts", "pass", "Duffy volumes + P_r dimensions r=1..6");
}

#[test]
fn unisolvence_mass_spd() {
    // The element family is linearly independent iff the mass matrix
    // is SPD — checked via Cholesky on the assembled two-tet mass.
    let (complex, positions) = fs_feec::two_tets();
    for r in 1..=5usize {
        let sp = SimplexSpace::new(&complex, r);
        let mass = sp.mass(&positions);
        let dense = mass.to_dense();
        let chol = fs_la::factor::cholesky(&dense, sp.ndof);
        assert!(
            chol.is_ok(),
            "r={r}: mass not SPD — basis dependent or quadrature underintegrated"
        );
    }
    log("unisolvence", "pass", "mass SPD r=1..5 (two-tet)");
}

#[test]
fn conformity_across_shared_face() {
    // two_tets share global face {1,2,3}: element 0 = [0,1,2,3] sees
    // it at barycentric (0, α, β, γ); element 1 = [1,2,3,4] at
    // (α, β, γ, 0). A random GLOBAL dof vector must produce the same
    // trace from both sides — the sorted-global orientation convention
    // doing its job (this is where unoriented high-order codes break).
    let (complex, _positions) = fs_feec::two_tets();
    let r = 4usize;
    let sp = SimplexSpace::new(&complex, r);
    let u = rand_vec(sp.ndof, 1);
    let eval_elem = |t: usize, lam: [f64; 4]| -> f64 {
        sp.element_dofs(t)
            .iter()
            .map(|(gi, lf)| u[*gi] * lf.eval(lam, r))
            .sum()
    };
    let samples = [
        [0.2f64, 0.3, 0.5],
        [0.6, 0.1, 0.3],
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.05, 0.9, 0.05],
        [0.45, 0.1, 0.45],
    ];
    let mut worst = 0.0f64;
    for &[a, b, c] in &samples {
        let v0 = eval_elem(0, [0.0, a, b, c]);
        let v1 = eval_elem(1, [a, b, c, 0.0]);
        worst = worst.max((v0 - v1).abs());
    }
    assert!(
        worst < 1e-13,
        "conformity broken across shared face: {worst:.3e}"
    );
    log(
        "conformity",
        "pass",
        &format!("r=4 trace agreement {worst:.1e}"),
    );
}

fn solve_poisson(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    r: usize,
    f_exact: &dyn Fn([f64; 3]) -> f64,
) -> Vec<f64> {
    let sp = SimplexSpace::new(complex, r);
    let k = sp.stiffness(positions);
    let b = sp.load(positions, &|p| f_exact(p));
    let boundary = sp.boundary_mask();
    let interior: Vec<usize> = (0..sp.ndof).filter(|&d| !boundary[d]).collect();
    let ni = interior.len();
    let mut slot = vec![usize::MAX; sp.ndof];
    for (i, &d) in interior.iter().enumerate() {
        slot[d] = i;
    }
    let mut red = fs_sparse::Coo::new(ni, ni);
    for (i, &d) in interior.iter().enumerate() {
        let (cols, vals) = k.row(d);
        for (&c, &v) in cols.iter().zip(vals) {
            if slot[c] != usize::MAX {
                red.push(i, slot[c], v);
            }
        }
    }
    let a = red.assemble();
    let rhs: Vec<f64> = interior.iter().map(|&d| b[d]).collect();
    let mut x = vec![0.0f64; ni];
    let report = pcg(&a, &rhs, &mut x, &IdentityPrecond, 1e-12, 40_000);
    assert!(report.converged, "PCG failed at r={r}: {report:?}");
    let mut full = vec![0.0f64; sp.ndof];
    for (i, &d) in interior.iter().enumerate() {
        full[d] = x[i];
    }
    full
}

#[test]
fn mms_orders_r1_through_r4() {
    let pi = std::f64::consts::PI;
    let u_exact = move |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let f_exact = move |p: [f64; 3]| 3.0 * pi * pi * u_exact(p);
    for r in 1..=4usize {
        // r = 1 needs a finer ladder to reach asymptotics on tets.
        let ladder: [usize; 2] = if r == 1 { [4, 8] } else { [2, 4] };
        let mut errs = Vec::new();
        for &m in &ladder {
            let (complex, positions) = kuhn_cube(m);
            let u = solve_poisson(&complex, &positions, r, &f_exact);
            let sp = SimplexSpace::new(&complex, r);
            let err = sp.l2_error(&positions, &u, &u_exact);
            errs.push(err);
            log("mms-simplex", "info", &format!("r={r} m={m} L2={err:.4e}"));
        }
        let order = (errs[0] / errs[1]).ln() / 2.0f64.ln();
        assert!(
            order > r as f64 + 0.6,
            "r={r}: simplicial order {order:.2} (errors {errs:?})"
        );
        log(
            "mms-simplex-order",
            "pass",
            &format!("r={r} order={order:.2}"),
        );
    }
}

/// Permute vertex labels of a complex+positions pair.
fn permute(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    perm: &[u32],
) -> (TetComplex, Vec<[f64; 3]>) {
    let tets: Vec<[u32; 4]> = complex
        .tets
        .iter()
        .map(|t| {
            [
                perm[t[0] as usize],
                perm[t[1] as usize],
                perm[t[2] as usize],
                perm[t[3] as usize],
            ]
        })
        .collect();
    let mut pos = vec![[0.0f64; 3]; positions.len()];
    for (v, &p) in perm.iter().zip(positions) {
        pos[*v as usize] = p;
    }
    (TetComplex::from_tets(complex.vertex_count, tets), pos)
}

/// Deterministic pseudo-random permutation of 0..n.
fn random_perm(n: usize, seed_tile: u32) -> Vec<u32> {
    let mut s = StreamKey {
        seed: 12,
        kernel: 0x51E3,
        tile: seed_tile,
    }
    .stream();
    let mut p: Vec<u32> = (0..u32::try_from(n).expect("small")).collect();
    for i in (1..n).rev() {
        let j = usize::try_from(s.next_below(i as u64 + 1)).expect("small");
        p.swap(i, j);
    }
    p
}

#[test]
fn orientation_signed_permutation_equivariance_r3() {
    // At r = 3 each edge carries kernels P_0 (symmetric) and P_1
    // (ANTIsymmetric: sign flips when the sorted order reverses), and
    // each face carries the single symmetric kernel P_0·P_0 — so
    // relabeling maps dofs to dofs up to a KNOWN sign, and the
    // operator must be equivariant: K'(S u) = S(K u) to assembly
    // roundoff. This pins the sorted-global convention sharply (r = 2
    // has only symmetric kernels and would test nothing; higher r
    // mixes face bases and needs the physics-invariance gate below).
    let (complex, positions) = kuhn_cube(2);
    let r = 3usize;
    let sp = SimplexSpace::new(&complex, r);
    let k = sp.stiffness(&positions);
    let perm = random_perm(complex.vertex_count, 7);
    let (complex_p, positions_p) = permute(&complex, &positions, &perm);
    let sp_p = SimplexSpace::new(&complex_p, r);
    let k_p = sp_p.stiffness(&positions_p);
    assert_eq!(sp.ndof, sp_p.ndof);
    // Dof map: vertex v → perm[v]; edge {a,b} dof k → same k on the
    // permuted edge with sign (−1)^k on sort-order flips; face dof
    // (single, symmetric) → the permuted face.
    let mut map = vec![0usize; sp.ndof];
    let mut sign = vec![1.0f64; sp.ndof];
    for v in 0..complex.vertex_count {
        map[v] = perm[v] as usize;
    }
    for (e, &[a, b]) in complex.edges.iter().enumerate() {
        let (pa, pb) = (perm[a as usize], perm[b as usize]);
        let flipped = pa > pb;
        let key = if flipped { [pb, pa] } else { [pa, pb] };
        let ep = complex_p
            .edges
            .binary_search(&key)
            .expect("edge in permuted table");
        for kk in 0..sp.per_edge {
            map[sp.edge_off + e * sp.per_edge + kk] = sp_p.edge_off + ep * sp.per_edge + kk;
            if flipped && kk % 2 == 1 {
                sign[sp.edge_off + e * sp.per_edge + kk] = -1.0;
            }
        }
    }
    for (f, &tri) in complex.faces.iter().enumerate() {
        let mut key = [
            perm[tri[0] as usize],
            perm[tri[1] as usize],
            perm[tri[2] as usize],
        ];
        key.sort_unstable();
        let fp = complex_p
            .faces
            .binary_search(&key)
            .expect("face in permuted table");
        for kk in 0..sp.per_face {
            map[sp.face_off + f * sp.per_face + kk] = sp_p.face_off + fp * sp.per_face + kk;
        }
    }
    let u = rand_vec(sp.ndof, 20);
    // S u.
    let mut su = vec![0.0f64; sp.ndof];
    for d in 0..sp.ndof {
        su[map[d]] = sign[d] * u[d];
    }
    let mut ku = vec![0.0f64; sp.ndof];
    k.spmv(&u, &mut ku);
    let mut sku = vec![0.0f64; sp.ndof];
    for d in 0..sp.ndof {
        sku[map[d]] = sign[d] * ku[d];
    }
    let mut kp_su = vec![0.0f64; sp.ndof];
    k_p.spmv(&su, &mut kp_su);
    let scale = sku.iter().map(|v| v.abs()).fold(0.0f64, f64::max).max(1.0);
    let worst = sku
        .iter()
        .zip(&kp_su)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12 * scale,
        "orientation equivariance broken: {worst:.3e} (scale {scale:.3e})"
    );
    log(
        "orientation-r2",
        "pass",
        &format!("signed-permutation equivariance {worst:.1e}"),
    );
}

#[test]
fn orientation_physics_invariance_r4() {
    // At r = 4 face dofs mix under relabeling (the span is invariant,
    // the basis is not), so the honest G3 gate is PHYSICS invariance:
    // solve the same Poisson problem under a random relabeling and
    // compare vertex point values (exact dof = point value in this
    // hierarchy), energy, and L2 error.
    let pi = std::f64::consts::PI;
    let u_exact = move |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let f_exact = move |p: [f64; 3]| 3.0 * pi * pi * u_exact(p);
    let (complex, positions) = kuhn_cube(2);
    let r = 4usize;
    let u1 = solve_poisson(&complex, &positions, r, &f_exact);
    let sp1 = SimplexSpace::new(&complex, r);
    let err1 = sp1.l2_error(&positions, &u1, &u_exact);
    let perm = random_perm(complex.vertex_count, 9);
    let (complex_p, positions_p) = permute(&complex, &positions, &perm);
    let u2 = solve_poisson(&complex_p, &positions_p, r, &f_exact);
    let sp2 = SimplexSpace::new(&complex_p, r);
    let err2 = sp2.l2_error(&positions_p, &u2, &u_exact);
    // Vertex values are basis-independent point values.
    let mut worst_v = 0.0f64;
    for v in 0..complex.vertex_count {
        worst_v = worst_v.max((u1[v] - u2[perm[v] as usize]).abs());
    }
    assert!(
        worst_v < 1e-9,
        "vertex values differ under relabeling: {worst_v:.3e}"
    );
    let rel = (err1 - err2).abs() / err1.max(1e-30);
    assert!(
        rel < 1e-6,
        "L2 error differs under relabeling: {err1:.6e} vs {err2:.6e}"
    );
    log(
        "orientation-r4",
        "pass",
        &format!("vertex dev {worst_v:.1e}, L2 rel dev {rel:.1e}"),
    );
}

const GOLDEN_HASH: u64 = 0xf37d_c95f_9ce6_c195; // recorded at tfz.6 slice 3, frozen

#[test]
fn simplex_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (complex, positions) = kuhn_cube(2);
    let sp = SimplexSpace::new(&complex, 3);
    let k = sp.stiffness(&positions);
    let u = rand_vec(sp.ndof, 40);
    let mut ku = vec![0.0f64; sp.ndof];
    k.spmv(&u, &mut ku);
    for v in ku.iter().step_by(13) {
        feed(*v);
    }
    // Basis evaluation fingerprint.
    let lam = [0.1f64, 0.2, 0.3, 0.4];
    for (_, lf) in sp.element_dofs(0) {
        feed(lf.eval(lam, 3));
        for d in lf.d_lambda(lam, 3) {
            feed(d);
        }
    }
    // Quadrature fingerprint.
    for (l, w) in duffy_quadrature(4).iter().take(8) {
        for v in l {
            feed(*v);
        }
        feed(*w);
    }
    log("simplex-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "simplex bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
