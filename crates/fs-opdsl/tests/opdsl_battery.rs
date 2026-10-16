//! fs-opdsl battery (tfz.4): the acceptance gates — generated primal
//! vs hand FEEC (materialized BITWISE, matrix-free to roundoff);
//! ⟨Av,w⟩ = ⟨v,Aᵀw⟩ mechanically for every fixture (nonsymmetric
//! advection makes it non-vacuous); JVP-vs-VJP transpose consistency;
//! generated derivatives vs fs-ad dual numbers through the pointwise
//! law; dd = 0 honored AT THE IR LEVEL; type errors structurally
//! rejected; generation determinism; DWR indicators; MMS order ≈ 2
//! THROUGH the generated operator; elasticity constitutive
//! integration (affine patch + symmetry); perf comparison logged; and
//! the cross-ISA golden hash.

use fs_feec::{element_geometry, kuhn_cube};
use fs_opdsl::fixtures::{
    VECTOR_DOFS, convection_diffusion, elasticity, poisson, reaction_diffusion,
};
use fs_opdsl::{Atom, Expr, OperatorDef, Space, TypeError, mms_poisson_study};
use fs_qty::Dims;
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-opdsl\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 4,
        kernel: 0x0D51,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn poisson_matches_hand_feec() {
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let (def, expr) = poisson(&complex, &geo);
    let op = def.lower(expr);
    // MATERIALIZED: same spgemm association as fs-feec's stiffness —
    // bitwise-equal dense forms.
    let generated = op.materialize().expect("linear").to_dense();
    let hand = fs_feec::stiffness(
        &fs_feec::incidence_to_csr(&complex.d0()),
        &fs_feec::mass_matrix(&complex, &geo, 1),
    )
    .to_dense();
    assert_eq!(generated.len(), hand.len());
    let mut worst_bits = 0u64;
    for (g, h) in generated.iter().zip(&hand) {
        worst_bits = worst_bits.max(g.to_bits() ^ h.to_bits());
    }
    assert_eq!(
        worst_bits, 0,
        "materialized Poisson must equal hand FEEC bitwise"
    );
    // MATRIX-FREE: different association (three SpMVs) — to roundoff.
    let u = rand_vec(complex.vertex_count, 1);
    let y_free = op.apply(&u);
    let mut y_hand = vec![0.0f64; u.len()];
    fs_feec::stiffness(
        &fs_feec::incidence_to_csr(&complex.d0()),
        &fs_feec::mass_matrix(&complex, &geo, 1),
    )
    .spmv(&u, &mut y_hand);
    let scale: f64 = y_hand.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let worst = y_free
        .iter()
        .zip(&y_hand)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-13 * scale.max(1.0),
        "matrix-free deviates {worst:.3e}"
    );
    log(
        "poisson-vs-hand",
        "pass",
        &format!("bitwise materialized; free dev {worst:.2e}"),
    );
}

#[test]
fn adjoint_identity_all_fixtures() {
    // ⟨Av, w⟩ = ⟨v, Aᵀw⟩ over random vectors — for Poisson
    // (symmetric), convection–diffusion (NONsymmetric: the gate does
    // real work), and elasticity (declared symmetric, verified here).
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let cases: Vec<(&str, OperatorDef, Expr)> = {
        let (pd, pe) = poisson(&complex, &geo);
        let (cd, ce) = convection_diffusion(&complex, &geo, 0.7, [1.0, -0.5, 0.25]);
        let (ed, ee) = elasticity(&complex, &geo, 200.0e9, 0.3);
        vec![
            ("poisson", pd, pe),
            ("convection-diffusion", cd, ce),
            ("elasticity", ed, ee),
        ]
    };
    for (name, def, expr) in &cases {
        let op = def.lower(expr.clone());
        let n = def.field_space.n;
        for trial in 0..3u32 {
            let v = rand_vec(n, 100 + trial);
            let w = rand_vec(n, 200 + trial);
            let av = op.jvp(&v);
            let atw = op.vjp(&w);
            let lhs = dot(&av, &w);
            let rhs = dot(&v, &atw);
            let scale = lhs.abs().max(rhs.abs()).max(1e-30);
            assert!(
                ((lhs - rhs) / scale).abs() < 1e-12,
                "{name}: adjoint identity broken: {lhs:.17e} vs {rhs:.17e}"
            );
        }
        log("adjoint-identity", "pass", name);
    }
}

#[test]
fn nonlinear_jvp_vjp_and_dual_check() {
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let (def, expr) = reaction_diffusion(&complex, &geo, 0.8);
    let mut op = def.lower(expr);
    let n = complex.vertex_count;
    let u0 = rand_vec(n, 10);
    op.linearize(&u0);
    // Transpose consistency of the LINEARIZED operator.
    for trial in 0..3u32 {
        let v = rand_vec(n, 300 + trial);
        let w = rand_vec(n, 400 + trial);
        let lhs = dot(&op.jvp(&v), &w);
        let rhs = dot(&v, &op.vjp(&w));
        let scale = lhs.abs().max(rhs.abs()).max(1e-30);
        assert!(
            ((lhs - rhs) / scale).abs() < 1e-12,
            "nonlinear transpose identity broken: {lhs:.17e} vs {rhs:.17e}"
        );
    }
    // Directional derivative vs fs-ad Dual64 through the WHOLE
    // residual: phi(t) = <c, R(u0 + t v)>; phi'(0) must equal
    // <c, J v>. Dual drives R elementwise through the law only (the
    // linear part is exact), so compare against the law's pointwise
    // chain: R(u) = K u + M0 N(u) => phi'(0) = <c, K v> + <c, M0 diag(N'(u0)) v>.
    let v = rand_vec(n, 500);
    let c = rand_vec(n, 501);
    let jv = op.jvp(&v);
    let analytic = dot(&c, &jv);
    // Dual check on the pointwise composition: for each dof,
    // N(u0 + eps v) through Dual64<1> gives dN = N'(u0)*v exactly.
    let law = fs_opdsl::CubicReaction { alpha: 0.8 };
    let mut dn = vec![0.0f64; n];
    for i in 0..n {
        let (_, d) = fs_ad::jvp([u0[i]], [v[i]], |x| {
            let a = fs_ad::Dual64::<1>::constant(0.8);
            a * x[0] * x[0] * x[0]
        });
        dn[i] = d;
        // The law's hand derivative must agree with the dual number.
        let hand = {
            use fs_opdsl::PointwiseLaw;
            law.derivative(u0[i]) * v[i]
        };
        assert!(
            (d - hand).abs() <= 1e-12 * hand.abs().max(1.0),
            "law derivative disagrees with Dual64 at dof {i}: {d} vs {hand}"
        );
    }
    // Rebuild <c, J v> from parts: K v + M0 dn.
    let (pdef, pexpr) = poisson(&complex, &geo);
    let pop = pdef.lower(pexpr);
    let kv = pop.apply(&v);
    let m0 = fs_feec::mass_matrix(&complex, &geo, 0);
    let mut m0dn = vec![0.0f64; n];
    m0.spmv(&dn, &mut m0dn);
    let reference = dot(&c, &kv) + dot(&c, &m0dn);
    let scale = analytic.abs().max(reference.abs()).max(1e-30);
    assert!(
        ((analytic - reference) / scale).abs() < 1e-12,
        "generated JVP disagrees with dual-number reference: {analytic:.17e} vs {reference:.17e}"
    );
    log(
        "nonlinear-dual",
        "pass",
        "chain rule == Dual64 through the cubic law",
    );
}

#[test]
fn dd_zero_at_the_ir_level() {
    let (complex, positions) = kuhn_cube(1);
    let geo = element_geometry(&complex, &positions);
    let dims = Dims::NONE;
    let space = Space {
        degree: 0,
        n: complex.vertex_count,
        dims,
    };
    let mut def = OperatorDef::new(space);
    let d0 = def.add_atom(Atom::d(&complex, 0, dims));
    let d1 = def.add_atom(Atom::d(&complex, 1, dims));
    let _ = geo;
    // d1(d0(u)) folds to Zero BEFORE any float exists.
    let e = def
        .apply(d0, def.field())
        .and_then(|e| def.apply(d1, e))
        .expect("well-typed");
    assert!(
        matches!(e, Expr::Zero(_)),
        "d after d must fold to Zero at the IR level"
    );
    // Through a scale too: d1(s·d0(u)) = 0.
    let e2 = def
        .apply(d0, def.field())
        .map(|e| def.scale(3.5, e))
        .and_then(|e| def.apply(d1, e))
        .expect("well-typed");
    assert!(
        matches!(e2, Expr::Zero(_)),
        "scale must not hide the exactness identity"
    );
    // The materialized zero applies as zero.
    let op = def.lower(e);
    let u = rand_vec(complex.vertex_count, 20);
    assert!(op.apply(&u).iter().all(|&v| v == 0.0));
    log("dd-zero-ir", "pass", "exactness folded symbolically");
}

#[test]
fn type_errors_are_structural() {
    let (complex, positions) = kuhn_cube(1);
    let geo = element_geometry(&complex, &positions);
    let dims = Dims::NONE;
    let secs = Dims([0, 0, 1, 0, 0, 0]);
    let space = Space {
        degree: 0,
        n: complex.vertex_count,
        dims,
    };
    let mut def = OperatorDef::new(space);
    let d0 = def.add_atom(Atom::d(&complex, 0, dims));
    let d1 = def.add_atom(Atom::d(&complex, 1, dims));
    let m0 = def.add_atom(Atom::mass(&complex, &geo, 0, dims));
    // Degree mismatch: d1 applied to a 0-cochain.
    let err = def.apply(d1, def.field()).unwrap_err();
    assert!(
        matches!(err, TypeError::ApplyMismatch { .. }),
        "degree mismatch must be caught"
    );
    // Dims mismatch: an atom expecting seconds fed a dimensionless field.
    let d0_secs = def.add_atom(Atom::d(&complex, 0, secs));
    let err = def.apply(d0_secs, def.field()).unwrap_err();
    assert!(
        matches!(err, TypeError::ApplyMismatch { .. }),
        "dims mismatch must be caught"
    );
    // Add mismatch: edge cochain + vertex cochain.
    let e_edges = def.apply(d0, def.field()).expect("typed");
    let e_verts = def.apply(m0, def.field()).expect("typed");
    let err = def.add(e_edges, e_verts).unwrap_err();
    assert!(matches!(err, TypeError::AddMismatch { .. }));
    log(
        "type-errors",
        "pass",
        "degree, dims, and add mismatches all structural",
    );
}

#[test]
fn generation_is_deterministic_and_reports_provenance() {
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let build = || {
        let (def, expr) = convection_diffusion(&complex, &geo, 0.7, [1.0, -0.5, 0.25]);
        let op = def.lower(expr);
        op.report().to_json()
    };
    let (r1, r2) = (build(), build());
    assert_eq!(r1, r2, "generation must be deterministic");
    println!("{r1}");
    assert!(
        r1.contains("\"provenance\":\"derived\""),
        "derived atoms marked"
    );
    assert!(r1.contains("\"name\":\"advection\",") && r1.contains("\"provenance\":\"hand\""));
    assert!(
        r1.contains("\"kernel\":\"scale_k\""),
        "tile-kernel metadata surfaced"
    );
    assert!(
        r1.contains("d0T(M1(d0(u)))"),
        "structural fingerprint present: {r1}"
    );
    log(
        "determinism",
        "pass",
        "reports identical across regeneration",
    );
}

#[test]
fn dwr_indicators_localize() {
    // Residual concentrated at one dof + uniform dual weight → the
    // indicator ranks that dof first.
    let r = vec![0.0, 0.1, -3.0, 0.2];
    let z = vec![1.0, 1.0, 0.5, 1.0];
    let eta = fs_opdsl::dwr_indicators(&r, &z);
    assert_eq!(eta.len(), 4);
    let argmax = eta
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("nonempty")
        .0;
    assert_eq!(argmax, 2);
    assert!((eta[2] - 1.5).abs() < 1e-15);
    log("dwr", "pass", "algebraic indicators |r_i z_i|");
}

#[test]
fn mms_poisson_through_generated_operator() {
    let pi = std::f64::consts::PI;
    let u_exact = move |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let f_exact = move |p: [f64; 3]| 3.0 * pi * pi * u_exact(p);
    let report = mms_poisson_study(
        &[4, 8, 16],
        |complex, geo| {
            let (def, expr) = poisson(complex, geo);
            let op = def.lower(expr);
            fs_opdsl::mms::materialize_dense(&op)
        },
        &u_exact,
        &f_exact,
    );
    for (n, e) in report.ns.iter().zip(&report.errors) {
        log("mms", "info", &format!("n={n} L2={e:.4e}"));
    }
    for o in &report.orders {
        assert!(
            (o - 2.0).abs() < 0.5,
            "MMS order {o:.2} (report {:?})",
            report.orders
        );
    }
    log("mms-order", "pass", &format!("orders {:?}", report.orders));
}

#[test]
fn elasticity_affine_patch_and_rigid_modes() {
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let (def, expr) = elasticity(&complex, &geo, 100.0, 0.25);
    assert_eq!(def.field_space.degree, VECTOR_DOFS);
    let op = def.lower(expr);
    let nv = positions.len();
    // (a) Rigid translation: zero strain → K u = 0 EXACTLY at interior
    // dofs (rows sum against constant fields; roundoff-level).
    let mut u_trans = vec![0.0f64; 3 * nv];
    for v in 0..nv {
        u_trans[3 * v] = 1.0;
        u_trans[3 * v + 1] = -2.0;
        u_trans[3 * v + 2] = 0.5;
    }
    let r = op.apply(&u_trans);
    let worst_t = r.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        worst_t < 1e-9,
        "translation must be a rigid mode: {worst_t:.3e}"
    );
    // (b) Infinitesimal rotation u = ω × x: antisymmetric gradient,
    // zero strain → K u = 0.
    let omega = [0.3f64, -0.7, 1.1];
    let mut u_rot = vec![0.0f64; 3 * nv];
    for (v, p) in positions.iter().enumerate() {
        u_rot[3 * v] = omega[1].mul_add(p[2], -(omega[2] * p[1]));
        u_rot[3 * v + 1] = omega[2].mul_add(p[0], -(omega[0] * p[2]));
        u_rot[3 * v + 2] = omega[0].mul_add(p[1], -(omega[1] * p[0]));
    }
    let r = op.apply(&u_rot);
    let worst_r = r.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        worst_r < 1e-9,
        "infinitesimal rotation must be a rigid mode: {worst_r:.3e}"
    );
    // (c) Uniform strain patch: u = B x (symmetric B) → constant
    // stress → interior residual rows vanish.
    let bmat: [[f64; 3]; 3] = [
        [1e-3, 2e-4, -1e-4],
        [2e-4, -5e-4, 3e-4],
        [-1e-4, 3e-4, 7e-4],
    ];
    let mut u_strain = vec![0.0f64; 3 * nv];
    for (v, p) in positions.iter().enumerate() {
        for i in 0..3 {
            u_strain[3 * v + i] =
                bmat[i][0].mul_add(p[0], bmat[i][1].mul_add(p[1], bmat[i][2] * p[2]));
        }
    }
    let r = op.apply(&u_strain);
    let mut worst_i = 0.0f64;
    for (v, p) in positions.iter().enumerate() {
        if !fs_feec::on_unit_cube_boundary(*p) {
            for i in 0..3 {
                worst_i = worst_i.max(r[3 * v + i].abs());
            }
        }
    }
    assert!(
        worst_i < 1e-12,
        "uniform-strain patch residual {worst_i:.3e}"
    );
    log(
        "elasticity",
        "pass",
        &format!("rigid {worst_t:.1e}/{worst_r:.1e}, patch {worst_i:.1e}"),
    );
}

#[test]
fn perf_generated_vs_hand_documented() {
    // MEASURED, not assumed: matrix-free (3 SpMV walk) and
    // materialized (single CSR) applies vs the hand-assembled
    // operator. Debug-build timings are logged for the record; only
    // sanity bounds are asserted (real numbers go in the bead close).
    let (complex, positions) = kuhn_cube(6);
    let geo = element_geometry(&complex, &positions);
    let (def, expr) = poisson(&complex, &geo);
    let op = def.lower(expr);
    let hand = fs_feec::stiffness(
        &fs_feec::incidence_to_csr(&complex.d0()),
        &fs_feec::mass_matrix(&complex, &geo, 1),
    );
    let materialized = op.materialize().expect("linear");
    let u = rand_vec(complex.vertex_count, 30);
    let reps = 50u32;
    let time = |f: &mut dyn FnMut()| -> f64 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            f();
        }
        t0.elapsed().as_secs_f64()
    };
    let mut sink = 0.0f64;
    let t_hand = time(&mut || {
        let mut y = vec![0.0f64; u.len()];
        hand.spmv(&u, &mut y);
        sink += y[0];
    });
    let t_mat = time(&mut || {
        let mut y = vec![0.0f64; u.len()];
        materialized.spmv(&u, &mut y);
        sink += y[0];
    });
    let t_free = time(&mut || {
        let y = op.apply(&u);
        sink += y[0];
    });
    assert!(sink.is_finite());
    let ratio_mat = t_hand / t_mat.max(1e-12);
    let ratio_free = t_hand / t_free.max(1e-12);
    log(
        "perf",
        "info",
        &format!(
            "hand={t_hand:.4}s materialized={t_mat:.4}s (x{ratio_mat:.2}) \
             matrix_free={t_free:.4}s (x{ratio_free:.2}) reps={reps}"
        ),
    );
    // Sanity only: materialized within 2x of hand (same structure);
    // matrix-free within 20x (three SpMVs + combinators, debug build).
    assert!(
        ratio_mat > 0.5,
        "materialized apply unexpectedly slow: x{ratio_mat:.2}"
    );
    assert!(
        ratio_free > 0.05,
        "matrix-free apply pathologically slow: x{ratio_free:.2}"
    );
}

const GOLDEN_HASH: u64 = 0x8b28_77cc_cb43_7cbc; // recorded at tfz.4 landing, frozen

#[test]
fn opdsl_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let n = complex.vertex_count;
    let u = rand_vec(n, 40);
    // Poisson apply + materialized sample.
    let (pd, pe) = poisson(&complex, &geo);
    let pop = pd.lower(pe);
    for v in pop.apply(&u).iter().step_by(3) {
        feed(*v);
    }
    // Convection–diffusion jvp + vjp.
    let (cd, ce) = convection_diffusion(&complex, &geo, 0.7, [1.0, -0.5, 0.25]);
    let cop = cd.lower(ce);
    for v in cop.jvp(&u).iter().step_by(3) {
        feed(*v);
    }
    for v in cop.vjp(&u).iter().step_by(3) {
        feed(*v);
    }
    // Reaction–diffusion linearized jvp.
    let (rd, re) = reaction_diffusion(&complex, &geo, 0.8);
    let mut rop = rd.lower(re);
    rop.linearize(&u);
    let v = rand_vec(n, 41);
    for x in rop.jvp(&v).iter().step_by(3) {
        feed(*x);
    }
    // Elasticity apply sample.
    let (ed, ee) = elasticity(&complex, &geo, 100.0, 0.25);
    let eop = ed.lower(ee);
    let uv = rand_vec(3 * n, 42);
    for x in eop.apply(&uv).iter().step_by(7) {
        feed(*x);
    }
    log("opdsl-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "opdsl bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
