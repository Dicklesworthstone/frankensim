//! fs-adjoint battery (tfz.24): the adjoint-vs-FD triangle on the
//! matrix-free IFT gradient (source AND density/SIMP parameters),
//! Sobolev smoothing demonstrably rescuing a mesh-noisy gradient
//! (cosine-similarity numbers, not vibes), Hadamard shape gradients
//! verified against perturb-and-resolve FD (volume: tight; Poisson
//! compliance: relative gate with the discretization error REPORTED),
//! revolve-checkpointed time-dependent adjoints (gradient w.r.t. the
//! initial condition vs FD, recompute counts logged), the
//! verification gate catching a corrupted gradient, and the golden.

use fs_adjoint::{
    DensityOp, DensityPoisson, GradientVerdict, HeatAdjoint, heat_initial_gradient,
    ift_gradient_matfree, sobolev_smooth, verify_gradient, volume_shape_gradient,
};
use fs_feec::{element_geometry, kuhn_cube};
use fs_rand::StreamKey;
use fs_solver::{CgState, CsrOp, LinearOp};
use fs_sparse::precond::IdentityPrecond;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn assert_gradient_refused(verdict: &GradientVerdict) {
    assert!(!verdict.pass, "invalid evidence passed the gradient gate");
    assert!(
        verdict.max_rel_err.is_infinite() && verdict.max_rel_err.is_sign_positive(),
        "refusal must use the deterministic +infinity sentinel: {verdict:?}"
    );
    assert!(
        verdict
            .pairs
            .iter()
            .all(|(analytic, fd)| analytic.is_finite() && fd.is_finite()),
        "refusal exposed non-finite directional evidence: {verdict:?}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 31,
        kernel: 0xAD10,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

/// Interior-reduced (M0, K0) pair on kuhn(m).
fn poisson_pair(m: usize) -> (fs_sparse::Csr, fs_sparse::Csr, usize) {
    let (complex, positions) = kuhn_cube(m);
    let geo = element_geometry(&complex, &positions);
    let k0 = fs_feec::stiffness(
        &fs_feec::incidence_to_csr(&complex.d0()),
        &fs_feec::mass_matrix(&complex, &geo, 1),
    );
    let m0 = fs_feec::mass_matrix(&complex, &geo, 0);
    let interior: Vec<usize> = (0..positions.len())
        .filter(|&v| !fs_feec::on_unit_cube_boundary(positions[v]))
        .collect();
    let mut slot = vec![usize::MAX; positions.len()];
    for (i, &v) in interior.iter().enumerate() {
        slot[v] = i;
    }
    let reduce = |a: &fs_sparse::Csr| -> fs_sparse::Csr {
        let mut red = fs_sparse::Coo::new(interior.len(), interior.len());
        for (i, &v) in interior.iter().enumerate() {
            let (cols, vals) = a.row(v);
            for (&c, &val) in cols.iter().zip(vals) {
                if slot[c] != usize::MAX {
                    red.push(i, slot[c], val);
                }
            }
        }
        red.assemble()
    };
    let n = interior.len();
    (reduce(&m0), reduce(&k0), n)
}

fn solve_spd(a: &CsrOp, b: &[f64]) -> Vec<f64> {
    let mut st = CgState::new(a, &IdentityPrecond, b);
    let rep = st.run(a, &IdentityPrecond, 1e-13, 20_000);
    assert!(rep.converged, "solve failed: {rep:?}");
    st.x
}

#[test]
fn ift_source_gradient_triangle() {
    // R(u; p) = K·u − M·p, J = ½(u−t)ᵀ(u−t). Adjoint gradient vs
    // central FD along random directions — the acceptance triangle
    // (the dual-number leg lives in fs-opdsl's battery; here the
    // solver is in the loop, which duals cannot reach).
    let (m0, k0, n) = poisson_pair(4);
    let k_op = CsrOp::symmetric(k0);
    let p0 = rand_vec(n, 1);
    let target = rand_vec(n, 2);
    let solve_u = |p: &[f64]| -> Vec<f64> {
        let mut b = vec![0.0f64; n];
        m0.spmv(p, &mut b);
        solve_spd(&k_op, &b)
    };
    let j = |p: &[f64]| -> f64 {
        let u = solve_u(p);
        0.5 * u
            .iter()
            .zip(&target)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
    };
    // Adjoint: Kᵀλ = (u − t); dJ/dp = +Mᵀλ (R = Ku − Mp ⇒ ∂R/∂p = −M).
    let u = solve_u(&p0);
    let djdu: Vec<f64> = u.iter().zip(&target).map(|(a, b)| a - b).collect();
    let (g, rep) = ift_gradient_matfree(
        &k_op,
        &djdu,
        &[],
        &|lam| {
            let mut out = vec![0.0f64; n];
            m0.spmv(lam, &mut out); // M symmetric
            for v in &mut out {
                *v = -*v;
            }
            out
        },
        1e-12,
        200,
    );
    assert!(rep.converged, "adjoint solve failed: {rep:?}");
    let dirs: Vec<Vec<f64>> = (0..3).map(|k| rand_vec(n, 10 + k)).collect();
    let verdict = verify_gradient(&j, &p0, &g, &dirs, 1e-5, 1e-6);
    assert!(
        verdict.pass,
        "IFT source gradient failed FD check: {:?}",
        verdict.max_rel_err
    );
    log(
        "ift-source",
        "pass",
        &format!(
            "rel_err={:.2e} adjoint_iters={}",
            verdict.max_rel_err, rep.iters
        ),
    );
}

#[test]
fn density_chain_rule_simp() {
    // K(ρ)u = b with per-cell densities: dJ/dρ_t = −λᵀK_t u — the
    // exact volumetric (SIMP) chain rule, FD-verified.
    let (complex, positions) = kuhn_cube(2);
    let nt = complex.tets.len();
    let rho0: Vec<f64> = rand_vec(nt, 20).iter().map(|v| 1.0 + 0.3 * v).collect();
    let problem = DensityPoisson::new(&complex, &positions, rho0.clone());
    let n = problem.n();
    let b = rand_vec(n, 21);
    let target = rand_vec(n, 22);
    let j = |rho: &[f64]| -> f64 {
        let pr = DensityPoisson::new(&complex, &positions, rho.to_vec());
        let op = DensityOp::new(&pr);
        let u = solve_spd_op(&op, &b);
        0.5 * u
            .iter()
            .zip(&target)
            .map(|(a, c)| (a - c) * (a - c))
            .sum::<f64>()
    };
    // Adjoint at rho0.
    let op = DensityOp::new(&problem);
    let u = solve_spd_op(&op, &b);
    let djdu: Vec<f64> = u.iter().zip(&target).map(|(a, c)| a - c).collect();
    let lambda = solve_spd_op(&op, &djdu); // symmetric K
    let mut g = problem.density_pullback(&lambda, &u);
    for v in &mut g {
        *v = -*v;
    }
    let dirs: Vec<Vec<f64>> = (0..3).map(|k| rand_vec(nt, 30 + k)).collect();
    let verdict = verify_gradient(&j, &rho0, &g, &dirs, 1e-6, 1e-5);
    assert!(
        verdict.pass,
        "SIMP density gradient failed FD check: {:?}",
        verdict.max_rel_err
    );
    log(
        "density-simp",
        "pass",
        &format!("rel_err={:.2e}", verdict.max_rel_err),
    );
}

fn solve_spd_op<A: LinearOp>(a: &A, b: &[f64]) -> Vec<f64> {
    let mut st = CgState::new(a, &IdentityPrecond, b);
    let rep = st.run(a, &IdentityPrecond, 1e-13, 20_000);
    assert!(rep.converged, "solve failed: {rep:?}");
    st.x
}

#[test]
fn sobolev_smoothing_rescues_noisy_gradient() {
    // A smooth "true" gradient plus grid-frequency noise: the H¹
    // Riesz representation must recover direction alignment the raw
    // vector has lost. Numbers, not vibes.
    let (m0, k0, n) = poisson_pair(6);
    let (complex, positions) = kuhn_cube(6);
    let _ = &complex;
    let interior: Vec<[f64; 3]> = positions
        .iter()
        .copied()
        .filter(|&p| !fs_feec::on_unit_cube_boundary(p))
        .collect();
    assert_eq!(interior.len(), n);
    let g_true: Vec<f64> = interior
        .iter()
        .map(|p| (2.0 * p[0] - 1.0) * (1.5 - p[1]) + 0.5 * p[2])
        .collect();
    // Grid-frequency noise: alternating-sign per index, same scale as
    // the signal.
    let noise: Vec<f64> = (0..n)
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    let g_raw: Vec<f64> = g_true.iter().zip(&noise).map(|(t, e)| t + e).collect();
    let cosine = |a: &[f64], b: &[f64]| -> f64 {
        let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
        let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
        let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
        dot / (na * nb)
    };
    let raw_cos = cosine(&g_raw, &g_true);
    // α = h² (the standard scaling). Measured on this fixture:
    // raw 0.49 → smoothed 0.93 (larger α clamps the signal toward
    // the zero-Dirichlet boundary of the interior-reduced K and
    // HURTS alignment — 0.89 at α = 0.2; the metric choice is a real
    // tradeoff, which is the point of making it configurable).
    let h = 1.0 / 6.0;
    let (g_smooth, iters) = sobolev_smooth(&m0, &k0, h * h, &g_raw, 1e-12);
    let smooth_cos = cosine(&g_smooth, &g_true);
    assert!(
        smooth_cos > 0.9,
        "smoothed gradient poorly aligned: cos={smooth_cos:.4} (raw {raw_cos:.4})"
    );
    assert!(
        smooth_cos > raw_cos + 0.3,
        "smoothing must materially improve alignment: raw={raw_cos:.4} smooth={smooth_cos:.4}"
    );
    log(
        "sobolev",
        "pass",
        &format!("raw_cos={raw_cos:.4} smooth_cos={smooth_cos:.4} iters={iters}"),
    );
}

#[test]
fn hadamard_volume_matches_fd_exactly() {
    // dVol[V] via the boundary integral vs central FD of the discrete
    // mesh volume under vertex perturbation x ← x + ε·V(x): both are
    // polynomial in ε, so agreement is tight.
    let (complex, positions) = kuhn_cube(2);
    // NONZERO divergence on purpose: a divergence-free field has
    // dVol = 0 exactly, which would make the relative comparison
    // meaningless (the first draft hit exactly that: two zeros).
    let velocity = |p: [f64; 3]| -> [f64; 3] {
        [
            0.3f64.mul_add(p[0], 0.2 * p[1]),
            0.1f64.mul_add(p[1], -0.4 * p[2]),
            0.5f64.mul_add(p[2], 0.25 * p[0]),
        ]
    };
    let analytic = volume_shape_gradient(&complex, &positions, &velocity);
    let vol = |pos: &[[f64; 3]]| -> f64 {
        let geo = element_geometry(&complex, pos);
        geo.vol_signed.iter().map(|v| v.abs()).sum()
    };
    let eps = 1e-6;
    let perturb = |sign: f64| -> Vec<[f64; 3]> {
        positions
            .iter()
            .map(|&p| {
                // Only boundary vertices move (interior motion does
                // not change the volume; keeping it fixed matches the
                // boundary-integral formula's domain).
                if fs_feec::on_unit_cube_boundary(p) {
                    let v = velocity(p);
                    [
                        (sign * eps).mul_add(v[0], p[0]),
                        (sign * eps).mul_add(v[1], p[1]),
                        (sign * eps).mul_add(v[2], p[2]),
                    ]
                } else {
                    p
                }
            })
            .collect()
    };
    let fd = (vol(&perturb(1.0)) - vol(&perturb(-1.0))) / (2.0 * eps);
    let rel = (analytic - fd).abs() / fd.abs().max(1e-30);
    assert!(
        rel < 1e-7,
        "volume Hadamard vs FD: {analytic:.10e} vs {fd:.10e} (rel {rel:.2e})"
    );
    log("hadamard-volume", "pass", &format!("rel={rel:.2e}"));
}

#[test]
fn heat_adjoint_gradient_vs_fd() {
    // Terminal-misfit gradient w.r.t. the initial condition through
    // the revolve-checkpointed reverse sweep; every reverse step is a
    // transposed solve. FD-verified; recompute count logged.
    let (m0, k0, n) = poisson_pair(3);
    let heat = HeatAdjoint::new(m0, &k0, 0.01, 12);
    let u0 = rand_vec(n, 40);
    let target = rand_vec(n, 41);
    let (g, fwd_evals) = heat_initial_gradient(&heat, &u0, &target);
    let j = |u0v: &[f64]| -> f64 {
        let un = heat.forward(u0v);
        0.5 * un
            .iter()
            .zip(&target)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
    };
    let dirs: Vec<Vec<f64>> = (0..2).map(|k| rand_vec(n, 50 + k)).collect();
    let verdict = verify_gradient(&j, &u0, &g, &dirs, 1e-5, 1e-6);
    assert!(
        verdict.pass,
        "heat adjoint failed FD check: {:?}",
        verdict.max_rel_err
    );
    log(
        "heat-adjoint",
        "pass",
        &format!(
            "rel_err={:.2e} forward_evals={fwd_evals} (steps=12, O(log N) memory)",
            verdict.max_rel_err
        ),
    );
}

#[test]
fn verification_gate_catches_corruption() {
    // The gate must FAIL a corrupted gradient — a gate that cannot
    // fail is not a gate.
    let j = |p: &[f64]| -> f64 { p.iter().map(|x| x * x * x).sum() };
    let p = rand_vec(20, 60);
    let good: Vec<f64> = p.iter().map(|x| 3.0 * x * x).collect();
    let mut bad = good.clone();
    bad[7] *= 1.5;
    let dirs: Vec<Vec<f64>> = (0..4).map(|k| rand_vec(20, 70 + k)).collect();
    let ok = verify_gradient(&j, &p, &good, &dirs, 1e-6, 1e-8);
    assert!(ok.pass, "correct gradient rejected: {:?}", ok.max_rel_err);
    let caught = verify_gradient(&j, &p, &bad, &dirs, 1e-6, 1e-8);
    assert!(!caught.pass, "corrupted gradient passed the gate");
    log("verify-gate", "pass", "accepts correct, rejects corrupted");
}

#[test]
fn verification_gate_fails_closed_on_non_finite_inputs() {
    let objective = |point: &[f64]| point.iter().map(|value| value * value).sum();
    let point = vec![1.0, -2.0];
    let gradient = vec![2.0, -4.0];
    let directions = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

    for non_finite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let invalid_objective = |_: &[f64]| non_finite;
        assert_gradient_refused(&verify_gradient(
            &invalid_objective,
            &point,
            &gradient,
            &directions,
            1e-6,
            1e-8,
        ));

        let mut invalid_point = point.clone();
        invalid_point[0] = non_finite;
        assert_gradient_refused(&verify_gradient(
            &objective,
            &invalid_point,
            &gradient,
            &directions,
            1e-6,
            1e-8,
        ));

        let mut invalid_gradient = gradient.clone();
        invalid_gradient[0] = non_finite;
        assert_gradient_refused(&verify_gradient(
            &objective,
            &point,
            &invalid_gradient,
            &directions,
            1e-6,
            1e-8,
        ));

        let mut invalid_directions = directions.clone();
        invalid_directions[0][0] = non_finite;
        assert_gradient_refused(&verify_gradient(
            &objective,
            &point,
            &gradient,
            &invalid_directions,
            1e-6,
            1e-8,
        ));
    }

    for invalid_eps in [0.0, -1e-6, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_gradient_refused(&verify_gradient(
            &objective,
            &point,
            &gradient,
            &directions,
            invalid_eps,
            1e-8,
        ));
    }
    for invalid_tol in [0.0, -1e-8, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_gradient_refused(&verify_gradient(
            &objective,
            &point,
            &gradient,
            &directions,
            1e-6,
            invalid_tol,
        ));
    }
}

#[test]
fn verification_gate_fails_closed_on_non_finite_intermediates() {
    let finite_objective = |point: &[f64]| point.iter().copied().sum();

    let analytic_overflow = verify_gradient(
        &finite_objective,
        &[0.0, 0.0],
        &[f64::MAX, f64::MAX],
        &[vec![1.0, 1.0]],
        1.0,
        1e-8,
    );
    assert_gradient_refused(&analytic_overflow);

    let step_overflow = verify_gradient(
        &finite_objective,
        &[0.0],
        &[0.0],
        &[vec![f64::MAX]],
        2.0,
        1e-8,
    );
    assert_gradient_refused(&step_overflow);

    let coordinate_overflow = verify_gradient(
        &finite_objective,
        &[f64::MAX],
        &[0.0],
        &[vec![1.0]],
        f64::MAX / 2.0,
        1e-8,
    );
    assert_gradient_refused(&coordinate_overflow);

    let subtraction_overflow = |point: &[f64]| {
        if point[0].is_sign_positive() {
            f64::MAX
        } else {
            -f64::MAX
        }
    };
    assert_gradient_refused(&verify_gradient(
        &subtraction_overflow,
        &[0.0],
        &[0.0],
        &[vec![1.0]],
        1.0,
        1e-8,
    ));

    assert_gradient_refused(&verify_gradient(
        &finite_objective,
        &[0.0],
        &[0.0],
        &[vec![f64::MIN_POSITIVE]],
        f64::MAX,
        1e-8,
    ));

    let division_overflow = |point: &[f64]| {
        if point[0].is_sign_positive() {
            1.0
        } else {
            0.0
        }
    };
    assert_gradient_refused(&verify_gradient(
        &division_overflow,
        &[0.0],
        &[0.0],
        &[vec![1.0]],
        f64::from_bits(1),
        1e-8,
    ));

    let relative_error_overflow = |point: &[f64]| -0.5 * f64::MAX * point[0];
    assert_gradient_refused(&verify_gradient(
        &relative_error_overflow,
        &[0.0],
        &[f64::MAX],
        &[vec![1.0]],
        1.0,
        1e-8,
    ));
}

#[test]
fn verification_gate_rejects_vacuous_or_unrepresentable_directions() {
    let objective = |point: &[f64]| point.iter().map(|value| value * value).sum();
    assert_gradient_refused(&verify_gradient(
        &objective,
        &[1.0],
        &[2.0],
        &[],
        1e-6,
        1e-8,
    ));
    assert_gradient_refused(&verify_gradient(
        &objective,
        &[1.0],
        &[2.0],
        &[vec![0.0]],
        1e-6,
        1e-8,
    ));
    assert_gradient_refused(&verify_gradient(
        &objective,
        &[f64::MAX],
        &[0.0],
        &[vec![f64::MIN_POSITIVE]],
        1.0,
        1e-8,
    ));
}

#[test]
fn verification_gate_is_invariant_to_paired_direction_step_rescaling() {
    let objective = |point: &[f64]| 2.0 * point[0] - 0.5 * point[1];
    let point = [0.0, 0.0];
    let gradient = [2.0, -0.5];
    let direction = vec![0.25, -0.5];
    let scaled_direction = direction.iter().map(|value| 8.0 * value).collect();
    let base = verify_gradient(
        &objective,
        &point,
        &gradient,
        &[direction],
        2.0_f64.powi(-10),
        1e-12,
    );
    let scaled = verify_gradient(
        &objective,
        &point,
        &gradient,
        &[scaled_direction],
        2.0_f64.powi(-13),
        1e-12,
    );

    assert!(base.pass && scaled.pass);
    assert_eq!(base.max_rel_err.to_bits(), scaled.max_rel_err.to_bits());
    assert_eq!(
        scaled.pairs[0].0.to_bits(),
        (8.0 * base.pairs[0].0).to_bits()
    );
    assert_eq!(
        scaled.pairs[0].1.to_bits(),
        (8.0 * base.pairs[0].1).to_bits()
    );
}

#[test]
fn verification_gate_tiny_direction_cannot_hide_a_wrong_gradient() {
    let verdict = verify_gradient(
        &|point: &[f64]| point[0],
        &[0.0],
        &[0.0],
        &[vec![f64::MIN_POSITIVE]],
        1.0,
        1e-8,
    );

    assert!(
        !verdict.pass,
        "tiny direction suppressed a unit gradient error"
    );
    assert_eq!(verdict.pairs, vec![(0.0, f64::MIN_POSITIVE)]);
    assert_eq!(verdict.max_rel_err.to_bits(), 1.0_f64.to_bits());
}

#[test]
fn verification_gate_retains_only_a_deterministic_finite_prefix() {
    let objective = |point: &[f64]| {
        if point[1] == 0.0 {
            point[0] * point[0]
        } else {
            f64::NAN
        }
    };
    let directions = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let first = verify_gradient(
        &objective,
        &[0.0, 0.0],
        &[0.0, 0.0],
        &directions,
        1e-6,
        1e-8,
    );
    let second = verify_gradient(
        &objective,
        &[0.0, 0.0],
        &[0.0, 0.0],
        &directions,
        1e-6,
        1e-8,
    );

    assert_gradient_refused(&first);
    assert_gradient_refused(&second);
    assert_eq!(first.pairs, vec![(0.0, 0.0)]);
    assert_eq!(first.pairs, second.pairs);
}

const GOLDEN_HASH: u64 = 0x0896_7e37_81b3_c044; // recorded at tfz.24 landing, frozen

#[test]
fn adjoint_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // IFT source gradient sample.
    let (m0, k0, n) = poisson_pair(3);
    let k_op = CsrOp::symmetric(k0.clone());
    let p0 = rand_vec(n, 80);
    let mut b = vec![0.0f64; n];
    m0.spmv(&p0, &mut b);
    let u = solve_spd(&k_op, &b);
    let djdu = u.clone();
    let (g, _) = ift_gradient_matfree(
        &k_op,
        &djdu,
        &[],
        &|lam| {
            let mut out = vec![0.0f64; n];
            m0.spmv(lam, &mut out);
            for v in &mut out {
                *v = -*v;
            }
            out
        },
        1e-11,
        200,
    );
    for v in g.iter().step_by(3) {
        feed(*v);
    }
    // Sobolev smoothing sample.
    let g_raw = rand_vec(n, 81);
    let (gs, _) = sobolev_smooth(&m0, &k0, 0.05, &g_raw, 1e-11);
    for v in gs.iter().step_by(5) {
        feed(*v);
    }
    // Heat adjoint sample.
    let heat = HeatAdjoint::new(m0, &k0, 0.02, 8);
    let u0 = rand_vec(n, 82);
    let target = rand_vec(n, 83);
    let (gh, _) = heat_initial_gradient(&heat, &u0, &target);
    for v in gh.iter().step_by(2) {
        feed(*v);
    }
    log("adjoint-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "adjoint bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn hadamard_compliance_agreement_improves_under_refinement() {
    // Compliance J = ∫ f·u (Dirichlet Poisson): the Hadamard boundary
    // form −∫(∂u/∂n)²·(V·n) dA on P1 solutions carries discretization
    // error, so the honest gate is not a tight tolerance but
    // CONSISTENCY: the relative gap to perturb-and-resolve FD must
    // SHRINK under mesh refinement, and be modest on the fine mesh.
    let velocity = |p: [f64; 3]| -> [f64; 3] {
        [
            0.2f64.mul_add(p[0], 0.1),
            0.15f64.mul_add(p[1], -0.05),
            0.1f64.mul_add(p[2], 0.05),
        ]
    };
    let f_src = |p: [f64; 3]| -> f64 { 1.0 + p[0] + 0.5 * p[1] };
    let gap_at = |m: usize| -> f64 {
        let (complex, positions) = kuhn_cube(m);
        let compliance = |pos: &[[f64; 3]]| -> (f64, Vec<f64>) {
            let geo = element_geometry(&complex, pos);
            let k0 = fs_feec::stiffness(
                &fs_feec::incidence_to_csr(&complex.d0()),
                &fs_feec::mass_matrix(&complex, &geo, 1),
            );
            let m0 = fs_feec::mass_matrix(&complex, &geo, 0);
            let interior: Vec<usize> = (0..pos.len())
                .filter(|&v| !fs_feec::on_unit_cube_boundary(positions[v]))
                .collect();
            let mut slot = vec![usize::MAX; pos.len()];
            for (i, &v) in interior.iter().enumerate() {
                slot[v] = i;
            }
            let mut red = fs_sparse::Coo::new(interior.len(), interior.len());
            for (i, &v) in interior.iter().enumerate() {
                let (cols, vals) = k0.row(v);
                for (&c, &val) in cols.iter().zip(vals) {
                    if slot[c] != usize::MAX {
                        red.push(i, slot[c], val);
                    }
                }
            }
            let a = CsrOp::symmetric(red.assemble());
            let fvals: Vec<f64> = pos.iter().map(|&p| f_src(p)).collect();
            let mut bfull = vec![0.0f64; pos.len()];
            m0.spmv(&fvals, &mut bfull);
            let b: Vec<f64> = interior.iter().map(|&v| bfull[v]).collect();
            let x = solve_spd(&a, &b);
            // J = Σ interior b_i x_i (boundary u = 0).
            let j: f64 = b.iter().zip(&x).map(|(bi, xi)| bi * xi).sum();
            let mut ufull = vec![0.0f64; pos.len()];
            for (i, &v) in interior.iter().enumerate() {
                ufull[v] = x[i];
            }
            (j, ufull)
        };
        let (_, u_full) = compliance(&positions);
        let geo = element_geometry(&complex, &positions);
        let analytic =
            fs_adjoint::compliance_shape_gradient(&complex, &positions, &geo, &u_full, &velocity);
        let eps = 1e-4;
        let perturb = |sign: f64| -> Vec<[f64; 3]> {
            positions
                .iter()
                .map(|&p| {
                    if fs_feec::on_unit_cube_boundary(p) {
                        let v = velocity(p);
                        [
                            (sign * eps).mul_add(v[0], p[0]),
                            (sign * eps).mul_add(v[1], p[1]),
                            (sign * eps).mul_add(v[2], p[2]),
                        ]
                    } else {
                        p
                    }
                })
                .collect()
        };
        let (jp, _) = compliance(&perturb(1.0));
        let (jm, _) = compliance(&perturb(-1.0));
        let fd = (jp - jm) / (2.0 * eps);
        (analytic - fd).abs() / fd.abs().max(1e-30)
    };
    let gap_coarse = gap_at(2);
    let gap_fine = gap_at(4);
    log(
        "hadamard-compliance",
        "info",
        &format!("rel gap: kuhn(2)={gap_coarse:.3} kuhn(4)={gap_fine:.3}"),
    );
    assert!(
        gap_fine < gap_coarse,
        "Hadamard/FD agreement must improve under refinement: {gap_coarse:.3} -> {gap_fine:.3}"
    );
    // The P1 one-sided normal trace squared is LOW-ORDER accurate, so
    // the boundary form converges slowly (measured 0.84 → 0.66 across
    // one refinement; the sign and trend are the verified claims). At
    // lowest order the exactly-FD-verified VOLUMETRIC form (the SIMP
    // density pullback above, rel err ~1e-6) is the production path;
    // the boundary form earns tight tolerances with high-order traces
    // (recorded in the CONTRACT).
    assert!(gap_fine < 0.75, "fine-mesh gap too large: {gap_fine:.3}");
    log(
        "hadamard-compliance",
        "pass",
        &format!("consistent, fine gap {gap_fine:.3}"),
    );
}

#[test]
fn exact_hvp_symmetry_fd_gap_and_tr_win() {
    use fs_adjoint::density_misfit_hvp;
    // Fixture: kuhn(3) density Poisson misfit (8 interior dofs, 162
    // cells... enough structure for a meaningful Hessian).
    let (complex, positions) = kuhn_cube(3);
    let nt = complex.tets.len();
    let rho0: Vec<f64> = rand_vec(nt, 90).iter().map(|v| 1.2 + 0.3 * v).collect();
    let problem = DensityPoisson::new(&complex, &positions, rho0.clone());
    let n = problem.n();
    let b = rand_vec(n, 91);
    let u_star = rand_vec(n, 92);
    // (1) SYMMETRY: vᵀH w == wᵀH v (the IFT Hessian of a scalar
    // objective is symmetric — a strong exactness gate).
    let v = rand_vec(nt, 93);
    let w = rand_vec(nt, 94);
    let hv = density_misfit_hvp(&problem, &b, &u_star, &v, 1e-13);
    let hw = density_misfit_hvp(&problem, &b, &u_star, &w, 1e-13);
    let vthw: f64 = v.iter().zip(&hw).map(|(a, c)| a * c).sum();
    let wthv: f64 = w.iter().zip(&hv).map(|(a, c)| a * c).sum();
    let sym_rel = (vthw - wthv).abs() / vthw.abs().max(1e-30);
    assert!(
        sym_rel < 1e-8,
        "Hessian not symmetric: {vthw:.10e} vs {wthv:.10e}"
    );
    // (2) FD-of-gradients agreement + the FD accuracy gap QUANTIFIED.
    let grad_at = |rho: &[f64]| -> Vec<f64> {
        let pr = DensityPoisson::new(&complex, &positions, rho.to_vec());
        let op = DensityOp::new(&pr);
        let u = solve_spd_op(&op, &b);
        let djdu: Vec<f64> = u.iter().zip(&u_star).map(|(a, c)| a - c).collect();
        let lambda = solve_spd_op(&op, &djdu);
        let mut g = pr.density_pullback(&lambda, &u);
        for x in &mut g {
            *x = -*x;
        }
        g
    };
    let mut fd_err_best = f64::INFINITY;
    let mut fd_err_worst = 0.0f64;
    for eps in [1e-4f64, 1e-6, 1e-8] {
        let rp: Vec<f64> = rho0
            .iter()
            .zip(&v)
            .map(|(r, vi)| eps.mul_add(*vi, *r))
            .collect();
        let rm: Vec<f64> = rho0
            .iter()
            .zip(&v)
            .map(|(r, vi)| eps.mul_add(-vi, *r))
            .collect();
        let gp = grad_at(&rp);
        let gm = grad_at(&rm);
        let fd: Vec<f64> = gp
            .iter()
            .zip(&gm)
            .map(|(a, c)| (a - c) / (2.0 * eps))
            .collect();
        let num: f64 = fd
            .iter()
            .zip(&hv)
            .map(|(a, c)| (a - c) * (a - c))
            .sum::<f64>()
            .sqrt();
        let den: f64 = hv.iter().map(|x| x * x).sum::<f64>().sqrt();
        let rel = num / den;
        fd_err_best = fd_err_best.min(rel);
        fd_err_worst = fd_err_worst.max(rel);
    }
    assert!(
        fd_err_best < 1e-5,
        "exact Hv disagrees with best-eps FD: {fd_err_best:.2e}"
    );
    assert!(
        fd_err_worst > 10.0 * fd_err_best,
        "the FD accuracy gap should be visible across eps: {fd_err_best:.2e}..{fd_err_worst:.2e}"
    );
    log(
        "exact-hvp",
        "pass",
        &format!("sym {sym_rel:.1e}, FD gap {fd_err_best:.1e}..{fd_err_worst:.1e}"),
    );
}

// ---- Claim-integrity regressions (E02 sweep) --------------------------

#[test]
fn verification_gate_refuses_bit_insensitive_directions() {
    // frankensim-extreal-program-f85xj.2.20. The silent-zero seam: a
    // forward path that drops the parameter produces BOTH a zero
    // cotangent and two bit-identical objective evaluations, so
    // analytic = fd = 0 and the relative error is 0/scale = 0 — the
    // PERFECT score. Before the fix this returned pass = true with
    // max_rel_err = 0.0; the refusal preflight did not fire because the
    // direction is finite, non-zero, and both perturbations move the
    // stored coordinates.
    let bit_insensitive = |_point: &[f64]| 7.0;
    let verdict = verify_gradient(
        &bit_insensitive,
        &[1.0, 2.0],
        &[0.0, 0.0],
        &[vec![1.0, -1.0]],
        1e-6,
        1e-6,
    );
    assert_eq!(
        verdict.pairs,
        vec![(0.0, 0.0)],
        "the no-signal pair must still be reported as diagnostics"
    );
    assert_eq!(
        verdict.informative_directions, 0,
        "a zero-vs-zero probe carries no signal"
    );
    assert_gradient_refused(&verdict);

    // The gate must not over-refuse: one no-signal direction alongside a
    // real one still passes, and the count says how much evidence there
    // actually was.
    let linear_in_first = |point: &[f64]| point[0];
    let mixed = verify_gradient(
        &linear_in_first,
        &[1.0, 2.0],
        &[1.0, 0.0],
        &[vec![0.0, 1.0], vec![1.0, 0.0]],
        1e-6,
        1e-6,
    );
    assert!(mixed.pass, "a genuinely corroborated gradient still passes");
    assert_eq!(
        mixed.informative_directions, 1,
        "only one of the two probes carried signal"
    );

    // And the silent-zero gradient is still CAUGHT on the informative
    // direction rather than being excused by the no-signal one.
    let zeroed = verify_gradient(
        &linear_in_first,
        &[1.0, 2.0],
        &[0.0, 0.0],
        &[vec![0.0, 1.0], vec![1.0, 0.0]],
        1e-6,
        1e-6,
    );
    assert!(!zeroed.pass, "a silent-zero gradient must not pass");
    assert_eq!(zeroed.informative_directions, 1);
}

#[test]
#[should_panic(expected = "parameter-space pullback length")]
fn ift_gradient_refuses_a_truncated_pullback() {
    // frankensim-extreal-program-f85xj.2.26. A masked/sliced pullback
    // shorter than ∂J/∂p used to be zipped in silently: the returned
    // vector had the expected length, its head carried the total
    // derivative and its tail only the explicit partial, and the report
    // still said converged.
    let (_, k0, n) = poisson_pair(2);
    let k_op = CsrOp::symmetric(k0);
    let djdu = rand_vec(n, 70);
    let djdp = vec![1.0f64; 4];
    let _ = ift_gradient_matfree(
        &k_op,
        &djdu,
        &djdp,
        &|_lambda| vec![0.5f64; 2], // a 2-component pullback against a 4-component ∂J/∂p
        1e-10,
        50,
    );
}

#[test]
#[should_panic(expected = "Hessian solve tolerance")]
fn hessian_refuses_a_vacuous_solve_tolerance() {
    // frankensim-extreal-program-f85xj.2.27. CG checks rel_residual < tol
    // at the top of the loop with x = 0 and rel_residual = 1, so any
    // tol > 1 breaks out at zero iterations and reports converged. All
    // four inner solves then returned the zero vector and the function
    // published exact zeros as the EXACT H·v, behind an
    // assert!(rep.converged) that passed vacuously.
    let (complex, positions) = kuhn_cube(2);
    let nt = complex.tets.len();
    let problem = DensityPoisson::new(&complex, &positions, vec![1.0; nt]);
    let n = problem.n();
    let b = rand_vec(n, 95);
    let u_star = rand_vec(n, 96);
    let v = rand_vec(nt, 97);
    let _ = fs_adjoint::density_misfit_hvp(&problem, &b, &u_star, &v, 2.0);
}
