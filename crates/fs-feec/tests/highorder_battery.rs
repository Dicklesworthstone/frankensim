//! High-order battery, slice 1 (tfz.6): Gauss–Legendre exactness,
//! Lobatto basis structure, sum-factorized apply vs the assembled
//! Kronecker reference to roundoff (the acceptance gate), Jacobi
//! diagonal exactness, G1 MMS spectral/order convergence for r = 1..6
//! through the matrix-free PCG path, and the golden hash.

use fs_feec::highorder::hex::{TensorSpace, pcg_matfree};
use fs_feec::highorder::quad1d::{gauss_legendre, legendre, lobatto_shapes};
use fs_math::det;
use fs_rand::StreamKey;

const SUITE: &str = "fs-feec/highorder";
const FIXED_INPUT_SEED: u64 = 0;
const STREAM_INPUT_SEED: u64 = 6;

fn verdict(case: &str, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass: true,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("high-order verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("high-order verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn measurement(identity: &str, name: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, identity);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: name.to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("high-order measurement must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("high-order measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn finite_json(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "null".to_string()
    }
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: STREAM_INPUT_SEED,
        kernel: 0x460,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

#[test]
fn gauss_legendre_exactness_and_symmetry() {
    for n in 1..=10usize {
        let (x, w) = gauss_legendre(n);
        // Weights sum to 2, nodes symmetric.
        let wsum: f64 = w.iter().sum();
        assert!((wsum - 2.0).abs() < 1e-14, "n={n}: weight sum {wsum}");
        for i in 0..n {
            assert!(
                (x[i] + x[n - 1 - i]).abs() < 1e-14,
                "n={n}: nodes not symmetric"
            );
        }
        // Exact for monomials up to degree 2n−1.
        for d in 0..2 * n {
            let quad: f64 = x
                .iter()
                .zip(&w)
                .map(|(&xi, &wi)| wi * det::powi(xi, i32::try_from(d).expect("small")))
                .sum();
            let exact = if d % 2 == 1 {
                0.0
            } else {
                2.0 / (d as f64 + 1.0)
            };
            assert!(
                (quad - exact).abs() < 1e-13,
                "n={n} degree {d}: {quad} vs {exact}"
            );
        }
    }
    verdict(
        "gauss-legendre",
        "n=1..10 exact to degree 2n-1",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn lobatto_bubbles_vanish_at_endpoints() {
    for r in 1..=6usize {
        for &x in &[-1.0f64, 1.0] {
            let (vals, _) = lobatto_shapes(r, x);
            // Vertex functions partition the endpoint values.
            assert!((vals[0] - f64::from(u8::from(x < 0.0))).abs() < 1e-15);
            assert!((vals[1] - f64::from(u8::from(x > 0.0))).abs() < 1e-15);
            for (k, v) in vals.iter().enumerate().skip(2) {
                assert!(
                    v.abs() < 1e-14,
                    "r={r}: bubble {k} nonzero at endpoint {x}: {v}"
                );
            }
        }
        // Legendre orthogonality spot check keeps the recurrence honest.
        let (qx, qw) = gauss_legendre(r + 2);
        let mut inner = 0.0f64;
        for (&x, &w) in qx.iter().zip(&qw) {
            let (la, _) = legendre(r, x);
            let (lb, _) = legendre(r.saturating_sub(1).max(1), x);
            if r != r.saturating_sub(1).max(1) {
                inner += w * la * lb;
            }
        }
        assert!(inner.abs() < 1e-13, "r={r}: Legendre orthogonality {inner}");
    }
    verdict("lobatto", "endpoint structure r=1..6", FIXED_INPUT_SEED);
}

#[test]
fn sum_factorized_matches_assembled_kronecker() {
    // The acceptance gate: sum-factorized apply == assembled operator
    // to roundoff. Reference: K1⊗M1⊗M1 + M1⊗K1⊗M1 + M1⊗M1⊗K1 over
    // assembled dense 1D operators (exactly the Galerkin operator on
    // a tensor grid).
    for &(m, r) in &[(2usize, 1usize), (2, 3), (1, 5), (3, 2)] {
        let sp = TensorSpace::new(m, r);
        let n1 = sp.n1();
        let (m1, k1) = sp.assembled_1d();
        let stream_tile = 30 + u32::try_from(10 * m + r).expect("small");
        let u = rand_vec(sp.ndof(), stream_tile);
        let y_fast = sp.apply_stiffness(&u);
        // Dense Kronecker apply (small fixtures only).
        let mut y_ref = vec![0.0f64; sp.ndof()];
        for i in 0..n1 {
            for j in 0..n1 {
                for k in 0..n1 {
                    let mut acc = 0.0f64;
                    for a in 0..n1 {
                        for b in 0..n1 {
                            for c in 0..n1 {
                                let kk = k1[i * n1 + a] * m1[j * n1 + b] * m1[k * n1 + c]
                                    + m1[i * n1 + a] * k1[j * n1 + b] * m1[k * n1 + c]
                                    + m1[i * n1 + a] * m1[j * n1 + b] * k1[k * n1 + c];
                                if kk != 0.0 {
                                    acc = kk.mul_add(u[sp.gid(a, b, c)], acc);
                                }
                            }
                        }
                    }
                    y_ref[sp.gid(i, j, k)] = acc;
                }
            }
        }
        let scale = y_ref
            .iter()
            .map(|v| v.abs())
            .fold(0.0f64, f64::max)
            .max(1.0);
        let worst = y_fast
            .iter()
            .zip(&y_ref)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(
            worst < 1e-12 * scale,
            "m={m} r={r}: sum-factorized deviates {worst:.3e} (scale {scale:.3e})"
        );
        let case = format!("sumfact-vs-assembled/m-{m}/r-{r}");
        measurement(
            &format!("{case}/measurement"),
            "sumfact-vs-assembled",
            format!(
                "{{\"m\":{m},\"r\":{r},\"deviation\":{},\"scale\":{},\
                 \"input_seed\":{STREAM_INPUT_SEED},\"stream_kernel\":\"0x460\",\
                 \"stream_tile\":{stream_tile}}}",
                finite_json(worst),
                finite_json(scale)
            ),
        );
        verdict(
            &case,
            &format!("m={m} r={r} dev={worst:.2e}"),
            STREAM_INPUT_SEED,
        );
    }
}

#[test]
fn jacobi_diagonal_is_exact() {
    let sp = TensorSpace::new(2, 3);
    let diag = sp.stiffness_diagonal();
    // Compare against applies to unit vectors on a sample of dofs.
    for &d in &[0usize, 17, 100, sp.ndof() - 1] {
        let mut e = vec![0.0f64; sp.ndof()];
        e[d] = 1.0;
        let col = sp.apply_stiffness(&e);
        assert!(
            (col[d] - diag[d]).abs() < 1e-12 * diag[d].abs().max(1.0),
            "diagonal mismatch at dof {d}: {} vs {}",
            col[d],
            diag[d]
        );
    }
    verdict(
        "jacobi-diag",
        "Kronecker diagonal == operator diagonal",
        FIXED_INPUT_SEED,
    );
}

#[test]
fn mms_orders_r1_through_r6() {
    // G1 slope gate: for each r, an h-ladder (m, 2m) must show
    // observed L2 order ≥ r + 0.6 (theory: r + 1). METRIC TRAP,
    // diagnosed during construction and kept as doctrine here: on a
    // SINGLE symmetric cell (m = 1) the sin-product solution aligns
    // with even-degree bubbles and superconverges, faking a shallow
    // slope at even r — so every ladder starts at m ≥ 2, where cells
    // see asymmetric arcs and the trap disappears (verified: both the
    // (1,1,1) and a mixed (2,1,3) eigenmode converge at order → r+1
    // on m ≥ 2 ladders; the probe survives as ho_probe.rs).
    let pi = std::f64::consts::PI;
    let u_exact = move |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let f_exact = move |p: [f64; 3]| 3.0 * pi * pi * u_exact(p);
    for r in 1..=6usize {
        let case = format!("mms-ho-order/r-{r}");
        // Low r needs finer ladders to reach asymptotics; high r must
        // stay above the 1e-13 solver floor.
        let ladder: [usize; 2] = if r <= 2 { [4, 8] } else { [2, 4] };
        let mut errs = Vec::new();
        for &m in &ladder {
            let sp = TensorSpace::new(m, r);
            let b = sp.load(&f_exact);
            let mask = sp.interior_mask();
            let diag = sp.stiffness_diagonal();
            let mut b_masked = b;
            for (bi, &mk) in b_masked.iter_mut().zip(&mask) {
                if !mk {
                    *bi = 0.0;
                }
            }
            let mut x = vec![0.0f64; sp.ndof()];
            let (iters, converged) = pcg_matfree(
                &|v| sp.apply_stiffness(v),
                &b_masked,
                &mut x,
                &mask,
                &diag,
                1e-13,
                20_000,
            );
            assert!(converged, "r={r} m={m}: PCG failed after {iters} iters");
            let err = sp.l2_error(&x, &u_exact);
            errs.push(err);
            measurement(
                &format!("{case}/measurement/m-{m}"),
                "mms-ho",
                format!(
                    "{{\"r\":{r},\"m\":{m},\"l2_error\":{},\"iterations\":{iters},\
                     \"input_seed\":{FIXED_INPUT_SEED}}}",
                    finite_json(err)
                ),
            );
        }
        let order = (errs[0] / errs[1]).ln() / (ladder[1] as f64 / ladder[0] as f64).ln();
        assert!(
            order > r as f64 + 0.6,
            "r={r}: observed order {order:.2} below gate {} (errors {errs:?})",
            r as f64 + 0.6
        );
        verdict(&case, &format!("r={r} order={order:.2}"), FIXED_INPUT_SEED);
    }
}

const GOLDEN_HASH: u64 = 0xaaf1_076a_196c_6902; // recorded at tfz.6 slice 1, frozen

#[test]
fn highorder_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // GL nodes/weights (n = 4, 7).
    for n in [4usize, 7] {
        let (x, w) = gauss_legendre(n);
        for v in x.iter().chain(&w) {
            feed(*v);
        }
    }
    // 1D element matrices at r = 4.
    let sp = TensorSpace::new(2, 4);
    for v in sp
        .element_mass_matrix()
        .iter()
        .chain(sp.element_stiffness_matrix())
    {
        feed(*v);
    }
    // Sum-factorized apply output sample.
    let u = rand_vec(sp.ndof(), 90);
    for v in sp.apply_stiffness(&u).iter().step_by(17) {
        feed(*v);
    }
    // Diagonal sample.
    for v in sp.stiffness_diagonal().iter().step_by(23) {
        feed(*v);
    }
    measurement(
        "ho-golden/measurement",
        "ho-golden",
        format!(
            "{{\"hash\":\"{acc:#018x}\",\"expected_hash\":\"{GOLDEN_HASH:#018x}\",\
             \"input_seed\":{STREAM_INPUT_SEED},\"stream_kernel\":\"0x460\",\
             \"stream_tile\":90}}"
        ),
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "high-order bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
