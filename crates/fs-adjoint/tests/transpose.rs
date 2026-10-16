//! Ledger-transposition conformance (the bk0o.1 bead; runs under the
//! `ledger-transpose` feature). Acceptance: a seam-crossing gradient
//! (control points → restriction → nonlinearity → SOLVE → functional)
//! matches the conditioning-aware FD falsifier; solver VJPs are
//! TRANSPOSED solves, never differentiation through iterations; a
//! missing VJP or a declared non-differentiable op BLOCKS the gradient
//! loudly; ⟨Av,w⟩=⟨v,Aᵀw⟩ holds for every registered linear op; revolve
//! checkpoints spilled through the REAL content-addressed ledger store
//! reproduce bit-equal gradients with and without spill.
#![cfg(feature = "ledger-transpose")]

use std::sync::Arc;

use fs_adjoint::transpose::{
    CheckpointStore, MemStore, Tape, TransposeError, TransposeProbeError, Vjp, VjpRegistry,
    check_transpose, fd_falsifier, spilled_adjoint,
};
use fs_adjoint::{HeatAdjoint, heat_initial_gradient};
use fs_sparse::Coo;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint/transpose\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

// ---- The seam ops -----------------------------------------------------

/// Linear restriction R (m×n smoothing kernel) — a Convert-like seam.
struct Restrict {
    m: usize,
    n: usize,
}

impl Restrict {
    fn weight(&self, r: usize, c: usize) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let x = (r * self.n + c) as f64;
        // A fixed smooth kernel: strictly positive, row-varying.
        0.25 / (1.0 + ((c as f64) - (r as f64) * (self.n as f64) / (self.m as f64)).powi(2))
            + 0.001 * (x * 0.37).sin()
    }

    fn apply(&self, x: &[f64]) -> Vec<f64> {
        (0..self.m)
            .map(|r| (0..self.n).map(|c| self.weight(r, c) * x[c]).sum())
            .collect()
    }

    fn apply_t(&self, y: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|c| (0..self.m).map(|r| self.weight(r, c) * y[r]).sum())
            .collect()
    }
}

impl Vjp for Restrict {
    fn vjp(&self, _primal: &[&[f64]], bar: &[f64]) -> Vec<Vec<f64>> {
        vec![self.apply_t(bar)]
    }
}

/// A smooth nonlinearity y_i = u_i + 0.1 u_i² (a blend-region stand-in).
struct SoftSquare;

impl SoftSquare {
    fn apply(u: &[f64]) -> Vec<f64> {
        u.iter().map(|v| v + 0.1 * v * v).collect()
    }
}

impl Vjp for SoftSquare {
    fn vjp(&self, primal: &[&[f64]], bar: &[f64]) -> Vec<Vec<f64>> {
        let u = primal[0];
        vec![
            bar.iter()
                .zip(u)
                .map(|(b, v)| b * 0.2f64.mul_add(*v, 1.0))
                .collect(),
        ]
    }
}

/// SPD tridiagonal solve y = A⁻¹x. The VJP is a TRANSPOSED SOLVE
/// (A symmetric ⇒ the same solve) on the cotangent — NEVER
/// differentiation through the iteration sequence.
struct SpdSolve {
    diag: f64,
    off: f64,
    n: usize,
}

impl SpdSolve {
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        // Deterministic CG, matrix-free tridiagonal apply.
        let apply = |x: &[f64]| -> Vec<f64> {
            (0..self.n)
                .map(|i| {
                    let mut v = self.diag * x[i];
                    if i > 0 {
                        v += self.off * x[i - 1];
                    }
                    if i + 1 < self.n {
                        v += self.off * x[i + 1];
                    }
                    v
                })
                .collect()
        };
        let mut x = vec![0.0f64; self.n];
        let mut r = b.to_vec();
        let mut p = r.clone();
        let mut rr: f64 = r.iter().map(|v| v * v).sum();
        for _ in 0..4 * self.n {
            if rr < 1e-28 {
                break;
            }
            let ap = apply(&p);
            let pap: f64 = p.iter().zip(&ap).map(|(a, c)| a * c).sum();
            let alpha = rr / pap;
            for i in 0..self.n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            let rr2: f64 = r.iter().map(|v| v * v).sum();
            let beta = rr2 / rr;
            rr = rr2;
            for i in 0..self.n {
                p[i] = r[i] + beta * p[i];
            }
        }
        x
    }
}

impl Vjp for SpdSolve {
    fn vjp(&self, _primal: &[&[f64]], bar: &[f64]) -> Vec<Vec<f64>> {
        vec![self.solve(bar)] // Aᵀ = A: the transposed solve.
    }
}

/// The lift-proxy functional J = g·y.
struct LiftProxy {
    g: Vec<f64>,
}

impl Vjp for LiftProxy {
    fn vjp(&self, _primal: &[&[f64]], bar: &[f64]) -> Vec<Vec<f64>> {
        vec![self.g.iter().map(|gi| gi * bar[0]).collect()]
    }
}

/// Forward pass of the whole seam chain, recording the tape.
fn run_chain(c: &[f64], registry_tape: Option<&mut Tape>) -> (f64, Option<usize>) {
    let restrict = Restrict { m: 12, n: c.len() };
    let solve = SpdSolve {
        diag: 4.0,
        off: -1.0,
        n: 12,
    };
    let g: Vec<f64> = (0..12).map(|i| 1.0 + 0.1 * f64::from(i as u8)).collect();
    let r = restrict.apply(c);
    let s = SoftSquare::apply(&r);
    let y = solve.solve(&s);
    let j: f64 = g.iter().zip(&y).map(|(a, b)| a * b).sum();
    if let Some(tape) = registry_tape {
        let leaf = tape.leaf(c.to_vec());
        let n1 = tape.apply("convert/restrict", &[leaf], r);
        let n2 = tape.apply("blend/soft-square", &[n1], s);
        let n3 = tape.apply("solver/spd", &[n2], y);
        let out = tape.apply("functional/lift-proxy", &[n3], vec![j]);
        (j, Some(out))
    } else {
        (j, None)
    }
}

fn registry() -> VjpRegistry {
    let mut reg = VjpRegistry::new();
    reg.register("convert/restrict", Arc::new(Restrict { m: 12, n: 8 }));
    reg.register("blend/soft-square", Arc::new(SoftSquare));
    reg.register(
        "solver/spd",
        Arc::new(SpdSolve {
            diag: 4.0,
            off: -1.0,
            n: 12,
        }),
    );
    reg.register(
        "functional/lift-proxy",
        Arc::new(LiftProxy {
            g: (0..12).map(|i| 1.0 + 0.1 * f64::from(i as u8)).collect(),
        }),
    );
    reg
}

#[test]
fn tr_001_seam_crossing_gradient_vs_fd_falsifier() {
    let c: Vec<f64> = (0..8).map(|i| 0.3 + 0.05 * f64::from(i as u8)).collect();
    let mut tape = Tape::new();
    let (_, out) = run_chain(&c, Some(&mut tape));
    let grads = tape
        .transpose(&registry(), out.expect("output"), &[1.0])
        .expect("full chain differentiates");
    let grad = grads.get(&0).expect("leaf gradient");
    assert_eq!(grad.len(), 8);
    // FD falsifier along 3 deterministic directions, conditioning-aware.
    let f = |x: &[f64]| run_chain(x, None).0;
    for k in 0..3 {
        let dir: Vec<f64> = (0..8)
            .map(|i| if i % 3 == k { 1.0 } else { 0.25 })
            .collect();
        let adjoint_dd: f64 = grad.iter().zip(&dir).map(|(g, d)| g * d).sum();
        let v = fd_falsifier(&f, &c, &dir, adjoint_dd, 1e-5, 1e-7);
        assert!(
            v.consistent,
            "seam-crossing gradient must satisfy FD: {v:?}"
        );
    }
    // Determinism: the transposed sweep is bit-equal on re-run.
    let mut tape2 = Tape::new();
    let (_, out2) = run_chain(&c, Some(&mut tape2));
    let grads2 = tape2
        .transpose(&registry(), out2.expect("output"), &[1.0])
        .expect("rerun");
    for (a, b) in grad.iter().zip(grads2.get(&0).expect("leaf")) {
        assert_eq!(a.to_bits(), b.to_bits(), "bit-equal gradients");
    }
    verdict(
        "tr-001",
        "control-points -> restrict -> blend -> SOLVE -> lift gradient passes the \
         conditioning-aware FD falsifier on 3 directions; re-runs bit-equal",
    );
}

#[test]
fn tr_002_transpose_consistency_battery() {
    let restrict = Restrict { m: 12, n: 8 };
    let report = check_transpose(&|x| restrict.apply(x), &|y| restrict.apply_t(y), 8, 12, 24)
        .expect("a well-posed probe request");
    assert_eq!(report.probes, 24);
    assert!(
        report.max_abs_residual < 1e-12,
        "restriction transpose exact: {report:?}"
    );
    assert!(
        report.scale > 0.0 && report.rel_residual < 1e-12,
        "the residual must also be small RELATIVE to the pairing scale: {report:?}"
    );
    let solve = SpdSolve {
        diag: 4.0,
        off: -1.0,
        n: 12,
    };
    let report_s = check_transpose(&|x| solve.solve(x), &|y| solve.solve(y), 12, 12, 24)
        .expect("a well-posed probe request");
    assert!(
        report_s.max_abs_residual < 1e-10,
        "symmetric solve is its own transpose to solver tolerance: {report_s:?}"
    );
    verdict(
        "tr-002",
        "the G0 battery: restriction transpose exact to 1e-12; solve self-transpose to \
         1e-10 over 24 probes each",
    );
}

#[test]
fn tr_003_missing_and_declared_vjps_fail_loud() {
    let c = vec![0.5; 8];
    // An op nobody registered lands mid-path.
    let mut tape = Tape::new();
    let leaf = tape.leaf(c.clone());
    // The lift-proxy VJP produces a 12-component cotangent (its `g`), so
    // the node it consumes must have 12 components: the sweep now CHECKS
    // cotangent lengths against the recorded primals instead of zipping
    // them into agreement.
    let v = vec![1.0; 12];
    let n1 = tape.apply("mystery/op", &[leaf], v.clone());
    let out = tape.apply("functional/lift-proxy", &[n1], vec![1.0]);
    let reg = registry();
    let err = tape
        .transpose(&reg, out, &[1.0])
        .expect_err("missing VJP blocks");
    assert!(
        matches!(&err, TransposeError::MissingVjp { op } if op == "mystery/op"),
        "the error names the op: {err}"
    );
    assert!(
        format!("{err}").contains("Goodhart"),
        "teaches the trap: {err}"
    );
    // A DECLARED non-differentiable op blocks with color consequences.
    let mut reg2 = registry();
    reg2.declare_non_differentiable(
        "mystery/op",
        "integer quantization has no useful derivative",
        "estimated at best",
    );
    let err2 = tape
        .transpose(&reg2, out, &[1.0])
        .expect_err("declared op still blocks");
    match &err2 {
        TransposeError::NonDifferentiableInPath {
            op,
            color_consequence,
            ..
        } => {
            assert_eq!(op, "mystery/op");
            assert_eq!(color_consequence, "estimated at best");
        }
        TransposeError::MissingVjp { .. } => panic!("wrong error kind: {err2}"),
    }
    // Coverage report shows both kinds.
    let (diff, nondiff) = reg2.coverage();
    assert_eq!(diff.len(), 4);
    assert_eq!(nondiff, vec!["mystery/op"]);
    verdict(
        "tr-003",
        "missing VJP blocks loudly naming the op and the Goodhart trap; declared \
         non-differentiable blocks with its color consequence; coverage reported",
    );
}

/// The REAL content-addressed store: fs-ledger artifacts (dev-dep —
/// the shared storage discipline with Proposal 2).
struct LedgerCas {
    ledger: fs_ledger::Ledger,
}

impl CheckpointStore for LedgerCas {
    fn put(&mut self, bytes: &[u8]) -> Vec<u8> {
        let receipt = self
            .ledger
            .put_artifact("revolve-checkpoint", bytes, None)
            .expect("cas put");
        receipt.hash.as_bytes().to_vec()
    }

    fn get(&self, key: &[u8]) -> Vec<u8> {
        let mut h = [0u8; 32];
        h.copy_from_slice(key);
        self.ledger
            .get_artifact(&fs_ledger::ContentHash(h))
            .expect("cas get")
            .expect("checkpoint present")
    }
}

/// A single backward-Euler step as its own problem (M = I, 1-D
/// Laplacian K) — the deterministic step function both sweeps share.
fn one_step_problem(n: usize) -> HeatAdjoint {
    let mut m = Coo::new(n, n);
    let mut k = Coo::new(n, n);
    for i in 0..n {
        m.push(i, i, 1.0);
        k.push(i, i, 2.0);
        if i > 0 {
            k.push(i, i - 1, -1.0);
        }
        if i + 1 < n {
            k.push(i, i + 1, -1.0);
        }
    }
    HeatAdjoint::new(m.assemble(), &k.assemble(), 0.05, 1)
}

#[test]
fn tr_004_cas_checkpoints_bit_equal_with_and_without_spill() {
    // A small heat problem (mass = I, 1-D Laplacian stiffness).
    let n = 24;
    let mut mass = Coo::new(n, n);
    let mut stiff = Coo::new(n, n);
    for i in 0..n {
        mass.push(i, i, 1.0);
        stiff.push(i, i, 2.0);
        if i > 0 {
            stiff.push(i, i - 1, -1.0);
        }
        if i + 1 < n {
            stiff.push(i, i + 1, -1.0);
        }
    }
    let problem = HeatAdjoint::new(mass.assemble(), &stiff.assemble(), 0.05, 12);
    let u0: Vec<f64> = (0..n).map(|i| (f64::from(i as u8) * 0.4).sin()).collect();
    let target = vec![0.0f64; n];
    let one_f = one_step_problem(n);
    let step_f = move |u: &[f64]| one_f.forward(u);
    let one_r = one_step_problem(n);
    // μ = M·A⁻ᵀ·λ with M = I: the single reverse step is the forward
    // solve on the cotangent (symmetric system).
    let step_r = move |bar: &[f64]| one_r.forward(bar);
    let terminal =
        |u_n: &[f64]| -> Vec<f64> { u_n.iter().zip(&target).map(|(a, b)| a - b).collect() };
    // No-spill baseline (MemStore) vs the REAL ledger CAS.
    let mut mem = MemStore::default();
    let (g_mem, cp_mem, _) = spilled_adjoint(&u0, 12, 4, &mut mem, &step_f, &step_r, &terminal);
    let dir = std::env::temp_dir().join(format!("fs-adjoint-cas-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let mut cas = LedgerCas {
        ledger: fs_ledger::Ledger::open(dir.join("cas.led").to_str().expect("utf8"))
            .expect("ledger"),
    };
    let (g_cas, cp_cas, fwd) = spilled_adjoint(&u0, 12, 4, &mut cas, &step_f, &step_r, &terminal);
    assert_eq!(cp_mem, 3);
    assert_eq!(cp_cas, 3, "3 checkpoints at stride 4 over 12 steps");
    for (a, b) in g_mem.iter().zip(&g_cas) {
        assert_eq!(a.to_bits(), b.to_bits(), "BIT-EQUAL with and without spill");
    }
    // And the whole thing agrees with the base crate's revolve gradient.
    let (g_revolve, _) = heat_initial_gradient(&problem, &u0, &target);
    for (a, b) in g_mem.iter().zip(&g_revolve) {
        assert!(
            (a - b).abs() < 1e-9,
            "uniform-checkpoint vs revolve gradients agree: {a} vs {b}"
        );
    }
    println!(
        "{{\"metric\":\"cas-checkpoints\",\"checkpoints\":{cp_cas},\"forward_evals\":{fwd},\
         \"steps\":12}}"
    );
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "tr-004",
        "checkpoints spilled through the real fs-ledger CAS reproduce BIT-EQUAL \
         gradients vs in-memory; agrees with the base revolve path to 1e-9",
    );
}

// ---- Claim-integrity regressions (E02 sweep) --------------------------

/// A VJP that returns a cotangent SHORTER than its input — the masking
/// or slicing bug the shape contract exists to catch.
struct TruncatingVjp {
    keep: usize,
}

impl Vjp for TruncatingVjp {
    fn vjp(&self, _primal: &[&[f64]], bar: &[f64]) -> Vec<Vec<f64>> {
        vec![bar.iter().take(self.keep).copied().collect()]
    }
}

#[test]
#[should_panic(expected = "VJP cotangent length for input 0")]
fn tr_005_truncated_cotangent_is_refused_not_zipped() {
    // frankensim-extreal-program-f85xj.2.23. Only ARITY was asserted;
    // individual cotangent lengths were never compared against the
    // recorded primal, and accumulation used `zip`, which truncates to
    // the shorter. A length-k cotangent was then presented as the
    // length-n leaf gradient with the dropped tail indistinguishable
    // from a genuine zero.
    let mut reg = VjpRegistry::new();
    reg.register("masking/op", Arc::new(TruncatingVjp { keep: 3 }));
    let mut tape = Tape::new();
    let leaf = tape.leaf(vec![1.0; 8]);
    let out = tape.apply("masking/op", &[leaf], vec![1.0; 8]);
    let _ = tape.transpose(&reg, out, &[1.0; 8]);
}

#[test]
#[should_panic(expected = "seed cotangent length")]
fn tr_006_short_seed_is_refused() {
    // Same shape, on the caller's side of the seam.
    let mut reg = VjpRegistry::new();
    reg.register("masking/op", Arc::new(TruncatingVjp { keep: 8 }));
    let mut tape = Tape::new();
    let leaf = tape.leaf(vec![1.0; 8]);
    let out = tape.apply("masking/op", &[leaf], vec![1.0; 8]);
    let _ = tape.transpose(&reg, out, &[1.0; 3]);
}

#[test]
#[should_panic(expected = "reserved leaf sentinel")]
fn tr_007_reserved_leaf_op_name_is_refused() {
    // The sweep skips any node whose op is "leaf" WITHOUT a registry
    // lookup, so recording a real op under that name terminated the
    // chain silently instead of raising MissingVjp.
    let mut tape = Tape::new();
    let leaf = tape.leaf(vec![1.0; 4]);
    let _ = tape.apply("leaf", &[leaf], vec![1.0; 4]);
}

#[test]
fn tr_008_transpose_probe_refuses_vacuous_requests() {
    // frankensim-extreal-program-f85xj.2.22. Zero probes returned 0.0 —
    // the STRONGEST possible score — having evaluated neither operator,
    // and an empty dimension made both inner products empty sums.
    let restrict = Restrict { m: 12, n: 8 };
    let apply = |x: &[f64]| restrict.apply(x);
    let apply_t = |y: &[f64]| restrict.apply_t(y);
    assert_eq!(
        check_transpose(&apply, &apply_t, 8, 12, 0),
        Err(TransposeProbeError::NoProbes),
        "zero probes must refuse, not score 0.0"
    );
    assert_eq!(
        check_transpose(&|_x| Vec::new(), &|_y| Vec::new(), 0, 0, 4),
        Err(TransposeProbeError::EmptyDimension { n_in: 0, n_out: 0 }),
        "empty dimensions agree vacuously"
    );
    let short = check_transpose(&|_x| vec![0.0; 3], &apply_t, 8, 12, 4)
        .expect_err("a wrong-length image must refuse");
    assert!(
        matches!(
            short,
            TransposeProbeError::ImageLength {
                probe: 0,
                operator: "apply",
                expected: 12,
                got: 3
            }
        ),
        "{short:?}"
    );
    let nan_apply = |x: &[f64]| x.iter().map(|_| f64::NAN).collect::<Vec<f64>>();
    let identity_t = |y: &[f64]| y.to_vec();
    assert!(matches!(
        check_transpose(&nan_apply, &identity_t, 8, 8, 4),
        Err(TransposeProbeError::NonFinite {
            probe: 0,
            operator: "apply"
        })
    ));

    // And the ABSOLUTE residual is no longer the whole story: a tiny
    // operator with a badly wrong transpose scores ~1e-40 in absolute
    // terms while its RELATIVE residual is O(1).
    let tiny = |x: &[f64]| x.iter().map(|v| 1e-20 * v).collect::<Vec<f64>>();
    let tiny_wrong_t = |y: &[f64]| y.iter().map(|v| -1e-20 * v).collect::<Vec<f64>>();
    let report = check_transpose(&tiny, &tiny_wrong_t, 6, 6, 8).expect("well-posed");
    assert!(
        report.max_abs_residual < 1e-19,
        "the absolute residual is vacuously small: {report:?}"
    );
    assert!(
        report.rel_residual > 1.0,
        "the relative residual exposes the sign-flipped transpose: {report:?}"
    );
    verdict(
        "tr-008",
        "the G0 probe refuses zero probes, empty dimensions, wrong-length and non-finite \
         images, and reports a scale so thresholds can be relative",
    );
}

#[test]
fn tr_009_fd_falsifier_band_is_scale_relative() {
    // frankensim-extreal-program-f85xj.2.21. With the production
    // constants of certs::fd_spot_checks (h = 1e-5, base_tol = 1e-7),
    // the old band was `base_tol * max(|dd|, |fd|, 1.0)` — an ABSOLUTE
    // floor. On a derivative smaller than base_tol every value in
    // [-base_tol, base_tol] passed, including 0 and the sign flip.
    let tiny_slope = |x: &[f64]| 1e-9 * x[0];
    let x = [1.0];
    let dir = [0.4];

    let silent_zero = fd_falsifier(&tiny_slope, &x, &dir, 0.0, 1e-5, 1e-7);
    assert!(
        silent_zero.falsifiable,
        "4e-10 is far above the difference-quotient noise floor: {silent_zero:?}"
    );
    assert!(
        !silent_zero.consistent,
        "a zero adjoint must not agree with a nonzero FD derivative: {silent_zero:?}"
    );

    let sign_flipped = fd_falsifier(&tiny_slope, &x, &dir, -4e-10, 1e-5, 1e-7);
    assert!(
        !sign_flipped.consistent,
        "an exact sign flip must be caught: {sign_flipped:?}"
    );

    let honest = fd_falsifier(&tiny_slope, &x, &dir, 4e-10, 1e-5, 1e-7);
    assert!(
        honest.consistent,
        "the correct adjoint must still pass: {honest:?}"
    );

    // Control from the bead: at unit scale the verdict is unchanged.
    let unit_slope = |x: &[f64]| 3.0 * x[0];
    assert!(!fd_falsifier(&unit_slope, &x, &dir, 0.0, 1e-5, 1e-7).consistent);
    assert!(fd_falsifier(&unit_slope, &x, &dir, 1.2, 1e-5, 1e-7).consistent);

    // A bit-insensitive objective gives NO SIGNAL — which is not
    // agreement, even though |0 - 0| <= tolerance holds trivially.
    let constant = |_x: &[f64]| 7.0;
    let no_signal = fd_falsifier(&constant, &x, &dir, 0.0, 1e-5, 1e-7);
    assert!(
        !no_signal.falsifiable,
        "a constant objective cannot falsify anything: {no_signal:?}"
    );
    assert!(
        !no_signal.consistent,
        "no signal must not be reported as consistent: {no_signal:?}"
    );
    verdict(
        "tr-009",
        "the falsifier band is homogeneous in the derivative scale: zero and sign-flipped \
         adjoints on a 1e-9 slope are caught, and a bit-insensitive probe reports no signal",
    );
}
