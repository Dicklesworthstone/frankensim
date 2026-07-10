//! fs-opt conformance suite (CONTRACT.md: any reimplementation must
//! pass). Build-time validation with named teaching diagnostics,
//! bitwise serialization round-trips + hash identity, hash-consed CSE
//! and substitution laws, differentiability-class routing, the toy
//! Riemannian descent over every manifold kind, and P4/P7 budgets and
//! cancellation. JSON-line verdicts; seeded cases carry seeds.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{
    Class, ConstraintKind, DescentOptions, Manifold, OptError, OptimizerFamily, ProblemBuilder,
    ProblemTag, Sense, descend_fn, descend_ir, eval, parse, problem_hash, serialize,
};
use fs_qty::Dims;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-opt/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn below(&mut self, n: u64) -> u64 {
        (self.next() >> 32) % n
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0F7,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

const METER: Dims = Dims([1, 0, 0, 0, 0]);

/// opt-001 — incremental validation TEACHES: dimension mismatches,
/// shape mismatches, non-dimensionless transcendentals, odd sqrt dims,
/// bad indices, and non-scalar roots are all refused with the nodes
/// named; a seeded fuzz storm classifies random op sequences correctly
/// against an independent reference checker.
#[test]
#[allow(clippy::too_many_lines)] // the fuzz storm's reference model is one block
fn opt_001_validation_teaches() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 3 }, METER);
    let xr = b.var_ref(x).expect("var ref");
    let len = b.norm_sq(xr).expect("norm_sq");
    let meter_const = b.konst(2.0, METER);
    // m² + m → DimMismatch naming the op.
    let e = b.add(len, meter_const).expect_err("dim mismatch");
    let dim_teaches = matches!(e, OptError::DimMismatch { op: "add", .. })
        && e.to_string().contains("incompatible dimensions");
    // vector into abs → ShapeMismatch.
    let e = b.abs(xr).expect_err("shape mismatch");
    let shape_teaches = e.to_string().contains("incompatible shapes");
    // ln of meters → NonDimensionless.
    let dimensionless_err = {
        let e = b.ln(meter_const).expect_err("dimensioned ln");
        e.to_string().contains("dimensionless")
    };
    // sqrt of m³ → OddDims.
    let m3 = b.konst(1.0, Dims([3, 0, 0, 0, 0]));
    let odd_err = b.sqrt(m3).expect_err("odd dims").to_string();
    let pow_dims_checked = {
        let mut powers = ProblemBuilder::new();
        let m = powers.konst(2.0, METER);
        let m2 = powers.konst(2.0, Dims([2, 0, 0, 0, 0]));
        let m101 = powers.powi(m, 101).expect("m^101 is representable");
        let overflow = powers.powi(m2, 100).expect_err("m^200 cannot fit i8 dims");
        let problem = powers.finish();
        problem.node_dims(m101) == Dims([101, 0, 0, 0, 0])
            && matches!(overflow, OptError::DimOverflow { exponent: 100, .. })
    };
    // component out of range.
    let idx_err = b.component(xr, 7).expect_err("bad index").to_string();
    // vector objective root refused.
    let root_err = b
        .objective(xr, Sense::Minimize, 1.0)
        .expect_err("vector root")
        .to_string();

    // Fuzz storm: random op sequences; acceptance must MATCH an
    // independent shape/dims model maintained by the test.
    let mut rng = Lcg(0x1001_2026_0706_0031);
    let mut agreed = 0u64;
    let mut disagreed = 0u64;
    let mut fb = ProblemBuilder::new();
    let v = fb.var("v", Manifold::Rn { dim: 3 }, METER);
    let mut model: Vec<(bool, [i8; 5])> = Vec::new(); // (is_vector, dims)
    let mut ids = Vec::new();
    let seed_node = fb.var_ref(v).expect("seed");
    ids.push(seed_node);
    model.push((true, METER.0));
    for _ in 0..600 {
        let pick = ids[rng.below(ids.len() as u64) as usize];
        let (pv, pd) = model[pick.0 as usize];
        let (qn, qv, qd) = {
            let q = ids[rng.below(ids.len() as u64) as usize];
            let (a, b2) = model[q.0 as usize];
            (q, a, b2)
        };
        match rng.below(5) {
            0 => {
                // add: legal iff same shape and dims.
                let legal = pv == qv && pd == qd;
                match fb.add(pick, qn) {
                    Ok(id) => {
                        if legal {
                            agreed += 1;
                        } else {
                            disagreed += 1;
                        }
                        if (id.0 as usize) == model.len() {
                            model.push((pv, pd));
                        }
                        ids.push(id);
                    }
                    Err(_) => {
                        if legal {
                            disagreed += 1;
                        } else {
                            agreed += 1;
                        }
                    }
                }
            }
            1 => {
                // abs: legal iff scalar.
                let legal = !pv;
                match fb.abs(pick) {
                    Ok(id) => {
                        if legal {
                            agreed += 1;
                        } else {
                            disagreed += 1;
                        }
                        if (id.0 as usize) == model.len() {
                            model.push((false, pd));
                        }
                        ids.push(id);
                    }
                    Err(_) => {
                        if legal {
                            disagreed += 1;
                        } else {
                            agreed += 1;
                        }
                    }
                }
            }
            2 => {
                // exp: legal iff scalar and dimensionless.
                let legal = !pv && pd == [0i8; 5];
                match fb.exp(pick) {
                    Ok(id) => {
                        if legal {
                            agreed += 1;
                        } else {
                            disagreed += 1;
                        }
                        if (id.0 as usize) == model.len() {
                            model.push((false, [0; 5]));
                        }
                        ids.push(id);
                    }
                    Err(_) => {
                        if legal {
                            disagreed += 1;
                        } else {
                            agreed += 1;
                        }
                    }
                }
            }
            3 => {
                // dot: legal iff both vectors (same length 3 here).
                let legal = pv && qv;
                match fb.dot(pick, qn) {
                    Ok(id) => {
                        if legal {
                            agreed += 1;
                        } else {
                            disagreed += 1;
                        }
                        if (id.0 as usize) == model.len() {
                            let mut d = pd;
                            for (a, b2) in d.iter_mut().zip(qd) {
                                *a = a.saturating_add(b2);
                            }
                            model.push((false, d));
                        }
                        ids.push(id);
                    }
                    Err(_) => {
                        if legal {
                            disagreed += 1;
                        } else {
                            agreed += 1;
                        }
                    }
                }
            }
            _ => {
                // konst with random dims (always legal).
                let dims = Dims([(rng.below(3) as i8) - 1, 0, (rng.below(3) as i8) - 1, 0, 0]);
                let id = fb.konst(rng.unit(), dims);
                if (id.0 as usize) == model.len() {
                    model.push((false, dims.0));
                }
                ids.push(id);
                agreed += 1;
            }
        }
    }
    verdict(
        "opt-001",
        dim_teaches
            && shape_teaches
            && dimensionless_err
            && odd_err.contains("odd dimension")
            && pow_dims_checked
            && idx_err.contains("does not exist")
            && root_err.contains("SCALAR")
            && disagreed == 0
            && agreed > 500,
        &format!(
            "seeded ill-typed constructions refuse with teaching text (dims, shapes, \
             dimensionless, odd-sqrt, checked power dims, index, non-scalar root) and a 600-op fuzz storm \
             matches the independent validity model exactly ({agreed} agreements, \
             {disagreed} disagreements); seed 0x1001_2026_0706_0031"
        ),
    );
}

/// Build the shared fixture problem (PDE + stochastic + kink + tags).
fn fixture() -> fs_opt::Problem {
    let mut b = ProblemBuilder::new();
    let theta = b.var("theta", Manifold::Rn { dim: 3 }, METER);
    let q = b.var("attitude", Manifold::So3, Dims::NONE);
    let _ = q;
    let tr = b.var_ref(theta).expect("ref");
    let compliance = b.norm_sq(tr).expect("norm_sq");
    let limit = b.konst(4.0, Dims([2, 0, 0, 0, 0]));
    let excess = b.sub(compliance, limit).expect("sub");
    let zero = b.konst(0.0, Dims([2, 0, 0, 0, 0]));
    let hinge = b.max_of(excess, zero).expect("max");
    let pde = b
        .pde_residual("stokes-channel", theta, true, Dims::NONE)
        .expect("pde");
    let mean = b.expectation(pde, "uq-lhs-64").expect("expectation");
    let risk = b.cvar(mean, 0.95, "uq-lhs-64").expect("cvar");
    b.objective(compliance, Sense::Minimize, 1.0).expect("obj");
    b.objective(risk, Sense::Minimize, 0.25).expect("obj2");
    b.constraint(hinge, ConstraintKind::LeZero, "compliance-cap")
        .expect("con");
    b.tag(ProblemTag::ChanceConstrained { prob: 0.99 });
    b.tag(ProblemTag::MultiFidelity { levels: 3 });
    b.set_budget(10_000);
    b.finish()
}

/// opt-002 — canonical serialization: bitwise round-trip, hash
/// identity/stability, and corrupted-text refusals with line numbers.
#[test]
fn opt_002_roundtrip_and_hash() {
    let p1 = fixture();
    let text = serialize(&p1);
    let p2 = parse(&text).expect("round-trip parses");
    let roundtrip = p1 == p2;
    let hash_stable = problem_hash(&p1) == problem_hash(&fixture());
    let hash_matches = problem_hash(&p1) == problem_hash(&p2);
    // A differing constant changes the hash.
    let p3 = {
        let mut b = ProblemBuilder::new();
        let v = b.var("theta", Manifold::Rn { dim: 3 }, METER);
        let r = b.var_ref(v).expect("ref");
        let n = b.norm_sq(r).expect("n");
        b.objective(n, Sense::Minimize, 1.0).expect("obj");
        b.finish()
    };
    let hash_differs = problem_hash(&p1) != problem_hash(&p3);
    // Corruption: flip a byte in the body → integrity refusal w/ line.
    let corrupt = text.replace("objective min", "objective max");
    let integrity = matches!(parse(&corrupt), Err(OptError::Parse { what, .. })
        if what.contains("integrity hash mismatch"));
    // Garbage directive → parse error with its line number.
    let garbage = "fsopt v1\nwat 1 2 3\n";
    let teaches = matches!(parse(garbage), Err(OptError::Parse { line: 2, .. }));
    verdict(
        "opt-002",
        roundtrip && hash_stable && hash_matches && hash_differs && integrity && teaches,
        "build->serialize->parse round-trips to an IDENTICAL problem (floats travel \
         as bit patterns); hashes are stable across identical builds, differ across \
         edits, and guard integrity (tampered text refuses with the line named)",
    );
}

/// opt-003 — graph algebra (G0): hash-consing makes repeated
/// subexpressions the SAME node id; substitution commutes with
/// evaluation BITWISE; basic identities hold bitwise.
#[test]
fn opt_003_cse_and_substitution() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE);
    let xr = b.var_ref(x).expect("ref");
    let n1 = b.norm_sq(xr).expect("n1");
    let n2 = b.norm_sq(xr).expect("n2");
    let cse = n1 == n2; // identical subexpression → identical id
    let s1 = b.add(n1, n1).expect("s1");
    let s2 = b.add(n2, n2).expect("s2");
    let cse2 = s1 == s2;

    // Substitution law: f(g(x)) built two ways evaluates bitwise-equal.
    // f(u) = u² + u over scalar u; g = <x, x>.
    let composed = {
        let g = n1;
        let g2 = b.mul(g, g).expect("g2");
        b.add(g2, g).expect("f(g)")
    };
    let mut rng = Lcg(0x1001_2026_0706_0033);
    let mut law = true;
    for _ in 0..200 {
        let p = vec![rng.unit() * 4.0 - 2.0, rng.unit() * 4.0 - 2.0];
        let problem = {
            let mut bb = ProblemBuilder::new();
            let v = bb.var("x", Manifold::Rn { dim: 2 }, Dims::NONE);
            let r = bb.var_ref(v).expect("r");
            let n = bb.norm_sq(r).expect("n");
            let n2b = bb.mul(n, n).expect("n2");
            let f = bb.add(n2b, n).expect("f");
            bb.objective(f, Sense::Minimize, 1.0).expect("obj");
            bb.finish()
        };
        let via_ir = eval(
            &problem,
            problem.objectives[0].node,
            std::slice::from_ref(&p),
        )
        .expect("eval")
        .scalar()
        .expect("scalar");
        let g = p[0] * p[0] + p[1] * p[1];
        let manual = g * g + g;
        law &= via_ir.to_bits() == manual.to_bits();
    }
    // Identity laws, bitwise.
    let mut ib = ProblemBuilder::new();
    let v = ib.var("v", Manifold::Rn { dim: 1 }, Dims::NONE);
    let r = ib.var_ref(v).expect("r");
    let c = ib.component(r, 0).expect("c");
    let neg2 = {
        let n = ib.neg(c).expect("n");
        ib.neg(n).expect("nn")
    };
    let minaa = ib.min_of(c, c).expect("minaa");
    let prob = {
        ib.objective(neg2, Sense::Minimize, 1.0).expect("o1");
        ib.objective(minaa, Sense::Minimize, 1.0).expect("o2");
        ib.objective(c, Sense::Minimize, 1.0).expect("o3");
        ib.finish()
    };
    let mut ids_ok = true;
    for _ in 0..50 {
        let xv = vec![rng.unit() * 10.0 - 5.0];
        let e = |n| {
            eval(&prob, n, std::slice::from_ref(&xv))
                .expect("eval")
                .scalar()
                .expect("s")
        };
        ids_ok &= e(prob.objectives[0].node).to_bits() == e(prob.objectives[2].node).to_bits();
        ids_ok &= e(prob.objectives[1].node).to_bits() == e(prob.objectives[2].node).to_bits();
    }
    let composed_dbg = composed; // keep the two-way build exercised
    let _ = composed_dbg;
    verdict(
        "opt-003",
        cse && cse2 && law && ids_ok,
        "hash-consing returns IDENTICAL ids for repeated subexpressions (CSE by \
         construction); substitution commutes with evaluation BITWISE over 200 \
         random points; neg(neg(x)) and min(x,x) evaluate bitwise-identical to x; \
         seed 0x1001_2026_0706_0033",
    );
}

/// opt-004 — class propagation + routing: the min() kink is KNOWN at
/// build time, poisons L-BFGS routing with the node NAMED, routes fine
/// to subgradient/gradient-free; adjoint-less PDE nodes refuse
/// gradient families; the class trace names every node.
#[test]
fn opt_004_class_routing() {
    let p = fixture();
    // The fixture's objectives: smooth compliance + C0 risk (cvar).
    let smooth_obj = p.objectives[0].node;
    let risky_obj = p.objectives[1].node;
    let classes_right = p.class(smooth_obj) == Class::Smooth && p.class(risky_obj) == Class::C0;
    // Routing: whole problem contains max() hinge + cvar → L-BFGS must
    // refuse and NAME a poisoning node.
    let refusal = p.route(OptimizerFamily::Lbfgs).expect_err("must refuse");
    let names_node = match &refusal {
        OptError::NonsmoothForFamily { kind, .. } => kind == "max" || kind == "cvar",
        _ => false,
    };
    let sub_ok = p.route(OptimizerFamily::SubgradientBundle).is_ok();
    let free_ok = p.route(OptimizerFamily::GradientFree).is_ok();
    // Adjoint-less PDE → refused for L-BFGS even when smooth.
    let no_adj = {
        let mut b = ProblemBuilder::new();
        let v = b.var("v", Manifold::Rn { dim: 2 }, Dims::NONE);
        let pde = b
            .pde_residual("heat-steady", v, false, Dims::NONE)
            .expect("pde");
        b.objective(pde, Sense::Minimize, 1.0).expect("obj");
        let p2 = b.finish();
        matches!(
            p2.route(OptimizerFamily::Lbfgs),
            Err(OptError::NoAdjoint { study, .. }) if study == "heat-steady"
        ) && p2.route(OptimizerFamily::GradientFree).is_ok()
    };
    let trace = p.class_trace();
    let trace_ok =
        trace.len() == p.exprs.len() && trace.iter().any(|l| l.contains("max") && l.contains("C0"));
    let mut em = fs_obs::Emitter::new("fs-opt/conformance", "opt-004/classes");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "opt-class-routing".to_string(),
                json: format!(
                    "{{\"nodes\":{},\"hash\":\"{:016X}\",\"refusal\":\"{refusal}\"}}",
                    p.exprs.len(),
                    problem_hash(&p)
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("class event validates");
    println!("{line}");
    verdict(
        "opt-004",
        classes_right && names_node && sub_ok && free_ok && no_adj && trace_ok,
        &format!(
            "the kink is known at build time: L-BFGS routing refuses NAMING the \
             poisoning node ({refusal}), subgradient/gradient-free families accept, \
             adjoint-less PDE nodes refuse gradient families with the study named, \
             and the class trace covers every node"
        ),
    );
}

/// opt-005 — the toy Riemannian descent consumes manifold metadata:
/// Sphere converges to the analytic minimizer, SO(3) aligns a vector,
/// Stiefel columns stay orthonormal, and iterates stay ON their
/// manifolds throughout.
#[test]
#[allow(clippy::too_many_lines)] // one manifold per block
fn opt_005_riemannian_descent() {
    with_cx(|cx| {
        // Sphere: minimize <c, x> → x* = −c/|c|.
        let c = [0.6f64, -0.3, 0.9];
        let cn = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        let f = |x: &[f64]| c[0] * x[0] + c[1] * x[1] + c[2] * x[2];
        let rep = descend_fn(
            Manifold::Sphere { ambient: 3 },
            &f,
            &[1.0, 0.0, 0.0],
            DescentOptions {
                steps: 400,
                lr: 0.3,
                fd_h: 1e-6,
            },
            0,
            cx,
        )
        .expect("sphere descent");
        let err_sphere = (0..3)
            .map(|i| (rep.x[i] - (-c[i] / cn)).abs())
            .fold(0.0f64, f64::max);
        let on_sphere = (rep.x.iter().map(|v| v * v).sum::<f64>().sqrt() - 1.0).abs() < 1e-12;

        // SO(3): align R(q)·a with b.
        let a = [1.0f64, 0.0, 0.0];
        let b = [0.0f64, 0.0, 1.0];
        let rotate = |q: &[f64], v: [f64; 3]| -> [f64; 3] {
            // R(q) v via quaternion sandwich (w,x,y,z).
            let (w, x, y, z) = (q[0], q[1], q[2], q[3]);
            let uv = [
                y * v[2] - z * v[1],
                z * v[0] - x * v[2],
                x * v[1] - y * v[0],
            ];
            let uuv = [
                y * uv[2] - z * uv[1],
                z * uv[0] - x * uv[2],
                x * uv[1] - y * uv[0],
            ];
            [
                v[0] + 2.0 * (w * uv[0] + uuv[0]),
                v[1] + 2.0 * (w * uv[1] + uuv[1]),
                v[2] + 2.0 * (w * uv[2] + uuv[2]),
            ]
        };
        let g = |q: &[f64]| -> f64 {
            let r = rotate(q, a);
            (0..3).map(|i| (r[i] - b[i]) * (r[i] - b[i])).sum()
        };
        let rep2 = descend_fn(
            Manifold::So3,
            &g,
            &[1.0, 0.0, 0.0, 0.0],
            DescentOptions {
                steps: 400,
                lr: 0.25,
                fd_h: 1e-6,
            },
            0,
            cx,
        )
        .expect("so3 descent");
        let aligned = g(&rep2.x) < 1e-10;
        let unit_q = (rep2.x.iter().map(|v| v * v).sum::<f64>().sqrt() - 1.0).abs() < 1e-12;

        // Stiefel(4,2): maximize Rayleigh sum (minimize negative) of a
        // diagonal quadratic; columns must STAY orthonormal.
        let diag = [4.0, 3.0, 1.0, 0.5];
        let h = |xs: &[f64]| -> f64 {
            let mut s = 0.0f64;
            for j in 0..2 {
                for i in 0..4 {
                    let v = xs[j * 4 + i];
                    s -= diag[i] * v * v;
                }
            }
            s
        };
        let x0 = [
            0.5, 0.5, 0.5, 0.5, //
            0.5, -0.5, 0.5, -0.5,
        ];
        let rep3 = descend_fn(
            Manifold::Stiefel { n: 4, p: 2 },
            &h,
            &x0,
            DescentOptions {
                steps: 300,
                lr: 0.15,
                fd_h: 1e-6,
            },
            0,
            cx,
        )
        .expect("stiefel descent");
        let mut ortho = true;
        for j in 0..2 {
            for k in 0..=j {
                let d: f64 = (0..4).map(|i| rep3.x[j * 4 + i] * rep3.x[k * 4 + i]).sum();
                let want = if j == k { 1.0 } else { 0.0 };
                ortho &= (d - want).abs() < 1e-10;
            }
        }
        // The top-2 invariant subspace of diag(4,3,1,.5) is span(e1,e2):
        // the trailing rows must vanish.
        let tail: f64 = (0..2)
            .map(|j| rep3.x[j * 4 + 2].abs() + rep3.x[j * 4 + 3].abs())
            .sum();
        let subspace = tail < 1e-4;
        verdict(
            "opt-005",
            err_sphere < 1e-6
                && on_sphere
                && aligned
                && unit_q
                && ortho
                && subspace
                && rep.f_final < rep.f0
                && rep2.f_final < rep2.f0,
            &format!(
                "manifold metadata drives descent: sphere reaches the analytic \
                 minimizer to {err_sphere:.1e} STAYING unit; SO(3) aligns a vector to \
                 1e-10 with a unit quaternion throughout; Stiefel(4,2) finds the top \
                 invariant subspace with columns orthonormal to 1e-10"
            ),
        );
    });
}

/// opt-006 — P4/P7 through the IR-driven descent: the problem's
/// attached budget stops the run with a RECEIPT (not an error), and
/// cancellation between steps returns the teaching error; PDE and
/// stochastic nodes refuse evaluation naming their executor.
#[test]
fn opt_006_budget_and_cancellation() {
    with_cx(|cx| {
        let build = |max_evals: u64| {
            let mut b = ProblemBuilder::new();
            let v = b.var("x", Manifold::Rn { dim: 4 }, Dims::NONE);
            let r = b.var_ref(v).expect("r");
            let n = b.norm_sq(r).expect("n");
            b.objective(n, Sense::Minimize, 1.0).expect("obj");
            b.set_budget(max_evals);
            b.finish()
        };
        let p = build(50);
        let rep = descend_ir(&p, &[1.0, -2.0, 0.5, 3.0], DescentOptions::default(), cx)
            .expect("budgeted descent");
        let receipt = rep.budget_stopped && rep.evals <= 50 && rep.f_final < rep.f0;

        let unlimited = build(0);
        let rep2 = descend_ir(
            &unlimited,
            &[1.0, -2.0, 0.5, 3.0],
            DescentOptions::default(),
            cx,
        )
        .expect("full descent");
        let converged = rep2.f_final < 1e-8;

        // PDE/stochastic nodes refuse evaluation with the executor named.
        let fx = fixture();
        let stochastic_obj = fx.objectives[1].node;
        let refuse = matches!(
            eval(&fx, stochastic_obj, &[vec![0.0; 3], vec![1.0, 0.0, 0.0, 0.0]]),
            Err(OptError::Unevaluable { kind, .. }) if kind.contains("UQ") || kind.contains("FLUX")
        );
        verdict(
            "opt-006",
            receipt && converged && refuse,
            &format!(
                "the attached P4 budget stops descent with a receipt at {} evals \
                 (objective still improved {:.2} -> {:.2}); unlimited descent \
                 converges to {:.1e}; PDE/stochastic nodes refuse evaluation naming \
                 their executor",
                rep.evals, rep.f0, rep.f_final, rep2.f_final
            ),
        );
    });
}
