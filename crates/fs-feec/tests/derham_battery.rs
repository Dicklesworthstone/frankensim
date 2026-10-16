//! Tensor de Rham battery (tfz.6 slice 2): curl∘grad = div∘curl = 0
//! to machine cancellation, exact-sequence dimension bookkeeping, the
//! COMMUTING DIAGRAM (1D d∘π_C = π_D∘d and its 3D tensor version on
//! product fields), G1 projection-convergence ladders for both 1D
//! factor families at r = 1..6 (these drive the four 3D space types'
//! rates by tensorization), Legendre-mass diagonality, and the golden
//! hash.

use fs_feec::highorder::derham::TensorDeRham;
use fs_rand::StreamKey;

const SUITE: &str = "fs-feec/derham";
const FIXED_INPUT_SEED: u64 = 0;
const STREAM_INPUT_SEED: u64 = 8;
const STREAM_KERNEL: u32 = 0xDE4A;

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
    fs_obs::lint_failure_record(&event).expect("de Rham verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("de Rham verdict must use the fs-obs wire schema");
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
    fs_obs::lint_failure_record(&event).expect("de Rham measurement must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("de Rham measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: STREAM_INPUT_SEED,
        kernel: STREAM_KERNEL,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

#[test]
fn dd_vanishes_to_machine_cancellation() {
    for &(m, r) in &[(2usize, 1usize), (2, 3), (1, 5), (3, 2)] {
        let dr = TensorDeRham::new(m, r);
        let [ns, ne, ..] = dr.dims();
        // Scale of intermediate entries: G entries are O(m·r), dofs
        // O(1), so second derivatives reach O((m·r)²); the residual
        // must sit at ε relative to that.
        let scale = (m * r * m * r) as f64;
        let scalar_tile = 100 + u32::try_from(10 * m + r).expect("small");
        let s = rand_vec(ns, scalar_tile);
        let cg = dr.curl(&dr.grad(&s));
        let worst = cg.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        assert!(
            worst < 1e-13 * scale,
            "m={m} r={r}: curl grad = {worst:.3e} (scale {scale})"
        );
        let edge_tile = 200 + u32::try_from(10 * m + r).expect("small");
        let e = rand_vec(ne, edge_tile);
        let dc = dr.div(&dr.curl(&e));
        let worst_dc = dc.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        assert!(
            worst_dc < 1e-13 * scale,
            "m={m} r={r}: div curl = {worst_dc:.3e}"
        );
        verdict(
            &format!("dd-highorder/m{m}-r{r}"),
            &format!(
                "m={m} r={r} cg={worst:.1e} dc={worst_dc:.1e}; input_root=8 \
                 kernel=0xDE4A scalar_tile={scalar_tile} edge_tile={edge_tile}"
            ),
            STREAM_INPUT_SEED,
        );
    }
}

#[test]
fn exact_sequence_dimensions() {
    // χ = dim S − dim E + dim F − dim W = 1 (the cube's Betti numbers)
    // at every (m, r); plus the raw counts against the closed forms.
    for m in 1..=3usize {
        for r in 1..=6usize {
            let dr = TensorDeRham::new(m, r);
            let [s, e, f, w] = dr.dims();
            let (c, d) = (m * r + 1, m * r);
            assert_eq!(s, c * c * c);
            assert_eq!(e, 3 * d * c * c);
            assert_eq!(f, 3 * c * d * d);
            assert_eq!(w, d * d * d);
            let chi = i64::try_from(s).expect("small") - i64::try_from(e).expect("small")
                + i64::try_from(f).expect("small")
                - i64::try_from(w).expect("small");
            assert_eq!(chi, 1, "m={m} r={r}: Euler characteristic {chi}");
        }
    }
    verdict("dims", "chi = 1 for m=1..3, r=1..6", FIXED_INPUT_SEED);
}

#[test]
fn commuting_diagram_1d_and_3d() {
    // 1D: G·π_C(f) == π_D(f′) — the derivative-moment construction
    // makes this hold by design; verifying it numerically pins both
    // the G entries and the projection scalings.
    let f1 = |x: f64| (2.5 * x).sin() + 0.3 * x * x;
    let df1 = |x: f64| 2.5 * (2.5 * x).cos() + 0.6 * x;
    for &(m, r) in &[(2usize, 2usize), (3, 4), (2, 6)] {
        let dr = TensorDeRham::new(m, r);
        let pc = dr.project_c_1d(&f1, &df1);
        let pd = dr.project_d_1d(&df1);
        // G·pc
        let mut gpc = vec![0.0f64; dr.nd];
        for (i, out) in gpc.iter_mut().enumerate() {
            let mut acc = 0.0f64;
            for (j, pcj) in pc.iter().enumerate() {
                acc = dr.g.a[i * dr.nc + j].mul_add(*pcj, acc);
            }
            *out = acc;
        }
        let worst = gpc
            .iter()
            .zip(&pd)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(
            worst < 1e-11,
            "m={m} r={r}: 1D commuting residual {worst:.3e}"
        );
        verdict(
            &format!("commute-1d/m{m}-r{r}"),
            &format!("m={m} r={r} res={worst:.1e}"),
            FIXED_INPUT_SEED,
        );
    }
    // 3D on a product field f(x)g(y)h(z): grad_x(Π_S f) must equal
    // Π_{E_x}(∂_x f); tensor dofs of product functions factor, so both
    // sides are computable from 1D projections alone.
    let dr = TensorDeRham::new(2, 3);
    let (fx, dfx) = (|x: f64| (1.7 * x).sin(), |x: f64| 1.7 * (1.7 * x).cos());
    let (fy, dfy) = (
        |y: f64| 1.0 / (1.0 + y),
        |y: f64| -1.0 / ((1.0 + y) * (1.0 + y)),
    );
    let (fz, dfz) = (|z: f64| z * z + 0.5, |z: f64| 2.0 * z);
    let (pcx, pcy, pcz) = (
        dr.project_c_1d(&fx, &dfx),
        dr.project_c_1d(&fy, &dfy),
        dr.project_c_1d(&fz, &dfz),
    );
    let pdx = dr.project_d_1d(&dfx);
    // Π_S f as a tensor of 1D dofs.
    let (c, d) = (dr.nc, dr.nd);
    let mut s_dofs = vec![0.0f64; c * c * c];
    for i in 0..c {
        for j in 0..c {
            for k in 0..c {
                s_dofs[(i * c + j) * c + k] = pcx[i] * pcy[j] * pcz[k];
            }
        }
    }
    let grad = dr.grad(&s_dofs);
    // Reference E_x block = pdx ⊗ pcy ⊗ pcz.
    let mut worst = 0.0f64;
    for i in 0..d {
        for j in 0..c {
            for k in 0..c {
                let reference = pdx[i] * pcy[j] * pcz[k];
                let got = grad[(i * c + j) * c + k];
                worst = worst.max((got - reference).abs());
            }
        }
    }
    assert!(worst < 1e-11, "3D commuting residual {worst:.3e}");
    verdict("commute-3d", &format!("res={worst:.1e}"), FIXED_INPUT_SEED);
}

#[test]
fn projection_convergence_both_families() {
    // G1 ladders for the two 1D factor families — these rates drive
    // all four 3D tensor spaces (C-factor: O(h^{r+1}); D-factor:
    // O(h^r)). Full 3D vector MMS solves join the solver-stack lane.
    let f = |x: f64| (2.2 * x).sin() * (1.0 + 0.5 * x);
    let df = |x: f64| 2.2 * (2.2 * x).cos() * (1.0 + 0.5 * x) + 0.5 * (2.2 * x).sin();
    for r in 1..=6usize {
        let ladder = [4usize, 8];
        let mut errs_c = Vec::new();
        let mut errs_d = Vec::new();
        for &m in &ladder {
            let dr = TensorDeRham::new(m, r);
            let pc = dr.project_c_1d(&f, &df);
            errs_c.push(dr.l2_error_c_1d(&pc, &f));
            let pd = dr.project_d_1d(&f);
            errs_d.push(dr.l2_error_d_1d(&pd, &f));
        }
        let order_c = (errs_c[0] / errs_c[1]).ln() / 2.0f64.ln();
        let order_d = (errs_d[0] / errs_d[1]).ln() / 2.0f64.ln();
        assert!(
            order_c > r as f64 + 0.6,
            "r={r}: C-family order {order_c:.2} (errors {errs_c:?})"
        );
        assert!(
            order_d > r as f64 - 0.4,
            "r={r}: D-family order {order_d:.2} (errors {errs_d:?})"
        );
        verdict(
            &format!("proj-orders/r{r}"),
            &format!(
                "r={r} C={order_c:.2} (gate {}) D={order_d:.2} (gate {})",
                r + 1,
                r
            ),
            FIXED_INPUT_SEED,
        );
    }
}

#[test]
fn legendre_masses_are_diagonal_and_positive() {
    let dr = TensorDeRham::new(3, 4);
    assert_eq!(dr.mass_d.len(), dr.nd);
    assert!(dr.mass_d.iter().all(|&v| v > 0.0));
    // Spot-check the closed form against quadrature: ∫ L_j² over a
    // cell of width h is (h/2)·2/(2j+1).
    let h = 1.0 / 3.0;
    for j in 0..4usize {
        let expect = (h / 2.0) * 2.0 / (2.0 * j as f64 + 1.0);
        assert!(
            (dr.mass_d[j] - expect).abs() < 1e-15,
            "mass_d[{j}] = {} vs {expect}",
            dr.mass_d[j]
        );
    }
    verdict(
        "legendre-mass",
        "diagonal, positive, closed form",
        FIXED_INPUT_SEED,
    );
}

const GOLDEN_HASH: u64 = 0x14f9_beb9_cec8_4078; // recorded at tfz.6 slice 2, frozen

#[test]
fn derham_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let dr = TensorDeRham::new(2, 4);
    for v in dr.g.a.iter().filter(|v| **v != 0.0) {
        feed(*v);
    }
    for v in &dr.mass_d {
        feed(*v);
    }
    let [ns, ne, nf, _] = dr.dims();
    let s = rand_vec(ns, 300);
    for v in dr.grad(&s).iter().step_by(7) {
        feed(*v);
    }
    let e = rand_vec(ne, 301);
    for v in dr.curl(&e).iter().step_by(11) {
        feed(*v);
    }
    let f = rand_vec(nf, 302);
    for v in dr.div(&f).iter().step_by(5) {
        feed(*v);
    }
    measurement(
        "derham-golden/measurement",
        "derham-golden",
        format!(
            "{{\"input_root\":8,\"kernel\":\"0xDE4A\",\"tiles\":{{\"scalar\":300,\
             \"edge\":301,\"face\":302}},\"actual_hash\":\"{acc:#018x}\",\
             \"expected_hash\":\"{GOLDEN_HASH:#018x}\"}}"
        ),
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "derham bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
