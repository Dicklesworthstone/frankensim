//! Per-mode convergence regression (kept from the tfz.6 probe that
//! separated pre-asymptotics from bugs): both the (1,1,1) and the
//! mixed (2,1,3) Laplacian eigenmodes must converge at order → r+1
//! on m ≥ 2 ladders at r = 2 — the diagnosis that pinned the MMS
//! battery's ladder policy (single-cell parity superconvergence and
//! coarse-mesh pre-asymptotics are metric traps, not method bugs).
use fs_feec::highorder::hex::{TensorSpace, pcg_matfree};

const SUITE: &str = "fs-feec/ho-probe";
const FIXED_INPUT_SEED: u64 = 0;

fn measurement(identity: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, identity);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: "mode-sweep".to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("mode-sweep measurement must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("mode-sweep measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn finite_json(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "null".to_string()
    }
}

fn optional_finite_json(value: Option<f64>) -> String {
    value.map_or_else(|| "null".to_string(), finite_json)
}

fn solve_err(
    m: usize,
    r: usize,
    u_exact: &dyn Fn([f64; 3]) -> f64,
    f_exact: &dyn Fn([f64; 3]) -> f64,
) -> f64 {
    let sp = TensorSpace::new(m, r);
    let b = sp.load(&|p| f_exact(p));
    let mask = sp.interior_mask();
    let diag = sp.stiffness_diagonal();
    let mut bm = b;
    for (bi, &mk) in bm.iter_mut().zip(&mask) {
        if !mk {
            *bi = 0.0;
        }
    }
    let mut x = vec![0.0f64; sp.ndof()];
    let (it, conv) = pcg_matfree(
        &|v| sp.apply_stiffness(v),
        &bm,
        &mut x,
        &mask,
        &diag,
        1e-13,
        40_000,
    );
    assert!(conv, "pcg failed m={m} r={r} it={it}");
    sp.l2_error(&x, &|p| u_exact(p))
}

#[test]
fn per_mode_orders_reach_asymptotics() {
    let pi = std::f64::consts::PI;
    // Mode A: (1,1,1); Mode B: (2,1,3).
    let ua = move |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let fa = move |p: [f64; 3]| 3.0 * pi * pi * ua(p);
    let ub =
        move |p: [f64; 3]| (2.0 * pi * p[0]).sin() * (pi * p[1]).sin() * (3.0 * pi * p[2]).sin();
    let fb = move |p: [f64; 3]| 14.0 * pi * pi * ub(p);
    for (mode_id, name, u, f) in [
        (
            "a-1-1-1",
            "A(1,1,1)",
            &ua as &dyn Fn([f64; 3]) -> f64,
            &fa as &dyn Fn([f64; 3]) -> f64,
        ),
        ("b-2-1-3", "B(2,1,3)", &ub, &fb),
    ] {
        let mut prev: Option<f64> = None;
        for m in [2usize, 4, 8] {
            let e = solve_err(m, 2, u, f);
            let slope = prev.map(|p: f64| (p / e).ln() / (2.0f64).ln());
            let detail = format!("{name} r=2 m={m} L2={e:.4e} slope={slope:?}");
            measurement(
                &format!("mode-sweep/{mode_id}/measurement/m-{m}"),
                format!(
                    "{{\"detail\":\"{detail}\",\"mode\":\"{name}\",\"r\":2,\"m\":{m},\
                     \"l2_error\":{},\"slope\":{},\"input_seed\":{FIXED_INPUT_SEED},\
                     \"execution_seed\":null}}",
                    finite_json(e),
                    optional_finite_json(slope)
                ),
            );
            if m == 8 {
                let s = slope.expect("ladder ran");
                assert!(s > 2.6, "mode {name}: asymptotic slope {s:.2} below gate");
            }
            prev = Some(e);
        }
    }
}
