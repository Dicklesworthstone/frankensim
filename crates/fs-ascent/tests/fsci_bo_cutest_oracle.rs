//! Bead 7tv.21.3: FrankenScipy oracle battery expansion — BO-standard and
//! CUTEst-class shared fixtures plus a seeded global-parity check.
//!
//! Protocol (inherited from `fsci_oracle.rs`): multimodal fixtures compare
//! from BASIN-DISCIPLINED shared starts (basin choice is not a parity
//! criterion); every analytic gradient is verified through
//! `fsci_opt::check_grad` before it is trusted (the §12 gradient-gate
//! crossover); tolerances are documented per fixture.

use fs_ascent::{LbfgsState, StopRule};
use fsci_opt::{
    DifferentialEvolutionOptions, MinimizeOptions, OptimizeMethod, check_grad,
    differential_evolution, minimize,
};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// Max-abs deviation between two points.
fn xdev(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max)
}

/// Run fs-ascent L-BFGS to a gradient-norm stop from one start.
fn ascend(
    f: &dyn Fn(&[f64]) -> f64,
    g: &dyn Fn(&[f64]) -> Vec<f64>,
    x0: &[f64],
    tol: f64,
) -> Vec<f64> {
    let mut fg = |x: &[f64]| -> (f64, Vec<f64>) { (f(x), g(x)) };
    let mut st = LbfgsState::new(x0, 10, &mut fg);
    st.run(&mut fg, &StopRule::GradNorm(tol), 8000);
    st.x.clone()
}

/// Run the fsci oracle from the same start.
fn oracle(f: impl Fn(&[f64]) -> f64 + Copy, x0: &[f64], method: OptimizeMethod) -> Vec<f64> {
    minimize(
        f,
        x0,
        MinimizeOptions {
            method: Some(method),
            tol: Some(1e-12),
            maxiter: Some(20_000),
            ..MinimizeOptions::default()
        },
    )
    .expect("fsci minimize runs")
    .x
}

/// Gradient-gate: the analytic gradient must agree with fsci's own
/// finite-difference check before any parity claim uses it.
fn gate_gradient(
    name: &str,
    f: impl Fn(&[f64]) -> f64 + Sync,
    g: impl Fn(&[f64]) -> Vec<f64>,
    probes: &[Vec<f64>],
) {
    for x0 in probes {
        let err = check_grad(&f, &g, x0).expect("check_grad runs");
        assert!(
            err < 1e-5,
            "{name}: analytic gradient fails fsci check_grad at {x0:?}: {err:.3e}"
        );
    }
}

// ---------------------------------------------------------------- Branin --

/// Branin-Hoo (BO standard): three global minima, f* = 0.39788735772973816.
fn branin(x: &[f64]) -> f64 {
    let (x1, x2) = (x[0], x[1]);
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = x2 - b * x1 * x1 + c * x1 - 6.0;
    inner * inner + 10.0 * (1.0 - t) * x1.cos() + 10.0
}

fn branin_grad(x: &[f64]) -> Vec<f64> {
    let (x1, x2) = (x[0], x[1]);
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = x2 - b * x1 * x1 + c * x1 - 6.0;
    vec![
        2.0 * inner * (-2.0 * b * x1 + c) - 10.0 * (1.0 - t) * x1.sin(),
        2.0 * inner,
    ]
}

#[test]
fn branin_bo_standard_parity() {
    const F_STAR: f64 = 0.397_887_357_729_738_16;
    // One shared in-basin start per named global minimum.
    let basins = [
        ([-3.5f64, 13.0], [-core::f64::consts::PI, 12.275]),
        ([3.0f64, 2.0], [core::f64::consts::PI, 2.275]),
        ([9.0f64, 3.0], [9.424_78, 2.475]),
    ];
    gate_gradient(
        "branin",
        branin,
        branin_grad,
        &basins.iter().map(|(s, _)| s.to_vec()).collect::<Vec<_>>(),
    );
    let mut worst_pair = 0.0f64;
    let mut worst_f = 0.0f64;
    for (start, literature) in &basins {
        let ours = ascend(&branin, &branin_grad, start, 1e-10);
        for method in [OptimizeMethod::Bfgs, OptimizeMethod::LBfgsB] {
            let theirs = oracle(branin, start, method);
            let dev = xdev(&ours, &theirs);
            worst_pair = worst_pair.max(dev);
            assert!(
                dev < 1e-4,
                "branin basin {literature:?} {method:?}: parity deviation {dev:.3e} \
                 (ours {ours:?} vs {theirs:?})"
            );
        }
        let fdev = (branin(&ours) - F_STAR).abs();
        worst_f = worst_f.max(fdev);
        assert!(
            fdev < 1e-8,
            "branin basin {literature:?}: f* off by {fdev:.3e}"
        );
        assert!(
            xdev(&ours, literature) < 1e-3,
            "branin: converged away from the literature minimum {literature:?}: {ours:?}"
        );
    }
    verdict(
        "7tv21-fsci-branin",
        true,
        &format!(
            "three global basins: fs-ascent L-BFGS vs fsci Bfgs+LBfgsB worst pair-dev \
             {worst_pair:.2e} (tol 1e-4); worst |f-f*| {worst_f:.2e} (tol 1e-8)"
        ),
    );
}

// -------------------------------------------------------------- Hartmann --

/// Hartmann-3 (BO standard): global minimum f* = -3.86278214782076 at
/// (0.114614, 0.555649, 0.852547).
const H3_ALPHA: [f64; 4] = [1.0, 1.2, 3.0, 3.2];
const H3_A: [[f64; 3]; 4] = [
    [3.0, 10.0, 30.0],
    [0.1, 10.0, 35.0],
    [3.0, 10.0, 30.0],
    [0.1, 10.0, 35.0],
];
const H3_P: [[f64; 3]; 4] = [
    [0.3689, 0.1170, 0.2673],
    [0.4699, 0.4387, 0.7470],
    [0.1091, 0.8732, 0.5547],
    [0.0381, 0.5743, 0.8828],
];

fn hartmann3(x: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..4 {
        let mut inner = 0.0;
        for j in 0..3 {
            let d = x[j] - H3_P[i][j];
            inner += H3_A[i][j] * d * d;
        }
        sum += H3_ALPHA[i] * (-inner).exp();
    }
    -sum
}

fn hartmann3_grad(x: &[f64]) -> Vec<f64> {
    let mut g = vec![0.0; 3];
    for i in 0..4 {
        let mut inner = 0.0;
        for j in 0..3 {
            let d = x[j] - H3_P[i][j];
            inner += H3_A[i][j] * d * d;
        }
        let e = H3_ALPHA[i] * (-inner).exp();
        for (j, gj) in g.iter_mut().enumerate() {
            *gj += e * 2.0 * H3_A[i][j] * (x[j] - H3_P[i][j]);
        }
    }
    g
}

#[test]
fn hartmann3_bo_standard_parity() {
    const F_STAR: f64 = -3.862_782_147_820_76;
    const X_STAR: [f64; 3] = [0.114_614, 0.555_649, 0.852_547];
    // Start inside the global basin (Hartmann-3 is multimodal; basin
    // discipline per the oracle protocol).
    let start = [0.2f64, 0.5, 0.8];
    gate_gradient("hartmann3", hartmann3, hartmann3_grad, &[start.to_vec()]);
    let ours = ascend(&hartmann3, &hartmann3_grad, &start, 1e-10);
    let mut worst_pair = 0.0f64;
    for method in [OptimizeMethod::Bfgs, OptimizeMethod::NelderMead] {
        let theirs = oracle(hartmann3, &start, method);
        let dev = xdev(&ours, &theirs);
        worst_pair = worst_pair.max(dev);
        assert!(
            dev < 5e-4,
            "hartmann3 {method:?}: parity deviation {dev:.3e} (ours {ours:?} vs {theirs:?})"
        );
    }
    // MEASURED FINDING kept on record: with the standard 4-decimal P
    // matrix, the true optimum of THIS analytic form is f = -3.86277978...
    // — the widely cited -3.862782147... presumes higher-precision P
    // entries, so the literature constant is only good to ~1e-5 here.
    let fdev = (hartmann3(&ours) - F_STAR).abs();
    assert!(fdev < 1e-5, "hartmann3: f* off by {fdev:.3e}");
    assert!(
        xdev(&ours, &X_STAR) < 1e-3,
        "hartmann3: converged away from the literature optimum: {ours:?}"
    );
    verdict(
        "7tv21-fsci-hartmann3",
        true,
        &format!(
            "global basin: L-BFGS vs fsci Bfgs+NelderMead worst pair-dev {worst_pair:.2e} \
             (tol 5e-4); |f-f*_literature| {fdev:.2e} (tol 1e-5: the 4-decimal P matrix \
             supports the cited constant only to ~1e-5)"
        ),
    );
}

// ---------------------------------------------------- CUTEst-class trio --

/// Himmelblau: four global minima, f* = 0.
fn himmelblau(x: &[f64]) -> f64 {
    let a = x[0] * x[0] + x[1] - 11.0;
    let b = x[0] + x[1] * x[1] - 7.0;
    a * a + b * b
}

fn himmelblau_grad(x: &[f64]) -> Vec<f64> {
    let a = x[0] * x[0] + x[1] - 11.0;
    let b = x[0] + x[1] * x[1] - 7.0;
    vec![4.0 * x[0] * a + 2.0 * b, 2.0 * a + 4.0 * x[1] * b]
}

/// Beale: f* = 0 at (3, 0.5).
fn beale(x: &[f64]) -> f64 {
    let (x1, x2) = (x[0], x[1]);
    let t1 = 1.5 - x1 + x1 * x2;
    let t2 = 2.25 - x1 + x1 * x2 * x2;
    let t3 = 2.625 - x1 + x1 * x2 * x2 * x2;
    t1 * t1 + t2 * t2 + t3 * t3
}

fn beale_grad(x: &[f64]) -> Vec<f64> {
    let (x1, x2) = (x[0], x[1]);
    let t1 = 1.5 - x1 + x1 * x2;
    let t2 = 2.25 - x1 + x1 * x2 * x2;
    let t3 = 2.625 - x1 + x1 * x2 * x2 * x2;
    vec![
        2.0 * t1 * (x2 - 1.0) + 2.0 * t2 * (x2 * x2 - 1.0) + 2.0 * t3 * (x2 * x2 * x2 - 1.0),
        2.0 * t1 * x1 + 2.0 * t2 * 2.0 * x1 * x2 + 2.0 * t3 * 3.0 * x1 * x2 * x2,
    ]
}

/// Booth: f* = 0 at (1, 3); convex quadratic.
fn booth(x: &[f64]) -> f64 {
    let a = x[0] + 2.0 * x[1] - 7.0;
    let b = 2.0 * x[0] + x[1] - 5.0;
    a * a + b * b
}

fn booth_grad(x: &[f64]) -> Vec<f64> {
    let a = x[0] + 2.0 * x[1] - 7.0;
    let b = 2.0 * x[0] + x[1] - 5.0;
    vec![2.0 * a + 4.0 * b, 4.0 * a + 2.0 * b]
}

#[test]
fn cutest_class_smooth_parity() {
    // (name, f, g, in-basin start, literature x*, x* tol)
    type Fixture = (
        &'static str,
        &'static (dyn Fn(&[f64]) -> f64 + Sync),
        &'static dyn Fn(&[f64]) -> Vec<f64>,
        [f64; 2],
        [f64; 2],
    );
    let fixtures: [Fixture; 4] = [
        (
            "himmelblau-basin-pp",
            &himmelblau,
            &himmelblau_grad,
            [2.5, 2.5],
            [3.0, 2.0],
        ),
        (
            "himmelblau-basin-mm",
            &himmelblau,
            &himmelblau_grad,
            [-3.5, -3.5],
            [-3.779_310, -3.283_186],
        ),
        ("beale", &beale, &beale_grad, [2.0, 0.0], [3.0, 0.5]),
        ("booth", &booth, &booth_grad, [0.0, 0.0], [1.0, 3.0]),
    ];
    let mut worst_pair = 0.0f64;
    let mut worst_f = 0.0f64;
    for (name, f, g, start, x_star) in &fixtures {
        gate_gradient(name, f, g, &[start.to_vec()]);
        let ours = ascend(f, g, start, 1e-10);
        for method in [OptimizeMethod::Bfgs, OptimizeMethod::LBfgsB] {
            let theirs = oracle(f, start, method);
            let dev = xdev(&ours, &theirs);
            worst_pair = worst_pair.max(dev);
            assert!(
                dev < 1e-4,
                "{name} {method:?}: parity deviation {dev:.3e} (ours {ours:?} vs {theirs:?})"
            );
        }
        let fdev = f(&ours).abs();
        worst_f = worst_f.max(fdev);
        assert!(fdev < 1e-8, "{name}: f* off by {fdev:.3e}");
        assert!(
            xdev(&ours, x_star) < 1e-3,
            "{name}: converged away from the literature minimum {x_star:?}: {ours:?}"
        );
    }
    verdict(
        "7tv21-fsci-cutest-trio",
        true,
        &format!(
            "himmelblau(2 basins)/beale/booth: L-BFGS vs fsci Bfgs+LBfgsB worst pair-dev \
             {worst_pair:.2e} (tol 1e-4); worst |f-f*| {worst_f:.2e} (tol 1e-8)"
        ),
    );
}

// --------------------------------------------------------- global parity --

/// Rastrigin (n=2): separable multimodal standard, f* = 0 at the origin.
fn rastrigin2(x: &[f64]) -> f64 {
    20.0 + x
        .iter()
        .map(|&v| v * v - 10.0 * (2.0 * core::f64::consts::PI * v).cos())
        .sum::<f64>()
}

fn rastrigin2_grad(x: &[f64]) -> Vec<f64> {
    x.iter()
        .map(|&v| 2.0 * v + 20.0 * core::f64::consts::PI * (2.0 * core::f64::consts::PI * v).sin())
        .collect()
}

#[test]
fn rastrigin_global_de_vs_multistart_parity() {
    // Seeded fsci differential_evolution as the global oracle vs the best
    // of a deterministic 5x5 multistart grid of fs-ascent L-BFGS.
    // Documented tolerance: DE is stochastic-but-seeded; both sides must
    // land the origin cell within 1e-3 in x and 1e-6 in f.
    gate_gradient(
        "rastrigin2",
        rastrigin2,
        rastrigin2_grad,
        &[vec![0.1, -0.1], vec![1.1, 0.9]],
    );
    let de = differential_evolution(
        rastrigin2,
        &[(-5.12, 5.12), (-5.12, 5.12)],
        DifferentialEvolutionOptions {
            seed: Some(42),
            ..DifferentialEvolutionOptions::default()
        },
    )
    .expect("fsci DE runs");

    let mut best_x = vec![f64::INFINITY; 2];
    let mut best_f = f64::INFINITY;
    for i in 0..5 {
        for j in 0..5 {
            let start = [-4.0 + 2.0 * f64::from(i), -4.0 + 2.0 * f64::from(j)];
            let x = ascend(&rastrigin2, &rastrigin2_grad, &start, 1e-10);
            let f = rastrigin2(&x);
            if f < best_f {
                best_f = f;
                best_x = x;
            }
        }
    }
    let de_f = de.fun.expect("DE reports an objective");
    let dev = xdev(&best_x, &de.x);
    assert!(
        best_f.abs() < 1e-6,
        "multistart L-BFGS missed the global optimum: f={best_f:.3e} at {best_x:?}"
    );
    assert!(
        de_f.abs() < 1e-6,
        "seeded DE missed the global optimum: f={de_f:.3e} at {:?}",
        de.x
    );
    assert!(
        dev < 1e-3,
        "global parity: multistart {best_x:?} vs DE {:?} deviate {dev:.3e}",
        de.x
    );
    verdict(
        "7tv21-fsci-rastrigin-global",
        true,
        &format!(
            "seeded DE (seed 42) and 25-start L-BFGS grid agree on the origin: \
             xdev {dev:.2e} (tol 1e-3), f_de {de_f:.2e}, f_ms {best_f:.2e} (tol 1e-6)"
        ),
    );
}
