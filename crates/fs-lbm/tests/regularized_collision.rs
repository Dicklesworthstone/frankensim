//! CollisionModel2::Regularized battery (bead frankensim-3zkcr): the
//! Latt–Chopard regularized collision must (1) agree with BGK
//! analytically on states whose non-equilibrium already lies in the
//! second-order Hermite subspace, (2) keep the exact shear-viscosity
//! law nu = (tau - 1/2)/3, (3) deliver the capability BGK cannot —
//! surviving high cell Reynolds number (the executed 9ok02 boundary:
//! plain BGK destabilizes at Re_cell = u/nu in (36, 48) on the
//! jet-labium family) — and (4) stay bitwise deterministic. The BGK
//! blow-up control in (3) doubles as the mode-discrimination
//! mutation catcher: if `collide_into` ignored `Grid::collision`,
//! the regularized run would destabilize identically.

use fs_lbm::core2::{Cell, CollisionModel2, Grid};
use fs_lbm::{CS2, Q, equilibrium};

const E: [(i32, i32); Q] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];
const W: [f64; Q] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// A population set whose non-equilibrium is EXACTLY the
/// second-order Hermite term `w_i (9/2) (Q_i : A)` for symmetric A:
/// zero mass/momentum perturbation (Sum w Q = 0, odd third moment),
/// so the collide path recovers `feq` bit-for-bit from the moments
/// and the regularized projection must reproduce the perturbation
/// identically (D2Q9's fourth velocity moment is exactly isotropic,
/// so Sum w Q_ab Q_cd = cs^4 (delta delta + delta delta)).
fn hermite_state(rho: f64, u: [f64; 2], a: [[f64; 2]; 2]) -> [f64; Q] {
    let mut f = equilibrium(rho, u[0], u[1]);
    for q in 0..Q {
        let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
        let qa = (ex * ex - CS2) * a[0][0] + (ey * ey - CS2) * a[1][1] + 2.0 * ex * ey * a[0][1];
        f[q] += 4.5 * W[q] * qa;
    }
    f
}

/// (1) Analytic BGK equivalence on Hermite states: both operators
/// applied to the same Hermite-subspace state must agree to
/// arithmetic roundoff (the derivation is exact; only the operation
/// order differs).
#[test]
fn regularized_equals_bgk_on_hermite_states() {
    let cases = [
        (1.0, [0.0, 0.0], [[3.0e-4, -1.5e-4], [-1.5e-4, -2.0e-4]]),
        (1.02, [0.05, -0.03], [[-2.0e-4, 1.0e-4], [1.0e-4, 4.0e-4]]),
        // With a trace component (Sum w Q_ab Q_cd has NO delta_ab
        // delta_cd term, so the trace part must still round-trip).
        (0.98, [0.02, 0.06], [[5.0e-4, 0.0], [0.0, 5.0e-4]]),
    ];
    for (rho, u, a) in cases {
        let mut bgk = Grid::uniform(4, 4, 0.62);
        let mut reg = Grid::uniform(4, 4, 0.62);
        reg.collision = CollisionModel2::Regularized;
        let state = hermite_state(rho, u, a);
        for i in 0..16 {
            bgk.f[i] = state;
            reg.f[i] = state;
        }
        let mut post_bgk = Vec::new();
        let mut post_reg = Vec::new();
        bgk.collide_into(&mut post_bgk);
        reg.collide_into(&mut post_reg);
        for (pb, pr) in post_bgk.iter().zip(&post_reg) {
            for q in 0..Q {
                let scale = pb[q].abs().max(1e-3);
                assert!(
                    (pb[q] - pr[q]).abs() / scale < 1e-13,
                    "Hermite-state equivalence violated at q={q}: bgk {} vs reg {}",
                    pb[q],
                    pr[q]
                );
            }
        }
    }
}

/// Project the ux field onto sin(2 pi y / ny) (single-k shear-wave
/// amplitude).
fn shear_amplitude(grid: &Grid) -> f64 {
    let ny = grid.ny;
    let mut acc = 0.0;
    for y in 0..ny {
        let m = grid.moments(grid.idx(0, y));
        let phase = 2.0 * std::f64::consts::PI * (y as f64) / (ny as f64);
        acc += m.u[0] * phase.sin();
    }
    2.0 * acc / ny as f64
}

/// (2) Viscosity identity: a periodic single-k shear wave decays at
/// exp(-nu k^2 t); the measured nu under BOTH operators must equal
/// (tau - 1/2)/3 (BGK is the control proving the fixture itself).
#[test]
fn shear_wave_viscosity_identity_both_modes() {
    for mode in [CollisionModel2::Bgk, CollisionModel2::Regularized] {
        let (nx, ny) = (8, 64);
        let tau = 0.56; // nu = 0.02
        let mut grid = Grid::uniform(nx, ny, tau);
        grid.collision = mode;
        let amp = 0.01;
        for y in 0..ny {
            let phase = 2.0 * std::f64::consts::PI * (y as f64) / (ny as f64);
            let ux = amp * phase.sin();
            for x in 0..nx {
                let i = grid.idx(x, y);
                grid.f[i] = equilibrium(1.0, ux, 0.0);
            }
        }
        let mut scratch = Vec::new();
        let (t0, t1) = (500usize, 3500usize);
        for _ in 0..t0 {
            grid.step(&mut scratch);
        }
        let a0 = shear_amplitude(&grid);
        for _ in t0..t1 {
            grid.step(&mut scratch);
        }
        let a1 = shear_amplitude(&grid);
        let k = 2.0 * std::f64::consts::PI / ny as f64;
        let nu_measured = (a0 / a1).ln() / (k * k * (t1 - t0) as f64);
        let nu_expected = (tau - 0.5) / 3.0;
        let dev = (nu_measured - nu_expected).abs() / nu_expected;
        println!(
            "{{\"suite\":\"fs-lbm\",\"case\":\"shear-viscosity\",\"mode\":\"{mode:?}\",\"nu_measured\":{nu_measured:.6},\"nu_expected\":{nu_expected:.6},\"deviation\":{dev:.5}}}"
        );
        assert!(
            dev < 0.02,
            "{mode:?}: measured nu {nu_measured:.6} vs (tau-1/2)/3 = {nu_expected:.6} ({dev:.3})"
        );
    }
}

/// Doubly periodic double shear layer (Minion–Brown style) at cell
/// Reynolds u/nu = 100 — far beyond the executed plain-BGK boundary
/// (36, 48).
fn shear_layer_grid(mode: CollisionModel2) -> Grid {
    let (nx, ny) = (96, 96);
    let u0 = 0.1;
    let tau = 0.503; // nu = 0.001 -> Re_cell = 100
    let mut grid = Grid::uniform(nx, ny, tau);
    grid.collision = mode;
    for y in 0..ny {
        let yf = y as f64 / ny as f64;
        let ux = if yf <= 0.5 {
            u0 * (80.0 * (yf - 0.25)).tanh()
        } else {
            u0 * (80.0 * (0.75 - yf)).tanh()
        };
        for x in 0..nx {
            let xf = x as f64 / nx as f64;
            let uy = 0.05 * u0 * (2.0 * std::f64::consts::PI * xf).sin();
            let i = grid.idx(x, y);
            grid.f[i] = equilibrium(1.0, ux, uy);
        }
    }
    grid
}

fn run_finite(mut grid: Grid, steps: usize) -> bool {
    let result = std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
        let mut scratch = Vec::new();
        for _ in 0..steps {
            grid.step(&mut scratch);
        }
        grid.f
            .iter()
            .zip(&grid.flags)
            .filter(|&(_, &fl)| fl == Cell::Fluid)
            .all(|(f, _)| f.iter().all(|q| q.is_finite()))
    }));
    result.unwrap_or(false)
}

/// (3) The capability claim + mode-discrimination control:
/// Regularized survives Re_cell = 100 where plain BGK destabilizes.
/// If BGK ever survives this config, the config no longer
/// discriminates the collision mode — strengthen it rather than
/// deleting the control.
#[test]
fn regularized_survives_cell_reynolds_100_where_bgk_dies() {
    let steps = 3000;
    let reg_ok = run_finite(shear_layer_grid(CollisionModel2::Regularized), steps);
    let bgk_ok = run_finite(shear_layer_grid(CollisionModel2::Bgk), steps);
    println!(
        "{{\"suite\":\"fs-lbm\",\"case\":\"shear-layer-re-cell-100\",\"regularized_finite\":{reg_ok},\"bgk_finite\":{bgk_ok}}}"
    );
    assert!(
        reg_ok,
        "Regularized destabilized at Re_cell 100 — the capability this mode exists for"
    );
    assert!(
        !bgk_ok,
        "plain BGK survived Re_cell 100: this config no longer discriminates the \
         collision mode — strengthen the control (do not delete it)"
    );
}

/// (4) Bitwise determinism of the regularized path.
#[test]
fn regularized_determinism_bitwise() {
    let run = || {
        let mut grid = shear_layer_grid(CollisionModel2::Regularized);
        let mut scratch = Vec::new();
        for _ in 0..50 {
            grid.step(&mut scratch);
        }
        grid.f
    };
    let (a, b) = (run(), run());
    for (fa, fb) in a.iter().zip(&b) {
        for q in 0..Q {
            assert_eq!(fa[q].to_bits(), fb[q].to_bits());
        }
    }
}
