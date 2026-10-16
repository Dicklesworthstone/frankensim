//! CollisionModel2 battery (bead frankensim-3zkcr): the D2Q9 core's
//! collision operators beyond plain BGK. Every claim here is an
//! EXECUTED result, not folklore — the 2026-08-10 stability scans
//! falsified two textbook expectations on this rig family (plain
//! projected regularization is WORSE than BGK on unresolved shear
//! layers, and TRT's "magic = 1/4" choice reduces the stable window
//! at tau near 1/2), so the pins below encode what actually
//! happened:
//!
//! - Projected regularized (`Regularized`) equals BGK analytically
//!   on second-order Hermite states and keeps the viscosity law,
//!   but DIES on an about-one-cell shear layer at Re_cell = 100
//!   where BGK survives (recorded limitation, and the
//!   mode-discrimination control: a collide path that ignored
//!   `Grid::collision` could not reproduce the asymmetry).
//! - Recursive regularized (`RecursiveRegularized`, RR3) repairs
//!   exactly that fragility and survives alongside BGK.
//! - TRT at the BGK-equivalent magic `(tau - 1/2)^2` reproduces BGK
//!   to roundoff; its viscosity law holds at magic = 1/4.
//! - Central moment (`CentralMoment`) survives everything in the
//!   scan, keeps the viscosity law, and is the operator that
//!   extended the jet-labium rig beyond Re_cell 190 (fs-aeroac
//!   probe) — the capability this bead exists for.
//! - All modes conserve mass and momentum exactly and are bitwise
//!   deterministic.

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

const ALL_MODES: [CollisionModel2; 5] = [
    CollisionModel2::Bgk,
    CollisionModel2::Regularized,
    CollisionModel2::RecursiveRegularized,
    CollisionModel2::Trt { magic: 0.25 },
    CollisionModel2::CentralMoment,
];

/// A population set whose non-equilibrium is EXACTLY the
/// second-order Hermite term `w_i (9/2) (Q_i : A)` for symmetric A:
/// zero mass/momentum perturbation (Sum w Q = 0, odd third moment),
/// so the collide path recovers `feq` from the moments and the
/// regularized projection must reproduce the perturbation
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

fn collide_uniform_state(mode: CollisionModel2, tau: f64, state: [f64; Q]) -> Vec<[f64; Q]> {
    let mut grid = Grid::uniform(4, 4, tau);
    grid.collision = mode;
    for i in 0..16 {
        grid.f[i] = state;
    }
    let mut post = Vec::new();
    grid.collide_into(&mut post);
    post
}

/// Analytic BGK equivalence of the PROJECTED regularized operator on
/// Hermite-subspace states (exact derivation; only operation order
/// differs). At u = 0 the recursive third-order terms vanish, so
/// RecursiveRegularized joins the identity there.
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
        let state = hermite_state(rho, u, a);
        let post_bgk = collide_uniform_state(CollisionModel2::Bgk, 0.62, state);
        let post_reg = collide_uniform_state(CollisionModel2::Regularized, 0.62, state);
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
        if u == [0.0, 0.0] {
            let post_rr3 =
                collide_uniform_state(CollisionModel2::RecursiveRegularized, 0.62, state);
            for (pb, pr) in post_bgk.iter().zip(&post_rr3) {
                for q in 0..Q {
                    assert!(
                        (pb[q] - pr[q]).abs() / pb[q].abs().max(1e-3) < 1e-13,
                        "RR3 must equal BGK at zero velocity"
                    );
                }
            }
        }
    }
}

/// TRT at the BGK-equivalent magic `(tau - 1/2)^2` (tau_odd = tau)
/// must reproduce BGK to roundoff on a generic state.
#[test]
fn trt_at_bgk_magic_equals_bgk() {
    let tau = 0.62;
    let state = hermite_state(1.01, [0.04, -0.02], [[2.0e-4, 1.0e-4], [1.0e-4, -3.0e-4]]);
    let post_bgk = collide_uniform_state(CollisionModel2::Bgk, tau, state);
    let magic = (tau - 0.5) * (tau - 0.5);
    let post_trt = collide_uniform_state(CollisionModel2::Trt { magic }, tau, state);
    for (pb, pt) in post_bgk.iter().zip(&post_trt) {
        for q in 0..Q {
            assert!(
                (pb[q] - pt[q]).abs() / pb[q].abs().max(1e-3) < 1e-13,
                "TRT at the BGK magic diverged from BGK at q={q}"
            );
        }
    }
}

/// Exact mass and momentum conservation of every operator on a
/// generic (non-equilibrium) state — the cheap catcher for
/// reconstruction bugs in the moment-space operators.
#[test]
fn all_modes_conserve_mass_and_momentum() {
    let state = hermite_state(1.03, [0.06, -0.04], [[4.0e-4, -2.0e-4], [-2.0e-4, 1.0e-4]]);
    let (m0, p0) = moments_of(&state);
    for mode in ALL_MODES {
        let post = collide_uniform_state(mode, 0.55, state);
        let (m1, p1) = moments_of(&post[0]);
        assert!(
            (m1 - m0).abs() < 1e-14,
            "{mode:?} mass drift {}",
            (m1 - m0).abs()
        );
        for d in 0..2 {
            assert!(
                (p1[d] - p0[d]).abs() < 1e-14,
                "{mode:?} momentum drift {}",
                (p1[d] - p0[d]).abs()
            );
        }
    }
}

fn moments_of(f: &[f64; Q]) -> (f64, [f64; 2]) {
    let mut mass = 0.0;
    let mut mom = [0.0; 2];
    for q in 0..Q {
        mass += f[q];
        mom[0] += f64::from(E[q].0) * f[q];
        mom[1] += f64::from(E[q].1) * f[q];
    }
    (mass, mom)
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

/// Viscosity identity: a periodic single-k shear wave decays at
/// exp(-nu k^2 t); the measured nu under EVERY operator must equal
/// (tau - 1/2)/3 (BGK is the control proving the fixture itself;
/// executed deviations were all under 0.1%).
#[test]
fn shear_wave_viscosity_identity_all_modes() {
    for mode in ALL_MODES {
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

/// Doubly periodic double shear layer (Minion–Brown style); `sharp`
/// controls the layer width (80 -> about 1.2 cells: unresolved).
fn shear_layer_grid(mode: CollisionModel2, sharp: f64, tau: f64, u0: f64) -> Grid {
    let (nx, ny) = (96, 96);
    let mut grid = Grid::uniform(nx, ny, tau);
    grid.collision = mode;
    for y in 0..ny {
        let yf = y as f64 / ny as f64;
        let ux = if yf <= 0.5 {
            u0 * (sharp * (yf - 0.25)).tanh()
        } else {
            u0 * (sharp * (0.75 - yf)).tanh()
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

/// The EXECUTED stability asymmetry on the unresolved shear layer at
/// Re_cell = 100 (executed 2026-08-10; example driver `stab_scan`):
/// BGK, RR3, and CentralMoment survive; plain PROJECTED
/// regularization destabilizes — its recorded limitation, and the
/// mode-discrimination control (a collide path ignoring
/// `Grid::collision` cannot reproduce the asymmetry). If the reg
/// column ever flips, the operator changed — re-derive, do not
/// delete.
#[test]
fn unresolved_shear_layer_stability_asymmetry() {
    let steps = 3000;
    let (sharp, tau, u0) = (80.0, 0.503, 0.1); // Re_cell = 100
    let bgk = run_finite(
        shear_layer_grid(CollisionModel2::Bgk, sharp, tau, u0),
        steps,
    );
    let reg = run_finite(
        shear_layer_grid(CollisionModel2::Regularized, sharp, tau, u0),
        steps,
    );
    let rr3 = run_finite(
        shear_layer_grid(CollisionModel2::RecursiveRegularized, sharp, tau, u0),
        steps,
    );
    let cm = run_finite(
        shear_layer_grid(CollisionModel2::CentralMoment, sharp, tau, u0),
        steps,
    );
    println!(
        "{{\"suite\":\"fs-lbm\",\"case\":\"shear-layer-re-cell-100\",\"bgk\":{bgk},\"reg\":{reg},\"rr3\":{rr3},\"cm\":{cm}}}"
    );
    assert!(bgk, "BGK lost the smooth-box baseline it held on record");
    assert!(
        rr3,
        "RecursiveRegularized destabilized at Re_cell 100 — the repair it exists for"
    );
    assert!(
        cm,
        "CentralMoment destabilized at Re_cell 100 — the capability operator regressed"
    );
    assert!(
        !reg,
        "plain projected Regularized survived the unresolved layer: the recorded \
         asymmetry (and mode-discrimination control) vanished — re-derive the pins"
    );
}

/// Bitwise determinism of every non-BGK path.
#[test]
fn new_modes_determinism_bitwise() {
    for mode in [
        CollisionModel2::Regularized,
        CollisionModel2::RecursiveRegularized,
        CollisionModel2::Trt { magic: 0.25 },
        CollisionModel2::CentralMoment,
    ] {
        let run = || {
            let mut grid = shear_layer_grid(mode, 20.0, 0.51, 0.05);
            let mut scratch = Vec::new();
            for _ in 0..50 {
                grid.step(&mut scratch);
            }
            grid.f
        };
        let (a, b) = (run(), run());
        for (fa, fb) in a.iter().zip(&b) {
            for q in 0..Q {
                assert_eq!(fa[q].to_bits(), fb[q].to_bits(), "{mode:?}");
            }
        }
    }
}
