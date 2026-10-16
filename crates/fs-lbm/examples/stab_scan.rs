//! Collision-operator stability scan (bead 3zkcr) — the executed
//! evidence driver behind the CollisionModel2 landscape recorded in
//! CONTRACT.md and tests/regularized_collision.rs: TRT-vs-BGK
//! equivalence at the BGK magic, viscosity-role checks, and
//! survive/destabilize verdicts on impulsively started double shear
//! layers across sharpness, tau, and operator. Deterministic;
//! re-run to re-derive the recorded landscape.
use fs_lbm::core2::{Cell, CollisionModel2, Grid};
use fs_lbm::equilibrium;

fn grid(mode: CollisionModel2, sharp: f64, tau: f64, u0: f64) -> Grid {
    let (nx, ny) = (96, 96);
    let mut g = Grid::uniform(nx, ny, tau);
    g.collision = mode;
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
            let i = g.idx(x, y);
            g.f[i] = equilibrium(1.0, ux, uy);
        }
    }
    g
}

fn survives(mut g: Grid, steps: usize) -> bool {
    std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
        let mut s = Vec::new();
        for _ in 0..steps {
            g.step(&mut s);
        }
        g.f.iter()
            .zip(&g.flags)
            .filter(|&(_, &fl)| fl == Cell::Fluid)
            .all(|(f, _)| f.iter().all(|q| q.is_finite()))
    }))
    .unwrap_or(false)
}

fn main() {
    // 1. implementation sanity: TRT at magic = (tau-1/2)^2 must be BGK.
    let tau = 0.62;
    let mut a = grid(CollisionModel2::Bgk, 20.0, tau, 0.05);
    let mut b = grid(
        CollisionModel2::Trt {
            magic: (tau - 0.5) * (tau - 0.5),
        },
        20.0,
        tau,
        0.05,
    );
    let mut s = Vec::new();
    for _ in 0..50 {
        a.step(&mut s);
        b.step(&mut s);
    }
    let max_dev =
        a.f.iter()
            .zip(&b.f)
            .flat_map(|(x, y)| x.iter().zip(y).map(|(p, q)| (p - q).abs()))
            .fold(0.0f64, f64::max);
    println!("bgk-equivalence max_dev={max_dev:.3e}");
    // 1b. viscosity role check under TRT (flipped even/odd would track tau_odd).
    for mode in [
        CollisionModel2::Bgk,
        CollisionModel2::Trt { magic: 0.25 },
        CollisionModel2::CentralMoment,
    ] {
        let (nx, ny) = (8usize, 64usize);
        let tau = 0.56;
        let mut g = Grid::uniform(nx, ny, tau);
        g.collision = mode;
        for y in 0..ny {
            let ph = 2.0 * std::f64::consts::PI * (y as f64) / (ny as f64);
            for x in 0..nx {
                let i = g.idx(x, y);
                g.f[i] = equilibrium(1.0, 0.01 * ph.sin(), 0.0);
            }
        }
        let mut sc = Vec::new();
        let amp = |g: &Grid| -> f64 {
            let mut acc = 0.0;
            for y in 0..ny {
                let m = g.moments(g.idx(0, y));
                acc += m.u[0] * (2.0 * std::f64::consts::PI * (y as f64) / (ny as f64)).sin();
            }
            2.0 * acc / ny as f64
        };
        for _ in 0..500 {
            g.step(&mut sc);
        }
        let a0 = amp(&g);
        for _ in 0..3000 {
            g.step(&mut sc);
        }
        let a1 = amp(&g);
        let k = 2.0 * std::f64::consts::PI / ny as f64;
        let nu = (a0 / a1).ln() / (k * k * 3000.0);
        println!(
            "mode={mode:?} nu_measured={nu:.6} nu_expected={:.6}",
            (tau - 0.5) / 3.0
        );
    }
    // 2. magic sweep on the hard sharp-layer case and the moderate one.
    for &(sharp, tau, u0) in &[
        (80.0, 0.503, 0.1),
        (20.0, 0.5015, 0.1),
        (20.0, 0.50075, 0.1),
    ] {
        let re_cell = u0 / ((tau - 0.5) / 3.0);
        let bgk = survives(grid(CollisionModel2::Bgk, sharp, tau, u0), 3000);
        let reg = survives(grid(CollisionModel2::Regularized, sharp, tau, u0), 3000);
        let rr3 = survives(
            grid(CollisionModel2::RecursiveRegularized, sharp, tau, u0),
            3000,
        );
        let cm = survives(grid(CollisionModel2::CentralMoment, sharp, tau, u0), 3000);
        println!(
            "sharp={sharp} tau={tau} Re_cell={re_cell:.0} bgk={bgk} reg={reg} rr3={rr3} cm={cm}"
        );
    }
}
