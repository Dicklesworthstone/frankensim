//! dynamics.rs — Tier-2 dynamical-systems demos:
//!
//! * [`symplectic_vs_euler`] — a Kepler 2-body orbit integrated with BOTH a
//!   symplectic Störmer–Verlet scheme and explicit Euler, exposing Euler's
//!   energy drift against the symplectic scheme's near-conservation.
//! * [`lorenz_points`]       — the Lorenz attractor via RK4.
//! * [`ga_motor_orbit`]      — a real geometric-algebra screw motion: a PGA
//!   Cl(3,0,1) MOTOR (rotation about an axis + translation along it) applied
//!   iteratively to a seed ring, via the `fs-ga` kernel.
//!
//! Every input is clamped and every loop capped.

use fs_ga::{Motor, Point};
use fs_math::det;

/* ----------------------------------------------------------------------- */
/*  symplectic_vs_euler — energy conservation, side by side                 */
/* ----------------------------------------------------------------------- */

/// Integrate a Kepler 2-body orbit (unit mass, GM = 1) with BOTH a symplectic
/// Störmer–Verlet (velocity-Verlet) integrator and explicit forward Euler,
/// recording each trajectory and its Hamiltonian `H = ½|v|² − 1/|r|` at every
/// step. Euler's energy drifts upward and its orbit spirals out; the symplectic
/// scheme keeps `H` bounded and the ellipse closed.
///
/// Output layout (length `1 + 6*steps`):
/// - `[0]`                        = `steps` (as an f64).
/// - `[1 .. 1+2*steps]`           = symplectic positions, interleaved `x,y`.
/// - `[1+2*steps .. 1+4*steps]`   = Euler positions, interleaved `x,y`.
/// - `[1+4*steps .. 1+5*steps]`   = symplectic energy `H` per step.
/// - `[1+5*steps .. 1+6*steps]`   = Euler energy `H` per step.
///
/// Step `k` (0-based) records the state BEFORE the `k`-th update, so index 0 is
/// the shared initial condition. `steps` clamped to `[2,4000]`, `dt` to
/// `[1e-4, 0.2]`.
pub fn symplectic_vs_euler(steps_in: usize, dt_in: f64) -> Vec<f64> {
    let steps = steps_in.clamp(2, 4000);
    let dt = dt_in.clamp(1.0e-4, 0.2);

    let accel = |x: f64, y: f64| -> (f64, f64) {
        let r2 = x * x + y * y;
        let r = r2.sqrt().max(1.0e-6);
        let inv_r3 = 1.0 / (r2 * r);
        (-x * inv_r3, -y * inv_r3)
    };
    let energy = |x: f64, y: f64, vx: f64, vy: f64| -> f64 {
        let r = (x * x + y * y).sqrt().max(1.0e-9);
        0.5 * (vx * vx + vy * vy) - 1.0 / r
    };

    // Shared elliptical initial condition.
    let (mut sx, mut sy, mut svx, mut svy) = (1.0f64, 0.0f64, 0.0f64, 0.8f64);
    let (mut ex, mut ey, mut evx, mut evy) = (1.0f64, 0.0f64, 0.0f64, 0.8f64);

    let mut symp_pos = Vec::with_capacity(2 * steps);
    let mut eul_pos = Vec::with_capacity(2 * steps);
    let mut symp_e = Vec::with_capacity(steps);
    let mut eul_e = Vec::with_capacity(steps);

    for _ in 0..steps {
        symp_pos.push(sx);
        symp_pos.push(sy);
        eul_pos.push(ex);
        eul_pos.push(ey);
        symp_e.push(energy(sx, sy, svx, svy));
        eul_e.push(energy(ex, ey, evx, evy));

        // Symplectic velocity-Verlet (a half-kick / drift / half-kick).
        let (ax, ay) = accel(sx, sy);
        svx += 0.5 * dt * ax;
        svy += 0.5 * dt * ay;
        sx += dt * svx;
        sy += dt * svy;
        let (ax2, ay2) = accel(sx, sy);
        svx += 0.5 * dt * ax2;
        svy += 0.5 * dt * ay2;

        // Explicit forward Euler (position from old velocity).
        let (eax, eay) = accel(ex, ey);
        ex += dt * evx;
        ey += dt * evy;
        evx += dt * eax;
        evy += dt * eay;
    }

    let mut out = Vec::with_capacity(1 + 6 * steps);
    out.push(steps as f64);
    out.extend_from_slice(&symp_pos);
    out.extend_from_slice(&eul_pos);
    out.extend_from_slice(&symp_e);
    out.extend_from_slice(&eul_e);
    out
}

/* ----------------------------------------------------------------------- */
/*  lorenz_points — deterministic chaos via RK4                             */
/* ----------------------------------------------------------------------- */

/// Integrate the Lorenz system (`σ=10`, `β=8/3`, and the supplied `ρ`) with a
/// fixed-step RK4, returning the 3D trajectory.
///
/// Output layout: `steps * 3` values, interleaved `x,y,z` per step; step 0 is
/// the initial condition `(0.1, 0, 0)`. Length `3*steps`.
///
/// `steps` clamped to `[1,200000]`, `dt` to `[1e-4, 0.02]`, `ρ` to `[0,200]`.
pub fn lorenz_points(steps_in: usize, dt_in: f64, rho_in: f64) -> Vec<f64> {
    let steps = steps_in.clamp(1, 200_000);
    let dt = dt_in.clamp(1.0e-4, 0.02);
    let rho = rho_in.clamp(0.0, 200.0);
    let sigma = 10.0f64;
    let beta = 8.0 / 3.0;
    let f = |s: [f64; 3]| -> [f64; 3] {
        [
            sigma * (s[1] - s[0]),
            s[0] * (rho - s[2]) - s[1],
            s[0] * s[1] - beta * s[2],
        ]
    };
    let mut s = [0.1f64, 0.0, 0.0];
    let mut out = Vec::with_capacity(3 * steps);
    for _ in 0..steps {
        out.push(s[0]);
        out.push(s[1]);
        out.push(s[2]);
        let k1 = f(s);
        let s2 = [s[0] + 0.5 * dt * k1[0], s[1] + 0.5 * dt * k1[1], s[2] + 0.5 * dt * k1[2]];
        let k2 = f(s2);
        let s3 = [s[0] + 0.5 * dt * k2[0], s[1] + 0.5 * dt * k2[1], s[2] + 0.5 * dt * k2[2]];
        let k3 = f(s3);
        let s4 = [s[0] + dt * k3[0], s[1] + dt * k3[1], s[2] + dt * k3[2]];
        let k4 = f(s4);
        for d in 0..3 {
            s[d] += dt / 6.0 * (k1[d] + 2.0 * k2[d] + 2.0 * k3[d] + k4[d]);
        }
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  ga_motor_orbit — real PGA screw motion (fs-ga Motor)                    */
/* ----------------------------------------------------------------------- */

/// Real geometric-algebra rigid motion: build a PGA Cl(3,0,1) **screw motor**
/// (rotation about the z-axis composed with translation along it) with `fs-ga`,
/// then apply its powers to a seed ring of points, sweeping out a helical coil.
///
/// The seed is a circle of radius 0.35 offset one unit from the screw axis, so
/// the screw sweeps it into a spring/slinky. Frame `s` (0-based) holds the
/// motor's `s`-th power applied to the seed (frame 0 is the seed itself).
///
/// Output layout: `[nPoints, steps, then steps*nPoints*3 xyz]`.
/// - `[0]` = `nPoints`, `[1]` = `steps`.
/// - then `steps * nPoints * 3` values, frame-major: for each frame `s`, for
///   each seed point `i`, its transformed `x,y,z`. Total `2 + steps*nPoints*3`.
///
/// `nPoints` clamped to `[3,200]`, `steps` to `[1,240]`.
pub fn ga_motor_orbit(npoints_in: usize, steps_in: usize) -> Vec<f64> {
    let npoints = npoints_in.clamp(3, 200);
    let steps = steps_in.clamp(1, 240);
    let two_pi = 2.0 * std::f64::consts::PI;
    let seed_r = 0.35f64;
    let offset = 1.0f64;

    // Seed ring, offset from the z screw-axis and standing in the radial–z plane.
    let mut seed: Vec<Point> = Vec::with_capacity(npoints);
    for i in 0..npoints {
        let a = two_pi * i as f64 / npoints as f64;
        seed.push(Point::new(offset + seed_r * det::cos(a), 0.0, seed_r * det::sin(a)));
    }

    // Screw motor: rotate about z by dθ AND translate along z by dz per step.
    let dtheta = two_pi / 60.0; // one full turn every 60 steps
    let dz = 0.06f64;
    let step_motor = Motor::translator(0.0, 0.0, dz).compose(&Motor::rotor([0.0, 0.0, 1.0], dtheta));

    let mut out = Vec::with_capacity(2 + steps * npoints * 3);
    out.push(npoints as f64);
    out.push(steps as f64);
    let mut acc = Motor::identity();
    for _s in 0..steps {
        for &pt in seed.iter() {
            match acc.transform_point(pt) {
                Ok(q) => {
                    out.push(q.x);
                    out.push(q.y);
                    out.push(q.z);
                }
                Err(_) => {
                    out.push(pt.x);
                    out.push(pt.y);
                    out.push(pt.z);
                }
            }
        }
        acc = acc.compose(&step_motor);
    }
    out
}
