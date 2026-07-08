//! fs-vpm — vortex particle method (2-D core). Layer: L3.
//!
//! Vorticity, not velocity, is the natural state variable for a wake: it is
//! compact (nonzero only where the fluid actually rotates) and it advects with
//! the flow. This v0 is the inviscid 2-D core — point vortices carrying
//! circulation, inducing velocity on each other by a DESINGULARIZED
//! BIOT–SAVART kernel, advanced with RK4.
//!
//! The load-bearing checks are analytic: a single vortex of circulation `Γ`
//! induces a purely TANGENTIAL field of magnitude `Γ/(2πr)` (and does not move
//! itself); a counter-rotating PAIR separated by `d` self-propels in a straight
//! line at exactly `Γ/(2πd)` (the 2-D analog of a vortex ring's translation);
//! and the total circulation is CONSERVED. Deterministic; no dependencies.

use core::f64::consts::PI;

/// A point vortex: a position carrying a circulation `Γ`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VortexParticle {
    /// The position.
    pub pos: [f64; 2],
    /// The circulation `Γ` (signed strength).
    pub circulation: f64,
}

impl VortexParticle {
    /// A vortex particle.
    #[must_use]
    pub fn new(pos: [f64; 2], circulation: f64) -> VortexParticle {
        VortexParticle { pos, circulation }
    }
}

/// The velocity induced at `point` by all `particles` via the desingularized
/// 2-D Biot–Savart kernel `u = Σ (Γⱼ/2π)·perp(r)/(|r|² + core²)`, with
/// `perp([a,b]) = [−b, a]`. A `core` of `0` skips coincident particles.
#[must_use]
pub fn induced_velocity(particles: &[VortexParticle], point: [f64; 2], core: f64) -> [f64; 2] {
    let mut u = [0.0, 0.0];
    let core2 = core * core;
    for p in particles {
        let r = [point[0] - p.pos[0], point[1] - p.pos[1]];
        let r2 = r[0] * r[0] + r[1] * r[1] + core2;
        if r2 <= 1e-30 {
            continue; // coincident (self) with no desingularization
        }
        let g = p.circulation / (2.0 * PI * r2);
        u[0] += -g * r[1];
        u[1] += g * r[0];
    }
    u
}

fn velocities(particles: &[VortexParticle], core: f64) -> Vec<[f64; 2]> {
    particles
        .iter()
        .map(|p| induced_velocity(particles, p.pos, core))
        .collect()
}

fn displaced(particles: &[VortexParticle], vel: &[[f64; 2]], dt: f64) -> Vec<VortexParticle> {
    particles
        .iter()
        .zip(vel)
        .map(|(p, v)| {
            VortexParticle::new([p.pos[0] + dt * v[0], p.pos[1] + dt * v[1]], p.circulation)
        })
        .collect()
}

/// Advance the vortex system by one RK4 step of size `dt` (circulations are
/// invariant under inviscid advection).
#[must_use]
pub fn advect(particles: &[VortexParticle], dt: f64, core: f64) -> Vec<VortexParticle> {
    let k1 = velocities(particles, core);
    let s2 = displaced(particles, &k1, 0.5 * dt);
    let k2 = velocities(&s2, core);
    let s3 = displaced(particles, &k2, 0.5 * dt);
    let k3 = velocities(&s3, core);
    let s4 = displaced(particles, &k3, dt);
    let k4 = velocities(&s4, core);
    particles
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let vx = k1[i][0] + 2.0 * k2[i][0] + 2.0 * k3[i][0] + k4[i][0];
            let vy = k1[i][1] + 2.0 * k2[i][1] + 2.0 * k3[i][1] + k4[i][1];
            VortexParticle::new(
                [p.pos[0] + dt / 6.0 * vx, p.pos[1] + dt / 6.0 * vy],
                p.circulation,
            )
        })
        .collect()
}

/// Advance the system by `steps` RK4 steps.
#[must_use]
pub fn simulate(
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
) -> Vec<VortexParticle> {
    let mut state = particles.to_vec();
    for _ in 0..steps {
        state = advect(&state, dt, core);
    }
    state
}

/// The total circulation `Σ Γᵢ` (a conserved invariant).
#[must_use]
pub fn total_circulation(particles: &[VortexParticle]) -> f64 {
    particles.iter().map(|p| p.circulation).sum()
}

/// The centroid of vorticity `Σ Γᵢ xᵢ / Σ Γᵢ` (the "center of vorticity";
/// invariant for an isolated inviscid system). Returns `None` when the total
/// circulation is zero (e.g. a symmetric pair — use per-particle tracking).
#[must_use]
pub fn vorticity_centroid(particles: &[VortexParticle]) -> Option<[f64; 2]> {
    let total = total_circulation(particles);
    if total.abs() <= 1e-30 {
        return None;
    }
    let mut c = [0.0, 0.0];
    for p in particles {
        c[0] += p.circulation * p.pos[0];
        c[1] += p.circulation * p.pos[1];
    }
    Some([c[0] / total, c[1] / total])
}
