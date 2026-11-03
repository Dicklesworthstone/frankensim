//! Symmetric convection of the existing clockwise-positive wake.
//!
//! Each mutual pair uses the original VPM algebraic-core flow. Its first and
//! quadratic circulation moments are preserved without a posteriori fixes.
//! Freestream translation and the fixed quarter-chord vortex are separate
//! exact subflows; the latter exchanges impulse with the supported airfoil.
//! This does not change the empirical bound-circulation or shedding closure.

use super::WakeVortex;
use fs_vpm::VortexParticle;
use fs_vpm::pairwise::{PairFlowError, advance_pair};

/// Evolve a private candidate. The caller retains ownership of publication
/// because a refused later pair may leave this candidate partially advanced.
pub(super) fn advance(
    wake: &mut [WakeVortex],
    freestream: [f64; 2],
    bound: f64,
    quarter: [f64; 2],
    dt: f64,
    core: f64,
) -> Result<(), PairFlowError> {
    translate(wake, freestream, 0.5 * dt)?;
    bound_flow(wake, bound, quarter, 0.5 * dt, core)?;
    for i in 0..wake.len() {
        for j in i + 1..wake.len() {
            pair_flow(wake, i, j, 0.5 * dt, core)?;
        }
    }
    for i in (0..wake.len()).rev() {
        for j in (i + 1..wake.len()).rev() {
            pair_flow(wake, i, j, 0.5 * dt, core)?;
        }
    }
    bound_flow(wake, bound, quarter, 0.5 * dt, core)?;
    translate(wake, freestream, 0.5 * dt)
}

fn pair_flow(
    wake: &mut [WakeVortex],
    i: usize,
    j: usize,
    dt: f64,
    core: f64,
) -> Result<(), PairFlowError> {
    // BEM circulation is clockwise-positive; VPM uses counterclockwise.
    let next = advance_pair(
        [
            VortexParticle::new(wake[i].pos, -wake[i].gamma),
            VortexParticle::new(wake[j].pos, -wake[j].gamma),
        ],
        dt,
        core,
    )?;
    wake[i].pos = next[0].pos;
    wake[j].pos = next[1].pos;
    Ok(())
}

fn bound_flow(
    wake: &mut [WakeVortex],
    bound: f64,
    quarter: [f64; 2],
    dt: f64,
    core: f64,
) -> Result<(), PairFlowError> {
    if bound == 0.0 {
        return Ok(());
    }
    let fixed = VortexParticle::new(quarter, -bound);
    for vortex in wake {
        // A zero-strength tracer leaves the bound source fixed and follows
        // its exact circular streamline at the existing regularized speed.
        vortex.pos = advance_pair([VortexParticle::new(vortex.pos, 0.0), fixed], dt, core)?[0].pos;
    }
    Ok(())
}

fn translate(wake: &mut [WakeVortex], velocity: [f64; 2], dt: f64) -> Result<(), PairFlowError> {
    for vortex in wake {
        for axis in 0..2 {
            vortex.pos[axis] += velocity[axis] * dt;
            if !vortex.pos[axis].is_finite() {
                return Err(PairFlowError::Unrepresentable);
            }
        }
    }
    Ok(())
}
