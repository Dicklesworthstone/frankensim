//! 6-DOF integrator core (bead wf-root-guzez.4.2.1, E3.2-i). The fixed
//! 120 Hz simulation spine: translation by velocity-Verlet with a force
//! callback, rotation by a Strang (kick–free–kick) split around
//! `fs_time::lie::rigid_body_step` — the torque-free Lie-group step
//! sandwiched between two half-tick body-moment kicks. Both halves are
//! second order; the composition is the DECLARED ORDER 2 of the spine,
//! measured (not asserted) by the battery's Richardson fixture.
//!
//! Determinism: pure f64 arithmetic + the det-routed Lie step; a per-tick
//! state digest (fs-blake3) seeds the E3.5 structured-checkpoint program.
//! The partitioned multi-rate integrator with time-scale certificates is
//! E3.2b and does NOT live here.

use crate::Refusal;
use fs_blake3::hash_domain;
use fs_time::lie::rigid_body_step;

/// The fixed simulation tick rate [Hz] (plan §4.1).
pub const TICK_HZ: f64 = 120.0;
/// Admitted timestep domain [s].
pub const MIN_DT_S: f64 = 1.0e-6;
/// Upper timestep bound [s].
pub const MAX_DT_S: f64 = 0.1;
/// Step-budget cap per advance call.
pub const MAX_STEPS: u32 = 1_000_000;
/// Identity domain for spine tick digests.
pub const TICK_DIGEST_DOMAIN: &str = "org.frankensim.fs-flyer.spine-tick.v1";

/// Full 6-DOF rigid state (frd-body-v1 / NED world, SI units).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SixDofState {
    /// World position [m], NED.
    pub pos_m: [f64; 3],
    /// World velocity [m/s], NED.
    pub vel_mps: [f64; 3],
    /// Unit quaternion (w, x, y, z), body→world.
    pub quat: [f64; 4],
    /// Body-frame angular velocity [rad/s].
    pub omega_body: [f64; 3],
}

/// External generalized loads at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Loads {
    /// World-frame force [N].
    pub force_n: [f64; 3],
    /// Body-frame moment [N·m].
    pub moment_nm: [f64; 3],
}

/// Rigid-body mass properties (diagonal principal inertia).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidBody {
    /// Mass [kg].
    pub mass_kg: f64,
    /// Principal inertia [kg·m²].
    pub inertia_kgm2: [f64; 3],
}

impl RigidBody {
    /// Validate mass properties.
    ///
    /// # Errors
    /// `non-finite-input`, `mass-outside-domain`, `inertia-outside-domain`.
    pub fn admit(&self) -> Result<(), Refusal> {
        if !self.mass_kg.is_finite() || !self.inertia_kgm2.iter().all(|v| v.is_finite()) {
            return Err(Refusal {
                code: "non-finite-input",
                message: "mass properties must be finite".into(),
                ranked_repairs: vec!["check the mass build-up output".into()],
            });
        }
        if self.mass_kg <= 0.0 {
            return Err(Refusal {
                code: "mass-outside-domain",
                message: format!("mass {} kg must be positive", self.mass_kg),
                ranked_repairs: vec!["use FlyerDesign::mass_build_up for the gross mass".into()],
            });
        }
        if self.inertia_kgm2.iter().any(|&i| i <= 0.0) {
            return Err(Refusal {
                code: "inertia-outside-domain",
                message: format!("principal inertia {:?} must be positive", self.inertia_kgm2),
                ranked_repairs: vec!["use the build-up inertias".into()],
            });
        }
        Ok(())
    }
}

fn admit_step(state: &SixDofState, dt_s: f64) -> Result<(), Refusal> {
    let finite = state.pos_m.iter().all(|v| v.is_finite())
        && state.vel_mps.iter().all(|v| v.is_finite())
        && state.quat.iter().all(|v| v.is_finite())
        && state.omega_body.iter().all(|v| v.is_finite())
        && dt_s.is_finite();
    if !finite {
        return Err(Refusal {
            code: "non-finite-input",
            message: "state and dt must be finite".into(),
            ranked_repairs: vec!["a NaN state upstream means a load model diverged".into()],
        });
    }
    if !(MIN_DT_S..=MAX_DT_S).contains(&dt_s) {
        return Err(Refusal {
            code: "timestep-outside-domain",
            message: format!("dt {dt_s:e} outside admitted [{MIN_DT_S:e}, {MAX_DT_S:e}]"),
            ranked_repairs: vec![format!("the spine tick is 1/{TICK_HZ} s")],
        });
    }
    Ok(())
}

/// One spine tick: velocity-Verlet translation (loads sampled at t and
/// t+dt via the callback) + Strang kick–free–kick rotation with the
/// body moment held over the tick.
///
/// The callback receives (time_s, &state) and returns the loads; it is
/// called exactly twice per tick (begin, end-predictor) — deterministic
/// call order is part of the contract.
///
/// # Errors
/// Admission refusals ([`RigidBody::admit`], state/dt checks).
pub fn step<F>(
    body: &RigidBody,
    state: &SixDofState,
    t_s: f64,
    dt_s: f64,
    mut loads: F,
) -> Result<SixDofState, Refusal>
where
    F: FnMut(f64, &SixDofState) -> Loads,
{
    body.admit()?;
    admit_step(state, dt_s)?;
    let l0 = loads(t_s, state);
    let inv_m = 1.0 / body.mass_kg;
    let a0 = [
        l0.force_n[0] * inv_m,
        l0.force_n[1] * inv_m,
        l0.force_n[2] * inv_m,
    ];
    // Drift: x+ = x + v·dt + a0·dt²/2.
    let mut pos = *state;
    for i in 0..3 {
        pos.pos_m[i] += state.vel_mps[i] * dt_s + 0.5 * a0[i] * dt_s * dt_s;
    }
    // Rotation: half-kick, torque-free Lie step, half-kick.
    let half = 0.5 * dt_s;
    let mut w = state.omega_body;
    for i in 0..3 {
        w[i] += l0.moment_nm[i] / body.inertia_kgm2[i] * half;
    }
    let (q_new, w_free) = rigid_body_step(state.quat, w, body.inertia_kgm2, dt_s);
    pos.quat = q_new;
    // End-of-tick loads at the predicted state (positions + rotated frame).
    let mut predictor = pos;
    for i in 0..3 {
        predictor.vel_mps[i] = state.vel_mps[i] + a0[i] * dt_s;
        predictor.omega_body[i] = w_free[i];
    }
    let l1 = loads(t_s + dt_s, &predictor);
    let a1 = [
        l1.force_n[0] * inv_m,
        l1.force_n[1] * inv_m,
        l1.force_n[2] * inv_m,
    ];
    for i in 0..3 {
        pos.vel_mps[i] = state.vel_mps[i] + 0.5 * (a0[i] + a1[i]) * dt_s;
        pos.omega_body[i] = w_free[i] + l1.moment_nm[i] / body.inertia_kgm2[i] * half;
    }
    Ok(pos)
}

/// Advance `steps` ticks from `state` at the fixed tick; returns the final
/// state and the per-tick digest trace (the E3.5 seed).
///
/// # Errors
/// `step-budget-exceeded` above [`MAX_STEPS`] (tested at cap AND cap+1);
/// per-tick refusals pass through.
pub fn advance<F>(
    body: &RigidBody,
    state: &SixDofState,
    t0_s: f64,
    dt_s: f64,
    steps: u32,
    mut loads: F,
) -> Result<(SixDofState, Vec<String>), Refusal>
where
    F: FnMut(f64, &SixDofState) -> Loads,
{
    if steps > MAX_STEPS {
        return Err(Refusal {
            code: "step-budget-exceeded",
            message: format!("steps {steps} exceeds cap {MAX_STEPS}"),
            ranked_repairs: vec!["advance long runs in bounded chunks".into()],
        });
    }
    let mut s = *state;
    let mut digests = Vec::with_capacity(steps as usize);
    for k in 0..steps {
        let t = t0_s + f64::from(k) * dt_s;
        s = step(body, &s, t, dt_s, &mut loads)?;
        digests.push(tick_digest(k, &s));
    }
    Ok((s, digests))
}

/// Content digest of one tick's full state (exact bits, versioned domain).
#[must_use]
pub fn tick_digest(tick: u32, state: &SixDofState) -> String {
    let mut payload = Vec::with_capacity(8 * 13 + 4);
    payload.extend_from_slice(&tick.to_le_bytes());
    for v in state
        .pos_m
        .iter()
        .chain(state.vel_mps.iter())
        .chain(state.quat.iter())
        .chain(state.omega_body.iter())
    {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    hash_domain(TICK_DIGEST_DOMAIN, &payload).to_hex()
}
