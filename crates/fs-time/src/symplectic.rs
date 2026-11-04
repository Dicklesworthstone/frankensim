//! Symplectic integrators for separable Hamiltonians H = T(p) + V(q)
//! with unit mass or a caller-owned constant inverse-mass operator.
//!
//! Störmer–Verlet IS the discrete-Lagrangian (variational) integrator
//! for L(q, q̇) = ½‖q̇‖² − V(q) with midpoint quadrature (Marsden–West):
//! extremizing the discrete action Σ L_d(q_k, q_{k+1}) yields exactly
//! the leapfrog update — the equivalence is tested, not just cited.
//! Symplecticity is what buys BOUNDED energy error over 10⁶ steps where
//! RK4 drifts secularly (the acceptance demo in the battery).

use fs_ad::revolve::checkpointed_adjoint;
use fs_solver::op::LinearOp;

/// One Störmer–Verlet step for q̈ = force(q): kick–drift–kick.
/// `force` writes F(q) into its output slice.
pub fn verlet_step<F: Fn(&[f64], &mut [f64])>(
    q: &mut [f64],
    p: &mut [f64],
    h: f64,
    force: &F,
    scratch: &mut [f64],
) {
    let n = q.len();
    force(q, scratch);
    for i in 0..n {
        p[i] = (0.5 * h).mul_add(scratch[i], p[i]);
    }
    for i in 0..n {
        q[i] = h.mul_add(p[i], q[i]);
    }
    force(q, scratch);
    for i in 0..n {
        p[i] = (0.5 * h).mul_add(scratch[i], p[i]);
    }
}

/// Discrete adjoint of an n-step Verlet trajectory for a terminal cost:
/// given the terminal cotangent `bar_n` = (∂J/∂q_N, ∂J/∂p_N), returns
/// (∂J/∂q_0, ∂J/∂p_0) — the adjoint OF THE STEPPER, propagated
/// backwards with fs-ad's checkpointed revolve (memory O(log N) instead
/// of O(N)). `force_jvp` computes the action of ∂F/∂q at `q` on a
/// direction `v` (exact user-supplied linearization — gradcheck'd in
/// the battery).
#[must_use]
pub fn verlet_adjoint<F, J>(
    q0: &[f64],
    p0: &[f64],
    h: f64,
    steps: usize,
    force: &F,
    force_jvp: &J,
    bar_n: (&[f64], &[f64]),
) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(&[f64], &mut [f64]),
    J: Fn(&[f64], &[f64], &mut [f64]),
{
    let n = q0.len();
    let state0: Vec<f64> = q0.iter().chain(p0.iter()).copied().collect();
    let fwd = |_i: usize, s: &Vec<f64>| -> Vec<f64> {
        let mut q = s[..n].to_vec();
        let mut p = s[n..].to_vec();
        let mut scratch = vec![0.0f64; n];
        verlet_step(&mut q, &mut p, h, force, &mut scratch);
        q.into_iter().chain(p).collect()
    };
    // Reverse of one kick–drift–kick step: transpose of the tangent map.
    // Forward maps (with F₀ = F(q₀), F₁ = F(q₁)):
    //   p_half = p + h/2·F(q)
    //   q'     = q + h·p_half
    //   p'     = p_half + h/2·F(q')
    // Adjoint (bar quantities pulled back in reverse order):
    let rev = |_i: usize, s: &Vec<f64>, bar: (Vec<f64>, Vec<f64>)| -> (Vec<f64>, Vec<f64>) {
        let (mut bar_q, bar_p) = bar;
        let q = &s[..n];
        let p = &s[n..];
        let mut scratch = vec![0.0f64; n];
        // Recompute intermediates.
        force(q, &mut scratch);
        let p_half: Vec<f64> = p
            .iter()
            .zip(&scratch)
            .map(|(pi, fi)| (0.5 * h).mul_add(*fi, *pi))
            .collect();
        let q_new: Vec<f64> = q
            .iter()
            .zip(&p_half)
            .map(|(qi, ph)| h.mul_add(*ph, *qi))
            .collect();
        // Step 3 adjoint: p' = p_half + h/2·F(q') ⇒
        //   bar_p_half += bar_p';  bar_q_new += h/2·(∂F/∂q')ᵀ·bar_p'.
        let mut jt = vec![0.0f64; n];
        force_jvp(&q_new, &bar_p, &mut jt); // symmetric ∂F/∂q assumed ⇒ Jᵀv = Jv
        for i in 0..n {
            bar_q[i] = (0.5 * h).mul_add(jt[i], bar_q[i]);
        }
        // Step 2 adjoint: q' = q + h·p_half ⇒ bar_p_half += h·bar_q'.
        let mut bar_p_half = bar_p.clone();
        for i in 0..n {
            bar_p_half[i] = h.mul_add(bar_q[i], bar_p_half[i]);
        }
        // Step 1 adjoint: p_half = p + h/2·F(q) ⇒
        //   bar_p = bar_p_half; bar_q += h/2·Jᵀ·bar_p_half.
        force_jvp(q, &bar_p_half, &mut jt);
        for i in 0..n {
            bar_q[i] = (0.5 * h).mul_add(jt[i], bar_q[i]);
        }
        (bar_q, bar_p_half)
    };
    let budget = fs_ad::revolve::min_budget(steps);
    let (bar, _stats) = checkpointed_adjoint(
        &state0,
        steps,
        budget,
        &fwd,
        &rev,
        (bar_n.0.to_vec(), bar_n.1.to_vec()),
    );
    bar
}

/// Refusal of a physical-momentum Verlet step or its discrete adjoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MomentumVerletError {
    /// State, operator, seed, workspace, or step dimensions are inadmissible.
    InvalidInput(&'static str),
    /// An input, operator action, force action, or derived value is nonfinite.
    NonFinite(&'static str),
    /// Cancellation was observed before publishing the state or gradient.
    Cancelled,
}

impl std::fmt::Display for MomentumVerletError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "momentum Verlet refused: {self:?}")
    }
}

impl std::error::Error for MomentumVerletError {}

fn momentum_poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MomentumVerletError> {
    if cancelled() {
        Err(MomentumVerletError::Cancelled)
    } else {
        Ok(())
    }
}

fn momentum_finite(
    values: &[f64],
    context: &'static str,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MomentumVerletError> {
    for chunk in values.chunks(256) {
        momentum_poll(cancelled)?;
        if chunk.iter().any(|x| !x.is_finite()) {
            return Err(MomentumVerletError::NonFinite(context));
        }
    }
    momentum_poll(cancelled)
}

fn momentum_inputs<M: LinearOp + ?Sized>(
    q: &[f64],
    p: &[f64],
    h: f64,
    inverse_mass: &M,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, MomentumVerletError> {
    momentum_poll(cancelled)?;
    if q.is_empty() || q.len() != p.len() || q.len() != inverse_mass.n() {
        return Err(MomentumVerletError::InvalidInput(
            "matching positive state and inverse-mass dimensions required",
        ));
    }
    if !h.is_finite() || h == 0.0 {
        return Err(MomentumVerletError::InvalidInput(
            "finite nonzero time step required",
        ));
    }
    let workspace = q
        .len()
        .checked_mul(3)
        .ok_or(MomentumVerletError::InvalidInput(
            "workspace length overflow",
        ))?;
    momentum_finite(q, "position", cancelled)?;
    momentum_finite(p, "momentum", cancelled)?;
    Ok(workspace)
}

fn momentum_action(
    output: &mut [f64],
    context: &'static str,
    action: impl FnOnce(&mut [f64]),
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MomentumVerletError> {
    for chunk in output.chunks_mut(256) {
        momentum_poll(cancelled)?;
        // An incompletely written action must not reuse an earlier stage.
        chunk.fill(f64::NAN);
    }
    momentum_poll(cancelled)?;
    action(output);
    momentum_finite(output, context, cancelled)
}

// Shared complete primal replay for the forward map and its pullback. The
// reverse sweep must also refuse an inadmissible endpoint force or momentum.
fn momentum_proposal<M: LinearOp + ?Sized, F: Fn(&[f64], &mut [f64])>(
    q: &[f64],
    p: &[f64],
    h: f64,
    inverse_mass: &M,
    force: &F,
    workspace: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MomentumVerletError> {
    let n = q.len();
    let (next_q, rest) = workspace.split_at_mut(n);
    let (half_p, work) = rest.split_at_mut(n);
    momentum_action(work, "initial force", |out| force(q, out), cancelled)?;
    for i in 0..n {
        if i % 256 == 0 {
            momentum_poll(cancelled)?;
        }
        half_p[i] = (0.5 * h).mul_add(work[i], p[i]);
    }
    momentum_finite(half_p, "half-step momentum", cancelled)?;
    momentum_action(
        work,
        "inverse-mass action",
        |out| inverse_mass.apply(half_p, out),
        cancelled,
    )?;
    for i in 0..n {
        if i % 256 == 0 {
            momentum_poll(cancelled)?;
        }
        next_q[i] = h.mul_add(work[i], q[i]);
    }
    momentum_finite(next_q, "drifted position", cancelled)?;
    momentum_action(work, "endpoint force", |out| force(next_q, out), cancelled)?;
    for i in 0..n {
        if i % 256 == 0 {
            momentum_poll(cancelled)?;
        }
        half_p[i] = (0.5 * h).mul_add(work[i], half_p[i]);
    }
    momentum_finite(half_p, "endpoint momentum", cancelled)
}

/// One kick–drift–kick step with physical canonical momentum `p = M velocity`.
///
/// `inverse_mass` supplies the constant action `M^-1 p`; it is never assembled
/// or inverted here. `force` writes physical generalized force, not acceleration.
/// Both actions must overwrite their complete outputs. The caller supplies at
/// least `3 * q.len()` scratch components; the step itself allocates nothing.
/// Input `q` and `p` remain unchanged on every returned refusal, while scratch
/// may change. Cancellation is polled around actions and every 256 components;
/// action implementations own their internal cancellation latency. Publication
/// is two final slice copies after the last poll.
///
/// For the Hamiltonian `H = p^T M^-1 p / 2 + V(q)`, symplecticity and the usual
/// bounded-energy interpretation require a constant symmetric positive-definite
/// inverse mass and `force = -grad V`. These properties are caller obligations,
/// not inferred or certified from an arbitrary [`LinearOp`]. Approximate mass
/// solves carry their own defects. The map is second order and reversible for
/// fixed smooth actions; negative `h` performs a reverse step. Configuration-
/// dependent mass, constraints, impulses and mass-parameter gradients are outside
/// this separable lane. The existing unit-mass [`verlet_step`] is unchanged.
pub fn verlet_momentum_step<M: LinearOp + ?Sized, F: Fn(&[f64], &mut [f64])>(
    q: &mut [f64],
    p: &mut [f64],
    h: f64,
    inverse_mass: &M,
    force: &F,
    scratch: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MomentumVerletError> {
    let required = momentum_inputs(q, p, h, inverse_mass, cancelled)?;
    if scratch.len() < required {
        return Err(MomentumVerletError::InvalidInput(
            "three state vectors of workspace required",
        ));
    }
    let workspace = &mut scratch[..required];
    momentum_proposal(q, p, h, inverse_mass, force, workspace, cancelled)?;
    let (next_q, rest) = workspace.split_at_mut(q.len());
    let (next_p, _) = rest.split_at_mut(p.len());
    q.copy_from_slice(next_q);
    p.copy_from_slice(next_p);
    Ok(())
}

/// Discrete adjoint of repeated [`verlet_momentum_step`] maps.
///
/// `force_vjp(q, seed, out)` writes `(d force / d q)^T seed`. The drift uses
/// [`LinearOp::apply_transpose`] of the SAME supplied inverse-mass action; no
/// symmetry shortcut is used. The returned vectors are the initial position
/// and physical-momentum cotangents. The inverse mass, force and derivatives
/// must remain pure and unchanged throughout replay. This differentiates the
/// supplied discrete map, not an exact trajectory or parameters of the mass.
///
/// Fallible binary checkpointing follows the existing Verlet adjoint's schedule:
/// O(n log steps) live state and O(steps log steps) forward work, with immediate
/// propagation of cancellation/refusal and no partial gradient. A zero-step
/// chain returns its seed after normal input admission. The Hamiltonian claim
/// requires the same physical assumptions as the forward step; the pullback
/// also remains valid for admitted nonsymmetric linear actions.
pub fn verlet_momentum_adjoint<M, F, J, C>(
    q0: &[f64],
    p0: &[f64],
    h: f64,
    steps: usize,
    inverse_mass: &M,
    force: &F,
    force_vjp: &J,
    bar_n: (&[f64], &[f64]),
    cancelled: &mut C,
) -> Result<(Vec<f64>, Vec<f64>), MomentumVerletError>
where
    M: LinearOp + ?Sized,
    F: Fn(&[f64], &mut [f64]),
    J: Fn(&[f64], &[f64], &mut [f64]),
    C: FnMut() -> bool,
{
    let required = momentum_inputs(q0, p0, h, inverse_mass, cancelled)?;
    let n = q0.len();
    if bar_n.0.len() != n || bar_n.1.len() != n {
        return Err(MomentumVerletError::InvalidInput(
            "matching adjoint seed required",
        ));
    }
    momentum_finite(bar_n.0, "position seed", cancelled)?;
    momentum_finite(bar_n.1, "momentum seed", cancelled)?;
    let initial: Vec<_> = q0.iter().chain(p0).copied().collect();
    let forward = |state: &[f64], cancelled: &mut C| {
        let mut next = state.to_vec();
        let (q, p) = next.split_at_mut(n);
        let mut scratch = vec![0.0; required];
        verlet_momentum_step(q, p, h, inverse_mass, force, &mut scratch, cancelled)?;
        Ok(next)
    };
    let reverse = |state: &[f64], bar: (Vec<f64>, Vec<f64>), cancelled: &mut C| {
        let (q, p) = state.split_at(n);
        momentum_inputs(q, p, h, inverse_mass, cancelled)?;
        let mut scratch = vec![0.0; required];
        momentum_proposal(q, p, h, inverse_mass, force, &mut scratch, cancelled)?;
        let (next_q, rest) = scratch.split_at_mut(n);
        let (_, work) = rest.split_at_mut(n);
        let (mut bar_q, mut bar_p) = bar;
        momentum_action(
            work,
            "endpoint force transpose",
            |out| {
                force_vjp(next_q, &bar_p, out);
            },
            cancelled,
        )?;
        for i in 0..n {
            if i % 256 == 0 {
                momentum_poll(cancelled)?;
            }
            bar_q[i] = (0.5 * h).mul_add(work[i], bar_q[i]);
        }
        momentum_finite(&bar_q, "drift position cotangent", cancelled)?;
        momentum_action(
            work,
            "inverse-mass transpose",
            |out| {
                inverse_mass.apply_transpose(&bar_q, out);
            },
            cancelled,
        )?;
        for i in 0..n {
            if i % 256 == 0 {
                momentum_poll(cancelled)?;
            }
            bar_p[i] = h.mul_add(work[i], bar_p[i]);
        }
        momentum_finite(&bar_p, "initial momentum cotangent", cancelled)?;
        momentum_action(
            work,
            "initial force transpose",
            |out| {
                force_vjp(q, &bar_p, out);
            },
            cancelled,
        )?;
        for i in 0..n {
            if i % 256 == 0 {
                momentum_poll(cancelled)?;
            }
            bar_q[i] = (0.5 * h).mul_add(work[i], bar_q[i]);
        }
        momentum_finite(&bar_q, "initial position cotangent", cancelled)?;
        Ok((bar_q, bar_p))
    };
    let seed = (bar_n.0.to_vec(), bar_n.1.to_vec());
    momentum_sweep(&initial, steps, seed, &forward, &reverse, cancelled)
}

fn momentum_sweep<F, R, C>(
    state: &[f64],
    steps: usize,
    bar: (Vec<f64>, Vec<f64>),
    forward: &F,
    reverse: &R,
    cancelled: &mut C,
) -> Result<(Vec<f64>, Vec<f64>), MomentumVerletError>
where
    F: Fn(&[f64], &mut C) -> Result<Vec<f64>, MomentumVerletError>,
    R: Fn(
        &[f64],
        (Vec<f64>, Vec<f64>),
        &mut C,
    ) -> Result<(Vec<f64>, Vec<f64>), MomentumVerletError>,
    C: FnMut() -> bool,
{
    momentum_poll(cancelled)?;
    if steps == 0 {
        return Ok(bar);
    }
    if steps == 1 {
        return reverse(state, bar, cancelled);
    }
    let mid = steps / 2 + steps % 2;
    let mut midpoint = state.to_vec();
    for _ in 0..mid {
        midpoint = forward(&midpoint, cancelled)?;
    }
    let bar = momentum_sweep(&midpoint, steps - mid, bar, forward, reverse, cancelled)?;
    drop(midpoint);
    momentum_sweep(state, mid, bar, forward, reverse, cancelled)
}
