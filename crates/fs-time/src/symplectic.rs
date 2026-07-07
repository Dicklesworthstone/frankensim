//! Symplectic integrators for separable Hamiltonians H = T(p) + V(q)
//! with T = ½‖p‖² (unit mass; scale q/p externally otherwise).
//!
//! Störmer–Verlet IS the discrete-Lagrangian (variational) integrator
//! for L(q, q̇) = ½‖q̇‖² − V(q) with midpoint quadrature (Marsden–West):
//! extremizing the discrete action Σ L_d(q_k, q_{k+1}) yields exactly
//! the leapfrog update — the equivalence is tested, not just cited.
//! Symplecticity is what buys BOUNDED energy error over 10⁶ steps where
//! RK4 drifts secularly (the acceptance demo in the battery).

use fs_ad::revolve::checkpointed_adjoint;

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
/// given ∂J/∂(q_N, p_N), returns ∂J/∂(q_0, p_0) — the adjoint OF THE
/// STEPPER, propagated backwards with fs-ad's checkpointed revolve
/// (memory O(log N) instead of O(N)). `force_jvp` computes the action
/// of ∂F/∂q at `q` on a direction `v` (exact user-supplied linearization
/// — gradcheck'd in the battery).
#[must_use]
pub fn verlet_adjoint<F, J>(
    q0: &[f64],
    p0: &[f64],
    h: f64,
    steps: usize,
    force: &F,
    force_jvp: &J,
    bar_q_n: &[f64],
    bar_p_n: &[f64],
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
        let (mut bar_q, mut bar_p) = bar;
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
        let q_new: Vec<f64> = q.iter().zip(&p_half).map(|(qi, ph)| h.mul_add(*ph, *qi)).collect();
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
        let _ = p_half;
        (bar_q, bar_p_half)
    };
    let budget = fs_ad::revolve::min_budget(steps);
    let (bar, _stats) = checkpointed_adjoint(
        &state0,
        steps,
        budget,
        &fwd,
        &rev,
        (bar_q_n.to_vec(), bar_p_n.to_vec()),
    );
    bar
}
