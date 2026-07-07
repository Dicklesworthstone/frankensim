//! Embedded-pair adaptivity: Dormand–Prince RK45 with a PI step-size
//! controller (smooth step sequences, deterministic rejection handling)
//! and a RESUMABLE state machine — checkpoint = clone, and split runs
//! are bitwise-equal to straight runs (the P7 obligation, tested).

/// PI controller settings (standard exponents).
#[derive(Debug, Clone)]
pub struct PiController {
    /// Proportional exponent (default 0.7/5).
    pub k_p: f64,
    /// Integral exponent (default 0.4/5).
    pub k_i: f64,
    /// Safety factor.
    pub safety: f64,
    /// Step growth clamp.
    pub max_growth: f64,
    /// Step shrink clamp.
    pub max_shrink: f64,
}

impl Default for PiController {
    fn default() -> PiController {
        PiController { k_p: 0.14, k_i: 0.08, safety: 0.9, max_growth: 5.0, max_shrink: 0.2 }
    }
}

/// Resumable integration state (plain data; `clone()` IS a checkpoint).
#[derive(Debug, Clone)]
pub struct AdaptiveState {
    /// Current time.
    pub t: f64,
    /// Current solution.
    pub u: Vec<f64>,
    /// Current step size.
    pub h: f64,
    /// Previous error ratio (the PI controller's integral memory).
    pub err_prev: f64,
    /// Accepted steps so far.
    pub accepted: usize,
    /// Rejected steps so far.
    pub rejected: usize,
}

impl AdaptiveState {
    /// Fresh state.
    #[must_use]
    pub fn new(t0: f64, u0: &[f64], h0: f64) -> AdaptiveState {
        AdaptiveState { t: t0, u: u0.to_vec(), h: h0, err_prev: 1.0, accepted: 0, rejected: 0 }
    }
}

/// Dormand–Prince 5(4) coefficients.
const A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
const C: [f64; 6] = [0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const B5: [f64; 7] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
    0.0,
];
const B4: [f64; 7] = [
    5179.0 / 57600.0,
    0.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    -92097.0 / 339_200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
];

/// Advance the state until `t_end` (or `max_steps` attempts) with the
/// PI controller at relative tolerance `rtol` (plus absolute `atol`).
/// Deterministic; resumable mid-flight.
pub fn rk45_adaptive<F: Fn(f64, &[f64], &mut [f64])>(
    state: &mut AdaptiveState,
    rhs: &F,
    t_end: f64,
    rtol: f64,
    atol: f64,
    pi: &PiController,
    max_steps: usize,
) {
    let n = state.u.len();
    let mut k: Vec<Vec<f64>> = vec![vec![0.0; n]; 7];
    let mut attempts = 0usize;
    while state.t < t_end && attempts < max_steps {
        attempts += 1;
        let h = state.h.min(t_end - state.t);
        rhs(state.t, &state.u, &mut k[0]);
        for stage in 0..6 {
            let mut u_stage = state.u.clone();
            for (j, kj) in k.iter().enumerate().take(stage + 1) {
                let a = A[stage][j];
                if a != 0.0 {
                    for i in 0..n {
                        u_stage[i] = (h * a).mul_add(kj[i], u_stage[i]);
                    }
                }
            }
            let (head, tail) = k.split_at_mut(stage + 1);
            let _ = head;
            rhs(state.t + C[stage] * h, &u_stage, &mut tail[0]);
        }
        // 5th-order solution + embedded error estimate.
        let mut u5 = state.u.clone();
        let mut err = 0.0f64;
        for i in 0..n {
            let mut du5 = 0.0f64;
            let mut du4 = 0.0f64;
            for (j, kj) in k.iter().enumerate() {
                du5 = B5[j].mul_add(kj[i], du5);
                du4 = B4[j].mul_add(kj[i], du4);
            }
            u5[i] = h.mul_add(du5, u5[i]);
            let scale = atol + rtol * state.u[i].abs().max(u5[i].abs());
            let e = h * (du5 - du4) / scale;
            err = err.max(e.abs());
        }
        let err = err.max(1e-300);
        if err <= 1.0 {
            // Accept; PI step update.
            state.t += h;
            state.u = u5;
            state.accepted += 1;
            let factor = pi.safety
                * fs_math::det::pow(err, -pi.k_p)
                * fs_math::det::pow(state.err_prev, pi.k_i);
            state.h = h * factor.clamp(pi.max_shrink, pi.max_growth);
            state.err_prev = err;
        } else {
            // Reject; shrink.
            state.rejected += 1;
            let factor = pi.safety * fs_math::det::pow(err, -0.2);
            state.h = h * factor.clamp(pi.max_shrink, 1.0);
        }
    }
}
