//! SE(3) Lie-group and variational integrator lanes on PGA motors
//! (bead `frankensim-ext-time-se3-lanes-3ol0`).
//!
//! Three obligations from the bead, each with an honest claim scope:
//!
//! 1. **Exp-map lane**: motor states updated by the fs-ga screw
//!    exponential, with the double cover (`M ~ −M`) fixed by a
//!    deterministic canonicalization and drift-controlled
//!    renormalization that returns receipts instead of silently
//!    patching the state.
//! 2. **Variational lane**: a discrete Euler–Poincaré step for the
//!    free rigid body in body-momentum form. Spatial angular momentum
//!    is conserved EXACTLY by construction (the update transports the
//!    momentum by the same group element that updates the attitude);
//!    energy behavior earns the conservative-theorem claim class only
//!    for declared smooth conservative fixtures at fixed step with a
//!    converged solve — everything else gets a measured balance
//!    receipt, never the theorem.
//! 3. **Discrete adjoint**: derived from the ACTUAL fixed-point
//!    residual of the variational step via the implicit-function
//!    theorem and verified against finite differences of the whole
//!    map. The 3×3 residual Jacobians are formed by central
//!    differences of that residual (a stated v1 boundary; analytic
//!    tangents are follow-up work).
//!
//! RATTLE-style constraint projection is exposed as a hook trait for
//! fs-mbd; the constrained lanes live there, not here.

use crate::lie::{quat_exp, quat_mul, quat_rotate};
use fs_ga::pga::{EVEN_BLADES, axis_bivector, ideal_bivector};
use fs_ga::{Motor, exp_bivector};
use fs_math::det;

/// Typed refusals for the SE(3) lanes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Se3Error {
    /// A state or parameter component is NaN or infinite.
    NonFinite,
    /// The motor has no usable even component to anchor the
    /// double-cover sign (all exactly zero).
    DegenerateMotor,
    /// A principal inertia or the step size is not positive.
    InvalidParameter,
    /// The variational fixed-point solve did not converge.
    SolverDiverged {
        /// Iterations spent.
        iters: u32,
        /// Final residual norm.
        residual: f64,
    },
    /// The IFT linear solve met a numerically singular Jacobian.
    SingularJacobian,
}

impl std::fmt::Display for Se3Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Se3Error::NonFinite => write!(f, "non-finite state or parameter"),
            Se3Error::DegenerateMotor => {
                write!(f, "motor has no nonzero even component to anchor the sign")
            }
            Se3Error::InvalidParameter => {
                write!(f, "step size and principal inertias must be positive and finite")
            }
            Se3Error::SolverDiverged { iters, residual } => write!(
                f,
                "variational fixed-point solve diverged after {iters} iterations \
                 (residual {residual:e})"
            ),
            Se3Error::SingularJacobian => {
                write!(f, "adjoint residual Jacobian is numerically singular")
            }
        }
    }
}

impl std::error::Error for Se3Error {}

/// A body-frame twist: angular velocity (rad/s) then linear velocity
/// (m/s).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Twist {
    /// Body-frame angular velocity.
    pub omega: [f64; 3],
    /// Body-frame linear velocity.
    pub vel: [f64; 3],
}

impl Twist {
    fn is_finite(&self) -> bool {
        self.omega.iter().chain(self.vel.iter()).all(|v| v.is_finite())
    }
}

/// Deterministically fix the double-cover representative: the first
/// even component (in fs-ga's `EVEN_BLADES` order) whose magnitude
/// exceeds zero decides the sign; `M` and `−M` map to the SAME
/// canonical representative bit-for-bit.
pub fn canonicalize_motor(m: &Motor) -> Result<(Motor, bool), Se3Error> {
    for &blade in &EVEN_BLADES {
        let c = m.0.0[blade];
        if !c.is_finite() {
            return Err(Se3Error::NonFinite);
        }
        if c != 0.0 {
            return if c < 0.0 {
                Ok((Motor(m.0.scale(-1.0)), true))
            } else {
                Ok((*m, false))
            };
        }
    }
    Err(Se3Error::DegenerateMotor)
}

/// One exp-map step of the motor kinematics for a body-frame twist:
/// `M ← M ∘ exp(−(h/2)·B(ω, v))`, the SE(3) analogue of
/// [`crate::lie::quat_exp_step`]. The state stays on the group to
/// roundoff by construction; the canonical double-cover representative
/// is returned.
pub fn se3_exp_step(m: &Motor, twist: &Twist, h: f64) -> Result<Motor, Se3Error> {
    if !twist.is_finite() || !h.is_finite() {
        return Err(Se3Error::NonFinite);
    }
    let b = axis_bivector(twist.omega[0], twist.omega[1], twist.omega[2])
        .add(&ideal_bivector(twist.vel[0], twist.vel[1], twist.vel[2]))
        .scale(-0.5 * h);
    let step = exp_bivector(&b);
    let (canonical, _) = canonicalize_motor(&m.compose(&step))?;
    Ok(canonical)
}

/// Renormalization policy for long exp-lane runs.
#[derive(Debug, Clone, Copy)]
pub struct RenormPolicy {
    /// Renormalize when the unit defect `‖M M̃ − 1‖∞` exceeds this.
    pub defect_threshold: f64,
}

impl Default for RenormPolicy {
    fn default() -> Self {
        RenormPolicy {
            defect_threshold: 1e-12,
        }
    }
}

/// What one renormalization decision actually did (ledger fodder —
/// drift is reported, never silently absorbed).
#[derive(Debug, Clone, Copy)]
pub struct RenormReceipt {
    /// Unit defect measured BEFORE the decision.
    pub defect_before: f64,
    /// Whether the versor residue was divided out.
    pub renormalized: bool,
    /// Drift magnitude reported by `Motor::renormalize` (0 when not
    /// renormalized).
    pub drift: f64,
}

/// [`se3_exp_step`] plus drift-controlled renormalization with a
/// receipt.
pub fn se3_exp_step_renorm(
    m: &Motor,
    twist: &Twist,
    h: f64,
    policy: &RenormPolicy,
) -> Result<(Motor, RenormReceipt), Se3Error> {
    let mut next = se3_exp_step(m, twist, h)?;
    let defect_before = next.unit_defect();
    let (renormalized, drift) = if defect_before > policy.defect_threshold {
        let drift = next.renormalize();
        // Renormalization scales uniformly; rescaling cannot flip the
        // anchor sign, but keep the invariant explicit.
        let (canonical, _) = canonicalize_motor(&next)?;
        next = canonical;
        (true, drift)
    } else {
        (false, 0.0)
    };
    Ok((
        next,
        RenormReceipt {
            defect_before,
            renormalized,
            drift,
        },
    ))
}

fn validate_inertia(inertia: [f64; 3]) -> Result<(), Se3Error> {
    if inertia.iter().all(|i| i.is_finite() && *i > 0.0) {
        Ok(())
    } else {
        Err(Se3Error::InvalidParameter)
    }
}

/// Free rigid-body dynamics on SE(3): Euler's equations for the
/// angular velocity plus the body-frame transport of a spatially
/// constant linear velocity (`v̇_b = v_b × ω`), midpoint (RK2) in the
/// algebra and one exp-map motor update at the midpoint twist.
/// Returns the canonical motor and the updated twist.
pub fn se3_rigid_body_step(
    m: &Motor,
    twist: &Twist,
    inertia: [f64; 3],
    h: f64,
) -> Result<(Motor, Twist), Se3Error> {
    validate_inertia(inertia)?;
    if !twist.is_finite() || !(h.is_finite() && h != 0.0) {
        return Err(Se3Error::InvalidParameter);
    }
    let torque_free = |w: [f64; 3]| -> [f64; 3] {
        let l = [inertia[0] * w[0], inertia[1] * w[1], inertia[2] * w[2]];
        [
            l[1].mul_add(w[2], -(l[2] * w[1])) / inertia[0],
            l[2].mul_add(w[0], -(l[0] * w[2])) / inertia[1],
            l[0].mul_add(w[1], -(l[1] * w[0])) / inertia[2],
        ]
    };
    let vel_dot = |v: [f64; 3], w: [f64; 3]| -> [f64; 3] {
        // v̇_b = v_b × ω  (spatially constant free velocity).
        [
            v[1].mul_add(w[2], -(v[2] * w[1])),
            v[2].mul_add(w[0], -(v[0] * w[2])),
            v[0].mul_add(w[1], -(v[1] * w[0])),
        ]
    };
    let k1w = torque_free(twist.omega);
    let k1v = vel_dot(twist.vel, twist.omega);
    let mid = Twist {
        omega: [
            (0.5 * h).mul_add(k1w[0], twist.omega[0]),
            (0.5 * h).mul_add(k1w[1], twist.omega[1]),
            (0.5 * h).mul_add(k1w[2], twist.omega[2]),
        ],
        vel: [
            (0.5 * h).mul_add(k1v[0], twist.vel[0]),
            (0.5 * h).mul_add(k1v[1], twist.vel[1]),
            (0.5 * h).mul_add(k1v[2], twist.vel[2]),
        ],
    };
    let k2w = torque_free(mid.omega);
    let k2v = vel_dot(mid.vel, mid.omega);
    let next = Twist {
        omega: [
            h.mul_add(k2w[0], twist.omega[0]),
            h.mul_add(k2w[1], twist.omega[1]),
            h.mul_add(k2w[2], twist.omega[2]),
        ],
        vel: [
            h.mul_add(k2v[0], twist.vel[0]),
            h.mul_add(k2v[1], twist.vel[1]),
            h.mul_add(k2v[2], twist.vel[2]),
        ],
    };
    let motor = se3_exp_step(m, &mid, h)?;
    Ok((motor, next))
}

/// Solve controls for the variational fixed point.
#[derive(Debug, Clone, Copy)]
pub struct DepSolveParams {
    /// Convergence tolerance on the midpoint angular velocity update
    /// (absolute, rad/s).
    pub tol: f64,
    /// Iteration cap.
    pub max_iters: u32,
}

impl Default for DepSolveParams {
    fn default() -> Self {
        DepSolveParams {
            tol: 1e-14,
            max_iters: 64,
        }
    }
}

/// What one variational step's solve actually did.
#[derive(Debug, Clone, Copy)]
pub struct DepStepReceipt {
    /// Fixed-point iterations spent.
    pub iters: u32,
    /// Final update norm (the convergence metric).
    pub residual: f64,
    /// Whether `residual <= tol` was reached within the cap.
    pub converged: bool,
}

fn norm3(v: [f64; 3]) -> f64 {
    det::sqrt(v[0].mul_add(v[0], v[1].mul_add(v[1], v[2] * v[2])))
}

fn quat_conj(q: [f64; 4]) -> [f64; 4] {
    [q[0], -q[1], -q[2], -q[3]]
}

/// One discrete Euler–Poincaré step for the FREE rigid body in
/// body-momentum form: find the midpoint velocity `ω_m` with
/// `F = exp(h·ω̂_m)`, transport the body momentum by `F⁻¹`
/// (`Π' = F⁻¹ · Π`), and update the attitude by the SAME `F`
/// (`q' = q · F`). Spatial angular momentum `R I ω` is conserved
/// EXACTLY by construction — the theorem-class claim for energy is
/// decided separately by [`claim_for`].
pub fn dep_free_step(
    q: [f64; 4],
    omega: [f64; 3],
    inertia: [f64; 3],
    h: f64,
    params: &DepSolveParams,
) -> Result<([f64; 4], [f64; 3], DepStepReceipt), Se3Error> {
    validate_inertia(inertia)?;
    if !omega.iter().all(|v| v.is_finite()) || !(h.is_finite() && h != 0.0) {
        return Err(Se3Error::InvalidParameter);
    }
    let pi_k = [
        inertia[0] * omega[0],
        inertia[1] * omega[1],
        inertia[2] * omega[2],
    ];
    let mut w_mid = omega;
    let mut residual = f64::INFINITY;
    let mut iters = 0u32;
    let mut w_next = omega;
    while iters < params.max_iters {
        iters += 1;
        let f = quat_exp([h * w_mid[0], h * w_mid[1], h * w_mid[2]]);
        let pi_next = quat_rotate(quat_conj(f), pi_k);
        w_next = [
            pi_next[0] / inertia[0],
            pi_next[1] / inertia[1],
            pi_next[2] / inertia[2],
        ];
        let candidate = [
            0.5 * (omega[0] + w_next[0]),
            0.5 * (omega[1] + w_next[1]),
            0.5 * (omega[2] + w_next[2]),
        ];
        residual = norm3([
            candidate[0] - w_mid[0],
            candidate[1] - w_mid[1],
            candidate[2] - w_mid[2],
        ]);
        w_mid = candidate;
        if residual <= params.tol {
            break;
        }
    }
    let receipt = DepStepReceipt {
        iters,
        residual,
        converged: residual <= params.tol,
    };
    if !receipt.converged {
        return Err(Se3Error::SolverDiverged {
            iters,
            residual,
        });
    }
    let f = quat_exp([h * w_mid[0], h * w_mid[1], h * w_mid[2]]);
    let q_next = quat_mul(q, f);
    Ok((q_next, w_next, receipt))
}

/// Claim classes for variational runs. Composition never upgrades:
/// one violated assumption anywhere in a run demotes the whole run to
/// measured receipts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Se3ClaimClass {
    /// The discrete Euler–Poincaré symplectic/momentum theorem class:
    /// valid only for a declared smooth, conservative,
    /// regular-constraint fixture at fixed step with every solve
    /// converged.
    ConservativeVariationalTheorem,
    /// Honest fallback: the run carries MEASURED balance receipts and
    /// claims nothing structural.
    MeasuredOnly,
}

/// What the caller declares about the fixture; the integrator cannot
/// infer smoothness or conservativity, so the declaration is part of
/// the claim's provenance.
#[derive(Debug, Clone, Copy)]
pub struct Se3FixtureDeclaration {
    /// No dissipation, no external forcing.
    pub conservative: bool,
    /// Smooth/analytic dynamics (no impacts or switching).
    pub smooth: bool,
    /// Fixed step over the whole horizon.
    pub fixed_step: bool,
    /// Constraints (if any) are regular on the horizon.
    pub regular_constraints: bool,
}

/// Decide the claim class from the declaration and solver behavior.
#[must_use]
pub fn claim_for(decl: &Se3FixtureDeclaration, all_solves_converged: bool) -> Se3ClaimClass {
    if decl.conservative
        && decl.smooth
        && decl.fixed_step
        && decl.regular_constraints
        && all_solves_converged
    {
        Se3ClaimClass::ConservativeVariationalTheorem
    } else {
        Se3ClaimClass::MeasuredOnly
    }
}

/// Measured balance receipt for a variational run: what was actually
/// observed, independent of what theorem class the run earned.
#[derive(Debug, Clone)]
pub struct BalanceReceipt {
    /// The earned claim class.
    pub claim: Se3ClaimClass,
    /// Kinetic energy at the start.
    pub energy_start: f64,
    /// Kinetic energy at the end.
    pub energy_end: f64,
    /// Max |E(t) − E(0)| over the sampled series.
    pub energy_max_abs_drift: f64,
    /// Max componentwise |L_spatial(t) − L_spatial(0)|.
    pub momentum_max_abs_drift: f64,
    /// Steps taken.
    pub steps: usize,
    /// Whether every fixed-point solve converged.
    pub all_solves_converged: bool,
    /// Worst per-step iteration count.
    pub max_solver_iters: u32,
}

fn rotational_energy(omega: [f64; 3], inertia: [f64; 3]) -> f64 {
    0.5 * (inertia[0] * omega[0] * omega[0]
        + inertia[1] * omega[1] * omega[1]
        + inertia[2] * omega[2] * omega[2])
}

fn spatial_momentum(q: [f64; 4], omega: [f64; 3], inertia: [f64; 3]) -> [f64; 3] {
    quat_rotate(
        q,
        [
            inertia[0] * omega[0],
            inertia[1] * omega[1],
            inertia[2] * omega[2],
        ],
    )
}

/// Per-step multiplicative velocity damping for the honesty fixture:
/// `damping = 0` is the conservative case. Any nonzero damping demotes
/// the claim to [`Se3ClaimClass::MeasuredOnly`] regardless of how
/// small the measured drift looks.
pub fn run_dep_free(
    q0: [f64; 4],
    omega0: [f64; 3],
    inertia: [f64; 3],
    h: f64,
    steps: usize,
    damping: f64,
    params: &DepSolveParams,
) -> Result<([f64; 4], [f64; 3], BalanceReceipt), Se3Error> {
    if !damping.is_finite() || damping < 0.0 {
        return Err(Se3Error::InvalidParameter);
    }
    let mut q = q0;
    let mut w = omega0;
    let e0 = rotational_energy(w, inertia);
    let l0 = spatial_momentum(q, w, inertia);
    let mut energy_drift = 0.0f64;
    let mut momentum_drift = 0.0f64;
    let mut max_iters = 0u32;
    let mut all_converged = true;
    for _ in 0..steps {
        let (q1, w1, receipt) = dep_free_step(q, w, inertia, h, params)?;
        q = q1;
        w = w1;
        if damping > 0.0 {
            for wi in &mut w {
                *wi *= 1.0 - damping;
            }
        }
        all_converged &= receipt.converged;
        max_iters = max_iters.max(receipt.iters);
        let e = rotational_energy(w, inertia);
        energy_drift = energy_drift.max((e - e0).abs());
        let l = spatial_momentum(q, w, inertia);
        for a in 0..3 {
            momentum_drift = momentum_drift.max((l[a] - l0[a]).abs());
        }
    }
    let decl = Se3FixtureDeclaration {
        conservative: damping == 0.0,
        smooth: true,
        fixed_step: true,
        regular_constraints: true,
    };
    let receipt = BalanceReceipt {
        claim: claim_for(&decl, all_converged),
        energy_start: e0,
        energy_end: rotational_energy(w, inertia),
        energy_max_abs_drift: energy_drift,
        momentum_max_abs_drift: momentum_drift,
        steps,
        all_solves_converged: all_converged,
        max_solver_iters: max_iters,
    };
    Ok((q, w, receipt))
}

// ---------------------------------------------------------------------
// Discrete adjoint of the variational momentum map, derived from the
// ACTUAL fixed-point residual via the implicit-function theorem.
// ---------------------------------------------------------------------

/// The converged step's residual, as a function of (ω_mid, ω_k):
/// `g(ω_m, ω_k) = ω_m − ½·(ω_k + I⁻¹·(exp(h ω̂_m)⁻¹ · I ω_k))`.
fn dep_residual(w_mid: [f64; 3], w_k: [f64; 3], inertia: [f64; 3], h: f64) -> [f64; 3] {
    let f = quat_exp([h * w_mid[0], h * w_mid[1], h * w_mid[2]]);
    let pi_k = [
        inertia[0] * w_k[0],
        inertia[1] * w_k[1],
        inertia[2] * w_k[2],
    ];
    let pi_next = quat_rotate(quat_conj(f), pi_k);
    [
        w_mid[0] - 0.5 * (w_k[0] + pi_next[0] / inertia[0]),
        w_mid[1] - 0.5 * (w_k[1] + pi_next[1] / inertia[1]),
        w_mid[2] - 0.5 * (w_k[2] + pi_next[2] / inertia[2]),
    ]
}

/// Central-difference 3×3 Jacobian of the residual in its first or
/// second argument. This IS the actual discrete residual being
/// differenced (the bead's requirement); replacing the stencils with
/// analytic tangents is tracked follow-up work.
fn residual_jacobian(
    w_mid: [f64; 3],
    w_k: [f64; 3],
    inertia: [f64; 3],
    h: f64,
    wrt_mid: bool,
) -> [[f64; 3]; 3] {
    let mut j = [[0.0f64; 3]; 3];
    let base = if wrt_mid { w_mid } else { w_k };
    let scale = 1.0 + norm3(base);
    let eps = 1e-7 * scale;
    for col in 0..3 {
        let mut plus = base;
        let mut minus = base;
        plus[col] += eps;
        minus[col] -= eps;
        let (gp, gm) = if wrt_mid {
            (
                dep_residual(plus, w_k, inertia, h),
                dep_residual(minus, w_k, inertia, h),
            )
        } else {
            (
                dep_residual(w_mid, plus, inertia, h),
                dep_residual(w_mid, minus, inertia, h),
            )
        };
        for row in 0..3 {
            j[row][col] = (gp[row] - gm[row]) / (2.0 * eps);
        }
    }
    j
}

/// Solve `Aᵀ·x = b` for a 3×3 matrix by Gaussian elimination with
/// partial pivoting (the transpose solve the adjoint recursion needs).
fn solve3_transpose(a: &[[f64; 3]; 3], b: [f64; 3]) -> Result<[f64; 3], Se3Error> {
    // Form Aᵀ explicitly (3×3: clarity beats cleverness).
    let mut m = [[0.0f64; 4]; 3];
    for r in 0..3 {
        for c in 0..3 {
            m[r][c] = a[c][r];
        }
        m[r][3] = b[r];
    }
    for col in 0..3 {
        let mut pivot = col;
        for r in (col + 1)..3 {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        if m[pivot][col].abs() < 1e-300 {
            return Err(Se3Error::SingularJacobian);
        }
        m.swap(col, pivot);
        for r in (col + 1)..3 {
            let factor = m[r][col] / m[col][col];
            for c in col..4 {
                m[r][c] -= factor * m[col][c];
            }
        }
    }
    let mut x = [0.0f64; 3];
    for r in (0..3).rev() {
        let mut acc = m[r][3];
        for c in (r + 1)..3 {
            acc -= m[r][c] * x[c];
        }
        x[r] = acc / m[r][r];
    }
    Ok(x)
}

fn mat_t_vec(a: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0f64; 3];
    for c in 0..3 {
        for r in 0..3 {
            out[c] += a[r][c] * v[r];
        }
    }
    out
}

/// Discrete adjoint of the `steps`-step variational momentum
/// trajectory `ω_0 ↦ ω_N`: pulls a terminal cotangent `bar_ω_N` back
/// to `bar_ω_0` through the transposed implicit-function tangent of
/// each step's ACTUAL residual. The forward trajectory is recomputed
/// and stored (O(N) memory; revolve checkpointing is the follow-up,
/// matching the Verlet template).
pub fn dep_momentum_adjoint(
    omega0: [f64; 3],
    inertia: [f64; 3],
    h: f64,
    steps: usize,
    params: &DepSolveParams,
    bar_omega_n: [f64; 3],
) -> Result<[f64; 3], Se3Error> {
    // Forward sweep: record (ω_k, ω_mid_k) pairs.
    let mut trajectory = Vec::with_capacity(steps);
    let mut q = [1.0, 0.0, 0.0, 0.0];
    let mut w = omega0;
    for _ in 0..steps {
        let (q1, w1, _) = dep_free_step(q, w, inertia, h, params)?;
        let w_mid = [
            0.5 * (w[0] + w1[0]),
            0.5 * (w[1] + w1[1]),
            0.5 * (w[2] + w1[2]),
        ];
        trajectory.push((w, w_mid));
        q = q1;
        w = w1;
    }
    // Reverse sweep. With g(ω_m, ω_k) = 0 defining ω_m(ω_k) and
    // ω_{k+1} = 2·ω_m − ω_k:
    //   dω_{k+1}/dω_k = −2·(∂g/∂ω_m)⁻¹·(∂g/∂ω_k) − 1
    // so the transposed pull-back of a cotangent v is
    //   bar_ω_k = −2·(∂g/∂ω_k)ᵀ·(∂g/∂ω_m)⁻ᵀ·v − v.
    let mut bar = bar_omega_n;
    for (w_k, w_mid) in trajectory.iter().rev() {
        let dg_dmid = residual_jacobian(*w_mid, *w_k, inertia, h, true);
        let dg_dk = residual_jacobian(*w_mid, *w_k, inertia, h, false);
        let y = solve3_transpose(&dg_dmid, bar)?;
        let z = mat_t_vec(&dg_dk, y);
        bar = [
            (-2.0f64).mul_add(z[0], -bar[0]),
            (-2.0f64).mul_add(z[1], -bar[1]),
            (-2.0f64).mul_add(z[2], -bar[2]),
        ];
    }
    Ok(bar)
}

/// RATTLE-style constraint projection hook. The constrained
/// variational lanes live in fs-mbd; this trait is the seam they plug
/// into so fs-time never learns multibody types. Implementations
/// return the constraint-violation magnitude they removed (a receipt,
/// not a claim).
pub trait RattleProjection {
    /// Project the configuration back onto the constraint manifold.
    ///
    /// # Errors
    /// Implementation-defined refusals (irregular constraint,
    /// non-convergent projection).
    fn project_position(&self, motor: &mut Motor) -> Result<f64, Se3Error>;

    /// Project the twist onto the constraint's tangent space.
    ///
    /// # Errors
    /// Implementation-defined refusals.
    fn project_velocity(&self, motor: &Motor, twist: &mut Twist) -> Result<f64, Se3Error>;
}

/// The trivial (unconstrained) projection: removes nothing, refuses
/// nothing.
#[derive(Debug, Clone, Copy, Default)]
pub struct Unconstrained;

impl RattleProjection for Unconstrained {
    fn project_position(&self, _motor: &mut Motor) -> Result<f64, Se3Error> {
        Ok(0.0)
    }

    fn project_velocity(&self, _motor: &Motor, _twist: &mut Twist) -> Result<f64, Se3Error> {
        Ok(0.0)
    }
}
