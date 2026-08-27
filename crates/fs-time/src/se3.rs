//! SE(3) Lie-group and variational integrator lanes on PGA motors
//! (bead `frankensim-ext-time-se3-lanes-3ol0`).
//!
//! Three obligations from the bead, each with an honest claim scope:
//!
//! 1. **Exp-map lane**: canonical `fs-ga` `Se3` states updated through
//!    explicit body/right or space/left group operations, with drift-controlled
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

use fs_ga::{GaError, Se3, So3, So3Tangent, Twist, Vec3};
use fs_math::det;

/// Typed refusals for the SE(3) lanes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Se3Error {
    /// Canonical group/tangent validation refused an input or result.
    Group(GaError),
    /// A time-integration control or physical parameter is outside its
    /// documented domain.
    InvalidParameter {
        /// Parameter family that was invalid.
        context: &'static str,
    },
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
            Se3Error::Group(error) => write!(f, "SE(3) group refusal: {error}"),
            Se3Error::InvalidParameter { context } => write!(f, "invalid {context}"),
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

impl From<GaError> for Se3Error {
    fn from(value: GaError) -> Self {
        Self::Group(value)
    }
}

/// One exp-map step for a body-frame twist:
/// `T_next = T * Exp(h * twist_body)`.
pub fn se3_exp_step(pose: Se3, twist: Twist, h: f64) -> Result<Se3, Se3Error> {
    if !h.is_finite() {
        return Err(Se3Error::InvalidParameter {
            context: "SE(3) step size",
        });
    }
    Ok(pose.body_plus(twist.scale(h))?)
}

/// One exp-map step for a space-frame twist:
/// `T_next = Exp(h * twist_space) * T`.
pub fn se3_space_exp_step(pose: Se3, twist: Twist, h: f64) -> Result<Se3, Se3Error> {
    if !h.is_finite() {
        return Err(Se3Error::InvalidParameter {
            context: "SE(3) step size",
        });
    }
    Ok(pose.space_plus(twist.scale(h))?)
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
    pose: Se3,
    twist: Twist,
    h: f64,
    policy: &RenormPolicy,
) -> Result<(Se3, RenormReceipt), Se3Error> {
    if !policy.defect_threshold.is_finite() || policy.defect_threshold < 0.0 {
        return Err(Se3Error::InvalidParameter {
            context: "SE(3) renormalization defect threshold",
        });
    }
    let mut next = se3_exp_step(pose, twist, h)?;
    let defect_before = next.as_motor().unit_defect();
    let (renormalized, drift) = if defect_before > policy.defect_threshold {
        let mut corrected = *next.as_motor();
        let drift = corrected.renormalize();
        next = Se3::try_from_motor(corrected)?;
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

fn validate_inertia(inertia: Vec3) -> Result<(), Se3Error> {
    if [inertia.x, inertia.y, inertia.z]
        .iter()
        .all(|value| value.is_finite() && *value > 0.0)
    {
        Ok(())
    } else {
        Err(Se3Error::InvalidParameter {
            context: "rigid-body principal inertia",
        })
    }
}

fn validate_dep_params(params: &DepSolveParams) -> Result<(), Se3Error> {
    if params.tol.is_finite() && params.tol > 0.0 && params.max_iters > 0 {
        Ok(())
    } else {
        Err(Se3Error::InvalidParameter {
            context: "variational solve controls",
        })
    }
}

fn validate_angular_step(omega: Vec3, h: f64) -> Result<(), Se3Error> {
    if [omega.x, omega.y, omega.z]
        .iter()
        .all(|value| value.is_finite())
        && h.is_finite()
        && h != 0.0
    {
        Ok(())
    } else {
        Err(Se3Error::InvalidParameter {
            context: "angular velocity or nonzero step size",
        })
    }
}

/// Free rigid-body dynamics on SE(3): Euler's equations for the
/// angular velocity plus the body-frame transport of a spatially
/// constant linear velocity (`v̇_b = v_b × ω`), midpoint (RK2) in the
/// algebra and one exp-map motor update at the midpoint twist.
/// Returns the canonical motor and the updated twist.
pub fn se3_rigid_body_step(
    pose: Se3,
    twist: Twist,
    inertia: Vec3,
    h: f64,
) -> Result<(Se3, Twist), Se3Error> {
    validate_inertia(inertia)?;
    if !(twist.to_array().iter().all(|value| value.is_finite()) && h.is_finite() && h != 0.0) {
        return Err(Se3Error::InvalidParameter {
            context: "rigid-body twist or nonzero step size",
        });
    }
    let torque_free = |omega: Vec3| -> Vec3 {
        let momentum = Vec3::new(
            inertia.x * omega.x,
            inertia.y * omega.y,
            inertia.z * omega.z,
        );
        Vec3::new(
            momentum.y.mul_add(omega.z, -(momentum.z * omega.y)) / inertia.x,
            momentum.z.mul_add(omega.x, -(momentum.x * omega.z)) / inertia.y,
            momentum.x.mul_add(omega.y, -(momentum.y * omega.x)) / inertia.z,
        )
    };
    let velocity_derivative = |linear: Vec3, angular: Vec3| -> Vec3 {
        // v̇_b = v_b × ω  (spatially constant free velocity).
        linear.cross(angular)
    };
    let first_angular_slope = torque_free(twist.angular);
    let first_linear_slope = velocity_derivative(twist.linear, twist.angular);
    let midpoint = Twist::new(
        Vec3::new(
            (0.5 * h).mul_add(first_angular_slope.x, twist.angular.x),
            (0.5 * h).mul_add(first_angular_slope.y, twist.angular.y),
            (0.5 * h).mul_add(first_angular_slope.z, twist.angular.z),
        ),
        Vec3::new(
            (0.5 * h).mul_add(first_linear_slope.x, twist.linear.x),
            (0.5 * h).mul_add(first_linear_slope.y, twist.linear.y),
            (0.5 * h).mul_add(first_linear_slope.z, twist.linear.z),
        ),
    );
    let second_angular_slope = torque_free(midpoint.angular);
    let second_linear_slope = velocity_derivative(midpoint.linear, midpoint.angular);
    let next = Twist::new(
        Vec3::new(
            h.mul_add(second_angular_slope.x, twist.angular.x),
            h.mul_add(second_angular_slope.y, twist.angular.y),
            h.mul_add(second_angular_slope.z, twist.angular.z),
        ),
        Vec3::new(
            h.mul_add(second_linear_slope.x, twist.linear.x),
            h.mul_add(second_linear_slope.y, twist.linear.y),
            h.mul_add(second_linear_slope.z, twist.linear.z),
        ),
    );
    let next_pose = se3_exp_step(pose, midpoint, h)?;
    Ok((next_pose, next))
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

fn norm3(value: Vec3) -> f64 {
    det::sqrt(
        value
            .x
            .mul_add(value.x, value.y.mul_add(value.y, value.z * value.z)),
    )
}

/// One discrete Euler–Poincaré step for the FREE rigid body in
/// body-momentum form: find the midpoint velocity `ω_m` with
/// `F = exp(h·ω̂_m)`, transport the body momentum by `F⁻¹`
/// (`Π' = F⁻¹ · Π`), and update the attitude by the SAME `F`
/// (`q' = q · F`). Spatial angular momentum `R I ω` is conserved
/// EXACTLY by construction — the theorem-class claim for energy is
/// decided separately by [`claim_for`].
pub fn dep_free_step(
    rotation: So3,
    omega: Vec3,
    inertia: Vec3,
    h: f64,
    params: &DepSolveParams,
) -> Result<(So3, Vec3, DepStepReceipt), Se3Error> {
    validate_inertia(inertia)?;
    validate_angular_step(omega, h)?;
    validate_dep_params(params)?;
    let momentum = Vec3::new(
        inertia.x * omega.x,
        inertia.y * omega.y,
        inertia.z * omega.z,
    );
    let mut w_mid = omega;
    let mut residual = f64::INFINITY;
    let mut iters = 0u32;
    while iters < params.max_iters {
        iters += 1;
        let increment = So3::exp(So3Tangent::new(w_mid.scale(h)))?;
        let next_momentum = increment.inverse().rotate(momentum)?;
        let next_velocity = Vec3::new(
            next_momentum.x / inertia.x,
            next_momentum.y / inertia.y,
            next_momentum.z / inertia.z,
        );
        let candidate = Vec3::new(
            0.5 * (omega.x + next_velocity.x),
            0.5 * (omega.y + next_velocity.y),
            0.5 * (omega.z + next_velocity.z),
        );
        residual = norm3(candidate - w_mid);
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
        return Err(Se3Error::SolverDiverged { iters, residual });
    }
    // Re-evaluate the momentum transport at the accepted midpoint. Without
    // this final evaluation, `w_next` belongs to the previous fixed-point
    // iterate while the attitude uses the accepted one, invalidating the
    // same-increment momentum-conservation construction by the solve tolerance.
    let increment = So3::exp(So3Tangent::new(w_mid.scale(h)))?;
    let next_momentum = increment.inverse().rotate(momentum)?;
    let w_next = Vec3::new(
        next_momentum.x / inertia.x,
        next_momentum.y / inertia.y,
        next_momentum.z / inertia.z,
    );
    let next_rotation = rotation.compose(increment)?;
    Ok((next_rotation, w_next, receipt))
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
#[allow(clippy::struct_excessive_bools)] // Four independent theorem hypotheses remain explicit.
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

fn rotational_energy(omega: Vec3, inertia: Vec3) -> f64 {
    0.5 * (inertia.x * omega.x * omega.x
        + inertia.y * omega.y * omega.y
        + inertia.z * omega.z * omega.z)
}

fn spatial_momentum(rotation: So3, omega: Vec3, inertia: Vec3) -> Result<Vec3, Se3Error> {
    Ok(rotation.rotate(Vec3::new(
        inertia.x * omega.x,
        inertia.y * omega.y,
        inertia.z * omega.z,
    ))?)
}

/// Per-step multiplicative velocity damping for the honesty fixture:
/// `damping = 0` is the conservative case. Any nonzero damping demotes
/// the claim to [`Se3ClaimClass::MeasuredOnly`] regardless of how
/// small the measured drift looks.
pub fn run_dep_free(
    rotation0: So3,
    omega0: Vec3,
    inertia: Vec3,
    h: f64,
    steps: usize,
    damping: f64,
    params: &DepSolveParams,
) -> Result<(So3, Vec3, BalanceReceipt), Se3Error> {
    validate_inertia(inertia)?;
    validate_angular_step(omega0, h)?;
    validate_dep_params(params)?;
    if !damping.is_finite() || !(0.0..=1.0).contains(&damping) {
        return Err(Se3Error::InvalidParameter {
            context: "per-step damping fraction",
        });
    }
    let mut rotation = rotation0;
    let mut w = omega0;
    let e0 = rotational_energy(w, inertia);
    let l0 = spatial_momentum(rotation, w, inertia)?;
    let mut energy_drift = 0.0f64;
    let mut momentum_drift = 0.0f64;
    let mut max_iters = 0u32;
    let mut all_converged = true;
    for _ in 0..steps {
        let (next_rotation, w1, receipt) = dep_free_step(rotation, w, inertia, h, params)?;
        rotation = next_rotation;
        w = w1;
        if damping > 0.0 {
            w = w.scale(1.0 - damping);
        }
        all_converged &= receipt.converged;
        max_iters = max_iters.max(receipt.iters);
        let e = rotational_energy(w, inertia);
        energy_drift = energy_drift.max((e - e0).abs());
        let momentum_delta = spatial_momentum(rotation, w, inertia)? - l0;
        momentum_drift = momentum_drift
            .max(momentum_delta.x.abs())
            .max(momentum_delta.y.abs())
            .max(momentum_delta.z.abs());
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
    Ok((rotation, w, receipt))
}

// ---------------------------------------------------------------------
// Discrete adjoint of the variational momentum map, derived from the
// ACTUAL fixed-point residual via the implicit-function theorem.
// ---------------------------------------------------------------------

/// The converged step's residual, as a function of (ω_mid, ω_k):
/// `g(ω_m, ω_k) = ω_m − ½·(ω_k + I⁻¹·(exp(h ω̂_m)⁻¹ · I ω_k))`.
fn dep_residual(w_mid: Vec3, w_k: Vec3, inertia: Vec3, h: f64) -> Result<Vec3, Se3Error> {
    let increment = So3::exp(So3Tangent::new(w_mid.scale(h)))?;
    let momentum = Vec3::new(inertia.x * w_k.x, inertia.y * w_k.y, inertia.z * w_k.z);
    let next_momentum = increment.inverse().rotate(momentum)?;
    Ok(Vec3::new(
        w_mid.x - 0.5 * (w_k.x + next_momentum.x / inertia.x),
        w_mid.y - 0.5 * (w_k.y + next_momentum.y / inertia.y),
        w_mid.z - 0.5 * (w_k.z + next_momentum.z / inertia.z),
    ))
}

/// Central-difference 3×3 Jacobian of the residual in its first or
/// second argument. This IS the actual discrete residual being
/// differenced (the bead's requirement); replacing the stencils with
/// analytic tangents is tracked follow-up work.
fn residual_jacobian(
    w_mid: Vec3,
    w_k: Vec3,
    inertia: Vec3,
    h: f64,
    wrt_mid: bool,
) -> Result<[[f64; 3]; 3], Se3Error> {
    let mut columns = [[0.0f64; 3]; 3];
    let base = if wrt_mid { w_mid } else { w_k };
    let scale = 1.0 + norm3(base);
    let eps = 1e-7 * scale;
    let bases = [
        Vec3::new(eps, 0.0, 0.0),
        Vec3::new(0.0, eps, 0.0),
        Vec3::new(0.0, 0.0, eps),
    ];
    for (column, basis) in columns.iter_mut().zip(bases) {
        let plus = base + basis;
        let minus = base - basis;
        let (gp, gm) = if wrt_mid {
            (
                dep_residual(plus, w_k, inertia, h)?,
                dep_residual(minus, w_k, inertia, h)?,
            )
        } else {
            (
                dep_residual(w_mid, plus, inertia, h)?,
                dep_residual(w_mid, minus, inertia, h)?,
            )
        };
        let gp = [gp.x, gp.y, gp.z];
        let gm = [gm.x, gm.y, gm.z];
        for row in 0..3 {
            column[row] = (gp[row] - gm[row]) / (2.0 * eps);
        }
    }
    Ok([
        [columns[0][0], columns[1][0], columns[2][0]],
        [columns[0][1], columns[1][1], columns[2][1]],
        [columns[0][2], columns[1][2], columns[2][2]],
    ])
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
        let pivot_row = m[col];
        for r in (col + 1)..3 {
            let factor = m[r][col] / m[col][col];
            for (entry, pivot_entry) in m[r][col..].iter_mut().zip(&pivot_row[col..]) {
                *entry -= factor * *pivot_entry;
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
    omega0: Vec3,
    inertia: Vec3,
    h: f64,
    steps: usize,
    params: &DepSolveParams,
    bar_omega_n: Vec3,
) -> Result<Vec3, Se3Error> {
    validate_inertia(inertia)?;
    validate_angular_step(omega0, h)?;
    validate_dep_params(params)?;
    if ![bar_omega_n.x, bar_omega_n.y, bar_omega_n.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(Se3Error::InvalidParameter {
            context: "terminal angular-velocity cotangent",
        });
    }
    // Forward sweep: record (ω_k, ω_mid_k) pairs.
    let mut trajectory = Vec::with_capacity(steps);
    let mut rotation = So3::identity();
    let mut w = omega0;
    for _ in 0..steps {
        let (next_rotation, w1, _) = dep_free_step(rotation, w, inertia, h, params)?;
        let w_mid = Vec3::new(0.5 * (w.x + w1.x), 0.5 * (w.y + w1.y), 0.5 * (w.z + w1.z));
        trajectory.push((w, w_mid));
        rotation = next_rotation;
        w = w1;
    }
    // Reverse sweep. With g(ω_m, ω_k) = 0 defining ω_m(ω_k) and
    // ω_{k+1} = 2·ω_m − ω_k:
    //   dω_{k+1}/dω_k = −2·(∂g/∂ω_m)⁻¹·(∂g/∂ω_k) − 1
    // so the transposed pull-back of a cotangent v is
    //   bar_ω_k = −2·(∂g/∂ω_k)ᵀ·(∂g/∂ω_m)⁻ᵀ·v − v.
    let mut bar = [bar_omega_n.x, bar_omega_n.y, bar_omega_n.z];
    for (w_k, w_mid) in trajectory.iter().rev() {
        let dg_dmid = residual_jacobian(*w_mid, *w_k, inertia, h, true)?;
        let dg_dk = residual_jacobian(*w_mid, *w_k, inertia, h, false)?;
        let y = solve3_transpose(&dg_dmid, bar)?;
        let z = mat_t_vec(&dg_dk, y);
        bar = [
            (-2.0f64).mul_add(z[0], -bar[0]),
            (-2.0f64).mul_add(z[1], -bar[1]),
            (-2.0f64).mul_add(z[2], -bar[2]),
        ];
    }
    Ok(Vec3::new(bar[0], bar[1], bar[2]))
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
    fn project_position(&self, pose: &mut Se3) -> Result<f64, Se3Error>;

    /// Project the twist onto the constraint's tangent space.
    ///
    /// # Errors
    /// Implementation-defined refusals.
    fn project_velocity(&self, pose: &Se3, twist: &mut Twist) -> Result<f64, Se3Error>;
}

/// The trivial (unconstrained) projection: removes nothing, refuses
/// nothing.
#[derive(Debug, Clone, Copy, Default)]
pub struct Unconstrained;

impl RattleProjection for Unconstrained {
    fn project_position(&self, _pose: &mut Se3) -> Result<f64, Se3Error> {
        Ok(0.0)
    }

    fn project_velocity(&self, _pose: &Se3, _twist: &mut Twist) -> Result<f64, Se3Error> {
        Ok(0.0)
    }
}
