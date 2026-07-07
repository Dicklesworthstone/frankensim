//! Lie-group integrators on SO(3) (unit quaternions): exponential-map
//! updates keep the state ON the group to roundoff by construction — no
//! renormalization hacks, no drift accumulation (measured over 10⁵
//! steps in the battery). SE(3) rods and fs-ga motor states build on
//! this once fs-ga lands (recorded).

use fs_math::det;

/// Quaternion product (w, x, y, z convention).
#[must_use]
pub fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0].mul_add(b[0], -(a[1] * b[1]) - a[2] * b[2] - a[3] * b[3]),
        a[0].mul_add(b[1], a[1].mul_add(b[0], a[2].mul_add(b[3], -(a[3] * b[2])))),
        a[0].mul_add(b[2], a[2].mul_add(b[0], a[3].mul_add(b[1], -(a[1] * b[3])))),
        a[0].mul_add(b[3], a[3].mul_add(b[0], a[1].mul_add(b[2], -(a[2] * b[1])))),
    ]
}

/// exp of a pure quaternion ½·h·ω (the SO(3) exponential in quaternion
/// form): cos(θ) + sin(θ)·(axis), θ = ‖h·ω‖/2. Small-angle branch uses
/// the series for sin(θ)/θ (no cancellation).
#[must_use]
pub fn quat_exp(hw: [f64; 3]) -> [f64; 4] {
    let half = [0.5 * hw[0], 0.5 * hw[1], 0.5 * hw[2]];
    let theta2 = half[0].mul_add(half[0], half[1].mul_add(half[1], half[2] * half[2]));
    let theta = det::sqrt(theta2);
    let (c, sinc) = if theta < 1e-6 {
        // cos θ ≈ 1 − θ²/2 + θ⁴/24; sin θ/θ ≈ 1 − θ²/6 + θ⁴/120.
        (
            theta2.mul_add(theta2 / 24.0, 1.0 - 0.5 * theta2),
            theta2.mul_add(theta2 / 120.0, 1.0 - theta2 / 6.0),
        )
    } else {
        (det::cos(theta), det::sin(theta) / theta)
    };
    [c, sinc * half[0], sinc * half[1], sinc * half[2]]
}

/// One exp-map step of the attitude kinematics q̇ = ½·q·ω̂ (body-frame
/// angular velocity ω): q ← q·exp(½hω). Norm preserved by construction.
#[must_use]
pub fn quat_exp_step(q: [f64; 4], omega_body: [f64; 3], h: f64) -> [f64; 4] {
    let step = quat_exp([h * omega_body[0], h * omega_body[1], h * omega_body[2]]);
    quat_mul(q, step)
}

/// One commutator-free CG2 (midpoint Lie-group) step of the FREE rigid
/// body with principal inertia `inertia`: Euler's equations for ω in the
/// body frame + exp-map attitude update at the midpoint angular
/// velocity. Returns (q', ω').
#[must_use]
pub fn rigid_body_step(
    q: [f64; 4],
    omega: [f64; 3],
    inertia: [f64; 3],
    h: f64,
) -> ([f64; 4], [f64; 3]) {
    let torque_free = |w: [f64; 3]| -> [f64; 3] {
        // ω̇ = I⁻¹·(Iω × ω)
        let l = [inertia[0] * w[0], inertia[1] * w[1], inertia[2] * w[2]];
        [
            l[1].mul_add(w[2], -(l[2] * w[1])) / inertia[0],
            l[2].mul_add(w[0], -(l[0] * w[2])) / inertia[1],
            l[0].mul_add(w[1], -(l[1] * w[0])) / inertia[2],
        ]
    };
    // Midpoint (RK2) for ω.
    let k1 = torque_free(omega);
    let w_mid = [
        (0.5 * h).mul_add(k1[0], omega[0]),
        (0.5 * h).mul_add(k1[1], omega[1]),
        (0.5 * h).mul_add(k1[2], omega[2]),
    ];
    let k2 = torque_free(w_mid);
    let omega_new = [
        h.mul_add(k2[0], omega[0]),
        h.mul_add(k2[1], omega[1]),
        h.mul_add(k2[2], omega[2]),
    ];
    // CG2 attitude update at the midpoint velocity (group-preserving).
    let q_new = quat_exp_step(q, w_mid, h);
    (q_new, omega_new)
}

/// Rotate a body-frame vector to the spatial frame by quaternion q.
#[must_use]
pub fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let qv = [0.0, v[0], v[1], v[2]];
    let qc = [q[0], -q[1], -q[2], -q[3]];
    let r = quat_mul(quat_mul(q, qv), qc);
    [r[1], r[2], r[3]]
}
