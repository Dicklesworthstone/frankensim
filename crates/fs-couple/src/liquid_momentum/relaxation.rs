//! Exact uniform-rate drag with a finite six-coordinate kinetic metric.
//!
//! With x the body's mass-whitened velocity and K_i its point-velocity map,
//! u_i' = gamma (K_i x-u_i), x' = gamma sum m_i K_i^T (u_i-K_i x).
//! C = sum m_i K_i^T K_i is positive semidefinite and only 6 by 6.
//! d = sum m_i K_i^T (u_i-K_i x) obeys d'=-gamma (I+C)d.
//! The remaining fluid slip modes decay at gamma. Resolving both sets gives
//! the exact frozen-position source map without a dense 3N by 3N operator.

use super::{LiquidMomentumError, finite_vec, invalid, poll};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::Vec3;

#[derive(Clone, Copy, Debug)]
pub(super) struct Cell {
    pub mass: f64,
    pub arm: Vec3,
    pub velocity: Vec3,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Body {
    pub velocity: Vec3,
    pub omega: Vec3,
    pub inverse_sqrt_mass: f64,
    pub inverse_sqrt_inertia: [f64; 3],
    /// Original inertia principal axes resolved in world coordinates.
    pub principal_axes: [Vec3; 3],
}

impl Body {
    pub(super) fn point_velocity(self, arm: Vec3) -> Vec3 {
        self.velocity.add(self.omega.cross(arm))
    }

    pub(super) fn columns(self, arm: Vec3) -> [Vec3; 6] {
        let m = self.inverse_sqrt_mass;
        [
            Vec3::new(m, 0.0, 0.0),
            Vec3::new(0.0, m, 0.0),
            Vec3::new(0.0, 0.0, m),
            self.principal_axes[0]
                .cross(arm)
                .scale(self.inverse_sqrt_inertia[0]),
            self.principal_axes[1]
                .cross(arm)
                .scale(self.inverse_sqrt_inertia[1]),
            self.principal_axes[2]
                .cross(arm)
                .scale(self.inverse_sqrt_inertia[2]),
        ]
    }

    pub(super) fn physical_impulse(self, increment: [f64; 6]) -> (Vec3, Vec3) {
        let linear =
            Vec3::new(increment[0], increment[1], increment[2]).scale(1.0 / self.inverse_sqrt_mass);
        let angular = (0..3).fold(Vec3::ZERO, |sum, k| {
            sum.add(self.principal_axes[k].scale(increment[k + 3] / self.inverse_sqrt_inertia[k]))
        });
        (linear, angular)
    }
}

pub(super) struct Plan {
    /// Total impulses ON the original liquid cells, in caller order.
    pub impulses: Vec<Vec3>,
    /// Independently integrated native-body recoil, about its fixed COM.
    pub body_linear: Vec3,
    pub body_angular: Vec3,
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn dot(x: Vec3, y: Vec3) -> Dd {
    dd(x.x) * dd(y.x) + dd(x.y) * dd(y.y) + dd(x.z) * dd(y.z)
}
fn vector(v: [Dd; 3]) -> Vec3 {
    Vec3::new(v[0].to_f64(), v[1].to_f64(), v[2].to_f64())
}

/// Integral weights for the body mode and the forced fluid mode.
/// Small z*(1+c) uses a divided-power series, so a fluid O(z²) response is
/// not obtained by subtracting two O(z) exponentials. c=0 has its true limit.
fn weights(z: f64, c: f64) -> Result<(f64, f64), LiquidMomentumError> {
    let rate = 1.0 + c;
    let t = z * rate;
    let body = -(-t).exp_m1() / rate;
    let fluid = if t <= 0.125 {
        let mut sum = 0.5;
        let mut divided = 1.0;
        let mut power = 1.0;
        let mut factorial = 2.0;
        for k in 1..=18 {
            power *= z;
            divided = t.mul_add(divided, power);
            factorial *= (k + 2) as f64;
            sum += if k % 2 == 0 {
                divided / factorial
            } else {
                -divided / factorial
            };
        }
        z * z * sum
    } else {
        let decay = (-z).exp();
        if decay == 0.0 {
            1.0 / rate
        } else {
            let argument = z * c;
            let phi = if argument == 0.0 {
                1.0
            } else {
                -(-argument).exp_m1() / argument
            };
            (-(-z).exp_m1() - decay * (z * phi)) / rate
        }
    };
    if !body.is_finite() || !fluid.is_finite() || body < 0.0 || fluid < 0.0 {
        return Err(invalid("exact liquid drag weights are unrepresentable"));
    }
    Ok((body, fluid))
}

pub(super) fn solve(
    cells: &[Cell],
    body: Body,
    exposure: f64,
    gate: &CancelGate,
) -> Result<Plan, LiquidMomentumError> {
    poll(gate)?;
    let mut gram = [Dd::ZERO; 36];
    let mut defect = [Dd::ZERO; 6];
    for (index, cell) in cells.iter().enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let columns = body.columns(cell.arm);
        let slip = cell.velocity.sub(body.point_velocity(cell.arm));
        finite_vec(slip)?;
        for j in 0..6 {
            defect[j] = defect[j] + dd(cell.mass) * dot(columns[j], slip);
            for k in j..6 {
                gram[j * 6 + k] = gram[j * 6 + k] + dd(cell.mass) * dot(columns[j], columns[k]);
            }
        }
    }
    let mut matrix = [0.0; 36];
    for j in 0..6 {
        for k in j..6 {
            matrix[j * 6 + k] = gram[j * 6 + k].to_f64();
            matrix[k * 6 + j] = matrix[j * 6 + k];
        }
    }
    let defect = defect.map(Dd::to_f64);
    if matrix.iter().chain(&defect).any(|v| !v.is_finite()) {
        return Err(invalid("liquid/rigid kinetic response is unrepresentable"));
    }
    let scale = matrix.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    let (values, vectors) = if scale == 0.0 {
        let mut identity = vec![0.0; 36];
        for i in 0..6 {
            identity[i * 6 + i] = 1.0;
        }
        (vec![0.0; 6], identity)
    } else {
        let normalized = matrix.map(|v| v / scale);
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&normalized, 6);
        let mut actual = Vec::with_capacity(6);
        for value in values {
            if !value.is_finite() || value < -256.0 * f64::EPSILON {
                return Err(invalid(
                    "liquid/rigid positive kinetic metric failed eigensolve",
                ));
            }
            // Only negative roundoff in a known PSD Gram form is removed;
            // every positive resolved eigenvalue, including small ones, stays.
            let value = value.max(0.0) * scale;
            if !value.is_finite() {
                return Err(invalid("liquid/rigid kinetic rate is unrepresentable"));
            }
            actual.push(value);
        }
        (actual, vectors)
    };
    poll(gate)?;
    let mut body_change = [Dd::ZERO; 6];
    let mut fluid_response = [Dd::ZERO; 6];
    for mode in 0..6 {
        let amplitude = (0..6).fold(Dd::ZERO, |sum, j| {
            sum + dd(vectors[j * 6 + mode]) * dd(defect[j])
        });
        let (body_weight, fluid_weight) = weights(exposure, values[mode])?;
        for j in 0..6 {
            body_change[j] =
                body_change[j] + dd(vectors[j * 6 + mode]) * dd(body_weight) * amplitude;
            fluid_response[j] =
                fluid_response[j] + dd(vectors[j * 6 + mode]) * dd(fluid_weight) * amplitude;
        }
    }
    let body_change = body_change.map(Dd::to_f64);
    let response = fluid_response.map(Dd::to_f64);
    let (body_linear, body_angular) = body.physical_impulse(body_change);
    finite_vec(body_linear)?;
    finite_vec(body_angular)?;
    let decay_change = -(-exposure).exp_m1();
    let mut impulses = Vec::new();
    impulses
        .try_reserve_exact(cells.len())
        .map_err(|_| invalid("compact liquid drag impulse allocation refused"))?;
    for (index, cell) in cells.iter().enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let columns = body.columns(cell.arm);
        let slip = cell.velocity.sub(body.point_velocity(cell.arm));
        let mut change = [
            dd(-decay_change) * dd(slip.x),
            dd(-decay_change) * dd(slip.y),
            dd(-decay_change) * dd(slip.z),
        ];
        for j in 0..6 {
            change[0] = change[0] + dd(columns[j].x) * dd(response[j]);
            change[1] = change[1] + dd(columns[j].y) * dd(response[j]);
            change[2] = change[2] + dd(columns[j].z) * dd(response[j]);
        }
        let impulse = vector(change.map(|v| dd(cell.mass) * v));
        finite_vec(impulse)?;
        impulses.push(impulse);
    }
    poll(gate)?;
    Ok(Plan {
        impulses,
        body_linear,
        body_angular,
    })
}

#[cfg(test)]
mod tests;
