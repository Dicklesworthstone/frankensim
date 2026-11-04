//! Vector-speed inertial porous drag on the native implicit resistance step.

use super::*;

/// Prescribed Darcy–Forchheimer law `F = b - D v - quadratic |v| v`,
/// where `v = u - skeleton_velocity`.
///
/// The linear tensor may be anisotropic. The quadratic term uses the full
/// Euclidean relative speed and is rotationally invariant; it is not a
/// separable sum of principal-component quadratic drags. `quadratic` is a
/// force-density/velocity-squared coefficient in lattice units (`dx=dt=1`),
/// already including any density factor from its physical constitutive law.
/// No extra density or factor one half is applied to the force.
///
/// Sample all coefficients, independent body force and prescribed skeleton
/// motion at the step midpoint and hold them for the tick. Zero quadratic
/// resistance takes the exact original [`LinearResistance3`] arithmetic.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DarcyForchheimer3 {
    /// Linear tensor, material axes, skeleton velocity and independent force.
    pub linear: LinearResistance3,
    /// Finite nonnegative vector-speed resistance coefficient.
    pub quadratic: f64,
}

impl From<LinearResistance3> for DarcyForchheimer3 {
    fn from(linear: LinearResistance3) -> Self {
        Self {
            linear,
            quadratic: 0.0,
        }
    }
}

impl DarcyForchheimer3 {
    pub(super) fn admit(self) -> Result<(), BoundaryStepError3> {
        self.linear.admit()?;
        if !self.quadratic.is_finite() || self.quadratic < 0.0 {
            return Err(BoundaryStepError3::EndpointForce {
                reason: "Forchheimer vector-speed resistance must be finite and nonnegative",
            });
        }
        Ok(())
    }

    pub(super) fn drag(self, velocity: [f64; 3]) -> [f64; 3] {
        let linear = self.linear.drag(velocity);
        if self.quadratic == 0.0 {
            return linear;
        }
        let relative = core::array::from_fn(|a| self.linear.skeleton_velocity[a] - velocity[a]);
        let coefficient = self.quadratic * norm(relative);
        core::array::from_fn(|a| linear[a] + coefficient * relative[a])
    }

    pub(super) fn force(self, velocity: [f64; 3]) -> [f64; 3] {
        if self.quadratic == 0.0 {
            return self.linear.force(velocity);
        }
        let drag = self.drag(velocity);
        core::array::from_fn(|a| self.linear.body_force[a] + drag[a])
    }

    pub(super) fn endpoint_force(
        self,
        rho: f64,
        momentum: [f64; 3],
    ) -> Option<([f64; 3], [f64; 3])> {
        if self.quadratic == 0.0 {
            return self.linear.endpoint_force(rho, momentum);
        }
        let law = self.linear;
        let rhs = core::array::from_fn(|a| {
            momentum[a] + 0.5 * law.body_force[a] - rho * law.skeleton_velocity[a]
        });
        let diagonal = law.principal.map(|d| rho + 0.5 * d);
        if !rhs.iter().chain(&diagonal).all(|v| v.is_finite()) {
            return None;
        }
        let c = 0.5 * self.quadratic;
        let drag = if law.principal.iter().all(|&d| d == law.principal[0]) {
            // Isotropic linear plus vector-speed drag is one scalar quadratic.
            // Avoid rotating an isotropic law through an arbitrary material frame.
            let speed = positive_root(diagonal[0], c, norm(rhs))?;
            let denominator = diagonal[0] + c * speed;
            if !denominator.is_finite() {
                return None;
            }
            rhs.map(|q| {
                let v = q / denominator;
                -law.principal[0] * v - (self.quadratic * speed) * v
            })
        } else {
            let q = law.axes.map(|axis| dot(axis, rhs));
            let speed = anisotropic_speed(diagonal, c, q)?;
            let principal = core::array::from_fn::<_, 3, _>(|a| {
                let v = q[a] / (diagonal[a] + c * speed);
                -law.principal[a] * v - (self.quadratic * speed) * v
            });
            core::array::from_fn(|a| (0..3).map(|k| law.axes[k][a] * principal[k]).sum())
        };
        let force = core::array::from_fn(|a| law.body_force[a] + drag[a]);
        force
            .iter()
            .chain(&drag)
            .all(|v| v.is_finite())
            .then_some((force, drag))
    }

    pub(super) fn nonlinear_dissipation(self, initial: [f64; 3], end: [f64; 3]) -> f64 {
        let v0 = core::array::from_fn(|a| initial[a] - self.linear.skeleton_velocity[a]);
        let v1 = core::array::from_fn(|a| end[a] - self.linear.skeleton_velocity[a]);
        let (r0, r1) = (norm(v0), norm(v1));
        let midpoint = core::array::from_fn(|a| 0.5 * v0[a] + 0.5 * v1[a]);
        let half_difference = 0.5 * r0 - 0.5 * r1;
        // Exactly the trapezoidal drag dotted with midpoint relative velocity:
        // beta/8 (r0+r1) (|v0+v1|² + (r0-r1)²). Every term is nonnegative.
        // beta |v_mid|³ is a different time quadrature and is not substituted.
        (self.quadratic * (0.5 * r0 + 0.5 * r1))
            * (dot(midpoint, midpoint) + half_difference * half_difference)
    }
}

// Stable positive solution of a*s + c*s² = q, with a>0 and c,q>=0.
fn positive_root(a: f64, c: f64, q: f64) -> Option<f64> {
    if !q.is_finite() {
        return None;
    }
    if q == 0.0 {
        return Some(0.0);
    }
    let half_a = 0.5 * a;
    let root = half_a.hypot(c.sqrt() * q.sqrt());
    let denominator = half_a + root;
    let speed = if denominator.is_finite() {
        q / denominator
    } else {
        (0.5 * q) / (0.5 * half_a + 0.5 * root)
    };
    (speed.is_finite() && speed > 0.0).then_some(speed)
}

fn anisotropic_speed(diagonal: [f64; 3], c: f64, q: [f64; 3]) -> Option<f64> {
    // If s_i solves a_i*s_i+c*s_i²=|q_i|, then |v_i|<=s_i<=|v|.
    // Thus max(s_i)<=|v|<=||s||, a bracket spanning at most sqrt(3),
    // independent of linear-resistance contrast or nonlinear stiffness.
    let mut component = [0.0; 3];
    for a in 0..3 {
        component[a] = positive_root(diagonal[a], c, q[a].abs())?;
    }
    let mut lower = component.into_iter().fold(0.0_f64, f64::max);
    let mut upper = norm(component);
    if !upper.is_finite() {
        return None;
    }
    if upper == 0.0 {
        return Some(0.0);
    }
    let mut speed = upper;
    for _ in 0..64 {
        let denominator = diagonal.map(|a| a + c * speed);
        if !denominator.iter().all(|v| v.is_finite()) {
            return None;
        }
        let velocity = core::array::from_fn(|a| q[a] / denominator[a]);
        let length = norm(velocity);
        let residual = speed - length;
        if residual.abs() <= 64.0 * f64::EPSILON * speed.max(length) {
            return Some(speed);
        }
        if residual > 0.0 {
            upper = speed;
        } else {
            lower = speed;
        }
        let derivative = 1.0
            + (0..3)
                .map(|a| {
                    let direction = velocity[a] / length;
                    direction * direction * ((c * length) / denominator[a])
                })
                .sum::<f64>();
        let newton = speed - residual / derivative;
        speed = if newton.is_finite() && newton > lower && newton < upper {
            newton
        } else {
            lower + 0.5 * (upper - lower)
        };
    }
    None
}

impl BoundaryGrid3 {
    /// Advance with implicit anisotropic Darcy and vector-speed Forchheimer drag.
    ///
    /// Uses the exact [`Self::step_with_linear_resistance`] owner, force-frame
    /// preparation, collision, pull stream and atomic publication. Isotropic
    /// linear resistance has an analytic endpoint; anisotropy uses one monotone
    /// speed root, with a tight bracket, safeguarded Newton and a 64-iteration
    /// cap. An unrepresentable or unconverged endpoint refuses atomically.
    /// Zero quadratic resistance reproduces the existing linear step exactly.
    ///
    /// The returned source work uses endpoint-averaged constitutive drag and
    /// midpoint relative velocity. Its nonnegative nonlinear dissipation is
    /// `beta/8 (r0+r1) (|v0+v1|²+(r0-r1)²)`, not `beta |v_mid|³`.
    /// Uniform resolved decay is second order. Stiff trapezoidal updates can
    /// reverse velocity; this is not an L-stable/exact relaxation integrator.
    /// Prescribed-skeleton, no-heat and native low-Mach limits
    /// are unchanged. This overload measures no wall links and uses zero origin.
    pub fn step_with_darcy_forchheimer(
        &mut self,
        resistance: impl FnMut([usize; 3]) -> DarcyForchheimer3,
        gate: &CancelGate,
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_resistance_inner::<fs_exec::TilePool>(None, gate, resistance, &[], [0.0; 3])
    }

    /// Same nonlinear step retaining the native exchange on selected walls.
    /// Skeleton, independent force and wall moments use `moment_origin`.
    pub fn step_with_darcy_forchheimer_and_wall_momentum(
        &mut self,
        resistance: impl FnMut([usize; 3]) -> DarcyForchheimer3,
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_resistance_inner::<fs_exec::TilePool>(
            None,
            gate,
            resistance,
            links,
            moment_origin,
        )
    }

    /// Pooled nonlinear step with identical numerical kernels and ordered receipts.
    pub fn step_pooled_with_darcy_forchheimer<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        resistance: impl FnMut([usize; 3]) -> DarcyForchheimer3,
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_resistance_inner(Some(pool), gate, resistance, &[], [0.0; 3])
    }

    /// Pooled nonlinear step with selected native wall impulse, torque and work.
    pub fn step_pooled_with_darcy_forchheimer_and_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        resistance: impl FnMut([usize; 3]) -> DarcyForchheimer3,
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
    ) -> Result<ResistanceMomentumExchange3, BoundaryStepError3> {
        self.step_resistance_inner(Some(pool), gate, resistance, links, moment_origin)
    }
}

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
