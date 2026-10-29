//! Positive viscous-work integrals for the same exact scalar transition.
use super::{ExactZohCoefficients, ModalAcousticMode, ModalAcousticState};
use fs_math::det;

/// The existing entire-series transition, shared with the heat preparation.
/// All fields are dimensionless on the current preparation interval.
pub(super) struct DampedTransition {
    pub angle: f64,
    pub drag: f64,
    pub b: f64,
    pub g: f64,
    pub position_retention: f64,
    pub velocity_retention: f64,
    pub position_decrement: f64,
    pub velocity_decrement: f64,
}
impl DampedTransition {
    pub fn small(angle: f64, drag: f64) -> Self {
        let stiffness = angle * angle;
        let mut previous = 0.0;
        let mut power = 1.0;
        let mut reciprocal_factorial = 1.0;
        let mut b = 1.0;
        let mut g = 0.5;
        // p_j = -drag*p_{j-1} - angle^2*p_{j-2};
        // B/h = sum p_j/(j+1)!, G/h^2 = sum p_j/(j+2)!.
        for j in 1..=18 {
            let next = -drag * power - stiffness * previous;
            previous = power;
            power = next;
            reciprocal_factorial /= f64::from(j + 1);
            b += power * reciprocal_factorial;
            g += power * reciprocal_factorial / f64::from(j + 2);
        }
        let position_decrement = stiffness * g;
        let velocity_decrement = position_decrement + drag * b;
        Self {
            angle,
            drag,
            b,
            g,
            position_decrement,
            velocity_decrement,
            position_retention: 1.0 - position_decrement,
            velocity_retention: 1.0 - velocity_decrement,
        }
    }

    pub fn double(&mut self) {
        let cross = (self.angle * self.b) * (self.angle * self.b);
        let next_position = self.position_retention * self.position_retention - cross;
        let next_velocity = self.velocity_retention * self.velocity_retention - cross;
        self.g = 0.25 * (self.g * (1.0 + self.position_retention) + self.b * self.b);
        self.b *= 0.5 * (self.position_retention + self.velocity_retention);
        self.angle *= 2.0;
        self.drag *= 2.0;
        self.position_decrement = (self.angle * self.angle) * self.g;
        self.velocity_decrement = self.position_decrement + self.drag * self.b;
        // Keep the original separate increments and retentions: subtraction
        // cannot recover both small attacks and exponentially small tails.
        self.position_retention = if self.position_decrement.abs() < 0.5 {
            1.0 - self.position_decrement
        } else {
            next_position
        };
        self.velocity_retention = if self.velocity_decrement.abs() < 0.5 {
            1.0 - self.velocity_decrement
        } else {
            next_velocity
        };
    }

    pub fn coefficients(&self, dt: f64) -> ExactZohCoefficients {
        ExactZohCoefficients::Damped {
            position_retention: self.position_retention,
            velocity_retention: self.velocity_retention,
            velocity_gain: dt * self.b,
            displacement_gain: dt * (dt * self.g),
            position_decrement: self.position_decrement,
            velocity_decrement: self.velocity_decrement,
        }
    }
}

/// Upper-triangular square root of the velocity observability Gramian.
/// Keeping a square root avoids subtractive determinants or clipping heat.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct ViscousIntegral {
    r00: f64,
    r01: f64,
    r11: f64,
    roundoff_factor: f64,
}
impl ViscousIntegral {
    pub fn new(mode: ModalAcousticMode, dt: f64) -> Self {
        if mode.damping_ratio == 0.0 || mode.angular_frequency_rad_s == 0.0 {
            return Self::default();
        }
        let original_angle = mode.angular_frequency_rad_s * dt;
        let mut angle = original_angle;
        let mut drag = 2.0 * (mode.damping_ratio * angle);
        if !angle.is_finite() || !drag.is_finite() {
            return Self {
                r00: f64::NAN,
                ..Self::default()
            };
        }
        // Use [v, (F-omega^2*q)/kappa], kappa=max(omega,1/dt).
        // This keeps very soft forced modes representable without F/omega
        // or a Gramian entry proportional to omega^3 underflowing first.
        let ratio = original_angle.min(1.0);
        let mut scaled_time = original_angle.max(1.0);
        let mut doublings = 0;
        while angle.abs().max(drag.abs()).max(scaled_time) > 0.25 {
            angle *= 0.5;
            drag *= 0.5;
            scaled_time *= 0.5;
            doublings += 1;
        }
        let mut flow = DampedTransition::small(angle, drag);
        let mut integral = Self::default();
        // Eight-point Gauss-Legendre on [0,1], exclusively on the scaled
        // interval. Its dimensionless generator has both induced 1/inf norms
        // <=1/2. The 16th derivative bound gives a quadrature remainder
        // <4.7e-23*drag*||[v,r]||^2, below f64 roundoff. Long/stiff/oscillatory
        // intervals compose this resolved integral; they are never sampled
        // by eight points over an unresolved full interval.
        const QUADRATURE: [(f64, f64); 8] = [
            (0.019855071751231884, 0.05061426814518813),
            (0.10166676129318664, 0.11119051722668724),
            (0.2372337950418355, 0.15685332293894366),
            (0.4082826787521751, 0.181341891689181),
            (0.5917173212478248, 0.181341891689181),
            (0.7627662049581645, 0.15685332293894366),
            (0.8983332387068134, 0.11119051722668724),
            (0.9801449282487681, 0.05061426814518813),
        ];
        for (node, weight) in QUADRATURE {
            let at = DampedTransition::small(angle * node, drag * node);
            let weight = det::sqrt(drag) * det::sqrt(weight);
            integral.append_row(
                weight * at.velocity_retention,
                weight * scaled_time * node * at.b,
            );
        }
        for _ in 0..doublings {
            // H(2h)=H(h)+E(h)^T H(h) E(h), as QR of stacked [R; R*E].
            let e01 = scaled_time * flow.b;
            let e10 = -(ratio * ratio) * e01;
            let a = integral.r00 * flow.velocity_retention + integral.r01 * e10;
            let b = integral.r00 * e01 + integral.r01 * flow.position_retention;
            let c = integral.r11 * e10;
            let d = integral.r11 * flow.position_retention;
            integral.append_row(a, b);
            integral.append_row(c, d);
            flow.double();
            scaled_time *= 2.0;
        }
        integral.roundoff_factor = (256.0 + 128.0 * f64::from(doublings)) * f64::EPSILON;
        integral
    }

    fn append_row(&mut self, a: f64, b: f64) {
        let r = hypot(self.r00, a);
        if r == 0.0 {
            self.r11 = hypot(self.r11, b);
            return;
        }
        let c = self.r00 / r;
        let s = a / r;
        let residual = -s * self.r01 + c * b;
        self.r01 = c * self.r01 + s * b;
        self.r00 = r;
        self.r11 = hypot(self.r11, residual);
    }

    pub fn evaluate(
        self,
        mode: ModalAcousticMode,
        state: ModalAcousticState,
        force: f64,
        dt: f64,
    ) -> (f64, f64) {
        if self.r00 == 0.0 && self.r01 == 0.0 && self.r11 == 0.0 {
            return (0.0, 0.0);
        }
        let omega = mode.angular_frequency_rad_s;
        let angle = omega * dt;
        let r = if angle > 1.0 {
            force / omega - omega * state.displacement_m_sqrt_kg
        } else {
            force * dt - (omega * state.displacement_m_sqrt_kg) * angle
        };
        let av = self.r00 * state.velocity_m_sqrt_kg_per_s;
        let br = self.r01 * r;
        let cr = self.r11 * r;
        let first = av + br;
        let heat = first * first + cr * cr;
        let magnitude = av.abs() + br.abs();
        let allowance =
            self.roundoff_factor * (magnitude * magnitude + cr * cr).max(f64::MIN_POSITIVE);
        (heat, allowance)
    }
}

fn hypot(a: f64, b: f64) -> f64 {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 {
        return 0.0;
    }
    let x = a / scale;
    let y = b / scale;
    scale * det::sqrt(x * x + y * y)
}

#[cfg(test)]
#[path = "heat/tests.rs"]
mod tests;
