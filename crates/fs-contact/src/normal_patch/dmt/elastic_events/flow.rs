//! Entire solution of the constant-sign, linearly driven scalar oscillator.
//! Cells obey h*(sqrt(k/m)+c/m)<=1/4. Velocity extrema therefore partition
//! every cell into monotone intervals: underdamped acceleration zeros are at
//! least pi/sqrt(k/m) apart; critical/overdamped acceleration has at most one.
use super::{DmtElasticEventError as Error, finite};

const ORDER: usize = 24;

pub(super) struct Flow {
    velocity: [f64; ORDER + 1],
    pub h: f64,
    pub truncation_speed: f64,
}

impl Flow {
    pub fn new(
        v: f64,
        acceleration: f64,
        stage_speed: f64,
        omega2: f64,
        drag: f64,
        h: f64,
    ) -> Result<Self, Error> {
        let w2 = finite(omega2 * h * h)?;
        let d = finite(drag * h)?;
        let mut velocity = [0.0; ORDER + 1];
        velocity[0] = v;
        velocity[1] = finite(acceleration * h)?;
        for j in 0..ORDER - 1 {
            let forcing = if j == 0 { stage_speed } else { 0.0 };
            velocity[j + 2] = finite(
                (w2 * (forcing - velocity[j]) - d * (j + 1) as f64 * velocity[j + 1])
                    / ((j + 1) * (j + 2)) as f64,
            )?;
        }
        // For [v-u,h*a], ||A||_infinity<=1 on this cell. The degree-24
        // exponential tail is < e/25! < 1.8e-25 times its initial norm.
        let truncation_speed =
            finite(1.8e-25 * (v - stage_speed).abs().max((acceleration * h).abs()))?;
        Ok(Self {
            velocity,
            h,
            truncation_speed,
        })
    }
    pub fn v(&self, z: f64) -> f64 {
        self.velocity
            .iter()
            .rev()
            .fold(0.0, |sum, coefficient| sum * z + coefficient)
    }
    /// h*acceleration, with no small-duration division.
    pub fn a(&self, z: f64) -> f64 {
        (1..=ORDER)
            .rev()
            .fold(0.0, |sum, i| sum * z + i as f64 * self.velocity[i])
    }
    /// h²*jerk, also a homogeneous damped-oscillator solution.
    pub fn jerk(&self, z: f64) -> f64 {
        (2..=ORDER).rev().fold(0.0, |sum, i| {
            sum * z + (i * (i - 1)) as f64 * self.velocity[i]
        })
    }
    pub fn dx(&self, z: f64) -> f64 {
        self.h
            * z
            * (0..=ORDER)
                .rev()
                .fold(0.0, |sum, i| sum * z + self.velocity[i] / (i + 1) as f64)
    }
    /// Integral v*t, integral displacement, integral displacement*t.
    pub fn moments(&self, z: f64) -> [f64; 3] {
        let mut tv = 0.0;
        let mut x = 0.0;
        let mut tx = 0.0;
        for i in (0..=ORDER).rev() {
            tv = tv * z + self.velocity[i] / (i + 2) as f64;
            x = x * z + self.velocity[i] / ((i + 1) * (i + 2)) as f64;
            tx = tx * z + self.velocity[i] / ((i + 1) * (i + 3)) as f64;
        }
        [
            self.h * self.h * z * z * tv,
            self.h * self.h * z * z * x,
            self.h * self.h * self.h * z * z * z * tx,
        ]
    }
    pub fn speed_roundoff(&self) -> f64 {
        self.truncation_speed + 64.0 * f64::EPSILON * self.speed_bound()
    }
    pub fn speed_bound(&self) -> f64 {
        self.velocity.iter().map(|x| x.abs()).sum::<f64>()
    }
    /// Positive quadrature of c*v²; never a mechanical energy residual.
    pub fn damping_loss(&self, c: f64, z: f64, stage_speed: f64) -> Result<(f64, f64), Error> {
        if c == 0.0 {
            return Ok((0.0, 0.0));
        }
        const GAUSS: [(f64, f64); 8] = [
            (0.019855071751231884, 0.05061426814518813),
            (0.10166676129318664, 0.11119051722668724),
            (0.2372337950418355, 0.15685332293894366),
            (0.4082826787521751, 0.181341891689181),
            (0.5917173212478248, 0.181341891689181),
            (0.7627662049581645, 0.15685332293894366),
            (0.8983332387068134, 0.11119051722668724),
            (0.9801449282487681, 0.05061426814518813),
        ];
        let square = GAUSS
            .iter()
            .map(|(node, weight)| weight * self.v(z * node).powi(2))
            .sum::<f64>();
        let loss = finite(c * self.h * z * square)?;
        let norm = (self.velocity[0] - stage_speed)
            .abs()
            .max(self.velocity[1].abs());
        // ||A||∞<=1: the Gauss-8 remainder constant is 1.698e-23.
        // Bounding the 16th derivative of (u+e_1 exp(A*t)y)^2 gives
        // <1e-16*(|u|+||y||∞)^2. Rounding is separately retained.
        let quadrature = c * self.h * z * 1e-16 * (stage_speed.abs() + norm).powi(2);
        let tail = c
            * self.h
            * z
            * self.truncation_speed
            * (2.0 * (stage_speed.abs() + 3.0 * norm) + self.truncation_speed);
        Ok((
            loss,
            finite(quadrature + tail + 256.0 * f64::EPSILON * loss)?,
        ))
    }
}
