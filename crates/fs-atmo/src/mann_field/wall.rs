//! Harmonic blocking of the represented Mann field above one plane.
//!
//! Hunt & Graham, *Free-stream turbulence near plane boundaries*, JFM 84
//! (1978), 209–235, <https://doi.org/10.1017/S0022112078000130>, separate
//! the outer irrotational blocking correction from the thin viscous layer.
//! For a Fourier amplitude `a`, horizontal wavenumber `q = hypot(kx, ky)`
//! and height `h = z - wall_z`, this module adds
//! `a_z (i kx/q, i ky/q, -1) exp(i k.wall_point - i k.U t - q h)`.
//! It is the gradient of a harmonic potential: normal velocity vanishes at
//! the plane, divergence is preserved, and the correction decays aloft.
//! The original two polarizations, quadratures and OU histories are retained.
//!
//! This is a slip-plane, outer-layer model, not a no-slip viscous layer,
//! wall-generated turbulence or a calibration of the near-ground spectrum.
//! Wall blocking and homogeneous shear are distinct effects in Mann (1994),
//! JFM 273, <https://doi.org/10.1017/S0022112094001886>.
//! [`filter`] integrates the same field over finite boxes, faces and lines;
//! the evanescent correction has its own analytic observation gain.

use super::*;

/// Exact finite-footprint observations of the wall-blocked field.
pub mod filter;

/// Plane `z = height_m` in the authored Mann coordinate frame.
/// Only points in its exact upper half-space are admitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannPlaneWall {
    height_m: f64,
}

impl MannPlaneWall {
    /// A finite plane height; negative heights are allowed.
    pub fn new(height_m: f64) -> Result<Self, Refusal> {
        if !height_m.is_finite() {
            return Err(invalid("Mann plane height must be finite"));
        }
        Ok(Self { height_m })
    }

    #[must_use]
    pub const fn height_m(self) -> f64 {
        self.height_m
    }

    fn distance(self, point_m: [f64; 3]) -> Result<f64, Refusal> {
        finite_vector(point_m)?;
        if point_m[2] < self.height_m {
            return Err(invalid("Mann wall sample lies below the plane"));
        }
        let distance = point_m[2] - self.height_m;
        if !distance.is_finite() {
            return Err(invalid("Mann wall distance is not representable"));
        }
        Ok(distance)
    }
}

// Geometry is shared by realization sampling and transformed-tensor covariance.
// Fixed-size stack values only; no additional modal arrays or random paths.
struct ModePoint {
    k: [f64; 3],
    q: f64,
    horizontal: [f64; 2],
    phase: [f64; 2],      // cosine, sine
    wall_phase: [f64; 2], // exp(-q h) times cosine, sine at the plane
    difference: [f64; 2],
}

impl ModePoint {
    fn new(
        mode: &SpectralMode,
        advection_mps: [f64; 3],
        time_s: f64,
        point_m: [f64; 3],
        wall: MannPlaneWall,
        height: f64,
    ) -> Result<Self, Refusal> {
        let k = mode.k;
        let position = core::array::from_fn(|i| point_m[i] - advection_mps[i] * time_s);
        let phase = dot(k, position);
        let q = det::hypot(k[0], k[1]);
        if !phase.is_finite() || !q.is_finite() {
            return Err(invalid("Mann wall sample phase is not representable"));
        }
        // A purely vertical solenoidal mode has no normal amplitude and
        // needs no correction. Do not replace q by an epsilon.
        if q == 0.0 {
            if mode.polarization.iter().any(|p| p[2] != 0.0) {
                return Err(invalid("vertical Mann mode has a normal amplitude"));
            }
            let phase = [det::cos(phase), det::sin(phase)];
            return Ok(Self {
                k,
                q,
                horizontal: [0.0; 2],
                phase,
                wall_phase: [0.0; 2],
                difference: phase,
            });
        }
        let attenuation = q * height;
        let decay = det::exp(-attenuation);
        if decay == 0.0 {
            // The represented correction has vanished. Do not evaluate an
            // irrelevant distant-plane phase that could overflow.
            let phase = [det::cos(phase), det::sin(phase)];
            return Ok(Self {
                k,
                q,
                horizontal: [k[0] / q, k[1] / q],
                phase,
                wall_phase: [0.0; 2],
                difference: phase,
            });
        }
        let wall_position = [
            position[0],
            position[1],
            wall.height_m - advection_mps[2] * time_s,
        ];
        let wall_angle = dot(k, wall_position);
        if !wall_angle.is_finite() {
            return Err(invalid("Mann wall trace phase is not representable"));
        }
        let wall_trig = [det::cos(wall_angle), det::sin(wall_angle)];
        let wall_phase = wall_trig.map(|v| decay * v);
        let delta = k[2] * height;
        let (phase, difference) = if delta.abs() <= 0.5 && attenuation <= 0.5 {
            // Evaluate e^(i delta) - e^(-q h) without subtracting nearly
            // equal cosines. This also makes w and its tangential derivatives
            // exactly zero at h=0, including a large common advection phase.
            let sin_delta = det::sin(delta);
            let sin_half = det::sin(0.5 * delta);
            let cos_minus_one = -2.0 * sin_half * sin_half;
            let cos_minus_decay = cos_minus_one - det::expm1(-attenuation);
            (
                [
                    wall_trig[0] * (1.0 + cos_minus_one) - wall_trig[1] * sin_delta,
                    wall_trig[1] * (1.0 + cos_minus_one) + wall_trig[0] * sin_delta,
                ],
                [
                    wall_trig[0] * cos_minus_decay - wall_trig[1] * sin_delta,
                    wall_trig[1] * cos_minus_decay + wall_trig[0] * sin_delta,
                ],
            )
        } else {
            let phase = [det::cos(phase), det::sin(phase)];
            (phase, core::array::from_fn(|i| phase[i] - wall_phase[i]))
        };
        Ok(Self {
            k,
            q,
            horizontal: [k[0] / q, k[1] / q],
            phase,
            wall_phase,
            difference,
        })
    }

    fn values(&self, a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        let wall_derivative = -a[2] * self.wall_phase[1] + b[2] * self.wall_phase[0];
        [
            a[0] * self.phase[0] + b[0] * self.phase[1] + self.horizontal[0] * wall_derivative,
            a[1] * self.phase[0] + b[1] * self.phase[1] + self.horizontal[1] * wall_derivative,
            a[2] * self.difference[0] + b[2] * self.difference[1],
        ]
    }

    fn add_sample(&self, a: [f64; 3], b: [f64; 3], sample: &mut FlowSample) {
        let value = self.values(a, b);
        let derivative: [f64; 3] =
            core::array::from_fn(|i| -a[i] * self.phase[1] + b[i] * self.phase[0]);
        let wall_value = a[2] * self.wall_phase[0] + b[2] * self.wall_phase[1];
        let wall_derivative = -a[2] * self.wall_phase[1] + b[2] * self.wall_phase[0];
        let normal_derivative = -a[2] * self.difference[1] + b[2] * self.difference[0];
        for i in 0..3 {
            sample.u[i] += value[i];
        }
        for i in 0..2 {
            for j in 0..2 {
                sample.grad[i][j] += (derivative[i] - self.horizontal[i] * wall_value) * self.k[j];
            }
            sample.grad[i][2] +=
                derivative[i] * self.k[2] - self.horizontal[i] * wall_derivative * self.q;
        }
        for j in 0..2 {
            sample.grad[2][j] += normal_derivative * self.k[j];
        }
        sample.grad[2][2] += derivative[2] * self.k[2] + self.q * wall_value;
    }
}

pub(super) fn sample_modes(
    model: &MannSpectralModel,
    time_s: f64,
    amplitudes: [&[f64]; 4],
    point_m: [f64; 3],
    wall: MannPlaneWall,
) -> Result<FlowSample, Refusal> {
    let height = wall.distance(point_m)?;
    let mut sample = FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    };
    for (index, mode) in model.modes.iter().enumerate() {
        let geometry = ModePoint::new(mode, model.advection_mps, time_s, point_m, wall, height)?;
        let a = core::array::from_fn(|i| {
            mode.polarization[0][i] * amplitudes[0][index]
                + mode.polarization[1][i] * amplitudes[2][index]
        });
        let b = core::array::from_fn(|i| {
            mode.polarization[0][i] * amplitudes[1][index]
                + mode.polarization[1][i] * amplitudes[3][index]
        });
        geometry.add_sample(a, b, &mut sample);
    }
    if !sample
        .u
        .iter()
        .chain(sample.grad.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(invalid("Mann wall sample is not representable"));
    }
    Ok(sample)
}

pub(super) fn covariance_modes(
    model: &MannSpectralModel,
    correlation_times_s: Option<&[f64]>,
    first_m: [f64; 3],
    second_m: [f64; 3],
    lag_s: f64,
    wall: MannPlaneWall,
) -> Result<[[f64; 3]; 3], Refusal> {
    let first_height = wall.distance(first_m)?;
    let second_height = wall.distance(second_m)?;
    if !lag_s.is_finite() {
        return Err(invalid("Mann wall covariance lag must be finite"));
    }
    // Remove the first point's common spatial phase, translating the plane
    // with it. This retains both wall distances and avoids cancellation of
    // large common phases when two nearby points are far above the plane.
    let first = [0.0; 3];
    let second = core::array::from_fn(|i| second_m[i] - first_m[i]);
    finite_vector(second)?;
    let origin_wall = MannPlaneWall {
        height_m: -first_height,
    };
    let mut covariance = [[0.0; 3]; 3];
    for (index, mode) in model.modes.iter().enumerate() {
        let a = ModePoint::new(
            mode,
            model.advection_mps,
            0.0,
            first,
            origin_wall,
            first_height,
        )?;
        let b = ModePoint::new(
            mode,
            model.advection_mps,
            lag_s,
            second,
            origin_wall,
            second_height,
        )?;
        let decay = correlation_times_s.map_or(1.0, |times| det::exp(-lag_s.abs() / times[index]));
        for &polarization in &mode.polarization {
            let ac = a.values(polarization, [0.0; 3]);
            let as_ = a.values([0.0; 3], polarization);
            let bc = b.values(polarization, [0.0; 3]);
            let bs = b.values([0.0; 3], polarization);
            for i in 0..3 {
                for j in 0..3 {
                    covariance[i][j] += decay * (ac[i] * bc[j] + as_[i] * bs[j]);
                }
            }
        }
    }
    finite_covariance(covariance)
}

#[cfg(test)]
mod tests;
