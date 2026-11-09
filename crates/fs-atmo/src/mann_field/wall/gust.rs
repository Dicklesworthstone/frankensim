//! Compact turbulent gust packets on the original wall-blocked realization.
//!
//! For the existing solenoidal wall field `v`, fix the vector-potential gauge
//! `A = (integral_wall^z v_y ds, -integral_wall^z v_x ds, 0)`. Then
//! `curl A = v`, because `v_z(wall) = 0`, and both tangential components of
//! `A` vanish at the plane. Localize the POTENTIAL: `u = curl(e A)` equals
//! `e v + grad(e) cross A`. This retains zero divergence and wall penetration.
//! Multiplying the velocity by `e` alone retains zero wall penetration but
//! generally violates mass continuity through `grad(e) dot v`.
//! The analytic gradient includes the envelope Hessian and the potential's
//! first derivatives. Fourier vertical integrals are analytic, not quadrature.
//!
//! Potential-based inhomogeneous synthesis is discussed by Yu & Bai (2014),
//! JCP 256, 234–253, <https://doi.org/10.1016/j.jcp.2013.08.055>. The compact
//! envelope and wall gauge here are explicitly prescribed kinematics, not
//! that paper's calibrated construction. The packet changes the original
//! stationary Mann statistics, supplies no momentum/pressure evolution, and
//! makes no claim to reproduce a measured gust or a boundary-layer transition.

use super::*;

/// A compact C2 potential envelope moving with the source's uniform advection.
///
/// At physical time `t`, its center is `center + U_adv * (t-reference_time)`.
/// Within the axis-aligned box, `e = strength * product((1-s_i²)^3)`, where
/// `s_i = (x_i-center_i(t))/half_width_i`; outside, `e` and its first two
/// derivatives are zero. The resulting velocity is C1 across the box boundary.
/// All geometry is in the original upward Mann frame. The fixed wall clips
/// the admitted domain, with exact no penetration even if the box intersects it.
///
/// Strength scales the potential, not a prescribed peak velocity or RMS.
/// Velocity at the packet center is `strength` times the original wall field;
/// away from the center the curl correction is essential. Keep this immutable
/// definition with the source checkpoint to reconstruct the same gust.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannWallGust {
    wall: MannPlaneWall,
    center_m: [f64; 3],
    half_width_m: [f64; 3],
    reference_time_s: f64,
    strength: f64,
}

fn error(message: &str) -> Refusal {
    refuse(
        "mann-wall-gust-invalid",
        message.into(),
        "use finite geometry, positive representable half widths, and nonnegative strength",
    )
}

#[derive(Clone, Copy)]
struct Envelope {
    value: f64,
    grad: [f64; 3],
    hessian: [[f64; 3]; 3],
}

impl MannWallGust {
    /// Define a packet without allocating or drawing another realization.
    /// Half widths must have representable reciprocal squares. The source
    /// determines the clock, advection and all stochastic amplitudes.
    pub fn new(
        wall: MannPlaneWall,
        center_m: [f64; 3],
        half_width_m: [f64; 3],
        reference_time_s: f64,
        strength: f64,
    ) -> Result<Self, Refusal> {
        if !center_m.iter().all(|x| x.is_finite())
            || !half_width_m.iter().all(|&w| {
                let inverse_square = (1.0 / w) * (1.0 / w);
                w.is_finite() && w > 0.0 && inverse_square.is_finite() && inverse_square > 0.0
            })
            || !reference_time_s.is_finite()
            || !strength.is_finite()
            || strength < 0.0
        {
            return Err(error("gust definition is not representable"));
        }
        Ok(Self {
            wall,
            center_m,
            half_width_m,
            reference_time_s,
            strength,
        })
    }

    #[must_use]
    pub const fn wall(&self) -> MannPlaneWall {
        self.wall
    }

    /// Authored center at the reference time, in metres.
    #[must_use]
    pub const fn center_m(&self) -> [f64; 3] {
        self.center_m
    }

    #[must_use]
    pub const fn half_width_m(&self) -> [f64; 3] {
        self.half_width_m
    }

    #[must_use]
    pub const fn reference_time_s(&self) -> f64 {
        self.reference_time_s
    }

    #[must_use]
    pub const fn strength(&self) -> f64 {
        self.strength
    }

    fn envelope(
        &self,
        point_m: [f64; 3],
        time_s: f64,
        advection_mps: [f64; 3],
    ) -> Result<Option<Envelope>, Refusal> {
        let elapsed = time_s - self.reference_time_s;
        let offset: [f64; 3] =
            core::array::from_fn(|i| (point_m[i] - self.center_m[i]) - advection_mps[i] * elapsed);
        if !elapsed.is_finite() || !offset.iter().all(|x| x.is_finite()) {
            return Err(error("advected gust position is not representable"));
        }
        // Reject/skip before division: distant points cannot overflow their
        // normalized coordinate, and outside samples do no spectral work.
        if self.strength == 0.0 || (0..3).any(|i| offset[i].abs() >= self.half_width_m[i]) {
            return Ok(None);
        }
        let mut f = [0.0; 3];
        let mut d = [0.0; 3];
        let mut dd = [0.0; 3];
        for i in 0..3 {
            let inv = 1.0 / self.half_width_m[i];
            let s = offset[i] * inv;
            let r = (1.0 - s) * (1.0 + s);
            f[i] = r * r * r;
            d[i] = -6.0 * s * r * r * inv;
            dd[i] = 6.0 * r * (5.0 * s * s - 1.0) * inv * inv;
        }
        let value = self.strength * f[0] * f[1] * f[2];
        let mut grad = [0.0; 3];
        let mut hessian = [[0.0; 3]; 3];
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            grad[i] = self.strength * d[i] * f[j] * f[k];
            hessian[i][i] = self.strength * dd[i] * f[j] * f[k];
            hessian[i][j] = self.strength * d[i] * d[j] * f[k];
            hessian[j][i] = hessian[i][j];
        }
        if !value.is_finite()
            || !grad
                .iter()
                .chain(hessian.iter().flatten())
                .all(|x| x.is_finite())
        {
            return Err(error("gust envelope derivatives are not representable"));
        }
        Ok(Some(Envelope {
            value,
            grad,
            hessian,
        }))
    }
}

// Integrate one wall mode vertically, using the same phase/polarization
// convention as ModePoint. The returned FlowSample stores A and grad(A),
// whose units are m²/s and m/s, rather than velocity and velocity gradient.
fn potential(
    geometry: &ModePoint,
    height: f64,
    wall_angle: f64,
    a: [f64; 3],
    b: [f64; 3],
) -> Result<FlowSample, Refusal> {
    let half_delta = 0.5 * geometry.k[2] * height;
    if !half_delta.is_finite() {
        return Err(error("gust potential phase is not representable"));
    }
    let sinc = if half_delta.abs() < 1e-4 {
        let x2 = half_delta * half_delta;
        1.0 + x2 * (-1.0 / 6.0 + x2 * (1.0 / 120.0 - x2 / 5040.0))
    } else {
        det::sin(half_delta) / half_delta
    };
    let bulk_gain = height * sinc;
    let (sw, cw) = det::sin_cos(wall_angle);
    // Angle addition retains a small vertical phase next to a large common
    // advected phase, matching the original wall owner's near-plane policy.
    let (sh, ch) = det::sin_cos(half_delta);
    let s = sw * ch + cw * sh;
    let c = cw * ch - sw * sh;
    let attenuation = geometry.q * height;
    let wall_gain = if attenuation == 0.0 {
        height
    } else if attenuation <= 0.5 {
        height * (-det::expm1(-attenuation) / attenuation)
    } else {
        // expm1 retains the h -> 0 limit without subtracting rounded ones.
        -det::expm1(-attenuation) / geometry.q
    };
    let wall_value = a[2] * cw + b[2] * sw;
    let wall_derivative = -a[2] * sw + b[2] * cw;
    let at_point = geometry.values(a, b);
    let mut result = FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    };
    for i in 0..2 {
        let target = 1 - i;
        let sign = if i == 0 { -1.0 } else { 1.0 };
        result.u[target] = sign
            * ((a[i] * c + b[i] * s) * bulk_gain
                + geometry.horizontal[i] * wall_derivative * wall_gain);
        let horizontal_derivative = sign
            * ((-a[i] * s + b[i] * c) * bulk_gain
                - geometry.horizontal[i] * wall_value * wall_gain);
        for j in 0..2 {
            result.grad[target][j] = horizontal_derivative * geometry.k[j];
        }
        result.grad[target][2] = sign * at_point[i];
    }
    if !result
        .u
        .iter()
        .chain(result.grad.iter().flatten())
        .all(|x| x.is_finite())
    {
        return Err(error("gust vector potential is not representable"));
    }
    Ok(result)
}

fn localize(flow: FlowSample, potential: FlowSample, envelope: Envelope) -> FlowSample {
    let correction = cross(envelope.grad, potential.u);
    let mut result = FlowSample {
        u: core::array::from_fn(|i| envelope.value * flow.u[i] + correction[i]),
        grad: [[0.0; 3]; 3],
    };
    for j in 0..3 {
        let curvature = cross(
            core::array::from_fn(|i| envelope.hessian[i][j]),
            potential.u,
        );
        let connection = cross(
            envelope.grad,
            core::array::from_fn(|i| potential.grad[i][j]),
        );
        for i in 0..3 {
            result.grad[i][j] = envelope.value * flow.grad[i][j]
                + envelope.grad[j] * flow.u[i]
                + curvature[i]
                + connection[i];
        }
    }
    result
}

pub(in crate::mann_field) fn sample_modes(
    model: &MannSpectralModel,
    time_s: f64,
    amplitudes: [&[f64]; 4],
    point_m: [f64; 3],
    gust: &MannWallGust,
) -> Result<FlowSample, Refusal> {
    let height = gust.wall.distance(point_m)?;
    let zero = FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    };
    let Some(envelope) = gust.envelope(point_m, time_s, model.advection_mps)? else {
        return Ok(zero);
    };
    let wall_position = [
        point_m[0] - model.advection_mps[0] * time_s,
        point_m[1] - model.advection_mps[1] * time_s,
        gust.wall.height_m - model.advection_mps[2] * time_s,
    ];
    let mut flow = zero;
    let mut vector_potential = zero;
    for (index, mode) in model.modes.iter().enumerate() {
        let geometry = ModePoint::new(
            mode,
            model.advection_mps,
            time_s,
            point_m,
            gust.wall,
            height,
        )?;
        let wall_angle = dot(mode.k, wall_position);
        if !wall_angle.is_finite() {
            return Err(error("gust wall trace phase is not representable"));
        }
        let a = core::array::from_fn(|i| {
            mode.polarization[0][i] * amplitudes[0][index]
                + mode.polarization[1][i] * amplitudes[2][index]
        });
        let b = core::array::from_fn(|i| {
            mode.polarization[0][i] * amplitudes[1][index]
                + mode.polarization[1][i] * amplitudes[3][index]
        });
        geometry.add_sample(a, b, &mut flow);
        let part = potential(&geometry, height, wall_angle, a, b)?;
        for i in 0..3 {
            vector_potential.u[i] += part.u[i];
            for j in 0..3 {
                vector_potential.grad[i][j] += part.grad[i][j];
            }
        }
    }
    let result = localize(flow, vector_potential, envelope);
    if !result
        .u
        .iter()
        .chain(result.grad.iter().flatten())
        .all(|x| x.is_finite())
    {
        return Err(error(
            "localized gust velocity or gradient is not representable",
        ));
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
