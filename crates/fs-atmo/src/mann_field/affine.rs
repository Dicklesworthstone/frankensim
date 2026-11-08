//! Local affine observations with explicit, realization-specific error budgets.
//!
//! The retained field supplies the value and analytic gradient at one authored
//! anchor. For a homogeneous mode `a cos(theta) + b sin(theta)`, let
//! `A = sqrt(|a|^2 + |b|^2)`. Its velocity Hessian and gradient Lipschitz
//! constants are bounded by `A |k|^2`. For the harmonic wall correction let
//! `R = hypot(a_z, b_z)` and `q = hypot(k_x, k_y)`. Its corresponding bounds
//! throughout the upper half-space are `R q^2` and `sqrt(2) R q^2`: both the
//! horizontal phase derivatives and the vertical exponential derivatives enter.
//! Thus a displacement of length r has velocity error <= H_v r^2/2 and
//! gradient Frobenius error <= H_g r. The modal constants use the actual OU
//! quadratures, not an ensemble RMS or a continuum spectral approximation.
//!
//! These are analytical bounds evaluated in ordinary floating-point arithmetic,
//! not interval certificates. They bound spatial truncation of the represented
//! realization, not roundoff, omitted spectrum, or temporal/trajectory error.
//! An observation is an immutable snapshot; rebuild it after changing time or
//! amplitudes. Preparation costs O(modes), uses fixed-size storage, and draws no
//! randomness. Each admitted observation then costs O(1).

use super::{FlowSample, MannPlaneWall, MannSpectralModel, Refusal, det, refuse};

/// Authored absolute spatial truncation tolerances; zero admits only a zero
/// estimated remainder. The gradient tolerance uses the Frobenius norm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannAffineAccuracy {
    pub maximum_velocity_error_mps: f64,
    pub maximum_gradient_error_per_s: f64,
}

/// Analytical upper-bound estimates for one point's spatial truncation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MannAffineError {
    pub velocity_mps: f64,
    pub gradient_per_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannAffineSample {
    pub flow: FlowSample,
    pub error: MannAffineError,
}

/// One fixed anchor and one accepted realization, including its physical time.
/// The affine gradient retains the source's zero trace up to roundoff. A wall
/// observation must be anchored on the plane, so its normal velocity and
/// tangential derivatives vanish exactly there. No piecewise patches, clipped
/// coordinates, or query-order-dependent anchor selection are introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannAffineField {
    anchor_m: [f64; 3],
    time_s: f64,
    wall: Option<MannPlaneWall>,
    anchor: FlowSample,
    velocity_curvature: f64,
    gradient_curvature: f64,
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use a finite fixed anchor, the same realization and explicit representable error budgets; use exact sampling when the affine budget is insufficient",
    )
}

fn numerical() -> Refusal {
    error(
        "mann-affine-numerical",
        "affine field or a nonzero remainder bound is not representable",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-affine-cancelled",
            "affine preparation cancelled",
        ))
    } else {
        Ok(())
    }
}

// A zero caused by multiplying two positive values is a refusal, not evidence
// of an exact field. Chained products avoid forming |k|^2 or r^2 prematurely.
fn product(a: f64, b: f64) -> Result<f64, Refusal> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        return Err(numerical());
    }
    Ok(value)
}

fn norm(values: impl IntoIterator<Item = f64>) -> Result<f64, Refusal> {
    let mut result = 0.0;
    let mut nonzero = false;
    for value in values {
        nonzero |= value != 0.0;
        result = det::hypot(result, value);
    }
    if !result.is_finite() || (nonzero && result == 0.0) {
        return Err(numerical());
    }
    Ok(result)
}

impl MannAffineAccuracy {
    pub fn admit(self) -> Result<(), Refusal> {
        if !self.maximum_velocity_error_mps.is_finite()
            || self.maximum_velocity_error_mps < 0.0
            || !self.maximum_gradient_error_per_s.is_finite()
            || self.maximum_gradient_error_per_s < 0.0
        {
            return Err(error(
                "mann-affine-accuracy",
                "affine error budgets must be finite and nonnegative",
            ));
        }
        Ok(())
    }
}

impl MannAffineField {
    pub const fn anchor_m(&self) -> [f64; 3] {
        self.anchor_m
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn wall(&self) -> Option<MannPlaneWall> {
        self.wall
    }

    /// Evaluate the same global affine field at every admitted point. The
    /// exact source value/gradient are returned unchanged at the anchor.
    /// A point outside either error budget refuses; there is no silent exact
    /// fallback or local reanchoring that could depend on query order.
    pub fn sample(
        &self,
        point_m: [f64; 3],
        accuracy: MannAffineAccuracy,
    ) -> Result<MannAffineSample, Refusal> {
        accuracy.admit()?;
        if !point_m.iter().all(|x| x.is_finite())
            || self.wall.is_some_and(|wall| point_m[2] < wall.height_m())
        {
            return Err(error(
                "mann-affine-domain",
                "affine point must be finite and above the declared plane",
            ));
        }
        let displacement: [f64; 3] = core::array::from_fn(|i| point_m[i] - self.anchor_m[i]);
        let radius = norm(displacement)?;
        let error = MannAffineError {
            velocity_mps: product(
                product(product(self.velocity_curvature, radius)?, radius)?,
                0.5,
            )?,
            gradient_per_s: product(self.gradient_curvature, radius)?,
        };
        if error.velocity_mps > accuracy.maximum_velocity_error_mps
            || error.gradient_per_s > accuracy.maximum_gradient_error_per_s
        {
            return Err(self::error(
                "mann-affine-accuracy",
                "point exceeds the authored affine spatial error budget",
            ));
        }
        let mut flow = self.anchor;
        if radius != 0.0 {
            for (velocity, gradient) in flow.u.iter_mut().zip(flow.grad) {
                *velocity += gradient[0] * displacement[0]
                    + gradient[1] * displacement[1]
                    + gradient[2] * displacement[2];
            }
        }
        if !flow.u.iter().all(|x| x.is_finite()) {
            return Err(numerical());
        }
        Ok(MannAffineSample { flow, error })
    }
}

pub(super) fn prepare_modes(
    model: &MannSpectralModel,
    time_s: f64,
    amplitudes: [&[f64]; 4],
    anchor_m: [f64; 3],
    wall: Option<MannPlaneWall>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MannAffineField, Refusal> {
    poll(&mut cancelled)?;
    if !anchor_m.iter().all(|x| x.is_finite())
        || wall.is_some_and(|wall| anchor_m[2] != wall.height_m())
    {
        return Err(error(
            "mann-affine-anchor",
            "anchor must be finite and lie on the declared slip plane",
        ));
    }
    let anchor = match wall {
        Some(wall) => super::wall::sample_modes(model, time_s, amplitudes, anchor_m, wall),
        None => super::sample_modes(model, time_s, amplitudes, anchor_m, None),
    }?;
    let mut velocity_curvature = 0.0;
    let mut gradient_curvature = 0.0;
    for (index, mode) in model.modes.iter().enumerate() {
        poll(&mut cancelled)?;
        let mut a = [0.0; 3];
        let mut b = [0.0; 3];
        for i in 0..3 {
            a[i] = product(mode.polarization[0][i], amplitudes[0][index])?
                + product(mode.polarization[1][i], amplitudes[2][index])?;
            b[i] = product(mode.polarization[0][i], amplitudes[1][index])?
                + product(mode.polarization[1][i], amplitudes[3][index])?;
        }
        let amplitude = norm(a.into_iter().chain(b))?;
        let k = norm(mode.k)?;
        let homogeneous = product(product(amplitude, k)?, k)?;
        let blocked = if wall.is_some() {
            let normal = norm([a[2], b[2]])?;
            let q = norm([mode.k[0], mode.k[1]])?;
            product(product(normal, q)?, q)?
        } else {
            0.0
        };
        velocity_curvature += homogeneous + blocked;
        gradient_curvature += homogeneous + product(core::f64::consts::SQRT_2, blocked)?;
    }
    if !velocity_curvature.is_finite() || !gradient_curvature.is_finite() {
        return Err(numerical());
    }
    poll(&mut cancelled)?;
    Ok(MannAffineField {
        anchor_m,
        time_s,
        wall,
        anchor,
        velocity_curvature,
        gradient_curvature,
    })
}

#[cfg(test)]
mod tests;
