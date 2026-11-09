//! Height-resolved Smith/GGX microsurface reflection geometry.
//!
//! This is the random walk of Heitz et al. (2016), sections 5--7:
//! https://jo.dreggn.org/home/2016_microfacets.pdf
//! We use depth `-ln(C1(height))`, so no arbitrary height distribution or
//! physical displacement is needed. A downward ray increases depth with rate
//! `1 + Lambda`; an upward ray decreases it with rate `Lambda` until escape.
//! Facets follow the FULL visible-normal distribution on both sides of the
//! macro horizon. The bounded reflection proposal used by the single-scatter
//! camera BSDF would incorrectly discard internal returns here.
//!
//! Geometry supplies probabilities, not attenuation: each collision's optical
//! weight belongs to the caller. No extra Smith G factor or albedo compensation
//! is appropriate. Returned densities describe individual hidden paths, not the
//! outgoing marginal BSDF proposal required by camera/BDPT MIS.

use core::f64::consts::PI;
use fs_geom::Vec3;
use fs_math::det;

pub mod connection;

/// Semantics of the explicit height and direction sampling operation tree.
pub const SMITH_GGX_WALK_SEMANTICS_VERSION: u32 = 1;

/// Invalid input or a floating-point representation failure, never absorption.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmithError {
    InvalidRoughness,
    InvalidDirection,
    InvalidDepth,
    InvalidSample,
    NumericalRange,
}
impl core::fmt::Display for SmithError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Smith GGX microsurface: {self:?}")
    }
}
impl core::error::Error for SmithError {}

/// A full, anisotropic GGX heightfield in a supplied local tangent frame.
/// The local +z axis points out of the opaque surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithGgx {
    alpha_x: f64,
    alpha_y: f64,
}

/// One geometric free-flight/phase sample. Direction is propagation direction,
/// not the BSDF convention pointing back toward the incident endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SmithEvent {
    /// The current direction leaves the microsurface, with its retained power.
    Escape,
    /// A real internal facet reflection, even if its direction remains below
    /// the macro horizon. The facet points toward local +z.
    Reflection {
        depth: f64,
        facet: Vec3,
        direction: Vec3,
    },
}

/// A sample in the mixed height/direction path measure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithStep {
    pub event: SmithEvent,
    /// Natural log of the conditional joint density: d(depth) d(omega_out)
    /// for a collision, or the discrete escape probability. At exact horizontal
    /// propagation, height is a Dirac atom and only d(omega_out) remains.
    /// This is NOT an outgoing marginal directional PDF.
    pub log_density: f64,
}

impl SmithGgx {
    /// Match the native conductor's admitted slope range without conflating
    /// a delta mirror with a very narrow continuous distribution.
    pub fn try_new(alpha_x: f64, alpha_y: f64) -> Result<Self, SmithError> {
        if [alpha_x, alpha_y]
            .iter()
            .any(|a| !a.is_finite() || !(1.0e-4..=1.0).contains(a))
        {
            return Err(SmithError::InvalidRoughness);
        }
        Ok(Self { alpha_x, alpha_y })
    }

    /// Sample the next interaction or escape from nonnegative CDF depth.
    /// Start an incident path at depth zero with direction.z < 0.
    /// `draw[0]` is an OPEN-unit-interval survival draw; the two full-VNDF
    /// coordinates are in [0,1). No retries or below-horizon rejections occur.
    pub fn sample(
        self,
        direction: Vec3,
        depth: f64,
        draw: [f64; 3],
    ) -> Result<SmithStep, SmithError> {
        if !depth.is_finite() || depth < 0.0 {
            return Err(SmithError::InvalidDepth);
        }
        if !draw[0].is_finite()
            || draw[0] <= 0.0
            || draw[0] >= 1.0
            || draw[1..]
                .iter()
                .any(|u| !u.is_finite() || !(0.0..1.0).contains(u))
        {
            return Err(SmithError::InvalidSample);
        }
        let direction = unit(direction)?;
        let wo = direction.scale(-1.0);
        let area = self.projected_area(wo);
        if direction.z > 0.0 && (area == 0.0 || depth == 0.0) {
            return Ok(SmithStep {
                event: SmithEvent::Escape,
                log_density: 0.0,
            });
        }
        if !area.is_finite() || area <= 0.0 {
            return Err(SmithError::NumericalRange);
        }
        let exponential = -det::ln(draw[0]);
        // Inverse rate avoids Lambda's division by a grazing cosine. At exact
        // horizontal propagation the height transition is a unit-mass atom.
        let inverse_rate = direction.z.abs() / area;
        let travel = exponential * inverse_rate;
        if direction.z > 0.0 && travel >= depth {
            let log_density = -depth / inverse_rate;
            if !log_density.is_finite() {
                return Err(SmithError::NumericalRange);
            }
            return Ok(SmithStep {
                event: SmithEvent::Escape,
                log_density,
            });
        }
        let next_depth = if direction.z > 0.0 {
            depth - travel
        } else {
            depth + travel
        };
        let height_log_density = if direction.z == 0.0 {
            0.0
        } else {
            if inverse_rate <= 0.0 || !inverse_rate.is_finite() {
                return Err(SmithError::NumericalRange);
            }
            -det::ln(inverse_rate) - exponential
        };
        if !next_depth.is_finite() || next_depth < 0.0 {
            return Err(SmithError::NumericalRange);
        }
        let facet = self.visible_normal(wo, draw[1], draw[2])?;
        let cosine = wo.dot(facet);
        if cosine <= 0.0 || facet.z <= 0.0 {
            return Err(SmithError::NumericalRange);
        }
        let reflected = unit(Vec3::new(
            2.0 * cosine * facet.x - wo.x,
            2.0 * cosine * facet.y - wo.y,
            2.0 * cosine * facet.z - wo.z,
        ))?;
        // VNDF = D(m) * dot(wo,m) / area; reflection contributes
        // 1 / (4 * dot(wo,m)). Cancel the dot product analytically.
        let log_density = height_log_density + det::ln(self.ndf(facet)) - det::ln(4.0 * area);
        if !log_density.is_finite() {
            return Err(SmithError::NumericalRange);
        }
        Ok(SmithStep {
            event: SmithEvent::Reflection {
                depth: next_depth,
                facet,
                direction: reflected,
            },
            log_density,
        })
    }

    pub(crate) fn projected_area(self, wo: Vec3) -> f64 {
        let x = self.alpha_x * wo.x;
        let y = self.alpha_y * wo.y;
        let radial_squared = x * x + y * y;
        let length = det::sqrt(radial_squared + wo.z * wo.z);
        if wo.z >= 0.0 {
            0.5 * (length + wo.z)
        } else {
            // Avoid cancellation for rays nearly escaping along the normal.
            0.5 * radial_squared / (length - wo.z)
        }
    }

    pub(crate) fn ndf(self, m: Vec3) -> f64 {
        let x = m.x / self.alpha_x;
        let y = m.y / self.alpha_y;
        let length_squared = x * x + y * y + m.z * m.z;
        1.0 / (PI * self.alpha_x * self.alpha_y * length_squared * length_squared)
    }

    // Dupuy & Benyoub (2023), sections 2--3, including below-horizon views:
    // https://arxiv.org/pdf/2306.05044
    fn visible_normal(self, wo: Vec3, u: f64, v: f64) -> Result<Vec3, SmithError> {
        let view = unit(Vec3::new(self.alpha_x * wo.x, self.alpha_y * wo.y, wo.z))?;
        let cap_height = if view.z >= 0.0 {
            1.0 + view.z
        } else {
            (view.x * view.x + view.y * view.y) / (1.0 - view.z)
        };
        let cap_drop = u * cap_height;
        let radius = det::sqrt((cap_drop * (2.0 - cap_drop)).max(0.0));
        let (sine, cosine) = det::sin_cos(2.0 * PI * v);
        // view.z + cap.z = (1 + view.z) * (1-u); this form retains a
        // positive facet even near the backfacing limit.
        unit(Vec3::new(
            self.alpha_x * (view.x + radius * cosine),
            self.alpha_y * (view.y + radius * sine),
            cap_height * (1.0 - u),
        ))
    }
}

pub(crate) fn unit(v: Vec3) -> Result<Vec3, SmithError> {
    if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) {
        return Err(SmithError::InvalidDirection);
    }
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if scale == 0.0 {
        return Err(SmithError::InvalidDirection);
    }
    let v = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
    Ok(v.scale(1.0 / det::sqrt(v.dot(v))))
}

#[cfg(test)]
mod tests;
