//! Fixed-direction connections for the same height-resolved Smith walk.
//!
//! Heitz et al. (2016), sections 6--7: at each sampled collision, connect to
//! the requested upper-hemisphere direction with phase density times escape
//! probability. Dividing by the exit cosine converts that exit density into
//! a BRDF. These are integrands, NEVER outgoing marginal proposal PDFs.
//!
//! The first collision can be integrated analytically. With incoming cosine
//! ci and outgoing co, its height-integrated denominator is
//! 4 * [ci*co + A(-wi)*co + A(-wo)*ci], where A is projected microfacet area.
//! This is the reciprocal height-correlated Smith single-scattering law.

use super::{SmithError, SmithGgx, unit};
use fs_geom::Vec3;
use fs_math::det;

/// A deterministic connection at a given hidden collision, before applying
/// the caller's complex Fresnel/Mueller response at this facet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithConnection {
    pub facet: Vec3,
    /// Geometric BRDF factor per steradian, with escape already included.
    /// It is not a sampling density and contains no optical attenuation.
    pub factor_per_sr: f64,
}

/// Positive connection retained in logarithmic form when a single hidden
/// contribution is outside the linear range. The complete sum may still be
/// representable. This is a BRDF integrand, never a directional proposal PDF.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithLogConnection {
    pub facet: Vec3,
    pub log_factor_per_sr: f64,
}

impl SmithGgx {
    /// Logarithmic companion to the linear connection query. It preserves a
    /// grazing positive escape contribution until optical multiplication and
    /// summation, rather than turning its isolated underflow into extinction.
    pub fn log_connection_at_depth(
        self,
        incoming: Vec3,
        depth: f64,
        outgoing: Vec3,
    ) -> Result<Option<SmithLogConnection>, SmithError> {
        if !depth.is_finite() || depth < 0.0 {
            return Err(SmithError::InvalidDepth);
        }
        let incoming = unit(incoming)?;
        let outgoing = unit(outgoing)?;
        if outgoing.z <= 0.0 {
            return Ok(None);
        }
        let Some(facet) = connection_facet(incoming, outgoing)? else {
            return Ok(None);
        };
        let area = self.projected_area(incoming.scale(-1.0));
        if area <= 0.0 || !area.is_finite() {
            return Err(SmithError::NumericalRange);
        }
        let escape_area = self.projected_area(outgoing.scale(-1.0));
        let optical_depth = if escape_area == 0.0 || depth == 0.0 {
            0.0
        } else {
            depth * (escape_area / outgoing.z)
        };
        let log_factor = det::ln(self.ndf(facet))
            - det::ln(4.0)
            - det::ln(area)
            - det::ln(outgoing.z)
            - optical_depth;
        if !log_factor.is_finite() {
            return Err(SmithError::NumericalRange);
        }
        Ok(Some(SmithLogConnection {
            facet,
            log_factor_per_sr: log_factor,
        }))
    }
    /// Integrate the first collision's height exactly. Input is propagation
    /// direction into the surface; output points out toward the light.
    /// No visible connection (including a macro horizon) returns `None`.
    pub fn single_scattering_connection(
        self,
        incoming: Vec3,
        outgoing: Vec3,
    ) -> Result<Option<SmithConnection>, SmithError> {
        let incoming = unit(incoming)?;
        let outgoing = unit(outgoing)?;
        if incoming.z >= 0.0 || outgoing.z <= 0.0 {
            return Ok(None);
        }
        let Some(facet) = connection_facet(incoming, outgoing)? else {
            return Ok(None);
        };
        let ci = -incoming.z;
        let co = outgoing.z;
        let denominator = 4.0
            * (ci * co
                + self.projected_area(incoming) * co
                + self.projected_area(outgoing.scale(-1.0)) * ci);
        let factor = self.ndf(facet) / denominator;
        Ok(Some(SmithConnection {
            facet,
            factor_per_sr: positive(factor)?,
        }))
    }

    /// Connect a sampled collision at nonnegative CDF depth to a fixed light
    /// direction. The caller samples arrival at the depth exactly once, then
    /// uses this factor before sampling/applying the continuing facet.
    pub fn connection_at_depth(
        self,
        incoming: Vec3,
        depth: f64,
        outgoing: Vec3,
    ) -> Result<Option<SmithConnection>, SmithError> {
        if !depth.is_finite() || depth < 0.0 {
            return Err(SmithError::InvalidDepth);
        }
        let incoming = unit(incoming)?;
        let outgoing = unit(outgoing)?;
        if outgoing.z <= 0.0 {
            return Ok(None);
        }
        let Some(facet) = connection_facet(incoming, outgoing)? else {
            return Ok(None);
        };
        let area = self.projected_area(incoming.scale(-1.0));
        if area <= 0.0 || !area.is_finite() {
            return Err(SmithError::NumericalRange);
        }
        let escape_area = self.projected_area(outgoing.scale(-1.0));
        let optical_depth = if escape_area == 0.0 || depth == 0.0 {
            0.0
        } else {
            depth * (escape_area / outgoing.z)
        };
        let factor = self.ndf(facet) / (4.0 * area * outgoing.z) * det::exp(-optical_depth);
        let factor = if factor.is_finite() && factor > 0.0 {
            factor
        } else {
            // Preserve representable results when either the phase density or
            // escape probability alone over/underflows. An unrepresentable
            // positive connection is a refusal, not exact optical extinction.
            positive(det::exp(
                det::ln(self.ndf(facet))
                    - det::ln(4.0)
                    - det::ln(area)
                    - det::ln(outgoing.z)
                    - optical_depth,
            ))?
        };
        Ok(Some(SmithConnection {
            facet,
            factor_per_sr: factor,
        }))
    }
}

fn connection_facet(incoming: Vec3, outgoing: Vec3) -> Result<Option<Vec3>, SmithError> {
    let half = Vec3::new(
        outgoing.x - incoming.x,
        outgoing.y - incoming.y,
        outgoing.z - incoming.z,
    );
    if half.z <= 0.0 {
        return Ok(None);
    }
    let facet = unit(half)?;
    if -incoming.dot(facet) <= 0.0 || outgoing.dot(facet) <= 0.0 {
        return Ok(None);
    }
    Ok(Some(facet))
}

fn positive(value: f64) -> Result<f64, SmithError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(SmithError::NumericalRange)
    }
}

#[cfg(test)]
mod tests;
