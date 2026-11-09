//! Exact positive absorption in a finite normal-depth interval.
//!
//! With `u = alpha*z/cos(theta)` and total optical depth `tau`, the
//! forward/backward intensities are `E exp(-u)` and
//! `E R_rear exp(-(2*tau-u))`, where
//! `E = T_front/(1-R_front*R_rear*exp(-2*tau))`. Integrating their positive
//! Beer losses gives the cell energy without subtracting boundary powers.
//! This is the incoherent intensity construction of Byrnes, *Multilayer
//! optical calculations*, §6, <https://arxiv.org/abs/1603.02720>.

use super::{DielectricError, IncoherentSlabResponse, det};

/// Original slab response and a prepared physical absorption profile.
///
/// The optical properties remain homogeneous and frozen. Depth cells partition
/// absorption, not the material: they introduce no new optical interfaces.
/// Each query sums all internal reflections exactly within the existing
/// incoherent real-index/Beer approximation. No allocation, stochastic sample,
/// spatial quadrature, or additional spectral material query is required.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedIncoherentSlabAbsorption {
    pub(super) optical: IncoherentSlabResponse,
    pub(super) thickness_m: f64,
    pub(super) alpha_per_m: f64,
    pub(super) entering_sum: [f64; 2],
    pub(super) rear_reflectance: [f64; 2],
}

impl PreparedIncoherentSlabAbsorption {
    /// Full original polarized reflection, transmission and absorption.
    #[must_use]
    pub const fn optical(&self) -> IncoherentSlabResponse {
        self.optical
    }

    /// Physical normal thickness, in metres.
    #[must_use]
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Positive `[s, p]` absorbed fractions of incident normal power in the
    /// interval `start_depth_m..end_depth_m` measured from the illuminated face.
    ///
    /// Depths must be finite and satisfy `0 <= start <= end <= thickness`.
    /// Empty intervals, transparent slabs and front-face total reflection
    /// absorb exactly zero. The full interval agrees with the original total
    /// absorptance; adjacent intervals add without introducing interface loss.
    ///
    /// For incident powers `[P_s, P_p]` in watts, this cell receives
    /// `A_s*P_s + A_p*P_p` watts. Unpolarized incidence splits total power equally
    /// between the modes. Power already transmitted through the front surface
    /// must not be used here: the original entrance transmission is included.
    /// Scattering, coherent standing waves and temperature-dependent index or
    /// absorption within the slab remain outside this prepared model.
    pub fn absorption_between(
        &self,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<[f64; 2], DielectricError> {
        self.absorption_interval(start_depth_m, end_depth_m, false)
    }

    /// Positive `[s, p]` fractions for an interval whose depths are measured
    /// from the opposite, unilluminated face. The prepared incident interface
    /// and ray remain unchanged; only the coordinate origin is reversed.
    ///
    /// This retains the original interval width even when subtracting both
    /// endpoints from a much larger thickness would round them to the same
    /// depth. Use this query when illumination reverses relative to an authored
    /// thermal partition. The same bounds and power convention as
    /// [`Self::absorption_between`] apply.
    pub fn absorption_between_from_back(
        &self,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<[f64; 2], DielectricError> {
        self.absorption_interval(start_depth_m, end_depth_m, true)
    }

    fn absorption_interval(
        &self,
        start_depth_m: f64,
        end_depth_m: f64,
        from_back: bool,
    ) -> Result<[f64; 2], DielectricError> {
        if !start_depth_m.is_finite()
            || !end_depth_m.is_finite()
            || start_depth_m < 0.0
            || end_depth_m < start_depth_m
            || end_depth_m > self.thickness_m
        {
            return Err(DielectricError::InvalidAbsorption);
        }
        if start_depth_m == end_depth_m
            || self.alpha_per_m == 0.0
            || self.optical.transmitted_cosine == 0.0
        {
            return Ok([0.0; 2]);
        }
        let interval_loss = -det::expm1(-self.depth(end_depth_m - start_depth_m));
        let (forward_distance, return_distance) = if from_back {
            (self.thickness_m - end_depth_m, start_depth_m)
        } else {
            (start_depth_m, self.thickness_m - end_depth_m)
        };
        let forward = det::exp(-self.depth(forward_distance));
        // Resolve the backward path as two positive distances. Computing
        // 2*tau-u_end would form infinity-infinity for an opaque slab and
        // would cancel away a very thin rear interval at finite large tau.
        let backward = det::exp(-(self.optical.optical_depth + self.depth(return_distance)));
        Ok([0, 1].map(|pol| {
            // Combine the bounded field weights before multiplying a
            // potentially subnormal loss. E*loss can round to zero even when
            // the forward-plus-return loss is still representable.
            interval_loss
                * (self.entering_sum[pol] * (forward + self.rear_reflectance[pol] * backward))
        }))
    }

    fn depth(&self, distance_m: f64) -> f64 {
        if distance_m == 0.0 {
            return 0.0;
        }
        let value = self.alpha_per_m * distance_m / self.optical.transmitted_cosine;
        if value > 0.0 && value.is_finite() {
            value
        } else {
            // Preserve a representable short-depth query even when the full
            // slab is opaque or alpha*distance loses range before division.
            det::exp(
                det::ln(self.alpha_per_m) + det::ln(distance_m)
                    - det::ln(self.optical.transmitted_cosine),
            )
        }
    }
}

#[cfg(test)]
mod tests;
