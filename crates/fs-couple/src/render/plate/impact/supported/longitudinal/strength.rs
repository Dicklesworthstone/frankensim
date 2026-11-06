//! Tensile demand on the original retained longitudinal strain polynomial.
//!
//! The native owner supplies `T_max = T0 + EA (epsilon_mean + sum |epsilon_j|)`.
//! Combining it with the relative two-plane curvature envelope gives each
//! perfectly bonded layer's bound `E_i (T_max/EA + r_i K)`. Affine endpoint
//! motion enters the exact original mean/harmonics; it adds no curvature.
//! Unretained axial channels remain statically equilibrated as in the native
//! model, not silently populated with fictitious wave modes.

use super::*;
use crate::render::plate::impact::string::StringTensileAssessment;
use fs_material::section::strength::AxialTensileCapacity;

impl LongitudinalPotential {
    pub(in crate::render::plate::impact::supported) fn with_tensile_capacity(
        mut self,
        capacity: AxialTensileCapacity,
        length_m: f64,
    ) -> Result<Self, ImpactError> {
        if self.tensile.is_some() || !length_m.is_finite() || length_m <= 0.0 {
            return Err(invalid(
                "longitudinal tensile admission needs its original span and an exclusive limit",
            ));
        }
        // These are the actual transverse inertia/slope factors. Axial modes
        // use solid inertia and may be sparse; neither may replace this basis.
        let curvature_scales = self
            .scales
            .iter()
            .enumerate()
            .map(|(i, scale)| (i + 1) as f64 * core::f64::consts::PI / length_m * scale)
            .collect();
        self.tensile = Some(Arc::new(StringTensileGuard::new(
            capacity,
            self.ea,
            curvature_scales,
        )?));
        Ok(self)
    }

    pub(in crate::render::plate::impact::supported) fn tensile_assessment(
        &self,
        q: &[f64],
    ) -> Result<Option<StringTensileAssessment>, ImpactError> {
        let (_, axial) = self.observe_kinematics(q)?;
        self.tensile
            .as_ref()
            .map(|guard| guard.assess(q, 2, axial.maximum_tension_n))
            .transpose()
    }
}
