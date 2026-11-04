//! One original generated-conduction layer embedded in world coordinates.
//!
//! The conduction owner supplies T, dT/dx and d²T/dx². This adapter only rotates
//! its derivatives: grad(T)=T' a and Hess(T)=T'' a a^T for the declared unit axis
//! a. The left-face origin belongs to the selected layer, not the whole stack.
//! Each contact/material jump needs a separate field and native boundary event.
//! Planar symmetry is explicit in transverse directions; the wall's scalar heat
//! transfer area does not imply a lateral shape or aperture.

use super::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureJet,
    TemperatureSecondOrderField,
};
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use fs_conduction::ConductionError;
use fs_conduction::layered::{GeneratedLayerTemperatureProfile, GeneratedLayeredWallSolution};
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Immutable borrow of one original solved conduction layer. No additional
/// thermal solve, averaged temperature, coefficient fit or state copy is used.
pub struct PlanarConductionTemperatureField<'a> {
    solution: &'a GeneratedLayeredWallSolution,
    profile: GeneratedLayerTemperatureProfile<'a>,
    layer_index: usize,
    left_face_origin_m: Point3,
    axis: Vec3,
    thickness_m: f64,
}

impl<'a> PlanarConductionTemperatureField<'a> {
    /// Embed a layer along a finite unit direction from its authored left face.
    /// Unit-length roundoff is normalized once; arbitrary scaling is refused.
    /// The layer's closed axial interval is retained even if its temperature
    /// happens to be constant, so uniform optical shortcuts cannot bypass it.
    pub fn try_new(
        solution: &'a GeneratedLayeredWallSolution,
        layer_index: usize,
        left_face_origin_m: Point3,
        axis: Vec3,
    ) -> Result<Self, ConductionTemperatureFieldError> {
        if ![
            left_face_origin_m.x,
            left_face_origin_m.y,
            left_face_origin_m.z,
            axis.x,
            axis.y,
            axis.z,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(ConductionTemperatureFieldError::InvalidPlacement);
        }
        let square = axis.dot(axis);
        if !square.is_finite() || (square - 1.0).abs() > 2.0e-10 {
            return Err(ConductionTemperatureFieldError::InvalidPlacement);
        }
        let profile = solution.temperature_profile(layer_index)?;
        profile.temperature_jet_at(0.0)?;
        let thickness_m = profile.thickness_m();
        profile.temperature_jet_at(thickness_m)?;
        let axis = axis.scale(1.0 / det::sqrt(square));
        Ok(Self {
            solution,
            profile,
            layer_index,
            left_face_origin_m,
            axis,
            thickness_m,
        })
    }

    pub const fn solution(&self) -> &GeneratedLayeredWallSolution {
        self.solution
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn left_face_origin_m(&self) -> Point3 {
        self.left_face_origin_m
    }
    pub const fn axis(&self) -> Vec3 {
        self.axis
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Original conduction refusals remain available here. The native optical
    /// traits map a spatial refusal to `OutsideFieldDomain` and invalid source
    /// state to `InvalidSample`, without returning a partial temperature jet.
    pub fn sample_temperature_jet(
        &self,
        position_m: Point3,
    ) -> Result<TemperatureJet, ConductionTemperatureFieldError> {
        let displacement = Vec3::new(
            position_m.x - self.left_face_origin_m.x,
            position_m.y - self.left_face_origin_m.y,
            position_m.z - self.left_face_origin_m.z,
        );
        if ![displacement.x, displacement.y, displacement.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(ConductionTemperatureFieldError::OutsideLayer);
        }
        let offset_m = self.axis.dot(displacement);
        if !offset_m.is_finite() || !(0.0..=self.thickness_m).contains(&offset_m) {
            return Err(ConductionTemperatureFieldError::OutsideLayer);
        }
        let jet = self.profile.temperature_jet_at(offset_m)?;
        let a = [self.axis.x, self.axis.y, self.axis.z];
        Ok(TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: jet.temperature_k,
                gradient_k_per_m: self.axis.scale(jet.gradient_k_per_m),
            },
            diagonal_k_per_m2: a.map(|component| component * component * jet.curvature_k_per_m2),
            off_diagonal_k_per_m2: [
                a[0] * a[1] * jet.curvature_k_per_m2,
                a[0] * a[2] * jet.curvature_k_per_m2,
                a[1] * a[2] * jet.curvature_k_per_m2,
            ],
        })
    }
}

impl TemperatureField for PlanarConductionTemperatureField<'_> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map(|jet| jet.sample.temperature_k)
            .map_err(ConductionTemperatureFieldError::into_gradient_error)
    }
}
impl TemperatureGradientField for PlanarConductionTemperatureField<'_> {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map(|jet| jet.sample)
            .map_err(ConductionTemperatureFieldError::into_gradient_error)
    }
}
impl TemperatureSecondOrderField for PlanarConductionTemperatureField<'_> {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map_err(ConductionTemperatureFieldError::into_gradient_error)
    }
}

/// Original source or spatial-placement refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ConductionTemperatureFieldError {
    InvalidPlacement,
    OutsideLayer,
    Conduction(ConductionError),
}
impl ConductionTemperatureFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::OutsideLayer => GradientIndexError::OutsideFieldDomain,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl From<ConductionError> for ConductionTemperatureFieldError {
    fn from(error: ConductionError) -> Self {
        Self::Conduction(error)
    }
}
impl core::fmt::Display for ConductionTemperatureFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar conduction temperature field: {self:?}")
    }
}
impl core::error::Error for ConductionTemperatureFieldError {}

#[cfg(test)]
mod tests;
