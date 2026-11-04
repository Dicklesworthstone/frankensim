//! One original temperature-dependent conduction layer in world coordinates.
//!
//! The original conduction owner inverts its retained conductivity segments and
//! supplies T, T' and (when unique) T''. This adapter rotates those derivatives
//! into the native optical field; it performs no thermal solve, spatial fit or
//! smoothing across conductivity knots, contacts or different material layers.
//! The left-face origin and unit axis place one physical layer. Planar symmetry
//! is explicit; a wall's scalar transfer area does not imply a lateral aperture.

use super::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureJet,
    TemperatureSecondOrderField,
};
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use fs_conduction::ConductionError;
use fs_conduction::layered::{
    NonlinearGeneratedLayerTemperatureProfile, NonlinearGeneratedLayeredWallSolution,
    NonlinearTemperatureProfileError,
};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Immutable original solved state and cancellation context. Each observation
/// reuses the source owner's bounded segment inversion without allocating.
/// First-order optical queries do not require second-order source smoothness.
pub struct PlanarNonlinearConductionTemperatureField<'a, 'cx> {
    profile: NonlinearGeneratedLayerTemperatureProfile<'a>,
    cx: &'a Cx<'cx>,
    left_face_origin_m: Point3,
    axis: Vec3,
}

impl<'a, 'cx> PlanarNonlinearConductionTemperatureField<'a, 'cx> {
    /// Embed one original physical layer along a finite unit axis. Unit-length
    /// roundoff is normalized once; arbitrary scaling is refused. The closed
    /// layer domain and source temperature support apply to every observation.
    pub fn try_new(
        cx: &'a Cx<'cx>,
        solution: &'a NonlinearGeneratedLayeredWallSolution,
        layer_index: usize,
        left_face_origin_m: Point3,
        axis: Vec3,
    ) -> Result<Self, NonlinearConductionTemperatureFieldError> {
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
            return Err(NonlinearConductionTemperatureFieldError::InvalidPlacement);
        }
        let square = axis.dot(axis);
        if !square.is_finite() || (square - 1.0).abs() > 2.0e-10 {
            return Err(NonlinearConductionTemperatureFieldError::InvalidPlacement);
        }
        let profile = solution.temperature_profile(cx, layer_index)?;
        profile.temperature_gradient_at(cx, 0.0)?;
        profile.temperature_gradient_at(cx, profile.thickness_m())?;
        Ok(Self {
            profile,
            cx,
            left_face_origin_m,
            axis: axis.scale(1.0 / det::sqrt(square)),
        })
    }

    pub const fn solution(&self) -> &NonlinearGeneratedLayeredWallSolution {
        self.profile.solution()
    }
    pub const fn layer_index(&self) -> usize {
        self.profile.layer_index()
    }
    pub const fn thickness_m(&self) -> f64 {
        self.profile.thickness_m()
    }
    pub const fn left_face_origin_m(&self) -> Point3 {
        self.left_face_origin_m
    }
    pub const fn axis(&self) -> Vec3 {
        self.axis
    }

    fn offset_m(&self, p: Point3) -> Result<f64, NonlinearConductionTemperatureFieldError> {
        let displacement = Vec3::new(
            p.x - self.left_face_origin_m.x,
            p.y - self.left_face_origin_m.y,
            p.z - self.left_face_origin_m.z,
        );
        if ![displacement.x, displacement.y, displacement.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(NonlinearConductionTemperatureFieldError::OutsideLayer);
        }
        let offset = self.axis.dot(displacement);
        if !offset.is_finite() || !(0.0..=self.thickness_m()).contains(&offset) {
            return Err(NonlinearConductionTemperatureFieldError::OutsideLayer);
        }
        Ok(offset)
    }

    /// T and its original Fourier gradient, including at continuous source
    /// knots. No curvature query is made by this first-order path.
    pub fn sample_temperature_gradient(
        &self,
        position_m: Point3,
    ) -> Result<TemperatureGradientSample, NonlinearConductionTemperatureFieldError> {
        let sample = self
            .profile
            .temperature_gradient_at(self.cx, self.offset_m(position_m)?)?;
        Ok(TemperatureGradientSample {
            temperature_k: sample.temperature_k,
            gradient_k_per_m: self.axis.scale(sample.gradient_k_per_m),
        })
    }

    /// Analytic Hessian T'' a a^T, available only when the original source
    /// admits that derivative. A source curvature jump is never smoothed.
    pub fn sample_temperature_jet(
        &self,
        position_m: Point3,
    ) -> Result<TemperatureJet, NonlinearConductionTemperatureFieldError> {
        let jet = self
            .profile
            .temperature_jet_at(self.cx, self.offset_m(position_m)?)?;
        let a = [self.axis.x, self.axis.y, self.axis.z];
        Ok(TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: jet.sample.temperature_k,
                gradient_k_per_m: self.axis.scale(jet.sample.gradient_k_per_m),
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

impl TemperatureField for PlanarNonlinearConductionTemperatureField<'_, '_> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_temperature_gradient(p)
            .map(|sample| sample.temperature_k)
            .map_err(NonlinearConductionTemperatureFieldError::into_gradient_error)
    }
}
impl TemperatureGradientField for PlanarNonlinearConductionTemperatureField<'_, '_> {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        self.sample_temperature_gradient(p)
            .map_err(NonlinearConductionTemperatureFieldError::into_gradient_error)
    }
}
impl TemperatureSecondOrderField for PlanarNonlinearConductionTemperatureField<'_, '_> {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map_err(NonlinearConductionTemperatureFieldError::into_gradient_error)
    }
}

/// Original source refusal, missing unique curvature or invalid placement.
#[derive(Clone, Debug, PartialEq)]
pub enum NonlinearConductionTemperatureFieldError {
    InvalidPlacement,
    OutsideLayer,
    Profile(NonlinearTemperatureProfileError),
}
impl NonlinearConductionTemperatureFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::OutsideLayer
            | Self::Profile(NonlinearTemperatureProfileError::Conduction(
                ConductionError::OutsideTemperatureSpan { .. },
            )) => GradientIndexError::OutsideFieldDomain,
            Self::Profile(NonlinearTemperatureProfileError::Conduction(
                ConductionError::Cancelled { .. },
            )) => GradientIndexError::Cancelled,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl From<ConductionError> for NonlinearConductionTemperatureFieldError {
    fn from(error: ConductionError) -> Self {
        Self::Profile(error.into())
    }
}
impl From<NonlinearTemperatureProfileError> for NonlinearConductionTemperatureFieldError {
    fn from(error: NonlinearTemperatureProfileError) -> Self {
        Self::Profile(error)
    }
}
impl core::fmt::Display for NonlinearConductionTemperatureFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar nonlinear conduction temperature field: {self:?}")
    }
}
impl core::error::Error for NonlinearConductionTemperatureFieldError {}

#[cfg(test)]
mod tests;
