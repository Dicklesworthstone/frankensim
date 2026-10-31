//! Analytic first variations across one smooth dielectric boundary.
//!
//! The central ray must already end at the admitted interface. The arrival
//! derivative includes the neighboring ray's different arc length, then the
//! normal Jacobian and BOTH spatial phase-index gradients differentiate Snell's
//! law. Reflection uses the same surface-arrival derivatives. Neither a new
//! intersection nor a finite-difference normal is inferred here.
//!
//! With incident direction `t`, face-forward normal `N`, `c=-N.t`, and
//! `eta=n_i/n_t`, transmission is `eta*t + (eta*c-c_t)*N`.
//! `dc_t=(-eta*(1-c*c)*deta + eta*eta*c*dc)/c_t`, including
//! `deta=eta*(grad(n_i).X/n_i-grad(n_t).X/n_t)`. These are local
//! differentials of a fixed smooth branch; grazing and critical transmission
//! require explicit refusal. A caustic inside a path does not itself prevent
//! differentiating the next interface; receiving-area admission remains a
//! separate operation.
//!
//! References: Igehy, *Tracing Ray Differentials* (1999),
//! <https://doi.org/10.1145/311535.311555>; PBRT 4e, section 10.1.3,
//! <https://www.pbr-book.org/4ed/Textures_and_Materials/Texture_Sampling_and_Antialiasing>.
//! The spatial index-gradient terms above extend the homogeneous-interface law.

use super::*;
use crate::dielectric::{DielectricEvent, fresnel_dielectric, reflect, refract};
use fs_geom::NormalDifferential;

/// Explicit local angular resolution at a smooth interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InterfaceOptions {
    pub min_incidence_cosine: f64,
    pub min_transmitted_cosine: f64,
}

impl Default for InterfaceOptions {
    fn default() -> Self {
        Self {
            min_incidence_cosine: 1e-6,
            min_transmitted_cosine: 1e-6,
        }
    }
}

impl InterfaceOptions {
    fn validate(self) -> Result<(), BundleInterfaceError> {
        if [self.min_incidence_cosine, self.min_transmitted_cosine]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0 || x > 1.0)
        {
            return Err(BundleInterfaceError::InvalidOptions);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BundleInterfaceError {
    Bundle(BundleError),
    Dielectric(DielectricError),
    InvalidOptions,
    InvalidNormalDifferential,
    GrazingIncidence,
    CriticalTransmission,
    TotalInternalReflection,
    UnavailableBranch,
    NumericalRange,
}

impl From<BundleError> for BundleInterfaceError {
    fn from(error: BundleError) -> Self {
        Self::Bundle(error)
    }
}
impl From<DielectricError> for BundleInterfaceError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}
impl core::fmt::Display for BundleInterfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "differentiated dielectric interface: {self:?}")
    }
}
impl core::error::Error for BundleInterfaceError {}

impl RayBundlePropagation {
    /// Continue a fixed reflection/transmission branch at the existing endpoint.
    /// `surface` is the analytic derivative of the geometric unit normal there;
    /// `target` is the actual optical sample across the boundary at that same
    /// position and vacuum wavelength. The caller owns contact localization,
    /// visibility, medium identity, and support of neighboring rays.
    ///
    /// Returned geometry retains the original launch solid-angle coordinates
    /// and past optical integrals. It carries NO Fresnel weight; use the native
    /// polarized bundle consumer for physical multi-interface irradiance.
    pub fn smooth_interface(
        &self,
        surface: NormalDifferential,
        target: GradientIndexSample,
        event: DielectricEvent,
        options: InterfaceOptions,
    ) -> Result<BundleContinuation, BundleInterfaceError> {
        options.validate()?;
        let (mut normal, mut jacobian) = admit_normal(surface)?;
        let tangent = self.propagation.ray.dir;
        if normal.dot(tangent) > 0.0 {
            normal = normal.scale(-1.0);
            jacobian = jacobian.map(|row| row.map(|v| -v));
        }
        let cosine = (-normal.dot(tangent)).clamp(0.0, 1.0);
        if cosine < options.min_incidence_cosine {
            return Err(BundleInterfaceError::GrazingIncidence);
        }
        let incident = self.propagation.final_sample;
        let eta_i = incident.phase_index();
        let eta_t = target.phase_index();
        let eta = eta_i / eta_t;
        let fresnel = fresnel_dielectric(cosine, eta_i, eta_t)?;
        let transmitted_cosine = fresnel.transmitted_cosine;
        if event == DielectricEvent::Transmission {
            if fresnel.total_internal_reflection {
                return Err(BundleInterfaceError::TotalInternalReflection);
            }
            if transmitted_cosine < options.min_transmitted_cosine {
                return Err(BundleInterfaceError::CriticalTransmission);
            }
        }

        // Central directions stay owned by the existing dielectric geometry.
        let wo = tangent.scale(-1.0);
        let raw_direction = match event {
            DielectricEvent::Reflection => reflect(wo, normal),
            DielectricEvent::Transmission => sub(
                tangent.scale(eta),
                normal.scale(transmitted_cosine - eta * cosine),
            ),
        };
        let raw_norm = raw_direction.norm();
        let direction = match event {
            DielectricEvent::Reflection => {
                unit(raw_direction).ok_or(BundleInterfaceError::NumericalRange)?
            }
            DielectricEvent::Transmission => refract(wo, normal, eta_i, eta_t, transmitted_cosine)?,
        };
        if !raw_norm.is_finite() || raw_norm <= 0.0 {
            return Err(BundleInterfaceError::NumericalRange);
        }
        let g = incident.phase_gradient_per_m().scale(1.0 / eta_i);
        let curvature = sub(g, tangent.scale(tangent.dot(g)));
        let mut columns = self.derivatives;
        for column in &mut columns {
            // Differentiating F(r(s(alpha),alpha))=0 needs only its normal.
            let ds = normal.dot(column.position_m_per_rad) / cosine;
            let x = sub(column.position_m_per_rad, tangent.scale(-ds));
            let v = sub(column.direction_per_rad, curvature.scale(-ds));
            let dn = apply_matrix(jacobian, x);
            let dc = -normal.dot(v) - tangent.dot(dn);
            let raw_v = match event {
                DielectricEvent::Reflection => {
                    sub(sub(v, normal.scale(-2.0 * dc)), dn.scale(-2.0 * cosine))
                }
                DielectricEvent::Transmission => {
                    let deta = eta
                        * (incident.phase_gradient_per_m().dot(x) / eta_i
                            - target.phase_gradient_per_m().dot(x) / eta_t);
                    let dct = (-eta * (1.0 - cosine * cosine) * deta + eta * eta * cosine * dc)
                        / transmitted_cosine;
                    let normal_factor = deta * cosine + eta * dc - dct;
                    sub(
                        sub(
                            sub(v.scale(eta), tangent.scale(-deta)),
                            normal.scale(-normal_factor),
                        ),
                        dn.scale(transmitted_cosine - eta * cosine),
                    )
                }
            };
            let v = sub(raw_v, direction.scale(direction.dot(raw_v))).scale(1.0 / raw_norm);
            if !ds.is_finite() || !finite_vec(x) || !finite_vec(v) {
                return Err(BundleInterfaceError::NumericalRange);
            }
            *column = AngularDerivative {
                position_m_per_rad: x,
                direction_per_rad: v,
            };
        }
        let s = cross(tangent, normal);
        let axis = if s == Vec3::new(0.0, 0.0, 0.0) {
            self.propagation.transverse_axis
        } else {
            unit(s).ok_or(BundleInterfaceError::NumericalRange)?
        };
        let axis = transverse(direction, axis).ok_or(BundleInterfaceError::NumericalRange)?;
        let outgoing_sample = if event == DielectricEvent::Transmission {
            target
        } else {
            incident
        };
        Ok(BundleContinuation::from_interface(
            self,
            Ray {
                origin: self.propagation.ray.origin,
                dir: direction,
            },
            axis,
            outgoing_sample,
            columns,
        )?)
    }
}

fn apply_matrix(matrix: [[f64; 3]; 3], value: Vec3) -> Vec3 {
    let v = [value.x, value.y, value.z];
    let product = matrix.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
    Vec3::new(product[0], product[1], product[2])
}

fn admit_normal(
    surface: NormalDifferential,
) -> Result<(Vec3, [[f64; 3]; 3]), BundleInterfaceError> {
    if !finite_vec(surface.normal)
        || surface
            .derivative_per_m
            .iter()
            .flatten()
            .any(|x| !x.is_finite())
    {
        return Err(BundleInterfaceError::InvalidNormalDifferential);
    }
    let length = vec_norm(surface.normal);
    if !length.is_finite() || (length - 1.0).abs() > 2e-12 {
        return Err(BundleInterfaceError::InvalidNormalDifferential);
    }
    let normal = surface.normal.scale(1.0 / length);
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        let column = Vec3::new(
            surface.derivative_per_m[0][i],
            surface.derivative_per_m[1][i],
            surface.derivative_per_m[2][i],
        );
        let norm = vec_norm(column);
        let parallel = normal.dot(column);
        if !norm.is_finite()
            || !parallel.is_finite()
            || parallel.abs() > 256.0 * f64::EPSILON * norm
        {
            return Err(BundleInterfaceError::InvalidNormalDifferential);
        }
        let column = sub(column, normal.scale(parallel)).scale(1.0 / length);
        result[0][i] = column.x;
        result[1][i] = column.y;
        result[2][i] = column.z;
    }
    Ok((normal, result))
}

#[cfg(test)]
mod tests;
