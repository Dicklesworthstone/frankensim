//! Axisymmetric mode-I/III Griffith equilibrium on the original JKR source.
//!
//! Full-stick torsion is elastically decoupled from axisymmetric normal
//! deformation: M=(16/3) Gbar a^3 phi, 1/Gbar=1/Ga+1/Gb. See Popov, Hess and
//! Willert (2019), DOI 10.1007/978-3-662-58709-6_5, Eqs. 5.2-5.3. The summed
//! mode-I/III release is also used by Giannakopoulos, Lindley and Suresh (1999),
//! DOI 10.1016/S1359-6454(99)00312-2, Eqs. 12-13. Here the EXISTING owner's
//! mode-independent, ideally brittle Griffith assumption explicitly gives
//! U=E*a[(delta-a^2/(3R))^2+4a^4/(45R^2)]
//!   +(8/3)Gbar*a^3*phi^2-pi*W*a^2,
//! E*(delta-a^2/R)^2+8Gbar*a^2*phi^2=2pi*W*a.
//!
//! Area, normal force and torque therefore relax together. The larger radial
//! local minimum is followed under prescribed (delta,phi). At negative delta
//! it ends at a genuine radial fold. At positive delta the opening-edge
//! condition a^2/R>=delta ends first: continuing the squared equation would
//! admit a closed/overlapping crack. Neither boundary invents a detached state,
//! frictional slip, heat or re-adhesion. At delta=0 there is no finite ideal
//! fold; original small-scale-yielding and authored twist limits still apply.
//! This reversible full-stick approximation is not local Coulomb partial slip
//! or experimentally fitted mixed-mode toughness.

mod radius;

use super::{JkrMixedModeError, SourceJkrMixedMode, finite, positive};
use crate::normal_patch::adhesive::JkrError;
use core::f64::consts::PI;
use fs_blake3::ContentHash;
use fs_tribo::InputAuthority;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrTorsionalLimitKind {
    RadialFold,
    /// Vanishing mode-I opening; not a radial instability or detachment event.
    OpeningEdge,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalLimit {
    pub kind: JkrTorsionalLimitKind,
    pub twist_angle_rad: f64,
    pub contact_radius_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrTorsionalError {
    Source(JkrMixedModeError),
    BeyondAttachedLimit {
        kind: JkrTorsionalLimitKind,
        twist_angle_rad: f64,
        limit_rad: f64,
    },
    NoConvergence,
}
impl From<JkrMixedModeError> for JkrTorsionalError {
    fn from(value: JkrMixedModeError) -> Self {
        Self::Source(value)
    }
}
impl From<JkrError> for JkrTorsionalError {
    fn from(value: JkrError) -> Self {
        Self::Source(value.into())
    }
}
impl From<radius::Unresolved> for JkrTorsionalError {
    fn from(_: radius::Unresolved) -> Self {
        Self::NoConvergence
    }
}
impl core::fmt::Display for JkrTorsionalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "axisymmetric I/III JKR equilibrium: {self:?}")
    }
}
impl std::error::Error for JkrTorsionalError {}

/// Diagnostics from the original immutable source, not a mutable state token.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalResponse {
    normal_identity: ContentHash,
    model_identity: ContentHash,
    pub boundary: Option<JkrTorsionalLimitKind>,
    pub indentation_m: f64,
    pub twist_angle_rad: f64,
    pub contact_radius_m: f64,
    pub contact_area_m2: f64,
    pub repulsive_normal_force_n: f64,
    pub normal_force_n: f64,
    /// M=dU/dphi enters storage. The torque resisting receiver rotation is -M.
    pub torsional_moment_n_m: f64,
    pub normal_elastic_energy_j: f64,
    pub torsional_elastic_energy_j: f64,
    pub surface_energy_j: f64,
    pub reversible_energy_j: f64,
    pub normal_energy_release_rate_j_m2: f64,
    pub torsional_energy_release_rate_j_m2: f64,
    pub griffith_residual_j_m2: f64,
    pub displacement_radius_curvature_n_per_m: f64,
    /// d(N,M)/d(delta,phi), symmetric after radius relaxation. Entry units are
    /// [[N/m,N/rad],[Nm/m,Nm/rad]]. None at the radial fold; at OpeningEdge
    /// this is the one-sided value on the admitted branch.
    pub coupled_stiffness: Option<[[f64; 2]; 2]>,
    pub authority: InputAuthority,
}
impl JkrTorsionalResponse {
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_identity
    }
    pub const fn model_identity(&self) -> ContentHash {
        self.model_identity
    }
}

impl SourceJkrMixedMode {
    /// Circular torsion compliance in series, without Mindlin's (2-nu) factors.
    pub fn torsional_reduced_shear_modulus_pa(&self) -> f64 {
        let elastic = self.normal.normal_state().elastic();
        let g = [elastic.surface_a_state(), elastic.surface_b_state()]
            .map(|s| s.elastic_state().shear_modulus_pa());
        g[0].min(g[1]) / (1.0 + g[0].min(g[1]) / g[0].max(g[1]))
    }

    /// Ideal attached-branch endpoint. Source/material envelopes can bind
    /// earlier. None means delta=0 has no finite radial or opening-edge limit.
    pub fn torsional_displacement_limit(
        &self,
        indentation_m: f64,
    ) -> Result<Option<JkrTorsionalLimit>, JkrTorsionalError> {
        let (_, boundary) = self.torsional_shape(indentation_m)?;
        boundary
            .map(|b| self.dimensional_torsional_limit(b))
            .transpose()
    }

    /// Prescribed twist of the receiver relative to its counterpart. Existing
    /// shear screens apply to rim displacement/a=|phi| and the mean absolute
    /// circumferential traction 2Gbar|phi|. They do not bound the singular edge.
    pub fn evaluate_torsion_displacement_controlled(
        &self,
        indentation_m: f64,
        twist_angle_rad: f64,
    ) -> Result<JkrTorsionalResponse, JkrTorsionalError> {
        finite(twist_angle_rad, "torsional angle")?;
        let g = self.torsional_reduced_shear_modulus_pa();
        let magnitude = twist_angle_rad.abs();
        for (ratio, value, limit) in [
            (
                "torsional rim displacement/radius",
                magnitude,
                self.envelope.maximum_tangential_displacement_to_radius,
            ),
            (
                "mean torsional shear/yield",
                2.0 * (g * magnitude)
                    / self
                        .normal
                        .normal_state()
                        .elastic()
                        .limiting_yield_stress_pa(),
                self.envelope.maximum_mean_shear_stress_to_yield,
            ),
        ] {
            if !value.is_finite() || value > limit {
                return Err(JkrMixedModeError::OutsideApplicability {
                    ratio,
                    value,
                    limit,
                }
                .into());
            }
        }
        let (d, shape) = self.torsional_shape(indentation_m)?;
        let mut boundary = None;
        let (y, offset) = if let Some(shape) = shape {
            let limit = self.dimensional_torsional_limit(shape)?;
            if magnitude > limit.twist_angle_rad {
                return Err(JkrTorsionalError::BeyondAttachedLimit {
                    kind: limit.kind,
                    twist_angle_rad,
                    limit_rad: limit.twist_angle_rad,
                });
            }
            if magnitude == limit.twist_angle_rad {
                boundary = Some(limit.kind);
                (shape.radius, 0.0)
            } else if magnitude == 0.0 {
                let y = self
                    .normal
                    .evaluate_attached(indentation_m)?
                    .contact_radius_m
                    / self.normal.spinodal_radius_m;
                (y, y - shape.radius)
            } else {
                let ratio = magnitude / limit.twist_angle_rad;
                let scaled = magnitude / self.torsional_twist_scale()?;
                radius::attached(
                    d,
                    scaled * scaled,
                    shape,
                    shape.twist_square * (1.0 - ratio) * (1.0 + ratio),
                )?
            }
        } else {
            let scaled = magnitude / self.torsional_twist_scale()?;
            (radius::zero_approach(scaled * scaled)?, 0.0)
        };
        let a = self.normal.spinodal_radius_m * y;
        self.normal.check_applicability(a)?;
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let hertz_depth = if boundary == Some(JkrTorsionalLimitKind::OpeningEdge) {
            indentation_m
        } else {
            a * (a / self.normal.radius_m)
        };
        let opening = hertz_depth - indentation_m;
        let mean = indentation_m - hertz_depth / 3.0;
        let area = PI * a * a;
        let stiffness = ((16.0 / 3.0 * a) * (g * a)) * a;
        let moment = stiffness * twist_angle_rad;
        let normal_energy =
            (e * a * mean) * mean + ((4.0 / 45.0) * e * a * hertz_depth) * hertz_depth;
        let torsion_energy = (0.5 * stiffness * twist_angle_rad) * twist_angle_rad;
        let surface_energy = -w * area;
        let normal_release = (e * opening) * (opening / a) / (2.0 * PI);
        let torsion_release = (4.0 / PI) * (g * twist_angle_rad) * (a * twist_angle_rad);
        let residual = normal_release + torsion_release - w;
        let curvature = if let Some(shape) = shape.filter(|s| s.radial_fold) {
            let f = shape.radius;
            (PI * w / 4.0) * (offset / y) * ((y + f) * (y * y + f * f) + 8.0)
        } else {
            let scaled_opening = if boundary.is_some() { 0.0 } else { y * y - d };
            (PI * w / 4.0) * (8.0 + scaled_opening * ((y * y + d) / y))
        };
        let diagonal = [2.0 * e * a, stiffness];
        let cross = [-2.0 * e * opening, (16.0 * (g * a) * a) * twist_angle_rad];
        let coupled_stiffness = if curvature > 0.0 {
            let off = -(cross[0] / curvature) * cross[1];
            Some([
                [diagonal[0] - (cross[0] / curvature) * cross[0], off],
                [off, diagonal[1] - (cross[1] / curvature) * cross[1]],
            ])
        } else {
            None
        };
        let normal_force = 2.0 * e * a * mean;
        let repulsive = (4.0 / 3.0) * e * a * hertz_depth;
        let energy = normal_energy + torsion_energy + surface_energy;
        for value in [
            a,
            area,
            normal_force,
            repulsive,
            moment,
            stiffness,
            normal_energy,
            torsion_energy,
            surface_energy,
            energy,
            normal_release,
            torsion_release,
            residual,
            curvature,
        ] {
            finite(value, "torsional equilibrium response")?;
        }
        for value in coupled_stiffness.into_iter().flatten().flatten() {
            finite(value, "torsional coupled stiffness")?;
        }
        if a <= 0.0
            || area <= 0.0
            || stiffness <= 0.0
            || opening < 0.0
            || curvature < 0.0
            || residual.abs() > 512.0 * f64::EPSILON * w
        {
            return Err(JkrTorsionalError::NoConvergence);
        }
        Ok(JkrTorsionalResponse {
            normal_identity: self.normal.normal_state().identity(),
            model_identity: self.identity(),
            boundary,
            indentation_m,
            twist_angle_rad,
            contact_radius_m: a,
            contact_area_m2: area,
            repulsive_normal_force_n: repulsive,
            normal_force_n: normal_force,
            torsional_moment_n_m: moment,
            normal_elastic_energy_j: normal_energy,
            torsional_elastic_energy_j: torsion_energy,
            surface_energy_j: surface_energy,
            reversible_energy_j: energy,
            normal_energy_release_rate_j_m2: normal_release,
            torsional_energy_release_rate_j_m2: torsion_release,
            griffith_residual_j_m2: residual,
            displacement_radius_curvature_n_per_m: curvature,
            coupled_stiffness,
            authority: self.normal.authority(),
        })
    }

    fn torsional_shape(
        &self,
        delta: f64,
    ) -> Result<(f64, Option<radius::Boundary>), JkrTorsionalError> {
        self.normal.validate_indentation(delta)?;
        let minimum = self.normal.minimum_attached_indentation_m();
        if delta < minimum {
            return Err(JkrError::NoAttachedEquilibrium {
                indentation_m: delta,
                minimum_indentation_m: minimum,
            }
            .into());
        }
        let d = delta / self.normal.displacement_scale_m;
        let distance = (delta - minimum) / self.normal.displacement_scale_m;
        Ok((d, radius::boundary(d, distance)?))
    }
    fn torsional_twist_scale(&self) -> Result<f64, JkrTorsionalError> {
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let scale = (self.normal.spinodal_radius_m / self.normal.radius_m)
            * (e / self.torsional_reduced_shear_modulus_pa() / 8.0).sqrt();
        positive(scale, "torsional twist scale")?;
        Ok(scale)
    }
    fn dimensional_torsional_limit(
        &self,
        shape: radius::Boundary,
    ) -> Result<JkrTorsionalLimit, JkrTorsionalError> {
        Ok(JkrTorsionalLimit {
            kind: if shape.radial_fold {
                JkrTorsionalLimitKind::RadialFold
            } else {
                JkrTorsionalLimitKind::OpeningEdge
            },
            twist_angle_rad: finite(
                self.torsional_twist_scale()? * shape.twist_square.sqrt(),
                "torsional boundary angle",
            )?,
            contact_radius_m: finite(
                self.normal.spinodal_radius_m * shape.radius,
                "torsional boundary radius",
            )?,
        })
    }
}
