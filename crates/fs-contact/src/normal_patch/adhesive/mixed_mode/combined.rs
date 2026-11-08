//! One relaxed circular contact under approach, vector shear AND twist.
//!
//! This extends the existing ideally brittle, mode-independent Griffith model:
//! U=E*a[(delta-a²/(3R))²+4a⁴/(45R²)]+4Gstar*a*|u|²
//!   +(8/3)Gbar*a³*phi²-pi*W*a².
//! Translation and twist retain their original full-stick compliances before
//! the SAME radius is equilibrated. All four conjugate efforts and every cross
//! derivative therefore respond to combined loading.
//!
//! The existing decoupled isotropic-halfspace and circular/perimeter-averaged
//! assumptions remain explicit. Tangential compliance: Popov, Hess and Willert
//! (2019), DOI 10.1007/978-3-662-58709-6_4, Eqs. 4.7-4.10 (their shear-modulus
//! normalization is four times Gstar here). No-slip torsion: chapter 5,
//! DOI 10.1007/978-3-662-58709-6_5, Eqs. 5.2-5.3. The original Savkoor-Briggs
//! Griffith choice and its limitations are documented by the parent module.
//! These are reversible full-stick branches, not Coulomb partial slip, fitted
//! mode-mixity toughness or a post-fold detachment/heat model.

mod radius;

/// Work-consistent four-coordinate exchange and whole-path source admission.
pub mod path;

use super::{
    JkrMixedModeError, SourceJkrMixedMode,
    torsion::{JkrTorsionalError, JkrTorsionalLimitKind},
};
use crate::normal_patch::adhesive::JkrError;
use core::f64::consts::PI;
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use fs_tribo::InputAuthority;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrCombinedBoundary {
    RadialFold,
    /// Vanishing normal opening; this is not automatically a detachment event.
    OpeningEdge,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrCombinedError {
    Source(JkrMixedModeError),
    Torsional(JkrTorsionalError),
    NoAttachedEquilibrium,
    /// The stationary branch would require a closed/reversed normal crack.
    OpeningEdge,
    /// Arithmetic cannot resolve the original stable opening branch.
    Unresolved,
}
impl From<JkrMixedModeError> for JkrCombinedError {
    fn from(error: JkrMixedModeError) -> Self {
        Self::Source(error)
    }
}
impl From<JkrTorsionalError> for JkrCombinedError {
    fn from(error: JkrTorsionalError) -> Self {
        Self::Torsional(error)
    }
}
impl From<JkrError> for JkrCombinedError {
    fn from(error: JkrError) -> Self {
        Self::Source(error.into())
    }
}
impl core::fmt::Display for JkrCombinedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "combined normal/shear/twist JKR refusal: {self:?}")
    }
}
impl std::error::Error for JkrCombinedError {}

/// Original attached equilibrium, including explicitly identified boundaries.
/// Public diagnostics cannot authenticate a caller-authored constitutive state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedResponse {
    normal_identity: ContentHash,
    model_identity: ContentHash,
    /// [normal approach m, shear1 m, shear2 m, twist rad].
    pub coordinates: [f64; 4],
    /// Work-conjugate [normal force N, shear1 N, shear2 N, twist moment Nm].
    /// These enter the potential; actual body reactions have the opposite sign.
    pub force: [f64; 4],
    pub boundary: Option<JkrCombinedBoundary>,
    pub contact_radius_m: f64,
    pub contact_area_m2: f64,
    pub repulsive_normal_force_n: f64,
    pub normal_elastic_energy_j: f64,
    pub tangential_elastic_energy_j: f64,
    pub torsional_elastic_energy_j: f64,
    pub surface_energy_j: f64,
    pub reversible_energy_j: f64,
    pub normal_energy_release_rate_j_m2: f64,
    pub shear_energy_release_rate_j_m2: f64,
    pub torsional_energy_release_rate_j_m2: f64,
    pub griffith_residual_j_m2: f64,
    /// Actual original U_a, retaining the equilibrium root's finite residual.
    pub radial_residual_n: f64,
    pub displacement_radius_curvature_n_per_m: f64,
    /// Symmetric relaxed d(force)/d(coordinates), with mixed SI units. Absent
    /// at the radial fold. Opening-edge applicability can end at finite tangent.
    pub coupled_stiffness: Option<[[f64; 4]; 4]>,
    pub authority: InputAuthority,
}
impl JkrCombinedResponse {
    pub const fn normal_state_identity(&self) -> ContentHash {
        self.normal_identity
    }
    pub const fn model_identity(&self) -> ContentHash {
        self.model_identity
    }
}

impl SourceJkrMixedMode {
    /// Original stable circular equilibrium with all four held coordinates.
    /// Radius stationarity is solved once; independent unchanged-area shear
    /// and torsion patches would not represent this potential.
    ///
    /// Joint rim displacement and mean-traction triangle bounds use the SAME
    /// authored envelope and original yield stress. They do not certify the
    /// singular edge traction. All normal geometry/material checks remain.
    /// Generic arithmetic is bounded and allocation-free. An unresolved fold,
    /// lost stable root or negative normal opening refuses without a successor.
    pub fn evaluate_combined_displacement_controlled(
        &self,
        q: [f64; 4],
    ) -> Result<JkrCombinedResponse, JkrCombinedError> {
        if q.iter().any(|x| !x.is_finite()) {
            return Err(JkrMixedModeError::Invalid {
                field: "combined contact coordinates",
            }
            .into());
        }
        // The original displacement-stable shear branch has y>=f and
        // f²-d=4/f>0, so this exact reduction never imports a closed crack.
        if q[3] == 0.0 {
            return self.combined_from_shear(q);
        }
        if q[1] == 0.0 && q[2] == 0.0 {
            return self.combined_from_torsion(q);
        }
        self.normal.validate_indentation(q[0])?;
        let minimum = self.normal.minimum_attached_indentation_m();
        if q[0] < minimum {
            return Err(JkrError::NoAttachedEquilibrium {
                indentation_m: q[0],
                minimum_indentation_m: minimum,
            }
            .into());
        }
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let gbar = self.torsional_reduced_shear_modulus_pa();
        // Normalize directly from original dimensional coefficients in Dd.
        // The existing rounded radius scale is an arbitrary conditioning scale;
        // it must not silently replace the original adhesion coefficient by 16.
        let radial_scale = dd(self.normal.spinodal_radius_m);
        let depth_scale = radial_scale * (radial_scale / dd(self.normal.radius_m));
        let stationary_scale = dd(e) * depth_scale * depth_scale;
        let shear_square = dd(q[1]) * dd(q[1]) + dd(q[2]) * dd(q[2]);
        let d = dd(q[0]) / depth_scale;
        let b = dd(8.0) * dd(gbar) * radial_scale * radial_scale * dd(q[3]) * dd(q[3])
            / stationary_scale;
        let t =
            dd(0.5) * dd(self.tangential_stiffness_per_radius) * shear_square / stationary_scale;
        let ell = dd(2.0) * dd(PI) * dd(w) * radial_scale / stationary_scale;
        let root = radius::solve(d, b, t, ell).map_err(|error| match error {
            radius::Error::NoEquilibrium => JkrCombinedError::NoAttachedEquilibrium,
            radius::Error::OpeningEdge => JkrCombinedError::OpeningEdge,
            radius::Error::Unresolved => JkrCombinedError::Unresolved,
        })?;
        let a = self.normal.spinodal_radius_m * root.radius;
        self.normal.check_applicability(a)?;
        let hertz = a * (a / self.normal.radius_m);
        let opening = hertz - q[0];
        if opening < 0.0 {
            return Err(JkrCombinedError::OpeningEdge);
        }
        let boundary = if opening == 0.0 {
            Some(JkrCombinedBoundary::OpeningEdge)
        } else if root.marginal {
            Some(JkrCombinedBoundary::RadialFold)
        } else {
            None
        };
        let magnitude = q[1].hypot(q[2]);
        self.admit_combined_nominal_bounds(magnitude, q[3].abs(), a)?;
        let k = self.tangential_stiffness_per_radius;
        let shear_stiffness = k * a;
        let twist_stiffness = self.combined_torsional_stiffness(a);
        let mean = q[0] - hertz / 3.0;
        let normal_energy = (e * a * mean) * mean + ((4.0 / 45.0) * e * a * hertz) * hertz;
        let shear_energy = (0.5 * shear_stiffness * magnitude) * magnitude;
        let twist_energy = (0.5 * twist_stiffness * q[3]) * q[3];
        let area = PI * a * a;
        let surface_energy = -w * area;
        let normal_release = (e * opening) * (opening / a) / (2.0 * PI);
        let shear_release = shear_energy / (2.0 * area);
        let twist_release = (4.0 / PI) * (gbar * q[3]) * (a * q[3]);
        let residual = normal_release + shear_release + twist_release - w;
        let curvature = ((stationary_scale / radial_scale) * dd(root.derivative)).to_f64();
        let radial_residual = self.combined_radial_residual(q, a);
        if root.marginal {
            let original_opening = dd(a) * (dd(a) / dd(self.normal.radius_m)) - dd(q[0]);
            let actual_curvature = (dd(4.0) * dd(e) * dd(a) / dd(self.normal.radius_m)
                * original_opening
                + dd(16.0) * dd(gbar) * dd(a) * dd(q[3]) * dd(q[3])
                - dd(2.0) * dd(PI) * dd(w))
            .to_f64();
            if radial_residual != 0.0 || actual_curvature != 0.0 {
                return Err(JkrCombinedError::Unresolved);
            }
        }
        let diagonal = [
            2.0 * e * a,
            shear_stiffness,
            shear_stiffness,
            twist_stiffness,
        ];
        let cross = [
            -2.0 * e * opening,
            k * q[1],
            k * q[2],
            (16.0 * (gbar * a) * a) * q[3],
        ];
        let coupled_stiffness = if curvature > 0.0 {
            let mut tangent = [[0.0; 4]; 4];
            for i in 0..4 {
                for j in i..4 {
                    let value = (if i == j { diagonal[i] } else { 0.0 })
                        - (cross[i] / curvature) * cross[j];
                    tangent[i][j] = value;
                    tangent[j][i] = value;
                }
            }
            Some(tangent)
        } else {
            None
        };
        self.checked_combined(JkrCombinedResponse {
            normal_identity: self.normal.normal_state().identity(),
            model_identity: self.identity(),
            coordinates: q,
            force: [
                2.0 * e * a * mean,
                shear_stiffness * q[1],
                shear_stiffness * q[2],
                twist_stiffness * q[3],
            ],
            boundary,
            contact_radius_m: a,
            contact_area_m2: area,
            repulsive_normal_force_n: (4.0 / 3.0) * e * a * hertz,
            normal_elastic_energy_j: normal_energy,
            tangential_elastic_energy_j: shear_energy,
            torsional_elastic_energy_j: twist_energy,
            surface_energy_j: surface_energy,
            reversible_energy_j: normal_energy + shear_energy + twist_energy + surface_energy,
            normal_energy_release_rate_j_m2: normal_release,
            shear_energy_release_rate_j_m2: shear_release,
            torsional_energy_release_rate_j_m2: twist_release,
            griffith_residual_j_m2: residual,
            radial_residual_n: radial_residual,
            displacement_radius_curvature_n_per_m: curvature,
            coupled_stiffness,
            authority: self.normal.authority(),
        })
    }

    fn combined_from_shear(&self, q: [f64; 4]) -> Result<JkrCombinedResponse, JkrCombinedError> {
        let original = self.evaluate_displacement_controlled(q[0], [q[1], q[2]])?;
        let a = original.contact_radius_m;
        let twist_stiffness = self.combined_torsional_stiffness(a);
        let tangent = original.coupled_stiffness_n_per_m.map(|old| {
            let mut out = [[0.0; 4]; 4];
            for i in 0..3 {
                for j in 0..3 {
                    out[i][j] = old[i][j];
                }
            }
            out[3][3] = twist_stiffness;
            out
        });
        self.checked_combined(JkrCombinedResponse {
            normal_identity: original.normal_state_identity(),
            model_identity: original.model_identity(),
            coordinates: q,
            force: [
                original.normal_force_n,
                original.tangential_force_n[0],
                original.tangential_force_n[1],
                twist_stiffness * q[3],
            ],
            boundary: original.marginal.then_some(JkrCombinedBoundary::RadialFold),
            contact_radius_m: a,
            contact_area_m2: original.contact_area_m2,
            repulsive_normal_force_n: original.repulsive_normal_force_n,
            normal_elastic_energy_j: original.normal_elastic_energy_j,
            tangential_elastic_energy_j: original.tangential_elastic_energy_j,
            torsional_elastic_energy_j: 0.0,
            surface_energy_j: original.surface_energy_j,
            reversible_energy_j: original.reversible_energy_j,
            normal_energy_release_rate_j_m2: original.normal_energy_release_rate_j_m2,
            shear_energy_release_rate_j_m2: original.shear_energy_release_rate_j_m2,
            torsional_energy_release_rate_j_m2: 0.0,
            griffith_residual_j_m2: original.griffith_residual_j_m2,
            radial_residual_n: self.combined_radial_residual(q, a),
            displacement_radius_curvature_n_per_m: original.displacement_radius_curvature_n_per_m,
            coupled_stiffness: tangent,
            authority: original.authority,
        })
    }

    fn combined_from_torsion(&self, q: [f64; 4]) -> Result<JkrCombinedResponse, JkrCombinedError> {
        let original = self.evaluate_torsion_displacement_controlled(q[0], q[3])?;
        let a = original.contact_radius_m;
        let shear_stiffness = self.tangential_stiffness_per_radius * a;
        let tangent = original.coupled_stiffness.map(|old| {
            let mut out = [[0.0; 4]; 4];
            for (i, row) in [0, 3].into_iter().enumerate() {
                for (j, column) in [0, 3].into_iter().enumerate() {
                    out[row][column] = old[i][j];
                }
            }
            out[1][1] = shear_stiffness;
            out[2][2] = shear_stiffness;
            out
        });
        self.checked_combined(JkrCombinedResponse {
            normal_identity: original.normal_state_identity(),
            model_identity: original.model_identity(),
            coordinates: q,
            force: [
                original.normal_force_n,
                shear_stiffness * q[1],
                shear_stiffness * q[2],
                original.torsional_moment_n_m,
            ],
            boundary: original.boundary.map(|boundary| match boundary {
                JkrTorsionalLimitKind::RadialFold => JkrCombinedBoundary::RadialFold,
                JkrTorsionalLimitKind::OpeningEdge => JkrCombinedBoundary::OpeningEdge,
            }),
            contact_radius_m: a,
            contact_area_m2: original.contact_area_m2,
            repulsive_normal_force_n: original.repulsive_normal_force_n,
            normal_elastic_energy_j: original.normal_elastic_energy_j,
            tangential_elastic_energy_j: 0.0,
            torsional_elastic_energy_j: original.torsional_elastic_energy_j,
            surface_energy_j: original.surface_energy_j,
            reversible_energy_j: original.reversible_energy_j,
            normal_energy_release_rate_j_m2: original.normal_energy_release_rate_j_m2,
            shear_energy_release_rate_j_m2: 0.0,
            torsional_energy_release_rate_j_m2: original.torsional_energy_release_rate_j_m2,
            griffith_residual_j_m2: original.griffith_residual_j_m2,
            radial_residual_n: self.combined_radial_residual(q, a),
            displacement_radius_curvature_n_per_m: original.displacement_radius_curvature_n_per_m,
            coupled_stiffness: tangent,
            authority: original.authority,
        })
    }

    fn combined_torsional_stiffness(&self, a: f64) -> f64 {
        ((16.0 / 3.0 * a) * (self.torsional_reduced_shear_modulus_pa() * a)) * a
    }

    fn combined_radial_residual(&self, q: [f64; 4], a: f64) -> f64 {
        let opening = dd(a) * (dd(a) / dd(self.normal.radius_m)) - dd(q[0]);
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let shear = dd(q[1]) * dd(q[1]) + dd(q[2]) * dd(q[2]);
        (dd(e) * opening * opening
            + dd(0.5) * dd(self.tangential_stiffness_per_radius) * shear
            + dd(8.0)
                * dd(self.torsional_reduced_shear_modulus_pa())
                * dd(a)
                * dd(a)
                * dd(q[3])
                * dd(q[3])
            - dd(2.0) * dd(PI) * dd(w) * dd(a))
        .to_f64()
    }

    fn admit_combined_nominal_bounds(
        &self,
        shear: f64,
        twist: f64,
        a: f64,
    ) -> Result<(), JkrCombinedError> {
        let yield_stress = self
            .normal
            .normal_state()
            .elastic()
            .limiting_yield_stress_pa();
        for (ratio, value, limit) in [
            (
                "combined rim displacement/radius",
                shear / a + twist,
                self.envelope.maximum_tangential_displacement_to_radius,
            ),
            (
                "combined mean shear/yield",
                (self.tangential_stiffness_per_radius * (shear / (PI * a))
                    + 2.0 * self.torsional_reduced_shear_modulus_pa() * twist)
                    / yield_stress,
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
        Ok(())
    }

    fn checked_combined(
        &self,
        response: JkrCombinedResponse,
    ) -> Result<JkrCombinedResponse, JkrCombinedError> {
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        if response.force.iter().any(|x| !x.is_finite())
            || [
                response.contact_radius_m,
                response.contact_area_m2,
                response.repulsive_normal_force_n,
                response.normal_elastic_energy_j,
                response.tangential_elastic_energy_j,
                response.torsional_elastic_energy_j,
                response.surface_energy_j,
                response.reversible_energy_j,
                response.normal_energy_release_rate_j_m2,
                response.shear_energy_release_rate_j_m2,
                response.torsional_energy_release_rate_j_m2,
                response.griffith_residual_j_m2,
                response.radial_residual_n,
                response.displacement_radius_curvature_n_per_m,
            ]
            .iter()
            .any(|x| !x.is_finite())
            || response
                .coupled_stiffness
                .iter()
                .flatten()
                .flatten()
                .any(|x| !x.is_finite())
            || response.contact_radius_m <= 0.0
            || response.contact_area_m2 <= 0.0
            || response.displacement_radius_curvature_n_per_m < 0.0
            || response.griffith_residual_j_m2.abs() > 512.0 * f64::EPSILON * w
        {
            return Err(JkrCombinedError::Unresolved);
        }
        Ok(response)
    }
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
