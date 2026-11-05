//! Graded original line mass in the existing Euler--Bernoulli inertia owner.

use super::{Frame2d, SolidError, invalid};
use crate::frame2d::SectionModel;

/// P1 original line mass and centroidal section spin inertia along a member.
///
/// Values follow the element's i-to-j reference coordinate. Line mass is
/// `integral rho dA` [kg/m]; spin inertia is `integral rho*y² dA` [kg m]. The
/// cross-section mass centroid lies on the declared centerline. These are
/// original inventories, not current density sampled on a deformed chord.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameInertiaProfile {
    line_mass_kg_m: [f64; 2],
    section_spin_kg_m: [f64; 2],
}

impl FrameInertiaProfile {
    /// Admit finite nonnegative endpoint moments. Section spin requires line
    /// mass at the same endpoint; zero moments retain the massless-member lane.
    pub fn new(line_mass_kg_m: [f64; 2], section_spin_kg_m: [f64; 2]) -> Result<Self, SolidError> {
        if (0..2).any(|i| {
            !line_mass_kg_m[i].is_finite()
                || line_mass_kg_m[i] < 0.0
                || !section_spin_kg_m[i].is_finite()
                || section_spin_kg_m[i] < 0.0
                || (line_mass_kg_m[i] == 0.0 && section_spin_kg_m[i] != 0.0)
        }) {
            return Err(invalid(
                "finite nonnegative endpoint line mass and section spin required; spin requires mass",
            ));
        }
        Ok(Self {
            line_mass_kg_m,
            section_spin_kg_m,
        })
    }

    /// Original endpoint line masses [kg/m], in member i-to-j order.
    #[must_use]
    pub const fn line_mass_kg_m(self) -> [f64; 2] {
        self.line_mass_kg_m
    }

    /// Original endpoint centroidal spin inertias per length [kg m].
    #[must_use]
    pub const fn section_spin_kg_m(self) -> [f64; 2] {
        self.section_spin_kg_m
    }

    pub(super) fn has_mass(self) -> bool {
        self.line_mass_kg_m.iter().any(|&v| v > 0.0)
    }

    fn endpoint_inventory(self, length: f64) -> [[f64; 2]; 2] {
        let half_product = |a: f64| {
            let product = a * length;
            if product.is_finite() {
                0.5 * product
            } else {
                (0.5 * a) * length
            }
        };
        core::array::from_fn(|i| {
            [
                half_product(self.line_mass_kg_m[i]),
                half_product(self.section_spin_kg_m[i]),
            ]
        })
    }

    /// Exact positive integral of each P1 shape function times each P1 moment.
    pub(in crate::frame2d) fn nodal_inventory(
        self,
        length: f64,
    ) -> Result<[[f64; 2]; 2], SolidError> {
        let end = self.endpoint_inventory(length);
        for k in 0..2 {
            let values = if k == 0 {
                self.line_mass_kg_m
            } else {
                self.section_spin_kg_m
            };
            if (0..2).any(|i| !end[i][k].is_finite() || (values[i] > 0.0 && end[i][k] <= 0.0))
                || !(end[0][k] + end[1][k]).is_finite()
            {
                return Err(invalid(
                    "graded reference member inventory must be representable",
                ));
            }
        }
        let nodal: [[f64; 2]; 2] = core::array::from_fn(|i| {
            core::array::from_fn(|k| (2.0 / 3.0) * end[i][k] + (1.0 / 3.0) * end[1 - i][k])
        });
        for k in 0..2 {
            if end.iter().any(|v| v[k] > 0.0) && nodal.iter().any(|v| v[k] <= 0.0) {
                return Err(invalid(
                    "graded nodal reference inertia must be representable",
                ));
            }
        }
        Ok(nodal)
    }

    /// Closed-form endpoint-weighted Gram matrices. They integrate P1 density
    /// times the original cubic Hermite velocities (degree seven), and P1 spin
    /// times their quadratic slopes (degree five), without runtime quadrature.
    pub(super) fn consistent_local(self, length: f64) -> [[f64; 6]; 6] {
        let end = self.endpoint_inventory(length);
        let mut local = [[0.0; 6]; 6];
        local[0][0] = end[0][0] / 2.0 + end[1][0] / 6.0;
        local[0][3] = (end[0][0] + end[1][0]) / 6.0;
        local[3][0] = local[0][3];
        local[3][3] = end[0][0] / 6.0 + end[1][0] / 2.0;
        let translation = [
            [
                [240.0, 30.0, 54.0, -14.0],
                [30.0, 5.0, 12.0, -3.0],
                [54.0, 12.0, 72.0, -14.0],
                [-14.0, -3.0, -14.0, 3.0],
            ],
            [
                [72.0, 14.0, 54.0, -12.0],
                [14.0, 3.0, 14.0, -3.0],
                [54.0, 14.0, 240.0, -30.0],
                [-12.0, -3.0, -30.0, 5.0],
            ],
        ];
        let slope = [
            [
                [36.0, 0.0, -36.0, 6.0],
                [0.0, 6.0, 0.0, -1.0],
                [-36.0, 0.0, 36.0, -6.0],
                [6.0, -1.0, -6.0, 2.0],
            ],
            [
                [36.0, 6.0, -36.0, 0.0],
                [6.0, 2.0, -6.0, -1.0],
                [-36.0, -6.0, 36.0, 0.0],
                [0.0, -1.0, 0.0, 6.0],
            ],
        ];
        let transverse = [1, 2, 4, 5];
        let h_scale = [1.0, length, 1.0, length];
        let d_scale = [1.0 / length, 1.0, 1.0 / length, 1.0];
        for i in 0..4 {
            for j in i..4 {
                let centerline = end[0][0] * (translation[0][i][j] / 420.0)
                    + end[1][0] * (translation[1][i][j] / 420.0);
                let spin =
                    end[0][1] * (slope[0][i][j] / 30.0) + end[1][1] * (slope[1][i][j] / 30.0);
                let mut value = centerline * h_scale[i] * h_scale[j];
                if spin != 0.0 {
                    value += spin * d_scale[i] * d_scale[j];
                }
                local[transverse[i]][transverse[j]] = value;
                local[transverse[j]][transverse[i]] = value;
            }
        }
        local
    }
}

impl Frame2d {
    /// Add a member retaining its original P1 line mass and section spin.
    ///
    /// The native consistent mass integrates the density--velocity correlation
    /// exactly; gravity, modes, impulses and Newmark all use that SAME owner.
    /// Explicit lumping uses `L*(2*value_i+value_j)/6` at endpoint i, preserving
    /// total inventory and its first axial moment. Both retain the undeformed
    /// reference geometry. No cross-section eccentricity, shear deformation or
    /// configuration-dependent inertia is introduced.
    ///
    /// Uniform profiles dispatch to the existing constant-coefficient path.
    /// Invalid geometry, section count or unrepresentable inventory refuses
    /// before changing the frame.
    pub fn add_element_with_inertia_profile(
        &mut self,
        i: usize,
        j: usize,
        np: usize,
        profile: FrameInertiaProfile,
        make_section: &dyn Fn() -> SectionModel,
    ) -> Result<usize, SolidError> {
        let mass = profile.line_mass_kg_m;
        let spin = profile.section_spin_kg_m;
        if mass[0] == mass[1] && spin[0] == spin[1] {
            return self.add_element_with_rotary_inertia(i, j, np, mass[0], spin[0], make_section);
        }
        self.add_element_with_reference_inertia(
            i,
            j,
            np,
            f64::midpoint(mass[0], mass[1]),
            f64::midpoint(spin[0], spin[1]),
            Some(profile),
            make_section,
        )
    }

    /// Complete original endpoint inertia, including the uniform-member limit.
    /// The ordering is the original member's i-to-j reference coordinate.
    pub fn reference_member_inertia(
        &self,
        element: usize,
    ) -> Result<FrameInertiaProfile, SolidError> {
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("reference member does not exist"))?;
        Ok(member
            .reference_inertia_profile
            .unwrap_or(FrameInertiaProfile {
                line_mass_kg_m: [member.mass_per_length; 2],
                section_spin_kg_m: [member.rotary_inertia_per_length; 2],
            }))
    }
}
