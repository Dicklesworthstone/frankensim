//! Graded original line mass in the existing Euler--Bernoulli inertia owner.

use super::{Frame2d, FrameMassFormulation, SolidError, invalid};
use crate::frame2d::SectionModel;

/// P1 original transverse density moments along a member.
///
/// Values follow the element's i-to-j reference coordinate. Line mass is
/// `integral rho dA` [kg/m], first moment is `integral rho*y dA` [kg], and
/// spin inertia is `integral rho*y² dA` [kg m]. Positive y is the reference
/// member's left normal. These are original inventories, not current density
/// sampled on a deformed chord.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameInertiaProfile {
    line_mass_kg_m: [f64; 2],
    section_first_kg: [f64; 2],
    section_spin_kg_m: [f64; 2],
}

impl FrameInertiaProfile {
    /// Admit finite nonnegative endpoint moments. Section spin requires line
    /// mass at the same endpoint; zero moments retain the massless-member lane.
    pub fn new(line_mass_kg_m: [f64; 2], section_spin_kg_m: [f64; 2]) -> Result<Self, SolidError> {
        Self::with_first_moment(line_mass_kg_m, [0.0; 2], section_spin_kg_m)
    }

    /// Admit original line mass, signed first moment and second moment about
    /// the declared reference centerline. Each is interpolated independently
    /// as P1; their ratio is the resulting mass-centroid offset, not an
    /// independently linear offset. Endpoint moment matrices must be positive
    /// semidefinite (`first² <= mass*spin`), which also guarantees positivity
    /// everywhere along the member. No density profile or material is inferred.
    ///
    /// Consistent reference inertia uses the Euler--Bernoulli velocity field
    /// `(u_dot-y*v_dot', v_dot)`, retaining the signed axial/slope coupling.
    /// Nonzero first moments require `ConsistentReference`; diagonal lumping
    /// cannot preserve that coupling and is refused on element insertion.
    pub fn with_first_moment(
        line_mass_kg_m: [f64; 2],
        section_first_kg: [f64; 2],
        section_spin_kg_m: [f64; 2],
    ) -> Result<Self, SolidError> {
        if (0..2).any(|i| {
            !line_mass_kg_m[i].is_finite()
                || line_mass_kg_m[i] < 0.0
                || !section_first_kg[i].is_finite()
                || !section_spin_kg_m[i].is_finite()
                || section_spin_kg_m[i] < 0.0
                || (line_mass_kg_m[i] == 0.0 && section_spin_kg_m[i] != 0.0)
                // Avoid forming mass*spin or first²: either can overflow for
                // finite admissible moments. Zero mass/spin requires zero first.
                || section_first_kg[i].abs()
                    > line_mass_kg_m[i].sqrt() * section_spin_kg_m[i].sqrt()
        }) {
            return Err(invalid(
                "finite positive-semidefinite endpoint density moments required; spin requires mass",
            ));
        }
        Ok(Self {
            line_mass_kg_m,
            section_first_kg,
            section_spin_kg_m,
        })
    }

    /// Original endpoint line masses [kg/m], in member i-to-j order.
    #[must_use]
    pub const fn line_mass_kg_m(self) -> [f64; 2] {
        self.line_mass_kg_m
    }

    /// Original endpoint signed transverse first density moments [kg].
    #[must_use]
    pub const fn section_first_kg(self) -> [f64; 2] {
        self.section_first_kg
    }

    /// Original endpoint spin inertias about the reference centerline [kg m].
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
        if self.section_first_kg != [0.0; 2] {
            // Integral_0^1 N_endpoint * N_axial * dH/ds ds. The negative
            // sign follows u_dot-y*v_dot'; H's rotational columns carry L.
            // Exact degree-four products require no new runtime quadrature.
            let coupling = [
                [
                    [-3.0 / 10.0, 1.0 / 10.0, 3.0 / 10.0, -1.0 / 15.0],
                    [-1.0 / 5.0, -1.0 / 60.0, 1.0 / 5.0, -1.0 / 60.0],
                ],
                [
                    [-1.0 / 5.0, -1.0 / 60.0, 1.0 / 5.0, -1.0 / 60.0],
                    [-3.0 / 10.0, -1.0 / 15.0, 3.0 / 10.0, 1.0 / 10.0],
                ],
            ];
            for (a, axial) in [0, 3].into_iter().enumerate() {
                for (b, transverse) in transverse.into_iter().enumerate() {
                    let value = -(self.section_first_kg[0] * coupling[0][a][b]
                        + self.section_first_kg[1] * coupling[1][a][b])
                        * h_scale[b];
                    local[axial][transverse] = value;
                    local[transverse][axial] = value;
                }
            }
        }
        local
    }
}

impl Frame2d {
    /// Add a member retaining its original P1 section density moments.
    ///
    /// The native consistent mass integrates the density--velocity correlation
    /// exactly; gravity, modes, impulses and Newmark all use that SAME owner.
    /// Explicit lumping uses `L*(2*value_i+value_j)/6` at endpoint i, preserving
    /// total inventory and its first axial moment. Both retain the undeformed
    /// reference geometry. Nonzero signed section first moments retain density
    /// eccentricity and require consistent inertia. Shear deformation and
    /// configuration-dependent inertia are not introduced.
    ///
    /// Uniform centered profiles dispatch to the constant-coefficient path.
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
        if profile.section_first_kg != [0.0; 2]
            && self.mass_formulation != FrameMassFormulation::ConsistentReference
        {
            return Err(invalid(
                "eccentric density moments require consistent reference inertia",
            ));
        }
        if mass[0] == mass[1] && spin[0] == spin[1] && profile.section_first_kg == [0.0; 2] {
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
                section_first_kg: [0.0; 2],
                section_spin_kg_m: [member.rotary_inertia_per_length; 2],
            }))
    }
}
