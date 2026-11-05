//! Original reference-area mass moments, retained independently of stiffness.

use super::{
    MitcBudget, MitcShellModel, NodalDirectors, PlateError, ShearShellSection, ShellConstraint,
    ShellMass, ShellMesh, assemble_mitc3p_impl,
};
use fs_material::laminate::LaminateSection;

/// Physical through-thickness inertia about the mesh reference surface.
///
/// `I0 = integral rho dz` [kg/m²], `I1 = integral rho*z dz` [kg/m], and
/// `I2 = integral rho*z² dz` [kg] retain mass, its offset and rotary inertia.
/// The strictly positive moment matrix `[I0 I1; I1 I2]` is required. Linear
/// interpolation of admitted corner matrices preserves this condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShearShellInertia {
    moments: [f64; 3],
}

impl ShearShellInertia {
    /// Admit a realizable positive reference inventory without changing its
    /// declared surface or replacing its signed first moment by zero.
    pub fn new(i0_kg_m2: f64, i1_kg_m: f64, i2_kg: f64) -> Result<Self, PlateError> {
        let out = Self {
            moments: [i0_kg_m2, i1_kg_m, i2_kg],
        };
        out.validate()?;
        Ok(out)
    }

    /// Copy the physical inertia already admitted by a native shell section.
    #[must_use]
    pub fn from_section(section: &ShearShellSection) -> Self {
        Self {
            moments: [
                section.areal_mass_kg_m2(),
                section.mass_first_moment_kg_m(),
                section.rotary_inertia_kg(),
            ],
        }
    }

    /// Preserve the source laminate's layer order and reference surface.
    pub fn from_laminate(section: &LaminateSection) -> Result<Self, PlateError> {
        Self::new(section.i0_kg_m2(), section.i1_kg_m(), section.i2_kg())
    }

    /// Zeroth thickness moment [kg/m²].
    #[must_use]
    pub const fn areal_mass_kg_m2(self) -> f64 {
        self.moments[0]
    }

    /// Signed first thickness moment [kg/m].
    #[must_use]
    pub const fn mass_first_moment_kg_m(self) -> f64 {
        self.moments[1]
    }

    /// Second thickness moment [kg].
    #[must_use]
    pub const fn rotary_inertia_kg(self) -> f64 {
        self.moments[2]
    }

    pub(super) fn validate(self) -> Result<(), PlateError> {
        let [i0, i1, i2] = self.moments;
        if self.moments.iter().any(|v| !v.is_finite())
            || i0 <= 0.0
            || i2 <= 0.0
            || !positive_determinant(i0, i1, i2)
        {
            return Err(PlateError::BadSection {
                what: "positive-definite original shell inertia moments required",
            });
        }
        Ok(())
    }
}

/// Compare I0*I2 > I1² exactly for the admitted finite f64 inputs. Products of
/// their binary significands fit in 106 bits; separate exponents avoid both
/// overflow and underflow. There is no artificial conditioning cutoff, and a
/// singular block cannot enter through rounded square roots or division.
fn positive_determinant(i0: f64, i1: f64, i2: f64) -> bool {
    if i1 == 0.0 {
        return true;
    }
    let binary = |value: f64| {
        let bits = value.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        if exponent == 0 {
            (u128::from(fraction), -1074)
        } else {
            (u128::from(fraction | (1_u64 << 52)), exponent - 1075)
        }
    };
    let (m0, e0) = binary(i0);
    let (m1, e1) = binary(i1);
    let (m2, e2) = binary(i2);
    let (a, b) = (m0 * m2, m1 * m1);
    let (a_bits, b_bits) = (128 - a.leading_zeros(), 128 - b.leading_zeros());
    let (a_top, b_top) = (e0 + e2 + a_bits as i32, 2 * e1 + b_bits as i32);
    if a_top != b_top {
        return a_top > b_top;
    }
    let width = a_bits.max(b_bits);
    (a << (width - a_bits)) > (b << (width - b_bits))
}

/// Assemble the native MITC shell with P1 original section inertia per facet.
///
/// Corner order follows `mesh.tris`; the moments are about the same reference
/// surface and oriented directors as `sections`. This field is the authoritative
/// inertia, replacing the sections' constant mass moments while retaining their
/// original elastic A/B/D/shear operators. It can describe spatially graded
/// density or layers without collapsing their mass centroids to a facet mean.
///
/// Consistent mass exactly integrates the degree-seven enriched velocity
/// products and retains the original bubble condensation. Explicit lumping uses
/// `integral N_i * [I0,I1,I2] dA`, preserving total inventory and its signed first
/// moment. Enriched transient dynamics retain the same field permanently;
/// changing temperature or stiffness cannot replace the original inventory.
#[allow(clippy::too_many_arguments)]
pub fn assemble_mitc3p_with_reference_inertia(
    mesh: &ShellMesh,
    directors: &NodalDirectors,
    sections: &[ShearShellSection],
    reference_inertia: &[[ShearShellInertia; 3]],
    constraints: &[ShellConstraint],
    mass: ShellMass,
    budget: MitcBudget,
    drilling_alpha: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<MitcShellModel, PlateError> {
    assemble_mitc3p_impl(
        mesh,
        directors,
        sections,
        Some(reference_inertia),
        constraints,
        mass,
        budget,
        drilling_alpha,
        cancelled,
    )
}

/// An element-local view; constant input keeps the original arithmetic path.
#[derive(Clone, Copy)]
pub(super) enum FacetInertia<'a> {
    Constant(ShearShellInertia),
    Linear(&'a [ShearShellInertia; 3]),
}

impl<'a> FacetInertia<'a> {
    pub(super) fn new(
        section: &ShearShellSection,
        corners: Option<&'a [ShearShellInertia; 3]>,
    ) -> Self {
        match corners {
            Some(corners) if corners[1] == corners[0] && corners[2] == corners[0] => {
                Self::Constant(corners[0])
            }
            Some(corners) => Self::Linear(corners),
            None => Self::Constant(ShearShellInertia::from_section(section)),
        }
    }

    pub(super) fn at(self, lambda: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Constant(value) => value.moments,
            Self::Linear(corners) => {
                core::array::from_fn(|k| (0..3).map(|i| lambda[i] * corners[i].moments[k]).sum())
            }
        }
    }

    /// Exact positive partition of each P1 moment, including the drilling
    /// gauge's I2 weight. The uniform limit preserves area/3 arithmetic.
    pub(super) fn at_node(self, node: usize, area: f64) -> [f64; 3] {
        match self {
            Self::Constant(value) => value.moments.map(|v| area / 3.0 * v),
            Self::Linear(corners) => core::array::from_fn(|k| {
                (0..3)
                    .map(|i| area / if i == node { 6.0 } else { 12.0 } * corners[i].moments[k])
                    .sum()
            }),
        }
    }
}
