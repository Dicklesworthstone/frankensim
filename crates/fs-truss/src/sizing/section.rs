//! Homothetic cross-section families for the existing axial/Euler sizing model.
//!
//! Physical material area A scales both transverse dimensions together, keeping
//! the declared aspect or diameter ratio fixed. The weak-axis second moment is
//! I_min = c A². These are geometric identities, not local-wall, torsional,
//! distortional, inelastic or code-based buckling resistances.

use crate::ground::TrussConstructionError;

/// Shape before admission; dimensions are represented by physical material area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SectionFamily {
    SolidSquare,
    SolidCircle,
    /// Major side divided by minor side, finite and at least one. The minor
    /// side is sqrt(A / aspect_ratio); Euler uses the weaker principal axis.
    SolidRectangle {
        aspect_ratio: f64,
    },
    /// Inner diameter divided by outer diameter, finite and in [0, 1).
    /// A = pi D² (1 - ratio²) / 4 is wall material area, not enclosed area.
    CircularTube {
        inner_diameter_ratio: f64,
    },
}

/// An admitted fixed shape, scaled homothetically by the selected material area.
///
/// The default is a solid square. For other families sqrt(A) is only an
/// equivalent-square side; it is not a diameter, width or wall thickness.
/// Shape ratios remain fixed when catalog snapping or self-weight changes A.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HomotheticSection {
    family: SectionFamily,
    inertia_coefficient: f64,
}

impl Default for HomotheticSection {
    fn default() -> Self {
        Self::solid_square()
    }
}

impl HomotheticSection {
    /// Original square-section geometry, with I_min / A² = 1/12.
    #[must_use]
    pub const fn solid_square() -> Self {
        Self {
            family: SectionFamily::SolidSquare,
            inertia_coefficient: 1.0 / 12.0,
        }
    }

    /// Admit a section family and its exact analytic weak-axis shape relation.
    /// Unit-aspect rectangles canonicalize to squares; zero-bore tubes to
    /// solid circles. No minimum wall thickness or fabrication limit is inferred.
    ///
    /// # Errors
    /// Nonfinite or out-of-domain ratios, or an unrepresentable positive
    /// coefficient, return a structured construction refusal.
    pub fn try_new(family: SectionFamily) -> Result<Self, TrussConstructionError> {
        let circle = || Self {
            family: SectionFamily::SolidCircle,
            inertia_coefficient: 1.0 / (4.0 * std::f64::consts::PI),
        };
        let inertia_coefficient = match family {
            SectionFamily::SolidSquare => return Ok(Self::solid_square()),
            SectionFamily::SolidCircle => return Ok(circle()),
            SectionFamily::SolidRectangle { aspect_ratio } => {
                if !aspect_ratio.is_finite() || aspect_ratio < 1.0 {
                    return Err(invalid(
                        "rectangle aspect ratio",
                        "requires finite major/minor >= 1",
                    ));
                }
                if aspect_ratio == 1.0 {
                    return Ok(Self::solid_square());
                }
                // A=b*h, b/h=r, I_min=b*h³/12 => I_min/A²=1/(12r).
                (1.0 / 12.0) / aspect_ratio
            }
            SectionFamily::CircularTube {
                inner_diameter_ratio: ratio,
            } => {
                if !ratio.is_finite() || !(0.0..1.0).contains(&ratio) {
                    return Err(invalid(
                        "tube diameter ratio",
                        "requires finite inner/outer in [0,1)",
                    ));
                }
                if ratio == 0.0 {
                    return Ok(circle());
                }
                // I=pi(D⁴-d⁴)/64. Factor 1-r² to retain thin-wall geometry
                // without subtracting two almost equal squared diameters.
                (1.0 + ratio * ratio) / (4.0 * std::f64::consts::PI * (1.0 - ratio) * (1.0 + ratio))
            }
        };
        if !inertia_coefficient.is_finite() || inertia_coefficient <= 0.0 {
            return Err(invalid(
                "section inertia coefficient",
                "must be finite and strictly positive",
            ));
        }
        Ok(Self {
            family,
            inertia_coefficient,
        })
    }

    /// Original admitted family, including its fixed shape parameter.
    #[must_use]
    pub const fn family(self) -> SectionFamily {
        self.family
    }

    /// Dimensionless weak-axis coefficient c = I_min / A².
    #[must_use]
    pub const fn weak_axis_inertia_coefficient(self) -> f64 {
        self.inertia_coefficient
    }

    /// Actual weakest principal second moment [m⁴] at material area A [m²].
    ///
    /// # Errors
    /// Nonpositive/nonfinite area or unrepresentable positive second moment.
    pub fn minimum_second_moment_m4(self, area_m2: f64) -> Result<f64, TrussConstructionError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "section area",
                "must be finite and strictly positive in m2",
            ));
        }
        let inertia = area_m2 * (self.inertia_coefficient * area_m2);
        if !inertia.is_finite() || inertia <= 0.0 {
            return Err(invalid(
                "section second moment",
                "must be representable and strictly positive in m4",
            ));
        }
        Ok(inertia)
    }

    pub(super) const fn name(self) -> &'static str {
        match self.family {
            SectionFamily::SolidSquare => "solid-square",
            SectionFamily::SolidCircle => "solid-circle",
            SectionFamily::SolidRectangle { .. } => "solid-rectangle",
            SectionFamily::CircularTube { .. } => "circular-tube",
        }
    }
}

fn invalid(field: &'static str, requirement: &'static str) -> TrussConstructionError {
    TrussConstructionError::InvalidInput { field, requirement }
}
