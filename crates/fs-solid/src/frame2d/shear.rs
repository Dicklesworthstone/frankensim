//! Source-bound transverse shear in the original force-based frame owner.
//!
//! The shear resultant is `Q = -(M_i + M_j)/L` in the member's local +y
//! direction. Its complementary energy `Q² L/(2 GA)` gives the exact basic
//! flexibility contribution; GA includes the source's authored shear correction.
//! Consistent inertia uses the same static Timoshenko interpolation, with
//! independent section rotation. Using the centerline slope for section spin
//! would incorrectly reintroduce zero shear.

use super::{Frame2d, SolidError};
use fs_material::section::ConcentricShearSection;

/// Accepted local transverse shear in an explicitly shear-deformable member.
/// These are corrected section resultants, not resolved three-dimensional
/// local stresses or a shear-failure criterion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameShearResponse {
    /// Section force in the member's local +y direction [N].
    pub force_n: f64,
    /// Effective engineering shear `dw/dx - theta = Q/(kappa GA)`.
    pub engineering_strain: f64,
    /// Reversible shear storage over the reference member [J].
    pub stored_energy_j: f64,
}

#[derive(Clone, Debug)]
pub(super) struct ShearElement {
    source: ConcentricShearSection,
    // Basic rotational compliance [1/(N m)], cached at source admission.
    pub(super) flexibility: f64,
    pub(super) bending_fraction: f64,
    pub(super) shear_fraction: f64,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

impl ShearElement {
    fn try_new(source: ConcentricShearSection, length: f64) -> Result<Self, SolidError> {
        for layer in 0..source.elastic().layer_integrals().len() {
            let properties =
                source
                    .source_properties(layer)
                    .map_err(|error| SolidError::MaterialRefused {
                        what: error.to_string(),
                    })?;
            for property in properties {
                let refused = || SolidError::MaterialRefused {
                    what: format!(
                        "frame shear layer {layer} property {} needs a frozen coefficient without a frequency restriction",
                        property.requirement().name()
                    ),
                };
                if !property.is_constant_at_fixed_temperature()
                    || property
                        .angular_frequency_band_rad_s()
                        .map_err(|_| refused())?
                        .is_some()
                {
                    return Err(refused());
                }
            }
        }
        let flexibility = 1.0 / source.shear_stiffness_n() / length;
        let bending = length / source.elastic().bending_stiffness_n_m2() / 12.0;
        if !length.is_finite()
            || length <= 0.0
            || !flexibility.is_finite()
            || flexibility <= 0.0
            || !bending.is_finite()
            || bending <= 0.0
        {
            return Err(invalid(
                "finite positive member length and representable shear/bending compliance required",
            ));
        }
        // beta=1/(1+12 EI/(GA L²)). Divide the smaller compliance by the
        // larger one so neither a large stiffness ratio nor its sum overflows.
        // Keep both fractions independently; 1-beta loses very small shear.
        let (bending_fraction, shear_fraction) = if flexibility <= bending {
            let ratio = flexibility / bending;
            (1.0 / (1.0 + ratio), ratio / (1.0 + ratio))
        } else {
            let ratio = bending / flexibility;
            (ratio / (1.0 + ratio), 1.0 / (1.0 + ratio))
        };
        Ok(Self {
            source,
            flexibility,
            bending_fraction,
            shear_fraction,
        })
    }

    pub(super) fn consistent_mass(
        &self,
        length: f64,
        line_mass: f64,
        section_spin: f64,
    ) -> [[f64; 6]; 6] {
        // Four positive Gauss points integrate cubic translation squared
        // and quadratic section rotation squared exactly (degrees 6 and 4).
        const POINTS: [(f64, f64); 4] = [
            (0.069_431_844_202_973_71, 0.173_927_422_568_726_93),
            (0.330_009_478_207_571_87, 0.326_072_577_431_273_05),
            (0.669_990_521_792_428_1, 0.326_072_577_431_273_05),
            (0.930_568_155_797_026_2, 0.173_927_422_568_726_93),
        ];
        let mut result = [[0.0; 6]; 6];
        let beta = self.bending_fraction;
        for (x, weight) in POINTS {
            let h = beta * x * x * (3.0 - 2.0 * x) + self.shear_fraction * x;
            let v = [
                0.0,
                1.0 - h,
                length * (x - 0.5 * x * x - 0.5 * h),
                0.0,
                h,
                0.5 * length * (x * x - h),
            ];
            let turning = 6.0 * beta * x * (1.0 - x);
            let theta = [
                0.0,
                -turning / length,
                1.0 - x - 0.5 * turning,
                0.0,
                turning / length,
                x - 0.5 * turning,
            ];
            let u = [1.0 - x, 0.0, 0.0, x, 0.0, 0.0];
            let mass = weight * (line_mass * length);
            let spin = weight * (section_spin * length);
            for i in 0..6 {
                for j in 0..=i {
                    result[i][j] +=
                        mass * (u[i] * u[j] + v[i] * v[j]) + spin * theta[i] * theta[j];
                }
            }
        }
        for i in 0..6 {
            for j in 0..i {
                result[j][i] = result[i][j];
            }
        }
        result
    }
}

impl Frame2d {
    /// Add an elastic Timoshenko member from the original layered material
    /// section, including its explicitly authored transverse shear correction.
    /// The existing force-based equilibrium owner adds `Q² L/(2 kappa GA)`;
    /// axial/bending material sources, density and section spin are unchanged.
    ///
    /// Native statics, modes, Newmark, impacts and harmonic/pole responses use
    /// this same stiffness and original inertia. Consistent mass integrates
    /// the exact unloaded static Timoshenko interpolation for translation AND
    /// section rotation. Lumped mass retains the original nodal inventory.
    /// This is a convergent beam discretization, not exact continuum dynamics.
    ///
    /// Sources must be constant at their retained temperatures and support all
    /// frequencies, as for `add_material_element`. No shear viscosity,
    /// warping, shear failure or three-dimensional stress recovery is inferred.
    /// Explicit shear viscosity enters through `FrameViscousSection::try_new_with_shear`
    /// and `add_material_viscous_element`, using this same interpolation.
    /// Longitudinal Maxwell/thermal adapters using Euler–Bernoulli strain
    /// traces refuse these members. Small corotational strains and fixed
    /// reference inertia retain their existing meanings.
    ///
    /// # Errors
    /// Unsupported source dependence, invalid geometry or section count,
    /// unrepresentable compliance or inertia. Refusal leaves `self` unchanged.
    pub fn add_material_shear_element(
        &mut self,
        i: usize,
        j: usize,
        np: usize,
        section: ConcentricShearSection,
    ) -> Result<usize, SolidError> {
        let [Some(a), Some(b)] = [self.nodes.get(i), self.nodes.get(j)] else {
            return Err(invalid("shear element nodes do not exist"));
        };
        let length = (b[0] - a[0]).hypot(b[1] - a[1]);
        let shear = ShearElement::try_new(section, length)?;
        let element = self.add_material_element(i, j, np, shear.source.elastic().clone())?;
        self.elements[element].transverse_shear = Some(shear);
        Ok(element)
    }

    /// Original shear binding, its correction and every layer's material source.
    /// Euler–Bernoulli members and unknown indices return `None`.
    #[must_use]
    pub fn material_shear_section(&self, element: usize) -> Option<&ConcentricShearSection> {
        Some(
            &self
                .elements
                .get(element)?
                .transverse_shear
                .as_ref()?
                .source,
        )
    }

    /// Observe accepted transverse shear through the same native member law.
    /// Local +y rotates with the corotational chord. The effective strain and
    /// energy use the declared corrected section stiffness, not an invented
    /// uniform local stress `Q/A`.
    pub fn element_shear_response(&self, element: usize) -> Result<FrameShearResponse, SolidError> {
        let member = self
            .elements
            .get(element)
            .ok_or_else(|| invalid("shear response member does not exist"))?;
        let shear = member
            .transverse_shear
            .as_ref()
            .ok_or_else(|| invalid("member has no explicit transverse shear law"))?;
        let force_n = -(member.commit.q[1] + member.commit.q[2]) / member.l0;
        let engineering_strain = force_n / shear.source.shear_stiffness_n();
        let stored_energy_j = 0.5 * force_n * engineering_strain * member.l0;
        if ![force_n, engineering_strain, stored_energy_j]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("accepted member shear response is unrepresentable"));
        }
        Ok(FrameShearResponse {
            force_n,
            engineering_strain,
            stored_energy_j,
        })
    }
}
