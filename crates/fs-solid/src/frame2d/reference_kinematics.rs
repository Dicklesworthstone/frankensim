//! Immutable reference strain interpolation for source-bound linear frames.
//!
//! Native Timoshenko stiffness, inertia and viscosity choose the same unloaded
//! static interpolation at construction. A thermal consumer may retain that
//! interpolation while E/G/eta change: this is a fixed Galerkin basis with the
//! original consistent mass, not a newly equilibrated shape with stale inertia.
//! The consumer owns the source potential, parameter work and admissible mass.

use super::{FrameViscousHeatFractions, FrameViscousSection, SolidError};

/// Original member interpolation, retained independently of later coefficients.
/// No mutable material history or allocation is held by this geometric value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameReferenceKinematics {
    length_m: f64,
    bending_fraction: f64,
    shear_fraction: f64,
    shear_correction: Option<f64>,
}

/// Exact constant axial/shear strain and curvature observations for a member.
/// RMS curvature integrates quadratic storage/loss exactly; its maximum is the
/// distinct endpoint quantity needed for an extreme-fibre strain bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameReferenceStrainMeasures {
    pub axial_strain: f64,
    pub transverse_shear_strain: f64,
    pub rms_curvature_per_m: f64,
    pub maximum_absolute_curvature_per_m: f64,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

impl FrameReferenceKinematics {
    pub(super) fn new(section: &FrameViscousSection, length_m: f64) -> Result<Self, SolidError> {
        if !length_m.is_finite() || length_m <= 0.0 {
            return Err(invalid("finite positive frame reference length required"));
        }
        let (bending_fraction, shear_fraction, shear_correction) = match section.transverse() {
            Some(transverse) => {
                let (_, bending, shear) =
                    super::shear::ShearElement::compliance_partition(transverse.shear(), length_m)?;
                (bending, shear, Some(transverse.shear().shear_correction()))
            }
            None => (1.0, 0.0, None),
        };
        Ok(Self {
            length_m,
            bending_fraction,
            shear_fraction,
            shear_correction,
        })
    }

    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Original beta=1/(1+12 EI/(kappa GA L²)), never recomputed at a new T.
    pub const fn bending_fraction(&self) -> f64 {
        self.bending_fraction
    }
    /// Independently retained 1-beta, preserving small representable shear.
    pub const fn shear_fraction(&self) -> f64 {
        self.shear_fraction
    }
    pub const fn has_transverse_shear(&self) -> bool {
        self.shear_correction.is_some()
    }

    /// Coordinates are the same native extensions in metres:
    /// x0=u_j-u_i, x1=2(v_i-v_j)+L(theta_i+theta_j), x2=L(theta_i-theta_j).
    /// Local curvature is -(3 beta(1-2x/L)x1+x2)/L² and gamma=-s*x1/(2L).
    pub fn strain_measures(&self, x: [f64; 3]) -> Result<FrameReferenceStrainMeasures, SolidError> {
        if x.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite original frame strain coordinates required"));
        }
        let l = self.length_m;
        let axial_strain = x[0] / l;
        let transverse_shear_strain = -0.5 * self.shear_fraction * (x[1] / l);
        let symmetric_curvature = 3.0 * self.bending_fraction * (x[1] / l / l);
        let antisymmetric_curvature = x[2] / l / l;
        let rms_curvature_per_m =
            (symmetric_curvature / 3.0_f64.sqrt()).hypot(antisymmetric_curvature);
        let maximum_absolute_curvature_per_m = (symmetric_curvature + antisymmetric_curvature)
            .abs()
            .max((symmetric_curvature - antisymmetric_curvature).abs());
        if [
            axial_strain,
            transverse_shear_strain,
            rms_curvature_per_m,
            maximum_absolute_curvature_per_m,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || (x[0] != 0.0 && axial_strain == 0.0)
            || (self.shear_fraction > 0.0 && x[1] != 0.0 && transverse_shear_strain == 0.0)
            || (x[1] != 0.0 && self.bending_fraction > 0.0 && symmetric_curvature == 0.0)
            || (x[2] != 0.0 && antisymmetric_curvature == 0.0)
        {
            return Err(invalid("frame reference strains are unrepresentable"));
        }
        Ok(FrameReferenceStrainMeasures {
            axial_strain,
            transverse_shear_strain,
            rms_curvature_per_m,
            maximum_absolute_curvature_per_m,
        })
    }

    /// One original layer's [axial, symmetric bending/shear, antisymmetric]
    /// spring coefficients [N/m]. Updated E/G act on the SAME fixed strains.
    /// Axial thermal rest extension is supplied by the original expansion law.
    pub fn layer_stiffness_n_m(
        &self,
        section: &FrameViscousSection,
        layer: usize,
    ) -> Result<[f64; 3], SolidError> {
        self.admit_section(section)?;
        let geometry = section
            .elastic()
            .layer_integrals()
            .get(layer)
            .ok_or_else(|| invalid("frame stiffness layer index out of range"))?;
        let shear = match section.transverse() {
            Some(source) => {
                let modulus = source
                    .shear()
                    .layer_shear_modulus_pa(layer)
                    .map_err(|error| SolidError::MaterialRefused {
                        what: error.to_string(),
                    })?;
                let value = (modulus * geometry.area_m2()) * source.shear().shear_correction();
                if !value.is_finite() || value <= 0.0 {
                    return Err(invalid("original layer shear stiffness is unrepresentable"));
                }
                value
            }
            None => 0.0,
        };
        let parts = self.parts(
            geometry.axial_stiffness_n(),
            geometry.bending_stiffness_n_m2(),
            shear,
        )?;
        Ok([parts[0], parts[1] + parts[2], parts[3]])
    }

    /// Exact native dashpots [N s/m] at these retained reference strains.
    pub fn damping_n_s_m(&self, section: &FrameViscousSection) -> Result<[f64; 3], SolidError> {
        let parts = self.damping_parts(section)?;
        Ok([parts[0], parts[1] + parts[2], parts[3]])
    }

    /// Original layer heat from its independent viscosities on this fixed basis.
    pub fn layer_heat_fractions(
        &self,
        section: &FrameViscousSection,
        layer: usize,
    ) -> Result<FrameViscousHeatFractions, SolidError> {
        let [_, bending, shear, _] = self.damping_parts(section)?;
        let fraction = |part: f64, whole: f64| {
            let value = if whole == 0.0 { 0.0 } else { part / whole };
            if !value.is_finite() || !(0.0..=1.0).contains(&value) || (part > 0.0 && value == 0.0) {
                return Err(invalid(
                    "original frame-layer loss fraction is unrepresentable",
                ));
            }
            Ok(value)
        };
        let axial = *section
            .axial()
            .layer_viscous_rigidity_n_s()
            .get(layer)
            .ok_or_else(|| invalid("frame viscosity layer index out of range"))?;
        let curvature = section.bending().layer_viscous_stiffness_n_m2_s()[layer];
        let axial = fraction(axial, section.axial().viscous_rigidity_n_s())?;
        let curvature = fraction(curvature, section.bending().viscous_stiffness_n_m2_s())?;
        let transverse = match section.transverse() {
            Some(source) => fraction(
                source.layer_viscous_rigidity_n_s()[layer],
                source.viscous_rigidity_n_s(),
            )?,
            None => 0.0,
        };
        let bending_share = fraction(bending, bending + shear)?;
        let shear_share = fraction(shear, bending + shear)?;
        let symmetric_bending = curvature * bending_share;
        let symmetric_shear = transverse * shear_share;
        if (curvature > 0.0 && bending_share > 0.0 && symmetric_bending == 0.0)
            || (transverse > 0.0 && shear_share > 0.0 && symmetric_shear == 0.0)
        {
            return Err(invalid("original frame-layer loss fraction underflow"));
        }
        Ok(FrameViscousHeatFractions {
            axial,
            symmetric_bending,
            symmetric_shear,
            antisymmetric_bending: curvature,
        })
    }

    fn admit_section(&self, section: &FrameViscousSection) -> Result<(), SolidError> {
        match (self.shear_correction, section.transverse()) {
            (None, None) => Ok(()),
            (Some(correction), Some(source))
                if correction.to_bits() == source.shear().shear_correction().to_bits() =>
            {
                Ok(())
            }
            _ => Err(invalid(
                "frame reference shear interpolation and source correction differ",
            )),
        }
    }

    fn damping_parts(&self, section: &FrameViscousSection) -> Result<[f64; 4], SolidError> {
        self.admit_section(section)?;
        self.parts(
            section.axial().viscous_rigidity_n_s(),
            section.bending().viscous_stiffness_n_m2_s(),
            section
                .transverse()
                .map_or(0.0, |source| source.viscous_rigidity_n_s()),
        )
    }

    // One exact geometric integral for elastic and viscous coefficients. Shear
    // rigidity already contains the source section's correction, exactly once.
    fn parts(
        &self,
        axial_rigidity: f64,
        bending_rigidity: f64,
        shear_rigidity: f64,
    ) -> Result<[f64; 4], SolidError> {
        let l = self.length_m;
        let axial = axial_rigidity / l;
        let curvature = bending_rigidity / l / l / l;
        let bending = (curvature * self.bending_fraction) * (3.0 * self.bending_fraction);
        let shear = (shear_rigidity * self.shear_fraction / l) * (0.25 * self.shear_fraction);
        if [axial, curvature, bending + shear]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || (axial_rigidity > 0.0 && axial == 0.0)
            || (bending_rigidity > 0.0 && (curvature == 0.0 || bending == 0.0))
            || (shear_rigidity > 0.0 && shear == 0.0)
        {
            return Err(invalid(
                "frame coefficient is unrepresentable on the original strain basis",
            ));
        }
        Ok([axial, bending, shear, curvature])
    }
}
