//! Transverse-force stress in a homogeneous isotropic circular section.
//!
//! This is the Saint-Venant flexure solution for a solid circle or concentric
//! circular annulus, including its Poisson effect and both traction-free
//! contours. Original E/G/nu roles must be explicitly isotropic. Subdivision
//! into physical layers is allowed only when their elastic coefficients agree;
//! different densities, strengths and source receipts remain distinct.
//!
//! The supplied force is a section resultant. Its local stress is NOT obtained
//! from the rod's corrected uniform shear strain. The authored kappa is neither
//! changed nor used here; the recovered local energy must not replace/add to
//! that reduced model's storage. This is a prismatic, small-strain, elastic
//! section observation, away from end/load-introduction effects. It supplies
//! no cross-sectional inertia, interlayer mismatch, plasticity or failure law.
//!
//! Reference: Di Paola, Pirrotta and Santoro, Acta Mechanica 217 (2011),
//! pp. 101–118, sections 2 and 5.2, Eq. (28),
//! <https://doi.org/10.1007/s00707-010-0376-8>.

use std::sync::Arc;

use super::{ConcentricShearSection, MaterialStatePointError, invalid};
use crate::state_point::ResolvedIsotropicElasticModuli;

/// Bonded concentric elastic layers with an explicitly common Poisson ratio.
pub mod layered;

fn original_moduli(
    section: &ConcentricShearSection,
    index: usize,
) -> Result<(f64, f64, f64), MaterialStatePointError> {
    let elastic = section.elastic();
    let bindings = elastic.isotropic_parameters().ok_or_else(|| {
        invalid("circular transverse stress requires explicit original isotropic source roles")
    })?;
    let moduli = ResolvedIsotropicElasticModuli::try_from_resolved(
        elastic.layers()[index].material.clone(),
        &bindings[index].parameters,
    )?;
    let coefficients = (
        moduli.young_modulus_pa(),
        moduli.shear_modulus_pa(),
        moduli.poisson_ratio(),
    );
    if coefficients.0 != elastic.layer_integrals()[index].young_modulus_pa()
        || coefficients.1 != section.layer_shear_modulus_pa(index)?
    {
        return Err(invalid(
            "circular transverse stress must retain the original E/G roles",
        ));
    }
    Ok(coefficients)
}

/// Prepared exact circular-section field with the complete original sources.
/// Admission is O(layers); point recovery is allocation-free and O(1).
#[derive(Clone, Debug, PartialEq)]
pub struct CircularTransverseShearRecovery {
    section: Arc<ConcentricShearSection>,
    poisson_ratio: f64,
    shear_modulus_pa: f64,
    scale: f64,
    bore_ratio_squared: f64,
    layer_peak_factors: Vec<f64>,
}

/// Physical transverse shear on the axial cut, in the original material frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseShearStressPoint {
    /// Undoubled [tau_xy,tau_xz] [Pa]; x is the member axis.
    pub shear_stress_pa: [f64; 2],
    /// Engineering [gamma_xy,gamma_xz], using the original uncorrected G.
    pub engineering_shear_strain: [f64; 2],
    /// Local elastic shear storage, tau dot gamma / 2 [J/m^3].
    pub stored_energy_density_j_m3: f64,
    /// Exact maximum |gamma| over the complete selected annulus for this force.
    /// A query on a stress-free contour cannot hide another point's overstrain.
    pub maximum_layer_engineering_shear: f64,
}

impl TransverseShearStressPoint {
    /// Tensor-Voigt [xx,yy,zz,xy,yz,zx], containing this contribution only.
    #[must_use]
    pub fn stress_tensor_pa(&self) -> [f64; 6] {
        [
            0.0,
            0.0,
            0.0,
            self.shear_stress_pa[0],
            0.0,
            self.shear_stress_pa[1],
        ]
    }
}

impl CircularTransverseShearRecovery {
    /// Admit the original explicit isotropic E/nu, E/G or K/G parameterization.
    /// Equal sampled coefficients allow elastic material subdivision; a legacy
    /// independent E/G bundle is not silently declared isotropic. Heterogeneous
    /// elastic layers refuse rather than receiving a homogenized stress field.
    pub fn try_new(section: Arc<ConcentricShearSection>) -> Result<Self, MaterialStatePointError> {
        let elastic = section.elastic();
        let mut reference = None;
        for index in 0..elastic.layers().len() {
            let coefficients = original_moduli(&section, index)?;
            if reference.is_some_and(|original| original != coefficients) {
                return Err(invalid(
                    "circular transverse stress requires one homogeneous original isotropic elastic law",
                ));
            }
            reference = Some(coefficients);
        }
        let (_, shear_modulus_pa, poisson_ratio) = reference.expect("nonempty original section");
        let bore_ratio = elastic.inner_radius_m() / elastic.outer_radius_m();
        let bore_ratio_squared = bore_ratio * bore_ratio;
        // Stress is expressed in Q/A and dimensionless geometry. Avoid fourth
        // powers/inverse radii, including the small-bore Laurent singularity.
        let scale = (0.5 / (1.0 + poisson_ratio)) / (1.0 + bore_ratio_squared);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid(
                "circular transverse stress normalization is unrepresentable",
            ));
        }
        let a = 3.0 + 2.0 * poisson_ratio;
        let b = 1.0 - 2.0 * poisson_ratio;
        let radial = |r: f64| {
            let s = (r / elastic.outer_radius_m()).powi(2);
            let inner = if r == 0.0 {
                0.0
            } else {
                (elastic.inner_radius_m() / r).powi(2)
            };
            a * (1.0 - s) * (1.0 - inner)
        };
        let tangential = |r: f64| {
            let s = (r / elastic.outer_radius_m()).powi(2);
            let inner = if r == 0.0 {
                0.0
            } else {
                (elastic.inner_radius_m() / r).powi(2)
            };
            a * (1.0 + bore_ratio_squared + inner) - b * s
        };
        let mut inner = elastic.inner_radius_m();
        let mut layer_peak_factors = Vec::with_capacity(elastic.layers().len());
        for layer in elastic.layers() {
            let outer = layer.outer_radius_m;
            // In s=(r/R)^2 the radial coefficient is concave with its maximum
            // at r=sqrt(R*Ri). The tangential coefficient decreases strictly
            // for every stable isotropic nu<1/2, so its absolute extrema are
            // endpoints. These are the radial/tangential eigenvalues of the
            // local force-to-stress map, hence an exact whole-annulus peak.
            let extremum = (elastic.outer_radius_m() * bore_ratio.sqrt()).clamp(inner, outer);
            let peak = radial(extremum)
                .max(tangential(inner).abs())
                .max(tangential(outer).abs());
            let factor = scale * peak;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(invalid("circular transverse shear peak is unrepresentable"));
            }
            layer_peak_factors.push(factor);
            inner = outer;
        }
        Ok(Self {
            section,
            poisson_ratio,
            shear_modulus_pa,
            scale,
            bore_ratio_squared,
            layer_peak_factors,
        })
    }

    /// Original immutable section, including all material layers and authored kappa.
    #[must_use]
    pub fn section(&self) -> &ConcentricShearSection {
        &self.section
    }

    /// Recover the exact section stress at a point in one named original layer.
    /// Both sides of an interface and both physical contours are valid. Points
    /// inside the bore or outside the exterior are refused.
    /// The finite bound in (0,1) applies to the WHOLE selected annulus, including
    /// directions different from the queried point. Resultant components are
    /// [Qy,Qz] in N, coordinates [y,z] in m, in the same material frame.
    pub fn evaluate(
        &self,
        layer: usize,
        coordinates_yz_m: [f64; 2],
        resultant_n: [f64; 2],
        maximum_engineering_shear: f64,
    ) -> Result<TransverseShearStressPoint, MaterialStatePointError> {
        let elastic = self.section.elastic();
        let geometry = elastic
            .layers()
            .get(layer)
            .ok_or_else(|| invalid("transverse stress layer out of range"))?;
        let inner = if layer == 0 {
            elastic.inner_radius_m()
        } else {
            elastic.layers()[layer - 1].outer_radius_m
        };
        let [y, z] = coordinates_yz_m;
        let r = y.hypot(z);
        if !y.is_finite()
            || !z.is_finite()
            || !r.is_finite()
            || r < inner
            || r > geometry.outer_radius_m
            || resultant_n.iter().any(|v| !v.is_finite())
            || !maximum_engineering_shear.is_finite()
            || maximum_engineering_shear <= 0.0
            || maximum_engineering_shear >= 1.0
        {
            return Err(invalid(
                "finite circular shear force, original material point and small-shear bound required",
            ));
        }
        let average = resultant_n.map(|q| q / elastic.area_m2());
        let magnitude = average[0].hypot(average[1]);
        let peak_stress = self.layer_peak_factors[layer] * magnitude;
        let peak_strain = peak_stress / self.shear_modulus_pa;
        if !peak_strain.is_finite()
            || peak_strain > maximum_engineering_shear
            || (resultant_n != [0.0; 2] && (magnitude == 0.0 || peak_strain == 0.0))
        {
            return Err(invalid(
                "selected annulus exceeds its circular transverse-shear domain",
            ));
        }
        let a = 3.0 + 2.0 * self.poisson_ratio;
        let b = 1.0 - 2.0 * self.poisson_ratio;
        let shear_stress_pa = if r == 0.0 {
            average.map(|q| (self.scale * a) * q)
        } else {
            let s = (r / elastic.outer_radius_m()).powi(2);
            let inner = (elastic.inner_radius_m() / r).powi(2);
            // Factored radial field gives zero traction on BOTH physical
            // contours without subtracting large, nearly equal Laurent terms.
            let radial = self.scale * (a * (1.0 - s) * (1.0 - inner));
            let tangential = self.scale * (a * (1.0 + self.bore_ratio_squared + inner) - b * s);
            let [ny, nz] = [y / r, z / r];
            let q_radial = ny * average[0] + nz * average[1];
            let q_tangential = -nz * average[0] + ny * average[1];
            [
                radial * q_radial * ny - tangential * q_tangential * nz,
                radial * q_radial * nz + tangential * q_tangential * ny,
            ]
        };
        let engineering_shear_strain = shear_stress_pa.map(|tau| tau / self.shear_modulus_pa);
        let energy = shear_stress_pa
            .iter()
            .zip(engineering_shear_strain)
            .map(|(tau, gamma)| (0.5 * tau) * gamma)
            .sum::<f64>();
        if shear_stress_pa
            .iter()
            .chain(&engineering_shear_strain)
            .any(|v| !v.is_finite())
            || !energy.is_finite()
            || energy < 0.0
            || shear_stress_pa
                .iter()
                .zip(engineering_shear_strain)
                .any(|(tau, gamma)| *tau != 0.0 && gamma == 0.0)
            || (shear_stress_pa != [0.0; 2] && energy == 0.0)
        {
            return Err(invalid(
                "circular transverse stress, strain or local storage is unrepresentable",
            ));
        }
        Ok(TransverseShearStressPoint {
            shear_stress_pa,
            engineering_shear_strain,
            stored_energy_density_j_m3: energy,
            maximum_layer_engineering_shear: peak_strain,
        })
    }
}
