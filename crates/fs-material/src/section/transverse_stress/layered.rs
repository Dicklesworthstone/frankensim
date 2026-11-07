//! Transverse shear in perfectly bonded, radially layered circular sections.
//!
//! Each layer retains its own G and E. One common, explicitly sourced isotropic
//! Poisson ratio admits the same in-plane bending/Poisson displacement field;
//! axial warping and radial shear traction are matched at every interface.
//! Different Poisson ratios require extra coupled in-plane fields and refuse.
//! This is not a homogeneous stress field evaluated using a local modulus.
//!
//! The analytic annular field uses one O(layers) tridiagonal preparation and
//! O(1), allocation-free point evaluation. Original density, authored kappa,
//! reference stiffness and source identity are not changed. Local stress energy
//! is not added to, or identified with, the reduced rod's stored energy.
//!
//! Governing Saint-Venant transmission conditions: Stavole et al., IJSS 317
//! (2025), 113407, section 2, doi:10.1016/j.ijsolstr.2025.113407. The coefficients
//! here solve its common-Poisson-ratio subcase through the exact annular harmonic
//! boundary map, retaining the heterogeneous shear traction distribution.

mod coefficients;

use super::{
    ConcentricShearSection, MaterialStatePointError, TransverseShearStressPoint, invalid,
    original_moduli,
};
use coefficients::LayerCoefficients;
use std::sync::Arc;

/// Original layered stress response with independently retained source moduli.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredCircularTransverseShearRecovery {
    section: Arc<ConcentricShearSection>,
    poisson_ratio: f64,
    shear_moduli_pa: Vec<f64>,
    layers: Vec<LayerCoefficients>,
}

impl LayeredCircularTransverseShearRecovery {
    /// Resolve the bounded original layers. All must explicitly share nu, but
    /// may have different elastic stiffness, density and calibration identities.
    /// An unsupported Poisson mismatch, missing role or unresolved numerical
    /// transmission problem refuses; no stiffness or nu is averaged or clipped.
    pub fn try_new(section: Arc<ConcentricShearSection>) -> Result<Self, MaterialStatePointError> {
        let elastic = section.elastic();
        let mut nu = None;
        let mut shear_moduli_pa = Vec::with_capacity(elastic.layers().len());
        let mut radii = Vec::with_capacity(elastic.layers().len() + 1);
        radii.push(elastic.inner_radius_m() / elastic.outer_radius_m());
        if elastic.inner_radius_m() > 0.0 && radii[0] == 0.0 {
            return Err(invalid(
                "original bore is unrepresentable in normalized circular geometry",
            ));
        }
        for index in 0..elastic.layers().len() {
            let (_, g, poisson) = original_moduli(&section, index)?;
            if nu.is_some_and(|value| value != poisson) {
                return Err(invalid(
                    "layered circular flexure requires a common original Poisson ratio",
                ));
            }
            nu = Some(poisson);
            shear_moduli_pa.push(g);
            radii.push(elastic.layers()[index].outer_radius_m / elastic.outer_radius_m());
        }
        let poisson_ratio = nu.expect("original nonempty section");
        let largest = shear_moduli_pa.iter().copied().fold(0.0_f64, f64::max);
        let relative: Vec<_> = shear_moduli_pa.iter().map(|g| g / largest).collect();
        let layers = coefficients::prepare(&radii, &relative, poisson_ratio).ok_or_else(|| {
            invalid("layered circular shear transmission is unrepresentable or unresolved")
        })?;
        Ok(Self {
            section,
            poisson_ratio,
            shear_moduli_pa,
            layers,
        })
    }

    /// The exact original section, including its independent authored correction.
    pub fn section(&self) -> &ConcentricShearSection {
        &self.section
    }

    /// Observe one original material point under section resultant [Qy,Qz] N.
    /// Either material side of a bonded interface can be selected explicitly:
    /// radial traction is continuous while tangential stress may jump with G.
    /// The bound covers the maximum engineering shear over the entire selected
    /// annulus. This prismatic small-strain observation excludes end effects,
    /// cross-sectional inertia, imperfect bonding and varying Poisson ratios.
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
        if !r.is_finite()
            || r < inner
            || r > geometry.outer_radius_m
            || resultant_n.iter().any(|q| !q.is_finite())
            || !maximum_engineering_shear.is_finite()
            || maximum_engineering_shear <= 0.0
            || maximum_engineering_shear >= 1.0
        {
            return Err(invalid(
                "finite original layered point, shear force and small-strain bound required",
            ));
        }
        let average = resultant_n.map(|q| q / elastic.area_m2());
        let g = self.shear_moduli_pa[layer];
        let magnitude = average[0].hypot(average[1]);
        let peak = (self.layers[layer].peak * magnitude) / g;
        if !peak.is_finite()
            || peak > maximum_engineering_shear
            || (resultant_n != [0.0; 2] && (magnitude == 0.0 || peak == 0.0))
        {
            return Err(invalid(
                "selected layered annulus exceeds its transverse-shear domain",
            ));
        }
        let [radial, tangential] =
            self.layers[layer].factors(r / elastic.outer_radius_m(), self.poisson_ratio);
        let stress = if r == 0.0 {
            average.map(|q| radial * q)
        } else {
            let [ny, nz] = [y / r, z / r];
            let qr = ny * average[0] + nz * average[1];
            let qt = -nz * average[0] + ny * average[1];
            [
                radial * qr * ny - tangential * qt * nz,
                radial * qr * nz + tangential * qt * ny,
            ]
        };
        let strain = stress.map(|tau| tau / g);
        let energy: f64 = stress
            .iter()
            .zip(strain)
            .map(|(tau, gamma)| (0.5 * tau) * gamma)
            .sum();
        if stress.iter().chain(&strain).any(|v| !v.is_finite())
            || !energy.is_finite()
            || energy < 0.0
            || stress
                .iter()
                .zip(strain)
                .any(|(tau, gamma)| *tau != 0.0 && gamma == 0.0)
            || (stress != [0.0; 2] && energy == 0.0)
        {
            return Err(invalid(
                "layered circular stress, strain or storage is unrepresentable",
            ));
        }
        Ok(TransverseShearStressPoint {
            shear_stress_pa: stress,
            engineering_shear_strain: strain,
            stored_energy_density_j_m3: energy,
            maximum_layer_engineering_shear: peak,
        })
    }
}
