//! Rankine tensile assessment of bonded circular layers under combined loads.
//!
//! The retained elastic section supplies `epsilon = N/EA` and `kappa = |M|/EI`.
//! Each annulus keeps its own E, G, radius and tensile-strength source. Circular
//! Saint-Venant torsion gives shear magnitude `tau_i = G_i r |T|/GJ`, with the
//! original `GJ = sum(G_i J_i)`, independent of transverse shear correction.
//! At a point the principal tensile stress is
//! `lambda+ = sigma/2 + sqrt((sigma/2)^2 + tau^2)`; the third principal value is
//! zero. The compression branch is rationalized to retain small positive
//! tension without cancellation or squaring dimensional stresses.
//!
//! For fixed r, lambda+ increases with axial normal stress, so the bending
//! tensile side governs. There `sigma_i = sigma_mean_i + E_i kappa r`, and both
//! this normal stress and the shear magnitude increase with r. The annular
//! maximum is therefore exactly at its OUTER radius. Internal layer boundaries
//! remain independent candidates; no angular/radial sampling or eigensolve is
//! needed. Thermal mismatch supplies the native per-layer mean stress.
//!
//! Calling these methods explicitly selects maximum-principal-tension
//! (Rankine) assessment using the capacity's named strength. It does not turn
//! a tensile yield/proof observation into a ductile shear or multiaxial yield
//! law. Transverse-force shear, radial/Poisson restraint, interfaces, warping,
//! plastic redistribution, buckling, fatigue and fracture evolution are absent.
//! The original material sources and small-strain circular-section assumptions
//! remain authoritative; a crossed limit changes no mechanical state.

use super::{AxialTensileCapacity, checked_utilization};
use crate::section::thermal::ConcentricThermoelasticSection;
use crate::section::{ConcentricShearSection, MaterialStatePointError, invalid};

/// Exact normal-stress extrema and combined principal tension on one annulus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerPrincipalTensileAssessment {
    /// Minimum axial normal stress over the layer [Pa], tension positive.
    pub minimum_axial_stress_pa: f64,
    /// Maximum axial normal stress over the layer [Pa], tension positive.
    pub maximum_axial_stress_pa: f64,
    /// Maximum magnitude of the circular torsional shear [Pa].
    pub maximum_torsional_shear_pa: f64,
    /// Maximum principal tensile stress over the entire annulus [Pa].
    pub principal_tensile_stress_pa: f64,
    /// Principal tensile stress divided by this layer's original strength.
    pub tensile_utilization: f64,
}

/// Combined Rankine observations in the capacity's original layer order.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricPrincipalTensileAssessment {
    layers: Vec<LayerPrincipalTensileAssessment>,
    governing_layer: usize,
}

impl ConcentricPrincipalTensileAssessment {
    /// All original layers, including internal material interfaces.
    #[must_use]
    pub fn layers(&self) -> &[LayerPrincipalTensileAssessment] {
        &self.layers
    }

    /// Largest utilization; the first layer wins an exact tie.
    #[must_use]
    pub const fn governing_layer(&self) -> usize {
        self.governing_layer
    }

    /// Maximum principal tensile demand over its own layer's sourced limit.
    #[must_use]
    pub fn utilization(&self) -> f64 {
        self.layers[self.governing_layer].tensile_utilization
    }

    /// Whether at least one layer reaches its limit, without a tolerance.
    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.utilization() >= 1.0
    }

    /// Whether at least one layer exceeds its limit, without a tolerance.
    #[must_use]
    pub fn exceeds_limit(&self) -> bool {
        self.utilization() > 1.0
    }
}

impl AxialTensileCapacity {
    /// Exact annular maximum principal tension under signed axial force [N],
    /// two centroidal bending moments [N m] and signed circular torque [N m].
    /// The shear section must retain the exact original elastic section and
    /// material states. Each layer's actual G supplies torsional load sharing;
    /// the transverse shear correction never enters this field.
    ///
    /// This explicitly selects a Rankine tensile criterion, not a ductile yield
    /// criterion. Pure torsion has principal tension `abs(tau)` even with zero
    /// axial resultant. Compression can reduce that tension but cannot erase it
    /// when torque is nonzero. Loads are assumed to remain in the declared
    /// small-strain linear constitutive regime; this overload supplies no
    /// experimental calibration envelope or inferred strain limit.
    ///
    /// # Errors
    /// Foreign/stale shear sections, nonfinite loads and unrepresentable stress
    /// or utilization refuse atomically. Strength admission and receipts remain
    /// with the original capacity, including explicitly bound source conditions.
    pub fn assess_axial_bending_torsion(
        &self,
        shear: &ConcentricShearSection,
        axial_force_n: f64,
        bending_moment_n_m: [f64; 2],
        torque_n_m: f64,
    ) -> Result<ConcentricPrincipalTensileAssessment, MaterialStatePointError> {
        let strain = quotient(
            axial_force_n,
            self.section().axial_stiffness_n(),
            "combined tensile assessment needs representable signed axial force and strain",
        )?;
        self.assess_combined(shear, bending_moment_n_m, torque_n_m, |i, _, _| {
            product(
                self.section().layer_integrals()[i].young_modulus_pa(),
                strain,
                "combined tensile axial stress is unrepresentable",
            )
        })
    }

    /// Combined principal tension using the original thermal layer stresses,
    /// retaining mismatch at zero resultant rather than replacing it by N/EA.
    /// Both the thermal and shear sections must match the capacity exactly.
    ///
    /// The native response retains its total/thermal/elastic axial strain
    /// checks. The SAME declared limit also bounds both axial fiber extrema
    /// after bending, and the maximum engineering torsional shear strain
    /// `gamma = r_outer * abs(T/GJ)`. This is a componentwise linear-regime
    /// assumption, not an inferred yield or equivalent-strain criterion.
    ///
    /// # Errors
    /// Stale/foreign states, native thermal refusal, strain-envelope escape,
    /// nonfinite loads and unrepresentable arithmetic refuse atomically.
    pub fn assess_thermoelastic_axial_bending_torsion(
        &self,
        thermal: &ConcentricThermoelasticSection,
        shear: &ConcentricShearSection,
        total_axial_strain: f64,
        bending_moment_n_m: [f64; 2],
        torque_n_m: f64,
        linear_strain_limit: f64,
    ) -> Result<ConcentricPrincipalTensileAssessment, MaterialStatePointError> {
        if thermal.section() != self.section() {
            return Err(invalid(
                "combined thermal tensile assessment needs the exact original section and material states",
            ));
        }
        let response = thermal.response_at_axial_strain(total_axial_strain, linear_strain_limit)?;
        self.assess_combined(shear, bending_moment_n_m, torque_n_m, |i, bending, torsion| {
            let native = &response.layers()[i];
            if total_axial_strain.abs() + bending > linear_strain_limit
                || native.elastic_strain().abs() + bending > linear_strain_limit
                || torsion > linear_strain_limit
            {
                return Err(invalid("combined thermal tensile fiber or torsional shear strain exceeds the original linear regime"));
            }
            Ok(native.stress_pa())
        })
    }

    fn assess_combined(
        &self,
        shear: &ConcentricShearSection,
        bending_moment_n_m: [f64; 2],
        torque_n_m: f64,
        mut mean_stress: impl FnMut(usize, f64, f64) -> Result<f64, MaterialStatePointError>,
    ) -> Result<ConcentricPrincipalTensileAssessment, MaterialStatePointError> {
        if shear.elastic() != self.section() {
            return Err(invalid(
                "combined tensile assessment needs the exact original shear section and material states",
            ));
        }
        if bending_moment_n_m.iter().any(|m| !m.is_finite()) {
            return Err(invalid(
                "combined tensile assessment requires finite bending moments",
            ));
        }
        let scale = bending_moment_n_m[0].abs().max(bending_moment_n_m[1].abs());
        let curvature = if scale == 0.0 {
            0.0
        } else {
            let [a, b] = bending_moment_n_m.map(|m| m / scale);
            product(
                quotient(
                    scale,
                    self.section().bending_stiffness_n_m2(),
                    "combined tensile curvature is unrepresentable",
                )?,
                fs_math::det::sqrt(a * a + b * b),
                "combined tensile curvature is unrepresentable",
            )?
        };
        let twist = quotient(
            torque_n_m,
            shear.torsional_stiffness_n_m2(),
            "combined tensile assessment needs finite torque and representable twist",
        )?
        .abs();
        let mut layers = Vec::with_capacity(self.section().layers().len());
        let mut governing_layer = 0;
        let mut maximum_utilization = 0.0;
        for (i, (layer, integral)) in self
            .section()
            .layers()
            .iter()
            .zip(self.section().layer_integrals())
            .enumerate()
        {
            let bending_strain = product(
                curvature,
                layer.outer_radius_m,
                "combined tensile bending strain is unrepresentable",
            )?;
            let torsional_strain = product(
                twist,
                layer.outer_radius_m,
                "combined tensile torsional shear strain is unrepresentable",
            )?;
            let mean = mean_stress(i, bending_strain, torsional_strain)?;
            let amplitude = product(
                integral.young_modulus_pa(),
                bending_strain,
                "combined tensile bending stress is unrepresentable",
            )?;
            let minimum = mean - amplitude;
            let maximum = mean + amplitude;
            if !minimum.is_finite() || !maximum.is_finite() {
                return Err(invalid(
                    "combined tensile normal stress extrema are unrepresentable",
                ));
            }
            let tau = product(
                shear.layer_shear_modulus_pa(i)?,
                torsional_strain,
                "combined tensile torsional shear stress is unrepresentable",
            )?;
            let principal = principal_tension(maximum, tau)?;
            let utilization = checked_utilization(principal, self.layer_limits_pa()[i])?;
            if utilization > maximum_utilization {
                maximum_utilization = utilization;
                governing_layer = i;
            }
            layers.push(LayerPrincipalTensileAssessment {
                minimum_axial_stress_pa: minimum,
                maximum_axial_stress_pa: maximum,
                maximum_torsional_shear_pa: tau,
                principal_tensile_stress_pa: principal,
                tensile_utilization: utilization,
            });
        }
        Ok(ConcentricPrincipalTensileAssessment {
            layers,
            governing_layer,
        })
    }
}

fn principal_tension(normal: f64, shear: f64) -> Result<f64, MaterialStatePointError> {
    if shear == 0.0 {
        return Ok(normal.max(0.0));
    }
    let scale = normal.abs().max(shear);
    let half = 0.5 * (normal.abs() / scale);
    let transverse = shear / scale;
    let radius = fs_math::det::sqrt(half * half + transverse * transverse);
    let tensile = if normal < 0.0 {
        // tau^2/(sqrt((sigma/2)^2+tau^2)+abs(sigma)/2), with no
        // dimensional square, overflowing denominator or subtractive loss.
        shear * (transverse / (radius + half))
    } else {
        scale * (half + radius)
    };
    if !tensile.is_finite() || tensile == 0.0 {
        return Err(invalid(
            "positive combined principal tension is unrepresentable",
        ));
    }
    Ok(tensile)
}

fn product(a: f64, b: f64, quantity: &'static str) -> Result<f64, MaterialStatePointError> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        Err(invalid(quantity))
    } else {
        Ok(value)
    }
}

fn quotient(a: f64, b: f64, quantity: &'static str) -> Result<f64, MaterialStatePointError> {
    let value = a / b;
    if !a.is_finite() || !value.is_finite() || (a != 0.0 && value == 0.0) {
        Err(invalid(quantity))
    } else {
        Ok(value)
    }
}
