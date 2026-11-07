//! Total normal stress of the original axial and bending Kelvin--Voigt layers.
//!
//! Each perfectly bonded layer retains its own E, axial eta and bending eta:
//! `sigma_mean_i = E_i*epsilon + eta_ai*epsilon_dot`, and the two bending
//! coefficients are `E_i*kappa + eta_bi*kappa_dot`. Their vector sum precedes
//! the annular maximum; opposing elastic and viscous stresses can cancel.
//! Summed section viscosities never redistribute stress through an averaged
//! eta/E ratio. Preparation retains source identities and caches coefficients;
//! a single-layer assessment allocates nothing and performs no source queries.

use super::{AxialTensileCapacity, ConcentricTensileAssessment, LayerTensileAssessment};
use crate::section::axial::ConcentricAxialViscosity;
use crate::section::{ConcentricKelvinVoigtSection, MaterialStatePointError, invalid};

/// Original fixed-state coefficients of one layer, in spatial layer order.
/// Constructed only by the source-bound capacity; these are not fitted averages.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KelvinVoigtTensileLayer {
    young_modulus_pa: f64,
    axial_viscosity_pa_s: f64,
    bending_viscosity_pa_s: f64,
}

impl KelvinVoigtTensileLayer {
    /// Original layer Young modulus [Pa].
    #[must_use]
    pub const fn young_modulus_pa(&self) -> f64 {
        self.young_modulus_pa
    }
    /// Independently selected axial viscosity [Pa s].
    #[must_use]
    pub const fn axial_viscosity_pa_s(&self) -> f64 {
        self.axial_viscosity_pa_s
    }
    /// Independently selected bending viscosity [Pa s].
    #[must_use]
    pub const fn bending_viscosity_pa_s(&self) -> f64 {
        self.bending_viscosity_pa_s
    }
}

/// Total elastic-plus-viscous normal stress against each original tensile limit.
/// The retained elastic capacity contains every geometry and material receipt,
/// including the viscosity sources selected here. No strength is inferred from
/// the viscosities, and no plasticity, rupture or hereditary memory is advanced.
#[derive(Clone, Debug, PartialEq)]
pub struct KelvinVoigtTensileCapacity {
    elastic_capacity: AxialTensileCapacity,
    axial_viscosity_property: String,
    bending_viscosity_property: String,
    layers: Vec<KelvinVoigtTensileLayer>,
    axial_viscous_rigidity_n_s: f64,
    bending_viscous_stiffness_n_m2_s: f64,
    omega_band_rad_s: (f64, f64),
}

impl AxialTensileCapacity {
    /// Prepare total normal-stress assessment from both original viscosity owners.
    ///
    /// Axial and bending owners must retain this exact elastic section, layer
    /// order, geometry, source selections and held material states. Their
    /// independent viscosity selectors and complete receipts remain available
    /// through the retained capacity. Explicit zero viscosity is supported;
    /// missing data is never interpreted as zero.
    ///
    /// The returned source band is the intersection of both owners' bands.
    /// The dynamics consumer must admit its actual strain/rate history against
    /// that band and retain the original small-strain section assumptions.
    /// An instantaneous rate alone cannot establish frequency applicability.
    ///
    /// # Errors
    /// Foreign or stale source sections and disjoint source bands refuse.
    pub fn with_kelvin_voigt(
        self,
        axial: &ConcentricAxialViscosity,
        bending: &ConcentricKelvinVoigtSection,
    ) -> Result<KelvinVoigtTensileCapacity, MaterialStatePointError> {
        if axial.elastic() != self.section() || bending.elastic() != self.section() {
            return Err(invalid(
                "viscous tensile assessment requires the exact original axial and bending sections and material sources",
            ));
        }
        let axial_band = axial.omega_band_rad_s();
        let bending_band = bending.omega_band_rad_s();
        let omega_band_rad_s = (
            axial_band.0.max(bending_band.0),
            axial_band.1.min(bending_band.1),
        );
        if omega_band_rad_s.0 > omega_band_rad_s.1 {
            return Err(invalid(
                "axial and bending tensile-assessment frequency bands do not intersect",
            ));
        }
        let mut layers = Vec::with_capacity(self.section().layers().len());
        for (layer, geometry) in self
            .section()
            .layers()
            .iter()
            .zip(self.section().layer_integrals())
        {
            let axial_source = layer
                .material
                .property(axial.viscosity_property())
                .ok_or_else(|| invalid("original axial viscosity source is missing"))?;
            let bending_source = layer
                .material
                .property(bending.viscosity_property())
                .ok_or_else(|| invalid("original bending viscosity source is missing"))?;
            layers.push(KelvinVoigtTensileLayer {
                young_modulus_pa: geometry.young_modulus_pa(),
                axial_viscosity_pa_s: axial_source.value_si(),
                bending_viscosity_pa_s: bending_source.value_si(),
            });
        }
        Ok(KelvinVoigtTensileCapacity {
            elastic_capacity: self,
            axial_viscosity_property: axial.viscosity_property().into(),
            bending_viscosity_property: bending.viscosity_property().into(),
            layers,
            axial_viscous_rigidity_n_s: axial.viscous_rigidity_n_s(),
            bending_viscous_stiffness_n_m2_s: bending.viscous_stiffness_n_m2_s(),
            omega_band_rad_s,
        })
    }
}

impl KelvinVoigtTensileCapacity {
    /// Original geometry, tensile limits and complete material source receipts.
    #[must_use]
    pub const fn elastic_capacity(&self) -> &AxialTensileCapacity {
        &self.elastic_capacity
    }
    /// Original per-layer E and both independently selected viscosities.
    #[must_use]
    pub fn layers(&self) -> &[KelvinVoigtTensileLayer] {
        &self.layers
    }
    /// Exact selector of the original axial viscosity claims.
    #[must_use]
    pub fn axial_viscosity_property(&self) -> &str {
        &self.axial_viscosity_property
    }
    /// Exact selector of the original bending viscosity claims.
    #[must_use]
    pub fn bending_viscosity_property(&self) -> &str {
        &self.bending_viscosity_property
    }
    /// Original sum of eta_axial_i*A_i [N s], without a second reduction.
    #[must_use]
    pub const fn axial_viscous_rigidity_n_s(&self) -> f64 {
        self.axial_viscous_rigidity_n_s
    }
    /// Original sum of eta_bending_i*I_i [N m2 s], without a second reduction.
    #[must_use]
    pub const fn bending_viscous_stiffness_n_m2_s(&self) -> f64 {
        self.bending_viscous_stiffness_n_m2_s
    }
    /// Common fixed-state source applicability [rad/s]; actual history admission
    /// remains the dynamics owner's responsibility.
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }

    /// Assess signed common axial engineering strain, its rate [1/s], and two
    /// centroidal bending curvatures [1/m] and curvature rates [1/(m s)].
    ///
    /// Every layer uses its own E, axial eta and bending eta. Elastic and
    /// viscous bending vectors are added before taking their magnitude, giving
    /// the exact normal-stress extrema on each complete concentric annulus.
    /// The supplied strain is mechanical strain; layer-specific thermal
    /// eigenstrains, residual stress and hereditary stress are not supplied by
    /// this method. This observation changes neither state nor force.
    ///
    /// # Errors
    /// Nonfinite kinematics and unrepresentable stresses or utilization refuse.
    pub fn assess_strain_rates(
        &self,
        axial_strain: f64,
        axial_strain_rate: f64,
        curvature: [f64; 2],
        curvature_rate: [f64; 2],
    ) -> Result<ConcentricTensileAssessment, MaterialStatePointError> {
        finite_kinematics(axial_strain, axial_strain_rate, curvature, curvature_rate)?;
        let mut layers = Vec::with_capacity(self.layers.len());
        let mut governing_layer = 0;
        let mut maximum_utilization = 0.0;
        for (index, coefficients) in self.layers.iter().enumerate() {
            let assessment = self.assess_layer(
                index,
                coefficients,
                axial_strain,
                axial_strain_rate,
                curvature,
                curvature_rate,
            )?;
            if assessment.tensile_utilization > maximum_utilization {
                maximum_utilization = assessment.tensile_utilization;
                governing_layer = index;
            }
            layers.push(assessment);
        }
        Ok(ConcentricTensileAssessment {
            layers,
            governing_layer,
        })
    }

    /// Allocation-free assessment of one original layer with the same signed
    /// kinematics as [`Self::assess_strain_rates`].
    ///
    /// # Errors
    /// An invalid layer index, nonfinite kinematics or unrepresentable result refuses.
    pub fn assess_layer_strain_rates(
        &self,
        layer: usize,
        axial_strain: f64,
        axial_strain_rate: f64,
        curvature: [f64; 2],
        curvature_rate: [f64; 2],
    ) -> Result<LayerTensileAssessment, MaterialStatePointError> {
        finite_kinematics(axial_strain, axial_strain_rate, curvature, curvature_rate)?;
        let coefficients = self.layers.get(layer).ok_or_else(|| {
            invalid("viscous tensile assessment layer is outside the original section")
        })?;
        self.assess_layer(
            layer,
            coefficients,
            axial_strain,
            axial_strain_rate,
            curvature,
            curvature_rate,
        )
    }

    fn assess_layer(
        &self,
        layer: usize,
        coefficients: &KelvinVoigtTensileLayer,
        axial_strain: f64,
        axial_strain_rate: f64,
        curvature: [f64; 2],
        curvature_rate: [f64; 2],
    ) -> Result<LayerTensileAssessment, MaterialStatePointError> {
        let radius = self.elastic_capacity.section().layers()[layer].outer_radius_m;
        let mean = product(coefficients.young_modulus_pa, axial_strain)?
            + product(coefficients.axial_viscosity_pa_s, axial_strain_rate)?;
        let mut bending = [0.0; 2];
        for axis in 0..2 {
            let strain = product(radius, curvature[axis])?;
            let viscous = if coefficients.bending_viscosity_pa_s == 0.0 {
                0.0
            } else {
                product(
                    coefficients.bending_viscosity_pa_s,
                    product(radius, curvature_rate[axis])?,
                )?
            };
            bending[axis] = product(coefficients.young_modulus_pa, strain)? + viscous;
        }
        let amplitude = fs_math::det::hypot(bending[0], bending[1]);
        self.elastic_capacity
            .assess_layer_normal_stress(layer, mean, amplitude)
    }
}

fn finite_kinematics(
    strain: f64,
    strain_rate: f64,
    curvature: [f64; 2],
    curvature_rate: [f64; 2],
) -> Result<(), MaterialStatePointError> {
    if [strain, strain_rate]
        .iter()
        .chain(&curvature)
        .chain(&curvature_rate)
        .any(|value| !value.is_finite())
    {
        return Err(invalid(
            "viscous tensile assessment requires finite signed strains and rates",
        ));
    }
    Ok(())
}

fn product(coefficient: f64, value: f64) -> Result<f64, MaterialStatePointError> {
    let result = coefficient * value;
    if !result.is_finite() || coefficient != 0.0 && value != 0.0 && result == 0.0 {
        Err(invalid(
            "viscous tensile stress or strain product is unrepresentable",
        ))
    } else {
        Ok(result)
    }
}
