//! Thermal derivatives of the original concentric Cosserat section energy.
//!
//! The axial owner supplies integrated expansion and its Helmholtz derivatives.
//! Added terms are E I |k_b|²/2 and G (kappa A |gamma_s|² + J k_t²)/2.
//! Affine E/G make these extra terms linear in temperature at fixed total
//! strain, so their entropy is minus their temperature slope and their fixed-
//! strain heat-capacity correction is zero. No thermal/motion integrator lives
//! here; geometry, source pins and physical shear correction remain unchanged.

use super::{ConcentricThermoelasticLaw, ConcentricThermoelasticState, finite, invalid};
use crate::section::ConcentricShearSection;
use crate::section_temperature::ConcentricSectionTemperatureLaw;
use crate::state_point::{
    MaterialStatePointError, affine_temperature_slope, pinned_temperature_line,
};
use fs_matdb::{InterpolationPolicy, PropertyValue};

/// Same concentric layer sources, with an admitted affine shear response.
#[derive(Clone, Debug)]
pub struct ConcentricRodThermoelasticLaw<'a> {
    axial: ConcentricThermoelasticLaw<'a>,
    shear_property: String,
    shear_correction: f64,
    shear_slopes_pa_k: Vec<f64>,
}

/// Integrated quantities per reference length, additional to stress-free heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodLayerThermoelasticState {
    pub temperature_k: f64,
    pub shear_modulus_pa: f64,
    pub shear_modulus_derivative_pa_k: f64,
    pub free_energy_j_m: f64,
    pub entropy_j_m_k: f64,
    pub internal_energy_j_m: f64,
    /// Can be negative; the consumer must admit caloric plus mechanical capacity.
    pub heat_capacity_j_m_k: f64,
    pub maximum_absolute_strain: f64,
}

/// Complete original source sections and per-layer thermodynamic observations.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricRodThermoelasticState {
    pub axial: ConcentricThermoelasticState,
    pub section: ConcentricShearSection,
    pub layers: Vec<RodLayerThermoelasticState>,
}

impl<'a> ConcentricRodThermoelasticLaw<'a> {
    /// Scalar rho/alpha and scalar or globally affine E/G are the admitted law.
    /// The existing section owner validates the original named shear coefficient
    /// and correction; every selected shear knot is checked before acceptance.
    pub fn try_new(
        temperature: ConcentricSectionTemperatureLaw<'a>,
        shear_property: &str,
        shear_correction: f64,
    ) -> Result<Self, MaterialStatePointError> {
        temperature
            .reference()
            .clone()
            .with_shear_modulus(shear_property, shear_correction)?;
        let axial = ConcentricThermoelasticLaw::try_new(temperature)?;
        let mut shear_slopes_pa_k = Vec::with_capacity(axial.reference().layers().len());
        for (binding, layer) in axial
            .temperature
            .bindings
            .iter()
            .zip(axial.reference().layers())
        {
            let value = layer
                .material
                .property(shear_property)
                .ok_or_else(|| invalid("rod thermodynamics requires its original shear modulus"))?;
            let claim = binding
                .card
                .claims()
                .claim(value.answer().receipt.selected)
                .ok_or_else(|| invalid("rod shear source claim is absent"))?;
            let slope = match &claim.value {
                PropertyValue::Scalar { .. } => 0.0,
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    knots,
                    ..
                } if abscissa == "T"
                    && *abscissa_dims == fs_qty::Temperature::DIMS
                    && claim.interpolation == InterpolationPolicy::LinearInside
                    && knots.len() <= crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS =>
                {
                    affine_temperature_slope(knots, binding.point.axes()["T"], value.value_si())
                        .map_err(invalid)?
                }
                _ => {
                    return Err(invalid(
                        "rod shear modulus must be scalar or one globally affine T curve",
                    ));
                }
            };
            shear_slopes_pa_k.push(slope);
        }
        Ok(Self {
            axial,
            shear_property: shear_property.into(),
            shear_correction,
            shear_slopes_pa_k,
        })
    }

    pub const fn axial_law(&self) -> &ConcentricThermoelasticLaw<'a> {
        &self.axial
    }

    /// Cosserat components follow the native rod: axial/shears in gamma and
    /// twist/bend/bend in curvature. Intrinsic geometry has already been removed
    /// by that owner; thermal axial strain is applied here exactly once.
    pub fn at_strain(
        &self,
        temperatures_k: &[f64],
        gamma: [f64; 3],
        curvature: [f64; 3],
        linear_strain_limit: f64,
    ) -> Result<ConcentricRodThermoelasticState, MaterialStatePointError> {
        if gamma.iter().chain(&curvature).any(|x| !x.is_finite()) {
            return Err(invalid("finite original Cosserat strains required"));
        }
        let axial = self
            .axial
            .at_temperatures_k(temperatures_k, gamma[0], linear_strain_limit)?;
        let section = axial
            .section
            .section()
            .clone()
            .with_shear_modulus(&self.shear_property, self.shear_correction)?;
        let bend_square = finite(curvature[1].powi(2) + curvature[2].powi(2))?;
        let shear_square = finite(gamma[1].powi(2) + gamma[2].powi(2))?;
        let mut layers = Vec::with_capacity(axial.layers.len());
        for (i, layer) in axial.layers.iter().enumerate() {
            let original = &section.elastic().layers()[i];
            let binding = &self.axial.temperature.bindings[i];
            let value = original
                .material
                .property(&self.shear_property)
                .ok_or_else(|| invalid("resolved rod lost its original shear modulus"))?;
            let claim = binding
                .card
                .claims()
                .claim(value.answer().receipt.selected)
                .ok_or_else(|| invalid("resolved rod lost its selected shear source"))?;
            pinned_temperature_line(
                binding.card,
                &claim.key,
                "T",
                &binding.point,
                &original.material.query(),
                value.answer().receipt.selected,
            )?;
            let maximum_absolute_strain = finite(
                (gamma[0].abs().max(layer.elastic_strain.abs())
                    + original.outer_radius_m * bend_square.sqrt())
                .max(shear_square.sqrt() + original.outer_radius_m * curvature[0].abs()),
            )?;
            if maximum_absolute_strain > linear_strain_limit {
                return Err(invalid(
                    "rod fibre/shear strain exceeds its thermodynamic envelope",
                ));
            }
            let bending_geometry = finite(layer.second_moment_m4 * bend_square)?;
            let shear_geometry = finite(
                self.shear_correction * layer.area_m2 * shear_square
                    + 2.0 * layer.second_moment_m4 * curvature[0].powi(2),
            )?;
            let free_energy_j_m = finite(
                layer.axial_free_energy_per_reference_length_j_m
                    + 0.5
                        * (layer.young_modulus_pa * bending_geometry
                            + value.value_si() * shear_geometry),
            )?;
            let entropy_j_m_k = finite(
                layer.entropy_per_reference_length_j_m_k
                    - 0.5
                        * (layer.young_modulus_derivative_pa_k * bending_geometry
                            + self.shear_slopes_pa_k[i] * shear_geometry),
            )?;
            layers.push(RodLayerThermoelasticState {
                temperature_k: layer.temperature_k,
                shear_modulus_pa: value.value_si(),
                shear_modulus_derivative_pa_k: self.shear_slopes_pa_k[i],
                free_energy_j_m,
                entropy_j_m_k,
                internal_energy_j_m: finite(free_energy_j_m + layer.temperature_k * entropy_j_m_k)?,
                heat_capacity_j_m_k: finite(
                    -layer.temperature_k * layer.free_energy_temperature_second_derivative_j_m_k2,
                )?,
                maximum_absolute_strain,
            });
        }
        Ok(ConcentricRodThermoelasticState {
            axial,
            section,
            layers,
        })
    }
}
