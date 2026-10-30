//! Shared thermodynamics of a source-bound, perfectly bonded axial section.
//!
//! Each layer retains its own temperature and material card. With elastic
//! strain `e = strain - integral(alpha dT)`, its Helmholtz energy per reference
//! length is `psi = E A e^2 / 2`. Entropy is `-d psi/dT` and the contribution to
//! internal energy is `psi + T s`; this energy is additional to the source's
//! stress-free caloric chart. The thermal capacity correction is `-T psi_TT`.
//! These are derivatives of ONE potential, shared by strings, supports and
//! other axial reductions, rather than independently authored force/heat laws.
//!
//! This rung admits scalar expansion, scalar reference density, and scalar or
//! globally affine temperature-dependent Young's modulus. It keeps the original
//! reference section and small-strain parallel-strain assumptions. Interlayer
//! slip, radial work, changing mass and phase transitions require other models.
//! Consumers own their actual displacement fields, thermal inventory, positive
//! combined heat-capacity check, integration and atomic publication.
//!
//! Thermodynamic convention: Abeyaratne, *Continuum Mechanics*, Vol. II, §8.4,
//! <https://web.mit.edu/abeyaratne/Volumes/RCA_Vol_II.pdf>.

/// Complete axial, bending, shear and torsional potentials of the same layers.
pub mod rod;

use fs_matdb::{InterpolationPolicy, PropertyValue};

use super::{ConcentricSectionTemperatureLaw, invalid};
use crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS;
use crate::section::{ConcentricElasticSection, thermal::ConcentricThermoelasticSection};
use crate::state_point::{
    DENSITY_PROPERTY, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
    MaterialStatePointError, YOUNG_MODULUS_PROPERTY,
    affine_temperature_slope, pinned_temperature_line,
};

/// Original section/card bindings and admitted constitutive derivatives.
/// This object does not contain mutable thermal or mechanical history.
#[derive(Clone, Debug)]
pub struct ConcentricThermoelasticLaw<'a> {
    temperature: ConcentricSectionTemperatureLaw<'a>,
    modulus_slopes_pa_k: Vec<f64>,
}

/// Mechanics and thermodynamics evaluated at the SAME ordered layer states.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricThermoelasticState {
    /// Original section owner, including source-integrated thermal strains.
    pub section: ConcentricThermoelasticSection,
    /// Per-layer thermodynamics in original spatial order, never material order.
    pub layers: Vec<LayerThermoelasticState>,
}

/// One layer's coefficients and derivatives at fixed total axial strain.
/// Energy, entropy and curvature are per unit reference length. A consumer
/// multiplies by its physical reference length exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerThermoelasticState {
    /// Absolute layer temperature [K].
    pub temperature_k: f64,
    /// Source-evaluated modulus [Pa].
    pub young_modulus_pa: f64,
    /// Exact admitted affine derivative [Pa/K].
    pub young_modulus_derivative_pa_k: f64,
    /// Source scalar linear expansion coefficient [1/K], possibly negative.
    pub alpha_per_k: f64,
    /// Fixed reference density [kg/m3], not an independently evolving mass.
    pub density_kg_m3: f64,
    /// Original annular area [m2].
    pub area_m2: f64,
    /// Original annular second moment [m4].
    pub second_moment_m4: f64,
    /// Total axial strain minus the original integrated thermal strain.
    pub elastic_strain: f64,
    /// `E A e^2 / 2` [J/m].
    pub axial_free_energy_per_reference_length_j_m: f64,
    /// `A (E alpha e - E' e^2 / 2)` [J/(m K)].
    pub entropy_per_reference_length_j_m_k: f64,
    /// `A (E alpha^2 - 2 E' alpha e)` [J/(m K2)].
    pub free_energy_temperature_second_derivative_j_m_k2: f64,
}

impl<'a> ConcentricThermoelasticLaw<'a> {
    /// Admit the original source laws without constructing an effective material.
    /// All knots of an affine modulus claim are checked, including knots not
    /// crossed by the first thermal step. A kink cannot silently create an
    /// entropy jump. Expansion must be explicit, including explicit zero.
    /// The original section supplies the common stress-free reference T.
    pub fn try_new(
        temperature: ConcentricSectionTemperatureLaw<'a>,
    ) -> Result<Self, MaterialStatePointError> {
        let temperatures = temperature.bindings.iter().map(|binding| {
            if binding.axis != "T" {
                return Err(invalid("live section thermodynamics requires the canonical absolute T coordinate"));
            }
            binding.point.axes().get("T").copied()
                .ok_or_else(|| invalid("live section thermodynamics needs a reference temperature"))
        }).collect::<Result<Vec<_>, _>>()?;
        // Reuse source-path admission and the common stress-free-reference
        // requirement of the existing thermal-section owner, even at zero span.
        temperature.thermoelastic_at_temperatures_k(&temperatures)?;
        let mut slopes = Vec::with_capacity(temperature.bindings.len());
        for ((binding, layer), &reference_k) in temperature.bindings.iter()
            .zip(temperature.reference.layers()).zip(&temperatures)
        {
            for name in [DENSITY_PROPERTY, ALPHA] {
                let property = layer.material.property(name)
                    .ok_or_else(|| invalid("thermoelastic layer lacks density or explicit expansion"))?;
                let claim = binding.card.claims().claim(property.answer().receipt.selected)
                    .ok_or_else(|| invalid("thermoelastic layer source claim is absent"))?;
                if !matches!(claim.value, PropertyValue::Scalar { .. }) {
                    return Err(invalid("live section thermodynamics requires scalar reference density and scalar expansion"));
                }
            }
            let modulus = layer.material.property(YOUNG_MODULUS_PROPERTY)
                .ok_or_else(|| invalid("thermoelastic layer lacks Young's modulus"))?;
            let claim = binding.card.claims().claim(modulus.answer().receipt.selected)
                .ok_or_else(|| invalid("thermoelastic elastic source claim is absent"))?;
            let slope = match &claim.value {
                PropertyValue::Scalar { .. } => 0.0,
                PropertyValue::Curve { abscissa, abscissa_dims, knots, .. }
                    if abscissa == "T" && *abscissa_dims == fs_qty::Temperature::DIMS
                        && claim.interpolation == InterpolationPolicy::LinearInside
                        && knots.len() <= MAX_HEAT_CAPACITY_ENTHALPY_KNOTS =>
                {
                    affine_temperature_slope(knots, reference_k, modulus.value_si())
                        .map_err(invalid)?
                }
                _ => return Err(invalid("live section modulus must be scalar or one bounded globally affine temperature curve")),
            };
            slopes.push(slope);
        }
        Ok(Self { temperature, modulus_slopes_pa_k: slopes })
    }

    /// Original fixed geometry, material inventory and complete property bundles.
    #[must_use]
    pub const fn reference(&self) -> &ConcentricElasticSection { self.temperature.reference() }

    /// Original source-pinned evaluator, including auxiliary transport properties.
    #[must_use]
    pub const fn temperature_law(&self) -> &ConcentricSectionTemperatureLaw<'a> { &self.temperature }

    /// Resolve all properties, continuous source paths and thermal derivatives.
    /// Every layer shares the supplied total strain, while its own temperature
    /// determines eigenstrain and stress. Compression remains admissible here;
    /// a taut-string consumer must separately admit positive resultant tension.
    ///
    /// Invalid source support, excessive total/thermal/elastic strain or
    /// unrepresentable arithmetic returns no partially resolved state.
    pub fn at_temperatures_k(
        &self, temperatures_k: &[f64], total_axial_strain: f64, linear_strain_limit: f64,
    ) -> Result<ConcentricThermoelasticState, MaterialStatePointError> {
        let section = self.temperature.thermoelastic_at_temperatures_k(temperatures_k)?;
        let response = section.response_at_axial_strain(total_axial_strain, linear_strain_limit)?;
        let mut layers = Vec::with_capacity(temperatures_k.len());
        for (index, ((layer, integral), axial)) in section.section().layers().iter()
            .zip(section.section().layer_integrals()).zip(response.layers()).enumerate()
        {
            let binding = &self.temperature.bindings[index];
            let point = layer.material.query();
            // Endpoint resolution alone cannot cross an unsupported source cell.
            // Alpha's complete path was already checked by the expansion owner.
            for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY] {
                let property = layer.material.property(name)
                    .ok_or_else(|| invalid("resolved thermoelastic layer lost its original property"))?;
                let claim = binding.card.claims().claim(property.answer().receipt.selected)
                    .ok_or_else(|| invalid("resolved thermoelastic layer lost its source claim"))?;
                pinned_temperature_line(binding.card, &claim.key, "T", &binding.point,
                    &point, property.answer().receipt.selected)?;
            }
            let value = |name| layer.material.property(name).map(|property| property.value_si())
                .ok_or_else(|| invalid("resolved thermoelastic layer is missing a coefficient"));
            let young = value(YOUNG_MODULUS_PROPERTY)?;
            let alpha = value(ALPHA)?;
            let rho = value(DENSITY_PROPERTY)?;
            let slope = self.modulus_slopes_pa_k[index];
            let area = integral.area_m2();
            let strain = axial.elastic_strain();
            // Use the native force evaluation for stored energy, keeping the
            // same force/energy convention as all existing section consumers.
            let energy = finite(0.5 * axial.force_n() * strain)?;
            let entropy = finite(area * (young * alpha * strain - 0.5 * slope * strain * strain))?;
            let curvature = finite(area * (young * alpha * alpha - 2.0 * slope * alpha * strain))?;
            layers.push(LayerThermoelasticState {
                temperature_k: temperatures_k[index], young_modulus_pa: young,
                young_modulus_derivative_pa_k: slope, alpha_per_k: alpha, density_kg_m3: rho,
                area_m2: area, second_moment_m4: integral.second_moment_m4(), elastic_strain: strain,
                axial_free_energy_per_reference_length_j_m: energy,
                entropy_per_reference_length_j_m_k: entropy,
                free_energy_temperature_second_derivative_j_m_k2: curvature,
            });
        }
        Ok(ConcentricThermoelasticState { section, layers })
    }
}

fn finite(value: f64) -> Result<f64, MaterialStatePointError> {
    if value.is_finite() { Ok(value) }
    else { Err(invalid("unrepresentable section Helmholtz derivative")) }
}

#[cfg(test)]
mod tests;
