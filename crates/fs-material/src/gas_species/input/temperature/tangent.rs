//! Analytic temperature response of the original retained gas transport.
//!
//! Differentiate the existing Sutherland/USSA/Eucken and Wilke/WMS equations
//! at fixed pressure and composition. Values still come from the original
//! gas evaluator. NASA-9 source gases retain their actual selected transport
//! claims and shared Wilke owner. This adds no source or accuracy authority.

mod thermochemical;

use super::{GasRecord, GasSpec, GasState, IsobaricGasLaw, MaterialError, Source, Temperature};
use crate::fluid_transport::FluidTransport;
use crate::gas::{ConductivityModel, R_USSA_1976};
use crate::gas_mixture::interaction;
use crate::gas_species::GasSpecies;
use fs_math::det;

/// One original admitted gas record and its local isobaric transport response.
/// All derivatives hold the original mole fractions fixed, including water.
/// They describe the declared constitutive models, not measured uncertainty.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricGasTransportTangent {
    pub gas: GasRecord,
    /// d rho / d T at fixed pressure/composition [kg/(m³ K)].
    pub density_derivative_kg_m3_k: f64,
    /// d Cp / d T [J/(kg K²)]; zero only for calorically perfect models.
    pub specific_heat_derivative_j_kg_k2: f64,
    /// d k / d T [W/(m K²)].
    pub conductivity_derivative_w_m_k2: f64,
    /// d mu / d T [Pa s/K].
    pub dynamic_viscosity_derivative_pa_s_k: f64,
}

impl IsobaricGasLaw {
    /// Evaluate the SAME gas and its analytic temperature transport tangent.
    ///
    /// Legacy dry/moist air and named species retain their original gas,
    /// mixture, temperature, pressure and phase admissions. Ideal-gas density
    /// gives `rho'=-rho/T`; Cp is frozen by these engineering models. Sutherland
    /// viscosity, USSA/Eucken conductivity and every temperature-dependent
    /// Wilke interaction denominator are differentiated without finite
    /// differences. The additional derivative work uses stack storage for at
    /// most six species; it does not clone/query a material database.
    ///
    /// Thermochemical gases differentiate original NASA-9 Cp and selected
    /// constant or linear-temperature transport claims. Their original Wilke
    /// evaluation supplies viscosity and WMS denominator derivatives. No source
    /// is reselected, numerically perturbed or replaced by an engineering fit.
    /// At a domain endpoint this is the inward derivative; interior source
    /// jumps/corners and discrete-only transport refuse. Unsupported derivative
    /// forms do not change value-only [`Self::at_temperature`] admission.
    pub fn transport_temperature_tangent(
        &self,
        temperature: Temperature,
    ) -> Result<IsobaricGasTransportTangent, MaterialError> {
        let gas = self.at_temperature(temperature)?;
        gas.transport_state()?;
        let t = temperature.value();
        let pressure = gas.state().pressure;
        let [specific_heat, viscosity, conductivity] = match &self.source {
            Source::MoistAir { .. } if self.reference.state().water_mole_fraction == 0.0 => {
                let logarithmic = component_log_derivatives(GasSpec::dry_air_ussa1976(), t);
                [
                    0.0,
                    gas.state().dynamic_viscosity * logarithmic[0],
                    gas.state().thermal_conductivity * logarithmic[1],
                ]
            }
            Source::MoistAir { composition } => {
                let [viscosity, conductivity] = mixture_derivatives(
                    &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
                    composition.fractions(),
                    t,
                    pressure,
                )?;
                [0.0, viscosity, conductivity]
            }
            Source::Mixture {
                species,
                composition,
            } => {
                let mut specs = [GasSpec::dry_air_ussa1976(); GasSpecies::ALL.len()];
                if species.len() > specs.len() {
                    return Err(invalid(
                        "bounded original gas species required for transport tangent",
                    ));
                }
                for (target, species) in specs.iter_mut().zip(species) {
                    *target = species.spec();
                }
                let [viscosity, conductivity] = mixture_derivatives(
                    &specs[..species.len()],
                    composition.fractions(),
                    t,
                    pressure,
                )?;
                [0.0, viscosity, conductivity]
            }
            Source::Thermochemical { model, .. } => {
                let GasRecord::Thermochemical(snapshot) = &gas else {
                    return Err(invalid(
                        "thermochemical derivative lost its original snapshot",
                    ));
                };
                thermochemical::derivatives(model, snapshot)?
            }
        };
        let density = -gas.state().density / t;
        if ![density, specific_heat, viscosity, conductivity]
            .into_iter()
            .all(f64::is_finite)
            || density == 0.0
        {
            return Err(invalid(
                "gas transport temperature tangent is unrepresentable",
            ));
        }
        Ok(IsobaricGasTransportTangent {
            gas,
            density_derivative_kg_m3_k: density,
            specific_heat_derivative_j_kg_k2: specific_heat,
            conductivity_derivative_w_m_k2: conductivity,
            dynamic_viscosity_derivative_pa_s_k: viscosity,
        })
    }
}

fn component_log_derivatives(spec: GasSpec, t: f64) -> [f64; 2] {
    // mu=beta*T^(3/2)/(T+S). The coefficients are the original GasSpec.
    let viscosity = 1.5 / t - 1.0 / (t + spec.sutherland_s);
    let conductivity = match spec.conductivity {
        ConductivityModel::Eucken => viscosity,
        ConductivityModel::Ussa1976AirFit => {
            // k=C*T^(3/2)/(T+v), v=245.4*exp(-12*ln(10)/T).
            let v = 245.4 * det::exp((-12.0 / t) * core::f64::consts::LN_10);
            1.5 / t - (1.0 + v * (12.0 * core::f64::consts::LN_10 / (t * t))) / (t + v)
        }
    };
    [viscosity, conductivity]
}

fn mixture_derivatives(
    specs: &[GasSpec],
    fractions: &[f64],
    t: f64,
    pressure: f64,
) -> Result<[f64; 2], MaterialError> {
    if specs.is_empty() || specs.len() > GasSpecies::ALL.len() || specs.len() != fractions.len() {
        return Err(invalid(
            "transport tangent needs the original bounded component fractions",
        ));
    }
    let mut samples = [None; GasSpecies::ALL.len()];
    for (slot, spec) in samples.iter_mut().zip(specs) {
        let state = GasState::try_new(spec, t, pressure)?;
        // Match the molar masses reconstructed by the original mixture owner.
        let mass = R_USSA_1976 / state.specific_gas_constant;
        *slot = Some((state, mass, component_log_derivatives(*spec, t)));
    }
    let mut derivative = [0.0; 2];
    for i in 0..specs.len() {
        let (left, mi, logarithmic) = samples[i].expect("original component evaluated");
        if specs.len() == 1 {
            return Ok([
                left.dynamic_viscosity * logarithmic[0],
                left.thermal_conductivity * logarithmic[1],
            ]);
        }
        let mut denominator = 0.0;
        let mut denominator_derivative = 0.0;
        for j in 0..specs.len() {
            let (right, mj, other) = samples[j].expect("original component evaluated");
            let (phi, phi_derivative) = if i == j {
                (1.0, 0.0)
            } else {
                let phi = interaction(&left, &right, mi, mj);
                let z = (left.dynamic_viscosity / right.dynamic_viscosity).sqrt()
                    * (mj / mi).sqrt().sqrt();
                // phi=(1+z)^2/B and z'/z=(log(mu_i)'-log(mu_j)')/2.
                (phi, phi * (z / (1.0 + z)) * (logarithmic[0] - other[0]))
            };
            denominator += fractions[j] * phi;
            denominator_derivative += fractions[j] * phi_derivative;
        }
        if !denominator.is_finite() || denominator <= 0.0 || !denominator_derivative.is_finite() {
            return Err(invalid(
                "gas mixture transport tangent denominator is unrepresentable",
            ));
        }
        let response = denominator_derivative / denominator;
        derivative[0] +=
            (fractions[i] * left.dynamic_viscosity / denominator) * (logarithmic[0] - response);
        derivative[1] +=
            (fractions[i] * left.thermal_conductivity / denominator) * (logarithmic[1] - response);
    }
    Ok(derivative)
}

fn invalid(what: &'static str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

#[cfg(test)]
mod tests;
