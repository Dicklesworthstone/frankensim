//! Thermal entropy on the already accepted finite body/gas endpoint.

use super::{Error, FiniteGasThermalEnvironment, gas_error};
use fs_material::{
    phase::{EquilibriumPhaseState, UniformEnthalpyStepInput},
    thermochemical_gas::energy::{GasEnergyError, IsochoricGasParcel},
};

/// Entropy of one accepted implicit thermal proposal, in J/K.
///
/// Every two-element array is ordered `[body, gas]`. Storage uses relative
/// entropy from the original material owners: no absolute entropy reference,
/// constant heat capacity, or separate latent-heat credit is introduced.
/// The body retains its chart's pressure/composition and energy convention;
/// the gas has fixed volume and composition, including its actual EOS pressure
/// change. Other mechanical, chemical and radiation-field entropy is not inferred.
///
/// Storage equals reservoir input + the two transfer productions + unclassified
/// heat-at-endpoint + temporal defects + source and algebraic residuals, up to
/// `balance_residual_j_k`. All terms retain their computed signs. This identity
/// is not a temporal-accuracy certificate or a global second-law proof.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasEntropyReport {
    /// Exact `[m_body integral(dh/T), m_gas delta_s(T,p)]` on accepted states.
    pub storage_change_j_k: [f64; 2],
    /// Q*(1/T_gas - 1/T_body), using the actual reciprocal convective heat Q.
    /// Physical finite-temperature transfer, separate from all solver defects.
    pub convection_production_j_k: f64,
    /// Q_rad/T_enclosure at the actual independent enclosure temperature.
    /// This is a THERMAL RESERVOIR net-heat account, not photon entropy flux.
    pub radiation_reservoir_input_j_k: f64,
    /// Q_rad*(1/T_body - 1/T_enclosure), the associated thermal-reservoir
    /// transfer production. No blackbody photon-field entropy is asserted.
    pub radiation_transfer_production_j_k: f64,
    /// Prescribed body/gas heat divided by its receiving endpoint temperature.
    /// Source owners must classify heat, work conversion and boundary entropy;
    /// these contributions are neither entropy-free sources nor certified
    /// physical production. Radiation and reciprocal convection are excluded.
    pub unclassified_heat_at_endpoint_j_k: [f64; 2],
    /// Exact storage entropy minus measured storage energy/T_endpoint for each
    /// owner. These backward-Euler temporal defects are not physical production.
    pub temporal_defect_j_k: [f64; 2],
    /// Interface/source-matching energy residual divided by T_body,end.
    pub source_residual_j_k: f64,
    /// Original owner energy residuals divided by their endpoint temperatures:
    /// `[(solid_residual-source_residual)/T_body, (gas_storage-gas_heat-Q)/T_gas]`.
    pub algebraic_residual_j_k: [f64; 2],
    /// Compensated storage-minus-account difference; never converted into heat.
    pub balance_residual_j_k: f64,
}

pub(super) fn evaluate(
    environment: &FiniteGasThermalEnvironment<'_>,
    input: UniformEnthalpyStepInput<'_>,
    next_gas: &IsochoricGasParcel,
    body: EquilibriumPhaseState,
    [body_storage, gas_storage]: [f64; 2],
    [heat_j, radiation_heat_j]: [f64; 2],
    [source_residual_j, solid_residual_j]: [f64; 2],
) -> Result<FiniteGasEntropyReport, Error> {
    let body_k = body.temperature_k();
    let gas_k = next_gas.state().state().temperature;
    let enclosure_k = environment.body.radiation_temperature_k;
    let storage_change_j_k = [
        input.mass_kg
            * input.curve.specific_entropy_change_j_kg_k(
                input.initial.specific_enthalpy_j_kg(),
                body.specific_enthalpy_j_kg(),
            )?,
        environment.gas.mass_kg()
            * environment
                .gas
                .model()
                .specific_entropy_change_j_kg_k(environment.gas.point(), next_gas.point())
                .map_err(GasEnergyError::from)
                .map_err(gas_error)?,
    ];
    let mut report = FiniteGasEntropyReport {
        storage_change_j_k,
        // Differences preserve the small force near equilibrium; the ratio is
        // bounded in either direction without a temperature product.
        convection_production_j_k: (heat_j / body_k.min(gas_k))
            * ((body_k - gas_k) / body_k.max(gas_k)),
        radiation_reservoir_input_j_k: radiation_heat_j / enclosure_k,
        radiation_transfer_production_j_k: (radiation_heat_j / body_k.min(enclosure_k))
            * ((enclosure_k - body_k) / body_k.max(enclosure_k)),
        unclassified_heat_at_endpoint_j_k: [
            input.internal_heat_j / body_k,
            environment.gas_heat_j / gas_k,
        ],
        temporal_defect_j_k: [
            storage_change_j_k[0] - body_storage / body_k,
            storage_change_j_k[1] - gas_storage / gas_k,
        ],
        source_residual_j_k: source_residual_j / body_k,
        algebraic_residual_j_k: [
            (solid_residual_j - source_residual_j) / body_k,
            sum([gas_storage, -environment.gas_heat_j, -heat_j])? / gas_k,
        ],
        balance_residual_j_k: 0.0,
    };
    report.balance_residual_j_k = sum([
        report.storage_change_j_k[0],
        report.storage_change_j_k[1],
        -report.convection_production_j_k,
        -report.radiation_reservoir_input_j_k,
        -report.radiation_transfer_production_j_k,
        -report.unclassified_heat_at_endpoint_j_k[0],
        -report.unclassified_heat_at_endpoint_j_k[1],
        -report.temporal_defect_j_k[0],
        -report.temporal_defect_j_k[1],
        -report.source_residual_j_k,
        -report.algebraic_residual_j_k[0],
        -report.algebraic_residual_j_k[1],
    ])?;
    Ok(report)
}

// Fixed-order Neumaier summation retains small residuals beside cancelling
// body/gas storage. Validation applies to every exposed term before publication.
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, Error> {
    let mut value: f64 = 0.0;
    let mut correction = 0.0;
    for term in values {
        let next = value + term;
        if !term.is_finite() || !next.is_finite() {
            return Err(Error::Invalid("finite gas/body entropy is unrepresentable"));
        }
        correction += if value.abs() >= term.abs() {
            (value - next) + term
        } else {
            (term - next) + value
        };
        value = next;
    }
    let total = value + correction;
    if !total.is_finite() {
        return Err(Error::Invalid("finite gas/body entropy is unrepresentable"));
    }
    Ok(total)
}
