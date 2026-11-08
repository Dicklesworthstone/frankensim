//! Exact fixed-temperature-bath exchange in the original energy coordinate.

use super::{CaloricConductionLayer, OriginalConduction};
use super::super::capacity::{CapacityCurve, CapacityExchange, CapacityKnot};
use fs_matdb::MaterialAnswer;
use fs_thermal::enthalpy::{DepositionError, ThermalStorageBasis};

/// One immutable exact exchange proposal and its represented native storage.
///
/// Heat, mean temperature and tangents use the original capacity integral and
/// retained relative temperature travel. They survive rounding of absolute
/// temperature or the stored energy. The explicit residual reports that storage
/// rounding; it is not source uncertainty or a heat-transfer-model error bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricHeatExchange {
    /// Original h, u, or declared sensible-storage basis; never relabeled.
    pub storage_basis: ThermalStorageBasis,
    pub initial_specific_energy_j_kg: f64,
    pub final_specific_energy_j_kg: f64,
    /// Original caloric inverse of the returned final native energy [K].
    /// This can differ by roundoff from the exchange's relative endpoint.
    pub temperature_k: f64,
    /// Conductance-weighted mean of the exact material trajectory [K].
    pub mean_temperature_k: f64,
    /// Signed integrated heat into the body [J/kg], before native-energy rounding.
    pub specific_heat_j_kg: f64,
    /// Returned final native energy minus the exact supplied initial energy.
    pub represented_energy_change_j_kg: f64,
    /// Integrated heat minus the represented native-energy change [J/kg].
    pub heat_residual_j_kg: f64,
    /// Exact exchange derivative d(e_out)/d(e_in), at fixed bath and exposure.
    pub inlet_energy_gain: f64,
    /// One minus inlet_energy_gain, retained directly for weak exchange.
    pub inlet_energy_loss: f64,
    /// d(e_out)/d(T_bath) at fixed inlet energy and exposure [J/(kg K)].
    pub reservoir_energy_derivative_j_kg_k: f64,
    /// d(mean_temperature)/d(T_bath) at fixed inlet energy and exposure.
    pub mean_reservoir_derivative: f64,
    /// T_bath minus the exchange endpoint, before absolute-temperature rounding.
    pub remaining_temperature_difference_k: f64,
    /// Exposure times d(mean_temperature)/d(exposure) [K].
    /// Equals the endpoint-minus-mean gap before absolute-temperature rounding.
    pub mean_log_exposure_derivative_k: f64,
}

// Original source pairs already retain their positive affine capacity. Borrow
// them directly: no copied/resampled Cp table and no runtime material-card query.
impl CapacityKnot for (f64, MaterialAnswer) {
    fn temperature(&self) -> f64 {
        self.0
    }
    fn capacity(&self) -> f64 {
        self.1.evidence.value.value
    }
}

impl CaloricConductionLayer<'_, '_> {
    /// Exchange with one fixed-temperature bath using every crossed original
    /// affine capacity segment. The caller declares a lumped uniform-temperature
    /// body and a constant conductance G over the window; exposure is G*dt/m
    /// [J/(kg K)]. Geometry, source uncertainty and the validity of that transfer
    /// model remain with the caller.
    ///
    /// Solves integral C_specific(T)/(T_bath-T) dT = exposure without freezing
    /// capacity. Cp retains enthalpy on its original isobar, Cv retains internal
    /// energy at its original density, and amount/volume storage retains its
    /// declared sensible basis. Free expansion retains the same Cp enthalpy;
    /// pressure work is not added a second time. Amount and volumetric capacities
    /// use the original reference mass or density exactly once.
    ///
    /// The bath may lie outside source support; the complete material path must
    /// remain inside it. Successful queries allocate nothing and query no card.
    /// Work is bounded by crossed source knots and the shared segment inversion.
    /// No timeline, cell, mass, source receipt or material object is mutated.
    ///
    /// Native energy advances by the relative heat even when absolute T cannot
    /// move. A supported source endpoint uses its canonical native energy bound.
    /// If the stored coordinate itself cannot resolve all the heat, the residual
    /// reports the difference. The final coordinate must pass its ORIGINAL
    /// caloric inverse; unsupported paths and unrepresentable arithmetic refuse.
    pub fn linear_heat_exchange(
        &self,
        initial_specific_energy_j_kg: f64,
        reservoir_temperature_k: f64,
        exposure_j_kg_k: f64,
    ) -> Result<CaloricHeatExchange, DepositionError> {
        if !reservoir_temperature_k.is_finite()
            || reservoir_temperature_k <= 0.0
            || !exposure_j_kg_k.is_finite()
            || exposure_j_kg_k < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "caloric exchange requires a positive bath and finite nonnegative exposure",
            ));
        }
        let caloric = self.caloric_law();
        let initial_temperature = caloric.temperature_k(initial_specific_energy_j_kg)?;
        let (exchange, temperatures) = match self.original {
            OriginalConduction::Isobaric(law) => solve(
                &law.caloric().source_line().knots,
                1.0, initial_temperature, reservoir_temperature_k, exposure_j_kg_k,
            )?,
            OriginalConduction::ExpandingIsobaric(law) => solve(
                &law.caloric().caloric().source_line().knots,
                1.0, initial_temperature, reservoir_temperature_k, exposure_j_kg_k,
            )?,
            OriginalConduction::Isochoric(law) => solve(
                &law.caloric().source_line().knots,
                1.0, initial_temperature, reservoir_temperature_k, exposure_j_kg_k,
            )?,
            OriginalConduction::FixedVolume(law) => solve(
                &law.path().source_lines()[0].knots,
                law.reference_density_kg_m3(),
                initial_temperature, reservoir_temperature_k, exposure_j_kg_k,
            )?,
            OriginalConduction::HeatStorage(law) => solve(
                law.caloric().storage().capacity_knots(),
                law.caloric().reference_mass_kg(),
                initial_temperature, reservoir_temperature_k, exposure_j_kg_k,
            )?,
        };
        let heat = exchange.energy_change;
        let bounds = self.energy_bounds_j_kg();
        let final_energy = if heat == 0.0 {
            initial_specific_energy_j_kg
        } else if heat < 0.0 && exchange.temperature == temperatures[0] {
            bounds[0]
        } else if heat > 0.0 && exchange.temperature == temperatures[1] {
            bounds[1]
        } else {
            initial_specific_energy_j_kg + heat
        };
        if !final_energy.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        let temperature_k = caloric.temperature_k(final_energy)?;
        let represented = final_energy - initial_specific_energy_j_kg;
        let residual = heat - represented;
        if !represented.is_finite() || !residual.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        Ok(CaloricHeatExchange {
            storage_basis: caloric.storage_basis(),
            initial_specific_energy_j_kg,
            final_specific_energy_j_kg: final_energy,
            temperature_k,
            mean_temperature_k: exchange.mean_temperature,
            specific_heat_j_kg: heat,
            represented_energy_change_j_kg: represented,
            heat_residual_j_kg: residual,
            inlet_energy_gain: exchange.inlet_energy_gain,
            inlet_energy_loss: exchange.inlet_energy_loss,
            reservoir_energy_derivative_j_kg_k: exchange.reservoir_energy_derivative,
            mean_reservoir_derivative: exchange.mean_reservoir_derivative,
            remaining_temperature_difference_k: exchange.remaining_temperature_difference,
            mean_log_exposure_derivative_k: exchange.mean_log_exposure_derivative,
        })
    }
}

// Capacity is J/(kg K), J/(m3 K), or J/K. Solve in that original source basis,
// multiplying exposure by rho_ref or m_ref, then return specific heat/tangents.
// Temperature means and dimensionless gains are invariant under this conversion.
fn solve<K: CapacityKnot>(
    knots: &[K],
    energy_scale: f64,
    initial: f64,
    reservoir: f64,
    specific_exposure: f64,
) -> Result<(CapacityExchange, [f64; 2]), DepositionError> {
    let exposure = specific_exposure * energy_scale;
    if !exposure.is_finite() || (specific_exposure != 0.0 && exposure == 0.0) {
        return Err(DepositionError::Arithmetic);
    }
    let mut exchange = CapacityCurve::new(knots)
        .linear_heat_exchange_preserving_heat(initial, reservoir, exposure)
        .map_err(|_| DepositionError::MaterialDomain)?;
    let native_heat = exchange.energy_change;
    let native_derivative = exchange.reservoir_energy_derivative;
    exchange.energy_change /= energy_scale;
    exchange.reservoir_energy_derivative /= energy_scale;
    if !exchange.energy_change.is_finite()
        || !exchange.reservoir_energy_derivative.is_finite()
        || (native_heat != 0.0 && exchange.energy_change == 0.0)
        || (native_derivative != 0.0 && exchange.reservoir_energy_derivative == 0.0)
    {
        return Err(DepositionError::Arithmetic);
    }
    Ok((exchange, [knots[0].temperature(), knots[knots.len() - 1].temperature()]))
}
