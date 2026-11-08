//! Exact fixed-temperature-bath exchange in the original energy coordinate.

use super::{CaloricConductionLayer, OriginalConduction};
use super::super::capacity::{CapacityCurve, CapacityDrivenExchange, CapacityExchange, CapacityKnot};
use fs_matdb::MaterialAnswer;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};
use crate::thermal::conduction::ResolvedConductionLayer;

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

/// Exact constant applied-source plus linear-bath response on native storage.
///
/// Inputs obey C_specific(T) dT/ds = Q0 - X*(T-T0), 0 <= s <= 1.
/// Q0 is the signed initial driving power integrated over the whole window;
/// X is its nonnegative conductance exposure. A virtual equilibrium may be
/// negative or unrepresentable: only the actual material trajectory is a
/// physical temperature and must remain inside the original source support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricDrivenHeatExchange {
    pub storage_basis: ThermalStorageBasis,
    pub initial_specific_energy_j_kg: f64,
    pub final_specific_energy_j_kg: f64,
    /// Original inverse of the returned final native energy [K].
    pub temperature_k: f64,
    /// Endpoint minus initial temperature, retained before absolute rounding.
    pub temperature_change_k: f64,
    /// Exact time mean, rounded to absolute kelvin.
    pub mean_temperature_k: f64,
    /// Time mean minus initial temperature, retained before absolute rounding.
    pub mean_temperature_change_k: f64,
    /// Signed integrated heat into the material [J/kg], before storage rounding.
    pub specific_heat_j_kg: f64,
    pub represented_energy_change_j_kg: f64,
    pub heat_residual_j_kg: f64,
    /// d(e_out)/d(e_in) holding the physical applied power, baths and X fixed.
    /// Q0 therefore changes by -X*dT0 under this perturbation.
    pub inlet_energy_gain: f64,
    /// One minus inlet_energy_gain, retained directly for weak feedback.
    pub inlet_energy_loss: f64,
    /// d(e_out)/d(Q0), holding inlet energy and X fixed; dimensionless.
    pub source_heat_gain: f64,
    /// d(T_mean)/d(Q0) at fixed inlet energy and X [kg K/J].
    pub mean_source_derivative_kg_k_j: f64,
    /// Q0 - X*(T_out-T0), retained before cancellation or absolute-T rounding.
    pub remaining_drive_j_kg: f64,
    /// d(T_mean)/d(log duration) when Q0 and X scale together [K].
    pub mean_log_duration_derivative_k: f64,
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
    ///
    /// Equilibrium-phase storage delegates to its original exact phase flow
    /// and response, retaining its response-domain refusals. The complete path
    /// must stay in this layer's prepared pure-phase/source span; this method
    /// does not infer a homogeneous mixed-phase conductivity profile.
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
            OriginalConduction::EquilibriumPhase(law) => {
                return phase_exchange(
                    law, initial_specific_energy_j_kg,
                    reservoir_temperature_k, exposure_j_kg_k,
                );
            }
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
        let (final_energy, temperature_k, represented, residual) = self.finish_energy(
            initial_specific_energy_j_kg, heat, exchange.temperature, temperatures,
        )?;
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

    /// Exact constant source plus held linear baths, expressed without an
    /// absolute weighted bath or virtual-equilibrium temperature.
    ///
    /// Let Q0 = dt*(P_source + sum G_i*(T_i-T0))/m [J/kg] and
    /// X = dt*sum G_i/m [J/(kg K)], where T0 is the original inverse of the
    /// input energy. This method integrates C(T)dT/ds=Q0-X*(T-T0) on 0..1.
    /// The caller retains the individual bath temperatures/conductances and
    /// applied-source account; their model validity and cell inventory remain
    /// with that owner. Positive and negative applied power are supported.
    ///
    /// Zero X is exact calorimetry with the capacity-weighted cubic time mean,
    /// not an arithmetic midpoint. Zero Q0 still retains the correct source
    /// and bath sensitivities. Cp, Cv, reference-volume and amount-bound laws
    /// retain their original storage basis and source support. All crossed
    /// capacities and analytic moments use the SAME shared exchange kernel.
    ///
    /// Equilibrium-phase layers retain their original phase exchange API and
    /// are not admitted by this sensible-storage response.
    ///
    /// Successful queries allocate nothing, query no material card and mutate
    /// no thermal state. Relative heat and temperature moments survive absolute
    /// rounding; the returned energy is admitted by the original native inverse.
    pub fn driven_heat_exchange(
        &self,
        initial_specific_energy_j_kg: f64,
        initial_drive_j_kg: f64,
        exposure_j_kg_k: f64,
    ) -> Result<CaloricDrivenHeatExchange, DepositionError> {
        if !initial_drive_j_kg.is_finite()
            || !exposure_j_kg_k.is_finite() || exposure_j_kg_k < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "driven caloric exchange requires finite initial drive and nonnegative exposure",
            ));
        }
        let caloric = self.caloric_law();
        let initial_temperature = caloric.temperature_k(initial_specific_energy_j_kg)?;
        let (exchange, temperatures) = match self.original {
            OriginalConduction::EquilibriumPhase(_) => {
                return Err(DepositionError::InvalidInput(
                    "driven caloric response is unsupported for equilibrium-phase storage",
                ));
            }
            OriginalConduction::Isobaric(law) => solve_driven(
                &law.caloric().source_line().knots, 1.0,
                initial_temperature, initial_drive_j_kg, exposure_j_kg_k,
            )?,
            OriginalConduction::ExpandingIsobaric(law) => solve_driven(
                &law.caloric().caloric().source_line().knots, 1.0,
                initial_temperature, initial_drive_j_kg, exposure_j_kg_k,
            )?,
            OriginalConduction::Isochoric(law) => solve_driven(
                &law.caloric().source_line().knots, 1.0,
                initial_temperature, initial_drive_j_kg, exposure_j_kg_k,
            )?,
            OriginalConduction::FixedVolume(law) => solve_driven(
                &law.path().source_lines()[0].knots, law.reference_density_kg_m3(),
                initial_temperature, initial_drive_j_kg, exposure_j_kg_k,
            )?,
            OriginalConduction::HeatStorage(law) => solve_driven(
                law.caloric().storage().capacity_knots(), law.caloric().reference_mass_kg(),
                initial_temperature, initial_drive_j_kg, exposure_j_kg_k,
            )?,
        };
        let (final_energy, temperature_k, represented, residual) = self.finish_energy(
            initial_specific_energy_j_kg, exchange.energy_change, exchange.temperature, temperatures,
        )?;
        Ok(CaloricDrivenHeatExchange {
            storage_basis: caloric.storage_basis(),
            initial_specific_energy_j_kg,
            final_specific_energy_j_kg: final_energy,
            temperature_k,
            temperature_change_k: exchange.temperature_change,
            mean_temperature_k: exchange.mean_temperature,
            mean_temperature_change_k: exchange.mean_temperature_change,
            specific_heat_j_kg: exchange.energy_change,
            represented_energy_change_j_kg: represented,
            heat_residual_j_kg: residual,
            inlet_energy_gain: exchange.inlet_energy_gain,
            inlet_energy_loss: exchange.inlet_energy_loss,
            source_heat_gain: exchange.source_gain,
            mean_source_derivative_kg_k_j: exchange.mean_source_derivative,
            remaining_drive_j_kg: exchange.remaining_drive,
            mean_log_duration_derivative_k: exchange.mean_log_duration_derivative,
        })
    }

    // Both exchange entry points retain the same native energy admission and
    // canonical support endpoints; neither reconstructs heat from rounded T.
    fn finish_energy(
        &self,
        initial: f64,
        heat: f64,
        solved_temperature: f64,
        temperatures: [f64; 2],
    ) -> Result<(f64, f64, f64, f64), DepositionError> {
        let bounds = self.energy_bounds_j_kg();
        let final_energy = if heat == 0.0 {
            initial
        } else if heat < 0.0 && solved_temperature == temperatures[0] {
            bounds[0]
        } else if heat > 0.0 && solved_temperature == temperatures[1] {
            bounds[1]
        } else {
            initial + heat
        };
        if !final_energy.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        let temperature = self.caloric_law().temperature_k(final_energy)?;
        let represented = final_energy - initial;
        let residual = heat - represented;
        if !represented.is_finite() || !residual.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        Ok((final_energy, temperature, represented, residual))
    }

}

// Reuse the original phase traversal and every stable response moment.
// The fixed-bath trajectory is monotone in h, so admitting its two endpoints
// admits its complete path on the prepared pure-phase/source interval.
fn phase_exchange(
    layer: &ResolvedConductionLayer<'_>,
    initial_h: f64,
    reservoir: f64,
    exposure: f64,
) -> Result<CaloricHeatExchange, DepositionError> {
    layer.temperature_at_enthalpy(initial_h).map_err(|_| DepositionError::MaterialDomain)?;
    let curve = layer.caloric_curve();
    let initial = curve.state_at_specific_enthalpy(initial_h)
        .map_err(|_| DepositionError::MaterialDomain)?;
    let response = curve.linear_heat_exchange_response(initial, reservoir, exposure)
        .map_err(|_| DepositionError::MaterialDomain)?;
    let (exchange, jacobian) = (response.exchange, response.jacobian);
    let final_h = exchange.state.specific_enthalpy_j_kg();
    layer.temperature_at_enthalpy(final_h).map_err(|_| DepositionError::MaterialDomain)?;
    Ok(CaloricHeatExchange {
        storage_basis: curve.storage_basis(),
        initial_specific_energy_j_kg: initial_h,
        final_specific_energy_j_kg: final_h,
        temperature_k: exchange.state.temperature_k(),
        mean_temperature_k: exchange.mean_temperature_k,
        specific_heat_j_kg: exchange.specific_heat_j_kg,
        represented_energy_change_j_kg: exchange.represented_enthalpy_change_j_kg,
        heat_residual_j_kg: exchange.heat_residual_j_kg,
        inlet_energy_gain: jacobian.inlet_enthalpy_gain,
        inlet_energy_loss: jacobian.inlet_enthalpy_loss,
        reservoir_energy_derivative_j_kg_k: jacobian.reservoir_enthalpy_derivative_j_kg_k,
        mean_reservoir_derivative: jacobian.mean_reservoir_derivative,
        remaining_temperature_difference_k: jacobian.exposure_enthalpy_derivative_k,
        mean_log_exposure_derivative_k: jacobian.mean_log_exposure_derivative_k,
    })
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


fn solve_driven<K: CapacityKnot>(
    knots: &[K],
    energy_scale: f64,
    initial: f64,
    specific_drive: f64,
    specific_exposure: f64,
) -> Result<(CapacityDrivenExchange, [f64; 2]), DepositionError> {
    let drive = specific_drive * energy_scale;
    let exposure = specific_exposure * energy_scale;
    if !drive.is_finite() || !exposure.is_finite()
        || (specific_drive != 0.0 && drive == 0.0)
        || (specific_exposure != 0.0 && exposure == 0.0)
    {
        return Err(DepositionError::Arithmetic);
    }
    let mut exchange = CapacityCurve::new(knots)
        .driven_heat_exchange(initial, drive, exposure)
        .map_err(|_| DepositionError::MaterialDomain)?;
    let native_heat = exchange.energy_change;
    let native_mean_response = exchange.mean_source_derivative;
    exchange.energy_change /= energy_scale;
    exchange.remaining_drive /= energy_scale;
    // Output heat and input drive scale together, so source_gain is invariant.
    // The derivative of a temperature with respect to specific input energy is
    // instead the native derivative multiplied by the original mass or density.
    exchange.mean_source_derivative *= energy_scale;
    if specific_exposure == 0.0 {
        // Preserve the exact specific source input through mass/volume scaling.
        exchange.energy_change = specific_drive;
        exchange.remaining_drive = specific_drive;
    }
    if !exchange.energy_change.is_finite()
        || !exchange.remaining_drive.is_finite()
        || !exchange.mean_source_derivative.is_finite()
        || (native_heat != 0.0 && exchange.energy_change == 0.0)
        || (native_mean_response != 0.0 && exchange.mean_source_derivative == 0.0)
    {
        return Err(DepositionError::Arithmetic);
    }
    Ok((exchange, [knots[0].temperature(), knots[knots.len() - 1].temperature()]))
}

