//! Prescribed-reservoir backward Euler on the original fixed-F inventory.

use super::{
    FiniteSolidCaloricState, FiniteThermoelasticHeatBalance, FiniteThermoelasticHeatError,
    FiniteThermoelasticHeatSolve, Result, finite, invalid, poll,
};
use fs_thermal::enthalpy::Reservoir;

/// Accepted finite-duration boundary exchange and separate entropy accounts.
/// The heat report retains actual Q=dt*G*(Tb-Tf), independently checked against
/// the original internal-energy increment. Backward Euler approximates the
/// finite-time evolution; its temporal entropy defect is not physical loss.
#[derive(Clone, Debug)]
pub struct FiniteThermoelasticReservoirBalance {
    heat: FiniteThermoelasticHeatBalance,
    reservoir: Reservoir,
    duration_s: f64,
    endpoint_power_into_solid_w: f64,
    endpoint_conduction_entropy_j_per_k: f64,
    temporal_entropy_defect_j_per_k: f64,
    combined_entropy_change_j_per_k: f64,
}

impl FiniteSolidCaloricState {
    /// Advance a restrained uniform inventory through one prescribed reservoir.
    ///
    /// The authored total conductance G [W/K] and actual reservoir temperature
    /// remain constant for this duration. Area is retained as metadata and is
    /// not multiplied into G again. The caller owns the uniform-temperature
    /// approximation and the physical conductance; no contact pressure, spatial
    /// gradients, bath heat capacity or other boundary channel is inferred.
    ///
    /// Solve ΔU(T)+dt*G*(T-Tb)=0 with the original c_F+dt*G tangent and the SAME
    /// bounded source sampling, inverse, heat/entropy closure and work budget as
    /// `add_heat`. Only material temperatures need source support; the bath may
    /// lie outside that interval. This is a backward-Euler step, not an exact
    /// time path, adaptive integrator, global stability certificate or CFL bound.
    ///
    /// Zero duration/conductance or equal temperatures preserve the exact state
    /// after descriptor/control admission. Unresolved nonzero heat, source exit,
    /// budget exhaustion or cancellation publishes no endpoint.
    pub fn exchange_with_reservoir(
        &self,
        reservoir: Reservoir,
        duration_s: f64,
        controls: FiniteThermoelasticHeatSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteThermoelasticReservoirBalance> {
        poll(&mut cancelled)?;
        if !reservoir.area_m2.is_finite()
            || reservoir.area_m2 < 0.0
            || !reservoir.conductance_w_k.is_finite()
            || reservoir.conductance_w_k < 0.0
            || !reservoir.temperature_k.is_finite()
            || reservoir.temperature_k <= 0.0
            || !duration_s.is_finite()
            || duration_s < 0.0
        {
            return Err(invalid(
                "reservoir needs nonnegative area, total conductance and duration, and positive absolute temperature",
            ));
        }
        let inactive = duration_s == 0.0
            || reservoir.conductance_w_k == 0.0
            || reservoir.temperature_k == self.temperature_k();
        let conductance_time = if inactive {
            0.0
        } else {
            finite(duration_s * reservoir.conductance_w_k)?
        };
        if !inactive && conductance_time == 0.0 {
            return Err(invalid("nonzero reservoir coupling is unrepresentable"));
        }
        let heat = self.solve_heat(
            controls,
            |temperature| {
                let difference = reservoir.temperature_k - temperature;
                let heat = finite(conductance_time * difference)?;
                if conductance_time > 0.0 && difference != 0.0 && heat == 0.0 {
                    return Err(invalid("nonzero reservoir heat is unrepresentable"));
                }
                Ok(heat)
            },
            -conductance_time,
            &mut cancelled,
        )?;
        let final_temperature = heat.final_state().temperature_k();
        if final_temperature < self.temperature_k().min(reservoir.temperature_k)
            || final_temperature > self.temperature_k().max(reservoir.temperature_k)
        {
            return Err(invalid(
                "reservoir exchange cannot cross its bath temperature",
            ));
        }
        let difference = reservoir.temperature_k - final_temperature;
        let endpoint_power_into_solid_w = finite(reservoir.conductance_w_k * difference)?;
        if reservoir.conductance_w_k > 0.0
            && difference != 0.0
            && endpoint_power_into_solid_w == 0.0
        {
            return Err(invalid("nonzero reservoir power is unrepresentable"));
        }
        let endpoint_heat_entropy = finite(heat.heat_j() / final_temperature)?;
        let bath_heat_entropy = finite(heat.heat_j() / reservoir.temperature_k)?;
        // Stable Q*(1/Tf-1/Tb), without subtracting nearly equal reciprocals.
        let endpoint_conduction_entropy_j_per_k =
            finite(endpoint_heat_entropy * (difference / reservoir.temperature_k))?;
        let temporal_entropy_defect_j_per_k =
            finite(heat.entropy_change_j_per_k() - endpoint_heat_entropy)?;
        let combined_entropy_change_j_per_k =
            finite(endpoint_conduction_entropy_j_per_k + temporal_entropy_defect_j_per_k)?;
        let allowance = self.entropy_allowance(
            heat.final_state(),
            controls,
            endpoint_heat_entropy.abs().max(bath_heat_entropy.abs()),
        ) + controls.thermal.energy_tolerance_j / reservoir.temperature_k;
        if !allowance.is_finite()
            || endpoint_conduction_entropy_j_per_k < -allowance
            || temporal_entropy_defect_j_per_k < -allowance
            || combined_entropy_change_j_per_k < -allowance
        {
            return Err(FiniteThermoelasticHeatError::ReservoirEntropy {
                endpoint_conduction_j_per_k: endpoint_conduction_entropy_j_per_k,
                temporal_defect_j_per_k: temporal_entropy_defect_j_per_k,
            });
        }
        poll(&mut cancelled)?;
        Ok(FiniteThermoelasticReservoirBalance {
            heat,
            reservoir,
            duration_s,
            endpoint_power_into_solid_w,
            endpoint_conduction_entropy_j_per_k,
            temporal_entropy_defect_j_per_k,
            combined_entropy_change_j_per_k,
        })
    }
}

impl FiniteThermoelasticReservoirBalance {
    pub const fn heat_balance(&self) -> &FiniteThermoelasticHeatBalance {
        &self.heat
    }
    pub const fn final_state(&self) -> &FiniteSolidCaloricState {
        self.heat.final_state()
    }
    pub const fn reservoir(&self) -> &Reservoir {
        &self.reservoir
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Actual final boundary power, positive into the solid [W].
    pub const fn endpoint_power_into_solid_w(&self) -> f64 {
        self.endpoint_power_into_solid_w
    }
    /// Endpoint conduction Q*(1/Tf-1/Tb), excluding temporal defect [J/K].
    pub const fn endpoint_conduction_entropy_j_per_k(&self) -> f64 {
        self.endpoint_conduction_entropy_j_per_k
    }
    /// Backward-Euler ΔS_body-Q/Tf; not a physical dissipation channel [J/K].
    pub const fn temporal_entropy_defect_j_per_k(&self) -> f64 {
        self.temporal_entropy_defect_j_per_k
    }
    /// ΔS_body-Q/Tb, including both endpoint conduction and temporal defect.
    pub const fn combined_entropy_change_j_per_k(&self) -> f64 {
        self.combined_entropy_change_j_per_k
    }
}
