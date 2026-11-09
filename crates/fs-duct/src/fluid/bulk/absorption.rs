//! Physical bulk-wave absorption of the original loaded traction cylinder.

use super::{BulkTractionCylinder, BulkTractionResponse, ComposeOptions, FluidDuctError, bad};
use fs_material::fluid_acoustics::thermoviscous::absorption::{
    BulkCoherentCellAbsorption, BulkCoherentWaveSegment, BulkThermoviscousEquilibrium,
    BulkWaveField,
};
use fs_math::c64::C64;

/// Loaded original cylinder, now carrying physical peak amplitudes and its
/// coherent bulk absorption. The same acoustic eigenmode supplies temperature
/// and heat-flux phasors; no independent thermal diffusion mode or thermal end
/// condition is introduced. Wall-layer and molecular losses remain absent.
///
/// Full acoustic flux includes stress work and thermal availability transport.
/// It need not equal pressure work or mechanical traction work alone. Direct
/// viscous heat and conductive entropy production remain separately reported.
#[derive(Clone, Copy, Debug)]
pub struct BulkTractionAbsorptionResponse {
    response: BulkTractionResponse,
    segment: BulkCoherentWaveSegment,
    inlet_flow_m3_s: C64,
    area_m2: f64,
}

impl BulkTractionCylinder {
    /// Solve the SAME traction load and bind the physical inlet volume flow to
    /// exact finite-distance absorption. Flow is a peak exp(-i omega t) phasor
    /// [m³/s], not RMS or intensity. `None` retains the original rigid outlet.
    ///
    /// The supplied thermal equilibrium must come from this cylinder's original
    /// material. Its entire checked acoustic/caloric/transport/bulk-viscosity
    /// tuple must match exactly; equal tuples do not invent source provenance.
    /// Signed expansion and absolute temperature are taken from that equilibrium.
    /// The original load solver and condition refusals remain authoritative.
    /// Endpoint-anchored modes never require a growing exponential. Once bound,
    /// field and cell queries cost O(1), allocate nothing and perform no solves.
    pub fn absorption_response(
        &self,
        equilibrium: &BulkThermoviscousEquilibrium,
        traction_load_pa_s_m3: Option<C64>,
        inlet_flow_m3_s: C64,
        options: ComposeOptions,
    ) -> Result<BulkTractionAbsorptionResponse, FluidDuctError> {
        if equilibrium.absorption() != &self.material {
            return Err(
                bad("bulk absorption requires the cylinder's original checked material").into(),
            );
        }
        if !inlet_flow_m3_s.re.is_finite() || !inlet_flow_m3_s.im.is_finite() {
            return Err(bad("bulk absorption requires finite peak inlet volume flow").into());
        }
        let absorption = equilibrium.prepare_wave_absorption(self.omega_rad_s)?;
        let response = self.response(traction_load_pa_s_m3, options)?;
        let traction = self.wave.traction_impedance_pa_s_m();
        let velocity = inlet_flow_m3_s / C64::from_re(self.area_m2);
        let inlet_traction = response.inlet_traction_impedance_pa_s_m3 * inlet_flow_m3_s;
        let forward = velocity.scale(0.5) + (inlet_traction / traction).scale(0.5);
        let backward = self.backward_at_outlet(&response, traction_load_pa_s_m3, inlet_flow_m3_s);
        let segment = absorption.coherent_segment(self.length_m, forward, backward)?;
        Ok(BulkTractionAbsorptionResponse {
            response,
            segment,
            inlet_flow_m3_s,
            area_m2: self.area_m2,
        })
    }

    fn backward_at_outlet(
        &self,
        response: &BulkTractionResponse,
        load: Option<C64>,
        inlet_flow: C64,
    ) -> C64 {
        let traction = self.wave.traction_impedance_pa_s_m();
        let velocity = (response.outlet_over_inlet_flow * inlet_flow) / C64::from_re(self.area_m2);
        if let Some(load) = load {
            // Match the original constructor's supplied characteristic bits.
            let characteristic = traction.scale(1.0 / self.area_m2);
            let difference = characteristic - load;
            if difference.abs() <= 0.5 * characteristic.abs() {
                // B=(1-Z_L/Z_c) U_out/(2 A). Near matching, subtract the
                // ORIGINAL impedances before division, not rounded endpoint
                // pressure and velocity. Exact matching then gives exact B=0;
                // a tiny reflected wave keeps the actual load difference.
                return (difference / characteristic) * velocity.scale(0.5);
            }
        }
        // General loads, including rigid flow=0 and active reference poles,
        // use the already resolved native endpoint traction. No load ratio
        // multiplies a tiny or rounded-to-zero outlet flow here.
        let outlet_traction = response.outlet_traction_over_inlet_flow_pa_s_m3 * inlet_flow;
        velocity.scale(0.5) - (outlet_traction / traction).scale(0.5)
    }
}

impl BulkTractionAbsorptionResponse {
    /// Original unit-inlet-flow mechanical load response and condition estimate.
    pub const fn response(&self) -> BulkTractionResponse {
        self.response
    }
    pub const fn inlet_volume_flow_m3_s(&self) -> C64 {
        self.inlet_flow_m3_s
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn length_m(&self) -> f64 {
        self.segment.length_m()
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.segment.absorption().angular_frequency_rad_s()
    }
    pub const fn equilibrium(&self) -> &BulkThermoviscousEquilibrium {
        self.segment.absorption().equilibrium()
    }

    /// Physical pressure, traction, velocity, temperature and conductive flux,
    /// plus signed +x acoustic energy flux and distinct positive local losses.
    pub fn at_distance(&self, distance_m: f64) -> Result<BulkWaveField, FluidDuctError> {
        Ok(self.segment.at_distance(distance_m)?)
    }

    /// Exact coherent interval absorption with the ORIGINAL physical area
    /// applied once. Endpoints lie within this cylinder. The independently
    /// retained balance residual includes thermal availability flux at both
    /// faces; only viscous_heat_w is a direct caloric heat source.
    pub fn absorption_between(
        &self,
        start_m: f64,
        end_m: f64,
    ) -> Result<BulkCoherentCellAbsorption, FluidDuctError> {
        Ok(self
            .segment
            .absorption_between(start_m, end_m, self.area_m2)?)
    }
}

#[cfg(test)]
mod tests;
