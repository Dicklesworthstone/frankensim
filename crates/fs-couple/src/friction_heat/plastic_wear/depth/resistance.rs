//! Original thermal resistance at the fully advanced depth-contact endpoint.
//!
//! The native equal-area pressure sampler and source binding are shared with
//! ordinary plastic sliding. Every caloric trial samples after temperature,
//! recession and exposed-hardness continuations. No beginning pressure, nominal
//! whole-face area, hardness-derived real area or extra heat partition is used.

use fs_exec::Cx;
use fs_material::state_point::ResolvedInterfaceStatePoint;

use super::super::resistance::{PlasticThermalResistance, PlasticThermalResistanceSample};
use super::super::{PlasticSlidingHeatError, checkpoint};
use super::{DepthSlidingHeatContact, SourceDepthSlidingWear};

impl PlasticThermalResistance {
    pub(super) fn sample_depth<E>(
        &self,
        source: &SourceDepthSlidingWear,
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<PlasticThermalResistanceSample, PlasticSlidingHeatError<E>> {
        let wear = source.wear_state();
        self.sample_native_pressure(
            wear.normal_source().elastic().interface().card(),
            wear.contact_state().pressure_pa(),
            interface_temperature_k,
            cx,
        )
    }
}

impl DepthSlidingHeatContact<'_> {
    /// Bind R(T_interface, mean_loaded_pressure, original fixed conditions).
    ///
    /// The canonical resistance claim must belong to the same ordered interface,
    /// with a distinct typed static-pressure axis and the original temperature
    /// role. The existing heat route must already carry its resolved split.
    /// `fraction_to_first` locates the massless generation plane in total contact
    /// resistance; the native finite-cell heat solve determines heat partition.
    ///
    /// Each trial uses final pressure after wear and exposed-depth hardening or
    /// softening. The source and split survive accepted steps, material export
    /// and adaptive continuation. A support refusal returns no partial state.
    pub fn with_sourced_contact_resistance(
        mut self,
        source: ResolvedInterfaceStatePoint,
        pressure_axis: String,
        fraction_to_first: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        checkpoint(cx)?;
        let resistance = PlasticThermalResistance::bind(
            source,
            pressure_axis,
            fraction_to_first,
            self.source.wear_state().normal_source().interface_state(),
            self.source.friction_source().axes().temperature.clone(),
            self.source.interface_temperature_k(),
        )?;
        resistance
            .sample_depth(&self.source, self.source.interface_temperature_k(), cx)?
            .admit_route(&self.route)?;
        self.resistance = Some(resistance);
        Ok(self)
    }

    pub const fn contact_resistance(&self) -> Option<&PlasticThermalResistance> {
        self.resistance.as_ref()
    }
}
