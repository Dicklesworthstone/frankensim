//! Physical friction entropy at the solved interface, with finite body storage.
//! The interface has no heat capacity. Its energy residual is therefore kept
//! separately from both irreversible production and the thermal time defect.

use super::{finite, invalid};
use fs_thermal::enthalpy::entropy::EntropyAudit;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

/// Entropy account of one coupled friction/thermal candidate, in J/K.
///
/// Storage equals reservoir input + network/interface conduction + friction
/// generation + other source-at-cell contribution + temporal defect + interface
/// residual/T, to floating-point summation. Acceptance still requires the
/// original cell, thermal, interface and total energy gates. Endpoint quadrature
/// is not an exact time integral, and other source watts need their own entropy
/// semantics before this can support a closed-system second-law claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionHeatEntropy {
    /// Exact original caloric integrals, including latent storage once.
    pub thermal_storage_change_j_k: f64,
    /// Actual ambient/support reservoir heat divided by its own temperature.
    pub reservoir_entropy_input_j_k: f64,
    /// Original network links and ambient/support contacts. The direct friction
    /// link is suppressed, so its heat cannot be counted here a second time.
    pub network_conduction_generation_j_k: f64,
    /// Sum Q_i*(1/T_i - 1/T_interface), using signed heat into each body.
    pub interface_conduction_generation_j_k: f64,
    /// Exported physical friction loss divided by the solved interface Kelvin.
    /// Stored mechanical energy and numerical damping are not thermalized here.
    pub friction_generation_j_k: f64,
    /// OTHER authored/callback source energy divided by receiving cell Kelvin.
    /// Known internal friction heat has been removed from this classification.
    pub unclassified_source_at_cell_j_k: f64,
    /// Whether any non-friction net cell source was nonzero.
    pub has_unclassified_sources: bool,
    /// Remaining signed temporal/algebraic/roundoff defect, never physical heat.
    pub temporal_defect_j_k: f64,
    /// (Sum Q_i - D)/T_interface; independently visible until energy acceptance.
    pub interface_residual_entropy_j_k: f64,
}
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let value = finite(value)?;
        let next = finite(self.value + value)?;
        let lost = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, DepositionError> {
        finite(self.value + self.correction)
    }
}
#[derive(Clone, Copy)]
pub(super) struct InterfaceEntropy {
    conduction: Sum,
    other_source_j_k: f64,
    has_other_source: bool,
}
impl InterfaceEntropy {
    // Read other sources BEFORE adding friction. Subtracting a large friction
    // contribution afterwards could erase a small, real external source.
    pub(super) fn new(
        state: &DepositionSimulation<'_>,
        sources: &BTreeMap<CellId, f64>,
        dt: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        let mut other = Sum::default();
        let mut has_other_source = false;
        for (index, (&id, cell)) in state.cells().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let watts = finite(cell.boundary.source_w + sources.get(&id).copied().unwrap_or(0.0))?;
            has_other_source |= watts != 0.0;
            other.add(finite(dt * watts)? / state.temperature_k(id)?)?;
        }
        Ok(Self {
            conduction: Sum::default(),
            other_source_j_k: other.total()?,
            has_other_source,
        })
    }
    pub(super) fn add_transfer(
        &mut self,
        heat_j: f64,
        body_temperature_k: f64,
        interface_temperature_k: f64,
    ) -> Result<(), DepositionError> {
        let production = finite(
            (heat_j / body_temperature_k)
                * ((interface_temperature_k - body_temperature_k) / interface_temperature_k),
        )?;
        if production < 0.0 {
            return Err(invalid(
                "friction interface entropy transfer is not passive",
            ));
        }
        self.conduction.add(production)
    }
}
impl FrictionHeatEntropy {
    pub(super) fn from_trial(
        thermal: &EntropyAudit,
        interface: InterfaceEntropy,
        dissipated_j: f64,
        interface_temperature_k: f64,
        interface_residual_j: f64,
    ) -> Result<Self, DepositionError> {
        let friction = finite(dissipated_j / interface_temperature_k)?;
        let conduction = interface.conduction.total()?;
        let interface_residual = finite(interface_residual_j / interface_temperature_k)?;
        let mut defect = Sum::default();
        for value in [
            thermal.thermal_storage_change_j_k,
            -thermal.reservoir_entropy_input_j_k,
            -thermal.conduction_generation_j_k,
            -conduction,
            -friction,
            -interface.other_source_j_k,
            -interface_residual,
        ] {
            defect.add(value)?;
        }
        Ok(Self {
            thermal_storage_change_j_k: thermal.thermal_storage_change_j_k,
            reservoir_entropy_input_j_k: thermal.reservoir_entropy_input_j_k,
            network_conduction_generation_j_k: thermal.conduction_generation_j_k,
            interface_conduction_generation_j_k: conduction,
            friction_generation_j_k: friction,
            unclassified_source_at_cell_j_k: interface.other_source_j_k,
            has_unclassified_sources: interface.has_other_source,
            temporal_defect_j_k: defect.total()?,
            interface_residual_entropy_j_k: interface_residual,
        })
    }
}
