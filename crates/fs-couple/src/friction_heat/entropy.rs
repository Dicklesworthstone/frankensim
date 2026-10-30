//! Physical friction entropy at the solved interface, with finite body storage.
//! The interface has no heat capacity. Its energy residual is therefore kept
//! separately from both irreversible production and the thermal time defect.

use super::{finite, invalid};
use fs_thermal::enthalpy::entropy::EntropyAudit;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

/// Entropy account of one coupled friction/thermal candidate, in J/K.
///
/// Caloric plus mechanical entropy storage equals reservoir input + network/interface conduction + friction
/// generation + other source-at-cell contribution + temporal defect + interface
/// residual/T, to floating-point summation. Acceptance still requires the
/// original cell, thermal, interface and total energy gates. Endpoint quadrature
/// is not an exact time integral, and other source watts need their own entropy
/// semantics before this can support a closed-system second-law claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionHeatEntropy {
    /// Exact original caloric integrals, including latent storage once.
    pub thermal_storage_change_j_k: f64,
    /// Mechanical entropy from an explicitly admitted temperature-dependent
    /// potential. Zero for isothermal potentials; this is storage, not production.
    pub mechanical_storage_change_j_k: f64,
    /// Signed reversible heat into the caloric owner divided by its endpoint
    /// temperature. This internal transfer is not friction generation or an
    /// unclassified external source; it is already inside caloric storage.
    pub reversible_caloric_exchange_j_k: f64,
    /// Actual ambient/support reservoir heat divided by its own temperature.
    pub reservoir_entropy_input_j_k: f64,
    /// Original network links and ambient/support contacts. The direct friction
    /// links are suppressed, so their heat cannot be counted here a second time.
    pub network_conduction_generation_j_k: f64,
    /// Sum Q_i*(1/T_i - 1/T_interface), using signed heat into each body.
    pub interface_conduction_generation_j_k: f64,
    /// Sum of each physical friction loss divided by its own interface Kelvin.
    /// Stored mechanical energy and numerical damping are not thermalized here.
    pub friction_generation_j_k: f64,
    /// OTHER authored/callback source energy divided by receiving cell Kelvin.
    /// Known internal friction heat has been removed from this classification.
    pub unclassified_source_at_cell_j_k: f64,
    /// Whether any non-friction net cell source was nonzero.
    pub has_unclassified_sources: bool,
    /// Remaining signed temporal/algebraic/roundoff defect, never physical heat.
    pub temporal_defect_j_k: f64,
    /// Sum of (Sum Q_i - D)/T_interface over the distinct physical interfaces;
    /// independently visible until energy acceptance.
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
    reservoir_correction: Sum,
    conduction_correction: Sum,
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
            if watts != 0.0 {
                has_other_source = true;
                other.add(finite(dt * watts)? / state.temperature_k(id)?)?;
            }
        }
        Ok(Self {
            conduction: Sum::default(),
            other_source_j_k: other.total()?,
            has_other_source,
            reservoir_correction: Sum::default(),
            conduction_correction: Sum::default(),
        })
    }
    // A moving footprint cancels the covered part of an existing Robin flux
    // through an equal-and-opposite source. That cancellation belongs to the
    // ORIGINAL reservoir, not to the unclassified source or friction accounts.
    pub(super) fn correct_reservoir_transfer(
        &mut self,
        correction_heat_j: f64,
        body_temperature_k: f64,
        reservoir_temperature_k: f64,
    ) -> Result<(), DepositionError> {
        if !body_temperature_k.is_finite() || body_temperature_k <= 0.0
            || !reservoir_temperature_k.is_finite() || reservoir_temperature_k <= 0.0
        {
            return Err(invalid("reservoir entropy needs positive finite temperatures"));
        }
        self.reservoir_correction.add(correction_heat_j / reservoir_temperature_k)?;
        self.conduction_correction.add((correction_heat_j / body_temperature_k)
            * ((reservoir_temperature_k - body_temperature_k) / reservoir_temperature_k))
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
    pub(super) fn with_mechanical_storage(mut self, entropy_change_j_k: f64,
        reversible_heat_j: f64, temperature_k: f64) -> Result<Self, DepositionError> {
        self.mechanical_storage_change_j_k = finite(entropy_change_j_k)?;
        self.reversible_caloric_exchange_j_k = finite(reversible_heat_j / temperature_k)?;
        // The original caloric defect contains the reversible transfer. Add
        // the separately owned mechanical entropy to close total storage;
        // neither term may be labeled irreversible friction production.
        self.temporal_defect_j_k = finite(self.temporal_defect_j_k + entropy_change_j_k)?;
        Ok(self)
    }

    pub(super) fn from_trial(
        thermal: &EntropyAudit,
        interface: InterfaceEntropy,
        dissipated_j: f64,
        interface_temperature_k: f64,
        interface_residual_j: f64,
    ) -> Result<Self, DepositionError> {
        Self::from_interfaces(
            thermal,
            interface,
            [(dissipated_j, interface_temperature_k, interface_residual_j)],
            || Ok(()),
        )
    }

    pub(super) fn from_interfaces(
        thermal: &EntropyAudit,
        interface: InterfaceEntropy,
        interfaces: impl IntoIterator<Item = (f64, f64, f64)>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        let mut friction = Sum::default();
        let mut residual = Sum::default();
        for (index, (dissipated_j, temperature_k, residual_j)) in
            interfaces.into_iter().enumerate()
        {
            if index % 256 == 0 {
                checkpoint()?;
            }
            friction.add(finite(dissipated_j / temperature_k)?)?;
            residual.add(finite(residual_j / temperature_k)?)?;
        }
        let friction = friction.total()?;
        let conduction = interface.conduction.total()?;
        let interface_residual = residual.total()?;
        let mut reservoir = interface.reservoir_correction;
        reservoir.add(thermal.reservoir_entropy_input_j_k)?;
        let reservoir = reservoir.total()?;
        let mut network = interface.conduction_correction;
        network.add(thermal.conduction_generation_j_k)?;
        let network = network.total()?;
        let mut defect = Sum::default();
        for value in [
            thermal.thermal_storage_change_j_k,
            -reservoir,
            -network,
            -conduction,
            -friction,
            -interface.other_source_j_k,
            -interface_residual,
        ] {
            defect.add(value)?;
        }
        Ok(Self {
            thermal_storage_change_j_k: thermal.thermal_storage_change_j_k,
            mechanical_storage_change_j_k: 0.0,
            reversible_caloric_exchange_j_k: 0.0,
            reservoir_entropy_input_j_k: reservoir,
            network_conduction_generation_j_k: network,
            interface_conduction_generation_j_k: conduction,
            friction_generation_j_k: friction,
            unclassified_source_at_cell_j_k: interface.other_source_j_k,
            has_unclassified_sources: interface.has_other_source,
            temporal_defect_j_k: defect.total()?,
            interface_residual_entropy_j_k: interface_residual,
        })
    }
}
