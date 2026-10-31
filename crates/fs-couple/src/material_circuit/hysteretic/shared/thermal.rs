//! One retained shared-core history on the existing finite enthalpy timeline.
//!
//! Every thermal trial restarts the SAME accepted winding linkages and common
//! magnetic memory. Each original resistor's native Joule heat enters its own
//! material receiver; the common core's actual pinning heat enters one whole
//! fixed solid caloric cell, exactly once. Leakage storage is reversible and
//! adds no core heat. Native winding voltages, currents and total linkages stay
//! available on the accepted shared circuit and its original step report.
//!
//! The existing source-pinned MaterialNetwork, complete HeatAssembly, finite
//! coolant, thermal-contact and radiation paths all use this same owner.
//! Temperature changes source-resolved winding resistance and hence the common
//! magnetic trajectory. The source-backed core fit remains constant over its
//! explicitly admitted temperature interval; no unmeasured thermal magnetic law
//! or empirical core-loss source is introduced.

use super::super::{
    HystereticCircuitError, HystereticCircuitSolve,
    thermal::{
        HystereticEnthalpyCandidate, HystereticThermalHistory, HystereticThermalStep,
        ThermalHystereticCircuit, sealed,
    },
};
use super::{
    SharedHystereticCircuit, SharedHystereticCircuitCandidate, SharedHystereticCircuitStep,
};
use crate::{
    circuit_enthalpy::material::MaterialNetworkSample,
    material_circuit::{MaterialCircuit, MaterialCircuitStep, MaterialResistor},
};
use fs_matdb::QueryPoint;
use fs_material::field_storage::hysteresis::UniformHystereticCore;

/// Shared-core electrical history, ONE core heat receiver and the original
/// source-pinned finite thermal inventory. Its bind/step/HeatAssembly methods
/// are the same native implementations used by the single-core owner.
pub type ThermalSharedHystereticCircuit<'n, 'a> =
    ThermalHystereticCircuit<'n, 'a, SharedHystereticCircuit>;

/// Publishable shared electrical/magnetic and finite thermal histories.
/// `report().electrical.magnetic.core.core_heat_j` is the single physical core
/// receipt. The report's legacy `deposited_joule_j` includes that core heat
/// together with the original resistor heat, never a numerical defect.
pub type SharedHystereticEnthalpyCandidate<'t, 'n, 'a, S = MaterialNetworkSample> =
    HystereticEnthalpyCandidate<'t, 'n, 'a, S, SharedHystereticCircuit>;

impl sealed::History for SharedHystereticCircuit {}
impl HystereticThermalHistory for SharedHystereticCircuit {
    type Step = SharedHystereticCircuitStep;

    fn uniform_core(&self) -> &UniformHystereticCore {
        self.core().core()
    }
    fn native_circuit(&self) -> &MaterialCircuit {
        SharedHystereticCircuit::native_circuit(self)
    }
    fn time_s(&self) -> f64 {
        SharedHystereticCircuit::time_s(self)
    }
    fn with_core_source_query(&self, query: &QueryPoint) -> Result<Self, HystereticCircuitError> {
        SharedHystereticCircuit::with_core_source_query(self, query)
    }
    fn step_for_heat(
        &self,
        target_s: f64,
        resistors: Vec<MaterialResistor>,
        external_inputs: &[f64],
        controls: HystereticCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Self::Step), HystereticCircuitError> {
        self.step_candidate(target_s, resistors, external_inputs, controls, cancelled)
            .map(SharedHystereticCircuitCandidate::into_parts)
    }
}

impl sealed::Step for SharedHystereticCircuitStep {}
impl HystereticThermalStep for SharedHystereticCircuitStep {
    fn electrical(&self) -> &MaterialCircuitStep {
        &self.circuit
    }
    fn core_heat_j(&self) -> f64 {
        self.magnetic.core.core_heat_j
    }
    fn heat_j(&self) -> f64 {
        self.heat_j
    }
    fn numerical_dissipation_j(&self) -> f64 {
        self.numerical_dissipation_j
    }
    fn stored_change_j(&self) -> f64 {
        self.stored_change_j
    }
    fn supplied_work_j(&self) -> f64 {
        self.supplied_work_j
    }
    fn circuit_trials(&self) -> usize {
        self.circuit_trials
    }
}
