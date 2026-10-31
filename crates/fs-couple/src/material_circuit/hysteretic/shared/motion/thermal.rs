//! Source-bound finite heating on the original shared hysteretic gap trajectory.
//!
//! The generic hysteretic enthalpy owner resolves every coil's original R(T)
//! and the common core's source query at the same proposed temperatures. It
//! reruns the native circuit/magnetic/mechanical step from the accepted linkage,
//! polarization, gap and momentum. One pinning-heat receipt and each original
//! resistor loss enter their existing finite caloric receivers exactly once.
//!
//! The motion receipt already includes magnetic/gap/leakage energy, original
//! circuit storage, body kinetic energy and frozen support elasticity. Magnetic
//! force work is internal to that account. Only held external force work and
//! external electrical work supply the combined motion/thermal transaction.
//! The original MaterialNetwork and complete HeatAssembly owners retain all
//! source pins, links, finite coolant, contacts and radiation semantics.
//!
//! Core fit coefficients remain constant over their explicitly admitted source
//! temperature interval. Heating changes the currents and resulting hysteresis
//! and motion through the original electrical transport law. Body/support mass,
//! gap-medium conditions and frozen support properties retain their original
//! sources; this adds no thermal expansion or unmeasured magnetic temperature law.
//! The native motion owner still refuses thermally live supports and unresolved
//! common-flux reversals. Endpoint thermal coupling remains first order.

use super::{
    SharedHystereticGapMotion, SharedHystereticGapMotionCandidate, SharedHystereticGapMotionSolve,
    SharedHystereticGapMotionStep,
};
use crate::{
    circuit_enthalpy::material::MaterialNetworkSample,
    material_circuit::{
        MaterialCircuit, MaterialCircuitStep, MaterialResistor,
        hysteretic::{
            HystereticCircuitError, HystereticCircuitSolve,
            thermal::{
                HystereticEnthalpyCandidate, HystereticThermalHistory, HystereticThermalStep,
                ThermalHystereticCircuit, sealed,
            },
        },
    },
};
use fs_matdb::QueryPoint;
use fs_material::field_storage::hysteresis::UniformHystereticCore;

/// Per-interval physical drive and the original native motion solve policy.
/// Force is positive toward opening and is held across all thermal trials.
/// It is never stored as a replacement material coefficient or heat source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedHystereticGapHeatInput {
    pub external_force_n: f64,
    pub controls: SharedHystereticGapMotionSolve,
}

/// One accepted magnetic/circuit/body history on its original finite heat owner.
/// `circuit()` returns the complete native SharedHystereticGapMotion. All bind,
/// bulk-step and HeatAssembly methods use the existing hysteretic transaction.
pub type ThermalSharedHystereticGapMotion<'n, 'a> =
    ThermalHystereticCircuit<'n, 'a, SharedHystereticGapMotion>;

/// Publishable native motion and finite thermal histories for the same interval.
/// `report().electrical` is the original full motion receipt; its circuit, core
/// heat, force work, kinetic and support energy accounts remain distinguishable.
pub type SharedHystereticGapMotionEnthalpyCandidate<'t, 'n, 'a, S = MaterialNetworkSample> =
    HystereticEnthalpyCandidate<'t, 'n, 'a, S, SharedHystereticGapMotion>;

impl sealed::History for SharedHystereticGapMotion {}
impl HystereticThermalHistory for SharedHystereticGapMotion {
    type Step = SharedHystereticGapMotionStep;
    type Input = SharedHystereticGapHeatInput;

    fn electrical_controls(input: Self::Input) -> HystereticCircuitSolve {
        input.controls.circuit
    }
    fn with_electrical_controls(
        mut input: Self::Input,
        controls: HystereticCircuitSolve,
    ) -> Self::Input {
        input.controls.circuit = controls;
        input
    }
    fn uniform_core(&self) -> &UniformHystereticCore {
        self.circuit().core().core()
    }
    fn native_circuit(&self) -> &MaterialCircuit {
        self.circuit().native_circuit()
    }
    fn time_s(&self) -> f64 {
        SharedHystereticGapMotion::time_s(self)
    }
    fn with_core_source_query(&self, query: &QueryPoint) -> Result<Self, HystereticCircuitError> {
        SharedHystereticGapMotion::with_core_source_query(self, query)
    }
    fn step_for_heat(
        &self,
        target_s: f64,
        resistors: Vec<MaterialResistor>,
        external_inputs: &[f64],
        input: Self::Input,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Self::Step), HystereticCircuitError> {
        self.step_candidate(
            target_s,
            resistors,
            external_inputs,
            input.external_force_n,
            input.controls,
            cancelled,
        )
        .map(SharedHystereticGapMotionCandidate::into_parts)
    }
}

impl sealed::Step for SharedHystereticGapMotionStep {}
impl HystereticThermalStep for SharedHystereticGapMotionStep {
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
        // The original native circuit and conservative mechanical descriptor
        // explicitly close their work; only resistor and pinning losses heat.
        0.0
    }
    fn stored_change_j(&self) -> f64 {
        self.stored_change_j
    }
    fn supplied_work_j(&self) -> f64 {
        self.electrical_supply_j + self.mechanical_supply_j
    }
    fn circuit_trials(&self) -> usize {
        self.circuit_trials
    }
}
