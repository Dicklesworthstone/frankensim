//! Shared statically selected native free motion for pair and graph constraints.

use crate::spatial_loads::{PreparedSpatialLoads, SpatialStepReceipt};
use crate::{
    DynamicsError, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt, Wrench,
};

pub(super) trait FreeFlow {
    type Receipt: Copy;

    fn len(&self) -> usize;
    fn torque_bound(&self, index: usize) -> f64;
    fn step(
        &self,
        index: usize,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        properties: MassProperties,
        duration_s: f64,
    ) -> Result<Self::Receipt, DynamicsError>;
    fn state_after(receipt: Self::Receipt) -> RigidBodyState;
}

impl FreeFlow for &[Wrench] {
    type Receipt = StepReceipt;

    fn len(&self) -> usize {
        (*self).len()
    }
    fn torque_bound(&self, index: usize) -> f64 {
        let v = self[index].torque_body;
        v.x.hypot(v.y).hypot(v.z)
    }
    fn step(
        &self,
        index: usize,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        properties: MassProperties,
        duration_s: f64,
    ) -> Result<StepReceipt, DynamicsError> {
        integrator.step(state, properties, self[index], duration_s)
    }
    fn state_after(receipt: StepReceipt) -> RigidBodyState {
        receipt.state_after
    }
}

// Bounds are prepared by the enclosing admitted transaction once, using stack
// storage for a pair or one budgeted bank for an arbitrary body slice.
pub(super) struct SpatialFlow<'a> {
    pub loads: &'a [PreparedSpatialLoads],
    pub bounds: &'a [f64],
}

impl FreeFlow for SpatialFlow<'_> {
    type Receipt = SpatialStepReceipt;

    fn len(&self) -> usize {
        self.loads.len()
    }
    fn torque_bound(&self, index: usize) -> f64 {
        self.bounds[index]
    }
    fn step(
        &self,
        index: usize,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        properties: MassProperties,
        duration_s: f64,
    ) -> Result<SpatialStepReceipt, DynamicsError> {
        integrator.step_spatial(state, properties, &self.loads[index], duration_s)
    }
    fn state_after(receipt: SpatialStepReceipt) -> RigidBodyState {
        receipt.step.state_after
    }
}
