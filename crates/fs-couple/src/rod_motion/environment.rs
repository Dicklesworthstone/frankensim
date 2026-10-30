//! Prescribed fluid refresh without resetting a native rod trajectory.
//!
//! The external reservoir supplies a new uniform hydrostatic state. At the
//! unchanged material pose its changed potential is external parameter work,
//! with zero impulse and no invented heat. Subsequent intervals use the same
//! native force, rotation, support and conservation owners as before.

use std::sync::Arc;

use fs_exec::CancelGate;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use super::{RodMotion, RodMotionError, RodMotionInertia, checkpoint, finite};

/// One accepted-time change of a prescribed ambient fluid, without a time step.
#[derive(Clone, Debug, PartialEq)]
pub struct RodFluidRefresh {
    /// Unchanged physical clock [s]; a refresh consumes no integration step.
    pub time_s: f64,
    /// Previous coherent material transport snapshot.
    pub before: FluidTransportState,
    /// Newly admitted coherent material transport snapshot.
    pub after: FluidTransportState,
    /// Previous hydrostatic potential at the actual dynamic pose [J].
    pub buoyancy_before_j: f64,
    /// New hydrostatic potential at the identical pose [J].
    pub buoyancy_after_j: f64,
    /// `buoyancy_after_j - buoyancy_before_j`, external parameter work [J].
    /// Positive means energy enters rod-plus-hydrostatic-potential storage.
    pub environmental_work_j: f64,
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Refresh a previously bound fully submerged environment at the accepted
    /// clock, retaining pose, linear/angular momenta, source inventory, stationary
    /// supports and step count. Fresh native forces replace every stale load.
    ///
    /// The potential reference remains the original unloaded rod pose. Its
    /// fixed-pose change is reported as work by the externally prescribed fluid
    /// parameter, separately from subsequent mechanical work. It is not heat or
    /// the thermodynamic cost of replacing a finite reservoir. This uniform,
    /// frozen-displaced-volume model does not acquire drag, added mass, partial
    /// immersion, free-surface flow or a compression transient through refresh.
    ///
    /// # Errors
    /// Absent fluid binding, invalid material transport, unsupported native pose,
    /// unrepresentable work or cancellation leaves the original state intact.
    pub fn rebind_submerged_fluid_candidate(
        &self,
        fluid: &impl FluidTransport,
        gate: &CancelGate,
    ) -> Result<(Self, RodFluidRefresh), RodMotionError> {
        checkpoint(gate)?;
        let model = self.model.rebind_submerged_fluid(fluid, gate)?;
        // Use the actual moving pose, never the source's construction pose.
        let old_loads = self.model.body_loads(&self.pose, gate)?;
        let new_loads = model.body_loads(&self.pose, gate)?;
        let response = model.evaluate(&self.pose, gate)?;
        let environmental_work_j = new_loads.buoyancy_potential_j - old_loads.buoyancy_potential_j;
        finite(&[environmental_work_j])?;
        let before = self
            .model
            .source()
            .submerged_fluid_load()
            .expect("refresh admitted an original fluid binding")
            .transport_state();
        let after = model
            .source()
            .submerged_fluid_load()
            .expect("refresh retains its admitted fluid binding")
            .transport_state();
        let mut candidate = self.clone();
        candidate.model = Arc::new(model);
        candidate.response = response;
        let report = RodFluidRefresh {
            time_s: self.time_s,
            before,
            after,
            buoyancy_before_j: old_loads.buoyancy_potential_j,
            buoyancy_after_j: new_loads.buoyancy_potential_j,
            environmental_work_j,
        };
        checkpoint(gate)?;
        Ok((candidate, report))
    }
}
