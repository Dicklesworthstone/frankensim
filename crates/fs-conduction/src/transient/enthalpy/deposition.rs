//! Shared constitutive and executor owners for the reduced deposition network.
//! The spatial P1 enthalpy solver remains stationary; this adapter uses its
//! existing material-chart owner, not a second phase or heat-capacity model.
//! The phase-curve implementation of the shared `ThermalLaw` lives in
//! `fs-material`, so it is available without loading this FEM solver.

/// Implicit, radiation-capable parcel transport using the existing nonlinear
/// enthalpy and conservative finite-contact solvers. Its aggregate Biot gate
/// refuses outside the isothermal regime; it is not a growing spatial P1 solve.
pub mod implicit;
/// Bond development driven by actual explicit-deposition substeps and contacts.
pub mod bonding;
/// Staggered cohesive opening/reclosure changes actual thermal contact.
pub mod fracture;

use fs_exec::Cx;

use crate::transient::deposition::{
    DepositionError, DepositionSimulation, EnergyAudit, StepConfig,
};

/// Advance within an existing executor scope. Cancellation discards the whole
/// candidate, including newly activated material and event cursor updates.
pub fn advance_to<'a>(
    cx: &Cx<'_>,
    simulation: &DepositionSimulation<'a>,
    target_s: f64,
    config: StepConfig,
) -> Result<(DepositionSimulation<'a>, EnergyAudit), DepositionError> {
    simulation.advance_to(target_s, config, || {
        cx.checkpoint().map_err(|_| DepositionError::Cancelled)
    })
}
