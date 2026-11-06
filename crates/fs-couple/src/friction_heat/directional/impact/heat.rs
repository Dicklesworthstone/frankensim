//! Same-time impact heating through the original mechanical and caloric owners.

use super::{
    SourceRigidImpactBinding, SourceRigidImpactError, SourceRigidImpactReport,
    apply_source_rigid_elliptic_impacts,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
};
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::elliptic::RigidEllipticSourceError;
use fs_mbd::frictional_impact::{
    RigidFrictionalContact, RigidImpactControls, RigidImpactKinematics,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

mod partition;
pub use partition::ImpactHeatReceipt;
use partition::{ImpactHeatPartition, compensated_sum};

/// Original mechanical admission and explicit represented-energy tolerances.
#[derive(Clone, Copy, Debug)]
pub struct ImpactHeatControls {
    pub mechanical: RigidImpactControls,
    pub cell_energy_tolerance_j: f64,
    /// Bounds summed absolute cell-storage error and heat-partition roundoff.
    pub energy_tolerance_j: f64,
    pub combined_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum ImpactHeatError<E> {
    Mechanical(RigidEllipticSourceError<SourceRigidImpactError<E>>),
    Thermal(DepositionError),
}
impl<E: core::fmt::Debug> core::fmt::Display for ImpactHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source rigid impact heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ImpactHeatError<E> {}
impl<E> From<DepositionError> for ImpactHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

#[derive(Debug)]
pub struct HeatedSourceRigidImpactReport {
    pub impact: SourceRigidImpactReport,
    pub heat: ImpactHeatReceipt,
    /// Delta K + original thermal-storage change + unallocated loss - support
    /// work. An H receiver retains its original pressure-work convention; this
    /// account does not relabel enthalpy or declared sensible storage as U.
    pub combined_energy_residual_j: f64,
}

/// Deposit explicitly declared fractions of one simultaneous event's total loss.
///
/// `heat_fractions` maps unique original cell IDs to fractions of the GLOBAL
/// native impact dissipation. Entries must be finite, nonnegative and sum to at
/// most one; the remainder stays explicitly unallocated. The loss can include
/// normal/restitution and tangential mechanisms: no contactwise attribution,
/// heat fraction, impulse pressure, event duration or contact area is inferred.
/// The original interface laws sample the incoming state once, as before.
///
/// The event occurs at the receiver's accepted time. The caller owns mechanical
/// event alignment; pending thermal events must be processed first. Original
/// cell masses, caloric identities, H/U/phase coordinates and the physical clock
/// remain with `apply_heat_j`. Boundary watts and contacts are not reintegrated.
/// Empty/zero fractions reproduce the original mechanical endpoint exactly.
/// Bodies, source histories and thermal storage publish together after the final
/// checkpoint. Callback-owned external side effects remain outside this transaction.
#[allow(clippy::too_many_arguments)]
pub fn apply_source_rigid_elliptic_impacts_with_heat<E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact<()>],
    histories: &mut [SourceThermalEllipticState],
    bindings: &[SourceRigidImpactBinding<'_>],
    thermal: &mut DepositionSimulation<'_>,
    heat_fractions: &BTreeMap<CellId, f64>,
    controls: ImpactHeatControls,
    conditions: impl FnMut(
        &SourceThermalEllipticState,
        &RigidImpactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    gate: &CancelGate,
) -> Result<HeatedSourceRigidImpactReport, ImpactHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    if bodies.len() > controls.mechanical.maximum_bodies
        || contacts.len() > controls.mechanical.maximum_contacts
        || bodies.len() != properties.len()
        || histories.len() != contacts.len()
        || bindings.len() != contacts.len()
        || !controls.combined_energy_tolerance_j.is_finite()
        || controls.combined_energy_tolerance_j < 0.0
    {
        return Err(DepositionError::InvalidInput(
            "bounded matching impact/heat inventories required",
        )
        .into());
    }
    let partition = ImpactHeatPartition::new(
        thermal,
        heat_fractions,
        controls.cell_energy_tolerance_j,
        controls.energy_tolerance_j,
        checkpoint,
    )?;
    let mut next_bodies = bodies.to_vec();
    let mut next_histories = histories.to_vec();
    let impact = apply_source_rigid_elliptic_impacts(
        &mut next_bodies,
        properties,
        contacts,
        &mut next_histories,
        bindings,
        controls.mechanical,
        conditions,
        gate,
    )
    .map_err(ImpactHeatError::Mechanical)?;
    let (next_thermal, heat) = partition.deposit(
        thermal,
        impact.mechanical.contact.dissipated_energy,
        checkpoint,
    )?;
    // Sum already admitted native defects rather than subtracting large total
    // kinetic/reference energies. No numerical defect becomes physical heat.
    let combined = compensated_sum([
        impact.mechanical.total_energy_residual_j,
        heat.thermal.residual_j,
        heat.partition_energy_residual_j,
    ]);
    if !combined.is_finite() {
        return Err(DepositionError::Arithmetic.into());
    }
    if combined.abs() > controls.combined_energy_tolerance_j {
        return Err(DepositionError::EnergyBalance(combined).into());
    }
    checkpoint()?;
    bodies.copy_from_slice(&next_bodies);
    histories.swap_with_slice(&mut next_histories);
    *thermal = next_thermal;
    Ok(HeatedSourceRigidImpactReport {
        impact,
        heat,
        combined_energy_residual_j: combined,
    })
}
