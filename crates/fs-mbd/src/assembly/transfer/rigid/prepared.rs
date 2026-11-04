//! Transfer an existing rigid payload without discarding its intrinsic spin.
//!
//! The prepared partition retains the donor and payload as finite rigid groups.
//! One native release applies their equal/opposite separation impulses; one
//! native latch attaches the whole payload to the receiver at unchanged poses.
//! The native owners retain full tensors, shared recoil and spin plus orbital
//! momentum. No retained group is disassembled and inelastically reconstructed.

use crate::assembly::partition::AssemblyPartition;
use crate::assembly::transfer::{
    PointMassTransferControls, attachment_arm_mismatch, combined_balance,
};
use crate::assembly::{
    AssemblyBalance, AssemblyControls, AssemblyError, AssemblyLatch, RigidAssembly,
    SeparationImpulse, maximum, poll, vector,
};
use crate::{MassProperties, RigidBodyState, Vec3};

/// Cached original donor geometry and one local receiving attachment.
///
/// Group zero is the nonempty retained donor; group one is the nonempty payload.
/// Both retain the original donor's axes and every original component exactly
/// once. Preparing O(parts) geometry avoids repeating tensor assembly for the
/// two outgoing groups at each event. The receiver remains one finite body.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedRigidPayloadTransfer {
    partition: AssemblyPartition,
    donor_properties: MassProperties,
    donor_attachment_body_m: Vec3,
    payload_attachment_body_m: Vec3,
    receiver_attachment_body_m: Vec3,
}

/// The two surviving endpoints of one admitted rigid handoff.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidPayloadTransferResult {
    /// Retained donor at its new COM, in the original donor's physical axes.
    pub donor_state: RigidBodyState,
    /// Retained donor mass and full central inertia.
    pub donor_properties: MassProperties,
    /// Retained COM relative to the original donor COM, in original axes [m].
    pub donor_center_shift_body_m: Vec3,
    /// Receiver plus captured payload. Its assembly parts are the original
    /// receiver (index zero) and whole payload (index one); it uses the original
    /// receiver orientation. The prepared partition retains payload leaf geometry.
    pub receiver: AssemblyLatch,
    /// Nonnegative work supplied by the complete separation impulse set [J].
    /// Debit this source independently of the receiver's sticking-energy loss.
    pub supplied_work_j: f64,
    /// Maximum absolute attachment mismatch component actually checked [m].
    pub maximum_attachment_mismatch_m: f64,
    /// Final minus initial combined stored mass [kg].
    pub mass_residual_kg: f64,
    /// Common world origin for the angular residual: original donor COM [m].
    pub balance_origin_world_m: Vec3,
    /// Actual combined momentum and common-translating-frame energy balance.
    /// Energy residual is delta kinetic - separation work + sticking loss.
    pub balance: AssemblyBalance,
}

impl PreparedRigidPayloadTransfer {
    /// Select a nonempty proper subset of original donor components as payload.
    /// Selection order does not change original component order within groups.
    /// Repeated or out-of-range indices refuse before preparing any partition.
    ///
    /// Both attachment arms use their ORIGINAL body's COM and axes. The donor
    /// arm identifies the point on its payload that meets the receiver arm.
    /// Connectivity, contact discovery, overlap, strength and inventory identity
    /// remain caller-owned. This preparation changes neither body's state.
    /// The part cap counts all original donor components plus the receiver.
    pub fn new(
        donor: &RigidAssembly,
        payload_parts: &[usize],
        donor_attachment_body_m: Vec3,
        receiver_attachment_body_m: Vec3,
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AssemblyError> {
        poll(&mut cancelled)?;
        let count = donor.parts().len();
        let parts = count
            .checked_add(1)
            .ok_or(AssemblyError::Invalid("rigid transfer part count overflow"))?;
        controls.validate(parts, 0)?;
        vector(donor_attachment_body_m)?;
        vector(receiver_attachment_body_m)?;
        if payload_parts.is_empty() || payload_parts.len() >= count {
            return Err(AssemblyError::Invalid(
                "nonempty payload and retained donor required",
            ));
        }
        let mut membership = vec![0; count];
        for &index in payload_parts {
            poll(&mut cancelled)?;
            let entry = membership.get_mut(index).ok_or(AssemblyError::Invalid(
                "payload endpoint is not an original donor part",
            ))?;
            if *entry != 0 {
                return Err(AssemblyError::Invalid("duplicate payload component"));
            }
            *entry = 1;
        }
        let partition = AssemblyPartition::new(donor, &membership, controls, &mut cancelled)?;
        let payload_attachment_body_m = vector(
            donor_attachment_body_m.sub(partition.groups()[1].center_of_mass_in_reference()),
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            partition,
            donor_properties: donor.properties(),
            donor_attachment_body_m,
            payload_attachment_body_m,
            receiver_attachment_body_m,
        })
    }

    /// Original component inventory and cached geometry: retained donor first,
    /// whole payload second. Original indices remain available for later events.
    pub fn partition(&self) -> &AssemblyPartition {
        &self.partition
    }

    /// Perform one instantaneous handoff from the ORIGINAL donor snapshot.
    ///
    /// Separation impulses use original donor component indices and parent-COM
    /// points; every pair must cross the retained-donor/payload interface. Their
    /// combined native recoil acts before the perfectly inelastic receiving
    /// latch. The latch retains payload intrinsic spin and full inertia, and
    /// reports all relative kinetic energy removed without assigning it to heat.
    ///
    /// The shared transfer envelope counts every original donor component plus
    /// the receiver, every separation descriptor (including zeros), and one
    /// latch impulse slot. Both original local points and the actual released
    /// payload's represented points must coincide. There is no free-flight or
    /// finite-duration grasp model. Per event work/storage is O(impulses), with
    /// only two finite bodies in each native solve and no dense multibody matrix.
    ///
    /// Actual endpoints must pass native phase and combined mass, momentum and
    /// energy gates. A refusal or cancellation publishes neither endpoint. The
    /// inputs and prepared inventory remain unchanged; no free payload state
    /// is returned alongside the receiving assembly.
    pub fn transfer(
        &self,
        donor: &RigidBodyState,
        receiver: &RigidBodyState,
        receiver_properties: MassProperties,
        impulses: &[SeparationImpulse],
        controls: PointMassTransferControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidPayloadTransferResult, AssemblyError> {
        poll(&mut cancelled)?;
        let parts = self
            .partition
            .membership()
            .len()
            .checked_add(1)
            .ok_or(AssemblyError::Invalid("rigid transfer part count overflow"))?;
        let actions = impulses.len().checked_add(1).ok_or(AssemblyError::Invalid(
            "rigid transfer impulse count overflow",
        ))?;
        controls.assembly.validate(parts, actions)?;
        if [controls.attachment_tolerance_m, controls.mass_tolerance_kg]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(AssemblyError::Invalid(
                "finite nonnegative transfer tolerances required",
            ));
        }
        let mut mismatch = attachment_arm_mismatch(
            donor,
            self.donor_attachment_body_m,
            receiver,
            self.receiver_attachment_body_m,
        )?;
        if mismatch > controls.attachment_tolerance_m {
            return Err(AssemblyError::Invalid(
                "noncoincident rigid transfer attachment points",
            ));
        }
        let released =
            self.partition
                .release(*donor, impulses, controls.assembly, &mut cancelled)?;
        let retained = &self.partition.groups()[0];
        let payload = &self.partition.groups()[1];
        let payload_state = released.bodies[1];
        mismatch = mismatch.max(attachment_arm_mismatch(
            &payload_state,
            self.payload_attachment_body_m,
            receiver,
            self.receiver_attachment_body_m,
        )?);
        let receiver_point = receiver
            .pose()
            .point_world_from_body(self.receiver_attachment_body_m)?;
        let payload_point = payload_state
            .pose()
            .point_world_from_body(self.payload_attachment_body_m)?;
        mismatch = mismatch.max(maximum(vector(receiver_point.sub(payload_point))?));
        if mismatch > controls.attachment_tolerance_m {
            return Err(AssemblyError::Invalid(
                "unrepresented rigid transfer attachment coincidence",
            ));
        }
        let captured = RigidAssembly::latch(
            &[*receiver, payload_state],
            &[receiver_properties, payload.properties()],
            receiver.pose().orientation(),
            controls.assembly,
            &mut cancelled,
        )?;
        let origin = donor.pose().position_world();
        let (mass_residual, balance) = combined_balance(
            [
                (*donor, self.donor_properties),
                (*receiver, receiver_properties),
            ],
            [
                (released.bodies[0], retained.properties()),
                (captured.state, captured.assembly.properties()),
            ],
            origin,
            donor.center_of_mass_velocity_world(self.donor_properties)?,
            released.supplied_work_j,
            captured.removed_relative_energy_j,
            &mut cancelled,
        )?;
        if mass_residual.abs() > controls.mass_tolerance_kg {
            return Err(AssemblyError::Invalid(
                "combined transferred mass balance exceeds tolerance",
            ));
        }
        controls.assembly.admit(balance)?;
        poll(&mut cancelled)?;
        Ok(RigidPayloadTransferResult {
            donor_state: released.bodies[0],
            donor_properties: retained.properties(),
            donor_center_shift_body_m: retained.center_of_mass_in_reference(),
            receiver: captured,
            supplied_work_j: released.supplied_work_j,
            maximum_attachment_mismatch_m: mismatch,
            mass_residual_kg: mass_residual,
            balance_origin_world_m: origin,
            balance,
        })
    }
}
