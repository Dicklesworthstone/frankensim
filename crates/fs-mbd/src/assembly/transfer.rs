//! Atomic point-material transfer between two finite moving rigid bodies.
//!
//! Release and capture share one instantaneous physical attachment location.
//! Both bodies recoil, both central inertias change, and canonical parcel
//! momenta cross the interface directly. No free parcel is published alongside
//! the receiving body, and no source work or sticking loss is hidden as heat.

use super::capture::{
    PointMassCaptureResult, PointMassMomentumCapture, capture_point_masses_with_momenta,
    relative_momentum,
};
use super::release::{PointMassRelease, release_point_masses};
use super::{
    AssemblyBalance, AssemblyControls, AssemblyError, Sum, VectorSum, finite, maximum, poll, vector,
};
use crate::{MassProperties, RigidBodyState, Vec3};

/// Simultaneous closed networks of material handoffs between finite bodies.
pub mod network;

/// Finite rigid-component networks retaining intrinsic spin and material geometry.
pub mod rigid;

/// One existing donor parcel transferred to one receiver material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassTransfer {
    /// Existing donor mass, original donor-COM arm and donor-frame separation
    /// impulse. The donor receives the equal/opposite recoil.
    pub release: PointMassRelease,
    /// Receiving attachment relative to its ORIGINAL COM, in receiver axes [m].
    pub receiver_attachment_arm_body_m: Vec3,
}

/// Bounded work and physical acceptance for one two-body transfer event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassTransferControls {
    /// Counts all original bodies plus EVERY descriptor, including zeros.
    /// Each positive capture and each nonzero separation impulse consumes one
    /// impulse slot. The same absolute balance gates also apply to each phase.
    pub assembly: AssemblyControls,
    /// Maximum absolute Cartesian component of attachment mismatch [m].
    /// Both the origin-relative geometry and represented world points must pass.
    pub attachment_tolerance_m: f64,
    /// Maximum actual combined mass-inventory residual [kg].
    pub mass_tolerance_kg: f64,
}

/// Both admitted endpoints, with separate work-source and loss ownership.
#[derive(Clone, Debug, PartialEq)]
pub struct PointMassTransferResult {
    /// Retained donor at its new COM, with original body orientation.
    pub donor_state: RigidBodyState,
    /// Retained donor mass and full central inertia.
    pub donor_properties: MassProperties,
    /// Donor COM shift in original donor axes [m]; rebase retained material arms.
    pub donor_center_shift_body_m: Vec3,
    /// Receiving native endpoint, aggregate recoil and nonnegative sticking loss.
    pub receiver: PointMassCaptureResult,
    /// Nonnegative work supplied by the donor's separation impulses [J].
    /// A source must debit this; receiver loss is a separate thermal destination.
    pub supplied_work_j: f64,
    /// Maximum absolute attachment mismatch component actually checked [m].
    pub maximum_attachment_mismatch_m: f64,
    /// Final minus initial combined stored mass [kg].
    pub mass_residual_kg: f64,
    /// Common world origin for the angular residual: the original donor COM [m].
    pub balance_origin_world_m: Vec3,
    /// Actual combined momentum and common-translating-frame kinetic balance.
    /// Energy residual is delta kinetic - supplied work + receiver sticking loss.
    pub balance: AssemblyBalance,
}

/// Transfer all parcels from the original donor snapshot to the original
/// receiver snapshot in one atomic physical event.
///
/// Every active donor point must coincide with its receiving point within the
/// declared tolerance. This is instantaneous local handoff, not transport across
/// a gap: free flight, collision detection, material identity, inventory geometry
/// and thermal/source debits remain caller-owned. All parcels use their SAME
/// original body frames; no sequential COM rebasing is needed within the batch.
///
/// The native release and capture owners each admit one full central inertia
/// and one shared recoil. Their canonical momenta are handed over directly.
/// O(N) compensated passes and bounded O(N) staging replace per-parcel solves.
/// Both actual endpoints must pass combined mass, linear/angular momentum and
/// kinetic/work/loss gates. Cancellation and all refusals return no endpoints
/// and leave both borrowed bodies unchanged.
#[allow(clippy::too_many_arguments)]
pub fn transfer_point_masses(
    donor: &RigidBodyState,
    donor_properties: MassProperties,
    receiver: &RigidBodyState,
    receiver_properties: MassProperties,
    transfers: &[PointMassTransfer],
    controls: PointMassTransferControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassTransferResult, AssemblyError> {
    poll(&mut cancelled)?;
    let parts = transfers
        .len()
        .checked_add(2)
        .ok_or(AssemblyError::Invalid("transfer part count overflow"))?;
    controls.assembly.validate(parts, 0)?;
    if [controls.attachment_tolerance_m, controls.mass_tolerance_kg]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(AssemblyError::Invalid(
            "finite nonnegative transfer tolerances required",
        ));
    }
    let mut impulses = 0_usize;
    let mut active = 0;
    let mut mismatch = 0.0_f64;
    for transfer in transfers {
        poll(&mut cancelled)?;
        let parcel = transfer.release;
        if !parcel.mass_kg.is_finite() || parcel.mass_kg < 0.0 {
            return Err(AssemblyError::Invalid(
                "finite nonnegative transferred mass required",
            ));
        }
        vector(parcel.attachment_arm_body_m)?;
        vector(parcel.impulse_on_parcel_body_ns)?;
        vector(transfer.receiver_attachment_arm_body_m)?;
        if parcel.mass_kg == 0.0 {
            if parcel.impulse_on_parcel_body_ns != Vec3::ZERO {
                return Err(AssemblyError::Invalid(
                    "zero transferred mass requires zero impulse",
                ));
            }
            continue;
        }
        active += 1;
        impulses = impulses
            .checked_add(1 + usize::from(parcel.impulse_on_parcel_body_ns != Vec3::ZERO))
            .ok_or(AssemblyError::Invalid("transfer impulse count overflow"))?;
        mismatch = mismatch.max(attachment_mismatch(donor, receiver, transfer)?);
    }
    controls.assembly.validate(parts, impulses)?;
    if mismatch > controls.attachment_tolerance_m {
        return Err(AssemblyError::Invalid(
            "noncoincident transfer attachment points",
        ));
    }

    // Release staging is dropped before incoming staging is allocated.
    let released = {
        let mut outgoing = Vec::with_capacity(transfers.len());
        for transfer in transfers {
            poll(&mut cancelled)?;
            outgoing.push(transfer.release);
        }
        release_point_masses(
            donor,
            donor_properties,
            &outgoing,
            controls.assembly,
            &mut cancelled,
        )?
    };
    let captured = {
        let mut incoming = Vec::with_capacity(active);
        let mut parcels = released.parcels.iter();
        for transfer in transfers {
            poll(&mut cancelled)?;
            if transfer.release.mass_kg == 0.0 {
                continue;
            }
            let parcel = parcels.next().ok_or(AssemblyError::Invalid(
                "native release parcel axis mismatch",
            ))?;
            let point = receiver
                .pose()
                .point_world_from_body(transfer.receiver_attachment_arm_body_m)?;
            mismatch = mismatch.max(maximum(vector(point.sub(parcel.position_world_m))?));
            if mismatch > controls.attachment_tolerance_m {
                return Err(AssemblyError::Invalid(
                    "unrepresented transfer attachment coincidence",
                ));
            }
            incoming.push(PointMassMomentumCapture {
                mass_kg: parcel.mass_kg,
                attachment_arm_body_m: transfer.receiver_attachment_arm_body_m,
                linear_momentum_world_ns: parcel.linear_momentum_world_ns,
            });
        }
        if parcels.next().is_some() {
            return Err(AssemblyError::Invalid(
                "native release parcel axis mismatch",
            ));
        }
        capture_point_masses_with_momenta(
            receiver,
            receiver_properties,
            &incoming,
            controls.assembly,
            &mut cancelled,
        )?
    };

    let origin = donor.pose().position_world();
    let (mass_residual, balance) = if active == 0 {
        (
            0.0,
            AssemblyBalance {
                linear_momentum_residual_world_ns: Vec3::ZERO,
                angular_momentum_residual_world_nms: Vec3::ZERO,
                energy_residual_j: 0.0,
            },
        )
    } else {
        combined_balance(
            [(*donor, donor_properties), (*receiver, receiver_properties)],
            [
                (released.state, released.properties),
                (captured.state, captured.properties),
            ],
            origin,
            donor.center_of_mass_velocity_world(donor_properties)?,
            released.supplied_work_j,
            captured.removed_relative_energy_j,
            &mut cancelled,
        )?
    };
    if mass_residual.abs() > controls.mass_tolerance_kg {
        return Err(AssemblyError::Invalid(
            "combined transferred mass balance exceeds tolerance",
        ));
    }
    controls.assembly.admit(balance)?;
    poll(&mut cancelled)?;
    Ok(PointMassTransferResult {
        donor_state: released.state,
        donor_properties: released.properties,
        donor_center_shift_body_m: released.center_shift_body_m,
        receiver: captured,
        supplied_work_j: released.supplied_work_j,
        maximum_attachment_mismatch_m: mismatch,
        mass_residual_kg: mass_residual,
        balance_origin_world_m: origin,
        balance,
    })
}

#[allow(clippy::too_many_arguments)]
fn combined_balance(
    before: impl IntoIterator<Item = (RigidBodyState, MassProperties)>,
    after: impl IntoIterator<Item = (RigidBodyState, MassProperties)>,
    origin: Vec3,
    boost: Vec3,
    supplied: f64,
    removed: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, AssemblyBalance), AssemblyError> {
    let mut mass = Sum::default();
    let mut linear = VectorSum::default();
    let mut angular = VectorSum::default();
    let mut first_moment = VectorSum::default();
    let mut energy = Sum::default();
    let mut before = before.into_iter();
    let mut after = after.into_iter();
    loop {
        poll(cancelled)?;
        let (old, old_properties, new, new_properties) = match (before.next(), after.next()) {
            (Some((old, old_properties)), Some((new, new_properties))) => {
                (old, old_properties, new, new_properties)
            }
            (None, None) => break,
            _ => {
                return Err(AssemblyError::Invalid(
                    "matching transfer endpoint axes required",
                ));
            }
        };
        if old == new && old_properties == new_properties {
            continue;
        }
        let m0 = old_properties.mass();
        let m1 = new_properties.mass();
        let dm = finite(m1 - m0)?;
        mass.add(dm);
        linear.add(new.linear_momentum_world());
        linear.add(old.linear_momentum_world().scale(-1.0));
        let p0 = relative_momentum(old.linear_momentum_world(), m0, boost)?;
        let p1 = relative_momentum(new.linear_momentum_world(), m1, boost)?;
        let r0 = vector(old.pose().position_world().sub(origin))?;
        let shift = vector(new.pose().position_world().sub(old.pose().position_world()))?;
        angular.add(vector(shift.cross(p1))?);
        angular.add(vector(r0.cross(vector(p1.sub(p0))?))?);
        first_moment.add(vector(shift.scale(m1))?);
        first_moment.add(vector(r0.scale(dm))?);
        let l0 = old.angular_momentum_body();
        let l1 = new.angular_momentum_body();
        angular.add(vector(
            old.pose()
                .orientation()
                .rotate_body_to_world(vector(l1.sub(l0))?),
        )?);
        let w0 = old_properties.angular_velocity_body_checked(l0)?;
        let w1 = new_properties.angular_velocity_body_checked(l1)?;
        let v0 = vector(Vec3::new(p0.x / m0, p0.y / m0, p0.z / m0))?;
        let v1 = vector(Vec3::new(p1.x / m1, p1.y / m1, p1.z / m1))?;
        // Polarized changes avoid subtracting large common spin or kinetic
        // totals. Both states retain their own original body orientation.
        energy.add(finite(0.5 * p1.sub(p0).dot(v1))?);
        energy.add(finite(0.5 * p0.dot(v1.sub(v0)))?);
        energy.add(finite(0.5 * l1.sub(l0).dot(w1))?);
        energy.add(finite(0.5 * l0.dot(w1.sub(w0)))?);
    }
    // Actual stored COM and mass changes also contribute under a common boost.
    angular.add(vector(vector(first_moment.value())?.cross(boost))?);
    energy.add(-supplied);
    energy.add(removed);
    Ok((
        finite(mass.value())?,
        AssemblyBalance {
            linear_momentum_residual_world_ns: vector(linear.value())?,
            angular_momentum_residual_world_nms: vector(angular.value())?,
            energy_residual_j: finite(energy.value())?,
        },
    ))
}

// Keep attachment admission common to pair and network handoffs. In particular,
// rounded equality of large world positions must not hide different local points.
fn attachment_mismatch(
    donor: &RigidBodyState,
    receiver: &RigidBodyState,
    transfer: &PointMassTransfer,
) -> Result<f64, AssemblyError> {
    let mut gap = VectorSum::default();
    gap.add(vector(
        receiver
            .pose()
            .position_world()
            .sub(donor.pose().position_world()),
    )?);
    gap.add(vector(
        receiver
            .pose()
            .orientation()
            .rotate_body_to_world(transfer.receiver_attachment_arm_body_m),
    )?);
    gap.add(vector(
        donor
            .pose()
            .orientation()
            .rotate_body_to_world(transfer.release.attachment_arm_body_m)
            .scale(-1.0),
    )?);
    Ok(maximum(vector(gap.value())?))
}
