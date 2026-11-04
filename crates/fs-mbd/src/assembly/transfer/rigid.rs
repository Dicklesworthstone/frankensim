//! Simultaneous transfer of finite rigid components between moving assemblies.
//!
//! Original components going from one donor to one receiver stay rigid as a
//! group. Each donor releases those groups once, with the native shared recoil.
//! Only after all releases do receivers latch their incoming groups. Cycles
//! cannot re-export new arrivals, and retained groups incur no fictitious
//! leaf-level release or sticking loss. Full intrinsic spin, orbital momentum
//! and changing central inertia participate in both original native owners.

use super::combined_balance;
use crate::assembly::partition::AssemblyPartition;
use crate::assembly::{
    AssemblyBalance, AssemblyControls, AssemblyError, AssemblyPart, RigidAssembly,
    SeparationImpulse, Sum, compose, finite, poll,
};
use crate::{Pose, RigidBodyState, Vec3};

/// One original rigid inventory and its authored simultaneous redistribution.
#[derive(Clone, Copy, Debug)]
pub struct RigidTransferInput<'a> {
    /// Original full-inertia inventory; its parts define material identity.
    pub assembly: &'a RigidAssembly,
    /// Original COM state, with axes matching the assembly.
    pub state: RigidBodyState,
    /// One destination body index for EVERY original component. Components with
    /// the same donor and destination stay rigid throughout the release.
    pub destinations: &'a [usize],
    /// Internal force/couple pairs between original component indices.
    /// Points and impulses use the original donor COM frame. Each pair must
    /// cross two different destination groups.
    pub impulses: &'a [SeparationImpulse],
}

/// Original identity of one component in a returned body's local part order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RigidComponentOrigin {
    /// Original input body index.
    pub body: usize,
    /// Original local component index.
    pub part: usize,
}

/// Aggregate allocation/work bounds and actual physical balance acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidTransferControls {
    /// The part cap counts all original bodies PLUS all original components;
    /// the impulse cap counts every authored separation pair. The same native
    /// momentum/work tolerances apply to each release and latch.
    pub assembly: AssemblyControls,
    /// Maximum actual final minus original stored mass magnitude [kg].
    pub mass_tolerance_kg: f64,
}

/// One final inventory in the original receiver's axes, at its new COM.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidTransferEndpoint {
    /// Flat ORIGINAL component geometry, preserving full properties and world
    /// poses. Parts are ordered by donor index, then original donor part index.
    /// This inventory can be redistributed again without losing leaf identity.
    pub assembly: RigidAssembly,
    /// Actual native momentum-preserving endpoint.
    pub state: RigidBodyState,
    /// Source identity for every corresponding final assembly part.
    pub origins: Vec<RigidComponentOrigin>,
    /// Work supplied by this ORIGINAL body's simultaneous separation [J].
    pub supplied_work_j: f64,
    /// Relative kinetic energy removed by THIS receiver's joint latch [J].
    pub removed_relative_energy_j: f64,
    /// Native outgoing group-release balance; zero when no separation occurs.
    pub release_balance: AssemblyBalance,
    /// Native incoming group-latch balance; zero when no latch is needed.
    pub capture_balance: AssemblyBalance,
}

/// All final bodies, with one closed-event mass/momentum/work account.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidTransferResult {
    /// Complete endpoints in original body-index order; none are empty.
    pub bodies: Vec<RigidTransferEndpoint>,
    /// Sum of simultaneous per-donor source work [J].
    pub supplied_work_j: f64,
    /// Sum of simultaneous per-receiver sticking losses [J].
    pub removed_relative_energy_j: f64,
    /// Actual final minus original total mass [kg].
    pub mass_residual_kg: f64,
    /// Common physical origin for the angular-momentum residual [m].
    pub balance_origin_world_m: Vec3,
    /// Actual global momentum and kinetic - source work + loss residuals.
    pub balance: AssemblyBalance,
}

/// Redistribute complete rigid components in one instantaneous closed event.
///
/// Each original component has exactly one destination. Every final body must
/// contain material, but a donor may transfer its ENTIRE original inventory if
/// other original bodies replace it during this event. Original world geometry
/// never moves; the caller owns whether the authored rigid attachments are
/// physically admissible. This owner does not discover contacts, teleport a
/// payload to a receiver, integrate free flight or infer fracture/actuator laws.
///
/// Native prepared partitions preserve all within-group rigidity during
/// separation, including shared force/couple recoil and work cross terms.
/// Native latching then includes intrinsic spin and relative group motion.
/// Source work and removed kinetic energy remain distinct: the caller owns
/// actuator debits and any thermal destination.
///
/// Routing and retained leaf geometry require O(bodies + parts + impulses)
/// work/storage, in fixed donor/part order. One reusable destination map avoids
/// a body-by-body dense table; no global mobility matrix is introduced.
/// Original inventories are borrowed. Any refusal or cancellation publishes
/// no endpoint, including cancellation after the final physical balance.
pub fn transfer_rigid_components(
    inputs: &[RigidTransferInput<'_>],
    controls: RigidTransferControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RigidTransferResult, AssemblyError> {
    poll(&mut cancelled)?;
    controls.assembly.validate(inputs.len(), 0)?;
    if !controls.mass_tolerance_kg.is_finite() || controls.mass_tolerance_kg < 0.0 {
        return Err(AssemblyError::Invalid(
            "finite nonnegative rigid transfer mass tolerance required",
        ));
    }
    let mut count = inputs.len();
    let mut impulse_count = 0_usize;
    for input in inputs {
        poll(&mut cancelled)?;
        let parts = input.assembly.parts().len();
        if parts != input.destinations.len() {
            return Err(AssemblyError::Invalid(
                "one destination per original rigid component required",
            ));
        }
        count = count.checked_add(parts).ok_or(AssemblyError::Invalid(
            "rigid transfer component count overflow",
        ))?;
        impulse_count =
            impulse_count
                .checked_add(input.impulses.len())
                .ok_or(AssemblyError::Invalid(
                    "rigid transfer impulse count overflow",
                ))?;
        input.assembly.properties().validate()?;
        input
            .state
            .center_of_mass_velocity_world(input.assembly.properties())?;
        input
            .assembly
            .properties()
            .angular_velocity_body_checked(input.state.angular_momentum_body())?;
    }
    // Admit aggregate bounds before allocating any network staging.
    controls.assembly.validate(count, impulse_count)?;
    let mut final_sizes = vec![0_usize; inputs.len()];
    let mut active = vec![false; inputs.len()];
    for (donor, input) in inputs.iter().enumerate() {
        for &destination in input.destinations {
            poll(&mut cancelled)?;
            if destination >= inputs.len() {
                return Err(AssemblyError::Invalid(
                    "rigid component destination is not an original body",
                ));
            }
            final_sizes[destination] += 1;
            if destination != donor {
                active[donor] = true;
                active[destination] = true;
            }
        }
        for impulse in input.impulses {
            poll(&mut cancelled)?;
            let Some((&a, &b)) = input
                .destinations
                .get(impulse.part_a)
                .zip(input.destinations.get(impulse.part_b))
            else {
                return Err(AssemblyError::Invalid(
                    "separation endpoint is not an original component",
                ));
            };
            if a == b
                || !impulse.point_assembly_m.is_finite()
                || !impulse.impulse_on_a_assembly_ns.is_finite()
                || !impulse.couple_on_a_assembly_nms.is_finite()
            {
                return Err(AssemblyError::Invalid(
                    "finite separation impulse must cross rigid destination groups",
                ));
            }
        }
    }
    if final_sizes.contains(&0) {
        return Err(AssemblyError::Invalid(
            "every final rigid transfer body must contain material",
        ));
    }

    let mut incoming: Vec<Vec<Packet>> = (0..inputs.len()).map(|_| Vec::new()).collect();
    let mut source_work = vec![0.0; inputs.len()];
    let mut release_balances = vec![zero_balance(); inputs.len()];
    let mut group_for_destination = vec![usize::MAX; inputs.len()];
    let mut membership = Vec::new();
    let mut destinations = Vec::new();
    let mut supplied = Sum::default();
    for (donor, input) in inputs.iter().enumerate() {
        poll(&mut cancelled)?;
        if !active[donor] {
            continue;
        }
        membership.clear();
        destinations.clear();
        for &destination in input.destinations {
            poll(&mut cancelled)?;
            let group = &mut group_for_destination[destination];
            if *group == usize::MAX {
                *group = destinations.len();
                destinations.push(destination);
            }
            membership.push(*group);
        }
        // One group undergoes no actual separation. Preserve its exact native
        // inventory/state, even when the whole assembly changes ownership.
        if destinations.len() == 1 {
            let mut origins = Vec::with_capacity(input.assembly.parts().len());
            for part in 0..input.assembly.parts().len() {
                poll(&mut cancelled)?;
                origins.push(RigidComponentOrigin { body: donor, part });
            }
            incoming[destinations[0]].push(Packet {
                assembly: input.assembly.clone(),
                state: input.state,
                origins,
            });
        } else {
            let partition = AssemblyPartition::new(
                input.assembly,
                &membership,
                controls.assembly,
                &mut cancelled,
            )?;
            let released = partition.release(
                input.state,
                input.impulses,
                controls.assembly,
                &mut cancelled,
            )?;
            source_work[donor] = released.supplied_work_j;
            release_balances[donor] = released.balance;
            supplied.add(released.supplied_work_j);
            for (group, (&destination, &state)) in
                destinations.iter().zip(&released.bodies).enumerate()
            {
                poll(&mut cancelled)?;
                let indices =
                    partition
                        .original_part_indices(group)
                        .ok_or(AssemblyError::Invalid(
                            "native rigid group identity missing",
                        ))?;
                let mut origins = Vec::with_capacity(indices.len());
                for &part in indices {
                    poll(&mut cancelled)?;
                    origins.push(RigidComponentOrigin { body: donor, part });
                }
                incoming[destination].push(Packet {
                    assembly: partition.groups()[group].clone(),
                    state,
                    origins,
                });
            }
        }
        for &destination in &destinations {
            poll(&mut cancelled)?;
            group_for_destination[destination] = usize::MAX;
        }
    }

    let mut endpoints = Vec::with_capacity(inputs.len());
    let mut states = Vec::new();
    let mut properties = Vec::new();
    let mut removed = Sum::default();
    for (receiver, mut packets) in incoming.into_iter().enumerate() {
        poll(&mut cancelled)?;
        let input = inputs[receiver];
        let (assembly, state, origins, loss, capture_balance) = if !active[receiver] {
            let mut origins = Vec::with_capacity(input.assembly.parts().len());
            for part in 0..input.assembly.parts().len() {
                poll(&mut cancelled)?;
                origins.push(RigidComponentOrigin {
                    body: receiver,
                    part,
                });
            }
            (
                input.assembly.clone(),
                input.state,
                origins,
                0.0,
                zero_balance(),
            )
        } else if packets.len() == 1
            && packets[0].state.pose().orientation() == input.state.pose().orientation()
        {
            let packet = packets
                .pop()
                .ok_or(AssemblyError::Invalid("missing rigid transfer packet"))?;
            (
                packet.assembly,
                packet.state,
                packet.origins,
                0.0,
                zero_balance(),
            )
        } else {
            states.clear();
            properties.clear();
            for packet in &packets {
                poll(&mut cancelled)?;
                states.push(packet.state);
                properties.push(packet.assembly.properties());
            }
            let capture = RigidAssembly::latch(
                &states,
                &properties,
                input.state.pose().orientation(),
                controls.assembly,
                &mut cancelled,
            )?;
            let mut parts = Vec::with_capacity(final_sizes[receiver]);
            let mut origins = Vec::with_capacity(final_sizes[receiver]);
            for (placed, packet) in capture.assembly.parts().iter().zip(&packets) {
                for (part, &origin) in packet.assembly.parts().iter().zip(&packet.origins) {
                    poll(&mut cancelled)?;
                    parts.push(AssemblyPart {
                        properties: part.properties,
                        pose_in_reference: Pose::new(
                            placed
                                .pose_in_reference
                                .point_world_from_body(part.pose_in_reference.position_world())?,
                            compose(
                                placed.pose_in_reference.orientation(),
                                part.pose_in_reference.orientation(),
                            )?,
                        )?,
                    });
                    origins.push(origin);
                }
            }
            // Flatten geometry only. Keep the ORIGINAL native latch inertia
            // and state; recomputing or relatching leaves would create another
            // physical event and could discard shared-group work/rigidity.
            let assembly = RigidAssembly {
                parts,
                ..capture.assembly
            };
            (
                assembly,
                capture.state,
                origins,
                capture.removed_relative_energy_j,
                capture.balance,
            )
        };
        if assembly.parts().len() != final_sizes[receiver] || origins.len() != final_sizes[receiver]
        {
            return Err(AssemblyError::Invalid(
                "rigid transfer leaf inventory mismatch",
            ));
        }
        removed.add(loss);
        endpoints.push(RigidTransferEndpoint {
            assembly,
            state,
            origins,
            supplied_work_j: source_work[receiver],
            removed_relative_energy_j: loss,
            release_balance: release_balances[receiver],
            capture_balance,
        });
    }
    let supplied = finite(supplied.value())?;
    let removed = finite(removed.value())?;
    let anchor = active.iter().position(|&a| a).unwrap_or(0);
    let origin = inputs[anchor].state.pose().position_world();
    let boost = inputs[anchor]
        .state
        .center_of_mass_velocity_world(inputs[anchor].assembly.properties())?;
    let (mass_residual, balance) = combined_balance(
        inputs.iter().map(|i| (i.state, i.assembly.properties())),
        endpoints.iter().map(|e| (e.state, e.assembly.properties())),
        origin,
        boost,
        supplied,
        removed,
        &mut cancelled,
    )?;
    if mass_residual.abs() > controls.mass_tolerance_kg {
        return Err(AssemblyError::Invalid(
            "rigid transfer stored mass balance exceeds tolerance",
        ));
    }
    controls.assembly.admit(balance)?;
    poll(&mut cancelled)?;
    Ok(RigidTransferResult {
        bodies: endpoints,
        supplied_work_j: supplied,
        removed_relative_energy_j: removed,
        mass_residual_kg: mass_residual,
        balance_origin_world_m: origin,
        balance,
    })
}

struct Packet {
    assembly: RigidAssembly,
    state: RigidBodyState,
    origins: Vec<RigidComponentOrigin>,
}

fn zero_balance() -> AssemblyBalance {
    AssemblyBalance {
        linear_momentum_residual_world_ns: Vec3::ZERO,
        angular_momentum_residual_world_nms: Vec3::ZERO,
        energy_residual_j: 0.0,
    }
}
