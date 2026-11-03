//! Simultaneous normal momentum exchange between one plate and finite bodies.
//!
//! Every site shares the original plate's consistent mass. Rigid cross terms
//! connect only sites belonging to the same supplied body. The original contact
//! law resolves that complete matrix before any native state is published.

use super::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use super::{
    PlateBodyTransfer, PlateImpulseBalance, PlateRigidImpactError, check, components, dd, dot,
    poll, relative_speed, vector,
};
use fs_contact::multicontact::{MultiImpactResult, NormalContactImpact, solve_normal_impacts};
use fs_math::dd::Dd;
use fs_mbd::{MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};
use fs_plate::PlateMesh;
use fs_plate_transient::{PlateDynamics, PlateLoad};
use fs_sparse::Csr;
use fs_time::galpha::SecondOrderState;
use fs_time::galpha::impulse::SecondOrderImpulse;

/// One continuing finite body with its original center-of-mass inertia.
/// Index identity is supplied by the event's slice, not inferred from equal
/// poses, masses or velocities. Unaffected entries retain their state bits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateImpactBody {
    pub state: RigidBodyState,
    pub properties: MassProperties,
}

/// Body attachment paired with one prepared plate site, in original site order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateImpactBodyWitness {
    pub body: usize,
    /// Original material arm FROM THE BODY CENTER OF MASS, in its body frame.
    pub arm_body: Vec3,
}

/// The original point/mass/contact budgets plus a finite participant limit.
/// The plate workspace cap also includes all live body resultants and native
/// receipts, using their actual Rust storage size converted to scalar slots.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateRigidImpactNetworkBudget {
    pub contacts: PlateRigidImpactSetBudget,
    pub max_bodies: usize,
}

/// One globally accepted event across the original native participants.
#[derive(Debug)]
pub struct PlateRigidImpactNetworkReport {
    /// Original determinacy and work receipt. For redundant contacts the
    /// impulses are representatives, not unique local pressure or heat shares.
    pub contact: MultiImpactResult,
    /// Complete row-major mobility in the original site order [kg^-1].
    pub delassus_kg_inverse: Vec<f64>,
    /// Measured from the actually stored candidate velocities at every site,
    /// including geometrically inactive sites affected through the shared plate.
    pub relative_normal_speeds_after: Vec<f64>,
    pub plate: Option<SecondOrderImpulse>,
    /// Original body order; absent for every unaffected participant.
    pub bodies: Vec<Option<MomentumImpulseReceipt>>,
    pub balance: Option<PlateImpulseBalance>,
}

#[derive(Clone, Copy)]
struct BodyResultant {
    linear: Dd,
    angular: [Dd; 3],
}

/// Prepared consistent-mass coupling between one plate and several bodies.
/// No body is bound at preparation: the same cached plate responses serve each
/// event's explicitly indexed participants and original material arms.
pub struct PlateRigidImpactNetwork<'a, L: ?Sized> {
    prepared: PlateRigidImpactSet<'a, L>,
    max_bodies: usize,
}

impl<'a, L: PlateLoad + ?Sized> PlateRigidImpactNetwork<'a, L> {
    /// Prepare each original point once, reusing the existing full mass/support
    /// admission and cached B M^-1 B^T owner. Body and point caps are admitted
    /// before allocation; no plate-sized dense matrix or new mass solve exists.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        dynamics: PlateDynamics<'a, L>,
        mesh: &'a PlateMesh,
        full_mass: &'a Csr,
        sites: &[PlateImpactSite],
        controls: PlateRigidImpactSetControls,
        budget: PlateRigidImpactNetworkBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateRigidImpactError> {
        poll(cancelled)?;
        if budget.max_bodies == 0 {
            return Err(PlateRigidImpactError::Budget("finite body count"));
        }
        let body_components = (core::mem::size_of::<BodyResultant>()
            + core::mem::size_of::<Option<MomentumImpulseReceipt>>())
        .div_ceil(core::mem::size_of::<f64>());
        let additional = budget
            .max_bodies
            .checked_mul(body_components)
            .ok_or(PlateRigidImpactError::Budget("body workspace overflow"))?;
        let mut contact_budget = budget.contacts;
        contact_budget.plate.max_workspace_components = contact_budget
            .plate
            .max_workspace_components
            .checked_sub(additional)
            .ok_or(PlateRigidImpactError::Budget("finite body workspace"))?;
        let prepared = PlateRigidImpactSet::prepare(
            dynamics,
            mesh,
            full_mass,
            sites,
            controls,
            contact_budget,
            cancelled,
        )?;
        Ok(Self {
            prepared,
            max_bodies: budget.max_bodies,
        })
    }

    /// Cached original plate contribution B M^-1 B^T [kg^-1].
    #[must_use]
    pub fn plate_cross_mobility(&self) -> &[f64] {
        self.prepared.plate_cross_mobility()
    }

    /// Resolve all sites together, then atomically publish the native plate and
    /// every affected finite body. Each body receives one resultant impulse;
    /// the shared plate receives one original impulse/acceleration restart.
    /// Stationary support inertia, full rigid tensors, actual endpoint speeds,
    /// transverse momentum and combined kinetic work must all be admitted.
    ///
    /// Contact assembly is O(sites^2), with O(sites) original body inertia
    /// actions. Cached plate responses are reused. There is no per-body plate
    /// solve, sequential impact ordering, time advance or thermal partition.
    pub fn apply(
        &self,
        plate: &mut SecondOrderState,
        bodies: &mut [PlateImpactBody],
        witnesses: &[PlateImpactBodyWitness],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateRigidImpactNetworkReport, PlateRigidImpactError> {
        poll(cancelled)?;
        let k = self.prepared.sites.len();
        if bodies.len() > self.max_bodies {
            return Err(PlateRigidImpactError::Budget("finite body count"));
        }
        if witnesses.len() != k || witnesses.iter().any(|w| w.body >= bodies.len()) {
            return Err(PlateRigidImpactError::Invalid(
                "one valid body witness per original plate site required",
            ));
        }
        let first = &self.prepared.first;
        first.admit_state(plate)?;
        let mut contacts = Vec::with_capacity(k);
        let mut torques = Vec::with_capacity(k);
        let mut rotations = Vec::with_capacity(k);
        for (i, witness) in witnesses.iter().enumerate() {
            poll(cancelled)?;
            let participant = bodies[witness.body];
            let (point, row, _) = self.prepared.row(i);
            let sample = point
                .sample(&plate.q, &plate.v)
                .map_err(PlateRigidImpactError::Point)?;
            let motion = participant
                .state
                .point_kinematics(participant.properties, witness.arm_body)
                .map_err(PlateRigidImpactError::Body)?;
            let [x, y] = point.point();
            check(
                "contact point x",
                motion.point_world.x - x,
                first.controls.point_tolerance_m,
            )?;
            check(
                "contact point y",
                motion.point_world.y - y,
                first.controls.point_tolerance_m,
            )?;
            let site = self.prepared.sites[i];
            contacts.push(NormalContactImpact {
                gap: motion.point_world.z - sample.displacement,
                activation_tolerance: site.activation_tolerance_m,
                relative_normal_speed: relative_speed(row, &motion, &plate.v),
                restitution: site.restitution,
            });
            let normal_body = participant
                .state
                .pose()
                .orientation()
                .rotate_world_to_body(Vec3::new(0.0, 0.0, 1.0));
            let torque = witness.arm_body.cross(normal_body);
            torques.push(torque);
            rotations.push(
                participant
                    .properties
                    .angular_velocity_body_checked(torque)
                    .map_err(PlateRigidImpactError::Body)?,
            );
        }
        let mut delassus = self.prepared.plate_cross_mobility().to_vec();
        for i in 0..k {
            poll(cancelled)?;
            for j in 0..k {
                if witnesses[i].body == witnesses[j].body {
                    let mass = bodies[witnesses[i].body].properties.mass();
                    delassus[i * k + j] = (dd(delassus[i * k + j])
                        + dd(1.0) / dd(mass)
                        + dd(dot(&components(torques[i]), &components(rotations[j]))))
                    .to_f64();
                }
            }
        }
        let contact =
            solve_normal_impacts(&contacts, &delassus, self.prepared.contact, &mut *cancelled)
                .map_err(PlateRigidImpactError::MultipleContact)?;
        let mut receipts = vec![None; bodies.len()];
        if contact.normal_impulses.iter().all(|&j| j == 0.0) {
            poll(cancelled)?;
            return Ok(PlateRigidImpactNetworkReport {
                relative_normal_speeds_after: contacts
                    .iter()
                    .map(|c| c.relative_normal_speed)
                    .collect(),
                contact,
                delassus_kg_inverse: delassus,
                plate: None,
                bodies: receipts,
                balance: None,
            });
        }

        let map = &first.dynamics.model().dof_map;
        let mut full_impulse = vec![0.0; map.len()];
        let mut impulse = vec![0.0; first.unit_impulse.len()];
        for (r, value) in full_impulse.iter_mut().enumerate() {
            if r % 256 == 0 {
                poll(cancelled)?;
            }
            let mut sum = Dd::ZERO;
            for (i, &j) in contact.normal_impulses.iter().enumerate() {
                sum = sum - dd(j) * dd(self.prepared.row(i).2[r]);
            }
            *value = sum.to_f64();
            if let Some(free) = map[r] {
                impulse[free] = *value;
            }
        }
        let mut resultants = vec![
            BodyResultant {
                linear: Dd::ZERO,
                angular: [Dd::ZERO; 3],
            };
            bodies.len()
        ];
        for (i, &j) in contact.normal_impulses.iter().enumerate() {
            let resultant = &mut resultants[witnesses[i].body];
            resultant.linear = resultant.linear + dd(j);
            for (sum, torque) in resultant.angular.iter_mut().zip(components(torques[i])) {
                *sum = *sum + dd(j) * dd(torque);
            }
        }
        let (trial, plate_receipt) = first.stage_plate_impulse(plate, &impulse, cancelled)?;
        for (index, resultant) in resultants.iter().enumerate() {
            poll(cancelled)?;
            if resultant.linear == Dd::ZERO && resultant.angular.iter().all(|x| *x == Dd::ZERO) {
                continue;
            }
            let participant = bodies[index];
            receipts[index] = Some(
                participant
                    .state
                    .apply_momentum_impulse(
                        participant.properties,
                        Vec3::new(0.0, 0.0, resultant.linear.to_f64()),
                        vector(resultant.angular.map(Dd::to_f64)),
                    )
                    .map_err(PlateRigidImpactError::Body)?,
            );
        }
        let mut after = Vec::with_capacity(k);
        let mut relative_work = Dd::ZERO;
        for (i, witness) in witnesses.iter().enumerate() {
            poll(cancelled)?;
            let participant = bodies[witness.body];
            let state = receipts[witness.body]
                .as_ref()
                .map_or(participant.state, |r| r.state_after);
            let motion = state
                .point_kinematics(participant.properties, witness.arm_body)
                .map_err(PlateRigidImpactError::Body)?;
            let speed = relative_speed(self.prepared.row(i).1, &motion, &trial.v);
            check(
                "coupled endpoint speed",
                speed - contact.relative_normal_speeds_after[i],
                first.controls.speed_tolerance_m_s
                    + first.controls.impulse.relative_impulse_tolerance
                        * contacts[i].relative_normal_speed.abs().max(speed.abs()),
            )?;
            relative_work = relative_work
                + dd(contact.normal_impulses[i])
                    * (dd(0.5 * contacts[i].relative_normal_speed) + dd(0.5 * speed));
            after.push(speed);
        }
        check(
            "original contact work",
            (relative_work - dd(contact.relative_work)).to_f64(),
            first.controls.energy_tolerance_j,
        )?;
        let transfers = bodies.iter().zip(&receipts).zip(&resultants).filter_map(
            |((body, receipt), resultant)| {
                receipt.as_ref().map(|r| PlateBodyTransfer {
                    before: body.state,
                    after: r.state_after,
                    properties: body.properties,
                    linear_impulse: [Dd::ZERO, Dd::ZERO, resultant.linear],
                    angular_impulse: resultant.angular,
                })
            },
        );
        let balance = first.account_transfers(
            plate,
            &trial,
            &plate_receipt,
            transfers,
            |i| full_impulse[i],
            relative_work,
            contact.dissipated_energy,
            cancelled,
        )?;
        poll(cancelled)?;
        plate.v = trial.v;
        plate.a = trial.a;
        for (body, receipt) in bodies.iter_mut().zip(&receipts) {
            if let Some(receipt) = receipt {
                body.state = receipt.state_after;
            }
        }
        Ok(PlateRigidImpactNetworkReport {
            contact,
            delassus_kg_inverse: delassus,
            relative_normal_speeds_after: after,
            plate: Some(plate_receipt),
            bodies: receipts,
            balance: Some(balance),
        })
    }
}
