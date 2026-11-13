//! Simultaneous normal exchange between an original rod and spatial plate.
//!
//! The original rod's full spatial inertia and supported degrees of freedom
//! add to the plate's cached consistent-mass response. A single contact law
//! drives one native point-impulse batch and one structural impulse/restart.

use super::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use super::plates::{PlateImpactSurface, requests};
use super::{PlateRigidImpactError, PlateTransferAudit, check, dd, poll, vector};
use crate::rod_motion::impulse::{RodImpulseControls, RodImpulseReport, RodPointImpulse};
use crate::rod_motion::rigid_impact::{RodRigidImpactError, RodTransferAudit};
use crate::rod_motion::{RodMotion, RodMotionInertia};
use fs_contact::multicontact::{MultiImpactResult, NormalContactImpact, solve_normal_impacts};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::Vec3;
use fs_plate_transient::PlateLoad;
use fs_time::galpha::SecondOrderState;
use fs_time::galpha::impulse::SecondOrderImpulse;

/// Original rod material point corresponding to one prepared plate site.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodPlatePoint {
    pub node: usize,
    /// Offset from the material node, in its rotating material axes [m].
    pub point_body_m: [f64; 3],
}

/// Admission before preparing plate responses or allocating event scratch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RodPlateImpactBudget {
    pub plate: PlateRigidImpactSetBudget,
    pub maximum_rod_nodes: usize,
    /// Includes the plate workspace ceiling, rod candidate and native impulse
    /// arrays, compensated nodal requests and contact response arrays. Shared
    /// immutable rod sources, borrowed inputs and the original contact solver's
    /// separately bounded scratch are excluded.
    pub max_workspace_components: usize,
}

#[derive(Debug)]
pub enum RodPlateImpactError {
    Plate(PlateRigidImpactError),
    Rod(RodRigidImpactError),
}
impl From<PlateRigidImpactError> for RodPlateImpactError {
    fn from(error: PlateRigidImpactError) -> Self {
        Self::Plate(error)
    }
}
impl From<RodRigidImpactError> for RodPlateImpactError {
    fn from(error: RodRigidImpactError) -> Self {
        Self::Rod(error)
    }
}
impl core::fmt::Display for RodPlateImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rod/plate impact refused: {self:?}")
    }
}
impl core::error::Error for RodPlateImpactError {}

/// Actual native changes about world zero, ordered [rod, plate]. Rod support
/// details remain in its native receipt; plate entries are full (w, wx, wy)
/// support impulses including constrained-row cross inertia.
#[derive(Debug)]
pub struct RodPlateImpulseBalance {
    pub plate_support_impulses: Vec<(usize, f64)>,
    pub linear_impulses_world: [Vec3; 2],
    pub angular_impulses_world: [Vec3; 2],
    pub support_linear_impulses_world: [Vec3; 2],
    pub support_angular_impulses_world: [Vec3; 2],
    pub linear_balance_residual: Vec3,
    pub angular_balance_residual: Vec3,
    pub kinetic_energy_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub energy_residual_j: f64,
}

/// Both original trajectories have been published when this report returns.
#[derive(Debug)]
pub struct RodPlateImpactReport {
    /// Original global admissibility and unique-impulse/unique-velocity result.
    /// Redundant impulses are not a unique local pressure or heat partition.
    pub contact: MultiImpactResult,
    /// Complete row-major contact mobility at the actual rod configuration.
    pub delassus_kg_inverse: Vec<f64>,
    pub gaps_m: Vec<f64>,
    pub relative_normal_speeds_before_m_s: Vec<f64>,
    pub relative_normal_speeds_after_m_s: Vec<f64>,
    pub rod: Option<RodImpulseReport>,
    pub plate: Option<SecondOrderImpulse>,
    /// Absent only for an all-zero contact event.
    pub balance: Option<RodPlateImpulseBalance>,
}

/// Prepared original plate sites. The rod lies on the positive (+z) side and
/// receives +z impulse; the plate receives its reciprocal -z action. The
/// declared face height is added to the native transverse displacement.
pub struct RodPlateImpact<'a, L: ?Sized> {
    plate: PlateRigidImpactSet<'a, L>,
    reference_height_m: f64,
    maximum_rod_nodes: usize,
}

impl<'a, L: PlateLoad + ?Sized> RodPlateImpact<'a, L> {
    pub fn prepare(
        surface: PlateImpactSurface<'a, L>,
        sites: &[PlateImpactSite],
        controls: PlateRigidImpactSetControls,
        budget: RodPlateImpactBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RodPlateImpactError> {
        poll(cancelled)?;
        if !surface.reference_height_m.is_finite()
            || [
                controls.linear_impulse_tolerance_ns,
                controls.angular_impulse_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(PlateRigidImpactError::Invalid(
                "finite face height and positive native rod impulse tolerances",
            )
            .into());
        }
        let k = sites.len();
        if budget.maximum_rod_nodes == 0 || k == 0 || k > 64 {
            return Err(PlateRigidImpactError::Budget("rod/plate participant count").into());
        }
        let required = (|| {
            budget
                .plate
                .plate
                .max_workspace_components
                .checked_add(budget.maximum_rod_nodes.checked_mul(256)?)?
                .checked_add(surface.dynamics.model().dof_map.len().checked_mul(4)?)?
                .checked_add(k.checked_mul(k)?)?
                .checked_add(k.checked_mul(64)?)
        })()
        .ok_or(PlateRigidImpactError::Budget(
            "rod/plate workspace overflow",
        ))?;
        if required > budget.max_workspace_components {
            return Err(PlateRigidImpactError::Budget("rod/plate live workspace").into());
        }
        let plate = PlateRigidImpactSet::prepare(
            surface.dynamics,
            surface.mesh,
            surface.full_mass,
            sites,
            controls,
            budget.plate,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            plate,
            reference_height_m: surface.reference_height_m,
            maximum_rod_nodes: budget.maximum_rod_nodes,
        })
    }

    /// One frozen-configuration event at the same native physical time. The
    /// original point witnesses must share x/y within the declared geometry
    /// bound; the signed z gap still passes the original nonpenetration law.
    /// Every active row participates together, including separating rows.
    ///
    /// Point rotation uses full central tensors and eccentric mass offsets.
    /// Native support reactions, actual endpoint speeds, momentum and kinetic
    /// work all pass before either trajectory changes. Common translation is
    /// subtracted before adding point rotation; raw compensated contributions
    /// remain unrounded until both participants enter the event balance.
    ///
    /// Plate responses are cached; an event adds O(k²) rod cross responses and
    /// makes one O(rod nodes + k) original impulse batch plus one native plate
    /// mass solve/restart. No native time step, collision search, persistent
    /// constraint, tangential law or mechanical-loss-to-heat route is implied.
    pub fn apply<I: RodMotionInertia>(
        &self,
        rod: &mut RodMotion<I>,
        plate: &mut SecondOrderState,
        points: &[RodPlatePoint],
        gate: &CancelGate,
    ) -> Result<RodPlateImpactReport, RodPlateImpactError> {
        let cancelled = &mut || gate.is_requested();
        poll(cancelled)?;
        self.plate.first.admit_state(plate)?;
        let k = self.plate.sites.len();
        if points.len() != k || rod.time_s() != plate.t {
            return Err(PlateRigidImpactError::Invalid(
                "rod/plate point count and shared native physical time",
            )
            .into());
        }
        if rod.linear_momenta_kg_m_s().len() > self.maximum_rod_nodes {
            return Err(PlateRigidImpactError::Budget("rod/plate original rod nodes").into());
        }
        let controls = self.plate.first.controls;
        let mut contacts = Vec::with_capacity(k);
        let mut directions = Vec::with_capacity(k);
        let mut before_dd = Vec::with_capacity(k);
        for (i, point) in points.iter().enumerate() {
            poll(cancelled)?;
            let native = rod.contact_point_kinematics(point.node, point.point_body_m)?;
            let site = self.plate.sites[i];
            let mismatch = (native.point_world_dd_m[0] - dd(site.point_xy[0]))
                .to_f64()
                .hypot((native.point_world_dd_m[1] - dd(site.point_xy[1])).to_f64());
            check(
                "rod/plate material-point alignment",
                mismatch,
                controls.point_tolerance_m,
            )?;
            let row = self.plate.row(i).1;
            let speed = native.translation_m_s[2] - projection(row, &plate.v, cancelled)?
                + dd(native.rotation_m_s[2]);
            before_dd.push(speed);
            contacts.push(NormalContactImpact {
                gap: (native.point_world_dd_m[2]
                    - dd(self.reference_height_m)
                    - projection(row, &plate.q, cancelled)?)
                .to_f64(),
                activation_tolerance: site.activation_tolerance_m,
                relative_normal_speed: speed.to_f64(),
                restitution: site.restitution,
            });
            directions.push(rod.contact_direction_response(
                point.node,
                native.central_arm_body_m,
                [0.0, 0.0, 1.0],
            )?);
        }
        let mut delassus = self.plate.plate_cross_mobility().to_vec();
        for i in 0..k {
            poll(cancelled)?;
            for j in 0..=i {
                if points[i].node == points[j].node {
                    let mobility = (dd(delassus[i * k + j])
                        + rod.contact_cross_mobility(
                            points[i].node,
                            directions[i],
                            directions[j],
                        )?)
                    .to_f64();
                    if !mobility.is_finite() {
                        return Err(PlateRigidImpactError::Invalid(
                            "representable rod/plate contact mobility",
                        )
                        .into());
                    }
                    delassus[i * k + j] = mobility;
                    delassus[j * k + i] = mobility;
                }
            }
        }
        let contact =
            solve_normal_impacts(&contacts, &delassus, self.plate.contact, &mut *cancelled)
                .map_err(PlateRigidImpactError::MultipleContact)?;
        let gaps_m = contacts.iter().map(|c| c.gap).collect();
        let before: Vec<_> = before_dd.iter().map(|v| v.to_f64()).collect();
        if contact.normal_impulses.iter().all(|&j| j == 0.0) {
            poll(cancelled)?;
            return Ok(RodPlateImpactReport {
                contact,
                delassus_kg_inverse: delassus,
                gaps_m,
                relative_normal_speeds_after_m_s: before.clone(),
                relative_normal_speeds_before_m_s: before,
                rod: None,
                plate: None,
                balance: None,
            });
        }
        let (plate_request, plate_impulse) =
            requests(&self.plate, &contact.normal_impulses, -1.0, cancelled)?;
        let rod_impulses: Vec<_> = points
            .iter()
            .zip(&contact.normal_impulses)
            .filter(|(_, j)| **j != 0.0)
            .map(|(point, &j)| RodPointImpulse {
                node: point.node,
                point_body_m: point.point_body_m,
                linear_world_n_s: [0.0, 0.0, j],
                couple_world_nm_s: [0.0; 3],
            })
            .collect();
        let (rod_trial, rod_receipt) = rod
            .apply_point_impulses(
                &rod_impulses,
                RodImpulseControls {
                    balance_origin_m: [0.0; 3],
                    linear_tolerance_n_s: controls.linear_impulse_tolerance_ns,
                    angular_tolerance_nm_s: controls.angular_impulse_tolerance_nms,
                    work_tolerance_j: controls.energy_tolerance_j,
                },
                gate,
            )
            .map_err(RodRigidImpactError::from)?;
        let (plate_trial, plate_receipt) =
            self.plate
                .first
                .stage_plate_impulse(plate, &plate_impulse, cancelled)?;
        let mut after = Vec::with_capacity(k);
        let mut work = Dd::ZERO;
        for (i, point) in points.iter().enumerate() {
            poll(cancelled)?;
            let native = rod_trial.contact_point_kinematics(point.node, point.point_body_m)?;
            let speed = native.translation_m_s[2]
                - projection(self.plate.row(i).1, &plate_trial.v, cancelled)?
                + dd(native.rotation_m_s[2]);
            check(
                "rod/plate endpoint speed",
                (speed - dd(contact.relative_normal_speeds_after[i])).to_f64(),
                controls.speed_tolerance_m_s
                    + controls.impulse.relative_impulse_tolerance
                        * before[i].abs().max(speed.to_f64().abs()),
            )?;
            work = work + dd(contact.normal_impulses[i]) * dd(0.5) * (before_dd[i] + speed);
            after.push(speed.to_f64());
        }
        check(
            "original rod/plate contact work",
            (work - dd(contact.relative_work)).to_f64(),
            controls.energy_tolerance_j,
        )?;
        let rod_transfer =
            rod.account_contact_transfer(&rod_trial, Some(&rod_receipt), [0.0; 3], gate)?;
        let plate_transfer = self.plate.first.account_plate_transfer(
            plate,
            &plate_trial,
            &plate_receipt,
            |i| plate_request[i],
            cancelled,
        )?;
        let balance = self.account_pair(rod_transfer, plate_transfer, &contact, work)?;
        poll(cancelled)?;
        *rod = rod_trial;
        plate.v = plate_trial.v;
        plate.a = plate_trial.a;
        Ok(RodPlateImpactReport {
            contact,
            delassus_kg_inverse: delassus,
            gaps_m,
            relative_normal_speeds_before_m_s: before,
            relative_normal_speeds_after_m_s: after,
            rod: Some(rod_receipt),
            plate: Some(plate_receipt),
            balance: Some(balance),
        })
    }

    fn account_pair(
        &self,
        rod: RodTransferAudit,
        plate: PlateTransferAudit,
        contact: &MultiImpactResult,
        work: Dd,
    ) -> Result<RodPlateImpulseBalance, RodPlateImpactError> {
        let controls = self.plate.first.controls;
        let linear: [Dd; 3] = std::array::from_fn(|c| {
            rod.linear[c] - rod.support_linear[c] + plate.linear[c] - plate.support_linear[c]
        });
        let angular: [Dd; 3] = std::array::from_fn(|c| {
            rod.angular[c] - rod.support_angular[c] + plate.angular[c] - plate.support_angular[c]
        });
        let scale = contact
            .normal_impulses
            .iter()
            .fold(Dd::ZERO, |sum, &j| sum + dd(j.abs()));
        for c in 0..3 {
            check(
                "rod/plate linear balance",
                linear[c].to_f64(),
                controls.linear_impulse_tolerance_ns
                    + controls.impulse.relative_impulse_tolerance * scale.to_f64(),
            )?;
            check(
                "rod/plate angular balance",
                angular[c].to_f64(),
                controls.angular_impulse_tolerance_nms,
            )?;
        }
        let energy = rod.kinetic_energy + plate.kinetic_energy;
        let work_residual = (energy - work).to_f64();
        let energy_residual = (energy + dd(contact.dissipated_energy)).to_f64();
        check(
            "rod/plate impact work",
            work_residual,
            controls.energy_tolerance_j,
        )?;
        check(
            "rod/plate impact energy",
            energy_residual,
            controls.energy_tolerance_j,
        )?;
        Ok(RodPlateImpulseBalance {
            plate_support_impulses: plate.support_impulses,
            linear_impulses_world: [
                vector(rod.linear.map(Dd::to_f64)),
                vector(plate.linear.map(Dd::to_f64)),
            ],
            angular_impulses_world: [
                vector(rod.angular.map(Dd::to_f64)),
                vector(plate.angular.map(Dd::to_f64)),
            ],
            support_linear_impulses_world: [
                vector(rod.support_linear.map(Dd::to_f64)),
                vector(plate.support_linear.map(Dd::to_f64)),
            ],
            support_angular_impulses_world: [
                vector(rod.support_angular.map(Dd::to_f64)),
                vector(plate.support_angular.map(Dd::to_f64)),
            ],
            linear_balance_residual: vector(linear.map(Dd::to_f64)),
            angular_balance_residual: vector(angular.map(Dd::to_f64)),
            kinetic_energy_change_j: energy.to_f64(),
            relative_impulse_work_j: work.to_f64(),
            work_energy_residual_j: work_residual,
            energy_residual_j: energy_residual,
        })
    }
}

fn projection(
    row: &[f64],
    state: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, PlateRigidImpactError> {
    let mut value = Dd::ZERO;
    for offset in (0..row.len()).step_by(256) {
        poll(cancelled)?;
        for i in offset..(offset + 256).min(row.len()) {
            value = value + dd(row[i]) * dd(state[i]);
        }
    }
    Ok(value)
}
