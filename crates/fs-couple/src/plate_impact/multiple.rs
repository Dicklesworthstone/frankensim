//! Simultaneous normal contacts on one original plate and one finite body.
//!
//! The cached plate block is B M^-1 B^T. Each event adds the same rigid body's
//! complete cross mobility and delegates the single coupled law to fs-contact.
//! All contact impulses enter one native structural update and one native
//! rigid resultant update; row order is never a sequence of physical impacts.

use super::{
    PlateImpulseBalance, PlateRigidImpact, PlateRigidImpactBudget, PlateRigidImpactControls,
    PlateRigidImpactError, check, components, dd, dot, poll, relative_speed,
};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactResult, NormalContactImpact, solve_normal_impacts,
};
use fs_math::dd::Dd;
use fs_mbd::{MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};
use fs_plate::PlateMesh;
use fs_plate::loading::PlatePointStencil;
use fs_plate_transient::{PlateDynamics, PlateLoad};
use fs_sparse::Csr;
use fs_time::galpha::SecondOrderState;
use fs_time::galpha::impulse::{ImpulseControls, SecondOrderImpulse, impulse_workspace_components};
use fs_time::galpha::initialization::second_order_acceleration;

/// An original plate material point and its declared normal law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateImpactSite {
    pub point_xy: [f64; 2],
    pub activation_tolerance_m: f64,
    pub restitution: f64,
}

/// Native transfer bounds and the original simultaneous contact solver controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateRigidImpactSetControls {
    pub impulse: ImpulseControls,
    pub contact: MultiImpactConfig,
    pub point_tolerance_m: f64,
    pub speed_tolerance_m_s: f64,
    pub linear_impulse_tolerance_ns: f64,
    pub angular_impulse_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Bounds before any point storage, matrix assembly or native solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateRigidImpactSetBudget {
    /// Workspace includes the cached point rows/responses and event arrays.
    /// The contact owner's separately bounded (at most 64-row) solver scratch
    /// and borrowed native model/telemetry are excluded.
    pub plate: PlateRigidImpactBudget,
    pub max_contacts: usize,
    /// Maximum scalar products in forming the full plate cross mobility: k*k*n.
    pub max_mobility_products: usize,
}

/// Both continuing states have been published when this receipt is returned.
#[derive(Debug)]
pub struct PlateRigidImpactSetReport {
    /// Retains UniqueImpulses versus UniqueVelocity. A redundant impulse vector
    /// is one representative, not a unique pressure or local heat partition.
    pub contact: MultiImpactResult,
    pub relative_normal_speeds_after: Vec<f64>,
    pub plate: Option<SecondOrderImpulse>,
    pub body: Option<MomentumImpulseReceipt>,
    pub balance: Option<PlateImpulseBalance>,
}

struct PointRow {
    point: PlatePointStencil,
    unit: Vec<f64>,
    full_unit: Vec<f64>,
}

/// Geometry and consistent plate cross mobility, reusable across impact events.
/// The original plate mass is solved once per point at preparation. An event
/// performs one structural impulse update/restart, regardless of contact count.
pub struct PlateRigidImpactSet<'a, L: ?Sized> {
    pub(super) first: PlateRigidImpact<'a, L>,
    additional: Vec<PointRow>,
    pub(super) sites: Vec<PlateImpactSite>,
    plate_delassus: Vec<f64>,
    pub(super) contact: MultiImpactConfig,
}

impl<'a, L: PlateLoad + ?Sized> PlateRigidImpactSet<'a, L> {
    /// All sites use the original flat plate's +z normal. Body material arms
    /// are supplied at each event, in matching site order. Full reference mass
    /// and stationary support admission are identical to the isolated owner.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        dynamics: PlateDynamics<'a, L>,
        mesh: &'a PlateMesh,
        full_mass: &'a Csr,
        sites: &[PlateImpactSite],
        controls: PlateRigidImpactSetControls,
        budget: PlateRigidImpactSetBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateRigidImpactError> {
        poll(cancelled)?;
        let k = sites.len();
        let n = dynamics.model().free;
        let full = dynamics.model().dof_map.len();
        if k == 0 || k > 64 || k > budget.max_contacts || k > controls.contact.solver.max_variables
        {
            return Err(PlateRigidImpactError::Budget("simultaneous contact count"));
        }
        let products = k.checked_mul(k).and_then(|x| x.checked_mul(n));
        if products.is_none_or(|x| x > budget.max_mobility_products) {
            return Err(PlateRigidImpactError::Budget(
                "plate cross-mobility products",
            ));
        }
        let native = impulse_workspace_components(n, controls.impulse.mass_solve).ok_or(
            PlateRigidImpactError::Budget("workspace dimension overflow"),
        )?;
        // Includes preparation's k native response columns and the maximum
        // simultaneous lifetime of cached rows, matrix and publication scratch.
        let workspace = (|| {
            native
                .checked_add(n.checked_mul(12)?)?
                .checked_add(full.checked_mul(7)?)?
                .checked_add((k - 1).checked_mul(n.checked_add(full)?)?)?
                .checked_add(k.checked_mul(n)?)?
                .checked_add(k.checked_mul(k)?.checked_mul(4)?)?
                .checked_add(k.checked_mul(40)?)
        })()
        .ok_or(PlateRigidImpactError::Budget(
            "workspace dimension overflow",
        ))?;
        if workspace > budget.plate.max_workspace_components {
            return Err(PlateRigidImpactError::Budget("simultaneous live workspace"));
        }
        for site in sites {
            poll(cancelled)?;
            if !site.activation_tolerance_m.is_finite()
                || site.activation_tolerance_m < 0.0
                || !site.restitution.is_finite()
                || !(0.0..=1.0).contains(&site.restitution)
            {
                return Err(PlateRigidImpactError::Invalid("normal site law"));
            }
        }
        let single = PlateRigidImpactControls {
            impulse: controls.impulse,
            restitution: sites[0].restitution,
            activation_tolerance_m: sites[0].activation_tolerance_m,
            point_tolerance_m: controls.point_tolerance_m,
            speed_tolerance_m_s: controls.speed_tolerance_m_s,
            linear_impulse_tolerance_ns: controls.linear_impulse_tolerance_ns,
            angular_impulse_tolerance_nms: controls.angular_impulse_tolerance_nms,
            energy_tolerance_j: controls.energy_tolerance_j,
        };
        let (first, first_response) = PlateRigidImpact::prepare_with_response(
            dynamics,
            mesh,
            full_mass,
            sites[0].point_xy,
            single,
            budget.plate,
            cancelled,
        )?;
        let mut additional = Vec::with_capacity(k - 1);
        let mut responses = Vec::with_capacity(k);
        responses.push(first_response);
        let zero = vec![0.0; n];
        for site in &sites[1..] {
            poll(cancelled)?;
            let point = PlatePointStencil::locate(
                mesh,
                first.dynamics.model(),
                site.point_xy,
                budget.plate.location,
            )
            .map_err(PlateRigidImpactError::Point)?;
            let mut unit = vec![0.0; n];
            let applied = point
                .add_load([1.0, 0.0, 0.0], &mut unit)
                .map_err(PlateRigidImpactError::Point)?;
            let mut full_unit = vec![0.0; full];
            for (i, r) in first.dynamics.model().dof_map.iter().enumerate() {
                if let Some(r) = r {
                    full_unit[i] = unit[*r];
                }
            }
            for (i, value) in applied.fixed_dof_loads {
                full_unit[i] = value;
            }
            let response = second_order_acceleration(
                &first.dynamics,
                &zero,
                &zero,
                &unit,
                &first.dynamics.mass_preconditioner(),
                controls.impulse.mass_solve,
                first.native_workspace,
                cancelled,
            )
            .map_err(PlateRigidImpactError::Mass)?;
            let mobility = dot(&unit, &response.value);
            if !mobility.is_finite()
                || mobility < 0.0
                || (mobility == 0.0 && unit.iter().any(|&x| x != 0.0))
            {
                return Err(PlateRigidImpactError::Invalid(
                    "nonpositive plate point mobility",
                ));
            }
            additional.push(PointRow {
                point,
                unit,
                full_unit,
            });
            responses.push(response.value);
        }
        let mut prepared = Self {
            first,
            additional,
            sites: sites.to_vec(),
            plate_delassus: vec![0.0; k * k],
            contact: controls.contact,
        };
        for i in 0..k {
            for (j, response) in responses.iter().enumerate() {
                let mut value = Dd::ZERO;
                let row = prepared.row(i).1;
                for offset in (0..n).step_by(256) {
                    poll(cancelled)?;
                    for r in offset..(offset + 256).min(n) {
                        value = value + dd(row[r]) * dd(response[r]);
                    }
                }
                let value = value.to_f64();
                if !value.is_finite() {
                    return Err(PlateRigidImpactError::Invalid(
                        "unrepresentable plate cross mobility",
                    ));
                }
                // Retain both independently solved cross responses. The native
                // contact owner admits symmetry; no averaging hides a bad solve.
                prepared.plate_delassus[i * k + j] = value;
            }
        }
        poll(cancelled)?;
        Ok(prepared)
    }

    /// Row-major original plate contribution B M^-1 B^T [1/kg].
    #[must_use]
    pub fn plate_cross_mobility(&self) -> &[f64] {
        &self.plate_delassus
    }

    /// Resolve the entire localized event, retaining original time/poses/history.
    /// Geometrically active separating sites participate in the coupled law:
    /// another contact can otherwise drive them into closing motion. Inactive
    /// rows still retain their predicted and measured cross response.
    pub fn apply(
        &self,
        plate: &mut SecondOrderState,
        body: &mut RigidBodyState,
        properties: MassProperties,
        arms_body: &[Vec3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateRigidImpactSetReport, PlateRigidImpactError> {
        poll(cancelled)?;
        self.first.admit_state(plate)?;
        let k = self.sites.len();
        if arms_body.len() != k {
            return Err(PlateRigidImpactError::Invalid(
                "one body witness per original site required",
            ));
        }
        let normal_body = body
            .pose()
            .orientation()
            .rotate_world_to_body(Vec3::new(0.0, 0.0, 1.0));
        let mut contacts = Vec::with_capacity(k);
        let mut torques = Vec::with_capacity(k);
        let mut rotations = Vec::with_capacity(k);
        for (i, (&arm, site)) in arms_body.iter().zip(&self.sites).enumerate() {
            poll(cancelled)?;
            let (point, row, _) = self.row(i);
            let sample = point
                .sample(&plate.q, &plate.v)
                .map_err(PlateRigidImpactError::Point)?;
            let motion = body
                .point_kinematics(properties, arm)
                .map_err(PlateRigidImpactError::Body)?;
            let [x, y] = point.point();
            check(
                "contact point x",
                motion.point_world.x - x,
                self.first.controls.point_tolerance_m,
            )?;
            check(
                "contact point y",
                motion.point_world.y - y,
                self.first.controls.point_tolerance_m,
            )?;
            contacts.push(NormalContactImpact {
                gap: motion.point_world.z - sample.displacement,
                activation_tolerance: site.activation_tolerance_m,
                relative_normal_speed: relative_speed(row, &motion, &plate.v),
                restitution: site.restitution,
            });
            let torque = arm.cross(normal_body);
            torques.push(torque);
            rotations.push(
                properties
                    .angular_velocity_body_checked(torque)
                    .map_err(PlateRigidImpactError::Body)?,
            );
        }
        let mut delassus = self.plate_delassus.clone();
        for i in 0..k {
            poll(cancelled)?;
            for j in 0..k {
                delassus[i * k + j] = (dd(delassus[i * k + j])
                    + dd(properties.mass().recip())
                    + dd(dot(&components(torques[i]), &components(rotations[j]))))
                .to_f64();
            }
        }
        let contact = solve_normal_impacts(&contacts, &delassus, self.contact, &mut *cancelled)
            .map_err(PlateRigidImpactError::MultipleContact)?;
        if contact.normal_impulses.iter().all(|&j| j == 0.0) {
            poll(cancelled)?;
            return Ok(PlateRigidImpactSetReport {
                relative_normal_speeds_after: contacts
                    .iter()
                    .map(|c| c.relative_normal_speed)
                    .collect(),
                contact,
                plate: None,
                body: None,
                balance: None,
            });
        }
        let map = &self.first.dynamics.model().dof_map;
        let mut full_impulse = vec![0.0; map.len()];
        let mut impulse = vec![0.0; self.first.unit_impulse.len()];
        for (r, value) in full_impulse.iter_mut().enumerate() {
            if r % 256 == 0 {
                poll(cancelled)?;
            }
            let mut sum = Dd::ZERO;
            for (i, &j) in contact.normal_impulses.iter().enumerate() {
                sum = sum - dd(j) * dd(self.row(i).2[r]);
            }
            *value = sum.to_f64();
            if let Some(free) = map[r] {
                impulse[free] = *value;
            }
        }
        let mut linear = Dd::ZERO;
        let mut angular = [Dd::ZERO; 3];
        for (&j, torque) in contact.normal_impulses.iter().zip(&torques) {
            linear = linear + dd(j);
            for (sum, value) in angular.iter_mut().zip(components(*torque)) {
                *sum = *sum + dd(j) * dd(value);
            }
        }
        let linear = Vec3::new(0.0, 0.0, linear.to_f64());
        let angular = Vec3::new(
            angular[0].to_f64(),
            angular[1].to_f64(),
            angular[2].to_f64(),
        );
        let (trial, plate_receipt) = self.first.stage_plate_impulse(plate, &impulse, cancelled)?;
        let rigid = body
            .apply_momentum_impulse(properties, linear, angular)
            .map_err(PlateRigidImpactError::Body)?;
        let mut after = Vec::with_capacity(k);
        for (i, &arm) in arms_body.iter().enumerate() {
            poll(cancelled)?;
            let motion = rigid
                .state_after
                .point_kinematics(properties, arm)
                .map_err(PlateRigidImpactError::Body)?;
            let speed = relative_speed(self.row(i).1, &motion, &trial.v);
            let tolerance = self.first.controls.speed_tolerance_m_s
                + self.first.controls.impulse.relative_impulse_tolerance
                    * contacts[i].relative_normal_speed.abs().max(speed.abs());
            check(
                "coupled endpoint speed",
                speed - contact.relative_normal_speeds_after[i],
                tolerance,
            )?;
            after.push(speed);
        }
        let balance = self.first.account_transfer(
            plate,
            &trial,
            &plate_receipt,
            *body,
            rigid.state_after,
            properties,
            linear,
            angular,
            |i| full_impulse[i],
            contact.relative_work,
            contact.dissipated_energy,
            cancelled,
        )?;
        let report = PlateRigidImpactSetReport {
            contact,
            relative_normal_speeds_after: after,
            plate: Some(plate_receipt),
            body: Some(rigid),
            balance: Some(balance),
        };
        poll(cancelled)?;
        plate.v = trial.v;
        plate.a = trial.a;
        *body = rigid.state_after;
        Ok(report)
    }

    pub(super) fn row(&self, i: usize) -> (&PlatePointStencil, &[f64], &[f64]) {
        if i == 0 {
            (
                &self.first.point,
                &self.first.unit_impulse,
                &self.first.full_unit_impulse,
            )
        } else {
            let row = &self.additional[i - 1];
            (&row.point, &row.unit, &row.full_unit)
        }
    }
}
