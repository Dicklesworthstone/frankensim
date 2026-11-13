//! Reciprocal normal contact between two original spatial plate trajectories.
//!
//! Each plate retains its own mesh, consistent mass, supports and time owner.
//! The shared sites couple through the SUM of the two native point responses;
//! opposite transpose impulses enter one original update on each participant.

use super::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use super::{PlateRigidImpactError, PlateTransferAudit, check, dd, poll, vector};
use fs_contact::multicontact::{MultiImpactResult, NormalContactImpact, solve_normal_impacts};
use fs_math::dd::Dd;
use fs_mbd::Vec3;
use fs_plate::PlateMesh;
use fs_plate_transient::{PlateDynamics, PlateLoad};
use fs_sparse::Csr;
use fs_time::galpha::SecondOrderState;
use fs_time::galpha::impulse::SecondOrderImpulse;

/// One original flat contact face. The material displacement is added to its
/// declared world height; thickness or face offsets are never inferred.
/// Mesh coordinates are world x/y, and transverse motion is along world +z.
pub struct PlateImpactSurface<'a, L: ?Sized> {
    pub dynamics: PlateDynamics<'a, L>,
    pub mesh: &'a PlateMesh,
    pub full_mass: &'a Csr,
    pub reference_height_m: f64,
}

/// Both original preparation budgets and a bound on their combined live work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlatePlateImpactBudget {
    pub positive: PlateRigidImpactSetBudget,
    pub negative: PlateRigidImpactSetBudget,
    /// Includes both admitted plate workspace ceilings, compensated nodal
    /// requests and the shared contact arrays. Borrowed original trajectories
    /// and the contact solver's separately bounded scratch are excluded.
    pub max_workspace_components: usize,
}

/// Actual impulse contributions, ordered [positive plate, negative plate].
/// Only the transverse Pz and bending Lx/Ly represented by this model enter.
#[derive(Debug)]
pub struct PlatePlateImpulseBalance {
    pub support_impulses: [Vec<(usize, f64)>; 2],
    pub plate_linear_impulses_world: [Vec3; 2],
    pub plate_angular_impulses_world: [Vec3; 2],
    pub support_linear_impulses_world: [Vec3; 2],
    pub support_angular_impulses_world: [Vec3; 2],
    pub linear_balance_residual: Vec3,
    pub angular_balance_residual: Vec3,
    pub kinetic_energy_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub energy_residual_j: f64,
}

/// Returned only after both continuing trajectories have been published.
#[derive(Debug)]
pub struct PlatePlateImpactReport {
    /// Original global-law work and unique-impulse/unique-velocity distinction.
    /// Redundant impulses do not define unique local pressure or heat shares.
    pub contact: MultiImpactResult,
    pub gaps_m: Vec<f64>,
    pub relative_normal_speeds_before_m_s: Vec<f64>,
    pub relative_normal_speeds_after_m_s: Vec<f64>,
    /// Original native receipts in [positive, negative] order. Both are absent
    /// for an all-zero event; existing time/step/history state is retained.
    pub plates: [Option<SecondOrderImpulse>; 2],
    pub balance: Option<PlatePlateImpulseBalance>,
}

/// Prepared simultaneous contact on two independently meshed spatial plates.
/// Each shared world-x/y location is located on BOTH original meshes. Positive
/// impulse acts along +z on the positive plate and -z on the negative plate.
pub struct PlatePlateImpact<'a, P: ?Sized, N: ?Sized> {
    positive: PlateRigidImpactSet<'a, P>,
    negative: PlateRigidImpactSet<'a, N>,
    reference_gap: Dd,
    delassus: Vec<f64>,
}

impl<'a, P: PlateLoad + ?Sized, N: PlateLoad + ?Sized> PlatePlateImpact<'a, P, N> {
    /// Reuse the original stencils, complete support-inertia admission and
    /// consistent-mass response caches. The complete fixed-configuration
    /// contact mobility is W_positive + W_negative and is assembled once.
    /// No plate is replaced by a modal, lumped or rigid surrogate.
    pub fn prepare(
        positive: PlateImpactSurface<'a, P>,
        negative: PlateImpactSurface<'a, N>,
        sites: &[PlateImpactSite],
        controls: PlateRigidImpactSetControls,
        budget: PlatePlateImpactBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateRigidImpactError> {
        poll(cancelled)?;
        let k = sites.len();
        if k == 0 || k > 64 || k > controls.contact.solver.max_variables {
            return Err(PlateRigidImpactError::Budget("plate-pair contact count"));
        }
        if !positive.reference_height_m.is_finite() || !negative.reference_height_m.is_finite() {
            return Err(PlateRigidImpactError::Invalid("finite plate face heights"));
        }
        let reference_gap = dd(positive.reference_height_m) - dd(negative.reference_height_m);
        if !reference_gap.to_f64().is_finite() {
            return Err(PlateRigidImpactError::Invalid(
                "representable reference face separation",
            ));
        }
        let required = (|| {
            budget
                .positive
                .plate
                .max_workspace_components
                .checked_add(budget.negative.plate.max_workspace_components)?
                .checked_add(
                    positive
                        .dynamics
                        .model()
                        .dof_map
                        .len()
                        .checked_add(negative.dynamics.model().dof_map.len())?
                        .checked_mul(4)?,
                )?
                .checked_add(k.checked_mul(k)?)?
                .checked_add(k.checked_mul(24)?)
        })()
        .ok_or(PlateRigidImpactError::Budget(
            "plate-pair workspace overflow",
        ))?;
        if required > budget.max_workspace_components {
            return Err(PlateRigidImpactError::Budget("plate-pair live workspace"));
        }
        let positive = PlateRigidImpactSet::prepare(
            positive.dynamics,
            positive.mesh,
            positive.full_mass,
            sites,
            controls,
            budget.positive,
            cancelled,
        )?;
        let negative = PlateRigidImpactSet::prepare(
            negative.dynamics,
            negative.mesh,
            negative.full_mass,
            sites,
            controls,
            budget.negative,
            cancelled,
        )?;
        let mut delassus = Vec::with_capacity(k * k);
        for (&a, &b) in positive
            .plate_cross_mobility()
            .iter()
            .zip(negative.plate_cross_mobility())
        {
            let value = (dd(a) + dd(b)).to_f64();
            if !value.is_finite() {
                return Err(PlateRigidImpactError::Invalid(
                    "representable plate-pair mobility",
                ));
            }
            delassus.push(value);
        }
        poll(cancelled)?;
        Ok(Self {
            positive,
            negative,
            reference_gap,
            delassus,
        })
    }

    /// Cached row-major complete original mobility [kg^-1].
    #[must_use]
    pub fn contact_mobility(&self) -> &[f64] {
        &self.delassus
    }

    /// Resolve the original simultaneous normal law and atomically update both
    /// native velocities and accelerations at their shared physical time.
    ///
    /// Every active row participates, including initially separating contacts.
    /// Inactive rows still retain their actual cross response. Complete support
    /// reactions, each native generalized transfer, measured contact speeds and
    /// aggregate kinetic work must all pass before either participant changes.
    /// The original displacement, clocks, accepted steps and telemetry remain.
    ///
    /// Cached response assembly is quadratic in sites and linear in the two
    /// original state sizes. There is one impulse solve/restart per plate and
    /// no per-contact trajectory clone. This fixed-normal event does not locate
    /// collision times, evolve persistent contact, rotate normals, apply
    /// tangential impulses or assign mechanical loss to thermal destinations.
    pub fn apply(
        &self,
        positive: &mut SecondOrderState,
        negative: &mut SecondOrderState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlatePlateImpactReport, PlateRigidImpactError> {
        poll(cancelled)?;
        self.positive.first.admit_state(positive)?;
        self.negative.first.admit_state(negative)?;
        if positive.t != negative.t {
            return Err(PlateRigidImpactError::Invalid(
                "plate contact requires the same native physical time",
            ));
        }
        let k = self.positive.sites.len();
        let mut contacts = Vec::with_capacity(k);
        for (i, site) in self.positive.sites.iter().enumerate() {
            poll(cancelled)?;
            contacts.push(NormalContactImpact {
                gap: (self.reference_gap
                    + difference(
                        self.positive.row(i).1,
                        &positive.q,
                        self.negative.row(i).1,
                        &negative.q,
                        cancelled,
                    )?)
                .to_f64(),
                activation_tolerance: site.activation_tolerance_m,
                relative_normal_speed: difference(
                    self.positive.row(i).1,
                    &positive.v,
                    self.negative.row(i).1,
                    &negative.v,
                    cancelled,
                )?
                .to_f64(),
                restitution: site.restitution,
            });
        }
        let contact = solve_normal_impacts(
            &contacts,
            &self.delassus,
            self.positive.contact,
            &mut *cancelled,
        )
        .map_err(PlateRigidImpactError::MultipleContact)?;
        let gaps = contacts.iter().map(|c| c.gap).collect();
        let before: Vec<_> = contacts.iter().map(|c| c.relative_normal_speed).collect();
        if contact.normal_impulses.iter().all(|&j| j == 0.0) {
            poll(cancelled)?;
            return Ok(PlatePlateImpactReport {
                contact,
                gaps_m: gaps,
                relative_normal_speeds_after_m_s: before.clone(),
                relative_normal_speeds_before_m_s: before,
                plates: [None, None],
                balance: None,
            });
        }
        let (positive_request, positive_impulse) =
            requests(&self.positive, &contact.normal_impulses, 1.0, cancelled)?;
        let (negative_request, negative_impulse) =
            requests(&self.negative, &contact.normal_impulses, -1.0, cancelled)?;
        let (positive_trial, positive_receipt) =
            self.positive
                .first
                .stage_plate_impulse(positive, &positive_impulse, cancelled)?;
        let (negative_trial, negative_receipt) =
            self.negative
                .first
                .stage_plate_impulse(negative, &negative_impulse, cancelled)?;
        let controls = self.positive.first.controls;
        let mut after = Vec::with_capacity(k);
        let mut work = Dd::ZERO;
        for i in 0..k {
            poll(cancelled)?;
            let speed = difference(
                self.positive.row(i).1,
                &positive_trial.v,
                self.negative.row(i).1,
                &negative_trial.v,
                cancelled,
            )?
            .to_f64();
            check(
                "plate-pair endpoint speed",
                speed - contact.relative_normal_speeds_after[i],
                controls.speed_tolerance_m_s
                    + controls.impulse.relative_impulse_tolerance
                        * before[i].abs().max(speed.abs()),
            )?;
            work = work + dd(contact.normal_impulses[i]) * (dd(0.5 * before[i]) + dd(0.5 * speed));
            after.push(speed);
        }
        check(
            "original plate-pair contact work",
            (work - dd(contact.relative_work)).to_f64(),
            controls.energy_tolerance_j,
        )?;
        let positive_transfer = self.positive.first.account_plate_transfer(
            positive,
            &positive_trial,
            &positive_receipt,
            |i| positive_request[i],
            cancelled,
        )?;
        let negative_transfer = self.negative.first.account_plate_transfer(
            negative,
            &negative_trial,
            &negative_receipt,
            |i| negative_request[i],
            cancelled,
        )?;
        let balance = self.account_pair(positive_transfer, negative_transfer, &contact, work)?;
        poll(cancelled)?;
        positive.v = positive_trial.v;
        positive.a = positive_trial.a;
        negative.v = negative_trial.v;
        negative.a = negative_trial.a;
        Ok(PlatePlateImpactReport {
            contact,
            gaps_m: gaps,
            relative_normal_speeds_before_m_s: before,
            relative_normal_speeds_after_m_s: after,
            plates: [Some(positive_receipt), Some(negative_receipt)],
            balance: Some(balance),
        })
    }

    fn account_pair(
        &self,
        positive: PlateTransferAudit,
        negative: PlateTransferAudit,
        contact: &MultiImpactResult,
        work: Dd,
    ) -> Result<PlatePlateImpulseBalance, PlateRigidImpactError> {
        let controls = self.positive.first.controls;
        let linear: [Dd; 3] = std::array::from_fn(|c| {
            positive.linear[c] - positive.support_linear[c] + negative.linear[c]
                - negative.support_linear[c]
        });
        let angular: [Dd; 3] = std::array::from_fn(|c| {
            positive.angular[c] - positive.support_angular[c] + negative.angular[c]
                - negative.support_angular[c]
        });
        let impulse_scale = contact
            .normal_impulses
            .iter()
            .fold(Dd::ZERO, |sum, &j| sum + dd(j.abs()));
        for c in 0..3 {
            check(
                "plate-pair linear balance",
                linear[c].to_f64(),
                controls.linear_impulse_tolerance_ns
                    + controls.impulse.relative_impulse_tolerance * impulse_scale.to_f64(),
            )?;
            check(
                "plate-pair angular balance",
                angular[c].to_f64(),
                controls.angular_impulse_tolerance_nms,
            )?;
        }
        let energy = positive.kinetic_energy + negative.kinetic_energy;
        let work_residual = (energy - work).to_f64();
        let energy_residual = (energy + dd(contact.dissipated_energy)).to_f64();
        check(
            "plate-pair impact work",
            work_residual,
            controls.energy_tolerance_j,
        )?;
        check(
            "plate-pair impact energy",
            energy_residual,
            controls.energy_tolerance_j,
        )?;
        Ok(PlatePlateImpulseBalance {
            plate_linear_impulses_world: [
                vector(positive.linear.map(Dd::to_f64)),
                vector(negative.linear.map(Dd::to_f64)),
            ],
            plate_angular_impulses_world: [
                vector(positive.angular.map(Dd::to_f64)),
                vector(negative.angular.map(Dd::to_f64)),
            ],
            support_linear_impulses_world: [
                vector(positive.support_linear.map(Dd::to_f64)),
                vector(negative.support_linear.map(Dd::to_f64)),
            ],
            support_angular_impulses_world: [
                vector(positive.support_angular.map(Dd::to_f64)),
                vector(negative.support_angular.map(Dd::to_f64)),
            ],
            support_impulses: [positive.support_impulses, negative.support_impulses],
            linear_balance_residual: vector(linear.map(Dd::to_f64)),
            angular_balance_residual: vector(angular.map(Dd::to_f64)),
            kinetic_energy_change_j: energy.to_f64(),
            relative_impulse_work_j: work.to_f64(),
            work_energy_residual_j: work_residual,
            energy_residual_j: energy_residual,
        })
    }
}

pub(super) fn difference(
    positive: &[f64],
    pv: &[f64],
    negative: &[f64],
    nv: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, PlateRigidImpactError> {
    let mut value = Dd::ZERO;
    for offset in (0..positive.len()).step_by(256) {
        poll(cancelled)?;
        for i in offset..(offset + 256).min(positive.len()) {
            value = value + dd(positive[i]) * dd(pv[i]);
        }
    }
    for offset in (0..negative.len()).step_by(256) {
        poll(cancelled)?;
        for i in offset..(offset + 256).min(negative.len()) {
            value = value - dd(negative[i]) * dd(nv[i]);
        }
    }
    Ok(value)
}

pub(super) fn requests<L: PlateLoad + ?Sized>(
    prepared: &PlateRigidImpactSet<'_, L>,
    impulses: &[f64],
    sign: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<Dd>, Vec<f64>), PlateRigidImpactError> {
    let map = &prepared.first.dynamics.model().dof_map;
    let mut full = vec![Dd::ZERO; map.len()];
    let mut reduced = vec![0.0; prepared.first.unit_impulse.len()];
    for (r, value) in full.iter_mut().enumerate() {
        if r % 256 == 0 {
            poll(cancelled)?;
        }
        for (i, &j) in impulses.iter().enumerate() {
            *value = *value + dd(sign * j) * dd(prepared.row(i).2[r]);
        }
        if let Some(free) = map[r] {
            reduced[free] = value.to_f64();
        }
    }
    Ok((full, reduced))
}
