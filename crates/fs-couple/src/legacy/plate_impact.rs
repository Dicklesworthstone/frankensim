//! One reciprocal normal impulse between an original spatial plate and a
//! finite rigid body. The point load/observation, complete plate mass, tensor
//! rigid inertia and restitution law retain their existing numerical owners.
//!
//! This is a localized event, at fixed configuration, in the flat plate's +z
//! direction. The caller supplies a body material point above the selected
//! plate point. [`multiple`] couples sites on one body; [`bodies`] couples
//! several finite bodies through the same original mass and momentum owners.
//! [`plates`] exchanges normal momentum between two original spatial plates.
//! These paths do not detect collisions, rotate the plate normal, apply
//! friction, or deposit impact loss as heat.
//! Momentum resultants are those represented by this transverse plate model:
//! Pz and the two bending angular momenta Lx/Ly. In-plane translation and
//! drilling rotation require the corresponding shell/solid dynamics owner.

use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_math::dd::Dd;
use fs_mbd::{
    DynamicsError, ImpulseReceipt, MassProperties, PointKinematics, RigidBodyState, Vec3,
};
use fs_plate::PlateMesh;
use fs_plate::loading::{PlateLoadBudget, PlateLoadError, PlatePointStencil};
use fs_plate_transient::{PlateDynamics, PlateLoad};
use fs_sparse::Csr;
use fs_time::galpha::SecondOrderState;
use fs_time::galpha::impulse::{
    ImpulseControls, ImpulseError, SecondOrderImpulse, impulse_workspace_components,
    second_order_impulse,
};
use fs_time::galpha::initialization::{InitialSolveError, second_order_acceleration};
use fs_time::galpha::second_order_adjoint::trajectory::StructuralTrajectoryModel;

#[path = "plate_impact/multiple.rs"]
pub mod multiple;

#[path = "plate_impact/bodies.rs"]
pub mod bodies;

#[path = "plate_impact/plates.rs"]
pub mod plates;

/// Simultaneous normal contact across independent native plate trajectories.
#[path = "plate_impact/network.rs"]
pub mod network;

/// Reciprocal contact with original rotating rod attachments and full inertia.
#[cfg(feature = "rod-dynamics")]
#[path = "plate_impact/rods.rs"]
pub mod rods;

/// Physical and numerical allowances. Background momentum never enlarges an
/// impulse allowance. All supports are stationary in the original plate frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateRigidImpactControls {
    pub impulse: ImpulseControls,
    pub restitution: f64,
    pub activation_tolerance_m: f64,
    /// Allowed xy mismatch of the original body and plate material points.
    pub point_tolerance_m: f64,
    pub speed_tolerance_m_s: f64,
    pub linear_impulse_tolerance_ns: f64,
    pub angular_impulse_tolerance_nms: f64,
    /// Combined kinetic change plus the original scalar-law dissipation.
    pub energy_tolerance_j: f64,
}

/// Pre-scan limits, including the original uneliminated mass needed to recover
/// actual dynamic support impulses. No dense matrix is formed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateRigidImpactBudget {
    pub location: PlateLoadBudget,
    pub max_full_mass_nonzeros: usize,
    /// Includes native solver scratch and this adapter's live scalar arrays;
    /// excludes borrowed inputs and pre-existing trajectory telemetry.
    pub max_workspace_components: usize,
}

#[derive(Debug)]
pub enum PlateRigidImpactError {
    Invalid(&'static str),
    Budget(&'static str),
    Point(PlateLoadError),
    Mass(InitialSolveError),
    Plate(ImpulseError),
    Body(DynamicsError),
    Contact(NormalImpactError),
    MultipleContact(fs_contact::multicontact::MultiImpactError),
    Load(String),
    Balance {
        what: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for PlateRigidImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "plate/rigid impact refused: {self:?}")
    }
}
impl core::error::Error for PlateRigidImpactError {}

/// Accepted actual changes; both original states continue from the same time
/// and configuration. Positive support impulse acts **on the plate**.
#[derive(Debug)]
pub struct PlateRigidImpactReport {
    pub contact: NormalImpactResult,
    pub plate: Option<SecondOrderImpulse>,
    pub body: Option<ImpulseReceipt>,
    /// Full `(w, wx, wy)` support rows, including constrained-row cross inertia.
    /// These are complete impulse reactions, not just applied fixed-DOF loads.
    pub support_impulses: Vec<(usize, f64)>,
    /// Represented transverse momentum `(0, 0, delta-Pz)`.
    pub plate_linear_impulse_world: Vec3,
    /// Represented bending angular momentum `(delta-Lx, delta-Ly, 0)`.
    pub plate_angular_impulse_world: Vec3,
    /// Transverse support resultant in the same flat-plate model.
    pub support_linear_impulse_world: Vec3,
    /// Bending support moment about the original world origin.
    pub support_angular_impulse_world: Vec3,
    /// Plate plus body, less the external support reaction.
    pub linear_balance_residual: Vec3,
    pub angular_balance_residual: Vec3,
    pub relative_normal_speed_after: f64,
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
}

/// Actual native transfer shared by isolated and simultaneous impacts.
/// For redundant contacts, support entries correspond to the contact owner's
/// checked impulse representative; they do not establish a unique local load.
#[derive(Debug)]
pub struct PlateImpulseBalance {
    pub support_impulses: Vec<(usize, f64)>,
    pub plate_linear_impulse_world: Vec3,
    pub plate_angular_impulse_world: Vec3,
    pub support_linear_impulse_world: Vec3,
    pub support_angular_impulse_world: Vec3,
    pub linear_balance_residual: Vec3,
    pub angular_balance_residual: Vec3,
    pub kinetic_energy_change_j: f64,
    pub work_energy_residual_j: f64,
    pub energy_residual_j: f64,
}

/// One actual native rigid update in the common plate momentum/work audit.
/// The iterator supplying these records does not copy an entire body history.
struct PlateBodyTransfer {
    before: RigidBodyState,
    after: RigidBodyState,
    properties: MassProperties,
    // Retain the unrounded aggregate request through the actual-state audit.
    linear_impulse: [Dd; 3],
    angular_impulse: [Dd; 3],
}

// Keep each original plate's contribution unrounded until every participant
// enters the common event balance, including constrained-row support inertia.
struct PlateTransferAudit {
    support_impulses: Vec<(usize, f64)>,
    linear: [Dd; 3],
    angular: [Dd; 3],
    support_linear: [Dd; 3],
    support_angular: [Dd; 3],
    kinetic_energy: Dd,
}

/// Prepared point response of the same spatial trajectory. Construction solves
/// the original scaled mass once for the point mobility. Each impacted event
/// uses the native impulse update and native acceleration restart; it never
/// substitutes lumped mass, a modal truncation or a fresh trajectory.
pub struct PlateRigidImpact<'a, L: ?Sized> {
    dynamics: PlateDynamics<'a, L>,
    mesh: &'a PlateMesh,
    full_mass: &'a Csr,
    point: PlatePointStencil,
    unit_impulse: Vec<f64>,
    full_unit_impulse: Vec<f64>,
    inverse_plate_mass: f64,
    controls: PlateRigidImpactControls,
    native_workspace: usize,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateRigidImpactError> {
    if cancelled() {
        Err(PlateRigidImpactError::Cancelled)
    } else {
        Ok(())
    }
}
fn check(what: &'static str, residual: f64, allowance: f64) -> Result<(), PlateRigidImpactError> {
    if !residual.is_finite() || !allowance.is_finite() || residual.abs() > allowance {
        Err(PlateRigidImpactError::Balance {
            what,
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}
fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .fold(Dd::ZERO, |s, (&a, &b)| s + dd(a) * dd(b))
        .to_f64()
}

impl<'a, L: PlateLoad + ?Sized> PlateRigidImpact<'a, L> {
    /// `full_mass` is the ORIGINAL uneliminated reference mass in mesh `(w,wx,wy)`
    /// order. Its supported principal submatrix must exactly equal the mass
    /// borrowed by `dynamics`; the dynamics' original mass scaling is reused.
    /// For a free plate this is simply `&dynamics.model().m`. For a supported
    /// plate retain the same assembly before support elimination. Missing
    /// constrained rows cannot be reconstructed from a reduced pencil.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        dynamics: PlateDynamics<'a, L>,
        mesh: &'a PlateMesh,
        full_mass: &'a Csr,
        point_xy: [f64; 2],
        controls: PlateRigidImpactControls,
        budget: PlateRigidImpactBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateRigidImpactError> {
        Self::prepare_with_response(
            dynamics, mesh, full_mass, point_xy, controls, budget, cancelled,
        )
        .map(|(prepared, _)| prepared)
    }

    // The simultaneous adapter reuses the first solved response; preparing a
    // set must not solve its first mass column twice or repeat mass admission.
    #[allow(clippy::too_many_arguments)]
    fn prepare_with_response(
        dynamics: PlateDynamics<'a, L>,
        mesh: &'a PlateMesh,
        full_mass: &'a Csr,
        point_xy: [f64; 2],
        controls: PlateRigidImpactControls,
        budget: PlateRigidImpactBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, Vec<f64>), PlateRigidImpactError> {
        poll(cancelled)?;
        let model = dynamics.model();
        let n = model.free;
        let full = model.dof_map.len();
        let config = controls.impulse.mass_solve;
        let native_workspace = impulse_workspace_components(n, config).ok_or(
            PlateRigidImpactError::Budget("workspace dimension overflow"),
        )?;
        let required = n
            .checked_mul(12)
            .and_then(|x| full.checked_mul(6).and_then(|y| x.checked_add(y)))
            .and_then(|x| x.checked_add(native_workspace))
            .ok_or(PlateRigidImpactError::Budget(
                "workspace dimension overflow",
            ))?;
        if required > budget.max_workspace_components
            || full_mass.nnz() > budget.max_full_mass_nonzeros
        {
            return Err(PlateRigidImpactError::Budget(
                "mass entries or live workspace",
            ));
        }
        let allowances = [
            controls.activation_tolerance_m,
            controls.point_tolerance_m,
            controls.speed_tolerance_m_s,
            controls.linear_impulse_tolerance_ns,
            controls.angular_impulse_tolerance_nms,
            controls.energy_tolerance_j,
            controls.impulse.absolute_impulse_tolerance,
            controls.impulse.work_tolerance,
        ];
        if allowances.iter().any(|v| !v.is_finite() || *v < 0.0)
            || !controls.restitution.is_finite()
            || !(0.0..=1.0).contains(&controls.restitution)
            || !controls.impulse.relative_impulse_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.impulse.relative_impulse_tolerance)
            || full_mass.nrows() != full
            || full_mass.ncols() != full
        {
            return Err(PlateRigidImpactError::Invalid(
                "controls or original full mass dimensions",
            ));
        }
        let point = PlatePointStencil::locate(mesh, model, point_xy, budget.location)
            .map_err(PlateRigidImpactError::Point)?;
        // Compare every retained nonzero, not just a diagonal or a random probe.
        // Explicit zeros are harmless. This scan also admits all support rows.
        for i in 0..full {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            let (cols, vals) = full_mass.row(i);
            if vals.iter().any(|v| !v.is_finite()) {
                return Err(PlateRigidImpactError::Invalid("nonfinite full mass"));
            }
            if let Some(r) = model.dof_map[i] {
                let (reduced_cols, reduced_vals) = model.m.row(r);
                let mut count = 0;
                for (&j, &value) in cols.iter().zip(vals) {
                    if let Some(c) = model.dof_map[j] {
                        let actual = reduced_cols
                            .binary_search(&c)
                            .map_or(0.0, |k| reduced_vals[k]);
                        if actual != value {
                            return Err(PlateRigidImpactError::Invalid(
                                "full/reduced mass mismatch",
                            ));
                        }
                        count += usize::from(value != 0.0);
                    }
                }
                if count != reduced_vals.iter().filter(|&&v| v != 0.0).count() {
                    return Err(PlateRigidImpactError::Invalid(
                        "full mass omits retained entries",
                    ));
                }
            }
        }
        let mut unit_impulse = vec![0.0; n];
        let applied = point
            .add_load([1.0, 0.0, 0.0], &mut unit_impulse)
            .map_err(PlateRigidImpactError::Point)?;
        let mut full_unit_impulse = vec![0.0; full];
        for (i, &r) in model.dof_map.iter().enumerate() {
            if let Some(r) = r {
                full_unit_impulse[i] = unit_impulse[r];
            }
        }
        for (i, value) in applied.fixed_dof_loads {
            full_unit_impulse[i] = value;
        }
        let zero = vec![0.0; n];
        let response = second_order_acceleration(
            &dynamics,
            &zero,
            &zero,
            &unit_impulse,
            &dynamics.mass_preconditioner(),
            config,
            native_workspace,
            cancelled,
        )
        .map_err(PlateRigidImpactError::Mass)?;
        let inverse_plate_mass = dot(&unit_impulse, &response.value);
        if !inverse_plate_mass.is_finite()
            || inverse_plate_mass < 0.0
            || (inverse_plate_mass == 0.0 && unit_impulse.iter().any(|&x| x != 0.0))
        {
            return Err(PlateRigidImpactError::Invalid(
                "nonpositive plate point mobility",
            ));
        }
        poll(cancelled)?;
        Ok((
            Self {
                dynamics,
                mesh,
                full_mass,
                point,
                unit_impulse,
                full_unit_impulse,
                inverse_plate_mass,
                controls,
                native_workspace,
            },
            response.value,
        ))
    }

    #[must_use]
    pub fn inverse_plate_effective_mass(&self) -> f64 {
        self.inverse_plate_mass
    }

    /// Apply a localized event to both original histories atomically. Body
    /// impulses point along +z; plate impulses are their exact transpose with
    /// the opposite sign. Native plate acceleration is recomputed using the
    /// original finite load, elasticity and damping at the unchanged time.
    /// Refusal, including a final cancellation, leaves both states unchanged.
    pub fn apply(
        &self,
        plate: &mut SecondOrderState,
        body: &mut RigidBodyState,
        properties: MassProperties,
        arm_body: Vec3,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateRigidImpactReport, PlateRigidImpactError> {
        poll(cancelled)?;
        self.admit_state(plate)?;
        let sample = self
            .point
            .sample(&plate.q, &plate.v)
            .map_err(PlateRigidImpactError::Point)?;
        let motion = body
            .point_kinematics(properties, arm_body)
            .map_err(PlateRigidImpactError::Body)?;
        let [x, y] = self.point.point();
        check(
            "contact point x",
            motion.point_world.x - x,
            self.controls.point_tolerance_m,
        )?;
        check(
            "contact point y",
            motion.point_world.y - y,
            self.controls.point_tolerance_m,
        )?;
        let relative_before = self.relative_speed(&motion, &plate.v);
        let inverse_body_mass = body
            .directional_effective_mass(properties, arm_body, Vec3::new(0.0, 0.0, 1.0))
            .map_err(PlateRigidImpactError::Body)?
            .inverse_mass_kg_inverse;
        let contact = solve_normal_impact(NormalImpactInput {
            gap: motion.point_world.z - sample.displacement,
            activation_tolerance: self.controls.activation_tolerance_m,
            relative_normal_speed: relative_before,
            inverse_effective_mass: self.inverse_plate_mass + inverse_body_mass,
            restitution: self.controls.restitution,
        })
        .map_err(PlateRigidImpactError::Contact)?;
        if !contact.impacted {
            poll(cancelled)?;
            return Ok(PlateRigidImpactReport {
                contact,
                plate: None,
                body: None,
                support_impulses: Vec::new(),
                plate_linear_impulse_world: Vec3::ZERO,
                plate_angular_impulse_world: Vec3::ZERO,
                support_linear_impulse_world: Vec3::ZERO,
                support_angular_impulse_world: Vec3::ZERO,
                linear_balance_residual: Vec3::ZERO,
                angular_balance_residual: Vec3::ZERO,
                relative_normal_speed_after: relative_before,
                kinetic_energy_change_j: 0.0,
                energy_residual_j: 0.0,
            });
        }
        let j = contact.normal_impulse;
        let impulse: Vec<_> = self.unit_impulse.iter().map(|&u| -j * u).collect();
        let (trial, plate_receipt) = self.stage_plate_impulse(plate, &impulse, cancelled)?;
        let rigid = body
            .apply_impulse_at_body_point(properties, arm_body, Vec3::new(0.0, 0.0, j))
            .map_err(PlateRigidImpactError::Body)?;
        let relative_after = self.relative_speed(&rigid.point_kinematics_after, &trial.v);
        let rel = self.controls.impulse.relative_impulse_tolerance;
        check(
            "restitution speed",
            relative_after - contact.relative_normal_speed_after,
            self.controls.speed_tolerance_m_s
                + rel * relative_before.abs().max(relative_after.abs()),
        )?;
        let balance = self.account_transfer(
            plate,
            &trial,
            &plate_receipt,
            *body,
            rigid.state_after,
            properties,
            rigid.impulse_world,
            rigid.angular_impulse_body,
            |i| -j * self.full_unit_impulse[i],
            -contact.dissipated_energy,
            contact.dissipated_energy,
            cancelled,
        )?;
        let report = PlateRigidImpactReport {
            contact,
            plate: Some(plate_receipt),
            body: Some(rigid),
            support_impulses: balance.support_impulses,
            plate_linear_impulse_world: balance.plate_linear_impulse_world,
            plate_angular_impulse_world: balance.plate_angular_impulse_world,
            support_linear_impulse_world: balance.support_linear_impulse_world,
            support_angular_impulse_world: balance.support_angular_impulse_world,
            linear_balance_residual: balance.linear_balance_residual,
            angular_balance_residual: balance.angular_balance_residual,
            relative_normal_speed_after: relative_after,
            kinetic_energy_change_j: balance.kinetic_energy_change_j,
            energy_residual_j: balance.energy_residual_j,
        };
        poll(cancelled)?;
        plate.v = trial.v;
        plate.a = trial.a;
        *body = rigid.state_after;
        Ok(report)
    }

    fn admit_state(&self, plate: &SecondOrderState) -> Result<(), PlateRigidImpactError> {
        let n = self.unit_impulse.len();
        if !plate.t.is_finite()
            || [&plate.q, &plate.v, &plate.a]
                .iter()
                .any(|v| v.len() != n || v.iter().any(|x| !x.is_finite()))
        {
            return Err(PlateRigidImpactError::Invalid(
                "finite matching plate checkpoint required",
            ));
        }
        Ok(())
    }

    fn stage_plate_impulse(
        &self,
        plate: &SecondOrderState,
        impulse: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(SecondOrderState, SecondOrderImpulse), PlateRigidImpactError> {
        let mut force = vec![f64::NAN; self.unit_impulse.len()];
        self.dynamics
            .forcing(plate.t, &mut force)
            .map_err(PlateRigidImpactError::Load)?;
        poll(cancelled)?;
        // Keep arbitrarily long telemetry in place. Only native v/a publish.
        let mut trial = SecondOrderState::new(plate.t, &plate.q, &plate.v, &plate.a);
        let receipt = second_order_impulse(
            &mut trial,
            &self.dynamics,
            impulse,
            &force,
            &self.dynamics.mass_preconditioner(),
            self.controls.impulse,
            self.native_workspace,
            cancelled,
        )
        .map_err(PlateRigidImpactError::Plate)?;
        Ok((trial, receipt))
    }

    // One owner for represented increments, complete support reactions and
    // combined compensated kinetic work, independent of the contact law.
    #[allow(clippy::too_many_arguments)]
    fn account_transfer(
        &self,
        plate: &SecondOrderState,
        trial: &SecondOrderState,
        plate_receipt: &SecondOrderImpulse,
        body: RigidBodyState,
        next_body: RigidBodyState,
        properties: MassProperties,
        linear_impulse: Vec3,
        angular_impulse: Vec3,
        full_impulse: impl Fn(usize) -> f64,
        relative_work: f64,
        dissipated_energy: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateImpulseBalance, PlateRigidImpactError> {
        self.account_transfers(
            plate,
            trial,
            plate_receipt,
            core::iter::once(PlateBodyTransfer {
                before: body,
                after: next_body,
                properties,
                linear_impulse: components(linear_impulse).map(dd),
                angular_impulse: components(angular_impulse).map(dd),
            }),
            full_impulse,
            dd(relative_work),
            dissipated_energy,
            cancelled,
        )
    }

    // Account every actual finite participant together with the ORIGINAL plate
    // and complete support inertia. Rounded per-body energies are not summed.
    #[allow(clippy::too_many_arguments)]
    fn account_transfers(
        &self,
        plate: &SecondOrderState,
        trial: &SecondOrderState,
        plate_receipt: &SecondOrderImpulse,
        bodies: impl IntoIterator<Item = PlateBodyTransfer>,
        full_impulse: impl Fn(usize) -> f64,
        relative_work: Dd,
        dissipated_energy: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateImpulseBalance, PlateRigidImpactError> {
        let rel = self.controls.impulse.relative_impulse_tolerance;
        let plate_transfer = self.account_plate_transfer(
            plate,
            trial,
            plate_receipt,
            |i| dd(full_impulse(i)),
            cancelled,
        )?;
        let mut linear_residual: [Dd; 3] =
            std::array::from_fn(|c| plate_transfer.linear[c] - plate_transfer.support_linear[c]);
        let mut angular_residual: [Dd; 3] =
            std::array::from_fn(|c| plate_transfer.angular[c] - plate_transfer.support_angular[c]);
        let mut linear_scale = Dd::ZERO;
        let mut energy = plate_transfer.kinetic_energy;
        for transfer in bodies {
            poll(cancelled)?;
            let PlateBodyTransfer {
                before: body,
                after: next_body,
                properties,
                linear_impulse,
                angular_impulse,
            } = transfer;
            let before_p = components(body.linear_momentum_world());
            let after_p = components(next_body.linear_momentum_world());
            let before_l = components(body.angular_momentum_body());
            let after_l = components(next_body.angular_momentum_body());
            let dp: [Dd; 3] = std::array::from_fn(|c| dd(after_p[c]) - dd(before_p[c]));
            let dl: [Dd; 3] = std::array::from_fn(|c| dd(after_l[c]) - dd(before_l[c]));
            for c in 0..3 {
                let commanded = linear_impulse[c];
                check(
                    "rigid linear transfer",
                    (dp[c] - commanded).to_f64(),
                    self.controls.linear_impulse_tolerance_ns
                        + rel * dp[c].to_f64().abs().max(commanded.to_f64().abs()),
                )?;
                let commanded = angular_impulse[c];
                check(
                    "rigid angular transfer",
                    (dl[c] - commanded).to_f64(),
                    self.controls.angular_impulse_tolerance_nms
                        + rel * dl[c].to_f64().abs().max(commanded.to_f64().abs()),
                )?;
            }
            linear_scale = linear_scale
                + dd(linear_impulse
                    .into_iter()
                    .map(Dd::to_f64)
                    .map(f64::abs)
                    .fold(0.0, f64::max));
            let spin_world = rotate_body_increment(body, dl);
            let orbital = cross_increment(body.pose().position_world(), dp);
            let before_omega = components(
                properties
                    .angular_velocity_body_checked(body.angular_momentum_body())
                    .map_err(PlateRigidImpactError::Body)?,
            );
            let after_omega = components(
                properties
                    .angular_velocity_body_checked(next_body.angular_momentum_body())
                    .map_err(PlateRigidImpactError::Body)?,
            );
            for c in 0..3 {
                linear_residual[c] = linear_residual[c] + dp[c];
                angular_residual[c] = angular_residual[c] + spin_world[c] + orbital[c];
                let mean_velocity =
                    (dd(0.5 * before_p[c]) + dd(0.5 * after_p[c])) / dd(properties.mass());
                energy = energy
                    + dp[c] * mean_velocity
                    + dl[c] * (dd(0.5 * before_omega[c]) + dd(0.5 * after_omega[c]));
            }
        }
        let linear_residual = vector(linear_residual.map(Dd::to_f64));
        let angular_residual = vector(angular_residual.map(Dd::to_f64));
        for residual in components(linear_residual) {
            check(
                "system linear balance",
                residual,
                self.controls.linear_impulse_tolerance_ns + rel * linear_scale.to_f64(),
            )?;
        }
        for residual in components(angular_residual) {
            check(
                "system angular balance",
                residual,
                self.controls.angular_impulse_tolerance_nms,
            )?;
        }
        let energy_residual = (energy + dd(dissipated_energy)).to_f64();
        let work_residual = (energy - relative_work).to_f64();
        check(
            "impact work",
            work_residual,
            self.controls.energy_tolerance_j,
        )?;
        check(
            "impact energy",
            energy_residual,
            self.controls.energy_tolerance_j,
        )?;
        Ok(PlateImpulseBalance {
            support_impulses: plate_transfer.support_impulses,
            plate_linear_impulse_world: vector(plate_transfer.linear.map(Dd::to_f64)),
            plate_angular_impulse_world: vector(plate_transfer.angular.map(Dd::to_f64)),
            support_linear_impulse_world: vector(plate_transfer.support_linear.map(Dd::to_f64)),
            support_angular_impulse_world: vector(plate_transfer.support_angular.map(Dd::to_f64)),
            linear_balance_residual: linear_residual,
            angular_balance_residual: angular_residual,
            kinetic_energy_change_j: energy.to_f64(),
            work_energy_residual_j: work_residual,
            energy_residual_j: energy_residual,
        })
    }

    fn account_plate_transfer(
        &self,
        plate: &SecondOrderState,
        trial: &SecondOrderState,
        receipt: &SecondOrderImpulse,
        full_impulse: impl Fn(usize) -> Dd,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateTransferAudit, PlateRigidImpactError> {
        let n = self.unit_impulse.len();
        let map = &self.dynamics.model().dof_map;
        let mut full_change = vec![0.0; map.len()];
        for (i, &r) in map.iter().enumerate() {
            if let Some(r) = r {
                full_change[i] = receipt.velocity_change[r];
            }
        }
        let mut momentum = vec![0.0; map.len()];
        self.full_mass.spmv(&full_change, &mut momentum);
        let scale = self.dynamics.parameters().mass_scale;
        let mut support_impulses = Vec::with_capacity(map.len() - n);
        for (i, value) in momentum.iter_mut().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            *value *= scale;
            let requested = full_impulse(i);
            if map[i].is_none() {
                support_impulses.push((i, (dd(*value) - requested).to_f64()));
            } else {
                // Retain the unrounded physical transpose request when an
                // assembled contact produces several contributions to a row.
                check(
                    "plate generalized transfer",
                    (dd(*value) - requested).to_f64(),
                    self.controls.impulse.absolute_impulse_tolerance
                        + self.controls.impulse.relative_impulse_tolerance
                            * value.abs().max(requested.to_f64().abs()),
                )?;
            }
        }
        let (linear, angular) = self.resultants(|i| dd(momentum[i]));
        let (support_linear, support_angular) = self.resultants(|i| {
            if map[i].is_none() {
                dd(momentum[i]) - full_impulse(i)
            } else {
                Dd::ZERO
            }
        });
        let mut kinetic_energy = Dd::ZERO;
        for i in 0..n {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            kinetic_energy = kinetic_energy
                + (dd(0.5 * plate.v[i]) + dd(0.5 * trial.v[i])) * dd(receipt.momentum_change[i]);
        }
        Ok(PlateTransferAudit {
            support_impulses,
            linear,
            angular,
            support_linear,
            support_angular,
            kinetic_energy,
        })
    }

    fn relative_speed(&self, motion: &PointKinematics, velocity: &[f64]) -> f64 {
        relative_speed(&self.unit_impulse, motion, velocity)
    }

    fn resultants(&self, impulse: impl Fn(usize) -> Dd) -> ([Dd; 3], [Dd; 3]) {
        // Rigid transverse fields are w=vz+omega_x*y-omega_y*x,
        // wx=-omega_y and wy=omega_x. Their exact transpose gives Pz,Lx,Ly.
        let mut p = Dd::ZERO;
        let mut lx = Dd::ZERO;
        let mut ly = Dd::ZERO;
        for (node, &(x, y)) in self.mesh.nodes.iter().enumerate() {
            p = p + impulse(3 * node);
            lx = lx + dd(y) * impulse(3 * node) + impulse(3 * node + 2);
            ly = ly - dd(x) * impulse(3 * node) - impulse(3 * node + 1);
        }
        ([Dd::ZERO, Dd::ZERO, p], [lx, ly, Dd::ZERO])
    }
}

fn relative_speed(row: &[f64], motion: &PointKinematics, velocity: &[f64]) -> f64 {
    // Subtract common translation before rounding the point's rotational flow.
    let mut relative = dd(motion.center_of_mass_velocity_world.z)
        + dd(motion.angular_velocity_world.x) * dd(motion.arm_world.y)
        - dd(motion.angular_velocity_world.y) * dd(motion.arm_world.x);
    for (&row, &v) in row.iter().zip(velocity) {
        relative = relative - dd(row) * dd(v);
    }
    relative.to_f64()
}

fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

fn cross_increment(arm: Vec3, impulse: [Dd; 3]) -> [Dd; 3] {
    [
        dd(arm.y) * impulse[2] - dd(arm.z) * impulse[1],
        dd(arm.z) * impulse[0] - dd(arm.x) * impulse[2],
        dd(arm.x) * impulse[1] - dd(arm.y) * impulse[0],
    ]
}

fn rotate_body_increment(body: RigidBodyState, increment: [Dd; 3]) -> [Dd; 3] {
    let orientation = body.pose().orientation();
    let high = components(orientation.rotate_body_to_world(vector(increment.map(|x| x.hi))));
    if increment.iter().all(|x| x.lo == 0.0) {
        return high.map(dd);
    }
    let low = components(orientation.rotate_body_to_world(vector(increment.map(|x| x.lo))));
    std::array::from_fn(|c| dd(high[c]) + dd(low[c]))
}
