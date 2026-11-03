//! One reciprocal normal impulse between an original spatial plate and a
//! finite rigid body. The point load/observation, complete plate mass, tensor
//! rigid inertia and restitution law retain their existing numerical owners.
//!
//! This is a localized event, at fixed configuration, in the flat plate's +z
//! direction. The caller supplies a body material point above the selected
//! plate point. It does not detect collisions, rotate the plate normal, solve
//! simultaneous contacts, apply friction, or deposit impact loss as heat.
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
        Ok(Self {
            dynamics,
            mesh,
            full_mass,
            point,
            unit_impulse,
            full_unit_impulse,
            inverse_plate_mass,
            controls,
            native_workspace,
        })
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
        let mut force = vec![f64::NAN; n];
        self.dynamics
            .forcing(plate.t, &mut force)
            .map_err(PlateRigidImpactError::Load)?;
        poll(cancelled)?;
        // Do not clone potentially long telemetry. The native event only changes
        // v/a; publication below keeps the entire original history in place.
        let mut trial = SecondOrderState::new(plate.t, &plate.q, &plate.v, &plate.a);
        let plate_receipt = second_order_impulse(
            &mut trial,
            &self.dynamics,
            &impulse,
            &force,
            &self.dynamics.mass_preconditioner(),
            self.controls.impulse,
            self.native_workspace,
            cancelled,
        )
        .map_err(PlateRigidImpactError::Plate)?;
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
        let dp = rigid
            .state_after
            .linear_momentum_world()
            .sub(body.linear_momentum_world());
        let dl = rigid
            .state_after
            .angular_momentum_body()
            .sub(body.angular_momentum_body());
        for (actual, commanded) in components(dp)
            .into_iter()
            .zip(components(rigid.impulse_world))
        {
            check(
                "rigid linear transfer",
                actual - commanded,
                self.controls.linear_impulse_tolerance_ns + rel * actual.abs().max(commanded.abs()),
            )?;
        }
        for (actual, commanded) in components(dl)
            .into_iter()
            .zip(components(rigid.angular_impulse_body))
        {
            check(
                "rigid angular transfer",
                actual - commanded,
                self.controls.angular_impulse_tolerance_nms
                    + rel * actual.abs().max(commanded.abs()),
            )?;
        }
        let map = &self.dynamics.model().dof_map;
        let mut full_change = vec![0.0; map.len()];
        for (i, &r) in map.iter().enumerate() {
            if let Some(r) = r {
                full_change[i] = plate_receipt.velocity_change[r];
            }
        }
        let mut momentum = vec![0.0; map.len()];
        self.full_mass.spmv(&full_change, &mut momentum);
        let scale = self.dynamics.parameters().mass_scale;
        let mut support = vec![0.0; map.len()];
        let mut support_impulses = Vec::with_capacity(map.len() - n);
        for (i, value) in momentum.iter_mut().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            *value *= scale;
            if map[i].is_none() {
                support[i] = *value + j * self.full_unit_impulse[i];
                support_impulses.push((i, support[i]));
            }
        }
        let (plate_linear, plate_angular) = self.resultants(&momentum);
        let (support_linear, support_angular) = self.resultants(&support);
        let body_angular = body
            .pose()
            .orientation()
            .rotate_body_to_world(dl)
            .add(body.pose().position_world().cross(dp));
        let linear_residual = plate_linear.add(dp).sub(support_linear);
        let angular_residual = plate_angular.add(body_angular).sub(support_angular);
        for residual in components(linear_residual) {
            check(
                "pair linear balance",
                residual,
                self.controls.linear_impulse_tolerance_ns + rel * j,
            )?;
        }
        for residual in components(angular_residual) {
            check(
                "pair angular balance",
                residual,
                self.controls.angular_impulse_tolerance_nms,
            )?;
        }
        // Accumulate actual changes from BOTH bodies before rounding. Separate
        // large energies or per-body work totals can erase a small impact loss.
        let mut energy = Dd::ZERO;
        for i in 0..n {
            energy = energy
                + (dd(0.5 * plate.v[i]) + dd(0.5 * trial.v[i]))
                    * dd(plate_receipt.momentum_change[i]);
        }
        for (increment, (before, after)) in components(dp).into_iter().zip(
            components(rigid.point_kinematics_before.center_of_mass_velocity_world)
                .into_iter()
                .zip(components(
                    rigid.point_kinematics_after.center_of_mass_velocity_world,
                )),
        ) {
            energy = energy + dd(increment) * (dd(0.5 * before) + dd(0.5 * after));
        }
        for (increment, (before, after)) in components(dl).into_iter().zip(
            components(rigid.point_kinematics_before.angular_velocity_body)
                .into_iter()
                .zip(components(
                    rigid.point_kinematics_after.angular_velocity_body,
                )),
        ) {
            energy = energy + dd(increment) * (dd(0.5 * before) + dd(0.5 * after));
        }
        let energy_residual = (energy + dd(contact.dissipated_energy)).to_f64();
        check(
            "impact energy",
            energy_residual,
            self.controls.energy_tolerance_j,
        )?;
        let report = PlateRigidImpactReport {
            contact,
            plate: Some(plate_receipt),
            body: Some(rigid),
            support_impulses,
            plate_linear_impulse_world: plate_linear,
            plate_angular_impulse_world: plate_angular,
            support_linear_impulse_world: support_linear,
            support_angular_impulse_world: support_angular,
            linear_balance_residual: linear_residual,
            angular_balance_residual: angular_residual,
            relative_normal_speed_after: relative_after,
            kinetic_energy_change_j: energy.to_f64(),
            energy_residual_j: energy_residual,
        };
        poll(cancelled)?;
        plate.v = trial.v;
        plate.a = trial.a;
        *body = rigid.state_after;
        Ok(report)
    }

    fn relative_speed(&self, motion: &PointKinematics, velocity: &[f64]) -> f64 {
        // Subtract common translation before rounding the contact flow. Forming
        // v_com + omega cross r in f64 first can erase real rotational approach
        // and incorrectly classify an approaching contact as an inactive no-op.
        let mut relative = dd(motion.center_of_mass_velocity_world.z)
            + dd(motion.angular_velocity_world.x) * dd(motion.arm_world.y)
            - dd(motion.angular_velocity_world.y) * dd(motion.arm_world.x);
        for (&row, &v) in self.unit_impulse.iter().zip(velocity) {
            relative = relative - dd(row) * dd(v);
        }
        relative.to_f64()
    }

    fn resultants(&self, impulse: &[f64]) -> (Vec3, Vec3) {
        // Rigid transverse fields are w=vz+omega_x*y-omega_y*x,
        // wx=-omega_y and wy=omega_x. Their exact transpose gives Pz,Lx,Ly.
        let mut p = Dd::ZERO;
        let mut lx = Dd::ZERO;
        let mut ly = Dd::ZERO;
        for (node, &(x, y)) in self.mesh.nodes.iter().enumerate() {
            p = p + dd(impulse[3 * node]);
            lx = lx + dd(y) * dd(impulse[3 * node]) + dd(impulse[3 * node + 2]);
            ly = ly - dd(x) * dd(impulse[3 * node]) - dd(impulse[3 * node + 1]);
        }
        (
            Vec3::new(0.0, 0.0, p.to_f64()),
            Vec3::new(lx.to_f64(), ly.to_f64(), 0.0),
        )
    }
}
