//! Native finite-solid and rigid free motion with one physical balance account.
//!
//! The original solid energy/momentum step and rigid endpoint-kick/free-rotation step
//! own every candidate state. Held loads are accounted at their actual material
//! points. Consistent nodal momentum and current coordinates retain the solid's
//! orbital angular momentum; rigid spin is evaluated from both actual attitudes.
//! Numerical energy defects remain signed integration errors, never impact heat.

use crate::finite_solid::{
    FiniteSolidBody, FiniteSolidControls, FiniteSolidError, FiniteSolidStep,
};
use fs_exec::Cx;
use fs_mbd::gravity::{GravityError, RigidGravityField, step_in_field_with_loads};
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt,
    UnitQuaternion, Vec3, Wrench,
};

/// Native solid solve limits and absolute complete-system balance gates.
#[derive(Clone, Copy, Debug)]
pub struct FreeFlightControls {
    pub maximum_rigid_bodies: usize,
    pub solid: FiniteSolidControls,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
}
impl Default for FreeFlightControls {
    fn default() -> Self {
        Self {
            maximum_rigid_bodies: 64,
            solid: FiniteSolidControls::default(),
            linear_tolerance_kg_m_s: 1e-9,
            angular_tolerance_kg_m2_s: 1e-9,
            energy_tolerance_j: 1e-8,
        }
    }
}

/// Native state receipts and actual increments over the same positive interval.
#[derive(Clone, Debug)]
pub struct FreeFlightReport {
    pub duration_s: f64,
    pub solid: FiniteSolidStep,
    pub rigid: Vec<StepReceipt>,
    /// Evaluated by polarization of actual momenta, not large energy subtraction.
    pub kinetic_change_j: f64,
    pub strain_change_j: f64,
    /// Solid load work, rigid force-displacement work and endpoint torque work.
    /// Uniform rigid gravity is included here; it is not counted again as storage.
    pub source_work_j: f64,
    /// Held rigid COM-force displacement work and body-couple kick work.
    pub rigid_held_work_j: f64,
    /// Original rigid field potential decrease (uniform gravity uses its exact
    /// displacement work). Included once in source work, not again as storage.
    pub rigid_gravity_work_j: f64,
    /// Stable solid kinetic + original strain change minus its native source work.
    pub solid_integration_defect_j: f64,
    /// Actual rotational energy change inside the native torque-free rigid maps.
    pub rigid_drift_energy_defect_j: f64,
    /// Gravity half-kick work minus the same field's potential decrease.
    pub rigid_gravity_work_defect_j: f64,
    /// Held COM-force kick work minus its actual displacement work. Spatial
    /// gravity can make the two endpoint linear half-kicks unequal.
    pub rigid_held_force_work_defect_j: f64,
    pub integration_defect_j: f64,
    /// Sum of absolute solid, individual rigid drift, gravity and held-force
    /// defects. Opposite signed errors cannot pass the integration gate.
    pub absolute_integration_defect_j: f64,
    /// Actual storage change minus source work and the labeled numerical defects.
    pub work_energy_residual_j: f64,
    pub linear_momentum_change_kg_m_s: Vec3,
    pub angular_momentum_change_kg_m2_s: Vec3,
    /// Held nodal/body loads, the solid's native midpoint gravity, and both
    /// native rigid force/torque endpoint kicks.
    pub external_impulse_kg_m_s: Vec3,
    pub external_angular_impulse_kg_m2_s: Vec3,
    /// External reaction ON the solid's stationary prescribed DOFs.
    pub support_impulse_kg_m_s: Vec3,
    pub support_angular_impulse_kg_m2_s: Vec3,
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Sum of contributor residual norms, including the supported solid and
    /// each rigid body. Opposite unresolved transfers cannot cancel globally.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    /// Sum of absolute individual rigid storage-minus-work-minus-drift errors.
    pub absolute_work_residual_j: f64,
}

#[derive(Debug)]
pub enum FreeFlightError {
    Invalid(&'static str),
    Solid(FiniteSolidError),
    Rigid(DynamicsError),
    Gravity(GravityError),
    Cancelled,
    /// Neither mechanical candidate is published after a failed physical gate.
    Audit(Box<FreeFlightReport>),
}
impl From<FiniteSolidError> for FreeFlightError {
    fn from(error: FiniteSolidError) -> Self {
        Self::Solid(error)
    }
}
impl From<DynamicsError> for FreeFlightError {
    fn from(error: DynamicsError) -> Self {
        Self::Rigid(error)
    }
}
impl From<GravityError> for FreeFlightError {
    fn from(error: GravityError) -> Self {
        match error {
            GravityError::Cancelled => Self::Cancelled,
            error => Self::Gravity(error),
        }
    }
}
impl core::fmt::Display for FreeFlightError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid/finite-solid free flight: {self:?}")
    }
}
impl std::error::Error for FreeFlightError {}

/// One declared native rigid gravity law for every trial in the shared driver.
#[derive(Clone, Copy)]
pub(crate) enum RigidGravity<'a> {
    Uniform(Gravity),
    Spatial {
        field: &'a dyn RigidGravityField,
        frame: u64,
    },
}

impl<'a> FiniteSolidBody<'a> {
    /// Advance both native owners under held loads and admit actual balances.
    ///
    /// The solid receives its explicit referential nodal/body loads and any
    /// retained StationaryGravity field. The `gravity` argument belongs to the
    /// rigid integrator. Held loads must not duplicate the solid's bound field;
    /// its potential work is already included once in the native source work.
    ///
    /// A wrench holds its world force and body torque at BOTH endpoints. Applied
    /// angular impulse uses those two rotated body torques and actual endpoint
    /// force moments. This is the native discrete torque convention, not an
    /// exact forced Euler-top flow. The solid uses its original objective
    /// energy/momentum discrete gradient, including the same algorithmic support
    /// reaction integrated at midpoint coordinates; fixed supports do zero work.
    ///
    /// The signed solid and rigid-drift energy defects are separately bounded
    /// through their absolute sum. Neither this gate nor free motion admits a
    /// collision path or proves trajectory accuracy; the enclosing event driver
    /// must admit geometry and perform timestep refinement. Pressure-bound
    /// solids retain their original explicit-pressure refusal on this entry.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn rigid_free_flight_candidate(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FreeFlightControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<RigidBodyState>, FreeFlightReport), FreeFlightError> {
        self.rigid_free_flight_with_source(
            bodies,
            properties,
            nodal_forces_n,
            body_force_n_m3,
            None,
            wrenches,
            RigidGravity::Uniform(gravity),
            duration_s,
            origin_m,
            controls,
            cx,
        )
    }

    /// Advance both owners with uniform pressure held on the solid's original
    /// closed material boundary. The solid uses its objective discrete gradient
    /// of material energy plus pressure-volume potential; native -p Delta V work
    /// remains in `report.solid.pressure_work_j` and total source work. Closed
    /// uniform pressure adds no net external linear/angular impulse. Stationary
    /// support reactions and all existing phase balance gates remain active.
    /// This does not add nonuniform pressure or resolve a contact path.
    #[allow(clippy::too_many_arguments)]
    pub fn rigid_free_flight_candidate_with_pressure(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_pa: f64,
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FreeFlightControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<RigidBodyState>, FreeFlightReport), FreeFlightError> {
        self.rigid_free_flight_with_source(
            bodies,
            properties,
            nodal_forces_n,
            body_force_n_m3,
            Some(pressure_pa),
            wrenches,
            RigidGravity::Uniform(gravity),
            duration_s,
            origin_m,
            controls,
            cx,
        )
    }

    /// Advance rigid participants in one explicitly declared stationary native
    /// field alongside the solid's own retained sources. The original rigid
    /// gravity owner samples force, torque and potential at both actual poses;
    /// its monopole or explicitly selected finite-body model is unchanged.
    ///
    /// `uniform_pressure_pa` supplies held pressure only to a pressure-bound
    /// solid. Held wrenches must exclude this field. Potential decrease is
    /// source work once; gravity, held-force and free-rotation numerical defects
    /// have separate noncancelling gates. Endpoint admission does not certify a
    /// field domain along the intervening path. Moving fields and sustained
    /// loaded contact require their own stage/trajectory owner.
    #[allow(clippy::too_many_arguments)]
    pub fn rigid_free_flight_candidate_in_field(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        uniform_pressure_pa: Option<f64>,
        wrenches: &[Wrench],
        field: &dyn RigidGravityField,
        frame: u64,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FreeFlightControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<RigidBodyState>, FreeFlightReport), FreeFlightError> {
        self.rigid_free_flight_with_source(
            bodies,
            properties,
            nodal_forces_n,
            body_force_n_m3,
            uniform_pressure_pa,
            wrenches,
            RigidGravity::Spatial { field, frame },
            duration_s,
            origin_m,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn rigid_free_flight_with_source(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        uniform_pressure_pa: Option<f64>,
        wrenches: &[Wrench],
        gravity: RigidGravity<'_>,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FreeFlightControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<RigidBodyState>, FreeFlightReport), FreeFlightError> {
        poll(cx)?;
        let origin = vector(origin_m);
        if !origin.is_finite()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || controls.maximum_rigid_bodies == 0
            || bodies.len() > controls.maximum_rigid_bodies
            || bodies.len() != properties.len()
            || bodies.len() != wrenches.len()
            || [
                controls.linear_tolerance_kg_m_s,
                controls.angular_tolerance_kg_m2_s,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(FreeFlightError::Invalid(
                "bounded matching bodies, positive interval and finite physical controls required",
            ));
        }
        admit_uniform_pressure(self, uniform_pressure_pa, cx)?;
        let (solid_next, solid) = match uniform_pressure_pa {
            Some(pressure_pa) => self.step_energy_momentum_candidate_with_pressure(
                duration_s,
                nodal_forces_n,
                body_force_n_m3,
                pressure_pa,
                controls.solid,
                cx,
            ),
            None => self.step_energy_momentum_candidate(
                duration_s,
                nodal_forces_n,
                body_force_n_m3,
                controls.solid,
                cx,
            ),
        }?;
        // Same prepared load owner as the native step: no second body-force
        // quadrature or lumped mass enters the momentum account. Uniform pressure
        // has zero total force/torque; its actual work and any constrained
        // support impulse already belong to the native solid receipt.
        let mut force = self
            .prepared
            .external_force_n(nodal_forces_n, body_force_n_m3, cx)
            .map_err(FiniteSolidError::from)?;
        // The solid may retain its own conservative point field. Its actual
        // native midpoint force must enter these impulses and moments exactly
        // once, including supported nodes. Its potential work is already in
        // solid.source_work_j; the rigid gravity argument remains independent.
        if let Some(gravity) = &solid.gravity {
            for (dof, (force, &field_force)) in
                force.iter_mut().zip(&gravity.midpoint_force_n).enumerate()
            {
                if dof % 4096 == 0 {
                    poll(cx)?;
                }
                *force += field_force;
                if !force.is_finite() {
                    return Err(FreeFlightError::Invalid(
                        "solid source force is not representable",
                    ));
                }
            }
        }
        let before_x = self.current_positions(cx)?;
        let after_x = solid_next.current_positions(cx)?;
        let mut delta_v = vec![0.0; force.len()];
        let mut average_v = vec![0.0; force.len()];
        for (row, &dof) in self.free.iter().enumerate() {
            if row % 4096 == 0 {
                poll(cx)?;
            }
            delta_v[dof] = solid_next.state.v[row] - self.state.v[row];
            average_v[dof] = 0.5 * self.state.v[row] + 0.5 * solid_next.state.v[row];
        }
        let mut delta_p = vec![0.0; force.len()];
        let mut average_p = vec![0.0; force.len()];
        self.mass.spmv(&delta_v, &mut delta_p);
        self.mass.spmv(&average_v, &mut average_p);
        let mut kinetic_change_j = 0.0;
        let mut linear_change = Vec3::ZERO;
        let mut angular_change = Vec3::ZERO;
        let mut applied = Vec3::ZERO;
        let mut applied_angular = Vec3::ZERO;
        for node in 0..before_x.len() {
            if node % 256 == 0 {
                poll(cx)?;
            }
            let r = midpoint_arm(before_x[node], after_x[node], origin);
            let dx = vector(after_x[node]).sub(vector(before_x[node]));
            let dp = triple(&delta_p, node);
            linear_change = linear_change.add(dp);
            // Exact bilinear identity for x1×p1-x0×p0, including mesh motion.
            angular_change = angular_change
                .add(r.cross(dp))
                .add(dx.cross(triple(&average_p, node)));
            let impulse = triple(&force, node).scale(duration_s);
            applied = applied.add(impulse);
            applied_angular = applied_angular.add(r.cross(impulse));
            for c in 0..3 {
                let dof = 3 * node + c;
                kinetic_change_j = average_v[dof].mul_add(delta_p[dof], kinetic_change_j);
            }
        }
        let mut support = Vec3::ZERO;
        let mut support_angular = Vec3::ZERO;
        for &(dof, reaction) in &solid.support_reactions_n {
            poll(cx)?;
            let mut components = [0.0; 3];
            components[dof % 3] = reaction * duration_s;
            let impulse = vector(components);
            support = support.add(impulse);
            support_angular = support_angular
                .add(midpoint_arm(before_x[dof / 3], after_x[dof / 3], origin).cross(impulse));
        }
        let strain_change_j = solid.strain_after_j - solid.strain_before_j;
        let solid_defect = kinetic_change_j + strain_change_j - solid.source_work_j;
        let mut source_work_j = solid.source_work_j;
        let mut rigid_held_work_j = 0.0;
        let mut rigid_gravity_work_j = 0.0;
        let mut drift_defect = 0.0;
        let mut gravity_defect = 0.0;
        let mut held_force_defect = 0.0;
        let mut absolute_defect = solid_defect.abs();
        let mut absolute_linear_residual = magnitude(linear_change.sub(applied).sub(support));
        let mut absolute_angular_residual =
            magnitude(angular_change.sub(applied_angular).sub(support_angular));
        let mut absolute_work_residual = 0.0;
        let mut rigid_next = Vec::with_capacity(bodies.len());
        let mut rigid = Vec::with_capacity(bodies.len());
        for ((&old, &property), &wrench) in bodies.iter().zip(properties).zip(wrenches) {
            let flight = match gravity {
                RigidGravity::Uniform(gravity) => {
                    rigid_free_flight(old, property, wrench, gravity, duration_s, origin, cx)?
                }
                _ => rigid_free_flight_with_source(
                    old, property, wrench, gravity, duration_s, origin, cx,
                )?,
            };
            kinetic_change_j += flight.kinetic_change_j;
            source_work_j += flight.source_work_j;
            rigid_held_work_j += flight.held_work_j;
            rigid_gravity_work_j += flight.gravity_work_j;
            drift_defect += flight.drift_defect_j;
            gravity_defect += flight.gravity_defect_j;
            held_force_defect += flight.held_force_defect_j;
            absolute_defect += flight.drift_defect_j.abs()
                + flight.gravity_defect_j.abs()
                + flight.held_force_defect_j.abs();
            linear_change = linear_change.add(flight.linear_change);
            angular_change = angular_change.add(flight.angular_change);
            applied = applied.add(flight.applied);
            applied_angular = applied_angular.add(flight.applied_angular);
            absolute_linear_residual += magnitude(flight.linear_residual);
            absolute_angular_residual += magnitude(flight.angular_residual);
            absolute_work_residual += flight.work_residual_j.abs();
            rigid_next.push(flight.receipt.state_after);
            rigid.push(flight.receipt);
        }
        let integration_defect_j = solid_defect + drift_defect + gravity_defect + held_force_defect;
        let report = FreeFlightReport {
            duration_s,
            solid,
            rigid,
            kinetic_change_j,
            strain_change_j,
            source_work_j,
            rigid_held_work_j,
            rigid_gravity_work_j,
            solid_integration_defect_j: solid_defect,
            rigid_drift_energy_defect_j: drift_defect,
            rigid_gravity_work_defect_j: gravity_defect,
            rigid_held_force_work_defect_j: held_force_defect,
            integration_defect_j,
            absolute_integration_defect_j: absolute_defect,
            work_energy_residual_j: kinetic_change_j + strain_change_j
                - source_work_j
                - integration_defect_j,
            linear_momentum_change_kg_m_s: linear_change,
            angular_momentum_change_kg_m2_s: angular_change,
            external_impulse_kg_m_s: applied,
            external_angular_impulse_kg_m2_s: applied_angular,
            support_impulse_kg_m_s: support,
            support_angular_impulse_kg_m2_s: support_angular,
            linear_momentum_residual_kg_m_s: linear_change.sub(applied).sub(support),
            angular_momentum_residual_kg_m2_s: angular_change
                .sub(applied_angular)
                .sub(support_angular),
            absolute_linear_residual_kg_m_s: absolute_linear_residual,
            absolute_angular_residual_kg_m2_s: absolute_angular_residual,
            absolute_work_residual_j: absolute_work_residual,
        };
        if [
            report.kinetic_change_j,
            report.strain_change_j,
            report.source_work_j,
            report.rigid_held_work_j,
            report.rigid_gravity_work_j,
            report.rigid_drift_energy_defect_j,
            report.rigid_gravity_work_defect_j,
            report.rigid_held_force_work_defect_j,
            report.integration_defect_j,
            report.absolute_integration_defect_j,
            report.work_energy_residual_j,
            report.absolute_linear_residual_kg_m_s,
            report.absolute_angular_residual_kg_m2_s,
            report.absolute_work_residual_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || [
                linear_change,
                angular_change,
                applied,
                applied_angular,
                support,
                support_angular,
                report.linear_momentum_residual_kg_m_s,
                report.angular_momentum_residual_kg_m2_s,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || magnitude(report.linear_momentum_residual_kg_m_s) > controls.linear_tolerance_kg_m_s
            || magnitude(report.angular_momentum_residual_kg_m2_s)
                > controls.angular_tolerance_kg_m2_s
            || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
            || report.absolute_integration_defect_j > controls.energy_tolerance_j
            || report.absolute_linear_residual_kg_m_s > controls.linear_tolerance_kg_m_s
            || report.absolute_angular_residual_kg_m2_s > controls.angular_tolerance_kg_m2_s
            || report.absolute_work_residual_j > controls.energy_tolerance_j
        {
            return Err(FreeFlightError::Audit(Box::new(report)));
        }
        poll(cx)?;
        Ok((solid_next, rigid_next, report))
    }
}

/// One original rigid endpoint-kick/free-rotation step and its actual account.
/// Consumers sum these receipts in body order alongside their native solid.
#[derive(Clone, Debug)]
pub(crate) struct RigidFreeFlight {
    pub receipt: StepReceipt,
    pub kinetic_change_j: f64,
    pub source_work_j: f64,
    pub held_work_j: f64,
    pub gravity_work_j: f64,
    pub drift_defect_j: f64,
    pub gravity_defect_j: f64,
    pub held_force_defect_j: f64,
    pub linear_change: Vec3,
    pub angular_change: Vec3,
    pub applied: Vec3,
    pub applied_angular: Vec3,
    pub linear_residual: Vec3,
    pub angular_residual: Vec3,
    pub work_residual_j: f64,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn rigid_free_flight(
    old: RigidBodyState,
    property: MassProperties,
    wrench: Wrench,
    gravity: Gravity,
    duration_s: f64,
    origin: Vec3,
    cx: &Cx<'_>,
) -> Result<RigidFreeFlight, FreeFlightError> {
    rigid_free_flight_with_source(
        old,
        property,
        wrench,
        RigidGravity::Uniform(gravity),
        duration_s,
        origin,
        cx,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn rigid_free_flight_with_source(
    old: RigidBodyState,
    property: MassProperties,
    wrench: Wrench,
    gravity: RigidGravity<'_>,
    duration_s: f64,
    origin: Vec3,
    cx: &Cx<'_>,
) -> Result<RigidFreeFlight, FreeFlightError> {
    poll(cx)?;
    let half = 0.5 * duration_s;
    let (receipt, drifted, gravity_before, gravity_after, potential_work) = match gravity {
        RigidGravity::Uniform(gravity) => {
            // Preserve the original uniform native map and its exact constant
            // force work. No subtraction of a large potential datum is needed.
            let (receipt, drifted) = RigidBodyIntegrator::new(gravity).step_verlet_with(
                old,
                property,
                wrench,
                duration_s,
                |drifted| {
                    poll(cx)?;
                    Ok::<_, FreeFlightError>((wrench, drifted))
                },
            )?;
            let load = Wrench {
                force_world: gravity.acceleration_world().scale(property.mass()),
                torque_body: Vec3::ZERO,
            };
            let work = load.force_world.dot(
                receipt
                    .state_after
                    .pose()
                    .position_world()
                    .sub(old.pose().position_world()),
            );
            (receipt, drifted, load, load, work)
        }
        RigidGravity::Spatial { field, frame } => {
            let native =
                step_in_field_with_loads(field, frame, old, property, wrench, duration_s, || {
                    cx.checkpoint().is_err()
                })?;
            (
                native.step,
                native.state_after_drift,
                native.load_before.wrench,
                native.load_after.wrench,
                native.load_before.potential_energy_j - native.load_after.potential_energy_j,
            )
        }
    };
    let new = receipt.state_after;
    let force_before = wrench.force_world.add(gravity_before.force_world);
    let force_after = wrench.force_world.add(gravity_after.force_world);
    let torque_before_half = wrench
        .torque_body
        .add(gravity_before.torque_body)
        .scale(half);
    let torque_after_half = wrench
        .torque_body
        .add(gravity_after.torque_body)
        .scale(half);
    // Exactly the original first half-kick arithmetic; accepted momenta/poses
    // are never reconstructed or projected by this shared account.
    let kicked = old.angular_momentum_body().add(torque_before_half);
    let drift_defect_j = rotational_change(property, kicked, drifted.angular_momentum_body())?;
    let first_omega = angular_average(property, old.angular_momentum_body(), kicked)?;
    let last_omega = angular_average(
        property,
        drifted.angular_momentum_body(),
        new.angular_momentum_body(),
    )?;
    let torque_work = wrench
        .torque_body
        .scale(half)
        .dot(first_omega.add(last_omega));
    let dx = new.pose().position_world().sub(old.pose().position_world());
    let dp = new.linear_momentum_world().sub(old.linear_momentum_world());
    let p_average = old
        .linear_momentum_world()
        .scale(0.5)
        .add(new.linear_momentum_world().scale(0.5));
    let kinetic_change_j = p_average.scale(property.mass().recip()).dot(dp)
        + rotational_change(
            property,
            old.angular_momentum_body(),
            new.angular_momentum_body(),
        )?;
    let held_force_work = wrench.force_world.dot(dx);
    let held_work_j = held_force_work + torque_work;
    // Uniform-only consumers retain the original constant-force account,
    // including its roundoff convention and native drift-only defect.
    let source_work_j = match gravity {
        RigidGravity::Uniform(_) => force_before.dot(dx) + torque_work,
        RigidGravity::Spatial { .. } => held_work_j + potential_work,
    };
    // These are commanded native half-kicks paired with the ACTUAL adjoining
    // momenta. Component decomposition error is checked in work_residual_j;
    // it cannot be silently assigned to gravity or free-rotation energy.
    let first_velocity = old
        .linear_momentum_world()
        .scale(0.5 / property.mass())
        .add(drifted.linear_momentum_world().scale(0.5 / property.mass()));
    let last_velocity = drifted
        .linear_momentum_world()
        .scale(0.5 / property.mass())
        .add(new.linear_momentum_world().scale(0.5 / property.mass()));
    let gravity_kick_work = gravity_before.force_world.scale(half).dot(first_velocity)
        + gravity_after.force_world.scale(half).dot(last_velocity)
        + gravity_before.torque_body.scale(half).dot(first_omega)
        + gravity_after.torque_body.scale(half).dot(last_omega);
    let held_force_kick_work = wrench
        .force_world
        .scale(half)
        .dot(first_velocity.add(last_velocity));
    let (gravity_defect_j, held_force_defect_j) = match gravity {
        RigidGravity::Uniform(_) => (0.0, 0.0),
        RigidGravity::Spatial { .. } => (
            gravity_kick_work - potential_work,
            held_force_kick_work - held_force_work,
        ),
    };
    let r_average = old
        .pose()
        .position_world()
        .sub(origin)
        .scale(0.5)
        .add(new.pose().position_world().sub(origin).scale(0.5));
    let dl = new.angular_momentum_body().sub(old.angular_momentum_body());
    let angular_change = r_average
        .cross(dp)
        .add(dx.cross(p_average))
        .add(old.pose().orientation().rotate_body_to_world(dl))
        .add(rotation_difference(
            old.pose().orientation(),
            new.pose().orientation(),
            new.angular_momentum_body(),
        ));
    let force_before_half = force_before.scale(half);
    let force_after_half = force_after.scale(half);
    let applied = force_before_half.add(force_after_half);
    let applied_angular = old
        .pose()
        .position_world()
        .sub(origin)
        .cross(force_before_half)
        .add(
            new.pose()
                .position_world()
                .sub(origin)
                .cross(force_after_half),
        )
        .add(
            old.pose()
                .orientation()
                .rotate_body_to_world(torque_before_half),
        )
        .add(
            new.pose()
                .orientation()
                .rotate_body_to_world(torque_after_half),
        );
    Ok(RigidFreeFlight {
        receipt,
        kinetic_change_j,
        source_work_j,
        held_work_j,
        gravity_work_j: potential_work,
        drift_defect_j,
        gravity_defect_j,
        held_force_defect_j,
        linear_change: dp,
        angular_change,
        applied,
        applied_angular,
        linear_residual: dp.sub(applied),
        angular_residual: angular_change.sub(applied_angular),
        work_residual_j: kinetic_change_j
            - source_work_j
            - drift_defect_j
            - gravity_defect_j
            - held_force_defect_j,
    })
}

// Pressure admission precedes contact work in the trajectory consumer. The
// original bound surface still owns closure/reference identity and finite load.
pub(super) fn admit_uniform_pressure(
    solid: &FiniteSolidBody<'_>,
    pressure_pa: Option<f64>,
    cx: &Cx<'_>,
) -> Result<(), FreeFlightError> {
    if let Some(pressure_pa) = pressure_pa {
        let surface = solid
            .pressure_surface
            .as_deref()
            .ok_or(FreeFlightError::Solid(FiniteSolidError::Invalid(
                "body was not initialized with a pressure surface",
            )))?;
        surface
            .validate_pressure(pressure_pa.into(), cx)
            .map_err(FiniteSolidError::from)?;
    }
    Ok(())
}

fn angular_average(
    property: MassProperties,
    before: Vec3,
    after: Vec3,
) -> Result<Vec3, DynamicsError> {
    Ok(property
        .angular_velocity_body_checked(before)?
        .scale(0.5)
        .add(property.angular_velocity_body_checked(after)?.scale(0.5)))
}
pub(super) fn rotational_change(
    property: MassProperties,
    before: Vec3,
    after: Vec3,
) -> Result<f64, DynamicsError> {
    Ok(angular_average(property, before, after)?.dot(after.sub(before)))
}

/// Polarize q(0,v)q* instead of subtracting two large rotated spins. This uses
/// actual quaternion components and keeps their representable norm residual;
/// it neither renormalizes nor assumes the native rotation conserved momentum.
pub(super) fn rotation_difference(before: UnitQuaternion, after: UnitQuaternion, v: Vec3) -> Vec3 {
    let a = before.components();
    let mut b = after.components();
    if a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
        b = b.map(|x| -x);
    }
    let mean = Vec3::new(
        0.5 * a[1] + 0.5 * b[1],
        0.5 * a[2] + 0.5 * b[2],
        0.5 * a[3] + 0.5 * b[3],
    );
    let delta = Vec3::new(b[1] - a[1], b[2] - a[2], b[3] - a[3]);
    let mean_w = 0.5 * a[0] + 0.5 * b[0];
    let delta_w = b[0] - a[0];
    v.scale(2.0 * (mean_w * delta_w - mean.dot(delta)))
        .add(delta.scale(2.0 * mean.dot(v)))
        .add(mean.scale(2.0 * delta.dot(v)))
        .add(mean.cross(v).scale(2.0 * delta_w))
        .add(delta.cross(v).scale(2.0 * mean_w))
}
fn midpoint_arm(a: [f64; 3], b: [f64; 3], origin: Vec3) -> Vec3 {
    vector(a)
        .sub(origin)
        .scale(0.5)
        .add(vector(b).sub(origin).scale(0.5))
}
fn triple(values: &[f64], node: usize) -> Vec3 {
    Vec3::new(values[3 * node], values[3 * node + 1], values[3 * node + 2])
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn poll(cx: &Cx<'_>) -> Result<(), FreeFlightError> {
    cx.checkpoint().map_err(|_| FreeFlightError::Cancelled)
}
