//! Finite-inertia feedback between native planar induction and one supported rotor.
//!
//! This composition represents a homogeneous isotropic circular/annular conductor
//! whose analytic material occupancy is invariant under rotation about z. The
//! fixed triangular mesh approximates its circles within an admitted chord error;
//! it is not a moving polygon, sliding mesh, or a moving material interface.
//! The complete field domain has uniform permeability, so omitted permeability
//! interface traction cannot masquerade as a transverse bearing-force balance.
//!
//! The native field midpoint and native constrained inverse-inertia map share
//! the same angular velocity after a bounded scalar iteration. No rejected field
//! trial advances history. Joule heat, shaft work, support impulse and numerical
//! work defects remain distinct; no state or energy account is repaired.

use fs_em::planar::eddy::{
    EddyError, PlanarEddy,
    moving::{
        MovingEddyControls, MovingEddyStep,
        mechanical::{PlanarLorentzLoadError, PlanarLorentzLoadLimits},
    },
};
use fs_ga::{GaError, Se3, So3, Twist, Vec3, Wrench};
use fs_math::{dd::Dd, det};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointKind, PrescribedBaseImpulseResponse,
    apply_prescribed_base_impulses, forward_kinematics,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

/// Analytic, rotationally invariant conductor section. The axis position comes
/// from the original link's zero transform, avoiding another position owner.
#[derive(Clone, Copy, Debug)]
pub struct CircularConductor {
    /// Zero denotes a solid disk; a positive value denotes an annulus [m].
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    /// Maximum radial departure of every boundary chord from its circle [m].
    /// This is geometry admission, not a field-error or convergence certificate.
    pub maximum_boundary_error_m: f64,
}

/// Coupling convergence and physical-accounting allowances, independent of the
/// existing field Krylov/energy/Peclet controls.
#[derive(Clone, Copy, Debug)]
pub struct RotorControls {
    pub maximum_iterations: usize,
    pub relaxation: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub energy_tolerance_j: f64,
    pub angular_impulse_tolerance_nms: f64,
    pub linear_impulse_tolerance_ns: f64,
}

#[derive(Debug)]
pub enum RotorError {
    Invalid(&'static str),
    Cancelled,
    Field(EddyError),
    Mechanical(ArticulatedError),
    Load(PlanarLorentzLoadError),
    Geometry(GaError),
    Nonconvergence {
        iterations: usize,
        velocity_residual_rad_s: f64,
    },
    Balance(&'static str),
}
impl core::fmt::Display for RotorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar induction rotor refused: {self:?}")
    }
}
impl std::error::Error for RotorError {}
impl From<EddyError> for RotorError {
    fn from(value: EddyError) -> Self {
        Self::Field(value)
    }
}
impl From<ArticulatedError> for RotorError {
    fn from(value: ArticulatedError) -> Self {
        Self::Mechanical(value)
    }
}
impl From<PlanarLorentzLoadError> for RotorError {
    fn from(value: PlanarLorentzLoadError) -> Self {
        Self::Load(value)
    }
}
impl From<GaError> for RotorError {
    fn from(value: GaError) -> Self {
        Self::Geometry(value)
    }
}
type Result<T> = core::result::Result<T, RotorError>;

/// One accepted field/mechanical history. The original model and field owner
/// are immutable through this API; accepting a candidate replaces both together.
pub struct PlanarInductionRotor<'m> {
    field: PlanarEddy<'m>,
    binding: Arc<RotorBinding>,
    angle_rad: f64,
    angular_velocity_rad_s: f64,
}

// Geometry and numerical-owner binding survive accepted steps by shared
// ownership; publishing a new history never clones the mesh-sized arrays.
struct RotorBinding {
    model: ArticulatedModel,
    triangles: Vec<usize>,
    unit_velocity: Vec<[f64; 2]>,
    load_limits: PlanarLorentzLoadLimits,
    boundary_error_m: f64,
}

/// Diagnostics of the actual accepted joint/field midpoint, never predictor loads.
#[derive(Clone, Copy, Debug)]
pub struct RotorBalance {
    pub iterations: usize,
    pub field_angular_velocity_rad_s: f64,
    pub mechanical_midpoint_velocity_rad_s: f64,
    pub velocity_residual_rad_s: f64,
    pub electromagnetic_torque_nm: f64,
    /// Actual midpoint speed times the held applied shaft torque and duration.
    pub shaft_work_j: f64,
    /// Native P1 field work minus duration * field speed * integrated torque.
    pub load_work_residual_j: f64,
    /// Actual kinetic increment minus native Lorentz work and shaft work.
    pub mechanical_work_residual_j: f64,
    /// Magnetic + kinetic change + heat - field source work - shaft work.
    pub total_energy_residual_j: f64,
}

pub struct RotorCandidate<'m> {
    state: PlanarInductionRotor<'m>,
    pub induction: MovingEddyStep,
    /// Native full body impulse and ideal-bearing reaction, not only z torque.
    pub mechanics: PrescribedBaseImpulseResponse,
    pub balance: RotorBalance,
}
impl<'m> RotorCandidate<'m> {
    pub fn state(&self) -> &PlanarInductionRotor<'m> {
        &self.state
    }
    pub fn into_state(self) -> PlanarInductionRotor<'m> {
        self.state
    }
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RotorError::Invalid("nonfinite rotor arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(RotorError::Cancelled)
    } else {
        Ok(())
    }
}
fn base() -> BaseState {
    BaseState::stationary(Se3::identity())
}

impl<'m> PlanarInductionRotor<'m> {
    /// Bind one original revolute-z root link and every conducting cell to a
    /// homogeneous circular rotor. The link origin is the circle axis at the
    /// extrusion mid-depth. Its COM must lie at that origin and its inertia
    /// must be axisymmetric there. Additional stator conductors, impressed
    /// current inside the rotor, varying material sources and permeability
    /// interfaces refuse; insulating impressed-current coils remain available.
    ///
    /// The boundary check rejects finite square/off-axis/sector interpretations.
    /// Its chord approximation does not certify exact material no-penetration
    /// on the polygon; refining the fixed circle mesh remains the caller's job.
    pub fn new(
        field: PlanarEddy<'m>,
        model: ArticulatedModel,
        section: CircularConductor,
        angle_rad: f64,
        angular_velocity_rad_s: f64,
        limits: PlanarLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let problem = field.problem();
        let mesh = problem.mesh();
        if model.link_count() != 1
            || model.dof_count() != 1
            || mesh.vertices().len() > limits.max_vertices
            || mesh.faces().len() > limits.max_cells
            || field.conductors().is_empty()
            || field.conductors().len() > limits.max_cells
        {
            return Err(RotorError::Invalid(
                "one bounded native rotor and nonempty conductor required",
            ));
        }
        let link = &model.links()[0];
        let axis = link.joint().motion_subspace();
        let origin = link.parent_from_child_zero().translation();
        let depth = match mesh.metric() {
            fs_rep_mesh::Metric2::Planar { thickness } => thickness,
            _ => return Err(RotorError::Invalid("planar extrusion required")),
        };
        let inertia = link.inertia();
        let tensor = inertia.inertia_com().m;
        if link.joint().kind() != JointKind::Revolute
            || axis != Twist::new(Vec3::new(0.0, 0.0, 1.0), Vec3::default())
            || link.parent_from_child_zero().rotation() != So3::identity()
            || origin.z != 0.5 * depth
            || inertia.center_of_mass() != Vec3::default()
            || tensor[0] != tensor[4]
            || [1, 2, 3, 5, 6, 7].into_iter().any(|i| tensor[i] != 0.0)
        {
            return Err(RotorError::Invalid(
                "centered axisymmetric native inertia and aligned revolute-z root required",
            ));
        }
        forward_kinematics(&model, base(), &[angle_rad], &[angular_velocity_rad_s])?;
        let first = &field.conductors()[0];
        let mu = problem.permeability_h_m()[first.triangle];
        if problem.permeability_h_m().iter().any(|value| *value != mu)
            || field.conductors().iter().any(|cell| {
                cell.region != first.region
                    || cell.conductor != first.conductor
                    || field.impressed_current_density_a_m2()[cell.triangle] != 0.0
            })
        {
            return Err(RotorError::Invalid(
                "uniform permeability and one homogeneous undriven conductor source required",
            ));
        }
        let triangles: Vec<_> = field.conductors().iter().map(|c| c.triangle).collect();
        let boundary_error_m = admit_circle(&field, &triangles, origin, section, cancelled)?;
        let unit_velocity = mesh
            .vertices()
            .iter()
            .enumerate()
            .map(|(i, point)| {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                Ok([finite(origin.y - point[1])?, finite(point[0] - origin.x)?])
            })
            .collect::<Result<Vec<_>>>()?;
        poll(cancelled)?;
        Ok(Self {
            field,
            binding: Arc::new(RotorBinding {
                model,
                triangles,
                unit_velocity,
                load_limits: limits,
                boundary_error_m,
            }),
            angle_rad,
            angular_velocity_rad_s,
        })
    }

    pub fn field(&self) -> &PlanarEddy<'m> {
        &self.field
    }
    pub fn model(&self) -> &ArticulatedModel {
        &self.binding.model
    }
    pub fn angle_rad(&self) -> f64 {
        self.angle_rad
    }
    pub fn angular_velocity_rad_s(&self) -> f64 {
        self.angular_velocity_rad_s
    }
    pub fn boundary_error_m(&self) -> f64 {
        self.binding.boundary_error_m
    }

    /// Implicit midpoint feedback using the original field solve and original
    /// ABA impulse map. Shaft torque is an explicit held external drive on the
    /// relative coordinate; electromagnetic force remains a full body wrench.
    /// All trials start from this SAME accepted field and velocity. Iteration,
    /// work/impulse budget or cancellation failure returns no partial history.
    pub fn step_candidate(
        &self,
        target_s: f64,
        impressed_jz_a_m2: &[f64],
        shaft_torque_nm: f64,
        field_controls: MovingEddyControls,
        controls: RotorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RotorCandidate<'m>> {
        poll(cancelled)?;
        let duration = finite(target_s - self.field.time_s())?;
        if duration <= 0.0
            || !shaft_torque_nm.is_finite()
            || !(1..=256).contains(&controls.maximum_iterations)
            || !controls.relaxation.is_finite()
            || controls.relaxation <= 0.0
            || controls.relaxation > 1.0
            || [
                controls.angular_velocity_tolerance_rad_s,
                controls.energy_tolerance_j,
                controls.angular_impulse_tolerance_nms,
                controls.linear_impulse_tolerance_ns,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || impressed_jz_a_m2.len() != self.field.problem().mesh().faces().len()
            || self
                .binding
                .triangles
                .iter()
                .any(|&i| impressed_jz_a_m2[i] != 0.0)
        {
            return Err(RotorError::Invalid(
                "positive finite coupling controls, interval and external-only impressed current required",
            ));
        }
        self.binding
            .model
            .validate_applied_effort(&[shaft_torque_nm])?;
        let old_speed = self.angular_velocity_rad_s;
        let origin = self.binding.model.links()[0]
            .parent_from_child_zero()
            .translation();
        let mut speed = old_speed;
        let mut velocity = vec![[0.0; 2]; self.binding.unit_velocity.len()];
        let mut last_residual = f64::INFINITY;
        for iteration in 1..=controls.maximum_iterations {
            poll(cancelled)?;
            for (i, (output, unit)) in velocity
                .iter_mut()
                .zip(&self.binding.unit_velocity)
                .enumerate()
            {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                *output = [finite(speed * unit[0])?, finite(speed * unit[1])?];
            }
            let field = self.field.step_with_eulerian_velocity_candidate(
                target_s,
                &velocity,
                impressed_jz_a_m2,
                field_controls,
                cancelled,
            )?;
            let loads = field.lorentz_loads(
                &self.binding.triangles,
                self.binding.load_limits,
                cancelled,
            )?;
            let world = loads.wrench(&Se3::identity(), [origin.x, origin.y, origin.z])?;
            let midpoint_angle = finite(self.angle_rad + 0.5 * duration * speed)?;
            let kinematics =
                forward_kinematics(&self.binding.model, base(), &[midpoint_angle], &[speed])?;
            let body = loads.wrench(&kinematics.world_from_link[0].inverse()?, [0.0; 3])?;
            let response = apply_prescribed_base_impulses(
                &self.binding.model,
                base(),
                &[midpoint_angle],
                &[old_speed],
                &[finite(duration * shaft_torque_nm)?],
                &[Wrench::new(
                    body.torque.scale(duration),
                    body.force.scale(duration),
                )],
            )?;
            poll(cancelled)?;
            let endpoint_speed = response.joint_velocity_after[0];
            let actual_midpoint = f64::midpoint(old_speed, endpoint_speed);
            last_residual = finite(actual_midpoint - speed)?.abs();
            let shaft_work_j = finite(duration * shaft_torque_nm * actual_midpoint)?;
            let load_work_residual_j = finite(
                (Dd::from_f64(field.step().mechanical_work_j)
                    - Dd::from_f64(duration) * Dd::from_f64(speed) * Dd::from_f64(world.torque.z))
                .to_f64(),
            )?;
            let mechanical_work_residual_j = finite(
                (Dd::from_f64(response.kinetic_energy_change)
                    - Dd::from_f64(field.step().mechanical_work_j)
                    - Dd::from_f64(shaft_work_j))
                .to_f64(),
            )?;
            let total_energy_residual_j = finite(
                (Dd::from_f64(field.step().energy_residual_j)
                    + Dd::from_f64(mechanical_work_residual_j))
                .to_f64(),
            )?;
            if load_work_residual_j.abs() > controls.energy_tolerance_j {
                return Err(RotorError::Balance(
                    "native Lorentz torque and native field work disagree",
                ));
            }
            if last_residual <= controls.angular_velocity_tolerance_rad_s
                && mechanical_work_residual_j.abs() <= controls.energy_tolerance_j
                && total_energy_residual_j.abs() <= controls.energy_tolerance_j
            {
                if response.angular_momentum_residual_world.norm()
                    > controls.angular_impulse_tolerance_nms
                    || response.linear_momentum_residual_world.norm()
                        > controls.linear_impulse_tolerance_ns
                    || response.kinetic_work_residual.abs() > controls.energy_tolerance_j
                {
                    return Err(RotorError::Balance(
                        "native mechanical impulse/work residual exceeds coupling budget",
                    ));
                }
                let angle = finite(self.angle_rad + duration * actual_midpoint)?;
                if duration * actual_midpoint != 0.0 && angle == self.angle_rad {
                    return Err(RotorError::Invalid(
                        "rotor angle increment disappeared in endpoint",
                    ));
                }
                forward_kinematics(&self.binding.model, base(), &[angle], &[endpoint_speed])?;
                let balance = RotorBalance {
                    iterations: iteration,
                    field_angular_velocity_rad_s: speed,
                    mechanical_midpoint_velocity_rad_s: actual_midpoint,
                    velocity_residual_rad_s: last_residual,
                    electromagnetic_torque_nm: world.torque.z,
                    shaft_work_j,
                    load_work_residual_j,
                    mechanical_work_residual_j,
                    total_energy_residual_j,
                };
                let (field, induction) = field.into_parts();
                let state = Self {
                    field,
                    binding: Arc::clone(&self.binding),
                    angle_rad: angle,
                    angular_velocity_rad_s: endpoint_speed,
                };
                poll(cancelled)?;
                return Ok(RotorCandidate {
                    state,
                    induction,
                    mechanics: response,
                    balance,
                });
            }
            let next = finite(speed + controls.relaxation * (actual_midpoint - speed))?;
            if next == speed {
                return Err(RotorError::Balance(
                    "coupled midpoint stalled before work/velocity acceptance",
                ));
            }
            speed = next;
        }
        Err(RotorError::Nonconvergence {
            iterations: controls.maximum_iterations,
            velocity_residual_rad_s: last_residual,
        })
    }
}

// One or two complete, consistently oriented circular boundary cycles. Positive
// cell orientation, opposed interior edges and one connected cell component
// bind the declared analytic section to the actual original conductor cells.
fn admit_circle(
    field: &PlanarEddy<'_>,
    triangles: &[usize],
    origin: Vec3,
    section: CircularConductor,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let CircularConductor {
        inner_radius_m: inner,
        outer_radius_m: outer,
        maximum_boundary_error_m: tolerance,
    } = section;
    if [inner, outer, tolerance]
        .into_iter()
        .any(|x| !x.is_finite())
        || inner < 0.0
        || outer <= inner
        || tolerance <= 0.0
        || tolerance >= 0.125 * (outer - inner)
        || tolerance >= 0.125 * outer
        || inner > 0.0 && tolerance >= 0.125 * inner
    {
        return Err(RotorError::Invalid(
            "bounded positive circular-section approximation required",
        ));
    }
    let mesh = field.problem().mesh();
    let relative = |v: u32| {
        let p = mesh.vertices()[v as usize];
        [p[0] - origin.x, p[1] - origin.y]
    };
    let mut edges = BTreeMap::<(u32, u32), Vec<(usize, u32, u32)>>::new();
    for (cell, &triangle) in triangles.iter().enumerate() {
        poll(cancelled)?;
        let vertices = mesh.faces()[triangle];
        let points = vertices.map(relative);
        let signed = finite(
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
                - (points[1][1] - points[0][1]) * (points[2][0] - points[0][0]),
        )?;
        if signed <= 0.0
            || points.iter().any(|p| {
                let r = p[0].hypot(p[1]);
                !r.is_finite() || r > outer + tolerance || r < inner - tolerance
            })
        {
            return Err(RotorError::Invalid(
                "conductor cells do not lie in the oriented circular section",
            ));
        }
        for i in 0..3 {
            let (a, b) = (vertices[i], vertices[(i + 1) % 3]);
            edges
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((cell, a, b));
        }
    }
    let mut adjacency = vec![Vec::new(); triangles.len()];
    let mut loops = [BTreeMap::<u32, u32>::new(), BTreeMap::<u32, u32>::new()];
    let mut incoming = [BTreeSet::new(), BTreeSet::new()];
    let mut angles = [Dd::ZERO; 2];
    let mut maximum_error: f64 = 0.0;
    for entries in edges.values() {
        poll(cancelled)?;
        match entries.as_slice() {
            &[(a, u, v), (b, x, y)] if u == y && v == x => {
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
            &[(_, u, v)] => {
                let a = relative(u);
                let b = relative(v);
                let ra = a[0].hypot(a[1]);
                let rb = b[0].hypot(b[1]);
                let which = if (ra - outer).abs().max((rb - outer).abs()) <= tolerance {
                    0
                } else if inner > 0.0 && (ra - inner).abs().max((rb - inner).abs()) <= tolerance {
                    1
                } else {
                    return Err(RotorError::Invalid(
                        "rotor boundary is not a declared circle",
                    ));
                };
                let radius = if which == 0 { outer } else { inner };
                let delta = [b[0] - a[0], b[1] - a[1]];
                let square = finite(delta[0] * delta[0] + delta[1] * delta[1])?;
                if square <= 0.0 {
                    return Err(RotorError::Invalid("zero circular boundary edge"));
                }
                let t = (-(a[0] * delta[0] + a[1] * delta[1]) / square).clamp(0.0, 1.0);
                let closest = (a[0] + t * delta[0]).hypot(a[1] + t * delta[1]);
                let error = finite(
                    (ra - radius)
                        .abs()
                        .max((rb - radius).abs())
                        .max((closest - radius).abs()),
                )?;
                let angle = finite(det::atan2(
                    a[0] * b[1] - a[1] * b[0],
                    a[0] * b[0] + a[1] * b[1],
                ))?;
                if error > tolerance
                    || angle == 0.0
                    || (angle > 0.0) != (which == 0)
                    || loops[which].insert(u, v).is_some()
                    || !incoming[which].insert(v)
                {
                    return Err(RotorError::Invalid(
                        "circular boundary chord, orientation or topology refused",
                    ));
                }
                angles[which] = angles[which] + Dd::from_f64(angle);
                maximum_error = maximum_error.max(error);
            }
            _ => {
                return Err(RotorError::Invalid(
                    "nonmanifold or multiply covered conductor edge",
                ));
            }
        }
    }
    for which in 0..2 {
        if which == 1 && inner == 0.0 {
            if !loops[which].is_empty() {
                return Err(RotorError::Invalid("disk has an undeclared hole"));
            }
            continue;
        }
        let Some(&start) = loops[which].keys().next() else {
            return Err(RotorError::Invalid("missing circular boundary"));
        };
        let mut vertex = start;
        for _ in 0..loops[which].len() {
            poll(cancelled)?;
            vertex = *loops[which]
                .get(&vertex)
                .ok_or(RotorError::Invalid("open circular boundary"))?;
        }
        let expected = if which == 0 {
            core::f64::consts::TAU
        } else {
            -core::f64::consts::TAU
        };
        if vertex != start || (angles[which].to_f64() - expected).abs() > 1e-10 {
            return Err(RotorError::Invalid(
                "circular boundary does not wind once about the rotor axis",
            ));
        }
    }
    let mut visited = vec![false; triangles.len()];
    let mut pending = vec![0];
    visited[0] = true;
    while let Some(cell) = pending.pop() {
        poll(cancelled)?;
        for &other in &adjacency[cell] {
            if !visited[other] {
                visited[other] = true;
                pending.push(other);
            }
        }
    }
    if visited.iter().any(|v| !v) {
        return Err(RotorError::Invalid(
            "conductor is not one connected circular section",
        ));
    }
    Ok(maximum_error)
}
