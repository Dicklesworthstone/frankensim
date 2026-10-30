//! Joint Hertz sliding/spin friction on two original finite rigid bodies.
//!
//! This is a held-geometry, held-effort midpoint impulse phase. Actual masses,
//! full inertia tensors, orientations and COM-relative arms determine the full
//! three-coordinate recoil. Normal load, patch radius, contact discovery and
//! pose advance remain with their original owners. Static midpoint selection
//! is not endpoint arrest: a coarse phase that would assist slip or spin at an
//! endpoint refuses and requires the trajectory owner to refine its interval.

use super::{RigidContactBinding, RigidFrictionError, components, magnitude, poll, totals, vector};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::TriboError;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchStep,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

/// One caller-owned closed normal patch between two finite dynamic bodies.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzContact<'a> {
    /// Original body indices and declared-body-frame COM-relative contact arms.
    pub binding: RigidContactBinding,
    /// Original normal and tangent gauge, held over the impulse phase.
    pub frame: TangentFrame,
    /// Original admitted normal load, circular radius and local coefficients.
    pub law: &'a HertzWrenchContact,
}

/// Independent native mechanics gates in addition to the original law budgets.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzControls {
    pub maximum_bodies: usize,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
    pub contact: HertzWrenchConfig,
}

/// Actual two-body observables and the unchanged native constitutive result.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzReport {
    pub contact: HertzWrenchStep,
    /// Full `J M^-1 J^T` for `[slip_x, slip_y, a*relative_spin]` [1/kg].
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    /// Original held-effort midpoint response `dt² W / 2` [m/N].
    pub scaled_compliance_m_per_n: [[f64; 3]; 3],
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    pub spin_before_rad_s: f64,
    pub spin_after_rad_s: f64,
    /// Exposed for the surrounding normal owner; this phase supplies no normal impulse.
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub actual_midpoint_displacement_m: [f64; 2],
    pub actual_midpoint_spin_angle_rad: f64,
    /// Contact impulse acting on A; B receives its exact opposite.
    pub force_impulse_world_ns: Vec3,
    /// Free normal-axis couple acting on A; B receives its exact opposite.
    pub spin_angular_impulse_world_nms: Vec3,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Actual point-force plus free-couple work on both finite bodies [J].
    pub midpoint_impulse_work_j: f64,
    /// Original Hertz local dissipation evaluated at actual native midpoint rates.
    pub dissipated_energy_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    pub contact_heat_residual_j: f64,
    pub native_force_residual_n: f64,
    pub native_spin_moment_residual_nm: f64,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_spin_residual_rad_s: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    /// Positive endpoint shear/spin power times dt; an unresolved reversal
    /// cannot be hidden by a passive midpoint or labeled endpoint stiction.
    pub maximum_assisting_work_j: f64,
}

/// Complete private proposal in the original binding's A/B order.
/// Publish both native endpoints and any coupled thermal/source state together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub report: RigidHertzReport,
}

/// Original finite pair and its full shear/spin response prepared once.
///
/// Preparation snapshots only the two bound bodies, their actual mass properties,
/// contact geometry and one fixed duration. Repeated trials change neither those
/// snapshots nor their native response. Normal force and radius remain exactly
/// fixed; a source owner may supply temperature-resolved static/kinetic friction
/// coefficients through the same original Hertz law. Source identity and thermal
/// applicability belong to that source owner, never this numerical adapter.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzPair {
    binding: RigidContactBinding,
    frame: TangentFrame,
    before: [RigidBodyState; 2],
    pair_properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    points_before: [PointKinematics; 2],
    normal_force_n: f64,
    radius: f64,
    dt_s: f64,
    controls: RigidHertzControls,
    mobility: [[f64; 3]; 3],
    input: HertzWrenchInput,
    rate_before: [f64; 3],
    totals_before: (f64, Vec3, Vec3),
}

/// Every refusal leaves the caller's body slice byte-for-byte unchanged.
#[derive(Debug)]
pub enum RigidHertzError {
    Invalid(&'static str),
    Native(HertzWrenchError),
    Law(TriboError),
    Rigid(RigidFrictionError),
    Dynamics(DynamicsError),
    Cancelled,
    /// The held wrench would assist endpoint slip/spin. Refine the physical
    /// interval; no kinetic override, endpoint snap, or fabricated heat is used.
    UnresolvedSlipEvent {
        assisting_work_j: f64,
    },
    Audit(Box<RigidHertzReport>),
}
impl From<DynamicsError> for RigidHertzError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<TriboError> for RigidHertzError {
    fn from(error: TriboError) -> Self {
        Self::Law(error)
    }
}
impl From<RigidFrictionError> for RigidHertzError {
    fn from(error: RigidFrictionError) -> Self {
        match error {
            RigidFrictionError::Cancelled => Self::Cancelled,
            _ => Self::Rigid(error),
        }
    }
}
impl From<HertzWrenchError> for RigidHertzError {
    fn from(error: HertzWrenchError) -> Self {
        match error {
            HertzWrenchError::Cancelled => Self::Cancelled,
            _ => Self::Native(error),
        }
    }
}
impl core::fmt::Display for RigidHertzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid Hertz friction: {self:?}")
    }
}
impl std::error::Error for RigidHertzError {}

/// Apply one joint force/spin candidate to both original bodies atomically.
///
/// The native Hertz solve uses `x_free = dt * relative_rate` and `C = dt² W/2`.
/// Native endpoint point kinematics independently determine actual midpoint
/// work and local friction heat. Only the original two bodies are staged;
/// unaffected bodies are neither copied nor altered. All original poses stay
/// fixed. There is no external drive, normal impulse, heat store or elastic
/// history reset. Exact kinetic stopping is retained when the held phase can
/// resolve it; a midpoint-static reversal refuses rather than inventing arrest.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_hertz(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    gate: &CancelGate,
) -> Result<RigidHertzReport, RigidHertzError> {
    let prepared = PreparedRigidHertzPair::new(bodies, properties, contact, dt_s, controls, gate)?;
    let candidate = prepared.advance(contact.law, gate)?;
    poll(gate)?;
    // The only publication boundary; no fallible operation follows it.
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

impl PreparedRigidHertzPair {
    /// Admit and snapshot one unchanged physical pair and interval. No body or
    /// constitutive solve is published; both endpoints remain caller-owned.
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzContact<'_>,
        dt_s: f64,
        controls: RigidHertzControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzError> {
        poll(gate)?;
        controls.contact.validate()?;
        let binding = contact.binding;
        if bodies.len() != properties.len()
            || bodies.len() > controls.maximum_bodies
            || binding.body_a >= bodies.len()
            || binding.body_b >= bodies.len()
            || binding.body_a == binding.body_b
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || [
                controls.geometry_tolerance_m,
                controls.velocity_tolerance_m_s,
                controls.angular_velocity_tolerance_rad_s,
                controls.linear_momentum_tolerance_ns,
                controls.angular_momentum_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(RigidHertzError::Invalid(
                "matching bounded native bodies, distinct pair, positive duration and finite gates required",
            ));
        }
        let before = [bodies[binding.body_a], bodies[binding.body_b]];
        let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
        let arms = [binding.arm_a_body, binding.arm_b_body];
        let axes = [
            vector(contact.frame.longitudinal()),
            vector(contact.frame.lateral()),
            vector(contact.frame.normal()),
        ];
        let points_before = points(before, pair_properties, arms)?;
        if magnitude(
            points_before[0]
                .point_world
                .sub(points_before[1].point_world),
        ) > controls.geometry_tolerance_m
        {
            return Err(RigidHertzError::Invalid(
                "original contact points do not coincide",
            ));
        }
        let radius = contact.law.parameters().contact_radius_m;
        let mobility = response(before, pair_properties, arms, axes, radius)?;
        let half_dt_squared = finite(0.5 * dt_s * dt_s)?;
        if half_dt_squared == 0.0 {
            return Err(RigidHertzError::Invalid(
                "unrepresentable midpoint duration",
            ));
        }
        let compliance = mobility.map(|row| row.map(|v| v * half_dt_squared));
        let rate_before = rates(points_before, axes, radius);
        let input = HertzWrenchInput {
            free_displacement_m: [dt_s * rate_before[0], dt_s * rate_before[1]],
            free_spin_angle_rad: dt_s * relative_spin(points_before).dot(axes[2]),
            scaled_compliance_m_per_n: compliance,
        };
        let totals_before = totals(&before, &pair_properties, gate)?;
        poll(gate)?;
        Ok(Self {
            binding,
            frame: contact.frame,
            before,
            pair_properties,
            arms,
            axes,
            points_before,
            normal_force_n: contact.law.parameters().normal_force_n,
            radius,
            dt_s,
            controls,
            mobility,
            input,
            rate_before,
            totals_before,
        })
    }

    /// Original caller indices. Candidate endpoints always follow this A/B order.
    #[must_use]
    pub const fn binding(&self) -> RigidContactBinding {
        self.binding
    }
    #[must_use]
    pub const fn frame(&self) -> TangentFrame {
        self.frame
    }
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }
    #[must_use]
    pub const fn normal_force_n(&self) -> f64 {
        self.normal_force_n
    }
    #[must_use]
    pub const fn contact_radius_m(&self) -> f64 {
        self.radius
    }
    #[must_use]
    pub const fn controls(&self) -> RigidHertzControls {
        self.controls
    }
    /// Original native world points, independently available for thermal geometry.
    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.points_before.map(|point| point.point_world)
    }
    #[must_use]
    pub const fn scaled_mobility_per_kg(&self) -> &[[f64; 3]; 3] {
        &self.mobility
    }
    #[must_use]
    pub const fn scaled_compliance_m_per_n(&self) -> &[[f64; 3]; 3] {
        &self.input.scaled_compliance_m_per_n
    }

    /// Propose original native recoil with temperature-resolved coefficients.
    ///
    /// The source must preserve the admitted normal force and circular radius
    /// exactly. It may vary only the original static/kinetic coefficients. Each
    /// trial starts from the same captured momenta; no previous proposal is
    /// accumulated. Original endpoint work, physical Hertz heat, momentum and
    /// coarse-reversal gates all run before returning the two owned endpoints.
    /// Any refusal or cancellation leaves preparation reusable for exact retry.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        law: &HertzWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidHertzCandidate, RigidHertzError> {
        poll(gate)?;
        let parameters = law.parameters();
        if parameters.normal_force_n.to_bits() != self.normal_force_n.to_bits()
            || parameters.contact_radius_m.to_bits() != self.radius.to_bits()
        {
            return Err(RigidHertzError::Invalid(
                "prepared Hertz normal force and radius must remain unchanged",
            ));
        }
        let Self {
            before,
            pair_properties,
            arms,
            axes,
            points_before,
            radius,
            dt_s,
            controls,
            mobility,
            input,
            rate_before,
            ..
        } = *self;
        let compliance = input.scaled_compliance_m_per_n;
        let resolved = law.advance(&input, controls.contact, || gate.is_requested())?;
        let force = axes[0]
            .scale(resolved.force_n[0])
            .add(axes[1].scale(resolved.force_n[1]));
        let couple = axes[2].scale(resolved.spin_moment_nm);
        let impulse = force.scale(dt_s);
        let angular_impulse = couple.scale(dt_s);
        let mut after = before;
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let forced = before[i]
                .apply_impulse_at_body_point(pair_properties[i], arms[i], impulse.scale(sign))?
                .state_after;
            let body_couple = forced
                .pose()
                .orientation()
                .rotate_world_to_body(angular_impulse.scale(sign));
            after[i] = RigidBodyState::new(
                forced.pose(),
                forced.linear_momentum_world(),
                forced.angular_momentum_body().add(body_couple),
            )?;
        }
        let points_after = points(after, pair_properties, arms)?;
        let relative_before = relative_point(points_before);
        let relative_after = relative_point(points_after);
        let spin_before = relative_spin(points_before);
        let spin_after = relative_spin(points_after);
        let midpoint = relative_before.add(relative_after).scale(0.5);
        let midpoint_spin = spin_before.add(spin_after).scale(0.5).dot(axes[2]);
        let actual_displacement = [dt_s * midpoint.dot(axes[0]), dt_s * midpoint.dot(axes[1])];
        let actual_spin_angle = dt_s * midpoint_spin;
        let midpoint_work =
            impulse.dot(midpoint) + angular_impulse.dot(spin_before.add(spin_after).scale(0.5));
        let mut heat = 0.0;
        let mut force_residual = 0.0;
        let mut moment_residual = 0.0;
        if resolved.mode != HertzWrenchMode::Sticking {
            let parameters = law.parameters();
            let native = HertzSlidingPatch::new(HertzSlidingPatchParameters {
                normal_force_n: parameters.normal_force_n,
                contact_radius_m: parameters.contact_radius_m,
                friction_coefficient: parameters.kinetic_mu,
            })?
            .evaluate(
                [midpoint.dot(axes[0]), midpoint.dot(axes[1])],
                midpoint_spin,
            )?;
            heat = finite(dt_s * native.dissipated_power_w)?;
            force_residual = (native.force_n[0] - resolved.force_n[0])
                .hypot(native.force_n[1] - resolved.force_n[1]);
            moment_residual = (native.spin_moment_nm - resolved.spin_moment_nm).abs();
        }
        let rate_after = rates(points_after, axes, radius);
        let expected_after = [
            2.0 * resolved.slip_displacement_m[0] / dt_s - rate_before[0],
            2.0 * resolved.slip_displacement_m[1] / dt_s - rate_before[1],
            2.0 * radius * resolved.spin_angle_rad / dt_s - rate_before[2],
        ];
        let endpoint_velocity_residual = (rate_after[0] - expected_after[0])
            .abs()
            .max((rate_after[1] - expected_after[1]).abs());
        let endpoint_spin_residual = (rate_after[2] - expected_after[2]).abs() / radius;
        let (energy_before, linear_before, angular_before) = self.totals_before;
        let (energy_after, linear_after, angular_after) = totals(&after, &pair_properties, gate)?;
        let kinetic_change = energy_after - energy_before;
        // A passive midpoint alone permits a spurious coarse-step reversal.
        // Check the individual shear/spin ports at BOTH actual endpoints. Their
        // held-effort powers vary linearly in this frozen-configuration phase.
        let assisting_work = [
            force.dot(relative_before),
            force.dot(relative_after),
            couple.dot(spin_before),
            couple.dot(spin_after),
        ]
        .into_iter()
        .fold(0.0, f64::max)
            * dt_s;
        let report = RigidHertzReport {
            contact: resolved,
            scaled_mobility_per_kg: mobility,
            scaled_compliance_m_per_n: compliance,
            slip_velocity_before_m_s: [rate_before[0], rate_before[1]],
            slip_velocity_after_m_s: [rate_after[0], rate_after[1]],
            spin_before_rad_s: spin_before.dot(axes[2]),
            spin_after_rad_s: spin_after.dot(axes[2]),
            normal_velocity_before_m_s: relative_before.dot(axes[2]),
            normal_velocity_after_m_s: relative_after.dot(axes[2]),
            actual_midpoint_displacement_m: actual_displacement,
            actual_midpoint_spin_angle_rad: actual_spin_angle,
            force_impulse_world_ns: impulse,
            spin_angular_impulse_world_nms: angular_impulse,
            kinetic_before_j: energy_before,
            kinetic_after_j: energy_after,
            midpoint_impulse_work_j: midpoint_work,
            dissipated_energy_j: heat,
            work_energy_residual_j: kinetic_change - midpoint_work,
            total_energy_residual_j: kinetic_change + heat,
            contact_work_residual_j: midpoint_work + resolved.contact_work_j,
            contact_heat_residual_j: heat - resolved.dissipated_energy_j,
            native_force_residual_n: force_residual,
            native_spin_moment_residual_nm: moment_residual,
            endpoint_velocity_residual_m_s: endpoint_velocity_residual,
            endpoint_spin_residual_rad_s: endpoint_spin_residual,
            linear_momentum_change_ns: linear_after.sub(linear_before),
            angular_momentum_change_nms: angular_after.sub(angular_before),
            maximum_assisting_work_j: assisting_work,
        };
        let energy_residuals = [
            report.work_energy_residual_j,
            report.total_energy_residual_j,
            report.contact_work_residual_j,
            report.contact_heat_residual_j,
        ];
        if [
            energy_before,
            energy_after,
            midpoint_work,
            heat,
            assisting_work,
            force_residual,
            moment_residual,
            endpoint_velocity_residual,
            endpoint_spin_residual,
            magnitude(report.linear_momentum_change_ns),
            magnitude(report.angular_momentum_change_nms),
        ]
        .into_iter()
        .any(|v| !v.is_finite())
            || energy_residuals
                .into_iter()
                .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
            || force_residual > controls.contact.force_tolerance_n
            || moment_residual > controls.contact.moment_tolerance_nm
            || endpoint_velocity_residual > controls.velocity_tolerance_m_s
            || endpoint_spin_residual > controls.angular_velocity_tolerance_rad_s
            || magnitude(report.linear_momentum_change_ns) > controls.linear_momentum_tolerance_ns
            || magnitude(report.angular_momentum_change_nms)
                > controls.angular_momentum_tolerance_nms
        {
            return Err(RigidHertzError::Audit(Box::new(report)));
        }
        if assisting_work > controls.energy_tolerance_j {
            return Err(RigidHertzError::UnresolvedSlipEvent {
                assisting_work_j: assisting_work,
            });
        }
        poll(gate)?;
        Ok(RigidHertzCandidate {
            bodies_after: after,
            report,
        })
    }
}

fn points(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
) -> Result<[PointKinematics; 2], DynamicsError> {
    Ok([
        bodies[0].point_kinematics(properties[0], arms[0])?,
        bodies[1].point_kinematics(properties[1], arms[1])?,
    ])
}
fn relative_point(points: [PointKinematics; 2]) -> Vec3 {
    points[0]
        .point_velocity_world
        .sub(points[1].point_velocity_world)
}
fn relative_spin(points: [PointKinematics; 2]) -> Vec3 {
    points[0]
        .angular_velocity_world
        .sub(points[1].angular_velocity_world)
}
fn rates(points: [PointKinematics; 2], axes: [Vec3; 3], radius: f64) -> [f64; 3] {
    let relative = relative_point(points);
    [
        relative.dot(axes[0]),
        relative.dot(axes[1]),
        radius * relative_spin(points).dot(axes[2]),
    ]
}
fn finite(value: f64) -> Result<f64, RigidHertzError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidHertzError::Invalid(
            "unrepresentable native Hertz mechanics",
        ))
    }
}

/// Same native Gram construction as the original finite-patch rigid owner,
/// specialized to one pair and its actual Hertz radius. All nine entries stay.
fn response(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    radius: f64,
) -> Result<[[f64; 3]; 3], RigidHertzError> {
    let mut factor = [[0.0; 12]; 3];
    for (body, sign) in [(0, 1.0), (1, -1.0)] {
        let inverse_mass_root = properties[body].mass().sqrt().recip();
        let orientation = bodies[body].pose().orientation();
        for axis in 0..3 {
            let direction = orientation.rotate_world_to_body(axes[axis]);
            let translation = if axis < 2 {
                components(axes[axis].scale(sign * inverse_mass_root))
            } else {
                [0.0; 3]
            };
            let rotation = if axis < 2 {
                arms[body].cross(direction).scale(sign)
            } else {
                direction.scale(sign * radius)
            };
            let rotation = components(properties[body].angular_mobility_factor(rotation));
            for j in 0..3 {
                factor[axis][6 * body + j] = translation[j];
                factor[axis][6 * body + 3 + j] = rotation[j];
            }
        }
    }
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            result[i][j] = finite(factor[i].iter().zip(&factor[j]).map(|(a, b)| a * b).sum())?;
            result[j][i] = result[i][j];
        }
    }
    Ok(result)
}
