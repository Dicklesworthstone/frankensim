//! Simultaneous Hertz shear/spin and viscoelastic sphere rolling on one pair.
//!
//! This is joint dynamics of the two original named laws under a declared
//! normal/shear decoupling approximation: the original symmetric Hertz surface
//! traction and the original first-order normal-pressure bulk rolling moment.
//! Traction/normal-field feedback is omitted. In particular, the O(A) pressure
//! asymmetry causing rolling can also perturb the friction integrals; this is
//! not a complete first-order mixed-contact continuum solution. Work/passivity
//! checks apply to this stated model, with all original applicability gates.
//!
//! Five conjugate coordinates retain every finite-body cross response:
//! `x = [du_x, du_y, a*dtheta_spin, R*dtheta_roll_x, R*dtheta_roll_y]` and
//! `w = [Fx, Fy, Mspin/a, Mroll_x/R, Mroll_y/R]`. The material rolling law
//! supplies `c=A*N*R`; exact two-coordinate Schur elimination of its viscous
//! moment precedes the original joint Hertz solve. One combined native impulse
//! supplies all recoil. Surface heat and sphere bulk heat remain distinct.
//!
//! Geometry, normal load, sphere material temperature and frequency are held.
//! Poses and the normal contact remain caller-owned. A held effort that assists
//! an endpoint shear, spin or rolling port refuses, requiring time refinement.
//! No endpoint snapping, artificial compliance, partial-slip memory or thermal
//! store is added. The association of the supplied Hertz law with the ordered
//! interface remains caller-owned unless a source-bound wrapper establishes it.

use super::hertz::{
    PreparedRigidHertzPair, RigidHertzContact, RigidHertzControls, RigidHertzError,
    apply_pair_wrench, points, relative_point, relative_spin, scaled_response,
};
use super::rolling::RigidRollingHistory;
use super::{RigidContactBinding, RigidFrictionError, magnitude, poll, totals};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::{DeformingSphereSurface, ViscoelasticSphereRollingCard};
use fs_tribo::rolling_loss::{
    RollingKinematics, RollingLossError, RollingLossLaw, RollingLossState, RollingLossStep,
    RollingPatchReceipt, RollingWorkOwnership,
};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchContact, HertzWrenchInput, HertzWrenchMode, HertzWrenchStep,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};
use fs_tribo::{InterfaceSystemRef, TriboError};

/// One actual sphere/rigid-plane patch. The frame normal points from B to A.
/// Hertz N/area and the original rolling receipt must describe that same patch.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub hertz: &'a HertzWrenchContact,
    pub patch: &'a RollingPatchReceipt,
    pub interface: &'a InterfaceSystemRef,
    /// The original material-derived `ViscoelasticSphere` variant only.
    pub rolling: &'a RollingLossLaw,
    pub temperature_k: f64,
    /// Original material frequency coordinate; never inferred from dt or spin.
    pub excitation_frequency_hz: f64,
    /// One dedicated tangent-vector rolling channel, with its native checkpoint.
    pub work: &'a RollingWorkOwnership,
}

/// Original Hertz budgets/native mechanics gates plus the rolling effort gate.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingControls {
    pub mechanical: RigidHertzControls,
    pub maximum_rolling_moment_nm: f64,
    pub rolling_moment_tolerance_nm: f64,
}

/// One accepted candidate's original constitutive responses and actual work.
/// Energy/momentum totals cover only the pair; every other body is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingReport {
    pub hertz: HertzWrenchStep,
    pub rolling: RollingLossStep,
    /// Original body-slice index of the sphere generating all bulk heat.
    pub deforming_body: usize,
    /// Actual native material-card coordinates used for this accepted trial.
    pub rolling_temperature_k: f64,
    pub rolling_excitation_frequency_hz: f64,
    pub scaled_mobility_per_kg: [[f64; 5]; 5],
    pub scaled_compliance_m_per_n: [[f64; 5]; 5],
    /// Eliminated numerical response; its free motion is not external work.
    pub effective_hertz_compliance_m_per_n: [[f64; 3]; 3],
    pub scaled_rates_before_m_s: [f64; 5],
    pub scaled_rates_after_m_s: [f64; 5],
    pub actual_midpoint_coordinates_m: [f64; 5],
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub force_impulse_world_ns: Vec3,
    pub spin_angular_impulse_world_nms: Vec3,
    pub rolling_angular_impulse_world_nms: Vec3,
    pub surface_heat_j: f64,
    pub bulk_heat_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub surface_midpoint_work_j: f64,
    pub rolling_midpoint_work_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub surface_work_residual_j: f64,
    pub surface_heat_residual_j: f64,
    pub rolling_work_residual_j: f64,
    pub history_heat_residual_j: f64,
    pub native_force_residual_n: f64,
    pub native_spin_moment_residual_nm: f64,
    pub native_rolling_moment_residual_nm: f64,
    pub maximum_compatibility_residual_m: f64,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_angular_velocity_residual_rad_s: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub maximum_assisting_work_j: f64,
}

/// All endpoints and the original rolling checkpoint are proposed together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub rolling_history: RigidRollingHistory,
    pub report: RigidHertzRollingReport,
}

/// Immutable original pair, fixed geometry/interval and full five-axis response.
/// Repeated coefficient trials start from the same accepted momenta and history.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzRollingPair<'a> {
    core: PreparedRigidHertzPair,
    contact: RigidHertzRollingContact<'a>,
    history: RigidRollingHistory,
    controls: RigidHertzRollingControls,
    radius_m: f64,
    deforming_body: usize,
    mobility: [[f64; 5]; 5],
    compliance: [[f64; 5]; 5],
    rate_before: [f64; 5],
    rolling_factor: Spd2,
    effective: HertzWrenchInput,
}

/// Every error leaves caller-owned body states and rolling history unchanged.
#[derive(Debug)]
pub enum RigidHertzRollingError {
    Invalid(&'static str),
    Hertz(RigidHertzError),
    Rolling(RollingLossError),
    /// A held component would assist an actual endpoint. Refine the interval.
    UnresolvedSlipEvent {
        assisting_work_j: f64,
    },
    RollingMomentBudget {
        required_nm: f64,
        maximum_nm: f64,
    },
    Audit(Box<RigidHertzRollingReport>),
}
impl From<RigidHertzError> for RigidHertzRollingError {
    fn from(e: RigidHertzError) -> Self {
        Self::Hertz(e)
    }
}
impl From<RigidFrictionError> for RigidHertzRollingError {
    fn from(e: RigidFrictionError) -> Self {
        Self::Hertz(e.into())
    }
}
impl From<DynamicsError> for RigidHertzRollingError {
    fn from(e: DynamicsError) -> Self {
        Self::Hertz(e.into())
    }
}
impl From<TriboError> for RigidHertzRollingError {
    fn from(e: TriboError) -> Self {
        Self::Hertz(e.into())
    }
}
impl From<RollingLossError> for RigidHertzRollingError {
    fn from(e: RollingLossError) -> Self {
        Self::Rolling(e)
    }
}
impl core::fmt::Display for RigidHertzRollingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid Hertz and sphere rolling: {self:?}")
    }
}
impl std::error::Error for RigidHertzRollingError {}

/// Apply a single joint wrench and original rolling history atomically.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_hertz_rolling(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut RigidRollingHistory,
    contact: RigidHertzRollingContact<'_>,
    dt_s: f64,
    controls: RigidHertzRollingControls,
    gate: &CancelGate,
) -> Result<RigidHertzRollingReport, RigidHertzRollingError> {
    let prepared = PreparedRigidHertzRollingPair::new(
        bodies, properties, history, contact, dt_s, controls, gate,
    )?;
    let candidate = prepared.advance(contact.hertz, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *history = candidate.rolling_history;
    Ok(candidate.report)
}

impl<'a> PreparedRigidHertzRollingPair<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        history: &RigidRollingHistory,
        contact: RigidHertzRollingContact<'a>,
        dt_s: f64,
        controls: RigidHertzRollingControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzRollingError> {
        poll(gate)?;
        if [
            controls.maximum_rolling_moment_nm,
            controls.rolling_moment_tolerance_nm,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(invalid(
                "finite nonnegative rolling effort and residual gates required",
            ));
        }
        let card = sphere_card(contact.rolling)?;
        admit_history(history, contact)?;
        let core = PreparedRigidHertzPair::new(
            bodies,
            properties,
            RigidHertzContact {
                binding: contact.binding,
                frame: contact.frame,
                law: contact.hertz,
            },
            dt_s,
            controls.mechanical,
            gate,
        )?;
        let parameters = contact.hertz.parameters();
        if parameters.normal_force_n != contact.patch.normal_load_n() {
            return Err(invalid(
                "Hertz and original rolling patch normal loads must agree",
            ));
        }
        // At release, the original normal patch has zero area; the positive
        // Hertz radius is then only an inactive work-coordinate scale.
        if parameters.normal_force_n > 0.0 {
            let area = finite(
                core::f64::consts::PI * parameters.contact_radius_m * parameters.contact_radius_m,
            )?;
            let original = contact.patch.contact_area_m2();
            if (area - original).abs() > 512.0 * f64::EPSILON * area.max(original) {
                return Err(invalid(
                    "Hertz radius must match the original rolling patch area",
                ));
            }
        }
        let radius_m = card.parameters().radius_m;
        let (slot, sign, deforming_body) = match card.parameters().sphere_surface {
            DeformingSphereSurface::SurfaceA => (0, -1.0, contact.binding.body_a),
            DeformingSphereSurface::SurfaceB => (1, 1.0, contact.binding.body_b),
        };
        if magnitude(
            core.points_before[slot]
                .arm_world
                .sub(core.axes[2].scale(sign * radius_m)),
        ) > controls.mechanical.geometry_tolerance_m
        {
            return Err(invalid(
                "original sphere arm must match its material radius and normal",
            ));
        }
        let rate_before = rates(
            core.points_before,
            core.axes,
            parameters.contact_radius_m,
            radius_m,
        );
        if parameters.normal_force_n > 0.0
            && relative_point(core.points_before).dot(core.axes[2]).abs()
                > controls.mechanical.velocity_tolerance_m_s
        {
            return Err(invalid(
                "loaded joint phase requires resolved normal contact",
            ));
        }
        // Check the original initial rate even when the midpoint will be slower.
        card.response(
            contact.patch,
            contact.interface,
            kinematics(
                contact,
                radius_m,
                norm([rate_before[3], rate_before[4]]) / radius_m,
                relative_spin(core.points_before).dot(core.axes[2]),
                dt_s,
            )?,
        )?;
        let axes = [
            core.axes[0],
            core.axes[1],
            core.axes[2],
            core.axes[0],
            core.axes[1],
        ];
        let mobility = scaled_response(
            core.before,
            core.pair_properties,
            core.arms,
            axes,
            [0.0, 0.0, parameters.contact_radius_m, radius_m, radius_m],
        )?;
        let half_dt_squared = finite(0.5 * dt_s * dt_s)?;
        let compliance = mobility.map(|row| row.map(|v| v * half_dt_squared));
        if compliance.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable joint compliance"));
        }
        let damping = finite((card.rolling_time_s() * radius_m) * parameters.normal_force_n)?;
        if damping == 0.0 && card.rolling_time_s() > 0.0 && parameters.normal_force_n > 0.0 {
            return Err(invalid("unrepresentable material rolling damping"));
        }
        let scaled_damping = finite((damping / radius_m) / radius_m)?;
        let half_damping = finite(0.5 * dt_s * scaled_damping)?;
        if damping > 0.0 && (scaled_damping == 0.0 || half_damping == 0.0) {
            return Err(invalid("unrepresentable scaled material rolling damping"));
        }
        // Algebraically h = (dt²/2)*d with d=c/(R²*dt). Working with h
        // avoids an artificial 1/dt overflow at a small physical timestep.
        let rolling_factor = Spd2::new([
            [
                1.0 + half_damping * mobility[3][3],
                half_damping * mobility[3][4],
            ],
            [
                half_damping * mobility[4][3],
                1.0 + half_damping * mobility[4][4],
            ],
        ])?;
        let effective = eliminate(&core, mobility, half_damping, rolling_factor, rate_before)?;
        poll(gate)?;
        Ok(Self {
            core,
            contact,
            history: history.clone(),
            controls,
            radius_m,
            deforming_body,
            mobility,
            compliance,
            rate_before,
            rolling_factor,
            effective,
        })
    }

    #[must_use]
    pub const fn binding(&self) -> RigidContactBinding {
        self.core.binding()
    }
    #[must_use]
    pub const fn frame(&self) -> TangentFrame {
        self.core.frame()
    }
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.core.duration_s()
    }
    #[must_use]
    pub const fn normal_force_n(&self) -> f64 {
        self.core.normal_force_n()
    }
    #[must_use]
    pub const fn contact_radius_m(&self) -> f64 {
        self.core.contact_radius_m()
    }
    #[must_use]
    pub const fn sphere_radius_m(&self) -> f64 {
        self.radius_m
    }
    #[must_use]
    pub const fn deforming_body(&self) -> usize {
        self.deforming_body
    }
    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.core.contact_points_world()
    }
    #[must_use]
    pub const fn controls(&self) -> RigidHertzRollingControls {
        self.controls
    }
    #[must_use]
    pub const fn scaled_mobility_per_kg(&self) -> &[[f64; 5]; 5] {
        &self.mobility
    }
    #[must_use]
    pub const fn scaled_compliance_m_per_n(&self) -> &[[f64; 5]; 5] {
        &self.compliance
    }

    /// Trial original Hertz coefficients with fixed load/radius and the same
    /// immutable rolling card, geometry, temperature, frequency and history.
    /// Every trial independently rechecks actual native work and applicability.
    pub fn advance(
        &self,
        hertz: &HertzWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidHertzRollingCandidate, RigidHertzRollingError> {
        self.advance_at_temperature(hertz, self.contact.temperature_k, gate)
    }

    /// Query the SAME constant-coefficient rolling card at a trial material
    /// temperature. Its complete original temperature/frequency/rate domain
    /// remains mandatory; no new material or viscosity is substituted.
    #[allow(clippy::too_many_lines)]
    pub fn advance_at_temperature(
        &self,
        hertz: &HertzWrenchContact,
        temperature_k: f64,
        gate: &CancelGate,
    ) -> Result<RigidHertzRollingCandidate, RigidHertzRollingError> {
        poll(gate)?;
        let contact = RigidHertzRollingContact {
            temperature_k,
            ..self.contact
        };
        let hp = hertz.parameters();
        if hp.normal_force_n.to_bits() != self.normal_force_n().to_bits()
            || hp.contact_radius_m.to_bits() != self.contact_radius_m().to_bits()
        {
            return Err(invalid(
                "prepared joint normal force and Hertz radius must remain unchanged",
            ));
        }
        let card = sphere_card(contact.rolling)?;
        let dt = self.duration_s();
        let radius = self.radius_m;
        let axes = self.core.axes;
        let controls = self.controls.mechanical;
        card.response(
            contact.patch,
            contact.interface,
            kinematics(
                contact,
                radius,
                norm([self.rate_before[3], self.rate_before[4]]) / radius,
                relative_spin(self.core.points_before).dot(axes[2]),
                dt,
            )?,
        )?;
        let resolved = hertz
            .advance(&self.effective, controls.contact, || gate.is_requested())
            .map_err(RigidHertzError::from)?;
        let hertz_effort = [
            resolved.force_n[0],
            resolved.force_n[1],
            resolved.spin_moment_nm / hp.contact_radius_m,
        ];
        let roll_mid_scaled = self.rolling_factor.solve(std::array::from_fn(|i| {
            self.rate_before[i + 3]
                + 0.5
                    * dt
                    * (0..3)
                        .map(|j| self.mobility[i + 3][j] * hertz_effort[j])
                        .sum::<f64>()
        }))?;
        let predicted_roll = roll_mid_scaled.map(|v| v / radius);
        let response = card.response(
            contact.patch,
            contact.interface,
            kinematics(
                contact,
                radius,
                norm(predicted_roll),
                resolved.spin_angle_rad / dt,
                dt,
            )?,
        )?;
        let rolling_moment = along(predicted_roll, response.rolling_moment_nm);
        let rolling_norm = norm(rolling_moment);
        if rolling_norm > self.controls.maximum_rolling_moment_nm {
            return Err(RigidHertzRollingError::RollingMomentBudget {
                required_nm: rolling_norm,
                maximum_nm: self.controls.maximum_rolling_moment_nm,
            });
        }
        let force = world(resolved.force_n, axes);
        let spin_couple = axes[2].scale(resolved.spin_moment_nm);
        let rolling_couple = world(rolling_moment, axes);
        let impulse = force.scale(dt);
        let spin_impulse = spin_couple.scale(dt);
        let rolling_impulse = rolling_couple.scale(dt);
        let after = apply_pair_wrench(
            self.core.before,
            self.core.pair_properties,
            self.core.arms,
            impulse,
            spin_impulse.add(rolling_impulse),
            gate,
        )?;
        let after_points = points(after, self.core.pair_properties, self.core.arms)?;
        let point_before = relative_point(self.core.points_before);
        let point_after = relative_point(after_points);
        let omega_before = relative_spin(self.core.points_before);
        let omega_after = relative_spin(after_points);
        let midpoint = point_before.add(point_after).scale(0.5);
        let omega_mid = omega_before.add(omega_after).scale(0.5);
        let actual_roll = [omega_mid.dot(axes[0]), omega_mid.dot(axes[1])];
        // Exactly ONE original history advance for the whole isotropic rolling
        // vector. Native spin is retained but generates no rolling-channel heat.
        let rolling = contact.rolling.advance(
            contact.patch,
            contact.interface,
            kinematics(
                contact,
                radius,
                norm(actual_roll),
                omega_mid.dot(axes[2]),
                dt,
            )?,
            contact.work,
            self.history.state(),
        )?;
        card.response(
            contact.patch,
            contact.interface,
            kinematics(
                contact,
                radius,
                norm([omega_after.dot(axes[0]), omega_after.dot(axes[1])]),
                omega_after.dot(axes[2]),
                dt,
            )?,
        )?;
        if rolling.wrench.contour_force_n != 0.0
            || rolling.wrench.spin_moment_nm != 0.0
            || rolling.stored_energy_j != 0.0
            || rolling.storage_change_j != 0.0
        {
            return Err(invalid(
                "original sphere rolling must retain its sole zero-storage moment channel",
            ));
        }
        let native_rolling = along(actual_roll, rolling.wrench.rolling_moment_nm);
        let rolling_residual = norm([
            native_rolling[0] - rolling_moment[0],
            native_rolling[1] - rolling_moment[1],
        ]);
        let mut surface_heat = 0.0;
        let mut force_residual = 0.0;
        let mut spin_residual = 0.0;
        if resolved.mode != HertzWrenchMode::Sticking {
            let native = HertzSlidingPatch::new(HertzSlidingPatchParameters {
                normal_force_n: hp.normal_force_n,
                contact_radius_m: hp.contact_radius_m,
                friction_coefficient: hp.kinetic_mu,
            })?
            .evaluate(
                [midpoint.dot(axes[0]), midpoint.dot(axes[1])],
                omega_mid.dot(axes[2]),
            )?;
            surface_heat = finite(dt * native.dissipated_power_w)?;
            force_residual = norm([
                native.force_n[0] - resolved.force_n[0],
                native.force_n[1] - resolved.force_n[1],
            ]);
            spin_residual = (native.spin_moment_nm - resolved.spin_moment_nm).abs();
        }
        let rate_after = rates(after_points, axes, hp.contact_radius_m, radius);
        let actual: [f64; 5] =
            std::array::from_fn(|i| 0.5 * dt * (self.rate_before[i] + rate_after[i]));
        let expected = [
            resolved.slip_displacement_m[0],
            resolved.slip_displacement_m[1],
            hp.contact_radius_m * resolved.spin_angle_rad,
            dt * roll_mid_scaled[0],
            dt * roll_mid_scaled[1],
        ];
        let compatibility = (0..5)
            .map(|i| (actual[i] - expected[i]).abs())
            .fold(0.0, f64::max);
        let endpoint_velocity = (0..2)
            .map(|i| (rate_after[i] - (2.0 * expected[i] / dt - self.rate_before[i])).abs())
            .fold(0.0, f64::max);
        let endpoint_angular = (2..5)
            .map(|i| {
                (rate_after[i] - (2.0 * expected[i] / dt - self.rate_before[i])).abs()
                    / if i == 2 { hp.contact_radius_m } else { radius }
            })
            .fold(0.0, f64::max);
        let surface_work = impulse.dot(midpoint) + spin_impulse.dot(omega_mid);
        let rolling_work = rolling_impulse.dot(omega_mid);
        let work = surface_work + rolling_work;
        let bulk_heat = rolling.dissipation.total_heat_j;
        let history_heat =
            rolling.next_state.cumulative_heat_j() - self.history.state().cumulative_heat_j();
        let (energy_before, linear_before, angular_before) = self.core.totals_before;
        let (energy_after, linear_after, angular_after) =
            totals(&after, &self.core.pair_properties, gate)?;
        let assisting_work = [
            force.dot(point_before),
            force.dot(point_after),
            spin_couple.dot(omega_before),
            spin_couple.dot(omega_after),
            rolling_couple.dot(omega_before),
            rolling_couple.dot(omega_after),
        ]
        .into_iter()
        .fold(0.0, f64::max)
            * dt;
        let rolling_history = RigidRollingHistory::from_native_step(&rolling);
        let report = RigidHertzRollingReport {
            surface_work_residual_j: surface_work + resolved.contact_work_j,
            surface_heat_residual_j: surface_heat - resolved.dissipated_energy_j,
            hertz: resolved,
            rolling,
            deforming_body: self.deforming_body,
            rolling_temperature_k: temperature_k,
            rolling_excitation_frequency_hz: contact.excitation_frequency_hz,
            scaled_mobility_per_kg: self.mobility,
            scaled_compliance_m_per_n: self.compliance,
            effective_hertz_compliance_m_per_n: self.effective.scaled_compliance_m_per_n,
            scaled_rates_before_m_s: self.rate_before,
            scaled_rates_after_m_s: rate_after,
            actual_midpoint_coordinates_m: actual,
            normal_velocity_before_m_s: point_before.dot(axes[2]),
            normal_velocity_after_m_s: point_after.dot(axes[2]),
            force_impulse_world_ns: impulse,
            spin_angular_impulse_world_nms: spin_impulse,
            rolling_angular_impulse_world_nms: rolling_impulse,
            surface_heat_j: surface_heat,
            bulk_heat_j: bulk_heat,
            kinetic_before_j: energy_before,
            kinetic_after_j: energy_after,
            surface_midpoint_work_j: surface_work,
            rolling_midpoint_work_j: rolling_work,
            midpoint_impulse_work_j: work,
            work_energy_residual_j: energy_after - energy_before - work,
            total_energy_residual_j: energy_after - energy_before + surface_heat + bulk_heat,
            rolling_work_residual_j: rolling_work + bulk_heat,
            history_heat_residual_j: history_heat - bulk_heat,
            native_force_residual_n: force_residual,
            native_spin_moment_residual_nm: spin_residual,
            native_rolling_moment_residual_nm: rolling_residual,
            maximum_compatibility_residual_m: compatibility,
            endpoint_velocity_residual_m_s: endpoint_velocity,
            endpoint_angular_velocity_residual_rad_s: endpoint_angular,
            linear_momentum_change_ns: linear_after.sub(linear_before),
            angular_momentum_change_nms: angular_after.sub(angular_before),
            maximum_assisting_work_j: assisting_work,
        };
        audit(&report, self.controls)?;
        if hp.normal_force_n > 0.0
            && report.normal_velocity_after_m_s.abs() > controls.velocity_tolerance_m_s
        {
            return Err(invalid("joint recoil induced unresolved normal motion"));
        }
        if assisting_work > controls.energy_tolerance_j {
            return Err(RigidHertzRollingError::UnresolvedSlipEvent {
                assisting_work_j: assisting_work,
            });
        }
        poll(gate)?;
        Ok(RigidHertzRollingCandidate {
            bodies_after: after,
            rolling_history,
            report,
        })
    }
}

fn eliminate(
    core: &PreparedRigidHertzPair,
    w: [[f64; 5]; 5],
    h: f64,
    factor: Spd2,
    rate_before: [f64; 5],
) -> Result<HertzWrenchInput, RigidHertzRollingError> {
    if h == 0.0 {
        return Ok(core.input);
    }
    let dt = core.duration_s();
    let radius = core.contact_radius_m();
    let free_roll = factor.solve([dt * rate_before[3], dt * rate_before[4]])?;
    let free_hertz = [
        core.input.free_displacement_m[0],
        core.input.free_displacement_m[1],
        radius * core.input.free_spin_angle_rad,
    ];
    let mut released = [0.0; 3];
    let mut columns = [[0.0; 2]; 3];
    for i in 0..3 {
        released[i] =
            finite(free_hertz[i] - h * (w[i][3] * free_roll[0] + w[i][4] * free_roll[1]))?;
        columns[i] = factor.solve([w[3][i], w[4][i]])?;
    }
    let mut compliance = [[0.0; 3]; 3];
    let half_dt_squared = 0.5 * dt * dt;
    for i in 0..3 {
        for j in 0..=i {
            let response =
                finite(w[i][j] - h * (w[i][3] * columns[j][0] + w[i][4] * columns[j][1]))?;
            compliance[i][j] = finite(half_dt_squared * response)?;
            compliance[j][i] = compliance[i][j];
        }
    }
    Ok(HertzWrenchInput {
        free_displacement_m: [released[0], released[1]],
        free_spin_angle_rad: finite(released[2] / radius)?,
        scaled_compliance_m_per_n: compliance,
    })
}

#[derive(Clone, Copy, Debug)]
struct Spd2 {
    scale: f64,
    l00: f64,
    l10: f64,
    l11: f64,
}
impl Spd2 {
    fn new(a: [[f64; 2]; 2]) -> Result<Self, RigidHertzRollingError> {
        let scale = a.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale <= 0.0 || a.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable rolling Schur factor"));
        }
        let l00 = finite((a[0][0] / scale).sqrt())?;
        let l10 = finite((a[1][0] / scale) / l00)?;
        let pivot = finite(a[1][1] / scale - l10 * l10)?;
        if l00 <= 0.0 || pivot <= 0.0 {
            return Err(invalid("unresolved positive rolling Schur factor"));
        }
        Ok(Self {
            scale,
            l00,
            l10,
            l11: pivot.sqrt(),
        })
    }
    fn solve(self, rhs: [f64; 2]) -> Result<[f64; 2], RigidHertzRollingError> {
        let y0 = finite((rhs[0] / self.scale) / self.l00)?;
        let y1 = finite((rhs[1] / self.scale - self.l10 * y0) / self.l11)?;
        let x1 = finite(y1 / self.l11)?;
        let x0 = finite((y0 - self.l10 * x1) / self.l00)?;
        Ok([x0, x1])
    }
}

fn admit_history(
    history: &RigidRollingHistory,
    c: RigidHertzRollingContact<'_>,
) -> Result<(), RigidHertzRollingError> {
    if let Some(checkpoint) = history.checkpoint() {
        let restored = c
            .rolling
            .restore_checkpoint(c.patch, c.interface, c.work, checkpoint)?;
        if &restored != history.state() {
            return Err(invalid(
                "rolling history differs from its original checkpoint",
            ));
        }
    } else if history.state() != &RollingLossState::zero() {
        return Err(invalid(
            "nonzero rolling history needs its original checkpoint",
        ));
    }
    Ok(())
}
fn sphere_card(
    law: &RollingLossLaw,
) -> Result<&ViscoelasticSphereRollingCard, RigidHertzRollingError> {
    let RollingLossLaw::ViscoelasticSphere(card) = law else {
        return Err(invalid(
            "joint consumer requires the original viscoelastic sphere law",
        ));
    };
    Ok(card)
}
fn kinematics(
    c: RigidHertzRollingContact<'_>,
    radius: f64,
    rate: f64,
    spin: f64,
    dt: f64,
) -> Result<RollingKinematics, RollingLossError> {
    RollingKinematics::new(
        radius * rate,
        rate,
        spin,
        c.temperature_k,
        c.excitation_frequency_hz,
        dt,
    )
}
fn rates(p: [PointKinematics; 2], axes: [Vec3; 3], a: f64, r: f64) -> [f64; 5] {
    let v = relative_point(p);
    let omega = relative_spin(p);
    [
        v.dot(axes[0]),
        v.dot(axes[1]),
        a * omega.dot(axes[2]),
        r * omega.dot(axes[0]),
        r * omega.dot(axes[1]),
    ]
}
fn world(v: [f64; 2], axes: [Vec3; 3]) -> Vec3 {
    axes[0].scale(v[0]).add(axes[1].scale(v[1]))
}
fn norm(v: [f64; 2]) -> f64 {
    v[0].hypot(v[1])
}
fn along(v: [f64; 2], scalar: f64) -> [f64; 2] {
    let n = norm(v);
    if n == 0.0 {
        [0.0; 2]
    } else {
        v.map(|x| (x / n) * scalar)
    }
}
fn finite(v: f64) -> Result<f64, RigidHertzRollingError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("unrepresentable joint Hertz/rolling mechanics"))
    }
}
fn invalid(message: &'static str) -> RigidHertzRollingError {
    RigidHertzRollingError::Invalid(message)
}

fn audit(
    r: &RigidHertzRollingReport,
    c: RigidHertzRollingControls,
) -> Result<(), RigidHertzRollingError> {
    let m = c.mechanical;
    let energies = [
        r.work_energy_residual_j,
        r.total_energy_residual_j,
        r.surface_work_residual_j,
        r.surface_heat_residual_j,
        r.rolling_work_residual_j,
        r.history_heat_residual_j,
    ];
    if energies
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > m.energy_tolerance_j)
        || [
            r.surface_heat_j,
            r.bulk_heat_j,
            r.kinetic_before_j,
            r.kinetic_after_j,
            r.surface_midpoint_work_j,
            r.rolling_midpoint_work_j,
            r.midpoint_impulse_work_j,
            r.maximum_assisting_work_j,
            r.normal_velocity_before_m_s,
            r.normal_velocity_after_m_s,
            r.native_force_residual_n,
            r.native_spin_moment_residual_nm,
            r.native_rolling_moment_residual_nm,
            r.maximum_compatibility_residual_m,
            r.endpoint_velocity_residual_m_s,
            r.endpoint_angular_velocity_residual_rad_s,
            magnitude(r.linear_momentum_change_ns),
            magnitude(r.angular_momentum_change_nms),
        ]
        .into_iter()
        .any(|v| !v.is_finite())
        || r.surface_heat_j < 0.0
        || r.bulk_heat_j < 0.0
        || r.native_force_residual_n > m.contact.force_tolerance_n
        || r.native_spin_moment_residual_nm > m.contact.moment_tolerance_nm
        || r.native_rolling_moment_residual_nm > c.rolling_moment_tolerance_nm
        || r.maximum_compatibility_residual_m > m.contact.displacement_tolerance_m
        || r.endpoint_velocity_residual_m_s > m.velocity_tolerance_m_s
        || r.endpoint_angular_velocity_residual_rad_s > m.angular_velocity_tolerance_rad_s
        || magnitude(r.linear_momentum_change_ns) > m.linear_momentum_tolerance_ns
        || magnitude(r.angular_momentum_change_nms) > m.angular_momentum_tolerance_nms
    {
        return Err(RigidHertzRollingError::Audit(Box::new(r.clone())));
    }
    Ok(())
}
