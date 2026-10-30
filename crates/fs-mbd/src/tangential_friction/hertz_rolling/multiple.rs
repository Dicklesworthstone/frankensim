//! Simultaneous material rolling and Hertz shear/spin on shared finite bodies.
//!
//! Five coordinates per original patch retain all shared-body cross responses.
//! The original material rolling moments are eliminated together, then the
//! original coupled Hertz owner solves the remaining surface problem. Each
//! participating body receives one summed impulse. Each original rolling
//! history is advanced once, from its captured accepted state and the actual
//! shared-body endpoint. Surface and sphere bulk heat remain separate.
//!
//! This retains the parent's normal/shear decoupling approximation. Independent
//! normal patch fields and their non-overlapping physical ownership are supplied
//! by the caller. Intersecting spherical patch footprints and contradictory
//! material/temperature assignments to the same sphere are refused. Long-range
//! deformation-field interaction is not resolved. Loads, geometry and frequency
//! are held. No normal impulse, pose
//! update, partial-slip memory, regularization or artificial heat is introduced.

use super::{
    PreparedRigidHertzRollingPair, RigidHertzRollingContact, RigidHertzRollingControls,
    RigidHertzRollingError, along, finite, kinematics, norm, rates, sphere_card, world,
};
use crate::tangential_friction::hertz::multiple::RigidHertzIslandControls;
use crate::tangential_friction::hertz::{accepted_phase, points, relative_point, relative_spin};
use crate::tangential_friction::rolling::RigidRollingHistory;
use crate::tangential_friction::{
    RigidContactBinding, RigidFrictionError, components, magnitude, poll, totals,
};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::TriboError;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::RollingLossStep;
use fs_tribo::sliding_patch::contact::multiple::MultipleHertzConfig;
use fs_tribo::sliding_patch::contact::multiple::rolling::{
    HertzRollingContactInput, MultipleHertzRollingConfig, MultipleHertzRollingError,
    MultipleHertzRollingStep, PreparedHertzRollingContacts,
};
use fs_tribo::sliding_patch::contact::{HertzWrenchContact, HertzWrenchMode, HertzWrenchStep};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

/// A uniquely identified original patch and its accepted rolling checkpoint.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingIslandContact<'a> {
    /// Must equal the original rolling patch identity.
    pub patch_id: &'a str,
    pub contact: RigidHertzRollingContact<'a>,
    pub history: &'a RigidRollingHistory,
}

/// Repeated source trial against the same initial bodies and histories.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingIslandLaw<'a> {
    pub patch_id: &'a str,
    pub law: &'a HertzWrenchContact,
    /// Actual temperature of the deforming material, not interface temperature.
    pub temperature_k: f64,
}

/// Physical tolerances and explicit full-island allocation/iteration bounds.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingIslandControls {
    pub mechanical: RigidHertzRollingControls,
    pub maximum_contacts: usize,
    /// Maximum entries in the full 5N-by-5N response, at most 409600.
    pub maximum_matrix_terms: usize,
    pub maximum_sweeps: usize,
    pub maximum_evaluations: usize,
}
impl RigidHertzRollingIslandControls {
    /// Admit controls before any source trials or per-patch allocations.
    pub fn validate(self) -> Result<(), RigidHertzRollingIslandError> {
        RigidHertzIslandControls {
            mechanical: self.mechanical.mechanical,
            maximum_contacts: self.maximum_contacts,
            maximum_matrix_terms: self.maximum_matrix_terms.min(147_456),
            maximum_sweeps: self.maximum_sweeps,
            maximum_evaluations: self.maximum_evaluations,
        }
        .validate()
        .map_err(|_| invalid("invalid original Hertz island controls"))?;
        if !(25..=409_600).contains(&self.maximum_matrix_terms)
            || [
                self.mechanical.maximum_rolling_moment_nm,
                self.mechanical.rolling_moment_tolerance_nm,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(invalid(
                "bounded full response and finite nonnegative rolling gates required",
            ));
        }
        Ok(())
    }
    fn hertz(self) -> MultipleHertzConfig {
        MultipleHertzConfig {
            maximum_contacts: self.maximum_contacts,
            maximum_matrix_terms: self.maximum_matrix_terms.min(147_456),
            maximum_sweeps: self.maximum_sweeps,
            maximum_evaluations: self.maximum_evaluations,
            local: self.mechanical.mechanical.contact,
            total_energy_tolerance_j: self.mechanical.mechanical.energy_tolerance_j,
        }
    }
}

/// Immutable native geometry in canonical original patch-ID order.
#[derive(Clone, Debug)]
pub struct RigidHertzRollingIslandPatch {
    pub patch_id: String,
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub contact_points_world: [Vec3; 2],
    pub normal_force_n: f64,
    pub contact_radius_m: f64,
    pub sphere_radius_m: f64,
    pub deforming_body: usize,
}

/// One original rolling checkpoint, keyed independently of caller ordering.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingIslandHistory {
    pub patch_id: String,
    pub history: RigidRollingHistory,
}

/// Original material responses re-evaluated on the actual common endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingIslandContactReport {
    pub patch_id: String,
    pub hertz: HertzWrenchStep,
    pub rolling: RollingLossStep,
    pub deforming_body: usize,
    pub rolling_temperature_k: f64,
    pub rolling_excitation_frequency_hz: f64,
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
    pub surface_midpoint_work_j: f64,
    pub rolling_midpoint_work_j: f64,
    pub surface_work_residual_j: f64,
    pub surface_heat_residual_j: f64,
    pub rolling_work_residual_j: f64,
    pub rolling_heat_residual_j: f64,
    pub history_heat_residual_j: f64,
    pub native_force_residual_n: f64,
    pub native_spin_moment_residual_nm: f64,
    pub native_rolling_moment_residual_nm: f64,
    pub maximum_compatibility_residual_m: f64,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_angular_velocity_residual_rad_s: f64,
    pub maximum_assisting_work_j: f64,
}

/// Each native body's energy and momentum occur exactly once in these totals.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingIslandReport {
    pub contacts: Vec<RigidHertzRollingIslandContactReport>,
    pub native: MultipleHertzRollingStep,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub surface_heat_j: f64,
    pub bulk_heat_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
}

/// Complete immutable proposal. No body or history is accepted separately.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingIslandCandidate {
    pub body_indices: Vec<usize>,
    pub bodies_after: Vec<RigidBodyState>,
    pub rolling_histories: Vec<RigidHertzRollingIslandHistory>,
    pub report: RigidHertzRollingIslandReport,
}

#[derive(Clone, Debug)]
struct PatchMechanics<'a> {
    contact: RigidHertzRollingContact<'a>,
    history: RigidRollingHistory,
    bodies: [usize; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    before: [PointKinematics; 2],
    rates_before: [f64; 5],
}

/// Reusable held geometry, full inertia, original histories and rolling factor.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzRollingIsland<'a> {
    patches: Vec<RigidHertzRollingIslandPatch>,
    mechanics: Vec<PatchMechanics<'a>>,
    body_indices: Vec<usize>,
    before: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    totals_before: (f64, Vec3, Vec3),
    dt_s: f64,
    controls: RigidHertzRollingIslandControls,
    mobility: Vec<f64>,
    native: PreparedHertzRollingContacts,
}

#[derive(Debug)]
pub enum RigidHertzRollingIslandError {
    Invalid(&'static str),
    Pair(RigidHertzRollingError),
    Native(MultipleHertzRollingError),
    UnresolvedSlipEvent {
        patch_id: String,
        assisting_work_j: f64,
    },
    RollingMomentBudget {
        patch_id: String,
        required_nm: f64,
        maximum_nm: f64,
    },
    ContactAudit(Box<RigidHertzRollingIslandContactReport>),
    Audit(Box<RigidHertzRollingIslandReport>),
}
impl From<RigidHertzRollingError> for RigidHertzRollingIslandError {
    fn from(e: RigidHertzRollingError) -> Self {
        Self::Pair(e)
    }
}
impl From<MultipleHertzRollingError> for RigidHertzRollingIslandError {
    fn from(e: MultipleHertzRollingError) -> Self {
        Self::Native(e)
    }
}
impl From<RigidFrictionError> for RigidHertzRollingIslandError {
    fn from(e: RigidFrictionError) -> Self {
        Self::Pair(e.into())
    }
}
impl From<DynamicsError> for RigidHertzRollingIslandError {
    fn from(e: DynamicsError) -> Self {
        Self::Pair(e.into())
    }
}
impl From<TriboError> for RigidHertzRollingIslandError {
    fn from(e: TriboError) -> Self {
        Self::Pair(e.into())
    }
}
impl From<fs_tribo::rolling_loss::RollingLossError> for RigidHertzRollingIslandError {
    fn from(e: fs_tribo::rolling_loss::RollingLossError) -> Self {
        Self::Pair(e.into())
    }
}
impl From<crate::tangential_friction::hertz::RigidHertzError> for RigidHertzRollingIslandError {
    fn from(e: crate::tangential_friction::hertz::RigidHertzError) -> Self {
        Self::Pair(e.into())
    }
}
impl core::fmt::Display for RigidHertzRollingIslandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid Hertz/rolling island: {self:?}")
    }
}
impl std::error::Error for RigidHertzRollingIslandError {}
fn invalid(message: &'static str) -> RigidHertzRollingIslandError {
    RigidHertzRollingIslandError::Invalid(message)
}

impl<'a> PreparedRigidHertzRollingIsland<'a> {
    /// Assemble sparse two-body row factors into one original full 5N response.
    #[allow(clippy::too_many_lines)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidHertzRollingIslandContact<'a>],
        dt_s: f64,
        controls: RigidHertzRollingIslandControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzRollingIslandError> {
        poll(gate)?;
        controls.validate()?;
        let count = contacts.len();
        let rows = count
            .checked_mul(5)
            .ok_or(invalid("contact dimension overflow"))?;
        let terms = rows
            .checked_mul(rows)
            .ok_or(invalid("matrix dimension overflow"))?;
        if count == 0
            || count > controls.maximum_contacts
            || terms > controls.maximum_matrix_terms
            || bodies.len() != properties.len()
            || !dt_s.is_finite()
            || dt_s <= 0.0
        {
            return Err(invalid(
                "matching bodies, positive interval and bounded nonempty contacts required",
            ));
        }
        let mut order: Vec<_> = (0..count).collect();
        order.sort_by(|&a, &b| contacts[a].patch_id.cmp(contacts[b].patch_id));
        if order.iter().any(|&i| {
            contacts[i].patch_id.is_empty()
                || contacts[i].patch_id != contacts[i].contact.patch.patch_id()
        }) || order
            .windows(2)
            .any(|p| contacts[p[0]].patch_id == contacts[p[1]].patch_id)
        {
            return Err(invalid("unique original rolling patch identities required"));
        }
        let mut body_indices = Vec::with_capacity(2 * count);
        for &i in &order {
            let b = contacts[i].contact.binding;
            if b.body_a >= bodies.len() || b.body_b >= bodies.len() || b.body_a == b.body_b {
                return Err(invalid("two distinct valid body indices required"));
            }
            body_indices.extend([b.body_a, b.body_b]);
        }
        body_indices.sort_unstable();
        body_indices.dedup();
        if body_indices.len() > controls.mechanical.mechanical.maximum_bodies {
            return Err(invalid("participating body budget exceeded"));
        }
        let before: Vec<_> = body_indices.iter().map(|&i| bodies[i]).collect();
        let captured_properties: Vec<_> = body_indices.iter().map(|&i| properties[i]).collect();
        let mut patches: Vec<RigidHertzRollingIslandPatch> = Vec::with_capacity(count);
        let mut mechanics: Vec<PatchMechanics<'a>> = Vec::with_capacity(count);
        let mut factors = Vec::with_capacity(count);
        let mut damping = Vec::with_capacity(count);
        let mut force_tolerance = f64::MAX;
        for &index in &order {
            poll(gate)?;
            let identified = contacts[index];
            let c = identified.contact;
            let binding = c.binding;
            let pair = [bodies[binding.body_a], bodies[binding.body_b]];
            let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
            let prepared = PreparedRigidHertzRollingPair::new(
                &pair,
                &pair_properties,
                identified.history,
                RigidHertzRollingContact {
                    binding: RigidContactBinding {
                        body_a: 0,
                        body_b: 1,
                        ..binding
                    },
                    ..c
                },
                dt_s,
                controls.mechanical,
                gate,
            )?;
            let core = &prepared.core;
            let local = [
                body_indices
                    .binary_search(&binding.body_a)
                    .expect("captured A"),
                body_indices
                    .binary_search(&binding.body_b)
                    .expect("captured B"),
            ];
            let radius = prepared.radius_m;
            let hertz_radius = core.contact_radius_m();
            let deforming_body = if prepared.deforming_body == 0 {
                binding.body_a
            } else {
                binding.body_b
            };
            let sphere_arm = core.points_before[prepared.deforming_body].arm_world;
            let card = sphere_card(c.rolling)?;
            for (previous, previous_mechanics) in patches.iter().zip(&mechanics) {
                if previous.deforming_body != deforming_body {
                    continue;
                }
                let previous_card = sphere_card(previous_mechanics.contact.rolling)?;
                // Ordered surface labels may differ; the same physical sphere
                // still has one immutable material and one held bulk state.
                let mut parameters = card.parameters();
                let previous_parameters = previous_card.parameters();
                parameters.sphere_surface = previous_parameters.sphere_surface;
                if parameters != previous_parameters
                    || card.source_id() != previous_card.source_id()
                    || card.model_id() != previous_card.model_id()
                    || card.authority() != previous_card.authority()
                    || c.temperature_k.to_bits()
                        != previous_mechanics.contact.temperature_k.to_bits()
                {
                    return Err(invalid(
                        "the same deforming sphere requires one material, radius and bulk temperature",
                    ));
                }
                if core.normal_force_n() > 0.0 && previous.normal_force_n > 0.0 {
                    let previous_slot = if previous.binding.body_a == deforming_body {
                        0
                    } else {
                        1
                    };
                    let previous_arm = previous_mechanics.before[previous_slot].arm_world;
                    let angular_width = (hertz_radius / radius).asin()
                        + (previous.contact_radius_m / radius).asin();
                    let footprint_separation = 2.0 * radius * (0.5 * angular_width).sin();
                    if !footprint_separation.is_finite()
                        || magnitude(sphere_arm.sub(previous_arm))
                            <= footprint_separation
                                + controls.mechanical.mechanical.geometry_tolerance_m
                    {
                        return Err(invalid(
                            "distinct loaded patches on one sphere require disjoint spherical footprints",
                        ));
                    }
                }
            }
            let original_damping =
                (sphere_card(c.rolling)?.rolling_time_s() * core.normal_force_n()) / radius;
            let d = finite(original_damping / dt_s)?;
            if sphere_card(c.rolling)?.rolling_time_s() > 0.0
                && core.normal_force_n() > 0.0
                && d == 0.0
            {
                return Err(invalid("unrepresentable original material rolling damping"));
            }
            damping.push(d);
            force_tolerance = force_tolerance.min(finite(
                controls.mechanical.rolling_moment_tolerance_nm / radius,
            )?);
            let mut factor = [[[0.0; 6]; 5]; 2];
            for (side, sign) in [(0, 1.0), (1, -1.0)] {
                let inverse_mass_root = pair_properties[side].mass().sqrt().recip();
                let orientation = pair[side].pose().orientation();
                for axis in 0..5 {
                    let direction_world = core.axes[if axis < 3 { axis } else { axis - 3 }];
                    let direction = orientation.rotate_world_to_body(direction_world);
                    let translation = if axis < 2 {
                        components(direction_world.scale(sign * inverse_mass_root))
                    } else {
                        [0.0; 3]
                    };
                    let torque = if axis < 2 {
                        core.arms[side].cross(direction).scale(sign)
                    } else {
                        direction.scale(sign * if axis == 2 { hertz_radius } else { radius })
                    };
                    let rotation =
                        components(pair_properties[side].angular_mobility_factor(torque));
                    factor[side][axis][..3].copy_from_slice(&translation);
                    factor[side][axis][3..].copy_from_slice(&rotation);
                }
            }
            factors.push(factor);
            patches.push(RigidHertzRollingIslandPatch {
                patch_id: identified.patch_id.to_owned(),
                binding,
                frame: c.frame,
                contact_points_world: core.contact_points_world(),
                normal_force_n: core.normal_force_n(),
                contact_radius_m: hertz_radius,
                sphere_radius_m: radius,
                deforming_body,
            });
            mechanics.push(PatchMechanics {
                contact: c,
                history: identified.history.clone(),
                bodies: local,
                arms: core.arms,
                axes: core.axes,
                before: core.points_before,
                rates_before: prepared.rate_before,
            });
        }
        let mut mobility = vec![0.0; terms];
        for i in 0..rows {
            poll(gate)?;
            for j in 0..=i {
                let mut value = 0.0;
                for side_i in 0..2 {
                    for side_j in 0..2 {
                        if mechanics[i / 5].bodies[side_i] == mechanics[j / 5].bodies[side_j] {
                            value += factors[i / 5][side_i][i % 5]
                                .iter()
                                .zip(&factors[j / 5][side_j][j % 5])
                                .map(|(a, b)| a * b)
                                .sum::<f64>();
                        }
                    }
                }
                let value = finite(value)?;
                mobility[i * rows + j] = value;
                mobility[j * rows + i] = value;
            }
        }
        let dt_factor = finite(0.5 * dt_s * dt_s)?;
        let native = PreparedHertzRollingContacts::new(
            count,
            mobility.iter().map(|v| v * dt_factor).collect(),
            damping,
            MultipleHertzRollingConfig {
                hertz: controls.hertz(),
                maximum_matrix_terms: controls.maximum_matrix_terms,
                rolling_force_tolerance_n: force_tolerance,
                total_energy_tolerance_j: controls.mechanical.mechanical.energy_tolerance_j,
            },
            || gate.is_requested(),
        )?;
        let totals_before = totals(&before, &captured_properties, gate)?;
        poll(gate)?;
        Ok(Self {
            patches,
            mechanics,
            body_indices,
            before,
            properties: captured_properties,
            totals_before,
            dt_s,
            controls,
            mobility,
            native,
        })
    }

    #[must_use]
    pub fn contacts(&self) -> &[RigidHertzRollingIslandPatch] {
        &self.patches
    }
    pub fn patch_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.patches.iter().map(|p| p.patch_id.as_str())
    }
    #[must_use]
    pub fn body_indices(&self) -> &[usize] {
        &self.body_indices
    }
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }
    #[must_use]
    pub fn scaled_mobility_per_kg(&self) -> &[f64] {
        &self.mobility
    }
    #[must_use]
    pub fn scaled_compliance_m_per_n(&self) -> &[f64] {
        self.native.scaled_compliance_m_per_n()
    }

    /// Solve every keyed law against the same captured states; no trial warms
    /// the next trial with an unaccepted history or body endpoint.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        laws: &[RigidHertzRollingIslandLaw<'_>],
        gate: &CancelGate,
    ) -> Result<RigidHertzRollingIslandCandidate, RigidHertzRollingIslandError> {
        poll(gate)?;
        if laws.len() != self.patches.len() {
            return Err(invalid("one law/temperature per original patch required"));
        }
        let mut ordered: Vec<_> = laws.iter().collect();
        ordered.sort_by(|a, b| a.patch_id.cmp(b.patch_id));
        let mut inputs = Vec::with_capacity(ordered.len());
        for (index, ((keyed, patch), m)) in ordered
            .iter()
            .zip(&self.patches)
            .zip(&self.mechanics)
            .enumerate()
        {
            poll(gate)?;
            let p = keyed.law.parameters();
            if keyed.patch_id != patch.patch_id
                || p.normal_force_n.to_bits() != patch.normal_force_n.to_bits()
                || p.contact_radius_m.to_bits() != patch.contact_radius_m.to_bits()
            {
                return Err(invalid(
                    "prepared patch identity, load and radius must remain unchanged",
                ));
            }
            for j in 0..index {
                if self.patches[j].deforming_body == patch.deforming_body
                    && ordered[j].temperature_k.to_bits() != keyed.temperature_k.to_bits()
                {
                    return Err(invalid(
                        "one bulk temperature per original deforming sphere required",
                    ));
                }
            }
            let c = RigidHertzRollingContact {
                temperature_k: keyed.temperature_k,
                ..m.contact
            };
            sphere_card(c.rolling)?.response(
                c.patch,
                c.interface,
                kinematics(
                    c,
                    patch.sphere_radius_m,
                    norm([m.rates_before[3], m.rates_before[4]]) / patch.sphere_radius_m,
                    m.rates_before[2] / patch.contact_radius_m,
                    self.dt_s,
                )?,
            )?;
            inputs.push(HertzRollingContactInput {
                law: keyed.law,
                free_coordinates_m: m.rates_before.map(|v| self.dt_s * v),
            });
        }
        let phases: Vec<_> = self
            .mechanics
            .iter()
            .map(|m| accepted_phase(&m.rates_before[..3]))
            .collect();
        let native = self
            .native
            .advance_phases(&inputs, &phases, || gate.is_requested())?;
        let mut linear = vec![Vec3::ZERO; self.before.len()];
        let mut angular = vec![Vec3::ZERO; self.before.len()];
        for (i, m) in self.mechanics.iter().enumerate() {
            poll(gate)?;
            let hertz = &native.hertz.contacts[i];
            let rolling = native.rolling[i]
                .scaled_force_n
                .map(|v| v * self.patches[i].sphere_radius_m);
            let required = norm(rolling);
            if required > self.controls.mechanical.maximum_rolling_moment_nm {
                return Err(RigidHertzRollingIslandError::RollingMomentBudget {
                    patch_id: self.patches[i].patch_id.clone(),
                    required_nm: required,
                    maximum_nm: self.controls.mechanical.maximum_rolling_moment_nm,
                });
            }
            let impulse = world(hertz.force_n, m.axes).scale(self.dt_s);
            let couple = m.axes[2]
                .scale(hertz.spin_moment_nm)
                .add(world(rolling, m.axes))
                .scale(self.dt_s);
            for (side, sign) in [(0, 1.0), (1, -1.0)] {
                let body = m.bodies[side];
                let orientation = self.before[body].pose().orientation();
                let body_impulse = orientation.rotate_world_to_body(impulse.scale(sign));
                let body_couple = orientation.rotate_world_to_body(couple.scale(sign));
                linear[body] = linear[body].add(impulse.scale(sign));
                angular[body] = angular[body]
                    .add(m.arms[side].cross(body_impulse))
                    .add(body_couple);
            }
        }
        let mut after = Vec::with_capacity(self.before.len());
        for (i, b) in self.before.iter().enumerate() {
            poll(gate)?;
            after.push(RigidBodyState::new(
                b.pose(),
                b.linear_momentum_world().add(linear[i]),
                b.angular_momentum_body().add(angular[i]),
            )?);
        }
        let mut reports = Vec::with_capacity(self.patches.len());
        let mut rolling_histories = Vec::with_capacity(self.patches.len());
        for (i, ((patch, m), keyed)) in self
            .patches
            .iter()
            .zip(&self.mechanics)
            .zip(&ordered)
            .enumerate()
        {
            poll(gate)?;
            let c = RigidHertzRollingContact {
                temperature_k: keyed.temperature_k,
                ..m.contact
            };
            let hertz = &native.hertz.contacts[i];
            let current = points(
                [after[m.bodies[0]], after[m.bodies[1]]],
                [self.properties[m.bodies[0]], self.properties[m.bodies[1]]],
                m.arms,
            )?;
            let point_before = relative_point(m.before);
            let point_after = relative_point(current);
            let omega_before = relative_spin(m.before);
            let omega_after = relative_spin(current);
            let midpoint = point_before.add(point_after).scale(0.5);
            let omega_mid = omega_before.add(omega_after).scale(0.5);
            let actual_roll = [omega_mid.dot(m.axes[0]), omega_mid.dot(m.axes[1])];
            let rolling = c.rolling.advance(
                c.patch,
                c.interface,
                kinematics(
                    c,
                    patch.sphere_radius_m,
                    norm(actual_roll),
                    omega_mid.dot(m.axes[2]),
                    self.dt_s,
                )?,
                c.work,
                m.history.state(),
            )?;
            sphere_card(c.rolling)?.response(
                c.patch,
                c.interface,
                kinematics(
                    c,
                    patch.sphere_radius_m,
                    norm([omega_after.dot(m.axes[0]), omega_after.dot(m.axes[1])]),
                    omega_after.dot(m.axes[2]),
                    self.dt_s,
                )?,
            )?;
            if rolling.wrench.contour_force_n != 0.0
                || rolling.wrench.spin_moment_nm != 0.0
                || rolling.stored_energy_j != 0.0
                || rolling.storage_change_j != 0.0
            {
                return Err(invalid(
                    "original sphere rolling requires its sole zero-storage moment channel",
                ));
            }
            let rolling_moment = native.rolling[i]
                .scaled_force_n
                .map(|v| v * patch.sphere_radius_m);
            let native_rolling = along(actual_roll, rolling.wrench.rolling_moment_nm);
            let rolling_residual = norm([
                native_rolling[0] - rolling_moment[0],
                native_rolling[1] - rolling_moment[1],
            ]);
            let force = world(hertz.force_n, m.axes);
            let spin_couple = m.axes[2].scale(hertz.spin_moment_nm);
            let rolling_couple = world(rolling_moment, m.axes);
            let impulse = force.scale(self.dt_s);
            let spin_impulse = spin_couple.scale(self.dt_s);
            let rolling_impulse = rolling_couple.scale(self.dt_s);
            let mut surface_heat = 0.0;
            let mut force_residual = 0.0;
            let mut spin_residual = 0.0;
            if hertz.mode != HertzWrenchMode::Sticking {
                let p = keyed.law.parameters();
                let response = HertzSlidingPatch::new(HertzSlidingPatchParameters {
                    normal_force_n: p.normal_force_n,
                    contact_radius_m: p.contact_radius_m,
                    friction_coefficient: p.kinetic_mu,
                })?
                .evaluate(
                    [midpoint.dot(m.axes[0]), midpoint.dot(m.axes[1])],
                    omega_mid.dot(m.axes[2]),
                )?;
                surface_heat = finite(self.dt_s * response.dissipated_power_w)?;
                force_residual = norm([
                    response.force_n[0] - hertz.force_n[0],
                    response.force_n[1] - hertz.force_n[1],
                ]);
                spin_residual = (response.spin_moment_nm - hertz.spin_moment_nm).abs();
            }
            let after_rates = rates(
                current,
                m.axes,
                patch.contact_radius_m,
                patch.sphere_radius_m,
            );
            let actual =
                std::array::from_fn(|j| 0.5 * self.dt_s * (m.rates_before[j] + after_rates[j]));
            let expected = native.scaled_motion_m[i];
            let compatibility = (0..5)
                .map(|j| (actual[j] - expected[j]).abs())
                .fold(0.0, f64::max);
            let endpoint_velocity = (0..2)
                .map(|j| {
                    (after_rates[j] - (2.0 * expected[j] / self.dt_s - m.rates_before[j])).abs()
                })
                .fold(0.0, f64::max);
            let endpoint_angular = (2..5)
                .map(|j| {
                    (after_rates[j] - (2.0 * expected[j] / self.dt_s - m.rates_before[j])).abs()
                        / if j == 2 {
                            patch.contact_radius_m
                        } else {
                            patch.sphere_radius_m
                        }
                })
                .fold(0.0, f64::max);
            let surface_work = impulse.dot(midpoint) + spin_impulse.dot(omega_mid);
            let rolling_work = rolling_impulse.dot(omega_mid);
            let bulk_heat = rolling.dissipation.total_heat_j;
            let history_heat =
                rolling.next_state.cumulative_heat_j() - m.history.state().cumulative_heat_j();
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
                * self.dt_s;
            let report = RigidHertzRollingIslandContactReport {
                patch_id: patch.patch_id.clone(),
                hertz: hertz.clone(),
                deforming_body: patch.deforming_body,
                rolling_temperature_k: keyed.temperature_k,
                rolling_excitation_frequency_hz: c.excitation_frequency_hz,
                scaled_rates_before_m_s: m.rates_before,
                scaled_rates_after_m_s: after_rates,
                actual_midpoint_coordinates_m: actual,
                normal_velocity_before_m_s: point_before.dot(m.axes[2]),
                normal_velocity_after_m_s: point_after.dot(m.axes[2]),
                force_impulse_world_ns: impulse,
                spin_angular_impulse_world_nms: spin_impulse,
                rolling_angular_impulse_world_nms: rolling_impulse,
                surface_heat_j: surface_heat,
                bulk_heat_j: bulk_heat,
                surface_midpoint_work_j: surface_work,
                rolling_midpoint_work_j: rolling_work,
                surface_work_residual_j: surface_work + hertz.contact_work_j,
                surface_heat_residual_j: surface_heat - hertz.dissipated_energy_j,
                rolling_work_residual_j: rolling_work + bulk_heat,
                rolling_heat_residual_j: bulk_heat - native.rolling[i].dissipated_energy_j,
                history_heat_residual_j: history_heat - bulk_heat,
                native_force_residual_n: force_residual,
                native_spin_moment_residual_nm: spin_residual,
                native_rolling_moment_residual_nm: rolling_residual,
                maximum_compatibility_residual_m: compatibility,
                endpoint_velocity_residual_m_s: endpoint_velocity,
                endpoint_angular_velocity_residual_rad_s: endpoint_angular,
                maximum_assisting_work_j: assisting_work,
                rolling,
            };
            audit_contact(&report, self.controls.mechanical)?;
            if patch.normal_force_n > 0.0
                && report.normal_velocity_after_m_s.abs()
                    > self.controls.mechanical.mechanical.velocity_tolerance_m_s
            {
                return Err(invalid(
                    "shared joint recoil induced unresolved normal motion",
                ));
            }
            if assisting_work > self.controls.mechanical.mechanical.energy_tolerance_j {
                return Err(RigidHertzRollingIslandError::UnresolvedSlipEvent {
                    patch_id: patch.patch_id.clone(),
                    assisting_work_j: assisting_work,
                });
            }
            rolling_histories.push(RigidHertzRollingIslandHistory {
                patch_id: patch.patch_id.clone(),
                history: RigidRollingHistory::from_native_step(&report.rolling),
            });
            reports.push(report);
        }
        let (energy_before, linear_before, angular_before) = self.totals_before;
        let (energy_after, linear_after, angular_after) = totals(&after, &self.properties, gate)?;
        let surface_heat = reports.iter().map(|r| r.surface_heat_j).sum::<f64>();
        let bulk_heat = reports.iter().map(|r| r.bulk_heat_j).sum::<f64>();
        let work = reports
            .iter()
            .map(|r| r.surface_midpoint_work_j + r.rolling_midpoint_work_j)
            .sum::<f64>();
        let report = RigidHertzRollingIslandReport {
            contacts: reports,
            native,
            kinetic_before_j: energy_before,
            kinetic_after_j: energy_after,
            surface_heat_j: surface_heat,
            bulk_heat_j: bulk_heat,
            midpoint_impulse_work_j: work,
            work_energy_residual_j: energy_after - energy_before - work,
            total_energy_residual_j: energy_after - energy_before + surface_heat + bulk_heat,
            linear_momentum_change_ns: linear_after.sub(linear_before),
            angular_momentum_change_nms: angular_after.sub(angular_before),
        };
        let m = self.controls.mechanical.mechanical;
        if [
            report.work_energy_residual_j,
            report.total_energy_residual_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > m.energy_tolerance_j)
            || [
                surface_heat,
                bulk_heat,
                energy_before,
                energy_after,
                work,
                magnitude(report.linear_momentum_change_ns),
                magnitude(report.angular_momentum_change_nms),
            ]
            .into_iter()
            .any(|v| !v.is_finite())
            || surface_heat < 0.0
            || bulk_heat < 0.0
            || magnitude(report.linear_momentum_change_ns) > m.linear_momentum_tolerance_ns
            || magnitude(report.angular_momentum_change_nms) > m.angular_momentum_tolerance_nms
        {
            return Err(RigidHertzRollingIslandError::Audit(Box::new(report)));
        }
        poll(gate)?;
        Ok(RigidHertzRollingIslandCandidate {
            body_indices: self.body_indices.clone(),
            bodies_after: after,
            rolling_histories,
            report,
        })
    }
}

fn audit_contact(
    r: &RigidHertzRollingIslandContactReport,
    c: RigidHertzRollingControls,
) -> Result<(), RigidHertzRollingIslandError> {
    let m = c.mechanical;
    if [
        r.surface_work_residual_j,
        r.surface_heat_residual_j,
        r.rolling_work_residual_j,
        r.rolling_heat_residual_j,
        r.history_heat_residual_j,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || v.abs() > m.energy_tolerance_j)
        || [
            r.surface_heat_j,
            r.bulk_heat_j,
            r.surface_midpoint_work_j,
            r.rolling_midpoint_work_j,
            r.maximum_assisting_work_j,
            r.normal_velocity_before_m_s,
            r.normal_velocity_after_m_s,
            r.native_force_residual_n,
            r.native_spin_moment_residual_nm,
            r.native_rolling_moment_residual_nm,
            r.maximum_compatibility_residual_m,
            r.endpoint_velocity_residual_m_s,
            r.endpoint_angular_velocity_residual_rad_s,
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
    {
        return Err(RigidHertzRollingIslandError::ContactAudit(Box::new(
            r.clone(),
        )));
    }
    Ok(())
}
