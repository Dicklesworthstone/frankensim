//! Simultaneous Hertz shear/spin patches on a bounded native body island.
//!
//! All patches see the same original body momenta and the full shared-body
//! response. A single summed impulse updates each participating body. Stable
//! patch identities determine assembly and reduction order; reordering caller
//! inputs does not change a solve. Geometry, normal loads and radii are held,
//! and the surrounding trajectory/normal owner must refine unresolved events.
//! The original constitutive owner admits positive-semidefinite full response,
//! including dependent wrench coordinates of distinct patches on the same pair.
//! No invented compliance or sequential accepted-patch sweep is introduced.

use super::{
    PreparedRigidHertzPair, RigidHertzContact, RigidHertzControls, RigidHertzError, finite, points,
    relative_point, relative_spin,
};
use crate::tangential_friction::{RigidContactBinding, components, magnitude, poll, totals};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::TriboError;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::multiple::{
    HertzContactInput, MultipleHertzConfig, MultipleHertzError, MultipleHertzStep,
    PreparedHertzContacts,
};
use fs_tribo::sliding_patch::contact::{HertzWrenchContact, HertzWrenchMode, HertzWrenchStep};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

/// One identified original patch. IDs are unique within the prepared island.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzIslandContact<'a> {
    pub patch_id: &'a str,
    pub contact: RigidHertzContact<'a>,
}

/// A source-resolved coefficient trial, keyed independently of input ordering.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzIslandLaw<'a> {
    pub patch_id: &'a str,
    pub law: &'a HertzWrenchContact,
}

/// Physical mechanics gates and bounded global constitutive work.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzIslandControls {
    /// `maximum_bodies` bounds participating bodies, not untouched spectators.
    pub mechanical: RigidHertzControls,
    pub maximum_contacts: usize,
    pub maximum_matrix_terms: usize,
    pub maximum_sweeps: usize,
    pub maximum_evaluations: usize,
}

impl RigidHertzIslandControls {
    /// Admit all allocation/work bounds before a coupled owner resolves sources
    /// or allocates per-contact trial data.
    pub fn validate(self) -> Result<(), RigidHertzIslandError> {
        self.constitutive().validate()?;
        if self.mechanical.maximum_bodies < 2
            || [
                self.mechanical.geometry_tolerance_m,
                self.mechanical.velocity_tolerance_m_s,
                self.mechanical.angular_velocity_tolerance_rad_s,
                self.mechanical.linear_momentum_tolerance_ns,
                self.mechanical.angular_momentum_tolerance_nms,
                self.mechanical.energy_tolerance_j,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err(RigidHertzIslandError::Invalid(
                "bounded native bodies and finite nonnegative mechanical gates required",
            ));
        }
        Ok(())
    }

    fn constitutive(self) -> MultipleHertzConfig {
        MultipleHertzConfig {
            maximum_contacts: self.maximum_contacts,
            maximum_matrix_terms: self.maximum_matrix_terms,
            maximum_sweeps: self.maximum_sweeps,
            maximum_evaluations: self.maximum_evaluations,
            local: self.mechanical.contact,
            total_energy_tolerance_j: self.mechanical.energy_tolerance_j,
        }
    }
}

/// Read-only original geometry and normal data in canonical patch-ID order.
#[derive(Clone, Debug)]
pub struct RigidHertzIslandPatch {
    pub patch_id: String,
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub contact_points_world: [Vec3; 2],
    pub normal_force_n: f64,
    pub contact_radius_m: f64,
}

/// Per-patch work evaluated from the same actual shared-body endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzIslandContactReport {
    pub patch_id: String,
    pub contact: HertzWrenchStep,
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    pub spin_before_rad_s: f64,
    pub spin_after_rad_s: f64,
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub actual_midpoint_displacement_m: [f64; 2],
    pub actual_midpoint_spin_angle_rad: f64,
    pub force_impulse_world_ns: Vec3,
    pub spin_angular_impulse_world_nms: Vec3,
    pub midpoint_impulse_work_j: f64,
    pub dissipated_energy_j: f64,
    pub contact_work_residual_j: f64,
    pub contact_heat_residual_j: f64,
    pub native_force_residual_n: f64,
    pub native_spin_moment_residual_nm: f64,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_spin_residual_rad_s: f64,
    pub maximum_assisting_work_j: f64,
}

/// Global totals count each participating native body exactly once.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzIslandReport {
    pub contacts: Vec<RigidHertzIslandContactReport>,
    pub native: MultipleHertzStep,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub dissipated_energy_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
}

/// All candidate momenta and receipts, without any accepted-state mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzIslandCandidate {
    /// Sorted original indices, corresponding exactly to `bodies_after`.
    pub body_indices: Vec<usize>,
    pub bodies_after: Vec<RigidBodyState>,
    pub report: RigidHertzIslandReport,
}

#[derive(Clone, Debug)]
struct PatchMechanics {
    bodies: [usize; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    before: [PointKinematics; 2],
    rates_before: [f64; 3],
}

/// Shared geometry/full-inertia response prepared once for coefficient trials.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzIsland {
    patches: Vec<RigidHertzIslandPatch>,
    mechanics: Vec<PatchMechanics>,
    body_indices: Vec<usize>,
    before: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    totals_before: (f64, Vec3, Vec3),
    dt_s: f64,
    controls: RigidHertzIslandControls,
    mobility: Vec<f64>,
    native: PreparedHertzContacts,
}

/// Every failure leaves all input bodies and reusable preparation unchanged.
#[derive(Debug)]
pub enum RigidHertzIslandError {
    Invalid(&'static str),
    Hertz(RigidHertzError),
    Native(MultipleHertzError),
    /// A held effort assists an actual endpoint: refine the trajectory interval.
    UnresolvedSlipEvent {
        patch_id: String,
        assisting_work_j: f64,
    },
    ContactAudit(Box<RigidHertzIslandContactReport>),
    Audit(Box<RigidHertzIslandReport>),
}
impl From<RigidHertzError> for RigidHertzIslandError {
    fn from(e: RigidHertzError) -> Self {
        Self::Hertz(e)
    }
}
impl From<DynamicsError> for RigidHertzIslandError {
    fn from(e: DynamicsError) -> Self {
        Self::Hertz(e.into())
    }
}
impl From<TriboError> for RigidHertzIslandError {
    fn from(e: TriboError) -> Self {
        Self::Hertz(e.into())
    }
}
impl From<crate::tangential_friction::RigidFrictionError> for RigidHertzIslandError {
    fn from(e: crate::tangential_friction::RigidFrictionError) -> Self {
        Self::Hertz(e.into())
    }
}
impl From<MultipleHertzError> for RigidHertzIslandError {
    fn from(e: MultipleHertzError) -> Self {
        Self::Native(e)
    }
}
impl core::fmt::Display for RigidHertzIslandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid Hertz island: {self:?}")
    }
}
impl std::error::Error for RigidHertzIslandError {}

impl PreparedRigidHertzIsland {
    /// Assemble the original shared-body Gram matrix in stable patch-ID order.
    /// Only the participating body states/properties are captured. The sparse
    /// two-body row factors avoid a dense contacts-by-all-bodies Jacobian.
    #[allow(clippy::too_many_lines)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidHertzIslandContact<'_>],
        dt_s: f64,
        controls: RigidHertzIslandControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzIslandError> {
        poll(gate)?;
        controls.validate()?;
        let count = contacts.len();
        let rows = count
            .checked_mul(3)
            .ok_or(RigidHertzIslandError::Invalid("contact dimension overflow"))?;
        let terms = rows
            .checked_mul(rows)
            .ok_or(RigidHertzIslandError::Invalid("matrix dimension overflow"))?;
        if count == 0
            || count > controls.maximum_contacts
            || terms > controls.maximum_matrix_terms
            || bodies.len() != properties.len()
            || !dt_s.is_finite()
            || dt_s <= 0.0
        {
            return Err(RigidHertzIslandError::Invalid(
                "matching native bodies, positive interval and bounded nonempty contacts required",
            ));
        }
        let mut order: Vec<_> = (0..count).collect();
        order.sort_by(|&a, &b| contacts[a].patch_id.cmp(contacts[b].patch_id));
        if order.iter().any(|&i| contacts[i].patch_id.is_empty())
            || order
                .windows(2)
                .any(|p| contacts[p[0]].patch_id == contacts[p[1]].patch_id)
        {
            return Err(RigidHertzIslandError::Invalid(
                "unique nonempty patch identities required",
            ));
        }
        let mut body_indices = Vec::with_capacity(count.saturating_mul(2));
        for &i in &order {
            let binding = contacts[i].contact.binding;
            if binding.body_a >= bodies.len()
                || binding.body_b >= bodies.len()
                || binding.body_a == binding.body_b
            {
                return Err(RigidHertzIslandError::Invalid(
                    "two distinct valid native body indices required",
                ));
            }
            body_indices.extend([binding.body_a, binding.body_b]);
        }
        body_indices.sort_unstable();
        body_indices.dedup();
        if body_indices.len() > controls.mechanical.maximum_bodies {
            return Err(RigidHertzIslandError::Invalid(
                "participating body budget exceeded",
            ));
        }
        let before: Vec<_> = body_indices.iter().map(|&i| bodies[i]).collect();
        let captured_properties: Vec<_> = body_indices.iter().map(|&i| properties[i]).collect();
        let mut patches = Vec::with_capacity(count);
        let mut mechanics = Vec::with_capacity(count);
        let mut factors = Vec::with_capacity(count);
        for &index in &order {
            poll(gate)?;
            let identified = contacts[index];
            let contact = identified.contact;
            let binding = contact.binding;
            let pair = [bodies[binding.body_a], bodies[binding.body_b]];
            let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
            // Reuse original native pair admission, without retaining duplicate
            // body snapshots or allowing spectator state into any energy sum.
            let pair_contact = RigidHertzContact {
                binding: RigidContactBinding {
                    body_a: 0,
                    body_b: 1,
                    ..binding
                },
                ..contact
            };
            let prepared = PreparedRigidHertzPair::new(
                &pair,
                &pair_properties,
                pair_contact,
                dt_s,
                controls.mechanical,
                gate,
            )?;
            let local = [
                body_indices
                    .binary_search(&binding.body_a)
                    .expect("captured A"),
                body_indices
                    .binary_search(&binding.body_b)
                    .expect("captured B"),
            ];
            let radius = prepared.contact_radius_m();
            let relative = relative_point(prepared.points_before);
            let normal_rate = relative.dot(prepared.axes[2]);
            if prepared.normal_force_n() > 0.0
                && normal_rate.abs() > controls.mechanical.velocity_tolerance_m_s
            {
                return Err(RigidHertzIslandError::Invalid(
                    "loaded held patch requires closed original normal motion",
                ));
            }
            let mut factor = [[[0.0; 6]; 3]; 2];
            for (side, sign) in [(0, 1.0), (1, -1.0)] {
                let inverse_mass_root = pair_properties[side].mass().sqrt().recip();
                let orientation = pair[side].pose().orientation();
                for axis in 0..3 {
                    let direction = orientation.rotate_world_to_body(prepared.axes[axis]);
                    let translation = if axis < 2 {
                        components(prepared.axes[axis].scale(sign * inverse_mass_root))
                    } else {
                        [0.0; 3]
                    };
                    let torque = if axis < 2 {
                        prepared.arms[side].cross(direction).scale(sign)
                    } else {
                        direction.scale(sign * radius)
                    };
                    let rotation =
                        components(pair_properties[side].angular_mobility_factor(torque));
                    factor[side][axis][..3].copy_from_slice(&translation);
                    factor[side][axis][3..].copy_from_slice(&rotation);
                }
            }
            factors.push(factor);
            mechanics.push(PatchMechanics {
                bodies: local,
                arms: prepared.arms,
                axes: prepared.axes,
                before: prepared.points_before,
                rates_before: [
                    relative.dot(prepared.axes[0]),
                    relative.dot(prepared.axes[1]),
                    radius * relative_spin(prepared.points_before).dot(prepared.axes[2]),
                ],
            });
            patches.push(RigidHertzIslandPatch {
                patch_id: identified.patch_id.to_owned(),
                binding,
                frame: contact.frame,
                contact_points_world: prepared.contact_points_world(),
                normal_force_n: prepared.normal_force_n(),
                contact_radius_m: radius,
            });
        }
        let mut mobility = vec![0.0; terms];
        for i in 0..rows {
            poll(gate)?;
            for j in 0..=i {
                let mut value = 0.0;
                for side_i in 0..2 {
                    for side_j in 0..2 {
                        if mechanics[i / 3].bodies[side_i] == mechanics[j / 3].bodies[side_j] {
                            value += factors[i / 3][side_i][i % 3]
                                .iter()
                                .zip(&factors[j / 3][side_j][j % 3])
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
        let factor = finite(0.5 * dt_s * dt_s)?;
        let compliance: Vec<_> = mobility.iter().map(|&w| w * factor).collect();
        let native =
            PreparedHertzContacts::new(count, compliance, controls.constitutive(), || {
                gate.is_requested()
            })?;
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
    pub fn contacts(&self) -> &[RigidHertzIslandPatch] {
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
    /// Full row-major 3N by 3N native shared-body mobility.
    #[must_use]
    pub fn scaled_mobility_per_kg(&self) -> &[f64] {
        &self.mobility
    }
    #[must_use]
    pub fn scaled_compliance_m_per_n(&self) -> &[f64] {
        self.native.scaled_compliance_m_per_n()
    }

    /// Solve all keyed coefficient trials against the same captured endpoints.
    /// Original N and radius must match exactly for every patch. No proposal is
    /// reused as the next trial's accepted state or as hidden warm-start memory.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        laws: &[RigidHertzIslandLaw<'_>],
        gate: &CancelGate,
    ) -> Result<RigidHertzIslandCandidate, RigidHertzIslandError> {
        poll(gate)?;
        if laws.len() != self.patches.len() {
            return Err(RigidHertzIslandError::Invalid(
                "one original law per prepared patch required",
            ));
        }
        let mut ordered: Vec<_> = laws.iter().collect();
        ordered.sort_by(|a, b| a.patch_id.cmp(b.patch_id));
        let mut inputs = Vec::with_capacity(ordered.len());
        for ((keyed, patch), mechanics) in ordered.iter().zip(&self.patches).zip(&self.mechanics) {
            let parameters = keyed.law.parameters();
            if keyed.patch_id != patch.patch_id
                || parameters.normal_force_n.to_bits() != patch.normal_force_n.to_bits()
                || parameters.contact_radius_m.to_bits() != patch.contact_radius_m.to_bits()
            {
                return Err(RigidHertzIslandError::Invalid(
                    "prepared patch identity, normal force and radius must remain unchanged",
                ));
            }
            inputs.push(HertzContactInput {
                law: keyed.law,
                free_displacement_m: [
                    self.dt_s * mechanics.rates_before[0],
                    self.dt_s * mechanics.rates_before[1],
                ],
                free_spin_angle_rad: self.dt_s * mechanics.rates_before[2] / patch.contact_radius_m,
            });
        }
        let native = self.native.advance(&inputs, || gate.is_requested())?;
        let mut linear = vec![Vec3::ZERO; self.before.len()];
        let mut angular = vec![Vec3::ZERO; self.before.len()];
        for (mechanics, contact) in self.mechanics.iter().zip(&native.contacts) {
            poll(gate)?;
            let impulse = mechanics.axes[0]
                .scale(self.dt_s * contact.force_n[0])
                .add(mechanics.axes[1].scale(self.dt_s * contact.force_n[1]));
            let couple = mechanics.axes[2].scale(self.dt_s * contact.spin_moment_nm);
            for (side, sign) in [(0, 1.0), (1, -1.0)] {
                let body = mechanics.bodies[side];
                let orientation = self.before[body].pose().orientation();
                let body_impulse = orientation.rotate_world_to_body(impulse.scale(sign));
                let body_couple = orientation.rotate_world_to_body(couple.scale(sign));
                linear[body] = linear[body].add(impulse.scale(sign));
                angular[body] = angular[body]
                    .add(mechanics.arms[side].cross(body_impulse))
                    .add(body_couple);
            }
        }
        let mut after = Vec::with_capacity(self.before.len());
        for (i, before) in self.before.iter().enumerate() {
            poll(gate)?;
            after.push(RigidBodyState::new(
                before.pose(),
                before.linear_momentum_world().add(linear[i]),
                before.angular_momentum_body().add(angular[i]),
            )?);
        }
        let mut reports = Vec::with_capacity(self.patches.len());
        for (((patch, mechanics), keyed), contact) in self
            .patches
            .iter()
            .zip(&self.mechanics)
            .zip(&ordered)
            .zip(&native.contacts)
        {
            poll(gate)?;
            let current = points(
                [after[mechanics.bodies[0]], after[mechanics.bodies[1]]],
                [
                    self.properties[mechanics.bodies[0]],
                    self.properties[mechanics.bodies[1]],
                ],
                mechanics.arms,
            )?;
            let relative_before = relative_point(mechanics.before);
            let relative_after = relative_point(current);
            let spin_before = relative_spin(mechanics.before);
            let spin_after = relative_spin(current);
            let midpoint = relative_before.add(relative_after).scale(0.5);
            let midpoint_spin = spin_before
                .add(spin_after)
                .scale(0.5)
                .dot(mechanics.axes[2]);
            let force = mechanics.axes[0]
                .scale(contact.force_n[0])
                .add(mechanics.axes[1].scale(contact.force_n[1]));
            let couple = mechanics.axes[2].scale(contact.spin_moment_nm);
            let impulse = force.scale(self.dt_s);
            let angular_impulse = couple.scale(self.dt_s);
            let midpoint_work =
                impulse.dot(midpoint) + self.dt_s * contact.spin_moment_nm * midpoint_spin;
            let slip_after = [
                relative_after.dot(mechanics.axes[0]),
                relative_after.dot(mechanics.axes[1]),
            ];
            let mut heat = 0.0;
            let mut force_residual = 0.0;
            let mut moment_residual = 0.0;
            if contact.mode != HertzWrenchMode::Sticking {
                let parameters = keyed.law.parameters();
                let response = HertzSlidingPatch::new(HertzSlidingPatchParameters {
                    normal_force_n: parameters.normal_force_n,
                    contact_radius_m: parameters.contact_radius_m,
                    friction_coefficient: parameters.kinetic_mu,
                })?
                .evaluate(
                    [
                        midpoint.dot(mechanics.axes[0]),
                        midpoint.dot(mechanics.axes[1]),
                    ],
                    midpoint_spin,
                )?;
                heat = finite(self.dt_s * response.dissipated_power_w)?;
                force_residual = (response.force_n[0] - contact.force_n[0])
                    .hypot(response.force_n[1] - contact.force_n[1]);
                moment_residual = (response.spin_moment_nm - contact.spin_moment_nm).abs();
            }
            let endpoint_velocity_residual = (slip_after[0]
                - (2.0 * contact.slip_displacement_m[0] / self.dt_s - mechanics.rates_before[0]))
                .abs()
                .max(
                    (slip_after[1]
                        - (2.0 * contact.slip_displacement_m[1] / self.dt_s
                            - mechanics.rates_before[1]))
                        .abs(),
                );
            let endpoint_spin_residual = (spin_after.dot(mechanics.axes[2])
                - (2.0 * contact.spin_angle_rad / self.dt_s
                    - mechanics.rates_before[2] / patch.contact_radius_m))
                .abs();
            let assisting_work = [
                force.dot(relative_before),
                force.dot(relative_after),
                couple.dot(spin_before),
                couple.dot(spin_after),
            ]
            .into_iter()
            .fold(0.0, f64::max)
                * self.dt_s;
            let report = RigidHertzIslandContactReport {
                patch_id: patch.patch_id.clone(),
                contact: *contact,
                slip_velocity_before_m_s: [mechanics.rates_before[0], mechanics.rates_before[1]],
                slip_velocity_after_m_s: slip_after,
                spin_before_rad_s: spin_before.dot(mechanics.axes[2]),
                spin_after_rad_s: spin_after.dot(mechanics.axes[2]),
                normal_velocity_before_m_s: relative_before.dot(mechanics.axes[2]),
                normal_velocity_after_m_s: relative_after.dot(mechanics.axes[2]),
                actual_midpoint_displacement_m: [
                    self.dt_s * midpoint.dot(mechanics.axes[0]),
                    self.dt_s * midpoint.dot(mechanics.axes[1]),
                ],
                actual_midpoint_spin_angle_rad: self.dt_s * midpoint_spin,
                force_impulse_world_ns: impulse,
                spin_angular_impulse_world_nms: angular_impulse,
                midpoint_impulse_work_j: midpoint_work,
                dissipated_energy_j: heat,
                contact_work_residual_j: midpoint_work + contact.contact_work_j,
                contact_heat_residual_j: heat - contact.dissipated_energy_j,
                native_force_residual_n: force_residual,
                native_spin_moment_residual_nm: moment_residual,
                endpoint_velocity_residual_m_s: endpoint_velocity_residual,
                endpoint_spin_residual_rad_s: endpoint_spin_residual,
                maximum_assisting_work_j: assisting_work,
            };
            let controls = self.controls.mechanical;
            if [
                midpoint_work,
                heat,
                force_residual,
                moment_residual,
                endpoint_velocity_residual,
                endpoint_spin_residual,
                assisting_work,
                report.normal_velocity_after_m_s,
            ]
            .into_iter()
            .any(|x| !x.is_finite())
                || [
                    report.contact_work_residual_j,
                    report.contact_heat_residual_j,
                ]
                .into_iter()
                .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j)
                || force_residual > controls.contact.force_tolerance_n
                || moment_residual > controls.contact.moment_tolerance_nm
                || endpoint_velocity_residual > controls.velocity_tolerance_m_s
                || endpoint_spin_residual > controls.angular_velocity_tolerance_rad_s
                || (patch.normal_force_n > 0.0
                    && report.normal_velocity_after_m_s.abs() > controls.velocity_tolerance_m_s)
            {
                return Err(RigidHertzIslandError::ContactAudit(Box::new(report)));
            }
            if assisting_work > controls.energy_tolerance_j {
                return Err(RigidHertzIslandError::UnresolvedSlipEvent {
                    patch_id: patch.patch_id.clone(),
                    assisting_work_j: assisting_work,
                });
            }
            reports.push(report);
        }
        let (energy_after, linear_after, angular_after) = totals(&after, &self.properties, gate)?;
        let (energy_before, linear_before, angular_before) = self.totals_before;
        let midpoint_work = reports
            .iter()
            .map(|r| r.midpoint_impulse_work_j)
            .sum::<f64>();
        let heat = reports.iter().map(|r| r.dissipated_energy_j).sum::<f64>();
        let report = RigidHertzIslandReport {
            contacts: reports,
            native,
            kinetic_before_j: energy_before,
            kinetic_after_j: energy_after,
            midpoint_impulse_work_j: midpoint_work,
            dissipated_energy_j: heat,
            work_energy_residual_j: energy_after - energy_before - midpoint_work,
            total_energy_residual_j: energy_after - energy_before + heat,
            linear_momentum_change_ns: linear_after.sub(linear_before),
            angular_momentum_change_nms: angular_after.sub(angular_before),
        };
        let cfg = self.controls.mechanical;
        if [
            energy_before,
            energy_after,
            midpoint_work,
            heat,
            magnitude(report.linear_momentum_change_ns),
            magnitude(report.angular_momentum_change_nms),
        ]
        .into_iter()
        .any(|x| !x.is_finite())
            || [
                report.work_energy_residual_j,
                report.total_energy_residual_j,
                heat - report.native.dissipated_energy_j,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x.abs() > cfg.energy_tolerance_j)
            || magnitude(report.linear_momentum_change_ns) > cfg.linear_momentum_tolerance_ns
            || magnitude(report.angular_momentum_change_nms) > cfg.angular_momentum_tolerance_nms
        {
            return Err(RigidHertzIslandError::Audit(Box::new(report)));
        }
        poll(gate)?;
        Ok(RigidHertzIslandCandidate {
            body_indices: self.body_indices.clone(),
            bodies_after: after,
            report,
        })
    }
}

/// Publish every participating endpoint only after the complete joint audit.
pub fn apply_rigid_hertz_island(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidHertzIslandContact<'_>],
    dt_s: f64,
    controls: RigidHertzIslandControls,
    gate: &CancelGate,
) -> Result<RigidHertzIslandReport, RigidHertzIslandError> {
    let prepared =
        PreparedRigidHertzIsland::new(bodies, properties, contacts, dt_s, controls, gate)?;
    let laws: Vec<_> = contacts
        .iter()
        .map(|c| RigidHertzIslandLaw {
            patch_id: c.patch_id,
            law: c.contact.law,
        })
        .collect();
    let candidate = prepared.advance(&laws, gate)?;
    poll(gate)?;
    // All indices were validated during preparation; nothing fallible follows.
    for (&i, &state) in candidate.body_indices.iter().zip(&candidate.bodies_after) {
        bodies[i] = state;
    }
    Ok(candidate.report)
}
