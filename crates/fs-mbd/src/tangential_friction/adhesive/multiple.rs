//! Source-owned adhesive junctions on one actual shared-body island.
//!
//! Canonical patch IDs fix assembly and reduction order. The original normal
//! response supplies each true attached disk; its original interface supplies
//! both junction strengths. Complete inverse-inertia Gram rows and simultaneous
//! impulses preserve all cross-contact recoil, including dependent coordinates.
//! No normal load or Hertz pressure profile is used to invent adhesive capacity.
//!
//! This is a frozen-normal, frozen-geometry impulse phase. Every attached patch
//! requires closed normal motion at both native endpoints, including at zero
//! or tensile net load. Detached zero-area sources remain identified but create
//! no traction, matrix rows or heat. An entirely detached group is an exact
//! no-op. Physical stopping events, normal evolution, poses and heat allocation
//! remain with the surrounding owners. Refusal publishes no body or source.

use super::RigidAdhesiveContact;
use crate::tangential_friction::hertz::multiple::RigidHertzIslandError;
use crate::tangential_friction::hertz::{
    RigidHertzError, accepted_phase, finite, points, relative_point, relative_spin,
};
use crate::tangential_friction::island::{self, PatchMechanics};
use crate::tangential_friction::{magnitude, poll, totals, vector};
use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_contact::interface_binding::adhesive_shear::spinning::temperature::SourceThermalAdhesivePatch;
use fs_exec::CancelGate;
use fs_tribo::adhesive::spinning::contact::multiple::{
    AdhesiveContactInput, MultipleAdhesiveError, MultipleAdhesiveStep, PreparedAdhesiveContacts,
};
use fs_tribo::adhesive::spinning::contact::{
    AdhesiveWrenchContact, AdhesiveWrenchError, AdhesiveWrenchParameters,
};

/// Law-neutral actual endpoint, wrench, work and local acceptance receipt.
pub use crate::tangential_friction::hertz::multiple::RigidHertzIslandContactReport as RigidAdhesiveIslandContactReport;
/// Numerical body/iteration/energy controls; no normal-load assumption is stored.
pub use crate::tangential_friction::hertz::multiple::RigidHertzIslandControls as RigidAdhesiveIslandControls;
/// Law-neutral common native-body energy/momentum and full-response receipt.
pub use crate::tangential_friction::hertz::multiple::RigidHertzIslandReport as RigidAdhesiveIslandMechanicsReport;

/// Stable physical identity and original borrowed normal/shear source patch.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveIslandContact<'a> {
    /// Unique nonempty identity; caller order does not affect the physical solve.
    pub patch_id: &'a str,
    /// Actual ordered body points, tangent frame and source-owned junction.
    pub contact: RigidAdhesiveContact<'a>,
}

/// Temperature-refreshed junction strengths for one unchanged original disk.
/// The stable patch ID and original source identities are checked before solving.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveIslandThermalLaw<'s, 'a> {
    /// Original attached mechanical patch identity.
    pub patch_id: &'s str,
    /// Only the original junction-strength claims have been refreshed.
    pub source: &'s SourceThermalAdhesivePatch<'a>,
}

/// Original source and geometry retained in canonical physical ID order.
#[derive(Clone, Debug)]
pub struct RigidAdhesiveIslandPatch<'a> {
    /// Original stable physical identity.
    pub patch_id: String,
    /// Original body binding, frame and complete borrowed interface sources.
    pub contact: RigidAdhesiveContact<'a>,
    /// Actual points, never moved to manufacture a common contact location.
    pub contact_points_world: [Vec3; 2],
}

/// Attached-junction work and explicit zero-area exclusions for one candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveIslandReport {
    /// Native totals count each body participating in attached patches once.
    /// If all sources are detached, contacts/native rows/totals are empty/zero.
    pub mechanical: RigidAdhesiveIslandMechanicsReport,
    /// Canonical identities whose original solid contact area is exactly zero.
    pub detached_patch_ids: Vec<String>,
}

/// Immutable proposal for every participating original body.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveIslandCandidate {
    /// Sorted original body indices matching `bodies_after`.
    pub body_indices: Vec<usize>,
    /// All endpoints after summed simultaneous physical impulses.
    pub bodies_after: Vec<RigidBodyState>,
    /// Actual contact and body energy/momentum accounts.
    pub report: RigidAdhesiveIslandReport,
}

/// Reusable frozen mechanics and unchanged original source-owned junctions.
#[derive(Clone, Debug)]
pub struct PreparedRigidAdhesiveIsland<'a> {
    patches: Vec<RigidAdhesiveIslandPatch<'a>>,
    detached: Vec<RigidAdhesiveIslandPatch<'a>>,
    laws: Vec<AdhesiveWrenchContact>,
    mechanics: Vec<PatchMechanics>,
    body_indices: Vec<usize>,
    before: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    totals_before: (f64, Vec3, Vec3),
    dt_s: f64,
    controls: RigidAdhesiveIslandControls,
    mobility: Vec<f64>,
    native: Option<PreparedAdhesiveContacts>,
}

/// Source geometry, native mechanics and original-junction failures are distinct.
#[derive(Debug)]
pub enum RigidAdhesiveIslandError {
    /// Invalid original identity, geometry or dimensional controls.
    Invalid(&'static str),
    /// Original rigid-body kinematics or numerical finite-arithmetic refusal.
    Mechanics(RigidHertzError),
    /// Original uniform junction law admission failed.
    Law(AdhesiveWrenchError),
    /// Simultaneous original supports did not produce an accepted candidate.
    Native(MultipleAdhesiveError),
    /// A held wrench assists an actual endpoint; refine the physical interval.
    UnresolvedSlipEvent {
        /// Physical junction whose event has not been resolved.
        patch_id: String,
        /// Positive assisting endpoint work [J].
        assisting_work_j: f64,
    },
    /// Actual body motion failed a junction's original physical closure gates.
    ContactAudit(Box<RigidAdhesiveIslandContactReport>),
    /// Actual shared endpoints failed whole-island energy or momentum closure.
    Audit(Box<RigidAdhesiveIslandMechanicsReport>),
}
impl From<RigidHertzError> for RigidAdhesiveIslandError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidAdhesiveIslandError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<crate::tangential_friction::RigidFrictionError> for RigidAdhesiveIslandError {
    fn from(error: crate::tangential_friction::RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<AdhesiveWrenchError> for RigidAdhesiveIslandError {
    fn from(error: AdhesiveWrenchError) -> Self {
        Self::Law(error)
    }
}
impl From<MultipleAdhesiveError> for RigidAdhesiveIslandError {
    fn from(error: MultipleAdhesiveError) -> Self {
        Self::Native(error)
    }
}
impl From<RigidHertzIslandError> for RigidAdhesiveIslandError {
    fn from(error: RigidHertzIslandError) -> Self {
        match error {
            RigidHertzIslandError::Invalid(field) => Self::Invalid(field),
            RigidHertzIslandError::Hertz(error) => Self::Mechanics(error),
            RigidHertzIslandError::Native(error) => Self::Native(error),
            RigidHertzIslandError::UnresolvedSlipEvent {
                patch_id,
                assisting_work_j,
            } => Self::UnresolvedSlipEvent {
                patch_id,
                assisting_work_j,
            },
            RigidHertzIslandError::ContactAudit(report) => Self::ContactAudit(report),
            RigidHertzIslandError::Audit(report) => Self::Audit(report),
        }
    }
}
impl core::fmt::Display for RigidAdhesiveIslandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid adhesive island: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveIslandError {}

impl<'a> PreparedRigidAdhesiveIsland<'a> {
    /// Capture real attached contacts once, retaining every original normal and
    /// shear source. Spectators and detached-only bodies do not enter the active
    /// energy sum. Allocation/work bounds apply before native factor assembly.
    #[allow(clippy::too_many_lines)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidAdhesiveIslandContact<'a>],
        dt_s: f64,
        controls: RigidAdhesiveIslandControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidAdhesiveIslandError> {
        poll(gate)?;
        controls.validate()?;
        if contacts.is_empty()
            || contacts.len() > controls.maximum_contacts
            || bodies.len() != properties.len()
            || !dt_s.is_finite()
            || dt_s <= 0.0
        {
            return Err(RigidAdhesiveIslandError::Invalid(
                "matching native bodies, positive interval and bounded nonempty sources required",
            ));
        }
        let mut order: Vec<_> = (0..contacts.len()).collect();
        order.sort_by(|&a, &b| contacts[a].patch_id.cmp(contacts[b].patch_id));
        if order.iter().any(|&i| contacts[i].patch_id.is_empty())
            || order
                .windows(2)
                .any(|p| contacts[p[0]].patch_id == contacts[p[1]].patch_id)
        {
            return Err(RigidAdhesiveIslandError::Invalid(
                "unique nonempty patch identities required",
            ));
        }
        let mut body_indices = Vec::with_capacity(contacts.len() * 2);
        let mut attached_count = 0;
        for &i in &order {
            poll(gate)?;
            let contact = contacts[i].contact;
            let binding = contact.binding;
            if binding.body_a >= bodies.len()
                || binding.body_b >= bodies.len()
                || binding.body_a == binding.body_b
            {
                return Err(RigidAdhesiveIslandError::Invalid(
                    "two distinct valid native body indices required",
                ));
            }
            if contact.source.patch().parameters().contact_radius_m > 0.0 {
                attached_count += 1;
                body_indices.extend([binding.body_a, binding.body_b]);
            }
        }
        if (3 * attached_count) * (3 * attached_count) > controls.maximum_matrix_terms {
            return Err(RigidAdhesiveIslandError::Invalid(
                "attached response exceeds matrix budget",
            ));
        }
        body_indices.sort_unstable();
        body_indices.dedup();
        if body_indices.len() > controls.mechanical.maximum_bodies {
            return Err(RigidAdhesiveIslandError::Invalid(
                "participating body budget exceeded",
            ));
        }
        let before: Vec<_> = body_indices.iter().map(|&i| bodies[i]).collect();
        let captured_properties: Vec<_> = body_indices.iter().map(|&i| properties[i]).collect();
        let mut patches = Vec::with_capacity(attached_count);
        let mut detached = Vec::with_capacity(contacts.len() - attached_count);
        let mut laws = Vec::with_capacity(attached_count);
        let mut mechanics = Vec::with_capacity(attached_count);
        let mut factors = Vec::with_capacity(attached_count);
        for &index in &order {
            poll(gate)?;
            let identified = contacts[index];
            let contact = identified.contact;
            let binding = contact.binding;
            let pair = [bodies[binding.body_a], bodies[binding.body_b]];
            let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
            let arms = [binding.arm_a_body, binding.arm_b_body];
            let axes = [
                vector(contact.frame.longitudinal()),
                vector(contact.frame.lateral()),
                vector(contact.frame.normal()),
            ];
            let original = points(pair, pair_properties, arms)?;
            let radius = contact.source.patch().parameters().contact_radius_m;
            let patch = RigidAdhesiveIslandPatch {
                patch_id: identified.patch_id.to_owned(),
                contact,
                contact_points_world: original.map(|point| point.point_world),
            };
            if radius == 0.0 {
                detached.push(patch);
                continue;
            }
            if magnitude(original[0].point_world.sub(original[1].point_world))
                > controls.mechanical.geometry_tolerance_m
            {
                return Err(RigidAdhesiveIslandError::Invalid(
                    "attached original contact points must coincide",
                ));
            }
            let relative = relative_point(original);
            if relative.dot(axes[2]).abs() > controls.mechanical.velocity_tolerance_m_s {
                return Err(RigidAdhesiveIslandError::Invalid(
                    "held adhesive patch requires closed original normal motion",
                ));
            }
            let patch_mechanics = PatchMechanics {
                bodies: [
                    body_indices
                        .binary_search(&binding.body_a)
                        .expect("captured A"),
                    body_indices
                        .binary_search(&binding.body_b)
                        .expect("captured B"),
                ],
                arms,
                axes,
                before: original,
                rates_before: [
                    relative.dot(axes[0]),
                    relative.dot(axes[1]),
                    radius * relative_spin(original).dot(axes[2]),
                ],
            };
            factors.push(island::factor(
                pair,
                pair_properties,
                &patch_mechanics,
                radius,
            ));
            mechanics.push(patch_mechanics);
            laws.push(AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
                contact_radius_m: radius,
                static_shear_strength_pa: contact
                    .source
                    .static_patch()
                    .parameters()
                    .shear_strength_pa,
                kinetic_shear_strength_pa: contact.source.patch().parameters().shear_strength_pa,
            })?);
            patches.push(patch);
        }
        let mobility = island::mobility(&mechanics, &factors, gate)?;
        let native = if attached_count == 0 {
            None
        } else {
            let factor = finite(0.5 * dt_s * dt_s)?;
            if factor == 0.0 {
                return Err(RigidAdhesiveIslandError::Invalid(
                    "unrepresentable midpoint duration",
                ));
            }
            let compliance = mobility.iter().map(|&w| w * factor).collect();
            Some(PreparedAdhesiveContacts::new(
                attached_count,
                compliance,
                controls.constitutive(),
                || gate.is_requested(),
            )?)
        };
        let totals_before = totals(&before, &captured_properties, gate)?;
        poll(gate)?;
        Ok(Self {
            patches,
            detached,
            laws,
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

    /// Original attached source patches in canonical identity order.
    pub fn contacts(&self) -> &[RigidAdhesiveIslandPatch<'a>] {
        &self.patches
    }
    /// Original detached source patches, retained without mechanical constraints.
    pub fn detached_contacts(&self) -> &[RigidAdhesiveIslandPatch<'a>] {
        &self.detached
    }
    /// Canonical attached identities matching the native matrix block order.
    pub fn patch_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.patches.iter().map(|p| p.patch_id.as_str())
    }
    /// Original participating body indices, excluding detached-only bodies.
    pub fn body_indices(&self) -> &[usize] {
        &self.body_indices
    }
    /// Original physical duration; no sweep advances this clock.
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }
    /// Complete original shared-body inverse-mass Gram matrix [1/kg].
    pub fn scaled_mobility_per_kg(&self) -> &[f64] {
        &self.mobility
    }
    /// Complete original midpoint compliance [m/N], or empty if all detached.
    pub fn scaled_compliance_m_per_n(&self) -> &[f64] {
        self.native
            .as_ref()
            .map_or(&[], |native| native.scaled_compliance_m_per_n())
    }

    /// Propose simultaneous original-source impulses from the same captured
    /// endpoints. Retrying never accumulates a prior candidate or hidden state.
    pub fn advance(
        &self,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveIslandCandidate, RigidAdhesiveIslandError> {
        self.advance_laws(&self.laws, gate)
    }

    /// Reuse the captured full body response with the original junction claims
    /// evaluated at explicitly supplied interface temperatures. Normal source,
    /// disk geometry, body binding and accepted velocities remain unchanged.
    /// Source order is immaterial; exactly one successor per attached patch is
    /// required. This does not evolve normal adhesion or allocate friction heat.
    pub fn advance_thermal(
        &self,
        sources: &[RigidAdhesiveIslandThermalLaw<'_, '_>],
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveIslandCandidate, RigidAdhesiveIslandError> {
        poll(gate)?;
        if sources.len() != self.patches.len() {
            return Err(RigidAdhesiveIslandError::Invalid(
                "one thermal junction source per original attached patch required",
            ));
        }
        let mut ordered: Vec<_> = sources.iter().collect();
        ordered.sort_by_key(|source| source.patch_id);
        let mut laws = Vec::with_capacity(ordered.len());
        for (source, patch) in ordered.iter().zip(&self.patches) {
            poll(gate)?;
            let original = source.source.original();
            let expected = patch.contact.source;
            if source.patch_id != patch.patch_id
                || original.normal_state().identity() != expected.normal_state().identity()
                || original.shear_source().resolved().identity()
                    != expected.shear_source().resolved().identity()
                || original.contact_area_m2() != expected.contact_area_m2()
                || original.patch().parameters().contact_radius_m
                    != expected.patch().parameters().contact_radius_m
            {
                return Err(RigidAdhesiveIslandError::Invalid(
                    "thermal junction changed its original patch, source or disk geometry",
                ));
            }
            laws.push(*source.source.law());
        }
        self.advance_laws(&laws, gate)
    }

    fn advance_laws(
        &self,
        laws: &[AdhesiveWrenchContact],
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveIslandCandidate, RigidAdhesiveIslandError> {
        poll(gate)?;
        let native = if let Some(prepared) = &self.native {
            let inputs: Vec<_> = laws
                .iter()
                .zip(&self.mechanics)
                .map(|(law, mechanics)| AdhesiveContactInput {
                    law,
                    free_displacement_m: [
                        self.dt_s * mechanics.rates_before[0],
                        self.dt_s * mechanics.rates_before[1],
                    ],
                    free_spin_angle_rad: self.dt_s * mechanics.rates_before[2]
                        / law.parameters().contact_radius_m,
                })
                .collect();
            let phases: Vec<_> = self
                .mechanics
                .iter()
                .map(|m| accepted_phase(&m.rates_before))
                .collect();
            prepared.advance_phases(&inputs, &phases, || gate.is_requested())?
        } else {
            MultipleAdhesiveStep {
                contacts: Vec::new(),
                sweeps: 0,
                evaluations: 0,
                maximum_compatibility_residual_m: 0.0,
                dissipated_energy_j: 0.0,
                contact_work_j: 0.0,
                work_residual_j: 0.0,
            }
        };
        let after = island::apply_wrenches(
            &self.before,
            &self.mechanics,
            &native.contacts,
            self.dt_s,
            gate,
        )?;
        let mut reports = Vec::with_capacity(self.patches.len());
        for (((patch, mechanics), contact), law) in self
            .patches
            .iter()
            .zip(&self.mechanics)
            .zip(&native.contacts)
            .zip(laws)
        {
            poll(gate)?;
            let original = law.kinetic_patch();
            let report = island::contact_report(
                &patch.patch_id,
                mechanics,
                &after,
                &self.properties,
                original.parameters().contact_radius_m,
                self.dt_s,
                contact,
                |slip, spin| {
                    let response = original.evaluate(slip, spin)?;
                    Ok(island::NativeTraction {
                        force_n: response.force_n,
                        spin_moment_nm: response.spin_moment_nm,
                        dissipated_power_w: response.dissipated_power_w,
                    })
                },
            )?;
            if !island::contact_admitted(&report, self.controls.mechanical, true) {
                return Err(RigidAdhesiveIslandError::ContactAudit(Box::new(report)));
            }
            if report.maximum_assisting_work_j > self.controls.mechanical.energy_tolerance_j {
                return Err(RigidAdhesiveIslandError::UnresolvedSlipEvent {
                    patch_id: patch.patch_id.clone(),
                    assisting_work_j: report.maximum_assisting_work_j,
                });
            }
            reports.push(report);
        }
        let mechanical = island::report(
            reports,
            native,
            self.totals_before,
            &after,
            &self.properties,
            gate,
        )?;
        if !island::report_admitted(&mechanical, self.controls.mechanical) {
            return Err(RigidAdhesiveIslandError::Audit(Box::new(mechanical)));
        }
        let report = RigidAdhesiveIslandReport {
            mechanical,
            detached_patch_ids: self.detached.iter().map(|p| p.patch_id.clone()).collect(),
        };
        poll(gate)?;
        Ok(RigidAdhesiveIslandCandidate {
            body_indices: self.body_indices.clone(),
            bodies_after: after,
            report,
        })
    }
}

/// Publish each participating body's summed endpoint once, after every native
/// junction, whole-island and cancellation check. Detached inputs are no-ops.
pub fn apply_rigid_adhesive_island(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidAdhesiveIslandContact<'_>],
    dt_s: f64,
    controls: RigidAdhesiveIslandControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveIslandReport, RigidAdhesiveIslandError> {
    let prepared =
        PreparedRigidAdhesiveIsland::new(bodies, properties, contacts, dt_s, controls, gate)?;
    let candidate = prepared.advance(gate)?;
    poll(gate)?;
    for (&index, &body) in candidate.body_indices.iter().zip(&candidate.bodies_after) {
        bodies[index] = body;
    }
    Ok(candidate.report)
}
