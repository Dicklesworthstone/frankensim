//! Original source histories along the shared, admitted mixed-contact trajectory.
//!
//! The selected material director convects an orthonormal contact frame. Native
//! slip and texture components follow that director by an isometry, not by an
//! affine stretch of a friction ellipse. Each interval samples source claims at
//! its actual post-force incoming speed, solves ONE simultaneous mechanical
//! event, and advances finite travel only after the complete path is admitted.
//! Rebound and unloaded encounters retain their instantaneous source sample but
//! accumulate no free-flight slip. Sustained travel uses the explicitly declared
//! first-order held post-contact slip quadrature, with conservative nearness of
//! the complete affine geometric path required in addition to nonpenetration.
//!
//! History belongs to stable original scene features. New features initialize
//! from their original source law; an old feature retains its own history across
//! separation/re-entry. No history is transferred to an unrelated adjacent face.
//! Every retry and every whole-call refusal discards all private source changes.

use super::{
    FiniteSolidSustainedAcceptedContext, FiniteSolidSustainedRefusal,
    FiniteSolidSustainedTrajectoryError, FiniteSolidSustainedTrajectoryReport,
};
use crate::finite_solid::contact::mixed::{
    FiniteSolidContactCounterpart, FiniteSolidMixedContact, FiniteSolidMixedContactKinematics,
    SourceFiniteSolidMixedContactBinding, SourceFiniteSolidMixedContactError,
    SourceFiniteSolidMixedContactReceipt,
};
use crate::finite_solid::contact::scene::{FiniteSolidSphereScene, FiniteSolidSphereSceneKey};
use crate::finite_solid::contact::step::FiniteSolidSphereTrajectoryControls;
use crate::finite_solid::contact::surface::FiniteSolidContactSurface;
use crate::finite_solid::contact::texture::{
    FiniteSolidContactTexture, FiniteSolidContactTextureError,
};
use crate::finite_solid::{FiniteSolidBody, FiniteSolidMaterial};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_contact::frictional_multicontact::elliptic::EllipticFrictionalContactResult;
use fs_contact::interface_binding::thermal_friction::elliptic::history::{
    SourcePostcontactMotionControls, SourcePostcontactMotionError, SourcePostcontactMotionReport,
};
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_exec::Cx;
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Wrench};
use fs_query::TriangleQueryError;
use fs_solid::hyper3::prepared::PreparedHyperTetProblem;
use fs_tribo::TangentialLaw;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::transport::{
    FramedTangentialState, TangentialHistoryRotation, TangentialTransportError,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// An explicit ordered source law, coefficient role and participating material.
/// The material director is the orthonormal frame's physical transport policy;
/// it is not a relative torsional-slip input or an average of both bodies' spin.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidSourceContactBinding<'law> {
    pub source: SourceFiniteSolidMixedContactBinding<'law>,
    pub texture: FiniteSolidContactTexture,
}

/// Independent bounded history storage and original source length gates.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidSourceTrajectoryControls {
    pub maximum_history_entries: usize,
    pub motion: SourcePostcontactMotionControls,
}
impl Default for FiniteSolidSourceTrajectoryControls {
    fn default() -> Self {
        Self {
            maximum_history_entries: 4096,
            motion: SourcePostcontactMotionControls {
                absolute_slip_tolerance_m: 1e-12,
                absolute_path_tolerance_m: 1e-12,
            },
        }
    }
}

/// One accepted feature's native source state, in its selected material frame.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteSolidSourceContactHistory {
    pub contact_id: u64,
    pub key: FiniteSolidSphereSceneKey,
    pub texture: FiniteSolidContactTexture,
    pub coefficient_role: SourceImpulseCoefficientRole,
    pub frame: TangentFrame,
    pub state: SourceThermalEllipticState,
    /// Last represented mechanical time at which this feature was observed.
    pub last_contact_time_s: f64,
}

/// A scene clock and bounded, feature-keyed original histories.
///
/// Checkpoints share the immutable scene, index and source entries. Free-flight
/// retries clone only Arc handles; contact changes copy the index once and only
/// active entries clone their original source states.
/// The native solid preparation and exact prepared scene are pinned. A changed
/// namespace or stale clock requires explicit rebinding, never a silent reset.
#[derive(Clone, Debug)]
pub struct FiniteSolidSourceTrajectoryHistory<'a> {
    prepared: Arc<PreparedHyperTetProblem<'a, FiniteSolidMaterial>>,
    scene: Arc<FiniteSolidSphereScene>,
    time_s: f64,
    entries: Arc<BTreeMap<u64, Arc<FiniteSolidSourceContactHistory>>>,
}
impl<'a> FiniteSolidSourceTrajectoryHistory<'a> {
    pub fn new(solid: &FiniteSolidBody<'a>, scene: &FiniteSolidSphereScene) -> Self {
        Self {
            prepared: Arc::clone(&solid.prepared),
            scene: Arc::new(scene.clone()),
            time_s: solid.native_state().t,
            entries: Arc::new(BTreeMap::new()),
        }
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    pub fn get(&self, contact_id: u64) -> Option<&FiniteSolidSourceContactHistory> {
        self.entries.get(&contact_id).map(Arc::as_ref)
    }

    pub fn entries(&self) -> impl ExactSizeIterator<Item = &FiniteSolidSourceContactHistory> {
        self.entries.values().map(Arc::as_ref)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteSolidSourceConditionPhase {
    /// Actual incoming kinematics after this attempt's force kick.
    PreContact,
    /// Accepted held post-contact slip, with explicit conditions at the end time.
    AcceptedMotion,
}

/// Conditions remain physical caller inputs; speed and distance are source-owned.
/// For AcceptedMotion, the kinematics retain the phase's contact position/frame
/// and carry the actual post-contact velocity used by the held-slip quadrature.
pub struct FiniteSolidSourceConditionContext<'s> {
    pub phase: FiniteSolidSourceConditionPhase,
    pub key: FiniteSolidSphereSceneKey,
    pub time_s: f64,
    pub duration_s: f64,
    pub kinematics: FiniteSolidMixedContactKinematics,
    pub state: &'s SourceThermalEllipticState,
}

/// A finite native source-motion receipt or an explicitly instantaneous event.
#[derive(Clone, Debug)]
pub struct FiniteSolidSourceContactMotion {
    pub contact_id: u64,
    pub key: FiniteSolidSphereSceneKey,
    /// None for rebound/unloaded events: their free flight is not contact travel.
    pub motion: Option<SourcePostcontactMotionReport>,
    /// Certified upper clearance for every fraction of a sustained affine path.
    pub maximum_clearance_m: Option<f64>,
    /// Physical selected-director spin relative to minimum normal transport.
    /// This is not relative torsional slip and carries no work or heat.
    pub material_spin_angle_rad: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct FiniteSolidSourceTrajectoryEvent {
    pub sources: Vec<SourceFiniteSolidMixedContactReceipt>,
    pub motion: Vec<FiniteSolidSourceContactMotion>,
}

pub type FiniteSolidSustainedSourceTrajectoryReport = FiniteSolidSustainedTrajectoryReport<
    EllipticFrictionalContactResult,
    FiniteSolidSourceTrajectoryEvent,
>;
pub type FiniteSolidSustainedSourceTrajectoryError<E> = FiniteSolidSustainedTrajectoryError<
    EllipticFrictionalContactResult,
    FiniteSolidSustainedSourceFailure<E>,
    FiniteSolidSourceTrajectoryEvent,
>;

#[derive(Debug)]
pub enum FiniteSolidSustainedSourceFailure<E> {
    Invalid(&'static str),
    Binding {
        contact_id: u64,
        error: E,
    },
    Conditions {
        contact_id: u64,
        error: E,
    },
    Texture {
        contact_id: u64,
        error: FiniteSolidContactTextureError,
    },
    Transport {
        contact_id: u64,
        error: TangentialTransportError,
    },
    Motion {
        contact_id: u64,
        error: SourcePostcontactMotionError,
    },
    Geometry(TriangleQueryError),
    Contact(SourceFiniteSolidMixedContactError<E>),
    HistoryBudget {
        maximum_entries: usize,
    },
    Cancelled,
}

impl<'a> FiniteSolidBody<'a> {
    /// Advance source-dependent directional friction on the complete mixed scene.
    ///
    /// Binding/condition resolvers observe private candidates and may be called
    /// again after refinement. Keep their external effects read-only. Each
    /// stable feature must retain its original law, coefficient role and material
    /// attachment. The original source owner checks pinned law identity. Accepted
    /// histories, body states and native time return together after all gates.
    ///
    /// The selected material line defines a CONVECTED ORTHONORMAL director model:
    /// R/F transports the line, and native history follows the resulting tangent
    /// frame isometrically. Coefficients are current physical coefficients; solid
    /// stretch does not rescale them or stored slip. No common spin is inferred
    /// from relative torsion. Antipodal normal transfers remain unresolved.
    ///
    /// Source travel advances only for loaded nonopening contacts whose same
    /// endpoint feature and entire affine path remain within activation distance.
    /// An unresolved sustained interval refines the complete mechanical attempt.
    /// The held-slip quadrature is first order; it is not an exact curved-slip
    /// integral. Normal impulse supplies neither area/pressure nor heat partition.
    #[allow(
        clippy::too_many_arguments,
        clippy::too_many_lines,
        clippy::type_complexity
    )]
    pub fn sustained_source_elliptic_sphere_trajectory_candidate<'law, E>(
        &self,
        surface: &FiniteSolidContactSurface,
        scene: &FiniteSolidSphereScene,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FiniteSolidSphereTrajectoryControls,
        history: &FiniteSolidSourceTrajectoryHistory<'a>,
        source_controls: FiniteSolidSourceTrajectoryControls,
        binding: impl Fn(
            &FiniteSolidSphereSceneKey,
            &FiniteSolidMixedContact,
        ) -> Result<FiniteSolidSourceContactBinding<'law>, E>,
        conditions: impl Fn(
            &FiniteSolidSourceConditionContext<'_>,
        ) -> Result<SourceFrictionConditions, E>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            FiniteSolidSourceTrajectoryHistory<'a>,
            FiniteSolidSustainedSourceTrajectoryReport,
        ),
        FiniteSolidSustainedSourceTrajectoryError<E>,
    > {
        use FiniteSolidSustainedSourceFailure as F;
        if !Arc::ptr_eq(&self.prepared, &history.prepared)
            || history.scene.as_ref() != scene
            || history.time_s != self.native_state().t
            || source_controls.maximum_history_entries == 0
            || source_controls.maximum_history_entries > 65_536
            || history.entries.len() > source_controls.maximum_history_entries
            || [
                source_controls.motion.absolute_slip_tolerance_m,
                source_controls.motion.absolute_path_tolerance_m,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(FiniteSolidSustainedTrajectoryError::Invalid(
                "matching native scene/history clock and bounded finite source controls required",
            ));
        }
        self.sustained_sphere_trajectory_with(
            surface, scene, bodies, properties, nodal_forces_n, body_force_n_m3,
            wrenches, gravity, duration_s, origin_m, controls, history,
            |context, auxiliary| {
                let mut contacts = Vec::with_capacity(context.contacts.len());
                let mut histories = Vec::with_capacity(context.contacts.len());
                let mut bindings = Vec::with_capacity(context.contacts.len());
                for (&contact, &key) in context.contacts.iter().zip(context.contact_keys) {
                    poll(context.cx)?;
                    let resolved = binding(&key, &contact)
                        .map_err(|error| F::Binding { contact_id: contact.id, error })?;
                    let texture = context.solid.contact_texture_frame(
                        context.bodies, &contact, resolved.texture, context.cx,
                    ).map_err(|error| F::Texture { contact_id: contact.id, error })?;
                    let state = if let Some(previous) = auxiliary.entries.get(&contact.id) {
                        check_binding(previous, key, resolved)?;
                        material_transport(resolved.source.law, &previous.state,
                            previous.frame, texture.frame, context.cx)
                            .map_err(|error| F::Transport { contact_id: contact.id, error })?.0
                    } else {
                        if auxiliary.entries.len() == source_controls.maximum_history_entries {
                            return Err(F::HistoryBudget {
                                maximum_entries: source_controls.maximum_history_entries,
                            });
                        }
                        let state = resolved.source.law.initial_state();
                        FramedTangentialState::new(resolved.source.law, state, texture.frame)
                            .map_err(|error| F::Transport { contact_id: contact.id, error })?
                            .into_state()
                    };
                    Arc::make_mut(&mut auxiliary.entries).insert(contact.id, Arc::new(FiniteSolidSourceContactHistory {
                        contact_id: contact.id, key, texture: resolved.texture,
                        coefficient_role: resolved.source.coefficient_role,
                        frame: texture.frame, state: state.clone(),
                        last_contact_time_s: context.solid.native_state().t,
                    }));
                    contacts.push(texture.bind(contact).with_friction(()));
                    histories.push(state);
                    bindings.push(resolved.source);
                }
                let (solid, bodies, report) = context.solid.mixed_source_elliptic_contact_candidate(
                    context.bodies, context.properties, &contacts, &histories, &bindings,
                    context.origin_m, context.controls,
                    |state, kinematics| conditions(&FiniteSolidSourceConditionContext {
                        phase: FiniteSolidSourceConditionPhase::PreContact,
                        key: context.contact_keys[kinematics.original_index],
                        time_s: context.solid.native_state().t,
                        duration_s: context.duration_s, kinematics: *kinematics, state,
                    }), context.cx,
                ).map_err(F::Contact)?;
                for receipt in &report.sources {
                    let entry = Arc::make_mut(&mut auxiliary.entries).get_mut(&receipt.kinematics.contact_id)
                        .ok_or(F::Invalid("missing staged source feature"))?;
                    Arc::make_mut(entry).state = receipt.sample.state.clone();
                }
                Ok((solid, bodies, report.mechanical, FiniteSolidSourceTrajectoryEvent {
                    sources: report.sources, motion: Vec::new(),
                }))
            },
            |context, auxiliary, event| {
                let Some(event) = event else {
                    auxiliary.time_s = context.endpoint_solid.native_state().t;
                    return Ok(None);
                };
                let mechanical = context.impact.ok_or(F::Invalid("source event lacks mechanical receipt"))?;
                for receipt in &event.sources {
                    poll(context.cx)?;
                    let id = receipt.kinematics.contact_id;
                    let row = receipt.kinematics.original_index;
                    let contact = context.contacts[row];
                    let key = context.contact_keys[row];
                    let resolved = binding(&key, &contact)
                        .map_err(|error| F::Binding { contact_id: id, error })?;
                    let entry = Arc::make_mut(&mut auxiliary.entries).get_mut(&id)
                        .ok_or(F::Invalid("missing admitted source feature"))?;
                    check_binding(entry, key, resolved)?;
                    let entry = Arc::make_mut(entry);
                    // The original mechanical result is sorted by contact ID;
                    // source rows retain their original geometry-array index.
                    let solved_row = mechanical.contact_ids.binary_search(&id)
                        .map_err(|_| F::Invalid("source ID absent from original contact solve"))?;
                    let velocity = mechanical.contact.relative_velocities_after[solved_row];
                    let loaded = mechanical.contact.impulses[solved_row][0] > 0.0;
                    // A positive physical restitution target is rebound even
                    // when smaller than the numerical velocity tolerance.
                    let nonopening = matches!(contact.normal_law, NormalVelocityLaw::Target { speed } if speed == 0.0)
                        && velocity[0].abs() <= controls.impact.contact.velocity_tolerance;
                    let endpoint_row = context.endpoint_contacts.binary_search_by_key(&id, |c| c.id).ok();
                    let mut motion = None;
                    let mut maximum_clearance_m = None;
                    if loaded && nonopening {
                        let Some(end_row) = endpoint_row else {
                            return Ok(Some(FiniteSolidSustainedRefusal::AuxiliaryPath {
                                contact_id: id, what: "sustained source feature leaves the active endpoint set",
                            }));
                        };
                        if context.endpoint_contact_keys[end_row] != key {
                            return Err(F::Invalid("source feature ID rebound within an interval"));
                        }
                        let upper = nearness(surface, scene, &context, contact,
                            context.endpoint_contacts[end_row], key)
                            .map_err(F::Geometry)?;
                        if upper > contact.activation_tolerance_m {
                            return Ok(Some(FiniteSolidSustainedRefusal::AuxiliaryPath {
                                contact_id: id, what: "whole sustained source path nearness unresolved",
                            }));
                        }
                        maximum_clearance_m = Some(upper);
                        let mut kinematics = receipt.kinematics;
                        kinematics.relative_velocity_m_s = velocity;
                        let boundary = conditions(&FiniteSolidSourceConditionContext {
                            phase: FiniteSolidSourceConditionPhase::AcceptedMotion, key,
                            time_s: context.endpoint_solid.native_state().t,
                            duration_s: context.duration_s, kinematics, state: &entry.state,
                        }).map_err(|error| F::Conditions { contact_id: id, error })?;
                        let report = resolved.source.law.accept_postcontact_motion(
                            &entry.state, [velocity[1], velocity[2]], context.duration_s,
                            &boundary, source_controls.motion, &mut || context.cx.checkpoint().is_err(),
                        ).map_err(|error| F::Motion { contact_id: id, error })?;
                        entry.state = report.state.clone();
                        motion = Some(report);
                    }
                    let mut material_spin_angle_rad = None;
                    if let Some(end_row) = endpoint_row {
                        if context.endpoint_contact_keys[end_row] != key {
                            return Err(F::Invalid("source endpoint identity changed"));
                        }
                        let frame = context.endpoint_solid.contact_texture_frame(
                            context.endpoint_bodies, &context.endpoint_contacts[end_row],
                            entry.texture, context.cx,
                        ).map_err(|error| F::Texture { contact_id: id, error })?.frame;
                        let (state, spin) = material_transport(resolved.source.law, &entry.state,
                            entry.frame, frame, context.cx)
                            .map_err(|error| F::Transport { contact_id: id, error })?;
                        entry.state = state;
                        entry.frame = frame;
                        entry.last_contact_time_s = context.endpoint_solid.native_state().t;
                        material_spin_angle_rad = Some(spin);
                    }
                    event.motion.push(FiniteSolidSourceContactMotion {
                        contact_id: id, key, motion, maximum_clearance_m, material_spin_angle_rad,
                    });
                }
                auxiliary.time_s = context.endpoint_solid.native_state().t;
                poll(context.cx)?;
                Ok(None)
            }, cx,
        )
    }
}

fn check_binding<E>(
    previous: &FiniteSolidSourceContactHistory,
    key: FiniteSolidSphereSceneKey,
    binding: FiniteSolidSourceContactBinding<'_>,
) -> Result<(), FiniteSolidSustainedSourceFailure<E>> {
    if previous.key != key
        || previous.texture != binding.texture
        || previous.coefficient_role != binding.source.coefficient_role
    {
        return Err(FiniteSolidSustainedSourceFailure::Invalid(
            "retained feature must keep its original material attachment and source role",
        ));
    }
    // Even an unloaded/rebounding event without an endpoint feature must keep
    // its original source identity. Such an event has no finite-motion call
    // on which to rely for the native owner's identity admission.
    FramedTangentialState::new(binding.source.law, previous.state.clone(), previous.frame)
        .map_err(|error| FiniteSolidSustainedSourceFailure::Transport {
            contact_id: previous.contact_id,
            error,
        })?;
    Ok(())
}

fn material_transport(
    law: &SourceThermalEllipticCoulomb,
    state: &SourceThermalEllipticState,
    previous: TangentFrame,
    current: TangentFrame,
    cx: &Cx<'_>,
) -> Result<(SourceThermalEllipticState, f64), TangentialTransportError> {
    // The selected physical director supplies the twist relative to minimum
    // normal transport. Applying that twist with the original isometry keeps
    // components convected with the director, instead of preserving an old
    // world texture by mistaking the physical turn for a mere tangent gauge.
    let [cosine, sine] = TangentialHistoryRotation::between(previous, current, 0.0)?.cosine_sine();
    let spin = -fs_math::det::atan2(sine, cosine);
    if cx.checkpoint().is_err() {
        return Err(TangentialTransportError::Cancelled);
    }
    // For this declared director convention, minimum-normal rotation followed
    // by the physical director twist is EXACTLY the identity in convected
    // coordinates. Preserve components exactly instead of a sin/atan2 round
    // trip that would manufacture tiny slip changes at every rotating step.
    // The original law still admits identity transport and its source pin.
    let transported = FramedTangentialState::new(law, state.clone(), current)?;
    Ok((transported.into_state(), spin))
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), FiniteSolidSustainedSourceFailure<E>> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidSustainedSourceFailure::Cancelled)
}

fn nearness(
    surface: &FiniteSolidContactSurface,
    scene: &FiniteSolidSphereScene,
    context: &FiniteSolidSustainedAcceptedContext<'_, '_, '_, EllipticFrictionalContactResult>,
    contact: FiniteSolidMixedContact,
    endpoint: FiniteSolidMixedContact,
    key: FiniteSolidSphereSceneKey,
) -> Result<f64, TriangleQueryError> {
    use fs_query::triangle::sweep::{linear_sphere_triangle_nearness, linear_spheres_nearness};
    let sphere = |id| {
        scene
            .spheres()
            .binary_search_by_key(&id, |s| s.id)
            .map(|row| scene.spheres()[row])
            .map_err(|_| TriangleQueryError::Invalid("source sphere is outside prepared scene"))
    };
    let center = |bodies: &[RigidBodyState], body: usize| {
        let p = bodies[body].pose().position_world();
        [p.x, p.y, p.z]
    };
    match key {
        FiniteSolidSphereSceneKey::Spheres { first, second } => {
            let a = sphere(first)?;
            let b = sphere(second)?;
            linear_spheres_nearness(
                center(context.post_contact_bodies, a.rigid_body),
                center(context.endpoint_bodies, a.rigid_body),
                a.radius_m,
                center(context.post_contact_bodies, b.rigid_body),
                center(context.endpoint_bodies, b.rigid_body),
                b.radius_m,
                context.cx,
            )
        }
        FiniteSolidSphereSceneKey::Solid {
            sphere: id, face, ..
        } => {
            let s = sphere(id)?;
            let nodes = *surface
                .faces()
                .get(face)
                .ok_or(TriangleQueryError::Invalid(
                    "source face is outside prepared surface",
                ))?;
            let weights = |c: FiniteSolidMixedContact| {
                let FiniteSolidContactCounterpart::Solid(point) = c.counterpart else {
                    return Err(TriangleQueryError::Invalid(
                        "source solid feature has rigid counterpart",
                    ));
                };
                let tet = context
                    .post_contact_solid
                    .problem()
                    .tetrahedra
                    .get(point.tetrahedron())
                    .ok_or(TriangleQueryError::Invalid(
                        "source point is outside original mesh",
                    ))?;
                let mut weights = [0.0; 3];
                for (row, node) in nodes.into_iter().enumerate() {
                    let local =
                        tet.iter()
                            .position(|&n| n == node)
                            .ok_or(TriangleQueryError::Invalid(
                                "source witness is not on selected material face",
                            ))?;
                    weights[row] = point.weights()[local];
                }
                Ok(weights)
            };
            // Match the original native current-position arithmetic, gathering
            // only these three nodes. The proof encloses convex interpolation
            // from the ORIGINAL node floats, never a rounded material point.
            let node_position = |solid: &FiniteSolidBody<'_>, node: usize| {
                std::array::from_fn(|axis| {
                    let dof = 3 * node + axis;
                    let u = solid
                        .free
                        .binary_search(&dof)
                        .map_or(solid.fixed[dof], |row| solid.state.q[row]);
                    solid.problem().nodes_m[node][axis] + u
                })
            };
            linear_sphere_triangle_nearness(
                center(context.post_contact_bodies, s.rigid_body),
                center(context.endpoint_bodies, s.rigid_body),
                s.radius_m,
                nodes.map(|n| node_position(context.post_contact_solid, n)),
                nodes.map(|n| node_position(context.endpoint_solid, n)),
                weights(contact)?,
                weights(endpoint)?,
                context.cx,
            )
        }
    }
}
