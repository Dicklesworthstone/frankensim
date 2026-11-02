//! Automatic native-geometry walks for one differentiated optical branch.
//!
//! A launch and deterministic reflection/transmission policy are sufficient:
//! the prepared scene finds every actual boundary, including blockers and the
//! opposite wall of the same object. The central ray and its two original
//! launch-angle columns use one augmented boundary solve per segment. Native
//! normal derivatives and local media then supply the smooth interface map.
//!
//! Physical Stokes transfer, log interface power, Beer depth and phase/group
//! path lengths remain separate and accumulate exactly once. Contact collars,
//! grazing/critical refusals and work budgets retain the native event owner's
//! semantics. No finite-difference Hessian, rough-scattering derivative or
//! complete enumeration of optical images is inferred. A limited result is
//! diagnostic output, not a resumable token with a missing departure witness.

use fs_blake3::DomainHasher;
use fs_geom::Vec3;

use super::*;
use crate::charts::Ray;
use crate::dielectric::fresnel_dielectric;
use crate::gradient_index::GradientIndexSample;
use crate::gradient_index::bundle::event::BundleEventOptions;
use crate::gradient_index::bundle::interface::InterfaceOptions;
use crate::gradient_index::bundle::{BundleError, BundleOptions, propagate_bundle};
use crate::polarization::dielectric::transport::bundle::{
    PolarizedBundleError, PolarizedRayBundle,
};

/// Bounded local source-to-receiver connections through these native paths.
pub mod connection;

/// One central integration envelope and two additional derivative tolerances.
/// The `trace.event.integration` controls are authoritative for every segment.
#[derive(Clone, Copy, Debug)]
pub struct CurvedBundleSettings {
    pub trace: CurvedTraceSettings,
    pub position_derivative_tolerance_m: f64,
    pub direction_derivative_tolerance: f64,
    pub interface: InterfaceOptions,
}

impl Default for CurvedBundleSettings {
    fn default() -> Self {
        let bundle = BundleOptions::default();
        Self {
            trace: CurvedTraceSettings {
                max_surface_events: 64,
                max_distance_m: 100.0,
                max_total_steps: 65_536,
                max_query_work: 1_000_000,
                event: GradientEventOptions::default(),
            },
            position_derivative_tolerance_m: bundle.position_derivative_tolerance_m,
            direction_derivative_tolerance: bundle.direction_derivative_tolerance,
            interface: InterfaceOptions::default(),
        }
    }
}

impl CurvedBundleSettings {
    fn validate(self) -> Result<(), CurvedBundleSceneError> {
        let trace = self.trace;
        if trace.max_surface_events == 0
            || trace.max_surface_events > MAX_POLARIZED_PATH_EVENTS
            || !trace.max_distance_m.is_finite()
            || trace.max_distance_m <= 0.0
            || trace.max_total_steps == 0
            || trace.max_total_steps > MAX_GRADIENT_INDEX_STEPS
            || trace.max_query_work == 0
            || [
                self.position_derivative_tolerance_m,
                self.direction_derivative_tolerance,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value <= 0.0)
            || [
                self.interface.min_incidence_cosine,
                self.interface.min_transmitted_cosine,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value <= 0.0 || value > 1.0)
        {
            return Err(SceneError::InvalidSettings.into());
        }
        Ok(())
    }
}

/// Actual local optics available to a deterministic branch policy. A callback
/// is invoked only for an admitted smooth dielectric contact, after topology,
/// field support and the native surface differential have been checked.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedBundleBranchContext {
    pub event_index: u32,
    pub interaction: Interaction,
    pub ray: Ray,
    pub incident_sample: GradientIndexSample,
    pub target_sample: GradientIndexSample,
    pub total_internal_reflection: bool,
}

/// Allocation-free result of a native differentiated branch walk.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedBundlePathResult {
    bundle: PolarizedRayBundle,
    medium: usize,
    medium_identity: Option<ContentHash>,
    termination: PathTermination,
    diffuse_receiver: bool,
    surface_events: u32,
    last_contact: Option<CurvedContact>,
    detector_reference_axis: Option<Vec3>,
    max_contact_resolution_m: f64,
    boundary_evaluations: u64,
    query_work: u64,
    interface_field_evaluations: u32,
    branch_signature: ContentHash,
}

impl CurvedBundlePathResult {
    pub const fn bundle(&self) -> &PolarizedRayBundle {
        &self.bundle
    }
    pub const fn medium_index(&self) -> usize {
        self.medium
    }
    pub const fn medium_snapshot_identity(&self) -> Option<ContentHash> {
        self.medium_identity
    }
    pub const fn termination(&self) -> PathTermination {
        self.termination
    }
    /// The selected native diffuse receiver was reached with its incoming
    /// bundle intact. No reflectance, BSDF continuation or detector frame was
    /// applied. An ordinary opaque blocker never sets this marker.
    pub const fn is_diffuse_receiver(&self) -> bool {
        self.diffuse_receiver
    }
    pub const fn surface_events(&self) -> u32 {
        self.surface_events
    }
    pub const fn last_contact(&self) -> Option<&CurvedContact> {
        self.last_contact.as_ref()
    }
    /// Authored receiver axis; it does not replace the bundle's transported
    /// frame or its immutable two launch coordinates.
    pub const fn detector_reference_axis(&self) -> Option<Vec3> {
        self.detector_reference_axis
    }
    /// Largest declared contact collar, not accumulated derivative uncertainty.
    pub const fn max_contact_resolution_m(&self) -> f64 {
        self.max_contact_resolution_m
    }
    /// Cumulative boundary/departure oracle calls from all segment solves.
    pub const fn boundary_evaluations(&self) -> u64 {
        self.boundary_evaluations
    }
    /// Native geometry work includes the analytic normal-differential queries.
    pub const fn query_work(&self) -> u64 {
        self.query_work
    }
    /// Field queries include launch, all integration trials and target media.
    pub fn field_evaluations(&self) -> u64 {
        u64::from(self.bundle.propagation().propagation().field_evaluations)
            + u64::from(self.interface_field_evaluations)
    }
    /// Ordered native object/geometry/frame/feature, medium and selected-branch
    /// identity. Numerical points, directions and step counts are excluded so
    /// neighboring launches on the same smooth branch can be compared. This is
    /// a local branch guard within one immutable scene, not a topology proof.
    pub const fn branch_signature(&self) -> ContentHash {
        self.branch_signature
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CurvedBundleSceneError {
    Scene(CurvedSceneError),
    Bundle(PolarizedBundleError),
    MissingSecondOrderField {
        medium: usize,
    },
    UnsupportedInteraction {
        object_id: u64,
    },
    /// An exact prescribed branch sequence omitted this dielectric encounter.
    /// The encounter is counted independently of conductors and detectors.
    BranchPlanExhausted {
        dielectric_event: usize,
    },
    /// Tracing starts at an original, zero-distance launch, not a partial path
    /// whose geometry identity and interface departure have been discarded.
    InvalidInitialBundle,
    NumericalRange,
}
impl From<CurvedSceneError> for CurvedBundleSceneError {
    fn from(error: CurvedSceneError) -> Self {
        Self::Scene(error)
    }
}
impl From<SceneError> for CurvedBundleSceneError {
    fn from(error: SceneError) -> Self {
        Self::Scene(error.into())
    }
}
impl From<ProximityError> for CurvedBundleSceneError {
    fn from(error: ProximityError) -> Self {
        Self::Scene(error.into())
    }
}
impl From<GradientIndexError> for CurvedBundleSceneError {
    fn from(error: GradientIndexError) -> Self {
        Self::Scene(error.into())
    }
}
impl From<PolarizedBundleError> for CurvedBundleSceneError {
    fn from(error: PolarizedBundleError) -> Self {
        Self::Bundle(error)
    }
}
impl From<BundleError> for CurvedBundleSceneError {
    fn from(error: BundleError) -> Self {
        Self::Bundle(error.into())
    }
}
impl From<BeamError> for CurvedBundleSceneError {
    fn from(error: BeamError) -> Self {
        Self::Scene(error.into())
    }
}
impl From<fs_exec::Cancelled> for CurvedBundleSceneError {
    fn from(_: fs_exec::Cancelled) -> Self {
        SceneError::Cancelled.into()
    }
}
impl core::fmt::Display for CurvedBundleSceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved polarized bundle scene: {self:?}")
    }
}
impl core::error::Error for CurvedBundleSceneError {}

impl CurvedPolarizedScene<'_> {
    /// Construct a fresh source-angle bundle in the explicitly selected local
    /// medium. Absolute spectral intensity remains with the calibrated source;
    /// `source_polarization` has unit intensity in the supplied transverse frame.
    pub fn launch_bundle(
        &self,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        source_polarization: Stokes,
        initial_medium: usize,
    ) -> Result<PolarizedRayBundle, CurvedBundleSceneError> {
        let homogeneous;
        let field: &dyn SecondOrderIndexField = match self
            .media
            .get(initial_medium)
            .ok_or(SceneError::InvalidInitialMedium)?
        {
            SceneMedium::Homogeneous(material) => {
                homogeneous = HomogeneousIndexField::new(material);
                &homogeneous
            }
            SceneMedium::Differentiable { field, .. } => *field,
            SceneMedium::GradientIndex { .. } => {
                return Err(CurvedBundleSceneError::MissingSecondOrderField {
                    medium: initial_medium,
                });
            }
        };
        Ok(PolarizedRayBundle::try_new(
            propagate_bundle(
                field,
                ray,
                transverse_axis,
                wavelength_nm,
                0.0,
                BundleOptions::default(),
            )?,
            source_polarization,
        )?)
    }

    /// Follow transmission whenever a propagating transmitted ray exists;
    /// follow reflection at total internal reflection. Smooth metals always
    /// reflect. This selects one physical branch, not a sum over all branches.
    pub fn trace_bundle_transmitted(
        &self,
        cx: &Cx<'_>,
        bundle: PolarizedRayBundle,
        initial_medium: usize,
        settings: CurvedBundleSettings,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        self.trace_bundle_with_branches(cx, bundle, initial_medium, settings, |_, context| {
            if context.total_internal_reflection {
                DielectricEvent::Reflection
            } else {
                DielectricEvent::Transmission
            }
        })
    }

    /// Locate actual native interfaces and propagate the original launch
    /// derivatives through one selected smooth branch. The callback supplies
    /// only a branch, never a distance, surface, medium or forged contact.
    pub fn trace_bundle_with_branches(
        &self,
        cx: &Cx<'_>,
        bundle: PolarizedRayBundle,
        initial_medium: usize,
        settings: CurvedBundleSettings,
        mut branches: impl FnMut(&CurvedContact, &CurvedBundleBranchContext) -> DielectricEvent,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        self.trace_bundle_with_checked_branches(
            cx,
            bundle,
            initial_medium,
            settings,
            false,
            |contact, context| Ok(branches(contact, context)),
        )
    }

    // One event walk also serves exact branch plans. A policy refusal must stop
    // at the encounter, before an invented fallback branch is propagated.
    // Connections continue zero-power geometry to validate the receiver and
    // complete plan; ordinary forward transport retains early extinction.
    pub(super) fn trace_bundle_with_checked_branches(
        &self,
        cx: &Cx<'_>,
        bundle: PolarizedRayBundle,
        initial_medium: usize,
        settings: CurvedBundleSettings,
        continue_extinguished: bool,
        branches: impl FnMut(
            &CurvedContact,
            &CurvedBundleBranchContext,
        ) -> Result<DielectricEvent, CurvedBundleSceneError>,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        self.trace_bundle_with_checked_branches_to_receiver(
            cx,
            bundle,
            initial_medium,
            settings,
            continue_extinguished,
            None,
            branches,
        )
    }

    // The connector may receive incident power at one original diffuse surface.
    // All other contacts retain the same native event and interaction policy.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_bundle_with_checked_branches_to_receiver(
        &self,
        cx: &Cx<'_>,
        bundle: PolarizedRayBundle,
        initial_medium: usize,
        settings: CurvedBundleSettings,
        continue_extinguished: bool,
        diffuse_receiver: Option<u64>,
        mut branches: impl FnMut(
            &CurvedContact,
            &CurvedBundleBranchContext,
        ) -> Result<DielectricEvent, CurvedBundleSceneError>,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        cx.checkpoint()?;
        settings.validate()?;
        self.media
            .get(initial_medium)
            .ok_or(SceneError::InvalidInitialMedium)?;
        let initial = bundle.propagation();
        if initial.propagation().distance_m != 0.0
            || initial.interface_events() != 0
            || initial.initial_ray() != initial.propagation().ray
        {
            return Err(CurvedBundleSceneError::InvalidInitialBundle);
        }
        let mut signature = DomainHasher::new("org.frankensim.curved-bundle-branch.v1");
        signature.update(&(initial_medium as u64).to_le_bytes());
        append_medium(&mut signature, self.medium_identity(initial_medium));
        let mut result = CurvedBundlePathResult {
            bundle,
            medium: initial_medium,
            medium_identity: self.medium_identity(initial_medium),
            termination: PathTermination::SurfaceLimit,
            diffuse_receiver: false,
            surface_events: 0,
            last_contact: None,
            detector_reference_axis: None,
            max_contact_resolution_m: 0.0,
            boundary_evaluations: 0,
            query_work: 0,
            interface_field_evaluations: 0,
            branch_signature: signature.finalize(),
        };
        let mut oracle = self.boundaries.oracle(
            cx,
            settings.trace.max_query_work,
            settings.trace.event.surface_tolerance_m,
        );
        let mut departure = None;
        let mut remaining_steps = settings.trace.max_total_steps;
        for event_index in 0..settings.trace.max_surface_events {
            cx.checkpoint()?;
            let remaining_distance = settings.trace.max_distance_m
                - result.bundle.propagation().propagation().distance_m;
            if remaining_distance <= 0.0 {
                result.termination = PathTermination::DistanceLimit;
                break;
            }
            if remaining_steps == 0 {
                return Err(GradientIndexError::StepBudgetExceeded.into());
            }
            let mut boundary = settings.trace.event;
            boundary.integration.max_steps = boundary.integration.max_steps.min(remaining_steps);
            let options = BundleEventOptions {
                boundary,
                position_derivative_tolerance_m: settings.position_derivative_tolerance_m,
                direction_derivative_tolerance: settings.direction_derivative_tolerance,
            };
            let homogeneous;
            let field: &dyn SecondOrderIndexField = match &self.media[result.medium] {
                SceneMedium::Homogeneous(material) => {
                    homogeneous = HomogeneousIndexField::new(material);
                    &homogeneous
                }
                SceneMedium::Differentiable { field, .. } => *field,
                SceneMedium::GradientIndex { .. } => {
                    return Err(CurvedBundleSceneError::MissingSecondOrderField {
                        medium: result.medium,
                    });
                }
            };
            let (next, segment) = result
                .bundle
                .continue_to_boundary_with_control(
                    field,
                    &mut oracle,
                    remaining_distance,
                    departure.take(),
                    options,
                    || cx.checkpoint().is_ok(),
                )
                .map_err(|error| {
                    if let Some(geometry) = oracle.last_error() {
                        CurvedBundleSceneError::from(geometry)
                    } else if event_index == 0
                        && error == PolarizedBundleError::Bundle(BundleError::InitialSampleMismatch)
                    {
                        SceneError::InvalidInitialMedium.into()
                    } else {
                        error.into()
                    }
                })?;
            remaining_steps = remaining_steps
                .checked_sub(segment.segment.accepted_steps)
                .and_then(|value| value.checked_sub(segment.segment.rejected_steps))
                .ok_or(CurvedBundleSceneError::NumericalRange)?;
            result.bundle = next;
            result.boundary_evaluations = result
                .boundary_evaluations
                .checked_add(u64::from(segment.boundary_evaluations))
                .ok_or(CurvedBundleSceneError::NumericalRange)?;
            result.max_contact_resolution_m = result
                .max_contact_resolution_m
                .max(segment.contact_resolution_m);
            let BoundaryTermination::Surface(witness) = segment.termination else {
                result.termination = if self.straight_ray_escapes(
                    result.medium,
                    result.bundle.propagation().propagation().ray,
                    settings.trace.event.surface_tolerance_m,
                )? {
                    PathTermination::Escaped
                } else {
                    PathTermination::DistanceLimit
                };
                break;
            };
            let instance = self
                .geometry
                .instances()
                .get(witness.boundary_id)
                .ok_or(CurvedSceneError::InvalidContact)?;
            let binding = self
                .bindings
                .get(witness.boundary_id)
                .ok_or(CurvedSceneError::InvalidContact)?;
            if instance.object_id() != binding.object_id {
                return Err(CurvedSceneError::InvalidContact.into());
            }
            let contact = CurvedContact {
                object_id: instance.object_id(),
                geometry_identity: instance.geometry_identity(),
                frame_identity: instance.frame_identity(),
                witness,
                resolution_m: segment.contact_resolution_m,
            };
            result.last_contact = Some(contact);
            result.surface_events += 1;
            let arriving = result.bundle.propagation().propagation().ray;
            let incoming_medium = result.medium;
            let outgoing_tag = match binding.interaction {
                Interaction::Detector { reference_axis } => {
                    PolarizationFrame::try_new(arriving.dir, reference_axis)?;
                    result.detector_reference_axis = Some(reference_axis);
                    result.termination = PathTermination::Detected;
                    0
                }
                Interaction::Opaque => {
                    result.termination = PathTermination::Opaque;
                    1
                }
                Interaction::DepolarizingDiffuse { medium, .. }
                    if diffuse_receiver == Some(instance.object_id()) =>
                {
                    require_medium(instance.object_id(), result.medium, medium)?;
                    result.termination = PathTermination::Opaque;
                    result.diffuse_receiver = true;
                    5
                }
                Interaction::SmoothConductor { optics, medium } => {
                    require_medium(instance.object_id(), result.medium, medium)?;
                    let surface = oracle.normal_differential(witness, arriving.origin)?;
                    result.bundle = result.bundle.smooth_conductor_interface(
                        surface.differential,
                        &self.conductors[optics],
                        settings.interface,
                    )?;
                    2
                }
                Interaction::RoughConductor { .. }
                | Interaction::BrushedConductor { .. }
                | Interaction::DepolarizingDiffuse { .. } => {
                    return Err(CurvedBundleSceneError::UnsupportedInteraction {
                        object_id: instance.object_id(),
                    });
                }
                interaction @ (Interaction::Dielectric {
                    outside, inside, ..
                }
                | Interaction::CoatedDielectric {
                    outside, inside, ..
                }) => {
                    if matches!(interaction, Interaction::Dielectric { surface, .. } if !surface.is_delta())
                    {
                        return Err(CurvedBundleSceneError::UnsupportedInteraction {
                            object_id: instance.object_id(),
                        });
                    }
                    let cosine = witness.normal.dot(arriving.dir);
                    if !cosine.is_finite() || cosine == 0.0 {
                        return Err(CurvedSceneError::InvalidContact.into());
                    }
                    let (expected, target) = if cosine < 0.0 {
                        (outside, inside)
                    } else {
                        (inside, outside)
                    };
                    require_medium(instance.object_id(), result.medium, expected)?;
                    let wavelength = result.bundle.propagation().wavelength_nm();
                    let target_sample =
                        self.bundle_medium_sample(target, arriving.origin, wavelength)?;
                    result.interface_field_evaluations += 1;
                    let surface = oracle.normal_differential(witness, arriving.origin)?;
                    let incident_sample = result.bundle.propagation().propagation().final_sample;
                    let total_internal_reflection = fresnel_dielectric(
                        surface
                            .differential
                            .normal
                            .dot(arriving.dir)
                            .abs()
                            .clamp(0.0, 1.0),
                        incident_sample.phase_index(),
                        target_sample.phase_index(),
                    )
                    .map_err(BeamError::from)?
                    .total_internal_reflection;
                    let context = CurvedBundleBranchContext {
                        event_index,
                        interaction,
                        ray: arriving,
                        incident_sample,
                        target_sample,
                        total_internal_reflection,
                    };
                    let branch = branches(&contact, &context)?;
                    cx.checkpoint()?;
                    result.bundle = match interaction {
                        Interaction::Dielectric {
                            surface: material, ..
                        } => result.bundle.smooth_surface_interface(
                            surface.differential,
                            target_sample,
                            &material,
                            branch,
                            settings.interface,
                        )?,
                        Interaction::CoatedDielectric { coating, .. } => {
                            result.bundle.smooth_coated_interface(
                                surface.differential,
                                target_sample,
                                &self.coatings[coating],
                                branch,
                                settings.interface,
                            )?
                        }
                        _ => return Err(CurvedSceneError::InvalidContact.into()),
                    };
                    if branch == DielectricEvent::Transmission {
                        result.medium = target;
                        result.medium_identity = self.medium_identity(target);
                        4
                    } else {
                        3
                    }
                }
            };
            append_contact(
                &mut signature,
                contact,
                incoming_medium,
                result.medium,
                witness.normal.dot(arriving.dir) > 0.0,
                outgoing_tag,
            );
            append_medium(&mut signature, result.medium_identity);
            if matches!(
                result.termination,
                PathTermination::Detected | PathTermination::Opaque
            ) {
                break;
            }
            if !continue_extinguished && result.bundle.interface_log_power() == f64::NEG_INFINITY {
                result.termination = PathTermination::Extinguished;
                break;
            }
            let outgoing = witness
                .normal
                .dot(result.bundle.propagation().propagation().ray.dir);
            if !outgoing.is_finite() || outgoing == 0.0 {
                return Err(CurvedSceneError::InvalidContact.into());
            }
            departure = Some(GradientDeparture {
                witness,
                side: if outgoing > 0.0 {
                    BoundarySide::Positive
                } else {
                    BoundarySide::Negative
                },
            });
        }
        cx.checkpoint()?;
        result.query_work = oracle.work_used();
        result.branch_signature = signature.finalize();
        Ok(result)
    }

    fn bundle_medium_sample(
        &self,
        index: usize,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, CurvedBundleSceneError> {
        match self.media[index] {
            SceneMedium::Homogeneous(material) => {
                Ok(HomogeneousIndexField::new(material).sample(point, wavelength_nm)?)
            }
            SceneMedium::Differentiable { field, .. } => {
                Ok(field.sample_second_order(point, wavelength_nm)?.optical())
            }
            SceneMedium::GradientIndex { .. } => {
                Err(CurvedBundleSceneError::MissingSecondOrderField { medium: index })
            }
        }
    }
}

fn require_medium(
    object_id: u64,
    actual: usize,
    expected: usize,
) -> Result<(), CurvedBundleSceneError> {
    if actual != expected {
        return Err(SceneError::MediumMismatch {
            object_id,
            expected,
            actual,
        }
        .into());
    }
    Ok(())
}

fn append_medium(signature: &mut DomainHasher, identity: Option<ContentHash>) {
    match identity {
        None => signature.update(&[0]),
        Some(identity) => {
            signature.update(&[1]);
            signature.update(&identity.0);
        }
    }
}

fn append_contact(
    signature: &mut DomainHasher,
    contact: CurvedContact,
    incoming_medium: usize,
    outgoing_medium: usize,
    outward: bool,
    branch: u8,
) {
    signature.update(&contact.object_id.to_le_bytes());
    signature.update(&contact.geometry_identity.0);
    signature.update(&contact.frame_identity.0);
    signature.update(&(contact.witness.boundary_id as u64).to_le_bytes());
    signature.update(&(contact.witness.feature_id as u64).to_le_bytes());
    signature.update(&(incoming_medium as u64).to_le_bytes());
    signature.update(&(outgoing_medium as u64).to_le_bytes());
    signature.update(&[u8::from(outward), branch]);
}

#[cfg(test)]
mod tests;
