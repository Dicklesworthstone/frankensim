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
use crate::gradient_index::event::{
    BoundaryQuery, BoundaryWitness, DepartureNeighborhood, DepartureRequest, SignedBoundaryDistance,
};
use crate::polarization::dielectric::transport::bundle::{
    PolarizedBundleError, PolarizedRayBundle,
};

/// Bounded local source-to-receiver connections through these native paths.
pub mod connection;

/// An explicit emission/departure intent at one original native surface.
/// This descriptor contains no caller-supplied geometric witness. The native
/// walk queries the unchanged scene at the actual launch point for every
/// trial and replay before admitting its selected geometric side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BundleSurfaceSource {
    object_id: u64,
    side: BoundarySide,
}

impl BundleSurfaceSource {
    pub fn try_new(
        object_id: u64,
        side: BoundarySide,
    ) -> Result<Self, CurvedBundleSceneError> {
        if side == BoundarySide::OnSurface {
            return Err(CurvedBundleSceneError::InvalidSurfaceSource { object_id });
        }
        Ok(Self { object_id, side })
    }

    pub const fn object_id(self) -> u64 {
        self.object_id
    }

    pub const fn side(self) -> BoundarySide {
        self.side
    }
}

/// A numerical stopping section is not a scene boundary or optical interaction.
/// Its medium is always an index into this exact immutable native scene.
#[derive(Clone, Copy, Debug)]
pub(crate) enum BundleReceiver {
    Diffuse(u64),
    Volume {
        medium: usize,
        position_m: Point3,
        normal: Vec3,
    },
}

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
    volume_receiver: bool,
    surface_events: u32,
    last_contact: Option<CurvedContact>,
    source_contact: Option<CurvedContact>,
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
    /// The numerical section for a point in the selected original medium was
    /// reached. No surface interaction or event was added at that endpoint.
    pub const fn is_volume_receiver(&self) -> bool {
        self.volume_receiver
    }
    pub const fn surface_events(&self) -> u32 {
        self.surface_events
    }
    pub const fn last_contact(&self) -> Option<&CurvedContact> {
        self.last_contact.as_ref()
    }
    /// Original native source admission, when emission begins on a surface.
    /// This records geometry/side resolution only; it is not a resume token.
    pub const fn source_contact(&self) -> Option<&CurvedContact> {
        self.source_contact.as_ref()
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
    /// The selected original surface, outgoing side or emission role was not
    /// admitted at the unsnapped launch point.
    InvalidSurfaceSource {
        object_id: u64,
    },
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
        receiver: Option<BundleReceiver>,
        branches: impl FnMut(
            &CurvedContact,
            &CurvedBundleBranchContext,
        ) -> Result<DielectricEvent, CurvedBundleSceneError>,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        self.trace_bundle_with_receiver_flights(
            cx,
            bundle,
            initial_medium,
            settings,
            continue_extinguished,
            receiver,
            branches,
            |_, _| Ok(None),
            |_, _, _, _, _| Ok(()),
        )
    }

    // Optional source-survival candidates share this exact augmented boundary
    // walk. The callbacks cannot replace geometry, interfaces, or derivatives.
    // A segment observation receives incident medium, wavelength, actual arc
    // length, native endpoint, and whether the endpoint is a tracking candidate.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_bundle_with_receiver_flights(
        &self,
        cx: &Cx<'_>,
        bundle: PolarizedRayBundle,
        initial_medium: usize,
        settings: CurvedBundleSettings,
        continue_extinguished: bool,
        receiver: Option<BundleReceiver>,
        branches: impl FnMut(
            &CurvedContact,
            &CurvedBundleBranchContext,
        ) -> Result<DielectricEvent, CurvedBundleSceneError>,
        flight: impl FnMut(usize, f64) -> Result<Option<f64>, CurvedBundleSceneError>,
        observe_segment: impl FnMut(
            usize,
            f64,
            f64,
            Point3,
            bool,
        ) -> Result<(), CurvedBundleSceneError>,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        self.trace_bundle_with_source_receiver_flights(
            cx,
            bundle,
            initial_medium,
            settings,
            continue_extinguished,
            None,
            receiver,
            branches,
            flight,
            observe_segment,
        )
    }

    // Surface emission and ordinary point launches share the entire native
    // event owner. Source admission seeds its original departure neighborhood;
    // it applies no source-surface optical interaction and invents no offset.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn trace_bundle_with_source_receiver_flights(
        &self,
        cx: &Cx<'_>,
        bundle: PolarizedRayBundle,
        initial_medium: usize,
        settings: CurvedBundleSettings,
        continue_extinguished: bool,
        source: Option<BundleSurfaceSource>,
        receiver: Option<BundleReceiver>,
        mut branches: impl FnMut(
            &CurvedContact,
            &CurvedBundleBranchContext,
        ) -> Result<DielectricEvent, CurvedBundleSceneError>,
        mut flight: impl FnMut(usize, f64) -> Result<Option<f64>, CurvedBundleSceneError>,
        mut observe_segment: impl FnMut(
            usize,
            f64,
            f64,
            Point3,
            bool,
        ) -> Result<(), CurvedBundleSceneError>,
    ) -> Result<CurvedBundlePathResult, CurvedBundleSceneError> {
        cx.checkpoint()?;
        settings.validate()?;
        self.media
            .get(initial_medium)
            .ok_or(SceneError::InvalidInitialMedium)?;
        if let Some(BundleReceiver::Volume { medium, .. }) = receiver {
            self.media
                .get(medium)
                .ok_or(SceneError::InvalidInitialMedium)?;
        }
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
            volume_receiver: false,
            surface_events: 0,
            last_contact: None,
            source_contact: None,
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
        let mut departure = if let Some(source) = source {
            // Query with THIS oracle: smooth departure admission depends on
            // its original last-contact identity and competing-surface bounds.
            let ray = result.bundle.propagation().initial_ray();
            let witness = oracle
                .query(ray.origin)
                .map_err(|error| {
                    oracle.last_error().map_or_else(
                        || CurvedBundleSceneError::from(error),
                        CurvedBundleSceneError::from,
                    )
                })?
                .witness
                .ok_or(CurvedBundleSceneError::InvalidSurfaceSource {
                    object_id: source.object_id,
                })?;
            result.boundary_evaluations = 1;
            let instance = self
                .geometry
                .instances()
                .get(witness.boundary_id)
                .ok_or(CurvedSceneError::InvalidContact)?;
            let binding = self
                .bindings
                .get(witness.boundary_id)
                .ok_or(CurvedSceneError::InvalidContact)?;
            let invalid = || CurvedBundleSceneError::InvalidSurfaceSource {
                object_id: source.object_id,
            };
            if instance.object_id() != source.object_id
                || binding.object_id != source.object_id
                || !matches!(source.side, BoundarySide::Positive | BoundarySide::Negative)
            {
                return Err(invalid());
            }
            let cosine = witness.normal.dot(ray.dir);
            let positive = source.side == BoundarySide::Positive;
            if !cosine.is_finite()
                || (positive && cosine <= 0.0)
                || (!positive && cosine >= 0.0)
            {
                return Err(invalid());
            }
            match binding.interaction {
                Interaction::Detector { .. } => return Err(invalid()),
                Interaction::Opaque => {}
                Interaction::DepolarizingDiffuse { medium, .. }
                | Interaction::SmoothConductor { medium, .. }
                | Interaction::RoughConductor { medium, .. }
                | Interaction::BrushedConductor { medium, .. } => {
                    require_medium(source.object_id, initial_medium, medium)?;
                }
                Interaction::Dielectric { outside, inside, .. }
                | Interaction::CoatedDielectric { outside, inside, .. } => {
                    require_medium(
                        source.object_id,
                        initial_medium,
                        if positive { outside } else { inside },
                    )?;
                }
            }
            let contact = CurvedContact {
                object_id: instance.object_id(),
                geometry_identity: instance.geometry_identity(),
                frame_identity: instance.frame_identity(),
                witness,
                resolution_m: witness.distance_upper_m / cosine.abs(),
            };
            result.max_contact_resolution_m = contact.resolution_m;
            result.source_contact = Some(contact);
            signature.update(b"surface-source");
            append_contact(
                &mut signature,
                contact,
                initial_medium,
                initial_medium,
                positive,
                6,
            );
            Some(GradientDeparture {
                witness,
                side: source.side,
            })
        } else {
            None
        };
        let mut remaining_steps = settings.trace.max_total_steps;
        let mut first_segment = true;
        while result.surface_events < settings.trace.max_surface_events {
            cx.checkpoint()?;
            let event_index = result.surface_events;
            let remaining_distance = settings.trace.max_distance_m
                - result.bundle.propagation().propagation().distance_m;
            if remaining_distance <= 0.0 {
                result.termination = PathTermination::DistanceLimit;
                break;
            }
            if remaining_steps == 0 {
                return Err(GradientIndexError::StepBudgetExceeded.into());
            }
            let wavelength = result.bundle.propagation().wavelength_nm();
            let tracking_flight = flight(result.medium, wavelength)?;
            if tracking_flight.is_some_and(|distance| distance.is_nan() || distance < 0.0) {
                return Err(CurvedBundleSceneError::NumericalRange);
            }
            let segment_distance = tracking_flight.map_or(remaining_distance, |distance| {
                distance.min(remaining_distance)
            });
            let candidate_due =
                tracking_flight.is_some_and(|distance| distance < remaining_distance);
            let mut boundary = settings.trace.event;
            if first_segment && source.is_some() {
                // The source query is part of this first segment's oracle
                // allowance, as well as the scene's cumulative geometry work.
                boundary.max_boundary_queries = boundary
                    .max_boundary_queries
                    .checked_sub(1)
                    .filter(|remaining| *remaining > 0)
                    .ok_or(GradientIndexError::BoundaryQueryBudgetExceeded)?;
            }
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
            let section = match receiver {
                Some(BundleReceiver::Volume {
                    medium,
                    position_m,
                    normal,
                }) if medium == result.medium => Some((position_m, normal)),
                _ => None,
            };
            let propagation = {
                let mut endpoint_oracle = TargetSectionOracle {
                    native: &mut oracle,
                    section,
                    surface_tolerance_m: settings.trace.event.surface_tolerance_m,
                };
                result.bundle.continue_to_boundary_with_control(
                    field,
                    &mut endpoint_oracle,
                    segment_distance,
                    departure.take(),
                    options,
                    || cx.checkpoint().is_ok(),
                )
            };
            let (next, segment) = propagation.map_err(|error| {
                if let Some(geometry) = oracle.last_error() {
                    CurvedBundleSceneError::from(geometry)
                } else if first_segment
                    && error == PolarizedBundleError::Bundle(BundleError::InitialSampleMismatch)
                {
                    SceneError::InvalidInitialMedium.into()
                } else {
                    error.into()
                }
            })?;
            first_segment = false;
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
            let endpoint = result.bundle.propagation().propagation().ray.origin;
            if candidate_due {
                match segment.termination {
                    BoundaryTermination::DistanceLimit => {
                        // A candidate inside an unresolved arrival/departure
                        // collar cannot be ordered against the real interface.
                        let clearance = oracle
                            .query(endpoint)
                            .map_err(|error| {
                                oracle.last_error().map_or_else(
                                    || CurvedBundleSceneError::from(error),
                                    CurvedBundleSceneError::from,
                                )
                            })?
                            .clearance_m;
                        if clearance <= 2.0 * settings.trace.event.surface_tolerance_m {
                            return Err(GradientIndexError::BoundaryResolutionLimit.into());
                        }
                        observe_segment(
                            result.medium,
                            wavelength,
                            segment.segment.distance_m,
                            endpoint,
                            true,
                        )?;
                        // Continue the original launch-angle columns from the
                        // exact endpoint; a null candidate is not a surface.
                        continue;
                    }
                    BoundaryTermination::Surface(_) => {
                        if tracking_flight.ok_or(CurvedBundleSceneError::NumericalRange)?
                            - segment.segment.distance_m
                            <= segment.contact_resolution_m
                        {
                            return Err(GradientIndexError::BoundaryResolutionLimit.into());
                        }
                    }
                }
            }
            observe_segment(
                result.medium,
                wavelength,
                segment.segment.distance_m,
                endpoint,
                false,
            )?;
            let BoundaryTermination::Surface(witness) = segment.termination else {
                result.termination = if tracking_flight.is_none()
                    && self.straight_ray_escapes(
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
            if witness.boundary_id == TARGET_SECTION_ID {
                // This is an explicit numerical endpoint, not a physical
                // surface. In particular it has no BSDF, detector, object ID,
                // interface derivative, or charge against surface events.
                result.volume_receiver = true;
                result.termination = PathTermination::DistanceLimit;
                signature.update(b"volume-receiver");
                signature.update(&(result.medium as u64).to_le_bytes());
                append_medium(&mut signature, result.medium_identity);
                break;
            }
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
                    if matches!(receiver, Some(BundleReceiver::Diffuse(id))
                        if id == instance.object_id()) =>
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

// A reserved event identity is consumed here before native instance lookup.
// It never enters the scene's object/geometry identity namespace.
const TARGET_SECTION_ID: usize = usize::MAX;

/// Compose a numerical endpoint section with the original geometry oracle.
/// Every query still visits that original oracle and its original work budget.
/// The section cannot mask a nearer object or lie inside an unexamined native
/// departure neighborhood. It supplies no interface or material conversion.
struct TargetSectionOracle<'a, O: BoundaryOracle + ?Sized> {
    native: &'a mut O,
    section: Option<(Point3, Vec3)>,
    surface_tolerance_m: f64,
}

impl<O: BoundaryOracle + ?Sized> TargetSectionOracle<'_, O> {
    fn distance(&self, position: Point3) -> Result<Option<(f64, f64, f64)>, GradientIndexError> {
        let Some((origin, normal)) = self.section else {
            return Ok(None);
        };
        let signed = position.delta_from(origin).dot(normal);
        // Include cancellation in world-coordinate subtraction as well as the
        // dot product. A section below this resolution is refused, not snapped.
        let roundoff = 16.0
            * f64::EPSILON
            * (normal.x.abs() * (position.x.abs() + origin.x.abs())
                + normal.y.abs() * (position.y.abs() + origin.y.abs())
                + normal.z.abs() * (position.z.abs() + origin.z.abs()));
        let upper = (signed.abs() + roundoff).next_up();
        if !signed.is_finite() || !upper.is_finite() {
            return Err(GradientIndexError::NumericalRange);
        }
        let lower = (signed.abs() - roundoff).next_down().max(0.0);
        Ok(Some((signed, lower, upper)))
    }

    fn clip_neighborhood(
        &self,
        center: Point3,
        neighborhood: DepartureNeighborhood,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        let Some((_, lower, _)) = self.distance(center)? else {
            return Ok(neighborhood);
        };
        let radius = 0.5 * lower;
        if radius <= 0.0 {
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        Ok(match neighborhood {
            DepartureNeighborhood::Planar { radius_m } => DepartureNeighborhood::Planar {
                radius_m: radius_m.min(radius),
            },
            DepartureNeighborhood::Smooth {
                radius_m,
                normal_axis,
                max_normal_deviation,
            } => DepartureNeighborhood::Smooth {
                radius_m: radius_m.min(radius),
                normal_axis,
                max_normal_deviation,
            },
        })
    }
}

impl<O: BoundaryOracle + ?Sized> BoundaryOracle for TargetSectionOracle<'_, O> {
    fn query(&mut self, position: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        let mut query = self.native.query(position)?;
        let Some((signed, lower, upper)) = self.distance(position)? else {
            return Ok(query);
        };
        if upper <= self.surface_tolerance_m && query.clearance_m <= 2.0 * self.surface_tolerance_m
        {
            // A volume endpoint cannot claim a side of a physical interface
            // inside its unresolved numerical contact collar.
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        if upper < query.clearance_m {
            let (_, normal) = self.section.ok_or(GradientIndexError::NumericalRange)?;
            let radius = if query.clearance_m.is_finite() {
                0.25 * (query.clearance_m - upper)
            } else {
                f64::MAX
            };
            query.witness = Some(BoundaryWitness {
                boundary_id: TARGET_SECTION_ID,
                feature_id: 0,
                surface_point: position.offset(normal.scale(-signed)),
                normal,
                distance_upper_m: upper,
                side: if lower == 0.0 {
                    BoundarySide::OnSurface
                } else if signed > 0.0 {
                    BoundarySide::Positive
                } else {
                    BoundarySide::Negative
                },
                departure_radius_m: (radius > 0.0).then_some(radius),
            });
        } else if let Some(witness) = query.witness.as_mut()
            && let Some(radius) = witness.departure_radius_m
        {
            let (_, separation, _) = self
                .distance(witness.surface_point)?
                .ok_or(GradientIndexError::NumericalRange)?;
            witness.departure_radius_m = (separation > 0.0).then_some(radius.min(0.5 * separation));
        }
        query.clearance_m = query.clearance_m.min(lower);
        Ok(query)
    }

    fn prepare_departure(
        &mut self,
        request: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        let neighborhood = self.native.prepare_departure(request)?;
        self.clip_neighborhood(request.contact.witness.surface_point, neighborhood)
    }

    fn admit_arrival(
        &mut self,
        request: DepartureRequest,
    ) -> Result<Option<DepartureNeighborhood>, GradientIndexError> {
        if request.contact.witness.boundary_id == TARGET_SECTION_ID {
            return Ok(Some(DepartureNeighborhood::planar(
                request.contact.witness,
            )?));
        }
        self.native
            .admit_arrival(request)?
            .map(|neighborhood| {
                self.clip_neighborhood(request.contact.witness.surface_point, neighborhood)
            })
            .transpose()
    }

    fn signed_departure(
        &mut self,
        contact: BoundaryWitness,
        position: Point3,
    ) -> Result<SignedBoundaryDistance, GradientIndexError> {
        self.native.signed_departure(contact, position)
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
