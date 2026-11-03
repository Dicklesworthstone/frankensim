//! Deterministic discovery of regular images from a bounded set of launch seeds.
//!
//! Different images can traverse the same ordered native surfaces. Each seed
//! runs the existing differentiated local connector, then its converged launch
//! direction is compared only with images having the same optical history.
//! Numerical clustering has separate duplicate and separation thresholds; an
//! unresolved pair refuses the catalog instead of silently losing or doubling
//! optical power. The inverse endpoint Jacobian supplies a local angular
//! resolution estimate, not a global uniqueness or integration-error bound.
//!
//! Every supplied seed is attempted in order. A successful catalog describes
//! those finite attempts only: missed basins, caustics and unseeded histories
//! may contain other images. An empty catalog is not a zero-light certificate.
//! Each seed must initially hit the native finite receiver on a regular branch
//! for the local connector to refine it. No virtual receiver or finite-
//! difference production derivative is introduced.

use super::{
    BundleBranchPolicy, BundleConnection, BundleConnectionError, BundleConnectionLaunch,
    BundleConnectionOptions, BundleConnectionTarget, BundleSurfaceSource, CurvedPolarizedScene,
    Work, cross, norm, unit,
};
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleError, ProjectionOptions};
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSceneError;
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

/// Maximum authored launch seeds in one discovery operation.
pub const MAX_BUNDLE_IMAGE_SEEDS: usize = 256;
/// Maximum distinct regular images retained in one catalog.
pub const MAX_BUNDLE_IMAGES: usize = 32;

/// A launch guess and deterministic optical choices. Directions are normalized
/// before tracing; source position, spectrum and polarization law remain fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleImageSeed<'a> {
    pub initial_direction: Vec3,
    pub policy: BundleBranchPolicy<'a>,
}

/// Hard work/storage bounds and explicit numerical image resolution.
/// `duplicate_tolerance_rad < min_separation_rad` leaves an ambiguity band.
/// Local angular resolution is included before either classification is made.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleImageOptions {
    pub max_images: usize,
    pub min_separation_rad: f64,
    pub duplicate_tolerance_rad: f64,
    pub max_trace_evaluations: u32,
}

impl Default for BundleImageOptions {
    fn default() -> Self {
        Self {
            max_images: MAX_BUNDLE_IMAGES,
            min_separation_rad: 1e-3,
            duplicate_tolerance_rad: 1e-4,
            max_trace_evaluations: 24_576,
        }
    }
}

impl BundleImageOptions {
    pub(crate) fn validate(self) -> Result<(), BundleImageError> {
        if self.max_images == 0
            || self.max_images > MAX_BUNDLE_IMAGES
            || self.max_trace_evaluations == 0
            || self.max_trace_evaluations > (MAX_BUNDLE_IMAGE_SEEDS as u32) * 4096
            || !self.min_separation_rad.is_finite()
            || self.min_separation_rad > core::f64::consts::PI
            || !self.duplicate_tolerance_rad.is_finite()
            || self.duplicate_tolerance_rad <= 0.0
            || self.duplicate_tolerance_rad >= self.min_separation_rad
        {
            return Err(BundleImageError::InvalidOptions);
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BundleImageError {
    EmptySeeds,
    TooManySeeds,
    InvalidOptions,
    AllocationFailure,
    WorkspaceCapacity {
        required_images: usize,
        required_seeds: usize,
        image_capacity: usize,
        seed_capacity: usize,
    },
    TraceBudgetExceeded {
        seed: usize,
    },
    ImageCapacity {
        seed: usize,
        capacity: usize,
    },
    AmbiguousImage {
        seed: usize,
        image: usize,
        separation_rad: f64,
    },
    /// Invalid inputs, cancellation, and native field/interface/geometry errors
    /// abort discovery. No partially enumerated catalog is returned.
    Seed {
        seed: usize,
        error: BundleConnectionError,
    },
    NumericalRange,
}

impl core::fmt::Display for BundleImageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native optical image discovery: {self:?}")
    }
}
impl core::error::Error for BundleImageError {}

/// One distinct numerical image and the full native optical solution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleImage {
    connection: BundleConnection,
    first_seed: usize,
    angular_resolution_rad: f64,
    projection_parity: i8,
}

impl BundleImage {
    pub const fn connection(&self) -> &BundleConnection {
        &self.connection
    }
    pub const fn first_seed(&self) -> usize {
        self.first_seed
    }
    /// Local linearized angular resolution from the receiving-plane footprint
    /// and native contact collar. This is not a certified root enclosure.
    pub const fn angular_resolution_rad(&self) -> f64 {
        self.angular_resolution_rad
    }
    /// Sign of the native receiving-area differential in the launch frame.
    /// Opposite parities cannot represent one regular image, even when a fold
    /// brings their launch directions inside the duplicate tolerance.
    pub const fn projection_parity(&self) -> i8 {
        self.projection_parity
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BundleImageSeedOutcome {
    NewImage {
        image: usize,
    },
    Duplicate {
        image: usize,
    },
    /// Completed geometric misses or unresolved bounded local connections.
    /// These do not prove that the corresponding basin contains no image.
    Miss(BundleConnectionError),
}

#[derive(Clone, Debug, PartialEq)]
pub struct BundleImageSeedResult {
    outcome: BundleImageSeedOutcome,
    work: Work,
}

impl BundleImageSeedResult {
    pub const fn outcome(&self) -> &BundleImageSeedOutcome {
        &self.outcome
    }
    pub const fn trace_evaluations(&self) -> u32 {
        self.work.traces
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.work.fields
    }
    pub const fn query_work(&self) -> u64 {
        self.work.geometry
    }
}

/// Successful completion of every supplied seed, including recorded misses.
/// Storage is bounded by the public limits and allocated once before tracing.
/// Images retain first-discovery order; duplicate seeds keep the representative
/// with the smallest endpoint residual, with earlier seeds winning exact ties.
#[derive(Clone, Debug, PartialEq)]
pub struct BundleImageCatalog {
    images: Vec<BundleImage>,
    seeds: Vec<BundleImageSeedResult>,
    work: Work,
}

impl BundleImageCatalog {
    pub fn images(&self) -> &[BundleImage] {
        &self.images
    }
    pub fn seed_results(&self) -> &[BundleImageSeedResult] {
        &self.seeds
    }
    /// Includes successful, duplicate and missed local connection trials.
    pub const fn trace_evaluations(&self) -> u32 {
        self.work.traces
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.work.fields
    }
    pub const fn query_work(&self) -> u64 {
        self.work.geometry
    }

    /// Admit this finite catalog as a light-sampling proposal paired with the
    /// original continuous camera continuation. Completed geometric misses
    /// have no light-technique support; an emitter hit absent from this SAME
    /// catalog must retain its full continuation weight. An empty catalog is
    /// valid here, but is not a certificate of zero physical illumination.
    ///
    /// This accepts only native source-side exclusion, proven escape, opaque
    /// or detector obstruction, a different receiver, and the wrong receiving
    /// hemisphere. Work limits, missing branch decisions, unresolved contact,
    /// failed Newton refinement and singular/grazing images remain errors.
    /// Standalone incident-irradiance consumers have no complementary camera
    /// estimator and must keep their stricter all-seeds-resolved admission.
    pub fn validate_proposal_support(&self) -> Result<(), BundleImageError> {
        for (seed, result) in self.seeds.iter().enumerate() {
            if let BundleImageSeedOutcome::Miss(error) = &result.outcome
                && !is_geometric_support_miss(error)
            {
                return Err(BundleImageError::Seed { seed, error: error.clone() });
            }
        }
        Ok(())
    }
}

/// Reusable storage for native image discovery, prepared outside the path loop.
/// Capacities are explicit admission limits even when an allocator supplies
/// more storage. Discovery never grows these buffers. A failed operation clears
/// the incomplete result while retaining both allocations for the next call.
///
/// Completed catalogs are borrowed only from a successful discovery call;
/// mutably reusing this workspace cannot coexist with a prior catalog borrow.
#[derive(Debug)]
pub struct BundleImageWorkspace {
    catalog: BundleImageCatalog,
    image_capacity: usize,
    seed_capacity: usize,
}

impl BundleImageWorkspace {
    pub fn try_new(max_images: usize, max_seeds: usize) -> Result<Self, BundleImageError> {
        if max_images == 0
            || max_images > MAX_BUNDLE_IMAGES
            || max_seeds == 0
            || max_seeds > MAX_BUNDLE_IMAGE_SEEDS
        {
            return Err(BundleImageError::InvalidOptions);
        }
        let mut catalog = BundleImageCatalog {
            images: Vec::new(),
            seeds: Vec::new(),
            work: Work::default(),
        };
        catalog
            .images
            .try_reserve_exact(max_images)
            .map_err(|_| BundleImageError::AllocationFailure)?;
        catalog
            .seeds
            .try_reserve_exact(max_seeds)
            .map_err(|_| BundleImageError::AllocationFailure)?;
        Ok(Self {
            catalog,
            image_capacity: max_images,
            seed_capacity: max_seeds,
        })
    }

    pub const fn image_capacity(&self) -> usize {
        self.image_capacity
    }

    pub const fn seed_capacity(&self) -> usize {
        self.seed_capacity
    }

    fn clear(&mut self) {
        self.catalog.images.clear();
        self.catalog.seeds.clear();
        self.catalog.work = Work::default();
    }
}

impl CurvedPolarizedScene<'_> {
    /// Discover distinct regular images from an authored finite seed list.
    /// A per-seed local trace limit is a recorded miss; exhausting the total
    /// discovery limit aborts. Native failures/cancellation and ambiguous image
    /// separation abort, so a partial search cannot look like completed work.
    pub fn discover_bundle_images(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>],
        options: BundleImageOptions,
    ) -> Result<BundleImageCatalog, BundleImageError> {
        validate_discovery_request(connection, seeds, options)?;
        let mut workspace = BundleImageWorkspace::try_new(options.max_images, seeds.len())?;
        self.discover_bundle_images_with_workspace(
            cx,
            launch,
            target,
            connection,
            seeds,
            options,
            &mut workspace,
        )?;
        Ok(workspace.catalog)
    }

    /// Discover through the same native shooting and image-classification
    /// owners using caller-prepared storage. Every authored seed is attempted
    /// in order, including duplicate seeds and recorded local misses. Their
    /// work remains part of the complete catalog's totals.
    ///
    /// Insufficient prepared capacity refuses before tracing; it never triggers
    /// a path-time allocation. Any failure invalidates the incomplete catalog
    /// while preserving storage for subsequent pixels or wavelength lanes.
    #[allow(clippy::too_many_arguments)]
    pub fn discover_bundle_images_with_workspace<'w>(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>],
        options: BundleImageOptions,
        workspace: &'w mut BundleImageWorkspace,
    ) -> Result<&'w BundleImageCatalog, BundleImageError> {
        self.discover_bundle_images_from_source_with_workspace(
            cx, launch, None, target, connection, seeds, options, workspace,
        )
    }

    /// Discover finite-seed images emitted from one original native surface.
    /// Each trial re-admits that surface through its actual departure oracle;
    /// the source's original geometry remains active for all later encounters.
    #[allow(clippy::too_many_arguments)]
    pub fn discover_bundle_images_from_surface_with_workspace<'w>(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        source: BundleSurfaceSource,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>],
        options: BundleImageOptions,
        workspace: &'w mut BundleImageWorkspace,
    ) -> Result<&'w BundleImageCatalog, BundleImageError> {
        self.discover_bundle_images_from_source_with_workspace(
            cx,
            launch,
            Some(source),
            target,
            connection,
            seeds,
            options,
            workspace,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn discover_bundle_images_from_source_with_workspace<'w>(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        source: Option<BundleSurfaceSource>,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>],
        options: BundleImageOptions,
        workspace: &'w mut BundleImageWorkspace,
    ) -> Result<&'w BundleImageCatalog, BundleImageError> {
        workspace.clear();
        validate_discovery_request(connection, seeds, options)?;
        if options.max_images > workspace.image_capacity || seeds.len() > workspace.seed_capacity {
            return Err(BundleImageError::WorkspaceCapacity {
                required_images: options.max_images,
                required_seeds: seeds.len(),
                image_capacity: workspace.image_capacity,
                seed_capacity: workspace.seed_capacity,
            });
        }
        if let Err(error) = self.discover_bundle_images_prepared(
            cx,
            launch,
            source,
            target,
            connection,
            seeds,
            options,
            &mut workspace.catalog,
        ) {
            workspace.clear();
            return Err(error);
        }
        Ok(&workspace.catalog)
    }

    #[allow(clippy::too_many_arguments)]
    fn discover_bundle_images_prepared(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        source: Option<BundleSurfaceSource>,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>],
        options: BundleImageOptions,
        catalog: &mut BundleImageCatalog,
    ) -> Result<(), BundleImageError> {
        for (seed_index, seed) in seeds.iter().enumerate() {
            cx.checkpoint().map_err(|_| BundleImageError::Seed {
                seed: seed_index,
                error: BundleConnectionError::Cancelled,
            })?;
            let remaining = options.max_trace_evaluations - catalog.work.traces;
            if remaining == 0 {
                return Err(BundleImageError::TraceBudgetExceeded { seed: seed_index });
            }
            let direction = unit(seed.initial_direction).ok_or(BundleImageError::Seed {
                seed: seed_index,
                error: BundleConnectionError::NumericalRange,
            })?;
            let local_options = BundleConnectionOptions {
                max_trace_evaluations: connection.max_trace_evaluations.min(remaining),
                ..connection
            };
            let mut work = Work::default();
            let result = self.connect_bundle_from_source_with_policy_work(
                cx,
                BundleConnectionLaunch {
                    ray: Ray {
                        origin: launch.ray.origin,
                        dir: direction,
                    },
                    ..launch
                },
                source,
                target,
                local_options,
                seed.policy,
                &mut work,
            );
            if let Err(BundleConnectionError::Scene(CurvedBundleSceneError::SourceSideExcluded {
                field_evaluations, query_work, ..
            })) = &result {
                // A supported exclusion happens before a native path result
                // exists. Charge its actual launch and geometry query once;
                // all preceding completed Newton trials are already in work.
                work.fields = work.fields.checked_add(*field_evaluations)
                    .ok_or(BundleImageError::NumericalRange)?;
                work.geometry = work.geometry.checked_add(*query_work)
                    .ok_or(BundleImageError::NumericalRange)?;
            }
            add_work(&mut catalog.work, work)?;
            let outcome = match result {
                Ok(solution) => {
                    let resolution =
                        connection_angular_resolution(&solution, target, connection.projection)
                            .map_err(|error| BundleImageError::Seed {
                                seed: seed_index,
                                error,
                            })?;
                    let projection = solution
                        .receiver_projection(target, connection.projection)
                        .map_err(|error| BundleImageError::Seed {
                            seed: seed_index,
                            error,
                        })?;
                    let candidate = BundleImage {
                        connection: solution,
                        first_seed: seed_index,
                        angular_resolution_rad: resolution,
                        projection_parity: if projection.signed_area_per_solid_angle_m2() < 0.0 {
                            -1
                        } else {
                            1
                        },
                    };
                    classify_candidate(&mut catalog.images, candidate, options)?
                }
                Err(BundleConnectionError::TraceBudgetExceeded)
                    if catalog.work.traces == options.max_trace_evaluations =>
                {
                    return Err(BundleImageError::TraceBudgetExceeded { seed: seed_index });
                }
                Err(error) if is_local_miss(&error) => BundleImageSeedOutcome::Miss(error),
                Err(error) => {
                    return Err(BundleImageError::Seed {
                        seed: seed_index,
                        error,
                    });
                }
            };
            catalog.seeds.push(BundleImageSeedResult { outcome, work });
        }
        Ok(())
    }
}

fn validate_discovery_request(
    connection: BundleConnectionOptions,
    seeds: &[BundleImageSeed<'_>],
    options: BundleImageOptions,
) -> Result<(), BundleImageError> {
    options.validate()?;
    connection
        .validate()
        .map_err(|error| BundleImageError::Seed { seed: 0, error })?;
    if seeds.is_empty() {
        return Err(BundleImageError::EmptySeeds);
    }
    if seeds.len() > MAX_BUNDLE_IMAGE_SEEDS {
        return Err(BundleImageError::TooManySeeds);
    }
    Ok(())
}

fn add_work(total: &mut Work, next: Work) -> Result<(), BundleImageError> {
    total.traces = total
        .traces
        .checked_add(next.traces)
        .ok_or(BundleImageError::NumericalRange)?;
    total.fields = total
        .fields
        .checked_add(next.fields)
        .ok_or(BundleImageError::NumericalRange)?;
    total.geometry = total
        .geometry
        .checked_add(next.geometry)
        .ok_or(BundleImageError::NumericalRange)?;
    Ok(())
}

fn is_local_miss(error: &BundleConnectionError) -> bool {
    matches!(
        error,
        BundleConnectionError::ReceiverNotReached(_)
            | BundleConnectionError::Scene(CurvedBundleSceneError::SourceSideExcluded { .. })
            | BundleConnectionError::WrongDetector
            | BundleConnectionError::BackFacingReceiver
            | BundleConnectionError::BranchChanged
            | BundleConnectionError::BranchPlanNotConsumed { .. }
            | BundleConnectionError::ContactResolution { .. }
            | BundleConnectionError::IterationLimit { .. }
            | BundleConnectionError::TraceBudgetExceeded
            | BundleConnectionError::NoDescent { .. }
            | BundleConnectionError::Projection(
                BundleError::SingularProjection | BundleError::GrazingProjection
            )
    )
}

fn is_geometric_support_miss(error: &BundleConnectionError) -> bool {
    matches!(error,
        BundleConnectionError::Scene(CurvedBundleSceneError::SourceSideExcluded { .. })
            | BundleConnectionError::ReceiverNotReached(
                PathTermination::Opaque | PathTermination::Detected | PathTermination::Escaped
            )
            | BundleConnectionError::WrongDetector
            | BundleConnectionError::BackFacingReceiver
    )
}

fn classify_candidate(
    images: &mut Vec<BundleImage>,
    candidate: BundleImage,
    options: BundleImageOptions,
) -> Result<BundleImageSeedOutcome, BundleImageError> {
    let mut duplicate = None;
    let direction = candidate
        .connection
        .path()
        .bundle()
        .propagation()
        .initial_ray()
        .dir;
    for (index, image) in images.iter().enumerate() {
        if image.connection.path().branch_signature()
            != candidate.connection.path().branch_signature()
        {
            continue;
        }
        let separation = angular_distance(
            direction,
            image
                .connection
                .path()
                .bundle()
                .propagation()
                .initial_ray()
                .dir,
        );
        match classify_image_pair(
            options,
            separation,
            candidate.angular_resolution_rad,
            image.angular_resolution_rad,
        ) {
            ImagePairClassification::Distinct => {}
            ImagePairClassification::Duplicate
                if duplicate.is_none()
                    && candidate.projection_parity == image.projection_parity =>
            {
                duplicate = Some(index);
            }
            ImagePairClassification::Duplicate | ImagePairClassification::Ambiguous => {
                return Err(BundleImageError::AmbiguousImage {
                    seed: candidate.first_seed,
                    image: index,
                    separation_rad: separation,
                });
            }
        }
    }
    if let Some(index) = duplicate {
        if candidate.connection.residual_m() < images[index].connection.residual_m() {
            let first_seed = images[index].first_seed;
            images[index] = BundleImage {
                first_seed,
                ..candidate
            };
        }
        return Ok(BundleImageSeedOutcome::Duplicate { image: index });
    }
    if images.len() == options.max_images {
        return Err(BundleImageError::ImageCapacity {
            seed: candidate.first_seed,
            capacity: options.max_images,
        });
    }
    let image = images.len();
    images.push(candidate);
    Ok(BundleImageSeedOutcome::NewImage { image })
}

/// Shared with spectral continuation so image-merger checks use the same
/// actual differentiated footprint as the discovery operation.
pub(crate) fn connection_angular_resolution(
    connection: &BundleConnection,
    target: BundleConnectionTarget,
    options: ProjectionOptions,
) -> Result<f64, BundleConnectionError> {
    let projection = connection.receiver_projection(target, options)?;
    // area = sigma_max*sigma_min and condition = sigma_max/sigma_min.
    // Separate square roots avoid underflow in area/condition.
    let sigma_min =
        det::sqrt(projection.area_per_solid_angle_m2()) / det::sqrt(projection.condition_number());
    let spatial_resolution =
        target.position_tolerance_m() + connection.path().max_contact_resolution_m();
    let resolution = spatial_resolution / sigma_min;
    if !resolution.is_finite() || resolution <= 0.0 {
        return Err(BundleConnectionError::NumericalRange);
    }
    Ok(resolution.max(64.0 * f64::EPSILON))
}

pub(crate) fn angular_distance(first: Vec3, second: Vec3) -> f64 {
    let (Some(first), Some(second)) = (unit(first), unit(second)) else {
        return f64::NAN;
    };
    det::atan2(
        norm(cross(first, second)),
        first.dot(second).clamp(-1.0, 1.0),
    )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ImagePairClassification {
    Distinct,
    Duplicate,
    Ambiguous,
}

pub(crate) fn classify_image_pair(
    options: BundleImageOptions,
    separation_rad: f64,
    first_resolution_rad: f64,
    second_resolution_rad: f64,
) -> ImagePairClassification {
    if !separation_rad.is_finite()
        || separation_rad < 0.0
        || !first_resolution_rad.is_finite()
        || first_resolution_rad < 0.0
        || !second_resolution_rad.is_finite()
        || second_resolution_rad < 0.0
    {
        return ImagePairClassification::Ambiguous;
    }
    let resolution = first_resolution_rad + second_resolution_rad;
    if separation_rad + resolution <= options.duplicate_tolerance_rad {
        ImagePairClassification::Duplicate
    } else if separation_rad - resolution >= options.min_separation_rad {
        ImagePairClassification::Distinct
    } else {
        ImagePairClassification::Ambiguous
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod support_tests;
