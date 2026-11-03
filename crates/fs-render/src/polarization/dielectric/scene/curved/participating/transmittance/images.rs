//! Molecular survival for each distinct discovered refractive image.
//!
//! All deterministic shooting and image deduplication finish before any source
//! tracking sample is drawn. Each retained image is replayed once on the same
//! original scene; duplicate seeds never cause another survival estimate or
//! another radiance contribution. The solved bundle remains authoritative.
//!
//! Prepared workspaces reuse native discovery and survival storage. Returned
//! views exist only after every image succeeds. Seed misses remain explicit in
//! the native catalog; a camera must not turn an empty or unresolved search
//! into certified darkness. This finite authored search is not image-complete.

use super::{
    BundleConnectionLaunch, BundleConnectionOptions, BundleConnectionTarget, CurvedRayleighScene,
    CurvedTransmittanceConnection, CurvedTransmittanceError, CurvedTransmittanceOptions,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageCatalog, BundleImageError, BundleImageOptions, BundleImageSeed, BundleImageWorkspace,
};
use fs_exec::Cx;

#[derive(Clone, Debug, PartialEq)]
pub enum CurvedTransmittanceImageError {
    Discovery(BundleImageError),
    Preparation(CurvedTransmittanceError),
    Replay {
        image: usize,
        error: CurvedTransmittanceError,
    },
    /// The complete discovery already consumed the remaining total trace cap.
    /// An image's tracked replay is charged before its native walk begins.
    TraceBudgetExceeded {
        image: usize,
    },
    AllocationFailure,
    NumericalRange,
    Cancelled,
}

impl From<BundleImageError> for CurvedTransmittanceImageError {
    fn from(error: BundleImageError) -> Self {
        Self::Discovery(error)
    }
}
impl From<fs_exec::Cancelled> for CurvedTransmittanceImageError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl core::fmt::Display for CurvedTransmittanceImageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "molecular optical images: {self:?}")
    }
}
impl core::error::Error for CurvedTransmittanceImageError {}

/// Storage prepared once per worker and reused across vertices/wavelengths.
/// Capacities are hard admission bounds; discovery/replay never grows them.
/// A failed call exposes no survival catalog; the next call clears prior work.
pub struct CurvedTransmittanceImageWorkspace {
    geometric: BundleImageWorkspace,
    paths: Vec<CurvedTransmittanceConnection>,
}

impl CurvedTransmittanceImageWorkspace {
    pub fn try_new(
        max_images: usize,
        max_seeds: usize,
    ) -> Result<Self, CurvedTransmittanceImageError> {
        let geometric = BundleImageWorkspace::try_new(max_images, max_seeds)?;
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(max_images)
            .map_err(|_| CurvedTransmittanceImageError::AllocationFailure)?;
        Ok(Self { geometric, paths })
    }

    pub fn image_capacity(&self) -> usize {
        self.geometric.image_capacity()
    }
    pub fn seed_capacity(&self) -> usize {
        self.geometric.seed_capacity()
    }
}

/// Complete borrowed survival estimates in the native catalog's image order.
/// Aggregate work includes every attempted seed, including misses/duplicates,
/// plus exactly one replay per retained image when scattering requires it.
/// Individual entries retain their original connection's own solve/replay work;
/// sum the aggregate getters, not those per-image connection totals.
#[derive(Clone, Copy, Debug)]
pub struct CurvedTransmittanceImageCatalog<'a> {
    geometric: &'a BundleImageCatalog,
    paths: &'a [CurvedTransmittanceConnection],
    trace_evaluations: u32,
    field_evaluations: u64,
    query_work: u64,
    free_flight_draws: u32,
    tracking_candidates: u32,
}

impl<'a> CurvedTransmittanceImageCatalog<'a> {
    pub const fn geometric_catalog(&self) -> &'a BundleImageCatalog {
        self.geometric
    }
    pub const fn images(&self) -> &'a [CurvedTransmittanceConnection] {
        self.paths
    }
    pub const fn trace_evaluations(&self) -> u32 {
        self.trace_evaluations
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.field_evaluations
    }
    pub const fn query_work(&self) -> u64 {
        self.query_work
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
    pub const fn tracking_candidates(&self) -> u32 {
        self.tracking_candidates
    }
}

impl CurvedRayleighScene<'_, '_, '_> {
    /// Discover/deduplicate before independently tracking every unique image.
    ///
    /// The total trace cap in image_options includes all Newton trials and
    /// source replays. Tracking draw/candidate limits apply cumulatively across
    /// all images. Zero tracking budgets admit clear or uniform source paths.
    /// The sampler receives (unique_image_ordinal, replay_flight_ordinal); image
    /// ordinals follow native first-discovery order, never the number of seeds.
    ///
    /// A miss remains a recorded native seed outcome. This method estimates
    /// survival on discovered paths only; consumers requiring every declared
    /// seed to resolve must inspect geometric_catalog().seed_results().
    #[allow(clippy::too_many_arguments)]
    pub fn discover_bundle_images_with_transmittance_samples<'w>(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        connection: BundleConnectionOptions,
        seeds: &[BundleImageSeed<'_>],
        image_options: BundleImageOptions,
        tracking: CurvedTransmittanceOptions,
        workspace: &'w mut CurvedTransmittanceImageWorkspace,
        mut samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<CurvedTransmittanceImageCatalog<'w>, CurvedTransmittanceImageError> {
        cx.checkpoint()?;
        workspace.paths.clear();
        let clear = self
            .scattering_known_clear(cx, launch.wavelength_nm)
            .map_err(CurvedTransmittanceImageError::Preparation)?;
        let geometric = self.scene.discover_bundle_images_with_workspace(
            cx,
            launch,
            target,
            connection,
            seeds,
            image_options,
            &mut workspace.geometric,
        )?;
        if geometric.images().len() > workspace.paths.capacity() {
            return Err(CurvedTransmittanceImageError::NumericalRange);
        }
        let mut trace_evaluations = geometric.trace_evaluations();
        let mut field_evaluations = geometric.field_evaluations();
        let mut query_work = geometric.query_work();
        let mut free_flight_draws = 0_u32;
        let mut tracking_candidates = 0_u32;
        for (image_index, image) in geometric.images().iter().enumerate() {
            cx.checkpoint()?;
            if !clear && trace_evaluations >= image_options.max_trace_evaluations {
                return Err(CurvedTransmittanceImageError::TraceBudgetExceeded {
                    image: image_index,
                });
            }
            let ordinal = u32::try_from(image_index)
                .map_err(|_| CurvedTransmittanceImageError::NumericalRange)?;
            let seed = seeds
                .get(image.first_seed())
                .ok_or(CurvedTransmittanceImageError::NumericalRange)?;
            let original = *image.connection();
            let remaining = CurvedTransmittanceOptions {
                max_free_flight_draws: tracking
                    .max_free_flight_draws
                    .checked_sub(free_flight_draws)
                    .ok_or(CurvedTransmittanceImageError::NumericalRange)?,
                max_tracking_candidates: tracking
                    .max_tracking_candidates
                    .checked_sub(tracking_candidates)
                    .ok_or(CurvedTransmittanceImageError::NumericalRange)?,
            };
            let tracked = self
                .replay_connected_bundle_with_transmittance(
                    cx,
                    launch,
                    target,
                    connection,
                    seed.policy,
                    remaining,
                    original,
                    clear,
                    |flight| samples(ordinal, flight),
                )
                .map_err(|error| CurvedTransmittanceImageError::Replay {
                    image: image_index,
                    error,
                })?;
            // Discovery counted every solve, including the representative's.
            // Only the replay delta belongs in the catalog aggregate.
            trace_evaluations = trace_evaluations
                .checked_add(
                    tracked
                        .trace_evaluations()
                        .checked_sub(original.trace_evaluations())
                        .ok_or(CurvedTransmittanceImageError::NumericalRange)?,
                )
                .ok_or(CurvedTransmittanceImageError::NumericalRange)?;
            field_evaluations = field_evaluations
                .checked_add(
                    tracked
                        .field_evaluations()
                        .checked_sub(original.field_evaluations())
                        .ok_or(CurvedTransmittanceImageError::NumericalRange)?,
                )
                .ok_or(CurvedTransmittanceImageError::NumericalRange)?;
            query_work = query_work
                .checked_add(
                    tracked
                        .query_work()
                        .checked_sub(original.query_work())
                        .ok_or(CurvedTransmittanceImageError::NumericalRange)?,
                )
                .ok_or(CurvedTransmittanceImageError::NumericalRange)?;
            free_flight_draws = free_flight_draws
                .checked_add(tracked.free_flight_draws())
                .ok_or(CurvedTransmittanceImageError::NumericalRange)?;
            tracking_candidates = tracking_candidates
                .checked_add(tracked.tracking_candidates())
                .ok_or(CurvedTransmittanceImageError::NumericalRange)?;
            workspace.paths.push(tracked);
        }
        cx.checkpoint()?;
        Ok(CurvedTransmittanceImageCatalog {
            geometric,
            paths: &workspace.paths,
            trace_evaluations,
            field_evaluations,
            query_work,
            free_flight_draws,
            tracking_candidates,
        })
    }
}

