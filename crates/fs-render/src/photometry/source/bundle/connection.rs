//! Calibrated point-source measurements through automatically connected paths.
//!
//! Every native spectral bin solves its own dispersive connecting ray and uses
//! that ray's actual source intensity, polarized interface power, integrated
//! absorption and receiving-area Jacobian. Success publishes all 80 bins and
//! native CIE quadrature together; a blocked or unresolved bin publishes none.
//! The preceding wavelength supplies the next launch guess. Ordered optical
//! branch identity must remain unchanged across the spectrum.
//!
//! This is a local, single-image measurement; [`multipath`] sums a bounded set
//! of separately connected optical histories. [`images`] discovers distinct
//! regular images within those histories and tracks their spectral power.
//! The supplied initial direction
//! must reach the finite detector in the first bin's convergence basin. Other
//! images, diffuse paths, finite-aperture caustics and wavelengths outside the
//! native quadrature require their own transport/integration owners.

use super::{
    BundleIrradianceMeter, BundleIrradianceReading, BundlePhotometryError, PhotometricSource,
    PlanarIrradianceReceiver, SPECTRAL_BINS, wavelength_nm,
};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionError, BundleConnectionLaunch, BundleConnectionOptions,
    BundleConnectionTarget,
};
use crate::polarization::{Stokes, photometry::PolarizedSource};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Vec3;

/// Incoherent spectral power from distinct reflected/transmitted branches.
pub mod multipath;

/// Spectral continuation and incoherent power from discovered regular images.
pub mod images;

/// One native planar receiver and a local initial guess for the first bin.
/// Source polarization is defined relative to the supplied world axis projected
/// into each launch direction. Detector orientation remains a separate input.
#[derive(Clone, Copy, Debug)]
pub struct ConnectedIrradianceRequest {
    pub detector_object_id: u64,
    pub receiver: PlanarIrradianceReceiver,
    pub initial_direction: Vec3,
    pub source_transverse_axis: Vec3,
    pub initial_medium: usize,
    pub connection: BundleConnectionOptions,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ConnectedIrradianceError {
    InvalidTarget(BundleConnectionError),
    Connection {
        bin: usize,
        error: BundleConnectionError,
    },
    Photometry {
        bin: usize,
        error: BundlePhotometryError,
    },
    /// Dispersion left the selected smooth optical branch. A partial spectrum
    /// must not masquerade as a complete measurement of that image.
    SpectralBranchChanged {
        bin: usize,
    },
    NumericalRange,
}
impl core::fmt::Display for ConnectedIrradianceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "connected spectral irradiance: {self:?}")
    }
}
impl core::error::Error for ConnectedIrradianceError {}

/// Complete visible-band reading, actual dispersive launch directions and work.
/// Fixed-size storage; no path list or per-trial allocation is retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ConnectedIrradianceSpectrum {
    bins_w_m2_nm: [f64; SPECTRAL_BINS],
    launch_directions: [Vec3; SPECTRAL_BINS],
    reading: BundleIrradianceReading,
    branch_signature: ContentHash,
    max_endpoint_error_m: f64,
    max_contact_resolution_m: f64,
    trace_evaluations: u64,
    field_evaluations: u64,
    query_work: u64,
}

impl ConnectedIrradianceSpectrum {
    pub const fn spectral_irradiance_w_m2_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.bins_w_m2_nm
    }
    pub const fn launch_directions(&self) -> &[Vec3; SPECTRAL_BINS] {
        &self.launch_directions
    }
    pub const fn reading(&self) -> BundleIrradianceReading {
        self.reading
    }
    pub const fn branch_signature(&self) -> ContentHash {
        self.branch_signature
    }
    pub const fn max_endpoint_error_m(&self) -> f64 {
        self.max_endpoint_error_m
    }
    /// Native numerical collar, not a bound on the reading's physical error.
    pub const fn max_contact_resolution_m(&self) -> f64 {
        self.max_contact_resolution_m
    }
    pub const fn trace_evaluations(&self) -> u64 {
        self.trace_evaluations
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.field_evaluations
    }
    pub const fn query_work(&self) -> u64 {
        self.query_work
    }
}

impl PhotometricSource {
    /// Connect and measure an initially unpolarized fixture. Every encountered
    /// smooth dielectric/coating/conductor retains its induced polarization.
    pub fn connect_irradiance(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ConnectedIrradianceRequest,
    ) -> Result<ConnectedIrradianceSpectrum, ConnectedIrradianceError> {
        self.connect_irradiance_with_policy(
            cx,
            scene,
            request,
            BundleBranchPolicy::TransmitOrReflectTir,
        )
    }

    /// Connect a selected reflection/transmission sequence at every wavelength.
    /// Every shooting trial and spectral bin starts the policy at its beginning;
    /// each native interface still determines its own position and local optics.
    pub fn connect_irradiance_with_policy(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ConnectedIrradianceRequest,
        policy: BundleBranchPolicy<'_>,
    ) -> Result<ConnectedIrradianceSpectrum, ConnectedIrradianceError> {
        connect_spectrum(self, cx, scene, request, policy, |_| {
            Stokes::unpolarized(1.0)
        })
    }
}

impl PolarizedSource {
    /// Preserve the fixture's wavelength-dependent source polarization and all
    /// upstream optical losses. The native calibrated intensity view already
    /// contains those losses; only its normalized surviving state is launched.
    pub fn connect_irradiance(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ConnectedIrradianceRequest,
    ) -> Result<ConnectedIrradianceSpectrum, ConnectedIrradianceError> {
        self.connect_irradiance_with_policy(
            cx,
            scene,
            request,
            BundleBranchPolicy::TransmitOrReflectTir,
        )
    }

    /// Connect a selected branch with the source's surviving spectral Stokes
    /// states. Upstream filtering and every selected interface loss apply once.
    pub fn connect_irradiance_with_policy(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ConnectedIrradianceRequest,
        policy: BundleBranchPolicy<'_>,
    ) -> Result<ConnectedIrradianceSpectrum, ConnectedIrradianceError> {
        connect_spectrum(self.intensity_source(), cx, scene, request, policy, |bin| {
            let state = self.relative_stokes()[bin];
            if state.i == 0.0 {
                Stokes::unpolarized(1.0)
            } else {
                Stokes {
                    i: 1.0,
                    q: state.q / state.i,
                    u: state.u / state.i,
                    v: state.v / state.i,
                }
            }
        })
    }
}

fn connect_spectrum(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    request: ConnectedIrradianceRequest,
    policy: BundleBranchPolicy<'_>,
    polarization: impl Fn(usize) -> Stokes,
) -> Result<ConnectedIrradianceSpectrum, ConnectedIrradianceError> {
    let target = BundleConnectionTarget::try_new(
        request.detector_object_id,
        request.receiver.position_m(),
        request.receiver.normal(),
        request.receiver.position_tolerance_m(),
    )
    .map_err(ConnectedIrradianceError::InvalidTarget)?;
    let mut meter =
        BundleIrradianceMeter::new(source, request.receiver, request.connection.projection);
    let mut bins_w_m2_nm = [0.0; SPECTRAL_BINS];
    let mut launch_directions = [Vec3::new(0.0, 0.0, 0.0); SPECTRAL_BINS];
    let mut direction = request.initial_direction;
    let mut signature = None;
    let mut max_endpoint_error_m = 0.0_f64;
    let mut max_contact_resolution_m = 0.0_f64;
    let mut trace_evaluations = 0_u64;
    let mut field_evaluations = 0_u64;
    let mut query_work = 0_u64;
    for bin in 0..SPECTRAL_BINS {
        let connected = scene
            .connect_bundle_with_policy(
                cx,
                BundleConnectionLaunch {
                    ray: Ray {
                        origin: source.position_m(),
                        dir: direction,
                    },
                    transverse_axis: request.source_transverse_axis,
                    wavelength_nm: wavelength_nm(bin),
                    normalized_source_stokes: polarization(bin),
                    initial_medium: request.initial_medium,
                },
                target,
                request.connection,
                policy,
            )
            .map_err(|error| ConnectedIrradianceError::Connection { bin, error })?;
        let path = connected.path();
        if signature.is_some_and(|previous| previous != path.branch_signature()) {
            return Err(ConnectedIrradianceError::SpectralBranchChanged { bin });
        }
        signature = Some(path.branch_signature());
        meter
            .record_refractive(bin, path.bundle())
            .map_err(|error| ConnectedIrradianceError::Photometry { bin, error })?;
        bins_w_m2_nm[bin] = meter
            .spectral_irradiance_w_m2_nm(bin)
            .ok_or(ConnectedIrradianceError::NumericalRange)?;
        direction = path.bundle().propagation().initial_ray().dir;
        launch_directions[bin] = direction;
        max_endpoint_error_m = max_endpoint_error_m.max(connected.residual_m());
        max_contact_resolution_m = max_contact_resolution_m.max(path.max_contact_resolution_m());
        trace_evaluations = trace_evaluations
            .checked_add(u64::from(connected.trace_evaluations()))
            .ok_or(ConnectedIrradianceError::NumericalRange)?;
        field_evaluations = field_evaluations
            .checked_add(connected.field_evaluations())
            .ok_or(ConnectedIrradianceError::NumericalRange)?;
        query_work = query_work
            .checked_add(connected.query_work())
            .ok_or(ConnectedIrradianceError::NumericalRange)?;
    }
    let reading = meter
        .reading()
        .map_err(|error| ConnectedIrradianceError::Photometry {
            bin: SPECTRAL_BINS - 1,
            error,
        })?;
    Ok(ConnectedIrradianceSpectrum {
        bins_w_m2_nm,
        launch_directions,
        reading,
        branch_signature: signature.ok_or(ConnectedIrradianceError::NumericalRange)?,
        max_endpoint_error_m,
        max_contact_resolution_m,
        trace_evaluations,
        field_evaluations,
        query_work,
    })
}

#[cfg(test)]
mod tests;
