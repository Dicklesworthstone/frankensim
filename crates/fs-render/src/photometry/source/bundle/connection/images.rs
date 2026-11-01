//! Calibrated spectral power from a bounded search for regular optical images.
//!
//! Several rays can reach one receiver after exactly the same ordered surface
//! interactions. The native image search distinguishes their launch directions,
//! and this owner continues each discovered image through the 80 wavelength
//! bins. Every contribution retains its own source intensity, polarized Fresnel
//! transfer, Beer absorption and receiving-area Jacobian. Powers add
//! incoherently; no extra cosine or refractive-index multiplier is introduced.
//!
//! Discovery uses the first spectral bin and the caller's finite seed list.
//! A missed seed is not evidence that its neighborhood contains no image, and
//! images born at other wavelengths are not enumerated. A successful reading
//! therefore describes the discovered, spectrally continued image set, not all
//! light reaching the receiver. The native finite-detector estimator remains
//! the appropriate owner for an unresolved caustic or aperture-integrated flux.
//!
//! Continuation refuses history changes, footprint parity changes, numerical
//! image collisions and excessive angular movement. Movement is bounded both
//! by the caller's limit and by one quarter of the nearest previous same-history
//! image separation. These local checks guard against jumping to another image;
//! they are not an interval proof of uniqueness between wavelength samples.

use super::super::{BundleIrradianceMeter, BundleIrradianceReading, Sum};
use super::{
    BundlePhotometryError, ConnectedIrradianceSpectrum, CurvedPolarizedScene, PhotometricSource,
    PlanarIrradianceReceiver, PolarizedSource, SPECTRAL_BINS, Stokes, wavelength_nm,
};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleConnection, BundleConnectionError, BundleConnectionLaunch, BundleConnectionOptions,
    BundleConnectionTarget,
    images::{
        BundleImageCatalog, BundleImageError, BundleImageOptions, BundleImageSeed,
        ImagePairClassification, MAX_BUNDLE_IMAGES, angular_distance, classify_image_pair,
        connection_angular_resolution,
    },
};
use fs_exec::Cx;
use fs_geom::Vec3;

/// Shared receiver, source frame and bounded numerical controls for image
/// discovery. Every seed supplies its own direction and optical branch policy.
/// `max_continuation_step_rad` must be finite, positive and at most one radian.
#[derive(Clone, Copy, Debug)]
pub struct ImageIrradianceRequest {
    pub detector_object_id: u64,
    pub receiver: PlanarIrradianceReceiver,
    pub source_transverse_axis: Vec3,
    pub initial_medium: usize,
    pub connection: BundleConnectionOptions,
    pub discovery: BundleImageOptions,
    pub max_continuation_step_rad: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ImageIrradianceError {
    InvalidContinuation,
    InvalidTarget(BundleConnectionError),
    Discovery(BundleImageError),
    /// No admitted image was found. This is not a zero-light measurement.
    NoImages,
    Connection {
        image: usize,
        bin: usize,
        error: BundleConnectionError,
    },
    Photometry {
        image: usize,
        bin: usize,
        error: BundlePhotometryError,
    },
    SumPhotometry(BundlePhotometryError),
    SpectralBranchChanged {
        image: usize,
        bin: usize,
    },
    SpectralParityChanged {
        image: usize,
        bin: usize,
    },
    ImageCollision {
        first: usize,
        second: usize,
        bin: usize,
    },
    AmbiguousImages {
        first: usize,
        second: usize,
        bin: usize,
    },
    ContinuationJump {
        image: usize,
        bin: usize,
        angle_rad: f64,
        allowed_rad: f64,
    },
    AllocationRefused,
    NumericalRange,
}

impl core::fmt::Display for ImageIrradianceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "discovered-image spectral irradiance: {self:?}")
    }
}
impl core::error::Error for ImageIrradianceError {}

/// Complete incoherent sum of the discovered images. Each member retains its
/// 80 launch directions and spectral contributions. The discovery catalog has
/// the same image ordering and reports every seed, including misses/duplicates.
///
/// Total work includes all discovery trials and spectral continuation. A member
/// spectrum's work describes its retained first-bin connection and continuation;
/// discarded discovery trials belong to the catalog and total, not that member.
#[derive(Clone, Debug, PartialEq)]
pub struct DiscoveredIrradianceSpectrum {
    bins_w_m2_nm: [f64; SPECTRAL_BINS],
    reading: BundleIrradianceReading,
    images: Vec<ConnectedIrradianceSpectrum>,
    discovery: BundleImageCatalog,
    max_endpoint_error_m: f64,
    max_contact_resolution_m: f64,
    trace_evaluations: u64,
    field_evaluations: u64,
    query_work: u64,
}

impl DiscoveredIrradianceSpectrum {
    pub const fn spectral_irradiance_w_m2_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.bins_w_m2_nm
    }
    pub const fn reading(&self) -> BundleIrradianceReading {
        self.reading
    }
    pub fn images(&self) -> &[ConnectedIrradianceSpectrum] {
        &self.images
    }
    pub const fn discovery(&self) -> &BundleImageCatalog {
        &self.discovery
    }
    pub const fn max_endpoint_error_m(&self) -> f64 {
        self.max_endpoint_error_m
    }
    /// Largest native contact collar; not a bound on physical irradiance error.
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
    /// Discover regular images at the first native wavelength and continue
    /// their calibrated powers across the complete visible quadrature. Distinct
    /// images may share the same surface/medium/reflection/transmission history.
    pub fn discover_irradiance_images(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ImageIrradianceRequest,
        seeds: &[BundleImageSeed<'_>],
    ) -> Result<DiscoveredIrradianceSpectrum, ImageIrradianceError> {
        discover_spectrum(self, cx, scene, request, seeds, |_| {
            Stokes::unpolarized(1.0)
        })
    }
}

impl PolarizedSource {
    /// Preserve the surviving source polarization per wavelength. The intensity
    /// source already contains upstream filtering; normalizing Stokes here does
    /// not restore lost power or apply that filtering a second time.
    pub fn discover_irradiance_images(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: ImageIrradianceRequest,
        seeds: &[BundleImageSeed<'_>],
    ) -> Result<DiscoveredIrradianceSpectrum, ImageIrradianceError> {
        discover_spectrum(self.intensity_source(), cx, scene, request, seeds, |bin| {
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

fn discover_spectrum(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    request: ImageIrradianceRequest,
    seeds: &[BundleImageSeed<'_>],
    polarization: impl Fn(usize) -> Stokes,
) -> Result<DiscoveredIrradianceSpectrum, ImageIrradianceError> {
    if !request.max_continuation_step_rad.is_finite()
        || request.max_continuation_step_rad <= 0.0
        || request.max_continuation_step_rad > 1.0
    {
        return Err(ImageIrradianceError::InvalidContinuation);
    }
    let target = BundleConnectionTarget::try_new(
        request.detector_object_id,
        request.receiver.position_m(),
        request.receiver.normal(),
        request.receiver.position_tolerance_m(),
    )
    .map_err(ImageIrradianceError::InvalidTarget)?;
    // The native discovery owner admits the seed count and substitutes each
    // seed's direction. The placeholder direction is never traced.
    let launch = BundleConnectionLaunch {
        ray: Ray {
            origin: source.position_m(),
            dir: Vec3::new(0.0, 0.0, 1.0),
        },
        transverse_axis: request.source_transverse_axis,
        wavelength_nm: wavelength_nm(0),
        normalized_source_stokes: polarization(0),
        initial_medium: request.initial_medium,
    };
    let discovery = scene
        .discover_bundle_images(
            cx,
            launch,
            target,
            request.connection,
            seeds,
            request.discovery,
        )
        .map_err(ImageIrradianceError::Discovery)?;
    if discovery.images().is_empty() {
        return Err(ImageIrradianceError::NoImages);
    }
    let mut images = Vec::new();
    images
        .try_reserve_exact(discovery.images().len())
        .map_err(|_| ImageIrradianceError::AllocationRefused)?;
    let mut resolutions = [0.0; MAX_BUNDLE_IMAGES];
    let mut parities = [false; MAX_BUNDLE_IMAGES];
    let mut trace_evaluations = u64::from(discovery.trace_evaluations());
    let mut field_evaluations = discovery.field_evaluations();
    let mut query_work = discovery.query_work();
    for (image, member) in discovery.images().iter().enumerate() {
        let connected = member.connection();
        let mut spectrum = ConnectedIrradianceSpectrum {
            bins_w_m2_nm: [0.0; SPECTRAL_BINS],
            launch_directions: [Vec3::new(0.0, 0.0, 0.0); SPECTRAL_BINS],
            reading: BundleIrradianceReading {
                radiant_irradiance_w_m2: 0.0,
                illuminance_lux: 0.0,
            },
            branch_signature: connected.path().branch_signature(),
            max_endpoint_error_m: 0.0,
            max_contact_resolution_m: 0.0,
            trace_evaluations: 0,
            field_evaluations: 0,
            query_work: 0,
        };
        record_connection(source, request, &mut spectrum, connected, image, 0)?;
        resolutions[image] = member.angular_resolution_rad();
        parities[image] = footprint_parity(connected, request, image, 0)?;
        images.push(spectrum);
    }
    for bin in 1..SPECTRAL_BINS {
        for image in 0..images.len() {
            cx.checkpoint()
                .map_err(|_| ImageIrradianceError::Connection {
                    image,
                    bin,
                    error: BundleConnectionError::Cancelled,
                })?;
            let previous = images[image].launch_directions[bin - 1];
            let connected = scene
                .connect_bundle_with_policy(
                    cx,
                    BundleConnectionLaunch {
                        ray: Ray {
                            origin: source.position_m(),
                            dir: previous,
                        },
                        wavelength_nm: wavelength_nm(bin),
                        normalized_source_stokes: polarization(bin),
                        ..launch
                    },
                    target,
                    request.connection,
                    seeds[discovery.images()[image].first_seed()].policy,
                )
                .map_err(|error| ImageIrradianceError::Connection { image, bin, error })?;
            if connected.path().branch_signature() != images[image].branch_signature {
                return Err(ImageIrradianceError::SpectralBranchChanged { image, bin });
            }
            if footprint_parity(&connected, request, image, bin)? != parities[image] {
                return Err(ImageIrradianceError::SpectralParityChanged { image, bin });
            }
            let direction = connected.path().bundle().propagation().initial_ray().dir;
            let resolution =
                connection_angular_resolution(&connected, target, request.connection.projection)
                    .map_err(|error| ImageIrradianceError::Connection { image, bin, error })?;
            let angle = angular_distance(previous, direction);
            let mut allowed = request.max_continuation_step_rad;
            for (other, spectrum) in images.iter().enumerate() {
                if other != image && spectrum.branch_signature == images[image].branch_signature {
                    allowed = allowed.min(
                        0.25 * angular_distance(previous, spectrum.launch_directions[bin - 1]),
                    );
                }
            }
            if !angle.is_finite() || angle + resolutions[image] + resolution > allowed {
                return Err(ImageIrradianceError::ContinuationJump {
                    image,
                    bin,
                    angle_rad: angle,
                    allowed_rad: allowed,
                });
            }
            for other in 0..image {
                if images[other].branch_signature != images[image].branch_signature {
                    continue;
                }
                let separation = angular_distance(direction, images[other].launch_directions[bin]);
                match classify_image_pair(
                    request.discovery,
                    separation,
                    resolution,
                    resolutions[other],
                ) {
                    ImagePairClassification::Distinct => {}
                    ImagePairClassification::Duplicate => {
                        return Err(ImageIrradianceError::ImageCollision {
                            first: other,
                            second: image,
                            bin,
                        });
                    }
                    ImagePairClassification::Ambiguous => {
                        return Err(ImageIrradianceError::AmbiguousImages {
                            first: other,
                            second: image,
                            bin,
                        });
                    }
                }
            }
            record_connection(source, request, &mut images[image], &connected, image, bin)?;
            resolutions[image] = resolution;
            trace_evaluations = trace_evaluations
                .checked_add(u64::from(connected.trace_evaluations()))
                .ok_or(ImageIrradianceError::NumericalRange)?;
            field_evaluations = field_evaluations
                .checked_add(connected.field_evaluations())
                .ok_or(ImageIrradianceError::NumericalRange)?;
            query_work = query_work
                .checked_add(connected.query_work())
                .ok_or(ImageIrradianceError::NumericalRange)?;
        }
    }
    let mut totals: [Sum; SPECTRAL_BINS] = core::array::from_fn(|_| Sum::default());
    let mut max_endpoint_error_m = 0.0_f64;
    let mut max_contact_resolution_m = 0.0_f64;
    for (image, spectrum) in images.iter_mut().enumerate() {
        spectrum.reading = integrate(source, request, &spectrum.bins_w_m2_nm).map_err(|error| {
            ImageIrradianceError::Photometry {
                image,
                bin: SPECTRAL_BINS - 1,
                error,
            }
        })?;
        for (total, value) in totals.iter_mut().zip(spectrum.bins_w_m2_nm) {
            total.add(value);
        }
        max_endpoint_error_m = max_endpoint_error_m.max(spectrum.max_endpoint_error_m);
        max_contact_resolution_m = max_contact_resolution_m.max(spectrum.max_contact_resolution_m);
    }
    let bins_w_m2_nm = totals.map(Sum::value);
    let reading =
        integrate(source, request, &bins_w_m2_nm).map_err(ImageIrradianceError::SumPhotometry)?;
    Ok(DiscoveredIrradianceSpectrum {
        bins_w_m2_nm,
        reading,
        images,
        discovery,
        max_endpoint_error_m,
        max_contact_resolution_m,
        trace_evaluations,
        field_evaluations,
        query_work,
    })
}

fn record_connection(
    source: &PhotometricSource,
    request: ImageIrradianceRequest,
    spectrum: &mut ConnectedIrradianceSpectrum,
    connected: &BundleConnection,
    image: usize,
    bin: usize,
) -> Result<(), ImageIrradianceError> {
    let path = connected.path();
    spectrum.bins_w_m2_nm[bin] = source
        .spectral_irradiance_from_refractive_bundle(
            path.bundle(),
            request.receiver,
            request.connection.projection,
        )
        .map_err(|error| ImageIrradianceError::Photometry { image, bin, error })?;
    spectrum.launch_directions[bin] = path.bundle().propagation().initial_ray().dir;
    spectrum.max_endpoint_error_m = spectrum.max_endpoint_error_m.max(connected.residual_m());
    spectrum.max_contact_resolution_m = spectrum
        .max_contact_resolution_m
        .max(path.max_contact_resolution_m());
    spectrum.trace_evaluations = spectrum
        .trace_evaluations
        .checked_add(u64::from(connected.trace_evaluations()))
        .ok_or(ImageIrradianceError::NumericalRange)?;
    spectrum.field_evaluations = spectrum
        .field_evaluations
        .checked_add(connected.field_evaluations())
        .ok_or(ImageIrradianceError::NumericalRange)?;
    spectrum.query_work = spectrum
        .query_work
        .checked_add(connected.query_work())
        .ok_or(ImageIrradianceError::NumericalRange)?;
    Ok(())
}

fn footprint_parity(
    connected: &BundleConnection,
    request: ImageIrradianceRequest,
    image: usize,
    bin: usize,
) -> Result<bool, ImageIrradianceError> {
    connected
        .path()
        .bundle()
        .propagation()
        .project_at_endpoint(request.receiver.normal(), request.connection.projection)
        .map(|projection| {
            projection
                .signed_area_per_solid_angle_m2()
                .is_sign_positive()
        })
        .map_err(|error| ImageIrradianceError::Connection {
            image,
            bin,
            error: BundleConnectionError::Projection(error),
        })
}

fn integrate(
    source: &PhotometricSource,
    request: ImageIrradianceRequest,
    values: &[f64; SPECTRAL_BINS],
) -> Result<BundleIrradianceReading, BundlePhotometryError> {
    let mut meter =
        BundleIrradianceMeter::new(source, request.receiver, request.connection.projection);
    for (bin, value) in values.iter().copied().enumerate() {
        meter.record_value(bin, 1_u128 << bin, value)?;
    }
    meter.reading()
}

#[cfg(test)]
mod tests;
