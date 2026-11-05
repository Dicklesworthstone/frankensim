//! Bounded incoherent sums of the original slab's transmitted/reflected images.

use super::*;

/// Truncation tolerances and hard work limits for one complete 80-bin reading.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabHistoryOptions {
    /// Absolute omitted irradiance allowed in each bin, in W/(m² nm).
    pub absolute_tolerance_w_m2_nm: f64,
    /// Additional allowed omitted irradiance as a fraction of that bin's sum.
    pub relative_tolerance: f64,
    /// Maximum evaluated images per wavelength, including the direct image
    /// (transmission) or front-surface mirror (reflection).
    pub max_images_per_wavelength: u32,
    /// Total reserved boundary interactions across all wavelengths. Order k
    /// reserves 2k+2 for transmission or 2k+1 for reflection before tracing,
    /// even if opacity or back-facing reception makes it finish early. The
    /// receiver visibility query is additional.
    pub max_interface_interactions: u64,
}

impl Default for SlabHistoryOptions {
    fn default() -> Self {
        Self {
            absolute_tolerance_w_m2_nm: 1.0e-12,
            relative_tolerance: 1.0e-7,
            max_images_per_wavelength: 128,
            max_interface_interactions: 1_048_576,
        }
    }
}

impl SlabHistoryOptions {
    fn validate(self) -> Result<(), SlabHistoryError> {
        if !self.absolute_tolerance_w_m2_nm.is_finite()
            || self.absolute_tolerance_w_m2_nm < 0.0
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || (self.absolute_tolerance_w_m2_nm == 0.0 && self.relative_tolerance == 0.0)
            || self.max_images_per_wavelength == 0
            || self.max_interface_interactions == 0
        {
            return Err(SlabHistoryError::InvalidOptions);
        }
        Ok(())
    }
}

/// A refusal never publishes an unconverged partial sum as a complete reading.
#[derive(Debug)]
pub enum SlabHistoryError {
    /// Invalid tolerance or zero work limit.
    InvalidOptions,
    /// Original source, medium, geometry, distance or cancellation refusal.
    Transport(PhotometryError),
    /// The next omitted-image bound still exceeds the requested tolerance.
    ImageLimitExceeded {
        /// Native wavelength bin that could not converge.
        bin: usize,
        /// Images already evaluated in this bin.
        images: u32,
        /// Current analytic omitted-image bound, in W/(m² nm).
        remaining_bound_w_m2_nm: f64,
    },
    /// The next complete history would exceed the total boundary-work cap.
    InterfaceBudgetExceeded {
        /// Native wavelength bin at the work limit.
        bin: usize,
        /// Number of internal round trips in the refused next image.
        round_trips: u32,
        /// Boundary interactions needed for that image.
        required: u64,
        /// Unreserved boundary interactions remaining in this job.
        remaining: u64,
    },
}

impl From<PhotometryError> for SlabHistoryError {
    fn from(error: PhotometryError) -> Self {
        Self::Transport(error)
    }
}

impl core::fmt::Display for SlabHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "slab optical history refused: {self:?}")
    }
}

impl std::error::Error for SlabHistoryError {}

/// Which broad-face family a complete reading contains.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlabHistoryFamily {
    /// T(RR)^kT across the slab, k >= 0.
    Transmitted,
    /// The front mirror R, followed by T R^(2k-1) T for k >= 1.
    Reflected,
}

/// Incoherent broad-face irradiance plus its omitted-history bound.
/// Bounds concern image truncation in the declared geometric-optics model;
/// they do not enclose floating-point, mesh, spectrum or calibration errors.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabHistorySpectrum {
    family: SlabHistoryFamily,
    bins_w_m2_nm: [f64; SPECTRAL_BINS],
    omitted_w_m2_nm: [f64; SPECTRAL_BINS],
    images: [u32; SPECTRAL_BINS],
    reserved_interfaces: u64,
    reading: BundleIrradianceReading,
    omitted_reading: BundleIrradianceReading,
}

/// Spectrum returned by the opposite-side transmission measurement.
pub type TransmittedSlabHistorySpectrum = SlabHistorySpectrum;

/// Spectrum returned by the same-side reflection measurement.
pub type ReflectedSlabHistorySpectrum = SlabHistorySpectrum;

impl SlabHistorySpectrum {
    /// Selected transmitted or reflected broad-face family.
    pub const fn family(&self) -> SlabHistoryFamily {
        self.family
    }

    /// Sum of evaluated images at native wavelength midpoints, in W/(m² nm).
    pub const fn spectral_irradiance_w_m2_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.bins_w_m2_nm
    }

    /// Nonnegative analytic upper bound on all unevaluated images in each bin.
    pub const fn omitted_spectral_irradiance_bound_w_m2_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.omitted_w_m2_nm
    }

    /// Evaluated images per bin; zero is possible for an absent or weak source.
    pub const fn images_per_wavelength(&self) -> &[u32; SPECTRAL_BINS] {
        &self.images
    }

    /// Conservative interface work reserved, including early terminated images.
    pub const fn reserved_interface_interactions(&self) -> u64 {
        self.reserved_interfaces
    }

    /// Visible radiant irradiance and lux in the evaluated selected histories.
    pub const fn reading(&self) -> BundleIrradianceReading {
        self.reading
    }

    /// Omitted radiant irradiance and lux bounds under the same native quadrature.
    pub const fn omitted_reading_bound(&self) -> BundleIrradianceReading {
        self.omitted_reading
    }
}

impl PhotometricSource {
    /// Sum incoherent transmitted images with zero, one, two, ... internal
    /// round trips through a named native smooth parallel glass slab.
    ///
    /// Each order independently solves its dispersive launch direction and
    /// receiving-area Jacobian, evaluates this source's actual angular law,
    /// and traces every folded physical segment through the original scene.
    /// Native s/p reflection, transmission, coherent surface coatings and
    /// actual internal Beer distance survive the complete history. Opacity,
    /// a zero source intensity or back-facing/grazing reception zeros only
    /// that order; later images are still considered.
    ///
    /// Completion requires each omitted-image bound to be no greater than
    /// `absolute_tolerance_w_m2_nm + relative_tolerance * accumulated_bin`.
    /// Bare glass uses a tight Fresnel geometric-series bound. Passive coated
    /// glass also works, using a more conservative inverse-square spreading
    /// bound that can need a larger budget or looser tolerance. Bounds use the
    /// source's global angular peak, so a dark/blocked order cannot imply that
    /// all later orders are dark. No sampled angular maximum is assumed.
    ///
    /// Geometry and media admission match the selected-image owner, including
    /// closed static mesh ownership and vacuum endpoints. Nonparallel/edge
    /// encounters, extra dielectrics and unresolved geometry refuse. Distance
    /// limits apply to each evaluated physical path; both image count and total
    /// interface work are bounded. Cancellation publishes no partial spectrum.
    /// This is the broad-face T(RR)^mT family, not total scene illumination or
    /// coherent interference between bulk faces. Existing general path-job
    /// point-source/glass restrictions remain in effect.
    pub fn transmitted_slab_history_irradiance(
        &self,
        cx: &Cx<'_>,
        scene: &Scene,
        slab_primitive: usize,
        receiver: IlluminanceProbe,
        transport: TransportOptions,
        history: SlabHistoryOptions,
    ) -> Result<TransmittedSlabHistorySpectrum, SlabHistoryError> {
        history.validate()?;
        let connection = SlabConnection::try_new(
            cx,
            scene,
            slab_primitive,
            self.position_m(),
            receiver,
            transport,
        )?;
        self.slab_history_irradiance(cx, &connection, receiver, transport, history)
    }

    /// Measure the front-surface mirror plus every incoherent internal return
    /// through the same broad face: R + T R (RR)^m T, m >= 0.
    ///
    /// The source and receiver are in vacuum on the same side of the selected
    /// static, closed, smooth parallel slab. `front_face_point_m` is a point
    /// on its original source-facing broad face. This witness selects the face
    /// pair and is validated against the original body; it does not prescribe
    /// a physical reflection point. Opaque geometry on that discovery segment
    /// does not occlude light. Each solved physical segment sees the complete
    /// scene. Colocated source and receiver positions are valid.
    ///
    /// Each wavelength and image independently determines its source direction,
    /// reflected receiver projection and area Jacobian. The native surface
    /// response retains s/p power, coatings seen from the actual side, and Beer
    /// loss over every folded internal segment. Order zero is the front mirror;
    /// order k >= 1 crosses 2k slab widths and has 2k-1 internal reflections.
    /// Its reserved work is 2k+1 boundary interactions. Dark, blocked or
    /// back-facing early images do not suppress later orders.
    ///
    /// The same per-bin truncation tolerances, image/interaction caps and
    /// complete-spectrum cancellation semantics as transmitted histories apply.
    /// The output excludes direct source-to-receiver light, side-face paths,
    /// other scene illumination and coherent interference between bulk faces.
    /// It does not relax general path jobs' point-source/glass admission.
    #[allow(clippy::too_many_arguments)]
    pub fn reflected_slab_history_irradiance(
        &self,
        cx: &Cx<'_>,
        scene: &Scene,
        slab_primitive: usize,
        front_face_point_m: fs_geom::Point3,
        receiver: IlluminanceProbe,
        transport: TransportOptions,
        history: SlabHistoryOptions,
    ) -> Result<ReflectedSlabHistorySpectrum, SlabHistoryError> {
        history.validate()?;
        let connection = SlabConnection::try_reflected(
            cx,
            scene,
            slab_primitive,
            self.position_m(),
            front_face_point_m,
            receiver,
            transport,
        )?;
        self.slab_history_irradiance(cx, &connection, receiver, transport, history)
    }

    fn slab_history_irradiance(
        &self,
        cx: &Cx<'_>,
        connection: &SlabConnection<'_>,
        receiver: IlluminanceProbe,
        transport: TransportOptions,
        history: SlabHistoryOptions,
    ) -> Result<SlabHistorySpectrum, SlabHistoryError> {
        let plane = PlanarIrradianceReceiver::try_new(
            receiver.position_m(),
            receiver.normal(),
            transport.ray_epsilon_m,
        )
        .map_err(|_| PhotometryError::InvalidProbe)?;
        let mut meter = BundleIrradianceMeter::new(self, plane, ProjectionOptions::default());
        let mut tail_meter = BundleIrradianceMeter::new(self, plane, ProjectionOptions::default());
        let zero = BundleIrradianceReading {
            radiant_irradiance_w_m2: 0.0,
            illuminance_lux: 0.0,
        };
        let mut result = SlabHistorySpectrum {
            family: if connection.is_reflected() {
                SlabHistoryFamily::Reflected
            } else {
                SlabHistoryFamily::Transmitted
            },
            bins_w_m2_nm: [0.0; SPECTRAL_BINS],
            omitted_w_m2_nm: [0.0; SPECTRAL_BINS],
            images: [0; SPECTRAL_BINS],
            reserved_interfaces: 0,
            reading: zero,
            omitted_reading: zero,
        };
        let peak_cd = self.emitted_flux_lm * self.distribution.peak_per_steradian();
        for bin in 0..SPECTRAL_BINS {
            cx.checkpoint().map_err(PhotometryError::from)?;
            let wavelength = wavelength_nm(bin);
            let peak = peak_cd * self.spectrum.value_w_per_lumen_nm(wavelength);
            let mut sum = 0.0;
            let mut images = 0;
            let bound = loop {
                cx.checkpoint().map_err(PhotometryError::from)?;
                if peak == 0.0 {
                    break 0.0;
                }
                // This solve is reused if the next image needs evaluation.
                let image = connection.image(cx, wavelength, images)?;
                let bound = connection.omitted_bound(cx, &image, peak)?;
                if bound <= history.absolute_tolerance_w_m2_nm + history.relative_tolerance * sum {
                    break bound;
                }
                if images >= history.max_images_per_wavelength {
                    return Err(SlabHistoryError::ImageLimitExceeded {
                        bin,
                        images,
                        remaining_bound_w_m2_nm: bound,
                    });
                }
                let required = image.interface_interactions(connection.is_reflected());
                let remaining = history.max_interface_interactions - result.reserved_interfaces;
                if required > remaining {
                    return Err(SlabHistoryError::InterfaceBudgetExceeded {
                        bin,
                        round_trips: images,
                        required,
                        remaining,
                    });
                }
                result.reserved_interfaces += required;
                let ray = connection.sample_image(cx, &image)?;
                let intensity = self
                    .spectral_intensity_w_sr_nm(ray.launch_direction, wavelength)
                    .map_err(PhotometryError::from)?;
                let value = ray.irradiance(intensity)?;
                sum += value;
                if !sum.is_finite() || value < 0.0 {
                    return Err(PhotometryError::NonFiniteResult.into());
                }
                images += 1;
            };
            meter
                .record_solved_irradiance(bin, sum)
                .map_err(|_| PhotometryError::NonFiniteResult)?;
            tail_meter
                .record_solved_irradiance(bin, bound)
                .map_err(|_| PhotometryError::NonFiniteResult)?;
            result.bins_w_m2_nm[bin] = sum;
            result.omitted_w_m2_nm[bin] = bound;
            result.images[bin] = images;
        }
        result.reading = meter
            .reading()
            .map_err(|_| PhotometryError::NonFiniteResult)?;
        result.omitted_reading = tail_meter
            .reading()
            .map_err(|_| PhotometryError::NonFiniteResult)?;
        cx.checkpoint().map_err(PhotometryError::from)?;
        Ok(result)
    }
}
