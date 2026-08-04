//! Deterministic shutter-time semantics and rays carrying explicit time.

use core::fmt;

/// Stable counter-domain separation for shutter samples.
pub const SHUTTER_TIME_SAMPLE_DOMAIN_V1: u64 = 0x6d6f_7469_6f6e_0001;

/// Placement of an exposure interval relative to a nominal frame time.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ShutterConvention {
    /// Exposure is centered on the nominal frame time.
    Centered,
    /// Exposure starts at the nominal frame time.
    FrontLoaded,
    /// Exposure ends at the nominal frame time.
    BackLoaded,
}

/// Deterministic distribution used to choose a time inside an exposure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ShutterDistribution {
    /// Named counter-separated uniform distribution, version 1.
    UniformCounterV1,
}

/// Shot-level interval within which frame shutters must remain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShotTimeBounds {
    /// Inclusive shot start [s].
    start_s: f64,
    /// Inclusive shot end [s].
    end_s: f64,
}

impl ShotTimeBounds {
    /// Admit finite, ordered shot bounds. A zero-duration shot is valid.
    pub fn try_new(start_s: f64, end_s: f64) -> Result<Self, MotionTimeError> {
        if !start_s.is_finite() || !end_s.is_finite() || start_s > end_s {
            return Err(MotionTimeError::InvalidShotBounds);
        }
        Ok(Self { start_s, end_s })
    }
}

/// Resolved absolute frame exposure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShutterInterval {
    /// Inclusive exposure start [s].
    open_s: f64,
    /// Inclusive exposure end [s].
    close_s: f64,
    /// Convention that produced the interval.
    convention: ShutterConvention,
    /// Deterministic sampling distribution.
    distribution: ShutterDistribution,
}

impl ShutterInterval {
    /// Resolve frame-relative shutter semantics inside declared shot bounds.
    pub fn resolve(
        frame_time_s: f64,
        exposure_duration_s: f64,
        convention: ShutterConvention,
        distribution: ShutterDistribution,
        shot: ShotTimeBounds,
    ) -> Result<Self, MotionTimeError> {
        if !frame_time_s.is_finite()
            || !exposure_duration_s.is_finite()
            || exposure_duration_s < 0.0
        {
            return Err(MotionTimeError::InvalidExposure);
        }
        let (open_s, close_s) = match convention {
            ShutterConvention::Centered => {
                let half = 0.5 * exposure_duration_s;
                (frame_time_s - half, frame_time_s + half)
            }
            ShutterConvention::FrontLoaded => (frame_time_s, frame_time_s + exposure_duration_s),
            ShutterConvention::BackLoaded => (frame_time_s - exposure_duration_s, frame_time_s),
        };
        if !open_s.is_finite()
            || !close_s.is_finite()
            || open_s < shot.start_s
            || close_s > shot.end_s
        {
            return Err(MotionTimeError::ExposureOutsideShot);
        }
        Ok(Self {
            open_s,
            close_s,
            convention,
            distribution,
        })
    }

    /// Exposure duration [s].
    #[must_use]
    pub fn duration_s(self) -> f64 {
        self.close_s - self.open_s
    }

    /// Inclusive exposure start [s].
    #[must_use]
    pub const fn open_s(self) -> f64 {
        self.open_s
    }

    /// Inclusive exposure end [s].
    #[must_use]
    pub const fn close_s(self) -> f64 {
        self.close_s
    }

    /// Frame-relative shutter convention.
    #[must_use]
    pub const fn convention(self) -> ShutterConvention {
        self.convention
    }

    /// Deterministic time-sampling distribution.
    #[must_use]
    pub const fn distribution(self) -> ShutterDistribution {
        self.distribution
    }

    /// Map an admitted normalized coordinate to absolute time. A zero-width
    /// shutter always returns its one exact endpoint.
    #[must_use]
    pub fn time_at(self, normalized_time: NormalizedShutterTime) -> f64 {
        self.duration_s()
            .mul_add(normalized_time.value(), self.open_s)
    }

    /// Deterministically sample a shutter coordinate from stable logical IDs.
    #[must_use]
    pub fn sample(self, pixel_identity: u64, sample_identity: u64) -> NormalizedShutterTime {
        match self.distribution {
            ShutterDistribution::UniformCounterV1 => {
                let counter = pixel_identity
                    .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                    .wrapping_add(sample_identity)
                    ^ SHUTTER_TIME_SAMPLE_DOMAIN_V1;
                NormalizedShutterTime::from_counter_bits(splitmix64(counter))
            }
        }
    }
}

/// Finite normalized shutter coordinate in `[0, 1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalizedShutterTime(f64);

impl NormalizedShutterTime {
    /// Admit an explicit normalized coordinate.
    pub fn try_new(value: f64) -> Result<Self, MotionTimeError> {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(MotionTimeError::InvalidNormalizedTime);
        }
        Ok(Self(value))
    }

    fn from_counter_bits(bits: u64) -> Self {
        let value = (bits >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0);
        Self(value)
    }

    /// Coordinate value.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0
    }
}

/// Any spatial ray paired with normalized and absolute shutter time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TimedRay<SpatialRay> {
    /// Existing backend-specific spatial ray.
    spatial: SpatialRay,
    /// Dimensionless exposure coordinate.
    normalized_time: NormalizedShutterTime,
    /// Absolute trajectory/camera time [s].
    absolute_time_s: f64,
}

impl<SpatialRay> TimedRay<SpatialRay> {
    /// Bind an existing spatial ray to one deterministic shutter sample.
    #[must_use]
    pub fn from_sample(
        spatial: SpatialRay,
        shutter: ShutterInterval,
        pixel_identity: u64,
        sample_identity: u64,
    ) -> Self {
        let normalized_time = shutter.sample(pixel_identity, sample_identity);
        Self {
            spatial,
            normalized_time,
            absolute_time_s: shutter.time_at(normalized_time),
        }
    }

    /// Bind a ray to an explicit admitted coordinate, useful for endpoints and
    /// deterministic reference accumulation.
    #[must_use]
    pub fn at_normalized(
        spatial: SpatialRay,
        shutter: ShutterInterval,
        normalized_time: NormalizedShutterTime,
    ) -> Self {
        Self {
            spatial,
            normalized_time,
            absolute_time_s: shutter.time_at(normalized_time),
        }
    }

    /// Borrow the backend-specific spatial ray.
    #[must_use]
    pub const fn spatial(&self) -> &SpatialRay {
        &self.spatial
    }

    /// Dimensionless coordinate inside the resolved exposure.
    #[must_use]
    pub const fn normalized_time(&self) -> NormalizedShutterTime {
        self.normalized_time
    }

    /// Absolute trajectory/camera time [s].
    #[must_use]
    pub const fn absolute_time_s(&self) -> f64 {
        self.absolute_time_s
    }

    /// Consume the envelope and recover the existing spatial ray.
    #[must_use]
    pub fn into_spatial(self) -> SpatialRay {
        self.spatial
    }
}

/// Structured shutter/ray-time refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MotionTimeError {
    /// Shot bounds were non-finite or decreasing.
    InvalidShotBounds,
    /// Frame time or exposure duration was invalid.
    InvalidExposure,
    /// Resolved exposure was non-finite or outside the shot.
    ExposureOutsideShot,
    /// Normalized shutter coordinate was not finite and inside `[0, 1]`.
    InvalidNormalizedTime,
}

impl fmt::Display for MotionTimeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for MotionTimeError {}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}
