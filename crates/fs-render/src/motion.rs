//! Deterministic shutter-time semantics and rays carrying explicit time.

use core::fmt;

/// Stable counter-domain separation for shutter samples.
pub const SHUTTER_TIME_SAMPLE_DOMAIN_V1: u64 = 0x6d6f_7469_6f6e_0001;

const SHUTTER_PIXEL_IDENTITY_DOMAIN_V1: u64 = 0x7069_7865_6c00_0001;
const SHUTTER_SAMPLE_IDENTITY_DOMAIN_V1: u64 = 0x7361_6d70_6c65_0001;
const SHUTTER_STRATUM_PERMUTATION_DOMAIN_V1: u64 = 0x7374_7261_7461_0001;
const SHUTTER_STREAM_IDENTITY_DOMAIN_V1: u64 = 0x7374_7265_616d_0001;
const GOLDEN_RATIO_CONJUGATE_FIXED_U64: u64 = 0x9e37_79b9_7f4a_7c15;

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
    /// Jittered strata, keyed by absolute logical sample identity. Every
    /// consecutive window of `strata` sample identities visits each stratum
    /// once, while retaining replay identity across progressive partitions.
    StratifiedCounterV1 {
        /// Number of equal-width temporal strata. Must be nonzero.
        strata: u32,
    },
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

    /// Inclusive shot start [s].
    #[must_use]
    pub const fn start_s(self) -> f64 {
        self.start_s
    }

    /// Inclusive shot end [s].
    #[must_use]
    pub const fn end_s(self) -> f64 {
        self.end_s
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
    /// Re-admit already-resolved canonical bounds from a validated binary
    /// artifact such as a renderer checkpoint.
    ///
    /// This crate-private seam avoids reconstructing a resolved interval from
    /// a nominal frame time, which can change endpoint bits through a second
    /// round of floating-point arithmetic. It enforces the same stored-field
    /// invariants as [`Self::resolve`] and canonicalizes the one physical zero.
    pub(crate) fn try_from_canonical_parts(
        open_s: f64,
        close_s: f64,
        convention: ShutterConvention,
        distribution: ShutterDistribution,
    ) -> Result<Self, MotionTimeError> {
        if !open_s.is_finite() || !close_s.is_finite() || open_s > close_s {
            return Err(MotionTimeError::InvalidExposure);
        }
        if matches!(
            distribution,
            ShutterDistribution::StratifiedCounterV1 { strata: 0 }
        ) {
            return Err(MotionTimeError::InvalidDistribution);
        }
        Ok(Self {
            open_s: canonical_zero(open_s),
            close_s: canonical_zero(close_s),
            convention,
            distribution,
        })
    }

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
        if matches!(
            distribution,
            ShutterDistribution::StratifiedCounterV1 { strata: 0 }
        ) {
            return Err(MotionTimeError::InvalidDistribution);
        }
        let (mut open_s, mut close_s) = match convention {
            ShutterConvention::Centered => {
                let half = 0.5 * exposure_duration_s;
                (frame_time_s - half, frame_time_s + half)
            }
            ShutterConvention::FrontLoaded => (frame_time_s, frame_time_s + exposure_duration_s),
            ShutterConvention::BackLoaded => (frame_time_s - exposure_duration_s, frame_time_s),
        };
        // There is one physical zero instant. Canonical storage makes exact
        // progressive-shutter provenance independent of an input zero's sign.
        open_s = canonical_zero(open_s);
        close_s = canonical_zero(close_s);
        let open_bits = open_s.to_bits();
        let close_bits = close_s.to_bits();
        if exposure_duration_s > 0.0
            && (open_bits == close_bits || (open_bits << 1 == 0 && close_bits << 1 == 0))
        {
            return Err(MotionTimeError::CollapsedExposure);
        }
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
        match normalized_time.value() {
            0.0 => self.open_s,
            1.0 => self.close_s,
            coordinate => self.duration_s().mul_add(coordinate, self.open_s),
        }
    }

    /// Deterministically sample a shutter coordinate from stable logical IDs
    /// in the compatibility stream with identity zero.
    #[must_use]
    pub fn sample(self, pixel_identity: u64, sample_identity: u64) -> NormalizedShutterTime {
        self.sample_for_stream(0, pixel_identity, sample_identity)
    }

    /// Deterministically sample a shutter coordinate from an explicit render
    /// stream plus stable pixel and absolute logical-sample identities.
    #[must_use]
    pub fn sample_for_stream(
        self,
        stream_identity: u64,
        pixel_identity: u64,
        sample_identity: u64,
    ) -> NormalizedShutterTime {
        match self.distribution {
            ShutterDistribution::UniformCounterV1 => {
                // Stream zero is the public V1 compatibility stream used by
                // `sample` and `TimedRay::from_sample`. Preserve its original
                // bit sequence; nonzero streams use the independently keyed
                // extension introduced with `sample_for_stream`.
                let bits = if stream_identity == 0 {
                    legacy_uniform_counter_bits(pixel_identity, sample_identity)
                } else {
                    shutter_counter_bits(stream_identity, pixel_identity, sample_identity)
                };
                NormalizedShutterTime::from_counter_bits(bits)
            }
            ShutterDistribution::StratifiedCounterV1 { strata } => {
                debug_assert_ne!(strata, 0, "distribution is validated by resolve");
                let strata_u64 = u64::from(strata);
                let stratum =
                    permuted_stratum(stream_identity, pixel_identity, sample_identity, strata_u64);
                let jitter = NormalizedShutterTime::from_counter_bits(shutter_counter_bits(
                    stream_identity,
                    pixel_identity,
                    sample_identity,
                ));
                stratified_coordinate(stratum, strata, jitter)
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
    /// Complete admitted shutter definition that produced this ray time.
    shutter: ShutterInterval,
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
        Self::from_stream_sample(spatial, shutter, 0, pixel_identity, sample_identity)
    }

    /// Bind an existing spatial ray to one deterministic shutter sample in an
    /// explicit render stream.
    #[must_use]
    pub fn from_stream_sample(
        spatial: SpatialRay,
        shutter: ShutterInterval,
        stream_identity: u64,
        pixel_identity: u64,
        sample_identity: u64,
    ) -> Self {
        let normalized_time =
            shutter.sample_for_stream(stream_identity, pixel_identity, sample_identity);
        Self {
            spatial,
            shutter,
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
            shutter,
            normalized_time,
            absolute_time_s: shutter.time_at(normalized_time),
        }
    }

    /// Borrow the backend-specific spatial ray.
    #[must_use]
    pub const fn spatial(&self) -> &SpatialRay {
        &self.spatial
    }

    /// Complete admitted shutter definition used to derive the ray time.
    #[must_use]
    pub const fn shutter(&self) -> ShutterInterval {
        self.shutter
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
    /// A positive requested duration was below the absolute-time resolution
    /// and collapsed to one representable endpoint.
    CollapsedExposure,
    /// Resolved exposure was non-finite or outside the shot.
    ExposureOutsideShot,
    /// Normalized shutter coordinate was not finite and inside `[0, 1]`.
    InvalidNormalizedTime,
    /// A shutter distribution had invalid parameters.
    InvalidDistribution,
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

fn legacy_uniform_counter_bits(pixel_identity: u64, sample_identity: u64) -> u64 {
    let counter = pixel_identity
        .wrapping_mul(GOLDEN_RATIO_CONJUGATE_FIXED_U64)
        .wrapping_add(sample_identity)
        ^ SHUTTER_TIME_SAMPLE_DOMAIN_V1;
    splitmix64(counter)
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn stratified_coordinate(
    stratum: u64,
    strata: u32,
    jitter: NormalizedShutterTime,
) -> NormalizedShutterTime {
    let denominator = f64::from(strata);
    let mut coordinate = (stratum as f64 + jitter.value()) / denominator;
    // Exact arithmetic lies in this stratum's half-open interval, but either
    // boundary can round into its neighbor (not only the global 1.0 endpoint).
    // Nudge only those exceptional results until the public floor-based stratum
    // recovery agrees with the permutation. For u32 strata the initial result
    // is at most a handful of ulps from the intended bucket.
    for _ in 0..8 {
        let recovered = (coordinate * denominator).floor() as u64;
        if recovered == stratum {
            return NormalizedShutterTime(coordinate);
        }
        coordinate = if recovered < stratum {
            next_up_unit(coordinate)
        } else {
            next_down_unit(coordinate)
        };
    }
    // The midpoint is separated from either boundary by half a stratum, which
    // dwarfs binary64 rounding for the admitted u32 stratum count.
    coordinate = (stratum as f64 + 0.5) / denominator;
    debug_assert_eq!((coordinate * denominator).floor() as u64, stratum);
    NormalizedShutterTime(coordinate)
}

fn next_up_unit(value: f64) -> f64 {
    if value == 0.0 {
        f64::from_bits(1)
    } else {
        f64::from_bits(value.to_bits() + 1)
    }
}

fn next_down_unit(value: f64) -> f64 {
    if value == 0.0 {
        0.0
    } else {
        f64::from_bits(value.to_bits() - 1)
    }
}

fn shutter_counter_bits(stream_identity: u64, pixel_identity: u64, sample_identity: u64) -> u64 {
    let stream = splitmix64(stream_identity ^ SHUTTER_STREAM_IDENTITY_DOMAIN_V1);
    let pixel = splitmix64(pixel_identity ^ SHUTTER_PIXEL_IDENTITY_DOMAIN_V1);
    let sample = splitmix64(sample_identity ^ SHUTTER_SAMPLE_IDENTITY_DOMAIN_V1);
    splitmix64(
        SHUTTER_TIME_SAMPLE_DOMAIN_V1
            ^ stream.rotate_left(7)
            ^ pixel.rotate_left(17)
            ^ sample.rotate_right(11),
    )
}

fn permuted_stratum(
    stream_identity: u64,
    pixel_identity: u64,
    sample_identity: u64,
    strata: u64,
) -> u64 {
    if strata == 1 {
        return 0;
    }
    let base_stratum = sample_identity % strata;
    let permutation_bits = shutter_counter_bits(
        stream_identity,
        pixel_identity ^ SHUTTER_STRATUM_PERMUTATION_DOMAIN_V1,
        0,
    );
    // A step near the golden-ratio conjugate gives useful prefix coverage even
    // when the render consumes far fewer samples than declared strata. Moving
    // to the next coprime retains a true permutation for arbitrary stratum
    // counts; the key controls the cyclic offset and within-stratum jitter.
    let scaled = (u128::from(strata) * u128::from(GOLDEN_RATIO_CONJUGATE_FIXED_U64)) >> 64;
    let mut multiplier = u64::try_from(scaled).unwrap_or(1).max(1);
    while greatest_common_divisor(multiplier, strata) != 1 {
        multiplier += 1;
        if multiplier == strata {
            multiplier = 1;
        }
    }
    let offset = permutation_bits.rotate_left(29) % strata;
    (multiplier * base_stratum + offset) % strata
}

fn greatest_common_divisor(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        (left, right) = (right, left % right);
    }
    left
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_rounding_stays_in_the_permuted_stratum() {
        let maximum_jitter = NormalizedShutterTime::from_counter_bits(u64::MAX);
        for (stratum, strata, jitter) in [
            (1, 8, maximum_jitter),
            (7, 8, maximum_jitter),
            (15, 22, NormalizedShutterTime(0.0)),
        ] {
            let coordinate = stratified_coordinate(stratum, strata, jitter).value();
            assert!((0.0..1.0).contains(&coordinate));
            assert_eq!((coordinate * f64::from(strata)).floor() as u64, stratum);
        }
    }
}
