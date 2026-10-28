//! Recorded earthquake histories through the existing seismic response owner.
//!
//! Records retain their own clock and duration. Ingestion converts explicit
//! acceleration units to SI without detrending, tapering, filtering, resampling,
//! or adding a free-decay tail. Recorded-suite IDA reuses the implicit bilinear
//! oscillator; PGA/PSA normalization compares records at one physical intensity.

use super::{FragilityPoint, bilinear_peak_ductility, sdof_coefficients, sdof_peak};

/// Units of the supplied acceleration samples before SI conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AccelerationUnit {
    /// Metres per second squared.
    MetresPerSecondSquared,
    /// Centimetres per second squared (Gal).
    CentimetresPerSecondSquared,
    /// Standard gravity, exactly 9.80665 m/s² per g.
    StandardGravity,
}

impl AccelerationUnit {
    fn to_si(self) -> f64 {
        match self {
            Self::MetresPerSecondSquared => 1.0,
            Self::CentimetresPerSecondSquared => 0.01,
            Self::StandardGravity => 9.80665,
        }
    }
}

/// Uniformly sampled acceleration in SI units.
///
/// N samples cover N-1 integration intervals. Response starts at rest at the
/// first sample and covers the supplied duration, including any explicitly
/// supplied zero-acceleration tail. Ingestion does not certify measurement
/// processing, usable frequency band, or temporal resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct GroundMotionRecord {
    dt_s: f64,
    acceleration_m_s2: Vec<f64>,
    description: Option<String>,
}

/// Recorded-motion admission refuses incomplete or ambiguous input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroundMotionRecordError {
    /// Invalid clock, duration, header, or unsupported quantity/units.
    Invalid(&'static str),
    /// Declared count exceeds the caller's ingestion budget.
    SampleLimit { requested: usize, maximum: usize },
    /// Encountered count disagrees with NPTS; stops at the first excess sample.
    SampleCount { expected: usize, actual: usize },
    /// Sample is malformed, nonfinite, or unrepresentable after SI conversion.
    Sample { index: usize },
}

impl core::fmt::Display for GroundMotionRecordError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-UQ-RECORD: {self:?}")
    }
}
impl std::error::Error for GroundMotionRecordError {}

impl GroundMotionRecord {
    /// Admit an acceleration history with explicit source units.
    ///
    /// # Errors
    /// Refuses fewer than two samples, invalid clock/duration, or any sample
    /// that is nonfinite or unrepresentable after conversion.
    #[allow(clippy::cast_precision_loss)]
    pub fn try_new(
        dt_s: f64,
        mut samples: Vec<f64>,
        units: AccelerationUnit,
    ) -> Result<Self, GroundMotionRecordError> {
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || samples.len() < 2
            || !((samples.len() - 1) as f64 * dt_s).is_finite()
        {
            return Err(GroundMotionRecordError::Invalid(
                "at least two samples and a positive finite clock/duration are required",
            ));
        }
        for (index, value) in samples.iter_mut().enumerate() {
            *value *= units.to_si();
            if !value.is_finite() {
                return Err(GroundMotionRecordError::Sample { index });
            }
        }
        Ok(Self {
            dt_s,
            acceleration_m_s2: samples,
            description: None,
        })
    }

    /// Read a four-header-line PEER acceleration record (AT2).
    ///
    /// Supports PEER and legacy Pacific Engineering title lines, acceleration
    /// HISTORY/SERIES in g, labelled NPTS/DT with seconds, and E/D exponents.
    /// The second header line is retained as the record description.
    ///
    /// Format: <https://opensees.berkeley.edu/OpenSees/manuals/usermanual/761.htm>.
    /// This reads supplied text only; it does not download or license records.
    ///
    /// # Errors
    /// Refuses unsupported quantity/units, missing or duplicate metadata,
    /// malformed samples, and disagreement with NPTS. Declared count is checked
    /// against maximum_samples BEFORE allocating sample storage.
    pub fn from_peer_at2(
        text: &str,
        maximum_samples: usize,
    ) -> Result<Self, GroundMotionRecordError> {
        use GroundMotionRecordError as Error;
        let mut lines = text.lines().filter(|line| !line.trim().is_empty());
        let title = lines
            .next()
            .ok_or(Error::Invalid("missing PEER title"))?
            .trim()
            .trim_start_matches('\u{feff}')
            .to_ascii_uppercase();
        if !title.starts_with("PEER ") && !title.starts_with("PACIFIC ENGINEERING") {
            return Err(Error::Invalid("unsupported PEER title"));
        }
        let description = lines
            .next()
            .ok_or(Error::Invalid("missing record description"))?;
        let quantity = lines
            .next()
            .ok_or(Error::Invalid("missing acceleration units"))?;
        let words: Vec<_> = quantity.split_whitespace().collect();
        if words.len() < 7
            || !words[0].eq_ignore_ascii_case("ACCELERATION")
            || !words[1].eq_ignore_ascii_case("TIME")
            || !(words[2].eq_ignore_ascii_case("SERIES")
                || words[2].eq_ignore_ascii_case("HISTORY"))
            || !words[3].eq_ignore_ascii_case("IN")
            || !words[4].eq_ignore_ascii_case("UNITS")
            || !words[5].eq_ignore_ascii_case("OF")
            || !words[6].trim_end_matches('.').eq_ignore_ascii_case("G")
        {
            return Err(Error::Invalid(
                "PEER input must declare acceleration in units of g",
            ));
        }
        let metadata = lines
            .next()
            .ok_or(Error::Invalid("missing NPTS/DT metadata"))?;
        let fields: Vec<_> = metadata
            .split(|c: char| c.is_ascii_whitespace() || c == ',' || c == '=')
            .filter(|field| !field.is_empty())
            .collect();
        let field = |key: &str| -> Result<(usize, &str), Error> {
            let mut matches = fields
                .iter()
                .enumerate()
                .filter(|(_, value)| value.eq_ignore_ascii_case(key));
            let (index, _) = matches.next().ok_or(Error::Invalid("missing NPTS or DT"))?;
            if matches.next().is_some() {
                return Err(Error::Invalid("duplicate NPTS or DT"));
            }
            let value = fields
                .get(index + 1)
                .copied()
                .ok_or(Error::Invalid("missing NPTS/DT value"))?;
            Ok((index, value))
        };
        let count = field("NPTS")?
            .1
            .parse::<usize>()
            .map_err(|_| Error::Invalid("invalid NPTS"))?;
        if count < 2 {
            return Err(Error::Invalid("NPTS must be at least two"));
        }
        if count > maximum_samples {
            return Err(Error::SampleLimit {
                requested: count,
                maximum: maximum_samples,
            });
        }
        let (dt_index, dt_value) = field("DT")?;
        if fields
            .get(dt_index + 2)
            .is_none_or(|unit| !unit.eq_ignore_ascii_case("SEC"))
        {
            return Err(Error::Invalid("DT must explicitly use SEC"));
        }
        let dt_s = number(dt_value).ok_or(Error::Invalid("invalid DT"))?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(Error::Invalid("DT must be positive and finite"));
        }
        let mut samples = Vec::with_capacity(count);
        for token in lines.flat_map(str::split_whitespace) {
            if samples.len() == count {
                return Err(Error::SampleCount {
                    expected: count,
                    actual: count + 1,
                });
            }
            let value = number(token).ok_or(Error::Sample {
                index: samples.len(),
            })?;
            samples.push(value);
        }
        if samples.len() != count {
            return Err(Error::SampleCount {
                expected: count,
                actual: samples.len(),
            });
        }
        let mut record = Self::try_new(dt_s, samples, AccelerationUnit::StandardGravity)?;
        record.description = Some(description.trim().to_owned());
        Ok(record)
    }

    /// Original sample interval [s].
    #[must_use]
    pub fn dt_s(&self) -> f64 {
        self.dt_s
    }

    /// Samples in m/s², with their original order and signs.
    #[must_use]
    pub fn acceleration_m_s2(&self) -> &[f64] {
        &self.acceleration_m_s2
    }

    /// AT2 event/station/component description, when supplied.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Duration between the first and final sample [s].
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn duration_s(&self) -> f64 {
        (self.acceleration_m_s2.len() - 1) as f64 * self.dt_s
    }

    /// Maximum absolute sampled acceleration [m/s²].
    #[must_use]
    pub fn peak_ground_acceleration(&self) -> f64 {
        self.acceleration_m_s2
            .iter()
            .fold(0.0_f64, |peak, a| peak.max(a.abs()))
    }

    /// Pseudo-spectral acceleration omega² * max|u_relative| [m/s²].
    ///
    /// Uses the existing linear Newmark response on this record's clock and
    /// duration. This differs from peak absolute oscillator acceleration.
    /// Resolve the response in time; ingestion cannot recover unmeasured peaks.
    ///
    /// # Panics
    /// Invalid period/damping, or unrepresentable response.
    #[must_use]
    pub fn pseudo_spectral_acceleration(&self, period_s: f64, damping_ratio: f64) -> f64 {
        assert!(
            period_s.is_finite() && period_s > 0.0,
            "spectral period must be positive and finite"
        );
        let omega = core::f64::consts::TAU / period_s;
        let value =
            omega * omega * sdof_peak(&self.acceleration_m_s2, self.dt_s, omega, damping_ratio);
        assert!(
            value.is_finite(),
            "spectral acceleration must be representable"
        );
        value
    }
}

fn number(text: &str) -> Option<f64> {
    text.replace(['D', 'd'], "E").parse::<f64>().ok()
}

/// Intensity shared by every record at an IDA level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IdaIntensity {
    /// Dimensionless multiplier of each original history.
    ScaleFactor,
    /// Target peak ground acceleration [m/s²].
    PeakGroundAcceleration,
    /// Target linear pseudo-spectral acceleration [m/s²] at this period/loss.
    PseudoSpectralAcceleration { period_s: f64, damping_ratio: f64 },
}

/// Physical oscillator and exceedance criterion for recorded-motion IDA.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BilinearIdaModel {
    /// Initial natural angular frequency [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Viscous damping ratio relative to the initial stiffness.
    pub damping_ratio: f64,
    /// Yield displacement [m].
    pub yield_displacement_m: f64,
    /// Post-yield/initial stiffness ratio in [0,1].
    pub post_yield_stiffness_ratio: f64,
    /// Exceedance is strict: peak ductility > this nonnegative limit.
    pub ductility_limit: f64,
}

/// Demands and empirical exceedances in caller-supplied order.
#[derive(Debug, Clone, PartialEq)]
pub struct RecordedIda {
    /// Meaning and units of the fragility ordinates.
    pub intensity_measure: IdaIntensity,
    /// Original intensity per record, or one for ScaleFactor.
    pub reference_intensities: Vec<f64>,
    /// Multipliers: outer index = intensity, inner index = record.
    pub record_scales: Vec<Vec<f64>>,
    /// Nonlinear peak demands in the same intensity-by-record order.
    pub peak_ductilities: Vec<Vec<f64>>,
    /// Empirical exceedance fractions; not population probabilities or CSs.
    pub fragility: Vec<FragilityPoint>,
}

/// Incremental dynamic analysis of supplied recorded motions, each with its
/// original clock, through the existing implicit bilinear response owner.
///
/// PGA/PSA references are computed once per record. Every level borrows the
/// original record: no scaled history is allocated, and no samples are aligned
/// across unrelated motions. A physical target I gives scale I/I_record.
/// ScaleFactor instead uses I itself. Signs and waveform are preserved.
///
/// Exceedance fractions do not establish a sampling law, confidence coverage,
/// or structural-model validity. This is the existing SDOF consumer, not a
/// fiber-frame fragility analysis. No monotonicity is imposed on demands.
///
/// # Panics
/// Empty suite, invalid oscillator/intensities, zero record intensity when
/// normalization is requested, or unrepresentable scaling/response.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn ida_fragility_from_records(
    records: &[GroundMotionRecord],
    intensity_measure: IdaIntensity,
    im_ladder: &[f64],
    model: BilinearIdaModel,
) -> RecordedIda {
    assert!(
        !records.is_empty(),
        "recorded IDA needs a nonempty motion suite"
    );
    let BilinearIdaModel {
        angular_frequency_rad_s: wn,
        damping_ratio: zeta,
        yield_displacement_m: uy,
        post_yield_stiffness_ratio: alpha,
        ductility_limit,
    } = model;
    assert!(
        uy.is_finite() && uy > 0.0,
        "yield displacement must be positive and finite"
    );
    assert!(
        (0.0..=1.0).contains(&alpha),
        "post-yield stiffness ratio must be in [0,1]"
    );
    assert!(
        ductility_limit.is_finite() && ductility_limit >= 0.0,
        "ductility limit must be finite and nonnegative"
    );
    assert!(
        im_ladder.iter().all(|im| im.is_finite() && *im >= 0.0),
        "IDA intensities must be finite and nonnegative"
    );
    let reference_intensities: Vec<f64> = records
        .iter()
        .map(|record| {
            sdof_coefficients(record.dt_s, wn, zeta);
            let reference = match intensity_measure {
                IdaIntensity::ScaleFactor => 1.0,
                IdaIntensity::PeakGroundAcceleration => record.peak_ground_acceleration(),
                IdaIntensity::PseudoSpectralAcceleration {
                    period_s,
                    damping_ratio,
                } => record.pseudo_spectral_acceleration(period_s, damping_ratio),
            };
            assert!(
                reference.is_finite() && reference > 0.0,
                "record intensity must be positive and finite for normalization"
            );
            reference
        })
        .collect();
    let mut report = RecordedIda {
        intensity_measure,
        reference_intensities,
        record_scales: Vec::with_capacity(im_ladder.len()),
        peak_ductilities: Vec::with_capacity(im_ladder.len()),
        fragility: Vec::with_capacity(im_ladder.len()),
    };
    for &im in im_ladder {
        let mut scales = Vec::with_capacity(records.len());
        let mut demands = Vec::with_capacity(records.len());
        let mut hits = 0;
        for (record, &reference) in records.iter().zip(&report.reference_intensities) {
            let scale = im / reference;
            assert!(scale.is_finite(), "record scaling must be representable");
            let demand = bilinear_peak_ductility(
                &record.acceleration_m_s2,
                record.dt_s,
                wn,
                zeta,
                uy,
                alpha,
                scale,
            );
            assert!(
                demand.is_finite(),
                "recorded IDA ductility must be representable"
            );
            hits += usize::from(demand > ductility_limit);
            scales.push(scale);
            demands.push(demand);
        }
        report.record_scales.push(scales);
        report.peak_ductilities.push(demands);
        report.fragility.push(FragilityPoint {
            im,
            p: hits as f64 / records.len() as f64,
            n: records.len(),
        });
    }
    report
}
