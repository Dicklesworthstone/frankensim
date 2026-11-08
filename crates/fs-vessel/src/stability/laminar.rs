//! Dimensional, locally parallel free-surface film screening.
//!
//! This is the physical alternative to the historical two-wall channel proxy
//! in the parent module. Each section declares its actual width and effective
//! acceleration components; no radius-profile slope is invented as gravity.
//! A section uses the native Nusselt film owner, including its long-wave
//! admission budgets. Surface-wave growth is not a turbulence certificate.
//! A laterally wide film (thickness small compared with width), gradual geometry
//! variation, constant Newtonian properties, no slip, negligible
//! gas shear and a supported free surface are assumptions of this local screen.

use fs_lbm::film::{FilmError, FilmMaterial, LongWaveLimits, NusseltFilm};

/// Finite-wavelength free-surface spectra using the native Chebyshev owner.
pub mod spectral;
pub use spectral::{LaminarFilmSpectralStudy, LaminarSpectralError, LaminarSpectralReport};

mod transit_band;
pub use transit_band::{LaminarFilmTransitBandStudy, LaminarTransitBandReport};

/// One locally parallel section of the actual liquid path, in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarFilmSection {
    /// Wetted span transverse to the flow; lateral-wall effects are omitted.
    pub width_m: f64,
    /// Effective acceleration in the downstream direction.
    pub tangential_acceleration_m_s2: f64,
    /// Effective acceleration toward the supporting wall.
    pub normal_acceleration_m_s2: f64,
}

/// Prescribed total liquid flow and constant properties at one operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarPourConditions {
    /// Total volume flow, conserved through every section, in m³/s.
    pub volume_flow_m3_s: f64,
    /// Actual liquid properties; no water/air constants are substituted.
    pub material: FilmMaterial,
}

/// A borrowed path and a common set of physical streamwise wavenumbers.
///
/// The caller chooses the sampled wavelengths and approximation budgets. A
/// negative result only means these admitted long-wave modes decay locally.
/// This does not maximize over all wavelengths, solve finite-length inlet/outlet
/// modes, or replace the transient free-surface pour.
#[derive(Debug, Clone, Copy)]
pub struct LaminarFilmStudy<'a> {
    /// Sections in path order.
    pub sections: &'a [LaminarFilmSection],
    /// Sampled angular wavenumbers, in rad/m.
    pub wavenumbers_rad_m: &'a [f64],
    /// Explicit long-wave/asymptotic admission budgets.
    pub limits: LongWaveLimits,
}

/// One locally parallel part of a finite lip path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarTransitSection {
    /// The same physical base-flow geometry as the local stability study.
    pub film: LaminarFilmSection,
    /// Downstream distance over which this local base state applies, in metres.
    pub length_m: f64,
}

/// Leading-order transport of one inlet disturbance along a slowly varying lip.
///
/// The steady path preserves the disturbance's angular frequency, not its
/// wavenumber: locally k = omega/c with the Nusselt kinematic-wave speed
/// c = d q/d h = 3 q/h. The weak-growth spatial exponent is s/c, so a section
/// contributes s L/c to the log amplitude. `max_growth_frequency_ratio` bounds
/// |s/omega|, in addition to the native long-wave restrictions. These are
/// applicability budgets, not a spatial Orr--Sommerfeld error estimate.
///
/// The transported amplitude is relative film thickness (equivalently relative
/// volume-flow perturbation). Since width*c*h = 3 Q, a gradual width change
/// introduces no extra relative-amplitude factor at this order. Abrupt section
/// jumps, reflections, turning points, nonlinear waves, and inlet/exit boundary
/// layers are outside this local transport model.
#[derive(Debug, Clone, Copy)]
pub struct LaminarFilmTransitStudy<'a> {
    /// Positive-length sections in path order, approximating gradual variation.
    pub sections: &'a [LaminarTransitSection],
    /// Positive inlet angular frequency, in rad/s, held fixed along the path.
    pub angular_frequency_rad_s: f64,
    /// Native local long-wave budgets, enforced at every section.
    pub limits: LongWaveLimits,
    /// Maximum |local temporal growth / angular frequency|, between zero and one.
    pub max_growth_frequency_ratio: f64,
}

/// Finite-path linear amplification; no exponential is formed or clipped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarTransitReport {
    /// Natural log of outlet relative amplitude divided by inlet amplitude.
    pub log_amplification: f64,
    /// Largest accumulated log amplification, including the inlet value zero.
    pub peak_log_amplification: f64,
    /// Number of completed sections at the first peak; zero denotes the inlet.
    pub peak_after_sections: usize,
    /// Transit time of the surface disturbance, sum L/(3*mean speed), in seconds.
    pub wave_transit_time_s: f64,
    /// Residence time of the liquid at its mean velocity, in seconds.
    pub liquid_residence_time_s: f64,
    /// Smallest and largest local angular wavenumbers, in rad/m.
    pub wavenumber_range_rad_m: [f64; 2],
    /// Minimum and maximum mean-speed film Reynolds numbers.
    pub reynolds_range: [f64; 2],
    /// Largest local kh.
    pub max_kh: f64,
    /// Largest local Re*kh.
    pub max_re_kh: f64,
    /// Largest admitted |s/omega| along the path.
    pub max_growth_frequency_ratio: f64,
}

/// Growth objective shared by sampled and continuous-band physical studies.
/// Risk and racing consume this owner without reimplementing film physics.
pub trait LaminarGrowthObjective {
    /// Worst admitted local temporal growth, in 1/s.
    fn growth_per_s(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<f64, LaminarFilmError>;
}

/// Physical study over every wavenumber in an explicit closed SI interval.
/// The native quartic dispersion needs at most three evaluations per section;
/// no sampled grid can miss an interior unstable long-wave mode.
#[derive(Debug, Clone, Copy)]
pub struct LaminarFilmBandStudy<'a> {
    /// Locally planar sections in path order.
    pub sections: &'a [LaminarFilmSection],
    /// Inclusive nonnegative angular-wavenumber interval, in rad/m.
    pub wavenumbers_rad_m: [f64; 2],
    /// Both limits must admit the entire interval at every section.
    pub limits: LongWaveLimits,
}

/// Maximum of the local Benney model over the declared continuous band.
/// This remains an Estimated model result, not a numerical enclosure or a
/// finite-length/global-stability certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarBandGrowthReport {
    /// Largest temporal growth rate, in 1/s.
    pub max_growth_per_s: f64,
    /// First section attaining the maximum.
    pub section: usize,
    /// Physical wavenumber attaining that maximum, in rad/m.
    pub wavenumber_rad_m: f64,
    /// Minimum and maximum mean-speed film Reynolds numbers.
    pub reynolds_range: [f64; 2],
    /// Maximum kh over every section and the complete interval.
    pub max_kh: f64,
    /// Maximum Re*kh over every section and the complete interval.
    pub max_re_kh: f64,
}

/// Worst sampled local temporal growth and the actual nondimensional range.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarGrowthReport {
    /// Largest sampled temporal growth rate, in 1/s.
    pub max_growth_per_s: f64,
    /// Path section attaining the maximum; first section wins exact ties.
    pub section: usize,
    /// Wavenumber attaining the maximum; first mode wins exact ties.
    pub mode: usize,
    /// Minimum and maximum mean-speed film Reynolds numbers, q/nu.
    pub reynolds_range: [f64; 2],
    /// Largest sampled k h.
    pub max_kh: f64,
    /// Largest sampled Re k h.
    pub max_re_kh: f64,
}

/// A physical film study that cannot report a growth value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LaminarFilmError {
    /// The native free-surface eigensolve was refused. The spectral study's
    /// direct evaluate method retains the full residual/resolution diagnosis.
    Spectral { section: usize, mode: usize },
    /// At least one section is needed.
    EmptyPath,
    /// At least one physical wavenumber is needed.
    EmptyModes,
    /// The width multiplier must be finite and strictly positive.
    InvalidWidthScale,
    /// Total volume flow must be finite and strictly positive.
    InvalidVolumeFlow,
    /// The transit frequency/band or the spatial-growth budget is invalid.
    InvalidTransitInput,
    /// A section length is not finite and positive.
    InvalidLength {
        /// Path section.
        section: usize,
    },
    /// Temporal-to-spatial conversion exceeds its declared weak-growth budget.
    OutsideSpatialGrowthRegime {
        /// Path section.
        section: usize,
    },
    /// A declared/scaled section width is not positive and representable.
    InvalidWidth {
        /// Path section.
        section: usize,
    },
    /// The native base flow refused its inputs.
    Film {
        /// Path section.
        section: usize,
        /// Native refusal, without substitution or clipping.
        source: FilmError,
    },
    /// A requested wave violates the native physical/asymptotic domain.
    Mode {
        /// Path section.
        section: usize,
        /// Wavenumber index.
        mode: usize,
        /// Native refusal, without substitution or clipping.
        source: FilmError,
    },
    /// The requested continuous interval was refused by the native owner.
    Band {
        /// Path section.
        section: usize,
        /// Original physical or long-wave admission refusal.
        source: FilmError,
    },
}

impl core::fmt::Display for LaminarFilmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Spectral { section, mode } => {
                write!(
                    f,
                    "free-surface spectrum refused at film section {section}, wavenumber {mode}"
                )
            }
            Self::EmptyPath => write!(f, "a laminar film study needs a path section"),
            Self::EmptyModes => write!(f, "a laminar film study needs a physical wavenumber"),
            Self::InvalidWidthScale => write!(f, "film width scale must be finite and positive"),
            Self::InvalidVolumeFlow => write!(f, "film volume flow must be finite and positive"),
            Self::InvalidTransitInput => {
                write!(f, "invalid film transit frequency or growth budget")
            }
            Self::InvalidLength { section } => {
                write!(f, "film section {section} needs a finite positive length")
            }
            Self::OutsideSpatialGrowthRegime { section } => {
                write!(
                    f,
                    "film section {section} exceeds the weak spatial-growth budget"
                )
            }
            Self::InvalidWidth { section } => {
                write!(f, "film section {section} has an invalid scaled width")
            }
            Self::Film { section, source } => write!(f, "film section {section}: {source}"),
            Self::Mode {
                section,
                mode,
                source,
            } => {
                write!(f, "film section {section}, mode {mode}: {source}")
            }
            Self::Band { section, source } => {
                write!(f, "film section {section}, wavenumber band: {source}")
            }
        }
    }
}

impl std::error::Error for LaminarFilmError {}

impl LaminarFilmTransitStudy<'_> {
    /// Integrate signed amplification along the finite path in O(sections)
    /// arithmetic and constant storage, retaining all native physical refusals.
    pub fn evaluate(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<LaminarTransitReport, LaminarFilmError> {
        if self.sections.is_empty() {
            return Err(LaminarFilmError::EmptyPath);
        }
        if !width_scale.is_finite() || width_scale <= 0.0 {
            return Err(LaminarFilmError::InvalidWidthScale);
        }
        if !conditions.volume_flow_m3_s.is_finite() || conditions.volume_flow_m3_s <= 0.0 {
            return Err(LaminarFilmError::InvalidVolumeFlow);
        }
        let omega = self.angular_frequency_rad_s;
        let budget = self.max_growth_frequency_ratio;
        if !omega.is_finite()
            || omega <= 0.0
            || !budget.is_finite()
            || budget <= 0.0
            || budget >= 1.0
        {
            return Err(LaminarFilmError::InvalidTransitInput);
        }
        let mut report = LaminarTransitReport {
            log_amplification: 0.0,
            peak_log_amplification: 0.0,
            peak_after_sections: 0,
            wave_transit_time_s: 0.0,
            liquid_residence_time_s: 0.0,
            wavenumber_range_rad_m: [f64::INFINITY, 0.0],
            reynolds_range: [f64::INFINITY, 0.0],
            max_kh: 0.0,
            max_re_kh: 0.0,
            max_growth_frequency_ratio: 0.0,
        };
        for (section, segment) in self.sections.iter().enumerate() {
            let length = segment.length_m;
            if !length.is_finite() || length <= 0.0 {
                return Err(LaminarFilmError::InvalidLength { section });
            }
            let film = section_film(segment.film, width_scale, conditions, section)?;
            let failure = || LaminarFilmError::Film {
                section,
                source: FilmError::NonRepresentable,
            };
            let speed = 3.0 * film.mean_velocity_m_s();
            let k = omega / speed;
            let liquid_time = length / film.mean_velocity_m_s();
            let wave_time = liquid_time / 3.0;
            if [speed, k, liquid_time, wave_time]
                .iter()
                .any(|&x| !x.is_finite() || x <= 0.0)
            {
                return Err(failure());
            }
            let wave =
                film.long_wave_mode(k, self.limits)
                    .map_err(|source| LaminarFilmError::Mode {
                        section,
                        mode: 0,
                        source,
                    })?;
            let ratio = wave.growth_per_s.abs() / omega;
            if ratio > budget {
                return Err(LaminarFilmError::OutsideSpatialGrowthRegime { section });
            }
            let increment = wave.growth_per_s * wave_time;
            if !increment.is_finite() || (increment == 0.0 && wave.growth_per_s != 0.0) {
                return Err(failure());
            }
            report.log_amplification += increment;
            report.wave_transit_time_s += wave_time;
            report.liquid_residence_time_s += liquid_time;
            if [
                report.log_amplification,
                report.wave_transit_time_s,
                report.liquid_residence_time_s,
            ]
            .iter()
            .any(|x| !x.is_finite())
            {
                return Err(failure());
            }
            if report.log_amplification > report.peak_log_amplification {
                report.peak_log_amplification = report.log_amplification;
                report.peak_after_sections = section + 1;
            }
            report.wavenumber_range_rad_m[0] = report.wavenumber_range_rad_m[0].min(k);
            report.wavenumber_range_rad_m[1] = report.wavenumber_range_rad_m[1].max(k);
            report.reynolds_range[0] = report.reynolds_range[0].min(film.reynolds());
            report.reynolds_range[1] = report.reynolds_range[1].max(film.reynolds());
            report.max_kh = report.max_kh.max(wave.kh);
            report.max_re_kh = report.max_re_kh.max(wave.re_kh);
            report.max_growth_frequency_ratio = report.max_growth_frequency_ratio.max(ratio);
        }
        Ok(report)
    }
}

impl LaminarFilmStudy<'_> {
    /// Evaluate a width candidate at one physical operating point.
    ///
    /// The multiplier changes every declared section width, while total volume
    /// flow stays fixed: q = Q / (width * scale), Re = q / nu. Every mode at
    /// every section must be admitted before a report is returned. Work is
    /// O(sections * modes), with no heap allocation or eigensolver per trial.
    pub fn evaluate(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<LaminarGrowthReport, LaminarFilmError> {
        if self.sections.is_empty() {
            return Err(LaminarFilmError::EmptyPath);
        }
        if self.wavenumbers_rad_m.is_empty() {
            return Err(LaminarFilmError::EmptyModes);
        }
        if !width_scale.is_finite() || width_scale <= 0.0 {
            return Err(LaminarFilmError::InvalidWidthScale);
        }
        if !conditions.volume_flow_m3_s.is_finite() || conditions.volume_flow_m3_s <= 0.0 {
            return Err(LaminarFilmError::InvalidVolumeFlow);
        }
        let mut report = LaminarGrowthReport {
            max_growth_per_s: f64::NEG_INFINITY,
            section: 0,
            mode: 0,
            reynolds_range: [f64::INFINITY, 0.0],
            max_kh: 0.0,
            max_re_kh: 0.0,
        };
        for (section, geometry) in self.sections.iter().enumerate() {
            let film = section_film(*geometry, width_scale, conditions, section)?;
            let re = film.reynolds();
            report.reynolds_range[0] = report.reynolds_range[0].min(re);
            report.reynolds_range[1] = report.reynolds_range[1].max(re);
            for (mode, &wavenumber) in self.wavenumbers_rad_m.iter().enumerate() {
                let wave = film
                    .long_wave_mode(wavenumber, self.limits)
                    .map_err(|source| LaminarFilmError::Mode {
                        section,
                        mode,
                        source,
                    })?;
                if wave.growth_per_s > report.max_growth_per_s {
                    report.max_growth_per_s = wave.growth_per_s;
                    report.section = section;
                    report.mode = mode;
                }
                report.max_kh = report.max_kh.max(wave.kh);
                report.max_re_kh = report.max_re_kh.max(wave.re_kh);
            }
        }
        Ok(report)
    }
}

impl LaminarGrowthObjective for LaminarFilmStudy<'_> {
    fn growth_per_s(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<f64, LaminarFilmError> {
        self.evaluate(width_scale, conditions)
            .map(|report| report.max_growth_per_s)
    }
}

impl LaminarFilmBandStudy<'_> {
    /// Evaluate all sections and the entire band in O(section count) work,
    /// without allocating or performing a wavenumber search.
    pub fn evaluate(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<LaminarBandGrowthReport, LaminarFilmError> {
        if self.sections.is_empty() {
            return Err(LaminarFilmError::EmptyPath);
        }
        if !width_scale.is_finite() || width_scale <= 0.0 {
            return Err(LaminarFilmError::InvalidWidthScale);
        }
        if !conditions.volume_flow_m3_s.is_finite() || conditions.volume_flow_m3_s <= 0.0 {
            return Err(LaminarFilmError::InvalidVolumeFlow);
        }
        let mut report = LaminarBandGrowthReport {
            max_growth_per_s: f64::NEG_INFINITY,
            section: 0,
            wavenumber_rad_m: 0.0,
            reynolds_range: [f64::INFINITY, 0.0],
            max_kh: 0.0,
            max_re_kh: 0.0,
        };
        for (section, &geometry) in self.sections.iter().enumerate() {
            let film = section_film(geometry, width_scale, conditions, section)?;
            let maximum = film
                .long_wave_band(self.wavenumbers_rad_m, self.limits)
                .map_err(|source| LaminarFilmError::Band { section, source })?;
            if maximum.mode.growth_per_s > report.max_growth_per_s {
                report.max_growth_per_s = maximum.mode.growth_per_s;
                report.section = section;
                report.wavenumber_rad_m = maximum.wavenumber_rad_m;
            }
            report.reynolds_range[0] = report.reynolds_range[0].min(film.reynolds());
            report.reynolds_range[1] = report.reynolds_range[1].max(film.reynolds());
            report.max_kh = report.max_kh.max(maximum.max_kh);
            report.max_re_kh = report.max_re_kh.max(maximum.max_re_kh);
        }
        Ok(report)
    }
}

impl LaminarGrowthObjective for LaminarFilmBandStudy<'_> {
    fn growth_per_s(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<f64, LaminarFilmError> {
        self.evaluate(width_scale, conditions)
            .map(|report| report.max_growth_per_s)
    }
}

fn section_film(
    geometry: LaminarFilmSection,
    width_scale: f64,
    conditions: LaminarPourConditions,
    section: usize,
) -> Result<NusseltFilm, LaminarFilmError> {
    let width = geometry.width_m * width_scale;
    if !geometry.width_m.is_finite()
        || geometry.width_m <= 0.0
        || !width.is_finite()
        || width <= 0.0
    {
        return Err(LaminarFilmError::InvalidWidth { section });
    }
    let flux = conditions.volume_flow_m3_s / width;
    if !flux.is_finite() || flux <= 0.0 {
        return Err(LaminarFilmError::Film {
            section,
            source: FilmError::NonRepresentable,
        });
    }
    NusseltFilm::new(
        conditions.material,
        flux,
        geometry.tangential_acceleration_m_s2,
        geometry.normal_acceleration_m_s2,
    )
    .map_err(|source| LaminarFilmError::Film { section, source })
}
