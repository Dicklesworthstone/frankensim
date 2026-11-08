//! SI finite-wavelength stability through the native free-surface spectral owner.

use super::{
    LaminarFilmError, LaminarFilmSection, LaminarGrowthObjective, LaminarPourConditions,
    section_film,
};
use fs_cheb::orr_sommerfeld::film::{
    FilmParameters, FilmSpectrumBudget, FilmSpectrumError, PreparedFilmSpectrum,
};
use fs_lbm::film::FilmError;

/// Sampled local free-surface spectra of the actual Newtonian pour path.
///
/// The native Nusselt base and all requested spectral modes are retained.
/// This admits finite kh and finite Re without the Benney expansion. Constant
/// Newtonian properties, a locally parallel planar film, no slip, negligible
/// gas stress and infinitesimal disturbances remain assumptions. A negative
/// sampled result is neither global stability nor a turbulence certificate.
#[derive(Clone, Copy, Debug)]
pub struct LaminarFilmSpectralStudy<'a> {
    pub sections: &'a [LaminarFilmSection],
    /// Positive angular wavenumbers in rad/m; no unsampled-band claim is made.
    pub wavenumbers_rad_m: &'a [f64],
    /// Both spectral orders, requested modes, work caps and accuracy tolerances.
    pub budget: FilmSpectrumBudget,
    /// Authored lateral-width applicability limit h/width, strictly in (0,1).
    /// This is an assumption budget, not a calibrated side-wall error estimate.
    pub max_thickness_width_ratio: f64,
}

/// Worst sampled temporal mode and the actual physical/numerical ranges.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarSpectralReport {
    pub max_growth_per_s: f64,
    pub section: usize,
    pub wavenumber_index: usize,
    pub spectral_mode: usize,
    pub wavenumber_rad_m: f64,
    /// Signed downstream angular frequency of the worst mode, in rad/s.
    pub angular_frequency_rad_s: f64,
    /// Complex eigenvalue change between orders at the worst mode, in 1/s.
    /// A numerical resolution indicator, not a rigorous error bound.
    pub order_difference_per_s: f64,
    /// Largest original-pencil residual over all returned modes and samples.
    pub max_relative_residual: f64,
    pub reynolds_range: [f64; 2],
    pub max_kh: f64,
    pub max_thickness_width_ratio: f64,
}

/// Original physical or spectral refusal, with the affected path/sample.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminarSpectralError {
    Film(LaminarFilmError),
    Spectrum {
        section: usize,
        wavenumber_index: usize,
        source: FilmSpectrumError,
    },
}
impl core::fmt::Display for LaminarSpectralError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Film(error) => error.fmt(f),
            Self::Spectrum {
                section,
                wavenumber_index,
                source,
            } => write!(
                f,
                "film section {section}, wavenumber {wavenumber_index}: {source}"
            ),
        }
    }
}
impl std::error::Error for LaminarSpectralError {}

impl LaminarFilmSpectralStudy<'_> {
    pub fn evaluate(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<LaminarSpectralReport, LaminarSpectralError> {
        self.evaluate_with_cancel(width_scale, conditions, || false)
    }

    /// Reuse one prepared pair of Chebyshev grids across every section and
    /// wavenumber. All results remain local until the complete study succeeds.
    pub fn evaluate_with_cancel(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminarSpectralReport, LaminarSpectralError> {
        let physical = LaminarSpectralError::Film;
        if self.sections.is_empty() {
            return Err(physical(LaminarFilmError::EmptyPath));
        }
        if self.wavenumbers_rad_m.is_empty() {
            return Err(physical(LaminarFilmError::EmptyModes));
        }
        if !positive(width_scale) {
            return Err(physical(LaminarFilmError::InvalidWidthScale));
        }
        if !positive(conditions.volume_flow_m3_s) {
            return Err(physical(LaminarFilmError::InvalidVolumeFlow));
        }
        if !positive(self.max_thickness_width_ratio) || self.max_thickness_width_ratio >= 1.0 {
            return Err(physical(LaminarFilmError::Film {
                section: 0,
                source: FilmError::InvalidInput("spectral film width budget"),
            }));
        }
        for (mode, &k) in self.wavenumbers_rad_m.iter().enumerate() {
            if !positive(k) {
                return Err(physical(LaminarFilmError::Mode {
                    section: 0,
                    mode,
                    source: FilmError::InvalidInput("positive spectral wavenumber"),
                }));
            }
        }
        let solver = PreparedFilmSpectrum::new(self.budget, &mut cancelled)
            .map_err(|source| spectral_error(0, 0, source))?;
        let mut report = LaminarSpectralReport {
            max_growth_per_s: f64::NEG_INFINITY,
            section: 0,
            wavenumber_index: 0,
            spectral_mode: 0,
            wavenumber_rad_m: 0.0,
            angular_frequency_rad_s: 0.0,
            order_difference_per_s: 0.0,
            max_relative_residual: 0.0,
            reynolds_range: [f64::INFINITY, 0.0],
            max_kh: 0.0,
            max_thickness_width_ratio: 0.0,
        };
        for (section, &geometry) in self.sections.iter().enumerate() {
            let film =
                section_film(geometry, width_scale, conditions, section).map_err(physical)?;
            let h = film.thickness_m();
            let speed = film.mean_velocity_m_s();
            let clock_rate = speed / h;
            let aspect = h / (geometry.width_m * width_scale);
            if !positive(clock_rate) || !positive(aspect) {
                return Err(physical(LaminarFilmError::Film {
                    section,
                    source: FilmError::NonRepresentable,
                }));
            }
            if aspect > self.max_thickness_width_ratio {
                return Err(physical(LaminarFilmError::Film {
                    section,
                    source: FilmError::InvalidInput("laterally wide spectral film"),
                }));
            }
            let normal_gravity = if geometry.normal_acceleration_m_s2 == 0.0 {
                0.0
            } else {
                geometry.normal_acceleration_m_s2 * h / (speed * speed)
            };
            let capillarity = if conditions.material.surface_tension_n_m == 0.0 {
                0.0
            } else {
                conditions.material.surface_tension_n_m
                    / (conditions.material.density_kg_m3 * speed * speed * h)
            };
            if !normal_gravity.is_finite()
                || !capillarity.is_finite()
                || (geometry.normal_acceleration_m_s2 > 0.0 && normal_gravity == 0.0)
                || (conditions.material.surface_tension_n_m > 0.0 && capillarity == 0.0)
            {
                return Err(physical(LaminarFilmError::Film {
                    section,
                    source: FilmError::NonRepresentable,
                }));
            }
            report.reynolds_range[0] = report.reynolds_range[0].min(film.reynolds());
            report.reynolds_range[1] = report.reynolds_range[1].max(film.reynolds());
            report.max_thickness_width_ratio = report.max_thickness_width_ratio.max(aspect);
            for (wavenumber_index, &k) in self.wavenumbers_rad_m.iter().enumerate() {
                let parameters = FilmParameters {
                    reynolds: film.reynolds(),
                    kh: k * h,
                    normal_gravity,
                    capillarity,
                };
                let modes = solver
                    .solve(parameters, &mut cancelled)
                    .map_err(|source| spectral_error(section, wavenumber_index, source))?;
                report.max_kh = report.max_kh.max(parameters.kh);
                for (spectral_mode, mode) in modes.iter().enumerate() {
                    let growth = clock_rate * mode.eigenvalue.re;
                    let frequency = -clock_rate * mode.eigenvalue.im;
                    let difference = clock_rate * mode.order_difference;
                    if !growth.is_finite() || !frequency.is_finite() || !difference.is_finite() {
                        return Err(physical(LaminarFilmError::Mode {
                            section,
                            mode: wavenumber_index,
                            source: FilmError::NonRepresentable,
                        }));
                    }
                    report.max_relative_residual =
                        report.max_relative_residual.max(mode.relative_residual);
                    if growth > report.max_growth_per_s {
                        report.max_growth_per_s = growth;
                        report.section = section;
                        report.wavenumber_index = wavenumber_index;
                        report.spectral_mode = spectral_mode;
                        report.wavenumber_rad_m = k;
                        report.angular_frequency_rad_s = frequency;
                        report.order_difference_per_s = difference;
                    }
                }
            }
        }
        if cancelled() {
            return Err(spectral_error(0, 0, FilmSpectrumError::Cancelled));
        }
        Ok(report)
    }
}

impl LaminarGrowthObjective for LaminarFilmSpectralStudy<'_> {
    /// The existing risk/racing error preserves the failing sample. Call
    /// evaluate/evaluate_with_cancel for the full native numerical diagnosis.
    fn growth_per_s(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<f64, LaminarFilmError> {
        self.evaluate(width_scale, conditions)
            .map(|report| report.max_growth_per_s)
            .map_err(|error| match error {
                LaminarSpectralError::Film(error) => error,
                LaminarSpectralError::Spectrum {
                    section,
                    wavenumber_index,
                    ..
                } => LaminarFilmError::Spectral {
                    section,
                    mode: wavenumber_index,
                },
            })
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn spectral_error(
    section: usize,
    wavenumber_index: usize,
    source: FilmSpectrumError,
) -> LaminarSpectralError {
    LaminarSpectralError::Spectrum {
        section,
        wavenumber_index,
        source,
    }
}
