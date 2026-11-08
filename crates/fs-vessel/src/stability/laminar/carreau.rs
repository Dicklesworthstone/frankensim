//! SI local stability of an original thermal Carreau-Yasuda film.
//!
//! Each actual section resolves thickness at its prescribed flux through the
//! native nonlinear film/slit owner. Both spectral grids sample that SAME law,
//! temperature, shear envelope and accepted profile. The original apparent and
//! differential viscosities enter different stress components. No fitted
//! viscosity, parabolic base, Benney expansion or temporal-to-spatial conversion
//! is introduced.

use super::{LaminarFilmError, LaminarFilmSection, LaminarSpectralReport};
use fs_cheb::orr_sommerfeld::film::generalized::FilmProfilePoint;
use fs_cheb::orr_sommerfeld::film::{
    FilmParameters, FilmSpectrumBudget, FilmSpectrumError, PreparedFilmSpectrum,
};
use fs_lbm::film::FilmError;
use fs_lbm::film::carreau::{CarreauFilm, CarreauFilmControls, CarreauFilmError, CarreauFilmSpec};
use fs_material::fluid_transport::rheology::ThermalCarreauYasuda;

/// Actual isothermal source and SI operating point, conserved along the path.
///
/// Density/capillarity are authored at this actual temperature. The original
/// constitutive declaration supplies viscosity and its shear derivative; it is
/// not a material-card provenance claim or a universal laminar-flow certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauPourConditions {
    pub volume_flow_m3_s: f64,
    pub density_kg_m3: f64,
    pub surface_tension_n_m: f64,
    pub law: ThermalCarreauYasuda,
    pub temperature_k: f64,
    /// Physical shear envelope, within the original source's actual T support.
    pub maximum_laminar_shear_rate_per_s: f64,
}

/// Sampled finite-wavelength temporal stability of a locally parallel film.
///
/// The strong variable-viscosity operator requires regular coefficient
/// derivatives at zero shear. Nonconstant Carreau sources currently admit
/// transition exponent a>=2; exact constant-law limits admit any source a.
/// Three-dimensional disturbances, inlet development, thermal coupling, gas
/// stress, side-wall corrections and contact lines are outside this local model.
#[derive(Clone, Copy, Debug)]
pub struct CarreauFilmSpectralStudy<'a> {
    pub sections: &'a [LaminarFilmSection],
    pub wavenumbers_rad_m: &'a [f64],
    pub budget: FilmSpectrumBudget,
    pub base_controls: CarreauFilmControls,
    /// Authored lateral-width applicability limit h/width, strictly in (0,1).
    pub max_thickness_width_ratio: f64,
}

/// Actual nonlinear base ranges and the worst accepted local temporal mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauFilmSpectralReport {
    /// The Reynolds range here uses the SAME original mu_0(T) as every
    /// normalized apparent/tangent viscosity. It is not a transition threshold.
    pub spectrum: LaminarSpectralReport,
    pub thickness_range_m: [f64; 2],
    /// Largest |actual-requested|/requested flux over physical sections.
    pub max_relative_flow_residual: f64,
    /// Largest native quadrature indicator in bases and point velocities.
    pub max_profile_quadrature_relative_change: f64,
}

/// Full original source/physics/spectral refusal, with no partial result.
#[derive(Debug)]
pub enum CarreauFilmSpectralError {
    Input(LaminarFilmError),
    UnsupportedSurfaceRegularity {
        transition_exponent: f64,
    },
    Base {
        section: usize,
        source: CarreauFilmError,
    },
    Spectrum {
        section: usize,
        wavenumber_index: usize,
        source: FilmSpectrumError,
    },
}
impl core::fmt::Display for CarreauFilmSpectralError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original-law film spectrum refused: {self:?}")
    }
}
impl std::error::Error for CarreauFilmSpectralError {}

impl CarreauFilmSpectralStudy<'_> {
    pub fn evaluate(
        &self,
        width_scale: f64,
        conditions: CarreauPourConditions,
    ) -> Result<CarreauFilmSpectralReport, CarreauFilmSpectralError> {
        self.evaluate_with_cancel(width_scale, conditions, || false)
    }

    /// Resolve each source profile once and reuse its original derivative
    /// grids across every wavelength. Results remain local until all physical,
    /// source, residual and paired-resolution checks have succeeded.
    pub fn evaluate_with_cancel(
        &self,
        width_scale: f64,
        conditions: CarreauPourConditions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CarreauFilmSpectralReport, CarreauFilmSpectralError> {
        let input = CarreauFilmSpectralError::Input;
        if self.sections.is_empty() {
            return Err(input(LaminarFilmError::EmptyPath));
        }
        if self.wavenumbers_rad_m.is_empty() {
            return Err(input(LaminarFilmError::EmptyModes));
        }
        if !positive(width_scale) {
            return Err(input(LaminarFilmError::InvalidWidthScale));
        }
        if !positive(conditions.volume_flow_m3_s) {
            return Err(input(LaminarFilmError::InvalidVolumeFlow));
        }
        if !positive(self.max_thickness_width_ratio) || self.max_thickness_width_ratio >= 1.0 {
            return Err(input(LaminarFilmError::Film {
                section: 0,
                source: FilmError::InvalidInput("Carreau spectral film width budget"),
            }));
        }
        for (mode, &k) in self.wavenumbers_rad_m.iter().enumerate() {
            if !positive(k) {
                return Err(input(LaminarFilmError::Mode {
                    section: 0,
                    mode,
                    source: FilmError::InvalidInput("positive spectral wavenumber"),
                }));
            }
        }
        let [mu0, mu_inf, time, index, transition] = conditions.law.curve().parameters();
        let constant = mu0 == mu_inf || time == 0.0 || index == 1.0;
        if !constant && transition < 2.0 {
            return Err(CarreauFilmSpectralError::UnsupportedSurfaceRegularity {
                transition_exponent: transition,
            });
        }
        let solver = PreparedFilmSpectrum::new(self.budget, &mut cancelled)
            .map_err(|source| spectral_error(0, 0, source))?;
        let mut report = CarreauFilmSpectralReport {
            spectrum: LaminarSpectralReport {
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
            },
            thickness_range_m: [f64::INFINITY, 0.0],
            max_relative_flow_residual: 0.0,
            max_profile_quadrature_relative_change: 0.0,
        };
        for (section, &geometry) in self.sections.iter().enumerate() {
            let width = geometry.width_m * width_scale;
            if !positive(geometry.width_m) || !positive(width) {
                return Err(input(LaminarFilmError::InvalidWidth { section }));
            }
            let spec = CarreauFilmSpec {
                density_kg_m3: conditions.density_kg_m3,
                surface_tension_n_m: conditions.surface_tension_n_m,
                flow_per_width_m2_s: conditions.volume_flow_m3_s / width,
                tangential_acceleration_m_s2: geometry.tangential_acceleration_m_s2,
                normal_acceleration_m_s2: geometry.normal_acceleration_m_s2,
            };
            let film = CarreauFilm::at_flow_with_cancel(
                conditions.law,
                conditions.temperature_k,
                spec,
                conditions.maximum_laminar_shear_rate_per_s,
                self.base_controls,
                &mut cancelled,
            )
            .map_err(|source| CarreauFilmSpectralError::Base { section, source })?;
            let h = film.thickness_m();
            let speed = film.mean_velocity_m_s();
            let clock_rate = speed / h;
            let aspect = h / width;
            if !positive(clock_rate) || !positive(aspect) {
                return Err(input(LaminarFilmError::Film {
                    section,
                    source: FilmError::NonRepresentable,
                }));
            }
            if aspect > self.max_thickness_width_ratio {
                return Err(input(LaminarFilmError::Film {
                    section,
                    source: FilmError::InvalidInput("laterally wide Carreau spectral film"),
                }));
            }
            let normal_gravity = if geometry.normal_acceleration_m_s2 == 0.0 {
                0.0
            } else {
                geometry.normal_acceleration_m_s2 * h / (speed * speed)
            };
            let capillarity = if conditions.surface_tension_n_m == 0.0 {
                0.0
            } else {
                conditions.surface_tension_n_m / (conditions.density_kg_m3 * speed * speed * h)
            };
            if !normal_gravity.is_finite()
                || !capillarity.is_finite()
                || (geometry.normal_acceleration_m_s2 > 0.0 && normal_gravity == 0.0)
                || (conditions.surface_tension_n_m > 0.0 && capillarity == 0.0)
            {
                return Err(input(LaminarFilmError::Film {
                    section,
                    source: FilmError::NonRepresentable,
                }));
            }
            let (coarse_nodes, fine_nodes) = solver.profile_coordinates();
            let (coarse, coarse_change) = sample_profile(&film, coarse_nodes, &mut cancelled)
                .map_err(|source| CarreauFilmSpectralError::Base { section, source })?;
            let (fine, fine_change) = sample_profile(&film, fine_nodes, &mut cancelled)
                .map_err(|source| CarreauFilmSpectralError::Base { section, source })?;
            // Every admitted nonconstant CY source has a>=2, so the exact
            // one-sided surface differential-viscosity gradient is zero.
            let prepared = solver
                .prepare_profile(&coarse, &fine, 0.0, &mut cancelled)
                .map_err(|source| spectral_error(section, 0, source))?;
            let spectrum = &mut report.spectrum;
            spectrum.reynolds_range[0] = spectrum.reynolds_range[0].min(film.reference_reynolds());
            spectrum.reynolds_range[1] = spectrum.reynolds_range[1].max(film.reference_reynolds());
            spectrum.max_thickness_width_ratio = spectrum.max_thickness_width_ratio.max(aspect);
            report.thickness_range_m[0] = report.thickness_range_m[0].min(h);
            report.thickness_range_m[1] = report.thickness_range_m[1].max(h);
            report.max_relative_flow_residual = report
                .max_relative_flow_residual
                .max((film.flow_residual_m2_s() / spec.flow_per_width_m2_s).abs());
            report.max_profile_quadrature_relative_change = report
                .max_profile_quadrature_relative_change
                .max(film.quadrature_relative_change())
                .max(coarse_change)
                .max(fine_change);
            for (wavenumber_index, &k) in self.wavenumbers_rad_m.iter().enumerate() {
                let parameters = FilmParameters {
                    reynolds: film.reference_reynolds(),
                    kh: k * h,
                    normal_gravity,
                    capillarity,
                };
                let modes = prepared
                    .solve(parameters, &mut cancelled)
                    .map_err(|source| spectral_error(section, wavenumber_index, source))?;
                spectrum.max_kh = spectrum.max_kh.max(parameters.kh);
                for (spectral_mode, mode) in modes.iter().enumerate() {
                    let growth = clock_rate * mode.eigenvalue.re;
                    let frequency = -clock_rate * mode.eigenvalue.im;
                    let difference = clock_rate * mode.order_difference;
                    if !growth.is_finite() || !frequency.is_finite() || !difference.is_finite() {
                        return Err(input(LaminarFilmError::Mode {
                            section,
                            mode: wavenumber_index,
                            source: FilmError::NonRepresentable,
                        }));
                    }
                    spectrum.max_relative_residual =
                        spectrum.max_relative_residual.max(mode.relative_residual);
                    if growth > spectrum.max_growth_per_s {
                        spectrum.max_growth_per_s = growth;
                        spectrum.section = section;
                        spectrum.wavenumber_index = wavenumber_index;
                        spectrum.spectral_mode = spectral_mode;
                        spectrum.wavenumber_rad_m = k;
                        spectrum.angular_frequency_rad_s = frequency;
                        spectrum.order_difference_per_s = difference;
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

fn sample_profile(
    film: &CarreauFilm,
    nodes: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<FilmProfilePoint>, f64), CarreauFilmError> {
    let h = film.thickness_m();
    let speed = film.mean_velocity_m_s();
    let mu0 = film.zero_shear_viscosity_pa_s();
    let mut profile = Vec::with_capacity(nodes.len());
    let mut change = 0.0_f64;
    for &y in nodes {
        let point = film.sample_with_cancel(y * h, &mut *cancelled)?;
        change = change.max(point.quadrature_relative_change);
        profile.push(FilmProfilePoint {
            velocity: point.velocity_m_s / speed,
            velocity_gradient: point.velocity_gradient_per_s * h / speed,
            velocity_curvature: point.velocity_curvature_per_m_s * h * h / speed,
            apparent_viscosity: point.dynamic_viscosity_pa_s / mu0,
            stress_tangent: point.stress_tangent_pa_s / mu0,
        });
    }
    Ok((profile, change))
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn spectral_error(
    section: usize,
    wavenumber_index: usize,
    source: FilmSpectrumError,
) -> CarreauFilmSpectralError {
    CarreauFilmSpectralError::Spectrum {
        section,
        wavenumber_index,
        source,
    }
}

#[cfg(test)]
mod tests;
