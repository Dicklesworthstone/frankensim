//! Coherent vector diffraction through the original source-backed planar films.
//!
//! Samples describe one PERIODIC, band-limited transverse cell, with x varying
//! fastest and exp(+i omega t) peak phasors in V/m. The supplied Ex/Ey are the
//! incident field at the first film face. Every Fourier mode keeps its complex
//! phase and its own TE/TM stack response. Exact kz propagation, Ez and all H
//! components follow Maxwell's equations; there is no paraxial approximation.
//!
//! The entire incident grid must lie strictly inside the incident light cone.
//! Evanescent transmission and tunnelling remain supported by the original
//! films. An exactly grazing exit refuses because tangential E alone does not
//! parameterize that mode. No spectral bins are discarded or windowed. Padding,
//! spatial bandwidth and periodic-image convergence are caller-owned physical
//! resolution choices, not a certificate for an isolated hard-edged aperture.
//!
//! Reference: ETH Photonics, Angular Spectrum Representation, equations 7.5–7.9:
//! <https://ethz.ch/content/dam/ethz/special-interest/itet/photonics-dam/documents/lectures/EandM/AngularSpectrumRepresentation.pdf>.
//! The phasor sign here is the conjugate of that reference, matching `layers`.

use super::{Polarization, ReflectionError, layers::PlanarDielectricLayers};
use fs_fft::{C64 as FftComplex, Fft};
use fs_material::field_storage::{linear_electric_permittivity, linear_magnetic_permeability};
use fs_math::{c64::C64, det};
use modes::{basis, cartesian, fields, finite, normal_flux, polarized, propagation};

mod modes;

/// Bounds the uninterruptible work of one existing FFT pencil.
pub const MAX_AXIS_POINTS: usize = 4096;

/// Periodic cell samples at (column*dx, row*dy), with no implicit centering.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumGrid {
    pub width: usize,
    pub height: usize,
    pub spacing_m: [f64; 2],
}

/// Logical numerical storage includes the retained modes, plans, one returned
/// six-component field and its transform scratch. Input/source storage and
/// allocator overhead are outside this bound. FFT work counts complex
/// point-stages; one field query costs 6*N*(log2(width)+log2(height)).
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumLimits {
    pub maximum_grid_points: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

#[derive(Debug)]
pub enum AngularSpectrumError {
    Source(ReflectionError),
    Invalid(&'static str),
    IncidentLightCone { column: usize, row: usize },
    GrazingExit { column: usize, row: usize },
    Cancelled,
}
impl From<ReflectionError> for AngularSpectrumError {
    fn from(value: ReflectionError) -> Self {
        Self::Source(value)
    }
}
impl core::fmt::Display for AngularSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar angular spectrum: {self:?}")
    }
}
impl std::error::Error for AngularSpectrumError {}

/// Integrated cycle-mean powers over the declared cell, in watts. Reflected
/// power is a positive outward magnitude. Transmission is at the exit face;
/// finite-film absorption excludes subsequent loss in the exit half-space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumPower {
    pub incident_w: f64,
    pub reflected_w: f64,
    pub transmitted_w: f64,
    pub film_absorbed_w: f64,
}
impl AngularSpectrumPower {
    pub fn balance_defect_w(self) -> f64 {
        self.incident_w - self.reflected_w - self.transmitted_w - self.film_absorbed_w
    }
}

/// Distance is measured OUTWARD from the corresponding stack face. Incident
/// and reflected fields at distance d occupy the same plane z=-d; transmitted
/// fields occupy z=last_face+d. TotalIncidentSide retains standing-wave terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AngularSpectrumSide {
    Incident,
    Reflected,
    TotalIncidentSide,
    Transmitted,
}

/// Complex peak fields on the original cell. Poynting flux uses E cross H*,
/// including all spatial/polarization interference, with its physical sign.
#[derive(Debug)]
pub struct AngularSpectrumField {
    grid: AngularSpectrumGrid,
    electric_v_per_m: Vec<[C64; 3]>,
    magnetic_a_per_m: Vec<[C64; 3]>,
    normal_power_w: f64,
}
impl AngularSpectrumField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }
    pub fn electric_v_per_m(&self) -> &[[C64; 3]] {
        &self.electric_v_per_m
    }
    pub fn magnetic_a_per_m(&self) -> &[[C64; 3]] {
        &self.magnetic_a_per_m
    }
    pub fn normal_power_w(&self) -> f64 {
        self.normal_power_w
    }
    /// Cycle-mean local normal Poynting flux in W/m², row-major like the field.
    pub fn normal_flux_w_m2(&self, index: usize) -> Option<f64> {
        Some(normal_flux(
            *self.electric_v_per_m.get(index)?,
            *self.magnetic_a_per_m.get(index)?,
        ))
    }
}

struct Mode {
    kx: f64,
    ky: f64,
    incident_kz: f64,
    exit_kz: C64,
    incident_e: [C64; 2],
    reflected_e: [C64; 2],
    transmitted_e: [C64; 2],
}

/// One prepared coherent beam and its original material/film owner. The two
/// input FFTs and all angle-dependent film solves happen once. Each additional
/// plane needs six separable FFTs and O(N) propagation/reconstruction work;
/// neither the material laws nor the film recursion are reevaluated there.
pub struct PlanarAngularSpectrum<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    grid: AngularSpectrumGrid,
    x_fft: Fft,
    y_fft: Fft,
    modes: Vec<Mode>,
    incident_omega_mu: f64,
    exit_omega_mu: f64,
    power: AngularSpectrumPower,
    workspace_bytes: usize,
    field_fft_work: u64,
}

impl<'s, 'm> PlanarAngularSpectrum<'s, 'm> {
    pub fn prepare(
        source: &'s PlanarDielectricLayers<'m>,
        grid: AngularSpectrumGrid,
        incident_tangential_e_v_per_m: &[[C64; 2]],
        limits: AngularSpectrumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(&mut cancelled)?;
        let (n, area, workspace_bytes, field_fft_work) = admit(grid, limits)?;
        if incident_tangential_e_v_per_m.len() != n {
            return Err(bad("finite Ex/Ey samples must cover the declared grid"));
        }
        for chunk in incident_tangential_e_v_per_m.chunks(256) {
            poll(&mut cancelled)?;
            if !chunk.iter().flatten().all(|&e| finite(e)) {
                return Err(bad("finite Ex/Ey samples must cover the declared grid"));
            }
        }
        let omega = source.angular_frequency_rad_s();
        let epsilon = linear_electric_permittivity(source.incident_material())
            .map_err(ReflectionError::from)?;
        let mu = linear_magnetic_permeability(source.incident_material())
            .map_err(ReflectionError::from)?;
        let exit_mu =
            linear_magnetic_permeability(source.exit_material()).map_err(ReflectionError::from)?;
        let k = omega * det::sqrt(epsilon) * det::sqrt(mu);
        let incident_omega_mu = omega * mu;
        let exit_omega_mu = omega * exit_mu;
        if [k, incident_omega_mu, exit_omega_mu]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(bad("wave number or magnetic admittance is unrepresentable"));
        }
        // Admit ALL bins before allocating FFT plans or evaluating the films.
        // Real transverse wave numbers are conserved across every interface.
        let max_kx = if grid.width == 1 {
            0.0
        } else {
            core::f64::consts::PI / grid.spacing_m[0]
        };
        let max_ky = if grid.height == 1 {
            0.0
        } else {
            core::f64::consts::PI / grid.spacing_m[1]
        };
        if det::hypot(max_kx / k, max_ky / k) >= 1.0 {
            return Err(AngularSpectrumError::IncidentLightCone {
                column: grid.width / 2,
                row: grid.height / 2,
            });
        }
        let x_fft = Fft::new(grid.width);
        poll(&mut cancelled)?;
        let y_fft = Fft::new(grid.height);
        let mut ex = Vec::with_capacity(n);
        let mut ey = Vec::with_capacity(n);
        for (index, sample) in incident_tangential_e_v_per_m.iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            ex.push(fft_complex(sample[0]));
            ey.push(fft_complex(sample[1]));
        }
        let mut line = vec![FftComplex::default(); grid.width.max(grid.height)];
        let mut scratch = line.clone();
        // Coefficients multiply exp(-i k_t.r_t), so fs-fft INVERSE analyzes the
        // samples, including its 1/N normalization. FORWARD synthesizes them.
        transform(
            grid,
            &x_fft,
            &y_fft,
            &mut ex,
            &mut line,
            &mut scratch,
            true,
            &mut cancelled,
        )?;
        transform(
            grid,
            &x_fft,
            &y_fft,
            &mut ey,
            &mut line,
            &mut scratch,
            true,
            &mut cancelled,
        )?;
        let mut modes = Vec::with_capacity(n);
        let mut power = AngularSpectrumPower {
            incident_w: 0.0,
            reflected_w: 0.0,
            transmitted_w: 0.0,
            film_absorbed_w: 0.0,
        };
        for row in 0..grid.height {
            for column in 0..grid.width {
                let index = row * grid.width + column;
                if index % 32 == 0 {
                    poll(&mut cancelled)?;
                }
                let kx = transverse_number(column, grid.width, grid.spacing_m[0]);
                let ky = transverse_number(row, grid.height, grid.spacing_m[1]);
                let sine = det::hypot(kx / k, ky / k);
                let cosine = det::sqrt((1.0 - sine) * (1.0 + sine));
                if !cosine.is_finite() || cosine <= 0.0 {
                    return Err(AngularSpectrumError::IncidentLightCone { column, row });
                }
                let te = source.response(cosine, Polarization::Te)?;
                let tm = source.response(cosine, Polarization::Tm)?;
                let exit_kz = te.exit_normal_wavenumber_rad_m;
                if exit_kz == C64::ZERO {
                    return Err(AngularSpectrumError::GrazingExit { column, row });
                }
                let incident_e = [complex(ex[index]), complex(ey[index])];
                let p = basis(kx, ky);
                let a = polarized(incident_e, p);
                let reflected_e = cartesian(
                    [
                        a[0] * te.reflected_tangential_e,
                        a[1] * tm.reflected_tangential_e,
                    ],
                    p,
                );
                let transmitted_e = cartesian(
                    [
                        a[0] * te.transmitted_tangential_e,
                        a[1] * tm.transmitted_tangential_e,
                    ],
                    p,
                );
                if !incident_e
                    .iter()
                    .chain(&reflected_e)
                    .chain(&transmitted_e)
                    .all(|&z| finite(z))
                {
                    return Err(bad("transformed beam amplitudes are unrepresentable"));
                }
                let incident_kz = te.incident_normal_wavenumber_rad_m;
                let te_power = (area * 0.5) * a[0].norm_sq() * (incident_kz / incident_omega_mu);
                let tm_power = (area * 0.5) * a[1].norm_sq() * (omega * epsilon / incident_kz);
                power.incident_w += te_power + tm_power;
                power.reflected_w +=
                    te_power * te.reflected_power_fraction + tm_power * tm.reflected_power_fraction;
                power.transmitted_w += te_power * te.transmitted_power_fraction
                    + tm_power * tm.transmitted_power_fraction;
                power.film_absorbed_w +=
                    te_power * te.absorbed_power_fraction + tm_power * tm.absorbed_power_fraction;
                modes.push(Mode {
                    kx,
                    ky,
                    incident_kz,
                    exit_kz,
                    incident_e,
                    reflected_e,
                    transmitted_e,
                });
            }
        }
        if [
            power.incident_w,
            power.reflected_w,
            power.transmitted_w,
            power.film_absorbed_w,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(bad("integrated beam power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            grid,
            x_fft,
            y_fft,
            modes,
            incident_omega_mu,
            exit_omega_mu,
            power,
            workspace_bytes,
            field_fft_work,
        })
    }

    pub fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.source
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }
    pub fn power(&self) -> AngularSpectrumPower {
        self.power
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn field_fft_work(&self) -> u64 {
        self.field_fft_work
    }

    /// Reconstruct a complete complex vector field. Cancellation or any
    /// unrepresentable mode leaves the prepared beam untouched and returns no
    /// partial field. Propagation into an absorbing exit only decays.
    pub fn field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !distance_from_face_m.is_finite() || distance_from_face_m < 0.0 {
            return Err(bad("field distance must be finite and nonnegative"));
        }
        let n = self.modes.len();
        let mut electric = Vec::with_capacity(n);
        let mut magnetic = Vec::with_capacity(n);
        for (index, mode) in self.modes.iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let mut e = [C64::ZERO; 3];
            let mut h = [C64::ZERO; 3];
            let mut add =
                |kz, omega_mu, transverse, distance| -> Result<(), AngularSpectrumError> {
                    let phase = propagation(kz, distance)
                        .ok_or_else(|| bad("propagation phase is unrepresentable"))?;
                    let (em, hm) = fields(mode.kx, mode.ky, kz, omega_mu, transverse)
                        .ok_or_else(|| bad("longitudinal or magnetic field is unrepresentable"))?;
                    for c in 0..3 {
                        e[c] = e[c] + em[c] * phase;
                        h[c] = h[c] + hm[c] * phase;
                    }
                    Ok(())
                };
            if matches!(
                side,
                AngularSpectrumSide::Incident | AngularSpectrumSide::TotalIncidentSide
            ) {
                add(
                    C64::from_re(mode.incident_kz),
                    self.incident_omega_mu,
                    mode.incident_e,
                    -distance_from_face_m,
                )?;
            }
            if matches!(
                side,
                AngularSpectrumSide::Reflected | AngularSpectrumSide::TotalIncidentSide
            ) {
                add(
                    C64::from_re(-mode.incident_kz),
                    self.incident_omega_mu,
                    mode.reflected_e,
                    -distance_from_face_m,
                )?;
            }
            if side == AngularSpectrumSide::Transmitted {
                add(
                    mode.exit_kz,
                    self.exit_omega_mu,
                    mode.transmitted_e,
                    distance_from_face_m,
                )?;
            }
            electric.push(e);
            magnetic.push(h);
        }
        let mut buffer = vec![FftComplex::default(); n];
        let mut line = vec![FftComplex::default(); self.grid.width.max(self.grid.height)];
        let mut scratch = line.clone();
        for field in [&mut electric, &mut magnetic] {
            for c in 0..3 {
                for (index, (dst, src)) in buffer.iter_mut().zip(field.iter()).enumerate() {
                    if index % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    *dst = fft_complex(src[c]);
                }
                transform(
                    self.grid,
                    &self.x_fft,
                    &self.y_fft,
                    &mut buffer,
                    &mut line,
                    &mut scratch,
                    false,
                    &mut cancelled,
                )?;
                for (index, (dst, src)) in field.iter_mut().zip(&buffer).enumerate() {
                    if index % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    dst[c] = complex(*src);
                }
            }
        }
        let mut flux = 0.0;
        for (index, (&e, &h)) in electric.iter().zip(&magnetic).enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            if !e.iter().chain(&h).all(|&z| finite(z)) {
                return Err(bad("synthesized field is unrepresentable"));
            }
            flux += normal_flux(e, h);
        }
        let normal_power_w = flux * self.grid.spacing_m[0] * self.grid.spacing_m[1];
        if !normal_power_w.is_finite() {
            return Err(bad("synthesized Poynting power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(AngularSpectrumField {
            grid: self.grid,
            electric_v_per_m: electric,
            magnetic_a_per_m: magnetic,
            normal_power_w,
        })
    }
}

fn bad(what: &'static str) -> AngularSpectrumError {
    AngularSpectrumError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AngularSpectrumError> {
    if cancelled() {
        Err(AngularSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}
fn complex(z: FftComplex) -> C64 {
    C64::new(z.re, z.im)
}
fn fft_complex(z: C64) -> FftComplex {
    FftComplex::new(z.re, z.im)
}

fn transverse_number(index: usize, count: usize, spacing: f64) -> f64 {
    let bin = if index <= (count - 1) / 2 {
        index as f64
    } else {
        index as f64 - count as f64
    };
    (core::f64::consts::TAU * (bin / count as f64)) / spacing
}

fn admit(
    grid: AngularSpectrumGrid,
    limits: AngularSpectrumLimits,
) -> Result<(usize, f64, usize, u64), AngularSpectrumError> {
    if [grid.width, grid.height]
        .iter()
        .any(|&n| !n.is_power_of_two() || n > MAX_AXIS_POINTS)
        || grid.spacing_m.iter().any(|&v| !v.is_finite() || v <= 0.0)
    {
        return Err(bad(
            "grid axes must be powers of two up to 4096, with positive finite spacing",
        ));
    }
    let n = grid
        .width
        .checked_mul(grid.height)
        .ok_or_else(|| bad("grid extent overflow"))?;
    if n > limits.maximum_grid_points {
        return Err(bad("grid exceeds its point budget"));
    }
    let area = (grid.width as f64 * grid.spacing_m[0]) * (grid.height as f64 * grid.spacing_m[1]);
    if !area.is_finite() || area <= 0.0 {
        return Err(bad("cell area is unrepresentable"));
    }
    // Query storage dominates preparation's two input component arrays. Plan
    // allowance is conservative; Fft stores at most (width+height)/2 twiddles.
    let per_point = core::mem::size_of::<Mode>()
        + 6 * core::mem::size_of::<C64>()
        + core::mem::size_of::<FftComplex>();
    let workspace_bytes = n
        .checked_mul(per_point)
        .and_then(|v| {
            (2 * grid.width.max(grid.height) + grid.width + grid.height)
                .checked_mul(core::mem::size_of::<FftComplex>())
                .and_then(|s| v.checked_add(s))
        })
        .ok_or_else(|| bad("workspace extent overflow"))?;
    let field_fft_work = (n as u64)
        .checked_mul(6 * u64::from(grid.width.ilog2() + grid.height.ilog2()))
        .ok_or_else(|| bad("FFT work overflow"))?;
    if workspace_bytes > limits.maximum_workspace_bytes || field_fft_work > limits.maximum_fft_work
    {
        return Err(bad(
            "beam preparation/field query exceeds its storage or FFT budget",
        ));
    }
    Ok((n, area, workspace_bytes, field_fft_work))
}

#[allow(clippy::too_many_arguments)]
fn transform(
    grid: AngularSpectrumGrid,
    x: &Fft,
    y: &Fft,
    data: &mut [FftComplex],
    line: &mut [FftComplex],
    scratch: &mut [FftComplex],
    inverse: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    for row in data.chunks_exact_mut(grid.width) {
        poll(cancelled)?;
        if inverse {
            x.inverse(row, &mut scratch[..grid.width]);
        } else {
            x.forward(row, &mut scratch[..grid.width]);
        }
    }
    for column in 0..grid.width {
        poll(cancelled)?;
        for row in 0..grid.height {
            line[row] = data[row * grid.width + column];
        }
        if inverse {
            y.inverse(&mut line[..grid.height], &mut scratch[..grid.height]);
        } else {
            y.forward(&mut line[..grid.height], &mut scratch[..grid.height]);
        }
        for row in 0..grid.height {
            data[row * grid.width + column] = line[row];
        }
    }
    Ok(())
}
