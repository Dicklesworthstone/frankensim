//! Optional pressure reaction for the same planar, infinite-baffle image used
//! by the primary piano receiver. All structural columns and area samples feed
//! the existing passive load fitter and the existing mechanical air exchange.
use super::{board_geometry::SurfaceSample, engine::radiation::Model, linear::Bank};
use fs_bem::helmholtz::{Medium, rayleigh_modal::{AreaSample, PreparedSurface, MAX_MODAL_RAYLEIGH_FIELDS}};
use fs_couple::render::plate::impact::radiation::fit;
use fs_math::det;
use std::f64::consts::TAU;

/// An offline ceiling on cell-mode-frequency visits, not a time guarantee.
/// The default 3918-cell preset fits even at 32 modes and 65 frequencies.
const MAX_SWEEP_TERMS: u64 = 64_000_000_000;

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub sides: u8,
    pub min_hz: f64,
    pub max_hz: f64,
    pub samples: usize,
}

impl Options {
    /// sides,min_hz,max_hz,samples. The shared fitter owns its fixed dictionary;
    /// there is no user-selected alternate fitting order or relaxed error gate.
    pub fn parse(text: &str) -> Result<Self, String> {
        let f: Vec<_> = text.split(',').map(str::trim).collect();
        let error = || "--rayleigh-loading expects sides,min_hz,max_hz,samples; sides 1|2 and odd 33..257 samples".to_string();
        if f.len() != 4 { return Err(error()); }
        let result = Self {
            sides: f[0].parse().map_err(|_| error())?,
            min_hz: f[1].parse().map_err(|_| error())?,
            max_hz: f[2].parse().map_err(|_| error())?,
            samples: f[3].parse().map_err(|_| error())?,
        };
        result.validate()?;
        Ok(result)
    }

    fn validate(&self) -> Result<(), String> {
        if !matches!(self.sides, 1 | 2) || !self.min_hz.is_finite() || self.min_hz <= 0.0
            || !self.max_hz.is_finite() || self.max_hz <= self.min_hz
            || !(33..=257).contains(&self.samples) || self.samples % 2 == 0 {
            return Err("Rayleigh loading needs 1 or 2 half-spaces, a finite increasing positive band and odd 33..257 frequencies".into());
        }
        Ok(())
    }

    pub fn validate_clock(&self, output_rate: u32, mechanics_rate: u32) -> Result<(), String> {
        self.validate()?;
        if output_rate < 8_000 || mechanics_rate < output_rate || mechanics_rate % output_rate != 0
            || self.max_hz >= 0.45 * f64::from(output_rate) {
            return Err("Rayleigh loading must lie below 45% of the output rate with an integer mechanical clock".into());
        }
        // The SAME dictionary and off-band storage guard as loaded BEM playback.
        if 8.0 * TAU * self.max_hz >= 0.9 * std::f64::consts::PI * f64::from(mechanics_rate) {
            return Err("Rayleigh loading's 8x off-band acoustic poles exceed the mechanical-rate guard".into());
        }
        Ok(())
    }
}

pub struct Report {
    pub sides: u8,
    pub surface_samples: usize,
    pub board_modes: usize,
    pub frequencies: usize,
    pub min_hz: f64,
    pub max_hz: f64,
    pub reference_min_hz: f64,
    pub reference_max_hz: f64,
    pub acoustic_poles: usize,
    pub minimum_equivalent_disc_ppw: f64,
    pub sweep_terms: u64,
    pub peak_error: f64,
    pub rms_error: f64,
    pub resistance_peak_error: f64,
    pub resistance_rms_error: f64,
}

impl std::fmt::Display for Report {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Rayleigh pressure reaction: {} half-space(s), {} complete loaded modes, {} area samples, {} acoustic poles; {}..{} Hz at {} frequencies; complex peak/RMS {:.5}/{:.5}, resistance {:.5}/{:.5}. Loaded reference oscillators {:.3}..{:.3} Hz; minimum equivalent-disc cells/wavelength {:.2}; {} cell-mode-frequency visits. Sampled-fit estimates, no out-of-band or quadrature-error certificate.",
            self.sides, self.board_modes, self.surface_samples, self.acoustic_poles,
            self.min_hz, self.max_hz, self.frequencies, self.peak_error, self.rms_error,
            self.resistance_peak_error, self.resistance_rms_error, self.reference_min_hz,
            self.reference_max_hz, self.minimum_equivalent_disc_ppw, self.sweep_terms)
    }
}

pub struct PreparedLoad { pub model: Model, pub report: Report }

/// Prepare before attaching the existing radiation owner or launching a note.
/// `SurfaceSample` shapes are bare-board shapes: apply the bank's actual loaded
/// map once, just as the primary microphone does. The source must be planar in
/// XY; a projected representation of a crowned source must also be refused by
/// its geometry frontdoor. Two sides means two independent equal half-spaces,
/// never an extra receiver gain and never an inferred lower cavity.
pub fn prepare(surface: &[SurfaceSample], bank: &Bank, options: &Options, medium: Medium,
    output_rate: u32) -> Result<PreparedLoad, String> {
    options.validate_clock(output_rate, bank.rate)?;
    let modes = bank.board_count;
    if !(1..=MAX_MODAL_RAYLEIGH_FIELDS).contains(&modes) || surface.is_empty()
        || surface.len() > fs_bem::helmholtz::MAX_DENSE_PANELS {
        return Err("Rayleigh loading retains 1..32 complete board modes and 1..8192 area samples; reduce the supplied physical image explicitly".into());
    }
    let count = surface.len() as u64;
    let sweep_terms = count.checked_mul(count).and_then(|v| v.checked_mul(modes as u64))
        .and_then(|v| v.checked_mul(options.samples as u64)).ok_or("Rayleigh sweep work overflow")?;
    if sweep_terms > MAX_SWEEP_TERMS {
        return Err(format!("Rayleigh loading requires {sweep_terms} cell-mode-frequency visits, above its {MAX_SWEEP_TERMS} offline ceiling; no samples or modes were discarded"));
    }
    let (_, references) = bank.board_pressure_basis()?;
    let reference_min_hz = references.iter().map(|v| det::sqrt(*v) / TAU).fold(f64::INFINITY, f64::min);
    let reference_max_hz = references.iter().map(|v| det::sqrt(*v) / TAU).fold(0.0_f64, f64::max);
    if options.min_hz > reference_min_hz || options.max_hz < reference_max_hz {
        return Err(format!("Rayleigh sampled band must cover all loaded board reference oscillators ({reference_min_hz:.6}..{reference_max_hz:.6} Hz); these are not full coupled-instrument poles"));
    }
    let z = surface[0].position_m[2];
    let scale = surface.iter().map(|p| p.position_m[2].abs()).fold(1.0_f64, f64::max);
    let tolerance = 256.0 * f64::EPSILON * scale;
    if surface.iter().any(|p| p.position_m.iter().any(|v| !v.is_finite())
        || (p.position_m[2] - z).abs() > tolerance) {
        return Err("Rayleigh loading requires one flat XY source plane; supply a body BEM for curved geometry".into());
    }
    let cells = surface.iter().map(|p| Ok(AreaSample {
        position_m: [p.position_m[0], p.position_m[1]], area_m2: p.area_m2,
        mode_shape: bank.project_board_shape(&p.mode_shape)?,
    })).collect::<Result<Vec<_>, String>>()?;
    let prepared = PreparedSurface::new(&cells).map_err(|e| e.to_string())?;
    let minimum_equivalent_disc_ppw = prepared.equivalent_disc_ppw(TAU * options.max_hz, medium)
        .map_err(|e| e.to_string())?;
    if minimum_equivalent_disc_ppw < fs_bem::helmholtz::MIN_PANELS_PER_WAVELENGTH {
        return Err(fs_bem::helmholtz::HelmholtzError::TooCoarse {
            panels_per_wavelength: minimum_equivalent_disc_ppw,
        }.to_string());
    }
    let low = det::ln(options.min_hz);
    let span = det::ln(options.max_hz) - low;
    let mut omega: Vec<_> = (0..options.samples).map(|i|
        TAU * det::exp(low + span * i as f64 / (options.samples - 1) as f64)).collect();
    omega[0] = TAU * options.min_hz;
    omega[options.samples - 1] = TAU * options.max_hz;
    let mut data = Vec::with_capacity(options.samples);
    for &w in &omega {
        let mut matrix = prepared.impedance(w, medium).map_err(|e| e.to_string())?;
        for value in &mut matrix { *value = value.scale(f64::from(options.sides)); }
        data.push(matrix);
    }
    let fitted = fit::fit(&omega, &data, modes)?;
    let report = Report {
        sides: options.sides, surface_samples: surface.len(), board_modes: modes,
        frequencies: options.samples, min_hz: options.min_hz, max_hz: options.max_hz,
        reference_min_hz, reference_max_hz, acoustic_poles: fitted.model.poles.len(),
        minimum_equivalent_disc_ppw, sweep_terms, peak_error: fitted.peak_error,
        rms_error: fitted.rms_error, resistance_peak_error: fitted.resistance_peak_error,
        resistance_rms_error: fitted.resistance_rms_error,
    };
    Ok(PreparedLoad { model: fitted.model, report })
}

#[cfg(test)]
#[path = "rayleigh_loading_tests.rs"]
mod tests;
