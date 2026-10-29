//! Continuous passive piston memory for a joint structural/acoustic solve.
//!
//! Geometry comes from the sole Bessel–Struve piston owner. The shared resolved
//! PSD-residue fitter supplies acoustic oscillators; this module performs no
//! integration and never reinterprets a digital filter's history as energy.
use crate::acoustic_realize::AcousticRealizeError;
use crate::render::plate::impact::radiation::{Model, fit};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::c64::C64;

/// Measured fit discrepancies, including frequency-local weak resistance.
/// Global errors use the reference peak or RMS on each disjoint sample grid.
/// Pointwise errors use the actual reference at each frequency. These are
/// sampled discrepancies, not continuum or closed-loop error certificates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonRadiationFit {
    pub complex_peak_relative_error: f64,
    pub complex_rms_relative_error: f64,
    pub resistance_peak_relative_error: f64,
    pub resistance_rms_relative_error: f64,
    pub complex_pointwise_relative_error: f64,
    pub resistance_pointwise_relative_error: f64,
    pub sample_points: usize,
    pub poles: usize,
}

/// Uniform circular piston in an infinite rigid baffle and homogeneous gas.
///
/// Its scalar port is VOLUME FLOW [m^3/s] -> pressure [Pa]. Connect a mechanical
/// velocity through its signed volume shape and apply the opposite pressure
/// force through the SAME shape. The point-force participation is independent.
/// The model has passive storage at every frequency; approximation accuracy is
/// checked only over the explicit source band. This is neither an unbaffled
/// source nor the full spatial radiation of a nonuniform vibrating face.
///
/// Distinct physical faces require distinct models unless their producer gives
/// a mutual impedance. Multiple coordinates of one uniform face may share its
/// ONE memory through summed volume flow. Equal face areas do not establish
/// common geometric identity. The caller owns closed-loop/band applicability.
#[derive(Clone, Debug)]
pub struct BaffledPistonRadiation {
    radius_m: f64,
    gas: GasState,
    band_rad_s: [f64; 2],
    model: Model,
    report: PistonRadiationFit,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl BaffledPistonRadiation {
    /// Fit the same continuous piston physics with bounded, geometry-resolving
    /// cold work. No coefficient is altered after the shared fit validates it.
    ///
    /// `sample_points` is odd, 33..=257. Even samples identify the model; odd
    /// samples only validate it. A quadratic frequency grid resolves weak
    /// low-frequency resistance as well as the diameter's travel-time detail.
    /// The shared `fit_resolved_relative` owner requires at most 2% complex and
    /// 2% resistive relative error at EVERY sample, in addition to its original
    /// global peak/RMS gates. These are explicit sampled approximation limits.
    ///
    /// The dictionary resolves at least eight intervals across the band and
    /// six per inverse diameter travel time, with logarithmic enrichment of
    /// its first interval. The training grid must resolve that spacing too.
    /// `maximum_poles`, 1..=128, bounds the COMPLETE resulting history; overflow
    /// refuses without truncation. Caller state and digital histories are never
    /// modified. The caller explicitly chooses initial acoustic q/p coordinates.
    /// Cancellation is checked between source samples and around one bounded fit.
    pub fn new(
        radius_m: f64,
        gas: &GasState,
        band_rad_s: [f64; 2],
        sample_points: usize,
        maximum_poles: usize,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        gas.validate_coherence()
            .map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
        let [low, high] = band_rad_s;
        if !radius_m.is_finite()
            || radius_m <= 0.0
            || !low.is_finite()
            || low <= 0.0
            || !high.is_finite()
            || high <= low
            || !(33..=257).contains(&sample_points)
            || sample_points % 2 == 0
            || !(1..=128).contains(&maximum_poles)
        {
            return Err(invalid(
                "piston radiation needs positive geometry, ordered positive band, odd 33..257 samples and 1..128 poles",
            ));
        }
        let width = high - low;
        let spacing =
            (width / 8.0).min(core::f64::consts::TAU * gas.sound_speed / (12.0 * radius_m));
        if !spacing.is_finite() || spacing <= 0.0 {
            return Err(invalid(
                "piston source travel-time spacing is unrepresentable",
            ));
        }
        let mut omega = Vec::with_capacity(sample_points);
        let mut data = Vec::with_capacity(sample_points);
        for i in 0..sample_points {
            if gate.is_requested() {
                return Err(AcousticRealizeError::Cancelled);
            }
            let fraction = i as f64 / (sample_points - 1) as f64;
            let w = if i + 1 == sample_points {
                high
            } else {
                low + width * fraction * fraction
            };
            let (r, x) =
                fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, radius_m, w)
                    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            omega.push(w);
            data.push(vec![C64::new(r, x)]);
        }
        let relative = fit::fit_resolved_relative(&omega, &data, 1, spacing)
            .map_err(AcousticRealizeError::Nonlinear)?;
        if gate.is_requested() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let fitted = relative.fit;
        if fitted.model.poles.len() > maximum_poles {
            return Err(invalid(
                "complete piston radiation history exceeds its declared pole capacity",
            ));
        }
        let report = PistonRadiationFit {
            complex_peak_relative_error: fitted.peak_error,
            complex_rms_relative_error: fitted.rms_error,
            resistance_peak_relative_error: fitted.resistance_peak_error,
            resistance_rms_relative_error: fitted.resistance_rms_error,
            complex_pointwise_relative_error: relative.complex_peak_relative_error,
            resistance_pointwise_relative_error: relative.resistance_peak_relative_error,
            sample_points,
            poles: fitted.model.poles.len(),
        };
        Ok(Self {
            radius_m,
            gas: *gas,
            band_rad_s,
            model: fitted.model,
            report,
        })
    }

    #[must_use]
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    #[must_use]
    pub const fn band_rad_s(&self) -> [f64; 2] {
        self.band_rad_s
    }
    #[must_use]
    pub const fn report(&self) -> PistonRadiationFit {
        self.report
    }
    /// Original scalar pressure/volume-flow model, with every fitted pole.
    #[must_use]
    pub const fn model(&self) -> &Model {
        &self.model
    }

    /// Fitted continuous p/U under exp(-i omega t). Queries outside the source
    /// band refuse, including DC. The passive continuation of the model alone
    /// supplies no geometry-fidelity evidence outside the fitted band.
    pub fn impedance(&self, omega_rad_s: f64) -> Result<C64, AcousticRealizeError> {
        if !omega_rad_s.is_finite()
            || omega_rad_s < self.band_rad_s[0]
            || omega_rad_s > self.band_rad_s[1]
        {
            return Err(invalid(
                "piston radiation frequency is outside the admitted source band",
            ));
        }
        self.model
            .impedance(omega_rad_s)
            .map(|z| z[0])
            .map_err(AcousticRealizeError::Nonlinear)
    }

    /// Realize the unchanged acoustic oscillators in the existing pHS owner.
    /// State order is interleaved acoustic q/p in `model().poles` order. Input
    /// is volume flow and output is pressure. The host's complete exact or
    /// Gonzalez solve owns their time evolution and stored-energy account.
    pub fn port_hamiltonian(&self) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
        let omega: Vec<_> = self.model.poles.iter().map(|p| p.omega).collect();
        let zeta: Vec<_> = self.model.poles.iter().map(|p| p.zeta).collect();
        let drive: Vec<_> = self.model.poles.iter().map(|p| p.coupling[0]).collect();
        fs_phs::modal_bank_ports(&omega, &zeta, &[&drive])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
    }
}

#[cfg(test)]
mod tests;
