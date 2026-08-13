//! Characteristic-port composition: a viscothermal TMM driving-point
//! becomes a [`fs_vfit::discretize::DelayedFilter`].
//!
//! This file is the coupling-layer filling of that nameless port. A
//! bore, a muffler, an HVAC run, and a pulse tube are the same object.
use fs_duct::{Duct, DuctError, LossModel, Termination, impedance_sweep, input_impedance};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_vfit::FitOptions;
use fs_vfit::discretize::{DelayedFilter, DiscretizeError, RealizeError, reflectance};

/// Typed refusal from a driving-point characteristic line.
#[derive(Debug, Clone, PartialEq)]
pub enum DrivingPointError {
    /// Geometry or delay is not realizable.
    Invalid {
        /// Which check failed.
        what: &'static str,
    },
    /// TMM refusal.
    Duct(DuctError),
    /// Vector-fit / bilinear refusal.
    Realize(RealizeError),
    /// Discrete runtime refusal.
    Discrete(DiscretizeError),
}

impl core::fmt::Display for DrivingPointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "FS-COUPLE-DRIVE: {what}"),
            Self::Duct(e) => write!(f, "{e}"),
            Self::Realize(e) => write!(f, "{e}"),
            Self::Discrete(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DrivingPointError {}

/// Build a causal characteristic line from a TMM duct.
///
/// Sweeps `Z_in` (acoustic `e^{-iωt}`), forms `R = (Z−Zc)/(Z+Zc)`,
/// conjugates into the vector-fitting convention, peels the geometric
/// delay `2L/c`, and realizes the residual as a [`DelayedFilter`].
///
/// # Errors
/// Empty geometry, a delay that does not fit `n` samples, TMM, or fit.
pub fn characteristic_line(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
) -> Result<DelayedFilter, DrivingPointError> {
    let length: f64 = physics
        .segments
        .iter()
        .map(|s| match *s {
            fs_duct::Segment::Cylinder { length, .. } | fs_duct::Segment::Cone { length, .. } => {
                length
            }
            fs_duct::Segment::ToneHole { .. } => 0.0,
        })
        .sum();
    if !(length > 0.0) {
        return Err(DrivingPointError::Invalid {
            what: "line needs positive axial length",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let delay_samples = 2.0 * length / gas.sound_speed / dt;
    if !(delay_samples >= 2.0 && delay_samples < n as f64 - 2.0) {
        return Err(DrivingPointError::Invalid {
            what: "round-trip delay does not fit the realized history",
        });
    }
    let omega0 = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
    let omega_lo = 0.25 * omega0;
    let omega_hi = (12.0 * omega0).min(0.45 * core::f64::consts::PI / dt);
    if !(omega_hi > omega_lo) {
        return Err(DrivingPointError::Invalid {
            what: "impedance sweep band is empty at this sample rate",
        });
    }
    let sweep = impedance_sweep(
        physics,
        gas,
        omega_lo,
        omega_hi,
        96,
        LossModel::AllRegime,
        termination,
    )
    .map_err(DrivingPointError::Duct)?;
    let omega: Vec<f64> = sweep.iter().map(|r| r.omega).collect();
    // Duct is e^{-iωt}; vfit is e^{+iωt}. Conjugate, then peel delay.
    let r_fit: Vec<C64> = sweep
        .iter()
        .map(|r| {
            let rac = reflectance(C64::new(r.impedance.re, r.impedance.im), zc);
            C64::new(rac.re, -rac.im)
        })
        .collect();
    let mut opts = FitOptions::new(4);
    opts.fit_e = false;
    opts.iterations = 12;
    let mut line = DelayedFilter::from_tabulated(&omega, &r_fit, delay_samples, dt, &opts, omega0)
        .map_err(DrivingPointError::Realize)?;
    // Admit loop gain at the quarter-wave: peel the same delay from
    // the TMM reflectance and pin |H(ω0)| to that magnitude.
    let z0 = input_impedance(physics, gas, omega0, LossModel::AllRegime, termination)
        .map_err(DrivingPointError::Duct)?;
    let r0 = reflectance(z0.impedance, zc);
    let r_fit0 = C64::new(r0.re, -r0.im);
    let tau = delay_samples * dt;
    let peeled = r_fit0
        * C64::new(
            fs_math::det::cos(omega0 * tau),
            fs_math::det::sin(omega0 * tau),
        );
    line.pin_magnitude_at(omega0, peeled.abs().clamp(0.05, 0.99));
    Ok(line)
}
