//! Characteristic-port composition: a viscothermal TMM driving-point
//! becomes a [`fs_vfit::discretize::DelayedFilter`].
//!
//! This file is the coupling-layer filling of that nameless port. A
//! bore, a muffler, an HVAC run, and a pulse tube are the same object.
use fs_duct::{Duct, DuctError, LossModel, Termination, input_impedance};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_math::c64::C64;
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
/// Samples `R(ω) = (Z−Zc)/(Z+Zc)` on the DFT grid of the run and
/// inverse-transforms it. The resulting FIR *is* the linear
/// scattering port: a tone hole, a cone, and a pulse tube change
/// `R` and therefore the impulse. There is no rational-fit residual
/// and no one-pole fallback.
///
/// # Errors
/// Empty geometry, a delay that does not fit `n` samples, TMM, or DFT.
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
    let geo_delay = 2.0 * length / gas.sound_speed / dt;
    if !(geo_delay >= 2.0 && geo_delay < n as f64 - 2.0) {
        return Err(DrivingPointError::Invalid {
            what: "round-trip delay does not fit the realized history",
        });
    }
    let n_fft = ((4.0 * geo_delay).ceil() as usize)
        .next_power_of_two()
        .clamp(256, 4096);
    let fft = Fft::new(n_fft);
    let mut buf = vec![FftC64::new(0.0, 0.0); n_fft];
    buf[0] = match termination {
        Termination::Closed => FftC64::new(1.0, 0.0),
        Termination::IdealOpen | Termination::UnflangedOpen | Termination::FlangedOpen => {
            FftC64::new(-1.0, 0.0)
        }
    };
    for k in 1..=n_fft / 2 {
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let response = input_impedance(physics, gas, omega, LossModel::AllRegime, termination)
            .map_err(DrivingPointError::Duct)?;
        let rac = reflectance(
            C64::new(response.impedance.re, response.impedance.im),
            zc,
        );
        buf[k] = FftC64::new(rac.re, rac.im);
        if k != n_fft / 2 {
            buf[n_fft - k] = FftC64::new(rac.re, -rac.im);
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.inverse(&mut buf, &mut scratch);
    let ir: Vec<f64> = buf.iter().map(|c| c.re).collect();
    DelayedFilter::from_impulse_response(dt, ir).map_err(DrivingPointError::Discrete)
}
