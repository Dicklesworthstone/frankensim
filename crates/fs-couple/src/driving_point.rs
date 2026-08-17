//! Characteristic-port composition: a viscothermal TMM driving-point
//! becomes a [`fs_vfit::discretize::DelayedFilter`].
//!
//! This file is the coupling-layer filling of that nameless port. A
//! bore, a muffler, an HVAC run, and a pulse tube are the same object.
use fs_duct::{Duct, DuctError, LossModel, Termination, input_impedance_wall};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_phs::WallPin;
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
    wall: Option<&WallPin>,
) -> Result<DelayedFilter, DrivingPointError> {
    characteristic_line_dense(physics, gas, termination, sample_rate_hz, n, zc, wall, None)
}

/// [`characteristic_line`] with an explicit realization length.
///
/// `n_fft_override` sets the reflectance sampling density (power of
/// two, clamped to [256, 8192]); `None` keeps the default 4 round
/// trips. A SHORT low-loss duct (a 17 cm vocal tract) rings far past
/// 4 round trips, and the wrapped/truncated FIR shifts its upper
/// resonances measurably (F2 of the /a/ tract measured 19% low at the
/// default) — densify for chart-fidelity consumers.
///
/// # Errors
/// As [`characteristic_line`].
#[allow(clippy::too_many_arguments)]
pub fn characteristic_line_dense(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
    n_fft_override: Option<usize>,
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
    let n_fft = n_fft_override.map_or_else(
        || {
            ((4.0 * geo_delay).ceil() as usize)
                .next_power_of_two()
                .clamp(256, 4096)
        },
        |v| v.next_power_of_two().clamp(256, 8192),
    );
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
        let response =
            input_impedance_wall(physics, gas, omega, LossModel::Bessel, termination, wall)
                .map_err(DrivingPointError::Duct)?;
        let rac = reflectance(C64::new(response.impedance.re, response.impedance.im), zc);
        // Acoustic e^{-iωt} → DFT e^{+iωt} so the IR is causal.
        buf[k] = FftC64::new(rac.re, -rac.im);
        if k != n_fft / 2 {
            buf[n_fft - k] = FftC64::new(rac.re, rac.im);
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.inverse(&mut buf, &mut scratch);
    let ir: Vec<f64> = buf.iter().map(|c| c.re).collect();
    DelayedFilter::from_impulse_response(dt, ir).map_err(DrivingPointError::Discrete)
}

/// Causal driving-point impedance as the same FIR port object.
///
/// Samples `Z(ω)` (closed) or the mouth transfer `p/u` (open) on the
/// DFT grid of the run and inverse-transforms it. Isolated linear
/// blow is then streaming convolution against that IR — the same
/// `DelayedFilter` a reed or a body uses for `R(ω)`. A tone hole
/// changes `Z` and therefore the impulse. There is no one-pole
/// fallback and no block-wise `U(ω) Z(ω)` special case.
///
/// `Z(0)` is taken as zero so a DC volume-velocity does not mint a
/// static pressure. TMM `Z` is stored on the DFT grid without a
/// further conjugate: the historical block `U(ω) Z(ω)` path and a
/// vented tone-hole period both require that convention.
///
/// # Errors
/// Empty geometry, TMM, or DFT.
pub fn impedance_line(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<&WallPin>,
) -> Result<DelayedFilter, DrivingPointError> {
    if physics.segments.is_empty() {
        return Err(DrivingPointError::Invalid {
            what: "line needs positive axial length",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let n_fft = n.next_power_of_two().clamp(256, 8192);
    let fft = Fft::new(n_fft);
    let mut buf = vec![FftC64::new(0.0, 0.0); n_fft];
    for k in 1..=n_fft / 2 {
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let response =
            input_impedance_wall(physics, gas, omega, LossModel::Bessel, termination, wall)
                .map_err(DrivingPointError::Duct)?;
        let z = match termination {
            Termination::Closed => response.impedance,
            Termination::IdealOpen | Termination::UnflangedOpen | Termination::FlangedOpen => {
                response.p_mouth_over_u_in
            }
        };
        // Same convention as the historical block `U(ω) Z(ω)` path:
        // TMM `Z` is already the DFT-side sample; conjugating it
        // time-reverses a vented IR and lengthens the measured period.
        buf[k] = FftC64::new(z.re, z.im);
        if k != n_fft / 2 {
            buf[n_fft - k] = FftC64::new(z.re, -z.im);
        }
    }
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.inverse(&mut buf, &mut scratch);
    let ir: Vec<f64> = buf.iter().map(|c| c.re).collect();
    DelayedFilter::from_impulse_response(dt, ir).map_err(DrivingPointError::Discrete)
}
