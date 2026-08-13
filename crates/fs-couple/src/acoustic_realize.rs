//! Realize an [`fs_scenario::AcousticAssembly`] into observer pressure.
//!
//! Strings become mass-normalized modal oscillators observed through a
//! compact first-principles radiation transfer that uses the assembly's
//! [`fs_material::gas::GasState`]. Ducts become an inlet-pressure history
//! `p_in(t) = IFFT[Z_in(ω) U(ω)]` from the viscothermal TMM. There is no
//! instrument crate: a guitar or clarinet is a description.

use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel,
};
use crate::pcm_wav::{WavError, encode_pcm16_wav};
use fs_duct::{
    Duct, DuctError, LossModel, MAX_RADIATION_KA, Segment, Termination, input_impedance,
};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_math::det;
use fs_scenario::{
    AcousticAssembly, AmbientGas, Pluck, PrestressedString, ViscothermalDuct, VolumeVelocityPulse,
    WaveguideEnd,
};

/// Typed realization refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum AcousticRealizeError {
    /// Description failed a physical domain check.
    InvalidDescription {
        /// Which field.
        what: &'static str,
    },
    /// Gas-state derivation failed.
    Ambient(String),
    /// Modal time-stepper refusal.
    Modal(ModalAcousticTimeError),
    /// Waveguide TMM refusal.
    Duct(fs_duct::DuctError),
    /// WAV encode refusal.
    Wav(WavError),
}

impl core::fmt::Display for AcousticRealizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDescription { what } => {
                write!(f, "FS-COUPLE-ASSEMBLY: {what}")
            }
            Self::Ambient(what) => write!(f, "FS-COUPLE-ASSEMBLY-GAS: {what}"),
            Self::Modal(e) => write!(f, "{e}"),
            Self::Duct(e) => write!(f, "{e}"),
            Self::Wav(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for AcousticRealizeError {}

/// Realized observer pressure and the gas state that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct RealizedAssembly {
    /// Sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Observer pressure [Pa], one sample per period.
    pub pressure_pa: Vec<f64>,
    /// Derived ambient state (rho, c, mu, …).
    pub gas: GasState,
}

/// Realize every described member and sum their observer pressures.
///
/// Linear acoustics: string radiation and duct inlet pressure occupy the
/// same timeline and add. Structure–bore coupling is a recorded no-claim.
///
/// # Errors
/// Empty assemblies, invalid events, gas/domain refusals, stepper or TMM
/// refusals.
pub fn realize_assembly(
    assembly: &AcousticAssembly,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    if assembly.string.is_none() && assembly.duct.is_none() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "assembly has neither a string nor a duct",
        });
    }
    if assembly.sample_rate_hz == 0
        || !(assembly.duration_s > 0.0 && assembly.duration_s.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "sample rate and duration must be positive and finite",
        });
    }
    if !(assembly.listener.distance_m > 0.0 && assembly.listener.distance_m.is_finite()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "listener distance must be positive and finite",
        });
    }
    let gas = gas_state(assembly.ambient)?;
    let n = sample_count(assembly.sample_rate_hz, assembly.duration_s)?;
    let mut pressure_pa = vec![0.0; n];
    if let Some(string) = assembly.string {
        let pluck = assembly
            .pluck
            .ok_or(AcousticRealizeError::InvalidDescription {
                what: "a string member requires a pluck event",
            })?;
        let hist = realize_plucked_string(
            string,
            pluck,
            &gas,
            assembly.listener.distance_m,
            assembly.sample_rate_hz,
            n,
        )?;
        add_in_place(&mut pressure_pa, &hist);
    }
    if let Some(duct) = assembly.duct.as_ref() {
        let blow = assembly
            .blow
            .ok_or(AcousticRealizeError::InvalidDescription {
                what: "a duct member requires a volume-velocity pulse",
            })?;
        let hist = realize_blown_duct(duct, blow, &gas, assembly.sample_rate_hz, n)?;
        add_in_place(&mut pressure_pa, &hist);
    }
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa,
        gas,
    })
}

/// Encode the realization as PCM16 WAV.
///
/// # Errors
/// WAV domain refusals.
pub fn assembly_wav(
    realized: &RealizedAssembly,
    full_scale_pa: f64,
) -> Result<(Vec<u8>, usize), AcousticRealizeError> {
    encode_pcm16_wav(
        &realized.pressure_pa,
        realized.sample_rate_hz,
        full_scale_pa,
    )
    .map_err(AcousticRealizeError::Wav)
}

fn gas_state(ambient: AmbientGas) -> Result<GasState, AcousticRealizeError> {
    GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        ambient.temperature_k,
        ambient.pressure_pa,
    )
    .map_err(|e| AcousticRealizeError::Ambient(e.to_string()))
}

fn sample_count(rate_hz: u32, duration_s: f64) -> Result<usize, AcousticRealizeError> {
    let n = f64::from(rate_hz) * duration_s;
    if !(n.is_finite() && n >= 2.0) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "realized history needs at least two samples",
        });
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    Ok(n.floor() as usize)
}

fn realize_plucked_string(
    string: PrestressedString,
    pluck: Pluck,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    if !(string.length_m > 0.0
        && string.tension_n > 0.0
        && string.lin_density_kg_m > 0.0
        && string.axial_stiffness_n >= 0.0
        && string.width_m > 0.0
        && string.n_modes > 0
        && string.damping_ratio >= 0.0
        && string.length_m.is_finite()
        && string.tension_n.is_finite()
        && string.lin_density_kg_m.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "string parameters must be physical and finite",
        });
    }
    if !(pluck.station_frac > 0.0 && pluck.station_frac < 1.0 && pluck.height_m.is_finite()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "pluck station must lie in (0, 1) with a finite height",
        });
    }
    let wave = det::sqrt(string.tension_n / string.lin_density_kg_m);
    let pi = core::f64::consts::PI;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let mut modes = Vec::with_capacity(string.n_modes);
    let mut states = Vec::with_capacity(string.n_modes);
    for k in 1..=string.n_modes {
        let kf = k as f64;
        let omega = kf * pi * wave / string.length_m;
        let q_phys = triangular_pluck_modal(pluck, k);
        let q = mass_scale * q_phys;
        let monopole_area = if k % 2 == 1 {
            string.width_m * 2.0 * string.length_m / (kf * pi)
        } else {
            0.0
        };
        let transfer = C64::new(
            0.0,
            omega * gas.density * monopole_area / (4.0 * pi * listener_m * mass_scale),
        );
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: string.damping_ratio,
            pressure_per_modal_velocity: transfer,
        });
        states.push(ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: 0.0,
        });
    }
    let mut model = ModalAcousticTimeModel::try_new(
        sample_rate_hz,
        modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .map_err(AcousticRealizeError::Modal)?;
    model
        .restore_states(&states)
        .map_err(AcousticRealizeError::Modal)?;
    let zeros = vec![0.0; string.n_modes];
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let frame = model.step(&zeros).map_err(AcousticRealizeError::Modal)?;
        out.push(frame.observer_pressure_pa);
    }
    Ok(out)
}

fn triangular_pluck_modal(pluck: Pluck, k: usize) -> f64 {
    let xi = pluck.station_frac;
    let kf = k as f64;
    let pi = core::f64::consts::PI;
    2.0 * pluck.height_m * det::sin(kf * pi * xi) / (kf * kf * pi * pi * xi * (1.0 - xi))
}

fn realize_blown_duct(
    duct: &ViscothermalDuct,
    blow: VolumeVelocityPulse,
    gas: &GasState,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    if duct.segments.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "duct has no segments",
        });
    }
    if !(blow.peak_m3_s.is_finite() && blow.duration_s > 0.0 && blow.duration_s.is_finite()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "blow pulse must have finite peak and positive duration",
        });
    }
    let physics = Duct {
        segments: duct
            .segments
            .iter()
            .map(|s| Segment::Cylinder {
                radius: s.radius_m,
                length: s.length_m,
            })
            .collect(),
    };
    let termination = match duct.termination {
        WaveguideEnd::Closed => Termination::Closed,
        WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
    };
    if matches!(duct.termination, WaveguideEnd::UnflangedOpen) {
        let mouth_r = duct
            .segments
            .last()
            .expect("non-empty checked above")
            .radius_m;
        let ka_nyquist =
            core::f64::consts::PI * f64::from(sample_rate_hz) * mouth_r / gas.sound_speed;
        if ka_nyquist > MAX_RADIATION_KA {
            return Err(AcousticRealizeError::Duct(DuctError::RadiationKaTooLarge {
                ka: ka_nyquist,
            }));
        }
    }
    let n_fft = n.next_power_of_two();
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut drive = vec![0.0; n_fft];
    for (i, sample) in drive.iter_mut().enumerate() {
        let t = i as f64 * dt;
        if t < blow.duration_s {
            *sample = blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s);
        }
    }
    let fft = Fft::new(n_fft);
    let mut buf: Vec<FftC64> = drive.iter().map(|&u| FftC64::new(u, 0.0)).collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.forward(&mut buf, &mut scratch);
    for (k, bin) in buf.iter_mut().enumerate().take(n_fft / 2 + 1) {
        if k == 0 {
            *bin = FftC64::new(0.0, 0.0);
            continue;
        }
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let response = duct_response(&physics, gas, omega, termination)?;
        // Closed: observe inlet pressure (reed/blow plane). Open: observe
        // mouth pressure from the ABCD mouth-flow ratio times Z_L.
        let h = match duct.termination {
            WaveguideEnd::Closed => response.impedance,
            WaveguideEnd::UnflangedOpen => response.p_mouth_over_u_in,
        };
        *bin = FftC64::new(bin.re * h.re - bin.im * h.im, bin.re * h.im + bin.im * h.re);
    }
    for k in 1..n_fft / 2 {
        let conj = buf[k];
        buf[n_fft - k] = FftC64::new(conj.re, -conj.im);
    }
    fft.inverse(&mut buf, &mut scratch);
    Ok(buf[..n].iter().map(|c| c.re).collect())
}

fn duct_response(
    physics: &Duct,
    gas: &GasState,
    omega: f64,
    termination: Termination,
) -> Result<fs_duct::DuctResponse, AcousticRealizeError> {
    match input_impedance(physics, gas, omega, LossModel::WideTube, termination) {
        Ok(response) => Ok(response),
        Err(DuctError::TooNarrow { .. }) => {
            // Wide-tube ZK is invalid at this shear number; lossless Z is
            // still defined. Losses at that bin are a recorded no-claim.
            input_impedance(physics, gas, omega, LossModel::Lossless, termination)
                .map_err(AcousticRealizeError::Duct)
        }
        Err(err) => Err(AcousticRealizeError::Duct(err)),
    }
}

fn add_in_place(acc: &mut [f64], addend: &[f64]) {
    for (a, b) in acc.iter_mut().zip(addend) {
        *a += *b;
    }
}
