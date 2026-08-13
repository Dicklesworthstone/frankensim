//! Realize an [`fs_scenario::AcousticAssembly`] into observer pressure.
//!
//! A guitar or clarinet is a description: Kirchhoff–Carrier strings,
//! Rayleigh damping, a radiating plate, a beating reed on a viscothermal
//! TMM bore, optional tone holes, and a regularized bow. There is no
//! instrument crate.

use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel,
};
use crate::pcm_wav::{WavError, encode_pcm16_wav};
use crate::reed_bore::realize_reed_bore;
use crate::stribeck_friction::StribeckFriction;
use fs_duct::{
    Duct, DuctError, HoleState, LossModel, MAX_RADIATION_KA, Segment, Termination, input_impedance,
};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::RayleighDamping;
use fs_math::c64::C64;
use fs_math::det;
use fs_nlmodal::{KcStringParams, assemble, kirchhoff_carrier_string, prestressed_beam_omega};
use fs_phs::step;
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, BowStroke, Pluck, PrestressedString, RadiatingPlate,
    RayleighParams, ViscothermalDuct, VolumeVelocityPulse, WaveguideEnd,
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
    /// Kirchhoff–Carrier / pHS refusal.
    Nonlinear(String),
    /// Reed–bore loop refusal.
    Reed {
        /// Which check failed.
        what: &'static str,
    },
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
            Self::Nonlinear(e) => write!(f, "FS-COUPLE-ASSEMBLY-NL: {e}"),
            Self::Reed { what } => write!(f, "FS-COUPLE-REED: {what}"),
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
        if assembly.pluck.is_none() && assembly.bow.is_none() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a string member requires a pluck or a bow",
            });
        }
        let hist = realize_string(
            string,
            assembly.pluck,
            assembly.bow,
            assembly.soundboard,
            &assembly.body_modes,
            &gas,
            assembly.listener.distance_m,
            assembly.sample_rate_hz,
            n,
        )?;
        add_in_place(&mut pressure_pa, &hist);
    }
    if let Some(duct) = assembly.duct.as_ref() {
        let hist = if let Some(reed) = assembly.reed {
            realize_reed_on_duct(duct, reed, &gas, assembly.sample_rate_hz, n)?
        } else {
            let blow = assembly
                .blow
                .ok_or(AcousticRealizeError::InvalidDescription {
                    what: "a duct member requires a volume-velocity pulse or a reed",
                })?;
            realize_blown_duct(duct, blow, &gas, assembly.sample_rate_hz, n)?
        };
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

/// Prestressed Euler–Bernoulli frequency of sine mode `k` (1-based).
#[must_use]
pub fn string_mode_omega(string: PrestressedString, k: usize) -> f64 {
    prestressed_beam_omega(
        string.length_m,
        string.tension_n,
        string.lin_density_kg_m,
        string.bending_stiffness_n_m2,
        k,
    )
}

fn mode_zeta(
    string: PrestressedString,
    omega: f64,
    gas: &GasState,
) -> Result<f64, AcousticRealizeError> {
    if let Some(RayleighParams {
        alpha_per_s,
        beta_s,
    }) = string.rayleigh
    {
        let r = RayleighDamping::new(alpha_per_s, beta_s).map_err(|_| {
            AcousticRealizeError::InvalidDescription {
                what: "Rayleigh coefficients must be finite and non-negative",
            }
        })?;
        return Ok(r.zeta_at(omega));
    }
    // Stokes air drag on a cylinder of radius ≈ width/2, plus the
    // caller-authored internal floor, plus a stiffness-proportional
    // term from bending (Valette–Cuesta / Chaigne class).
    let radius = 0.5 * string.width_m;
    let stokes = core::f64::consts::PI
        * radius
        * det::sqrt(2.0 * gas.dynamic_viscosity * gas.density / omega.max(1.0))
        / (2.0 * string.lin_density_kg_m * det::sqrt(omega.max(1.0)));
    let bend = if string.bending_stiffness_n_m2 > 0.0 {
        2.0e-7 * omega
    } else {
        0.0
    };
    Ok(string.damping_ratio + stokes + bend)
}

#[allow(clippy::too_many_arguments)]
fn realize_string(
    string: PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
    soundboard: Option<RadiatingPlate>,
    body_modes: &[RadiatingPlate],
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    validate_string(string, pluck, bow)?;
    let mut hist = if string.axial_stiffness_n > 0.0 {
        realize_kc_string(
            string,
            pluck,
            bow,
            soundboard,
            body_modes,
            gas,
            listener_m,
            sample_rate_hz,
            n,
            1.0,
        )?
    } else {
        realize_linear_string(
            string,
            pluck,
            bow,
            soundboard,
            body_modes,
            gas,
            listener_m,
            sample_rate_hz,
            n,
            1.0,
        )?
    };
    if string.polarization_detune > 0.0 {
        let mut twin = string;
        twin.polarization_detune = 0.0;
        twin.tension_n *= (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
        let other = if string.axial_stiffness_n > 0.0 {
            realize_kc_string(
                twin,
                pluck,
                bow,
                None,
                &[],
                gas,
                listener_m,
                sample_rate_hz,
                n,
                0.85,
            )?
        } else {
            realize_linear_string(
                twin,
                pluck,
                bow,
                None,
                &[],
                gas,
                listener_m,
                sample_rate_hz,
                n,
                0.85,
            )?
        };
        add_in_place(&mut hist, &other);
    }
    Ok(hist)
}

fn validate_string(
    string: PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
) -> Result<(), AcousticRealizeError> {
    if !(string.length_m > 0.0
        && string.tension_n > 0.0
        && string.lin_density_kg_m > 0.0
        && string.axial_stiffness_n >= 0.0
        && string.width_m > 0.0
        && string.n_modes > 0
        && string.damping_ratio >= 0.0
        && string.bending_stiffness_n_m2 >= 0.0
        && string.polarization_detune >= 0.0
        && string.length_m.is_finite()
        && string.tension_n.is_finite()
        && string.lin_density_kg_m.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "string parameters must be physical and finite",
        });
    }
    if let Some(pluck) = pluck
        && !(pluck.station_frac > 0.0 && pluck.station_frac < 1.0 && pluck.height_m.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "pluck station must lie in (0, 1) with a finite height",
        });
    }
    if let Some(bow) = bow
        && !(bow.station_frac > 0.0
            && bow.station_frac < 1.0
            && bow.normal_force_n > 0.0
            && bow.stribeck_m_s > 0.0
            && bow.mu_static >= bow.mu_dynamic
            && bow.mu_dynamic >= 0.0)
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "bow station must lie in (0, 1) with physical friction",
        });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn realize_linear_string(
    string: PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
    soundboard: Option<RadiatingPlate>,
    body_modes: &[RadiatingPlate],
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    radiation_scale: f64,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let pi = core::f64::consts::PI;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let mut modes = Vec::with_capacity(string.n_modes);
    let mut states = Vec::with_capacity(string.n_modes);
    let mut phi_bow = vec![0.0; string.n_modes];
    for k in 1..=string.n_modes {
        let kf = k as f64;
        let omega = string_mode_omega(string, k);
        let q_phys = pluck.map_or(0.0, |p| triangular_pluck_modal(p, k));
        let q = mass_scale * q_phys;
        let monopole_area = if k % 2 == 1 {
            string.width_m * 2.0 * string.length_m / (kf * pi)
        } else {
            0.0
        };
        let transfer = C64::new(
            0.0,
            radiation_scale * omega * gas.density * monopole_area
                / (4.0 * pi * listener_m * mass_scale),
        );
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: mode_zeta(string, omega, gas)?,
            pressure_per_modal_velocity: transfer,
        });
        states.push(ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: 0.0,
        });
        if let Some(bow) = bow {
            phi_bow[k - 1] = det::sin(kf * pi * bow.station_frac) / mass_scale;
        }
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
    let mut plates = plate_bank(soundboard, body_modes)?;
    let mut out = Vec::with_capacity(n);
    let dt = 1.0 / f64::from(sample_rate_hz);
    for _ in 0..n {
        let mut force = vec![0.0; string.n_modes];
        if let Some(bow) = bow {
            let v_string: f64 = model
                .states()
                .iter()
                .zip(phi_bow.iter())
                .map(|(s, phi)| s.velocity_m_sqrt_kg_per_s * phi)
                .sum();
            let f_bow = bow_force(bow, v_string);
            for (f, phi) in force.iter_mut().zip(phi_bow.iter()) {
                *f = f_bow * *phi;
            }
        }
        let frame = model.step(&force).map_err(AcousticRealizeError::Modal)?;
        let mut p = frame.observer_pressure_pa;
        if !plates.is_empty() {
            let q_phys: Vec<f64> = model
                .states()
                .iter()
                .map(|s| s.displacement_m_sqrt_kg / mass_scale)
                .collect();
            let fb = bridge_force(string, &q_phys);
            for plate in &mut plates {
                p += plate.drive_and_radiate(fb, dt, gas.density, listener_m);
            }
        }
        out.push(p);
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn realize_kc_string(
    string: PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
    soundboard: Option<RadiatingPlate>,
    body_modes: &[RadiatingPlate],
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    radiation_scale: f64,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut storage = kirchhoff_carrier_string(
        &KcStringParams {
            length: string.length_m,
            tension: string.tension_n,
            lin_density: string.lin_density_kg_m,
            ea: string.axial_stiffness_n,
        },
        string.n_modes,
    )
    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    if string.bending_stiffness_n_m2 > 0.0 {
        for (k, w) in storage.omegas.iter_mut().enumerate() {
            *w = string_mode_omega(string, k + 1);
        }
    }
    let zetas: Result<Vec<f64>, _> = storage
        .omegas
        .iter()
        .map(|&w| mode_zeta(string, w, gas))
        .collect();
    let zetas = zetas?;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let pi = core::f64::consts::PI;
    let mut strike = vec![0.0; string.n_modes];
    if let Some(bow) = bow {
        for k in 1..=string.n_modes {
            strike[k - 1] = det::sin(k as f64 * pi * bow.station_frac) / mass_scale;
        }
    }
    let sys = assemble(storage, &zetas, &strike)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let mut x = vec![0.0; 2 * string.n_modes];
    if let Some(pluck) = pluck {
        for k in 1..=string.n_modes {
            x[2 * (k - 1)] = mass_scale * triangular_pluck_modal(pluck, k);
        }
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut plates = plate_bank(soundboard, body_modes)?;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let u = if let Some(bow) = bow {
            let v_string: f64 = (0..string.n_modes)
                .map(|k| {
                    let phi = det::sin((k + 1) as f64 * pi * bow.station_frac) / mass_scale;
                    x[2 * k + 1] * phi
                })
                .sum();
            vec![bow_force(bow, v_string)]
        } else {
            vec![0.0]
        };
        let rec =
            step(&sys, &x, &u, dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        x = rec.x;
        let mut p = 0.0;
        let mut q_phys = vec![0.0; string.n_modes];
        for k in 0..string.n_modes {
            let kf = (k + 1) as f64;
            let q = x[2 * k];
            q_phys[k] = q / mass_scale;
            if k % 2 == 0 {
                let omega = string_mode_omega(string, k + 1);
                let area = string.width_m * 2.0 * string.length_m / (kf * pi);
                let h_im = radiation_scale * omega * gas.density * area
                    / (4.0 * pi * listener_m * mass_scale);
                p += h_im * omega * q;
            }
        }
        if !plates.is_empty() {
            let fb = bridge_force(string, &q_phys);
            for plate in &mut plates {
                p += plate.drive_and_radiate(fb, dt, gas.density, listener_m);
            }
        }
        out.push(p);
    }
    Ok(out)
}

fn triangular_pluck_modal(pluck: Pluck, k: usize) -> f64 {
    let xi = pluck.station_frac;
    let kf = k as f64;
    let pi = core::f64::consts::PI;
    2.0 * pluck.height_m * det::sin(kf * pi * xi) / (kf * kf * pi * pi * xi * (1.0 - xi))
}

fn bow_force(bow: BowStroke, v_string: f64) -> f64 {
    StribeckFriction {
        mu_static: bow.mu_static,
        mu_dynamic: bow.mu_dynamic,
        stiction_m_s: bow.stribeck_m_s,
    }
    .traction(bow.velocity_m_s - v_string, bow.normal_force_n)
}

fn plate_bank(
    soundboard: Option<RadiatingPlate>,
    extra: &[RadiatingPlate],
) -> Result<Vec<PlateState>, AcousticRealizeError> {
    let mut out = Vec::with_capacity(extra.len() + 1);
    if let Some(spec) = soundboard {
        out.push(PlateState::new(spec)?);
    }
    for &spec in extra {
        out.push(PlateState::new(spec)?);
    }
    Ok(out)
}

fn bridge_force(string: PrestressedString, q_phys: &[f64]) -> f64 {
    let pi = core::f64::consts::PI;
    let mut slope = 0.0;
    let mut stretch = 0.0;
    for (k, &q) in q_phys.iter().enumerate() {
        let kn = (k + 1) as f64 * pi / string.length_m;
        slope += q * kn;
        stretch += kn * kn * q * q;
    }
    let integral = 0.5 * string.length_m * stretch;
    let t_eff = string.tension_n + string.axial_stiffness_n / (2.0 * string.length_m) * integral;
    t_eff * slope
}

struct PlateState {
    spec: RadiatingPlate,
    omega: f64,
    y: f64,
    v: f64,
}

impl PlateState {
    fn new(spec: RadiatingPlate) -> Result<Self, AcousticRealizeError> {
        if !(spec.area_m2 > 0.0
            && spec.mass_kg > 0.0
            && spec.frequency_hz > 0.0
            && spec.damping_ratio >= 0.0)
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "soundboard parameters must be physical",
            });
        }
        Ok(Self {
            spec,
            omega: core::f64::consts::TAU * spec.frequency_hz,
            y: 0.0,
            v: 0.0,
        })
    }

    fn drive_and_radiate(&mut self, force_n: f64, dt: f64, rho: f64, listener_m: f64) -> f64 {
        let acc = force_n / self.spec.mass_kg
            - 2.0 * self.spec.damping_ratio * self.omega * self.v
            - self.omega * self.omega * self.y;
        self.v += dt * acc;
        self.y += dt * self.v;
        rho * self.spec.area_m2 * acc / (4.0 * core::f64::consts::PI * listener_m)
    }
}

fn realize_reed_on_duct(
    duct: &ViscothermalDuct,
    reed: BeatingReed,
    gas: &GasState,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let physics = physics_duct(duct)?;
    let termination = match duct.termination {
        WaveguideEnd::Closed => Termination::Closed,
        WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
    };
    refuse_open_nyquist(duct, gas, sample_rate_hz)?;
    realize_reed_bore(&physics, gas, reed, termination, sample_rate_hz, n)
}

fn realize_blown_duct(
    duct: &ViscothermalDuct,
    blow: VolumeVelocityPulse,
    gas: &GasState,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    if !(blow.peak_m3_s.is_finite() && blow.duration_s > 0.0 && blow.duration_s.is_finite()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "blow pulse must have finite peak and positive duration",
        });
    }
    let physics = physics_duct(duct)?;
    let termination = match duct.termination {
        WaveguideEnd::Closed => Termination::Closed,
        WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
    };
    refuse_open_nyquist(duct, gas, sample_rate_hz)?;
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
        let response = input_impedance(&physics, gas, omega, LossModel::AllRegime, termination)
            .map_err(AcousticRealizeError::Duct)?;
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

fn physics_duct(duct: &ViscothermalDuct) -> Result<Duct, AcousticRealizeError> {
    if duct.segments.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "duct has no segments",
        });
    }
    let mut segments = Vec::new();
    for (i, s) in duct.segments.iter().enumerate() {
        segments.push(Segment::Cylinder {
            radius: s.radius_m,
            length: s.length_m,
        });
        for hole in &duct.tone_holes {
            if hole.after_segment == i {
                segments.push(Segment::ToneHole {
                    hole_radius: hole.radius_m,
                    chimney_height: hole.chimney_m,
                    bore_radius: s.radius_m,
                    state: if hole.open {
                        HoleState::Open
                    } else {
                        HoleState::Closed
                    },
                });
            }
        }
    }
    Ok(Duct { segments })
}

fn refuse_open_nyquist(
    duct: &ViscothermalDuct,
    gas: &GasState,
    sample_rate_hz: u32,
) -> Result<(), AcousticRealizeError> {
    if !matches!(duct.termination, WaveguideEnd::UnflangedOpen) {
        return Ok(());
    }
    let mouth_r = duct
        .segments
        .last()
        .expect("non-empty checked above")
        .radius_m;
    let ka_nyquist = core::f64::consts::PI * f64::from(sample_rate_hz) * mouth_r / gas.sound_speed;
    if ka_nyquist > MAX_RADIATION_KA {
        return Err(AcousticRealizeError::Duct(DuctError::RadiationKaTooLarge {
            ka: ka_nyquist,
        }));
    }
    Ok(())
}

fn add_in_place(acc: &mut [f64], addend: &[f64]) {
    for (a, b) in acc.iter_mut().zip(addend) {
        *a += *b;
    }
}
