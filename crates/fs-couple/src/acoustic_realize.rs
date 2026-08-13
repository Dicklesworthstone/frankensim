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
use crate::thin_plate::{PlateBank, VkBody, certified_radiators};
use crate::unilateral_contact::modal_contact_forces;
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
    AcousticAssembly, AmbientGas, BeatingReed, BowStroke, ContactTexture, Pluck, PrestressedString,
    RadiatingPlate, RayleighParams, ThinPlate, UnilateralObstacle, ViscothermalDuct,
    VolumeVelocityPulse, WaveguideEnd,
};
use fs_tribo::{
    InputAuthority, InterfaceMedium, InterfaceSystemRef,
    surface_excitation::{
        PeriodicHarmonicSurface, SelfAffinePeriodicProfileSpectrum, SurfaceTraceMotion,
        UniformSurfaceTrace, evaluate_point_surface_pair,
    },
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
    let mut bodies = plate_bank(assembly.soundboard, &assembly.body_modes, assembly.plate)?;
    bodies.attach_radiation_loads(&gas, assembly.sample_rate_hz);
    if assembly.string.is_some() && assembly.duct.is_some() {
        return realize_coupled(assembly, &gas, &mut bodies, n);
    }
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
            assembly.contact_texture,
            &mut bodies,
            &assembly.obstacles,
            &gas,
            assembly.listener.distance_m,
            assembly.sample_rate_hz,
            n,
        )?;
        add_in_place(&mut pressure_pa, &hist);
    }
    if let Some(duct) = assembly.duct.as_ref() {
        let hist = if let Some(reed) = assembly.reed {
            realize_reed_on_duct(
                duct,
                reed,
                &mut bodies,
                &gas,
                assembly.listener.distance_m,
                assembly.sample_rate_hz,
                n,
            )?
        } else {
            let blow = assembly
                .blow
                .ok_or(AcousticRealizeError::InvalidDescription {
                    what: "a duct member requires a volume-velocity pulse or a reed",
                })?;
            realize_blown_duct(
                duct,
                blow,
                &mut bodies,
                &gas,
                assembly.listener.distance_m,
                assembly.sample_rate_hz,
                n,
            )?
        };
        add_in_place(&mut pressure_pa, &hist);
    }
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa,
        gas,
    })
}

/// One shared clock: string, plate, and duct exchange force and flow
/// every sample. Sequential string-then-duct is not this function.
fn realize_coupled(
    assembly: &AcousticAssembly,
    gas: &GasState,
    plates: &mut PlateBank,
    n: usize,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    use crate::driving_point::characteristic_line;
    use crate::reed_bore::solve_reed_wave;
    use crate::traveling_wave_line::TravelingWaveLine;
    let string = assembly.string.expect("checked");
    let duct = assembly.duct.as_ref().expect("checked");
    if assembly.pluck.is_none() && assembly.bow.is_none() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "a string member requires a pluck or a bow",
        });
    }
    if assembly.reed.is_none() && assembly.blow.is_none() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "a coupled duct requires a reed or a volume-velocity pulse",
        });
    }
    validate_string(string, assembly.pluck, assembly.bow)?;
    let physics = physics_duct(duct)?;
    let termination = match duct.termination {
        WaveguideEnd::Closed => Termination::Closed,
        WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
    };
    refuse_open_nyquist(duct, gas, assembly.sample_rate_hz)?;
    let inlet_r = physics
        .segments
        .first()
        .ok_or(AcousticRealizeError::InvalidDescription {
            what: "duct has no segments",
        })?
        .outlet_radius();
    let area = core::f64::consts::PI * inlet_r * inlet_r;
    let zc = gas.density * gas.sound_speed / area;
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    let listener_m = assembly.listener.distance_m;
    let mut fitted =
        characteristic_line(&physics, gas, termination, assembly.sample_rate_hz, n, zc).ok();
    let mut delay = if fitted.is_none() {
        Some(
            TravelingWaveLine::from_duct(
                &physics,
                gas,
                termination,
                assembly.sample_rate_hz,
                n,
                zc,
            )
            .map_err(|e| match e {
                crate::traveling_wave_line::TravelingWaveError::Invalid { what } => {
                    AcousticRealizeError::Reed { what }
                }
                crate::traveling_wave_line::TravelingWaveError::Duct(d) => {
                    AcousticRealizeError::Duct(d)
                }
            })?,
        )
    } else {
        None
    };
    let mut texture = TextureDrive::try_new(assembly.contact_texture)?;
    if string.axial_stiffness_n > 0.0 {
        return realize_coupled_kc(
            assembly,
            gas,
            plates,
            n,
            fitted,
            delay,
            zc,
            area,
            &mut texture,
        );
    }
    let pi = core::f64::consts::PI;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let mut modes = Vec::new();
    let mut states = Vec::new();
    let mut phi_bow = vec![0.0; string.n_modes];
    for k in 1..=string.n_modes {
        let kf = k as f64;
        let omega = string_mode_omega(string, k);
        let q_phys = assembly.pluck.map_or(0.0, |p| triangular_pluck_modal(p, k));
        let monopole_area = if k % 2 == 1 {
            string.width_m * 2.0 * string.length_m / (kf * pi)
        } else {
            0.0
        };
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: mode_zeta(string, omega, gas)?,
            pressure_per_modal_velocity: C64::new(
                0.0,
                omega * gas.density * monopole_area / (4.0 * pi * listener_m * mass_scale),
            ),
        });
        states.push(ModalAcousticState {
            displacement_m_sqrt_kg: mass_scale * q_phys,
            velocity_m_sqrt_kg_per_s: 0.0,
        });
        if let Some(bow) = assembly.bow {
            phi_bow[k - 1] = det::sin(kf * pi * bow.station_frac) / mass_scale;
        }
    }
    let mut model = ModalAcousticTimeModel::try_new(
        assembly.sample_rate_hz,
        modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .map_err(AcousticRealizeError::Modal)?;
    model
        .restore_states(&states)
        .map_err(AcousticRealizeError::Modal)?;
    let mut p_plus_prev = 5.0;
    let mut out = vec![0.0; n];
    for i in 0..n {
        let mut force = vec![0.0; string.n_modes];
        if let Some(bow) = assembly.bow {
            let v_string: f64 = model
                .states()
                .iter()
                .zip(phi_bow.iter())
                .map(|(s, phi)| s.velocity_m_sqrt_kg_per_s * phi)
                .sum();
            let f_bow = bow_force(bow, v_string, texture.delta_n(bow.velocity_m_s, dt));
            for (f, phi) in force.iter_mut().zip(phi_bow.iter()) {
                *f = f_bow * *phi;
            }
        }
        if !assembly.obstacles.is_empty() {
            let x: Vec<f64> = model
                .states()
                .iter()
                .flat_map(|s| [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s])
                .collect();
            let extra = modal_contact_forces(string, &assembly.obstacles, &x)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            for (f, c) in force.iter_mut().zip(extra) {
                *f += c;
            }
        }
        let frame = model.step(&force).map_err(AcousticRealizeError::Modal)?;
        let q_phys: Vec<f64> = model
            .states()
            .iter()
            .map(|s| s.displacement_m_sqrt_kg / mass_scale)
            .collect();
        let fb = bridge_force(string, &q_phys);
        let p_minus = if let Some(line) = fitted.as_ref() {
            line.incoming()
        } else {
            delay.as_ref().expect("line").incoming()
        };
        let u_body = plates.volume_velocity();
        let t = i as f64 * dt;
        let p_plus = if let Some(reed) = assembly.reed {
            let p_m = if reed.attack_s <= 0.0 {
                reed.blowing_pressure_pa
            } else if t >= reed.attack_s {
                reed.blowing_pressure_pa
            } else {
                let x = t / reed.attack_s;
                reed.blowing_pressure_pa * 0.5 * (1.0 - det::cos(core::f64::consts::PI * x))
            };
            solve_reed_wave(
                reed,
                gas.density,
                zc,
                0.0,
                p_minus,
                p_m,
                p_plus_prev,
                u_body,
            )?
        } else {
            let blow = assembly.blow.expect("checked");
            let u_blow = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
            p_minus + zc * (u_blow + u_body)
        };
        p_plus_prev = p_plus;
        let p_minus_now = if let Some(line) = fitted.as_mut() {
            line.push(p_plus).map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?
        } else {
            delay.as_mut().expect("line").push(p_plus)
        };
        let p_bore = p_plus + p_minus_now;
        let mut p = frame.observer_pressure_pa + p_bore;
        p += plates.drive_and_radiate(fb + p_bore * area, dt, gas.density, listener_m)?;
        out[i] = p;
    }
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa: out,
        gas: gas.clone(),
    })
}

/// Shared clock with a Kirchhoff–Carrier string (EA > 0).
#[allow(clippy::too_many_arguments)]
fn realize_coupled_kc(
    assembly: &AcousticAssembly,
    gas: &GasState,
    plates: &mut PlateBank,
    n: usize,
    mut fitted: Option<fs_vfit::discretize::DelayedFilter>,
    mut delay: Option<crate::traveling_wave_line::TravelingWaveLine>,
    zc: f64,
    area: f64,
    texture: &mut TextureDrive,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    use crate::reed_bore::solve_reed_wave;
    let string = assembly.string.expect("checked");
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
    if let Some(bow) = assembly.bow {
        for k in 1..=string.n_modes {
            strike[k - 1] = det::sin(k as f64 * pi * bow.station_frac) / mass_scale;
        }
    }
    let sys = assemble(storage, &zetas, &strike)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let mut x = vec![0.0; 2 * string.n_modes];
    if let Some(pluck) = assembly.pluck {
        for k in 1..=string.n_modes {
            x[2 * (k - 1)] = mass_scale * triangular_pluck_modal(pluck, k);
        }
    }
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    let listener_m = assembly.listener.distance_m;
    let mut p_plus_prev = 5.0;
    let mut out = vec![0.0; n];
    for i in 0..n {
        let u = if let Some(bow) = assembly.bow {
            let v_string: f64 = (0..string.n_modes)
                .map(|k| {
                    let phi = det::sin((k + 1) as f64 * pi * bow.station_frac) / mass_scale;
                    x[2 * k + 1] * phi
                })
                .sum();
            vec![bow_force(
                bow,
                v_string,
                texture.delta_n(bow.velocity_m_s, dt),
            )]
        } else {
            vec![0.0]
        };
        let rec =
            step(&sys, &x, &u, dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        x = rec.x;
        if !assembly.obstacles.is_empty() {
            let extra = modal_contact_forces(string, &assembly.obstacles, &x)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            for k in 0..string.n_modes {
                x[2 * k + 1] += dt * extra[k];
            }
        }
        let mut q_phys = vec![0.0; string.n_modes];
        let mut p_string = 0.0;
        for k in 0..string.n_modes {
            let kf = (k + 1) as f64;
            let q = x[2 * k];
            q_phys[k] = q / mass_scale;
            if k % 2 == 0 {
                let omega = string_mode_omega(string, k + 1);
                let monopole = string.width_m * 2.0 * string.length_m / (kf * pi);
                let h_im = omega * gas.density * monopole / (4.0 * pi * listener_m * mass_scale);
                p_string += h_im * omega * q;
            }
        }
        let fb = bridge_force(string, &q_phys);
        let p_minus = if let Some(line) = fitted.as_ref() {
            line.incoming()
        } else {
            delay.as_ref().expect("line").incoming()
        };
        let u_body = plates.volume_velocity();
        let t = i as f64 * dt;
        let p_plus = if let Some(reed) = assembly.reed {
            let p_m = if reed.attack_s <= 0.0 {
                reed.blowing_pressure_pa
            } else if t >= reed.attack_s {
                reed.blowing_pressure_pa
            } else {
                let frac = t / reed.attack_s;
                reed.blowing_pressure_pa * 0.5 * (1.0 - det::cos(core::f64::consts::PI * frac))
            };
            solve_reed_wave(
                reed,
                gas.density,
                zc,
                0.0,
                p_minus,
                p_m,
                p_plus_prev,
                u_body,
            )?
        } else {
            let blow = assembly.blow.expect("checked");
            let u_blow = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
            p_minus + zc * (u_blow + u_body)
        };
        p_plus_prev = p_plus;
        let p_minus_now = if let Some(line) = fitted.as_mut() {
            line.push(p_plus).map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?
        } else {
            delay.as_mut().expect("line").push(p_plus)
        };
        let p_bore = p_plus + p_minus_now;
        let mut p = p_string + p_bore;
        p += plates.drive_and_radiate(fb + p_bore * area, dt, gas.density, listener_m)?;
        out[i] = p;
    }
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa: out,
        gas: gas.clone(),
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
    texture: Option<ContactTexture>,
    plates: &mut PlateBank,
    obstacles: &[UnilateralObstacle],
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
            texture,
            plates,
            obstacles,
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
            texture,
            plates,
            obstacles,
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
        let mut unused_bodies = PlateBank::default();
        let other = if string.axial_stiffness_n > 0.0 {
            realize_kc_string(
                twin,
                pluck,
                bow,
                None,
                &mut unused_bodies,
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
                &mut unused_bodies,
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
    texture: Option<ContactTexture>,
    plates: &mut PlateBank,
    obstacles: &[UnilateralObstacle],
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    radiation_scale: f64,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut texture = TextureDrive::try_new(texture)?;
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
            let f_bow = bow_force(bow, v_string, texture.delta_n(bow.velocity_m_s, dt));
            for (f, phi) in force.iter_mut().zip(phi_bow.iter()) {
                *f = f_bow * *phi;
            }
        }
        if !obstacles.is_empty() {
            let x: Vec<f64> = model
                .states()
                .iter()
                .flat_map(|s| [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s])
                .collect();
            let extra = modal_contact_forces(string, obstacles, &x)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            for (f, c) in force.iter_mut().zip(extra) {
                *f += c;
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
            p += plates.drive_and_radiate(fb, dt, gas.density, listener_m)?;
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
    texture: Option<ContactTexture>,
    plates: &mut PlateBank,
    obstacles: &[UnilateralObstacle],
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    radiation_scale: f64,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut texture = TextureDrive::try_new(texture)?;
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
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let u = if let Some(bow) = bow {
            let v_string: f64 = (0..string.n_modes)
                .map(|k| {
                    let phi = det::sin((k + 1) as f64 * pi * bow.station_frac) / mass_scale;
                    x[2 * k + 1] * phi
                })
                .sum();
            vec![bow_force(
                bow,
                v_string,
                texture.delta_n(bow.velocity_m_s, dt),
            )]
        } else {
            vec![0.0]
        };
        let rec =
            step(&sys, &x, &u, dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        x = rec.x;
        if !obstacles.is_empty() {
            let extra = modal_contact_forces(string, obstacles, &x)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            for k in 0..string.n_modes {
                x[2 * k + 1] += dt * extra[k];
            }
        }
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
            p += plates.drive_and_radiate(fb, dt, gas.density, listener_m)?;
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

fn bow_force(bow: BowStroke, v_string: f64, normal_delta_n: f64) -> f64 {
    let law = fs_tribo::FrictionLaw::Stribeck {
        static_mu: bow.mu_static,
        kinetic_mu: bow.mu_dynamic,
        characteristic_speed: bow.stribeck_m_s,
        viscous_per_speed: 0.0,
    };
    let normal = (bow.normal_force_n + normal_delta_n).max(0.0);
    // Driven-body sign: + when the driver is faster.
    law.regularized_traction_1d(bow.velocity_m_s - v_string, normal, bow.stribeck_m_s)
        .map(|f| -f)
        .unwrap_or(0.0)
}

/// Declared surface-height drive of the contact normal.
struct TextureDrive {
    inner: Option<TextureInner>,
}

struct TextureInner {
    iface: InterfaceSystemRef,
    moving: UniformSurfaceTrace,
    fixed: UniformSurfaceTrace,
    kn: f64,
    path_m: f64,
}

impl TextureDrive {
    fn try_new(spec: Option<ContactTexture>) -> Result<Self, AcousticRealizeError> {
        let Some(spec) = spec else {
            return Ok(Self { inner: None });
        };
        if !(spec.rms_height_m > 0.0
            && spec.track_length_m > 0.0
            && spec.tangent_stiffness_n_m > 0.0)
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "contact texture must have positive height, track, and stiffness",
            });
        }
        let spectrum = SelfAffinePeriodicProfileSpectrum::new(
            spec.rms_height_m,
            spec.hurst_exponent,
            spec.min_cycles,
            spec.max_cycles,
            spec.phase_seed,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let harmonics = spectrum
            .realize_harmonics()
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let samples = (spec.max_cycles as usize).saturating_mul(8).max(32);
        let moving = PeriodicHarmonicSurface::new(
            "contact.moving",
            "fs-scenario.ContactTexture",
            InputAuthority::CallerDeclared,
            spec.track_length_m,
            samples,
            harmonics,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
        .realize()
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let fixed = PeriodicHarmonicSurface::new(
            "contact.fixed",
            "fs-scenario.ContactTexture",
            InputAuthority::CallerDeclared,
            spec.track_length_m,
            samples,
            Vec::new(),
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
        .realize()
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let iface = InterfaceSystemRef::new(
            "assembly.contact.a-b",
            "assembly.contact.history",
            "fs-scenario.ContactTexture",
            InputAuthority::CallerDeclared,
            InterfaceMedium::Dry,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        Ok(Self {
            inner: Some(TextureInner {
                iface,
                moving,
                fixed,
                kn: spec.tangent_stiffness_n_m,
                path_m: 0.0,
            }),
        })
    }

    fn delta_n(&mut self, speed_m_s: f64, dt: f64) -> f64 {
        let Some(inner) = self.inner.as_mut() else {
            return 0.0;
        };
        inner.path_m += speed_m_s.abs() * dt;
        let motion_a = SurfaceTraceMotion {
            trace: &inner.moving,
            path_coordinate_m: inner.path_m,
            path_speed_m_per_s: speed_m_s,
        };
        let motion_b = SurfaceTraceMotion {
            trace: &inner.fixed,
            path_coordinate_m: 0.0,
            path_speed_m_per_s: 0.0,
        };
        evaluate_point_surface_pair(&inner.iface, motion_a, motion_b)
            .map(|r| inner.kn * r.combined_effective_height_m)
            .unwrap_or(0.0)
    }
}

fn plate_bank(
    soundboard: Option<RadiatingPlate>,
    extra: &[RadiatingPlate],
    plate: Option<ThinPlate>,
) -> Result<PlateBank, AcousticRealizeError> {
    let mut out = PlateBank::default();
    if let Some(mesh) = plate {
        if mesh.geometric_nonlinearity {
            out.vk = Some(VkBody::from_plate(mesh)?);
        } else {
            out.linear.extend(certified_radiators(mesh)?);
        }
    }
    if let Some(spec) = soundboard {
        out.linear
            .push(crate::thin_plate::CompactBody::from_radiator(spec)?);
    }
    for &spec in extra {
        out.linear
            .push(crate::thin_plate::CompactBody::from_radiator(spec)?);
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

fn realize_reed_on_duct(
    duct: &ViscothermalDuct,
    reed: BeatingReed,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let physics = physics_duct(duct)?;
    let termination = match duct.termination {
        WaveguideEnd::Closed => Termination::Closed,
        WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
    };
    refuse_open_nyquist(duct, gas, sample_rate_hz)?;
    realize_reed_bore(
        &physics,
        gas,
        reed,
        termination,
        plates,
        listener_m,
        sample_rate_hz,
        n,
    )
}

fn realize_blown_duct(
    duct: &ViscothermalDuct,
    blow: VolumeVelocityPulse,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
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
    if !plates.is_empty() {
        return realize_blown_with_body(
            &physics,
            blow,
            plates,
            gas,
            termination,
            listener_m,
            sample_rate_hz,
            n,
        );
    }
    // Isolated linear blow: IFFT of the TMM driving point is the
    // exact linear response, including tone-hole shortening. The
    // characteristic line is the time port when a body or valve
    // exchanges volume velocity — it cannot yet replace this oracle.
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

fn realize_blown_with_body(
    physics: &Duct,
    blow: VolumeVelocityPulse,
    plates: &mut PlateBank,
    gas: &GasState,
    termination: Termination,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    use crate::driving_point::characteristic_line;
    use crate::traveling_wave_line::TravelingWaveLine;
    let inlet_r = physics
        .segments
        .first()
        .ok_or(AcousticRealizeError::InvalidDescription {
            what: "duct has no segments",
        })?
        .outlet_radius();
    let area = core::f64::consts::PI * inlet_r * inlet_r;
    let zc = gas.density * gas.sound_speed / area;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut fitted = characteristic_line(physics, gas, termination, sample_rate_hz, n, zc).ok();
    let mut delay = if fitted.is_none() {
        Some(
            TravelingWaveLine::from_duct(physics, gas, termination, sample_rate_hz, n, zc)
                .map_err(|e| match e {
                    crate::traveling_wave_line::TravelingWaveError::Invalid { what } => {
                        AcousticRealizeError::Reed { what }
                    }
                    crate::traveling_wave_line::TravelingWaveError::Duct(d) => {
                        AcousticRealizeError::Duct(d)
                    }
                })?,
        )
    } else {
        None
    };
    let mut out = vec![0.0; n];
    for i in 0..n {
        let t = i as f64 * dt;
        let u_blow = if t < blow.duration_s {
            blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
        } else {
            0.0
        };
        let u_body = plates.volume_velocity();
        let p_minus = if let Some(line) = fitted.as_ref() {
            line.incoming()
        } else {
            delay.as_ref().expect("line").incoming()
        };
        let p_plus = p_minus + zc * (u_blow + u_body);
        let p_minus_now = if let Some(line) = fitted.as_mut() {
            line.push(p_plus).map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?
        } else {
            delay.as_mut().expect("line").push(p_plus)
        };
        let p_bore = p_plus + p_minus_now;
        let mut p_obs = p_bore;
        p_obs += plates.drive_and_radiate(p_bore * area, dt, gas.density, listener_m)?;
        out[i] = p_obs;
    }
    Ok(out)
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
