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
use crate::reed_bore::{blowing_envelope, realize_reed_bore, reed_structural};
use crate::thin_plate::{PlateBank, VkBody, certified_radiators, vk_plate_phs};
use crate::unilateral_contact::{
    modal_contact_forces, modal_friction_forces, modal_hunt_crossley_forces, slit_contact_force,
    slit_lay,
};
/// Sections-per-metre budget for slicing a duct into lumped cells.
const SECTION_BUDGET: f64 = 8.0;

use fs_duct::{Duct, DuctError, HoleState, MAX_RADIATION_KA, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::{GeneralizedMaxwell, RayleighDamping};
use fs_math::c64::C64;
use fs_math::det;
use fs_nlmodal::{
    KcStringParams, assemble_storage, kirchhoff_carrier_moving_end, kirchhoff_carrier_string,
    prestressed_beam_omega,
};
use fs_phs::{
    AcousticSection, AcousticTap, MouthFlange, ViscothermalPin, WallPin, acoustic_chain_mouth_wall,
    bernoulli_volume_flow, join_port, mass_spring_damper, modal_bank, modal_bank_ports,
    quasistatic_aperture_opening, step, step_descriptor, transformer,
};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, BowStroke, ContactTexture, HelmholtzCavity,
    LocallyReactingWall, Pluck, PrestressedString, RadiatingPlate, RayleighParams, ThinPlate,
    UnilateralObstacle, ViscothermalDuct, VolumeVelocityPulse, WaveguideEnd,
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
#[allow(clippy::too_many_lines)] // one coherent realization dispatcher
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
    let dirac_base = assembly.string.is_some_and(|s| s.moving_end);
    let dirac_string_only = dirac_base && assembly.duct.is_none();
    let dirac_string_duct = dirac_base
        && assembly.duct.is_some()
        && assembly.cavity.is_none()
        && assembly
            .duct
            .as_ref()
            .is_some_and(|d| bore_spec(d).is_some());
    if let Some(cavity) = assembly.cavity
        && !dirac_string_only
    {
        bodies.attach_cavity(cavity, &gas)?;
    }
    if dirac_string_duct {
        let string = assembly.string.expect("checked");
        if assembly.pluck.is_none() && assembly.bow.is_none() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a string member requires a pluck or a bow",
            });
        }
        let hist = realize_dirac_join(
            string,
            assembly.pluck,
            assembly.bow,
            &assembly.obstacles,
            assembly.contact_texture,
            &bodies,
            assembly.plate,
            None,
            assembly.duct.as_ref(),
            assembly.reed,
            assembly.blow,
            &gas,
            assembly.listener.distance_m,
            assembly.sample_rate_hz,
            n,
        )?;
        add_in_place(&mut pressure_pa, &hist);
        return Ok(finish_observer(assembly, &gas, pressure_pa));
    }
    if assembly.string.is_some() && assembly.duct.is_some() {
        if assembly
            .duct
            .as_ref()
            .is_some_and(|d| bore_spec(d).is_some())
        {
            return realize_coupled_ode(assembly, &gas, &mut bodies, n);
        }
        return realize_coupled(assembly, &gas, &mut bodies, n);
    }
    if let Some(string) = assembly.string {
        if assembly.pluck.is_none() && assembly.bow.is_none() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a string member requires a pluck or a bow",
            });
        }
        if dirac_string_only {
            let hist = realize_dirac_join(
                string,
                assembly.pluck,
                assembly.bow,
                &assembly.obstacles,
                assembly.contact_texture,
                &bodies,
                assembly.plate,
                assembly.cavity,
                None,
                assembly.reed,
                assembly.blow,
                &gas,
                assembly.listener.distance_m,
                assembly.sample_rate_hz,
                n,
            )?;
            add_in_place(&mut pressure_pa, &hist);
            return Ok(finish_observer(assembly, &gas, pressure_pa));
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
    Ok(finish_observer(assembly, &gas, pressure_pa))
}

/// One-way string (fixed-fixed, `φ(0)=0`) plus an
/// [`acoustic_chain`]. With a plate the bridge force is a leftover
/// port of `transformer(plate, chain)`. Without a plate the
/// support does no work (`φ(0)=0`) so the members share a clock
/// but not a force — the chain is driven only by blow/reed.
/// Shared clock, not FIR.
#[allow(clippy::too_many_lines)] // one coherent coupled realization
#[allow(clippy::needless_range_loop)] // the sample index is the time axis
fn realize_coupled_ode(
    assembly: &AcousticAssembly,
    gas: &GasState,
    plates: &mut PlateBank,
    n: usize,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    let string = assembly.string.expect("checked");
    let duct = assembly.duct.as_ref().expect("checked");
    if assembly.pluck.is_none() && assembly.bow.is_none() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "a string member requires a pluck or a bow",
        });
    }
    validate_string(string, assembly.pluck, assembly.bow)?;
    refuse_open_nyquist(duct, gas, assembly.sample_rate_hz)?;
    let (sections, taps) = bore_spec(duct).ok_or(AcousticRealizeError::InvalidDescription {
        what: "ODE coupled duct expected a cylindrical chain",
    })?;
    let has_plate = !plates.linear.is_empty() || plates.vk.is_some();
    let need_drive = assembly.reed.is_some() || assembly.blow.is_some();
    let inlets = if has_plate && need_drive { 2 } else { 1 };
    let chain = ode_bore_chain(
        &sections,
        duct_mouth(duct.termination),
        inlets,
        &taps,
        gas,
        duct.wall,
    )?;
    let body = if has_plate {
        let plate = if plates.vk.is_some() {
            vk_plate_phs(
                assembly
                    .plate
                    .ok_or(AcousticRealizeError::InvalidDescription {
                        what: "von Karman coupled ODE needs the plate description",
                    })?,
                true,
            )?
        } else {
            plate_modal_ports(plates, true)?.ok_or(AcousticRealizeError::InvalidDescription {
                what: "ODE coupled duct expected a linear plate bank",
            })?
        };
        let face = usize::from(inlets == 2);
        transformer(plate, chain, 1, face, 1.0)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
    } else {
        chain
    };
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    let listener_m = assembly.listener.distance_m;
    let mut texture = TextureDrive::try_new(assembly.contact_texture)?;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let n_plate = if !has_plate {
        0
    } else if let Some(vk) = plates.vk.as_ref() {
        2 * vk.areas.len()
    } else {
        2 * plates.linear.len()
    };
    let mut x_body = vec![0.0; body.state_dim()];
    let mut out = vec![0.0; n];
    if string.axial_stiffness_n > 0.0 {
        let mut members = vec![kc_string_member(
            string,
            assembly.pluck,
            assembly.bow,
            &assembly.obstacles,
            gas,
            listener_m,
            1.0,
        )?];
        if string.polarization_detune > 0.0 {
            let mut twin = string;
            twin.polarization_detune = 0.0;
            twin.tension_n *=
                (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
            members.push(kc_string_member(
                twin,
                assembly.pluck,
                assembly.bow,
                &assembly.obstacles,
                gas,
                listener_m,
                0.85,
            )?);
        }
        for i in 0..n {
            let t = i as f64 * dt;
            let mut p_string = 0.0;
            let mut fb = 0.0;
            for (idx, member) in members.iter_mut().enumerate() {
                let member_string = if idx == 0 {
                    string
                } else {
                    let mut twin = string;
                    twin.polarization_detune = 0.0;
                    twin.tension_n *=
                        (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
                    twin
                };
                let (p_m, q_phys) = step_kc_member(
                    member,
                    assembly.bow,
                    &mut texture,
                    &assembly.obstacles,
                    member_string,
                    gas,
                    dt,
                )?;
                p_string += p_m;
                fb += bridge_force(member_string, &q_phys);
            }
            let rec = step_string_chain(
                &body,
                &x_body,
                fb,
                n_plate,
                assembly.reed,
                assembly.blow,
                gas,
                t,
                dt,
            )?;
            out[i] = p_string
                + plate_chain_radiation(plates, &rec.x, &x_body, n_plate, gas, listener_m, dt)
                + chain_inlet_pressure(&body, &rec.x, n_plate);
            x_body = rec.x;
        }
        return Ok(finish_observer(assembly, gas, out));
    }
    let mut members = vec![linear_string_member(
        string,
        assembly.pluck,
        assembly.bow,
        gas,
        listener_m,
        assembly.sample_rate_hz,
        1.0,
    )?];
    if string.polarization_detune > 0.0 {
        let mut twin = string;
        twin.polarization_detune = 0.0;
        twin.tension_n *= (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
        members.push(linear_string_member(
            twin,
            assembly.pluck,
            assembly.bow,
            gas,
            listener_m,
            assembly.sample_rate_hz,
            0.85,
        )?);
    }
    for i in 0..n {
        let t = i as f64 * dt;
        let mut p_string = 0.0;
        let mut fb = 0.0;
        for (idx, member) in members.iter_mut().enumerate() {
            let member_string = if idx == 0 {
                string
            } else {
                let mut twin = string;
                twin.polarization_detune = 0.0;
                twin.tension_n *=
                    (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
                twin
            };
            p_string += step_linear_member(
                member,
                assembly.bow,
                &mut texture,
                &assembly.obstacles,
                member_string,
                dt,
            )?;
            let q_phys: Vec<f64> = member
                .model
                .states()
                .iter()
                .map(|s| s.displacement_m_sqrt_kg / mass_scale)
                .collect();
            fb += bridge_force(member_string, &q_phys);
        }
        let rec = step_string_chain(
            &body,
            &x_body,
            fb,
            n_plate,
            assembly.reed,
            assembly.blow,
            gas,
            t,
            dt,
        )?;
        out[i] = p_string
            + plate_chain_radiation(plates, &rec.x, &x_body, n_plate, gas, listener_m, dt)
            + chain_inlet_pressure(&body, &rec.x, n_plate);
        x_body = rec.x;
    }
    Ok(finish_observer(assembly, gas, out))
}

#[allow(clippy::too_many_arguments)] // one coherent chain step record
fn step_string_chain(
    body: &fs_phs::PortHamiltonian,
    x: &[f64],
    f_bridge: f64,
    n_plate: usize,
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
) -> Result<fs_phs::StepRecord, AcousticRealizeError> {
    if n_plate == 0 {
        return step_chain_inlet(body, x, reed, blow, gas, t, dt);
    }
    step_plate_chain(body, x, f_bridge, reed, blow, gas, t, dt)
}

fn chain_inlet_pressure(body: &fs_phs::PortHamiltonian, x: &[f64], n_plate: usize) -> f64 {
    if n_plate != 0 {
        return 0.0;
    }
    body.output(x).first().copied().unwrap_or(0.0)
}

fn step_chain_inlet(
    body: &fs_phs::PortHamiltonian,
    x: &[f64],
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
) -> Result<fs_phs::StepRecord, AcousticRealizeError> {
    let m = body.port_dim();
    let y = body.output(x);
    let mut u = vec![0.0; m];
    if m > 0 {
        if let Some(reed) = reed {
            let p_m = blowing_envelope(reed, t);
            let p_bore = y.first().copied().unwrap_or(0.0);
            let h = quasistatic_aperture_opening(
                reed.rest_opening_m,
                reed.closing_pressure_pa,
                p_m - p_bore,
            );
            u[0] = bernoulli_volume_flow(reed.width_m, h, p_m - p_bore, gas.density);
        } else if let Some(blow) = blow {
            u[0] = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
        }
    }
    step(body, x, &u, dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
}

#[allow(clippy::too_many_arguments)] // one coherent chain step record
fn step_plate_chain(
    body: &fs_phs::PortHamiltonian,
    x: &[f64],
    f_bridge: f64,
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
) -> Result<fs_phs::StepRecord, AcousticRealizeError> {
    let m = body.port_dim();
    let y = body.output(x);
    let mut u = vec![0.0; m];
    if m == 0 {
        return step(body, x, &u, dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()));
    }
    u[0] = f_bridge;
    if m > 1 {
        if let Some(reed) = reed {
            let p_m = blowing_envelope(reed, t);
            let p_bore = y.get(1).copied().unwrap_or(0.0);
            let h = quasistatic_aperture_opening(
                reed.rest_opening_m,
                reed.closing_pressure_pa,
                p_m - p_bore,
            );
            u[1] = bernoulli_volume_flow(reed.width_m, h, p_m - p_bore, gas.density);
        } else if let Some(blow) = blow {
            u[1] = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
        }
    }
    step(body, x, &u, dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
}

fn plate_chain_radiation(
    plates: &PlateBank,
    x1: &[f64],
    x0: &[f64],
    n_plate: usize,
    gas: &GasState,
    listener_m: f64,
    dt: f64,
) -> f64 {
    let mut p = 0.0;
    if let Some(vk) = plates.vk.as_ref() {
        for (k, &area) in vk.areas.iter().enumerate() {
            let acc = (x1.get(2 * k + 1).copied().unwrap_or(0.0)
                - x0.get(2 * k + 1).copied().unwrap_or(0.0))
                / dt;
            p += gas.density * area * acc / (2.0 * core::f64::consts::PI * listener_m);
        }
        return p;
    }
    let n = n_plate / 2;
    for (k, body) in plates.linear.iter().take(n).enumerate() {
        let acc = (x1.get(2 * k + 1).copied().unwrap_or(0.0)
            - x0.get(2 * k + 1).copied().unwrap_or(0.0))
            / dt;
        p += body.radiate(acc, gas.density, listener_m);
    }
    p
}

/// One shared clock: string, plate, and duct exchange force and flow
/// every sample via the characteristic FIR. Used when the duct is
/// not a cylindrical chain.
#[allow(clippy::too_many_lines)] // one coherent coupled realization
#[allow(clippy::needless_range_loop)] // the sample index is the time axis
fn realize_coupled(
    assembly: &AcousticAssembly,
    gas: &GasState,
    plates: &mut PlateBank,
    n: usize,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    use crate::driving_point::characteristic_line;
    use crate::reed_bore::solve_reed_wave;
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
        WaveguideEnd::FlangedOpen => Termination::FlangedOpen,
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
    let mut fitted = characteristic_line(
        &physics,
        gas,
        termination,
        assembly.sample_rate_hz,
        n,
        zc,
        wall_pin(duct.wall).as_ref(),
    )
    .map_err(map_drive)?;
    let mut texture = TextureDrive::try_new(assembly.contact_texture)?;
    if string.axial_stiffness_n > 0.0 {
        return realize_coupled_kc(assembly, gas, plates, n, fitted, zc, area, &mut texture);
    }
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let mut members = vec![linear_string_member(
        string,
        assembly.pluck,
        assembly.bow,
        gas,
        listener_m,
        assembly.sample_rate_hz,
        1.0,
    )?];
    if string.polarization_detune > 0.0 {
        let mut twin = string;
        twin.polarization_detune = 0.0;
        twin.tension_n *= (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
        members.push(linear_string_member(
            twin,
            assembly.pluck,
            assembly.bow,
            gas,
            listener_m,
            assembly.sample_rate_hz,
            0.85,
        )?);
    }
    let mut p_plus_prev = 5.0;
    let mut out = vec![0.0; n];
    #[allow(clippy::needless_range_loop)] // the sample index is the time axis
    for i in 0..n {
        let mut p_string = 0.0;
        let mut fb = 0.0;
        for (idx, member) in members.iter_mut().enumerate() {
            let member_string = if idx == 0 {
                string
            } else {
                let mut twin = string;
                twin.polarization_detune = 0.0;
                twin.tension_n *=
                    (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
                twin
            };
            p_string += step_linear_member(
                member,
                assembly.bow,
                &mut texture,
                &assembly.obstacles,
                member_string,
                dt,
            )?;
            let q_phys: Vec<f64> = member
                .model
                .states()
                .iter()
                .map(|s| s.displacement_m_sqrt_kg / mass_scale)
                .collect();
            fb += bridge_force(member_string, &q_phys);
        }
        let p_minus = fitted.incoming();
        let u_body = plates.volume_velocity();
        let t = i as f64 * dt;
        let p_plus = if let Some(reed) = assembly.reed {
            // No ramp configured OR ramp finished: full pressure.
            let p_m = if reed.attack_s <= 0.0 || t >= reed.attack_s {
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
        let p_minus_now = fitted
            .push(p_plus)
            .map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?;
        let p_bore = p_plus + p_minus_now;
        let mut p = p_string + p_bore;
        p += plates.drive_and_radiate(fb + p_bore * area, dt, gas.density, listener_m)?;
        out[i] = p;
    }
    Ok(finish_observer(assembly, gas, out))
}

/// Shared clock with a Kirchhoff–Carrier string (EA > 0).
#[allow(clippy::too_many_arguments)]
fn realize_coupled_kc(
    assembly: &AcousticAssembly,
    gas: &GasState,
    plates: &mut PlateBank,
    n: usize,
    mut fitted: fs_vfit::discretize::DelayedFilter,
    zc: f64,
    area: f64,
    texture: &mut TextureDrive,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    use crate::reed_bore::solve_reed_wave;
    let string = assembly.string.expect("checked");
    let listener_m = assembly.listener.distance_m;
    let mut members = vec![kc_string_member(
        string,
        assembly.pluck,
        assembly.bow,
        &assembly.obstacles,
        gas,
        listener_m,
        1.0,
    )?];
    if string.polarization_detune > 0.0 {
        let mut twin = string;
        twin.polarization_detune = 0.0;
        twin.tension_n *= (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
        members.push(kc_string_member(
            twin,
            assembly.pluck,
            assembly.bow,
            &assembly.obstacles,
            gas,
            listener_m,
            0.85,
        )?);
    }
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    let mut p_plus_prev = 5.0;
    let mut out = vec![0.0; n];
    #[allow(clippy::needless_range_loop)] // the sample index is the time axis
    for i in 0..n {
        let mut p_string = 0.0;
        let mut fb = 0.0;
        for (idx, member) in members.iter_mut().enumerate() {
            let member_string = if idx == 0 {
                string
            } else {
                let mut twin = string;
                twin.polarization_detune = 0.0;
                twin.tension_n *=
                    (1.0 + string.polarization_detune) * (1.0 + string.polarization_detune);
                twin
            };
            let (p_m, q_phys) = step_kc_member(
                member,
                assembly.bow,
                texture,
                &assembly.obstacles,
                member_string,
                gas,
                dt,
            )?;
            p_string += p_m;
            fb += bridge_force(member_string, &q_phys);
        }
        let p_minus = fitted.incoming();
        let u_body = plates.volume_velocity();
        let t = i as f64 * dt;
        let p_plus = if let Some(reed) = assembly.reed {
            // No ramp configured OR ramp finished: full pressure.
            let p_m = if reed.attack_s <= 0.0 || t >= reed.attack_s {
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
        let p_minus_now = fitted
            .push(p_plus)
            .map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?;
        let p_bore = p_plus + p_minus_now;
        let mut p = p_string + p_bore;
        p += plates.drive_and_radiate(fb + p_bore * area, dt, gas.density, listener_m)?;
        out[i] = p;
    }
    Ok(finish_observer(assembly, gas, out))
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

fn map_drive(err: crate::driving_point::DrivingPointError) -> AcousticRealizeError {
    match err {
        crate::driving_point::DrivingPointError::Invalid { what } => {
            AcousticRealizeError::Reed { what }
        }
        crate::driving_point::DrivingPointError::Duct(d) => AcousticRealizeError::Duct(d),
        crate::driving_point::DrivingPointError::Realize(_) => AcousticRealizeError::Reed {
            what: "characteristic realization refused",
        },
        crate::driving_point::DrivingPointError::Discrete(_) => AcousticRealizeError::Reed {
            what: "characteristic line left the finite set",
        },
    }
}

fn gas_state(ambient: AmbientGas) -> Result<GasState, AcousticRealizeError> {
    if !(ambient.relative_humidity >= 0.0 && ambient.relative_humidity <= 1.0) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "relative humidity must be an explicit fraction in [0, 1]",
        });
    }
    GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        ambient.temperature_k,
        ambient.pressure_pa,
    )
    .map_err(|e| AcousticRealizeError::Ambient(e.to_string()))
}

fn finish_observer(
    assembly: &AcousticAssembly,
    gas: &GasState,
    mut pressure_pa: Vec<f64>,
) -> RealizedAssembly {
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    crate::air_path::absorb_pressure_history(
        &mut pressure_pa,
        dt,
        assembly.listener.distance_m,
        gas,
        assembly.ambient.relative_humidity,
    );
    RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa,
        gas: *gas,
    }
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

/// Compact observer transfer `Im H` for sine mode `k` (1-based).
///
/// Odd modes are free-space monopoles (`p ∝ ρ A ÿ / 4π r`). Even
/// modes have vanishing monopole area; they radiate as compact
/// dipoles (`p ∝ ρ Π̈ / 4π r c`). Both are the same Green's-function
/// family — not a BEM field and not a 3-D jet.
#[allow(clippy::too_many_arguments)] // one coherent field record
fn string_mode_h_im(
    string: PrestressedString,
    k: usize,
    omega: f64,
    mass_scale: f64,
    radiation_scale: f64,
    rho: f64,
    sound_speed: f64,
    listener_m: f64,
) -> f64 {
    let kf = k as f64;
    let pi = core::f64::consts::PI;
    if k % 2 == 1 {
        let area = string.width_m * 2.0 * string.length_m / (kf * pi);
        radiation_scale * omega * rho * area / (4.0 * pi * listener_m * mass_scale)
    } else {
        // ∫_0^L sin(kπx/L) (x − L/2) dx = −L²/(kπ) for even k.
        let moment = string.width_m * string.length_m * string.length_m / (kf * pi);
        radiation_scale * omega * rho * moment / (4.0 * pi * listener_m * sound_speed * mass_scale)
    }
}

fn assemble_kc(
    string: PrestressedString,
    storage: fs_nlmodal::SosModalStorage,
    zetas: &[f64],
    obstacles: &[UnilateralObstacle],
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    use crate::unilateral_contact::wrap_modal_contact;
    let n = string.n_modes;
    let omegas = storage.omegas.clone();
    let wrapped = wrap_modal_contact(Box::new(storage), string, obstacles)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let dim = 2 * n;
    let mut g = vec![0.0; dim * n];
    for k in 0..n {
        g[(2 * k + 1) * n + k] = 1.0;
    }
    assemble_storage(n, &omegas, zetas, n, g, wrapped)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
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
    let internal = prony_internal_zeta(string, omega);
    Ok(internal + stokes + bend)
}

/// Authored `ζ` at the fundamental becomes one Prony branch; higher
/// modes see `η(ω)/2` from that branch, not a constant ratio.
fn prony_internal_zeta(string: PrestressedString, omega: f64) -> f64 {
    let z0 = string.damping_ratio;
    if !(z0 > 0.0 && z0 < 0.49) {
        return z0.max(0.0);
    }
    let omega1 = string_mode_omega(string, 1);
    let eta0 = (2.0 * z0).min(0.99);
    GeneralizedMaxwell::matching_loss(1.0, omega1.max(1.0), eta0)
        .map_or(z0, |gm| 0.5 * gm.loss_factor(omega))
}

/// Two-way Dirac realize: moving-end waveguide ⊕ plate ⊕ optional
/// Helmholtz cavity or [`acoustic_chain`] duct. `EA > 0` uses
/// [`kirchhoff_carrier_moving_end`] on the same free-fixed port.
/// Bow, obstacle stations, and blow/reed are leftover ports of
/// [`join_port`]. Von Karman is the same plate pHS with a quartic
/// storage. Hunt–Crossley and Stribeck stay port forces, not terms
/// in `H`.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)] // one coherent realization
fn realize_dirac_join(
    string: PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
    obstacles: &[UnilateralObstacle],
    texture: Option<ContactTexture>,
    plates: &PlateBank,
    plate_spec: Option<ThinPlate>,
    cavity: Option<HelmholtzCavity>,
    duct: Option<&ViscothermalDuct>,
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    validate_string(string, pluck, bow)?;
    let n_s = string.n_modes.max(1);
    let zetas: Vec<f64> = (0..n_s)
        .map(|k| {
            let omega = moving_end_omega(string, k);
            mode_zeta(string, omega, gas)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let obs_stations: Vec<f64> = obstacles
        .iter()
        .flat_map(|o| o.stations.iter().copied())
        .collect();
    let waveguide = moving_end_string_phs(string, &zetas, bow, &obs_stations)?;
    let n_string = waveguide.state_dim();
    let need_area = cavity.is_some() || duct.is_some();
    let body = if plates.vk.is_some() {
        let spec = plate_spec.ok_or(AcousticRealizeError::InvalidDescription {
            what: "von Karman Dirac join needs the plate description",
        })?;
        Some(vk_plate_phs(spec, need_area)?)
    } else {
        plate_modal_ports(plates, need_area)?
    };
    let drive_inlet = reed.is_some() || blow.is_some();
    let (waveguide, sys, ode_start) = match (body, cavity, duct) {
        (None, _, Some(line)) => {
            refuse_open_nyquist(line, gas, sample_rate_hz)?;
            let (sections, taps) =
                bore_spec(line).ok_or(AcousticRealizeError::InvalidDescription {
                    what: "Dirac duct expected a cylindrical chain",
                })?;
            let inlets = if drive_inlet { 2 } else { 1 };
            let chain = ode_bore_chain(
                &sections,
                duct_mouth(line.termination),
                inlets,
                &taps,
                gas,
                line.wall,
            )?;
            let chain_face = usize::from(inlets == 2);
            let radius = sections.first().map_or(0.0, |s| s.radius);
            let area = core::f64::consts::PI * radius * radius;
            (
                Some(
                    transformer(waveguide, chain, 0, chain_face, area)
                        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
                ),
                None,
                0,
            )
        }
        (None, _, None) => (Some(waveguide), None, 1),
        (Some(plate), None, None) => (
            None,
            Some(
                join_port(waveguide, plate, 0, 0)
                    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
            ),
            0,
        ),
        (Some(plate), Some(cav), None) => {
            let flow = helmholtz_flow_cavity(cav, gas)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            let plate_cav = transformer(plate, flow, 1, 0, 1.0)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            (
                None,
                Some(
                    join_port(waveguide, plate_cav, 0, 0)
                        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
                ),
                0,
            )
        }
        (Some(plate), None, Some(line)) => {
            refuse_open_nyquist(line, gas, sample_rate_hz)?;
            let (sections, taps) =
                bore_spec(line).ok_or(AcousticRealizeError::InvalidDescription {
                    what: "Dirac duct expected a cylindrical chain",
                })?;
            let inlets = if drive_inlet { 2 } else { 1 };
            let chain = ode_bore_chain(
                &sections,
                duct_mouth(line.termination),
                inlets,
                &taps,
                gas,
                line.wall,
            )?;
            let chain_face = usize::from(inlets == 2);
            let plate_line = transformer(plate, chain, 1, chain_face, 1.0)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            (
                None,
                Some(
                    join_port(waveguide, plate_line, 0, 0)
                        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
                ),
                0,
            )
        }
        (Some(_), Some(_), Some(_)) => {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "Dirac join takes a cavity or a duct, not both",
            });
        }
    };
    let mut x = match &sys {
        Some(s) => vec![0.0; s.state_dim()],
        None => vec![
            0.0;
            waveguide
                .as_ref()
                .map_or(n_string, fs_phs::PortHamiltonian::state_dim)
        ],
    };
    if let Some(pluck) = pluck {
        for k in 0..n_s {
            x[2 * k] = free_fixed_pluck_modal(pluck, string, k);
        }
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut texture = TextureDrive::try_new(texture)?;
    let mut out = Vec::with_capacity(n);
    let phi0 = det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m));
    let pi = core::f64::consts::PI;
    for i in 0..n {
        let t = i as f64 * dt;
        let rec_x = if let Some(s) = &sys {
            let u = leftover_u(
                s,
                &x,
                n_s,
                string,
                bow,
                obstacles,
                reed,
                blow,
                gas,
                t,
                dt,
                &mut texture,
            )?;
            step_descriptor(s, &x, &u, dt)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
                .x
        } else {
            let wg = waveguide.as_ref().expect("waveguide-only clock");
            let u = leftover_u_ode(
                wg,
                &x,
                n_s,
                string,
                bow,
                obstacles,
                reed,
                blow,
                gas,
                t,
                dt,
                &mut texture,
                ode_start,
            )?;
            step(wg, &x, &u, dt)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
                .x
        };
        let mut p = 0.0;
        for k in 0..n_s {
            let acc = (rec_x[2 * k + 1] - x[2 * k + 1]) / dt;
            let kappa = (k as f64 + 0.5) * pi / string.length_m;
            let integral = phi0 * (if k % 2 == 0 { 1.0 } else { -1.0 }) / kappa;
            p += gas.density * string.width_m * acc * integral / (4.0 * pi * listener_m);
        }
        if !plates.linear.is_empty() {
            let base = n_string;
            for (i, body) in plates.linear.iter().enumerate() {
                let acc = (rec_x[base + 2 * i + 1] - x[base + 2 * i + 1]) / dt;
                p += body.radiate(acc, gas.density, listener_m);
            }
        }
        if let Some(vk) = plates.vk.as_ref() {
            let base = n_string;
            for (k, &area) in vk.areas.iter().enumerate() {
                let acc = (rec_x[base + 2 * k + 1] - x[base + 2 * k + 1]) / dt;
                p += gas.density * area * acc / (2.0 * pi * listener_m);
            }
        }
        if plates.linear.is_empty()
            && plates.vk.is_none()
            && duct.is_some()
            && let Some(wg) = waveguide.as_ref()
        {
            p += wg.effort(&rec_x).get(n_string).copied().unwrap_or(0.0);
        }
        if !p.is_finite() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "Dirac observer pressure left the finite set",
            });
        }
        x = rec_x;
        out.push(p);
    }
    Ok(out)
}

fn moving_end_string_phs(
    string: PrestressedString,
    zetas: &[f64],
    bow: Option<BowStroke>,
    obs_stations: &[f64],
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    let n = string.n_modes.max(1);
    let phi0 = det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m));
    let m = 1 + usize::from(bow.is_some()) + obs_stations.len();
    let mut g = vec![0.0; (2 * n) * m];
    for k in 0..n {
        let kappa = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
        g[(2 * k + 1) * m] = phi0;
        let mut col = 1;
        if let Some(bow) = bow {
            g[(2 * k + 1) * m + col] = phi0 * det::cos(kappa * bow.station_frac * string.length_m);
            col += 1;
        }
        for &s in obs_stations {
            g[(2 * k + 1) * m + col] = phi0 * det::cos(kappa * s.clamp(0.0, 1.0) * string.length_m);
            col += 1;
        }
    }
    if string.axial_stiffness_n > 0.0 {
        let mut storage = kirchhoff_carrier_moving_end(
            &KcStringParams {
                length: string.length_m,
                tension: string.tension_n,
                lin_density: string.lin_density_kg_m,
                ea: string.axial_stiffness_n,
            },
            n,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        if string.bending_stiffness_n_m2 > 0.0 {
            for (k, w) in storage.omegas.iter_mut().enumerate() {
                *w = moving_end_omega(string, k);
            }
        }
        let omegas = storage.omegas.clone();
        assemble_storage(n, &omegas, zetas, m, g, Box::new(storage))
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
    } else {
        let omegas: Vec<f64> = (0..n).map(|k| moving_end_omega(string, k)).collect();
        let mut drives: Vec<Vec<f64>> = Vec::with_capacity(m);
        for p in 0..m {
            drives.push((0..n).map(|k| g[(2 * k + 1) * m + p]).collect());
        }
        let refs: Vec<&[f64]> = drives.iter().map(Vec::as_slice).collect();
        modal_bank_ports(&omegas, zetas, &refs)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
    }
}

#[allow(clippy::too_many_arguments)] // one coherent solver probe
#[allow(clippy::unnecessary_wraps)] // uniform Result surface with the ODE twin
fn leftover_u(
    sys: &fs_phs::DescriptorPortHamiltonian,
    x: &[f64],
    n_s: usize,
    string: PrestressedString,
    bow: Option<BowStroke>,
    obstacles: &[UnilateralObstacle],
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
    texture: &mut TextureDrive,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let m = sys.port_dim();
    if m == 0 {
        return Ok(Vec::new());
    }
    let y = sys.output(x);
    let mut u = vec![0.0; m];
    let mut col = 0;
    if let Some(bow) = bow {
        let v = y.get(col).copied().unwrap_or(0.0);
        u[col] = bow_force(bow, v, texture.delta_n(bow.velocity_m_s - v, dt));
        col += 1;
    }
    for spec in obstacles {
        for (i, &s) in spec.stations.iter().enumerate() {
            let opening = moving_end_opening(string, x, n_s, s);
            let v = y.get(col).copied().unwrap_or(0.0);
            let gap = spec.gaps_m.get(i).copied().unwrap_or(0.0);
            let pen = (opening - gap).max(0.0);
            let fel = -spec.stiffness * det::pow(pen, spec.alpha);
            u[col] = fel - spec.internal_loss * spec.stiffness * det::pow(pen, spec.alpha) * v;
            col += 1;
        }
    }
    if col < m {
        if let Some(reed) = reed {
            let p_m = blowing_envelope(reed, t);
            let p_bore = y.get(col).copied().unwrap_or(0.0);
            let h = quasistatic_aperture_opening(
                reed.rest_opening_m,
                reed.closing_pressure_pa,
                p_m - p_bore,
            );
            u[col] = bernoulli_volume_flow(reed.width_m, h, p_m - p_bore, gas.density);
        } else if let Some(blow) = blow {
            u[col] = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
        }
    }
    Ok(u)
}

#[allow(clippy::too_many_arguments)] // one coherent solver probe
#[allow(clippy::unnecessary_wraps)] // uniform Result surface with the DAE twin
fn leftover_u_ode(
    sys: &fs_phs::PortHamiltonian,
    x: &[f64],
    n_s: usize,
    string: PrestressedString,
    bow: Option<BowStroke>,
    obstacles: &[UnilateralObstacle],
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
    texture: &mut TextureDrive,
    start_col: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let m = sys.port_dim();
    let y = sys.output(x);
    let mut u = vec![0.0; m];
    let mut col = start_col;
    if col >= m {
        return Ok(u);
    }
    if let Some(bow) = bow {
        let v = y.get(col).copied().unwrap_or(0.0);
        u[col] = bow_force(bow, v, texture.delta_n(bow.velocity_m_s - v, dt));
        col += 1;
    }
    for spec in obstacles {
        for (i, &s) in spec.stations.iter().enumerate() {
            if col >= m {
                break;
            }
            let opening = moving_end_opening(string, x, n_s, s);
            let v = y.get(col).copied().unwrap_or(0.0);
            let gap = spec.gaps_m.get(i).copied().unwrap_or(0.0);
            let pen = (opening - gap).max(0.0);
            let fel = -spec.stiffness * det::pow(pen, spec.alpha);
            u[col] = fel - spec.internal_loss * spec.stiffness * det::pow(pen, spec.alpha) * v;
            col += 1;
        }
    }
    if col < m {
        if let Some(reed) = reed {
            let p_m = blowing_envelope(reed, t);
            let p_bore = y.get(col).copied().unwrap_or(0.0);
            let h = quasistatic_aperture_opening(
                reed.rest_opening_m,
                reed.closing_pressure_pa,
                p_m - p_bore,
            );
            u[col] = bernoulli_volume_flow(reed.width_m, h, p_m - p_bore, gas.density);
        } else if let Some(blow) = blow {
            u[col] = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
        }
    }
    Ok(u)
}

fn moving_end_opening(string: PrestressedString, x: &[f64], n_s: usize, station: f64) -> f64 {
    let phi0 = det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m));
    let mut y = 0.0;
    for k in 0..n_s {
        let kappa = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
        let phi = phi0 * det::cos(kappa * station.clamp(0.0, 1.0) * string.length_m);
        y += x.get(2 * k).copied().unwrap_or(0.0) * phi;
    }
    y
}

fn moving_end_omega(string: PrestressedString, k_zero: usize) -> f64 {
    let n = k_zero as f64 + 0.5;
    let wave = det::sqrt(string.tension_n / string.lin_density_kg_m);
    let omega = n * core::f64::consts::PI * wave / string.length_m;
    if !(string.bending_stiffness_n_m2 > 0.0) {
        return omega;
    }
    let inharm = core::f64::consts::PI * core::f64::consts::PI * string.bending_stiffness_n_m2
        / (string.tension_n * string.length_m * string.length_m);
    omega * det::sqrt(1.0 + inharm * n * n)
}

fn free_fixed_pluck_modal(pluck: Pluck, string: PrestressedString, k_zero: usize) -> f64 {
    let l = string.length_m;
    let s = pluck.station_frac * l;
    let h = pluck.height_m;
    let mu = string.lin_density_kg_m;
    let kappa = (k_zero as f64 + 0.5) * core::f64::consts::PI / l;
    let phi0 = det::sqrt(2.0 / (mu * l));
    let nq = 48usize;
    let mut acc = 0.0;
    for i in 0..nq {
        let x = (i as f64 + 0.5) * l / nq as f64;
        let w = if x < s {
            h * x / s.max(1.0e-18)
        } else {
            h * (l - x) / (l - s).max(1.0e-18)
        };
        acc += w * det::cos(kappa * x);
    }
    mu * phi0 * acc * (l / nq as f64)
}

fn plate_modal_ports(
    plates: &PlateBank,
    with_area_port: bool,
) -> Result<Option<fs_phs::PortHamiltonian>, AcousticRealizeError> {
    if plates.linear.is_empty() {
        return Ok(None);
    }
    let omegas: Vec<f64> = plates.linear.iter().map(|b| b.omega).collect();
    let zetas: Vec<f64> = plates.linear.iter().map(|b| b.zeta).collect();
    let drives: Vec<f64> = plates
        .linear
        .iter()
        .map(|b| 1.0 / b.mass_kg.sqrt().max(1.0e-18))
        .collect();
    let sys = if with_area_port {
        let areas: Vec<f64> = plates.linear.iter().map(|b| b.area_m2).collect();
        modal_bank_ports(&omegas, &zetas, &[&drives, &areas])
    } else {
        modal_bank(&omegas, &zetas, &drives)
    }
    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    Ok(Some(sys))
}

fn helmholtz_flow_cavity(
    cavity: HelmholtzCavity,
    gas: &GasState,
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    use fs_phs::{MouthFlange, compact_radiation_impedance, helmholtz_resonator_flow};
    let pi = core::f64::consts::PI;
    let neck_area = pi * cavity.neck_radius_m * cavity.neck_radius_m;
    let l_eff = cavity.neck_length_m + 2.0 * (8.0 / (3.0 * pi)) * cavity.neck_radius_m;
    let omega0 = gas.sound_speed * (neck_area / (cavity.volume_m3 * l_eff)).sqrt();
    let r_rad = compact_radiation_impedance(
        gas.density,
        gas.sound_speed,
        cavity.neck_radius_m,
        omega0,
        MouthFlange::Unflanged,
    )
    .map_or(0.0, |(r, _)| r);
    helmholtz_resonator_flow(
        cavity.volume_m3,
        cavity.neck_radius_m,
        cavity.neck_length_m,
        gas.density,
        gas.sound_speed,
        r_rad,
    )
    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
}

#[allow(clippy::too_many_arguments)] // one coherent realization record
#[allow(clippy::needless_range_loop)] // modal index spans state and shape arrays
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
    let twin_tension = if string.polarization_detune > 0.0 {
        Some((1.0 + string.polarization_detune) * (1.0 + string.polarization_detune))
    } else {
        None
    };
    if string.axial_stiffness_n > 0.0 {
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
            twin_tension,
        )
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
            twin_tension,
        )
    }
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
    twin_tension: Option<f64>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut texture = TextureDrive::try_new(texture)?;
    let mut members = vec![linear_string_member(
        string,
        pluck,
        bow,
        gas,
        listener_m,
        sample_rate_hz,
        1.0,
    )?];
    if let Some(scale) = twin_tension {
        let mut twin = string;
        twin.polarization_detune = 0.0;
        twin.tension_n *= scale;
        members.push(linear_string_member(
            twin,
            pluck,
            bow,
            gas,
            listener_m,
            sample_rate_hz,
            0.85,
        )?);
    }
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let mut out = Vec::with_capacity(n);
    let dt = 1.0 / f64::from(sample_rate_hz);
    for _ in 0..n {
        let mut p = 0.0;
        let mut fb = 0.0;
        for (idx, member) in members.iter_mut().enumerate() {
            let member_string = if idx == 0 {
                string
            } else {
                let mut twin = string;
                twin.polarization_detune = 0.0;
                if let Some(scale) = twin_tension {
                    twin.tension_n *= scale;
                }
                twin
            };
            p += step_linear_member(member, bow, &mut texture, obstacles, member_string, dt)?;
            if !plates.is_empty() {
                let q_phys: Vec<f64> = member
                    .model
                    .states()
                    .iter()
                    .map(|s| s.displacement_m_sqrt_kg / mass_scale)
                    .collect();
                fb += bridge_force(member_string, &q_phys);
            }
        }
        if !plates.is_empty() {
            p += plates.drive_and_radiate(fb, dt, gas.density, listener_m)?;
        }
        out.push(p);
    }
    Ok(out)
}

struct LinearMember {
    model: ModalAcousticTimeModel,
    phi_bow: Vec<f64>,
}

fn linear_string_member(
    string: PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    radiation_scale: f64,
) -> Result<LinearMember, AcousticRealizeError> {
    let pi = core::f64::consts::PI;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let mut modes = Vec::with_capacity(string.n_modes);
    let mut states = Vec::with_capacity(string.n_modes);
    let mut phi_bow = vec![0.0; string.n_modes];
    for k in 1..=string.n_modes {
        let omega = string_mode_omega(string, k);
        let q_phys = pluck.map_or(0.0, |p| triangular_pluck_modal(p, k));
        let h_im = string_mode_h_im(
            string,
            k,
            omega,
            mass_scale,
            radiation_scale,
            gas.density,
            gas.sound_speed,
            listener_m,
        );
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: mode_zeta(string, omega, gas)?,
            pressure_per_modal_velocity: C64::new(0.0, h_im),
        });
        states.push(ModalAcousticState {
            displacement_m_sqrt_kg: mass_scale * q_phys,
            velocity_m_sqrt_kg_per_s: 0.0,
        });
        if let Some(bow) = bow {
            phi_bow[k - 1] = det::sin(k as f64 * pi * bow.station_frac) / mass_scale;
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
    Ok(LinearMember { model, phi_bow })
}

fn step_linear_member(
    member: &mut LinearMember,
    bow: Option<BowStroke>,
    texture: &mut TextureDrive,
    obstacles: &[UnilateralObstacle],
    string: PrestressedString,
    dt: f64,
) -> Result<f64, AcousticRealizeError> {
    let mut force = vec![0.0; string.n_modes];
    if let Some(bow) = bow {
        let v_string: f64 = member
            .model
            .states()
            .iter()
            .zip(member.phi_bow.iter())
            .map(|(s, phi)| s.velocity_m_sqrt_kg_per_s * phi)
            .sum();
        let f_bow = bow_force(
            bow,
            v_string,
            texture.delta_n(bow.velocity_m_s - v_string, dt),
        );
        for (f, phi) in force.iter_mut().zip(member.phi_bow.iter()) {
            *f = f_bow * *phi;
        }
    }
    if !obstacles.is_empty() {
        let x: Vec<f64> = member
            .model
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
    let frame = member
        .model
        .step(&force)
        .map_err(AcousticRealizeError::Modal)?;
    Ok(frame.observer_pressure_pa)
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
    twin_tension: Option<f64>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut texture = TextureDrive::try_new(texture)?;
    let mut members = vec![kc_string_member(
        string, pluck, bow, obstacles, gas, listener_m, 1.0,
    )?];
    if let Some(scale) = twin_tension {
        let mut twin = string;
        twin.polarization_detune = 0.0;
        twin.tension_n *= scale;
        members.push(kc_string_member(
            twin, pluck, bow, obstacles, gas, listener_m, 0.85,
        )?);
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let mut p = 0.0;
        let mut fb = 0.0;
        for (idx, member) in members.iter_mut().enumerate() {
            let member_string = if idx == 0 {
                string
            } else {
                let mut twin = string;
                twin.polarization_detune = 0.0;
                if let Some(scale) = twin_tension {
                    twin.tension_n *= scale;
                }
                twin
            };
            let (p_m, q_phys) =
                step_kc_member(member, bow, &mut texture, obstacles, member_string, gas, dt)?;
            p += p_m;
            if !plates.is_empty() {
                fb += bridge_force(member_string, &q_phys);
            }
        }
        if !plates.is_empty() {
            p += plates.drive_and_radiate(fb, dt, gas.density, listener_m)?;
        }
        out.push(p);
    }
    Ok(out)
}

struct KcMember {
    sys: fs_phs::PortHamiltonian,
    x: Vec<f64>,
    radiation_scale: f64,
    listener_m: f64,
}

fn kc_string_member(
    string: PrestressedString,
    pluck: Option<Pluck>,
    _bow: Option<BowStroke>,
    obstacles: &[UnilateralObstacle],
    gas: &GasState,
    listener_m: f64,
    radiation_scale: f64,
) -> Result<KcMember, AcousticRealizeError> {
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
    let sys = assemble_kc(string, storage, &zetas, obstacles)?;
    let mut x = vec![0.0; 2 * string.n_modes];
    if let Some(pluck) = pluck {
        for k in 1..=string.n_modes {
            x[2 * (k - 1)] = mass_scale * triangular_pluck_modal(pluck, k);
        }
    }
    Ok(KcMember {
        sys,
        x,
        radiation_scale,
        listener_m,
    })
}

fn step_kc_member(
    member: &mut KcMember,
    bow: Option<BowStroke>,
    texture: &mut TextureDrive,
    obstacles: &[UnilateralObstacle],
    string: PrestressedString,
    gas: &GasState,
    dt: f64,
) -> Result<(f64, Vec<f64>), AcousticRealizeError> {
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let pi = core::f64::consts::PI;
    let mut u = vec![0.0; string.n_modes];
    if let Some(bow) = bow {
        let v_string: f64 = (0..string.n_modes)
            .map(|k| {
                let phi = det::sin((k + 1) as f64 * pi * bow.station_frac) / mass_scale;
                member.x[2 * k + 1] * phi
            })
            .sum();
        let f_bow = bow_force(
            bow,
            v_string,
            texture.delta_n(bow.velocity_m_s - v_string, dt),
        );
        #[allow(clippy::needless_range_loop)] // the modal index spans u and shapes
        for k in 0..string.n_modes {
            let phi = det::sin((k + 1) as f64 * pi * bow.station_frac) / mass_scale;
            u[k] = f_bow * phi;
        }
    }
    if obstacles.iter().any(|o| o.mu_kinetic > 0.0) {
        let extra = modal_friction_forces(string, obstacles, &member.x)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        for (uk, fk) in u.iter_mut().zip(extra) {
            *uk += fk;
        }
    }
    if obstacles.iter().any(|o| o.internal_loss > 0.0) {
        let extra = modal_hunt_crossley_forces(string, obstacles, &member.x)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        for (uk, fk) in u.iter_mut().zip(extra) {
            *uk += fk;
        }
    }
    let rec = step(&member.sys, &member.x, &u, dt)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    member.x = rec.x;
    let mut p = 0.0;
    let mut q_phys = vec![0.0; string.n_modes];
    #[allow(clippy::needless_range_loop)] // the modal index spans state and outputs
    for k in 0..string.n_modes {
        let q = member.x[2 * k];
        q_phys[k] = q / mass_scale;
        let omega = string_mode_omega(string, k + 1);
        let h_im = string_mode_h_im(
            string,
            k + 1,
            omega,
            mass_scale,
            member.radiation_scale,
            gas.density,
            gas.sound_speed,
            member.listener_m,
        );
        p += h_im * omega * q;
    }
    Ok((p, q_phys))
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
        .map_or(0.0, |f| -f)
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
            .map_or(0.0, |r| inner.kn * r.combined_effective_height_m)
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
        WaveguideEnd::FlangedOpen => Termination::FlangedOpen,
    };
    refuse_open_nyquist(duct, gas, sample_rate_hz)?;
    if let Some((sections, taps)) = bore_spec(duct)
        && plates.vk.is_none()
    {
        return realize_reed_ode(
            &sections,
            duct_mouth(duct.termination),
            &taps,
            reed,
            plates,
            gas,
            listener_m,
            sample_rate_hz,
            n,
            duct.wall,
        );
    }
    realize_reed_bore(
        &physics,
        gas,
        reed,
        termination,
        plates,
        listener_m,
        sample_rate_hz,
        n,
        wall_pin(duct.wall).as_ref(),
    )
}

fn viscothermal_pin(gas: &GasState) -> ViscothermalPin {
    ViscothermalPin {
        dynamic_viscosity: gas.dynamic_viscosity,
        gamma: gas.gamma,
        prandtl: gas.prandtl,
        foster_branches: 3,
    }
}

fn wall_pin(wall: Option<LocallyReactingWall>) -> Option<WallPin> {
    wall.map(|w| WallPin {
        surface_density: w.surface_density_kg_m2,
        stiffness_per_area: w.stiffness_pa_per_m,
        resistance: w.resistance_pa_s_per_m,
    })
}

fn ode_bore_chain(
    sections: &[AcousticSection],
    mouth: Option<MouthFlange>,
    inlets: usize,
    taps: &[AcousticTap],
    gas: &GasState,
    wall: Option<LocallyReactingWall>,
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    let wall = wall_pin(wall);
    acoustic_chain_mouth_wall(
        sections,
        gas.density,
        gas.sound_speed,
        mouth,
        inlets,
        taps,
        Some(&viscothermal_pin(gas)),
        wall.as_ref(),
    )
    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
}

fn bore_spec(duct: &ViscothermalDuct) -> Option<(Vec<AcousticSection>, Vec<AcousticTap>)> {
    if duct.segments.is_empty()
        || duct
            .segments
            .iter()
            .any(|s| !(s.radius_m > 0.0 && s.length_m > 0.0 && s.outlet_radius_m > 0.0))
    {
        return None;
    }
    let total: f64 = duct.segments.iter().map(|s| s.length_m).sum();
    if !(total > 0.0) {
        return None;
    }
    let mut sections = Vec::new();
    for s in &duct.segments {
        let n = ((SECTION_BUDGET * s.length_m / total).round() as usize).max(2);
        sections.push(AcousticSection {
            length: s.length_m,
            radius: s.radius_m,
            outlet_radius: s.outlet_radius_m,
            cells: n,
        });
    }
    let mut taps = Vec::new();
    let mut acc = 0.0;
    for (i, s) in duct.segments.iter().enumerate() {
        acc += s.length_m;
        for hole in &duct.tone_holes {
            if hole.after_segment == i && hole.radius_m > 0.0 {
                taps.push(AcousticTap {
                    station: (acc / total).clamp(0.0, 1.0),
                    neck_length: hole.chimney_m.max(0.0),
                    neck_radius: hole.radius_m,
                    open_fraction: hole.open_fraction.clamp(0.0, 1.0),
                });
            }
        }
    }
    Some((sections, taps))
}

#[allow(clippy::too_many_arguments)]
fn realize_blown_ode(
    sections: &[AcousticSection],
    mouth: Option<MouthFlange>,
    taps: &[AcousticTap],
    blow: VolumeVelocityPulse,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<LocallyReactingWall>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let radius = sections.first().map_or(0.0, |s| s.radius);
    let area = core::f64::consts::PI * radius * radius;
    let inlets = if plates.linear.is_empty() { 1 } else { 2 };
    let line = ode_bore_chain(sections, mouth, inlets, taps, gas, wall)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    if plates.linear.is_empty() {
        let mut x = vec![0.0; line.state_dim()];
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 * dt;
            let u = if t < blow.duration_s {
                blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
            } else {
                0.0
            };
            let rec = step(&line, &x, &[u], dt)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            out.push(rec.y[0]);
            x = rec.x;
        }
        return Ok(out);
    }
    let plate =
        plate_modal_ports(plates, false)?.ok_or(AcousticRealizeError::InvalidDescription {
            what: "ODE duct expected a linear plate bank",
        })?;
    let sys = transformer(plate, line, 0, 1, area)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let mut x = vec![0.0; sys.state_dim()];
    let mut out = Vec::with_capacity(n);
    let n_plate = 2 * plates.linear.len();
    for i in 0..n {
        let t = i as f64 * dt;
        let u = if t < blow.duration_s {
            blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
        } else {
            0.0
        };
        let rec =
            step(&sys, &x, &[u], dt).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let mut p = rec.y[0];
        for (k, body) in plates.linear.iter().enumerate() {
            let acc = (rec.x[2 * k + 1] - x[2 * k + 1]) / dt;
            p += body.radiate(acc, gas.density, listener_m);
        }
        debug_assert!(rec.x.len() >= n_plate);
        out.push(p);
        x = rec.x;
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)] // one coherent realization
fn realize_reed_ode(
    sections: &[AcousticSection],
    mouth: Option<MouthFlange>,
    taps: &[AcousticTap],
    reed: BeatingReed,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<LocallyReactingWall>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    if !(reed.rest_opening_m > 0.0
        && reed.width_m > 0.0
        && reed.closing_pressure_pa > 0.0
        && reed.blowing_pressure_pa >= 0.0
        && reed.mass_kg >= 0.0
        && reed.mass_kg.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "reed parameters must be physical and finite",
        });
    }
    let radius = sections.first().map_or(0.0, |s| s.radius);
    let area = core::f64::consts::PI * radius * radius;
    let inlets = if plates.linear.is_empty() { 1 } else { 2 };
    let line = ode_bore_chain(sections, mouth, inlets, taps, gas, wall)?;
    let n_line = line.state_dim();
    let (line, joined) = if plates.linear.is_empty() {
        (Some(line), None)
    } else {
        let plate =
            plate_modal_ports(plates, false)?.ok_or(AcousticRealizeError::InvalidDescription {
                what: "ODE reed expected a linear plate bank",
            })?;
        (
            None,
            Some(
                transformer(plate, line, 0, 1, area)
                    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
            ),
        )
    };
    let massive = if reed.mass_kg > 0.0 {
        let (k, r_damp) = reed_structural(reed);
        let phs = mass_spring_damper(reed.mass_kg, k, r_damp)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let face = reed.width_m * 0.025;
        let k_lay = 1.0e7 * reed.width_m;
        let chi = r_damp / (k_lay * reed.rest_opening_m * reed.rest_opening_m).max(1.0e-18);
        let lay = slit_lay(k_lay, 2.0)
            .and_then(|o| o.with_internal_loss(chi))
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        Some((phs, face, lay))
    } else {
        None
    };
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut x = vec![
        0.0;
        joined
            .as_ref()
            .map_or(n_line, fs_phs::PortHamiltonian::state_dim)
    ];
    let mut x_reed = vec![0.0; 2];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 * dt;
        let p_m = blowing_envelope(reed, t);
        let p_bore = if let Some(sys) = &joined {
            sys.output(&x).first().copied().unwrap_or(0.0)
        } else {
            line.as_ref()
                .map_or(0.0, |s| s.output(&x).first().copied().unwrap_or(0.0))
        };
        let u = if let Some((reed_phs, face, lay)) = &massive {
            let q = x_reed[0];
            let opening = (reed.rest_opening_m + q).max(0.0);
            let v = reed_phs.output(&x_reed).first().copied().unwrap_or(0.0);
            let dp = p_m - p_bore;
            let u_jet = bernoulli_volume_flow(reed.width_m, opening, dp, gas.density);
            let y = reed.rest_opening_m + q;
            let mut f = -face * dp;
            let contact = slit_contact_force(lay, y)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            f += contact;
            f += lay.dissipative_modal_forces(1, &[y, v], &[v])[0];
            let rec_r = step(reed_phs, &x_reed, &[f], dt)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            if !rec_r.x.iter().all(|xi| xi.is_finite()) {
                return Err(AcousticRealizeError::Reed {
                    what: "massive reed left the finite set",
                });
            }
            x_reed = rec_r.x;
            u_jet + face * v
        } else {
            let h = quasistatic_aperture_opening(
                reed.rest_opening_m,
                reed.closing_pressure_pa,
                p_m - p_bore,
            );
            bernoulli_volume_flow(reed.width_m, h, p_m - p_bore, gas.density)
        };
        let rec = if let Some(sys) = &joined {
            step(sys, &x, &[u], dt)
        } else {
            step(line.as_ref().expect("waveguide-only reed"), &x, &[u], dt)
        }
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let mut p = rec.y[0];
        if !plates.linear.is_empty() {
            for (k, body) in plates.linear.iter().enumerate() {
                let acc = (rec.x[2 * k + 1] - x[2 * k + 1]) / dt;
                p += body.radiate(acc, gas.density, listener_m);
            }
        }
        if !p.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "ODE reed pressure left the finite set",
            });
        }
        x = rec.x;
        out.push(p);
    }
    Ok(out)
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
        WaveguideEnd::FlangedOpen => Termination::FlangedOpen,
    };
    refuse_open_nyquist(duct, gas, sample_rate_hz)?;
    if let Some((sections, taps)) = bore_spec(duct)
        && plates.vk.is_none()
    {
        return realize_blown_ode(
            &sections,
            duct_mouth(duct.termination),
            &taps,
            blow,
            plates,
            gas,
            listener_m,
            sample_rate_hz,
            n,
            duct.wall,
        );
    }
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
            wall_pin(duct.wall).as_ref(),
        );
    }
    // Isolated linear blow: the same DelayedFilter port as a reed or
    // a body, filled with IFFT[Z] rather than IFFT[R]. A vented
    // reflectance FIR does not ring a measurable period; the
    // impedance FIR does, so tone-hole shortening stays TMM-emergent.
    let mut line = crate::driving_point::impedance_line(
        &physics,
        gas,
        termination,
        sample_rate_hz,
        n,
        wall_pin(duct.wall).as_ref(),
    )
    .map_err(map_drive)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 * dt;
        let u = if t < blow.duration_s {
            blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
        } else {
            0.0
        };
        let p = line.push(u).map_err(|_| AcousticRealizeError::Reed {
            what: "impedance line left the finite set",
        })?;
        out.push(p);
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)] // one coherent realization record
#[allow(clippy::needless_range_loop)] // the sample index is the time axis
fn realize_blown_with_body(
    physics: &Duct,
    blow: VolumeVelocityPulse,
    plates: &mut PlateBank,
    gas: &GasState,
    termination: Termination,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<&WallPin>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    use crate::driving_point::characteristic_line;
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
    let mut line = characteristic_line(physics, gas, termination, sample_rate_hz, n, zc, wall)
        .map_err(map_drive)?;
    let mut out = vec![0.0; n];
    for i in 0..n {
        let t = i as f64 * dt;
        let u_blow = if t < blow.duration_s {
            blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
        } else {
            0.0
        };
        let u_body = plates.volume_velocity();
        let p_minus = line.incoming();
        let p_plus = p_minus + zc * (u_blow + u_body);
        let p_minus_now = line.push(p_plus).map_err(|_| AcousticRealizeError::Reed {
            what: "characteristic line left the finite set",
        })?;
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
        if s.is_taper() {
            segments.push(Segment::Cone {
                inlet_radius: s.radius_m,
                outlet_radius: s.outlet_radius_m,
                length: s.length_m,
            });
        } else {
            segments.push(Segment::Cylinder {
                radius: s.radius_m,
                length: s.length_m,
            });
        }
        for hole in &duct.tone_holes {
            if hole.after_segment == i {
                segments.push(Segment::ToneHole {
                    hole_radius: hole.radius_m,
                    chimney_height: hole.chimney_m,
                    bore_radius: s.radius_m,
                    state: hole_state(hole.open_fraction),
                });
            }
        }
    }
    Ok(Duct { segments })
}

fn hole_state(sigma: f64) -> HoleState {
    let s = sigma.clamp(0.0, 1.0);
    if s <= 0.0 {
        HoleState::Closed
    } else if s >= 1.0 {
        HoleState::Open
    } else {
        HoleState::Vent(s)
    }
}

fn duct_mouth(end: WaveguideEnd) -> Option<MouthFlange> {
    match end {
        WaveguideEnd::Closed => None,
        WaveguideEnd::UnflangedOpen => Some(MouthFlange::Unflanged),
        WaveguideEnd::FlangedOpen => Some(MouthFlange::Flanged),
    }
}

fn refuse_open_nyquist(
    duct: &ViscothermalDuct,
    gas: &GasState,
    sample_rate_hz: u32,
) -> Result<(), AcousticRealizeError> {
    if !matches!(duct.termination, WaveguideEnd::UnflangedOpen) {
        // Flanged mouths use the Rayleigh piston above ka = 1.
        return Ok(());
    }
    let mouth_r = duct
        .segments
        .last()
        .expect("non-empty checked above")
        .outlet_radius_m;
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
