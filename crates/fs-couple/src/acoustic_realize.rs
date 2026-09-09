//! Realize an [`fs_scenario::AcousticAssembly`] into observer pressure.
//!
//! A guitar or clarinet is a description: Kirchhoff–Carrier strings,
//! Rayleigh damping, a radiating plate, a beating reed on a viscothermal
//! TMM bore, optional tone holes, and a regularized bow. There is no
//! instrument crate.

use crate::modal_acoustic_time::{
    MAX_TIME_DOMAIN_ACOUSTIC_MODES, ModalAcousticFrame, ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeError, ModalAcousticTimeModel,
};
use crate::pcm_wav::{WavError, encode_pcm16_wav};
use crate::reed_bore::{blowing_envelope, realize_reed_bore, reed_pressure_face, reed_structural};
use crate::thin_plate::{
    PlateBank, PlateChartRadiation, VkBody, certified_chart_radiators, certified_radiators,
    vk_plate_phs,
};
use crate::unilateral_contact::{
    modal_contact_forces, modal_friction_forces, modal_hunt_crossley_forces, slit_contact_force,
    slit_lay,
};
/// Sections-per-metre budget for slicing a duct into lumped cells.
const SECTION_BUDGET: f64 = 8.0;

use fs_duct::{Duct, DuctError, HoleState, MAX_RADIATION_KA, Segment, Termination};
use fs_material::gas::GasState;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::visco::RayleighDamping;
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
    /// The execution scope refused work before candidate publication.
    Cancelled,
    /// Description failed a physical domain check.
    InvalidDescription {
        /// Which field.
        what: &'static str,
    },
    /// Gas-state derivation failed.
    Ambient(String),
    /// Modal time-stepper refusal.
    Modal(ModalAcousticTimeError),
    /// Plate mesh, section, or eigenproblem refusal, retaining the owner diagnostic.
    Plate(fs_plate::PlateError),
    /// A named plate-material assignment cannot be compiled.
    PlateAssignment {
        /// Region label, when the failure belongs to a particular assignment.
        region: Option<String>,
        /// Mesh triangle index, when applicable.
        element: Option<usize>,
        /// Failed assignment condition.
        what: &'static str,
    },
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
            Self::Cancelled => write!(f, "FS-COUPLE-ASSEMBLY-CANCELLED"),
            Self::InvalidDescription { what } => {
                write!(f, "FS-COUPLE-ASSEMBLY: {what}")
            }
            Self::Ambient(what) => write!(f, "FS-COUPLE-ASSEMBLY-GAS: {what}"),
            Self::Modal(e) => write!(f, "{e}"),
            Self::Plate(e) => write!(f, "{e}"),
            Self::PlateAssignment {
                region,
                element,
                what,
            } => write!(
                f,
                "FS-COUPLE-PLATE-ASSIGNMENT: region {region:?}, element {element:?}: {what}"
            ),
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
    realize_assembly_inner(assembly, None)
}

/// Realize an assembly with a linear plate chart, including regional materials
/// and an explicit force footprint. The chart replaces `assembly.plate`, which
/// must be absent. It feeds the existing body ports for strings, ducts, cavities
/// and the same pressure observer; extra authored lumped bodies still compose.
///
/// # Errors
/// Conflicting plate descriptions, invalid chart/controls, or the ordinary
/// assembly admission and stepping errors refuse.
pub fn realize_assembly_with_plate_chart(
    assembly: &AcousticAssembly,
    chart: &fs_plate::PlateChart,
    options: &PlateChartRadiation,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    if assembly.plate.is_some() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "chart assembly requires the uniform plate description to be absent",
        });
    }
    realize_assembly_inner(assembly, Some((chart, options)))
}

#[allow(clippy::too_many_lines)] // one coherent realization dispatcher
fn realize_assembly_inner(
    assembly: &AcousticAssembly,
    chart: Option<(&fs_plate::PlateChart, &PlateChartRadiation)>,
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
    if let Some(string) = assembly
        .string
        .as_ref()
        .filter(|s| s.relaxation_bending.is_some())
    {
        validate_relaxation_modes(string, assembly.sample_rate_hz)?;
    }
    let n = sample_count(assembly.sample_rate_hz, assembly.duration_s)?;
    let mut pressure_pa = vec![0.0; n];
    let mut bodies = plate_bank(assembly.soundboard, &assembly.body_modes, assembly.plate)?;
    if let Some((chart, options)) = chart {
        bodies
            .linear
            .extend(certified_chart_radiators(chart, options)?);
    }
    bodies.attach_radiation_loads(&gas, assembly.sample_rate_hz);
    let dirac_base = assembly.string.as_ref().is_some_and(|s| s.moving_end);
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
        let string = assembly.string.as_ref().expect("checked");
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
    if let Some(string) = &assembly.string {
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
    let string = assembly.string.as_ref().expect("checked");
    let twin = secondary_string(string, (1.0 + string.polarization_detune).powi(2));
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
    if string.axial_stiffness_n > 0.0 || string.relaxation_bending.is_some() {
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
            members.push(kc_string_member(
                &twin,
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
                let member_string = if idx == 0 { string } else { &twin };
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
        members.push(linear_string_member(
            &twin,
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
            let member_string = if idx == 0 { string } else { &twin };
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
    // Structural port slot: a shorter output vector is an internal
    // invariant break, never physics — refuse loudly instead of zeroing.
    body.output(x)[0]
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
        p += body.radiate_mass_normalized(acc, gas.density, listener_m);
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
    let string = assembly.string.as_ref().expect("checked");
    let twin = secondary_string(string, (1.0 + string.polarization_detune).powi(2));
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
    if string.axial_stiffness_n > 0.0 || string.relaxation_bending.is_some() {
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
        members.push(linear_string_member(
            &twin,
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
            let member_string = if idx == 0 { string } else { &twin };
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
    let string = assembly.string.as_ref().expect("checked");
    let twin = secondary_string(string, (1.0 + string.polarization_detune).powi(2));
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
        members.push(kc_string_member(
            &twin,
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
            let member_string = if idx == 0 { string } else { &twin };
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

pub(crate) fn map_drive(err: crate::driving_point::DrivingPointError) -> AcousticRealizeError {
    match err {
        crate::driving_point::DrivingPointError::Invalid { what } => {
            AcousticRealizeError::Reed { what }
        }
        crate::driving_point::DrivingPointError::Duct(d) => AcousticRealizeError::Duct(d),
        crate::driving_point::DrivingPointError::Realize(_) => AcousticRealizeError::Reed {
            what: "characteristic realization refused",
        },
        crate::driving_point::DrivingPointError::Discrete(e) => AcousticRealizeError::Nonlinear(
            format!("characteristic-line discretization refused: {e}"),
        ),
    }
}

fn gas_state(ambient: AmbientGas) -> Result<GasState, AcousticRealizeError> {
    if !(ambient.relative_humidity >= 0.0 && ambient.relative_humidity <= 1.0) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "relative humidity must be an explicit fraction in [0, 1]",
        });
    }
    GasState::try_new_moist_air(
        ambient.temperature_k,
        ambient.pressure_pa,
        ambient.relative_humidity,
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
pub fn string_mode_omega(string: &PrestressedString, k: usize) -> f64 {
    prestressed_beam_omega(
        string.length_m,
        string.tension_n,
        string.lin_density_kg_m,
        string.bending_stiffness_n_m2,
        k,
    )
}

/// Compact strip observer from actual modal accelerations, not `-omega² q`.
/// Expanding distributed retarded monopoles gives
/// `p = rho/(4 pi r) [integral(w a dx) + integral(w (x-L/2) jerk dx)/c]`.
/// See Wood, Euphonics 11.7.1, eq. (6):
/// <https://euphonics.org/11-7-1-compact-sound-sources-monopoles-dipoles-quadrupoles/>.
/// The first moment assumes an observer on the positive string axis; the
/// scenario has no directional listener. This remains a compact strip model,
/// not the radiation solution for a cylindrical wire. Common propagation
/// delay is omitted. Jerk uses backward differences of accepted accelerations
/// (first-order); instantaneous release impulses are not reconstructed.
#[derive(Clone)]
struct StringObserver {
    weights: Vec<(f64, f64)>,
    previous_acceleration: Option<Vec<f64>>,
}

impl StringObserver {
    fn new(string: &PrestressedString, gas: &GasState, range: f64, scale: f64) -> Self {
        let pi = core::f64::consts::PI;
        let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
        let factor = scale * gas.density * string.width_m / (4.0 * pi * range * mass_scale);
        let weights = (1..=string.n_modes)
            .map(|k| {
                let kpi = k as f64 * pi;
                if k % 2 == 1 {
                    (factor * 2.0 * string.length_m / kpi, 0.0)
                } else {
                    // Signed first moment of the sine about the string midpoint.
                    (
                        0.0,
                        -factor * string.length_m.powi(2) / (kpi * gas.sound_speed),
                    )
                }
            })
            .collect();
        Self {
            weights,
            previous_acceleration: None,
        }
    }

    fn observe(
        &mut self,
        acceleration: Vec<f64>,
        initial_acceleration: Option<Vec<f64>>,
        dt: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let previous = self
            .previous_acceleration
            .as_ref()
            .or(initial_acceleration.as_ref())
            .expect("first step supplies its initial acceleration");
        let pressure: f64 = self
            .weights
            .iter()
            .zip(&acceleration)
            .zip(previous)
            .map(|((&(area, moment), &a), &before)| area * a + moment * ((a - before) / dt))
            .sum();
        let limit = ModalAcousticTimeBudget::audible_reference().maximum_abs_pressure_pa;
        if !pressure.is_finite() || pressure.abs() > limit {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "compact string observer pressure is nonfinite or exceeds its pressure budget",
            });
        }
        self.previous_acceleration = Some(acceleration);
        Ok(pressure)
    }
}

fn assemble_kc(
    string: &PrestressedString,
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
    string: &PrestressedString,
    omega: f64,
    wave_number: f64,
    gas: &GasState,
) -> Result<f64, AcousticRealizeError> {
    Ok(mode_loss_parts(string, omega, wave_number, gas)?.total())
}

#[derive(Clone, Copy, Default)]
struct StringModeLossParts {
    material: f64,
    air: f64,
    authored: f64,
}

impl StringModeLossParts {
    fn total(self) -> f64 {
        self.material + self.authored + self.air
    }
}

fn mode_loss_parts(
    string: &PrestressedString,
    omega: f64,
    wave_number: f64,
    gas: &GasState,
) -> Result<StringModeLossParts, AcousticRealizeError> {
    if string.relaxation_bending.is_some()
        && (string.kelvin_voigt_bending.is_some()
            || string.rayleigh.is_some()
            || string.damping_ratio != 0.0)
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "relaxation bending cannot be combined with other internal/Rayleigh losses",
        });
    }
    let material_bending = string
        .kelvin_voigt_bending
        .map(|law| {
            let (lo, hi) = law.omega_band_rad_s;
            if string.rayleigh.is_some()
                || string.damping_ratio != 0.0
                || !law.viscous_stiffness_n_m2_s.is_finite()
                || law.viscous_stiffness_n_m2_s < 0.0
                || !(lo.is_finite() && hi.is_finite() && lo >= 0.0 && hi >= lo)
                || !(omega.is_finite() && omega > 0.0 && omega >= lo && omega <= hi)
            {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "Kelvin-Voigt bending needs nonnegative viscosity, an admitted modal frequency, and no authored internal/Rayleigh loss",
                });
            }
            // Galerkin projection of (eta I y_xxt)_xx gives c_k = eta I k^4.
            // With mu q_tt + c_k q_t + (T k² + EI k⁴) q = 0,
            // zeta_k = c_k/(2 mu omega_k). Tension energy is lossless: using
            // eta*omega/(2E) here would omit the bending-energy dilution.
            // Sakthivel et al. (2023), https://arxiv.org/abs/2301.07931.
            let zeta = 0.5 * law.viscous_stiffness_n_m2_s
                * wave_number.powi(4)
                / string.lin_density_kg_m
                / omega;
            if !zeta.is_finite() || zeta < 0.0 {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "Kelvin-Voigt bending damping is unrepresentable",
                });
            }
            Ok(zeta)
        })
        .transpose()?;
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
        return Ok(StringModeLossParts {
            authored: r.zeta_at(omega),
            ..StringModeLossParts::default()
        });
    }
    // The compact string path interprets width as circular diameter for air
    // drag. R [N s/m²] / (2 mu_linear omega) is dimensionless. Material
    // bindings derive this diameter and mu_linear from the same cross-section.
    let resistance = crate::air_path::oscillating_cylinder_air_resistance_per_length(
        0.5 * string.width_m,
        omega,
        gas,
    )?;
    let stokes = 0.5 * resistance / string.lin_density_kg_m / omega;
    if !stokes.is_finite() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "cylinder air damping ratio is unrepresentable",
        });
    }
    if let Some(bending) = material_bending {
        return Ok(StringModeLossParts {
            material: bending,
            air: stokes,
            authored: 0.0,
        });
    }
    if string.relaxation_bending.is_some() {
        // The memory arms supply their own storage and dissipation in pHS.
        return Ok(StringModeLossParts {
            air: stokes,
            ..StringModeLossParts::default()
        });
    }
    // A caller-authored modal ratio is an explicit reduced model. Elastic EI
    // supplies storage, not an additional dissipation law, and one ratio does
    // not identify a relaxation spectrum. Physical bending loss must select
    // the Kelvin-Voigt or causal memory path above.
    Ok(StringModeLossParts {
        authored: string.damping_ratio,
        air: stokes,
        material: 0.0,
    })
}

fn relaxing_mode_stiffness(string: &PrestressedString, stiffness_n_m2: f64, k: usize) -> f64 {
    let order = k as f64 + if string.moving_end { 0.5 } else { 1.0 };
    let wave_number = order * core::f64::consts::PI / string.length_m;
    stiffness_n_m2 * wave_number.powi(4) / string.lin_density_kg_m
}

fn validate_relaxation_modes(
    string: &PrestressedString,
    sample_rate: u32,
) -> Result<(), AcousticRealizeError> {
    let law = string
        .relaxation_bending
        .as_ref()
        .expect("selected relaxation law");
    let refuse = || AcousticRealizeError::InvalidDescription {
        what: "relaxation bending needs finite coefficients, exclusive losses, and admitted relaxed/instantaneous modal frequencies below Nyquist",
    };
    let (lo, hi) = law.omega_band_rad_s;
    if !(law.branches.iter().all(|branch| {
        branch.relaxing_stiffness_n_m2.is_finite()
            && branch.relaxing_stiffness_n_m2 >= 0.0
            && branch.relaxation_time_s.is_finite()
            && branch.relaxation_time_s > 0.0
    }) && lo.is_finite()
        && hi.is_finite()
        && lo >= 0.0
        && hi >= lo)
        || string.kelvin_voigt_bending.is_some()
        || string.rayleigh.is_some()
        || string.damping_ratio != 0.0
    {
        return Err(refuse());
    }
    if string.moving_end && string.polarization_detune != 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "moving-end relaxation does not support a second polarization",
        });
    }
    // The fixed sample clock uses the midpoint discrete gradient. For an
    // isolated pole, h <= 0.1 keeps both phase-rate and relaxation-rate error
    // below 0.1%; Nyquist alone admits severely warped/stiff transients.
    // This reference check is not a bound on nonlinear or coupled spectra.
    let reference_rate_limit = 0.1 * f64::from(sample_rate);
    if law.branches.iter().any(|branch| {
        branch.relaxing_stiffness_n_m2 > 0.0
            && 1.0 / branch.relaxation_time_s > reference_rate_limit
    }) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "sample rate does not resolve bending relaxation: require dt/tau <= 0.1; increase sample rate",
        });
    }
    let guard = ModalAcousticTimeBudget::audible_reference().nyquist_guard_fraction
        * core::f64::consts::PI
        * f64::from(sample_rate);
    let mut member = string.clone();
    let total_stiffness = law.branches.iter().map(|b| b.relaxing_stiffness_n_m2).sum();
    for polarization in 0..=usize::from(string.polarization_detune > 0.0) {
        if polarization == 1 {
            member.tension_n *= (1.0 + string.polarization_detune).powi(2);
        }
        for k in 0..string.n_modes {
            let relaxed = if string.moving_end {
                moving_end_omega(&member, k)
            } else {
                string_mode_omega(&member, k + 1)
            };
            let increment = relaxing_mode_stiffness(&member, total_stiffness, k);
            let instant = det::sqrt(relaxed * relaxed + increment);
            if !(relaxed.is_finite()
                && relaxed > 0.0
                && relaxed >= lo
                && increment.is_finite()
                && increment >= 0.0
                && instant.is_finite()
                && instant <= hi
                && instant <= guard)
            {
                return Err(refuse());
            }
            if instant > reference_rate_limit {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "sample rate does not resolve bending reference phase: require dt*omega_instant <= 0.1; increase sample rate",
                });
            }
        }
    }
    Ok(())
}

fn with_bending_relaxation(
    string: &PrestressedString,
    sys: fs_phs::PortHamiltonian,
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    let Some(law) = &string.relaxation_bending else {
        return Ok(sys);
    };
    let mut branches = Vec::new();
    // Branch-major, then mode-major: the initializer uses this same order.
    for branch in law
        .branches
        .iter()
        .filter(|b| b.relaxing_stiffness_n_m2 > 0.0)
    {
        for k in 0..string.n_modes {
            let mut projection = vec![0.0; sys.state_dim()];
            projection[2 * k] = 1.0;
            branches.push(fs_phs::RelaxationBranch {
                projection,
                stiffness: relaxing_mode_stiffness(string, branch.relaxing_stiffness_n_m2, k),
                relaxation_time_s: branch.relaxation_time_s,
            });
        }
    }
    sys.with_relaxation_branches(branches)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
}

fn initialize_relaxed_bending(string: &PrestressedString, x: &mut [f64]) {
    let Some(law) = &string.relaxation_bending else {
        return;
    };
    for (j, branch) in law
        .branches
        .iter()
        .filter(|b| b.relaxing_stiffness_n_m2 > 0.0)
        .enumerate()
    {
        for k in 0..string.n_modes {
            x[(2 + j) * string.n_modes + k] = det::sqrt(relaxing_mode_stiffness(
                string,
                branch.relaxing_stiffness_n_m2,
                k,
            )) * x[2 * k];
        }
    }
}

fn secondary_string(string: &PrestressedString, tension_scale: f64) -> PrestressedString {
    let mut twin = string.clone();
    twin.polarization_detune = 0.0;
    twin.tension_n *= tension_scale;
    twin
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
    string: &PrestressedString,
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
            let wave_number = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
            mode_zeta(string, omega, wave_number, gas)
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
    initialize_relaxed_bending(string, &mut x);
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
                p += body.radiate_mass_normalized(acc, gas.density, listener_m);
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
    string: &PrestressedString,
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
    let base = if string.axial_stiffness_n > 0.0 {
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
    }?;
    with_bending_relaxation(string, base)
}

#[allow(clippy::too_many_arguments)] // one coherent solver probe
#[allow(clippy::unnecessary_wraps)] // uniform Result surface with the ODE twin
fn leftover_u(
    sys: &fs_phs::DescriptorPortHamiltonian,
    x: &[f64],
    n_s: usize,
    string: &PrestressedString,
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
        u[col] = bow_force(bow, v, texture.delta_n(bow.velocity_m_s - v, dt))
            .map_err(AcousticRealizeError::Nonlinear)?;
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
    string: &PrestressedString,
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
        u[col] = bow_force(bow, v, texture.delta_n(bow.velocity_m_s - v, dt))
            .map_err(AcousticRealizeError::Nonlinear)?;
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

fn moving_end_opening(string: &PrestressedString, x: &[f64], n_s: usize, station: f64) -> f64 {
    let phi0 = det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m));
    let mut y = 0.0;
    for k in 0..n_s {
        let kappa = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
        let phi = phi0 * det::cos(kappa * station.clamp(0.0, 1.0) * string.length_m);
        y += x.get(2 * k).copied().unwrap_or(0.0) * phi;
    }
    y
}

fn moving_end_omega(string: &PrestressedString, k_zero: usize) -> f64 {
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

fn free_fixed_pluck_modal(pluck: Pluck, string: &PrestressedString, k_zero: usize) -> f64 {
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
        .map(|b| b.drive_participation / b.mass_kg.sqrt())
        .collect();
    let sys = if with_area_port {
        let areas: Vec<f64> = plates
            .linear
            .iter()
            .map(|b| b.area_m2 / b.mass_kg.sqrt())
            .collect();
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
    string: &PrestressedString,
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
    if string.axial_stiffness_n > 0.0 || string.relaxation_bending.is_some() {
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
    string: &PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
) -> Result<(), AcousticRealizeError> {
    if !(string.damping_ratio.is_finite() && string.damping_ratio >= 0.0) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "authored string modal damping ratio must be finite and nonnegative",
        });
    }
    if !(string.length_m > 0.0
        && string.tension_n > 0.0
        && string.lin_density_kg_m > 0.0
        && string.axial_stiffness_n >= 0.0
        && string.width_m > 0.0
        && string.n_modes > 0
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
            && bow.velocity_m_s.is_finite()
            && bow.normal_force_n.is_finite()
            && bow.stribeck_m_s > 0.0
            && bow.mu_static.is_finite()
            && bow.mu_dynamic.is_finite()
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
    string: &PrestressedString,
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
    let twin = secondary_string(string, twin_tension.unwrap_or(1.0));
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
    if twin_tension.is_some() {
        members.push(linear_string_member(
            &twin,
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
            let member_string = if idx == 0 { string } else { &twin };
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

#[derive(Clone)]
struct LinearMember {
    model: ModalAcousticTimeModel,
    phi_bow: Vec<f64>,
    observer: StringObserver,
    loss_parts: Vec<StringModeLossParts>,
}

/// Incremental small-amplitude vibration of one source-bound, pinned string.
///
/// This uses the same exact-ZOH modes, loss admission and compact strip observer
/// as the assembly realizer. The explicit linearization omits axial stretching
/// induced by transverse motion. Ambient gas and listener distance stay fixed;
/// no atmospheric path filter, common propagation delay or radiation loading is
/// included. The observer is a diagnostic, not an additional energy debit.
///
/// Rebinding keeps the exact length, linear density and sine-mode count, hence
/// the same mass-normalized basis. Changing geometry/supports, Prony memory,
/// nonlinear response or phase requires a different state-transfer owner.
/// Material temperature is supplied by the specimen, never copied from air.
/// The caller owns heat transport and must account for returned parameter work;
/// this runtime alone does not implement a closed thermomechanical system.
#[derive(Clone)]
pub struct LinearMaterialStringRuntime {
    specimen: crate::string_specimen::ResolvedStringSpecimen,
    member: LinearMember,
    gas: GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    epoch: u64,
    accepted_samples: u64,
}

/// Energy change when replacing stiffness at fixed modal displacement/momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringMaterialUpdate {
    /// Shared material/operator/observer epoch after publication.
    pub epoch: u64,
    /// Modal mechanical energy before the update [J].
    pub energy_before_j: f64,
    /// Modal mechanical energy after the update [J].
    pub energy_after_j: f64,
    /// Signed work into the retained vibration: `sum((omega_new²-omega_old²) q²/2)`.
    /// The thermal or external parameter owner owes the opposite energy transfer.
    /// This excludes the axial prestress reservoir and dissipated heat.
    pub parameter_work_j: f64,
}

/// One accepted vibration/pressure step and its physical loss destinations.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialStringFrame {
    /// Epoch of mechanics, observer and loss publication.
    pub epoch: u64,
    /// The underlying exact-ZOH acoustic frame, including work and modal energy.
    pub acoustic: ModalAcousticFrame,
    /// Energy removed by the separately admitted damping mechanisms.
    pub dissipation: StringDissipation,
}

/// Losses during one sample, excluding parameter work and diagnostic radiation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringDissipation {
    /// Kelvin-Voigt bending dissipation inside the solid [J].
    pub material_heat_j: f64,
    /// Energy removed by the admitted oscillating-cylinder air resistance [J].
    /// It leaves the solid; this reduced law does not resolve its fluid heating.
    pub air_loss_j: f64,
    /// Authored modal/Rayleigh damping with no identified thermal destination [J].
    pub authored_loss_j: f64,
    /// Aggregate work-minus-energy loss minus the three assigned channels [J].
    /// Includes negative mode-local roundoff and undamped numerical drift;
    /// it is never deposited as material heat.
    pub roundoff_residual_j: f64,
    /// Scale-aware allowance for the residual and numerical loss estimates [J].
    /// This is not a constitutive-model or temporal-discretization error bound.
    pub roundoff_tolerance_j: f64,
}

fn partition_string_dissipation(
    frame: &ModalAcousticFrame,
    parts: &[StringModeLossParts],
) -> Result<StringDissipation, AcousticRealizeError> {
    let mut loss = StringDissipation {
        material_heat_j: 0.0,
        air_loss_j: 0.0,
        authored_loss_j: 0.0,
        roundoff_residual_j: 0.0,
        roundoff_tolerance_j: frame.dissipation_roundoff_tolerance_j,
    };
    for (mode, parts) in frame.modal_viscous_losses.iter().zip(parts) {
        loss.roundoff_tolerance_j += mode.roundoff_tolerance_j;
        let total = parts.total();
        if total > 0.0 {
            // With constant c_i during the step, each mechanism removes
            // c_i integral(v^2 dt), so its share is c_i / sum(c_i).
            // Only admitted negative roundoff is excluded, with its signed
            // difference retained below; the mechanical state is untouched.
            let dissipated = mode.energy_j.max(0.0);
            loss.material_heat_j += dissipated * (parts.material / total);
            loss.air_loss_j += dissipated * (parts.air / total);
            loss.authored_loss_j += dissipated * (parts.authored / total);
        }
    }
    let assigned = loss.material_heat_j + loss.air_loss_j + loss.authored_loss_j;
    loss.roundoff_residual_j = frame.viscous_dissipation_j - assigned;
    loss.roundoff_tolerance_j += 8.0 * f64::EPSILON * parts.len() as f64 * assigned;
    if !loss.roundoff_residual_j.is_finite()
        || !loss.roundoff_tolerance_j.is_finite()
        || loss.roundoff_residual_j.abs() > loss.roundoff_tolerance_j
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "string damping partition does not close within its roundoff allowance",
        });
    }
    Ok(loss)
}

impl LinearMaterialStringRuntime {
    /// Initialize zero vibration or a released triangular pluck.
    ///
    /// All admission and cancellation checks precede returning the runtime.
    /// The retained model is bounded by the ordinary audible-reference budgets.
    pub fn try_new(
        cx: &fs_exec::Cx<'_>,
        specimen: crate::string_specimen::ResolvedStringSpecimen,
        pluck: Option<Pluck>,
        ambient: AmbientGas,
        listener_m: f64,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        let string = specimen.string();
        validate_incremental_string(&string)?;
        validate_string(&string, pluck, None)?;
        if !(listener_m.is_finite() && listener_m > 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "incremental string listener distance must be finite and positive",
            });
        }
        let gas = gas_state(ambient)?;
        let member =
            linear_string_member(&string, pluck, None, &gas, listener_m, sample_rate_hz, 1.0)?;
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        Ok(Self {
            specimen,
            member,
            gas,
            listener_m,
            sample_rate_hz,
            epoch: 0,
            accepted_samples: 0,
        })
    }

    /// Source receipts and physical coefficients at the accepted epoch.
    #[must_use]
    pub const fn specimen(&self) -> &crate::string_specimen::ResolvedStringSpecimen {
        &self.specimen
    }

    /// One epoch covers the specimen, operators, vibration and observer history.
    #[must_use]
    pub const fn epoch(&self) -> u64 {
        self.epoch
    }

    /// Number of accepted audio steps; rebinding never consumes a sample.
    #[must_use]
    pub const fn accepted_samples(&self) -> u64 {
        self.accepted_samples
    }

    /// Current states in the unchanged mass-normalized sine basis.
    #[must_use]
    pub fn states(&self) -> &[ModalAcousticState] {
        self.member.model.states()
    }

    /// Frequencies and complete admitted damping at the same epoch.
    #[must_use]
    pub fn modes(&self) -> &[ModalAcousticMode] {
        self.member.model.modes()
    }

    /// Advance one sample under held modal forces [N/sqrt(kg)].
    ///
    /// Mechanics, accepted acceleration history and sample position publish
    /// together after pressure/energy limits and a final cancellation poll.
    /// A refusal preserves the entire runtime for retry. The returned pressure
    /// is the compact strip observation; work and loss refer to retained modes.
    pub fn step(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        forces: &[f64],
    ) -> Result<MaterialStringFrame, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        let epoch = next_string_epoch(self.epoch)?;
        let samples = next_string_epoch(self.accepted_samples)?;
        let mut candidate = self.member.clone();
        let initial = candidate
            .observer
            .previous_acceleration
            .is_none()
            .then(|| linear_accelerations(&candidate.model, forces));
        let mut frame = candidate
            .model
            .step(forces)
            .map_err(AcousticRealizeError::Modal)?;
        frame.observer_pressure_pa = candidate.observer.observe(
            linear_accelerations(&candidate.model, forces),
            initial,
            candidate.model.sample_period_s(),
        )?;
        let dissipation = partition_string_dissipation(&frame, &candidate.loss_parts)?;
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        self.member = candidate;
        self.epoch = epoch;
        self.accepted_samples = samples;
        Ok(MaterialStringFrame {
            epoch,
            acoustic: frame,
            dissipation,
        })
    }

    /// Replace the material state without resetting vibration or its observer.
    ///
    /// The exact card, length, linear density and mode count must be unchanged.
    /// No mass tolerance or mode-number matching hides a basis transfer. Loss,
    /// tension, EI and observer weights rebuild together using existing owners.
    /// Previous accepted acceleration is retained for the observer's backward
    /// difference; an abrupt parameter jump can therefore emit a transient.
    /// Updates are piecewise constant, with no crossfade or adiabatic claim.
    ///
    /// `maximum_abs_parameter_work_j` is a caller-owned energy-transfer budget,
    /// not a rescaling target. A refusal leaves all state/receipts unchanged.
    pub fn rebind(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        specimen: crate::string_specimen::ResolvedStringSpecimen,
        maximum_abs_parameter_work_j: f64,
    ) -> Result<StringMaterialUpdate, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        let epoch = next_string_epoch(self.epoch)?;
        let old = self.specimen.string();
        let new = specimen.string();
        validate_incremental_string(&new)?;
        if self.specimen.material().card_identity() != specimen.material().card_identity()
            || old.length_m.to_bits() != new.length_m.to_bits()
            || old.lin_density_kg_m.to_bits() != new.lin_density_kg_m.to_bits()
            || old.n_modes != new.n_modes
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "incremental material update requires the same card, length, linear density and sine basis; use an explicit state transfer for other changes",
            });
        }
        if !(maximum_abs_parameter_work_j.is_finite() && maximum_abs_parameter_work_j >= 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "parameter-work budget must be finite and nonnegative",
            });
        }
        let mut candidate = linear_string_member(
            &new,
            None,
            None,
            &self.gas,
            self.listener_m,
            self.sample_rate_hz,
            1.0,
        )?;
        candidate
            .model
            .restore_states(self.states())
            .map_err(AcousticRealizeError::Modal)?;
        candidate
            .observer
            .previous_acceleration
            .clone_from(&self.member.observer.previous_acceleration);
        let energy = |model: &ModalAcousticTimeModel| -> f64 {
            model
                .modes()
                .iter()
                .zip(model.states())
                .map(|(mode, state)| {
                    0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                        + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2))
                })
                .sum()
        };
        let energy_before_j = energy(&self.member.model);
        let energy_after_j = energy(&candidate.model);
        // Difference the stiffness terms mode by mode; do not subtract large
        // kinetic energies when the parameter work is small.
        let parameter_work_j: f64 = self
            .modes()
            .iter()
            .zip(candidate.model.modes())
            .zip(self.states())
            .map(|((before, after), state)| {
                let old_wq = before.angular_frequency_rad_s * state.displacement_m_sqrt_kg;
                let new_wq = after.angular_frequency_rad_s * state.displacement_m_sqrt_kg;
                0.5 * (new_wq - old_wq) * (new_wq + old_wq)
            })
            .sum();
        if !parameter_work_j.is_finite() || parameter_work_j.abs() > maximum_abs_parameter_work_j {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "material update exceeds its parameter-work budget",
            });
        }
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        self.specimen = specimen;
        self.member = candidate;
        self.epoch = epoch;
        Ok(StringMaterialUpdate {
            epoch,
            energy_before_j,
            energy_after_j,
            parameter_work_j,
        })
    }
}

/// Coupled refusal; no mechanical, thermal or observer state is published.
#[derive(Debug)]
pub enum ThermalStringError {
    /// Unsupported physical model or invalid coupling budget.
    Admission(&'static str),
    /// The caller has already consumed this epoch, or is using stale state.
    Epoch {
        /// Epoch the caller intended to advance.
        expected: u64,
        /// Currently accepted epoch.
        actual: u64,
    },
    /// Mechanical or acoustic owner refused the candidate.
    Acoustic(AcousticRealizeError),
    /// Thermal owner refused the candidate.
    Phase(fs_material::phase::PhaseStateError),
    /// The supplied thermal transport refused, retaining its diagnostic.
    Transport(String),
    /// The material-backed uniform DC conductor refused the candidate.
    Electrical(fs_material::conductor::ConductorError),
    /// Material query coordinates were invalid.
    Query(fs_matdb::MatDbError),
    /// The pinned material claims cannot supply the candidate state.
    Material(fs_material::state_point::MaterialStatePointError),
}

impl core::fmt::Display for ThermalStringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-THERMAL-STRING: {self:?}")
    }
}

impl std::error::Error for ThermalStringError {}

/// Uniform solid temperature coupled to material-resolved viscous vibration.
///
/// Density and Young's modulus must be validity-wide constants, with prescribed
/// radius or mass and fixed tension: no changing elastic/prestress potential or
/// expansion work is hidden in heat. Viscosity may depend on temperature.
/// The admitted thermal chart has constant solid density and uses total
/// specific enthalpy; neither latent heat nor air loss is credited twice.
///
/// Each sample freezes viscosity at its initial temperature, evolves the exact
/// mechanical substep, deposits its material loss, and resolves viscosity at
/// the new temperature. This is first-order partitioned feedback, not an exact
/// continuous-temperature solution. The explicit temperature-increment budget
/// limits a step; time refinement is still required to assess coupling error.
/// Spatial transport, thermal expansion, melting and acoustic radiation loading are
/// outside this rung. Ambient gas remains independent of the solid temperature.
#[derive(Clone)]
pub struct ThermalMaterialStringRuntime {
    mechanical: LinearMaterialStringRuntime,
    card: fs_matdb::MaterialCard,
    point: fs_matdb::QueryPoint,
    requirements: Vec<ScalarPropertyRequirement>,
    selection: MaterialPropertySelection,
    curve: EquilibriumEnthalpyPhaseCurve,
    thermal: EquilibriumPhaseState,
    maximum_temperature_increment_k: f64,
}

/// One accepted thermal/mechanical sample, with a single publication epoch.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalMaterialStringFrame {
    /// Mechanical work/loss over the sample and final-state diagnostic pressure.
    pub vibration: MaterialStringFrame,
    /// Thermal state at the same epoch as `vibration.epoch`.
    pub thermal: EquilibriumPhaseState,
    /// Signed externally supplied heat [J]; its source lies outside this owner.
    pub external_heat_j: f64,
    /// Actual increase in total thermal enthalpy [J].
    pub enthalpy_change_j: f64,
    /// Work + external heat - mechanical/thermal storage changes - outgoing
    /// air/authored loss - the separately reported mechanical numerical loss [J].
    pub energy_balance_residual_j: f64,
    /// Arithmetic allowance [J], including the absolute enthalpy reference.
    /// This does not bound coupling truncation or constitutive-model error.
    pub energy_roundoff_tolerance_j: f64,
    /// Separately declared thermal-solve residual allowance [J], zero for
    /// directly supplied heat. Neither this nor roundoff bounds timestep error.
    pub thermal_solve_tolerance_j: f64,
}

/// Atomic vibration/temperature result and the actual thermal solver trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalStringTransportFrame<R> {
    /// Coupled state and first-law balance at the newly accepted epoch.
    pub coupled: ThermalMaterialStringFrame,
    /// Concrete report returned by the transport owner, without type erasure.
    pub transport: R,
}

/// One accepted ideal-source sample of a uniform conducting string.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectrothermalStringFrame<R = ()> {
    /// Accepted mechanics, temperature and pressure at one publication epoch.
    pub coupled: ThermalMaterialStringFrame,
    /// Authored ideal-source control retained independently of resolved current.
    pub drive: fs_material::conductor::OhmicDrive,
    /// Signed current evaluated with initial resistance and frozen for the sample [A].
    pub current_a: f64,
    /// Current demanded by the same source at the accepted final resistance [A].
    pub final_current_a: f64,
    /// Electrical source work converted to internal Joule heat [J].
    pub joule_heat_j: f64,
    /// Signed thermal boundary heat [J], excluding Joule and vibration heat.
    pub boundary_heat_j: f64,
    /// Actual thermal transport report, or `()` for insulated boundaries.
    pub transport: R,
    /// Source-resolved resistance at the beginning of the sample.
    /// This frozen law supplies the sample's Joule heat.
    pub conductor_before: fs_material::conductor::ResolvedConductor,
    /// Source-resolved resistance at the accepted final temperature.
    pub conductor_after: fs_material::conductor::ResolvedConductor,
}

impl ThermalMaterialStringRuntime {
    /// Attach the same immutable card and thermal chart to an accepted string.
    /// Existing material coordinates and selected claims are retained exactly.
    pub fn try_new(
        cx: &fs_exec::Cx<'_>,
        mechanical: LinearMaterialStringRuntime,
        card: fs_matdb::MaterialCard,
        curve: EquilibriumEnthalpyPhaseCurve,
        specific_enthalpy_j_kg: f64,
        maximum_temperature_increment_k: f64,
    ) -> Result<Self, ThermalStringError> {
        use crate::string_specimen::StringPrestress;
        let refuse = ThermalStringError::Admission;
        cx.checkpoint()
            .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
        let specimen = mechanical.specimen();
        let material = specimen.material();
        let thermal = curve
            .state_at_specific_enthalpy(specific_enthalpy_j_kg)
            .map_err(ThermalStringError::Phase)?;
        if !(maximum_temperature_increment_k.is_finite() && maximum_temperature_increment_k > 0.0)
            || card.content_hash() != material.card_identity()
            || thermal.material_card_identity() != material.card_identity()
            || !matches!(specimen.prestress(), StringPrestress::FixedTension(_))
        {
            return Err(refuse(
                "thermal string requires one card, fixed tension and a positive temperature-increment budget",
            ));
        }
        for name in ["density", "young_modulus"] {
            let property = material
                .property(name)
                .ok_or(refuse("thermal string needs density and Young's modulus"))?;
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(refuse("thermal string selected claim is absent"))?;
            if !matches!(claim.value, fs_matdb::PropertyValue::Scalar { .. }) {
                return Err(refuse(
                    "changing density or elastic stiffness requires thermodynamic work and state transfer",
                ));
            }
        }
        for property in material.properties() {
            cx.checkpoint()
                .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(refuse("thermal string selected claim is absent"))?;
            // Endpoint queries of a TabulatedOnly curve do not admit the path
            // between them. Constant claims and linear T-only curves do.
            if !property.is_constant_at_fixed_temperature()
                || matches!(claim.value, fs_matdb::PropertyValue::Curve { .. })
                    && claim.interpolation != fs_matdb::InterpolationPolicy::LinearInside
            {
                return Err(refuse(
                    "thermal string requires constant or linearly interpolated temperature-only claims",
                ));
            }
            if property.requirement().name()
                == crate::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY
                && let fs_matdb::PropertyValue::Curve { knots, .. } = &claim.value
                && knots.iter().any(|(_, viscosity)| *viscosity < 0.0)
            {
                return Err(refuse(
                    "thermal string viscosity must be nonnegative over its whole source curve",
                ));
            }
        }
        let density = material
            .property("density")
            .ok_or(refuse("thermal string needs density"))?
            .value_si();
        if curve.knots().iter().any(|knot| {
            knot.liquid_mass_fraction == 0.0
                && knot.bulk_density_kg_m3.to_bits() != density.to_bits()
        }) {
            return Err(refuse(
                "solid thermal density must equal the constant mechanical density",
            ));
        }
        let mut point = fs_matdb::QueryPoint::new();
        for (axis, value) in material.query_point() {
            point = match material.axis_quantities().get(axis) {
                Some(quantity) => point.with_quantity(axis, *quantity, *value),
                None => point.with(axis, *value),
            }
            .map_err(ThermalStringError::Query)?;
        }
        let absolute_temperature =
            fs_qty::QuantitySpec::semantic(fs_qty::semantic::SemanticType::new(
                fs_qty::semantic::QuantityKind::AbsoluteTemperature,
                fs_qty::semantic::ValueForm::Static,
            ));
        if point.axes().get("T").copied() != Some(thermal.temperature_k())
            || point.axis_quantities().get("T").is_some_and(|q| {
                *q != absolute_temperature
                    && *q != fs_qty::QuantitySpec::dimensional(fs_qty::Temperature::DIMS)
            })
            || thermal.phase() != SolidLiquidPhase::Solid
        {
            return Err(refuse(
                "initial solid thermal temperature must match the material's absolute T coordinate",
            ));
        }
        let requirements = material
            .properties()
            .iter()
            .map(|p| p.requirement().clone())
            .collect();
        let selection = MaterialPropertySelection::PinnedByProperty(
            material
                .properties()
                .iter()
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let result = Self {
            mechanical,
            card,
            point,
            requirements,
            selection,
            curve,
            thermal,
            maximum_temperature_increment_k,
        };
        // Validate the initial law against its original receipt. Pinning an
        // already selected claim changes the query-policy receipt identity;
        // that is not a physical mismatch. Future temperature queries are
        // pinned, but attaching this owner preserves the accepted source state.
        let rebound = result
            .mechanical
            .specimen()
            .clone()
            .with_kelvin_voigt_bending_loss()
            .map_err(ThermalStringError::Acoustic)?;
        if rebound.string() != result.mechanical.specimen().string() {
            return Err(refuse(
                "initial string must already use its sourced Kelvin-Voigt bending law",
            ));
        }
        cx.checkpoint()
            .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
        Ok(result)
    }

    /// Mechanics, material receipts, observer history and epoch, read only.
    #[must_use]
    pub const fn mechanical(&self) -> &LinearMaterialStringRuntime {
        &self.mechanical
    }

    /// Accepted uniform thermal state.
    #[must_use]
    pub const fn thermal(&self) -> EquilibriumPhaseState {
        self.thermal
    }

    fn conductor(&self) -> Result<fs_material::conductor::ResolvedConductor, ThermalStringError> {
        use fs_material::conductor::{ELECTRICAL_RESISTIVITY_PROPERTY, resolve_uniform_conductor};
        let specimen = self.mechanical.specimen();
        let property = specimen
            .material()
            .property(ELECTRICAL_RESISTIVITY_PROPERTY)
            .ok_or(ThermalStringError::Admission(
                "current-driven string requires a resolved electrical_resistivity claim",
            ))?;
        let point = self
            .point
            .clone()
            .with("T", self.thermal.temperature_k())
            .map_err(ThermalStringError::Query)?;
        resolve_uniform_conductor(
            &self.card,
            &point,
            MaterialPropertySelection::PinnedByProperty(vec![(
                ELECTRICAL_RESISTIVITY_PROPERTY.to_owned(),
                property.answer().receipt.selected,
            )]),
            specimen.string().length_m,
            specimen.area_m2(),
        )
        .map_err(ThermalStringError::Electrical)
    }

    /// Deposit `I² R(T_initial) dt` from an ideal prescribed current source.
    /// Resistance uses the same immutable card, pinned resistivity claim,
    /// temperature and actual length/area as mechanics. Both electrical states
    /// are resolved before publishing mechanics, temperature or sound.
    ///
    /// This is a first-order, uniform quasi-static Ohmic model, with insulated
    /// thermal boundaries. The external heat account is electrical source work
    /// converted once to Joule heat; material vibration loss stays separate.
    /// No circuit dynamics, skin/contact resistance, charge field, Lorentz
    /// force, motional emf or changing-geometry electrical work is represented.
    pub fn step_with_current(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        current_a: f64,
    ) -> Result<ElectrothermalStringFrame, ThermalStringError> {
        self.step_with_drive(
            cx,
            expected_epoch,
            forces,
            fs_material::conductor::OhmicDrive::Current(current_a),
        )
    }

    /// Prescribe voltage instead of current, with `I_initial = V/R(T_initial)`.
    /// The same first-order frozen-law and insulated-boundary model applies.
    /// Final current is resolved at the accepted temperature before publication.
    pub fn step_with_voltage(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        voltage_v: f64,
    ) -> Result<ElectrothermalStringFrame, ThermalStringError> {
        self.step_with_drive(
            cx,
            expected_epoch,
            forces,
            fs_material::conductor::OhmicDrive::Voltage(voltage_v),
        )
    }

    fn step_with_drive(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: fs_material::conductor::OhmicDrive,
    ) -> Result<ElectrothermalStringFrame, ThermalStringError> {
        let (conductor_before, current_a, heat) = self.prepare_drive(cx, expected_epoch, drive)?;
        let mut candidate = self.clone();
        let coupled = candidate.step(cx, expected_epoch, forces, heat)?;
        let conductor_after = candidate.conductor()?;
        let final_current_a = conductor_after
            .current_a_for_drive(drive)
            .map_err(ThermalStringError::Electrical)?;
        cx.checkpoint()
            .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
        *self = candidate;
        Ok(ElectrothermalStringFrame {
            coupled,
            drive,
            current_a,
            final_current_a,
            joule_heat_j: heat,
            boundary_heat_j: 0.0,
            transport: (),
            conductor_before,
            conductor_after,
        })
    }

    fn prepare_drive(
        &self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        drive: fs_material::conductor::OhmicDrive,
    ) -> Result<(fs_material::conductor::ResolvedConductor, f64, f64), ThermalStringError> {
        cx.checkpoint()
            .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        let conductor_before = self.conductor()?;
        let current_a = conductor_before
            .current_a_for_drive(drive)
            .map_err(ThermalStringError::Electrical)?;
        let power = conductor_before
            .joule_power_w(current_a)
            .map_err(ThermalStringError::Electrical)?;
        let heat = power * self.mechanical.member.model.sample_period_s();
        if !heat.is_finite() || (power > 0.0 && heat == 0.0) {
            return Err(ThermalStringError::Admission(
                "Joule heat is not representable in joules",
            ));
        }
        Ok((conductor_before, current_a, heat))
    }

    /// Couple prescribed-current heating to an existing thermal transport owner.
    /// The transport receives vibration plus Joule heat as its internal source;
    /// it must return only actual boundary heat in `external_heat_j`. This owner
    /// separately retains both inputs and includes electrical source work once
    /// in the coupled external-energy account. For a prescribed ambient use
    /// `|input| environment.advance(cx, input)`.
    ///
    /// The callback proposes state without publishing or debiting a source.
    /// Initial resistance is frozen for the sample, so this retains the same
    /// first-order electrical/material partition and uniform-body limits as
    /// `step_with_current`, while transport owns its thermal discretization.
    pub fn step_with_current_and_thermal_transport<R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        current_a: f64,
        transport: impl FnOnce(
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<fs_material::phase::UniformEnthalpyStep<R>, E>,
    ) -> Result<ElectrothermalStringFrame<R>, ThermalStringError> {
        self.step_with_drive_and_thermal_transport(
            cx,
            expected_epoch,
            forces,
            fs_material::conductor::OhmicDrive::Current(current_a),
            transport,
        )
    }

    /// Prescribed-voltage counterpart of `step_with_current_and_thermal_transport`.
    /// Joule heat uses the initial `V²/R`; temperature-dependent current, thermal
    /// transport, material laws and sound publish at the same accepted epoch.
    pub fn step_with_voltage_and_thermal_transport<R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        voltage_v: f64,
        transport: impl FnOnce(
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<fs_material::phase::UniformEnthalpyStep<R>, E>,
    ) -> Result<ElectrothermalStringFrame<R>, ThermalStringError> {
        self.step_with_drive_and_thermal_transport(
            cx,
            expected_epoch,
            forces,
            fs_material::conductor::OhmicDrive::Voltage(voltage_v),
            transport,
        )
    }

    /// Couple an ideal current or voltage source to caller-owned thermal state.
    /// The callback has the same proposal-only contract as
    /// [`Self::step_with_thermal_state`]. Both owners publish only after final
    /// electrical resolution as well as thermal/mechanical acceptance.
    /// Joule work is counted once by the existing electrical transport path.
    pub fn step_with_drive_and_thermal_state<S, R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: fs_material::conductor::OhmicDrive,
        state: &mut S,
        transport: impl FnOnce(
            &S,
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<(S, fs_material::phase::UniformEnthalpyStep<R>), E>,
    ) -> Result<ElectrothermalStringFrame<R>, ThermalStringError> {
        let frame = self.step_with_drive_and_thermal_transport(
            cx,
            expected_epoch,
            forces,
            drive,
            |input| {
                let (candidate, proposal) = transport(state, input)?;
                Ok::<_, E>(fs_material::phase::UniformEnthalpyStep {
                    state: proposal.state,
                    external_heat_j: proposal.external_heat_j,
                    energy_residual_tolerance_j: proposal.energy_residual_tolerance_j,
                    report: (candidate, proposal.report),
                })
            },
        )?;
        let _previous = core::mem::replace(state, frame.transport.0);
        Ok(ElectrothermalStringFrame {
            coupled: frame.coupled,
            drive: frame.drive,
            current_a: frame.current_a,
            final_current_a: frame.final_current_a,
            joule_heat_j: frame.joule_heat_j,
            boundary_heat_j: frame.boundary_heat_j,
            transport: frame.transport.1,
            conductor_before: frame.conductor_before,
            conductor_after: frame.conductor_after,
        })
    }

    fn step_with_drive_and_thermal_transport<R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: fs_material::conductor::OhmicDrive,
        transport: impl FnOnce(
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<fs_material::phase::UniformEnthalpyStep<R>, E>,
    ) -> Result<ElectrothermalStringFrame<R>, ThermalStringError> {
        let (conductor_before, current_a, heat) = self.prepare_drive(cx, expected_epoch, drive)?;
        let mut candidate = self.clone();
        let result =
            candidate.step_with_thermal_transport(cx, expected_epoch, forces, |mut input| {
                input.internal_heat_j += heat;
                if !input.internal_heat_j.is_finite() {
                    return Err("combined Joule and vibration heat is not finite".to_owned());
                }
                let proposed = transport(input).map_err(|error| error.to_string())?;
                Ok(fs_material::phase::UniformEnthalpyStep {
                    state: proposed.state,
                    external_heat_j: proposed.external_heat_j + heat,
                    energy_residual_tolerance_j: proposed.energy_residual_tolerance_j,
                    report: (proposed.external_heat_j, proposed.report),
                })
            })?;
        let conductor_after = candidate.conductor()?;
        let final_current_a = conductor_after
            .current_a_for_drive(drive)
            .map_err(ThermalStringError::Electrical)?;
        cx.checkpoint()
            .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
        *self = candidate;
        Ok(ElectrothermalStringFrame {
            coupled: result.coupled,
            drive,
            current_a,
            final_current_a,
            joule_heat_j: heat,
            boundary_heat_j: result.transport.0,
            transport: result.transport.1,
            conductor_before,
            conductor_after,
        })
    }

    fn specimen_at(
        &self,
        thermal: EquilibriumPhaseState,
    ) -> Result<crate::string_specimen::ResolvedStringSpecimen, ThermalStringError> {
        let point = self
            .point
            .clone()
            .with("T", thermal.temperature_k())
            .map_err(ThermalStringError::Query)?;
        let state = resolve_material_state_point(
            &self.card,
            &point,
            &self.requirements,
            self.selection.clone(),
        )
        .map_err(ThermalStringError::Material)?;
        let old = self.mechanical.specimen();
        crate::string_specimen::with_uniform_circular_material_and_constraints(
            old.string(),
            old.geometry_constraint(),
            &state,
            old.prestress(),
        )
        .and_then(|s| s.with_kelvin_voigt_bending_loss())
        .map_err(ThermalStringError::Acoustic)
    }

    /// Advance mechanics, material heating and viscosity feedback as one commit.
    /// `expected_epoch` prevents replaying an accepted heat/source application.
    /// Signed external heat is supplied as energy [J], not an inferred ambient
    /// flux. A refusal consumes neither this input nor an audio sample.
    pub fn step(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        external_heat_j: f64,
    ) -> Result<ThermalMaterialStringFrame, ThermalStringError> {
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        if !external_heat_j.is_finite() {
            return Err(ThermalStringError::Admission(
                "external heat must be finite joules",
            ));
        }
        let mut candidate = self.mechanical.clone();
        let vibration = candidate
            .step(cx, forces)
            .map_err(ThermalStringError::Acoustic)?;
        let heat = vibration.dissipation.material_heat_j + external_heat_j;
        let thermal = self
            .curve
            .advance_specific_energy(self.thermal, heat / candidate.specimen().mass_kg())
            .map_err(ThermalStringError::Phase)?;
        let result = self.finish_thermal_step(
            cx,
            forces,
            candidate,
            vibration,
            fs_material::phase::UniformEnthalpyStep {
                state: thermal,
                external_heat_j,
                energy_residual_tolerance_j: 0.0,
                report: (),
            },
        )?;
        Ok(result.coupled)
    }

    /// Advance material heating and boundary heat exchange through a supplied
    /// transport owner. For fs-conduction's `LumpedThermalEnvironment`, pass
    /// `|input| environment.advance(cx, input)`; no reverse crate dependency
    /// from this coupled owner to a spatial solver is required.
    ///
    /// Area is `2 pi r (L + r)` and the Biot length is the actual volume/area.
    /// The transport must justify its lumped approximation and return actual
    /// boundary heat, separately from internal heating and solver residual.
    /// It proposes a result; it must not publish state or debit a mutable source.
    pub fn step_with_thermal_transport<R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        transport: impl FnOnce(
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<fs_material::phase::UniformEnthalpyStep<R>, E>,
    ) -> Result<ThermalStringTransportFrame<R>, ThermalStringError> {
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        let specimen = self.mechanical.specimen();
        let radius = specimen.radius_m();
        let length = specimen.string().length_m;
        let area = 2.0 * core::f64::consts::PI * radius * (length + radius);
        let volume = specimen.area_m2() * length;
        let mut candidate = self.mechanical.clone();
        let vibration = candidate
            .step(cx, forces)
            .map_err(ThermalStringError::Acoustic)?;
        let dt = candidate.member.model.sample_period_s();
        let proposed = transport(fs_material::phase::UniformEnthalpyStepInput {
            curve: &self.curve,
            initial: self.thermal,
            mass_kg: specimen.mass_kg(),
            volume_m3: volume,
            surface_area_m2: area,
            internal_heat_j: vibration.dissipation.material_heat_j,
            duration_s: dt,
        })
        .map_err(|error| ThermalStringError::Transport(error.to_string()))?;
        if proposed.state.phase_curve_identity() != self.curve.identity()
            || !proposed.external_heat_j.is_finite()
            || !proposed.energy_residual_tolerance_j.is_finite()
            || proposed.energy_residual_tolerance_j < 0.0
        {
            return Err(ThermalStringError::Admission(
                "thermal transport returned a foreign chart, invalid heat or invalid residual budget",
            ));
        }
        self.finish_thermal_step(cx, forces, candidate, vibration, proposed)
    }

    /// Advance with caller-owned thermal state, such as a finite support's
    /// enthalpy, publishing that state only after coupled acceptance.
    ///
    /// The callback reads accepted state and returns an owned candidate with
    /// the usual thermal proposal. It must not mutate shared interior state,
    /// debit external resources or emit observations. All transport, material,
    /// energy, epoch and cancellation refusals preserve both accepted owners.
    /// This method does not validate the opaque state's physical meaning;
    /// that remains the transport owner's responsibility. Publication uses
    /// exclusive borrows, not a cross-thread or durable transaction.
    pub fn step_with_thermal_state<S, R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        state: &mut S,
        transport: impl FnOnce(
            &S,
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<(S, fs_material::phase::UniformEnthalpyStep<R>), E>,
    ) -> Result<ThermalStringTransportFrame<R>, ThermalStringError> {
        let frame = self.step_with_thermal_transport(cx, expected_epoch, forces, |input| {
            let (candidate, proposal) = transport(state, input)?;
            Ok::<_, E>(fs_material::phase::UniformEnthalpyStep {
                state: proposal.state,
                external_heat_j: proposal.external_heat_j,
                energy_residual_tolerance_j: proposal.energy_residual_tolerance_j,
                report: (candidate, proposal.report),
            })
        })?;
        // No fallible operation or cancellation point separates publication.
        // Move the old value out so even its destructor runs after both commits.
        let _previous = core::mem::replace(state, frame.transport.0);
        Ok(ThermalStringTransportFrame {
            coupled: frame.coupled,
            transport: frame.transport.1,
        })
    }

    fn finish_thermal_step<R>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        forces: &[f64],
        mut candidate: LinearMaterialStringRuntime,
        mut vibration: MaterialStringFrame,
        proposed: fs_material::phase::UniformEnthalpyStep<R>,
    ) -> Result<ThermalStringTransportFrame<R>, ThermalStringError> {
        let fs_material::phase::UniformEnthalpyStep {
            state: thermal,
            external_heat_j,
            energy_residual_tolerance_j: thermal_solve_tolerance_j,
            report,
        } = proposed;
        let energy_before_j: f64 = self
            .mechanical
            .modes()
            .iter()
            .zip(self.mechanical.states())
            .map(|(mode, state)| {
                0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                    + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2))
            })
            .sum();
        let mass = candidate.specimen().mass_kg();
        let heat = vibration.dissipation.material_heat_j + external_heat_j;
        if thermal.phase() != SolidLiquidPhase::Solid
            || (thermal.temperature_k() - self.thermal.temperature_k()).abs()
                > self.maximum_temperature_increment_k
        {
            return Err(ThermalStringError::Admission(
                "thermal step left the solid phase or exceeded its temperature-increment budget",
            ));
        }
        let specimen = self.specimen_at(thermal)?;
        // The internal rebind shares this sample's publication epoch. Start
        // from the accepted epoch so the last representable epoch is usable.
        candidate.epoch = self.mechanical.epoch();
        candidate
            .rebind(cx, specimen, 0.0)
            .map_err(ThermalStringError::Acoustic)?;
        // Rebinding is internal to this sample's transaction, not a second
        // accepted instant. Observe the final acceleration with the previous
        // accepted history, so pressure and its retained history use the new law.
        candidate.member.observer = self.mechanical.member.observer.clone();
        let initial = candidate
            .member
            .observer
            .previous_acceleration
            .is_none()
            .then(|| linear_accelerations(&self.mechanical.member.model, forces));
        vibration.acoustic.observer_pressure_pa = candidate
            .member
            .observer
            .observe(
                linear_accelerations(&candidate.member.model, forces),
                initial,
                candidate.member.model.sample_period_s(),
            )
            .map_err(ThermalStringError::Acoustic)?;
        let enthalpy_change_j =
            mass * (thermal.specific_enthalpy_j_kg() - self.thermal.specific_enthalpy_j_kg());
        let a = &vibration.acoustic;
        let d = &vibration.dissipation;
        let energy_balance_residual_j = a.input_work_j + external_heat_j
            - ((a.total_modal_energy_j - energy_before_j)
                + enthalpy_change_j
                + d.air_loss_j
                + d.authored_loss_j
                + d.roundoff_residual_j);
        let scale = a.input_work_j.abs()
            + external_heat_j.abs()
            + energy_before_j
            + a.total_modal_energy_j
            + mass
                * (thermal.specific_enthalpy_j_kg().abs()
                    + self.thermal.specific_enthalpy_j_kg().abs())
            + heat.abs()
            + d.air_loss_j
            + d.authored_loss_j;
        let energy_roundoff_tolerance_j = d.roundoff_tolerance_j + 32.0 * f64::EPSILON * scale;
        if !energy_balance_residual_j.is_finite()
            || !energy_roundoff_tolerance_j.is_finite()
            || energy_balance_residual_j.abs()
                > energy_roundoff_tolerance_j + thermal_solve_tolerance_j
        {
            return Err(ThermalStringError::Admission(
                "coupled thermal/mechanical energy balance exceeds its roundoff and solve allowances",
            ));
        }
        cx.checkpoint()
            .map_err(|_| ThermalStringError::Acoustic(AcousticRealizeError::Cancelled))?;
        self.mechanical = candidate;
        self.thermal = thermal;
        Ok(ThermalStringTransportFrame {
            coupled: ThermalMaterialStringFrame {
                vibration,
                thermal,
                external_heat_j,
                enthalpy_change_j,
                energy_balance_residual_j,
                energy_roundoff_tolerance_j,
                thermal_solve_tolerance_j,
            },
            transport: report,
        })
    }
}

fn next_string_epoch(epoch: u64) -> Result<u64, AcousticRealizeError> {
    epoch
        .checked_add(1)
        .ok_or(AcousticRealizeError::InvalidDescription {
            what: "incremental string epoch or sample counter overflow",
        })
}

fn validate_incremental_string(string: &PrestressedString) -> Result<(), AcousticRealizeError> {
    validate_string(string, None, None)?;
    if string.moving_end
        || string.polarization_detune != 0.0
        || string.relaxation_bending.is_some()
        || string.n_modes > MAX_TIME_DOMAIN_ACOUSTIC_MODES
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "incremental linear string needs one fixed-end sine basis without relaxation memory, within the modal budget",
        });
    }
    Ok(())
}

fn linear_string_member(
    string: &PrestressedString,
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
    let mut loss_parts = Vec::with_capacity(string.n_modes);
    for k in 1..=string.n_modes {
        let omega = string_mode_omega(string, k);
        let q_phys = pluck.map_or(0.0, |p| triangular_pluck_modal(p, k));
        let loss = mode_loss_parts(string, omega, k as f64 * pi / string.length_m, gas)?;
        modes.push(ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: loss.total(),
            // The time-domain observer below includes actual applied forces.
            pressure_per_modal_velocity: C64::ZERO,
        });
        loss_parts.push(loss);
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
    Ok(LinearMember {
        model,
        phi_bow,
        observer: StringObserver::new(string, gas, listener_m, radiation_scale),
        loss_parts,
    })
}

fn linear_accelerations(model: &ModalAcousticTimeModel, force: &[f64]) -> Vec<f64> {
    model
        .modes()
        .iter()
        .zip(model.states())
        .zip(force)
        .map(|((mode, state), &f)| {
            let omega = mode.angular_frequency_rad_s;
            f - omega * omega * state.displacement_m_sqrt_kg
                - 2.0 * mode.damping_ratio * omega * state.velocity_m_sqrt_kg_per_s
        })
        .collect()
}

fn step_linear_member(
    member: &mut LinearMember,
    bow: Option<BowStroke>,
    texture: &mut TextureDrive,
    obstacles: &[UnilateralObstacle],
    string: &PrestressedString,
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
        )
        .map_err(AcousticRealizeError::Nonlinear)?;
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
    advance_linear_member(member, &force, dt)
}

fn advance_linear_member(
    member: &mut LinearMember,
    force: &[f64],
    dt: f64,
) -> Result<f64, AcousticRealizeError> {
    let initial = member
        .observer
        .previous_acceleration
        .is_none()
        .then(|| linear_accelerations(&member.model, force));
    member
        .model
        .step(force)
        .map_err(AcousticRealizeError::Modal)?;
    member
        .observer
        .observe(linear_accelerations(&member.model, force), initial, dt)
}

#[allow(clippy::too_many_arguments)]
fn realize_kc_string(
    string: &PrestressedString,
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
    let twin = secondary_string(string, twin_tension.unwrap_or(1.0));
    let mut texture = TextureDrive::try_new(texture)?;
    let mut members = vec![kc_string_member(
        string, pluck, bow, obstacles, gas, listener_m, 1.0,
    )?];
    if twin_tension.is_some() {
        members.push(kc_string_member(
            &twin, pluck, bow, obstacles, gas, listener_m, 0.85,
        )?);
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let mut p = 0.0;
        let mut fb = 0.0;
        for (idx, member) in members.iter_mut().enumerate() {
            let member_string = if idx == 0 { string } else { &twin };
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
    observer: StringObserver,
}

fn kc_string_member(
    string: &PrestressedString,
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
        .enumerate()
        .map(|(k, &w)| {
            let wave_number = (k + 1) as f64 * core::f64::consts::PI / string.length_m;
            mode_zeta(string, w, wave_number, gas)
        })
        .collect();
    let zetas = zetas?;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let sys = with_bending_relaxation(string, assemble_kc(string, storage, &zetas, obstacles)?)?;
    let mut x = vec![0.0; sys.state_dim()];
    if let Some(pluck) = pluck {
        for k in 1..=string.n_modes {
            x[2 * (k - 1)] = mass_scale * triangular_pluck_modal(pluck, k);
        }
    }
    initialize_relaxed_bending(string, &mut x);
    Ok(KcMember {
        sys,
        x,
        observer: StringObserver::new(string, gas, listener_m, radiation_scale),
    })
}

/// Momentum rows of `(J-R) grad H + G u` are accelerations in the
/// mass-normalized string basis. Read the admitted structure and storage so
/// nonlinear elasticity and conservative contact reach the same observer.
fn kc_accelerations(sys: &fs_phs::PortHamiltonian, x: &[f64], u: &[f64]) -> Vec<f64> {
    let effort = sys.effort(x);
    let (j, r, g) = sys.structure();
    let n = sys.state_dim();
    let m = sys.port_dim();
    (0..m)
        .map(|k| {
            let row = 2 * k + 1;
            let internal: f64 = effort
                .iter()
                .enumerate()
                .map(|(col, e)| (j[row * n + col] - r[row * n + col]) * e)
                .sum();
            internal
                + u.iter()
                    .enumerate()
                    .map(|(col, f)| g[row * m + col] * f)
                    .sum::<f64>()
        })
        .collect()
}

fn step_kc_member(
    member: &mut KcMember,
    bow: Option<BowStroke>,
    texture: &mut TextureDrive,
    obstacles: &[UnilateralObstacle],
    string: &PrestressedString,
    _gas: &GasState,
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
        )
        .map_err(AcousticRealizeError::Nonlinear)?;
        #[allow(clippy::needless_range_loop)] // the modal index spans u and shapes
        for k in 0..string.n_modes {
            let phi = det::sin((k + 1) as f64 * pi * bow.station_frac) / mass_scale;
            u[k] = f_bow * phi;
        }
    }
    if obstacles.iter().any(|o| o.mu_kinetic > 0.0) {
        let extra = modal_friction_forces(string, obstacles, &member.x[..2 * string.n_modes])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        for (uk, fk) in u.iter_mut().zip(extra) {
            *uk += fk;
        }
    }
    if obstacles.iter().any(|o| o.internal_loss > 0.0) {
        let extra = modal_hunt_crossley_forces(string, obstacles, &member.x[..2 * string.n_modes])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        for (uk, fk) in u.iter_mut().zip(extra) {
            *uk += fk;
        }
    }
    let initial = member
        .observer
        .previous_acceleration
        .is_none()
        .then(|| kc_accelerations(&member.sys, &member.x, &u));
    let rec = step(&member.sys, &member.x, &u, dt)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let p = member
        .observer
        .observe(kc_accelerations(&member.sys, &rec.x, &u), initial, dt)?;
    member.x = rec.x;
    let mut q_phys = vec![0.0; string.n_modes];
    #[allow(clippy::needless_range_loop)] // the modal index spans state and outputs
    for k in 0..string.n_modes {
        let q = member.x[2 * k];
        q_phys[k] = q / mass_scale;
    }
    Ok((p, q_phys))
}

fn triangular_pluck_modal(pluck: Pluck, k: usize) -> f64 {
    let xi = pluck.station_frac;
    let kf = k as f64;
    let pi = core::f64::consts::PI;
    2.0 * pluck.height_m * det::sin(kf * pi * xi) / (kf * kf * pi * pi * xi * (1.0 - xi))
}

fn bow_force(bow: BowStroke, v_string: f64, normal_delta_n: f64) -> Result<f64, String> {
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
        .map_err(|e| e.to_string())
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

fn bridge_force(string: &PrestressedString, q_phys: &[f64]) -> f64 {
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

/// Fail-closed admission for tone holes: every realization lane must agree
/// on whether a hole exists, so malformed entries refuse the WHOLE duct
/// instead of being dropped by one image and kept by the other. This also
/// closes the mis-indexed-hole escape hatch: `after_segment` beyond the
/// last segment previously matched no realization loop at all.
fn validate_tone_holes(duct: &ViscothermalDuct) -> Result<(), &'static str> {
    if duct.tone_holes.len() > duct.segments.len().saturating_mul(64) {
        return Err("unreasonable tone-hole count for the duct");
    }
    for hole in &duct.tone_holes {
        if hole.after_segment >= duct.segments.len() {
            return Err("tone hole references a segment index past the duct end");
        }
        if !(hole.radius_m.is_finite() && hole.radius_m > 0.0) {
            return Err("tone-hole radius must be finite and positive");
        }
        if !(hole.chimney_m.is_finite() && hole.chimney_m >= 0.0) {
            return Err("tone-hole chimney length must be finite and nonnegative");
        }
        if !hole.open_fraction.is_finite() {
            return Err("tone-hole open fraction must be finite");
        }
    }
    Ok(())
}

fn bore_spec(duct: &ViscothermalDuct) -> Option<(Vec<AcousticSection>, Vec<AcousticTap>)> {
    if validate_tone_holes(duct).is_err() {
        return None;
    }
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
            p += body.radiate_mass_normalized(acc, gas.density, listener_m);
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
    // Parity with ReedBoreVoice::new (render.rs): the ODE path must
    // reject non-finite/nonphysical attack ramps and reed stiffnesses
    // up front, or NaN slips past typed refusals into generic per-step
    // Nonlinear noise mid-integration (bead frankensim-dtb76).
    if !(reed.rest_opening_m > 0.0
        && reed.width_m > 0.0
        && reed.closing_pressure_pa > 0.0
        && reed.blowing_pressure_pa >= 0.0
        && reed.attack_s >= 0.0
        && reed.mass_kg >= 0.0
        && reed.stiffness_n_m >= 0.0
        && reed.damping_ratio >= 0.0
        && reed.damping_ratio.is_finite()
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
        let face = reed_pressure_face(reed);
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
                p += body.radiate_mass_normalized(acc, gas.density, listener_m);
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
    validate_tone_holes(duct).map_err(|what| AcousticRealizeError::InvalidDescription { what })?;
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

#[cfg(test)]
mod string_observer_tests {
    use super::*;

    #[test]
    fn g1_plate_ports_and_observer_use_the_same_modal_mass() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let mut base = crate::thin_plate::CompactBody::from_radiator(RadiatingPlate {
            area_m2: 0.3,
            mass_kg: 4.0,
            frequency_hz: 100.0,
            damping_ratio: 0.0,
        })
        .unwrap();
        base.area_m2 = -0.12;
        base.drive_participation = 0.25;
        let (force, face_pressure, dt) = (2.0, 3.0, 1.0e-5);
        let expected = gas.density * -0.12 / (2.0 * core::f64::consts::PI)
            * (0.25 * force - 0.12 * face_pressure)
            / 4.0
            / (1.0 + (base.omega * dt / 2.0).powi(2));
        for scale in [1.0, -4.0, 0.125] {
            let mut body = base.clone();
            body.mass_kg *= scale * scale;
            body.area_m2 *= scale;
            body.drive_participation *= scale;
            let mut plates = PlateBank::default();
            plates.linear.push(body);
            let sys = plate_modal_ports(&plates, true).unwrap().unwrap();
            let zero = vec![0.0; sys.state_dim()];
            let rec = step(&sys, &zero, &[force, face_pressure], dt).unwrap();
            let p = plate_chain_radiation(&plates, &rec.x, &zero, 2, &gas, 1.0, dt);
            assert!(
                (p / expected - 1.0).abs() < 1.0e-10,
                "mass-normalized pressure {p} vs {expected}"
            );
            let physical_v = rec.x[1] / plates.linear[0].mass_kg.sqrt();
            let output = sys.output(&rec.x);
            assert!(
                (output[0] - plates.linear[0].drive_participation * physical_v).abs() < 1.0e-14
            );
            assert!((output[1] - plates.linear[0].area_m2 * physical_v).abs() < 1.0e-14);
        }
    }

    #[test]
    fn g1_string_observer_held_static_load_is_silent() {
        let string = PrestressedString {
            length_m: 0.5,
            tension_n: 20.0,
            lin_density_kg_m: 0.006,
            axial_stiffness_n: 0.0,
            width_m: 0.001,
            n_modes: 2,
            damping_ratio: 0.0,
            rayleigh: Some(RayleighParams {
                alpha_per_s: 2.0,
                beta_s: 0.0,
            }),
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
        };
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let mut member = linear_string_member(&string, None, None, &gas, 1.0, 8_000, 1.0).unwrap();
        let force = [1.0, 2.0];
        member.model.initialize_static_equilibrium(&force).unwrap();
        assert!(
            member
                .model
                .states()
                .iter()
                .all(|s| s.displacement_m_sqrt_kg > 0.0)
        );
        for _ in 0..8 {
            let pressure = advance_linear_member(&mut member, &force, 1.0 / 8_000.0).unwrap();
            assert!(
                pressure.abs() < 1.0e-12,
                "stationary loaded string radiated {pressure} Pa"
            );
        }
    }
}

#[cfg(test)]
mod tone_hole_admission_tests {
    use super::*;
    use fs_scenario::{CylinderSegment, ToneHole};

    fn base_duct(hole: Option<ToneHole>) -> ViscothermalDuct {
        ViscothermalDuct {
            segments: vec![CylinderSegment::cylinder(0.01, 0.3)],
            tone_holes: hole.into_iter().collect(),
            termination: WaveguideEnd::Closed,
            wall: None,
        }
    }

    fn malformed_variants() -> Vec<ToneHole> {
        vec![
            ToneHole {
                after_segment: 0,
                radius_m: f64::NAN,
                chimney_m: 0.0,
                open_fraction: 1.0,
            },
            ToneHole {
                after_segment: 5,
                radius_m: 0.004,
                chimney_m: 0.0,
                open_fraction: 1.0,
            },
            ToneHole {
                after_segment: 0,
                radius_m: -0.004,
                chimney_m: 0.0,
                open_fraction: 1.0,
            },
            ToneHole {
                after_segment: 0,
                radius_m: 0.004,
                chimney_m: 0.0,
                open_fraction: f64::NAN,
            },
        ]
    }

    #[test]
    fn bore_spec_refuses_every_malformed_tone_hole_shape() {
        for hole in malformed_variants() {
            assert!(
                bore_spec(&base_duct(Some(hole))).is_none(),
                "bore_spec must refuse duct with malformed hole {hole:?}"
            );
        }
        assert!(bore_spec(&base_duct(None)).is_some(), "control must pass");
    }

    #[test]
    fn physics_duct_refuses_with_typed_error() {
        for hole in malformed_variants() {
            let error = physics_duct(&base_duct(Some(hole)))
                .expect_err("physics_duct must refuse malformed holes");
            let what = match &error {
                AcousticRealizeError::InvalidDescription { what } => *what,
                other => panic!("unexpected refusal shape: {other:?}"),
            };
            assert!(what.contains("tone"), "{what}");
        }
        assert!(physics_duct(&base_duct(None)).is_ok(), "control must pass");
    }
}
