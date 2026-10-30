//! Realize an [`fs_scenario::AcousticAssembly`] into observer pressure.
//!
//! A guitar or clarinet is a description: Kirchhoff–Carrier strings,
//! Rayleigh damping, a radiating plate, a beating reed on a viscothermal
//! TMM bore, optional tone holes, and a regularized bow. There is no
//! instrument crate.

/// Finite source-backed gas storage coupled to the live material string.
pub mod finite_gas;
/// Source-backed fixed-support expansion and elastic storage during vibration.
pub mod thermoelastic;
/// Layer-resolved heat and continuing vibration on the same material sources.
pub mod concentric;
mod thermal_step;

use crate::modal_acoustic_time::{
    MAX_TIME_DOMAIN_ACOUSTIC_MODES, ModalAcousticFrame, ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeError, ModalAcousticTimeModel,
};
use crate::pcm_wav::{WavError, encode_pcm16_wav};
use crate::reed_bore::{
    blowing_envelope, reed_pressure_face, reed_structural, reed_swept_flow,
};
use crate::thin_plate::{
    PlateBank, PlateBankObservation, PlateChartRadiation, ResolvedPlateChart, VkBody, certified_chart_radiators,
    certified_material_chart_radiators, certified_radiators, vk_plate_phs,
};
use crate::unilateral_contact::{
    modal_contact_forces, modal_friction_forces, modal_hunt_crossley_forces, slit_contact_force,
    slit_lay,
};
/// Legacy sections-per-metre budget for default non-reed ODE duct paths.
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
    AcousticSection, AcousticTap, ChainOutlet, MouthFlange, ViscothermalPin, WallPin,
    acoustic_chain_outlet_wall, bernoulli_volume_flow, join_port, mass_spring_damper, modal_bank,
    modal_bank_ports, quasistatic_aperture_opening, step, step_descriptor, transformer,
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
    /// A non-exact block exceeded the clock's explicit <=1-cent linear band.
    TimeResolution(fs_phs::linear_exact::ResolvedBandError),
    /// The execution scope refused work before candidate publication.
    Cancelled,
    /// Description failed a physical domain check.
    InvalidDescription {
        /// Which field.
        what: &'static str,
    },
    /// Gas-state derivation failed.
    Ambient(String),
    /// ISO 9613-1 was asked of a gas not admitted as terrestrial air; the
    /// refusal names the species.
    NonTerrestrialAir(fs_material::iso9613::Iso9613Refusal),
    /// A receiver too close for its observer filter. Its own propagation delay
    /// at the nearest outlet element cannot absorb the filter's low-frequency
    /// delay (dispersion plus declared latency) and keep the receiver's
    /// one-sample minimum.
    ObserverLatency {
        /// Receiver distance from the outlet centre [m].
        range_m: f64,
        /// Propagation delay at the nearest outlet element at `c_0` [samples].
        available_samples: f64,
        /// Delay the filter and the receiver need there [samples].
        needed_samples: f64,
        /// The certificate the filter was designed to: `min(1e-9, 1e-3 x`
        /// what ISO 9613-1's +-10 % on alpha is in response energy).
        budget_energy: f64,
    },
    /// A render that applies no observer law, or Stokes-Kirchhoff only, would
    /// leave out more absorption than the tolerance, or an unbounded amount
    /// (a polyatomic listener gas without sourced relaxation rates), without
    /// declaring a law. The receipt carries the figure.
    ObserverOmission(Box<crate::air_path::ObserverOmission>),
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
    /// Algebraic work/flow refusal with the original quantitative certificate.
    AlgebraicSolve(fs_phs::PhsError),
    /// Lumped LC ladder requested outside its admissible pitch-error band.
    LumpedLadderBand {
        /// Requested total cell count.
        cells: usize,
        /// Smallest cell count whose spatial dispersion alone meets the
        /// tolerance. This is a necessary condition, not an admission.
        min_admissible_cells: usize,
        /// Refusing pitch error at `cells` [cents]: the dispersion bound below
        /// `min_admissible_cells`, otherwise the measured passive error
        /// (infinite when the ladder has no resonance within ±386 cents of
        /// the line's).
        pitch_error_cents: f64,
        /// Declared pitch-error tolerance [cents].
        max_pitch_error_cents: f64,
    },
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
            Self::TimeResolution(error) => core::fmt::Display::fmt(error, f),
            Self::Cancelled => write!(f, "FS-COUPLE-ASSEMBLY-CANCELLED"),
            Self::InvalidDescription { what } => {
                write!(f, "FS-COUPLE-ASSEMBLY: {what}")
            }
            Self::Ambient(what) => write!(f, "FS-COUPLE-ASSEMBLY-GAS: {what}"),
            Self::NonTerrestrialAir(refusal) => write!(f, "FS-COUPLE-ASSEMBLY-ISO9613: {refusal}"),
            Self::ObserverLatency { range_m, available_samples, needed_samples, budget_energy } => write!(
                f,
                "FS-COUPLE-ASSEMBLY-OBSERVER-LATENCY: a receiver at {range_m} m has \
                 {available_samples:.3} samples of propagation delay; its observer filter, \
                 certified to {budget_energy:e} of its response energy, and the receiver \
                 need {needed_samples:.3}. This law's dispersion cannot be certified this \
                 close: declare observer_law none or move the receiver farther"
            ),
            Self::ObserverOmission(omission) => write!(
                f,
                "FS-COUPLE-ASSEMBLY-OBSERVER-OMISSION: {omission}. Declare the observer law \
                 (observer_law none or stokes-kirchhoff) to accept the omission, or move the \
                 listener closer"
            ),
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
            Self::AlgebraicSolve(e) => write!(f, "FS-COUPLE-ASSEMBLY-BALANCE: {e}"),
            Self::LumpedLadderBand {
                cells,
                min_admissible_cells,
                pitch_error_cents,
                max_pitch_error_cents,
            } => write!(
                f,
                "FS-COUPLE-LUMPED-LADDER-BAND: requested {cells} cells has model-form pitch error \
                 {pitch_error_cents:.4} cents > {max_pitch_error_cents:.4} cents (minimum admissible \
                 cells: {min_admissible_cells})"
            ),
            Self::Reed { what } => write!(f, "FS-COUPLE-REED: {what}"),
            Self::Wav(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for AcousticRealizeError {}

/// Bore time-domain discretization selector for [`RealizeOptions`].
///
/// The default is the convergence-sized TMM characteristic line
/// ([`Self::CharacteristicLine`]). The lumped LC ladder ([`Self::LumpedLadder`])
/// is retained only as an explicitly requested model-form rung with a refinable
/// cell-count axis and an enforced pitch-error band.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoreRealization {
    /// Convergence-sized TMM characteristic line (`M2.2-W2`, default).
    CharacteristicLine(crate::driving_point::LineLadder),
    /// Explicit lumped LC ladder model-form rung (`realize_reed_ode`).
    LumpedLadder {
        /// Total lumped LC cells across the bore (refinable axis, no 8-cell budget clamp).
        cells: usize,
        /// Maximum admitted model-form pitch error against the line [cents].
        max_pitch_error_cents: f64,
    },
}

impl Default for BoreRealization {
    fn default() -> Self {
        Self::CharacteristicLine(crate::driving_point::LineLadder::default())
    }
}

impl BoreRealization {
    /// Request an explicit `cells`-section lumped LC ladder with the default
    /// 1.0-cent model-form pitch-error ceiling.
    #[must_use]
    pub const fn lumped_ladder(cells: usize) -> Self {
        Self::LumpedLadder {
            cells,
            max_pitch_error_cents: 1.0,
        }
    }

    /// Request an explicit `cells`-section lumped LC ladder with a caller-supplied
    /// pitch-error ceiling [cents].
    #[must_use]
    pub const fn lumped_ladder_with_tolerance(cells: usize, max_pitch_error_cents: f64) -> Self {
        Self::LumpedLadder {
            cells,
            max_pitch_error_cents,
        }
    }
}

/// Provenance and convergence/model-form receipt for a realized bore.
#[derive(Debug, Clone, PartialEq)]
pub enum BoreRealizationReceipt {
    /// Realized via the convergence-sized TMM characteristic line.
    CharacteristicLine(crate::driving_point::LineRealization),
    /// Realized via the explicit lumped LC ladder model-form rung.
    LumpedLadder {
        /// Total lumped LC cells assembled across the bore.
        cells: usize,
        /// Smallest cell count whose spatial dispersion alone meets the
        /// tolerance. This is a necessary condition, not an admission.
        min_admissible_cells: usize,
        /// Spatial dispersion model-form pitch error at `cells` [cents].
        dispersion_error_cents: f64,
        /// Measured, signed model-form pitch error against the characteristic
        /// line [cents]. It is taken from the lowest passive input-impedance
        /// resonance of each: the ladder's port-Hamiltonian `Z(iω)` against the
        /// line's own TMM. Positive means the ladder is sharp.
        pitch_error_cents: f64,
        /// Declared pitch-error tolerance [cents].
        max_pitch_error_cents: f64,
    },
}

/// Non-IR realization options (`scenario.ir` remains frozen per `schema-policy.json`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct RealizeOptions<'a> {
    /// Bore discretization choice for reed/duct realizations. Nondefault
    /// requests refuse on joined string/duct and prescribed-flow paths whose
    /// propagators do not implement this selector.
    pub bore: BoreRealization,
    /// Explicit, coherently resolved surrounding gas. Every excitation,
    /// propagation, loading and pressure-observation consumer receives this
    /// same state. `None` resolves the scenario's ambient-air declaration.
    /// Solid material temperature and prestress remain independent inputs.
    pub gas: Option<GasState>,
    /// Explicit observer absorption law, independent of gas-state provenance
    /// and the mechanical discretization. `None` preserves the scenario-air
    /// path when `gas` is absent, and selects Stokes-Kirchhoff for supplied gas.
    /// A supplied NASA-9 air state therefore needs an explicit ISO declaration
    /// to include that air-specific molecular-relaxation model.
    pub observer_absorption: Option<crate::air_path::GasAbsorptionModel>,
    /// The gas's molecular relaxation, evaluated at exactly the realized gas
    /// state. It is required by, and only by,
    /// [`crate::air_path::GasAbsorptionModel::MolecularRelaxation`]. A state
    /// evaluated for another gas, or a state supplied with another law, refuses
    /// instead of being ignored.
    pub observer_relaxation: Option<&'a fs_material::relaxation::RelaxationState>,
}

/// The observer law and, for molecular relaxation, its evaluated state.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct ObserverLaw<'a> {
    model: Option<crate::air_path::GasAbsorptionModel>,
    relaxation: Option<&'a fs_material::relaxation::RelaxationState>,
}

impl From<Option<crate::air_path::GasAbsorptionModel>> for ObserverLaw<'_> {
    fn from(model: Option<crate::air_path::GasAbsorptionModel>) -> Self {
        Self { model, relaxation: None }
    }
}

/// Realized observer pressure and the gas state that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct RealizedAssembly {
    /// Sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Observer pressure [Pa], one sample per period.
    pub pressure_pa: Vec<f64>,
    /// Derived ambient state (rho, c, mu, …).
    pub gas: GasState,
    /// Bore realization receipt when a reed bore was realized.
    pub bore_receipt: Option<BoreRealizationReceipt>,
    /// The observer absorption law applied, and any molecular relaxation it
    /// left out, with an upper bound on its effect at the listener.
    pub observer_absorption: crate::air_path::ObserverAbsorptionReceipt,
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
    realize_assembly_with_options(assembly, &RealizeOptions::default())
}

/// Realize an assembly with explicit non-IR [`RealizeOptions`] (such as an
/// explicit [`BoreRealization`] model-form rung).
///
/// # Errors
/// Ordinary assembly refusals plus [`AcousticRealizeError::LumpedLadderBand`]
/// when an explicit lumped LC ladder is requested outside its admissible band.
pub fn realize_assembly_with_options(
    assembly: &AcousticAssembly,
    options: &RealizeOptions,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    realize_assembly_inner(assembly, None, options)
}

/// Realize the same assembly in an explicitly supplied gas state.
///
/// This is the general material/composition ingress. Density, sound speed,
/// viscosity, thermal conductivity, heat capacities and impedance all come
/// from the supplied state, so a caller can swap an admitted gas or mixture
/// without rewriting instrument or solver code. The scenario AmbientGas is
/// ignored for constitutive coefficients on this path.
///
/// ISO 9613-1 is an atmospheric-air model and is not applied to arbitrary
/// mixtures. Free-field propagation therefore uses the existing generic
/// Stokes-Kirchhoff model rather than pretending the external gas has an
/// atmospheric relative-humidity coordinate.
///
/// # Errors
/// Ordinary assembly refusals plus incoherent/nonphysical supplied gas states.
pub fn realize_assembly_with_gas(
    assembly: &AcousticAssembly,
    gas: &GasState,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    realize_assembly_with_options(
        assembly,
        &RealizeOptions {
            gas: Some(*gas),
            ..RealizeOptions::default()
        },
    )
}

/// Realize a caller-resolved gas with an explicit observer constitutive law.
/// NASA-9 producers use this path: ISO applicability is a caller assertion,
/// never inferred from numerical equality with a legacy atmosphere. Invalid
/// declared ISO conditions refuse instead of selecting a different law.
pub fn realize_assembly_with_gas_and_absorption(
    assembly: &AcousticAssembly,
    gas: &GasState,
    absorption: crate::air_path::GasAbsorptionModel,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    realize_assembly_with_options(
        assembly,
        &RealizeOptions {
            gas: Some(*gas),
            observer_absorption: Some(absorption),
            ..RealizeOptions::default()
        },
    )
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
    realize_assembly_with_plate_chart_and_options(
        assembly,
        chart,
        options,
        &RealizeOptions::default(),
    )
}

/// Plate-chart variant of realize_assembly_with_gas. Geometry/material
/// charting and acoustic-medium selection remain independent inputs.
///
/// # Errors
/// Conflicting plate descriptions, invalid gas/chart controls, or ordinary
/// assembly admission and stepping errors refuse.
pub fn realize_assembly_with_plate_chart_and_gas(
    assembly: &AcousticAssembly,
    chart: &fs_plate::PlateChart,
    options: &PlateChartRadiation,
    gas: &GasState,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    realize_assembly_with_plate_chart_and_options(
        assembly,
        chart,
        options,
        &RealizeOptions {
            gas: Some(*gas),
            ..RealizeOptions::default()
        },
    )
}

/// Realize a plate chart with the same medium, absorption and bore controls as
/// an ordinary assembly. Regional material assignment does not select or
/// discard any gas constitutive law or numerical-resolution request.
///
/// # Errors
/// Conflicting plate descriptions, incoherent gas, invalid observer laws, and
/// ordinary chart/assembly/numerical refusals propagate before publication.
pub fn realize_assembly_with_plate_chart_and_options(
    assembly: &AcousticAssembly,
    chart: &fs_plate::PlateChart,
    radiation: &PlateChartRadiation,
    options: &RealizeOptions,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    if assembly.plate.is_some() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "chart assembly requires the uniform plate description to be absent",
        });
    }
    realize_assembly_inner(
        assembly,
        Some((AssemblyPlateChart::Authored(chart), radiation)),
        options,
    )
}

/// Realize a source-bound regional plate without erasing its selected loss law.
///
/// Elastic modes, force ports and any material bending loss are reduced once
/// by the existing material-chart owner. The resulting bodies join the same
/// string/duct/cavity coupling, gas loading and pressure observer as a bare
/// chart. Source-band admission and the prohibition on double-counted damping
/// survive this boundary. The retained regional loss is the existing diagonal
/// weak-loss image; this does not introduce full nonproportional modal damping.
///
/// # Errors
/// Conflicting plate descriptions, material-loss/source-domain refusals and
/// ordinary chart, medium, numerical or runtime refusals.
pub fn realize_assembly_with_material_plate_chart_and_options(
    assembly: &AcousticAssembly,
    specimen: &ResolvedPlateChart,
    radiation: &PlateChartRadiation,
    options: &RealizeOptions,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    if assembly.plate.is_some() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "chart assembly requires the uniform plate description to be absent",
        });
    }
    realize_assembly_inner(
        assembly,
        Some((AssemblyPlateChart::Material(specimen), radiation)),
        options,
    )
}

enum AssemblyPlateChart<'a> {
    Authored(&'a fs_plate::PlateChart),
    Material(&'a ResolvedPlateChart),
}

#[allow(clippy::too_many_lines)] // one coherent realization dispatcher
fn realize_assembly_inner(
    assembly: &AcousticAssembly,
    chart: Option<(AssemblyPlateChart<'_>, &PlateChartRadiation)>,
    realize_options: &RealizeOptions,
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
    if realize_options.bore != BoreRealization::default()
        && (assembly.reed.is_none() || assembly.duct.is_none() || assembly.string.is_some())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "explicit bore discretization requires a reed-driven duct without a string member",
        });
    }
    let gas = match realize_options.gas {
        Some(gas) => {
            validate_external_gas(&gas)?;
            gas
        }
        None => gas_state(assembly.ambient)?,
    };
    // Constitutive selection follows the input path, not equality of floating
    // coefficients with a freshly resolved air preset. This also lets sourced
    // plates, strings and reeds carry the same explicit observer law.
    let absorption = realize_options.observer_absorption.or_else(|| {
        realize_options
            .gas
            .map(|_| crate::air_path::GasAbsorptionModel::StokesKirchhoff)
    });
    match (absorption, realize_options.observer_relaxation) {
        (Some(crate::air_path::GasAbsorptionModel::MolecularRelaxation), relaxation) => {
            admit_relaxation(relaxation, &gas)?;
        }
        (_, Some(_)) => {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "an observer relaxation state needs the MolecularRelaxation observer law",
            });
        }
        (Some(model), None) => crate::air_path::absorb_pressure_history_with_model(
            &mut [],
            1.0 / f64::from(assembly.sample_rate_hz),
            assembly.listener.distance_m,
            &gas,
            model,
        )?,
        (None, None) => {}
    }
    let absorption =
        ObserverLaw { model: absorption, relaxation: realize_options.observer_relaxation };
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
        bodies.linear.extend(match chart {
            AssemblyPlateChart::Authored(chart) => certified_chart_radiators(chart, options)?,
            AssemblyPlateChart::Material(specimen) => {
                certified_material_chart_radiators(specimen, options)?
            }
        });
    }
    bodies.declare_radiation_loads(&gas, assembly.sample_rate_hz)?;
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
        return finish_observer(assembly, &gas, pressure_pa, absorption);
    }
    if assembly.string.is_some() && assembly.duct.is_some() {
        if assembly
            .duct
            .as_ref()
            .is_some_and(|d| bore_spec(d).is_some())
        {
            return realize_coupled_ode(assembly, &gas, &mut bodies, n, absorption);
        }
        return realize_coupled(assembly, &gas, &mut bodies, n, absorption);
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
            return finish_observer(assembly, &gas, pressure_pa, absorption);
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
    let mut bore_receipt = None;
    if let Some(duct) = assembly.duct.as_ref() {
        let hist = if let Some(reed) = assembly.reed {
            let (hist, receipt) = realize_reed_on_duct(
                duct,
                reed,
                &mut bodies,
                &gas,
                assembly.listener.distance_m,
                assembly.sample_rate_hz,
                n,
                realize_options.bore,
            )?;
            bore_receipt = Some(receipt);
            hist
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
    let mut realized = finish_observer(assembly, &gas, pressure_pa, absorption)?;
    realized.bore_receipt = bore_receipt;
    Ok(realized)
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
    absorption: ObserverLaw<'_>,
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
        duct_outlet(duct.termination),
        inlets,
        &taps,
        gas,
        duct.wall,
    )?;
    let (body, mut x_body, plate_observation) = compose_plate_chain(
        plates, chain, usize::from(inlets == 2), gas, assembly.sample_rate_hz,
    )?;
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    let listener_m = assembly.listener.distance_m;
    let mut texture = TextureDrive::try_new(assembly.contact_texture)?;
    let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
    let n_plate = plate_observation.as_ref().map_or(0, PlateBankObservation::state_dim);
    let body_linear = prepare_render_linear(&body)?;
    let body_exact = prepare_chain_exact(&body, &body_linear, dt)?;
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
                &body_linear,
                body_exact.as_ref(),
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
                + plate_chain_radiation(plate_observation.as_ref(), &rec.x, &x_body, gas, listener_m, dt)?
                + chain_inlet_pressure(&body, &rec.x, n_plate);
            x_body = rec.x;
        }
        return finish_observer(assembly, gas, out, absorption);
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
            &body_linear,
            body_exact.as_ref(),
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
            + plate_chain_radiation(plate_observation.as_ref(), &rec.x, &x_body, gas, listener_m, dt)?
            + chain_inlet_pressure(&body, &rec.x, n_plate);
        x_body = rec.x;
    }
    finish_observer(assembly, gas, out, absorption)
}

#[allow(clippy::too_many_arguments)] // one coherent chain step record
fn step_string_chain(
    body: &fs_phs::PortHamiltonian,
    linear: &fs_phs::linear_exact::LinearSplit,
    exact: Option<&fs_phs::linear_exact::PreparedExactZoh<'_>>,
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
        return step_chain_inlet(body, linear, exact, x, reed, blow, gas, t, dt);
    }
    step_plate_chain(body, linear, exact, x, f_bridge, reed, blow, gas, t, dt)
}

fn chain_inlet_pressure(body: &fs_phs::PortHamiltonian, x: &[f64], n_plate: usize) -> f64 {
    if n_plate != 0 {
        return 0.0;
    }
    // Structural port slot: a shorter output vector is an internal
    // invariant break, never physics — refuse loudly instead of zeroing.
    body.output(x)[0]
}

#[allow(clippy::too_many_arguments)] // one coherent chain step record
fn step_chain_inlet(
    body: &fs_phs::PortHamiltonian,
    linear: &fs_phs::linear_exact::LinearSplit,
    exact: Option<&fs_phs::linear_exact::PreparedExactZoh<'_>>,
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
    step_chain_flow(body, linear, exact, x, &u, dt)
}

#[allow(clippy::too_many_arguments)] // one coherent chain step record
fn step_plate_chain(
    body: &fs_phs::PortHamiltonian,
    linear: &fs_phs::linear_exact::LinearSplit,
    exact: Option<&fs_phs::linear_exact::PreparedExactZoh<'_>>,
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
        return step_chain_flow(body, linear, exact, x, &u, dt);
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
    step_chain_flow(body, linear, exact, x, &u, dt)
}

// Reuse the same complete quadratic propagator and independently integrated
// heat/work for every sample. Nonlinear storage keeps the existing shared
// Gonzalez owner and its explicit temporal admission.
fn prepare_chain_exact<'a>(
    system: &'a fs_phs::PortHamiltonian,
    linear: &fs_phs::linear_exact::LinearSplit,
    dt: f64,
) -> Result<Option<fs_phs::linear_exact::PreparedExactZoh<'a>>, AcousticRealizeError> {
    if !linear.is_exact_quadratic() { return Ok(None); }
    fs_phs::linear_exact::PreparedExactZoh::new(system, dt,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: system.state_dim(), max_omega_dt: core::f64::consts::PI,
        })
        .map(Some).map_err(render_phs_error)
}

fn step_chain_flow(
    system: &fs_phs::PortHamiltonian,
    linear: &fs_phs::linear_exact::LinearSplit,
    exact: Option<&fs_phs::linear_exact::PreparedExactZoh<'_>>,
    state: &[f64],
    input: &[f64],
    dt: f64,
) -> Result<fs_phs::StepRecord, AcousticRealizeError> {
    if let Some(exact) = exact {
        if !core::ptr::eq(exact.system(),system) || exact.duration().to_bits()!=dt.to_bits() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "prepared duct clock and complete physical system must match",
            });
        }
        exact.step(state,input).map_err(render_phs_error)
    } else {
        fs_phs::linear_exact::step_render(system,linear,state,input,dt).map_err(render_phs_error)
    }
}

fn prepare_render_linear(
    sys: &fs_phs::PortHamiltonian,
) -> Result<fs_phs::linear_exact::LinearSplit, AcousticRealizeError> {
    fs_phs::linear_exact::split_quadratic(
        sys,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: sys.state_dim(),
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .map_err(render_phs_error)
}

// The complete original structural/material/cavity/exterior system occupies
// the leading block. The duct starts after that COMPLETE extent, not after an
// inferred number of q/p pairs. Its positive inlet flow compresses the fluid;
// outward surface flow therefore enters with the opposite sign, while pressure
// acts on the same area map. External point force remains its independent port.
fn compose_plate_chain(
    plates: &mut PlateBank,
    chain: fs_phs::PortHamiltonian,
    face: usize,
    gas: &GasState,
    sample_rate_hz: u32,
) -> Result<(fs_phs::PortHamiltonian, Vec<f64>, Option<PlateBankObservation>), AcousticRealizeError> {
    let Some(image) = plates.physical_image(gas, sample_rate_hz, &fs_exec::CancelGate::default())? else {
        let state = vec![0.0; chain.state_dim()];
        return Ok((chain, state, None));
    };
    let (plate, mut state, observation) = image.into_parts();
    let system = transformer(plate, chain, 1, face, 1.0)
        .map_err(render_phs_error)?;
    state.resize(system.state_dim(), 0.0);
    Ok((system, state, Some(observation)))
}

fn plate_chain_radiation(
    observation: Option<&PlateBankObservation>,
    x1: &[f64],
    x0: &[f64],
    gas: &GasState,
    listener_m: f64,
    dt: f64,
) -> Result<f64, AcousticRealizeError> {
    let Some(observation) = observation else { return Ok(0.0); };
    let extent = observation.state_dim();
    let before = x0.get(..extent).ok_or(AcousticRealizeError::InvalidDescription {
        what: "coupled duct lost the original plate state prefix",
    })?;
    let after = x1.get(..extent).ok_or(AcousticRealizeError::InvalidDescription {
        what: "coupled duct lost the original plate state prefix",
    })?;
    observation.pressure_between(before, after, dt, gas.density, listener_m)
}

#[cfg(test)]
mod plate_chain_tests;
#[cfg(test)]
mod blown_ode_tests;

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
    absorption: ObserverLaw<'_>,
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
        .inlet_radius();
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
        return realize_coupled_kc(
            assembly,
            gas,
            plates,
            n,
            fitted,
            zc,
            area,
            &mut texture,
            absorption,
        );
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
    finish_observer(assembly, gas, out, absorption)
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
    absorption: ObserverLaw<'_>,
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
    finish_observer(assembly, gas, out, absorption)
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
        crate::driving_point::DrivingPointError::LadderUnconverged { .. } => {
            AcousticRealizeError::Reed {
                what: "characteristic line ladder reached max_taps ceiling before convergence",
            }
        }
        crate::driving_point::DrivingPointError::Duct(d) => AcousticRealizeError::Duct(d),
        crate::driving_point::DrivingPointError::Material(e) => AcousticRealizeError::Ambient(e.to_string()),
        crate::driving_point::DrivingPointError::Realize(_) => AcousticRealizeError::Reed {
            what: "characteristic realization refused",
        },
        crate::driving_point::DrivingPointError::Discrete(e) => AcousticRealizeError::Nonlinear(
            format!("characteristic-line discretization refused: {e}"),
        ),
    }
}

pub(crate) fn validate_external_gas(gas: &GasState) -> Result<(), AcousticRealizeError> {
    // Use the same EOS/caloric/transport admission as heat and flow. Checking
    // only rho*c and Pr would admit a pressure-only edit with stale density,
    // or a detuned sound speed unrelated to the declared gas composition.
    gas.validate_coherence().map_err(|error| AcousticRealizeError::Ambient(error.to_string()))
}

pub(crate) fn gas_state(ambient: AmbientGas) -> Result<GasState, AcousticRealizeError> {
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

/// The relaxation state that molecular relaxation needs, evaluated at exactly
/// `gas`.
fn admit_relaxation<'a>(
    relaxation: Option<&'a fs_material::relaxation::RelaxationState>,
    gas: &GasState,
) -> Result<&'a fs_material::relaxation::RelaxationState, AcousticRealizeError> {
    let state = relaxation.ok_or(AcousticRealizeError::InvalidDescription {
        what: "the MolecularRelaxation observer law needs RealizeOptions::observer_relaxation",
    })?;
    if state.gas() != gas {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "the observer relaxation state was evaluated for a different gas state",
        });
    }
    Ok(state)
}

fn finish_observer(
    assembly: &AcousticAssembly,
    gas: &GasState,
    mut pressure_pa: Vec<f64>,
    absorption: ObserverLaw<'_>,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    use crate::air_path::{GasAbsorptionModel, ObserverAbsorptionReceipt, relaxation_omission};
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    // The law actually applied. Only the scenario-air path arrives without one
    // (supplied gases and explicitly chosen laws were admitted before
    // stepping). Its legacy moist air is terrestrial air by construction, so it
    // makes its own admission, bound to the state. The law is ISO 9613 when ISO
    // admits the conditions, and otherwise the documented Stokes-Kirchhoff
    // fallback. Stokes-Kirchhoff leaves molecular relaxation out, which is
    // never silent: the receipt carries the labelled figure over the band.
    let law = absorption.model.unwrap_or_else(|| {
        let ambient = assembly.ambient;
        let relative_humidity = ambient.relative_humidity;
        match fs_material::iso9613::TerrestrialAir::moist_air(
            ambient.temperature_k,
            ambient.pressure_pa,
            relative_humidity,
        ) {
            Ok((_, air))
                if fs_material::iso9613::iso9613_absorption(&air, gas, relative_humidity, 1.0)
                    .is_ok() =>
            {
                GasAbsorptionModel::Iso9613Air { relative_humidity, air }
            }
            _ => GasAbsorptionModel::StokesKirchhoff,
        }
    });
    let relaxation_omitted = match law {
        GasAbsorptionModel::Iso9613Air { .. } | GasAbsorptionModel::MolecularRelaxation => None,
        GasAbsorptionModel::StokesKirchhoff | GasAbsorptionModel::Lossless => Some(relaxation_omission(
            assembly.listener.distance_m,
            0.5 * f64::from(assembly.sample_rate_hz),
        )?),
    };
    if law == GasAbsorptionModel::MolecularRelaxation {
        // The bulk delay r / c_inf is returned, not applied: no observer law
        // here places the arrival.
        crate::air_path::propagate_pressure_history_relaxing(
            &mut pressure_pa,
            dt,
            assembly.listener.distance_m,
            admit_relaxation(absorption.relaxation, gas)?,
        )?;
    } else {
        crate::air_path::absorb_pressure_history_with_model(
            &mut pressure_pa,
            dt,
            assembly.listener.distance_m,
            gas,
            law,
        )?;
    }
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa,
        gas: *gas,
        bore_receipt: None,
        observer_absorption: ObserverAbsorptionReceipt { law, relaxation_omitted },
    })
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

pub(crate) fn mode_zeta(
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
    let mut loss = intrinsic_mode_loss_parts(string, omega, wave_number)?;
    // A supplied Rayleigh pair is already a complete reduced loss model.
    if string.rayleigh.is_some() {
        return Ok(loss);
    }
    // The compact string path interprets width as circular diameter for air
    // drag. R [N s/m²] / (2 mu_linear omega) is dimensionless. Material
    // bindings derive this diameter and mu_linear from the same cross-section.
    let resistance = crate::air_path::oscillating_cylinder_air_resistance_per_length(
        0.5 * string.width_m,
        omega,
        gas,
    )?;
    loss.air = 0.5 * resistance / string.lin_density_kg_m / omega;
    if !loss.air.is_finite() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "cylinder air damping ratio is unrepresentable",
        });
    }
    Ok(loss)
}

/// The string's own loss, excluding its surrounding fluid. A causal fluid
/// realization must replace only the air term, keeping tension and material
/// bending losses in the actual (possibly fluid-loaded) modal mass basis.
pub(crate) fn intrinsic_string_mode_zeta(
    string: &PrestressedString,
    omega: f64,
    wave_number: f64,
) -> Result<f64, AcousticRealizeError> {
    Ok(intrinsic_mode_loss_parts(string, omega, wave_number)?.total())
}

fn intrinsic_mode_loss_parts(
    string: &PrestressedString,
    omega: f64,
    wave_number: f64,
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
    let tension = tension_loss_zeta(string, omega, wave_number)?;
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
            // zeta_k = c_k/(2 mu omega_k). This bending law does not dissipate
            // tension energy: eta*omega/(2E) would omit the bending dilution.
            // Any independently declared tension loss is added separately.
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
    if let Some(bending) = material_bending {
        return Ok(StringModeLossParts {
            material: bending + tension,
            ..StringModeLossParts::default()
        });
    }
    if string.relaxation_bending.is_some() {
        // The memory arms supply their own bending storage and dissipation in
        // pHS; the tension's loss is not theirs and stays modal.
        return Ok(StringModeLossParts {
            material: tension,
            ..StringModeLossParts::default()
        });
    }
    // A caller-authored modal ratio is an explicit reduced model. Elastic EI
    // supplies storage, not an additional dissipation law, and one ratio does
    // not identify a relaxation spectrum. Physical bending loss must select
    // the Kelvin-Voigt or causal memory path above.
    Ok(StringModeLossParts {
        authored: string.damping_ratio,
        material: tension,
        ..StringModeLossParts::default()
    })
}

/// Modal damping ratio of the tension's own loss factor eta_F.
///
/// The constant internal-friction term of the string loss model of Valette
/// and Cuesta, Mécanique de la corde vibrante (1993):
/// Q_F^-1 = eta_F T k² / (T k² + EI k⁴). A complex tension T (1 + i eta_F)
/// loses eta_F of the tension's share of the mode's stiffness per radian, so
/// zeta = eta_F T k² / (2 mu omega²) at the mode's own admitted frequency,
/// held as a modal viscous rate. It belongs to the string's construction
/// (core and winding), not to its bending law, and is material heat.
/// `None` contributes nothing. An authored ratio or Rayleigh pair is already a
/// complete reduced loss, so it cannot be combined with a physical eta_F.
fn tension_loss_zeta(
    string: &PrestressedString,
    omega: f64,
    wave_number: f64,
) -> Result<f64, AcousticRealizeError> {
    let Some(eta) = string.tension_loss_factor else {
        return Ok(0.0);
    };
    if !eta.is_finite() || eta < 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "tension loss factor must be finite and nonnegative",
        });
    }
    if string.rayleigh.is_some() || string.damping_ratio != 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "tension loss factor cannot be combined with authored internal/Rayleigh loss",
        });
    }
    let zeta = 0.5 * eta * string.tension_n * wave_number * wave_number
        / (string.lin_density_kg_m * omega * omega);
    if !zeta.is_finite() || zeta < 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "tension loss damping is unrepresentable",
        });
    }
    Ok(zeta)
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
                duct_outlet(line.termination),
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
                duct_outlet(line.termination),
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
    let descriptor_rate = sys
        .as_ref()
        .map(fs_phs::linear_exact::descriptor_pole_rate_bound)
        .transpose()
        .map_err(render_phs_error)?;
    if let Some(rate) = descriptor_rate {
        fs_phs::linear_exact::admit_midpoint_band(rate, dt)
            .map_err(AcousticRealizeError::TimeResolution)?;
    }
    let waveguide_linear = waveguide.as_ref().map(prepare_render_linear).transpose()?;
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
            fs_phs::linear_exact::step_render(
                wg,
                waveguide_linear.as_ref().expect("prepared waveguide clock"),
                &x,
                &u,
                dt,
            )
            .map_err(render_phs_error)?
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
    // One owner (DEF-09): the declared ends set both end corrections and the
    // outer radiation load; the Bessel neck loss joins it at the resonance.
    let resistance = crate::cavity_modal::resonance_neck_resistance(&cavity, gas)?;
    fs_phs::helmholtz_resonator_flow(
        cavity.volume_m3,
        cavity.neck_radius_m,
        cavity.neck_length_m,
        gas.density,
        gas.sound_speed,
        resistance,
        crate::cavity_modal::neck_ends(&cavity),
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
/// induced by transverse motion. The surrounding gas can be rebound without
/// resetting vibration; listener distance stays fixed;
/// no atmospheric path filter, common propagation delay or radiation loading is
/// included. The observer is a diagnostic, not an additional energy debit.
///
/// Ordinary rebinding keeps the exact length, linear density and sine-mode
/// count, hence the same mass-normalized basis. Explicit material substitution
/// at fixed geometry transfers physical displacement and conjugate momentum to
/// the new mass normalization and reports the external parameter work. Changing
/// geometry/supports, Prony memory, nonlinear response or phase requires a
/// different state-transfer owner.
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

/// Prescribed transverse force at a fixed material station for one sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringPointForce {
    /// Fraction of loaded string length, in `[0, 1]`.
    pub station_frac: f64,
    /// Signed force in the string's transverse vibration direction [N].
    pub force_n: f64,
}

/// Physical motion at a station, reconstructed from the accepted retained modes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringPointMotion {
    /// Accepted mechanical/material epoch of this observation.
    pub epoch: u64,
    /// Transverse displacement from the pre-stressed reference line [m].
    pub displacement_m: f64,
    /// Transverse velocity [m/s].
    pub velocity_m_s: f64,
}

/// Energy change when replacing a material at fixed physical displacement and
/// conjugate momentum. Ordinary rebinding also keeps the mass normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringMaterialUpdate {
    /// Shared material/operator/observer epoch after publication.
    pub epoch: u64,
    /// Modal mechanical energy before the update [J].
    pub energy_before_j: f64,
    /// Modal mechanical energy after the update [J].
    pub energy_after_j: f64,
    /// Signed work into the retained vibration. With unchanged mass this is
    /// `sum((omega_new²-omega_old²) q²/2)`; material substitution also includes
    /// the change in kinetic energy at fixed physical conjugate momentum.
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
        let gas = gas_state(ambient)?;
        Self::try_new_with_gas(cx, specimen, pluck, &gas, listener_m, sample_rate_hz)
    }

    /// Initialize in one complete caller-resolved gas or mixture.
    ///
    /// The material owner checks EOS, caloric, transport and acoustic coherence.
    /// Drag and the compact observer use this SAME snapshot; no air preset or
    /// atmospheric absorption law is inferred from its numerical coefficients.
    /// Retain the originating material model/receipts when changing this gas.
    pub fn try_new_with_gas(
        cx: &fs_exec::Cx<'_>,
        specimen: crate::string_specimen::ResolvedStringSpecimen,
        pluck: Option<Pluck>,
        gas: &GasState,
        listener_m: f64,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        validate_external_gas(gas)?;
        let gas = *gas;
        let string = specimen.string();
        validate_incremental_string(&string)?;
        validate_string(&string, pluck, None)?;
        if !(listener_m.is_finite() && listener_m > 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "incremental string listener distance must be finite and positive",
            });
        }
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

    /// Current states in the accepted material's mass-normalized sine basis.
    #[must_use]
    pub fn states(&self) -> &[ModalAcousticState] {
        self.member.model.states()
    }

    /// Frequencies and complete admitted damping at the same epoch.
    #[must_use]
    pub fn modes(&self) -> &[ModalAcousticMode] {
        self.member.model.modes()
    }

    /// Accepted surrounding gas used by air drag and the compact observer.
    #[must_use]
    pub const fn ambient(&self) -> &GasState {
        &self.gas
    }

    /// Project physical point forces into this epoch's retained modal basis.
    ///
    /// With `phi_k = sin(k pi x/L)/sqrt(mu L/2)`, generalized force is
    /// `sum(phi_k F)` and point velocity is `sum(phi_k qdot_k)`. Using the same
    /// map in both directions preserves instantaneous power and held-load work.
    /// Pass the returned `[N/sqrt(kg)]` vector to the existing mechanical or
    /// thermal step before any rebind/refinement; do not cache it across epochs.
    /// Loads are held at fixed stations during the step. End-support loads do
    /// no work in this pinned model; support reaction dynamics are not modeled.
    /// This projection supplies no contact law or opposite-body reaction.
    pub fn project_point_forces(
        &self,
        cx: &fs_exec::Cx<'_>,
        loads: &[StringPointForce],
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        let mut forces = vec![0.0; self.states().len()];
        let string = self.specimen.string();
        let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
        for load in loads {
            cx.checkpoint()
                .map_err(|_| AcousticRealizeError::Cancelled)?;
            validate_string_station(load.station_frac)?;
            if !load.force_n.is_finite() {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "point force must be finite newtons",
                });
            }
            for (index, force) in forces.iter_mut().enumerate() {
                if index % 64 == 0 {
                    cx.checkpoint()
                        .map_err(|_| AcousticRealizeError::Cancelled)?;
                }
                *force +=
                    string_point_participation(index, load.station_frac, mass_scale) * load.force_n;
                if !force.is_finite() {
                    return Err(AcousticRealizeError::InvalidDescription {
                        what: "projected string point force is not finite",
                    });
                }
            }
        }
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        Ok(forces)
    }

    /// Initialize an unexcited string in equilibrium under settled point loads.
    ///
    /// The loading history is assumed to precede the simulation window: this
    /// establishes stored elastic energy, not a measured actuator-work transfer.
    /// Continue stepping the same loads to hold the shape, or release them to
    /// excite vibration. Only a zero-state runtime with no accepted samples or
    /// observer history can be initialized; existing motion is never reset.
    /// Material state and sample count are unchanged. Publication advances the
    /// epoch and records zero previous acceleration for the stationary history.
    /// Cancellation or a force/state/energy refusal leaves the runtime untouched.
    pub fn initialize_point_load_equilibrium(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        loads: &[StringPointForce],
    ) -> Result<u64, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        if self.accepted_samples != 0
            || self.member.observer.previous_acceleration.is_some()
            || self.states().iter().any(|state| {
                state.displacement_m_sqrt_kg != 0.0 || state.velocity_m_sqrt_kg_per_s != 0.0
            })
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "settled point-load initialization requires an unexcited string before stepping",
            });
        }
        let epoch = next_string_epoch(self.epoch)?;
        let forces = self.project_point_forces(cx, loads)?;
        let mut candidate = self.member.clone();
        candidate
            .model
            .initialize_static_equilibrium(&forces)
            .map_err(AcousticRealizeError::Modal)?;
        candidate.observer.previous_acceleration = Some(vec![0.0; forces.len()]);
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        self.member = candidate;
        self.epoch = epoch;
        Ok(epoch)
    }

    /// Observe physical displacement and conjugate velocity at a fixed station.
    /// Reads accepted state only; it does not advance time or alter acoustics.
    pub fn point_motion(
        &self,
        cx: &fs_exec::Cx<'_>,
        station_frac: f64,
    ) -> Result<StringPointMotion, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        validate_string_station(station_frac)?;
        let string = self.specimen.string();
        let mass_scale = det::sqrt(string.lin_density_kg_m * string.length_m / 2.0);
        let mut motion = StringPointMotion {
            epoch: self.epoch,
            displacement_m: 0.0,
            velocity_m_s: 0.0,
        };
        for (index, state) in self.states().iter().enumerate() {
            if index % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| AcousticRealizeError::Cancelled)?;
            }
            let phi = string_point_participation(index, station_frac, mass_scale);
            motion.displacement_m += phi * state.displacement_m_sqrt_kg;
            motion.velocity_m_s += phi * state.velocity_m_sqrt_kg_per_s;
        }
        if !motion.displacement_m.is_finite() || !motion.velocity_m_s.is_finite() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "reconstructed string point motion is not finite",
            });
        }
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        Ok(motion)
    }

    /// Replace surrounding gas without restarting vibration or advancing time.
    /// Air damping and acoustic observation rebuild through their shared owners;
    /// material receipts, modal coordinates and previous acceleration survive.
    /// The update advances the epoch, but consumes no sample or material heat.
    /// This prescribed gas change does not model fluid-loading work, boundary
    /// heat transfer or a propagating change of the surrounding atmosphere.
    /// A refusal or cancellation preserves the complete accepted runtime.
    pub fn rebind_ambient(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        ambient: AmbientGas,
    ) -> Result<StringMaterialUpdate, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        self.rebind_gas(cx, &gas_state(ambient)?)
    }

    /// Replace the complete surrounding gas through its shared material owner.
    /// This has the same state-preserving, prescribed-medium semantics as
    /// `rebind_ambient`, and admits source-backed gases and mixtures.
    /// No thermal energy is inferred from this instantaneous parameter update.
    pub fn rebind_gas(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        gas: &GasState,
    ) -> Result<StringMaterialUpdate, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        validate_external_gas(gas)?;
        let mut candidate = self.clone();
        candidate.gas = *gas;
        let update = candidate.rebind(cx, self.specimen.clone(), 0.0)?;
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        *self = candidate;
        Ok(update)
    }

    /// Add higher sine modes at rest without advancing the physical clock.
    /// Existing coordinates and accepted acceleration history are retained
    /// exactly. New modes have zero prior acceleration (no prior forcing).
    /// This is explicit fixed-basis refinement, not reconstruction of omitted
    /// history, a tail-error estimate or transfer across changed geometry.
    /// Returns the new epoch; smaller/equal counts or inadmissible frequencies
    /// refuse without changing the accepted runtime.
    pub fn refine_modes(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        mode_count: usize,
    ) -> Result<u64, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        let epoch = next_string_epoch(self.epoch)?;
        let old_count = self.states().len();
        if mode_count <= old_count {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "mode refinement requires a strictly larger retained sine basis",
            });
        }
        let specimen = self.specimen.clone().with_retained_modes(mode_count);
        let string = specimen.string();
        validate_incremental_string(&string)?;
        let mut member = linear_string_member(
            &string,
            None,
            None,
            &self.gas,
            self.listener_m,
            self.sample_rate_hz,
            1.0,
        )?;
        let mut states = member.model.states().to_vec();
        states[..old_count].clone_from_slice(self.states());
        member
            .model
            .restore_states(&states)
            .map_err(AcousticRealizeError::Modal)?;
        member.observer.previous_acceleration = self
            .member
            .observer
            .previous_acceleration
            .as_ref()
            .map(|previous| {
                let mut extended = previous.clone();
                extended.resize(mode_count, 0.0);
                extended
            });
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        self.specimen = specimen;
        self.member = member;
        self.epoch = epoch;
        Ok(epoch)
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
        self.replace_material_state(cx, specimen, maximum_abs_parameter_work_j, false)
    }

    /// Substitute a source-bound material without restarting the vibrating string.
    ///
    /// Length, circular radius, pinned supports and retained sine modes must be
    /// unchanged. The new specimen supplies its own density, elasticity,
    /// prestress and admitted loss law. No material name selects solver code.
    /// The unchanged sine shapes allow an explicit canonical transfer: with
    /// `s = sqrt(mu_new / mu_old)`, `q_new = s q_old` and `p_new = p_old / s`.
    /// Physical displacement and transverse momentum therefore survive; velocity
    /// generally changes. This is not a velocity-preserving replacement.
    ///
    /// The caller imposes an external parameter intervention and owes the signed
    /// change in kinetic plus elastic energy returned as `parameter_work_j`.
    /// It is not chemical conversion, melting, mass transport or material heat.
    /// The axial prestress reservoir remains outside the retained transverse
    /// vibration. This operation is unavailable on the thermal wrapper, whose
    /// immutable material/enthalpy binding must not silently change identity.
    ///
    /// Force projection, damping and pressure weights rebuild together. Accepted
    /// physical acceleration history survives its change in normalization, so
    /// the next observer sample includes the change in endpoint acceleration.
    /// Distributional radiation from an instantaneous velocity jump is not
    /// reconstructed; that needs a time-resolved transition. Time does not
    /// advance and no crossfade is inserted. Reproject held physical forces
    /// after this operation; old mass-normalized force vectors are stale.
    ///
    /// # Errors
    /// Changed geometry/basis, unsupported loss or material domains, unrepresentable
    /// transferred state, exceeded parameter-work budget or cancellation refuse
    /// without publishing any material, motion, observer or clock change.
    pub fn substitute_material(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        specimen: crate::string_specimen::ResolvedStringSpecimen,
        maximum_abs_parameter_work_j: f64,
    ) -> Result<StringMaterialUpdate, AcousticRealizeError> {
        self.replace_material_state(cx, specimen, maximum_abs_parameter_work_j, true)
    }

    fn replace_material_state(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        specimen: crate::string_specimen::ResolvedStringSpecimen,
        maximum_abs_parameter_work_j: f64,
        substitute: bool,
    ) -> Result<StringMaterialUpdate, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        let epoch = next_string_epoch(self.epoch)?;
        let old = self.specimen.string();
        let new = specimen.string();
        validate_incremental_string(&new)?;
        let same_shapes =
            old.length_m.to_bits() == new.length_m.to_bits() && old.n_modes == new.n_modes;
        if substitute {
            if !same_shapes || self.specimen.radius_m().to_bits() != specimen.radius_m().to_bits() {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "material substitution requires the same length, circular radius and pinned sine basis",
                });
            }
        } else if !same_shapes
            || self.specimen.material().card_identity() != specimen.material().card_identity()
            || old.lin_density_kg_m.to_bits() != new.lin_density_kg_m.to_bits()
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
        // Take roots first so a representable normalization need not overflow
        // in the intermediate density ratio. Both must be strictly positive.
        let mass_scale = if old.lin_density_kg_m.to_bits() == new.lin_density_kg_m.to_bits() {
            1.0
        } else {
            det::sqrt(new.lin_density_kg_m) / det::sqrt(old.lin_density_kg_m)
        };
        if !(mass_scale.is_finite() && mass_scale > 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "material substitution mass normalization is unrepresentable",
            });
        }
        let transferred: Vec<_> = self
            .states()
            .iter()
            .map(|state| ModalAcousticState {
                displacement_m_sqrt_kg: state.displacement_m_sqrt_kg * mass_scale,
                velocity_m_sqrt_kg_per_s: state.velocity_m_sqrt_kg_per_s / mass_scale,
            })
            .collect();
        if transferred.iter().zip(self.states()).any(|(next, previous)| {
            (previous.displacement_m_sqrt_kg != 0.0 && next.displacement_m_sqrt_kg == 0.0)
                || (previous.velocity_m_sqrt_kg_per_s != 0.0
                    && next.velocity_m_sqrt_kg_per_s == 0.0)
        }) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "material substitution underflows nonzero motion",
            });
        }
        candidate
            .model
            .restore_states(&transferred)
            .map_err(AcousticRealizeError::Modal)?;
        candidate.observer.previous_acceleration = self
            .member
            .observer
            .previous_acceleration
            .as_ref()
            .map(|history| history.iter().map(|value| value * mass_scale).collect());
        if candidate
            .observer
            .previous_acceleration
            .as_ref()
            .zip(self.member.observer.previous_acceleration.as_ref())
            .is_some_and(|(history, previous)| {
                history.iter().zip(previous).any(|(value, before)| {
                    !value.is_finite() || (*before != 0.0 && *value == 0.0)
                })
            })
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "material substitution acceleration history is unrepresentable",
            });
        }
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
        // Difference kinetic and stiffness terms mode by mode; subtracting
        // total energies would erase small work against unchanged large terms.
        let parameter_work_j: f64 = self
            .modes()
            .iter()
            .zip(candidate.model.modes())
            .zip(self.states())
            .zip(&transferred)
            .map(|(((before, after), state), transferred)| {
                let old_wq = before.angular_frequency_rad_s * state.displacement_m_sqrt_kg;
                let new_wq = after.angular_frequency_rad_s * transferred.displacement_m_sqrt_kg;
                let old_v = state.velocity_m_sqrt_kg_per_s;
                let new_v = transferred.velocity_m_sqrt_kg_per_s;
                0.5 * (new_wq - old_wq) * (new_wq + old_wq)
                    + 0.5 * (new_v - old_v) * (new_v + old_v)
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
    /// The declared reversible-heat bracket does not contain a consistent endpoint.
    ThermoelasticBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    /// The bounded thermal/mechanical energy solve exhausted its work allowance.
    ThermoelasticConvergence {
        iterations: usize,
        residual_j: f64,
    },
}

impl core::fmt::Display for ThermalStringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-THERMAL-STRING: {self:?}")
    }
}

impl std::error::Error for ThermalStringError {}

/// Uniform solid temperature coupled to material-resolved viscous vibration.
///
/// The ordinary constructor uses constant density/Young's modulus and fixed
/// tension. `try_new_thermoelastic` also admits sourced fixed-support expansion
/// and affine temperature-dependent stiffness, with axial/modal elastic storage
/// and its entropy contribution exchanged through the same thermal owner.
/// Viscosity may depend on temperature in both cases.
/// The admitted thermal chart has constant solid density and uses total
/// specific enthalpy; neither latent heat nor air loss is credited twice.
///
/// Each sample freezes viscosity at its initial temperature, evolves the exact
/// mechanical substep, deposits its material loss, and resolves viscosity at
/// the new temperature. This is first-order partitioned feedback, not an exact
/// continuous-temperature solution. The explicit temperature-increment budget
/// limits a step; time refinement is still required to assess coupling error.
/// Spatial gradients, moving geometry, melting and acoustic radiation loading
/// are outside this rung. Ambient gas remains independent of the solid temperature
/// unless explicitly coupled through the finite-gas transaction.
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
    thermoelastic: Option<thermoelastic::ThermoelasticString>,
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
    /// Change in fixed-support axial elastic free energy [J].
    pub axial_free_energy_change_j: f64,
    /// Change in `T*s_mechanical`, including axial and retained modal entropy [J].
    /// Together with free energy this forms the mechanical internal energy.
    pub mechanical_entropy_energy_change_j: f64,
    /// Elastic parameter work into the retained modes at their unchanged endpoint
    /// coordinates [J]. Included in the final modal energy, not added twice.
    pub modal_parameter_work_j: f64,
    /// Signed reversible heat supplied to the original caloric/transport owner [J],
    /// separate from the physical Kelvin-Voigt heat in `vibration.dissipation`.
    pub reversible_heat_j: f64,
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
        Self::try_new_impl(cx, mechanical, card, curve, specific_enthalpy_j_kg,
            maximum_temperature_increment_k, None)
    }

    #[allow(clippy::too_many_arguments)]
    fn try_new_impl(
        cx: &fs_exec::Cx<'_>,
        mechanical: LinearMaterialStringRuntime,
        card: fs_matdb::MaterialCard,
        curve: EquilibriumEnthalpyPhaseCurve,
        specific_enthalpy_j_kg: f64,
        maximum_temperature_increment_k: f64,
        thermoelastic: Option<thermoelastic::ThermoelasticString>,
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
            || (thermoelastic.is_none()
                && !matches!(specimen.prestress(), StringPrestress::FixedTension(_)))
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
            if !matches!(claim.value, fs_matdb::PropertyValue::Scalar { .. })
                && (name == "density" || thermoelastic.is_none()) {
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
            thermoelastic,
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

    /// Update prescribed surrounding gas at the expected accepted epoch.
    /// This changes drag/observation only; the material enthalpy is unchanged.
    /// Use the thermal transport step for boundary heat exchange.
    pub fn rebind_ambient(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        ambient: AmbientGas,
    ) -> Result<StringMaterialUpdate, ThermalStringError> {
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        self.mechanical
            .rebind_ambient(cx, ambient)
            .map_err(ThermalStringError::Acoustic)
    }

    /// Rebind a complete source-resolved gas at the expected accepted epoch.
    /// Solid enthalpy and modal coordinates survive unchanged. For finite gas
    /// heating and reciprocal drag-energy storage, use `step_with_isochoric_gas`.
    pub fn rebind_gas(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        gas: &GasState,
    ) -> Result<StringMaterialUpdate, ThermalStringError> {
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        self.mechanical
            .rebind_gas(cx, gas)
            .map_err(ThermalStringError::Acoustic)
    }

    /// Accepted uniform thermal state.
    ///
    /// Mode refinement and ambient rebinding do not alter this state.
    #[must_use]
    pub const fn thermal(&self) -> EquilibriumPhaseState {
        self.thermal
    }

    /// Refine the fixed sine basis while preserving material enthalpy.
    /// New modes begin at rest; no omitted vibration history is reconstructed.
    pub fn refine_modes(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        mode_count: usize,
    ) -> Result<u64, ThermalStringError> {
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        self.mechanical
            .refine_modes(cx, mode_count)
            .map_err(ThermalStringError::Acoustic)
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
        transport: impl FnMut(
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
        transport: impl FnMut(
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
        mut transport: impl FnMut(
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
        mut transport: impl FnMut(
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
        if let Some(law) = &self.thermoelastic {
            return law.specimen_at(self, thermal);
        }
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

    /// Advance mechanics, material heating and thermoelastic feedback as one commit.
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
        let proposed = self.propose_thermal_step(cx, &candidate, &vibration, |input| {
            let heat = input.internal_heat_j + external_heat_j;
            let thermal = input.curve.advance_specific_energy(input.initial, heat / input.mass_kg)
                .map_err(ThermalStringError::Phase)?;
            Ok::<_, ThermalStringError>(fs_material::phase::UniformEnthalpyStep {
                state: thermal,
                external_heat_j,
                energy_residual_tolerance_j: 0.0,
                report: (),
            })
        })?;
        let result = self.finish_thermal_step(
            cx,
            forces,
            candidate,
            vibration,
            proposed,
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
    /// Thermoelastic feedback can call it repeatedly from the SAME accepted
    /// state with different signed internal heat. Only the accepted proposal
    /// survives; the callback must be repeatable and retain its actual boundary
    /// transfer. Ordinary fixed-tension strings still call it once.
    pub fn step_with_thermal_transport<R, E: core::fmt::Display>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        mut transport: impl FnMut(
            fs_material::phase::UniformEnthalpyStepInput<'_>,
        ) -> Result<fs_material::phase::UniformEnthalpyStep<R>, E>,
    ) -> Result<ThermalStringTransportFrame<R>, ThermalStringError> {
        if expected_epoch != self.mechanical.epoch() {
            return Err(ThermalStringError::Epoch {
                expected: expected_epoch,
                actual: self.mechanical.epoch(),
            });
        }
        let mut candidate = self.mechanical.clone();
        let vibration = candidate
            .step(cx, forces)
            .map_err(ThermalStringError::Acoustic)?;
        let proposed = self.propose_thermal_step(cx, &candidate, &vibration, |input| {
            transport(input).map_err(|error| ThermalStringError::Transport(error.to_string()))
        })?;
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
        mut transport: impl FnMut(
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
        candidate: LinearMaterialStringRuntime,
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
        let (mut candidate, exchange) = self.prepare_thermal_rebind(cx, &candidate, thermal)?;
        // The accepted frame reports the new elastic law's actual energy.
        // Losses/work above still belong to the frozen-temperature mechanical
        // substep; parameter work is the separately balanced thermal exchange.
        if self.thermoelastic.is_some() {
            vibration.acoustic.modal_energy_j = candidate.modes().iter().zip(candidate.states())
                .map(|(mode, state)| 0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
                    + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2)))
                .collect();
            vibration.acoustic.total_modal_energy_j = vibration.acoustic.modal_energy_j.iter().sum();
        }
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
                + exchange.axial_free_energy_change_j
                + exchange.mechanical_entropy_energy_change_j
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
        let energy_roundoff_tolerance_j = d.roundoff_tolerance_j
            + exchange.roundoff_tolerance_j + 32.0 * f64::EPSILON * scale;
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
                axial_free_energy_change_j: exchange.axial_free_energy_change_j,
                mechanical_entropy_energy_change_j: exchange.mechanical_entropy_energy_change_j,
                modal_parameter_work_j: exchange.modal_parameter_work_j,
                reversible_heat_j: exchange.reversible_heat_j,
                energy_balance_residual_j,
                energy_roundoff_tolerance_j,
                thermal_solve_tolerance_j,
            },
            transport: report,
        })
    }
}

fn validate_string_station(station_frac: f64) -> Result<(), AcousticRealizeError> {
    if !station_frac.is_finite() || !(0.0..=1.0).contains(&station_frac) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "string station must be a finite fraction in [0, 1]",
        });
    }
    Ok(())
}

fn string_point_participation(index: usize, station_frac: f64, mass_scale: f64) -> f64 {
    // Preserve exact pinned endpoints instead of taking a rounded sin(k*pi).
    if station_frac == 0.0 || station_frac == 1.0 {
        0.0
    } else {
        det::sin((index + 1) as f64 * core::f64::consts::PI * station_frac) / mass_scale
    }
}

fn next_string_epoch(epoch: u64) -> Result<u64, AcousticRealizeError> {
    epoch
        .checked_add(1)
        .ok_or(AcousticRealizeError::InvalidDescription {
            what: "incremental string epoch or sample counter overflow",
        })
}

pub(crate) fn validate_incremental_string(
    string: &PrestressedString,
) -> Result<(), AcousticRealizeError> {
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
    linear: fs_phs::linear_exact::LinearSplit,
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
    let linear = fs_phs::linear_exact::split_quadratic(
        &sys,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: sys.state_dim(),
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .map_err(render_phs_error)?;
    let mut x = vec![0.0; sys.state_dim()];
    if let Some(pluck) = pluck {
        for k in 1..=string.n_modes {
            x[2 * (k - 1)] = mass_scale * triangular_pluck_modal(pluck, k);
        }
    }
    initialize_relaxed_bending(string, &mut x);
    Ok(KcMember {
        sys,
        linear,
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
    let rec = fs_phs::linear_exact::step_render(&member.sys, &member.linear, &member.x, &u, dt)
        .map_err(render_phs_error)?;
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

fn render_phs_error(error: fs_phs::PhsError) -> AcousticRealizeError {
    match error {
        fs_phs::PhsError::TimeResolution(band) => AcousticRealizeError::TimeResolution(band),
        other => AcousticRealizeError::Nonlinear(other.to_string()),
    }
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

#[allow(clippy::too_many_arguments)]
fn realize_reed_on_duct(
    duct: &ViscothermalDuct,
    reed: BeatingReed,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    bore: BoreRealization,
) -> Result<(Vec<f64>, BoreRealizationReceipt), AcousticRealizeError> {
    let physics = physics_duct(duct)?;
    let termination = match duct.termination {
        WaveguideEnd::Closed => Termination::Closed,
        WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
        WaveguideEnd::FlangedOpen => Termination::FlangedOpen,
    };
    refuse_open_nyquist(duct, gas, sample_rate_hz)?;
    match bore {
        BoreRealization::CharacteristicLine(ladder) => {
            let (samples, receipt) = crate::reed_bore::realize_reed_bore_with_ladder(
                &physics,
                gas,
                reed,
                termination,
                plates,
                listener_m,
                sample_rate_hz,
                n,
                wall_pin(duct.wall).as_ref(),
                Some(&ladder),
            )?;
            Ok((samples, BoreRealizationReceipt::CharacteristicLine(receipt)))
        }
        BoreRealization::LumpedLadder {
            cells,
            max_pitch_error_cents,
        } => {
            if cells < 2 {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "lumped ladder requires at least 2 cells",
                });
            }
            if !(max_pitch_error_cents.is_finite() && max_pitch_error_cents > 0.0) {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "lumped ladder pitch tolerance must be finite and positive",
                });
            }
            if plates.vk.is_some() {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "lumped ladder ODE reed requires a linear plate bank",
                });
            }
            let dispersion_error_cents = lumped_ladder_dispersion_cents(cells);
            let min_admissible_cells = min_admissible_lumped_cells(max_pitch_error_cents);
            if dispersion_error_cents > max_pitch_error_cents {
                return Err(AcousticRealizeError::LumpedLadderBand {
                    cells,
                    min_admissible_cells,
                    pitch_error_cents: dispersion_error_cents,
                    max_pitch_error_cents,
                });
            }
            let (sections, taps) = bore_spec_with_cells(duct, cells).ok_or(
                AcousticRealizeError::InvalidDescription {
                    what: "duct cannot be sliced into a lumped ladder",
                },
            )?;
            // Dispersion alone is a necessary condition. The admission is the
            // ladder's MEASURED passive pitch against the line's own TMM.
            let pitch_error_cents = lumped_ladder_pitch_error_cents(
                &physics,
                &sections,
                &taps,
                duct,
                termination,
                gas,
            )?;
            if pitch_error_cents.abs() > max_pitch_error_cents {
                return Err(AcousticRealizeError::LumpedLadderBand {
                    cells,
                    min_admissible_cells,
                    pitch_error_cents,
                    max_pitch_error_cents,
                });
            }
            let samples = realize_reed_ode(
                &sections,
                duct_outlet(duct.termination),
                &taps,
                reed,
                plates,
                gas,
                listener_m,
                sample_rate_hz,
                n,
                duct.wall,
            )?;
            Ok((
                samples,
                BoreRealizationReceipt::LumpedLadder {
                    cells,
                    min_admissible_cells,
                    dispersion_error_cents,
                    pitch_error_cents,
                    max_pitch_error_cents,
                },
            ))
        }
    }
}

/// Quarter-wave fundamental spatial dispersion error (in cents) of an `N`-cell
/// lumped LC ladder: `|1200 * log2(sin(pi / (2N)) / (pi / (2N)))|`.
#[must_use]
pub fn lumped_ladder_dispersion_cents(cells: usize) -> f64 {
    if cells == 0 {
        return f64::INFINITY;
    }
    let theta = core::f64::consts::PI / (2.0 * cells as f64);
    let ratio = det::sin(theta) / theta;
    (1200.0 * det::ln(ratio) / core::f64::consts::LN_2).abs()
}

/// Minimum cell count `N >= 2` whose quarter-wave fundamental spatial
/// dispersion error is at most `max_pitch_error_cents`.
#[must_use]
pub fn min_admissible_lumped_cells(max_pitch_error_cents: f64) -> usize {
    if !(max_pitch_error_cents.is_finite() && max_pitch_error_cents > 0.0) {
        return usize::MAX;
    }
    let mut cells = 2usize;
    while cells < 4096 && lumped_ladder_dispersion_cents(cells) > max_pitch_error_cents {
        cells += 1;
    }
    cells
}

/// Passive input impedance of a linear single-port bore chain,
/// `Z(iω) = Gᵀ Q (iω I − (J − R) Q)⁻¹ G`, read from its port-Hamiltonian
/// structure. `Q` comes column by column from the (linear) effort map.
struct ChainImpedance {
    n: usize,
    a: Vec<f64>,
    g: Vec<f64>,
    gq: Vec<f64>,
}

impl ChainImpedance {
    fn new(line: &fs_phs::PortHamiltonian) -> Result<Self, AcousticRealizeError> {
        let n = line.state_dim();
        if line.port_dim() != 1 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "passive ladder impedance needs exactly one inlet port",
            });
        }
        let (j, r, g) = line.structure();
        let mut q = vec![0.0; n * n];
        let mut unit = vec![0.0; n];
        for col in 0..n {
            unit[col] = 1.0;
            for (row, value) in line.effort(&unit).into_iter().enumerate() {
                q[row * n + col] = value;
            }
            unit[col] = 0.0;
        }
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for k in 0..n {
                a[i * n + k] = (0..n)
                    .map(|l| (j[i * n + l] - r[i * n + l]) * q[l * n + k])
                    .sum();
            }
        }
        let gq = (0..n)
            .map(|k| (0..n).map(|i| g[i] * q[i * n + k]).sum())
            .collect();
        Ok(Self {
            n,
            a,
            g: g.to_vec(),
            gq,
        })
    }

    fn magnitude(&self, omega: f64) -> Result<f64, AcousticRealizeError> {
        let n = self.n;
        let mut m: Vec<C64> = self.a.iter().map(|&v| C64::new(-v, 0.0)).collect();
        for i in 0..n {
            m[i * n + i] = m[i * n + i] + C64::new(0.0, omega);
        }
        let lu = fs_la::eigen_complex::lu_complex(&m, n).map_err(|_| {
            AcousticRealizeError::InvalidDescription {
                what: "lumped ladder chain is singular at a sampled frequency",
            }
        })?;
        let mut x: Vec<C64> = self.g.iter().map(|&v| C64::new(v, 0.0)).collect();
        lu.solve(&mut x);
        let z = self
            .gq
            .iter()
            .zip(&x)
            .fold(C64::ZERO, |acc, (&w, &xv)| acc + xv.scale(w));
        Ok(z.abs())
    }
}

/// Angular frequency of the lowest passive input-impedance resonance on
/// `[lo, hi]`: the first interior local maximum of `|Z|` on a `scan`-cell
/// uniform grid, refined by golden-section search across its two neighbouring
/// cells. `None` when the grid has no interior maximum: a window edge is never
/// a resonance.
fn first_impedance_peak(
    lo: f64,
    hi: f64,
    scan: usize,
    mut magnitude: impl FnMut(f64) -> Result<f64, AcousticRealizeError>,
) -> Result<Option<f64>, AcousticRealizeError> {
    let step = (hi - lo) / scan as f64;
    let (mut before, mut here) = (magnitude(lo)?, magnitude(lo + step)?);
    let mut peak = None;
    for i in 2..=scan {
        let after = magnitude(lo + step * i as f64)?;
        if here > before && here >= after {
            peak = Some(lo + step * (i - 1) as f64);
            break;
        }
        (before, here) = (here, after);
    }
    let Some(center) = peak else {
        return Ok(None);
    };
    let (mut a, mut b) = (center - step, center + step);
    let ratio = 0.5 * (det::sqrt(5.0) - 1.0);
    let (mut c, mut d) = (b - ratio * (b - a), a + ratio * (b - a));
    let (mut fc, mut fd) = (magnitude(c)?, magnitude(d)?);
    for _ in 0..48 {
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - ratio * (b - a);
            fc = magnitude(c)?;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + ratio * (b - a);
            fd = magnitude(d)?;
        }
    }
    Ok(Some(f64::midpoint(a, b)))
}

/// Signed model-form pitch error [cents] of the lumped ladder against the
/// characteristic line's own TMM (exact Bessel losses, the same termination
/// and wall), compared at the line's lowest passive input-impedance resonance
/// ([`lumped_ladder_resonances`]). Positive means the ladder is sharp; infinite
/// when the ladder has no resonance within ±386 cents of the line's.
fn lumped_ladder_pitch_error_cents(
    physics: &Duct,
    sections: &[AcousticSection],
    taps: &[AcousticTap],
    duct: &ViscothermalDuct,
    termination: Termination,
    gas: &GasState,
) -> Result<f64, AcousticRealizeError> {
    let (ladder, line) = lumped_ladder_resonances(physics, sections, taps, duct, termination, gas)?;
    Ok(ladder.map_or(f64::INFINITY, |ladder| {
        1200.0 * det::ln(ladder / line) / core::f64::consts::LN_2
    }))
}

/// Lowest passive input-impedance resonances `(ladder, line)` [rad/s]. The
/// line is scanned from a quarter of the closed-open quarter-wave frequency to
/// six times it, so a closed far end (half-wave), a cone or an open tone hole
/// moves the resonance instead of leaving the scan on a window edge. Each
/// ladder sample is a dense chain solve, so the ladder is scanned only within
/// ±386 cents of the line's resonance; `None` when it has none there.
fn lumped_ladder_resonances(
    physics: &Duct,
    sections: &[AcousticSection],
    taps: &[AcousticTap],
    duct: &ViscothermalDuct,
    termination: Termination,
    gas: &GasState,
) -> Result<(Option<f64>, f64), AcousticRealizeError> {
    let length: f64 = physics.segments.iter().map(Segment::axial_length).sum();
    let quarter = core::f64::consts::TAU * gas.sound_speed / (4.0 * length);
    let wall = wall_pin(duct.wall);
    let line = first_impedance_peak(0.25 * quarter, 6.0 * quarter, 920, |omega| {
        fs_duct::input_impedance_wall(
            physics,
            gas,
            omega,
            fs_duct::LossModel::Bessel,
            termination,
            wall.as_ref(),
        )
        .map(|response| response.impedance.abs())
        .map_err(AcousticRealizeError::Duct)
    })?
    .ok_or(AcousticRealizeError::InvalidDescription {
        what: "no passive line resonance inside the lumped-ladder scan",
    })?;
    let chain = ChainImpedance::new(&ode_bore_chain(
        sections,
        duct_outlet(duct.termination),
        1,
        taps,
        gas,
        duct.wall,
    )?)?;
    let ladder = first_impedance_peak(line / 1.25, line * 1.25, 160, |omega| {
        chain.magnitude(omega)
    })?;
    Ok((ladder, line))
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
    outlet: ChainOutlet,
    inlets: usize,
    taps: &[AcousticTap],
    gas: &GasState,
    wall: Option<LocallyReactingWall>,
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    let wall = wall_pin(wall);
    acoustic_chain_outlet_wall(
        sections,
        gas.density,
        gas.sound_speed,
        outlet,
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
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    bore_spec_inner(duct, SECTION_BUDGET as usize, false)
}

fn bore_spec_with_cells(
    duct: &ViscothermalDuct,
    total_cells: usize,
) -> Option<(Vec<AcousticSection>, Vec<AcousticTap>)> {
    bore_spec_inner(duct, total_cells, true)
}

fn bore_spec_inner(
    duct: &ViscothermalDuct,
    total_cells: usize,
    staggered_half_cell_align: bool,
) -> Option<(Vec<AcousticSection>, Vec<AcousticTap>)> {
    if validate_tone_holes(duct).is_err() || total_cells < 2 {
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
    let budget = total_cells as f64;
    // The half-cell alignment offsets the LC chain's (N + 1/2) dx effective
    // length behind an open outlet. fs-phs realizes a lone widening taper as
    // a psi line instead, which has no such offset (lossless -1.3% at 8 cells
    // against -6% for a cylinder). A rigid outlet cuts the last flux, and N
    // charges with N - 1 fluxes span exactly N dx. Neither is shortened.
    let psi_cone = matches!(
        duct.segments.as_slice(),
        [s] if s.is_taper() && s.outlet_radius_m > s.radius_m
    ) && duct.termination != WaveguideEnd::Closed;
    let exact_length = psi_cone || duct.termination == WaveguideEnd::Closed;
    let length_scale = if staggered_half_cell_align && !exact_length {
        budget / (budget + 0.5)
    } else {
        1.0
    };
    let mut sections = Vec::new();
    for s in &duct.segments {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = ((budget * s.length_m / total).round() as usize).max(2);
        sections.push(AcousticSection {
            length: s.length_m * length_scale,
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

struct PreparedLinearZoh {
    n: usize,
    m: usize,
    phi_x: Vec<f64>,
    gamma_x: Vec<f64>,
    phi_y: Vec<f64>,
    gamma_y: Vec<f64>,
}

impl PreparedLinearZoh {
    fn from_phs_scattering(
        sys: &fs_phs::PortHamiltonian,
        zc: f64,
        dt: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let n = sys.state_dim();
        let m = sys.port_dim();
        let (j, r, g) = sys.structure();
        let mut q = vec![0.0; n * n];
        let mut unit = vec![0.0; n];
        for col in 0..n {
            unit[col] = 1.0;
            let e = sys.effort(&unit);
            unit[col] = 0.0;
            for row in 0..n {
                q[row * n + col] = e[row];
            }
        }
        let inv_zc = if zc > 0.0 && zc.is_finite() {
            1.0 / zc
        } else {
            0.0
        };
        // Scattering state matrix A_scat = (J - R - (1/Zc) g_0 g_0^T) * Q
        // and output matrix C = G^T * Q.
        let mut a = vec![0.0; n * n];
        for row in 0..n {
            let g_row0 = if m > 0 { g[row * m] } else { 0.0 };
            for k in 0..n {
                let g_k0 = if m > 0 { g[k * m] } else { 0.0 };
                let r_scat = r[row * n + k] + inv_zc * g_row0 * g_k0;
                let jrk = j[row * n + k] - r_scat;
                if jrk == 0.0 {
                    continue;
                }
                for col in 0..n {
                    a[row * n + col] += jrk * q[k * n + col];
                }
            }
        }
        let mut c = vec![0.0; m * n];
        for p in 0..m {
            for k in 0..n {
                let gkp = g[k * m + p];
                if gkp == 0.0 {
                    continue;
                }
                for col in 0..n {
                    c[p * n + col] += gkp * q[k * n + col];
                }
            }
        }
        let mut d = vec![1.0; n];
        let mut inv_d = vec![1.0; n];
        for i in 0..n {
            let qii = q[i * n + i];
            if qii > 0.0 && qii.is_finite() {
                let s = qii.sqrt();
                d[i] = s;
                inv_d[i] = 1.0 / s;
            }
        }
        let mut a_bal = vec![0.0; n * n];
        for row in 0..n {
            let d_row = d[row];
            for col in 0..n {
                a_bal[row * n + col] = d_row * a[row * n + col] * inv_d[col] * dt;
            }
        }
        let (phi_bal, f1_bal, f2_bal) = dense_phi_f1_f2(&a_bal, n)?;
        let mut phi_x = vec![0.0; n * n];
        let mut f1 = vec![0.0; n * n];
        let mut f2 = vec![0.0; n * n];
        for row in 0..n {
            let inv_d_row = inv_d[row];
            for col in 0..n {
                let scale = inv_d_row * d[col];
                let idx = row * n + col;
                phi_x[idx] = phi_bal[idx] * scale;
                f1[idx] = f1_bal[idx] * scale;
                f2[idx] = f2_bal[idx] * scale;
            }
        }
        let mut gamma_x = vec![0.0; n * m];
        let mut f2_g_dt = vec![0.0; n * m];
        for row in 0..n {
            for k in 0..n {
                let f1_rk = f1[row * n + k] * dt;
                let f2_rk = f2[row * n + k] * dt;
                for p in 0..m {
                    let g_scale = if p == 0 && inv_zc > 0.0 {
                        2.0 * inv_zc
                    } else {
                        1.0
                    };
                    let g_kp = g[k * m + p] * g_scale;
                    if g_kp != 0.0 {
                        gamma_x[row * m + p] += f1_rk * g_kp;
                        f2_g_dt[row * m + p] += f2_rk * g_kp;
                    }
                }
            }
        }
        let mut phi_y = vec![0.0; m * n];
        let mut gamma_y = vec![0.0; m * m];
        for p in 0..m {
            for k in 0..n {
                let c_pk = c[p * n + k];
                if c_pk == 0.0 {
                    continue;
                }
                for col in 0..n {
                    phi_y[p * n + col] += c_pk * f1[k * n + col];
                }
                for col_p in 0..m {
                    gamma_y[p * m + col_p] += c_pk * f2_g_dt[k * m + col_p];
                }
            }
        }
        Ok(Self {
            n,
            m,
            phi_x,
            gamma_x,
            phi_y,
            gamma_y,
        })
    }

    fn step(&self, x: &[f64], u: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut x_next = vec![0.0; self.n];
        for row in 0..self.n {
            let base = row * self.n;
            let mut sum = 0.0;
            for col in 0..self.n {
                sum += self.phi_x[base + col] * x[col];
            }
            for p in 0..self.m {
                sum += self.gamma_x[row * self.m + p] * u[p];
            }
            x_next[row] = sum;
        }
        let mut y = vec![0.0; self.m];
        for p in 0..self.m {
            let base = p * self.n;
            let mut sum = 0.0;
            for col in 0..self.n {
                sum += self.phi_y[base + col] * x[col];
            }
            for col_p in 0..self.m {
                sum += self.gamma_y[p * self.m + col_p] * u[col_p];
            }
            y[p] = sum;
        }
        (x_next, y)
    }
}

#[allow(clippy::type_complexity)]
fn dense_phi_f1_f2(
    a: &[f64],
    dim: usize,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), AcousticRealizeError> {
    let mut norm1 = 0.0_f64;
    for col in 0..dim {
        let mut col_sum = 0.0;
        for row in 0..dim {
            col_sum += a[row * dim + col].abs();
        }
        if !col_sum.is_finite() {
            return Err(AcousticRealizeError::Nonlinear(
                "non-finite linear ZOH generator".to_string(),
            ));
        }
        norm1 = norm1.max(col_sum);
    }
    let mut s = 0u32;
    let mut scaled_norm = norm1;
    while scaled_norm > 0.5 && s < 30 {
        scaled_norm *= 0.5;
        s += 1;
    }
    let scale = 0.5_f64.powi(s as i32);
    let mut a_s = vec![0.0; dim * dim];
    for (dst, &src) in a_s.iter_mut().zip(a.iter()) {
        *dst = src * scale;
    }
    let mut phi = vec![0.0; dim * dim];
    let mut f1 = vec![0.0; dim * dim];
    let mut f2 = vec![0.0; dim * dim];
    let mut term = vec![0.0; dim * dim];
    for i in 0..dim {
        let diag = i * dim + i;
        phi[diag] = 1.0;
        f1[diag] = 1.0;
        f2[diag] = 0.5;
        term[diag] = 1.0;
    }
    let mut next_term = vec![0.0; dim * dim];
    for k in 1..=24 {
        let inv_k = 1.0 / f64::from(k);
        next_term.fill(0.0);
        for i in 0..dim {
            let row_i = i * dim;
            for p in 0..dim {
                let t_ip = term[row_i + p];
                if t_ip == 0.0 {
                    continue;
                }
                let scaled = t_ip * inv_k;
                let row_p = p * dim;
                for j in 0..dim {
                    next_term[row_i + j] += scaled * a_s[row_p + j];
                }
            }
        }
        term.copy_from_slice(&next_term);
        let c1 = 1.0 / f64::from(k + 1);
        let c2 = 1.0 / (f64::from(k + 1) * f64::from(k + 2));
        let mut max_delta = 0.0_f64;
        for idx in 0..(dim * dim) {
            let t = term[idx];
            phi[idx] += t;
            f1[idx] += t * c1;
            f2[idx] += t * c2;
            max_delta = max_delta.max(t.abs());
        }
        if max_delta <= 1.0e-16 {
            break;
        }
    }
    let mut one_plus_phi = vec![0.0; dim * dim];
    let mut prod = vec![0.0; dim * dim];
    for _ in 0..s {
        one_plus_phi.copy_from_slice(&phi);
        for i in 0..dim {
            one_plus_phi[i * dim + i] += 1.0;
        }
        // f2_next = 0.25 * ((I + phi) * f2 + f1)
        dense_matmul(&one_plus_phi, &f2, &mut prod, dim);
        for idx in 0..(dim * dim) {
            f2[idx] = 0.25 * (prod[idx] + f1[idx]);
        }
        // f1_next = 0.5 * (I + phi) * f1
        dense_matmul(&one_plus_phi, &f1, &mut prod, dim);
        for idx in 0..(dim * dim) {
            f1[idx] = 0.5 * prod[idx];
        }
        // phi_next = phi * phi
        dense_matmul(&phi, &phi, &mut prod, dim);
        phi.copy_from_slice(&prod);
    }
    Ok((phi, f1, f2))
}

fn dense_matmul(a: &[f64], b: &[f64], out: &mut [f64], dim: usize) {
    out.fill(0.0);
    for i in 0..dim {
        let row_i = i * dim;
        for p in 0..dim {
            let a_ip = a[row_i + p];
            if a_ip == 0.0 {
                continue;
            }
            let row_p = p * dim;
            for j in 0..dim {
                out[row_i + j] += a_ip * b[row_p + j];
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn realize_blown_ode(
    sections: &[AcousticSection],
    outlet: ChainOutlet,
    taps: &[AcousticTap],
    blow: VolumeVelocityPulse,
    plates: &mut PlateBank,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<LocallyReactingWall>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let inlets = if plates.is_empty() { 1 } else { 2 };
    let line = ode_bore_chain(sections, outlet, inlets, taps, gas, wall)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let (sys, mut x, observation) = compose_plate_chain(plates, line, 1, gas, sample_rate_hz)?;
    let linear = prepare_render_linear(&sys)?;
    let exact = prepare_chain_exact(&sys, &linear, dt)?;
    // A complete bank retains its independent external force port at slot zero.
    // Only the remaining duct inlet receives this prescribed volume source;
    // pressure reaction must never depend on an unrelated point-drive footprint.
    let source = usize::from(observation.is_some());
    let mut input = vec![0.0; sys.port_dim()];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 * dt;
        input[source] = if t < blow.duration_s {
            blow.peak_m3_s * det::sin(core::f64::consts::PI * t / blow.duration_s)
        } else {
            0.0
        };
        let rec = step_chain_flow(&sys, &linear, exact.as_ref(), &x, &input, dt)?;
        let p = rec.y[source]
            + plate_chain_radiation(observation.as_ref(), &rec.x, &x, gas, listener_m, dt)?;
        if !p.is_finite() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "blown ODE pressure left the finite set",
            });
        }
        out.push(p);
        x = rec.x;
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)] // one coherent realization
fn realize_reed_ode(
    sections: &[AcousticSection],
    outlet: ChainOutlet,
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
    if !crate::reed_bore::reed_parameters_valid(reed) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "reed parameters must be physical and finite",
        });
    }
    let radius = sections.first().map_or(0.0, |s| s.radius);
    let area = core::f64::consts::PI * radius * radius;
    let zc = gas.density * gas.sound_speed / area;
    let inlets = if plates.linear.is_empty() { 1 } else { 2 };
    let line = ode_bore_chain(sections, outlet, inlets, taps, gas, wall)?;
    let bore_sys = if plates.linear.is_empty() {
        line
    } else {
        let plate =
            plate_modal_ports(plates, false)?.ok_or(AcousticRealizeError::InvalidDescription {
                what: "ODE reed expected a linear plate bank",
            })?;
        transformer(plate, line, 0, 1, area)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
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
        let linear = prepare_render_linear(&phs)?;
        Some((phs, linear, face, lay))
    } else {
        None
    };
    let dt = 1.0 / f64::from(sample_rate_hz);
    let block_fits_phs =
        2 * (bore_sys.state_dim() + bore_sys.port_dim() + 1)
            <= fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;
    let bore_linear = if block_fits_phs && n <= 4 {
        Some(prepare_render_linear(&bore_sys)?)
    } else {
        None
    };
    let zoh = if bore_linear.is_none() {
        let zc_scat = if massive.is_none() { zc } else { 0.0 };
        Some(PreparedLinearZoh::from_phs_scattering(&bore_sys, zc_scat, dt)?)
    } else {
        None
    };
    let mut x = vec![0.0; bore_sys.state_dim()];
    let mut x_reed = vec![0.0; 2];
    let mut p_plus_prev = 5.0;
    let mut p_minus_prev = 0.0;
    if massive.is_none()
        && let Some(z) = &zoh
    {
        let (nx, ny) = z.step(&x, &[p_plus_prev]);
        p_minus_prev = ny[0] - p_plus_prev;
        x = nx;
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 * dt;
        let p_m = blowing_envelope(reed, t);
        let p_bore = bore_sys.output(&x).first().copied().unwrap_or(0.0);
        if let Some((reed_phs, reed_linear, face, lay)) = &massive {
            let q = x_reed[0];
            let opening = (reed.rest_opening_m + q).max(0.0);
            let v = reed_phs.output(&x_reed).first().copied().unwrap_or(0.0);
            let dp = p_m - p_bore;
            let u_jet = bernoulli_volume_flow(reed.width_m, opening, dp, gas.density);
            let u_swept = reed_swept_flow(*face, v);
            let y = reed.rest_opening_m + q;
            let mut f = -face * dp;
            let contact = slit_contact_force(lay, y)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            f += contact;
            f += lay
                .dissipative_modal_forces(1, &[y, v], &[v])
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?[0];
            let rec_r = fs_phs::linear_exact::step_render(reed_phs, reed_linear, &x_reed, &[f], dt)
                .map_err(render_phs_error)?;
            if !rec_r.x.iter().all(|xi| xi.is_finite()) {
                return Err(AcousticRealizeError::Reed {
                    what: "massive reed left the finite set",
                });
            }
            x_reed = rec_r.x;
            let u = u_jet + u_swept;
            let (next_x, y0) = if let Some(linear) = &bore_linear {
                let rec = fs_phs::linear_exact::step_render(&bore_sys, linear, &x, &[u], dt)
                    .map_err(render_phs_error)?;
                (rec.x, rec.y[0])
            } else {
                let (nx, ny) = zoh.as_ref().expect("precomputed ZOH").step(&x, &[u]);
                (nx, ny[0])
            };
            let mut p = y0;
            if !plates.linear.is_empty() {
                for (k, body) in plates.linear.iter().enumerate() {
                    let acc = (next_x[2 * k + 1] - x[2 * k + 1]) / dt;
                    p += body.radiate_mass_normalized(acc, gas.density, listener_m);
                }
            }
            if !p.is_finite() {
                return Err(AcousticRealizeError::Reed {
                    what: "ODE reed pressure left the finite set",
                });
            }
            x = next_x;
            out.push(p);
            continue;
        }
        let (next_x, y0) = if let Some(linear) = &bore_linear {
            let h = quasistatic_aperture_opening(
                reed.rest_opening_m,
                reed.closing_pressure_pa,
                p_m - p_bore,
            );
            let u = bernoulli_volume_flow(reed.width_m, h, p_m - p_bore, gas.density);
            let rec = fs_phs::linear_exact::step_render(&bore_sys, linear, &x, &[u], dt)
                .map_err(render_phs_error)?;
            (rec.x, rec.y[0])
        } else {
            let z = zoh.as_ref().expect("precomputed ZOH");
            let p_minus = p_minus_prev;
            let p_plus = crate::reed_bore::solve_reed_wave(
                reed,
                gas.density,
                zc,
                0.0,
                p_minus,
                p_m,
                p_plus_prev,
                0.0,
            )?;
            p_plus_prev = p_plus;
            let (nx, ny) = z.step(&x, &[p_plus]);
            let p_minus_now = ny[0] - p_plus;
            p_minus_prev = p_minus_now;
            (nx, p_plus + p_minus_now)
        };
        let mut p = y0;
        if !plates.linear.is_empty() {
            for (k, body) in plates.linear.iter().enumerate() {
                let acc = (next_x[2 * k + 1] - x[2 * k + 1]) / dt;
                p += body.radiate_mass_normalized(acc, gas.density, listener_m);
            }
        }
        if !p.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "ODE reed pressure left the finite set",
            });
        }
        x = next_x;
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
    if let Some((sections, taps)) = bore_spec(duct) {
        return realize_blown_ode(
            &sections,
            duct_outlet(duct.termination),
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
    // Isolated linear blow on non-ODE geometry (e.g. conical taper):
    // the same causal DelayedFilter port as a reed or a body, filled
    // with the conjugated IFFT[Z] under the shared e^{-iωt} → DFT
    // convention (DEF-01).
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
        .inlet_radius();
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

/// Chain outlet of a scenario termination. `Closed` is a rigid closure
/// (`U = 0`), never the chain's pressure-release end (frankensim-gs17m).
fn duct_outlet(end: WaveguideEnd) -> ChainOutlet {
    match end {
        WaveguideEnd::Closed => ChainOutlet::Rigid,
        WaveguideEnd::UnflangedOpen => ChainOutlet::Radiating(MouthFlange::Unflanged),
        WaveguideEnd::FlangedOpen => ChainOutlet::Radiating(MouthFlange::Flanged),
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
            let image = plates.physical_image(&gas, 48_000, &fs_exec::CancelGate::default())
                .unwrap().unwrap();
            let (sys, zero, observation) = image.into_parts();
            let rec = step(&sys, &zero, &[force, face_pressure], dt).unwrap();
            let p = plate_chain_radiation(Some(&observation), &rec.x, &zero, &gas, 1.0, dt).unwrap();
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
            tension_loss_factor: None,
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

    #[test]
    fn kc_render_refuses_unresolved_linear_poles_without_advancing() {
        let string = PrestressedString {
            length_m: 0.5,
            tension_n: 20.0,
            lin_density_kg_m: 0.006,
            axial_stiffness_n: 1000.0,
            width_m: 0.001,
            n_modes: 2,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor: None,
        };
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let mut member = kc_string_member(&string, None, None, &[], &gas, 1.0, 1.0).unwrap();
        let mut texture = TextureDrive::try_new(None).unwrap();
        let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
        let dt = cap / member.linear.pole_rate_bound_rad_s();
        assert!(
            step_kc_member(&mut member, None, &mut texture, &[], &string, &gas, dt * (1.0 - 1e-12))
                .is_ok()
        );
        let before = member.x.clone();
        assert!(matches!(
            step_kc_member(&mut member, None, &mut texture, &[], &string, &gas, dt * (1.0 + 1e-10)),
            Err(AcousticRealizeError::TimeResolution(_))
        ));
        assert_eq!(member.x, before);
    }

    #[test]
    fn coupled_ode_linear_step_uses_the_complete_transformer() {
        let plate = fs_phs::modal_bank_ports(&[1.0, 1.7], &[0.0, 0.0], &[&[1.0, 0.5]]).unwrap();
        let line = fs_phs::mass_spring_damper(1.0, 2.0, 0.0).unwrap();
        let sys = transformer(plate, line, 0, 0, 0.7).unwrap();
        let linear = prepare_render_linear(&sys).unwrap();
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let x = [0.1, 0.02, -0.03, 0.01, 0.04, -0.02];
        let full = step_plate_chain(&sys, &linear, None, &x, 0.0, None, None, &gas, 0.0, 0.5).unwrap();
        let half = step_plate_chain(&sys, &linear, None, &x, 0.0, None, None, &gas, 0.0, 0.25).unwrap();
        let halves =
            step_plate_chain(&sys, &linear, None, &half.x, 0.0, None, None, &gas, 0.25, 0.25).unwrap();
        for (a, b) in full.x.iter().zip(halves.x) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(full.balance_residual().abs() < 1e-12);
    }

    #[test]
    fn blown_ode_observes_the_exact_held_flow_step() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let sections = [AcousticSection {
            length: 0.5,
            radius: 0.01,
            outlet_radius: 0.0,
            cells: 2,
        }];
        for rate in [48_000, 96_000, 192_000] {
            let dt = 1.0 / f64::from(rate);
            let out = realize_blown_ode(
                &sections,
                ChainOutlet::PressureRelease,
                &[],
                VolumeVelocityPulse {
                    peak_m3_s: 1e-4,
                    duration_s: 2.0 * dt,
                },
                &mut PlateBank::default(),
                &gas,
                1.0,
                rate,
                2,
                None,
            )
            .unwrap();
            let line = ode_bore_chain(&sections, ChainOutlet::PressureRelease, 1, &[], &gas, None).unwrap();
            let split = prepare_render_linear(&line).unwrap();
            let expected = fs_phs::linear_exact::step_split(
                &line,
                &split,
                fs_phs::linear_exact::LinearPropagator::ExactZoh,
                &vec![0.0; line.state_dim()],
                &[1e-4],
                dt,
            )
            .unwrap();
            assert_eq!(out[0], 0.0);
            assert!((out[1] - expected.y[0]).abs() <= 1e-12 * expected.y[0].abs());
            let legacy = fs_phs::step(&line, &vec![0.0; line.state_dim()], &[1e-4], dt).unwrap();
            assert!((out[1] - legacy.y[0]).abs() > 1e-9 * expected.y[0].abs());
        }
    }

    #[test]
    fn reed_ode_observes_exact_held_flow_for_both_valve_selections() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let sections = [AcousticSection {
            length: 0.5,
            radius: 0.01,
            outlet_radius: 0.0,
            cells: 2,
        }];
        for rate in [48_000, 96_000, 192_000] {
            for mass in [0.0, 1e-4] {
                let dt = 1.0 / f64::from(rate);
                let reed = BeatingReed {
                    rest_opening_m: 1e-3,
                    width_m: 0.01,
                    closing_pressure_pa: 1e6,
                    blowing_pressure_pa: 100.0,
                    attack_s: 2.0 * dt,
                    mass_kg: mass,
                    stiffness_n_m: 100.0,
                    damping_ratio: 0.01,
                };
                let out = realize_reed_ode(
                    &sections,
                    ChainOutlet::PressureRelease,
                    &[],
                    reed,
                    &mut PlateBank::default(),
                    &gas,
                    1.0,
                    rate,
                    2,
                    None,
                )
                .unwrap();
                let pressure = blowing_envelope(reed, dt);
                let opening = if mass == 0.0 {
                    quasistatic_aperture_opening(
                        reed.rest_opening_m,
                        reed.closing_pressure_pa,
                        pressure,
                    )
                } else {
                    reed.rest_opening_m
                };
                let flow = bernoulli_volume_flow(reed.width_m, opening, pressure, gas.density);
                let line = ode_bore_chain(&sections, ChainOutlet::PressureRelease, 1, &[], &gas, None).unwrap();
                let split = prepare_render_linear(&line).unwrap();
                let expected = fs_phs::linear_exact::step_split(
                    &line,
                    &split,
                    fs_phs::linear_exact::LinearPropagator::ExactZoh,
                    &vec![0.0; line.state_dim()],
                    &[flow],
                    dt,
                )
                .unwrap();
                assert_eq!(out[0], 0.0);
                assert!((out[1] - expected.y[0]).abs() <= 1e-12 * expected.y[0].abs());
            }
        }
    }

    #[test]
    fn actual_dirac_renderer_refuses_an_unresolved_common_clock() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        for rate in [48_000, 96_000, 192_000] {
            let mut plates = PlateBank::default();
            plates.linear.push(
                crate::thin_plate::CompactBody::from_radiator(RadiatingPlate {
                    area_m2: 0.1,
                    mass_kg: 1.0,
                    frequency_hz: 1.0,
                    damping_ratio: 0.0,
                })
                .unwrap(),
            );
            let make = |frequency: f64| PrestressedString {
                length_m: 0.5,
                tension_n: (2.0 * frequency).powi(2),
                lin_density_kg_m: 1.0,
                axial_stiffness_n: 0.0,
                width_m: 0.001,
                n_modes: 1,
                damping_ratio: 0.0,
                rayleigh: None,
                bending_stiffness_n_m2: 0.0,
                kelvin_voigt_bending: None,
                relaxation_bending: None,
                polarization_detune: 0.0,
                moving_end: true,
                tension_loss_factor: None,
            };
            let render = |string: &PrestressedString| {
                realize_dirac_join(
                    string,
                    Some(Pluck {
                        station_frac: 0.5,
                        height_m: 1e-8,
                    }),
                    None,
                    &[],
                    None,
                    &plates,
                    None,
                    None,
                    None,
                    None,
                    None,
                    &gas,
                    1.0,
                    rate,
                    1,
                )
            };
            assert!(render(&make(1.0)).is_ok());
            assert!(matches!(
                render(&make(10_000.0)),
                Err(AcousticRealizeError::TimeResolution(_))
            ));
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

#[cfg(test)]
mod lumped_ladder_resonance_tests {
    use super::*;
    use fs_scenario::{CylinderSegment, ToneHole};

    fn resonances(duct: &ViscothermalDuct, cells: usize) -> (Option<f64>, f64) {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let termination = match duct.termination {
            WaveguideEnd::Closed => Termination::Closed,
            WaveguideEnd::UnflangedOpen => Termination::UnflangedOpen,
            WaveguideEnd::FlangedOpen => Termination::FlangedOpen,
        };
        let (sections, taps) = bore_spec_with_cells(duct, cells).unwrap();
        let physics = physics_duct(duct).unwrap();
        lumped_ladder_resonances(&physics, &sections, &taps, duct, termination, &gas).unwrap()
    }

    #[test]
    fn a_window_edge_is_never_a_resonance() {
        assert_eq!(
            first_impedance_peak(1.0, 2.0, 920, |w| Ok(w)).unwrap(),
            None
        );
        assert_eq!(
            first_impedance_peak(1.0, 2.0, 920, |w| Ok(-w)).unwrap(),
            None
        );
        // A compliance-like edge far stronger than the resonance at 1.7.
        let peak = first_impedance_peak(1.0, 2.0, 920, |w| {
            Ok(1.0e5 * det::exp(-50.0 * (w - 1.0)) + 1.0 / ((w - 1.7) * (w - 1.7) + 1.0e-4))
        })
        .unwrap()
        .unwrap();
        assert!((peak - 1.7).abs() < 1.0e-9, "{peak}");
    }

    #[test]
    fn closed_end_cone_and_open_hole_compare_their_lowest_resonance() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let length = 0.5;
        let quarter = core::f64::consts::TAU * gas.sound_speed / (4.0 * length);
        let duct = |segments: Vec<CylinderSegment>, holes: Vec<ToneHole>, end| ViscothermalDuct {
            segments,
            tone_holes: holes,
            termination: end,
            wall: None,
        };
        let closed = duct(
            vec![CylinderSegment::cylinder(0.0075, length)],
            Vec::new(),
            WaveguideEnd::Closed,
        );
        let cone = duct(
            vec![CylinderSegment::taper(0.003, 0.02, length)],
            Vec::new(),
            WaveguideEnd::UnflangedOpen,
        );
        let holed = duct(
            vec![
                CylinderSegment::cylinder(0.0075, 0.5 * length),
                CylinderSegment::cylinder(0.0075, 0.5 * length),
            ],
            vec![ToneHole {
                after_segment: 0,
                radius_m: 0.004,
                chimney_m: 0.003,
                open_fraction: 1.0,
            }],
            WaveguideEnd::UnflangedOpen,
        );
        // Lossless truncated cone driven at its narrow end, p = 0 one unflanged
        // end correction past its wide end: k L' + atan(k x1) = pi, x1 the
        // apex distance of the inlet.
        let (x1, wide) = (length * 0.003 / (0.02 - 0.003), length + 0.6133 * 0.02);
        let (mut lo, mut hi) = (1.0, 20.0);
        for _ in 0..200 {
            let k = f64::midpoint(lo, hi);
            if k * wide + det::atan(k * x1) < core::f64::consts::PI {
                lo = k;
            } else {
                hi = k;
            }
        }
        let cone_oracle = f64::midpoint(lo, hi) * gas.sound_speed;
        // The hole lands in a cell, so coarse ladders misplace it (b151: -114
        // cents at 6 cells, +185 at 12, -25 at 24, -12 at 48).
        let rungs: [(&str, ViscothermalDuct, &[usize]); 3] = [
            ("closed", closed, &[6, 12]),
            ("cone", cone, &[6, 12]),
            ("hole", holed, &[24]),
        ];
        for (name, duct, cell_counts) in rungs {
            for &cells in cell_counts {
                let (ladder, line) = resonances(&duct, cells);
                let cents = ladder.map(|l| 1200.0 * det::ln(l / line) / core::f64::consts::LN_2);
                println!(
                    "{{\"case\":\"{name}\",\"cells\":{cells},\"line_over_quarter\":{},\"ladder_over_quarter\":{:?},\"cents\":{cents:?}}}",
                    line / quarter,
                    ladder.map(|l| l / quarter)
                );
                // Every case's lowest resonance sits outside the old
                // [0.5, 1.5] quarter-wave window.
                assert!(line / quarter > 1.5, "{name}: {}", line / quarter);
                match name {
                    "closed" => {
                        // The half-wave, lowered only by viscothermal dispersion.
                        assert!((line / (2.0 * quarter) - 1.0).abs() < 2.0e-2);
                        // The rigid chain outlet (frankensim-gs17m) is admitted
                        // and measured; as pressure release it sat an octave low
                        // (0.987 x c/4L) and was refused.
                        assert!(cents.is_some_and(|c| c.abs() < 100.0), "{name} {cells}");
                    }
                    "cone" => {
                        assert!(
                            (line / cone_oracle - 1.0).abs() < 3.0e-2,
                            "{}",
                            line / cone_oracle
                        );
                        // The psi cone with its one inlet near-field shunt is
                        // admitted and measured (frankensim-d1j6s); with a
                        // shunt on every cell it ran away to 2.6 and 3.3 x c/4L.
                        assert!(cents.is_some_and(|c| c.abs() < 100.0), "{name} {cells}");
                    }
                    _ => {
                        // The open hole at mid-bore shortens the speaking length.
                        assert!(line / quarter < 2.0, "{}", line / quarter);
                        assert!(cents.is_some_and(|c| c.abs() < 100.0), "{name} {cells}");
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod arbitrary_gas_realization_tests {
    use super::*;
    use fs_material::gas::{ConductivityModel, GasSpec};
    use fs_scenario::{Listener, Pluck};

    #[test]
    fn observer_absorption_law_is_explicit_for_thermochemical_gas() {
        use crate::air_path::{GasAbsorptionModel, absorb_pressure_history_with_model};
        let description = fixture();
        let mut gas = gas_state(description.ambient).unwrap();
        // Perturb a coherent coefficient set to distinguish an externally
        // resolved state from the legacy atmosphere. This unit tests routing;
        // source-backed NASA-9 numerical values are checked in matdb_pack_cli.
        gas.specific_heat_cp *= 1.001;
        gas.gamma = gas.specific_heat_cp / (gas.specific_heat_cp - gas.specific_gas_constant);
        gas.sound_speed = det::sqrt(gas.gamma * gas.specific_gas_constant * gas.temperature);
        gas.characteristic_impedance = gas.density * gas.sound_speed;
        gas.prandtl = gas.dynamic_viscosity * gas.specific_heat_cp / gas.thermal_conductivity;
        gas.validate_coherence().unwrap();
        let history: Vec<_> = (0..128).map(|i| det::sin(i as f64 * 0.29)).collect();
        let model = GasAbsorptionModel::Iso9613Air {
            relative_humidity: 0.5,
            air: admitted_air(gas),
        };
        // The legacy atmosphere's admission does not transfer to the external
        // state, even though both are air: ISO needs this state's own.
        let (_, legacy_air) = legacy_admission(&description);
        assert_eq!(
            finish_observer(&description, &gas, history.clone(),
                Some(GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5, air: legacy_air }).into()),
            Err(AcousticRealizeError::NonTerrestrialAir(
                fs_material::iso9613::Iso9613Refusal::DifferentState
            ))
        );
        let mut expected = history.clone();
        absorb_pressure_history_with_model(
            &mut expected,
            1.0 / 16000.0,
            description.listener.distance_m,
            &gas,
            model,
        )
        .unwrap();
        let result = finish_observer(&description, &gas, history.clone(), Some(model).into()).unwrap();
        assert_eq!(result.pressure_pa, expected);
        let classical = finish_observer(
            &description,
            &gas,
            history.clone(),
            Some(GasAbsorptionModel::StokesKirchhoff).into(),
        )
        .unwrap();
        assert_ne!(result.pressure_pa, classical.pressure_pa);
        assert!(
            finish_observer(
                &description,
                &gas,
                history.clone(),
                Some(GasAbsorptionModel::Iso9613Air {
                    relative_humidity: f64::NAN,
                    air: admitted_air(gas),
                })
                .into()
            )
            .is_err()
        );
        let (legacy, legacy_air) = legacy_admission(&description);
        assert_eq!(legacy, gas_state(description.ambient).unwrap());
        let scenario_law = GasAbsorptionModel::Iso9613Air {
            relative_humidity: description.ambient.relative_humidity,
            air: legacy_air,
        };
        let mut original = history.clone();
        absorb_pressure_history_with_model(
            &mut original,
            1.0 / 16000.0,
            description.listener.distance_m,
            &legacy,
            scenario_law,
        )
        .unwrap();
        assert_eq!(
            finish_observer(&description, &legacy, history, ObserverLaw::default())
                .unwrap()
                .pressure_pa,
            original
        );
        assert_eq!(
            realize_assembly_with_gas_and_absorption(&description, &legacy, scenario_law)
                .unwrap()
                .pressure_pa,
            realize_assembly(&description).unwrap().pressure_pa
        );
    }

    /// The scenario's legacy moist air with its own terrestrial-air admission.
    fn legacy_admission(description: &AcousticAssembly) -> (GasState, fs_material::iso9613::TerrestrialAir) {
        let ambient = description.ambient;
        fs_material::iso9613::TerrestrialAir::moist_air(
            ambient.temperature_k,
            ambient.pressure_pa,
            ambient.relative_humidity,
        )
        .unwrap()
    }

    /// `gas` with an explicit U.S. Standard Atmosphere 1976 dry composition
    /// plus the state's own water.
    fn ussa_medium(gas: GasState) -> fs_material::relaxation::RelaxingMedium {
        let dry_parts = [("N2", 0.780_84), ("O2", 0.209_476), ("Ar", 0.009_34), ("CO2", 0.000_314)];
        let total: f64 = dry_parts.iter().map(|(_, x)| x).sum();
        let water = gas.water_mole_fraction;
        let mut components: Vec<_> = dry_parts
            .iter()
            .map(|(name, x)| (fs_thermochem::SpeciesId::new(*name).unwrap(), x / total * (1.0 - water)))
            .collect();
        if water > 0.0 {
            components.push((fs_thermochem::SpeciesId::new("H2O").unwrap(), water));
        }
        fs_material::relaxation::RelaxingMedium::new(gas, components).unwrap()
    }

    /// Terrestrial-air admission of `gas` from its explicit USSA composition.
    fn admitted_air(gas: GasState) -> fs_material::iso9613::TerrestrialAir {
        fs_material::iso9613::TerrestrialAir::admit(&ussa_medium(gas)).unwrap()
    }

    /// MolecularRelaxation applies the evaluated law as the causal transfer and
    /// names itself in the receipt, with nothing omitted. It refuses:
    /// - a missing state;
    /// - a state evaluated for another gas;
    /// - a state supplied with another law, or with no law at all.
    #[test]
    fn molecular_relaxation_observer_is_the_causal_law_and_refuses_mismatches() {
        use crate::air_path::{GasAbsorptionModel, propagate_pressure_history_relaxing};
        let description = fixture();
        let (legacy, _) = legacy_admission(&description);
        let state = fs_material::relaxation::RelaxationAbsorption::iso9613_terrestrial_air()
            .unwrap()
            .evaluate(&ussa_medium(legacy))
            .unwrap();
        let relaxing = GasAbsorptionModel::MolecularRelaxation;
        let history: Vec<_> = (0..128).map(|i| det::sin(i as f64 * 0.29)).collect();
        let law = ObserverLaw { model: Some(relaxing), relaxation: Some(&state) };
        let realized = finish_observer(&description, &legacy, history.clone(), law).unwrap();
        let mut expected = history.clone();
        propagate_pressure_history_relaxing(&mut expected, 1.0 / 16000.0,
            description.listener.distance_m, &state).unwrap();
        assert_eq!(realized.pressure_pa, expected);
        assert_eq!(realized.observer_absorption.law, relaxing);
        assert_eq!(realized.observer_absorption.relaxation_omitted, None);
        assert!(finish_observer(&description, &helium_like(), history, law).is_err());
        let options = RealizeOptions {
            gas: Some(legacy),
            observer_absorption: Some(relaxing),
            observer_relaxation: Some(&state),
            ..RealizeOptions::default()
        };
        let rendered = realize_assembly_with_options(&description, &options).unwrap();
        assert_eq!(rendered.observer_absorption.law, relaxing);
        assert_ne!(rendered.pressure_pa, realize_assembly(&description).unwrap().pressure_pa);
        for refused in [
            RealizeOptions { observer_relaxation: None, ..options },
            RealizeOptions { gas: Some(helium_like()), ..options },
            RealizeOptions { observer_absorption: Some(GasAbsorptionModel::StokesKirchhoff), ..options },
            RealizeOptions { gas: None, observer_absorption: None, ..options },
        ] {
            assert!(matches!(
                realize_assembly_with_options(&description, &refused),
                Err(AcousticRealizeError::InvalidDescription { .. })
            ));
        }
    }

    /// Every realized observer names its absorption law. Scenario air carries
    /// ISO 9613 with nothing omitted. A supplied gas on Stokes-Kirchhoff, by
    /// default or by choice, carries the molecular-relaxation omission. Its
    /// figure at the listener distance over the band up to Nyquist is labelled
    /// the standard-air convention, not a bound. ISO asked of the helium-like
    /// gas refuses. Numerics are unchanged.
    #[test]
    fn observer_absorption_receipt_names_the_law_and_labels_any_omitted_relaxation() {
        use crate::air_path::{GasAbsorptionModel, OmissionBound, relaxation_omission};
        let description = fixture();
        let (legacy, legacy_air) = legacy_admission(&description);
        let history: Vec<_> = (0..128).map(|i| det::sin(i as f64 * 0.29)).collect();
        let scenario = finish_observer(&description, &legacy, history.clone(), ObserverLaw::default()).unwrap();
        assert_eq!(
            scenario.observer_absorption.law,
            GasAbsorptionModel::Iso9613Air {
                relative_humidity: description.ambient.relative_humidity,
                air: legacy_air,
            }
        );
        assert_eq!(scenario.observer_absorption.relaxation_omitted, None);
        let supplied = helium_like();
        let classical = finish_observer(&description, &supplied, history.clone(),
            Some(GasAbsorptionModel::StokesKirchhoff).into()).unwrap();
        let expected = relaxation_omission(description.listener.distance_m,
            0.5 * f64::from(description.sample_rate_hz)).unwrap();
        assert_eq!(classical.observer_absorption.law, GasAbsorptionModel::StokesKirchhoff);
        assert_eq!(classical.observer_absorption.relaxation_omitted, Some(expected));
        let OmissionBound::StandardAirConvention { db } = expected.bound;
        assert!(db > 0.0);
        let realized = realize_assembly_with_gas(&description, &supplied).unwrap();
        assert_eq!(realized.observer_absorption.relaxation_omitted, Some(expected));
        let mut reference = history.clone();
        crate::air_path::absorb_pressure_history_with_model(&mut reference, 1.0 / 16000.0,
            description.listener.distance_m, &supplied, GasAbsorptionModel::StokesKirchhoff).unwrap();
        assert_eq!(classical.pressure_pa, reference);
        assert!(matches!(
            realize_assembly_with_gas_and_absorption(&description, &supplied,
                GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5, air: legacy_air }),
            Err(AcousticRealizeError::NonTerrestrialAir(_))
        ));
        println!(
            "{{\"bead\":\"M0.1-W2\",\"test\":\"observer_absorption_receipt\",\"range_m\":{},\"band_hz\":{},\"omitted_relaxation_standard_air_convention_db\":{db}}}",
            expected.range_m, expected.band_hz
        );
    }

    fn fixture() -> AcousticAssembly {
        AcousticAssembly {
            ambient: AmbientGas::sea_level(),
            string: Some(PrestressedString {
                length_m: 0.65, tension_n: 80.0, lin_density_kg_m: 0.006,
                axial_stiffness_n: 0.0, width_m: 0.001, n_modes: 6,
                damping_ratio: 0.0, rayleigh: None, bending_stiffness_n_m2: 0.0,
                kelvin_voigt_bending: None, relaxation_bending: None,
                polarization_detune: 0.0, moving_end: false, tension_loss_factor: None,
            }),
            duct: None,
            pluck: Some(Pluck { station_frac: 0.25, height_m: 1.0e-4 }),
            bow: None, blow: None, reed: None, soundboard: None, body_modes: vec![],
            plate: None, cavity: None, obstacles: vec![], contact_texture: None,
            listener: Listener { distance_m: 10.0 }, sample_rate_hz: 16_000, duration_s: 0.02,
        }
    }

    fn helium_like() -> GasState {
        GasState::try_new(&GasSpec {
            molar_mass: 4.002_602e-3, gamma: 5.0 / 3.0,
            sutherland_beta: 1.05e-6, sutherland_s: 79.4,
            conductivity: ConductivityModel::Eucken,
        }, 288.15, 101_325.0).expect("admitted helium-like ideal-gas fixture")
    }

    #[test]
    fn gas_data_only_swap_reaches_actual_pressure_realization() {
        let description = fixture();
        let air = realize_assembly(&description).expect("air realization");
        let helium = helium_like();
        let changed = realize_assembly_with_gas(&description, &helium).expect("gas realization");
        assert_eq!(changed.gas, helium);
        assert!(changed.gas.sound_speed > air.gas.sound_speed);
        assert_eq!(changed.pressure_pa.len(), air.pressure_pa.len());
        assert!(changed.pressure_pa.iter().all(|value| value.is_finite()));
        assert_ne!(changed.pressure_pa, air.pressure_pa);
        assert_eq!(description.ambient, AmbientGas::sea_level());
    }

    #[test]
    fn explicit_gas_law_does_not_depend_on_equality_with_an_ambient_preset() {
        use crate::air_path::GasAbsorptionModel;
        let description = fixture();
        let gas = gas_state(description.ambient).unwrap();
        let options = RealizeOptions {
            gas: Some(gas),
            ..RealizeOptions::default()
        };
        let supplied = realize_assembly_with_options(&description, &options).unwrap();
        let classical = realize_assembly_with_gas_and_absorption(
            &description,
            &gas,
            GasAbsorptionModel::StokesKirchhoff,
        )
        .unwrap();
        assert_eq!(supplied, classical);
        let mut unrelated_ambient = description.clone();
        unrelated_ambient.ambient.relative_humidity = 0.75;
        assert_eq!(
            supplied,
            realize_assembly_with_options(&unrelated_ambient, &options).unwrap()
        );
        let explicit_air = realize_assembly_with_options(
            &description,
            &RealizeOptions {
                observer_absorption: Some(GasAbsorptionModel::Iso9613Air {
                    relative_humidity: 0.0,
                    air: legacy_admission(&description).1,
                }),
                ..options
            },
        )
        .unwrap();
        assert_ne!(supplied.pressure_pa, explicit_air.pressure_pa);
        assert!(matches!(
            realize_assembly_with_options(
                &description,
                &RealizeOptions {
                    bore: BoreRealization::lumped_ladder(1),
                    ..options
                }
            ),
            Err(AcousticRealizeError::InvalidDescription {
                what: "explicit bore discretization requires a reed-driven duct without a string member"
            })
        ));
    }

    #[test]
    fn independently_mutated_gas_derived_fields_refuse() {
        let mut gas = helium_like();
        gas.characteristic_impedance *= 2.0;
        assert!(matches!(realize_assembly_with_gas(&fixture(), &gas), Err(AcousticRealizeError::Ambient(_))));
        let mut gas = helium_like();
        gas.prandtl = f64::NAN;
        assert!(matches!(realize_assembly_with_gas(&fixture(), &gas), Err(AcousticRealizeError::Ambient(_))));
    }

    #[test]
    fn acoustic_admission_cannot_disagree_with_the_material_eos() {
        let changes: [fn(&mut GasState); 5] = [
            |gas| gas.pressure *= 2.0,
            |gas| gas.temperature *= 1.1,
            |gas| gas.gamma = 1.3,
            |gas| gas.specific_gas_constant *= 1.1,
            |gas| {
                gas.sound_speed *= 1.1;
                gas.characteristic_impedance = gas.density * gas.sound_speed;
            },
        ];
        for change in changes {
            let mut gas = helium_like();
            change(&mut gas);
            assert!(matches!(realize_assembly_with_gas(&fixture(), &gas),
                Err(AcousticRealizeError::Ambient(_))), "accepted split physical state: {gas:?}");
        }
    }
}

#[cfg(test)]
mod tension_loss_tests {
    use super::*;
    use fs_scenario::acoustic::{BendingRelaxationBranch, KelvinVoigtBending, PronyBending};

    fn wire(tension_loss_factor: Option<f64>) -> PrestressedString {
        PrestressedString {
            length_m: 0.328,
            tension_n: 70.0,
            lin_density_kg_m: 4.0e-4,
            axial_stiffness_n: 0.0,
            width_m: 2.5e-4,
            n_modes: 24,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 4.0e-4,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor,
        }
    }

    fn mode(string: &PrestressedString, k: usize) -> (f64, f64) {
        (
            string_mode_omega(string, k),
            k as f64 * core::f64::consts::PI / string.length_m,
        )
    }

    #[test]
    fn tension_loss_is_eta_f_on_the_tension_share_of_every_mode_under_each_bending_law() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let eta = 6.0e-4;
        let kelvin_voigt = KelvinVoigtBending {
            viscous_stiffness_n_m2_s: 2.0e-8,
            omega_band_rad_s: (0.0, 1.0e6),
            material_state_identity: None,
        };
        let prony = PronyBending {
            branches: vec![BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 1.0e-4,
                relaxation_time_s: 1.0e-3,
            }],
            omega_band_rad_s: (0.0, 1.0e6),
            material_state_identity: None,
            source_properties: None,
        };
        let laws: [(Option<KelvinVoigtBending>, Option<PronyBending>); 3] = [
            (None, None),
            (Some(kelvin_voigt), None),
            (None, Some(prony)),
        ];
        for (kv, relaxation) in laws {
            let law = |factor| PrestressedString {
                kelvin_voigt_bending: kv,
                relaxation_bending: relaxation.clone(),
                ..wire(factor)
            };
            let (plain, zero, lossy) = (law(None), law(Some(0.0)), law(Some(eta)));
            let mut last = f64::INFINITY;
            for k in 1..=plain.n_modes {
                let (omega, wave_number) = mode(&plain, k);
                let base = mode_loss_parts(&plain, omega, wave_number, &gas).unwrap();
                // A zero factor is the absent one, bit for bit.
                let same = mode_loss_parts(&zero, omega, wave_number, &gas).unwrap();
                assert_eq!(same.total().to_bits(), base.total().to_bits());
                let parts = mode_loss_parts(&lossy, omega, wave_number, &gas).unwrap();
                assert_eq!(parts.air.to_bits(), base.air.to_bits());
                assert_eq!(parts.authored.to_bits(), base.authored.to_bits());
                let tension = plain.tension_n * wave_number * wave_number;
                let share =
                    tension / (tension + plain.bending_stiffness_n_m2 * wave_number.powi(4));
                // Q_F^-1 = 2 zeta_F, the tension's share of the mode's stiffness times eta_F.
                let inverse_q = 2.0 * (parts.material - base.material);
                assert!(
                    (inverse_q / (eta * share) - 1.0).abs() < 1.0e-9,
                    "mode {k}: 1/Q_F {inverse_q:e} vs eta_F x share {:e}",
                    eta * share
                );
                // Bending takes a growing share up the series, so the tension's loss dilutes.
                assert!(
                    inverse_q < last,
                    "mode {k}: {inverse_q:e} did not dilute below {last:e}"
                );
                last = inverse_q;
            }
        }
    }

    #[test]
    fn tension_loss_refuses_nonphysical_factors_and_authored_reduced_losses() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let (omega, wave_number) = mode(&wire(None), 3);
        for bad in [f64::NAN, f64::INFINITY, -1.0e-6] {
            assert!(
                mode_zeta(&wire(Some(bad)), omega, wave_number, &gas).is_err(),
                "accepted eta_F {bad}"
            );
        }
        let authored = PrestressedString {
            damping_ratio: 1.0e-3,
            ..wire(Some(1.0e-4))
        };
        let rayleigh = PrestressedString {
            rayleigh: Some(RayleighParams {
                alpha_per_s: 1.0,
                beta_s: 0.0,
            }),
            ..wire(Some(1.0e-4))
        };
        for reduced in [authored, rayleigh] {
            assert!(mode_zeta(&reduced, omega, wave_number, &gas).is_err());
            let without = PrestressedString {
                tension_loss_factor: None,
                ..reduced
            };
            assert!(mode_zeta(&without, omega, wave_number, &gas).is_ok());
        }
    }
}
