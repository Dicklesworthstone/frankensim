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
/// Regional reed-bore media and a physical exterior microphone on the same assembly.
pub mod duct_observer;
/// Native bell and tone-hole harmonics through the original Rayleigh receiver.
pub mod duct_harmonic;
mod thermal_step;
mod dirac_body;
mod reed_ode;
mod kc_observer;

#[cfg(test)]
mod polarization_tests;
#[cfg(test)]
mod moving_contact_tests;
#[cfg(test)]
mod kc_step_tests;
#[cfg(test)]
#[path = "acoustic_realize/string_observer_tests.rs"]
mod moving_string_observer_tests;

#[cfg(test)]
use reed_ode::realize_reed_ode;

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
    certified_material_chart_radiators, certified_radiators,
};
use crate::unilateral_contact::{
    modal_contact_forces, modal_friction_forces, modal_hunt_crossley_forces, slit_contact_force,
    slit_lay, station::StationContact,
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
    acoustic_chain_outlet_wall, bernoulli_volume_flow, mass_spring_damper,
    modal_bank_ports, quasistatic_aperture_opening, step, step_descriptor, transformer,
};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, BowStroke, ContactTexture,
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
    /// Original terminal-transfer or retarded-pressure observer refusal.
    DuctObserver(crate::driving_point::regional::outlet::OutletError),
    /// Kirchhoff–Carrier / pHS refusal.
    Nonlinear(String),
    /// Algebraic work/flow refusal with the original quantitative certificate.
    AlgebraicSolve(fs_phs::PhsError),
    /// Lumped LC ladder requested outside its admissible pitch-error band.
    LumpedLadderBand {
        /// Requested total cell count.
        cells: usize,
        /// Smallest cell count whose spatial screening estimate meets the
        /// tolerance. This does not admit material-interface eigenmode error.
        min_admissible_cells: usize,
        /// Refusing pitch error at `cells` [cents]: the spatial estimate below
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
    /// A fixed supply passage is so heavy against the moving reed that the three-level coupling (the passage steps
    /// on the reed's step-mean flow) leaves a parasitic alternating mode that decays too slowly. The mode's root is
    /// -mu/(mu+2) per step; it must fall to the declared decay within one period of the fastest admitted mode.
    SupplyParasiticRoot {
        /// Added-mass ratio mu = (rho L / A_t) sum(A_i^2 / m_i) of the passage against the moving masses.
        mass_ratio: f64,
        /// Magnitude of the parasitic root, mu / (mu + 2).
        root: f64,
        /// Largest admitted magnitude, decay^(max angular step / tau).
        maximum_root: f64,
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
            Self::DuctObserver(e) => write!(f, "{e}"),
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
            Self::SupplyParasiticRoot { mass_ratio, root, maximum_root } => write!(
                f,
                "FS-COUPLE-SUPPLY-PARASITIC-ROOT: supply channel added-mass ratio {mass_ratio:.4} gives a \
                 parasitic step root of magnitude {root:.6}, above {maximum_root:.6} (the -60 dB decay \
                 within one period of the fastest admitted mode). Shorten the passage, widen it, or use a \
                 lighter reed"
            ),
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
        /// Smallest cell count whose spatial screening estimate meets the
        /// tolerance. This does not admit material-interface eigenmode error.
        min_admissible_cells: usize,
        /// Spatial screening estimate at `cells` [cents]. Nonuniform interiors
        /// use the largest cell travel time relative to total regional travel
        /// time; their measured TMM pitch error remains the admission.
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
    /// Explicit, coherently resolved gas shared by ordinary excitation,
    /// propagation, loading and pressure observation. `None` resolves the
    /// scenario's ambient-air declaration. `baffled_duct` can instead bind
    /// individual interior regions and an independent exterior explicitly.
    /// Solid material temperature and prestress remain independent inputs.
    pub gas: Option<GasState>,
    /// Explicit observer absorption law, independent of gas-state provenance
    /// and the mechanical discretization. `None` preserves the scenario-air
    /// path when `gas` is absent, and selects Stokes-Kirchhoff for supplied gas
    /// in the ordinary assembly. Physical `baffled_duct` material overrides
    /// require an explicit law, including an explicit choice to omit relaxation.
    /// A supplied NASA-9 air state therefore needs an explicit ISO declaration
    /// to include that air-specific molecular-relaxation model.
    pub observer_absorption: Option<crate::air_path::GasAbsorptionModel>,
    /// The gas's molecular relaxation, evaluated at exactly the realized gas
    /// state. It is required by, and only by,
    /// [`crate::air_path::GasAbsorptionModel::MolecularRelaxation`]. A state
    /// evaluated for another gas, or a state supplied with another law, refuses
    /// instead of being ignored.
    pub observer_relaxation: Option<&'a fs_material::relaxation::RelaxationState>,
    /// Physical mouth microphones for a reed-driven duct with a flanged outlet.
    /// Its regional gas, exterior load and causal observation share the existing
    /// owners. A lumped ladder admits and advances the same original regional
    /// gas/material-wall frustum model, including for a single cone, and exports
    /// its terminal flow. Uniform-gas rigid-wall ladders also observe each native
    /// open/vented side mouth through its individual parallel-baffle path. Other coupled
    /// members and unsupported medium/wall combinations refuse explicitly.
    /// The observer law/relaxation above applies to this path's exterior gas;
    /// the receiver position replaces `assembly.listener.distance_m`.
    pub baffled_duct: Option<duct_observer::BaffledDuctObserver<'a>>,
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

/// Realized pressure, its observation binding and the gas state that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct RealizedAssembly {
    /// Sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Pressure [Pa], one sample per period. A bound `duct_observer` supplies
    /// exterior microphone pressure. Legacy duct images without that binding
    /// retain their source-time/local pressure diagnostics; sharing Pa units
    /// does not make those diagnostics an exterior acoustic observation.
    pub pressure_pa: Vec<f64>,
    /// Derived ambient state (rho, c, mu, …). On the physical duct path this is
    /// the exterior observer gas; `duct_observer` retains the interior regions.
    pub gas: GasState,
    /// Bore realization receipt when a reed bore was realized.
    pub bore_receipt: Option<BoreRealizationReceipt>,
    /// The observer absorption law applied, and any molecular relaxation it
    /// left out, with an upper bound on its effect at the listener.
    pub observer_absorption: crate::air_path::ObserverAbsorptionReceipt,
    /// Physical media, receiver and causal-delay binding when explicitly selected.
    pub duct_observer: Option<duct_observer::BaffledDuctObservation>,
}

/// Realize every described member and sum their observer pressures.
///
/// # Errors
/// Empty assemblies, invalid events, gas/domain refusals, stepper or TMM
/// refusals.
#[allow(clippy::too_many_lines)]// one coherent realization dispatcher
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

#[allow(clippy::too_many_lines)]// one coherent realization dispatcher
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
    if assembly.string.as_ref().is_some_and(|s| s.moving_end && s.polarization_detune > 0.0) {
        // One scalar body port cannot specify both transverse attachment
        // motions. Do not silently drop the second plane or constrain it to
        // the first plane's velocity.
        return Err(AcousticRealizeError::InvalidDescription {
            what: "moving-end second polarization requires explicit transverse attachment geometry and is not supported by this assembly",
        });
    }
    if let Some(observer) = realize_options.baffled_duct {
        return duct_observer::realize(assembly, chart.is_some(), realize_options, observer);
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
        && assembly
            .duct
            .as_ref()
            .is_some_and(|d| bore_spec(d).is_some());
    // Every assembly path consumes the same complete body image. A cavity is
    // an internal load of that body, so adding a duct must not replace the
    // moving-end constraint with the fixed-fixed, one-way bridge path.
    if let Some(cavity) = assembly.cavity {
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
            &mut bodies,
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
                &mut bodies,
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
#[allow(clippy::too_many_lines)]// one coherent coupled realization
#[allow(clippy::needless_range_loop)]// the sample index is the time axis
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
            1.0, dt,
        )?];
        if string.polarization_detune > 0.0 {
            members.push(kc_string_member(
                &twin,
                assembly.pluck,
                assembly.bow,
                &assembly.obstacles,
                gas,
                listener_m,
                0.85, dt,
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

#[allow(clippy::too_many_arguments)]// one coherent chain step record
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

#[allow(clippy::too_many_arguments)]// one coherent chain step record
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

#[allow(clippy::too_many_arguments)]// one coherent chain step record
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
#[allow(clippy::too_many_lines)]// one coherent coupled realization
#[allow(clippy::needless_range_loop)]// the sample index is the time axis
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
    let dt = 1.0 / f64::from(assembly.sample_rate_hz);
    let mut members = vec![kc_string_member(
        string,
        assembly.pluck,
        assembly.bow,
        &assembly.obstacles,
        gas,
        listener_m,
        1.0, dt,
    )?];
    if string.polarization_detune > 0.0 {
        members.push(kc_string_member(
            &twin,
            assembly.pluck,
            assembly.bow,
            &assembly.obstacles,
            gas,
            listener_m,
            0.85, dt,
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
    let law = selected_observer_law(assembly.ambient, gas, absorption.model);
    let relaxation_omitted = match law {
        GasAbsorptionModel::Iso9613Air { .. } | GasAbsorptionModel::MolecularRelaxation => None,
        GasAbsorptionModel::StokesKirchhoff | GasAbsorptionModel::Lossless => Some(relaxation_omission(
            assembly.listener.distance_m,
            0.5 * f64::from(assembly.sample_rate_hz),
        )?),
    };
    // Every law uses its original causal response, including ISO's molecular
    // phase and the classical law's declared minimum-phase realization.
    // The remaining receiver delay is returned, not applied: this is a
    // source-time history including FIR latency, not an absolute arrival.
    crate::air_path::propagate_pressure_history_with_model(
        &mut pressure_pa,
        dt,
        assembly.listener.distance_m,
        gas,
        law,
        absorption.relaxation,
    )?;
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa,
        gas: *gas,
        bore_receipt: None,
        observer_absorption: ObserverAbsorptionReceipt { law, relaxation_omitted },
        duct_observer: None,
    })
}

// Explicit constitutive selection is shared by offline and physical causal
// observers. Only a scenario-air input is eligible for the ambient default.
fn selected_observer_law(
    ambient: AmbientGas,
    gas: &GasState,
    selected: Option<crate::air_path::GasAbsorptionModel>,
) -> crate::air_path::GasAbsorptionModel {
    use crate::air_path::GasAbsorptionModel;
    selected.unwrap_or_else(|| {
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
/// delay is omitted. Linear held-force modes use their exact endpoint
/// acceleration and jerk. Fixed-end nonlinear modes differentiate their full
/// original storage law at the accepted endpoint under the same held input.
/// The joined moving-end path supplies
/// interval-mean accelerations from its accepted momentum increments, retaining
/// that clock instead of inferring endpoint reactions from a descriptor step's
/// interval multiplier. Distributional impulses at ideal force, release or
/// parameter jumps are not reconstructed by either observation.
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
                if string.moving_end {
                    let kpi = (k as f64 - 0.5) * pi;
                    let sign = if k % 2 == 1 { 1.0 } else { -1.0 };
                    // The original half-integer cosine shape has both signed
                    // area and first moment about the same string midpoint.
                    return (
                        factor * sign * string.length_m / kpi,
                        factor * string.length_m.powi(2)
                            * (0.5 * sign / kpi - 1.0 / kpi.powi(2))
                            / gas.sound_speed,
                    );
                }
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
        Self::check_pressure(pressure)?;
        self.previous_acceleration = Some(acceleration);
        Ok(pressure)
    }

    /// Smooth endpoint limit under this sample's held force and coefficients.
    /// From q'' + d q' + omega^2 q = f, a = f - omega^2 q - d v and
    /// j = -omega^2 v - d a. In particular the even sine modes radiate through
    /// j alone: differencing a would attenuate them by sinc(omega dt/2) and
    /// introduce a half-sample phase lag. No sampled history or allocation is
    /// needed, including after a material or gas rebind. An ideal jump's delta
    /// impulse is outside this smooth, held-parameter observation.
    fn observe_linear(
        &self,
        model: &ModalAcousticTimeModel,
        force: &[f64],
    ) -> Result<f64, AcousticRealizeError> {
        let pressure = self
            .weights
            .iter()
            .zip(model.modes())
            .zip(model.states())
            .zip(model.damping_rates_per_s())
            .zip(force)
            .map(|((((&(area, moment), mode), state), &damping), &f)| {
                let omega_sq = mode.angular_frequency_rad_s.powi(2);
                let acceleration = f
                    - omega_sq * state.displacement_m_sqrt_kg
                    - damping * state.velocity_m_sqrt_kg_per_s;
                let jerk = -omega_sq * state.velocity_m_sqrt_kg_per_s - damping * acceleration;
                area * acceleration + moment * jerk
            })
            .sum();
        Self::check_pressure(pressure)
    }

    fn check_pressure(pressure: f64) -> Result<f64, AcousticRealizeError> {
        let limit = ModalAcousticTimeBudget::audible_reference().maximum_abs_pressure_pa;
        if !pressure.is_finite() || pressure.abs() > limit {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "compact string observer pressure is nonfinite or exceeds its pressure budget",
            });
        }
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
    if string.bending_loss_factor.is_some()
        && (string.kelvin_voigt_bending.is_some() || string.relaxation_bending.is_some())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "a bending loss factor cannot be combined with another bending loss law",
        });
    }
    let tension = tension_loss_zeta(string, omega, wave_number)?;
    let factor_bending = bending_loss_zeta(string, omega, wave_number)?;
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
    // the Kelvin-Voigt or causal memory path above, or a constant loss factor.
    Ok(StringModeLossParts {
        authored: string.damping_ratio,
        material: tension + factor_bending,
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

/// Modal damping ratio of the bending stiffness's constant loss factor eta_B.
///
/// The bending term of Woodhouse, Plucked guitar transients, Acta Acustica
/// united with Acustica 90 (2004) 945-965, eq. (8): a complex flexural rigidity
/// EI (1 + i eta_B) loses eta_B of bending's share of the mode's stiffness per
/// radian, so zeta = eta_B EI k⁴ / (2 mu omega²) at the mode's own admitted
/// frequency, held as a modal viscous rate like eta_F. With
/// mu omega² = T k² + EI k⁴ this and [`tension_loss_zeta`] are that equation's
/// friction and bending terms over its one denominator. eta_B is constant in
/// frequency: Kelvin-Voigt's viscosity has a loss factor that grows with omega
/// and is a different law. It is material heat. `None` contributes nothing.
/// An authored ratio or Rayleigh pair is already a complete reduced loss, so it
/// cannot be combined with a physical eta_B.
fn bending_loss_zeta(
    string: &PrestressedString,
    omega: f64,
    wave_number: f64,
) -> Result<f64, AcousticRealizeError> {
    let Some(eta) = string.bending_loss_factor else {
        return Ok(0.0);
    };
    if !eta.is_finite() || eta < 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "bending loss factor must be finite and nonnegative",
        });
    }
    if string.rayleigh.is_some() || string.damping_ratio != 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "bending loss factor cannot be combined with authored internal/Rayleigh loss",
        });
    }
    let zeta = 0.5 * eta * string.bending_stiffness_n_m2 * wave_number.powi(4)
        / (string.lin_density_kg_m * omega * omega);
    if !zeta.is_finite() || zeta < 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "bending loss damping is unrepresentable",
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
    // Fixed-end members propagate the complete quadratic block, INCLUDING
    // these Maxwell coordinates, through PreparedNonlinearZoh. Imposing the
    // old midpoint limits there would reject audible partials and fast
    // material relaxation despite exact linear propagation. Moving-end joins
    // can still use midpoint for their coupled/nonlinear storage, so retain
    // that path's reference limits in addition to its own coupled-rate gate.
    // Neither policy bounds nonlinear/contact-generated frequencies.
    let reference_rate_limit = 0.1 * f64::from(sample_rate);
    if string.moving_end
        && law.branches.iter().any(|branch| {
            branch.relaxing_stiffness_n_m2 > 0.0
                && 1.0 / branch.relaxation_time_s > reference_rate_limit
        })
    {
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
            if string.moving_end && instant > reference_rate_limit {
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

/// Two-way Dirac realize: moving-end waveguide ⊕ complete plate/body image
/// ⊕ optional [`acoustic_chain`] duct. The image retains its cavity, material
/// and exterior histories together. `EA > 0` uses
/// [`kirchhoff_carrier_moving_end`] on the same free-fixed port.
/// Bow, obstacle stations, and blow/reed are leftover ports of
/// [`join_port`]. Von Karman is the same plate pHS with a quartic
/// storage. Hunt–Crossley and Stribeck stay port forces, not terms
/// in `H`.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)]// one coherent realization
fn realize_dirac_join(
    string: &PrestressedString,
    pluck: Option<Pluck>,
    bow: Option<BowStroke>,
    obstacles: &[UnilateralObstacle],
    texture: Option<ContactTexture>,
    plates: &mut PlateBank,
    duct: Option<&ViscothermalDuct>,
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    validate_string(string, pluck, bow)?;
    if duct.is_none() && (reed.is_some() || blow.is_some()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "a Dirac fluid exciter requires an independent duct source port",
        });
    }
    let n_s = string.n_modes.max(1);
    let zetas: Vec<f64> = (0..n_s)
        .map(|k| {
            let omega = moving_end_omega(string, k);
            let wave_number = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
            mode_zeta(string, omega, wave_number, gas)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let contacts = StationContact::prepare(obstacles)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let waveguide = moving_end_string_phs(string, &zetas, bow, &contacts)?;
    let n_string = waveguide.state_dim();
    let mut string_state = vec![0.0; n_string];
    if let Some(pluck) = pluck {
        for k in 0..n_s {
            string_state[2 * k] = free_fixed_pluck_modal(pluck, string, k);
        }
    }
    initialize_relaxed_bending(string, &mut string_state);
    let body = plates.physical_image(gas, sample_rate_hz, &fs_exec::CancelGate::default())?;
    let drive_inlet = reed.is_some() || blow.is_some();
    let (
        waveguide,
        sys,
        mut x,
        ode_start,
        body_observation,
        fluid_port,
    ) = match (body, duct) {
        (None, Some(line)) => {
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
            let fluid_port = drive_inlet.then_some(waveguide.port_dim() - 1);
            let system = transformer(waveguide, chain, 0, chain_face, area)
                .map_err(render_phs_error)?;
            string_state.resize(system.state_dim(), 0.0);
            (Some(system), None, string_state, 0, None, fluid_port)
        }
        (None, None) => (Some(waveguide), None, string_state, 1, None, None),
        (Some(plate), None) => {
            let joined = dirac_body::join_complete_body(waveguide, string_state, plate, None)?;
            (
                None,
                Some(joined.system),
                joined.state,
                0,
                Some(joined.observation),
                joined.fluid_port,
            )
        }
        (Some(plate), Some(line)) => {
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
            let load = dirac_body::DiracBodyLoad {
                state: vec![0.0; chain.state_dim()],
                system: chain,
                face: chain_face,
                drive_port: drive_inlet.then_some(0),
            };
            let joined = dirac_body::join_complete_body(waveguide, string_state, plate, Some(load))?;
            (
                None,
                Some(joined.system),
                joined.state,
                0,
                Some(joined.observation),
                joined.fluid_port,
            )
        }
    };
    let dt = 1.0 / f64::from(sample_rate_hz);
    let (sys, mut descriptor_exact) = match sys {
        Some(system) if system.storage().quadratic_part().is_some_and(|part| part.exact) => {
            let budget = fs_phs::linear_exact::LinearSplitBudget {
                max_modes: system.state_dim(),
                max_omega_dt: core::f64::consts::PI,
            };
            let prepared = fs_phs::linear_exact::PreparedCommonFlowZoh::new_owned(
                system, dt, budget,
            )
            .map_err(render_phs_error)?;
            (None, Some(prepared))
        }
        original => (original, None),
    };
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
    let waveguide_exact = match (&waveguide, &waveguide_linear) {
        (Some(system), Some(linear)) => prepare_chain_exact(system, linear, dt)?,
        _ => None,
    };
    let mut texture = TextureDrive::try_new(texture)?;
    let mut out = Vec::with_capacity(n);
    let mut rec_x = vec![0.0; x.len()];
    let mut average_output = vec![0.0; descriptor_exact.as_ref().map_or(0, |s| s.port_dim())];
    let mut string_observer = StringObserver::new(string, gas, listener_m, 1.0);
    for i in 0..n {
        let t = i as f64 * dt;
        if let Some(prepared) = descriptor_exact.as_mut() {
            let u = leftover_u(
                prepared.system(),
                &x,
                n_s,
                string,
                bow,
                &contacts,
                reed,
                blow,
                gas,
                t,
                dt,
                &mut texture,
                fluid_port,
            )?;
            prepared.step_into(&x, &u, &mut rec_x, &mut average_output)
                .map_err(render_phs_error)?;
        } else if let Some(s) = &sys {
            let u = leftover_u(
                s,
                &x,
                n_s,
                string,
                bow,
                &contacts,
                reed,
                blow,
                gas,
                t,
                dt,
                &mut texture,
                fluid_port,
            )?;
            rec_x = step_descriptor(s, &x, &u, dt)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
                .x;
        } else {
            let wg = waveguide.as_ref().expect("waveguide-only clock");
            let u = leftover_u_ode(
                wg,
                &x,
                n_s,
                string,
                bow,
                &contacts,
                reed,
                blow,
                gas,
                t,
                dt,
                &mut texture,
                ode_start,
                fluid_port,
            )?;
            rec_x = step_chain_flow(
                wg,
                waveguide_linear.as_ref().expect("prepared waveguide clock"),
                waveguide_exact.as_ref(),
                &x,
                &u,
                dt,
            )?
            .x;
        }
        let acceleration: Vec<f64> = (0..n_s)
            .map(|k| (rec_x[2 * k + 1] - x[2 * k + 1]) / dt)
            .collect();
        // Keep the original accepted interval-mean source. Its first interval
        // has no preceding mean: prime with itself, without a fabricated jerk
        // or release impulse. All later means share the same observer history.
        let initial = string_observer.previous_acceleration.is_none()
            .then(|| acceleration.clone());
        let mut p = string_observer.observe(acceleration, initial, dt)?;
        if let Some(observation) = &body_observation {
            p += observation.pressure_between(&x, &rec_x, dt, gas.density, listener_m)?;
        }
        if body_observation.is_none()
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
        core::mem::swap(&mut x, &mut rec_x);
        out.push(p);
    }
    Ok(out)
}

fn moving_end_string_phs(
    string: &PrestressedString,
    zetas: &[f64],
    bow: Option<BowStroke>,
    contacts: &[StationContact],
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    let n = string.n_modes.max(1);
    let phi0 = det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m));
    let m = 1 + usize::from(bow.is_some()) + contacts.len();
    let mut g = vec![0.0; (2 * n) * m];
    for k in 0..n {
        let kappa = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
        g[(2 * k + 1) * m] = phi0;
        let mut col = 1;
        if let Some(bow) = bow {
            g[(2 * k + 1) * m + col] = phi0 * det::cos(kappa * bow.station_frac * string.length_m);
            col += 1;
        }
        for contact in contacts {
            g[(2 * k + 1) * m + col] = phi0 * det::cos(kappa * contact.station_frac() * string.length_m);
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

#[allow(clippy::too_many_arguments)]// one coherent solver probe
#[allow(clippy::unnecessary_wraps)]// uniform Result surface with the ODE twin
fn leftover_u(
    sys: &fs_phs::DescriptorPortHamiltonian,
    x: &[f64],
    n_s: usize,
    string: &PrestressedString,
    bow: Option<BowStroke>,
    contacts: &[StationContact],
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
    texture: &mut TextureDrive,
    fluid_port: Option<usize>,
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
    col = moving_end_contact_inputs(string, x, n_s, contacts, &y, &mut u, col)?;
    if let Some(port) = fluid_port {
        if port != col || port >= m {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "Dirac fluid source must follow the original string force ports",
            });
        }
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

#[allow(clippy::too_many_arguments)]// one coherent solver probe
#[allow(clippy::unnecessary_wraps)]// uniform Result surface with the DAE twin
fn leftover_u_ode(
    sys: &fs_phs::PortHamiltonian,
    x: &[f64],
    n_s: usize,
    string: &PrestressedString,
    bow: Option<BowStroke>,
    contacts: &[StationContact],
    reed: Option<BeatingReed>,
    blow: Option<VolumeVelocityPulse>,
    gas: &GasState,
    t: f64,
    dt: f64,
    texture: &mut TextureDrive,
    start_col: usize,
    fluid_port: Option<usize>,
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
    col = moving_end_contact_inputs(string, x, n_s, contacts, &y, &mut u, col)?;
    if let Some(port) = fluid_port {
        if port != col || port >= m {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "ODE fluid source must follow the original string force ports",
            });
        }
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

// Both the free/duct ODE and common-velocity descriptor use these same physical
// force ports, including original quadrature and the actual conjugate velocity.
#[allow(clippy::too_many_arguments)]
fn moving_end_contact_inputs(
    string: &PrestressedString, x: &[f64], n_s: usize, contacts: &[StationContact],
    velocities: &[f64], input: &mut [f64], start: usize,
) -> Result<usize, AcousticRealizeError> {
    let end = start.checked_add(contacts.len()).ok_or(AcousticRealizeError::InvalidDescription {
        what: "moving-end contact port count overflows",
    })?;
    if n_s.checked_mul(2).is_none_or(|n| n > x.len())
        || end > input.len() || velocities.len() != input.len()
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "moving-end contact lost its original string state or force/velocity ports",
        });
    }
    for (i, contact) in contacts.iter().enumerate() {
        let opening = moving_end_opening(string, x, n_s, contact.station_frac());
        input[start + i] = contact.force(opening, velocities[start + i])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    }
    Ok(end)
}

fn moving_end_opening(string: &PrestressedString, x: &[f64], n_s: usize, station: f64) -> f64 {
    let phi0 = det::sqrt(2.0 / (string.lin_density_kg_m * string.length_m));
    let mut y = 0.0;
    for k in 0..n_s {
        let kappa = (k as f64 + 0.5) * core::f64::consts::PI / string.length_m;
        let phi = phi0 * det::cos(kappa * station * string.length_m);
        y += x[2 * k] * phi;
    }
    y
}

fn moving_end_omega(string: &PrestressedString, k_zero: usize) -> f64 {
    // This image owns its half-integer cosine basis; the material dispersion
    // remains the same law used by pinned strings, fibers and impact wires.
    fs_nlmodal::prestressed_beam_omega_at_wavenumber(
        (k_zero as f64 + 0.5) * core::f64::consts::PI / string.length_m,
        string.tension_n,
        string.lin_density_kg_m,
        string.bending_stiffness_n_m2,
    )
}

#[cfg(test)]
mod moving_end_dispersion_tests {
    use super::*;

    // G1/G3: both boundary images substitute their own spatial eigenvalue
    // in the same tension+bending operator, including the bending limit.
    #[test]
    fn moving_end_and_pinned_images_share_constitutive_dispersion() {
        let mut string = PrestressedString {
            length_m: 0.65,
            tension_n: 70.0,
            lin_density_kg_m: 0.004,
            axial_stiffness_n: 0.0,
            width_m: 0.001,
            n_modes: 8,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.006,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: true,
            tension_loss_factor: None,
            bending_loss_factor: None,
        };
        for tension in [70.0, 1.0e-300, 0.0] {
            string.tension_n = tension;
            for k_zero in 0..string.n_modes {
                let kappa = (k_zero as f64 + 0.5) * core::f64::consts::PI / string.length_m;
                let expected = ((tension * kappa.powi(2)
                    + string.bending_stiffness_n_m2 * kappa.powi(4))
                    / string.lin_density_kg_m).sqrt();
                let actual = moving_end_omega(&string, k_zero);
                assert!(actual.is_finite());
                assert!((actual - expected).abs() <= 4.0e-14 * expected);
                // A pinned member of twice this span has the same odd spatial
                // eigenvalues; this compares the two public consumer paths.
                let mut pinned = string.clone();
                pinned.length_m *= 2.0;
                pinned.moving_end = false;
                assert_eq!(actual, string_mode_omega(&pinned, 2 * k_zero + 1));
            }
        }
        string.bending_stiffness_n_m2 = f64::NAN;
        assert!(moving_end_omega(&string, 0).is_nan());
    }
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

#[allow(clippy::too_many_arguments)]// one coherent realization record
#[allow(clippy::needless_range_loop)]// modal index spans state and shape arrays
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
    settled_load_initialized: bool,
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
    /// Dissipation inside the solid: Kelvin-Voigt bending and the constant
    /// tension and bending loss factors [J].
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
            settled_load_initialized: false,
        })
    }

    /// Source receipts and physical coefficients at the accepted epoch.
    #[must_use]
    pub const fn specimen(&self) -> &crate::string_specimen::ResolvedStringSpecimen {
        &self.specimen
    }

    /// One epoch covers the specimen, operators, vibration and observer.
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
    /// prior load initialization can be initialized; existing motion is never
    /// reset. Material state and sample count are unchanged. Publication
    /// advances the epoch and records the one-time initialization.
    /// Cancellation or a force/state/energy refusal leaves the runtime untouched.
    pub fn initialize_point_load_equilibrium(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        loads: &[StringPointForce],
    ) -> Result<u64, AcousticRealizeError> {
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        if self.accepted_samples != 0
            || self.settled_load_initialized
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
        cx.checkpoint()
            .map_err(|_| AcousticRealizeError::Cancelled)?;
        self.member = candidate;
        self.settled_load_initialized = true;
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
    /// material receipts, modal coordinates and settled-load initialization survive.
    /// The update advances the epoch, but consumes no sample or material heat.
    /// This prescribed gas change does not model fluid-load…33935 tokens truncated…, +185 at 12, -25 at 24, -12 at 48).
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
        use crate::air_path::{GasAbsorptionModel, propagate_pressure_history_with_model};
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
        propagate_pressure_history_with_model(
            &mut expected,
            1.0 / 16000.0,
            description.listener.distance_m,
            &gas,
            model,
            None,
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
        propagate_pressure_history_with_model(
            &mut original,
            1.0 / 16000.0,
            description.listener.distance_m,
            &legacy,
            scenario_law,
            None,
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

    /// G3: the public offline default and a supplied gas's classical observer
    /// produce the same phased pressure as the existing streaming law.
    #[test]
    fn offline_observer_defaults_use_the_streaming_complex_transfer() {
        use crate::air_path::{GasAbsorptionModel, ObserverFilter, absorb_pressure_history_with_model};
        let mut description = fixture();
        description.listener.distance_m = 30.0;
        let (air_gas, air) = legacy_admission(&description);
        let iso = GasAbsorptionModel::Iso9613Air {
            relative_humidity: description.ambient.relative_humidity, air,
        };
        let history: Vec<_> = (0..513).map(|i| det::sin(i as f64 * 0.29)).collect();
        for (gas, selected, expected_law) in [
            (air_gas, None, iso),
            (helium_like(), Some(GasAbsorptionModel::StokesKirchhoff), GasAbsorptionModel::StokesKirchhoff),
        ] {
            let realized = finish_observer(&description, &gas, history.clone(), selected.into()).unwrap();
            let mut filter = ObserverFilter::for_law(
                expected_law, &gas, None, description.listener.distance_m, description.sample_rate_hz,
            ).unwrap();
            let expected: Vec<_> = history.iter().map(|&x| filter.process(x)).collect();
            assert!(realized.pressure_pa.iter().zip(&expected).all(|(a, b)| (a - b).abs() < 1e-12));
            assert_eq!(realized.observer_absorption.law, expected_law);
            let mut magnitude_only = history.clone();
            absorb_pressure_history_with_model(&mut magnitude_only, 1.0 / f64::from(description.sample_rate_hz),
                description.listener.distance_m, &gas, expected_law).unwrap();
            assert!(realized.pressure_pa.iter().zip(&magnitude_only).any(|(a, b)| (a - b).abs() > 1e-6));
        }
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
    /// gas refuses. The source-time history uses that law's causal response.
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
        crate::air_path::propagate_pressure_history_with_model(&mut reference, 1.0 / 16000.0,
            description.listener.distance_m, &supplied, GasAbsorptionModel::StokesKirchhoff, None).unwrap();
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
                bending_loss_factor: None,
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
            bending_loss_factor: None,
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

#[cfg(test)]
mod bending_loss_tests {
    use super::*;
    use core::f64::consts::{PI, TAU};
    use fs_scenario::acoustic::{BendingRelaxationBranch, KelvinVoigtBending, PronyBending};

    /// Woodhouse 2004 Table I, the E2 string, on the guitar's 0.68 m pinned span.
    fn e2(eta_f: Option<f64>, eta_b: Option<f64>) -> PrestressedString {
        PrestressedString {
            length_m: 0.68,
            tension_n: 71.6,
            lin_density_kg_m: 0.00624,
            axial_stiffness_n: 0.0,
            width_m: 0.001118,
            n_modes: 24,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 57.0e-6,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor: eta_f,
            bending_loss_factor: eta_b,
        }
    }

    /// (mode, omega, zeta_F, zeta_B) of Woodhouse's eq. (8), p. 956, as zeta = eta / 2.
    /// Worked from the printed form in 50-digit decimal by a separate script, not
    /// by this crate, for eta_F 2e-5 and eta_B 0.02.
    const E2_LAW: [(usize, f64, f64, f64); 3] = [
        (1, 4.948_902_138_948_980_65e2, 9.999_830_083_179_121_68e-6, 1.699_168_208_783_223_04e-7),
        (12, 5.945_893_126_065_888_00e3, 9.975_591_286_343_197_58e-6, 2.440_871_365_680_241_96e-5),
        (24, 1.193_524_631_775_046_69e4, 9.903_074_890_545_638_05e-6, 9.692_510_945_436_195_23e-5),
    ];

    fn wave_number(string: &PrestressedString, j: usize) -> f64 {
        j as f64 * PI / string.length_m
    }

    #[test]
    fn bending_loss_is_eq_8s_term_and_moves_nothing_else() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let (plain, friction) = (e2(None, None), e2(Some(2.0e-5), None));
        let (both, zero) = (e2(Some(2.0e-5), Some(0.02)), e2(None, Some(0.0)));
        for (j, omega, zeta_f, zeta_b) in E2_LAW {
            let (actual, k) = (string_mode_omega(&both, j), wave_number(&both, j));
            assert!((actual / omega - 1.0).abs() < 1.0e-14, "mode {j}: omega {actual:e}");
            let parts = |string: &PrestressedString| mode_loss_parts(string, actual, k, &gas).unwrap();
            let (base, tension, full) = (parts(&plain), parts(&friction), parts(&both));
            assert!((tension.material / zeta_f - 1.0).abs() < 1.0e-13, "mode {j}");
            assert!((full.material / (zeta_f + zeta_b) - 1.0).abs() < 1.0e-13, "mode {j}");
            // The air stays the one computed term: this law adds nothing to it.
            assert_eq!(full.air.to_bits(), base.air.to_bits());
            assert_eq!(full.authored.to_bits(), base.authored.to_bits());
            // A zero factor is the absent one, bit for bit.
            assert_eq!(parts(&zero).total().to_bits(), base.total().to_bits());
        }
    }

    #[test]
    fn a_kelvin_voigt_viscosity_matched_at_one_mode_is_not_this_law_at_another() {
        // Planted negative. The viscosity that equals eta_B at mode 12 has a loss
        // factor growing with omega, so it overshoots mode 24 by the frequency ratio.
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let (_, omega_12, _, zeta_12) = E2_LAW[1];
        let (_, omega_24, _, zeta_24) = E2_LAW[2];
        let viscous = PrestressedString {
            kelvin_voigt_bending: Some(KelvinVoigtBending {
                viscous_stiffness_n_m2_s: 0.02 * 57.0e-6 / omega_12,
                omega_band_rad_s: (0.0, 1.0e6),
                material_state_identity: None,
            }),
            ..e2(None, None)
        };
        let rate = |j, omega| mode_loss_parts(&viscous, omega, wave_number(&viscous, j), &gas).unwrap().material;
        assert!((rate(12, omega_12) / zeta_12 - 1.0).abs() < 1.0e-12);
        let overshoot = rate(24, omega_24) / zeta_24;
        assert!((overshoot / (omega_24 / omega_12) - 1.0).abs() < 1.0e-12, "overshoot {overshoot}");
        assert!(overshoot > 1.9);
    }

    #[test]
    fn bending_loss_refuses_other_bending_laws_reduced_losses_and_nonphysical_factors() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let refused = |string: &PrestressedString, expected: &str| {
            let (omega, k) = (string_mode_omega(string, 3), wave_number(string, 3));
            match mode_zeta(string, omega, k, &gas) {
                Err(AcousticRealizeError::InvalidDescription { what }) => assert_eq!(what, expected),
                other => panic!("admitted or misreported {other:?}, wanted {expected}"),
            }
        };
        let kelvin_voigt = Some(KelvinVoigtBending {
            viscous_stiffness_n_m2_s: 2.0e-8,
            omega_band_rad_s: (0.0, 1.0e6),
            material_state_identity: None,
        });
        let prony = Some(PronyBending {
            branches: vec![BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 1.0e-6,
                relaxation_time_s: 1.0e-3,
            }],
            omega_band_rad_s: (0.0, 1.0e6),
            material_state_identity: None,
            source_properties: None,
        });
        // Two bending laws would book one stiffness's loss twice, whatever the factor.
        for eta in [0.02, 0.0] {
            let other = "a bending loss factor cannot be combined with another bending loss law";
            refused(&PrestressedString { kelvin_voigt_bending: kelvin_voigt, ..e2(None, Some(eta)) }, other);
            refused(&PrestressedString { relaxation_bending: prony.clone(), ..e2(None, Some(eta)) }, other);
        }
        let reduced = "bending loss factor cannot be combined with authored internal/Rayleigh loss";
        refused(&PrestressedString { damping_ratio: 1.0e-3, ..e2(None, Some(0.02)) }, reduced);
        let rayleigh = Some(RayleighParams { alpha_per_s: 1.0, beta_s: 0.0 });
        refused(&PrestressedString { rayleigh, ..e2(None, Some(0.02)) }, reduced);
        for bad in [f64::NAN, f64::INFINITY, -1.0e-6] {
            refused(&e2(None, Some(bad)), "bending loss factor must be finite and nonnegative");
        }
    }

    #[test]
    fn the_bending_term_is_booked_as_material_heat_and_the_ledger_fails_without_it() {
        let gas = gas_state(AmbientGas::sea_level()).unwrap();
        let (lossy, omitted) = (e2(Some(2.0e-5), Some(0.02)), e2(Some(2.0e-5), None));
        let quiet = vec![0.0; lossy.n_modes];
        // The same admitted frequencies, booked as if the bending term did not exist.
        let without_term: Vec<_> = (1..=lossy.n_modes)
            .map(|k| {
                mode_loss_parts(&omitted, string_mode_omega(&lossy, k), wave_number(&lossy, k), &gas)
                    .unwrap()
            })
            .collect();
        let (rate, samples) = (48_000_u32, 9_600_u32);
        let time = f64::from(samples) / f64::from(rate);
        for (j, omega, zeta_f, zeta_b) in E2_LAW {
            // The drag's printed form, written out here rather than called.
            let radius = 0.5 * lossy.width_m;
            let resistance = TAU
                * gas.dynamic_viscosity
                * (1.0 + radius * det::sqrt(2.0 * gas.density * omega / gas.dynamic_viscosity));
            let zeta_air = resistance / (2.0 * lossy.lin_density_kg_m * omega);
            let zeta = zeta_f + zeta_b + zeta_air;

            let mut member = linear_string_member(&lossy, None, None, &gas, 1.0, rate, 1.0).unwrap();
            let mut states = member.model.states().to_vec();
            let start = 1.0e-3;
            states[j - 1].displacement_m_sqrt_kg = start;
            member.model.restore_states(&states).unwrap();
            let initial = 0.5 * omega * omega * start * start;

            let (mut heat, mut air, mut heat_without, mut air_without) = (0.0, 0.0, 0.0, 0.0);
            let mut energy = initial;
            for _ in 0..samples {
                let frame = member.model.step(&quiet).unwrap();
                let booked = partition_string_dissipation(&frame, &member.loss_parts).unwrap();
                assert_eq!(booked.authored_loss_j, 0.0);
                heat += booked.material_heat_j;
                air += booked.air_loss_j;
                let short = partition_string_dissipation(&frame, &without_term).unwrap();
                heat_without += short.material_heat_j;
                air_without += short.air_loss_j;
                energy = frame.total_modal_energy_j;
            }

            // The motion itself carries the term: exact free decay at the full rate.
            let damped = omega * det::sqrt(1.0 - zeta * zeta);
            let decay = det::exp(-zeta * omega * time);
            let (sine, cosine) = (det::sin(damped * time), det::cos(damped * time));
            let displacement = start * decay * (cosine + zeta * omega / damped * sine);
            let velocity = -start * decay * omega * omega / damped * sine;
            let expected = 0.5 * velocity * velocity + 0.5 * omega * omega * displacement * displacement;
            assert!((energy / expected - 1.0).abs() < 1.0e-8, "mode {j}: {energy:e} vs {expected:e}");
            let drop = initial - energy;
            assert!(drop > 0.0, "mode {j} lost nothing");

            // Sources and sinks: the drop must equal the booked heat plus the air's
            // share taken from the constants above, not from the booking under test.
            let air_share = zeta_air / zeta;
            let residual = (heat + drop * air_share - drop) / drop;
            assert!(residual.abs() < 1.0e-9, "mode {j}: ledger open by {residual:e}");
            assert!((air / drop / air_share - 1.0).abs() < 1.0e-9, "mode {j}");

            // Planted negative: left out of the booking, the same ledger stays open
            // by the term's whole share of the heat.
            let short = (heat_without + drop * air_share - drop) / drop;
            let share_lost = (zeta_f + zeta_b) / zeta - zeta_f / (zeta_f + zeta_air);
            assert!((short / -share_lost - 1.0).abs() < 1.0e-9, "mode {j}: {short:e}");
            assert!(short.abs() > 1.0e-3, "mode {j}: omission is invisible at {short:e}");
            // Comparing the booking only with itself cannot see that: the shares
            // always sum to one, so the omitted term moves into air and still closes.
            assert!(((heat_without + air_without) / drop - 1.0).abs() < 1.0e-9, "mode {j}");
        }
    }
}
