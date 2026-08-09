//! Physical structural modes and contact-force participation for resolved
//! three-dimensional specimens.
//!
//! This module is an integration boundary, not an object-specific sound
//! preset. Geometry comes from [`crate::specimen`], elastic properties come
//! from the evidence-bearing `fs-material` state point, the body-fitted mesh
//! comes from `fs-mesh`, `(K,M)` comes from `fs-solid`, and eigenpairs come
//! from `fs-modal`. A material name never selects a frequency, decay time, or
//! gain.
//!
//! The retained modes are mass normalized: `phi^T M phi = 1`. Consequently a
//! nodal mode-shape component has units `kg^-1/2`, while projecting a physical
//! point force produces a generalized modal force in `N kg^-1/2`. These units
//! matter later when acoustic radiation maps modal velocity to pressure.
//!
//! Applicability is explicit. The current volume producer admits sharp and
//! circularly filleted solid cylinders. Other already-supported Euler profile
//! families refuse until a conforming volume producer for their exact
//! meridians exists. Small-strain elasticity also refuses by construction once
//! thermal softening, yield, finite strain, phase change, or evolving topology
//! requires a higher constitutive rung.

use fs_bem::BemError;
use fs_bem::helmholtz::{
    DirectivityTable, Formulation as HelmholtzFormulation, HelmholtzError, MAX_SH_DEGREE, Medium,
    directivity_sh_table, far_field, solve_radiation,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::{ContentHash, DomainHasher};
use fs_couple::modal_acoustic_time::{
    ModalAcousticFrame, ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel,
};
use fs_exec::Cx;
use fs_material::gas::GasState;
use fs_material::visco::{LoweredModel, RayleighDamping, ViscoError, loss_factor_to_zeta};
use fs_math::c64::C64;
use fs_math::det;
use fs_mbd::{Pose, Vec3};
use fs_mesh::{
    RoundedCylinderMeshError, RoundedCylinderMeshSpec, RoundedCylinderTetMesh,
    rounded_cylinder_tet_mesh,
};
use fs_modal::{ModalError, SliceOptions, SliceStats, slice_window};
use fs_rep_frep::SquatDiscEdgeTreatment;
use fs_solid::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetLinearElasticProblem,
    TetMaterialField,
};

use crate::specimen::{DiscProfileSpec, ResolvedElasticDiscProfile};
use crate::timeline_resampling::{EventEvaluationSide, TimelineResampler, TimelineResamplingError};
use crate::{ChannelControl, EulerControlStream};

/// Schema version of the integrated structural modal artifact.
pub const STRUCTURAL_MODAL_BASIS_SCHEMA_VERSION: u32 = 1;
/// Maximum number of simultaneous physical pressure observers in one pass.
pub const MAX_PHYSICAL_PRESSURE_OBSERVERS: usize = 64;
const STRUCTURAL_MODAL_BASIS_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-modal-basis.v1";
const MODAL_ACOUSTIC_RADIATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.modal-acoustic-radiation.v1";
const MODAL_ACOUSTIC_DIRECTIVITY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.modal-acoustic-directivity.v1";

/// Resolution and resource controls, deliberately separate from physical
/// dimensions so callers cannot accidentally restate geometry differently
/// for mechanics and sound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StructuralMeshControls {
    /// Radial intervals from the axis through the planar-cap region.
    pub core_radial_segments: u32,
    /// Radial intervals through each outer circular fillet. Must be zero for
    /// a sharp rim and positive for a circular fillet.
    pub fillet_radial_segments: u32,
    /// Periodic angular intervals.
    pub azimuthal_segments: u32,
    /// Axial intervals.
    pub axial_segments: u32,
    /// Maximum admitted vertices.
    pub maximum_vertices: usize,
    /// Maximum admitted tetrahedra.
    pub maximum_tetrahedra: usize,
}

impl StructuralMeshControls {
    /// A bounded starting point for modal refinement studies.
    #[must_use]
    pub const fn modal_default(has_fillet: bool) -> Self {
        Self {
            core_radial_segments: 6,
            fillet_radial_segments: if has_fillet { 3 } else { 0 },
            azimuthal_segments: 24,
            axial_segments: 2,
            maximum_vertices: 100_000,
            maximum_tetrahedra: 600_000,
        }
    }
}

/// One physical modal-basis request.
pub struct StructuralModeRequest<'a> {
    /// Resolved geometry and complete material state.
    pub specimen: &'a ResolvedElasticDiscProfile,
    /// Discretization and work envelope.
    pub mesh: StructuralMeshControls,
    /// Strictly positive lower edge of the requested band [Hz].
    pub minimum_frequency_hz: f64,
    /// Upper edge of the requested band [Hz].
    pub maximum_frequency_hz: f64,
    /// Maximum number of modes the caller will retain.
    pub maximum_modes: usize,
    /// Certified sparse modal-slice controls.
    pub slice: SliceOptions,
    /// Element assembly and quality envelope.
    pub assembly: TetAssemblyBudget,
}

/// One certified, mass-normalized three-dimensional elastic mode.
#[derive(Clone, Debug)]
pub struct StructuralMode {
    /// Eigenvalue `lambda = omega^2` [s^-2].
    pub eigenvalue_s2: f64,
    /// Natural angular frequency [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Natural frequency [Hz].
    pub frequency_hz: f64,
    /// Certified eigenvalue interval [s^-2].
    pub eigenvalue_interval_s2: (f64, f64),
    /// Residual-derived eigenvalue distance bound [s^-2].
    pub eigenvalue_residual_s2: f64,
    /// Three displacement components at each volume node [kg^-1/2].
    pub nodal_shape_per_sqrt_kg: Vec<[f64; 3]>,
    /// Outward-normal displacement sampled at each boundary triangle
    /// centroid by P1 interpolation [kg^-1/2].
    pub panel_normal_shape_per_sqrt_kg: Vec<f64>,
}

/// Complete geometry/mechanics/modal artifact used by contact and acoustics.
#[derive(Clone, Debug)]
pub struct StructuralModalBasis {
    /// Integrated artifact schema.
    pub schema_version: u32,
    /// Identity of the resolved geometry plus material state.
    pub specimen_identity: ContentHash,
    /// Exact geometry/density profile identity shared with trajectory
    /// metadata, before additional elastic-state properties are attached.
    pub profile_identity: ContentHash,
    /// Exact evidence-bearing material state used to assemble `(K,M)`.
    pub material_state_identity: ContentHash,
    /// Exact body-fitted volume and boundary panelization.
    pub mesh: RoundedCylinderTetMesh,
    /// Assembled physical mass and stiffness operators.
    pub assembly: TetElasticAssembly,
    /// Certified in-band modes, ascending by frequency.
    pub modes: Vec<StructuralMode>,
    /// Certified requested eigenvalue window [s^-2].
    pub eigenvalue_window_s2: (f64, f64),
    /// Inertia-certified number of modes in the window.
    pub certified_mode_count: usize,
    /// Sparse eigensolver work accounting.
    pub slice_stats: SliceStats,
    /// Content identity binding all physical and numerical inputs and outputs.
    pub identity: ContentHash,
}

/// Projection of a physical point force onto every retained structural mode.
#[derive(Clone, Debug)]
pub struct PointForceProjection {
    /// Boundary triangle selected by closest-point distance.
    pub boundary_triangle: usize,
    /// Closest point on the piecewise-planar boundary [m].
    pub closest_point_m: [f64; 3],
    /// Distance from requested application point to that boundary [m].
    pub distance_to_boundary_m: f64,
    /// Barycentric weights on the selected triangle.
    pub barycentric: [f64; 3],
    /// Generalized force per retained mode [N kg^-1/2].
    pub modal_force_n_per_sqrt_kg: Vec<f64>,
}

/// An evidence-bound gas state used by the exterior acoustic solve.
///
/// `GasState` carries every derived transport/acoustic scalar but not the
/// originating species/model identity, so the latter remains explicit rather
/// than being guessed from density or sound speed.
#[derive(Clone, Copy, Debug)]
pub struct ResolvedAcousticMedium<'a> {
    /// Complete thermodynamic/transport state.
    pub gas: &'a GasState,
    /// Identity of the gas species, EOS, and transport model that produced it.
    pub gas_model_identity: ContentHash,
}

/// One microphone location expressed in the undeformed specimen body frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticObserver {
    /// Microphone point relative to the specimen origin [m].
    pub position_m: [f64; 3],
}

/// One microphone location in the inertial world frame.
///
/// Unlike [`AcousticObserver`], this point does not rotate with the radiating
/// body. It is the appropriate observation frame for a real stationary
/// microphone while a rigid specimen tumbles beneath it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticWorldObserver {
    /// Microphone position in world coordinates [m].
    pub position_world_m: [f64; 3],
}

/// Accuracy controls for a reusable body-frame acoustic directivity field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticDirectivityControls {
    /// Maximum spherical-harmonic degree retained in the far-field pattern.
    pub maximum_spherical_harmonic_degree: usize,
    /// Minimum admitted fraction of quadrature-estimated far-field power.
    pub minimum_captured_fraction: f64,
}

/// Exterior radiation of one mass-normalized mode, independent of observer.
#[derive(Clone, Debug)]
pub struct AcousticModeDirectivity {
    /// Index into [`StructuralModalBasis::modes`].
    pub structural_mode: usize,
    /// Evaluation angular frequency [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Wavenumber [1/m].
    pub wavenumber_rad_m: f64,
    /// Boundary-integral formulation used for this mode.
    pub formulation: HelmholtzFormulation,
    /// Body-frame far-field amplitude `F(direction)` where
    /// `p = F exp(i k r) / r` per unit generalized modal velocity.
    pub directivity: DirectivityTable,
    /// Radiated power per squared generalized modal velocity
    /// `[W s^2 / (m^2 kg)]`.
    pub radiated_power_per_modal_velocity_squared: f64,
    /// BEM panels per wavelength.
    pub panels_per_wavelength: f64,
    /// Probe-based lower bound on the BEM matrix condition number.
    pub condition_lower_bound: f64,
    /// Minimum distance admitted by the far-field approximation [m].
    pub minimum_far_field_distance_m: f64,
}

/// Observer-independent, gas-state-dependent modal radiation field.
///
/// The table is solved once in the undeformed body frame. Any world-fixed
/// observer can then be evaluated against the current rigid pose without
/// repeating the boundary-integral solve.
#[derive(Clone, Debug)]
pub struct ModalAcousticDirectivity {
    /// Structural artifact consumed by this solve.
    pub structural_basis_identity: ContentHash,
    /// Gas model/species identity supplied by the caller.
    pub gas_model_identity: ContentHash,
    /// Temperature [K].
    pub temperature_k: f64,
    /// Absolute pressure [Pa].
    pub ambient_pressure_pa: f64,
    /// Derived acoustic density [kg/m3].
    pub density_kg_m3: f64,
    /// Derived sound speed [m/s].
    pub sound_speed_m_s: f64,
    /// One reusable directivity table per retained structural mode.
    pub modes: Vec<AcousticModeDirectivity>,
    /// Identity binding structure, gas state/model, controls, and tables.
    pub identity: ContentHash,
}

/// Exterior radiation of one mass-normalized structural mode at its natural
/// frequency.
#[derive(Clone, Debug)]
pub struct AcousticModeRadiation {
    /// Index into [`StructuralModalBasis::modes`].
    pub structural_mode: usize,
    /// Evaluation angular frequency [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Boundary-integral formulation selected from nondimensional acoustic
    /// size. Plain CBIE avoids the documented low-`ka` Burton--Miller
    /// resistance artifact; Burton--Miller protects the higher-frequency arm
    /// from fictitious interior resonances.
    pub formulation: HelmholtzFormulation,
    /// Complex pressure per unit generalized modal velocity at the observer
    /// `[Pa s / (m sqrt(kg))]`, under the shared `exp(-i omega t)` convention.
    pub observer_pressure_per_modal_velocity: C64,
    /// Radiated power per squared generalized modal velocity
    /// `[W s^2 / (m^2 kg)]`.
    pub radiated_power_per_modal_velocity_squared: f64,
    /// BEM panels per wavelength.
    pub panels_per_wavelength: f64,
    /// Probe-based lower bound on the BEM matrix condition number.
    pub condition_lower_bound: f64,
    /// Minimum distance required by the declared far-field approximation [m].
    pub minimum_far_field_distance_m: f64,
}

/// Gas-state-dependent modal radiation transfer at one observer.
#[derive(Clone, Debug)]
pub struct ModalAcousticRadiation {
    /// Structural artifact consumed by this solve.
    pub structural_basis_identity: ContentHash,
    /// Gas model/species identity supplied by the caller.
    pub gas_model_identity: ContentHash,
    /// Temperature [K].
    pub temperature_k: f64,
    /// Absolute pressure [Pa].
    pub ambient_pressure_pa: f64,
    /// Derived acoustic density [kg/m3].
    pub density_kg_m3: f64,
    /// Derived sound speed [m/s].
    pub sound_speed_m_s: f64,
    /// Body-frame observer location [m].
    pub observer: AcousticObserver,
    /// One transfer value for every retained structural mode.
    pub modes: Vec<AcousticModeRadiation>,
    /// Identity binding the structural basis, gas state/model, observer, and
    /// computed SI transfer values.
    pub identity: ContentHash,
}

/// Frequency-dependent material loss values evaluated on one exact structural
/// basis. This is the neutral handoff accepted from any constitutive damping
/// producer; material names never enter modal time integration.
#[derive(Clone, Debug)]
pub struct ModalLossSpectrum {
    /// Structural basis whose frequencies were evaluated.
    pub structural_basis_identity: ContentHash,
    /// Material state to which the damping model applies.
    pub material_state_identity: ContentHash,
    /// Identity of the constitutive damping model and its parameter evidence.
    pub damping_model_identity: ContentHash,
    /// Loss factor `eta(omega_k)` for every retained mode.
    pub loss_factors: Vec<f64>,
}

/// One physical contact-force transition and its SI pressure observation.
#[derive(Clone, Debug)]
pub struct PhysicalModalPressureFrame {
    /// Actual point-force projection used to drive the retained modes.
    pub force_projection: PointForceProjection,
    /// Exact-ZOH modal transition and physical pressure in pascals.
    pub acoustic: ModalAcousticFrame,
}

/// Discrete semantics used to map mechanics interval controls into body-frame
/// structural forcing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhysicalContactForceSampling {
    /// The authoritative interval-mean force is held over its exact mechanics
    /// interval. Its spatial point and body transform use the closing contact
    /// endpoint when available, otherwise the opening endpoint. Audio samples
    /// are split at every mechanics boundary.
    IntervalMeanAtClosingElseOpeningEndpointZohV1,
}

/// Coordinate frame and exact location of a physical pressure observer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PhysicalPressureObserver {
    /// A legacy observer fixed in the undeformed specimen body frame.
    BodyFixed(AcousticObserver),
    /// A stationary or externally animated observer in the inertial frame.
    WorldFixed(AcousticWorldObserver),
}

/// Unmastered physical pressure signal at one observer.
#[derive(Clone, Debug)]
pub struct PhysicalPressureSignal {
    /// First sample-interval boundary on the source trajectory clock [s].
    pub start_time_s: f64,
    /// Fixed output sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Pressure at each closing sample boundary [Pa].
    pub pressure_pa: Vec<f64>,
    /// Largest absolute pressure in the signal [Pa].
    pub peak_abs_pressure_pa: f64,
    /// Exact mechanics-to-audio sampling convention.
    pub contact_force_sampling: PhysicalContactForceSampling,
    /// Exact coordinate frame and location at which pressure was evaluated.
    pub observer: PhysicalPressureObserver,
    /// Structural basis consumed by synthesis.
    pub structural_basis_identity: ContentHash,
    /// Acoustic observer transfer consumed by synthesis.
    pub radiation_identity: ContentHash,
    /// Constitutive damping model consumed by synthesis.
    pub damping_model_identity: ContentHash,
    /// Identity binding source controls, all physical models, the sampling
    /// convention, and the resulting SI pressure samples.
    pub identity: ContentHash,
}

/// Integrated physical runtime: one structural basis, state-dependent modal
/// damping, and one BEM-derived observer transfer.
pub struct PhysicalModalAudioModel<'basis> {
    basis: &'basis StructuralModalBasis,
    runtime: ModalAcousticTimeModel,
    sample_rate_hz: u32,
    static_observer: Option<AcousticObserver>,
    /// Identity of the acoustic radiation artifact.
    pub radiation_identity: ContentHash,
    /// Identity of the damping model and its evidence.
    pub damping_model_identity: ContentHash,
}

/// Typed refusal from structural-mode construction or force projection.
#[derive(Debug)]
pub enum StructuralModalBasisError {
    /// A scalar or count in the request is outside its physical domain.
    InvalidRequest {
        /// Failed invariant.
        what: &'static str,
    },
    /// The resolved profile is real geometry, but its conforming volume
    /// producer has not yet been implemented.
    UnsupportedProfile {
        /// Exact missing producer.
        what: &'static str,
    },
    /// The body-fitted volume mesh refused.
    Mesh(RoundedCylinderMeshError),
    /// The physical finite-element assembly refused.
    Elastic(TetElasticError),
    /// The certified modal solve refused.
    Modal(ModalError),
    /// The generic boundary-panel carrier refused.
    BemSurface(BemError),
    /// The exterior Helmholtz solve refused.
    Acoustic(HelmholtzError),
    /// The state-dependent viscoelastic model refused.
    Viscoelastic(ViscoError),
    /// The exact time-domain modal runtime refused.
    ModalTime(ModalAcousticTimeError),
    /// Pose reconstruction for a world-fixed observer refused.
    Timeline(TimelineResamplingError),
    /// The requested band contains no elastic modes.
    NoModesInBand,
    /// The certified count exceeds the caller's retained-mode budget.
    ModeBudgetExceeded {
        /// Certified in-band count.
        requested: usize,
        /// Caller-declared cap.
        maximum: usize,
    },
    /// A returned interval is not entirely within the stable positive branch.
    NonPositiveCertifiedMode {
        /// Zero-based retained mode index.
        mode: usize,
    },
    /// A contact point lies farther from the discrete boundary than allowed.
    ContactOutsideTolerance {
        /// Measured closest distance [m].
        distance_m: f64,
        /// Caller-declared tolerance [m].
        tolerance_m: f64,
    },
    /// The observer is not far enough for the selected far-field evaluator.
    ObserverOutsideFarField {
        /// Observer radius [m].
        distance_m: f64,
        /// Frequency-dependent minimum [m].
        minimum_m: f64,
        /// Structural mode at which the gate failed.
        mode: usize,
    },
    /// A nominally passive exterior solve returned negative outgoing power.
    NegativeRadiatedPower {
        /// Structural mode at which the gate failed.
        mode: usize,
        /// Returned power coefficient.
        power: f64,
    },
    /// A spherical-harmonic table omitted more far-field power than the
    /// caller admitted.
    DirectivityTruncation {
        /// Structural mode at which the gate failed.
        mode: usize,
        /// Quadrature-estimated captured fraction.
        captured_fraction: f64,
        /// Caller-required minimum fraction.
        minimum_fraction: f64,
    },
    /// Two independently produced physical artifacts do not share an exact
    /// structural basis or material state.
    IdentityMismatch {
        /// Failed identity relationship.
        what: &'static str,
    },
    /// A nonzero mechanics force has no retained application point.
    MissingContactLocation {
        /// Source sample closing the affected interval.
        source_sample: usize,
    },
    /// A contact-active mechanics interval has no authoritative force.
    MissingContactForce {
        /// Source sample closing the affected interval.
        source_sample: usize,
    },
    /// Source intervals cannot form one contiguous fixed-rate signal.
    ControlTimeline {
        /// Failed timeline invariant.
        what: &'static str,
    },
    /// Output pressure allocation refused.
    PressureCapacity {
        /// Requested sample count.
        requested: usize,
    },
    /// Explicit cancellation at a bounded audio-frame checkpoint.
    Cancelled,
}

impl core::fmt::Display for StructuralModalBasisError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidRequest { what } => {
                write!(formatter, "FS-EULER-STRUCTURAL-MODE-REQUEST: {what}")
            }
            Self::UnsupportedProfile { what } => {
                write!(formatter, "FS-EULER-STRUCTURAL-MODE-PROFILE: {what}")
            }
            Self::Mesh(source) => write!(formatter, "structural volume mesh refused: {source}"),
            Self::Elastic(source) => write!(formatter, "structural assembly refused: {source}"),
            Self::Modal(source) => write!(formatter, "structural modal solve refused: {source}"),
            Self::BemSurface(source) => {
                write!(formatter, "structural acoustic surface refused: {source}")
            }
            Self::Acoustic(source) => {
                write!(formatter, "structural acoustic radiation refused: {source}")
            }
            Self::Viscoelastic(source) => {
                write!(formatter, "structural damping refused: {source}")
            }
            Self::ModalTime(source) => {
                write!(
                    formatter,
                    "structural modal time integration refused: {source}"
                )
            }
            Self::Timeline(source) => {
                write!(
                    formatter,
                    "physical audio pose reconstruction refused: {source}"
                )
            }
            Self::NoModesInBand => write!(formatter, "FS-EULER-STRUCTURAL-MODE-EMPTY-BAND"),
            Self::ModeBudgetExceeded { requested, maximum } => write!(
                formatter,
                "FS-EULER-STRUCTURAL-MODE-BUDGET: certified {requested} modes exceeds {maximum}"
            ),
            Self::NonPositiveCertifiedMode { mode } => write!(
                formatter,
                "FS-EULER-STRUCTURAL-MODE-NONPOSITIVE: mode {mode} interval reaches zero"
            ),
            Self::ContactOutsideTolerance {
                distance_m,
                tolerance_m,
            } => write!(
                formatter,
                "FS-EULER-STRUCTURAL-CONTACT-DISTANCE: {distance_m:.6e} m exceeds {tolerance_m:.6e} m"
            ),
            Self::ObserverOutsideFarField {
                distance_m,
                minimum_m,
                mode,
            } => write!(
                formatter,
                "FS-EULER-ACOUSTIC-FAR-FIELD: observer {distance_m:.6e} m is below {minimum_m:.6e} m for mode {mode}"
            ),
            Self::NegativeRadiatedPower { mode, power } => write!(
                formatter,
                "FS-EULER-ACOUSTIC-NONPASSIVE: mode {mode} returned {power:.6e} W per squared modal velocity"
            ),
            Self::DirectivityTruncation {
                mode,
                captured_fraction,
                minimum_fraction,
            } => write!(
                formatter,
                "FS-EULER-ACOUSTIC-DIRECTIVITY-TRUNCATION: mode {mode} captured {captured_fraction:.6e}, below required {minimum_fraction:.6e}"
            ),
            Self::IdentityMismatch { what } => {
                write!(formatter, "FS-EULER-STRUCTURAL-IDENTITY: {what}")
            }
            Self::MissingContactLocation { source_sample } => write!(
                formatter,
                "FS-EULER-PHYSICAL-AUDIO-CONTACT-LOCATION: nonzero interval force at source sample {source_sample} has no retained contact point"
            ),
            Self::MissingContactForce { source_sample } => write!(
                formatter,
                "FS-EULER-PHYSICAL-AUDIO-CONTACT-FORCE: contact-active source sample {source_sample} has no authoritative interval force"
            ),
            Self::ControlTimeline { what } => {
                write!(formatter, "FS-EULER-PHYSICAL-AUDIO-TIMELINE: {what}")
            }
            Self::PressureCapacity { requested } => write!(
                formatter,
                "FS-EULER-PHYSICAL-AUDIO-CAPACITY: allocation of {requested} pressure samples refused"
            ),
            Self::Cancelled => formatter.write_str("FS-EULER-PHYSICAL-AUDIO-CANCELLED"),
        }
    }
}

impl std::error::Error for StructuralModalBasisError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Mesh(source) => Some(source),
            Self::Elastic(source) => Some(source),
            Self::Modal(source) => Some(source),
            Self::BemSurface(source) => Some(source),
            Self::Acoustic(source) => Some(source),
            Self::Viscoelastic(source) => Some(source),
            Self::ModalTime(source) => Some(source),
            Self::Timeline(source) => Some(source),
            _ => None,
        }
    }
}

impl From<RoundedCylinderMeshError> for StructuralModalBasisError {
    fn from(source: RoundedCylinderMeshError) -> Self {
        Self::Mesh(source)
    }
}

impl From<TetElasticError> for StructuralModalBasisError {
    fn from(source: TetElasticError) -> Self {
        Self::Elastic(source)
    }
}

impl From<ModalError> for StructuralModalBasisError {
    fn from(source: ModalError) -> Self {
        Self::Modal(source)
    }
}

impl From<BemError> for StructuralModalBasisError {
    fn from(source: BemError) -> Self {
        Self::BemSurface(source)
    }
}

impl From<HelmholtzError> for StructuralModalBasisError {
    fn from(source: HelmholtzError) -> Self {
        Self::Acoustic(source)
    }
}

impl From<ViscoError> for StructuralModalBasisError {
    fn from(source: ViscoError) -> Self {
        Self::Viscoelastic(source)
    }
}

impl From<ModalAcousticTimeError> for StructuralModalBasisError {
    fn from(source: ModalAcousticTimeError) -> Self {
        Self::ModalTime(source)
    }
}

impl From<TimelineResamplingError> for StructuralModalBasisError {
    fn from(source: TimelineResamplingError) -> Self {
        Self::Timeline(source)
    }
}

/// Assemble a body-fitted structural basis directly from a resolved specimen.
///
/// # Errors
/// Refuses invalid requests, unsupported volume-profile families, mesh or
/// elasticity failures, unresolved spectrum slices, empty bands, or mode caps.
pub fn build_structural_modal_basis(
    request: &StructuralModeRequest<'_>,
    cx: &Cx<'_>,
) -> Result<StructuralModalBasis, StructuralModalBasisError> {
    validate_request(request)?;
    let (outer_radius_m, thickness_m, fillet_radius_m) =
        rounded_cylinder_dimensions(request.specimen)?;
    let mesh_spec = RoundedCylinderMeshSpec {
        outer_radius_m,
        thickness_m,
        fillet_radius_m,
        core_radial_segments: request.mesh.core_radial_segments,
        fillet_radial_segments: request.mesh.fillet_radial_segments,
        azimuthal_segments: request.mesh.azimuthal_segments,
        axial_segments: request.mesh.axial_segments,
        maximum_vertices: request.mesh.maximum_vertices,
        maximum_tetrahedra: request.mesh.maximum_tetrahedra,
    };
    let mesh = rounded_cylinder_tet_mesh(mesh_spec, cx)?;
    let assembly = TetLinearElasticProblem {
        nodes_m: &mesh.nodes_m,
        tetrahedra: &mesh.tetrahedra,
        materials: TetMaterialField::Uniform(&request.specimen.elastic_material),
        fixed_dofs: &[],
        budget: request.assembly,
    }
    .assemble(cx)?;

    let angular_min = core::f64::consts::TAU * request.minimum_frequency_hz;
    let angular_max = core::f64::consts::TAU * request.maximum_frequency_hz;
    let eigenvalue_window_s2 = (angular_min * angular_min, angular_max * angular_max);
    let report = slice_window(
        &assembly.stiffness,
        &assembly.mass,
        eigenvalue_window_s2,
        &request.slice,
    )?;
    if report.expected == 0 {
        return Err(StructuralModalBasisError::NoModesInBand);
    }
    if report.expected > request.maximum_modes {
        return Err(StructuralModalBasisError::ModeBudgetExceeded {
            requested: report.expected,
            maximum: request.maximum_modes,
        });
    }

    let mut modes = Vec::with_capacity(report.modes.len());
    for (mode_index, pair) in report.modes.iter().enumerate() {
        if !(pair.lambda > 0.0 && pair.interval.0 > 0.0 && pair.interval.1.is_finite()) {
            return Err(StructuralModalBasisError::NonPositiveCertifiedMode { mode: mode_index });
        }
        let mut nodal_shape = vec![[0.0; 3]; mesh.nodes_m.len()];
        for (reduced_dof, &full_dof) in assembly.free_dofs.iter().enumerate() {
            nodal_shape[full_dof / 3][full_dof % 3] = pair.phi[reduced_dof];
        }
        let panel_normal_shape = mesh
            .boundary
            .triangles
            .iter()
            .zip(&mesh.boundary.normals)
            .map(|(triangle, normal)| {
                let displacement = [
                    (nodal_shape[triangle[0]][0]
                        + nodal_shape[triangle[1]][0]
                        + nodal_shape[triangle[2]][0])
                        / 3.0,
                    (nodal_shape[triangle[0]][1]
                        + nodal_shape[triangle[1]][1]
                        + nodal_shape[triangle[2]][1])
                        / 3.0,
                    (nodal_shape[triangle[0]][2]
                        + nodal_shape[triangle[1]][2]
                        + nodal_shape[triangle[2]][2])
                        / 3.0,
                ];
                dot(displacement, *normal)
            })
            .collect();
        let omega = pair.lambda.sqrt();
        modes.push(StructuralMode {
            eigenvalue_s2: pair.lambda,
            angular_frequency_rad_s: omega,
            frequency_hz: omega / core::f64::consts::TAU,
            eigenvalue_interval_s2: pair.interval,
            eigenvalue_residual_s2: pair.residual,
            nodal_shape_per_sqrt_kg: nodal_shape,
            panel_normal_shape_per_sqrt_kg: panel_normal_shape,
        });
    }

    let identity = basis_identity(request, mesh_spec, &mesh, &assembly, &modes);
    Ok(StructuralModalBasis {
        schema_version: STRUCTURAL_MODAL_BASIS_SCHEMA_VERSION,
        specimen_identity: request.specimen.identity,
        profile_identity: request.specimen.profile.content_identities().profile,
        material_state_identity: request.specimen.material_state_identity,
        mesh,
        assembly,
        modes,
        eigenvalue_window_s2,
        certified_mode_count: report.expected,
        slice_stats: report.stats,
        identity,
    })
}

/// Evaluate a certified generalized-Maxwell material model at every retained
/// structural frequency.
///
/// This adapter intentionally accepts an explicit model identity. The
/// parameter/evidence author, not the material's display name, owns that
/// identity and the certified frequency band.
///
/// # Errors
/// Refuses a foreign specimen/material binding, a zero model identity, or any
/// modal frequency outside the lowered model's certified band.
pub fn modal_loss_spectrum_from_prony(
    basis: &StructuralModalBasis,
    specimen: &ResolvedElasticDiscProfile,
    model: &LoweredModel,
    damping_model_identity: ContentHash,
) -> Result<ModalLossSpectrum, StructuralModalBasisError> {
    if basis.specimen_identity != specimen.identity
        || basis.material_state_identity != specimen.material_state_identity
    {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "damping specimen does not match the structural basis",
        });
    }
    if damping_model_identity == ContentHash([0; 32]) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "damping_model_identity must not be zero",
        });
    }
    let mut loss_factors = Vec::with_capacity(basis.modes.len());
    for mode in &basis.modes {
        let loss = model.loss_factor_checked(mode.angular_frequency_rad_s)?;
        if !(loss >= 0.0 && loss.is_finite()) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "constitutive modal loss factor must be finite and non-negative",
            });
        }
        loss_factors.push(loss);
    }
    Ok(ModalLossSpectrum {
        structural_basis_identity: basis.identity,
        material_state_identity: basis.material_state_identity,
        damping_model_identity,
        loss_factors,
    })
}

/// Evaluate a caller-identified Rayleigh damping law on one structural basis.
///
/// This is the inexpensive constitutive rung for materials whose damping data
/// are available as mass- and stiffness-proportional coefficients. The
/// coefficients are numerical material/process state, not a material-name
/// preset; changing either coefficient changes the retained loss spectrum and
/// all downstream pressure samples.
///
/// # Errors
/// Refuses a foreign specimen/material binding, a zero damping identity, or a
/// non-finite/negative modal loss factor.
pub fn modal_loss_spectrum_from_rayleigh(
    basis: &StructuralModalBasis,
    specimen: &ResolvedElasticDiscProfile,
    model: RayleighDamping,
    damping_model_identity: ContentHash,
) -> Result<ModalLossSpectrum, StructuralModalBasisError> {
    if basis.specimen_identity != specimen.identity
        || basis.material_state_identity != specimen.material_state_identity
    {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "damping specimen does not match the structural basis",
        });
    }
    if damping_model_identity == ContentHash([0; 32]) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "damping_model_identity must not be zero",
        });
    }
    let mut loss_factors = Vec::with_capacity(basis.modes.len());
    for mode in &basis.modes {
        let loss = 2.0 * model.zeta_at(mode.angular_frequency_rad_s);
        if !(loss >= 0.0 && loss.is_finite()) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "Rayleigh modal loss factor must be finite and non-negative",
            });
        }
        loss_factors.push(loss);
    }
    Ok(ModalLossSpectrum {
        structural_basis_identity: basis.identity,
        material_state_identity: basis.material_state_identity,
        damping_model_identity,
        loss_factors,
    })
}

impl<'basis> PhysicalModalAudioModel<'basis> {
    /// Bind structural, damping, and radiation artifacts into one exact-ZOH
    /// physical-pressure runtime.
    ///
    /// # Errors
    /// Refuses foreign identities, wrong cardinalities, malformed loss
    /// factors, frequencies above the Nyquist guard, or invalid budgets.
    pub fn try_new(
        basis: &'basis StructuralModalBasis,
        loss: &ModalLossSpectrum,
        radiation: &ModalAcousticRadiation,
        sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, StructuralModalBasisError> {
        if loss.structural_basis_identity != basis.identity
            || radiation.structural_basis_identity != basis.identity
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "damping or acoustic artifact does not match the structural basis",
            });
        }
        if loss.material_state_identity != basis.material_state_identity {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "damping material state does not match structural assembly material state",
            });
        }
        if loss.damping_model_identity == ContentHash([0; 32]) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "damping_model_identity must not be zero",
            });
        }
        if radiation.identity == ContentHash([0; 32])
            || radiation.gas_model_identity == ContentHash([0; 32])
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "acoustic radiation and gas model identities must not be zero",
            });
        }
        if loss.loss_factors.len() != basis.modes.len()
            || radiation.modes.len() != basis.modes.len()
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "loss, radiation, and structural mode counts must agree",
            });
        }
        let mut modes = Vec::with_capacity(basis.modes.len());
        for (mode_index, ((structural, loss_factor), acoustic)) in basis
            .modes
            .iter()
            .zip(&loss.loss_factors)
            .zip(&radiation.modes)
            .enumerate()
        {
            if acoustic.structural_mode != mode_index
                || acoustic.angular_frequency_rad_s.to_bits()
                    != structural.angular_frequency_rad_s.to_bits()
                || !(loss_factor.is_finite() && *loss_factor >= 0.0)
                || acoustic.formulation == HelmholtzFormulation::BurtonMillerWrongAlphaSign
                || !(acoustic.radiated_power_per_modal_velocity_squared >= 0.0
                    && acoustic
                        .radiated_power_per_modal_velocity_squared
                        .is_finite()
                    && acoustic.panels_per_wavelength > 0.0
                    && acoustic.panels_per_wavelength.is_finite()
                    && acoustic.condition_lower_bound >= 1.0
                    && acoustic.condition_lower_bound.is_finite()
                    && acoustic.minimum_far_field_distance_m > 0.0
                    && acoustic.minimum_far_field_distance_m.is_finite())
            {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "per-mode structural, damping, and radiation rows are misaligned",
                });
            }
            modes.push(ModalAcousticMode {
                angular_frequency_rad_s: structural.angular_frequency_rad_s,
                damping_ratio: loss_factor_to_zeta(*loss_factor),
                pressure_per_modal_velocity: acoustic.observer_pressure_per_modal_velocity,
            });
        }
        let runtime = ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?;
        Ok(Self {
            basis,
            runtime,
            sample_rate_hz,
            static_observer: Some(radiation.observer),
            radiation_identity: radiation.identity,
            damping_model_identity: loss.damping_model_identity,
        })
    }

    /// Bind structural modes and constitutive damping to an
    /// observer-independent directivity artifact.
    ///
    /// The modal oscillator is identical to [`Self::try_new`], but pressure is
    /// deliberately initialized with zero static transfers. Callers must use
    /// [`Self::synthesize_control_stream_world_observers`] so every sample is
    /// observed through an explicit current pose and world microphone.
    ///
    /// # Errors
    /// Refuses foreign identities, wrong cardinalities, malformed loss or
    /// directivity rows, frequencies above Nyquist, or invalid budgets.
    pub fn try_new_with_directivity(
        basis: &'basis StructuralModalBasis,
        loss: &ModalLossSpectrum,
        directivity: &ModalAcousticDirectivity,
        sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, StructuralModalBasisError> {
        validate_modal_acoustic_directivity(directivity)?;
        if loss.structural_basis_identity != basis.identity
            || directivity.structural_basis_identity != basis.identity
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "damping or acoustic directivity artifact does not match the structural basis",
            });
        }
        if loss.material_state_identity != basis.material_state_identity {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "damping material state does not match structural assembly material state",
            });
        }
        if loss.damping_model_identity == ContentHash([0; 32])
            || directivity.identity == ContentHash([0; 32])
            || directivity.gas_model_identity == ContentHash([0; 32])
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "damping, directivity, and gas model identities must not be zero",
            });
        }
        if loss.loss_factors.len() != basis.modes.len()
            || directivity.modes.len() != basis.modes.len()
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "loss, directivity, and structural mode counts must agree",
            });
        }
        let mut modes = Vec::with_capacity(basis.modes.len());
        for (mode_index, ((structural, loss_factor), acoustic)) in basis
            .modes
            .iter()
            .zip(&loss.loss_factors)
            .zip(&directivity.modes)
            .enumerate()
        {
            if acoustic.structural_mode != mode_index
                || acoustic.angular_frequency_rad_s.to_bits()
                    != structural.angular_frequency_rad_s.to_bits()
                || !(loss_factor.is_finite() && *loss_factor >= 0.0)
            {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "per-mode structural, damping, and directivity rows are misaligned",
                });
            }
            modes.push(ModalAcousticMode {
                angular_frequency_rad_s: structural.angular_frequency_rad_s,
                damping_ratio: loss_factor_to_zeta(*loss_factor),
                pressure_per_modal_velocity: C64::ZERO,
            });
        }
        let runtime = ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?;
        Ok(Self {
            basis,
            runtime,
            sample_rate_hz,
            static_observer: None,
            radiation_identity: directivity.identity,
            damping_model_identity: loss.damping_model_identity,
        })
    }

    /// Current physical sample-boundary modal states.
    #[must_use]
    pub fn states(&self) -> &[fs_couple::modal_acoustic_time::ModalAcousticState] {
        self.runtime.states()
    }

    /// Project one body-frame point force and advance the physical pressure
    /// runtime by one audio sample.
    ///
    /// # Errors
    /// Refuses an off-boundary contact, invalid force, or a transactional
    /// modal-time failure. A refusal leaves all modal states unchanged.
    pub fn step_point_force(
        &mut self,
        point_body_m: [f64; 3],
        force_body_n: [f64; 3],
        maximum_distance_m: f64,
    ) -> Result<PhysicalModalPressureFrame, StructuralModalBasisError> {
        let force_projection =
            self.basis
                .project_point_force(point_body_m, force_body_n, maximum_distance_m)?;
        let acoustic = self
            .runtime
            .step(&force_projection.modal_force_n_per_sqrt_kg)?;
        Ok(PhysicalModalPressureFrame {
            force_projection,
            acoustic,
        })
    }

    /// Produce one fixed-rate, unmastered SI-pressure signal from an admitted
    /// Euler control stream.
    ///
    /// Mechanics force changes are never rounded onto the audio grid: an
    /// output sample crossing a mechanics boundary is integrated as multiple
    /// exact-ZOH modal substeps. When a full contact wrench is unavailable,
    /// an explicitly admitted interval-mean normal reaction is reconstructed
    /// along the retained contact normal. No impulse is invented from a
    /// timing-only contact event.
    ///
    /// # Errors
    /// Refuses discontinuous clocks, nonintegral output duration, missing
    /// force/location authority, projection/runtime failures, capacity, or
    /// cancellation. A runtime refusal may have committed earlier complete
    /// audio frames; callers seeking all-or-nothing publication must discard
    /// this model instance together with the returned error.
    pub fn synthesize_control_stream(
        &mut self,
        controls: &EulerControlStream<'_>,
        maximum_contact_distance_m: f64,
        cx: &Cx<'_>,
    ) -> Result<PhysicalPressureSignal, StructuralModalBasisError> {
        let observer = self
            .static_observer
            .ok_or(StructuralModalBasisError::InvalidRequest {
                what: "world-observer runtime requires world-observer synthesis",
            })?;
        let intervals = controls.audio();
        if controls.source().metadata().specimen_profile_identity != self.basis.profile_identity {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "control trajectory specimen profile does not match structural basis profile",
            });
        }
        let first = intervals
            .first()
            .ok_or(StructuralModalBasisError::ControlTimeline {
                what: "control stream has no positive-duration audio intervals",
            })?;
        let last = intervals
            .last()
            .expect("nonempty interval slice has a last item");
        for pair in intervals.windows(2) {
            if pair[0].end_time_s.to_bits() != pair[1].start_time_s.to_bits() {
                return Err(StructuralModalBasisError::ControlTimeline {
                    what: "mechanics audio intervals are not exactly contiguous",
                });
            }
        }
        let duration_s = last.end_time_s - first.start_time_s;
        let exact_frames = duration_s * f64::from(self.sample_rate_hz);
        let rounded_frames = exact_frames.round();
        let frame_tolerance = 128.0 * f64::EPSILON * exact_frames.abs().max(1.0);
        if !(duration_s > 0.0
            && duration_s.is_finite()
            && exact_frames.is_finite()
            && rounded_frames >= 1.0
            && (exact_frames - rounded_frames).abs() <= frame_tolerance
            && rounded_frames <= usize::MAX as f64)
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "control horizon is not a positive integral number of audio samples",
            });
        }
        let frame_count = rounded_frames as usize;
        let mut pressure_pa = Vec::new();
        pressure_pa.try_reserve_exact(frame_count).map_err(|_| {
            StructuralModalBasisError::PressureCapacity {
                requested: frame_count,
            }
        })?;
        let mut modal_forces = Vec::with_capacity(intervals.len());
        for interval in intervals {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
            modal_forces.push(self.modal_force_for_control_interval(
                controls,
                interval,
                maximum_contact_distance_m,
            )?);
        }

        let sample_period_s = self.runtime.sample_period_s();
        let mut interval_index = 0usize;
        let mut peak_abs_pressure_pa = 0.0_f64;
        let mut sample_start = first.start_time_s;
        for frame in 0..frame_count {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            // Anchor both horizon endpoints to the authoritative mechanics
            // clock.  Repeated fixed-rate arithmetic is used only for
            // interior boundaries, so the final sample cannot overshoot the
            // source horizon by a floating-point ulp.
            let sample_end = if frame + 1 == frame_count {
                last.end_time_s
            } else {
                first.start_time_s + (frame + 1) as f64 * sample_period_s
            };
            let mut time = sample_start;
            let mut final_pressure = 0.0;
            while time < sample_end {
                while interval_index + 1 < intervals.len()
                    && time >= intervals[interval_index].end_time_s
                {
                    interval_index += 1;
                }
                let interval = &intervals[interval_index];
                let segment_end = sample_end.min(interval.end_time_s);
                let segment_duration = segment_end - time;
                if !(segment_duration > 0.0 && segment_duration.is_finite()) {
                    return Err(StructuralModalBasisError::ControlTimeline {
                        what: "audio sample subdivision made no forward progress",
                    });
                }
                final_pressure = self
                    .runtime
                    .step_duration(&modal_forces[interval_index], segment_duration)?
                    .observer_pressure_pa;
                time = segment_end;
            }
            peak_abs_pressure_pa = peak_abs_pressure_pa.max(final_pressure.abs());
            pressure_pa.push(final_pressure);
            sample_start = sample_end;
        }
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let observer = PhysicalPressureObserver::BodyFixed(observer);
        let identity = physical_pressure_signal_identity(
            self,
            controls,
            first.start_time_s,
            observer,
            &pressure_pa,
        );
        Ok(PhysicalPressureSignal {
            start_time_s: first.start_time_s,
            sample_rate_hz: self.sample_rate_hz,
            pressure_pa,
            peak_abs_pressure_pa,
            contact_force_sampling:
                PhysicalContactForceSampling::IntervalMeanAtClosingElseOpeningEndpointZohV1,
            observer,
            structural_basis_identity: self.basis.identity,
            radiation_identity: self.radiation_identity,
            damping_model_identity: self.damping_model_identity,
            identity,
        })
    }

    /// Produce simultaneous physical-pressure signals at world-fixed
    /// observers from one admitted control stream.
    ///
    /// Mechanics forcing advances one shared modal state exactly once. At
    /// each closing audio boundary the source pose is reconstructed from the
    /// authoritative trajectory, each world line of sight is rotated into the
    /// body-frame BEM directivity field, and each observer gets its own SI
    /// pressure sample. This preserves inter-channel phase and avoids running
    /// two independently drifting oscillator copies.
    ///
    /// Propagation is the narrow-band Helmholtz phase of each retained mode at
    /// the current sample-boundary pose. It does not claim a broadband moving-
    /// boundary retarded-time solution or room reflections.
    ///
    /// # Errors
    /// Refuses the static-observer constructor, foreign directivity, empty or
    /// excessive observer sets, invalid timelines/forces/poses, far-field
    /// violations, capacity, cancellation, or modal-time failure.
    pub fn synthesize_control_stream_world_observers(
        &mut self,
        controls: &EulerControlStream<'_>,
        directivity: &ModalAcousticDirectivity,
        observers: &[AcousticWorldObserver],
        maximum_contact_distance_m: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<PhysicalPressureSignal>, StructuralModalBasisError> {
        if self.static_observer.is_some() {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "static-observer runtime cannot synthesize world observers",
            });
        }
        if directivity.identity != self.radiation_identity
            || directivity.structural_basis_identity != self.basis.identity
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "world-observer directivity does not match the admitted audio runtime",
            });
        }
        if observers.is_empty() || observers.len() > MAX_PHYSICAL_PRESSURE_OBSERVERS {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "world observer count must be in 1..=64",
            });
        }
        let intervals = controls.audio();
        if controls.source().metadata().specimen_profile_identity != self.basis.profile_identity {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "control trajectory specimen profile does not match structural basis profile",
            });
        }
        let first = intervals
            .first()
            .ok_or(StructuralModalBasisError::ControlTimeline {
                what: "control stream has no positive-duration audio intervals",
            })?;
        let last = intervals
            .last()
            .expect("nonempty interval slice has a last item");
        for pair in intervals.windows(2) {
            if pair[0].end_time_s.to_bits() != pair[1].start_time_s.to_bits() {
                return Err(StructuralModalBasisError::ControlTimeline {
                    what: "mechanics audio intervals are not exactly contiguous",
                });
            }
        }
        let duration_s = last.end_time_s - first.start_time_s;
        let exact_frames = duration_s * f64::from(self.sample_rate_hz);
        let rounded_frames = exact_frames.round();
        let frame_tolerance = 128.0 * f64::EPSILON * exact_frames.abs().max(1.0);
        if !(duration_s > 0.0
            && duration_s.is_finite()
            && exact_frames.is_finite()
            && rounded_frames >= 1.0
            && (exact_frames - rounded_frames).abs() <= frame_tolerance
            && rounded_frames <= usize::MAX as f64)
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "control horizon is not a positive integral number of audio samples",
            });
        }
        let frame_count = rounded_frames as usize;
        let mut pressure_channels = vec![Vec::new(); observers.len()];
        for pressure in &mut pressure_channels {
            pressure.try_reserve_exact(frame_count).map_err(|_| {
                StructuralModalBasisError::PressureCapacity {
                    requested: frame_count,
                }
            })?;
        }
        let mut modal_forces = Vec::with_capacity(intervals.len());
        for interval in intervals {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
            modal_forces.push(self.modal_force_for_control_interval(
                controls,
                interval,
                maximum_contact_distance_m,
            )?);
        }

        let timeline = TimelineResampler::new(controls.source());
        let mut transfer_scratch = vec![Vec::new(); observers.len()];
        let mut peaks = vec![0.0_f64; observers.len()];
        let sample_period_s = self.runtime.sample_period_s();
        let mut interval_index = 0usize;
        let mut sample_start = first.start_time_s;
        for frame in 0..frame_count {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let sample_end = if frame + 1 == frame_count {
                last.end_time_s
            } else {
                first.start_time_s + (frame + 1) as f64 * sample_period_s
            };
            let mut time = sample_start;
            while time < sample_end {
                while interval_index + 1 < intervals.len()
                    && time >= intervals[interval_index].end_time_s
                {
                    interval_index += 1;
                }
                let segment_end = sample_end.min(intervals[interval_index].end_time_s);
                let segment_duration = segment_end - time;
                if !(segment_duration > 0.0 && segment_duration.is_finite()) {
                    return Err(StructuralModalBasisError::ControlTimeline {
                        what: "audio sample subdivision made no forward progress",
                    });
                }
                self.runtime
                    .step_duration(&modal_forces[interval_index], segment_duration)?;
                time = segment_end;
            }
            let pose = timeline
                .sample(sample_end, EventEvaluationSide::RightLimit)?
                .state
                .pose();
            for (observer_index, observer) in observers.iter().copied().enumerate() {
                directivity.write_observer_transfers_at_pose(
                    pose,
                    observer,
                    &mut transfer_scratch[observer_index],
                )?;
                let pressure = self
                    .runtime
                    .observer_pressure_with_transfers(&transfer_scratch[observer_index])?;
                peaks[observer_index] = peaks[observer_index].max(pressure.abs());
                pressure_channels[observer_index].push(pressure);
            }
            sample_start = sample_end;
        }
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;

        Ok(observers
            .iter()
            .copied()
            .zip(pressure_channels)
            .zip(peaks)
            .map(|((observer, pressure_pa), peak_abs_pressure_pa)| {
                let observer = PhysicalPressureObserver::WorldFixed(observer);
                let identity = physical_pressure_signal_identity(
                    self,
                    controls,
                    first.start_time_s,
                    observer,
                    &pressure_pa,
                );
                PhysicalPressureSignal {
                    start_time_s: first.start_time_s,
                    sample_rate_hz: self.sample_rate_hz,
                    pressure_pa,
                    peak_abs_pressure_pa,
                    contact_force_sampling:
                        PhysicalContactForceSampling::IntervalMeanAtClosingElseOpeningEndpointZohV1,
                    observer,
                    structural_basis_identity: self.basis.identity,
                    radiation_identity: self.radiation_identity,
                    damping_model_identity: self.damping_model_identity,
                    identity,
                }
            })
            .collect())
    }

    fn modal_force_for_control_interval(
        &self,
        controls: &EulerControlStream<'_>,
        interval: &crate::AudioControlInterval,
        maximum_contact_distance_m: f64,
    ) -> Result<Vec<f64>, StructuralModalBasisError> {
        let visualization = controls.visualization();
        let end = visualization
            .get(interval.visual_coverage.end_visualization_index)
            .ok_or(StructuralModalBasisError::ControlTimeline {
                what: "audio interval closing visualization index is out of bounds",
            })?;
        let start = interval
            .visual_coverage
            .start_visualization_index
            .and_then(|index| visualization.get(index));
        let (point, orientation, normal_world) = if let Some(contact) = end.contact {
            (
                Some(contact.point_body_m),
                end.disc_pose.orientation(),
                Some(contact.normal_world),
            )
        } else if let Some((start, contact)) =
            start.and_then(|start| start.contact.map(|c| (start, c)))
        {
            (
                Some(contact.point_body_m),
                start.disc_pose.orientation(),
                Some(contact.normal_world),
            )
        } else {
            (None, end.disc_pose.orientation(), None)
        };
        let force_world = match interval.channels.contact {
            ChannelControl::Available(contact) => contact.mean_force_world_n,
            ChannelControl::Unavailable => {
                if let (Some(normal_force), Some(normal)) =
                    (interval.mean_base_normal_contact_force_n, normal_world)
                {
                    normal.scale(normal_force)
                } else if interval.interval_contact_active {
                    return Err(StructuralModalBasisError::MissingContactForce {
                        source_sample: interval.source_sample_index,
                    });
                } else {
                    fs_mbd::Vec3::ZERO
                }
            }
        };
        if force_world.norm_squared() == 0.0 {
            return Ok(vec![0.0; self.basis.modes.len()]);
        }
        let point = point.ok_or(StructuralModalBasisError::MissingContactLocation {
            source_sample: interval.source_sample_index,
        })?;
        let force_body = orientation.rotate_world_to_body(force_world);
        self.basis
            .project_point_force(
                [point.x, point.y, point.z],
                [force_body.x, force_body.y, force_body.z],
                maximum_contact_distance_m,
            )
            .map(|projection| projection.modal_force_n_per_sqrt_kg)
    }
}

fn physical_pressure_signal_identity(
    model: &PhysicalModalAudioModel<'_>,
    controls: &EulerControlStream<'_>,
    start_time_s: f64,
    observer: PhysicalPressureObserver,
    pressure_pa: &[f64],
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-disc.physical-pressure-signal.v1");
    hasher.update(model.basis.identity.as_bytes());
    hasher.update(model.radiation_identity.as_bytes());
    hasher.update(model.damping_model_identity.as_bytes());
    hasher.update(
        controls
            .source()
            .metadata()
            .configuration_identity
            .as_bytes(),
    );
    hasher.update(&model.sample_rate_hz.to_le_bytes());
    hasher.update(&start_time_s.to_bits().to_le_bytes());
    match observer {
        PhysicalPressureObserver::BodyFixed(observer) => {
            hasher.update(&[0]);
            for value in observer.position_m {
                hasher.update(&value.to_bits().to_le_bytes());
            }
        }
        PhysicalPressureObserver::WorldFixed(observer) => {
            hasher.update(&[1]);
            for value in observer.position_world_m {
                hasher.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    hasher.update(&[0]);
    hasher.update(
        &u64::try_from(pressure_pa.len())
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    for pressure in pressure_pa {
        hasher.update(&pressure.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

impl StructuralModalBasis {
    /// Project one physical force at a point on (or close to) the discrete
    /// boundary onto every retained mode using the closest triangle and P1
    /// barycentric interpolation.
    ///
    /// # Errors
    /// Refuses non-finite values, negative tolerance, or a point farther from
    /// the discrete boundary than `maximum_distance_m`.
    pub fn project_point_force(
        &self,
        point_m: [f64; 3],
        force_n: [f64; 3],
        maximum_distance_m: f64,
    ) -> Result<PointForceProjection, StructuralModalBasisError> {
        if point_m
            .iter()
            .chain(force_n.iter())
            .any(|value| !value.is_finite())
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "contact point and force must be finite",
            });
        }
        if !(maximum_distance_m.is_finite() && maximum_distance_m >= 0.0) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "maximum contact-to-boundary distance must be finite and non-negative",
            });
        }
        let mut best = None;
        for (boundary_triangle, triangle) in self.mesh.boundary.triangles.iter().enumerate() {
            let vertices = triangle.map(|node| self.mesh.nodes_m[node]);
            let (closest, barycentric) = closest_point_on_triangle(point_m, vertices);
            let distance_squared = norm_squared(sub(point_m, closest));
            if best.as_ref().is_none_or(
                |(_, _, _, best_distance): &(usize, [f64; 3], [f64; 3], f64)| {
                    distance_squared < *best_distance
                },
            ) {
                best = Some((boundary_triangle, closest, barycentric, distance_squared));
            }
        }
        let (boundary_triangle, closest_point_m, barycentric, distance_squared) =
            best.expect("a structural mesh always has a non-empty boundary");
        let distance_to_boundary_m = distance_squared.sqrt();
        if distance_to_boundary_m > maximum_distance_m {
            return Err(StructuralModalBasisError::ContactOutsideTolerance {
                distance_m: distance_to_boundary_m,
                tolerance_m: maximum_distance_m,
            });
        }
        let triangle = self.mesh.boundary.triangles[boundary_triangle];
        let modal_force_n_per_sqrt_kg = self
            .modes
            .iter()
            .map(|mode| {
                let mut shape = [0.0; 3];
                for corner in 0..3 {
                    let nodal = mode.nodal_shape_per_sqrt_kg[triangle[corner]];
                    for component in 0..3 {
                        shape[component] =
                            barycentric[corner].mul_add(nodal[component], shape[component]);
                    }
                }
                dot(shape, force_n)
            })
            .collect();
        Ok(PointForceProjection {
            boundary_triangle,
            closest_point_m,
            distance_to_boundary_m,
            barycentric,
            modal_force_n_per_sqrt_kg,
        })
    }

    /// Compute exterior acoustic radiation at every retained natural
    /// frequency from the same boundary-normal mode shapes used by contact.
    ///
    /// The current BEM exposes an asymptotic far-field evaluator, so this
    /// method enforces a frequency-dependent Fraunhofer distance instead of
    /// silently using it at a near-field microphone. Returned pressure is
    /// physical SI pressure per generalized modal velocity; it has no digital
    /// full-scale gain or loudness mastering folded into it.
    ///
    /// # Errors
    /// Refuses malformed or identity-free media, invalid observers, far-field
    /// violations, BEM resolution/work failures, or negative radiated power.
    pub fn modal_acoustic_radiation(
        &self,
        medium: ResolvedAcousticMedium<'_>,
        observer: AcousticObserver,
    ) -> Result<ModalAcousticRadiation, StructuralModalBasisError> {
        if medium.gas_model_identity == ContentHash([0; 32]) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "gas_model_identity must not be zero",
            });
        }
        validate_acoustic_medium(medium.gas)?;
        if observer.position_m.iter().any(|value| !value.is_finite()) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "observer position must be finite",
            });
        }
        let observer_distance = norm_squared(observer.position_m).sqrt();
        if !(observer_distance > 0.0 && observer_distance.is_finite()) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "observer must not coincide with the specimen origin",
            });
        }
        let surface = SpherePanels::new(
            self.mesh.boundary.centroids_m.clone(),
            self.mesh.boundary.normals.clone(),
            self.mesh.boundary.areas_m2.clone(),
        )?;
        let acoustic_medium = Medium {
            density: medium.gas.density,
            sound_speed: medium.gas.sound_speed,
        };
        let direction = [
            observer.position_m[0] / observer_distance,
            observer.position_m[1] / observer_distance,
            observer.position_m[2] / observer_distance,
        ];
        let diameter_m = self
            .mesh
            .nodes_m
            .iter()
            .map(|node| 2.0 * norm_squared(*node).sqrt())
            .fold(0.0_f64, f64::max);
        let mut radiations = Vec::with_capacity(self.modes.len());
        for (mode_index, mode) in self.modes.iter().enumerate() {
            let k = mode.angular_frequency_rad_s / acoustic_medium.sound_speed;
            let wavelength_m = core::f64::consts::TAU / k;
            // Fraunhofer aperture criterion plus an enclosing-body distance.
            // This is an applicability rule for the asymptotic evaluator, not
            // a claim that its truncation error is rigorously enclosed.
            let minimum_far_field_distance_m =
                (2.0 * diameter_m * diameter_m / wavelength_m).max(2.0 * diameter_m);
            if observer_distance < minimum_far_field_distance_m {
                return Err(StructuralModalBasisError::ObserverOutsideFarField {
                    distance_m: observer_distance,
                    minimum_m: minimum_far_field_distance_m,
                    mode: mode_index,
                });
            }
            let velocity: Vec<C64> = mode
                .panel_normal_shape_per_sqrt_kg
                .iter()
                .map(|value| C64::from_re(*value))
                .collect();
            let acoustic_radius_m = 0.5 * diameter_m;
            let formulation = if k * acoustic_radius_m < 0.5 {
                HelmholtzFormulation::PlainCbie
            } else {
                HelmholtzFormulation::BurtonMiller
            };
            let solution = solve_radiation(&surface, k, acoustic_medium, &velocity, formulation)?;
            if solution.radiated_power < 0.0 {
                return Err(StructuralModalBasisError::NegativeRadiatedPower {
                    mode: mode_index,
                    power: solution.radiated_power,
                });
            }
            let far = far_field(&surface, &solution, acoustic_medium, &[direction])[0];
            let phase = C64::new(
                det::cos(k * observer_distance),
                det::sin(k * observer_distance),
            );
            radiations.push(AcousticModeRadiation {
                structural_mode: mode_index,
                angular_frequency_rad_s: mode.angular_frequency_rad_s,
                formulation,
                observer_pressure_per_modal_velocity: (far * phase)
                    .scale(observer_distance.recip()),
                radiated_power_per_modal_velocity_squared: solution.radiated_power,
                panels_per_wavelength: solution.panels_per_wavelength,
                condition_lower_bound: solution.condition_lower_bound,
                minimum_far_field_distance_m,
            });
        }
        let identity = acoustic_radiation_identity(self, medium, observer, &radiations);
        Ok(ModalAcousticRadiation {
            structural_basis_identity: self.identity,
            gas_model_identity: medium.gas_model_identity,
            temperature_k: medium.gas.temperature,
            ambient_pressure_pa: medium.gas.pressure,
            density_kg_m3: medium.gas.density,
            sound_speed_m_s: medium.gas.sound_speed,
            observer,
            modes: radiations,
            identity,
        })
    }

    /// Solve observer-independent far-field directivity for every retained
    /// structural mode.
    ///
    /// The BEM velocity boundary condition comes directly from the
    /// mass-normalized elastic mode shape. The resulting body-frame spherical
    /// harmonic tables can be evaluated for any rigid pose and world-fixed
    /// microphone; neither material names nor hand-authored pan/gain curves
    /// enter the transfer.
    ///
    /// # Errors
    /// Refuses malformed media/controls, BEM work or resolution failures,
    /// negative radiated power, and directivity truncation below the explicit
    /// captured-power floor.
    pub fn modal_acoustic_directivity(
        &self,
        medium: ResolvedAcousticMedium<'_>,
        controls: AcousticDirectivityControls,
    ) -> Result<ModalAcousticDirectivity, StructuralModalBasisError> {
        if medium.gas_model_identity == ContentHash([0; 32]) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "gas_model_identity must not be zero",
            });
        }
        validate_acoustic_medium(medium.gas)?;
        if !(controls.minimum_captured_fraction > 0.0
            && controls.minimum_captured_fraction <= 1.0
            && controls.minimum_captured_fraction.is_finite())
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "minimum directivity captured fraction must be finite and in (0, 1]",
            });
        }
        let surface = SpherePanels::new(
            self.mesh.boundary.centroids_m.clone(),
            self.mesh.boundary.normals.clone(),
            self.mesh.boundary.areas_m2.clone(),
        )?;
        let acoustic_medium = Medium {
            density: medium.gas.density,
            sound_speed: medium.gas.sound_speed,
        };
        let diameter_m = self
            .mesh
            .nodes_m
            .iter()
            .map(|node| 2.0 * norm_squared(*node).sqrt())
            .fold(0.0_f64, f64::max);
        let mut modes = Vec::with_capacity(self.modes.len());
        for (mode_index, mode) in self.modes.iter().enumerate() {
            let k = mode.angular_frequency_rad_s / acoustic_medium.sound_speed;
            let wavelength_m = core::f64::consts::TAU / k;
            let minimum_far_field_distance_m =
                (2.0 * diameter_m * diameter_m / wavelength_m).max(2.0 * diameter_m);
            let velocity: Vec<C64> = mode
                .panel_normal_shape_per_sqrt_kg
                .iter()
                .map(|value| C64::from_re(*value))
                .collect();
            let formulation = if k * (0.5 * diameter_m) < 0.5 {
                HelmholtzFormulation::PlainCbie
            } else {
                HelmholtzFormulation::BurtonMiller
            };
            let solution = solve_radiation(&surface, k, acoustic_medium, &velocity, formulation)?;
            if solution.radiated_power < 0.0 {
                return Err(StructuralModalBasisError::NegativeRadiatedPower {
                    mode: mode_index,
                    power: solution.radiated_power,
                });
            }
            let directivity = directivity_sh_table(
                &surface,
                &solution,
                acoustic_medium,
                controls.maximum_spherical_harmonic_degree,
            )?;
            if solution.radiated_power > 0.0
                && directivity.captured_fraction < controls.minimum_captured_fraction
            {
                return Err(StructuralModalBasisError::DirectivityTruncation {
                    mode: mode_index,
                    captured_fraction: directivity.captured_fraction,
                    minimum_fraction: controls.minimum_captured_fraction,
                });
            }
            modes.push(AcousticModeDirectivity {
                structural_mode: mode_index,
                angular_frequency_rad_s: mode.angular_frequency_rad_s,
                wavenumber_rad_m: k,
                formulation,
                directivity,
                radiated_power_per_modal_velocity_squared: solution.radiated_power,
                panels_per_wavelength: solution.panels_per_wavelength,
                condition_lower_bound: solution.condition_lower_bound,
                minimum_far_field_distance_m,
            });
        }
        let identity = acoustic_directivity_identity(self, medium, controls, &modes);
        Ok(ModalAcousticDirectivity {
            structural_basis_identity: self.identity,
            gas_model_identity: medium.gas_model_identity,
            temperature_k: medium.gas.temperature,
            ambient_pressure_pa: medium.gas.pressure,
            density_kg_m3: medium.gas.density,
            sound_speed_m_s: medium.gas.sound_speed,
            modes,
            identity,
        })
    }
}

impl ModalAcousticDirectivity {
    /// Evaluate one physical pressure transfer per mode for a world-fixed
    /// microphone and the specimen's current rigid pose.
    ///
    /// Translation controls spherical spreading/propagation phase. Rotation
    /// maps the world line of sight into the body-frame directivity field.
    /// The observer is never implicitly attached to the moving body.
    ///
    /// # Errors
    /// Refuses non-finite/coincident observers, malformed mode tables, or a
    /// pose-observer distance below any mode's far-field applicability bound.
    pub fn observer_transfers_at_pose(
        &self,
        pose: Pose,
        observer: AcousticWorldObserver,
    ) -> Result<Vec<C64>, StructuralModalBasisError> {
        let mut transfers = Vec::with_capacity(self.modes.len());
        self.write_observer_transfers_at_pose(pose, observer, &mut transfers)?;
        Ok(transfers)
    }

    /// Evaluate into a reusable caller-owned buffer.
    ///
    /// This is equivalent to [`Self::observer_transfers_at_pose`] but avoids
    /// per-sample allocation in fixed-rate audio synthesis.
    pub fn write_observer_transfers_at_pose(
        &self,
        pose: Pose,
        observer: AcousticWorldObserver,
        transfers: &mut Vec<C64>,
    ) -> Result<(), StructuralModalBasisError> {
        validate_modal_acoustic_directivity(self)?;
        if observer
            .position_world_m
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "world observer position must be finite",
            });
        }
        let source_world = pose.position_world();
        let relative_world = Vec3::new(
            observer.position_world_m[0] - source_world.x,
            observer.position_world_m[1] - source_world.y,
            observer.position_world_m[2] - source_world.z,
        );
        let distance_m = relative_world.norm_squared().sqrt();
        if !(distance_m > 0.0 && distance_m.is_finite()) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "world observer must not coincide with the specimen origin",
            });
        }
        let relative_body = pose.orientation().rotate_world_to_body(relative_world);
        let direction_body = [relative_body.x, relative_body.y, relative_body.z];
        transfers.clear();
        transfers.try_reserve(self.modes.len()).map_err(|_| {
            StructuralModalBasisError::PressureCapacity {
                requested: self.modes.len(),
            }
        })?;
        for (mode_index, mode) in self.modes.iter().enumerate() {
            if distance_m < mode.minimum_far_field_distance_m {
                return Err(StructuralModalBasisError::ObserverOutsideFarField {
                    distance_m,
                    minimum_m: mode.minimum_far_field_distance_m,
                    mode: mode_index,
                });
            }
            let far = mode.directivity.evaluate(direction_body);
            let phase = C64::new(
                det::cos(mode.wavenumber_rad_m * distance_m),
                det::sin(mode.wavenumber_rad_m * distance_m),
            );
            transfers.push((far * phase).scale(distance_m.recip()));
        }
        Ok(())
    }
}

fn validate_modal_acoustic_directivity(
    directivity: &ModalAcousticDirectivity,
) -> Result<(), StructuralModalBasisError> {
    if directivity.identity == ContentHash([0; 32])
        || directivity.structural_basis_identity == ContentHash([0; 32])
        || directivity.gas_model_identity == ContentHash([0; 32])
        || !(directivity.temperature_k > 0.0
            && directivity.temperature_k.is_finite()
            && directivity.ambient_pressure_pa > 0.0
            && directivity.ambient_pressure_pa.is_finite()
            && directivity.density_kg_m3 > 0.0
            && directivity.density_kg_m3.is_finite()
            && directivity.sound_speed_m_s > 0.0
            && directivity.sound_speed_m_s.is_finite())
        || directivity.modes.is_empty()
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "acoustic directivity header is malformed",
        });
    }
    for (mode_index, mode) in directivity.modes.iter().enumerate() {
        let expected_coefficients = mode
            .directivity
            .l_max
            .checked_add(1)
            .and_then(|side| side.checked_mul(side));
        let expected_wavenumber = mode.angular_frequency_rad_s / directivity.sound_speed_m_s;
        if mode.structural_mode != mode_index
            || !(mode.angular_frequency_rad_s > 0.0 && mode.angular_frequency_rad_s.is_finite())
            || !(mode.wavenumber_rad_m > 0.0 && mode.wavenumber_rad_m.is_finite())
            || mode.wavenumber_rad_m.to_bits() != expected_wavenumber.to_bits()
            || mode.wavenumber_rad_m.to_bits() != mode.directivity.k.to_bits()
            || mode.formulation == HelmholtzFormulation::BurtonMillerWrongAlphaSign
            || mode.directivity.l_max > MAX_SH_DEGREE
            || expected_coefficients != Some(mode.directivity.coefficients.len())
            || !(mode.directivity.captured_fraction >= 0.0
                && mode.directivity.captured_fraction.is_finite()
                && mode.radiated_power_per_modal_velocity_squared >= 0.0
                && mode.radiated_power_per_modal_velocity_squared.is_finite()
                && mode.panels_per_wavelength > 0.0
                && mode.panels_per_wavelength.is_finite()
                && mode.condition_lower_bound >= 1.0
                && mode.condition_lower_bound.is_finite()
                && mode.minimum_far_field_distance_m > 0.0
                && mode.minimum_far_field_distance_m.is_finite())
            || mode
                .directivity
                .coefficients
                .iter()
                .any(|coefficient| !(coefficient.re.is_finite() && coefficient.im.is_finite()))
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "directivity mode rows are malformed or misaligned",
            });
        }
    }
    Ok(())
}

fn validate_acoustic_medium(gas: &GasState) -> Result<(), StructuralModalBasisError> {
    for (value, what) in [
        (
            gas.temperature,
            "gas temperature must be positive and finite",
        ),
        (gas.pressure, "gas pressure must be positive and finite"),
        (gas.density, "gas density must be positive and finite"),
        (
            gas.sound_speed,
            "gas sound speed must be positive and finite",
        ),
        (
            gas.dynamic_viscosity,
            "gas dynamic viscosity must be positive and finite",
        ),
        (
            gas.thermal_conductivity,
            "gas thermal conductivity must be positive and finite",
        ),
        (
            gas.specific_gas_constant,
            "gas specific constant must be positive and finite",
        ),
        (
            gas.specific_heat_cp,
            "gas heat capacity must be positive and finite",
        ),
        (
            gas.prandtl,
            "gas Prandtl number must be positive and finite",
        ),
        (
            gas.characteristic_impedance,
            "gas characteristic impedance must be positive and finite",
        ),
    ] {
        if !(value > 0.0 && value.is_finite()) {
            return Err(StructuralModalBasisError::InvalidRequest { what });
        }
    }
    if !(gas.gamma > 1.0 && gas.gamma.is_finite()) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "gas heat-capacity ratio must be finite and greater than one",
        });
    }
    let expected_impedance = gas.density * gas.sound_speed;
    if (gas.characteristic_impedance - expected_impedance).abs()
        > 32.0 * f64::EPSILON * expected_impedance
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "gas characteristic impedance disagrees with density times sound speed",
        });
    }
    Ok(())
}

fn acoustic_radiation_identity(
    basis: &StructuralModalBasis,
    medium: ResolvedAcousticMedium<'_>,
    observer: AcousticObserver,
    modes: &[AcousticModeRadiation],
) -> ContentHash {
    let mut hasher = DomainHasher::new(MODAL_ACOUSTIC_RADIATION_IDENTITY_DOMAIN);
    hasher.update(basis.identity.as_bytes());
    hasher.update(medium.gas_model_identity.as_bytes());
    for value in [
        medium.gas.temperature,
        medium.gas.pressure,
        medium.gas.density,
        medium.gas.sound_speed,
        medium.gas.dynamic_viscosity,
        medium.gas.thermal_conductivity,
        medium.gas.gamma,
        medium.gas.specific_gas_constant,
        medium.gas.specific_heat_cp,
        medium.gas.prandtl,
        medium.gas.characteristic_impedance,
        observer.position_m[0],
        observer.position_m[1],
        observer.position_m[2],
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    for mode in modes {
        hasher.update(&[match mode.formulation {
            HelmholtzFormulation::PlainCbie => 0,
            HelmholtzFormulation::BurtonMiller => 1,
            HelmholtzFormulation::BurtonMillerWrongAlphaSign => 2,
        }]);
        for value in [
            mode.angular_frequency_rad_s,
            mode.observer_pressure_per_modal_velocity.re,
            mode.observer_pressure_per_modal_velocity.im,
            mode.radiated_power_per_modal_velocity_squared,
            mode.panels_per_wavelength,
            mode.condition_lower_bound,
            mode.minimum_far_field_distance_m,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn acoustic_directivity_identity(
    basis: &StructuralModalBasis,
    medium: ResolvedAcousticMedium<'_>,
    controls: AcousticDirectivityControls,
    modes: &[AcousticModeDirectivity],
) -> ContentHash {
    let mut hasher = DomainHasher::new(MODAL_ACOUSTIC_DIRECTIVITY_IDENTITY_DOMAIN);
    hasher.update(basis.identity.as_bytes());
    hasher.update(medium.gas_model_identity.as_bytes());
    hasher.update(
        &u64::try_from(controls.maximum_spherical_harmonic_degree)
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    for value in [
        controls.minimum_captured_fraction,
        medium.gas.temperature,
        medium.gas.pressure,
        medium.gas.density,
        medium.gas.sound_speed,
        medium.gas.dynamic_viscosity,
        medium.gas.thermal_conductivity,
        medium.gas.gamma,
        medium.gas.specific_gas_constant,
        medium.gas.specific_heat_cp,
        medium.gas.prandtl,
        medium.gas.characteristic_impedance,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    for mode in modes {
        hasher.update(
            &u64::try_from(mode.structural_mode)
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        hasher.update(&[match mode.formulation {
            HelmholtzFormulation::PlainCbie => 0,
            HelmholtzFormulation::BurtonMiller => 1,
            HelmholtzFormulation::BurtonMillerWrongAlphaSign => 2,
        }]);
        hasher.update(
            &u64::try_from(mode.directivity.l_max)
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        for value in [
            mode.angular_frequency_rad_s,
            mode.wavenumber_rad_m,
            mode.directivity.captured_fraction,
            mode.radiated_power_per_modal_velocity_squared,
            mode.panels_per_wavelength,
            mode.condition_lower_bound,
            mode.minimum_far_field_distance_m,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
        for coefficient in &mode.directivity.coefficients {
            hasher.update(&coefficient.re.to_bits().to_le_bytes());
            hasher.update(&coefficient.im.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn validate_request(request: &StructuralModeRequest<'_>) -> Result<(), StructuralModalBasisError> {
    if !(request.minimum_frequency_hz.is_finite()
        && request.minimum_frequency_hz > 0.0
        && request.maximum_frequency_hz.is_finite()
        && request.maximum_frequency_hz > request.minimum_frequency_hz)
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "frequency band must satisfy 0 < minimum < maximum",
        });
    }
    if request.maximum_modes == 0 {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "maximum_modes must be positive",
        });
    }
    Ok(())
}

fn rounded_cylinder_dimensions(
    specimen: &ResolvedElasticDiscProfile,
) -> Result<(f64, f64, f64), StructuralModalBasisError> {
    match specimen.profile.spec {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m,
            thickness_m,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        } => Ok((outer_radius_m, thickness_m, 0.0)),
        DiscProfileSpec::SolidCylinder {
            outer_radius_m,
            thickness_m,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius },
        } => Ok((outer_radius_m, thickness_m, radius)),
        _ => Err(StructuralModalBasisError::UnsupportedProfile {
            what: "exact axisymmetric tetrahedralization is not yet available for this profile family",
        }),
    }
}

fn basis_identity(
    request: &StructuralModeRequest<'_>,
    mesh_spec: RoundedCylinderMeshSpec,
    mesh: &RoundedCylinderTetMesh,
    assembly: &TetElasticAssembly,
    modes: &[StructuralMode],
) -> ContentHash {
    let mut hasher = DomainHasher::new(STRUCTURAL_MODAL_BASIS_IDENTITY_DOMAIN);
    hasher.update(&STRUCTURAL_MODAL_BASIS_SCHEMA_VERSION.to_le_bytes());
    hasher.update(request.specimen.identity.as_bytes());
    for value in [
        mesh_spec.outer_radius_m,
        mesh_spec.thickness_m,
        mesh_spec.fillet_radius_m,
        request.minimum_frequency_hz,
        request.maximum_frequency_hz,
        assembly.total_mass_kg,
        assembly.minimum_scaled_jacobian,
        mesh.maximum_meridian_chord_error_m,
        mesh.maximum_azimuthal_chord_error_m,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    for count in [
        u64::from(mesh_spec.core_radial_segments),
        u64::from(mesh_spec.fillet_radial_segments),
        u64::from(mesh_spec.azimuthal_segments),
        u64::from(mesh_spec.axial_segments),
        u64::try_from(mesh.nodes_m.len()).unwrap_or(u64::MAX),
        u64::try_from(mesh.tetrahedra.len()).unwrap_or(u64::MAX),
        u64::try_from(modes.len()).unwrap_or(u64::MAX),
    ] {
        hasher.update(&count.to_le_bytes());
    }
    for mode in modes {
        for value in [
            mode.eigenvalue_s2,
            mode.eigenvalue_interval_s2.0,
            mode.eigenvalue_interval_s2.1,
            mode.eigenvalue_residual_s2,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
        for displacement in &mode.nodal_shape_per_sqrt_kg {
            for value in displacement {
                hasher.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    hasher.finalize()
}

// Closest point and barycentric coordinates, following the Voronoi-region
// construction from Real-Time Collision Detection. The selected triangle is
// non-degenerate because fs-mesh validates every boundary panel.
fn closest_point_on_triangle(point: [f64; 3], triangle: [[f64; 3]; 3]) -> ([f64; 3], [f64; 3]) {
    let [a, b, c] = triangle;
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(point, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, [1.0, 0.0, 0.0]);
    }
    let bp = sub(point, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, [0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (add_scaled(a, ab, v), [1.0 - v, v, 0.0]);
    }
    let cp = sub(point, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, [0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (add_scaled(a, ac, w), [1.0 - w, 0.0, w]);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let bc = sub(c, b);
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (add_scaled(b, bc, w), [0.0, 1.0 - w, w]);
    }
    let denominator = (va + vb + vc).recip();
    let v = vb * denominator;
    let w = vc * denominator;
    (add_scaled(add_scaled(a, ab, v), ac, w), [1.0 - v - w, v, w])
}

fn add_scaled(a: [f64; 3], direction: [f64; 3], scale: f64) -> [f64; 3] {
    [
        direction[0].mul_add(scale, a[0]),
        direction[1].mul_add(scale, a[1]),
        direction[2].mul_add(scale, a[2]),
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm_squared(value: [f64; 3]) -> f64 {
    dot(value, value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, resolve_isotropic_elastic_state_point,
        resolve_orthotropic_elastic_state_point,
    };
    use fs_mbd::UnitQuaternion;
    use fs_qty::{Density, Dims, Pressure};

    fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x5354_5255_4354_4143,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn material_card(young_modulus_pa: f64, density_kg_m3: f64) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, density_kg_m3),
            ("young_modulus", Pressure::DIMS, young_modulus_pa),
            ("poisson_ratio", Dims::NONE, 0.29),
            ("yield_stress", Pressure::DIMS, 250.0e6),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: format!("structural-acoustics test {name}"),
                        license: "CC0-1.0".to_owned(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test-isotropic-solid".to_owned(),
                phase: "solid".to_owned(),
                process: "synthetic".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn orthotropic_material_card() -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, 8_000.0),
            ("young_modulus_1", Pressure::DIMS, 220.0e9),
            ("young_modulus_2", Pressure::DIMS, 80.0e9),
            ("young_modulus_3", Pressure::DIMS, 30.0e9),
            ("poisson_ratio_12", Dims::NONE, 0.20),
            ("poisson_ratio_13", Dims::NONE, 0.10),
            ("poisson_ratio_23", Dims::NONE, 0.15),
            ("shear_modulus_12", Pressure::DIMS, 50.0e9),
            ("shear_modulus_23", Pressure::DIMS, 20.0e9),
            ("shear_modulus_31", Pressure::DIMS, 25.0e9),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: format!("structural-acoustics orthotropic test {name}"),
                        license: "CC0-1.0".to_owned(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "test-orthotropic-solid".to_owned(),
                phase: "solid".to_owned(),
                process: "synthetic-principal-axis-data".to_owned(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn specimen() -> ResolvedElasticDiscProfile {
        let point = QueryPoint::new().with("T", 293.15).unwrap();
        let material = resolve_isotropic_elastic_state_point(
            &material_card(193.0e9, 8_000.0),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        with_cx(|cx| {
            DiscProfileSpec::SolidCylinder {
                outer_radius_m: 0.038,
                thickness_m: 0.006,
                edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
            }
            .resolve_with_isotropic_elastic_state(&material, cx)
            .unwrap()
        })
    }

    fn coarse_modal_request(specimen: &ResolvedElasticDiscProfile) -> StructuralModeRequest<'_> {
        StructuralModeRequest {
            specimen,
            mesh: StructuralMeshControls {
                core_radial_segments: 2,
                fillet_radial_segments: 1,
                azimuthal_segments: 8,
                axial_segments: 1,
                maximum_vertices: 1_000,
                maximum_tetrahedra: 10_000,
            },
            minimum_frequency_hz: 100.0,
            maximum_frequency_hz: 100_000.0,
            maximum_modes: 64,
            slice: SliceOptions::default(),
            assembly: TetAssemblyBudget::standard(),
        }
    }

    #[test]
    fn closest_point_preserves_barycentric_reconstruction() {
        let triangle = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        for point in [[0.5, 0.25, 1.0], [-1.0, -1.0, 0.0], [2.0, 2.0, 0.0]] {
            let (closest, barycentric) = closest_point_on_triangle(point, triangle);
            let reconstructed = [
                barycentric[0] * triangle[0][0]
                    + barycentric[1] * triangle[1][0]
                    + barycentric[2] * triangle[2][0],
                barycentric[0] * triangle[0][1]
                    + barycentric[1] * triangle[1][1]
                    + barycentric[2] * triangle[2][1],
                barycentric[0] * triangle[0][2]
                    + barycentric[1] * triangle[1][2]
                    + barycentric[2] * triangle[2][2],
            ];
            assert!(norm_squared(sub(closest, reconstructed)) < 1.0e-28);
            assert!((barycentric.iter().sum::<f64>() - 1.0).abs() < 1.0e-14);
            assert!(barycentric.iter().all(|weight| *weight >= 0.0));
        }
    }

    #[test]
    fn g1_resolved_specimen_produces_mass_normalized_modes_and_contact_forces() {
        let specimen = specimen();
        let request = StructuralModeRequest {
            specimen: &specimen,
            mesh: StructuralMeshControls {
                core_radial_segments: 2,
                fillet_radial_segments: 1,
                azimuthal_segments: 8,
                axial_segments: 1,
                maximum_vertices: 1_000,
                maximum_tetrahedra: 10_000,
            },
            minimum_frequency_hz: 100.0,
            maximum_frequency_hz: 100_000.0,
            maximum_modes: 64,
            slice: SliceOptions::default(),
            assembly: TetAssemblyBudget::standard(),
        };
        let basis = with_cx(|cx| build_structural_modal_basis(&request, cx)).unwrap();
        assert_eq!(basis.certified_mode_count, basis.modes.len());
        assert!(!basis.modes.is_empty());
        assert!(basis.assembly.total_mass_kg > 0.0);
        for mode in &basis.modes {
            assert_eq!(
                mode.panel_normal_shape_per_sqrt_kg.len(),
                basis.mesh.boundary.triangles.len()
            );
            let reduced: Vec<f64> = basis
                .assembly
                .free_dofs
                .iter()
                .map(|dof| mode.nodal_shape_per_sqrt_kg[dof / 3][dof % 3])
                .collect();
            let mass = basis.assembly.mass.to_dense();
            let modal_mass = (0..reduced.len())
                .map(|row| {
                    (0..reduced.len())
                        .map(|column| {
                            reduced[row] * mass[row * reduced.len() + column] * reduced[column]
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
            assert!((modal_mass - 1.0).abs() < 1.0e-8, "{modal_mass}");
        }

        let rayleigh = RayleighDamping::new(1.25, 2.5e-8).unwrap();
        let loss = modal_loss_spectrum_from_rayleigh(
            &basis,
            &specimen,
            rayleigh,
            ContentHash([0x48; 32]),
        )
        .unwrap();
        assert_eq!(loss.loss_factors.len(), basis.modes.len());
        for (loss_factor, mode) in loss.loss_factors.iter().zip(&basis.modes) {
            let expected = rayleigh.alpha / mode.angular_frequency_rad_s
                + rayleigh.beta * mode.angular_frequency_rad_s;
            assert_eq!(loss_factor.to_bits(), expected.to_bits());
        }

        let point = [0.038, 0.0, 0.0];
        let one = basis
            .project_point_force(point, [0.0, 0.0, 1.0], 1.0e-12)
            .unwrap();
        let two = basis
            .project_point_force(point, [0.0, 0.0, 2.0], 1.0e-12)
            .unwrap();
        assert!(one.distance_to_boundary_m < 1.0e-14);
        for (one, two) in one
            .modal_force_n_per_sqrt_kg
            .iter()
            .zip(two.modal_force_n_per_sqrt_kg)
        {
            assert!((two - 2.0 * one).abs() < 1.0e-12);
        }
    }

    #[test]
    fn g1_material_axis_orientation_changes_the_actual_modal_and_acoustic_basis() {
        let point = QueryPoint::new().with("T", 293.15).unwrap();
        let material = resolve_orthotropic_elastic_state_point(
            &orthotropic_material_card(),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            1.0e-3,
        )
        .unwrap();
        let profile = DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        };
        let identity_orientation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let principal_one_into_axial = [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
        let (radial_stiff, axial_stiff) = with_cx(|cx| {
            (
                profile
                    .resolve_with_orthotropic_elastic_state(
                        &material,
                        identity_orientation,
                        ContentHash([0x31; 32]),
                        cx,
                    )
                    .unwrap(),
                profile
                    .resolve_with_orthotropic_elastic_state(
                        &material,
                        principal_one_into_axial,
                        ContentHash([0x32; 32]),
                        cx,
                    )
                    .unwrap(),
            )
        });
        assert_ne!(radial_stiff.identity, axial_stiff.identity);
        assert_ne!(
            radial_stiff.elastic_material.stiffness_mandel_pa(),
            axial_stiff.elastic_material.stiffness_mandel_pa()
        );

        let radial_basis =
            with_cx(|cx| build_structural_modal_basis(&coarse_modal_request(&radial_stiff), cx))
                .unwrap();
        let axial_basis =
            with_cx(|cx| build_structural_modal_basis(&coarse_modal_request(&axial_stiff), cx))
                .unwrap();
        assert_ne!(radial_basis.identity, axial_basis.identity);
        assert!(
            radial_basis
                .modes
                .iter()
                .zip(&axial_basis.modes)
                .any(|(radial, axial)| {
                    radial.frequency_hz.to_bits() != axial.frequency_hz.to_bits()
                }),
            "rotating a strongly orthotropic tensor out of the disc plane must alter at least one retained physical mode"
        );
    }

    #[test]
    fn g1_point_force_drives_identity_bound_si_pressure_runtime() {
        let specimen = specimen();
        let request = StructuralModeRequest {
            specimen: &specimen,
            mesh: StructuralMeshControls {
                core_radial_segments: 2,
                fillet_radial_segments: 1,
                azimuthal_segments: 8,
                axial_segments: 1,
                maximum_vertices: 1_000,
                maximum_tetrahedra: 10_000,
            },
            minimum_frequency_hz: 100.0,
            maximum_frequency_hz: 100_000.0,
            maximum_modes: 64,
            slice: SliceOptions::default(),
            assembly: TetAssemblyBudget::standard(),
        };
        let basis = with_cx(|cx| build_structural_modal_basis(&request, cx)).unwrap();
        let loss = ModalLossSpectrum {
            structural_basis_identity: basis.identity,
            material_state_identity: basis.material_state_identity,
            damping_model_identity: ContentHash([0x5a; 32]),
            loss_factors: vec![0.02; basis.modes.len()],
        };
        // This manufactured transfer isolates the structural/contact/time
        // composition. The BEM radiation solver has independent tests; no
        // synthetic value is shipped by the production constructor.
        let radiation = ModalAcousticRadiation {
            structural_basis_identity: basis.identity,
            gas_model_identity: ContentHash([0x6b; 32]),
            temperature_k: 293.15,
            ambient_pressure_pa: 101_325.0,
            density_kg_m3: 1.204,
            sound_speed_m_s: 343.0,
            observer: AcousticObserver {
                position_m: [1.0, 0.0, 0.0],
            },
            modes: basis
                .modes
                .iter()
                .enumerate()
                .map(|(index, mode)| AcousticModeRadiation {
                    structural_mode: index,
                    angular_frequency_rad_s: mode.angular_frequency_rad_s,
                    formulation: HelmholtzFormulation::PlainCbie,
                    observer_pressure_per_modal_velocity: if index == 0 {
                        C64::from_re(1.0)
                    } else {
                        C64::ZERO
                    },
                    radiated_power_per_modal_velocity_squared: 0.0,
                    panels_per_wavelength: 10.0,
                    condition_lower_bound: 1.0,
                    minimum_far_field_distance_m: 0.1,
                })
                .collect(),
            identity: ContentHash([0x7c; 32]),
        };
        let mut runtime = PhysicalModalAudioModel::try_new(
            &basis,
            &loss,
            &radiation,
            500_000,
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();

        let first_mode = &basis.modes[0];
        let mut selected = None;
        for triangle in &basis.mesh.boundary.triangles {
            for &node in triangle {
                let shape = first_mode.nodal_shape_per_sqrt_kg[node];
                let norm = norm_squared(shape);
                if selected
                    .as_ref()
                    .is_none_or(|(_, _, best): &([f64; 3], [f64; 3], f64)| norm > *best)
                {
                    selected = Some((basis.mesh.nodes_m[node], shape, norm));
                }
            }
        }
        let (point, shape, norm) = selected.unwrap();
        assert!(norm > 0.0);
        let force = shape.map(|component| component / norm.sqrt());
        let frame = runtime.step_point_force(point, force, 1.0e-12).unwrap();
        assert!(frame.force_projection.modal_force_n_per_sqrt_kg[0] > 0.0);
        assert_ne!(frame.acoustic.observer_pressure_pa, 0.0);
        assert!(frame.acoustic.total_modal_energy_j > 0.0);
        assert!(
            frame.acoustic.viscous_dissipation_j
                >= -frame.acoustic.dissipation_roundoff_tolerance_j
        );
    }

    #[test]
    fn g0_world_observer_uses_current_body_pose_and_spherical_spreading() {
        // A pure Y_10 table has a body-axis null in the equatorial plane. A
        // quarter-turn about body y therefore rotates a world-z microphone
        // from the polar lobe into that null without any pan/gain scripting.
        let directivity = ModalAcousticDirectivity {
            structural_basis_identity: ContentHash([0x11; 32]),
            gas_model_identity: ContentHash([0x22; 32]),
            temperature_k: 293.15,
            ambient_pressure_pa: 101_325.0,
            density_kg_m3: 1.204,
            sound_speed_m_s: 343.0,
            modes: vec![AcousticModeDirectivity {
                structural_mode: 0,
                angular_frequency_rad_s: 343.0,
                wavenumber_rad_m: 1.0,
                formulation: HelmholtzFormulation::PlainCbie,
                directivity: DirectivityTable {
                    k: 1.0,
                    l_max: 1,
                    coefficients: vec![C64::ZERO, C64::ZERO, C64::from_re(1.0), C64::ZERO],
                    captured_fraction: 1.0,
                },
                radiated_power_per_modal_velocity_squared: 1.0,
                panels_per_wavelength: 10.0,
                condition_lower_bound: 1.0,
                minimum_far_field_distance_m: 1.0,
            }],
            identity: ContentHash([0x33; 32]),
        };
        let observer = AcousticWorldObserver {
            position_world_m: [0.0, 0.0, 10.0],
        };
        let polar = directivity
            .observer_transfers_at_pose(Pose::identity(), observer)
            .unwrap()[0];
        assert!(polar.norm_sq() > 1.0e-6);

        let quarter_turn =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.5 * core::f64::consts::PI)
                .unwrap();
        let rotated_pose = Pose::new(Vec3::ZERO, quarter_turn).unwrap();
        let equatorial = directivity
            .observer_transfers_at_pose(rotated_pose, observer)
            .unwrap()[0];
        assert!(equatorial.norm_sq() < polar.norm_sq() * 1.0e-24);

        let nearer_pose = Pose::new(Vec3::new(0.0, 0.0, 5.0), UnitQuaternion::IDENTITY).unwrap();
        let nearer = directivity
            .observer_transfers_at_pose(nearer_pose, observer)
            .unwrap()[0];
        assert!((nearer.norm_sq() / polar.norm_sq() - 4.0).abs() < 1.0e-12);

        let mut forged_medium = directivity.clone();
        forged_medium.sound_speed_m_s = 344.0;
        assert!(matches!(
            forged_medium.observer_transfers_at_pose(Pose::identity(), observer),
            Err(StructuralModalBasisError::InvalidRequest { .. })
        ));

        let mut forged_far_field = directivity;
        forged_far_field.modes[0].minimum_far_field_distance_m = f64::NAN;
        assert!(matches!(
            forged_far_field.observer_transfers_at_pose(Pose::identity(), observer),
            Err(StructuralModalBasisError::InvalidRequest { .. })
        ));
    }
}
