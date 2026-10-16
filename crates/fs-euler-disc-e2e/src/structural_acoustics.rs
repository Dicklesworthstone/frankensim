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
    DirectivityTable, Formulation as HelmholtzFormulation, HelmholtzError,
    MAX_RADIATION_FIELDS_PER_BATCH, MAX_SH_DEGREE, Medium, RadiationSolution, directivity_sh_table,
    far_field, solve_radiation, solve_radiation_batch,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::{ContentHash, DomainHasher};
use fs_couple::broadband_radiation::{
    BroadbandRadiationArtifact, BroadbandRadiationAuthority, BroadbandRadiationControls,
    BroadbandRadiationError, BroadbandRadiationRuntime, ComplexShTrainingSample,
    DirectFarFieldHeldOutSample, HarmonicTimeConvention, MAX_BROADBAND_FREQUENCIES,
    MAX_VALIDATION_DIRECTIONS, RadiationSampleDiagnostics, RealTesseralChannel,
    SampledRadiationData, build_broadband_radiation_artifact, evaluate_real_tesseral,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticFrame, ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel,
};
use fs_exec::Cx;
use fs_material::gas::GasState;
use fs_material::state_point::IsotropicElasticStatePoint;
use fs_material::visco::{LoweredModel, RayleighDamping, ViscoError, loss_factor_to_zeta};
use fs_math::c64::C64;
use fs_math::det;
use fs_mbd::{Pose, Vec3};
use fs_mesh::{
    RoundedCylinderMeshError, RoundedCylinderMeshSpec, RoundedCylinderTetMesh,
    rounded_cylinder_tet_mesh,
};
use fs_modal::{ModalError, ModePair, SliceOptions, SliceStats, slice_window};
use fs_plate::{
    AssemblyOptions as PlateAssemblyOptions, EdgeSupport as PlateEdgeSupport, PlateError,
    PlateMesh, PlateModel, PlateSection, assemble as assemble_plate,
};
use fs_rep_frep::SquatDiscEdgeTreatment;
use fs_solid::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetLinearElasticProblem,
    TetMaterialField,
};

use crate::audio_resampling::{
    DecimatedModalAcceleration, fixed_rate_frame_count_with_roundoff_bound,
};
use crate::render_trajectory::RenderTrajectory;
use crate::specimen::{DiscProfileSpec, ResolvedElasticDiscProfile};
use crate::timeline_resampling::{EventEvaluationSide, TimelineResampler, TimelineResamplingError};
use crate::{
    AudioResamplingError, ChannelControl, EulerControlStream, GeneralizedForceMeasureInterval,
    GeneralizedForceReconstructionInput, ReconstructedGeneralizedForce,
    reconstruct_generalized_force_measures,
};

/// Schema version of the integrated structural modal artifact.
pub const STRUCTURAL_MODAL_BASIS_SCHEMA_VERSION: u32 = 1;
/// Schema version of the count-certified residual-flexibility estimate.
pub const STRUCTURAL_RESIDUAL_FLEXIBILITY_SCHEMA_VERSION: u32 = 1;
/// Exact limitation attached to every truncated residual-flexibility basis.
pub const STRUCTURAL_RESIDUAL_FLEXIBILITY_NO_CLAIM: &str = "Estimate-only: inertia certifies eigenvalue counts, not an enclosure of the eigenvector-derived compliance; eigenpair residuals, modes above the declared enrichment band, mesh and constitutive error, dynamic correction, and broadband radiation/audio are not bounded";
/// Limitation carried by the broadband body-frame source artifact and stem.
pub const STRUCTURAL_BROADBAND_SOURCE_NO_CLAIM: &str = "Estimate-only linear source stem about the undeformed stationary body: sampled BEM, SH truncation, rational fitting, discretization, structural truncation, and constitutive damping are not certified; static residual flexibility is deliberately excluded, and no 1/r propagation, delay, Doppler, listener/room response, impact, air-film sound, mastering, or calibrated SPL is claimed";
/// Limitation carried by the rigid-disc acceleration radiation artifact.
pub const RIGID_DISC_BROADBAND_SOURCE_NO_CLAIM: &str = "Estimate-only linear low-Mach radiation from rigid center-of-mass translation of the undeformed reference disc: sampled BEM, SH truncation, rational fitting, discretization, moving-boundary, rotational, convective, near-field, room, calibration, and backreaction errors are not certified; rotational coordinates are excluded because their low-frequency boundary-work estimate is not passivity-admissible on the production mesh";
/// Limitation carried by the retarded rigid-source far-field observer.
pub const RETARDED_FAR_FIELD_OBSERVER_NO_CLAIM: &str = "Estimate-only rigid low-Mach far-field observer: causal retarded delay, emission-pose direction, 1/r spreading, and deterministic multi-observer timing are modeled; no moving-boundary/FW-H, convective Green/Jacobian or exact Doppler amplitude, near field, deformation, room/support/head scattering, absorption, backreaction, impact radiation, calibration, or certified far-field enclosure is claimed";
/// Maximum number of simultaneous physical pressure observers in one pass.
pub const MAX_PHYSICAL_PRESSURE_OBSERVERS: usize = 64;
const STRUCTURAL_MODAL_BASIS_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-modal-basis.v1";
const STRUCTURAL_RESIDUAL_FLEXIBILITY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-residual-flexibility.v1";
const STRUCTURAL_RESIDUAL_RESPONSE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-residual-response.v1";
const STRUCTURAL_BROADBAND_SOURCE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-broadband-source.v1";
const STRUCTURAL_BROADBAND_ARTIFACT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-broadband-artifact.v1";
const RIGID_DISC_BROADBAND_SOURCE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.rigid-disc-broadband-source.v1";
const RIGID_DISC_BROADBAND_ARTIFACT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.rigid-disc-broadband-artifact.v1";
const RETARDED_FAR_FIELD_SIGNAL_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.retarded-far-field-signal.v1";
const MODAL_ACOUSTIC_RADIATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.modal-acoustic-radiation.v1";
const MODAL_ACOUSTIC_DIRECTIVITY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.modal-acoustic-directivity.v1";

/// Authority ceiling carried by residual-flexibility artifacts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StructuralResidualFlexibilityAuthority {
    /// Truncated numerical estimate under [`STRUCTURAL_RESIDUAL_FLEXIBILITY_NO_CLAIM`].
    EstimateOnly,
}

impl StructuralResidualFlexibilityAuthority {
    /// Stable artifact spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        "estimate-only"
    }

    const fn tag(self) -> u8 {
        1
    }
}

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

/// Spectral controls for the residual-flexibility complement above a forcing band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructuralResidualFlexibilityControls {
    /// Highest elastic frequency included in the static sum [Hz].
    pub maximum_enrichment_frequency_hz: f64,
    /// Hard cap on certified positive-frequency enrichment modes.
    pub maximum_enrichment_modes: usize,
}

/// Count-certified truncated sum `u_H = sum(phi phi^T f / lambda)` with an
/// Estimate-only authority ceiling under [`STRUCTURAL_RESIDUAL_FLEXIBILITY_NO_CLAIM`].
#[derive(Clone, Debug)]
pub struct StructuralResidualFlexibilityEstimateBasis {
    /// Integrated artifact schema.
    pub schema_version: u32,
    /// Explicit authority ceiling; never inferred from a type name.
    pub authority: StructuralResidualFlexibilityAuthority,
    /// Identity of resolved geometry plus material state.
    pub specimen_identity: ContentHash,
    /// Geometry/density profile identity.
    pub profile_identity: ContentHash,
    /// Elastic material state used to assemble `(K,M)`.
    pub material_state_identity: ContentHash,
    /// Exact body-fitted volume and boundary panelization.
    pub mesh: RoundedCylinderTetMesh,
    /// Physical mass and stiffness operators whose positive spectrum was used.
    pub assembly: TetElasticAssembly,
    /// Exact finite-element assembly budget and quality threshold.
    pub assembly_budget: TetAssemblyBudget,
    /// Exact identity of mesh, panels, reduced operators, and DOF semantics.
    pub operator_identity: ContentHash,
    /// Six analytic free-body modes after verified M-orthonormalization.
    pub rigid_modes_per_sqrt_kg: Vec<Vec<f64>>,
    /// Largest relative stiffness-null residual over the rigid basis.
    pub maximum_rigid_stiffness_relative_residual: f64,
    /// Largest M-inner-product defect among rigid and retained elastic modes.
    pub maximum_mass_orthogonality_error: f64,
    /// Certified mass-normalized enrichment modes with full nodal/panel shapes.
    pub enrichment_modes: Vec<StructuralMode>,
    /// Declared forcing band whose upper edge begins the residual sum [Hz].
    pub forcing_frequency_band_hz: (f64, f64),
    /// Half-open certified enrichment band `(low, high]` [Hz].
    pub enrichment_frequency_band_hz: (f64, f64),
    /// Inertia-certified number of modes in the enrichment window.
    pub certified_enrichment_mode_count: usize,
    /// Inertia-certified number of modes inside the forcing band.
    pub certified_in_band_mode_count: usize,
    /// Exact inertia count over `(forcing_min, enrichment_max]`.
    pub certified_partition_mode_count: usize,
    /// Sparse eigensolver work accounting.
    pub slice_stats: SliceStats,
    /// Content identity binding physical, spectral, and numerical inputs.
    pub identity: ContentHash,
}

/// Static elastic response to one arbitrary moving boundary point force.
#[derive(Clone, Debug)]
pub struct StructuralResidualFlexibilityEstimateResponse {
    /// Residual-flexibility basis used for the response.
    pub basis_identity: ContentHash,
    /// Authority inherited exactly from the basis.
    pub authority: StructuralResidualFlexibilityAuthority,
    /// Requested body-frame force application point [m].
    pub requested_point_m: [f64; 3],
    /// Applied body-frame point force [N].
    pub applied_force_n: [f64; 3],
    /// Boundary interpolation and generalized forces for every enrichment mode.
    pub force_projection: PointForceProjection,
    /// Full nodal load after removal of all six rigid generalized loads [N].
    pub inertia_relieved_nodal_force_n: Vec<[f64; 3]>,
    /// Largest relative rigid generalized-force remainder.
    pub maximum_rigid_force_relative_residual: f64,
    /// Static coordinates `q_k = (phi_k^T f)/lambda_k` [m sqrt(kg)].
    pub modal_displacement_m_sqrt_kg: Vec<f64>,
    /// Full nodal elastic displacement reconstructed from the truncated sum [m].
    pub nodal_displacement_m: Vec<[f64; 3]>,
    /// Boundary-panel normal displacement for broadband radiation [m].
    pub panel_normal_displacement_m: Vec<f64>,
    /// Independently interpolated physical contact work `f^T u_H` [J].
    pub elastic_work_j: f64,
    /// Independently evaluated strain energy `u_H^T K u_H / 2` [J].
    pub recoverable_strain_energy_j: f64,
    /// Closure residual `f^T u_H - u_H^T K u_H` [J].
    pub energy_closure_residual_j: f64,
    /// Identity binding the basis, physical force, projection, and response.
    pub identity: ContentHash,
}

/// Difference metrics for one nested residual-flexibility enrichment study.
#[derive(Clone, Debug)]
pub struct StructuralResidualFlexibilityEstimateComparison {
    /// Coarser response identity.
    pub coarse_response_identity: ContentHash,
    /// Finer response identity.
    pub fine_response_identity: ContentHash,
    /// Signed added static work from the larger spectral window [J].
    pub elastic_work_increment_j: f64,
    /// Absolute work difference normalized by the finer work scale.
    pub relative_elastic_work_difference: f64,
    /// Full-nodal displacement L2 difference normalized by the finer response.
    pub relative_nodal_displacement_l2_difference: f64,
    /// Panel-normal displacement L2 difference normalized by the finer response.
    pub relative_panel_normal_l2_difference: f64,
}

/// Offline BEM sampling and causal-fit request for the residual basis's fixed
/// scalar modal-acceleration coordinates.
pub struct StructuralBroadbandRadiationRequest<'a> {
    /// Count-certified source basis; only `enrichment_modes` are dynamic.
    pub basis: &'a StructuralResidualFlexibilityEstimateBasis,
    /// Loss factors evaluated in the exact enrichment-mode order.
    pub loss: &'a StructuralResidualModalLossSpectrum,
    /// Evidence-bound exterior gas state.
    pub medium: ResolvedAcousticMedium<'a>,
    /// Strictly increasing positive training frequencies [Hz].
    pub training_frequency_hz: &'a [f64],
    /// Strictly increasing, disjoint validation frequencies [Hz].
    pub held_out_frequency_hz: &'a [f64],
    /// Independent nonzero body-frame directions for direct held-out BEM data.
    pub held_out_directions_body: &'a [[f64; 3]],
    /// Spherical-harmonic truncation controls.
    pub directivity: AcousticDirectivityControls,
    /// Rational-fit, validation, and audio-clock controls.
    pub fit: BroadbandRadiationControls,
}

/// Constitutive damping samples bound specifically to one residual basis's
/// enrichment-mode order.
#[derive(Clone, Debug, PartialEq)]
pub struct StructuralResidualModalLossSpectrum {
    /// Exact residual basis identity.
    pub structural_basis_identity: ContentHash,
    /// Exact elastic material-state identity.
    pub material_state_identity: ContentHash,
    /// Constitutive damping model/evidence identity.
    pub damping_model_identity: ContentHash,
    /// Loss factor `eta(omega_k)` in enrichment-mode order.
    pub loss_factors: Vec<f64>,
}

/// Identity-bound body-frame radiation bank from modal acceleration to
/// far-field source amplitude. Each filter input has units
/// `m sqrt(kg) / s^2`; each output has units `Pa m` after runtime evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct StructuralBroadbandRadiationArtifact {
    /// Complete identity of the fitted, discretized physical source bank.
    pub identity: ContentHash,
    /// Evidence ceiling inherited from the neutral broadband bridge.
    pub authority: BroadbandRadiationAuthority,
    /// Exact residual/enrichment basis identity.
    pub structural_basis_identity: ContentHash,
    /// Gas species/EOS/transport identity used by BEM.
    pub gas_model_identity: ContentHash,
    /// Constant homogeneous propagation speed [m/s].
    pub sound_speed_m_s: f64,
    /// Maximum undeformed source extent [m].
    pub source_diameter_m: f64,
    /// Exact constitutive damping identity.
    pub damping_model_identity: ContentHash,
    /// Audio sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Loss factors in enrichment-mode/input order.
    pub modal_loss_factors: Vec<f64>,
    /// Solver-neutral fitted real-tesseral filter bank.
    pub radiation: BroadbandRadiationArtifact,
    /// Explicit applicability boundary.
    pub no_claims: &'static str,
}

/// Admitted body-frame translational coordinates for the rigid disc.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RigidDiscAcousticCoordinate {
    /// Center-of-mass translation along body x [m/s2].
    TranslationX,
    /// Center-of-mass translation along body y [m/s2].
    TranslationY,
    /// Center-of-mass translation along body z [m/s2].
    TranslationZ,
}

const RIGID_DISC_ACOUSTIC_COORDINATES: [RigidDiscAcousticCoordinate; 3] = [
    RigidDiscAcousticCoordinate::TranslationX,
    RigidDiscAcousticCoordinate::TranslationY,
    RigidDiscAcousticCoordinate::TranslationZ,
];

/// Offline BEM sampling and causal-fit request for rigid-disc acceleration.
pub struct RigidDiscBroadbandRadiationRequest<'a> {
    /// Exact disc mesh and profile identity used for the boundary surface.
    pub basis: &'a StructuralResidualFlexibilityEstimateBasis,
    /// Evidence-bound exterior gas state.
    pub medium: ResolvedAcousticMedium<'a>,
    /// Strictly increasing positive training frequencies [Hz].
    pub training_frequency_hz: &'a [f64],
    /// Strictly increasing disjoint validation frequencies [Hz].
    pub held_out_frequency_hz: &'a [f64],
    /// Independent nonzero body-frame validation directions.
    pub held_out_directions_body: &'a [[f64; 3]],
    /// Spherical-harmonic truncation controls.
    pub directivity: AcousticDirectivityControls,
    /// Rational-fit, validation, and audio-clock controls.
    pub fit: BroadbandRadiationControls,
}

/// Identity-bound body-frame bank from rigid generalized acceleration to
/// far-field source amplitude in pascal-metres.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidDiscBroadbandRadiationArtifact {
    /// Complete fitted-source identity.
    pub identity: ContentHash,
    /// Estimate-only authority inherited from the broadband bridge.
    pub authority: BroadbandRadiationAuthority,
    /// Exact disc mesh/basis identity.
    pub structural_basis_identity: ContentHash,
    /// Gas model used by the Helmholtz solves.
    pub gas_model_identity: ContentHash,
    /// Constant homogeneous propagation speed [m/s].
    pub sound_speed_m_s: f64,
    /// Maximum undeformed source extent [m].
    pub source_diameter_m: f64,
    /// Audio sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Fixed input order consumed by the runtime.
    pub coordinates: Vec<RigidDiscAcousticCoordinate>,
    /// Solver-neutral fitted real-tesseral filter bank.
    pub radiation: BroadbandRadiationArtifact,
    /// Explicit applicability boundary.
    pub no_claims: &'static str,
}

/// One real body-frame far-field source coefficient in pascal-metres.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FarFieldSourceCoefficientPaM(
    /// Coefficient value [Pa m].
    pub f64,
);

/// Fixed-rate body-frame real-tesseral source stem, before listener
/// propagation. This is intentionally not a [`PhysicalPressureSignal`].
#[derive(Clone, Debug, PartialEq)]
pub struct StructuralBroadbandSourceStem {
    /// Opening boundary of the first sample cell [s]. Frame zero is at
    /// `start_time_s + 1/sample_rate_hz`.
    pub start_time_s: f64,
    /// Fixed sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Canonical coefficient channel order.
    pub channels: Vec<RealTesseralChannel>,
    /// Frame-major coefficients [Pa m].
    pub coefficients: Vec<FarFieldSourceCoefficientPaM>,
    /// Evidence ceiling.
    pub authority: BroadbandRadiationAuthority,
    /// Exact broadband radiation artifact that produced the coefficients.
    pub source_identity: ContentHash,
    /// Exact structural basis that supplied generalized acceleration.
    pub structural_basis_identity: ContentHash,
}

impl StructuralBroadbandSourceStem {
    /// Number of complete coefficient frames.
    #[must_use]
    pub fn frame_count(&self) -> usize {
        self.coefficients.len() / self.channels.len()
    }

    /// One complete coefficient frame in canonical channel order.
    #[must_use]
    pub fn frame(&self, index: usize) -> Option<&[FarFieldSourceCoefficientPaM]> {
        let start = index.checked_mul(self.channels.len())?;
        self.coefficients.get(start..start + self.channels.len())
    }
}

/// Controls for the deterministic rigid-source retarded far-field evaluator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedFarFieldObserverControls {
    /// Maximum admitted source-surface speed divided by sound speed.
    pub maximum_surface_mach: f64,
    /// Absolute emission-time bisection tolerance [s].
    pub root_time_tolerance_s: f64,
    /// Maximum deterministic bisection iterations.
    pub maximum_root_iterations: u32,
    /// Lanczos radius; version one requires exactly eight (16 taps).
    pub interpolation_radius_frames: u8,
    /// Maximum common output frames across all observers.
    pub maximum_output_frames: usize,
}

/// Geometric support constraint for a rectangular thin plate.
///
/// Coordinates are expressed in the centered plate frame, with the plate
/// occupying `[-width/2,width/2] x [-depth/2,depth/2]`. Point supports are
/// snapped to the nearest structured-mesh nodes under an explicit tolerance;
/// the resolved nodes and maximum snap error become part of the modal-basis
/// artifact and its identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RectangularPlateSupport {
    /// Apply one idealized support condition to every perimeter node.
    Perimeter(PlateEdgeSupport),
    /// Pin transverse displacement at three non-collinear support locations.
    /// Plate rotations remain free at each pin.
    ThreePointPinned {
        /// Requested support locations in centered plate coordinates [m].
        points_centered_m: [[f64; 2]; 3],
        /// Largest permitted request-to-mesh-node distance [m].
        maximum_snap_distance_m: f64,
    },
}

impl RectangularPlateSupport {
    /// Stable, human-readable support-family code for artifact manifests.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Perimeter(PlateEdgeSupport::SimplySupported) => "simply-supported-perimeter",
            Self::Perimeter(PlateEdgeSupport::Clamped) => "clamped-perimeter",
            Self::ThreePointPinned { .. } => "three-point-pinned",
        }
    }

    /// Validate support geometry and mesh snapping for a rectangular grid.
    ///
    /// # Errors
    /// Refuses non-finite dimensions, an empty grid, out-of-bounds points,
    /// excessive snap error, duplicate resolved nodes, or collinear pins.
    pub fn validate_for_rectangular_grid(
        self,
        width_m: f64,
        depth_m: f64,
        cells_x: usize,
        cells_y: usize,
    ) -> Result<(), StructuralModalBasisError> {
        if !(width_m.is_finite() && width_m > 0.0)
            || !(depth_m.is_finite() && depth_m > 0.0)
            || cells_x == 0
            || cells_y == 0
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "plate support validation requires positive dimensions and grid cells",
            });
        }
        let mesh = PlateMesh::rectangle(width_m, depth_m, cells_x, cells_y);
        resolve_rectangular_plate_support(self, &mesh, width_m, depth_m, cells_x, cells_y)
            .map(|_| ())
    }
}

#[derive(Clone, Debug)]
struct ResolvedRectangularPlateSupport {
    boundary_nodes: Vec<usize>,
    condition: PlateEdgeSupport,
    maximum_snap_error_m: f64,
}

fn resolve_rectangular_plate_support(
    support: RectangularPlateSupport,
    mesh: &PlateMesh,
    width_m: f64,
    depth_m: f64,
    cells_x: usize,
    cells_y: usize,
) -> Result<ResolvedRectangularPlateSupport, StructuralModalBasisError> {
    match support {
        RectangularPlateSupport::Perimeter(condition) => Ok(ResolvedRectangularPlateSupport {
            boundary_nodes: PlateMesh::rectangle_boundary(cells_x, cells_y),
            condition,
            maximum_snap_error_m: 0.0,
        }),
        RectangularPlateSupport::ThreePointPinned {
            points_centered_m,
            maximum_snap_distance_m,
        } => {
            if !(maximum_snap_distance_m.is_finite() && maximum_snap_distance_m >= 0.0)
                || points_centered_m
                    .iter()
                    .flatten()
                    .any(|coordinate| !coordinate.is_finite())
            {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "plate point supports and snap tolerance must be finite and nonnegative",
                });
            }
            let half_width = 0.5 * width_m;
            let half_depth = 0.5 * depth_m;
            if points_centered_m.iter().any(|point| {
                point[0] < -half_width
                    || point[0] > half_width
                    || point[1] < -half_depth
                    || point[1] > half_depth
            }) {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "plate point support lies outside the rectangular plate",
                });
            }

            let mut boundary_nodes = Vec::with_capacity(3);
            let mut maximum_snap_error_m = 0.0_f64;
            for point in points_centered_m {
                let local = [point[0] + half_width, point[1] + half_depth];
                let mut nearest = None::<(usize, f64)>;
                for (node, &(x, y)) in mesh.nodes.iter().enumerate() {
                    let dx = x - local[0];
                    let dy = y - local[1];
                    let distance_squared = dx.mul_add(dx, dy * dy);
                    if nearest.is_none_or(|(_, best)| distance_squared < best) {
                        nearest = Some((node, distance_squared));
                    }
                }
                let (node, distance_squared) =
                    nearest.ok_or(StructuralModalBasisError::InvalidRequest {
                        what: "plate point support cannot resolve against an empty mesh",
                    })?;
                let distance_m = distance_squared.sqrt();
                if distance_m > maximum_snap_distance_m {
                    return Err(StructuralModalBasisError::InvalidRequest {
                        what: "plate point support exceeds its mesh-snap tolerance",
                    });
                }
                maximum_snap_error_m = maximum_snap_error_m.max(distance_m);
                boundary_nodes.push(node);
            }
            let mut distinct_nodes = boundary_nodes.clone();
            distinct_nodes.sort_unstable();
            distinct_nodes.dedup();
            if distinct_nodes.len() != 3 {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "plate point supports must resolve to three distinct mesh nodes",
                });
            }
            let resolved_points =
                [boundary_nodes[0], boundary_nodes[1], boundary_nodes[2]].map(|node| {
                    let (x, y) = mesh.nodes[node];
                    [x - half_width, y - half_depth]
                });
            let [a, b, c] = resolved_points;
            let twice_area = (b[0] - a[0]).mul_add(c[1] - a[1], -(b[1] - a[1]) * (c[0] - a[0]));
            let area_scale = det::powi(width_m.max(depth_m), 2);
            if twice_area.abs() <= 64.0 * f64::EPSILON * area_scale {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "plate point supports must be non-collinear",
                });
            }
            Ok(ResolvedRectangularPlateSupport {
                boundary_nodes,
                condition: PlateEdgeSupport::SimplySupported,
                maximum_snap_error_m,
            })
        }
    }
}

/// Physical and numerical request for a supported rectangular thin plate.
///
/// Geometry, boundary support, and the resolved material state are all
/// explicit. A display name cannot select stiffness, density, dimensions, or
/// support conditions.
pub struct RectangularPlateModeRequest<'a> {
    /// Plate span along the base-frame x axis [m].
    pub width_m: f64,
    /// Plate span along the base-frame y axis [m].
    pub depth_m: f64,
    /// Uniform plate thickness [m].
    pub thickness_m: f64,
    /// Exact evidence-bearing isotropic elastic state.
    pub elastic: &'a IsotropicElasticStatePoint,
    /// Explicit geometric support constraint admitted by the plate model.
    pub support: RectangularPlateSupport,
    /// Structured cells along x.
    pub cells_x: usize,
    /// Structured cells along y.
    pub cells_y: usize,
    /// Maximum admitted mesh nodes.
    pub maximum_nodes: usize,
    /// Strictly positive lower edge of the requested band [Hz].
    pub minimum_frequency_hz: f64,
    /// Upper edge of the requested band [Hz].
    pub maximum_frequency_hz: f64,
    /// Maximum certified in-band modes retained.
    pub maximum_modes: usize,
    /// Certified sparse modal-slice controls.
    pub slice: SliceOptions,
}

/// One mass-normalized supported-plate bending mode.
#[derive(Clone, Debug)]
pub struct RectangularPlateMode {
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
    /// Transverse displacement at every full plate node [kg^-1/2].
    pub nodal_transverse_shape_per_sqrt_kg: Vec<f64>,
    /// Plate-x slope at every full plate node [m^-1 kg^-1/2].
    pub nodal_slope_x_per_m_sqrt_kg: Vec<f64>,
    /// Plate-y slope at every full plate node [m^-1 kg^-1/2].
    pub nodal_slope_y_per_m_sqrt_kg: Vec<f64>,
    /// Shared nodal mixed derivative used by the C1 rectangular Hermite
    /// reconstruction [m^-2 kg^-1/2].
    pub nodal_mixed_slope_per_m2_sqrt_kg: Vec<f64>,
}

/// Certified supported-plate basis used for contact-force projection and
/// acoustic radiation.
#[derive(Clone, Debug)]
pub struct RectangularPlateModalBasis {
    /// Width [m].
    pub width_m: f64,
    /// Depth [m].
    pub depth_m: f64,
    /// Thickness [m].
    pub thickness_m: f64,
    /// Exact material state used by the section.
    pub material_state_identity: ContentHash,
    /// Requested support geometry and condition.
    pub support: RectangularPlateSupport,
    /// Mesh nodes to which support constraints were actually applied.
    pub support_node_indices: Vec<usize>,
    /// Largest request-to-node support snap error [m].
    pub maximum_support_snap_error_m: f64,
    /// Structured cells along x.
    pub cells_x: usize,
    /// Structured cells along y.
    pub cells_y: usize,
    /// Structured plate mesh in `[0,width] x [0,depth]` coordinates.
    pub mesh: PlateMesh,
    /// Assembled reduced `(K,M)` pencil and full-to-reduced DOF map.
    pub model: PlateModel,
    /// Certified retained bending modes.
    pub modes: Vec<RectangularPlateMode>,
    /// Certified requested eigenvalue window [s^-2].
    pub eigenvalue_window_s2: (f64, f64),
    /// Inertia-certified count in the requested window.
    pub certified_mode_count: usize,
    /// Sparse eigensolver work accounting.
    pub slice_stats: SliceStats,
    /// Identity binding geometry, support, material state, mesh, and modes.
    pub identity: ContentHash,
}

/// Modal projection of one transverse point force on the supported plate.
#[derive(Clone, Debug, PartialEq)]
pub struct PlatePointForceProjection {
    /// Zero-based structured cell containing the application point.
    pub cell: [usize; 2],
    /// Normalized coordinates within the selected rectangular cell.
    pub cell_coordinates: [f64; 2],
    /// Generalized force for every retained mode [N kg^-1/2].
    pub modal_force_n_per_sqrt_kg: Vec<f64>,
}

/// Frequency-domain and causal retarded-time Rayleigh transfers for fixed
/// observers above one rigidly baffled supported plate.
#[derive(Clone, Debug)]
pub struct BaffledPlateObserverRadiation {
    /// Structural plate basis consumed by the transfer.
    pub structural_basis_identity: ContentHash,
    /// Gas model/species identity.
    pub gas_model_identity: ContentHash,
    /// World-space observer positions [m].
    pub observers: Vec<AcousticWorldObserver>,
    /// Per-observer, per-mode pressure divided by generalized modal velocity.
    pub pressure_per_modal_velocity: Vec<Vec<C64>>,
    /// Per-observer, per-mode causal FIR mapping sampled generalized modal
    /// acceleration to pressure [Pa / (m sqrt(kg) s^-2)]. Lag zero is the
    /// current sample boundary; later entries are progressively older states.
    pub pressure_per_modal_acceleration_fir: Vec<Vec<Vec<f64>>>,
    /// Sampling rate at which the retarded-time FIR was constructed [Hz].
    pub retarded_sample_rate_hz: u32,
    /// Longest retained propagation delay including its interpolation tap.
    pub maximum_retarded_delay_frames: usize,
    /// Smallest acoustic quadrature panels-per-wavelength over retained modes.
    pub minimum_panels_per_wavelength: f64,
    /// Identity binding the basis, gas, observers, quadrature, and transfers.
    pub identity: ContentHash,
}

/// Physical time-domain runtime for one fixed supported plate.
pub struct BaffledPlateModalAudioModel<'basis> {
    basis: &'basis RectangularPlateModalBasis,
    runtime: ModalAcousticTimeModel,
    modal_damping_ratios: Vec<f64>,
    sample_rate_hz: u32,
    maximum_abs_pressure_pa: f64,
    radiation_identity: ContentHash,
    damping_model_identity: ContentHash,
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
    /// Authoritative interval force-time measures are conservatively
    /// rasterized onto the audio clock and low-pass reconstructed before
    /// structural integration. Spatial projection still uses the closing
    /// contact endpoint, falling back to the opening endpoint when required.
    IntervalMeasureAtClosingElseOpeningEndpointBandLimitedV1,
    /// Actual mechanics-cadence moving-contact modal acceleration was filtered
    /// before each factor-two sample removal, then causal delay was compensated
    /// using real mechanics postroll. No second structural integration occurs.
    MechanicsModalAccelerationAntiAliasedDecimatedV1,
}

/// Explicit modal state at the first retained mechanics boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhysicalModalInitialState {
    /// Every retained displacement and velocity begins at zero. This is
    /// appropriate only when the modeled excitation truly begins at the
    /// source horizon.
    Zero,
    /// Begin at static equilibrium under the first held generalized load.
    /// This avoids inventing an impact when a causal preroll starts during an
    /// already-running contact experiment; subsequent force motion still
    /// excites the modes dynamically.
    StaticEquilibriumAtFirstHeldForce,
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

impl PhysicalPressureSignal {
    /// Extract a contiguous pressure interval and rebase its first boundary.
    ///
    /// This is the generic causal-preroll publication seam: modal state may be
    /// warmed on a longer trajectory, while the published signal remains an
    /// exact sample crop. No resampling, gain, filtering, or fade is applied.
    pub fn try_crop_rebased(
        &self,
        first_frame: usize,
        end_frame: usize,
        rebased_start_time_s: f64,
    ) -> Result<Self, StructuralModalBasisError> {
        if first_frame >= end_frame
            || end_frame > self.pressure_pa.len()
            || !rebased_start_time_s.is_finite()
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "physical pressure crop must be finite, nonempty, and in bounds",
            });
        }
        let pressure_pa = self.pressure_pa[first_frame..end_frame].to_vec();
        let peak_abs_pressure_pa = pressure_pa
            .iter()
            .fold(0.0_f64, |peak, pressure| peak.max(pressure.abs()));
        let mut identity = DomainHasher::new("org.frankensim.euler-disc.physical-pressure-crop.v1");
        identity.update(self.identity.as_bytes());
        identity.update(&(first_frame as u64).to_le_bytes());
        identity.update(&(end_frame as u64).to_le_bytes());
        identity.update(&rebased_start_time_s.to_bits().to_le_bytes());
        for pressure in &pressure_pa {
            identity.update(&pressure.to_bits().to_le_bytes());
        }
        Ok(Self {
            start_time_s: rebased_start_time_s,
            sample_rate_hz: self.sample_rate_hz,
            pressure_pa,
            peak_abs_pressure_pa,
            contact_force_sampling: self.contact_force_sampling,
            observer: self.observer,
            structural_basis_identity: self.structural_basis_identity,
            radiation_identity: self.radiation_identity,
            damping_model_identity: self.damping_model_identity,
            identity: identity.finalize(),
        })
    }
}

/// Superpose simultaneous physical pressure fields at one observer.
///
/// Linear acoustic pressure is additive. This operation stays in pascals and
/// therefore belongs before any pressure-to-digital listening gain. Inputs
/// are ordered by their complete signal identities so the result is invariant
/// to caller iteration order; the composite structural, damping, and
/// radiation identities bind every contributing model.
///
/// # Errors
/// Refuses an empty set, mismatched clocks/observers/sampling conventions,
/// non-finite source or summed pressure, allocation failure, or cancellation.
pub fn superpose_pressure_signals(
    signals: &[&PhysicalPressureSignal],
    cx: &Cx<'_>,
) -> Result<PhysicalPressureSignal, StructuralModalBasisError> {
    let first = signals
        .first()
        .copied()
        .ok_or(StructuralModalBasisError::InvalidRequest {
            what: "pressure superposition requires at least one signal",
        })?;
    let mut ordered = Vec::new();
    ordered.try_reserve_exact(signals.len()).map_err(|_| {
        StructuralModalBasisError::PressureCapacity {
            requested: signals.len(),
        }
    })?;
    ordered.extend_from_slice(signals);
    ordered.sort_by(|left, right| left.identity.as_bytes().cmp(right.identity.as_bytes()));

    for signal in &ordered {
        if signal.start_time_s.to_bits() != first.start_time_s.to_bits()
            || signal.sample_rate_hz != first.sample_rate_hz
            || signal.pressure_pa.len() != first.pressure_pa.len()
            || signal.contact_force_sampling != first.contact_force_sampling
            || signal.observer != first.observer
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "superposed pressure fields must share clock, observer, length, and force-sampling semantics",
            });
        }
    }

    let mut pressure_pa = Vec::new();
    pressure_pa
        .try_reserve_exact(first.pressure_pa.len())
        .map_err(|_| StructuralModalBasisError::PressureCapacity {
            requested: first.pressure_pa.len(),
        })?;
    let mut peak_abs_pressure_pa = 0.0_f64;
    for frame in 0..first.pressure_pa.len() {
        if frame % 4_096 == 0 {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
        }
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        for signal in &ordered {
            let value = signal.pressure_pa[frame];
            if !value.is_finite() {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "pressure superposition source sample must be finite",
                });
            }
            let corrected = value - correction;
            let next = sum + corrected;
            correction = (next - sum) - corrected;
            sum = next;
        }
        if !sum.is_finite() {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "pressure superposition produced a non-finite sample",
            });
        }
        peak_abs_pressure_pa = peak_abs_pressure_pa.max(sum.abs());
        pressure_pa.push(sum);
    }
    cx.checkpoint()
        .map_err(|_| StructuralModalBasisError::Cancelled)?;

    let composite_identity =
        |domain: &'static str, select: fn(&PhysicalPressureSignal) -> ContentHash| {
            let mut hasher = DomainHasher::new(domain);
            let mut identities: Vec<_> = ordered.iter().map(|signal| select(signal)).collect();
            identities.sort_by(|left, right| left.as_bytes().cmp(right.as_bytes()));
            for identity in identities {
                hasher.update(identity.as_bytes());
            }
            hasher.finalize()
        };
    let structural_basis_identity = composite_identity(
        "org.frankensim.euler-disc.pressure-superposition.structural-bases.v1",
        |signal| signal.structural_basis_identity,
    );
    let radiation_identity = composite_identity(
        "org.frankensim.euler-disc.pressure-superposition.radiation.v1",
        |signal| signal.radiation_identity,
    );
    let damping_model_identity = composite_identity(
        "org.frankensim.euler-disc.pressure-superposition.damping.v1",
        |signal| signal.damping_model_identity,
    );
    let mut identity =
        DomainHasher::new("org.frankensim.euler-disc.pressure-superposition.signal.v1");
    for signal in &ordered {
        identity.update(signal.identity.as_bytes());
    }
    identity.update(structural_basis_identity.as_bytes());
    identity.update(radiation_identity.as_bytes());
    identity.update(damping_model_identity.as_bytes());
    for pressure in &pressure_pa {
        identity.update(&pressure.to_bits().to_le_bytes());
    }

    Ok(PhysicalPressureSignal {
        start_time_s: first.start_time_s,
        sample_rate_hz: first.sample_rate_hz,
        pressure_pa,
        peak_abs_pressure_pa,
        contact_force_sampling: first.contact_force_sampling,
        observer: first.observer,
        structural_basis_identity,
        radiation_identity,
        damping_model_identity,
        identity: identity.finalize(),
    })
}

/// Integrated physical runtime: one structural basis, state-dependent modal
/// damping, and one BEM-derived observer transfer.
pub struct PhysicalModalAudioModel<'basis> {
    basis: &'basis StructuralModalBasis,
    runtime: ModalAcousticTimeModel,
    static_pressure_transfers: Vec<C64>,
    sample_rate_hz: u32,
    static_observer: Option<AcousticObserver>,
    /// Identity of the acoustic radiation artifact.
    pub radiation_identity: ContentHash,
    /// Identity of the damping model and its evidence.
    pub damping_model_identity: ContentHash,
}

/// Persistent enrichment-mode oscillator plus causal broadband radiation bank.
pub struct StructuralBroadbandSourceRuntime<'artifact> {
    basis: &'artifact StructuralResidualFlexibilityEstimateBasis,
    source: &'artifact StructuralBroadbandRadiationArtifact,
    modal_runtime: ModalAcousticTimeModel,
    modal_damping_ratios: Vec<f64>,
    radiation_runtime: BroadbandRadiationRuntime<'artifact>,
    modal_acceleration: Vec<f64>,
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
    /// The supported thin-plate assembly or modal solve refused.
    Plate(PlateError),
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
    /// Broadband fitting or persistent radiation filtering refused.
    BroadbandRadiation(BroadbandRadiationError),
    /// Conservative generalized-force reconstruction refused.
    ForceReconstruction(AudioResamplingError),
    /// Pose reconstruction for a world-fixed observer refused.
    Timeline(TimelineResamplingError),
    /// The requested band contains no elastic modes.
    NoModesInBand,
    /// A rigid-subspace, inertia-relief, or energy-closure check failed.
    ResidualFlexibilityVerification {
        /// Failed invariant.
        what: &'static str,
        /// Dimensionless or SI residual named by `what`.
        residual: f64,
        /// Corresponding admitted tolerance.
        tolerance: f64,
    },
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
    /// A broadband sample failed passivity after its diagnostic far field was evaluated.
    BroadbandNegativeRadiatedPower(String),
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
    /// The plate surface quadrature is too coarse for its shortest retained
    /// acoustic wavelength.
    PlateRadiationUnderresolved {
        /// Smallest cells per wavelength over both plate axes.
        panels_per_wavelength: f64,
        /// Required minimum.
        minimum: f64,
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
            Self::Plate(source) => write!(formatter, "supported plate refused: {source}"),
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
            Self::BroadbandRadiation(source) => {
                write!(
                    formatter,
                    "structural broadband radiation refused: {source:?}"
                )
            }
            Self::ForceReconstruction(source) => {
                write!(
                    formatter,
                    "physical generalized-force reconstruction refused: {source}"
                )
            }
            Self::Timeline(source) => {
                write!(
                    formatter,
                    "physical audio pose reconstruction refused: {source}"
                )
            }
            Self::NoModesInBand => write!(formatter, "FS-EULER-STRUCTURAL-MODE-EMPTY-BAND"),
            Self::ResidualFlexibilityVerification {
                what,
                residual,
                tolerance,
            } => write!(
                formatter,
                "FS-EULER-STRUCTURAL-RESIDUAL-VERIFY: {what} residual {residual:.6e} exceeds {tolerance:.6e}"
            ),
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
            Self::BroadbandNegativeRadiatedPower(diagnostic) => {
                write!(formatter, "FS-EULER-ACOUSTIC-NONPASSIVE: {diagnostic}")
            }
            Self::DirectivityTruncation {
                mode,
                captured_fraction,
                minimum_fraction,
            } => write!(
                formatter,
                "FS-EULER-ACOUSTIC-DIRECTIVITY-TRUNCATION: mode {mode} captured {captured_fraction:.6e}, below required {minimum_fraction:.6e}"
            ),
            Self::PlateRadiationUnderresolved {
                panels_per_wavelength,
                minimum,
            } => write!(
                formatter,
                "FS-EULER-PLATE-RADIATION-RESOLUTION: {panels_per_wavelength:.6e} panels/wavelength is below {minimum:.6e}"
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
            Self::Plate(source) => Some(source),
            Self::Modal(source) => Some(source),
            Self::BemSurface(source) => Some(source),
            Self::Acoustic(source) => Some(source),
            Self::Viscoelastic(source) => Some(source),
            Self::ModalTime(source) => Some(source),
            Self::BroadbandRadiation(_) => None,
            Self::ForceReconstruction(source) => Some(source),
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

impl From<PlateError> for StructuralModalBasisError {
    fn from(source: PlateError) -> Self {
        Self::Plate(source)
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

impl From<BroadbandRadiationError> for StructuralModalBasisError {
    fn from(source: BroadbandRadiationError) -> Self {
        Self::BroadbandRadiation(source)
    }
}

impl From<AudioResamplingError> for StructuralModalBasisError {
    fn from(source: AudioResamplingError) -> Self {
        Self::ForceReconstruction(source)
    }
}

impl From<TimelineResamplingError> for StructuralModalBasisError {
    fn from(source: TimelineResamplingError) -> Self {
        Self::Timeline(source)
    }
}

fn reconstruct_projected_modal_forces(
    intervals: &[crate::AudioControlInterval],
    modal_forces: &[Vec<f64>],
    input: GeneralizedForceReconstructionInput,
    cx: &Cx<'_>,
) -> Result<ReconstructedGeneralizedForce, StructuralModalBasisError> {
    if intervals.len() != modal_forces.len() {
        return Err(StructuralModalBasisError::ControlTimeline {
            what: "projected modal forces do not match mechanics intervals",
        });
    }
    let mut measures = Vec::new();
    measures.try_reserve_exact(intervals.len()).map_err(|_| {
        StructuralModalBasisError::PressureCapacity {
            requested: intervals.len(),
        }
    })?;
    for (interval, mean_force) in intervals.iter().zip(modal_forces) {
        let mut force_time_measure = Vec::new();
        force_time_measure
            .try_reserve_exact(mean_force.len())
            .map_err(|_| StructuralModalBasisError::PressureCapacity {
                requested: mean_force.len(),
            })?;
        force_time_measure.extend(mean_force.iter().map(|force| force * interval.duration_s));
        measures.push(GeneralizedForceMeasureInterval {
            start_time_s: interval.start_time_s,
            end_time_s: interval.end_time_s,
            force_time_measure,
        });
    }
    reconstruct_generalized_force_measures(&measures, input, cx).map_err(Into::into)
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
    let (mesh_spec, mesh, assembly) = assemble_structural_operators(request, cx)?;

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

    let modes = structural_modes_from_pairs(&mesh, &assembly, &report.modes)?;

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

/// Build a count-certified positive-spectrum complement above the forcing band.
///
/// # Errors
/// Refuses malformed inputs, incomplete partitions/rigid proofs, empty enrichment, caps, or non-positive modes.
pub fn build_structural_residual_flexibility_estimate_basis(
    request: &StructuralModeRequest<'_>,
    controls: StructuralResidualFlexibilityControls,
    cx: &Cx<'_>,
) -> Result<StructuralResidualFlexibilityEstimateBasis, StructuralModalBasisError> {
    validate_request(request)?;
    if !(controls.maximum_enrichment_frequency_hz.is_finite()
        && controls.maximum_enrichment_frequency_hz > request.maximum_frequency_hz)
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "residual-flexibility enrichment maximum must be finite and exceed the forcing-band maximum",
        });
    }
    if controls.maximum_enrichment_modes == 0 {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "residual-flexibility enrichment mode cap must be positive",
        });
    }
    let (_mesh_spec, mesh, assembly) = assemble_structural_operators(request, cx)?;
    let angular_low = core::f64::consts::TAU * request.minimum_frequency_hz;
    let angular_partition = core::f64::consts::TAU * request.maximum_frequency_hz;
    let angular_high = core::f64::consts::TAU * controls.maximum_enrichment_frequency_hz;
    let eigenvalue_window_s2 = (angular_low * angular_low, angular_high * angular_high);
    let report = slice_window(
        &assembly.stiffness,
        &assembly.mass,
        eigenvalue_window_s2,
        &request.slice,
    )?;
    verify_residual(
        "exactly six eigenvalues below the positive-frequency window",
        report.below_low.abs_diff(6) as f64,
        0.0,
    )?;
    let partition_s2 = angular_partition * angular_partition;
    let mut in_band_count = 0;
    for pair in &report.modes {
        if pair.interval.0 <= eigenvalue_window_s2.0 || pair.interval.1 > eigenvalue_window_s2.1 {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "residual-flexibility eigenvalue interval escapes its requested window",
            });
        }
        if pair.interval.1 <= partition_s2 {
            in_band_count += 1;
        } else if pair.interval.0 <= partition_s2 {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "residual-flexibility eigenvalue interval crosses the forcing-band edge",
            });
        }
    }
    let enrichment_count = report.expected - in_band_count;
    if enrichment_count == 0 {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "residual-flexibility enrichment window contains no elastic modes",
        });
    }
    if in_band_count > request.maximum_modes {
        return Err(StructuralModalBasisError::ModeBudgetExceeded {
            requested: in_band_count,
            maximum: request.maximum_modes,
        });
    }
    if enrichment_count > controls.maximum_enrichment_modes {
        return Err(StructuralModalBasisError::ModeBudgetExceeded {
            requested: enrichment_count,
            maximum: controls.maximum_enrichment_modes,
        });
    }
    let mut all_modes = structural_modes_from_pairs(&mesh, &assembly, &report.modes)?;
    let (rigid_modes_per_sqrt_kg, rigid_null_residual, mass_orthogonality_error) =
        verified_rigid_subspace(&mesh, &assembly, all_modes.iter())?;
    let operator_identity = structural_operator_identity(&mesh, &assembly);
    let enrichment_modes = all_modes.split_off(in_band_count);
    let mut basis = StructuralResidualFlexibilityEstimateBasis {
        schema_version: STRUCTURAL_RESIDUAL_FLEXIBILITY_SCHEMA_VERSION,
        authority: StructuralResidualFlexibilityAuthority::EstimateOnly,
        specimen_identity: request.specimen.identity,
        profile_identity: request.specimen.profile.content_identities().profile,
        material_state_identity: request.specimen.material_state_identity,
        mesh,
        assembly,
        assembly_budget: request.assembly,
        operator_identity,
        rigid_modes_per_sqrt_kg,
        maximum_rigid_stiffness_relative_residual: rigid_null_residual,
        maximum_mass_orthogonality_error: mass_orthogonality_error,
        enrichment_modes,
        forcing_frequency_band_hz: (request.minimum_frequency_hz, request.maximum_frequency_hz),
        enrichment_frequency_band_hz: (
            request.maximum_frequency_hz,
            controls.maximum_enrichment_frequency_hz,
        ),
        certified_enrichment_mode_count: enrichment_count,
        certified_in_band_mode_count: in_band_count,
        certified_partition_mode_count: report.expected,
        slice_stats: report.stats,
        identity: ContentHash([0; 32]),
    };
    basis.identity = basis.recomputed_identity();
    Ok(basis)
}

fn assemble_structural_operators(
    request: &StructuralModeRequest<'_>,
    cx: &Cx<'_>,
) -> Result<
    (
        RoundedCylinderMeshSpec,
        RoundedCylinderTetMesh,
        TetElasticAssembly,
    ),
    StructuralModalBasisError,
> {
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
    Ok((mesh_spec, mesh, assembly))
}

fn structural_modes_from_pairs(
    mesh: &RoundedCylinderTetMesh,
    assembly: &TetElasticAssembly,
    pairs: &[ModePair],
) -> Result<Vec<StructuralMode>, StructuralModalBasisError> {
    let mut modes = Vec::with_capacity(pairs.len());
    for (mode_index, pair) in pairs.iter().enumerate() {
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
                let displacement = core::array::from_fn(|component| {
                    (nodal_shape[triangle[0]][component]
                        + nodal_shape[triangle[1]][component]
                        + nodal_shape[triangle[2]][component])
                        / 3.0
                });
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
    Ok(modes)
}

fn verified_rigid_subspace<'a>(
    mesh: &RoundedCylinderTetMesh,
    assembly: &TetElasticAssembly,
    elastic_modes: impl Iterator<Item = &'a StructuralMode>,
) -> Result<(Vec<Vec<f64>>, f64, f64), StructuralModalBasisError> {
    let n = assembly.free_dofs.len();
    let mut rigid = Vec::<Vec<f64>>::with_capacity(6);
    for kind in 0..6 {
        let mut candidate = Vec::with_capacity(n);
        for &full_dof in &assembly.free_dofs {
            let node = full_dof / 3;
            let component = full_dof % 3;
            let [x, y, z] = mesh.nodes_m[node];
            candidate.push(match kind {
                0..=2 => {
                    if component == kind {
                        1.0
                    } else {
                        0.0
                    }
                }
                3 => [0.0, -z, y][component],
                4 => [z, 0.0, -x][component],
                _ => [-y, x, 0.0][component],
            });
        }
        for _ in 0..2 {
            for prior in &rigid {
                let mut mass_candidate = vec![0.0; n];
                assembly.mass.spmv(&candidate, &mut mass_candidate);
                let coefficient = dot_slice(prior, &mass_candidate);
                for (value, basis_value) in candidate.iter_mut().zip(prior) {
                    *value = coefficient.mul_add(-basis_value, *value);
                }
            }
        }
        let mut mass_candidate = vec![0.0; n];
        assembly.mass.spmv(&candidate, &mut mass_candidate);
        let norm = dot_slice(&candidate, &mass_candidate).sqrt();
        if !(norm.is_finite() && norm > 0.0) {
            return Err(StructuralModalBasisError::ResidualFlexibilityVerification {
                what: "analytic rigid basis rank",
                residual: norm,
                tolerance: f64::MIN_POSITIVE,
            });
        }
        candidate.iter_mut().for_each(|value| *value /= norm);
        rigid.push(candidate);
    }

    let stiffness_scale = (0..assembly.stiffness.nrows())
        .map(|row| assembly.stiffness.row(row).1.iter().map(|v| v.abs()).sum())
        .fold(0.0_f64, f64::max);
    let mut maximum_null_residual = 0.0_f64;
    let mut maximum_orthogonality_error = 0.0_f64;
    for (i, mode) in rigid.iter().enumerate() {
        let mut stiffness_mode = vec![0.0; n];
        assembly.stiffness.spmv(mode, &mut stiffness_mode);
        maximum_null_residual = maximum_null_residual.max(
            maximum_abs(&stiffness_mode)
                / (stiffness_scale * maximum_abs(mode)).max(f64::MIN_POSITIVE),
        );
        let mut mass_mode = vec![0.0; n];
        assembly.mass.spmv(mode, &mut mass_mode);
        for (j, other) in rigid.iter().enumerate() {
            maximum_orthogonality_error = maximum_orthogonality_error
                .max((dot_slice(other, &mass_mode) - if i == j { 1.0 } else { 0.0 }).abs());
        }
    }
    for mode in elastic_modes {
        let reduced: Vec<f64> = assembly
            .free_dofs
            .iter()
            .map(|dof| mode.nodal_shape_per_sqrt_kg[dof / 3][dof % 3])
            .collect();
        let mut mass_mode = vec![0.0; n];
        assembly.mass.spmv(&reduced, &mut mass_mode);
        for rigid_mode in &rigid {
            maximum_orthogonality_error =
                maximum_orthogonality_error.max(dot_slice(rigid_mode, &mass_mode).abs());
        }
    }
    verify_residual(
        "analytic rigid stiffness nullspace",
        maximum_null_residual,
        1.0e-10,
    )?;
    verify_residual(
        "rigid/elastic M orthogonality",
        maximum_orthogonality_error,
        1.0e-8,
    )?;
    Ok((rigid, maximum_null_residual, maximum_orthogonality_error))
}

fn verify_residual(
    what: &'static str,
    residual: f64,
    tolerance: f64,
) -> Result<(), StructuralModalBasisError> {
    if residual.is_finite() && residual <= tolerance {
        Ok(())
    } else {
        Err(StructuralModalBasisError::ResidualFlexibilityVerification {
            what,
            residual,
            tolerance,
        })
    }
}

fn dot_slice(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn maximum_abs(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |max, value| max.max(value.abs()))
}

/// Assemble a certified bending basis for one resolved supported plate.
///
/// # Errors
/// Refuses invalid geometry/discretization, an unsupported material state,
/// plate assembly or modal failures, empty bands, or a certified count above
/// the caller's mode budget.
pub fn build_rectangular_plate_modal_basis(
    request: &RectangularPlateModeRequest<'_>,
    cx: &Cx<'_>,
) -> Result<RectangularPlateModalBasis, StructuralModalBasisError> {
    for (value, what) in [
        (request.width_m, "plate width must be finite and positive"),
        (request.depth_m, "plate depth must be finite and positive"),
        (
            request.thickness_m,
            "plate thickness must be finite and positive",
        ),
        (
            request.minimum_frequency_hz,
            "plate minimum frequency must be finite and positive",
        ),
        (
            request.maximum_frequency_hz,
            "plate maximum frequency must be finite and positive",
        ),
    ] {
        if !(value.is_finite() && value > 0.0) {
            return Err(StructuralModalBasisError::InvalidRequest { what });
        }
    }
    if request.maximum_frequency_hz <= request.minimum_frequency_hz
        || request.cells_x < 2
        || request.cells_y < 2
        || request.maximum_modes == 0
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "plate frequency band, mesh cells, and mode cap must be ordered and nonzero",
        });
    }
    let node_count = request
        .cells_x
        .checked_add(1)
        .and_then(|x| {
            request
                .cells_y
                .checked_add(1)
                .and_then(|y| x.checked_mul(y))
        })
        .ok_or(StructuralModalBasisError::InvalidRequest {
            what: "plate mesh node count overflowed",
        })?;
    if node_count > request.maximum_nodes {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "plate mesh exceeds the caller node budget",
        });
    }
    cx.checkpoint()
        .map_err(|_| StructuralModalBasisError::Cancelled)?;
    let section = PlateSection::isotropic(
        request.elastic.young_modulus_pa(),
        request.elastic.poisson_ratio(),
        request.thickness_m,
        request.elastic.density_kg_m3(),
    )?;
    let mesh = PlateMesh::rectangle(
        request.width_m,
        request.depth_m,
        request.cells_x,
        request.cells_y,
    );
    let support = resolve_rectangular_plate_support(
        request.support,
        &mesh,
        request.width_m,
        request.depth_m,
        request.cells_x,
        request.cells_y,
    )?;
    let model = assemble_plate(
        &mesh,
        &section,
        &support.boundary_nodes,
        &[],
        &PlateAssemblyOptions {
            pretension: 0.0,
            support: support.condition,
        },
    )?;
    let angular_min = core::f64::consts::TAU * request.minimum_frequency_hz;
    let angular_max = core::f64::consts::TAU * request.maximum_frequency_hz;
    let eigenvalue_window_s2 = (angular_min * angular_min, angular_max * angular_max);
    let report = fs_plate::modes(&model, eigenvalue_window_s2, &request.slice)?;
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
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        if !(pair.lambda > 0.0 && pair.interval.0 > 0.0 && pair.interval.1.is_finite()) {
            return Err(StructuralModalBasisError::NonPositiveCertifiedMode { mode: mode_index });
        }
        let mut nodal_transverse = vec![0.0; mesh.nodes.len()];
        let mut nodal_slope_x = vec![0.0; mesh.nodes.len()];
        let mut nodal_slope_y = vec![0.0; mesh.nodes.len()];
        for node in 0..mesh.nodes.len() {
            if let Some(reduced) = model.dof_map[3 * node] {
                nodal_transverse[node] = pair.phi[reduced];
            }
            if let Some(reduced) = model.dof_map[3 * node + 1] {
                nodal_slope_x[node] = pair.phi[reduced];
            }
            if let Some(reduced) = model.dof_map[3 * node + 2] {
                nodal_slope_y[node] = pair.phi[reduced];
            }
        }
        let nodal_mixed_slope = rectangular_mixed_slopes(
            &nodal_slope_x,
            &nodal_slope_y,
            request.cells_x,
            request.cells_y,
            request.width_m,
            request.depth_m,
        );
        let angular_frequency_rad_s = pair.lambda.sqrt();
        modes.push(RectangularPlateMode {
            eigenvalue_s2: pair.lambda,
            angular_frequency_rad_s,
            frequency_hz: angular_frequency_rad_s / core::f64::consts::TAU,
            eigenvalue_interval_s2: pair.interval,
            eigenvalue_residual_s2: pair.residual,
            nodal_transverse_shape_per_sqrt_kg: nodal_transverse,
            nodal_slope_x_per_m_sqrt_kg: nodal_slope_x,
            nodal_slope_y_per_m_sqrt_kg: nodal_slope_y,
            nodal_mixed_slope_per_m2_sqrt_kg: nodal_mixed_slope,
        });
    }
    let material_state_identity = request.elastic.resolved().identity();
    let mut identity =
        DomainHasher::new("org.frankensim.euler-disc.supported-rectangular-plate-modal-basis.v2");
    identity.update(material_state_identity.as_bytes());
    for value in [request.width_m, request.depth_m, request.thickness_m] {
        identity.update(&value.to_bits().to_le_bytes());
    }
    match request.support {
        RectangularPlateSupport::Perimeter(condition) => {
            identity.update(&[0]);
            identity.update(&[match condition {
                PlateEdgeSupport::SimplySupported => 0,
                PlateEdgeSupport::Clamped => 1,
            }]);
        }
        RectangularPlateSupport::ThreePointPinned {
            points_centered_m,
            maximum_snap_distance_m,
        } => {
            identity.update(&[1]);
            for point in points_centered_m {
                for coordinate in point {
                    identity.update(&coordinate.to_bits().to_le_bytes());
                }
            }
            identity.update(&maximum_snap_distance_m.to_bits().to_le_bytes());
        }
    }
    for node in &support.boundary_nodes {
        identity.update(&u64::try_from(*node).unwrap_or(u64::MAX).to_le_bytes());
    }
    identity.update(&support.maximum_snap_error_m.to_bits().to_le_bytes());
    identity.update(
        &u64::try_from(request.cells_x)
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    identity.update(
        &u64::try_from(request.cells_y)
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    for &(x, y) in &mesh.nodes {
        identity.update(&x.to_bits().to_le_bytes());
        identity.update(&y.to_bits().to_le_bytes());
    }
    for triangle in &mesh.tris {
        for node in triangle {
            identity.update(&u64::try_from(*node).unwrap_or(u64::MAX).to_le_bytes());
        }
    }
    for mode in &modes {
        for value in [
            mode.eigenvalue_s2,
            mode.eigenvalue_interval_s2.0,
            mode.eigenvalue_interval_s2.1,
            mode.eigenvalue_residual_s2,
        ] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        for value in &mode.nodal_transverse_shape_per_sqrt_kg {
            identity.update(&value.to_bits().to_le_bytes());
        }
        for values in [
            &mode.nodal_slope_x_per_m_sqrt_kg,
            &mode.nodal_slope_y_per_m_sqrt_kg,
            &mode.nodal_mixed_slope_per_m2_sqrt_kg,
        ] {
            for value in values {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    Ok(RectangularPlateModalBasis {
        width_m: request.width_m,
        depth_m: request.depth_m,
        thickness_m: request.thickness_m,
        material_state_identity,
        support: request.support,
        support_node_indices: support.boundary_nodes,
        maximum_support_snap_error_m: support.maximum_snap_error_m,
        cells_x: request.cells_x,
        cells_y: request.cells_y,
        mesh,
        model,
        modes,
        eigenvalue_window_s2,
        certified_mode_count: report.expected,
        slice_stats: report.stats,
        identity: identity.finalize(),
    })
}

fn rectangular_mixed_slopes(
    slope_x: &[f64],
    slope_y: &[f64],
    cells_x: usize,
    cells_y: usize,
    width_m: f64,
    depth_m: f64,
) -> Vec<f64> {
    let row = cells_x + 1;
    let dx = width_m / cells_x as f64;
    let dy = depth_m / cells_y as f64;
    let derivative = |values: &[f64], i: usize, j: usize, along_x: bool| {
        if along_x {
            let (left, right, span) = if i == 0 {
                (j * row, j * row + 1, dx)
            } else if i == cells_x {
                (j * row + i - 1, j * row + i, dx)
            } else {
                (j * row + i - 1, j * row + i + 1, 2.0 * dx)
            };
            (values[right] - values[left]) / span
        } else {
            let (lower, upper, span) = if j == 0 {
                (i, row + i, dy)
            } else if j == cells_y {
                ((j - 1) * row + i, j * row + i, dy)
            } else {
                ((j - 1) * row + i, (j + 1) * row + i, 2.0 * dy)
            };
            (values[upper] - values[lower]) / span
        }
    };
    let mut mixed = vec![0.0; slope_x.len()];
    for j in 0..=cells_y {
        for i in 0..=cells_x {
            let index = j * row + i;
            // Both estimates represent d2w/(dx dy). Averaging keeps the
            // reconstruction symmetric in the two coordinate directions and
            // gives every adjacent cell the exact same nodal derivative.
            mixed[index] =
                0.5 * (derivative(slope_x, i, j, false) + derivative(slope_y, i, j, true));
        }
    }
    mixed
}

fn cubic_hermite_basis(t: f64) -> [f64; 4] {
    let t2 = t * t;
    let t3 = t2 * t;
    [
        2.0 * t3 - 3.0 * t2 + 1.0,
        t3 - 2.0 * t2 + t,
        -2.0 * t3 + 3.0 * t2,
        t3 - t2,
    ]
}

fn rectangular_mode_value(
    mode: &RectangularPlateMode,
    lower_left: usize,
    row: usize,
    u: f64,
    v: f64,
    cell_width_m: f64,
    cell_depth_m: f64,
) -> f64 {
    let nodes = [
        lower_left,
        lower_left + 1,
        lower_left + row,
        lower_left + row + 1,
    ];
    let hu = cubic_hermite_basis(u);
    let hv = cubic_hermite_basis(v);
    let interpolate_y = |lower: usize, upper: usize, values: &[f64], derivatives: &[f64]| {
        hv[0] * values[lower]
            + hv[1] * cell_depth_m * derivatives[lower]
            + hv[2] * values[upper]
            + hv[3] * cell_depth_m * derivatives[upper]
    };
    let value_left = interpolate_y(
        nodes[0],
        nodes[2],
        &mode.nodal_transverse_shape_per_sqrt_kg,
        &mode.nodal_slope_y_per_m_sqrt_kg,
    );
    let value_right = interpolate_y(
        nodes[1],
        nodes[3],
        &mode.nodal_transverse_shape_per_sqrt_kg,
        &mode.nodal_slope_y_per_m_sqrt_kg,
    );
    let slope_left = interpolate_y(
        nodes[0],
        nodes[2],
        &mode.nodal_slope_x_per_m_sqrt_kg,
        &mode.nodal_mixed_slope_per_m2_sqrt_kg,
    );
    let slope_right = interpolate_y(
        nodes[1],
        nodes[3],
        &mode.nodal_slope_x_per_m_sqrt_kg,
        &mode.nodal_mixed_slope_per_m2_sqrt_kg,
    );
    hu[0] * value_left
        + hu[1] * cell_width_m * slope_left
        + hu[2] * value_right
        + hu[3] * cell_width_m * slope_right
}

impl RectangularPlateModalBasis {
    /// Project a base-frame transverse point force through a C1 rectangular
    /// Hermite reconstruction of the DKT nodal displacement and slope DOFs.
    ///
    /// Using the retained slope DOFs is essential for a moving load: the old
    /// P1 projection was continuous but had artificial slope jumps at every
    /// triangle edge, which injected mesh-crossing harmonics into both the
    /// production support dynamics and the acoustic replay.
    pub fn project_transverse_point_force(
        &self,
        point_base_m: [f64; 3],
        transverse_force_n: f64,
        maximum_surface_distance_m: f64,
    ) -> Result<PlatePointForceProjection, StructuralModalBasisError> {
        if point_base_m.iter().any(|value| !value.is_finite())
            || !transverse_force_n.is_finite()
            || !(maximum_surface_distance_m.is_finite() && maximum_surface_distance_m >= 0.0)
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "plate point force and projection tolerance must be finite",
            });
        }
        if point_base_m[2].abs() > maximum_surface_distance_m {
            return Err(StructuralModalBasisError::ContactOutsideTolerance {
                distance_m: point_base_m[2].abs(),
                tolerance_m: maximum_surface_distance_m,
            });
        }
        let x = point_base_m[0] + 0.5 * self.width_m;
        let y = point_base_m[1] + 0.5 * self.depth_m;
        let scale = self.width_m.max(self.depth_m).max(1.0);
        let tolerance = 128.0 * f64::EPSILON * scale;
        if x < -tolerance
            || x > self.width_m + tolerance
            || y < -tolerance
            || y > self.depth_m + tolerance
        {
            return Err(StructuralModalBasisError::ContactOutsideTolerance {
                distance_m: (-x)
                    .max(x - self.width_m)
                    .max(-y)
                    .max(y - self.depth_m)
                    .max(0.0),
                tolerance_m: tolerance,
            });
        }
        let gx = (x.clamp(0.0, self.width_m) / self.width_m) * self.cells_x as f64;
        let gy = (y.clamp(0.0, self.depth_m) / self.depth_m) * self.cells_y as f64;
        let cell_x = (gx.floor() as usize).min(self.cells_x - 1);
        let cell_y = (gy.floor() as usize).min(self.cells_y - 1);
        let u = (gx - cell_x as f64).clamp(0.0, 1.0);
        let v = (gy - cell_y as f64).clamp(0.0, 1.0);
        let row = self.cells_x + 1;
        let lower_left = cell_y * row + cell_x;
        let cell_width_m = self.width_m / self.cells_x as f64;
        let cell_depth_m = self.depth_m / self.cells_y as f64;
        let modal_force_n_per_sqrt_kg = self
            .modes
            .iter()
            .map(|mode| {
                let shape =
                    rectangular_mode_value(mode, lower_left, row, u, v, cell_width_m, cell_depth_m);
                transverse_force_n * shape
            })
            .collect();
        Ok(PlatePointForceProjection {
            cell: [cell_x, cell_y],
            cell_coordinates: [u, v],
            modal_force_n_per_sqrt_kg,
        })
    }

    /// Evaluate the baffled Rayleigh surface integral at fixed world observers
    /// for every retained structural mode.
    ///
    /// The plate lies in world `z = 0`, centered at the origin, and radiates
    /// into the upper half-space. Triangle-centroid quadrature is admitted only
    /// with at least six cells per shortest retained wavelength.
    pub fn baffled_observer_radiation(
        &self,
        medium: ResolvedAcousticMedium<'_>,
        observers: &[AcousticWorldObserver],
        minimum_panels_per_wavelength: f64,
        retarded_sample_rate_hz: u32,
        cx: &Cx<'_>,
    ) -> Result<BaffledPlateObserverRadiation, StructuralModalBasisError> {
        validate_acoustic_medium(medium.gas)?;
        if medium.gas_model_identity == ContentHash([0; 32])
            || observers.is_empty()
            || observers.len() > MAX_PHYSICAL_PRESSURE_OBSERVERS
            || retarded_sample_rate_hz == 0
            || !(minimum_panels_per_wavelength.is_finite() && minimum_panels_per_wavelength >= 2.0)
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "plate radiation medium, observer count, and resolution floor must be valid",
            });
        }
        if observers.iter().any(|observer| {
            observer
                .position_world_m
                .iter()
                .any(|value| !value.is_finite())
                || observer.position_world_m[2] <= 0.0
        }) {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "baffled-plate observers must be finite and above the plate plane",
            });
        }
        let largest_cell_m =
            (self.width_m / self.cells_x as f64).max(self.depth_m / self.cells_y as f64);
        let maximum_omega = self
            .modes
            .last()
            .expect("admitted plate basis has modes")
            .angular_frequency_rad_s;
        let shortest_wavelength_m = core::f64::consts::TAU * medium.gas.sound_speed / maximum_omega;
        let observed_resolution = shortest_wavelength_m / largest_cell_m;
        if observed_resolution < minimum_panels_per_wavelength {
            return Err(StructuralModalBasisError::PlateRadiationUnderresolved {
                panels_per_wavelength: observed_resolution,
                minimum: minimum_panels_per_wavelength,
            });
        }
        let mut pressure_per_modal_velocity = vec![Vec::new(); observers.len()];
        let mut pressure_per_modal_acceleration_fir =
            vec![vec![Vec::new(); self.modes.len()]; observers.len()];
        for transfers in &mut pressure_per_modal_velocity {
            transfers.try_reserve_exact(self.modes.len()).map_err(|_| {
                StructuralModalBasisError::PressureCapacity {
                    requested: self.modes.len(),
                }
            })?;
        }
        let mut maximum_retarded_delay_frames = 0_usize;
        for (mode_index, mode) in self.modes.iter().enumerate() {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
            let k = mode.angular_frequency_rad_s / medium.gas.sound_speed;
            for (observer_index, observer) in observers.iter().enumerate() {
                let mut integral = C64::ZERO;
                let acceleration_fir =
                    &mut pressure_per_modal_acceleration_fir[observer_index][mode_index];
                for triangle in &self.mesh.tris {
                    let a = self.mesh.nodes[triangle[0]];
                    let b = self.mesh.nodes[triangle[1]];
                    let c = self.mesh.nodes[triangle[2]];
                    let twice_area = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
                    let area_m2 = 0.5 * twice_area.abs();
                    let centroid = [
                        (a.0 + b.0 + c.0) / 3.0 - 0.5 * self.width_m,
                        (a.1 + b.1 + c.1) / 3.0 - 0.5 * self.depth_m,
                        0.0,
                    ];
                    let dx = observer.position_world_m[0] - centroid[0];
                    let dy = observer.position_world_m[1] - centroid[1];
                    let dz = observer.position_world_m[2];
                    let distance_m = (dx * dx + dy * dy + dz * dz).sqrt();
                    let shape = (mode.nodal_transverse_shape_per_sqrt_kg[triangle[0]]
                        + mode.nodal_transverse_shape_per_sqrt_kg[triangle[1]]
                        + mode.nodal_transverse_shape_per_sqrt_kg[triangle[2]])
                        / 3.0;
                    let phase = C64::new(det::cos(k * distance_m), det::sin(k * distance_m));
                    integral = integral + phase.scale(area_m2 * shape / distance_m);

                    let delay_frames =
                        distance_m / medium.gas.sound_speed * f64::from(retarded_sample_rate_hz);
                    let delay_floor = delay_frames.floor();
                    if !(delay_floor >= 0.0 && delay_floor <= (usize::MAX - 2) as f64) {
                        return Err(StructuralModalBasisError::InvalidRequest {
                            what: "baffled-plate propagation delay exceeds addressable history",
                        });
                    }
                    let lag = delay_floor as usize;
                    let fraction = delay_frames - delay_floor;
                    let required = lag + 2;
                    if acceleration_fir.len() < required {
                        acceleration_fir.resize(required, 0.0);
                    }
                    maximum_retarded_delay_frames = maximum_retarded_delay_frames.max(required);
                    // Under this module's exp(-i omega t) convention the
                    // existing frequency transfer is +i rho omega times
                    // velocity. Therefore its causal time-domain equivalent
                    // is -rho times retarded normal acceleration.
                    let coefficient = -medium.gas.density * area_m2 * shape
                        / (2.0 * core::f64::consts::PI * distance_m);
                    acceleration_fir[lag] += coefficient * (1.0 - fraction);
                    acceleration_fir[lag + 1] += coefficient * fraction;
                }
                // Rayleigh I for a velocity-prescribed source in an infinite
                // rigid baffle: p = i rho omega/(2 pi) INT v exp(i k R)/R dS.
                let coefficient = medium.gas.density * mode.angular_frequency_rad_s
                    / (2.0 * core::f64::consts::PI);
                let transfer = C64::new(-integral.im, integral.re).scale(coefficient);
                if !transfer.re.is_finite() || !transfer.im.is_finite() {
                    return Err(StructuralModalBasisError::InvalidRequest {
                        what: "baffled-plate Rayleigh transfer became non-finite",
                    });
                }
                pressure_per_modal_velocity[observer_index].push(transfer);
            }
        }
        if maximum_retarded_delay_frames == 0
            || pressure_per_modal_acceleration_fir.iter().any(|observer| {
                observer.iter().any(|fir| {
                    fir.is_empty() || fir.iter().any(|coefficient| !coefficient.is_finite())
                })
            })
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "baffled-plate retarded-time radiation kernel is empty or non-finite",
            });
        }
        let mut identity =
            DomainHasher::new("org.frankensim.euler-disc.baffled-plate-observer-radiation.v1");
        identity.update(self.identity.as_bytes());
        identity.update(medium.gas_model_identity.as_bytes());
        for value in [
            medium.gas.temperature,
            medium.gas.pressure,
            medium.gas.density,
            medium.gas.sound_speed,
            minimum_panels_per_wavelength,
            observed_resolution,
        ] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        identity.update(&retarded_sample_rate_hz.to_le_bytes());
        identity.update(&(maximum_retarded_delay_frames as u64).to_le_bytes());
        for observer in observers {
            for value in observer.position_world_m {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        for transfers in &pressure_per_modal_velocity {
            for transfer in transfers {
                identity.update(&transfer.re.to_bits().to_le_bytes());
                identity.update(&transfer.im.to_bits().to_le_bytes());
            }
        }
        for observer in &pressure_per_modal_acceleration_fir {
            for fir in observer {
                identity.update(&(fir.len() as u64).to_le_bytes());
                for coefficient in fir {
                    identity.update(&coefficient.to_bits().to_le_bytes());
                }
            }
        }
        Ok(BaffledPlateObserverRadiation {
            structural_basis_identity: self.identity,
            gas_model_identity: medium.gas_model_identity,
            observers: observers.to_vec(),
            pressure_per_modal_velocity,
            pressure_per_modal_acceleration_fir,
            retarded_sample_rate_hz,
            maximum_retarded_delay_frames,
            minimum_panels_per_wavelength: observed_resolution,
            identity: identity.finalize(),
        })
    }
}

impl<'basis> BaffledPlateModalAudioModel<'basis> {
    /// Bind one supported-plate basis, Rayleigh loss law, and fixed-observer
    /// Rayleigh radiation artifact into an exact-ZOH modal runtime. The
    /// frequency-domain transfers remain diagnostic; published samples use
    /// the artifact's causal retarded-acceleration FIR.
    pub fn try_new(
        basis: &'basis RectangularPlateModalBasis,
        damping: RayleighDamping,
        damping_model_identity: ContentHash,
        radiation: &BaffledPlateObserverRadiation,
        sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
    ) -> Result<Self, StructuralModalBasisError> {
        if damping_model_identity == ContentHash([0; 32])
            || radiation.identity == ContentHash([0; 32])
            || radiation.structural_basis_identity != basis.identity
            || radiation.pressure_per_modal_velocity.len() != radiation.observers.len()
            || radiation.pressure_per_modal_acceleration_fir.len() != radiation.observers.len()
            || radiation.retarded_sample_rate_hz != sample_rate_hz
            || radiation.maximum_retarded_delay_frames == 0
            || radiation
                .pressure_per_modal_velocity
                .iter()
                .any(|transfers| transfers.len() != basis.modes.len())
            || radiation
                .pressure_per_modal_acceleration_fir
                .iter()
                .any(|kernels| {
                    kernels.len() != basis.modes.len()
                        || kernels.iter().any(|kernel| {
                            kernel.is_empty()
                                || kernel.len() > radiation.maximum_retarded_delay_frames
                                || kernel.iter().any(|coefficient| !coefficient.is_finite())
                        })
                })
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "plate damping/radiation artifacts do not match the structural basis",
            });
        }
        let mut modes = Vec::with_capacity(basis.modes.len());
        let mut modal_damping_ratios = Vec::with_capacity(basis.modes.len());
        for mode in &basis.modes {
            let damping_ratio = damping.zeta_at(mode.angular_frequency_rad_s);
            if !(damping_ratio.is_finite() && damping_ratio >= 0.0) {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "plate modal damping ratio must be finite and non-negative",
                });
            }
            modes.push(ModalAcousticMode {
                angular_frequency_rad_s: mode.angular_frequency_rad_s,
                damping_ratio,
                // Observer transfers are applied after one shared state step,
                // keeping stereo phase exact without duplicating oscillators.
                pressure_per_modal_velocity: C64::ZERO,
            });
            modal_damping_ratios.push(damping_ratio);
        }
        let runtime = ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?;
        Ok(Self {
            basis,
            runtime,
            modal_damping_ratios,
            sample_rate_hz,
            maximum_abs_pressure_pa: budget.maximum_abs_pressure_pa,
            radiation_identity: radiation.identity,
            damping_model_identity,
        })
    }

    /// Synthesize simultaneous SI-pressure signals at fixed observers from the
    /// equal-and-opposite base contact reaction.
    ///
    /// Mechanics force-time measures are conservatively reconstructed on the
    /// audio clock before the structural state advances once per audio cell.
    /// Modal acceleration is propagated from every surface triangle at its
    /// physical sound-travel delay using a linear fractional-sample kernel.
    /// Only the transverse force admitted by the current DKT bending model is
    /// projected; no in-plane or housing response is invented.
    pub fn synthesize_control_stream_observers(
        &mut self,
        controls: &EulerControlStream<'_>,
        radiation: &BaffledPlateObserverRadiation,
        force_reconstruction: GeneralizedForceReconstructionInput,
        maximum_contact_surface_distance_m: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<PhysicalPressureSignal>, StructuralModalBasisError> {
        if radiation.identity != self.radiation_identity
            || radiation.structural_basis_identity != self.basis.identity
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "plate observer radiation does not match the audio runtime",
            });
        }
        let intervals = controls.audio();
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
        let frame_count = fixed_rate_frame_count_with_roundoff_bound(
            first.start_time_s,
            last.end_time_s,
            self.sample_rate_hz,
            force_reconstruction.clock_roundoff_operation_count,
        )
        .ok_or(StructuralModalBasisError::ControlTimeline {
            what: "control horizon is not an integral number of plate-audio samples within declared clock roundoff",
        })?;
        let mut pressure_channels = vec![Vec::new(); radiation.observers.len()];
        for channel in &mut pressure_channels {
            channel.try_reserve_exact(frame_count).map_err(|_| {
                StructuralModalBasisError::PressureCapacity {
                    requested: frame_count,
                }
            })?;
        }
        let mut modal_forces = Vec::with_capacity(intervals.len());
        for interval in intervals {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
            modal_forces.push(self.modal_force_for_interval(
                controls,
                interval,
                maximum_contact_surface_distance_m,
            )?);
        }
        let reconstructed =
            reconstruct_projected_modal_forces(intervals, &modal_forces, force_reconstruction, cx)?;
        if reconstructed.sample_rate_hz != self.sample_rate_hz
            || reconstructed.frame_count() != frame_count
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "reconstructed plate force clock does not match modal audio clock",
            });
        }
        // This source trajectory is a causal preroll cropped from an already
        // running spin, not a force-application experiment. Begin from the
        // static deflection under its first reconstructed contact load so the output
        // cannot invent a plate impact at the arbitrary retained horizon.
        self.runtime.initialize_static_equilibrium(
            reconstructed
                .frame(0)
                .expect("positive reconstructed frame count has a first force row"),
        )?;
        let sample_period_s = self.runtime.sample_period_s();
        let mut peaks = vec![0.0_f64; radiation.observers.len()];
        let history_length = radiation.maximum_retarded_delay_frames;
        let history_cells = self.basis.modes.len().checked_mul(history_length).ok_or(
            StructuralModalBasisError::PressureCapacity {
                requested: usize::MAX,
            },
        )?;
        let mut modal_acceleration_history = vec![0.0_f64; history_cells];
        let mut history_head = history_length - 1;
        for frame in 0..frame_count {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let applied_force = reconstructed
                .frame(frame)
                .expect("reconstructed frame count matches plate pressure horizon");
            self.runtime.step_duration(applied_force, sample_period_s)?;
            history_head = (history_head + 1) % history_length;
            for (mode_index, ((mode, damping_ratio), state)) in self
                .basis
                .modes
                .iter()
                .zip(&self.modal_damping_ratios)
                .zip(self.runtime.states())
                .enumerate()
            {
                let omega = mode.angular_frequency_rad_s;
                let acceleration = applied_force[mode_index]
                    - 2.0 * damping_ratio * omega * state.velocity_m_sqrt_kg_per_s
                    - omega * omega * state.displacement_m_sqrt_kg;
                if !acceleration.is_finite() {
                    return Err(StructuralModalBasisError::InvalidRequest {
                        what: "plate modal acceleration became non-finite",
                    });
                }
                modal_acceleration_history[mode_index * history_length + history_head] =
                    acceleration;
            }
            for (observer_index, mode_kernels) in radiation
                .pressure_per_modal_acceleration_fir
                .iter()
                .enumerate()
            {
                let mut pressure = 0.0_f64;
                let mut correction = 0.0_f64;
                for (mode_index, kernel) in mode_kernels.iter().enumerate() {
                    let history_offset = mode_index * history_length;
                    for (lag, coefficient) in kernel.iter().copied().enumerate() {
                        let history_index = (history_head + history_length - lag) % history_length;
                        let contribution = coefficient
                            * modal_acceleration_history[history_offset + history_index];
                        let corrected = contribution - correction;
                        let next = pressure + corrected;
                        correction = (next - pressure) - corrected;
                        pressure = next;
                    }
                }
                if !pressure.is_finite() || pressure.abs() > self.maximum_abs_pressure_pa {
                    return Err(StructuralModalBasisError::InvalidRequest {
                        what: "retarded baffled-plate observer pressure exceeded its finite safety envelope",
                    });
                }
                peaks[observer_index] = peaks[observer_index].max(pressure.abs());
                pressure_channels[observer_index].push(pressure);
            }
        }
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        Ok(radiation
            .observers
            .iter()
            .copied()
            .zip(pressure_channels)
            .zip(peaks)
            .map(|((observer, pressure_pa), peak_abs_pressure_pa)| {
                let observer = PhysicalPressureObserver::WorldFixed(observer);
                let identity = baffled_plate_pressure_signal_identity(
                    self,
                    controls,
                    first.start_time_s,
                    observer,
                    reconstructed.identity,
                    &pressure_pa,
                );
                PhysicalPressureSignal {
                    start_time_s: first.start_time_s,
                    sample_rate_hz: self.sample_rate_hz,
                    pressure_pa,
                    peak_abs_pressure_pa,
                    contact_force_sampling:
                        PhysicalContactForceSampling::IntervalMeasureAtClosingElseOpeningEndpointBandLimitedV1,
                    observer,
                    structural_basis_identity: self.basis.identity,
                    radiation_identity: self.radiation_identity,
                    damping_model_identity: self.damping_model_identity,
                    identity,
                }
            })
            .collect())
    }

    /// Apply the fixed-plate retarded radiation FIR directly to modal
    /// acceleration from the accepted coupled-mechanics state.
    pub fn synthesize_decimated_acceleration_observers(
        &self,
        acceleration: &DecimatedModalAcceleration,
        radiation: &BaffledPlateObserverRadiation,
        cx: &Cx<'_>,
    ) -> Result<Vec<PhysicalPressureSignal>, StructuralModalBasisError> {
        if radiation.identity != self.radiation_identity
            || radiation.structural_basis_identity != self.basis.identity
            || acceleration.plate_model_identity
                != [self.basis.identity, self.damping_model_identity]
            || acceleration.sample_rate_hz != self.sample_rate_hz
            || acceleration.coordinate_count() != self.basis.modes.len()
            || acceleration.frame_count() == 0
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "decimated mechanics acceleration does not match plate radiation",
            });
        }
        let frame_count = acceleration.frame_count();
        let history_length = radiation.maximum_retarded_delay_frames;
        let history_cells = self.basis.modes.len().checked_mul(history_length).ok_or(
            StructuralModalBasisError::PressureCapacity {
                requested: usize::MAX,
            },
        )?;
        let mut history = vec![0.0_f64; history_cells];
        let mut history_head = history_length - 1;
        let mut channels = vec![Vec::new(); radiation.observers.len()];
        let mut peaks = vec![0.0_f64; radiation.observers.len()];
        for channel in &mut channels {
            channel.try_reserve_exact(frame_count).map_err(|_| {
                StructuralModalBasisError::PressureCapacity {
                    requested: frame_count,
                }
            })?;
        }
        for frame in 0..frame_count {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let row = acceleration
                .frame(frame)
                .expect("admitted acceleration frame");
            history_head = (history_head + 1) % history_length;
            for (mode, value) in row.iter().copied().enumerate() {
                history[mode * history_length + history_head] = value;
            }
            for (observer, mode_kernels) in radiation
                .pressure_per_modal_acceleration_fir
                .iter()
                .enumerate()
            {
                let mut pressure = 0.0_f64;
                let mut correction = 0.0_f64;
                for (mode, kernel) in mode_kernels.iter().enumerate() {
                    for (lag, coefficient) in kernel.iter().copied().enumerate() {
                        let index = (history_head + history_length - lag) % history_length;
                        let value =
                            coefficient * history[mode * history_length + index] - correction;
                        let next = pressure + value;
                        correction = (next - pressure) - value;
                        pressure = next;
                    }
                }
                if !pressure.is_finite() || pressure.abs() > self.maximum_abs_pressure_pa {
                    return Err(StructuralModalBasisError::InvalidRequest {
                        what: "retarded plate pressure exceeded its finite safety envelope",
                    });
                }
                peaks[observer] = peaks[observer].max(pressure.abs());
                channels[observer].push(pressure);
            }
        }
        Ok(radiation
            .observers
            .iter()
            .copied()
            .zip(channels)
            .zip(peaks)
            .map(|((observer, pressure_pa), peak_abs_pressure_pa)| {
                let observer = PhysicalPressureObserver::WorldFixed(observer);
                let identity = baffled_plate_decimated_pressure_signal_identity(
                    self,
                    acceleration,
                    observer,
                    &pressure_pa,
                );
                PhysicalPressureSignal {
                    start_time_s: acceleration.start_time_s,
                    sample_rate_hz: self.sample_rate_hz,
                    pressure_pa,
                    peak_abs_pressure_pa,
                    contact_force_sampling:
                        PhysicalContactForceSampling::MechanicsModalAccelerationAntiAliasedDecimatedV1,
                    observer,
                    structural_basis_identity: self.basis.identity,
                    radiation_identity: self.radiation_identity,
                    damping_model_identity: self.damping_model_identity,
                    identity,
                }
            })
            .collect())
    }

    fn modal_force_for_interval(
        &self,
        controls: &EulerControlStream<'_>,
        interval: &crate::AudioControlInterval,
        maximum_contact_surface_distance_m: f64,
    ) -> Result<Vec<f64>, StructuralModalBasisError> {
        let visualization = controls.visualization();
        let end = visualization
            .get(interval.visual_coverage.end_visualization_index)
            .ok_or(StructuralModalBasisError::ControlTimeline {
                what: "plate audio closing visualization index is out of bounds",
            })?;
        let start = interval
            .visual_coverage
            .start_visualization_index
            .and_then(|index| visualization.get(index));
        let selected = if end.contact.is_some() {
            Some(end)
        } else {
            start.filter(|point| point.contact.is_some())
        };
        let Some(selected) = selected else {
            if interval.interval_contact_active {
                return Err(StructuralModalBasisError::MissingContactLocation {
                    source_sample: interval.source_sample_index,
                });
            }
            return Ok(vec![0.0; self.basis.modes.len()]);
        };
        let contact = selected.contact.expect("selected point has contact");
        let transverse_force_on_plate_n = match interval.channels.contact {
            ChannelControl::Available(contact_control) => {
                let force_on_disc_base = selected
                    .orientation_base_to_world
                    .rotate_world_to_body(contact_control.mean_force_world_n);
                -force_on_disc_base.z
            }
            ChannelControl::Unavailable => {
                if let Some(normal_force_n) = interval.mean_base_normal_contact_force_n {
                    -normal_force_n * contact.normal_base.z
                } else if interval.interval_contact_active {
                    return Err(StructuralModalBasisError::MissingContactForce {
                        source_sample: interval.source_sample_index,
                    });
                } else {
                    0.0
                }
            }
        };
        if transverse_force_on_plate_n == 0.0 {
            return Ok(vec![0.0; self.basis.modes.len()]);
        }
        let plate_reference_point_m = plate_reference_contact_point(
            [
                contact.point_base_m.x,
                contact.point_base_m.y,
                contact.point_base_m.z,
            ],
            selected.signed_gap_m,
        )?;
        self.basis
            .project_transverse_point_force(
                plate_reference_point_m,
                transverse_force_on_plate_n,
                maximum_contact_surface_distance_m,
            )
            .map(|projection| projection.modal_force_n_per_sqrt_kg)
    }
}

fn plate_reference_contact_point(
    disc_contact_point_base_m: [f64; 3],
    signed_gap_m: f64,
) -> Result<[f64; 3], StructuralModalBasisError> {
    if disc_contact_point_base_m
        .iter()
        .any(|coordinate| !coordinate.is_finite())
        || !signed_gap_m.is_finite()
    {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "disc-side base-frame contact point or retained signed gap is non-finite",
        });
    }

    // The retained point is the actual disc material point expressed relative
    // to the reduced base-mode origin. Its z coordinate equals the smooth
    // profile gap, not necessarily `signed_gap_m`: the production contact law
    // may additionally resolve a local surface-height field that is not a
    // rigid displacement of the whole plate. Plate modes are defined on the
    // undeformed z=0 reference midsurface, so only the material point's exact
    // in-plane coordinates are projected onto that surface. The resolved gap
    // is independently retained for unilateral branch/contact validation.
    Ok([
        disc_contact_point_base_m[0],
        disc_contact_point_base_m[1],
        0.0,
    ])
}

fn baffled_plate_pressure_signal_identity(
    model: &BaffledPlateModalAudioModel<'_>,
    controls: &EulerControlStream<'_>,
    start_time_s: f64,
    observer: PhysicalPressureObserver,
    force_reconstruction_identity: ContentHash,
    pressure_pa: &[f64],
) -> ContentHash {
    let mut hasher =
        DomainHasher::new("org.frankensim.euler-disc.baffled-plate-pressure-signal.v2");
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
    hasher.update(force_reconstruction_identity.as_bytes());
    hash_pressure_observer(&mut hasher, observer);
    for pressure in pressure_pa {
        hasher.update(&pressure.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

fn baffled_plate_decimated_pressure_signal_identity(
    model: &BaffledPlateModalAudioModel<'_>,
    acceleration: &DecimatedModalAcceleration,
    observer: PhysicalPressureObserver,
    pressure_pa: &[f64],
) -> ContentHash {
    let mut hasher =
        DomainHasher::new("org.frankensim.euler-disc.baffled-plate-pressure-signal.v3");
    hasher.update(model.basis.identity.as_bytes());
    hasher.update(model.radiation_identity.as_bytes());
    hasher.update(model.damping_model_identity.as_bytes());
    hasher.update(acceleration.identity.as_bytes());
    hasher.update(&acceleration.start_time_s.to_bits().to_le_bytes());
    hash_pressure_observer(&mut hasher, observer);
    for pressure in pressure_pa {
        hasher.update(&pressure.to_bits().to_le_bytes());
    }
    hasher.finalize()
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
        let mut static_pressure_transfers = Vec::with_capacity(basis.modes.len());
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
            static_pressure_transfers.push(acoustic.observer_pressure_per_modal_velocity);
        }
        let runtime = ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?;
        Ok(Self {
            basis,
            runtime,
            static_pressure_transfers,
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
            static_pressure_transfers: vec![C64::ZERO; basis.modes.len()],
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
    /// Mechanics force-time measures are conservatively rasterized onto the
    /// exact audio cells and band-limited before modal integration. When a full
    /// contact wrench is unavailable, an explicitly admitted interval-mean
    /// normal reaction is reconstructed along the retained contact normal. No
    /// impulse is invented from a timing-only contact event.
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
        initial_state: PhysicalModalInitialState,
        force_reconstruction: GeneralizedForceReconstructionInput,
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
        let frame_count = fixed_rate_frame_count_with_roundoff_bound(
            first.start_time_s,
            last.end_time_s,
            self.sample_rate_hz,
            force_reconstruction.clock_roundoff_operation_count,
        )
        .ok_or(StructuralModalBasisError::ControlTimeline {
            what: "control horizon is not a positive integral number of audio samples within declared clock roundoff",
        })?;
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
        let reconstructed =
            reconstruct_projected_modal_forces(intervals, &modal_forces, force_reconstruction, cx)?;
        if reconstructed.sample_rate_hz != self.sample_rate_hz
            || reconstructed.frame_count() != frame_count
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "reconstructed specimen force clock does not match modal audio clock",
            });
        }
        if initial_state == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce {
            self.runtime.initialize_static_equilibrium(
                reconstructed
                    .frame(0)
                    .expect("positive reconstructed frame count has a first force row"),
            )?;
        }

        let sample_period_s = self.runtime.sample_period_s();
        let mut peak_abs_pressure_pa = 0.0_f64;
        for frame in 0..frame_count {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let applied_force = reconstructed
                .frame(frame)
                .expect("reconstructed frame count matches specimen pressure horizon");
            self.runtime.step_duration(applied_force, sample_period_s)?;
            let final_pressure =
                if initial_state == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce {
                    self.runtime
                        .observer_pressure_with_transfers_about_static_equilibrium(
                            &self.static_pressure_transfers,
                            applied_force,
                        )?
                } else {
                    self.runtime
                        .observer_pressure_with_transfers(&self.static_pressure_transfers)?
                };
            peak_abs_pressure_pa = peak_abs_pressure_pa.max(final_pressure.abs());
            pressure_pa.push(final_pressure);
        }
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let observer = PhysicalPressureObserver::BodyFixed(observer);
        let identity = physical_pressure_signal_identity(
            self,
            controls,
            first.start_time_s,
            observer,
            reconstructed.identity,
            &pressure_pa,
        );
        Ok(PhysicalPressureSignal {
            start_time_s: first.start_time_s,
            sample_rate_hz: self.sample_rate_hz,
            pressure_pa,
            peak_abs_pressure_pa,
            contact_force_sampling:
                PhysicalContactForceSampling::IntervalMeasureAtClosingElseOpeningEndpointBandLimitedV1,
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
        initial_state: PhysicalModalInitialState,
        force_reconstruction: GeneralizedForceReconstructionInput,
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
        let frame_count = fixed_rate_frame_count_with_roundoff_bound(
            first.start_time_s,
            last.end_time_s,
            self.sample_rate_hz,
            force_reconstruction.clock_roundoff_operation_count,
        )
        .ok_or(StructuralModalBasisError::ControlTimeline {
            what: "control horizon is not a positive integral number of audio samples within declared clock roundoff",
        })?;
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
        let reconstructed =
            reconstruct_projected_modal_forces(intervals, &modal_forces, force_reconstruction, cx)?;
        if reconstructed.sample_rate_hz != self.sample_rate_hz
            || reconstructed.frame_count() != frame_count
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "reconstructed specimen force clock does not match world-observer audio clock",
            });
        }
        if initial_state == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce {
            self.runtime.initialize_static_equilibrium(
                reconstructed
                    .frame(0)
                    .expect("positive reconstructed frame count has a first force row"),
            )?;
        }

        let timeline = TimelineResampler::new(controls.source());
        let mut transfer_scratch = vec![Vec::new(); observers.len()];
        let mut peaks = vec![0.0_f64; observers.len()];
        let sample_period_s = self.runtime.sample_period_s();
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
            let applied_force = reconstructed
                .frame(frame)
                .expect("reconstructed frame count matches world-observer pressure horizon");
            self.runtime.step_duration(applied_force, sample_period_s)?;
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
                let pressure = if initial_state
                    == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce
                {
                    self.runtime
                        .observer_pressure_with_transfers_about_static_equilibrium(
                            &transfer_scratch[observer_index],
                            applied_force,
                        )?
                } else {
                    self.runtime
                        .observer_pressure_with_transfers(&transfer_scratch[observer_index])?
                };
                peaks[observer_index] = peaks[observer_index].max(pressure.abs());
                pressure_channels[observer_index].push(pressure);
            }
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
                    reconstructed.identity,
                    &pressure_pa,
                );
                PhysicalPressureSignal {
                    start_time_s: first.start_time_s,
                    sample_rate_hz: self.sample_rate_hz,
                    pressure_pa,
                    peak_abs_pressure_pa,
                    contact_force_sampling:
                        PhysicalContactForceSampling::IntervalMeasureAtClosingElseOpeningEndpointBandLimitedV1,
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
        modal_force_for_control_interval(
            &self.basis.mesh,
            &self.basis.modes,
            controls,
            interval,
            maximum_contact_distance_m,
        )
    }
}

fn modal_force_for_control_interval(
    mesh: &RoundedCylinderTetMesh,
    modes: &[StructuralMode],
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
    } else if let Some((start, contact)) = start.and_then(|start| start.contact.map(|c| (start, c)))
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
            if let (Some(force), Some(normal)) =
                (interval.mean_base_normal_contact_force_n, normal_world)
            {
                normal.scale(force)
            } else if interval.interval_contact_active {
                return Err(StructuralModalBasisError::MissingContactForce {
                    source_sample: interval.source_sample_index,
                });
            } else {
                Vec3::ZERO
            }
        }
    };
    if force_world.norm_squared() == 0.0 {
        return Ok(vec![0.0; modes.len()]);
    }
    let point = point.ok_or(StructuralModalBasisError::MissingContactLocation {
        source_sample: interval.source_sample_index,
    })?;
    let force_body = orientation.rotate_world_to_body(force_world);
    project_point_force_on_modes(
        mesh,
        modes,
        [point.x, point.y, point.z],
        [force_body.x, force_body.y, force_body.z],
        maximum_contact_distance_m,
    )
    .map(|projection| projection.modal_force_n_per_sqrt_kg)
}

fn hash_pressure_observer(hasher: &mut DomainHasher, observer: PhysicalPressureObserver) {
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
}

fn physical_pressure_signal_identity(
    model: &PhysicalModalAudioModel<'_>,
    controls: &EulerControlStream<'_>,
    start_time_s: f64,
    observer: PhysicalPressureObserver,
    force_reconstruction_identity: ContentHash,
    pressure_pa: &[f64],
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-disc.physical-pressure-signal.v2");
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
    hasher.update(force_reconstruction_identity.as_bytes());
    hash_pressure_observer(&mut hasher, observer);
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
        project_point_force_on_modes(
            &self.mesh,
            &self.modes,
            point_m,
            force_n,
            maximum_distance_m,
        )
    }
}

impl StructuralResidualFlexibilityEstimateBasis {
    /// Recompute the content identity from every stored semantic field.
    #[must_use]
    pub fn recomputed_identity(&self) -> ContentHash {
        residual_flexibility_basis_identity(self)
    }

    /// Evaluate truncated static flexibility at any admitted moving boundary point.
    ///
    /// # Errors
    /// Refuses malformed/off-boundary force inputs, non-finite reconstruction,
    /// or a negative compliance-work defect beyond roundoff.
    pub fn evaluate_point_force(
        &self,
        point_m: [f64; 3],
        force_n: [f64; 3],
        maximum_distance_m: f64,
    ) -> Result<StructuralResidualFlexibilityEstimateResponse, StructuralModalBasisError> {
        let (force_projection, inertia_relieved_nodal_force_n, rigid_force_residual) =
            inertia_relieved_point_force(self, point_m, force_n, maximum_distance_m)?;
        let mut modal_displacement_m_sqrt_kg = Vec::with_capacity(self.enrichment_modes.len());
        let mut nodal_displacement_m = vec![[0.0; 3]; self.mesh.nodes_m.len()];
        let mut panel_normal_displacement_m = vec![0.0; self.mesh.boundary.triangles.len()];
        for (mode, generalized_force) in self
            .enrichment_modes
            .iter()
            .zip(&force_projection.modal_force_n_per_sqrt_kg)
        {
            let coordinate = generalized_force / mode.eigenvalue_s2;
            if !coordinate.is_finite() {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "residual-flexibility modal displacement became non-finite",
                });
            }
            modal_displacement_m_sqrt_kg.push(coordinate);
            for (displacement, shape) in nodal_displacement_m
                .iter_mut()
                .zip(&mode.nodal_shape_per_sqrt_kg)
            {
                for component in 0..3 {
                    displacement[component] =
                        coordinate.mul_add(shape[component], displacement[component]);
                }
            }
            for (displacement, shape) in panel_normal_displacement_m
                .iter_mut()
                .zip(&mode.panel_normal_shape_per_sqrt_kg)
            {
                *displacement = coordinate.mul_add(*shape, *displacement);
            }
        }
        if nodal_displacement_m
            .iter()
            .flatten()
            .chain(&panel_normal_displacement_m)
            .any(|value| !value.is_finite())
        {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "residual-flexibility reconstruction became non-finite",
            });
        }
        let triangle = self.mesh.boundary.triangles[force_projection.boundary_triangle];
        let contact_displacement_m = core::array::from_fn(|component| {
            (0..3).fold(0.0, |value, corner| {
                force_projection.barycentric[corner]
                    .mul_add(nodal_displacement_m[triangle[corner]][component], value)
            })
        });
        let elastic_work_j = dot(contact_displacement_m, force_n);
        let reduced_displacement: Vec<f64> = self
            .assembly
            .free_dofs
            .iter()
            .map(|dof| nodal_displacement_m[dof / 3][dof % 3])
            .collect();
        let mut stiffness_displacement = vec![0.0; reduced_displacement.len()];
        self.assembly
            .stiffness
            .spmv(&reduced_displacement, &mut stiffness_displacement);
        let recoverable_strain_energy_j =
            0.5 * dot_slice(&reduced_displacement, &stiffness_displacement);
        let energy_closure_residual_j = elastic_work_j - 2.0 * recoverable_strain_energy_j;
        let energy_scale_j = elastic_work_j
            .abs()
            .max((2.0 * recoverable_strain_energy_j).abs())
            .max(f64::MIN_POSITIVE);
        let work_tolerance_j = 1.0e-8 * energy_scale_j;
        verify_residual(
            "non-negative residual-flexibility work [J]",
            (-elastic_work_j.min(recoverable_strain_energy_j)).max(0.0),
            work_tolerance_j,
        )?;
        verify_residual(
            "contact-work/strain-energy closure [J]",
            energy_closure_residual_j.abs(),
            work_tolerance_j,
        )?;
        let mut response = StructuralResidualFlexibilityEstimateResponse {
            basis_identity: self.identity,
            authority: self.authority,
            requested_point_m: point_m,
            applied_force_n: force_n,
            force_projection,
            inertia_relieved_nodal_force_n,
            maximum_rigid_force_relative_residual: rigid_force_residual,
            modal_displacement_m_sqrt_kg,
            nodal_displacement_m,
            panel_normal_displacement_m,
            elastic_work_j,
            recoverable_strain_energy_j,
            energy_closure_residual_j,
            identity: ContentHash([0; 32]),
        };
        response.identity = response.recomputed_identity();
        Ok(response)
    }
}

impl StructuralResidualFlexibilityEstimateResponse {
    /// Recompute the content identity from every stored semantic field.
    #[must_use]
    pub fn recomputed_identity(&self) -> ContentHash {
        residual_flexibility_response_identity(self)
    }
}

/// Compare one physical point-force response across nested enrichment bands.
///
/// Caller owns the QoI tolerance; values bind only this force/point/mesh/cutoff pair.
///
/// # Errors
/// Refuses non-nested or physically mismatched bases and every point-force
/// evaluation error.
pub fn compare_structural_residual_flexibility_estimates(
    coarse: &StructuralResidualFlexibilityEstimateBasis,
    fine: &StructuralResidualFlexibilityEstimateBasis,
    point_m: [f64; 3],
    force_n: [f64; 3],
    maximum_distance_m: f64,
) -> Result<StructuralResidualFlexibilityEstimateComparison, StructuralModalBasisError> {
    if coarse.specimen_identity != fine.specimen_identity
        || coarse.profile_identity != fine.profile_identity
        || coarse.material_state_identity != fine.material_state_identity
        || coarse.operator_identity != fine.operator_identity
        || !pair_bits_equal(
            coarse.forcing_frequency_band_hz,
            fine.forcing_frequency_band_hz,
        )
        || coarse.enrichment_frequency_band_hz.0.to_bits()
            != fine.enrichment_frequency_band_hz.0.to_bits()
        || fine.enrichment_frequency_band_hz.1 <= coarse.enrichment_frequency_band_hz.1
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "residual-flexibility comparison needs nested bands on identical physical operators",
        });
    }
    let coarse_response = coarse.evaluate_point_force(point_m, force_n, maximum_distance_m)?;
    let fine_response = fine.evaluate_point_force(point_m, force_n, maximum_distance_m)?;
    let elastic_work_increment_j = fine_response.elastic_work_j - coarse_response.elastic_work_j;
    let work_scale = fine_response
        .elastic_work_j
        .abs()
        .max(coarse_response.elastic_work_j.abs())
        .max(f64::MIN_POSITIVE);
    let relative_elastic_work_difference = elastic_work_increment_j.abs() / work_scale;
    let relative_nodal_displacement_l2_difference = relative_vec3_l2_difference(
        &coarse_response.nodal_displacement_m,
        &fine_response.nodal_displacement_m,
    );
    let relative_panel_normal_l2_difference = relative_l2_difference(
        &coarse_response.panel_normal_displacement_m,
        &fine_response.panel_normal_displacement_m,
    );
    Ok(StructuralResidualFlexibilityEstimateComparison {
        coarse_response_identity: coarse_response.identity,
        fine_response_identity: fine_response.identity,
        elastic_work_increment_j,
        relative_elastic_work_difference,
        relative_nodal_displacement_l2_difference,
        relative_panel_normal_l2_difference,
    })
}

fn inertia_relieved_point_force(
    basis: &StructuralResidualFlexibilityEstimateBasis,
    point_m: [f64; 3],
    force_n: [f64; 3],
    maximum_distance_m: f64,
) -> Result<(PointForceProjection, Vec<[f64; 3]>, f64), StructuralModalBasisError> {
    let mut projection = project_point_force_on_modes(
        &basis.mesh,
        &basis.enrichment_modes,
        point_m,
        force_n,
        maximum_distance_m,
    )?;
    let full_dof_count = basis.mesh.nodes_m.len() * 3;
    let mut reduced_by_full = vec![usize::MAX; full_dof_count];
    for (reduced, &full) in basis.assembly.free_dofs.iter().enumerate() {
        reduced_by_full[full] = reduced;
    }
    let mut relieved = vec![0.0; basis.assembly.free_dofs.len()];
    let triangle = basis.mesh.boundary.triangles[projection.boundary_triangle];
    for corner in 0..3 {
        for component in 0..3 {
            let reduced = reduced_by_full[3 * triangle[corner] + component];
            if reduced == usize::MAX {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "residual flexibility requires a free-body nodal load",
                });
            }
            relieved[reduced] += projection.barycentric[corner] * force_n[component];
        }
    }
    let rigid_loads: Vec<f64> = basis
        .rigid_modes_per_sqrt_kg
        .iter()
        .map(|mode| dot_slice(mode, &relieved))
        .collect();
    for (mode, generalized_load) in basis.rigid_modes_per_sqrt_kg.iter().zip(&rigid_loads) {
        let mut inertial_load = vec![0.0; relieved.len()];
        basis.assembly.mass.spmv(mode, &mut inertial_load);
        for (load, inertial) in relieved.iter_mut().zip(inertial_load) {
            *load = generalized_load.mul_add(-inertial, *load);
        }
    }
    let rigid_scale = maximum_abs(&rigid_loads).max(f64::MIN_POSITIVE);
    let rigid_residual = basis
        .rigid_modes_per_sqrt_kg
        .iter()
        .map(|mode| dot_slice(mode, &relieved).abs() / rigid_scale)
        .fold(0.0_f64, f64::max);
    verify_residual(
        "inertia-relieved rigid generalized force",
        rigid_residual,
        1.0e-10,
    )?;
    // Overwrite the bare geometric projection with the inertia-relieved
    // generalized force: downstream residual-flexibility displacements and
    // response identity consume the elastic forcing with rigid-body
    // equilibrium removed. This field therefore intentionally differs from a
    // fresh `project_point_force_on_modes` result by the relief residue.
    projection.modal_force_n_per_sqrt_kg = basis
        .enrichment_modes
        .iter()
        .map(|mode| {
            basis
                .assembly
                .free_dofs
                .iter()
                .zip(&relieved)
                .map(|(dof, load)| mode.nodal_shape_per_sqrt_kg[dof / 3][dof % 3] * load)
                .sum()
        })
        .collect();
    let mut full_relieved = vec![[0.0; 3]; basis.mesh.nodes_m.len()];
    for (&full, load) in basis.assembly.free_dofs.iter().zip(relieved) {
        full_relieved[full / 3][full % 3] = load;
    }
    Ok((projection, full_relieved, rigid_residual))
}

fn project_point_force_on_modes(
    mesh: &RoundedCylinderTetMesh,
    modes: &[StructuralMode],
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
    for (boundary_triangle, triangle) in mesh.boundary.triangles.iter().enumerate() {
        let vertices = triangle.map(|node| mesh.nodes_m[node]);
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
    let triangle = mesh.boundary.triangles[boundary_triangle];
    let modal_force_n_per_sqrt_kg = modes
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

/// Build a causal broadband body-frame source bank from the exact residual
/// basis, material loss evidence, gas state, and disjoint BEM grids.
///
/// At drive frequency `omega`, unit generalized modal acceleration implies
/// `q = -1/omega^2` and therefore `v_n = +i phi_n/omega` under
/// `exp(-i omega t)`. One shared BEM factorization serves every modal input at
/// each frequency. Static residual flexibility is never an input to this path.
pub fn build_structural_broadband_radiation_artifact(
    request: &StructuralBroadbandRadiationRequest<'_>,
    cx: &Cx<'_>,
) -> Result<StructuralBroadbandRadiationArtifact, StructuralModalBasisError> {
    let basis = request.basis;
    let loss = request.loss;
    let sample_rate_hz = request.fit.sample_rate_hz as u32;
    let forcing = basis.forcing_frequency_band_hz;
    let enrichment = basis.enrichment_frequency_band_hz;
    let grids_valid = |grid: &[f64]| {
        !grid.is_empty()
            && grid.len() <= MAX_BROADBAND_FREQUENCIES
            && grid.iter().enumerate().all(|(index, &frequency)| {
                frequency > 0.0
                    && frequency <= enrichment.1
                    && frequency < 0.5 * request.fit.sample_rate_hz
                    && frequency.is_finite()
                    && (index == 0 || grid[index - 1] < frequency)
            })
    };
    if basis.schema_version != STRUCTURAL_RESIDUAL_FLEXIBILITY_SCHEMA_VERSION
        || basis.authority != StructuralResidualFlexibilityAuthority::EstimateOnly
        || basis.recomputed_identity() != basis.identity
        || basis.operator_identity == ContentHash([0; 32])
        || basis.enrichment_modes.is_empty()
        || basis.enrichment_modes.len() > MAX_RADIATION_FIELDS_PER_BATCH
        || basis.certified_enrichment_mode_count != basis.enrichment_modes.len()
        || basis.certified_partition_mode_count
            != basis
                .certified_in_band_mode_count
                .checked_add(basis.certified_enrichment_mode_count)
                .unwrap_or(usize::MAX)
        || forcing.1.to_bits() != enrichment.0.to_bits()
        || !(forcing.0 > 0.0 && forcing.0 < forcing.1 && enrichment.0 < enrichment.1)
        || sample_rate_hz == 0
        || f64::from(sample_rate_hz).to_bits() != request.fit.sample_rate_hz.to_bits()
        || !grids_valid(request.training_frequency_hz)
        || !grids_valid(request.held_out_frequency_hz)
        || request.held_out_frequency_hz.iter().any(|held| {
            request
                .training_frequency_hz
                .iter()
                .any(|training| training.to_bits() == held.to_bits())
        })
        || !(8..=MAX_VALIDATION_DIRECTIONS).contains(&request.held_out_directions_body.len())
        || request.held_out_directions_body.iter().any(|direction| {
            let norm = norm_squared(*direction);
            !(norm > 0.0 && norm.is_finite())
        })
        || request.directivity.maximum_spherical_harmonic_degree > MAX_SH_DEGREE
        || !(request.directivity.minimum_captured_fraction > 0.0
            && request.directivity.minimum_captured_fraction <= 1.0
            && request.directivity.minimum_captured_fraction.is_finite())
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "broadband basis partition, grids, directions, SH controls, or sample clock are inconsistent",
        });
    }
    validate_acoustic_medium(request.medium.gas)?;
    if request.medium.gas_model_identity == ContentHash([0; 32])
        || loss.structural_basis_identity != basis.identity
        || loss.material_state_identity != basis.material_state_identity
        || loss.damping_model_identity == ContentHash([0; 32])
        || loss.loss_factors.len() != basis.enrichment_modes.len()
        || loss
            .loss_factors
            .iter()
            .any(|value| !(value.is_finite() && *value >= 0.0))
        || basis.enrichment_modes.iter().any(|mode| {
            mode.panel_normal_shape_per_sqrt_kg.len() != basis.mesh.boundary.triangles.len()
                || !(mode.angular_frequency_rad_s > 0.0
                    && mode.angular_frequency_rad_s.is_finite()
                    && mode.eigenvalue_interval_s2.0 > (core::f64::consts::TAU * forcing.1).powi(2)
                    && mode.frequency_hz <= enrichment.1
                    && mode.frequency_hz < 0.5 * request.fit.sample_rate_hz)
        })
    {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "broadband damping, medium, or enrichment modes do not match the exact residual basis",
        });
    }

    let surface = SpherePanels::from_triangles(
        basis
            .mesh
            .boundary
            .triangles
            .iter()
            .map(|triangle| triangle.map(|node| basis.mesh.nodes_m[node]))
            .collect(),
    )?;
    let medium = Medium {
        density: request.medium.gas.density,
        sound_speed: request.medium.gas.sound_speed,
    };
    let radius_m = basis
        .mesh
        .nodes_m
        .iter()
        .map(|point| norm_squared(*point).sqrt())
        .fold(0.0_f64, f64::max);
    let mut training = Vec::with_capacity(request.training_frequency_hz.len());
    let mut held_out = Vec::with_capacity(request.held_out_frequency_hz.len());
    for &frequency_hz in request.training_frequency_hz {
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        let (_, solutions) = solve_modal_acceleration_radiation_batch(
            &surface,
            &basis.enrichment_modes,
            omega,
            medium,
            radius_m,
        )?;
        let tables = checked_directivity_tables(
            &surface,
            &solutions,
            &basis.enrichment_modes,
            medium,
            request.directivity.maximum_spherical_harmonic_degree,
            request.directivity.minimum_captured_fraction,
            "training",
            frequency_hz,
        )?;
        let captured_fraction = captured_fraction(&solutions, &tables);
        training.push(ComplexShTrainingSample {
            omega_rad_s: omega,
            coefficients_by_input: tables.into_iter().map(|table| table.coefficients).collect(),
            diagnostics: radiation_sample_diagnostics(&solutions, captured_fraction),
        });
    }
    for &frequency_hz in request.held_out_frequency_hz {
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        let (_, solutions) = solve_modal_acceleration_radiation_batch(
            &surface,
            &basis.enrichment_modes,
            omega,
            medium,
            radius_m,
        )?;
        let tables = checked_directivity_tables(
            &surface,
            &solutions,
            &basis.enrichment_modes,
            medium,
            request.directivity.maximum_spherical_harmonic_degree,
            request.directivity.minimum_captured_fraction,
            "held-out",
            frequency_hz,
        )?;
        let captured_fraction = captured_fraction(&solutions, &tables);
        held_out.push(DirectFarFieldHeldOutSample {
            omega_rad_s: omega,
            directions: request.held_out_directions_body.to_vec(),
            far_field_by_input: solutions
                .iter()
                .map(|solution| {
                    far_field(&surface, solution, medium, request.held_out_directions_body)
                })
                .collect(),
            diagnostics: radiation_sample_diagnostics(&solutions, captured_fraction),
        });
    }
    let input_ids: Vec<String> = (0..basis.enrichment_modes.len())
        .map(|mode| format!("enrichment-mode-{mode:04}-qddot[m*sqrt(kg)/s^2]"))
        .collect();
    let source_identity = broadband_sample_identity(request);
    let samples = SampledRadiationData {
        source_id: format!("euler-structural-modal-acceleration-pa-m-v1:{source_identity}"),
        harmonic_time_convention: HarmonicTimeConvention::ExpNegativeIOmegaT,
        l_max: request.directivity.maximum_spherical_harmonic_degree,
        input_ids,
        training,
        held_out,
    };
    cx.checkpoint()
        .map_err(|_| StructuralModalBasisError::Cancelled)?;
    let radiation = build_broadband_radiation_artifact(&samples, request.fit)?;
    let source_diameter_m = 2.0 * radius_m;
    let identity = structural_broadband_artifact_identity(
        source_identity,
        basis.identity,
        request.medium.gas_model_identity,
        request.medium.gas.sound_speed,
        source_diameter_m,
        loss.damping_model_identity,
        &loss.loss_factors,
        &radiation,
    );
    Ok(StructuralBroadbandRadiationArtifact {
        identity,
        authority: BroadbandRadiationAuthority::EstimateOnly,
        structural_basis_identity: basis.identity,
        gas_model_identity: request.medium.gas_model_identity,
        sound_speed_m_s: request.medium.gas.sound_speed,
        source_diameter_m,
        damping_model_identity: loss.damping_model_identity,
        sample_rate_hz,
        modal_loss_factors: loss.loss_factors.clone(),
        radiation,
        no_claims: STRUCTURAL_BROADBAND_SOURCE_NO_CLAIM,
    })
}

/// Build the compact rigid-disc radiation bank used for the audible wobble
/// source. Each input is one unit body-frame generalized acceleration; under
/// `exp(-i omega t)` its boundary velocity is `+i shape/omega`.
pub fn build_rigid_disc_broadband_radiation_artifact(
    request: &RigidDiscBroadbandRadiationRequest<'_>,
    cx: &Cx<'_>,
) -> Result<RigidDiscBroadbandRadiationArtifact, StructuralModalBasisError> {
    let basis = request.basis;
    let sample_rate_hz = request.fit.sample_rate_hz as u32;
    let grids_valid = |grid: &[f64]| {
        !grid.is_empty()
            && grid.len() <= MAX_BROADBAND_FREQUENCIES
            && grid.iter().enumerate().all(|(index, &frequency)| {
                frequency > 0.0
                    && frequency < 0.5 * request.fit.sample_rate_hz
                    && frequency.is_finite()
                    && (index == 0 || grid[index - 1] < frequency)
            })
    };
    if basis.schema_version != STRUCTURAL_RESIDUAL_FLEXIBILITY_SCHEMA_VERSION
        || basis.recomputed_identity() != basis.identity
        || basis.mesh.boundary.triangles.is_empty()
        || sample_rate_hz == 0
        || f64::from(sample_rate_hz).to_bits() != request.fit.sample_rate_hz.to_bits()
        || !grids_valid(request.training_frequency_hz)
        || !grids_valid(request.held_out_frequency_hz)
        || request.held_out_frequency_hz.iter().any(|held| {
            request
                .training_frequency_hz
                .iter()
                .any(|training| training.to_bits() == held.to_bits())
        })
        || !(8..=MAX_VALIDATION_DIRECTIONS).contains(&request.held_out_directions_body.len())
        || request.held_out_directions_body.iter().any(|direction| {
            let norm = norm_squared(*direction);
            !(norm > 0.0 && norm.is_finite())
        })
        || request.directivity.maximum_spherical_harmonic_degree > MAX_SH_DEGREE
        || !(request.directivity.minimum_captured_fraction > 0.0
            && request.directivity.minimum_captured_fraction <= 1.0
            && request.directivity.minimum_captured_fraction.is_finite())
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "rigid-disc broadband mesh, grids, directions, SH controls, or sample clock are inconsistent",
        });
    }
    validate_acoustic_medium(request.medium.gas)?;
    if request.medium.gas_model_identity == ContentHash([0; 32]) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "rigid-disc broadband gas identity must not be zero",
        });
    }

    let surface = SpherePanels::from_triangles(
        basis
            .mesh
            .boundary
            .triangles
            .iter()
            .map(|triangle| triangle.map(|node| basis.mesh.nodes_m[node]))
            .collect(),
    )?;
    let radius_m = basis
        .mesh
        .nodes_m
        .iter()
        .map(|point| norm_squared(*point).sqrt())
        .fold(0.0_f64, f64::max);
    let mut shapes =
        vec![vec![0.0; surface.centroids().len()]; RIGID_DISC_ACOUSTIC_COORDINATES.len()];
    for (panel, &normal) in surface.normals().iter().enumerate() {
        shapes[0][panel] = normal[0];
        shapes[1][panel] = normal[1];
        shapes[2][panel] = normal[2];
    }
    let medium = Medium {
        density: request.medium.gas.density,
        sound_speed: request.medium.gas.sound_speed,
    };
    let solve_frequency = |frequency_hz: f64| {
        let omega = core::f64::consts::TAU * frequency_hz;
        let velocity: Vec<Vec<C64>> = shapes
            .iter()
            .map(|shape| {
                shape
                    .iter()
                    .map(|value| C64::new(0.0, value / omega))
                    .collect()
            })
            .collect();
        let fields: Vec<&[C64]> = velocity.iter().map(Vec::as_slice).collect();
        let k = omega / medium.sound_speed;
        let formulation = if k * radius_m < 0.5 {
            HelmholtzFormulation::PlainCbie
        } else {
            HelmholtzFormulation::BurtonMiller
        };
        solve_radiation_batch(&surface, k, medium, &fields, formulation)
            .map(|solutions| (omega, solutions))
    };
    let mut training = Vec::with_capacity(request.training_frequency_hz.len());
    let mut held_out = Vec::with_capacity(request.held_out_frequency_hz.len());
    for &frequency_hz in request.training_frequency_hz {
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let (omega_rad_s, solutions) = solve_frequency(frequency_hz)?;
        let tables = checked_rigid_directivity_tables(
            &surface,
            &solutions,
            medium,
            request.directivity.maximum_spherical_harmonic_degree,
            request.directivity.minimum_captured_fraction,
            "training",
            frequency_hz,
        )?;
        let captured_fraction = captured_fraction(&solutions, &tables);
        training.push(ComplexShTrainingSample {
            omega_rad_s,
            coefficients_by_input: tables.into_iter().map(|table| table.coefficients).collect(),
            diagnostics: radiation_sample_diagnostics(&solutions, captured_fraction),
        });
    }
    for &frequency_hz in request.held_out_frequency_hz {
        cx.checkpoint()
            .map_err(|_| StructuralModalBasisError::Cancelled)?;
        let (omega_rad_s, solutions) = solve_frequency(frequency_hz)?;
        let tables = checked_rigid_directivity_tables(
            &surface,
            &solutions,
            medium,
            request.directivity.maximum_spherical_harmonic_degree,
            request.directivity.minimum_captured_fraction,
            "held-out",
            frequency_hz,
        )?;
        let captured_fraction = captured_fraction(&solutions, &tables);
        held_out.push(DirectFarFieldHeldOutSample {
            omega_rad_s,
            directions: request.held_out_directions_body.to_vec(),
            far_field_by_input: solutions
                .iter()
                .map(|solution| {
                    far_field(&surface, solution, medium, request.held_out_directions_body)
                })
                .collect(),
            diagnostics: radiation_sample_diagnostics(&solutions, captured_fraction),
        });
    }
    let source_identity = rigid_disc_broadband_sample_identity(request);
    let samples = SampledRadiationData {
        source_id: format!("euler-rigid-disc-acceleration-pa-m-v1:{source_identity}"),
        harmonic_time_convention: HarmonicTimeConvention::ExpNegativeIOmegaT,
        l_max: request.directivity.maximum_spherical_harmonic_degree,
        input_ids: vec![
            "translation-x[m/s2]".to_owned(),
            "translation-y[m/s2]".to_owned(),
            "translation-z[m/s2]".to_owned(),
        ],
        training,
        held_out,
    };
    cx.checkpoint()
        .map_err(|_| StructuralModalBasisError::Cancelled)?;
    let radiation = build_broadband_radiation_artifact(&samples, request.fit)?;
    let source_diameter_m = 2.0 * radius_m;
    let identity = rigid_disc_broadband_artifact_identity(
        source_identity,
        basis.identity,
        request.medium.gas_model_identity,
        request.medium.gas.sound_speed,
        source_diameter_m,
        &radiation,
    );
    Ok(RigidDiscBroadbandRadiationArtifact {
        identity,
        authority: BroadbandRadiationAuthority::EstimateOnly,
        structural_basis_identity: basis.identity,
        gas_model_identity: request.medium.gas_model_identity,
        sound_speed_m_s: request.medium.gas.sound_speed,
        source_diameter_m,
        sample_rate_hz,
        coordinates: RIGID_DISC_ACOUSTIC_COORDINATES.to_vec(),
        radiation,
        no_claims: RIGID_DISC_BROADBAND_SOURCE_NO_CLAIM,
    })
}

fn structural_broadband_artifact_identity(
    source_identity: ContentHash,
    structural_basis_identity: ContentHash,
    gas_model_identity: ContentHash,
    sound_speed_m_s: f64,
    source_diameter_m: f64,
    damping_model_identity: ContentHash,
    modal_loss_factors: &[f64],
    radiation: &BroadbandRadiationArtifact,
) -> ContentHash {
    let mut h = DomainHasher::new(STRUCTURAL_BROADBAND_ARTIFACT_IDENTITY_DOMAIN);
    for identity in [
        source_identity,
        structural_basis_identity,
        gas_model_identity,
        damping_model_identity,
    ] {
        h.update(identity.as_bytes());
    }
    hash_f64s(
        &mut h,
        [
            sound_speed_m_s,
            source_diameter_m,
            radiation.sample_interval_s,
        ],
    );
    hash_usizes(
        &mut h,
        [
            radiation.l_max,
            radiation.inputs.len(),
            radiation.channels.len(),
        ],
    );
    hash_f64_slice(&mut h, modal_loss_factors);
    h.update(radiation.report.source_id.as_bytes());
    for channel in &radiation.channels {
        hash_usizes(&mut h, [channel.l]);
        h.update(&channel.signed_m.to_le_bytes());
    }
    for input in &radiation.inputs {
        h.update(input.id.as_bytes());
        for filter in &input.filters {
            hash_usizes(&mut h, [filter.n]);
            hash_f64_slice(&mut h, &filter.a);
            hash_f64_slice(&mut h, &filter.b);
            hash_f64_slice(&mut h, &filter.c);
            hash_f64s(&mut h, [filter.d, filter.e_leftover, filter.t_s]);
        }
    }
    h.finalize()
}

fn rigid_disc_broadband_sample_identity(
    request: &RigidDiscBroadbandRadiationRequest<'_>,
) -> ContentHash {
    let mut h = DomainHasher::new(RIGID_DISC_BROADBAND_SOURCE_IDENTITY_DOMAIN);
    for identity in [
        request.basis.identity,
        request.basis.operator_identity,
        request.medium.gas_model_identity,
    ] {
        h.update(identity.as_bytes());
    }
    hash_usizes(
        &mut h,
        [
            RIGID_DISC_ACOUSTIC_COORDINATES.len(),
            request.directivity.maximum_spherical_harmonic_degree,
            request.fit.fit_order,
            request.fit.fit_iterations,
        ],
    );
    h.update(request.fit.fit_weights.label().as_bytes());
    h.update(&[u8::from(request.fit.fit_d)]);
    hash_f64s(
        &mut h,
        [
            request.fit.sample_rate_hz,
            request.fit.minimum_captured_fraction,
            request.fit.far_field_signal_floor,
            request.fit.maximum_normalized_error,
            request.fit.rms_normalized_error,
            request.directivity.minimum_captured_fraction,
            request.medium.gas.temperature,
            request.medium.gas.pressure,
            request.medium.gas.density,
            request.medium.gas.sound_speed,
        ],
    );
    hash_f64_slice(&mut h, request.training_frequency_hz);
    hash_f64_slice(&mut h, request.held_out_frequency_hz);
    hash_vec3_slice(&mut h, request.held_out_directions_body);
    h.finalize()
}

fn rigid_disc_broadband_artifact_identity(
    source_identity: ContentHash,
    structural_basis_identity: ContentHash,
    gas_model_identity: ContentHash,
    sound_speed_m_s: f64,
    source_diameter_m: f64,
    radiation: &BroadbandRadiationArtifact,
) -> ContentHash {
    let mut h = DomainHasher::new(RIGID_DISC_BROADBAND_ARTIFACT_IDENTITY_DOMAIN);
    for identity in [
        source_identity,
        structural_basis_identity,
        gas_model_identity,
    ] {
        h.update(identity.as_bytes());
    }
    hash_f64s(
        &mut h,
        [
            sound_speed_m_s,
            source_diameter_m,
            radiation.sample_interval_s,
        ],
    );
    hash_usizes(
        &mut h,
        [
            radiation.l_max,
            radiation.inputs.len(),
            radiation.channels.len(),
        ],
    );
    h.update(radiation.report.source_id.as_bytes());
    for channel in &radiation.channels {
        hash_usizes(&mut h, [channel.l]);
        h.update(&channel.signed_m.to_le_bytes());
    }
    for input in &radiation.inputs {
        h.update(input.id.as_bytes());
        for filter in &input.filters {
            hash_usizes(&mut h, [filter.n]);
            hash_f64_slice(&mut h, &filter.a);
            hash_f64_slice(&mut h, &filter.b);
            hash_f64_slice(&mut h, &filter.c);
            hash_f64s(&mut h, [filter.d, filter.e_leftover, filter.t_s]);
        }
    }
    h.finalize()
}

fn checked_rigid_directivity_tables(
    surface: &SpherePanels,
    solutions: &[RadiationSolution],
    medium: Medium,
    l_max: usize,
    minimum_captured_fraction: f64,
    grid: &'static str,
    frequency_hz: f64,
) -> Result<Vec<DirectivityTable>, StructuralModalBasisError> {
    solutions
        .iter()
        .enumerate()
        .map(|(coordinate, solution)| {
            let table = directivity_sh_table(surface, solution, medium, l_max)?;
            if solution.radiated_power < 0.0 {
                let d = solution.power_diagnostics;
                let sh_norm: f64 = table.coefficients.iter().map(|value| value.norm_sq()).sum();
                let sh_power = sh_norm / (2.0 * medium.density * medium.sound_speed);
                let direct_power = if table.captured_fraction > 0.0 {
                    sh_power / table.captured_fraction
                } else {
                    0.0
                };
                return Err(StructuralModalBasisError::BroadbandNegativeRadiatedPower(format!(
                    "grid={grid} frequency_hz={frequency_hz:.17e} rigid_coordinate={:?} Ssurf=({:.17e},{:.17e})W interval=[{:.17e},{:.17e}]W positive={:.17e}W negative={:.17e}W apparent={:.17e}W Pref={:.17e}W signed_efficiency={:.17e} Pinf_sh={sh_power:.17e}W Pinf_direct={direct_power:.17e}W captured_fraction={:.17e} ppw={:.17e} condition_lower_bound={:.17e}",
                    RIGID_DISC_ACOUSTIC_COORDINATES[coordinate],
                    d.surface_power.re,
                    d.surface_power.im,
                    solution.radiated_power_roundoff_interval.0,
                    solution.radiated_power_roundoff_interval.1,
                    d.positive_real_power,
                    d.negative_real_power,
                    d.apparent_power,
                    d.plane_wave_reference_power,
                    d.surface_power.re / d.plane_wave_reference_power,
                    table.captured_fraction,
                    solution.panels_per_wavelength,
                    solution.condition_lower_bound,
                )));
            }
            if solution.radiated_power > 0.0
                && table.captured_fraction < minimum_captured_fraction
            {
                return Err(StructuralModalBasisError::DirectivityTruncation {
                    mode: coordinate,
                    captured_fraction: table.captured_fraction,
                    minimum_fraction: minimum_captured_fraction,
                });
            }
            Ok(table)
        })
        .collect()
}

fn checked_directivity_tables(
    surface: &SpherePanels,
    solutions: &[RadiationSolution],
    modes: &[StructuralMode],
    medium: Medium,
    l_max: usize,
    minimum_captured_fraction: f64,
    grid: &'static str,
    frequency_hz: f64,
) -> Result<Vec<DirectivityTable>, StructuralModalBasisError> {
    solutions
        .iter()
        .enumerate()
        .map(|(mode, solution)| {
            let table = directivity_sh_table(surface, solution, medium, l_max)?;
            if solution.radiated_power < 0.0 {
                let d = solution.power_diagnostics;
                let sh_norm: f64 = table.coefficients.iter().map(|value| value.norm_sq()).sum();
                let sh_power = sh_norm / (2.0 * medium.density * medium.sound_speed);
                let direct_power = if table.captured_fraction > 0.0 {
                    sh_power / table.captured_fraction
                } else { 0.0 };
                return Err(StructuralModalBasisError::BroadbandNegativeRadiatedPower(format!(
                    "grid={grid} frequency_hz={frequency_hz:.17e} shape={mode} shape_frequency_hz={:.17e} Ssurf=({:.17e},{:.17e})W interval=[{:.17e},{:.17e}]W positive={:.17e}W negative={:.17e}W apparent={:.17e}W Pref={:.17e}W signed_efficiency={:.17e} Pinf_sh={sh_power:.17e}W Pinf_direct={direct_power:.17e}W captured_fraction={:.17e} ppw={:.17e} condition_lower_bound={:.17e}",
                    modes[mode].frequency_hz, d.surface_power.re, d.surface_power.im,
                    solution.radiated_power_roundoff_interval.0, solution.radiated_power_roundoff_interval.1,
                    d.positive_real_power, d.negative_real_power, d.apparent_power,
                    d.plane_wave_reference_power, d.surface_power.re / d.plane_wave_reference_power,
                    table.captured_fraction, solution.panels_per_wavelength, solution.condition_lower_bound,
                )));
            }
            if solution.radiated_power > 0.0 && table.captured_fraction < minimum_captured_fraction
            {
                return Err(StructuralModalBasisError::DirectivityTruncation {
                    mode,
                    captured_fraction: table.captured_fraction,
                    minimum_fraction: minimum_captured_fraction,
                });
            }
            Ok(table)
        })
        .collect()
}

fn captured_fraction(solutions: &[RadiationSolution], tables: &[DirectivityTable]) -> f64 {
    solutions
        .iter()
        .zip(tables)
        .filter(|(solution, _)| solution.radiated_power > 0.0)
        .map(|(_, table)| table.captured_fraction)
        .fold(1.0_f64, f64::min)
}

fn solve_modal_acceleration_radiation_batch(
    surface: &SpherePanels,
    modes: &[StructuralMode],
    omega_rad_s: f64,
    medium: Medium,
    radius_m: f64,
) -> Result<(HelmholtzFormulation, Vec<RadiationSolution>), HelmholtzError> {
    let velocity: Vec<Vec<C64>> = modes
        .iter()
        .map(|mode| {
            mode.panel_normal_shape_per_sqrt_kg
                .iter()
                .map(|shape| C64::new(0.0, shape / omega_rad_s))
                .collect()
        })
        .collect();
    let fields: Vec<&[C64]> = velocity.iter().map(Vec::as_slice).collect();
    let k = omega_rad_s / medium.sound_speed;
    let formulation = if k * radius_m < 0.5 {
        HelmholtzFormulation::PlainCbie
    } else {
        HelmholtzFormulation::BurtonMiller
    };
    solve_radiation_batch(surface, k, medium, &fields, formulation)
        .map(|solutions| (formulation, solutions))
}

fn radiation_sample_diagnostics(
    solutions: &[RadiationSolution],
    captured_fraction: f64,
) -> RadiationSampleDiagnostics {
    RadiationSampleDiagnostics {
        captured_fraction,
        panels_per_wavelength: solutions
            .iter()
            .map(|solution| solution.panels_per_wavelength)
            .fold(f64::INFINITY, f64::min),
        condition_lower_bound: solutions
            .iter()
            .map(|solution| solution.condition_lower_bound)
            .fold(1.0_f64, f64::max),
    }
}

fn broadband_sample_identity(request: &StructuralBroadbandRadiationRequest<'_>) -> ContentHash {
    let mut h = DomainHasher::new(STRUCTURAL_BROADBAND_SOURCE_IDENTITY_DOMAIN);
    for identity in [
        request.basis.identity,
        request.basis.operator_identity,
        request.basis.material_state_identity,
        request.loss.damping_model_identity,
        request.medium.gas_model_identity,
    ] {
        h.update(identity.as_bytes());
    }
    hash_usizes(
        &mut h,
        [
            request.basis.enrichment_modes.len(),
            request.directivity.maximum_spherical_harmonic_degree,
            request.fit.fit_order,
            request.fit.fit_iterations,
        ],
    );
    h.update(request.fit.fit_weights.label().as_bytes());
    h.update(&[u8::from(request.fit.fit_d)]);
    hash_f64s(
        &mut h,
        [
            request.fit.sample_rate_hz,
            request.fit.minimum_captured_fraction,
            request.fit.far_field_signal_floor,
            request.fit.maximum_normalized_error,
            request.fit.rms_normalized_error,
            request.directivity.minimum_captured_fraction,
            request.medium.gas.temperature,
            request.medium.gas.pressure,
            request.medium.gas.density,
            request.medium.gas.sound_speed,
        ],
    );
    hash_f64_slice(&mut h, &request.loss.loss_factors);
    hash_f64_slice(&mut h, request.training_frequency_hz);
    hash_f64_slice(&mut h, request.held_out_frequency_hz);
    hash_vec3_slice(&mut h, request.held_out_directions_body);
    h.finalize()
}

impl StructuralBroadbandRadiationArtifact {
    fn recomputed_identity(&self) -> Result<ContentHash, StructuralModalBasisError> {
        if self.no_claims == RIGID_DISC_BROADBAND_SOURCE_NO_CLAIM {
            return Ok(rigid_disc_broadband_artifact_identity(
                broadband_source_id_with_prefix(
                    &self.radiation,
                    "euler-rigid-disc-acceleration-pa-m-v1:",
                )?,
                self.structural_basis_identity,
                self.gas_model_identity,
                self.sound_speed_m_s,
                self.source_diameter_m,
                &self.radiation,
            ));
        }
        Ok(structural_broadband_artifact_identity(
            broadband_source_id_from_report(&self.radiation)?,
            self.structural_basis_identity,
            self.gas_model_identity,
            self.sound_speed_m_s,
            self.source_diameter_m,
            self.damping_model_identity,
            &self.modal_loss_factors,
            &self.radiation,
        ))
    }

    /// Bind this exact source bank to its residual basis and a modal budget.
    pub fn try_runtime<'a>(
        &'a self,
        basis: &'a StructuralResidualFlexibilityEstimateBasis,
        budget: ModalAcousticTimeBudget,
    ) -> Result<StructuralBroadbandSourceRuntime<'a>, StructuralModalBasisError> {
        if self.authority != BroadbandRadiationAuthority::EstimateOnly
            || self.no_claims != STRUCTURAL_BROADBAND_SOURCE_NO_CLAIM
            || self.recomputed_identity()? != self.identity
            || basis.recomputed_identity() != basis.identity
            || self.structural_basis_identity != basis.identity
            || self.gas_model_identity == ContentHash([0; 32])
            || !(self.sound_speed_m_s > 0.0 && self.sound_speed_m_s.is_finite())
            || !(self.source_diameter_m > 0.0 && self.source_diameter_m.is_finite())
            || self.damping_model_identity == ContentHash([0; 32])
            || self.modal_loss_factors.len() != basis.enrichment_modes.len()
            || self.radiation.inputs.len() != basis.enrichment_modes.len()
            || !self
                .radiation
                .report
                .source_id
                .starts_with("euler-structural-modal-acceleration-pa-m-v1:")
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "broadband source artifact does not match the exact residual basis",
            });
        }
        let modes = basis
            .enrichment_modes
            .iter()
            .zip(&self.modal_loss_factors)
            .map(|(mode, loss)| ModalAcousticMode {
                angular_frequency_rad_s: mode.angular_frequency_rad_s,
                damping_ratio: loss_factor_to_zeta(*loss),
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect();
        Ok(StructuralBroadbandSourceRuntime {
            basis,
            source: self,
            modal_runtime: ModalAcousticTimeModel::try_new(self.sample_rate_hz, modes, budget)?,
            modal_damping_ratios: self
                .modal_loss_factors
                .iter()
                .map(|loss| loss_factor_to_zeta(*loss))
                .collect(),
            radiation_runtime: self.radiation.try_runtime()?,
            modal_acceleration: vec![0.0; basis.enrichment_modes.len()],
        })
    }
}

impl RigidDiscBroadbandRadiationArtifact {
    fn recomputed_identity(&self) -> Result<ContentHash, StructuralModalBasisError> {
        let source_identity = broadband_source_id_with_prefix(
            &self.radiation,
            "euler-rigid-disc-acceleration-pa-m-v1:",
        )?;
        Ok(rigid_disc_broadband_artifact_identity(
            source_identity,
            self.structural_basis_identity,
            self.gas_model_identity,
            self.sound_speed_m_s,
            self.source_diameter_m,
            &self.radiation,
        ))
    }

    /// Apply the fitted rigid-disc radiation bank directly to acceleration
    /// sampled at accepted mechanics boundaries and anti-aliased before each
    /// factor-two decimation.
    pub fn synthesize_decimated_acceleration(
        &self,
        acceleration: &DecimatedModalAcceleration,
        basis: &StructuralResidualFlexibilityEstimateBasis,
        cx: &Cx<'_>,
    ) -> Result<StructuralBroadbandSourceStem, StructuralModalBasisError> {
        if self.authority != BroadbandRadiationAuthority::EstimateOnly
            || self.no_claims != RIGID_DISC_BROADBAND_SOURCE_NO_CLAIM
            || self.recomputed_identity()? != self.identity
            || basis.recomputed_identity() != basis.identity
            || self.structural_basis_identity != basis.identity
            || self.gas_model_identity == ContentHash([0; 32])
            || !(self.sound_speed_m_s > 0.0 && self.sound_speed_m_s.is_finite())
            || !(self.source_diameter_m > 0.0 && self.source_diameter_m.is_finite())
            || self.coordinates != RIGID_DISC_ACOUSTIC_COORDINATES
            || self.radiation.inputs.len() != self.coordinates.len()
            || acceleration.plate_model_identity != [basis.identity, self.identity]
            || acceleration.sample_rate_hz != self.sample_rate_hz
            || acceleration.coordinate_count() != self.coordinates.len()
            || acceleration.frame_count() == 0
        {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "decimated rigid acceleration does not match the disc radiation artifact",
            });
        }
        let channel_count = self.radiation.channels.len();
        let coefficient_count = acceleration
            .frame_count()
            .checked_mul(channel_count)
            .ok_or(StructuralModalBasisError::PressureCapacity {
                requested: usize::MAX,
            })?;
        let mut coefficients = Vec::new();
        coefficients
            .try_reserve_exact(coefficient_count)
            .map_err(|_| StructuralModalBasisError::PressureCapacity {
                requested: coefficient_count,
            })?;
        let mut runtime = self.radiation.try_runtime()?;
        for frame in 0..acceleration.frame_count() {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let outputs = runtime.step(
                acceleration
                    .frame(frame)
                    .expect("admitted rigid acceleration frame"),
            )?;
            coefficients.extend(outputs.iter().copied().map(FarFieldSourceCoefficientPaM));
        }
        Ok(StructuralBroadbandSourceStem {
            start_time_s: acceleration.start_time_s,
            sample_rate_hz: acceleration.sample_rate_hz,
            channels: self.radiation.channels.clone(),
            coefficients,
            authority: self.authority,
            source_identity: self.identity,
            structural_basis_identity: basis.identity,
        })
    }
}

impl StructuralBroadbandSourceRuntime<'_> {
    /// Project the moving contact, conservatively reconstruct its modal force
    /// measures, advance the oscillators, and emit closing-boundary Pa-m
    /// source coefficients. A listener/pressure/WAV stage is intentionally absent.
    pub fn synthesize_control_stream(
        &mut self,
        controls: &EulerControlStream<'_>,
        initial_state: PhysicalModalInitialState,
        force_reconstruction: GeneralizedForceReconstructionInput,
        maximum_contact_distance_m: f64,
        cx: &Cx<'_>,
    ) -> Result<StructuralBroadbandSourceStem, StructuralModalBasisError> {
        if controls.source().metadata().specimen_profile_identity != self.basis.profile_identity {
            return Err(StructuralModalBasisError::IdentityMismatch {
                what: "control trajectory profile does not match broadband structural source",
            });
        }
        let intervals = controls.audio();
        let first = intervals
            .first()
            .ok_or(StructuralModalBasisError::ControlTimeline {
                what: "control stream has no positive-duration audio intervals",
            })?;
        let last = intervals
            .last()
            .expect("nonempty interval slice has a last item");
        if intervals
            .windows(2)
            .any(|pair| pair[0].end_time_s.to_bits() != pair[1].start_time_s.to_bits())
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "mechanics audio intervals are not exactly contiguous",
            });
        }
        let frame_count = fixed_rate_frame_count_with_roundoff_bound(
            first.start_time_s,
            last.end_time_s,
            self.source.sample_rate_hz,
            force_reconstruction.clock_roundoff_operation_count,
        )
        .ok_or(StructuralModalBasisError::ControlTimeline {
            what: "control horizon is not an integral number of broadband audio samples",
        })?;
        let mut modal_forces = Vec::with_capacity(intervals.len());
        for interval in intervals {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
            modal_forces.push(modal_force_for_control_interval(
                &self.basis.mesh,
                &self.basis.enrichment_modes,
                controls,
                interval,
                maximum_contact_distance_m,
            )?);
        }
        let reconstructed =
            reconstruct_projected_modal_forces(intervals, &modal_forces, force_reconstruction, cx)?;
        if reconstructed.sample_rate_hz != self.source.sample_rate_hz
            || reconstructed.frame_count() != frame_count
        {
            return Err(StructuralModalBasisError::ControlTimeline {
                what: "reconstructed force clock does not match broadband source clock",
            });
        }
        if initial_state == PhysicalModalInitialState::StaticEquilibriumAtFirstHeldForce {
            self.modal_runtime.initialize_static_equilibrium(
                reconstructed
                    .frame(0)
                    .expect("positive frame count has a force row"),
            )?;
        }
        let channel_count = self.source.radiation.channels.len();
        let capacity = frame_count.checked_mul(channel_count).ok_or(
            StructuralModalBasisError::PressureCapacity {
                requested: usize::MAX,
            },
        )?;
        let mut coefficients = Vec::new();
        coefficients.try_reserve_exact(capacity).map_err(|_| {
            StructuralModalBasisError::PressureCapacity {
                requested: capacity,
            }
        })?;
        for frame in 0..frame_count {
            if frame % 64 == 0 {
                cx.checkpoint()
                    .map_err(|_| StructuralModalBasisError::Cancelled)?;
            }
            let force = reconstructed
                .frame(frame)
                .expect("validated force frame count");
            self.modal_runtime
                .step_duration(force, self.modal_runtime.sample_period_s())?;
            write_closing_modal_acceleration(
                &self.basis.enrichment_modes,
                &self.modal_damping_ratios,
                self.modal_runtime.states(),
                force,
                &mut self.modal_acceleration,
            )?;
            coefficients.extend(
                self.radiation_runtime
                    .step(&self.modal_acceleration)?
                    .iter()
                    .copied()
                    .map(FarFieldSourceCoefficientPaM),
            );
        }
        Ok(StructuralBroadbandSourceStem {
            start_time_s: first.start_time_s,
            sample_rate_hz: self.source.sample_rate_hz,
            channels: self.source.radiation.channels.clone(),
            coefficients,
            authority: BroadbandRadiationAuthority::EstimateOnly,
            source_identity: self.source.identity,
            structural_basis_identity: self.basis.identity,
        })
    }
}

fn write_closing_modal_acceleration(
    modes: &[StructuralMode],
    damping: &[f64],
    states: &[fs_couple::modal_acoustic_time::ModalAcousticState],
    force: &[f64],
    out: &mut [f64],
) -> Result<(), StructuralModalBasisError> {
    if modes.len() != damping.len()
        || modes.len() != states.len()
        || modes.len() != force.len()
        || modes.len() != out.len()
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "closing modal-acceleration rows must have identical cardinality",
        });
    }
    for ((((mode, damping_ratio), state), applied), acceleration) in
        modes.iter().zip(damping).zip(states).zip(force).zip(out)
    {
        let omega = mode.angular_frequency_rad_s;
        let equilibrium = applied / (omega * omega);
        // Algebraically Q - 2*zeta*omega*qdot - omega^2*q. Centering the
        // elastic term preserves an exactly held static equilibrium.
        *acceleration = (-2.0 * damping_ratio * omega).mul_add(
            state.velocity_m_sqrt_kg_per_s,
            -omega * omega * (state.displacement_m_sqrt_kg - equilibrium),
        );
        if *acceleration == 0.0 {
            *acceleration = 0.0;
        }
        if !acceleration.is_finite() {
            return Err(StructuralModalBasisError::InvalidRequest {
                what: "closing modal acceleration became non-finite",
            });
        }
    }
    Ok(())
}

/// Synthesize simultaneous world-fixed pressure from the rigid-disc
/// acceleration bank. The internal carrier reuses the common retarded source
/// evaluator but is never exposed as an elastic modal artifact.
pub fn synthesize_rigid_disc_retarded_far_field_world_observers(
    source: &RigidDiscBroadbandRadiationArtifact,
    stem: &StructuralBroadbandSourceStem,
    basis: &StructuralResidualFlexibilityEstimateBasis,
    trajectory: &RenderTrajectory,
    observers: &[AcousticWorldObserver],
    controls: RetardedFarFieldObserverControls,
    cx: &Cx<'_>,
) -> Result<Vec<PhysicalPressureSignal>, StructuralModalBasisError> {
    if source.recomputed_identity()? != source.identity
        || source.structural_basis_identity != basis.identity
        || source.coordinates != RIGID_DISC_ACOUSTIC_COORDINATES
    {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "rigid-disc source does not match the exact structural mesh",
        });
    }
    let carrier = StructuralBroadbandRadiationArtifact {
        identity: source.identity,
        authority: source.authority,
        structural_basis_identity: source.structural_basis_identity,
        gas_model_identity: source.gas_model_identity,
        sound_speed_m_s: source.sound_speed_m_s,
        source_diameter_m: source.source_diameter_m,
        damping_model_identity: source.identity,
        sample_rate_hz: source.sample_rate_hz,
        modal_loss_factors: Vec::new(),
        radiation: source.radiation.clone(),
        no_claims: source.no_claims,
    };
    synthesize_retarded_far_field_world_observers(
        &carrier, stem, basis, trajectory, observers, controls, cx,
    )
}

/// Synthesize simultaneous world-fixed pressure signals from body-frame
/// broadband far-field coefficients using physical retarded time.
///
/// Version one is interior-only: every accepted emission query has a complete
/// 16-tap Lanczos-8 stencil. All observers are validated and synthesized into
/// private candidates before any signal is returned.
pub fn synthesize_retarded_far_field_world_observers(
    source: &StructuralBroadbandRadiationArtifact,
    stem: &StructuralBroadbandSourceStem,
    basis: &StructuralResidualFlexibilityEstimateBasis,
    trajectory: &RenderTrajectory,
    observers: &[AcousticWorldObserver],
    controls: RetardedFarFieldObserverControls,
    cx: &Cx<'_>,
) -> Result<Vec<PhysicalPressureSignal>, StructuralModalBasisError> {
    validate_retarded_observer_inputs(source, stem, basis, trajectory, observers, controls, cx)?;
    let fs = f64::from(stem.sample_rate_hz);
    let period = fs.recip();
    let first_emission_s = stem.start_time_s + 8.0 * period;
    let last_emission_s = stem.start_time_s + (stem.frame_count() as f64 - 8.0) * period;
    let timeline = TimelineResampler::new(trajectory);
    let mut arrival_bounds = Vec::with_capacity(observers.len());
    for observer in observers {
        let first = arrival_time(
            &timeline,
            *observer,
            first_emission_s,
            source.sound_speed_m_s,
        )?;
        let last = arrival_time(
            &timeline,
            *observer,
            last_emission_s,
            source.sound_speed_m_s,
        )?;
        arrival_bounds.push((first, last));
    }
    let common_start = arrival_bounds
        .iter()
        .map(|bound| bound.0)
        .fold(f64::NEG_INFINITY, f64::max);
    let common_end = arrival_bounds
        .iter()
        .map(|bound| bound.1)
        .fold(f64::INFINITY, f64::min);
    let first_frame = ((common_start - stem.start_time_s) * fs).ceil().max(1.0) as usize;
    let end_frame = ((common_end - stem.start_time_s) * fs).floor() as usize;
    let frame_count = end_frame
        .checked_sub(first_frame)
        .and_then(|span| span.checked_add(1))
        .unwrap_or(0);
    if frame_count == 0 || frame_count > controls.maximum_output_frames {
        return Err(StructuralModalBasisError::ControlTimeline {
            what: "retarded observers have no bounded common interior arrival horizon",
        });
    }
    let output_start_s = stem.start_time_s + (first_frame - 1) as f64 * period;
    let capacity = frame_count.checked_mul(observers.len()).ok_or(
        StructuralModalBasisError::PressureCapacity {
            requested: usize::MAX,
        },
    )?;
    let mut pressure = Vec::new();
    pressure.try_reserve_exact(observers.len()).map_err(|_| {
        StructuralModalBasisError::PressureCapacity {
            requested: capacity,
        }
    })?;
    for _ in observers {
        let mut observer_pressure = Vec::new();
        observer_pressure
            .try_reserve_exact(frame_count)
            .map_err(|_| StructuralModalBasisError::PressureCapacity {
                requested: capacity,
            })?;
        pressure.push(observer_pressure);
    }
    let mut coefficient_scratch = vec![0.0; stem.channels.len()];
    let mut complex_scratch = vec![C64::ZERO; stem.channels.len()];
    let minimum_far_field_m = minimum_broadband_far_field_distance(source, fs);
    for frame in 0..frame_count {
        if frame % 64 == 0 {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
        }
        let arrival_s = stem.start_time_s + (first_frame + frame) as f64 * period;
        for (observer_index, observer) in observers.iter().enumerate() {
            let (emission_s, sample, direction_world, radius_m) = retarded_emission(
                &timeline,
                *observer,
                arrival_s,
                first_emission_s,
                last_emission_s,
                source.sound_speed_m_s,
                controls,
            )?;
            validate_retarded_state(
                source,
                trajectory,
                &sample,
                radius_m,
                minimum_far_field_m,
                controls.maximum_surface_mach,
            )?;
            interpolate_stem_lanczos8(stem, emission_s, &mut coefficient_scratch)?;
            let direction_body = sample
                .state
                .pose()
                .orientation()
                .rotate_world_to_body(direction_world);
            for (complex, real) in complex_scratch.iter_mut().zip(&coefficient_scratch) {
                *complex = C64::from_re(*real);
            }
            let far = evaluate_real_tesseral(
                source.radiation.l_max,
                &complex_scratch,
                [direction_body.x, direction_body.y, direction_body.z],
            )?;
            let value = far.re / radius_m;
            if !value.is_finite() {
                return Err(StructuralModalBasisError::InvalidRequest {
                    what: "retarded far-field pressure became non-finite",
                });
            }
            pressure[observer_index].push(value);
        }
    }
    Ok(observers
        .iter()
        .copied()
        .zip(pressure)
        .map(|(observer, pressure_pa)| {
            let peak_abs_pressure_pa = pressure_pa
                .iter()
                .fold(0.0_f64, |peak, value| peak.max(value.abs()));
            let observer = PhysicalPressureObserver::WorldFixed(observer);
            let identity = retarded_pressure_identity(
                source,
                stem,
                trajectory,
                observer,
                output_start_s,
                controls,
                &pressure_pa,
            );
            PhysicalPressureSignal {
                start_time_s: output_start_s,
                sample_rate_hz: stem.sample_rate_hz,
                pressure_pa,
                peak_abs_pressure_pa,
                contact_force_sampling:
                    PhysicalContactForceSampling::IntervalMeasureAtClosingElseOpeningEndpointBandLimitedV1,
                observer,
                structural_basis_identity: basis.identity,
                radiation_identity: source.identity,
                damping_model_identity: source.damping_model_identity,
                identity,
            }
        })
        .collect())
}

fn validate_retarded_observer_inputs(
    source: &StructuralBroadbandRadiationArtifact,
    stem: &StructuralBroadbandSourceStem,
    basis: &StructuralResidualFlexibilityEstimateBasis,
    trajectory: &RenderTrajectory,
    observers: &[AcousticWorldObserver],
    controls: RetardedFarFieldObserverControls,
    cx: &Cx<'_>,
) -> Result<(), StructuralModalBasisError> {
    drop(source.radiation.try_runtime()?);
    if source.identity == ContentHash([0; 32])
        || source.authority != BroadbandRadiationAuthority::EstimateOnly
        || !matches!(
            source.no_claims,
            STRUCTURAL_BROADBAND_SOURCE_NO_CLAIM | RIGID_DISC_BROADBAND_SOURCE_NO_CLAIM
        )
        || source.structural_basis_identity != basis.identity
        || basis.recomputed_identity() != basis.identity
        || trajectory.metadata().specimen_profile_identity != basis.profile_identity
        || source.gas_model_identity == ContentHash([0; 32])
        || !(source.sound_speed_m_s > 0.0 && source.sound_speed_m_s.is_finite())
        || !(source.source_diameter_m > 0.0 && source.source_diameter_m.is_finite())
        || stem.authority != BroadbandRadiationAuthority::EstimateOnly
        || stem.source_identity != source.identity
        || stem.structural_basis_identity != basis.identity
        || stem.sample_rate_hz != source.sample_rate_hz
        || source.radiation.sample_interval_s.to_bits()
            != f64::from(source.sample_rate_hz).recip().to_bits()
        || stem.channels != source.radiation.channels
        || stem.coefficients.len() % stem.channels.len() != 0
        || stem.frame_count() < 16
        || !stem.start_time_s.is_finite()
        || stem.coefficients.iter().any(|value| !value.0.is_finite())
    {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "retarded observer source, stem, trajectory, gas, or basis identity is inconsistent",
        });
    }
    if observers.is_empty()
        || observers.len() > MAX_PHYSICAL_PRESSURE_OBSERVERS
        || observers.iter().any(|observer| {
            observer
                .position_world_m
                .iter()
                .any(|value| !value.is_finite())
        })
        || !(controls.maximum_surface_mach > 0.0 && controls.maximum_surface_mach <= 0.1)
        || !(controls.root_time_tolerance_s > 0.0 && controls.root_time_tolerance_s.is_finite())
        || controls.maximum_root_iterations == 0
        || controls.maximum_root_iterations > 256
        || controls.interpolation_radius_frames != 8
        || controls.maximum_output_frames == 0
    {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "retarded observer set or controls are outside the version-one domain",
        });
    }
    if source.recomputed_identity()? != source.identity {
        return Err(StructuralModalBasisError::IdentityMismatch {
            what: "retarded observer broadband artifact identity does not recompute",
        });
    }
    let timeline = TimelineResampler::new(trajectory);
    let minimum = minimum_broadband_far_field_distance(source, f64::from(stem.sample_rate_hz));
    for (index, sample) in trajectory.samples().iter().enumerate() {
        if index % 64 == 0 {
            cx.checkpoint()
                .map_err(|_| StructuralModalBasisError::Cancelled)?;
        }
        let state = sample.state();
        let reconstructed =
            timeline.sample(sample.input().time_s, EventEvaluationSide::RightLimit)?;
        for observer in observers {
            let radius = observer_vector(*observer, state.pose().position_world()).1;
            validate_retarded_state(
                source,
                trajectory,
                &reconstructed,
                radius,
                minimum,
                controls.maximum_surface_mach,
            )?;
        }
    }
    Ok(())
}

fn broadband_source_id_from_report(
    radiation: &BroadbandRadiationArtifact,
) -> Result<ContentHash, StructuralModalBasisError> {
    broadband_source_id_with_prefix(radiation, "euler-structural-modal-acceleration-pa-m-v1:")
}

fn broadband_source_id_with_prefix(
    radiation: &BroadbandRadiationArtifact,
    prefix: &str,
) -> Result<ContentHash, StructuralModalBasisError> {
    let hex = radiation.report.source_id.strip_prefix(prefix).ok_or(
        StructuralModalBasisError::IdentityMismatch {
            what: "broadband producer source id has an unsupported schema",
        },
    )?;
    ContentHash::from_hex(hex).ok_or(StructuralModalBasisError::IdentityMismatch {
        what: "broadband producer source id has an invalid digest",
    })
}

fn minimum_broadband_far_field_distance(
    source: &StructuralBroadbandRadiationArtifact,
    sample_rate_hz: f64,
) -> f64 {
    let wavelength = 2.0 * source.sound_speed_m_s / sample_rate_hz;
    (2.0 * source.source_diameter_m)
        .max(2.0 * source.source_diameter_m * source.source_diameter_m / wavelength)
}

fn validate_retarded_state(
    source: &StructuralBroadbandRadiationArtifact,
    trajectory: &RenderTrajectory,
    sample: &crate::timeline_resampling::ResampledTimelineSample,
    radius_m: f64,
    minimum_far_field_m: f64,
    maximum_surface_mach: f64,
) -> Result<(), StructuralModalBasisError> {
    if !(radius_m >= minimum_far_field_m && radius_m.is_finite()) {
        return Err(StructuralModalBasisError::ObserverOutsideFarField {
            distance_m: radius_m,
            minimum_m: minimum_far_field_m,
            mode: 0,
        });
    }
    let properties = trajectory.metadata().mass_properties.properties;
    let kinematics_error = |_| StructuralModalBasisError::InvalidRequest {
        what: "retarded trajectory kinematics are invalid",
    };
    let linear = sample
        .state
        .center_of_mass_velocity_world(properties)
        .map_err(kinematics_error)?;
    let angular = properties
        .angular_velocity_body_checked(sample.state.angular_momentum_body())
        .map_err(kinematics_error)?;
    let speed = det::sqrt(linear.norm_squared())
        + 0.5 * source.source_diameter_m * det::sqrt(angular.norm_squared());
    if !speed.is_finite() || speed / source.sound_speed_m_s > maximum_surface_mach {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "retarded far-field source exceeds the admitted surface Mach number",
        });
    }
    Ok(())
}

fn observer_vector(observer: AcousticWorldObserver, source_world: Vec3) -> (Vec3, f64) {
    let offset = Vec3::new(
        observer.position_world_m[0] - source_world.x,
        observer.position_world_m[1] - source_world.y,
        observer.position_world_m[2] - source_world.z,
    );
    let radius = det::sqrt(offset.norm_squared());
    (offset, radius)
}

fn arrival_time(
    timeline: &TimelineResampler<'_>,
    observer: AcousticWorldObserver,
    emission_s: f64,
    sound_speed_m_s: f64,
) -> Result<f64, StructuralModalBasisError> {
    let sample = timeline.sample(emission_s, EventEvaluationSide::RightLimit)?;
    let (_, radius) = observer_vector(observer, sample.state.pose().position_world());
    let arrival = emission_s + radius / sound_speed_m_s;
    if !(radius > 0.0 && radius.is_finite() && arrival.is_finite()) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "retarded observer radius or arrival time is invalid",
        });
    }
    Ok(arrival)
}

fn retarded_emission(
    timeline: &TimelineResampler<'_>,
    observer: AcousticWorldObserver,
    arrival_s: f64,
    low: f64,
    high: f64,
    sound_speed_m_s: f64,
    controls: RetardedFarFieldObserverControls,
) -> Result<
    (
        f64,
        crate::timeline_resampling::ResampledTimelineSample,
        Vec3,
        f64,
    ),
    StructuralModalBasisError,
> {
    let emission_s = bisect_retarded_emission_time(
        arrival_s,
        low,
        high,
        controls.root_time_tolerance_s,
        controls.maximum_root_iterations,
        |time_s| arrival_time(timeline, observer, time_s, sound_speed_m_s),
    )?;
    let sample = timeline.sample(emission_s, EventEvaluationSide::RightLimit)?;
    let (offset, radius_m) = observer_vector(observer, sample.state.pose().position_world());
    if !(radius_m > 0.0 && radius_m.is_finite()) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "retarded observer coincides with the source",
        });
    }
    Ok((emission_s, sample, offset.scale(radius_m.recip()), radius_m))
}

fn bisect_retarded_emission_time(
    arrival_s: f64,
    mut low: f64,
    mut high: f64,
    tolerance_s: f64,
    maximum_iterations: u32,
    mut arrival_at: impl FnMut(f64) -> Result<f64, StructuralModalBasisError>,
) -> Result<f64, StructuralModalBasisError> {
    let mut low_residual = arrival_at(low)? - arrival_s;
    let high_residual = arrival_at(high)? - arrival_s;
    if low_residual > 0.0 || high_residual < 0.0 {
        return Err(StructuralModalBasisError::ControlTimeline {
            what: "arrival time has no emission root inside the complete source stencil",
        });
    }
    for _ in 0..maximum_iterations {
        let mid = 0.5 * (low + high);
        let residual = arrival_at(mid)? - arrival_s;
        if residual <= 0.0 {
            low = mid;
            low_residual = residual;
        } else {
            high = mid;
        }
        if high - low <= tolerance_s || residual == 0.0 {
            break;
        }
    }
    if high - low > tolerance_s && low_residual != 0.0 {
        return Err(StructuralModalBasisError::ControlTimeline {
            what: "retarded emission root exceeded the deterministic bisection budget",
        });
    }
    Ok(if low_residual == 0.0 {
        low
    } else {
        0.5 * (low + high)
    })
}

fn interpolate_stem_lanczos8(
    stem: &StructuralBroadbandSourceStem,
    emission_s: f64,
    out: &mut [f64],
) -> Result<(), StructuralModalBasisError> {
    if out.len() != stem.channels.len() {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "retarded interpolation scratch has the wrong channel count",
        });
    }
    let u = (emission_s - stem.start_time_s) * f64::from(stem.sample_rate_hz) - 1.0;
    let base = u.floor() as isize;
    if base < 7 || base + 8 >= stem.frame_count() as isize {
        return Err(StructuralModalBasisError::ControlTimeline {
            what: "retarded coefficient query lacks a complete Lanczos-8 stencil",
        });
    }
    out.fill(0.0);
    for tap in -7_isize..=8 {
        let index = usize::try_from(base + tap).expect("prevalidated nonnegative stem index");
        let x = u - index as f64;
        let weight = lanczos8(x);
        for (sum, coefficient) in out
            .iter_mut()
            .zip(stem.frame(index).expect("complete stem stencil").iter())
        {
            *sum = weight.mul_add(coefficient.0, *sum);
        }
    }
    if out.iter().any(|value| !value.is_finite()) {
        return Err(StructuralModalBasisError::InvalidRequest {
            what: "retarded Lanczos interpolation became non-finite",
        });
    }
    Ok(())
}

fn lanczos8(x: f64) -> f64 {
    if x == 0.0 {
        1.0
    } else if x.abs() >= 8.0 {
        0.0
    } else {
        let pix = core::f64::consts::PI * x;
        (det::sin(pix) / pix) * (det::sin(pix / 8.0) / (pix / 8.0))
    }
}

fn retarded_pressure_identity(
    source: &StructuralBroadbandRadiationArtifact,
    stem: &StructuralBroadbandSourceStem,
    trajectory: &RenderTrajectory,
    observer: PhysicalPressureObserver,
    start_time_s: f64,
    controls: RetardedFarFieldObserverControls,
    pressure_pa: &[f64],
) -> ContentHash {
    let mut h = DomainHasher::new(RETARDED_FAR_FIELD_SIGNAL_IDENTITY_DOMAIN);
    h.update(source.identity.as_bytes());
    h.update(stem.source_identity.as_bytes());
    h.update(trajectory.metadata().configuration_identity.as_bytes());
    hash_pressure_observer(&mut h, observer);
    hash_f64s(
        &mut h,
        [
            start_time_s,
            controls.maximum_surface_mach,
            controls.root_time_tolerance_s,
        ],
    );
    hash_usizes(
        &mut h,
        [
            controls.maximum_root_iterations as usize,
            controls.interpolation_radius_frames as usize,
            pressure_pa.len(),
        ],
    );
    hash_f64_slice(&mut h, pressure_pa);
    h.finalize()
}

impl StructuralModalBasis {
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
        let surface = SpherePanels::from_triangles(
            self.mesh
                .boundary
                .triangles
                .iter()
                .map(|triangle| triangle.map(|node| self.mesh.nodes_m[node]))
                .collect(),
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
        let surface = SpherePanels::from_triangles(
            self.mesh
                .boundary
                .triangles
                .iter()
                .map(|triangle| triangle.map(|node| self.mesh.nodes_m[node]))
                .collect(),
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

fn residual_flexibility_basis_identity(
    basis: &StructuralResidualFlexibilityEstimateBasis,
) -> ContentHash {
    let mut hasher = DomainHasher::new(STRUCTURAL_RESIDUAL_FLEXIBILITY_IDENTITY_DOMAIN);
    hasher.update(&basis.schema_version.to_le_bytes());
    hasher.update(&[basis.authority.tag()]);
    hasher.update(basis.specimen_identity.as_bytes());
    hasher.update(basis.profile_identity.as_bytes());
    hasher.update(basis.material_state_identity.as_bytes());
    hasher.update(basis.operator_identity.as_bytes());
    hasher.update(structural_operator_identity(&basis.mesh, &basis.assembly).as_bytes());
    hash_f64s(
        &mut hasher,
        [
            basis.assembly_budget.minimum_scaled_jacobian,
            basis.maximum_rigid_stiffness_relative_residual,
            basis.maximum_mass_orthogonality_error,
            basis.forcing_frequency_band_hz.0,
            basis.forcing_frequency_band_hz.1,
            basis.enrichment_frequency_band_hz.0,
            basis.enrichment_frequency_band_hz.1,
            basis.slice_stats.shift,
        ],
    );
    hash_usizes(
        &mut hasher,
        [
            basis.assembly_budget.maximum_nodes,
            basis.assembly_budget.maximum_tetrahedra,
            basis.assembly_budget.maximum_free_dofs,
            basis.certified_enrichment_mode_count,
            basis.certified_in_band_mode_count,
            basis.certified_partition_mode_count,
            basis.slice_stats.factorizations,
            basis.slice_stats.lanczos_iters,
            basis.slice_stats.restarts,
            basis.slice_stats.factor_nnz_l,
            basis.slice_stats.factor_peak_bytes,
            basis.slice_stats.pivots_delayed,
            basis.rigid_modes_per_sqrt_kg.len(),
            basis.enrichment_modes.len(),
        ],
    );
    for mode in &basis.rigid_modes_per_sqrt_kg {
        hash_f64_slice(&mut hasher, mode);
    }
    for mode in &basis.enrichment_modes {
        hash_f64s(
            &mut hasher,
            [
                mode.eigenvalue_s2,
                mode.angular_frequency_rad_s,
                mode.frequency_hz,
                mode.eigenvalue_interval_s2.0,
                mode.eigenvalue_interval_s2.1,
                mode.eigenvalue_residual_s2,
            ],
        );
        hash_vec3_slice(&mut hasher, &mode.nodal_shape_per_sqrt_kg);
        hash_f64_slice(&mut hasher, &mode.panel_normal_shape_per_sqrt_kg);
    }
    hasher.finalize()
}

fn structural_operator_identity(
    mesh: &RoundedCylinderTetMesh,
    assembly: &TetElasticAssembly,
) -> ContentHash {
    let mut hasher = DomainHasher::new(
        "org.frankensim.fs-euler-disc-e2e.structural-operator-and-p1-panel-normal.v1",
    );
    hash_usizes(
        &mut hasher,
        [
            mesh.nodes_m.len(),
            mesh.tetrahedra.len(),
            mesh.boundary.triangles.len(),
            assembly.free_dofs.len(),
            assembly.element_volumes_m3.len(),
        ],
    );
    hash_f64s(&mut hasher, mesh.nodes_m.iter().flatten().copied());
    for tetrahedron in &mesh.tetrahedra {
        hash_usizes(&mut hasher, tetrahedron.iter().copied());
    }
    for index in 0..mesh.boundary.triangles.len() {
        hash_usizes(&mut hasher, mesh.boundary.triangles[index]);
        hash_f64s(
            &mut hasher,
            mesh.boundary.centroids_m[index]
                .into_iter()
                .chain(mesh.boundary.normals[index])
                .chain([mesh.boundary.areas_m2[index]]),
        );
    }
    for matrix in [&assembly.stiffness, &assembly.mass] {
        hash_usizes(&mut hasher, [matrix.nrows(), matrix.ncols()]);
        for row in 0..matrix.nrows() {
            let (columns, values) = matrix.row(row);
            hash_usizes(&mut hasher, [columns.len()]);
            hash_usizes(&mut hasher, columns.iter().copied());
            hash_f64s(&mut hasher, values.iter().copied());
        }
    }
    hash_usizes(&mut hasher, assembly.free_dofs.iter().copied());
    hash_f64s(
        &mut hasher,
        [
            mesh.maximum_meridian_chord_error_m,
            mesh.maximum_azimuthal_chord_error_m,
            assembly.total_mass_kg,
            assembly.minimum_scaled_jacobian,
        ]
        .into_iter()
        .chain(assembly.element_volumes_m3.iter().copied()),
    );
    hasher.finalize()
}

fn residual_flexibility_response_identity(
    response: &StructuralResidualFlexibilityEstimateResponse,
) -> ContentHash {
    let mut hasher = DomainHasher::new(STRUCTURAL_RESIDUAL_RESPONSE_IDENTITY_DOMAIN);
    hasher.update(response.basis_identity.as_bytes());
    hasher.update(&[response.authority.tag()]);
    hash_f64s(
        &mut hasher,
        response
            .requested_point_m
            .into_iter()
            .chain(response.applied_force_n)
            .chain(response.force_projection.closest_point_m)
            .chain(response.force_projection.barycentric)
            .chain([
                response.force_projection.distance_to_boundary_m,
                response.maximum_rigid_force_relative_residual,
                response.elastic_work_j,
                response.recoverable_strain_energy_j,
                response.energy_closure_residual_j,
            ]),
    );
    hash_usizes(&mut hasher, [response.force_projection.boundary_triangle]);
    hash_f64_slice(
        &mut hasher,
        &response.force_projection.modal_force_n_per_sqrt_kg,
    );
    hash_vec3_slice(&mut hasher, &response.inertia_relieved_nodal_force_n);
    hash_f64_slice(&mut hasher, &response.modal_displacement_m_sqrt_kg);
    hash_vec3_slice(&mut hasher, &response.nodal_displacement_m);
    hash_f64_slice(&mut hasher, &response.panel_normal_displacement_m);
    hasher.finalize()
}

fn hash_f64s(hasher: &mut DomainHasher, values: impl IntoIterator<Item = f64>) {
    for value in values {
        hasher.update(&value.to_bits().to_le_bytes());
    }
}

fn hash_f64_slice(hasher: &mut DomainHasher, values: &[f64]) {
    hash_usizes(hasher, [values.len()]);
    hash_f64s(hasher, values.iter().copied());
}

fn hash_vec3_slice(hasher: &mut DomainHasher, values: &[[f64; 3]]) {
    hash_usizes(hasher, [values.len()]);
    hash_f64s(hasher, values.iter().flatten().copied());
}

fn hash_usizes(hasher: &mut DomainHasher, values: impl IntoIterator<Item = usize>) {
    for value in values {
        hasher.update(&u64::try_from(value).unwrap_or(u64::MAX).to_le_bytes());
    }
}

fn pair_bits_equal(left: (f64, f64), right: (f64, f64)) -> bool {
    left.0.to_bits() == right.0.to_bits() && left.1.to_bits() == right.1.to_bits()
}

fn relative_l2_difference(coarse: &[f64], fine: &[f64]) -> f64 {
    let difference = coarse
        .iter()
        .zip(fine)
        .fold(0.0_f64, |norm, (left, right)| norm.hypot(left - right));
    let scale = fine.iter().fold(0.0_f64, |norm, value| norm.hypot(*value));
    difference / scale.max(f64::MIN_POSITIVE)
}

fn relative_vec3_l2_difference(coarse: &[[f64; 3]], fine: &[[f64; 3]]) -> f64 {
    let difference = coarse
        .iter()
        .flatten()
        .zip(fine.iter().flatten())
        .fold(0.0_f64, |norm, (left, right)| norm.hypot(left - right));
    let scale = fine
        .iter()
        .flatten()
        .fold(0.0_f64, |norm, value| norm.hypot(*value));
    difference / scale.max(f64::MIN_POSITIVE)
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

    #[test]
    fn g0_rectangular_hermite_projection_reproduces_nodal_slopes_without_edge_kinks() {
        let width_m = 2.0;
        let depth_m = 3.0;
        let nodes = [
            (0.0, 0.0),
            (width_m, 0.0),
            (0.0, depth_m),
            (width_m, depth_m),
        ];
        let value = |x: f64, y: f64| 1.0 + 2.0 * x + 3.0 * y + 4.0 * x * y;
        let slope_x = |_x: f64, y: f64| 2.0 + 4.0 * y;
        let slope_y = |x: f64, _y: f64| 3.0 + 4.0 * x;
        let mode = RectangularPlateMode {
            eigenvalue_s2: 1.0,
            angular_frequency_rad_s: 1.0,
            frequency_hz: 1.0 / core::f64::consts::TAU,
            eigenvalue_interval_s2: (1.0, 1.0),
            eigenvalue_residual_s2: 0.0,
            nodal_transverse_shape_per_sqrt_kg: nodes.iter().map(|&(x, y)| value(x, y)).collect(),
            nodal_slope_x_per_m_sqrt_kg: nodes.iter().map(|&(x, y)| slope_x(x, y)).collect(),
            nodal_slope_y_per_m_sqrt_kg: nodes.iter().map(|&(x, y)| slope_y(x, y)).collect(),
            nodal_mixed_slope_per_m2_sqrt_kg: vec![4.0; 4],
        };
        for (u, v) in [(0.0, 0.0), (1.0, 1.0), (0.37, 0.61), (1.0, 0.43)] {
            let observed = rectangular_mode_value(&mode, 0, 2, u, v, width_m, depth_m);
            let expected = value(u * width_m, v * depth_m);
            assert!((observed - expected).abs() <= 2.0e-14 * expected.abs().max(1.0));
        }
    }

    #[test]
    fn g0_compliant_disc_contact_maps_to_plate_reference_surface() {
        let penetration_m = -4.005_401e-7;
        assert_eq!(
            plate_reference_contact_point([0.012, -0.007, penetration_m], penetration_m)
                .expect("consistent disc-side contact and signed gap"),
            [0.012, -0.007, 0.0]
        );
        assert_eq!(
            plate_reference_contact_point([0.012, -0.007, 3.0e-6], penetration_m,)
                .expect("surface-resolved gap may differ from the smooth disc-point height"),
            [0.012, -0.007, 0.0]
        );
        assert!(plate_reference_contact_point([0.012, -0.007, f64::NAN], penetration_m).is_err());
        assert!(plate_reference_contact_point([0.012, -0.007, 0.0], f64::NAN).is_err());
    }

    fn pressure_signal(identity_byte: u8, pressure_pa: Vec<f64>) -> PhysicalPressureSignal {
        PhysicalPressureSignal {
            start_time_s: 0.0,
            sample_rate_hz: 48_000,
            peak_abs_pressure_pa: pressure_pa
                .iter()
                .fold(0.0_f64, |peak, value| peak.max(value.abs())),
            pressure_pa,
            contact_force_sampling:
                PhysicalContactForceSampling::IntervalMeanAtClosingElseOpeningEndpointZohV1,
            observer: PhysicalPressureObserver::WorldFixed(AcousticWorldObserver {
                position_world_m: [0.25, -0.1, 0.3],
            }),
            structural_basis_identity: ContentHash([identity_byte; 32]),
            radiation_identity: ContentHash([identity_byte.wrapping_add(1); 32]),
            damping_model_identity: ContentHash([identity_byte.wrapping_add(2); 32]),
            identity: ContentHash([identity_byte.wrapping_add(3); 32]),
        }
    }

    #[test]
    fn g0_physical_pressure_superposition_is_si_additive_and_order_invariant() {
        let disc = pressure_signal(0x20, vec![1.0, -2.0, 0.5]);
        let plate = pressure_signal(0x40, vec![0.25, 3.0, -0.5]);
        let (forward, reversed) = with_cx(|cx| {
            (
                superpose_pressure_signals(&[&disc, &plate], cx).unwrap(),
                superpose_pressure_signals(&[&plate, &disc], cx).unwrap(),
            )
        });
        assert_eq!(forward.pressure_pa, vec![1.25, 1.0, 0.0]);
        assert_eq!(forward.pressure_pa, reversed.pressure_pa);
        assert_eq!(forward.peak_abs_pressure_pa.to_bits(), 1.25_f64.to_bits());
        assert_eq!(forward.identity, reversed.identity);
        assert_eq!(
            forward.structural_basis_identity,
            reversed.structural_basis_identity
        );
        assert_ne!(
            forward.structural_basis_identity,
            disc.structural_basis_identity
        );
        assert_ne!(
            forward.structural_basis_identity,
            plate.structural_basis_identity
        );
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

    #[test]
    fn g0_three_point_plate_support_resolves_geometry_and_builds_physical_modes() {
        let point = QueryPoint::new().with("T", 293.15).unwrap();
        let material = resolve_isotropic_elastic_state_point(
            &material_card(70.0e9, 2_500.0),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let support = RectangularPlateSupport::ThreePointPinned {
            points_centered_m: [[-0.0675, -0.0675], [0.0675, -0.0675], [0.0, 0.0675]],
            maximum_snap_distance_m: 1.0e-12,
        };
        let basis = with_cx(|cx| {
            build_rectangular_plate_modal_basis(
                &RectangularPlateModeRequest {
                    width_m: 0.18,
                    depth_m: 0.18,
                    thickness_m: 0.010,
                    elastic: &material,
                    support,
                    cells_x: 16,
                    cells_y: 16,
                    maximum_nodes: 2_048,
                    minimum_frequency_hz: 10.0,
                    maximum_frequency_hz: 5_000.0,
                    maximum_modes: 64,
                    slice: SliceOptions::default(),
                },
                cx,
            )
            .unwrap()
        });
        assert_eq!(basis.support, support);
        assert_eq!(basis.support_node_indices, vec![36, 48, 246]);
        assert!(basis.maximum_support_snap_error_m <= 1.0e-12);
        assert!(!basis.modes.is_empty());
        assert!(basis.modes.iter().all(|mode| {
            mode.frequency_hz >= 10.0
                && mode.frequency_hz <= 5_000.0
                && mode.frequency_hz.is_finite()
        }));
    }

    #[test]
    fn g0_modal_base_port_moves_the_actual_contact_and_closes_energy() {
        use std::sync::Arc;

        use crate::modal_base_response::{
            RectangularModalBaseIdentity, RectangularModalBasePort, RectangularModalBaseStepInput,
        };

        let point = QueryPoint::new().with("T", 293.15).unwrap();
        let material = resolve_isotropic_elastic_state_point(
            &material_card(70.0e9, 2_500.0),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let support = RectangularPlateSupport::ThreePointPinned {
            points_centered_m: [[-0.0675, -0.0675], [0.0675, -0.0675], [0.0, 0.0675]],
            maximum_snap_distance_m: 1.0e-12,
        };
        let basis = with_cx(|cx| {
            build_rectangular_plate_modal_basis(
                &RectangularPlateModeRequest {
                    width_m: 0.18,
                    depth_m: 0.18,
                    thickness_m: 0.010,
                    elastic: &material,
                    support,
                    cells_x: 8,
                    cells_y: 8,
                    maximum_nodes: 128,
                    minimum_frequency_hz: 10.0,
                    maximum_frequency_hz: 5_000.0,
                    maximum_modes: 16,
                    slice: SliceOptions::default(),
                },
                cx,
            )
            .unwrap()
        });
        let port = RectangularModalBasePort::try_new(
            RectangularModalBaseIdentity {
                model_id: "test/rectangular-modal-base".into(),
                configuration_id: "test/moving-contact".into(),
            },
            Arc::new(basis),
            RayleighDamping::new(0.15, 2.0e-7).unwrap(),
            48_000,
            ModalAcousticTimeBudget::audible_reference(),
            2,
            1.0e-12,
        )
        .unwrap();
        let start = [0.0, 0.0, 0.0];
        let checkpoint = port.initial_static_contact_checkpoint(start, 2.0).unwrap();
        let initial_surface = port.surface_state(&checkpoint, start).unwrap();
        assert!(initial_surface.displacement_m < 0.0);
        assert_eq!(
            initial_surface.velocity_m_per_s.to_bits(),
            0.0_f64.to_bits()
        );

        let proposal = port
            .propose(
                &checkpoint,
                &RectangularModalBaseStepInput {
                    step_id: "moving-contact".into(),
                    expected_version: 0,
                    duration_s: 1.0 / 48_000.0,
                    contact_point_start_base_m: start,
                    contact_point_force_base_m: [0.001, 0.0, 0.0],
                    contact_point_end_base_m: [0.002, 0.0, 0.0],
                    compressive_normal_force_on_base_n: 2.0,
                },
            )
            .unwrap();
        assert!(proposal.receipt().surface_end.displacement_m.is_finite());
        assert!(proposal.receipt().energy_closure_residual_j.abs() <= 1.0e-12);
        let next = port.accept(&checkpoint, proposal).unwrap();
        assert_eq!(next.accepted_version(), 1);
        assert_ne!(
            next.accepted_step_lineage_root(),
            checkpoint.accepted_step_lineage_root()
        );
        assert!(
            port.propose(
                &next,
                &RectangularModalBaseStepInput {
                    step_id: "stale".into(),
                    expected_version: 0,
                    duration_s: 1.0 / 48_000.0,
                    contact_point_start_base_m: [0.002, 0.0, 0.0],
                    contact_point_force_base_m: [0.003, 0.0, 0.0],
                    contact_point_end_base_m: [0.004, 0.0, 0.0],
                    compressive_normal_force_on_base_n: 2.0,
                }
            )
            .is_err()
        );
    }

    #[test]
    fn g0_three_point_plate_support_refuses_degenerate_resolved_constraints() {
        let duplicate = RectangularPlateSupport::ThreePointPinned {
            points_centered_m: [[0.0, 0.0], [1.0e-9, 0.0], [0.04, 0.04]],
            maximum_snap_distance_m: 0.01,
        };
        assert!(
            duplicate
                .validate_for_rectangular_grid(0.18, 0.18, 16, 16)
                .is_err()
        );

        let collinear = RectangularPlateSupport::ThreePointPinned {
            points_centered_m: [[-0.045, 0.0], [0.0, 0.0], [0.045, 0.0]],
            maximum_snap_distance_m: 1.0e-12,
        };
        assert!(
            collinear
                .validate_for_rectangular_grid(0.18, 0.18, 16, 16)
                .is_err()
        );
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
        eprintln!(
            "coarse structural mode frequencies_hz={:?}",
            basis
                .modes
                .iter()
                .map(|mode| mode.frequency_hz)
                .collect::<Vec<_>>()
        );
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
        let loss =
            modal_loss_spectrum_from_rayleigh(&basis, &specimen, rayleigh, ContentHash([0x48; 32]))
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
    fn g0_g3_residual_flexibility_retains_stiff_body_response_above_empty_forcing_band() {
        let specimen = specimen();
        let mut request = coarse_modal_request(&specimen);
        request.maximum_frequency_hz = 12_000.0;
        assert!(matches!(
            with_cx(|cx| build_structural_modal_basis(&request, cx)),
            Err(StructuralModalBasisError::NoModesInBand)
        ));

        let build = |maximum_enrichment_frequency_hz| {
            with_cx(|cx| {
                build_structural_residual_flexibility_estimate_basis(
                    &request,
                    StructuralResidualFlexibilityControls {
                        maximum_enrichment_frequency_hz,
                        maximum_enrichment_modes: 64,
                    },
                    cx,
                )
                .unwrap()
            })
        };
        let mut fine = build(100_000.0);
        assert_eq!(
            fine.authority,
            StructuralResidualFlexibilityAuthority::EstimateOnly
        );
        assert_eq!(fine.authority.code(), "estimate-only");
        assert_eq!(fine.recomputed_identity(), fine.identity);
        let original_budget = fine.assembly_budget.maximum_free_dofs;
        fine.assembly_budget.maximum_free_dofs ^= 1;
        assert_ne!(fine.recomputed_identity(), fine.identity);
        fine.assembly_budget.maximum_free_dofs = original_budget;
        let original_rigid = fine.rigid_modes_per_sqrt_kg[0][0];
        fine.rigid_modes_per_sqrt_kg[0][0] = f64::from_bits(original_rigid.to_bits() ^ 1);
        assert_ne!(fine.recomputed_identity(), fine.identity);
        fine.rigid_modes_per_sqrt_kg[0][0] = original_rigid;
        assert_eq!(fine.rigid_modes_per_sqrt_kg.len(), 6);
        assert_eq!(
            fine.certified_enrichment_mode_count,
            fine.enrichment_modes.len()
        );
        assert_eq!(
            fine.certified_partition_mode_count,
            fine.certified_in_band_mode_count + fine.certified_enrichment_mode_count
        );
        assert!(fine.maximum_rigid_stiffness_relative_residual <= 1.0e-10);
        assert!(fine.maximum_mass_orthogonality_error <= 1.0e-8);
        let cutoff_s2 = (core::f64::consts::TAU * request.maximum_frequency_hz).powi(2);
        assert!(fine.enrichment_modes.iter().all(|mode| {
            mode.eigenvalue_interval_s2.0 > cutoff_s2
                && mode.frequency_hz <= fine.enrichment_frequency_band_hz.1
        }));

        let point = [0.038, 0.0, 0.0];
        let mut unit = fine
            .evaluate_point_force(point, [0.0, 0.0, 1.0], 1.0e-12)
            .unwrap();
        let doubled = fine
            .evaluate_point_force(point, [0.0, 0.0, 2.0], 1.0e-12)
            .unwrap();
        assert!(unit.elastic_work_j > 0.0);
        assert!(unit.recoverable_strain_energy_j > 0.0);
        assert!(unit.maximum_rigid_force_relative_residual <= 1.0e-10);
        assert!(unit.energy_closure_residual_j.abs() <= 1.0e-8 * unit.elastic_work_j);
        assert_eq!(
            unit.energy_closure_residual_j.to_bits(),
            (unit.elastic_work_j - 2.0 * unit.recoverable_strain_energy_j).to_bits()
        );
        assert_eq!(unit.authority, fine.authority);
        assert_eq!(unit.recomputed_identity(), unit.identity);
        let original_panel = unit.panel_normal_displacement_m[0];
        unit.panel_normal_displacement_m[0] = f64::from_bits(original_panel.to_bits() ^ 1);
        assert_ne!(unit.recomputed_identity(), unit.identity);
        unit.panel_normal_displacement_m[0] = original_panel;
        assert_eq!(unit.nodal_displacement_m.len(), fine.mesh.nodes_m.len());
        assert_eq!(
            unit.panel_normal_displacement_m.len(),
            fine.mesh.boundary.triangles.len()
        );
        assert!(
            unit.panel_normal_displacement_m
                .iter()
                .any(|value| *value != 0.0)
        );
        assert!((doubled.elastic_work_j / unit.elastic_work_j - 4.0).abs() < 1.0e-12);
        assert_ne!(
            unit.identity,
            fine.evaluate_point_force([0.0, 0.038, 0.0], [0.0, 0.0, 1.0], 1.0e-12)
                .unwrap()
                .identity
        );
        // The stored projection is the inertia-relieved generalized force, so
        // exact replay must go through the admitted relieved path, not the
        // bare geometric projection: the two arithmetic paths differ by the
        // rigid-relief residue and by summation order (and may diverge across
        // ISAs because only one is written in terms of fused multiply-adds).
        let (relieved_replay, _, _) =
            inertia_relieved_point_force(&fine, point, [0.0, 0.0, 1.0], 1.0e-12).unwrap();
        assert_eq!(
            relieved_replay.modal_force_n_per_sqrt_kg,
            unit.force_projection.modal_force_n_per_sqrt_kg
        );
        let geometric = project_point_force_on_modes(
            &fine.mesh,
            &fine.enrichment_modes,
            point,
            [0.0, 0.0, 1.0],
            1.0e-12,
        )
        .unwrap();
        for (geometric_force, relieved_force) in geometric
            .modal_force_n_per_sqrt_kg
            .iter()
            .zip(&unit.force_projection.modal_force_n_per_sqrt_kg)
        {
            // Inertia relief may move an elastic generalized force only by
            // the mass-orthogonality residue of this basis, already gated to
            // 1.0e-8 during construction above; the bare geometric projection
            // must therefore agree inside a band well below that gate.
            let scale = geometric_force.abs().max(relieved_force.abs());
            assert!(
                (relieved_force - geometric_force).abs() <= 1.0e-9 * scale.max(f64::MIN_POSITIVE)
            );
        }

        let split = fine
            .enrichment_modes
            .windows(2)
            .position(|modes| modes[1].frequency_hz > modes[0].frequency_hz * (1.0 + 1.0e-10))
            .map(|index| index + 1)
            .expect("the coarse specimen must expose at least two distinct elastic frequencies");
        let coarse_cutoff_hz = 0.5
            * (fine.enrichment_modes[split - 1].frequency_hz
                + fine.enrichment_modes[split].frequency_hz);
        let coarse = build(coarse_cutoff_hz);
        assert!(coarse.enrichment_modes.len() < fine.enrichment_modes.len());
        assert_eq!(coarse.operator_identity, fine.operator_identity);
        let comparison = compare_structural_residual_flexibility_estimates(
            &coarse,
            &fine,
            point,
            [0.0, 0.0, 1.0],
            1.0e-12,
        )
        .unwrap();
        assert!(comparison.elastic_work_increment_j >= -1.0e-12 * unit.elastic_work_j);
        for metric in [
            comparison.relative_elastic_work_difference,
            comparison.relative_nodal_displacement_l2_difference,
            comparison.relative_panel_normal_l2_difference,
        ] {
            assert!(metric.is_finite());
        }
    }

    #[test]
    fn g0_g3_broadband_acceleration_sign_and_static_closure() {
        let surface = SpherePanels::icosphere(1.0, 1).unwrap();
        let mode = StructuralMode {
            eigenvalue_s2: 16.0,
            angular_frequency_rad_s: 4.0,
            frequency_hz: 4.0 / core::f64::consts::TAU,
            eigenvalue_interval_s2: (15.0, 17.0),
            eigenvalue_residual_s2: 0.0,
            nodal_shape_per_sqrt_kg: Vec::new(),
            panel_normal_shape_per_sqrt_kg: vec![2.0; surface.centroids().len()],
        };
        let (_, solved) = solve_modal_acceleration_radiation_batch(
            &surface,
            core::slice::from_ref(&mode),
            10.0,
            Medium::air(),
            1.0,
        )
        .unwrap();
        assert!(solved[0].velocity.iter().all(|v| *v == C64::new(0.0, 0.2)));

        let mut acceleration = [f64::NAN];
        let state = fs_couple::modal_acoustic_time::ModalAcousticState {
            displacement_m_sqrt_kg: 0.5,
            velocity_m_sqrt_kg_per_s: -0.25,
        };
        write_closing_modal_acceleration(
            &[mode.clone()],
            &[0.1],
            &[state],
            &[3.0],
            &mut acceleration,
        )
        .unwrap();
        assert!((acceleration[0] - (3.0 - 2.0 * 0.1 * 4.0 * -0.25 - 16.0 * 0.5)).abs() < 1.0e-15);
        let static_state = fs_couple::modal_acoustic_time::ModalAcousticState {
            displacement_m_sqrt_kg: 3.0 / 16.0,
            velocity_m_sqrt_kg_per_s: 0.0,
        };
        write_closing_modal_acceleration(
            &[mode],
            &[0.1],
            &[static_state],
            &[3.0],
            &mut acceleration,
        )
        .unwrap();
        assert_eq!(acceleration[0].to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn g0_g3_retarded_root_lanczos_clock_and_inverse_radius_are_analytic() {
        let arrival_s = 3.0;
        let sound_speed_m_s = 4.0;
        let emission_s =
            bisect_retarded_emission_time(arrival_s, 0.0, 4.0, 1.0e-13, 64, |time_s| {
                Ok(time_s + (10.0 - 0.25 * time_s) / sound_speed_m_s)
            })
            .unwrap();
        let expected_emission_s =
            (arrival_s - 10.0 / sound_speed_m_s) / (1.0 - 0.25 / sound_speed_m_s);
        let radius_m = 10.0 - 0.25 * emission_s;
        assert!((emission_s - expected_emission_s).abs() < 2.0e-13);
        assert!((emission_s + radius_m / sound_speed_m_s - arrival_s).abs() < 2.0e-13);
        let stationary_emission =
            bisect_retarded_emission_time(arrival_s, 0.0, 4.0, 1.0e-13, 64, |time_s| {
                Ok(time_s + 10.0 / sound_speed_m_s)
            })
            .unwrap();
        assert!((stationary_emission - 0.5).abs() < 2.0e-13);

        let stem = StructuralBroadbandSourceStem {
            start_time_s: 0.0,
            sample_rate_hz: 16,
            channels: vec![RealTesseralChannel { l: 0, signed_m: 0 }],
            coefficients: (0..64)
                .map(|frame| FarFieldSourceCoefficientPaM((frame + 1) as f64))
                .collect(),
            authority: BroadbandRadiationAuthority::EstimateOnly,
            source_identity: ContentHash([1; 32]),
            structural_basis_identity: ContentHash([2; 32]),
        };
        let mut out = [0.0];
        interpolate_stem_lanczos8(&stem, 17.0 / 16.0, &mut out).unwrap();
        assert!((out[0] - 17.0).abs() < 1.0e-13);
        let y00 = 0.5 / det::sqrt(core::f64::consts::PI);
        assert!(((out[0] * y00 / radius_m) * radius_m - out[0] * y00).abs() < 1.0e-13);
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
