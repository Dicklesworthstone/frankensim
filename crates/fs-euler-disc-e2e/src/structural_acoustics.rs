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
    Formulation as HelmholtzFormulation, HelmholtzError, Medium, far_field, solve_radiation,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::{ContentHash, DomainHasher};
use fs_couple::modal_acoustic_time::{
    ModalAcousticFrame, ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel,
};
use fs_exec::Cx;
use fs_material::gas::GasState;
use fs_material::visco::{LoweredModel, ViscoError, loss_factor_to_zeta};
use fs_math::c64::C64;
use fs_math::det;
use fs_mesh::{
    RoundedCylinderMeshError, RoundedCylinderMeshSpec, RoundedCylinderTetMesh,
    rounded_cylinder_tet_mesh,
};
use fs_modal::{ModalError, SliceOptions, SliceStats, slice_window};
use fs_rep_frep::SquatDiscEdgeTreatment;
use fs_solid::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial,
    TetLinearElasticProblem, TetMaterialField,
};

use crate::specimen::{DiscProfileSpec, ResolvedMaterialDiscProfile};

/// Schema version of the integrated structural modal artifact.
pub const STRUCTURAL_MODAL_BASIS_SCHEMA_VERSION: u32 = 1;
const STRUCTURAL_MODAL_BASIS_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.structural-modal-basis.v1";
const MODAL_ACOUSTIC_RADIATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.modal-acoustic-radiation.v1";

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
    pub specimen: &'a ResolvedMaterialDiscProfile,
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

/// Exterior radiation of one mass-normalized structural mode at its natural
/// frequency.
#[derive(Clone, Debug)]
pub struct AcousticModeRadiation {
    /// Index into [`StructuralModalBasis::modes`].
    pub structural_mode: usize,
    /// Evaluation angular frequency [rad/s].
    pub angular_frequency_rad_s: f64,
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

/// Integrated physical runtime: one structural basis, state-dependent modal
/// damping, and one BEM-derived observer transfer.
pub struct PhysicalModalAudioModel<'basis> {
    basis: &'basis StructuralModalBasis,
    runtime: ModalAcousticTimeModel,
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
    /// Two independently produced physical artifacts do not share an exact
    /// structural basis or material state.
    IdentityMismatch {
        /// Failed identity relationship.
        what: &'static str,
    },
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
            Self::IdentityMismatch { what } => {
                write!(formatter, "FS-EULER-STRUCTURAL-IDENTITY: {what}")
            }
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
    let material = TetElasticMaterial::from_resolved_state(&request.specimen.material);
    let assembly = TetLinearElasticProblem {
        nodes_m: &mesh.nodes_m,
        tetrahedra: &mesh.tetrahedra,
        materials: TetMaterialField::Uniform(&material),
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
        material_state_identity: request.specimen.material.resolved().identity(),
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
    specimen: &ResolvedMaterialDiscProfile,
    model: &LoweredModel,
    damping_model_identity: ContentHash,
) -> Result<ModalLossSpectrum, StructuralModalBasisError> {
    if basis.specimen_identity != specimen.identity
        || basis.material_state_identity != specimen.material.resolved().identity()
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
            radiation_identity: radiation.identity,
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
            let solution = solve_radiation(
                &surface,
                k,
                acoustic_medium,
                &velocity,
                HelmholtzFormulation::BurtonMiller,
            )?;
            let velocity_area_norm = solution
                .velocity
                .iter()
                .zip(surface.areas())
                .map(|(value, area)| value.norm_sq() * area)
                .sum::<f64>();
            let plane_wave_power_scale =
                0.5 * acoustic_medium.density * acoustic_medium.sound_speed * velocity_area_norm;
            let power_tolerance = 1.0e-11 * plane_wave_power_scale.max(f64::MIN_POSITIVE);
            if solution.radiated_power < -power_tolerance {
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
    specimen: &ResolvedMaterialDiscProfile,
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
        MaterialPropertySelection, resolve_isotropic_solid_state_point,
    };
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

    fn specimen() -> ResolvedMaterialDiscProfile {
        let point = QueryPoint::new().with("T", 293.15).unwrap();
        let material = resolve_isotropic_solid_state_point(
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
            .resolve_with_material_state(&material, cx)
            .unwrap()
        })
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
}
