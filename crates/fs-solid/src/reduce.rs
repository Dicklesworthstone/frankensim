//! Offline lip/reed reduce lab (music bead
//! `frankensim-music-v8-root-3ez8g.3.3`): a 3-D leaflet mesh + a
//! material card become a VALVE CARD — reduced masses, stiffnesses,
//! damping, and rest-orifice geometry measured from the mesh — with the
//! provenance chain and the retained-compliance disclosure that make
//! "that reed" different from "a plausible reed". The mesh itself never
//! runs at audio rate; the card is what
//! `fs_phs::mass_spring_damper` + `bernoulli_volume_flow` islands (and
//! `fs_scenario::BeatingReed`) consume.
//!
//! REDUCTION RECIPE (no new FEM — D23): [`crate::TetLinearElasticProblem`]
//! assembles the ordinary `(K, M)` pencil; [`fs_modal::slice_window`]
//! certifies the low modes; each retained mode is RE-NORMALIZED from the
//! solver's `phi^T M phi = 1` convention to UNIT MEAN TRANSLATION of the
//! orifice face along the opening axis (mass normalization would destroy
//! the units — the recorded fs-solid lesson), giving physical
//! `m_eff = phi^T M phi` [kg] and `k_eff = lambda m_eff` [N/m].
//!
//! RETAINED-COMPLIANCE HONESTY: the modal compliance
//! `sum Gamma_n^2 / lambda_n` of the retained modes is compared against
//! the EXACT static compliance `r^T K^{-1} r` (sparse LDLT solve) for
//! the face-opening load pattern; a reduction whose retained fraction
//! falls below the authored floor REFUSES rather than shipping a stiff
//! lie.
//!
//! DAMPING HONESTY: compatibility reduction uses an explicitly authored loss
//! factor, with `c = eta sqrt(k m)`. The [`visco`] path projects a selected
//! tensor-Prony/Kelvin-Voigt material law on the SAME modes, with explicit
//! frozen-frequency and coupled-damping admission. The per-mode damping
//! coefficient is authoritative for runtime lowering and survives the codec.
//! No measured quality factor is inferred from a material name.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::state_point::IsotropicElasticStatePoint;
use fs_modal::{SliceOptions, slice_window};
use std::fmt::Write as _;

use crate::linear3::tet10::{QuadraticTetProblem, TetQuadrature};
use crate::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
};

/// Anisotropic and spatially varying material fields through the same reducer.
pub mod material_field;
/// Thermal equilibrium and prestress through the same pressure-loaded reducer.
pub mod thermal;
/// Source tensor-Prony/Kelvin-Voigt loss through the common valve reduction.
pub mod visco;
mod pressure_compliance;
pub use pressure_compliance::{PressureCompliance, SingleModeValveAdmission};

/// Domain-separated identity for valve cards.
pub const VALVE_CARD_HASH_DOMAIN: &str = "org.frankensim.fs-solid.valve-card.v1";
/// Card schema version (bumps refuse old bytes, never reinterpret).
pub const VALVE_CARD_SCHEMA_VERSION: u32 = 2;

/// Required discretization identity, retained by the card and its codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValveDiscretization {
    /// Existing affine four-node path, with its original assembly arithmetic.
    Tet4,
    /// Quadratic isoparametric path, with explicit quadrature.
    Tet10(TetQuadrature),
}

/// Required request choice; quadratic connectivity is supplied, never guessed.
#[derive(Debug, Clone, Copy)]
pub enum ValveElementOrder<'a> {
    /// Existing P1 connectivity in `ValveCardRequest::tetrahedra`.
    Tet4,
    /// Full connectivity in canonical TET10 order. Its corner sets must agree
    /// with the request's per-element P1 topology; midsides may be curved.
    Tet10 {
        /// Actual quadratic topology, including the pressure/observation nodes.
        elements: &'a [[usize; 10]],
        /// Explicit integration choice.
        quadrature: TetQuadrature,
    },
}

impl ValveElementOrder<'_> {
    fn discretization(self) -> ValveDiscretization {
        match self {
            Self::Tet4 => ValveDiscretization::Tet4,
            Self::Tet10 { quadrature, .. } => ValveDiscretization::Tet10(quadrature),
        }
    }
}

/// Typed refusals from the reduce lab.
#[derive(Debug)]
pub enum ReduceError {
    /// A request parameter is unusable.
    Invalid {
        /// Diagnosis.
        what: &'static str,
    },
    /// The elasticity assembly refused.
    Assembly(crate::linear3::TetElasticError),
    /// The source-bound viscoelastic element or frequency query refused.
    Visco(crate::linear3::visco::ViscoElementError),
    /// The modal solve refused.
    Modal(fs_modal::ModalError),
    /// The static-compliance factorization refused.
    Factor(fs_sparse::LdltError),
    /// The retained modes carry too little of the opening compliance.
    RetainedComplianceTooLow {
        /// Measured retained fraction.
        fraction: f64,
        /// The authored floor.
        floor: f64,
    },
    /// Card bytes failed decode/verify.
    Card {
        /// Diagnosis.
        what: &'static str,
    },
}

impl core::fmt::Display for ReduceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ReduceError::Invalid { what } => write!(f, "FS-SOLID-REDUCE: {what}"),
            ReduceError::Assembly(e) => write!(f, "FS-SOLID-REDUCE assembly: {e:?}"),
            ReduceError::Visco(e) => write!(f, "FS-SOLID-REDUCE visco: {e}"),
            ReduceError::Modal(e) => write!(f, "FS-SOLID-REDUCE modal: {e:?}"),
            ReduceError::Factor(e) => write!(f, "FS-SOLID-REDUCE factor: {e:?}"),
            ReduceError::RetainedComplianceTooLow { fraction, floor } => write!(
                f,
                "FS-SOLID-REDUCE: retained modes carry {fraction:.4} of the opening \
                 compliance, below the authored floor {floor:.4} — a reduction that \
                 discards the compliance refuses rather than shipping a stiff lie"
            ),
            ReduceError::Card { what } => write!(f, "FS-SOLID-REDUCE card: {what}"),
        }
    }
}

impl core::error::Error for ReduceError {}

impl From<crate::linear3::TetElasticError> for ReduceError {
    fn from(e: crate::linear3::TetElasticError) -> Self {
        ReduceError::Assembly(e)
    }
}

impl From<crate::linear3::visco::ViscoElementError> for ReduceError {
    fn from(e: crate::linear3::visco::ViscoElementError) -> Self {
        Self::Visco(e)
    }
}

impl From<fs_modal::ModalError> for ReduceError {
    fn from(e: fs_modal::ModalError) -> Self {
        ReduceError::Modal(e)
    }
}

impl From<fs_sparse::LdltError> for ReduceError {
    fn from(e: fs_sparse::LdltError) -> Self {
        ReduceError::Factor(e)
    }
}

/// The orifice: which mesh nodes form the moving valve face, which axis
/// opens, and where the opposing surface sits.
#[derive(Debug, Clone)]
pub struct OrificeSpec<'a> {
    /// Node indices of the moving face (the free lip/reed edge surface).
    pub face_nodes: &'a [usize],
    /// Unit opening axis (the face moves along this to close the gap).
    pub opening_axis: [f64; 3],
    /// Unit width axis (the slit's long direction), perpendicular to the
    /// opening axis within 1e-12 in normalized dot product.
    pub width_axis: [f64; 3],
    /// Offset of the opposing surface along the opening axis [m]; the
    /// rest gap is measured FROM THE MESH as the CLOSEST approach
    /// `min(face . axis) - offset`. Face nodes are the surface strip
    /// FACING the opposing plane (the underside of a lip/reed near its
    /// free edge); the third axis (opening x width) is the flow
    /// direction, whose face extent is the channel's effective
    /// thickness.
    pub opposing_plane_offset_m: f64,
}

/// A reduce-lab request.
pub struct ValveCardRequest<'a, M = &'a IsotropicElasticStatePoint> {
    /// Leaflet vertex coordinates [m].
    pub nodes_m: &'a [[f64; 3]],
    /// Conforming P1 tetrahedra.
    pub tetrahedra: &'a [[usize; 4]],
    /// Required element order; no default or implicit promotion.
    pub element_order: ValveElementOrder<'a>,
    /// Strongly fixed DOFs (`3*node + component`) — the clamped root.
    pub fixed_dofs: &'a [usize],
    /// Resolved isotropic state on the compatibility path, or the common
    /// tetrahedral material field with [`material_field`]. Geometry, loads,
    /// constraints and reduction controls do not depend on this choice.
    pub material: M,
    /// The orifice geometry spec.
    pub orifice: OrificeSpec<'a>,
    /// Modal window [Hz] handed to `slice_window`.
    pub window_hz: (f64, f64),
    /// Refuse when the retained modal compliance falls below this
    /// fraction of the exact static compliance.
    pub retained_compliance_floor: f64,
    /// AUTHORED structural loss factor eta (Estimate until a measured
    /// tissue/cane pack exists); `c = eta sqrt(k m)`.
    pub loss_factor: f64,
    /// Caller-supplied source identity (mesh + geometry provenance).
    pub source_id: &'a str,
}

/// One reduced mode of the card.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedMode {
    /// Natural frequency [Hz].
    pub frequency_hz: f64,
    /// Effective mass with the face-unit-translation normalization [kg].
    pub mass_kg: f64,
    /// `lambda * mass` [N/m].
    pub stiffness_n_m: f64,
    /// Modal damping [N s/m]. The legacy path uses `eta * sqrt(k m)`;
    /// the material-loss path projects the source viscosity on this same mode.
    pub damping_n_s_m: f64,
    /// Opening-force projection `Gamma = phi^T r` for the unit face
    /// force pattern (phi M-normalized; N/kg^(1/2) per unit pattern).
    pub participation_sqrt_kg: f64,
}

/// The minted valve card.
#[derive(Debug, Clone, PartialEq)]
pub struct ValveCard {
    /// Content identity over the canonical bytes (domain-separated).
    pub identity: ContentHash,
    /// Caller-supplied source identity.
    pub source_id: String,
    /// Material state receipt identity (fs-material resolution chain).
    pub material_identity: ContentHash,
    /// Mesh digest (nodes + tets, bit-exact).
    pub mesh_digest: ContentHash,
    /// Order and quadrature bound into canonical bytes and card identity.
    pub discretization: ValveDiscretization,
    /// Reduced modes, ascending frequency (index 0 drives the 1-DOF
    /// island; 0..2 drive the two-mass realization).
    pub modes: Vec<ReducedMode>,
    /// Rest gap measured from the mesh [m].
    pub rest_gap_m: f64,
    /// Slit width measured from the mesh [m].
    pub width_m: f64,
    /// Effective thickness measured from the mesh [m].
    pub effective_thickness_m: f64,
    /// Retained modal compliance / exact static compliance.
    pub retained_compliance_fraction: f64,
    /// Authored loss factor on the compatibility path; zero when modal
    /// damping is supplied entirely by the selected source material law.
    /// Runtime lowering always reads each mode's actual `damping_n_s_m`.
    pub loss_factor: f64,
    /// Total leaflet mass [kg].
    pub total_mass_kg: f64,
}

/// A valve reduction with a work-conjugate distributed pressure port.
/// Load vectors describe force per pascal toward opening; the reed's positive
/// closing pressure acts against them. Areas may be signed for higher modes.
#[derive(Debug, Clone)]
pub struct PressureLoadedValve {
    card: ValveCard,
    effective_areas_m2: Vec<f64>,
    source_modes: Vec<fs_modal::ModePair>,
    mass_normalized_pressure_load: Vec<f64>,
    pressure_compliance: PressureCompliance,
    identity: ContentHash,
}

impl PressureLoadedValve {
    /// Structural reduction with its material and mesh identities.
    #[must_use]
    pub fn card(&self) -> &ValveCard {
        &self.card
    }

    /// Generalized force per pascal, also volume flow per modal velocity.
    #[must_use]
    pub fn effective_areas_m2(&self) -> &[f64] {
        &self.effective_areas_m2
    }

    /// Original mass-normalized shapes in card order and in the source
    /// assembly's reduced free-DOF ordering. Sorted mode indices do not identify
    /// the same physical displacement field when the source changes.
    #[must_use]
    pub fn source_modes(&self) -> &[fs_modal::ModePair] {
        &self.source_modes
    }

    /// Exact `phi^T load_per_pa` [m² / sqrt(kg)] on each source mode.
    /// The transpose maps modal velocity to work-conjugate swept volume flow.
    #[must_use]
    pub fn mass_normalized_pressure_load(&self) -> &[f64] {
        &self.mass_normalized_pressure_load
    }

    /// Identity of the structural card, nodal pressure load and projections.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact static responses and the response carried by the retained modes.
    /// These diagnose truncation of the actual pressure load, independently of
    /// the opening-load compliance gate used to choose structural modes.
    #[must_use]
    pub const fn pressure_compliance(&self) -> &PressureCompliance {
        &self.pressure_compliance
    }

    /// Admit the ACTUAL first-mode image before converting to a playing reed.
    ///
    /// A multi-mode reduction passing its opening-compliance floor does not
    /// establish that mode zero represents either opening or pressure work.
    /// The caller explicitly bounds all three static truncation measures.
    /// This is not a nonlinear or broadband acoustic accuracy certificate.
    ///
    /// # Errors
    /// Refuses inadequate first-mode compliance, an unrepresentable closing
    /// direction, invalid thresholds, or any ordinary reed-conversion refusal.
    pub fn beating_reed_checked(
        &self,
        pressure_pa: f64,
        attack_s: f64,
        admission: SingleModeValveAdmission,
    ) -> Result<fs_scenario::BeatingReed, ReduceError> {
        self.pressure_compliance.admit_single_mode(admission)?;
        self.beating_reed(pressure_pa, attack_s)
    }

    /// Convert the first mode using its projected area, without authored area.
    /// This compatibility path does NOT admit first-mode adequacy. New compiled
    /// physical assemblies should use [`Self::beating_reed_checked`].
    /// # Errors
    /// Refuses nonpositive first-mode area or invalid reed inputs. Signed or
    /// zero ports cannot be represented by this positive-closing reed model.
    pub fn beating_reed(
        &self,
        pressure_pa: f64,
        attack_s: f64,
    ) -> Result<fs_scenario::BeatingReed, ReduceError> {
        self.card
            .beating_reed(self.effective_areas_m2[0], pressure_pa, attack_s)
    }
}

fn unit(v: [f64; 3]) -> Result<[f64; 3], ReduceError> {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if !(n > 0.0 && n.is_finite()) {
        return Err(ReduceError::Invalid {
            what: "axes must be nonzero finite vectors",
        });
    }
    Ok([v[0] / n, v[1] / n, v[2] / n])
}

fn mesh_digest(nodes: &[[f64; 3]], tets: &[[usize; 4]]) -> ContentHash {
    let mut h = DomainHasher::new("org.frankensim.fs-solid.valve-mesh.v1");
    h.update(&(nodes.len() as u64).to_le_bytes());
    for n in nodes {
        for c in n {
            h.update(&c.to_bits().to_le_bytes());
        }
    }
    h.update(&(tets.len() as u64).to_le_bytes());
    for t in tets {
        for &i in t {
            h.update(&(i as u64).to_le_bytes());
        }
    }
    h.finalize()
}

// The face is a set of nodes (validated by the caller). Prescribed zero
// components contribute zero to its mean, but their nodes remain in the
// denominator. Counting free x components would change normalization under
// unrelated constraints or a permutation of coordinate axes.
fn mean_face_translation(
    phi: &[f64],
    free_dofs: &[usize],
    face_nodes: &[usize],
    axis: [f64; 3],
) -> f64 {
    let mut sum = 0.0;
    for (&value, &full) in phi.iter().zip(free_dofs) {
        if face_nodes.contains(&(full / 3)) {
            sum += value * axis[full % 3];
        }
    }
    sum / face_nodes.len() as f64
}

/// Run the reduction: assemble, certify modes, renormalize, disclose the
/// retained compliance, measure the orifice, and mint the card.
///
/// # Errors
/// [`ReduceError`] on any refusal (see the module doc for the
/// retained-compliance law).
#[allow(clippy::too_many_lines)] // one reduction pipeline, kept whole
pub fn reduce_valve(request: &ValveCardRequest<'_>, cx: &Cx<'_>) -> Result<ValveCard, ReduceError> {
    let material = TetElasticMaterial::from_resolved_elastic_state(request.material);
    reduce_valve_impl(
        request,
        Some(TetMaterialField::Uniform(&material)),
        request.material.identity(),
        cx,
        None,
        None,
    )
    .map(|reduced| reduced.card)
}

/// Reduce with explicit nodal force-per-pascal vectors [m²], in mesh axes.
/// Each vector is the consistent nodal integral of the chosen pressure shape
/// over the caller's loaded surface. Fixed DOFs do no virtual work. The load
/// may cover a different surface than the orifice's opening observation.
/// For a unit-mean-opening mode psi, A = psi^T load, so F = p A and Q = A v.
/// This uses the existing structural modes; it does not solve pressure fields
/// or infer pressure-load adequacy from the opening-compliance test. Exact
/// static pressure responses are retained for explicit first-mode admission.
/// # Errors
/// Refuses mismatched/nonfinite load fields, nonfinite projections, or any
/// ordinary reduction refusal. Negative/zero areas are retained, not clamped.
pub fn reduce_pressure_loaded_valve(
    request: &ValveCardRequest<'_>,
    nodal_load_per_pa_m2: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<PressureLoadedValve, ReduceError> {
    let material = TetElasticMaterial::from_resolved_elastic_state(request.material);
    reduce_pressure_loaded_valve_impl(
        request,
        TetMaterialField::Uniform(&material),
        request.material.identity(),
        nodal_load_per_pa_m2,
        cx,
        None,
    )
}

fn reduce_pressure_loaded_valve_impl<M>(
    request: &ValveCardRequest<'_, M>,
    materials: TetMaterialField<'_>,
    material_identity: ContentHash,
    nodal_load_per_pa_m2: &[[f64; 3]],
    cx: &Cx<'_>,
    prepared: Option<PreparedValve<'_>>,
) -> Result<PressureLoadedValve, ReduceError> {
    if nodal_load_per_pa_m2.len() != request.nodes_m.len()
        || !nodal_load_per_pa_m2.iter().flatten().all(|v| v.is_finite())
    {
        return Err(ReduceError::Invalid {
            what: "pressure load must contain one finite vector per mesh node",
        });
    }
    let reduced = reduce_valve_impl(
        request,
        Some(materials),
        material_identity,
        cx,
        Some(nodal_load_per_pa_m2),
        prepared,
    )?;
    pressure_loaded_from_reduction(reduced, nodal_load_per_pa_m2)
}

fn pressure_loaded_from_reduction(
    reduced: ReducedValve,
    nodal_load_per_pa_m2: &[[f64; 3]],
) -> Result<PressureLoadedValve, ReduceError> {
    let card = reduced.card;
    let effective_areas_m2 = reduced.effective_areas_m2;
    let pressure_compliance = reduced.pressure_compliance.ok_or(ReduceError::Invalid {
        what: "pressure-loaded reduction did not retain its static compliance",
    })?;
    let mut hash = DomainHasher::new("org.frankensim.fs-solid.pressure-loaded-valve.v2");
    hash.update(card.identity.as_bytes());
    hash.update(&(nodal_load_per_pa_m2.len() as u64).to_le_bytes());
    for value in nodal_load_per_pa_m2
        .iter()
        .flatten()
        .chain(effective_areas_m2.iter())
    {
        hash.update(&value.to_bits().to_le_bytes());
    }
    pressure_compliance.append_identity(&mut hash);
    Ok(PressureLoadedValve {
        card,
        effective_areas_m2,
        source_modes: reduced.source_modes,
        mass_normalized_pressure_load: reduced.mass_normalized_pressure_load,
        pressure_compliance,
        identity: hash.finalize(),
    })
}

fn validate_element_topology<M>(request: &ValveCardRequest<'_, M>) -> Result<(), ReduceError> {
    if let ValveElementOrder::Tet10 { elements, .. } = request.element_order {
        if elements.len() != request.tetrahedra.len()
            || elements.iter().zip(request.tetrahedra).any(|(q, p)| {
                let mut a = [q[0], q[1], q[2], q[3]];
                let mut b = *p;
                a.sort_unstable();
                b.sort_unstable();
                a != b
            })
        {
            return Err(ReduceError::Invalid {
                what: "TET10 corners must match the request topology element by element",
            });
        }
    }
    Ok(())
}

fn pressure_area(phi: &[f64], free_dofs: &[usize], load: &[[f64; 3]], mean_opening: f64) -> f64 {
    phi.iter()
        .zip(free_dofs)
        .map(|(value, full)| value * load[full / 3][full % 3])
        .sum::<f64>()
        / mean_opening
}

// Prepared operators stay with their source owner. Thermal observation geometry
// and viscoelastic material loss reuse this reducer's modes and pressure ports.
struct PreparedValve<'a> {
    assembly: &'a TetElasticAssembly,
    observation_nodes_m: &'a [[f64; 3]],
    /// Equilibrium K for static compliance when the modal pencil is frozen K'(omega).
    static_stiffness: Option<&'a fs_sparse::Csr>,
    material_loss: Option<&'a visco::Projection<'a>>,
}

struct ReducedValve {
    card: ValveCard,
    effective_areas_m2: Vec<f64>,
    source_modes: Vec<fs_modal::ModePair>,
    mass_normalized_pressure_load: Vec<f64>,
    pressure_compliance: Option<PressureCompliance>,
    material_loss: Option<visco::Approximation>,
    // Only the explicit complete material image owns all selected modes,
    // including modes orthogonal to the mean opening observation.
    complete_material_modes: Option<Vec<fs_modal::ModePair>>,
}

fn reduce_valve_impl<M>(
    request: &ValveCardRequest<'_, M>,
    materials: Option<TetMaterialField<'_>>,
    material_identity: ContentHash,
    cx: &Cx<'_>,
    pressure_load: Option<&[[f64; 3]]>,
    prepared: Option<PreparedValve<'_>>,
) -> Result<ReducedValve, ReduceError> {
    if request.face_nodes_invalid() {
        return Err(ReduceError::Invalid {
            what: "orifice face nodes must be non-empty, unique and in range",
        });
    }
    if !(request.window_hz.0 > 0.0 && request.window_hz.1 > request.window_hz.0) {
        return Err(ReduceError::Invalid {
            what: "modal window needs 0 < lo < hi",
        });
    }
    if !(request.retained_compliance_floor > 0.0 && request.retained_compliance_floor <= 1.0) {
        return Err(ReduceError::Invalid {
            what: "retained-compliance floor must be in (0, 1]",
        });
    }
    if !(request.loss_factor >= 0.0 && request.loss_factor.is_finite()) {
        return Err(ReduceError::Invalid {
            what: "loss factor must be finite and non-negative",
        });
    }
    let axis = unit(request.orifice.opening_axis)?;
    let width_axis = unit(request.orifice.width_axis)?;
    let alignment: f64 = axis.iter().zip(width_axis).map(|(a, w)| a * w).sum();
    if alignment.abs() > 1e-12 {
        return Err(ReduceError::Invalid {
            what: "orifice opening and width axes must be perpendicular",
        });
    }
    if !request.orifice.opposing_plane_offset_m.is_finite() {
        return Err(ReduceError::Invalid {
            what: "orifice opposing-plane offset must be finite",
        });
    }

    validate_element_topology(request)?;
    let thermal_tangent = prepared.as_ref().is_some_and(|p| p.material_loss.is_none());
    let prepared_pencil = prepared.is_some();
    let assembled;
    let (assembly, observation_nodes_m, static_stiffness, material_loss) = if let Some(prepared) = prepared {
        (
            prepared.assembly,
            prepared.observation_nodes_m,
            prepared.static_stiffness.unwrap_or(&prepared.assembly.stiffness),
            prepared.material_loss,
        )
    } else {
        let materials = materials.ok_or(ReduceError::Invalid {
            what: "unprepared valve reduction requires an elastic material field",
        })?;
        // The original unheated path retains its assembly and arithmetic.
        assembled = match request.element_order {
            ValveElementOrder::Tet4 => TetLinearElasticProblem {
                nodes_m: request.nodes_m,
                tetrahedra: request.tetrahedra,
                materials,
                fixed_dofs: request.fixed_dofs,
                budget: TetAssemblyBudget::default(),
            }
            .assemble(cx)?,
            ValveElementOrder::Tet10 { elements, quadrature } => QuadraticTetProblem {
                nodes_m: request.nodes_m,
                elements,
                materials,
                fixed_dofs: request.fixed_dofs,
                budget: TetAssemblyBudget::default(),
                quadrature,
            }
            .assemble(cx)?,
        };
        (&assembled, request.nodes_m, &assembled.stiffness, None)
    };

    // 2. Certified low modes.
    let tau = core::f64::consts::TAU;
    let window = (
        (tau * request.window_hz.0).powi(2),
        (tau * request.window_hz.1).powi(2),
    );
    let report = slice_window(
        &assembly.stiffness,
        &assembly.mass,
        window,
        &SliceOptions::default(),
    )?;
    if report.modes.is_empty() {
        return Err(ReduceError::Invalid {
            what: "no modes in the requested window",
        });
    }

    let retain_complete_material_basis = material_loss.is_some_and(|p| p.complete_basis);
    let material_loss = material_loss
        .map(|projection| projection.project(&report.modes, cx))
        .transpose()?;

    // 3. Opening-direction load pattern on the reduced DOFs: the unit
    //    translation of the FACE nodes along the opening axis.
    let n_red = assembly.free_dofs.len();
    let mut r = vec![0.0f64; n_red];
    let mut face_in_free = 0usize;
    for (red, &full) in assembly.free_dofs.iter().enumerate() {
        let node = full / 3;
        let comp = full % 3;
        if request.orifice.face_nodes.contains(&node) {
            r[red] = axis[comp];
            face_in_free += 1;
        }
    }
    if face_in_free == 0 {
        return Err(ReduceError::Invalid {
            what: "every orifice face node is fixed; the valve cannot move",
        });
    }

    // 4. Exact static compliance r^T K^{-1} r (sparse LDLT).
    let symbolic =
        fs_sparse::SymbolicLdlt::analyze(static_stiffness, fs_sparse::DirectOrdering::Amd)?;
    let factor = symbolic.factor(static_stiffness, &fs_sparse::LdltOptions::default())?;
    // A positive modal window cannot hide an unstable direction elsewhere.
    // Singular tangents already refuse in the existing factorization owner.
    if prepared_pencil
        && (factor.inertia().negative != 0 || factor.inertia().positive != n_red)
    {
        return Err(ReduceError::Invalid {
            what: if thermal_tangent {
                "thermal prestress tangent must be positive definite before valve reduction"
            } else {
                "material valve equilibrium stiffness must be positive definite before reduction"
            },
        });
    }
    let k_inv_r = factor.solve(&r);
    let static_compliance: f64 = r.iter().zip(&k_inv_r).map(|(a, b)| a * b).sum();
    if !(static_compliance > 0.0 && static_compliance.is_finite()) {
        return Err(ReduceError::Invalid {
            what: "static opening compliance is not positive finite",
        });
    }

    // Reuse the SAME factorization for the actual pressure load. Reciprocity
    // gives mean opening per pascal from K^-1 r, with ALL face nodes in the
    // denominator (including fixed nodes). The pressure energy compliance is
    // l^T K^-1 l; an opening-compliance gate cannot stand in for this quantity.
    let exact_pressure = pressure_load.map(|load| {
        let reduced_load: Vec<f64> = assembly
            .free_dofs
            .iter()
            .map(|&full| load[full / 3][full % 3])
            .collect();
        let displacement = factor.solve(&reduced_load);
        let volume: f64 = reduced_load
            .iter()
            .zip(&displacement)
            .map(|(a, b)| a * b)
            .sum();
        let opening: f64 = reduced_load
            .iter()
            .zip(&k_inv_r)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            / request.orifice.face_nodes.len() as f64;
        (volume, opening)
    });

    // 5. Participations + retained-compliance fraction. For M-orthonormal
    //    modes, the EXPORTED image compliance of the FORCE pattern r is
    //    sum (phi^T r)^2 / lambda. For the material-visco path its frozen
    //    K differs from the equilibrium K used above: this remains the actual
    //    image/source comparison, not an exact equilibrium spectral partition.
    //    The projection is phi^T r directly (phi^T M r is the
    //    base-excitation participation, a different animal; the first
    //    build used it and retained 1e-11 of the compliance).
    let mut modes = Vec::new();
    let mut effective_areas_m2 = Vec::new();
    let mut source_modes = Vec::new();
    let mut mass_normalized_pressure_load = Vec::new();
    let mut retained = 0.0f64;
    for (mode_index, pair) in report.modes.iter().enumerate() {
        let gamma: f64 = pair.phi.iter().zip(&r).map(|(a, b)| a * b).sum();
        retained += gamma * gamma / pair.lambda;
        let u_face = mean_face_translation(
            &pair.phi,
            &assembly.free_dofs,
            request.orifice.face_nodes,
            axis,
        );
        if u_face.abs() < 1e-30 {
            // A mode that does not move the face contributes no valve
            // dynamics; keep it out of the card's realizations.
            continue;
        }
        let scale = 1.0 / u_face;
        // phi_hat = phi * scale; m_eff = phi_hat^T M phi_hat = scale^2
        // (solver convention phi^T M phi = 1).
        let mass_kg = scale * scale;
        let stiffness_n_m = pair.lambda * mass_kg;
        let damping_n_s_m = if let Some(loss) = &material_loss {
            loss.modal_damping_per_s[mode_index] * mass_kg
        } else {
            request.loss_factor * (stiffness_n_m * mass_kg).sqrt()
        };
        if let Some(load) = pressure_load {
            let modal_load = pressure_area(&pair.phi, &assembly.free_dofs, load, 1.0);
            let area = modal_load / u_face;
            if !area.is_finite() {
                return Err(ReduceError::Invalid {
                    what: "projected pressure area is not finite",
                });
            }
            effective_areas_m2.push(area);
            mass_normalized_pressure_load.push(modal_load);
        }
        source_modes.push(pair.clone());
        modes.push(ReducedMode {
            frequency_hz: pair.lambda.sqrt() / tau,
            mass_kg,
            stiffness_n_m,
            damping_n_s_m,
            participation_sqrt_kg: gamma,
        });
    }
    let fraction = retained / static_compliance;
    if fraction < request.retained_compliance_floor {
        return Err(ReduceError::RetainedComplianceTooLow {
            fraction,
            floor: request.retained_compliance_floor,
        });
    }
    if modes.is_empty() {
        return Err(ReduceError::Invalid {
            what: "no retained mode moves the orifice face",
        });
    }

    // 6. Orifice geometry measured from the accepted observation mesh.
    let thickness_axis = [
        axis[1] * width_axis[2] - axis[2] * width_axis[1],
        axis[2] * width_axis[0] - axis[0] * width_axis[2],
        axis[0] * width_axis[1] - axis[1] * width_axis[0],
    ];
    let mut gap_min = f64::INFINITY;
    let mut w_lo = f64::INFINITY;
    let mut w_hi = f64::NEG_INFINITY;
    let mut t_lo = f64::INFINITY;
    let mut t_hi = f64::NEG_INFINITY;
    for &node in request.orifice.face_nodes {
        let p = observation_nodes_m[node];
        let a = p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2];
        let w = p[0] * width_axis[0] + p[1] * width_axis[1] + p[2] * width_axis[2];
        let t = p[0] * thickness_axis[0] + p[1] * thickness_axis[1] + p[2] * thickness_axis[2];
        gap_min = gap_min.min(a - request.orifice.opposing_plane_offset_m);
        w_lo = w_lo.min(w);
        w_hi = w_hi.max(w);
        t_lo = t_lo.min(t);
        t_hi = t_hi.max(t);
    }
    let rest_gap_m = gap_min;
    if !(rest_gap_m > 0.0) {
        return Err(ReduceError::Invalid {
            what: "measured rest gap is not positive (face behind the opposing plane)",
        });
    }
    let width_m = w_hi - w_lo;
    let effective_thickness_m = t_hi - t_lo;

    let mut card = ValveCard {
        identity: ContentHash([0u8; 32]),
        source_id: request.source_id.to_string(),
        material_identity,
        mesh_digest: {
            let base = mesh_digest(request.nodes_m, request.tetrahedra);
            match request.element_order {
                ValveElementOrder::Tet4 => base,
                ValveElementOrder::Tet10 { elements, .. } => {
                    let mut h = DomainHasher::new("org.frankensim.fs-solid.valve-mesh.tet10.v1");
                    h.update(base.as_bytes());
                    for tet in elements {
                        for node in tet {
                            h.update(&(*node as u64).to_le_bytes());
                        }
                    }
                    h.finalize()
                }
            }
        },
        discretization: request.element_order.discretization(),
        modes,
        rest_gap_m,
        width_m,
        effective_thickness_m,
        retained_compliance_fraction: fraction,
        loss_factor: request.loss_factor,
        total_mass_kg: assembly.total_mass_kg,
    };
    card.identity = card.recomputed_identity();
    let pressure_compliance = exact_pressure
        .map(|(volume, opening)| {
            let face_count = request.orifice.face_nodes.len() as f64;
            PressureCompliance::from_modes(
                &card.modes,
                &effective_areas_m2,
                static_compliance / face_count / face_count,
                volume,
                opening,
            )
        })
        .transpose()?;
    Ok(ReducedValve {
        card,
        effective_areas_m2,
        source_modes,
        mass_normalized_pressure_load,
        pressure_compliance,
        material_loss,
        complete_material_modes: retain_complete_material_basis.then(|| report.modes.clone()),
    })
}

impl<M> ValveCardRequest<'_, M> {
    fn face_nodes_invalid(&self) -> bool {
        self.orifice.face_nodes.is_empty()
            || self
                .orifice
                .face_nodes
                .iter()
                .enumerate()
                .any(|(i, &n)| n >= self.nodes_m.len() || self.orifice.face_nodes[..i].contains(&n))
    }
}

impl ValveCard {
    /// The card's identity over its canonical fields (schema-versioned,
    /// domain-separated; upstream identities first, then f64 bits, then
    /// length-prefixed mode payload).
    #[must_use]
    pub fn recomputed_identity(&self) -> ContentHash {
        let mut h = DomainHasher::new(VALVE_CARD_HASH_DOMAIN);
        h.update(&VALVE_CARD_SCHEMA_VERSION.to_le_bytes());
        h.update(self.material_identity.as_bytes());
        h.update(self.mesh_digest.as_bytes());
        match self.discretization {
            ValveDiscretization::Tet4 => h.update(&[1]),
            ValveDiscretization::Tet10(q) => {
                h.update(&[2]);
                h.update(&(q.stiffness_points_per_dir as u64).to_le_bytes());
                h.update(&(q.mass_points_per_dir as u64).to_le_bytes());
            }
        }
        h.update(self.source_id.as_bytes());
        h.update(&[0u8]);
        for v in [
            self.rest_gap_m,
            self.width_m,
            self.effective_thickness_m,
            self.retained_compliance_fraction,
            self.loss_factor,
            self.total_mass_kg,
        ] {
            h.update(&v.to_bits().to_le_bytes());
        }
        h.update(&(self.modes.len() as u64).to_le_bytes());
        for m in &self.modes {
            for v in [
                m.frequency_hz,
                m.mass_kg,
                m.stiffness_n_m,
                m.damping_n_s_m,
                m.participation_sqrt_kg,
            ] {
                h.update(&v.to_bits().to_le_bytes());
            }
        }
        h.finalize()
    }

    /// Canonical line-oriented bytes (schema line first, fixed field
    /// order, `{:e}` floats, identity LAST so decode can verify).
    #[must_use]
    pub fn to_canonical_bytes(&self) -> Vec<u8> {
        let mut s = String::new();
        s.push_str("frankensim-valve-card-v2\n");
        let _ = writeln!(s, "source_id\t{}", self.source_id);
        let _ = writeln!(s, "material_identity\t{}", self.material_identity.to_hex());
        let _ = writeln!(s, "mesh_digest\t{}", self.mesh_digest.to_hex());
        match self.discretization {
            ValveDiscretization::Tet4 => {
                let _ = writeln!(s, "element_order\ttet4");
            }
            ValveDiscretization::Tet10(q) => {
                let _ = writeln!(
                    s,
                    "element_order\ttet10\t{}\t{}",
                    q.stiffness_points_per_dir, q.mass_points_per_dir
                );
            }
        }
        for (key, v) in [
            ("rest_gap_m", self.rest_gap_m),
            ("width_m", self.width_m),
            ("effective_thickness_m", self.effective_thickness_m),
            (
                "retained_compliance_fraction",
                self.retained_compliance_fraction,
            ),
            ("loss_factor", self.loss_factor),
            ("total_mass_kg", self.total_mass_kg),
        ] {
            let _ = writeln!(s, "{key}\t{v:e}");
        }
        let _ = writeln!(s, "modes\t{}", self.modes.len());
        for m in &self.modes {
            let _ = writeln!(
                s,
                "mode\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}",
                m.frequency_hz,
                m.mass_kg,
                m.stiffness_n_m,
                m.damping_n_s_m,
                m.participation_sqrt_kg
            );
        }
        let _ = writeln!(s, "identity\t{}", self.identity.to_hex());
        s.into_bytes()
    }

    /// Decode + VERIFY canonical bytes (identity recomputed and matched;
    /// tampered bytes refuse).
    ///
    /// # Errors
    /// [`ReduceError::Card`] naming the violated rule.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<ValveCard, ReduceError> {
        fn take<'a>(
            lines: &[&'a str],
            cursor: &mut usize,
            key: &str,
        ) -> Result<&'a str, ReduceError> {
            let line = lines.get(*cursor).copied().ok_or(ReduceError::Card {
                what: "truncated card",
            })?;
            *cursor += 1;
            let (k, v) = line.split_once('\t').ok_or(ReduceError::Card {
                what: "malformed card line",
            })?;
            if k != key {
                return Err(ReduceError::Card {
                    what: "card field order violated",
                });
            }
            Ok(v)
        }
        fn take_num(lines: &[&str], cursor: &mut usize, key: &str) -> Result<f64, ReduceError> {
            take(lines, cursor, key)?
                .parse::<f64>()
                .map_err(|_| ReduceError::Card {
                    what: "bad numeric field",
                })
        }
        let text = core::str::from_utf8(bytes).map_err(|_| ReduceError::Card {
            what: "card bytes are not UTF-8",
        })?;
        let lines: Vec<&str> = text.lines().collect();
        let mut cursor = 0usize;
        if lines.first().copied() != Some("frankensim-valve-card-v2") {
            return Err(ReduceError::Card {
                what: "schema line mismatch",
            });
        }
        cursor += 1;
        let source_id = take(&lines, &mut cursor, "source_id")?.to_string();
        let material_identity =
            ContentHash::from_hex(take(&lines, &mut cursor, "material_identity")?).ok_or(
                ReduceError::Card {
                    what: "bad material identity hex",
                },
            )?;
        let mesh_digest = ContentHash::from_hex(take(&lines, &mut cursor, "mesh_digest")?).ok_or(
            ReduceError::Card {
                what: "bad mesh digest hex",
            },
        )?;
        let order = take(&lines, &mut cursor, "element_order")?;
        let discretization = if order == "tet4" {
            ValveDiscretization::Tet4
        } else {
            let fields: Vec<_> = order.split('\t').collect();
            if fields.len() != 3 || fields[0] != "tet10" {
                return Err(ReduceError::Card {
                    what: "unknown element order or missing quadrature",
                });
            }
            let parse = |s: &str| {
                s.parse::<usize>().map_err(|_| ReduceError::Card {
                    what: "bad quadrature count",
                })
            };
            let q = TetQuadrature {
                stiffness_points_per_dir: parse(fields[1])?,
                mass_points_per_dir: parse(fields[2])?,
            };
            q.validate().map_err(|_| ReduceError::Card {
                what: "inadmissible quadrature",
            })?;
            ValveDiscretization::Tet10(q)
        };
        let rest_gap_m = take_num(&lines, &mut cursor, "rest_gap_m")?;
        let width_m = take_num(&lines, &mut cursor, "width_m")?;
        let effective_thickness_m = take_num(&lines, &mut cursor, "effective_thickness_m")?;
        let retained_compliance_fraction =
            take_num(&lines, &mut cursor, "retained_compliance_fraction")?;
        let loss_factor = take_num(&lines, &mut cursor, "loss_factor")?;
        let total_mass_kg = take_num(&lines, &mut cursor, "total_mass_kg")?;
        let n_modes = take(&lines, &mut cursor, "modes")?
            .parse::<usize>()
            .map_err(|_| ReduceError::Card {
                what: "bad mode count",
            })?;
        // Every mode needs an input row, followed by one identity row.
        // Do not allocate from an unchecked file count. A zero-mode card
        // also cannot supply the first mode required by beating_reed().
        if n_modes == 0 || n_modes > lines.len().saturating_sub(cursor).saturating_sub(1) {
            return Err(ReduceError::Card {
                what: "mode count must be nonzero and fit the remaining input rows",
            });
        }
        let mut modes = Vec::with_capacity(n_modes);
        for _ in 0..n_modes {
            let line = lines.get(cursor).copied().ok_or(ReduceError::Card {
                what: "truncated mode row",
            })?;
            cursor += 1;
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() != 6 || cols[0] != "mode" {
                return Err(ReduceError::Card {
                    what: "malformed mode row",
                });
            }
            let mut vals = [0.0f64; 5];
            for (slot, col) in vals.iter_mut().zip(&cols[1..]) {
                *slot = col.parse::<f64>().map_err(|_| ReduceError::Card {
                    what: "bad mode value",
                })?;
            }
            modes.push(ReducedMode {
                frequency_hz: vals[0],
                mass_kg: vals[1],
                stiffness_n_m: vals[2],
                damping_n_s_m: vals[3],
                participation_sqrt_kg: vals[4],
            });
        }
        let identity = ContentHash::from_hex(take(&lines, &mut cursor, "identity")?).ok_or(
            ReduceError::Card {
                what: "bad identity hex",
            },
        )?;
        if cursor != lines.len() {
            return Err(ReduceError::Card {
                what: "trailing records after card identity",
            });
        }
        let card = ValveCard {
            identity,
            source_id,
            material_identity,
            mesh_digest,
            discretization,
            modes,
            rest_gap_m,
            width_m,
            effective_thickness_m,
            retained_compliance_fraction,
            loss_factor,
            total_mass_kg,
        };
        if card.recomputed_identity() != card.identity {
            return Err(ReduceError::Card {
                what: "identity mismatch: card bytes were altered",
            });
        }
        Ok(card)
    }

    /// The card's 1-DOF island as an `fs_scenario::BeatingReed` (massive
    /// branch: mode 0's mass and stiffness). The caller supplies the effective
    /// pressure area [m²], conjugate to this mode's mean-opening coordinate;
    /// `P_c = k H / A` makes the pressure force and swept flow use that area.
    /// The reduction does not yet derive this distributed-load projection.
    /// Supplying an authored area does not upgrade the card's authority.
    ///
    /// # Errors
    /// Refuses empty cards, nonphysical converted parameters and unrepresentable
    /// closing pressure. No default face length or pressure area is substituted.
    pub fn beating_reed(
        &self,
        effective_pressure_area_m2: f64,
        blowing_pressure_pa: f64,
        attack_s: f64,
    ) -> Result<fs_scenario::BeatingReed, ReduceError> {
        let m0 = self.modes.first().ok_or(ReduceError::Card {
            what: "reed conversion requires a retained mode",
        })?;
        if ![
            effective_pressure_area_m2,
            self.rest_gap_m,
            self.width_m,
            m0.mass_kg,
            m0.stiffness_n_m,
        ]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0)
            || ![self.loss_factor, m0.damping_n_s_m, blowing_pressure_pa, attack_s]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0)
        {
            return Err(ReduceError::Card {
                what: "reed conversion requires finite physical area, mechanics and excitation",
            });
        }
        let closing_pressure_pa = m0.stiffness_n_m * self.rest_gap_m / effective_pressure_area_m2;
        // Retain the legacy bit pattern only when its authored coefficient
        // reproduces this mode's ACTUAL damping exactly. Source-derived or
        // independently edited modal damping always determines the ratio.
        let damping_ratio = if m0.damping_n_s_m
            == self.loss_factor * (m0.stiffness_n_m * m0.mass_kg).sqrt()
        {
            0.5 * self.loss_factor
        } else {
            0.5 * (m0.damping_n_s_m / m0.stiffness_n_m.sqrt()) / m0.mass_kg.sqrt()
        };
        if !closing_pressure_pa.is_finite() || closing_pressure_pa <= 0.0
            || !damping_ratio.is_finite()
        {
            return Err(ReduceError::Card {
                what: "reed closing pressure is not positive finite",
            });
        }
        Ok(fs_scenario::BeatingReed {
            rest_opening_m: self.rest_gap_m,
            width_m: self.width_m,
            closing_pressure_pa,
            blowing_pressure_pa,
            attack_s,
            mass_kg: m0.mass_kg,
            stiffness_n_m: m0.stiffness_n_m,
            damping_ratio,
        })
    }

    /// JSON-lines log rows: the modal ladder with participations and the
    /// retained-compliance disclosure. The legacy `loss_factor` field is the
    /// authored contribution only; each mode's `damping_n_s_m` is authoritative.
    #[must_use]
    pub fn debug_lines(&self) -> Vec<String> {
        let mut out = vec![format!(
            "{{\"suite\":\"fs-solid\",\"case\":\"valve-card\",\"source_id\":\"{}\",\
             \"rest_gap_m\":{:.4e},\"width_m\":{:.4e},\"thickness_m\":{:.4e},\
             \"retained_compliance_fraction\":{:.4},\"loss_factor\":{:.3},\
             \"total_mass_kg\":{:.4e}}}",
            self.source_id,
            self.rest_gap_m,
            self.width_m,
            self.effective_thickness_m,
            self.retained_compliance_fraction,
            self.loss_factor,
            self.total_mass_kg
        )];
        for (i, m) in self.modes.iter().enumerate() {
            out.push(format!(
                "{{\"suite\":\"fs-solid\",\"case\":\"valve-card-mode\",\"index\":{i},\
                 \"f_hz\":{:.3},\"mass_kg\":{:.4e},\"stiffness_n_m\":{:.4e},\
                 \"damping_n_s_m\":{:.4e},\"participation_sqrt_kg\":{:.4e}}}",
                m.frequency_hz,
                m.mass_kg,
                m.stiffness_n_m,
                m.damping_n_s_m,
                m.participation_sqrt_kg
            ));
        }
        out
    }
}

#[cfg(test)]
#[path = "reduce/tests.rs"]
mod reduce_tests;
