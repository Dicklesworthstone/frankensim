//! Stationary reference-mass enthalpy transport on the existing P1 mesh.
//!
//! With specific enthalpy `h` [J/kg], frozen reference density `rho0` and
//! `q = -k(T,f) grad(T)`, the weak balance is
//! `integral(rho0 N_i dh/dt) + integral(grad(N_i) k grad(T)) = load_i`.
//! Neumann flux is positive OUTWARD; Robin contributes `htc*(T-T_ref)`
//! outward. Existing source, boundary and contact assemblers own all spatial
//! quadrature. Energetic internal variables and geometry are frozen here.
//!
//! Nodal constitutive evaluation gives `T_i=T(h_i)` and the existing P1
//! temperature interpolant. Row-sum storage gives invariant reference masses
//! `m_i=sum_e integral(rho0 N_i)`, equal to `sum_e rho0 V_e/4` for
//! cellwise density. Positive P1 nodal density retains its exact row moments.
//! Backward Euler solves, in joules,
//! `R(h)=m*(h-h_old)+dt*(A(T(h),f(h))*T(h)-b)=0`. The exact spatial
//! Jacobian includes both temperature-column scaling and the direct phase
//! conductivity derivative. An optional caller-declared law multiplies each
//! element's existing tensor by `s(mean_nodal_liquid_mass_fraction)`. Its
//! derivative can remain nonzero on a latent plateau where `T'=0`; reference
//! mass still supplies storage there. No apparent heat capacity or extra
//! latent account is introduced. The shared conduction assembly retains the
//! exact base `k'(T)` contribution, natural boundaries and contacts.
//!
//! The source-bound melting adapter instead supplies phase-specific nodal
//! conductivity and its enthalpy derivative. P1 interpolation gives the exact
//! element mean of that admitted conductivity field, through this same spatial
//! assembly. Each phase table stays within its original temperature support.
//!
//! Uniform and explicitly assigned heterogeneous reference materials share this
//! transport kernel, natural flux/Robin boundaries and finite-resistance
//! contacts. Distinct materials require distinct interface vertices; no nodal
//! mixture law is inferred. Prescribed endpoint temperatures are inverted only
//! where the original chart gives one unique enthalpy. A latent-plateau
//! temperature remains ambiguous and is refused. Reservoir reaction includes
//! fixed-node storage and transport, including an initial temperature jump.
//! Equilibrium density from the chart never changes the stored reference mass.
//! This does not model expansion, material motion, pressure work or remapping.

/// Prescribed spectral heating on the original spatial enthalpy balance.
pub mod interference;

/// Implicit endpoint derivatives through the accepted physical balance.
pub mod adjoint;
/// Error-controlled time integration of the same nonlinear enthalpy transport.
pub mod adaptive;
/// Second-order SDIRK stages on the same nonlinear enthalpy residual.
pub mod sdirk;
/// Implicit closed-enclosure radiation with complete physical enthalpy feedback.
pub mod enclosure;
/// Explicit heterogeneous phase charts and frozen reference densities.
pub mod heterogeneous;
/// Paired source caloric storage and explicitly isotropic conductivity.
pub mod material;
/// Growing-cell deposition adapters for existing enthalpy charts and Cx.
pub mod deposition;
/// Implicit ambient radiation with immutable enthalpy history and joule gates.
pub mod radiation;
pub use radiation::adjoint::{EnthalpyRadiationStepGradient, EnthalpyRadiationStepLinearization};
mod endpoint;
pub use endpoint::EnthalpyEndpointEvaluation;
mod contact;
mod pointwise;
mod reference_mass;
use pointwise::{PointwiseBoundary, PointwiseBoundaryAssembly};

use std::{cell::RefCell, fmt};

use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NewtonStallDiagnosis,
    NonlinearProblem, SolverCallback, SolverRunProgress,
};
use fs_sparse::Csr;
use crate::radiation::{
    PointwiseAmbientRadiation, PointwiseRadiationReport,
    PointwiseEnclosureRadiation, PointwiseEnclosureReport,
};

use crate::{
    ConductionError, ConductionMesh, ConductionProblem, InterfaceFlux, RobinFlux,
    ThermalInterfaces,
    assemble::{
        ASSEMBLY_TILE, AssembledSystem, DofMap,
        assemble_enthalpy_jacobian_with_optional_interfaces,
        assemble_operator_scaled_with_interfaces,
    },
    material::LiquidMassFractionConductivity,
    solve::energy_balance,
};

/// Spatial work and allocation admission for one prepared reference mass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnthalpyBudget {
    /// Maximum nodal enthalpy unknowns, including boundary nodes.
    pub max_vertices: usize,
    /// Maximum tetrahedra used by storage and transport assembly.
    pub max_elements: usize,
}

/// Explicit nonlinear/Krylov work and absolute energy-closure policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnthalpyStepConfig {
    /// Shared Newton controls; absolute residual tolerance is in joules.
    /// Restart and cycle bounds apply to each bounded Newton attempt.
    pub newton: NewtonKrylovConfig,
    /// Maximum attempts for this endpoint; no unconverged endpoint is returned.
    pub max_newton_iterations: usize,
    /// Maximum absolute storage-minus-external-input defect, in joules.
    pub energy_tolerance_j: f64,
}

/// A refused construction or step returns no new physical state.
#[derive(Debug, Clone, PartialEq)]
pub enum EnthalpyError {
    /// Invalid shape, nonpositive physical input, or unrepresentable size/work.
    InvalidInput(&'static str),
    /// The caller's spatial admission cap is insufficient.
    Budget {
        /// Which admitted resource exceeded its cap.
        resource: &'static str,
        /// Actual mesh requirement.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
    /// The requesting consumer does not support prescribed-temperature rows.
    UnsupportedDirichlet,
    /// A fixed temperature selects more than one enthalpy on a latent plateau.
    AmbiguousPrescribedTemperature {
        /// Boundary vertex whose chart contains the plateau.
        vertex: usize,
        /// Requested absolute temperature, K.
        temperature_k: f64,
    },
    /// A prescribed temperature is outside the original constitutive chart.
    PrescribedTemperatureOutsideDomain {
        /// Boundary vertex with insufficient constitutive support.
        vertex: usize,
        /// Requested absolute temperature, K.
        temperature_k: f64,
        /// Lowest chart temperature, K.
        lower_k: f64,
        /// Highest chart temperature, K.
        upper_k: f64,
    },
    /// Existing spatial/material/interface refusal, including cancellation.
    Conduction(ConductionError),
    /// An exact nodal chart or derivative query was refused.
    Phase {
        /// Nodal state that could not be resolved.
        vertex: usize,
        /// Original constitutive refusal.
        source: PhaseStateError,
    },
    /// A finite nodal enthalpy leaves the original conductivity source support.
    /// Newton trial proposals may backtrack; initial/direct states are refused.
    SourceEnthalpyOutsideDomain {
        /// Vertex whose source cannot resolve this enthalpy.
        vertex: usize,
        /// Unsupported specific enthalpy, J/kg.
        specific_enthalpy_j_kg: f64,
        /// Lower original source bound, J/kg.
        lower_j_kg: f64,
        /// Upper original source bound, J/kg.
        upper_j_kg: f64,
    },
    /// Shared Newton configuration or finite-residual admission failed.
    Newton(NewtonError),
    /// Shared Newton/Krylov work or globalization failed to converge.
    NotConverged(NewtonReport),
    /// Independent physical nodal residual exceeds its unchanged local gate.
    /// This is not a Newton stall or a whole-domain energy defect.
    EndpointResidual {
        residual_j: f64,
        tolerance_j: f64,
    },
    /// A nonlinear callback unwound; the incomplete Newton attempt was discarded.
    SolverCallbackPanicked(SolverCallback),
    /// Independently integrated input did not balance the published storage.
    EnergyBalance {
        /// Actual storage minus net external input, joules.
        residual_j: f64,
        /// Declared absolute allowance, joules.
        tolerance_j: f64,
    },
    /// Finite physical inputs produced an unrepresentable numerical quantity.
    NonFiniteArithmetic,
}

impl fmt::Display for EnthalpyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(reason) => write!(formatter, "enthalpy input: {reason}"),
            Self::Budget {
                resource,
                required,
                limit,
            } => {
                write!(
                    formatter,
                    "enthalpy {resource} requires {required}, cap {limit}"
                )
            }
            Self::UnsupportedDirichlet => write!(
                formatter,
                "this enthalpy consumer does not support prescribed-temperature boundaries"
            ),
            Self::AmbiguousPrescribedTemperature { vertex, temperature_k } => write!(
                formatter,
                "enthalpy vertex {vertex}: prescribed {temperature_k} K selects a latent plateau, not a unique enthalpy"
            ),
            Self::PrescribedTemperatureOutsideDomain { vertex, temperature_k, lower_k, upper_k } => write!(
                formatter,
                "enthalpy vertex {vertex}: prescribed {temperature_k} K is outside [{lower_k}, {upper_k}] K"
            ),
            Self::Conduction(error) => error.fmt(formatter),
            Self::Phase { vertex, source } => {
                write!(formatter, "enthalpy vertex {vertex}: {source}")
            }
            Self::SourceEnthalpyOutsideDomain {
                vertex, specific_enthalpy_j_kg, lower_j_kg, upper_j_kg,
            } => write!(
                formatter,
                "enthalpy vertex {vertex}: {specific_enthalpy_j_kg} J/kg is outside the conductivity source [{lower_j_kg}, {upper_j_kg}] J/kg"
            ),
            Self::Newton(error) => error.fmt(formatter),
            Self::NotConverged(report) => write!(
                formatter,
                "enthalpy Newton refused after {} attempts: {:?}; residual {} J",
                report.iterations, report.diagnosis, report.residual_norm
            ),
            Self::EndpointResidual { residual_j, tolerance_j } => write!(
                formatter,
                "physical enthalpy endpoint residual {residual_j} J exceeds {tolerance_j} J"
            ),
            Self::SolverCallbackPanicked(callback) => {
                write!(
                    formatter,
                    "enthalpy Newton {callback:?} callback panicked; step discarded"
                )
            }
            Self::EnergyBalance {
                residual_j,
                tolerance_j,
            } => write!(
                formatter,
                "enthalpy energy defect {residual_j} J exceeds {tolerance_j} J"
            ),
            Self::NonFiniteArithmetic => {
                write!(formatter, "enthalpy arithmetic is not representable")
            }
        }
    }
}

impl std::error::Error for EnthalpyError {}

impl From<ConductionError> for EnthalpyError {
    fn from(error: ConductionError) -> Self {
        Self::Conduction(error)
    }
}

/// Accepted nodal enthalpy, temperature and phase with a direct energy audit.
#[derive(Debug, Clone)]
pub struct EnthalpyStepSolution {
    /// Published specific enthalpy at each vertex, J/kg; use as the next history.
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// Chart-resolved nodal absolute temperature, K.
    pub temperature: Vec<f64>,
    /// Chart-resolved liquid mass fractions in [0,1].
    pub liquid_mass_fraction: Vec<f64>,
    /// Endpoint volumetric generation integrated by the existing source rule, W.
    pub source_w: f64,
    /// Endpoint outward Neumann heat, W; negative means imposed heating.
    pub neumann_out_w: f64,
    /// Endpoint outward Robin heat, W.
    pub robin_out_w: f64,
    /// Spatially integrated nonlinear radiation, when explicitly requested.
    /// Its heat is already included in the energy balance and is separate from
    /// `robin_out_w`; it must not be handed to a convection/coolant heat port.
    pub pointwise_radiation: Option<PointwiseRadiationReport>,
    /// Spatial local emission with reflected, patch-uniform irradiation.
    /// Net physical-surface heat is included once in the energy balance, apart
    /// from authored Neumann/convection. An open model separately records its
    /// black bath; that bath is never a second debit against solid storage.
    pub pointwise_enclosure_radiation: Option<PointwiseEnclosureReport>,
    /// Prescribed-temperature reservoir heat, W inward, including fixed-node
    /// storage changes. A first-step boundary jump is charged exactly once.
    pub dirichlet_in_w: f64,
    /// `sum_i m_i*(h_new_i-h_old_i)` of the actual returned state, J.
    pub stored_energy_change_j: f64,
    /// Storage minus dt times net external input, J.
    pub energy_residual_j: f64,
    /// Existing per-region endpoint convection reports.
    pub robin_fluxes: Vec<RobinFlux>,
    /// Internal contact exchanges, excluded from net external input.
    pub contact_fluxes: Vec<InterfaceFlux>,
    /// Shared Newton's true residual and bounded inner-solve telemetry.
    pub newton: NewtonReport,
}

/// Fixed reference nodal masses and one immutable equilibrium enthalpy chart.
#[derive(Debug)]
pub struct EnthalpyBackwardEuler<'m, 'c> {
    mesh: &'m ConductionMesh,
    curve: &'c EquilibriumEnthalpyPhaseCurve,
    // Only the heterogeneous wrapper constructs this private override. The
    // public uniform API, including phase_curve(), keeps its original meaning.
    nodal_curves: Option<Vec<&'c EquilibriumEnthalpyPhaseCurve>>,
    masses: Vec<f64>,
    phase_conductivity: Option<PhaseConductivity>,
    // A source-bound melting adapter supplies nodal k(h), using the existing
    // element-scale assembly with its retained unit tensor as numerical basis.
    sourced_melting: Option<material::melting::MeltingFieldTransport<'c>>,
    // Original multiaxis source integrated over P1 h, with no inferred phase mix.
    sourced_phase: Option<material::phase::PhaseFieldTransport<'c>>,
}

#[derive(Debug)]
struct PhaseConductivity {
    element_laws: Vec<LiquidMassFractionConductivity>,
    // A fraction-slope corner matters only where an incident element actually
    // changes conductivity with that fraction. Constant multipliers do not.
    variable_vertices: Vec<bool>,
}

fn reference_mass_node_count(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    budget: EnthalpyBudget,
) -> Result<usize, EnthalpyError> {
    poll(cx, 0)?;
    let n = mesh.vertex_count();
    for (resource, required, limit) in [
        ("vertices", n, budget.max_vertices),
        ("elements", mesh.element_count(), budget.max_elements),
    ] {
        if required > limit {
            return Err(EnthalpyError::Budget {
                resource,
                required,
                limit,
            });
        }
    }
    if n == 0 || mesh.element_count() == 0 {
        return Err(EnthalpyError::InvalidInput(
            "nonempty mesh and positive finite reference density required",
        ));
    }
    // Mesh/contact assembly is at most dense in its admitted vertex set.
    // Check both dense-index and element-staging products before allocation.
    n.checked_mul(n)
        .and_then(|v| v.checked_mul(24))
        .and_then(|_| mesh.element_count().checked_mul(32 * 24))
        .ok_or(EnthalpyError::InvalidInput(
            "spatial allocation size overflow",
        ))?;
    Ok(n)
}

impl<'m, 'c> EnthalpyBackwardEuler<'m, 'c> {
    /// Prepare invariant P1 nodal masses from declared reference density [kg/m3].
    /// This is mass quadrature, not a heat-capacity declaration. The curve's
    /// state-dependent density is never substituted into these masses.
    pub fn uniform(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        curve: &'c EquilibriumEnthalpyPhaseCurve,
        reference_density_kg_m3: f64,
        budget: EnthalpyBudget,
    ) -> Result<Self, EnthalpyError> {
        let n = reference_mass_node_count(cx, mesh, budget)?;
        if !reference_density_kg_m3.is_finite() || reference_density_kg_m3 <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "nonempty mesh and positive finite reference density required",
            ));
        }
        let mut masses = vec![0.0; n];
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let contribution =
                finite(reference_density_kg_m3 * mesh.element_volume(element) / 4.0)?;
            if contribution <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "positive reference mass is not representable",
                ));
            }
            for &vertex in vertices {
                masses[vertex as usize] = finite(masses[vertex as usize] + contribution)?;
            }
        }
        Self::with_reference_masses(cx, mesh, curve, masses)
    }

    /// Prepare the native enthalpy balance with a positive P1 reference density
    /// supplied once per original mesh vertex, in kg/m³.
    ///
    /// Each tetrahedron contributes `V*(sum(rho)+rho_a)/20` to vertex a.
    /// These exact row moments preserve total mass and the integral of a P1
    /// specific-enthalpy field. A cell mean loses the correlation between
    /// density and enthalpy, and produces incorrect local heating/phase rates.
    /// Globally uniform input retains the original uniform-constructor bits.
    ///
    /// The resulting reference masses remain fixed during every nonlinear
    /// trial and time step, including latent plateaus and prescribed-temperature
    /// reservoir accounting. Existing history/source adjoints use these same
    /// masses; this entry does not add density-control derivatives or reinterpret
    /// heterogeneous material-density controls. The caloric chart's equilibrium
    /// density never replaces this declared inventory. Geometry remains fixed.
    ///
    /// # Errors
    /// Wrong nodal count, nonpositive/nonfinite density, original mesh/work
    /// admission, unrepresentable mass or cancellation.
    pub fn from_nodal_reference_density(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        curve: &'c EquilibriumEnthalpyPhaseCurve,
        reference_nodal_density_kg_m3: &[f64],
        budget: EnthalpyBudget,
    ) -> Result<Self, EnthalpyError> {
        let n = reference_mass_node_count(cx, mesh, budget)?;
        if reference_nodal_density_kg_m3.len() != n
            || reference_nodal_density_kg_m3
                .iter()
                .any(|rho| !rho.is_finite() || *rho <= 0.0)
        {
            return Err(EnthalpyError::InvalidInput(
                "one positive finite reference density per mesh vertex required",
            ));
        }
        if reference_nodal_density_kg_m3
            .iter()
            .all(|rho| rho.to_bits() == reference_nodal_density_kg_m3[0].to_bits())
        {
            return Self::uniform(cx, mesh, curve, reference_nodal_density_kg_m3[0], budget);
        }
        let mut masses = vec![0.0; n];
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let density = vertices.map(|vertex| reference_nodal_density_kg_m3[vertex as usize]);
            let local = reference_mass::p1_nodal_masses(mesh.element_volume(element), density)
                .ok_or(EnthalpyError::InvalidInput(
                    "positive graded reference mass is not representable",
                ))?;
            for (a, &vertex) in vertices.iter().enumerate() {
                masses[vertex as usize] = finite(masses[vertex as usize] + local[a])?;
            }
        }
        Self::with_reference_masses(cx, mesh, curve, masses)
    }

    fn with_reference_masses(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        curve: &'c EquilibriumEnthalpyPhaseCurve,
        masses: Vec<f64>,
    ) -> Result<Self, EnthalpyError> {
        for (vertex, &mass) in masses.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if mass <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "every vertex requires positive reference mass",
                ));
            }
        }
        poll(cx, masses.len())?;
        Ok(Self {
            mesh,
            curve,
            nodal_curves: None,
            masses,
            phase_conductivity: None,
            sourced_melting: None,
            sourced_phase: None,
        })
    }

    /// Bind an explicit phase multiplier to every tetrahedron, replacing any
    /// previous manual assignment on this consumed owner. A source-bound
    /// melting law already owns complete conductivity and refuses this override.
    ///
    /// Each law is evaluated at the arithmetic mean of the element's nodal
    /// liquid **mass** fractions. It scales the existing, independently
    /// assigned conductivity tensor at mean nodal temperature. Boundary and
    /// contact conductances are not scaled. Reference mass and geometry remain
    /// fixed. The same law participates in residuals, Newton and adjoints.
    pub fn with_phase_conductivity(
        mut self,
        cx: &Cx<'_>,
        element_laws: Vec<LiquidMassFractionConductivity>,
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        if self.sourced_melting.is_some() || self.sourced_phase.is_some() {
            return Err(EnthalpyError::InvalidInput(
                "source-bound phase transport already owns the complete conductivity law",
            ));
        }
        if element_laws.len() != self.mesh.element_count() {
            return Err(EnthalpyError::InvalidInput(
                "one phase conductivity law per tetrahedron required",
            ));
        }
        let mut variable_vertices = vec![false; self.masses.len()];
        for (element, (vertices, law)) in self
            .mesh
            .complex()
            .tets
            .iter()
            .zip(&element_laws)
            .enumerate()
        {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            if law.derivative() != 0.0 {
                for &vertex in vertices {
                    variable_vertices[vertex as usize] = true;
                }
            }
        }
        poll(cx, element_laws.len())?;
        self.phase_conductivity = Some(PhaseConductivity {
            element_laws,
            variable_vertices,
        });
        Ok(self)
    }

    /// Invariant reference masses in nodal mesh order, kg.
    #[must_use]
    pub fn reference_nodal_masses_kg(&self) -> &[f64] {
        &self.masses
    }

    /// Immutable constitutive chart retained by this spatial storage model.
    #[must_use]
    pub const fn phase_curve(&self) -> &EquilibriumEnthalpyPhaseCurve {
        self.curve
    }

    // Callers have already admitted the nodal shape. Keep constitutive lookup
    // common to residuals, Newton tangents, phase reporting and adjoint guards.
    fn curve_for_vertex(&self, vertex: usize) -> &EquilibriumEnthalpyPhaseCurve {
        self.nodal_curves
            .as_ref()
            .map_or(self.curve, |curves| curves[vertex])
    }

    /// Advance from immutable specific-enthalpy history, with current endpoint
    /// source and boundary data. Prescribed temperatures may change from the
    /// history, but must select unique enthalpies on their original charts.
    /// Conductivity is evaluated at every
    /// actual nonlinear trial. No returned field escapes before both Newton's
    /// residual gate and the independent absolute joule-balance gate pass.
    ///
    /// The mesh caps bound spatial buffers. Before Newton allocation, checked
    /// products bound the restart basis/Hessenberg sizes and the total maximum
    /// Krylov columns `max_newton_iterations*restart*cycles`. Each shared Newton
    /// attempt uses its declared globalization bounds; Cx is polled between
    /// attempts and inside constitutive, assembly and sparse-action tiles.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.advance_radiative(cx, problem, interfaces, old_specific_h, dt_s, config, None)
    }

    /// Advance with the integral of local `epsilon sigma (T_h^4 - T_a^4)`
    /// on each prepared surface triangle, including its exact Newton tangent.
    ///
    /// Conduction, physical source, contact, convection, fixed-node reaction
    /// and latent storage are solved together in the original enthalpy system.
    /// Radiation has its own power account; hot spots are never replaced by a
    /// patch mean. The prepared trace must borrow this exact mesh and boundary.
    /// The existing Newton work and joule gates apply; no outer radiation fixed
    /// point, finite-difference Jacobian or artificial heat capacity is used.
    /// This is first-order time integration with fixed emissivity and geometry,
    /// not an emissivity-law derivative, enclosure or spatial-error certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &PointwiseAmbientRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        poll(cx, 0)?;
        if !radiation.matches(problem.mesh, problem.boundary) {
            return Err(EnthalpyError::InvalidInput(
                "pointwise radiation must retain the exact transport mesh and boundary",
            ));
        }
        self.advance_radiative(
            cx, problem, interfaces, old_specific_h, dt_s, config,
            Some(PointwiseBoundary::Ambient(radiation)),
        )
    }

    /// Advance with original-material nodal emissivity on the same P1 trace.
    /// The full h tangent includes phase-dependent emissivity even when dT/dh
    /// vanishes. Original Robin heat stays separate from radiative heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_material_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_specific_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseAmbientRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        poll(cx, 0)?;
        if !std::ptr::eq(self.mesh, problem.mesh)
            || !radiation.matches(problem.mesh, problem.boundary)
            || !radiation.matches_curves(cx, |vertex| self.curve_for_vertex(vertex))? {
            return Err(EnthalpyError::InvalidInput(
                "material radiation must retain the exact mesh, boundary and original nodal caloric objects",
            ));
        }
        self.advance_radiative(cx, problem, interfaces, old_specific_h, dt_s, config,
            Some(PointwiseBoundary::MaterialAmbient(radiation)))
    }

    /// Advance original-material local emission and freshly updated reflected
    /// irradiation in the same h solve, including isothermal phase response.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_material_enclosure_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], dt_s: f64,
        config: EnthalpyStepConfig, radiation: &crate::MaterialPointwiseEnclosureRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        poll(cx, 0)?;
        if !std::ptr::eq(self.mesh, problem.mesh)
            || !radiation.matches(problem.mesh, problem.boundary)
            || !radiation.matches_curves(cx, |vertex| self.curve_for_vertex(vertex))? {
            return Err(EnthalpyError::InvalidInput(
                "material enclosure requires the original mesh, boundary and nodal caloric objects",
            ));
        }
        self.advance_radiative(cx, problem, interfaces, old_h, dt_s, config,
            Some(PointwiseBoundary::MaterialEnclosure(radiation)))
    }

    /// Advance local P1 emission and reflected enclosure exchange in the SAME
    /// nonlinear enthalpy solve. The complete generally nonsymmetric radiation
    /// tangent uses the existing Newton-FGMRES lane. No outer frozen-radiation
    /// iteration or patch-mean temperature replaces a local hot spot.
    ///
    /// Original convection, prescribed flux, contact, nonlinear conductivity,
    /// latent storage and prescribed-temperature reactions remain active.
    /// An owner from `PointwiseEnclosureRadiation::new_open` also exchanges
    /// heat with its explicit black bath. Finite storage supplies the transient
    /// anchor without an artificial Dirichlet/Robin boundary or heat capacity.
    /// Irradiation is uniform within each declared view-factor patch; this
    /// does not infer spatial visibility, emissivity variation or moving geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &PointwiseEnclosureRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        poll(cx, 0)?;
        if !radiation.matches(problem.mesh, problem.boundary) {
            return Err(EnthalpyError::InvalidInput(
                "pointwise enclosure must retain the exact transport mesh and boundary",
            ));
        }
        self.advance_radiative(
            cx, problem, interfaces, old_specific_h, dt_s, config,
            Some(PointwiseBoundary::Enclosure(radiation)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_radiative(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        pointwise_radiation: Option<PointwiseBoundary<'_>>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.admit_step(cx, problem, old_specific_h, dt_s, config)?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation,
        };
        let (h, report) = context.solve_newton(config, None)?;
        context.finish(h, report, config.energy_tolerance_j)
    }

    fn admit_step(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<(), EnthalpyError> {
        poll(cx, 0)?;
        if !std::ptr::eq(self.mesh, problem.mesh) || old_specific_h.len() != self.masses.len() {
            return Err(EnthalpyError::InvalidInput(
                "reference mass and history must match the exact transport mesh",
            ));
        }
        if let Some(source) = self.sourced_melting.as_ref() {
            if !source.matches_numerical_basis(problem) {
                return Err(EnthalpyError::InvalidInput(
                    "source-bound melting requires its original problem tensor and no element override",
                ));
            }
        }
        if let Some(source) = self.sourced_phase.as_ref() {
            if !source.matches_numerical_basis(problem) {
                return Err(EnthalpyError::InvalidInput(
                    "phase source requires its original problem tensor and no element override",
                ));
            }
            source.admit_enthalpies(cx, old_specific_h)?;
        }
        validate_work(self.masses.len(), dt_s, config)?;
        for condition in problem.boundary.conditions() {
            condition.validate(self.masses.len())?;
        }
        for (index, &(vertex, temperature)) in problem.boundary.dirichlet().iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            if vertex >= self.masses.len() {
                return Err(EnthalpyError::InvalidInput(
                    "prescribed vertex is outside the reference mesh",
                ));
            }
            self.prescribed_enthalpy(vertex, temperature)?;
        }
        Ok(())
    }

    // Invert the admitted piecewise-linear chart, never a separate Cp model.
    // Monotone temperatures permit binary search; both ends of every plateau
    // are ambiguous too, regardless of the history's liquid fraction.
    fn prescribed_enthalpy(&self, vertex: usize, temperature_k: f64) -> Result<f64, EnthalpyError> {
        let knots = self.curve_for_vertex(vertex).knots();
        let lower_k = knots[0].temperature_k;
        let upper_k = knots[knots.len() - 1].temperature_k;
        if !temperature_k.is_finite() || temperature_k < lower_k || temperature_k > upper_k {
            return Err(EnthalpyError::PrescribedTemperatureOutsideDomain {
                vertex,
                temperature_k,
                lower_k,
                upper_k,
            });
        }
        let upper = knots.partition_point(|knot| knot.temperature_k < temperature_k);
        if knots[upper].temperature_k == temperature_k {
            if knots
                .get(upper + 1)
                .is_some_and(|next| next.temperature_k == temperature_k)
            {
                return Err(EnthalpyError::AmbiguousPrescribedTemperature {
                    vertex,
                    temperature_k,
                });
            }
            return Ok(knots[upper].specific_enthalpy_j_kg);
        }
        let a = knots[upper - 1];
        let b = knots[upper];
        let fraction = (temperature_k - a.temperature_k) / (b.temperature_k - a.temperature_k);
        finite(a.specific_enthalpy_j_kg * (1.0 - fraction) + b.specific_enthalpy_j_kg * fraction)
    }
}

struct StepContext<'a, 'm, 'c, 'cx> {
    storage: &'a EnthalpyBackwardEuler<'m, 'c>,
    cx: &'a Cx<'cx>,
    problem: ConductionProblem<'a>,
    interfaces: Option<&'a ThermalInterfaces>,
    old: &'a [f64],
    dt: f64,
    pointwise_radiation: Option<PointwiseBoundary<'a>>,
}

impl StepContext<'_, '_, '_, '_> {
    // Both source modes use the original Newton/Krylov driver and tangent.
    // Prescribed interval joules are independent of h and change only the RHS.
    fn solve_newton<'s>(
        &'s self,
        config: EnthalpyStepConfig,
        source_j: Option<&'s [fs_math::dd::Dd]>,
    ) -> Result<(Vec<f64>, NewtonReport), EnthalpyError> {
        self.solve_newton_from(config, source_j, None)
    }

    fn solve_newton_from<'s>(
        &'s self, config: EnthalpyStepConfig,
        source_j: Option<&'s [fs_math::dd::Dd]>, initial_guess: Option<&[f64]>,
    ) -> Result<(Vec<f64>, NewtonReport), EnthalpyError> {
        // Resolve history and the initial tangent before shared Newton admits it.
        // A chart/domain refusal here retains its original typed error.
        let initial_h = self.initial_values_from(initial_guess)?;
        let mut stage = self.stage(&initial_h)?;
        stage.source_j = source_j;
        let mut newton =
            match NewtonKrylovState::new(&stage, initial_h, config.newton) {
                Ok(state) => state,
                Err(error) => {
                    return Err(stage.take_failure().unwrap_or(EnthalpyError::Newton(error)));
                }
            };
        poll(self.cx, 0)?;
        for attempt in 0..config.max_newton_iterations {
            poll(self.cx, attempt)?;
            let progress = newton.run_cancellable(&stage, 1, self.cx);
            match progress.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => {
                    return Err(ConductionError::Cancelled {
                        stage: "spatial-enthalpy-solver",
                        at: attempt,
                    }
                    .into());
                }
                SolverRunProgress::CallbackPanicked(callback) => {
                    return Err(EnthalpyError::SolverCallbackPanicked(callback));
                }
                SolverRunProgress::DimensionMismatch { expected, actual } => {
                    return Err(EnthalpyError::Newton(NewtonError::Dimension {
                        problem: actual,
                        state: expected,
                    }));
                }
            }
            if let Some(error) = stage.fatal.take() {
                return Err(error);
            }
            poll(self.cx, attempt)?;
            let report = progress.report;
            if report.converged {
                return Ok((newton.x, report));
            }
            if report.diagnosis != Some(NewtonStallDiagnosis::BudgetExhausted)
                || attempt + 1 == config.max_newton_iterations
            {
                return Err(stage
                    .take_failure()
                    .unwrap_or(EnthalpyError::NotConverged(report)));
            }
            // One shared Newton attempt uses one exact tangent point. Rebuild
            // after acceptance rather than freezing k(T) across Newton updates.
            stage = self.stage(&newton.x)?;
            stage.source_j = source_j;
        }
        unreachable!("positive Newton attempt budget admitted above")
    }

    fn is_prescribed(&self, vertex: usize) -> bool {
        self.problem
            .boundary
            .dirichlet()
            .binary_search_by_key(&vertex, |&(v, _)| v)
            .is_ok()
    }

    fn initial_values(&self) -> Result<Vec<f64>, EnthalpyError> {
        self.initial_values_from(None)
    }

    fn initial_values_from(&self, guess: Option<&[f64]>) -> Result<Vec<f64>, EnthalpyError> {
        // Validate even replaced boundary history before accepting its storage.
        if !self.problem.boundary.dirichlet().is_empty() {
            self.temperatures(self.old)?;
        }
        let values = guess.unwrap_or(self.old);
        if values.len() != self.old.len() {
            return Err(EnthalpyError::InvalidInput("initial enthalpy guess has the wrong nodal shape"));
        }
        let mut initial = values.to_vec();
        for (index, &(vertex, temperature)) in self.problem.boundary.dirichlet().iter().enumerate()
        {
            if index % ASSEMBLY_TILE == 0 {
                poll(self.cx, index)?;
            }
            initial[vertex] = self.storage.prescribed_enthalpy(vertex, temperature)?;
        }
        Ok(initial)
    }

    fn check_prescribed_values(&self, h: &[f64]) -> Result<(), EnthalpyError> {
        for (index, &(vertex, temperature)) in self.problem.boundary.dirichlet().iter().enumerate()
        {
            if index % ASSEMBLY_TILE == 0 {
                poll(self.cx, index)?;
            }
            if h[vertex] != self.storage.prescribed_enthalpy(vertex, temperature)? {
                return Err(EnthalpyError::InvalidInput(
                    "candidate enthalpy violates prescribed endpoint temperature",
                ));
            }
        }
        Ok(())
    }

    fn prescribed_heat_input_w(
        &self,
        h: &[f64],
        temperature: &[f64],
        system: &AssembledSystem,
    ) -> Result<f64, EnthalpyError> {
        let radiation = if self.problem.boundary.dirichlet().is_empty() {
            None
        } else {
            self.pointwise_radiation.map(|radiation| radiation.load_and_report(self.cx, temperature, h)).transpose()?
        };
        let mut heat = 0.0;
        for (index, &(vertex, _)) in self.problem.boundary.dirichlet().iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(self.cx, index)?;
            }
            let (columns, entries) = system.operator.row(vertex);
            let transport: f64 = columns
                .iter()
                .zip(entries)
                .map(|(&column, &entry)| entry * temperature[column])
                .sum();
            let transport = match &radiation {
                Some((load, _)) => finite(transport + load[vertex])?,
                None => transport,
            };
            heat = finite(
                heat + finite(
                    self.storage.masses[vertex] * (h[vertex] - self.old[vertex]) / self.dt
                        + transport
                        - system.load[vertex],
                )?,
            )?;
        }
        Ok(heat)
    }

    fn physical_prescribed_heat_input_w(&self, h: &[f64]) -> Result<f64, EnthalpyError> {
        if self.problem.boundary.dirichlet().is_empty() {
            return Ok(0.0);
        }
        let temperature = self.temperatures(h)?;
        let scales = self.phase_scales(h)?;
        let system = self.assemble(&temperature, scales.as_deref())?;
        self.prescribed_heat_input_w(h, &temperature, &system)
    }

    fn temperatures(&self, h: &[f64]) -> Result<Vec<f64>, EnthalpyError> {
        let mut temperature = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            let state = self
                .storage
                .curve_for_vertex(vertex)
                .state_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            temperature.push(state.temperature_k());
        }
        Ok(temperature)
    }

    fn phase_scales(&self, h: &[f64]) -> Result<Option<Vec<f64>>, EnthalpyError> {
        if let Some(source) = self.storage.sourced_phase.as_ref() {
            return source.assemble(self.cx, self.storage.mesh, h, false)
                .map(|(values, _)| Some(values));
        }
        if let Some(source) = self.storage.sourced_melting.as_ref() {
            let mut nodal = Vec::with_capacity(h.len());
            for (vertex, &value) in h.iter().enumerate() {
                if vertex % ASSEMBLY_TILE == 0 {
                    poll(self.cx, vertex)?;
                }
                nodal.push(source.sample_for_vertex(value, vertex)?.conductivity_w_m_k);
            }
            let mut scales = Vec::with_capacity(self.storage.mesh.element_count());
            for (element, vertices) in self.storage.mesh.complex().tets.iter().enumerate() {
                if element % ASSEMBLY_TILE == 0 {
                    poll(self.cx, element)?;
                }
                // Exactly integrate the P1 interpolant of nodal conductivity;
                // do not query either phase at an averaged off-phase state.
                let mean = vertices.iter().map(|&v| nodal[v as usize] / 4.0).sum();
                scales.push(finite(mean)?);
            }
            return Ok(Some(scales));
        }
        let Some(policy) = &self.storage.phase_conductivity else {
            return Ok(None);
        };
        let mut scales = Vec::with_capacity(policy.element_laws.len());
        for (element, (vertices, law)) in self
            .storage
            .mesh
            .complex()
            .tets
            .iter()
            .zip(&policy.element_laws)
            .enumerate()
        {
            if element % ASSEMBLY_TILE == 0 {
                poll(self.cx, element)?;
            }
            let mut fraction = 0.0;
            for &vertex in vertices {
                let vertex = vertex as usize;
                fraction += self
                    .storage
                    .curve_for_vertex(vertex)
                    .state_at_specific_enthalpy(h[vertex])
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .liquid_mass_fraction()
                    / 4.0;
            }
            scales.push(law.multiplier_at(fraction)?);
        }
        Ok(Some(scales))
    }

    fn phase_scale_derivatives(&self, h: &[f64]) -> Result<Option<Vec<[f64; 4]>>, EnthalpyError> {
        if let Some(source) = self.storage.sourced_melting.as_ref() {
            let mut nodal = Vec::with_capacity(h.len());
            for (vertex, &value) in h.iter().enumerate() {
                if vertex % ASSEMBLY_TILE == 0 {
                    poll(self.cx, vertex)?;
                }
                let derivative = source.sample_for_vertex(value, vertex)?.derivative_per_j_kg;
                let scaled = finite(derivative / 4.0)?;
                if derivative != 0.0 && scaled == 0.0 {
                    return Err(EnthalpyError::InvalidInput(
                        "source conductivity tangent is not representable in P1 quadrature",
                    ));
                }
                nodal.push(scaled);
            }
            let mut derivatives = Vec::with_capacity(self.storage.mesh.element_count());
            for (element, vertices) in self.storage.mesh.complex().tets.iter().enumerate() {
                if element % ASSEMBLY_TILE == 0 {
                    poll(self.cx, element)?;
                }
                derivatives.push(vertices.map(|v| nodal[v as usize]));
            }
            return Ok(Some(derivatives));
        }
        let Some(policy) = &self.storage.phase_conductivity else {
            return Ok(None);
        };
        let mut derivatives = Vec::with_capacity(policy.element_laws.len());
        for (element, (vertices, law)) in self
            .storage
            .mesh
            .complex()
            .tets
            .iter()
            .zip(&policy.element_laws)
            .enumerate()
        {
            if element % ASSEMBLY_TILE == 0 {
                poll(self.cx, element)?;
            }
            let mut row = [0.0; 4];
            let contrast = law.derivative();
            if contrast != 0.0 {
                for (column, &vertex) in vertices.iter().enumerate() {
                    let vertex = vertex as usize;
                    let fraction_derivative = self
                        .storage
                        .curve_for_vertex(vertex)
                        .liquid_mass_fraction_derivative_at_specific_enthalpy(h[vertex])
                        .map_err(|source| EnthalpyError::Phase { vertex, source })?;
                    // Divide the larger factor first: a subnormal contrast
                    // can still produce a representable derivative on a
                    // narrow latent interval with a large fraction slope.
                    row[column] = finite(if contrast.abs() >= fraction_derivative {
                        (contrast / 4.0) * fraction_derivative
                    } else {
                        contrast * (fraction_derivative / 4.0)
                    })?;
                    if fraction_derivative != 0.0 && row[column] == 0.0 {
                        return Err(EnthalpyError::InvalidInput(
                            "nonzero phase conductivity derivative is not representable",
                        ));
                    }
                }
            }
            derivatives.push(row);
        }
        Ok(Some(derivatives))
    }

    fn assemble(
        &self,
        temperature: &[f64],
        phase_scales: Option<&[f64]>,
    ) -> Result<AssembledSystem, EnthalpyError> {
        Ok(assemble_operator_scaled_with_interfaces(
            self.cx,
            self.storage.mesh,
            self.problem.boundary,
            self.problem.material,
            self.problem.source,
            temperature,
            phase_scales,
            self.interfaces,
            self.problem.element_materials,
        )?)
    }

    fn stage(&self, h: &[f64]) -> Result<NewtonStage<'_, '_, '_, '_, '_>, EnthalpyError> {
        let temperature = self.temperatures(h)?;
        let radiation = self.pointwise_radiation
            .map(|radiation| radiation.assemble(self.cx, &temperature, h)).transpose()?;
        // The exact source moment pass returns values and derivatives together;
        // do not walk crossed source knots twice during a Newton linearization.
        let (phase_scales, phase_derivatives) = if let Some(source) = self.storage.sourced_phase.as_ref() {
            let (values, derivatives) = source.assemble(self.cx, self.storage.mesh, h, true)?;
            (Some(values), derivatives)
        } else {
            (self.phase_scales(h)?, self.phase_scale_derivatives(h)?)
        };
        let picard = self.assemble(&temperature, phase_scales.as_deref())?;
        let mut slope = Vec::with_capacity(h.len());
        let mut inverse = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            let derivative = self
                .storage
                .curve_for_vertex(vertex)
                .temperature_derivative_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            let diagonal = if self.is_prescribed(vertex) {
                self.storage.masses[vertex]
            } else {
                let transport_diagonal = match &radiation {
                    Some(radiation) => finite(picard.operator.get(vertex, vertex)
                        + radiation.preconditioner_diagonal(vertex))?,
                    None => picard.operator.get(vertex, vertex),
                };
                finite(self.storage.masses[vertex] + self.dt * transport_diagonal * derivative)?
            };
            if diagonal <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "positive enthalpy Picard diagonal required",
                ));
            }
            slope.push(derivative);
            inverse.push(finite(1.0 / diagonal)?);
        }
        let tangent = assemble_enthalpy_jacobian_with_optional_interfaces(
            self.cx,
            self.storage.mesh,
            self.problem.boundary,
            self.problem.material,
            &temperature,
            &slope,
            phase_scales.as_deref(),
            phase_derivatives.as_deref(),
            self.interfaces,
            self.problem.element_materials,
        )?;
        poll(self.cx, h.len())?;
        Ok(NewtonStage {
            context: self,
            source_j: None,
            tangent,
            slope,
            inverse,
            fatal: RefCell::new(None),
            trial_refusal: RefCell::new(None),
            radiation,
        })
    }

    fn residual(&self, h: &[f64], out: &mut [f64]) -> Result<(), EnthalpyError> {
        self.residual_with_source(h, out, None)
    }

    fn residual_with_source(
        &self,
        h: &[f64],
        out: &mut [f64],
        source_j: Option<&[fs_math::dd::Dd]>,
    ) -> Result<(), EnthalpyError> {
        let temperature = self.temperatures(h)?;
        let phase_scales = self.phase_scales(h)?;
        let system = self.assemble(&temperature, phase_scales.as_deref())?;
        self.assembled_residual_with_source(h, &temperature, &system, out, source_j)
    }

    fn assembled_residual(
        &self,
        h: &[f64],
        temperature: &[f64],
        system: &AssembledSystem,
        out: &mut [f64],
    ) -> Result<(), EnthalpyError> {
        self.assembled_residual_with_source(h, temperature, system, out, None)
    }

    fn assembled_residual_with_source(
        &self,
        h: &[f64],
        temperature: &[f64],
        system: &AssembledSystem,
        out: &mut [f64],
        source_j: Option<&[fs_math::dd::Dd]>,
    ) -> Result<(), EnthalpyError> {
        let radiation = self.pointwise_radiation
            .map(|radiation| radiation.load_and_report(self.cx, temperature, h)).transpose()?;
        for (vertex, value) in out.iter_mut().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            let (columns, entries) = system.operator.row(vertex);
            let transport: f64 = columns
                .iter()
                .zip(entries)
                .map(|(&column, &entry)| entry * temperature[column])
                .sum();
            let transport = match &radiation {
                Some((load, _)) => finite(transport + load[vertex])?,
                None => transport,
            };
            *value = match source_j {
                Some(source) => interference::physical_residual_j(
                    self.storage.masses[vertex],
                    h[vertex],
                    self.old[vertex],
                    self.dt,
                    transport,
                    system.load[vertex],
                    source[vertex],
                )?,
                None => finite(
                    self.storage.masses[vertex] * (h[vertex] - self.old[vertex])
                        + self.dt * (transport - system.load[vertex]),
                )?,
            };
        }
        for (index, &(vertex, temperature)) in self.problem.boundary.dirichlet().iter().enumerate()
        {
            if index % ASSEMBLY_TILE == 0 {
                poll(self.cx, index)?;
            }
            let prescribed = self.storage.prescribed_enthalpy(vertex, temperature)?;
            out[vertex] = finite(self.storage.masses[vertex] * (h[vertex] - prescribed))?;
        }
        Ok(())
    }

    fn finish(
        &self,
        h: Vec<f64>,
        newton: NewtonReport,
        tolerance_j: f64,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.check_prescribed_values(&h)?;
        let temperature = self.temperatures(&h)?;
        let phase_scales = self.phase_scales(&h)?;
        let system = self.assemble(&temperature, phase_scales.as_deref())?;
        let mut liquid_mass_fraction = Vec::with_capacity(h.len());
        let mut storage = 0.0;
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            storage = finite(storage + self.storage.masses[vertex] * (value - self.old[vertex]))?;
            liquid_mass_fraction.push(
                self.storage
                    .curve_for_vertex(vertex)
                    .state_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .liquid_mass_fraction(),
            );
        }
        // Independent physical source/face integration, not the sum of Newton
        // residuals or an inferred heat input from storage. This shared owner
        // pass and its contact reporting are bounded by the admitted mesh.
        let dofs = DofMap::new_allow_fully_prescribed(self.problem.boundary, h.len())?;
        let (energy, robin_fluxes) = energy_balance(
            self.storage.mesh,
            self.problem.boundary,
            self.problem.source,
            &system,
            &dofs,
            &temperature,
        );
        poll(self.cx, h.len())?;
        let dirichlet_in_w = self.prescribed_heat_input_w(&h, &temperature, &system)?;
        let pointwise_radiation = self.pointwise_radiation
            .map(|radiation| radiation.load_and_report(self.cx, &temperature, &h).map(|(_, report)| report)).transpose()?;
        let net =
            finite(energy.source_w + dirichlet_in_w - energy.neumann_out_w - energy.robin_out_w)?;
        let net = match &pointwise_radiation {
            Some(radiation) => finite(net - radiation.outward_heat_w())?,
            None => net,
        };
        let energy_residual_j = finite(storage - self.dt * net)?;
        if energy_residual_j.abs() > tolerance_j {
            return Err(EnthalpyError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j,
            });
        }
        let contact_fluxes = match self.interfaces {
            Some(interfaces) => interfaces.fluxes(&temperature)?,
            None => Vec::new(),
        };
        poll(self.cx, h.len())?;
        let (pointwise_radiation, pointwise_enclosure_radiation) = pointwise_radiation
            .map_or((None, None), pointwise::PointwiseBoundaryReport::into_reports);
        Ok(EnthalpyStepSolution {
            specific_enthalpy_j_kg: h,
            temperature,
            liquid_mass_fraction,
            source_w: energy.source_w,
            neumann_out_w: energy.neumann_out_w,
            robin_out_w: energy.robin_out_w,
            pointwise_radiation,
            pointwise_enclosure_radiation,
            dirichlet_in_w,
            stored_energy_change_j: storage,
            energy_residual_j,
            robin_fluxes,
            contact_fluxes,
            newton,
        })
    }
}

// Private adapter valid for precisely one shared Newton attempt. All line-
// search residuals evaluate the full nonlinear model; its Jacobian and Jacobi
// are the exact current Newton point's operators, reused across Krylov columns.
struct NewtonStage<'s, 'a, 'm, 'c, 'cx> {
    context: &'s StepContext<'a, 'm, 'c, 'cx>,
    source_j: Option<&'s [fs_math::dd::Dd]>,
    tangent: Csr,
    slope: Vec<f64>,
    inverse: Vec<f64>,
    fatal: RefCell<Option<EnthalpyError>>,
    trial_refusal: RefCell<Option<EnthalpyError>>,
    radiation: Option<PointwiseBoundaryAssembly>,
}

impl NewtonStage<'_, '_, '_, '_, '_> {
    fn take_failure(&self) -> Option<EnthalpyError> {
        self.fatal.take().or_else(|| self.trial_refusal.take())
    }

    fn poll_action(&self, at: usize) -> bool {
        if self.fatal.borrow().is_some() {
            return false;
        }
        if let Err(error) = poll(self.context.cx, at) {
            self.fatal.replace(Some(error));
            return false;
        }
        true
    }
}

impl NonlinearProblem for NewtonStage<'_, '_, '_, '_, '_> {
    fn dimension(&self) -> usize {
        self.slope.len()
    }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.poll_action(0) {
            return;
        }
        match self.context.residual_with_source(x, output, self.source_j) {
            Ok(()) => {
                self.trial_refusal.take();
            }
            Err(error) => {
                // Out-of-domain trial states are rejected by the shared line
                // search. Constitutive/assembly cancellation is always fatal.
                let recoverable = matches!(
                    error,
                    EnthalpyError::Phase {
                        source: PhaseStateError::OutsideEnthalpyDomain { .. },
                        ..
                    } | EnthalpyError::SourceEnthalpyOutsideDomain { .. }
                        | EnthalpyError::Conduction(ConductionError::OutsideTemperatureSpan { .. })
                ) || matches!(
                    &error,
                    EnthalpyError::Conduction(ConductionError::MaterialQuery { property, .. })
                        if self.context.pointwise_radiation.is_some()
                            && property == crate::SURFACE_EMISSIVITY_PROPERTY
                );
                if recoverable {
                    self.trial_refusal.replace(Some(error));
                } else {
                    self.fatal.replace(Some(error));
                }
                output.fill(f64::NAN);
            }
        }
    }

    fn jacobian_apply(&self, _x: &[f64], direction: &[f64], output: &mut [f64]) {
        for (row, value) in output.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 && !self.poll_action(row) {
                output.fill(f64::NAN);
                return;
            }
            if self.context.is_prescribed(row) {
                *value = self.context.storage.masses[row] * direction[row];
                continue;
            }
            let (columns, entries) = self.tangent.row(row);
            let applied: f64 = columns
                .iter()
                .zip(entries)
                .map(|(&column, &entry)| entry * direction[column])
                .sum();
            *value = self.context.storage.masses[row] * direction[row] + self.context.dt * applied;
        }
        if let Some(radiation) = &self.radiation
            && let Err(error) = radiation.add_scaled_tangent(
                self.context.cx, direction, Some(&self.slope), self.context.dt,
                |row| !self.context.is_prescribed(row), output,
            )
        {
            self.fatal.replace(Some(error.into()));
            output.fill(f64::NAN);
        }
    }

    fn preconditioner_apply(
        &self,
        _x: &[f64],
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        for (row, value) in output.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 && !self.poll_action(row) {
                output.fill(f64::NAN);
                return;
            }
            *value = residual[row] * self.inverse[row];
        }
    }
}

fn validate_work(n: usize, dt: f64, config: EnthalpyStepConfig) -> Result<(), EnthalpyError> {
    if !dt.is_finite()
        || dt <= 0.0
        || config.max_newton_iterations == 0
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j <= 0.0
    {
        return Err(EnthalpyError::InvalidInput(
            "positive finite step, energy tolerance and Newton budget required",
        ));
    }
    let restart = config.newton.linear_restart;
    // Conservative numerical-buffer size: Arnoldi and flexible bases,
    // Hessenberg, rotations, nonlinear trials and physical nodal scratch.
    restart
        .checked_add(1)
        .and_then(|r| r.checked_mul(restart))
        .and_then(|h| {
            n.checked_mul(restart.checked_mul(2)?.checked_add(32)?)?
                .checked_add(h)
        })
        .and_then(|v| v.checked_add(restart.checked_mul(8)?))
        .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or(EnthalpyError::InvalidInput(
            "Newton workspace size overflow",
        ))?;
    config
        .max_newton_iterations
        .checked_mul(restart)
        .and_then(|v| v.checked_mul(config.newton.max_linear_cycles))
        .ok_or(EnthalpyError::InvalidInput("Newton work budget overflow"))?;
    Ok(())
}

fn finite(value: f64) -> Result<f64, EnthalpyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EnthalpyError::NonFiniteArithmetic)
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), EnthalpyError> {
    cx.checkpoint().map_err(|_| {
        EnthalpyError::Conduction(ConductionError::Cancelled {
            stage: "spatial-enthalpy",
            at,
        })
    })
}
