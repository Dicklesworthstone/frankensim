//! Thermal equilibrium and initial stress in the existing valve reduction.
//!
//! Reference stiffness, mass, thermal strain, equilibrium and geometric
//! stiffness remain in the original solid owners. The common reducer uses the
//! resulting tangent for BOTH its modal solve and static pressure compliance.
//! Only aperture observation uses the accepted displaced geometry.

/// Retained pressure-driven motion and source-mesh modal state transfer.
pub mod motion;

use super::{
    PreparedValve, PressureLoadedValve, ReduceError, SingleModeValveAdmission,
    ValveCardRequest, ValveElementOrder, material_field, validate_element_topology,
};
use crate::linear3::{
    QuadraticPrestressedAssembly, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig,
    TetThermalDisplacementSolution, TetThermalStrainField,
    tet10::{QuadraticTetProblem, tet10_point},
};
use crate::prestress::PrestressedTetAssembly;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;

/// Explicit equilibrium solve and small-displacement applicability limits.
#[derive(Clone, Copy, Debug)]
pub struct ThermalValveControls {
    /// Original bounded static-solve controls.
    pub solve: TetStaticSolveConfig,
    /// Nonzero frame of the reference mesh and all Mandel components.
    pub mesh_frame: ContentHash,
    /// Original total/elastic Mandel component ceiling.
    pub maximum_absolute_mandel_strain: f64,
    /// Maximum Euclidean nodal displacement before aperture observation [m].
    /// This caller-declared limit does not establish finite-strain accuracy.
    pub maximum_nodal_displacement_m: f64,
}

/// Complete accepted source-bound equilibrium, without a second stress law.
#[derive(Clone, Debug)]
pub enum ThermalValveEquilibrium {
    /// Affine tetrahedra retain both their displacement and signed tangent.
    Tet4 {
        /// Original thermal/dead-force equilibrium.
        displacement: TetThermalDisplacementSolution,
        /// Original initial-stress assembly and actual element stress.
        tangent: PrestressedTetAssembly,
    },
    /// Quadratic tetrahedra retain actual quadrature stress in their tangent.
    Tet10(QuadraticPrestressedAssembly),
}

impl ThermalValveEquilibrium {
    /// Compatible displacement at every reference node [m].
    #[must_use]
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        match self {
            Self::Tet4 { displacement, .. } => displacement.displacement_m(),
            Self::Tet10(tangent) => tangent.response().displacement_m(),
        }
    }

    /// Reference-mass pencil with signed initial-stress stiffness.
    #[must_use]
    pub const fn pencil(&self) -> &TetElasticAssembly {
        match self {
            Self::Tet4 { tangent, .. } => tangent.pencil(),
            Self::Tet10(tangent) => tangent.pencil(),
        }
    }

    /// Original identity binding equilibrium to geometry, loads and material.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        match self {
            Self::Tet4 { tangent, .. } => tangent.identity(),
            Self::Tet10(tangent) => tangent.identity(),
        }
    }
}

/// Thermally preloaded valve, its accepted geometry and complete source chain.
/// The underlying card codec remains unchanged; persist this thermal receipt
/// separately when retaining the equilibrium inputs is required.
#[derive(Clone, Debug)]
pub struct ThermalPressureLoadedValve {
    valve: PressureLoadedValve,
    equilibrium: ThermalValveEquilibrium,
    observation_nodes_m: Vec<[f64; 3]>,
    mesh_frame: ContentHash,
    identity: ContentHash,
}

impl ThermalPressureLoadedValve {
    /// Existing pressure-port/card owner, with tangent modes and compliance.
    #[must_use]
    pub const fn valve(&self) -> &PressureLoadedValve {
        &self.valve
    }

    /// Original equilibrium and stress, rather than an authored preload scalar.
    #[must_use]
    pub const fn equilibrium(&self) -> &ThermalValveEquilibrium {
        &self.equilibrium
    }

    /// Accepted displaced nodes used to measure gap, width and thickness.
    /// The vibration pencil retains reference mass and reference gradients.
    #[must_use]
    pub fn observation_nodes_m(&self) -> &[[f64; 3]] {
        &self.observation_nodes_m
    }

    /// Coordinate frame of the retained source mesh and tensor fields.
    #[must_use]
    pub const fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }

    /// Binds the valve, equilibrium, displacement limit and observed geometry.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Compile the same physical pressure port after original first-mode
    /// compliance admission. There is no independently tuned reed image.
    pub fn beating_reed_checked(
        &self,
        blowing_pressure_pa: f64,
        attack_s: f64,
        admission: SingleModeValveAdmission,
    ) -> Result<fs_scenario::BeatingReed, ReduceError> {
        self.valve
            .beating_reed_checked(blowing_pressure_pa, attack_s, admission)
    }
}

/// Solve equilibrium and reduce its thermal/dead-load tangent to a valve.
///
/// Material tensors and total thermal eigenstrains are already resolved by
/// their original source owners. The same reference mesh and homogeneous
/// supports govern equilibrium and incremental vibration. The original mass
/// is conserved; density is not independently rescaled after expansion.
///
/// `equilibrium_nodal_forces_n` uses the original sorted full-DOF dead-force
/// convention. `nodal_load_per_pa_m2` is a separate, caller-declared FIXED
/// pressure-load pattern on the reference nodes. Its work-conjugate projected
/// areas use the actual tangent modes. No follower stiffness, updated surface
/// normal, pressure traction or pressure-dependent static equilibrium is
/// inferred from this field.
///
/// The existing displaced-geometry checks admit the aperture observation.
/// Quadratic maps use the original stiffness AND mass quadrature points; this
/// does not certify injectivity between them. The opposing plane and observed
/// axes stay fixed. Contact/closed gaps, excessive displacement or strain, and
/// mismatched thermal/material states refuse rather than selecting a contact model.
///
/// The whole tangent must pass the existing LDLT positive-inertia test. A
/// positive requested frequency window cannot conceal an unstable direction
/// elsewhere. Exact static pressure/opening compliance and modal truncation
/// admission use that SAME tangent. This is the conservative small-displacement
/// initial-stress approximation, not finite-strain or live modal-state transfer.
///
/// Nothing is published on equilibrium, geometry, modal, compliance or
/// cancellation refusal. Loss remains the request's explicitly authored scalar.
#[allow(clippy::too_many_arguments)]
pub fn reduce_pressure_loaded_valve(
    request: &ValveCardRequest<'_, TetMaterialField<'_>>,
    thermal: TetThermalStrainField<'_>,
    equilibrium_nodal_forces_n: &[(usize, f64)],
    nodal_load_per_pa_m2: &[[f64; 3]],
    controls: ThermalValveControls,
    cx: &Cx<'_>,
) -> Result<ThermalPressureLoadedValve, ReduceError> {
    poll(cx)?;
    if !controls.maximum_nodal_displacement_m.is_finite()
        || controls.maximum_nodal_displacement_m <= 0.0
    {
        return Err(TetElasticError::InvalidGeometryUpdate {
            what: "maximum_nodal_displacement_m must be finite and positive",
        }
        .into());
    }
    validate_element_topology(request)?;
    let material_identity =
        material_field::field_identity(request.material, request.tetrahedra.len())?;
    let budget = TetAssemblyBudget::default();
    let (equilibrium, observation_nodes_m) = match request.element_order {
        ValveElementOrder::Tet4 => {
            let problem = TetLinearElasticProblem {
                nodes_m: request.nodes_m,
                tetrahedra: request.tetrahedra,
                materials: request.material,
                fixed_dofs: request.fixed_dofs,
                budget,
            };
            let displacement = problem.solve_thermal_displacement_with_nodal_forces(
                thermal,
                equilibrium_nodal_forces_n,
                controls.solve,
                cx,
            )?;
            let tangent = problem.assemble_thermal_tangent_with_nodal_forces(
                &displacement,
                thermal,
                equilibrium_nodal_forces_n,
                controls.mesh_frame,
                controls.maximum_absolute_mandel_strain,
                cx,
            )?;
            let geometry = problem.update_geometry_from_displacement(
                &displacement,
                controls.maximum_nodal_displacement_m,
                cx,
            )?;
            (
                ThermalValveEquilibrium::Tet4 {
                    displacement,
                    tangent,
                },
                geometry.nodes_m().to_vec(),
            )
        }
        ValveElementOrder::Tet10 {
            elements,
            quadrature,
        } => {
            let problem = QuadraticTetProblem {
                nodes_m: request.nodes_m,
                elements,
                materials: request.material,
                fixed_dofs: request.fixed_dofs,
                budget,
                quadrature,
            };
            let tangent = problem.solve_thermal_tangent_with_nodal_forces(
                thermal,
                equilibrium_nodal_forces_n,
                controls.solve,
                controls.mesh_frame,
                controls.maximum_absolute_mandel_strain,
                cx,
            )?;
            let nodes = displaced_quadratic_nodes(
                &problem,
                tangent.response().displacement_m(),
                controls.maximum_nodal_displacement_m,
                cx,
            )?;
            (ThermalValveEquilibrium::Tet10(tangent), nodes)
        }
    };
    poll(cx)?;
    let valve = super::reduce_pressure_loaded_valve_impl(
        request,
        request.material,
        material_identity,
        nodal_load_per_pa_m2,
        cx,
        Some(PreparedValve {
            assembly: equilibrium.pencil(),
            observation_nodes_m: &observation_nodes_m,
            static_stiffness: None,
            material_loss: None,
            complete_basis: false,
        }),
    )?;
    let mut identity =
        DomainHasher::new("org.frankensim.fs-solid.thermal-pressure-loaded-valve.v1");
    identity.update(valve.identity().as_bytes());
    identity.update(equilibrium.identity().as_bytes());
    identity.update(&controls.maximum_nodal_displacement_m.to_bits().to_le_bytes());
    for (i, node) in observation_nodes_m.iter().enumerate() {
        if i % 4096 == 0 {
            poll(cx)?;
        }
        for value in node {
            identity.update(&value.to_bits().to_le_bytes());
        }
    }
    poll(cx)?;
    Ok(ThermalPressureLoadedValve {
        valve,
        equilibrium,
        observation_nodes_m,
        mesh_frame: controls.mesh_frame,
        identity: identity.finalize(),
    })
}

fn displaced_quadratic_nodes(
    problem: &QuadraticTetProblem<'_>,
    displacement: &[[f64; 3]],
    maximum_displacement_m: f64,
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>, ReduceError> {
    let mut nodes = Vec::with_capacity(problem.nodes_m.len());
    for (i, (point, delta)) in problem.nodes_m.iter().zip(displacement).enumerate() {
        if i % 4096 == 0 {
            poll(cx)?;
        }
        let norm = delta[0].hypot(delta[1]).hypot(delta[2]);
        if !norm.is_finite() || norm > maximum_displacement_m {
            return Err(TetElasticError::InvalidGeometryUpdate {
                what: "nodal displacement exceeds its admitted bound",
            }
            .into());
        }
        let moved: [f64; 3] = core::array::from_fn(|axis| point[axis] + delta[axis]);
        if moved.iter().any(|value| !value.is_finite()) {
            return Err(TetElasticError::InvalidGeometryUpdate {
                what: "updated coordinate is non-finite",
            }
            .into());
        }
        nodes.push(moved);
    }
    let stiffness = duffy_quadrature(problem.quadrature.stiffness_points_per_dir);
    let mass = duffy_quadrature(problem.quadrature.mass_points_per_dir);
    for (element, connectivity) in problem.elements.iter().enumerate() {
        poll(cx)?;
        let points = connectivity.map(|node| nodes[node]);
        for (q, (lambda, _)) in stiffness.iter().chain(&mass).enumerate() {
            poll(cx)?;
            tet10_point(
                points,
                *lambda,
                problem.budget.minimum_scaled_jacobian,
                element,
                q,
            )?;
        }
    }
    Ok(nodes)
}

fn poll(cx: &Cx<'_>) -> Result<(), ReduceError> {
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    Ok(())
}

impl fs_scenario::reed_source::BeatingReedSource for ThermalPressureLoadedValve {
    type Admission = SingleModeValveAdmission;
    type Error = ReduceError;

    fn compile_beating_reed(
        &self,
        blowing_pressure_pa: f64,
        attack_s: f64,
        admission: &Self::Admission,
    ) -> Result<fs_scenario::BeatingReed, Self::Error> {
        self.beating_reed_checked(blowing_pressure_pa, attack_s, *admission)
    }
}
