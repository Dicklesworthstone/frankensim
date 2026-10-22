//! Reuse the EXACT accepted thermal operator inputs, not a second source query.
use fs_exec::Cx;
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState};
use fs_solid::prestress::PrestressedTetAssembly;
use super::{ThermalBody, ThermalBodyError, ThermalSolveConfig};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Linearization {
    pub reference_nodes_m: Vec<[f64; 3]>,
    pub materials: Vec<TetElasticMaterial>,
    pub strains: Vec<TetThermalStrainState>,
    pub fixed_dofs: Vec<usize>,
    pub nodal_forces_n: Vec<(usize, f64)>,
    pub config: ThermalSolveConfig,
}
impl ThermalBody<'_> {
    /// Assemble the initial-stress vibration pencil of THIS accepted thermal state.
    ///
    /// Current source-resolved E/nu, reference density, thermal strains, mesh,
    /// supports, mechanical dead forces and actual displacement are retained at
    /// the successful shape solve. No new card lookup, density query, stress multiplier or hand-entered
    /// modal coefficient is needed. The original solid owner rechecks solution
    /// binding and supplies the ordinary (K_e + K_g, M) carrier for fs-modal.
    /// The source-backed heat path exposes this through
    /// `ConstrainedHeatedBody::configuration()` at every accepted endpoint.
    ///
    /// This is explicitly the small-displacement REFERENCE-mesh initial-stress
    /// approximation. It is not assembly on independently rescaled elements,
    /// a full deformed-configuration tangent, a new equilibrium solve, or live
    /// vibration-state transfer. Mass and geometry updates remain available in
    /// `body()` and `geometry()` with their separate meanings. Followers and
    /// large-rotation effects need another linearization, not frequency retuning.
    ///
    /// The original assembly/strain budgets remain in force. An all-fixed body
    /// has no vibration DOFs and refuses, rather than silently releasing supports.
    /// Compressive/unstable eigenvalues remain signed; the caller must not turn
    /// a nonpositive eigenvalue into a real acoustic frequency.
    pub fn initial_stress_pencil(&self, cx: &Cx<'_>)
        -> Result<PrestressedTetAssembly, ThermalBodyError> {
        let input = &self.linearization;
        let problem = TetLinearElasticProblem {
            nodes_m: &input.reference_nodes_m,
            tetrahedra: self.geometry.tetrahedra(),
            materials: TetMaterialField::PerElement(&input.materials),
            fixed_dofs: &input.fixed_dofs,
            budget: input.config.assembly,
        };
        let field = TetThermalStrainField::PerElement(&input.strains);
        if input.nodal_forces_n.is_empty() {
            Ok(problem.assemble_thermal_tangent(&self.displacement, field,
                input.config.mesh_frame, input.config.maximum_absolute_mandel_strain, cx)?)
        } else {
            Ok(problem.assemble_thermal_tangent_with_nodal_forces(&self.displacement,
                field, &input.nodal_forces_n, input.config.mesh_frame,
                input.config.maximum_absolute_mandel_strain, cx)?)
        }
    }
}
