//! Source-bound thermal loads and equilibrium on the HEX8 strain formulation.

use super::{
    Hex8ElasticProblem, Hex8PrestressedAssembly, Hex8ThermalResponse, isoparametric, tensor_rule,
};
use crate::linear3::{
    TetElasticError, TetStaticSolveConfig, TetThermalLoad, TetThermalStrainField,
    add_nodal_forces_to_thermal_load,
};
use fs_blake3::ContentHash;
use fs_exec::Cx;

impl Hex8ElasticProblem<'_> {
    /// Assemble integral B^T C epsilon_free on the same Standard or B-bar
    /// columns as stiffness. Full rows include all supports; reduced rows use
    /// the original pencil's DOF map. The state must bind the exact current
    /// elastic material, including its full mesh-frame tensor.
    ///
    /// Integrated expansion is uniform or per element. This supports distinct
    /// temperatures and anisotropic expansion across elements, without defining
    /// an intra-element thermal-pressure projection. All density and geometry
    /// remain in the original reference configuration.
    pub fn assemble_thermal_load(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        self.validate(false, cx)?;
        let stiffness_rule = tensor_rule(self.stiffness_points_per_dir);
        let extra_mass_rule = (self.stiffness_points_per_dir < 3).then(|| tensor_rule(3));
        let mass_rule = extra_mass_rule.as_deref().unwrap_or(&stiffness_rule);
        isoparametric::thermal_load_with_formulation::<8, 24, _>(
            self.nodes_m,
            self.hexahedra,
            self.materials,
            thermal_strains.into(),
            self.fixed_dofs,
            self.budget,
            self.formulation,
            &stiffness_rule,
            mass_rule,
            |nodes, reference, element, point| self.point(nodes, reference, element, point),
            cx,
        )
    }

    /// Add sorted, unique full node-major xyz dead loads [N]. Forces on fixed
    /// DOFs retain their support-reaction meaning and do not create motion.
    pub fn assemble_thermal_load_with_nodal_forces(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        nodal_forces_n: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        add_nodal_forces_to_thermal_load(
            self.assemble_thermal_load(thermal_strains, cx)?,
            self.nodes_m.len(),
            nodal_forces_n,
            cx,
        )
    }

    /// Solve thermal equilibrium through the existing checked SPD owner and
    /// recover stress, storage and full support forces on the same columns.
    /// B-bar observations use its assumed strain; compatible and assumed
    /// strains independently obey the supplied finite small-strain ceiling.
    ///
    /// Every connected component must have its six rigid motions restrained.
    /// A fully fixed solid still returns its physical stress and reactions.
    /// This is frozen-reference linear equilibrium, without geometric
    /// prestress stiffness, changing temperature or finite-rotation kinematics.
    pub fn solve_thermal_response(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        config: TetStaticSolveConfig,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Hex8ThermalResponse, TetElasticError> {
        self.solve_thermal_response_impl(
            thermal_strains,
            None,
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// The same thermal equilibrium with explicit full dead nodal forces.
    /// Returned force balance is internal minus applied force: support
    /// reactions on fixed DOFs and the actual residual on free DOFs.
    pub fn solve_thermal_response_with_nodal_forces(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        nodal_forces_n: &[(usize, f64)],
        config: TetStaticSolveConfig,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Hex8ThermalResponse, TetElasticError> {
        self.solve_thermal_response_impl(
            thermal_strains,
            Some(nodal_forces_n),
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// Solve thermal equilibrium and assemble its signed initial-stress pencil.
    ///
    /// The same compatible or B-bar stress accepted by equilibrium enters
    /// `K_g[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma grad(N_b) dV`.
    /// Actual quadrature stress is retained, including anisotropic coupling and
    /// within-element variation; an element-mean stress is not substituted.
    /// The original material stiffness, consistent mass, reference geometry
    /// and homogeneous incremental supports feed the existing spectral owner.
    ///
    /// This conservative initial-stress approximation is not a finite-strain
    /// tangent or stability certificate. Compression may make the tangent
    /// indefinite; nonpositive eigenvalues must remain signed. A fully fixed
    /// body has a static response but no vibration pencil and returns NoFreeDofs.
    pub fn solve_thermal_tangent(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        config: TetStaticSolveConfig,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Hex8PrestressedAssembly, TetElasticError> {
        self.solve_thermal_tangent_impl(
            thermal_strains,
            None,
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// The same signed tangent about thermal and dead-nodal-force equilibrium.
    ///
    /// Solve and recover from these exact material, thermal and force inputs
    /// before assembling the original geometric-stiffness owner. Fixed-DOF
    /// loads retain their support-reaction meaning. Dead forces add no follower
    /// stiffness; no independently supplied or stale stress field is accepted.
    pub fn solve_thermal_tangent_with_nodal_forces(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        nodal_forces_n: &[(usize, f64)],
        config: TetStaticSolveConfig,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Hex8PrestressedAssembly, TetElasticError> {
        self.solve_thermal_tangent_impl(
            thermal_strains,
            Some(nodal_forces_n),
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    fn solve_thermal_tangent_impl(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        nodal_forces_n: Option<&[(usize, f64)]>,
        config: TetStaticSolveConfig,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Hex8PrestressedAssembly, TetElasticError> {
        let response = self.solve_thermal_response_impl(
            thermal_strains,
            nodal_forces_n,
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )?;
        let pencil = self.assemble(cx)?;
        isoparametric::thermal_response::prestress::assemble::<8, 24, _>(
            self.nodes_m,
            self.hexahedra,
            self.materials,
            thermal_strains.into(),
            pencil,
            response,
            &tensor_rule(self.stiffness_points_per_dir),
            |nodes, reference, element, point| self.point(nodes, reference, element, point),
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    fn solve_thermal_response_impl(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        nodal_forces_n: Option<&[(usize, f64)]>,
        config: TetStaticSolveConfig,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<Hex8ThermalResponse, TetElasticError> {
        self.validate(false, cx)?;
        let stiffness_rule = tensor_rule(self.stiffness_points_per_dir);
        let extra_mass_rule = (self.stiffness_points_per_dir < 3).then(|| tensor_rule(3));
        let mass_rule = extra_mass_rule.as_deref().unwrap_or(&stiffness_rule);
        isoparametric::thermal_response::solve_with_formulation::<8, 24, _>(
            self.nodes_m,
            self.hexahedra,
            self.materials,
            thermal_strains.into(),
            nodal_forces_n,
            self.fixed_dofs,
            self.budget,
            self.formulation,
            &stiffness_rule,
            mass_rule,
            |nodes, reference, element, point| self.point(nodes, reference, element, point),
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }
}
