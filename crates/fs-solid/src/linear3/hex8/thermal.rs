//! Source-bound thermal loads and equilibrium on the HEX8 strain formulation.

use super::{Hex8ElasticProblem, Hex8ThermalResponse, isoparametric, tensor_rule};
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
