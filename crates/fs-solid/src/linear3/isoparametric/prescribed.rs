//! P1 Dirichlet lifting using linear3's original element, solve and stress owners.
//! The private assembly module keeps those implementation details private; the
//! public methods attach directly to TetLinearElasticProblem, not to HEX/TET10.
use super::super::*;

fn invalid(what: &'static str) -> TetElasticError {
    TetElasticError::InvalidStaticSolveConfig { what }
}

impl TetLinearElasticProblem<'_> {
    /// Solve `K_ff u_f = f_f + f_thermal,f - K_fc g_c` with nonzero support
    /// displacements in metres. Entries are strictly increasing full xyz DOFs,
    /// each already present in `self.fixed_dofs`; omitted fixed DOFs stay zero.
    /// This call never changes which components are constrained. Prescriptions
    /// are absolute displacements from the reference mesh, not increments.
    ///
    /// Thermal strain and dead nodal loads retain their existing meaning. The
    /// same P1 element matrices form the lift; the original bounded PCG owns the
    /// reduced solve, including rigid-mode refusal. The reported true relative
    /// residual is normalized by the LIFTED right-hand side. Recovery through
    /// `recover_thermal_stress_with_prescribed_dofs` uses actual total motion.
    /// Fully prescribed meshes require no matrix solve. No thermal surrogate,
    /// penalty stiffness, follower force or transient inertia is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_thermal_displacement_with_prescribed_dofs(
        &self, thermal: TetThermalStrainField<'_>, nodal_forces_n: &[(usize, f64)],
        prescribed_m: &[(usize, f64)], maximum_prescribed_displacement_m: f64,
        config: TetStaticSolveConfig, cx: &Cx<'_>,
    ) -> Result<TetThermalDisplacementSolution, TetElasticError> {
        let (mut load, g) = self.prescribed_load(thermal, nodal_forces_n, prescribed_m,
            maximum_prescribed_displacement_m, cx)?;
        if !load.free_dofs.is_empty() && prescribed_m.iter().any(|(_, value)| *value != 0.0) {
            let mut reduced_of = vec![None; g.len()];
            for (i, &dof) in load.free_dofs.iter().enumerate() {
                if i % 256 == 0 { cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?; }
                reduced_of[dof] = Some(i);
            }
            for (element, tet) in self.tetrahedra.iter().enumerate() {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                let local: [f64; 12] = core::array::from_fn(|d| g[3 * tet[d / 3] + d % 3]);
                if local.iter().all(|value| *value == 0.0) { continue; }
                let geometry = TetGeometry::try_new(element, tet.map(|node| self.nodes_m[node]),
                    self.budget.minimum_scaled_jacobian)?;
                let (ke, _) = element_matrices(&geometry, self.materials.at(element));
                for row in 0..12 {
                    let Some(i) = reduced_of[3 * tet[row / 3] + row % 3] else { continue; };
                    let lift = (0..12).fold(0.0, |sum, col| ke[row][col].mul_add(local[col], sum));
                    load.reduced_force_n[i] -= lift;
                    if !lift.is_finite() || !load.reduced_force_n[i].is_finite() {
                        return Err(TetElasticError::NonFiniteThermalLoad {
                            element, local_dof: Some(row),
                        });
                    }
                }
            }
        }
        // This temporary load is never published as a physical force artifact:
        // only its reduced RHS includes -K_fc g. Its identity binds physical f
        // and the complete prescription, identically for subsequent recovery.
        let mut solution = self.solve_thermal_displacement_from_load(load, config, cx)?;
        for &(dof, _) in prescribed_m {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            solution.displacement_m[dof / 3][dof % 3] = g[dof];
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-prescribed-solve.v1");
        identity.update(solution.identity.as_bytes());
        for (node, values) in solution.displacement_m.iter().enumerate() {
            if node % 256 == 0 { cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?; }
            for value in values { identity.update(&value.to_bits().to_le_bytes()); }
        }
        solution.identity = identity.finalize();
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(solution)
    }

    /// Recover `C (B u - epsilon_thermal)`, stored energy and nodal internal
    /// forces from the SAME load, material and prescribed-support state. A
    /// different prescription, including its admission bound, refuses rather
    /// than rebinding a solved field. All strain/frame checks are the original
    /// stress owner's. The original geometry-update method accepts this full u.
    #[allow(clippy::too_many_arguments)]
    pub fn recover_thermal_stress_with_prescribed_dofs(
        &self, solution: &TetThermalDisplacementSolution, thermal: TetThermalStrainField<'_>,
        nodal_forces_n: &[(usize, f64)], prescribed_m: &[(usize, f64)],
        maximum_prescribed_displacement_m: f64, mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64, cx: &Cx<'_>,
    ) -> Result<TetThermalStressReport, TetElasticError> {
        let (load, _) = self.prescribed_load(thermal, nodal_forces_n, prescribed_m,
            maximum_prescribed_displacement_m, cx)?;
        self.recover_thermal_stress_from_load(solution, EigenstrainFieldView::Thermal(thermal), load, mesh_frame,
            maximum_absolute_mandel_strain, cx)
    }

    fn prescribed_load(
        &self, thermal: TetThermalStrainField<'_>, nodal_forces_n: &[(usize, f64)],
        prescribed_m: &[(usize, f64)], maximum: f64, cx: &Cx<'_>,
    ) -> Result<(TetThermalLoad, Vec<f64>), TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        if !maximum.is_finite() || maximum <= 0.0 {
            return Err(invalid("maximum prescribed displacement must be finite and positive"));
        }
        if prescribed_m.len() > self.fixed_dofs.len() {
            return Err(invalid("prescribed entries exceed the fixed-DOF count"));
        }
        // The canonical owner admits all mesh, material, thermal, fixed-DOF and
        // load data before g allocation or connectivity indexing below.
        let mut load = self.assemble_thermal_load_with_nodal_forces(thermal, nodal_forces_n, cx)?;
        let mut g = vec![0.0; load.full_equivalent_force_n.len()];
        let mut previous = None;
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-prescribed-load.v1");
        identity.update(load.identity.as_bytes());
        identity.update(&maximum.to_bits().to_le_bytes());
        hash_usize(&mut identity, prescribed_m.len())?;
        for &(dof, value) in prescribed_m {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            if dof >= g.len() {
                return Err(TetElasticError::FixedDofOutOfRange { dof, full_dofs: g.len() });
            }
            if previous.is_some_and(|last| dof <= last) {
                return Err(invalid("prescribed DOFs must be strictly increasing and unique"));
            }
            if load.free_dofs.binary_search(&dof).is_ok() {
                return Err(invalid("prescribed displacement DOF is not a declared fixed component"));
            }
            if !value.is_finite() || value.abs() > maximum {
                return Err(invalid("prescribed displacement is nonfinite or exceeds its metre bound"));
            }
            g[dof] = value;
            previous = Some(dof);
            hash_usize(&mut identity, dof)?;
            identity.update(&value.to_bits().to_le_bytes());
        }
        load.identity = identity.finalize();
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok((load, g))
    }
}
