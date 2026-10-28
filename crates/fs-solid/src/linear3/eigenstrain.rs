//! Generic initial-strain entry to the original tetrahedral solid kernels.
//!
//! Hygroscopic, manufacturing or other source-integrated small strains do not
//! acquire a fictitious temperature. Thermal inputs keep their own temperature
//! and source-path semantics; both paths use one load, solve, stress and geometry
//! implementation. This adapter does not supply a chemical/thermal free energy.

use super::*;

/// Source-integrated symmetric engineering strain in the mesh Mandel frame.
/// The source identity owns the physical datum and integration convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TetEigenstrainState {
    free_strain_mandel: [f64; 6],
    elastic_material_state_identity: [u8; 32],
    source_identity: ContentHash,
    maximum_absolute_mandel_strain: f64,
    identity: ContentHash,
}
impl TetEigenstrainState {
    /// Bind a finite small-strain source to the exact current elastic state.
    /// This never interprets strain as temperature or creates an energy source.
    pub fn try_new(
        free_strain_mandel: [f64; 6],
        elastic_material_state_identity: [u8; 32],
        source_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, TetElasticError> {
        if !maximum_absolute_mandel_strain.is_finite()
            || maximum_absolute_mandel_strain <= 0.0
            || free_strain_mandel
                .iter()
                .any(|x| !x.is_finite() || x.abs() > maximum_absolute_mandel_strain)
            || elastic_material_state_identity == [0; 32]
            || source_identity == ContentHash([0; 32])
        {
            return Err(TetElasticError::InvalidEigenstrain {
                what: "finite bounded strain and nonzero elastic/source identities required",
            });
        }
        let mut h = DomainHasher::new("org.frankensim.fs-solid.tet-eigenstrain.v1");
        h.update(&elastic_material_state_identity);
        h.update(source_identity.as_bytes());
        h.update(&maximum_absolute_mandel_strain.to_bits().to_le_bytes());
        for x in free_strain_mandel {
            h.update(&x.to_bits().to_le_bytes());
        }
        Ok(Self {
            free_strain_mandel,
            elastic_material_state_identity,
            source_identity,
            maximum_absolute_mandel_strain,
            identity: h.finalize(),
        })
    }

    /// Rotate three principal integrated strains with the existing solid
    /// second-order Mandel transform. Columns of Q are source axes in mesh
    /// coordinates; source and target frame labels must be nonzero. Equal frame
    /// labels require the identity rotation. Elastic symmetry is independent.
    #[allow(clippy::too_many_arguments)]
    pub fn try_from_principal_strains(
        material: &TetElasticMaterial,
        source_principal_strains: [f64; 3],
        principal_to_mesh: [[f64; 3]; 3],
        principal_frame: ContentHash,
        mesh_frame: ContentHash,
        source_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, TetElasticError> {
        validate_rotation(principal_to_mesh)?;
        let unit = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        if principal_frame == ContentHash([0; 32])
            || mesh_frame == ContentHash([0; 32])
            || (principal_frame == mesh_frame && principal_to_mesh != unit)
            || source_principal_strains.iter().any(|x| !x.is_finite())
        {
            return Err(TetElasticError::InvalidEigenstrain {
                what: "finite principal strain and consistent nonzero frames required",
            });
        }
        let transform = mandel_rotation(principal_to_mesh);
        let strain = core::array::from_fn(|i| {
            (0..3).fold(0.0, |a, j| {
                transform[i][j].mul_add(source_principal_strains[j], a)
            })
        });
        let mut result = Self::try_new(
            strain,
            material.material_state_identity,
            source_identity,
            maximum_absolute_mandel_strain,
        )?;
        let mut h = DomainHasher::new("org.frankensim.fs-solid.principal-eigenstrain.v1");
        h.update(result.identity.as_bytes());
        h.update(principal_frame.as_bytes());
        h.update(mesh_frame.as_bytes());
        for x in source_principal_strains
            .iter()
            .chain(principal_to_mesh.iter().flatten())
        {
            h.update(&x.to_bits().to_le_bytes());
        }
        result.identity = h.finalize();
        Ok(result)
    }
    /// Integrated mesh-frame strain, with sqrt(2)-scaled physical shear.
    pub const fn free_strain_mandel(&self) -> &[f64; 6] {
        &self.free_strain_mandel
    }
    /// Exact elastic state this strain loads.
    pub const fn elastic_material_state_identity(&self) -> [u8; 32] {
        self.elastic_material_state_identity
    }
    /// Original physical source/datum identity.
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    /// Authored small-strain ceiling.
    pub const fn maximum_absolute_mandel_strain(&self) -> f64 {
        self.maximum_absolute_mandel_strain
    }
    /// Complete state including orientation.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Uniform or per-tetrahedron initial strain, in original connectivity order.
#[derive(Clone, Copy, Debug)]
pub enum TetEigenstrainField<'a> {
    /// One source state throughout the body.
    Uniform(&'a TetEigenstrainState),
    /// One source state per tetrahedron.
    PerElement(&'a [TetEigenstrainState]),
}

/// Internal read-only view; physical source types remain distinct.
#[derive(Clone, Copy)]
pub(super) enum EigenstrainFieldView<'a> {
    Thermal(TetThermalStrainField<'a>),
    Generic(TetEigenstrainField<'a>),
}
pub(super) struct EigenstrainValue {
    pub(super) free_strain_mandel: [f64; 6],
    pub(super) elastic_material_state_identity: [u8; 32],
    pub(super) identity: ContentHash,
}
impl EigenstrainFieldView<'_> {
    pub(super) fn validate(&self, elements: usize) -> Result<(), TetElasticError> {
        match self {
            Self::Thermal(field) => field.validate(elements),
            Self::Generic(TetEigenstrainField::PerElement(states)) if states.len() != elements => {
                Err(TetElasticError::EigenstrainCountMismatch {
                    expected: elements,
                    got: states.len(),
                })
            }
            _ => Ok(()),
        }
    }
    pub(super) fn at(&self, element: usize) -> EigenstrainValue {
        match self {
            Self::Thermal(field) => {
                let state = field.at(element);
                EigenstrainValue {
                    free_strain_mandel: state.free_strain_mandel,
                    elastic_material_state_identity: state.elastic_material_state_identity,
                    identity: state.identity,
                }
            }
            Self::Generic(field) => {
                let state = match field {
                    TetEigenstrainField::Uniform(s) => *s,
                    TetEigenstrainField::PerElement(s) => &s[element],
                };
                EigenstrainValue {
                    free_strain_mandel: state.free_strain_mandel,
                    elastic_material_state_identity: state.elastic_material_state_identity,
                    identity: state.identity,
                }
            }
        }
    }
    pub(super) fn load_domain(&self) -> &'static str {
        match self {
            Self::Thermal(_) => "org.frankensim.fs-solid.tet-thermal-load.v1",
            Self::Generic(_) => "org.frankensim.fs-solid.tet-eigenstrain-load.v1",
        }
    }
    pub(super) fn material_mismatch(&self, element: usize) -> TetElasticError {
        match self {
            Self::Thermal(_) => TetElasticError::ThermalElasticStateMismatch { element },
            Self::Generic(_) => TetElasticError::EigenstrainElasticStateMismatch { element },
        }
    }
}

impl TetLinearElasticProblem<'_> {
    /// Assemble B^T C epsilon_free on the same mesh and supports as elasticity.
    pub fn assemble_eigenstrain_load(
        &self,
        strains: TetEigenstrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetEigenstrainLoad, TetElasticError> {
        self.assemble_eigenstrain_load_inner(EigenstrainFieldView::Generic(strains), cx)
    }
    /// Existing checked displacement solve, now accepting a nonthermal source.
    /// Failed source/material/solve/cancellation admission returns no state.
    pub fn solve_eigenstrain_displacement(
        &self,
        strains: TetEigenstrainField<'_>,
        config: TetStaticSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<TetEigenstrainDisplacementSolution, TetElasticError> {
        let load = self.assemble_eigenstrain_load(strains, cx)?;
        self.solve_thermal_displacement_from_load(load, config, cx)
    }
    /// Existing C(Bu-epsilon_free) stress/energy recovery with exact input binding.
    pub fn recover_eigenstrain_stress(
        &self,
        solution: &TetEigenstrainDisplacementSolution,
        strains: TetEigenstrainField<'_>,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<TetEigenstrainStressReport, TetElasticError> {
        let load = self.assemble_eigenstrain_load(strains, cx)?;
        self.recover_thermal_stress_from_load(
            solution,
            EigenstrainFieldView::Generic(strains),
            load,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// Assemble generic initial strain and explicit dead nodal forces together.
    ///
    /// Forces are newtons paired with full node-major xyz DOFs, strictly
    /// increasing and unique. Fixed-DOF forces remain in the full load for
    /// support reactions. All finite/range/ordering checks and load identities
    /// are shared with the original solid force owner. An empty force list
    /// preserves the original eigenstrain-only load identity exactly.
    pub fn assemble_eigenstrain_load_with_nodal_forces(
        &self,
        strains: TetEigenstrainField<'_>,
        nodal_forces_n: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<TetEigenstrainLoad, TetElasticError> {
        let load = self.assemble_eigenstrain_load(strains, cx)?;
        if nodal_forces_n.is_empty() {
            return Ok(load);
        }
        add_nodal_forces_to_thermal_load(load, self.nodes_m.len(), nodal_forces_n, cx)
    }

    /// Solve mechanical equilibrium with generic initial strain and dead load.
    ///
    /// This uses the same stiffness, force signs, checked linear solve and
    /// input binding as thermal loading, without creating a temperature.
    pub fn solve_eigenstrain_displacement_with_nodal_forces(
        &self,
        strains: TetEigenstrainField<'_>,
        nodal_forces_n: &[(usize, f64)],
        config: TetStaticSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<TetEigenstrainDisplacementSolution, TetElasticError> {
        let load = self.assemble_eigenstrain_load_with_nodal_forces(strains, nodal_forces_n, cx)?;
        self.solve_thermal_displacement_from_load(load, config, cx)
    }

    /// Recover stress and internal forces from the same initial-strain/load solve.
    ///
    /// Changed force bits, source strain, material, mesh or supports refuse
    /// through the existing physical input identity check. Actual reactions
    /// are internal minus external force at constrained DOFs.
    pub fn recover_eigenstrain_stress_with_nodal_forces(
        &self,
        solution: &TetEigenstrainDisplacementSolution,
        strains: TetEigenstrainField<'_>,
        nodal_forces_n: &[(usize, f64)],
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<TetEigenstrainStressReport, TetElasticError> {
        let load = self.assemble_eigenstrain_load_with_nodal_forces(strains, nodal_forces_n, cx)?;
        self.recover_thermal_stress_from_load(
            solution,
            EigenstrainFieldView::Generic(strains),
            load,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }
}

#[cfg(test)]
mod nodal_force_tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx(f: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 812,
                    kernel_id: 812,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx);
        });
    }

    #[test]
    fn g1_initial_strain_plus_dead_load_has_exact_affine_response() {
        with_cx(|cx| {
            // Numerical patch fixture: 1/6 m^3 tetrahedron, E=120 GPa,
            // nu=1/4, uniform initial strain 1e-4, sigma_xx=120 MPa.
            // The only nonzero free nodal force is sigma_xx*V=20 MN.
            let nodes = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            let tets = [[0, 1, 2, 3]];
            let material = TetElasticMaterial::try_new(2700.0, 120.0e9, 0.25, [21; 32]).unwrap();
            let initial = TetEigenstrainState::try_new(
                [1.0e-4, 1.0e-4, 1.0e-4, 0.0, 0.0, 0.0],
                [21; 32],
                ContentHash([22; 32]),
                0.01,
            )
            .unwrap();
            let field = TetEigenstrainField::Uniform(&initial);
            let problem = TetLinearElasticProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                materials: TetMaterialField::Uniform(&material),
                fixed_dofs: &[0, 1, 2, 4, 5, 8],
                budget: TetAssemblyBudget::standard(),
            };
            let forces = [(0, 7.0), (3, 20.0e6)];
            let solution = problem
                .solve_eigenstrain_displacement_with_nodal_forces(
                    field,
                    &forces,
                    TetStaticSolveConfig::default(),
                    cx,
                )
                .unwrap();
            let u = solution.displacement_m();
            assert!((u[1][0] - 0.0011).abs() < 2.0e-15);
            assert!((u[2][1] + 0.00015).abs() < 2.0e-15);
            assert!((u[3][2] + 0.00015).abs() < 2.0e-15);
            let report = problem
                .recover_eigenstrain_stress_with_nodal_forces(
                    &solution,
                    field,
                    &forces,
                    ContentHash([23; 32]),
                    0.01,
                    cx,
                )
                .unwrap();
            let stress = report.elements()[0].stress_mandel_pa();
            assert!((stress[0] / 120.0e6 - 1.0).abs() < 2.0e-12);
            assert!(stress[1..].iter().all(|x| x.abs() < 1.0e-5));
            // 1/2 * V * sigma_xx * elastic_epsilon_xx = 10 kJ.
            assert!((report.elastic_energy_j() / 10_000.0 - 1.0).abs() < 2.0e-12);
            let support_x = report.nodal_internal_force_n()[0][0] - forces[0].1;
            assert!((support_x + 20.0e6 + 7.0).abs() < 1.0e-5);

            let without_support_load = problem
                .solve_eigenstrain_displacement_with_nodal_forces(
                    field,
                    &[(3, 20.0e6)],
                    TetStaticSolveConfig::default(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                solution.displacement_m(),
                without_support_load.displacement_m()
            );
            assert_ne!(solution.identity(), without_support_load.identity());
            assert!(
                problem
                    .recover_eigenstrain_stress_with_nodal_forces(
                        &solution,
                        field,
                        &[(0, 7.0), (3, 19.0e6)],
                        ContentHash([23; 32]),
                        0.01,
                        cx,
                    )
                    .is_err()
            );
        });
    }

    #[test]
    fn g0_empty_forces_preserve_state_and_malformed_forces_refuse() {
        with_cx(|cx| {
            let nodes = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            let tets = [[0, 1, 2, 3]];
            let material = TetElasticMaterial::try_new(1000.0, 1.0e9, 0.25, [31; 32]).unwrap();
            let initial = TetEigenstrainState::try_new(
                [1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0],
                [31; 32],
                ContentHash([32; 32]),
                0.01,
            )
            .unwrap();
            let field = TetEigenstrainField::Uniform(&initial);
            let problem = TetLinearElasticProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                materials: TetMaterialField::Uniform(&material),
                fixed_dofs: &[0, 1, 2, 4, 5, 8],
                budget: TetAssemblyBudget::standard(),
            };
            let original = problem
                .solve_eigenstrain_displacement(field, TetStaticSolveConfig::default(), cx)
                .unwrap();
            let empty = problem
                .solve_eigenstrain_displacement_with_nodal_forces(
                    field,
                    &[],
                    TetStaticSolveConfig::default(),
                    cx,
                )
                .unwrap();
            assert_eq!(original, empty);
            let original_stress = problem
                .recover_eigenstrain_stress(&original, field, ContentHash([33; 32]), 0.01, cx)
                .unwrap();
            let empty_stress = problem
                .recover_eigenstrain_stress_with_nodal_forces(
                    &empty,
                    field,
                    &[],
                    ContentHash([33; 32]),
                    0.01,
                    cx,
                )
                .unwrap();
            assert_eq!(original_stress, empty_stress);
            for bad in [
                vec![(3, f64::NAN)],
                vec![(3, f64::INFINITY)],
                vec![(12, 1.0)],
                vec![(3, 1.0), (3, 2.0)],
                vec![(3, 1.0), (0, 1.0)],
            ] {
                assert!(matches!(
                    problem.assemble_eigenstrain_load_with_nodal_forces(field, &bad, cx),
                    Err(TetElasticError::InvalidNodalForce { .. })
                ));
            }
        });
    }
}
