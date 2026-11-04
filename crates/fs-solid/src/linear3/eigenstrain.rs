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
}
