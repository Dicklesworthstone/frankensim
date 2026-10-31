//! Native elastic material fields in the existing density-design operator.
//!
//! Constitutive tensors and reference densities keep their original native
//! material-state identities. SIMP weights multiply the reference stiffness;
//! the consistent mass blocks remain unit-density so existing mass interpolation
//! and its derivative can use the same operator without hidden double weighting.

use super::DensityElasticity;
use fs_feec::whitney::{ElementGeometryError, try_element_geometry};
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{TetElasticError, TetElasticMaterial, TetMaterialField};

/// Admission failures from the existing geometry/material owners, or a finite
/// input whose element operator or reference mass cannot be represented.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialElasticityError {
    /// Invalid coordinates, connectivity or finite-element metrics.
    Geometry(ElementGeometryError),
    /// Invalid native material, field cardinality or mesh inventory.
    Solid(TetElasticError),
    /// Nonfinite stiffness/stress/mass, or a positive physical quantity that
    /// underflowed to zero in this cell's stored operators.
    UnrepresentableElement { cell: usize },
}
impl core::fmt::Display for MaterialElasticityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry(error) => write!(f, "material elasticity geometry: {error:?}"),
            Self::Solid(error) => write!(f, "{error}"),
            Self::UnrepresentableElement { cell } => {
                write!(f, "material elasticity element {cell} is not representable")
            }
        }
    }
}
impl std::error::Error for MaterialElasticityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            _ => None,
        }
    }
}
impl From<ElementGeometryError> for MaterialElasticityError {
    fn from(error: ElementGeometryError) -> Self {
        Self::Geometry(error)
    }
}
impl From<TetElasticError> for MaterialElasticityError {
    fn from(error: TetElasticError) -> Self {
        Self::Solid(error)
    }
}

/// One retained uniform state or an ordered field; never one clone per cell
/// for the uniform case.
pub(super) enum OwnedMaterialField {
    Uniform(TetElasticMaterial),
    PerElement(Vec<TetElasticMaterial>),
}
impl OwnedMaterialField {
    fn new(field: TetMaterialField<'_>, cells: usize) -> Result<Self, TetElasticError> {
        // Public density and identity fields can be edited after initial native
        // construction. Re-admit the complete state at this consumer boundary.
        let admit = |material: &TetElasticMaterial| {
            TetElasticMaterial::try_new_mandel(
                material.density_kg_m3,
                *material.stiffness_mandel_pa(),
                material.material_state_identity,
            )
        };
        match field {
            TetMaterialField::Uniform(material) => Ok(Self::Uniform(admit(material)?)),
            TetMaterialField::PerElement(materials) => {
                if materials.len() != cells {
                    return Err(TetElasticError::MaterialCountMismatch {
                        expected: cells,
                        got: materials.len(),
                    });
                }
                Ok(Self::PerElement(
                    materials.iter().map(admit).collect::<Result<_, _>>()?,
                ))
            }
        }
    }

    fn at(&self, cell: usize) -> &TetElasticMaterial {
        match self {
            Self::Uniform(material) => material,
            Self::PerElement(materials) => &materials[cell],
        }
    }

    fn borrowed(&self) -> TetMaterialField<'_> {
        match self {
            Self::Uniform(material) => TetMaterialField::Uniform(material),
            Self::PerElement(materials) => TetMaterialField::PerElement(materials),
        }
    }
}

impl DensityElasticity {
    /// Build the existing design operator from complete native material states.
    ///
    /// The material field is copied once and remains immutable. Its full global
    /// Mandel tensors, densities and state identities survive; no isotropic fit,
    /// scalar modulus extraction or material-name lookup occurs. Per-element
    /// states follow the complex's tetrahedron order.
    ///
    /// The existing engineering-shear B matrix needs the work-conjugate change
    /// of basis C_eng = D C_M D, D = diag(1,1,1,1/sqrt(2),1/sqrt(2),1/sqrt(2)).
    /// Both rows and columns are scaled, preserving shear stress and energy.
    /// Initial stiffness weights are one; later weights are dimensionless.
    /// Consistent mass blocks remain unit-density. Consumers obtain physical
    /// density from material_density_kg_m3 before their mass interpolation.
    ///
    /// # Errors
    /// Empty, malformed or unrepresentable meshes; incomplete or invalid
    /// material fields; nonfinite/underflowed element operators or cell masses.
    pub fn with_material_field(
        complex: &TetComplex,
        positions: &[[f64; 3]],
        materials: TetMaterialField<'_>,
        fixed: &dyn Fn([f64; 3]) -> bool,
    ) -> Result<Self, MaterialElasticityError> {
        if complex.vertex_count == 0 || complex.tets.is_empty() {
            return Err(TetElasticError::EmptyMesh.into());
        }
        complex
            .vertex_count
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        complex
            .tets
            .len()
            .checked_mul(144)
            .ok_or(TetElasticError::SizeOverflow)?;
        let materials = OwnedMaterialField::new(materials, complex.tets.len())?;
        let geometry = try_element_geometry(complex, positions)?;
        let mut referenced = vec![false; complex.vertex_count];
        for tet in &complex.tets {
            for &vertex in tet {
                referenced[vertex as usize] = true;
            }
        }
        if let Some(node) = referenced.iter().position(|used| !used) {
            return Err(TetElasticError::UnreferencedNode { node }.into());
        }
        let mut masses = Vec::with_capacity(complex.tets.len());
        for (cell, volume) in geometry.vol_signed.iter().enumerate() {
            let mass = volume.abs() * materials.at(cell).density_kg_m3;
            if !mass.is_finite() || mass <= 0.0 {
                return Err(MaterialElasticityError::UnrepresentableElement { cell });
            }
            masses.push(mass);
        }
        let mut result = Self::from_engineering_field(
            complex,
            positions,
            &geometry,
            |cell| {
                let mandel = materials.at(cell).stiffness_mandel_pa();
                let factors = [
                    1.0,
                    1.0,
                    1.0,
                    core::f64::consts::FRAC_1_SQRT_2,
                    core::f64::consts::FRAC_1_SQRT_2,
                    core::f64::consts::FRAC_1_SQRT_2,
                ];
                core::array::from_fn(|row| {
                    core::array::from_fn(|column| {
                        mandel[row][column] * factors[row] * factors[column]
                    })
                })
            },
            fixed,
        );
        for cell in 0..result.cells() {
            if result.ke[cell]
                .iter()
                .chain(&result.cb[cell])
                .chain(&result.me[cell])
                .any(|value| !value.is_finite())
                || (0..12).any(|dof| result.ke[cell][13 * dof] <= 0.0)
                || result.me[cell][0] <= 0.0
            {
                return Err(MaterialElasticityError::UnrepresentableElement { cell });
            }
        }
        result.materials = Some(materials);
        result.reference_cell_masses_kg = Some(masses);
        Ok(result)
    }

    /// Original native field in connectivity order; absent for the legacy
    /// scalar constructor. Borrowing cannot mutate the retained source states.
    #[must_use]
    pub fn material_field(&self) -> Option<TetMaterialField<'_>> {
        self.materials.as_ref().map(OwnedMaterialField::borrowed)
    }

    /// Physical reference density [kg/m³], or None for an invalid cell index
    /// or a legacy model with no explicit density/material-state binding.
    #[must_use]
    pub fn material_density_kg_m3(&self, cell: usize) -> Option<f64> {
        if cell >= self.cells() {
            return None;
        }
        self.materials
            .as_ref()
            .map(|materials| materials.at(cell).density_kg_m3)
    }

    /// Full reference-cell masses [kg], before design interpolation and without
    /// eliminating mass attached to fixed DOFs. Absent for the scalar constructor.
    #[must_use]
    pub fn reference_cell_masses_kg(&self) -> Option<&[f64]> {
        self.reference_cell_masses_kg.as_deref()
    }
}
