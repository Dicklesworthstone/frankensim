//! Native elastic material fields in the existing density-design operator.
//!
//! Constitutive tensors and reference densities keep their original native
//! material-state identities. SIMP weights multiply the reference stiffness.
//! An explicitly authored nodal reference density retains its exact cubic
//! mass moments, normalized by each cell's mean density. Both modal front doors
//! and their derivatives use that same mean and mass tile.

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

pub(super) struct NodalReferenceDensity {
    values_kg_m3: Vec<f64>,
    means_kg_m3: Vec<f64>,
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
    /// Retain an original P1 nodal reference density alongside native stiffness.
    ///
    /// The supplied positive densities [kg/m³] follow the original mesh vertex
    /// order and are interpolated with its barycentric basis. They determine
    /// the closed reference inventory; material states still determine stiffness
    /// and keep their source identities. Neither a later design interpolation
    /// nor prescribed thermal strain requeries or changes this density field.
    /// This is an authored mass distribution, not a thermal expansion law.
    ///
    /// Integrating `rho N_a N_b` exactly retains center-of-mass and rotational
    /// inertia effects that replacing density by its centroid value would lose.
    /// Preparation is O(cells); all existing dense and matrix-free mass actions
    /// reuse the same 12×12 tiles without extra work during a solver iteration.
    /// Low-level mass weights multiply tiles normalized by
    /// [`Self::reference_density_kg_m3`], rather than the nominal material density.
    pub fn with_material_field_and_nodal_reference_density(
        complex: &TetComplex,
        positions: &[[f64; 3]],
        materials: TetMaterialField<'_>,
        reference_density_kg_m3: &[f64],
        fixed: &dyn Fn([f64; 3]) -> bool,
    ) -> Result<Self, MaterialElasticityError> {
        if reference_density_kg_m3.len() != complex.vertex_count
            || reference_density_kg_m3
                .iter()
                .any(|rho| !rho.is_finite() || *rho <= 0.0)
        {
            return Err(TetElasticError::InvalidMaterial {
                what: "one positive finite reference density per mesh vertex required",
            }
            .into());
        }
        let mut result = Self::with_material_field(complex, positions, materials, fixed)?;
        let mut means = Vec::with_capacity(result.cells());
        let mut masses = Vec::with_capacity(result.cells());
        for (cell, tet) in result.tets.iter().enumerate() {
            let rho = tet.map(|vertex| reference_density_kg_m3[vertex as usize]);
            let scale = rho.into_iter().fold(0.0_f64, f64::max);
            let scaled = rho.map(|value| value / scale);
            let sum: f64 = scaled.iter().sum();
            let uniform = rho.iter().all(|value| value.to_bits() == rho[0].to_bits());
            let mean = if uniform {
                rho[0]
            } else {
                scale * (0.25 * sum)
            };
            let mass = result.volumes[cell] * mean;
            if !mean.is_finite() || mean <= 0.0 || !mass.is_finite() || mass <= 0.0 {
                return Err(MaterialElasticityError::UnrepresentableElement { cell });
            }
            // Exact barycentric moments: diagonal V*(sum(rho)+2*rho_a)/60,
            // off-diagonal V*(sum(rho)+rho_a+rho_b)/120. Normalize before
            // summing so the mean does not overflow for large finite samples.
            // Uniform cells retain the original unit-density tile bit for bit.
            if !uniform {
                let fractions = scaled.map(|value| value / sum);
                for a in 0..4 {
                    for b in a..4 {
                        let value = if a == b {
                            (result.volumes[cell] / 15.0) * (1.0 + 2.0 * fractions[a])
                        } else {
                            (result.volumes[cell] / 30.0) * (1.0 + (fractions[a] + fractions[b]))
                        };
                        if !value.is_finite() || value <= 0.0 {
                            return Err(MaterialElasticityError::UnrepresentableElement { cell });
                        }
                        for component in 0..3 {
                            let row = 3 * a + component;
                            let column = 3 * b + component;
                            result.me[cell][12 * row + column] = value;
                            result.me[cell][12 * column + row] = value;
                        }
                    }
                }
            }
            // A representable cell total does not guarantee representable
            // consistent entries, including for a uniform subnormal density.
            for a in 0..4 {
                for b in a..4 {
                    let physical = result.me[cell][36 * a + 3 * b] * mean;
                    if !physical.is_finite() || physical <= 0.0 {
                        return Err(MaterialElasticityError::UnrepresentableElement { cell });
                    }
                }
            }
            means.push(mean);
            masses.push(mass);
        }
        result.nodal_reference_density = Some(NodalReferenceDensity {
            values_kg_m3: reference_density_kg_m3.to_vec(),
            means_kg_m3: means,
        });
        result.reference_cell_masses_kg = Some(masses);
        Ok(result)
    }

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
    /// density from reference_density_kg_m3 before their mass interpolation.
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

    /// Nominal density [kg/m³] of the retained material state, or None for an
    /// invalid cell or an unbound legacy model. An authored nodal inventory
    /// leaves this source datum unchanged; mass consumers use
    /// [`Self::reference_density_kg_m3`].
    #[must_use]
    pub fn material_density_kg_m3(&self, cell: usize) -> Option<f64> {
        if cell >= self.cells() {
            return None;
        }
        self.materials
            .as_ref()
            .map(|materials| materials.at(cell).density_kg_m3)
    }

    /// Mean physical reference density [kg/m³] for the prepared mass tile.
    /// An authored nodal field takes precedence over the nominal material
    /// density. Multiplying this mean by the normalized tile recovers the
    /// exact original density-weighted mass, including its spatial moments.
    /// Returns None for an invalid cell or a legacy unit-density model.
    #[must_use]
    pub fn reference_density_kg_m3(&self, cell: usize) -> Option<f64> {
        match &self.nodal_reference_density {
            Some(field) => field.means_kg_m3.get(cell).copied(),
            None => self.material_density_kg_m3(cell),
        }
    }

    /// Original immutable nodal reference inventory, if explicitly supplied.
    #[must_use]
    pub fn nodal_reference_density_kg_m3(&self) -> Option<&[f64]> {
        self.nodal_reference_density
            .as_ref()
            .map(|field| field.values_kg_m3.as_slice())
    }

    /// Full reference-cell masses [kg], before design interpolation and without
    /// eliminating mass attached to fixed DOFs. Absent for the scalar constructor.
    #[must_use]
    pub fn reference_cell_masses_kg(&self) -> Option<&[f64]> {
        self.reference_cell_masses_kg.as_deref()
    }
}
