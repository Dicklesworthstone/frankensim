//! Existing tetrahedral strain/quadrature maps projected into caller coordinates.
//! No mass lumping, new element kernel, material averaging or modal damping fit.

use super::*;
use fs_exec::CancelGate;
use fs_material::visco_tensor::{ProjectedTensorProny, TensorProjectionBudget};

/// Explicit total cold-projection capacities. These do not truncate a basis,
/// quadrature rule or constitutive memory when the request exceeds a limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SolidProjectionBudget {
    /// Displacement coordinates; in 1..=256, as admitted by the material owner.
    pub max_coordinates: usize,
    /// Total stored free-DOF basis entries.
    pub max_basis_scalars: usize,
    /// Total material quadrature points across all elements.
    pub max_material_points: usize,
    /// Total scalar multiply-add terms for mass, strain and material projections.
    pub max_scalar_terms: usize,
}

/// Immutable projection bound to the original geometry, mass, material laws and
/// fixed support map. Columns map displacement coordinates [m] into nodal
/// displacements [m]; they are dimensionless, not implicitly mass-normalized.
pub struct ProjectedViscoSolid<'a> {
    source: &'a ViscoElements,
    basis: Vec<Vec<f64>>,
    mass_kg: Vec<f64>,
    points: Vec<ProjectedTensorProny>,
    point_elements: Vec<usize>,
}
impl<'s> ProjectedViscoSolid<'s> {
    /// Original mesh and reference assembly. After a material rebuild, current
    /// constitutive operators are `material_points()`, not this old stiffness.
    pub fn source(&self) -> &ViscoElements { self.source }
    pub fn coordinates(&self) -> usize { self.basis.len() }
    /// Original columns in the source assembly's canonical free-DOF order.
    pub fn basis(&self) -> &[Vec<f64>] { &self.basis }
    /// Complete symmetric Phi^T M Phi, row-major [kg]; never lumped.
    pub fn mass_matrix_kg(&self) -> &[f64] { &self.mass_kg }
    /// Original tensor material projections, in element/quadrature order.
    pub fn material_points(&self) -> &[ProjectedTensorProny] { &self.points }
    /// Element index for each retained quadrature point.
    pub fn point_elements(&self) -> &[usize] { &self.point_elements }
    pub fn element_count(&self) -> usize { self.source.elements.len() }

    /// Rebuild material operators on the SAME reference geometry, support map,
    /// complete basis, positive quadrature and consistent mass. No basis change
    /// or internal-variable transfer is inferred. The coupling owner must bind
    /// the replacement laws to their original source and transfer live history.
    pub fn with_element_laws(
        &self, laws: &[TensorProny], budget: SolidProjectionBudget, gate: &CancelGate,
    ) -> Result<ProjectedViscoSolid<'s>, ViscoElementError> {
        poll_map(gate)?;
        let n = self.coordinates();
        if laws.len() != self.element_count() || n > budget.max_coordinates
            || !(1..=256).contains(&budget.max_coordinates)
            || self.points.len() > budget.max_material_points
            || checked_mul(n, self.source.equilibrium.free_dofs.len())? > budget.max_basis_scalars {
            return Err(invalid("material rebuild requires complete element laws and original projection capacity"));
        }
        let mut work = 0;
        for &element in &self.point_elements {
            work = checked_add(work, material_work(n, &laws[element])?)?;
        }
        if work > budget.max_scalar_terms {
            return Err(invalid("material rebuild exceeds projection work capacity"));
        }
        let mut points = Vec::with_capacity(self.points.len());
        for (old, &element) in self.points.iter().zip(&self.point_elements) {
            poll_map(gate)?;
            points.push(laws[element].project_strain(
                old.strain_columns(), old.volume_m3(),
                TensorProjectionBudget {
                    max_coordinates: budget.max_coordinates,
                    max_scalar_terms: material_work(n, &laws[element])?,
                },
                || gate.is_requested(),
            ).and_then(|point| point.with_eigenstrain_mandel(*old.free_strain_mandel()))
                .map_err(ViscoElementError::Projection)?);
        }
        poll_map(gate)?;
        Ok(ProjectedViscoSolid {
            source: self.source, basis: self.basis.clone(), mass_kg: self.mass_kg.clone(),
            points, point_elements: self.point_elements.clone(),
        })
    }

    /// Apply the original solid owner's admitted per-element thermal strains
    /// to the SAME quadrature, basis and reference mass. The source coupling
    /// owner must bind material/frame identities before this lowering.
    pub fn with_element_thermal_strains(
        &self, strains: &[crate::linear3::TetThermalStrainState], gate: &CancelGate,
    ) -> Result<Self, ViscoElementError> {
        poll_map(gate)?;
        if strains.len() != self.element_count() {
            return Err(invalid("one admitted thermal strain per original element required"));
        }
        let mut points = Vec::with_capacity(self.points.len());
        for (point, &element) in self.points.iter().zip(&self.point_elements) {
            poll_map(gate)?;
            points.push(point.with_eigenstrain_mandel(strains[element].free_strain_mandel())
                .map_err(ViscoElementError::Projection)?);
        }
        poll_map(gate)?;
        Ok(Self { source: self.source, basis: self.basis.clone(), mass_kg: self.mass_kg.clone(),
            points, point_elements: self.point_elements.clone() })
    }

    /// Apply the SAME basis transpose to a full node-major force field.
    /// Fixed rows do no virtual work. The returned generalized force is
    /// conjugate to the caller's displacement coordinates and their velocities.
    pub fn project_nodal_force(
        &self, force_n: &[[f64; 3]], gate: &CancelGate,
    ) -> Result<Vec<f64>, ViscoElementError> {
        poll_map(gate)?;
        if force_n.len().checked_mul(3) != Some(self.source.full_dofs)
            || force_n.iter().flatten().any(|x| !x.is_finite())
        { return Err(invalid("one finite force vector per original node required")); }
        let mut force = Vec::with_capacity(self.coordinates());
        for column in &self.basis {
            poll_map(gate)?;
            let value: f64 = column.iter().zip(&self.source.equilibrium.free_dofs)
                .map(|(a, &d)| a*force_n[d/3][d%3]).sum();
            if !value.is_finite() { return Err(ViscoElementError::NonFiniteResult); }
            force.push(value);
        }
        Ok(force)
    }

    /// Expand displacement or velocity through the same basis and support map.
    /// Fixed components remain zero; no second mode or observer fit is used.
    pub fn expand_nodal(
        &self, coordinates: &[f64], gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, ViscoElementError> {
        poll_map(gate)?;
        if coordinates.len() != self.coordinates() || coordinates.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite coordinates matching the retained solid basis required"));
        }
        let mut nodal = vec![[0.0; 3]; self.source.full_dofs/3];
        for (r, &d) in self.source.equilibrium.free_dofs.iter().enumerate() {
            poll_map(gate)?;
            let value: f64 = self.basis.iter().zip(coordinates).map(|(column, q)| column[r]*q).sum();
            if !value.is_finite() { return Err(ViscoElementError::NonFiniteResult); }
            nodal[d/3][d%3] = value;
        }
        Ok(nodal)
    }
}

fn poll_map(gate: &CancelGate) -> Result<(), ViscoElementError> {
    if gate.is_requested() { Err(TetElasticError::Cancelled.into()) } else { Ok(()) }
}
fn invalid(what: &'static str) -> ViscoElementError {
    ViscoElementError::InvalidHistory { what }
}
fn checked_add(a: usize, b: usize) -> Result<usize, ViscoElementError> {
    a.checked_add(b).ok_or_else(|| invalid("solid projection capacity arithmetic overflow"))
}
fn checked_mul(a: usize, b: usize) -> Result<usize, ViscoElementError> {
    a.checked_mul(b).ok_or_else(|| invalid("solid projection capacity arithmetic overflow"))
}
fn material_work(n: usize, law: &TensorProny) -> Result<usize, ViscoElementError> {
    let tensors = 1 + usize::from(law.kv_pa_s().is_some());
    checked_add(
        checked_mul(checked_mul(checked_add(law.branches().len(), tensors)?, 36)?, n)?,
        checked_mul(checked_mul(6, tensors)?, checked_mul(n, n)?)?,
    )
}

impl ViscoElements {
    /// Project complete 3D source laws and CONSISTENT inertia on a supplied
    /// finite basis. Full branch factors and Kelvin cross terms survive.
    ///
    /// The source's original B matrices, physical quadrature volumes and free
    /// DOF map are reused for both Tet4 and Tet10. A dynamic consumer must admit
    /// positive-definite projected inertia; dependent columns are not repaired.
    /// Spatial/modal convergence and source frequency applicability remain
    /// explicit obligations, as with the original time-domain history owner.
    pub fn project_basis(
        &self, columns: &[Vec<f64>], budget: SolidProjectionBudget, cx: &Cx<'_>,
    ) -> Result<ProjectedViscoSolid<'_>, ViscoElementError> {
        checkpoint(cx)?;
        let n = columns.len();
        let free = self.equilibrium.free_dofs.len();
        if !(1..=256).contains(&budget.max_coordinates) || n == 0 || n > budget.max_coordinates
            || budget.max_material_points == 0
            || checked_mul(n, free)? > budget.max_basis_scalars
            || columns.iter().any(|v| v.len() != free || v.iter().any(|x| !x.is_finite()))
        { return Err(invalid("bounded complete finite solid basis required")); }
        let mut count = 0;
        let mut work = checked_add(
            checked_mul(n, self.equilibrium.mass.nnz())?,
            checked_mul(checked_mul(n, n)?, free)?,
        )?;
        for element in &self.elements {
            checkpoint(cx)?;
            count = checked_add(count, element.points.len())?;
            for point in &element.points {
                work = checked_add(work, checked_add(
                    checked_mul(checked_mul(point.b.len(), 6)?, n)?,
                    material_work(n, &element.law)?,
                )?)?;
            }
        }
        if count > budget.max_material_points || work > budget.max_scalar_terms {
            return Err(invalid("complete solid projection exceeds point/work capacity"));
        }

        let mut mass = vec![0.0; n*n];
        let mut product = vec![0.0; free];
        for j in 0..n {
            checkpoint(cx)?;
            self.equilibrium.mass.spmv(&columns[j], &mut product);
            for i in 0..=j {
                let value: f64 = columns[i].iter().zip(&product).map(|(a,b)| a*b).sum();
                if !value.is_finite() { return Err(ViscoElementError::NonFiniteResult); }
                mass[i*n+j] = value; mass[j*n+i] = value;
            }
        }
        let mut points = Vec::with_capacity(count);
        let mut point_elements = Vec::with_capacity(count);
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let reduced: Vec<_> = element.dofs.iter()
                .map(|d| self.equilibrium.free_dofs.binary_search(d).ok()).collect();
            for point in &element.points {
                checkpoint(cx)?;
                let mut strains = Vec::with_capacity(n);
                for column in columns {
                    let mandel: [f64; 6] = core::array::from_fn(|i|
                        point.b.iter().zip(&reduced).map(|(b, r)|
                            r.map_or(0.0, |r| b[i]*column[r])).sum());
                    // The element uses Mandel shear; project_strain takes
                    // tensor shear. This conversion retains stress:strain work.
                    strains.push(core::array::from_fn(|i| if i < 3 { mandel[i] }
                        else { mandel[i]*core::f64::consts::FRAC_1_SQRT_2 }));
                }
                let projected = element.law.project_strain(
                    &strains, point.volume_m3,
                    TensorProjectionBudget {
                        max_coordinates: budget.max_coordinates,
                        max_scalar_terms: material_work(n, &element.law)?,
                    },
                    || cx.checkpoint().is_err(),
                ).map_err(ViscoElementError::Projection)?;
                points.push(projected);
                point_elements.push(e);
            }
        }
        checkpoint(cx)?;
        Ok(ProjectedViscoSolid {
            source: self, basis: columns.to_vec(), mass_kg: mass, points, point_elements,
        })
    }
}
