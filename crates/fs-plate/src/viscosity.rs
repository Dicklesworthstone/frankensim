//! Physical bending viscosity on the original spatial DKT coordinates.
//!
//! Each element integrates `C_e = integral B^T V_b B dA` with the same
//! curvature quadrature as its elastic owner. Geometric prestress, stiffeners,
//! rigid motion and numerical stabilization acquire no material loss here.

use crate::{PlateError, PlateMesh, PlateModel, dkt_stiffness};
use fs_material::plane_stress_viscosity::{PlaneStressViscosity, factor_passive_resultant};
use fs_sparse::{Coo, Csr};
use std::sync::Arc;

fn bad(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// One element's physical bending viscosity in mesh coordinates [N m s].
#[derive(Clone, Debug)]
pub struct ViscousPlateSection {
    bending: [f64; 9],
    source: Option<Arc<PlaneStressViscosity>>,
    resolved_band_rad_s: Option<(f64, f64)>,
}

impl ViscousPlateSection {
    /// Admit an explicitly supplied passive resultant. Directional loss and
    /// exact zero channels are legal; no positive loss floor is introduced.
    pub fn from_resultant(bending: [f64; 9]) -> Result<Self, PlateError> {
        factor_passive_resultant(&bending).map_err(|_| {
            bad("plate bending viscosity must be finite symmetric positive semidefinite")
        })?;
        Ok(Self {
            bending,
            source: None,
            resolved_band_rad_s: None,
        })
    }

    /// Reduce a retained material law at the supplied thickness/orientation.
    /// `frequency_band` is the caller's declared resolved transient band [rad/s];
    /// this admission does not establish that a trajectory stays in that band.
    pub fn from_material(
        source: PlaneStressViscosity,
        thickness_m: f64,
        angle_rad: f64,
        frequency_band: (f64, f64),
    ) -> Result<Self, PlateError> {
        if !frequency_band.0.is_finite()
            || !frequency_band.1.is_finite()
            || frequency_band.0 < 0.0
            || frequency_band.1 < frequency_band.0
        {
            return Err(bad("plate viscosity needs a finite ordered resolved frequency band"));
        }
        for omega in [frequency_band.0, frequency_band.1] {
            source.validate_frequency(omega).map_err(|_| {
                bad("plate viscosity resolved band exceeds its material source band")
            })?;
        }
        let bending = source
            .bending_matrix_n_m_s(thickness_m, angle_rad)
            .map_err(|_| bad("plate material viscosity section reduction failed"))?;
        let mut section = Self::from_resultant(bending)?;
        section.source = Some(Arc::new(source));
        section.resolved_band_rad_s = Some(frequency_band);
        Ok(section)
    }

    /// Original selected material law, when supplied through `from_material`.
    pub fn source(&self) -> Option<&PlaneStressViscosity> {
        self.source.as_deref()
    }

    /// Declared resolved transient band [rad/s], for a source-derived section.
    pub const fn resolved_band_rad_s(&self) -> Option<(f64, f64)> {
        self.resolved_band_rad_s
    }

    /// Physical mesh-frame bending tensor [N m s].
    pub const fn bending_n_m_s(&self) -> &[f64; 9] {
        &self.bending
    }
}

/// Cold assembly ceilings, checked before matrix or support-map allocation.
#[derive(Clone, Copy, Debug)]
pub struct PlateViscosityBudget {
    /// Maximum triangle count.
    pub max_elements: usize,
    /// Maximum number of free spatial coordinates.
    pub max_dofs: usize,
    /// Maximum number of coordinates before support elimination.
    pub max_full_dofs: usize,
    /// Upper bound of 81 element entries per triangle before support elimination.
    pub max_assembly_entries: usize,
}

/// A passive material-loss operator tied to one original supported plate model.
/// Its sparse matrix and all original selected section laws are immutable.
#[derive(Debug)]
pub struct PlateBendingViscosity<'a> {
    model: &'a PlateModel,
    matrix: Csr,
    sections: Vec<ViscousPlateSection>,
}

impl<'a> PlateBendingViscosity<'a> {
    /// Assemble complete element coverage on the supplied original support map.
    /// No modal projection, Rayleigh fit, new inertia, or altered stiffness is used.
    /// The caller supplies the same mesh geometry/connectivity used to assemble
    /// `model`; the model does not itself retain a mesh for comparison.
    pub fn assemble<C: FnMut() -> bool>(
        mesh: &PlateMesh,
        model: &'a PlateModel,
        sections: Vec<ViscousPlateSection>,
        budget: PlateViscosityBudget,
        cancelled: &mut C,
    ) -> Result<Self, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let full = mesh
            .nodes
            .len()
            .checked_mul(3)
            .ok_or_else(|| bad("plate viscosity DOF overflow"))?;
        let entries = mesh
            .tris
            .len()
            .checked_mul(81)
            .ok_or_else(|| bad("plate viscosity assembly overflow"))?;
        for (what, requested, limit) in [
            ("plate viscosity elements", mesh.tris.len(), budget.max_elements),
            ("plate viscosity DOFs", model.free, budget.max_dofs),
            ("plate viscosity full DOFs", full, budget.max_full_dofs),
            ("plate viscosity assembly entries", entries, budget.max_assembly_entries),
        ] {
            if requested > limit {
                return Err(PlateError::Budget { what, requested, limit });
            }
        }
        if sections.len() != mesh.tris.len() {
            return Err(PlateError::SectionCount {
                expected: mesh.tris.len(),
                actual: sections.len(),
            });
        }
        if model.dof_map.len() != full || model.free > full {
            return Err(bad("plate viscosity requires the original mesh support map"));
        }
        let mut seen = vec![false; model.free];
        for &index in model.dof_map.iter().flatten() {
            if index >= model.free || seen[index] {
                return Err(bad("invalid plate viscosity support map"));
            }
            seen[index] = true;
        }
        if seen.iter().any(|present| !present) {
            return Err(bad("incomplete plate viscosity support map"));
        }
        mesh.validate()?;
        let mut coo = Coo::new(model.free, model.free);
        for (element, (tri, section)) in mesh.tris.iter().zip(&sections).enumerate() {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let x = tri.map(|node| mesh.nodes[node].0);
            let y = tri.map(|node| mesh.nodes[node].1);
            let (local, _) = dkt_stiffness(&x, &y, &section.bending, element)?;
            if local.iter().any(|entry| !entry.is_finite()) {
                return Err(bad("plate viscosity element overflow"));
            }
            for row in 0..9 {
                for column in 0..9 {
                    let value = local[9 * row + column];
                    if let (Some(i), Some(j)) = (
                        model.dof_map[3 * tri[row / 3] + row % 3],
                        model.dof_map[3 * tri[column / 3] + column % 3],
                    ) {
                        if value != 0.0 {
                            coo.push(i, j, value);
                        }
                    }
                }
            }
        }
        let matrix = coo.assemble();
        if (0..model.free).any(|row| matrix.row(row).1.iter().any(|x| !x.is_finite())) {
            return Err(bad("plate viscosity assembly overflow"));
        }
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        Ok(Self { model, matrix, sections })
    }

    /// Exact original model identity; a different support/pencil needs reassembly.
    pub const fn model(&self) -> &'a PlateModel {
        self.model
    }

    /// Full sparse spatial damping matrix in the original reduced coordinates.
    pub const fn matrix(&self) -> &Csr {
        &self.matrix
    }

    /// Original per-element material/resultant selections in mesh order.
    pub fn sections(&self) -> &[ViscousPlateSection] {
        &self.sections
    }
}
