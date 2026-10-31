//! Steady anisotropic Stokes–Brinkman flow on the original BDM1–P0 mesh.
//!
//! The superficial velocity solves
//! `-div(2 nu_eff D(u)) + nu K^-1 u + grad(pi) = f`, `div(u) = 0`.
//! Intrinsic `K` is in m², `nu` is the actual fluid's molecular viscosity
//! divided by the constant reference density, and `pi = p/rho`. The effective
//! Brinkman shear viscosity is a separately declared closure, never inferred
//! from porosity or substituted for the molecular viscosity in Darcy drag.
//!
//! Reference: Sandia Aria, *Brinkman Momentum*,
//! <https://www.sandia.gov/files/sierra/Aria_Users_5_28/theory/gov_eq/brinkman_mom.html>.
//! The rigid, stationary skeleton is fully saturated by one incompressible
//! fluid. There is no pore inertia/storage, Forchheimer/slip correction,
//! moving skeleton, unsaturated mobility, thermal deposition or Biot wave.
//! This steady entry inherits the existing synchronous Stokes solver's
//! iteration budget and cancellation limitations; it adds no cancellation claim.

use super::boundary::{StokesAssembly, StokesBoundary, StokesBoundaryError, StokesSolution};
use super::{FluxParams, FluxSystem};
use crate::bdm::{CellBasis, eval_basis, tri_quad};
use fs_material::porous::permeability::IntrinsicPermeability;
use fs_sparse::Csr;

/// One cell's fixed geometry and explicitly selected Brinkman shear closure.
/// The original permeability tensor and its world axes remain intact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrinkmanCell {
    permeability: Option<IntrinsicPermeability>,
    effective_viscosity_ratio: f64,
}

impl BrinkmanCell {
    /// A resolved free-fluid cell: molecular stress and exactly zero pore drag.
    #[must_use]
    pub const fn free_fluid() -> Self {
        Self {
            permeability: None,
            effective_viscosity_ratio: 1.0,
        }
    }

    /// Admit a porous cell with the stated `mu_eff / mu` closure.
    ///
    /// This two-dimensional owner requires the xy plane to be an invariant
    /// subspace of the permeability tensor. Rotations within that plane are
    /// admitted; coupling to an unresolved z velocity is refused. A sealed
    /// principal direction needs an explicit wall/constraint formulation and
    /// is refused here, without replacing zero permeability by an epsilon.
    ///
    /// # Errors
    /// Refuses a nonpositive/nonfinite viscosity ratio, sealed axes or
    /// out-of-plane permeability coupling.
    pub fn try_new(
        permeability: IntrinsicPermeability,
        effective_viscosity_ratio: f64,
    ) -> Result<Self, StokesBoundaryError> {
        if !positive(effective_viscosity_ratio)
            || permeability.principal_m2().iter().any(|&k| !positive(k))
        {
            return Err(invalid(
                "Brinkman requires positive permeability and an explicit positive effective-viscosity ratio",
            ));
        }
        // Each omitted xz/yz entry must cancel relative to the terms that
        // actually form that entry. A large resistance on an unrelated axis
        // cannot hide real coupling by setting a global normalization scale.
        for a in 0..2 {
            let mut products = [0.0; 3];
            let mut smallest = f64::INFINITY;
            for (i, axis) in permeability.material_axes_in_world().iter().enumerate() {
                products[i] = axis[a] * axis[2];
                if axis[a] != 0.0 && axis[2] != 0.0 {
                    if products[i] == 0.0 {
                        return Err(invalid("out-of-plane tensor coupling is not representable"));
                    }
                    smallest = smallest.min(permeability.principal_m2()[i]);
                }
            }
            if smallest.is_infinite() {
                continue;
            }
            let mut cross = 0.0;
            let mut magnitude = 0.0;
            for (i, product) in products.into_iter().enumerate() {
                if product == 0.0 {
                    continue;
                }
                let term = (smallest / permeability.principal_m2()[i]) * product;
                cross += term;
                magnitude += term.abs();
            }
            if cross.abs() > 4.0e-12 * magnitude {
                return Err(invalid(
                    "two-dimensional Brinkman refuses out-of-plane permeability coupling",
                ));
            }
        }
        Ok(Self {
            permeability: Some(permeability),
            effective_viscosity_ratio,
        })
    }

    /// Original intrinsic tensor; `None` identifies a resolved free-fluid cell.
    #[must_use]
    pub const fn permeability(&self) -> Option<&IntrinsicPermeability> {
        self.permeability.as_ref()
    }

    /// Declared effective shear viscosity divided by molecular fluid viscosity.
    #[must_use]
    pub const fn effective_viscosity_ratio(&self) -> f64 {
        self.effective_viscosity_ratio
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(reason: &'static str) -> StokesBoundaryError {
    StokesBoundaryError::InvalidInput(reason)
}

/// Spectral factors preserve nonnegative dissipation even for rotated tensors.
/// Conversion from geometry and fluid viscosity occurs once per cell.
#[derive(Clone, Copy)]
pub(super) struct DarcyDrag {
    rates: [f64; 3],
    axes: [[f64; 2]; 3],
}

impl DarcyDrag {
    fn compile(cell: BrinkmanCell, nu: f64) -> Result<Self, StokesBoundaryError> {
        let Some(permeability) = cell.permeability else {
            return Ok(Self {
                rates: [0.0; 3],
                axes: [[0.0; 2]; 3],
            });
        };
        let rates = permeability.principal_m2().map(|k| nu / k);
        if rates.iter().any(|&value| !positive(value)) {
            return Err(invalid(
                "molecular viscosity divided by intrinsic permeability is not representable",
            ));
        }
        let axes = permeability
            .material_axes_in_world()
            .map(|axis| [axis[0], axis[1]]);
        Ok(Self { rates, axes })
    }

    pub(super) fn local_matrix(
        &self,
        basis: &CellBasis,
        vertices: [[f64; 2]; 3],
        area: f64,
    ) -> [[f64; 6]; 6] {
        let mut matrix = [[0.0; 6]; 6];
        if self.rates == [0.0; 3] {
            return matrix;
        }
        for (point, weight) in tri_quad(vertices, area) {
            let projected: [[f64; 3]; 6] = core::array::from_fn(|i| {
                let value = eval_basis(basis, i, point);
                self.axes
                    .map(|axis| axis[0] * value[0] + axis[1] * value[1])
            });
            for i in 0..6 {
                for j in i..6 {
                    let value: f64 = (0..3)
                        .map(|a| self.rates[a] * projected[i][a] * projected[j][a])
                        .sum();
                    matrix[i][j] += weight * value;
                }
            }
        }
        for i in 0..6 {
            for j in 0..i {
                matrix[i][j] = matrix[j][i];
            }
        }
        matrix
    }

    fn power(&self, velocity: [f64; 2]) -> f64 {
        (0..3)
            .map(|a| {
                let projection = self.axes[a][0] * velocity[0] + self.axes[a][1] * velocity[1];
                self.rates[a] * projection * projection
            })
            .sum()
    }
}

/// The native Stokes operator plus the exact compiled pore-drag coefficients.
/// Assembly and checked solves reuse its boundary, continuity and pressure owner.
pub struct BrinkmanAssembly<'a, 'm> {
    system: &'a FluxSystem<'m>,
    stokes: StokesAssembly<'a, 'm>,
    drag: Vec<DarcyDrag>,
}

impl BrinkmanAssembly<'_, '_> {
    /// Full saddle matrix; the existing discrete transpose solver can consume it.
    #[must_use]
    pub fn matrix(&self) -> &Csr {
        self.stokes.matrix()
    }
    /// Native volume, velocity and total-traction load.
    #[must_use]
    pub fn rhs(&self) -> &[f64] {
        self.stokes.rhs()
    }

    /// Solve with the original residual checks and declared pressure datum.
    ///
    /// # Errors
    /// Preserves the native Stokes factorization and convergence refusals.
    pub fn solve(&self) -> Result<StokesSolution, StokesBoundaryError> {
        self.stokes.solve()
    }

    /// Integrate `nu u.K^-1.u` in each triangle using the assembled coefficients.
    /// Values have units m⁴/s³ for this 2-D mesh; multiply by the SAME reference
    /// density to obtain physical Darcy dissipation per out-of-plane depth [W/m].
    /// These losses are diagnostics, never automatically added to a heat model.
    ///
    /// # Errors
    /// Refuses incorrect velocity dimensions, nonfinite states or overflow.
    pub fn cell_darcy_power_per_density(
        &self,
        velocity: &[f64],
    ) -> Result<Vec<f64>, StokesBoundaryError> {
        if velocity.len() != self.system.n_u || velocity.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "Darcy dissipation requires exactly the finite native velocity DOFs",
            ));
        }
        let mesh = self.system.mesh;
        let mut powers = Vec::with_capacity(mesh.tris.len());
        for (cell, drag) in self.drag.iter().enumerate() {
            let vertices = core::array::from_fn(|i| mesh.verts[mesh.tris[cell][i]]);
            let mut power = 0.0;
            for (point, weight) in tri_quad(vertices, mesh.areas[cell]) {
                let mut value = [0.0; 2];
                for i in 0..6 {
                    let basis = eval_basis(&self.system.bases[cell], i, point);
                    let dof = velocity[self.system.dof(cell, i)];
                    value[0] += dof * basis[0];
                    value[1] += dof * basis[1];
                }
                power += weight * drag.power(value);
            }
            if !power.is_finite() || power < 0.0 {
                return Err(invalid("Darcy dissipation is not representable"));
            }
            powers.push(power);
        }
        Ok(powers)
    }
}

impl<'m> FluxSystem<'m> {
    /// Assemble stationary, fully saturated Stokes–Brinkman momentum.
    ///
    /// `params.nu` is molecular fluid viscosity/density. Each cell provides an
    /// intrinsic tensor and separate effective shear-viscosity ratio; free-fluid
    /// cells have no drag. No permeability component is replaced by a diagonal
    /// approximation: all in-plane couplings enter the velocity block together.
    /// Cellwise coefficients are compiled once outside all quadrature loops.
    /// Strong normal moments, total-traction boundaries, SIP interfaces, exact
    /// cell continuity and both true residual checks remain owned by Stokes.
    /// A positive porous cell removes all rigid velocity modes of this single
    /// connected domain, admitting all-traction boundaries without an artificial
    /// velocity anchor. Disconnected meshes still refuse. Without pore drag the
    /// original requirement for a velocity boundary remains in force.
    ///
    /// # Errors
    /// Refuses the native Stokes admission failures, mismatched cell counts,
    /// nonpositive molecular viscosity and unrepresentable material coefficients.
    pub fn assemble_brinkman_mixed<'a>(
        &'a self,
        params: FluxParams,
        cells: &[BrinkmanCell],
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
    ) -> Result<BrinkmanAssembly<'a, 'm>, StokesBoundaryError> {
        if cells.len() != self.mesh.tris.len() || !positive(params.nu) {
            return Err(invalid(
                "Brinkman requires one cell law per triangle and positive molecular viscosity",
            ));
        }
        let mut drag = Vec::with_capacity(cells.len());
        let mut effective_nu = Vec::with_capacity(cells.len());
        for &cell in cells {
            drag.push(DarcyDrag::compile(cell, params.nu)?);
            let effective = params.nu * cell.effective_viscosity_ratio;
            if !positive(effective) {
                return Err(invalid("effective Brinkman viscosity is not representable"));
            }
            effective_nu.push(effective);
        }
        let admitted = self.admit_stokes_problem_with_drag(
            params,
            boundaries,
            cells.iter().any(|cell| cell.permeability.is_some()),
        )?;
        let stokes = self.assemble_stokes_admitted_viscosity_and_drag(
            params,
            &|cell, _| effective_nu[cell],
            force,
            &admitted,
            Some(&drag),
        )?;
        Ok(BrinkmanAssembly {
            system: self,
            stokes,
            drag,
        })
    }
}
