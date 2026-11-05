//! Anisotropic Stokes–Brinkman flow on the original BDM1–P0 mesh.
//!
//! The steady macroscopic velocity solves
//! `-div(2 nu_eff D(u)) + nu K^-1 (u-u_s) + grad(pi) = f`, `div(u) = 0`.
//! Intrinsic `K` is in m², `nu` is the actual fluid's molecular viscosity
//! divided by the constant reference density, and `pi = p/rho`. The effective
//! Brinkman shear viscosity is a separately declared closure, never inferred
//! from porosity or substituted for the molecular viscosity in Darcy drag.
//!
//! Reference: Sandia Aria, *Brinkman Momentum*,
//! <https://www.sandia.gov/files/sierra/Aria_Users_5_28/theory/gov_eq/brinkman_mom.html>.
//! Prescribed skeleton motion uses the one-domain Brinkman approximation,
//! with the phi=1 velocity convention of Sherwood (2020), equations 2.4–2.5,
//! <https://doi.org/10.1017/jfm.2020.265>. It does not identify finite-porosity
//! laboratory superficial velocity with intrinsic pore velocity: that would
//! require the phi factor in equation 2.3 and a different storage/flux model.
//! Skeleton geometry and coefficients are frozen during each solve. The BDF1
//! entry adds native fluid inertia, with full tensor drag implicit at the endpoint.
//! There is no finite-porosity storage or pore-scale added-mass correction,
//! Forchheimer/slip correction, moving mesh, deformation, unsaturated mobility,
//! thermal deposition or Biot wave.
//! These entries inherit the existing synchronous Stokes solver's
//! iteration budget and cancellation limitations; they add no cancellation claim.

use super::boundary::{StokesAssembly, StokesBoundary, StokesBoundaryError, StokesSolution};
use super::{FluxParams, FluxSystem};
use crate::bdm::{CellBasis, eval_basis, tri_quad};
use fs_material::porous::permeability::IntrinsicPermeability;
use fs_sparse::Csr;

mod transient;

/// One cell's fixed geometry and explicitly selected Brinkman shear closure.
/// The original permeability tensor and its world axes remain intact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrinkmanCell {
    permeability: Option<IntrinsicPermeability>,
    effective_viscosity_ratio: f64,
    skeleton_velocity: [f64; 2],
}

impl BrinkmanCell {
    /// A resolved free-fluid cell: molecular stress and exactly zero pore drag.
    #[must_use]
    pub const fn free_fluid() -> Self {
        Self {
            permeability: None,
            effective_viscosity_ratio: 1.0,
            skeleton_velocity: [0.0; 2],
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
            skeleton_velocity: [0.0; 2],
        })
    }

    /// Prescribe this cell's uniform skeleton velocity [m/s] in world axes.
    ///
    /// Drag acts on `u-u_s` and its reaction can perform mechanical work. This
    /// selects the one-domain/high-porosity Brinkman velocity convention
    /// (phi=1), not a finite-porosity fluid-only superficial flux. The supplied
    /// geometry is a frozen snapshot: this does not move interfaces, evolve
    /// porosity, deform the solid or solve for the skeleton's dynamics.
    /// Zero velocity retains the original stationary arithmetic exactly.
    ///
    /// # Errors
    /// Refuses nonfinite velocity or nonzero skeleton motion in a free-fluid
    /// cell, which has no skeleton to prescribe.
    pub fn with_skeleton_velocity(
        mut self,
        velocity: [f64; 2],
    ) -> Result<Self, StokesBoundaryError> {
        if !velocity.iter().all(|value| value.is_finite())
            || self.permeability.is_none() && velocity != [0.0; 2]
        {
            return Err(invalid(
                "Brinkman skeleton motion requires finite velocity in a porous cell",
            ));
        }
        self.skeleton_velocity = velocity;
        Ok(self)
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

    /// Authored skeleton velocity [m/s] in the same world frame as the flow.
    #[must_use]
    pub const fn skeleton_velocity(&self) -> [f64; 2] {
        self.skeleton_velocity
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
    skeleton_velocity: [f64; 2],
}

impl DarcyDrag {
    fn compile(cell: BrinkmanCell, nu: f64) -> Result<Self, StokesBoundaryError> {
        let Some(permeability) = cell.permeability else {
            return Ok(Self {
                rates: [0.0; 3],
                axes: [[0.0; 2]; 3],
                skeleton_velocity: [0.0; 2],
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
        Ok(Self {
            rates,
            axes,
            skeleton_velocity: cell.skeleton_velocity,
        })
    }

    pub(super) fn local_matrix_and_load(
        &self,
        basis: &CellBasis,
        vertices: [[f64; 2]; 3],
        area: f64,
    ) -> ([[f64; 6]; 6], [f64; 6]) {
        let mut matrix = [[0.0; 6]; 6];
        let mut load = [0.0; 6];
        if self.rates == [0.0; 3] {
            return (matrix, load);
        }
        let moving = self.skeleton_velocity != [0.0; 2];
        let skeleton = self
            .axes
            .map(|axis| axis[0] * self.skeleton_velocity[0] + axis[1] * self.skeleton_velocity[1]);
        for (point, weight) in tri_quad(vertices, area) {
            let projected: [[f64; 3]; 6] = core::array::from_fn(|i| {
                let value = eval_basis(basis, i, point);
                self.axes
                    .map(|axis| axis[0] * value[0] + axis[1] * value[1])
            });
            for i in 0..6 {
                if moving {
                    // The same spectral factors and quadrature supply D*u_s
                    // to the original load/reaction rows, once per cell.
                    let value: f64 = (0..3)
                        .map(|a| self.rates[a] * projected[i][a] * skeleton[a])
                        .sum();
                    load[i] += weight * value;
                }
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
        (matrix, load)
    }

    fn power(&self, velocity: [f64; 2]) -> f64 {
        let relative = [
            velocity[0] - self.skeleton_velocity[0],
            velocity[1] - self.skeleton_velocity[1],
        ];
        (0..3)
            .map(|a| {
                let projection = self.axes[a][0] * relative[0] + self.axes[a][1] * relative[1];
                self.rates[a] * projection * projection
            })
            .sum()
    }

    fn exchange(&self, velocity: [f64; 2]) -> BrinkmanDragExchange {
        let relative = [
            velocity[0] - self.skeleton_velocity[0],
            velocity[1] - self.skeleton_velocity[1],
        ];
        let mut exchange = BrinkmanDragExchange::default();
        for (rate, axis) in self.rates.into_iter().zip(self.axes) {
            let projection = axis[0] * relative[0] + axis[1] * relative[1];
            let effort = rate * projection;
            for a in 0..2 {
                exchange.skeleton_force[a] += effort * axis[a];
            }
            exchange.fluid_power -= effort * (axis[0] * velocity[0] + axis[1] * velocity[1]);
            exchange.skeleton_power += effort
                * (axis[0] * self.skeleton_velocity[0] + axis[1] * self.skeleton_velocity[1]);
            exchange.dissipation += effort * projection;
        }
        exchange
    }
}

/// Area-integrated Darcy exchange divided by the same reference density.
///
/// Force is per out-of-plane depth/density [m³/s²], and powers are per
/// depth/density [m⁴/s³]. Positive power is work ON the named subsystem.
/// `fluid_power + skeleton_power + dissipation = 0` to quadrature/roundoff;
/// a driven skeleton may give the fluid positive work while dissipating heat.
/// These are constitutive exchanges, not automatic thermal/body updates or
/// a complete energy balance for the spatial discretization.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BrinkmanDragExchange {
    /// Force exerted by the fluid ON the skeleton; fluid force is its negative.
    pub skeleton_force: [f64; 2],
    /// Signed work rate of Darcy force on the fluid.
    pub fluid_power: f64,
    /// Signed work rate of the opposite reaction on the prescribed skeleton.
    pub skeleton_power: f64,
    /// Nonnegative `integral(nu (u-u_s).K^-1.(u-u_s))`.
    pub dissipation: f64,
}

impl BrinkmanDragExchange {
    fn add_scaled(&mut self, value: Self, weight: f64) {
        for a in 0..2 {
            self.skeleton_force[a] += weight * value.skeleton_force[a];
        }
        self.fluid_power += weight * value.fluid_power;
        self.skeleton_power += weight * value.skeleton_power;
        self.dissipation += weight * value.dissipation;
    }

    fn finite(self) -> bool {
        self.skeleton_force.iter().all(|value| value.is_finite())
            && self.fluid_power.is_finite()
            && self.skeleton_power.is_finite()
            && self.dissipation.is_finite()
            && self.dissipation >= 0.0
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

    /// Integrate `nu (u-u_s).K^-1.(u-u_s)` using the assembled coefficients.
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

    /// Integrate the reaction, signed subsystem work and positive relative
    /// dissipation in each original triangle. Uses the same spectral factors,
    /// velocity reconstruction and quadrature as the assembled Darcy operator.
    /// A common frame translation preserves reaction and dissipation; signed
    /// subsystem powers change by equal and opposite force-times-velocity work.
    ///
    /// # Errors
    /// Refuses incorrect velocity dimensions, nonfinite states or overflow.
    pub fn cell_darcy_exchange_per_density(
        &self,
        velocity: &[f64],
    ) -> Result<Vec<BrinkmanDragExchange>, StokesBoundaryError> {
        if velocity.len() != self.system.n_u || velocity.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "Darcy exchange requires exactly the finite native velocity DOFs",
            ));
        }
        let mesh = self.system.mesh;
        let mut exchanges = Vec::with_capacity(mesh.tris.len());
        for (cell, drag) in self.drag.iter().enumerate() {
            let vertices = core::array::from_fn(|i| mesh.verts[mesh.tris[cell][i]]);
            let mut exchange = BrinkmanDragExchange::default();
            for (point, weight) in tri_quad(vertices, mesh.areas[cell]) {
                let mut value = [0.0; 2];
                for i in 0..6 {
                    let basis = eval_basis(&self.system.bases[cell], i, point);
                    let dof = velocity[self.system.dof(cell, i)];
                    value[0] += dof * basis[0];
                    value[1] += dof * basis[1];
                }
                exchange.add_scaled(drag.exchange(value), weight);
            }
            if !exchange.finite() {
                return Err(invalid("Darcy force or work is not representable"));
            }
            exchanges.push(exchange);
        }
        Ok(exchanges)
    }
}

// Compile material coefficients once per solve, outside the quadrature loops.
// All native Brinkman entries share this admission and molecular/shear split.
fn compile_cells(
    cells: &[BrinkmanCell],
    expected_cells: usize,
    nu: f64,
) -> Result<(Vec<DarcyDrag>, Vec<f64>), StokesBoundaryError> {
    if cells.len() != expected_cells || !positive(nu) {
        return Err(invalid(
            "Brinkman requires one cell law per triangle and positive molecular viscosity",
        ));
    }
    let mut drag = Vec::with_capacity(cells.len());
    let mut effective_nu = Vec::with_capacity(cells.len());
    for &cell in cells {
        drag.push(DarcyDrag::compile(cell, nu)?);
        let effective = nu * cell.effective_viscosity_ratio;
        if !positive(effective) {
            return Err(invalid("effective Brinkman viscosity is not representable"));
        }
        effective_nu.push(effective);
    }
    Ok((drag, effective_nu))
}

impl<'m> FluxSystem<'m> {
    /// Assemble steady, fully saturated Stokes–Brinkman momentum.
    ///
    /// `params.nu` is molecular fluid viscosity/density. Each cell provides an
    /// intrinsic tensor and separate effective shear-viscosity ratio; free-fluid
    /// cells have no drag. Prescribed skeleton velocity contributes the matching
    /// `nu K^-1 u_s` load; it does not change the matrix or stress law.
    /// No permeability component is replaced by a diagonal
    /// approximation: all in-plane couplings enter the velocity block together.
    /// Cellwise coefficients are compiled once outside all quadrature loops.
    /// Strong normal moments, total-traction boundaries, SIP interfaces, exact
    /// cell continuity and both true residual checks remain owned by Stokes.
    /// A positive porous cell removes all rigid velocity modes of its connected
    /// domain, admitting all-traction boundaries there without an artificial
    /// velocity anchor. Every component without pore drag independently
    /// requires a velocity boundary.
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
        let (drag, effective_nu) = compile_cells(cells, self.mesh.tris.len(), params.nu)?;
        let mut anchors = vec![false; self.pressure_component_count()];
        for (cell, law) in cells.iter().enumerate() {
            anchors[self.pressure_component[cell]] |= law.permeability.is_some();
        }
        let admitted =
            self.admit_stokes_problem_with_velocity_anchors(params, boundaries, &anchors)?;
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
