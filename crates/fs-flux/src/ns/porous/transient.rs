//! Native backward Euler with fully implicit anisotropic skeleton drag.

use super::{BrinkmanAssembly, BrinkmanCell, compile_cells};
use crate::ns::boundary::{StokesBoundary, StokesBoundaryError};
use crate::ns::{FluxParams, FluxSystem};

impl<'m> FluxSystem<'m> {
    /// Assemble one backward-Euler step of incompressible Stokes–Brinkman:
    /// `(u-u_previous)/dt - div(2 nu_eff D(u)) + nu K^-1 (u-u_s) + grad(pi) = f`.
    ///
    /// Previous velocity contains exactly the original BDM1 velocity moments.
    /// The full anisotropic resistance is implicit; cell coefficients, prescribed
    /// skeleton velocities, force and mixed boundaries specify the new time level
    /// and remain frozen during this linear solve. No drag timestep restriction
    /// is introduced by explicit splitting. Time accuracy is first order.
    ///
    /// Native mass, pressure/continuity, SIP interfaces, boundary reactions and
    /// both true-residual checks retain their original owners. Finite positive
    /// inertia anchors even free-fluid all-traction components, without adding
    /// artificial steady resistance. Cell Darcy diagnostics are endpoint powers;
    /// backward-Euler numerical energy loss is not physical drag heating.
    ///
    /// Geometry is fixed, density constant, and the macroscopic velocity uses
    /// the one-domain/phi=1 convention of [`BrinkmanCell`]. This adds fluid inertia
    /// to linear Stokes, without convection, pore-scale added mass, finite-porosity
    /// storage, moving interfaces or a skeleton dynamics solver. The caller owns
    /// the clock and accepted history; neither is changed by assembly or refusal.
    ///
    /// # Errors
    /// Refuses invalid native velocity history or timestep, mismatched cell laws,
    /// unrepresentable coefficients and every original mixed-boundary admission
    /// failure. Invalid history and coefficients are rejected before callbacks.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_brinkman_mixed_bdf1<'a>(
        &'a self,
        params: FluxParams,
        cells: &[BrinkmanCell],
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<BrinkmanAssembly<'a, 'm>, StokesBoundaryError> {
        let inv_dt = self.admit_mixed_bdf1_history(previous_velocity, timestep_s)?;
        let (drag, effective_nu) = compile_cells(cells, self.mesh.tris.len(), params.nu)?;
        let inertia_anchors = vec![true; self.pressure_component_count()];
        let admitted =
            self.admit_stokes_problem_with_velocity_anchors(params, boundaries, &inertia_anchors)?;
        let stokes = self.assemble_stokes_admitted_operator(
            params,
            &|cell, _| effective_nu[cell],
            force,
            &admitted,
            Some(&drag),
            inv_dt,
            Some(previous_velocity),
            None,
        )?;
        Ok(BrinkmanAssembly {
            system: self,
            stokes,
            drag,
        })
    }
}
