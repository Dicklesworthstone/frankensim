//! Checked linear unsteady Stokes on the original mixed-boundary operator.

use super::{StokesAssembly, StokesBoundary, StokesBoundaryError, StokesSolution};
use crate::ns::{FluxParams, FluxSystem};

impl<'m> FluxSystem<'m> {
    /// Assemble one backward-Euler step of linear incompressible Stokes:
    /// `(u - u_previous)/dt - div(2 nu D(u)) + grad(p) = force`.
    ///
    /// Force, velocity and TOTAL-traction callbacks specify the new time level.
    /// Capture that time in the callbacks to change applied pressure or stress
    /// between steps. Previous velocity contains exactly the original BDM1
    /// velocity moments; no previous pressure is needed. This uses the same
    /// mass quadrature, symmetric stress and strong-normal/SIP-tangent boundary
    /// treatment as the existing transient operator, without convection.
    ///
    /// The positive mass term supplies physical inertia for entirely natural
    /// components, so a velocity boundary is unnecessary at finite dt. Every
    /// traction component retains all continuity equations and its physical
    /// traction datum. Sealed components retain independent mean-zero gauges
    /// and must satisfy their own endpoint flux compatibility.
    ///
    /// Geometry is metres, time seconds, force acceleration [m/s²], and stress
    /// and pressure are divided by the same constant reference density. This
    /// is fixed-mesh linear Stokes, not nonlinear open-boundary Navier–Stokes.
    /// The original synchronous solver's cancellation limitations still apply.
    ///
    /// # Errors
    /// Refuses invalid viscosity, history dimensions/values or a nonpositive,
    /// nonfinite or noninvertible timestep before querying callbacks. All
    /// mixed-boundary coefficient, ownership and component-flux checks remain.
    pub fn assemble_stokes_mixed_bdf1<'a>(
        &'a self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Stokes viscosity",
            ));
        }
        self.assemble_stokes_mixed_bdf1_variable_viscosity(
            params,
            &|_, _| params.nu,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
        )
    }

    /// The same BDF1 mixed-Stokes step with positive endpoint kinematic
    /// viscosity `viscosity(cell, point)` [m²/s]; `params.nu` is unused.
    /// The callback must be pure and use the cell index for one-sided traces.
    /// Its field is frozen during this linear step, with constant reference
    /// density at both time levels. No implicit rheological update is inferred.
    ///
    /// # Errors
    /// Inherits [`Self::assemble_stokes_mixed_bdf1`] admission and checked
    /// variable-viscosity sampling. Refusal never changes the history input.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_stokes_mixed_bdf1_variable_viscosity<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        self.assemble_stokes_mixed_bdf1_inertia(
            params,
            None,
            viscosity,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
        )
    }

    /// Linear unsteady Stokes with a fixed, spatially varying mass density:
    /// `r (u - u_previous)/dt - div(2 nu_reference D(u)) + grad(pi) = f`.
    ///
    /// Choose ONE constant positive reference density rho_reference [kg/m³].
    /// `relative_density(cell, point)` is r = rho(x)/rho_reference, and
    /// `viscosity(cell, point)` is mu(x)/rho_reference [m²/s], NOT mu/rho(x).
    /// Force density [N/m³], pressure [Pa], total traction [Pa] and any Navier
    /// slip resistance must be divided by that SAME reference density.
    /// Returned pressure and reactions retain this common normalization.
    /// `params.nu` is unused. Callbacks must be pure; the cell index admits
    /// one-sided material jumps without smearing density across an interface.
    ///
    /// Density is sampled once at each positive mass quadrature point. The
    /// consistent matrix and immutable previous-velocity load use exactly the
    /// same stored weights, including supported reaction rows. BDM1 velocity
    /// is affine, so the degree-four rule integrates density polynomials of
    /// degree at most two per cell exactly. General density fields retain the
    /// ordinary spatial-quadrature/refinement obligation. The assembly's
    /// [`StokesAssembly::kinetic_energy`] uses this same metric without resampling.
    ///
    /// The density field is fixed at both time levels. This covers linear
    /// perturbations about a stationary stratification and shear parallel to
    /// density layers; it does not transport density or material interfaces,
    /// add nonlinear momentum advection, or model mixing/compressibility.
    ///
    /// # Errors
    /// Refuses invalid history/timestep, nonpositive/nonfinite density samples,
    /// unrepresentable weighted mass, and the original mixed-boundary failures.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_stokes_mixed_bdf1_variable_density<'a>(
        &'a self,
        params: FluxParams,
        relative_density: &dyn Fn(usize, [f64; 2]) -> f64,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        self.assemble_stokes_mixed_bdf1_inertia(
            params,
            Some(relative_density),
            viscosity,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_stokes_mixed_bdf1_inertia<'a>(
        &'a self,
        params: FluxParams,
        relative_density: Option<&dyn Fn(usize, [f64; 2]) -> f64>,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        let inv_dt = timestep_s.recip();
        if !timestep_s.is_finite()
            || timestep_s <= 0.0
            || !inv_dt.is_finite()
            || inv_dt <= 0.0
            || previous_velocity.len() != self.n_u
            || previous_velocity.iter().any(|value| !value.is_finite())
        {
            return Err(StokesBoundaryError::InvalidInput(
                "BDF1 needs finite velocity moments and a representable positive timestep in seconds",
            ));
        }
        // M/dt is positive definite on every velocity component, including
        // rigid translation/rotation. This is physical transient momentum,
        // not an artificial steady pin or porous resistance.
        let inertia_anchors = vec![true; self.pressure_roots.len()];
        let admitted =
            self.admit_stokes_problem_with_velocity_anchors(params, boundaries, &inertia_anchors)?;
        self.assemble_stokes_admitted_density(
            params,
            viscosity,
            force,
            &admitted,
            None,
            inv_dt,
            Some(previous_velocity),
            None,
            None,
            relative_density,
        )
    }

    /// Advance one checked constant-viscosity linear mixed-Stokes step.
    ///
    /// This assembles [`Self::assemble_stokes_mixed_bdf1`] and applies the
    /// existing dense-LU/sparse-FGMRES selection. Both original and row-scaled
    /// residuals must converge before any candidate is returned. To inspect
    /// the matrix or use a spatial viscosity field, use the assembly entry.
    ///
    /// # Errors
    /// Refuses invalid inputs, unrepresentable assembly, singular systems or
    /// failed linear convergence without modifying `previous_velocity`.
    pub fn bdf1_stokes_mixed_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesSolution, StokesBoundaryError> {
        self.assemble_stokes_mixed_bdf1(params, force, boundaries, previous_velocity, timestep_s)?
            .solve()
    }
}
