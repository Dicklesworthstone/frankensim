//! Sparse spectral windows about the same finitely deformed equilibrium.
//!
//! The original hyper3 owner supplies equilibrium, supports, reference inertia
//! and the complete material/follower tangent. The existing fs-modal window
//! owner supplies shared sparse factors, signed inertia counts, cancellation,
//! residual intervals and complete-in-window certificates. No dense n-by-n
//! allocation or 128-DOF limit is introduced by this adapter.
//!
//! Windows are in signed squared-frequency units [s^-2], with the native
//! (low, high] convention. Negative modes remain negative. A positive window
//! does not establish stability outside that window; omitted modes do not
//! become a certified dynamics reduction. Singular spectral endpoints refuse
//! without being shifted. The native residual/abs(lambda) gate can refuse a
//! window containing rigid or nearly neutral modes; the bounded complete dense
//! route remains available for those small pencils. Neither route establishes
//! continuum stability or resolves the original P1 volumetric locking.
//!
//! Equilibrium/assembly retain the original TetAssemblyBudget. WindowPlan and
//! SliceOptions retain their separate spectral payload, factor, iteration and
//! worker budgets, including indivisible sparse symbolic/factorization leaves.
//! These are native numerical payload limits, not a process-memory guarantee.

use super::{ReferenceDensityInput, RetainedReferenceDensity, VibrationError, invalid};
use crate::hyper3::pressure::ClosedTetPressureSurface;
use crate::hyper3::{HyperTetError, HyperTetProblem, HyperTetSolution, checkpoint};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_modal::{MultiSliceReport, SliceOptions, WindowPlan, slice_windows};
use fs_sparse::Csr;

/// One accepted equilibrium, its full sparse free-coordinate pencil, and every
/// mode in the requested windows. No eigenvalue outside those windows is
/// inferred, clipped, or replaced by a real frequency.
#[derive(Debug)]
pub struct HyperTetVibrationWindows<'a, L: HyperelasticMaterial = Hyperelastic> {
    problem: HyperTetProblem<'a, L>,
    held_pressure: Option<(&'a ClosedTetPressureSurface, f64)>,
    reference_density: RetainedReferenceDensity,
    equilibrium: HyperTetSolution,
    free_dofs: Vec<usize>,
    stiffness_n_m: Csr,
    mass_kg: Csr,
    plan: WindowPlan,
    options: SliceOptions,
    spectrum: MultiSliceReport,
}

// Rows/columns are already canonical and the free map is increasing. Merge
// the native material/load rows directly: no COO sort, dense intermediate,
// stiffness reconstruction or second initial-stress matrix is needed.
fn restrict(
    matrix: &Csr,
    subtract: Option<&Csr>,
    free: &[usize],
    reduced: &[usize],
    cx: &Cx<'_>,
) -> Result<Csr, HyperTetError> {
    let n = free.len();
    let mut rows = Vec::with_capacity(n + 1);
    let mut columns = Vec::new();
    let mut values = Vec::new();
    let mut scale = 0.0_f64;
    rows.push(0);
    for &row in free {
        checkpoint(cx)?;
        let (ac, av) = matrix.row(row);
        let (bc, bv) = subtract.map_or((&[][..], &[][..]), |load| load.row(row));
        let (mut a, mut b) = (0, 0);
        while a < ac.len() || b < bc.len() {
            let col = ac
                .get(a)
                .copied()
                .unwrap_or(usize::MAX)
                .min(bc.get(b).copied().unwrap_or(usize::MAX));
            let original = if ac.get(a) == Some(&col) {
                let value = av[a];
                a += 1;
                value
            } else {
                0.0
            };
            let value = if bc.get(b) == Some(&col) {
                let value = original - bv[b];
                b += 1;
                value
            } else {
                original
            };
            if reduced[col] != usize::MAX {
                if !value.is_finite() {
                    return Err(invalid("nonfinite finite-strain vibration pencil"));
                }
                columns.push(reduced[col]);
                values.push(value);
                scale = scale.max(value.abs());
            }
        }
        rows.push(values.len());
    }
    let result =
        Csr::try_from_parts_with_checkpoint(n, n, rows, columns, values, || checkpoint(cx))?
            .ok_or_else(|| invalid("noncanonical native vibration restriction"))?;
    let tolerance = (64.0 * f64::EPSILON * n as f64) * scale;
    for row in 0..n {
        checkpoint(cx)?;
        let (columns, values) = result.row(row);
        for (&col, &value) in columns.iter().zip(values) {
            if (value - result.get(col, row)).abs() > tolerance {
                return Err(invalid(
                    "native vibration pencil is not symmetric to roundoff",
                ));
            }
        }
    }
    Ok(result)
}

impl<'a, L: HyperelasticMaterial> HyperTetProblem<'a, L> {
    /// Solve the original dead-load equilibrium and certify all requested
    /// sparse spectral windows. The original assembly budget bounds free DOFs;
    /// this route has no dense-spectrum cap.
    ///
    /// Supply one positive reference density per original tetrahedron. Window
    /// coordinates, when supplied for nested dissection, use ascending free
    /// node-major xyz DOF order (repeat a node position for its free components).
    /// Material/source/reference adapters remain in the equilibrium and tangent.
    ///
    /// # Errors
    /// Native equilibrium, mass, material, symmetry, spectral budget, unresolved
    /// window and cancellation refusals. No partial window report is returned.
    pub fn solve_vibration_windows(
        &self,
        reference_density_kg_m3: &[f64],
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        self.solve_vibration_windows_impl(
            ReferenceDensityInput::PerElement(reference_density_kg_m3),
            None,
            None,
            plan,
            options,
            cx,
        )
    }

    /// The same signed sparse windows with a positive P1 reference density on
    /// the original mesh nodes. Native cubic moment integration supplies the
    /// complete consistent mass, without centroid sampling or cell averaging.
    /// The retained nodal field and all original spectral budgets are preserved.
    ///
    /// # Errors
    /// Nodal-density admission and all errors of [`Self::solve_vibration_windows`].
    pub fn solve_vibration_windows_with_nodal_density(
        &self,
        reference_nodal_density_kg_m3: &[f64],
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        self.solve_vibration_windows_impl(
            ReferenceDensityInput::Nodal(reference_nodal_density_kg_m3),
            None,
            None,
            plan,
            options,
            cx,
        )
    }

    /// Re-equilibrate one dead-load endpoint from an initial iterate, then
    /// certify the selected signed windows. The original solve_from one-step
    /// rule applies; the iterate is never treated as an accepted equilibrium.
    /// Unstable branches may be queried with negative eigenvalue windows.
    pub fn solve_vibration_windows_from(
        &self,
        initial_displacement_m: &[f64],
        reference_density_kg_m3: &[f64],
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        self.solve_vibration_windows_impl(
            ReferenceDensityInput::PerElement(reference_density_kg_m3),
            Some(initial_displacement_m),
            None,
            plan,
            options,
            cx,
        )
    }

    /// Re-equilibrate from an initial iterate while retaining the original P1
    /// nodal density in every requested signed window. The initial iterate and
    /// one-load-step rules of [`Self::solve_vibration_windows_from`] apply.
    ///
    /// # Errors
    /// Nodal-density and native initial-iterate/spectral refusals.
    pub fn solve_vibration_windows_from_with_nodal_density(
        &self,
        initial_displacement_m: &[f64],
        reference_nodal_density_kg_m3: &[f64],
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        self.solve_vibration_windows_impl(
            ReferenceDensityInput::Nodal(reference_nodal_density_kg_m3),
            Some(initial_displacement_m),
            None,
            plan,
            options,
            cx,
        )
    }

    /// Solve with held uniform closed pressure, then certify the windows of
    /// K = K_material - d(f_pressure)/du = K_material + p*d²V/du².
    /// Positive pressure is inward; the exact surface/value remain bound to
    /// the report. A finite gas reservoir or nonconservative pressure field is
    /// not inferred. All original reference/pressure gates remain in force.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_vibration_windows_with_pressure(
        &self,
        surface: &'a ClosedTetPressureSurface,
        pressure_pa: f64,
        reference_density_kg_m3: &[f64],
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        self.solve_vibration_windows_impl(
            ReferenceDensityInput::PerElement(reference_density_kg_m3),
            None,
            Some((surface, pressure_pa)),
            plan,
            options,
            cx,
        )
    }

    /// Original closed-pressure equilibrium and signed sparse windows with a
    /// positive P1 nodal reference density. Pressure contributes its original
    /// follower tangent exactly once; finite deformation never refreshes mass.
    ///
    /// # Errors
    /// Nodal-density and native pressure/equilibrium/spectral refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_vibration_windows_with_nodal_density_and_pressure(
        &self,
        surface: &'a ClosedTetPressureSurface,
        pressure_pa: f64,
        reference_nodal_density_kg_m3: &[f64],
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        self.solve_vibration_windows_impl(
            ReferenceDensityInput::Nodal(reference_nodal_density_kg_m3),
            None,
            Some((surface, pressure_pa)),
            plan,
            options,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_vibration_windows_impl(
        &self,
        density: ReferenceDensityInput<'_>,
        initial: Option<&[f64]>,
        held_pressure: Option<(&'a ClosedTetPressureSurface, f64)>,
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibrationWindows<'a, L>, VibrationError> {
        checkpoint(cx)?;
        self.validate_inputs()?;
        if let Some((surface, pressure)) = held_pressure {
            if !pressure.is_finite() || !surface.matches_reference(self.nodes_m, self.tetrahedra) {
                return Err(invalid(
                    "finite pressure and the exact reference pressure surface are required",
                )
                .into());
            }
        }
        let mut fixed = vec![false; 3 * self.nodes_m.len()];
        for &(dof, _) in self.prescribed_m {
            fixed[dof] = true;
        }
        let free_dofs: Vec<_> = (0..fixed.len()).filter(|&dof| !fixed[dof]).collect();
        if free_dofs.is_empty() {
            return Err(invalid("finite-strain vibration requires free displacement DOFs").into());
        }
        let mut reduced = vec![usize::MAX; fixed.len()];
        for (index, &dof) in free_dofs.iter().enumerate() {
            reduced[dof] = index;
        }
        let mass_kg = restrict(
            &density.mass_matrix(self, cx)?,
            None,
            &free_dofs,
            &reduced,
            cx,
        )?;
        let equilibrium = match held_pressure {
            Some((surface, pressure)) => self.solve_with_pressure(surface, pressure, cx)?,
            None => match initial {
                None => self.solve(cx)?,
                Some(displacement) => self.solve_from(displacement, cx)?,
            },
        };
        let displacement: Vec<_> = equilibrium
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        let pressure_tangent = held_pressure
            .filter(|(_, pressure)| *pressure != 0.0)
            .map(|(surface, pressure)| surface.force_tangent(&displacement, pressure, cx))
            .transpose()?;
        let stiffness_n_m = restrict(
            &self.prepare(cx)?.tangent_matrix(&displacement, cx)?,
            pressure_tangent.as_ref(),
            &free_dofs,
            &reduced,
            cx,
        )?;
        // Free the full load operator and kinematic scratch before the spectral
        // owner's independent payload admission. No dense matrix is retained.
        drop(pressure_tangent);
        drop(displacement);
        drop(reduced);
        drop(fixed);
        let spectrum = slice_windows(&stiffness_n_m, &mass_kg, plan, options, cx)
            .map_err(VibrationError::Modal)?;
        let reference_density = density.retain(self.tetrahedra, cx)?;
        let result = HyperTetVibrationWindows {
            problem: *self,
            held_pressure,
            reference_density,
            equilibrium,
            free_dofs,
            stiffness_n_m,
            mass_kg,
            plan: plan.clone(),
            options: *options,
            spectrum,
        };
        checkpoint(cx)?;
        Ok(result)
    }
}

impl<L: HyperelasticMaterial> HyperTetVibrationWindows<'_, L> {
    /// Original mesh, material/source adapters, supports and dead loads.
    pub const fn problem(&self) -> &HyperTetProblem<'_, L> {
        &self.problem
    }

    /// Exact closed surface and held signed pressure [Pa]; positive is inward.
    pub const fn held_pressure(&self) -> Option<(&ClosedTetPressureSurface, f64)> {
        self.held_pressure
    }

    /// Conserved mass per original volume, in tetrahedron order. For a graded
    /// field these are inventory means; the retained mass uses all nodal moments.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density.cell_means_kg_m3
    }

    /// Complete original P1 nodal density, or None for cellwise constant density.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_density.nodal_kg_m3.as_deref()
    }

    /// Accepted native equilibrium; no stress reset or foreign state injection.
    pub const fn equilibrium(&self) -> &HyperTetSolution {
        &self.equilibrium
    }

    /// Ascending reduced-to-full node-major xyz DOFs; fixed increments are zero.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Full sparse free-coordinate material minus held-pressure tangent [N/m].
    pub const fn stiffness_n_m(&self) -> &Csr {
        &self.stiffness_n_m
    }

    /// Full sparse free-coordinate consistent original reference mass [kg].
    pub const fn mass_kg(&self) -> &Csr {
        &self.mass_kg
    }

    /// Exact requested coverage, spectral budgets and residual gate.
    pub const fn plan(&self) -> &WindowPlan {
        &self.plan
    }

    /// Original sparse factorization, ordering and Lanczos controls.
    pub const fn options(&self) -> &SliceOptions {
        &self.options
    }

    /// Native signed modes, endpoint inertia counts, residual intervals and
    /// numerical work receipt. Complete only within the requested windows;
    /// positive selected modes alone never certify the whole equilibrium stable.
    pub const fn spectrum(&self) -> &MultiSliceReport {
        &self.spectrum
    }
}
