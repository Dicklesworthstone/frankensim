//! Steady and BDF1 isotropic form drag in the original anisotropic Brinkman owner.
//!
//! The added physical resistance is `rho beta |w| w`, `w = u-u_s`, beta [1/m].
//! The original molecular-viscosity Darcy tensor is retained independently.
//! This convention agrees with `fs_material::porous::permeability`'s directional
//! Forchheimer law. Its gas/specimen-specific calibration is deliberately not
//! converted to an isotropic two-dimensional closure without that authority.
//! See https://doc.comsol.com/6.3/doc/com.comsol.help.cfd/cfd_ug_fluidflow_porous.10.004.html.
//!
//! These are saturated, incompressible, fixed-mesh closures with the parent's
//! one-domain/phi=1 velocity convention. BDF1 adds native fluid inertia only.
//! There is no convective momentum, pore-scale added mass, finite-porosity
//! storage, turbulent stress, moving geometry, solid evolution or automatic
//! heat deposition. Authored beta alone is not a material calibration or an
//! applicability certificate. The caller owns validation over its flow regime.

use super::{
    BrinkmanAssembly, BrinkmanCell, BrinkmanDragExchange, CellDragTerms, CellResistance, DarcyDrag,
    compile_cells, invalid,
};
use crate::bdm::{eval_basis, tri_quad};
use crate::ns::boundary::{StokesAssembly, StokesBoundary, StokesBoundaryError, StokesSolution};
use crate::ns::{FluxParams, FluxRheologyError, FluxSolution, FluxSystem};

mod law;

/// One original Brinkman cell plus an explicitly authored isotropic form-drag
/// coefficient [1/m], with no hidden factor of 1/2 or permeability conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForchheimerCell {
    brinkman: BrinkmanCell,
    beta_per_m: f64,
}

impl ForchheimerCell {
    /// Retain all anisotropic Darcy and moving-skeleton data. Positive beta is
    /// meaningful only in a porous cell; zero beta retains the original law.
    ///
    /// # Errors
    /// Refuses negative/nonfinite beta or form drag in a resolved free fluid.
    pub fn try_new(brinkman: BrinkmanCell, beta_per_m: f64) -> Result<Self, StokesBoundaryError> {
        if !beta_per_m.is_finite()
            || beta_per_m < 0.0
            || (beta_per_m > 0.0 && brinkman.permeability().is_none())
        {
            return Err(invalid(
                "Forchheimer beta must be finite, nonnegative and belong to a porous cell",
            ));
        }
        Ok(Self {
            brinkman,
            beta_per_m,
        })
    }

    /// Original permeability, molecular/effective-viscosity split and skeleton.
    #[must_use]
    pub const fn brinkman(&self) -> &BrinkmanCell {
        &self.brinkman
    }

    /// Authored isotropic inertial resistance [1/m].
    #[must_use]
    pub const fn beta_per_m(&self) -> f64 {
        self.beta_per_m
    }
}

/// Finite budgets for a consistent Newton solve with residual backtracking.
/// Linear solves retain `FluxParams` and must have at least this accuracy.
#[derive(Clone, Copy, Debug)]
pub struct ForchheimerControls {
    /// Maximum Newton linear solves, in 1..=256.
    pub max_iterations: usize,
    /// Maximum halvings per Newton direction, in 0..=60; zero tries only the
    /// full step. Each trial is checked against the actual nonlinear equations.
    pub max_backtracks: usize,
    /// Required original and row-scaled nonlinear residuals AND unrelaxed
    /// relative velocity correction, using the fixed references documented by
    /// the solve entry. Finite and strictly between zero and one.
    pub relative_tolerance: f64,
}

impl Default for ForchheimerControls {
    fn default() -> Self {
        Self {
            max_iterations: 40,
            max_backtracks: 24,
            relative_tolerance: 1e-9,
        }
    }
}

impl ForchheimerControls {
    fn admit(self, params: FluxParams) -> Result<(), StokesBoundaryError> {
        if !(1..=256).contains(&self.max_iterations)
            || self.max_backtracks > 60
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance <= 0.0
            || self.relative_tolerance >= 1.0
            || !params.tol.is_finite()
            || params.tol <= 0.0
            || params.tol > self.relative_tolerance
        {
            return Err(invalid(
                "invalid Forchheimer Newton budgets/tolerance or insufficient linear tolerance",
            ));
        }
        Ok(())
    }
}

/// Accepted native flow with separate actual-state Darcy and inertial work.
/// No last-iterate frozen resistance is used for its residuals or reactions.
pub struct ForchheimerSolution {
    /// Original pressure/continuity/boundary owner. Both residual fields refer
    /// to the accepted NONLINEAR equations, including constrained rows.
    pub flow: StokesSolution,
    /// Newton solves performed; zero identifies an accepted initial state or
    /// the exact zero-beta linear path.
    pub iterations: usize,
    /// Total rejected trial steps before acceptance.
    pub line_search_reductions: usize,
    /// Last FULL Newton velocity correction relative to the fixed maximum of
    /// initial and first Newton candidate velocity norms. No velocity floor or
    /// small relaxed step can alone establish convergence.
    pub relative_update: f64,
    /// Actual accepted Darcy exchange, divided by reference density, per cell.
    pub cell_darcy_exchange: Vec<BrinkmanDragExchange>,
    /// Actual accepted `beta |w| w` exchange, divided by reference density,
    /// per cell. Dissipation is the quadrature of `beta |w|^3`.
    pub cell_forchheimer_exchange: Vec<BrinkmanDragExchange>,
}

// Compile the native Darcy matrix/load once per cell. Only the nonlinear
// tangent changes between Newton steps; no dense global matrix is formed here.
fn darcy_terms(system: &FluxSystem<'_>, drag: &[DarcyDrag]) -> Vec<CellDragTerms> {
    drag.iter()
        .enumerate()
        .map(|(cell, coefficient)| {
            let vertices = core::array::from_fn(|i| system.mesh.verts[system.mesh.tris[cell][i]]);
            coefficient.local_matrix_and_load(
                &system.bases[cell],
                vertices,
                system.mesh.areas[cell],
            )
        })
        .collect()
}

fn linearized_terms(
    system: &FluxSystem<'_>,
    cells: &[ForchheimerCell],
    base: &[CellDragTerms],
    velocity: &[f64],
) -> Result<Vec<CellDragTerms>, StokesBoundaryError> {
    let mut terms = base.to_vec();
    for (cell, coefficient) in cells.iter().enumerate() {
        if coefficient.beta_per_m == 0.0 {
            continue;
        }
        let vertices = core::array::from_fn(|i| system.mesh.verts[system.mesh.tris[cell][i]]);
        let (matrix, load) = &mut terms[cell];
        for (point, weight) in tri_quad(vertices, system.mesh.areas[cell]) {
            let values: [[f64; 2]; 6] =
                core::array::from_fn(|i| eval_basis(&system.bases[cell], i, point));
            let mut u = [0.0; 2];
            for (i, basis) in values.iter().enumerate() {
                let dof = velocity[system.dof(cell, i)];
                u[0] += dof * basis[0];
                u[1] += dof * basis[1];
            }
            let response = law::evaluate(
                coefficient.beta_per_m,
                u,
                coefficient.brinkman.skeleton_velocity(),
            )
            .map_err(invalid)?;
            // Form J*v once per column; retain its directional coupling even
            // though beta itself is isotropic.
            let products = values.map(|trial| {
                [
                    response.tangent[0][0] * trial[0] + response.tangent[0][1] * trial[1],
                    response.tangent[1][0] * trial[0] + response.tangent[1][1] * trial[1],
                ]
            });
            for (i, test) in values.iter().enumerate() {
                load[i] += weight * (test[0] * response.load[0] + test[1] * response.load[1]);
                for (j, jv) in products.iter().enumerate().skip(i) {
                    let value = weight * (test[0] * jv[0] + test[1] * jv[1]);
                    matrix[i][j] += value;
                    if i != j {
                        matrix[j][i] += value;
                    }
                }
            }
        }
        if !matrix
            .iter()
            .flatten()
            .chain(load.iter())
            .all(|value| value.is_finite())
        {
            return Err(invalid(
                "Forchheimer cell tangent or load is not representable",
            ));
        }
    }
    Ok(terms)
}

fn forchheimer_exchanges(
    system: &FluxSystem<'_>,
    cells: &[ForchheimerCell],
    velocity: &[f64],
) -> Result<Vec<BrinkmanDragExchange>, StokesBoundaryError> {
    let mut exchanges = Vec::with_capacity(cells.len());
    for (cell, coefficient) in cells.iter().enumerate() {
        let mut exchange = BrinkmanDragExchange::default();
        if coefficient.beta_per_m > 0.0 {
            let vertices = core::array::from_fn(|i| system.mesh.verts[system.mesh.tris[cell][i]]);
            for (point, weight) in tri_quad(vertices, system.mesh.areas[cell]) {
                let mut u = [0.0; 2];
                for i in 0..6 {
                    let basis = eval_basis(&system.bases[cell], i, point);
                    let dof = velocity[system.dof(cell, i)];
                    u[0] += dof * basis[0];
                    u[1] += dof * basis[1];
                }
                let response = law::evaluate(
                    coefficient.beta_per_m,
                    u,
                    coefficient.brinkman.skeleton_velocity(),
                )
                .map_err(invalid)?;
                exchange.add_scaled(
                    BrinkmanDragExchange {
                        skeleton_force: response.force,
                        fluid_power: response.fluid_power,
                        skeleton_power: response.skeleton_power,
                        dissipation: response.dissipation,
                    },
                    weight,
                );
            }
        }
        if !exchange.finite() {
            return Err(invalid(
                "Forchheimer integrated force or work is not representable",
            ));
        }
        exchanges.push(exchange);
    }
    Ok(exchanges)
}

fn residual_norms(
    assembly: &StokesAssembly<'_, '_>,
    x: &[f64],
    fixed_row_scales: &[f64],
    work: &mut [f64],
) -> (f64, f64) {
    assembly.matrix().spmv(x, work);
    // Weights remain fixed for the complete solve: a larger trial tangent
    // must not manufacture a smaller line-search merit or convergence defect.
    work.iter().zip(assembly.rhs()).zip(fixed_row_scales).fold(
        (0.0_f64, 0.0_f64),
        |(norm, scaled), ((left, right), scale)| {
            let defect = left - right;
            (norm.hypot(defect), scaled.hypot(defect / scale))
        },
    )
}

fn relative_norm(norm: f64, reference: f64) -> f64 {
    if norm == 0.0 {
        0.0
    } else if reference > 0.0 {
        norm / reference
    } else {
        f64::INFINITY
    }
}

fn nonlinear_failure(iterations: usize, update: f64, residual: f64) -> StokesBoundaryError {
    StokesBoundaryError::Rheology(FluxRheologyError::NonlinearConvergence {
        iterations,
        relative_update: update,
        residual,
    })
}

impl<'m> FluxSystem<'m> {
    /// Solve steady `-div(2 nu_eff D(u)) + grad(p) + nu K^-1 w + beta |w| w = f`
    /// with `w = u-u_s` in the original mixed-boundary BDM1/P0 owner.
    ///
    /// Pressure and force are divided by one constant reference density;
    /// beta is explicitly authored in 1/m. Molecular `params.nu` sets Darcy
    /// resistance independently of the effective shear viscosity. No material
    /// fit is inferred from permeability. Boundary/load callbacks must be pure.
    ///
    /// A consistent symmetric tangent and the native positive cell quadrature
    /// drive bounded Newton/backtracking. Both true original and row-scaled
    /// residuals plus the full Newton correction must converge. Residual row
    /// weights and reference loads come from the zero-velocity operator; only
    /// a homogeneous zero-load norm falls back to the initial true defect.
    /// Thus the artificial Newton affine load cannot change convergence scales.
    /// The velocity reference is fixed to the larger norm of the projected
    /// initial state and first Newton candidate, without an artificial floor.
    /// Strong normal
    /// moments, continuity, traction pressure datums, slip and SIP terms retain
    /// their original owner. Reactions use a fresh accepted-state tangent/load.
    ///
    /// All-zero beta takes the original Darcy solve and diagnostics exactly.
    /// Positive beta follows the regime and no-claim limits of this module;
    /// quadrature of a general velocity norm still requires mesh convergence.
    ///
    /// # Errors
    /// Refuses malformed laws, controls/initial states, native boundary errors,
    /// unrepresentable force/work, unconverged linear solves and exhausted
    /// Newton/backtracking budgets. No partial solution is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_brinkman_forchheimer_mixed(
        &self,
        params: FluxParams,
        cells: &[ForchheimerCell],
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: ForchheimerControls,
    ) -> Result<ForchheimerSolution, StokesBoundaryError> {
        self.solve_forchheimer(
            params,
            cells,
            force,
            boundaries,
            initial_velocity,
            controls,
            None,
        )
    }

    /// Implicit BDF1 fluid inertia with the SAME nonlinear endpoint drag and
    /// Newton loop as [`Self::solve_brinkman_forchheimer_mixed`]. The accepted
    /// previous velocity supplies the original native mass/history operator;
    /// it is immutable, and the caller owns time and subsequent state commit.
    ///
    /// Coefficients, skeleton velocities and boundaries belong to the endpoint
    /// and remain frozen within the step. Rates in the result are physical
    /// endpoint exchanges, not heat increments or BDF1 numerical energy loss.
    /// Inertia anchors all-traction free-fluid components. This fixed-domain
    /// unsteady Stokes closure does not add convective or pore-scale inertia.
    ///
    /// # Errors
    /// Also refuses nonpositive/nonrepresentable dt and malformed history
    /// before callbacks. All steady physical and convergence refusals remain.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_brinkman_forchheimer_mixed_step(
        &self,
        params: FluxParams,
        cells: &[ForchheimerCell],
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: ForchheimerControls,
    ) -> Result<ForchheimerSolution, StokesBoundaryError> {
        let inv_dt = self.admit_mixed_bdf1_history(previous_velocity, timestep_s)?;
        self.solve_forchheimer(
            params,
            cells,
            force,
            boundaries,
            previous_velocity,
            controls,
            Some((previous_velocity, inv_dt)),
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn solve_forchheimer(
        &self,
        params: FluxParams,
        cells: &[ForchheimerCell],
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: ForchheimerControls,
        history: Option<(&[f64], f64)>,
    ) -> Result<ForchheimerSolution, StokesBoundaryError> {
        controls.admit(params)?;
        if initial_velocity.len() != self.n_u || initial_velocity.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "Forchheimer initial state requires exactly the finite native velocity DOFs",
            ));
        }
        let brinkman: Vec<_> = cells.iter().map(|cell| cell.brinkman).collect();
        let (drag, effective_nu) = compile_cells(&brinkman, self.mesh.tris.len(), params.nu)?;
        let mut anchors = vec![history.is_some(); self.pressure_component_count()];
        for (cell, law) in cells.iter().enumerate() {
            anchors[self.pressure_component[cell]] |= law.brinkman.permeability().is_some();
        }
        let admitted =
            self.admit_stokes_problem_with_velocity_anchors(params, boundaries, &anchors)?;
        let (previous, inv_dt) =
            history.map_or((None, 0.0), |(velocity, rate)| (Some(velocity), rate));
        let assemble = |resistance| {
            self.assemble_stokes_admitted_cell_resistance(
                params,
                &|cell, _| effective_nu[cell],
                force,
                &admitted,
                Some(resistance),
                inv_dt,
                previous,
                None,
                None,
                None,
            )
        };
        if cells.iter().all(|cell| cell.beta_per_m == 0.0) {
            let assembly = BrinkmanAssembly {
                system: self,
                stokes: assemble(CellResistance::Darcy(&drag))?,
                drag,
            };
            let flow = assembly.solve()?;
            let cell_darcy_exchange =
                assembly.cell_darcy_exchange_per_density(&flow.field.x[..self.n_u])?;
            return Ok(ForchheimerSolution {
                flow,
                iterations: 0,
                line_search_reductions: 0,
                relative_update: 0.0,
                cell_darcy_exchange,
                cell_forchheimer_exchange: vec![BrinkmanDragExchange::default(); cells.len()],
            });
        }
        let base = darcy_terms(self, &drag);
        let mut field = FluxSolution {
            x: initial_velocity
                .iter()
                .copied()
                .chain(std::iter::repeat_n(0.0, self.n - self.n_u))
                .collect(),
            n_u: self.n_u,
            iters: 0,
            rel_residual: f64::INFINITY,
            scaled_rel_residual: f64::INFINITY,
        };
        // The zero-velocity RHS is the physical forcing, including prescribed
        // skeleton motion, boundary data and old-time fluid inertia. It is not
        // the changing J(u)u-F(u) Newton affine source.
        let zero = vec![0.0; self.n_u];
        let reference = assemble(CellResistance::Linearized(&linearized_terms(
            self, cells, &base, &zero,
        )?))?;
        let fixed_row_scales: Vec<_> = (0..self.n)
            .map(|row| {
                reference
                    .matrix()
                    .row(row)
                    .1
                    .iter()
                    .fold(0.0_f64, |scale, entry| scale.max(entry.abs()))
                    .max(f64::MIN_POSITIVE)
            })
            .collect();
        // Ownership tags select constraints; never infer them from a
        // coincidentally diagonal physical row. Values come from the original
        // normal-moment quadrature and pressure-pin RHS without re-integration.
        let mut prescribed = Vec::new();
        for edge in 0..self.mesh.edges.len() {
            if self.mesh.is_boundary(edge) && !admitted.is_traction(edge) {
                prescribed.push((2 * edge, reference.rhs()[2 * edge]));
                prescribed.push((2 * edge + 1, reference.rhs()[2 * edge + 1]));
            }
        }
        for (component, &root) in self.pressure_roots.iter().enumerate() {
            if !admitted.has_traction[component] {
                prescribed.push((self.n_u + root, reference.rhs()[self.n_u + root]));
            }
        }
        for &(row, value) in &prescribed {
            field.x[row] = value;
        }
        let (load_norm, scaled_load_norm) = reference
            .rhs()
            .iter()
            .zip(&fixed_row_scales)
            .fold((0.0_f64, 0.0_f64), |(norm, scaled), (value, scale)| {
                (norm.hypot(*value), scaled.hypot(value / scale))
            });
        let mut assembly = if field.x[..self.n_u].iter().all(|value| *value == 0.0) {
            reference
        } else {
            assemble(CellResistance::Linearized(&linearized_terms(
                self, cells, &base, &field.x,
            )?))?
        };
        let mut work = vec![0.0; self.n];
        let (initial_norm, initial_scaled_norm) =
            residual_norms(&assembly, &field.x, &fixed_row_scales, &mut work);
        let residual_reference = if load_norm > 0.0 {
            load_norm
        } else {
            initial_norm
        };
        let scaled_reference = if scaled_load_norm > 0.0 {
            scaled_load_norm
        } else {
            initial_scaled_norm
        };
        if !residual_reference.is_finite() || !scaled_reference.is_finite() {
            return Err(invalid(
                "Forchheimer residual reference is not representable",
            ));
        }
        let mut velocity_reference = field.x[..self.n_u]
            .iter()
            .fold(0.0_f64, |norm, value| norm.hypot(*value));
        let (mut relative, mut scaled, mut merit) = (
            relative_norm(initial_norm, residual_reference),
            relative_norm(initial_scaled_norm, scaled_reference),
            initial_scaled_norm,
        );
        let mut iterations = 0;
        let mut reductions = 0;
        let mut update = 0.0;
        while relative > controls.relative_tolerance
            || scaled > controls.relative_tolerance
            || update > controls.relative_tolerance
        {
            if !relative.is_finite() || !scaled.is_finite() || !merit.is_finite() {
                return Err(invalid(
                    "Forchheimer nonlinear residual is not representable",
                ));
            }
            if iterations == controls.max_iterations {
                return Err(nonlinear_failure(iterations, update, relative.max(scaled)));
            }
            let linear = self.solve_linear_result(assembly.matrix(), assembly.rhs(), params)?;
            if !linear.x.iter().all(|value| value.is_finite()) || !linear.converged(params.tol) {
                return Err(StokesBoundaryError::NotConverged {
                    iterations: linear.iters,
                    relative_residual: linear.rel_residual,
                    scaled_relative_residual: linear.scaled_rel_residual,
                });
            }
            let direction: Vec<_> = linear
                .x
                .iter()
                .zip(&field.x)
                .map(|(next, old)| next - old)
                .collect();
            let numerator = direction[..self.n_u]
                .iter()
                .fold(0.0_f64, |norm, value| norm.hypot(*value));
            let denominator = linear.x[..self.n_u]
                .iter()
                .fold(0.0_f64, |norm, value| norm.hypot(*value));
            if iterations == 0 {
                velocity_reference = velocity_reference.max(denominator);
            }
            if !velocity_reference.is_finite() {
                return Err(invalid(
                    "Forchheimer velocity reference is not representable",
                ));
            }
            update = relative_norm(numerator, velocity_reference);
            iterations += 1;
            let mut step = 1.0;
            let mut accepted = None;
            for backtrack in 0..=controls.max_backtracks {
                // Preserve the native exact constrained moments on a full
                // step; old + (new-old) can lose a small new value by cancellation.
                let mut candidate: Vec<_> = if backtrack == 0 {
                    linear.x.clone()
                } else {
                    field
                        .x
                        .iter()
                        .zip(&linear.x)
                        .map(|(old, next)| (1.0 - step) * old + step * next)
                        .collect()
                };
                for &(row, value) in &prescribed {
                    candidate[row] = value;
                }
                if candidate.iter().all(|value| value.is_finite()) {
                    // Overflow in a trial constitutive state can be reduced;
                    // immutable native load/admission errors still propagate.
                    if let Ok(terms) = linearized_terms(self, cells, &base, &candidate) {
                        let endpoint = assemble(CellResistance::Linearized(&terms))?;
                        let (next_norm, next_merit) =
                            residual_norms(&endpoint, &candidate, &fixed_row_scales, &mut work);
                        let next_relative = relative_norm(next_norm, residual_reference);
                        let next_scaled = relative_norm(next_merit, scaled_reference);
                        let converged = next_relative <= controls.relative_tolerance
                            && next_scaled <= controls.relative_tolerance
                            && update <= controls.relative_tolerance;
                        if next_merit.is_finite()
                            && next_relative.is_finite()
                            && next_scaled.is_finite()
                            && ((next_merit < merit
                                && next_merit <= (1.0 - 1e-4 * step) * merit)
                                || converged)
                        {
                            accepted = Some((
                                candidate,
                                endpoint,
                                next_relative,
                                next_scaled,
                                next_merit,
                                backtrack,
                            ));
                            break;
                        }
                    }
                }
                step *= 0.5;
            }
            let Some((candidate, endpoint, next_relative, next_scaled, next_merit, backtrack)) =
                accepted
            else {
                return Err(nonlinear_failure(iterations, update, relative.max(scaled)));
            };
            field.x = candidate;
            field.iters = linear.iters;
            assembly = endpoint;
            relative = next_relative;
            scaled = next_scaled;
            merit = next_merit;
            reductions += backtrack;
        }
        if !relative.is_finite() || !scaled.is_finite() || !update.is_finite() {
            return Err(invalid("Forchheimer accepted state is not representable"));
        }
        field.rel_residual = relative;
        field.scaled_rel_residual = scaled;
        // The endpoint assembly contains J(u) and J(u)u-F(u), including its
        // constrained rows. Evaluating them reports ACTUAL nonlinear reactions.
        let flow = assembly.report_solution(field)?;
        let cell_forchheimer_exchange =
            forchheimer_exchanges(self, cells, &flow.field.x[..self.n_u])?;
        let linear_owner = BrinkmanAssembly {
            system: self,
            stokes: assembly,
            drag,
        };
        let cell_darcy_exchange =
            linear_owner.cell_darcy_exchange_per_density(&flow.field.x[..self.n_u])?;
        Ok(ForchheimerSolution {
            flow,
            iterations,
            line_search_reductions: reductions,
            relative_update: update,
            cell_darcy_exchange,
            cell_forchheimer_exchange,
        })
    }
}
